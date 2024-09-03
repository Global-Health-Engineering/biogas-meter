import os
import numpy as np
from numpy.polynomial.polynomial import Polynomial
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import CoolProp as CP
from plot_errors_sbg import MinorSymLogLocator
from get_git_root import get_git_root

def set_mpl():
    fontsize = 15
    mpl.rcParams.update({"font.family": "Times New Roman",
                         "mathtext.fontset": "dejavuserif",
                         "font.size": fontsize,
                         "axes.linewidth": 1.,
                         "axes.labelsize": fontsize,
                         "axes.titlesize" : fontsize,
                         "legend.fontsize": fontsize,
                         "xtick.top": True,
                         "xtick.bottom": True,
                         "ytick.left": True,
                         "ytick.right": True,
                         "xtick.direction": "in",
                         "ytick.direction": "in",
                         "xtick.major.pad": 7,
                         "ytick.major.pad": 7})


def fit_func(x, coefs):
    y = 0 if type(x) is int else [0]*len(x)
    x, y = np.array(x), np.array(y)
    for power, coef in enumerate(coefs):
        y = y + coef * x**power
    return y


def plot_calibration_for_3_humidities(csv_file):
    df = pd.read_csv(csv_file)

    fit_deg = 1

    df["CH4 Vol. fraction round"] = df["CH4 Vol. fraction"].round(1)
    df["YF201 vol. flow (ls/min)"] = pd.NA

    poly_coefs = {}

    # --- fit to all data ---

    # for i, group in df.groupby("CH4 Vol. fraction round"):
    #     p, (resid, rank, sv, rcond) = Polynomial.fit(group["YF201 pulses"],
    #                                                  group["Total vol. flow (ls/min)"],
    #                                                  deg = fit_deg,
    #                                                  full = True)
    #     poly_coefs[i] = p.convert().coef
    #     df.loc[group.index, "YF201 vol. flow (ls/min)"] = fit_func(group["YF201 pulses"], poly_coefs[i])
    #     print(f"{100*i}% CH4 SSQ = {resid[0]:.2e}")

    # --- fit to the highest humidity only ---

    for i, group in df[(df["AHT Humidity (%rH)"] >= 61) & (df["AHT Humidity (%rH)"] <= 80)].groupby("CH4 Vol. fraction round"):
        p, (resid, rank, sv, rcond) = Polynomial.fit(group["YF201 pulses"],
                                                     group["Total vol. flow (ls/min)"],
                                                     deg = fit_deg,
                                                     full = True)
        poly_coefs[i] = p.convert().coef
    for i, group in df.groupby("CH4 Vol. fraction round"):
        df.loc[group.index, "YF201 vol. flow (ls/min)"] = fit_func(group["YF201 pulses"], poly_coefs[i])
    
    # ---

    df["YF201 vol. flow rel-err"] = (df["Total vol. flow (ls/min)"] - df["YF201 vol. flow (ls/min)"]) / df["Total vol. flow (ls/min)"]
    print(f"Mean flow measurement error: {100*df['YF201 vol. flow rel-err'].abs().mean():.2f} %")

    fig, axes = plt.subplots(2, 3, figsize=(9, 5.5), sharex="col", sharey="row")
    [ax0, ax1, ax2], [ax3, ax4, ax5] = axes

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cmap = plt.get_cmap('viridis')
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    rh_ranges = ((0, 15), (41, 60), (61, 80))
    ax_pairs = ((ax0, ax3), (ax1, ax4), (ax2, ax5))
    x = np.arange(0, 150.01, 0.01)
    y1 = [-5]*len(x)
    y2 = [5]*len(x)

    for rh, ax_pair in zip(rh_ranges, ax_pairs):
        _df = df[(df["AHT Humidity (%rH)"] >= rh[0]) & (df["AHT Humidity (%rH)"] <= rh[1])].copy()
        ax_pair[0].errorbar(_df["YF201 pulses"], _df["Total vol. flow (ls/min)"],
                            xerr = 1.96*_df["YF201 pulses STUNC"],
                            yerr = 1.96*_df["Total vol. flow (ls/min) COMB-UNC"],
                            fmt = "o",
                            mfc = "None",
                            mec = "None",
                            ecolor = "k",
                            capsize = 2,
                            elinewidth = 0.8,
                            zorder = -10)
        ax_pair[0].scatter(_df["YF201 pulses"], _df["Total vol. flow (ls/min)"],
                           c = _df["CH4 Vol. fraction"],
                           cmap = cmap,
                           norm = norm,
                           marker = "o",
                           s = 70,
                           alpha = 0.9,
                           edgecolors = "k",
                           zorder = 10)
        
        ax_pair[1].scatter(_df["YF201 pulses"], 100*_df["YF201 vol. flow rel-err"],
                           c = _df["CH4 Vol. fraction"],
                           cmap = cmap,
                           norm = norm,
                           marker = "o",
                           s = 70,
                           alpha = 0.9,
                           edgecolors = "k")

        for k, v in poly_coefs.items():
            line_colors = cmap(norm(k))
            ax_pair[0].plot(x, fit_func(x, poly_coefs[k]),
                            zorder = 0,
                            color = line_colors)

        ax_pair[0].set_xlim((0, 150))
        ax_pair[0].set_ylim((0, 8))

        ax_pair[1].set_ylim((-1e3, 1e3))
        ax_pair[1].set_yscale("symlog")
        ax_pair[1].fill_between(x, y1, y2, color="lightgrey", zorder=-10)
        ax_pair[1].plot(x, y1, color="k", linewidth=1., zorder=0)
        ax_pair[1].plot(x, y2, color="k", linewidth=1., zorder=0)
        ax_pair[1].hlines(0, 0, 150, linestyle=":", linewidth=1., color="k", zorder=0)

    ax0.set_title(r"$0~{\rm \%rh}$")
    ax1.set_title(r"$50 \pm 5~{\rm \%rh}$")
    ax2.set_title(r"$72 \pm 5~{\rm \%rh}$")

    ax0.set_ylabel(r"$\dot{V}_{\rm ref}~/~{\rm l_{\rm s}~min^{-1}}$")
    ax4.set_xlabel(r"$N$")
    ax3.set_ylabel(r"$100~(\dot{V}_{\rm ref} - \dot{V}_{\rm m})~\dot{V}_{\rm ref}^{-1}$")
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(1e-1))

    fig.subplots_adjust(left=0.1, right=1.03, bottom=0.11, top=0.93)

    cbar = fig.colorbar(colormap, ax=axes)
    cbar.set_label(r"$\dot{V}~/~{\rm l_{\rm s}~min^{-1}}$")
    cbar.ax.tick_params(axis='y', direction='out')

    fig.savefig(os.path.join(get_git_root(os.getcwd()), "img", f"{os.path.basename(csv_file).split('.')[0]}_6figs.pdf"))

if __name__ == "__main__":
    set_mpl()
    csv_file = os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "yf201_calibration.csv")
    plot_calibration_for_3_humidities(csv_file)
    
    plt.show()



