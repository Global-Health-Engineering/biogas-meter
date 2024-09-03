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
    y[y < 0] = 0
    y[y > 1] = 1
    return y


def plot_all_at_once(csv_file):
    df = pd.read_csv(csv_file)

    fit_deg = 3

    df["Total vol. flow round (ls/min)"] = df["Total vol. flow (ls/min)"].round(1)
    poly_coefs = {}
    for i, group in df.groupby("Total vol. flow round (ls/min)"):
        p = Polynomial.fit(group["DC voltage (VDC)"], group["CH4 Vol. fraction"], deg=fit_deg)
        poly_coefs[i] = p.convert().coef

    p = Polynomial.fit(df["DC voltage (VDC)"], df["CH4 Vol. fraction"], deg=fit_deg)
    poly_coefs_all = p.convert().coef

    df["CH4 measured vol. fraction"] = fit_func(df["DC voltage (VDC)"], poly_coefs_all)
    df["CH4 vol. fraction rel-err"] = (df["CH4 Vol. fraction"] - df["CH4 measured vol. fraction"]) / df["CH4 Vol. fraction"]

    fig, axes = plt.subplots(2, 1, figsize=(5, 5.5), sharex=True)
    ax0, ax1 = axes

    norm = mpl.colors.Normalize(vmin=0, vmax=max(df["Total vol. flow (ls/min)"]))
    cmap = plt.get_cmap('viridis')
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    ax0.errorbar(_df["DC voltage (VDC)"], _df["CH4 Vol. fraction"],
                 xerr = 1.96*_df["DC voltage STUNC (VDC)"],
                 yerr = 1.96*_df["CH4 Flow (ls/min) STUNC"],
                 fmt = "o",
                 mfc = "None",
                 mec = "None",
                 ecolor = "k",
                 capsize = 2,
                 elinewidth = 0.8,
                 zorder = -10)
    ax0.scatter(df["DC voltage (VDC)"], df["CH4 Vol. fraction"],
                c = df["Total vol. flow (ls/min)"],
                cmap = cmap,
                norm = norm,
                s = 70,
                alpha=0.9,
                edgecolors = "k",
                zorder=10)

    x = np.arange(-0.1, 1.01, 0.01)
    for k, v in poly_coefs.items():
        line_colors = cmap(norm(k))
        ax0.plot(x, fit_func(x, poly_coefs[k]),
                 zorder = 0,
                 color = line_colors)
    ax0.plot(x, fit_func(x, poly_coefs_all),
             zorder = 100,
             color = "k",
             linestyle = "-",
             linewidth = 2)
    ax1.scatter(df["DC voltage (VDC)"], 100*df["CH4 vol. fraction rel-err"],
                c = df["Total vol. flow (ls/min)"],
                cmap = cmap,
                norm = norm,
                s = 70,
                alpha=0.9,
                edgecolors = "k")
    x = np.arange(-0.1, 0.31, 0.01)
    y1 = [-5]*len(x)
    y2 = [5]*len(x)
    ax1.fill_between(x, y1, y2, color="lightgrey", zorder=-10)
    ax1.plot(x, y1, color="k", linewidth=1., zorder=0)
    ax1.plot(x, y2, color="k", linewidth=1., zorder=0)
    ax1.hlines(0, -0.1, 0.31, linestyle=":", linewidth=1., color="k", zorder=0)

    ax0.set_xlim((-0.1, 0.3))
    ax0.set_ylim((-0.1, 1.1))
    ax0.set_ylabel(r"$\phi_{\rm CH_4}$")

    ax1.set_ylim((-1e3, 1e3))
    ax1.set_xlabel(r"$U~/~{\rm V_{DC}}$")
    ax1.set_ylabel(r"$(\phi_{\rm CH_4}^{\rm ref} - \phi_{\rm CH_4}^{\rm m})~(\phi_{\rm CH_4}^{\rm ref})^{-1}$")
    ax1.set_yscale("symlog")
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(1e-1))

    fig.subplots_adjust(left=0.19, right=0.95, bottom=0.15, top=0.95)
    cbar = fig.colorbar(colormap, ax=axes)
    cbar.set_label(r"$\dot{V}~/~{\rm l_{\rm s}~min^{-1}}$")
    cbar.ax.tick_params(axis='y', direction='out')


def plot_all_at_once_with_humidity_markers(csv_file):
    df = pd.read_csv(csv_file)

    fit_deg = 3

    df["Total vol. flow round (ls/min)"] = df["Total vol. flow (ls/min)"].round(1)
    poly_coefs = {}
    for i, group in df.groupby("Total vol. flow round (ls/min)"):
        p = Polynomial.fit(group["DC voltage (VDC)"], group["CH4 Vol. fraction"], deg=fit_deg)
        poly_coefs[i] = p.convert().coef

    p = Polynomial.fit(df["DC voltage (VDC)"], df["CH4 Vol. fraction"], deg=fit_deg)
    poly_coefs_all = p.convert().coef

    df["CH4 measured vol. fraction"] = fit_func(df["DC voltage (VDC)"], poly_coefs_all)
    df["CH4 vol. fraction rel-err"] = (df["CH4 Vol. fraction"] - df["CH4 measured vol. fraction"]) / df["CH4 Vol. fraction"]

    rh_ranges = ((0, 15), (41, 60), (61, 80))
    markers = ("o", "s", "v")
    
    fig, axes = plt.subplots(2, 1, figsize=(5, 5.5), sharex=True)
    ax0, ax1 = axes

    norm = mpl.colors.Normalize(vmin=0, vmax=max(df["Total vol. flow (ls/min)"]))
    cmap = plt.get_cmap('viridis')
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    for rh, marker in zip(rh_ranges, markers):
        _df = df[(df["AHT Humidity (%rH)"] >= rh[0]) & (df["AHT Humidity (%rH)"] <= rh[1])].copy()
        ax0.errorbar(_df["DC voltage (VDC)"], _df["CH4 Vol. fraction"],
                     xerr = 1.96*_df["DC voltage STUNC (VDC)"],
                     yerr = 1.96*_df["CH4 Flow (ls/min) STUNC"],
                     fmt = "o",
                     mfc = "None",
                     mec = "None",
                     ecolor = "k",
                     capsize = 2,
                     elinewidth = 0.8,
                     zorder = -10)
        ax0.scatter(_df["DC voltage (VDC)"], _df["CH4 Vol. fraction"],
                    c = _df["Total vol. flow (ls/min)"],
                    cmap = cmap,
                    norm = norm,
                    marker = marker,
                    s = 70,
                    alpha = 0.9,
                    edgecolors = "k",
                    zorder = 10)
        
        ax1.scatter(_df["DC voltage (VDC)"], 100*_df["CH4 vol. fraction rel-err"],
                    c = _df["Total vol. flow (ls/min)"],
                    cmap = cmap,
                    norm = norm,
                    marker = marker,
                    s = 70,
                    alpha = 0.9,
                    edgecolors = "k")

    x = np.arange(-0.1, 1.01, 0.01)
    ax0.plot(x, fit_func(x, poly_coefs_all),
                 zorder = 100,
                 color = "k",
                 linestyle = "-",
                 linewidth = 2)
    for k, v in poly_coefs.items():
        line_colors = cmap(norm(k))
        ax0.plot(x, fit_func(x, poly_coefs[k]),
                 zorder = 0,
                 color = line_colors)

    x = np.arange(-0.1, 0.31, 0.01)
    y1 = [-5]*len(x)
    y2 = [5]*len(x)
    ax1.fill_between(x, y1, y2, color="lightgrey", zorder=-10)
    ax1.plot(x, y1, color="k", linewidth=1., zorder=0)
    ax1.plot(x, y2, color="k", linewidth=1., zorder=0)
    ax1.hlines(0, -0.1, 0.31, linestyle=":", linewidth=1., color="k", zorder=0)

    ax0.set_xlim((-0.1, 0.3))
    ax0.set_ylim((-0.1, 1.1))
    ax0.set_yticks(np.arange(0, 1.25, 0.25))
    ax0.set_ylabel(r"$\phi_{\rm CH_4}^{\rm ref}$")

    ax1.set_ylim((-1e3, 1e3))
    ax1.set_xlabel(r"$U~/~{\rm V_{DC}}$")
    ax1.set_ylabel(r"$(\phi_{\rm CH_4}^{\rm ref} - \phi_{\rm CH_4}^{\rm m})~(\phi_{\rm CH_4}^{\rm ref})^{-1}$")
    ax1.set_yscale("symlog")
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(1e-1))

    fig.subplots_adjust(left=0.19, right=0.95, bottom=0.15, top=0.95)
    cbar = fig.colorbar(colormap, ax=axes)
    cbar.set_label(r"$\dot{V}~/~{\rm l_{\rm s}~min^{-1}}$")
    cbar.ax.tick_params(axis='y', direction='out')

    fig.savefig(os.path.join(get_git_root(os.getcwd()), "img", f"{os.path.basename(csv_file).split('.')[0]}.pdf"))


def plot_calibration_for_3_humidities(csv_file):
    df = pd.read_csv(csv_file)

    fit_deg = 4

    df["Total vol. flow round (ls/min)"] = df["Total vol. flow (ls/min)"].round(1)
    poly_coefs = {}
    for i, group in df.groupby("Total vol. flow round (ls/min)"):
        p, (resid, rank, sv, rcond) = Polynomial.fit(group["DC voltage (VDC)"],
                                                     group["CH4 Vol. fraction"],
                                                     deg = fit_deg,
                                                     full = True)
        poly_coefs[i] = p.convert().coef
        print(f"{i} lspm SSQ = {resid[0]:.2e}")

    p, (resid, rank, sv, rcond) = Polynomial.fit(df["DC voltage (VDC)"],
                                                 df["CH4 Vol. fraction"],
                                                 deg = fit_deg,
                                                 full = True)
    poly_coefs_all = p.convert().coef
    print(f"Total SSQ = {resid[0]:.2e}")

    df["CH4 measured vol. fraction"] = fit_func(df["DC voltage (VDC)"], poly_coefs_all)
    df["CH4 vol. fraction rel-err"] = (df["CH4 Vol. fraction"] - df["CH4 measured vol. fraction"]) / df["CH4 Vol. fraction"]
    print(f"Mean composition error: {100*df['CH4 vol. fraction rel-err'][df['CH4 Vol. fraction'] > 0.1].abs().mean():.2f}%")

    fig, axes = plt.subplots(2, 3, figsize=(9, 5.5), sharex="col", sharey="row")
    [ax0, ax1, ax2], [ax3, ax4, ax5] = axes

    norm = mpl.colors.Normalize(vmin=0, vmax=max(df["Total vol. flow (ls/min)"]))
    cmap = plt.get_cmap('viridis')
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    rh_ranges = ((0, 15), (41, 60), (61, 80))
    ax_pairs = ((ax0, ax3), (ax1, ax4), (ax2, ax5))
    x = np.arange(-0.1, 1.01, 0.01)
    y1 = [-5]*len(x)
    y2 = [5]*len(x)

    for rh, ax_pair in zip(rh_ranges, ax_pairs):
        _df = df[(df["AHT Humidity (%rH)"] >= rh[0]) & (df["AHT Humidity (%rH)"] <= rh[1])].copy()
        ax_pair[0].errorbar(_df["DC voltage (VDC)"], _df["CH4 Vol. fraction"],
                            xerr = 1.96*_df["DC voltage STUNC (VDC)"],
                            yerr = 1.96*_df["CH4 Flow (ls/min) STUNC"],
                            fmt = "o",
                            mfc = "None",
                            mec = "None",
                            ecolor = "k",
                            capsize = 2,
                            elinewidth = 0.8,
                            zorder = -10)
        ax_pair[0].scatter(_df["DC voltage (VDC)"], _df["CH4 Vol. fraction"],
                           c = _df["Total vol. flow (ls/min)"],
                           cmap = cmap,
                           norm = norm,
                           marker = "o",
                           s = 70,
                           alpha = 0.9,
                           edgecolors = "k",
                           zorder = 10)
        
        ax_pair[1].scatter(_df["DC voltage (VDC)"], 100*_df["CH4 vol. fraction rel-err"],
                           c = _df["Total vol. flow (ls/min)"],
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

        ax_pair[0].plot(x, fit_func(x, poly_coefs_all),
                        zorder = 100,
                        color = "k",
                        linestyle = "-",
                        linewidth = 1)

        ax_pair[0].set_xlim((-0.1, 0.3))
        ax_pair[0].set_ylim((-0.1, 1.1))
        ax_pair[0].set_xticks(np.arange(-0.1, 0.35, 0.1))
        ax_pair[0].set_yticks(np.arange(0, 1.25, 0.25))

        ax_pair[1].set_ylim((-1e3, 1e3))
        ax_pair[1].set_yscale("symlog")
        ax_pair[1].fill_between(x, y1, y2, color="lightgrey", zorder=-10)
        ax_pair[1].plot(x, y1, color="k", linewidth=1., zorder=0)
        ax_pair[1].plot(x, y2, color="k", linewidth=1., zorder=0)
        ax_pair[1].hlines(0, -0.1, 0.31, linestyle=":", linewidth=1., color="k", zorder=0)

    ax0.set_title(r"$0~{\rm \%rh}$")
    ax1.set_title(r"$50 \pm 5~{\rm \%rh}$")
    ax2.set_title(r"$72 \pm 5~{\rm \%rh}$")

    ax0.set_ylabel(r"$\phi_{\rm CH_4}$")
    ax4.set_xlabel(r"$U~/~{\rm V_{DC}}$")
    ax3.set_ylabel(r"$(\phi_{\rm CH_4}^{\rm ref} - \phi_{\rm CH_4}^{\rm m})~(\phi_{\rm CH_4}^{\rm ref})^{-1}$")
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(1e-1))

    fig.subplots_adjust(left=0.1, right=1.05, bottom=0.13, top=0.93)

    cbar = fig.colorbar(colormap, ax=axes)
    cbar.set_label(r"$\dot{V}~/~{\rm l_{\rm s}~min^{-1}}$")
    cbar.ax.tick_params(axis='y', direction='out')

    fig.savefig(os.path.join(get_git_root(os.getcwd()), "img", f"{os.path.basename(csv_file).split('.')[0]}_6figs.pdf"))


if __name__ == "__main__":
    set_mpl()
    csv_file = os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "composition_calibration.csv")
    
    # plot_all_at_once(csv_file)
    # plot_all_at_once_with_humidity_markers(csv_file)
    plot_calibration_for_3_humidities(csv_file)
    
    plt.show()


