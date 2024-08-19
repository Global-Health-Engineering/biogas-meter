import os
import git
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib.ticker import Locator


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        self.linthresh = linthresh

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a %s type.' % type(self))

    def __call__(self):
        majorlocs = self.axis.get_majorticklocs()
        minorlocs = []
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)
        return self.raise_if_exceeds(np.array(minorlocs))


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


def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir


def plot_rel_err(csv_file):
    df = pd.read_csv(csv_file, index_col="date/yyyy-mm-dd")
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    norm = mpl.colors.Normalize(vmin = min(df["Total vol. flow (ln/min)"]),
                                vmax = max(df["Total vol. flow (ln/min)"]))
    cmap = "viridis"
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    ax.scatter(df["CH4 Vol. fraction"], 100*df["SBG rel-err"],
               c = df["Total vol. flow (ln/min)"],
               cmap = "viridis",
               norm = norm,
               s = 70,
               alpha=0.9,
               edgecolors = "k")
    
    cbar = plt.colorbar(colormap, ax=ax)
    cbar.set_label(r"$\dot{V}~/~{\rm l_n~min^{-1}}$")
    cbar.ax.tick_params(axis='y', direction='out')

    ax.set_xlabel(r"$\phi_{\rm CH_4} = \dot{V}_{\rm CH_4}~\dot{V}^{-1}$")
    ax.set_ylabel(r"$100~(\dot{V}_{\rm ref} - \dot{V}_{\rm sbg})~\dot{V}_{\rm ref}^{-1}$")
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    # ax.set_xlim((0.2, 0.8))
    # ax.set_ylim((-60, 0))
    plt.tight_layout(pad=0.2)


def plot_flow_vals_unc(csv_file):
    df = pd.read_csv(csv_file, index_col="date/yyyy-mm-dd")

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    norm = mpl.colors.Normalize(vmin = 0, vmax = 1)
    cmap = "viridis"
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    
    x = np.arange(0, 12.1, 0.1)
    ax.plot(x, x, linestyle=":", color="k")

    ax.errorbar(df["Total vol. flow (ln/min)"], df["Flow SBG [ln/min]"],
                xerr = 2*df["Total vol. flow (ln/min) COMB-UNC"],
                yerr = 2*df["Flow SBG [ln/min] STUNC"],
                fmt = "o",
                mfc = "None",
                mec = "None",
                ecolor = "k",
                capsize = 2,
                elinewidth = 0.8,
                zorder = -10)
    ax.scatter(df["Total vol. flow (ln/min)"], df["Flow SBG [ln/min]"],
               c = df["CH4 Vol. fraction"],
               cmap = "viridis",
               norm = norm,
               s = 70,
               alpha=0.9,
               edgecolors = "k",
               zorder = 10)

    cbar = plt.colorbar(colormap, ax=ax)
    cbar.set_label(r"$\phi_{\rm CH_4}$")
    cbar.ax.tick_params(axis='y', direction='out')
    
    ax.set_xlabel(r"$\dot{V}_{\rm ref}~/~{\rm l_n~min^{-1}}$")
    ax.set_ylabel(r"$\dot{V}_{\rm sbg}~/~{\rm l_n~min^{-1}}$")
    ax.set_xlim((0, 10))
    ax.set_ylim((0, 10))
    ax.set_xticks(np.arange(0, 11, 2))
    
    plt.tight_layout(pad=0.2)


def plot_rel_and_abs_unc(csv_file):
    df = pd.read_csv(csv_file, index_col="date/yyyy-mm-dd")
    
    fig, axes = plt.subplots(2, 1, figsize=(5, 5), sharex=True)
    ax0, ax1 = axes

    norm = mpl.colors.Normalize(vmin = 0, vmax = 1)
    cmap = "viridis"
    colormap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    x = np.arange(0, 12.1, 0.1)
    ax0.plot(x, x, linestyle=":", color="k")
    ax0.errorbar(df["Total vol. flow (ln/min)"], df["Flow SBG [ln/min]"],
                xerr = 2*df["Total vol. flow (ln/min) COMB-UNC"],
                yerr = 2*df["Flow SBG [ln/min] STUNC"],
                fmt = "o",
                mfc = "None",
                mec = "None",
                ecolor = "k",
                capsize = 2,
                elinewidth = 0.8,
                zorder = -10)
    ax0.scatter(df["Total vol. flow (ln/min)"], df["Flow SBG [ln/min]"],
               c = df["CH4 Vol. fraction"],
               cmap = "viridis",
               norm = norm,
               s = 70,
               alpha=0.9,
               edgecolors = "k",
               zorder = 10)

    ax0.set_ylabel(r"$\dot{V}_{\rm sbg}~/~{\rm l_n~min^{-1}}$")
    ax0.set_xlim((0, 10))
    ax0.set_ylim((0, 10))
    ax0.set_xticks(np.arange(0, 11, 2))
    ax0.set_yticks(np.arange(0, 11, 2))

    ax1.scatter(df["Total vol. flow (ln/min)"], 100*df["SBG rel-err"],
               c = df["CH4 Vol. fraction"],
               cmap = "viridis",
               norm = norm,
               s = 70,
               alpha=0.9,
               edgecolors = "k")
    
    fig.subplots_adjust(left=0.18, right=0.9, bottom=0.13, top=0.97)
    cbar = fig.colorbar(colormap, ax=axes)
    cbar.set_label(r"$\phi_{\rm CH_4} = \dot{V}_{\rm CH_4}~\dot{V}^{-1}$")
    cbar.ax.tick_params(axis='y', direction='out')

    ax1.set_xlabel(r"$\dot{V}~/~{\rm l_n~min^{-1}}$")
    ax1.set_ylabel(r"$100~(\dot{V}_{\rm ref} - \dot{V}_{\rm sbg})~\dot{V}_{\rm ref}^{-1}$")
    ax1.set_xticks(np.arange(0, 11, 2))
    ax1.set_xlim((0, 8))
    ax1.set_ylim((-1e3, 0))
    
    ax1.set_yscale("symlog")
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(1e-1))


def main():
    csv_file = os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "full_sgb_dry_run.csv")

    # plot_rel_err(csv_file)
    # plot_flow_vals_unc(csv_file)
    plot_rel_and_abs_unc(csv_file)


if __name__ == "__main__":
    set_mpl()
    main()
    plt.show()