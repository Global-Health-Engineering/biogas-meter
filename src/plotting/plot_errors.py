import os
import git
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx


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
    ax.set_xlabel(r"$\phi_{\rm CH_4} = \dot{V}_{\rm CH_4}~\dot{V}^{-1}$")
    ax.set_ylabel(r"$100~(\dot{V}_{\rm ref} - \dot{V}_{\rm sbg})~\dot{V}_{\rm ref}^{-1}$")

    cbar.ax.tick_params(axis='y', direction='out')
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    # ax.set_xlim((0.2, 0.8))
    # ax.set_ylim((-60, 0))
    plt.tight_layout(pad=0.2)


def plot_flow_vals_unc(csv_file):
    df = pd.read_csv(csv_file, index_col="date/yyyy-mm-dd")

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))


def main():
    csv_file = os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "full_sgb_dry_run.csv")

    plot_rel_err(csv_file)
    # plot_flow_vals_unc(csv_file)


if __name__ == "__main__":
    set_mpl()
    main()
    plt.show()
