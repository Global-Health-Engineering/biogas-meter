import os
import git
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir


def plot_flow_over_time():
    df = pd.read_csv(os.path.join(get_git_root(os.getcwd()), "data", "derived_data", "mfc-only-2024-08-13.csv"), index_col=[0])
    df.index = pd.to_datetime(df.index)
    t1 = "2024-08-13 10:09:00"
    t2 = "2024-08-13 10:15:00"
    t1 = datetime.strptime(t1, "%Y-%m-%d %H:%M:%S")
    t2 = datetime.strptime(t2, "%Y-%m-%d %H:%M:%S")
    _df = df[(df.index >= t1) & (df.index <= t2)].copy()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.plot(_df.index, _df["CH4 Flow (ln/min)"])
    ax.set_xlabel("time")
    ax.set_ylabel("CH4 flow (ln/min)")
    fig.tight_layout(pad=0.2)


def main():
    plot_flow_over_time()


if __name__ == "__main__":
    main()
    plt.show()
