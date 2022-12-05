import logging
import os
import json
import sys
import re
import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from typing import List
from pathlib import Path
from parse_neighbours import read_neighbours


def main():
    logging.basicConfig(level=logging.INFO)

    neighpath = Path(sys.argv[1]).resolve()
    datafile = Path(sys.argv[2]).resolve()
    outfile = Path(sys.argv[3]).resolve()
    start_season = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    n_seasons = int(sys.argv[5]) if len(sys.argv) > 5 else None

    plot_adhesion(neighpath,
                  datafile,
                  outfile,
                  start_season,
                  n_seasons)


def plot_adhesion(neighpath, datafile, outfile, start_season, n_seasons):
    logging.info("Reading neighbours file")
    if os.path.isfile(neighpath):
        with open(neighpath, "r") as file:
            neighs = json.load(file)
        season_filter = build_season_filter(neighs, start_season, n_seasons)
        neighs = {k: v for k, v in neighs.items() if int(k) in season_filter}
    else:
        logging.warning("Neighbours file not found")
        logging.info(f"Try to build neighbours from '{neighpath}'")
        seasons = []
        for filepath in Path(neighpath).iterdir():
            m = re.search(r"\d+", filepath.name)
            if m is not None:
                seasons.append(m.group())
        season_filter = build_season_filter(seasons, start_season, n_seasons)
        neighs = read_neighbours(neighpath, season_filter)

    logging.info(f"Reading adhesion information from {datafile}")
    adhesion_info = parse_adhesion_info(datafile, season_filter)

    logging.info(f"Plotting figure to {outfile}")
    fig = make_plots(neighs, adhesion_info)

    if "html" in outfile.name:
        fig.write_html(outfile)
    else:
        fig.write_image(outfile)
    logging.info("Finished")


def make_plots(neighs, adhesion_info):
    fig = make_subplots(3, 1, subplot_titles=["Prevalence of multicellularity",
                                              "Medium gamma of neighbouring cells",
                                              "Medium gamma of all x all cells"])

    logging.info("Plotting prevalence of multicellularity")
    x1 = []
    y1 = []
    for season, clusters in sorted(neighs.items(), key=lambda kv: int(kv[0])):
        x1.append(int(season))
        unicells = 0
        total_pop = 0
        for c in clusters:
            if len(c) == 1:
                unicells += 1
            total_pop += len(c)
        y1.append(1 - unicells / total_pop)

    fig.add_trace(go.Scatter(
        x=x1,
        y=y1,
        mode="lines"
    ), row=1, col=1)
    fig.update_yaxes(title="% of adhering cells", range=[0, 1], row=1, col=1)

    sorted_adhesion = sorted(adhesion_info.items(), key=lambda kv: int(kv[0]))

    logging.info("Plotting medium gamma of neighbouring cells")
    x2 = []
    y2 = []
    for season, cells in sorted_adhesion:
        x2.append(int(season))
        gammas = []
        for c_info in cells.values():
            gammas += [c_info["medJ"] - cJ/2 for cJ in c_info["cellJs"].values()]
        y2.append(np.median(gammas))

    fig.add_trace(go.Scatter(
        x=x2,
        y=y2,
        mode="lines"
    ), row=2, col=1)
    fig.update_yaxes(title="gamma", row=2, col=1)

    logging.info("Plotting medium gamma of all x all cells")
    x3 = []
    y3 = []
    for season, cells in sorted_adhesion:
        x3.append(int(season))

        keylock_matrix = np.array([[c_info["key"], c_info["lock"]] for c_info in cells.values()])
        y3.append(get_medium_gamma(keylock_matrix, 4, 8, np.array([4, 3, 2, 1, 1, 1]), 50))
    fig.add_trace(go.Scatter(
        x=x3,
        y=y3,
        mode="lines",
    ), row=3, col=1)
    fig.update_yaxes(title="gamma", row=3, col=1)

    return fig.update_layout(showlegend=False)


def parse_adhesion_info(datafile, season_filter: List[int] = None):
    seasons = {}
    with open(datafile) as file:
        text = file.read()

    matches = re.finditer(
        r"^(\d+) (\d+) [\d. -]+ ([01]{24}) ([01]{24}) [\d.]+ 0 (\d+) ((?:\d+ \d+ )*)",
        text,
        flags=re.MULTILINE
    )
    for m in matches:
        print(m.group())
        season_i = int(m.group(1))
        if season_i not in season_filter:
            continue

        season = seasons.setdefault(season_i, {})

        key = np.array(list(m.group(3)), dtype=int)
        lock = np.array(list(m.group(4)), dtype=int)
        cellJs = {}
        for c_cJ in re.finditer(r"(\d+) (\d+)", m.group(6)):
            cellJs[int(c_cJ.group(1))] = int(c_cJ.group(2))

        season[int(m.group(2))] = {
            "key": key,
            "lock": lock,
            "medJ": int(m.group(5)),
            "cellJs": cellJs
        }

    return seasons


def build_season_filter(seasons, start_season, n_seasons=None):
    seasons = set(int(season) for season in seasons if int(season) >= start_season)
    step = round(len(seasons) / n_seasons) if n_seasons is not None else 1
    return sorted(seasons)[::max(1, step)]


def get_medium_gamma(key_locks, ja, jam, f_arr, n=None):
    """Calculates the medium gamma from all x all cells.

    :param key_locks: Numpy matrix of (key, lock) pairs
    :param ja:
    :param jam:
    :param f_arr:
    :param n: Samples n key-lock pairs from total population to speed up calculations
    :return:
    """
    if n is not None:
        np.random.seed(0)
        indexes = np.random.randint(0, len(key_locks), n)
        key_locks = key_locks[indexes]

    gammas = []
    for k1, l1 in key_locks:
        for k2, l2 in key_locks:
            # print(np.all(k1 == k2) and np.all(l1 == l2))
            #print(l1, k1, l2, k2)
            gammas.append(gamma(ja, jam, f_arr, l1, k1, l2, k2))

    return np.mean(gammas)


def cell_cell_j(ja, l1, k1, l2, k2):
    v = len(k1)
    return ja + 2 * v - np.count_nonzero(l1 != k2) - np.count_nonzero(l2 != k1)


def cell_med_j(jam, f_arr, key):
    return jam + np.sum(f_arr * key[:len(f_arr)])


def gamma(ja, jam, f_seq, l1, k1, l2, k2):
    # print(cell_med_j(jam, f_seq, k1), cell_cell_j(ja, l1, k1, l2, k2))
    return cell_med_j(jam, f_seq, k1) - cell_cell_j(ja, l1, k1, l2, k2) / 2


if __name__ == "__main__":
    main()
