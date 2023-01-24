import argparse
import logging
import sys
import warnings
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from typing import List, Set
from plotly.subplots import make_subplots
from enlighten import Counter
from parse import *

# Type alias
CellCluster = List[Set[int]]
logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="plot_adhesion",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("datadir", help="directory containing the cell CSV files")
    parser.add_argument("outputfile", help="output HTML or SVG file")
    parser.add_argument("-t",
                        "--time-step",
                        help="first timestep to plot",
                        default=0,
                        type=int)
    parser.add_argument("-n",
                        help="number of time-steps to plot (can be used to speed up plotting)",
                        default=None,
                        type=int)
    args = parser.parse_args()

    t_filter = build_time_filter(get_time_points(args.datadir), start=args.time_step, n=args.n)
    celldf = parse_cell_data(args.datadir, time_filter=t_filter)
    plot_adhesion(celldf, args.outputfile)

    logger.info("Finished")


def plot_adhesion(celldf: pd.DataFrame, outputfile):
    logger.info("Making adhesion plots")

    fig = make_subplots(2,
                        1,
                        shared_xaxes=True,
                        subplot_titles=["Prevalence of multicellularity",
                                        "Median gamma of neighbouring cells"])
    x = []
    multicel_perc = []
    gamma_plots = {
        "all x all": [],
        "mig x mig": [],
        "div x div": [],
        "mig x div": []
    }

    time_steps = celldf["time"].unique()
    pbar = Counter(total=len(time_steps), desc='Time-steps')
    for time_step in time_steps:
        x.append(int(time_step))
        sdf = celldf[celldf["time"] == time_step]

        clusters = get_adhering_clusters(sdf)
        unicells = 0
        total_pop = 0
        for c in clusters:
            if len(c) == 1:
                unicells += 1
            total_pop += len(c)
        multicel_perc.append(1 - unicells / total_pop)

        neigh_gammas = []
        mig_mig_gammas = []
        div_div_gammas = []
        mig_div_gammas = []
        for sigma, tau, Jmed, neighs, Jneighs in zip(sdf["sigma"],
                                                     sdf["tau"],
                                                     sdf["Jmed"],
                                                     sdf["neighbour_list"],
                                                     sdf["Jneighbour_list"]):
            for neigh, neighJ in zip(np.fromstring(neighs, sep=' ', dtype=int),
                                     np.fromstring(Jneighs, sep=' ')):
                if neigh == 0:
                    continue
                gamma = Jmed - neighJ / 2
                if tau == sdf.loc[(time_step, neigh), "tau"]:
                    if tau == 1:
                        mig_mig_gammas.append(gamma)
                    else:
                        div_div_gammas.append(gamma)
                else:
                    mig_div_gammas.append(gamma)
                neigh_gammas.append(gamma)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gamma_plots["all x all"].append(np.median(neigh_gammas))
            gamma_plots["mig x mig"].append(np.median(mig_mig_gammas))
            gamma_plots["div x div"].append(np.median(div_div_gammas))
            gamma_plots["mig x div"].append(np.median(mig_div_gammas))

        pbar.update()

    colors = px.colors.qualitative.Plotly
    fig.add_trace(go.Scatter(
        x=x,
        y=multicel_perc,
        mode="lines",
        showlegend=False
    ), row=1, col=1)
    for i, (name, ys) in enumerate(gamma_plots.items()):
        fig.add_trace(go.Scatter(
            x=x,
            y=ys,
            mode="lines",
            name=name,
            line_color=colors[i]
        ), row=2, col=1)
    fig.update_yaxes(title="% of adhering cells", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title="gamma", row=2, col=1)

    logger.info(f"Writing plot to: {outputfile}")
    if outputfile[-5:] == ".html":
        fig.write_html(outputfile)
    else:
        fig.write_image(outputfile)


def get_adhering_neighbours(celldf: pd.DataFrame):
    """Return a list of lists with the neighbours that each cell is adhering to."""
    if not celldf["time"].nunique() == 1:
        raise ValueError("make sure this function is called for a single time point")

    adh_neigh_lists = []
    for neighs, Jneighs, Jmed in zip(celldf["neighbour_list"],
                                     celldf["Jneighbour_list"],
                                     celldf["Jmed"]):
        ad_neigh_list = []
        for neigh, neighJ in zip(np.fromstring(neighs, sep=' ', dtype=int),
                                 np.fromstring(Jneighs, sep=' ')):
            if neigh != 0 and Jmed - neighJ / 2 > 0:
                ad_neigh_list.append(neigh)
        adh_neigh_lists.append(ad_neigh_list)
    return adh_neigh_lists


def get_adhering_clusters(celldf: pd.DataFrame) -> CellCluster:
    """Returns a list of list representing the different clusters of adhering cells in the
    system."""

    adh_neigh_lists = get_adhering_neighbours(celldf)
    clusters = []
    neigh_set = set(celldf["sigma"])
    neigh_dict = dict(zip(celldf["sigma"], adh_neigh_lists))
    while neigh_set:
        cluster = set()
        qneighs = {neigh_set.pop()}
        while qneighs:
            neigh = qneighs.pop()
            # Prevent infinite recursion
            if neigh not in cluster:
                cluster.add(neigh)
                qneighs.update(neigh_dict[neigh])

        neigh_set -= cluster
        clusters.append(set(cluster))

    return clusters


def get_median_gamma(key_locks, Jmed, Ja, weights, n=None):
    """Calculates the medium gamma from all x all cells.

    :param key_locks: Numpy matrix of (key, lock) pairs or each cell
    :param n: Samples n key-lock pairs from total population to speed up calculations
    :param Ja:
    :param Jmed:
    :param weights:
    """
    if n is not None:
        rng = np.random.default_rng(0)
        indexes = rng.choice(len(key_locks), min(n, len(key_locks)), replace=False)
        key_locks = key_locks[indexes]

    gammas = []
    for k1, l1 in key_locks:
        for k2, l2 in key_locks:
            gammas.append(get_gamma(k1, l1, k2, l2, Jmed, Ja, weights))

    return np.median(gammas)


def get_cell_cell_j(k1, l1, k2, l2, Ja, weights):
    return Ja + np.sum((k1 == l2) * weights) + np.sum((k2 == l1) * weights)


def get_gamma(k1, l1, k2, l2, Jmed, Ja, weights):
    return Jmed - get_cell_cell_j(k1, l1, k2, l2, Ja, weights) / 2


if __name__ == "__main__":
    main()
