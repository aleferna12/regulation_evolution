import logging
import sys

import numpy as np
import plotly.graph_objects as go
from typing import List, Set
from plotly.subplots import make_subplots
from enlighten import Counter
from parse import *

# Type alias
CellCluster = List[Set[int]]
logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    datadir = sys.argv[1]
    outfile = sys.argv[2]
    start_time_step = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    n_time_steps = int(sys.argv[4]) if len(sys.argv) > 4 else None

    t_filter = build_time_filter(get_time_points(datadir), start=start_time_step, n=n_time_steps)
    celldf = parse_cell_data(datadir, time_filter=t_filter)
    fig = make_plots(celldf)

    logger.info(f"Writing files to: {outfile}")
    if outfile[-5:] == ".html":
        fig.write_html(outfile)
    else:
        fig.write_image(outfile)

    logger.info("Finished")


def make_plots(celldf: pd.DataFrame):
    logger.info("Making adhesion plots")

    fig = make_subplots(5, 1, subplot_titles=[
        "Prevalence of multicellularity",
        "Median gamma of neighbouring cells",
        "Median gamma of neighbouring migrating cells",
        "Median gamma of neighbouring dividing cells",
        "Median gamma of neighbouring migrating-dividing cells"
    ])

    x = []
    prev_y = []
    neigh_gammas_med = []
    mig_mig_gammas_med = []
    div_div_gammas_med = []
    mig_div_gammas_med = []
    
    time_steps = celldf["time"].unique()
    pbar = Counter(total=len(time_steps), desc='Time-steps')
    for season in time_steps:
        x.append(int(season))
        sdf = celldf[celldf["time"] == season]

        clusters = get_adhering_clusters(sdf)
        unicells = 0
        total_pop = 0
        for c in clusters:
            if len(c) == 1:
                unicells += 1
            total_pop += len(c)
        prev_y.append(1 - unicells / total_pop)

        neigh_gammas = []
        mig_mig_gammas = []
        div_div_gammas = []
        mig_div_gammas = []
        for sigma, tau, Jmed, neighs, neighJs in zip(sdf["sigma"],
                                                     sdf["tau"],
                                                     sdf["Jmed"],
                                                     sdf["neighbour_list"],
                                                     sdf["Jneighbour_list"]):
            for neigh, neighJ in zip(neighs, neighJs):
                if neigh == 0:
                    continue
                gamma = Jmed - neighJ / 2
                if tau == sdf.loc[(season, neigh), "tau"]:
                    if tau == 1:
                        mig_mig_gammas.append(gamma)
                    else:
                        div_div_gammas.append(gamma)
                else:
                    mig_div_gammas.append(gamma)
                neigh_gammas.append(gamma)

        neigh_gammas_med.append(np.median(neigh_gammas))
        mig_mig_gammas_med.append(np.median(mig_mig_gammas))
        div_div_gammas_med.append(np.median(div_div_gammas))
        mig_div_gammas_med.append(np.median(mig_div_gammas))

        pbar.update()

    for i, ys in enumerate(
            [prev_y, neigh_gammas_med, mig_mig_gammas_med, div_div_gammas_med, mig_div_gammas_med],
            1):
        fig.add_trace(go.Scatter(
            x=x,
            y=ys,
            mode="lines"
        ), row=i, col=1)
    fig.update_yaxes(title="gamma")
    fig.update_yaxes(title="% of adhering cells", range=[0, 1], row=1, col=1)
    return fig.update_layout(showlegend=False)


def get_adhering_neighbours(celldf: pd.DataFrame):
    """Return a list of lists with the neighbours that each cell is adhering to."""
    for key in ["neighbour_list", "Jneighbour_list"]:
        if isinstance(celldf[key].iat[0], str):
            celldf[key] = celldf[key].apply(lambda s: np.fromstring(s, sep=' ', dtype=int))

    adh_neigh_lists = []
    for neigh_list, neighJ_list, Jmed in zip(celldf["neighbour_list"],
                                             celldf["Jneighbour_list"],
                                             celldf["Jmed"]):
        ad_neigh_list = []
        for neigh, neighJ in zip(neigh_list, neighJ_list):
            if neigh != 0 and Jmed - neighJ / 2 > 0:
                ad_neigh_list.append(neigh)
        adh_neigh_lists.append(ad_neigh_list)
    return adh_neigh_lists


def get_adhering_clusters(celldf: pd.DataFrame) -> CellCluster:
    """Returns a list of list representing the different clusters of adhering cells in the
    system."""
    if not celldf["time"].nunique() == 1:
        raise ValueError("make sure this function is called for a single time point")

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
