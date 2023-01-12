import logging
import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from typing import List, Set
from plotly.subplots import make_subplots
from parse import parse_cell_data, build_time_filter, get_time_points

# Type alias
CellCluster = List[Set[int]]


def main():
    logging.basicConfig(level=logging.INFO)

    datadir = sys.argv[1]
    outfile = sys.argv[2]
    start_season = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    n_seasons = int(sys.argv[4]) if len(sys.argv) > 4 else None
    n_samples = int(sys.argv[5]) if len(sys.argv) > 5 else None

    t_filter = build_time_filter(get_time_points(datadir), start=start_season, n=n_seasons)
    celldf = parse_cell_data(datadir, time_filter=t_filter)
    fig = make_plots(celldf, n_samples)

    if outfile[-5:] == ".html":
        fig.write_html(outfile)
    else:
        fig.write_image(outfile)


def make_plots(celldf: pd.DataFrame, n_samples):
    fig = make_subplots(3, 1, subplot_titles=["Prevalence of multicellularity",
                                              "Medium gamma of neighbouring cells",
                                              "Medium gamma of all x all cells"])

    x = []
    prev_y = []
    neigh_gammas = []
    all_gammas = []

    weights = np.array([1, 2, 3, 4, 5, 6])

    def array_from_dec(j_dec):
        return np.array(list(np.binary_repr(j_dec, width=len(weights))), dtype=int)

    celldf["jkey"] = celldf["jkey_dec"].apply(array_from_dec)
    celldf["jlock"] = celldf["jlock_dec"].apply(array_from_dec)
    celldf["jkey_jlock"] = list(np.stack([celldf["jkey"], celldf["jlock"]], axis=1))

    for season in celldf["time"].unique():
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

        for sigma, medJ, neighs, neighJs in zip(sdf["sigma"],
                                                sdf["medJ"],
                                                sdf["neighbour_list"],
                                                sdf["neighbourJ_list"]):
            for neigh, neighJ in zip(neighs, neighJs):
                if neigh != 0:
                    neigh_gammas.append(medJ - neighJ / 2)

        all_gammas.append(get_medium_gamma(
            sdf["jkey_jlock"].to_numpy(),
            14,
            7,
            weights,
            n_samples
        ))

    fig.add_trace(go.Scatter(
        x=x,
        y=prev_y,
        mode="lines"
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=x,
        y=neigh_gammas,
        mode="lines"
    ), row=2, col=1)
    fig.add_trace(go.Scatter(
        x=x,
        y=all_gammas,
        mode="lines",
    ), row=3, col=1)
    fig.update_yaxes(title="% of adhering cells", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title="gamma", row=2, col=1)
    fig.update_yaxes(title="gamma", row=3, col=1)
    return fig.update_layout(showlegend=False)


def get_adhering_neighbours(celldf: pd.DataFrame):
    """Return a list of lists with the neighbours that each cell is adhering to."""
    adh_neigh_lists = []
    for neigh_list, neighJ_list, medJ in zip(celldf["neighbour_list"],
                                             celldf["neighbourJ_list"],
                                             celldf["medJ"]):
        ad_neigh_list = []
        for neigh, neighJ in zip(neigh_list, neighJ_list):
            if neigh != 0 and medJ - neighJ / 2 > 0:
                ad_neigh_list.append(neigh)
        adh_neigh_lists.append(ad_neigh_list)
    return adh_neigh_lists


def get_adhering_clusters(celldf: pd.DataFrame) -> CellCluster:
    """Returns a list of list representing the different clusters of adhering cells in the
    system."""
    if not celldf["sigma"].is_unique:
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


def get_medium_gamma(key_locks, jmed, ja, weights, n=None):
    """Calculates the medium gamma from all x all cells.

    :param key_locks: Numpy matrix of (key, lock) pairs or each cell
    :param n: Samples n key-lock pairs from total population to speed up calculations
    :param ja:
    :param jmed:
    :param weights:
    """
    if n is not None:
        rng = np.random.default_rng(0)
        indexes = rng.choice(len(key_locks), min(n, len(key_locks)), replace=False)
        key_locks = key_locks[indexes]

    gammas = []
    for k1, l1 in key_locks:
        for k2, l2 in key_locks:
            gammas.append(get_gamma(k1, l1, k2, l2, jmed, ja, weights))

    return np.mean(gammas)


def get_cell_cell_j(k1, l1, k2, l2, ja, weights):
    return ja + np.sum((k1 == l2) * weights) + np.sum((k2 == l1) * weights)


def get_gamma(k1, l1, k2, l2, jmed, ja, weights):
    return jmed - get_cell_cell_j(k1, l1, k2, l2, ja, weights) / 2


if __name__ == "__main__":
    main()
