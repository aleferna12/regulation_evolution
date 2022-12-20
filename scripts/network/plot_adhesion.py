import logging
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from parse import parse_cell_data, build_time_filter, get_time_points


def main():
    logging.basicConfig(level=logging.INFO)

    datapath = sys.argv[1]
    outfile = sys.argv[2]
    start_season = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    n_seasons = int(sys.argv[4]) if len(sys.argv) > 4 else None
    n_samples = int(sys.argv[5]) if len(sys.argv) > 5 else None

    t_filter = build_time_filter(get_time_points(datapath), start=start_season, n=n_seasons)
    cell_df = parse_cell_data(datapath, time_filter=t_filter)
    fig = make_plots(cell_df, n_samples)

    if outfile[-5:] == ".html":
        fig.write_html(outfile)
    else:
        fig.write_image(outfile)


def make_plots(cell_df, n_samples):
    fig = make_subplots(3, 1, subplot_titles=["Prevalence of multicellularity",
                                              "Medium gamma of neighbouring cells",
                                              "Medium gamma of all x all cells"])

    x = []
    prev_y = []
    neigh_gammas = []
    all_gammas = []
    for key in ["jkey", "jlock"]:
        cell_df[key] = cell_df[key].apply(lambda string: np.array(list(string), dtype=int))
    cell_df["jkey_jlock"] = list(np.stack([cell_df["jkey"], cell_df["jlock"]], axis=1))

    for season in set(cell_df["time"]):
        x.append(int(season))
        sdf = cell_df[cell_df["time"] == season]

        clusters = get_adhering_clusters(
            sdf["sigma"],
            get_adhering_neighbours(
                sdf["neighbour_list"],
                sdf["neighbourJ_list"],
                sdf["medJ"]
            )
        )
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
            sdf["jkey_jlock"],
            4,
            8,
            np.array([4, 3, 2, 1, 1, 1]),
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


def get_adhering_neighbours(neigh_lists, neighJ_lists, medJs):
    """Return a list of lists with the neighbours that each cell is adhering to."""
    adh_neigh_lists = []
    for neigh_list, neighJ_list, medJ in zip(neigh_lists, neighJ_lists, medJs):
        ad_neigh_list = []
        for neigh, neighJ in zip(neigh_list, neighJ_list):
            if neigh != 0 and medJ - neighJ / 2 > 0:
                ad_neigh_list.append(neigh)
        adh_neigh_lists.append(ad_neigh_list)
    return adh_neigh_lists


def get_adhering_clusters(sigmas, adh_neigh_lists):
    """Returns a list of list representing the different clusters of adhering cells in the
    system."""
    neigh_set = set(sigmas)
    if len(neigh_set) != len(sigmas):
        raise ValueError("make sure this function is called for a single time point")
    clusters = []
    neigh_dict = dict(zip(sigmas, adh_neigh_lists))
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
        clusters.append(list(cluster))

    return clusters


def get_medium_gamma(key_locks, ja, jam, f_arr, n=None):
    """Calculates the medium gamma from all x all cells.

    :param key_locks: Numpy matrix of (key, lock) pairs or each cell
    :param ja:
    :param jam:
    :param f_arr:
    :param n: Samples n key-lock pairs from total population to speed up calculations
    :return:
    """
    if n is not None:
        rng = np.random.default_rng(0)
        indexes = rng.choice(len(key_locks), n, replace=False)
        key_locks = key_locks[indexes]

    gammas = []
    for k1, l1 in key_locks:
        for k2, l2 in key_locks:
            gammas.append(get_gamma(ja, jam, f_arr, l1, k1, l2, k2))

    return np.mean(gammas)


def get_cell_cell_j(ja, l1, k1, l2, k2):
    v = len(k1)
    return ja + 2 * v - np.count_nonzero(l1 != k2) - np.count_nonzero(l2 != k1)


def get_cell_med_j(jam, f_arr, key):
    return jam + np.sum(f_arr * key[:len(f_arr)])


def get_gamma(ja, jam, f_seq, l1, k1, l2, k2):
    return get_cell_med_j(jam, f_seq, k1) - get_cell_cell_j(ja, l1, k1, l2, k2) / 2


if __name__ == "__main__":
    main()
