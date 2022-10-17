import re
import sys
import logging
import json
from pathlib import Path
from typing import List


def main():
    logging.basicConfig(level=logging.INFO)
    neighpath = sys.argv[1]
    outpath = sys.argv[2]
    logging.info("Reading neighbour information")
    neighs = read_neighbours(neighpath)
    logging.info(f"Writing neighbours to: '{outpath}'")
    with open(outpath, 'w') as file:
        json.dump(neighs, file, indent=4)
    logging.info("Finished")


def exhaust_cluster(neigh, neights, cluster):
    if neigh not in cluster:
        cluster.add(neigh)
        for next_neigh in neights[neigh]:
            exhaust_cluster(next_neigh, neights, cluster)
    return cluster


def get_clusters(neighs):
    clusters = []
    neigh_list = list(neighs)
    while neigh_list:
        cluster = list(exhaust_cluster(neigh_list[0], neighs, set()))
        for found_neigh in cluster:
            if found_neigh in neigh_list:
                neigh_list.remove(found_neigh)
        clusters.append(cluster)

    return clusters


def read_neighbours(path, season_filter: List[int] = None):
    seasons = {}
    path = Path(path)
    for filepath in path.glob("neigh*.txt"):
        season_i = int(re.search(r"(?<=t)\d+", filepath.name).group())
        if season_filter is not None and season_i not in season_filter:
            continue

        neighbours = {}
        with open(filepath) as file:
            lines = file.read().split("\n")
        for line in lines:
            if not line:
                continue
            info = line.split(" ")
            if info[0] == '0':
                continue
            neighbours[int(info[0])] = [int(n) for n in info[1:]]

        clusters = get_clusters(neighbours)
        seasons[season_i] = clusters
    return seasons


if __name__ == "__main__":
    main()