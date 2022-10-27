import re
import sys
import logging
import json
from pathlib import Path
from typing import List


def main(neighpath, outpath):
    logging.info("Reading neighbour information")
    neighs = read_neighbours(neighpath)
    logging.info(f"Writing neighbours to: '{outpath}'")
    with open(outpath, 'w') as file:
        json.dump(neighs, file, indent=4)
    logging.info("Finished")


def get_clusters(neighs):
    clusters = []
    neigh_set = set(neighs)
    while neigh_set:
        cluster = set()
        qneighs = {neigh_set.pop()}
        while qneighs:
            neigh = qneighs.pop()
            # Prevent infitine recursion
            if neigh not in cluster:
                cluster.add(neigh)
                qneighs.update(neighs[neigh])

        neigh_set -= cluster
        # Sets are not JSON hashable
        clusters.append(list(cluster))

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
    logging.basicConfig(level=logging.INFO)
    neighpath = sys.argv[1]
    outpath = sys.argv[2]

    main(neighpath, outpath)