from typing import Set, Dict, Tuple
import re
from pathlib import Path
from dataclasses import dataclass


NETWORKS_PATH = "../runs/neigh_info/networkdir_"


def main():
    pass


@dataclass
class InNode:
    gene_id: int
    scaling: float

    def __hash__(self):
        return self.gene_id


@dataclass
class RegNode:
    gene_id: int
    threshold: float

    def __hash__(self):
        return self.gene_id


@dataclass
class OutNode:
    gene_id: int
    threshold: float

    def __hash__(self):
        return self.gene_id


@dataclass
class Cell:
    cell_id: int
    in_nr: int
    reg_nr: int
    out_nr: int
    id_map: dict  # Points to Nodes
    network: Dict[int, Set[Tuple[int, float]]]

    def __hash__(self):
        return self.cell_id


def read_networks(path):
    path = Path(path)
    networks = {}
    for filepath in path.glob("t*.txt"):
        cell_i = int(re.search(r"(?<=c)\d+", filepath.name).group())
        if cell_i == 0:
            continue
        season_i = int(re.search(r"(?<=t)\d+", filepath.name).group())

        with open(filepath) as file:
            raw = file.read()
        lines = [line for line in raw.split("\n") if line]
        node_nrs = [int(nr) for nr in lines[0].split(" ") if nr]
        scalings = [float(s) for s in lines[1].split(" ") if s]
        if len(scalings) != node_nrs[0]:
            raise ValueError("Parsing scalings failed")

        id_map = {i: InNode(i, scalings[i]) for i in range(node_nrs[0])}
        network = {}
        for line in lines[2:]:
            gene_info = [info for info in line.split(" ") if info]
            if gene_info[0] == "1":
                gene_type = RegNode
            elif gene_info[0] == "2":
                gene_type = OutNode
            else:
                raise ValueError("Unidentified gene type")
            gene = gene_type(
                gene_id=int(gene_info[1]),
                threshold=float(gene_info[2])
            )
            id_map[gene.gene_id] = gene
            network[gene.gene_id] = set()
            weights = [float(w) for w in gene_info[3:]]
            for i, w in enumerate(weights):
                if gene_type == RegNode:
                    network[gene.gene_id].add((i, w))
                else:
                    network[gene.gene_id].add((i + node_nrs[0], w))

        cell = Cell(
            cell_i,
            node_nrs[0],
            node_nrs[1],
            node_nrs[2],
            id_map,
            network
        )
        season = networks.setdefault(season_i, {})
        season[cell_i] = cell

    return networks


def exhaust_cluster(neigh, neights, cluster):
    if neigh not in cluster:
        cluster.add(neigh)
        for next_neigh in neights[neigh]:
            exhaust_cluster(next_neigh, neights, cluster)
    return cluster


def get_clusters(neighs):
    clusters = []
    neigh_list = list(neighs)
    for neigh in neigh_list:
        cluster = exhaust_cluster(neigh, neighs, set())
        for found_neigh in cluster:
            neigh_list.remove(found_neigh)
        clusters.append(cluster)

    return clusters


def read_neighbours(path):
    seasons = {}
    path = Path(path)
    for filepath in path.glob("neigh*.txt"):
        neighbours = {}
        with open(filepath) as file:
            lines = file.read().split("\n")

        for line in lines:
            if not line:
                continue
            info = line.split(" ")
            neighbours[int(info[0])] = [int(n) for n in info[1:]]

        get_clusters(neighbours)
        season_i = int(re.search(r"(?<=t)\d+", filepath.name).group())
        seasons[season_i] = neighbours
    return seasons


if __name__ == "__main__":
    main()