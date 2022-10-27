import re
import sys
import pickle
import logging
from typing import Set, Dict, Tuple, List
from pathlib import Path
from dataclasses import dataclass


def main(netpath, outpath):
    logging.info("Reading networks")
    networks = read_networks(netpath)
    logging.info(f"Writing networks to: '{outpath}'")
    with open(outpath, 'wb') as file:
        pickle.dump(networks, file)
    logging.info("Finished")


@dataclass
class InNode:
    __slots__ = ("gene_id", "scaling")
    gene_id: int
    scaling: float

    def __hash__(self):
        return self.gene_id


@dataclass
class RegNode:
    __slots__ = ("gene_id", "threshold")
    gene_id: int
    threshold: float

    def __hash__(self):
        return self.gene_id


@dataclass
class OutNode:
    __slots__ = ("gene_id", "threshold")
    gene_id: int
    threshold: float

    def __hash__(self):
        return self.gene_id


@dataclass
class Cell:
    __slots__ = ("cell_sigma", "in_nr", "reg_nr", "out_nr", "id_map", "network")
    cell_sigma: int
    in_nr: int
    reg_nr: int
    out_nr: int
    id_map: dict  # Points to Nodes
    network: Dict[int, Set[Tuple[int, float]]]

    def __hash__(self):
        return self.cell_sigma


def read_networks(path, season_filter: List[int] = None):
    path = Path(path)
    seasons = {}
    for filepath in path.glob("t*.txt"):
        season_i = int(re.search(r"(?<=t)\d+", filepath.name).group())
        if season_filter is not None and season_i not in season_filter:
            continue
        cell_sigma = int(re.search(r"(?<=c)\d+", filepath.name).group())
        if cell_sigma == 0:
            continue

        with open(filepath) as file:
            raw = file.read()
        lines = [line for line in raw.split("\n") if line]
        node_nrs = [int(nr) for nr in lines[0].split(" ") if nr]
        scalings = [float(s) for s in lines[1].split(" ") if s]
        if len(scalings) != node_nrs[0]:
            raise ValueError("parsing scalings failed")

        id_map = {i: InNode(i, scalings[i]) for i in range(node_nrs[0])}
        network = {}
        for line in lines[2:]:
            gene_info = [info for info in line.split(" ") if info]
            if gene_info[0] == "1":
                gene_type = RegNode
            elif gene_info[0] == "2":
                gene_type = OutNode
            else:
                raise ValueError("unidentified gene type")
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
            cell_sigma,
            node_nrs[0],
            node_nrs[1],
            node_nrs[2],
            id_map,
            network
        )
        season = seasons.setdefault(season_i, {})
        season[cell_sigma] = cell
    return seasons


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    netpath = sys.argv[1]
    outpath = sys.argv[2]

    main(netpath, outpath)