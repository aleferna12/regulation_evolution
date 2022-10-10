import numpy as np
from ete3 import Tree, TreeStyle
from skbio import DistanceMatrix
from skbio.tree import nj
from io import StringIO
from parse_networks import read_networks, RegNode, InNode, OutNode

NETWORKS_PATH = "../runs/evolreg_false/networks"


def main():
    networks = read_networks(NETWORKS_PATH)
    # Build tree using cells from last season
    last_season = list(networks.values())[-1]
    tree = build_tree(last_season)
    tree_style = TreeStyle()
    tree_style.mode = 'c'
    tree.show(tree_style=tree_style)


def build_tree(season):
    # Not converting to index confuses nj builder
    ids = [str(i) for i in season]
    dists = dist_matrix(list(season.values()))
    dm = DistanceMatrix(dists, ids)
    sk_tree = nj(dm)
    newick = sk_tree.write(StringIO()).getvalue()
    return Tree(newick)


# This is a bit redundant as it mostly reconstructs the info as its saved in the file
# Still its good to be able to reconstruct this info from the network
def ordered_parameters(cell):
    weights = []
    thresholds = []
    scalings = []
    for gene_id in sorted(cell.network.keys()):
        for _, weigth in sorted(cell.network[gene_id]):
            weights.append(weigth)
    for gene_id, gene in sorted(cell.id_map.items()):
        if isinstance(gene, (RegNode, OutNode)):
            thresholds.append(gene.threshold)
        elif isinstance(gene, InNode):
            scalings.append(gene.scaling)
    return {
        "weights": weights,
        "thresholds": thresholds,
        "scalings": scalings
    }


def distance(pars1, pars2):
    if len(pars1) != len(pars2):
        raise ValueError("Parameters arrays have different lengths")

    dist = 0
    for x1, x2 in zip(pars1, pars2):
        dist += abs(x1 - x2)
    return dist


def dist_matrix(cells: list):
    array = np.zeros([len(cells), len(cells)])
    for i in range(len(cells)):
        for j in range(len(cells)):
            # Right now we only look at the weights so no transformation is necessary
            # to correct for differences in the variances of the parameter types
            array[i, j] = distance(
                ordered_parameters(cells[i])["weights"],
                ordered_parameters(cells[j])["weights"]
            )
            array[j, i] = array[i, j]
    return array


if __name__ == "__main__":
    main()