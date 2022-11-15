import os.path
import sys
import pickle
import logging
import numpy as np
from pathlib import Path
# Need to import Cell for pickle deserialization
from parse_networks import Cell, RegNode, InNode, OutNode, read_networks


def main():
    logging.basicConfig(level=logging.INFO)
    season = int(sys.argv[1])
    netpath = Path(sys.argv[2]).resolve()
    outpath = Path(sys.argv[3]).resolve()

    construct_dist_matrix(season, netpath, outpath)


def construct_dist_matrix(season, netpath, outpath):
    logging.info("Reading network file")
    if os.path.isfile(netpath):
        with open(netpath, "rb") as file:
            networks = pickle.load(file)
    else:
        logging.warning("Network file not found")
        logging.info(f"Try to build networks from '{netpath}'")
        networks = read_networks(netpath, season_filter=[season])

    logging.info(f"This season has {len(networks[season])} cells")
    logging.info("Retrieving regulatory parameters from networks")
    param_season = construct_param_season(networks[season])

    logging.info("Creating distance matrix")
    dm, ids = dist_matrix(param_season)
    logging.info(f"Writing distance matrix to '{outpath}'")
    write_dist_matrix(dm, ids, outpath)
    logging.info("Finished")


def write_dist_matrix(dm, ids, outpath):
    # Convert floats to string
    dm = dm.astype("<U16")
    with open(outpath, "w") as file:
        file.write(str(len(ids)) + '\n')
        for row, cell_sigma in zip(dm, ids):
            file.write(str(cell_sigma) + '\t')
            file.write('\t'.join(row) + '\n')


def construct_param_season(season):
    param_dict = {}
    for cell_sigma, cell in season.items():
        param_dict[cell_sigma] = ordered_parameters(cell)
    return param_dict


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
        "weights": np.array(weights),
        "thresholds": np.array(thresholds),
        "scalings": np.array(scalings)
    }


def dist_matrix(season):
    # Not converting to index confuses nj builder
    array = np.zeros([len(season), len(season)])
    for i, cell_sigmai in enumerate(season):
        for j, cell_sigmaj in enumerate(season):
            # Right now we only look at the weights so no transformation is necessary
            # to correct for differences in the variances of the parameter types
            array[i, j] = np.abs(
                season[cell_sigmai]["weights"] - season[cell_sigmaj]["weights"]
            ).sum()
            array[j, i] = array[i, j]
    return array, [str(cell_sigma) for cell_sigma in season]


if __name__ == "__main__":
    main()