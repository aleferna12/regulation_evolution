"""Makes newick trees from cell ancestry data.

Be mindful of the fact that not all generations are captured, since data is only saved every X
MCSs. This means that some branches will be collapsed into multifurcations, resulting in some loss
of phylogenetic information. This can be improved (with some amount of work) by saving a vector of
ancestors rather than a single one.
"""
import argparse
import logging
import sys
import pandas as pd
from pathlib import Path
from enlighten import Counter
from ete3 import Tree
from fileio import *

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="make_trees",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("datadir", help="directory containing the cell data as CSV files")
    parser.add_argument("outdir", help="directory where the trees will be created")
    parser.add_argument("-e", 
                        "--extinct", 
                        action="store_true",
                        help="include extinct lineages")
    parser.add_argument("-n", 
                        "--no-names",
                        action="store_true",
                        help="do not include the sigmas of cells as branch names")
    parser.add_argument(
        "-s",
        "--single-lineage",
        action="store_true",
        help="allow trees that consist of a single lineage to be created (may cause problems with "
             "visualization in some programs)"
    )
    parser.add_argument("-x",
                        "--no-nhx",
                        action="store_true",
                        help="do not include NHX attributes in the trees (such as time)")
    parser.add_argument("-m",
                        "--stop-mrca",
                        action="store_true",
                        help="stop computing when the MRCA of the living cells is found")
    parser.add_argument(
        "-c",
        "--no-collapse",
        action="store_true",
        help="do not collapse lineages with a single descendent into a single branch"
    )
    parser.add_argument("-r",
                        "--no-root",
                        action="store_true",
                        help="do not include root attributes and information for compatibility")
    parser.add_argument("-f",
                        "--format",
                        default=5,
                        type=int,
                        help="tree output format (following etetoolkit conventions)")
    parser.add_argument("-t",
                        "--last-time",
                        default=float("inf"),
                        type=int,
                        help="last time point to look at when making the trees (i.e. the present)")
    args = parser.parse_args()

    t_filter = build_time_filter(get_time_points(args.datadir), end=args.last_time)
    celldf = parse_cell_data(args.datadir, time_filter=t_filter)
    trees = make_trees(celldf,
                       extinct=args.extinct,
                       names=not args.no_names,
                       single_lineage=args.single_lineage,
                       stop_mrca=args.stop_mrca,
                       collapse_branches=not args.no_collapse)
    longest_time, longest = write_trees(trees,
                                        args.outdir,
                                        not args.no_nhx,
                                        args.format,
                                        not args.no_root)

    logger.info("Longest lineages are: " + ", ".join(longest))
    logger.info(f"These lineages survived for {longest_time} MCSs)")
    logger.info("Finished")


def write_trees(trees, outdir, nhx=True, fmt=5, root_info=True):
    logger.info(f"Writing trees to '{outdir}'")
    outdir = Path(outdir)

    longest = []
    longest_time = 0
    for i, tree in enumerate(trees):
        tree_name = "tree" + str(i)
        time = tree.get_farthest_leaf()[0].time - tree.time
        if time > longest_time:
            longest = [tree_name]
            longest_time = time
        elif time == longest_time:
            longest.append(tree_name)
        tree.write(
            outfile=str(outdir / (tree_name + ".newick")),
            format=fmt,
            features=["time"] if nhx else None,
            format_root_node=root_info
        )

    return longest_time, longest


# Remember that although the sigmas are added as names, the same id can refer to DIFFERENT cells
def make_trees(celldf: pd.DataFrame,
               extinct=False,
               names=True,
               single_lineage=True,
               stop_mrca=False,
               collapse_branches=True):
    logger.info("Building trees from cell ancestry data")
    # This dict holds the descendents of the current cells
    prev_anc_child = {}
    # Order by time backwards
    # This is needed to deal with dead-ends in a nice way and to terminate early if MRCA is found
    timesteps = celldf["time"].sort_values(ascending=False).unique()
    pbar = Counter(total=len(timesteps), desc="Time-steps")
    time = -1  # Just in case there are no timesteps
    for time in timesteps:
        if stop_mrca and len(prev_anc_child) == 1:
            pbar.close()
            logger.info("Found MRCA, stopping analysis")
            break
        # This dict holds the ancestors of the current cells
        next_anc_child = {}
        tdf = celldf.loc[time]
        for sigma, anc_sigma in zip(tdf["sigma"], tdf["ancestor"]):
            # Keep dead ends from being added to the trees
            if not extinct and prev_anc_child and sigma not in prev_anc_child:
                continue

            node = _build_node(
                sigma,
                prev_anc_child.get(sigma, []),
                time,
                names,
                collapse_branches
            )
            if anc_sigma not in next_anc_child:
                next_anc_child[anc_sigma] = []
            next_anc_child[anc_sigma].append(node)
        prev_anc_child = dict(next_anc_child)
        pbar.update()

    # Construct root trees
    trees = []
    for anc_sigma, children in prev_anc_child.items():
        root = _build_node(anc_sigma, children, time, names, collapse_branches)
        # Exclude trees that only have a single leaf due to being part of an exclusively
        # migrating lineage or because the lineage only has one survivor and extinct=False
        # These single-leaf lineages can cause problems in some newick parsers
        if single_lineage or len(root.get_leaves()) > 1:
            trees.append(root)
    return trees


def _build_node(sigma, children, time, names, collapse_branches):
    if collapse_branches and len(children) == 1:
        children[0].dist += 1
        return children[0]
    node = Tree(name=str(sigma) if names else "")
    node.add_feature("time", time)
    for child in children:
        node.add_child(child)
    return node


if __name__ == "__main__":
    main()
