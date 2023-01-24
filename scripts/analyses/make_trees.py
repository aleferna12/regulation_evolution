"""Makes newick trees from cell ancestry data.

Be mindful of the fact that not all generations are captured, since data is only saved every X
MCSs. This means that some branches will be collapsed into multifurcations, resulting in some loss
of phylogenetic information. This can be improved (with some amount of work) by saving a vector of
ancestors rather than a single one.
"""
import logging
import sys
import pandas as pd
from pathlib import Path
from enlighten import Counter
from ete3 import Tree
from parse import *

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    datadir = sys.argv[1]
    outdir = Path(sys.argv[2])
    if not outdir.is_dir():
        raise ValueError("second argument is not a valid existing directory")
    dead_ends = False if len(sys.argv) > 3 and sys.argv[3] in ["0", "false"] else True
    names = False if len(sys.argv) > 4 and sys.argv[4] in ["0", "false"] else True
    single_leaf = False if len(sys.argv) > 5 and sys.argv[5] in ["0", "false"] else True
    nhx = False if len(sys.argv) > 6 and sys.argv[6] in ["0", "false"] else True
    fmt = int(sys.argv[7]) if len(sys.argv) > 7 else 5

    t_filter = build_time_filter(get_time_points(datadir))
    celldf = parse_cell_data(datadir, time_filter=t_filter)
    trees = make_trees(celldf, dead_ends, names, single_leaf)
    longest_time, longest = write_trees(trees, outdir, nhx, fmt)

    logger.info("Longest lineages are: " + ", ".join(longest))
    logger.info(f"These lineages survived for {longest_time} MCSs)")
    logger.info("Finished")


def write_trees(trees, outdir, nhx=True, fmt=5):
    logger.info(f"Writing trees to '{outdir}'")
    outdir = Path(outdir)

    longest = []
    longest_time = 0
    for i, tree in enumerate(trees):
        tree_name = "tree" + str(i)
        time = tree.get_farthest_leaf()[0].time
        if time > longest_time:
            longest = [tree_name]
            longest_time = time
        elif time == longest_time:
            longest.append(tree_name)
        tree.write(
            outfile=str(outdir / (tree_name + ".newick")),
            format=fmt,
            features=["time"] if nhx else None
        )

    return longest_time, longest


# Remember that although the sigmas are added as names, the same id can refer to DIFFERENT cells
def make_trees(celldf: pd.DataFrame, dead_ends=True, names=True, single_leaf=True):
    logger.info("Building trees from cell ancestry data")
    # This dict holds the descendents of the current cells
    prev_anc_child = {}
    # Order by time backwards
    # This is needed to deal with dead-ends in a nice way
    timesteps = celldf["time"].sort_values(ascending=False).unique()
    pbar = Counter(total=len(timesteps), desc="Time-steps")
    for time in timesteps:
        # This dict holds the ancestors of the current cells
        next_anc_child = {}
        tdf = celldf.loc[time]
        for sigma, anc_sigma in zip(tdf["sigma"], tdf["ancestor"]):
            # Keep dead ends from being added to the trees
            if not dead_ends and prev_anc_child and sigma not in prev_anc_child:
                continue
            node = Tree(name=str(sigma) if names else "")
            node.add_feature("time", time)
            if sigma in prev_anc_child:
                for child in prev_anc_child[sigma]:
                    node.add_child(child)
            if anc_sigma not in next_anc_child:
                next_anc_child[anc_sigma] = []
            next_anc_child[anc_sigma].append(node)
        prev_anc_child = dict(next_anc_child)

        pbar.update()

    # Construct root trees
    trees = []
    for anc_sigma, children in prev_anc_child.items():
        # The roots date from the start of the simulation
        # Unfortunately, root attributes can't be saved in NHX, so adding it as a feature is
        # futile
        root = Tree(name=anc_sigma if names else "")
        for child in children:
            root.add_child(child)
        # Exclude trees that only have a single leaf due to being part of an exclusively
        # migrating lineage or because the lineage only has one survivor and dead_ends=False
        # These single-leaf lineages can cause problems in some newick parsers
        if single_leaf or len(root.get_leaves()) > 1:
            trees.append(root)
    return trees


if __name__ == "__main__":
    main()
