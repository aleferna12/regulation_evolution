"""Makes newick trees from cell ancestry data.

Be mindful of the fact that not all generations are captured, since data is only saved every X
MCSs. This means that some branches will be collapsed into multifurcations, resulting in some loss
of phylogenetic information.
"""
import logging
import sys
import pandas as pd
from pathlib import Path
from ete3 import Tree
from parse import parse_cell_data


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

    celldf = parse_cell_data(datadir)
    trees = make_trees(celldf, dead_ends, names, single_leaf)
    longest_time, longest = write_trees(trees, outdir, nhx, fmt)

    logging.info("Longest lineages are: " + ", ".join(longest))
    logging.info(f"These lineages survived for {longest_time} MCSs)")
    logging.info("Finished")


def write_trees(trees, outdir, nhx=True, fmt=5):
    logging.info(f"Writing trees to '{outdir}'")
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
    logging.info("Building trees from cell ancestry data")
    # This dict holds the descendents of the current cells
    prev_anc_child = {}
    # Order by time backwards
    # This is needed to deal with dead-ends in a nice way
    for time in celldf["time"].sort_values(ascending=False).unique():
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
