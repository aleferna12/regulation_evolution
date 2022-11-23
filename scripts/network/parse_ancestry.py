"""Parse ancestry files into newick trees.

Be mindful of the fact that not all generations are captured, since data is only saved every X
MCSs. This means that some branches will be collapsed into multifurcations, resulting in some loss
of phylogenetic information.
"""
import re
import logging
import sys
import os
from pathlib import Path
from ete3 import Tree


def main():
    logging.basicConfig(level=logging.INFO)
    netpath = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve()
    if not os.path.isdir(outdir):
        raise ValueError("second argument is not a valid existing directory")
    dead_ends = False if len(sys.argv) > 3 and sys.argv[3] in ["0", "false"] else True
    names = False if len(sys.argv) > 4 and sys.argv[4] in ["0", "false"] else True
    single_leaf = False if len(sys.argv) > 5 and sys.argv[5] in ["0", "false"] else True
    nhx = False if len(sys.argv) > 6 and sys.argv[6] in ["0", "false"] else True
    fmt = int(sys.argv[7]) if len(sys.argv) > 7 else 5

    parse_ancestry(netpath, outdir, dead_ends, names, single_leaf, nhx, fmt)


def parse_ancestry(netpath, outdir, dead_ends, names, single_leaf, nhx, fmt):
    logging.info(f"Writing trees to '{outdir}'")
    longest = []
    longest_time = 0
    for i, tree in enumerate(read_ancestry(netpath, dead_ends, names, single_leaf)):
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
    logging.info("Longest lineages are: " + ", ".join(longest))
    logging.info(f"These lineages survived for {longest_time} MCSs)")
    logging.info("Finished")

    return longest


def read_ancestry(path, dead_ends=True, names=True, single_leaf=True, up_to=-1):
    path = Path(path)
    timepoints = {}
    for filepath in path.glob("anc_*.txt"):
        time = int(re.search(r"(?<=t)\d+", filepath.name).group())
        if up_to != -1 and time > up_to:
            continue
            
        cells = []
        with open(filepath) as file:
            lines = file.read().split("\n")
        for line in lines:
            if not line:
                continue
            cell_sigma, anc_sigma = line.split(" ")
            if cell_sigma == "0":
                continue
            cells.append((int(cell_sigma), int(anc_sigma)))
        timepoints[time] = cells
    return make_trees(timepoints, dead_ends, names, single_leaf)


# Remember that although the sigmas are added as names, the same id can refer to DIFFERENT cells
def make_trees(timepoints, dead_ends=True, names=True, single_leaf=True):
    # Order by time backwards
    # This is needed to deal with dead-ends
    timepoints = {k: timepoints[k] for k in sorted(timepoints, reverse=True)}
    # This dict holds the descendents of the current cells
    prev_anc_child = {}
    for time, cell_batch in timepoints.items():
        # This dict holds the ancestors of the current cells
        next_anc_child = {}
        for cell_sigma, anc_sigma in cell_batch:
            # Keep dead ends from being added to the trees
            if not dead_ends and prev_anc_child and cell_sigma not in prev_anc_child:
                continue
            node = Tree(name=str(cell_sigma) if names else "")
            node.add_feature("time", time)
            if cell_sigma in prev_anc_child:
                for child in prev_anc_child[cell_sigma]:
                    node.add_child(child)
            if anc_sigma not in next_anc_child:
                next_anc_child[anc_sigma] = []
            next_anc_child[anc_sigma].append(node)
        prev_anc_child = dict(next_anc_child)

    # Construct root trees
    trees = []
    for anc_sigma, children in prev_anc_child.items():
        root = Tree(name=anc_sigma if names else "")
        # The roots date from the start of the simulation
        # Unfortunately, root attributes can't be saved in NHX, so adding it as a feature is
        # futile
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