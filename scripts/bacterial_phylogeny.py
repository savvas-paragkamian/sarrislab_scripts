#!/usr/bin/env python3

###############################################################################
# script name: bacterial_phylogeny.py
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:

###############################################################################
# usage:./bacterial_phylogeny.py
###############################################################################
import os,sys
import dendropy

print("Loading the tree")

path = "/home1/s.paragkamian/bacillus/543_de_novo/"
filename = "gtdbtk.bac120.decorated.tree"
tree = dendropy.Tree.get(
    path=os.path.join(path, filename),
    schema="newick")

print(type(tree))

print("Writing the tree")
output = "tree_ascii.txt"
tree_ascii = open(os.path.join(path, output),"w")

tree_ascii.write(tree.as_ascii_plot())
tree_ascii.close()
#sys.exit()
