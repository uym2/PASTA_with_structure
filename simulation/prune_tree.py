#!/usr/bin/env python

from dendropy import Tree
from basic_utils import binary_search
from tree_lib import tree_as_newick
from sys import argv,exit
import os
import copy
import os.path

if __name__ == '__main__':
    if len(argv) < 3: 
        print("USAGE: treefile pruning_list_file [output]")
        exit(1)

    treefile = argv[1]
    listfile = argv[2]
    if len(argv) > 3:
	    outfile = argv[3]
    else:
	    outfile = None
    sample = sorted([ s.rstrip() for s in open(listfile).readlines() ])
    tree = Tree.get_from_path(treefile, 'newick',preserve_underscores=True)
    filt = lambda node: True if ((node.taxon is not None and not binary_search(sample,node.taxon.label)) or (node.label is not None and not binary_search(sample,node.label))) else False
    tree.filter_leaf_nodes(filt,recursive=True)
    #trees.write(file=open(sys.argv[-1],'w'), schema='newick', suppress_rooting=True)
    tree_as_newick(tree,outfile=outfile)
