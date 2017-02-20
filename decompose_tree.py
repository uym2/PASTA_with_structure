#! /usr/bin/env python

from dendropy import Tree
from sepp_tool.tree import PhylogeneticTree
import sys
from basic_utils import  trl_zero
import os

path = sys.argv[1]
out_path = sys.argv[2]
M = int(sys.argv[3])
if len(sys.argv) > 4:
	m = int(sys.argv[4])
else:
	m = None

t = Tree()
t.read_from_path(path,'newick')
T = PhylogeneticTree(t)

print('computing treeMap ... ')
treeMap = T.decompose_tree(M,'centroid',minSize=m,decomp_strategy='normal')

treeNum = len(treeMap.keys())
#print treeNum

print('writing tree ... ')

try:
    os.stat(out_path)
except:
    os.mkdir(out_path)

for (k,T) in treeMap.items():
	#print k
	T.write_newick_to_path( out_path + '/tree_' + trl_zero(k,len(str(treeNum))) + '.tre')
