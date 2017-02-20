#! /usr/bin/env python

from sys import argv

def gap_propagate(cons_seq,targ_seq):
# propagate gaps from cons_seq to targ_seq; this is used after the cons_seq was aligned with another sequence and we want to
# propage the alignment to targ_seq. This is a similar idea with transitivity used in PASTA; the ultimate goal is to merge 2 alignments,
# but in this case we have a consensus sequence for each alignments and we also have a good way (properly using seconday structure) to align them.
# NOTE: gap_propagate is NOT symmetric: only propagate from consensus to target, not in reverse; careful consider which sequence is the cons_seq and which is targ_seq!

	out_seq = ''
	i = 0
	for c in cons_seq:
		if c != '-':
			out_seq += targ_seq[i]
			i += 1
		else:
			out_seq += '-'

	return out_seq

# open input files
cons_seq = argv[1]
targ_aln = argv[2]

with open(cons_seq,'r') as f1:
	cons_seq = f1.readline().rstrip()
	with open(targ_aln,'r') as f:
			for line in f:
				if line[0] == '>':
					print(line.rstrip())
				else:
					targ_seq = line.rstrip()
					print(gap_propagate(cons_seq,targ_seq))
