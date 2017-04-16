from smpl_merger import smplMerger
from sys import argv,stdout
from sequence_lib import read_fasta, write_fasta,gap_propagate
from math import log

aln1_file = argv[1]
aln2_file = argv[2]

if len(argv) > 3:
	outfile = argv[3]
else:	
	outfile = None

#taxon_name1, aln1 = read_fasta(aln1_file)
#taxon_name2, aln2 = read_fasta(aln2_file)

MGR = smplMerger(aln1_file,aln2_file)
score,cons1,cons2 = MGR.smpl_merge(nsmpl=10,n1=20,n2=20)

if outfile:
	fout = open(outfile,'w')
else:
	fout = stdout

for i in range(len(MGR.aln1)):
	fout.write(">"+MGR.name1[i]+"\n")
	fout.write(gap_propagate(cons1,MGR.aln1[i])+"\n")
for i in range(len(MGR.aln2)):
	fout.write(">"+MGR.name2[i]+"\n")
	fout.write(gap_propagate(cons2,MGR.aln2[i])+"\n")

if outfile:
	fout.close()
