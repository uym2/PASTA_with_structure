from merger_backup import merger
from sys import argv,stdout
from sequence_lib import read_fasta, write_fasta,gap_propagate
from rand_aln_prob import rand_Aln
from math import log

aln1_file = argv[1]
aln2_file = argv[2]
ref_aln_file = argv[3]
if len(argv) > 4:
	outfile = argv[4]
else:	
	outfile = None

taxon_name1, aln1 = read_fasta(aln1_file)
taxon_name2, aln2 = read_fasta(aln2_file)
#taxon_name_all, ref_aln = read_fasta(ref_aln_file)


rand_P = rand_Aln(max(len(aln1[0]),len(aln2[0])))
#print(rand_P.size())
#print(len(aln1[0]))
#print(len(aln2[0]))
#print(rand_P.prob(len(aln1[0]),len(aln2[0]),100,100))
#print(heuristic_score(aln1,aln2,ref_aln))
MGR = merger(ref_aln_file)
score,cons1,cons2 = MGR.logodd_merge(aln1,taxon_name1,aln2,taxon_name2,rand_P)

if outfile:
	fout = open(outfile,'w')
else:
	fout = stdout

for i in range(len(aln1)):
	fout.write(">"+taxon_name1[i]+"\n")
	fout.write(gap_propagate(cons1,aln1[i])+"\n")
for i in range(len(aln2)):
	fout.write(">"+taxon_name2[i]+"\n")
	fout.write(gap_propagate(cons2,aln2[i])+"\n")

if outfile:
	fout.close()
