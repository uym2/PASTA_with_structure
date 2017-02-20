from sys import argv
from sequence_lib import read_fasta, write_fasta


def ref_matching(aln1,aln2,ref_aln):
	# for now, assume everything is matched: taxon_name_all == taxon_name1 + taxon_name2 and the order is matched perfectly as well
	# also assume that everything in each alignment is perfectly aligned (all seqs have the same length)
	m = len(ref_aln[0])
	n = len(ref_aln)
	matching = [[-1 for x in range(m)] for y in range(n)]


	# match aln1 to ref_aln
	for j in range(len(aln1)):
		i1 = 0
		for i in range(m):
			if ref_aln[j][i] == '-':
				continue
			while aln1[j][i1] == '-':
				i1 += 1
			if aln1[j][i1] == ref_aln[j][i]:
				matching[j][i] = i1
				i1 += 1
				#print(matching)
			else:
				print("reference alignment and alignment 1 are not matched at sequence #" + str(j))



	# match aln2 to ref_aln

	for j in range(len(aln2)):
		i2 = 0
		for i in range(m):
			if ref_aln[j+len(aln1)][i] == '-':
				continue
			while aln2[j][i2] == '-':
				i2 += 1
			if aln2[j][i2] == ref_aln[j+len(aln1)][i]:
				matching[j+len(aln1)][i] = i2
				i2 += 1
			else:
				print("reference alignment and alignment 2 are not matched at sequence #" + str(j))
	return matching

def scoring(aln1,aln2,ref_aln):
	matching = ref_matching(aln1,aln2,ref_aln)
	m = len(matching[0])
	n = len(matching)
	scoring = {}
	
	for i in range(m):
		L1 = []
		d1 = {}
		L2 = []
		d2 = {}
		for j in range(n):
		# below L and d are used as "references" (just as pointers in C++): they are not deep copy, but a shallow copy of L1/L2 and d1/d2
			if j < len(aln1):
				L = L1
				d = d1
			else:
				L = L2
				d = d2
			if j == len(aln1) and L1 == []:
				L1 = L1 + [-1]		
				d1[-1] = len(aln1)
			if matching[j][i] >= 0:
				all_gap = False
				if matching[j][i] in d:
					d[matching[j][i]] += 1
				else:
					L += [matching[j][i]]
					d[matching[j][i]] = 1
		if L2 == []:
			L2 = L2 + [-1]
			d2[-1] = len(aln2)
		for x in L1:
			for y in L2:
				scoring[(x,y)] = d1[x]*d2[y]
	return scoring
	

aln1_file = argv[1]
aln2_file = argv[2]
ref_aln_file = argv[3]

taxon_name1, aln1 = read_fasta(aln1_file)
taxon_name2, aln2 = read_fasta(aln2_file)
taxon_name_all, ref_aln = read_fasta(ref_aln_file)

print(scoring(aln1,aln2,ref_aln))
