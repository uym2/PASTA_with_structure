from sys import argv,stdout
from sequence_lib import read_fasta, write_fasta
from math import log

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

def expected_gap(l1,l2,n,m,ref_aln):
#l1 : len(aln1) ; l2 : len(aln2)
	eg1 = 0
	eg2 = 0
	for i in range(l1):
		for x in ref_aln[i]:
			if x == '-':
				eg1 += 1
	
	for i in range(l1,l1+l2):
		for x in ref_aln[i]:
			if x == '-':
				eg2 += 1

	return eg1/n, eg2/m

def estimate_score(aln1,aln2,ref_aln):
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
			#if j == len(aln1) and L1 == []:
				#L1 = L1 + [-1]		
				#d1[-1] = len(aln1)
			if matching[j][i] >= 0:
				all_gap = False
				if matching[j][i] in d:
					d[matching[j][i]] += 1
				else:
					L += [matching[j][i]]
					d[matching[j][i]] = 1
		#if L2 == []:
		#	L2 = L2 + [-1]
		#	d2[-1] = len(aln2)
		for x in L1:
			for y in L2:
				if (x,y) not in scoring:
					scoring[(x,y)] = 0
				scoring[(x,y)] += d1[x]*d2[y]
	return scoring

def merge(aln1,aln2,ref_aln):
	scoring = estimate_score(aln1,aln2,ref_aln)
	m = len(aln2[0])
	n = len(aln1[0])

	l1 = len(aln1)
	l2 = len(aln2)

	eg1,eg2 = expected_gap(l1,l2,n,m,ref_aln)
	gap_penalty = (l1-eg1)*(l2-eg2)/10
	#print(gap_penalty)

	aln_score = [[0 for i in range(m+1)] for j in range(n+1)]
	backtrack = [['-' for i in range(m+1)] for j in range(n+1)]
	for i in range(1,m+1):
		backtrack[0][i] = 'L'
		#aln_score[0][i] = aln_score[0][i-1] + scoring[(-1,i-1)] if (-1,i-1) in scoring else aln_score[0][i-1]
	for j in range(1,n+1):
		backtrack[j][0] = 'U'
		#aln_score[j][0] = aln_score[j-1][0] + scoring[(j-1,-1)] if (j-1,-1) in scoring else aln_score[j-1][0]

	for j in range(1,n+1):
		for i in range(1,m+1):
			ms  = aln_score[j-1][i-1] + scoring[(j-1,i-1)] if (j-1,i-1) in scoring else aln_score[j-1][i-1]
			g1 = aln_score[j][i-1] #- gap_penalty
			g2 = aln_score[j-1][i] #- gap_penalty

			if ms >= g1 and ms >= g2:
				aln_score[j][i] = ms
				backtrack[j][i] = 'D'
			elif g1 >= ms and g1 >= g2:
				aln_score[j][i] = g1
				backtrack[j][i] = 'L'
			else:
				aln_score[j][i] = g2
				backtrack[j][i] = 'U'

	i = m
	j = n
	M1 = ""
	M2 = ""
	while (i > 0 or j > 0):
		#print(aln_score[j][i])
		if backtrack[j][i] == 'D':
			#print('D')
			#M1 = str(j-1) + M1
			#M2 = str(i-1) + M2
			M1 = "." + M1
			M2 = "." + M2
			i -= 1
			j -= 1	
		elif backtrack[j][i] == 'L':
			#print('L')
			#M2 = str(i-1) + M2
			M2 = "." + M2
			M1 = "-" + M1
			i -= 1
		else:
			#print('U')
			M2 = "-" + M2
			#M1 = str(j-1) + M1
			M1 = "." + M1
			j -= 1
	return aln_score[n][m], M1, M2		

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

aln1_file = argv[1]
aln2_file = argv[2]
ref_aln_file = argv[3]
if len(argv) > 4:
	outfile = argv[4]
else:	
	outfile = None

taxon_name1, aln1 = read_fasta(aln1_file)
taxon_name2, aln2 = read_fasta(aln2_file)
taxon_name_all, ref_aln = read_fasta(ref_aln_file)

for i in range(len(aln1)):
	aln1[i] = aln1[i].upper()
for i in range(len(aln2)):
	aln2[i] = aln2[i].upper()
for i in range(len(ref_aln)):
	ref_aln[i] = ref_aln[i].upper()


#print(estimate_score(aln1,aln2,ref_aln))
score,cons1,cons2 = merge(aln1,aln2,ref_aln)

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

