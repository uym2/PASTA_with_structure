from sequence_lib import read_fasta
from rand_aln_prob import rand_Aln
from math import log

class merger:
	def __init__(self,ref_aln_file):
		taxa_names, ref_aln = read_fasta(ref_aln_file)
		self.ref_aln = ref_aln	
		self.tax2seqidx = {}

		for i in range(len(taxa_names)):
			self.tax2seqidx[taxa_names[i]] = i

	def show_taxa(self):
		return self.tax2seqidx.keys()

	def seqidx(self,tax):
		return self.tax2seqidx[tax]

	def ref_matching(self,aln1,taxa1,aln2,taxa2):
		m = len(self.ref_aln[0])
		n = len(aln1) + len(aln2)
		matching = [[-1 for x in range(m)] for y in range(n)]

		#taxa1, aln1 = read_fasta(aln1_file)
		#taxa2, aln2 = read_fasta(aln2_file)

		# match aln1 to ref_aln
		for j1 in range(len(aln1)):
			tax = taxa1[j1]
			j = self.seqidx(tax)
			i1 = 0
			gap1 = 0
			for i in range(m):
				if self.ref_aln[j][i] == '-':
					gap1 += 1
					continue
				while aln1[j1][i1] == '-':
					i1 += 1
				if aln1[j1][i1] == self.ref_aln[j][i]:
					matching[j1][i] = i1
					i1 += 1
					#print(matching)
				else:
					print("reference alignment and alignment 1 are not matched at taxon " + tax)
		gap_rate1 = float(gap1)/len(self.ref_aln[0])/len(aln1)
		# match aln2 to ref_aln
		for j2 in range(len(aln2)):
			tax = taxa2[j2]
			j = self.seqidx(tax)
			i2 = 0
			gap2 = 0
			for i in range(m):
				if self.ref_aln[j][i] == '-':
					gap2 += 1
					continue
				while aln2[j2][i2] == '-':
					i2 += 1
				if aln2[j2][i2] == self.ref_aln[j][i]:
					matching[j2+len(aln1)][i] = i2
					i2 += 1
				else:
					print("reference alignment and alignment 2 are not matched at taxon " + tax)
		gap_rate2 = float(gap2)/len(self.ref_aln[0])/len(aln2)
		return matching,gap_rate1,gap_rate2

	def heuristic_score(self,aln1,taxa1,aln2,taxa2,norm=False):
		matching,gap_rate1,gap_rate2 = self.ref_matching(aln1,taxa1,aln2,taxa2)
		m = len(matching[0])
		n = len(matching)
		pairing = {}
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
					if matching[j][i] in d:
						d[matching[j][i]] += 1
					else:
						L += [matching[j][i]]
						d[matching[j][i]] = 1
			for x in L1:
				for y in L2:
					if (x,y) not in pairing:
						pairing[(x,y)] = 0
					pairing[(x,y)] += d1[x]*d2[y]
		if norm:
			for i in range(len(aln1[0])):
				for j in range(len(aln2[0])):
					if (i,j) in pairing:
							
						count = 0
						for k in range(len(aln2[0])):
							count = count + pairing[(i,k)] if (i,k) in pairing else count
						for k in range(len(aln1[0])):
							count = count + pairing[(k,j)] if (k,j) in pairing else count
						count -= pairing[(i,j)]
						scoring[(i,j)] = float(pairing[(i,j)])/count
		else:					
			scoring = pairing
		
		return scoring,gap_rate1,gap_rate2
	
	def logodd_score(self,aln1,taxa1,aln2,taxa2,rand_P):
		scoring,gap_rate1,gap_rate2 = self.heuristic_score(aln1,taxa1,aln2,taxa2,norm=True)
		#print(gap_rate1)
		#print(gap_rate2)
		for key in scoring:
			#print(rand_P.prob(len(aln1[0]),len(aln2[0]),key[0]+1,key[1]+1))
			scoring[key] = log(scoring[key]/rand_P.prob(len(aln1[0]),len(aln2[0]),key[0]+1,key[1]+1))
		#print(scoring)
		del_score = log(gap_rate2/rand_P.del_rate(len(aln1[0]),len(aln2[0]),1))
		ins_score = log(gap_rate1/rand_P.ins_rate(len(aln1[0]),len(aln2[0]),1))
		#print(del_score)
		#print(ins_score)
		return scoring,del_score,ins_score
	
	
	def merge(self,n,m,scoring,default=0,ins_score=0,del_score=0):
		aln_score = [[0 for i in range(m+1)] for j in range(n+1)]
		backtrack = [['-' for i in range(m+1)] for j in range(n+1)]
		for i in range(1,m+1):
			backtrack[0][i] = 'L'
			aln_score[0][i] = aln_score[0][i-1] + ins_score
		for j in range(1,n+1):
			backtrack[j][0] = 'U'
			aln_score[j][0] = aln_score[j-1][0] + del_score

		for j in range(1,n+1):
			for i in range(1,m+1):
				ms  = aln_score[j-1][i-1] + scoring[(j-1,i-1)] if (j-1,i-1) in scoring else aln_score[j-1][i-1] + default
				g1 = aln_score[j][i-1] + ins_score
				g2 = aln_score[j-1][i] + del_score

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

	def heuristic_merge(self,aln1,taxa1,aln2,taxa2):
		scoring = self.heuristic_score(aln1,taxa1,aln2,taxa2)
		m = len(aln2[0])
		n = len(aln1[0])
		return self.merge(n,m,scoring)

	def logodd_merge(self,aln1,taxa1,aln2,taxa2,rand_P):
		scoring,del_score,ins_score = self.logodd_score(aln1,taxa1,aln2,taxa2,rand_P)
		print(del_score)
		print(ins_score)
		m = len(aln2[0])
		n = len(aln1[0])
		#return self.merge(n,m,scoring,default=0,ins_score=0,del_score=0)
		return self.merge(n,m,scoring)
