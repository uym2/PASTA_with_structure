from sequence_lib import read_fasta
from rand_aln_prob import rand_Aln
from math import log,sqrt
from FastSP import FastSP

class merger:
	def __init__(self,ref_aln_file):
		taxa_names, ref_aln = read_fasta(ref_aln_file)
		self.ref_aln = ref_aln	
		self.tax2seqidx = {}
		self.wp = 1 # P(true|ref_aln): probability that a homology is true given it is found in ref_aln (modeler)
		self.wn = 0 # P(true|~ref_aln): probability that a homology not found in ref_aln is true
			    # = P(true|~ref_aln) = P(~ref_aln|true)*P(true)/P(~ref_aln) = SPFN*P(true)/P(~ref_aln)

		for i in range(len(taxa_names)):
			self.tax2seqidx[taxa_names[i]] = i

	def show_taxa(self):
		return self.tax2seqidx.keys()

	def seqidx(self,tax):
		return self.tax2seqidx[tax]

	def __eval(self,aln,taxa):
		M = []
		count = 0
		subtr = 0

		for i,tax in enumerate(taxa):
			M.append(self.ref_aln[self.tax2seqidx[tax]])
			l = len([x for x in aln[i] if x != '-'])
			count += l
			subtr += l*(l-1)/2

		hshared,hA,hM = FastSP(aln,M)
		h_all=count*(count-1)/2 - subtr
		wp = 1.0*hshared/hM
		wn = (1.0-hshared/hA)*hA/(h_all-hM)
		
		return wp,wn

	def evalRef(self,aln1,taxa1,aln2=None,taxa2=None):
		wp1,wn1 = self.__eval(aln1,taxa1)

		if aln2 is None or taxa2 is None:
			self.wp = wp1
			self.wn = wn1
		else:
			wp2,wn2 = self.__eval(aln2,taxa2)
			self.wp = sqrt(wp1*wp2)
			self.wn = sqrt(wn1*wn2)
		

	def ref_matching(self,aln1,taxa1,aln2,taxa2):
		m = len(self.ref_aln[0])
		n = len(aln1) + len(aln2)
		matching = [[-1 for x in range(m)] for y in range(n)]
		#match1 = [[-1 for x in range(len(aln1[0])+1)] for y in range(len(aln1))] # the last column stores the length of each sequence
		#match2 = [[-1 for x in range(len(aln2[0])+1)] for y in range(len(aln2))] # the last column stores the length of each sequence

		# match aln1 to ref_aln
		for j1 in range(len(aln1)):
			tax = taxa1[j1]
			j = self.seqidx(tax)
			i1 = 0
			k = 0
			#gap1 = 0
			for i in range(m):
				if self.ref_aln[j][i] == '-':
					#gap1 += 1
					continue
				while aln1[j1][i1] == '-':
					i1 += 1
				#match1[j1][i1] = k
				k += 1
				if aln1[j1][i1] == self.ref_aln[j][i]:
					matching[j1][i] = i1
					i1 += 1
				else:
					print("reference alignment and alignment 1 are not matched at taxon " + tax)
			#match1[j1][len(aln1[0])] = k
		# match aln2 to ref_aln
		for j2 in range(len(aln2)):
			tax = taxa2[j2]
			j = self.seqidx(tax)
			i2 = 0
			k = 0
			gap2 = 0
			for i in range(m):
				if self.ref_aln[j][i] == '-':
					#gap2 += 1
					continue
				while aln2[j2][i2] == '-':
					i2 += 1
				#match2[j2][i2] = k
				k += 1
				if aln2[j2][i2] == self.ref_aln[j][i]:
					matching[j2+len(aln1)][i] = i2
					i2 += 1
				else:
					print("reference alignment and alignment 2 are not matched at taxon " + tax)
			#match2[j2][len(aln2[0])] = k
		return matching #,match1,match2

	def residue_count(self,aln):
		R = []
		for j in range(len(aln[0])):
			R += [0]
			for i in range(len(aln)):
				if aln[i][j] != '-':
					R[j] += 1
		return R

	def heuristic_score(self,aln1,taxa1,aln2,taxa2,w=1):
		R1 = self.residue_count(aln1)
		R2 = self.residue_count(aln2)		
		matching = self.ref_matching(aln1,taxa1,aln2,taxa2)
		m = len(matching[0])
		n = len(matching)
		score_tab = {}
		
		for i in range(len(aln1[0])):
			score_tab[(i,-1)] = 0

		for j in range(len(aln2[0])):
			score_tab[(-1,j)] = 0
		
		for i in range(len(aln1[0])):
			for j in range(len(aln2[0])):
				score_tab[(i,j)] = (w-1)*R1[i]*R2[j]
			
		for i in range(m):
			L1 = []
			d1 = {}
			L2 = []
			d2 = {}
			for j in range(n):
			# below L and d are used as "references" (just as pointers in C++): they are not deep copy, but a shallow copy of L1|L2 and d1|d2
				if j < len(aln1):
					L = L1
					d = d1
				else:
					L = L2
					d = d2
				if matching[j][i] >= 0:
					if matching[j][i] in d:
						d[matching[j][i]] += 1
					else:
						L += [matching[j][i]]
						d[matching[j][i]] = 1
			for x in L1:
				for y in L2:
					score_tab[(x,y)] += (2*w-1)*d1[x]*d2[y]


		return score_tab,R1,R2 #,m1,m2

	def TP_score(self,aln1,taxa1,aln2,taxa2):
		return self.heuristic_score(aln1,taxa1,aln2,taxa2,w=1)

	def llh_score(self,aln1,taxa1,aln2,taxa2):
		TP_score,R1,R2 = self.TP_score(aln1,taxa1,aln2,taxa2)
		score_tab = {}

		self.evalRef(aln1,taxa1,aln2,taxa2)
		print(self.wp,self.wn)

		for j in range(len(aln2[0])):
			score_tab[(-1,j)] = 0
			for i in range(len(aln1[0])):
				H_ij = R1[i]*R2[j]
				p = 1.0-(self.wp*TP_score[(i,j)]+self.wn*(H_ij-TP_score[(i,j)]))/H_ij
				score_tab[(-1,j)] += log(p)

		for i in range(len(aln1[0])):
			score_tab[(i,-1)] = 0
			for j in range(len(aln2[0])):
				H_ij = R1[i]*R2[j]
				p = 1.0-(self.wp*TP_score[(i,j)]+self.wn*(H_ij-TP_score[(i,j)]))/H_ij
				#print(p)
				score_tab[(i,-1)] += log(p)

		for (i,j) in TP_score:
			H_ij = R1[i]*R2[j]
			score_tab[(i,j)] = log((self.wp*TP_score[(i,j)]+self.wn*(H_ij-TP_score[(i,j)]))/H_ij)

		return score_tab
	
	def logodd_score(self,aln1,taxa1,aln2,taxa2,rand_P):
		#score_tab,match1,match2= self.heuristic_score(aln1,taxa1,aln2,taxa2)
		score_tab,R1,R2 = self.heuristic_score(aln1,taxa1,aln2,taxa2)
		#print(gap_rate1)
		#print(gap_rate2)
		for key in score_tab:
			'''
			print(key)
			p = 0
			for s1 in match1:
				for s2 in match2:
					if s1[key[0]] >= 0 and s2[key[1]] >= 0 :
						p += rand_P.prob(s1[-1],s2[-1],s1[key[0]]+1,s2[key[1]]+1)
						#print(p)
			#score_tab[key] = log(score_tab[key]/rand_P.prob(len(aln1[0]),len(aln2[0]),key[0]+1,key[1]+1))
			p = p/len(aln1)/len(aln2)
			#print(p)
			score_tab[key] = log(score_tab[key]/p)
			#print(score_tab[key])
		#print(score_tab)
			'''
			if score_tab[key]:
				score_tab[key] = log(score_tab[key]/rand_P.prob(len(aln1[0]),len(aln2[0]),key[0]+1,key[1]+1))
		#del_score = log(gap_rate2/rand_P.del_rate(len(aln1[0]),len(aln2[0]),1))
		#ins_score = log(gap_rate1/rand_P.ins_rate(len(aln1[0]),len(aln2[0]),1))
		return score_tab #,del_score,ins_score
	
	
	def merge(self,aln1,aln2,score_tab,default=0):
		n = len(aln1[0])
		m = len(aln2[0])

		aln_score = [[0 for i in range(m+1)] for j in range(n+1)]
		backtrack = [['-' for i in range(m+1)] for j in range(n+1)]
		for i in range(1,m+1):
			backtrack[0][i] = 'L'
			aln_score[0][i] = aln_score[0][i-1] + score_tab[(-1,i-1)]#+ ins_score
		for j in range(1,n+1):
			backtrack[j][0] = 'U'
			aln_score[j][0] = aln_score[j-1][0] + score_tab[(j-1,-1)]#+ del_score

		for j in range(1,n+1):
			for i in range(1,m+1):
				ms = aln_score[j-1][i-1] + score_tab[(j-1,i-1)] 
				g1 = aln_score[j][i-1] + score_tab[(-1,i-1)] #+ ins_score
				g2 = aln_score[j-1][i] + score_tab[(j-1,-1)] #+ del_score

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
			if backtrack[j][i] == 'D':
				M1 = "." + M1
				M2 = "." + M2
				i -= 1
				j -= 1	
			elif backtrack[j][i] == 'L':
				M2 = "." + M2
				M1 = "-" + M1
				i -= 1
			else:
				M2 = "-" + M2
				M1 = "." + M1
				j -= 1
		return aln_score[n][m], M1, M2		

	def heuristic_merge(self,aln1,taxa1,aln2,taxa2):
		score_tab,R1,R2 = self.heuristic_score(aln1,taxa1,aln2,taxa2)
		return self.merge(aln1,aln2,score_tab)

	def ML_merge(self,aln1,taxa1,aln2,taxa2):
		score_tab = self.llh_score(aln1,taxa1,aln2,taxa2)
		return self.merge(aln1,aln2,score_tab)

	def logodd_merge(self,aln1,taxa1,aln2,taxa2,rand_P):
		score_tab= self.logodd_score(aln1,taxa1,aln2,taxa2,rand_P)
		return self.merge(aln1,aln2,score_tab)
