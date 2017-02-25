class rand_Aln:
	def __init__(self,n=0):
		self.n = n
		self.Q = {}

		self.Q[(0,0)] = 1
		for i in range(1,n+1):
			self.Q[(0,i)] = self.Q[(0,i-1)]*0.5

		for j in range(1,n+1):
			self.Q[(j,j)] = 0.25*self.Q[(j-1,j-1)] + self.Q[(j-1,j)]
			for i in range(j+1,n+1):
				self.Q[(j,i)] = 0.25*self.Q[(j-1,i-1)] + 0.5*self.Q[(j,i-1)] + 0.5*self.Q[(j-1,i)]
	
	def __q(self,m1=None,n1=None):
		m1 = self.n if m1 is None else m1
		n1 = self.n if n1 is None else n1
		
		if max(m1,n1) > self.n:
			print("Query number exceeds current size. Please call extend(max(" + str(m1) + "," + str(n1) + ")) then try again")
			return None
	
		return self.Q[(min(m1,n1),max(m1,n1))]	

	def extend(self,n1):
		if (n1 < self.n):
			print ("Already up-to-date with the specified size")
		else:
			for j in range(self.n+1,n1+1):
				self.Q[(j,j)] = 0.25*self.Q[(j-1,j-1)] + self.Q[(j-1,j)]
				for i in range(j+1,n+1):
					self.Q[(j,i)] = 0.25*self.Q[(j-1,i-1)] + 0.5*Q[(j,i-1)] + 0.5*Q[(j-1,i)]
			self.n = n1
			print ("Matrix Q extended successfully")

	def size(self):
		return self.n

	def prob(self,m1,n1,i,j):
	# probability of column i of the first seq of length m is matched with column j of the second seq of length n
		if (i > m1 or j > n1):
			print("Invalid input")
			return -1
		else:
			print(self.__q(i-1,j-1))
			print(self.__q(m1-i,n1-j))
			print(self.__q(m1,n1))
			return 0.25*self.__q(i-1,j-1)*self.__q(m1-i,n1-j)/float(self.__q(m1,n1))

	def del_rate(self,m1,n1,i):
	# deletion rate (prob of aligning ith column of seq 1 with a gap)
		p = 0
		for j in range(n1+1):
			p += 0.5*self.__q(i-1,j)*self.__q(m1-i,n1-j)/float(self.__q(m1,n1))
		return p	


	def ins_rate(self,m1,n1,j):
	# insertion rate (prob of aligning jth column of seq 2 with a gap)
		p = 0
		for i in range(m1+1):
			p += 0.5*self.__q(i,j-1)*self.__q(m1-i,n1-j)/float(self.__q(m1,n1))
		return p	
