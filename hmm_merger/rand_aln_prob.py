class rand_Aln:
	def __init__(self):
		self.n = 0
		self.P = {}
		self.S = {}

	def __init__(self,n):
		self.n = n
		self.setup(n)

	def setup(self,n):
		self.n = n
		self.P = {}
		self.S = {}

		for i in range(n+1):
			self.P[(0,i)] = 1
			self.S[(0,i)] = 1

		for j in range(1,n+1):
			self.P[(j,j)] = 2*self.S[(j-1,j-1)] + self.P[(j-1,j)]
			self.S[(j,j)] = self.S[(j-1,j)] + self.P[(j,j)]
			for i in range(j+1,n+1):
				self.P[(j,i)] = 2*self.S[(j-1,i-1)] + self.P[(j,i-1)]
				self.S[(j,i)] = self.S[(j-1,i)] + self.P[(j,i)]
	
	def count(self,m1=None,n1=None):
		m1 = self.n if m1 is None else m1
		n1 = self.n if n1 is None else n1
		
		if max(m1,n1) > self.n:
			print("Query number exceeds current size. Please call update(max(" + str(m1) + "," + str(n1) + ")) then try again")
			return None
	
		return self.P[(min(m1,n1),max(m1,n1))]	

	def update(self,n1):
		if (n1 < self.n):
			print ("Already up-to-date with the specified size")
		else:
			for j in range(self.n+1,n1+1):
				self.P[(j,j)] = 2*self.S[(j-1,j-1)] + self.P[(j-1,j)]
				self.S[(j,j)] = self.S[(j-1,j)] + self.P[(j,j)]
				for i in range(j+1,n+1):
					self.P[(j,i)] = 2*self.S[(j-1,i-1)] + self.P[(j,i-1)]
					self.S[(j,i)] = self.S[(j-1,i)] + self.P[(j,i)]
			self.n = n1
			print ("Updated successfully")

	def size(self):
		return self.n


	def prob(self,m1,n1,i,j):
	# probability of column i of the first seq of length m is matched with column j of the second seq of length n
		if (i > m1 or j > n1):
			print("Invalid input")
			return -1
		else:
			return self.count(i-1,j-1)*self.count(m1-i,n1-j)/float(self.count(m1,n1))

	def in_rate(self,m1,n1,i):
	# insertion rate (prob of aligning ith column of seq 1 with a gap)
		p = 0
		for j in range(n1+1):
			p += self.count(i-1,j)*self.count(m1-i,n1-j)/float(self.count(m1,n1))
		return p	


	def del_rate(self,m1,n1,j):
	# deletion rate (prob of aligning jth column of seq 2 with a gap)
		p = 0
		for i in range(m1+1):
			p += self.count(i,j-1)*self.count(m1-i,n1-j)/float(self.count(m1,n1))
		return p	
