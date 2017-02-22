from rand_aln_prob import rand_Aln
from sys import argv

m = int(argv[1])
n = int(argv[2])
i = int(argv[3])
j = int(argv[4])


R = rand_Aln(max(m,n))

print(R.prob(m,n,i,j))
print(R.in_rate(m,n,i))
print(R.del_rate(m,n,j))
