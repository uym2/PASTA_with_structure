#! /usr/bin/env python

import sys
#import random
#import dendropy
import subprocess


list1 = ['ACGT','TAGC','GGTC','ACCG','TGGA','ATTG','AGTC','CTGT']
list2 = ['AAAA','CCCC','GGGG','TTTT','TGCA','ACGT','ACTG','CAGT']

#The Fastsp function
def Fastsp(list1,list2):
	name_list= []
	n = len(list1)
	m = len(list2)
	if n != m:
		print "Input lists have different sizes"
	f = open('ref.aln', 'w')
	for i in range(n):
		name_list.append('>'+'L'+ str(i)+'0273514')
#	print name_list
	for i in range(n):
		f.write(name_list[i])
		f.write('\n')
		f.write(list1[i])
		f.write('\n')
	f.close()
	f = open('est.aln', 'w')
	for i in range(n):
                f.write(name_list[i])
                f.write('\n')
                f.write(list2[i])
                f.write('\n')
	f.close()
	s = subprocess.check_output(["java","-jar","/home/uym2/Packages_N_Libraries/FastSP/FastSP.jar","-r","ref.aln","-e","est.aln"]).split()
#	print(s)
	return float(s[1]), float(s[3])

SP, Modeler = Fastsp(list1,list2)
print SP
print Modeler
