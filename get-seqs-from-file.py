#!/usr/bin/python

from Bio import SeqIO
import sys

listfile1 = open(sys.argv[1],'r') 

inputfile = open(sys.argv[2],'r') 

g = open(sys.argv[3],'w') 

infmt = sys.argv[4]
outfmt=sys.argv[5]

x = SeqIO.to_dict(SeqIO.parse(inputfile,infmt)) 
n=0
c=0
for i in listfile1:
	if i[0]<>"#":
		
		iname = i.split(" ")[0].rstrip("\n")
		try:
			SeqIO.write(x[iname],g,outfmt)
			c=c+1
		except:
			n=n+1

print "got %d sequences, could not find %d " %(c,n)








	

		






