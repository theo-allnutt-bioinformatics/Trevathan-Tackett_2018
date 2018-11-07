#!/usr/bin/python
import sys


f = open(sys.argv[1],'r') #otu table in txt
ref = open(sys.argv[2],'r') #gg taxonomy file
g = open(sys.argv[3],'w') #output
h=open("chl.hits",'w')#chlorplast otus
data={}
terms=sys.argv[4].split(",")
c=0
for x in ref:

	data[x.split("\t")[0]]=x.split("\t")[1].rstrip("\n")
	

for i in f:
	
	if i[0]<>"#":
	
		k = i.split("\t")
		
		gg=k[0]
		
		try:
			tax=data[gg]
		except:
			tax="unknown"
			
		g.write(i.rstrip("\n")+"\t"+tax+"\n")
		
		for t in terms:
			
			if t.lower() in tax.lower():
				c=c+1
				print gg,tax
				h.write(gg+"\t"+tax+"\n")
				break
				
	else:
		g.write(i.rstrip("\n")+"\ttaxonomy\n")
	
print "green hits=",c
