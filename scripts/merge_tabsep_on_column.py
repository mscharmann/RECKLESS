##
## merges two tab-separated files based on the strings in a column, fills NA if line not exists!

import sys

infile1 = sys.argv[1]
column_1 = int( sys.argv[2] )
infile2 = sys.argv[3]
column_2 = int( sys.argv[4] )



with open(infile1, "r") as INF:
	indict1 = { line.strip("\n").split("\t")[column_1] : line.strip("\n") for line in INF }

with open(infile2, "r") as INF:
	indict2 = { line.strip("\n").split("\t")[column_2] : line.strip("\n") for line in INF }

all_keys = set.union( set( indict1.keys() ), set( indict2.keys() ) )

outlines = []
for k in all_keys:	
	outl = ""	
	try:
		outl += indict1[k]
	except KeyError:
		outl += "\t".join( ["NA"]*len( next(iter(indict1.values())).split( "\t"   ) ) )
	outl += "\t"	
	try:
		outl += indict2[k]
	except KeyError:
		outl += "\t".join( ["NA"]*len( next(iter(indict2.values())).split( "\t"   ) ) )
	outlines.append( outl )
	

sys.stdout.write("\n".join(outlines)+"\n")

