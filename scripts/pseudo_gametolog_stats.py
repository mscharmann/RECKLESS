## pseudo_gametolog_stats.py

"""
INPUT: VCF called by bcftools. pipe it.
"""


import sys

sys.stdout.write( "\t".join(["chrom","XY_like","ZW_like"]) + "\n" )

#print(popdict_vcf_idx)

chrom_current = "dummy"
XY_like = 0
ZW_like = 0
for line in sys.stdin:
	if line.startswith("#"):
		if line.startswith("#CHROM"):
			header_line = line.lstrip("#").strip("\n").split("\t")
			samples = sorted(header_line[9:])
			continue
		continue
	if len(line) < 2: # empty lines or so
		continue
	fields = line.strip("\n").split("\t")
	chrom = fields[0]
	if chrom != chrom_current:
		# triggers report and reset of counters
		sys.stdout.write( "\t".join([chrom_current,str(XY_like),str(ZW_like)]) + "\n" )
		chrom_current = chrom
		XY_like = 0
		ZW_like = 0
		
	pos = fields[1]
	refallele = fields[3]
	altallele = fields[4]
	
	# p1 are the 'males', p2 are the 'females'
	
	gts_p1 = set( fields[9].split(":")[0].replace("/","") )		
	gts_p2 = set( fields[10].split(":")[0].replace("/","") )
	
	if not gts_p1 == set("."):
		if not gts_p2 == set("."):
			if len(gts_p1) > 1:
				if len(gts_p2) == 1:
					allelic_depths_p1 = [int(x) for x in fields[9].split(":")[3].split(",")]
					if 0.4 <= allelic_depths_p1[0]/float(sum(allelic_depths_p1)) <= 0.6:		
						# print("XY-like", allelic_depths_p1[0]/float(sum(allelic_depths_p1)), chrom, fields[9],fields[10])
						XY_like += 1
			elif len(gts_p2) > 1:
				if len(gts_p1) == 1:
					allelic_depths_p2 = [int(x) for x in fields[10].split(":")[3].split(",")]
					if 0.4 <= allelic_depths_p2[0]/float(sum(allelic_depths_p2)) <= 0.6:		
						# print("ZW-like", allelic_depths_p2[0]/float(sum(allelic_depths_p2)), chrom, fields[9],fields[10])
						ZW_like += 1



# final chrom report
sys.stdout.write( "\t".join([chrom_current,str(XY_like),str(ZW_like)]) + "\n" )







