# RECKLESS - gRoup spEcifiC Kmers to pLant dEnovo aSSembly pipeline
RECKLESS - gRoup spEcifiC Kmers to pLant dEnovo aSSembly pipeline



inputs: short-reads (PE) of two groups

databases required: NCBI nt, lists of taxids under Viridiplantae subtrees etc

steps:
- get group-specific kmers | parameters for KMC
- retrieve reads containing these kmers: from each read-pair unit separately, and only from 
- assemble with Megahit
- filter for contigs containing original kmers
- BLASTN against nt: do not limit number of target seqs
- classify taxonomy
- filter for plant-contigs
- blastx against TAIR proteins, rename contigs if hit
- map reads to contigs
- filter contigs: drop those without mapped reads
- output: fasta of the final contigs, BAM files for reads mapped to those contigs
		
tools required:
- BLAST
- snakemake-minimal
- samtools
- bwa
- megahit
- seqtk
- kmc
- pandas

## simulate test data

A test dataset is composed as follows:

- 100kb from the first chromosome of Arabidopsis thaliana, Genbank CP002684.1
- 100kb from Pseudomonas syringiae, Genbank CP068034.2
- 100kb from the linear chromosome of Agrobacterium tumefaciens, CP033032.1

group-divergence exists both in the degree and species of contaminants AND in the plant data.


import random

with open("Ath_100kb.fasta", "r") as I:
	I.readline()
	inseq = I.readline().strip()
	

div_region = list(inseq[1000:11000])
snpsites = random.sample(range(0, 10000), 500)
print (len(snpsites))
for s in snpsites:
	s = int(s)
	isnuc = div_region[s]
	newnuc = random.choice([ x for x in ["A","C","T","G"] if not x == isnuc])
	div_region[s] = newnuc 
	
outseq = inseq[:1000] + "".join(div_region) + inseq[11000:]

with open("Ath_100kb.2nd_allele.fasta", "w") as OUT:
	OUT.write(">second_allele"+"\n")
	OUT.write(outseq+"\n")
	
### now simulate reads and compose the groups.
We wanted 100x coverage for the plant, i.e. c. 33000 PE150 reads, BUT there are contaminants which make up c. 10% of the reads.

wgsim -N 15000 -1 150 -2 150 Ath_100kb.fasta tmp1.1.fastq tmp1.2.fastq
wgsim -N 15000 -1 150 -2 150 Ath_100kb.2nd_allele.fasta tmp2.1.fastq tmp2.2.fastq
wgsim -N 3000 -1 150 -2 150 Pseu_100kb.fasta tmp3.1.fastq tmp3.2.fastq

cat tmp1.1.fastq tmp2.1.fastq tmp3.1.fastq | sed 's/\// /g' > sample_1.1.fastq
cat tmp1.2.fastq tmp2.2.fastq tmp3.2.fastq | sed 's/\// /g' > sample_1.2.fastq
rm tmp*

wgsim -N 30000 -1 150 -2 150 Ath_100kb.fasta tmp1.1.fastq tmp1.2.fastq
wgsim -N 3000 -1 150 -2 150 Agrob_100kb.fasta tmp3.1.fastq tmp3.2.fastq

cat tmp1.1.fastq tmp3.1.fastq | sed 's/\// /g' > sample_2.1.fastq
cat tmp1.2.fastq tmp3.2.fastq | sed 's/\// /g' > sample_2.2.fastq
rm tmp*

gzip *.fastq

#########

MUST use absolute paths in table data/pools_units_readfiles.txt
	
