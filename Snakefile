configfile: "config.yaml"

samples_units_fqs_map = config["samples_units_fqs_map"]


import pandas as pd

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) ) 
#print(samples_units_fqs)
#print(SAMPLES)

groups = ["group1","group2"]


# ci_param_combs = [[20,2],[50,2],[100,2],[250,2],[20,10],[50,10],[100,10],[250,10],[50,20],[100,20],[250,20]]
ci_param_combs = [[50,2],[100,2],[250,2],[50,10],[100,10],[250,10],[100,20],[250,20]]



ruleorder: KMC_substract_and_report_group1 > KMC_substract_and_report_group2 > KMC_count

def get_sample_spRead_fqs(wildcards):
	"""Get all fq files with sp reads of given sample."""
	return expand(
		"mapped_reads_per_unit/{sample}-{unit}.sorted.bam",
		sample=wildcards.sample,
		unit=samples_units_fqs.loc[wildcards.sample].unit,
	)


def get_fastq_sample_unit(wildcards):
	"""Get fastq files of given sample-unit."""
	fastqs = samples_units_fqs.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]
	if fastqs.fq2.isnull().values.any():
		return [ fastqs.fq1.item() ]
	return [ fastqs.fq1.item(), fastqs.fq2.item() ]


def get_fastq_sample_ALL(wildcards):
	"""Get list of fastq files of given sample including ALL units."""
	fastqs_pd = samples_units_fqs.loc[(wildcards.sample), ["fq1", "fq2"]]
	fastqs = set( fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist() )
	fastqs_clean = list( {x for x in fastqs if pd.notna(x)} )
	return fastqs_clean


rule all:
	input:
		expand("results/asm.{group}.cir{param_comb[0]}.cie{param_comb[1]}.stats_per_contig.txt", group=groups, param_comb = ci_param_combs),
		"working_dir/all_credible_contigs.fa.taxonomy_report.txt",
		"results/report_plant_contig_stats.txt"


rule softlink_inputs_and_protocol_renaming:
	output:
		r11="working_dir/group1.1.fq.gz",
		r12="working_dir/group1.2.fq.gz",
		r21="working_dir/group2.1.fq.gz",
		r22="working_dir/group2.2.fq.gz",	
		report="results/renaming_protocol.txt"
	run:
		# softlink read files: stable, hardcoded names.
		fastqs_pd = samples_units_fqs.loc[SAMPLES[0], ["fq1", "fq2"]]
		fastqs_clean = fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist()
		shell("ln -s {fastqs_clean[0]} {output.r11}")
		shell("ln -s {fastqs_clean[1]} {output.r12}")
	
		fastqs_pd = samples_units_fqs.loc[SAMPLES[1], ["fq1", "fq2"]]
		fastqs_clean = fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist()
		shell("ln -s {fastqs_clean[0]} {output.r21}")
		shell("ln -s {fastqs_clean[1]} {output.r22}")
	
		with open("results/renaming_protocol.txt", "w") as O:
			O.write(SAMPLES[0] + " renamed as group1" + "\n")
			O.write(SAMPLES[1] + " renamed as group2"+ "\n")
			O.write("these names apply throughout the analyses"+ "\n")
		
		shell("echo group min_ci_to_retain min_ci_to_exclude number_of_specific_kmers > results/kmer_spec_report.txt")

	

rule KMC_count:
	input:
		expand("working_dir/{group}.1.fq.gz", group=groups)
	output:
		"results/kmer_{group}.kmc_suf"
	params:
		ks = config["kmer_size"]
	threads: 14
	shell:
		"""
		# count kmers that occur at least 2 times (-ci), and count them up to 1 million (-cs)
		echo "working_dir/{wildcards.group}.1.fq.gz" "working_dir/{wildcards.group}.2.fq.gz" | tr " " "\n" >input_file_names.{wildcards.group}
		mkdir -p ./kmc_tempdir.{wildcards.group}
		kmc -m5 -sm -k{params.ks} -fq -ci2 -cs1000000 -t{threads} @input_file_names.{wildcards.group} results/kmer_{wildcards.group} ./kmc_tempdir.{wildcards.group}
		rm input_file_names.{wildcards.group}
		rm -r ./kmc_tempdir.{wildcards.group}
		"""


rule KMC_substract_and_report_group1:
	input:
		"results/kmer_group1.kmc_suf",
		"results/kmer_group2.kmc_suf"
	output:
		"results/kmer_diff.group1_minus_group2.cir{p1}.cie{p2}.kmc_suf"
	threads: 14
	shell:
		"""
		kmc_tools -t{threads} simple results/kmer_group1 -ci{wildcards.p1} results/kmer_group2 -ci{wildcards.p2} kmers_subtract results/kmer_diff.group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2}
		echo group1 {wildcards.p1} {wildcards.p2} $( kmc_tools transform results/kmer_diff.group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2} dump /dev/stdout | wc -l ) >> "results/kmer_spec_report.txt"			
		"""

rule KMC_substract_and_report_group2:
	input:
		"results/kmer_group1.kmc_suf",
		"results/kmer_group2.kmc_suf"
	output:
		"results/kmer_diff.group2_minus_group1.cir{p1}.cie{p2}.kmc_suf"
	threads: 14
	shell:
		"""
		kmc_tools -t{threads} simple results/kmer_group2 -ci{wildcards.p1} results/kmer_group1 -ci{wildcards.p2} kmers_subtract results/kmer_diff.group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2}
		echo group2 {wildcards.p1} {wildcards.p2} $( kmc_tools transform results/kmer_diff.group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2} dump /dev/stdout | wc -l ) >> "results/kmer_spec_report.txt"			
		"""


rule get_reads_containing_specific_kmers_group1:
	input:
		"results/kmer_diff.group1_minus_group2.cir{p1}.cie{p2}.kmc_suf"
	output:
		r1=temp("working_dir/target_reads_group1_minus_group2.cir{p1}.cie{p2}.1.fq.gz"),
		r2=temp("working_dir/target_reads_group1_minus_group2.cir{p1}.cie{p2}.2.fq.gz")
	threads: 14
	shell:
		"""
		echo working_dir/group1.1.fq.gz working_dir/group1.2.fq.gz | tr " " "\n" >input_file_names.getreads_group1.cir{wildcards.p1}.cie{wildcards.p2}
		
		kmc_tools -t{threads} filter results/kmer_diff.group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2} -ci1 @input_file_names.getreads_group1.cir{wildcards.p1}.cie{wildcards.p2} -ci1 tmp.group1.cir{wildcards.p1}.cie{wildcards.p2}.fq
 		cat tmp.group1.cir{wildcards.p1}.cie{wildcards.p2}.fq | awk '/^@/ {{print $1}}' | sort | uniq | sed 's/@//g' | sed 's/\/.//g' > working_dir/target_reads_group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2}.txt
		rm tmp.group1.cir{wildcards.p1}.cie{wildcards.p2}.fq input_file_names.getreads_group1.cir{wildcards.p1}.cie{wildcards.p2}

		seqtk subseq working_dir/group1.1.fq.gz working_dir/target_reads_group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2}.txt | gzip -c > {output.r1}
		seqtk subseq working_dir/group1.2.fq.gz working_dir/target_reads_group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2}.txt | gzip -c > {output.r2}
		
		rm working_dir/target_reads_group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2}.txt
		
		"""
		
		
rule get_reads_containing_specific_kmers_group2:
	input:
		"results/kmer_diff.group2_minus_group1.cir{p1}.cie{p2}.kmc_suf"
	output:
		r1=temp("working_dir/target_reads_group2_minus_group1.cir{p1}.cie{p2}.1.fq.gz"),
		r2=temp("working_dir/target_reads_group2_minus_group1.cir{p1}.cie{p2}.2.fq.gz")
	threads: 14
	shell:
		"""
		echo working_dir/group2.1.fq.gz working_dir/group2.2.fq.gz | tr " " "\n" >input_file_names.getreads_group2.cir{wildcards.p1}.cie{wildcards.p2}
		
		kmc_tools -t{threads} filter results/kmer_diff.group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2} -ci1 @input_file_names.getreads_group2.cir{wildcards.p1}.cie{wildcards.p2} -ci1 tmp.group2.cir{wildcards.p1}.cie{wildcards.p2}.fq
 		cat tmp.group2.cir{wildcards.p1}.cie{wildcards.p2}.fq | awk '/^@/ {{print $1}}' | sort | uniq | sed 's/@//g' | sed 's/\/.//g' > working_dir/target_reads_group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2}.txt
		rm tmp.group2.cir{wildcards.p1}.cie{wildcards.p2}.fq input_file_names.getreads_group2.cir{wildcards.p1}.cie{wildcards.p2}
		
		seqtk subseq working_dir/group2.1.fq.gz working_dir/target_reads_group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2}.txt | gzip -c > {output.r1}
		seqtk subseq working_dir/group2.2.fq.gz working_dir/target_reads_group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2}.txt | gzip -c > {output.r2}
		
		rm working_dir/target_reads_group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2}.txt
		
		"""


rule asm_megahit_group1:
	input:
		r1="working_dir/target_reads_group1_minus_group2.cir{p1}.cie{p2}.1.fq.gz",
		r2="working_dir/target_reads_group1_minus_group2.cir{p1}.cie{p2}.2.fq.gz"		
	output:
		"working_dir/asm.group1.cir{p1}.cie{p2}.raw.fa"
	threads: 12
	shell:
		"""
		# -m : fraction of all memory on the machine to use
		megahit -m 0.2 -t {threads} -1 {input.r1} -2 {input.r2} --no-mercy -o working_dir/megahit1_cir{wildcards.p1}.cie{wildcards.p2}
		mv working_dir/megahit1_cir{wildcards.p1}.cie{wildcards.p2}/final.contigs.fa {output}
		rm -r working_dir/megahit1_cir{wildcards.p1}.cie{wildcards.p2}		
		"""
	

rule asm_megahit_group2:
	input:
		r1="working_dir/target_reads_group2_minus_group1.cir{p1}.cie{p2}.1.fq.gz",
		r2="working_dir/target_reads_group2_minus_group1.cir{p1}.cie{p2}.2.fq.gz"		
	output:
		"working_dir/asm.group2.cir{p1}.cie{p2}.raw.fa"
	threads: 12
	shell:
		"""
		# -m : fraction of all memory on the machine to use
		megahit -m 0.2 -t {threads} -1 {input.r1} -2 {input.r2} --no-mercy -o working_dir/megahit2_cir{wildcards.p1}.cie{wildcards.p2}
		mv working_dir/megahit2_cir{wildcards.p1}.cie{wildcards.p2}/final.contigs.fa {output}
		rm -r working_dir/megahit2_cir{wildcards.p1}.cie{wildcards.p2}		
		"""


rule sanity_check_asms_group1:
	input:
		"working_dir/asm.group1.cir{p1}.cie{p2}.raw.fa",
		"results/kmer_diff.group1_minus_group2.cir{p1}.cie{p2}.kmc_suf"
	output:
		"working_dir/asm.group1.cir{p1}.cie{p2}.credible.fa"
	threads: 6
	shell:
		"""
		# filter ASMs for only those contigs containing the kmers that were used to retreive the reads / sanity-check
		kmc_tools -t{threads} filter results/kmer_diff.group1_minus_group2.cir{wildcards.p1}.cie{wildcards.p2} -ci1 working_dir/asm.group1.cir{wildcards.p1}.cie{wildcards.p2}.raw.fa -fa -ci1 working_dir/asm.group1.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa
		
		echo group_1_specific {wildcards.p1} {wildcards.p2} $( grep ">" working_dir/asm.group1.cir{wildcards.p1}.cie{wildcards.p2}.raw.fa | wc -l ) $( grep ">" working_dir/asm.group1.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa | wc -l ) >> working_dir/asm_sanity_check_report.txt	
		"""

rule sanity_check_asms_group2:
	input:
		"working_dir/asm.group2.cir{p1}.cie{p2}.raw.fa",
		"results/kmer_diff.group2_minus_group1.cir{p1}.cie{p2}.kmc_suf"
	output:
		"working_dir/asm.group2.cir{p1}.cie{p2}.credible.fa"
	threads: 6
	shell:
		"""
		# filter ASMs for only those contigs containing the kmers that were used to retreive the reads / sanity-check
		kmc_tools -t{threads} filter results/kmer_diff.group2_minus_group1.cir{wildcards.p1}.cie{wildcards.p2} -ci1 working_dir/asm.group2.cir{wildcards.p1}.cie{wildcards.p2}.raw.fa -fa -ci1 working_dir/asm.group2.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa
		
		echo group_2_specific {wildcards.p1} {wildcards.p2} $( grep ">" working_dir/asm.group2.cir{wildcards.p1}.cie{wildcards.p2}.raw.fa | wc -l ) $( grep ">" working_dir/asm.group2.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa | wc -l ) >> working_dir/asm_sanity_check_report.txt	
		"""


rule rename_and_merge_fastas:
	#
	# this rule elicits iteration over the parameter combinations
	#
	input:
		expand("working_dir/asm.{group}.cir{param_comb[0]}.cie{param_comb[1]}.credible.fa", group=groups, param_comb = ci_param_combs)
	output:
		temp( "working_dir/all_credible_contigs.fa" )
	shell:
		"""
		echo "" > {output}
		for f in {input} ; do
 			fname=$(echo $f | sed 's/working_dir\///g')
			cat $f | sed "s/flag.*/$fname/g" | sed 's/ /-/g' >> {output}
		done	
		"""
		
rule blastn_nt:
	input:
		"working_dir/all_credible_contigs.fa"
	output:
		"working_dir/all_credible_contigs.fa.blastn_nt.txt"
	threads: 12
	params:
		ntpath = config["nt_path"]
	shell:
		"""
		blastn -db {params.ntpath} -query {input} -evalue 1e-5 -max_hsps 1 -max_target_seqs 20 -num_threads {threads} -outfmt '6 qseqid sseqid qlen slen length pident evalue bitscore staxid' > {output}
		"""	
		
rule classify_taxonomy:
	input:
		"working_dir/all_credible_contigs.fa.blastn_nt.txt"
	output:
		"working_dir/all_credible_contigs.fa.taxonomy_report.txt"
	params:
		p = config["taxid_list_Viridiplantae"],
		f = config["taxid_list_Fungi"],
		m = config["taxid_list_Metazoa"],
		b = config["taxid_list_Bacteria"]
	run:

		plantset = []
		with open(params.p, "r") as INF:
			for line in INF:
				plantset.append(line.strip("\n"))


		bactset = []
		with open(params.b, "r") as INF:
			for line in INF:
				bactset.append(line.strip("\n"))


		fungiset = []
		with open(params.f, "r") as INF:
			for line in INF:
				fungiset.append(line.strip("\n"))


		metazoaset = []
		with open(params.b, "r") as INF:
			for line in INF:
				metazoaset.append(line.strip("\n"))

		plantset = set(plantset)
		fungiset = set(fungiset)
		bactset = set(bactset)
		metazoaset = set(metazoaset)

		plantcnt = 0
		bactcnt = 0
		fungicnt = 0
		metazoacnt = 0
		othercnt = 0
		humancnt = 0
		contig_score_dict = {}
		with open("working_dir/all_credible_contigs.fa.blastn_nt.txt","r") as INF:
			for line in INF:
				fields = line.strip("\n").split("\t")
				if len(fields) >= 9:
					contig = fields[0].split(":")[0]
					contig_score_dict[contig] = [0,0,0,0,0,0]


		with open("working_dir/all_credible_contigs.fa.blastn_nt.txt","r") as INF:
			for line in INF:
				fields = line.strip("\n").split("\t")
				if len(fields) >= 9:
					contig = fields[0].split(":")[0]
					bitscore = float(fields[7])
					taxid = fields[8]
					if taxid in plantset:
						contig_score_dict[contig][0] += bitscore
					elif taxid in bactset:
						contig_score_dict[contig][1] += bitscore
					elif taxid in fungiset:
						contig_score_dict[contig][2] += bitscore
					elif taxid == "9606": # Human
						contig_score_dict[contig][3] += bitscore
					elif taxid in metazoaset:
						contig_score_dict[contig][4] += bitscore
					else:
						contig_score_dict[contig][5] += bitscore




		classes = ["Viridiplantae","Bacteria","Fungi","human","Metazoa","other"]
		for k,v in contig_score_dict.items():
			contig_score_dict[k] = v + [ classes[ v.index(max(v)) ] ]


		# export a report:
		outlines =  ["\t".join(["contig","Viridiplantae","Bacteria","Fungi","human","Metazoa","other","classification_result"])]
		for k,v in contig_score_dict.items():
			outlines.append( k + "\t" + "\t".join( [str(x) for x in v] ) )


		with open("working_dir/all_credible_contigs.fa.taxonomy_report.txt", "w") as O:
			O.write("\n".join(outlines) + "\n")
		##################


rule filter_for_plant_contigs_and_split:
	input:
		"working_dir/all_credible_contigs.fa.taxonomy_report.txt",
		"working_dir/asm.{group}.cir{p1}.cie{p2}.credible.fa"
	output:
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.fa" )
	shell:
		"""
		set +o pipefail;
		tail -n +2 working_dir/all_credible_contigs.fa.taxonomy_report.txt | grep Viridiplantae | grep asm.{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa | cut -f1 | tr "-" " " | cut -d " " -f1 > working_dir/asm.{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}.TARGETS
		
		seqtk subseq working_dir/asm.{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}.credible.fa working_dir/asm.{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}.TARGETS > {output}
		rm working_dir/asm.{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}.TARGETS	
		"""


rule rename_contigs_by_TAIR_hit:
	input:
		"working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.fa"
	output:
		"working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa"
	params:
		ATpath = config["TAIR_path"]
	threads: 6
	run:
		shell("""
		blastx -db {params.ATpath} -query {input} -num_threads {threads} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen length pident evalue' | awk 'seen[$1]++ == 0' > {input}.blastx_AT.txt
		""")
			
		namedict = {}
		with open(str(input) + ".blastx_AT.txt", "r") as I:
			for line in I:
				fields = line.strip().split()
				namedict[fields[0]] = fields[1]
		
		outlines = []		
		with open(str(input), "r") as I:
			for line in I:
				if line.startswith(">"):				
					oname = line.strip(">").strip("\n").split()[0]
					try:
						outlines.append(">" + oname + "__" + namedict[oname])
					except KeyError:
						outlines.append(">" + oname + "__noAthit")
				else:
					outlines.append(line.strip())
					
		with open(str(output), "w") as O:
			O.write("\n".join(outlines)+"\n")






	
rule score_plant_contig_stats:
	input:
		expand("working_dir/asm.{group}.cir{param_comb[0]}.cie{param_comb[1]}.credible.plant.fa", group=groups, param_comb = ci_param_combs)
	output:
		"results/report_plant_contig_stats.txt"			
	shell:
		"""
		for f in {input} ; do
			echo $f $( seqtk comp $f | awk 'BEGIN{{sum=0}} {{sum+=$2;}} END{{print sum;}}' ) >> {output}
		done
		"""

		


rule get_cov_ratio:
	input:
		bam1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam",
		bai1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam.bai",
		bam2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam",
		bai2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam.bai"
	output:
		temp( "results/asm.{group}.cir{p1}.cie{p2}.read_coverage_proportions.txt" )
	shell:
		"""	
		samtools idxstats {input.bam1} | head -n -1 > {input.bam1}.idxstats.txt
		samtools idxstats {input.bam2} | head -n -1 > {input.bam2}.idxstats.txt
		
		paste <(cut -f3 {input.bam1}.idxstats.txt) <(cut -f3 {input.bam2}.idxstats.txt) | awk '{{if(($1+$2)>0) {{print ($1+$2)"\t"($1/($1+$2)) }} else {{print "0\tNA"}} }}' > propgr1_{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2} 
		
		paste <(cut -f1,2 {input.bam1}.idxstats.txt) propgr1_{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2} > {output}
		rm propgr1_{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2}
		rm {input.bam1}.idxstats.txt {input.bam2}.idxstats.txt
		"""
			

rule samtools_index:
	input:
		bam1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam",
		bam2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam"
	output:
		bai1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam.bai",
		bai2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam.bai"
	shell:
		"""
		samtools index {input.bam1}
		samtools index {input.bam2}		
		"""


rule samtools_sort:
	input:
		bam1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.bam",
		bam2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.bam"
	output:
		bam1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam",
		bam2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam"
	shell:
		"""
		samtools sort -T mapped_reads/{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2} -O bam {input.bam1} > {output.bam1}
		samtools sort -T mapped_reads/{wildcards.group}.cir{wildcards.p1}.cie{wildcards.p2} -O bam {input.bam2} > {output.bam2}
		"""


rule bwa_map:
	input:
		fa="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa",
		gidx1="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.bwt",
		gidx2="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.amb",
		gidx3="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.ann",
		gidx4="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.pac",
		gidx5="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.sa",
		r11="working_dir/group1.1.fq.gz",
		r21="working_dir/group2.1.fq.gz"
	output:
		bam1=temp( "mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.bam"),
		bam2=temp( "mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.bam")
	threads: 14
	shell:
		"""
		set +o pipefail;
		
		# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
		# -F 256 == -F 0x0100 == NOT not primary alignment
		# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
		# -F 2048 == -F 0x800 == NOT supplementary alignment
		# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
		# filtering alignments to be "properly paired": -f 2
		bwa mem -t {threads} -a {input.fa} working_dir/group1.1.fq.gz working_dir/group1.2.fq.gz -R "@RG\\tID:group1\\tSM:group1\\tPL:Illumina" | samtools view -F 2304 -f 2 -b -@ 2 - > {output.bam1}
		bwa mem -t {threads} -a {input.fa} working_dir/group2.1.fq.gz working_dir/group2.2.fq.gz -R "@RG\\tID:group2\\tSM:group2\\tPL:Illumina" | samtools view -F 2304 -f 2 -b -@ 2 - > {output.bam2}  
		"""

			
rule bwa_idx:
	input:
	   "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa"
	output:
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.bwt"),
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.amb"),
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.ann"),
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.pac"),
		temp( "working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.sa")
	shell:
		"""
		if [[ ! $( grep ">" {input} ) =~ "|" ]]; then
			bwa index {input}
		else
			echo "refusing to run, fasta headers contain pipe '|' character, dying"
		fi
		"""
		



		
rule pseudo_phase_gametologs:
	input:
		"results_varcalls/asm.{group}.cir{p1}.cie{p2}.vcf.gz"
	output:
		temp( "results/asm.{group}.cir{p1}.cie{p2}.XY_and_ZW_stats.txt" )
	shell:
		"""
		zcat {input} | python scripts/pseudo_gametolog_stats.py > {output}.temp
		
		tail -n +3 {output}.temp > {output}
		rm {output}.temp
		"""

	
rule bcftools_call_variants:
	input:
		ref="working_dir/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa",
		bam1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam",
		bai1="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group1.sorted.bam.bai",
		bam2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam",
		bai2="mapped_reads/asm.{group}.cir{p1}.cie{p2}.credible.plant.AT.fa.reads_group2.sorted.bam.bai"
	output:
		"results_varcalls/asm.{group}.cir{p1}.cie{p2}.vcf.gz"
	threads: 2
	shell:
		"""
		set +o pipefail;
		
		# --max-depth 250: use at most 250 reads per input BAM file, apparently these are sampled RANDOMLY!?
		# --min-MQ 15: minimum mapping quality of an alignment, otherwise skip
		# --no-BAQ : do NOT re-calculate mapping quality (which involves re-aligning). Instead, will use MAPQ as stored in BAM file.
		# --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
		# -a DP,AD : annotate for each sample genotype DP = total depth and AD = allelic depth, count of reads separately for each allele.
		
		bcftools mpileup -Ou -f {input.ref} {input.bam1} {input.bam2} --max-depth 10000 --min-MQ 20 --min-BQ 15 --no-BAQ -a DP,AD | bcftools call -m --skip-variants indels -Ov -v | bgzip -c > {output}
		"""


rule combine_cov_stats_and_gametolog_stats:
	input:
		gametolog_stats="results/asm.{group}.cir{p1}.cie{p2}.XY_and_ZW_stats.txt",
		covstats="results/asm.{group}.cir{p1}.cie{p2}.read_coverage_proportions.txt"
	output:
		"results/asm.{group}.cir{p1}.cie{p2}.stats_per_contig.txt"
	shell:
		"""
		set +o pipefail;
		echo -e "contig\tlength\taligned_reads_total\tprop_reads_group1\tSNPs_XY_like\tSNPs_ZW_like" > {output}
		python scripts/merge_tabsep_on_column.py {input.covstats} 0 {input.gametolog_stats} 0 | cut -f 1,2,3,4,6,7 >> {output}
		
		"""




