rule wham_call:
	input:
		unpack(get_bams),
		ref = "results/{refGenome}/data/genome/{refGenome}.fna",
	output:
		wham_vcf=temp("results/{refGenome}/SV/wham/{sample}.raw.vcf"),
	log:
		"logs/{refGenome}/SV/wham/{sample}.log"
	conda:
		"../envs/wham.yaml"
	params:
		contigs = lambda wc: read_contig_file(config["include_contigs"])
	threads: 10 #TODO Make this dynamic based on user provided threads
	shell:
		"""
		whamg -x {threads} -c {params.contigs} -a {input.ref} -f {input.bam} > {output.wham_vcf} 2> {log}
		"""
rule fix_whamcall_header:
	"""
	Add contig ID and length to wham vcf file
	"""
	input:
		wham_vcf = "results/{refGenome}/SV/wham/{sample}.raw.vcf",
		fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
	output:
		"results/{refGenome}/SV/wham/{sample}.vcf",
	conda:
		"../envs/wham.yaml"
	shell:
		"""
		bcftools reheader -f {input.fai} -o {output} {input.wham_vcf}
		"""