# Checkpoint rule to ensure all SV caller outputs are available before processing
checkpoint ensure_sv_calls_complete:
    input:
        get_sv_caller_outputs
    output:
        touch("results/{refGenome}/SV/postprocess/.sv_calls_complete")
    shell:
        """
        echo "All SV caller outputs are available, proceeding with post-processing"
        """

# Sort vcf files from delly, lumpy and wham
rule sort_vcfs:
	input:
		delly_vcf="results/{refGenome}/SV/delly/{sample}.vcf",
		lumpy_vcf="results/{refGenome}/SV/lumpy/{sample}.vcf",
		wham_vcf="results/{refGenome}/SV/wham/{sample}.vcf",
		checkpoint="results/{refGenome}/SV/postprocess/.sv_calls_complete",
	output:
		delly_vcf_sorted="results/{refGenome}/SV/delly/{sample}.sorted.vcf",
		lumpy_vcf_sorted="results/{refGenome}/SV/lumpy/{sample}.sorted.vcf",
		wham_vcf_sorted="results/{refGenome}/SV/wham/{sample}.sorted.vcf",
	conda:
		"../envs/wham.yaml"
	log:
		"logs/{refGenome}/SV/postprocess/sort_{sample}.log"
	shell:
		"""
		bcftools sort {input.delly_vcf} -o {output.delly_vcf_sorted} 2> {log}
		bcftools sort {input.lumpy_vcf} -o {output.lumpy_vcf_sorted} 2>> {log}
		bcftools sort {input.wham_vcf} -o {output.wham_vcf_sorted} 2>> {log}
		"""
# Rule to merge SVS from all methods
checkpoint sample_sv_call_merge:
	input:
		delly_vcf="results/{refGenome}/SV/delly/{sample}.sorted.vcf",
		lumpy_vcf="results/{refGenome}/SV/lumpy/{sample}.sorted.vcf",
		wham_vcf="results/{refGenome}/SV/wham/{sample}.sorted.vcf",
	output:
		merged_vcf="results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf",
	conda:
		"../envs/survivor.yaml"
	log:
		"logs/{refGenome}/SV/postprocess/{sample}.log"
	params:
		sample_file = "results/{refGenome}/SV/postprocess/raw_merge/{sample}.txt",
		merge_param=lambda wc: config.get("sv_merge", 1),
	shell:
		"""
		# Create a file with the paths to the vcf files
		# This is needed for survivor merge
		ls {input.delly_vcf} {input.lumpy_vcf} {input.wham_vcf} > {params.sample_file}
		SURVIVOR merge {params.sample_file} 1000 {params.merge_param} 1 1 0 50 {output.merged_vcf} 2> {log}
		"""

checkpoint filter_sv_calls:
    input:
        merged_vcf="results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf",
    output:
        filtered_vcf="results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf",

    conda:
        "../envs/survivor.yaml"
    log:
        "logs/{refGenome}/SV/postprocess/{sample}.filter.log"
    shell:
        """
		# Filter out SV calls with length > 10kb
        bcftools view -i 'abs(SVLEN) <= 10000' {input.merged_vcf} -o {output.filtered_vcf} 2> {log}
        """
# Rule to get paths of all processed VCF files
rule get_processed_vcf_paths:
    input:
        processed_vcfs=lambda wc: expand(
            "results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf",
            refGenome=wc.refGenome,
            sample=samples["BioSample"].unique().tolist()
        ),
    output:
        path_list="results/{refGenome}/SV/postprocess/processed/all_samples.txt",
    run:
        import os
        # Create directory if needed
        os.makedirs(os.path.dirname(output.path_list), exist_ok=True)
        # Write sorted VCF paths to list
        with open(output.path_list, "w") as out_f:
            for p in sorted(input.processed_vcfs):
                out_f.write(p + "\n")

# Update rules to use the default value
rule all_sample_merge:
    input:
        # Use the file containing the list of VCF paths as input
        path_list="results/{refGenome}/SV/postprocess/processed/all_samples.txt",
    output:
        merged_vcf="results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf",
    conda:
        "../envs/survivor.yaml",
    log:
        "logs/{refGenome}/SV/postprocess/all_samples.log",
    shell:
        """
        # SURVIVOR merge takes the file with paths as the first argument
        SURVIVOR merge {input.path_list} 1000 1 1 1 0 50 {output.merged_vcf} 2> {log}
        """
rule alt_metadata_save:
    input:
        vcf_file="results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf",
    output:
        metadata_file="results/{refGenome}/SV/sv_metadata/metadata.tsv"
    shell:
        """
        mkdir -p $(dirname {output.metadata_file})
        grep -v "^#" {input.vcf_file} | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1, $3, $5}}' > {output.metadata_file}
        """
rule post_merge_process:
	input:
		merged_vcf="results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf",
	output:
		final_vcf="results/{refGenome}/SV/postprocess/processed/all_samples_final.vcf",
	shell:
		"""
		awk '
		BEGIN {{ FS = OFS = "\t" }}
		# Print header lines unchanged
		/^#/ {{ print; next }}
		# Process records where the 5th field is NOT <DEL>, <DUP>, or <INV>
		$5 != "<DEL>" && $5 != "<DUP>" && $5 != "<INV>" {{
			# Find SVTYPE=value in the INFO field ($8)
			if ($8 ~ /SVTYPE=[^;]+/) {{
				# Extract the SVTYPE value
				split($8, info_fields, ";")
				for (i in info_fields) {{
					if (info_fields[i] ~ /^SVTYPE=/) {{
						svtype_value = substr(info_fields[i], 8)
						break
					}}
				}}
				# Replace ALT field ($5) with the extracted SVTYPE value enclosed in <>
				$5 = "<" svtype_value ">"
			}}
			# Print the modified line
			print
			next
		}}
		# Print other records unchanged
		{{ print }}
		' {input.merged_vcf} > {output.final_vcf}
		"""