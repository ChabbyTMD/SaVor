wildcard_constraints:
    sample="[a-zA-Z0-9_]+"

ruleorder:
    # fix_lumpycall_header > lumpy_call
    lumpy_call > fix_lumpycall_header

rule discordant_extract:
    input:
        unpack(get_bams),
    output:
        unsorted_discordant_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.discordant.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.discordant_extract.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools view -b -F 1294 {input.bam} > {output.unsorted_discordant_bam} 2> {log}
        """

rule split_read_extract:
    input:
        unpack(get_bams),
    output:
        unsorted_split_read_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.split_read.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.split_read_extract.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools view -h {input.bam} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb - \
            > {output.unsorted_split_read_bam} 2> {log}
        """

rule sort_splits_discordants:
    input:
        unsorted_discordant_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.discordant.bam",
        unsorted_split_read_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.split_read.bam",
    output:
        discordant_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.discordant.bam"),
        split_read_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.split_read.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.sort_splits_discordants.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools sort -o {output.discordant_bam} {input.unsorted_discordant_bam} 2> {log}
        samtools sort -o {output.split_read_bam} {input.unsorted_split_read_bam} 2>> {log}
        """

rule lumpy_call:
    input:
        unpack(get_bams),
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
        sorted_discordant_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.discordant.bam",
        sorted_split_read_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.split_read.bam",
    output:
        lumpy_vcf = temp("results/{refGenome}/SV/lumpy/{sample}.raw.vcf"),
    log:
        "logs/{refGenome}/SV/lumpy/{sample}.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        lumpyexpress \
            -B {input.bam} \
            -S {input.sorted_split_read_bam} \
            -D {input.sorted_discordant_bam} \
            -o {output.lumpy_vcf} 2> {log}
        """

rule fix_lumpycall_header:
    input:
        vcf = "results/{refGenome}/SV/lumpy/{sample}.raw.vcf",
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
    output:
        "results/{refGenome}/SV/lumpy/{sample}.vcf",
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        bcftools reheader -f {input.fai} -o {output} {input.vcf}
        """