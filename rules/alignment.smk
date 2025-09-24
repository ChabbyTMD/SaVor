rule bwa_map:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        r1 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_1.fastq.gz",
        r2 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_2.fastq.gz",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["0123", "pac", "bwt.2bit.64", "ann", "amb", "fai"]),
    output: 
        bam = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam"),
        bai = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam.bai"),
    params:
        rg = get_read_group
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/bwa_mem/{sample}/{run}.txt"
    benchmark:
        "benchmarks/{refGenome}/bwa_mem/{sample}_{run}.txt"
    threads: 8
    shell:
        "bwa-mem2 mem -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam} - && samtools index {output.bam} {output.bai}"

rule merge_bams:
    input:
        merge_bams_input
    output:
        bam = temp("results/{refGenome}/bams/postMerge/{sample}.bam"),
        bai = temp("results/{refGenome}/bams/postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_bams/{sample}.txt"
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam} > {log}"

rule dedup:
    input:
        unpack(dedup_input)
    output:
        dedupBam = "results/{refGenome}/bams/{sample}_final.bam",
        dedupBai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    conda:
        "../envs/sambamba.yml"
    log:
        "logs/{refGenome}/sambamba_dedup/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/sambamba_dedup/{sample}.txt"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"

rule download_reference:
    """Download reference genome if not provided by user."""
    output:
        ref="results/{refGenome}/data/genome/{refGenome}.fna",
    params:
        refGenome="{refGenome}"
    conda:
        "../envs/fastq2bam.yml"
    shell:
        """
        datasets download genome accession {params.refGenome} --filename {params.refGenome}.zip
        unzip -j {params.refGenome}.zip -d results/{wildcards.refGenome}/data/genome/
        mv results/{wildcards.refGenome}/data/genome/*.fna {output.ref}
        rm {params.refGenome}.zip
        """

rule index_reference:
    """Create BWA and samtools indices for reference genome."""
    input:
        ref="results/{refGenome}/data/genome/{refGenome}.fna"
    output:
        expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["0123", "pac", "bwt.2bit.64", "ann", "amb", "fai"])
    conda:
        "../envs/fastq2bam.yml"
    shell:
        """
        bwa-mem2 index {input.ref}
        samtools faidx {input.ref}
        """