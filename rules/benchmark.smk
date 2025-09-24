# ruleorder: sort_benchmark_vcfs > index_vcfs
# Rule to separate DEL DUPS and INV for benchmarking
rule generate_sv_benchmark_test_vcfs:
    input:
        filtered_vcf="results/{refGenome}/SV/postprocess/processed/all_samples_final.vcf",
    output:
        del_vcf=temp("results/{refGenome}/SV/benchmark/DEL/all_samples.del.vcf"),
        dup_vcf=temp("results/{refGenome}/SV/benchmark/DUP/all_samples.dup.vcf"),
        inv_vcf=temp("results/{refGenome}/SV/benchmark/INV/all_samples.inv.vcf"),
    conda:
        "../envs/survivor.yaml"
    shell:
        """
        bcftools view -i 'INFO/SVTYPE=="DEL"' {input.filtered_vcf} > {output.del_vcf}
        bcftools view -i 'INFO/SVTYPE=="DUP"' {input.filtered_vcf} > {output.dup_vcf}
        bcftools view -i 'INFO/SVTYPE=="INV"' {input.filtered_vcf} > {output.inv_vcf}
        """

rule sort_and_index_vcfs:
    input:
        del_vcf="results/{refGenome}/SV/benchmark/DEL/all_samples.del.vcf",
        dup_vcf="results/{refGenome}/SV/benchmark/DUP/all_samples.dup.vcf",
        inv_vcf="results/{refGenome}/SV/benchmark/INV/all_samples.inv.vcf",
    output:
        del_vcf_sorted="results/{refGenome}/SV/benchmark/DEL/all_samples.del.sorted.vcf.gz",
        dup_vcf_sorted="results/{refGenome}/SV/benchmark/DUP/all_samples.dup.sorted.vcf.gz",
        inv_vcf_sorted="results/{refGenome}/SV/benchmark/INV/all_samples.inv.sorted.vcf.gz",
    conda:
        "../envs/survivor.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/sort.index.log"
    shell:
        """
        bcftools sort {input.del_vcf} -Oz -o {output.del_vcf_sorted} 2> {log}
        bcftools sort {input.dup_vcf} -Oz -o {output.dup_vcf_sorted} 2>> {log}
        bcftools sort {input.inv_vcf} -Oz -o {output.inv_vcf_sorted} 2>> {log}
        tabix -p vcf {output.del_vcf_sorted} 2>> {log}
        tabix -p vcf {output.dup_vcf_sorted} 2>> {log}
        tabix -p vcf {output.inv_vcf_sorted} 2>> {log}
        """
# Rule to benchamrk SV calls against Arabidopsis Gotkay set

rule del_truvari_bench:
    input:
        del_vcf="results/{refGenome}/SV/benchmark/DEL/all_samples.del.sorted.vcf.gz",
        reference = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        truvari_report="results/{refGenome}/SV/benchmark/DEL/DEL_METRICS/summary.json",
    conda:
        "../envs/truvari.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/DEL/truvari_del.log"
    params:
        del_truth = config["deletions"],
        bench_dir = "results/{refGenome}/SV/benchmark/DEL/DEL_METRICS/",
    shell:
        """
        rm -rf {params.bench_dir}
        truvari bench -b {params.del_truth} -c {input.del_vcf} --reference {input.reference} -o {params.bench_dir} -p 0.6 -P 0.6 -r 1000 2> {log}
        """

rule inv_truvari_bench:
    input:
        inv_vcf="results/{refGenome}/SV/benchmark/INV/all_samples.inv.sorted.vcf.gz",
        reference = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        truvari_report="results/{refGenome}/SV/benchmark/INV/INV_METRICS/summary.json",
    conda:
        "../envs/truvari.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/INV/truvari_inv.log"
    params:
        inv_truth = config["inversions"],
        bench_dir = "results/{refGenome}/SV/benchmark/INV/INV_METRICS/",
    shell:
        """
        rm -rf {params.bench_dir}
        truvari bench -b {params.inv_truth} -c {input.inv_vcf} --reference {input.reference} -o {params.bench_dir} -p 0.6 -P 0.6 -r 1000 2> {log}
        """

rule dup_truvari_bench:
    input:
        dup_vcf="results/{refGenome}/SV/benchmark/DUP/all_samples.dup.sorted.vcf.gz",
        reference = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        truvari_report="results/{refGenome}/SV/benchmark/DUP/DUP_METRICS/summary.json",
    conda:
        "../envs/truvari.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/DUP/truvari_dup.log"
    params:
        dup_truth = config["duplications"],
        bench_dir = "results/{refGenome}/SV/benchmark/DUP/DUP_METRICS/",
    shell:
        """
        rm -rf {params.bench_dir}
        truvari bench -b {params.dup_truth} -c {input.dup_vcf} --reference {input.reference} -o {params.bench_dir} -p 0.6 -P 0.6 -r 1000 2> {log}
        """