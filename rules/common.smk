import sys
import os
import tempfile
import random
import string
import statistics
import json
from pathlib import Path
from collections import defaultdict
import pandas as pd
try:
    from importlib.metadata import version
except ImportError:
    # Fallback for Python < 3.8
    from importlib_metadata import version

def parse_version(v):
    """Simple version parsing function"""
    return tuple(map(int, v.split('.')))

# Can't be less than 7 cuz of min version in snakefile
SNAKEMAKE_VERSION = 8 if parse_version(snakemake.__version__) >= parse_version("8.0.0") else 7
logger.warning(f"svArcher: Using Snakemake {snakemake.__version__}")
if SNAKEMAKE_VERSION >= 8:
    DEFAULT_STORAGE_PREFIX = StorageSettings.default_storage_prefix if StorageSettings.default_storage_prefix is not None else ""
else:
    # backwards compatibility w/ snakemake <= 7
    DEFAULT_STORAGE_PREFIX = workflow.default_remote_prefix
    if config.get("remote_reads", False):
        from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
        GS = GSRemoteProvider()
        GS_READS_PREFIX = config['remote_reads_prefix']

def parse_sample_sheet(config):
    """Parse sample sheet CSV file."""
    return pd.read_csv(config["samples"])

def has_custom_reference(refGenome):
    """Check if a custom reference path exists for the given refGenome.
    
    Validates that all samples with the same refGenome have consistent refPath values.
    Returns True if exactly one distinct non-empty refPath exists.
    Returns False if no non-empty refPath exists.
    Raises WorkflowError if multiple distinct non-empty refPath values exist.
    """
    if not config.get("use_custom_reference", False) or "refPath" not in samples.columns:
        return False
    
    # Get all rows with this refGenome
    ref_rows = samples[samples["refGenome"] == refGenome]
    if ref_rows.empty:
        return False
    
    # Collect all refPath values and normalize (treat NaN/empty/whitespace as missing)
    ref_paths = []
    for _, row in ref_rows.iterrows():
        ref_path = row["refPath"]
        if pd.notna(ref_path) and str(ref_path).strip():
            ref_paths.append(str(ref_path).strip())
    
    # Check for consistency
    unique_ref_paths = list(set(ref_paths))
    
    if len(unique_ref_paths) == 0:
        # No non-empty refPath exists
        return False
    elif len(unique_ref_paths) == 1:
        # Exactly one distinct non-empty refPath exists
        return True
    else:
        # Multiple distinct non-empty refPath values exist
        raise WorkflowError(
            f"Inconsistent refPath values found for refGenome '{refGenome}': "
            f"{unique_ref_paths}. All samples with the same refGenome must have "
            f"the same refPath value or be empty/missing."
        )

def get_custom_reference_path(refGenome):
    """Get the custom reference path for the given refGenome.
    
    Returns the consistent refPath value if one exists, None otherwise.
    """
    if has_custom_reference(refGenome):
        # Get all rows with this refGenome
        ref_rows = samples[samples["refGenome"] == refGenome]
        
        # Find the first non-empty refPath (they're all the same due to validation in has_custom_reference)
        for _, row in ref_rows.iterrows():
            ref_path = row["refPath"]
            if pd.notna(ref_path) and str(ref_path).strip():
                return str(ref_path).strip()
    return None

def get_bams(wc):
    """Get BAM files for SV calling - use final deduplicated BAMs if mark_duplicates is enabled."""
    out = {"bam": None, "bai": None}
    if config.get("mark_duplicates", True):
        out["bam"] = "results/{refGenome}/bams/{sample}_final.bam"
        out["bai"] = "results/{refGenome}/bams/{sample}_final.bam.bai"
        return out
    else:
        return dedup_input(wc)

def dedup_input(wc):
    """Get input for deduplication step - either single BAM or merged BAM."""
    runs = samples.loc[samples["BioSample"] == wc.sample]["Run"].tolist()

    if len(runs) == 1:
        bam = expand("results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam", run=runs)
        bai = expand("results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam.bai", run=runs)
    else:
        bam = "results/{refGenome}/bams/postMerge/{sample}.bam"
        bai = "results/{refGenome}/bams/postMerge/{sample}.bam.bai"
    return {"bam": bam, "bai": bai}

def merge_bams_input(wc):
    """Get all BAM files for a sample that need to be merged."""
    return expand(
        "results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam",
        run=samples.loc[samples["BioSample"] == wc.sample]["Run"].tolist(),
    )

def get_reads(wc):
    """Returns local read files if present. Defaults to SRR if no local reads in sample sheet."""
    row = samples.loc[samples["Run"] == wc.run]
    r1 = f"results/data/fastq/{wc.refGenome}/{wc.sample}/{wc.run}_1.fastq.gz"
    r2 = f"results/data/fastq/{wc.refGenome}/{wc.sample}/{wc.run}_2.fastq.gz"
    if "fq1" in samples.columns and "fq2" in samples.columns:
        if row["fq1"].notnull().any() and row["fq2"].notnull().any():
            r1 = row.fq1.item()
            r2 = row.fq2.item()
            if config.get("remote_reads", False):
                if SNAKEMAKE_VERSION >= 8:
                    # remote read path must have full remote prefix, eg: gs://reads_bucket/sample1/...
                    # depends on snakemake>8 to figure out proper remote provider from prefix using storage()
                    return {"r1": storage(r1), "r2": storage(r2)}
                else:
                    return get_remote_reads(wc)
            if os.path.exists(row.fq1.item()) and os.path.exists(row.fq2.item()):
                return {"r1": r1, "r2": r2}
            else:
                raise WorkflowError(
                    f"fq1 and fq2 specified for {wc.sample}, but files were not found."
                )
        else:
            # this allows mixed srr and user-specified paths for reads
            return {"r1": r1, "r2": r2}
    else:
        return {"r1": r1, "r2": r2}

def get_reads_fastp(wc):
    """Get reads for fastp - either sorted or original reads."""
    if config.get("sort_reads", False):
        return {"r1":"results/{refGenome}/sorted_reads/{sample}/{run}_1.fastq.gz", "r2":"results/{refGenome}/sorted_reads/{sample}/{run}_2.fastq.gz"}
    else:
        return get_reads(wc)

def get_remote_reads(wildcards):
    """Use this for reads on a different remote bucket than the default. For backwards compatibility."""
    row = samples.loc[samples["Run"] == wildcards.run]
    r1 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq1.item()))
    r2 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq2.item()))
    return {"r1": r1, "r2": r2}

def get_read_group(wc):
    """Denote sample name and library_id in read group."""
    libname = samples.loc[samples["Run"] == wc.run]["LibraryName"].tolist()[0]
    return r"'@RG\tID:{lib}\tSM:{sample}\tLB:{lib}\tPL:ILLUMINA'".format(
        sample=wc.sample, lib=libname
    )

def get_big_temp(wildcards):
    """Sets a temp dir for rules that need more temp space that is typical on some cluster environments. Defaults to system temp dir."""
    if config.get("bigtmp"):
        if config["bigtmp"].endswith("/"):
            return (
                config["bigtmp"]
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
        else:
            return (
                config["bigtmp"]
                + "/"
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
    else:
        return tempfile.gettempdir()

def read_contig_file(contig_file: str) -> list[str]:
    """Read contig file for wham SV caller."""
    try:
        with open(contig_file, "r") as file:
            contigs = [line.strip() for line in file if line.strip()]
        if not contigs:
            logger.warning(f"Contig file {contig_file} is empty")
        return contigs
    except FileNotFoundError:
        logger.error(f"Contig file {contig_file} not found")
        raise
    except Exception as e:
        logger.error(f"Error reading contig file {contig_file}: {e}")
        raise

def get_sv_caller_outputs(wildcards):
    """Get all SV caller outputs that must be completed before post-processing"""
    return expand(
        "results/{refGenome}/SV/{method}/{sample}.vcf",
        refGenome=REFGENOME,
        method=["delly", "lumpy", "wham"],
        sample=samples["BioSample"].unique().tolist()
    )

def merge_config(wildcards):
    """
    Pass merge parameter to SURVIVOR merge from config file.
    """
    try:
        if config.get("sv_merge")>=1 and config.get("sv_merge")<=3:
            return config.get("sv_merge")
        else:
            raise WorkflowError(f"Invalid sv_merge value: {config.get('sv_merge')}. Must be 1, 2, or 3.")
    except TypeError:
        sv_merge_value = config.get('sv_merge')
        logger.error(f"Invalid sv_merge value provided: '{sv_merge_value}'. Must be an integer 1, 2, or 3, not {type(sv_merge_value).__name__} type")
        raise

def svArcher_output(wildcards):
    """Define all final outputs for svArcher pipeline."""
    output = []
    # Return SV calls from all methods
    output.extend(get_sv_caller_outputs(wildcards))
    # Return merged SV calls per sample
    output.extend(expand("results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf", refGenome=REFGENOME, sample=samples["BioSample"].unique().tolist()))
    # Return processed SV calls per sample
    output.extend(expand("results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf", refGenome=REFGENOME, sample=samples["BioSample"].unique().tolist()))
    # Return an all sample merged vcf
    output.extend(expand("results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf", refGenome=REFGENOME))
    output.extend(expand("results/{refGenome}/SV/postprocess/processed/all_samples_final.vcf", refGenome=REFGENOME))
    output.extend(expand("results/{refGenome}/SV/sv_metadata/metadata.tsv", refGenome=REFGENOME))
    if config.get("svBenchmark", False):
        output.extend(expand("results/{refGenome}/SV/benchmark/{svType}/{svType}_METRICS/summary.json", refGenome=REFGENOME, svType=["DEL", "INV", "DUP"]))
    return output