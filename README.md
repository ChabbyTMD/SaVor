# SaVor: A Structural Variant Calling and Benchmarking Workflow for Short-Read Sequence Data.

<div align="center">
<img src="img/savor_logo.png" alt="SaVor Logo" width="350" height="350">
</div>

*Savor* is a standalone reproducible Snakemake pipeline designed to call structural variants (SVs) from short read genomic sequencing data. The pipeline performs:

1. **Data ingestion**: FASTQ download/processing and quality control
2. **Read alignment**: BWA-MEM2 alignment with duplicate marking
3. **Structural variant calling**: Using an ensemble of three SV callers (DELLY, Lumpy, and WHAM)
4. **Post-processing**: Merging, filtering, and consensus calling with SURVIVOR

## Features

- **Complete standalone pipeline**: No dependencies on external workflows
- **Ensemble SV calling**: Uses DELLY, Lumpy, and WHAM for robust SV detection
- **Automated data ingestion**: Downloads FASTQ files from SRA or processes local files
- **Quality control**: Integrated fastp for read filtering and quality assessment
- **Consensus calling**: SURVIVOR-based merging requiring agreement from multiple callers
- **Benchmarking support**: Optional comparison against truth sets using Truvari
- **Conda environments**: All dependencies managed through conda/mamba
- **Flexible input**: Support for FASTQ files (SRA or local) or user-provided BAM files

## Directory Structure

```{text}
SaVor_Standalone/
├── Snakefile                 # Main workflow file
├── config/
│   ├── config.yaml          # Configuration file
│   ├── samples.csv          # Sample metadata
│   └── include_contigs.csv  # Contigs for WHAM caller
├── rules/                   # Snakemake rule files
│   ├── common.smk           # Common functions
│   ├── fastq_ingestion.smk  # Data download and QC
│   ├── alignment.smk        # Read alignment
│   ├── lumpy.smk            # Lumpy SV caller
│   ├── delly.smk            # DELLY SV caller
│   ├── wham.smk             # WHAM SV caller
│   ├── svcallprocess.smk    # Post-processing
│   └── benchmark.smk        # Benchmarking (optional)
└── envs/                    # Conda environment files
    ├── fastq2bam.yml
    ├── sambamba.yml
    ├── lumpy.yaml
    ├── delly.yaml
    ├── wham.yaml
    ├── survivor.yaml
    └── truvari.yaml
```

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥9.0)

### Setup

1. Clone or download this repository
2. Activate the existing snakemake environment:

```bash
conda activate snakemake
```

**Note**: This pipeline assumes you have a conda environment named "snakemake" with Snakemake installed. If you don't have it, create one:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake>=9.0
conda activate snakemake
```

## Usage

### Basic Run

```bash
# Activate environment
conda activate snakemake

# Dry run to check workflow
snakemake -n

# Run with at least 8 cores
snakemake --cores <num_of_cores> --use-conda

# Run on cluster (SLURM example)
snakemake --workflow-profile workflow-profile/cluster-generic-slurm/
```

### User Guide

For detailed guidance on workflow setup, please consult the [user guide](docs/README.md).
### Output Files

The pipeline generates several key outputs:

- **Individual SV calls**: `results/{refGenome}/SV/{caller}/{sample}.vcf`
- **Merged per-sample calls**: `results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf`
- **Processed calls**: `results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf`
- **Final consensus**: `results/{refGenome}/SV/postprocess/processed/all_samples_final.vcf`
- **SV Call Metadata**: `results/{refGenome}/SV/sv_metadata/metadata.tsv`
- **BAM files**: `results/{refGenome}/bams/{sample}_final.bam`

### Benchmarking (Optional)

To enable benchmarking against truth sets:

1. Set `svBenchmark: true` in config
2. Provide paths to truth VCF files for deletions, duplications, and inversions
3. Results will be in `results/{refGenome}/SV/benchmark/{svType}/`

## Workflow Steps

1. **Data Ingestion**:

   - Download FASTQ files from SRA (`get_fastq_pe`)
   - Optional read sorting (`sort_reads`)
   - Quality control with fastp (`fastp`)
2. **Reference Preparation**:

   - Download reference genome if needed (`download_reference`)
   - Create BWA and samtools indices (`index_reference`)
3. **Read Alignment**:

   - BWA-MEM2 alignment (`bwa_map`)
   - BAM merging if multiple runs per sample (`merge_bams`)
   - Duplicate marking with sambamba (`dedup`)
4. **Structural Variant Calling**:

   - DELLY calling (`delly_call`, `delly_vcf`)
   - Lumpy calling (`discordant_extract`, `split_read_extract`, `lumpy_call`)
   - WHAM calling (`wham_call`)
5. **Post-processing**:

   - VCF sorting (`sort_vcfs`)
   - Per-sample merging (`sample_sv_call_merge`)
   - Filtering (`filter_sv_calls`)
   - All-sample merging (`all_sample_merge`)
   - Final processing (`post_merge_process`)

## Resource Requirements

- **Memory**: 8-16GB per sample for alignment, 4-8GB for SV calling
- **Storage**: ~10-50GB per sample depending on coverage
- **Time**: 2-8 hours per sample depending on coverage and resources

## Workflow Tutorial
Consult the documentation directory for an in-depth [tutorial](docs/savor_tutorial.md).

## Citation

If you use SaVor, please cite in addition to all the tools implemented:

```{text}
SaVor - A Reproducible Structural Variant Calling and Benchmarking Platform from Short-Read Data
[In Prep]
```

## License

[TBD]

## Support

For questions and issues, please submit an issue on the repo.
