# svArcher Standalone Pipeline

<div align="center">
<img src="https://raw.githubusercontent.com/ChabbyTMD/snpArcher/main/workflow/modules/svArcher/img/svArcher_Logo.png" alt="svArcher Logo" width="150" height="190">
</div>

*svArcher* is a standalone reproducible Snakemake pipeline designed to call structural variants (SVs) from short read genomic sequencing data. The pipeline performs:

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

## Directory Structure

```
svArcher_Standalone/
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
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥7.0.0)

### Setup
1. Clone or download this repository
2. Activate the existing snakemake environment:

```bash
conda activate snakemake
```

**Note**: This pipeline assumes you have a conda environment named "snakemake" with Snakemake installed. If you don't have it, create one:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake>=7.0.0
conda activate snakemake
```

## Configuration

### 1. Sample Metadata (`config/samples.csv`)

Create a CSV file with your sample information:

```csv
BioSample,Run,LibraryName,refGenome,fq1,fq2
sample1,SRR12345678,lib1,GCF_000001735.4,,
sample2,SRR12345679,lib2,GCF_000001735.4,,
```

**Required columns:**
- `BioSample`: Sample identifier
- `Run`: SRA run accession (e.g., SRR12345678) or unique identifier for local files
- `LibraryName`: Library name for read groups
- `refGenome`: Reference genome accession or identifier

**Optional columns:**
- `fq1`, `fq2`: Paths to local FASTQ files (leave empty for SRA download)

### 2. Reference Contigs (`config/include_contigs.csv`)

List the chromosomes/contigs to include in WHAM analysis:

```
Chr1
Chr2
Chr3
Chr4
Chr5
```

### 3. Main Configuration (`config/config.yaml`)

```yaml
# Required settings
samples: "config/samples.csv"
include_contigs: "config/include_contigs.csv"

# Optional settings
sort_reads: false
mark_duplicates: true
remote_reads: false
svBenchmark: false

# Benchmarking (if enabled)
deletions: "benchmark_truth/deletions.vcf.gz"
duplications: "benchmark_truth/duplications.vcf.gz"
inversions: "benchmark_truth/inversions.vcf.gz"
```

## Usage

### Basic Run

```bash
# Activate environment
conda activate snakemake

# Dry run to check workflow
snakemake -n

# Run with 8 cores
snakemake --cores 8 --use-conda

# Run on cluster (SLURM example)
snakemake --cores 100 --use-conda --cluster "sbatch -p partition -c {threads} --mem={resources.mem_mb}MB -t 24:00:00"
```

### Output Files

The pipeline generates several key outputs:

- **Individual SV calls**: `results/{refGenome}/SV/{caller}/{sample}.vcf`
- **Merged per-sample calls**: `results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf`
- **Processed calls**: `results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf`
- **Final consensus**: `results/{refGenome}/SV/postprocess/processed/all_samples_final.vcf`
- **Metadata**: `results/{refGenome}/SV/sv_metadata/metadata.tsv`
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

## Citation

If you use svArcher, please cite:

```
svArcher: A Structural Variant Calling Extension module for snpArcher
[Citation information to be added]
```

## License

[License information to be added]

## Support

For questions and issues, please visit the [GitHub repository](https://github.com/ChabbyTMD/snpArcher).