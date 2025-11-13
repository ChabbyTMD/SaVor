# SaVor User Guide

# 1. Prerequisites

## 1.1 Conda/Mamba

Ensure you have conda/mamba installed on the system you intend to run the workflow on. A detailed mamba installation guide can be found [here](https://www.howtogeek.com/how-to-set-up-a-development-environment-with-mamba/). For cluster environments, please contact your system administrator for conda/mamba usage and environment setup.

Verify your conda/mamba installation with:

```bash
# Verify conda works with
conda --version

# Verify mamba works with
mamba --version
```

### 1.1.1 Snakemake Conda/Mamba environment

Once your conda/mamba environment is working, you shall create a new conda/mamba environment and install snakemake with the following command.

> [!NOTE]
> For the rest of this tutorial, we'll use mamba commands. However, conda commands work the same way—simply substitute "conda" for "mamba."
<br>

```bash
mamba create -n snakemake -c bioconda snakemake>=9.11.0
```

Once the environment is created, test your snakemake installation.

```bash
mamba activate snakemake
snakemake --version
# This should report a snakemake version greater than 9.11.0
```

# 2. Workflow Download and Setup

Clone the repository from github onto a suitable directory on your local machine or cluster with;

```bash
git clone https://github.com/ChabbyTMD/svArcher_Standalone.git
```

Once the SaVor repository has been cloned onto your machine. You may set up your analysis in one of two ways depending on the location of your reads.

> [!NOTE]
> 1. SaVor functions with both local sequence reads or sequence reads available from the NCBI Sequence Read Archive (NCBI).
> 
> 2. SaVor only runs with paired-end sequence data. Please check all your sequence reads prior to running your workflow.
<br>

> [!TIP]
> Check the fastq report and pay special attention to the Sequence length distribution section. The length of your R1 and R2 files match.
<br>

Detail all sequence reads you intend to analyze in the `.csv` sample sheet located in the `config/` directory, including the appropriate metadata. Refer to the next section for a description of sample sheet columns.

## 2.1. Description of sample sheet columns

### 2.1.1 Mandatory Fields

1. `BioSample`. This is the unique sample name for your local paired-end read data. For NCBI SRA data this will be the BioSample identifier.
2. `Run`. This can either be the sample name, if it was sequenced on a single lane. For multi-lane samples this column should be the sequence lane number for the R1 and R2 of the particular sequence record. For NCBI SRA data this is the run identifier listed on the samples’ SRA run selector record.
3. `LibraryName`: For single and multi-lane local data, this can be the name of the sample prefixed with `lib_` . For multi-lane local fastq files, the library name should be the same for each sample. The following example represents a single sample, sample1, sequenced on two lanes, 1 and 2;

```bash
BioSample,Run,LibraryName,refGenome,fq1,fq2
sample1,1,lib_sample1,GCF_XXXXXXXXX.X,/path/to/sample1_lane1_1.fastq.gz,/path/to/sample1_lane1_2.fastq.gz
sample1,2,lib_sample1,GCF_XXXXXXXXX.X,/path/to/sample1_lane2_1.fastq.gz,/path/to/sample1_lane2_2.fastq.gz
```

1. `refGenome`: Provide the NCBI RefSeq accession of the reference genome for your analysis. All samples must use the same RefSeq identifier. For custom assemblies, use a consistent identifier.

### 2.1.2 Optional fields

1. `refPath`: For a custom assembly or reference genome e.g. masked reference, input the absolute or relative (to the repository root directory) path here. This column should be omitted if using an NCBI refSeq genome.
2. `fq1` and `fq2`: These are the R1 and R2 absolute or relative paths for your local samples’ paired end data. Note that for samples sequenced on multiple lanes, each lanes R1 and R2 must be present in the same record.

## 2.2. Workflow configuration

In the sample config.yaml you’ll find a number of parameters you’ll need to change specific to your analysis.

### Mandatory Options:

1. `samples` : The relative or absolute path to your samples.csv file.
2. `include_contigs` : A simple text file with the chromosome/contig names of your reference genome/assembly. Each entry should exist on its own line.
3. `use_custom_reference` : Set to `True` if you are using a locally available reference genome. Set to `False` if you want your reference downloaded from NCBI’s refSeq database.

### Optional options:

1. `sort_reads`: Set to `True` for the workflow to sort .bam files
2. `mark_duplicates`: Set to `True` for the workflow to mark duplicates
3. `svBenchmark`: Set to `True` to enable the workflows’ benchmark sub-module.
4. Path to SV benchmark base sets, `deletions`, `duplications`, `inversions` : These paths must be defined for your ground truth SV calls if the benchmark module is activated.

# 3. Running SaVor

Use the following as general guide to determine the columns you need to include in your sample.csv metadata file to set up your structural variant calling run.

> [!NOTE]
> For a quick tutorial with Arabidopsis data from NCBI follow this [link](savor_tutorial.md).

## 3.1. Choosing your sample sheet layout

### 3.1.1. Single lane samples from NCBI SRA using a reference genome from NCBI refSeq

```csv
Run,BioSample,LibraryName,refGenome
SRRXXXXXXX,SAMNXXXXXXXX,SAMNXXXXXXXX,GCF_XXXXXXXXX.X
```

`use_custom_reference`: False

### 3.1.2. Single lane local reads and NCBI refSeq reference genome

```csv
Run,BioSample,LibraryName,refGenome,fq1,fq1
sample1,sample1,lib_sample1,GCF_XXXXXXXXX.X,/path/to/sample1_lane1_1.fastq.gz,/path/to/sample1_lane1_2.fastq.gz
```

`use_custom_reference`: False

### 3.1.3. Single lane local reads and a custom reference

```csv
Run,BioSample,LibraryName,refGenome,refPath,fq1,fq1
sample1,sample1,lib_sample1,GCF_XXXXXXXXX.X,/path/to/custom/reference/REF_NAME.fna,/path/to/sample1_lane1_1.fastq.gz,/path/to/sample1_lane1_2.fastq.gz
```

`use_custom_reference`: True

### 3.1.4. Mutli lane local reads and NCBI refSeq reference genome

```csv
Run,BioSample,LibraryName,refGenome,fq1,fq1
sample1,1,lib_sample1,GCF_XXXXXXXXX.X,/path/to/sample1_lane1_1.fastq.gz,/path/to/sample1_lane1_2.fastq.gz
sample1,2,lib_sample1,GCF_XXXXXXXXX.X,/path/to/sample1_lane2_1.fastq.gz,/path/to/sample1_lane2_2.fastq.gz
sample2,1,lib_sample2,GCF_XXXXXXXXX.X,/path/to/sample2_lane1_1.fastq.gz,/path/to/sample2_lane1_2.fastq.gz
sample2,2,lib_sample2,GCF_XXXXXXXXX.X,/path/to/sample2_lane2_1.fastq.gz,/path/to/sample2_lane2_2.fastq.gz
```

`use_custom_reference`: False

### 3.1.5. Mutli lane local reads and a custom reference

```csv
Run,BioSample,LibraryName,refGenome,refPath,fq1,fq1
sample1,1,lib_sample1,GCF_XXXXXXXXX.X,/path/to/custom/reference/REF_NAME.fna,/path/to/sample1_lane1_1.fastq.gz,/path/to/sample1_lane1_2.fastq.gz
sample1,2,lib_sample1,GCF_XXXXXXXXX.X,/path/to/custom/reference/REF_NAME.fna,/path/to/sample1_lane2_1.fastq.gz,/path/to/sample1_lane2_2.fastq.gz
sample2,1,lib_sample2,GCF_XXXXXXXXX.X,/path/to/custom/reference/REF_NAME.fna,/path/to/sample2_lane1_1.fastq.gz,/path/to/sample2_lane1_2.fastq.gz
sample2,2,lib_sample2,GCF_XXXXXXXXX.X,/path/to/custom/reference/REF_NAME.fna,/path/to/sample2_lane2_1.fastq.gz,/path/to/sample2_lane2_2.fastq.gz
```

`use_custom_reference`: True

### 3.1.6 Single lane with custom reference and user provided BAM files

```csv
Run,BioSample,LibraryName,refGenome,refPath,bamPath,baiPath
sample1,1,lib_sample1,GCF_XXXXXXXXX.X,REF_NAME.fna,/path/to/your/sample1.bam,/path/to/your/sample1.bam.bai
sample2,1,lib_sample2,GCF_XXXXXXXXX.X,REF_NAME.fna,/path/to/your/sample2.bam,/path/to/your/sample2.bam.bai
```

For further details on executing the workflow from BAM files refer to the guide [here](user_provided_bams.md).


## 3.2. Executing SaVor

After setting up your sample sheet and workflow configuration in the `config.yaml` file, its good to perform a dry run to confirm whether you’ve set everything up correctly.

### 3.2.1. Local Execution

While in the workflow root directory with the snakmake mamba environment activated, perform a dry run with the following command:

```bash
snakemake -np --cores 1 --workflow-profile workflow-profiles/default
```

The options specify;

`-np` : Perform a dry run and print out shell commands to stdout for each rule file

`--cores`: Specify how many cores to provide to snakemake. Specifying `--all` will use all available cores on your system.

`--workflow-profile`: A path to a configuration file that species the number of threads to dedicate to certain rules. Scale up and down these threads based on the capacity of your system. Current defaults assume a system with at least 16 cores/threads.

Assuming the workflow dry run completed with no errors, execute the pipeline with:

```bash
snakemake -p --cores all --workflow-profile workflow-profiles/default
```

### 3.2.2 Cluster Execution

To execute on a SLURM cluster, you need to first ensure the right executor plugin are installed in the snakemake environment. At the moment, SaVor only supports SLURM execution with the `cluster-generic` [snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html). Depending on the resources present on your cluster, scale the number of jobs based on how many nodes you can dedicate to your job. For example, if you have access to 2 nodes with 16 cores each, the `jobs` parameter ought to be set to 32. This will effectively utilize two entire nodes for your snakemake jobs.

Execute the workflow dry run from the head-node by running the following command:

```bash
snakemake -np --workflow-profile workflow-profiles/cluster-generic-slurm
```

Once the dry run completes with no errors, proceed to execute the workflow with:

```bash
snakemake -p --workflow-profile workflow-profiles/cluster-generic-slurm
```

# 4. Troubleshooting

## 4.1. Restarting SaVor

If the workflow is killed or interrupted, append the `--rerun-incomplete` flag to the snakemake command. This allows regeneration of output files for rules that were in progress during the interruption.

If the user had to manually interupt the workflow, the following command must first be executed before restarting the SaVor:

```bash
snakemake -np --cores all --unlock
```
