
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

In the workflow `config.yaml` located in the `config/` directory you’ll find a number of parameters you’ll need to change specific to your analysis.

### Mandatory Options:

1. `samples` : The relative or absolute path to your `samples.csv` file.
2. `include_contigs` : A simple text file with the chromosome/contig names of your reference genome/assembly. Each entry should exist on its own line.
3. `use_custom_reference` : Set to `True` if you are using a locally available reference genome. Set to `False` if you want your reference downloaded from NCBI’s refSeq database.

### Optional options:

1. `sort_reads`: Set to `True` for the workflow to sort .bam files
2. `mark_duplicates`: Set to `True` for the workflow to mark duplicates
3. `svBenchmark`: Set to `True` to enable the workflows’ benchmark sub-module.
4. Path to SV benchmark base sets, `deletions`, `duplications`, `inversions` : These paths must be defined for your ground truth SV calls if the benchmark module is activated.