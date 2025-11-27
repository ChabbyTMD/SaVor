# Savor Tutorial

In this quick tutorial, we shall run SaVor to generate a union set of structural variants from a small set of Arabidopsis whole genome samples from NCBI.

For brevity, we shall use BAM files, their associated indices of 4 *Arabidopsis thaliana* and the TAIR10 reference sequence contained in the SaVor tutorial [Google Drive](https://drive.google.com/drive/folders/1ZZhS3nRH0XJbpVIZpgbx_NQjUkCvnI4U?usp=sharing) directory.

This tutorial assumes you have;

1. Set up and configured your snakemake mamba/conda environment and verified snakemake is installed and functional. Instructions for which are in the user guide [here](README.md)

2. You have cloned the repo onto your local machine.

## Step 0: Download Reference Sequence and BAM Files.

Create a new directory and name it `REFERENCE`. Download the TAIR10 reference genome from this [link](https://drive.google.com/file/d/1X9S57uUYM7A1Vxc8C2XKeR4Hk-tkFU6l/view?usp=drive_link) and move it into the `REFERENCE` directory. Additionally, create an empty directory named `bams/` and download all files from this [Google Drive](https://drive.google.com/drive/folders/1v5DyYavAvpMPYInNnZxpix4m9zvk7jkr?usp=sharing) folder into it.


## Step 1: Setup the workflow `config.yaml`


Change into the `config` directory and edit the following options in the `config.yaml` file with a text editor of your choice.

1. `samples:` - Edit this path to point to the `test.csv` sample sheet. 


>[!NOTE]
> Create a new file in the `config` directory with the following content and name it `test.csv`
> ```csv
>Run,BioSample,LibraryName,refGenome,refPath,bamPath,baiPath
>SRR1945442,SAMN03326286,SAMN03326286,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326286_final.bam,bams/SAMN03326286_final.bam.bai
>SRR1945443,SAMN03326287,SAMN03326287,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326287_final.bam,bams/SAMN03326287_final.bam.bai
>SRR1945444,SAMN03326288,SAMN03326288,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326288_final.bam,bams/SAMN03326288_final.bam.bai
>SRR1945445,SAMN03326289,SAMN03326289,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326289_final.bam,bams/SAMN03326289_final.bam.bai
>```


2. `include_contigs:` - Edit this path to point to the file `include_contigs.csv` containing the contig/chromosome IDs of the TAIR10 Reference.

While in the `config` directory, create an empty file with all the contig IDs of your reference genome. Below are the contig IDs for the TAIR10 reference you downloaded earlier.

```
Chr1
Chr2
Chr3
Chr4
Chr5
```
Alternatively, you can easily determine the contig IDs from the reference file itself using the `grep` command.

```bash
grep "^>" ../REFERENCE/TAIR10_REF.fna > include_contigs.csv
```

3. `sv_merge` - Change the value to `1` to generate a union set of SVs of the 4 samples.


## Step 3: Perform a Dry Run

Change out of the `config` directory and into the SaVor root directory. Next, activate the snakemake conda/mamba environment unless its already active. Run the following Snakemake command to generate a list of jobs SaVor will execute. This should execute in a couple of seconds. If you encounter any errors, ensure you are in the repository root directory and verify whether you have all necessary configuration and data files present. Refer to the [user guide](README.md) if needed.

```bash
snakemake --cores 1 -np --workflow-profile workflow-profiles/default/
```

If your workflow is configured correctly, you should see the following message print out.


```{image} ../img/savor_dry_run.png
:width: 1200
:height: 400
:align: center
```

## Step 4: Perform a Wet Run

Once your dry run completes successfully, you may remove the `n` option above to run the workflow.

```bash
# This assumes you have at least 8 cores on your local machines CPU.
snakemake --cores 8 -np --workflow-profile workflow-profiles/default/

```

>[!NOTE]
> Note this tutorial will not run the benchmark module.
