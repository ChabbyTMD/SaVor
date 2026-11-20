# Savor Tutorial

In this quick tutorial, we shall run SaVor to generate a union set of structural variants from a small set of Arabidopsis whole genome samples from NCBI.

For brevity, we shall use BAM files, their associated indices of 4 *Arabidopsis thaliana* and the TAIR10 reference sequence contained in the SaVor tutorial [Google Drive](https://drive.google.com/drive/folders/1ZZhS3nRH0XJbpVIZpgbx_NQjUkCvnI4U?usp=sharing) directory.

This tutorial assumes you have;

1. Set up and configured your snakemake mamba/conda environment and verified snakemake is installed and functional. For instructions, please refer to the prerequisites [here](prerequisites.md)

2. You have cloned the repo onto your local machine.

## Step 0: Download Reference Sequence and BAM Files.

Change directory into the `SaVor` repo root directory. Create a new directory and name it `REFERENCE` and another directory named `bams/`. In order to download the tutorial data, run the following commands in your terminal.

```bash
# Install gdown, a command line tool to download files from Google Drive.
pip install gdown

# Download the tutorial data
gdown https://drive.google.com/drive/folders/1ZZhS3nRH0XJbpVIZpgbx_NQjUkCvnI4U

# Copy over tutorial data into the bams and REFERENCE directories you created.
mv SaVor_Tutorial_Files/bams/* bams/
mv SaVor_Tutorial_Files/REFERENCE/* REFERENCE/

#Clean up download folder
rm -rf SaVor_Tutorial_Files/
```


## Step 1: Setup the workflow `config.yaml`


A tutorial configuration file has been provided at `config/tutorial-config.yaml`. Below are the explanations of the attributes.

1. `samples:` - A relative or absolute path to a workflow sample sheet CSV file. A sample sheet has been provided at `config/bam_samples.csv`. Below is a preview of the tutorial sample sheet.


>[!NOTE]
> ```csv
>Run,BioSample,LibraryName,refGenome,refPath,bamPath,baiPath
>SRR1945442,SAMN03326286,SAMN03326286,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326286_final.bam,bams/SAMN03326286_final.bam.bai
>SRR1945443,SAMN03326287,SAMN03326287,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326287_final.bam,bams/SAMN03326287_final.bam.bai
>SRR1945444,SAMN03326288,SAMN03326288,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326288_final.bam,bams/SAMN03326288_final.bam.bai
>SRR1945445,SAMN03326289,SAMN03326289,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,bams/SAMN03326289_final.bam,bams/SAMN03326289_final.bam.bai
>```


2. `include_contigs:` - A path pointing to the `include_contigs.csv` file containing the contig/chromosome IDs of the TAIR10 Reference.

Alternatively, you can easily determine the contig IDs from the reference file itself using the `grep` command.

```bash
grep "^>" REFERENCE/TAIR10_REF.fna > include_contigs.csv
```

3. `use_custom_reference:` = Set to `True` to enable the tutorial reference genome.

4. `user_provided_bams:` = Set to `True` to enable the use of user-provided BAM files.

5. `mark_duplicates:` = Set to `False` since we are using user-provided BAM files that are already marked for duplicates.

6. `sv_merge` - Change the value to `1` to generate a union set of SVs of the 4 samples.

>[!NOTE]
> A tutorial configuration file has been provided in the `config` directory as `tutorial-config.yaml`.


## Step 3: Perform a Dry Run

Next, activate the snakemake conda/mamba environment unless its already active. Run the following Snakemake command to generate a list of jobs SaVor will execute. This should execute in a couple of seconds. If you encounter any errors, ensure you are in the repository root directory and verify whether you have all necessary configuration and data files present.

```bash
snakemake --cores 1 -np --workflow-profile workflow-profiles/default/ --configfile config/tutorial-config.yaml
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
snakemake --cores 8 -np --workflow-profile workflow-profiles/default/ --configfile config/tutorial-config.yaml

```

>[!NOTE]
> Note this tutorial will not run the benchmark module.
