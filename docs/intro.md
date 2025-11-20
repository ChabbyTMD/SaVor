# SaVor: A Structural Variant Calling and Benchmarking Workflow for Short-Read Sequence Data.

```{image} ../img/savor_logo.png
:width: 350px
:align: center
```

Structural variant calling from short-read sequence data is a complex task often requiring multiple tools and various file types. With SaVor, we have streamlined the process of predicting structural variants from whole-genome resequencing projects by designing a flexible workflow that can utilize NCBI SRA sequence data or your local Illumina sequence reads from one or multiple lanes. In case you have on hand BAM files that you may have generated as part of a SNP analysis project, those can be used by SaVor to save on computational time.

SaVor can additionally benchmark SV callset against ground truth sets that you provide. For now, benchmarking is limited to deletion, inversion and duplication SV types. This  feature is disabled by default but can be enabled in the workflow `config.yaml`. You can find more information on how to enable it [here](workflow-setup.md#optional-options).

## Workflow Requirements

1. A working conda/mamba installation.
2. Reference genome (a local .fasta/.fna file or an NCBI RefSeq accession)
3. Illumina paired-end sequence reads.
4. Optionally, ground truth deletion, duplication and inversion SV sets.

## Resource Requirements

- **Memory**: 8-16GB per sample for alignment, 4-8GB for SV calling
- **Storage**: ~10-50GB per sample depending on coverage
- **Time**: 2-8 hours per sample depending on coverage and resources

## Workflow Tutorial

For a quick example on how to setup and execute the workflow, please follow the guide [here](savor_tutorial.md).

## User Guide

For a more detailed user guide, you can jump straight in [here](prerequisites.md)

## Citation


If you use SaVor, please cite in addition to all the tools implemented:

```text
SaVor - A Reproducible Structural Variant Calling and Benchmarking Platform from Short-Read Data
[In Prep]
```

## Support

For questions and issues, please submit an issue on the [repo](https://github.com/ChabbyTMD/svArcher_Standalone/issues/new).