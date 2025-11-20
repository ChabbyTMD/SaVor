# 3. Sample Sheet Setup

Use the following as general guide to determine the columns you need to include in your sample.csv metadata file to set up your structural variant calling run.


## 3.1. Choosing your sample sheet layout

Depending on what scenario best describes your data implement the recommended columns. Below each example is a corresponding boolean value for what the `use_custom_reference` directive should be in the `config.yaml` file.

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
`use_custom_reference`: True

For further details on executing the workflow from BAM files refer to the guide [here](user_provided_bams.md).

