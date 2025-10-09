# User-Provided BAM Files

The workflow now supports using your own BAM files instead of generating them from FASTQ files. To use this feature:

1. Add two new columns to your sample sheet (`samples.csv`):
   - `bamPath`: Full path to your BAM file
   - `baiPath`: Full path to your BAI file (index)

Example sample sheet with user-provided BAM files:

```csv
Run,BioSample,LibraryName,refGenome,refPath,bamPath,baiPath
SRR1945442,SAMN03326286,SAMN03326286,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,/path/to/your/sample1.bam,/path/to/your/sample1.bam.bai
SRR1945443,SAMN03326287,SAMN03326287,GCF_000001735.3,REFERENCE/TAIR10_REF.fna,/path/to/your/sample2.bam,/path/to/your/sample2.bam.bai
```

When user-provided BAM files are detected, the workflow will:

1. Skip the FASTQ download/ingestion steps
2. Skip the read alignment steps
3. Copy your BAM files to the workflow's expected location (`results/{refGenome}/bams/{sample}_final.bam`)
4. Use these copied BAM files for SV calling

**Note**:

- Ensure your BAM files are properly sorted and indexed.
- The workflow will copy your BAM files to maintain a consistent structure for downstream steps.
- You can mix samples with user-provided BAMs and samples that need alignment in the same run.
- The `mark_duplicates` setting is ignored for user-provided BAMs, as they are used as-is. Ensure BAM files provided to the workflow are properly mark duplicated. One can do this with [sambamba](https://github.com/biod/sambamba) for example;
 ```bash
  sambamba markdup -t 8 -p input.bam output.bam
  ```
- The workflow performs file existence and permission checks at the time of use, providing clear error messages if files are missing or unreadable.
- Logs for user-provided BAM file operations are saved in `logs/{refGenome}/user_bams/`.
