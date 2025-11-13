from snakemake.exceptions import WorkflowError

rule link_user_bam:
    """Link or copy user-provided BAM files to the expected location in the workflow.
    This rule will only be triggered when user-provided BAM paths are specified in the sample sheet.
    The rule includes robust error handling for file not found and access errors."""
    output:
        bam = "results/{refGenome}/bams/{sample}_final.bam",
        bai = "results/{refGenome}/bams/{sample}_final.bam.bai"
    params:
        user_bams = lambda wc: get_user_bams(wc.sample)
    log:
        "logs/{refGenome}/user_bams/{sample}.log"
    run:
        import os
        import shutil
        import sys
        from pathlib import Path
        
        # Create output directories for both output files and logs
        Path(output.bam).parent.mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        
        with open(log[0], "w") as log_file:
            try:
                bam_path = params.user_bams["bam"]
                bai_path = params.user_bams["bai"]
                
                # Verify that the BAM file exists and is readable
                if not os.path.exists(bam_path):
                    raise FileNotFoundError(f"User-provided BAM file not found: {bam_path}")
                if not os.access(bam_path, os.R_OK):
                    raise PermissionError(f"Cannot read user-provided BAM file (check permissions): {bam_path}")
                
                # Verify that the BAI file exists and is readable
                if not os.path.exists(bai_path):
                    raise FileNotFoundError(f"User-provided BAI file not found: {bai_path}")
                if not os.access(bai_path, os.R_OK):
                    raise PermissionError(f"Cannot read user-provided BAI file (check permissions): {bai_path}")
                
                # Remove existing files if they exist
                if os.path.lexists(output.bam):
                    os.remove(output.bam)
                if os.path.lexists(output.bai):
                    os.remove(output.bai)
                
                # Copy the files
                shutil.copy2(bam_path, output.bam)
                shutil.copy2(bai_path, output.bai)
                
                message = f"Successfully copied user-provided BAM file from {bam_path} to {output.bam}"
                print(message)
                log_file.write(message + "\n")
                
            except (FileNotFoundError, PermissionError, OSError) as e:
                error_message = f"ERROR: {str(e)}"
                log_file.write(error_message + "\n")
                print(error_message, file=sys.stderr)
                raise WorkflowError(error_message)