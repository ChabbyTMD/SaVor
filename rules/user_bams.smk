rule link_user_bam:
    """Link or copy user-provided BAM files to the expected location in the workflow"""
    output:
        bam = "results/{refGenome}/bams/{sample}_final.bam",
        bai = "results/{refGenome}/bams/{sample}_final.bam.bai"
    params:
        user_bams = lambda wc: get_user_bams(wc.sample)
    run:
        import os
        import shutil
        from pathlib import Path
        
        # Create output directory
        Path(output.bam).parent.mkdir(parents=True, exist_ok=True)
        
        # Link or copy the BAM file
        if os.path.lexists(output.bam):
            os.remove(output.bam)
        shutil.copy2(params.user_bams["bam"], output.bam)
        
        # Link or copy the BAI file
        if os.path.lexists(output.bai):
            os.remove(output.bai)
        shutil.copy2(params.user_bams["bai"], output.bai)

# Define a conditional workflow based on whether user BAMs are provided
def get_alignment_input(wildcards):
    """Determine if we need to run the alignment workflow or use user-provided BAMs"""
    if has_user_bams(wildcards.sample):
        # If user-provided BAMs exist, point to the output of link_user_bam
        return {
            "bam": f"results/{wildcards.refGenome}/bams/{wildcards.sample}_final.bam",
            "bai": f"results/{wildcards.refGenome}/bams/{wildcards.sample}_final.bam.bai"
        }
    else:
        # Otherwise, use the normal alignment workflow inputs
        if config.get("mark_duplicates", True):
            return {
                "bam": f"results/{wildcards.refGenome}/bams/{wildcards.sample}_final.bam",
                "bai": f"results/{wildcards.refGenome}/bams/{wildcards.sample}_final.bam.bai"
            }
        else:
            return dedup_input(wildcards)