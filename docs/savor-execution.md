
# 4. Executing SaVor

After setting up your sample sheet and workflow configuration in the `config.yaml` file, it's a good idea to perform a dry run to confirm whether you've set everything up correctly.

### 4.1. Local Execution

While in the workflow root directory with the snakemake mamba environment activated, perform a dry run with the following command:

```bash
snakemake -np --cores 1 --workflow-profile workflow-profiles/default
```

The options above specify the following;

`-np` : Perform a dry run and print out shell commands to stdout for each rule file

`--cores`: Specify how many cores to provide to snakemake. Specifying `--all` will use all available cores on your system.

`--workflow-profile`: A path to a configuration file that species the number of threads to dedicate to certain rules. Scale up and down these threads based on the capacity of your system. Current defaults assume a system with at least 16 cores/threads.

Assuming the workflow dry run completed with no errors, execute the pipeline with:

```bash
snakemake -p --cores all --workflow-profile workflow-profiles/default
```

### 4.2 Cluster Execution

To execute on a SLURM cluster, you need to first ensure the right executor plugin are installed in the snakemake environment. At the moment, SaVor only supports SLURM execution with the `cluster-generic` [snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html). Depending on the resources present on your cluster, scale the number of jobs based on how many nodes you can dedicate to your job. For example, if you have access to 2 nodes with 16 cores each, the `jobs` parameter ought to be set to 32. This will effectively utilize two entire nodes for your snakemake jobs.

Execute the workflow dry run from the head-node by running the following command:

```bash
snakemake -np --workflow-profile workflow-profiles/cluster-generic-slurm
```

Once the dry run completes with no errors, proceed to execute the workflow with:

```bash
snakemake -p --workflow-profile workflow-profiles/cluster-generic-slurm
```