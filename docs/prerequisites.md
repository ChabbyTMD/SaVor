
# 1. Prerequisites

## 1.1 Conda/Mamba

Ensure you have conda/mamba installed on the system you intend to run the workflow on. A detailed mamba installation guide can be found [here](https://www.howtogeek.com/how-to-set-up-a-development-environment-with-mamba/). For cluster environments, please contact your system administrator for conda/mamba usage and environment setup.

Verify your conda/mamba installation with:

```{code} bash
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