#!/bin/bash
# Test script for svArcher Standalone Pipeline

echo "=== svArcher Standalone Test Script ==="
echo ""

# Check if we're in the right directory
if [ ! -f "Snakefile" ]; then
    echo "Error: Snakefile not found. Please run this script from the svArcher_Standalone directory."
    exit 1
fi

# Check if conda/mamba is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "Error: Neither conda nor mamba found. Please install one of them."
    exit 1
fi

echo "Using $CONDA_CMD for environment management"

# Activate snakemake environment
echo "Activating snakemake environment..."
eval "$(conda shell.bash hook)"
conda activate snakemake

# Check if snakemake is available in the environment
if ! command -v snakemake &> /dev/null; then
    echo "Error: Snakemake not found in the 'snakemake' environment."
    echo "Please ensure snakemake is installed in the environment:"
    echo "$CONDA_CMD install -c conda-forge -c bioconda -n snakemake snakemake>=7.0.0"
    exit 1
fi

echo "Snakemake version: $(snakemake --version)"

# Validate configuration files
echo ""
echo "=== Validating Configuration ==="
if [ ! -f "config/config.yaml" ]; then
    echo "Error: config/config.yaml not found"
    exit 1
fi

if [ ! -f "config/samples.csv" ]; then
    echo "Error: config/samples.csv not found"
    exit 1
fi

if [ ! -f "config/include_contigs.csv" ]; then
    echo "Error: config/include_contigs.csv not found"
    exit 1
fi

echo "✓ Configuration files found"

# Dry run test
echo ""
echo "=== Running Snakemake Dry Run ==="
snakemake -n --use-conda

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Dry run successful! The workflow is ready to execute."
    echo ""
    echo "To run the full pipeline:"
    echo "  snakemake --cores 8 --use-conda"
    echo ""
    echo "For cluster execution:"
    echo "  snakemake --cores 100 --use-conda --cluster 'your-cluster-command'"
else
    echo ""
    echo "✗ Dry run failed. Please check the error messages above."
    exit 1
fi