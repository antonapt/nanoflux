#!/bin/bash

set -e

ENV_NAME="nanoflux"
REF_DIR="./refs"

echo "Starting setup for $ENV_NAME..."

# 1. Conda Check
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed. Please install Conda before running this script or check the README for alternative installation methods."
    exit 1
fi

echo " -------- Environment Creation ($ENV_NAME with Python 3.12) --------"
conda create -n $ENV_NAME python=3.12 -y
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo " -------- Installing Bioinformatics Tools (Bioconda) --------"
conda install -y -c conda-forge -c bioconda uv
if [[ "$OSTYPE" == darwin* ]]; then
    conda install -y -c conda-forge -c bioconda samtools==1.21 bedtools==2.31.1 ont-modkit minimap2
elif [[ "$OSTYPE" == linux-gnu* ]]; then
    conda install -y -c conda-forge -c bioconda samtools==1.21 bedtools==2.31.1 ont-modkit==0.2.5 minimap2==2.26
else
    echo "Error: Unsupported OS type ($OSTYPE). This script only supports macOS and Linux."
    exit 1
fi

echo " -------- Install Python dependencies --------"
if [ -f "pyproject.toml" ]; then
    echo "Installing package and dependencies into conda environment..."
    if ! uv pip install -e .; then
        echo "Error: Failed to install dependencies with uv!"
        echo "Reverting environment setup..."
        conda deactivate
        conda remove -n $ENV_NAME --all -y
        echo "Environment $ENV_NAME has been removed."
        exit 1
    fi
else
    echo "Warning: pyproject.toml not found."
fi


echo " -------- Download Reference Genome --------"
mkdir -p "$REF_DIR"
if [ ! -f "$REF_DIR/chm13v2.fa" ]; then
    echo "Downloading CHM13 reference genome..."
    curl -L https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o "$REF_DIR/chm13v2.fa.gz"
    gunzip -c "$REF_DIR/chm13v2.fa.gz" > "$REF_DIR/chm13v2.fa"
    rm "$REF_DIR/chm13v2.fa.gz"
fi

echo " -------- Model Download --------"
if [ -f "src/data/download_models.py" ]; then
    python -m src.data.download_models
else
    echo "Warning: download_models.py not found at src/data/."
fi

echo "------------------------------------------------"
echo ""
echo "Setup complete..."
echo ""
echo "Run: conda activate $ENV_NAME"
echo ""
echo "------------------------------------------------"
