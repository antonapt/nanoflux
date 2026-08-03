#!/bin/bash

set -e

ENV_NAME="nanoflux"
REF_DIR="./data/refs"
REFERENCE="chm13v2"

usage() {
    echo "Usage: $0 [--hg38 | --ch13v2 | --chm13v2]"
    echo "  --hg38      Use hg38 reference genome (default)"
    echo "  --chm13v2   Use CHM13 v2 reference genome"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hg38)
            REFERENCE="hg38"
            shift
            ;;
        --hg38-as)
            REFERENCE="hg38-as"
            shift
            ;;
        --chm13v2)
            REFERENCE="chm13v2"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option '$1'"
            usage
            ;;
    esac
done

echo "Starting setup for $ENV_NAME..."
echo "Selected reference genome: $REFERENCE"

# 1. Conda Check
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed. Please install Conda before running this script or check the README for alternative installation methods."
    exit 1
fi

echo " -------- Environment Creation ($ENV_NAME with Python 3.12) --------"
conda create -n $ENV_NAME python=3.12 -y
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: Conda environment is not active."
    exit 1
fi

CONDA_PYTHON="$CONDA_PREFIX/bin/python"

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

echo " -------- Download Reference Genome --------"
mkdir -p "$REF_DIR"

if [[ "$REFERENCE" == "hg38" ]]; then
    if [ ! -f "$REF_DIR/hg38.fa" ]; then
        echo "Downloading hg38 reference genome..."
        curl -L --progress-bar https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -o "$REF_DIR/hg38.fa.gz"
        gunzip -c "$REF_DIR/hg38.fa.gz" > "$REF_DIR/hg38.fa"
        rm "$REF_DIR/hg38.fa.gz"
    fi
elif [[ "$REFERENCE" == "hg38-as" ]]; then
    if [ ! -f "$REF_DIR/hg38-as.fna" ]; then
        echo "Downloading hg38 analysis set reference genome..."
        curl -L --progress-bar https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
        -o "$REF_DIR/hg38-as.fna.gz"
        gunzip -c "$REF_DIR/hg38-as.fna.gz" > "$REF_DIR/hg38-as.fa"
        rm "$REF_DIR/hg38-as.fna.gz"
    fi
else
    if [ ! -f "$REF_DIR/chm13v2.fa" ]; then
        echo "Downloading CHM13 reference genome..."
        curl -L --progress-bar https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o "$REF_DIR/chm13v2.fa.gz"
        gunzip -c "$REF_DIR/chm13v2.fa.gz" > "$REF_DIR/chm13v2.fa"
        rm "$REF_DIR/chm13v2.fa.gz"
    fi
fi


echo " -------- Model Download --------"
if [ -f "data/download_models.py" ]; then
    uv run --python "$CONDA_PYTHON" --no-project --with zenodo_get data/download_models.py
else
    echo "Warning: download_models.py not found at data/."
fi

echo " -------- Install Python dependencies --------"
if [ -f "pyproject.toml" ]; then
    echo "Installing package and dependencies into conda environment..."
    if ! uv pip install --python "$CONDA_PYTHON" -e .; then
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

echo "------------------------------------------------"
echo ""
echo "Setup complete..."
echo ""
echo "Run: conda activate $ENV_NAME"
echo ""
echo "------------------------------------------------"
