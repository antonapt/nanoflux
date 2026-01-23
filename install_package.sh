#!/bin/bash

set -e

ENV_NAME="nanoflux_env"
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
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -y uv
if [[ "$OSTYPE" == darwin* ]]; then
    conda install -y samtools==1.21 bedtools==2.31.1 ont-modkit minimap2
elif [[ "$OSTYPE" == linux-gnu* ]]; then
    conda install -y samtools==1.21 bedtools==2.31.1 ont-modkit==0.2.5 minimap2==2.26
else
    echo "Error: Unsupported OS type ($OSTYPE). This script only supports macOS and Linux."
    exit 1
fi

echo " -------- Install Python dependencies with Fallback for OpenMP --------"
if [ -f "pyproject.toml" ]; then
    echo "Attempting to sync dependencies..."

    export UV_PYTHON=$(which python)

    if ! uv sync; then
    echo "------------------------------------------------"
    echo "Sync failed. This is likely due to liblinear-multicore on your system."
    echo "Switching to standard liblinear for compatibility..."
    echo "------------------------------------------------"

python << 'EOF'
import pathlib
import shutil

original = pathlib.Path('pyproject.toml')
backup = pathlib.Path('pyproject.toml.bak')
orig_lock = pathlib.Path('uv.lock')
backup_lock = pathlib.Path('uv.lock.bak')

if original.exists():
    shutil.copy2(original, backup)

if orig_lock.exists():
    shutil.copy2(orig_lock, backup_lock)

lines = original.read_text().splitlines(keepends=True)

with open(original, 'w') as f:
    for line in lines:
        if 'liblinear-multicore' in line:
            f.write('    "liblinear-official",\n')
        else:
            f.write(line)
EOF
    cat pyproject.toml
    echo "Retrying sync with liblinear-official..."
    uv sync
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
    uv run python -m src.data.download_models
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
