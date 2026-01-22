# 🧬 NanoFlux

DNA methylation ensemble classifier for Oxford Nanopore sequencing data.

## 📋 Requirements

- **Conda** installed on your computer: [Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

## 🔧 Installation

Run the installation script to set up the environment and download necessary dependencies:

```bash
bash install_package.sh
```

This script will:
- Create a Conda environment with Python 3.12
- Install bioinformatics tools (samtools, bedtools, ont-modkit, minimap2)
- Install Python dependencies
- Download the CHM13 reference genome
- Download pre-trained ensemble models

After installation, activate the environment:

```bash
conda activate nanoflux_env
```

## 🚀 Usage

### 1. 📁 Prepare BAM Files

[TODO]

### 2. 🤖 Infer with NanoFlux

Run the NanoFlux inference pipeline on your `.bed` file (output from previous step):

```bash
uv run python -m src.main infer -i "$input_file" -o "$output_dir" -c
```

The inference will:
- Extract methylation features from the prepared `.bed` file (`input_file`)
- Apply the ensemble classifier
- Output predictions and probability scores to the specified `output_dir`
- The `-c` flag indicates that output directories will be created if they do not exist

## 🧪 Development

Run tests to verify dependencies:

```bash
uv run pytest
```