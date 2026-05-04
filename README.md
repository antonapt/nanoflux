# 🧬 NanoFlux

DNA methylation ensemble classifier for Oxford Nanopore sequencing data. This python package accompanies the work:
*"Advancing methylation-based brain tumor classification through Nanopore sequencing data of cerebrospinal fluid-derived cell-free DNA"*  by Julia Freese, Anton Appelt, Laure Ciernik et al.

## 📋 Requirements

- **Conda** installed on your computer: [Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- A **basecaller** of your choice. We recommend ONT Dorado, see "Usage" for model recommendation. [Installation Guide](https://software-docs.nanoporetech.com/dorado/latest/)

## 🔧 Installation

Clone the repository and run the installation script to set up the environment and download necessary dependencies:

```console
git clone https://github.com/antonapt/nanoflux.git
cd nanoflux
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
conda activate nanoflux
```

## 🚀 Usage

### 1. 🔎 Base calling

We recommend you perform base calling with `dorado` and it's `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` model or newer:

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 "$pod5_path" --modified-bases 5mCG_5hmCG > "$out_BAM"
```

### 2. 📁 Prepare BAM Files

Nanoflux will handle both SAM and BAM files obtained from base calling. As these are often already aligned to some reference sequence, Nanoflux will automatically re-align your files to the T2T CHM13v2.0 reference. If your data happens to be already aligned to this reference, you can pass the `--skip-alignment` flag to the `prepare` command in this step.

Run the NanoFlux file preparation pipeline on a SAM/BAM file:

```bash
nanoflux prepare -i "$input_file" -o "$output_directory" 
```

Additional flags and arguments:
```bash
-c/--create-dir: "create "$output_directory" if it does not exist"
-t/--threads "$n_threads": "number of threads to use (default: 4)"
--skip-alignment: "skip minimap2 alignment if SAM/BAM is already aligned to T2T CHM13v2.0"
```

This will:
- align SAM/BAM with `minimap2` (optional)
- sort, index and convert to BAM with `samtools`
- extract and and tabulate methylation information with `modkit`
- annotate methylation calls with `bedtools`

### 3. 🤖 Infer with NanoFlux

Run the NanoFlux inference pipeline on your `methylation.bed` file (output from previous step):

```bash
nanoflux infer -i "$input_file" -o "$output_directory" -c
```

By default, the brain tumor classifier is used. Use the `-m`/`--model` flag to select a different model:

```bash
nanoflux infer -i "$input_file" -o "$output_directory" -c -m pancancer_coarse_lbls
```

#### Available models

| Model | Classes | Description |
|---|---|---|
| `brain_only` *(default)* | 21 | Brain tumor methylation classifier (ATRT, DMG, EPN subtypes, GBM, IDH GLM, MB subtypes, PA, PLEX, and controls). Corresponds to the model described in the accompanying publication. |
| `pancancer_coarse_lbls` | 33 | Brain tumor classifier extended with 12 broad extra-cranial cancer categories (bladder, breast, colon, kidney, liver, lung, lymphoma/leukemia, melanoma, ovary, pancreas, prostate, sarcoma). Use when the tumor origin is unknown. |
| `pancancer_fine_lbls` | 107 | Fine-grained pan-cancer classifier covering brain tumor subtypes and detailed histological subtypes across all major cancer types. |

The inference will:
- Extract methylation features from the prepared `.bed` file (`input_file`)
- Apply the ensemble classifier
- Output predictions and probability scores to the specified `output_directory`
- The `-c` flag indicates that output directories will be created if they do not exist

## 🧪 Development

Run tests to verify dependencies:

```bash
uv run pytest
```

## 📫 Contact 
If you have any questions, comments, or contribution ideas, don't hesitate to reach out.
- Anton Appelt, Appelt [at] kinderkrebs-forschung.de
- Laure Ciernik, ciernik [at] tu-berlin.de
- Michael Bockmayr, m.bockmayr [at] uke.de
 
  


> [!IMPORTANT]
> License: PolyForm Noncommercial License 1.0.0
> <https://polyformproject.org/licenses/noncommercial/1.0.0>
