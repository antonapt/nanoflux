# 🧬 NanoFlux

DNA methylation ensemble classifier for Oxford Nanopore sequencing data. This python package accompanies the work:
*"Advancing methylation-based brain tumor classification through Nanopore sequencing data of cerebrospinal fluid-derived cell-free DNA"*  by Julia Freese, Anton Appelt, Laure Ciernik et al.

## 📋 Requirements

- **Pixi** (recommended), a cross-platform package manager for the conda ecosystem or a **Conda** distribution of your choice.

    [Pixi Installation Guide](https://pixi.prefix.dev/latest/installation/)

    [Miniforge Installation Guide](https://conda-forge.org/download/)

- A **basecaller** of your choice. We recommend ONT Dorado, see "Usage" for model recommendation. 

    [Doarao Installation Guide](https://software-docs.nanoporetech.com/dorado/latest/)

## 🔧 Installation
### Global:
```bash
pixi global install \
    -c conda-forge \
    -c bioconda \
    --git https://github.com/antonapt/nanoflux
```

Try it with:
```bash
nanoflux --help
```
### In a conda environment:
```bash
conda create -n nanoflux python=3.11 pixi -c conda-forge
conda activate nanoflux
pixi install \
    -c conda-forge \
    -c bioconda \
    --git https://github.com/antonapt/nanoflux
```
Try it with:
```bash
conda activate nanoflux
nanoflux --help
```

## 🚀 Usage

### 1. 🔎 Base calling

We recommend you perform base calling with `dorado` and it's `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` model or newer:

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 "$pod5_path" --modified-bases 5mCG_5hmCG > "$out_BAM"
```

### 2. 📁 Prepare BAM Files

Nanoflux will handle both SAM and BAM files obtained from base calling. As these are often already aligned to some reference sequence, Nanoflux will automatically re-align your files to the T2T CHM13v2.0 reference. If your data happens to be already aligned to this reference, you can pass the `--skip-alignment` flag to the `prepare` command in this step.

⚠️ If you are using a reference for the first time, nanoflux will automatically download it. This can take a few minutes depending on your internet connection. 

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

⚠️ If you are using ```infer``` for the first time, nanoflux will automatically download the model data. This can take a few minutes depending on your internet connection. 

The inference will:
- Extract methylation features from the prepared `.bed` file (`input_file`)
- Apply the ensemble classifier
- Output predictions and probability scores to the specified `output_directory`
- The `-c` flag indicates that output directories will be created if they do not exist

## 📫 Contact 
If you have any questions, comments, or contribution ideas, don't hesitate to reach out.
- Anton Appelt, Appelt [at] kinderkrebs-forschung.de
- Laure Ciernik, ciernik [at] tu-berlin.de
- Michael Bockmayr, m.bockmayr [at] uke.de
 
  


> [!IMPORTANT]
> License: PolyForm Noncommercial License 1.0.0
> <https://polyformproject.org/licenses/noncommercial/1.0.0>
