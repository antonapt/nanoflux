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

The inference will:
- Extract methylation features from the prepared `.bed` file (`input_file`)
- Apply the ensemble classifier
- Output predictions and probability scores to the specified `output_directory`
- The `-c` flag indicates that output directories will be created if they do not exist

### 4. ⚖️ Score reads against a pooled CpG atlas (optional)

`nanoflux score` compares every read with a pooled reference of control samples and turns the result into a per-read weight: reads that look like the reference get a low weight, reads that do not get a high weight. Weights are keyed by read name, so they can be computed on one alignment (the atlas genome build, hg38) and applied to another.

Input is a per-read call table aligned to the atlas build, either the `reads.tsv` written by `nanoflux prepare --extract-reads` (modkit) or a nanopolish/f5c `call-methylation` table (`--format nanopolish`). The atlas is a feather file with one row per CpG: chromosome, 1-based plus-strand C position, pooled methylated count and pooled total count.

```bash
nanoflux score -i "$reads_tsv" --atlas "$atlas_feather" -o "$output_directory" -c
```

Shallow atlas sites are read with a credibility weight: a Beta prior `p = (m + alpha) / (n + alpha + beta)` pools the real reads with `alpha` methylated and `beta` unmethylated pseudo-reads. By default `alpha` and `beta` are fitted to the well-covered sites of the atlas (`--prior auto`); pass two numbers (`--prior 1 1`) to set them by hand. The caller's behaviour at methylated and unmethylated sites (four moments) is fitted on the sample; reuse one `moments.json` across a cohort with `--moments`.

Outputs in `output_directory`:
- `read_scores.parquet`: one row per read with `n_cpg` (atlas CpGs on the read), `n_bin`, `z`, `llr_per_call` and `weight`
- `moments.json`: the prior and the fitted caller moments
- `score_summary.json`: counts, per-bin cutoffs and z quantiles

An optional calibration table (`--fit-calibration`, then `--calibration calibration.json`) replaces the prior with an empirical lookup per atlas fraction and coverage. It must be fitted on reads that are not part of the atlas, otherwise the correction cancels out.

To use the weights, run `prepare` as usual on the same sample with `--read-weights read_scores.parquet`. Prepare then builds the pileup itself from the extracted read calls, counting every read by its weight (restricted to the annotation sites, strands combined, calls below the pass threshold counted as failed), and writes `methylation.bed` in the layout `infer` expects. Weights are matched by read name, so the score can come from an hg38 alignment while prepare runs on chm13v2. Reads without a score get `--unscored-weight` (default 1). With all weights equal to 1 the percentages equal modkit's pileup.

```bash
nanoflux prepare -i "$input_file" -o "$output_directory" -c --read-weights "$scores_dir/read_scores.parquet"
nanoflux infer -i "$output_directory/methylation.bed" -o "$results_dir" -c
```

The weight is a plain decreasing function of z with no data-derived cutoff: `w = w_min + (1 - w_min) * sigmoid((center - z) / temperature)`. With the defaults (`--weight-center 0`, `--weight-temperature 1`) a read that fits the atlas gets about 0.5, a read at z = -2 about 0.88 and at z = -4 about 0.98; a smaller temperature sharpens the contrast and `--weight-min` keeps atlas-like reads from vanishing. Only the ratio of weights between reads at the same site matters in the pileup. `--weight-mode hard` with `--weight-thresholds '1:-3,2:-3.8,3-4:-4.6,5-9:-5.8,10+:-8'` gives 1 below the per-bin threshold and `--weight-min` above.

Two scripts run this for a whole cohort laid out as `<samples_dir>/<sample>/{aligned_to_13v2.bam,extracted.tsv}`: [scripts/score_samples.sh](scripts/score_samples.sh) only scores every sample (and documents each step of the scoring), [scripts/run_weighted_cohort.sh](scripts/run_weighted_cohort.sh) scores, runs the weighted prepare and infer, runs an unweighted baseline and writes a `summary.csv` with the top prediction per sample and run.

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
