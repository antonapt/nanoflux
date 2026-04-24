import logging
import tempfile
from importlib.resources import files
from pathlib import Path
from subprocess import PIPE, Popen, run

import pandas as pd

from src.utils.filehandling import prepare_location, with_tmpfile
from src.utils.log import logger


@with_tmpfile
def minimap2_align(
    *,  # force keyword args to ensure decorator compatibility
    input_path: Path,
    output_path: Path,
    reference: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    cmd_sam = [
        "samtools",
        "fastq",
        "-@",
        str(threads),
        "-T",
        "MM,ML",
        str(input_path),
    ]
    cmd_mm2 = [
        "minimap2",
        str(reference),
        "-",
        "-o",
        str(output_path),
        "-ayx",
        "map-ont",
        "-t",
        str(threads),
    ]
    logger.info("Running [green bold]minimap2[/]")
    logger.debug(f"Running command: {' '.join(cmd_sam)} | {' '.join(cmd_mm2)}")

    if not dry_run:
        sam = Popen(cmd_sam, stdout=PIPE, stderr=PIPE)
        mm2 = Popen(cmd_mm2, stdin=sam.stdout, stdout=PIPE, stderr=PIPE, text=True)
        if sam.stdout is not None:
            sam.stdout.close()  # allow sam to get SIGPIPE if mm2 exits
        out, err = mm2.communicate()
        sam.wait()
        return mm2.returncode, out, err
    return 0, "", ""


@with_tmpfile
def samtools_sort(
    *,
    input_path: Path,
    output_path: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    cmd = [
        "samtools",
        "sort",
        "-O",
        "BAM",
        "-o",
        str(output_path),
        "-@",
        str(threads),
        str(input_path),
    ]
    logger.info("Running [green bold]samtools sort[/]")
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


def _load_thresholds_row(thresholds_file: Path, quantile: float) -> dict:
    df = pd.read_csv(thresholds_file, index_col=0)
    return df.loc[quantile].to_dict()


def _assign_bin_thresholds(covered_cpgs: pd.Series, thresholds_row: dict) -> pd.Series:
    """Assign per-read threshold based on CpG coverage bins derived from sorted column lower bounds.
    Column keys are integer lower bounds; bins are [col_i, col_{i+1}) with [last_col, inf) for the last.
    """
    boundaries = sorted((int(k), float(v)) for k, v in thresholds_row.items())
    result = pd.Series(float("nan"), index=covered_cpgs.index, dtype=float)
    for i, (lo, thresh) in enumerate(boundaries):
        hi = boundaries[i + 1][0] if i + 1 < len(boundaries) else None
        mask = (
            (covered_cpgs >= lo) & (covered_cpgs < hi)
            if hi is not None
            else (covered_cpgs >= lo)
        )
        result[mask] = thresh
    return result


def get_readnames_from_feather(
    feather_path: Path,
    short_read_min_covered_cpgs: int | None = None,
    non_tumor_score_threshold: float | None = None,
    score_filter_min_cpgs: int | None = None,
    thresholds_file: Path | None = None,
    quantile: float | None = None,
    read_name_col: str = "read_name",
    score_col: str = "nontumor_score_sum",
    cpg_col: str = "covered_cpgs",
) -> pd.Series:
    """Extracts read names from a feather file, applying up to three filter modes:
    1. Short-read filter: drop reads with covered_cpgs < short_read_min_covered_cpgs.
    2. Score filter (single threshold): keep reads with covered_cpgs < score_filter_min_cpgs
       OR score < threshold. Threshold comes from non_tumor_score_threshold directly, or
       is looked up from thresholds_file at the given quantile and score_filter_min_cpgs column.
    3. Bin-based score filter: when thresholds_file + quantile are given without
       score_filter_min_cpgs, each read is assigned a threshold based on its cpg coverage bin.
       Reads with score < bin_threshold are kept (NaN bin → keep).
    Modes 2 and 3 are mutually exclusive.
    """
    feather_path = Path(feather_path).resolve()
    logger.info(f"Loading read filter feather: {feather_path}")
    df = pd.read_feather(feather_path)
    nr_reads_orig = len(df)

    if short_read_min_covered_cpgs is not None:
        df = df.loc[df[cpg_col] >= short_read_min_covered_cpgs]

    if thresholds_file is not None and quantile is not None:
        row = _load_thresholds_row(Path(thresholds_file), quantile)
        if score_filter_min_cpgs is not None:
            # Approach 2: single threshold from CSV column indexed by min-CpG value
            col_key = str(score_filter_min_cpgs)
            if col_key not in row:
                raise ValueError(f"Column '{col_key}' not found in {thresholds_file}")
            threshold = float(row[col_key])
            logger.info(
                f"Score filter (approach 2): threshold={threshold} (quantile={quantile}, min_cpgs={score_filter_min_cpgs})"
            )
            likely_tumor = (df[cpg_col] < score_filter_min_cpgs) | (
                df[score_col] < threshold
            )
        else:
            # Approach 3: per-read threshold from coverage bins
            per_read_thresh = _assign_bin_thresholds(df[cpg_col], row)
            n_binned = per_read_thresh.notna().sum()
            logger.info(
                f"Score filter (approach 3): bin-based thresholds (quantile={quantile}), {n_binned}/{len(df)} reads matched a bin"
            )
            likely_tumor = per_read_thresh.isna() | (df[score_col] < per_read_thresh)
        df = df.loc[likely_tumor]
    elif non_tumor_score_threshold is not None and score_filter_min_cpgs is not None:
        # Approach 2: direct threshold (no CSV lookup)
        logger.info(
            f"Score filter (approach 2): threshold={non_tumor_score_threshold} (direct, min_cpgs={score_filter_min_cpgs})"
        )
        likely_tumor = (df[cpg_col] < score_filter_min_cpgs) | (
            df[score_col] < non_tumor_score_threshold
        )
        df = df.loc[likely_tumor]

    read_names = df.loc[:, read_name_col]
    logger.info(
        f"Keeping {len(read_names)} / {nr_reads_orig} reads after feather-based filtering"
    )
    return read_names


def samtools_filter_ids_from_names(
    *,
    read_names,
    input_path: Path,
    output_path: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    """Filters a SAM/BAM file using an iterable of read names, via a temporary file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        tmp.write("\n".join(read_names))
    try:
        return samtools_filter_ids(
            input_path=input_path,
            output_path=output_path,
            read_list=tmp_path,
            threads=threads,
            dry_run=dry_run,
        )
    finally:
        tmp_path.unlink(missing_ok=True)


@with_tmpfile
def samtools_filter_ids(
    *,
    input_path: Path,
    output_path: Path,
    read_list: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    """Filters a SAM/BAM file using a list of read IDs."""
    if not read_list.exists() or read_list.stat().st_size == 0:
        logger.error(
            f"Read filter list {read_list} does not exist or is empty. Skipping filtering step."
        )
        raise FileNotFoundError(
            f"Read filter list {read_list} does not exist or is empty."
        )

    cmd = [
        "samtools",
        "view",
        "-h",
        "-N",
        str(read_list),
        "-@",
        str(threads),
        "-o",
        str(output_path),
        str(input_path),
    ]
    logger.info("Running [green bold]samtools filter reads[/]")
    logger.debug(f"Running command: {' '.join(cmd)}")
    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


@with_tmpfile
def samtools_index(
    *,
    input_path: Path,
    output_path: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    cmd = [
        "samtools",
        "index",
        "-o",
        str(output_path),
        "-@",
        str(threads),
        str(input_path),
    ]
    logger.info("Running [green bold]samtools index[/]")
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


@with_tmpfile
def modkit_pileup(
    *,
    input_path: Path,
    output_path: Path,
    reference: Path,
    threads: int,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    cmd = [
        "modkit",
        "pileup",
        str(input_path),
        str(output_path),
        "-t",
        str(threads),
        "--preset",
        "traditional",
        "--ref",
        str(reference),
        "--suppress-progress",  # prevent huge stdout
        # "--log-file" # TODO: add a path to this
    ]
    logger.info("Running [green bold]modkit pileup[/]")
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


@with_tmpfile
def bedtools_intersect(
    *,
    input_path: Path,
    output_path: Path,
    anno_path: Path,
    dry_run: bool = False,
) -> tuple[int, str, str]:
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(input_path),
        "-b",
        str(anno_path),
        "-wa",
        "-wb",
    ]
    logger.info("Running [green bold]bedtools intersect[/]")
    logger.debug(f"Running command: {' '.join(cmd)} > {output_path}")

    if not dry_run:
        with open(output_path, "w") as outfile:
            proc = run(cmd, stdout=outfile, stderr=PIPE, text=True)
            return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


def main(args):
    if args.debug or args.dry_run:
        logger.setLevel(logging.DEBUG)

    logger.info("Running [blue bold]nanoflux prepare[/]")

    input_file = args.input.resolve()
    output_dir = args.output.resolve()
    reference = args.ref.resolve()

    anno_file = files("data") / "features" / "mapping_EPIC.bed"

    methyl_file = output_dir / "methylation.bed"
    prepare_location(methyl_file, args.create_dir)

    with tempfile.TemporaryDirectory(prefix="nanoflux_prepare_") as tmpdir:
        work_dir = Path(tmpdir)

        if not args.skip_alignment:
            output_file = work_dir / "aligned_to_CHM13v2.sam"
            minimap2_align(
                input_path=input_file,
                output_path=output_file,
                reference=reference,
                threads=args.threads,
                dry_run=args.dry_run,
            )
            input_file = output_file

        if args.read_filter_list:
            filter_list = Path(args.read_filter_list).resolve()
            filtered_bam = work_dir / "aligned_to_CHM13v2_and_filtered_reads.bam"
            samtools_filter_ids(
                input_path=input_file,
                output_path=filtered_bam,
                read_list=filter_list,
                threads=args.threads,
                dry_run=args.dry_run,
            )
            input_file = filtered_bam

        elif args.read_filter_feather:
            read_name_col, score_col, cpg_col = args.feather_cols
            read_names = get_readnames_from_feather(
                feather_path=Path(args.read_filter_feather).resolve(),
                short_read_min_covered_cpgs=args.short_read_min_covered_cpgs,
                non_tumor_score_threshold=args.non_tumor_score_threshold,
                score_filter_min_cpgs=args.score_filter_min_cpgs,
                thresholds_file=args.thresholds_file,
                quantile=args.quantile,
                read_name_col=read_name_col,
                score_col=score_col,
                cpg_col=cpg_col,
            )
            filtered_bam = work_dir / "aligned_to_CHM13v2_and_filtered_reads.bam"
            samtools_filter_ids_from_names(
                input_path=input_file,
                output_path=filtered_bam,
                read_names=read_names,
                threads=args.threads,
                dry_run=args.dry_run,
            )
            input_file = filtered_bam

        output_file = work_dir / "aligned_to_CHM13v2.bam"
        samtools_sort(
            input_path=input_file,
            output_path=output_file,
            threads=args.threads,
            dry_run=args.dry_run,
        )
        input_file = output_file

        output_file = work_dir / "aligned_to_CHM13v2.bam.bai"
        samtools_index(
            input_path=input_file,
            output_path=output_file,
            threads=args.threads,
            dry_run=args.dry_run,
        )

        pileup_file = work_dir / "pileup.bed"
        modkit_pileup(
            input_path=input_file,
            output_path=pileup_file,
            reference=reference,
            threads=args.threads,
            dry_run=args.dry_run,
        )

        bedtools_intersect(
            input_path=pileup_file,
            output_path=methyl_file,
            anno_path=anno_file,  # type: ignore
            dry_run=args.dry_run,
        )

    logger.info("[blue bold]nanoflux prepare[/] has successfully run!")
    logger.info(f"The intermediate file has been saved to: {methyl_file}")
