import logging
from importlib.resources import files
from pathlib import Path
from subprocess import PIPE, Popen, run

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
def samtools_filter(
    *,
    input_path: Path,
    output_path: Path,
    threads: int,
    min_mapq: int = 0,  # mapq > 0 discards multimapped reads
    dry_run: bool = False,
) -> tuple[int, str, str]:

    cmd = [
        "samtools",
        "view",
        "-b",
        "-o",
        str(output_path),
        "-@",
        str(threads),
        "-q",
        str(min_mapq),
        str(input_path),
    ]
    logger.info("Running [green bold]samtools view[/] (mapq filtering)")
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
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

    input_file = args.input.resolve()
    output_dir = args.output.resolve()

    reference = Path(str(files("data") / "refs" / f"{args.ref}.fa"))
    annotation = Path(str(files("data") / "features" / f"EPIC_{args.ref}.bed"))

    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found at {reference}")

    logger.info("Running [blue bold]nanoflux prepare[/]")

    if not args.skip_alignment:
        output_file = output_dir / f"aligned_to_{args.ref}.sam"
        prepare_location(output_file, args.create_dir)
        minimap2_align(
            input_path=input_file,
            output_path=output_file,
            reference=reference,
            threads=args.threads,
            dry_run=args.dry_run,
        )
        input_file = output_file

    # if sort and index are skipped mapq filtering can not be applied. This gets checked in the parser
    if not args.skip_sort_index:
        filename = f"aligned_to_{args.ref}_q{args.min_mapq}.bam"

        if args.min_mapq > 0:
            prepare_location(output_dir / filename, args.create_dir)
            samtools_filter(
                input_path=input_file,
                output_path=output_dir / filename,
                threads=args.threads,
                min_mapq=args.min_mapq,
                dry_run=args.dry_run,
            )
            input_file = output_dir / filename
        elif args.min_mapq == 0:
            logger.warning(
                "[red bold]min_mapq is set to 0[/], which means multimapped reads will be included. This may lead to inaccurate methylation calls in repetitive regions."
            )

        output_file = output_dir / filename
        prepare_location(
            output_file, args.create_dir, allow_overwrite=(args.min_mapq > 0)
        )
        samtools_sort(
            input_path=input_file,
            output_path=output_file,
            threads=args.threads,
            dry_run=args.dry_run,
        )
        input_file = output_file

        output_file = output_dir / f"{filename}.bai"
        prepare_location(output_file, args.create_dir)
        samtools_index(
            input_path=input_file,
            output_path=output_file,
            threads=args.threads,
            dry_run=args.dry_run,
        )

    pileup_file = output_dir / "pileup.bed"
    prepare_location(pileup_file, args.create_dir)
    modkit_pileup(
        input_path=input_file,
        output_path=pileup_file,
        reference=reference,
        threads=args.threads,
        dry_run=args.dry_run,
    )

    methyl_file = output_dir / "methylation.bed"
    prepare_location(methyl_file, args.create_dir)
    bedtools_intersect(
        input_path=pileup_file,
        output_path=methyl_file,
        anno_path=annotation,  # type: ignore
        dry_run=args.dry_run,
    )

    logger.info("[blue bold]nanoflux prepare[/] has successfully run!")
    logger.info(f"The intermediate file has been saved to: {methyl_file}")
