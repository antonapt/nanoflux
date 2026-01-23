import shutil
import tempfile
from functools import partial, wraps
from pathlib import Path
from subprocess import PIPE, Popen, run

from src.utils.filehandling import prepare_location
from src.utils.log import logger


# intercept output_path and redirect to tmpfile
# if returncode == 0 move tmpfile to output_path, else raise exception
def with_tmpfile(fun):
    @wraps(fun)
    def wrapper(**kwargs):
        with tempfile.TemporaryDirectory() as tmpdir:
            orig_path = Path(kwargs.pop("output_path"))
            tmp_path = Path(tmpdir) / orig_path.name
            logger.debug(f"Intercepting {orig_path} and redirecting to {tmp_path}")

            p = partial(fun, output_path=tmp_path)
            returncode, stdout, stderr = p(**kwargs)

            if returncode == 0:
                logger.debug(stdout)
                logger.debug(stderr)
                logger.debug(f"Writing {tmp_path} to {orig_path}")
                shutil.move(tmp_path, orig_path)
            else:
                if tmp_path.exists():
                    tmp_path.unlink()
                raise Exception(
                    f"Errorcode: {returncode}\nStdout: {stdout}\nStderr: {stderr}"
                )

    return wrapper


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
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        proc = run(cmd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


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
    logger.debug(f"Running command: {' '.join(cmd)} > {output_path}")

    if not dry_run:
        with open(output_path, "w") as outfile:
            proc = run(cmd, stdout=outfile, stderr=PIPE, text=True)
            return proc.returncode, proc.stdout, proc.stderr
    return 0, "", ""


def main(args):
    input_file = args.input.resolve()
    output_dir = args.output.resolve()
    reference = args.ref.resolve()

    anno_file = Path("src/data/mapping_EPIC.bed").resolve()

    if not args.skip_alignment:
        output_file = output_dir / "aligned_to_CHM13v2.sam"
        prepare_location(output_file, args.create_dir)
        minimap2_align(
            input_path=input_file,
            output_path=output_file,
            reference=reference,
            threads=args.threads,
            dry_run=args.dry_run,
        )
        input_file = output_file

    output_file = output_dir / "aligned_to_CHM13v2.bam"
    prepare_location(output_file, args.create_dir)
    samtools_sort(
        input_path=input_file,
        output_path=output_file,
        threads=args.threads,
        dry_run=args.dry_run,
    )
    input_file = output_file

    output_file = output_dir / "aligned_to_CHM13v2.bam.bai"
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
        anno_path=anno_file,
        dry_run=args.dry_run,
    )
