import os
import shutil
import tempfile
from functools import partial, wraps
from pathlib import Path
from subprocess import run

from src.utils.filehandling import prepare_location
from src.utils.log import logger


# seems not to work yet
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
                logger.debug(f"Writing {tmp_path} to {orig_path}")
                shutil.move(tmp_path, orig_path)
            else:
                # attempt to remove tmpfile?
                raise Exception(stderr)

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
    cmd = [
        "minimap2",
        str(reference),
        str(input_path),
        "-x",
        "map-ont",
        "-t",
        str(threads),
    ]
    logger.debug(f"Running command: {' '.join(cmd)} > {output_path}")

    if not dry_run:
        with open(output_path, "w") as outfile:  # Make this a tempfile
            res = run(cmd, stdout=outfile, text=True)
        return res.returncode, res.stdout, res.stderr
    return 0, "", ""


### Change below -----------------------------------------------------------------


def samtools_sort(input_path: Path, threads: int, dry_run: bool = False) -> None:
    cmd = [
        "samtools",
        "sort",
        str(input_path),
        "-@",
        str(threads),
    ]
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        run(
            cmd,
            check=True,
        )


def samtools_index(input_path: Path, threads: int, dry_run: bool = False) -> None:
    cmd = [
        "samtools",
        "index",
        str(input_path),
        "-@",
        str(threads),
    ]
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        run(
            cmd,
            check=True,
        )


def modkit_pileup(
    input_path: Path,
    output_path: Path,
    reference: Path,
    threads: int,
    dry_run: bool = False,
) -> None:
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
    ]
    logger.debug(f"Running command: {' '.join(cmd)}")

    if not dry_run:
        run(
            cmd,
            check=True,
        )
        # if nonzero exit remove file


def bedtools_intersect(
    input_path: Path, output_path: Path, anno_path: Path, dry_run: bool = False
) -> None:
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
            run(
                cmd,
                stdout=outfile,
            )
            # same as for minimap


def main(args):
    input_file = args.input.resolve()
    output_dir = args.output.resolve()
    reference = args.ref.resolve()

    anno_file = Path("src/data/mapping_EPIC.bed").resolve()

    logger.debug(f"Running nanoflux prepare: \n{os.getcwd()}")

    if not args.skip_alignment:
        output_file = output_dir / "aligned_to_CHM13v2.bam"
        prepare_location(output_file, args.create_dir)
        minimap2_align(
            input_path=input_file,
            output_path=output_file,
            reference=reference,
            threads=args.threads,
            dry_run=args.dry_run,
        )
        input_file = output_file

    # samtools_sort(input_file, args.threads, args.dry_run)
    # samtools_index(input_file, args.threads, args.dry_run)

    # pileup_file = output_dir / "pileup.bed"
    # prepare_location(pileup_file, args.create_dir)
    # modkit_pileup(input_file, pileup_file, reference, args.threads, args.dry_run)

    # methyl_file = output_dir / "methylation.bed"
    # prepare_location(methyl_file, args.create_dir)
    # bedtools_intersect(pileup_file, methyl_file, anno_file, args.dry_run)
