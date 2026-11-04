import os
from argparse import ArgumentTypeError
from collections.abc import Callable
from pathlib import Path


def min_max_float(min_: float, max_: float) -> Callable:
    def _func(arg) -> float:
        try:
            f = float(arg)
        except ValueError:
            raise ArgumentTypeError("Must be of type float")

        if f <= min_ or f >= max_:
            msg = f"value needs to be between {min_} and {max_}"
            raise ArgumentTypeError(msg)

        return f

    return _func


def output_path_type(arg: str) -> Path:
    """
    Custom type for argparse to validate output path.
    Ensures the path is writable or creatable.
    """

    path = Path(arg)
    if path.exists() and not path.is_dir():
        raise ArgumentTypeError(f"{arg} exists and is not a directory.")
    if path.exists() and not os.access(path, os.W_OK):
        raise ArgumentTypeError(f"{arg} is not writable.")
    return path


def bam_sam_path_type(arg: str) -> Path:
    """
    Custom type for argparse to validate input path.
    Ensures the path exists and is readable.
    """

    path = Path(arg)
    if not path.exists():
        raise ArgumentTypeError(f"{arg} does not exist.")
    if not path.is_file():
        raise ArgumentTypeError(f"{arg} is not a file.")
    if not os.access(path, os.R_OK):
        raise ArgumentTypeError(f"{arg} is not readable.")
    if not path.suffix.lower() in [".sam", ".bam"]:
        raise ArgumentTypeError(f"{arg} must have a .sam or .bam extension.")
    return path


def bed_path_type(arg: str) -> Path:
    """
    Custom type for argparse to validate BED file path.
    Ensures the path exists and is readable.
    """

    path = Path(arg)
    if not path.exists():
        raise ArgumentTypeError(f"{arg} does not exist.")
    if not path.is_file():
        raise ArgumentTypeError(f"{arg} is not a file.")
    if not os.access(path, os.R_OK):
        raise ArgumentTypeError(f"{arg} is not readable.")
    if not path.suffix.lower() == ".bed":
        raise ArgumentTypeError(f"{arg} must have a .bed extension.")
    return path
