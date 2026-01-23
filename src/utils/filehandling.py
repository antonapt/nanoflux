import shutil
import tempfile
from functools import partial, wraps
from pathlib import Path

from src.utils.log import logger


def prepare_location(path: Path, create: bool = False) -> None:
    path = path.resolve()
    parent = path.parent

    if path.exists():
        raise FileExistsError(f"Target already exists: {path}")
    elif not parent.exists():
        if create:
            parent.mkdir(parents=True)
        else:
            raise NotADirectoryError(
                "Target directory not found. Specify -c/--create_dirs to resolve any given path."
            )


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
