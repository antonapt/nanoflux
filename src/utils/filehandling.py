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
            logger.debug(f"Creating dir: {parent}")
            parent.mkdir(parents=True)
        else:
            raise NotADirectoryError(
                "Target directory not found. Specify -c/--create_dirs to resolve any given path."
            )


def with_tmpfile(fun):
    @wraps(fun)
    def wrapper(**kwargs) -> tuple[int, str, str]:
        if kwargs.get("dry_run"):
            logger.debug("Dry-running, not attempting to intercept save path")
            returncode, stdout, stderr = fun(**kwargs)

        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                orig_path = Path(kwargs.pop("output_path"))
                tmp_path = Path(tmpdir) / orig_path.name
                logger.debug(f"Intercepting {orig_path} and redirecting to {tmp_path}")

                p = partial(fun, output_path=tmp_path)
                returncode, stdout, stderr = p(**kwargs)

                logger.debug(stdout)
                logger.debug(stderr)

                if returncode == 0:
                    logger.debug(f"Writing {tmp_path} to {orig_path}")
                    shutil.move(tmp_path, orig_path)
                else:
                    logger.debug("Nonzero return, did not write file")
                    if tmp_path.exists():
                        tmp_path.unlink()

                    raise Exception(
                        f"Tool call failed!\nTool: {fun.__name__}\nStdout: {stdout}\nStderr:{stderr}"
                    )

        return returncode, stdout, stderr

    return wrapper
