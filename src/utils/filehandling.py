from pathlib import Path


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
