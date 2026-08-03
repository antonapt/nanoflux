from pathlib import Path

from zenodo_get import download

if __name__ == "__main__":
    download(
        "10.5281/zenodo.18743210",
        output_dir=Path(__file__).parent / "models_test",
    )
