import gzip
import shutil
from xml.etree import ElementTree

import requests
from platformdirs import user_data_path
from rich.live import Live
from rich.spinner import Spinner

share_token = "y5ZpkFFbDRirBog"
base_url = "https://tubcloud.tu-berlin.de"
dav_url = f"{base_url}/public.php/webdav/"
password = ""


def download_models():
    target_dir = user_data_path("nanoflux") / "models"

    # Adding 'Depth: 2' allows us to see inside the subfolders
    response = requests.request(
        "PROPFIND", dav_url, auth=(share_token, password), headers={"Depth": "2"}
    )

    if response.status_code != 207:
        print(f"Error: {response.status_code}")
        return

    tree = ElementTree.fromstring(response.content)
    ns = {"d": "DAV:"}

    for resp in tree.findall("d:response", ns):
        href = resp.find("d:href", ns).text

        # If it ends in a slash, it's a folder. We skip folders for downloading.
        if href.endswith("/"):
            print(f"Found folder: {href}")
            continue

        file_name = href.split("/")[-1]
        dest_path = target_dir / file_name

        if dest_path.exists():
            print(f"Skipping existing file: {dest_path}")
            continue

        download_url = f"{base_url}{href}"

        print(f"Downloading: {file_name}...")
        r = requests.get(download_url, auth=(share_token, password))
        with open(dest_path, "wb") as f:
            f.write(r.content)


def download_genome(name: str) -> None:
    match name:
        case "chm13v2.fa":
            genome_url = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
        case "hg38.fa":
            genome_url = (
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
            )
        case "hg38-as.fa":
            genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
        case _:
            raise ValueError(f"Unknown genome name: {name}")

    target_dir = user_data_path("nanoflux") / "refs"
    target_dir.mkdir(parents=True, exist_ok=True)

    target_path = target_dir / name
    temporary_path = target_path.with_suffix(target_path.suffix + ".part")

    spinner = Spinner(
        "dots",
        text=f"Downloading {name} to {target_path}...",
    )

    try:
        with (
            requests.get(genome_url, stream=True, allow_redirects=True) as response,
            Live(spinner, refresh_per_second=10),
        ):
            response.raise_for_status()
            response.raw.decode_content = False

            with (
                gzip.GzipFile(fileobj=response.raw) as decompressed,
                temporary_path.open("wb") as output,
            ):
                shutil.copyfileobj(decompressed, output, length=1024 * 1024)

        temporary_path.replace(target_path)
    finally:
        temporary_path.unlink(missing_ok=True)
