from pathlib import Path
from xml.etree import ElementTree

import requests

share_token = "y5ZpkFFbDRirBog"
base_url = "https://tubcloud.tu-berlin.de"
dav_url = f"{base_url}/public.php/webdav/"
password = ""


def download():
    target_dir = Path(__file__).parent / "models"
    target_dir.mkdir(parents=True, exist_ok=True)

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


if __name__ == "__main__":
    download()
