#!/usr/bin/env python3
"""
Download MPC orbital element files.

Download MPCORB.DAT and NEAm00.txt from the Minor Planet Center
and save them in the data directory.
"""

from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlretrieve


def main():
    # Locate and create the data directory.
    data_dir = Path(__file__).resolve().parent.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Current UTC date.
    date = datetime.now(timezone.utc).strftime("%Y%m%d")

    # MPC download URLs and output filenames.
    files = [
        (
            "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT",
            f"MPCORB{date}.DAT",
            "MPCORB.DAT",
        ),
        (
            "https://www.minorplanetcenter.net/iau/MPCORB/NEAm00.txt",
            f"NEAm00{date}.txt",
            "NEAm00.txt",
        ),
    ]

    for url, dated_name, standard_name in files:
        dated_path = data_dir / dated_name
        standard_path = data_dir / standard_name

        print(f"Downloading {url}")
        urlretrieve(url, dated_path)

        print(f"Copy {dated_name} as {standard_name}")
        standard_path.write_bytes(dated_path.read_bytes())

    print(f"Files saved in {data_dir}")


if __name__ == "__main__":
    main()

