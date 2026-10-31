#!/usr/bin/env python3
"""
Prepare MPC orbital element files.

Copy the date-stamped MPCORB.DAT and NEAm00.txt files
into the standard filenames in the data directory.
"""

from datetime import datetime, timezone
from pathlib import Path
import shutil


def main():
    # Locate the data directory.
    data_dir = Path(__file__).resolve().parent.parent / "data"

    # Current UTC date.
    date = datetime.now(timezone.utc).strftime("%Y%m%d")

    # MPCORB.DAT
    mpcorb = data_dir / f"MPCORB{date}.DAT"
    print(f"Copy {mpcorb} as MPCORB.DAT")
    shutil.copy2(mpcorb, data_dir / "MPCORB.DAT")

    # NEAm00.txt
    nea = data_dir / f"NEAm00{date}.txt"
    print(f"Copy {nea} as NEAm00.txt")
    shutil.copy2(nea, data_dir / "NEAm00.txt")


if __name__ == "__main__":
    main()

