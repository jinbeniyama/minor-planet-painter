#!/bin/bash
# Obtain orbital elements of minor planets and NEAs,
# MPCORB.DAT (~280 MB) and NEAm00.txt (8 MB), from MPC

SCRIPT_PATH="$(command -v "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
DATA_DIR="$(realpath "$SCRIPT_DIR/../data")"
DATE=$(date -u +%Y%m%d)
MPCORB=${DATA_DIR}/MPCORB${DATE}.DAT
echo "Save MPCORB.DAT as ${MPCORB}"
#wget https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT -O ${MPCORB}
echo "Copy ${MPCORB} as MPCORB.DAT"
cp ${MPCORB} ${DATA_DIR}/MPCORB.DAT

NEA=${DATA_DIR}/NEAm00${DATE}.txt
echo "Save NEAm00.txt as ${NEA}"
#wget https://www.minorplanetcenter.net/iau/MPCORB/NEAm00.txt -O ${NEA}
echo "Copy ${NEA} as NEAm00.txt"
cp ${NEA} ${DATA_DIR}/NEAm00.txt
