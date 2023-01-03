#!/bin/bash

# This script downloads the latest data (as a JSON file)
# of the orbits of all asteroids from the Minor Planet Center
# (MPC) database

SCRIPT_PATH="$(dirname "$(realpath -s "$0")")"
DATA_DIR="${SCRIPT_PATH}/data"
DATA_ARCHIVE="${DATA_DIR}/data.json.gz"

# Make sure the data directory exists
mkdir -p "${SCRIPT_PATH}/data"

# Download the latest archive
curl -o "${DATA_ARCHIVE}" "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz"

# Unpack the archive
gzip -d "${DATA_ARCHIVE}"
