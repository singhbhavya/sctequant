#!/usr/bin/env bash
# Guess the 10x chemistry of one sample.
# Usage:  ./check_chemistry.sh /path/to/sample_fastq_folder
set -euo pipefail
cd "$(dirname "$0")"
if [ $# -ne 1 ]; then
    echo "Usage: ./check_chemistry.sh /path/to/one/sample/fastq/folder"
    exit 1
fi
ref_dir=$(python -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['reference_dir'])")
python workflow/scripts/check_chemistry.py "$1" "$ref_dir/whitelists"
