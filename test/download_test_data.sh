#!/usr/bin/env bash
# Download a public 10x Genomics mouse dataset (1k brain cells, E18,
# 3' v3) and keep only the first 2 million reads per lane, so a test run
# takes about an hour instead of a day.
set -euo pipefail
cd "$(dirname "$0")"

URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar"
READS=2000000
out="data/neuron_1k_v3"

if ls "$out"/*R1_001.fastq.gz >/dev/null 2>&1; then
    echo "Test data already present in test/$out"
    exit 0
fi

mkdir -p data/_download "$out"
echo "Downloading test data (~7 GB, a few minutes) ..."
curl -fL --retry 3 "$URL" | tar -x -C data/_download

echo "Keeping the first $READS reads of each file ..."
for f in data/_download/*/*_R[12]_001.fastq.gz; do
    { gzip -dc "$f" || true; } | head -n $((READS * 4)) | gzip > "$out/$(basename "$f")"
done
rm -rf data/_download
echo "Done. Test FASTQs are in test/$out"
ls -lh "$out"
