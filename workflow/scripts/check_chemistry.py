"""Guess the 10x chemistry of a library from its R1 (barcode) reads.

Usage: python check_chemistry.py <folder or R1.fastq.gz> <whitelist folder>
"""

import glob
import gzip
import os
import subprocess
import sys
from collections import Counter

import yaml

N_READS = 200_000
TSO = "TTTCTTATATGGG"  # template switch oligo, follows barcode+UMI in 5' libraries
HERE = os.path.dirname(os.path.abspath(__file__))


def find_r1(path):
    if os.path.isfile(path):
        return path
    for pattern in ("*_R1_*.fastq.gz", "*_R1.fastq.gz", "*_1.fastq.gz"):
        hits = sorted(glob.glob(os.path.join(path, pattern)))
        if hits:
            return hits[0]
    sys.exit(f"No R1 FASTQ found in {path}")


def read_r1(path):
    seqs = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= N_READS:
                    break
    return seqs


def whitelist_hits(whitelist, barcodes):
    opener = gzip.open if whitelist.endswith(".gz") else open
    found = 0
    with opener(whitelist, "rt") as fh:
        for line in fh:
            if line.strip() in barcodes:
                found += barcodes[line.strip()]
    return found


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    r1_path, wl_dir = find_r1(sys.argv[1]), sys.argv[2]
    with open(os.path.join(HERE, "..", "chemistries.yaml")) as fh:
        chems = yaml.safe_load(fh)

    print(f"Reading the first {N_READS:,} reads of {r1_path}")
    seqs = read_r1(r1_path)
    lengths = Counter(len(s) for s in seqs)
    print(f"R1 read length(s): {', '.join(f'{k} bp' for k, _ in lengths.most_common(3))}")

    barcodes = Counter(s[:16] for s in seqs)
    os.makedirs(wl_dir, exist_ok=True)
    scores = {}
    for wl_name, url in sorted({(c["whitelist"], c["whitelist_url"]) for c in chems.values()}):
        wl = os.path.join(wl_dir, wl_name)
        if not os.path.exists(wl):
            print(f"Downloading whitelist {wl_name} ...")
            subprocess.run(f"curl -fsSL {url} | gunzip -c > {wl}.tmp && mv {wl}.tmp {wl}", shell=True, check=True)
        scores[wl_name] = whitelist_hits(wl, barcodes) / len(seqs)

    print("\nFraction of reads whose first 16 bp is a valid 10x barcode:")
    for name, frac in sorted(scores.items(), key=lambda x: -x[1]):
        print(f"  {name:28s} {100 * frac:5.1f}%")

    tso = {umi: sum(s[16 + umi:16 + umi + len(TSO)] == TSO for s in seqs) / len(seqs) for umi in (10, 12)}
    best = max(scores, key=scores.get)
    is_5prime = max(tso.values()) > 0.2

    print(f"\nTemplate-switch oligo right after barcode+UMI: {100 * max(tso.values()):.1f}% of reads")
    if scores[best] < 0.5:
        print("\nRESULT: unclear. Fewer than half the reads match any 10x whitelist.")
        print("Is this really 10x data, and is this the barcode read (R1)?")
        return

    candidates = [n for n, c in chems.items() if c["whitelist"] == best]
    if len(candidates) > 1:
        candidates = [n for n in candidates if ("5prime" in n) == is_5prime] or candidates
    print(f"\nRESULT: chemistry looks like  {candidates[0]}")
    if candidates[0] == "10x_3prime_v2" and max(lengths) < 40:
        print("(3' v2 and 5' v1/v2 share a whitelist; R1 is too short to see the 5' TSO,")
        print(" so check your kit if this sample might be 5'.)")


if __name__ == "__main__":
    main()
