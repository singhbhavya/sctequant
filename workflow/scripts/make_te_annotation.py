"""Turn a UCSC RepeatMasker .fa.out.gz file into the TE annotations each tool needs.

Outputs
  --subfamily-bed  6-column BED, name = TE subfamily (for scTE and IRescue)
  --solote-bed     BED in SoloTE's format (same as SoloTE_RepeatMasker_to_BED.py -i)
  --families       table: subfamily, family, class (used to annotate AnnData .var)
"""

import argparse
import gzip

# Repeat classes counted as transposable elements.
TE_CLASSES = {"LINE", "SINE", "LTR", "DNA", "RC", "Retroposon"}


def main_chromosome(chrom):
    return not any(tag in chrom for tag in ("_random", "chrUn", "_alt", "_fix"))


def parse(args):
    families = {}
    with gzip.open(args.rmsk, "rt") as rmsk, \
            open(args.subfamily_bed, "w") as sub_bed, \
            open(args.solote_bed, "w") as solo_bed:
        for line in rmsk:
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue  # header or blank line
            div, chrom, start, end, strand = fields[1], fields[4], fields[5], fields[6], fields[8]
            subfamily, class_family = fields[9], fields[10]
            if not main_chromosome(chrom):
                continue
            strand = "-" if strand == "C" else "+"
            parts = class_family.split("/")
            te_class = parts[0].rstrip("?")
            te_family = parts[1].rstrip("?") if len(parts) > 1 else "-"
            if te_class not in TE_CLASSES:
                continue

            sub_bed.write(f"{chrom}\t{int(start) - 1}\t{end}\t{subfamily}\t0\t{strand}\n")

            # SoloTE's own format: 1-based start and a composite name.
            solo_family = parts[1] if len(parts) > 1 else "-"
            te_name = f"{chrom}|{start}|{end}|{subfamily}:{solo_family}:{parts[0]}|{div}|{strand}"
            solo_bed.write(f"{chrom}\t{start}\t{end}\t{te_name}\t{div}\t{strand}\n")

            families.setdefault(subfamily, (te_family, te_class))

    with open(args.families, "w") as out:
        out.write("subfamily\tfamily\tclass\n")
        for subfamily in sorted(families):
            out.write(f"{subfamily}\t{families[subfamily][0]}\t{families[subfamily][1]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rmsk", required=True)
    parser.add_argument("--subfamily-bed", required=True)
    parser.add_argument("--solote-bed", required=True)
    parser.add_argument("--families", required=True)
    parse(parser.parse_args())
