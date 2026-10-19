######################################################################
## Shared setup: read the config, chemistry presets and sample sheet,
## find each sample's FASTQs, and stop early with a clear message if
## anything looks wrong.
######################################################################

import glob
import os
import re

import pandas as pd
import yaml

OUT = config["output_dir"].rstrip("/")
REF = config["reference_dir"].rstrip("/")
GENOME = config["genome"]["name"]

with open(os.path.join(workflow.basedir, "chemistries.yaml")) as fh:
    CHEMISTRIES = yaml.safe_load(fh)

TOOLS = [tool for tool, on in config["tools"].items() if on]


def check_placeholders():
    left = [key for key, value in (
        ("fastq_dir", config["fastq_dir"]),
        ("reference_dir", config["reference_dir"]),
        ("lsf: account", config["lsf"]["account"]),
    ) if "CHANGE_ME" in str(value)]
    if left:
        fail("Edit config/config.yaml first. Still set to CHANGE_ME: " + ", ".join(left))


def fail(message):
    raise ValueError("\n\n" + "=" * 70 + "\n  PROBLEM: " + message + "\n" + "=" * 70 + "\n")


def load_samples():
    sheet = config["sample_sheet"]
    if not os.path.exists(sheet):
        fail(f"Sample sheet not found: {sheet}\n  Check 'sample_sheet' in config/config.yaml.")

    df = pd.read_csv(sheet, dtype=str, sep=None, engine="python").fillna("")
    df.columns = df.columns.str.strip()
    cols = config["sample_sheet_columns"]

    if cols["sample"] not in df.columns:
        fail(
            f"The sample sheet has no '{cols['sample']}' column.\n"
            f"  Columns found: {', '.join(df.columns)}\n"
            "  Rename a column, or set sample_sheet_columns -> sample in config/config.yaml."
        )

    for column, wanted in (config.get("sample_filter") or {}).items():
        if column not in df.columns:
            fail(f"sample_filter uses column '{column}', which is not in the sample sheet.")
        wanted = [wanted] if isinstance(wanted, str) else [str(w) for w in wanted]
        df = df[df[column].isin(wanted)]

    if df.empty:
        fail("No samples left to run. Check the sample sheet and sample_filter.")

    df = df.rename(columns={cols["sample"]: "sample"})
    df["sample"] = df["sample"].str.strip()

    bad = df.loc[~df["sample"].str.match(r"^[A-Za-z0-9._-]+$"), "sample"].tolist()
    if bad:
        fail(f"Sample names may only use letters, numbers, '.', '_' and '-'. Fix: {bad}")
    dups = df.loc[df["sample"].duplicated(), "sample"].unique().tolist()
    if dups:
        fail(f"Sample names must be unique. Repeated: {dups}")

    if "chemistry" not in df.columns or (df["chemistry"] == "").all():
        df["chemistry"] = config["chemistry"]
    df.loc[df["chemistry"] == "", "chemistry"] = config["chemistry"]
    unknown = sorted(set(df["chemistry"]) - set(CHEMISTRIES))
    if unknown:
        fail(f"Unknown chemistry {unknown}. Choose from: {', '.join(CHEMISTRIES)}")

    folder_col = cols.get("folder")
    if folder_col and folder_col in df.columns and folder_col != "folder":
        df = df.rename(columns={folder_col: "folder"})
    if "folder" not in df.columns:
        df["folder"] = ""

    return df.set_index("sample", drop=False)


R1_PATTERNS = ("*_R1_*.fastq.gz", "*_R1.fastq.gz", "*_1.fastq.gz", "*_R1_*.fq.gz", "*_R1.fq.gz", "*_1.fq.gz")


def find_fastqs(sample, folder):
    """Return matched (R1 list, R2 list) for one sample, sorted by lane."""
    base = config["fastq_dir"]
    if folder:
        search_dir, prefix = os.path.join(base, folder), ""
    elif os.path.isdir(os.path.join(base, sample)):
        search_dir, prefix = os.path.join(base, sample), ""
    else:
        search_dir, prefix = base, sample + "_"

    if not os.path.isdir(search_dir):
        fail(f"FASTQ folder for sample '{sample}' does not exist:\n  {search_dir}")

    r1 = sorted({f for p in R1_PATTERNS for f in glob.glob(os.path.join(search_dir, prefix + p))})
    if not r1:
        fail(
            f"No R1 FASTQ files found for sample '{sample}' in:\n  {search_dir}\n"
            "  Expected names like  <name>_S1_L001_R1_001.fastq.gz"
        )

    r2 = []
    for f in r1:
        name = os.path.basename(f)
        mate = re.sub(r"_R1(_|\.)", r"_R2\1", name, count=1)
        if mate == name:
            mate = re.sub(r"_1(\.f(ast)?q\.gz)$", r"_2\1", name)
        mate = os.path.join(os.path.dirname(f), mate)
        if not os.path.exists(mate):
            fail(f"Found {f}\n  but not its partner R2 file:\n  {mate}")
        r2.append(mate)
    return r1, r2


check_placeholders()
SAMPLES = load_samples()
FASTQS = {s: find_fastqs(s, SAMPLES.loc[s, "folder"]) for s in SAMPLES.index}


def chem(sample):
    return CHEMISTRIES[SAMPLES.loc[sample, "chemistry"]]


def whitelist_path(sample):
    return f"{REF}/whitelists/{chem(sample)['whitelist']}"


wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES.index),
    tool="|".join(TOOLS) if TOOLS else "none",
