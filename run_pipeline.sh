#!/usr/bin/env bash
# Run the sctequant pipeline on Minerva.
#
#   ./run_pipeline.sh check     show what will run, without running anything
#   ./run_pipeline.sh           install software, download references, then
#                               submit the pipeline as an LSF job (you can log out)
#   ./run_pipeline.sh here      run in this terminal (use inside tmux/screen)
#   ./run_pipeline.sh unlock    fix "Directory cannot be locked" after a crash
#
# Add  --test  to any of these to use the small test dataset instead of your
# own data (download it first with ./test/download_test_data.sh).
set -euo pipefail
cd "$(dirname "$0")"

mode="submit"
extra=()
for arg in "$@"; do
    case "$arg" in
        check|submit|here|unlock|prepare) mode="$arg" ;;
        --test) extra+=(--configfile test/config_test.yaml) ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        *) extra+=("$arg") ;;   # anything else is passed to snakemake
    esac
done

die() { echo -e "\nERROR: $*\n" >&2; exit 1; }

if [[ " ${extra[*]:-} " == *" test/config_test.yaml "* ]] && ! ls test/data/*/*_R1_*.fastq.gz >/dev/null 2>&1; then
    die "The test data isn't downloaded yet. Run this first:\n  ./test/download_test_data.sh"
fi

command -v snakemake >/dev/null || die "snakemake not found. Activate your environment first:  conda activate sctequant"
python - <<'PY' || die "Snakemake 8.6 or newer plus the LSF plugin are needed. Install with:\n  pip install 'snakemake>=8.6' snakemake-executor-plugin-lsf"
import importlib.util, sys
from snakemake import __version__ as v
major, minor = (int(x) for x in v.split(".")[:2])
sys.exit(0 if (major, minor) >= (8, 6) and importlib.util.find_spec("snakemake_executor_plugin_lsf") else 1)
PY

cfg() { python -c "import yaml,sys; c=yaml.safe_load(open('config/config.yaml')); print(eval('c'+sys.argv[1]) or '')" "$1"; }
account=$(cfg "['lsf']['account']")
queue=$(cfg "['lsf']['queue']")
long_queue=$(cfg "['lsf']['long_queue']")
max_jobs=$(cfg "['lsf']['max_jobs']")
conda_dir=$(cfg "['conda_envs_dir']")

conda_args=()
[ -n "$conda_dir" ] && conda_args=(--conda-prefix "$conda_dir")

local_args=(--cores 1 --software-deployment-method conda ${conda_args[@]+"${conda_args[@]}"} ${extra[@]+"${extra[@]}"})
lsf_args=(--profile profiles/minerva --jobs "$max_jobs" ${conda_args[@]+"${conda_args[@]}"}
          --default-resources "lsf_project=$account" "lsf_queue=$queue" mem_mb=4000 runtime=120
          ${extra[@]+"${extra[@]}"})

prepare() {
    echo ">>> Step 1/2: installing software (first time only, can take 30+ minutes)"
    snakemake ${local_args[@]+"${local_args[@]}"} --conda-create-envs-only
    echo ">>> Step 2/2: downloading reference files (first time only)"
    snakemake prepare ${local_args[@]+"${local_args[@]}"}
}

case "$mode" in
    check)
        snakemake --cores 1 --dry-run --quiet rules ${extra[@]+"${extra[@]}"}
        echo -e "\nLooks good. Start the pipeline with:  ./run_pipeline.sh ${extra[*]:-}"
        ;;
    unlock)
        snakemake ${local_args[@]+"${local_args[@]}"} --unlock
        ;;
    prepare)
        prepare
        echo -e "\nSoftware and reference downloads are ready. Start the pipeline with:  ./run_pipeline.sh ${extra[*]:-}"
        ;;
    here)
        prepare
        snakemake ${lsf_args[@]+"${lsf_args[@]}"}
        ;;
    submit)
        command -v bsub >/dev/null || die "bsub not found. Are you logged in to Minerva?"
        prepare
        mkdir -p logs
        log="logs/sctequant_$(date +%Y%m%d_%H%M%S).log"
        bsub -J sctequant -P "$account" -q "$long_queue" -n 1 -W 144:00 \
             -R "rusage[mem=4000]" -oo "$log" \
             snakemake ${lsf_args[@]+"${lsf_args[@]}"}
        echo
        echo "Submitted! The pipeline is now running as an LSF job named 'sctequant'."
        echo "  Check on it:   bjobs        (or: tail -f $log)"
        echo "  Stop it:       bkill -J sctequant"
        ;;
esac
