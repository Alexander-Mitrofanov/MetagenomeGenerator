#!/bin/bash
# Same 5×5 matrix as run_benchmark_25_matrix.sh (5 lengths × 5 split seeds) but:
#   - No genome pool / NCBI: creates one tiny FASTA per category under mock_input_genomes/
#   - Tiny READS_PER_ORG by default for quick timing checks
#
# Still requires a real viral BLAST DB (--viral-db) for blastn-filter.
#
# Env overrides:
#   BENCH_ROOT (default: output/benchmark_25_matrix_mock)
#   MOCK_GENOME_LEN (default: 20000 bp per genome, single record each)
#   READS_PER_ORG (default: 5)
#   CHUNK_SEED, TRAIN_PCT, VIRAL_DB, SNAP unused for pool (no snapshot check)

set -euo pipefail
export PYTHONUNBUFFERED=1
if command -v stdbuf >/dev/null 2>&1; then
  MG() { stdbuf -oL -eL metagenome-generator "$@"; }
else
  MG() { metagenome-generator "$@"; }
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VIRAL_DB="${VIRAL_DB:-${ROOT}/viral_reference/viral_db_2026-03-10/blastn_db/viral_db}"
BENCH_ROOT="${BENCH_ROOT:-${ROOT}/output/benchmark_25_matrix_mock}"
CHUNK_SEED="${CHUNK_SEED:-1}"
READS_PER_ORG="${READS_PER_ORG:-5}"
TRAIN_PCT="${TRAIN_PCT:-80}"
MOCK_GENOME_LEN="${MOCK_GENOME_LEN:-20000}"

LENGTHS=(150 250 500 1000 3000)
SPLIT_SEEDS=(1 2 3 4 5)

LOG_DIR="${BENCH_ROOT}/logs"
mkdir -p "${LOG_DIR}"
RUN_LOG="${LOG_DIR}/benchmark_25_mock_$(date +%Y%m%d_%H%M%S).log"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "${RUN_LOG}"; }

TOTAL_START=$SECONDS
log "=== benchmark 25 MOCK: 5 lengths × 5 split seeds (1 FASTA per category) ==="
log "BENCH_ROOT=${BENCH_ROOT}"
log "MOCK_GENOME_LEN=${MOCK_GENOME_LEN} | reads-per-organism=${READS_PER_ORG} | train-test-split=${TRAIN_PCT}%"

if [[ ! -e "${VIRAL_DB}" && ! -e "${VIRAL_DB}.nhr" ]]; then
  log "ERROR: viral DB not found: ${VIRAL_DB}"
  exit 1
fi

GENOME_DIR="${BENCH_ROOT}/mock_input_genomes"
BLAST_DIR="${BENCH_ROOT}/blastn_shared"
EVE_QUERY_STORE="${BENCH_ROOT}/eve_query_store"
mkdir -p "${BENCH_ROOT}" "${EVE_QUERY_STORE}"

# --- 1) One synthetic genome per category (no NCBI) ---
STEP_START=$SECONDS
log "MOCK | writing 4 FASTAs -> ${GENOME_DIR}"
python3 - "${GENOME_DIR}" "${MOCK_GENOME_LEN}" << 'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1])
n = int(sys.argv[2])
seq = ("ACGT" * ((n // 4) + 1))[:n]
specs = [
    ("bacteria", "mock_bacteria"),
    ("virus", "mock_virus"),
    ("archaea", "mock_archaea"),
    ("plasmid", "mock_plasmid"),
]
for sub, acc in specs:
    d = root / sub
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{acc}.fasta").write_text(f">{acc}\n{seq}\n")
print("OK", root)
PY
log "MOCK | genomes ready in $((SECONDS - STEP_START))s"

# --- 2) BLAST / EVE once ---
STEP_START=$SECONDS
log "BLASTN | ${BLAST_DIR} | eve-query-store=${EVE_QUERY_STORE}"
MG blastn-filter \
  --genome-dir "${GENOME_DIR}" \
  --out-dir "${BLAST_DIR}" \
  --viral-db "${VIRAL_DB}" \
  --eve-query-store "${EVE_QUERY_STORE}" 2>&1 | tee -a "${RUN_LOG}"
log "BLASTN | done in $((SECONDS - STEP_START))s"

EVE_JSON="${BLAST_DIR}/eve_intervals.json"

# --- 3) Five chunk passes; each length then five splits ---
for L in "${LENGTHS[@]}"; do
  LDIR="${BENCH_ROOT}/length_${L}nt"
  mkdir -p "${LDIR}"
  STEP_START=$SECONDS
  log "CHUNK | length=${L}nt -> ${LDIR}/metagenome.fasta"
  MG chunk \
    --input "${GENOME_DIR}" \
    --output metagenome.fasta \
    --output-dir "${LDIR}" \
    --sequence-length "${L}" \
    --reads-per-organism "${READS_PER_ORG}" \
    --balanced \
    --eve-intervals "${EVE_JSON}" \
    --seed "${CHUNK_SEED}" 2>&1 | tee -a "${RUN_LOG}"
  log "CHUNK | length=${L}nt | done in $((SECONDS - STEP_START))s"

  META="${LDIR}/metagenome.fasta"
  if [[ ! -f "${META}" ]]; then
    log "ERROR: missing ${META}"
    exit 1
  fi

  for S in "${SPLIT_SEEDS[@]}"; do
    SDIR="${LDIR}/split_seed_${S}"
    mkdir -p "${SDIR}"
    STEP_START=$SECONDS
    log "SPLIT | length=${L}nt split_seed=${S} -> ${SDIR}"
    MG split-metagenome-train-test \
      --input "${META}" \
      --output-dir "${SDIR}" \
      --train-test-split "${TRAIN_PCT}" \
      --seed "${S}" 2>&1 | tee -a "${RUN_LOG}"
    log "SPLIT | length=${L}nt seed=${S} | done in $((SECONDS - STEP_START))s"
  done
done

log "=== ALL 25 MOCK RUNS COMPLETE ==="
log "Total wall time: $(( (SECONDS - TOTAL_START) / 60 ))m $(( (SECONDS - TOTAL_START) % 60 ))s"
log "Log: ${RUN_LOG}"
