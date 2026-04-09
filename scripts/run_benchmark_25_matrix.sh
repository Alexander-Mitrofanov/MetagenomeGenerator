#!/bin/bash
set -euo pipefail
export PYTHONUNBUFFERED=1
if command -v stdbuf >/dev/null 2>&1; then
  MG() { stdbuf -oL -eL metagenome-generator "$@"; }
else
  MG() { metagenome-generator "$@"; }
fi

# 25 runs = 5 read lengths × 5 train/test split seeds.
# Shared work (once per full execution of this script):
#   - genome-pool prepare (NCBI download into pool, skipped if pool_manifest.json exists)
#   - genome-pool materialize → one genome dir (symlinks into pool)
#   - blastn-filter once → eve_intervals.json + shared --eve-query-store (incremental across reruns)
# Per length: one chunk → metagenome.fasta
# Per (length, split_seed): split-metagenome-train-test → metagenome_train.fasta / metagenome_test.fasta
#
# Outputs:
#   ${BENCH_ROOT}/length_${L}nt/metagenome.fasta
#   ${BENCH_ROOT}/length_${L}nt/split_seed_${S}/metagenome_{train,test}.fasta
#
# Env overrides:
#   BENCH_ROOT, MATERIALIZE_SEED (genome subset RNG; default 1), CHUNK_SEED (default 1),
#   READS_PER_ORG (default 1000), TRAIN_PCT (default 80), SNAP, VIRAL_DB

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAP="${SNAP:-${ROOT}/snapshots/accession_snapshot_2026-03-10.json}"
VIRAL_DB="${VIRAL_DB:-${ROOT}/viral_reference/viral_db_2026-03-10/blastn_db/viral_db}"
BENCH_ROOT="${BENCH_ROOT:-${ROOT}/output/benchmark_25_matrix}"
MATERIALIZE_SEED="${MATERIALIZE_SEED:-1}"
CHUNK_SEED="${CHUNK_SEED:-1}"
READS_PER_ORG="${READS_PER_ORG:-1000}"
TRAIN_PCT="${TRAIN_PCT:-80}"

MAX_B="${MAX_BACTERIA:-500}"
MAX_V="${MAX_VIRUS:-500}"
MAX_A="${MAX_ARCHAEA:-100}"
MAX_P="${MAX_PLASMID:-100}"
POOL_SEED="${POOL_SEED:-0}"

LENGTHS=(150 250 500 1000 3000)
SPLIT_SEEDS=(1 2 3 4 5)

LOG_DIR="${BENCH_ROOT}/logs"
mkdir -p "${LOG_DIR}"
RUN_LOG="${LOG_DIR}/benchmark_25_matrix_$(date +%Y%m%d_%H%M%S).log"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "${RUN_LOG}"; }

TOTAL_START=$SECONDS
log "=== benchmark 25 matrix: 5 lengths × 5 split seeds ==="
log "BENCH_ROOT=${BENCH_ROOT}"
log "POOL max: bacteria=${MAX_B} virus=${MAX_V} archaea=${MAX_A} plasmid=${MAX_P} pool-seed=${POOL_SEED}"
log "Materialize sample-seed=${MATERIALIZE_SEED} (same genome set for all 25 runs)"
log "Chunk seed=${CHUNK_SEED} | reads-per-organism=${READS_PER_ORG} | train-test-split=${TRAIN_PCT}%"

if [[ ! -f "${SNAP}" ]]; then
  log "ERROR: snapshot not found: ${SNAP}"
  exit 1
fi
if [[ ! -e "${VIRAL_DB}" && ! -e "${VIRAL_DB}.nhr" ]]; then
  log "ERROR: viral DB not found: ${VIRAL_DB}"
  exit 1
fi

POOL_DIR="${BENCH_ROOT}/genome_pool"
GENOME_DIR="${BENCH_ROOT}/shared_genomes"
BLAST_DIR="${BENCH_ROOT}/blastn_shared"
EVE_QUERY_STORE="${BENCH_ROOT}/eve_query_store"
mkdir -p "${BENCH_ROOT}" "${POOL_DIR}" "${EVE_QUERY_STORE}"

# --- 1) Pool (one download) ---
STEP_START=$SECONDS
log "POOL | prepare | ${POOL_DIR}"
if [[ -f "${POOL_DIR}/pool_manifest.json" ]]; then
  log "POOL | prepare | skipped"
else
  MG genome-pool prepare \
    --accessions-file "${SNAP}" \
    --pool-dir "${POOL_DIR}" \
    --max-bacteria "${MAX_B}" --max-virus "${MAX_V}" \
    --max-archaea "${MAX_A}" --max-plasmid "${MAX_P}" \
    --pool-seed "${POOL_SEED}" 2>&1 | tee -a "${RUN_LOG}"
fi
log "POOL | prepare | done in $((SECONDS - STEP_START))s"

# --- 2) One materialized genome dir (symlinks) ---
STEP_START=$SECONDS
log "MATERIALIZE | ${GENOME_DIR} (sample-seed=${MATERIALIZE_SEED})"
MG genome-pool materialize \
  --pool-dir "${POOL_DIR}" \
  --output-dir "${GENOME_DIR}" \
  --sample-seed "${MATERIALIZE_SEED}" \
  --max-bacteria "${MAX_B}" --max-virus "${MAX_V}" \
  --max-archaea "${MAX_A}" --max-plasmid "${MAX_P}" 2>&1 | tee -a "${RUN_LOG}"
log "MATERIALIZE | done in $((SECONDS - STEP_START))s"

# --- 3) BLAST / EVE once for this genome dir ---
STEP_START=$SECONDS
log "BLASTN | ${BLAST_DIR} | eve-query-store=${EVE_QUERY_STORE}"
MG blastn-filter \
  --genome-dir "${GENOME_DIR}" \
  --out-dir "${BLAST_DIR}" \
  --viral-db "${VIRAL_DB}" \
  --eve-query-store "${EVE_QUERY_STORE}" 2>&1 | tee -a "${RUN_LOG}"
log "BLASTN | done in $((SECONDS - STEP_START))s"

EVE_JSON="${BLAST_DIR}/eve_intervals.json"

# --- 4) Five chunk passes; each length then five splits ---
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

log "=== ALL 25 RUNS COMPLETE ==="
log "Total wall time: $(( (SECONDS - TOTAL_START) / 60 ))m $(( (SECONDS - TOTAL_START) % 60 ))s"
log "Log: ${RUN_LOG}"
