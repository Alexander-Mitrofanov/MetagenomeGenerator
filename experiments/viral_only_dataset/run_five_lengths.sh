#!/bin/bash
# Build two "viral-heavy" metagenomes (prok vs euk inferred host slices) at 5 read lengths.
# Requires: metagenome-generator on PATH, ENTREZ_EMAIL, optional ENTREZ_API_KEY.
# CHIMERA requires ≥1 non-viral FASTA — we add a long synthetic bacteria/ decoy (~250 kb)
# and use --abundance-profile bacteria=0,virus=1 so the decoy gets 1 read per file, viruses scale.
#
# Env (defaults tuned so viral read count ≥ 5000 per metagenome when MAX_VIRUS is large enough):
#   REPO_ROOT, OUT, MAX_VIRUS_PROK, MAX_VIRUS_EUK, READS_PER_ORG, SAMPLE_SEED, SNAP_PROK, SNAP_EUK

set -euo pipefail
export PYTHONUNBUFFERED=1
if command -v stdbuf >/dev/null 2>&1; then
  MG() { stdbuf -oL -eL metagenome-generator "$@"; }
else
  MG() { metagenome-generator "$@"; }
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${OUT:-${REPO_ROOT}/output/viral_only_dataset}"
MAX_VIRUS_PROK="${MAX_VIRUS_PROK:-2500}"
MAX_VIRUS_EUK="${MAX_VIRUS_EUK:-1876}"
READS_PER_ORG="${READS_PER_ORG:-3}"
SAMPLE_SEED="${SAMPLE_SEED:-42}"
SNAP_PROK="${SNAP_PROK:-${REPO_ROOT}/experiments/viral_host_stratify/out/accession_snapshot_viral_prok.json}"
SNAP_EUK="${SNAP_EUK:-${REPO_ROOT}/experiments/viral_host_stratify/out/accession_snapshot_viral_euk.json}"
DECOY_LEN="${DECOY_LEN:-250000}"

LENGTHS=(150 250 500 1000 3000)

write_decoy() {
  local dest_dir="$1"
  mkdir -p "${dest_dir}/bacteria"
  python3 - "${dest_dir}/bacteria/synthetic_nonviral_decoy.fasta" "${DECOY_LEN}" << 'PY'
import sys
from pathlib import Path
path, n = Path(sys.argv[1]), int(sys.argv[2])
seq = ("ACGT" * ((n // 4) + 1))[:n]
path.write_text(">synthetic_nonviral_decoy\n" + seq + "\n")
print("Wrote", path)
PY
}

run_one_panel() {
  local label="$1" snap="$2" maxv="$3" genomes="$4"
  echo "=== ${label}: download max-virus=${maxv} -> ${genomes} ==="
  mkdir -p "${genomes}/virus" "${genomes}/bacteria" "${genomes}/archaea" "${genomes}/plasmid"
  write_decoy "${genomes}"
  MG download \
    --accessions-file "${snap}" \
    --max-bacteria 0 \
    --max-archaea 0 \
    --max-plasmid 0 \
    --max-virus "${maxv}" \
    --sample-seed "${SAMPLE_SEED}" \
    --output-dir "${genomes}"

  nvir="$(find "${genomes}/virus" -maxdepth 1 -name '*.fasta' 2>/dev/null | wc -l)"
  min_reads=$((nvir * READS_PER_ORG))
  echo "${label}: viral genomes=${nvir} x reads-per-organism=${READS_PER_ORG} => ~${min_reads} viral reads per chunk (plus 1 decoy read)."
  if [[ "${min_reads}" -lt 5000 ]]; then
    echo "WARNING: ${label}: n_virus * READS_PER_ORG = ${min_reads} < 5000. Increase max virus count or READS_PER_ORG."
  fi

  for L in "${LENGTHS[@]}"; do
    odir="${OUT}/${label}/length_${L}nt"
    mkdir -p "${odir}"
    echo "=== ${label} CHUNK length=${L}nt -> ${odir}/metagenome.fasta ==="
    MG chunk \
      --input "${genomes}" \
      --output metagenome.fasta \
      --output-dir "${odir}" \
      --sequence-length "${L}" \
      --reads-per-organism "${READS_PER_ORG}" \
      --abundance-profile "bacteria=0,virus=1,archaea=0,plasmid=0" \
      --seed "${SAMPLE_SEED}"
  done
}

mkdir -p "${OUT}"

run_one_panel "prokaryotic_viruses" "${SNAP_PROK}" "${MAX_VIRUS_PROK}" "${OUT}/prokaryotic_viruses/genomes"
run_one_panel "eukaryotic_viruses" "${SNAP_EUK}" "${MAX_VIRUS_EUK}" "${OUT}/eukaryotic_viruses/genomes"

echo "Done. Outputs under: ${OUT}/{prokaryotic_viruses,eukaryotic_viruses}/length_*nt/metagenome.fasta"
