# Viral-only-style metagenomes (5 lengths × prok / euk slices)

CHIMERA requires one non-viral FASTA; this script adds a synthetic `bacteria/` decoy and uses `--abundance-profile bacteria=0,virus=1` so each chunk run emits ~`n_virus × reads-per-organism` viral reads plus **one** decoy read.

## Run

```bash
cd /path/to/MetagenomeGenerator
export ENTREZ_EMAIL="you@example.com"
bash experiments/viral_only_dataset/run_five_lengths.sh
```

Outputs: `output/viral_only_dataset/{prokaryotic_viruses,eukaryotic_viruses}/length_{150,250,500,1000,3000}nt/metagenome.fasta`

## Scale (≥ 5000 reads)

Defaults: `READS_PER_ORG=3`, `MAX_VIRUS_PROK=2500`, `MAX_VIRUS_EUK=1876` → ~7500 / ~5628 viral reads per metagenome.

Override:

```bash
MAX_VIRUS_PROK=4000 READS_PER_ORG=2 OUT=/tmp/vod bash experiments/viral_only_dataset/run_five_lengths.sh
```

Requires prior stratified snapshots: `experiments/viral_host_stratify/out/accession_snapshot_viral_{prok,euk}.json`.
