# AGENTS Guide for SaVor

## Overview
SaVor is a Snakemake‑based pipeline for structural variant (SV) calling from short‑read data. Agents should understand the Snakemake workflow, configuration, and the conda environments that provide tool dependencies.

## Essential Commands
| Task | Command |
|------|---------|
| Create conda env (snakemake) | `mamba create -c conda-forge -c bioconda -n snakemake snakemake>=9.0` |
| Activate env | `mamba activate snakemake` |
| Dry‑run (local) | `snakemake -np --cores 1 --workflow-profile workflow-profiles/default --use-conda --rerun-incomplete` |
| Full run (local) | `snakemake -p --cores all --workflow-profile workflow-profiles/default` |
| Full run (SLURM cluster) | `snakemake --workflow-profile workflow-profiles/cluster-generic-slurm/ --cores <n>` |
| Resume after interruption | add `--rerun-incomplete` to the above commands |
| View config values | `snakemake --print-configfile` |

## Directory Layout & Architecture
- **Snakefile** – includes all rule files in `rules/`.
- **rules/** – individual Snakemake rule modules:
  - `fastq_ingestion.smk` – download/quality‑control reads.
  - `alignment.smk` – minibwa mapping, merging, duplicate marking.
  - `delly.smk`, `lumpy.smk`, `wham.smk` – SV callers.
  - `svcallprocess.smk` – VCF sorting, merging, filtering, consensus.
  - `benchmark.smk` – optional Truvari benchmark.
  - `user_bams.smk` – handling of user‑provided BAM/BAI paths.
- **config/** – `config.yaml` (pipeline options) plus CSVs for samples and contigs.
- **envs/** – conda environment YAML files for each tool (e.g., `fastq2bam.yml`, `lumpy.yaml`).
- **workflow‑profiles/** – `default/` (local) and `cluster-generic-slurm/` (SLURM) with `config.yaml` that sets thread counts.
- **docs/** – user guide, tutorial, troubleshooting.
- **results/** (generated) – hierarchical output:
  - `data/fastq/` – raw/filtered reads.
  - `bams/` – pre‑merge, post‑merge, final BAMs.
  - `SV/` – caller VCFs, merged VCFs, final consensus.
  - `benchmark/` – optional benchmark reports.

## Configuration (`config/config.yaml`)
- `samples`: path to sample sheet CSV (default `config/bam_samples.csv`).
- `include_contigs`: contig list for WHAM.
- `use_custom_reference`: `true` enables `refPath` column handling.
- `sort_reads`, `mark_duplicates`: toggle read sorting and duplicate marking (ignored for user‑provided BAMs).
- `remote_reads`: set to `true` + `remote_reads_prefix` for SRA reads stored remotely.
- `sv_merge`: 1/2/3 control caller consensus (default 3).
- `svBenchmark`: enable benchmarking; required truth‑set paths.

## Gotchas & Non‑Obvious Patterns
- **Snakemake version check** – `common.smk` forces version ≥8; pipeline expects ≥9 (see README).
- **Custom reference handling** – `has_custom_reference()` validates that all samples sharing a `refGenome` have the same non‑empty `refPath`. Inconsistent entries raise a `WorkflowError`.
- **User‑provided BAMs** – ruleorder `link_user_bam > dedup` ensures that when BAM paths are present, the pipeline links them instead of running deduplication. `mark_duplicates` flag is ignored in this case.
- **Remote reads** – if `remote_reads: true`, reads are wrapped with Snakemake `storage()` (Snakemake ≥8) or fallback `get_remote_reads`. Ensure `remote_reads_prefix` ends with a trailing slash.
- **Cluster profile** – `workflow-profiles/cluster-generic-slurm/config.yaml` sets per‑rule thread counts. The `status_control.sh` script reports SLURM job status; it expects the `--parsable` flag on `sbatch` calls (handled by Snakemake automatically).
- **Temporary directories** – many rules use `temp()` for intermediate files (e.g., BAMs before merging). Clean up occurs automatically on successful runs.
- **Benchmarking** – `benchmark.smk` only runs when `svBenchmark: true`. Truth‑set VCF paths must exist; otherwise the rule fails.
- **Dry‑run validation** – verify the sample sheet and config before heavy compute; missing columns (`fq1`, `fq2`, `bamPath`, etc.) will cause rule failures.
- **Conda environments** – each rule points to a specific environment file relative to the rule (`../envs/*.yml`). Do **not** modify these paths unless you relocate the `envs/` folder.

## Naming Conventions & Style
- Snakefile wildcards: `{refGenome}`, `{sample}`, `{run}` – always lower‑case, alphanumeric with underscores.
- Rule names are verb‑noun (`bwa_map`, `fastp`, `lumpy_call`). Follow the same pattern for new rules.
- Output files use the pattern `results/{refGenome}/...` to keep everything under a single root.
- Temporary files use `temp()` to ensure automatic cleanup.

## Testing Approach
- No unit tests are provided; validation is performed via Snakemake dry‑runs and by inspecting final VCFs.
- Use `--rerun-incomplete` after interruptions (see `docs/troubleshooting.md`).
- Verify configuration changes with `snakemake --print-configfile`.

## Frequently Used Scripts
- **`status_control.sh`** – SLURM job status helper used by the cluster profile.
- **Conda env files** – located in `envs/`; install with `mamba env create -f envs/<name>.yml` if needed to test isolated tool commands.

## Adding New Rules
1. Create `<name>.smk` in `rules/`.
2. Define rule(s) with `input`, `output`, `params`, `log`, `benchmark`, `conda` (point to appropriate env file).
3. Include the new file in `Snakefile` via `include: "rules/<name>.smk"`.
4. Update `workflow-profiles/*/config.yaml` with thread settings if needed.
5. Follow existing naming and temporary‑file conventions.

---
*Generated for agents interacting with the SaVor pipeline.*