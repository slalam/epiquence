# EPIquence usage guide

Everything needed to install, configure, run and maintain EPIquence. For an overview of the project see the [README](../README.md).

**Contents:**
[Requirements](#requirements) |
[Installation](#installation) |
[Containers](#containers-and-caches) |
[Reference databases](#reference-databases) |
[Configuration](#configuration) |
[Inputs](#inputs) |
[Running](#running) |
[Examples](#examples) |
[Parameters](#parameters) |
[Outputs](#outputs) |
[Species identification](#species-identification) |
[Surveillance](#surveillance) |
[Upgrading](#upgrading) |
[Troubleshooting](#troubleshooting)

---

## Requirements

EPIquence runs on a single Linux host with Nextflow's local executor. Choose one of two compute profiles. The profile sets the size of every task and selects the species engine and Bakta database.

| | `-profile low_capacity` | `-profile high_capacity` |
|---|---|---|
| Intended host | 8 to 32 threads, 16 to 64 GB RAM (for example 16 threads, 20 GB) | 64+ threads, 256+ GB RAM |
| Species identification | sylph, c1000 database (about 6 GB) | GTDB-Tk (110 GB per task, two at a time) |
| Bakta (core genome, plasmid maps) | light database, bundled in its container | full database |
| Reference databases on disk | about 60 GB | about 300 GB |

Without `-profile` the original settings apply (largest tasks request 50 CPUs / 60 GB).

**How sizing works.** Host size comes from `max_cpus` / `max_memory` in `conf/env.config`. Nextflow starts a task only when its reservation fits in what is free, and no task can request more than the host has. Per-task sizes match what each tool needs for a bacterial genome, so the machine fills up by running several tasks in parallel. A task killed for memory is retried up to twice with double the memory, capped at `max_memory`.

**Software:**

| Component | Version | Notes |
|---|---|---|
| Linux x86-64 | any recent distribution | container images are amd64 |
| Java | 17 or newer | for Nextflow |
| Nextflow | 25.10.2 | `manifest.nextflowVersion` |
| Apptainer or SingularityCE | Apptainer 1.2+ / SingularityCE 3.11+ | every analysis tool runs in a container |
| Conda | Miniforge recommended | a few Python steps and the database loaders |
| Python | 3.9+ on the host | a few small host-run steps |
| procps, tzdata | | `/usr/bin/ps` is used inside two containers, and tzdata handles time zones |
| marimo | 0.20+ | only for the web launcher |
| PostgreSQL | 14+ | Optional. SQLite is the default |

Network access is needed on first use to pull container images (Docker Hub, quay.io) and build Conda environments (conda-forge, bioconda), and once to download the reference databases. Analysis then runs offline. Plasmid comparison maps load two JavaScript libraries from `cdn.jsdelivr.net` when viewed.

---

## Installation

**Debian / Ubuntu**

```bash
sudo apt-get update
sudo apt-get install -y openjdk-17-jre-headless python3 procps tzdata curl git
# Apptainer:
# On Ubuntu use the PPA below. On Debian, install the .deb from https://github.com/apptainer/apptainer/releases
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt-get update && sudo apt-get install -y apptainer
# SingularityCE:
# Based on your Linux distribution, download the .deb from https://github.com/sylabs/singularity/releases, then
sudo apt-get install -y ./singularity-ce-4.6.0.deb
```

**RHEL / Rocky Linux / AlmaLinux**

```bash
sudo dnf install -y java-17-openjdk-headless python3 procps-ng tzdata curl git
# Apptainer:
sudo dnf install -y epel-release
sudo dnf install -y apptainer
# SingularityCE:
# Based on your Linux distribution, download the .rpm from https://github.com/sylabs/singularity/releases, then
sudo dnf install -y ./singularity-ce-4.6.0.rpm
```

Install either Apptainer or SingularityCE, not both. See [Containers and caches](#containers-and-caches).

**Both distributions** (Nextflow, Conda, folders and the pipeline):

```bash
# Nextflow 25.10.2
export NXF_VER=25.10.2
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
nextflow -version
apptainer --version || singularity --version

# Miniforge (Conda)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3 && /opt/miniforge3/bin/conda init bash
conda config --set channel_priority strict

# Folders used by the default configuration
sudo mkdir -p /data/epiquence /data/apps/singularity/{cache,layers} /data/apps/conda/envs /data/tmp
sudo chown -R "$USER" /data/epiquence /data/apps /data/tmp

# The pipeline
cd /data/epiquence
git clone https://github.com/slalam/epiquence.git pipeline
cd pipeline && chmod +x bin/*.py
```

**Web launcher (optional):**

```bash
python3 -m venv /opt/epiquence-launcher
/opt/epiquence-launcher/bin/pip install "marimo>=0.20"
/opt/epiquence-launcher/bin/marimo run /data/epiquence/pipeline/epiquence_launcher.py --host 0.0.0.0 --port 2718
```

The launcher calls `nextflow`, so Nextflow, Java and Apptainer must be on the `PATH` of the user running it.

Suggested layout:

```text
/data/
├── epiquence/      pipeline/ (this repository), reads/, assemblies/, results/, runs/, backup/
├── databases/      reference databases, results database, surveillance state
├── apps/           singularity/cache (images Nextflow runs), singularity/layers (runtime layer cache), Conda environments
└── tmp/            scratch, bound to /tmp inside containers
```

---

## Containers and caches

Every analysis tool runs in a container. Nextflow's `singularity` settings work with both **Apptainer** and **SingularityCE**, because Apptainer also installs a `singularity` command. Images are pulled from Docker Hub and quay.io and converted on first use.

Three locations matter. Put all of them on a large local disk:

| What | Set with | Template value | Default if unset |
|---|---|---|---|
| Images Nextflow runs (`.img` files, 20 to 40 GB) | `singularity.cacheDir` in `conf/env.config` (or `NXF_SINGULARITY_CACHEDIR`) | `/data/apps/singularity/cache` | `work/singularity` of each run |
| Layer cache used while pulling | `APPTAINER_CACHEDIR` or `SINGULARITY_CACHEDIR` in the shell | | `~/.apptainer/cache` or `~/.singularity/cache` (can fill your home directory) |
| Temporary build space while pulling | `APPTAINER_TMPDIR` or `SINGULARITY_TMPDIR` in the shell | | `/tmp` (often too small) |

The layer cache and build space are read by the runtime itself, not by Nextflow, so set them in the environment of every user who runs the pipeline or the launcher, for example in `/etc/profile.d/epiquence.sh`:

```bash
# Apptainer:
export APPTAINER_CACHEDIR=/data/apps/singularity/layers
export APPTAINER_TMPDIR=/data/tmp
# SingularityCE:
export SINGULARITY_CACHEDIR=/data/apps/singularity/layers
export SINGULARITY_TMPDIR=/data/tmp
```

In `conf/env.config`, the `singularity` block sets the image cache, a pull timeout of one hour (large images can take longer than Nextflow's default of 20 minutes) and the folders bound into containers:

```groovy
singularity {
    enabled     = true
    autoMounts  = true
    cacheDir    = "/data/apps/singularity/cache"
    pullTimeout = '1h'
    runOptions  = '--bind /data/databases:/data/databases --bind /data/tmp:/tmp'
}
```

Add a `--bind` for any other folder the containers must read (reads, assemblies or databases on another disk). On a shared cluster filesystem you can also point `singularity.libraryDir` at a read-only folder of pre-pulled images. Nextflow then uses those images without pulling.

To reclaim space from the layer cache later: `apptainer cache clean` or `singularity cache clean`.

### Prefetching containers (optional)

`scripts/prefetch_containers.sh` pulls every image the pipeline uses into the image cache once, before the first run, so runs start immediately and a slow download cannot hit the pull timeout mid-run. Images are named exactly as Nextflow names them, and images already present are skipped, so it is safe to re-run after an upgrade.

```bash
cd /data/epiquence/pipeline
scripts/prefetch_containers.sh --dry-run                  # list images and whether each is present
scripts/prefetch_containers.sh                            # all images (43)
scripts/prefetch_containers.sh --profile low_capacity     # skip GTDB-Tk and the full Bakta image (41)
scripts/prefetch_containers.sh --profile high_capacity    # skip the Bakta light image (42)
scripts/prefetch_containers.sh --runtime singularity --jobs 4
```

The image cache is read from `conf/env.config` (or `NXF_SINGULARITY_CACHEDIR`, or `--cache-dir`). The runtime is Apptainer if installed, otherwise SingularityCE (or `--runtime`). If you run the low capacity profile with `--taxonomy_tool gtdbtk` or `both`, prefetch with the default `--profile all`. Conda environments are not prefetched. Nextflow builds them on first use.

---

## Reference databases

Run downloads inside the pipeline's own containers so database and tool versions match. Commands are shown for Apptainer and for SingularityCE, so use the one you installed.

| Database | Setting in `conf/env.config` | Needed for | Size |
|---|---|---|---:|
| BUSCO `bacteria_odb12` | `busco_db` | always | under 1 GB |
| sylph GTDB R226 + GTDB taxonomy files | `sylph_db`, `sylph_taxonomy` | low capacity | 3.7 GB (c1000) |
| GTDB-Tk R226 package | `gtdbtk_db` | high capacity | about 110 GB |
| Bakta full database | `bakta_db` | high capacity with `--core_genome` or `--hgt_maps` | about 70 GB |
| Results database | `surveillance_db_url` | `--surveillance` | grows over time |

All other databases (AMR, MLST, typing, plasmids, virulence, the Bakta light database) are inside the container images.

**BUSCO (always):**

```bash
mkdir -p /data/databases/busco
# Apptainer:
apptainer exec -B /data/databases docker://staphb/busco:6.0.0 \
  busco --download bacteria_odb12 --download_path /data/databases/busco
# SingularityCE:
singularity exec -B /data/databases docker://staphb/busco:6.0.0 \
  busco --download bacteria_odb12 --download_path /data/databases/busco
```

**sylph and GTDB taxonomy (low capacity).** Use the same GTDB release (R226) as GTDB-Tk, so species names match between engines. Use the c1000 database on small hosts. The c200 database (18.4 GB, about 24 GB RAM) is better suited to servers.

```bash
mkdir -p /data/databases/sylph && cd /data/databases/sylph
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c1000-dbv1.syldb
# Optional: c200 database (18.4 GB on disk, needs about 24 GB RAM), then point sylph_db at it
# wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_taxonomy_r226.tsv.gz
```

**GTDB-Tk R226 (high capacity):**

```bash
mkdir -p /data/databases/gtdb && cd /data/databases/gtdb
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
tar -xzf gtdbtk_r226_data.tar.gz && rm gtdbtk_r226_data.tar.gz        # creates release226/
# Check the installation
# Apptainer:
apptainer exec -B /data/databases --env GTDBTK_DATA_PATH=/data/databases/gtdb/release226 \
  docker://ecogenomic/gtdbtk:2.6.1 gtdbtk check_install
# SingularityCE:
singularity exec -B /data/databases --env GTDBTK_DATA_PATH=/data/databases/gtdb/release226 \
  docker://ecogenomic/gtdbtk:2.6.1 gtdbtk check_install
```

**Bakta full database (high capacity):**

```bash
mkdir -p /data/databases/bakta/1.11.4
# Apptainer:
apptainer exec -B /data/databases docker://staphb/bakta:1.11.4 \
  bakta_db download --output /data/databases/bakta/1.11.4 --type full
# SingularityCE:
singularity exec -B /data/databases docker://staphb/bakta:1.11.4 \
  bakta_db download --output /data/databases/bakta/1.11.4 --type full
```

If a download URL has moved, check the GTDB release folder for R226 and the sylph pre-built databases page.

**Results database (used with `--surveillance`).** Tables are created and extended automatically, so there are no migrations to run.

- SQLite (default, one host): `surveillance_db_url = "sqlite:////data/databases/epiquence/epiquence_db.db"` (four slashes for an absolute path). One writer at a time, so avoid overlapping runs.
- PostgreSQL (recommended for production, concurrent runs or Metabase):

    ```sql
    CREATE ROLE epiquence LOGIN PASSWORD 'change-me';
    CREATE DATABASE epiquence OWNER epiquence;
    ```

    `surveillance_db_url = "postgresql+psycopg2://epiquence:change-me@db.example.org:5432/epiquence"`

**Surveillance state.** Create the folders below. Back up the sketch folder and the results database together: the sketch folder keeps cluster names stable between runs.

```bash
mkdir -p /data/databases/epiquence/{surveillance,surveillance-results,analysis_reports}
```

---

## Configuration

Host-specific settings live in `conf/env.config`, which is gitignored. The pipeline will not start without it.

```bash
cp conf/env.config.template conf/env.config
chmod 600 conf/env.config
```

| Setting | Example | Purpose |
|---|---|---|
| `max_cpus`, `max_memory` | `16`, `18.GB` | Host size: all threads, and RAM minus about 2 GB (small host) or 8 GB (server). Drives the executor and the per-task ceiling. |
| `busco_db`, `gtdbtk_db`, `bakta_db` | paths | Reference databases |
| `sylph_db`, `sylph_taxonomy` | `.syldb` path, and comma-separated taxonomy files | sylph engine |
| `taxonomy_tool` | `sylph` | Command-line runs only: `gtdbtk`, `sylph` or `both`. The launcher always follows the profile |
| `surveillance_db_url` | `sqlite:////...` or `postgresql+psycopg2://user:pass@host/db` | Results database |
| `surveillance_sketch_dir` | `/data/databases/epiquence/surveillance` | Cumulative surveillance state |
| `surveillance_results_dir` | `/data/databases/epiquence/surveillance-results` | Surveillance trees and distance matrices |
| `analysis_report_dir` | `/data/databases/epiquence/analysis_reports` | Archive of every run's report |
| `report_timezone` | `America/Chicago` | IANA zone for all timestamps |
| `metabase_base_url`, `metabase_db_id`, `metabase_api_key` | | Optional Metabase refresh after each load |
| `ai_summary`, `ai_summary_full`, `ai_model`, `ai_endpoint`, `ai_context_window` | `true`, `true`, `'gemma4'`, `'http://localhost:11434/v1/chat/completions'`, `131072` | Optional AI summary in the report through an OpenAI-compatible endpoint (for example Ollama). `ai_context_window` must match the model's context length |

The `executor` block in the template reads the host size: `executor { name = 'local'; cpus = params.max_cpus; memory = params.max_memory }`. Also set the container image cache, pull timeout and bind paths (see [Containers and caches](#containers-and-caches)) and the Conda cache to your folders.

**Credentials.** Keep the database password and Metabase key in `conf/env.config` with `chmod 600`. At start-up the pipeline writes them to a private per-run file (mode 600) that only the database and Metabase steps read, and deletes it when the run ends. The values are never written to `work/`, the reports or the logs.

**Launcher settings** (`epiquence_config.json`, next to the launcher):

| Key | Meaning |
|---|---|
| `pipeline_dir` | Folder containing `main.nf`. `"."` means the launcher's own folder |
| `output_base` | Each launch creates `<output_base>/W<MMDDYYHHMMSS>/` |
| `reads_base` | Its sub-folders are offered as read and assembly folders |
| `reads_sort`, `reads_sort_desc`, `reads_show_date` | Folder ordering and labels |
| `groups` | Group IDs and descriptions. Groups added in the launcher are saved here |
| `default_profile` | `"high_capacity"` or `"low_capacity"`, preselected in the launcher |

---

## Inputs

A run needs a **sample sheet** of reads (`--sample_sheet`), a **folder of assemblies** (`--assemblies_dir`), or both.

**Sample sheet:** a CSV with exactly these columns. Leave a cell empty when a read type is not available. Use absolute paths and gzipped FASTQ.

```csv
sample_id,long_R,short_R1,short_R2
SAMPLE-01,,/data/epiquence/reads/batch1/SAMPLE-01_R1.fastq.gz,/data/epiquence/reads/batch1/SAMPLE-01_R2.fastq.gz
SAMPLE-02,/data/epiquence/reads/batch1/SAMPLE-02_L.fastq.gz,,
SAMPLE-03,/data/epiquence/reads/batch1/SAMPLE-03_L.fastq.gz,/data/epiquence/reads/batch1/SAMPLE-03_R1.fastq.gz,/data/epiquence/reads/batch1/SAMPLE-03_R2.fastq.gz
```

Each assembly type only runs for rows with the reads it needs.

**Assembly types** (`--assembly_type`, comma-separated):

| Code | Meaning | Reads needed |
|---|---|---|
| `S` | short-read assembly | Illumina R1 + R2 |
| `L` | long-read assembly with polishing | Nanopore |
| `HL` | hybrid, long-read first, short-read polished | both |
| `HS` | hybrid, short-read first | both |
| `ALL` | all four | both |

Every sample is analyzed once per assembly type and named `<sample_id>_<type>` (for example `SAMPLE-01_HL`). Imported assemblies are named `<name>_A`.

**Sample ID rules:** use letters, digits and hyphens (for example `SAMPLE-01`), and keep each ID unique within the sheet. Do not use spaces, slashes or dots, and avoid IDs ending in `_S`, `_L`, `_HL`, `_HS` or `_A`. Never use the suffix `_O<n>`, which marks superseded surveillance versions. With `--surveillance`, re-using an ID in a later run means the isolate was re-sequenced: the earlier version is archived and replaced.

**Validation:** before analysis, every read file is checked for gzip integrity, corrupted bytes and valid FASTQ. Any failure stops the run and lists the files. Nanopore reads shorter than 1,000 bp are discarded.

**Assemblies folder:** every `*.fasta`, `*.fna` or `*.fa` file becomes a sample (`<file name>_A`). Each file needs at least one contig and 100 bp. Imported assemblies skip read QC, assembly and coverage.

**Metadata (optional, `--metadata`):** a CSV (or a TSV with a `.tsv` name) with a `sample_id` column. Other columns (for example `collection_date`, `facility`, `specimen_source`) are shown on outbreak trees and tables. Rows match the full sample name or the bare isolate ID. In the launcher, put the file in the Reads Directory and select it from the **Metadata File** dropdown.

**Launcher file names:** the launcher recognizes `<id>_R1.fastq.gz`, `<id>_R2.fastq.gz` and `<id>_L.fastq.gz` (also `.fq.gz`).

---

## Running

**Command line:**

```bash
nextflow run /data/epiquence/pipeline/main.nf -profile low_capacity \
    --sample_sheet samples.csv --assembly_type S \
    --outdir results \
    --workflow_id "W$(date +%m%d%y%H%M%S)" --workflow_description "Weekly batch" \
    --group_id GP001 --group_description "Healthcare-Associated Infections" \
    --outbreak --plasmid_hgt --surveillance
```

Always pass `--workflow_id`: it labels the report, every table and every database row.

| Nextflow option | Effect |
|---|---|
| `-profile low_capacity` / `high_capacity` | Task sizing, species engine and Bakta database |
| `-resume` | Re-use finished steps from the cache |
| `-bg` | Run in the background (output in `.nextflow.log`) |
| `-w <dir>` | Work directory (default `./work`) |
| `-params-file params.yaml` | Read parameters from a file |
| `-preview` | Check configuration and databases without running tasks |

The start-up banner shows the version, species engine and profile. Invalid settings or missing databases stop the run at this point with a clear message.

**Web launcher:** start it (see [Installation](#installation)) and open it in a browser.

1. Enter a *Workflow Description* and choose a *Group* (or add one).
2. Tick the assembly types and the long-read subsampling option.
3. Choose a *Reads Directory* (and optionally an *Assemblies Directory*) and click **Scan**. You can deselect or rename samples in the table.
4. Choose the **Compute Profile**: `High capacity | 64+ threads, 256+ GB RAM` or `Low capacity | 8 to 32 threads, 16 to 64 GB RAM`. The host size from `conf/env.config` is shown underneath.
5. Set outbreak detection, core genome, plasmid HGT and comparison maps, and surveillance. For metadata, choose a **Metadata File** from the dropdown, which lists the CSV and TSV files in the selected Reads Directory (place your metadata file there).
6. Click **Launch Workflow**.

Before launching, the launcher checks the setup: a missing `gtdbtk_db` (high capacity), `sylph_db` or taxonomy file (low capacity), or `conf/env.config` blocks the launch and shows the reason. A missing `bakta_db`, an undersized host or missing `max_cpus` / `max_memory` show as warnings. Each launch creates `<output_base>/W<MMDDYYHHMMSS>/` with `raw_reads/` (symlinks), `sample_sheet.csv`, `run.sh` (the exact command), `nextflow_launch.log`, `work/` and `results/`.

The launcher turns on outbreak detection, plasmid HGT, comparison maps and surveillance by default and skips long-read subsampling. On the command line, all modules are off unless you add the flags.

**Resume.** Re-run the same command with `-resume`, for example after a failure or to add a module. In the launcher, *Resume Previous Workflow* lists earlier runs with their status, presets the profile they used, applies the module toggles you choose and resumes in the background. Changing the profile of an existing run re-runs most steps. Report and database steps always re-run.

**Partial runs:** `--stop_at qc`, `assembly` or `polish` stops after that stage, for example to inspect assemblies. Continue later with the full command and `-resume`.

**Monitoring:** follow progress with `tail -f .nextflow.log`. Per-task resources are in `results/pipeline_info/trace.txt`, and timeline and resource reports are in `results/pipeline_info/`. For a failed task, look in its work directory (shown in the error) at `.command.sh`, `.command.err` and `.command.log`.

---

## Examples

```bash
EPQ=/data/epiquence/pipeline/main.nf
WID="W$(date +%m%d%y%H%M%S)"
```

**Nanopore isolates, with subsampling to 100x**

```bash
nextflow run $EPQ -profile low_capacity --sample_sheet ont.csv --assembly_type L \
    --rasusa_coverage 100 --rasusa_genome_size 5.5mb \
    --outdir results --workflow_id $WID --workflow_description "ONT batch"
```

**Compare all assembly strategies**

```bash
nextflow run $EPQ -profile high_capacity --sample_sheet hybrid.csv --assembly_type ALL \
    --outdir results --workflow_id $WID --workflow_description "Assembler comparison"
```

**Outbreak investigation with core-genome confirmation**

Core-genome confirmation requires at least 5 samples in an outbreak cluster.

```bash
nextflow run $EPQ -profile high_capacity --sample_sheet batch.csv --assembly_type S \
    --outbreak --outbreak_threshold 0.05 --core_genome --metadata metadata.csv \
    --outdir results --workflow_id $WID --workflow_description "Ward 4 investigation"
```

**Resistance plasmid spread**

```bash
nextflow run $EPQ -profile high_capacity --sample_sheet batch.csv --assembly_type HL \
    --plasmid_hgt --hgt_maps --amr_priority_genes 'blaKPC*,blaNDM*,blaOXA-48*' \
    --outdir results --workflow_id $WID --workflow_description "Carbapenemase plasmids"
```

**Validate sylph against GTDB-Tk on known isolates**

This only works with the `high_capacity` profile, because it runs GTDB-Tk.

```bash
nextflow run $EPQ -profile high_capacity --taxonomy_tool both \
    --assemblies_dir /data/epiquence/assemblies/known_isolates \
    --outdir results --workflow_id $WID --workflow_description "sylph validation"
```

---

## Parameters

Set with `--name value` on the command line or in a `-params-file`. Defaults are in `nextflow.config`.

| Parameter | Default | Description |
|---|---|---|
| `--sample_sheet` | none | Reads sample sheet (CSV) |
| `--assemblies_dir` | none | Folder of assemblies |
| `--outdir` | `results` | Output folder |
| `--metadata` | none | Metadata CSV/TSV for outbreak trees and tables |
| `--workflow_id` | none | Run identifier. Always set it |
| `--workflow_description` | empty | Free text for the report and database |
| `--group_id`, `--group_description` | `GP001`, `Healthcare-Associated Infections` | Program group |
| `--assembly_type` | `ALL` | `L`, `S`, `HL`, `HS` (comma-separated) or `ALL` |
| `--medaka_model` | `r1041_e82_400bps_sup_v5.0.0` | Long-read polishing model. Match it to your chemistry |
| `--skip_rasusa` | `false` | Use all long reads (no subsampling) |
| `--rasusa_coverage`, `--rasusa_genome_size` | `100`, `6mb` | Long-read subsampling target |
| `--stop_at` | `all` | `qc`, `assembly`, `polish` or `all` |
| `--taxonomy_tool` | from `conf/env.config` or the profile | `gtdbtk`, `sylph` or `both` |
| `--sylph_species_ani`, `--sylph_borderline_ani`, `--sylph_min_margin` | `95`, `94`, `0.5` | sylph call thresholds |
| `--outbreak` | `false` | Clonal outbreak detection (2+ samples) |
| `--outbreak_threshold` | `0.05` | Clustering distance. A lower value gives tighter clusters |
| `--core_genome` | `false` | Core-genome confirmation for clusters of 5+ |
| `--plasmid_hgt` | `false` | Plasmid HGT detection (2+ samples) |
| `--hgt_maps` | `false` | Plasmid comparison maps for AMR-bearing plasmid groups |
| `--hgt_mash_threshold` | `0.005` | Plasmid clustering distance |
| `--amr_priority_genes` | empty | Genes that rank the reference plasmid, for example `'blaKPC*,blaNDM*'` |
| `--comparison_map_palette` | `cool` | `warm`, `cool` or `earth` |
| `--surveillance` | `false` | Store results and add isolates to cumulative surveillance |
| `--surveillance_group_views` | `true` | Also write per-group surveillance pages |
| `--surveillance_min_chromosome_bp` | `500000` | Skip isolates whose chromosome is smaller (likely incomplete) |
| `--multiqc_title`, `--multiqc_subtitle`, `--multiqc_intro` | | Report text |

Species routing (which organisms get point-mutation calling, which MLST scheme is forced, and which species go to each typing tool) is set by the maps and `*_species` lists in `nextflow.config`.

---

## Outputs

```text
results/
├── multiqc/multiqc_report.html     the interactive report (start here)
├── summaries/                      combined tables and sequences
├── <sample>/                       per-sample QC, assemblies, stats, annotation, typing
├── outbreak/                       with --outbreak
├── plasmid_hgt/<assembly type>/    with --plasmid_hgt
└── pipeline_info/                  Nextflow report, timeline, trace
```

| What | Where |
|---|---|
| Interactive report | `results/multiqc/multiqc_report.html` (also archived as `<analysis_report_dir>/<workflow_id>_report.html`) |
| Final assemblies | `results/<sample>/assemblies/`, `results/summaries/assemblies/` |
| Chromosome and plasmid sequences | `results/summaries/chromosomes/`, `results/summaries/plasmids/<type>/` |
| Combined tables | `results/summaries/*.tsv`: general statistics, taxonomy, MLST, AMR genes, point mutations, resistance phenotypes, plasmid and chromosomal resistance markers, virulence, typing, software versions |
| Outbreak clusters | `results/outbreak/`: cluster list, SNP matrices, trees, pairwise distances |
| Plasmid HGT | `results/plasmid_hgt/<type>/`: plasmid clusters, pairwise ANI, AMR plasmid tree and matrix, UMAP plot, comparison maps |
| Run diagnostics | `results/pipeline_info/` |
| Surveillance pages | `<surveillance_results_dir>/all/`, `/<group_id>/`, `/ungrouped/` |
| Database tables | `surveillance_db_url` |

---

## Species identification

The species call drives the MLST scheme, point-mutation calling, species-specific typing, outbreak grouping and the surveillance population an isolate joins. Two engines produce identical species names (both use GTDB R226):

| | GTDB-Tk | sylph |
|---|---|---|
| Method | marker genes and phylogenetic placement | containment ANI of the assembly against GTDB species representatives |
| Memory | about 90+ GB | about the database size (c1000: about 6 GB) |
| Default for | `high_capacity`, no profile | `low_capacity` |

**Which engine runs:** in the launcher, the profile decides. On the command line: `--taxonomy_tool`, then `taxonomy_tool` in `conf/env.config`, then the profile default, then GTDB-Tk.

**sylph call rule** (top hit containment ANI): 95 or more gives the species, and 94 to 95 gives the species flagged `borderline`. Below 94, only the genus is reported. If another species is within 0.5 ANI, the call is flagged `ambiguous`, and with no hit it is `unclassified`. Containment ANI on a single isolate runs slightly below true ANI, which is why the borderline band exists.

**Validation mode** (`--taxonomy_tool both`): GTDB-Tk stays authoritative and sylph is reported beside it with an Agreement column. Run known isolates this way before relying on sylph alone, particularly for close species complexes.

GTDB names some traditional subspecies as species (for example *Legionella pneumophila* subsp. *fraseri* is *Legionella fraseri*). Both engines report GTDB names.

---

## Surveillance

With `--surveillance`, results are stored in the database and each isolate's chromosome joins its species' cumulative population in `surveillance_sketch_dir`. Every run recomputes distances for the whole species and clusters at 10, 20, 50 and 100 SNPs, keeping cluster names stable and recording changes.

Pages are published per species in `surveillance_results_dir`:

```text
surveillance-results/
├── all/          species-trees/<Species>_surveillance_tree.html, distance-matrix/<Species>_distance_matrix.{html,tsv}
├── GP001/        the same, for isolates of group GP001
└── ungrouped/    isolates with no group recorded
```

Group pages use the species-wide distances and cluster names, so values match `all/`. The tree and the matrix order are recomputed for the group. An isolate belongs to the group of the latest run that processed it. A group with fewer than 2 isolates of a species gets no page. Turn group pages off with `--surveillance_group_views false`.

Back up `surveillance_sketch_dir` and the results database together, and restore them together.

---

## Upgrading

1. Let running workflows finish. Back up `conf/env.config`, `epiquence_config.json`, the surveillance sketch folder and the results database.
2. Update the code on a branch (`git pull`, or copy a release over the repository without overwriting `conf/env.config` and `epiquence_config.json`), then `chmod +x bin/*.py`.
3. Compare `conf/env.config` with `conf/env.config.template` and add new settings (for this release: `max_cpus`, `max_memory`, the `executor` block reading them, the sylph settings for low-capacity hosts, and optionally the AI settings).
4. Add `"default_profile"` to `epiquence_config.json`. Optionally prefetch new container images: `scripts/prefetch_containers.sh` (only missing images are pulled).
5. Verify, then run a few known isolates:

    ```bash
    nextflow config -profile low_capacity . | grep resourceLimits
    nextflow run main.nf -preview -profile low_capacity --assemblies_dir /path/to/test/assemblies --outdir /tmp/p --workflow_id W_PREVIEW
    ```

6. Restart the launcher.

To roll back, check out the previous version and restore the backed-up files. If a new run used `--surveillance`, restore the sketch folder and the database together.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Pipeline does not start: `conf/env.config` not found | Not created | `cp conf/env.config.template conf/env.config` |
| `Either --sample_sheet or --assemblies_dir ... must be provided` | Wrong path | Use absolute paths |
| No QC runs, all reads treated as missing | Sample sheet header misspelt | Header must be `sample_id,long_R,short_R1,short_R2` |
| `VALIDATION FAILED` | Corrupt or truncated read file | Re-transfer the listed files |
| `taxonomy_tool 'sylph' needs sylph_db` or similar | Database missing or path wrong | See [Reference databases](#reference-databases) |
| `Process requirement exceeds available CPUs/memory` | Host size not set, no profile | Set `max_cpus` / `max_memory` and use a profile |
| GTDB-Tk killed (exit 137) | Out of memory | Use `high_capacity` on a 256+ GB host, or sylph |
| `command not found` (exit 127) | Container not used: Apptainer not enabled, or an edited profile replaced a container setting | Check `singularity.enabled` in `conf/env.config`. Profile selectors must keep the `'^(NAME)$'` form |
| Container cannot find a database path | Folder not bound | Add `--bind /your/path:/your/path` to `singularity.runOptions` |
| `Failed to pull singularity image ... Try and increase singularity.pullTimeout` | Slow download of a large image | Keep `pullTimeout = '1h'` (or longer), or run `scripts/prefetch_containers.sh` first |
| Home directory full, or `no space left on device` while pulling | Layer cache or build space on a small disk | Set `APPTAINER_CACHEDIR` / `SINGULARITY_CACHEDIR` and `APPTAINER_TMPDIR` / `SINGULARITY_TMPDIR` to a large disk |
| Every run pulls images again | No image cache set, so each run uses its own `work/singularity` | Set `singularity.cacheDir` in `conf/env.config` |
| `database is locked` | SQLite with overlapping runs | Do not overlap runs, or use PostgreSQL |
| `permission denied for schema public` | PostgreSQL 15+ privileges | Make the pipeline role the database owner |
| Unexpected species, `borderline` or `ambiguous` | Close species or isolate near the species boundary | Review the call, and compare engines with `--taxonomy_tool both` |
| No point mutations for a species | Species not in the AMR organism maps | Add it to the maps in `nextflow.config` |
| Isolate missing from surveillance | Chromosome below `--surveillance_min_chromosome_bp` | Check the task log, and lower the threshold only for small genomes |
| Warning that `conf/env.config` is readable by others | Permissions | `chmod 600 conf/env.config` |
| Resume re-runs everything | Different launch folder or changed parameters | Resume from the same folder with the same command |
