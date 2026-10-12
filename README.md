# Frailty WGBS: BSBolt Index and Alignment on Hoffman2

Workflow for building a bisulfite index and aligning the Human_Blood_Frailty samples (96 trimmed paired-end FASTQ pairs) with BSBolt 1.6.0 on the UCLA Hoffman2 cluster.

## Overview

1. Set up a Python environment with BSBolt (once).
2. Build the hg38 BSBolt index (once, as a queued job).
3. Make a sample list.
4. Align each sample (SGE array job), then sort and index the BAMs.
5. Check results, then scale from 1 test sample to all 96.

| Item | Value |
|---|---|
| Cluster | Hoffman2 (SGE scheduler) |
| Raw data (read-only) | `/u/project/pellegrini_archive/data/TBS2/fastq/Human_Blood_Frailty` |
| Reference | `/u/project/pellegrini_archive/data/Charlotte/run_program/run_program/refdata-gex-GRCh38-2024-A/fasta/genome.fa` |
| Work directory | `$SCRATCH/frailty_project` (`/u/scratch/j/jzeng`) |
| Tools | BSBolt 1.6.0 (in a venv), samtools 1.15 (module) |

## Key concepts

- **Login node (head node):** only for light work such as `ls`, editing, and submitting jobs. Never compute here.
- **`qsub`:** submits a script to the queue. The job keeps running even if your SSH session drops.
- **`qrsh`:** gives an interactive shell on a compute node, useful for testing a command.
- **Archive:** `/u/project/pellegrini_archive/...` is shared storage for raw data.
- **Scratch:** `$SCRATCH` is where large outputs (indexes, BAMs) go.
- **Job states in `qstat`:** `qw` waiting, `hqw` held until another job finishes, `r` running, nothing listed means finished or failed.

## Step 0: One-time setup

```bash
ssh jzeng@hoffman2.idre.ucla.edu

module load python/3.9.6
python3 -m venv ~/bsbolt_env
source ~/bsbolt_env/bin/activate
pip install --upgrade pip
pip install BSBolt
bsbolt --help

mkdir -p $SCRATCH/frailty_project/{index,logs,alignments}
```

Use `$SCRATCH`, not `/u/scratch/jzeng`. The real path is `/u/scratch/j/jzeng`, and `mkdir` on the wrong one fails with permission denied.

## Step 1: Index generation

The `hg38.fa.gz` in Mike's archive is only a text list of chromosome sizes, not a genome. Use the 10x GRCh38 FASTA instead.

```bash
cat > $SCRATCH/frailty_project/build_index.sh <<'EOF'
#!/bin/bash
#$ -cwd
#$ -j y

. /u/local/Modules/default/init/bash
module load python/3.9.6
source ~/bsbolt_env/bin/activate

REF="/u/project/pellegrini_archive/data/Charlotte/run_program/run_program/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
DB="$SCRATCH/frailty_project/index/hg38_bsbolt"

mkdir -p "$(dirname "$DB")"

bsbolt Index -G "$REF" -DB "$DB"
EOF

qsub -l h_data=16G,h_rt=24:00:00 -pe shared 4 \
  -o $SCRATCH/frailty_project/logs/index.log \
  $SCRATCH/frailty_project/build_index.sh
```

`qsub` prints the job ID. Note it for Step 3. The index took several hours in this run. Watch progress with `tail -2 $SCRATCH/frailty_project/logs/index.log`.

## Step 2: Sample list

```bash
ls /u/project/pellegrini_archive/data/TBS2/fastq/Human_Blood_Frailty/*_R1_trimmed.fastq.gz \
  | xargs -n1 basename | sed 's/_R1_trimmed.fastq.gz//' > $SCRATCH/frailty_project/samples.txt
wc -l $SCRATCH/frailty_project/samples.txt   # expect 96
```

## Step 3: Alignment

`#$ -t 1-1` runs only the first sample as a test.

```bash
cat > $SCRATCH/frailty_project/align.sh <<'EOF'
#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
#$ -t 1-1
#$ -o /u/scratch/j/jzeng/frailty_project/logs/

. /u/local/Modules/default/init/bash
module load python/3.9.6
source ~/bsbolt_env/bin/activate
module load samtools/1.15

FQ=/u/project/pellegrini_archive/data/TBS2/fastq/Human_Blood_Frailty
OUT=$SCRATCH/frailty_project/alignments
DB=$SCRATCH/frailty_project/index/hg38_bsbolt
S=$(sed -n "${SGE_TASK_ID}p" $SCRATCH/frailty_project/samples.txt)

bsbolt Align -F1 $FQ/${S}_R1_trimmed.fastq.gz -F2 $FQ/${S}_R2_trimmed.fastq.gz \
  -DB $DB -O $OUT/$S -t $NSLOTS

samtools sort -@ $NSLOTS -o $OUT/$S.sorted.bam $OUT/$S.bam
samtools index $OUT/$S.sorted.bam
EOF
```

Submit it held until the index job finishes. Replace `INDEX_JOB_ID` with the real number (no `<` or `>`):

```bash
qsub -hold_jid INDEX_JOB_ID $SCRATCH/frailty_project/align.sh
```

## Step 4: Check results

```bash
qstat -u $USER
tail -20 $SCRATCH/frailty_project/logs/index.log
ls -lh $SCRATCH/frailty_project/index/hg38_bsbolt/
ls -lh $SCRATCH/frailty_project/alignments/
tail -20 $SCRATCH/frailty_project/logs/align.sh.o*

module load samtools/1.15
samtools flagstat $SCRATCH/frailty_project/alignments/130010.sorted.bam
```

In this run the test sample (130010) produced `130010.bam` (812M), `130010.sorted.bam` (628M), and `130010.sorted.bam.bai` (5.2M), and the index folder contained per-contig `.pkl` files. Check the mapped percentage from `flagstat` before scaling up (not yet verified).

## Step 5: Run all 96 samples

```bash
sed -i 's/#\$ -t 1-1/#$ -t 1-96/' $SCRATCH/frailty_project/align.sh
qsub $SCRATCH/frailty_project/align.sh
```

## Troubleshooting log

| Symptom | Cause | Fix |
|---|---|---|
| `libpython3.9.so.1.0: cannot open shared object file` | Compute node did not load the Python module the venv depends on | Load Python in the script before activating the venv |
| `module: command not found` | Batch jobs do not have the `module` command | Add `. /u/local/Modules/default/init/bash` at the top of every script |
| `argument -t: expected one argument` | `$NSLOTS` is empty outside a queued job | Ignore on the login node; it is set inside jobs |
| Alignment job shows `qw` instead of `hqw` | `-hold_jid` was given a placeholder or wrong job ID | `qdel` it and resubmit with the real index job ID |
| `bash: new_index_job_id: No such file` | `<...>` placeholder was read as file redirection | Use the actual number |
| `qstat` shows nothing | Job finished or failed | Check the log and output folders to tell which |

## Notes

- Every script that uses Python or samtools needs the module init line first.
- Do heavy work only through `qsub` or `qrsh`, never on the login node.
- Keep scripts (not data) in GitHub: `build_index.sh`, `align.sh`, the sample-list command, and this README.
- Methylation calling (`bsbolt CallMethylation`) is not covered here.
