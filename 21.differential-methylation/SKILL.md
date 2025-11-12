---
name: differential-methylation
description: This skill performs differential DNA methylation analysis (DMRs and DMCs) between experimental conditions using WGBS methylation tracks (BED/BedGraph). It standardizes input files into per-sample five-column Metilene tables, constructs a merged methylation matrix, runs Metilene for DMR detection, filters the results, and generates quick visualizations. The user must specify the column indices corresponding to methylation_fraction and coverage values.
---

# WGBS Differential Methylation with metilene

## Overview

- **Validate inputs**: Ask the user to specify the column indices corresponding to methylation_fraction and coverage values. **Never dicide by yourself**.
- **Standardize**: convert heterogeneous inputs to a **per‑sample 5‑column Metilene table** (chrom, start, end, methylation_fraction, coverage). 
- **Define groups**: create group files mapping samples to group A/B for metilene. Make sure that each group has at least 2 replicates. 
- **Merge**: build a **merged methylation matrix** (coordinates + one methylation column per sample) for metilene using `bedtools unionbedg`. Remember to add headers ("chorm\tstart\end\t\g1_rep1\g2_rep2...") to the merged file.
- **Run metilene**: call DMRs with tunable parameters. The `-a` and `-b` parameter should be the group name (like cell type name or control/treated and so on)
- **Filter & export**: retain significant hits and emit BED/TSV summaries plus QC stats.  
- **Optional DMC mode**: treat single CpGs as regions (min‑CpG=1, max‑dist=1).  
- **Visualize**: quick plots (Δmethylation vs –log10(q), length histograms) and browser tracks.

---

## Inputs

- **Per-sample methylation tracks** at CpG resolution (one file per sample), supported raw formats:
  - **BedGraph** with `%methylation`, **methylated_count**
  - **Bismark coverage/bedGraph-like** with `%methylation`, **methylated_count**, **unmethylated_count**.
  - **BED-like** with `%methylation`, **methylated_count**


- **Groups file** (`groups.tsv`): whitespace‑separated lines: `<sample_id>  <group_name>`  
  Example:
  ```text
  K562_rep1   case
  K562_rep2   case
  Ery_rep1    control
  Ery_rep2    control
  ```

> **Assumptions**: All samples share the same reference genome build and chromosome naming scheme.

---

## Outputs

- `results/` directory containing:
  - `dmr_results.txt` – raw metilene output.
  - `significant_dmrs.txt` – filtered significant DMRs (TSV).
  - `significant_dmrs.bed` – BED for genome browser.
  - `dmr_summary.txt` – counts and length statistics.
  - `plots/` – `volcano.png`, `length_hist.png`.
  - Optional DMC outputs mirroring the above.
- `work/` directory (intermediate):
  - `per_sample/*.metilene.tsv` – **5‑column per‑sample Metilene tables**.
  - `merged_methylation.bed` – merged methylation matrix (coords + 1 column per sample).

---

## Workflow Decision Tree

### 1) Validate & detect input modes

```bash
mkdir -p work/per_sample results plots
```

**Heuristics** (used in conversion below):
- Values in `[0,1]` → **fraction**; values in `[0,100]` → **percent**.
- If two count columns exist (methylated & unmethylated), coverage = sum(counts); fraction = methylated/coverage.
- If only one count/coverage column exists, treat it as coverage.

### 2) Convert to per-sample 5‑column Metilene tables

**Target schema per sample (TSV):**
```
chrom    start    end    methylation_fraction(0-1)    coverage(int)
```

> The 5‑column per‑sample tables are retained for provenance but metilene itself works on a **merged matrix of methylation fractions** (next step).

### 3) Build the merged methylation matrix (fractions per sample)

Create a manifest from `groups.tsv` and check matching files:
```bash
cut -f1 groups.tsv > work/samples.list

while read s; do
  test -s "work/per_sample/${s}.metilene.tsv" || { echo "Missing converted table for sample: ${s}"; exit 1; }
done < work/samples.list
```

Extract only the **fraction column (4th)** from each per-sample table and merge on genomic intervals using `bedtools unionbedg`:
```bash
(echo -e "#chrom\tstart\tend\t$(paste -sd '\t' work/samples.list)";  cat work/merged_methylation.bed) > work/merged_with_header.tsv
```


### 4) Run metilene (DMR mode)

```bash
metilene -a ${name of groupA} -b ${name of groupB} -t 1 work/merged_methylation.bed > results/dmr_results.txt 2> results/metilene.stderr.log
```

### 5) Filter significant DMRs and export BED

```bash
# Keep FDR q<=0.05 and |Delta methylation| >= 0.25
awk 'BEGIN{OFS="\t"}
     $4!~/^#/ && $4<=0.05 && ( $5>=0.25 || $5<=-0.25 ) {print $1,$2,$3,$4,$5}'      results/dmr_results.txt > results/significant_dmrs.txt

cut -f1-3 results/significant_dmrs.txt > results/significant_dmrs.bed

# Summary
total=$(grep -vc '^#' results/dmr_results.txt || true)
sig=$(wc -l < results/significant_dmrs.txt || echo 0)
hypo=$(awk '$5<0' results/significant_dmrs.txt | wc -l | tr -d ' ')
hyper=$(awk '$5>0' results/significant_dmrs.txt | wc -l | tr -d ' ')
awk 'BEGIN{print "DMR length stats (bp):";}
     {len=$3-$2; if(min==""||len<min)min=len; if(len>max)max=len; sum+=len; n++}
     END{if(n>0){printf "min=%d\nmax=%d\nmean=%.2f\n",min,max,sum/n} else {print "No DMRs"}}'      results/significant_dmrs.bed > results/dmr_summary.txt
printf "Total DMRs: %d\nSignificant (q<=0.05 & |Δ|>=0.25): %d\nHyper: %s\nHypo: %s\n"        "$total" "$sig" "$hyper" "$hypo" >> results/dmr_summary.txt
```

### 6) Optional: DMC mode (single‑CpG tests)

```bash
metilene -a ${name of groupA} -b ${name of groupB} -t 1 work/merged_methylation.bed > results/dmc_results.txt 2> results/metilene_dmc.stderr.log

awk 'BEGIN{OFS="\t"} $4!~/^#/ && $4<=0.05 && ( $5>=0.10 || $5<=-0.10 ) {print $1,$2,$3,$4,$5}' results/dmc_results.txt > results/significant_dmcs.txt
cut -f1-3 results/significant_dmcs.txt > results/significant_dmcs.bed
```

### 7) Visualization (quick, optional)

**Volcano-like (Δmethylation vs –log10(q))**
```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt
df = pd.read_csv('results/dmr_results.txt', sep='\t', header=None, comment='#')
df = df.iloc[:, :5]; df.columns=['chr','start','end','q','diff']
df = df[(df['q']>0) & np.isfinite(df['q'])]
plt.figure(figsize=(6,5))
plt.scatter(df['diff'], -np.log10(df['q']), s=4, alpha=0.5)
plt.xlabel('Delta methylation (case - control)')
plt.ylabel('-log10(q)')
plt.tight_layout(); plt.savefig('plots/volcano.png', dpi=200)
```

**DMR length histogram**
```python
import pandas as pd, matplotlib.pyplot as plt
s = pd.read_csv('work/dmr.lengths', header=None)[0]
plt.figure(figsize=(6,4))
plt.hist(s, bins=50)
plt.xlabel('DMR length (bp)'); plt.ylabel('Count')
plt.tight_layout(); plt.savefig('plots/length_hist.png', dpi=200)
```

---

## Parameters & Tuning

- **`-m / --min-cpg`**: minimum CpGs per DMR (use 1 for DMC mode).
- **`-d / --max-dist`**: maximum distance between adjacent CpGs within a DMR.
- **FDR threshold** (`q<=0.05`) and **effect size** (|Δ|≥0.25) are common starting points; adapt to sample size/biology.
- Balance group sizes; extremely unbalanced designs reduce power.

---


## Troubleshooting

- **"Input not sorted"**: re‑sort each per-sample table and the prepared fraction BedGraphs (`sort -k1,1 -k2,2n`).
- **Missing samples**: ensure per-sample file basenames match sample IDs in `groups.tsv`.
- **All zeros/ones**: confirm the fraction vs percent conversion and that invalid rows were filtered.
- **Chromosome naming mismatches**: standardize to a single scheme (`chr1` vs `1`) across all samples.
