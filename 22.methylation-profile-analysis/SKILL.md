---
name: methylation-profile-analysis
description: Analyze DNA methylation profiles around genomic regions using methylKit. This single document embeds concise R code hints (no external .R files) for preprocessing, profiling, and comparing profiles. It asks the user which columns are methylation fraction and coverage, converts inputs to a 5-column BED for methylKit, and highlights key functions/parameters only.
---

# Methylation Profile Analysis

## Overview
- **Confirm columns**: Ask which columns are methylation fraction/percent and coverage.
- **Preprocess**: Convert input → 5-column BED: `chr, start, end, methyl_percent, coverage` (percent 0–100).
- **Profile**: Bin methylation around regions (±flank, fixed bin size), aggregate mean±SE.
- **Visualize**: Plot mean profile with ribbon and center line.
- **Compare (optional)**: Merge two profiles, compute differences, paired t-test across bins.

## Decision Tree

### 1) Preprocess input → 5-column BED (for methylKit), and 3-column BED (for target regions)
**Goal**: Standardize any input to `chr  start  end  methyl_percent  coverage` for methylation and any input of target regions to `chr start end`.

**Ask the user**: Which column is **methylation** (fraction 0–1 or percent 0–100)? Which column is **coverage**?

**Key R hints**:
```r
library(data.table)
DT <- fread("INPUT.tsv", header = TRUE)
# choose by user indices
methyl <- as.numeric(DT[[methyl_col_i]])
cov    <- as.integer(DT[[coverage_col_i]])
looks_fraction <- median(methyl, na.rm=TRUE) <= 1
methyl_percent <- if (looks_fraction) methyl*100 else methyl
OUT <- data.table(chr=DT[[chr_i]], start=DT[[start_i]], end=DT[[end_i]],
                  methyl_percent=pmin(pmax(methyl_percent,0),100), coverage=cov)[coverage>=min_coverage]
fwrite(OUT, "OUTPUT.5col.bed", sep="\t", col.names=FALSE)
```
**Notes**: Ensure percent in [0,100]; strip any `track` line if present; tolerate header/no-header.

---

### 2) Build methylation profiles around regions
**Goal**: For each region, collect CpGs within ±`flank_size`, bin by `bin_size`, compute mean per bin; aggregate across regions.

**Key parameters**: `flank_size = 2000`, `bin_size = 50`, `min_coverage = 10`.

**Key R hints**:
```r
library(data.table); library(ggplot2)
M <- fread("METH.5col.bed", header=FALSE); setnames(M, c("chr","start","end","methyl_percent","coverage"))
R <- fread("REGIONS.bed", header=FALSE); setnames(R, c("chr","start","end"))
R[, center := (start+end)%/%2]
# per region → relative positions and binning
# bin: floor((CpG_start - center)/bin_size)*bin_size
# aggregate per region, then across regions: mean and SE
# minimal plotting
p <- ggplot(overall, aes(bin, mean_methyl)) + geom_line() +
     geom_ribbon(aes(ymin=mean_methyl-se_methyl, ymax=mean_methyl+se_methyl), alpha=.2) +
     geom_vline(xintercept=0, linetype="dashed")
```
**Tips**: Split methylation by chromosome for speed; ensure consistent genome build (chr vs 1).

---

### 3) Compare two profiles (optional)
**Goal**: Merge by `bin`, compute `diff = mean_c2 - mean_c1`, propagate SE, run paired t-test across bins.

**Key R hints**:
```r
library(data.table); library(ggplot2)
A <- fread("cond1_profile.tsv"); B <- fread("cond2_profile.tsv")
X <- merge(A, B, by="bin", suffixes=c("_c1","_c2"))
X[, `:=`(methyl_diff = mean_methyl_c2 - mean_methyl_c1,
         se_diff = sqrt(se_methyl_c1^2 + se_methyl_c2^2))]
T <- t.test(X$mean_methyl_c1, X$mean_methyl_c2, paired=TRUE)
# quick plots
ggplot(X, aes(bin)) + geom_line(aes(y=mean_methyl_c1)) + geom_line(aes(y=mean_methyl_c2))
ggplot(X, aes(bin, methyl_diff)) + geom_line() +
  geom_ribbon(aes(ymin=methyl_diff-se_diff, ymax=methyl_diff+se_diff), alpha=.2)
```
**Report**: save merged table, plots, and a small `comparison_stats.tsv` with `t_statistic`, `p_value`, `n_bins`.

---

## Parameter Guidelines
| Context   | Flank | Bin  | Min cov |
|-----------|-------|------|---------|
| TF peaks  | ±2 kb | 50bp | 10x     |
| Promoters | ±1 kb | 50bp | 10x     |
| Enhancers | ±5 kb | 100bp| 5x      |
| Motifs    | ±0.5kb| 10–20| 10x     |

## Quality Control
- Consistent chromosome naming (`chr1` vs `1`).
- Methylation stored as **percent (0–100)** after preprocessing.
- Prefer ≥100 regions with data; ≥5 CpGs per region on average.

## Output


## Troubleshooting
- **No CpGs in window** → check genome build, increase flank, lower coverage.
- **All 0/100** → input likely fraction not scaled; re-run preprocessing.
- **Asymmetric profile** → confirm region centering (peak summit vs motif center).

## Notes
- Snippets are *usage hints* and must be adapted to your paths and column indices.
