---
name: functional-enrichment
description: Perform GO and KEGG functional enrichment using HOMER from genomic regions (BED/narrowPeak/broadPeak) or gene lists, and produce R-based barplot/dotplot visualizations.

---

# Functional Enrichment (HOMER + R)

## Overview

- **Validate input**: Accept BED/peak files with genomic coordinates or gene lists; check format and genome assembly.
- **Map regions to genes**: Convert regions to a unique gene set using HOMER `annotatePeaks.pl`.
- **Run GO enrichment**: Use HOMER `findGO.pl` (or `annotatePeaks.pl -go`) for BP/MF/CC.
- **Run KEGG enrichment**: Use HOMER `findGO.pl -kegg` (or `annotatePeaks.pl -kegg`).
- **Collect outputs**: Save tidy tables for downstream plotting and a compact summary of top terms.
- **Visualize in R**: Create barplots and dotplots (GO/KEGG) with `ggplot2` from standardized outputs.
- **QC & troubleshooting**: Provide checks for genome mismatch, chromosome naming, and low-signal inputs.

## Inputs & Outputs

**Inputs (choose one):**
- `*.bed|*.narrowPeak|*.broadPeak` with at least 3 columns (chrom, start, end); extra columns allowed.
- `gene_list.txt` with one official gene symbol per line (no header). And an optional `gene_list_background.txt` with one official gene symbol per line (no header).

**Required parameters:**
- `--genome <hg38|hg19|mm10|mm9|...>` HOMER genome key (must be installed).
- Optional: `--background <bg_gene_list.txt>` custom background set.

**Outputs (directory layout):**
```
results/
  {run}/
    input/
      peaks.bed or gene_list.txt
    tables/
      genes.txt
      go_bp.tsv, go_mf.tsv, go_cc.tsv
      kegg.tsv
    plots/
      GO_barplot.pdf, GO_dotplot.pdf
      KEGG_barplot.pdf, KEGG_dotplot.pdf
    logs/
      annotate.log, findGO.log
```

## Requirements

- HOMER installed and target genome added (e.g., `configureHomer.pl -install hg38`).
- R (≥3.6) with packages: `ggplot2`, `dplyr`, `readr`, `stringr`, `forcats`.

## Decision Tree (key examples only)

### A) Input is a region file (BED/narrowPeak/broadPeak)

**A1. Sanity check**
```bash
head -n 3 peaks.bed
```

**A2. Map regions → genes**
```bash
annotatePeaks.pl peaks.bed hg38 -annStats results/{run}/logs/annStats.txt   > results/{run}/tables/annotated.tsv 2> results/{run}/logs/annotate.log

awk 'NR>1{print $2}' results/{run}/tables/annotated.tsv | sed 's/[;,].*//' | sort -u   > results/{run}/tables/genes.txt
```

**A3. (Optional) Build a background**
```bash
annotatePeaks.pl all_regions.bed hg38 > results/{run}/tables/all_annotated.tsv
awk 'NR>1{print $2}' results/{run}/tables/all_annotated.tsv | sed 's/[;,].*//' | sort -u   > results/{run}/tables/background.txt
```

**A4. GO & KEGG enrichment**
```bash
findGO.pl results/{run}/tables/genes.txt human results/{run}/tables > results/{run}/tables/go_results.tsv ## if no background gene list provided
```

> **Alternative direct from BED**  
> `annotatePeaks.pl peaks.bed hg38 -go results/{run}/tables/go_dir -genomeOntology`  
> `annotatePeaks.pl peaks.bed hg38 -kegg results/{run}/tables/kegg_dir`

### B) Input is a gene list

**B1. GO & KEGG enrichment**
```bash
findGO.pl results/{run}/input/genes_list.txt human results/{run}/tables > results/{run}/tables/go_results.tsv ## if no background gene list provided
```

**B2. Optional background**
```bash
findGO.pl results/{run}/input/genes_list.txt human results/{run}/tables -bp results/{run}/input/genes_list_background.txt > results/{run}/tables/go_results.tsv ## if background gene list provided
```

## Visualization in R (barplot & dotplot)

> The plotting assumes `findGO.pl`-style columns (e.g., `Term`, `BH/FDR` or `P-value`, `NumDEInCat`, `NumInCat`). Adjust column names if needed.

**Minimal example (GO BP KEGG barplot)**

only check the first line of the file to get column name if any bug occurs
```r
library(readr); library(dplyr); library(ggplot2); library(forcats)

df <- read_tsv("results/{run}/tables/go_results.tsv", comment = "#", show_col_types = FALSE)
fdr <- intersect(c("BH/FDR","FDR","Q-value","Q"), names(df)); p <- intersect(c("P-value","P","P.value"), names(df))
df$FDR2 <- if (length(fdr)) df[[fdr[1]]] else df[[p[1]]]
df$score <- -log10(pmax(df$FDR2, 1e-300))

p <- df %>% slice_min(FDR2, n = 15) %>% mutate(Term = forcats::fct_reorder(Term, score)) %>%
  ggplot(aes(x = Term, y = score)) + geom_col() + coord_flip() +
  labs(x = NULL, y = "-log10(FDR or P)", title = "GO BP enrichment") + theme_bw(12)

ggsave("results/{run}/plots/GO_barplot.pdf", p, width = 7, height = 5)
```

## Quick Run Template

```bash
run=my_run; genome=hg38
mkdir -p results/${run}/{input,tables,plots,logs}

# from regions:
cp peaks.bed results/${run}/input/
annotatePeaks.pl results/${run}/input/peaks.bed ${genome} -annStats results/${run}/logs/annStats.txt > results/${run}/tables/annotated.tsv
awk 'NR>1{print $2}' results/${run}/tables/annotated.tsv | sed 's/[;,].*//' | sort -u > results/${run}/tables/genes.txt

# enrichment:
findGO.pl results/{run}/input/genes_list.txt human results/{run}/tables > results/{run}/tables/go_results.tsv ## if no background gene list provided

# then run the minimal R snippets above to create plots
```

## Notes & Best Practices

- **Genome & naming**: Ensure the HOMER genome key matches the species; chromosome naming must be consistent (`chr1` vs `1`).
- **BED format**: Tab-delimited, ≥3 columns, 0-based coordinates, no header.
- **Multiple testing**: Prefer FDR (BH) if provided; otherwise fallback to P-value.
- **Background set**: `-bg` helps reduce bias; choose a reasonable universe (e.g., all expressed or all accessible regions → genes).
- **Direct-from-BED**: `annotatePeaks.pl -go/-kegg` is convenient; the gene-list route yields uniform TSVs for plotting.

## Troubleshooting

- **Many NAs after annotation**: Check genome version, chromosome naming, BED formatting, and headers.
- **Empty/weak enrichment**: Ensure sufficient genes (suggest ≥50), verify species of symbols, tune thresholds or background.
- **Column name drift**: HOMER versions may differ; adjust R column mappings if needed.
