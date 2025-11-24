---  
name: global-methylation-profile 
description: This skill performs genome-wide DNA methylation profiling. It supports single-sample and multi-sample workflows to compute methylation density distributions, genomic feature distribution of the methylation profile, and sample-level clustering/PCA. Use it when you want to systematically characterize global methylation patterns from WGBS or similar per-CpG methylation call files.  
---

# Global DNA Methylation Profiling

## Overview

Main steps include:

- Refer to the **Inputs & Outputs** section to check available inputs and design the output structure.  
- **Always prompt user** for genome assembly used.
- **Always prompt user** for which columns are methylation fraction/percent.
- Prepare the raw methylation BED into 6-column standard format BED file received by HOMER.
- Annoate the CpG sites to genomic features with HOMER.
- Analyze and visualize the genomic feature distribution of methylations for each sample.
- For multi-sample datasets, prepare the matrix of methylation data.
- Compute and visualize genome-wide methylation density distributions.
- Perform PCA and hierarchical clustering to assess sample similarity based on global methylation.

---

## When to use this skill

Use the **global-methylation-profiling** skill when you want to:

- Characterize **global DNA methylation status** of one or multiple samples (e.g. normal vs tumor, different cell types).  
- Compare broad methylation patterns across samples:  
  - Are some samples globally hypo-/hyper-methylated?  
  - Are certain chromosomes or genomic regions more strongly affected?  
- Explore genomic feature of your methylation dataset (e.g. promoter hypomethylation, gene body hypermethylation).  
- Perform **unsupervised clustering/PCA** to see if samples separate by condition based on genome-wide methylation patterns.

---

## Inputs & Outputs

### Inputs
`<sample.bed`

### Outputs

```bash
global_methylation_profile/
  stats/
    ...
  plots/
    sample1_genomic_feature_pie.pdf
    sample2_genomic_feature_pie.pdf
    ... # Other samples
    allSamples_methylation_density_overlay.pdf
    PCA_scatterplot.pdf
    sample_correlation_heatmap.pdf
    ...
  logs/
  temp/ # all the temp files
```

---

## Decision Tree

### Step 1: Prepare the raw methylation BED into 6-column standard format BED file received by HOMER.

```bash
for sample in samples;do
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, "peak_"NR, "*", "+"}' sample.bed > homer_sample.bed
done
```

### Step 2: Annoate the CpG sites to genomic features with HOMER.

- for each homer_sample.bed

Call:

- mcp__homer-tools__homer_annotate_peaks

with:
- `peaks_path`: 6-column standard format BED file from Step 1.
- `genome`: Provide by user.
- `output_path`: Output path of the annotated file


### Step 3: Analyze and visualize Genomic Feature Distribution of CpGs

Call:
- mcp__methyl-tools_global_methyl_genomic_features

with:

- `sample_paths`: List of sample path. (HOMER output for each sample from Step 2)
- `stats_dir`: Stats output directory.
- `plots_dir`: Plots output directory.

### Step 4: For multi-sample datasets, prepare the matrix of methylation data

(1) Generate a sample meta TSV file containing with columns: sample_id,file_path,treatment

(2) Call: mcp__methyl-tools_prepare_methylation_matrix
with:
- `sample_table`: sample meta TSV file
- `output_dir`: Output directory
- `freqC_col`: 1-based column index for methylation value (fraction or percent of Cm, provided by user)

### Step 5: Visualization

Call: 
- mcp__methyl-tools_global_methyl_pca_clustering

with:
methylation_matrix: Path to `methylation_matrix.tsv` (generated in Step 4)
sample_table: sample meta TSV file (generated in  Step 4)
out_dir: Output directory

