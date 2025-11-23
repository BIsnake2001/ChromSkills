---
name: loop-annotation-homer
description: This skill annotates chromatin loops using HOMER’s annotateInteractions.pl, including enhancer/promoter assignments, CTCF-peak overlap. It automatically constructs enhancer and promoter sets when missing and outputs standardized loop categories.
---

# Loop Annotation with HOMER annotateInteractions.pl

## Overview
This skill performs loop annotation for Hi-C/HiChIP/ChIA-PET interaction data using HOMER’s annotateInteractions.pl. It identifies regulatory and structural loop types (E–E, E–P, P–P, CTCF-CTCF) and computes CTCF motif orientation. Steps:

- Refer to **Inputs & Outputs** to verify necessary files.
- **Always prompt user** if required files are missing.
- **Always prompt user** for genome assembly used.
- Build enhancers.bed if absent (ATAC + H3K27ac).
- Build promoters.bed if absent (.tss annotation). The .tss annotation are located in `${HOMER_data}/genomes/${genome/${genome}.tss`
- Standardize the format of the `.bedpe` file as the input of `annotateInteractions.pl`
- Run `annotateInteractions.pl` with feature sets.
- Visualization

## When to use this skill
Use when you need:
- Regulatory loop landscape analysis.
- Enhancer–promoter mapping from chromatin loops.
- Structural loop analysis via CTCF orientation.
- Integration with ATAC/H3K27ac/TSS/CTCF datasets.
- Upstream to differential loop testing or expression integration.

## Inputs & Outputs

### Inputs
Required:
- loops.bedpe
- ctcf_peaks.bed
- CTCF.motif
- genome version (user must provide)

Optional:
- enhancers.bed
- promoters.bed
- ATAC_peaks.bed
- H3K27ac_peaks.bed
- .tss or .gtf gene annotation

### Outputs

```bash
loop_annotation/
    logs/
        annotateInteractions.log
    annotations/
        interactionAnnotation.txt
        lengthDist.txt
        featureEnrichment.txt
        pairwiseFeatureEnrichment.txt # assign feature pairs to 0x0 0x1 and so on, represent the feature pairs like CTCF-CTCF, E-P
        ... # other outputs by annotateInteractions.pl
    features/
        enhancers.bed
        promoters.bed
    plots:
        loop_type.pdf
        lengthDist.pdf
```

## Decision Tree

### Step 1 — Validate inputs

If enhancers.bed or promoters.bed missing, go on to Step 2, otherwise go on to step 4 directly.

### Step 2 — Build enhancers.bed (if missing)

```bash
bedtools intersect -a ATAC_peaks.bed -b H3K27ac_peaks.bed > enhancers_raw.bed
sort -k1,1 -k2,2n enhancers_raw.bed | bedtools merge -i - > enhancers.bed
```

### Step 3 — Build promoters.bed (if missing)

```bash
awk 'BEGIN{OFS="\t"} {
 chr=$2; tss=$3+2000; gene=$1;
 start=tss-1000; if(start<0)start=0;
 end=tss+100;
 print chr,start,end,gene
}' input.tss > promoters.bed
```
### Step 4 — Standardize the input file of annotateInteractions.pl

Standardize the loop file of `annotateInteractions.pl` to contain the following columns, fill the column with 0 if not provided in the vailable loop file.

1) Interaction ID (must be unique)
2) Peak ID for region 1
3) chr for region 1
4) start position for region 1
5) end position for region 1
6) strand for region 1
7) total reads for region 1
8) Peak ID for region 2
9) chr for region 2
10) start position for region 2
11) end position for region 2
12) strand for region 2
13) total reads for region 2
14) Distance between regions (or "interchromosomal")
15) Interaction Reads (total Hi-C reads connecting the regions)
16) Expected Interaction Reads (total expected Hi-C reads based on background model)
17) Modified Z-score
18) Natural log of the p-value for the interaction (binomial)
19) False Discovery Rate (based on Benjamini correction)
20) Circos Thickness (used for visualization by Circos)

- **example code**
```python
input_file = "loop.bedpe"
output_file = "loops_standardized.bedpe"
df = pd.read_csv(input_file, sep="\t", skiprows=[1], header=0)

columns_out = [
    "InteractionID",
    "PeakID_1", "Chr1", "Start1", "End1", "Strand1", "Reads1",
    "PeakID_2", "Chr2", "Start2", "End2", "Strand2", "Reads2",
    "Distance",
    "InteractionReads",
    "ExpectedReads",
    "ModifiedZ",
    "ln_pvalue",
    "FDR",
    "CircosThickness"
]

out = pd.DataFrame(index=df.index)

# Interaction ID
out["InteractionID"] = ["int{}".format(i+1) for i in range(len(df))]

# Region 1
out["PeakID_1"] = df["name"].astype(str) + "_1"
out["Chr1"] = df["#chr1"]
out["Start1"] = df["x1"]
out["End1"] = df["x2"]
out["Strand1"] = df["strand1"]
out["Reads1"] = df["observed"]

# Region 2
out["PeakID_2"] = df["name"].astype(str) + "_2"
out["Chr2"] = df["chr2"]
out["Start2"] = df["y1"]
out["End2"] = df["y2"]
out["Strand2"] = df["strand2"]
out["Reads2"] = df["observed"]

# Distance
same_chr = df["#chr1"] == df["chr2"]
out["Distance"] = same_chr * (abs(df["x1"] - df["y1"]))
out.loc[~same_chr, "Distance"] = "interchromosomal"

# InteractionReads
out["InteractionReads"] = df["observed"]
out["ExpectedReads"] = df["expectedBL"]
out["FDR"] = df['fdrBL']

# 以下列在原文件中没有 → 填 NA
out["ModifiedZ"] = 0
out["ln_pvalue"] = 0

out["CircosThickness"] = 0
out = out[columns_out]
out.to_csv(output_file, sep="\t", index=False)
```
Output: `loops_standardized.bedpe`

### Step 5 — Run annotateInteractions.pl
```bash
# adjust the -res parameter if user provide custom resolution
annotateInteractions.pl <loops_standardized> <genome> annotations -res 10000 -p ctcf_peaks.bed enhancers.bed promoters.bed -cpu 4 > annotateInteractions.log
```

### Step 6 — Classify and visualize loop types
Python example:
```python
df_anno = pd.csv("interactionAnnotation.txt", sep='\t')
# The peak Link column contains the 0x0 0x1 and so on, represent the feature pairs like CTCF-CTCF, E-P
df_anno_target =[['Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', 'Peak Links']]

# draw the bar plot of the 'Peak Links' columns
# plot the length distribution from the lengthDist.txt
...
```
Output: `barplot_loop_type.pdf` `lengthDist.pdf`


## Advanced Usage
- Use -minDist/-maxDist to restrict loop ranges.
