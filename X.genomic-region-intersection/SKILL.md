---
name: genomic-region-intersection
description: This skill performs intersection and overlap analysis between genomic interval files (e.g., BED, narrowPeak, or GTF). It identifies overlapping or nearby regions such as TF binding sites overlapping gene annotations, promoters, or other regulatory regions using bedtools or pybedtools.

---

# Genomic Interval Intersection and Overlap Analysis

## Overview

1. Prepare two or more genomic files in BED, GTF, or narrowPeak format.  
2. Use `bedtools intersect` or `bedtools window` to identify overlapping or nearby genomic features.  
3. Optionally calculate overlap statistics (e.g., percent overlap, number of shared regions).  
4. Generate summary reports and visualizations of overlaps.  

---

## Decision Tree

### 1. Input File Preparation

Accepts any standard genomic interval files containing at least three columns: `chromosome`, `start`, and `end`.

**Supported formats:**
- BED, narrowPeak, broadPeak, GTF, GFF
- If using GTF, convert to BED first for simpler intersection:
  ```bash
  awk '$3=="gene"{print $1"\t"$4-1"\t"$5"\t"$9}' genes.gtf > genes.bed
  ```

**Check format validity:**

```bash
head -n 5 TF_peaks.bed
```

---

### 2. Intersection and Overlap Calculation

Use **bedtools** to identify overlapping regions.

#### Example 1. Simple overlap between two BED files

```bash
bedtools intersect -a TF_peaks.bed -b genes.bed -wa -wb > TF_gene_overlap.bed
```

**Flags explanation:**

* `-a`: primary file (e.g., TF peaks)
* `-b`: secondary file (e.g., gene annotations)
* `-wa`: write original entries from file A
* `-wb`: write original entries from file B that overlap

#### Example 2. Require minimum overlap fraction (e.g., 50%)

```bash
bedtools intersect -a TF_peaks.bed -b enhancers.bed -f 0.5 -wa -wb > TF_enhancer_overlap.bed
```

#### Example 3. Count number of overlaps per region

```bash
bedtools intersect -a TF_peaks.bed -b genes.bed -c > TF_gene_overlap_counts.bed
```

#### Example 4. Find nearby (not necessarily overlapping) features

```bash
bedtools window -a TF_peaks.bed -b TSS.bed -w 1000 > TF_near_TSS.bed
```

This finds TF peaks within ±1 kb of a TSS.

---

### 3. Visualization of Overlap Results

After generating overlap files, visualize the extent of intersection.

**Example in Python:**

```python
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("TF_gene_overlap_counts.bed", sep="\t", header=None)
plt.hist(data.iloc[:, -1], bins=30)
plt.xlabel('Number of overlaps per region')
plt.ylabel('Count')
plt.title('Distribution of TF-gene overlaps')
plt.show()
```

**Example in R:**

```R
data <- read.table("TF_gene_overlap_counts.bed", header=FALSE)
hist(data$V4, breaks=30, main="TF-gene overlap count distribution", xlab="Overlap count")
```

---

### 4. Summary and Statistics

Summarize results by counting overlaps and calculating overlap fractions.

```bash
# Count total overlapping regions
wc -l TF_gene_overlap.bed

# Unique peaks overlapping at least one gene
cut -f1-3 TF_gene_overlap.bed | sort | uniq | wc -l
```

**Optional overlap fraction calculation (Python):**

```python
import pandas as pd
df = pd.read_csv("TF_gene_overlap.bed", sep="\t", header=None)
overlap_fraction = len(df) / len(open("TF_peaks.bed").readlines())
print(f"Overlap fraction: {overlap_fraction:.2%}")
```

---

### 5. Example Use Cases

**(1) TF peaks overlapping gene promoters:**

```bash
bedtools intersect -a TF_peaks.bed -b promoters.bed -wa -wb > TF_promoter_overlap.bed
```

**(2) Histone mark peaks overlapping genes:**

```bash
bedtools intersect -a H3K4me3.bed -b genes.gtf -wa -wb > H3K4me3_gene_overlap.bed
```

**(3) Enhancer-promoter proximity (within 2 kb):**

```bash
bedtools window -a enhancers.bed -b promoters.bed -w 2000 > enhancer_promoter_pairs.bed
```

---

### 6. Troubleshooting

| Issue                | Description                           | Solution                             |
| -------------------- | ------------------------------------- | ------------------------------------ |
| File format mismatch | GTF/GFF cannot be parsed directly     | Convert to BED using `awk`           |
| Empty overlap        | Chromosome naming differs (chr1 vs 1) | Normalize chromosome names           |
| Missing bedtools     | Not installed                         | `conda install -c bioconda bedtools` |

---

### 7. Output Files

| File                  | Description                                    |
| --------------------- | ---------------------------------------------- |
| `*_overlap.bed`       | Overlapping genomic regions from both files    |
| `*_counts.bed`        | File A with overlap counts                     |
| `*_nearby.bed`        | File A regions within given distance of file B |
| `overlap_summary.txt` | Optional summary statistics                    |

---

### 8. Best Practices

* Ensure consistent genome assemblies and chromosome naming.
* Filter peaks or annotations by confidence before intersection.
* Use strandedness options (`-s`, `-S`) if strand-specific.
* Save parameters and version info for reproducibility.

---

## Example Workflow Summary

```bash
# Step 1: Intersect TF peaks and gene annotations
bedtools intersect -a TF_peaks.bed -b genes.bed -wa -wb > TF_gene_overlap.bed

# Step 2: Summarize overlap statistics
wc -l TF_gene_overlap.bed

# Step 3: Visualize results
python plot_overlap_hist.py TF_gene_overlap.bed
```
---

