
# BEDTools Intersect Command Reference

## Core Command

### bedtools intersect

Find overlapping regions between two BED/GTF/VCF files.

**Basic Syntax:**

```bash
bedtools intersect -a <fileA> -b <fileB> [options] > output.bed
```

* **-a**: Primary file (e.g., peaks, features)
* **-b**: Secondary file(s) to intersect with
* **Output:** intervals in `-a` that overlap with intervals in `-b` (by default)

---

## Essential Options

| Option       | Description                                                                              | Example                                                            |
| :----------- | :--------------------------------------------------------------------------------------- | :----------------------------------------------------------------- |
| `-wa`        | Write the original entry in A for each overlap                                           | `bedtools intersect -a a.bed -b b.bed -wa > overlap_in_a.bed`      |
| `-wb`        | Write the original entry in B for each overlap                                           | `bedtools intersect -a a.bed -b b.bed -wb > overlap_pairs.bed`     |
| `-u`         | Write each entry in A once if it overlaps B                                              | `bedtools intersect -a a.bed -b b.bed -u > unique_hits.bed`        |
| `-v`         | Report entries in A **with no overlap** in B                                             | `bedtools intersect -a a.bed -b b.bed -v > nonoverlap.bed`         |
| `-c`         | For each entry in A, report the number of overlaps with B                                | `bedtools intersect -a a.bed -b b.bed -c > count_overlaps.bed`     |
| `-wo`        | Write the original A and B entries plus the number of overlapping bases                  | `bedtools intersect -a a.bed -b b.bed -wo > overlap_basepairs.bed` |
| `-f <float>` | Minimum overlap fraction of A required (0–1)                                             | `bedtools intersect -a a.bed -b b.bed -f 0.5`                      |
| `-r`         | Require reciprocal overlap (both A and B must overlap by at least `-f` fraction)         | `bedtools intersect -a a.bed -b b.bed -f 0.5 -r`                   |
| `-s`         | Require same strand                                                                      | `bedtools intersect -a a.bed -b b.bed -s`                          |
| `-S`         | Require opposite strand                                                                  | `bedtools intersect -a a.bed -b b.bed -S`                          |
| `-sorted`    | Use for pre-sorted inputs (faster and memory-efficient)                                  | `bedtools intersect -a a.sorted.bed -b b.sorted.bed -sorted`       |
| `-loj`       | Left outer join: report all A intervals, even if no overlap in B (adds `-1` when no hit) | `bedtools intersect -a a.bed -b b.bed -loj > left_join.bed`        |

---

## Common Use Cases

### 1. Intersect TF Peaks with Gene Annotations

```bash
bedtools intersect -a tf_peaks.bed -b genes.gtf -wa -wb > peaks_with_genes.bed
```

> Returns TF peaks and their overlapping gene annotation entries.

### 2. Find Unique TF Peaks (Non-overlapping)

```bash
bedtools intersect -a tf_peaks.bed -b histone_marks.bed -v > unique_tf_peaks.bed
```

> Returns TF peaks that **do not** overlap with histone mark regions.

### 3. Count Overlaps per Feature

```bash
bedtools intersect -a genes.bed -b peaks.bed -c > gene_overlap_counts.bed
```

> Appends a column with the number of overlapping peaks for each gene.

### 4. Reciprocal 50% Overlap

```bash
bedtools intersect -a human_peaks.bed -b mouse_peaks.bed -f 0.5 -r > conserved_peaks.bed
```

> Reports only those peaks that overlap at least 50% in both species.

### 5. Calculate Basepair Overlap Length

```bash
bedtools intersect -a a.bed -b b.bed -wo > overlap_lengths.bed
```

> Adds a final column with the basepair overlap length between each pair.

---

## Advanced Options

| Option         | Description                                                                                   |
| :------------- | :-------------------------------------------------------------------------------------------- |
| `-filenames`   | Prepend filenames (A and B) to each output line                                               |
| `-split`       | Treat “split” BED12 features (e.g., exons in transcripts) as separate intervals               |
| `-nonamecheck` | Disable name checking when chromosome names differ between A and B                            |
| `-e`           | Require at least one base overlap                                                             |
| `-header`      | Include header lines from input BED/GTF/VCF files                                             |
| `-wao`         | Similar to `-wo`, but reports all A records even when no overlap (0 in overlap length column) |

---

## Example Workflows

### Gene Annotation Overlap

```bash
bedtools intersect -a peaks.bed -b genes.gtf -wa -wb | awk '$NF=="exon"' > peaks_in_exons.bed
```

### TF Co-binding Analysis

```bash
bedtools intersect -a TF1_peaks.bed -b TF2_peaks.bed -u > cobound_sites.bed
```

### Promoter Overlap (TSS ±2kb)

```bash
bedtools intersect -a promoters_2kb.bed -b histone_H3K4me3.bed -c > promoter_H3K4me3_counts.bed
```

### Using Multiple B Files

```bash
bedtools intersect -a peaks.bed -b mark1.bed mark2.bed mark3.bed -c > peak_overlap_summary.bed
```

---

## Output Explanation

| Column                     | Description                               |
| :------------------------- | :---------------------------------------- |
| 1–3                        | Chromosome, Start, End of A feature       |
| 4+                         | Fields from A and/or B depending on flags |
| Last Column (`-wo`/`-wao`) | Overlap length (bp)                       |
| `-c`                       | Count of overlapping features             |
| `-v`                       | Only A records with **no** overlap        |

---

## Performance Tips

* **Sort inputs first** using:

  ```bash
  sort -k1,1 -k2,2n file.bed > file.sorted.bed
  ```
* Use `-sorted` for large datasets.
* Use `-c` or `-u` when only counts or existence of overlap is needed (saves memory).
* Use `bgzip` and `tabix` for compressed, indexed inputs to speed up queries.

---

## Troubleshooting

| Issue                  | Cause                                                  | Solution                                                      |
| :--------------------- | :----------------------------------------------------- | :------------------------------------------------------------ |
| No overlaps found      | Coordinates differ (e.g., different genome assemblies) | Ensure both files use same assembly (e.g., hg38 vs hg19)      |
| “chromosome not found” | Name mismatch (chr1 vs 1)                              | Use `-nonamecheck` or rename chromosomes                      |
| Memory/Speed issues    | Very large files                                       | Sort inputs and use `-sorted`                                 |
| Missing columns        | Input format not recognized                            | Ensure BED/GTF files are tab-delimited and properly formatted |

---
