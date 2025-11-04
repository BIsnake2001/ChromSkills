---
name: track-generation
description: This skill generates normalized BigWig (.bw) tracks from BAM files for ATAC-seq and ChIP-seq visualization. It handles genome detection, chromosome size matching, normalization (RPM or fold-change), and Tn5 offset correction automatically.
---

## Overview

This skill converts filtered BAM files into normalized signal tracks (BigWig) for genome browser visualization.  
It supports both ATAC-seq and ChIP-seq datasets, automatically detecting genome assembly and chromosome size files.

**Steps:**

1. **Identify BAM files** in the current directory or user-specified input.  
2. **Detect genome assembly** from BAM header and locate `chrom.sizes`.  The location of `chrom.sizes` is  `$(dirname $(which homer))/../data/genomes/${GENOME}/chrom.sizes` 
3. **Ensure chromosome naming consistency** between BAM and `chrom.sizes`, for example, with or without "chr". If not consistent, then generate a `chrom.size` file with the consistent name according to the existing `chrom.sizes` file.
4. **Normalize all tracks** to 1 million mapped reads (RPM normalization).  
5. **For ATAC-seq**, apply Tn5 offset correction (+4/−5) and generate normalized BigWig (RPM).  
6. **For ChIP-seq**, generate two BigWig files:  
   - RPM-normalized track (signal intensity) using `bedtools` and `bedGraphToBigWig`
   - log₂(IP/Input) fold-change track using `bamCompare` [optional, depends on user, default don't generate] 
7. **Save all BigWig files** to the `tracks` folder; intermediate files to `temp/`.
8. **Always use the filtered BAM file** if available

---

## Decision Tree

### 1. Identify Input BAM Files

Detect BAM files automatically:

- For paired-end data, ensure that BAM files are coordinate-sorted and indexed (`.bai`).  index sample.sorted.bam

---

### 2. Detect Genome Assembly

Extract genome build from BAM header:
```bash
samtools view -H sample.bam | grep "@SQ" | head -n 1
```

Determine assembly:
- `chr1` → `hg38` or `mm10`
- `1` → `hg19` or `mm9`

Locate `chrom.sizes`:
```bash
GENOME=hg38
CHROMSIZES=$(dirname $(which homer))/../data/genomes/${GENOME}/chrom.sizes
```

If chromosome naming is inconsistent, regenerate a consistent version of `chrom.sizes` to match BAM.

---

### 3. Generate BigWig for ATAC-seq

**Goal:** Create RPM-normalized BigWig scaled to 1M mapped reads with Tn5 offset correction.

**Steps:**

1. **Shift Tn5 cut sites:**
   ```bash
   bamToBed -i sample.bam |    awk 'BEGIN{OFS="\t"} $6=="+"{$2=$2+4} $6=="-"{$3=$3-5} {print $0}' > temp/sample_shifted.bed
   ```

2. **Generate bedGraph normalized to 1M mapped reads:**
   ```bash
   SCALE=$(echo "1000000/$(samtools view -c -F 260 sample.bam)" | bc -l)
   bedtools genomecov -bg -scale $SCALE -i temp/sample_shifted.bed -g $CHROMSIZES > temp/sample.bg
   ```

3. **Convert bedGraph to BigWig:**
   ```bash
   bedGraphToBigWig temp/sample.bg $CHROMSIZES tracks/sample.RPM.bw
   ```

✅ **Output:** `tracks/sample.RPM.bw`  
✅ **Normalization:** RPM (Reads per Million mapped reads)  
✅ **Correction:** Tn5 (+4/−5)

---

### 4. Generate BigWig for ChIP-seq

#### (A) RPM-normalized signal track (same as ATAC-seq normalization)

**Goal:** Generate signal intensity BigWig scaled to 1M mapped reads.

```bash
SCALE=$(echo "1000000/$(samtools view -c -F 260 IP.bam)" | bc -l)
bedtools genomecov -bg -scale $SCALE -ibam IP.bam -g $CHROMSIZES > temp/IP.bg
bedGraphToBigWig temp/IP.bg $CHROMSIZES tracks/IP.RPM.bw
```

✅ **Output:** `tracks/IP.RPM.bw`  
✅ **Normalization:** RPM (1M mapped reads)

---

#### (B) Fold-change track using `bamCoverage`

**Goal:** Compute log₂(IP/Input) normalized BigWig for enrichment visualization.

**Command:**
```bash
bamCompare   -b IP.bam   -c Input.bam   --outFileName tracks/IP_vs_Input.foldchange.bw   --outFileFormat bigwig   --scaleFactorsMethod readCount   --operation log2   --binSize 25   --extendReads   --normalizeUsing None
```

✅ **Output:** `tracks/IP_vs_Input.foldchange.bw`  
✅ **Normalization:** log₂(IP/Input) fold-change  
✅ **Bin size:** 25 bp (adjustable)

---

### 5. Notes

- Always use **filtered BAMs** (duplicates/blacklist removed).  
- Check for **consistent chromosome naming** before `bedGraphToBigWig`.  
- **Genome autodetection** can be overridden by setting `$GENOME` manually.  
- All intermediate files stored in `temp/`, final outputs in `tracks/`.  
- Requires: `samtools`, `bedtools`, `bedGraphToBigWig`, `deepTools`.  
