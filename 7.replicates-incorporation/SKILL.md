name: replicates-incorporation
description: This skill handles replicate incorporation for ATAC-seq and ChIP-seq data, including pooled BAM generation, normalized bigWig creation, and reproducible peak merging based on IDR or overlap results.
---

# Replicates Incorporation Skill

## Overview

This skill integrates replicate-level ChIP-seq and ATAC-seq data to generate reproducible tracks and optimal peak sets:

- **Identify genome build** from BAM header and load proper `chrom.sizes`. The location of `chrom.sizes` is  `$(dirname $(which homer))/../data/genomes/${GENOME}/chrom.sizes` 
- **Ensure chromosome naming consistency** between BAM and `chrom.sizes`, for example, with or without "chr". If not consistent, then generate a `chrom.size` file with the consistent name according to the existing `chrom.sizes` file.
- **Pool replicate BAM files** into a single BAM file per sample.
- **Generate normalized bigWig files** (RPM-scaled) using `bedtools` and `bedGraphToBigWig`.
- **For ATAC-seq**, perform Tn5 cut-site adjustment (+4/−5 bp).
- **For ChIP-seq**, generate both normalized signal (`RPM.bw`) and log₂(IP/Input) fold change (`foldchange.bw`) tracks. Among which the `foldchange.bw` is optional by the user, default is not perform this step.
- **Integrate peaks from replicates** using IDR or overlap-based filtering (keep IDR<0.05).
- **Save all track outputs** to the `track/` folder and merged peaks to `qc_results/`.

---

## Workflow Decision Tree

### 1. Identify Genome Assembly

Automatically extract genome information from BAM headers:
```bash
samtools view -H sample_rep1.bam | grep "@SQ" | head -n 1
# Example output: @SQ SN:chr1 LN:248956422 → genome=hg38
GENOME=$(infer_genome_from_header.sh sample_rep1.bam)
CHROM_SIZES="$(dirname $(which homer))/../data/genomes/${GENOME}/chrom.sizes"
```

---

### 2. Pool Replicate BAMs

Merge BAM files from multiple replicates:
```bash
samtools merge -f temp/{sample}.pooled.bam rep1.bam rep2.bam
samtools index temp/{sample}.pooled.bam
```

---

### 3. Generate Normalized bigWig (RPM) Tracks

#### For ATAC-seq:
Adjust Tn5 cut-site (+4/−5 bp) and normalize to 1M mapped reads:
```bash
bedtools genomecov -bg -scale $(echo "1000000 / $(samtools view -c -F 260 temp/{sample}.pooled.bam)" | bc -l)   -ibam temp/{sample}.pooled.bam | awk '{OFS="\t"; if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5}; print}' | bedGraphToBigWig -trackName {sample}.RPM -color=0,0,255 -lineWidth=2   -clip stdin ${CHROM_SIZES} track/{sample}.pooled.RPM.bw
```

#### For ChIP-seq:
Normalize IP signal to 1M mapped reads:
```bash
bedtools genomecov -bg -scale $(echo "1000000 / $(samtools view -c -F 260 temp/{sample}_IP.pooled.bam)" | bc -l)   -ibam temp/{sample}_IP.pooled.bam | bedGraphToBigWig stdin ${CHROM_SIZES} track/{sample}.pooled.RPM.bw
```

Generate log₂(IP/Input) fold change:
```bash
bamCoverage -b temp/{sample}_IP.pooled.bam -c temp/{sample}_Input.pooled.bam   -o track/{sample}.foldchange.bw --scaleFactorsMethod readCount   --normalizeUsing log2 --operation ratio --binSize 25
```

---

### 4. Merge Replicate Peaks

#### For Narrow ChIP-seq peaks and ATAC-seq:
Use IDR results from `qc_results` to select reproducible peaks, names with `{sample}.optimal_peak.narrowPeak`, ouput file should be located in `peaks` directory. Note that: the .txt file ontains the scaled IDR value, min(int(log2(-125IDR), 1000). e.g. peaks with an IDR of 0 have a score of 1000, idr 0.05 have a score of int(-125log2(0.05)) = 540, and idr 1.0 has a score of 0. So you should filter the fifth columns by `awk '{if($5 >= 540) print $0}' sample-idr`. Make sure the output format is a valid narrowPeak format with proper column information.

#### For Broad Histone Marks:
Merge replicate peaks by ≥50% overlap:
```bash
bedtools intersect -a rep1.broadPeak -b rep2.broadPeak -f 0.5 -r | sort -k1,1 -k2,2n | bedtools merge -i - > qc_results/{sample}.optimal_peak.broadPeak
```
---

### 5. Output Structure

| Output Type | Location | Naming |
|--------------|-----------|--------|
| Pooled BAM | `temp/{sample}.pooled.bam` | Intermediate |
| Normalized bigWig | `track/{sample}.pooled.RPM.bw` | Scaled to 1M reads |
| Fold Change bigWig (ChIP-seq only) | `track/{sample}.foldchange.bw` | log₂(IP/Input) |
| Optimal Peaks | `qc_results/{sample}.optimal_peak.{narrow|broad}Peak` | IDR < 0.05 or overlap ≥50% |

---

## Example Usage

### Example 1: ATAC-seq Replicates
```bash
samtools merge -f temp/K562_ATAC.pooled.bam K562_ATAC_rep1.bam K562_ATAC_rep2.bam

bedtools genomecov -bg -scale $(echo "1000000 / $(samtools view -c -F 260 temp/K562_ATAC.pooled.bam)" | bc -l)   -ibam temp/K562_ATAC.pooled.bam | awk '{OFS="\t"; if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5}; print}' | bedGraphToBigWig stdin $(dirname $(which homer))/../data/genomes/hg38/chrom.sizes track/K562_ATAC.pooled.RPM.bw
```

### Example 2: ChIP-seq Replicates
```bash
samtools merge -f temp/CTCF_IP.pooled.bam CTCF_rep1.bam CTCF_rep2.bam
samtools merge -f temp/CTCF_Input.pooled.bam Input_rep1.bam Input_rep2.bam

bedtools genomecov -bg -scale $(echo "1000000 / $(samtools view -c -F 260 temp/CTCF_IP.pooled.bam)" | bc -l)   -ibam temp/CTCF_IP.pooled.bam | bedGraphToBigWig stdin $(dirname $(which homer))/../data/genomes/mm10/chrom.sizes track/CTCF.pooled.RPM.bw

bamCompare -b temp/CTCF_IP.pooled.bam -c temp/CTCF_Input.pooled.bam   -o track/CTCF.foldchange.bw --normalizeUsing log2 --operation ratio --binSize 25
```

---

## Best Practices

- **Use pooled tracks** for visualization and differential analysis.
- **Keep individual replicate tracks** for QC and reproducibility evaluation.
- **Use IDR ≤ 0.05** for reproducible narrow ChIP-seq peaks and ATAC-seq.
- **Use overlap ≥50% ** for broad histone mark peaks.
- **Always normalize** to 1 million mapped reads for consistency.

---

**Note:** Temporary files are stored in `temp/`, normalized tracks in `track/`, and reproducible peak files in `peaks/`.
