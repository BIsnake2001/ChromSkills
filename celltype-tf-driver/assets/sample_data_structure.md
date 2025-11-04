# Sample Data Structure for Testing

This document describes the expected file structure and format for testing the cell type transcription factor driver analysis.

## Directory Structure

```
test_data/
├── input/
│   ├── ATAC_peak_fibroblast.bed
│   ├── ATAC_peak_cardiomyocyte.bed
│   ├── H3K27ac_peak_fibroblast.bed
│   ├── H3K27ac_peak_cardiomyocyte.bed
│   ├── hg38.fa
│   └── human.meme
├── config.yaml
└── expected_output/
    ├── key_driver_tfs.csv
    ├── motif_enrichment_results.csv
    ├── differential_accessibility_results.bed
    └── analysis_summary.md
```

## File Formats

### BED Files (ATAC-seq and H3K27ac peaks)

**Format:** Standard 3-column BED format
```
chr1	1000	2000	peak_1	100	+
chr1	5000	6000	peak_2	50	-
chr2	3000	4000	peak_3	75	.
```

**Columns:**
1. Chromosome name
2. Start position (0-based)
3. End position (0-based, exclusive)
4. Peak name (optional)
5. Score (optional)
6. Strand (optional)

**Requirements:**
- Minimum 1000 peaks per cell type
- Peak sizes between 100-1000 bp
- Valid chromosome names matching genome assembly

### Genome FASTA File

**Format:** Standard FASTA format
```
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
>chr2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

**Requirements:**
- Must include .fai index file (created automatically)
- Chromosome names must match BED files
- Should be the same genome assembly used for peak calling

### MEME Motif Database

**Format:** MEME motif database format
```
MEME version 5

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000

MOTIF MA0001.1 Ahr::Arnt

letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
  0.200000  0.300000  0.300000  0.200000
  0.400000  0.100000  0.100000  0.400000
  0.300000  0.200000  0.200000  0.300000
  0.200000  0.300000  0.300000  0.200000
  0.400000  0.100000  0.100000  0.400000
  0.300000  0.200000  0.200000  0.300000
```

**Requirements:**
- Must include transcription factor names
- Should be species-specific (human.meme, mouse.meme, etc.)
- Can be from JASPAR, HOCOMOCO, or custom databases

## Sample Configuration File

```yaml
# test_data/config.yaml
cell_type_a: "fibroblast"
cell_type_b: "cardiomyocyte"
genome_assembly: "hg38"
species: "human"

input_directory: "./input"
output_directory: "./results"

atac_peak_a: "ATAC_peak_fibroblast.bed"
atac_peak_b: "ATAC_peak_cardiomyocyte.bed"
h3k27ac_peak_a: "H3K27ac_peak_fibroblast.bed"
h3k27ac_peak_b: "H3K27ac_peak_cardiomyocyte.bed"
genome_fasta: "hg38.fa"
meme_database: "human.meme"

# Analysis parameters
fdr_threshold: 0.05
fold_change_threshold: 1.5
fimo_pvalue_threshold: 1e-4

# Integration weights
weight_motif_enrichment: 0.4
weight_accessibility_change: 0.3
weight_enhancer_activity: 0.2
weight_expression_change: 0.1
```

## Expected Output Files

### Key Driver TFs (`key_driver_tfs.csv`)
```csv
TF_Name,Composite_Score,Motif_Enrichment,Accessibility_Change,Enhancer_Activity,Expression_Change
GATA4,3.2,2.1,0.8,0.3,0.0
TBX5,2.8,1.8,0.7,0.3,0.0
NKX2-5,2.5,1.5,0.6,0.4,0.0
MEF2C,2.1,1.2,0.5,0.4,0.0
```

### Motif Enrichment Results (`motif_enrichment_results.csv`)
```csv
motif_id,tf_name,foreground_count,background_count,enrichment_score,p_value,fdr,odds_ratio
MA0035.2,GATA4,45,12,1.91,1.2e-5,0.001,3.75
MA0080.3,TBX5,38,15,1.34,2.3e-4,0.012,2.53
MA0074.1,NKX2-5,32,14,1.19,5.6e-4,0.025,2.29
```

### Differential Accessibility (`differential_accessibility_results.bed`)
```
chr1	1000	2000	peak_1	fibroblast	-2.1	1.3e-6	1.2e-4	lost
chr1	5000	6000	peak_2	cardiomyocyte	1.8	2.1e-5	8.7e-4	gained
chr2	3000	4000	peak_3	cardiomyocyte	2.3	4.5e-7	3.1e-5	gained
```

## Generating Test Data

### Small Test Dataset
For quick testing, you can create minimal test files:

**Create minimal BED files:**
```bash
# ATAC peaks for fibroblast (100 peaks)
for i in {1..100}; do
  echo -e "chr1\t$((i*1000))\t$((i*1000+500))\tpeak_fib_$i\t100\t+"
done > ATAC_peak_fibroblast.bed

# ATAC peaks for cardiomyocyte (100 peaks with some overlap)
for i in {1..100}; do
  if [ $i -le 50 ]; then
    # Different peaks for first 50
    echo -e "chr2\t$((i*1000))\t$((i*1000+500))\tpeak_cm_$i\t100\t+"
  else
    # Same peaks for last 50 (shared regions)
    echo -e "chr1\t$(((i-50)*1000))\t$(((i-50)*1000+500))\tpeak_shared_$i\t100\t+"
  fi
done > ATAC_peak_cardiomyocyte.bed
```

### Using Public Data
For more realistic testing, use publicly available datasets:

1. **ENCODE Project**: ATAC-seq and H3K27ac data for various cell types
2. **Roadmap Epigenomics**: Epigenomic data across multiple cell types
3. **JASPAR**: Public motif databases
4. **UCSC Genome Browser**: Genome FASTA files

## Validation Commands

```bash
# Validate configuration
python scripts/validate_inputs.py --config test_data/config.yaml

# Run analysis
python scripts/tf_driver_analysis.py \
  --atac-a test_data/input/ATAC_peak_fibroblast.bed \
  --atac-b test_data/input/ATAC_peak_cardiomyocyte.bed \
  --h3k27ac-a test_data/input/H3K27ac_peak_fibroblast.bed \
  --h3k27ac-b test_data/input/H3K27ac_peak_cardiomyocyte.bed \
  --genome-fasta test_data/input/hg38.fa \
  --meme-db test_data/input/human.meme \
  --output-dir test_data/results \
  --cell-type-a fibroblast \
  --cell-type-b cardiomyocyte
```

## Troubleshooting Test Data

### Common Issues
- **BED file coordinates invalid**: Ensure start < end and coordinates are 0-based
- **Genome mismatch**: Verify chromosome names match between BED and FASTA files
- **MEME database empty**: Check database contains valid motifs
- **Low peak count**: Ensure each cell type has sufficient peaks (>1000 recommended)

### Solutions
- Use `bedtools validate` to check BED files
- Create genome index with `samtools faidx`
- Download pre-validated MEME databases from JASPAR
- Use public datasets with sufficient data coverage