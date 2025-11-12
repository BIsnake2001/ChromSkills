---
name: integrative-methylation-expression
description: This skill should be used when users need to integrate differential DNA methylation and differential gene expression data to identify coordinated epigenetic regulation patterns. It provides workflows for correlation analysis, promoter methylation-expression relationships, pattern classification, and functional enrichment of genes with coordinated changes.
---

# Integrative Methylation-Expression Analysis

## Overview

This skill enables comprehensive integration of differential DNA methylation and differential gene expression data to identify genes with coordinated epigenetic regulation. It focuses on promoter methylation changes and their relationship to gene expression, providing statistical correlation analysis, pattern classification, and functional enrichment.

## When to Use This Skill

Use this skill when:
- Investigating epigenetic regulation of gene expression
- Identifying genes with coordinated methylation-expression changes
- Analyzing promoter methylation effects on transcription
- Performing functional enrichment of epigenetically regulated genes
- Integrating WGBS/RRBS methylation data with RNA-seq expression data

## Quick Start

To begin integrative methylation-expression analysis:

1. **Preprocess differential data** (if needed):
   ```bash
   # Preprocess methylation data
   Rscript scripts/preprocess_differential_data.R \
     methyl_results.tsv \
     processed_methyl.tsv \
     --type methyl \
     --source methylkit

   # Preprocess expression data
   Rscript scripts/preprocess_differential_data.R \
     expr_results.csv \
     processed_expr.tsv \
     --type expr \
     --source deseq2
   ```

2. **Perform integrative analysis**:
   ```bash
   Rscript scripts/integrative_analysis.R \
     processed_methyl.tsv \
     processed_expr.tsv \
     output_prefix \
     --promoter_window 2000 \
     --methyl_pval 0.05 \
     --methyl_diff 10 \
     --expr_pval 0.05 \
     --expr_fc 1.5 \
     --genome hg19
   ```

3. **Run pathway enrichment** (optional):
   ```bash
   Rscript scripts/pathway_enrichment.R \
     output_prefix_integrated_results.tsv \
     enrichment_output \
     --pattern both
   ```

## Data Preprocessing

Convert various differential data formats to standardized format using `scripts/preprocess_differential_data.R`.

### Supported Input Formats

**Differential Methylation:**
- methylKit output (`meth.diff`, `qvalue` columns)
- DSS output (`diff.Methy`, `fdr` columns)
- Custom BED format with methylation differences

**Differential Expression:**
- DESeq2 output (`log2FoldChange`, `padj` columns)
- edgeR output (`logFC`, `FDR` columns)
- Custom tables with gene symbols and fold changes

### Preprocessing Parameters

- `--type`: Specify data type (`methyl` or `expr`)
- `--source`: Data source (`deseq2`, `methylkit`, `edger`, `dss`, `custom`)
- Column mapping parameters for custom formats

## Integrative Analysis

Perform correlation and pattern analysis using `scripts/integrative_analysis.R`.

### Analysis Parameters

- `--promoter_window`: Promoter region size from TSS (default: 2000bp)
- `--methyl_pval`: Methylation p-value threshold (default: 0.05)
- `--methyl_diff`: Methylation difference threshold (default: 10%)
- `--expr_pval`: Expression p-value threshold (default: 0.05)
- `--expr_fc`: Expression fold change threshold (default: 1.5x)
- `--genome`: Genome assembly (hg19, hg38, mm10)
- `--gene_annotation`: Custom gene annotation file

### Output Files

- `*_integrated_results.tsv`: Genes with coordinated changes
- `*_correlation_plot.pdf/png`: Methylation vs expression scatter plot
- `*_pattern_plot.pdf/png`: Pattern classification visualization
- `*_pattern_counts.tsv`: Counts of methylation-expression patterns
- `*_summary_stats.tsv`: Statistical summary

## Pathway Enrichment

Analyze functional enrichment for coordinated changes using `scripts/pathway_enrichment.R`.

### Enrichment Parameters

- `--pattern`: Pattern to analyze (`hyper_down`, `hypo_up`, `both`)
- `--pvalue_cutoff`: P-value cutoff (default: 0.05)
- `--qvalue_cutoff`: Q-value cutoff (default: 0.1)
- `--min_genes`: Minimum genes per category (default: 5)
- `--max_categories`: Maximum categories to show (default: 20)
- `--organism`: Organism (`hsapiens`, `mmusculus`)

### Enrichment Outputs

- `*_go_bp.tsv`: GO Biological Process enrichment
- `*_go_mf.tsv`: GO Molecular Function enrichment
- `*_go_cc.tsv`: GO Cellular Component enrichment
- `*_kegg.tsv`: KEGG pathway enrichment
- `*_*_dotplot.pdf/png`: Enrichment dotplots
- `*_enrichment_summary.tsv`: Summary statistics

## Pattern Classification

### Expected Biological Patterns

1. **Hypermethylation + Downregulation**:
   - Classical epigenetic silencing
   - Tumor suppressor genes in cancer
   - Developmental gene regulation

2. **Hypomethylation + Upregulation**:
   - Epigenetic activation
   - Oncogene activation in cancer
   - Cell type-specific gene expression

3. **Unexpected Patterns**:
   - Hypermethylation + Upregulation: Indirect effects
   - Hypomethylation + Downregulation: Alternative regulation

### Statistical Analysis

- **Spearman correlation**: Non-parametric correlation analysis
- **Pattern counts**: Statistical significance of pattern enrichment
- **Functional enrichment**: Pathway analysis of coordinated changes

## Parameter Guidelines

### Significance Thresholds

**Methylation:**
- p-value: 0.05 (default)
- Methylation difference: 10% (default)
- For WGBS: 5-15% difference
- For targeted: 10-25% difference

**Expression:**
- p-value: 0.05 (default)
- Fold change: 1.5x (default)
- For RNA-seq: 1.5-2x fold change
- For microarray: 1.2-1.5x fold change

### Promoter Window Sizes

- **Standard**: ±2kb from TSS (default)
- **Narrow**: ±1kb for focused analysis
- **Broad**: ±5kb for extended regulatory regions

## Quality Control

### Input Data Validation

- Verify gene symbol consistency between datasets
- Check chromosome naming conventions
- Validate statistical significance thresholds
- Ensure proper genomic coordinate systems

### Analysis Quality Metrics

- Minimum overlapping genes: 50
- Correlation significance: p < 0.05
- Enrichment FDR: q < 0.1
- Pattern consistency across biological replicates

## Troubleshooting

### Common Issues

1. **No overlapping genes**:
   - Check gene symbol mapping
   - Verify promoter window size
   - Consider alternative gene annotations

2. **Weak correlation**:
   - Adjust significance thresholds
   - Consider different genomic contexts
   - Check for batch effects

3. **Limited enrichment**:
   - Lower p-value thresholds
   - Increase minimum gene set size
   - Try different pattern combinations

### Performance Optimization

- Pre-filter differential data by significance
- Use efficient data structures (data.table)
- Consider parallel processing for enrichment
- Monitor memory usage with large gene sets

## Resources

### Scripts

- `scripts/integrative_analysis.R`: Main integrative analysis script
- `scripts/preprocess_differential_data.R`: Data format conversion
- `scripts/pathway_enrichment.R`: Functional enrichment analysis

### References

- `references/workflow_guide.md`: Complete workflow documentation
- `references/statistical_methods.md`: Statistical methods reference

## Biological Interpretation

### Cancer Studies
- Focus on tumor suppressor and oncogene patterns
- Consider tissue-specific methylation patterns
- Integrate with clinical data when available

### Developmental Studies
- Analyze stage-specific regulation
- Consider dynamic methylation changes
- Integrate with chromatin accessibility data

### Disease Studies
- Focus on disease-relevant pathways
- Consider cell type-specific effects
- Validate with independent datasets
