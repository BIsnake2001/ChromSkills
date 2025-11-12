---
name: integrative-methylation-expression
description: This skill performs correlation analysis between differential methylation and differential gene expression, identifying genes with coordinated epigenetic regulation. It provides preprocessing, integration, and enrichment workflows, using promoter-level methylation–expression relationships.

---

# Integrative Methylation–Expression Correlation Analysis

## Overview

This skill integrates **differential methylation** and **differential expression** datasets to reveal coordinated epigenetic regulation patterns.

- **Preprocess** differential methylation and expression tables into a standard format.
- **Integrate** methylation and expression data by promoter proximity.
- **Detect the gene name version** (Entrez ID or Ensembl ID), and make sure the consistence when excute merging.
- **Calculate correlation** between methylation change and expression fold change.
- **Classify patterns** such as hypermethylation–downregulation or hypomethylation–upregulation.
- **Perform pathway enrichment** for genes with coordinated changes.

---

## Decision Tree

### 1. Preprocessing Differential Data

Convert heterogeneous methylation or expression outputs (from tools like `methylKit`, `DSS`, `DESeq2`, or `edgeR`) into a standard tab-delimited table.

#### Input
- **Differential methylation table** (from methylKit/DSS/custom)
- **Differential expression table** (from DESeq2/edgeR/custom)

#### Command Example
```bash
# Methylation data preprocessing
Rscript -e '
library(dplyr); library(data.table);
data <- fread("methyl_diff.tsv");
data <- data %>%
  select(chr, start, end, qvalue, meth.diff) %>%
  rename(pvalue = qvalue, meth_diff = meth.diff);
fwrite(data, "processed_methyl.tsv", sep="\t");
'

# Expression data preprocessing
Rscript -e '
library(dplyr);
expr <- read.csv("DESeq2_results.csv");
expr <- expr %>%
  select(gene, padj, log2FoldChange) %>%
  rename(pvalue = padj, log2fc = log2FoldChange);
write.table(expr, "processed_expr.tsv", sep="\t", row.names=FALSE, quote=FALSE);
'
```

#### Key Parameters
- `--type`: `methyl` or `expr`
- `--source`: `methylkit`, `dss`, `deseq2`, `edger`, or `custom`
- Columns standardized as:
  - Methylation: `chr`, `start`, `end`, `pvalue`, `meth_diff`
  - Expression: `gene`, `pvalue`, `log2fc`

#### Output
- `processed_methyl.tsv`
- `processed_expr.tsv`

---

### 2. Integrative Analysis

Join methylation and expression datasets to detect genes whose promoter methylation correlates with expression change.

#### Concept
- Promoter: ±2kb window around TSS.
- Genes overlapped by significantly methylated promoter CpGs are paired with their expression data.
- Spearman correlation measures association.

#### Command Example
```bash
Rscript -e '
library(dplyr); library(ggplot2); library(GenomicRanges); library(clusterProfiler); library(org.Hs.eg.db);
meth <- read.table("processed_methyl.tsv", header=TRUE);
expr <- read.table("processed_expr.tsv", header=TRUE);

sig_meth <- subset(meth, pvalue < 0.05 & abs(meth_diff) > 10)
sig_expr <- subset(expr, pvalue < 0.05 & abs(log2fc) > log2(1.5))

# Example promoter annotation (pseudo)
# promoters <- makeGRangesFromDataFrame(annotation_df, keep.extra.columns=TRUE)
# methyl_gr <- makeGRangesFromDataFrame(sig_meth)
# overlaps <- findOverlaps(methyl_gr, promoters)

# Merge methylation and expression data by gene symbol
integrated <- merge(sig_expr, sig_meth, by="gene", all=FALSE)

# Spearman correlation
cor_res <- cor.test(integrated$meth_diff, integrated$log2fc, method="spearman")

# Scatter plot
p <- ggplot(integrated, aes(x=meth_diff, y=log2fc)) +
  geom_point(alpha=0.6, color="blue") +
  geom_smooth(method="lm", se=TRUE, color="red") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(title="Promoter Methylation vs Gene Expression",
       subtitle=paste("Spearman rho =", round(cor_res$estimate,3),
                      "p =", round(cor_res$p.value,4)))
ggsave("integrative_correlation_plot.pdf", p, width=8, height=6)
'
```

#### Output
- `*_integrated_results.tsv`
- `*_correlation_plot.pdf/png`
- `*_pattern_counts.tsv`
- `*_summary_stats.tsv`

#### Key Parameters
| Parameter | Description | Default |
|------------|--------------|----------|
| `--promoter_window` | Promoter region size | 2000 bp |
| `--methyl_pval` | Methylation significance cutoff | 0.05 |
| `--methyl_diff` | Methylation difference threshold (%) | 10 |
| `--expr_pval` | Expression significance cutoff | 0.05 |
| `--expr_fc` | Expression FC threshold | 1.5 |
| `--genome` | Genome assembly | hg38 |

---

### 3. Pattern Classification

Genes are categorized by methylation–expression relationships:

| Pattern | Biological Meaning |
|----------|--------------------|
| **Hypermethylation + Downregulation** | Epigenetic silencing (tumor suppressor-like) |
| **Hypomethylation + Upregulation** | Activation (oncogene-like) |
| **Hypermethylation + Upregulation** | Indirect effect |
| **Hypomethylation + Downregulation** | Alternative regulation |

#### Example Code
```r
integrated$pattern <- "No significant pattern"
integrated$pattern[integrated$meth_diff > 0 & integrated$log2fc < 0] <- "Hypermethylation + Downregulation"
integrated$pattern[integrated$meth_diff < 0 & integrated$log2fc > 0] <- "Hypomethylation + Upregulation"

table(integrated$pattern)
```

---

### 4. Pathway Enrichment

Enrichment identifies functional pathways for genes with coordinated methylation–expression changes.

#### Command Example
```bash
Rscript -e '
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
data <- read.table("integrated_results.tsv", header=TRUE)
genes <- subset(data, pattern %in% c("Hypermethylation + Downregulation","Hypomethylation + Upregulation"))$gene
ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ego <- enrichGO(ids$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1)
p <- dotplot(ego, showCategory=15) + labs(title="GO Enrichment of Coordinated Genes")
ggsave("go_enrichment.pdf", p, width=8, height=6)
'
```

#### Output
- `*_go_bp.tsv`, `*_kegg.tsv`
- `*_dotplot.pdf/png`
- `*_enrichment_summary.tsv`

#### Typical Parameters
| Option | Default | Description |
|--------|----------|-------------|
| `--pattern` | both | Target pattern (`hyper_down`, `hypo_up`, `both`) |
| `--pvalue_cutoff` | 0.05 | P-value cutoff |
| `--qvalue_cutoff` | 0.1 | Adjusted q-value cutoff |
| `--organism` | hsapiens | Organism database |

---

## Parameter Guidelines

**Methylation**
- `p < 0.05`, difference ≥ 10% recommended for WGBS  
**Expression**
- `p < 0.05`, `|log2FC| > 0.58` (≈1.5×) typical

**Promoter window:** ±2kb around TSS  
**Minimum genes for enrichment:** ≥50 overlap recommended

---

## Quality Control

- Ensure consistent gene IDs between methylation and expression tables.
- Harmonize chromosome naming (`chr1` vs `1`).
- Verify significant overlap (>50 genes).
- Correlation significance: *p* < 0.05.
- FDR threshold for enrichment: *q* < 0.1.

---

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|---------------|-----------|
| No overlapping genes | Wrong gene mapping | Check annotation file and promoter window |
| Weak correlation | Low-quality input | Adjust cutoffs or verify batch correction |
| Few enrichment terms | Small gene set | Relax significance or combine both patterns |

---

## Biological Interpretation

- **Cancer:** Hypermethylation silences tumor suppressors; hypomethylation activates oncogenes.  
- **Development:** Dynamic promoter methylation drives stage-specific expression.  
- **Disease:** Identify methylation–expression-linked pathways in relevant tissues.
