# Chromatin State Interpretation Guide

## Standard Chromatin State Definitions

### 15-State Model (ENCODE Standard)

| State | Name | Histone Marks | Genomic Features | Biological Function |
|-------|------|---------------|------------------|---------------------|
| 1 | Active Promoter | H3K4me3++, H3K27ac++ | TSS, CpG islands | Active transcription initiation |
| 2 | Weak Promoter | H3K4me3+ | TSS | Poised or weak promoters |
| 3 | Inactive/Repressed Promoter | H3K27me3+ | TSS | Polycomb-repressed genes |
| 4 | Strong Enhancer | H3K27ac++, H3K4me1+ | Intergenic, intronic | Active enhancer elements |
| 5 | Strong Enhancer 2 | H3K27ac++, H3K4me1++ | Intergenic | Very active enhancers |
| 6 | Weak/Poised Enhancer | H3K4me1+ | Intergenic | Poised or weak enhancers |
| 7 | Weak/Poised Enhancer 2 | H3K4me1++ | Intergenic | Developmental enhancers |
| 8 | Insulator | CTCF+ | Boundary elements | Chromatin domain boundaries |
| 9 | Transcribed | H3K36me3+ | Gene bodies | Active transcription elongation |
| 10 | Transcribed 2 | H3K36me3++ | Gene bodies | Highly transcribed genes |
| 11 | Repressed | H3K27me3++ | Intergenic | Polycomb-repressed regions |
| 12 | Repressed 2 | H3K27me3+++ | Intergenic | Strongly repressed regions |
| 13 | Heterochromatin | Low signals | Repeat regions | Constitutive heterochromatin |
| 14 | Repetitive/CNV | Variable | Repeat regions | Copy number variations |
| 15 | Quiescent/Low | Low signals | Intergenic | Inactive chromatin |

### 8-State Model (Simplified)

| State | Name | Histone Marks | Genomic Features |
|-------|------|---------------|------------------|
| 1 | Active Promoter | H3K4me3, H3K27ac | TSS regions |
| 2 | Weak Promoter | H3K4me3 | TSS regions |
| 3 | Strong Enhancer | H3K27ac, H3K4me1 | Intergenic regions |
| 4 | Weak Enhancer | H3K4me1 | Intergenic regions |
| 5 | Insulator | CTCF | Boundary elements |
| 6 | Transcribed | H3K36me3 | Gene bodies |
| 7 | Repressed | H3K27me3 | Various regions |
| 8 | Quiescent | Low signals | Intergenic regions |

## Histone Mark Combinations for State Identification

### Promoter States

**Active Promoter (State 1):**
- High H3K4me3 (>0.8)
- High H3K27ac (>0.7)
- Low H3K27me3 (<0.2)
- Enriched at transcription start sites

**Weak Promoter (State 2):**
- Moderate H3K4me3 (0.4-0.8)
- Low H3K27ac (<0.4)
- Low H3K27me3 (<0.2)
- Found at weakly expressed genes

**Repressed Promoter (State 3):**
- Low H3K4me3 (<0.3)
- High H3K27me3 (>0.6)
- Associated with polycomb target genes

### Enhancer States

**Strong Enhancer (States 4-5):**
- High H3K27ac (>0.7)
- Moderate to high H3K4me1 (0.5-0.9)
- Low H3K4me3 (<0.3)
- Enriched in intergenic and intronic regions

**Weak/Poised Enhancer (States 6-7):**
- Low H3K27ac (<0.4)
- Moderate H3K4me1 (0.3-0.6)
- May represent developmental or poised enhancers

### Transcription States

**Transcribed Regions (States 9-10):**
- High H3K36me3 (>0.7)
- Low promoter marks (H3K4me3 <0.3)
- Enriched throughout gene bodies
- Correlates with gene expression levels

### Repressive States

**Polycomb Repressed (States 11-12):**
- High H3K27me3 (>0.7)
- Low active marks (<0.2)
- Associated with developmental repression

**Heterochromatin (State 13):**
- Low signals across all marks (<0.2)
- Enriched in repeat-rich regions
- Constitutive heterochromatin

## Genomic Feature Enrichment Patterns

### Gene Regions
- **Promoters**: States 1-3
- **5' UTR**: States 1-2
- **Exons**: States 9-10
- **Introns**: States 4-7, 9-10
- **3' UTR**: States 9-10

### Regulatory Elements
- **Enhancers**: States 4-7
- **Insulators**: State 8
- **CTCF sites**: State 8
- **Super-enhancers**: States 4-5 with large domains

### Repetitive Elements
- **LINE elements**: States 13-15
- **SINE elements**: States 13-15
- **LTR elements**: States 13-15
- **Satellite repeats**: State 13

## Biological Context Interpretation

### Cell Type Specificity
- **Ubiquitous states**: States 1, 9, 13 (housekeeping functions)
- **Cell-type specific**: States 4-7 (enhancer states)
- **Developmental**: States 3, 11-12 (polycomb states)

### Developmental Dynamics
- **Stable states**: Promoters, transcribed regions
- **Dynamic states**: Enhancers, repressed regions
- **Plastic states**: Poised enhancers, bivalent promoters

### Disease Associations
- **Cancer**: Altered enhancer states, gained/lost promoter states
- **Neurodevelopmental**: Altered polycomb states
- **Autoimmune**: Altered enhancer states in immune cells

## Validation Methods

### Expression Correlation
- Promoter states (1-2) should correlate with gene expression
- Transcribed states (9-10) should correlate with expression levels
- Repressed states (3, 11-12) should anti-correlate with expression

### Functional Validation
- Enhancer states should validate with reporter assays
- Promoter states should match RNA-seq expression
- Insulator states should match CTCF ChIP-seq

### Comparative Analysis
- Compare with ENCODE chromatin states
- Validate with orthogonal methods (ATAC-seq, DNAse-seq)
- Check conservation across species

## Common Interpretation Pitfalls

### Over-interpretation
- Avoid assigning precise biological functions to every state
- Some states may represent technical artifacts
- Consider the limitations of the histone mark panel

### Context Dependence
- State definitions vary by cell type
- Developmental stage affects state interpretation
- Disease states may have unique chromatin patterns

### Technical Considerations
- Input data quality affects state calling
- Number of states should match biological complexity
- Bin size affects state resolution and interpretation