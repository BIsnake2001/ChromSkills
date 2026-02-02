# Integrative Analysis Report
## Identification of Genomic Regions with Increased Methylation, Decreased Accessibility, and Down-regulated Gene Expression

### Analysis Overview
This analysis integrates WGBS methylation data, ATAC-seq chromatin accessibility data, and RNA-seq differential expression results to identify genomic regions where:
1. DNA methylation increases (hyper-methylated DMRs)
2. Chromatin accessibility decreases (down-regulated DARs)
3. Target gene expression is down-regulated (DEGs)

All data are from chromosome 19 comparing K562 vs GM12878 conditions.

### Summary Statistics
| Analysis Step | Count |
|---------------|-------|
| Hyper-methylated DMRs (K562 vs GM12878) | 1,108 |
| Hypo-methylated DMRs (K562 vs GM12878) | 8,628 |
| Total significant DMRs (q<0.05, |Δ|>0.25) | 9,736 |
| Up-regulated DARs (K562 vs GM12878) | 1,400 |
| Down-regulated DARs (K562 vs GM12878) | 1,404 |
| Total significant DARs (padj<0.05, |log2FC|>1) | 3,449 |
| Intersected hyperDMR & downDAR regions | 492 |
| Unique genes annotated to intersected regions | 257 |
| Down-regulated DEGs (padj<0.05, log2FC<0) | 10,904 |
| Intersected genes (annotated & down-regulated) | 78 |
| Final integrated regions (peaks) | 158 |

### Analysis Pipeline Details

#### Step 1: Differential Methylation Analysis (Skill 21)

**Objective**: Identify differentially methylated regions (DMRs) between K562 and GM12878.

**Commands executed**:
1. Standardize WGBS BED files to 4-column Metilene format:
   ```bash
   awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $11/100}' WGBS_K562_rep1_chr19.bed | sort -V -k1,1 -k2,2n > DMR_DMC_detection/temp/WGBS_K562_rep1_chr19.sorted.bed
   # Repeated for WGBS_K562_rep2_chr19.bed, WGBS_GM12878_rep1_chr19.bed, WGBS_GM12878_rep2_chr19.bed
   ```

2. Generate merged methylation matrix for metilene:
   ```bash
   metilene_input.pl --in1 DMR_DMC_detection/temp/WGBS_K562_rep1_chr19.sorted.bed,DMR_DMC_detection/temp/WGBS_K562_rep2_chr19.sorted.bed --in2 DMR_DMC_detection/temp/WGBS_GM12878_rep1_chr19.sorted.bed,DMR_DMC_detection/temp/WGBS_GM12878_rep2_chr19.sorted.bed --out DMR_DMC_detection/temp/merged_input.bed --h1 K562 --h2 GM12878
   ```

3. Run metilene for DMR detection:
   ```bash
   metilene -a K562 -b GM12878 -t 1 -f 1 DMR_DMC_detection/temp/merged_input.bed > DMR_DMC_detection/stats/dmr_results.txt
   ```

4. Filter significant DMRs (q<0.05, |Δ|>0.25):
   ```bash
   # Internal filtering by metilene_tools
   # Output: significant_dmrs.tsv, significant_dmrs.bed
   ```

5. Extract hyper-methylated DMRs (Δ>0):
   ```bash
   awk -F'\t' '$5 > 0' DMR_DMC_detection/stats/significant_dmrs.tsv | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, "hyperDMR_"NR, $5, "."}' > DMR_DMC_detection/stats/hyper_dmrs.bed
   ```

**Output files**:
- `DMR_DMC_detection/temp/merged_input.bed` - Merged methylation matrix
- `DMR_DMC_detection/stats/dmr_results.txt` - Raw DMR results
- `DMR_DMC_detection/stats/significant_dmrs.tsv` - Significant DMRs (TSV)
- `DMR_DMC_detection/stats/significant_dmrs.bed` - Significant DMRs (BED)
- `DMR_DMC_detection/stats/hyper_dmrs.bed` - Hyper-methylated DMRs (1,108 regions)

---

#### Step 2: ATAC-seq Peak Calling (Skill 3)

**Objective**: Call peaks from ATAC-seq BAM files for each replicate.

**Commands executed**:
1. MACS2 peak calling for ATAC-seq (paired-end, nomodel, shift=-100, extsize=200):
   ```bash
   macs2 callpeak -t ATAC_K562_rep1_chr19.filtered.bam -f BAMPE -g hs --nomodel --shift -100 --extsize 200 -q 0.05 -n ATAC_K562_rep1 --outdir all_peak_calling/peaks
   # Repeated for ATAC_K562_rep2, ATAC_GM12878_rep1, ATAC_GM12878_rep2
   ```

**Output files**:
- `all_peak_calling/peaks/ATAC_K562_rep1_peaks.narrowPeak` - K562 replicate 1 peaks
- `all_peak_calling/peaks/ATAC_K562_rep2_peaks.narrowPeak` - K562 replicate 2 peaks
- `all_peak_calling/peaks/ATAC_GM12878_rep1_peaks.narrowPeak` - GM12878 replicate 1 peaks
- `all_peak_calling/peaks/ATAC_GM12878_rep2_peaks.narrowPeak` - GM12878 replicate 2 peaks

---

#### Step 3: Differential Accessibility Analysis (Skill 8)

**Objective**: Identify differentially accessible regions (DARs) between K562 and GM12878.

**Commands executed**:
1. Generate consensus peak set:
   ```bash
   # bedtools merge across all replicate peaks
   # Output: consensus_peaks.bed, consensus_peaks.saf
   ```

2. Generate read count matrix using featureCounts:
   ```bash
   featureCounts -a K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.saf -F SAF -o K562_vs_GM12878_DAR_analysis/tables/atac_counts.txt -T 4 -p -B -C ATAC_K562_rep1_chr19.filtered.bam ATAC_K562_rep2_chr19.filtered.bam ATAC_GM12878_rep1_chr19.filtered.bam ATAC_GM12878_rep2_chr19.filtered.bam
   ```

3. Prepare sample metadata:
   ```csv
   sample,condition,replicate
   ATAC_K562_rep1_chr19.filtered.bam,K562,1
   ATAC_K562_rep2_chr19.filtered.bam,K562,2
   ATAC_GM12878_rep1_chr19.filtered.bam,GM12878,1
   ATAC_GM12878_rep2_chr19.filtered.bam,GM12878,2
   ```

4. Run DESeq2 differential analysis:
   ```python
   # PyDESeq2 analysis: design = ~condition, contrast: K562 vs GM12878
   # Output: DAR_results.csv
   ```

5. Filter significant DARs (padj<0.05, |log2FC|>1):
   ```bash
   # Filtering and export to BED
   # Output: DAR_sig.bed, DAR_up.bed, DAR_down.bed
   ```

**Output files**:
- `K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.bed` - Consensus peak set
- `K562_vs_GM12878_DAR_analysis/tables/atac_counts.txt` - Read count matrix
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv` - DESeq2 results
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_down.bed` - Down-regulated DARs (1,404 regions)
- `K562_vs_GM12878_DAR_analysis/plots/Volcano.pdf` - Volcano plot
- `K562_vs_GM12878_DAR_analysis/plots/PCA.pdf` - PCA plot

---

#### Step 4: Intersect Hyper-methylated DMRs and Down-regulated DARs

**Objective**: Find regions with both increased methylation and decreased accessibility.

**Commands executed**:
```bash
bedtools intersect -a DMR_DMC_detection/stats/hyper_dmrs.bed -b K562_vs_GM12878_DAR_analysis/DARs/DAR_down.bed > integrated_analysis/hyperDMR_downDAR_intersect.bed
```

**Output files**:
- `integrated_analysis/hyperDMR_downDAR_intersect.bed` - 492 intersected regions

---

#### Step 5: Genomic Feature Annotation (Skill 10)

**Objective**: Annotate intersected regions to nearby genes using HOMER.

**Commands executed**:
1. Check genome installation:
   ```bash
   # Check if hg38 genome is installed in HOMER
   ```

2. Annotate peaks with HOMER:
   ```bash
   annotatePeaks.pl integrated_analysis/hyperDMR_downDAR_intersect.bed hg38 -annStats hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features_stats.txt -size given > hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features.txt
   ```

3. Visualize annotation statistics:
   ```python
   # Generate pie chart of genomic feature distribution
   # Output: hyperDMR_downDAR.anno_genomic_features_stats.pdf
   ```

**Output files**:
- `hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features.txt` - Annotated regions (full annotation)
- `hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features_stats.txt` - Annotation statistics
- `hyperDMR_downDAR_genomic_feature_annotation/plots/hyperDMR_downDAR.anno_genomic_features_stats.pdf` - Pie chart visualization

---

#### Step 6: Filter Down-regulated DEGs from RNA-seq

**Objective**: Identify significantly down-regulated genes from RNA-seq data.

**Commands executed**:
```bash
awk -F',' 'NR==1 {next} {if ($7+0 < 0.05 && $3+0 < 0) {gsub(/"/, "", $1); print $1, $3, $7}}' K562_vs_GM12878_deseq2.csv > DEG_analysis/down_deg_genes.tsv

# Strip version numbers from Ensembl IDs
awk '{split($1,a,"."); print a[1]}' DEG_analysis/down_deg_genes.tsv | sort -u > DEG_analysis/down_deg_genes_noversion.txt
```

**Output files**:
- `DEG_analysis/down_deg_genes.tsv` - Down-regulated DEGs with log2FC and padj (10,904 genes)
- `DEG_analysis/down_deg_genes_noversion.txt` - Gene IDs without version numbers

---

#### Step 7: Integrate with Down-regulated Genes

**Objective**: Identify which annotated genes are also down-regulated.

**Commands executed**:
1. Extract unique gene IDs from HOMER annotation:
   ```bash
   tail -n +2 hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features.txt | awk -F'\t' 'length($15) > 0 {print $15}' | sort -u > integrated_analysis/annotated_genes.txt
   ```

2. Find intersecting genes:
   ```bash
   comm -12 <(sort integrated_analysis/annotated_genes.txt) <(sort DEG_analysis/down_deg_genes_noversion.txt) > integrated_analysis/intersected_genes.txt
   ```

3. Create peak-gene mapping:
   ```bash
   tail -n +2 hyperDMR_downDAR_genomic_feature_annotation/results/hyperDMR_downDAR.anno_genomic_features.txt | awk -F'\t' 'length($15) > 0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $15 "\t" $16 "\t" $8 "\t" $10}' > integrated_analysis/peak_gene_mapping.tsv
   ```

4. Filter peaks for intersected genes:
   ```bash
   awk 'NR==FNR {genes[$1]=1; next} FNR==1 {print} $5 in genes' integrated_analysis/intersected_genes.txt integrated_analysis/peak_gene_mapping_header.tsv > integrated_analysis/final_integrated_regions.tsv
   ```

5. Create final BED file:
   ```bash
   awk 'NR>1 {print $2 "\t" $3 "\t" $4 "\t" $1 "\t1000\t."}' integrated_analysis/final_integrated_regions.tsv > integrated_analysis/final_integrated_regions.bed
   ```

**Output files**:
- `integrated_analysis/annotated_genes.txt` - 257 unique genes annotated to intersected regions
- `integrated_analysis/intersected_genes.txt` - 78 genes both annotated and down-regulated
- `integrated_analysis/peak_gene_mapping.tsv` - Peak to gene mapping
- `integrated_analysis/final_integrated_regions.tsv` - Final integrated regions with annotations
- `integrated_analysis/final_integrated_regions.bed` - Final integrated regions in BED format (158 regions)

---

### Key Findings

1. **78 genes** show coordinated epigenetic and transcriptional regulation:
   - Increased DNA methylation in their regulatory regions
   - Decreased chromatin accessibility
   - Down-regulated expression in K562 vs GM12878

2. **158 genomic regions** (peaks) are associated with these genes, representing candidate regulatory elements where epigenetic changes may drive gene expression changes.

3. The integration successfully identified regions meeting all three criteria: hyper-methylation, decreased accessibility, and target gene down-regulation.

### Biological Interpretation

#### 1. **Genomic Context of Identified Regions**
Based on HOMER annotation of the 492 intersected hyper-methylated/decreased accessibility regions:

| Genomic Feature | Number of Peaks | Enrichment (LogP) | Biological Significance |
|-----------------|-----------------|-------------------|-------------------------|
| **Promoter** | 48 | -20.297 (enriched) | Regions near transcription start sites; methylation here directly silences genes |
| **5' UTR** | 11 | -16.278 (enriched) | Untranslated regions near TSS; methylation can affect translation initiation |
| **CpG Islands** | 33 | -36.908 (enriched) | GC-rich regions normally protected from methylation; aberrant methylation here is a hallmark of epigenetic dysregulation |
| **Intron** | 168 | -34.606 (enriched) | Intronic regulatory elements, potential enhancers or silencers |
| **Exon** | 41 | -11.912 (enriched) | Coding regions; methylation can affect splicing |
| **Intergenic** | 28 | +10.364 (depleted) | Distal regulatory elements less represented |
| **LINE/SINE** | 117 | +13.496/+9.042 (depleted) | Repetitive elements show decreased representation |

**Interpretation**: The strong enrichment in promoters and CpG islands suggests these regions are **key regulatory elements** where DNA methylation changes likely directly impact gene expression. CpG island hyper-methylation is particularly significant as these regions are normally unmethylated in healthy cells. The distribution is visualized in `hyperDMR_downDAR_genomic_feature_annotation/plots/hyperDMR_downDAR.anno_genomic_features_stats.pdf`.

#### 2. **Potential Biological Mechanisms**
- **DNA methylation → Chromatin compaction → Reduced accessibility**: Increased CpG methylation recruits methyl-binding proteins (e.g., MeCP2) that promote heterochromatin formation, reducing ATAC-seq signal.
- **Transcription factor exclusion**: Methylation of CpG dinucleotides in transcription factor binding sites prevents TF binding, leading to decreased chromatin accessibility.
- **Coordinated epigenetic silencing**: The correlation between hyper-methylation and decreased accessibility suggests a **reinforcing epigenetic silencing loop** common in cell fate decisions and disease states.

#### 3. **Cell Type Context: K562 vs GM12878**
- **K562**: Chronic myelogenous leukemia (CML) erythroid precursor cells
- **GM12878**: Lymphoblastoid B-cell line (Epstein-Barr virus transformed)

The identified regions likely represent:
1. **Lymphoid-specific enhancers/promoters** that become silenced in erythroid lineage (K562)
2. **Epigenetic reprogramming** during leukemic transformation
3. **Cell identity genes** differentially regulated between hematopoietic lineages

#### 4. **Functional Implications**
- **Candidate silenced tumor suppressor genes**: In cancer (K562), promoter hyper-methylation often silences tumor suppressors.
- **Lineage-specific regulation**: Genes involved in B-cell function may be epigenetically silenced in erythroid cells.
- **Therapeutic targets**: Regions showing coordinated epigenetic changes could be targeted with epigenetic therapies (DNMT inhibitors, HDAC inhibitors).

#### 5. **Validation and Follow-up**
The 158 final integrated regions prioritize loci for:
- **ChIP-seq validation**: Check histone modifications (H3K4me3, H3K27ac) to confirm active/silenced states
- **Motif analysis**: Identify transcription factors whose binding is disrupted by methylation
- **CRISPR epigenetic editing**: Test causal relationships by targeted demethylation
- **Pathway analysis**: Determine biological pathways enriched among the 78 target genes

### Directory Structure

```
.
├── DMR_DMC_detection/          # Differential methylation results
├── all_peak_calling/           # ATAC-seq peak calling results
├── K562_vs_GM12878_DAR_analysis/ # Differential accessibility results
├── hyperDMR_downDAR_genomic_feature_annotation/ # HOMER annotation
├── DEG_analysis/               # Differential expression results
├── integrated_analysis/        # Integrated analysis results
│   ├── hyperDMR_downDAR_intersect.bed      # 492 intersected regions
│   ├── annotated_genes.txt                 # 257 annotated genes
│   ├── intersected_genes.txt               # 78 intersected genes
│   ├── final_integrated_regions.tsv        # Final regions with annotations
│   └── final_integrated_regions.bed        # Final regions in BED format
└── integration_analysis_report.md          # This report
```

### Conclusion
The analysis successfully identified 158 genomic regions on chromosome 19 where DNA methylation increases, chromatin accessibility decreases, and target gene expression is down-regulated in K562 compared to GM12878. These regions represent strong candidates for epigenetic regulation of gene expression and can be prioritized for further functional validation.

**Key biological insights**:
1. **Regulatory element identification**: Promoter and CpG island enrichment suggests these are key transcriptional control regions.
2. **Epigenetic silencing mechanism**: Coordinated hyper-methylation and chromatin closing indicates a reinforcing silencing loop.
3. **Cell type specificity**: Likely represent lymphoid-specific regulatory elements silenced in erythroid leukemia cells.
4. **Disease relevance**: Aberrant CpG island methylation is a hallmark of cancer epigenetics, suggesting potential tumor suppressor silencing.

**Limitations and future directions**:
- Analysis limited to chromosome 19; genome-wide analysis would provide complete picture
- Correlation does not prove causality; functional validation needed
- Time-course or treatment data could establish directionality of changes

**Recommended next steps**:
1. **Pathway analysis**: Use Skill 11 (functional enrichment) on the 78 target genes to identify affected biological processes
2. **Motif discovery**: Analyze sequences for enriched transcription factor binding sites
3. **Validation experiments**: ChIP-seq for histone modifications, CRISPR epigenetic editing
4. **Integration with additional data**: Hi-C for 3D chromatin structure, additional histone marks

The integrated approach demonstrates how multi-omics data can pinpoint regulatory elements where epigenetic changes correlate with transcriptional output, providing testable hypotheses for mechanistic studies.