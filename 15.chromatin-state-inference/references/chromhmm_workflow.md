# ChromHMM Detailed Workflow Guide

## Complete Analysis Pipeline

### Step 1: Software Installation and Setup

Install chromHMM and required dependencies:

```bash
# Download chromHMM
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip

# Set Java memory allocation
export JAVA_OPTS="-Xmx8g"
```

### Step 2: Input File Preparation

#### Cell Mark File Format
Create a tab-delimited file defining histone modifications:

```
H3K4me3	CellType1	/path/to/H3K4me3.bam
H3K27ac	CellType1	/path/to/H3K27ac.bam
H3K27me3	CellType1	/path/to/H3K27me3.bam
H3K36me3	CellType1	/path/to/H3K36me3.bam
```

#### Chromosome Sizes File
Create genome chromosome sizes file:

```
chr1	248956422
chr2	242193529
chr3	198295559
...
```

### Step 3: Data Binarization

Convert BAM files to binary presence/absence format:

```bash
java -jar ChromHMM.jar BinarizeBed \
  -b 200 \
  chromSizes.txt \
  bam_files/ \
  cellmarkfile.txt \
  binarized_data/
```

**Output files:**
- `binarized_data/*_binary.txt`: Binarized data files
- `binarized_data/*.bed`: BED files with binarized regions

### Step 4: Model Training

Train chromatin state model with different state numbers:

```bash
# Train 15-state model
java -jar ChromHMM.jar LearnModel \
  -p 8 \
  binarized_data/ \
  model_15states/ \
  15 \
  hg38

# Train 8-state model for comparison
java -jar ChromHMM.jar LearnModel \
  -p 8 \
  binarized_data/ \
  model_8states/ \
  8 \
  hg38
```

**Key output files:**
- `emissions_*.txt`: Emission probabilities for each state
- `transitions_*.txt`: Transition probabilities between states
- `model_*.txt`: Complete model parameters

### Step 5: Genome Segmentation

Generate chromatin state assignments:

```bash
java -jar ChromHMM.jar MakeSegmentation \
  model_15states/model_15.txt \
  binarized_data/ \
  segmentation_15states/
```

**Output files:**
- `*_segments.bed`: Genome segmentation in BED format
- `*_dense.bed`: Dense segmentation format

### Step 6: State Annotation and Enrichment

#### Genomic Feature Enrichment
```bash
# Gene annotations
java -jar ChromHMM.jar OverlapEnrichment \
  segmentation_15states/*_segments.bed \
  genes.bed \
  enrichment_genes/

# Promoter regions
java -jar ChromHMM.jar OverlapEnrichment \
  segmentation_15states/*_segments.bed \
  promoters.bed \
  enrichment_promoters/

# Enhancer annotations
java -jar ChromHMM.jar OverlapEnrichment \
  segmentation_15states/*_segments.bed \
  enhancers.bed \
  enrichment_enhancers/
```

### Step 7: Visualization

#### Heatmap Generation
```bash
# TSS-centered heatmaps
java -jar ChromHMM.jar MakeHeatmap \
  segmentation_15states/*_segments.bed \
  tss_coordinates.bed \
  heatmaps_tss/

# Enhancer-centered heatmaps
java -jar ChromHMM.jar MakeHeatmap \
  segmentation_15states/*_segments.bed \
  enhancer_coordinates.bed \
  heatmaps_enhancers/
```

#### Genome Browser Tracks
```bash
java -jar ChromHMM.jar MakeBrowserFiles \
  segmentation_15states/*_segments.bed \
  hg38 \
  browser_tracks/
```

## Multi-sample Analysis

### Comparative Analysis Between Conditions

```bash
# Compare two models
java -jar ChromHMM.jar CompareModels \
  model_condition1/ \
  model_condition2/ \
  model_comparison/

# Differential enrichment
java -jar ChromHMM.jar DifferentialEnrichment \
  segmentation_condition1/*_segments.bed \
  segmentation_condition2/*_segments.bed \
  genes.bed \
  diff_enrichment/
```

## Parameter Optimization

### Determining Optimal State Number

Use Bayesian Information Criterion (BIC) to select optimal state number:

```bash
# Train models with different state numbers
for states in 8 10 12 15 18 20 25; do
  java -jar ChromHMM.jar LearnModel \
    -p 8 \
    binarized_data/ \
    model_${states}states/ \
    $states \
    hg38

  # Calculate BIC
  java -jar ChromHMM.jar BIC \
    model_${states}states/model_${states}.txt \
    binarized_data/ \
    bic_${states}.txt
done
```

### Bin Size Selection

Test different bin sizes:

```bash
for binsize in 100 200 500 1000; do
  java -jar ChromHMM.jar BinarizeBed \
    -b $binsize \
    chromSizes.txt \
    bam_files/ \
    cellmarkfile.txt \
    binarized_${binsize}bp/

  java -jar ChromHMM.jar LearnModel \
    -p 8 \
    binarized_${binsize}bp/ \
    model_${binsize}bp/ \
    15 \
    hg38
done
```

## Quality Control Metrics

### Model Convergence
- Check log-likelihood convergence in learning output
- Verify emission parameters stabilize across iterations
- Assess transition probability consistency

### Biological Validation
- Compare with known chromatin state annotations (ENCODE)
- Validate promoter states with gene expression data
- Check enhancer states with functional validation data
- Assess reproducibility across replicates

### Technical Validation
- Verify input data quality (ChIP-seq quality metrics)
- Check binarization coverage across genome
- Assess model performance on held-out data
- Compare with alternative segmentation methods

## Advanced Workflows

### Custom State Definitions
Define custom state models based on specific biological contexts:

```bash
# Create custom emission priors
java -jar ChromHMM.jar LearnModel \
  -i custom_emissions.txt \
  -p 8 \
  binarized_data/ \
  custom_model/ \
  12 \
  hg38
```

### Integration with Other Data Types
Combine chromHMM with:
- RNA-seq for expression correlation
- ATAC-seq for accessibility validation
- TF ChIP-seq for factor binding enrichment
- Hi-C for 3D chromatin organization
- DNA methylation for epigenetic state correlation

### Time-series Analysis
Track chromatin state dynamics:

```bash
# Analyze multiple timepoints
for timepoint in T0 T6 T12 T24 T48; do
  java -jar ChromHMM.jar LearnModel \
    -p 8 \
    binarized_${timepoint}/ \
    model_${timepoint}/ \
    15 \
    hg38
done

# Compare timepoint models
java -jar ChromHMM.jar CompareModels \
  model_T0/ \
  model_T48/ \
  time_comparison/
```