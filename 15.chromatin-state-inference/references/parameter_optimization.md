# ChromHMM Parameter Optimization Guide

## Determining Optimal Number of States

### Bayesian Information Criterion (BIC)

Use BIC to select the optimal number of chromatin states:

```bash
# Calculate BIC for different state numbers
for states in 6 8 10 12 15 18 20 25; do
  java -jar ChromHMM.jar BIC \
    model_${states}states/model_${states}.txt \
    binarized_data/ \
    bic_${states}.txt
done

# Compare BIC values
echo "State Number\tBIC Value"
for states in 6 8 10 12 15 18 20 25; do
  bic_value=$(grep "BIC" bic_${states}.txt | awk '{print $2}')
  echo "${states}\t${bic_value}"
done
```

**Interpretation:**
- Lower BIC values indicate better model fit
- Choose the state number where BIC plateaus or shows an "elbow"
- Consider biological interpretability alongside statistical fit

### Model Selection Criteria

1. **Statistical Fit**: BIC, AIC, log-likelihood
2. **Biological Interpretability**: Clear state definitions
3. **Stability**: Consistent results across runs
4. **Reproducibility**: Similar states across replicates

## Bin Size Selection

### Impact of Bin Size

| Bin Size | Resolution | Memory Usage | Computation Time | Use Cases |
|----------|------------|--------------|------------------|-----------|
| 50bp | Very High | Very High | Very High | High-resolution studies |
| 100bp | High | High | High | Detailed chromatin analysis |
| 200bp | Standard | Moderate | Moderate | Most applications |
| 500bp | Low | Low | Low | Genome-wide screening |
| 1000bp | Very Low | Very Low | Very Low | Large-scale comparisons |

### Testing Different Bin Sizes

```bash
# Test multiple bin sizes
for binsize in 50 100 200 500 1000; do
  # Binarize data
  java -jar ChromHMM.jar BinarizeBed \
    -b $binsize \
    chromSizes.txt \
    bam_files/ \
    cellmarkfile.txt \
    binarized_${binsize}bp/

  # Train model
  java -jar ChromHMM.jar LearnModel \
    -p 8 \
    binarized_${binsize}bp/ \
    model_${binsize}bp/ \
    15 \
    hg38

  # Calculate BIC
  java -jar ChromHMM.jar BIC \
    model_${binsize}bp/model_15.txt \
    binarized_${binsize}bp/ \
    bic_${binsize}bp.txt
done
```

## Learning Algorithm Parameters

### Convergence Criteria

```bash
# More iterations for better convergence
java -jar ChromHMM.jar LearnModel \
  -i 200 \
  -p 8 \
  binarized_data/ \
  model_converged/ \
  15 \
  hg38

# Check convergence in output
# Look for: "Converged after X iterations"
# Or stable log-likelihood values
```

### Memory and Performance Optimization

```bash
# Increase Java memory
export JAVA_OPTS="-Xmx16g -Xms4g"

# Use more processors
java -jar ChromHMM.jar LearnModel \
  -p 16 \
  binarized_data/ \
  model_fast/ \
  15 \
  hg38

# Reduce memory by filtering chromosomes
java -jar ChromHMM.jar BinarizeBed \
  -b 200 \
  -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
  chromSizes.txt \
  bam_files/ \
  cellmarkfile.txt \
  binarized_autosomes/
```

## Input Data Quality Assessment

### BAM File Quality Metrics

Before running chromHMM, ensure:

1. **Read Quality**:
   - Mapping quality > 30
   - Remove duplicates
   - Proper pair alignment

2. **Coverage**:
   - Sufficient depth (>10 million reads per mark)
   - Even coverage across genome
   - No extreme GC bias

3. **Signal-to-Noise**:
   - High enrichment over input
   - Clean peak calls
   - Reproducible across replicates

### Binarization Quality Control

Check binarization output:

```bash
# Check binarization statistics
grep "Total" binarized_data/*_stats.txt

# Expected outputs:
# - High coverage regions should be binarized
# - Low coverage regions should be excluded
# - Consistent binarization across marks
```

## Model Quality Assessment

### Emission Parameter Analysis

Evaluate state emission parameters:

1. **Distinct States**: Each state should have unique mark combinations
2. **Biological Plausibility**: States should match known chromatin patterns
3. **Consistency**: Similar states across different runs

### Transition Probability Analysis

Check state transition patterns:

1. **Spatial Coherence**: Similar states should cluster
2. **Biological Transitions**: Logical state transitions (e.g., promoter → transcribed)
3. **Boundary Definition**: Clear transitions at domain boundaries

### Segmentation Quality

Assess genome segmentation:

1. **Segment Size Distribution**: Appropriate for chromatin domains
2. **State Proportions**: Reasonable distribution across states
3. **Genomic Feature Overlap**: Expected enrichment patterns

## Advanced Parameter Tuning

### Custom Emission Priors

Define custom emission priors for specific biological contexts:

```bash
# Create custom emission file
cat > custom_emissions.txt << EOF
# State1: Active Promoter
1	H3K4me3	0.9
1	H3K27ac	0.8
# State2: Enhancer
2	H3K27ac	0.7
2	H3K4me1	0.6
EOF

# Use custom emissions
java -jar ChromHMM.jar LearnModel \
  -i custom_emissions.txt \
  -p 8 \
  binarized_data/ \
  custom_model/ \
  8 \
  hg38
```

### Multi-resolution Analysis

Combine different bin sizes for comprehensive analysis:

```bash
# Coarse resolution for large domains
java -jar ChromHMM.jar LearnModel \
  -p 8 \
  binarized_1000bp/ \
  model_coarse/ \
  8 \
  hg38

# Fine resolution for detailed patterns
java -jar ChromHMM.jar LearnModel \
  -p 8 \
  binarized_100bp/ \
  model_fine/ \
  15 \
  hg38
```

## Troubleshooting Common Issues

### Memory Problems

**Symptoms:** Java out of memory errors, slow performance

**Solutions:**
- Reduce bin size (e.g., 200bp → 500bp)
- Filter to autosomes only
- Increase Java memory allocation
- Use fewer states

### Convergence Issues

**Symptoms:** Model doesn't converge, unstable parameters

**Solutions:**
- Increase iterations (`-i 200`)
- Check input data quality
- Verify binarization coverage
- Try different random seeds

### Uninterpretable States

**Symptoms:** States with mixed or unclear mark patterns

**Solutions:**
- Adjust number of states
- Check histone mark combinations
- Verify input data quality
- Compare with known chromatin states