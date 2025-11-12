---
name: hic-compartments-calling
description: This skill performs PCA-based A/B compartments calling on Hi-C .mcool datasets using cooltools, with automatic genome detection, chromosome naming adjustment.
---

# Hi-C Compartments Calling

## Overview

This skill provides an automated workflow for compartments calling using cooltools on .mcool Hi-C data.

- Always ask the user for the genome assembly of the HiC data. **Never detect the genome assembly automatically**.
- Load the genome FASTA file from `$(dirname $(which homer))/../data/genomes/${genome}/genome.fa` according to the user input and detect the chromosome name.
- Allow users to assign the resolution to call compartments on (default is 100kb resolution)
- Change chromosome name (e.g., "chr1" vs "1") in the .mcool or .cool if it is not consistent with the chromosome name used in the .fa file with `cooler.rename_chroms` (see step2 in decision tree for details).
- Generate chromosome-arm view files for compartments calling **after changing the chromosome name used in the .mcool or .cool file**.
- Perform PCA-based compartment analysis (`cooltools eigs-cis`) and extract the first principal component (PC1).
- Generate compartment interaction saddle plots and bigWig outputs for visualization.
- All the generated file should located in `${sample}_compartments_calling`, run `mkdir -p ${sample}_compartments_calling` first.

---

## When to Use This Skill

Use this skill when:
- You want to identify A/B compartments from Hi-C .mcool files.
- You need compartment scores and tracks for visualization.
- You want to perform reproducible, normalized, and automated compartment calling with minimal manual adjustment.

---

## Quick Start

Example commands to harmonize chromosome naming

```python
import cooler
clr = cooler.Cooler('input.mcool::/resolutions/100000')

rename_dict = {chrom:f"chr{chrom}" for chrom in clr.chromnames if not chrom.startswith('chr')}
cooler.rename_chroms(clr, rename_dict)
```

Example commands to perform compartment calling at 100 kb resolution:

```bash
cooler dump --header -t bins $cool_file::resolutions/100000 | cut -f1-3 > bins.100000.tsv

cooltools genome gc bins.100000.tsv ${genone}.fa > gc.100000.tsv

# Run compartment calling
cooltools eigs-cis -o sample.eigs.100000   --view view_hg38.tsv   --phasing-track gc.100000.tsv   --n-eigs 1   --clr-weight-name weight --bigwig input.mcool::resolutions/100000

# Generate saddle plot
cooltools saddle --qrange 0.02 0.98 -o sample.saddle.cis.100000   --view view_hg38.tsv   --clr-weight-name weight input.mcool::resolutions/100000   sample.eigs.100000.cis.vecs.tsv   sample.expected.cis.100000.tsv
```

```python
# Generate saddle plot
saddle = np.load('outputs/test.saddle.cis.100000.saddledump.npz', allow_pickle=True)
plt.figure(figsize=(6,6))
norm = LogNorm(    vmin=10**(-1), vmax=10**1)
im = plt.imshow(
    saddle['saddledata'],
    cmap='RdBu_r',
    norm = norm
);
plt.xlabel("saddle category")
plt.ylabel("saddle category")
plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
plt.savefig("${sample}.saddle.cis.100000.pdf")
```
---

## Decision Tree

### Step 1: Detect and Load Genome FASTA File

Automatically locate the genome file used for sequence-based reference.

```bash
GENOME=hg38
GENOME_FA=$(dirname $(which homer))/../data/genomes/${GENOME}/genome.fa
```

Verify the file exists:
```bash
if [ ! -f "$GENOME_FA" ]; then
  echo "Genome FASTA not found: $GENOME_FA"
  exit 1
fi
```

---

### Step 2: Harmonize Chromosome Naming

Change the chromosome name in the .mcool if not consistent with the chromosome name in the .fa file.

**Python Example:**
```python
import cooler
clr = cooler.Cooler('input.mcool::/resolutions/100000')

rename_dict = {chrom.replace('chr',''):chrom for chrom in clr.chromnames if chrom.startswith('chr')}
cooler.rename_chroms(clr, rename_dict)
```

---

### Step 3: Create Chromosome-Arm View

Generate view file defining chromosome arms (based on centromere positions). Change chromosome name (e.g., "chr1" vs "1") in the .mcool or .cool if it is not consistent to those in the .fa file with `cooler.rename_chroms`.

**Python Example:**
```python
import bioframe, cooler

hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
view_hg38 = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

clr = cooler.Cooler('input.mcool::/resolutions/100000')
view_hg38 = view_hg38[view_hg38.chrom.isin(clr.chromnames)]
view_hg38.to_csv('view_hg38.tsv', sep='\t', header=False, index=False)
```

---

### Step 4: Perform PCA-Based Compartment Calling

Run PCA on the observed/expected matrix and extract PC1.

**Command:**
```bash
cooltools expected-cis \
  --view view_hg38.tsv \
  --clr-weight-name weight \
  input.mcool::/resolutions/100000 \
  -o ${sample}.expected.cis.100kb.tsv

cooltools eigs-cis   -o ${sample}.eigs.100000   --view view_hg38.tsv   --phasing-track gc.100000.tsv   --n-eigs 1   --clr-weight-name weight --bigwig input.mcool::resolutions/100000
```

---

### Step 5: Generate Saddle Plot

Visualize the compartment interactions

**Command:**
```bash
cooltools saddle   --qrange 0.02 0.98 -o ${sample}.saddle.cis.100000 --view view_hg38.tsv --clr-weight-name weight input.mcool::resolutions/100000 ${sample}.eigs.100000.cis.vecs.tsv ${sample}.expected.cis.100000.tsv
```

- See `scripts/compartment_visualization.py` for the example of visualization of the A/B compartments. 

```python
saddle = np.load('outputs/test.saddle.cis.100000.saddledump.npz', allow_pickle=True)
plt.figure(figsize=(6,6))
norm = LogNorm(    vmin=10**(-1), vmax=10**1)
im = plt.imshow(
    saddle['saddledata'],
    cmap='RdBu_r',
    norm = norm
);
plt.xlabel("saddle category")
plt.ylabel("saddle category")
plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
plt.savefig("${sample}.saddle.cis.100000.pdf")
```

---

## Output

- `${sample}.eigs.100000.cis.vecs.tsv`: PC1 compartment scores  
- `${sample}.eigs.100000.bigwig`: PC1 track for genome browser  
- `${sample}.saddle.cis.100000.pdf`: Saddle plot visualization 
- `${sample}.compartment_distribution`: Compartment_distribution
- `view_${genome}.tsv`: Chromosome-arm view definition  

---

## Best Practices

- Always check balancing column names (`clr.bins().columns`) before proceeding.
- Use consistent genome versions for both .mcool and .fa files.
- Validate compartment polarity using GC content or gene density as phasing tracks.
- Start with 100 kb resolution; adjust for dataset size and sequencing depth.

---

## Resources

- See `scripts/compartment_visualization.py` for the example of visualization of the A/B compartments. 
