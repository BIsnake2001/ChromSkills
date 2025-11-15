---
name: hic-compartments-calling
description: This skill performs PCA-based A/B compartments calling on Hi-C .mcool or .cool datasets using cooltools.
---

# Hi-C Compartments Calling

## Overview

This skill provides an automated workflow for compartments calling using **cooltools** on .mcool or .cool Hi-C data.

Main steps include:
- Refer to the **Inputs & Outputs** section to verify required files and output structure.
- **Always prompt user** for genome assembly used.
- **Always prompt user** for resolution used to call compartments.
- Load the genome FASTA file from `$(dirname $(which homer))/../data/genomes/${genome}/genome.fa` based on user input and detect the chromosome name.
- **Modify chromosome names** in the .mcool or .cool file.
- Generate chromosome-arm view files for compartment calling **after changing the chromosome name**.
- Perform **PCA-based compartment analysis** (`cooltools eigs-cis`) and extract the first principal component (PC1).
- Generate compartment interaction saddle plots and BigWig outputs for visualization.

## When to use this skill

Use this skill when:
- You want to identify A/B compartments from Hi-C .mcool or .cool files.
- You need compartment scores and tracks for visualization.
- You want to perform reproducible, normalized, and automated compartment calling with minimal manual adjustment.

## Inputs & Outputs

### Inputs

- **File format:** .mcool or .cool (Hi-C data file) data.
- **Genome assembly:** Prompt the user for genome assembly used.
- **Resolution:** Prompt the user for resolution used to call compartments.

### Outputs

```bash
compartments_calling/
    compartments/
      ${sample}.eigs.${resolution}.cis.vecs.tsv    # PC1 compartment scores  
      ${sample}.eigs.${resolution}.bigwig
      ${sample}.expected.${resolution}.cis.csv
      ... # other saddle outputs
    plots/         # PC1 track for genome browser  
      ${sample}.saddle.cis.${resolution}.pdf      # Saddle plot visualization 
      ${sample}.chr1_compartments.pdf # Compartment Scores along Chromosome 1
      ${sample}.compartment_scores_distribution.pdf
    temp/
      view_${genome}.tsv # Chromosome-arm view definition
      bins.${res}.tsv
      ... # other temp files                  
```

## Decision Tree

### Step 1: Detect and Load Genome FASTA File

Automatically locate the genome file used for sequence-based reference.

```bash
GENOME=hg38 # provided by user
GENOME_FA=$(dirname $(which homer))/../data/genomes/${GENOME}/genome.fa
```

Verify the file exists:

```bash
if [ ! -f "$GENOME_FA" ]; then
  echo "Genome FASTA not found: $GENOME_FA"
  exit 1
fi
```

### Step 2: Modify the chromosome name in .mcool file

**Python Example:**
```python
import cooler
clr = cooler.Cooler(f'input.mcool::/resolutions/{resoluton}')
rename_dict = {chrom:f"chr{chrom}" for chrom in clr.chromnames if not chrom.startswith('chr')}
cooler.rename_chroms(clr, rename_dict)
```

### Step 3: Create Chromosome-Arm View

Generate a view file defining chromosome arms (based on centromere positions). Change chromosome name in the .mcool or .cool file if not consistent with those in the .fa file using `cooler.rename_chroms`.

**Python Example:**
```python
import bioframe, cooler

hg38_chromsizes = bioframe.fetch_chromsizes('hg38') # change the genome name according to user feedback
hg38_cens = bioframe.fetch_centromeres('hg38')
view_hg38 = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

clr = cooler.Cooler(f'input.mcool::/resolutions/{res}')
view_hg38 = view_hg38[view_hg38.chrom.isin(clr.chromnames)]
view_hg38.to_csv('view_hg38.tsv', sep='	', header=False, index=False)
```

### Step 4: Perform PCA-Based Compartment Calling

Run PCA on the observed/expected matrix and extract PC1.

**example command:**

```python
import cooler
clr = cooler.Cooler(f'input.mcool::/resolutions/{resoluton}')
bins_df = clr.bins()[:][['chrom', 'start', 'end']]
bins_df.to_csv(f'bins.{res}.tsv', sep='\t', index=False, header=True)
```

```bash
cooltools genome gc bins.<res>.tsv <genome>.fa > gc.<res>.tsv
cooltools expected-cis --view view_hg38.tsv --clr-weight-name weight input.mcool::/resolutions/<res> -o <sample>.expected.cis.<res>.tsv
cooltools eigs-cis -o <sample>.eigs.<res> --view view_hg38.tsv --phasing-track gc.<res>.tsv --n-eigs 1 --clr-weight-name weight --bigwig input.mcool::resolutions/<res>
```

### Step 5: Generate Saddle Plot

Visualize the compartment interactions.

**Command:**
```bash
cooltools saddle --qrange 0.02 0.98 -o <sample>.saddle.cis.<res> --view view_hg38.tsv --clr-weight-name weight input.mcool::resolutions/<res> <sample>.eigs.<res>.cis.vecs.tsv <sample>.expected.cis.<res>.tsv
```

- See `scripts/compartment_visualization.py` for the example of visualization of the A/B compartments. 

```python
saddle = np.load(f'{sample}.saddle.cis.{res}.saddledump.npz', allow_pickle=True)
plt.figure(figsize=(6,6))
norm = LogNorm(vmin=10**(-1), vmax=10**1)
im = plt.imshow(
    saddle['saddledata'],
    cmap='RdBu_r',
    norm = norm
);
plt.xlabel("saddle category")
plt.ylabel("saddle category")
plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
plt.savefig(f"{sample}.saddle.cis.{res}.pdf")
```

# Best Practices

- Always check balancing column names (`clr.bins().columns`) before proceeding.
- Use consistent genome versions for both `.mcool` and `.fa` files.
- Validate compartment polarity using GC content or gene density as phasing tracks.

# Resources

- See `scripts/compartment_visualization.py` for the example of visualization of the A/B compartments. 
