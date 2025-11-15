---
name: hic-tad-calling
description: This skill should be used when users need to identify topologically associating domains (TADs) from Hi-C data in .mcools (or .cool) files or when users want to visualize the TAD in target genome loci. It provides workflows for TAD calling and visualization.
---

# TADs Calling with HiCExplorer and Cooltools

## Overview

This skill enables comprehensive identification and analysis of topologically associating domains (TADs) from Hi-C data stored in .mcool (or .cool) files. It integrates **HiCExplorer** for robust TAD calling and visualization capabilities.

Main steps include:

- Refer to the **Inputs & Outputs** section to verify required files and output structure.
- **Data Preparation**: Ensure .mcool files are formatted correctly and resolutions are verified.
- **Always prompt user** for resolution used to call TADs.
- **TAD Calling**: Use **HiCExplorer** to call TADs with customizable parameters.
- **Always prompt user** for target genomic loci for visualization.
- **Visualization**: Generate contact maps with TAD boundaries overlayed, for specific regions of the genome.

---

## When to use this skill

Use this skill when:

- You need to identify TADs in Hi-C data stored in .mcool (or .cool) files.
- You want to visualize TADs in a specific region of the genome.
- You need to perform automated TAD calling with HiCExplorer, including statistical corrections.

---

## Inputs & Outputs

### Inputs

- **File format:** .mcool or .cool (Hi-C data file).
- **Resolution:** Provided by user
- **Target region:** Genome region provided by user to visualize TADs (e.g., `"chr1:1000000-5000000"`).

### Outputs

```bash
TAD_calling/
    TADs/
        ${sample}_TAD_boundaries.bed  # Called TADs in BED format
        ${sample}_TAD_boundaries.gff
        ${sample}_TAD_domains.bed
        ... # other files output by the hicFindTADs
    plots/
        ${sample}_TADs_${genome_loci}.pdf  # TADs visualization (contact map)
    temp/
        track.ini            # Configuration file for visualization
```
---

## Decision Tree

### Step 1: Data Preparation

Verify that .mcool files are properly formatted and accessible. Check file integrity and resolution:

```bash
cooler ls <input.mcool>
```

### Step 2: HiCExplorer TAD Calling

Use `hicFindTADs` for comprehensive TAD identification. Customize parameters to suit the resolution and depth of your Hi-C data:

```bash
hicFindTADs --matrix input.mcool::/resolutions/<res> --outPrefix <output_prefix> --minDepth 75000 --maxDepth 200000 --step 25000 --correctForMultipleTesting fdr --thresholdComparisons 0.05 --delta 0.01
```

Key parameters:
- `--minDepth`: Minimum window size, should be at least 3 times as large as the bin size of the Hi-C matrix..
- `--maxDepth`: Maximum window size, should around 6-10 times as large as the bin size of the Hi-C matrix.
- `--step`: Step size for sliding window.
- `--correctForMultipleTesting`: Multiple testing correction method.
- `--thresholdComparisons`: FDR threshold for significant TADs.

## Step 3: Visualization
Contact Maps with TAD Overlays
- Ask the user for the target region, like `"chr1:1000000-5000000"`.
- Generate the `<track.ini>` file first for visualization, using the template provided in `reference/template.ini`.
- Generate contact maps with TAD boundaries.
- Adjust the parameters in the `<track.ini>` to make the plot more precise.

```bash
hicPlotTADs --tracks <track.ini> --region <genome_loci> --outFileName <outname>.pdf --dpi 300
```
---

## Resources

- See `reference/template.ini` for an example of how to generate the `<track.ini>` file.
