---
name: hic-tad-calling
description: This skill should be used when users need to identify topologically associating domains (TADs) from Hi-C data in .mcools files using HiCExplorer and cooltools. It provides workflows for TAD calling and visualization.
---

# TADs Calling with HiCExplorer and Cooltools

## Overview

This skill enables comprehensive identification and analysis of topologically associating domains (TADs) from Hi-C data stored in .mcools files. It integrates both HiCExplorer provide robust TAD calling and visualization capabilities. All the generated file should located in `${sample}_TAD_calling`, run `mkdir -p ${sample}_TAD_calling` first.

## Quick Start

To begin TADs analysis, first check for available .mcools files in the current directory:

```bash
find . -name "*.mcool" -o -name "*.cool"
```

Then proceed with the appropriate workflow based on your analysis goals.

## TAD Calling Workflow

### Step 1: Data Preparation

Verify that .mcools files are properly formatted and accessible. Check file integrity and resolution:

```bash
h5ls <input.mcool>resolutions
```

### Step 2: HiCExplorer TAD Calling

Use `hicFindTADs` for comprehensive TAD identification:

```bash
hicFindTADs --matrix <input.mcool::/resolutions/25000> --outPrefix <output_prefix> --minDepth 75000 --maxDepth 200000 --step 25000 --correctForMultipleTesting fdr --thresholdComparisons 0.05 --delta 0.01
```

Key parameters:
- `--minDepth`: Minimum window size (default: 75000)
- `--maxDepth`: Maximum window size (default: 200000)
- `--step`: Step size for sliding window (default: 25000)
- `--correctForMultipleTesting`: Multiple testing correction method
- `--thresholdComparisons`: FDR threshold for significant TADs


## Visualization

### Contact Maps with TAD Overlays

Ask the user for the target region, like `"chr1:1000000-5000000"`

Generate the `<track.ini>` file first for visualization, see `reference/template.ini` for example.

Generate contact maps with TAD boundaries:

```bash
hicPlotTADs --tracks <track.ini> --region chr1:1000000-5000000 --outFileName <outname>.pdf --dpi 300
```

## Resource

See `reference/template.ini` as an example for generating the `<track.ini>` file.
