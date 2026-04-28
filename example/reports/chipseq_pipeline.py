#!/usr/bin/env python3
"""
End-to-end ChIP-seq analysis pipeline: H3K4me3 and H3K27me3 peak calling & track generation
============================================================================================

Reproduces the complete analysis workflow:
  1. BAM filtration (remove artifacts, duplicates, MT reads, blacklist)
  2. Pool BAM replicates
  3. Peak calling with MACS2 (narrow for H3K4me3, broad for H3K27me3)
  4. BigWig track generation for IGV visualization

All functions imported from the MCP tool modules with exact parameters used.

Usage:
    python3 chipseq_pipeline.py

Requirements:
    - samtools, bedtools, picard, macs2, bedGraphToBigWig in PATH
    - MCP tool modules available
"""

import os
import sys
sys.path.insert(0, "/root")
from pathlib import Path
import asyncio

# ---------------------------------------------------------------------------
# MCP tool imports
# ---------------------------------------------------------------------------
from MCPs.qc_tools import bam_artifacts
from MCPs.bw_tools import (
    pool_bams,
    calculate_scaling_factor,
    generate_chrom_sizes,
    bam_to_bigwig,
)
from MCPs.macs2_tools import run_macs2

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR = Path("/root/ChromOmics/check_cost/ChIPseq")
FILTERED_DIR = BASE_DIR / "all_bam_filtration" / "filtered_bam"
BLACKLIST = BASE_DIR / "fixed_blacklist.bed"
PEAK_DIR = BASE_DIR / "all_peak_calling" / "peaks"
TRACK_DIR = BASE_DIR / "tracks"
TEMP_DIR = BASE_DIR / "all_bam_filtration" / "temp"

# Genome parameters (mouse mm10, no chr prefix)
GENOME_SIZE = "mm"
GENOME = "mm10"

# Sample definitions
SAMPLES = {
    "H3K4me3": {
        "rep1": "wt_H3K4me3_rep1",
        "rep2": "wt_H3K4me3_rep2",
        "peak_type": "narrow",
    },
    "H3K27me3": {
        "rep1": "wt_H3K27me3_rep1",
        "rep2": "wt_H3K27me3_rep2",
        "peak_type": "broad",
    },
    "input": {
        "rep1": "wt_input_rep1",
        "rep2": "wt_input_rep2",
    },
}


async def main():
    # =======================================================================
    # Step 1: BAM Filtration
    # =======================================================================
    print("=" * 70)
    print("STEP 1: BAM Filtration (remove artifacts)")
    print("=" * 70)

    os.makedirs(FILTERED_DIR, exist_ok=True)
    os.makedirs(TEMP_DIR, exist_ok=True)

    for mark, info in SAMPLES.items():
        for rep_key in ["rep1", "rep2"]:
            if rep_key not in info:
                continue
            sample_name = info[rep_key]
            raw_bam = BASE_DIR / f"{sample_name}.bam"
            filtered_bam = FILTERED_DIR / f"{sample_name}.filtered.bam"

            if filtered_bam.exists():
                print(f"  [SKIP] {sample_name} already filtered")
                continue

            print(f"  [FILTERING] {sample_name}")
            await bam_artifacts(
                bam_file=str(raw_bam),
                output_bam=str(filtered_bam),
                temp_dir=str(TEMP_DIR),
                blacklist_bed=str(BLACKLIST),
            )

    print("\n  BAM filtration complete.\n")

    # =======================================================================
    # Step 2: Pool replicates
    # =======================================================================
    print("=" * 70)
    print("STEP 2: Pool BAM replicates")
    print("=" * 70)

    pooled_bams = {}

    for mark, info in SAMPLES.items():
        if mark == "input":
            pooled_name = "wt_input_pooled"
            bam_list = [
                str(FILTERED_DIR / f"{info['rep1']}.filtered.bam"),
                str(FILTERED_DIR / f"{info['rep2']}.filtered.bam"),
            ]
        else:
            pooled_name = f"wt_{mark}_pooled"
            bam_list = [
                str(FILTERED_DIR / f"{info['rep1']}.filtered.bam"),
                str(FILTERED_DIR / f"{info['rep2']}.filtered.bam"),
            ]

        pooled_bam = FILTERED_DIR / f"{pooled_name}.filtered.bam"
        pooled_bams[mark] = str(pooled_bam)

        if pooled_bam.exists():
            print(f"  [SKIP] {pooled_name} already pooled")
            continue

        print(f"  [POOLING] {pooled_name}")
        await pool_bams(bam_files=bam_list, output_bam=str(pooled_bam))

    print("\n  Pooling complete.\n")

    # =======================================================================
    # Step 3: Peak Calling with MACS2
    # =======================================================================
    print("=" * 70)
    print("STEP 3: Peak Calling")
    print("=" * 70)

    os.makedirs(PEAK_DIR, exist_ok=True)

    input_pooled = pooled_bams["input"]

    for mark in ["H3K4me3", "H3K27me3"]:
        info = SAMPLES[mark]
        treatment_bam = pooled_bams[mark]
        peak_name = f"wt_{mark}_pooled"
        is_broad = info["peak_type"] == "broad"

        output_peak = PEAK_DIR / f"{peak_name}_peaks.{'broadPeak' if is_broad else 'narrowPeak'}"

        if output_peak.exists():
            print(f"  [SKIP] {peak_name} peaks already called")
            continue

        print(f"  [CALLING] {peak_name} ({'broad' if is_broad else 'narrow'})")
        await run_macs2(
            treatment_file=treatment_bam,
            control_file=input_pooled,
            name=peak_name,
            out_dir=str(PEAK_DIR),
            genome_size=GENOME_SIZE,
            format="BAMPE",              # paired-end
            broad=is_broad,
            broad_cutoff=0.1 if is_broad else None,
            qvalue=0.05,
            nomodel=False,
            shift=0,
            extsize=200,
        )

    print("\n  Peak calling complete.\n")

    # =======================================================================
    # Step 4: Generate BigWig tracks for IGV
    # =======================================================================
    print("=" * 70)
    print("STEP 4: BigWig track generation")
    print("=" * 70)

    os.makedirs(TRACK_DIR, exist_ok=True)
    bw_temp = str(TRACK_DIR / "temp")

    # Generate chrom.sizes from BAM header
    chrom_sizes = str(TRACK_DIR / "chrom.sizes")
    first_pooled = pooled_bams["H3K4me3"]
    print(f"  Generating chrom.sizes from {first_pooled}")
    await generate_chrom_sizes(bam_file=first_pooled, output_path=chrom_sizes)

    # Collect all BAMs to convert (pooled + individual replicates)
    bams_to_convert = {}

    # Pooled BAMs
    for mark in ["H3K4me3", "H3K27me3", "input"]:
        bams_to_convert[f"wt_{mark}_pooled"] = pooled_bams[mark]

    # Individual replicates
    for mark, info in SAMPLES.items():
        for rep_key in ["rep1", "rep2"]:
            if rep_key not in info:
                continue
            sample_name = info[rep_key]
            bams_to_convert[sample_name] = str(
                FILTERED_DIR / f"{sample_name}.filtered.bam"
            )

    for track_name, bam_path in bams_to_convert.items():
        output_bw = TRACK_DIR / f"{track_name}.bw"

        if output_bw.exists():
            print(f"  [SKIP] {track_name}.bw already exists")
            continue

        # Calculate RPM scaling factor
        scale = await calculate_scaling_factor(bam_file=bam_path)
        print(f"  [TRACK] {track_name}  scale_factor={scale:.4f}")

        await  bam_to_bigwig(
            bam_file=bam_path,
            chrom_sizes=chrom_sizes,
            output_bw=str(output_bw),
            scale_factor=scale,
            shift_tn5=False,            # ChIP-seq: no Tn5 shifting
            temp_dir=bw_temp,
        )

    print("\n  Track generation complete.\n")

    # =======================================================================
    # Summary
    # =======================================================================
    print("=" * 70)
    print("PIPELINE COMPLETE — Summary")
    print("=" * 70)

    h3k4me3_peaks = PEAK_DIR / "wt_H3K4me3_pooled_peaks.narrowPeak"
    h3k27me3_peaks = PEAK_DIR / "wt_H3K27me3_pooled_peaks.broadPeak"

    for peak_file, label in [
        (h3k4me3_peaks, "H3K4me3 (narrowPeak)"),
        (h3k27me3_peaks, "H3K27me3 (broadPeak)"),
    ]:
        if peak_file.exists():
            with open(peak_file) as f:
                count = sum(1 for _ in f)
            print(f"  {label}: {count} peaks  →  {peak_file}")

    print(f"\n  BigWig tracks  →  {TRACK_DIR}/")
    for bw in sorted(TRACK_DIR.glob("*.bw")):
        size_mb = bw.stat().st_size / 1e6
        print(f"    {bw.name}  ({size_mb:.1f} MB)")

    print("\n  Parameter logs  →  all_peak_calling/logs/")
    print("\n  Load .bw files and .narrowPeak/.broadPeak into IGV for visualization.")
    print("=" * 70)


if __name__ == "__main__":
    asyncio.run(main())
