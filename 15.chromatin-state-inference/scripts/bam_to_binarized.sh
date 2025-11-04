#!/bin/bash

# ChromHMM BAM to Binarized Conversion Script
# This script automates the conversion of BAM files to chromHMM binarized format

set -e  # Exit on error

# Default parameters
BIN_SIZE=200
NUM_PROCESSORS=8
GENOME_ASSEMBLY="hg38"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -b|--bin-size)
      BIN_SIZE="$2"
      shift 2
      ;;
    -p|--processors)
      NUM_PROCESSORS="$2"
      shift 2
      ;;
    -g|--genome)
      GENOME_ASSEMBLY="$2"
      shift 2
      ;;
    -i|--input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    -o|--output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -c|--cell-mark-file)
      CELL_MARK_FILE="$2"
      shift 2
      ;;
    -s|--chrom-sizes)
      CHROM_SIZES="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [options]"
      echo ""
      echo "Options:"
      echo "  -b, --bin-size SIZE        Bin size in bp (default: 200)"
      echo "  -p, --processors NUM       Number of processors (default: 8)"
      echo "  -g, --genome ASSEMBLY      Genome assembly (default: hg38)"
      echo "  -i, --input-dir DIR        Input directory with BAM files"
      echo "  -o, --output-dir DIR       Output directory for binarized data"
      echo "  -c, --cell-mark-file FILE  Cell mark file"
      echo "  -s, --chrom-sizes FILE     Chromosome sizes file"
      echo "  -h, --help                 Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$CELL_MARK_FILE" || -z "$CHROM_SIZES" ]]; then
  echo "Error: Missing required arguments"
  echo "Please provide: --input-dir, --output-dir, --cell-mark-file, --chrom-sizes"
  exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if ChromHMM JAR exists
if [[ ! -f "ChromHMM.jar" ]]; then
  echo "Error: ChromHMM.jar not found in current directory"
  echo "Please download ChromHMM from http://compbio.mit.edu/ChromHMM/"
  exit 1
fi

# Validate input files
echo "Validating input files..."
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input directory not found: $INPUT_DIR"
  exit 1
fi

if [[ ! -f "$CELL_MARK_FILE" ]]; then
  echo "Error: Cell mark file not found: $CELL_MARK_FILE"
  exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "Error: Chromosome sizes file not found: $CHROM_SIZES"
  exit 1
fi

# Check BAM files referenced in cell mark file
echo "Checking BAM files in cell mark file..."
while IFS=$'\t' read -r mark cell_type bam_file; do
  if [[ ! -f "$bam_file" ]]; then
    echo "Error: BAM file not found: $bam_file"
    exit 1
  fi
  # Check if BAM file is indexed
  if [[ ! -f "${bam_file}.bai" ]]; then
    echo "Warning: BAM file not indexed: $bam_file"
    echo "Creating index..."
    samtools index "$bam_file"
  fi
done < "$CELL_MARK_FILE"

# Run binarization
echo "Starting binarization with bin size: $BIN_SIZE bp"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Cell mark file: $CELL_MARK_FILE"
echo "Chromosome sizes: $CHROM_SIZES"

java -jar ChromHMM.jar BinarizeBed \
  -b "$BIN_SIZE" \
  "$CHROM_SIZES" \
  "$INPUT_DIR" \
  "$CELL_MARK_FILE" \
  "$OUTPUT_DIR"

# Check if binarization was successful
if [[ $? -eq 0 ]]; then
  echo "Binarization completed successfully!"
  echo "Output files created in: $OUTPUT_DIR"

  # Display binarization statistics
  echo ""
  echo "Binarization Statistics:"
  for stats_file in "$OUTPUT_DIR"/*_stats.txt; do
    if [[ -f "$stats_file" ]]; then
      echo "File: $(basename "$stats_file")"
      grep -E "(Total|Covered)" "$stats_file" | head -5
      echo ""
    fi
  done
else
  echo "Error: Binarization failed"
  exit 1
fi

echo "Binarization workflow completed!"
echo "Next step: Run model learning with:"
echo "java -jar ChromHMM.jar LearnModel -p $NUM_PROCESSORS $OUTPUT_DIR model_output/ 15 $GENOME_ASSEMBLY"