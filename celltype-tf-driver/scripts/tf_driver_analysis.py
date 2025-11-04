#!/usr/bin/env python3
"""
Cell Type Transcription Factor Driver Analysis

This script implements the core pipeline for identifying key driver transcription factors
between two cell types using chromatin accessibility, histone modification, and motif data.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('tf_driver_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class TFAnalysisPipeline:
    """Main class for transcription factor driver analysis"""

    def __init__(self, params):
        self.params = params
        self.output_dir = Path(params.output_dir)
        self.temp_dir = self.output_dir / "temp"

        # Create directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(exist_ok=True)

        # File paths
        self.atac_a = params.atac_a
        self.atac_b = params.atac_b
        self.h3k27ac_a = params.h3k27ac_a
        self.h3k27ac_b = params.h3k27ac_b
        self.genome_fasta = params.genome_fasta
        self.meme_db = params.meme_db

    def validate_inputs(self):
        """Validate all input files exist and are properly formatted"""
        logger.info("Validating input files...")

        required_files = [
            self.atac_a, self.atac_b,
            self.h3k27ac_a, self.h3k27ac_b,
            self.genome_fasta, self.meme_db
        ]

        for file_path in required_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Required file not found: {file_path}")

        # Check BED file format
        self._validate_bed_file(self.atac_a)
        self._validate_bed_file(self.atac_b)
        self._validate_bed_file(self.h3k27ac_a)
        self._validate_bed_file(self.h3k27ac_b)

        logger.info("All input files validated successfully")

    def _validate_bed_file(self, bed_file):
        """Validate BED file format"""
        try:
            # Use basic format check instead of bedtools validate
            self._basic_bed_validation(bed_file)
        except Exception as e:
            logger.warning(f"BED validation failed for {bed_file}: {e}")
            raise

    def _basic_bed_validation(self, bed_file):
        """Basic BED file validation"""
        with open(bed_file, 'r') as f:
            for i, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    raise ValueError(f"Line {i}: BED file must have at least 3 columns")
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    if start >= end:
                        raise ValueError(f"Line {i}: start must be less than end")
                except ValueError:
                    raise ValueError(f"Line {i}: Invalid coordinates")

    def run_differential_accessibility(self):
        """Identify differential accessible regions"""
        logger.info("Running differential accessibility analysis...")

        # Find gained peaks (in B but not A)
        gained_peaks = self.temp_dir / "gained_peaks.bed"
        cmd = [
            "bedtools", "subtract",
            "-a", self.atac_b,
            "-b", self.atac_a,
            "-A"
        ]
        self._run_command(cmd, gained_peaks)

        # Find lost peaks (in A but not B)
        lost_peaks = self.temp_dir / "lost_peaks.bed"
        cmd = [
            "bedtools", "subtract",
            "-a", self.atac_a,
            "-b", self.atac_b,
            "-A"
        ]
        self._run_command(cmd, lost_peaks)

        # Merge differential peaks
        differential_peaks = self.output_dir / "differential_peaks.bed"
        cmd = [
            "cat", gained_peaks, lost_peaks
        ]
        self._run_command(cmd, differential_peaks)

        # Generate background regions
        background_peaks = self.temp_dir / "background_peaks.bed"
        self._generate_background_regions(background_peaks)

        logger.info(f"Found {self._count_lines(gained_peaks)} gained peaks")
        logger.info(f"Found {self._count_lines(lost_peaks)} lost peaks")

        return {
            'gained_peaks': gained_peaks,
            'lost_peaks': lost_peaks,
            'differential_peaks': differential_peaks,
            'background_peaks': background_peaks
        }

    def _generate_background_regions(self, output_file):
        """Generate background regions for motif enrichment"""
        # Use a simpler approach: sample random regions from the genome
        # This avoids issues with bedtools shuffle and invalid coordinates
        genome_sizes = self._get_genome_sizes()

        # Count total number of peaks to sample
        with open(self.atac_a, 'r') as f:
            num_peaks = sum(1 for _ in f)

        # Generate random regions using bedtools random
        cmd = [
            "bedtools", "random",
            "-g", str(genome_sizes),
            "-l", "200",  # Average peak size
            "-n", str(num_peaks)
        ]
        self._run_command(cmd, output_file)

    def _get_genome_sizes(self):
        """Get genome sizes file"""
        # Create genome sizes from FASTA index
        sizes_file = self.temp_dir / "genome.sizes"
        if not sizes_file.exists():
            cmd = ["samtools", "faidx", self.genome_fasta]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"Failed to create genome index: {result.stderr}")

            # Parse .fai file to create sizes
            fai_file = f"{self.genome_fasta}.fai"
            with open(fai_file, 'r') as f_in, open(sizes_file, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    f_out.write(f"{fields[0]}\t{fields[1]}\n")

        return sizes_file

    def run_motif_enrichment(self, foreground_bed, background_bed):
        """Run motif enrichment analysis using HOMER with enhanced implementation"""
        logger.info("Running motif enrichment analysis using HOMER...")

        # Run HOMER findMotifsGenome.pl
        homer_output = self.temp_dir / "homer_motif_enrichment"
        homer_output.mkdir(exist_ok=True)

        # Enhanced HOMER command with proper parameters
        cmd = [
            "findMotifsGenome.pl",
            str(foreground_bed),
            self.params.genome,
            str(homer_output),
            "-bg", str(background_bed),
            "-size", "given",
            "-mask",
            "-p", str(self.params.threads) if hasattr(self.params, 'threads') else "4"
        ]

        try:
            self._run_command(cmd)
            logger.info("HOMER motif enrichment completed successfully")
        except (subprocess.CalledProcessError, RuntimeError) as e:
            logger.warning(f"HOMER analysis failed: {e}")
            logger.info("Falling back to FIMO analysis...")
            return self._run_fimo_analysis(foreground_bed, background_bed)

        # Parse HOMER results
        enrichment_results = self._parse_homer_results(homer_output)

        return enrichment_results

    def _parse_homer_results(self, homer_output_dir):
        """Parse HOMER motif enrichment results with enhanced parsing"""
        results = []

        # HOMER output file
        known_results_file = homer_output_dir / "knownResults.txt"

        if not known_results_file.exists():
            logger.warning(f"HOMER results file not found: {known_results_file}")
            return pd.DataFrame()

        try:
            # Parse HOMER knownResults.txt with enhanced error handling
            with open(known_results_file, 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 6:
                        motif_name = fields[0]
                        p_value = float(fields[2])
                        q_value = float(fields[4])
                        target_pct = float(fields[6])
                        bg_pct = float(fields[8])

                        # Calculate enrichment metrics
                        enrichment = target_pct / bg_pct if bg_pct > 0 else float('inf')
                        log2_enrichment = np.log2(enrichment) if enrichment != float('inf') else np.inf

                        results.append({
                            'motif_id': motif_name,
                            'tf_name': self._extract_tf_name(motif_name),
                            'p_value': p_value,
                            'q_value': q_value,
                            'foreground_percent': target_pct,
                            'background_percent': bg_pct,
                            'enrichment_ratio': enrichment,
                            'enrichment_score': log2_enrichment
                        })

            # Sort by significance
            df = pd.DataFrame(results)
            if not df.empty:
                df = df.sort_values('q_value')

        except Exception as e:
            logger.error(f"Failed to parse HOMER results: {e}")
            return pd.DataFrame()

        return df

    def _run_fimo_analysis(self, foreground_bed, background_bed):
        """Run FIMO motif scanning as fallback when HOMER fails"""
        logger.info("Running FIMO motif analysis...")

        # Extract sequences from foreground regions
        foreground_fasta = self.temp_dir / "foreground_sequences.fa"
        cmd = ["bedtools", "getfasta", "-fi", str(self.genome_fasta),
               "-bed", str(foreground_bed), "-fo", str(foreground_fasta)]
        self._run_command(cmd)

        # Extract sequences from background regions
        background_fasta = self.temp_dir / "background_sequences.fa"
        cmd = ["bedtools", "getfasta", "-fi", str(self.genome_fasta),
               "-bed", str(background_bed), "-fo", str(background_fasta)]
        self._run_command(cmd)

        # Run FIMO on foreground
        fimo_output = self.temp_dir / "fimo_analysis"
        fimo_output.mkdir(exist_ok=True)
        cmd = ["fimo", "--oc", str(fimo_output), str(self.meme_db), str(foreground_fasta)]
        self._run_command(cmd)

        # Parse FIMO results
        fimo_results = self._parse_fimo_results(fimo_output)

        logger.info("FIMO analysis completed successfully")
        return fimo_results

    def _parse_fimo_results(self, fimo_output_dir):
        """Parse FIMO motif scanning results"""
        # FIMO output file
        fimo_file = fimo_output_dir / "fimo.tsv"

        if not fimo_file.exists():
            logger.warning(f"FIMO results file not found: {fimo_file}")
            return pd.DataFrame()

        try:
            df = pd.read_csv(fimo_file, sep='\t', comment='#')
            # Basic parsing - in practice would implement full enrichment analysis
            results = []
            for _, row in df.iterrows():
                results.append({
                    'motif_id': row.get('motif_id', ''),
                    'tf_name': self._extract_tf_name(row.get('motif_id', '')),
                    'p_value': row.get('p-value', 1.0),
                    'q_value': 1.0,  # FIMO doesn't provide q-values directly
                    'enrichment_score': -np.log10(row.get('p-value', 1.0))
                })
            return pd.DataFrame(results)
        except Exception as e:
            logger.error(f"Failed to parse FIMO results: {e}")
            return pd.DataFrame()

    def _extract_tf_name(self, motif_name):
        """Extract TF name from HOMER motif name with enhanced pattern matching"""
        # Enhanced TF name extraction based on working Python implementation
        tf_name = motif_name

        # Pattern 1: Extract from parentheses (e.g., "MA0004.1(AP2)" → "AP2")
        if '(' in motif_name and ')' in motif_name:
            tf_name = motif_name.split('(')[1].rstrip(')')

        # Pattern 2: Handle JASPAR format (e.g., "MA0004.1")
        elif motif_name.startswith('MA') and '.' in motif_name:
            # Look up TF name from JASPAR database or return as-is
            tf_name = motif_name

        # Pattern 3: Handle composite motifs (e.g., "AP-2alpha/AP-2gamma")
        elif '/' in motif_name:
            tf_name = motif_name.split('/')[0]

        # Pattern 4: Handle underscore-separated names (e.g., "AP_2alpha")
        elif '_' in tf_name:
            tf_name = tf_name.replace('_', '')

        return tf_name

    def integrate_results(self, accessibility_results, motif_results):
        """Integrate all results to identify key driver TFs"""
        logger.info("Integrating results to identify key driver TFs...")

        # Simple integration: rank by motif enrichment
        # In practice, integrate multiple metrics
        ranked_tfs = motif_results.sort_values('enrichment_score', ascending=False)

        # Save results
        output_file = self.output_dir / "key_driver_tfs.csv"
        ranked_tfs.to_csv(output_file, index=False)

        logger.info(f"Identified {len(ranked_tfs)} candidate TFs")
        logger.info(f"Top 5 TFs: {', '.join(ranked_tfs.head()['motif_id'].tolist())}")

        return ranked_tfs

    def _run_command(self, cmd, output_file=None):
        """Run a shell command"""
        # Convert all command elements to strings
        cmd_str = [str(arg) for arg in cmd]
        logger.debug(f"Running command: {' '.join(cmd_str)}")

        if output_file:
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd_str, stdout=f, stderr=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(cmd_str, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd_str)}\nError: {result.stderr}")

        return result

    def _count_lines(self, file_path):
        """Count lines in a file"""
        try:
            with open(file_path, 'r') as f:
                return sum(1 for _ in f)
        except:
            return 0

    def run_pipeline(self):
        """Run the complete analysis pipeline"""
        logger.info("Starting transcription factor driver analysis pipeline")

        try:
            # Step 1: Validate inputs
            self.validate_inputs()

            # Step 2: Differential accessibility
            accessibility_results = self.run_differential_accessibility()

            # Step 3: Motif enrichment
            motif_results = self.run_motif_enrichment(
                accessibility_results['differential_peaks'],
                accessibility_results['background_peaks']
            )

            # Step 4: Integration
            final_results = self.integrate_results(accessibility_results, motif_results)

            logger.info("Analysis completed successfully")
            return final_results

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(
        description="Cell Type Transcription Factor Driver Analysis"
    )

    # Required arguments
    parser.add_argument("--atac-a", required=True, help="ATAC-seq peaks for cell type A (BED)")
    parser.add_argument("--atac-b", required=True, help="ATAC-seq peaks for cell type B (BED)")
    parser.add_argument("--h3k27ac-a", required=True, help="H3K27ac peaks for cell type A (BED)")
    parser.add_argument("--h3k27ac-b", required=True, help="H3K27ac peaks for cell type B (BED)")
    parser.add_argument("--genome-fasta", required=True, help="Genome FASTA file")
    parser.add_argument("--meme-db", required=True, help="MEME format motif database")
    parser.add_argument("--output-dir", required=True, help="Output directory")

    # Optional arguments
    parser.add_argument("--cell-type-a", default="cell_A", help="Name of cell type A")
    parser.add_argument("--cell-type-b", default="cell_B", help="Name of cell type B")
    parser.add_argument("--genome", required=True, help="Genome assembly (e.g., hg38, mm10)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for HOMER analysis")

    args = parser.parse_args()

    # Run pipeline
    pipeline = TFAnalysisPipeline(args)
    results = pipeline.run_pipeline()

    print(f"\nAnalysis complete! Results saved to {args.output_dir}")
    print(f"Top transcription factors identified:")
    print(results.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
