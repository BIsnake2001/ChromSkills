#!/usr/bin/env python3
"""
Input Validation and Error Handling for TF Driver Analysis

This script provides comprehensive validation of input files and parameters
before running the main analysis pipeline.
"""

import os
import sys
import argparse
import yaml
import pandas as pd
from pathlib import Path
import subprocess
import logging

class InputValidator:
    """Class for validating all input files and parameters"""

    def __init__(self, config_file=None, **kwargs):
        self.config = {}
        self.errors = []
        self.warnings = []

        if config_file:
            self.load_config(config_file)
        else:
            self.config.update(kwargs)

    def load_config(self, config_file):
        """Load configuration from YAML file"""
        try:
            with open(config_file, 'r') as f:
                self.config = yaml.safe_load(f)
        except Exception as e:
            self.errors.append(f"Failed to load config file: {e}")

    def validate_all(self):
        """Run all validation checks"""
        logger.info("Starting comprehensive input validation...")

        # Check required parameters
        self._validate_required_parameters()

        # Validate file existence
        self._validate_file_existence()

        # Validate file formats
        self._validate_bed_files()
        self._validate_fasta_file()
        self._validate_meme_database()

        # Validate optional files
        if self.config.get('rna_seq_counts'):
            self._validate_rna_seq_files()

        # Validate analysis parameters
        self._validate_analysis_parameters()

        # Report results
        self._report_validation_results()

        return len(self.errors) == 0

    def _validate_required_parameters(self):
        """Validate that all required parameters are present"""
        required_params = [
            'cell_type_a', 'cell_type_b', 'genome_assembly', 'species',
            'atac_peak_a', 'atac_peak_b', 'h3k27ac_peak_a', 'h3k27ac_peak_b',
            'genome_fasta', 'meme_database', 'output_directory'
        ]

        for param in required_params:
            if param not in self.config or not self.config[param]:
                self.errors.append(f"Missing required parameter: {param}")

    def _validate_file_existence(self):
        """Validate that all required files exist"""
        required_files = [
            'atac_peak_a', 'atac_peak_b',
            'h3k27ac_peak_a', 'h3k27ac_b',
            'genome_fasta', 'meme_database'
        ]

        input_dir = self.config.get('input_directory', '')

        for file_param in required_files:
            file_path = self.config[file_param]
            if input_dir:
                file_path = os.path.join(input_dir, file_path)

            if not os.path.exists(file_path):
                self.errors.append(f"File not found: {file_path}")
            elif os.path.getsize(file_path) == 0:
                self.errors.append(f"File is empty: {file_path}")

    def _validate_bed_files(self):
        """Validate BED file format and content"""
        bed_files = [
            ('atac_peak_a', 'ATAC-seq peaks for cell type A'),
            ('atac_peak_b', 'ATAC-seq peaks for cell type B'),
            ('h3k27ac_peak_a', 'H3K27ac peaks for cell type A'),
            ('h3k27ac_peak_b', 'H3K27ac peaks for cell type B')
        ]

        input_dir = self.config.get('input_directory', '')

        for file_param, description in bed_files:
            file_path = self.config[file_param]
            if input_dir:
                file_path = os.path.join(input_dir, file_path)

            try:
                # Basic format validation
                df = pd.read_csv(file_path, sep='\t', header=None, comment='#')

                # Check minimum columns
                if df.shape[1] < 3:
                    self.errors.append(f"{description} ({file_path}): BED file must have at least 3 columns")

                # Check coordinate validity
                if df.shape[1] >= 3:
                    invalid_coords = df[(df[1] >= df[2]) | (df[1] < 0) | (df[2] < 0)]
                    if len(invalid_coords) > 0:
                        self.errors.append(f"{description} ({file_path}): Contains invalid coordinates")

                # Check peak count
                if len(df) < self.config.get('min_peak_count', 1000):
                    self.warnings.append(f"{description} ({file_path}): Low peak count ({len(df)}), may affect statistical power")

                # Check peak sizes
                if df.shape[1] >= 3:
                    peak_sizes = df[2] - df[1]
                    invalid_sizes = peak_sizes[(peak_sizes < self.config.get('min_peak_size', 100)) |
                                              (peak_sizes > self.config.get('max_peak_size', 1000))]
                    if len(invalid_sizes) > 0:
                        self.warnings.append(f"{description} ({file_path}): Contains peaks outside recommended size range")

            except Exception as e:
                self.errors.append(f"{description} ({file_path}): Failed to parse BED file - {e}")

    def _validate_fasta_file(self):
        """Validate genome FASTA file"""
        file_path = self.config['genome_fasta']
        input_dir = self.config.get('input_directory', '')

        if input_dir:
            file_path = os.path.join(input_dir, file_path)

        try:
            # Check if FASTA index exists
            fai_file = file_path + '.fai'
            if not os.path.exists(fai_file):
                self.warnings.append(f"Genome FASTA index not found: {fai_file}. Creating index...")
                self._create_fasta_index(file_path)

            # Basic FASTA validation
            cmd = ["samtools", "faidx", file_path]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.errors.append(f"Invalid FASTA file: {file_path}")

        except Exception as e:
            self.errors.append(f"Failed to validate FASTA file {file_path}: {e}")

    def _validate_meme_database(self):
        """Validate MEME motif database"""
        file_path = self.config['meme_database']
        input_dir = self.config.get('input_directory', '')

        if input_dir:
            file_path = os.path.join(input_dir, file_path)

        try:
            # Basic MEME format check
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('MEME version'):
                    self.errors.append(f"Invalid MEME database format: {file_path}")

            # Check if database contains motifs
            cmd = ["meme-get-motif", "-db", file_path]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.warnings.append(f"MEME database may be empty or corrupted: {file_path}")

        except Exception as e:
            self.errors.append(f"Failed to validate MEME database {file_path}: {e}")

    def _validate_rna_seq_files(self):
        """Validate RNA-seq files if provided"""
        input_dir = self.config.get('input_directory', '')

        # Validate count matrix
        counts_file = self.config['rna_seq_counts']
        if input_dir:
            counts_file = os.path.join(input_dir, counts_file)

        if not os.path.exists(counts_file):
            self.errors.append(f"RNA-seq count matrix not found: {counts_file}")
        else:
            try:
                # Check if it's a valid count matrix
                df = pd.read_csv(counts_file, sep='\t', index_col=0)
                if df.empty:
                    self.errors.append(f"RNA-seq count matrix is empty: {counts_file}")
            except Exception as e:
                self.errors.append(f"Failed to parse RNA-seq count matrix {counts_file}: {e}")

        # Validate metadata if provided
        if self.config.get('rna_seq_metadata'):
            metadata_file = self.config['rna_seq_metadata']
            if input_dir:
                metadata_file = os.path.join(input_dir, metadata_file)

            if not os.path.exists(metadata_file):
                self.errors.append(f"RNA-seq metadata not found: {metadata_file}")

    def _validate_analysis_parameters(self):
        """Validate analysis parameter ranges"""
        # FDR threshold
        fdr = self.config.get('fdr_threshold', 0.05)
        if not (0 < fdr <= 1):
            self.errors.append(f"Invalid FDR threshold: {fdr}. Must be between 0 and 1.")

        # Fold change threshold
        fc = self.config.get('fold_change_threshold', 1.5)
        if fc <= 0:
            self.errors.append(f"Invalid fold change threshold: {fc}. Must be positive.")

        # Weights
        weights = [
            self.config.get('weight_motif_enrichment', 0.4),
            self.config.get('weight_accessibility_change', 0.3),
            self.config.get('weight_enhancer_activity', 0.2),
            self.config.get('weight_expression_change', 0.1)
        ]

        if abs(sum(weights) - 1.0) > 0.01:
            self.errors.append(f"Integration weights must sum to 1.0. Current sum: {sum(weights)}")

        # Performance parameters
        threads = self.config.get('num_threads', 4)
        if threads < 1:
            self.errors.append(f"Invalid number of threads: {threads}. Must be at least 1.")

        memory = self.config.get('memory_limit_gb', 16)
        if memory < 1:
            self.errors.append(f"Invalid memory limit: {memory} GB. Must be at least 1 GB.")

    def _create_fasta_index(self, fasta_file):
        """Create FASTA index if it doesn't exist"""
        try:
            cmd = ["samtools", "faidx", fasta_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"Created FASTA index for {fasta_file}")
            else:
                self.errors.append(f"Failed to create FASTA index: {result.stderr}")
        except Exception as e:
            self.errors.append(f"Failed to create FASTA index: {e}")

    def _report_validation_results(self):
        """Report validation results"""
        if self.errors:
            logger.error("VALIDATION FAILED - Errors found:")
            for error in self.errors:
                logger.error(f"  • {error}")
        else:
            logger.info("✓ All required validations passed")

        if self.warnings:
            logger.warning("Validation warnings:")
            for warning in self.warnings:
                logger.warning(f"  • {warning}")

        # Summary
        logger.info(f"Validation summary: {len(self.errors)} errors, {len(self.warnings)} warnings")

    def get_validation_report(self):
        """Get comprehensive validation report"""
        report = {
            'valid': len(self.errors) == 0,
            'errors': self.errors,
            'warnings': self.warnings,
            'summary': f"{len(self.errors)} errors, {len(self.warnings)} warnings"
        }
        return report

def main():
    parser = argparse.ArgumentParser(description="Validate input files for TF driver analysis")
    parser.add_argument("--config", required=True, help="Configuration YAML file")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    global logger
    logger = logging.getLogger(__name__)

    # Run validation
    validator = InputValidator(args.config)
    is_valid = validator.validate_all()

    # Exit with appropriate code
    sys.exit(0 if is_valid else 1)

if __name__ == "__main__":
    main()