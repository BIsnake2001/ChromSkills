#!/usr/bin/env python3
"""
Generate compartment visualizations from existing Hi-C analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

print("Loading compartment analysis data...")

# Load eigenvectors
eigs_df = pd.read_csv('K562.eigs.100000.cis.vecs.tsv', sep='\t')
print(f"Loaded eigenvectors for {len(eigs_df)} bins")

# Load expected values
expected_df = pd.read_csv('K562.expected.cis.100000.tsv', sep='\t')
print(f"Loaded expected values for {len(expected_df)} regions")

# Load saddle plot data
saddle_data = np.load('K562.saddle.cis.100000.saddledump.npz', allow_pickle=True)
print("Loaded saddle plot data")

# Load eigenvalues
eigenvalues = pd.read_csv('K562.eigs.100000.cis.lam.txt', sep='\t')
print(f"Loaded eigenvalues for {len(eigenvalues)} chromosome arms")

print("\nGenerating compartment visualizations...")

# 1. Plot compartment scores along chromosomes
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Plot PC1 scores for first few chromosomes
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4']
for i, chrom in enumerate(chromosomes):
    if i >= 4:
        break
    chrom_data = eigs_df[eigs_df['chrom'] == chrom]
    if len(chrom_data) > 0:
        ax = axes[i//2, i%2]
        positions = (chrom_data['start'] + chrom_data['end']) / 2 / 1e6  # Convert to Mb
        ax.plot(positions, chrom_data['E1'], linewidth=1)
        ax.set_title(f'{chrom} - PC1 Compartment Scores')
        ax.set_xlabel('Position (Mb)')
        ax.set_ylabel('PC1 Score')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('K562.compartment_scores.pdf', dpi=300, bbox_inches='tight')
print("✓ Saved compartment scores plot: K562.compartment_scores.pdf")

# 2. Plot saddle plot
plt.figure(figsize=(8, 6))
saddledata = saddle_data['saddledata']

# Create saddle plot
norm = LogNorm(vmin=10**(-1), vmax=10**1)
im = plt.imshow(
    saddledata,
    cmap='RdBu_r',
    norm=norm,
    aspect='equal'
)

plt.xlabel("Compartment Category (A→B)")
plt.ylabel("Compartment Category (A→B)")
plt.title("K562 - Compartment Interaction Saddle Plot (100kb)")
plt.colorbar(im, label='Observed/Expected', pad=0.025, shrink=0.8)

# Add compartment labels
n_bins = saddledata.shape[0]
plt.xticks([0, n_bins-1], ['A', 'B'])
plt.yticks([0, n_bins-1], ['A', 'B'])

plt.savefig('K562.saddle_plot.pdf', dpi=300, bbox_inches='tight')
print("✓ Saved saddle plot: K562.saddle_plot.pdf")

# 3. Distribution of compartment scores
plt.figure(figsize=(10, 6))

plt.subplot(1, 2, 1)
plt.hist(eigs_df['E1'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
plt.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Zero line')
plt.xlabel('PC1 Compartment Score')
plt.ylabel('Frequency')
plt.title('Distribution of Compartment Scores')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
# Calculate percentage of A vs B compartments
A_compartments = len(eigs_df[eigs_df['E1'] > 0])
B_compartments = len(eigs_df[eigs_df['E1'] < 0])
total = len(eigs_df)

labels = ['A Compartments', 'B Compartments']
sizes = [A_compartments, B_compartments]
colors = ['lightcoral', 'lightskyblue']

plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
plt.title('A/B Compartment Distribution')

plt.tight_layout()
plt.savefig('K562.compartment_distribution.pdf', dpi=300, bbox_inches='tight')
print("✓ Saved compartment distribution: K562.compartment_distribution.pdf")

# 4. Create summary statistics
print("\n=== COMPARTMENT ANALYSIS SUMMARY ===")
print(f"Total genomic bins analyzed: {len(eigs_df):,}")
print(f"A compartments (PC1 > 0): {A_compartments:,} ({A_compartments/total*100:.1f}%)")
print(f"B compartments (PC1 < 0): {B_compartments:,} ({B_compartments/total*100:.1f}%)")
# Calculate average eigenvalue (excluding empty values)
valid_eigenvalues = pd.to_numeric(eigenvalues['eigval1'], errors='coerce')
avg_eigenvalue = valid_eigenvalues.mean()
print(f"Average eigenvalue: {avg_eigenvalue:.4f}")
print(f"Mean PC1 score: {eigs_df['E1'].mean():.4f}")
print(f"PC1 score std: {eigs_df['E1'].std():.4f}")

# Calculate compartment sizes
chrom_sizes = eigs_df.groupby('chrom').agg({
    'E1': ['count', 'mean', 'std']
}).round(4)
chrom_sizes.columns = ['num_bins', 'mean_pc1', 'std_pc1']

print("\nTop chromosomes by compartment strength:")
print(chrom_sizes.nlargest(5, 'mean_pc1'))

# Save summary statistics
summary_stats = {
    'total_bins': len(eigs_df),
    'A_compartments': A_compartments,
    'B_compartments': B_compartments,
    'A_percentage': A_compartments/total*100,
    'B_percentage': B_compartments/total*100,
    'avg_eigenvalue': avg_eigenvalue,
    'mean_pc1': eigs_df['E1'].mean(),
    'std_pc1': eigs_df['E1'].std()
}

summary_df = pd.DataFrame([summary_stats])
summary_df.to_csv('K562.compartment_summary.tsv', sep='\t', index=False)
print("\n✓ Saved summary statistics: K562.compartment_summary.tsv")

print("\n=== VISUALIZATION COMPLETE ===")
print("Generated files:")
print("  - K562.compartment_scores.pdf: PC1 scores along chromosomes")
print("  - K562.saddle_plot.pdf: Compartment interaction matrix")
print("  - K562.compartment_distribution.pdf: A/B compartment distribution")
print("  - K562.compartment_summary.tsv: Summary statistics")

print("\nThe A/B compartment analysis is complete!")