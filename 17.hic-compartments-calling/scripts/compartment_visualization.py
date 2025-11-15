#!/usr/bin/env python3
"""
Visualize compartment calling results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns

# Parameters
sample_name = "K562"
resolution = 100000

# Load compartment scores
print("Loading compartment scores...")
eigs_file = f"compartments_calling/compartments/{sample_name}.eigs.{resolution}.cis.vecs.tsv"
eigs_df = pd.read_csv(eigs_file, sep='\t')

# Load saddle data
print("Loading saddle data...")
saddle_file = f"compartments_calling/plots/{sample_name}.saddle.cis.{resolution}.saddledump.npz"
saddle_data = np.load(saddle_file, allow_pickle=True)

# Plot 1: Compartment scores distribution
print("Plotting compartment scores distribution...")
plt.figure(figsize=(10, 6))
plt.hist(eigs_df['E1'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
plt.xlabel('PC1 (Compartment Score)')
plt.ylabel('Frequency')
plt.title(f'{sample_name} - Compartment Scores Distribution at {resolution}bp')
plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='Zero (A/B boundary)')
plt.legend()
plt.grid(alpha=0.3)
plt.savefig(f"compartments_calling/plots/{sample_name}.compartment_scores_distribution.pdf", bbox_inches='tight')
plt.close()

# Plot 2: Saddle plot
print("Plotting saddle plot...")
plt.figure(figsize=(8, 8))
norm = LogNorm(vmin=10**(-1), vmax=10**1)
im = plt.imshow(
    saddle_data['saddledata'],
    cmap='RdBu_r',
    norm=norm
)
plt.xlabel("Saddle Category")
plt.ylabel("Saddle Category")
plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
plt.title(f'{sample_name} - Saddle Plot at {resolution}bp')
plt.savefig(f"compartments_calling/plots/{sample_name}.saddle.cis.{resolution}.pdf", bbox_inches='tight')
plt.close()

# Plot 3: Compartment scores along chromosome 1 (example)
print("Plotting compartment scores along chromosome 1...")
chr1_data = eigs_df[eigs_df['chrom'] == 'chr1']
plt.figure(figsize=(12, 4))
plt.plot(chr1_data['start'], chr1_data['E1'], linewidth=1, color='blue')
plt.fill_between(chr1_data['start'], chr1_data['E1'], where=chr1_data['E1']>0,
                 alpha=0.3, color='red', label='A compartment (E1 > 0)')
plt.fill_between(chr1_data['start'], chr1_data['E1'], where=chr1_data['E1']<0,
                 alpha=0.3, color='blue', label='B compartment (E1 < 0)')
plt.xlabel('Genomic Position (bp)')
plt.ylabel('PC1 Score')
plt.title(f'{sample_name} - Compartment Scores along Chromosome 1')
plt.legend()
plt.grid(alpha=0.3)
plt.savefig(f"compartments_calling/plots/{sample_name}.chr1_compartments.pdf", bbox_inches='tight')
plt.close()

print("Visualization completed!")
print(f"Generated plots:")
print(f"- Compartment scores distribution: compartments_calling/plots/{sample_name}.compartment_scores_distribution.pdf")
print(f"- Saddle plot: compartments_calling/plots/{sample_name}.saddle.cis.{resolution}.pdf")
print(f"- Chromosome 1 compartments: compartments_calling/plots/{sample_name}.chr1_compartments.pdf")