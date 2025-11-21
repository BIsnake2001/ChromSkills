---
name: hic-matrix-qc
description: This skill performs standardized quality control (QC) on Hi-C contact matrices stored in .mcool or .cool format. It computes coverage and cis/trans ratios, distance-dependent contact decay (P(s) curves), coverage uniformity, and replicate correlation at a chosen resolution using cooler and cooltools. Use it to assess whether Hi-C data are of sufficient quality for downstream analyses such as TAD calling, loop detection, and compartment analysis.
---

# Hi-C Contact Matrix QC for .mcool Files

## Overview

This skill performs QC on Hi-C matrices stored in .cool or .mcool format files at a user-selected resolution.

Main steps include:

- Refer to the **Inputs & Outputs** section to check required inputs and set up the output directory structure.  
- **Always wait the user feedback** if required files are not available in the current working directory by asking  
  `"${files} not available, provide required files or skip and proceed ?"`
- Inspect the `.mcool` file to list available resolutions and confirm the analysis resolution with the user.
- Compute coverage and cis/trans ratios using `cooltools coverage`.
- Assess coverage uniformity across bins from coverage tables.
- Compute cis expected contact frequency and distance-dependent contact decay (P(s) curves) with `cooltools expected-cis`.
- Visualize contact decay and P(s) scaling curves using Python.
- If multiple Hi-C replicates are provided, compute pairwise correlation of balanced matrices at the chosen resolution.
- Summarize QC metrics and plots into a structured output directory.

---

## When to use this skill

Use the hic-matrix-qc skill when you need to evaluate the quality of Hi-C contact matrices that are already stored in .cool or mcool format.

---

## Inputs & Outputs

### Inputs (choose one)

- sample.mcool
- sample.cool
User must choose one of the resolutions listed with cooler ls.

Optional: Multiple Hi-C matrices for replicate QC
`rep1.mcool` `rep2.mcool` `rep3.mcool`

### Outputs

```bash
hic_qc/
  logs/
    hic_qc.log               # Commands, parameters, and software versions
  metrics/
    coverage.${RES}.tsv               # Per-bin cis/total coverage from cooltools coverage
    cis_trans_summary.${RES}.txt      # Summarized cis, total, trans counts, and ratios
    expected_cis.${RES}.tsv           # Expected cis contacts vs distance from expected-cis
    ps_scaling_summary.${RES}.txt     # Optional table with P(s) slope(s) in defined distance ranges
    replicate_correlation.${RES}.tsv  # Pairwise correlation coefficients between replicates
  plots/
    coverage_histogram.${RES}.pdf     # Coverage uniformity plot
    ps_curve.${RES}.pdf               # P(s) curve (contact probability vs distance)
    decay_curve.${RES}.pdf            # Contact decay curve (raw/normalized)
    replicate_correlation_heatmap.${RES}.pdf  # Correlation matrix heatmap (if multiple replicates)
  comparison/
    replicate_vectors_${RES}.npz      # (Optional) Stored vectors used for replicate correlations
```
---

## Decision Tree

### Step 1: Check required inputs and list resolutions
```bash
MCool=sample.mcool
cooler ls "$MCool" # List groups/resolutions inside .mcool
RES=25000  # chosen by user
COOL="${MCool}::/resolutions/${RES}"
```

### Step 2: Compute coverage and cis/trans ratio

- Quantify cis and total coverage and derive cis/trans ratio at the chosen resolution.
- If the cooler is unbalanced or has a different weight column name, ask the user for the correct weight name or whether to use raw counts (empty --clr_weight_name).

- **example code**

```bash
# Compute coverage (cis and total) for each bin at the chosen resolution
cooltools coverage -o coverage.${RES}.tsv \
  --clr_weight_name weight \
   "$COOL"

# Summarize cis, total, and trans counts.
awk 'BEGIN{FS=OFS="\t"}
     NR>1 {
       cis += $NF-1;  # adjust columns according to actual header; verify first
       tot += $NF;    # example; inspect header before finalizing
     }
     END {
       trans = tot - cis;
       print "cis", cis;
       print "total", tot;
       print "trans", trans;
       print "cis_ratio", cis / tot;
     }' coverage.${RES}.tsv > cis_trans_summary.${RES}.txt
```
Output: `coverage.${RES}.tsv` `cis_trans_summary.${RES}.txt`

### Step 3: Assess coverage uniformity

- Draw the histogrm in python with `cis_trans_summary.${RES}.txt` into `coverage_histogram.${RES}.pdf`
- A reasonably broad distribution is expected; a long tail is common.
- Many zero-coverage bins may indicate insufficient depth at this resolution.
- A few bins with extremely high coverage may indicate local artifacts (e.g. centromeres, rDNA, mapping issues).

Output: `coverage_histogram.${RES}.pdf`

### Step 4: Compute cis expected and P(s) contact decay curve

- This computes expected cis interactions per diagonal (distance), averaged over all regions or view if specified.
- Plot the P(s) curve (log–log distance vs expected contacts).

- **example code**

```bash
cooltools expected-cis -o expected_cis.${RES}.tsv \
  --clr-weight-name weight \
  --ignore-diags 2 \
  "$COOL"
```
```python
exp_cis = pd.read_csv(f"expected_cis.{res}.tsv", sep="\t")
ps = exp_cis.groupby("dist")["balanced.avg"].mean().reset_index()
plt.loglog(ps["dist"], ps["balanced.avg"], marker=".")
```
Output: `ps_scaling_summary.${RES}.txt` `ps_curve.${RES}.pdf`

### Step 5: Contact decay curve from raw or balanced counts

- Visualize decay of interaction frequency vs distance, optionally comparing raw and balanced counts.
- This is often similar to the P(s) curve; you may reuse expected_cis.${RES}.tsv and plot:
  - dist vs n_valid (raw counts per distance)
  - dist vs balanced.avg (normalized P(s))

Output: `decay_curve.${RES}.pdf  `

### Step 6: Replicate correlation of Hi-C matrices (optional)

- Quantify similarity between Hi-C replicates at matrix level.
- Assumes:
  - At least two .mcool files (e.g. rep1.mcool, rep2.mcool, etc.).
  - Same genome assembly and resolution.

- **example code**

```python
mcools = {
    "rep1": f"rep1.mcool::/resolutions/{res}",
    "rep2": f"rep2.mcool::/resolutions/{res}",
    # add more replicates as needed
}
# load coolers
clrs = {name: cooler.Cooler(uri) for name, uri in mcools.items()}
# collect genome-wide vectors for each replicate
vectors = {name: [] for name in mcools.keys()}
for chrom in chroms:
    mats = {}
    for name, clr in clrs.items():
        mats[name] = clr.matrix(balance=True).fetch(chrom)

    # upper triangle indices (excluding diagonal)
    tri = np.triu_indices_from(next(iter(mats.values())), k=1)

    for name, mat in mats.items():
        v = mat[tri]
        # filter NaNs and infs later
        vectors[name].append(v)

# concatenate chromosome-wise vectors per replicate
for name in vectors:
    vectors[name] = np.concatenate(vectors[name])

# build correlation matrix
names = list(vectors.keys())
corr_mat = pd.DataFrame(index=names, columns=names, dtype=float)

for i, ni in enumerate(names):
    for j, nj in enumerate(names):
        v1 = vectors[ni]
        v2 = vectors[nj]
        mask = np.isfinite(v1) & np.isfinite(v2)
        r, _ = pearsonr(v1[mask], v2[mask])
        corr_mat.loc[ni, nj] = r

corr_mat.to_csv(f"replicate_correlation.{res}.tsv", sep="\t")

# plot heatmap
plt.figure()
im = plt.imshow(corr_mat.values.astype(float), vmin=0, vmax=1)
plt.xticks(range(len(names)), names, rotation=45, ha="right")
plt.yticks(range(len(names)), names)
plt.colorbar(im, label="Pearson correlation")
plt.title(f"Hi-C replicate correlation at {res}-bp resolution")
plt.tight_layout()
plt.savefig(f"replicate_correlation_heatmap.{res}.pdf")
```
Output: `replicate_correlation.{res}.tsv` `replicate_correlation_heatmap.{res}.pdf`

## Notes & troubleshooting

- If balancing weights are missing or correlation is calculated on raw counts, explicitly record this in logs and interpret with caution.
- Very low correlations (<0.7) between supposed biological replicates may indicate experimental issues or mismatched samples.