---

```markdown
name: gene-id-conversion
description: This skill performs gene identifier conversion between Ensembl IDs, Gene Symbols, Entrez IDs, and Protein IDs (e.g., UniProt). It supports multiple species (human, mouse, zebrafish, etc.) and can use online databases such as Ensembl BioMart or local annotation tables for offline mapping.

---

# Gene Identifier Conversion and Mapping

## Overview

1. Prepare a list or file of gene identifiers to convert (e.g., Ensembl IDs, Gene Symbols).  
2. Specify the source and target ID types (e.g., Ensembl → Symbol, Symbol → UniProt).  
3. Use R (`biomaRt`) or Python (`mygene`, `gprofiler-official`) for conversion.  
4. Save the converted results as a table for downstream analysis.  

---

## Decision Tree

### 1. Input Preparation

Input can be:
- A **text file** with one gene ID per line  
- A **column** in a CSV/TSV file  
- A **Python list** or **R vector**

**Example input file (`genes.txt`):**
```

ENSG00000139618
ENSG00000157764
ENSG00000141510

```

---

### 2. Choose Conversion Method

You can use one of the following tools:

| Tool | Language | Strength |
|------|-----------|-----------|
| `biomaRt` | R | Highly reliable for Ensembl → Symbol / Entrez / UniProt |
| `mygene` | Python | Fast, flexible for Symbol ↔ Ensembl ↔ UniProt |
| `gprofiler-official` | Python | Supports ortholog mapping (cross-species) |
| Local GTF annotation | shell/R | Offline mapping possible |

---

## Example 1. Ensembl → Gene Symbol (R + biomaRt)

```R
# Load library
library(biomaRt)

# Select dataset (example: human)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Convert Ensembl IDs to gene symbols
genes <- read.table("genes.txt", header=FALSE)$V1
mapping <- getBM(
  attributes=c("ensembl_gene_id", "hgnc_symbol"),
  filters="ensembl_gene_id",
  values=genes,
  mart=mart
)

write.table(mapping, "ensembl_to_symbol.txt", sep="\t", row.names=FALSE, quote=FALSE)
````

**Output (`ensembl_to_symbol.txt`):**

| ensembl_gene_id | hgnc_symbol |
| --------------- | ----------- |
| ENSG00000139618 | BRCA2       |
| ENSG00000157764 | BRAF        |
| ENSG00000141510 | TP53        |

---

## Example 2. Gene Symbol → UniProt ID (Python + mygene)

```python
from mygene import MyGeneInfo
mg = MyGeneInfo()

# Read gene symbols
genes = [line.strip() for line in open("symbols.txt")]

# Query UniProt and Ensembl mappings
res = mg.querymany(genes, scopes='symbol', fields='uniprot,ensembl.gene', species='human')

import pandas as pd
df = pd.DataFrame(res)
df[['query', 'uniprot', 'ensembl']].to_csv('symbol_to_uniprot.csv', index=False)
```

**Output:**

| query | uniprot | ensembl         |
| ----- | ------- | --------------- |
| TP53  | P04637  | ENSG00000141510 |
| BRCA2 | P51587  | ENSG00000139618 |

---

## Example 3. Cross-Species Gene Mapping (Python + gprofiler-official)

```python
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)

res = gp.convert(
    organism='hsapiens',
    query=['TP53', 'BRCA2', 'MYC'],
    target='drerio'  # zebrafish
)

res.to_csv('human_to_zebrafish_homologs.csv', index=False)
```

**Output example:**

| incoming | converted | target_organism |
| -------- | --------- | --------------- |
| TP53     | tp53      | drerio          |
| BRCA2    | brca2     | drerio          |
| MYC      | myca      | drerio          |

---

## Example 4. Offline Mapping Using GTF

If network access is unavailable, extract mappings from a local GTF file:

```bash
awk '$3=="gene"{split($9,a,";"); for(i in a){if(a[i]~/gene_id/) gid=a[i]; if(a[i]~/gene_name/) gname=a[i]} print gid"\t"gname}' annotation.gtf | sed 's/"//g; s/gene_id //; s/gene_name //' > gtf_mapping.tsv
```

Then use `join` or `merge` to convert IDs:

```bash
join -1 1 -2 1 <(sort genes.txt) <(sort gtf_mapping.tsv) > converted.txt
```

---

## Output Files

| File                              | Description                         |
| --------------------------------- | ----------------------------------- |
| `ensembl_to_symbol.txt`           | Ensembl → Gene Symbol mapping       |
| `symbol_to_uniprot.csv`           | Gene Symbol → UniProt mapping       |
| `human_to_zebrafish_homologs.csv` | Cross-species homolog mapping       |
| `gtf_mapping.tsv`                 | Local annotation-derived ID mapping |

---

## Best Practices

* Always specify the correct **species** (e.g., `hsapiens`, `mmusculus`, `drerio`).
* For batch conversion, clean duplicates and handle missing values.
* For reproducibility, save database version or GTF release.
* Prefer `biomaRt` for Ensembl-based pipelines, `mygene` for REST flexibility.

---

## Example Workflow Summary

```bash
# Step 1: Ensembl → Symbol
Rscript convert_ensembl_to_symbol.R

# Step 2: Symbol → UniProt
python convert_symbol_to_uniprot.py

# Step 3: Cross-species mapping (optional)
python human_to_zebrafish.py
```

---

## Troubleshooting

| Issue                 | Cause                        | Solution                                 |
| --------------------- | ---------------------------- | ---------------------------------------- |
| “Dataset not found”   | Wrong species dataset        | Check `listDatasets(useMart("ensembl"))` |
| Empty mapping results | ID version numbers included  | Remove `.1`, `.2` suffixes               |
| Inconsistent results  | Different Ensembl releases   | Fix version (e.g., Ensembl 110)          |
| Network timeout       | Remote database inaccessible | Use local GTF-based mapping              |

---
