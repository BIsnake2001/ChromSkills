# ChromSkills

ChromSkills is a curated library of domain-specific Claude Code Skills for **agentic, interpretable, and reproducible chromatin data analysis** inside a Docker environment.

ChromSkills translates **natural-language analysis intents** into structured, tool-guided workflows across common epigenomic assays, while enforcing assay-aware decision logic and reducing unstable, free-form command generation.

---
![image](https://github.com/BIsnake2001/ChromSkills/blob/master/img/ChromSKills_framework.png)

## Scope

### What ChromSkills is for

ChromSkills focuses on **lightweight, workstation-friendly downstream analysis** that is practical on a desktop or laptop and supports reproducible reporting.

Typical use cases include:

- ATAC-seq / ChIP-seq peak calling and QC
- Differential accessibility / binding analysis
- Motif discovery and enrichment
- Hi-C compartment, TAD, and loop analysis
- WGBS methylation profiling and DMR analysis
- Integrative multi-omics analysis (e.g. ATAC + WGBS + RNA)

### 🚫 What ChromSkills does *NOT* do

To avoid incorrect expectations:

- ❌ **No read mapping** (FASTQ → BAM)
- ❌ **No heavy HPC pipelines**

ChromSkills **starts from processed inputs**, such as BAM, cool/mcool, or methylation tables.

---

## 📊 Can I use ChromSkills with my data?

Use the table below to quickly determine whether your data can be analyzed **directly** with ChromSkills.

| Assay | Required starting input | Supported analyses |
|------|------------------------|--------------------|
| **ATAC-seq** | BAM | QC, filtering, peak calling, replicate handling, footprinting, differential accessibility, motif analysis, track generation |
| **ChIP-seq** | BAM | QC, peak calling, genomic annotation, motif discovery/enrichment, differential binding |
| **WGBS** | Per-CpG methylation table (BED / BedGraph / TSV) | Global/local methylation, DMR/DMC, methylation variability, UMR/LMR/PMD detection |
| **Hi-C** | `.cool` / `.mcool` | Matrix QC, normalization, compartments, compartment shifts, TADs (including nested), loops, differential TADs |
| **Multi-omics** | Any combination above | ATAC–WGBS correlation, DMR–DEG integration, regulatory feature association |

---

## Installation

### 1. Pull Docker image

```bash
docker pull yuxuan2001/chromskills
```

---

## 2. Run the Docker container

Run ChromSkills in interactive mode:

```bash
docker run -it -p 8080:8080 --name chromskills yuxuan2001/chromskills /bin/bash
```

### Recommended: mount a local project directory

```bash
docker run -it -p 8080:8080 \
  -v /path/to/your/project:/work \
  --name chromskills \
  yuxuan2001/chromskills /bin/bash
```

Inside the container, work under:

```
/work
```

---

## 3. Configure HOMER (required for motif & annotation tasks)

```bash
cd /mnt/softwares/homer
perl configureHomer.pl -list
perl configureHomer.pl -install hg38
```

Install other species as needed.

---

## 4. Configure Claude Code and build MCP tools

Inside the container:

```bash
cd ~
./build_mcps.sh
```

This builds the structured MCP tool interfaces required by ChromSkills.

---

## 5. Choose your LLM backend

ChromSkills works with **any Claude Code–compatible model backend**.

### Option A — Claude (default)

```bash
export ANTHROPIC_API_KEY=your_key_here
```

### Option B — DeepSeek

```bash
export ANTHROPIC_BASE_URL=https://api.deepseek.com/anthropic
export ANTHROPIC_AUTH_TOKEN=your_key_here
```

### Option C — MiniMax

```bash
export ANTHROPIC_BASE_URL=https://api.minimax.io/anthropic
export ANTHROPIC_AUTH_TOKEN=your_key_here
```

---

## 6. Initialize a project

```bash
cd ${path_to_project_dir}
claude /mcp
```

Wait until initialization completes, then exit.

---

## 7. Start a chat session

```bash
chat
```

You can now perform analyses through natural language.

---

## Example usage

### Input data
```bash
fixed_blacklist.bed
wt_input_rep1.bam
wt_input_rep2.bam   
wt_H3K27me3_rep2.bam  
wt_H3K4me3_rep2.bam
wt_H3K27me3_rep1.bam  
wt_H3K4me3_rep1.bam   
```

### Prompt
```
Identify H3K4me3 and H3K27me3 peaks, annotate their genomic features, evaluate the data quality, and generate genome-wide signal tracks for visualization in IGV with available skills.
```
### Output

[example output reports from ChromSkills](example/reports/Task1.md)
---

## Skill architecture

Each Skill is:

- Markdown-based
- Human-readable
- Decision-tree driven
- Tool-constrained via MCP

Directory structure:

```
~/.claude/skills/SkillName/SKILL.md
~/MCPs/ToolName.py
```

Execution flow:

1.Claude reads Skill metadata<br>
2.Matches Skills to user intent<br>
3.Loads full Skill logic on demand<br>
4.Executes tools with context-aware parameters<br>

---

## Writing your own Skills

1. Create a new Skill.md file:

```bash
~/.claude/skills/MySkill/SKILL.md
```

2. Write instructions + decision logic in Markdown
3. Implement required MCP tools:

```bash
~/MCPs/MyTool.py
```
4. Register the MCP tool:
```bash
claude mcp add MyTool -s user -- python /root/MCPs/MyTool.py
```
5. Restart Claude Code by executing:
```bash
claude /mcp
```
---

## Design principles

ChromSkills:

- separates reasoning from execution
- constrains tools via structured MCP interfaces
- improves reproducibility
- reduces hallucinated commands
- enables laptop-friendly epigenomic analysis

## Contributing

We welcome:

- new Skills
- improved decision trees
- additional assay support
- documentation improvements

Feel free to open issues or pull requests.

## Contact us

If you have any questions or suggestions, feel free to reach out to us at 2211289@tongji.edu.cn