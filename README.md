# ChromSkills

ChromSkills is a curated library of domain-specific Claude Code Skills for **agentic, interpretable, and reproducible chromatin data analysis** inside a Docker environment.

ChromSkills translates **natural-language analysis intents** into structured, tool-guided workflows across common epigenomic assays, while enforcing assay-aware decision logic and reducing unstable, free-form command generation.

---
![image](https://github.com/BIsnake2001/ChromSkills/blob/master/img/ChromSKills_framework.png)

## Scope

### What ChromSkills is for

ChromSkills focuses on **lightweight, workstation-friendly downstream analysis** that is practical on a desktop or laptop and supports reproducible reporting.

### What ChromSkills does NOT do

ChromSkills **does not include raw read mapping** (e.g., FASTQ → BAM) and is **not intended as a full HPC pipeline replacement**.

Users should provide **aligned or processed inputs** (see below).

---

## Supported assays and starting inputs

| Assay / data type | Start from | Examples of supported analysis |
|------------------|-----------|--------------------------------|
| ChIP-seq / ATAC-seq | BAM | alignment-level QC, filtering, peak calling, track generation, replicate handling, differential binding/accessibility, motif discovery/enrichment, peak annotation |
| WGBS | per-CpG methylation calls (BED/BedGraph/TSV) | DMR/DMC, global/local methylation profiling, variability/heterogeneity, UMR/LMR/PMD detection, integration with chromatin features, DMR–DEG integration |
| Hi-C | .cool / .mcool | normalization, matrix QC, A/B compartments and shifts, TADs (including nested TADs), loops, loop annotation, differential TAD analysis, regulatory community analysis |

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

Each Skill is stored as `SkillName/SKILL.md` in `~/.claude/skills` directory
Each MCP tools is stored as `MCPName.py` in `~/MCPs` directory

Skills are:

- Markdown-based
- Human-readable
- Domain-aware
- Deterministic (tool-constrained)

When needed:

1. Claude reads Skill metadata
2. Selects relevant Skills
3. Loads full Skill content
4. Executes tools with context-aware parameters

---

## Adding new Skills

1. Create a new folder in `~/.claude/skills` directory:

```
MySkill/SKILL.md
```

2. Write instructions + decision logic in Markdown
3. Create requiered `MCPName.py` in `~/MCPs` directory
4. Add MCPs into Claude Code by executing:
```bash
claude mcp add <MCPName> -s user -- python /root/MCPs/MCPName.py
```
5. Make sure the <MCPName> invoked in the Skills is the same as those in the command line above
6. Restart Claude Code by executing:

```bash
claude /mcp
```

---

## Design philosophy

ChromSkills:

- separates reasoning from execution
- constrains tools via structured MCP interfaces
- improves reproducibility
- reduces hallucinated commands
- enables laptop-friendly epigenomic analysis

