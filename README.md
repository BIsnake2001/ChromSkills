# ChromSkills README

ChromSkills provides domain-specific biological analysis capabilities powered by Claude skills inside a Docker-based environment.

This guide walks you through installing prerequisites, preparing your workspace, running the Docker container, and starting a ChromSkills analysis session.

---

## 1. Install HOMER Locally (Required)

ChromSkills depends on the **HOMER** toolkit. Install HOMER and the required species packages on your local machine:

Installation guide: http://homer.ucsd.edu/homer/introduction/install.html

> **Note:** Downloading HOMER reference packages may take a long time.

After installation, locate your HOMER data directory:

```bash
$(dirname $(which homer))/../data/
```

This directory is referred to as **`$path_local_homer_data`** and must be mounted into the Docker container.

---

## 2. Prepare a Working Directory

Create a local directory containing the files you want to analyze with ChromSkills:

```
/path/to/my_work_dir/
```

This directory will later be mounted into the Docker container.

---

## 3. Pull the ChromSkills Docker Image

```bash
docker pull yuxuan2001/chromskills
```

---

## 4. Start the Docker Container

Replace the placeholders with your actual paths:

```bash
docker run -it   -p 8080:8080   -v <path_local_homer_data>:/mnt/softwares/homer/data   -v <path_local_work_dir>:/root/ChromOmics   --name chromskills   yuxuan2001/chromskills   /bin/bash
```

---

## 5. Configure Claude MCP for ChromSkills

Inside the container:

```bash
cd ~/ChromOmics/project_dir
claude /mcp
```

When the MCP initialization finishes, exit Claude.

---

## 6. Start a Chat Session

Inside the container:

```bash
chat
```

Now you can use ChromSkills to run biological analyses directly through the chat interface.

---

## You're Ready!

You now have a working ChromSkills environment capable of performing HOMER‑powered biological analyses with Claude MCP skills.

