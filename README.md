# ChromSkills README

ChromSkills provides domain-specific biological analysis capabilities powered by Claude skills inside a Docker-based environment.

This guide walks you through installing prerequisites, preparing your workspace, running the Docker container, and starting a ChromSkills analysis session.

---

## 1. Pull the ChromSkills Docker Image

```bash
docker pull yuxuan2001/chromskills
```

---

## 2. Start the Docker Container

```bash
docker run -it -p 8080:8080 --name chromskills yuxuan2001/chromskills /bin/bash
```

---

## 3. Configure HOMER (Required)

ChromSkills depends on the **HOMER** toolkit. Install the required species packages according to the [guideline](http://homer.ucsd.edu/homer/introduction/configure.html):

```bash
cd /mnt/softwares/homer
perl configureHomer.pl -list # check available species
perl configureHomer.pl -install hg38 # this is an example
```


## 4. Configure Claude and tools required for ChromSkills

Inside the container:

- First run `claude`, log into your account and configure your claude code.

- Then build required MCPs by running:

```bash
cd ~
./build_mcps.sh
```

---

## 5. Initiate a project

Inside the container:

```bash
cd ${path_to_project_dir}
claude /mcp
```

When the project initialization finishes, exit Claude.

---

## 6. Start a Chat Session

Inside the container:

```bash
chat
```

Now you can use ChromSkills to run biological analyses directly through the chat interface.

---

## You're Ready!

You now have a working ChromSkills environment capable of performing biological analyses with Claude skills.

## Contact us

If you have any questions or suggestions, feel free to reach out to us at 2211289@tongji.edu.cn

