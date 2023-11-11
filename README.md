# nextflow-scripts

This repository contains various nextflow workflows.

## Orientation

Each subdirectory contains nextflow scripts, parameter templates, and other files for running the nextflow process. Generally, you run the process by executing:

```shell
nextflow run <script.nf> -params-file <JSON param file> -c nextflow.config
```

See comments in the nextflow files which describe the parameters contained in the template parameter JSON file.

**docker:/**

This directory contains general `Dockerfile`s required to execute the workflows. They will generally be named like `Dockerfile.<foo>` where `<foo>` is some descriptor. Note that these should relatively
general in nature and not include process-specific scripts, files, etc. If you require this, you can extend other containers or add process-specific `Dockerfile`s in the appropriate subdirectory for your
process.

**star_create_index/:**

This directory contains nextflow scripts and files for creating the STAR index
