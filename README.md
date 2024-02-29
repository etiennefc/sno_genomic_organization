# SnoRNA genomic organization pipeline

__Author__ : Etienne Fafard-Couture

__Email__ :  _<etienne.fafard-couture@usherbrooke.ca>_

## License
SnoRNA genomic organization pipeline. Copyright (C) 2021  Étienne Fafard-Couture

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details at http://www.gnu.org/licenses/.

## Software to install
Conda (Miniconda3) and Mamba need to be installed (https://docs.conda.io/en/latest/miniconda.html)

For Linux users :
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Answer `yes` to `Do you wish the installer to initialize Miniconda3?`

To install mamba:
```bash
conda install -n base -c conda-forge mamba
```

To create the Snakemake environment used to launch Snakemake, run the following. The `conda create` command can appear to be stuck on `Solving environment`. While we are actually arguably [never going to solve the environment](https://www.ipcc.ch/sr15/chapter/spm/), the command is probably not stuck. Just be patient. If it does not work, use `mamba create` instead.

```bash
exec bash
conda config --set auto_activate_base False
conda create --name snakemake -c bioconda -c conda-forge snakemake=7.18
```

Before running Snakemake, you have to initialize the environment
```bash
conda activate snakemake
```

This workflow was tested with the following versions:
```bash
snakemake version == 7.18.2
conda version == 4.12.0
mamba version == 0.15.3
```

## Run

Run firstly the tasks requiring the internet locally (i.e all downloads) from the workflow directory with :
```bash
snakemake all_downloads --profile ../profile_local/
```


To run the workflow locally simply run the following command in the Snakemake conda environment from the workflow directory.
```bash
snakemake --profile ../profile_local/
```


To look at your entire workflow in svg, use the following commands (combined or separate conditions respectively):
```bash
snakemake --rulegraph | dot -Tsvg | display
snakemake --dag | dot -Tsvg | display
```

