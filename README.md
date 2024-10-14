# Can Parameterizations Reproduce Gravity Waves Momentum Fluxes and Drag Simulated by a Global High-Resolution Model?

This repository contains the necessary files and instructions for generating data and figures used in the paper *"Can parameterizations reproduce the gravity waves momentum fluxes and drag simulated by a global high-resolution model?"*. It includes raw ICON data, parameterization schemes, and Python code for figure generation.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)

## Introduction

- The `ICON_data/` directory contains ICON meteorological data.
- The `offline/` directory contains code to run the parameterization schemes offline using ICON meteorological data.
- The `paper_figures.py` script generates the figures used in the paper, based on data from the parameterization and ICON simulations.

## Features

- **Figure Generation**: Python scripts to produce plots and figures from the processed data.
- **Data Handling**: Fortran routines to preprocess and analyze data from parameterization and ICON models.
- **Documentation**: Includes scripts and instructions for reproducing the paper's figures and data.

## Installation

Follow these steps to set up the project on your local machine:

### 1. Clone the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/imantoghraei/python_codes_all_figures.git
cd python_codes_all_figures
```

### 2. Running Offline Parameterization Schemes

1. **Download ICON Data**: The data is linked in the `ICON.txt` file. Download the data as instructed in that file.
2. **Offline Parameterization Setup**: The required files for running the parameterization offline are in the `offline/` directory.

- **`run/` subdirectory**: Contains scripts for compiling the programs, linking input datasets, and generating outputs. The `Makefile` may need to be adapted for your specific system. To launch the offline parameterization, run:

  ```bash
  ./loop_gwd_healpix_chopin_2d.sh
  ```

- **`prog/` subdirectory**: Contains the Fortran routines that execute the parameterization schemes:
  - `laun_gwd_healpix.f90`: Loads ICON data in NetCDF format and calculates drag and momentum fluxes.
  - `acama_gwd_rando_QBOi.f90`: Calculates fluxes due to fronts.
  - `flott_gwd_rando_QBOi.f90`: Calculates fluxes due to precipitation.
  - `orografi_strato_QBOi.f90`: Calculates fluxes due to orography.

### 3. Generating Figures

Once the parameterization is completed, you can generate the figures by running the following Python script:

```bash
python paper_figures.py
```

## Requirements

- **Fortran**: Required for running the offline parameterization schemes.
- **Python**: Required for generating figures using `paper_figures.py`. Install any necessary Python packages listed in the `requirements.txt` (if available).
