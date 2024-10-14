# Can parameterizations reproduce the gravity waves momentum fluxes and drag simulated by a global high resolution model?

This document contains how to use or generate the data in the paper "Can parameterizations reproduce the gravity waves momentum fluxes and drag simulated by a global high resolution model?" including the following files: raw ICON data, parameterization schemes and a python code for generating figures. 

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)


## Introduction
... directory contains ICON meteorological data.
.... directory contains code for run parameterization offline using ICON meteorological data.
paper_figures.py contains code for generating figures used in the paper, based on data obtained from parameterization and ICON simulations.

## Features

- **Figure Generation**: Scripts to produce plots and figures from data.
- **Data Handling**: Code to preprocess and analyze data from parameterization and ICON models.
- **Documentation**: Includes scripts and documentation relevant to the paper's figures and data.


## Installation

To set up the project on your local machine, follow the steps below:

### 1. Clone the repository

Start by cloning this repository to your local machine. Open a terminal and run the following command:
```bash
git clone https://github.com//imantoghraei/python_codes_all_figures.git
cd python_codes_all_figures
```

### 2. Running offline parameterization schemes 

"offline/run/" directory contains the scripts that compile the programs, link to the input dataset and produce various outputs. The Makefile
certainly needs to be adapted to the computer. To launch the offline parameterization, launch: ./laun.sh

### 3. Generating figures
To generate figures, simply run
python paper_figures.py


## Requirements

You need Fortran to run the parameterization offline and Python to generate figuers from paper_fugures.py
