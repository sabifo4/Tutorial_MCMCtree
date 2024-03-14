# Timetree inference with PAML: a step-by-step tutorial

**DISCLAIMER**: This tutorial is based on a phylogenetics tool that I am working on at the moment, which I am still developing and is yet to be published. While some of the scripts/tools that you will find here have been validated and used in published research ([Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1)), I am actively implementing new features as part of the current workflow of this pipeline as well as developing new scripts/tools. In other words, the code is not stable and I am still validating the new features. If you want to use the tools that I have developed as part of this tutorial/pipeline, please first <a href="mailto:sandra.ac93@gmail.com"><b>contact me</b></a>. Thank you :)

## Introduction

In this repository, you will find step-by-step guidelines to run `MCMCtree` using the approximate likelihood.

For this tutorial, we assume that the user has already (1) collected their data, (2) inferred the corresponding gene/genome alignment/s, and (3) inferred the corresponding phylogeny/ies. We will focus on the following:

1. Getting the data ready (i.e., correct format to run `MCMCtree`).
2. Setting a prior for the rates using R.
3. Running `BASEML` to calculate the Hessian and the gradient so we can use the approximate likelihood implemented in `MCMCtree` for timetree inference.
4. Using the former to estimate the divergence times with `MCMCtree`.

Specific `README.md` files and scripts have been generated so that users can follow this tutorial regardless of their operating system:

* `00_data_formatting`
  * Linux and WSL users can follow the instructions detailed in the [`REAMDE.md` file](00_data_formatting/README.md).
  * Mac OSX users can follow the instruction detailed in the [`REAMDE_MacOSX.md` file](00_data_formatting/README_MacOSX.md).
* `01_PAML/00_BASEML`
  * Linux and WSL users can follow the instructions detailed in the [`README.md` file](01_PAML/00_BASEML/README.md).
  * Mac OSX users can follow the instructions detailed in the [`README_MacOSX.md` file](01_PAML/00_BASEML/README_MacOSX.md).
  * Users that want to submit scripts to a High-Performance Computing server (SGE scheduler) can follow the instructions detailed in the [`README_HPC_SGE.md` file](01_PAML/00_BASEML/README_HPC_SGE.md).
* `01_PAML/01_MCMCtree`
  * Linux and WSL users can follow the instructions detailed in the [`README.md` file](01_PAML/01_MCMCtree/README.md).
  * Mac OSX users can follow the instructions detailed in the [`README_MacOSX.md` file](01_PAML/01_MCMCtree/README_MacOSX.md).
  * Users that want to submit scripts to a High-Performance Computing server (SGE scheduler) can follow the instructions detailed in the [`README_HPC_SGE.md` file](01_PAML/01_MCMCtree/README_HPC_SGE.md).

> NOTE: Everyone needs to be familiar with their dataset/s before proceeding with timetree inference: how were the data collected? How were the alignments generated? How are the files going to be organised? In this tutorial, we are not addressing these questions as we are just focusing on the subsequent steps. For a summary on how to approach them, however, we suggest reading [Álvarez-Carretero & dos Reis, 2022](https://link.springer.com/chapter/10.1007/978-3-030-60181-2_13).
