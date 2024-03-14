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

## Software requirements

Before you start this practical session, please make sure you have the following software installed on your PCs:

* **`PAML`**: you will be using the latest `PAML` release ([at the time of writing, v4.10.7](https://github.com/abacus-gene/paml/releases/tag/4.10.7)), available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml). If you do not want to install the software from the source code, then follow (A). If you want to install `PAML` from the source code, then follow (B). If you have a Mac with M1/M2 chips, please follow (C):

  * Installation (A): if you have problems installing `PAML` from the source code or you do not have the tools required to compile the source code, then you can [download the pre-compiled binaries available from the latest release by following this link](https://github.com/abacus-gene/paml/releases/tag/4.10.7). Please choose the pre-compiled binaries you need according to your OS, download the corresponding compressed file, and save it in your preferred directory. Then, after decompressing the file, please give executable permissions, export the path to this binary file so you can execute it from a terminal, and you should be ready to go!
    > **Windows users**: we highly recommend you install the Windows Subsystem for Linux (i.e., WSL) on your PCs to properly follow this tutorial -- otherwise, you may experience problems with the Windows Command Prompt. Once you have the WSL installed, then you can download the binaries for Linux.
  * Installation (B): to install `PAML` from the source code, please follow the instructions given in the code snippet below:

    ```sh
    # Clone to the `PAML` GitHub repository to get the latest `PAML` version
    # You can go to "https://github.com/abacus-gene/paml" and manually clone
    # the repository or continue below from the command line
    git clone https://github.com/abacus-gene/paml
    # Change name of cloned directory to keep track of version
    mv paml paml4.10.7
    # Move to `src` directory and compile programs
    cd paml4.10.7/src
    make -f Makefile
    rm *o
    # Move the new executable files to the `bin` directory and give executable
    # permissions
    mkdir ../bin
    mv baseml basemlg chi2 codeml evolver infinitesites mcmctree pamp yn00 ../bin
    chmod 775 ../bin/*
    ```
  
    Now, you just need to export the path to the `bin` directory where you have saved the executable file. If you want to automatically export this path to your `./bashrc` or your `~/.bash_profile`, you can run the following commands **AFTER ADAPTING** the absolute paths written in the code snippet below to those in your filesystem:

    ```sh
    # Run from any location. Change `~/.bashrc` if you are 
    # using another file
    printf "\n# Export path to PAML\n" >> ~/.bashrc
    # Replace "/c/Users/Bioinfor_tools/" with the path
    # that leads to the location where you have saved the
    # `paml4.10.7` directory. Modify any other part of the
    # absolute path if you have made other changes to the 
    # name of the directory where you have downloaded `PAML`
    printf "export PATH=/c/usr/bioinfo_tools/paml4.10.7/bin:\$PATH\n" >>  ~/.bashrc
    # Now, source the `~/.bashrc` file (or the file you are 
    # using) to update the changes
    source ~/.bashrc
    ```

    Alternatively, you can edit this file using your preferred text editor (e.g., `vim`, `nano`, etc.).
    >> **Windows users**: we highly recommend you install the Windows Subsystem for Linux (i.e., WSL) on your PCs to properly follow this tutorial -- otherwise, you may experience problems with the Windows Command Prompt. Once you have the WSL installed, then download the source code and follow the instructions listed above.

* Installation (C) for M1/M2 chips (Mac OSX): you will need to download the `dev` branch on the PAML GitHub repository and compile the binaries from the `dev` source code. Please [follow this link](https://github.com/abacus-gene/paml/tree/dev) and click the green button [`<> Code`] to start the download. You will see that a compressed file called `paml-dev.zip` will start to download. Once you decompress this file, you can go to directory `src` and follo the instructions in (B) to compile the binaries from the source code.

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) as we will be using them throughout the practical. The packages we will be using should work with R versions that are either newer than or equal to v4.1.2. If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages we will be using
    # Note that it may take a while if you have not 
    # installed all these software before
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan', 'colorBlindness'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned.
    ```

* **`FigTree`**: you can use this graphical interface to display tree topologies with/without branch lengths and with/without additional labels. You can then decide what you want to be displayed by selecting the buttons and options that you require for that to happen. You can [download the latest pre-compiled binaries, `FigTree v1.4.4` at the time of writing, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2` at the time of writing, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

* **Visual Studio Code**: for best experience with this tutorial, we highly recommend you install Visual Studio Code and run the tutorial from this IDE to keep everything tidy, organised, and self-contained. You can download it from [this website](https://code.visualstudio.com/). If you are new to VSC, you can check their webinars to learn about its various features and how to make the most out of it. You may also want to install the following extensions:
  * Markdown PDF -- developer: yzane
  * markdownlint -- developer: David Anson
  * Spell Right -- developer: Bartosz Antosik
  * vscode-pdf -- developer: tomoki1207

## Are you ready?

If you have gone through the previous sections, have a clear understanding of the how this repository is organised, and have installed the required software... Then you are ready to go!

You can start the tutorial by jumping on the [`00_data_formatting` directory](00_data_formatting), choose the `README.md` file that aligns with your OS, and...

Happy timetree inference! :)

----

ⓒ Dr Sandra Álvarez-Carretero | [`@sabifo4`](https://github.com/sabifo4/)