# MCMCtree tutorial

## Introduction

In this repository, you will find step-by-step guidelines to run `MCMCtree` using the approximate likelihood.

For this tutorial, we assume that the user has already (1) collected their data and (2) inferred the corresponding gene/genome alignment/s, and (3) inferred the corresponding phylogeny/ies. We will focus on the following:

1. Getting the data ready (i.e., correct format to run `MCMCtree`).
2. Setting a prior for the rates using R.
3. Running `BASEML` to calculate the Hessian and the gradient so we can use the approximate likelihood implemented in `MCMCtree` for timetree inference.
4. Using the former to estimate the divergence times with `MCMCtree`.

> NOTE: Note that the first thing that one needs to do is to get familiar with their dataset/s before proceeding with timetree inference: how were the data collected? How were the alignments generated? How are the files going to be organised? In this tutorial, we are not addressing these questions as we are just focusing on the subsequent steps (but we are working on other tutorials to provide some guidance to help the user navigate the aforementioned questions about their dataset/s). For a summary on how to approach them, however, we suggest reading [Alvarez-Carretero & dos Reis, 2022](https://link.springer.com/chapter/10.1007/978-3-030-60181-2_13).
