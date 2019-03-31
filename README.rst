.. |date| date::

*******************************
PPTC PDX Oncoprint Generation
*******************************

:authors: Jo Lynne Rokita, Alvin Farrel, Khushbu Patel
:contact: Jo Lynne Rokita (rokita@email.chop.edu)
:organization: CHOP
:status: In-process
:date: |date|

.. meta::
   :keywords: pdx, mouse, WES, RNA-Seq, Fusions, SNP array, TMB, 2019
   :description: code to create PPTC PDX oncoprints by histology using WES mutations, RNA Fusions, and Copy Number data

Introduction
============

Here, we provide scripts to enable reproducible generation of Manuscript Figure 2: oncoprints by PDX histology. This repo contains code for:

1. Create expanded oncoprint matrices used to reproduce Figure 2

Details
=======

- data-fetch-pptc-pdx-oncoprints.sh
- RUN-create-full-oncoprint.R
- install-packages.R		
- demog-color-function.R
- mutation-color-function.R
- reformat-fusion-as-matrix.R
- load-maf-maftools.R
- create-mut-sigs-matrix.R		
- create-mut-matrices.R
- create-CN-matrices.R			
- merge-mut-CN-matrices.R
- merge-CN-gistic-matrices-new.R
- create-complexheat-oncoprint-all.R	
- merge-mut-CN-fusion-matrices.R


Software Requirements
=====================

R 3.4.3

Pipeline
========

.. code-block:: bash

         # How to run:
         # Download github repository in your home directory (~/)
         # Make sure to not clone the repository inside any other repository in your home directory
         git clone https://github.com/marislab/create-pptc-pdx-oncoprints.git
         
         # Run script to fetch data files
         ./data-fetch-pptc-pdx-oncoprints.sh
         
         # Run script to generate oncoprints
         Rscript RUN-create-full-oncoprint.R

