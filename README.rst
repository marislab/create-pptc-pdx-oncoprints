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

Here, we provide scripts to enable reproducible generation of Manuscript Figure 2: oncoprints by PDX histology and Figure 3: co-oncoplots for diagnosis/relapse samples in neuroblastoma, BCP-ALL, T-ALL, and osteosarcoma. This repo contains code for:

1. Create expanded oncoprint matrices used to reproduces Figure 2 and 3

Details
=======

- RUN-create-full-oncoprint-revision.R    
- load-maf-maftools.R
- cluster-scores.R                        
- merge-mut-CN-fusion-matrices.R
- co-oncoplots.R                          
- merge-mut-CN-matrices.R
- create-CN-matrices.R                    
- mutation-color-function.R
- create-complexheat-oncoprint-revision.R 
- reformat-CN-matrix-to-df.R
- create-mut-matrices.R                   
- reformat-clinical.R
- create-mut-sigs-matrix.R                
- reformat-fusion-as-matrix.R
- demog-color-function.R                  
- reformat-fusion-for-maf-revision.R
- install-packages.R



Software Requirements
=====================

R 3.4.3

Pipeline
========

.. code-block:: bash

         # How to run:
         # Download github repository in your home directory (~/)
         git clone https://github.com/marislab/create-pptc-pdx-oncoprints.git
        
         #
         brew update && brew install llvm
         # Run script to create pie chart
         Rscript ~/create-pptc-pdx-oncoprints/R/RUN-create-full-oncoprint-revision.R 

