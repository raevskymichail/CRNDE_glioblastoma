# CRNDE in glioblastoma

<!-- <h1 align="center">
  CRNDE in glioblastoma
  <br>
</h1> -->

[![](https://img.shields.io/github/languages/code-size/raevskymichail/CRNDE_glioblastoma)](https://img.shields.io/github/languages/code-size/raevskymichail/CRNDE_glioblastoma)
[![](https://img.shields.io/github/languages/top/raevskymichail/CRNDE_glioblastoma)](https://img.shields.io/github/languages/top/raevskymichail/CRNDE_glioblastoma)
[![](https://img.shields.io/github/issues/raevskymichail/CRNDE_glioblastoma)](https://img.shields.io/github/issues/raevskymichail/CRNDE_glioblastoma)
[![](https://img.shields.io/github/license/raevskymichail/CRNDE_glioblastoma)](https://img.shields.io/github/license/raevskymichail/CRNDE_glioblastoma)


## Introduction

This repository contains primary source code for *"Overexpression of CRNDE in glioblastoma is a poor survival prognosis biomarker"* manuscript. 

Glioblastoma (GBM) is the most common and malignant brain malignancy worldwide with a 10-year survival of only *0.7%*. Aggressive multimodal treatment is not enough to increase life expectancy and provide good quality of life for glioblastoma patients. In addition, despite decades of research, there are no established biomarkers for early disease diagnosis and monitoring of patient response to treatment. High throughput sequencing technologies allow for identification of molecules from large clinically annotated datasets. Thus, the aim of our study was to identify significant molecular changes between short- and long-term glioblastoma survivors by transcriptome RNA sequencing profiling using previous data from the publicly available repositories The Cancer Genome Atlas (TCGA; number of annotated cases = 135) and Chinese Glioma Genome Atlas (CGGA; number of annotated cases = 218), and using experimental clinically annotated GBM tissue samples from the Institute of Pathology, Faculty of medicine in Ljubljana corresponding to 2-50 months overall survival (*n = 16*). We identified overlapping sets of congruently regulated differential genes involved in cell growth, division and migration, structure and dynamics of extracellular matrix, DNA methylation, and regulation through noncoding RNAs. We found one differential gene for long noncoding RNA CRNDE that showed the best potential to be used as a negative prognostic biomarker for glioblastoma and was significant in Kaplan-Meier overall survival analyses in both literature and experimental datasets investigated (*p=0.013 - 0.046*) and in Kaplan-Meier progression-free survival in literature dataset (*p=0.0016*) 

## 📝 Requirements

Main dependencies are:
* TCGAbiolinks
* GEOquery
* preprocessCore
* ROCR
* DESeq2
* biomaRt
* pacman

Other (minor) dependecies will be automatically installed if they are missing by `pacman` package within an execution of source scripts.

## 🚀 Quick start

Source scripts that can be used for a reproduction is located at `./scripts`.
Annotation of clinical biosamples with *time-to-progression (TTP)* and *over-survival (OS)* as well supplementary data files used across analyses can be found at `./annotations`

## 🆘 Help

Please feel free to contact Mikhail Raevskiy (raevskii.mm@phystech.edu) if you have any questions about the software.

## 📃 License

This project is [Apache 2.0](https://github.com/raevskymichail/CRNDE_glioblastoma/blob/main/LICENSE) licensed.