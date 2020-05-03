# R Script to join and analyze output from CiliaQ

This template markdown aims to facilitate and automate the analysis of tabular output data generated by CiliaQ. It searches, fetches and joins all CiliaQ output files ('*_CQs.txt*') allowing downstream analysis in R or to export the collected data as an .xslx (Excel) file. Additionally, the script implements basic quality control steps and aims to create a brief overview of the retrieved data.

## Installation
Downlaod the [repository](https://github.com/sRassmann/ciliaQ-analysis/archive/master.zip) or directly clone the project [using Git](https://happygitwithr.com/rstudio-git-github.html). The script was written and tested in R version 4.0.0, however, earlier versions might be functional, as well.  
__Note:__ You should clone the whole project structure, as it also contains the blueprints of additional documents (exclude lists, annotation lists) and custom functions imported into the actual script.

## Usage
Before running the script only the paths, names, and experimental conditions need to be modified (marked with the tag '*TODO*'). Hence, you on't need (higher) skills in R programming. All the steps are explained in the markdown parts of the documents.
