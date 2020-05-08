
# Codon Usage Analysis

The purpose of this project is to analyzes codon usage of selected genes
based on microarray samples. for microarray analysis the polysomal
fraction was separated to small and large polysomes, which was then
analyzed by microarray analysis. To understand whether there is a codon
usage preferrence associated with the polysomal fraction of the gene
level. The following genes are analyzed:

1.  The top 10 percent of genes with the highest large polysomal
    fractions (LP)

2.  The top 10 percent of genes with the highest small polysomal
    fractions (SP).

3.  The bottom 10 percent of genes with the highest large polysomal
    fractions (LP)

4.  The bottom 10 percent of genes with the highest small polysomal
    fractions (SP).

The results are used to guide downstream experiments.

Data and results are not included in this repo. Some sample plots are
shown in the `results` file.

## Usage

To replicate the analysis, clone this GitHub repository, install the
[dependencies](#dependencies) listed below, and run the following
commands at the command line/terminal from the root directory of this
project:

    make all

To reset the repo to a clean state, with no intermediate or results
files, run the following command at the command line/terminal from the
root directory of this project:

    make clean

## Dependencies

  - Python 3.7.4 and Python packages:
      - numpy==1.17.2
      - scipy==1.4.1
      - matplotlib==3.1.1
