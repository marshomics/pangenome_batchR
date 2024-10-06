# pangenome_batchR
An increasing number of programs are available for microbial pangenome analysis, including Roary, Panaroo, Pandora, Pagoo, Panakeia, and GET_HOMOLOGUES. However, few offer standardised output files that can be used for more complex downstream analyses, especially for multiple pangenomes at once. In response, this small program was developed to take multiple standard gene presence-absence files that are produced by all pangenome programs and provide a number of output plots and data tables to inform a first-pass analysis of each pangenome. By leveragint the capabilities of Pagoo and other R packages, the output includes core, shell, and cloud pangenome matrices, PCA plots, core gene alignments, random forest-based feature selection, Heaps' law (power law regression) calculations, gene clusters under positive, negative, and neutral selection, and population structure analyses (with associated phylogenetic trees). See below for a detailed overview of these analyses and their output files. The priority of this program is to quickly produce a large amount of data associated with each input pangenome as a starting point for a more focused and refined analysis by the user.

## Installation

```git clone https://github.com/marshomics/pangenome_batchR.git```

## Input preparation




