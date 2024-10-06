# pangenome_batchR
An increasing number of programs are available for microbial pangenome analysis, including Roary, Panaroo, Pandora, Pagoo, Panakeia, and GET_HOMOLOGUES. However, few offer standardised output files that can be used for more complex downstream analyses, especially for multiple pangenomes at once. In response, this small program was developed to take multiple standard gene presence-absence files that are produced by all pangenome programs and provide a number of output plots and data tables to inform a first-pass analysis of each pangenome. By leveragint the capabilities of Pagoo and other R packages, the output includes core, shell, and cloud pangenome matrices, PCA plots, core gene alignments, random forest-based feature selection, Heaps' law (power law regression) calculations, gene clusters under positive, negative, and neutral selection, and population structure analyses (with associated phylogenetic trees). See below for a detailed overview of these analyses and their output files. The priority of this program is to quickly produce a large amount of data associated with each input pangenome as a starting point for a more focused and refined analysis by the user.

### How to cite
Coming soon

## Installation

```
git clone https://github.com/marshomics/pangenome_batchR.git
```

## Dependencies

`pangenome_batchR` has been tested with the following package versions. When running the program it will attempt to install any missing packages.

```
RColorBrewer: 1.1.3 
pagoo: 0.3.18 
ggplot2: 3.5.1 
Biostrings: 2.72.1 
data.table: 1.16.0 
dplyr: 1.1.4 
NbClust: 3.0.1 
randomForest: 4.7.1.1 
magrittr: 2.0.3 
DECIPHER: 3.0.0 
phangorn: 2.11.1 
ggtree: 3.12.0 
rhierbaps: 1.1.4 
optparse: 1.7.5 
yaml: 2.3.10 
tidyr: 1.3.1 
ape: 5.8 
pegas: 1.3 
patchwork: 1.2.0 
```

## Input preparation

`pangenome_batchR` takes a `config.yaml` file as the input which specifies the data directory, the categorical variables to be used for plotting, and the percentage pangenome definition. The data directory can either be a single directory containing the *.gff files used for pangenome construction, the presence-absence file, and metadata file, or a parent directory containing multiple sub-directories that each contain these files (i.e. when processing multiple pangenomes at once).

### Input file formats

The *.gff files should be those that were used by the pangenome construction program.

The presence-absence file should be a comma-separated pangenome matrix similar to what is produced by Roary.

The tab-separated metadata table should list the genome names in a column called `filename`, and these should match the assemblies that are listed in the gene presence-absence file and match the *.gff filenames. The remaining columns should include all categorical metadata associated with each genome that might be relevant when plotting (see below).

### `config.yaml` specifics

The `config.yaml` enables global parameters to be set prior to running the program. The input directory should be specified, as well as the column names in the metdata table that list the categorical variables. Additionally, the percentage of genomes a gene should be associated with for it to be considered "core" for both pangenome and phylogenetic analyses can be specified (core_level_pangenome: 95, core_level_phylogeny: 100 are the defaults).

## Running the program

```
Rscript \
    pagoo_batch_analysis.R \
    --config config.yaml
```

## Output files and analyses










