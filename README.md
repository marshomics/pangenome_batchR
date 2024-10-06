# pangenome_batchR
An increasing number of programs are available for microbial pangenome analysis, including Roary, Panaroo, Pandora, Pagoo, Panakeia, and GET_HOMOLOGUES. However, few offer standardised output files that can be used for more complex downstream analyses, especially for multiple pangenomes at once. In response, this small program was developed to take multiple standard gene presence-absence files that are produced by all pangenome programs and provide a number of output plots and data tables to inform a first-pass analysis of each pangenome. By leveraging the capabilities of Pagoo and other R packages, the output includes core, shell, and cloud pangenome matrices, PCA plots, core gene alignments, random forest-based feature selection, Heaps' law (power law regression) calculations, gene clusters under positive, negative, and neutral selection, and population structure analyses (with associated phylogenetic trees). See below for a detailed overview of these analyses and their output files. The priority of this program is to quickly produce a large amount of data associated with each input pangenome as a starting point for a more focused and refined analysis by the user.

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

The `config.yaml` enables global parameters to be set prior to running the program. The input directory should be specified, as well as the column names in the metdata table that list the categorical variables. Additionally, the percentage of genomes a gene should be associated with for it to be considered "core" for both pangenome and phylogenetic analyses can be specified (`core_level_pangenome`: 95, `core_level_phylogeny`: 100 are the defaults).

## Running the program

```
Rscript \
    pagoo_batch_analysis.R \
    --config config.yaml
```

## Output files and analyses

The following files and plots will be produced:

1. **Summary statistics:** Pangenome summary statistics listing the number of core, shell, and cloud genes.
2. **Extract gene metadata:** Table listing the genes and contig associated with each pangenome cluster, their genome coordinates, strand, and gene products (if included in the *.gff files). Additionally, separated tables for core, cloud, and shell genes are provided.
3. **Extract cluster metadata:** Tables listing the pangenome clusters and their functional annotation for core, shell, and cloud genes.
4. **Core, cloud, and shell gene sequences:** Fasta files of the genes associated with each pangenome section.
5. **Core and accessory genome PCA plots:** A separate plot is produced for each categorical variable (column) in the metadata, with the points coloured according to metadata category.
6. **Pangenome curves:** Curves plotting the number of gene clusters that are identified as more genomes are added to the dataset. As implemented in Pagoo, the Power-law distribution is used to fit the pangenome size and the Exponential decay function to fit the core genome size.
7. **Heaps' law:** Numeric estimate of the Heaps' law value indicating whether the pangenome is open or closed.
8. **Pangenome pie chart:** Basic pie chart to show core, shell, and cloud genome proportions.
9. **Filtered pangenome matrices:** The pangenome matrix subset to include only core, shell, or cloud genes.
10. **Filtered pangenome presence-absence tables:** Presence-absence tables subset to include only core, shell, or cloud genes.
11. **PCA partitioning and k-means clustering:** Optimal clusters are identified in the core and shell PCA plots and the genomes associated with each cluster are written to file. Additional PCA plots are produced with the points coloured according to cluster membership.
12. **Random forest feature selection:** Identification of the features (pangenome gene clusters) that are driving the separation between PCA clusters, well as features differentiating the categorical groups provided in the metadata table. The user should determine whether these are relevant comparisons by manually inspecting the PCA clusters. Corresponding tables and bar charts are produced.
13. **Core gene codon-aware alignment:** Alignment fasta is produced.
14. **Core genes under positive, negative, or neutral selection:** Tajima's D is calculated for the neutral core genes. Neutral core gene alignments with and without gaps as well as a concatenated alignment with and without gaps is produced, as well as tables listing each gene and whether it is likely to be under positive, negative, or neutral selection.
15. **Pangenome population structure:** Neutral core genes are hierarchically clustered to reveal nested population structures for each pangenome. Phylogenetic tree plots are produced, colored according to the predicted evolutionary lineaage as well as each categorical variable listed in the metadata table.
16. **Neutral core gene phylogenetic trees:** Newick phylogenetic tree files are produced for both the neutral core genes and the neutral core genes with gaps removed from the alignment.

## Next steps

## References










