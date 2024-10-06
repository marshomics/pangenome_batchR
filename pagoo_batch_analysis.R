# Ensure 'remotes' package is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Ensure 'devtools' package is installed (needed to install pagoo from GitHub)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Required CRAN packages and their versions
cran_packages <- list(
  RColorBrewer = "1.1-3",
  ggplot2 = "3.5.1",
  data.table = "1.16.0",
  dplyr = "1.1.4",
  NbClust = "3.0",
  randomForest = "4.7-1.1",
  magrittr = "2.0.3",
  optparse = "1.7.5",
  yaml = "2.3.10",
  tidyr = "1.3.1",
  ape = "5.8",
  pegas = "1.3",
  patchwork = "1.1.2"
)

# Function to check and install CRAN packages if not installed
install_cran_if_needed <- function(pkg, version) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Package", pkg, "is not installed. Installing version", version, "from CRAN...\n")
    remotes::install_version(pkg, version = version)
  } else {
    cat(pkg, "is already installed.\n")
  }
  # Load the package
  library(pkg, character.only = TRUE)
}

# Install and load CRAN packages
for (pkg in names(cran_packages)) {
  version <- cran_packages[[pkg]]
  install_cran_if_needed(pkg, version)
}

# Install pagoo via GitHub if not installed
if (!requireNamespace("pagoo", quietly = TRUE)) {
  cat("Installing pagoo version 0.3.18 from GitHub...\n")
  devtools::install_github('iferres/pagoo', ref = 'v0.3.18')
} else {
  cat("pagoo is already installed.\n")
}
# Load pagoo
library(pagoo)

# Note about Bioconductor packages
cat("\nNote: Installing specific versions of Bioconductor packages can be challenging due to dependencies and compatibility issues.\n")
cat("For Bioconductor packages, the latest compatible versions will be installed if they are not already installed.\n")
cat("If specific versions are required, consider using package management tools like 'renv' or 'packrat' for reproducibility.\n\n")

# Ensure 'BiocManager' is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if not installed
bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", "ggtree", "rhierbaps")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from Bioconductor...\n")
    BiocManager::install(pkg)
  } else {
    cat(pkg, "is already installed.\n")
  }
  # Load the package
  library(pkg, character.only = TRUE)
}

# Load other required libraries
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dplyr)
library(NbClust)
library(randomForest)
library(magrittr)
library(optparse)
library(yaml)
library(tidyr)
library(ape)
library(pegas)
library(patchwork)

# Parse command-line arguments
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL,
              help = "Path to the config file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if config file is provided
if (is.null(opt$config)) {
  stop("Please provide a config file using -c or --config option")
}

# Read the config file
config <- yaml::read_yaml(opt$config)

# Extract input directory and parameters from the config
input_directory <- config$input_directory
categorical_variables <- config$categorical_variables
core_level_pangenome <- config$core_level_pangenome
core_level_phylogeny <- config$core_level_phylogeny

# Function to find directories containing the required files
find_input_directories <- function(base_dir) {
  # Required filenames
  required_files <- c("metadata.txt", "gene_presence_absence.csv")

  # Find all subdirectories
  dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  dirs <- c(base_dir, dirs)  # Include base directory

  # Filter directories that contain all required files and at least one .gff file
  valid_dirs <- dirs[sapply(dirs, function(d) {
    files_in_dir <- list.files(d)
    has_required_files <- all(required_files %in% files_in_dir)
    has_gff_files <- length(list.files(d, pattern = "\\.gff$", full.names = TRUE)) > 0
    has_required_files && has_gff_files
  })]

  return(valid_dirs)
}

# Find all valid input directories
input_dirs <- find_input_directories(input_directory)

# Check if any valid directories are found
if (length(input_dirs) == 0) {
  stop("No directories with the required files were found.")
}

# Process each input directory separately
for (dir_idx in seq_along(input_dirs)) {
  # Current input directory
  input_dir <- input_dirs[dir_idx]

  # Create output directory
  output_dir <- file.path(input_dir, "pagoo_output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Set up logging
  log_file <- file.path(output_dir, "analysis_log.txt")
  sink(log_file, append = TRUE, split = TRUE)

  # Start logging
  cat("Starting analysis for directory:", input_dir, "\n")

  # List of .gff files
  gffs <- list.files(path = input_dir, pattern = "\\.gff$", full.names = TRUE)

  # Path to gene presence/absence CSV file
  csv <- file.path(input_dir, "gene_presence_absence.csv")

  # Read metadata
  meta_path <- file.path(input_dir, "metadata.txt")
  meta <- read.table(file = meta_path, sep = "\t", header = TRUE)

  # Rename 'filename' column to 'org' if necessary
  colnames(meta)[colnames(meta) == "filename"] <- "org"

  # Load pangenome data using pagoo
  p <- roary_2_pagoo(gene_presence_absence_csv = csv, gffs = gffs)

  # Add metadata to pangenome object
  p$add_metadata(map = "org", meta)

  # Set core level percentage for pangenome plots
  p$core_level <- core_level_pangenome

  # Get summary statistics
  p_summary <- as.data.frame(p$summary_stats)

  # Write summary statistics to file
  write.table(p_summary, file = file.path(output_dir, "summary_stats.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write gene metadata
  p_gene_metadata <- as.data.frame(p$genes)
  write.table(p_gene_metadata, file = file.path(output_dir, "gene_metadata.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write cluster information
  p_clusters <- as.data.frame(p$clusters)
  write.table(p_clusters, file = file.path(output_dir, "clusters.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write shell clusters
  p_shell_clusters <- as.data.frame(p$shell_clusters)
  write.table(p_shell_clusters, file = file.path(output_dir, "shell_clusters.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write cloud clusters
  p_cloud_clusters <- as.data.frame(p$cloud_clusters)
  write.table(p_cloud_clusters, file = file.path(output_dir, "cloud_clusters.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write core clusters
  p_core_clusters <- as.data.frame(p$core_clusters)
  write.table(p_core_clusters, file = file.path(output_dir, "core_clusters.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write shell genes
  p_shell_genes <- as.data.frame(p$shell_genes)
  write.table(p_shell_genes, file = file.path(output_dir, "shell_genes.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write cloud genes
  p_cloud_genes <- as.data.frame(p$cloud_genes)
  write.table(p_cloud_genes, file = file.path(output_dir, "cloud_genes.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Extract and write core genes
  p_core_genes <- as.data.frame(p$core_genes)
  write.table(p_core_genes, file = file.path(output_dir, "core_genes.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Combine and write shell sequences to FASTA file
  combined_dna_shell <- unlist(p$shell_sequences, use.names = TRUE)
  writeXStringSet(combined_dna_shell, filepath = file.path(output_dir, "shell_genes.fasta"))

  # Combine and write core sequences to FASTA file
  combined_dna_core <- unlist(p$core_sequences, use.names = TRUE)
  writeXStringSet(combined_dna_core, filepath = file.path(output_dir, "core_genes.fasta"))

  # Combine and write cloud sequences to FASTA file
  combined_dna_cloud <- unlist(p$cloud_sequences, use.names = TRUE)
  writeXStringSet(combined_dna_cloud, filepath = file.path(output_dir, "cloud_genes.fasta"))

  # Perform PCA and generate plots for each specified categorical variable
  for (cat_var in categorical_variables) {
    if (cat_var %in% colnames(meta)) {
      pca <- p$pan_pca()
      p$gg_pca(colour = cat_var, size = 4) +
        theme_bw(base_size = 15) +
        scale_color_brewer(palette = "Set2")

      # Save PCA plot
      ggsave(filename = file.path(output_dir, paste0("pca_", cat_var, ".png")))
    } else {
      warning(paste("Categorical variable", cat_var, "not found in metadata."))
    }
  }

  # Generate pan-genome curves plot
  p$gg_curves(size = 2) +
    geom_point() +
    facet_wrap(~Category, scales = 'free_y') +
    theme_bw(base_size = 15) +
    scale_color_brewer(palette = "Accent")

  # Save curves plot
  ggsave(filename = file.path(output_dir, "curves.png"), width = 16, height = 8)

  # Fit pan-genome to power law and write parameters to file
  fit_result <- p$pg_power_law_fit()
  params <- fit_result$params
  alpha <- attr(fit_result, "alpha")
  result_df <- data.frame(Parameter = c("K", "delta", "alpha"),
                          Value = c(params["K"], params["delta"], alpha))
  write.table(result_df, file = file.path(output_dir, "power.txt"),
              sep = "\t", col.names = NA)

  # Generate pie chart of pan-genome composition
  p$gg_pie() + theme_bw(base_size = 15) +
    scale_fill_brewer(palette = "Blues") +
    scale_x_discrete(breaks = c(0, 25, 50, 75)) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank())

  # Save pie chart
  ggsave(filename = file.path(output_dir, "pie.png"))

  # Additional analysis steps
  # Prepare pan-genome matrix and filter clusters
  pan_matrix <- p$pan_matrix

  # Extract clusters based on their categories
  core_clusters <- unique(as.character(p$core_clusters$cluster))
  shell_clusters <- unique(as.character(p$shell_clusters$cluster))
  cloud_clusters <- unique(as.character(p$cloud_clusters$cluster))

  # Filter pan-genome matrix for core, shell, and cloud clusters
  filtered_pan_matrix_core <- pan_matrix[, colnames(pan_matrix) %in% core_clusters]
  filtered_pan_matrix_shell <- pan_matrix[, colnames(pan_matrix) %in% shell_clusters]
  filtered_pan_matrix_cloud <- pan_matrix[, colnames(pan_matrix) %in% cloud_clusters]

  # Remove rows with all zeros
  non_zero_row_indices_core <- which(rowSums(filtered_pan_matrix_core != 0) > 0)
  non_zero_row_indices_shell <- which(rowSums(filtered_pan_matrix_shell != 0) > 0)
  non_zero_row_indices_cloud <- which(rowSums(filtered_pan_matrix_cloud != 0) > 0)

  final_matrix_core <- filtered_pan_matrix_core[non_zero_row_indices_core, ]
  final_matrix_shell <- filtered_pan_matrix_shell[non_zero_row_indices_shell, ]
  final_matrix_cloud <- filtered_pan_matrix_cloud[non_zero_row_indices_cloud, ]

  # Write filtered pan-genome matrices to files
  write.table(final_matrix_core, file = file.path(output_dir, "core_pan_matrix.txt"),
              sep = "\t", col.names = NA, quote = FALSE)
  write.table(final_matrix_shell, file = file.path(output_dir, "shell_pan_matrix.txt"),
              sep = "\t", col.names = NA, quote = FALSE)
  write.table(final_matrix_cloud, file = file.path(output_dir, "cloud_pan_matrix.txt"),
              sep = "\t", col.names = NA, quote = FALSE)

  # Read gene presence/absence CSV
  csv_pres_abs <- read.table(file = csv, sep = ",", header = TRUE, check.names = FALSE)

  # Filter CSV for core, shell, and cloud clusters
  final_matrix_core_colnames <- colnames(final_matrix_core)
  final_matrix_shell_colnames <- colnames(final_matrix_shell)
  final_matrix_cloud_colnames <- colnames(final_matrix_cloud)

  filtered_csv_core <- csv_pres_abs[csv_pres_abs$Gene %in% final_matrix_core_colnames, ]
  filtered_csv_shell <- csv_pres_abs[csv_pres_abs$Gene %in% final_matrix_shell_colnames, ]
  filtered_csv_cloud <- csv_pres_abs[csv_pres_abs$Gene %in% final_matrix_cloud_colnames, ]

  # Write filtered CSVs to files
  write.table(filtered_csv_core, file = file.path(output_dir, "presence_absence_core.csv"),
              sep = ",", row.names = FALSE)
  write.table(filtered_csv_shell, file = file.path(output_dir, "presence_absence_shell.csv"),
              sep = ",", row.names = FALSE)
  write.table(filtered_csv_cloud, file = file.path(output_dir, "presence_absence_cloud.csv"),
              sep = ",", row.names = FALSE)

  # Create pagoo objects for core and shell genes
  gffs_core <- gffs
  csv_core <- file.path(output_dir, "presence_absence_core.csv")
  meta_core <- meta

  p_core <- roary_2_pagoo(gene_presence_absence_csv = csv_core, gffs = gffs_core)
  p_core$add_metadata(map = "org", meta_core)

  gffs_shell <- gffs
  csv_shell <- file.path(output_dir, "presence_absence_shell.csv")
  meta_shell <- meta

  p_shell <- roary_2_pagoo(gene_presence_absence_csv = csv_shell, gffs = gffs_shell)
  p_shell$add_metadata(map = "org", meta_shell)

  # Perform PCA on core and shell genes and generate plots dynamically
  for (cat_var in categorical_variables) {
    if (cat_var %in% colnames(meta_core)) {
      # PCA for core genes
      pca_core <- p_core$pan_pca()
      p_core$gg_pca(colour = cat_var, size = 4) +
        theme_bw(base_size = 15) +
        scale_color_brewer(palette = "Set2")

      # Save PCA plot for core genes
      ggsave(filename = file.path(output_dir, paste0("pca_core_", cat_var, ".png")))
    }
    if (cat_var %in% colnames(meta_shell)) {
      # PCA for shell genes
      pca_shell <- p_shell$pan_pca()
      p_shell$gg_pca(colour = cat_var, size = 4) +
        theme_bw(base_size = 15) +
        scale_color_brewer(palette = "Set2")

      # Save PCA plot for shell genes
      ggsave(filename = file.path(output_dir, paste0("pca_shell_", cat_var, ".png")))
    }
  }

  # Convert PCA results to data frames
  pca_df_core <- data.frame(pca_core$x)
  pca_df_shell <- data.frame(pca_shell$x)

  # Determine optimal number of clusters using NbClust
  set.seed(123)
  nb_core <- NbClust(pca_df_core[, c("PC1", "PC2")], distance = "euclidean",
                     min.nc = 2, max.nc = 10, method = "complete")
  nb_shell <- NbClust(pca_df_shell[, c("PC1", "PC2")], distance = "euclidean",
                      min.nc = 2, max.nc = 10, method = "complete")

  num_clusters_core <- length(unique(nb_core$Best.partition))
  num_clusters_shell <- length(unique(nb_shell$Best.partition))

  # K-means clustering on PCA coordinates
  set.seed(123)
  kmeans_result_core <- kmeans(pca_df_core[, 1:2], centers = num_clusters_core)
  kmeans_result_shell <- kmeans(pca_df_shell[, 1:2], centers = num_clusters_shell)

  # Add cluster assignments to PCA data frames
  pca_df_core$cluster <- as.factor(kmeans_result_core$cluster)
  pca_df_shell$cluster <- as.factor(kmeans_result_shell$cluster)

  # Save cluster assignments
  sample_clusters_core <- data.frame(sample = rownames(pca_df_core), cluster = pca_df_core$cluster)
  sample_clusters_shell <- data.frame(sample = rownames(pca_df_shell), cluster = pca_df_shell$cluster)

  write.csv(sample_clusters_core, file = file.path(output_dir, "sample_clusters_core.txt"),
            row.names = FALSE)
  write.csv(sample_clusters_shell, file = file.path(output_dir, "sample_clusters_shell.txt"),
            row.names = FALSE)

  # Plot PCA with clusters
  ggplot(pca_df_core, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 4) +
    theme_bw(base_size = 15) +
    scale_color_brewer(palette = "Set2") +
    labs(title = "PCA Plot with Clusters (Core Genes)")

  ggsave(filename = file.path(output_dir, "pca_core_clusters.png"))

  ggplot(pca_df_shell, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 4) +
    theme_bw(base_size = 15) +
    scale_color_brewer(palette = "Set2") +
    labs(title = "PCA Plot with Clusters (Shell Genes)")

  ggsave(filename = file.path(output_dir, "pca_shell_clusters.png"))

  # Prepare pan-genome matrices as data frames and add sample names
  pan_matrix_core_df <- as.data.frame(p_core$pan_matrix)
  pan_matrix_core_df$sample <- rownames(pan_matrix_core_df)

  pan_matrix_shell_df <- as.data.frame(p_shell$pan_matrix)
  pan_matrix_shell_df$sample <- rownames(pan_matrix_shell_df)

  # Ensure 'meta' contains the necessary columns and 'sample' column
  colnames(meta)[colnames(meta) == "org"] <- "sample"

  # Function to perform pairwise feature selection
  perform_pairwise_feature_selection <- function(pan_matrix_df, response_df, response_var, output_prefix) {
  # Get unique categories or clusters
  categories <- unique(response_df[[response_var]])
  
  # Generate all unique pairs of categories
  category_pairs <- combn(categories, 2, simplify = FALSE)
  
  # Loop over each pair of categories
  for (pair in category_pairs) {
    category_a <- pair[1]
    category_b <- pair[2]
    pair_label <- paste0(response_var, "_", category_a, "_vs_", category_b)
    cat("Processing", pair_label, "\n")
    
    # Subset data for the two categories
    samples_in_pair <- response_df$sample[response_df[[response_var]] %in% c(category_a, category_b)]
    pan_matrix_pair <- pan_matrix_df[pan_matrix_df$sample %in% samples_in_pair, ]
    response_df_pair <- response_df[response_df$sample %in% samples_in_pair, ]
    
    # Prepare response variable
    response <- as.factor(response_df_pair[[response_var]])
    
    # Check if there are enough samples for both categories
    if (length(unique(response)) < 2) {
      cat("Not enough samples for both categories in", pair_label, ". Skipping...\n")
      next
    }
    
    # Prepare predictors (exclude 'sample' column)
    predictors <- pan_matrix_pair %>%
      select(-sample)
    
    # Check if predictors data frame is empty
    if (ncol(predictors) == 0) {
      cat("No predictors available for", pair_label, ". Skipping...\n")
      next
    }
    
    # Ensure predictors have variance
    variances <- apply(predictors, 2, var, na.rm = TRUE)
    
    # Replace NA variances with zero
    variances[is.na(variances)] <- 0
    
    variable_columns <- variances != 0
    
    # Check if lengths match
    if (length(variable_columns) != ncol(predictors)) {
      cat("Mismatch in variable_columns length for", pair_label, ". Skipping...\n")
      next
    }
    
    # Check if there are any variable predictors
    if (!any(variable_columns)) {
      cat("No variable predictors for", pair_label, ". Skipping...\n")
      next
    }
    
    # Subset predictors to variable columns
    predictors <- predictors[, variable_columns, drop = FALSE]
    
    # Train Random Forest model
    set.seed(123)
    rf_model <- randomForest(x = predictors, y = response, importance = TRUE)
    
    # Extract feature importance
    feature_importance <- importance(rf_model, type = 1)
    feature_importance_df <- data.frame(Feature = rownames(feature_importance),
                                        Importance = feature_importance[, 1])
    
    # Sort by importance
    feature_importance_df <- feature_importance_df %>%
      arrange(desc(Importance))
    
    # Save feature importance results
    output_file <- file.path(output_dir, paste0("feature_importance_", output_prefix, "_", pair_label, ".txt"))
    write.table(feature_importance_df, file = output_file, sep = "\t", row.names = FALSE)
    
    # Plot top N important features
    top_n <- min(50, nrow(feature_importance_df))
    if (top_n > 0) {
      ggplot(feature_importance_df[1:top_n, ], aes(x = reorder(Feature, Importance), y = Importance)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        xlab("Genes") +
        ylab("Importance") +
        ggtitle(paste("Feature Importance for", pair_label, "(", output_prefix, "Genes)")) +
        theme_bw(base_size = 15)
      
      # Save plot
      plot_file <- file.path(output_dir, paste0("feature_importance_", output_prefix, "_", pair_label, "_top", top_n, ".png"))
      ggsave(filename = plot_file)
    } else {
      cat("No features to plot for", pair_label, ". Skipping plot.\n")
    }
  }
}

  # Perform pairwise feature selection for each categorical variable
  for (cat_var in categorical_variables) {
    if (cat_var %in% colnames(meta)) {
      cat("Performing feature selection for categorical variable:", cat_var, "\n")

      # Prepare response data frame with 'sample' and the categorical variable
      response_df <- meta[, c("sample", cat_var)]

      # Remove samples with missing values in the categorical variable
      response_df <- response_df[!is.na(response_df[[cat_var]]), ]

      # Merge response data with pan-genome matrices
      pan_matrix_core_df_response <- merge(pan_matrix_core_df, response_df, by = "sample")
      pan_matrix_shell_df_response <- merge(pan_matrix_shell_df, response_df, by = "sample")

      # Perform pairwise feature selection for core genes
      perform_pairwise_feature_selection(
        pan_matrix_df = pan_matrix_core_df_response,
        response_df = response_df,
        response_var = cat_var,
        output_prefix = paste0("core_", cat_var)
      )

      # Perform pairwise feature selection for shell genes
      perform_pairwise_feature_selection(
        pan_matrix_df = pan_matrix_shell_df_response,
        response_df = response_df,
        response_var = cat_var,
        output_prefix = paste0("shell_", cat_var)
      )
    } else {
      warning(paste("Categorical variable", cat_var, "not found in metadata."))
    }
  }

  # Align individual core genes at specified core level for phylogeny
  p$core_level <- core_level_phylogeny
  
  # Extract sequences for alignment
  seqs <- p$core_seqs_4_phylo()
  
  # Align sequences using DECIPHER
  ali <- seqs %>%
    lapply(DECIPHER::AlignTranslation)
  
  # Calculate Tajima's D for each alignment
  tajD <- ali %>%
    lapply(ape::as.DNAbin) %>%
    lapply(pegas::tajima.test) %>%
    sapply("[[", "D")
  
  # Identify genes under neutral selection
  neutral <- which(tajD <= 2 & tajD >= -2)
  
  # Create data frame of Tajima's D results
  tajD_df <- data.frame(
    cluster = names(tajD),
    tajD = as.numeric(tajD)
  )
  
  tajD_df$class <- ifelse(tajD_df$tajD >= -2 & tajD_df$tajD <= 2, "neutral",
                          ifelse(tajD_df$tajD > 2, "positive", "negative"))
  
  # Write Tajima's D results to file
  write.table(tajD_df, file = file.path(output_dir, "selection_commongenes.txt"),
              sep = "\t", row.names = FALSE)
  
  # Extract neutral alignments
  neutral_ali <- ali[neutral]
  
  # Define directory for individual alignments
  dir_path <- file.path(output_dir, "neutral_gene_alignment_commongenes")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Write individual alignments to FASTA files
  for (i in seq_along(neutral_ali)) {
    alignment <- neutral_ali[[i]]
    file_name <- file.path(dir_path, paste0("alignment", i, ".fasta"))
    writeXStringSet(alignment, filepath = file_name, format = "fasta")
  }
  
  # Function to remove columns with gaps in 50% or more sequences
  remove_gaps_manual <- function(dna_set) {
    dna_matrix <- as.matrix(dna_set)
    gap_proportion <- colMeans(dna_matrix == "-")
    dna_matrix_nogaps <- dna_matrix[, gap_proportion < 0.50]
    return(DNAStringSet(apply(dna_matrix_nogaps, 1, paste0, collapse = "")))
  }
  
  # Apply gap removal function
  filtered_dnaStringSetList <- DNAStringSetList(lapply(neutral_ali, remove_gaps_manual))
  
  # Define directory for gap-filtered alignments
  dir_path <- file.path(output_dir, "neutral_gene_alignment_commongenes_nogaps")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Write gap-filtered alignments to FASTA files
  for (i in seq_along(filtered_dnaStringSetList)) {
    alignment <- filtered_dnaStringSetList[[i]]
    file_name <- file.path(dir_path, paste0("alignment", i, ".fasta"))
    writeXStringSet(alignment, filepath = file_name, format = "fasta")
  }
  
  # Concatenate neutral core gene clusters
  concatenated_sequences <- do.call(Biostrings::xscat, neutral_ali)
  names(concatenated_sequences) <- p$organisms$org
  output_alignment <- DNAStringSet(concatenated_sequences)
  writeXStringSet(output_alignment, filepath = file.path(output_dir, "concatenated_alignment.fasta"))
  
  # Concatenate gap-filtered sequences
  concatenated_sequences_nogaps <- do.call(Biostrings::xscat, filtered_dnaStringSetList)
  names(concatenated_sequences_nogaps) <- p$organisms$org
  output_alignment_nogaps <- DNAStringSet(concatenated_sequences_nogaps)
  writeXStringSet(output_alignment_nogaps, filepath = file.path(output_dir, "concatenated_alignment_nogaps.fasta"))
  
  # Prepare data for rhierbaps
  concat_neu <- do.call(Biostrings::xscat, neutral_ali) %>%
    setNames(p$organisms$org) %>%
    as("matrix") %>%
    tolower()
  
  concat_neu_nogaps <- do.call(Biostrings::xscat, filtered_dnaStringSetList) %>%
    setNames(p$organisms$org) %>%
    as("matrix") %>%
    tolower()
  
  # Run rhierbaps for population structure
  set.seed(123)
  rhb <- hierBAPS(snp.matrix = concat_neu, n.pops = 10,
                  max.depth = 1, n.extra.rounds = 5, n.cores = 10)
  
  set.seed(123)
  rhb_nogaps <- hierBAPS(snp.matrix = concat_neu_nogaps, n.pops = 10,
                         max.depth = 1, n.extra.rounds = 5, n.cores = 10)
  
  rhb_df <- as.data.frame(rhb$partition.df)
  rhb_df_nogaps <- as.data.frame(rhb_nogaps$partition.df)
  
  # Write rhierbaps results to files
  write.table(rhb_df, file = file.path(output_dir, "tree_clusters_commongenes.txt"),
              sep = "\t", row.names = FALSE)
  write.table(rhb_df_nogaps, file = file.path(output_dir, "tree_clusters_commongenes_nogaps.txt"),
              sep = "\t", row.names = FALSE)
  
  # Add lineage information to organisms metadata
  res <- rhb$partition.df
  res_nogaps <- rhb_nogaps$partition.df
  
  lin <- data.frame(org = as.character(res[, 1]), lineage = as.factor(res[, 2]))
  lin_nogaps <- data.frame(org = as.character(res_nogaps[, 1]), lineage = as.factor(res_nogaps[, 2]))
  
  # Create copies of pangenome objects for 'nogaps' versions
  p_nogaps <- p
  p_nogaps$add_metadata(map = "org", data = lin_nogaps)
  
  p$add_metadata(map = "org", data = lin)
  
  # Compute phylogenetic trees
  tre <- concat_neu %>%
    phangorn::phyDat(type = "DNA") %>%
    phangorn::dist.ml() %>%
    phangorn::NJ()
  
  tre_nogaps <- concat_neu_nogaps %>%
    phangorn::phyDat(type = "DNA") %>%
    phangorn::dist.ml() %>%
    phangorn::NJ()
  
  # Adjust branch lengths by adding 0.1 to avoid zero-length branches
  tre$edge.length <- tre$edge.length + 0.1
  tre_nogaps$edge.length <- tre_nogaps$edge.length + 0.1
  
  # Function to generate and save trees
  generate_and_save_trees <- function(tre, pangenome_object, layout_type, version_label) {
    # Always generate the 'lineage' tree
    gg_lineage <- ggtree(tre, ladderize = TRUE, layout = layout_type) %<+%
      as.data.frame(pangenome_object$organisms) +
      geom_tippoint(aes(color = as.factor(lineage))) +
      labs(subtitle = paste("Lineage -", version_label, "-", layout_type)) +
      scale_color_discrete("Lineage")
    
    # Save 'lineage' tree plot
    ggsave(filename = file.path(output_dir, paste0("tree_lineage_", version_label, "_", layout_type, ".png")))
    
    # Generate trees for each categorical variable
    for (cat_var in categorical_variables) {
      if (cat_var %in% colnames(pangenome_object$organisms)) {
        gg_tree <- ggtree(tre, ladderize = TRUE, layout = layout_type) %<+%
          as.data.frame(pangenome_object$organisms) +
          geom_tippoint(aes_string(colour = cat_var)) +
          labs(subtitle = paste(cat_var, "-", version_label, "-", layout_type)) +
          scale_colour_discrete(cat_var)
        
        # Save tree plot
        ggsave(filename = file.path(output_dir, paste0("tree_", cat_var, "_", version_label, "_", layout_type, ".png")))
        
        # Tree with tip labels
        gg_tree_tips <- gg_tree +
          geom_tiplab()
        
        # Save tree plot with tip labels
        ggsave(filename = file.path(output_dir, paste0("tree_", cat_var, "_", version_label, "_", layout_type, "_tips.png")),
               plot = gg_tree_tips, height = 30, width = 10)
      } else {
        warning(paste("Variable", cat_var, "not found in organisms metadata."))
      }
    }
  }
  
  # Generate and save trees for 'gg' version (with gaps)
  for (layout in c("slanted", "circular")) {
    generate_and_save_trees(tre, p, layout_type = layout, version_label = "commongenes")
  }
  
  # Generate and save trees for 'gg_nogaps' version (without gaps)
  for (layout in c("slanted", "circular")) {
    generate_and_save_trees(tre_nogaps, p_nogaps, layout_type = layout, version_label = "commongenes_nogaps")
  }
  
  # Save phylogenetic trees in Newick format
  write.tree(tre, file = file.path(output_dir, "population_structure_tree_commongenes.nwk"))
  write.tree(tre_nogaps, file = file.path(output_dir, "population_structure_tree_commongenes_nogaps.nwk"))

  # Close the log file
  sink()

  cat("Analysis completed for directory:", input_dir, "\n")
}