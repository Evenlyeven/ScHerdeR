suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(ggthemes)))
suppressMessages(suppressWarnings(library(writexl)))

## ===== define a function of the analysis ===== ##
ScHerdeR <- function(seurat_obj,
                     output_dir,
                     npc,
                     resolution,
                     reduction,
                     sample_ident,
                     perc_calculation = TRUE,
                     cols) {
  ## ----- handle input file ----- ##
  # User may input .rds or RData of the seurat_obj, we will read it in
  if (grepl("\\.rds$", seurat_obj, ignore.case = TRUE)) {
    seurat_obj <- readRDS(seurat_obj)
  } else if (grepl("\\.RData$", seurat_obj, ignore.case = TRUE)) {
    temp_env <- new.env()
    load(seurat_obj, envir = temp_env)
    
    # Look for the first Seurat object in the loaded environment
    seurat_vars <- ls(temp_env)
    seurat_objs <- Filter(function(x)
      inherits(temp_env[[x]], "Seurat"), seurat_vars)
    
    if (length(seurat_objs) == 0) {
      stop("No Seurat object found in the RData file.")
    } else if (length(seurat_objs) > 1) {
      warning("Multiple Seurat objects found. Using the first one: ",
              seurat_objs[1])
    }
    
    seurat_obj <- temp_env[[seurat_objs[1]]]
    message("Using Seurat object from .RData: ", seurat_objs[1])
  } else {
    stop("The input file must be either .rds or .RData format.")
  }
  
  # Check if the input is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The input must be a Seurat object.")
  }
  
  # check if the required metadata columns exist
  if (!sample_ident %in% colnames(seurat_obj@meta.data)) {
    stop(
      paste(
        "The sample identity column",
        sample_ident,
        "does not exist in the metadata of the Seurat object."
      )
    )
  }
  
  ## ----- handle output directory ----- ##
  # Wrap output in a timestamped folder, so repeated runs do not overwrite results unless desired
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(output_dir, paste0("ScHerdeR_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----- load color vector if provided ----- #
  if (!is.null(cols)) {
    if (grepl("\\.rds$", cols, ignore.case = TRUE)) {
      cols <- readRDS(cols)
    } else if (grepl("\\.RData$", cols, ignore.case = TRUE)) {
      temp_env <- new.env()
      load(cols, envir = temp_env)
      
      # Look for the first vector or factor
      color_vars <- ls(temp_env)
      vector_objs <- Filter(function(x) is.atomic(temp_env[[x]]), color_vars)
      
      if (length(vector_objs) == 0) {
        stop("No color vector found in the RData file.")
      } else if (length(vector_objs) > 1) {
        warning("Multiple objects found in color file. Using the first one: ", vector_objs[1])
      }
      
      cols <- temp_env[[vector_objs[1]]]
      message("Using color vector from .RData: ", vector_objs[1])
    } else {
      stop("The color file must be either .rds or .RData format.")
    }
  }
  
  ## ----- run the analysis ----- ##
  # Dimensional Reduction
  set.seed(10086)
  
  i <- as.numeric(npc)
  j <- as.numeric(resolution)
  
  seu_m.c <- seurat_obj %>%
    FindNeighbors(dims = 1:i, reduction = reduction) %>%
    FindClusters(resolution = j) %>%
    RunUMAP(dims = 1:i, reduction = reduction)
  
  # save the dimensional reduced seurat object
  saveRDS(seu_m.c,
          file = file.path(output_dir, paste0("seurat_obj_PC", i, "_res", j, ".rds")))
  message(
    "Saved the Seurat object with PCA and clustering results to: ",
    file.path(output_dir, paste0("seurat_obj_PC", i, "_res", j, ".rds"))
  )
  
  # UMAPs
  # Helper function to generate and save DimPlot
  make_umap_plot <- function(seu,
                             filename,
                             group.by = NULL,
                             split.by = NULL,
                             label = TRUE,
                             repel = TRUE,
                             label.box = TRUE,
                             ncol = 2,
                             width = 1000,
                             height = 800,
                             cols = NULL) {
    plot_args <- list(
      object = seu,
      reduction = "umap",
      raster = FALSE,
      label = label,
      repel = repel,
      label.box = label.box
    )
    
    if (!is.null(group.by)) {
      plot_args$group.by <- group.by
    }
    if (!is.null(split.by)) {
      plot_args$split.by <- split.by
      plot_args$ncol <- ncol
    }
    if (!is.null(cols)) {
      plot_args$cols <- cols
    }
    
    p <- do.call(DimPlot, plot_args) + coord_fixed()
    
    if (!is.null(split.by)) {
      n_samps <- length(unique(seu@meta.data[[sample_ident]]))
      nrow <- ceiling(n_samps / ncol)
      height <- height * nrow * 0.6  # adjust height based on rows
      width <- width * 0.6 * 2
    } else {
      width <- width * 0.75
      height <- height * 0.75
    }
    
    png(file.path(output_dir, filename),
        width = width,
        height = height)
    print(p)
    dev.off()
    
    message("UMAP plot saved to: ", file.path(output_dir, filename))
  }
  
  make_umap_plot(seu_m.c, "umap.png", cols = cols, group.by = "seurat_clusters")
  
  make_umap_plot(seu_m.c, "umap_grouped_by_sample.png", group.by = sample_ident, label.box = FALSE, label = FALSE)
  
  make_umap_plot(
    seu_m.c,
    "umap_split_by_sample.png",
    split.by = sample_ident,
    group.by = "seurat_clusters",
    cols = cols
  )
  
  # generate a cell count by sample table
  ## We will check if any sample contains zero counts in any cluster
  #plot for each sample, how many cells in each cluster
  df_c <- seu_m.c@meta.data
  df_c$cell_ID <- rownames(df_c)
  
  df_c_keep <- df_c %>%
    select(seurat_clusters, !!as.name(sample_ident), cell_ID)
  rownames(df_c_keep) <- NULL
  
  df_perc <- df_c_keep %>%
    group_by(!!as.name(sample_ident), seurat_clusters) %>%
    summarize(cell_count = n()) %>%
    group_by(!!as.name(sample_ident)) %>%
    mutate(perc = cell_count / sum(cell_count) * 100) %>%
    tidyr::complete(seurat_clusters, fill = list(cell_count = 0, perc = 0))
  
  df_perc$seurat_clusters <- paste0("Cluster_", df_perc$seurat_clusters)
  
  saveRDS(df_perc, file = file.path(output_dir, "cell_count_by_sample.rds"))
  
  message(
    "Saved the cell count by sample table to: ",
    file.path(output_dir, "cell_count_by_sample.rds")
  )
  
  # Optionally output excel file of the cell count table and visualization
  # Actually we already calculated this above
  if (perc_calculation) {
    # save the table
    writexl::write_xlsx(df_perc, path = file.path(output_dir, "cell_count_by_sample.xlsx"))
    
    message(
      "Saved the cell count by sample table to: ",
      file.path(output_dir, "cell_count_by_sample.xlsx")
    )
    
    # plot the cell count by sample and save the plot
    pc <- ggplot(data = df_perc,
                 mapping = aes(
                   x = seurat_clusters,
                   y = perc,
                   fill = !!as.name(sample_ident)
                 )) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(
        aes(label = cell_count),
        hjust = -0.25,
        color = "black",
        position = position_dodge(width = 1),
        size = 4
      ) +
      xlab("Cell cluster") +
      ylab("Percentage of cell counts in each cluster (%)") +
      ggthemes::theme_few() +
      theme(
        axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "right",
        text = element_text(size = 20)
      ) +
      coord_flip()
    
    if (!is.null(cols)) {
      pcc <- pc + scale_fill_manual(values = cols)
      
      png(
        file.path(output_dir, "cell_count_by_sample.png"),
        width = 1000,
        height = 200 + 50 * length(unique(seu_m.c$seurat_clusters))
      )
      print(pcc)
      dev.off()
    } else {
      png(
        file.path(output_dir, "cell_count_by_sample.png"),
        width = 1000,
        height = 200 + 50 * length(unique(seu_m.c$seurat_clusters))
      )
      print(pc)
      dev.off()
    }
    
    message(
      "Saved the cell count by sample plot to: ",
      file.path(output_dir, "cell_count_by_sample.png")
    )
  }
}

## ===== define options for the script ===== ##
description_text <- '
This script takes a Seurat object (.rds/.RData), performs dimensionality reduction, clustering, and UMAP visualization using user-specified number of principal components and clustering resolution. It evaluates clustering consistency across samples, and outputs UMAPs grouped and split by sample identity. It also generates a summary table of cell counts and percentages per cluster per sample, and optionally visualizes the percentage breakdown as a bar plot.

ScHerdeR expects a normalized Seurat object and is designed to help check for batch effects or sample-specific dropout in clusters. It writes out:
  • A Seurat object with neighbors, clustering, and UMAP results (`seurat_obj_PC<pc>_res<res>.rds`)  
  • UMAP plots grouped and split by samples (`umap.png`, `umap_grouped_by_sample.png`, `umap_split_by_sample.png`)  
  • A table summarizing cell counts and percentages by sample and cluster (`cell_count_by_sample.rds`, `.xlsx`)  
  • An optional bar plot of cell percentages by sample and cluster (`cell_count_by_sample.png`)  

Usage:
  Rscript ScHerdeR.R \\
    --seurat_obj <Path to a Seurat .rds|.RData file> \\
    --output_dir <output_directory> (default: current directory) \\
    --npc <number_of_PCs> (default: 30) \\
    --resolution <clustering_resolution> (default: 0.6) \\
    --reduction <reduction_method> (default: "pca") \\
    --sample_ident <metadata_column_for_sample_ID> (default: orig.ident) \\
    --perc_calculation <TRUE/FALSE> (default: FALSE) \\
    --cols <Path to color vector in .rds| .RData format> (optional)
'

option_list <- list(
  make_option(
    c("--seurat_obj"),
    type = "character",
    default = NULL,
    help = "Path to the Seurat object file (.rds or .RData)",
    metavar = "character"
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = ".",
    help = "Directory to save output files (default: current directory)",
    metavar = "character"
  ),
  make_option(
    c("--npc"),
    type = "integer",
    default = 30,
    help = "Number of principal components to use (default: 30)",
    metavar = "integer"
  ),
  make_option(
    c("--resolution"),
    type = "numeric",
    default = 0.6,
    help = "Clustering resolution (default: 0.6)",
    metavar = "numeric"
  ),
  make_option(
    c("--reduction"),
    type = "character",
    default = "pca",
    help = "Dimensionality reduction method (default: 'pca')",
    metavar = "character"
  ),
  make_option(
    c("--sample_ident"),
    type = "character",
    default = "orig.ident",
    help = "Metadata column for sample identity (default: 'orig.ident')",
    metavar = "character"
  ),
  make_option(
    c("--perc_calculation"),
    type = "logical",
    default = FALSE,
    help = "Calculate and output cell percentage by sample and cluster (default: FALSE)",
    metavar = "logical"
  ),
  make_option(
    c("--cols"),
    type = "character",
    default = NULL,
    help = "Path to a color vector file (.rds or .RData) for UMAP plots (optional)",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ===== check the input parameters ===== ##
if (is.null(opt$seurat_obj)) {
  print_help(opt_parser)
  stop("Please provide a Seurat object file using --seurat_obj option.")
}
if (!file.exists(opt$seurat_obj)) {
  stop("The specified Seurat object file does not exist: ", opt$seurat_obj)
}

ScHerdeR(
  seurat_obj = opt$seurat_obj,
  output_dir = opt$output_dir,
  npc = opt$npc,
  resolution = opt$resolution,
  reduction = opt$reduction,
  sample_ident = opt$sample_ident,
  perc_calculation = opt$perc_calculation,
  cols = opt$cols
)