############################################################
# Single-sample spatial transcriptomics pipeline for NOD_1
# REVISED FIGURE 8 VERSION: run this script from the top.
# Before running, verify root_dir below points to the folder containing NOD_1.
# Focus: spatial context within one disease-model lacrimal gland section
#
# Recommended sample:
#   Folder: NOD_1
#   Analysis label: Case_1
#   Biological label: NOD.H-2b_1
#
# What this script does:
# 1. Loads one 10x Visium sample from NOD_1 / processed data
# 2. Performs QC and optional spot filtering
# 3. Normalizes with SCTransform
# 4. Runs PCA, UMAP, graph-based clustering
# 5. Maps clusters back to tissue space
# 6. Finds cluster marker genes within this single section
# 7. Plots biologically selected marker genes
# 8. Scores biological gene programs
# 9. Finds spatially variable genes
# 10. Fits spatial gradients with GAM models
# 11. Performs simple boundary and neighborhood analysis
# 12. Runs optional GO enrichment for cluster markers
#
# Important:
# This is a single-sample spatial-context analysis.
# It should NOT be described as Case versus Control DEG.
############################################################


############################################################
# Module 0. Editable settings
############################################################

rm(list = ls())
gc()

# EDIT THIS PATH if your folder is different.
# This folder should contain BALBC_1, BALBC_2, NOD_1, NOD_2.
root_dir <- "C:/Users/sawro/OneDrive/Desktop/spatial transcriptomics"

# Main sample selected for one-sample downstream analysis.
sample_folder <- "NOD_1"
sample_id <- "Case_1"
original_sample_id <- "NOD.H-2b_1"
condition_label <- "Case"

# Output folder.
out_dir <- file.path(root_dir, "NOD_1_single_sample_spatial_context")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis settings.
set.seed(1234)
APPLY_QC_FILTER <- TRUE

# Lenient filtering because your earlier QC showed good mitochondrial signal.
min_features <- 200
min_counts <- 500
max_percent_mt <- 5

# PCA and clustering settings.
npcs_run <- 50
dims_use <- 1:30
cluster_resolution <- 0.45

# Spatially variable gene settings.
n_svg_features <- 2000
svg_method_preferred <- "moransi"       # faster when supported
svg_method_fallback <- "markvariogram"  # slower but commonly available

# Figure settings.
FIG_DPI <- 300
FIG_FONT_FAMILY <- "sans"
FIG_BASE_SIZE <- 14
SPATIAL_POINT_SIZE <- 2.2

# GO enrichment.
RUN_GO_ENRICHMENT <- TRUE

# This single-sample pipeline does not run SPOTlight, STdeconvolve, or Giotto cell-cell communication by default.
# Those require a reliable single-cell reference or additional curated cell-type labels.
RUN_DECONVOLUTION <- FALSE


############################################################
# Module 1. Install and load packages
############################################################

install_if_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing)
  }
}

cran_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "stringr",
  "ggrepel",
  "mgcv",
  "Matrix",
  "FNN",
  "viridis",
  "scales",
  "future"
)

install_if_missing_cran(cran_pkgs)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  # Explicit dplyr:: namespaces are used below to avoid Bioconductor function masking.
  library(tidyr)
  library(stringr)
  library(ggrepel)
  library(mgcv)
  library(Matrix)
  library(FNN)
  library(viridis)
  library(scales)
  library(future)
})

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

if (RUN_GO_ENRICHMENT) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot")
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(enrichplot)
  })
}


############################################################
# Module 2. Helper functions
############################################################

message_section <- function(txt) {
  message("\n############################################################")
  message(txt)
  message("############################################################")
}

save_plot <- function(filename, plot_obj, width = 8, height = 6, dpi = FIG_DPI) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

save_table <- function(filename, table_obj) {
  write.csv(
    table_obj,
    file = file.path(out_dir, filename),
    row.names = FALSE
  )
}

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("_$", "", x)
  x
}

k_format <- function(x) {
  dplyr::case_when(
    x >= 1e6 ~ paste0(round(x / 1e6, 1), "M"),
    x >= 1e3 ~ paste0(round(x / 1e3, 1), "k"),
    TRUE ~ as.character(x)
  )
}

manuscript_theme <- theme(
  text = element_text(family = FIG_FONT_FAMILY, size = FIG_BASE_SIZE),
  plot.title = element_text(family = FIG_FONT_FAMILY, size = 16, face = "bold", hjust = 0.5),
  axis.title = element_text(family = FIG_FONT_FAMILY, size = 13),
  axis.text = element_text(family = FIG_FONT_FAMILY, size = 11),
  legend.title = element_text(family = FIG_FONT_FAMILY, size = 11),
  legend.text = element_text(family = FIG_FONT_FAMILY, size = 10),
  strip.text = element_text(family = FIG_FONT_FAMILY, size = 12, face = "bold")
)

find_spaceranger_dir <- function(root_dir, folder_name) {
  direct_dir <- file.path(root_dir, folder_name)
  processed_dir <- file.path(root_dir, folder_name, "processed data")

  direct_h5 <- file.path(direct_dir, "filtered_feature_bc_matrix.h5")
  processed_h5 <- file.path(processed_dir, "filtered_feature_bc_matrix.h5")

  if (file.exists(processed_h5)) {
    return(processed_dir)
  }
  if (file.exists(direct_h5)) {
    return(direct_dir)
  }

  stop(
    "Could not find filtered_feature_bc_matrix.h5.\nChecked:\n",
    processed_h5, "\n",
    direct_h5
  )
}

get_assay_matrix_safe <- function(obj, assay = NULL, layer = "data") {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  mat <- tryCatch(
    {
      SeuratObject::LayerData(obj, assay = assay, layer = layer)
    },
    error = function(e1) {
      tryCatch(
        {
          Seurat::GetAssayData(obj, assay = assay, slot = layer)
        },
        error = function(e2) {
          stop("Could not extract assay matrix for assay = ", assay, " and layer or slot = ", layer)
        }
      )
    }
  )

  return(mat)
}

get_feature_values <- function(obj, feature, assay = NULL, layer = "data") {
  if (feature %in% colnames(obj@meta.data)) {
    return(as.numeric(obj@meta.data[[feature]]))
  }

  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  if (feature %in% rownames(obj)) {
    mat <- get_assay_matrix_safe(obj, assay = assay, layer = layer)
    return(as.numeric(mat[feature, colnames(obj)]))
  }

  warning("Feature not found: ", feature)
  return(rep(NA_real_, ncol(obj)))
}

get_spatial_coords <- function(obj) {
  coords <- tryCatch(
    {
      as.data.frame(GetTissueCoordinates(obj))
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(coords) || nrow(coords) == 0) {
    if (length(obj@images) > 0) {
      coords <- as.data.frame(obj@images[[1]]@coordinates)
    } else {
      stop("No spatial coordinates found.")
    }
  }

  if ("cell" %in% colnames(coords)) {
    coords$barcode <- as.character(coords$cell)
  } else if ("barcode" %in% colnames(coords)) {
    coords$barcode <- as.character(coords$barcode)
  } else {
    coords$barcode <- rownames(coords)
  }

  x_candidates <- c("imagecol", "col", "x", "pxl_col_in_fullres", "tissuecol", "array_col")
  y_candidates <- c("imagerow", "row", "y", "pxl_row_in_fullres", "tissuerow", "array_row")

  x_col <- intersect(x_candidates, colnames(coords))[1]
  y_col <- intersect(y_candidates, colnames(coords))[1]

  if (is.na(x_col) || is.na(y_col)) {
    numeric_cols <- colnames(coords)[vapply(coords, is.numeric, logical(1))]
    if (length(numeric_cols) >= 2) {
      x_col <- numeric_cols[1]
      y_col <- numeric_cols[2]
    } else {
      stop("Could not identify spatial x and y columns.")
    }
  }

  out <- data.frame(
    barcode = as.character(coords$barcode),
    x = as.numeric(coords[[x_col]]),
    y = as.numeric(coords[[y_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$x) & !is.na(out$y), ]
  out <- out[out$barcode %in% colnames(obj), ]
  rownames(out) <- out$barcode
  out <- out[colnames(obj), , drop = FALSE]

  return(out)
}

plot_spatial_feature_simple <- function(obj, feature, title = NULL, assay = NULL, layer = "data") {
  coords <- get_spatial_coords(obj)
  val <- get_feature_values(obj, feature, assay = assay, layer = layer)

  df <- data.frame(
    barcode = colnames(obj),
    x = coords[colnames(obj), "x"],
    y = coords[colnames(obj), "y"],
    value = val,
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = x, y = y, color = value)) +
    geom_point(size = 1.0, alpha = 0.9) +
    scale_color_viridis_c(option = "C", na.value = "gray90") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
    labs(title = ifelse(is.null(title), feature, title), color = "Expression")
}

select_pc_number <- function(obj, max_pc = 50) {
  stdev <- obj[["pca"]]@stdev
  pct <- stdev / sum(stdev) * 100
  cumu <- cumsum(pct)

  pc1 <- which(cumu > 90 & pct < 5)[1]
  pc2_candidates <- which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1)
  pc2 <- ifelse(length(pc2_candidates) > 0, max(pc2_candidates) + 1, 30)

  pc_use <- suppressWarnings(min(pc1, pc2, max_pc, na.rm = TRUE))
  if (!is.finite(pc_use) || pc_use < 10) {
    pc_use <- 30
  }

  return(seq_len(pc_use))
}


############################################################
# Module 3. Load NOD_1 Visium sample
############################################################

message_section("Loading NOD_1 Visium sample")

sample_dir <- find_spaceranger_dir(root_dir, sample_folder)
message("Using sample directory: ", sample_dir)

hires_file <- file.path(sample_dir, "spatial", "tissue_hires_image.png")
lowres_file <- file.path(sample_dir, "spatial", "tissue_lowres_image.png")

if (file.exists(hires_file)) {
  message("High-resolution tissue image found. Loading with Read10X_Image.")
  image_obj <- Read10X_Image(
    image.dir = file.path(sample_dir, "spatial"),
    image.name = "tissue_hires_image.png",
    assay = "Spatial"
  )

  obj <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample_id,
    image = image_obj,
    filter.matrix = TRUE
  )
} else {
  message("High-resolution image not found. Loading with default Load10X_Spatial behavior.")
  obj <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample_id,
    filter.matrix = TRUE
  )
}

obj$sample_id <- sample_id
obj$original_sample_id <- original_sample_id
obj$condition <- condition_label
obj$orig.ident <- sample_id

# Mouse mitochondrial genes usually begin with mt-. Human genes usually begin with MT-.
if (length(grep("^mt-", rownames(obj), value = TRUE)) > 0) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
} else if (length(grep("^MT-", rownames(obj), value = TRUE)) > 0) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
} else {
  obj$percent.mt <- 0
}

saveRDS(obj, file.path(out_dir, "01_NOD_1_raw_spatial_object.rds"))

message("Loaded spots: ", ncol(obj))
message("Loaded genes: ", nrow(obj))


############################################################
# Module 4. Quality control before filtering
############################################################

message_section("Quality control before filtering")

qc_summary_before <- data.frame(
  sample_id = sample_id,
  original_sample_id = original_sample_id,
  condition = condition_label,
  n_spots = ncol(obj),
  median_nCount_Spatial = median(obj$nCount_Spatial, na.rm = TRUE),
  median_nFeature_Spatial = median(obj$nFeature_Spatial, na.rm = TRUE),
  median_percent_mt = median(obj$percent.mt, na.rm = TRUE),
  mean_nCount_Spatial = mean(obj$nCount_Spatial, na.rm = TRUE),
  mean_nFeature_Spatial = mean(obj$nFeature_Spatial, na.rm = TRUE),
  mean_percent_mt = mean(obj$percent.mt, na.rm = TRUE)
)

save_table("QC_summary_before_filtering.csv", qc_summary_before)

p_qc_vln_before <- VlnPlot(
  obj,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) +
  plot_annotation(title = "NOD_1 QC before filtering") &
  manuscript_theme

p_count_spatial <- SpatialFeaturePlot(
  obj,
  features = "nCount_Spatial",
  pt.size.factor = SPATIAL_POINT_SIZE,
  image.scale = ifelse(file.exists(hires_file), "hires", "lowres")
) +
  ggtitle("Spatial distribution of total counts") &
  manuscript_theme

p_feature_spatial <- SpatialFeaturePlot(
  obj,
  features = "nFeature_Spatial",
  pt.size.factor = SPATIAL_POINT_SIZE,
  image.scale = ifelse(file.exists(hires_file), "hires", "lowres")
) +
  ggtitle("Spatial distribution of detected genes") &
  manuscript_theme

p_mt_spatial <- SpatialFeaturePlot(
  obj,
  features = "percent.mt",
  pt.size.factor = SPATIAL_POINT_SIZE,
  image.scale = ifelse(file.exists(hires_file), "hires", "lowres")
) +
  ggtitle("Spatial distribution of mitochondrial percent") &
  manuscript_theme

save_plot("QC_violin_before_filtering.png", p_qc_vln_before, width = 12, height = 4)
save_plot("QC_spatial_counts_before_filtering.png", p_count_spatial, width = 7, height = 6)
save_plot("QC_spatial_features_before_filtering.png", p_feature_spatial, width = 7, height = 6)
save_plot("QC_spatial_percent_mt_before_filtering.png", p_mt_spatial, width = 7, height = 6)


############################################################
# Module 5. Optional spot filtering
############################################################

message_section("Optional spot filtering")

obj$pass_qc <- obj$nFeature_Spatial >= min_features &
  obj$nCount_Spatial >= min_counts &
  obj$percent.mt <= max_percent_mt

qc_filter_table <- data.frame(
  filter = c("Total spots before filtering", "Spots passing filter", "Spots removed"),
  n_spots = c(ncol(obj), sum(obj$pass_qc), ncol(obj) - sum(obj$pass_qc))
)

save_table("QC_filtering_summary.csv", qc_filter_table)

if (APPLY_QC_FILTER) {
  obj <- subset(
    obj,
    subset = nFeature_Spatial >= min_features &
      nCount_Spatial >= min_counts &
      percent.mt <= max_percent_mt
  )
}

qc_summary_after <- data.frame(
  sample_id = sample_id,
  original_sample_id = original_sample_id,
  condition = condition_label,
  n_spots = ncol(obj),
  median_nCount_Spatial = median(obj$nCount_Spatial, na.rm = TRUE),
  median_nFeature_Spatial = median(obj$nFeature_Spatial, na.rm = TRUE),
  median_percent_mt = median(obj$percent.mt, na.rm = TRUE),
  mean_nCount_Spatial = mean(obj$nCount_Spatial, na.rm = TRUE),
  mean_nFeature_Spatial = mean(obj$nFeature_Spatial, na.rm = TRUE),
  mean_percent_mt = mean(obj$percent.mt, na.rm = TRUE)
)

save_table("QC_summary_after_filtering.csv", qc_summary_after)

p_qc_vln_after <- VlnPlot(
  obj,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) +
  plot_annotation(title = "NOD_1 QC after filtering") &
  manuscript_theme

save_plot("QC_violin_after_filtering.png", p_qc_vln_after, width = 12, height = 4)
saveRDS(obj, file.path(out_dir, "02_NOD_1_after_QC.rds"))


############################################################
# Module 6. SCTransform normalization
############################################################

message_section("SCTransform normalization")

DefaultAssay(obj) <- "Spatial"

obj <- SCTransform(
  obj,
  assay = "Spatial",
  verbose = FALSE,
  return.only.var.genes = FALSE
)

saveRDS(obj, file.path(out_dir, "03_NOD_1_SCT_normalized.rds"))


############################################################
# Module 7. PCA, UMAP, and clustering within NOD_1
############################################################

message_section("PCA, UMAP, and clustering")

DefaultAssay(obj) <- "SCT"

obj <- RunPCA(obj, npcs = npcs_run, verbose = FALSE)

dims_auto <- select_pc_number(obj, max_pc = npcs_run)
message("Automatically suggested PCs: 1:", max(dims_auto))

# Use fixed PCs for reproducibility with your previous analysis unless you want auto.
dims_use <- seq_len(min(max(dims_use), max(dims_auto), npcs_run))
message("Using PCs: ", paste(range(dims_use), collapse = " to "))

p_elbow <- ElbowPlot(obj, ndims = npcs_run) +
  geom_vline(xintercept = max(dims_use), color = "red", linetype = "dashed") +
  ggtitle("PCA elbow plot for NOD_1") +
  manuscript_theme

save_plot("PCA_elbow_NOD_1.png", p_elbow, width = 7, height = 5)

obj <- RunUMAP(obj, reduction = "pca", dims = dims_use, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)

Idents(obj) <- "seurat_clusters"

cluster_counts <- as.data.frame(table(obj$seurat_clusters))
colnames(cluster_counts) <- c("cluster", "n_spots")
cluster_counts$percent_spots <- round(cluster_counts$n_spots / sum(cluster_counts$n_spots) * 100, 2)

save_table("cluster_spot_counts_NOD_1.csv", cluster_counts)

p_umap_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("NOD_1 UMAP by Seurat cluster") +
  theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
  manuscript_theme

save_plot("UMAP_clusters_NOD_1.png", p_umap_cluster, width = 7, height = 6)

p_spatial_cluster <- SpatialDimPlot(
  obj,
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size.factor = SPATIAL_POINT_SIZE,
  image.scale = ifelse(file.exists(hires_file), "hires", "lowres")
) +
  ggtitle("NOD_1 spatial clusters") &
  manuscript_theme

save_plot("Spatial_clusters_NOD_1.png", p_spatial_cluster, width = 7, height = 7)

saveRDS(obj, file.path(out_dir, "04_NOD_1_clustered.rds"))


############################################################
# Module 8. Cluster marker genes within one sample
############################################################

message_section("Finding cluster marker genes")

DefaultAssay(obj) <- "SCT"
Idents(obj) <- "seurat_clusters"

obj_prepared <- tryCatch(
  {
    PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
  },
  error = function(e) {
    message("PrepSCTFindMarkers failed. Continuing without recorrection.")
    message(e$message)
    obj
  }
)
obj <- obj_prepared

cluster_markers <- FindAllMarkers(
  obj,
  assay = "SCT",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

cluster_markers <- cluster_markers %>%
  dplyr::arrange(cluster, p_val_adj, desc(avg_log2FC))

save_table("cluster_markers_all_NOD_1.csv", cluster_markers)

top10_markers <- cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

save_table("cluster_markers_top10_NOD_1.csv", top10_markers)

# Heatmap can fail if many marker genes are not present in scale.data. This block handles that safely.
heatmap_features <- unique(top10_markers$gene)
heatmap_features <- intersect(heatmap_features, rownames(obj))

p_heatmap <- tryCatch(
  {
    DoHeatmap(obj, features = heatmap_features, group.by = "seurat_clusters") +
      NoLegend() +
      ggtitle("Top cluster markers in NOD_1") +
      manuscript_theme
  },
  error = function(e) {
    message("DoHeatmap failed: ", e$message)
    NULL
  }
)

if (!is.null(p_heatmap)) {
  save_plot("Cluster_marker_heatmap_NOD_1.png", p_heatmap, width = 12, height = 14)
}


############################################################
# Module 9. Biological gene panels for lacrimal gland spatial context
############################################################

message_section("Defining biological gene panels")

gene_panels <- list(
  epithelial_secretory = c("Scgb2b17", "Ltf", "Aqp5", "Krt8", "Krt18", "Krt19", "Muc1", "Pip"),
  antigen_presentation = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "B2m"),
  b_cell_plasma_cell = c("Cd79a", "Ms4a1", "Jchain", "Mzb1", "Ighm", "Igha"),
  t_cell = c("Cd3d", "Cd3g", "Cd4", "Cd8a", "Trac"),
  myeloid_stromal = c("Gsn", "Lyz2", "Cd68", "C1qa", "C1qb", "Csf1r", "Col1a1", "Col3a1"),
  interferon_chemokine = c("Cxcl9", "Cxcl10", "Ccl5", "Ifit1", "Ifit3", "Isg15", "Stat1", "Irf1", "Irf7"),
  metabolism_core = c("Hk1", "Hk2", "Gpi1", "Pfkp", "Aldoa", "Gapdh", "Pgk1", "Eno1", "Pkm", "Ldha", "Pdha1", "Cs", "Aco2", "Idh3a", "Sdha", "Mdh2", "Ndufa1", "Ndufb8", "Uqcrc1", "Cox4i1", "Cox5a", "Atp5f1a", "Atp5f1b"),
  lipid_metabolism = c("Acaca", "Fasn", "Scd1", "Cpt1a", "Cpt2", "Acadm", "Ppara", "Pparg", "Fabp4", "Fabp5"),
  oxidative_stress = c("Sod1", "Sod2", "Gpx1", "Gpx3", "Cat", "Prdx1", "Prdx2", "Hmox1", "Nqo1")
)

panel_presence <- dplyr::bind_rows(lapply(names(gene_panels), function(panel_name) {
  data.frame(
    panel = panel_name,
    gene = gene_panels[[panel_name]],
    present = gene_panels[[panel_name]] %in% rownames(obj),
    stringsAsFactors = FALSE
  )
}))

save_table("gene_panel_presence_NOD_1.csv", panel_presence)

# Add module scores. AddModuleScore creates columns named panel1.
for (panel_name in names(gene_panels)) {
  genes_present <- intersect(gene_panels[[panel_name]], rownames(obj))
  if (length(genes_present) >= 2) {
    obj <- AddModuleScore(
      obj,
      features = list(genes_present),
      name = paste0(panel_name, "_score"),
      assay = "SCT",
      ctrl = min(50, length(genes_present))
    )
  } else {
    message("Skipping module score for ", panel_name, ": fewer than 2 genes present.")
  }
}

module_score_cols <- grep("_score1$", colnames(obj@meta.data), value = TRUE)
module_score_table <- data.frame(
  module_score_column = module_score_cols,
  biological_program = gsub("_score1$", "", module_score_cols),
  stringsAsFactors = FALSE
)

save_table("module_score_columns_NOD_1.csv", module_score_table)

saveRDS(obj, file.path(out_dir, "05_NOD_1_with_module_scores.rds"))


############################################################
# Module 10. Marker gene and module-score plots
############################################################

message_section("Plotting marker genes and module scores")

core_marker_genes <- c(
  "Scgb2b17", "Ltf", "Jchain", "Cd79a", "Cd3g", "Gsn",
  "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1",
  "Cxcl9", "Cxcl10", "Isg15", "Stat1",
  "Ldha", "Sod2"
)

core_marker_genes_present <- intersect(core_marker_genes, rownames(obj))

if (length(core_marker_genes_present) > 0) {
  p_feature_umap <- FeaturePlot(
    obj,
    features = core_marker_genes_present,
    reduction = "umap",
    ncol = 4,
    order = TRUE,
    raster = FALSE
  ) &
    manuscript_theme

  save_plot("UMAP_core_marker_genes_NOD_1.png", p_feature_umap, width = 14, height = 12)

  p_dot_cluster <- DotPlot(
    obj,
    features = core_marker_genes_present,
    group.by = "seurat_clusters",
    assay = "SCT"
  ) +
    RotatedAxis() +
    ggtitle("Core marker genes by NOD_1 cluster") +
    manuscript_theme

  save_plot("DotPlot_core_marker_genes_by_cluster_NOD_1.png", p_dot_cluster, width = 12, height = 5)

  # Spatial plots in chunks so the text and legends remain readable.
  marker_chunks <- split(core_marker_genes_present, ceiling(seq_along(core_marker_genes_present) / 6))

  for (i in seq_along(marker_chunks)) {
    p_spatial_markers <- SpatialFeaturePlot(
      obj,
      features = marker_chunks[[i]],
      ncol = 2,
      pt.size.factor = SPATIAL_POINT_SIZE,
      image.scale = ifelse(file.exists(hires_file), "hires", "lowres"),
      alpha = c(0.1, 1)
    ) &
      manuscript_theme

    save_plot(
      paste0("Spatial_core_marker_genes_chunk_", i, "_NOD_1.png"),
      p_spatial_markers,
      width = 10,
      height = 12
    )
  }
}

if (length(module_score_cols) > 0) {
  p_module_dot <- DotPlot(
    obj,
    features = module_score_cols,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    scale_x_discrete(labels = gsub("_score1$", "", module_score_cols)) +
    ggtitle("Biological module scores by cluster") +
    manuscript_theme

  save_plot("DotPlot_module_scores_by_cluster_NOD_1.png", p_module_dot, width = 12, height = 5)

  module_chunks <- split(module_score_cols, ceiling(seq_along(module_score_cols) / 4))
  for (i in seq_along(module_chunks)) {
    plots_now <- lapply(module_chunks[[i]], function(score_col) {
      plot_spatial_feature_simple(
        obj,
        feature = score_col,
        title = gsub("_score1$", "", score_col),
        assay = "SCT"
      )
    })

    p_module_spatial <- wrap_plots(plots_now, ncol = 2)

    save_plot(
      paste0("Spatial_module_scores_chunk_", i, "_NOD_1.png"),
      p_module_spatial,
      width = 10,
      height = 10
    )
  }
}


############################################################
# Module 11. Spatially variable genes
############################################################

message_section("Spatially variable gene analysis")

DefaultAssay(obj) <- "SCT"

candidate_features <- VariableFeatures(obj)
if (length(candidate_features) < 100) {
  candidate_features <- rownames(obj)
}
candidate_features <- intersect(candidate_features, rownames(obj))
candidate_features <- candidate_features[seq_len(min(length(candidate_features), n_svg_features))]

obj_svg <- tryCatch(
  {
    FindSpatiallyVariableFeatures(
      obj,
      assay = "SCT",
      features = candidate_features,
      selection.method = svg_method_preferred,
      verbose = TRUE
    )
  },
  error = function(e1) {
    message("Preferred SVG method failed: ", svg_method_preferred)
    message(e1$message)
    message("Trying fallback method: ", svg_method_fallback)

    tryCatch(
      {
        FindSpatiallyVariableFeatures(
          obj,
          assay = "SCT",
          features = candidate_features,
          selection.method = svg_method_fallback,
          verbose = TRUE
        )
      },
      error = function(e2) {
        message("Fallback SVG method also failed: ", e2$message)
        NULL
      }
    )
  }
)

if (!is.null(obj_svg)) {
  obj <- obj_svg

  svg_genes <- tryCatch(
    {
      SpatiallyVariableFeatures(obj, assay = "SCT", selection.method = svg_method_preferred)
    },
    error = function(e) {
      tryCatch(
        {
          SpatiallyVariableFeatures(obj, assay = "SCT", selection.method = svg_method_fallback)
        },
        error = function(e2) {
          character(0)
        }
      )
    }
  )

  svg_table <- data.frame(
    rank = seq_along(svg_genes),
    gene = svg_genes,
    stringsAsFactors = FALSE
  )

  save_table("spatially_variable_genes_NOD_1.csv", svg_table)

  top_svg_genes <- head(svg_genes, 12)
  top_svg_genes <- intersect(top_svg_genes, rownames(obj))

  if (length(top_svg_genes) > 0) {
    svg_chunks <- split(top_svg_genes, ceiling(seq_along(top_svg_genes) / 6))

    for (i in seq_along(svg_chunks)) {
      p_svg <- SpatialFeaturePlot(
        obj,
        features = svg_chunks[[i]],
        ncol = 2,
        pt.size.factor = SPATIAL_POINT_SIZE,
        image.scale = ifelse(file.exists(hires_file), "hires", "lowres"),
        alpha = c(0.1, 1)
      ) &
        manuscript_theme

      save_plot(
        paste0("Spatial_top_SVG_chunk_", i, "_NOD_1.png"),
        p_svg,
        width = 10,
        height = 12
      )
    }
  }
} else {
  message("Spatially variable gene analysis did not produce results.")
}


############################################################
# Module 12. Spatial gradient modeling with GAM
############################################################

message_section("Spatial gradient analysis with GAM")

coords <- get_spatial_coords(obj)

gradient_features <- unique(c(
  "Cd74", "Scgb2b17", "H2-Eb1", "H2-Aa", "Jchain", "Cd79a", "Ltf", "Gsn",
  "Cxcl9", "Cxcl10", "Isg15", "Stat1", "Ldha", "Sod2",
  module_score_cols
))

gradient_features <- gradient_features[
  gradient_features %in% rownames(obj) | gradient_features %in% colnames(obj@meta.data)
]

fit_spatial_gradient <- function(obj, feature, coords, assay = "SCT", layer = "data") {
  val <- get_feature_values(obj, feature, assay = assay, layer = layer)

  df <- data.frame(
    barcode = colnames(obj),
    x = coords[colnames(obj), "x"],
    y = coords[colnames(obj), "y"],
    value = val,
    stringsAsFactors = FALSE
  )

  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$value), ]

  if (nrow(df) < 50 || length(unique(df$value)) < 5) {
    return(NULL)
  }

  k_use <- min(80, max(20, floor(nrow(df) / 20)))

  fit <- tryCatch(
    {
      mgcv::gam(value ~ s(x, y, k = k_use), data = df, method = "REML")
    },
    error = function(e) {
      message("GAM failed for ", feature, ": ", e$message)
      NULL
    }
  )

  if (is.null(fit)) {
    return(NULL)
  }

  df$fitted <- as.numeric(predict(fit, newdata = df))
  fit_sum <- summary(fit)

  out <- list(
    data = df,
    fit = fit,
    summary = data.frame(
      feature = feature,
      n_spots = nrow(df),
      deviance_explained_percent = round(fit_sum$dev.expl * 100, 2),
      adjusted_r_squared = round(fit_sum$r.sq, 4),
      stringsAsFactors = FALSE
    )
  )

  return(out)
}

gradient_results <- list()

for (feature in gradient_features) {
  message("Fitting spatial gradient for: ", feature)
  res <- fit_spatial_gradient(obj, feature, coords)
  if (!is.null(res)) {
    gradient_results[[feature]] <- res
  }
}

if (length(gradient_results) > 0) {
  gradient_summary <- dplyr::bind_rows(lapply(gradient_results, function(x) x$summary)) %>%
    dplyr::arrange(desc(deviance_explained_percent))

  save_table("GAM_spatial_gradient_summary_NOD_1.csv", gradient_summary)

  top_gradient_features <- head(gradient_summary$feature, 12)

  for (feature in top_gradient_features) {
    df <- gradient_results[[feature]]$data

    p_raw <- ggplot(df, aes(x = x, y = y, color = value)) +
      geom_point(size = 1.0, alpha = 0.9) +
      scale_color_viridis_c(option = "C", na.value = "gray90") +
      scale_y_reverse() +
      coord_fixed() +
      theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
      labs(title = paste0(feature, " raw expression"), color = "Expression")

    p_fit <- ggplot(df, aes(x = x, y = y, color = fitted)) +
      geom_point(size = 1.0, alpha = 0.9) +
      scale_color_viridis_c(option = "C", na.value = "gray90") +
      scale_y_reverse() +
      coord_fixed() +
      theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
      labs(title = paste0(feature, " fitted spatial gradient"), color = "GAM fitted")

    p_pair <- p_raw | p_fit

    save_plot(
      paste0("GAM_gradient_", clean_filename(feature), "_NOD_1.png"),
      p_pair,
      width = 10,
      height = 5
    )
  }
}


############################################################
# Module 13. Boundary and neighborhood analysis
# Revised: native Visium lattice adjacency + sensitivity analyses
############################################################

message_section("Boundary and neighborhood analysis: native Visium adjacency")

# ----------------------------------------------------------
# 13A. Read native Space Ranger array coordinates
# ----------------------------------------------------------
read_visium_positions <- function(sample_dir) {
  pos_csv <- file.path(sample_dir, "spatial", "tissue_positions.csv")
  pos_old <- file.path(sample_dir, "spatial", "tissue_positions_list.csv")

  if (file.exists(pos_csv)) {
    pos <- read.csv(pos_csv, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (file.exists(pos_old)) {
    pos <- read.csv(pos_old, header = FALSE, stringsAsFactors = FALSE)
    colnames(pos)[1:6] <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
  } else {
    stop("Could not find tissue_positions.csv or tissue_positions_list.csv in: ",
         file.path(sample_dir, "spatial"))
  }

  required_cols <- c("barcode", "array_row", "array_col")
  if (!all(required_cols %in% colnames(pos))) {
    stop("Tissue-position file is missing required columns: ",
         paste(setdiff(required_cols, colnames(pos)), collapse = ", "))
  }

  pos$barcode <- as.character(pos$barcode)
  pos$array_row <- as.integer(pos$array_row)
  pos$array_col <- as.integer(pos$array_col)
  pos
}

coords <- get_spatial_coords(obj)
cluster_vec <- as.character(obj$seurat_clusters)
names(cluster_vec) <- colnames(obj)

positions <- read_visium_positions(sample_dir)
positions <- positions[match(colnames(obj), positions$barcode), , drop = FALSE]

if (anyNA(positions$barcode)) {
  stop("Some retained Seurat spots could not be matched to Space Ranger tissue positions.")
}

boundary_df <- data.frame(
  barcode = colnames(obj),
  x = coords[colnames(obj), "x"],
  y = coords[colnames(obj), "y"],
  array_row = positions$array_row,
  array_col = positions$array_col,
  cluster = cluster_vec[colnames(obj)],
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 13B. Construct exact native Visium adjacency
# ----------------------------------------------------------
# Standard Visium uses an "orange-crate" hexagonal lattice.
# Immediate array neighbors are:
#   same row: (row, col - 2), (row, col + 2)
#   adjacent rows: (row +/- 1, col +/- 1)
# Crucially, only spots that remain in the current filtered object are connected.
# Missing/non-tissue/QC-removed spots are NOT replaced by farther spots.

build_native_visium_edges <- function(df) {
  offsets <- data.frame(
    dr = c(0, 0, -1, -1, 1, 1),
    dc = c(-2, 2, -1, 1, -1, 1)
  )

  key <- paste(df$array_row, df$array_col, sep = "_")
  index_by_key <- setNames(seq_len(nrow(df)), key)

  edge_list <- lapply(seq_len(nrow(df)), function(i) {
    candidate_keys <- paste(
      df$array_row[i] + offsets$dr,
      df$array_col[i] + offsets$dc,
      sep = "_"
    )
    j <- unname(index_by_key[candidate_keys])
    j <- j[!is.na(j)]
    if (length(j) == 0) return(NULL)
    data.frame(i = i, j = j)
  })

  edges <- dplyr::bind_rows(edge_list)
  if (nrow(edges) == 0) stop("No native Visium adjacency edges were found.")

  edges %>%
    dplyr::mutate(i_min = pmin(i, j), j_max = pmax(i, j)) %>%
    transmute(i = i_min, j = j_max) %>%
    dplyr::filter(i != j) %>%
    dplyr::distinct()
}

compute_graph_metrics <- function(cluster_labels, edges, n_spots) {
  stopifnot(length(cluster_labels) == n_spots)

  degree <- tabulate(c(edges$i, edges$j), nbins = n_spots)
  is_diff_edge <- cluster_labels[edges$i] != cluster_labels[edges$j]
  n_diff <- tabulate(
    c(edges$i[is_diff_edge], edges$j[is_diff_edge]),
    nbins = n_spots
  )

  neighbor_diff_fraction <- ifelse(degree > 0, n_diff / degree, NA_real_)

  # "Any-different" is a graph-theoretic interface label and is kept as a
  # descriptive secondary measure, not the primary coherence statistic.
  is_boundary_any <- ifelse(degree > 0, n_diff > 0, NA)

  # Stricter sensitivity definition: at least half of available native
  # neighbors have a different cluster label.
  is_boundary_majority <- ifelse(
    degree > 0,
    neighbor_diff_fraction >= 0.5,
    NA
  )

  spot_table <- data.frame(
    degree = degree,
    n_different_neighbors = n_diff,
    neighbor_diff_fraction = neighbor_diff_fraction,
    is_boundary_any = is_boundary_any,
    is_boundary_majority = is_boundary_majority
  )

  overall <- data.frame(
    n_spots = n_spots,
    n_isolated_spots = sum(degree == 0),
    mean_degree = mean(degree),
    mean_neighbor_diff_fraction = mean(neighbor_diff_fraction, na.rm = TRUE),
    boundary_any_fraction = mean(is_boundary_any, na.rm = TRUE),
    boundary_majority_fraction = mean(is_boundary_majority, na.rm = TRUE)
  )

  list(spot_table = spot_table, overall = overall)
}

native_edges <- build_native_visium_edges(boundary_df)
native_metrics <- compute_graph_metrics(
  cluster_labels = boundary_df$cluster,
  edges = native_edges,
  n_spots = nrow(boundary_df)
)

boundary_df <- bind_cols(boundary_df, native_metrics$spot_table)

save_table("boundary_spot_table_native_adjacency_NOD_1.csv", boundary_df)
save_table("boundary_summary_overall_native_adjacency_NOD_1.csv", native_metrics$overall)

boundary_summary <- boundary_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_spots = n(),
    mean_native_degree = round(mean(degree), 3),
    mean_neighbor_diff_fraction = round(mean(neighbor_diff_fraction, na.rm = TRUE), 4),
    boundary_any_fraction = round(mean(is_boundary_any, na.rm = TRUE), 4),
    boundary_majority_fraction = round(mean(is_boundary_majority, na.rm = TRUE), 4),
    .groups = "drop"
  )

save_table("boundary_summary_by_cluster_native_adjacency_NOD_1.csv", boundary_summary)

# Primary Figure 8 map: continuous local mixing measure.
p_boundary_native <- ggplot(
  boundary_df,
  aes(x = x, y = y, color = cluster, alpha = neighbor_diff_fraction)
) +
  geom_point(size = 1.1) +
  scale_alpha_continuous(range = c(0.15, 1), na.value = 0.15) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
  labs(
    title = "NOD_1 local cluster coherence: native Visium adjacency",
    color = "Cluster",
    alpha = "Neighbor\ndifference"
  )

save_plot("Boundary_map_native_adjacency_NOD_1.png", p_boundary_native, width = 8, height = 7)

# Native neighboring-cluster adjacency table.
adjacency_table_native <- data.frame(
  cluster_a = boundary_df$cluster[native_edges$i],
  cluster_b = boundary_df$cluster[native_edges$j],
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  dplyr::mutate(
    pair_a = sort(c(cluster_a, cluster_b))[1],
    pair_b = sort(c(cluster_a, cluster_b))[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::count(pair_a, pair_b, name = "native_neighbor_pair_count") %>%
  dplyr::arrange(desc(native_neighbor_pair_count))

save_table("cluster_neighbor_adjacency_native_NOD_1.csv", adjacency_table_native)

# ----------------------------------------------------------
# 13C. Sensitivity analysis 1: distance-limited graph
# ----------------------------------------------------------
# Estimate the physical first-neighbor spacing from the native edges, then
# retain only coordinate-neighbor links within 1.25x that spacing. This prevents
# kNN from bridging tissue holes or disconnected fragments.

native_edge_distance <- sqrt(
  (boundary_df$x[native_edges$i] - boundary_df$x[native_edges$j])^2 +
  (boundary_df$y[native_edges$i] - boundary_df$y[native_edges$j])^2
)

native_spacing <- median(native_edge_distance, na.rm = TRUE)
radius_multiplier <- 1.25
radius_cutoff <- radius_multiplier * native_spacing

k_radius <- min(12, nrow(boundary_df) - 1)
knn_radius <- FNN::get.knn(
  as.matrix(boundary_df[, c("x", "y")]),
  k = k_radius
)

radius_edge_list <- lapply(seq_len(nrow(boundary_df)), function(i) {
  keep <- which(knn_radius$nn.dist[i, ] <= radius_cutoff)
  if (length(keep) == 0) return(NULL)
  data.frame(i = i, j = knn_radius$nn.index[i, keep])
})

radius_edges <- dplyr::bind_rows(radius_edge_list) %>%
  dplyr::mutate(i_min = pmin(i, j), j_max = pmax(i, j)) %>%
  transmute(i = i_min, j = j_max) %>%
  dplyr::filter(i != j) %>%
  dplyr::distinct()

radius_metrics <- compute_graph_metrics(
  cluster_labels = boundary_df$cluster,
  edges = radius_edges,
  n_spots = nrow(boundary_df)
)

graph_sensitivity <- dplyr::bind_rows(
  native_metrics$overall %>% dplyr::mutate(graph = "native_array_adjacency"),
  radius_metrics$overall %>% dplyr::mutate(graph = "distance_limited")
) %>%
  dplyr::mutate(
    native_spacing_pixels = native_spacing,
    radius_multiplier = ifelse(graph == "distance_limited", radius_multiplier, NA_real_),
    radius_cutoff_pixels = ifelse(graph == "distance_limited", radius_cutoff, NA_real_)
  ) %>%
  dplyr::select(graph, everything())

native_edge_key <- paste(native_edges$i, native_edges$j, sep = "_")
radius_edge_key <- paste(radius_edges$i, radius_edges$j, sep = "_")
edge_jaccard <- length(intersect(native_edge_key, radius_edge_key)) /
  length(union(native_edge_key, radius_edge_key))

spot_graph_agreement <- data.frame(
  spearman_neighbor_diff = suppressWarnings(cor(
    native_metrics$spot_table$neighbor_diff_fraction,
    radius_metrics$spot_table$neighbor_diff_fraction,
    use = "complete.obs",
    method = "spearman"
  )),
  boundary_any_agreement = mean(
    native_metrics$spot_table$is_boundary_any == radius_metrics$spot_table$is_boundary_any,
    na.rm = TRUE
  ),
  edge_jaccard = edge_jaccard
)

save_table("Figure8_graph_sensitivity_NOD_1.csv", graph_sensitivity)
save_table("Figure8_graph_agreement_NOD_1.csv", spot_graph_agreement)

# ----------------------------------------------------------
# 13D. Sensitivity analysis 2: clustering resolution
# ----------------------------------------------------------
# The manuscript's primary clustering resolution is 0.45. Bracket it with
# lower and higher resolutions to assess whether local-coherence conclusions
# depend strongly on cluster granularity.

cluster_resolutions_sensitivity <- sort(unique(c(0.30, cluster_resolution, 0.60)))

resolution_sensitivity <- dplyr::bind_rows(lapply(cluster_resolutions_sensitivity, function(res) {
  obj_tmp <- FindClusters(obj, resolution = res, verbose = FALSE)
  cl_tmp <- as.character(obj_tmp$seurat_clusters)

  met_tmp <- compute_graph_metrics(
    cluster_labels = cl_tmp,
    edges = native_edges,
    n_spots = nrow(boundary_df)
  )

  data.frame(
    resolution = res,
    n_clusters = length(unique(cl_tmp)),
    mean_neighbor_diff_fraction = met_tmp$overall$mean_neighbor_diff_fraction,
    boundary_any_fraction = met_tmp$overall$boundary_any_fraction,
    boundary_majority_fraction = met_tmp$overall$boundary_majority_fraction
  )
}))

save_table("Figure8_cluster_resolution_sensitivity_NOD_1.csv", resolution_sensitivity)

p_res_mixing <- ggplot(
  resolution_sensitivity,
  aes(x = resolution, y = mean_neighbor_diff_fraction)
) +
  geom_line() +
  geom_point(size = 2.5) +
  theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
  labs(
    title = "Sensitivity to clustering resolution",
    x = "Seurat clustering resolution",
    y = "Mean neighbor-difference fraction"
  )

p_res_boundary <- ggplot(
  resolution_sensitivity,
  aes(x = resolution, y = boundary_any_fraction)
) +
  geom_line() +
  geom_point(size = 2.5) +
  theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
  labs(
    title = "Binary interface sensitivity",
    x = "Seurat clustering resolution",
    y = "Fraction with >=1 different native neighbor"
  )

p_resolution_sensitivity <- p_res_mixing | p_res_boundary
save_plot(
  "Figure8_cluster_resolution_sensitivity_NOD_1.png",
  p_resolution_sensitivity,
  width = 11,
  height = 5
)

message("Native adjacency analysis complete.")
message("Native edges: ", nrow(native_edges))
message("Distance-limited edges: ", nrow(radius_edges))
message("Graph edge Jaccard: ", round(edge_jaccard, 4))
message("See Figure8_graph_sensitivity_NOD_1.csv and Figure8_cluster_resolution_sensitivity_NOD_1.csv")


############################################################
# Module 14. Functional enrichment of cluster marker genes
############################################################

if (RUN_GO_ENRICHMENT) {

  message_section("GO Biological Process enrichment for cluster markers")

  run_go_for_gene_list <- function(gene_symbols, output_prefix, title_text) {
    gene_symbols <- unique(na.omit(gene_symbols))

    if (length(gene_symbols) < 10) {
      message("Skipping GO for ", output_prefix, ": fewer than 10 genes.")
      return(NULL)
    }

    gene_map <- tryCatch(
      {
        bitr(
          gene_symbols,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Mm.eg.db
        )
      },
      error = function(e) {
        message("Gene ID conversion failed for ", output_prefix, ": ", e$message)
        NULL
      }
    )

    if (is.null(gene_map) || nrow(gene_map) < 10) {
      message("Skipping GO for ", output_prefix, ": fewer than 10 mapped genes.")
      return(NULL)
    }

    ego <- enrichGO(
      gene = unique(gene_map$ENTREZID),
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )

    ego_df <- as.data.frame(ego)

    if (nrow(ego_df) == 0) {
      message("No enriched GO terms for ", output_prefix)
      return(NULL)
    }

    save_table(paste0(output_prefix, "_GO_BP.csv"), ego_df)

    top_n <- min(10, nrow(ego_df))
    ego_plot_df <- ego_df %>%
      dplyr::arrange(p.adjust) %>%
      slice_head(n = top_n) %>%
      dplyr::mutate(
        Description_wrapped = stringr::str_wrap(Description, width = 35),
        Description_wrapped = factor(Description_wrapped, levels = rev(Description_wrapped))
      )

    ego_plot_df$GeneRatio_numeric <- sapply(
      ego_plot_df$GeneRatio,
      function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }
    )

    p_go <- ggplot(
      ego_plot_df,
      aes(x = GeneRatio_numeric, y = Description_wrapped, size = Count, color = p.adjust)
    ) +
      geom_point(alpha = 0.9) +
      scale_color_viridis_c(option = "C", direction = -1) +
      scale_size(range = c(3, 9)) +
      theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
      manuscript_theme +
      labs(
        title = title_text,
        x = "Gene ratio",
        y = NULL,
        color = "Adjusted\np-value",
        size = "Gene\ncount"
      )

    save_plot(paste0(output_prefix, "_GO_BP_dotplot.png"), p_go, width = 9, height = 6)

    return(ego_df)
  }

  significant_cluster_markers <- cluster_markers %>%
    dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.25)

  cluster_ids <- sort(unique(as.character(significant_cluster_markers$cluster)))

  for (cl in cluster_ids) {
    genes_cl <- significant_cluster_markers %>%
      dplyr::filter(as.character(cluster) == cl) %>%
      dplyr::pull(gene)

    run_go_for_gene_list(
      gene_symbols = genes_cl,
      output_prefix = paste0("cluster_", clean_filename(cl), "_marker_genes"),
      title_text = paste0("GO BP enrichment: cluster ", cl, " marker genes")
    )
  }
}


############################################################
# Module 15. Save final object and output manifest
############################################################

message_section("Saving final outputs")

saveRDS(obj, file.path(out_dir, "final_NOD_1_single_sample_spatial_context.rds"))

output_manifest <- data.frame(
  output_group = c(
    "RDS",
    "QC",
    "Clustering",
    "Marker genes",
    "Biological panels",
    "Spatially variable genes",
    "Spatial gradients",
    "Boundary analysis",
    "GO enrichment"
  ),
  key_files = c(
    "final_NOD_1_single_sample_spatial_context.rds",
    "QC_summary_before_filtering.csv; QC_summary_after_filtering.csv; QC_violin_after_filtering.png",
    "UMAP_clusters_NOD_1.png; Spatial_clusters_NOD_1.png; cluster_spot_counts_NOD_1.csv",
    "cluster_markers_all_NOD_1.csv; cluster_markers_top10_NOD_1.csv",
    "gene_panel_presence_NOD_1.csv; DotPlot_core_marker_genes_by_cluster_NOD_1.png; Spatial_core_marker_genes_chunk_*.png",
    "spatially_variable_genes_NOD_1.csv; Spatial_top_SVG_chunk_*.png",
    "GAM_spatial_gradient_summary_NOD_1.csv; GAM_gradient_*.png",
    "boundary_summary_overall_native_adjacency_NOD_1.csv; boundary_summary_by_cluster_native_adjacency_NOD_1.csv; Boundary_map_native_adjacency_NOD_1.png; Figure8_graph_sensitivity_NOD_1.csv; Figure8_graph_agreement_NOD_1.csv; Figure8_cluster_resolution_sensitivity_NOD_1.csv; Figure8_cluster_resolution_sensitivity_NOD_1.png",
    "cluster_*_marker_genes_GO_BP.csv; cluster_*_marker_genes_GO_BP_dotplot.png"
  ),
  interpretation = c(
    "Final processed object for NOD_1 only",
    "Checks count depth, gene detection, and mitochondrial signal",
    "Finds transcriptional spot groups and maps them to tissue space",
    "Finds genes that define each cluster within this sample",
    "Connects clusters and tissue regions to epithelial, immune, stromal, metabolic, and oxidative-stress signals",
    "Finds genes whose expression is spatially organized",
    "Quantifies how much expression pattern is explained by x and y tissue location",
    "Describes where cluster domains meet or intermix",
    "Summarizes biological processes represented by cluster marker genes"
  ),
  stringsAsFactors = FALSE
)

save_table("output_manifest_NOD_1_single_sample_pipeline.csv", output_manifest)

message("Pipeline complete.")
message("All outputs saved in: ", out_dir)
message("Remember: this is a single-sample spatial-context analysis, not Case versus Control DEG.")


