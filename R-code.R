# NOD_1 single-sample Visium analysis
# Mauduit et al. (2022), GSE210380
# Run from a directory containing NOD_1, or set the NOD1_ROOT environment variable.

rm(list = ls())
gc()

config <- list(
  root_dir = Sys.getenv("NOD1_ROOT", unset = getwd()),
  sample_folder = "NOD_1",
  sample_id = "Case_1",
  original_sample_id = "NOD.H-2b_1",
  condition = "Case",
  output_folder = "NOD_1_single_sample_spatial_context",
  seed = 1234L,
  apply_qc_filter = TRUE,
  min_features = 200L,
  min_counts = 500L,
  max_percent_mt = 5,
  n_pcs = 50L,
  dims = 1:30,
  cluster_resolution = 0.45,
  svg_n_features = 2000L,
  svg_methods = c("moransi", "markvariogram"),
  figure_dpi = 300L,
  figure_base_size = 14,
  spatial_point_size = 2.2,
  run_go = TRUE,
  resolution_sensitivity = c(0.30, 0.45, 0.60),
  radius_multiplier = 1.25
)

config$root_dir <- normalizePath(
  config$root_dir,
  winslash = "/",
  mustWork = TRUE
)
output_dir <- file.path(config$root_dir, config$output_folder)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(config$seed)

cran_packages <- c(
  "Seurat", "SeuratObject", "ggplot2", "patchwork", "dplyr", "tidyr",
  "stringr", "ggrepel", "mgcv", "Matrix", "FNN", "viridis", "scales",
  "future"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_cran, collapse = ", "),
    ". Install them before running this script."
  )
}

if (config$run_go) {
  bioc_packages <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot")
  missing_bioc <- bioc_packages[
    !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_bioc) > 0) {
    stop(
      "Missing Bioconductor packages: ",
      paste(missing_bioc, collapse = ", "),
      ". Install them before running this script."
    )
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
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

if (config$run_go) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(enrichplot)
  })
}

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

theme_publication <- ggplot2::theme(
  text = ggplot2::element_text(size = config$figure_base_size),
  plot.title = ggplot2::element_text(
    size = 16,
    face = "bold",
    hjust = 0.5
  ),
  axis.title = ggplot2::element_text(size = 13),
  axis.text = ggplot2::element_text(size = 11),
  legend.title = ggplot2::element_text(size = 11),
  legend.text = ggplot2::element_text(size = 10),
  strip.text = ggplot2::element_text(size = 12, face = "bold")
)

log_section <- function(label) {
  message("\n", label)
  message(strrep("-", nchar(label)))
}

save_plot <- function(filename, plot, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = config$figure_dpi,
    bg = "white"
  )
}

save_table <- function(filename, x) {
  utils::write.csv(
    x,
    file = file.path(output_dir, filename),
    row.names = FALSE
  )
}

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("_$", "", x)
}

find_spaceranger_dir <- function(root_dir, sample_folder) {
  candidates <- c(
    file.path(root_dir, sample_folder, "processed data"),
    file.path(root_dir, sample_folder)
  )

  has_matrix <- file.exists(
    file.path(candidates, "filtered_feature_bc_matrix.h5")
  )

  if (!any(has_matrix)) {
    stop(
      "filtered_feature_bc_matrix.h5 was not found under ",
      file.path(root_dir, sample_folder)
    )
  }

  candidates[which(has_matrix)[1]]
}

get_assay_matrix <- function(object, assay = NULL, layer = "data") {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  tryCatch(
    SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer
    ),
    error = function(e) {
      Seurat::GetAssayData(
        object,
        assay = assay,
        slot = layer
      )
    }
  )
}

get_feature_values <- function(object, feature, assay = NULL, layer = "data") {
  if (feature %in% colnames(object@meta.data)) {
    return(as.numeric(object@meta.data[[feature]]))
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  if (!feature %in% rownames(object[[assay]])) {
    return(rep(NA_real_, ncol(object)))
  }

  matrix <- get_assay_matrix(object, assay = assay, layer = layer)
  as.numeric(matrix[feature, colnames(object)])
}

get_spatial_coordinates <- function(object) {
  coordinates <- tryCatch(
    as.data.frame(Seurat::GetTissueCoordinates(object)),
    error = function(e) NULL
  )

  if (is.null(coordinates) || nrow(coordinates) == 0) {
    if (length(object@images) == 0) {
      stop("No spatial coordinates are available in the Seurat object.")
    }
    coordinates <- as.data.frame(object@images[[1]]@coordinates)
  }

  if ("cell" %in% colnames(coordinates)) {
    coordinates$barcode <- as.character(coordinates$cell)
  } else if ("barcode" %in% colnames(coordinates)) {
    coordinates$barcode <- as.character(coordinates$barcode)
  } else {
    coordinates$barcode <- rownames(coordinates)
  }

  x_candidates <- c(
    "imagecol", "col", "x", "pxl_col_in_fullres", "tissuecol", "array_col"
  )
  y_candidates <- c(
    "imagerow", "row", "y", "pxl_row_in_fullres", "tissuerow", "array_row"
  )

  x_col <- intersect(x_candidates, colnames(coordinates))[1]
  y_col <- intersect(y_candidates, colnames(coordinates))[1]

  if (is.na(x_col) || is.na(y_col)) {
    numeric_cols <- colnames(coordinates)[
      vapply(coordinates, is.numeric, logical(1))
    ]

    if (length(numeric_cols) < 2) {
      stop("Unable to identify spatial x and y coordinates.")
    }

    x_col <- numeric_cols[1]
    y_col <- numeric_cols[2]
  }

  out <- data.frame(
    barcode = as.character(coordinates$barcode),
    x = as.numeric(coordinates[[x_col]]),
    y = as.numeric(coordinates[[y_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[
    is.finite(out$x) &
      is.finite(out$y) &
      out$barcode %in% colnames(object),
    ,
    drop = FALSE
  ]

  rownames(out) <- out$barcode
  missing_barcodes <- setdiff(colnames(object), rownames(out))

  if (length(missing_barcodes) > 0) {
    stop(
      "Spatial coordinates are missing for ",
      length(missing_barcodes),
      " retained spots."
    )
  }

  out[colnames(object), , drop = FALSE]
}

plot_spatial_feature <- function(
    object,
    feature,
    title = feature,
    assay = NULL,
    layer = "data") {

  coordinates <- get_spatial_coordinates(object)
  values <- get_feature_values(
    object,
    feature,
    assay = assay,
    layer = layer
  )

  data <- data.frame(
    x = coordinates[colnames(object), "x"],
    y = coordinates[colnames(object), "y"],
    value = values
  )

  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, color = value)) +
    ggplot2::geom_point(size = 1, alpha = 0.9) +
    viridis::scale_color_viridis(
      option = "C",
      na.value = "grey90"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = config$figure_base_size) +
    ggplot2::labs(title = title, color = "Expression")
}

select_pc_range <- function(object, max_pc = 50L) {
  stdev <- object[["pca"]]@stdev
  percent <- stdev / sum(stdev) * 100
  cumulative <- cumsum(percent)

  pc_cumulative <- which(cumulative > 90 & percent < 5)[1]
  drop_candidates <- which(
    (percent[-length(percent)] - percent[-1]) > 0.1
  )
  pc_drop <- if (length(drop_candidates) > 0) {
    max(drop_candidates) + 1L
  } else {
    30L
  }

  candidates <- c(pc_cumulative, pc_drop, max_pc)
  candidates <- candidates[is.finite(candidates)]
  n_pc <- min(candidates)

  if (!is.finite(n_pc) || n_pc < 10) {
    n_pc <- min(30L, max_pc)
  }

  seq_len(n_pc)
}

read_visium_positions <- function(sample_dir) {
  current_file <- file.path(
    sample_dir,
    "spatial",
    "tissue_positions.csv"
  )
  legacy_file <- file.path(
    sample_dir,
    "spatial",
    "tissue_positions_list.csv"
  )

  if (file.exists(current_file)) {
    positions <- utils::read.csv(
      current_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (file.exists(legacy_file)) {
    positions <- utils::read.csv(
      legacy_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(positions)[1:6] <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  } else {
    stop("No Space Ranger tissue-position file was found.")
  }

  required <- c("barcode", "array_row", "array_col")
  missing <- setdiff(required, colnames(positions))

  if (length(missing) > 0) {
    stop(
      "Missing tissue-position columns: ",
      paste(missing, collapse = ", ")
    )
  }

  positions$barcode <- as.character(positions$barcode)
  positions$array_row <- as.integer(positions$array_row)
  positions$array_col <- as.integer(positions$array_col)
  positions
}

build_native_visium_edges <- function(spots) {
  offsets <- data.frame(
    dr = c(0L, 0L, -1L, -1L, 1L, 1L),
    dc = c(-2L, 2L, -1L, 1L, -1L, 1L)
  )

  spot_key <- paste(spots$array_row, spots$array_col, sep = "_")
  index <- stats::setNames(seq_len(nrow(spots)), spot_key)

  edge_list <- lapply(seq_len(nrow(spots)), function(i) {
    neighbor_key <- paste(
      spots$array_row[i] + offsets$dr,
      spots$array_col[i] + offsets$dc,
      sep = "_"
    )

    j <- unname(index[neighbor_key])
    j <- j[!is.na(j)]

    if (length(j) == 0) {
      return(NULL)
    }

    data.frame(i = i, j = j)
  })

  edges <- dplyr::bind_rows(edge_list)

  if (nrow(edges) == 0) {
    stop("No native Visium adjacency edges were found.")
  }

  edges |>
    dplyr::mutate(
      i_min = pmin(i, j),
      j_max = pmax(i, j)
    ) |>
    dplyr::transmute(i = i_min, j = j_max) |>
    dplyr::filter(i != j) |>
    dplyr::distinct()
}

compute_graph_metrics <- function(cluster_labels, edges, n_spots) {
  if (length(cluster_labels) != n_spots) {
    stop("Cluster labels and spot count are inconsistent.")
  }

  if (nrow(edges) == 0) {
    stop("Graph contains no edges.")
  }

  degree <- tabulate(c(edges$i, edges$j), nbins = n_spots)
  different <- cluster_labels[edges$i] != cluster_labels[edges$j]

  n_different <- tabulate(
    c(edges$i[different], edges$j[different]),
    nbins = n_spots
  )

  neighbor_diff_fraction <- ifelse(
    degree > 0,
    n_different / degree,
    NA_real_
  )

  spot_table <- data.frame(
    degree = degree,
    n_different_neighbors = n_different,
    neighbor_diff_fraction = neighbor_diff_fraction,
    is_boundary_any = ifelse(
      degree > 0,
      n_different > 0,
      NA
    ),
    is_boundary_majority = ifelse(
      degree > 0,
      neighbor_diff_fraction >= 0.5,
      NA
    )
  )

  overall <- data.frame(
    n_spots = n_spots,
    n_isolated_spots = sum(degree == 0),
    mean_degree = mean(degree),
    mean_neighbor_diff_fraction = mean(
      neighbor_diff_fraction,
      na.rm = TRUE
    ),
    boundary_any_fraction = mean(
      spot_table$is_boundary_any,
      na.rm = TRUE
    ),
    boundary_majority_fraction = mean(
      spot_table$is_boundary_majority,
      na.rm = TRUE
    )
  )

  list(spot_table = spot_table, overall = overall)
}

fit_spatial_gradient <- function(
    object,
    feature,
    coordinates,
    assay = "SCT",
    layer = "data") {

  data <- data.frame(
    barcode = colnames(object),
    x = coordinates[colnames(object), "x"],
    y = coordinates[colnames(object), "y"],
    value = get_feature_values(
      object,
      feature,
      assay = assay,
      layer = layer
    ),
    stringsAsFactors = FALSE
  )

  data <- data[
    is.finite(data$x) &
      is.finite(data$y) &
      is.finite(data$value),
    ,
    drop = FALSE
  ]

  if (nrow(data) < 50 || length(unique(data$value)) < 5) {
    return(NULL)
  }

  k <- min(80L, max(20L, floor(nrow(data) / 20)))

  fit <- tryCatch(
    mgcv::gam(
      value ~ s(x, y, k = k),
      data = data,
      method = "REML"
    ),
    error = function(e) {
      message("GAM failed for ", feature, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(NULL)
  }

  fit_summary <- summary(fit)
  data$fitted <- as.numeric(
    stats::predict(fit, newdata = data)
  )

  list(
    data = data,
    fit = fit,
    summary = data.frame(
      feature = feature,
      n_spots = nrow(data),
      deviance_explained_percent = round(
        fit_summary$dev.expl * 100,
        2
      ),
      adjusted_r_squared = round(fit_summary$r.sq, 4)
    )
  )
}

run_spatial_variable_analysis <- function(object, features) {
  for (method in config$svg_methods) {
    result <- tryCatch(
      Seurat::FindSpatiallyVariableFeatures(
        object,
        assay = "SCT",
        features = features,
        selection.method = method,
        verbose = TRUE
      ),
      error = function(e) {
        message(
          "Spatial-variable analysis failed with ",
          method,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(result)) {
      genes <- tryCatch(
        Seurat::SpatiallyVariableFeatures(
          result,
          assay = "SCT",
          selection.method = method
        ),
        error = function(e) character(0)
      )

      return(list(object = result, genes = genes, method = method))
    }
  }

  NULL
}

run_go_enrichment <- function(gene_symbols, output_prefix, title) {
  genes <- unique(stats::na.omit(gene_symbols))

  if (length(genes) < 10) {
    return(NULL)
  }

  gene_map <- tryCatch(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    ),
    error = function(e) NULL
  )

  if (is.null(gene_map) || nrow(gene_map) < 10) {
    return(NULL)
  }

  enrichment <- clusterProfiler::enrichGO(
    gene = unique(gene_map$ENTREZID),
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  results <- as.data.frame(enrichment)

  if (nrow(results) == 0) {
    return(NULL)
  }

  save_table(
    paste0(output_prefix, "_GO_BP.csv"),
    results
  )

  plot_data <- results |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = min(10, nrow(results))) |>
    dplyr::mutate(
      Description_wrapped = stringr::str_wrap(
        Description,
        width = 35
      ),
      Description_wrapped = factor(
        Description_wrapped,
        levels = rev(Description_wrapped)
      ),
      GeneRatio_numeric = vapply(
        GeneRatio,
        function(x) {
          parts <- strsplit(x, "/", fixed = TRUE)[[1]]
          as.numeric(parts[1]) / as.numeric(parts[2])
        },
        numeric(1)
      )
    )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = Description_wrapped,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    viridis::scale_color_viridis(
      option = "C",
      direction = -1
    ) +
    ggplot2::scale_size(range = c(3, 9)) +
    ggplot2::theme_classic(
      base_size = config$figure_base_size
    ) +
    theme_publication +
    ggplot2::labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted\np-value",
      size = "Gene\ncount"
    )

  save_plot(
    paste0(output_prefix, "_GO_BP_dotplot.png"),
    plot,
    width = 9,
    height = 6
  )

  results
}


# Load data -------------------------------------------------------------------

log_section("Load Visium data")

sample_dir <- find_spaceranger_dir(
  config$root_dir,
  config$sample_folder
)

hires_file <- file.path(
  sample_dir,
  "spatial",
  "tissue_hires_image.png"
)
image_scale <- if (file.exists(hires_file)) "hires" else "lowres"

if (file.exists(hires_file)) {
  image <- Seurat::Read10X_Image(
    image.dir = file.path(sample_dir, "spatial"),
    image.name = "tissue_hires_image.png",
    assay = "Spatial"
  )

  object <- Seurat::Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = config$sample_id,
    image = image,
    filter.matrix = TRUE
  )
} else {
  object <- Seurat::Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = config$sample_id,
    filter.matrix = TRUE
  )
}

object$sample_id <- config$sample_id
object$original_sample_id <- config$original_sample_id
object$condition <- config$condition
object$orig.ident <- config$sample_id

if (any(grepl("^mt-", rownames(object)))) {
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object,
    pattern = "^mt-"
  )
} else if (any(grepl("^MT-", rownames(object)))) {
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object,
    pattern = "^MT-"
  )
} else {
  object$percent.mt <- 0
}

saveRDS(
  object,
  file.path(output_dir, "01_NOD_1_raw_spatial_object.rds")
)

message("Spots: ", ncol(object))
message("Genes: ", nrow(object))


# Quality control -------------------------------------------------------------

log_section("Quality control")

qc_summary <- function(object) {
  data.frame(
    sample_id = config$sample_id,
    original_sample_id = config$original_sample_id,
    condition = config$condition,
    n_spots = ncol(object),
    median_nCount_Spatial = stats::median(
      object$nCount_Spatial,
      na.rm = TRUE
    ),
    median_nFeature_Spatial = stats::median(
      object$nFeature_Spatial,
      na.rm = TRUE
    ),
    median_percent_mt = stats::median(
      object$percent.mt,
      na.rm = TRUE
    ),
    mean_nCount_Spatial = mean(
      object$nCount_Spatial,
      na.rm = TRUE
    ),
    mean_nFeature_Spatial = mean(
      object$nFeature_Spatial,
      na.rm = TRUE
    ),
    mean_percent_mt = mean(
      object$percent.mt,
      na.rm = TRUE
    )
  )
}

save_table(
  "QC_summary_before_filtering.csv",
  qc_summary(object)
)

qc_violin_before <- Seurat::VlnPlot(
  object,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 3
) +
  patchwork::plot_annotation(title = "NOD_1 QC before filtering") &
  theme_publication

qc_counts <- Seurat::SpatialFeaturePlot(
  object,
  features = "nCount_Spatial",
  pt.size.factor = config$spatial_point_size,
  image.scale = image_scale
) +
  ggplot2::ggtitle("Spatial distribution of total counts") &
  theme_publication

qc_features <- Seurat::SpatialFeaturePlot(
  object,
  features = "nFeature_Spatial",
  pt.size.factor = config$spatial_point_size,
  image.scale = image_scale
) +
  ggplot2::ggtitle("Spatial distribution of detected genes") &
  theme_publication

qc_mt <- Seurat::SpatialFeaturePlot(
  object,
  features = "percent.mt",
  pt.size.factor = config$spatial_point_size,
  image.scale = image_scale
) +
  ggplot2::ggtitle("Spatial distribution of mitochondrial percentage") &
  theme_publication

save_plot(
  "QC_violin_before_filtering.png",
  qc_violin_before,
  width = 12,
  height = 4
)
save_plot(
  "QC_spatial_counts_before_filtering.png",
  qc_counts,
  width = 7,
  height = 6
)
save_plot(
  "QC_spatial_features_before_filtering.png",
  qc_features,
  width = 7,
  height = 6
)
save_plot(
  "QC_spatial_percent_mt_before_filtering.png",
  qc_mt,
  width = 7,
  height = 6
)

object$pass_qc <-
  object$nFeature_Spatial >= config$min_features &
  object$nCount_Spatial >= config$min_counts &
  object$percent.mt <= config$max_percent_mt

save_table(
  "QC_filtering_summary.csv",
  data.frame(
    filter = c(
      "Total spots before filtering",
      "Spots passing filter",
      "Spots removed"
    ),
    n_spots = c(
      ncol(object),
      sum(object$pass_qc),
      sum(!object$pass_qc)
    )
  )
)

if (config$apply_qc_filter) {
  object <- subset(object, subset = pass_qc)
}

save_table(
  "QC_summary_after_filtering.csv",
  qc_summary(object)
)

qc_violin_after <- Seurat::VlnPlot(
  object,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 3
) +
  patchwork::plot_annotation(title = "NOD_1 QC after filtering") &
  theme_publication

save_plot(
  "QC_violin_after_filtering.png",
  qc_violin_after,
  width = 12,
  height = 4
)

saveRDS(
  object,
  file.path(output_dir, "02_NOD_1_after_QC.rds")
)


# Normalization and clustering ------------------------------------------------

log_section("Normalization and clustering")

Seurat::DefaultAssay(object) <- "Spatial"

object <- Seurat::SCTransform(
  object,
  assay = "Spatial",
  verbose = FALSE,
  return.only.var.genes = FALSE
)

saveRDS(
  object,
  file.path(output_dir, "03_NOD_1_SCT_normalized.rds")
)

Seurat::DefaultAssay(object) <- "SCT"

object <- Seurat::RunPCA(
  object,
  npcs = config$n_pcs,
  verbose = FALSE
)

pc_auto <- select_pc_range(
  object,
  max_pc = config$n_pcs
)
dims_use <- seq_len(
  min(
    max(config$dims),
    max(pc_auto),
    config$n_pcs
  )
)

pca_elbow <- Seurat::ElbowPlot(
  object,
  ndims = config$n_pcs
) +
  ggplot2::geom_vline(
    xintercept = max(dims_use),
    linetype = "dashed"
  ) +
  ggplot2::ggtitle("PCA elbow plot for NOD_1") +
  theme_publication

save_plot(
  "PCA_elbow_NOD_1.png",
  pca_elbow,
  width = 7,
  height = 5
)

object <- Seurat::RunUMAP(
  object,
  reduction = "pca",
  dims = dims_use,
  seed.use = config$seed,
  verbose = FALSE
)

object <- Seurat::FindNeighbors(
  object,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

object <- Seurat::FindClusters(
  object,
  resolution = config$cluster_resolution,
  random.seed = config$seed,
  verbose = FALSE
)

Seurat::Idents(object) <- "seurat_clusters"

cluster_counts <- as.data.frame(
  table(object$seurat_clusters),
  stringsAsFactors = FALSE
)
colnames(cluster_counts) <- c("cluster", "n_spots")
cluster_counts$percent_spots <- round(
  100 * cluster_counts$n_spots / sum(cluster_counts$n_spots),
  2
)

save_table(
  "cluster_spot_counts_NOD_1.csv",
  cluster_counts
)

umap_clusters <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggplot2::ggtitle("NOD_1 UMAP by Seurat cluster") +
  ggplot2::theme_classic(
    base_size = config$figure_base_size
  ) +
  theme_publication

save_plot(
  "UMAP_clusters_NOD_1.png",
  umap_clusters,
  width = 7,
  height = 6
)

spatial_clusters <- Seurat::SpatialDimPlot(
  object,
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size.factor = config$spatial_point_size,
  image.scale = image_scale
) +
  ggplot2::ggtitle("NOD_1 spatial clusters") &
  theme_publication

save_plot(
  "Spatial_clusters_NOD_1.png",
  spatial_clusters,
  width = 7,
  height = 7
)

saveRDS(
  object,
  file.path(output_dir, "04_NOD_1_clustered.rds")
)


# Cluster markers -------------------------------------------------------------

log_section("Cluster markers")

prepared <- tryCatch(
  Seurat::PrepSCTFindMarkers(
    object,
    assay = "SCT",
    verbose = TRUE
  ),
  error = function(e) {
    message(
      "PrepSCTFindMarkers failed: ",
      conditionMessage(e)
    )
    object
  }
)

object <- prepared
Seurat::Idents(object) <- "seurat_clusters"

cluster_markers <- Seurat::FindAllMarkers(
  object,
  assay = "SCT",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
) |>
  dplyr::arrange(
    cluster,
    p_val_adj,
    dplyr::desc(avg_log2FC)
  )

save_table(
  "cluster_markers_all_NOD_1.csv",
  cluster_markers
)

top10_markers <- cluster_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

save_table(
  "cluster_markers_top10_NOD_1.csv",
  top10_markers
)

heatmap_features <- intersect(
  unique(top10_markers$gene),
  rownames(object)
)

if (length(heatmap_features) > 0) {
  marker_heatmap <- tryCatch(
    Seurat::DoHeatmap(
      object,
      features = heatmap_features,
      group.by = "seurat_clusters"
    ) +
      Seurat::NoLegend() +
      ggplot2::ggtitle("Top cluster markers in NOD_1") +
      theme_publication,
    error = function(e) {
      message(
        "Marker heatmap was not generated: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(marker_heatmap)) {
    save_plot(
      "Cluster_marker_heatmap_NOD_1.png",
      marker_heatmap,
      width = 12,
      height = 14
    )
  }
}


# Biological programs ---------------------------------------------------------

log_section("Biological programs")

gene_panels <- list(
  epithelial_secretory = c(
    "Scgb2b17", "Ltf", "Aqp5", "Krt8",
    "Krt18", "Krt19", "Muc1", "Pip"
  ),
  antigen_presentation = c(
    "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "B2m"
  ),
  b_cell_plasma_cell = c(
    "Cd79a", "Ms4a1", "Jchain", "Mzb1", "Ighm", "Igha"
  ),
  t_cell = c(
    "Cd3d", "Cd3g", "Cd4", "Cd8a", "Trac"
  ),
  myeloid_stromal = c(
    "Gsn", "Lyz2", "Cd68", "C1qa",
    "C1qb", "Csf1r", "Col1a1", "Col3a1"
  ),
  interferon_chemokine = c(
    "Cxcl9", "Cxcl10", "Ccl5", "Ifit1",
    "Ifit3", "Isg15", "Stat1", "Irf1", "Irf7"
  ),
  metabolism_core = c(
    "Hk1", "Hk2", "Gpi1", "Pfkp", "Aldoa", "Gapdh",
    "Pgk1", "Eno1", "Pkm", "Ldha", "Pdha1", "Cs",
    "Aco2", "Idh3a", "Sdha", "Mdh2", "Ndufa1",
    "Ndufb8", "Uqcrc1", "Cox4i1", "Cox5a",
    "Atp5f1a", "Atp5f1b"
  ),
  lipid_metabolism = c(
    "Acaca", "Fasn", "Scd1", "Cpt1a", "Cpt2",
    "Acadm", "Ppara", "Pparg", "Fabp4", "Fabp5"
  ),
  oxidative_stress = c(
    "Sod1", "Sod2", "Gpx1", "Gpx3", "Cat",
    "Prdx1", "Prdx2", "Hmox1", "Nqo1"
  )
)

panel_presence <- dplyr::bind_rows(
  lapply(names(gene_panels), function(panel_name) {
    data.frame(
      panel = panel_name,
      gene = gene_panels[[panel_name]],
      present = gene_panels[[panel_name]] %in% rownames(object)
    )
  })
)

save_table(
  "gene_panel_presence_NOD_1.csv",
  panel_presence
)

for (panel_name in names(gene_panels)) {
  genes <- intersect(
    gene_panels[[panel_name]],
    rownames(object)
  )

  if (length(genes) < 2) {
    next
  }

  object <- Seurat::AddModuleScore(
    object,
    features = list(genes),
    name = paste0(panel_name, "_score"),
    assay = "SCT",
    ctrl = min(50, length(genes)),
    seed = config$seed
  )
}

module_score_cols <- grep(
  "_score1$",
  colnames(object@meta.data),
  value = TRUE
)

save_table(
  "module_score_columns_NOD_1.csv",
  data.frame(
    module_score_column = module_score_cols,
    biological_program = sub(
      "_score1$",
      "",
      module_score_cols
    )
  )
)

saveRDS(
  object,
  file.path(
    output_dir,
    "05_NOD_1_with_module_scores.rds"
  )
)


# Marker and module-score figures --------------------------------------------

log_section("Marker and module-score figures")

core_marker_genes <- c(
  "Scgb2b17", "Ltf", "Jchain", "Cd79a",
  "Cd3g", "Gsn", "Cd74", "H2-Aa",
  "H2-Ab1", "H2-Eb1", "Cxcl9", "Cxcl10",
  "Isg15", "Stat1", "Ldha", "Sod2"
)

core_marker_genes <- intersect(
  core_marker_genes,
  rownames(object)
)

if (length(core_marker_genes) > 0) {
  marker_umap <- Seurat::FeaturePlot(
    object,
    features = core_marker_genes,
    reduction = "umap",
    ncol = 4,
    order = TRUE,
    raster = FALSE
  ) &
    theme_publication

  save_plot(
    "UMAP_core_marker_genes_NOD_1.png",
    marker_umap,
    width = 14,
    height = 12
  )

  marker_dotplot <- Seurat::DotPlot(
    object,
    features = core_marker_genes,
    group.by = "seurat_clusters",
    assay = "SCT"
  ) +
    Seurat::RotatedAxis() +
    ggplot2::ggtitle("Core marker genes by NOD_1 cluster") +
    theme_publication

  save_plot(
    "DotPlot_core_marker_genes_by_cluster_NOD_1.png",
    marker_dotplot,
    width = 12,
    height = 5
  )

  marker_chunks <- split(
    core_marker_genes,
    ceiling(seq_along(core_marker_genes) / 6)
  )

  for (i in seq_along(marker_chunks)) {
    plot <- Seurat::SpatialFeaturePlot(
      object,
      features = marker_chunks[[i]],
      ncol = 2,
      pt.size.factor = config$spatial_point_size,
      image.scale = image_scale,
      alpha = c(0.1, 1)
    ) &
      theme_publication

    save_plot(
      sprintf(
        "Spatial_core_marker_genes_chunk_%d_NOD_1.png",
        i
      ),
      plot,
      width = 10,
      height = 12
    )
  }
}

if (length(module_score_cols) > 0) {
  module_dotplot <- Seurat::DotPlot(
    object,
    features = module_score_cols,
    group.by = "seurat_clusters"
  ) +
    Seurat::RotatedAxis() +
    ggplot2::scale_x_discrete(
      labels = sub(
        "_score1$",
        "",
        module_score_cols
      )
    ) +
    ggplot2::ggtitle("Biological module scores by cluster") +
    theme_publication

  save_plot(
    "DotPlot_module_scores_by_cluster_NOD_1.png",
    module_dotplot,
    width = 12,
    height = 5
  )

  module_chunks <- split(
    module_score_cols,
    ceiling(seq_along(module_score_cols) / 4)
  )

  for (i in seq_along(module_chunks)) {
    plots <- lapply(
      module_chunks[[i]],
      function(feature) {
        plot_spatial_feature(
          object,
          feature = feature,
          title = sub("_score1$", "", feature),
          assay = "SCT"
        )
      }
    )

    save_plot(
      sprintf(
        "Spatial_module_scores_chunk_%d_NOD_1.png",
        i
      ),
      patchwork::wrap_plots(plots, ncol = 2),
      width = 10,
      height = 10
    )
  }
}


# Spatially variable genes ----------------------------------------------------

log_section("Spatially variable genes")

Seurat::DefaultAssay(object) <- "SCT"

svg_candidates <- Seurat::VariableFeatures(object)

if (length(svg_candidates) < 100) {
  svg_candidates <- rownames(object)
}

svg_candidates <- intersect(
  svg_candidates,
  rownames(object)
)
svg_candidates <- head(
  svg_candidates,
  config$svg_n_features
)

svg_result <- run_spatial_variable_analysis(
  object,
  svg_candidates
)

if (!is.null(svg_result)) {
  object <- svg_result$object
  svg_genes <- svg_result$genes

  save_table(
    "spatially_variable_genes_NOD_1.csv",
    data.frame(
      rank = seq_along(svg_genes),
      gene = svg_genes,
      method = svg_result$method
    )
  )

  top_svg <- intersect(
    head(svg_genes, 12),
    rownames(object)
  )

  if (length(top_svg) > 0) {
    svg_chunks <- split(
      top_svg,
      ceiling(seq_along(top_svg) / 6)
    )

    for (i in seq_along(svg_chunks)) {
      plot <- Seurat::SpatialFeaturePlot(
        object,
        features = svg_chunks[[i]],
        ncol = 2,
        pt.size.factor = config$spatial_point_size,
        image.scale = image_scale,
        alpha = c(0.1, 1)
      ) &
        theme_publication

      save_plot(
        sprintf(
          "Spatial_top_SVG_chunk_%d_NOD_1.png",
          i
        ),
        plot,
        width = 10,
        height = 12
      )
    }
  }
}


# Spatial gradients -----------------------------------------------------------

log_section("Spatial gradients")

coordinates <- get_spatial_coordinates(object)

gradient_features <- unique(
  c(
    "Cd74", "Scgb2b17", "H2-Eb1", "H2-Aa",
    "Jchain", "Cd79a", "Ltf", "Gsn",
    "Cxcl9", "Cxcl10", "Isg15", "Stat1",
    "Ldha", "Sod2",
    module_score_cols
  )
)

gradient_features <- gradient_features[
  gradient_features %in% rownames(object) |
    gradient_features %in% colnames(object@meta.data)
]

gradient_results <- lapply(
  gradient_features,
  function(feature) {
    fit_spatial_gradient(
      object,
      feature,
      coordinates
    )
  }
)
names(gradient_results) <- gradient_features
gradient_results <- Filter(
  Negate(is.null),
  gradient_results
)

if (length(gradient_results) > 0) {
  gradient_summary <- dplyr::bind_rows(
    lapply(
      gradient_results,
      function(x) x$summary
    )
  ) |>
    dplyr::arrange(
      dplyr::desc(deviance_explained_percent)
    )

  save_table(
    "GAM_spatial_gradient_summary_NOD_1.csv",
    gradient_summary
  )

  top_gradient_features <- head(
    gradient_summary$feature,
    12
  )

  for (feature in top_gradient_features) {
    data <- gradient_results[[feature]]$data

    raw_plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = x, y = y, color = value)
    ) +
      ggplot2::geom_point(size = 1, alpha = 0.9) +
      viridis::scale_color_viridis(option = "C") +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_fixed() +
      ggplot2::theme_void(
        base_size = config$figure_base_size
      ) +
      ggplot2::labs(
        title = paste(feature, "raw expression"),
        color = "Expression"
      )

    fitted_plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = x, y = y, color = fitted)
    ) +
      ggplot2::geom_point(size = 1, alpha = 0.9) +
      viridis::scale_color_viridis(option = "C") +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_fixed() +
      ggplot2::theme_void(
        base_size = config$figure_base_size
      ) +
      ggplot2::labs(
        title = paste(feature, "fitted spatial gradient"),
        color = "GAM fitted"
      )

    save_plot(
      sprintf(
        "GAM_gradient_%s_NOD_1.png",
        clean_filename(feature)
      ),
      patchwork::wrap_plots(
        raw_plot,
        fitted_plot,
        ncol = 2
      ),
      width = 10,
      height = 5
    )
  }
}


# Spatial neighborhood analysis ----------------------------------------------

log_section("Spatial neighborhood analysis")

coordinates <- get_spatial_coordinates(object)
positions <- read_visium_positions(sample_dir)

positions <- positions[
  match(colnames(object), positions$barcode),
  ,
  drop = FALSE
]

if (anyNA(positions$barcode)) {
  stop(
    "Some retained spots could not be matched to ",
    "Space Ranger array coordinates."
  )
}

boundary_data <- data.frame(
  barcode = colnames(object),
  x = coordinates[colnames(object), "x"],
  y = coordinates[colnames(object), "y"],
  array_row = positions$array_row,
  array_col = positions$array_col,
  cluster = as.character(object$seurat_clusters),
  stringsAsFactors = FALSE
)

native_edges <- build_native_visium_edges(
  boundary_data
)

native_metrics <- compute_graph_metrics(
  boundary_data$cluster,
  native_edges,
  nrow(boundary_data)
)

boundary_data <- dplyr::bind_cols(
  boundary_data,
  native_metrics$spot_table
)

save_table(
  "boundary_spot_table_native_adjacency_NOD_1.csv",
  boundary_data
)

save_table(
  "boundary_summary_overall_native_adjacency_NOD_1.csv",
  native_metrics$overall
)

boundary_summary <- boundary_data |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n_spots = dplyr::n(),
    mean_native_degree = round(mean(degree), 3),
    mean_neighbor_diff_fraction = round(
      mean(neighbor_diff_fraction, na.rm = TRUE),
      4
    ),
    boundary_any_fraction = round(
      mean(is_boundary_any, na.rm = TRUE),
      4
    ),
    boundary_majority_fraction = round(
      mean(is_boundary_majority, na.rm = TRUE),
      4
    ),
    .groups = "drop"
  )

save_table(
  "boundary_summary_by_cluster_native_adjacency_NOD_1.csv",
  boundary_summary
)

boundary_map <- ggplot2::ggplot(
  boundary_data,
  ggplot2::aes(
    x = x,
    y = y,
    color = cluster,
    alpha = neighbor_diff_fraction
  )
) +
  ggplot2::geom_point(size = 1.1) +
  ggplot2::scale_alpha_continuous(
    range = c(0.15, 1),
    na.value = 0.15
  ) +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(
    base_size = config$figure_base_size
  ) +
  ggplot2::labs(
    title = "NOD_1 local cluster coherence",
    color = "Cluster",
    alpha = "Neighbor\ndifference"
  )

save_plot(
  "Boundary_map_native_adjacency_NOD_1.png",
  boundary_map,
  width = 8,
  height = 7
)

cluster_pairs <- data.frame(
  cluster_a = boundary_data$cluster[native_edges$i],
  cluster_b = boundary_data$cluster[native_edges$j],
  stringsAsFactors = FALSE
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    pair_a = sort(c(cluster_a, cluster_b))[1],
    pair_b = sort(c(cluster_a, cluster_b))[2]
  ) |>
  dplyr::ungroup() |>
  dplyr::count(
    pair_a,
    pair_b,
    name = "native_neighbor_pair_count"
  ) |>
  dplyr::arrange(
    dplyr::desc(native_neighbor_pair_count)
  )

save_table(
  "cluster_neighbor_adjacency_native_NOD_1.csv",
  cluster_pairs
)

native_edge_distance <- sqrt(
  (
    boundary_data$x[native_edges$i] -
      boundary_data$x[native_edges$j]
  )^2 +
    (
      boundary_data$y[native_edges$i] -
        boundary_data$y[native_edges$j]
    )^2
)

native_spacing <- stats::median(
  native_edge_distance,
  na.rm = TRUE
)
radius_cutoff <- config$radius_multiplier * native_spacing
k_radius <- min(12L, nrow(boundary_data) - 1L)

knn <- FNN::get.knn(
  as.matrix(boundary_data[, c("x", "y")]),
  k = k_radius
)

radius_edge_list <- lapply(
  seq_len(nrow(boundary_data)),
  function(i) {
    keep <- which(
      knn$nn.dist[i, ] <= radius_cutoff
    )

    if (length(keep) == 0) {
      return(NULL)
    }

    data.frame(
      i = i,
      j = knn$nn.index[i, keep]
    )
  }
)

radius_edges <- dplyr::bind_rows(
  radius_edge_list
) |>
  dplyr::mutate(
    i_min = pmin(i, j),
    j_max = pmax(i, j)
  ) |>
  dplyr::transmute(
    i = i_min,
    j = j_max
  ) |>
  dplyr::filter(i != j) |>
  dplyr::distinct()

radius_metrics <- compute_graph_metrics(
  boundary_data$cluster,
  radius_edges,
  nrow(boundary_data)
)

graph_sensitivity <- dplyr::bind_rows(
  native_metrics$overall |>
    dplyr::mutate(graph = "native_array_adjacency"),
  radius_metrics$overall |>
    dplyr::mutate(graph = "distance_limited")
) |>
  dplyr::mutate(
    native_spacing_pixels = native_spacing,
    radius_multiplier = dplyr::if_else(
      graph == "distance_limited",
      config$radius_multiplier,
      NA_real_
    ),
    radius_cutoff_pixels = dplyr::if_else(
      graph == "distance_limited",
      radius_cutoff,
      NA_real_
    )
  ) |>
  dplyr::select(graph, dplyr::everything())

native_edge_key <- paste(
  native_edges$i,
  native_edges$j,
  sep = "_"
)
radius_edge_key <- paste(
  radius_edges$i,
  radius_edges$j,
  sep = "_"
)

edge_union <- union(
  native_edge_key,
  radius_edge_key
)
edge_jaccard <- if (length(edge_union) == 0) {
  NA_real_
} else {
  length(
    intersect(
      native_edge_key,
      radius_edge_key
    )
  ) / length(edge_union)
}

graph_agreement <- data.frame(
  spearman_neighbor_diff = suppressWarnings(
    stats::cor(
      native_metrics$spot_table$neighbor_diff_fraction,
      radius_metrics$spot_table$neighbor_diff_fraction,
      use = "complete.obs",
      method = "spearman"
    )
  ),
  boundary_any_agreement = mean(
    native_metrics$spot_table$is_boundary_any ==
      radius_metrics$spot_table$is_boundary_any,
    na.rm = TRUE
  ),
  edge_jaccard = edge_jaccard
)

save_table(
  "Figure8_graph_sensitivity_NOD_1.csv",
  graph_sensitivity
)

save_table(
  "Figure8_graph_agreement_NOD_1.csv",
  graph_agreement
)

resolution_values <- sort(
  unique(
    c(
      config$resolution_sensitivity,
      config$cluster_resolution
    )
  )
)

resolution_sensitivity <- dplyr::bind_rows(
  lapply(
    resolution_values,
    function(resolution) {
      temp_object <- Seurat::FindClusters(
        object,
        resolution = resolution,
        random.seed = config$seed,
        verbose = FALSE
      )

      labels <- as.character(
        temp_object$seurat_clusters
      )

      metrics <- compute_graph_metrics(
        labels,
        native_edges,
        nrow(boundary_data)
      )

      data.frame(
        resolution = resolution,
        n_clusters = length(unique(labels)),
        mean_neighbor_diff_fraction =
          metrics$overall$mean_neighbor_diff_fraction,
        boundary_any_fraction =
          metrics$overall$boundary_any_fraction,
        boundary_majority_fraction =
          metrics$overall$boundary_majority_fraction
      )
    }
  )
)

save_table(
  "Figure8_cluster_resolution_sensitivity_NOD_1.csv",
  resolution_sensitivity
)

resolution_mixing_plot <- ggplot2::ggplot(
  resolution_sensitivity,
  ggplot2::aes(
    x = resolution,
    y = mean_neighbor_diff_fraction
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::theme_classic(
    base_size = config$figure_base_size
  ) +
  ggplot2::labs(
    title = "Sensitivity to clustering resolution",
    x = "Seurat clustering resolution",
    y = "Mean neighbor-difference fraction"
  )

resolution_boundary_plot <- ggplot2::ggplot(
  resolution_sensitivity,
  ggplot2::aes(
    x = resolution,
    y = boundary_any_fraction
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::theme_classic(
    base_size = config$figure_base_size
  ) +
  ggplot2::labs(
    title = "Binary interface sensitivity",
    x = "Seurat clustering resolution",
    y = "Fraction with at least one different neighbor"
  )

save_plot(
  "Figure8_cluster_resolution_sensitivity_NOD_1.png",
  patchwork::wrap_plots(
    resolution_mixing_plot,
    resolution_boundary_plot,
    ncol = 2
  ),
  width = 11,
  height = 5
)


# Gene Ontology enrichment ----------------------------------------------------

if (config$run_go) {
  log_section("Gene Ontology enrichment")

  significant_markers <- cluster_markers |>
    dplyr::filter(
      p_val_adj < 0.05,
      avg_log2FC > 0.25
    )

  cluster_ids <- sort(
    unique(
      as.character(significant_markers$cluster)
    )
  )

  for (cluster_id in cluster_ids) {
    genes <- significant_markers |>
      dplyr::filter(
        as.character(cluster) == cluster_id
      ) |>
      dplyr::pull(gene)

    run_go_enrichment(
      genes,
      output_prefix = paste0(
        "cluster_",
        clean_filename(cluster_id),
        "_marker_genes"
      ),
      title = paste(
        "GO BP enrichment: cluster",
        cluster_id,
        "marker genes"
      )
    )
  }
}


# Final outputs ---------------------------------------------------------------

log_section("Final outputs")

saveRDS(
  object,
  file.path(
    output_dir,
    "final_NOD_1_single_sample_spatial_context.rds"
  )
)

analysis_parameters <- data.frame(
  parameter = c(
    "sample_folder",
    "sample_id",
    "original_sample_id",
    "condition",
    "seed",
    "apply_qc_filter",
    "min_features",
    "min_counts",
    "max_percent_mt",
    "n_pcs",
    "dims_used",
    "cluster_resolution",
    "svg_n_features",
    "svg_methods",
    "figure_dpi",
    "resolution_sensitivity",
    "radius_multiplier"
  ),
  value = c(
    config$sample_folder,
    config$sample_id,
    config$original_sample_id,
    config$condition,
    config$seed,
    config$apply_qc_filter,
    config$min_features,
    config$min_counts,
    config$max_percent_mt,
    config$n_pcs,
    paste(range(dims_use), collapse = "-"),
    config$cluster_resolution,
    config$svg_n_features,
    paste(config$svg_methods, collapse = ","),
    config$figure_dpi,
    paste(resolution_values, collapse = ","),
    config$radius_multiplier
  )
)

save_table(
  "analysis_parameters.csv",
  analysis_parameters
)

output_manifest <- data.frame(
  output_group = c(
    "RDS",
    "QC",
    "Clustering",
    "Marker genes",
    "Biological programs",
    "Spatially variable genes",
    "Spatial gradients",
    "Spatial neighborhood",
    "GO enrichment",
    "Reproducibility"
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
    "cluster_*_marker_genes_GO_BP.csv; cluster_*_marker_genes_GO_BP_dotplot.png",
    "analysis_parameters.csv; sessionInfo.txt"
  ),
  stringsAsFactors = FALSE
)

save_table(
  "output_manifest_NOD_1_single_sample_pipeline.csv",
  output_manifest
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("Analysis complete: ", output_dir)
