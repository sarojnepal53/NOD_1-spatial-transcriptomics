[README(6).md](https://github.com/user-attachments/files/32173587/README.6.md)
# NOD_1 Spatial Transcriptomics Analysis

This repository contains a reproducible **single-sample 10x Genomics Visium spatial transcriptomics workflow** for the lacrimal gland sample **NOD.H-2b_1**, analyzed as `Case_1`.

The analysis is implemented in **R using Seurat** and focuses on spatial organization within one disease-model tissue section. It is **not** a case-versus-control differential-expression analysis.

## Repository contents

```text
NOD_1-spatial-transcriptomics/
├── README.md
├── R_code.R
├── Supplementary_File.docx
├── Data and Result.rar
├── .gitattributes
└── .gitignore
```

### `Data and Result.rar`

The archive contains both the original/processed Visium input for `NOD_1` and the complete analysis output folder.

```text
Data and Result.rar
├── NOD_1/
│   ├── filtered_feature_bc_matrix.h5
│   ├── filtered_feature_bc_matrix/
│   ├── spatial/
│   │   ├── scalefactors_json.json
│   │   ├── tissue_hires_image.png
│   │   ├── tissue_lowres_image.png
│   │   └── tissue_positions*.csv
│   └── SLIDE2C1.tif
│
└── NOD_1_single_sample_spatial_context/
    ├── processed Seurat objects (.rds)
    ├── QC tables and plots
    ├── clustering results
    ├── marker-gene tables and figures
    ├── biological module-score results
    ├── spatially variable gene results
    ├── GAM spatial-gradient results
    ├── native Visium boundary/neighborhood results
    ├── Figure 8 sensitivity analyses
    └── GO enrichment outputs
```

Because the archive is large, `Data and Result.rar` is stored using **Git LFS**.

## Analysis workflow

`R_code.R` performs the following steps:

1. Load the `NOD_1` Visium sample
2. Quality control and optional spot filtering
3. SCTransform normalization
4. PCA, UMAP, and graph-based clustering
5. Spatial visualization of clusters
6. Cluster marker-gene analysis
7. Selected marker-gene visualization
8. Biological program scoring
9. Spatially variable gene analysis
10. Spatial-gradient modeling using generalized additive models (GAMs)
11. Native Visium boundary and neighborhood analysis
12. Graph and clustering-resolution sensitivity analyses
13. Optional Gene Ontology enrichment
14. Save the final processed object and output manifest

## Main outputs

### Processed Seurat objects

```text
01_NOD_1_raw_spatial_object.rds
02_NOD_1_after_QC.rds
03_NOD_1_SCT_normalized.rds
04_NOD_1_clustered.rds
05_NOD_1_with_module_scores.rds
final_NOD_1_single_sample_spatial_context.rds
```

### Main result tables

```text
QC_summary_before_filtering.csv
QC_summary_after_filtering.csv
cluster_spot_counts_NOD_1.csv
cluster_markers_all_NOD_1.csv
cluster_markers_top10_NOD_1.csv
gene_panel_presence_NOD_1.csv
spatially_variable_genes_NOD_1.csv
GAM_spatial_gradient_summary_NOD_1.csv
boundary_spot_table_native_adjacency_NOD_1.csv
boundary_summary_overall_native_adjacency_NOD_1.csv
boundary_summary_by_cluster_native_adjacency_NOD_1.csv
cluster_neighbor_adjacency_native_NOD_1.csv
Figure8_graph_sensitivity_NOD_1.csv
Figure8_graph_agreement_NOD_1.csv
Figure8_cluster_resolution_sensitivity_NOD_1.csv
output_manifest_NOD_1_single_sample_pipeline.csv
```

### Main figures

```text
QC_violin_before_filtering.png
QC_violin_after_filtering.png
QC_spatial_counts_before_filtering.png
QC_spatial_features_before_filtering.png
QC_spatial_percent_mt_before_filtering.png

PCA_elbow_NOD_1.png
UMAP_clusters_NOD_1.png
Spatial_clusters_NOD_1.png

Cluster_marker_heatmap_NOD_1.png
DotPlot_core_marker_genes_by_cluster_NOD_1.png
Spatial_core_marker_genes_chunk_*_NOD_1.png

DotPlot_module_scores_by_cluster_NOD_1.png
Spatial_module_scores_chunk_*_NOD_1.png

Spatial_top_SVG_chunk_*_NOD_1.png
GAM_gradient_*_NOD_1.png

Boundary_map_native_adjacency_NOD_1.png
Figure8_cluster_resolution_sensitivity_NOD_1.png
```

## Boundary and Figure 8 analysis

The revised boundary analysis uses the **native Visium hexagonal array geometry** rather than unrestricted nearest-neighbor connections.

Immediate neighbors are defined from the Visium array coordinates, so missing, non-tissue, or QC-removed spots are not replaced by more distant spots.

The primary local-mixing measure is the fraction of native neighboring spots assigned to a different cluster.

Two sensitivity analyses are included:

- **Graph sensitivity:** native Visium adjacency is compared with a distance-limited graph.
- **Clustering-resolution sensitivity:** the analysis is repeated at Seurat resolutions `0.30`, `0.45`, and `0.60`.

Relevant files include:

```text
Figure8_graph_sensitivity_NOD_1.csv
Figure8_graph_agreement_NOD_1.csv
Figure8_cluster_resolution_sensitivity_NOD_1.csv
Figure8_cluster_resolution_sensitivity_NOD_1.png
```

## Biological signals examined

The workflow evaluates selected marker genes and gene programs associated with:

- epithelial and secretory signals
- antigen presentation
- B-cell and plasma-cell signals
- T-cell signals
- myeloid and stromal signals
- interferon and chemokine signaling
- core metabolism
- lipid metabolism
- oxidative stress

Because Visium spots can contain transcripts from multiple cells, these results should be interpreted as **spatial expression signals**, not definitive single-cell identities.

## How to run

Clone the repository:

```bash
git clone https://github.com/sarojnepal53/NOD_1-spatial-transcriptomics.git
cd NOD_1-spatial-transcriptomics
```

Download the Git LFS files:

```bash
git lfs pull
```

Extract:

```text
Data and Result.rar
```

Then open `R_code.R` in RStudio and update:

```r
root_dir <- "C:/Users/yourname/path/to/spatial transcriptomics"
```

`root_dir` should point to the directory containing the extracted `NOD_1` folder.

Run the script from top to bottom. New results will be written to:

```text
NOD_1_single_sample_spatial_context/
```

## Main R packages

The workflow uses:

```text
Seurat
SeuratObject
ggplot2
patchwork
dplyr
tidyr
stringr
ggrepel
mgcv
Matrix
FNN
viridis
scales
future
```

Optional Gene Ontology enrichment uses:

```text
clusterProfiler
org.Mm.eg.db
enrichplot
```

## Data source

This secondary analysis uses publicly available data from:

**Mauduit, O., Delcroix, V., Umazume, T., de Paiva, C. S., Dartt, D. A., & Makarenkova, H. P. (2022).**  
*Spatial transcriptomics of the lacrimal gland features macrophage activity and epithelium metabolism as key alterations during chronic inflammation.*  
**Frontiers in Immunology, 13**, 1011125.  
https://doi.org/10.3389/fimmu.2022.1011125

GEO accessions:

- Bulk RNA-seq: `GSE210332`
- Visium spatial gene-expression dataset: `GSE210380`

## Notes

- This repository analyzes one sample only.
- Cluster-marker analysis compares clusters within `NOD_1`.
- Cluster identities and spatial patterns can change with QC thresholds, PCA dimensions, clustering resolution, or software versions.
- Marker genes and module scores support biological interpretation but are not definitive cell-type annotations.
- Large data and result files are distributed through Git LFS in `Data and Result.rar`.
