[README(7).md](https://github.com/user-attachments/files/32173701/README.7.md)
# NOD_1 Spatial Transcriptomics Analysis

A reproducible single-sample **10x Genomics Visium** spatial transcriptomics workflow for the lacrimal gland sample **NOD.H-2b_1**, analyzed as `Case_1`.

The analysis is implemented in **R with Seurat** and is designed to characterize transcriptional and spatial organization within one disease-model tissue section. It includes quality control, normalization, clustering, marker analysis, biological program scoring, spatially variable gene analysis, spatial-gradient modeling, neighborhood analysis, sensitivity analyses, and optional Gene Ontology enrichment.

> **Scope:** This is a single-sample spatial-context analysis. It should not be interpreted as a case-versus-control differential-expression analysis.

## Repository structure

```text
NOD_1-spatial-transcriptomics/
├── README.md
├── R_code.R
├── Supplementary_File.docx
├── Data and Result.rar
├── .gitattributes
└── .gitignore
```

`Data and Result.rar` is stored with **Git LFS** because it contains large input and output files.

The archive contains:

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
    ├── processed Seurat objects
    ├── QC tables and figures
    ├── clustering outputs
    ├── marker-gene results
    ├── biological program scores
    ├── spatially variable gene results
    ├── GAM spatial-gradient results
    ├── spatial neighborhood and boundary results
    ├── sensitivity-analysis outputs
    └── Gene Ontology enrichment results
```

## Analysis workflow

The main script, `R_code.R`, performs the following steps:

1. Load the `NOD_1` Visium dataset.
2. Calculate quality-control metrics and optionally filter low-quality spots.
3. Normalize expression using `SCTransform()`.
4. Perform PCA, UMAP, nearest-neighbor analysis, and Seurat graph-based clustering.
5. Map transcriptional clusters back to tissue space.
6. Identify positive marker genes for each cluster.
7. Visualize selected marker genes in UMAP and spatial coordinates.
8. Calculate biological program scores.
9. Identify spatially variable genes.
10. Model spatial expression gradients using generalized additive models (GAMs).
11. Quantify local cluster organization using Visium array neighbors.
12. Evaluate sensitivity to graph definition and clustering resolution.
13. Perform optional Gene Ontology Biological Process enrichment.
14. Save processed objects, summary tables, figures, and an output manifest.

## Key analysis settings

| Parameter | Setting |
|---|---|
| Sample folder | `NOD_1` |
| Analysis label | `Case_1` |
| Original sample ID | `NOD.H-2b_1` |
| QC: minimum detected features | 200 |
| QC: minimum counts | 500 |
| QC: maximum mitochondrial percentage | 5% |
| Normalization | SCTransform |
| PCA components calculated | 50 |
| Clustering resolution | 0.45 |
| Figure resolution | 300 dpi |
| GO enrichment | Enabled by default |

The clustering-resolution sensitivity analysis additionally evaluates resolutions `0.30` and `0.60`.

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

### Core result tables

```text
QC_summary_before_filtering.csv
QC_summary_after_filtering.csv
QC_filtering_summary.csv

cluster_spot_counts_NOD_1.csv
cluster_markers_all_NOD_1.csv
cluster_markers_top10_NOD_1.csv

gene_panel_presence_NOD_1.csv
module_score_columns_NOD_1.csv
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

### Representative figures

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

## Spatial neighborhood analysis

Local spatial organization is evaluated using the **native Visium array geometry**. Immediate neighbors are defined from the array row and column coordinates rather than by unrestricted nearest-neighbor assignment.

For each retained spot, the analysis calculates the fraction of immediate Visium neighbors assigned to a different cluster. This provides a local measure of cluster mixing and spatial coherence.

Two sensitivity analyses are included:

- **Graph-definition sensitivity:** native Visium adjacency is compared with a distance-limited coordinate graph.
- **Clustering-resolution sensitivity:** neighborhood statistics are recalculated at Seurat resolutions `0.30`, `0.45`, and `0.60`.

The main outputs are:

```text
Boundary_map_native_adjacency_NOD_1.png
Figure8_graph_sensitivity_NOD_1.csv
Figure8_graph_agreement_NOD_1.csv
Figure8_cluster_resolution_sensitivity_NOD_1.csv
Figure8_cluster_resolution_sensitivity_NOD_1.png
```

## Biological programs examined

The workflow evaluates selected genes and gene sets representing:

- epithelial and secretory programs
- antigen presentation
- B-cell and plasma-cell signals
- T-cell signals
- myeloid and stromal signals
- interferon and chemokine signaling
- core metabolism
- lipid metabolism
- oxidative stress

Because a Visium spot can contain transcripts from multiple cells, these measurements should be interpreted as **spatial expression signals**, not definitive single-cell identities.

## Running the analysis

### 1. Clone the repository

```bash
git clone https://github.com/sarojnepal53/NOD_1-spatial-transcriptomics.git
cd NOD_1-spatial-transcriptomics
```

### 2. Retrieve Git LFS files

Install Git LFS if necessary, then run:

```bash
git lfs install
git lfs pull
```

### 3. Extract the data archive

Extract:

```text
Data and Result.rar
```

The extracted directory should contain both `NOD_1/` and `NOD_1_single_sample_spatial_context/`.

### 4. Set the project path

Open `R_code.R` in RStudio and edit:

```r
root_dir <- "C:/path/to/spatial transcriptomics"
```

`root_dir` must point to the directory containing the extracted `NOD_1` folder.

### 5. Run the script

Run `R_code.R` from top to bottom.

New outputs are written automatically to:

```text
NOD_1_single_sample_spatial_context/
```

## R dependencies

Core packages:

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

The script checks for missing packages and installs them when needed.

## Interpretation and limitations

This repository contains a **single tissue-section analysis**. Therefore:

- cluster-marker tests compare transcriptional clusters within `NOD_1`;
- no disease-versus-control statistical inference is performed;
- cluster identities can vary with QC thresholds, PCA dimensions, clustering resolution, and software versions;
- marker genes and module scores support biological interpretation but are not definitive cell-type annotations;
- spatial neighborhood statistics describe organization of the analyzed section and should not be generalized to biological replicates without additional samples.

## Data source

The analysis uses publicly available lacrimal gland spatial transcriptomics data from:

**Mauduit, O., Delcroix, V., Umazume, T., de Paiva, C. S., Dartt, D. A., & Makarenkova, H. P. (2022).**  
*Spatial transcriptomics of the lacrimal gland features macrophage activity and epithelium metabolism as key alterations during chronic inflammation.*  
**Frontiers in Immunology, 13**, 1011125.  
https://doi.org/10.3389/fimmu.2022.1011125

GEO accessions:

- Bulk RNA-seq: `GSE210332`
- Visium spatial gene expression: `GSE210380`

## Citation

If you use the underlying dataset, cite the original study above.

If you reuse or adapt this analysis workflow, please also cite or link to this repository so that the computational workflow can be traced and reproduced.
