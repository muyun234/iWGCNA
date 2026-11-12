# iWGCNA

An interactive Shiny application for weighted gene co-expression network
analysis (WGCNA), visualization, and functional enrichment analysis.

iWGCNA provides a graphical interface for exploring gene co-expression
modules, investigating module–trait relationships, identifying candidate
hub genes, and interpreting biological functions without writing R code.

## Web Application

Access the online application: [Launch iWGCNA](https://shiny.cgu.edu.tw/app/iWGCNA)

## Documentation

For data preparation, parameter settings, and analysis workflows,
please refer to the [iWGCNA User Guide](https://appliedbioinformatics.notion.site/ebd//2cd160e47a6980e3a0abe769e96eed7b).

## Main Features

- **Data input:** Upload gene expression matrices and sample trait data.
- **Data preprocessing:** Filter genes, transform expression values,
  select variable genes, and inspect sample clustering.
- **Network construction:** Evaluate soft-thresholding powers and
  identify co-expression modules.
- **Module–trait analysis:** Explore correlations between module
  eigengenes and sample traits.
- **Hub gene and network analysis:** Examine module membership and
  gene significance, identify candidate hub genes, and visualize
  co-expression and STRING protein–protein interaction networks.
- **Functional enrichment analysis:** Perform Gene Ontology (GO)
  and KEGG pathway enrichment analyses.
- **Visualization and export:** Generate and download analysis
  results and figures.

## Input Data

- **Gene expression matrix:** A `.txt`, `.tsv`, or `.csv` file with
  genes in rows and samples in columns. The first column should
  contain gene identifiers.
- **Sample trait data:** A table containing sample information for
  module–trait analysis. Sample identifiers must match those in
  the expression matrix.

Please refer to the user guide for supported formats and example datasets.

## Requirements

The online application can be accessed without installing R locally.

For local use, iWGCNA was developed using **R version 4.5.1**.
The package versions used in the development environment are listed below.

| Package | Version |
| --- | --- |
| shiny | 1.14.0 |
| bs4Dash | 2.3.4 |
| shinyWidgets | 0.9.0 |
| shinyjs | 2.1.0 |
| shinyalert | 3.1.0 |
| waiter | 0.2.5 |
| summaryBox | 0.1.0 |
| DT | 0.33 |
| reactable | 0.4.5 |
| htmltools | 0.5.9 |
| shinyscreenshot | 0.2.1 |
| dplyr | 1.2.1 |
| tibble | 3.3.1 |
| magrittr | 2.0.5 |
| readr | 2.1.5 |
| purrr | 1.2.2 |
| fastDummies | 1.7.5 |
| GEOquery | 2.76.0 |
| Biobase | 2.68.0 |
| limma | 3.64.3 |
| DESeq2 | 1.48.1 |
| WGCNA | 1.73 |
| AnnotationDbi | 1.70.0 |
| ClusterGVis | 0.1.4 |
| clusterProfiler | 4.16.0 |
| ggplot2 | 4.0.3 |
| ggdendro | 0.2.0 |
| ggprism | 1.0.7 |
| ggrepel | 0.9.6 |
| patchwork | 1.3.1 |
| plotly | 4.12.1 |
| ComplexHeatmap | 2.24.1 |
| circlize | 0.4.16 |
| scales | 1.4.0 |
| igraph | 2.1.4 |
| tidygraph | 1.3.1 |
| ggraph | 2.2.1 |

## Feedback

For bug reports and feature requests, please open an issue in this
repository. Include a description of the problem, steps to reproduce
it, and relevant error messages.

## License

This project is licensed under the MIT License.
See the [LICENSE](LICENSE) file for details.
