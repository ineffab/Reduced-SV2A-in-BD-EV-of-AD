# Reduced SV2A in brain-derived extracellular vesicles in Alzheimer's disease

Proteomic network analysis code for:

> Nussbaumer J., **Barve A.**, Zufferey V., Espourteille J., Kirabali T., Konietzko U., Razansky D., Nordberg A., Buée L., Colin M., Nitsch R.M., Hock C., Richetin K., Ni R. **Reduced synaptic vesicle protein 2A in extracellular vesicles and brains of Alzheimer's disease: associations with Aβ, tau, synaptic proteins and APOE ε4.** *Translational Neurodegeneration* (2025). [doi:10.1186/s40035-025-00508-2](https://doi.org/10.1186/s40035-025-00508-2) · [preprint](https://www.biorxiv.org/content/10.1101/2025.01.21.634009)

This repository holds the brain-derived extracellular vesicle (BD-EV) proteomics analysis I contributed to the paper: co-expression networks, functional and cell-type annotation, synaptic protein correlations and protein–protein interaction networks.

## Summary

Synaptic vesicle protein 2A (SV2A) is the target of synaptic-density PET imaging in Alzheimer's disease (AD). We asked whether SV2A and other synaptic proteins change in vesicles released by brain cells.

BD-EVs were isolated from post-mortem prefrontal cortex of AD cases across Braak stages and non-demented controls, and profiled by label-free LC-MS/MS. Co-expression analysis identified eight protein modules. A **synaptic-vesicle module (M3)**, with SV2A as a central hub, was reduced in AD and was the only module altered from early Braak stages onward. GAP43 and SYT1 were reduced in BD-EVs from early stages, and SV2A at late stages.

| | |
|---|---|
| **Samples** | Post-mortem prefrontal cortex: 4 controls, 17 AD (Braak 0–II, III–IV, V–VI) |
| **Proteomics** | Label-free LC-MS/MS of BD-EVs; riBAQ protein intensities |
| **Analyses** | WGCNA, GO enrichment, cell-type enrichment, Spearman correlations with a synaptic protein panel, SIGNOR protein–protein interaction networks |

## Repository structure

```
├── WGCNA/          Co-expression modules, differential module analysis, functional and cell-type annotation
├── Correlations/   SV2A correlations with a 24-protein CSF synaptic panel (GraphPad Prism files)
└── Cytoscape/      Protein–protein interaction networks (SV2A network, module hub proteins)
```

## Analysis workflow

The R scripts in `WGCNA/` run as a chain. Each one calls the previous one with `source()`, so run them in this order:

| Step | Script | What it does |
|---|---|---|
| 1 | `generate_modules_WGCNA_AD_class.R` | Sample outlier check, soft-threshold power 16, topological overlap, dynamic tree cut (minimum module size 30, deepSplit 2), module merging (cut height 0.15), module eigengenes and kME, limma tests of module eigengenes and proteins between groups, hub-protein selection |
| 2 | `Functional_annotation_WGCNA_AD_class.R` | GO, KEGG, Reactome and WikiPathways enrichment of each module with g:Profiler (FDR correction) |
| 3 | `ORA_pruning_WGCNA_AD_class.R` | Reduces redundant GO terms by semantic similarity (GOSemSim, rrvgo) |
| 4 | `Cell_type_FET_WGCNA_AD_class.R` | Cell-type enrichment of modules against brain cell-type marker lists with Fisher's exact test, using `geneListFET.R` |

**Correlations:** Spearman correlations between SV2A and the synaptic protein panel (21 of 24 proteins detected in BD-EVs), in `Correlations/`.

**Networks:** protein–protein interactions from SIGNOR, visualised in Cytoscape 3.10 (`Cytoscape/Sv2a_PPI.cys`).

## How to run

Requires R ≥ 4.2 with `WGCNA`, `flashClust`, `limma`, `gprofiler2`, `GOSemSim`, `rrvgo`, `org.Hs.eg.db`, `clusterProfiler`, `biomaRt`, `HGNChelper`, `readxl` and `tidyverse` (CRAN and Bioconductor).

Before running:
1. Replace the `setwd(...)` and `source("C:/...")` paths at the top of each script with your local paths.
2. Place the protein intensity matrix and sample metadata in `WGCNA/`.

## Acknowledgements

`geneListFET.R` and the cell-type Fisher's exact test wrapper are by Eric Dammer and Divya Nandakumar (Emory University), as published in Seyfried et al., *Cell Systems* 2017 and Johnson et al., *Nature Medicine* 2020, and are used here with attribution.

## Citation

If you use this code, please cite the paper above.

## Contact

Aatmika Barve · [ORCID 0000-0003-2331-5825](https://orcid.org/0000-0003-2331-5825) · [LinkedIn](https://www.linkedin.com/in/aatmika-b-504aaa212)
