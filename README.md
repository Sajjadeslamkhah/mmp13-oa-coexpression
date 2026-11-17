# MMP13 co-expression in osteoarthritic joint tissues

Code, results and supplementary material for the manuscript

> Aydoğdu, E.D.; Eslamkhah, S.; Gülmez, H.; Aslan, E.S.; Şahin, H. *Beyond the case–control signal: MMP13 is embedded in a reproducible hypertrophic and fibro-stromal co-expression programme within osteoarthritic joint tissues.* Under review at *Current Issues in Molecular Biology* (MDPI).

This repository is the data source cited in the manuscript's Data Availability Statement; the version referred to there is release `v1.0.0`.

## What the study does

- Re-analyses seven public GEO microarray datasets from four joint tissues (cartilage, meniscus, synovium, subchondral bone; 123 OA and 66 control samples).
- Estimates the OA-versus-control difference in MMP13 with a random-effects meta-analysis and tests its sensitivity to donor age and sex.
- Builds a consensus MMP13 co-expression neighbourhood within OA tissue, locks the top partners, and tests them in held-out cartilage (GSE114007) and synovium (GSE82107) and in single-cell data (GSE104782).
- Revision analyses: comparison with other seed genes and expression-matched random genes, network centrality, cell-type composition adjustment, co-expression in control tissue, correlation-aware gene-set tests (CAMERA, sample-label permutation) and single-cell sensitivity analysis.

No new samples were generated. All input data are public in GEO. The datasets were taken from a list compiled by the authors before analysis; no systematic database search was performed (see `data/assessed_series.csv` and Table S5).

## Repository layout

| Path | Content |
|---|---|
| `notebooks/01_discovery_analysis.ipynb` | Discovery analysis: data loading, case–control meta-analysis, co-expression meta-analysis, paired analysis, programme coupling, GSEA, enrichment, drug–gene records, age sensitivity |
| `notebooks/02_validation_analysis.ipynb` | Locked signature; held-out bulk validation; single-cell localisation |
| `notebooks/03_revision_analyses.ipynb` | Analyses added during peer review |
| `src/revision_analyses.py` | Functions used by the revision notebook |
| `scripts/make_figures.py` | Draws Figure 8 (with panel e) and Figures S5–S7 from `results/revision` |
| `scripts/make_supplementary_tables.py` | Builds `supplementary/Supplementary_Tables_S5-S12.xlsx` |
| `data/assessed_series.csv` | Every GEO series assessed for inclusion, with the decision and reason |
| `results/discovery/` | Discovery tables (`tables/`: S0–S4 dataset audit, duplicate-array check, donor age, age-adjusted models and age-controlled co-expression; T1–T9 effect sizes, meta-analysis, AUCs, consensus co-expression for all genes, programme coupling, GSEA, enrichment, DGIdb records, paired analysis) and Figures 1–7 and S1–S3 with source data (`figures/`) |
| `results/validation/` | Locked signature (`V0_locked_signature.json`), validation tables and figures with source data |
| `results/revision/` | Gene-level and summary outputs of the revision analyses (see `results/README.md`) |
| `supplementary/` | Supplementary Figures S5–S7 (PNG and Word file with legends), Supplementary Tables S5–S12 (Excel), Figure 8 with panel e |

## Running the notebooks

The notebooks were run on Google Colab with Google Drive mounted. They expect a working folder `MyDrive/mmp13_oa/`; change `BASE` in the first configuration cell to use another location. Run them in order (01, 02, 03): the validation notebook reads the discovery table, and the revision notebook reads both. GEO series are downloaded with GEOparse on first use.

Package versions used for the published results are listed in `requirements.txt`. The reproduction gate in `03_revision_analyses.ipynb` compares the rerun with the published 1,510 / 407 robust partners; the rerun used for the revision gave 1,513 / 407 (a three-gene difference that was not traced further).

## Licence

Code (`notebooks/`, `src/`, `scripts/`): MIT licence (`LICENSE`).
Results, tables and figures (`results/`, `supplementary/`, `figures/`, `data/`): Creative Commons Attribution 4.0 International (CC BY 4.0).

## Citation

Please cite the article once published; `CITATION.cff` will be updated with the DOI.
