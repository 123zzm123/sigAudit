# Data

Raw expression matrices and survival tables are **not redistributed** in this repository.

## TCGA

cBioPortal study identifiers (all suffix `_tcga_pan_can_atlas_2018`):

`brca_tcga_pan_can_atlas_2018`, `coadread_tcga_pan_can_atlas_2018`,
`kirc_tcga_pan_can_atlas_2018`, `lgg_tcga_pan_can_atlas_2018`,
`luad_tcga_pan_can_atlas_2018`, `meso_tcga_pan_can_atlas_2018`,
`paad_tcga_pan_can_atlas_2018`, `stad_tcga_pan_can_atlas_2018`

- Molecular profile: `{study}_rna_seq_v2_mrna`
- Sample list: `{study}_all`
- API: `https://www.cbioportal.org/api`

## METABRIC

- Study identifier: `brca_metabric`
- Molecular profile: `brca_metabric_mrna`
- Same API.

## Gene sets

MSigDB v2024.1, six collections: `h.all`, `c2.cp`, `c2.cp.kegg_medicus`,
`c2.cp.reactome`, `c5.go.bp`, `c7.immunesigdb`
(`*.v2024.1.Hs.symbols.gmt`, downloaded from `https://www.gsea-msigdb.org/`).

## Signature gene lists

`genefu` 2.44.0 (Bioconductor). The six objects used and their MD5 checksums are recorded in
`third_party/genefu/README.md`; the objects themselves are not redistributed.

## Where each artefact lives

| Artefact | Location |
|---|---|
| Raw expression and survival tables | fetched by the user into `data/cohorts/` |
| Per-patient predictions | produced by `reproduce_current.sh train`; not in this repository |
| Derived result tables | `results/tables/`, `results/manifest/` |
| Analysis code | `scripts/` |
