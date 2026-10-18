# Data

The raw expression matrices and survival tables are **not redistributed** here.
TCGA Pan-Cancer Atlas 2018 and METABRIC data are available from cBioPortal:

- TCGA: `https://www.cbioportal.org/api` — eight study identifiers are listed in Table 1 of
  the manuscript (`brca_tcga_pan_can_atlas_2018`, `coadread_...`, `kirc_...`, `lgg_...`,
  `luad_...`, `meso_...`, `paad_...`, `stad_...`), profile `*_rna_seq_v2_mrna`.
- METABRIC: study `brca_metabric`, profile `brca_metabric_mrna`.

Gene sets: MSigDB v2024.1 (`https://www.gsea-msigdb.org/`), six collections listed in
Supplementary Methods S1.

Signature gene lists were read from `genefu` 2.44.0 (Bioconductor) and are not redistributed;
`third_party/genefu/README.md` records the objects used and their MD5 checksums.

After download, place files under `data/cohorts/` as described in `reproduce_current.sh`.
