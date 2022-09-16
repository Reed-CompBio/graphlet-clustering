# COSMIC Cancer Database: Cancer Gene Consensus (CGC)

The Catalogue of Somatic Mutations in Cancer (COSMIC) database contains a well-curated list of genes whose mutations are associated with different cancers.  This folder processes the downloaded CGC file into a format that can be used for graphlet clustering.

### 1. Get the CGC File

Go to the [GCG Website](https://cancer.sanger.ac.uk/census) and download the `.csv` file that includes both tiers of annoataions (Tier 1 & Tier 2). Note that you must make an account in order to download this file for academic purposes. The file should be named `cancer_gene_census.csv` and placed in this directory.

_The file used for experiments was downloaded on Sept 16, 2022._

### 2. Process the CGC File

The `csv` module is required.

```
python process_cgc.py
```

The output files include: `CGC_tier1.csv`,`CGC_tier2.csv`, and `CGC_bothtiers.csv`. They are in the same format as the SNAP dataset.

### TODOs

1. Map the abbreviated names to full names
2. Filter cancer types that have very few genes.

```
cut -f 1 CGC_bothtiers.csv | sort | uniq -c | sort -grk1
  39 AML
  19 T-ALL
  14 prostate
  12 ALL
  11 melanoma
  11 NHL
   9 Spitzoid tumour
   9 AL
   8 colorectal
   8 NSCLC
   7 papillary thyroid
   7 DLBCL
   7 ALCL
   6 MM
   6 AML*
   5 synovial sarcoma
   5 lung cancer
   5 glioma
   5 breast
   5 aneurysmal bone cyst
   5 B-NHL
   ...
```
