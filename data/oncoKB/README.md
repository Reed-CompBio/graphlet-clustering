# OncoKB

OncoKB is a knowledgebase of actionable genes and cancer genes developed at Memorial Sloan Kettering Cancer Center. This folder processes the downloaded OncoKB files into a format that can be used for graphlet clustering.

### 1. Get the OncoKB File

Go to the [OncoKB Website](https://www.oncokb.org/) and make an account to access the data. Note that it takes 1-2 business days to verify your account.  Download the Associations file from [Actionable Genes](https://www.oncokb.org/actionableGenes#sections=Tx), which includes 5 levels (listed at the top of that webpage). The file should be named `oncokb_biomarker_drug_associations.tsv` and placed in this directory.

_The file used for experiments was downloaded on Sept 19, 2022._

### 2. Process the OncoKB File

The `csv` module is required.

```
python process_oncokb.py
```

The output files include: `oncokb.out` (same format as SNAP dataset) and `cancer_ids.txt` (maps the cancer ID to the cancer name).

### TODOs

- Filter cancer types that have very few genes.
- There is also a `cancerGeneList.txt` file, of 1,067 cancer-associated genes (not broken down by cancer types). Maybe this is a good 'PanCancer' example?

```
cut -f 2 oncokb.out | sort | uniq -c | sort -grk1

28 NOS
14 Prostate Cancer
11 Diffuse Large B-Cell Lymphoma
 9 BCR-ABL1 Like
 9 B-Lymphoblastic Leukemia/Lymphoma
 3 with MYC and BCL2 and/or BCL6 Rearrangements
 3 Pancreatic Adenocarcinoma
 3 Ovarian Cancer
 3 High-Grade B-Cell Lymphoma
 3 Acinar Cell Carcinoma of the Pancreas
 ...
```
