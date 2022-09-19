# OncoVar

[OncoVar](https://oncovar.org/) is an integrated database and analysis platform for oncogenic driver variants in cancers.  It has mutation-level, gene-level, and pathway-level data for TCGA and ICGC mutations. Within each dataset, there is also a list of driver mutations and all mutations.  

Suggestion: Use driver mutations for enrichment; then evaluate all mutations to find other genes nearby that harbor mutations. Could also use TCGA to find clusters and validate with ICGC (or vice versa). Note that TCGA and ICGC record different cancers (though many cancers overlap).

### 1. Input Files

OncoVar is freely-accessible, and all download files are available at [this page](https://oncovar.org/welcome/download). In their [Terms of Use](https://oncovar.org/welcome/about), "Users may freely use the OncoVar database for non-commercial purposes as long as they properly cite it."  We have included the relevant input files this repo.

The input files include **gene-level** files for **All cancers** (TCGA all genes, TCGA driver genes, ICGC all genes, and ICGC driver genes). Download these four files, extract them (e.g., `tar -xvzf <FILENAME>`), and place them in this directory. They should be four directories with the names:

- `All_genes_OncoVar_ICGC`  
- `Onco_genes_OncoVar_ICGC`
- `All_genes_OncoVar_TCGA`  
- `Onco_genes_OncoVar_TCGA`

Then, use `gunzip` to unzip the compressed files in each directory:

```
gunzip */*.tsv.gz
```

_The files used in experiments were downloaded on Sept 19, 2022._

### 2. Process the OncoVar Files

```
python process_oncovar.py
```

The output files include:

- `All_ICGC.out`  
- `Onco_ICGC.out`  
- `All_TCGA.out`  
- `Onco_TCGA.out`  

All files are in the same format as the SNAP dataset.

### TODOs

- Filter cancer types that have very few genes.
- Consider PanCancer separately (or ignore it completely).

```
cut -f 1 Onco_ICGC.out | sort | uniq -c | sort -grk1 | head

686 PanCancer
 66 PRAD
 62 BRCA
 35 MALY
 32 PACA
 29 CLLE
 27 EOPC
 19 LICA
 14 PBCA
  8 OV
   ...
```
