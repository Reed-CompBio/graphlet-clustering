## Mapping Namespace

This directory has a `map-files.py` script that maps any SNAP-formatted file into three namespaces: `symbol`,`entrez`, and `ensembl`.  The original file **must** be in one of these three namespaces.

It uses the `biomart` module to retrieve this mapping from Ensembl, in the case that the file `mapped_genes.txt` does not exist.

### Usage:

```
python map_files.py ../snap/ppi_DG.tsv
```

This will produce the following files:
- `../snap/ppi_DG-droppedids.tsv`
- `../snap/ppi_DG-ensembl.tsv`    
-  `../snap/ppi_DG-entrez.tsv`
- `../snap/ppi_DG-symbol.tsv`

Note: you can add any number of SNAP-formatted gene sets as arguments, in any of the three namespaces allowed.

### BioMart attributes

The attributes in the BioMart human ensembl database are listed in `BioMart-hsapiens_gene_ensembl-attributes.txt`.  Note that the mirror is set to `useast` when accessing BioMart; the `uswest` mirror was down.

BioMart mapping completed on September 29, 2022.
