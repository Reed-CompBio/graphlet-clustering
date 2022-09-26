# MSigDB

GSEA's [Molecular Signatures Database](http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) contains collections of gene sets.  

### Get the Datasets

[Datasets are available for download](http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) after registration.  Download the "Gene Symbols" and the "Entrez Gene IDs" GMT files and place them in the respective sub-directories:

- `MSigDB/hallmark/`: hallmark gene sets from MSigDB
- `MSigDB/canonical_BioCarta/`: canonical pathway gene sets from BioCarta
- `MSigDB/canonical_KEGG/`: canonical pathway gene sets from KEGG
- `MSigDB/canonical_PID/`: canonical pathway gene sets from NCI-PID
- `MSigDB/oncogenic/`: gene sets of pathway signatures dysregulated in cancer


### Process Files

Within _this_ directory, call

```
python gmt-to-snap.py
```

This will convert **every** `*.gmt` file in a subdirectory into an `*.out` file that is formatted like the SNAP dataset.  The GeneSet ID column (column 1) includes the shortened abbreviation for each curated subset:

- H:hallmark
- C2: curated gene sets (BioCarta, KEGG, and PID are subsets)
- C6: oncongenic gene sets
