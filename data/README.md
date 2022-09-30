# Datasets

Datasets with asterisks require registration to access the data or have some other data limitations (licensing, size, etc.).  Instructions for accessing the data are provided in each subdirectory.

### Interactomes

These are undirected but potentially weighted protein-protein interaction networks.

- `interactomes/All_Pathway_Commons.txt` [namespace: *Symbol*]
- `interactomes/PP-Pathways_ppi.csv`: protein interaction network from SNAP ([PP-Pathways](https://snap.stanford.edu/biodata/datasets/10000/10000-PP-Pathways.html)) [namespace: *Entrez*]
- `interactomes/huri.tsv`: Human Reference Interactome ([HuRI](http://www.interactome-atlas.org/)) [namespace: *Ensembl*]

Additionally, we use the following two interactomes provided with the DREAM challenge:

- `DREAM/1_networks/original/1_ppi_string_v2.txt` [namespace: *Symbol*]
- `DREAM/1_networks/original/2_ppi_inweb_v2.txt` [namespace: *Symbol*]

- **TODO fill in**

### Biological Process Gene Sets

- `MSigDB/hallmark/`: hallmark gene sets from MSigDB*
- `MSigDB/canonical_BioCarta/`: canonical pathway gene sets from BioCarta*
- `MSigDB/canonical_KEGG/`: canonical pathway gene sets from KEGG*
- `MSigDB/canonical_PID/`: canonical pathway gene sets from NCI-PID*

Once these are processed, see instructions in the `mappers/` directory to map to different namespaces.

### Disease Gene Sets

Once these are processed, see instructions in the `mappers/` directory to map to different namespaces.

#### Multiple Diseases

- `snap/`: disease-gene association network from [DG-Association Miner](http://snap.stanford.edu/biodata/datasets/10012/10012-DG-AssocMiner.html) **TODO license info?? can we commit this to the repo??**

#### Cancer

- `cosmic/`: cancer gene sets (level 1 and 2) from the COSMIC dataset (very few gene sets)*
- `MSigDB/oncogenic/`: gene sets of pathway signatures dysregulated in cancer*
- `oncoKB/`: cancer gene sets from the OncoKB dataset (very few gene sets)*
- `OncoVar/`: cancer gene sets from TCGA and ICGC curated by OncoVar*
