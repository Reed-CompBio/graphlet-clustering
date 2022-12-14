# graphlet-clustering

This repo contains code for graphlet-based community detection in Protein-Protein Interaction (PPI) netowrks. 

*Identification of Disease Modules Using Higher-Order Network Structure,*
Pramesh Singh, Hannah Kuder, Anna Ritz
https://www.biorxiv.org/content/10.1101/2022.12.24.521876v1

## Directory Contents

- **src**: source code.

- **data**: Interactomes, biological pathway gene sets, disease gene sets, DREAM challenge datasets and scoring tools.


## Orca

Our method code uses ORbit Counting Algorithm (ORCA) for computing edge-orbits. To use ORCA run the following commands:

```
git submodule init
git submodule update
cd src/clustering/orca/
g++ -O2 -std=c++11 -o orca.exe orca.cpp
```

See the `README.md` file in `orca/` for more info.
