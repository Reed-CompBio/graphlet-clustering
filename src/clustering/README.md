# graphlet-clustering

This code finds clusters in a given PPI network with respect to specific graphlets. There are 30 graphlets of size 5 or smaller. Thus using the code we can find 30 different partitions of the same network.

The code `graphlet_clustering.py` can be run with 5 arguments:

- path to the PPI network. This is file should be a tab separated 2 (unweighted) or 3 (weighted) column edgelist.
- path to output directory. This is where all the output files will be saved. If output directory does not exist, it will be created by this code.
- Graphlet size. This argument is either `4` or `5`.
- Enter `1` for weighted or `0` for unweighted.
- Inflation parameter. This controls the granularity of clusters.

### Example:
`python graphlet_clustering.py ../../data/DREAM/1_networks/original/2_ppi_inweb.txt ../../out/2_ppi_inweb/ 5 1 4.0`

## Dependencies
### ORCA
Orca is an efficient algorith for computed node and edge orbits (See https://doi.org/10.1093/bioinformatics/btt717 ). The code makes use of orca code, which needs to be compiled to generate an executable `orca.exe` in `orca/`. Here orca is added as a submodule. To compile:

`cd orca`
`g++ -O2 -std=c++11 -o orca.exe orca.cpp`

See the `README.md` file in `orca/` for more information.

### Markov Clustering Algorithm (MCL) 
MCL is random-walk based clustering algorithm that our code uses (See http://micans.org/mcl/ ). `mcl` needs to be installed before graphlet clustering can be used ( https://github.com/micans/mcl ). 
On macOS, it can be installed easily through Mac Ports ( https://ports.macports.org/port/mcl/ ).
