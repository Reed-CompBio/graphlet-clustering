import biomart
import sys
import os

NAMESPACES = ['symbol','entrez','ensembl']
ATTRIBUTES = ['uniprot_gn_symbol','entrezgene_id','ensembl_gene_id']
MAPPING_FILE = 'mapped_genes.txt'

'''
This script takes a SNAP-formatted gene set file and converts the file
into three attributes: symbol (common name), entrez ID, and ensembl.
The file must be in one of these three namespaces.

You can add any number of files, in any namespace, to process at once.
'''
def main(files):

    print('mapping %d files' % (len(files)))

    ## If the mapping file doesn't exist, use biomart to build it.
    if not os.path.isfile(MAPPING_FILE):
        print('making mapping file...')
        make_mapping_file()

    ## Convert mapping file into a dictionary of dictionaries and a
    ## dictionary of node sets.
    print('reading mapping file...')
    mapper,node_sets = read_mapping_file()

    ## process each file.
    print('processing files...')
    for f in files:
        print('  processing %s...' % (f))
        process(f,mapper,node_sets)

    return

'''
Takes a file, the mapping dictionary, and the dictionary of node sets
and writes four files: three namespace mapping files and one file of ids
dropped from the mapping.

The namespace mapping file for the SAME attribute as the original file
is different in that unmapped IDs are dropped. All namespace mapping files
have the same number of genes in each gene set.
'''
def process(file,mapper,node_sets):
    this_attribute = None
    prefix = '.'.join(file.split('.')[:-1])
    suffix = file.split('.')[-1]

    # make one outfile per namespace.
    outfiles = []
    for a in NAMESPACES:
        print('opening file:',prefix+'-'+a+'.'+suffix)
        outfiles.append(open(prefix+'-'+a+'.'+suffix,'w'))

    dropped_set = set()
    with open(file) as fin:
        for line in fin:
            row = line.strip().split('\t')
            # Write headers if this is the first line.
            if row[0][0]=='#':
                for i in range(len(NAMESPACES)):
                    outfiles[i].write('%s\t%s\t%s\n' % (row[0],row[1],NAMESPACES[i]))
                continue

            # For the first row w/ data, determine which namespace the
            # original file is in.
            if this_attribute== None:
                for a in ATTRIBUTES:
                    if row[2] in node_sets[a]:
                        this_attribute = a
            if this_attribute==None:
                for a in ATTRIBUTES:
                    print(len(node_sets[a]),row[2],row[2] in node_sets[a])
                sys.exit('ERROR: gene %s is not in a namespace.' % (row[2]))

            # for each attribute, write the appropriate line in to teh mapped outfile.
            for i in range(len(ATTRIBUTES)):
                # if this gene is not in the mapper, include it in the dropped set.
                if row[2] not in mapper[this_attribute][ATTRIBUTES[i]]:
                    dropped_set.add(row[2])
                    #print('WARNING SKIPPING BECAUSE:',row[2],'not in ',this_attribute,'->',ATTRIBUTES[i])
                else:
                    outfiles[i].write('%s\t%s\t%s\n' % (row[0],row[1],mapper[this_attribute][ATTRIBUTES[i]].get(row[2],row[2])))
    for out in outfiles:
        out.close()

    # write genes that were dropped to a droppedids file.
    dropped_out = open(prefix+'-droppedids.'+suffix,'w')
    for d in dropped_set:
        dropped_out.write(d+'\n')
    dropped_out.close()

    print('done writing mapping files.')
    return

## Use BioMart to get mapping among the three namespaces.
def make_mapping_file():
    # lines come from these two code snippets
    # https://pypi.org/project/biomart/
    # https://gist.github.com/ben-heil/cffbebf8865795fe2efbbfec041da969
    ## Values for the mirror are: useast, uswest, asia, and www.
    server = biomart.BiomartServer( "http://useast.ensembl.org/biomart" )
    server.verbose = True

    mart = server.datasets['hsapiens_gene_ensembl']
    response = mart.search({'attributes':ATTRIBUTES})
    data = response.raw.data.decode('ascii')

    out = open(MAPPING_FILE,'w')
    out.write('#'+'\t'.join([NAMESPACES[n] for n in ATTRIBUTES])+'\n')
    print('\n')
    print(ATTRIBUTES)
    t = 1
    for line in data.splitlines():
        row = line.strip().split('\t')
        if len(row) < 3:
            continue
        if row[0] == "" or row[1]=="" or row[2] == "": # keep genes that have have full mapping
            continue
        if 'MT-' in row[0]: # skip mitochondrial (MT-) genes
            continue

        out.write(line+'\n')
    print('Done writing mapping file')
    return

## If mapping file exists, read it and build dictionaries of mappers.
def read_mapping_file():
    mapper = {}
    node_sets = {}
    for a in ATTRIBUTES:
        mapper[a] = {}
        node_sets[a] = set()
        for b in ATTRIBUTES:
            mapper[a][b] = {}

    with open(MAPPING_FILE) as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            for i in range(len(ATTRIBUTES)):
                node_sets[ATTRIBUTES[i]].add(row[i])
                for j in range(len(ATTRIBUTES)):
                    mapper[ATTRIBUTES[i]][ATTRIBUTES[j]][row[i]] = row[j]
    print('finished reading mapped_genes.txt as dictionary of dictionaries')
    for k1 in mapper:
        for k2 in mapper[k1]:
            print(k1,'->',k2)
    return mapper,node_sets

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('ARGUMENTS:',sys.argv)
        sys.exit('USAGE: python map_files.py <FILE1> [<FILE2> <FILE3> ...]\n')
    main(sys.argv[1:])
