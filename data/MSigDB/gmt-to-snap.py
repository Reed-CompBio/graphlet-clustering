import glob

def main():

    files = glob.glob('*/*.gmt')
    for f in files:
        outf = f.replace('.gmt','.out')
        print('converting %s to %s...'% (f,outf))
        out = open(outf,'w')
        out.write('#GeneSet ID\tGeneSet Name\tGene\n')

        idprefix = f.split('.')[0].split('/')[-1]
        with open(f) as fin:
            id = 1
            for line in fin:
                row = line.strip().split('\t')
                geneset_name = row[0]
                for gene in row[2:]:
                    out.write('%s-%d\t%s\t%s\n' % (idprefix,id,geneset_name,gene))
                id+=1
        out.close()

if __name__ == '__main__':
    main()
