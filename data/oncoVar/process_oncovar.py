import sys
import os
import glob

DIR_INFIX = '_genes_OncoVar_'

def main():

    for scope in ['Onco','All']:
        for dataset in ['ICGC','TCGA']:
            indir = '%s%s%s' % (scope,DIR_INFIX,dataset)
            outfile = '%s_%s.out' % (scope,dataset)

            if not os.path.exists(indir):
                sys.exit("ERROR! Directory %s does not exist in this directory. Refer to the README." % (indir))

            process_files(indir,outfile)



def process_files(indir,outfile):
    out = open(outfile,'w')
    out.write('#CancerCode\tCancerName\tGene\n')
    for file in glob.glob(indir+'/*'):
        print(file)
        header = None
        with open(file) as fin:
            for line in fin:
                if header == None:
                    header = line.strip().split('\t')
                    continue

                row = line.strip().split('\t')
                if 'PanCancer' in row[0]:
                    out.write('%s\t%s\t%s\n'% (row[0],row[0],row[1]))
                else:
                    split_name = row[0].split('[')
                    full = split_name[0].strip()
                    abbrev = split_name[1][:-1]
                    out.write('%s\t%s\t%s\n'% (abbrev,full,row[1]))
    out.close()
    print('Wrote to %s' % (outfile))
    return

if __name__ == '__main__':
    main()
