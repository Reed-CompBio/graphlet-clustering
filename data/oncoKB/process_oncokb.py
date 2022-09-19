import sys
import csv
import os

INFILE = 'oncokb_biomarker_drug_associations.tsv'

def main():

    if not os.path.exists(INFILE):
        sys.exit("ERROR! %s does not exist in this directory. Refer to the README." % (INFILE))

    cancer_dict = read_file()

    # write to outfiles
    out = open('oncokb.out','w')
    outid = open('cancer_ids.txt','w')

    out.write('#CancerID\tCancer Name\tGene Symbol\n')
    outid.write('#CancerID\tCancerName\n')

    num = 1
    for c in cancer_dict:
        name = 'oncoKB_%d' %(num)
        num+=1
        for n in cancer_dict[c]:
            out.write('%s\t%s\t%s\n' % (name,c,n))
        outid.write('%s\t%s\n' % (name,c))


    out.close()
    outid.close()

def read_file():
    # read file
    header = None
    cancer_id = 3
    gene_id = 1
    cancer_dict = {}
    with open(INFILE) as csvfile:
        filereader = csv.reader(csvfile,delimiter='\t',quotechar='"')
        for row in filereader:
            if header == None:
                header = row
                continue
            #for i in range(len(header)):
            #    print(i,header[i])
            #print(row)

            if ',' in row[cancer_id]:
                for c in row[cancer_id].split(','):
                    c = c.strip()
                    if c not in cancer_dict:
                        cancer_dict[c] = set()
                    cancer_dict[c].add(row[gene_id])
            else:
                c = row[cancer_id].strip()
                if c not in cancer_dict:
                    cancer_dict[c] = set()
            #print()
    print('Done reading file: %d cancers mapped' % (len(cancer_dict)))

    return cancer_dict


if __name__ == '__main__':
    main()
