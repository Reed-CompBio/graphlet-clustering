import sys
import csv
import os

INFILE = 'cancer_gene_census.csv'
ABBREVFILE = 'abbrevs.txt'

def main():

    if not os.path.exists(INFILE):
        sys.exit("ERROR! %s does not exist in this directory. Refer to the README." % (INFILE))

    abbrevs = read_abbrevs()
    mapping_dict,tumor_dict,cancers = read_file()
    #for c in cancers:
    #    print('*',c,'*')
    cancer_dict = map_by_cancer_types(tumor_dict,cancers)

    # write to outfiles
    out1 = open('CGC_tier1.csv','w')
    out2 = open('CGC_tier2.csv','w')
    outboth = open('CGC_bothtiers.csv','w')
    outid = open('cancer_ids.txt','w')

    for out in [out1,out2,outboth]:
        out.write('#CancerID\tCancer Name\tGene Symbol\n')

    outid.write('#CancerID\tCancerName\n')

    num = 1
    for c in cancer_dict:
        name = 'CGC_%d' %(num)
        num+=1
        for tier in ['1','2']:
            for n in cancer_dict[c][tier]:
                if tier == '1':
                    out1.write('%s\t%s\t%s\n' % (name,abbrevs.get(c,c),n))
                elif tier == '2':
                    out2.write('%s\t%s\t%s\n' % (name,abbrevs.get(c,c),n))
                outboth.write('%s\t%s\t%s\n' % (name,abbrevs.get(c,c),n))
        outid.write('%s\t%s\n' % (name,abbrevs.get(c,c)))

    for out in [out1,out2,outboth,outid]:
        out.close()

def read_abbrevs():
    abbrevs = {}
    with open(ABBREVFILE) as fin:
        for line in fin:
            row = line.strip().split()
            abbrevs[row[0]] = ' '.join(row[1:])
    return abbrevs

def map_by_cancer_types(tumor_dict,cancers):
    cancer_dict = {}
    for c in cancers:
        cancer_dict[c] = {'1':set(),'2':set()}
    for n in tumor_dict:
        for tier in ['1','2']:
            for x in tumor_dict[n][tier]:
                cancer_dict[x][tier].add(n)
    return cancer_dict

def read_file():
    # read file
    header = None
    map_inds = [1,2,19,21]
    tumor_ind = 9
    tier_ind = 4
    mapping_dict = {}
    tumor_dict = {}
    with open(INFILE) as csvfile:
        filereader = csv.reader(csvfile,delimiter=',',quotechar='"')
        for row in filereader:
            if header == None:
                header = row
                continue
            #for i in range(len(header)):
                #print(i,header[i])
            #print(row)
            mapping_dict[row[0]] = [row[i] for i in map_inds]
            if row[0] not in tumor_dict:
                tumor_dict[row[0]] = {'1':[],'2':[]}
            if ',' in row[tumor_ind]:
                tumor_dict[row[0]][row[tier_ind]] = [a for a in ','.split(row[tumor_ind]) if a not in ['',' ',',']]
            elif row[tumor_ind] not in ['',' ',',']:
                tumor_dict[row[0]][row[tier_ind]] = [row[tumor_ind]]
            #print()
    print('Done reading file: %d genes mapped; %d genes associated with cancers' % (len(mapping_dict),len(tumor_dict)))

    cancers = set()
    for n in tumor_dict:
        cancers.update(tumor_dict[n]['1'])
        cancers.update(tumor_dict[n]['2'])
    print('%d distinct cancer types that harbor somatic mutations' % (len(cancers)))

    return mapping_dict,tumor_dict,cancers


if __name__ == '__main__':
    main()
