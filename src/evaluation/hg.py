import numpy as np
import networkx as nx
import pandas as pd
import os
#import os.path
import sys
from csv import reader
from networkx.algorithms import approximation
from matplotlib import pyplot as plt
#from sklearn.metrics.cluster import adjusted_rand_score
#from sklearn.preprocessing import MultiLabelBinarizer
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection


def get_gene_population(outdir):
	nodefile = outdir+'g_network_nodeid.txt'
	genes = pd.read_csv(nodefile,sep='\t',header=None,dtype=str)
	gene_pop = list(set(genes[1]))
	return(gene_pop)

def make_disease_gene_dictionary(interactome_dg_file, gene_pop):
	G_initial = interactome_dg_file
	ith_clusters = []
	ith_lengths = []

	Gi = open(G_initial, "r")
	ppi_dg_dict = {}
	lines = Gi.readlines()[1:]
	for line in lines:
		stripped_line = line.strip()
		line_list = stripped_line.split()
		disease_id = line_list[0]
		gene_id = line_list[-1]
		if gene_id in gene_pop:
			if disease_id in ppi_dg_dict:
				updated_list = ppi_dg_dict[disease_id]
				updated_list.append(gene_id)  
				ppi_dg_dict[disease_id] = updated_list
			else:
				ppi_dg_dict[disease_id] = [gene_id]
	Gi.close()
	return(ppi_dg_dict)

def make_disease_name_dictionary(interactome_dg_file):
	Gi = pd.read_csv(interactome_dg_file,sep='\t')
	Gi.drop(Gi.iloc[:,[2]],axis = 1,inplace=True)
	Gi.drop_duplicates(inplace=True)
	dg_name_dict = dict(zip(Gi.iloc[:,0],Gi.iloc[:,1]))
	return(dg_name_dict)

def convert_to_integers(unconverted,interactome_mcl_dir):
	converted = []
	G_initial = interactome_mcl_dir+'g_network_nodeid.txt'
	Gi = open(G_initial, "r")
	for line in Gi:
		stripped_line = line.strip()
		pair_list = stripped_line.split()

		if pair_list[1] in unconverted:
			converted.append(pair_list[0])
	Gi.close()
	return(converted)

def get_all_clusters(outdir,inflation):
	all_clusters = []
	cluster_lengths = []
	
	for i in range(1):

		G_initial = outdir+'clusters_'+str(inflation)+'_G'+str(i)+'.txt' #'inflation_tests/G'+str(i)+'_I3.0'
		ith_clusters = []
		ith_lengths = []
	
		Gi = open(G_initial, "r")
		for line in Gi:
			
			stripped_line = line.strip()
			line_list = stripped_line.split()
			
			ith_clusters.append(line_list)
			ith_lengths.append(len(line_list))

		Gi.close()

		#add all the nodes taken off by the graphlet step
		all_clusters.append(ith_clusters)
		cluster_lengths.append(ith_lengths)
	return(all_clusters,cluster_lengths)


# Given all of the clusters and the target nodes, it will return a list of lists of the number of target nodes
# in each cluster
def get_target_counts(all_clusters,target_nodes):

	all_target_counts = []


	for i in range(len(all_clusters)):

		ith_target_counts = []
	

		for j in range(len(all_clusters[i])):

			count = 0
			for target in target_nodes:
				if target in all_clusters[i][j]:
					count += 1
					
			
			#if count>0:
			ith_target_counts.append(count)

		
		if len(ith_target_counts)<0:
			ith_target_counts = [None]
	

		all_target_counts.append(ith_target_counts)

	return(all_target_counts)

'''
def list_of_strings(list,prefix='',suffix=''):
	return_list = []
	for item in list:
		return_list.append(prefix+str(item)+suffix)
	return return_list

# a max function that also returns the index of the max value in the list
def max_returns_index(ranlis):
	mx = 0
	mx_index = 0
	for i in range(len(ranlis)):
		if ranlis[i] > mx:
			mx = ranlis[i]
			mx_index = i
	return(mx,mx_index)
'''

def perform_hg_test(dg_file, outdir, inflation, dg_file_name):
	all_clusters,cluster_lengths = get_all_clusters(outdir,inflation)


	gene_pop = get_gene_population(outdir)
	total_genes = len(gene_pop)

	disease_gene_dict = make_disease_gene_dictionary(dg_file,gene_pop) #key: disease id, value: list of genes
	disease_name_dict = make_disease_name_dictionary(dg_file) #key: disease id, value: disease name

	ids_lis = []
	for ids in disease_gene_dict:
		ids_lis.append(ids)

	all_counts = []

	for disease_id in disease_gene_dict:
		target_nodes_unconverted = disease_gene_dict[disease_id]
		target_nodes = convert_to_integers(target_nodes_unconverted,outdir)
		target_counts = get_target_counts(all_clusters,target_nodes)
		all_counts.append(target_counts)


	dg_list_sizes = []
	for ids in disease_gene_dict:
		dg_list_sizes.append(len(disease_gene_dict[ids]))

	

	n_d = len(dg_list_sizes)

	for i in range(1):
		f = open(outdir+'enrichment_modules_{}_G{}_{}.txt'.format(inflation,i,dg_file_name),'w')
		for j in range(n_d):
			#print(j)
			target_counts = all_counts[j]

			num_of_poss_dg = dg_list_sizes[j]
				#total_dg_by_graphlet[i][j]
		   
			if len(target_counts[i]) > 0:
					#mx, mx_index = max_returns_index(target_counts[i])
					#cluster_index_array[i][j] = mx_index
					#P[i][j] = max(target_counts[i])
					for k in range(len(target_counts[i])):
						correct_dg=target_counts[i][k]
						chosen_cluster_size = len(all_clusters[i][k])
						
						#if (correct_dg>0): # this will bias the results
						if (chosen_cluster_size>2): # include only clusters >= 3
							#dis_id #genes in disease cluster #graphlet #cluster_id #cluster size #genes in common #p-value
							#f.write('j '+str(j)+' '+str(num_of_poss_dg)+' i '+str(i)+' k '+str(k)+' common '+str(correct_dg)+' '+str(hypergeom.sf(correct_dg-1, total_genes, num_of_poss_dg, chosen_cluster_size))+'\n')
							f.write(disease_name_dict[ids_lis[j]]+'\t'+str(num_of_poss_dg)+'\t'+str(i)+'\t'+str(k)+'\t'+str(chosen_cluster_size)+'\t'+str(correct_dg)+'\t'+str(hypergeom.sf(correct_dg-1, total_genes, num_of_poss_dg, chosen_cluster_size))+'\n')
		f.close()

def perform_enrichment(outdir, dg_file_name, inflation):

	t=0.05 #significance threshold

	with open(outdir+'NS_{}.txt'.format(dg_file_name), 'w') as fp1:
		for i in range(1):
			df = pd.read_csv(outdir+'enrichment_modules_{}_G{}_{}.txt'.format(inflation,i,dg_file_name),
							 sep='\t',header=None)
			df[7] = fdrcorrection(df[6],t)[1]
			df[8] = fdrcorrection(df[6],t)[0]
			df.sort_values(by=7,inplace=True)
			df=df[df[8]==True]
			df.drop(columns={2,6,8},inplace=True) #2: graphlet
			
			x=len(np.unique(df[3]))#sig. clusters
			y=len(np.unique(df[0]))#sig. diseases

			fp1.write('G'+str(i)+'\t'+str(x)+'\t'+str(y)+'\n')

			df=df.rename(columns={0: 'Disease', 1:'DiseaseGenes',3:'ClusterID',4:'ClusterSize',5:'Common',7:'P-value'})
			df.to_csv(outdir+'enrichment_modules_{}_G{}_{}.csv'.format(inflation,i,dg_file_name),sep='\t',index=False)



	df = pd.read_csv(outdir+'NS_{}.txt'.format(dg_file_name),sep='\t',header=None)

	plt.clf()
	plt.bar(df[0],df[1])
	plt.xlabel('Graphlet')
	plt.xticks(rotation=90)
	plt.ylabel("No. of significant modules")
	plt.savefig(outdir+'sig_modules_{}.pdf'.format(dg_file_name),format='pdf',bbox_inches='tight')

	plt.clf()
	plt.bar(df[0],df[2])
	plt.xlabel('Graphlet')
	plt.xticks(rotation=90)
	plt.ylabel("No. of significantly associated diseases")
	plt.savefig(outdir+'sig_diseases_{}.pdf'.format(dg_file_name),format='pdf',bbox_inches='tight')


def main(argv):
	outdir = argv[1] 
	dg_file = argv[2]
	inflation = argv[3]

	
	#if not os.path.exists(outdir):
	#	os.makedirs(outdir)

	#dg_file = '../../data/snap/ppi_DG.tsv'
	#outdir = '../../out/PP-Pathways_ppi/'

	dg_file_name = dg_file.split('/')[-1]#.split('.')[0]

	perform_hg_test(dg_file,outdir,inflation,dg_file_name)
	perform_enrichment(outdir, dg_file_name,inflation)


if __name__ == "__main__":
	main(sys.argv)