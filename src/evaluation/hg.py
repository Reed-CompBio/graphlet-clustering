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


def make_disease_gene_dictionary(interactome_dg_file):
	Gi = open(interactome_dg_file, "r")
	ppi_dg_dict = {}
	for line in Gi:
		stripped_line = line.strip()
		line_list = stripped_line.split()
		disease_id = line_list[0]
		gene_id = line_list[-1]
		
		if disease_id in ppi_dg_dict:
			updated_list = ppi_dg_dict[disease_id]
			updated_list.append(gene_id)  
			ppi_dg_dict[disease_id] = updated_list
		else:
			ppi_dg_dict[disease_id] = [gene_id]
	Gi.close()
	ppi_dg_dict.pop('#');
	return(ppi_dg_dict)

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

def get_all_clusters(outdir):
	all_clusters = []
	cluster_lengths = []
	
	for i in range(30):

		G_initial = outdir+'clusters_G'+str(i)+'.txt' #'inflation_tests/G'+str(i)+'_I3.0'
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


def main(argv):
	outdir = argv[1] 
	dg_file = argv[2]
	inflation = 4.0

	
	#if not os.path.exists(outdir):
	#	os.makedirs(outdir)

	#dg_file = '../../data/snap/ppi_DG.tsv'
	#outdir = '../../out/PP-Pathways_ppi/'

	disease_gene_dict = make_disease_gene_dictionary(dg_file) #key: disease, value: list of genes
	all_clusters,cluster_lengths = get_all_clusters(outdir)

	ids_lis = []
	for ids in disease_gene_dict:
		ids_lis.append(ids)

	all_counts = []

	for disease_id in disease_gene_dict:
		target_nodes_unconverted = disease_gene_dict[disease_id]
		target_nodes = convert_to_integers(target_nodes_unconverted,outdir)
		target_counts = get_target_counts(all_clusters,target_nodes)
		all_counts.append(target_counts)


	P = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	significant_array = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	cluster_index_array = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	significant_labels = []
	sig_counts = []

	for j in range(519):
		#print(j)
		disease_counts_by_graphlet = []
		target_counts = all_counts[j]
		for i in range(30):
			if len(target_counts[i]) > 0:
					mx, mx_index = max_returns_index(target_counts[i])
					cluster_index_array[i][j] = mx_index
					P[i][j] = max(target_counts[i])
					
					disease_counts_by_graphlet.append(max(target_counts[i]))
		if max(disease_counts_by_graphlet) > 4:
			sig_counts.append(target_counts)
			significant_labels.append(list(disease_gene_dict)[j])
			for i in range(30):
				if len(target_counts[i]) > 0:
					significant_array[i][j] = max(target_counts[i])

	H = np.zeros(519*30).reshape(30,519)

	dg_list_sizes = []
	for ids in disease_gene_dict:
		dg_list_sizes.append(len(disease_gene_dict[ids]))

	for j in range(519):
		for i in range(30):

			num_of_poss_dg = dg_list_sizes[j]
			#total_dg_by_graphlet[i][j]
			
			total_genes = 21538
			#total_genes_by_graphlet[i]

			cor_cluster = int(cluster_index_array[i][j]) #corresponding cluster
			chosen_cluster_size = len(all_clusters[i][cor_cluster])

			correct_dg = P[i][j]

			H[i][j] = hypergeom.sf(correct_dg-1, total_genes, num_of_poss_dg, chosen_cluster_size)

	row_names = list_of_strings(list(range(30)),'G')
	#row_names.append('DG_size')
	#H = np.vstack([H,dg_list_sizes])    
	column_names = ids_lis
	pd.DataFrame(H,columns=column_names,index = row_names).to_csv(outdir+'ppi_dg_hypergeo3_new.csv')

if __name__ == "__main__":
	main(sys.argv)