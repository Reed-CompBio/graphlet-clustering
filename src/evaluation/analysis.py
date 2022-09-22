# argument 1: path to the directory with clusters_G0 - G29.txt files 
# relative to the current directory.
#
# argument 2: path to save the output 
#
import pandas as pd
import numpy as np
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score


def get_ari_matrix(clusterdir, t):

	all_clusters=None
	for g in range(0,30):

		clusters={}
		cluster_file = clusterdir+'clusters_G'+str(g)+'.txt'
		f2=open(cluster_file,'r')
		lines = f2.read().splitlines()
		i=0
		for line in lines:
			#if (len(line.split('\t'))>=3): #clusters larger than 2
			for l in line.split('\t'):
				clusters[l]=i
			i=i+1
		f2.close()
		
		
		if (all_clusters is None):
			all_clusters=pd.DataFrame([clusters.keys(),clusters.values()]).T
			all_clusters.rename(columns={1:'G'+str(g),0:'node'},inplace=True)
			all_clusters['freq''G'+str(g)] = all_clusters.groupby('G'+str(g))['G'+str(g)].transform('count')
			all_clusters.reset_index(drop=True)
			#print(len(all_clusters))
		else:
			c2=pd.DataFrame([clusters.keys(),clusters.values()]).T
			c2.rename(columns={1:'G'+str(g),0:'node'},inplace=True)
			c2['freq''G'+str(g)] = c2.groupby('G'+str(g))['G'+str(g)].transform('count')
			c2.reset_index(drop=True)
			all_clusters = pd.merge(all_clusters,c2,on=['node'],how='inner')
			all_clusters=all_clusters.drop_duplicates(subset=['node'])
			#print(len(all_clusters))

	ari_matrix = np.zeros(30*30).reshape(30,30)
	for i in range(30):
		for j in range(30):
			x = all_clusters[(all_clusters['freqG'+str(i)]>=t) & (all_clusters['freqG'+str(j)]>=t)]['G'+str(i)]
			y = all_clusters[(all_clusters['freqG'+str(i)]>=t) & (all_clusters['freqG'+str(j)]>=t)]['G'+str(j)]
			ari_matrix[i][j]=adjusted_rand_score(x,y)
	return ari_matrix

def plot_ari(savepath,ari_matrix,plottitle):
	plt.clf()
	plt.matshow(ari_matrix)
	plt.colorbar()
	plt.title(plottitle)
	plt.savefig(savepath)


# Hannah's functions



# Parameter Sweep functions


# Runs MCL for multiple inflation values
def run_mcl_multiple_inflations(outdir, inflations, graphlet_list):
	# outdir = '../../out/1_ppi_string/'
	# inflations = [1.2, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]
	# graphlet_list = [0,1,2]

	for graphlet in graphlet_list:
		for inflation in inflations:

			str_inflation = str(inflation)
			new_file_name = "{}clusters_G{}_.txt".format(outdir,graphlet,str_inflation)

			if os.path.exists(new_file_name) == False:
				print("Running mcl for G"+str(graphlet))
				cmd = "mcl {}network_G{}.txt --abc -I {} -o {}clusters_G{}.txt".format(outdir,graphlet,inflation,outdir,graphlet)
				os.system(cmd)

# Returns the largest cluster size and number of clusters for the chosen graphlets and inflations
def parameter_lists(outdir, inflations, graphlet_list):
	# outdir = '../../out/1_ppi_string/'
	# inflations = [1.2, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]
	# graphlet_list = [0,1,2]

	ppi_largest_cluster_lengths = []
	ppi_number_of_clusters = []

	for graphlet in graphlet_list:
		ith_graphlet_lengths = []
		ith_cluster_lengths = []
		for inflation in inflations:

			str_inflation = str(inflation)
			new_file_name = "{}clusters_G{}_.txt".format(outdir,graphlet,str_inflation)
			

			list_of_lists = []
			cluster_num = 0
			Gi = open(new_file_name, "r")
			for line in Gi:
				cluster_num += 1
				stripped_line = line.strip()
				line_list = stripped_line.split()
				list_of_lists.append(line_list)
			Gi.close()

			largest_cluster_lenth = len(list_of_lists[0])
			ith_graphlet_lengths.append(largest_cluster_lenth)

			ith_cluster_lengths.append(cluster_num)
		ppi_largest_cluster_lengths.append(ith_graphlet_lengths)
		ppi_number_of_clusters.append(ith_cluster_lengths)

def parameter_sweep_plot(title, num_of_graphlets, parameter_list, ylabel, inflations):

	fig = plt.figure()
	ax1 = fig.add_subplot(111)

	inflations = [1.2, 1.5, 2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 5.0]
	colors = ['k','dimgray','lime','rosybrown','maroon','gold',
				'chocolate', 'fuchsia','orangered','saddlebrown',
				'darkorange','darkkhaki','olive','darkolivegreen',
				'chartreuse','palegreen','turquoise','cornflowerblue','blue',
				'blueviolet','violet','crimson','deeppink',
				'darkslategray','yellow','gray','cyan','brown','olive',
				'deepskyblue'
				]
	#,'pink','lightcoral','darkgray',

	for graphlet in range(int(num_of_graphlets)):
		ax1.scatter(inflations, parameter_list[graphlet], s=10, c=colors[graphlet], marker=".", label='G'+str(graphlet))
			


	#ax1.set_yscale('log')

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5));
	ax1.set_ylabel(ylabel)
	ax1.set_xlabel('MCL inflation paramter')
	ax1.set_title(title)
	plt.show()

# Prints the edges missing by each graphlet, see how different the networks are after you remove non graphlet edges
# outdir is where the network_G files are
def edges_missing(outdir):
	total_edges = 0
	for i in range(30):
		ith_edgelist = []
		intermed = open("{}network_G{}.txt".format(outdir,str(i)), "r")
		for line in intermed:
			intermed_stripped_line = line.strip()
			intermed_line_list = intermed_stripped_line.split()
			ith_edgelist.append(intermed_line_list)
		intermed.close()
		ith_edges = len(ith_edgelist)
		if i == 0:
			total_edges = ith_edges
		edge_difference = total_edges-ith_edges
		edge_percent = round(float(ith_edges)/float(total_edges)*100,2)
		print('Graphlet ' + str(i) + ': ' + str(edge_difference) + ' edges removed or ' + str(edge_percent) + ' percent left.')




# Disease Gene Functions


# preliminary functions


# 'interactome_dg_file' has disease id as the first column and the gene associated with it
# as the last column
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


# 'interactome_mcl_dir' is where the 5 node orca file 'g_network_nodeid.txt' is
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

# 'interactome_mcl_dir' is also where all of the graphlet cluster files are
# this function returns both a list of the clusters from all the graphlets and another
# list with just the lengths of those clusters
def get_all_clusters(interactome_mcl_dir):
	all_clusters = []
	cluster_lengths = []
	
	for i in range(30):

		G_initial = interactome_mcl_dir
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

# this is just useful for making the DataFrame's labels
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


def convert_from_integers(integers,interactome_mcl_dir):
	original = []
	G_initial = interactome_mcl_dir
	Gi = open(G_initial, "r")
	for line in Gi:
		stripped_line = line.strip()
		pair_list = stripped_line.split()

		if pair_list[0] in integers:
			original.append(pair_list[1])
	Gi.close()
	return(original)





# total

# Gets all of the disease gene counts for every cluster of every graphlet.
# interactome_dg_file is the disease gene file
# outdir is where the cluster files are
def get_all_counts(interactome_dg_file,outdir):
	disease_gene_dict = make_disease_gene_dictionary(interactome_dg_file)
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
	return all_counts


# returns an array of the number of disease genes for the max disease clusters for each graphlet
# also returns an array of the original cluster index for each max disease cluster
# also returns an array P
def get_sig_and_index(all_counts, disease_gene_dict):
	P = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	significant_array = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	cluster_index_array = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))
	significant_labels = []
	sig_counts = []

	for j in range(len(disease_gene_dict)):
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
		step1 = significant_array.transpose()
		step2 = step1[~np.all(step1 == 0, axis=1)]
		significant_array = step2.transpose()
	return(significant_array,cluster_index_array,P)

# creates .csv file of hypergeometric scores for each of the 'best' clusters
def hypergeo_dataframe(disease_gene_dict, cluster_index_array, all_clusters,P):
	H = np.zeros(len(disease_gene_dict)*30).reshape(30,len(disease_gene_dict))

	dg_list_sizes = []
	for ids in disease_gene_dict:
		dg_list_sizes.append(len(disease_gene_dict[ids]))

	for j in range(len(disease_gene_dict)):
		for i in range(30):

			num_of_poss_dg = dg_list_sizes[j]
			
			total_genes = 0
			for disease in disease_gene_dict:
				for gene in disease_gene_dict[disease]:
					total_genes+=1
	
			cor_cluster = int(cluster_index_array[i][j]) #corresponding cluster
			chosen_cluster_size = len(all_clusters[i][cor_cluster])

			correct_dg = P[i][j]

			H[i][j] = hypergeom.sf(correct_dg-1, total_genes, num_of_poss_dg, chosen_cluster_size)

	row_names = list_of_strings(list(range(30)),'G')
	row_names.append('DG_size')

	
	ids_lis = []
	for ids in disease_gene_dict:
		ids_lis.append(ids)
	column_names = ids_lis

	pd.DataFrame(H,columns=column_names,index = row_names).to_csv('ppi_dg_hypergeo.csv')





def main(argv):
	clusterdir = argv[1] #'../../out/1_ppi_string/'
	outdir = argv[2] #'../../out/'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	savepath = outdir+clusterdir.split('/')[-2]+'_ari.pdf'
	plottitle = clusterdir.split('/')[-2]

	ari_matrix = get_ari_matrix(clusterdir,3) # 2nd argument is the min. cluster size.
	plot_ari(savepath,ari_matrix,plottitle)

if __name__ == "__main__":
	main(sys.argv)