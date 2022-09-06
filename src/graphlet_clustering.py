# Motif/Graphlet enriched clusters. 
# Clustering with respect to a given graphlet (up to size 5)

import numpy as np
import math
import networkx as nx
import pandas as pd
import os
import random
import matplotlib.pyplot as plt
import time
import os.path


def perform_clustering(G, outdir_beginning, orca_direct, orca_num):
	# make sure the outdir doesn't contain other large network files.
	nx.write_edgelist(G,outdir_beginning+orca_direct+'g_network.txt',data=False,delimiter=" ")
	os.system('sh files_for_orca.sh '+outdir_beginning+orca_direct)
	os.system('sh execute_orca.sh edge '+orca_num+' '+outdir_beginning+orca_direct)

	print("1) reading csv files")
	df = pd.read_csv(outdir_beginning+orca_direct+'g_network.txt.orca.ocount',sep=" ",header=None)

	edges = pd.read_csv(outdir_beginning+orca_direct+'g_network.txt.orca',sep=" ",header=None,skiprows=1)

	edgeorbit = {1: [0], 2: [1], 3: [2,3], 4:[4], 5:[5], 6:[6,7,8], 7:[9,10], 8:[11],
			9: [12,13], 10:[14,15,16], 11:[17], 12:[18,19,20], 13:[21,22,23,24], 14:[25,26,27],15:[28],
			16:[29,30,31], 17:[32,33,34,35], 18:[36,37], 19:[38,39,40,41], 20:[42], 21:[43,44,45,46],
			22:[47,48],23:[49,50,51],24:[52,53,54,55], 25:[56,57,58,59], 26:[60,61,62], 27:[63,64],
			28:[65,66], 29:[67]} # edge orbits counted from G1 (will manually add G0)
	
	n = G.number_of_nodes()
	A = nx.to_pandas_adjacency(G)

	

	for graphlet in graphlet_num: # use 30 for up to five node graphlets.


		# makes a directory for the specific graphlet, ie '../huri_out/G2'
		outdir = outdir_beginning + 'G' + str(graphlet) + '/'
		if os.path.isdir(outdir) == False:
			directory = 'G' + str(graphlet)
			path = os.path.join(outdir_beginning, directory)
			os.makedirs(path) 
			print("Directory '% s' created" % directory) 


		print("2) Making probability matrix for G"+str(graphlet))

		P = np.zeros(n*n).reshape(n,n)
		if (graphlet==0):
			for index, i in edges.iterrows():
				P[i[0]][i[1]]=A[i[0]][i[1]]
				P[i[1]][i[0]]=A[i[1]][i[0]]
				if (i[0]==i[1]):
					P[i][j]=0.0

		elif(graphlet > 0):
			if type == 'unweighted':
				print("made it here")
				for index, i in edges.iterrows():
					if (df.loc[index,edgeorbit[graphlet]].sum()>0):
						P[i[0]][i[1]] = A[i[0]][i[1]]
						P[i[1]][i[0]] = A[i[1]][i[0]]
					if (i[0]==i[1]):
						
						P[i][j]=0.0
			elif type == 'weighted':
				for index, i in edges.iterrows():
					orbit_total = df.loc[index,edgeorbit[graphlet]].sum()
					if (orbit_total>0):
						P[i[0]][i[1]] = orbit_total
						P[i[1]][i[0]] = orbit_total
					if (i[0]==i[1]):
						P[i][j]=0.0

		

		step1 = nx.to_networkx_graph(P)

		step2 = nx.to_edgelist(step1)

		
		
		f = open(outdir+"intermed_G"+str(graphlet)+".txt", "w")
		for line in step2:
			weight_text = str(line[2])
			weight_list = weight_text.split()
			last_weight = weight_list[1]
			weight = last_weight[:-1]
			f.write(str(line[0])+"	"+str(line[1])+"	"+str(weight))
			f.write("\n")
		f.close()

		print("3) running mcl & making "+outdir_beginning + "inflation_tests/G" +str(graphlet)+ "_I3.25")
		# run mcl using the os module (like a command line)
		os.system("mcl " +outdir+ "intermed_G" +str(graphlet)+ ".txt -I 3.25 --abc -o " +outdir_beginning + "inflation_tests/G" +str(graphlet)+ "_I3.25")
		
		
#specify the output directory.
interact = input('Interactome: ')
orca_num = input('Orca: ')


graphlet_num = [28,29]
#[7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
#[int(input('Graphlet: '))]


# Graphs must be converted to simple txt files of the format 
# node1		node2
if interact == 'apc':
	G_initial = '../data/interactomes/All_Pathway_Commons.txt'
	outdir_beginning = '../apc_out/'
	type = 'unweighted'
if interact == 'weighted apc':
	G_initial = '../data/interactomes/All_Pathway_Commons.txt'
	outdir_beginning = '../weighted_apc_out/'
	type = 'weighted'
if interact == 'huri':
	G_initial = '../data/interactomes/huri.tsv'
	outdir_beginning = '../huri_out/'
	type = 'unweighted'
if interact == 'weighted huri':
	G_initial = '../data/interactomes/huri.tsv'
	outdir_beginning = '../weighted_huri_out/'
	type = 'weighted'
if interact == 'random six':
	G_initial = '../data/interactomes/six_graph.txt'
	outdir_beginning = '../random_six_out/'
	type = 'unweighted'
if interact == 'random one':
	G_initial = '../data/interactomes/one_graph.txt'
	outdir_beginning = '../random_one_out/'
	type = 'unweighted'
if interact == 'ppi':
	G_initial = '../data/interactomes/PP-Pathways_ppi.csv'
	outdir_beginning = '../ppi_out/'
	type = 'unweighted'


if orca_num == '4':
	orca_direct = 'orca4/'
if orca_num == '5':
	orca_direct = 'orca5/'
if os.path.isdir(outdir_beginning+orca_direct) == False:
		directory = orca_direct
		path = os.path.join(outdir_beginning, directory)
		os.makedirs(path) 
		print("Directory '% s' created" % directory)

if interact != 'ppi':
	list_of_lists = []
	Gi = open(G_initial, "r")
	for line in Gi:
		stripped_line = line.strip()
		line_list = stripped_line.split()
		list_of_lists.append(line_list)
	Gi.close()
	G = nx.from_edgelist(list_of_lists)
if interact == 'ppi':
	list_of_lists = []
	Gi = open(G_initial, "r")
	for line in Gi:
		stripped_line = line.strip()
		line_list = stripped_line.split(',')
		list_of_lists.append(line_list)
	Gi.close()
	G = nx.from_edgelist(list_of_lists)



# if node ids are not 0-n.
G = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default', label_attribute='ID')
pd.DataFrame(list(G.nodes('ID'))).to_csv(outdir_beginning+orca_direct+'g_network_nodeid.txt',sep='\t',header=False,index=False)
G.remove_edges_from(nx.selfloop_edges(G))


#perform_clustering(G, outdir_beginning, orca_direct, orca_num)




nx.write_edgelist(G,outdir_beginning+orca_direct+'g_network.txt',data=False,delimiter=" ")
os.system('sh files_for_orca.sh '+outdir_beginning+orca_direct)
os.system('sh execute_orca.sh edge '+orca_num+' '+outdir_beginning+orca_direct)

