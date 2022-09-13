# Motif/Graphlet enriched clusters. 
# Clustering with respect to a given graphlet (up to size 5)

import numpy as np
import math
import networkx as nx
import pandas as pd
import os
#import markov_clustering as mc
import random
import matplotlib.pyplot as plt
import time

def prepare_network(G_initial,outdir,weighted=True):

	if weighted:
	    G=nx.read_weighted_edgelist(G_initial)
	else:
	    G=nx.read_edgelist(G_initial)


	# write to outdir.
	G = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default', label_attribute='ID')
	pd.DataFrame(list(G.nodes('ID'))).to_csv('{}g_network_nodeid.txt'.format(outdir),sep='\t',header=False,index=False)
	G.remove_edges_from(nx.selfloop_edges(G))

	return G

def final_modules(outdir):
	
	nodeid = pd.read_csv('{}/g_network_nodeid.txt'.format(outdir), sep = '\t', header=None)
	nodeid = dict(zip(nodeid[0],nodeid[1]))
	
	for graphlet in range(0,30):
		with open('{}/module_G{}.txt'.format(outdir,graphlet), 'w') as f1:    
			with open('{}/clusters_G{}.txt'.format(outdir,graphlet)) as f2:
				lines = f2.read().splitlines()
				i=0
				for line in lines:
					if (len(line.split('\t'))>=3): #clusters larger than 2
						f1.write(str(i)+'\t'+'1.0'+'\t')
						for l in line.split('\t'):
							f1.write(nodeid[int(l)]+'\t')
						f1.write('\n')
						i=i+1
						

def perform_clustering(G, outdir, weighted=True):
	# make sure the outdir doesn't contain other large network files.
	nx.write_edgelist(G,outdir+'g_network.txt',data=False,delimiter=" ")
	os.system('sh files_for_orca.sh '+outdir)
	os.system('sh execute_orca.sh edge 5 '+outdir)

	df = pd.read_csv(outdir+'g_network.txt.orca.ocount',sep=" ",header=None)

	edges = pd.read_csv(outdir+'g_network.txt.orca',sep=" ",header=None,skiprows=1)

	edgeorbit = {1: [0], 2: [1], 3: [2,3], 4:[4], 5:[5], 6:[6,7,8], 7:[9,10], 8:[11],
			9: [12,13], 10:[14,15,16], 11:[17], 12:[18,19,20], 13:[21,22,23,24], 14:[25,26,27],15:[28],
			16:[29,30,31], 17:[32,33,34,35], 18:[36,37], 19:[38,39,40,41], 20:[42], 21:[43,44,45,46],
			22:[47,48],23:[49,50,51],24:[52,53,54,55], 25:[56,57,58,59], 26:[60,61,62], 27:[63,64],
			28:[65,66], 29:[67]} # edge orbits counted from G1 (will manually add G0)
	
	n = G.number_of_nodes()
	if weighted:
		A = nx.to_pandas_adjacency(G, weight='weight')
	else:
		A = nx.to_pandas_adjacency(G)
		
	for graphlet in range(0,30): # use 30 for upto five node graphlets.

		P = np.zeros(n*n).reshape(n,n)
		if (graphlet==0):
			for index, i in edges.iterrows():
				P[i[0]][i[1]]=A[i[0]][i[1]]
				P[i[1]][i[0]]=A[i[1]][i[0]]
				if (i[0]==i[1]):
					P[i][j]=0.0

		elif(graphlet > 0):
			for index, i in edges.iterrows():
				if (df.loc[index,edgeorbit[graphlet]].sum()>0):
					P[i[0]][i[1]]=A[i[0]][i[1]]
					P[i[1]][i[0]]=A[i[1]][i[0]]
				if (i[0]==i[1]):
					P[i][j]=0.0

		step1 = nx.to_networkx_graph(P)
		nx.info(step1)
		step2 = nx.to_edgelist(step1)
		print(len(step2))

		
		f = open("{}network_G{}.txt".format(outdir,graphlet), "w")
		for line in step2:
			#f.write(str(line[0])+"	"+str(line[1])) #unweighted
			f.write(str(line[0])+" "+str(line[1])+" "+str(line[2]['weight'])) #weighted
			f.write("\n")
		f.close()

		# run mcl using the os module (like a command line)
		# I is the inflation parameter - set to 4
		cmd = "mcl {}network_G{}.txt --abc -I 4 -o {}clusters_G{}.txt".format(outdir,graphlet,outdir,graphlet)
		os.system(cmd)
		




#specify the output directory.
outdir = '../../out/1_ppi_string/'

#interactome 
G_initial = '../../data/DREAM/1_networks/original/1_ppi_string_cutoff_8.txt'

weighted = True

G = prepare_network(G_initial,outdir,weighted)

perform_clustering(G, outdir, weighted)

final_modules(outdir)