import pandas as pd
import numpy as np
import os
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

def plot_ari(savepath,ari_matrix):
	plt.clf()
	plt.matshow(ari_matrix)
	plt.colorbar()
	plt.title('1_ppi_string (weights > 0.8)')
	plt.savefig(savepath)


clusterdir = '../../../../GitHub/graphlet-clustering/out/1_ppi_string/'
outdir = '../../../../GitHub/graphlet-clustering/out/'
savepath = outdir+'1_ppi_string'+'_ari.pdf'

ari_matrix = get_ari_matrix(clusterdir,3)
plot_ari(savepath,ari_matrix)
