import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
from IPython.core.display import Image
from networkx.drawing.nx_pydot import to_pydot
from itertools import combinations


df = pd.read_csv('%s'%sys.argv[1])
df.index=df['Family']
del df['Family']
del df['Unnamed: 0']
#################################################
################## SNP NMF H heat map ############
fig = plt.figure(figsize=(12,12))
sns.clustermap(df,cmap='RdYlGn_r',linewidths=.7)
plt.title('NMF H and SNP Cluster', fontsize=20)
plt.savefig('%s.SNP_BACTERIAL_BETA_heat_map.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
#################################
#################### SNP correlation heat map and  SNP interaction network 

plt.figure(figsize=(12,12))
corr =df.corr('pearson')
sns.clustermap(data = df.corr(),
        fmt = '.2f', linewidths=.5, cmap='RdYlGn_r')
plt.title('SNP and SNP Beta correlation',fontsize=20)
plt.savefig('%s.SNP_SNP_correaltion_map.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
links = corr.stack().reset_index()
links.columns = ['var1', 'var2','value']
links_filtered=links.loc[ (links['value'] > np.percentile(links['value'],75)) & (links['var1'] != links['var2']) ]
G=nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
plt.figure(figsize=(12,12))
nx.draw(G, with_labels=True, node_color='orange', node_size=400 ,edge_color='black', linewidths=1, font_size=12)
plt.savefig('%s.SNP_SNP_correaltion_network.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()


NMF_W = pd.read_csv(sys.argv[2],sep='\t',index_col=0)
Family=[]
name_family =[] 
cnt = 0
for i in NMF_W.index : 
    if sys.argv[3] in i :
        if i[i.find(sys.argv[3]):] != sys.argv[3]:
            Family.append(i)
            if i[i.find(sys.argv[3]):] in  name_family :
                name_family.append(i[i.find(i.split(';')[-2]):])
            else :
                name_family.append(i[i.find(sys.argv[3]):])
        cnt = cnt +1
NMF_W =NMF_W.transpose()[Family].transpose()#()[T_Names]
NMF_W.index= name_family
fig = plt.figure(figsize=(12,12))
sns.clustermap(NMF_W,cmap='RdYlGn_r')
plt.title('NMF W and bacteria correlation map', fontsize=20)
plt.savefig('%s.W_Bacterial_correlation.map.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
W_m = nx.Graph()
### NMF K and K interation 
corr_W = NMF_W.corr()
links3 = corr_W.stack().reset_index()
links3.columns = ['var1', 'var2','value']
links_filtered3=links3.loc[ (links3['value'] >  np.mean(links3['value'])) & (links['var1'] != links['var2']) ]
W_m=nx.from_pandas_edgelist(links_filtered3, 'var1', 'var2')

############ NMF k and SNP interation ###########
BETA_col = list(df.columns)
BETA_IDX = list(df.index)
positive_edge = [] 
negative_edge= [] 
bacterial_cnt =0
NODES = [] 
for i in np.array(df) :
    snp_cnt = 0
    for j in i :
        EDGE=[BETA_col[snp_cnt],BETA_IDX[bacterial_cnt],j]
        if j > np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) and BETA_col[snp_cnt]in G.nodes :
            positive_edge.append(EDGE)
            NODES.append(BETA_IDX[bacterial_cnt])
            
        elif  j< -np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) and BETA_col[snp_cnt] in G.nodes :
            negative_edge.append(EDGE)
            NODES.append(BETA_IDX[bacterial_cnt])
        snp_cnt =snp_cnt +1
    bacterial_cnt = bacterial_cnt +1

labels = {}
for i in list(df.index):
        labels[i]=i.split(';')[1]
G2 = nx.Graph()
G2.add_nodes_from(df.index)
G2.add_edges_from(G.edges)
G2.add_weighted_edges_from(positive_edge)
G2.add_weighted_edges_from(negative_edge)
G2.add_edges_from(W_m.edges)

pos= nx.spring_layout(G2,k=0.3,threshold= 1e-8) 
fig = plt.figure()
ax = plt.subplot(1, 1, 1)
fig.set_size_inches(20,15)
nx.draw_networkx_nodes(G2,pos,nodelist=list(df.index),node_color='c')
nx.draw_networkx_nodes(G2,pos,nodelist=list(G.nodes),node_color='orange')
nx.draw_networkx_edges(G2,pos,edgelist=positive_edge,edge_color='b',width=[1,3,5],alpha=0.4)
nx.draw_networkx_edges(G2,pos,edgelist=negative_edge,edge_color='r',width=[1,3,5],alpha=0.4)
nx.draw_networkx_edges(G2,pos,edgelist=list(G.edges)+list(W_m.edges))
nx.draw_networkx_labels(G2, pos,labels=labels,font_size=17,)
plt.title('SNP and NMF cluster Network',fontsize= 20)
if len(positive_edge) == 0 : 
    ax.legend(['Bacteria','SNP','negative interaction','Correlation'],loc='best',fontsize =('xx-large'))
elif len(negative_edge) == 0:
    ax.legend(['Bacteria','SNP','positive correation','Correlation'],loc='best',fontsize =('xx-large'))
else : 
     ax.legend(['Bacteria','SNP','positive correation','negative interaction','Correlation'],loc='best',fontsize =('xx-large'))
plt.savefig('%s_SNP_bacteria_network.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
############# SAVE NETWORK ##########
DO = []
DA = [] 
df.index.name ='Bacterial_cluster'
for i in list(G2.node) : 
    if i in df.index : 
        DO.append(i)
    else : 
        DA.append(i)
df.loc[DO][DA].to_csv('%s_SIG_NODE_Beta.csv'%sys.argv[1])
NMF_W.index.name ='Bacterial_cluster'
NMF_W[NMF_W>np.percentile(NMF_W,90)].to_csv('%s.NMF_90percentile.bacteria.csv'%sys.argv[2])
pd.DataFrame(list(G2.edges),columns=['Node1','Node2']).to_csv('%s.EDGEs.csv'%sys.argv[1])
