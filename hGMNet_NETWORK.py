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
from sklearn.metrics.pairwise import cosine_similarity
### PHEWAS summary file 
df = pd.read_csv('%s'%sys.argv[1])
df.index =df['Family']
del df['Family']
del df['Unnamed: 0']
###########
Family=[]
name_family =[] 
for i in df.index :
    if sys.argv[3] in i :
        if i[i.find(sys.argv[3]):] != sys.argv[3]:
            Family.append(i)
            if i[i.find('%s'%sys.argv[3]):] in name_family :
                name_family.append(i[i.find(i.split(';')[-2]):])
            else:
                name_family.append(i[i.find('%s'%sys.argv[3]):])
df.index = name_family


fig = plt.figure(figsize=(20,20))
sns.clustermap(df,cmap='RdYlGn_r',linewidths=.7)
plt.title('colormap of cmap=RdYlGn_r', fontsize=20)
plt.savefig('%s.SNP_BACTERIAL_BETA_heat_map.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
fig = plt.figure(figsize=(20,20))
G = ''
corr =df.corr('pearson')
links = corr.stack().reset_index()
links.columns = ['var1', 'var2','value']
links_filtered=links.loc[ (links['value'] > 0.8) & (links['var1'] != links['var2']) ]
plt.figure(figsize=(20,20))
sns.clustermap(data = df.corr(),
        fmt = '.2f', linewidths=.5, cmap='RdYlGn_r',vmin=0 , vmax=1 )
plt.savefig('%s.SNP_SNP_pearson_heatmap.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
fig = plt.figure(figsize=(20,20))
G=nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
plt.figure(figsize=(8,8))
nx.draw(G, with_labels=True, node_color='orange', node_size=400 ,edge_color='black', linewidths=1, font_size=12)
plt.savefig('%s.SNP_SNP_pearson_network.png'%sys.argv[1],dpi =300,bbox_inches = 'tight')
plt.close()

G2= ''
df2 = pd.read_table('%s'%sys.argv[2])
df2.index=df2["SAMPLES"]
del df2['SAMPLES']
#del df2['IID']
cnt = 0 
Family =  []
name_family = [] 
for i in df2.columns : 
    if sys.argv[3] in i :
        if i[i.find(sys.argv[3]):] != sys.argv[3]:
            Family.append(i)
            if i[i.find('%s'%sys.argv[3]):] in name_family :
                name_family.append(i[i.find(i.split(';')[-2]):])
            else: 
                name_family.append(i[i.find('%s'%sys.argv[3]):])        
        cnt = cnt +1 
df2 = df2[Family]
COS=pd.DataFrame(cosine_similarity(df2.transpose()))
#COS = COS[Family]
COS.index=list(name_family)
COS.columns=list(name_family)
COS.index.name='Family'
plt.figure(figsize=(20,20))
sns.clustermap((COS), cmap='RdYlGn_r')
plt.savefig('%s_BACTERIA_BACTERIA_clusterMap.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()
AA = COS.stack()
cnt =0
Nal = []
for i in  AA.index :
    Nal.append([i[0],i[1],AA[i]])
links2= pd.DataFrame(Nal)
Nal = []
fig = plt.figure(figsize=(20,20))
links2.columns = ['var1', 'var2','value']
if sys.argv[3] == 'f__' or sys.argv[3] == 'o__': 
    percut= 95
elif sys.argv[3] == 'g__' or sys.argv[3] == 's__':
    percut= 99
else : 
    percut= 90
links_filtered2=links2.loc[ (links2['value'] > np.percentile(links2['value'] ,percut)) & (links2['var1'] != links2['var2']) ]

G2=nx.from_pandas_edgelist(links_filtered2, 'var1', 'var2')
nx.draw(G2, with_labels=True, node_color='cyan', node_size=400 ,edge_color='black', linewidths=1, font_size=12,alpha= 0.8)
plt.savefig('%s_bac_bac_network.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()

BETA_col = list(df.columns)
BETA_IDX = list(df.index)
## Beta cutoff ? 0.03 ? 
positive_edge = [] 
negative_edge= [] 
NODES = [] 
bacterial_cnt =0
for i in np.array(df) :
    snp_cnt = 0
        #    print i[COUNT]
    for j in i : 
        EDGE=[BETA_col[snp_cnt],BETA_IDX[bacterial_cnt],j]
        if j >  np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) and BETA_col[snp_cnt] in G.nodes and BETA_IDX[bacterial_cnt] in G2.nodes :
            positive_edge.append(EDGE)
        elif j < -np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) and BETA_col[snp_cnt] in G.nodes and BETA_IDX[bacterial_cnt] in G2.nodes:
            negative_edge.append(EDGE)
        elif  j >  np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) :
            positive_edge.append(EDGE)
            NODES.append(BETA_IDX[bacterial_cnt])
            G.add_node(list(df.columns)[snp_cnt])
        elif  j< -np.percentile(df[list(df.columns)[snp_cnt]].abs(),95) :
            negative_edge.append(EDGE)
            NODES.append(BETA_IDX[bacterial_cnt])
            G.add_node(list(df.columns)[snp_cnt])
        snp_cnt =snp_cnt +1 
    bacterial_cnt = bacterial_cnt +1    
labels = {}

G3 = nx.from_pandas_edgelist(links_filtered2, 'var1', 'var2')
#G3 = nx.Graph()
#G3.add_nodes_from(G2.nodes)
#G3.add_edges_from(G2.edges)
G3.add_nodes_from(G.nodes)
G3.add_edges_from(G.edges)
G3.add_weighted_edges_from(positive_edge)
G3.add_weighted_edges_from(negative_edge)
G3.add_nodes_from(NODES)
#G3.add_edges_from(positive_edge)
#G3.add_edges_from(negative_edge)
pos= nx.spring_layout(G3,k=0.3,threshold= 1e-8) 
Faile = []
for i in G.nodes: 
    p=nx.shortest_path_length(G3,source=i)
    for j in G3.nodes : 
        try:
            p[j]
        except KeyError:
            Faile.append(j)
G3.remove_nodes_from(Faile)
G2.remove_nodes_from(Faile)
for i in list(G3.nodes): 
        labels[i] = i
fig = plt.figure()
fig.set_size_inches(20,15)
ax = plt.subplot(1, 1, 1)
nx.draw_networkx_nodes(G3, pos,nodelist=list(G2.nodes)+NODES,node_color='c',node_size=800,alpha=1)
nx.draw_networkx_nodes(G3, pos , nodelist=G.nodes, node_color='orange',node_size=800,alpha=1)
nx.draw_networkx_edges(G3, pos,edgelist=list(G2.edges)+list(G.edges))
#nx.draw_networkx_edges(G3,pos,edgelist =G.edges)
nx.draw_networkx_edges(G3,pos,edgelist=positive_edge,edge_color='b',width=[1,3,5],alpha=0.4)
nx.draw_networkx_edges(G3,pos,edgelist=negative_edge,edge_color='r',width=[1,3,5],alpha=0.4)
nx.draw_networkx_labels(G3, pos,labels=labels,font_size=15)
if len(list(G2.edges)+list(G.edges)) != 0 and len(positive_edge) == 0: 
    ax.legend(['Bacteria','SNP','Correlation','negative interaction'],loc='best',fontsize =('xx-large'))
elif len(list(G2.edges)+list(G.edges)) != 0 and len(positive_edge) != 0:
    ax.legend(['Bacteria','SNP','Correlation','positive correation','negative interaction'],loc='best',fontsize =('xx-large'))
elif len(list(G2.edges)+list(G.edges)) == 0 and len(positive_edge) == 0  : 
    ax.legend(['Bacteria','SNP','negative interaction'],loc='best',fontsize =('xx-large'))
else : 
    ax.legend(['Bacteria','SNP','positive correation','negative interaction'],loc='best',fontsize =('xx-large'))

#ax.legend(['Bacteria','SNP','Correlation','positive interaction','negative interaction'],loc='best',fontsize =('xx-large'))
plt.tight_layout()
plt.ylim(-1.1,1.1)
plt.xlim(-1.4,1.4)
plt.savefig('%s_snp_bacteria_network.png'%sys.argv[1],dpi = 300,bbox_inches = 'tight')
plt.close()

### SAVE SIG SNP and Bacteria in PHEWAS results
nod_idx=list(G3.nodes())
in_bac=[]
for i in nod_idx:
    if i in list(df.index) :
        in_bac.append(i)
ADS= df.loc[in_bac]
ADS.index.name = 'Bactera_%s'%sys.argv[3][0]
ADS.to_csv('%s_Networked_SNP_Bac.csv'%sys.argv[1][:sys.argv[1].find('.csv')],sep='\t')

