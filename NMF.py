import sys
import numpy as np
from sklearn.decomposition import NMF
import pandas as pd 
FILE= sys.argv[1]
print FILE
df =pd.read_csv(FILE,delimiter="\t",index_col=0).transpose()

model = NMF(n_components=int(sys.argv[2]), init='random', random_state=12,)
W = model.fit_transform(df)
H = model.components_

k_index=[]
for i in range(int(sys.argv[2])): 
        k_index.append('k%s;k%s'%(str(i+1),str(i+1)))
AA= pd.DataFrame(W,index=df.index,columns=k_index)
BB = pd.DataFrame(H,index=k_index,columns=df.columns)
BB = BB.transpose()
AA.to_csv("%s_W"%FILE,sep='\t')
BB.to_csv("%s_H"%FILE,sep='\t')
