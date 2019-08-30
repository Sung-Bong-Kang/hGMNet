import sys
otu = open('%s'%sys.argv[1])
idx=otu.readline().strip('\n').split('\t')
idx2 = idx[1:]
cnt = -1
## kingdom:1 ,phylum:2, class:3, order:4, family:5,genus:6, species:7 
Class=int(sys.argv[2])
print len(idx2)
##same_column= []
same_column={}
save= []
#order = [] 
for i in idx2 : 
    cnt = cnt +1
    idx = i.split(';')
    if ';'.join(idx[:Class]) in same_column : 
        I=same_column[';'.join(idx[:Class])] 
        I.append(cnt)
        same_column[';'.join(idx[:Class])]= I
      #  save = []
    else :
        save.append(cnt)
        same_column[';'.join(idx[:Class])] = save
#        print save,same_column
     #   order.append(same_column.values())
        save =[]
#print same_column,'\n',len(same_column),'\n', same_column.values() 
#print same_column
VALUES= same_column.values()
cnt0 = 0
fail = [] 
fail_name = []
names = same_column.keys()
if sys.argv[2] != '7' :
    for i in same_column.keys() : 
        idx = i.split(';')
        if i[-2:] =='__' or '__;' in i : 
            fail.append(cnt0+1)
            fail_name.append(i)
        cnt0 =cnt0 +1 
    for i in fail_name :
        names.remove(i)
if '.txt' in sys.argv[1] : 
    file_idx= sys.argv[1].find('.txt')
else : 
    file_idx = len(sys.argv[1])
NEW_file=open('%s_cluster_%s_level.txt'%(sys.argv[1][:file_idx],sys.argv[2]),'w') 
Header = 'SAMPLES\t'+'\t'.join(names)+'\n'
NEW_file.write(Header)
array= []
#print same_column.keys()
import numpy as np 
import pandas as pd
sample=[]
DIGO = []
for i in otu :
    idx = i.strip('\n').split()[1:] 
    sample.append(i.split()[0])
    array.append(idx)
##    DIGO.append(i.split())
data =pd.DataFrame(np.array(array,float))
DD = 0
DF = []

for i in VALUES :
    DF.append(data[i].sum(axis=1))
data2 = pd.DataFrame(np.array(DF,str))
cnt = 0
array = [] 
data = []
fail.reverse()
for i in sample :
    line =i+'\t'+'\t'.join(data2[cnt])+'\n'#'\t'+'\t'.join(DIGO[cnt])+'\n'
    idx = line.split()
    for j in fail :
        del idx[j]
    NEW_file.write('\t'.join(idx)+'\n')
    cnt = cnt+1 

