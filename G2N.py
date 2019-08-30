import sys

map_gene = open('%s'%sys.argv[1])
net=open('%s'%sys.argv[2])
w = open('%s_GENE.csv'%sys.argv[2][:sys.argv[2].find('.csv')],'w')
cnt = 0 
dic = {}
for i in map_gene:
    idx = i.strip('\n').split()
    dic[idx[2]]= idx[0]
head =net.readline().strip('\n').split(',')
save =[]
new_head=',Family'
for i in head : 
    if i in dic.keys() : 
        save.append(cnt)
        new_head = new_head+',%s'%dic[i]
    cnt = cnt + 1 
new_head = new_head+'\n'
w.write(new_head)
for i in net : 
    idx = i.strip('\n').split(',')
    line = '%s,%s'%(idx[0],idx[1]) 
    for j in save : 
        line= line +',%s'%idx[j] 
    w.write(line+'\n')

