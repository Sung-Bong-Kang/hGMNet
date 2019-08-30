import sys

otu= open('%s'%sys.argv[1])
fam= open('%s'%sys.argv[2]) 
dim_idx_otu=sys.argv[1].split('/')[-1]
print dim_idx_otu
dim_idx_fam='/'.join(sys.argv[2].split('/')[:-1])
print dim_idx_fam
AA = sys.argv[1].find(str(dim_idx_otu))+len(dim_idx_otu)+1
AAA= '%s/%s.pheno'%(dim_idx_fam,dim_idx_otu)

cov = open('%s'%AAA,'w')
        #cov =open('%s/%s.pheno'%(dim_idx_fam,sys.argv[1][AA:sys.argv[1].find('.txt')],'w'))
dic={}
line = otu.readline().strip('\n').split('\t')[1:]
cov.write('FID\tIID\t'+'\t'.join(line)+'\n')
for i in otu :
    D = i.strip('\n')
    idx=D.split('\t')
    dic[idx[0]]=i[len(idx[0]):].strip('\n')#i.find(idx[-2])]

for i in fam : 
    i.strip('\n')
    idx = i.split()[:2]
    if idx[0] in dic.keys(): 
        line= '\t'.join(idx)+dic[idx[0]] +'\n'
        cov.write(line)
