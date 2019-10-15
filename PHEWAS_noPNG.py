import sys 
import os
import matplotlib as mpl
#if os.environ.get('DISPLAY','') == '':
 #   print('no display found. Using non-interactive Agg backend')
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
## check string 
def is_digit(str):
    try:
        tmp=float(str)
        return True
    except ValueError:
        return False
Color = {}
## b g r c m y k 
for i in range(400):
    if i%4 ==0 :
        Color[i] = 'b'
    elif i%4 ==1 :
        Color[i] = 'r'
    elif i%4 ==2 : 
        Color[i] = 'y'
    elif i%4 ==3 :
        Color[i] = 'm'
## merged : merged and sorted bacterial Qassoc  
merged=open('%s'%sys.argv[1],'r')
w = open('%s.PASSed.SNP'%sys.argv[1],'w')
## mean of beta range 
P_value_cutoff=str(sys.argv[2])
## P-value count Cutoff P< 0.05
P_count_cutoff=str(sys.argv[3])
DIREACTORU= sys.argv[4]
## column  ID ,SNP , BETA, P
ID_IDX= [] ## ID_IDX : DIR/BACTERIAL full ID 
SNP_ID= []
BETA_IDX = {}
P_IDX = {}
cnt =0
PRE_SNP=''
Class = {}
Family= {}
Family_indx = []
for i in merged : 
    idx = i.strip('\n').split('\t')
    BACTERIAL_ID=idx[0].split('.OTU.')[1].split(';')
    ######################## fix
    Family_indx.append(BACTERIAL_ID[-2])
    if 'linear' in BACTERIAL_ID[-1] : 
        Family_indx.append(idx[0].split('.OTU.')[1][:idx[0].split('.OTU.')[1].find('.assoc.linear')])
    else :
        Family_indx.append(idx[0].split('.OTU.')[1][:idx[0].split('.OTU.')[1].find('.qassoc')])
    #print Family_indx,BACTERIAL_ID
    #Family_indx.append('%s;%s'%(BACTERIAL_ID[-2],inside))
    #Family_indx.append(BACTERIAL_ID[-1][:BACTERIAL_ID[-1].find('.qassoc')])
    SNP= idx[1]
    P_value= idx[3]
    Beta=idx[2]
    if P_value == 'NA' or Beta =='NA' : 
        continue
    #if SNP not in BETA_IDX.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
    if SNP not in Family.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
        BETA_IDX[SNP]=[float(Beta)]

        Family[SNP] =[Family_indx]
        Family_indx =[]
        ###################### fix
        Class[SNP]=[BACTERIAL_ID[-2]]
        P_IDX[SNP]=[float(P_value)]
        if len(PRE_SNP) !=0 :
            MEAN_BETA= np.mean(BETA_IDX[PRE_SNP]) 
            cnt_comfim_P = 0
            for j in P_IDX[PRE_SNP] :
                if float(j) < float(P_value_cutoff) : ## P is fixing 
                    cnt_comfim_P = cnt_comfim_P +1  
            if cnt_comfim_P > int(P_count_cutoff):
                w.write('%s\n'%PRE_SNP)
                print "P-vlause cut off counting:",cnt_comfim_P ,PRE_SNP,len(P_IDX[PRE_SNP]),len(Family[PRE_SNP]),len(Class[PRE_SNP])
                OO = -1
                X = 0
                Cindx = {}
                CLASS = []
                for DD in Class[PRE_SNP] :
                    if float(P_IDX[PRE_SNP][X]) == 0 : 
                        Y = 1
                    else:
                        Y= -math.log10(float(P_IDX[PRE_SNP][X]))
                    #Name= Family[PRE_SNP][X][1][3:].split(';')[1]
                    Name =Family[PRE_SNP][X][1].split(';')[-1]
                    if DD[3:] not in Cindx.keys() : 
                        CLASS.append(DD[3:])
                        OO = OO + 1
                        Cindx[DD[3:]]= OO
                    else: 
                        Cindx[DD[3:]]= Cindx[DD[3:]]+ 0.05
                        Col_idx = Color[int(Cindx[DD[3:]])%4]
                    X = X +1 
                name_idx= []
                for kd in Family[PRE_SNP] :
                    name_idx.append(kd[1])
                mydict  = pd.DataFrame({'Family':name_idx,'%s'%PRE_SNP:BETA_IDX[PRE_SNP]})
                if 'Network' not in globals().keys():
                    Network=mydict
                else :
                    Network = pd.merge(Network, mydict,on ='Family')#),inplace= True)
                print Network.shape
                del Family[PRE_SNP]
                del Class[PRE_SNP]
                del P_IDX[PRE_SNP]
            else:
                del BETA_IDX[PRE_SNP]
                del Family[PRE_SNP] 
                del Class[PRE_SNP] 
                del P_IDX[PRE_SNP]
    #elif SNP in BETA_IDX.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
    elif SNP in Family.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
    ## BETA
        BETA_info= BETA_IDX[SNP]
        BETA_info.append(float(Beta))
        BETA_IDX[SNP]= BETA_info
        #print Family[SNP]
        FFF = Family[SNP]
        FFF.append(Family_indx)
        Family[SNP] = FFF
        Family_indx=[]
        Class_ii = Class[SNP]
        ########## fix############
        Class_ii.append(BACTERIAL_ID[-2])
        Class[SNP]=Class_ii
        ############P-value
        P_info= P_IDX[SNP]
        P_info.append(float(P_value))
        P_IDX[SNP]= P_info
    cnt = cnt +1
    PRE_SNP=SNP
Network.to_csv("%sfor_network.csv"%(sys.argv[1])) 
