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
#w = open('%sDDD'%sys.argv[1],'w')
## mean of beta range 
P_values_cutoff=str(sys.argv[2])
## P-value count Cutoff P< 0.05
P_count_cutoff=str(sys.argv[3])
DRECTORY=sys.argv[4]
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

    Family_indx.append(BACTERIAL_ID[-2])
    if 'linear' in BACTERIAL_ID[-1] :
        Family_indx.append(idx[0].split('.OTU.')[1][:idx[0].split('.OTU.')[1].find('.assoc.linear')])
    else :
        Family_indx.append(idx[0].split('.OTU.')[1][:idx[0].split('.OTU.')[1].find('.qassoc')])
    SNP= idx[1]
    P_value= idx[3]
    Beta=idx[2]
    if P_value == 'NA' or Beta =='NA' : 
        continue
    if SNP not in Family.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
    #if SNP not in BETA_IDX.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
        BETA_IDX[SNP]=[float(Beta)]
        Family[SNP] =[Family_indx]
        Family_indx =[]
        Class[SNP]=[BACTERIAL_ID[-2]]
        P_IDX[SNP]=[float(P_value)]
        ### mean of beta analysis
        if len(PRE_SNP) !=0 :
            MEAN_BETA= np.mean(BETA_IDX[PRE_SNP]) 
            cnt_comfim_P = 0
            #print Family[PRE_SNP]
            for j in P_IDX[PRE_SNP] :
                if float(j) < float(P_values_cutoff) : ## P is fixing 
                    cnt_comfim_P = cnt_comfim_P +1  
            if cnt_comfim_P > int(P_count_cutoff):
                print "P-vlause cut off counting:",cnt_comfim_P ,PRE_SNP,len(P_IDX[PRE_SNP])
                fig = plt.figure()
                fig.set_size_inches(10,10) 
                ax = plt.subplot(1, 1, 1)
                plt.ylabel("-log(P-value)",size=15)
                plt.title("Phe-WAS:%s"%PRE_SNP,size= 15)
#               plt.ylim(0,15)
                #P< 1e-4
                plt.axhline(y=4, color='b', linewidth=1,linestyle=':')
                plt.axhline(y=2, color='g', linewidth=1,linestyle=':')
                plt.axhline(y=-math.log10(5e-8), 
                        color='r', linewidth=1,linestyle=':')
                X = 0
                OO = -1
                Cindx = {}
                print PRE_SNP,Class[PRE_SNP][0], Family[PRE_SNP][0]
                CLASS = []
                HI = Class[PRE_SNP]
                for DD in Class[PRE_SNP]:
                    if float(P_IDX[PRE_SNP][X]) == 0 : 
                        Y = 0
                    else:
                        Y= -math.log10(float(P_IDX[PRE_SNP][X]))
                    #Name= Family[PRE_SNP][X][1][3:].split(';')[1]
                    Name =Family[PRE_SNP][X][1].split(';')[-1]
                   # print Name
                    if DD[3:] not in Cindx.keys() : 
                        CLASS.append(DD[3:])
                        OO = OO + 1
                        Cindx[DD[3:]]= OO
                        #print Name
                        if Y > 2 :
                            plt.plot(OO,Y,ls='', marker='o',color= Color[OO%4])
                            plt.text(OO,Y+0.001,'%s'%Name,size= 12)
                        else : 
                            plt.plot(OO,Y,ls='', marker='o',color= Color[OO%4])
                    else: 
                        Cindx[DD[3:]]= Cindx[DD[3:]]+ 0.05
                        Col_idx = Color[int(Cindx[DD[3:]])%4]
                        if Y > 2 :
                            plt.plot(Cindx[DD[3:]],Y,ls='', marker='o',
                                    color= Col_idx)
                            plt.text(Cindx[DD[3:]],Y+0.001,'%s'%Name)
                        else:
                            plt.plot(Cindx[DD[3:]],Y,ls='', 
                                    marker='o',color= Col_idx)
                    X = X +1 
                plt.xticks([GG  for GG in range(len(CLASS))], CLASS)
                for label in ax.xaxis.get_ticklabels() :
                    label.set_rotation(90)
                plt.xlabel("Microbial Class",size=10)
                plt.savefig('%s/%s_Pcount.png'%(DRECTORY,PRE_SNP),
                        dpi = 500,bbox_inches = 'tight')
                plt.close()
                ##################################################################
                ######################## NetWork analysis ########################
                ##################################################################
                #print len(Family)
                #print len(P_IDX)
                name_idx= []
                for kd in Family[PRE_SNP] :
                    name_idx.append(kd[1])
                mydict  = pd.DataFrame({'Family':name_idx,'%s'%PRE_SNP:BETA_IDX[PRE_SNP]})
                if 'Network' not in globals().keys():    
                    Network=mydict
                else :
                    Network = pd.merge(Network, mydict,on ='Family')#),inplace= True)
                print Network.shape
               # print Network
                #print mydict['0']
                # for k in P_IDX[PRE_SNP] :
                #    print k
                 #   display(
                    
                ##################################################################
                ######################## End of Analysis #########################
                ##################################################################
                del BETA_IDX[PRE_SNP]
                del Family[PRE_SNP]
                del Class[PRE_SNP]
                del P_IDX[PRE_SNP]
            else:
                del BETA_IDX[PRE_SNP]
                del Family[PRE_SNP] 
                del Class[PRE_SNP] 
                del P_IDX[PRE_SNP]
    elif SNP in Family.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
    #elif SNP in BETA_IDX.keys() and is_digit(P_value) ==True and is_digit(Beta)==True:
        ## BETA
        BETA_info= BETA_IDX[SNP]
        BETA_info.append(float(Beta))
        BETA_IDX[SNP]= BETA_info
        FFF = Family[SNP]
        FFF.append(Family_indx)
        Family[SNP] = FFF
        Family_indx=[]
        Class_ii = Class[SNP]
        Class_ii.append(BACTERIAL_ID[-2])
        Class[SNP]=Class_ii
        ############P-value
        P_info= P_IDX[SNP]
        P_info.append(float(P_value))
        P_IDX[SNP]= P_info
    cnt = cnt +1
    PRE_SNP=SNP
Network.to_csv("%sfor_network.csv"%sys.argv[1])
print BETA_IDX.keys()
