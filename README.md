# Automatic pipeline of microbiome PheWAS for the discovery of host-microbe interaction networks
## hGMNet : host Genetics and Microbe interaction Networks 
  
 
![image](https://user-images.githubusercontent.com/51352117/64014502-3b5f2900-cb5d-11e9-9ad0-e333ade4e442.png)

 Gene mode is not available 
 
# Required Libraries and Tools

`Python 2.7` version from https://www.python.org/download/releases/2.7

###python libraries : 

 `matplotlib`  :  `pip install matplotlib`
 
 `math`
 
 `pandas 0.23.4`  :  `pip install pandas` 
 
 `numpy 1.15.4`  :  `pip install numpy`
 
 `seaborn 0.9.0` :  `pip install seaborn`
 
 `networkx 2.2`  :  `pip install networkx`
 
 `sklearn` from  

`Plink 1.9` version for genotype analysis from `https://www.cog-genomics.org/plink2`


`R` version 3.4 or higher 

# How to use

### Analysis Mode 1 : make PheWAS figures, no covariate , analysis Class Family
  /downloaded/hGMNet/Path/hGMNet.sh \
        --OTU_ID your_OTU.txt \
        --Bacterial_class F \
        --OTU_DIR /data/rowi007/SCRIPT/PHEWAS/EXAMPLE \
        --Input_prefix KB \
        --DIR /data/rowi007/SCRIPT/PHEWAS/EXAMPLE \
        --Analysis Linear \
        --P_cut 5e-6 \
        --P_count 10 \
        --Gene_mode N \
        --PHEWAS_image_mode Y
        
### Analysis Mode 2 : non make PheWAS figures, 

### Analysis Mode 3 : Non-negative matrix factorization(NMF) analysis 
