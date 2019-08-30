# Automatic pipeline of microbiome PheWAS for the discovery of host-microbe interaction networks -BETA Versions
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

# How to use ?

### Analysis Mode 1 : make PheWAS figures, no covariate , analysis bacterial-Class:Family
  /downloaded/hGMNet/Path/hGMNet.sh \
    --OTU_ID your_OTU.txt \
    --Bacterial_class F # choise S,G,F,O,C,P,K \
    --OTU_DIR /your/OTU/path \
    --Input_prefix plink_file_id \
    --DIR /your/plink/.bed.bim.fam/path \
    --Analysis Linear \
    --P_cut 5e-6 \
    --P_count 10 \
    --PHEWAS_image_mode Y
        
### Analysis Mode 2 : non make PheWAS figures, Covariate, analysis bacterial-Class:Species

/downloaded/hGMNet/Path/hGMNet.sh \
    --OTU_ID your_OTU.txt \ 
    --Bacterial_class S \
    --OTU_DIR /your/OTU/path \
    --Input_prefix plink_file_id \
    --DIR /your/plink/.bed.bim.fam/path \
    --Analysis Linear \
    --P_cut 5e-6 \
    --P_count 4 \
    --Cov /covariate/path/covariate.txt \
    --Cov_names age,sex,bmi,.etc
### Analysis Mode 3 : Non-negative matrix factorization(NMF) analysis , non PheWAS immage 
  /downloaded/hGMNet/Path/hGMNet.sh \
    --OTU_ID your_OTU.txt \
    --Bacterial_class F # choise S,G,F,O,C,P,K \
    --OTU_DIR /your/OTU/path \
    --Input_prefix plink_file_id \
    --DIR /your/plink/.bed.bim.fam/path \
    --Analysis NMF \
    --NMF_K 8 \
    --P_cut 5e-6 \
    --P_count 10 
    
    
    
# Reference 
 

