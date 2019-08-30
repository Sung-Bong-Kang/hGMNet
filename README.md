# Automatic pipeline of microbiome PheWAS for the discovery of host-microbe interaction networks -BETA Versions
## hGMNet : host Genetics and Microbe interaction Networks 
 
 The assay, which proceeds to measure changes in microbe abundance according to a single genotype, does not take into account the interactions of bacteria and genetics.
 
 Therefore, we have created a tool to discover new host genetics and microbe interaction networks using microbiome PheWAS.
 
![image](https://user-images.githubusercontent.com/51352117/64018100-6c902700-cb66-11e9-8d2f-c99823794e8e.png)

 
## Required Libraries and Tools

`Python 2.7` version from https://www.python.org/download/releases/2.7

### python libraries : 

 `matplotlib`  :  `pip install matplotlib`
 

 
 `pandas 0.23.4`  : pip install pandas
 
 `numpy 1.15.4`  :  pip install numpy
 
 `seaborn 0.9.0` :  pip install seaborn
 
 `networkx 2.2`  :  pip install networkx
 
 `sklearn` , `math`

`Plink 1.9` version for genotype analysis from https://www.cog-genomics.org/plink2


`R` version 3.4 or higher 

`db19_20k.gz` for Gene mode from  https://drive.google.com/open?id=1hEUdViceUQIO-_-zSShxUqW6W4qashXu

### How to use ?
옵션설명

#### Analysis Mode 1 : make PheWAS figures, no covariate , analysis bacterial-Class:Family
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
        
#### Analysis Mode 2 : non make PheWAS figures, Covariate, analysis bacterial-Class:Species
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
    
#### Analysis Mode 3 : Non-negative matrix factorization(NMF) analysis , non PheWAS immage 
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
    
    
### Example run     


#### Results
![ALL_chr resultsfor_network csv_snp_bacteria_network](https://user-images.githubusercontent.com/51352117/64018795-5aaf8380-cb68-11e9-984e-d138e42b7ca3.png)

![image](https://user-images.githubusercontent.com/51352117/64018954-d3aedb00-cb68-11e9-9506-c16b18a1b5de.png)

## Reference 
 [1]

