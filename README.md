# Automatic pipeline of microbiome PheWAS for the discovery of host-microbe interaction networks --BETA Versions
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
#### options dicription 
##### Require Options
 `--DIR`              : path of your Plink file format data 
 
 `--Input_prefix`     : Plink files(.bed,.bim,.fam) ID
 
 `--OTU_ID`           : File ID of OTU file format       
 
 `--Bacterial_class`  : Choose bacterial taxonomic level such as Species(S),Genus(G),Family(F),Order(O),Class(C),Pylumn(P)
 
 `--Analysis`         : Choose Analysis Mode such as Linear, NMF(non-negative matrix factorization), (Logistic is not yet available) \
                       see http://zzz.bwh.harvard.edu/plink/anal.shtml 
                       
 `--P_cut`            : Cut off of Single SNP P-value base on linear Quantitative Trait Loci Wald Test.\
                       see http://zzz.bwh.harvard.edu/plink/anal.shtml .qassoc
                       
 `--P_count`          : Set the number of bacteria that exceed significance P. This is to find SNPs that control several bacteria.
 
 #### Selective Options
 
 `--PHEWAS_image_mode`: PheWAS results make image like fig.1     `Choose Y or YES make image default None`
 
 `--NMF_K`            : If the `--Analysis` option is `NMF`, set the NMF component K
 
 `--Gene_mode`        : Y or YES : SNP in the gene region 20kb are generated by gene name. default None
                        require db19_20k.gz from https://drive.google.com/open?id=1hEUdViceUQIO-_-zSShxUqW6W4qashXu
 
 `--Cov`              : Covariate File 
 
 `--Cov_names`        : "," sparate Covariate name such as age,sex,bmi, .etc
 
 `--Corr`             : Analysis method of SNP-SNP Beta correlation (defalt pearson, bray-cutis, .etc )  *not yet available option
 
 `--Corr_cut`         : Set correlation coefficient cutoff  *not yet available option

### Analysis Examples
#### Analysis example 1 : make PheWAS figures, no covariate , analysis bacterial-Class:Family
/downloaded/hGMNet/Path/hGMNet.sh \
    --OTU_ID your_OTU.txt \
    --Bacterial_class F \
    --OTU_DIR /your/OTU/path \
    --Input_prefix plink_file_id \
    --DIR /your/plink/.bed.bim.fam/path \
    --Analysis Linear \
    --P_cut 5e-6 \
    --P_count 10 \
    --PHEWAS_image_mode Y
        
#### Analysis example 2 : non make PheWAS figures, Covariate, analysis bacterial-Class:Species
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
    
#### Analysis example 3 : Non-negative matrix factorization(NMF) analysis , non PheWAS image 
/downloaded/hGMNet/Path/hGMNet.sh \
    --OTU_ID your_OTU.txt \
    --Bacterial_class F \
    --OTU_DIR /your/OTU/path \
    --Input_prefix plink_file_id \
    --DIR /your/plink/.bed.bim.fam/path \
    --Analysis NMF \
    --NMF_K 8 \
    --P_cut 5e-6 \
    --P_count 10 
    
    
### Example run     


#### --Results figure

![image](https://user-images.githubusercontent.com/51352117/64018954-d3aedb00-cb68-11e9-9506-c16b18a1b5de.png)
[fig.1 microbiome PheWAS image mode result]

![ALL_chr resultsfor_network csv_snp_bacteria_network](https://user-images.githubusercontent.com/51352117/64018795-5aaf8380-cb68-11e9-984e-d138e42b7ca3.png)
[fig.2 Bacteria and host Genotype interaction network

## Reference 
 [1]
 
 [2]  Cronin, Robert M.; Field, Julie R.; Bradford, Yuki; Shaffer, Christian M.; Carroll, Robert J.; Mosley, Jonathan D.; Bastarache,         Lisa; Edwards, Todd L.; Hebbring, Scott J. (2014). "Phenome-wide association studies demonstrating pleiotropy of genetic variants       within FTO with and without adjustment for body mass index". Frontiers in Genetics. 5: 250. doi:10.3389/fgene.2014.00250. ISSN           1664-8021. PMC 4134007. PMID 25177340.
 
 [3] Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
      PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
  
 [4] Fevotte, C., & Idier, J. (2011). Algorithms for nonnegative matrix factorization with the beta-divergence. Neural Computation, 23(9).
