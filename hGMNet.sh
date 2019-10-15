#!/bin/bash 
script_DIR="$( cd "$( dirname "$0" )" && pwd -P )"
if ! options=$(getopt -o h -l help,DIR:,Input_prefix:,OTU_DIR:,OTU_ID:,Bacterial_class:,P_cut:,P_count:,PHEWAS_image_mode:,Corr:,Corr_cut:,Analysis:,Beta_cut:,Gene_mode:,NMF_K:,Cov:,Cov_names:,Norm: -- "$@")
then
	echo "ERROR: print usage"
	exit 1
fi
eval set -- "$options"
while true; do
	case "$1" in
		-h|--help)
			$script_DIR/Help
			break
			shift ;;
		--DIR)
			DIR="$2"
			echo --DIR == $DIR
		shift 2 ;;
		--Input_prefix)
			Input_prefix="$2"
			echo --Input_prefix == $Input_prefix
		shift 2 ;;	
		--OTU_DIR)
			OTU_DIR="$2"
			echo --OTU_DIR == $OTU_DIR
		shift 2 ;;	
		--OTU_ID)
			OTU_ID="$2"
			echo --OTU_ID == $OTU_ID
		shift 2 ;;

		--Bacterial_class)
			Bacterial_class="$2"
			echo --Bacterial_class == $Bacterial_class
		shift 2 ;;
		--P_cut)
			P_cut="$2"
			echo --P_cut == $P_cut
		shift 2 ;;
		--P_count)
			P_count="$2"
			echo --P_count == $P_count
		shift 2 ;;
		--PHEWAS_image_mode)
			PHEWAS_image_mode="$2"
			echo --PHEWAS_image_mode == $PHEWAS_image_mode
		shift 2 ;;
		--Corr)
			Corr="$2"
			echo --Corr == $Corr
		shift 2 ;;
		--Corr_cut)
			Corr_cut="$2"
			echo --Corr_cut == $Corr_cut
		shift 2 ;;
		--Analysis)
			Analysis="$2"
			echo --Analysis == $Analysis
		shift 2 ;;	
		--Beta_cut)
			Beta_cut="$2"
			echo --Beta_cut == $Beta_cut
		shift 2 ;;
		--Gene_mode)
			Gene_mode="$2"
			echo --Gene_mode == $Gene_mode
			
		shift 2 ;;
		--NMF_K)
			NMF_K="$2"
			echo --NMF-K == $NMF_K
		shift 2 ;;
		--Cov) 
			Cov="$2"
			echo --Cov == $Cov
		shift 2 ;;
		############ comma split Covar names #####
		--Cov_names)
			Cov_names="$2"
			echo --Cov_name == $Cov_names
		shift 2 ;;
		--Norm)
			Norm="$2"
			echo --Norm == $Norm
		shift 2 ;;
		--)
		shift
		break
	esac
done
echo "$@"
################################################   OTU level Down  ########################################################### 
if [ "$Bacterial_class" == '' ];
then 
	echo "put in --Bacterial_class { K(kingdom), P(phylum), C(class), O(order), F(family), G(genus), S(species)}"
else
	if [ ${Bacterial_class} == 'K' ];then Bacterial_class=1 Bacteria_idx='k__'
	elif [ ${Bacterial_class} == 'P' ];then Bacterial_class=2 Bacteria_idx='p__'
	elif [ ${Bacterial_class} == 'C' ];then	Bacterial_class=3 Bacteria_idx='c__'
	elif [ ${Bacterial_class} == 'O' ];then	Bacterial_class=4 Bacteria_idx='o__'
	elif [ ${Bacterial_class} == 'F' ];then	Bacterial_class=5 Bacteria_idx='f__'
	elif [ ${Bacterial_class} == 'G' ];then	Bacterial_class=6 Bacteria_idx='g__'
	elif [ ${Bacterial_class} == 'S' ];then	Bacterial_class=7 Bacteria_idx='s__'
	else echo 'put  K(kingdom), P(phylum), C(class), O(order), F(family), G(genus), S(species)' 
	fi;
fi
echo $Bacterial_class $Bacteria_idx
#######################  Script directory ######################################
if [ "$Norm" == "CSS" ];
then
	python $script_DIR/otu_level_down.py $OTU_DIR/$OTU_ID $Bacterial_class
	Rscript $script_DIR/CSS_norm.R $OTU_DIR/${OTU_ID%%".txt"}"_cluster_"$Bacterial_class"_level.txt" $OTU_DIR/${OTU_ID%%".txt"}"_cluster_"$Bacterial_class"_level."CSS.txt
	OTU_ID=${OTU_ID%%".txt"}"_cluster_"$Bacterial_class"_level."CSS.txt
	OTU_T=$OTU_DIR/${OTU_ID%%".txt"}"_cluster_"$Bacterial_class"_level."CSS.txt

	
else
	## Default 
	Rscript $script_DIR/Log_TSS_normalization_minmax.R $OTU_DIR/$OTU_ID $OTU_DIR/${OTU_ID%%".txt"}'.LogTSSMinMax.txt' $OTU_DIR/${OTU_ID%%".txt"}'.TSS.txt'
	python $script_DIR/otu_level_down.py $OTU_DIR/${OTU_ID%%".txt"}'.TSS.txt' $Bacterial_class
	OTU_T=$OTU_DIR/${OTU_ID%%".txt"}'.TSS_cluster_'$Bacterial_class'_level.txt'
	echo $OTU_T
	OTU_ID=${OTU_ID%%".txt"}'.LogTSSMinMax.txt'
	python $script_DIR/otu_level_down.py $OTU_DIR/$OTU_ID  $Bacterial_class 
	OTU_ID=${OTU_ID%%".txt"}'_cluster_'$Bacterial_class'_level.txt'
sed -e 's/ /_/' $OTU_DIR/$OTU_ID  >$OTU_DIR/temp
mv $OTU_DIR/temp $OTU_DIR/$OTU_ID
fi
echo $OTU_ID
if [ "$Analysis" == "NMF" ];
then
	echo $Analysis
	python $script_DIR/NMF.py $OTU_DIR/$OTU_ID $NMF_K
	python $script_DIR/OTU_to_PHENO.py $OTU_DIR/$OTU_ID'_H' $DIR/$Input_prefix.fam
	cp $OTU_DIR/$OTU_ID'_W' $DIR/
	phenos=$DIR/$OTU_ID'_H'.pheno
else 
	python $script_DIR/OTU_to_PHENO.py $OTU_DIR/$OTU_ID $DIR/$Input_prefix.fam
	phenos=$DIR/$OTU_ID.pheno 
	
fi
## microbial abundance table phenotype format
## phenome File : FID IID Bacteria1 Bacteria2 ....
###################### Linear QTL analysis ####################################
##############################################################################
if [ ! -d $DIR/Qassoc ];
then
	mkdir $DIR/Qassoc
else
	rm $DIR/Qassoc/*
fi
## Linear or Logistic or NMF	
if [ "$Analysis" == "Linear" -a "$Cov" == "" ];
then
	echo Analysis mode is Linear models 
	for idx in {1..22};
	do
		plink --bfile $DIR/$Input_prefix --chr $idx --pheno $phenos --all-pheno --allow-no-sex --assoc --out $DIR/Qassoc/chr$idx.OTU 
		for B_idx in `ls $DIR/Qassoc |grep chr$idx.OTU.*$Bacteria_idx|grep -v $Bacteria_idx'.qassoc'`;
		do
			B_idx=${B_idx##/*/}
			awk -v DDD=$B_idx '{print DDD"\t"$2"\t"$5"\t"$9}' $DIR/Qassoc/$B_idx|grep -v NA
		done | sort -k 2 |grep -v -E "SNP.*P|P.*SNP" >$DIR/chr$idx.results &
		echo Chr$idx Qassoc Test complet
	done
	wait 
	echo END
elif [ "$Analysis" == 'NMF' -a "$Cov" == "" ];
then
	for idx in {1..22};
	do
		plink --bfile $DIR/$Input_prefix --chr $idx --pheno $phenos --all-pheno --allow-no-sex --assoc --out $DIR/Qassoc/chr$idx.OTU
		for B_idx in `ls $DIR/Qassoc |grep chr$idx.OTU.*.qassoc`;
		do
			B_idx=${B_idx##/*/}
			awk -v DDD=$B_idx '{print DDD"\t"$2"\t"$5"\t"$9}' $DIR/Qassoc/$B_idx|grep -v NA
		done | sort -k 2 |grep -v -E "SNP.*P|P.*SNP" >$DIR/chr$idx.results &
	done
	wait
	echo END
elif [ "$Analysis" == "Linear" -a "$Cov" != "" ];
then
	VARS=$(echo $Cov_names | sed -e 's/,/ /g')
	for idx in {1..22};
	do
		plink --bfile $DIR/$Input_prefix --chr $idx --pheno $phenos --all-pheno --allow-no-sex --covar $Cov --covar-name $VARS --parameters 1 --linear --out $DIR/Qassoc/chr$idx.OTU
		for B_idx in `ls $DIR/Qassoc |grep chr$idx.OTU.*$Bacteria_idx|grep -v $Bacteria_idx'.assoc.linear'`;
		do
			B_idx=${B_idx##/*/}
			awk -v DDD=$B_idx '{print DDD"\t"$2"\t"$7"\t"$9}' $DIR/Qassoc/$B_idx|grep -v NA
		done | sort -k 2 |grep -v -E "SNP.*P|P.*SNP" >$DIR/chr$idx.results &
		echo Chr$idx --COV Qassoc Test complet
	done
	wait
	echo END
elif [ "$Analysis" == 'NMF' -a "$Cov" != "" ];
then
	VARS=$(echo $Cov_names | sed -e 's/,/ /g')
	for idx in {1..22};
	do
		plink --bfile $DIR/$Input_prefix --chr $idx --pheno $phenos --all-pheno --allow-no-sex --covar $Cov --covar-name $VARS --assoc --out $DIR/Qassoc/chr$idx.OTU
		for B_idx in `ls $DIR/Qassoc |grep chr$idx.OTU.*.qassoc`;
		do
			B_idx=${B_idx##/*/}
			awk -v DDD=$B_idx '{print DDD"\t"$2"\t"$5"\t"$9}' $DIR/Qassoc/$B_idx|grep -v NA
		done | sort -k 2 |grep -v -E "SNP.*P|P.*SNP" >$DIR/chr$idx.results &
	done
	wait
	echo END
fi
cat $DIR/chr*.results > $DIR/ALL_chr.results
echo -e "chrEND.OTU.END;END.qassoc\trsEND\t10000.0\t10000.0" >>$DIR/ALL_chr.results
rm $DIR/chr*.results
##########################################################PHEAWS#########################################
if [ "$PHEWAS_image_mode" != 'Y' -a "$PHEWAS_image_mode" != 'YES' -o "$Analysis" == 'NMF' ]; 
then
	python $script_DIR/PHEWAS_noPNG.py $DIR/ALL_chr.results $P_cut $P_count $DIR
else 
	if [ ! -d $DIR/PHEWAS_fig ];
	then
		mkdir $DIR/PHEWAS_fig
	else
		rm $DIR/PHEWAS_fig/*
	fi
	python $script_DIR/PHEWAS.py $DIR/ALL_chr.results $P_cut $P_count $DIR 
	mv $DIR/rs*.png $DIR/PHEWAS_fig/
fi
################################ NETWORK Analysis ############################################################ 
echo $Gene_mode
if [ "$Gene_mode" != 'TRUE' -a "$Analysis" != 'NMF' ];
then	
	python $script_DIR/hGMNet_NETWORK.py $DIR/ALL_chr.resultsfor_network.csv $OTU_T $Bacteria_idx
elif [ "$Gene_mode" != 'TRUE' -a "$Analysis" == 'NMF' ];
then
	python $script_DIR/hGMNet_NMF_network.py $DIR/ALL_chr.resultsfor_network.csv $DIR/$OTU_ID'_W' $Bacteria_idx
else
	echo Gene mdoe start! 	
	head -n 1 $DIR/ALL_chr.resultsfor_network.csv |sed -e 's/,/\n/g' |sed 's/Family//g'|sed '/^$/d'|awk '{print $1"\t"0.01}' > $DIR/temp 
	$script_DIR/GeneP.pl $script_DIR/db19_20k.gz 1 < $DIR/temp > $DIR/temp.gene 
	python $script_DIR/G2N.py $DIR/temp.gene $DIR/ALL_chr.resultsfor_network.csv
	python $script_DIR/hGMNet_NETWORK.py $DIR/ALL_chr.resultsfor_network_GENE.csv $OTU_T $Bacteria_idx
	rm $DIR/temp*	
fi
