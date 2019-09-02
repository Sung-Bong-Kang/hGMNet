#!/bin/bash
../hGMNet.sh \
	--OTU_ID Example_OTU.txt \
	--Bacterial_class F \
	--OTU_DIR ../SetUP_example \
	--Input_prefix Example \
	--DIR ../SetUP_example \
	--Analysis Linear \
	--P_cut 5e-4 \
	--P_count 31 \
	--Gene_mode N \
	--PHEWAS_image_mode N 
