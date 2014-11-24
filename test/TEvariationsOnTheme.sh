#!/usr/bin/env bash

if [ ! -e run_TEtest_data.sh ]
then
    echo "I want to be run from the directory containing run_TEtest_data.sh.  I'm not finding that script here, so I'm exiting"
    exit
fi

if [ -z teclust ]
then
    >&2 echo "teclust not found"
    exit
else
    echo "found "`which teclust`
fi    

if [ ! -s TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz ]
then
    echo "Looks like you need to run run_TEtest_data.sh first!"
exit
fi

if [ ! -s TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz ]
then
    echo "Looks like you need to run run_TEtest_data.sh first!"
exit
fi

##Let's redo what run_TE_test_data.sh does, but with a different TE positions file
teclust -b TEtestData/pecnv_output/pecnv_bamfile_sorted.bam -t TEtestData/TE_position_r5.1_KRT.bed -o method1b.out.gz -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99`

##Method 2: using just the UMU/UMM output.  No TE annotation file, no bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method2.out.gz

##Method 3: using the UMU/UMM files, the TE annotation file, and require that Ms in UMM overlap known TEs.  No bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1 --ummHitTE -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method3.out.gz
##Method 3b: different TE positions file
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1_KRT.bed --ummHitTE -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method3b.out.gz

##Method 4: using the UMU/UMM files, the TE annotation file, and require that Ms in UMM overlap known TEs.  Include bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1 --ummHitTE -b TEtestData/pecnv_output/pecnv_bamfile_sorted.bam -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method4.out.gz
##Method 4b: different positions files
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1_KRT.bed --ummHitTE -b TEtestData/pecnv_output/pecnv_bamfile_sorted.bam -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method4b.out.gz

##Method 1 = Cridland et al approach
./checkTEoutput.R TEtestData/teclust_output.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method1.compare.out
./checkTEoutput.R method1b.out.gz TEtestData/TE_position_r5.1_KRT.bed TEtestData/line99_truth method1b.compare.out
##Methods 2 = Cridland approach w/o agressive scanning for more data
./checkTEoutput.R method2.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method2.compare.out
#Method 3 = What we did for DGRP (Mackay et al. 2012, Nature), more or less
./checkTEoutput.R method3.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method3.compare.out
./checkTEoutput.R method3b.out.gz TEtestData/TE_position_r5.1_KRT.bed TEtestData/line99_truth method3b.compare.out
#Method 4 = Method 3 + bam file scanning
./checkTEoutput.R method4.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method4.compare.out
./checkTEoutput.R method4b.out.gz TEtestData/TE_position_r5.1_KRT.bed TEtestData/line99_truth method4b.compare.out