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

##Method 2: using just the UMU/UMM output.  No TE annotation file, no bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method2.out.gz

##Method 3: using the UMU/UMM files, the TE annotation file, and require that Ms in UMM overlap known TEs.  No bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1 --ummHitTE -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method3.out.gz

##Method 4: using the UMU/UMM files, the TE annotation file, and require that Ms in UMM overlap known TEs.  Include bam file scanning
teclust -u TEtestData/pecnv_output/pecnv_bamfile.um_u.csv.gz -m TEtestData/pecnv_output/pecnv_bamfile.um_m.csv.gz -t TEtestData/TE_position_r5.1 --ummHitTE -b TEtestData/pecnv_output/pecnv_bamfile_sorted.bam -i `pecnv_insert_qtile TEtestData/pecnv_output/pecnv_bamfile.mdist.gz 0.99` -o method4.out.gz

./checkTEoutput method2.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method2.compare.out
./checkTEoutput method3.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method3.compare.out
./checkTEoutput method4.out.gz TEtestData/TE_position_r5.1 TEtestData/line99_truth method4.compare.out