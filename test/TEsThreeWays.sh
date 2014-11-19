#!/usr/bin/env bash

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

