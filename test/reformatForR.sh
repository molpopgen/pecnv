#!/usr/bin/env bash

echo "chromo pfirst plast chrom2 mfirst mlast name score strand1 strand2 nplus nminus pdist pin mdist min"
zcat $1 | sed 's/:/\t/g' 