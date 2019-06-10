#!/bin/bash


cut sampleTableV5.tsv -f5 | sort -u | while read t
do
	RANDOM=46
	cut sampleTableV5.tsv -f-5| grep $t - |awk '$3 == "y" {print $1}' | shuf -n 5 | grep -Ff - sampleTableV5.tsv 
done | grep -Ff ref/subtissues_dev.txt - > sampleTableDev.tsv


