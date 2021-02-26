#!/bin/bash

MENDLIST=$1

echo "period,numtrios,num_consistent,num_inf_consistent,num_inconsistent,mi_all,mi_inf" | sed 's/,/\t/g'
for period in $(seq 1 6)
do
    cat $(echo $MENDLIST | sed 's/,//g') | \
	awk -v "period=$period" '(length($3)==period)' | \
	datamash sum 4 sum 5 sum 6 sum 7 | \
	awk -v"period=$period" '{print period "\t" $0 "\t" $2/$1 "\t" $3/($3+$4)}'
done
