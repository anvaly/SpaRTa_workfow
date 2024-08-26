#!/bin/bash

if [ $1 == nhops ] ; then
    for i in 1 2 3 4 5 6 7 8 9 ;
    do
	python ./nhops_subgraphs.py -n $i
    done
fi

if [ $1 == nsize ] ; then
    for i in 1 2 3 5 7 9 ;
    do
	for j in 20 30 40 50 60 ;
	do
	    python ./nhops_subgraphs.py -n $i -s $j
	done
    done
fi

	
