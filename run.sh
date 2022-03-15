#!/bin/bash

kmc=./kmc # set kmc path if needed
kmc_tools=./kmc_tools #set kmc_tools if needed

#build
g++ *.cpp -std=c++20 -O3 -o matrixer

# count for example data
$kmc -k28 -ci1 -fa test1.fa o1 .
$kmc_tools transform o1 sort o1s
echo "o1s" > input.txt

$kmc -k28 -ci1 -fa test2.fa o2 .
$kmc_tools transform o2 sort o2s
echo "o2s" >> input.txt

$kmc -k28 -ci1 -fa test3.fa o3 .
$kmc_tools transform o3 sort o3s
echo "o3s" >> input.txt

# run
./matrixer input.txt > matrix.txt