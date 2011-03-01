#!/usr/bin/env bash

echo "creating grafic2csv"

gcc -lm romberg.c fio.c grafic2csv.c -o grafic2csv


echo "testing grafic2csv"
./grafic2csv /project/s201/corbett/GlobularClusters/testbox_32/ic_files/ > testsim32cubed.csv