#!/usr/bin/env bash

echo "creating library fio"
#static
gcc -c romberg.c fio.c
ar cr libFio.a fio.o romberg.o 
#ranlib libFio.a
#shared
#gcc -c -fPIC romberg.c fio.c 
#gcc -shared -o libFio.so fio.o romberg.o
FIO_LIB_DIR=`pwd`

echo "testing library and demoing fio functionality"
gcc test_lib.c -lm -lFio -L$FIO_LIB_DIR  -o test_lib

echo "creating grafic2csv"
gcc grafic2csv.c -lm -lFio -L$FIO_LIB_DIR  -o grafic2csv

echo "testing grafic2csv"
./grafic2csv /Users/corbett/Local/BigData/GraficIcs/ic_files > testsim32cubed.csv