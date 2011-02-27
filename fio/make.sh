#!/usr/bin/env bash

echo "creating library fio"
#static
gcc -c romberg.c fio.c
ar rvs libFio.a fio.o romberg.o 
#shared
gcc -c -fPIC romberg.c fio.c 
gcc -shared -o libFio.so fio.o romberg.o

echo "testing library and demoing fio functionality"
gcc test_lib.c -lfio -L/Users/corbett/Documents/Projects/pvaddons/ParaViz/ParaViz_src/fio/test  -o test_lib
