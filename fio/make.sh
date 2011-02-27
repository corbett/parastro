#!/usr/bin/env bash
echo "building fio"
#gcc romberg.c fio.c -o fio
gcc test_lib.c -lGraficHelpers -L/Users/corbett/Documents/Projects/pvaddons/ParaViz/lib/ -o test_lib
