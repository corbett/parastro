#!/usr/bin/env sh

for runtype in Nohydro Hydro
do
	datadir=$runtype
	cd $datadir
	~/Projects/ccprojects/Ramses2Tipsy/output2tipsy_v2 -inp output_00036 -out output00036.ascii
	cd -
done
