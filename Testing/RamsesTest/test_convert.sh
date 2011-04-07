#!/usr/bin/env sh

for runtype in Nohydro Hydro
do
	datadir=$runtype/output_00036/
	cd $datadir
	infofile=`ls info_*.txt`
  /Users/corbett/Documents/Projects/ccprojects/RamsesRead++/ramses2tipsy $runtype/output_00036/$infofile
	cd -
done
