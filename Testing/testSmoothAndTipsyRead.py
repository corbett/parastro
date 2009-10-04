import sys
from collections import defaultdict
#define the header files
stadelSmoothHeader=["smoothed density"]
tipsyAsciiHeader=["Points:0","Points:1","Points:2","velocity:0","velocity:1",\
						 "velocity:2","mass","potential","eps"]
paraviewHeader=["potential","mass","eps","rho","hsmooth","temperature",\
								"metals","tform","velocity:0","velocity:1","velocity:2",\
								"smoothed mass","smoothed density",\
								"Points:0","Points:1","Points:2"]

def updateDict(dict,keys,vals):
	for k,v in zip(keys,vals):
		dict[k].append(v)

def compareResults(tipsyFileAscii,paraviewFileAscii):
	if ".den" in tipsyFileAscii:
		tipsyHeader=stadelSmoothHeader
		skiplines=1
	else:
		tipsyHeader=tipsyAsciiHeader
		skiplines=3
	#read in files
	#skipping first skiplines lines
	tipsyFileAscii=open(tipsyFileAscii).readlines()[skiplines:] 
	#skipping first line
	paraviewFileAscii=open(paraviewFileAscii).readlines()[1:] 
	#make the header files dictionaries
	tipsyRead=defaultdict(list)
	paraviewRead=defaultdict(list)
  #constructing the dictionary
	for line in tipsyFileAscii:
		updateDict(tipsyRead,tipsyHeader,map(float,line.strip().split(' ')))
	for line in paraviewFileAscii:
		updateDict(paraviewRead,paraviewHeader,\
							 map(float,line.strip().split(',')))
	for var in tipsyHeader:
		tipsyList=tipsyRead[var]
		paraviewList=paraviewRead[var]
		assert(len(tipsyList)==len(paraviewList))
		for i in range(len(tipsyList)):
			assert tipsyList[i]==paraviewList[i],\
			"%f != %f"% (tipsyList[i],paraviewList[i])

if __name__ == '__main__':
	compareResults('b1.00300.d0-1000.ascii','b1.00300.d0-1000.plussmooth.csv')
	compareResults('b1.00300.d0-1000.den','b1.00300.d0-1000.plussmooth.csv')
