import sys
def calculateWeightedMass(mass,point):
	return map(lambda x: mass*x,point)

def testCOM(csvFile):
	f=open(csvFile).readlines()[1:] #skipping first line
	totalMass=0
	totalWeightedMass=[0,0,0]
	for line in f:
		line=map(float,line.strip().split(','))
		mass=line[1]
		point=line[13:16]
		weightedMass=calculateWeightedMass(mass,point)
		totalMass+=mass	
		for i in range(3):
			totalWeightedMass[i]+=weightedMass[i]
	centerOfMass=map(lambda x: x/totalMass,totalWeightedMass)
	print "total mass=" + str(totalMass)
	print "total weighted mass="+str(totalWeightedMass)
	print "center of mass="+str(centerOfMass)

if __name__ == '__main__':
	testCOM('b1.00300.d0-1000.plussmooth.csv')
	
	