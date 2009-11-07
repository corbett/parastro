#include "ProfileHelpers.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "DataSetHelpers.h"
#include <assert.h>
#include <cmath>
#include "vtkMath.h"
//----------------------------------------------------------------------------
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,\
											double r,double s,double xacc,double yacc,\
											int *pnIter) 
{
	// This code copied from Doug Potter's and Joachim Stadel's
	// pkdgrav, class master.c, method illinois
	// Changed to return -1 if there are problems finding root
  const int maxIter = 100;
  double t,fr,fs,ft,phis,phir,gamma;
  int i;

  fr = func(r,ctx);
  fs = func(s,ctx);
	if(fr*fs>0)
		{
			return -1;
		}
  t = (s*fr - r*fs)/(fr - fs);

  for(i=0; i<maxIter && fabs(t-s) > xacc; ++i) 
		{
		ft = func(t,ctx);
		if (fabs(ft)<=yacc)
		 {
		 break;
		 }
		if(ft*fs<0) 
			{
	   	/*
	   	** Unmodified step.
	   	*/
	   	r = s;
	   	s = t;
	   	fr = fs;
	   	fs = ft;
			}
		else 
			{
	   	/*
	   	** Modified step to make sure we do not retain the 
	   	** endpoint r indefinitely.
	   	*/
			phis = ft/fs;
	    phir = ft/fr;
	    gamma = 1 - (phis/(1-phir));  /* method 3, for true illinois gamma=0.5*/
	    if(gamma < 0) 
				{
				gamma = 0.5;
	   		}
	   	fr *= gamma;
	   	s = t;
	   	fs = ft;
			}
		t = (s*fr - r*fs)/(fr - fs);
		}	
		
  if(pnIter)
		{
		*pnIter = i;
		}
		
  return(t);
}

//----------------------------------------------------------------------------

double ComputeMaxR(vtkPointSet* input,double point[])
{
	double bounds[6]; //xmin,xmax,ymin,ymax,zmin,zmax
	input->GetPoints()->ComputeBounds();
	input->GetPoints()->GetBounds(bounds);
	double maxR=0;
	double testR=0;
	// for each of the 8 corners of the bounding box, compute the 
	// distance to the point. maxR is the max distance.
	for(int x = 0; x < 2; ++x)
		{
		for(int y = 2; y <4; ++y)
			{
			for(int z = 4; z < 6; ++z)
				{
				double testCorner[3] = {bounds[x],bounds[y],bounds[z]};
				testR = sqrt(vtkMath::Distance2BetweenPoints(testCorner,point));
				// only if our test R is greater than the current max do we update
				maxR=std::max(maxR,testR);
				}
			}
		}
	return maxR;
}

//----------------------------------------------------------------------------
double OverDensityInSphere(double r,void* inputVirialRadiusInfo)
{
	VirialRadiusInfo* virialRadiusInfo = \
	 											static_cast<VirialRadiusInfo*>(inputVirialRadiusInfo);
	vtkSmartPointer<vtkIdList> pointsInRadius = \
																vtkSmartPointer<vtkIdList>::New();
	virialRadiusInfo->locator->FindPointsWithinRadius(r,
		virialRadiusInfo->center,
		pointsInRadius);
	// calculating the average mass, dividing this by the volume of the sphere
	// to get the density
	double totalMass=0;
	vtkPointSet* dataSet=\
		vtkPointSet::SafeDownCast(
		virialRadiusInfo->locator->GetDataSet());
	for(int pointLocalId = 0; 
			pointLocalId < pointsInRadius->GetNumberOfIds(); 
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		double* nextPoint=GetPoint(dataSet,pointGlobalId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have float method
		double* mass=GetDataValue(dataSet,"mass",pointGlobalId);
		totalMass+=mass[0];
		// Finally, some memory management
		delete [] mass;
		delete [] nextPoint;
		}
	// Returning the density minus the critical density. Density is defined
	// as zero if the number points within the radius is zero
	double density = (pointsInRadius->GetNumberOfIds() > 0) ? \
		totalMass/(4./3*M_PI*pow(r,3)) : 0;
	double overdensity = density - 	virialRadiusInfo->criticalValue;
	return overdensity;
}

//----------------------------------------------------------------------------
double OverNumberInSphere(double r,void* inputVirialRadiusInfo)
{
	VirialRadiusInfo* virialRadiusInfo = \
		static_cast<VirialRadiusInfo*>(inputVirialRadiusInfo);
	vtkSmartPointer<vtkIdList> pointsInRadius = \
		vtkSmartPointer<vtkIdList>::New();
	virialRadiusInfo->locator->FindPointsWithinRadius(r,
		virialRadiusInfo->center,
		pointsInRadius);
	// Returning the number minus the critical number
	return pointsInRadius->GetNumberOfIds() - \
	 	virialRadiusInfo->criticalValue;
}

//----------------------------------------------------------------------------
VirialRadiusInfo ComputeVirialRadius(vtkPointSet* input,
	double softening,double overdensity,double maxR,double center[])
{
		// Building the point locator and the struct to use as an 
		// input to the rootfinder.
		// 1. Building the point locator
		vtkPointLocator* locator = vtkPointLocator::New();
			locator->SetDataSet(input);
			locator->BuildLocator();
		// 2. Building the struct to use as argument to root finder and density
		// functions. Contains locator, center, softening info and stores virial
		// radius info for output
		VirialRadiusInfo virialRadiusInfo;
		virialRadiusInfo.locator=locator;
		for(int i = 0; i < 3; ++i)
		{
			virialRadiusInfo.center[i]=center[i];
		}
		virialRadiusInfo.softening=softening;
		virialRadiusInfo.virialRadius = -1; // if stays -1 means not found
		// but IllinoisRootFinder takes in a void pointer
		void* pntrVirialRadiusInfo = &virialRadiusInfo;
		// 3. Define necessary variables to find virial radius, then search for 
		// it
		int numIter=0; // don't ever use this info, but root finder needs it
	 // keeps track of our guesses and their associated overdensities
	 /// initial guess is the softening
		double guessR[3]={softening,softening,softening};
		double denGuessR[3]={0,0,0}; // keeps track of the density within each R
		int fib[2]={1,1};
		while(guessR[2]<maxR)
			{
			// if our last three guesses have been monotonically decreasing
			// in density, then try to calculate the root
			if(denGuessR[0]>denGuessR[1]>denGuessR[2])
				{
				virialRadiusInfo.criticalValue=overdensity;
				virialRadiusInfo.virialRadius = \
					IllinoisRootFinder(OverDensityInSphere,
					pntrVirialRadiusInfo,
					guessR[0],guessR[2],
					softening,softening,
					&numIter);
				// we are done trying to find the root if the virial radius found is 
				// greater than zero, as rootfinder returns -1 if there were problems 
				if(virialRadiusInfo.virialRadius>0)
					{
					break;
					}
				}
			int nextFib=fib[0]+fib[1];
			// updating the fibonacci sequence
			shiftLeftUpdate(fib,2,nextFib);
			// Updating guessR
			double nextR=nextFib*virialRadiusInfo.softening;
			shiftLeftUpdate(guessR,3,nextR);
			// Updating density estimates
			// Means that OverDensityInSphere will just return DensityInSphere
			virialRadiusInfo.criticalValue=0; 
			shiftLeftUpdate(denGuessR,3,OverDensityInSphere(nextR,
				pntrVirialRadiusInfo));
		}
  	return virialRadiusInfo;
}

//----------------------------------------------------------------------------
template <class T> void shiftLeftUpdate(T* array,int size, T updateValue)
{
	// for everything but the last, value is equal to item one to right
	for(int i = 0; i < size-1; ++i)
		{
		array[i]=array[i+1];
		}
	// for last item, value is equal to updateValue
	array[size-1]=updateValue;
}
//----------------------------------------------------------------------------
vtkPolyData* CopyPolyPointsAndData(vtkPolyData* dataSet, vtkIdList*
 	pointsInRadius)
{
	// TODO: I was using CopyCells method of vtkPolyData
	// but this wasn't working so I decided to do manually
	// go back to finding the way using the VTK api to do this
	int numNewPoints=pointsInRadius->GetNumberOfIds();
	// Initilizing
	vtkPolyData* newDataSet = vtkPolyData::New(); // this memory must be managed
		// Initializing points and verts
	  newDataSet->SetPoints(vtkSmartPointer<vtkPoints>::New());
		newDataSet->SetVerts(vtkSmartPointer<vtkCellArray>::New());
	  // Initializing data
	vtkSmartPointer<vtkDataArray> nextArray;
	for(int i = 0; i < dataSet->GetPointData()->GetNumberOfArrays(); ++i)
		{
		nextArray = dataSet->GetPointData()->GetArray(i);
		AllocateDataArray(newDataSet,
			nextArray->GetName(),
			nextArray->GetNumberOfComponents(),
			numNewPoints);
		}
	// Copying

	for(int pointLocalId = 0; 
			pointLocalId < pointsInRadius->GetNumberOfIds(); 
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		double* dbNextPoint=GetPoint(dataSet,pointGlobalId);
		float* nextPoint=DoublePointToFloat(dbNextPoint);
		vtkIdType newId =SetPointValue(newDataSet,nextPoint);
		// adding this to the newDataSet
		// adding this point's data to the newDataSet, for each data array
			for(int i = 0; i < dataSet->GetPointData()->GetNumberOfArrays(); ++i)
				{
				nextArray = dataSet->GetPointData()->GetArray(i);
				double* nextData = GetDataValue(dataSet,
					nextArray->GetName(),pointGlobalId);
				SetDataValue(newDataSet,nextArray->GetName(),newId,nextData);
				delete [] nextData;
				}
			delete [] nextPoint;
			delete [] dbNextPoint;
			}
	return newDataSet;
}



//----------------------------------------------------------------------------
vtkPolyData* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo)
{

	vtkSmartPointer<vtkIdList> pointsInRadius = \
		vtkSmartPointer<vtkIdList>::New();
	virialRadiusInfo.locator->FindPointsWithinRadius(
		virialRadiusInfo.virialRadius,
		virialRadiusInfo.center,
		pointsInRadius);
  vtkPolyData* dataSet = \
		vtkPolyData::SafeDownCast(virialRadiusInfo.locator->GetDataSet());	
	// Creating a new dataset
	// first allocating
	vtkPolyData* newDataSet = \
		CopyPolyPointsAndData(dataSet,pointsInRadius);
	return newDataSet;
}


//----------------------------------------------------------------------------
double* ComputeRadialVelocity(double v[],double r[])
{
	return ComputeProjection(v,r);
}

//----------------------------------------------------------------------------
double* ComputeTangentialVelocity(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
  double* vTan=PointVectorDifference(v,vRad);
	delete [] vRad;
	return vTan;	
}
//----------------------------------------------------------------------------
double* ComputeAngularMomentum(double v[], double r[])
{
	double* angularMomentum = new double[3];
	vtkMath::Cross(v,r,angularMomentum);
	return angularMomentum;
}

//----------------------------------------------------------------------------
double* ComputeVelocitySquared(double v[],double r[])
{
	double* velocitySquared = new double[1];
	velocitySquared[0]=vtkMath::Dot(v,v);
	return velocitySquared;
}

//----------------------------------------------------------------------------
double* ComputeRadialVelocitySquared(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
	double* vRadSquared=new double[1];
	vRadSquared[0]=vtkMath::Dot(vRad,vRad);
	delete [] vRad;
	return vRadSquared;
}

//----------------------------------------------------------------------------
double* ComputeTangentialVelocitySquared(double v[],double r[])
{
	double* vRad=ComputeRadialVelocity(v,r);
	double* vTan=ComputeTangentialVelocity(v,r);
	double* vTanSquared=new double[1];
	vTanSquared[0]=vtkMath::Dot(vTan,vTan);
	delete [] vRad;
	delete [] vTan;
	return vTanSquared;	
}

//----------------------------------------------------------------------------
double* ComputeVelocityDispersion(vtkVariant vSquaredAve, vtkVariant vAve)
{
	// vSquared ave required to be a variant which holds a double,
	// vAve required to be a variant which holds a double array with 3 
	// components
	double* velocityDispersion = new double[3];
	for(int comp = 0; comp < 3; ++comp)
		{
		velocityDispersion[comp] = sqrt(fabs(vSquaredAve.ToDouble() -
			pow(vAve.ToArray()->GetVariantValue(comp).ToDouble(),2)));
		}
	return velocityDispersion;
}
//----------------------------------------------------------------------------
double* ComputeCircularVelocity(vtkVariant cumulativeMass, 
	vtkVariant binRadius)
{
	double* circularVelocity = new double[1];
	circularVelocity[0]=cumulativeMass.ToDouble()/binRadius.ToDouble();
	return circularVelocity;
}

//----------------------------------------------------------------------------
double* ComputeDensity(vtkVariant cumulativeMass, 
	vtkVariant binRadius)
{
	double* density = new double[1];
	density[0] = cumulativeMass.ToDouble()/(4./3*vtkMath::Pi()*pow(
		binRadius.ToDouble(),3));
	return density;
}

//----------------------------------------------------------------------------
double* ComputeProjection(double  vectorOne[],double vectorTwo[])
{
	double normVectorTwo = vtkMath::Norm(vectorTwo);
	double projectionMagnitude = \
		vtkMath::Dot(vectorOne,vectorTwo)/normVectorTwo;
	double* projectionVector = new double[3];
	for(int i = 0; i < 3; ++i)
	{
		projectionVector[i] = projectionMagnitude * vectorTwo[i] /normVectorTwo;
	}
	return projectionVector;
}

//----------------------------------------------------------------------------
double* PointVectorDifference(double vectorOne[], double vectorTwo[])
{
	double* pointVectorDifference = new double[3];
	for(int i = 0; i < 3; ++i)
	{
		pointVectorDifference[i] = vectorOne[i] - vectorTwo[i];
	}
	return pointVectorDifference;
}

//----------------------------------------------------------------------------
double* ComputeMidpoint(double pointOne[], double pointTwo[])
{
	double* midpoint = new double[3];
	for(int i = 0; i < 3; ++i)
	{
	midpoint[i] = (pointOne[i] + pointTwo[i])/2;
	}
	return midpoint;
}

//----------------------------------------------------------------------------
void VecMultConstant(double vector[],double constant)
{
	for(int i = 0; i < 3; ++i)
	{
		vector[i] *= constant;
	}
}


//----------------------------------------------------------------------------
	
	
	