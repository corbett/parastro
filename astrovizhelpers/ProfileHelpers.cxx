#include "ProfileHelpers.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
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
	// NOTE: Only use for positive roots. 
  const int maxIter = 100;
  double t,fr,fs,ft,phis,phir,gamma;
  int i;

  fr = func(r,ctx);
  fs = func(s,ctx);
	if(fr*fs>0)
		{
			// used to be an assert, but removing 
			throw "something went wrong with the root finding";
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
	virialRadiusInfo->locator->FindPointsWithinRadius(r,\
																										virialRadiusInfo->center,\
																										pointsInRadius);
	// calculating the average mass, dividing this by the volume of the sphere
	// to get the density
	double totalMass=0;
	vtkPointSet* dataSet=\
								vtkPointSet::SafeDownCast(\
															virialRadiusInfo->locator->GetDataSet());
	for(int pointLocalId = 0; \
			pointLocalId < pointsInRadius->GetNumberOfIds(); \
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		double* nextPoint=GetPoint(dataSet,pointGlobalId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		double* mass=GetDataValue(dataSet,\
															"mass",pointGlobalId);
		totalMass+=mass[0];
		// Finally, some memory management
		delete [] mass;
		delete [] nextPoint;
		}
	// Returning the density minus the critical density
	return totalMass/(4./3*M_PI*pow(r,3)) - virialRadiusInfo->criticalDensity;
}

VirialRadiusInfo ComputeVirialRadius(vtkPointSet* input,\
																		double overdensity,double center[])
{
		// calculating the the max an min r of this pointset
		double maxR=ComputeMaxR(input,center);
		// Building the point locator and the struct to use as an 
		// input to the rootfinder.
		// 1. Building the point locator
		vtkPointLocator* locator = vtkPointLocator::New();
			locator->SetDataSet(input);
			locator->BuildLocator();
		// 2. Building the struct
		VirialRadiusInfo virialRadiusInfo;
		virialRadiusInfo.locator=locator;
		// copies the contents of center to virialRadiusInfo's arg center
		for(int i = 0; i < 3; ++i)
		{
			virialRadiusInfo.center[i]=center[i];
		}
		virialRadiusInfo.criticalDensity=overdensity;
		// but IllinoisRootFinder takes in a void pointer
		void* pntrVirialRadiusInfo = &virialRadiusInfo;
		// 3. Now we are ready to run the root finder
		int numIter=0;
		try
			{
			virialRadiusInfo.virialRadius=IllinoisRootFinder(OverDensityInSphere,\
																				pntrVirialRadiusInfo,\
																				maxR,1e-11f,//minR is almost zero
																				0.0,0.0,
																			  &numIter);
			}
		catch (const char* e)
			{
				// This indicates that something has gone wrong with the root finding
				virialRadiusInfo.virialRadius=-1; 
			}
  	return virialRadiusInfo;
}

//----------------------------------------------------------------------------
vtkPolyData* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo)
{

	vtkSmartPointer<vtkIdList> pointsInRadius = \
																vtkSmartPointer<vtkIdList>::New();
	virialRadiusInfo.locator->FindPointsWithinRadius(\
															virialRadiusInfo.virialRadius,\
															virialRadiusInfo.center,\
															pointsInRadius);
  vtkPolyData* dataSet = \
					vtkPolyData::SafeDownCast(virialRadiusInfo.locator->GetDataSet());	
	// Creating a new dataset

	vtkSmartPointer<vtkPolyData> newDataSet = \
																	vtkSmartPointer<vtkPolyData>::New();
	  newDataSet->SetPoints(vtkSmartPointer<vtkPoints>::New());
		newDataSet->SetVerts(vtkSmartPointer<vtkCellArray>::New());
	// Copy cells listed in idList from pd, including points, point data, 
	// and cell data. This method assumes that point and cell data have been
	// allocated.
	newDataSet->CopyCells(dataSet,pointsInRadius);
	return newDataSet;
}
	
	
	
	
	
	
	
	
	