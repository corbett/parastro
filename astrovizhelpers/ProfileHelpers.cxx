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
//----------------------------------------------------------------------------
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,\
											double r,double s,double xacc,double yacc,\
											int *pnIter) 
{
	// This code copied from Doug Potter's and Joachim Stadel's
	// pkdgrav, class master.c, method illinois
  const int maxIter = 100;
  double t,fr,fs,ft,phis,phir,gamma;
  int i;

  fr = func(r,ctx);
  fs = func(s,ctx);
	cout << "fr and fs are " << fr << " " << fs << "\n";
  assert(fr*fs < 0);
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

double OverDensityInSphere(double r,void* inputLocatorInfo)
{
	cout << "overdensity called with r " << r << "\n";
	LocatorInfo* locatorInfo = static_cast<LocatorInfo*>(inputLocatorInfo);
	cout <<"locator info says delta is " \
								<< locatorInfo->criticalDensity << "\n";
	vtkSmartPointer<vtkIdList> pointsInRadius = \
																vtkSmartPointer<vtkIdList>::New();
	locatorInfo->locator->FindPointsWithinRadius(r,locatorInfo->center,\
																							pointsInRadius);
	// calculating the average mass, dividing this by the volume of the sphere
	// to get the density
	double totalMass;
	for(int pointLocalId = 0; \
			pointLocalId < pointsInRadius->GetNumberOfIds(); \
			++pointLocalId)
		{
		vtkIdType pointGlobalId = pointsInRadius->GetId(pointLocalId);
		vtkPointSet* dataSet=\
								vtkPointSet::SafeDownCast(locatorInfo->locator->GetDataSet());
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
		return totalMass/(4./3*M_PI*pow(r,3)) - locatorInfo->criticalDensity;
}
	
	
	
	
	
	
	
	
	