#include "vtkPolyData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkDataSetAttributes.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkPointLocator.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
// Description:
// Uses the Illinois root finding method to find the root of the function
// func. The root must lie between r and s. Root is returned when it is found 
// within the accuracy xacc, yacc. pnIter indicates how many iterations the
// algorithm took to converge.
// have to wrap as it is a template function.
// Returns -1 if there is a problem, thus is designed to ONLY find positive 
// values as root
//BTX
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter);
//ETX
// Description:
// This assumes locatorStruct is a LocatorInfo struct, which contains
// a locator for a given vtkdataset, a center from which to calculate
// the volume 
// Given a radius, a center, calculates the density of within the sphere
// of radius r around the center and subtracts the critical density of 
// the user's choice.
// 
double OverDensityInSphere(double r, void* locatorStruct);

// Description:
// The locator struct is an pointLocatorStruct containing:
// .locator which is a vtkPointLocator
// .center  which is a double[3]
// .criticalDensity which is a double
struct LocatorInfo 
{
	vtkPointLocator* locator;
	double center[3];
	double criticalDensity;
};
