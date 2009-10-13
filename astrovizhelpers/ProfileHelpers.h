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
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter);

// Description:
// Given a point and an input data set, computes the maximum distance from the
// point to the data set boundaries.
// This is the maximum distance from the point to the 8 corners of the 
// bounding cube of this dataset
/*    	  
*										             (xmin,ymax,zmax)
*														       / |\
*														      /  | \
*														     /   |  \
*														   //    |   \\                							^  ^
*														  /      |     \              						z| /y
*														 /       \(xmin,ymax,zmin)    						|/
*			    (xmin,ymin,zmax)  |       / \      | (xmax,ymax,zmax)       \
*														|\     /   \    /|             						x\
*														| \   /  .p \\ / |               						v
*														|  \ /       //  |
*														|  /\       /  \ |
*														| /  \\    /    \|
*				  (xmin,ymin,zmin)	|/     \  /      | (xmax,ymax,zmin)
*														/   (xmax,ymin,zmax)  
*														 \       |      /
*														  \      |     /
*														   \\    |   //
*														     \   |  /
*														      \  | /
*														       \  /
*										            (xmax,ymin,zmin)
*
*/
double ComputeMaxR(vtkPointSet* input,double point[]);

// Description:
// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
// a locator for a given vtkdataset, a center from which to calculate
// the volume 
// Given a radius, a center, calculates the density of within the sphere
// of radius r around the center and subtracts the critical density of 
// the user's choice.
// 
double OverDensityInSphere(double r, void* locatorStruct);

// Description:
// The VirialRadiusInfo struct is an containing:
// .locator which is a vtkPointLocator
// .center  which is a double[3]
// .criticalDensity which is a double
// .virialRadius
struct VirialRadiusInfo 
{
	vtkPointLocator* locator;
	double center[3];
	double criticalDensity;
	double virialRadius;
};

// Description:
// Computes the virial radius >=0 base upon the user defined 
// overdensity and center. Returns -1 if there is a problem. 
VirialRadiusInfo ComputeVirialRadius(vtkPointSet* input,double overdensity,double center[]);

// Description:
// From dataSet, initializes a newDataSet with point and cell arrays (empty),
// and corresponding  data arrays (empty) with identical names, and number of 
// components to dataSet. Then copys, if listed in the vtkIdList,
// from the old data set to the new data set the points and their cell data.
vtkPolyData* CopyPolyPointsAndData(vtkPolyData* dataSet, vtkIdList*
 	pointsInRadius);

// Description:
// Given a populated virialradiusinfo struct, returns a dataset corresponding
// to only those points within the virial radius.
// This method only works if input was vtkPolyData...
vtkPolyData* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo);
