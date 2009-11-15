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
#include "vtkTable.h"
#include "vtkPointSet.h"
#include <vtkstd/string>
// Description:
// Uses the Illinois root finding method to find the root of the function
// func. The root must lie between r and s. Root is returned when it is found 
// within the accuracy xacc, yacc. pnIter indicates how many iterations the
// algorithm took to converge. Returns -1 if there are problems finding the 
// root, thus can only be used to find positive roots.
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
// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
// a locator for a given vtkdataset, a center from which to calculate
// the volume 
// Given a radius, a center, calculates the number density of within the
// sphere of radius r around the center and subtracts the critical density of 
// the user's choice.
double OverNumberInSphere(double r, void* locatorStruct);

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
	double criticalValue;
	double virialRadius;
	double softening;
	vtkstd::string massArrayName;
};


// Description:
// Computes the virial radius >=0 base upon the user defined 
// overdensity and center. Returns -1 if there is a problem. 
VirialRadiusInfo ComputeVirialRadius(vtkPointSet* input,
	vtkstd::string massArrayName, double softening,double overdensity,
	double maxR,double center[]);

// Description:
// shifts every item in array one to left (the first element is thrown away)
// then sets inserts updateValue in the last, free slot
template <class T> void shiftLeftUpdate(T* array,int size, T updateValue);

// Description:
// From dataSet, initializes a newDataSet with point and cell arrays (empty),
// and corresponding  data arrays (empty) with identical names, and number of 
// components to dataSet. Then copys, if listed in the vtkIdList,
// from the old data set to the new data set the points and their cell data.
vtkPointSet* CopyPointsAndData(vtkPointSet* dataSet, vtkIdList*
 	pointsInRadius);

// Description:
// Given a populated virialradiusinfo struct, returns a dataset corresponding
// to only those points within the virial radius.
// This method only works if input was vtkPolyData...
vtkPointSet* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo);

// Description
// Given an input data set, the bin number, a list of points in the relevant
// bin,  and the output table, computes the average radial velocity in the
// bin, the radial velocity dispersion in the bin, the tangential velocity
// in the bin, the tangentical velocity dispersion in the bin, and the angular
// momentum in the bin and adds it to the output table in the appropriate bin.
//
// Average velocity: 
// \vec v_{ave} = \frac{\sum_{i=1}^{N} \vec v_i}{N}
//
// Velocity Dispersion:
// \vec \sigma_v = \frac{\sum_{i=1}^{N} 
// \sqrt{(\vec v_i)^2-(\vec v_{ave}^2)}}{N}
//
// Specific Angular Momentum:
// \vec j = \frac{\sum_{i=1}^{N} \vec r_i x \vec v_i}{N}
//
void ComputeStatisticsInBin(vtkPolyData* inputDataSet, double center[],
	int binNum, vtkIdList* pointsInBin, vtkTable* output);
	
// Description:
// given 3-vector vectorOne and 3-vector vectorTwo, computes the 
// projection of vectorOnew onto vector two.
// Can be used to calculate e.g.
// radialVelocity = ComputeProjection(v,r);
//  or velocitySquared = ComputeProjection(v,v);
double* ComputeProjection(double  vectorOne[],double vectorTwo[]);

// Description:
// Computes the vector difference between two 3-vectors
// Can be used to calculate e.g. tangential velocity
// tangentialVelocity = PointVectorDifference(v,radialVelocity);
double* PointVectorDifference(double vectorOne[], double vectorTwo[]);

// Description
// Give a 3 vector v and a three vector r computes the specific angular 
// momentum = r x v
double* ComputeAngularMomentum(double v[], double r[]);

// Description:
// Multiplies in place a 3-vector by a constant
void VecMultConstant(double vector[],double constant);

// Description
// Given a vSquaredAve and a vAve calculates the velocity dispersion
// placing it in the output vector velocityDispersion
double* ComputeVelocityDispersion(vtkVariant vSquaredAve, vtkVariant vAve);
	
// Description:
// helper function to compute radial velocity
double* ComputeRadialVelocity(double v[],double r[]);

// Description:
// helper function to compute tangential velocity
double* ComputeTangentialVelocity(double v[],double r[]);

// Description:
// helper function to compute angular momentum
double* ComputeAngularMomentum(double v[], double r[]);

// Description:
// helper function to compute velocity squared
double* ComputeVelocitySquared(double v[],double r[]);

// Description:
// helper function to compute radial velocity squared
double* ComputeRadialVelocitySquared(double v[],double r[]);

// Description:
// helper function to compute tangential velocity squared
double* ComputeTangentialVelocitySquared(double v[],double r[]);

// Description:
// helper function to compute circular velocity
double* ComputeCircularVelocity(vtkVariant cumulativeMass, 
	vtkVariant binRadius);

// Description:
// helper function to compute density
double* ComputeDensity(vtkVariant cumulativeMass, 
	vtkVariant binRadius);
	
// Description:
// Helper function to compute the midpoint between two points
double* ComputeMidpoint(double pointOne[], double pointTwo[]);

