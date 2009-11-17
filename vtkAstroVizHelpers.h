// TODO: this class will replace DataSetHelpers and ProfileHelpers,
// wrapping them in a class in the vtk coding style. will also contain
// all mpi enums used throughout astroviz for communication. Currently in
// progress and unused, will make the swap when I have time.
/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkAstroVizHelpers.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkAstroVizHelpers - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkAstroVizHelpers 
// Contains many static helper functions used throughout astroviz

#ifndef __vtkAstroVizHelpers_h
#define __vtkAstroVizHelpers_h
#include "vtkObject.h"
#include <vtkstd/string> // some functions use this as argument
#include <iostream> // needed for cstring/vtkstring conversion
#include <sstream> // needed for cstring/vtkstring conversion
// Forward declarations in the vtk style
class vtkCell;
class vtkCellArray;
class vtkDataSetAttributes;
class vtkDataArray;
class vtkFloatArray;
class vtkFieldData;
class vtkIdTypeArray;
class vtkIdList;
class vtkPolyData;
class vtkPointData;
class vtkPointLocator;
class vtkPointSet;
class vtkTable;
class vtkInformationVector;
class vtkTable;
class vtkMultiProcessController;
class VTK_COMMON_EXPORT vtkAstroVizHelpers : public vtkObject
{
public:
	static vtkAstroVizHelpers *New();
	vtkTypeRevisionMacro(vtkAstroVizHelpers, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
	// TODO: include all MPI enums here in a single enum
	enum PointsInRadiusMPIData
	{
		TOTAL_MASS_IN_SPHERE,
		TOTAL_NUMBER_IN_SPHERE,
		MAX_R
	};	
	
	// Description:
	// The VirialRadiusInfo struct is an containing:
	// .locator which is a vtkPointLocator
	// .center  which is a double[3]
	// .criticalDensity which is a double
	// .virialRadius
	struct VirialRadiusInfo 
	{
		vtkPointLocator* locator;
		vtkMultiProcessController* controller;
		double center[3];
		double criticalValue;
		double virialRadius;
		double softening;
		vtkstd::string massArrayName;
	};
	
	/*
	* The following methods take and modify vtkPolyData
	*/
	// Description:
	// sets the point vertices in the output vector, assigning the point 
	// a unique id and places one point per cell. Works only with PolyData
	static vtkIdType SetPointValue(vtkPolyData* output,float pos[]);

	// Description:
	// creates a sphere of radius R, center center in the output.
	static void CreateSphere(vtkPolyData* output,double radius,double center[]);

	/*
	* The following methods take and modify vtkPointSet data
	*/
	// Description:
	// sets the point vertices in the output vector, assigning the point 
	// a unique. 
	static vtkIdType SetPointValue(vtkPointSet* output,float pos[]);

	// Description:
	// sets the data value in the output vector in array arrayName at 
	// position id to data.
	static void SetDataValue(vtkPointSet* output, const char* arrayName,
		vtkIdType id,float data[]);
	static void SetDataValue(vtkPointSet* output, const char* arrayName,
		vtkIdType id,double data[]);

	// Description:
	// points the double* data to the value in the output vector in 
	// array arrayName at position id would like to be float, but the 
	// version of vtk I am working with does not support 
	static double* GetDataValue(vtkPointSet* output, const char* arrayName,
		vtkIdType id);

	// Description:
	// create a vtkDataArray of floats with the  name arrayName, number /
	// of components. place it in the vtkPointSet
	static void AllocateDataArray(vtkPointSet* output, const char* arrayName,
		int numComponents, int numTuples);

	// Description:
	// create a vtkDataArray with the  name arrayName, number of components 
	// numComponents and number of tuples numTuples of type T. place it in
	// the vtkTable
	static void AllocateDataArray(vtkTable* output, const char* arrayName,
		int numComponents, int numTuples);

	// Description:
	// create a vtkDataArray with the  name arrayName, number of components 
	// numComponents and number of tuples numTuples of type T
	// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
	// scalar double densities
	static void AllocateDoubleDataArray(vtkPointSet* output, 
		const char* arrayName, int numComponents, int numTuples);

	// Description:
	// create a vtkDataArray with the  name arrayName, number of components 
	// numComponents and number of tuples numTuples of type T
	// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
	// scalar double densities
	static void AllocateIntDataArray(vtkPointSet* output, const char* arrayName,
		int numComponents, int numTuples);

	// Description:
	// Given a VTK array, sets it to have arrayName, numComponents and numTuples
	// unlike the above methods, does not create the array, nor add it to the 
	// output	
	static void InitializeDataArray(vtkDataArray* dataArray, 
		const char* arrayName, int numComponents, int numTuples);

	// Description:
	// returns a pointer to the point's coordinates in output which corresponds 
	// to this id
	static double* GetPoint(vtkPointSet* output,vtkIdType id);

	// Description:
	// takes in a double array of size three representing a point
	// and converts it to an array of the same size but in float precision
	static float* DoublePointToFloat(double point[]);

	// Description:
	// Creates a new information vector that is a deep copy of the old one
	// Note, returned vector of pointers to vtkInformationVectors must be 
	// deleted by caller. See DeleteDeepCopyInput below.
	static vtkInformationVector** DeepCopyInputVector(
		vtkInformationVector** inputVector,int inputVectorSize); 

	// Description: 
	// Deletes an array of pointers to vtkInformationVectors
	static void DeleteDeepCopyInputVector(vtkInformationVector** inputVector, 
		int inputVectorSize);

	// Description:
	// converts a double to a vtkstd::string
	static inline vtkstd::string ToString(double x)
	{
	  std::ostringstream o;
		o << x;
	  return o.str();
	}

	// Description:
	// returns true if this process should be run in parallel
	// (i.e. we have  non-null controller and more than one process to 
	// work with)
	static bool RunInParallel(vtkMultiProcessController* controller);

	// Description:
	// Uses the Illinois root finding method to find the root of the function
	// func. The root must lie between r and s. Root is returned when it 
	// is found  within the accuracy xacc, yacc. pnIter indicates how 
	// many iterations the algorithm took to converge. Returns -1 if there are
	// problems finding the  root, thus can only be used to find positive 
	// roots
	static double IllinoisRootFinder(double (*func)(double,void *),void *ctx,
		double r,double s,double xacc,double yacc,int *pnIter);

	// Description:
	// Given a point and an input data set, computes the maximum distance 
	// from the point to the data set boundaries.
	// This is the maximum distance from the point to the 8 corners of the 
	// bounding cube of this dataset
	//   	  
	//									             (xmin,ymax,zmax)
	//													       / |\
	//													      /  | \
	//													     /   |  \
	//													   //    |   \\                							^  ^
	//													  /      |     \              						z| /y
	//													 /       \(xmin,ymax,zmin)    						|/
	//		    (xmin,ymin,zmax)  |       / \      | (xmax,ymax,zmax)       \
	//													|\     /   \    /|             						x\
	//													| \   /  .p \\ / |               						v
	//													|  \ /       //  |
	//													|  /\       /  \ |
	//													| /  \\    /    \|
	//			  (xmin,ymin,zmin)	|/     \  /      | (xmax,ymax,zmin)
	//													/   (xmax,ymin,zmax)  
	//													 \       |      /
	//													  \      |     /
	//													   \\    |   //
	//													     \   |  /
	//													      \  | /
	//													       \  /
	//									            (xmax,ymin,zmin)
	//
	static double ComputeMaxR(vtkPointSet* input,double point[]);

	// Description:
	// Runs ComputeMaxR on each process, then syncronizes results taking the
	// over all global maximum. This result will be slightly different that that
	// if run serially, as the bounding boxes may be tighter on subprocesses
	// than an overall bounding box. But if anything, the parallel 
	// version should be an even better approximation to the maximum radius 
	// (as we use the bounding box, not the individual points therein to
	// compute) than the serial version, and is always guarenteed to be
	// greater than or equal to the true maximum (if one used the points 
	// within the data set)
	static double ComputeMaxRadiusInParallel(
		vtkMultiProcessController* controller,vtkPointSet* input,double point[]);
	
	// Description:
	// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
	// a locator for a given vtkdataset, a center from which to calculate
	// the volume 
	// Given a radius, a center, calculates the density of within the sphere
	// of radius r around the center and subtracts the critical density of 
	// the user's choice.
	// 
	static double OverDensityInSphere(double r, void* locatorStruct);

	// Description:
	// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
	// a locator for a given vtkdataset, a center from which to calculate
	// the volume 
	// Given a radius, a center, calculates the number density of within the
	// sphere of radius r around the center and subtracts the critical d
	// density of the user's choice.
	static double OverNumberInSphere(double r, void* locatorStruct);

	// Descriptions:
	// If this is called by each process in the domain of the controller, finds
	// all points within a given radius on this process
	static vtkIdList* FindPointsWithinRadius(double r, double* center,
		vtkPointLocator* locatorOfThisProcess);

	// Description:
	// Computes the virial radius >=0 base upon the user defined 
	// overdensity and center. Returns -1 if there is a problem. 
	// Works in parallel if a controller is specified not equal to null and if 
	// the number of processors is > 1
	static VirialRadiusInfo ComputeVirialRadius(
		vtkMultiProcessController* controller, vtkPointLocator* locator,
		vtkstd::string massArrayName, double softening,double overdensity,
		double maxR,double center[]);

	// Description:
	// From dataSet, initializes a newDataSet with point and cell 
	// arrays (empty), and corresponding  data arrays (empty) with 
	// identical names, and number of components to dataSet. Then copys, if
	// listed in the vtkIdList, from the old data set to the new data set the
	//  points and their cell data.
	static vtkPointSet* CopyPointsAndData(vtkPointSet* dataSet, vtkIdList*
	 	pointsInRadius);

	// Description:
	// Given a populated virialradiusinfo struct, returns a dataset
	// corresponding to only those points within the virial radius.
	// This method only works if input was vtkPolyData...
	static vtkPointSet* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo);

	// Description:
	// helper function to calculate the center based upon the source.
	// either the point, or the midpoint of a line
	static double* CalculateCenter(vtkDataSet* source);

	// Description
	// Given an input data set, the bin number, a list of points in the relevant
	// bin,  and the output table, computes the average radial velocity in the
	// bin, the radial velocity dispersion in the bin, the tangential velocity
	// in the bin, the tangentical velocity dispersion in the bin, and 
	// the angular momentum in the bin and adds it to the output table in 
	// the appropriate bin.
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
	static void ComputeStatisticsInBin(vtkPolyData* inputDataSet,
	 	double center[], int binNum, vtkIdList* pointsInBin, vtkTable* output);
	
	// Description:
	// given 3-vector vectorOne and 3-vector vectorTwo, computes the 
	// projection of vectorOnew onto vector two.
	// Can be used to calculate e.g.
	// radialVelocity = ComputeProjection(v,r);
	//  or velocitySquared = ComputeProjection(v,v);
	static double* ComputeProjection(double  vectorOne[],double vectorTwo[]);

	// Description:
	// Computes the vector difference between two 3-vectors
	// Can be used to calculate e.g. tangential velocity
	// tangentialVelocity = PointVectorDifference(v,radialVelocity);
	static double* PointVectorDifference(double vectorOne[], 
		double vectorTwo[]);

	// Description
	// Give a 3 vector v and a three vector r computes the specific angular 
	// momentum = r x v
	static double* ComputeAngularMomentum(double v[], double r[]);

	// Description:
	// Multiplies in place a 3-vector by a constant
	static void VecMultConstant(double vector[],double constant);

	// Description
	// Given a vSquaredAve and a vAve calculates the velocity dispersion
	// placing it in the output vector velocityDispersion
	static double* ComputeVelocityDispersion(vtkVariant vSquaredAve, 
		vtkVariant vAve);
	
	// Description:
	// helper function to compute radial velocity
	static double* ComputeRadialVelocity(double v[],double r[]);

	// Description:
	// helper function to compute tangential velocity
	static double* ComputeTangentialVelocity(double v[],double r[]);

	// Description:
	// helper function to compute velocity squared
	static double* ComputeVelocitySquared(double v[],double r[]);

	// Description:
	// helper function to compute radial velocity squared
	static double* ComputeRadialVelocitySquared(double v[],double r[]);

	// Description:
	// helper function to compute tangential velocity squared
	static double* ComputeTangentialVelocitySquared(double v[],double r[]);

	// Description:
	// helper function to compute circular velocity
	static double* ComputeCircularVelocity(vtkVariant cumulativeMass, 
		vtkVariant binRadius);

	// Description:
	// helper function to compute density
	static double* ComputeDensity(vtkVariant cumulativeMass, 
		vtkVariant binRadius);
	
	// Description:
	// Helper function to compute the midpoint between two points
	static double* ComputeMidpoint(double pointOne[], double pointTwo[]);
	
	// Don't wrap methods which use templates in python
	//BTX
	// Description:
	// shifts every item in array one to left (the first element is thrown away)
	// then sets inserts updateValue in the last, free slot
	static template <class T> void shiftLeftUpdate(T* array,
		int size, T updateValue);
private:
	vtkAstroVizHelpers(const vtkAstroVizHelpers&);  // Not implemented.
  void operator=(const vtkAstroVizHelpers&);  // Not implemented.
	//ETX
};
#endif
