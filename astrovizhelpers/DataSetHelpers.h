#include "vtkPolyData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkDataSetAttributes.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkTable.h"
#include "vtkFloatArray.h"
#include "vtkInformationVector.h"
#include <iostream>
#include <sstream>

class vtkMultiProcessController;
/*
* The following methods take and modify vtkPolyData
*/
// Description:
// sets the point vertices in the output vector, assigning the point 
// a unique id and places one point per cell. Works only with PolyData
vtkIdType SetPointValue(vtkPolyData* output,float pos[]);

// Description:
// creates a sphere of radius R, center center in the output.
void CreateSphere(vtkPolyData* output,double radius,double center[]);

/*
* The following methods take and modify vtkPointSet data
*/
// Description:
// sets the point vertices in the output vector, assigning the point 
// a unique. 
vtkIdType SetPointValue(vtkPointSet* output,float pos[]);

// Description:
// sets the data value in the output vector in array arrayName at 
// position id to data.
void SetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id,float data[]);
void SetDataValue(vtkPointSet* output, const char* arrayName,vtkIdType id,double data[]);

// Description:
// points the double* data to the value in the output vector in 
// array arrayName at position id would like to be float, but the 
// version of vtk I am working with does not support 
double* GetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id);

// Description:
// create a vtkDataArray of floats with the  name arrayName, number /
// of components. place it in the vtkPointSet
void AllocateDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, int numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T. place it in
// the vtkTable
void AllocateDataArray(vtkTable* output, const char* arrayName,
	int numComponents, int numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
// scalar double densities
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, int numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
// scalar double densities
void AllocateIntDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, int numTuples);
// Description:
// Given a VTK array, sets it to have arrayName, numComponents and numTuples
// unlike the above methods, does not create the array, nor add it to the 
// output	
void InitializeDataArray(vtkDataArray* dataArray, const char* arrayName,
	int numComponents, int numTuples);

// Description:
// returns a pointer to the point's coordinates in output which corresponds 
// to this id
double* GetPoint(vtkPointSet* output,vtkIdType id);

// Description:
// takes in a double array of size three representing a point
// and converts it to an array of the same size but in float precision
float* DoublePointToFloat(double point[]);

// Description:
// Creates a new information vector that is a deep copy of the old one
// Note, returned vector of pointers to vtkInformationVectors must be 
// deleted by caller. See DeleteDeepCopyInput below.
vtkInformationVector** DeepCopyInputVector(
	vtkInformationVector** inputVector,int inputVectorSize); 

// Description: 
// Deletes an array of pointers to vtkInformationVectors
void DeleteDeepCopyInputVector(vtkInformationVector** inputVector, 
	int inputVectorSize);

// Description:
// converts a double to a vtkstd::string
inline vtkstd::string ToString(double x)
{
  std::ostringstream o;
	o << x;
  return o.str();
}
// Description:
// returns true if this process should be run in parallel
// (i.e. we have  non-null controller and more than one process to work with)
bool RunInParallel(vtkMultiProcessController* controller);




