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
#include "vtkFloatArray.h"
// Description:
// sets the point vertices in the output vector, assigning the point 
// a unique id and places one point per cell
vtkIdType SetPointValue(vtkPolyData* output,float pos[]);
// Description:
// sets the data value in the output vector in array arrayName at 
// position id to data.
void SetDataValue(vtkPointSet* output, const char* arrayName,\
			vtkIdType id,float data[]);
void SetDataValue(vtkPointSet* output, const char* arrayName,vtkIdType id,double data[]);
// Description:
// points the double* data to the value in the output vector in array arrayName at position id
// would like to be float, but the version of vtk I am working with does not support 
double* GetDataValue(vtkPointSet* output, const char* arrayName,\
 					vtkIdType id);
// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateFloatArray("density",1,100) creates a array of 100 
// scalar float densities
// AllocateFloatArray("velocity",3,100) creates a array of 100 vector 
// float velocities
void AllocateDataArray(vtkPointSet* output, const char* arrayName,\
 			int numComponents, int numTuples);
// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
// scalar double densities
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,\
	int numComponents, int numTuples);

// Description:
// returns a pointer to the point's coordinates in output which corresponds 
// to this id
double* GetPoint(vtkPointSet* output,vtkIdType id);
// Description:
// Computes the radial distance between two two dimensional points
double ComputeRadialDistance(double pointOne[],double pointTwo[]);



