#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
// Description:
// sets the point vertices in the output vector, assigning the point a unique id
// and places one point per cell
vtkIdType SetPointValue(vtkPolyData* output,float pos[3]);
// Description:
// sets the data value in the output vector in array arrayName at position id to data.
void SetDataValue(vtkPolyData* output, const char* arrayName, vtkIdType& id,float* data);
// Description:
// points the double* data to the value in the output vector in array arrayName at position id
// would like to be float, but the version of vtk I am working with does not support 
void GetDataValue(vtkPolyData* output, const char* arrayName, vtkIdType& id,double* data);
// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateFloatArray("density",1,100) creates a array of 100 scalar float densities
// AllocateFloatArray("velocity",3,100) creates a array of 100 vector float velocities
void AllocateDataArray(vtkPolyData* output, const char* arrayName, int numComponents, int numTuples);
