#include "DataSetHelpers.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkSphereSource.h"
#include "vtkSmartPointer.h"
#include <cmath>
/*----------------------------------------------------------------------------
*
* Work with vtkDataArray
*
*---------------------------------------------------------------------------*/
void InitializeDataArray(vtkDataArray* dataArray, const char* arrayName,\
 			int numComponents, int numTuples)
{
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	dataArray->SetNumberOfTuples(numTuples);
}

/*----------------------------------------------------------------------------
*
* Work with vtkTable
*
*---------------------------------------------------------------------------*/
void AllocateDataArray(vtkTable* output, const char* arrayName,\
 			int numComponents, int numTuples)
{
	vtkSmartPointer<vtkFloatArray> dataArray=\
		vtkSmartPointer<vtkFloatArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);		
  output->AddColumn(dataArray);
}
/*----------------------------------------------------------------------------
*
* Work with VtkPolyData (a derived class from vtkPointSet)
*
*---------------------------------------------------------------------------*/

//----------------------------------------------------------------------------
vtkIdType SetPointValue(vtkPolyData* output,float pos[])
{
	vtkIdType id=output->GetPoints()->InsertNextPoint(pos);
	output->GetVerts()->InsertNextCell(1, &id);
	return id;
}

//----------------------------------------------------------------------------
float* DoublePointToFloat(double point[])
{
	float* floatPoint = new float[3];
	for(int i = 0; i < 3; ++i)
	{
		floatPoint[i]=static_cast<float>(point[i]);
	}
	return floatPoint;
}

//----------------------------------------------------------------------------
void CreateSphere(vtkPolyData* output,double radius,double center[])
{
	vtkSmartPointer<vtkSphereSource> sphere = \
	 															vtkSmartPointer<vtkSphereSource>::New();
	sphere->SetRadius(radius);
	sphere->SetCenter(center);
	sphere->Update();
	//Setting the points in the output to be those of the sphere
	output->SetPoints(sphere->GetOutput()->GetPoints());
	output->SetVerts(sphere->GetOutput()->GetVerts());
	output->SetPolys(sphere->GetOutput()->GetPolys());
}

/*----------------------------------------------------------------------------
*
* Work with VtkPointSet
*
*---------------------------------------------------------------------------*/

//----------------------------------------------------------------------------
void AllocateDataArray(vtkPointSet* output, const char* arrayName,\
 			int numComponents, int numTuples)
{
	vtkSmartPointer<vtkFloatArray> dataArray=\
		vtkSmartPointer<vtkFloatArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,\
 			int numComponents, int numTuples)
{
	vtkSmartPointer<vtkDoubleArray> dataArray=\
		vtkSmartPointer<vtkDoubleArray>::New();
	InitializeDataArray(dataArray,arrayName,numComponents,numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
double* GetPoint(vtkPointSet* output,vtkIdType id)
{
	double* nextPoint=new double[3]; 
	output->GetPoints()->GetPoint(id,nextPoint);
	return nextPoint;
}

//----------------------------------------------------------------------------
void SetDataValue(vtkPointSet* output, const char* arrayName,\
			vtkIdType id,float data[])
{
	output->GetPointData()->GetArray(arrayName)->SetTuple(id,data);
}

//----------------------------------------------------------------------------
void SetDataValue(vtkPointSet* output, const char* arrayName,\
			vtkIdType id,double data[])
{
	output->GetPointData()->GetArray(arrayName)->SetTuple(id,data);
}

//----------------------------------------------------------------------------
double* GetDataValue(vtkPointSet* output, const char* arrayName,\
 					vtkIdType id)
{
	double* data=new double[3];
	output->GetPointData()->GetArray(arrayName)->GetTuple(id,data);
	return data;
}

/*----------------------------------------------------------------------------
*
* Work with vtkInformationVector and vtkInformation objets
*
*---------------------------------------------------------------------------*/
//----------------------------------------------------------------------------
vtkInformationVector** DeepCopyInputVector(vtkInformationVector** inputVector,
	int inputVectorSize)
{
	vtkInformationVector** newInputVector = \
		new vtkInformationVector *[inputVectorSize];
	for(int i = 0; i < inputVectorSize; ++i)
		{
		vtkInformationVector* newInput = vtkInformationVector::New();
		// performas a deep copy of inputVector
		newInput->Copy(inputVector[i],1);
		newInputVector[i]=newInput;
		}
	return newInputVector; // note caller must manage this memory
												// see helper function DeleteDeepCopyInput below
}

//----------------------------------------------------------------------------
