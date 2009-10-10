#include "DataSetHelpers.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include <cmath>
//----------------------------------------------------------------------------
void AllocateDataArray(vtkPointSet* output, const char* arrayName,\
 			int numComponents, int numTuples)
{
	vtkSmartPointer<vtkFloatArray> dataArray=\
		vtkSmartPointer<vtkFloatArray>::New();
  	dataArray->SetNumberOfComponents(numComponents);
  	dataArray->SetName(arrayName);
		dataArray->SetNumberOfTuples(numTuples);
  output->GetPointData()->AddArray(dataArray);
}
//----------------------------------------------------------------------------
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,\
 			int numComponents, int numTuples)
{
	vtkSmartPointer<vtkDoubleArray> dataArray=\
		vtkSmartPointer<vtkDoubleArray>::New();
  	dataArray->SetNumberOfComponents(numComponents);
  	dataArray->SetName(arrayName);
		dataArray->SetNumberOfTuples(numTuples);
  output->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
vtkIdType SetPointValue(vtkPolyData* output,float pos[])
{
	vtkIdType id=output->GetPoints()->InsertNextPoint(pos);
	output->GetVerts()->InsertNextCell(1, &id);
	return id;
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
//----------------------------------------------------------------------------

