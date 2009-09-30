#include "PolyDataHelpers.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkSmartPointer.h"
//----------------------------------------------------------------------------
void AllocateDataArray(vtkPointSet* output, const char* arrayName, int numComponents, int numTuples)
{
	vtkSmartPointer<vtkFloatArray> dataArray=vtkSmartPointer<vtkFloatArray>::New();
  	dataArray->SetNumberOfComponents(numComponents);
  	dataArray->SetName(arrayName);
		dataArray->SetNumberOfTuples(numTuples);
  output->GetPointData()->AddArray(dataArray);
}
//----------------------------------------------------------------------------
vtkIdType SetPointValue(vtkPolyData* output,float pos[3])
{
	vtkIdType id=output->GetPoints()->InsertNextPoint(pos);
	output->GetVerts()->InsertNextCell(1, &id);
	return id;
}

//----------------------------------------------------------------------------
void SetDataValue(vtkPolyData* output, const char* arrayName, vtkIdType& id,float* data)
{
	output->GetPointData()->GetArray(arrayName)->SetTuple(id,data);
}

void GetDataValue(vtkPolyData* output, const char* arrayName, vtkIdType& id,double* data)
{
	output->GetPointData()->GetArray(arrayName)->GetTuple(id,data);
}
