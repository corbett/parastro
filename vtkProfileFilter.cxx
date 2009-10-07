/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkProfileFilter.cxx,v $
=========================================================================*/
#include "vtkProfileFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "vtkCellData.h"
#include "vtkTable.h"
#include "vtkSortDataArray.h"

vtkCxxRevisionMacro(vtkProfileFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkProfileFilter);

//----------------------------------------------------------------------------
vtkProfileFilter::vtkProfileFilter()
{
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
}

//----------------------------------------------------------------------------
void vtkProfileFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkProfileFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}




//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkTable* output = this->GetOutput();

	vtkSmartPointer<vtkFloatArray> XArray=vtkSmartPointer<vtkFloatArray>::New();
		XArray->DeepCopy(input->GetPointData()->GetArray("mass"));
		XArray->SetName("mass");
		vtkSortDataArray::Sort(XArray);
	vtkSmartPointer<vtkIntArray> dataValues=vtkSmartPointer<vtkIntArray>::New();
  	dataValues->SetNumberOfComponents(1);
		dataValues->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());	
		dataValues->SetName("N(particles>=mass)");		
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		dataValues->InsertValue(nextPointId,nextPointId);
		}
	// Updating the output
	output->AddColumn(XArray);
	output->AddColumn(dataValues);
	return 1;
}
