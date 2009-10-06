/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMassFunctionFilter.cxx,v $
=========================================================================*/
#include "vtkMassFunctionFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "vtkCellData.h"
#include <cmath>

vtkCxxRevisionMacro(vtkMassFunctionFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkMassFunctionFilter);

//----------------------------------------------------------------------------
vtkMassFunctionFilter::vtkMassFunctionFilter()
{
}

//----------------------------------------------------------------------------
vtkMassFunctionFilter::~vtkMassFunctionFilter()
{
}

//----------------------------------------------------------------------------
void vtkMassFunctionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkMassFunctionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}




//----------------------------------------------------------------------------
int vtkMassFunctionFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkRectilinearGrid* output = vtkRectilinearGrid::GetData(outputVector);
	// Setting the dimensions of this to be equal to our number of points,
	// in X, and equal to 1 in the Y and Z directions,
	// as these are dummy arrays.
		output->SetDimensions(input->GetPoints()->GetNumberOfPoints(),1,1);
		output->SetWholeExtent(0,input->GetPoints()->GetNumberOfPoints(),\
													 0,0,\
													 0,0);
	// Allocate the arrays for the X,Y, and Z coordinates, and inserts a 
	// single value for the dummy arrays, as well as allocate the array  
	// for the scalar data
	vtkSmartPointer<vtkDoubleArray> XArray=vtkSmartPointer<vtkDoubleArray>::New();
  	XArray->SetNumberOfComponents(1);
		XArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray> DummyYArray=vtkSmartPointer<vtkDoubleArray>::New();
  	DummyYArray->SetNumberOfComponents(1);
		DummyYArray->SetNumberOfTuples(1);
		DummyYArray->InsertValue(0,0.0);
	vtkSmartPointer<vtkDoubleArray> DummyZArray=vtkSmartPointer<vtkDoubleArray>::New();
  	DummyZArray->SetNumberOfComponents(1);	
		DummyZArray->SetNumberOfTuples(1);
		DummyZArray->InsertValue(0,0.0);
	vtkSmartPointer<vtkIntArray> dataValues=vtkSmartPointer<vtkIntArray>::New();
  	dataValues->SetNumberOfComponents(1);
		dataValues->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());	
		dataValues->SetName("data values");		
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		XArray->InsertValue(nextPointId,nextPointId);
		dataValues->InsertValue(nextPointId,pow(nextPointId,2));
		}
	// Updating the output
	output->SetXCoordinates(XArray);
	output->SetYCoordinates(DummyYArray);
	output->SetZCoordinates(DummyZArray);
	output->GetPointData()->SetScalars(dataValues);
	return 1;
}
