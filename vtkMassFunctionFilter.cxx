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
	// setting the dimensions of this to be equal to our number of points,
	// in X and Y each direction, and equal to 1 in the z direction,
	// as this is a dummy array
		output->SetDimensions(input->GetPoints()->GetNumberOfPoints(),\
						1,\
						1);
	// Allocate the arrays for the X,Y, and Z data
	vtkSmartPointer<vtkDoubleArray> XArray=\
		vtkSmartPointer<vtkDoubleArray>::New();
  	XArray->SetNumberOfComponents(1);
		XArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	vtkSmartPointer<vtkIntArray> DummyYArray=\
		vtkSmartPointer<vtkIntArray>::New();
  	DummyYArray->SetNumberOfComponents(1);
		DummyYArray->SetNumberOfTuples(1);
	vtkSmartPointer<vtkDoubleArray> DummyZArray=\
		vtkSmartPointer<vtkDoubleArray>::New();
  	DummyZArray->SetNumberOfComponents(1);	
		DummyZArray->SetNumberOfTuples(1);
	// if we first sorted the input points by mass (increasing)
	// then this would produce a scatter plot of the mass function
	// TODO: implement
	// right now just putting the architecture in place
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		// GetTuple function which operates with float
		double* mass=GetDataValue(input,"mass",nextPointId);
		XArray->InsertValue(nextPointId,mass[0]);
		//finally some memory management
		delete [] nextPoint;
		delete [] mass;
		}
	// Updating the output
	output->SetXCoordinates(XArray);
	output->SetYCoordinates(DummyYArray);
	output->SetZCoordinates(DummyZArray);  
	vtkWarningMacro("output has points: "\
									<< output->GetNumberOfPoints());
	vtkWarningMacro("output while the x coordinates: "\
									<< output->GetXCoordinates()->GetNumberOfTuples());
	
  return 1;
}
