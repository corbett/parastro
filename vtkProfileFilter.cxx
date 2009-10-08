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
vtkProfileFilter::vtkProfileFilter():vtkExtractHistogram()
{
	this->SetCalculateAverages(1); // no longer taking this in as an option
																 // may later actually disable
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
}

//----------------------------------------------------------------------------
void vtkProfileFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkExtractHistogram::PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation *request,\
																	vtkInformationVector **inputVector,\
																	vtkInformationVector *outputVector)
{
	// If we should bin by radius, first calculated add the radii to outputdata
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkTable* output = this->GetOutput();		
	//Just calling the superclass as a test
	vtkExtractHistogram::RequestData(request,inputVector,outputVector);
	/*
  // Get input and output data.
	// Allocate data structures
	// TODO: dummy allocation
	vtkSmartPointer<vtkFloatArray> XArray=vtkSmartPointer<vtkFloatArray>::New();
		XArray->DeepCopy(input->GetPointData()->GetArray("mass"));
		XArray->SetName("mass");
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
			double* nextPoint=GetPoint(input,nextPointId);
			int binNum=this->GetBinNum(nextPoint,\
																 binLowerBound,binUpperBound,binWidth);
			// Finally some memory management
			delete [] nextPoint;
		}
	// Updating the output
	// TODO: dummy output
	output->AddColumn(XArray);
	return 1;
	*/
}
