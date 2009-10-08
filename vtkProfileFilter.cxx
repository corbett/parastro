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
	// If we should bin by radius, first calculate and add the radii 
	// to outputdata
	if(this->BinByRadius)
		{
		vtkErrorMacro("bin by radius=true");
		// getting the input and output
		// according to paraview's pipeline architecture
		// I am *not* supposed to modify the input
		// breaking the law here 
		// TODO: change the strategy. modify the inputVector instead
  	vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
		// intializing the radii array
		vtkSmartPointer<vtkFloatArray> radiiArray=\
																		vtkSmartPointer<vtkFloatArray>::New();
			radiiArray->SetName("radii from center");
			radiiArray->SetNumberOfComponents(1);
			radiiArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
		double* nextPoint; // need for the loop
		for(int nextPointId = 0;\
		 		nextPointId < input->GetPoints()->GetNumberOfPoints();
		 		++nextPointId)
			{
				nextPoint=GetPoint(input,nextPointId);
				float radius=\
						static_cast<float>(ComputeRadialDistance(nextPoint,this->Center));
				radiiArray->InsertValue(nextPointId,radius);
			}
		// finally adding the new radius vector to our output 
		input->GetPointData()->AddArray(radiiArray);
		// and finally finally some memory management
		// setting the input array to process to be radii
		this->SetInputArrayToProcess(0,0,0,
						vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
						"radii from center");
		delete [] nextPoint;
		}
	// Just calling the superclass for now
	vtkExtractHistogram::RequestData(request,inputVector,outputVector);
	return 1;
}
