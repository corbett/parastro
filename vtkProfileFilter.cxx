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
 	vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// If we should bin by radius, first calculate and add the radii 
	// to outputdata
	if(this->BinByRadius)
		{
		// getting the input and output
		// according to paraview's pipeline architecture
		// I am *not* supposed to modify the input
		// breaking the law here 
		// TODO: change the strategy. modify the inputVector instead
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
		// setting the input array to process to be radii
		this->SetInputArrayToProcess(0,0,0,
						vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
						"radii from center");
		// and finally finally some memory management
		delete [] nextPoint;
		}
	// Just calling the superclass for now
	vtkExtractHistogram::RequestData(request,inputVector,outputVector);
	// Getting the output data
	vtkTable* const output = vtkTable::GetData(outputVector, 0);
	// Modifying the output from vtkExtractHistogram
		// intializing the cumulative mass array
	vtkSmartPointer<vtkFloatArray> cumulativeMassArray=\
																	vtkSmartPointer<vtkFloatArray>::New();
		cumulativeMassArray->SetName("cumulative mass");
		cumulativeMassArray->SetNumberOfComponents(1);
		cumulativeMassArray->SetNumberOfTuples(output->GetNumberOfRows());
	output->AddColumn(cumulativeMassArray);

	// adding to the table, just to check. that this strategy works
	float totalMass=0;
	for(int rowId = 0; rowId < output->GetNumberOfColumns(); ++rowId)
		{
			// TODO: make this calculation correct.
			vtkVariant binMassTotal = output->GetValueByName(rowId,"mass_total");
			totalMass+=binMassTotal.ToFloat();
			vtkErrorMacro(" bin/row id: " << rowId \
			 							<< " total_mass of bin: " << binMassTotal.ToFloat()\
										<< " total mass: " << totalMass);
			output->SetValueByName(rowId,"cumulative mass",totalMass);
		}	
	return 1;
}









