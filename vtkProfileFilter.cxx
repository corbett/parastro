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
#include "astrovizhelpers/ProfileHelpers.h"
#include "vtkCellData.h"
#include "vtkTable.h"
#include "vtkSortDataArray.h"
#include "vtkMath.h"
#include "vtkInformationDataObjectKey.h"
#include <cmath>

vtkCxxRevisionMacro(vtkProfileFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkProfileFilter);

//----------------------------------------------------------------------------
vtkProfileFilter::vtkProfileFilter():vtkExtractHistogram()
{
	this->SetCalculateAverages(1); // no longer taking this in as an option 
	 															// may later actually disable
  this->SetNumberOfInputPorts(2);
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
void vtkProfileFilter::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::CalculateAndSetCenter(vtkDataSet* source)
{
	//TODO: this can later be done as in the XML documentation for this filter; 	  
	// for now, only getting the first point. this is the point selected in the
	// GUI, or the first end of the line selected in the GUI
	double* selectedCenter=source->GetPoint(0);
	this->SetCenter(selectedCenter);
	delete [] selectedCenter;

	/*
  delete this->PointList;
  delete this->CellList;

  this->PointList = new vtkDataSetAttributes::FieldList(1);
  this->PointList->InitializeFieldList(source->GetPointData());

  this->CellList = new vtkDataSetAttributes::FieldList(1);
  this->CellList->InitializeFieldList(source->GetCellData());
	*/
}

//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation *request,
																	vtkInformationVector **inputVector,
																	vtkInformationVector *outputVector)
{
 	vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);
	vtkInformationVector* newInputVector = vtkInformationVector::New();
	newInputVector->Copy(*inputVector);
	/*
	vtkPolyData* newDataSet = vtkPolyData::New();
	newDataSet->DeepCopy(oldInput);
	cout << "number of points "<< newDataSet->GetNumberOfPoints() << "\n\n";
	double* nextData = GetDataValue(newDataSet,
		"potential",1000);
	cout << "potential is "<< nextData[0];
	*/
/*
	inputVector[0]->GetInformationObject(0)->Set(
		vtkDataObject::DATA_OBJECT(),newDataSet);

	vtkDataSet* pointInfo = vtkDataSet::GetData(inputVector[1]);
	this->CalculateAndSetCenter(pointInfo);
*/
//for testing only 
	//TODO: change
	// temporary hack swapping inputs so that superclass doesn't get confused
	// If we should cutoff at the virial radius, then calculating where
	// the cutoff is and then redefining the data set to be only the points
	// up to this cutoff
	if(this->CutOffAtVirialRadius)
		{
		VirialRadiusInfo virialRadiusInfo = \
		 	ComputeVirialRadius(input,this->Delta,this->Center);
		vtkErrorMacro("virial radius is " << virialRadiusInfo.virialRadius);
		// note that if there was an error finding the virialRadius the 
		// radius returned is < 0
		//setting the input to this newInput
		if(virialRadiusInfo.virialRadius>0)
			{
//			vtkPolyData* newDataSet = \
				GetDatasetWithinVirialRadius(virialRadiusInfo);	// note vtkPolyData 
																						// must now manually be deleted!
			//setting the input to this newInput
			
//			inputVector[0]->GetInformationObject(0)->Set(
//				vtkDataObject::DATA_OBJECT(),newDataSet);
			// TODO: debug this--- GetDatasetWithinVirialRadius is currently
			// one point only, moreover even if it's all 1001 points, the new 		
			// histogram doesn't find the input data arrays
			}
		else
			{
			vtkErrorMacro("Something has gone wrong with the virial radius finding. Perhaps change your delta, or your center, or if you are truely puzzled check out ProfileHelpers.cxx. For now binning out to the max radius instead of the virial.");
			}
		}

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
		AllocateDataArray(input,
			"radii from center",1,input->GetPoints()->GetNumberOfPoints());
		double* nextPoint = new double[3]; // need for the loop
		for(int nextPointId = 0;\
		 		nextPointId < input->GetPoints()->GetNumberOfPoints();
		 		++nextPointId)
			{
				nextPoint = GetPoint(input,nextPointId);
				float* radius=new float[1];
				radius[0] = \
					static_cast<float>(
					sqrt(vtkMath::Distance2BetweenPoints(nextPoint,
					this->Center)));
				SetDataValue(input,"radii from center",nextPointId,radius);
				// and finally some memory management
				delete [] radius;
			}
		// setting the input array to process to be radii
		this->SetInputArrayToProcess(0,0,0,
			vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
			"radii from center");
		// and finally finally some memory management
		delete [] nextPoint;
		}
	// Just calling the superclass for now
	//TODO: changing to newInputVector, for test
	vtkExtractHistogram::RequestData(request,&newInputVector,outputVector);
	// Getting the output data
	vtkTable* const output = vtkTable::GetData(outputVector, 0);
	// Modifying the output from vtkExtractHistogram
	AllocateDataArray(output,
		"cumulative mass",1,output->GetNumberOfRows());
	AllocateDataArray(output,
		"cumulative number",1,output->GetNumberOfRows());
		// we can only compute circular velocity, vtan
		// vrad (and their associated dispersions), angular momentum
		// and density, if we bin by radius
	if(this->BinByRadius)
		{
		AllocateDataArray(output,
			"circular velocity",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"density",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"radial velocity",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"radial velocity dispersion",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"tangential velocity",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"tangential velocity dispersion",1,output->GetNumberOfRows());
		AllocateDataArray(output,
			"angular momentum",1,output->GetNumberOfRows());
		}

	float totalMass=0;
	float totalNumber=0;
	for(int rowId = 0; rowId < output->GetNumberOfRows(); ++rowId)
		{
		// computing the cumulative mass
		float binMassTotal = \
			output->GetValueByName(rowId,"mass_total").ToFloat();
		totalMass+=binMassTotal;
		output->SetValueByName(rowId,"cumulative mass",totalMass);

		// computing the cumulative number
		float binNumberTotal = \
		 	output->GetValueByName(rowId,"bin_values").ToFloat();
		totalNumber+=binNumberTotal;
		output->SetValueByName(rowId,"cumulative number",totalNumber);
		
		// we can only compute circular velocity, vtan
		// vrad (and their associated dispersions), angular momentum
		// and density, if we bin by radius
		if(this->BinByRadius)
			{
			// we will need the bin radius for both calculations
			float binRadius = \
				output->GetValueByName(rowId,"radii from center").ToFloat();
			// computing the circular velocity (sqrt(M(<r)/r))
			float binCircularVelocity = sqrt(totalMass/binRadius);
			output->SetValueByName(rowId,"circular velocity",binCircularVelocity);
	
			// computing the density M(<r)/(4/3 pi r^3)
			float binDensity=totalMass/(4./3*M_PI*pow(binRadius,3));
			output->SetValueByName(rowId,"density",binDensity);
			// computing the radial velocity
// TODO also keep track of radial vector
//			double* velocity = output->GetValueByName(rowId,"density");
//			double* radius	= output->GetValueByName(rowId,"radius");
//			double radialVelocity = Dot(velocity,radius)/Norm(radius);
//			output->SetValueByName(rowId,"radial velocity",radialVelocity);
			// computing the tangential velocity
			// computing the specific angular momentum
			}
		}	
	return 1;
}









