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
#include <cmath>

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

	// If we should cutoff at the virial radius, then calculating where
	// the cutoff is and then redefining the data set to be only the points
	// up to this cutoff
	if(this->CutOffAtVirialRadius)
		{
		// calculating the bounds of this pointset
		double bounds[6];
		double upperBound[3];
		double lowerBound[3];
		input->GetPoints()->ComputeBounds();
		input->GetPoints()->GetBounds(bounds);
		upperBound[0]=bounds[0];
		upperBound[1]=bounds[2];
		upperBound[2]=bounds[4];
		lowerBound[0]=bounds[1];
		lowerBound[1]=bounds[3];
		lowerBound[2]=bounds[5];
		double maxR = sqrt(vtkMath::Distance2BetweenPoints(upperBound,\
																											 this->Center));
		double minR = sqrt(vtkMath::Distance2BetweenPoints(lowerBound,\
																											 this->Center));
		// Building the point locator and the struct to use as an 
		// input to the rootfinder.
		// 1. Building the point locator
		vtkSmartPointer<vtkPointLocator> locator = \
		 																	vtkSmartPointer<vtkPointLocator>::New();
		locator->SetDataSet(input);
		locator->BuildLocator();
		// 2. Building the struct
		LocatorInfo locatorInfo;
		locatorInfo.locator=locator;
		// copies the contents of this->Center to locatorInfo's arg center
		for(int i = 0; i < 3; ++i)
		{
			locatorInfo.center[i]=this->Center[i];
		}
		locatorInfo.criticalDensity=this->Delta;
		// but IllinoisRootFinder takes in a void pointer
		void* pntrLocatorInfo = &locatorInfo;
		// Now we are ready to run the root finder
		int numIter=0;
		double virialRadius=IllinoisRootFinder(OverDensityInSphere,\
																					pntrLocatorInfo,\
																					maxR,minR,
																					0.0,0.0,
																				  &numIter);
		vtkDebugMacro("virial radius is " << virialRadius);
		// Great, now we build a new dataset consisting only of points
		// which are within the virial radius. note that if there was an error
		// finding the virialRadius the radius returned is < 0
		if(virialRadius>0)
			{
				
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
						static_cast<float>(sqrt(vtkMath::Distance2BetweenPoints(nextPoint,\
																												this->Center)));
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

	// intializing the cumulative number array
	vtkSmartPointer<vtkFloatArray> cumulativeNumberArray=\
																		vtkSmartPointer<vtkFloatArray>::New();
		cumulativeNumberArray->SetName("cumulative number");
		cumulativeNumberArray->SetNumberOfComponents(1);
		cumulativeNumberArray->SetNumberOfTuples(output->GetNumberOfRows());
	output->AddColumn(cumulativeNumberArray);
	
	// we can only compute circular velocity and density if we bin by radius
	if(this->BinByRadius)
		{
		// intializing the circular velocity array
		vtkSmartPointer<vtkFloatArray> circularVelocityArray=\
																			vtkSmartPointer<vtkFloatArray>::New();
			circularVelocityArray->SetName("circular velocity");
			circularVelocityArray->SetNumberOfComponents(1);
			circularVelocityArray->SetNumberOfTuples(output->GetNumberOfRows());
		output->AddColumn(circularVelocityArray);

		// initializing the density in bin array
		vtkSmartPointer<vtkFloatArray> densityInBinArray=\
																			vtkSmartPointer<vtkFloatArray>::New();
			densityInBinArray->SetName("density");
			densityInBinArray->SetNumberOfComponents(1);
			densityInBinArray->SetNumberOfTuples(output->GetNumberOfRows());
		output->AddColumn(densityInBinArray);
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
		
		// we can only compute circular velocity and density if we bin by radius
		if(this->BinByRadius)
			{
			// we will need the bin radius for both calculations
			float binRadius = \
							output->GetValueByName(rowId,"radii from center").ToFloat();
			// computing the circular velocity (sqrt(M(<r)/r))
			float binCircularVelocity = sqrt(totalMass/binRadius);
			output->SetValueByName(rowId,"circular velocity",binCircularVelocity);
	
			// computing the density M(<r) / (4/3 pi r^3)
			float binDensity=totalMass/(4./3*M_PI*pow(binRadius,3));
			output->SetValueByName(rowId,"density",binDensity);
			}
		}	
	return 1;
}









