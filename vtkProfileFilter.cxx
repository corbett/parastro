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
	// TODO: doesn't actually initialize like this, but shorthand for now
	// for what I have in mind
	this->AdditionalProfileQuantities = \
		{"number in bin","radii from center","cumulative mass","cumulative number","circular velocity","density","radial velocity","radial velocity dispersion","tangential velocity","tangential velocity dispersion","angular momentum"}; // ALWAYS need at least "number in bin"
	// TODO: doesn't actually initialize like this, but shorthand for now
	// for what I have in mind
	this->CumulativeQuantities = \
		{"mass","number in bin"};
	
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
	this->Center=source->GetPoint(0);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::GenerateProfile(vtkPolyData* input,vtkTable* output)
{	
	this->InitializeBins(input,output);
	this->ComputeStatistics(input,output);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::InitializeBins(vtkPolyData input,
	vtkTable* output)
{
	this->CalculateAndSetBinExtents(input);
	vtkSmartPointer<vtkDataArray> nextArray;
	for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
		{
		nextArray = input->GetPointData()->GetArray(i);
		vtkstd::string totalName = nextArray->GetName() + "_total";
		// Allocating an column for the total sum of the existing quantities
		AllocateDataArray(output,totalName.c_str(),
			nextArray->GetNumberOfComponents(),this->BinNumber);
		// Allocating an column for the averages of the existing quantities
		vtkstd::string averageName = nextArray->GetName() + "_average";
		AllocateDataArray(output,averageName.c_str(),
			nextArray->GetNumberOfComponents(),this->BinNumber);
		if(this->CumulativeQuantities->LookupValue(nextArray->GetName())>=0)
			{
			// we should also consider this a cumulative quantity
			vtkstd::string cumulativeName = nextArray->GetName() + "_cumulative";
			AllocateDataArray(output,cumulativeName.c_str(),1,this->BinNumber);
			}
		}
	// For our additional quantities, allocating a column for the average 
	// and the sum. These are restricted to be scalars
	for(int i = 0; 
		i < this->AdditionalProfileQuantities.GetNumberOfValues();
	 	++i)
		{
		vtkstd::string nextName=this->AdditionalProfileQuantities.GetValue(i);
		vtkstd::string totalName = nextName + "_total";
		// Allocating an column for the total sum of the existing quantities
		AllocateDataArray(output,totalName.c_str(),1,this->BinNumber);
		// Allocating an column for the averages of the existing quantities
		vtkstd::string averageName = nextName + "_average";
		AllocateDataArray(output,averageName.c_str(),1,this->BinNumber);
		if(this->CumulativeQuantities->LookupValue(nextName)>=0)
			{
			// we should also consider this a cumulative quantity
			vtkstd::string cumulativeName = nextName + "_cumulative";
			AllocateDataArray(output,cumulativeName.c_str(),1,this->BinNumber);
			}
		}
}


//----------------------------------------------------------------------------
void vtkProfileFilter::ComputeStatistics(vtkPolyData* inputDataSet,vtkTable* output)
{
	double* nextPoint = new double[3]; // need for the loop
	for(int nextPointId = 0;
	 		nextPointId < dataSet->GetPoints()->GetNumberOfPoints();
	 		++nextPointId)
		{
			nextPoint = GetPoint(dataSet,nextPointId);
			UpdateBinStatistics(inputDataSet,nextPoint,nextPointId,output);
		}
	// Updating averages
	UpdateBinAverages(this->BinNumber,output);
	//Finally some memory management
	delete [] nextPoint;
}

//----------------------------------------------------------------------------
void 	UpdateBinAverages(vtkTable* output)
{
	/*
	for(bin in this->BinNumber)
	{
		double binSize=
		//done with for loop divide everything by N
		VecMultConstant(vAve,1./binSize);
		VecMultConstant(vRadAve,1./binSize);
		VecMultConstant(vTanAve,1./binSize);
		VecMultConstant(jAve,1./binSize);	
		VecMultConstant(vSquaredAve,1./binSize);	
		VecMultConstant(vRadSquaredAve,1./binSize);		
		VecMultConstant(vTanSquaredAve,1./binSize);		
		// calculate velocity dispersions, taking the necessary square roots	
		ComputeVelocityDispersion(vSquaredAve,vAve,vDisp);
		ComputeVelocityDispersion(vRadSquaredAve,vRadAve,vRadDisp);
		ComputeVelocityDispersion(vTanSquaredAve,vTanAve,vTanDisp);
		// add vAve, vRadAve, vTanAve,, vAveDisp, vRadAveDisp, vTanAveDisp and j
		// to the ouput table.
	
		//
	/*
		 finally getting and computing
			// computing the circular velocity (sqrt(M(<r)/r))
			float binCircularVelocity = sqrt(totalMass/binRadius);
			output->SetValueByName(rowId,"circular velocity",binCircularVelocity);

			// computing the density M(<r)/(4/3 pi r^3)
			float binDensity=totalMass/(4./3*M_PI*pow(binRadius,3));
			output->SetValueByName(rowId,"density",binDensity);
			// computing the radial velocity
			// computing the tangential velocity
			// computing the specific angular momentum
			
			}
	*/
}



//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBin(int binNum, int totalBins, 
	BinUpdateType updateType, char* attributeName, int attributeNumComponents,
	double* dataToAdd[], vtkTable* output)
{
	double* data=output->GetValueByName(binNum,attributeName).ToDouble();
	for(int i = 0; i < numComponents; ++i)
		{
			switch(updateType)
			{
			case add:
				data[i]+=dataToAdd[i];
				break;
			case multiply:
				data[i]*=dataToAdd[i];
				break;
			}
		}
	output.SetValueByName(binNum,attributeName,data);
}

void vtkProfileFilter::UpdateCumulativeBins(int binNum, int totalBins,
	BinUpdateType updateType, char* attributeName, int attributeNumComponents,
	double* dataToAdd[], vtkTable* output)
{
		for(int bin = binNum; bin < totalBins; ++bin)
		{
			UpdateBin(bin,totalBins,updateType,attributeName,\
				attributeNumComponents,dataToAdd[],output);
		}
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBinStatistics(vtkPolyData* inputDataSet,
 	double x[], vtkIdType pointGlobalId,vtkTable* output)
{
	//TODO: add later cumulative capabilities.
	// As we bin by radius allways need this
	double* r=PointVectorDifference(x,center);
	// Many of the quantities explicitely require the velocity
	double* v=GetDataValue(inputDataSet,"velocity",pointGlobalId);
	double binNum=GetBinNumber(r,this->BinNumber,this->BinSpacing);
	
	for(int i = 0; i < dataSet->GetPointData()->GetNumberOfArrays(); ++i)
		{
		vtkSmartPointer<vtkDataArray> nextArray = \
		 	input->GetPointData()->GetArray(i);
		// getting the data for this point
		double* nextData = GetDataValue(inputDataSet,
			nextArray->GetName(),pointGlobalId);
		// Updating the total bin
		vtkstd::string totalName = nextArray->GetName() + "_total";		
		UpdateBin(binNum, totalBins, BinUpdateType.add,totalName,
			nextArray->GetNumberOfComponents(), nextData, output);
		if(this->CumulativeQuantities->LookupValue(nextArray->GetName())>=0)
			{
			// we should also consider this a cumulative quantity
			vtkstd::string cumulativeName = nextArray->GetName() + "_cumulative";
			UpdateCumulativeBins(binNum,totalBins,BinUpdateType.add,
				cumulativeName,1,additionalData,output);
			}
		}
	// For our additional quantities, allocating a column for the average 
	// and the sum. These are restricted to be scalars
	for(int i = 0; i < 
		this->AdditionalProfileQuantities.GetNumberOfValues(); 
		++i)
		{
		vtkstd::string nextName=this->AdditionalProfileQuantities.GetValue(i);
		double* additionalData = CalculateAdditionalProfileQuantity(nextName,v,r);
		// updating the totalbin
		vtkstd::string totalName = nextName + "_total";
		UpdateBin(binNum,totalBins,BinUpdateType.add,totalName,\
			1,additionalData,output);
		if(this->CumulativeQuantities->LookupValue(nextName)>=0)
			{
			// we should also consider this a cumulative quantity
			vtkstd::string cumulativeName = nextName + "_cumulative";
			UpdateCumulativeBins(binNum,totalBins,BinUpdateType.add,
				cumulativeName,1,additionalData,output);
			}
		}
}

//----------------------------------------------------------------------------
double* vtkProfileFilter::CalculateAdditionalProfileQuantity(
	vtkstd::string nextName)
{
	// Some inefficiency by recomputing quantities, but paid for with 
	// flexibility, i.e don't have to call in a certain order or can compute
	// one quantity without storing the other
	if(nextName == "radial velocity")
		{
		return ComputeRadialVelocity(v,r);
		}
	else if(nextName == "tangential velocity")
		{
		return ComputeTangentialVelocity(v,r);
		}
	else if(nextName == "angular momentum")
		{
		return ComputeAngularMomentum(v,r);
		}
	else if(nextName == "velocity dispersion")
		{
		return ComputeVelocitySquared(v,r);
		}
	else if(nextName == "radial velocity dispersion")
		{
		return ComputeRadialVelocitySquared(v,r);
		}
	else if(nextName == "tangential velocity dispersion")
		{
		return ComputeTangentialVelocitySquared(v,r);
		}
	else if(nextName=="number in bin")
		{
			return {1.0};
		}
	else
		{
		vtkWarningMacro("input arrray requested not found, quantity returned as \
			array of zero");
		return {0.0};
		}
}

//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation *request,
																	vtkInformationVector **inputVector,
																	vtkInformationVector *outputVector)
{
	// Now we can get the dataSet with which we want to work
 	vtkPolyData* dataSet = vtkPolyData::GetData(inputVector[0]);
	// Setting the center based upon the selection in the GUI
	vtkDataSet* pointInfo = vtkDataSet::GetData(inputVector[1]);
	vtkTable* output = vtkTable::GetData(outputVector[0]);
	this->CalculateAndSetCenter(pointInfo);
	// If we want to cut off at the virial radius, compute this, and remove the
	// portion of the data set we don't care about
	if(this->CutOffAtVirialRadius)
		{
		VirialRadiusInfo virialRadiusInfo = \
		 	ComputeVirialRadius(dataSet,this->Delta,this->Center);
		vtkErrorMacro("virial radius is " << virialRadiusInfo.virialRadius);
		// note that if there was an error finding the virialRadius the 
		// radius returned is < 0
		//setting the input to this newInput
		if(virialRadiusInfo.virialRadius>0)
			{
			dataSet = \
				GetDatasetWithinVirialRadius(virialRadiusInfo);	
			this->GenerateProfile(dataSet,output);
			delete [] dataSet;
			return 1;
			}
		else
			{
			vtkErrorMacro("Something has gone wrong with the virial radius finding. Perhaps change your delta, or your center, or if you are truely puzzled check out ProfileHelpers.cxx. For now binning out to the max radius instead of the virial.");
			}
		}
	this->GenerateProfile(dataSet,output);
	return 1;
}









