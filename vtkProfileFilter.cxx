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
using vtkstd::string;

vtkCxxRevisionMacro(vtkProfileFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkProfileFilter);

//----------------------------------------------------------------------------
vtkProfileFilter::vtkProfileFilter()
{
  this->SetNumberOfInputPorts(2);
	// TODO: doesn't actually initialize like this, but shorthand for now
	// for what I have in mind
	this->AdditionalProfileQuantities = vtkStringArray::New();
//		{"number in bin","radii from center","cumulative mass","cumulative number","circular velocity","density","radial velocity","radial velocity dispersion","tangential velocity","tangential velocity dispersion","angular momentum"}; // ALWAYS need at least "number in bin"
	// TODO: doesn't actually initialize like this, but shorthand for now
	// for what I have in mind
	this->CumulativeQuantities = vtkStringArray::New();
//		{"mass","number in bin"};	
	this->MaxR=1.0;
	this->Delta=0.0;
	this->BinNumber=30;
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
	this->CumulativeQuantities->Delete();
	this->AdditionalProfileQuantities->Delete();
}

//----------------------------------------------------------------------------
void vtkProfileFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	// TODO: finish
  os << indent << "overdensity: "
     << this->Delta << "\n"
		 << indent << "bin number: "
     << this->BinNumber << "\n";
}

//----------------------------------------------------------------------------
void vtkProfileFilter::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkProfileFilter::FillInputPortInformation (int port, 
                                                   vtkInformation *info)
{
  this->Superclass::FillInputPortInformation(port, info);
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation *request,
																	vtkInformationVector **inputVector,
																	vtkInformationVector *outputVector)
{
	// Now we can get the input with which we want to work
 	vtkPolyData* dataSet = vtkPolyData::GetData(inputVector[0]);
	// Setting the center based upon the selection in the GUI
	vtkDataSet* pointInfo = vtkDataSet::GetData(inputVector[1]);
	vtkTable* const output = vtkTable::GetData(outputVector,0);
	output->Initialize();
	this->CalculateAndSetBounds(dataSet,pointInfo);
	// If we want to cut off at the virial radius, compute this, and remove the
	// portion of the data set we don't care about
	if(this->CutOffAtVirialRadius)
		{
		VirialRadiusInfo virialRadiusInfo = \
		 	ComputeVirialRadius(dataSet,this->Delta,this->MaxR,this->Center);
		vtkErrorMacro("virial radius is " << virialRadiusInfo.virialRadius);
		// note that if there was an error finding the virialRadius the 
		// radius returned is < 0
		//setting the dataSet to this newInput
		if(virialRadiusInfo.virialRadius>0)
			{
			dataSet = \
				GetDatasetWithinVirialRadius(virialRadiusInfo);	
			this->GenerateProfile(dataSet,output);
			dataSet->Delete();
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

//----------------------------------------------------------------------------
void vtkProfileFilter::CalculateAndSetBounds(vtkPolyData* input, 
	vtkDataSet* source)
{
	//TODO: this can later be done as in the XML documentation for this filter; 	  
	// for now, only getting the first point. this is the point selected in the
	// GUI, or the first end of the line selected in the GUI
	double* center = source->GetPoint(0);
	for(int i = 0; i < 3; ++i)
	{
		this->Center[i]=center[i];
	}
	// calculating the the max R
	this->MaxR=ComputeMaxR(input,this->Center);
	delete [] center;
}

//----------------------------------------------------------------------------
int vtkProfileFilter::GetBinNumber(double x[])
{
	double distanceToCenter = \
		sqrt(vtkMath::Distance2BetweenPoints(this->Center,x));
	return floor(distanceToCenter/this->BinSpacing);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::GenerateProfile(vtkPolyData* input,vtkTable* output)
{	
	this->InitializeBins(input,output);
	this->ComputeStatistics(input,output);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::InitializeBins(vtkPolyData* input,
	vtkTable* output)
{
	this->CalculateAndSetBinExtents(input);
	vtkSmartPointer<vtkDataArray> nextArray;
	for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
		{
		nextArray = input->GetPointData()->GetArray(i);
		string baseName = nextArray->GetName();
		for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
			{
			string totalName = GetColumnName(baseName,TOTAL,comp); 
			// Allocating an column for the total sum of the existing quantities
			AllocateDataArray(output,totalName.c_str(),1,this->BinNumber);
			// Allocating an column for the averages of the existing quantities
			string averageName = GetColumnName(baseName,AVERAGE,comp); 
			AllocateDataArray(output,averageName.c_str(),1,this->BinNumber);
			if(this->CumulativeQuantities->LookupValue(baseName)>=0)
				{
				// we should also consider this a cumulative quantity
				string cumulativeName = GetColumnName(baseName,CUMULATIVE,comp);
				AllocateDataArray(output,cumulativeName.c_str(),1,this->BinNumber);
				}
			}
		}
	// For our additional quantities, allocating a column for the average 
	// and the sum. These are restricted to be scalars
	for(int i = 0; 
		i < this->AdditionalProfileQuantities->GetNumberOfValues();
	 	++i)
		{
		string baseName=this->AdditionalProfileQuantities->GetValue(i);
		string totalName = GetColumnName(baseName,TOTAL,0); 
		// Allocating an column for the total sum of the existing quantities
		AllocateDataArray(output,totalName.c_str(),1,this->BinNumber);
		// Allocating an column for the averages of the existing quantities
		string averageName = GetColumnName(baseName,AVERAGE,0); 
		AllocateDataArray(output,averageName.c_str(),1,this->BinNumber);
		if(this->CumulativeQuantities->LookupValue(baseName)>=0)
			{
			// we should also consider this a cumulative quantity
			string cumulativeName = GetColumnName(baseName,CUMULATIVE,0);
			AllocateDataArray(output,cumulativeName.c_str(),1,this->BinNumber);
			}
		}
}

//----------------------------------------------------------------------------
void vtkProfileFilter::CalculateAndSetBinExtents(vtkPolyData* input)
{
	// TODO: only supporting non-log bin spacing for now
	this->BinSpacing=this->MaxR/this->BinNumber;
}

//----------------------------------------------------------------------------
void vtkProfileFilter::ComputeStatistics(vtkPolyData* input,vtkTable* output)
{
	for(int nextPointId = 0;
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();
	 		++nextPointId)
		{
			this->UpdateBinStatistics(input,nextPointId,output);
		}
	// Updating averages and doing relevant postprocessing
	this->BinAveragesAndPostprocessing(output);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBinStatistics(vtkPolyData* input,
 	vtkIdType pointGlobalId,vtkTable* output)
{
	double* x = GetPoint(input,pointGlobalId);
	// As we bin by radius always need
	double* r=PointVectorDifference(x,this->Center);
	// Many of the quantities explicitely require the velocity
	double* v=GetDataValue(input,"velocity",pointGlobalId);
	// TODO: debug this function, for now placing everything in
	// bin 0 as this is returning results > number of bins!
	// DO AN ASSERT TO MAKE SURE THIS DOESN'T HAPPEN!
	int binNum=this->GetBinNumber(x);
	assert(0<=binNum<=this->BinNum);
	// Updating quanties for the input data arrays
	for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
		{
		vtkSmartPointer<vtkDataArray> nextArray = \
		 	input->GetPointData()->GetArray(i);
		// getting the data for this point
		double* nextData = GetDataValue(input,
			nextArray->GetName(),pointGlobalId);
		// Updating the total bin
		for(int comp = 0; comp <	nextArray->GetNumberOfComponents(); ++comp)
			{
			string baseName = nextArray->GetName();
			string totalName =GetColumnName(baseName,TOTAL,comp);		
//			cout << "updating bin " << binNum << " total name " << totalName \
				<< "with data " << nextData[comp] << "\n";
			this->UpdateBin(binNum,ADD,totalName, nextData[comp], output);
			if(this->CumulativeQuantities->LookupValue(baseName)>=0)
				{
				// we should also consider this a cumulative quantity
				string cumulativeName =GetColumnName(baseName,CUMULATIVE,comp); 
				this->UpdateCumulativeBins(binNum,ADD,cumulativeName,
					nextData[comp],output);
				}
			}
		//Finally some memory management
		delete [] nextData;
		}
	// For our additional quantities, allocating a column for the average 
	// and the sum. These are restricted to be scalars
	for(int i = 0; i < 
		this->AdditionalProfileQuantities->GetNumberOfValues(); 
		++i)
		{
		string baseName=this->AdditionalProfileQuantities->GetValue(i);
		double* additionalData = \
			this->CalculateAdditionalProfileQuantity(baseName,v,r);
		// updating the totalbin
		string totalName=GetColumnName(baseName,TOTAL,0);
		this->UpdateBin(binNum,ADD,totalName,additionalData[0],output);
		if(this->CumulativeQuantities->LookupValue(baseName)>=0)
			{
			// we should also consider this a cumulative quantity
			string cumulativeName =GetColumnName(baseName,CUMULATIVE,0); 
			this->UpdateCumulativeBins(binNum,ADD,cumulativeName,
				additionalData[0],output);
			}
			//Finally some memory management
			delete [] additionalData;
		}
	// Finally some memory management
	delete [] x;
	delete [] r;
	delete [] v;
}




//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBin(int binNum, BinUpdateType updateType,
 	string attributeName, double dataToAdd, vtkTable* output)
{
	double data=output->GetValueByName(binNum,attributeName.c_str()).ToDouble();
	switch(updateType)
		{
		case ADD:
			data+=dataToAdd;
			break;
		case MULTIPLY:
			data*=dataToAdd;
			break;
		}
	output->SetValueByName(binNum,attributeName.c_str(),data);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateCumulativeBins(int binNum, BinUpdateType 		
	updateType, string attributeName, double dataToAdd, vtkTable* output)
{
	for(int bin = binNum; bin < this->BinNumber; ++bin)
		{
		this->UpdateBin(bin,updateType,attributeName,dataToAdd,output);
		}
}

//----------------------------------------------------------------------------
double* vtkProfileFilter::CalculateAdditionalProfileQuantity(
	string additionalQuantityName, double v[], double r[])
{
	// Some inefficiency by recomputing quantities, but paid for with 
	// flexibility, i.e don't have to call in a certain order or can compute
	// one quantity without storing the other
	// TODO: probably not string == here, just a proxy for
	// the function I should use throwing unimplemented until
	// I implement
	if(additionalQuantityName == "radial velocity")
		{
		return ComputeRadialVelocity(v,r);
		}
	else if(additionalQuantityName == "tangential velocity")
		{
		return ComputeTangentialVelocity(v,r);
		}
	else if(additionalQuantityName == "angular momentum")
		{
		return ComputeAngularMomentum(v,r);
		}
	else if(additionalQuantityName == "velocity dispersion")
		{
		return ComputeVelocitySquared(v,r);
		}
	else if(additionalQuantityName == "radial velocity dispersion")
		{
		return ComputeRadialVelocitySquared(v,r);
		}
	else if(additionalQuantityName == "tangential velocity dispersion")
		{
		return ComputeTangentialVelocitySquared(v,r);
		}
	else if(additionalQuantityName=="number in bin")
		{
			double* numberInBinDummy = new double[1];
			numberInBinDummy[0]=1.0;
			return numberInBinDummy;
		}
	else
		{
		vtkWarningMacro("input arrray requested not found, quantity returned as \
			array of zero");
		double* numberInBinDummy = new double[1];
		numberInBinDummy[0]=1.0;
		return numberInBinDummy;
		}
}

//----------------------------------------------------------------------------
void 	vtkProfileFilter::BinAveragesAndPostprocessing(vtkTable* output)
{
	// TODO: implement
//	throw "unimplemented"
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

string vtkProfileFilter::GetColumnName(string baseName, 
	ColumnType columnType, int dataIndex)
{
	switch(columnType)
		{
		case AVERAGE:
			return baseName+"_"+ToString(dataIndex)+"_average";
		case TOTAL:
			return baseName+"_"+ToString(dataIndex)+"_total";
		case CUMULATIVE:
			return baseName+"_"+ToString(dataIndex)+"_cumulative";
		default:
			vtkDebugMacro("columnType not found for function GetColumnName, returning error string");
			return "error";
		}
}














