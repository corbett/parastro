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
	this->AdditionalProfileQuantities = vtkStringArray::New();
		this->AdditionalProfileQuantities->InsertNextValue("circular velocity");
		this->AdditionalProfileQuantities->InsertNextValue("density");
		this->AdditionalProfileQuantities->InsertNextValue("radial velocity");
		this->AdditionalProfileQuantities->InsertNextValue("tangential velocity");
		this->AdditionalProfileQuantities->InsertNextValue("angular momentum");
		this->AdditionalProfileQuantities->InsertNextValue("velocity squared");
		this->AdditionalProfileQuantities->InsertNextValue("velocity dispersion");
		this->AdditionalProfileQuantities->InsertNextValue("radial velocity squared");
		this->AdditionalProfileQuantities->InsertNextValue("radial velocity dispersion");
	 	this->AdditionalProfileQuantities->InsertNextValue("tangential velocity squared");
	 	this->AdditionalProfileQuantities->InsertNextValue("tangential velocity dispersion");
	
	this->CumulativeQuantities = vtkStringArray::New();
		this->CumulativeQuantities->InsertNextValue("number in bin");
		this->CumulativeQuantities->InsertNextValue("mass");

	this->MaxR=1.0;
	this->Delta=0.0;
	this->BinNumber=30;
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
/*	this->CumulativeQuantities->Delete();
	this->AdditionalProfileQuantities->Delete(); */ 
	// removed this--get	paraview(54834,0xa01ef500) malloc: *** error for object 0x21c8f8c0: incorrect checksum for freed object - object was probably modified after being freed. *** set a breakpoint in malloc_error_break to debug
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
	this->CalculateAndSetBinExtents(input,output);
	vtkSmartPointer<vtkDataArray> nextArray;
	// always need this for averages
	AllocateDataArray(output,GetColumnName("number in bin",TOTAL,0).c_str(),
		1,this->BinNumber);
	AllocateDataArray(output,GetColumnName("number in bin",
		CUMULATIVE,0).c_str(),1,this->BinNumber);
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
	// and the sum. These are restricted to be 3-vectors
	for(int i = 0; 
		i < this->AdditionalProfileQuantities->GetNumberOfValues();
	 	++i)
		{
		string baseName=this->AdditionalProfileQuantities->GetValue(i);
		for(int comp = 0; comp < 3; ++comp)
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
}

//----------------------------------------------------------------------------
void vtkProfileFilter::CalculateAndSetBinExtents(vtkPolyData* input,
	vtkTable* output)
{
	this->BinSpacing=this->MaxR/this->BinNumber;
	// the first column will be the bin radius
	string binRadiusColumnName=this->GetColumnName("bin radius",
		TOTAL,0);
	AllocateDataArray(output,binRadiusColumnName.c_str(),1,this->BinNumber);
	// setting the bin radii in the output
	for(int binNum = 0; binNum < this->BinNumber; ++binNum)
	{
	this->UpdateBin(binNum,SET,binRadiusColumnName,(binNum+1)*this->BinSpacing,
		output);
	}
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
	this->BinAveragesAndPostprocessing(input,output);
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
	int binNum=this->GetBinNumber(x);
	this->UpdateBin(binNum,ADD,GetColumnName(
		"number in bin",TOTAL,0).c_str(),1,output);	
	this->UpdateCumulativeBins(binNum,ADD,GetColumnName(
		"number in bin",CUMULATIVE,0).c_str(),1,output);
	
	assert(0<=binNum<=this->BinNumber);
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
	for(int i = 0; i < 
		this->AdditionalProfileQuantities->GetNumberOfValues(); 
		++i)
		{
		string baseName=this->AdditionalProfileQuantities->GetValue(i);
		double* additionalData = \
			this->CalculateAdditionalProfileQuantity(baseName,v,r);
		// updating the totalbin
		for(int comp = 0; comp < 3; ++comp)
			{
			string totalName=GetColumnName(baseName,TOTAL,comp);
			this->UpdateBin(binNum,ADD,totalName,additionalData[comp],output);
			if(this->CumulativeQuantities->LookupValue(baseName)>=comp)
				{
				// we should also consider this a cumulative quantity
				string cumulativeName =GetColumnName(baseName,CUMULATIVE,comp); 
				this->UpdateCumulativeBins(binNum,ADD,cumulativeName,
					additionalData[comp],output);
				}
			}
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
		case SET:
			data=dataToAdd;
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
	else if(additionalQuantityName == "velocity squared")
		{
		return ComputeVelocitySquared(v,r);
		}
	else if(additionalQuantityName == "radial velocity squared")
		{
		return ComputeRadialVelocitySquared(v,r);
		}
	else if(additionalQuantityName == "tangential velocity squared")
		{
		return ComputeTangentialVelocitySquared(v,r);
		}
	else
		{
		vtkDebugMacro("input arrray requested not found, quantity returned as \
			array of zero");
		// sometimes the quantities should be 0, only updated at post-processing
		// final step
		double* emptyBin = new double[3];
		emptyBin[0]=0.0;
		emptyBin[1]=0.0;
		emptyBin[2]=0.0;
		return emptyBin;
		}
}

//----------------------------------------------------------------------------
void 	vtkProfileFilter::BinAveragesAndPostprocessing(
	vtkPolyData* input,vtkTable* output)
{
	// TODO: finish implementation
	for(int binNum = 0; binNum < this->BinNumber; ++binNum)
		{
		int binSize=output->GetValueByName(binNum,
			GetColumnName("number in bin",TOTAL,0).c_str()).ToInt();
		// For each input array, update its average column by getting
		// the total from the total column then dividing by 
		// the number in the bin. Only do this if the number in the bin
		// is greater than zero
			if(binSize>0)
				{
				vtkSmartPointer<vtkDataArray> nextArray;
				for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
					{
					nextArray = input->GetPointData()->GetArray(i);
					string baseName = nextArray->GetName();
					for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
						{
						string totalName = GetColumnName(baseName,TOTAL,comp); 
						double totalData=output->GetValueByName(binNum,
							totalName.c_str()).ToDouble();
						string averageName = GetColumnName(baseName,AVERAGE,comp); 
						this->UpdateBin(binNum,SET,averageName.c_str(),
							totalData/binSize,output);
						}	
					}
					// For each additional quantity we also divide by N in the average, 
					for(int i = 0; 
						i < this->AdditionalProfileQuantities->GetNumberOfValues();
					 	++i)
						{
						string baseName=this->AdditionalProfileQuantities->GetValue(i);
						for(int comp = 0; comp < 3; ++comp)
							{
							string totalName = GetColumnName(baseName,TOTAL,comp); 
							double totalData=output->GetValueByName(binNum,
								totalName.c_str()).ToDouble();
							string averageName = GetColumnName(baseName,AVERAGE,comp); 
							this->UpdateBin(binNum,SET,averageName.c_str(),
								totalData/binSize,output);
							}
						}
				}

		// TODO: code repetition, clean this up (perhaps keep track of
		// which quantities to compute dispersions for)
		// Special handling should come last
		/* Dispersions */
/*

		vDispColumnName=GetColumnName("velocity dispersion",AVERAGE,0);
		vRadDispColumnName=GetColumnName("radial velocity dispersion",AVERAGE,0);
		vTanDispColumnName=GetColumnName("tangential velocity dispersion",AVERAGE,0);		
		// column data
		double vAve=output->GetValueByName(binNum,
			GetColumnName("velocity",AVERAGE,0).c_str()).ToDouble();
		double vSquaredAve=output->GetValueByName(binNum,
			GetColumnName("velocity squared",AVERAGE,0).c_str()).ToDouble();

		double vRadAve=output->GetValueByName(binNum,
			GetColumnName("radial velocity",AVERAGE,0).c_str()).ToDouble();
		double vRadSquaredAve=output->GetValueByName(binNum,
			GetColumnName("radial velocity squared",AVERAGE,0).c_str()).ToDouble();

		double vTanAve=output->GetValueByName(binNum,
			GetColumnName("tangential velocity",AVERAGE,0).c_str()).ToDouble();
		double vTanSquaredAve=output->GetValueByName(binNum,
			GetColumnName("tangential velocity squared",
			AVERAGE,0).c_str()).ToDouble();

		// computing quantities
		double vDisp=vtkMath::Norm(ComputeVelocityDispersion(vSquaredAve,vAve));
		double vRadDisp = \
			vtkMath::Norm(ComputeVelocityDispersion(vRadSquaredAve,vRadAve));
		double vTanDisp = vtkMath::Norm(ComputeVelocityDispersion(vTanSquaredAve,\
			vTanAve,vTanDisp));
		// updating output
		this->UpdateBin(binNum,SET,vDispColumnName,vDisp,output);
		this->UpdateBin(binNum,SET,vRadDispColumnName,vRadDisp,output);
		this->UpdateBin(binNum,SET,vTanDispColumnName,vTanDisp,output);
*/
		/* Circular velocity and density */
		// computing the circular velocity (sqrt(M(<r)/r))
		// Column names
		// Column data
		double cumulativeMass=output->GetValueByName(binNum,
			GetColumnName("mass",CUMULATIVE,0).c_str()).ToDouble();
		double binRadius=output->GetValueByName(binNum,
			GetColumnName("bin radius",TOTAL,0).c_str()).ToDouble();
		// Computation and updating
		double circularVelocity = sqrt(cumulativeMass/binRadius);
		this->UpdateBin(binNum,SET,
			GetColumnName("circular velocity",AVERAGE,0),circularVelocity,
			output);
		// computing the density M(<r)/(4/3 pi r^3)
		double density=cumulativeMass/(4./3*vtkMath::Pi()*pow(binRadius,3));
		this->UpdateBin(binNum,SET,GetColumnName("density",AVERAGE,0),
			density,output);
		}
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














