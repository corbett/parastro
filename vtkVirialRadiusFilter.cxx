/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVirialRadiusFilter.cxx,v $
=========================================================================*/
#include "vtkVirialRadiusFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"
#include "vtkCellData.h"
#include "vtkSortDataArray.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkInformationDataObjectKey.h"
#include "vtkPolyData.h" // helper functions take this as argument
#include <cmath>
using vtkstd::string;

vtkCxxRevisionMacro(vtkVirialRadiusFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkVirialRadiusFilter);

//----------------------------------------------------------------------------
vtkVirialRadiusFilter::vtkVirialRadiusFilter()
{
  this->SetNumberOfInputPorts(2);
	// Defaults for quantities which will be computed based on user's
	// later input
	this->MaxR=1.0;
	this->Delta=0.0;
}

//----------------------------------------------------------------------------
vtkVirialRadiusFilter::~vtkVirialRadiusFilter()
{
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "overdensity: " << this->Delta << "\n"
		<< "softening :" << this->Softening << "\n";
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkVirialRadiusFilter::FillInputPortInformation (int port, 
	vtkInformation *info)
{
  this->Superclass::FillInputPortInformation(port, info);
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkVirialRadiusFilter::RequestData(vtkInformation *request,
																	vtkInformationVector **inputVector,
																	vtkInformationVector *outputVector)
{
	// Now we can get the input with which we want to work
 	vtkPolyData* dataSet = vtkPolyData::GetData(inputVector[0]);
	// Setting the center based upon the selection in the GUI
	vtkDataSet* pointInfo = vtkDataSet::GetData(inputVector[1]);
	vtkPointSet* output = vtkPointSet::GetData(outputVector,0);
	output->ShallowCopy(dataSet);
	
/*	
	this->CalculateAndSetBounds(dataSet,pointInfo);
	VirialRadiusInfo virialRadiusInfo = \
	 	ComputeVirialRadius(dataSet,this->Softening,
		this->Delta,this->MaxR,this->Center);
		// note that if there was an error finding the virialRadius the 
		// radius returned is < 0
	if(virialRadiusInfo.virialRadius>0)
		{
		vtkWarningMacro("virial radius is " << virialRadiusInfo.virialRadius);
		//setting the dataSet to this newInput
		vtkPolyData* newDataSet = \
			GetDatasetWithinVirialRadius(virialRadiusInfo);
		// TODO: add back in
		//output->DeepCopy(newDataSet);
		//newDataSet->Delete();
		}
	else	
		{
		vtkErrorMacro("Unable to find virial radius: considering changing your delta or selecting a different point around which to search. For now simply copying input");
		}
	*/
	return 1;
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::CalculateAndSetBounds(vtkPolyData* input, 
	vtkDataSet* source)
{
	//TODO: this can later be done as in the XML documentation for this filter; 	  
	// for now, only getting the first point. this is the point selected in the
 // GUI, or the midpoint of the line selected in the GUI
	double* center;
	if(source->GetNumberOfPoints()==1)
		{
		// we are dealing with a point
		center = source->GetPoint(0);
		}
	else
		{
		// we are dealing with a line
		double* pointOne=source->GetPoint(0);
		double* pointTwo=source->GetPoint(source->GetNumberOfPoints()-1);
		// TODO: fix this is currently == pointTwo (for some reason p1=p2?)
		center=ComputeMidpoint(pointOne,pointTwo);
		}
	for(int i = 0; i < 3; ++i)
		{
		this->Center[i]=center[i];
		}
 	// calculating the the max R
	this->MaxR=ComputeMaxR(input,this->Center);
}

