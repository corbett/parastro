/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMostBoundFilter.cxx,v $
=========================================================================*/
#include "vtkMostBoundFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"


vtkCxxRevisionMacro(vtkMostBoundFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkMostBoundFilter);

//----------------------------------------------------------------------------
vtkMostBoundFilter::vtkMostBoundFilter()
{
	this->Overdensity = 0; // this file is also optional
}

//----------------------------------------------------------------------------
vtkMostBoundFilter::~vtkMostBoundFilter()
{
}

//----------------------------------------------------------------------------
void vtkMostBoundFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkMostBoundFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

double* vtkMostBoundFilter::GetMostBoundParticle(vtkPointSet* input)
{
	double* mostBoundParticle = new double[3];
	//TODO: implement
	mostBoundParticle[0]=0;
	mostBoundParticle[1]=0;
	mostBoundParticle[2]=0;
	return mostBoundParticle;
}

//----------------------------------------------------------------------------
int vtkMostBoundFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);	
	// Calculate the most bound particle
	double* dbMostBound = GetMostBoundParticle(input); // db to calc virial r
	float mostBound[3]; // float to store in data structure
	for(int i = 0; i < 3; ++i)
		{
		mostBound[i]=static_cast<float>(dbMostBound[i]);
		}
	// Initialize the output structure. If we are able to calculate
	// a virial radius from the user defined overdensity, we will
	// display a sphere about the mostBoundParticle at the virial radius.
	// If not, we will only display the most bound particle
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	if(this->Overdensity>0)
		{
			double maxR=ComputeMaxR(input,dbMostBound);
			VirialRadiusInfo virialRadiusInfo=\
								ComputeVirialRadius(input,this->Overdensity,maxR,dbMostBound);
			if(virialRadiusInfo.virialRadius>0)
				{
				//Here is where we create the sphere around the COM to display
				vtkWarningMacro("the virial radius is " 
												<< virialRadiusInfo.virialRadius);
				// Creating the sphere
				CreateSphere(output,\
										virialRadiusInfo.virialRadius,dbMostBound);
				}
			else
				{
				vtkWarningMacro("unable to find the virial radius from over density you specified. Perhaps it is too high. For now displaying only the center of mass");
				// Placing the point's data in the output
				SetPointValue(output,mostBound);					
				}
		}
	else
		{
			// Placing the point's data in the output
			SetPointValue(output,mostBound);
		}
	// finally, some memory management
	delete [] dbMostBound;
  return 1;
}
