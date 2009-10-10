/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCenterOfMassFilter.cxx,v $
=========================================================================*/
#include "vtkCenterOfMassFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"


vtkCxxRevisionMacro(vtkCenterOfMassFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkCenterOfMassFilter);

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::vtkCenterOfMassFilter()
{
	this->Overdensity = 0; // this file is also optional
}

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::~vtkCenterOfMassFilter()
{
}

//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

double* vtkCenterOfMassFilter::CalculateWeightedMass(double& mass,\
																double* point)
{
	double* weightedMass = new double[3];
	for(int i = 0; i < 3; ++i)
	{
		weightedMass[i]=mass*point[i];
	}
	return weightedMass;
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	double totalMass=0;
	double totalWeightedMass[3];
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		// GetTuple function which operates with float
		double* mass=GetDataValue(input,"mass",nextPointId);
		//calculating the weighted mass
		double* weightedMass=CalculateWeightedMass(mass[0],nextPoint);
		// updating the mass and the weighted mass
		totalMass+=mass[0];
		for(int i = 0; i < 3; ++i)
			{
			totalWeightedMass[i]+=weightedMass[i];
			}
		// Finally, some memory management
		delete [] weightedMass;
		delete [] mass;
		delete [] nextPoint;
		}
	// calculating the result
	// our final data is in float, as Tipsy's data is stored in float
	float* centerOfMass=new float[3];
	if(totalMass!=0)
		{
		for(int i = 0; i < 3; ++i)
			{
			// casting to float to be the same precision as other Tipsy
			// variables.	
			centerOfMass[i]=\
				static_cast<float>(totalWeightedMass[i]/totalMass);
			}
		}
	else
		{
		vtkErrorMacro("total mass is zero, cannot calculate center of mass");
		return 0;
		}			
	// we will create one point in the output: the center of mass point
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	// if the Overdensity is non zero and we are able to find a
	// virial radius then we set the output to the sphere
	// around the COM at the virial radius.
	if(this->Overdensity>0)
		{
			VirialRadiusInfo virialRadiusInfo=\
										ComputeVirialRadius(input,this->Overdensity,centerOfMass);
			if(virialRadiusInfo.virialRadius>0)
				{
				//Here is where we create the sphere around the COM to display
				// TODO: do that
				vtkWarningMacro("the virial radius is " 
												<< virialRadiusInfo.virialRadius);
				}
			else
				{
				vtkWarningMacro("unable to find the virial radius from over density you specified. Perhaps it is too high. For now displaying only the center of mass");
				// Placing the point's data in the output
				SetPointValue(output,centerOfMass);					
				}
		}
	else
		{
			// Placing the point's data in the output
			SetPointValue(output,centerOfMass);
		}
	// finally, some memory management
	delete [] centerOfMass;
  return 1;
}
