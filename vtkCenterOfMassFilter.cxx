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


vtkCxxRevisionMacro(vtkCenterOfMassFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkCenterOfMassFilter);

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::vtkCenterOfMassFilter()
{
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
	double totalMass,totalWeightedMass[3];
  vtkDebugMacro("2. Calculating the quantities \
									we are interested in.");
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		vtkDebugMacro("next point is " << nextPoint[0] << "," \
									<< nextPoint[1] << ","<< nextPoint[2]);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		// GetTuple function which operates with float
		double* mass=GetDataValue(input,"mass",nextPointId);
		//calculating the weighted mass
		double* weightedMass=CalculateWeightedMass(mass[0],nextPoint);
		// updating the mass and the weighted mass
		for(int i = 0; i < 3; ++i)
		{
			totalMass+=mass[0];
			totalWeightedMass[i]+=weightedMass[i];
		}
		// Finally, some memory management
		delete [] weightedMass;
		delete [] mass;
		delete [] nextPoint;
		}
	// calculating the result
	// our final data is in float, as Tipsy's data is stored in float
	float centerOfMass[3];
	if(totalMass!=0)
		{
		for(int i = 0; i < 3; ++i)
			{
			centerOfMass[i]=static_cast<float>(totalWeightedMass[i]/totalMass);
			vtkDebugMacro("center of mass is " << centerOfMass[i]);
			}
		}
	else
		{
		vtkErrorMacro("total mass is zero, cannot calculate center of mass");
		return 0;
		}
	// Displaying/storing the result
	// we will create one point in the output: the center of mass point
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	// Placing the point's data in the output
	SetPointValue(output,centerOfMass);
  return 1;
}
