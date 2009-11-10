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
#include "vtkSphereSource.h"
#include "vtkMultiProcessController.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"


vtkCxxRevisionMacro(vtkCenterOfMassFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkCenterOfMassFilter);

vtkCxxSetObjectMacro(vtkCenterOfMassFilter,Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::vtkCenterOfMassFilter()
{
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::~vtkCenterOfMassFilter()
{
 	this->SetController(0);
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
//----------------------------------------------------------------------------
double* vtkCenterOfMassFilter::ComputeCenterOfMassFinal(
	vtkPointSet* input,double& totalMass,double totalWeightedMass[])
{
	// calculating the result
	// our final data is in float, as Tipsy's data is stored in float
	double* dbCenterOfMass=new double[3]; // this is needed for the virial calc
	if(totalMass!=0)
		{
		for(int i = 0; i < 3; ++i)
			{
			dbCenterOfMass[i]=totalWeightedMass[i]/totalMass;
			}
		}
	else
		{
		 vtkErrorMacro("total mass is zero, cannot calculate center of mass, setting center to 0,0,0");
		for(int i = 0; i < 3; ++i)
			{
			dbCenterOfMass[i]=0;	
			}
		}
	return dbCenterOfMass;
}

//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::UpdateCenterOfMassVariables(
	vtkPointSet* input,double& totalMass,double totalWeightedMass[])
{
	for(int nextPointId = 0;
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		// extracting the mass
		// has to be double as this version of VTK doesn't have 
		// GetTuple function which operates with float
		double* mass=GetDataValue(input,"mass",nextPointId);
		//calculating the weighted mass
		double* weightedMass=this->ComputeWeightedMass(mass[0],nextPoint);
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
}
//----------------------------------------------------------------------------	
double* vtkCenterOfMassFilter::ComputeWeightedMass(double& mass,double* point)
{
	double* weightedMass = new double[3];
	for(int i = 0; i < 3; ++i)
	{
	weightedMass[i]=mass*point[i];
	}
	return weightedMass;
}

//----------------------------------------------------------------------------
double* vtkCenterOfMassFilter::ComputeCenterOfMass(vtkPointSet* input,
	vtkstd::string massArrayName)
{

	// TODO: use the massArrayName
	// Allocating data arrays and setting to zero
	double totalMass[1];
	double totalWeightedMass[3];
	// testing to make sure we can get to work with D3
	if(RunInParallel(this->Controller))
		{
		int procId=this->Controller->GetLocalProcessId();
		int numProc=this->Controller->GetNumberOfProcesses();
		if(procId!=0)
			{
			// We are at non-root process so simply update and move on
			// Private variables to aid computation of COM
			this->UpdateCenterOfMassVariables(input,totalMass[0],totalWeightedMass);
			// Sending to root
			this->Controller->Send(totalMass,1,0,TOTAL_MASS);
			this->Controller->Send(totalWeightedMass,3,0,TOTAL_WEIGHTED_MASS);
			return NULL;
			}
		else
			{
			// We are at root process so update results from root process 
			this->UpdateCenterOfMassVariables(input,totalMass[0],totalWeightedMass);
			// Now gather results from each process other than this one
			for(int proc = 1; proc < numProc; ++proc)
				{
				double recTotalMass[1];
				double recTotalWeightedMass[3];
				// Receiving
				this->Controller->Receive(recTotalMass,
					1,proc,TOTAL_MASS);
				this->Controller->Receive(recTotalWeightedMass,
					3,proc,TOTAL_WEIGHTED_MASS);
				// Updating
				totalMass[0]+=recTotalMass[0];
				for(int i = 0; i < 3; ++i)
					{
					totalWeightedMass[i]+=recTotalWeightedMass[i];
					
					}
				}
			double* centerOfMassFinal = \
			this->ComputeCenterOfMassFinal(input,totalMass[0],
				totalWeightedMass);
			return centerOfMassFinal;
			}
		}
	else
		{
		// we aren't using MPI or have only one process
		this->UpdateCenterOfMassVariables(input,totalMass[0],totalWeightedMass);
		double* centerOfMassFinal = \
			this->ComputeCenterOfMassFinal(input,totalMass[0],totalWeightedMass);
		return centerOfMassFinal;
		}
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// Place result in output
	// need to output point set, otherwise D3 is VERY unhappy
	vtkPointSet* output = vtkPointSet::GetData(outputVector);
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	// TODO: actually use array name
	double* dbCenterOfMass=this->ComputeCenterOfMass(input,"mass");
	if(dbCenterOfMass!=NULL)
		{
		cout << "not null COM is " << dbCenterOfMass[0] << ","
		 	<< dbCenterOfMass[1] << ","
			<< dbCenterOfMass[2] << "\n";
		// we are in serial or at process 0
		float* centerOfMass = DoublePointToFloat(dbCenterOfMass);
		// Placing the point's data in the output
		SetPointValue(output,centerOfMass); 
		// finally, some memory management
		delete [] dbCenterOfMass;
		delete [] centerOfMass;
		}
  return 1;
}
