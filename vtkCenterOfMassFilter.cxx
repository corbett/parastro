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
int vtkCenterOfMassFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	// Allocating data arrays and setting to zero
	double* totalMass =new double[0];
	totalMass[0]=0;
	double* totalWeightedMass = new double[3];
	for(int i = 0; i < 3; ++i)
		{
		totalWeightedMass[i]=0;
		}
	if (this->Controller != NULL && 
		this->Controller->GetNumberOfProcesses() > 1)
		{
		int procId=this->Controller->GetLocalProcessId();
		int numProc=this->Controller->GetNumberOfProcesses();
		if(procId!=0)
			{
			// We are at non-root process so simply update and move on
			// Private variables to aid computation of COM
			UpdateCOMVars(input,totalMass[0],totalWeightedMass);
			// Sending to root
			this->Controller->Send(totalMass,1,0,TOTAL_MASS);
			this->Controller->Send(totalWeightedMass,3,0,TOTAL_WEIGHTED_MASS);
			return 1;
			}
		else
			{
			// We are at root process so update results from root process 
			UpdateCOMVars(input,totalMass[0],totalWeightedMass);
			// Now gather results from each process other than this one
			// TODO: add back in
			for(int proc = 1; proc < numProc; ++proc)
				{
				double* recTotalMass;
				double* recTotalWeightedMass;
				// Receiving
				cout << " receiving data from proc " << proc << "\n";
				this->Controller->Receive(recTotalMass,
					1,proc,TOTAL_MASS);
				this->Controller->Receive(recTotalWeightedMass,
					3,proc,TOTAL_WEIGHTED_MASS);
				// TODO: add back in
				/*
				// Updating
				totalMass[0]+=recTotalMass[0];
				for(int i = 0; i < 3; ++i)
					{
					totalWeightedMass[i]+=recTotalWeightedMass[i];
					}
				delete [] recTotalMass;
				delete [] recTotalWeightedMass;
				*/
				}
			}
		}
	else
		{
		// we aren't using MPI or have only one process
		UpdateCOMVars(input,totalMass[0],totalWeightedMass);
		}

	// Place result in output
	double* dbCenterOfMass = ComputeCOM(input,
		totalMass[0],totalWeightedMass);
	float* centerOfMass = new float[3];
	for(int i = 0; i < 3; ++i)
		{
		// TODO: add back in
		//centerOfMass[i]=static_cast<float>(dbCenterOfMass[i]);
		centerOfMass[i]=0;
		}
	// we will create one point in the output: the center of mass point
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	// Placing the point's data in the output
	SetPointValue(output,centerOfMass);
	// finally, some memory management
	delete [] totalMass;
	delete [] totalWeightedMass;
	delete [] dbCenterOfMass;
	delete [] centerOfMass;
  return 1;
}
