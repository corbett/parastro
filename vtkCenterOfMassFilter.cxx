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
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);
//  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// Place result in output
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	// TODO: trying this out
//	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
//	output->SetVerts(vtkSmartPointer<vtkCellArray>::New());
	output->ShallowCopy(input);
	// Allocating data arrays and setting to zero
	//TODO: add back in
/*
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
			for(int proc = 1; proc < numProc; ++proc)
				{
				double* recTotalMass;
				double* recTotalWeightedMass;
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
			}
		}
	else
		{
		// we aren't using MPI or have only one process
		UpdateCOMVars(input,totalMass[0],totalWeightedMass);
		}
	// we will create one point in the output: the center of mass point
	double* dbCenterOfMass = ComputeCOM(input,totalMass[0],totalWeightedMass);
	float* centerOfMass = new float[3];
	for(int i = 0; i < 3; ++i)
		{
		cout << " com " << i << " is " << dbCenterOfMass[i] << "\n";
		centerOfMass[i]=static_cast<float>(dbCenterOfMass[i]);
		}
	// Placing the point's data in the output
	// TODO: add back in
//	SetPointValue(output,centerOfMass); 
	// finally, some memory management
	delete [] totalMass;
	delete [] totalWeightedMass;
	delete [] dbCenterOfMass;
	delete [] centerOfMass;
	*/
	cout << " totally done\n";
  return 1;
}
