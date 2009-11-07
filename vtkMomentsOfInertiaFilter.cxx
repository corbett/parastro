/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMomentsOfInertiaFilter.cxx,v $
=========================================================================*/
#include "vtkMomentsOfInertiaFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkMultiProcessController.h"
#include "vtkCenterOfMassFilter.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"
#include "vtkMath.h"


vtkCxxRevisionMacro(vtkMomentsOfInertiaFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkMomentsOfInertiaFilter);
vtkCxxSetObjectMacro(vtkMomentsOfInertiaFilter,Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::vtkMomentsOfInertiaFilter()
{
	this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::~vtkMomentsOfInertiaFilter()
{
 	this->SetController(0);
}

//----------------------------------------------------------------------------
void vtkMomentsOfInertiaFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkMomentsOfInertiaFilter::FillInputPortInformation(int, 
	vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}


//----------------------------------------------------------------------------
int vtkMomentsOfInertiaFilter::RequestData(vtkInformation*,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get input and output data
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	// computing the center of mass, works in parallel if necessary
	vtkSmartPointer<vtkCenterOfMassFilter> centerOfMassFilter = \
		vtkSmartPointer<vtkCenterOfMassFilter>::New();
	centerOfMassFilter->SetController(this->Controller);
	// will be != null only for root process or serial
	double* centerOfMass = new double[3];
	centerOfMass=centerOfMassFilter->ComputeCenterOfMass(input,"mass"); 
	if(RunInParallel(this->Controller))
		{
		// syncs the value of centerofmass among all processes
		this->Controller->Broadcast(centerOfMass,3,0);
		cout << "COM on proc " << this->Controller->GetLocalProcessId() 
			<< " is synced to be " << centerOfMass[0] << ","
			<< centerOfMass[1] << "," << centerOfMass[2] << "\n";
		}
	double inertiaTensor[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];
	if(centerOfMass!=NULL)
		{
		// TODO: finish implementation
		// computing the moment of inertia tensor 3x3 matrix, and its
		// eigenvalues and eigenvectors
		ComputeInertiaTensor(input,centerOfMass,inertiaTensor);
		// if we are not at proc 0 send results
		// if we are on proc 0 & running in parallel receive results 
		// from other processes
		// finally perform final computation
		vtkMath::Diagonalize3x3(inertiaTensor,eigenvalues,eigenvectors);
		// displaying eigenvectors
		DisplayVectorsAsLines(input,output,eigenvectors,centerOfMass);
		// memory management
		delete [] centerOfMass;
		}
  return 1;
}
