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
	double* centerOfMass = \
	 	centerOfMassFilter->ComputeCenterOfMass(input,"mass");
	double inertiaTensor[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];
	if(centerOfMass!=NULL)
		{
		// TODO: implement
		// we are at process 0 or running in serial
		// first broadcast to other processes if 
		// necessary
		if(RunInParallel(this->Controller))
			{
			
			}
		// computing the moment of inertia tensor 3x3 matrix, and its
		// eigenvalues and eigenvectors
		ComputeInertiaTensor(input,centerOfMass,inertiaTensor);
		// receive values from other processes if necessary
		if(RunInParallel(this->Controller))
			{
			
			}
		// finally perform final computation
		vtkMath::Diagonalize3x3(inertiaTensor,eigenvalues,eigenvectors);
		// displaying eigenvectors
		DisplayVectorsAsLines(input,output,eigenvectors,centerOfMass);
		// memory management
		delete [] centerOfMass;
		}
	else
		{
		// TODO: implement
		// we are running in parallel

		}
  return 1;
}
