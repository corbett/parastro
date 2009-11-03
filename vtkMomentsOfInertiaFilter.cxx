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
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"
#include "vtkMath.h"


vtkCxxRevisionMacro(vtkMomentsOfInertiaFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkMomentsOfInertiaFilter);

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::vtkMomentsOfInertiaFilter()
{
}

//----------------------------------------------------------------------------
vtkMomentsOfInertiaFilter::~vtkMomentsOfInertiaFilter()
{
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
	// computing the center of mass
	double* centerOfMass = ComputeCOM(input);
	// computing the moment of inertia tensor 3x3 matrix, and its
	// eigenvalues and eigenvectors
	double inertiaTensor[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];
	ComputeInertiaTensor(input,centerOfMass,inertiaTensor);
	vtkMath::Diagonalize3x3(inertiaTensor,eigenvalues,eigenvectors);
	// displaying eigenvectors
	//TODO: displaying eigenvectors
	cout << "displaying eigenvectors\n";
	DisplayVectorsAsLines(input,output,eigenvectors,centerOfMass);
	delete [] centerOfMass;
  return 1;
}
