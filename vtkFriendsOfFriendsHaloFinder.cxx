/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFriendsOfFriendsHaloFinder.cxx,v $
=========================================================================*/
#include "vtkFriendsOfFriendsHaloFinder.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkPKdTree.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkGenericPointIterator.h"
#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkPKdTree.h"
#include "vtkDistributedDataFilter.h"
#include "vtkCallbackCommand.h"
#include <vtkstd/vector>
#include "astrovizhelpers/DataSetHelpers.h"


vtkCxxRevisionMacro(vtkFriendsOfFriendsHaloFinder, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkFriendsOfFriendsHaloFinder);
vtkCxxSetObjectMacro(vtkFriendsOfFriendsHaloFinder,Controller,
	vtkMultiProcessController);
vtkCxxSetObjectMacro(vtkFriendsOfFriendsHaloFinder,PKdTree,vtkPKdTree);
vtkCxxSetObjectMacro(vtkFriendsOfFriendsHaloFinder,D3,
	vtkDistributedDataFilter);

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::vtkFriendsOfFriendsHaloFinder()
{
  this->LinkingLength = 1e-6; //default
	this->PKdTree  = NULL;
	this->Controller = NULL;
	this->D3 = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::~vtkFriendsOfFriendsHaloFinder()
{
  this->SetPKdTree(NULL);
  this->SetController(NULL);
  this->SetD3(NULL);
}

//----------------------------------------------------------------------------
void vtkFriendsOfFriendsHaloFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Linking Length: " << this->LinkingLength << "\n";
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

int vtkFriendsOfFriendsHaloFinder::FindHaloes(vtkPointSet* input,
	vtkPointSet* output)
{
	// Building the Kd tree, should already be built
	vtkSmartPointer<vtkPKdTree> pointTree = vtkSmartPointer<vtkPKdTree>::New();
		pointTree->BuildLocatorFromPoints(input);
	// calculating the initial haloes
	vtkSmartPointer<vtkIdTypeArray> haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");
	/*
	int haloId=0;
	for(int nextPointId = 0;
		nextPointId < input->GetPoints()->GetNumberOfPoints();
	 	++nextPointId)
		{
		if(haloIdArray->GetValue(nextPointId)>0)
			{
			// this means this point is already in a unique halo
			continue;
			}
		double* nextPoint=GetPoint(input,nextPointId);
		// finding the points within a linking length of this point
		vtkSmartPointer<vtkIdList> pointsInLinkingLength = \
			vtkSmartPointer<vtkIdList>::New();
		pointTree->FindPointsWithinRadius(
			this->LinkingLength,
			nextPoint,
			pointsInLinkingLength);
		// TODO: change, just the prototype
		// increment halo id if necessary
		haloId+=(pointsInLinkingLength->GetNumberOfIds()>0) ? 1 : 0;
		for(int neighborPointLocalId = 0;
	 		neighborPointLocalId < pointsInLinkingLength->GetNumberOfIds();
			++neighborPointLocalId)
			{
			vtkIdType neighborPointGlobalId  = \
				pointsInLinkingLength->GetId(neighborPointLocalId);
			haloIdArray->SetValue(neighborPointLocalId,haloId);
			}
		// find bounding box 		
		// for each point in bounding box (8), find points within radius tau
			// merge with this id list,
			// compute new bounding box
			// continue until bounding box converges			
		// when done, mark all these particles as halo and as done in the 
		// vector by assigning them unique ID
		
		// Finally, some memory management
		delete [] nextPoint;
		}
	// merging the haloes
	// recording the results
	*/
	output->GetPointData()->AddArray(haloIdArray);
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet* output = vtkPointSet::GetData(outputVector);
  output->ShallowCopy(input);
	this->FindHaloes(input,output);
	// Finally, some memory management
  output->Squeeze();
  return 1;
}
