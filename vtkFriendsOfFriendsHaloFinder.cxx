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
	// Allocating array to store the number of the halo each point belongs to
 	AllocateIntDataArray(output,"halo number", 
		1,output->GetPoints()->GetNumberOfPoints());
	// calculating the initial haloes
	vtkstd::vector<vtkIdList> initialHaloes;
	vtkstd::vector<int> done;
	haloId=1;
	for(int nextPointId = 0;
		nextPointId < input->GetPoints()->GetNumberOfPoints();
	 	++nextPointId)
		{
		if(done[nextPoint]>0)
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
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet* output = vtkPointSet::GetData(outputVector);
	this->FindHaloes(input,output);
  output->ShallowCopy(input);
	// Finally, some memory management
  output->Squeeze();
  return 1;
}
