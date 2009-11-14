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
#include "vtkCallbackCommand.h"
#include <vtkstd/vector>
#include <vtkstd/map>
#include "astrovizhelpers/DataSetHelpers.h"


vtkCxxRevisionMacro(vtkFriendsOfFriendsHaloFinder, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkFriendsOfFriendsHaloFinder);

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::vtkFriendsOfFriendsHaloFinder():vtkDistributedDataFilter()
{
  this->LinkingLength = 1e-6; //default
	this->MinimumNumberOfParticles = 50; // default
}

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::~vtkFriendsOfFriendsHaloFinder()
{
}

//----------------------------------------------------------------------------
void vtkFriendsOfFriendsHaloFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Linking Length: " << this->LinkingLength 
		<<	indent << "Minimum Number Of Particles: " 
		<<  this->MinimumNumberOfParticles << "\n";
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FillInputPortInformation(int, 
  vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FindHaloes(vtkPKdTree* pointTree,
	vtkPointSet* output)
{
	if(this->MinimumNumberOfParticles < 2)
		{
		vtkWarningMacro("setting minimum number of particles to 2, a minimum number of particles below this makes no sense.");
		this->MinimumNumberOfParticles=2;
		}
	// building a locator
	pointTree->BuildLocatorFromPoints(output);
	// calculating the initial haloes- yes it really is done in just this 
	// one line.
	vtkSmartPointer<vtkIdTypeArray> haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(output->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");

	vtkstd::vector<vtkIdTypeArray*> allHaloIdArrays;
	if(RunInParallel(this->GetController()))
		{
		int procId=this->Controller->GetLocalProcessId();
		int numProc=this->Controller->GetNumberOfProcesses();
		if(procId!=0)
			{
			this->GetController()->Send(haloIdArray,0,HALO_ID_ARRAY_INITIAL);
			// waiting to recieve the final result, as computed by root
			this->GetController()->Receive(haloIdArray,0,
				HALO_ID_ARRAY_FINAL);
			// setting output
			output->GetPointData()->AddArray(haloIdArray);
			// returning
			return 1;
			}
		else
			{
			// Syncing all hashtables if process 0
			// Filling it first with process zero's info
			allHaloIdArrays.push_back(haloIdArray);
			vtkSmartPointer<vtkIdTypeArray> recHaloIdArray = \
				vtkSmartPointer<vtkIdTypeArray>::New();
			for(int proc = 1; proc < numProc; ++proc)
				{
				this->GetController()->Receive(recHaloIdArray,proc,
					HALO_ID_ARRAY_INITIAL);
				allHaloIdArrays.push_back(recHaloIdArray);
				}
			// don't return, proc 0 should execute code after if statement
			}
		}
	else
		{
		// running in serial
		allHaloIdArrays.push_back(haloIdArray);
		}
	// Now assign halos, if this point has at least one other pair,
	// it is a halo, if not it is not (set to 0)
	// first building map of id to count of that id, O(N)
	vtkstd::map<vtkIdType,int> haloCount;
	int uniqueId=1;
	for(int procHaloIdArrayIndex = 0; 
		procHaloIdArrayIndex < allHaloIdArrays.size(); 
		++procHaloIdArrayIndex)
		{
		vtkSmartPointer<vtkIdTypeArray> nextHaloIdArray = \
			allHaloIdArrays[procHaloIdArrayIndex];
		// use negatives to figure differentiate between unique id assignment 
		// (positive) and count (negative) in the same map
		for(int nextHaloId = 0;
			nextHaloId < nextHaloIdArray->GetNumberOfTuples();
		 	++nextHaloId)
			{
			vtkIdType haloId = nextHaloIdArray->GetValue(nextHaloId);
			if(haloCount[haloId]==-1*this->MinimumNumberOfParticles)
				{
				// we have seen the id minimum number of particles times,
				// so we assign it a unique halo id, considering it a halo
				haloCount[haloId]=uniqueId;
				uniqueId+=1;
				}
			else if(haloCount[haloId]<1)
				{
				// this counts each time we see the id, until we reach the
				// minimum number of particles to count as a halo
				haloCount[haloId]-=1;
				}
			}
		}
	// finally setting to zero points which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id
	// process 0 or running in serial
	for(int procHaloIdArrayIndex = 0; 
		procHaloIdArrayIndex < allHaloIdArrays.size(); 
		++procHaloIdArrayIndex)
		{
		vtkSmartPointer<vtkIdTypeArray> nextHaloIdArray = \
			allHaloIdArrays[procHaloIdArrayIndex];
		// use negatives to figure differentiate between unique id assignment 
		// (positive) and count (negative) in the same map
		for(int nextHaloId = 0;
			nextHaloId < nextHaloIdArray->GetNumberOfTuples();
		 	++nextHaloId)
			{
			vtkIdType haloId = nextHaloIdArray->GetValue(nextHaloId);
			if(haloCount[haloId]<1)
				{
				// we only saw it less than requisite number of times
				nextHaloIdArray->SetValue(nextHaloId,0);
				}
			else
				{
				// we saw it more than once, assign it to its unique id
				nextHaloIdArray->SetValue(nextHaloId,haloCount[haloId]);
				}
			}
		if(RunInParallel(this->GetController()) && procHaloIdArrayIndex > 0)
			{
			// if running in parallel and if we are not dealing with our own 
			// halo array on process 0
			// dispatch it to its process processes
			this->GetController()->Send(haloIdArray,
				procHaloIdArrayIndex,HALO_ID_ARRAY_FINAL);
			}
		}
	output->GetPointData()->AddArray(haloIdArray);
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	vtkPointSet* output;
	vtkSmartPointer<vtkPKdTree> pointTree;
	if(RunInParallel(this->GetController()))
		{
		// call D3, setting retain PKTree to 1; this can be accessed by later
		// methods
		this->RetainKdtreeOn();
		// Just calling the superclass' method to distribute data and build
		// PKdTree
	  this->Superclass::RequestData(request,inputVector,outputVector);
		output = vtkPointSet::GetData(outputVector);
		// setting the KdTree to the output from D3
		pointTree=this->GetKdtree();
		}
	else
		{		
		output = vtkPointSet::GetData(outputVector);
  	output->ShallowCopy(input);
		// Building the Kd tree, should already be built
		pointTree = vtkSmartPointer<vtkPKdTree>::New();
		}
	this->FindHaloes(pointTree,output);
  return 1;
}
