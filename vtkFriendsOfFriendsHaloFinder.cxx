/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFriendsOfFriendsHaloFinder.cxx,v $
=========================================================================*/
#include "vtkFriendsOfFriendsHaloFinder.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkKdTree.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkGenericPointIterator.h"
#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkCallbackCommand.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <vtkstd/vector>
#include <vtkstd/map>


vtkCxxRevisionMacro(vtkFriendsOfFriendsHaloFinder, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkFriendsOfFriendsHaloFinder);
vtkCxxSetObjectMacro(vtkFriendsOfFriendsHaloFinder,Controller,
	vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::vtkFriendsOfFriendsHaloFinder()
{
	this->SetInputArrayToProcess(
    0,
    0,
    0,
    vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
    vtkDataSetAttributes::SCALARS);
  this->LinkingLength = 1e-6; //default
	this->MinimumNumberOfParticles = 50; // default
	this->Controller = NULL;
	this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::~vtkFriendsOfFriendsHaloFinder()
{
	this->SetController(NULL);
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
vtkIdType vtkFriendsOfFriendsHaloFinder::GetUniqueId(
	int index, vtkIdTypeArray* globalIdArray)
{
	if(RunInParallel(this->GetController()))
		{
		return globalIdArray->GetValue(index)+1;
		}
	else
		{
		return index+1;
		}
}		

//----------------------------------------------------------------------------
vtkIdTypeArray* vtkFriendsOfFriendsHaloFinder::FindHaloes(
	vtkKdTree* pointTree, vtkIdTypeArray* globalIdArray, vtkPointSet* input)
{
	if(this->MinimumNumberOfParticles < 2)
		{
		vtkWarningMacro("setting minimum number of particles to 2, a minimum number of particles below this makes no sense.");
		this->MinimumNumberOfParticles=2;
		}
	// calculating the initial haloes
	// TODO: manage memory
	vtkIdTypeArray* haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");
	// Now assign halos, if this point has at least one other pair,
	// it is a halo, if not it is not (set to 0)
	// first building map of id to count of that id, O(N)
	vtkstd::map<vtkIdType,int> haloCount;
	vtkstd::map<vtkIdType,int> haloUniqueId;	
	// Will only use the following quantities if running in parallel
	vtkstd::map<vtkIdType,int> isHaloSpitAcrossProcessors;
	vtkSmartPointer<vtkPointSet> ghostPoints = \
		vtkSmartPointer<vtkPolyData>::New();
		ghostPoints->Initialize();
	vtkSmartPointer<vtkIdTypeArray> ghostPointLocalHaloIdArray=\
		vtkSmartPointer<vtkIdTypeArray>::New();
		ghostPointLocalHaloIdArray->Initialize();
		ghostPointLocalHaloIdArray->SetName("local halo id");
	for(int nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		if(RunInParallel(this->GetController()))
			{
			// TODO: remove, just testing ghost level functionality
			// TODO:
			// this is where we check if the point for which the nextHaloIdIndex
			// is recorded is a ghost cell. 
			//if it is:
			// o set isHaloSpitAcrossProcessors[haloId] = 1
			// o add it to the ghost point set
			// ghostPoints->SetNextPoint(nextHaloIdIndex,output->getPoint(nextHaloIdIndex))
			// record this ghostPoint's local haloID in  ghostPointLocalHaloIdArray
			if(input->GetPointData()->GetArray("vtkGhostLevels")->GetTuple(
					nextHaloIdIndex)[0]==1)
				{
				cout << "id " << nextHaloIdIndex << " has ghost level 1 ";
				isHaloSpitAcrossProcessors[haloId] = 1;
				ghostPoints->GetPoints()->InsertNextPoint(
					input->GetPoint(nextHaloIdIndex));
				ghostPointLocalHaloIdArray->InsertNextValue(haloId);
				}
			}
		haloCount[haloId]+=1;
		}

	if(RunInParallel(this->GetController()))
		{
		ghostPoints->GetPointData()->AddArray(ghostPointLocalHaloIdArray);
		int procId=this->GetController()->GetLocalProcessId();
		int numProc=this->GetController()->GetNumberOfProcesses();
		if(procId!=0)
			{
			// Sending	ghostPoints to root 
			this->GetController()->Send(ghostPoints,0,
				GHOST_POINTS_AND_LOCAL_HALO_IDS);
			// Wait for roots response with same point set, only containing
			// additional "global halo id" array that root has assigned
			ghostPoints->Initialize();
			this->GetController()->Receive(ghostPoints,0,
				GHOST_POINTS_AND_LOCAL_HALO_IDS_TO_GLOBAL);
			}
		else
			{
			// if we are root:
			vtkSmartPointer<vtkKdTree> ghostPointTree = \
			 	vtkSmartPointer<vtkKdTree>::New();
			vtkPointSet** allGhostPointSetArrays = new vtkPointSet*[numProc];
			vtkPoints** allGhostPointArrays = new vtkPoints*[numProc];
			// first add process 0's ghost points
			allGhostPointSetArrays[0] = ghostPoints;
			allGhostPointArrays[0] = ghostPoints->GetPoints();
			// Receiving ghostPoints from all processors if we are root
			// adding these data sets to pointTree one by one
			for(int proc = 1; proc < numProc; ++proc)
				{
				// this memory will be managed below
				vtkPointSet* recGhostPointSet = vtkPolyData::New();
				recGhostPointSet->Initialize();
				this->GetController()->Receive(recGhostPointSet,proc,
					GHOST_POINTS_AND_LOCAL_HALO_IDS);
				allGhostPointSetArrays[numProc] = recGhostPointSet;
				allGhostPointArrays[numProc] = recGhostPointSet->GetPoints();
				}
			// building a locator from all the points we have received
 		  ghostPointTree->BuildLocatorFromPoints(allGhostPointArrays,numProc);
			// merging these point ids within the linking length across processors
			vtkIdTypeArray* mergeGhostPointHalosIdArray = \
				pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
			// TODO:
			// Calculating a global ID for each of the points, as below

			// Sending result set to the process that root received pointset from
			for(int proc = 1; proc < numProc; ++proc)
				{
				vtkPointSet* sendGhostPointSet = allGhostPointSetArrays[proc];
				this->GetController()->Send(sendGhostPointSet,proc,
				GHOST_POINTS_AND_LOCAL_HALO_IDS_TO_GLOBAL);
				sendGhostPointSet->Delete();
				}
			}
		}

	// finally setting to zero points which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id
	// process 0 or running in serial
	for(int nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		// TODO: check this
		/*
		if(RunInParallel(this->GetController()) && \
			input->GetPointData()->GetArray("vtkGhostLevels")->GetTuple(
				haloId)[0]==1)
			{
			// This point is a ghost point, ignore in output
			continue;
			}
		*/
		if(isHaloSpitAcrossProcessors[haloId])
			{
			// TODO:
			// we use the unique id sent to us by the root process, instead of
			// one that we calculate on our own
			}
		else if(haloCount[haloId]<this->GetMinimumNumberOfParticles())
			{
			// we only saw it less than requisite number of times
			haloIdArray->SetValue(nextHaloIdIndex,0);
			}
		else
			{
			// we saw it requisite number of times, getting global id
			int uniqueId = haloUniqueId[haloId];
			if(!uniqueId)
				{
				haloUniqueId[haloId] = \
				 	this->GetUniqueId(nextHaloIdIndex,globalIdArray);
				}
			haloIdArray->SetValue(nextHaloIdIndex,uniqueId);
			}
		}
	return haloIdArray;
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkPointSet* output = vtkPointSet::GetData(outputVector);
	output->ShallowCopy(input);
	// Get name of data array containing global ids
	// Simply ignored if we are not running in parallel
	vtkSmartPointer<vtkIdTypeArray> globalIdArray = \
		vtkSmartPointer<vtkIdTypeArray>::New();
	if (RunInParallel(this->GetController()))  
		{
		// Requesting one level of ghost cells
		inInfo->Set(
			vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),1);
		// TODO: also assert here that we have a ghost cell array
		vtkSmartPointer<vtkDataArray> globalIdArrayGeneric = \
		 	this->GetInputArrayToProcess(0, inputVector);
		if(!globalIdArray || !globalIdArray->IsA("vtkIdTypeArray"))
	    {
	    vtkErrorMacro("Failed to locate global ID array, this is required if running in parallel. Generate by using Tipsy Reader to read in data, by running D3 with ghost cell generation, or by loading in with the original data in your preferred reader.");
	    return 0;
			}
		else
			{
			globalIdArray = vtkIdTypeArray::SafeDownCast(globalIdArrayGeneric);
			}
		}
	// Building a local KdTree for locator purposes
	vtkSmartPointer<vtkKdTree> pointTree = vtkSmartPointer<vtkKdTree>::New();
	// building a locator
	// changing back as input didn't find my ghost cells
	pointTree->BuildLocatorFromPoints(output);	
	vtkIdTypeArray* haloIdArray = \
		this->FindHaloes(pointTree,globalIdArray,output);
		
	// TODO add back in
	// output->GetPointData()->AddArray(haloIdArray);
	// TODO: manage memory
	// TODO: haloIdArray may be longer than number of points in output
	// if output doesn't copy ghost cells as I expect, should remove
	// ghost point
  return 1;
}
