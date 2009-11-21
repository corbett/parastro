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
	// Initializing maps and sets to keep track of haloes, unique ids, and
	// if running in parallel ghost cells, whether a halo is split accross 
	// processors, and a map from the local halo id to the global
	vtkstd::map<vtkIdType,int> haloCount;
	vtkstd::map<vtkIdType,int> haloUniqueId;	
	// Will only use the following quantities if running in parallel
	vtkstd::map<vtkIdType,int> isHaloSpitAcrossProcessors;
	vtkstd::map<vtkIdType,vtkIdType> ghostHaloLocalToGlobalHaloIds;
	vtkSmartPointer<vtkPointSet> ghostPoints = \
		vtkSmartPointer<vtkPolyData>::New();
		ghostPoints->SetPoints(vtkSmartPointer<vtkPoints>::New());
	vtkSmartPointer<vtkIdTypeArray> ghostPointLocalHaloIdArray=\
		vtkSmartPointer<vtkIdTypeArray>::New();
		ghostPointLocalHaloIdArray->Initialize();
		ghostPointLocalHaloIdArray->SetName("local halo id");

	// Calculating the initial haloes
	vtkIdTypeArray* haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");
	for(int nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		haloCount[haloId]+=1;
		if(RunInParallel(this->GetController()))
			{
			// If we are running in parallel and 
			// dealing with a ghost cell then consider this whole halo is considered
			// to be split across processors, thus the root process instead
			// of this process will be assigning it its unique id
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
		}
	// We are done counting the haloes locally, but if we are running in 
	// parallel we need to consider and merge if necessary the subparts of
	// haloes split accross processes. We do this here
	if(RunInParallel(this->GetController()))
		{
		int procId=this->GetController()->GetLocalProcessId();
		int numProc=this->GetController()->GetNumberOfProcesses();
		// Initializing/adding the arrays keeping track of the local and global
		// halo ids as well as the total halo count. 
		ghostPoints->GetPointData()->AddArray(ghostPointLocalHaloIdArray);
		AllocateIdTypeDataArray(ghostPoints,"global halo id",1,
			ghostPoints->GetPoints()->GetNumberOfPoints());
		AllocateIdTypeDataArray(ghostPoints,"global halo count",1,
			ghostPoints->GetPoints()->GetNumberOfPoints());

		// updating each ghost point with the count of its elements
		for(int i = 0; 
			i < ghostPoints->GetPoints()->GetNumberOfPoints(); ++i)
			{
			vtkIdType localGhostId = vtkIdTypeArray::SafeDownCast(
				ghostPoints->GetPointData()->GetArray("local halo id"))->GetValue(i);
			vtkIdTypeArray::SafeDownCast(ghostPoints->GetPointData()->GetArray(
				"global halo count"))->SetValue(i,haloCount[localGhostId]);
			}
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
			// We are not root
			vtkSmartPointer<vtkKdTree> ghostPointTree = \
			 	vtkSmartPointer<vtkKdTree>::New();
			ghostPointTree->Initialize();
			vtkPointSet** allGhostPointSetArrays = new vtkPointSet*[numProc];
			vtkPoints** allGhostPointArrays = new vtkPoints*[numProc];
			// first add process 0's ghost points
			allGhostPointSetArrays[0] = ghostPoints;
			allGhostPointArrays[0] = ghostPoints->GetPoints();
			// Receiving ghostPoints from all other processors,
			// adding these data sets to pointTree one by one
			int totalNumberOfPoints = 0;
			for(int proc = 1; proc < numProc; ++proc)
				{
				// this memory will be managed below
				vtkPointSet* recGhostPointSet = vtkPolyData::New();
				recGhostPointSet->Initialize();
				this->GetController()->Receive(recGhostPointSet,proc,
					GHOST_POINTS_AND_LOCAL_HALO_IDS);
				allGhostPointSetArrays[proc] = recGhostPointSet;
				allGhostPointArrays[proc] = recGhostPointSet->GetPoints();
				totalNumberOfPoints += \
				 	recGhostPointSet->GetPoints()->GetNumberOfPoints();
				}
			if(totalNumberOfPoints>0)
				{
				// building a locator from all the points we have received
	 		  ghostPointTree->BuildLocatorFromPoints(allGhostPointArrays,numProc);
				// merging these point ids within the linking length across processors
				vtkIdTypeArray* mergeGhostPointHalosIdArray = \
					pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
				// TODO: actually count results here and assign global ids
				// TODO
				// TODO
				}
			// Sending result set to the process that root received pointset from
			for(int proc = 1; proc < numProc; ++proc)
				{
				vtkPointSet* sendGhostPointSet = allGhostPointSetArrays[proc];
				this->GetController()->Send(sendGhostPointSet,proc,
				GHOST_POINTS_AND_LOCAL_HALO_IDS_TO_GLOBAL);
				sendGhostPointSet->Delete();
				}
			}
		// Both root and other processes need to do this step, mapping local
		// to global ids based on the ghostPoint's array data assigned by root.
		for(int i = 0; 
			i < ghostPoints->GetPoints()->GetNumberOfPoints(); ++i)
			{
			vtkIdType localId = \
				vtkIdTypeArray::SafeDownCast(
				ghostPoints->GetPointData()->GetArray("local halo id"))->GetValue(i);
			vtkIdType globalId = \
				vtkIdTypeArray::SafeDownCast(
				ghostPoints->GetPointData()->GetArray("global halo id"))->GetValue(i);
			int globalHaloCount = \
				vtkIntArray::SafeDownCast(ghostPoints->GetPointData()->GetArray(
				"global halo count"))->GetValue(i);
			ghostHaloLocalToGlobalHaloIds[localId]=globalId;
			haloCount[localId] = globalHaloCount;
			}
		// Done with the parallel part.
		}

	// Finally, whether in serial or parallel, setting to zero points which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id.
	for(int nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		// Switching to the global ID if this halo is split accross processsors
		haloId = (isHaloSpitAcrossProcessors[haloId]) ? \
		 	ghostHaloLocalToGlobalHaloIds[haloId] : haloId;
		if(haloCount[haloId]<this->GetMinimumNumberOfParticles())
			{
			// we only saw it less than requisite number of times
			haloIdArray->SetValue(nextHaloIdIndex,0);
			}
		else
			{
			// we saw it requisite number of times, getting unique id if it 
			// has already been assigned, otherwise assigning it.
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
	output->GetPointData()->AddArray(haloIdArray);
	// Managing ßmemory
	haloIdArray->Delete();
  return 1;
}
