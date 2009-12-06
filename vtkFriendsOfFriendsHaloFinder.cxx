/*=========================================================================

  Program:   AstroViz plugin for ParaView
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
	unsigned long index, vtkIdTypeArray* globalIdArray)
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
	vtkstd::map<vtkIdType,unsigned long> haloCount;
	vtkstd::map<vtkIdType,unsigned long> haloUniqueId;	
	// 1.  Calculating the initial haloes
	vtkIdTypeArray* haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");
	for(unsigned long nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		haloCount[haloId]+=1;
		}
	// Setting proto-halos which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id, meaning we have truly found a halo.
	for(unsigned long nextHaloIdIndex = 0;
		nextHaloIdIndex < haloIdArray->GetNumberOfTuples();
	 	++nextHaloIdIndex)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloIdIndex);
		if(haloCount[haloId]<this->GetMinimumNumberOfParticles())
			{
			// we only saw it less than requisite number of times
			haloIdArray->SetValue(nextHaloIdIndex,0);
			}
		else
			{
			// we saw it requisite number of times, getting unique id if it 
			// has already been assigned, otherwise assigning it.
			unsigned long uniqueId = haloUniqueId[haloId];
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
	// Outline of this filter:
	// 1. Build Kd tree
	// 2. Go through each point in output
	// 		o calculate points within linking length
	// 		o these form a halo
	// 3. Cutoff by particle count
	// 		o if proto-halo doesn't have minimum particle count, 
	// 		it is not considered a halo. if it does, it is given a unique id.
	
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	// TODO: trying some ghost level stuff
	// Requesting one level of ghost cells
	inInfo->Set(
		vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),1);
	vtkPointSet* output = vtkPointSet::GetData(outputVector);
	output->ShallowCopy(input);
	// Get name of data array containing global ids
	// Simply ignored if we are not running in parallel
	vtkSmartPointer<vtkIdTypeArray> globalIdArray = \
		vtkSmartPointer<vtkIdTypeArray>::New();
	if (RunInParallel(this->GetController()))  
		{
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
	// Managing memory
	haloIdArray->Delete();
  return 1;
}
