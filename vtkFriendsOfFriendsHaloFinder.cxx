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
int vtkFriendsOfFriendsHaloFinder::FindHaloes(vtkKdTree* pointTree,
	vtkIdTypeArray* globalIdArray, vtkPointSet* output)
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
	haloIdArray->SetNumberOfTuples(output->GetPoints()->GetNumberOfPoints());
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
	for(int nextHaloId = 0;
		nextHaloId < haloIdArray->GetNumberOfTuples();
	 	++nextHaloId)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloId);
		if(RunInParallel(this->GetController()))
			{
			// TODO: remove, just testing ghost level functionality
			cout << "id " << haloId << " has ghost level " \
			 << output->GetPointData()->GetArray("vtkGhostLevels")->GetTuple(
				nextHaloId)[0]
			 << "\n";
			// TODO:
			// this is where we check if the point for which the nextHaloId
			// is recorded is a ghost cell. 
			//if it is:
			// o set isHaloSpitAcrossProcessors[haloId] = 1
			// o add it to the ghost point set
			// ghostPoints->SetNextPoint(nextHaloId,output->getPoint(nextHaloId))
			// record this ghostPoint's local haloID in  ghostPointLocalHaloIdArray
			}
		haloCount[haloId]+=1;
		}

	if(RunInParallel(this->GetController()))
		{
		// TODO:
		// if we are not root:
			// Sending	ghostPoints to root 
			// waiting for root's response with global ID
			// waiting for root's response with total globalIds already assigned
		// if we are root:
		// Receiving ghostPoints from all processors if we are root
		// Adding these data sets to a point locator, one by one
		// Merging duplicate points in this point locator
		// Calculating a global ID for each of the points, as below

		}

	// finally setting to zero points which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id
	// process 0 or running in serial
	for(int nextHaloId = 0;
		nextHaloId < haloIdArray->GetNumberOfTuples();
	 	++nextHaloId)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloId);	
		if(isHaloSpitAcrossProcessors[haloId])
			{
			// TODO:
			// we use the unique id sent to us by the root process, instead of
			// one that we calculate on our own
			}
		else if(haloCount[haloId]<this->GetMinimumNumberOfParticles())
			{
			// we only saw it less than requisite number of times
			haloIdArray->SetValue(nextHaloId,0);
			}
		else
			{
			// we saw it requisite number of times, getting global id
			int uniqueId = haloUniqueId[haloId];
			if(!uniqueId)
				{
				haloUniqueId[haloId] = this->GetUniqueId(nextHaloId,globalIdArray);
				}
			haloIdArray->SetValue(nextHaloId,uniqueId);
			}
		}
	output->GetPointData()->AddArray(haloIdArray);
	return 1;
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
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
	pointTree->BuildLocatorFromPoints(output);	
	this->FindHaloes(pointTree,globalIdArray,output);
  return 1;
}
