/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNSmoothFilter.cxx,v $
=========================================================================*/
#include "vtkNSmoothFilter.h"
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
#include "astrovizhelpers/DataSetHelpers.h"

using vtkstd::string;

vtkCxxRevisionMacro(vtkNSmoothFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkNSmoothFilter);
vtkCxxSetObjectMacro(vtkNSmoothFilter,Controller,vtkMultiProcessController);
vtkCxxSetObjectMacro(vtkNSmoothFilter,PKdTree,vtkPKdTree);
vtkCxxSetObjectMacro(vtkNSmoothFilter,D3,vtkDistributedDataFilter);



//----------------------------------------------------------------------------
vtkNSmoothFilter::vtkNSmoothFilter()
{
  this->NeighborNumber = 50; //default
  this->Controller = NULL;
  this->PKdTree = NULL;
  this->D3 = NULL;

  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkNSmoothFilter::~vtkNSmoothFilter()
{
  this->SetPKdTree(NULL);
  this->SetController(NULL);
  this->SetD3(NULL);
}

//----------------------------------------------------------------------------
void vtkNSmoothFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Neighbor Number: " << this->NeighborNumber << "\n";
}

//----------------------------------------------------------------------------
int vtkNSmoothFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
double vtkNSmoothFilter::CalculateDensity(double pointOne[],\
													double pointTwo[], double smoothedMass)
{
	// now calculating the radial distance from the last point to the
  // center point to which it is a neighbor
	double radialDistance=sqrt(vtkMath::Distance2BetweenPoints(pointOne,\
																														pointTwo));
	// the volume is a sphere around nextPoint with radius of the 
	// last in the list of the closestNpoints
	// so 4/3 pi r^3 where 
	double neighborhoodVolume=4./3 * M_PI * pow(radialDistance,3);
	double smoothedDensity=0.0;
	if(neighborhoodVolume!=0)
		{
		smoothedDensity=smoothedMass/neighborhoodVolume;
		}
	return smoothedDensity;
}

//----------------------------------------------------------------------------
int vtkNSmoothFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet* output = vtkPointSet::GetData(outputVector);
  // Outline of this filter:
	// 1. Build Kd tree
	// 2. Go through each point in output
	// 		o calculate N nearest neighbors
	//		o calculate smoothed quantities
	// 		o add to their respective arrays.
	// 3. Add the arrays of smoothed varaibles to the output
  // First copying the input to the output, as the output will be 
	// identical to the input, with some additional properties
	// copies the point positions
  output->CopyStructure(input);
	// copies the point attributes
  output->CopyAttributes(input);
	// TODO: remove-just seeing if I understand how this filter accepts input
	if (this->D3 == NULL)
    {
		cout << "NO D3!!\n";
	  this->D3 = vtkDistributedDataFilter::New();
    }
  else
		{
			cout<<"D3!!\n";
		}
  this->D3->SetBoundaryModeToSplitBoundaryCells();
  this->D3->SetInput(input);
  this->D3->SetController(this->Controller);
  this->D3->Modified();
  this->D3->Update();

	if (this->PKdTree == NULL)
    {
		this->PKdTree = this->D3->GetKdtree();
		this->PKdTree->SetController(this->Controller);
		cout << "now Pkdtree\n";
    }
  else
		{
			cout<<"alread Pkdtree!!\n";
		}
	// Building the Kd tree, should already be built
	// TODO: remove
//	vtkSmartPointer<vtkPkdtree> pointTree = vtkSmartPointer<vtkPkdtree>::New();
	cout << "building locator from input \n";
	this->PKdTree->BuildLocatorFromPoints(input);
	cout << " done with locator \n";
	// Allocating arrays to store our smoothed values
	// smoothed density
 	AllocateDoubleDataArray(output,"smoothed density", 
		1,output->GetPoints()->GetNumberOfPoints());
	// smoothing each quantity in the input
	for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
		{
		vtkSmartPointer<vtkDataArray> nextArray = \
		 	input->GetPointData()->GetArray(i);
		string baseName = nextArray->GetName();
		for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
			{
			string totalName = GetSmoothedArrayName(baseName,comp);
			// Allocating an column for the total sum of the existing quantities
			AllocateDoubleDataArray(output,totalName.c_str(),
				1,output->GetPoints()->GetNumberOfPoints());
			}
		}
	for(int nextPointId = 0;
		nextPointId < input->GetPoints()->GetNumberOfPoints();
	 	++nextPointId)
		{
		double* nextPoint=GetPoint(input,nextPointId);
		// finding the closest N points
		vtkSmartPointer<vtkIdList> closestNPoints = \
			vtkSmartPointer<vtkIdList>::New();
		// plus one as the first point returned by locator is always one's self, 
		// and the user expects specifying 1 neighbor will actually find
		// one neighbor 
		this->PKdTree->FindClosestNPoints(this->NeighborNumber+1,
																	nextPoint,closestNPoints);
		// looping over the closestNPoints, 
		// only if we have more neighbors than ourselves
		if(closestNPoints->GetNumberOfIds()>0)
			{
			for(int neighborPointLocalId = 0;
		 		neighborPointLocalId < closestNPoints->GetNumberOfIds();
				++neighborPointLocalId)
				{
				vtkIdType neighborPointGlobalId = \
										closestNPoints->GetId(neighborPointLocalId);
				double* neighborPoint=GetPoint(input,neighborPointGlobalId);
			// keeps track of the totals for each quantity inside
			// the output, only dividing by N at the end
				for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
					{
					vtkSmartPointer<vtkDataArray> nextArray = \
						input->GetPointData()->GetArray(i);
					string baseName = nextArray->GetName();
					double* data=GetDataValue(input,baseName.c_str(),
						neighborPointGlobalId);
					for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
						{
						string totalName = GetSmoothedArrayName(baseName,comp);
						// adds data[comp] to value in column totalName
						double* total = \
							GetDataValue(output,totalName.c_str(),nextPointId);
						total[0]+=data[comp];
						SetDataValue(output,totalName.c_str(),
							nextPointId,total);
						// memory management
						delete [] total;
						}
					delete [] data;
					}
				// Finally, some memory management
				delete [] neighborPoint;
				}
			// dividing by N at the end
			double numberPoints = closestNPoints->GetNumberOfIds();
			for(int i = 0; i < input->GetPointData()->GetNumberOfArrays(); ++i)
				{
				vtkSmartPointer<vtkDataArray> nextArray = \
					input->GetPointData()->GetArray(i);
				string baseName = nextArray->GetName();
				for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
					{
					string totalName = GetSmoothedArrayName(baseName,comp);
					// adds data[comp] to value in column totalName
					double* smoothedData = \
						GetDataValue(output,totalName.c_str(),nextPointId);
					smoothedData[0]/=numberPoints;
					SetDataValue(output,totalName.c_str(),nextPointId,smoothedData);
					// memory management
					delete [] smoothedData;
					}
				}
			// for the smoothed Density we need the identity of the 
			// last neighbor point, as this is farthest from the original point
			// we use this to calculate the volume over which to smooth
			vtkIdType lastNeighborPointGlobalId = \
									closestNPoints->GetId(closestNPoints->GetNumberOfIds()-1);
			double* lastNeighborPoint=GetPoint(output,lastNeighborPointGlobalId);		
			// there MUST be an array named mass in the input.
			double* smoothedMass = \
				GetDataValue(output,GetSmoothedArrayName("mass",0).c_str(),
				nextPointId);
			double* smoothedDensity=new double[1];
			smoothedDensity[0]=\
			 	CalculateDensity(nextPoint,lastNeighborPoint,smoothedMass[0]);
			//storing the smooth density
			SetDataValue(output,"smoothed density",nextPointId,smoothedDensity);
			// Finally, some memory management
			delete [] lastNeighborPoint;
			delete [] smoothedMass;
			delete [] smoothedDensity;
			}
		else
			{
			// This point has no neighbors, so smoothed mass is identicle to 
			// this point's mass, and smoothed density is meaningless, set to -1
			// to indicate it is useless
			double* density=new double[1];
			density[0]=-1;
			SetDataValue(output,"smoothed density",nextPointId,density);
			// Finally, some memory management
			delete [] density;
			}
		// Finally, some memory management
		delete [] nextPoint;
		}
	// Finally, some memory management
  output->Squeeze();
  return 1;
}

string vtkNSmoothFilter::GetSmoothedArrayName(string baseName, int comp){
	return "smoothed_"+baseName+"_"+ToString(comp);
}
