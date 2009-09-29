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
#include "vtkKdTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkGenericPointIterator.h"

vtkCxxRevisionMacro(vtkNSmoothFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkNSmoothFilter);

//----------------------------------------------------------------------------
vtkNSmoothFilter::vtkNSmoothFilter()
{
  this->NeighborNumber = 50; //default
//  this->GetInformation()->Set(vtkAlgorithm::PRESERVES_RANGES(), 1);
// this->GetInformation()->Set(vtkAlgorithm::PRESERVES_BOUNDS(), 1);
}

//----------------------------------------------------------------------------
vtkNSmoothFilter::~vtkNSmoothFilter()
{
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
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkNSmoothFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  // Outline of this filter:
	// 1. Build Kd tree
	// 2. Go through each point in output
	// 		o calculate N nearest neighbors
	//		o calculate smoothed quantities
	// 		o add to their respective arrays.
	// 3. Add the arrays of smoothed varaibles to the output
  // First copying the input to the output, as the output will be identical to the input, with some additional properties
	// copies the point positions
  output->CopyStructure(input);
	// copies the point attributes
  output->CopyAttributes(input);
	// Building the Kd tree
  vtkDebugMacro("1. Building Kd tree.");
	vtkSmartPointer<vtkKdTree> pointTree = vtkSmartPointer<vtkKdTree>::New();
		pointTree->BuildLocatorFromPoints(output);
  vtkDebugMacro("2. Calculating the smoothed quantities we are interested in.");
	for(int id = 0; id < output->GetPoints()->GetNumberOfPoints(); ++id)
		{
		double nextPoint[3]; 
		output->GetPoints()->GetPoint(id,nextPoint);
		vtkDebugMacro("next point is " << nextPoint[0] << ","<< nextPoint[1] << ","<< nextPoint[2]);
		//finding the closest N points
		vtkDebugMacro("2. Finding the closeset N points");
		vtkSmartPointer<vtkIdList> closestNPoints = vtkSmartPointer<vtkIdList>::New();
		pointTree->FindClosestNPoints(this->NeighborNumber,nextPoint,closestNPoints);
		vtkErrorMacro("found " << closestNPoints->GetNumberOfIds() << " closest points");
		//looping over the closestNPoints
		for(int j = 0; j < closestNPoints->GetNumberOfIds(); ++j)
			{
			double neighborPoint[3];
			vtkIdType neighborPointId = closestNPoints->GetId(j);
			output->GetPoints()->GetPoint(neighborPointId,neighborPoint);
			vtkErrorMacro("the " << j <<"th nearest point coordiates are (" << neighborPoint[0] << "," << neighborPoint[1] << "," << neighborPoint[2] << ")");
			//extracting the mass
			double mass[1];
			output->GetPointData()->SetActiveScalars("mass");
			output->GetPointData()->GetScalars()->GetTuple(neighborPointId,mass);
			vtkErrorMacro("it has mass " << mass[0]);
			}
		//finding the average of each property we are interested in by dividing by #closestNPoints
		//the volume is a sphere around nextPoint with radius of the last in the list of the closestNpoints
		//so 4/3 pi r^3 where r=sqrt((nextPoint->x-nextPoint->x)^2+(nextPoint->y-nextPoint->y)^2+(nextPoint->z-nextPoint->z)^2)
		}
	vtkDebugMacro("3. Storing smoothed quantities in output.");
	// Finally, some memory management
  output->Squeeze();
  // Memory management of the Kd-tree we have built. 
	// TODO: Could later consider outputting this structure for further analysis.
  vtkDebugMacro("4. Smoothing calculation sucessful, deleting Kd tree.");
  return 1;
}
