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
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include <cmath>

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
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
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
  // First copying the input to the output, as the output will be identical to the input, with some additional properties
	// copies the point positions
  output->CopyStructure(input);
	// copies the point attributes
  output->CopyAttributes(input);
	// Building the Kd tree
  vtkDebugMacro("1a. Building Kd tree.");
	vtkSmartPointer<vtkPKdTree> pointTree = vtkSmartPointer<vtkPKdTree>::New();
		pointTree->BuildLocatorFromPoints(output);
	vtkDebugMacro("1b. Allocating arrays to store our smoothed values.");
	// allocating an array to store the smoothed mass
	vtkSmartPointer<vtkFloatArray> smoothedMassArray = vtkSmartPointer<vtkFloatArray>::New();
  	smoothedMassArray->SetNumberOfComponents(1);
  	smoothedMassArray->SetName("smoothed mass");
		smoothedMassArray->SetNumberOfTuples(output->GetPoints()->GetNumberOfPoints());
	vtkSmartPointer<vtkFloatArray> smoothedDensityArray = vtkSmartPointer<vtkFloatArray>::New();
  	smoothedDensityArray->SetNumberOfComponents(1);
  	smoothedDensityArray->SetName("smoothed density");
		smoothedDensityArray->SetNumberOfTuples(output->GetPoints()->GetNumberOfPoints());
  vtkDebugMacro("2. Calculating the smoothed quantities we are interested in.");
	for(int id = 0; id < output->GetPoints()->GetNumberOfPoints(); ++id)
		{
		double nextPoint[3]; 
		output->GetPoints()->GetPoint(id,nextPoint);
		vtkDebugMacro("next point is " << nextPoint[0] << ","<< nextPoint[1] << ","<< nextPoint[2]);
		// finding the closest N points
		vtkDebugMacro("2. Finding the closeset N points");
		vtkSmartPointer<vtkIdList> closestNPoints = vtkSmartPointer<vtkIdList>::New();
		pointTree->FindClosestNPoints(this->NeighborNumber,nextPoint,closestNPoints);
		vtkDebugMacro("found " << closestNPoints->GetNumberOfIds() << " closest points");
		// looping over the closestNPoints, only if we have more neighbors than ourselves
		if(closestNPoints->GetNumberOfIds()>0)
			{
			float totalMass=0;
			// closestNPoints is ordered by distance, so the last is the radius we want to calculate the volume with
			// thus, only loop up to the second to last
			for(int j = 0; j < closestNPoints->GetNumberOfIds() -1 ; ++j)
				{
				double neighborPoint[3];
				vtkIdType neighborPointId = closestNPoints->GetId(j);
				output->GetPoints()->GetPoint(neighborPointId,neighborPoint);
				vtkDebugMacro("the " << j <<"th nearest point coordiates are (" << neighborPoint[0] << "," << neighborPoint[1] << "," << neighborPoint[2] << ")");
				// extracting the mass
				// has to be double as this version of VTK doesn't have GetTuple function which operates with float
				double mass[1]; 
				output->GetPointData()->SetActiveScalars("mass");
				output->GetPointData()->GetScalars()->GetTuple(neighborPointId,mass);
				// taking log to help stay off loss of precision
				totalMass+=static_cast<float>(log(mass[0])); 
				}
			//now for the last point
			double neighborPoint[3];
			vtkIdType neighborPointId = closestNPoints->GetId(closestNPoints->GetNumberOfIds()-1);
			output->GetPoints()->GetPoint(neighborPointId,neighborPoint);
			vtkDebugMacro("the " << closestNPoints->GetNumberOfIds()-1 <<"th nearest point coordiates are (" << neighborPoint[0] << "," << neighborPoint[1] << "," << neighborPoint[2] << ")");
			// extracting the mass
			// has to be double as this version of VTK doesn't have GetTuple function which operates with float
			double mass[1]; 
			output->GetPointData()->SetActiveScalars("mass");
			output->GetPointData()->GetScalars()->GetTuple(neighborPointId,mass);
			// taking log to help stay off loss of precision
			totalMass+=log(mass[0]); 
			//done with the last point			
			// storing the smoothed mass in the output vector
			double smoothedMass=totalMass/(closestNPoints->GetNumberOfIds());
			smoothedMassArray->SetValue(neighborPointId,static_cast<float>(smoothedMass));
/*
			// finding the average of each property we are interested in by dividing by #closestNPoints
			// now calculating the radial distance from the last point to the center point to which it is a neighbor
			double radialDistance=sqrt(pow(nextPoint[0]-neighborPoint[0],2)+pow(nextPoint[1]-neighborPoint[1],2)+pow(nextPoint[2]-neighborPoint[2],2));
			// the volume is a sphere around nextPoint with radius of the last in the list of the closestNpoints
			// so 4/3 pi r^3 where r=sqrt((nextPoint->x-nextPoint->x)^2+(nextPoint->y-nextPoint->y)^2+(nextPoint->z-nextPoint->z)^2)
			double neighborhoodVolume=4./3 * M_PI * pow(radialDistance,3);
			double smoothedDensity;
			if(neighborhoodVolume!=0)
				{
				smoothedDensity=smoothedMass/neighborhoodVolume;
				}
			vtkDebugMacro("smoothed density is " << static_cast<float>(smoothedDensity)); 
			smoothedDensityArray->SetValue(neighborPointId,static_cast<float>(smoothedDensity));
*/ 			
			// taking the exp to reverse the log, and also converting back to float precision
			vtkDebugMacro("smoothed mass is " << static_cast<float>(exp(smoothedMass))); 
		  smoothedMassArray->SetValue(neighborPointId, static_cast<float>(exp(smoothedMass))); 
			
			}
		else
			{
				//This point has no neighbors, so smoothed mass is identicle to this point's mass, and smoothed density is meaningless
				double mass[1];
				output->GetPointData()->SetActiveScalars("mass");
				output->GetPointData()->GetScalars()->GetTuple(id,mass);
			  smoothedMassArray->SetValue(id,static_cast<float>(mass[0]));
			}
		}
	vtkDebugMacro("3. Storing smoothed quantities in output.");
	//storing the output vector in the output
	output->GetPointData()->AddArray(smoothedMassArray);
	output->GetPointData()->AddArray(smoothedDensityArray);
	// Finally, some memory management
  output->Squeeze();
  vtkDebugMacro("4. Smoothing calculation sucessful, deleting Kd tree.");
  return 1;
}
