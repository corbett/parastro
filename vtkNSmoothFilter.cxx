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
#include "astrovizhelpers/DataSetHelpers.h"

vtkCxxRevisionMacro(vtkNSmoothFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkNSmoothFilter);

//----------------------------------------------------------------------------
vtkNSmoothFilter::vtkNSmoothFilter()
{
  this->NeighborNumber = 50; //default
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
	double smoothedDensity;
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
	// Building the Kd tree
  vtkDebugMacro("1a. Building Kd tree.");
	vtkSmartPointer<vtkPKdTree> pointTree = vtkSmartPointer<vtkPKdTree>::New();
	pointTree->BuildLocatorFromPoints(output);
	vtkDebugMacro("1b. Allocating arrays to store our smoothed values.");
	// allocating an arrays for each of our smoothed values
 	AllocateDoubleDataArray(output,"smoothed mass", \
										1,output->GetPoints()->GetNumberOfPoints());
	AllocateDoubleDataArray(output,"smoothed density", \
										1,output->GetPoints()->GetNumberOfPoints());
	/*
	//TODO: add these later
	AllocateDoubleDataArray(output,"smoothed velocity", \
		3,output->GetPoints()->GetNumberOfPoints());
	AllocateDoubleDataArray(output,"smoothed  speed",\
		1,output->GetPoints()->GetNumberOfPoints());
	*/
  vtkDebugMacro("2. Calculating the smoothed quantities \
									we are interested in.");
	for(int nextPointId = 0;\
	 		nextPointId < output->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
		double* nextPoint=GetPoint(output,nextPointId);
		vtkDebugMacro("next point is " << nextPoint[0] << "," \
									<< nextPoint[1] << ","<< nextPoint[2]);
		// finding the closest N points
		vtkDebugMacro("2. Finding the closeset N points");
		vtkSmartPointer<vtkIdList> closestNPoints = \
			vtkSmartPointer<vtkIdList>::New();
		pointTree->FindClosestNPoints(this->NeighborNumber,\
																	nextPoint,closestNPoints);
		vtkDebugMacro("found " << closestNPoints->GetNumberOfIds() \
									<< " closest points");
		// looping over the closestNPoints, 
		// only if we have more neighbors than ourselves
		if(closestNPoints->GetNumberOfIds()>0)
			{
				double totalMass;
				//we do want to count this point's mass 
				double* smoothedMass=new double[1];
				double* smoothedDensity=new double[1];
				for(int neighborPointLocalId = 0; \
			 		neighborPointLocalId < closestNPoints->GetNumberOfIds()-1; \
					++neighborPointLocalId)
				{
				vtkIdType neighborPointGlobalId = \
										closestNPoints->GetId(neighborPointLocalId);
				double* neighborPoint=GetPoint(output,neighborPointGlobalId);
				vtkDebugMacro("the " << neighborPointLocalId \
											<< "th nearest point coordiates are (" \
				 							<< neighborPoint[0] << "," << neighborPoint[1] \
				 							<< "," << neighborPoint[2] << ")");
				// extracting the mass
				// has to be double as this version of VTK doesn't have 
				double* mass=GetDataValue(output,"mass",neighborPointGlobalId);
				totalMass+=mass[0];
				// Finally, some memory management
				delete [] mass;
				delete [] neighborPoint;
				}
			// storing the smoothed mass in the output vector
			smoothedMass[0]=totalMass/(closestNPoints->GetNumberOfIds());
			SetDataValue(output,"smoothed mass",nextPointId,smoothedMass);
			//for the smoothed Density we need the identity of the 
			// last neighbor point, as this is farthest from the original point
			vtkIdType lastNeighborPointGlobalId = \
									closestNPoints->GetId(closestNPoints->GetNumberOfIds()-1);
			double* lastNeighborPoint=GetPoint(output,lastNeighborPointGlobalId);		
			smoothedDensity[0]=\
			 				CalculateDensity(nextPoint,lastNeighborPoint,smoothedMass[0]);
			vtkDebugMacro("smoothed density is " << smoothedDensity[0]); 
			//storing the smooth density
			SetDataValue(output,"smoothed density",nextPointId,smoothedDensity);
			// Finally, some memory management
			delete [] lastNeighborPoint;
			delete [] smoothedDensity;
			delete [] smoothedMass;
			}
		else
			{
			// This point has no neighbors, so smoothed mass is identicle to 
			// this point's mass, and smoothed density is meaningless, set to -1
			// to indicate it is useless
			double* density=new double[1];
			density[0]=-1;
			double* mass=GetDataValue(output,"mass",nextPointId);
			SetDataValue(output,"smoothed mass",nextPointId,mass);
			SetDataValue(output,"smoothed density",nextPointId,density);
			// Finally, some memory management
			delete [] density;
			delete [] mass;
			}
		// Finally, some memory management
		delete [] nextPoint;
		}
	// Finally, some memory management
  output->Squeeze();
  return 1;
}
