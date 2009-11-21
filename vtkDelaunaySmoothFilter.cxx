/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkDelaunaySmoothFilter.cxx,v $
=========================================================================*/
#include "vtkDelaunaySmoothFilter.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkMultiProcessController.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkGenericPointIterator.h"
#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkCallbackCommand.h"
#include "vtkPointLocator.h"
#include "vtkDelaunay3D.h"

using vtkstd::string;
vtkCxxRevisionMacro(vtkDelaunaySmoothFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkDelaunaySmoothFilter);
//----------------------------------------------------------------------------
vtkDelaunaySmoothFilter::vtkDelaunaySmoothFilter():vtkPointSetAlgorithm()
{
	this->SetInputArrayToProcess(
    0,
    0,
    0,
    vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
    vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
vtkDelaunaySmoothFilter::~vtkDelaunaySmoothFilter()
{
}

//----------------------------------------------------------------------------
void vtkDelaunaySmoothFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkDelaunaySmoothFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkDelaunaySmoothFilter::RequestData(vtkInformation *request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// smoothing each quantity in the output, we will need a record of the
	// original number of input arrays
	int numberOriginalArrays = input->GetPointData()->GetNumberOfArrays();
	vtkDataArray* massArray = this->GetInputArrayToProcess(0, inputVector);
  if (!massArray)
    {
    vtkErrorMacro("Failed to locate mass array");
    return 0;
    }
  vtkPointSet* output = vtkPointSet::GetData(outputVector);
  output->ShallowCopy(input);
	// 1. Allocating arrays to store our smoothed values
 	AllocateDoubleDataArray(output,"smoothed density", 
		1,output->GetPoints()->GetNumberOfPoints());
	for(int i = 0; i < numberOriginalArrays; ++i)
		{
		vtkSmartPointer<vtkDataArray> nextArray = \
		 	output->GetPointData()->GetArray(i);
		string baseName = nextArray->GetName();
		for(int comp = 0; comp < nextArray->GetNumberOfComponents(); ++comp)
			{
			string totalName = GetSmoothedArrayName(baseName,comp);
			// Allocating an column for the total sum of the existing quantities
			AllocateDoubleDataArray(output,totalName.c_str(),
				1,output->GetPoints()->GetNumberOfPoints());
			}
		}
	// 2. Performing the Delaunay tessellation, local to the process. Heart
	// of the filter here.
  vtkSmartPointer<vtkDelaunay3D> delaunay = \
 		vtkSmartPointer<vtkDelaunay3D>::New();
  delaunay->AddInput(output);
  delaunay->Update();
	// 3. For each point, smoothing over delaunay cells neighboring that point's
	// delaunay cell. The smoothed density is the smoothed mass of the point
	// divided by the volume of the delaunay cell in which it resides.
	for(int nextPointId = 0;
		nextPointId < output->GetPoints()->GetNumberOfPoints();
	 	++nextPointId)
		{
		double* nextPoint=GetPoint(output,nextPointId);
		// 1. Finding the volume in delaunay cell of this point
		// TODO:, do this
		double delaunayCellVolume = 0;
		// 2. finding the points in neighboring cells
		// TODO: do this
		vtkSmartPointer<vtkIdList> neighborCellPoints = \
			vtkSmartPointer<vtkIdList>::New();
		// summing quantities of neighbor delaunay cell points for smoothing, 
		// only if we have more neighbors than ourselves
		if(neighborCellPoints->GetNumberOfIds()>0)
			{
			for(int neighborPointLocalId = 0;
		 		neighborPointLocalId < neighborCellPoints->GetNumberOfIds();
				++neighborPointLocalId)
				{
				vtkIdType neighborPointGlobalId = \
										neighborCellPoints->GetId(neighborPointLocalId);
				double* neighborPoint=GetPoint(output,neighborPointGlobalId);
			// keeps track of the totals for each quantity inside
			// the output, only dividing by N at the end
			for(int i = 0; i < numberOriginalArrays; ++i)
				{
				vtkSmartPointer<vtkDataArray> nextArray = \
					output->GetPointData()->GetArray(i);
				string baseName = nextArray->GetName();
				double* data=GetDataValue(output,baseName.c_str(),
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
			double numberPoints = neighborCellPoints->GetNumberOfIds();
			for(int i = 0; i < numberOriginalArrays; ++i)
				{
				vtkSmartPointer<vtkDataArray> nextArray = \
					output->GetPointData()->GetArray(i);
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
			double* smoothedMass = \
				GetDataValue(output,GetSmoothedArrayName(massArray->GetName(),
				0).c_str(),
				nextPointId);
			double smoothedDensity=smoothedMass[0]/delaunayCellVolume;
			//storing the smooth density
			SetDataValue(output,"smoothed density",nextPointId,&smoothedDensity);
			// Finally, some memory management
			delete [] smoothedMass;
			}
		else
			{
			// This point has no neighbors, so smoothed mass is identicle to 
			// this point's mass, and smoothed density is meaningless, set to -1
			// to indicate it is useless
			double density=-1;
			SetDataValue(output,"smoothed density",nextPointId,&density);
			}
		// Finally, some memory management
		delete [] nextPoint;
		}
  return 1;
}

string vtkDelaunaySmoothFilter::GetSmoothedArrayName(string baseName, 
	int comp)
{
	return "smoothed_"+baseName+"_"+ToString(comp);
}
