/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkProfileFilter.cxx,v $
=========================================================================*/
#include "vtkProfileFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "vtkCellData.h"
#include "vtkTable.h"
#include "vtkSortDataArray.h"

vtkCxxRevisionMacro(vtkProfileFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkProfileFilter);

//----------------------------------------------------------------------------
vtkProfileFilter::vtkProfileFilter()
{
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
}

//----------------------------------------------------------------------------
void vtkProfileFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkProfileFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkProfileFilter::GetBinInfo(vtkPointSet* input, double* binLowerBound,\
												double* binUpperBound, double* binWidth)
{
	// TODO: implement
	double bounds[6];
	// calculating the bounds of this pointset
	// later this will be the virial radius calculation
	input->GetPoints()->ComputeBounds();
	input->GetPoints()->GetBounds(bounds);
	double xMax=bounds[1];
	double yMax=bounds[3];
	double zMax=bounds[5];
	vtkErrorMacro("bounds info is "\
								<< xMax << " " \
								<< yMax << " "\
								<< zMax);
	// setting the upper bounds appropriately
	double upperBounds[3]={xMax,yMax,zMax};
	*binLowerBound=1;
	*binUpperBound=log(ComputeRadialDistance(this->Center,upperBounds));
	*binWidth=.33;
}

//----------------------------------------------------------------------------
int vtkProfileFilter::GetBinNum(double point[],
																double binLowerBound,\
																double binUpperBound,\
																double binWidth)
{
	//TODO: implement
	return 0;
}

//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkTable* output = this->GetOutput();	
	// Allocate data structures
	double binLowerBound,binUpperBound,binWidth;
	this->GetBinInfo(input,&binLowerBound,&binUpperBound,&binWidth);
	vtkErrorMacro("bin info is "\
								<< binLowerBound << " " \
								<< binUpperBound << " "\
								<< binWidth);
	// TODO: dummy allocation
	vtkSmartPointer<vtkFloatArray> XArray=vtkSmartPointer<vtkFloatArray>::New();
		XArray->DeepCopy(input->GetPointData()->GetArray("mass"));
		XArray->SetName("mass");
	for(int nextPointId = 0;\
	 		nextPointId < input->GetPoints()->GetNumberOfPoints();\
	 		++nextPointId)
		{
			double* nextPoint=GetPoint(input,nextPointId);
			int binNum=this->GetBinNum(nextPoint,\
																 binLowerBound,binUpperBound,binWidth);
			// Finally some memory management
			delete [] nextPoint;
		}
	// Updating the output
	// TODO: dummy output
	output->AddColumn(XArray);
	return 1;
}
