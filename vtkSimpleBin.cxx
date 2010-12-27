/*=========================================================================

TODO: change output so something useful..


  Program:   
  Module:    vtkSimpleBin.cxx

  Copyright (c) Rafael Küng

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSimpleBin.h"

#include "AstroVizHelpersLib/AstroVizHelpers.h"

#include "vtkPolyData.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "vtkTimerLog.h"
#include "vtkSmartPointer.h"
#include "vtkDataSetAttributes.h"

#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkPolyLine.h"
#include "vtkDataArray.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"
#include "vtkMultiBlockDataSet.h"

#include <vector>


vtkCxxRevisionMacro(vtkSimpleBin, "$Revision: 0.1 $");
vtkStandardNewMacro(vtkSimpleBin);

//----------------------------------------------------------------------------
vtkSimpleBin::vtkSimpleBin()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->SetInputArrayToProcess(
		0,
		0,
		0,
		vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
		vtkDataSetAttributes::SCALARS);

	this->DoStdDerr = false;
}

//----------------------------------------------------------------------------
vtkSimpleBin::~vtkSimpleBin()
{
}

//----------------------------------------------------------------------------
void vtkSimpleBin::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkSimpleBin::FillInputPortInformation(int port, vtkInformation* info)
{
	this->Superclass::FillInputPortInformation(port, info);
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}
//----------------------------------------------------------------------------
int vtkSimpleBin::FillOutputPortInformation(
	int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
	return 1;
}

//----------------------------------------------------------------------------
int vtkSimpleBin::RequestData(vtkInformation*,
								vtkInformationVector** inputVector,
								vtkInformationVector* outputVector)
{
	vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
	timer->StartTimer();

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkTable* output = vtkTable::GetData(outputVector);
	//vtkMultiBlockDataSet * multiOutput = vtkMultiBlockDataSet::GetData(outputVector);

	vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);
	vtkPointData * pData = input->GetPointData();
	int nArr = pData->GetNumberOfArrays();

	vtkstd::vector<int> counter; // stores, how many entries a bin has

	// calculate nBin (simple, only works for int values..
	double range [2];
	filterArray->GetRange(range);
	int nBin = range[1] - range[0] +1;



	// create table structre
	// mean values
	for (int i = 0; i<nArr; i++)
	{
		vtkAbstractArray * inArr = pData->GetArray(i);
		vtkAbstractArray * arr_mean = vtkAbstractArray::CreateArray(inArr->GetDataType());
		vtkstd::string name = inArr->GetName();
		name.append("_mean");
		arr_mean->SetName(name.c_str());
		output->AddColumn(arr_mean);
		arr_mean->Delete();
	}
	// Standart derrivation fields
	if (this->DoStdDerr)
	{
		for (int i = 0; i<nArr; i++)
		{
			vtkAbstractArray * inArr = pData->GetArray(i);
			vtkAbstractArray * arr_StdDerr = vtkDoubleArray::New();
			vtkstd::string name = inArr->GetName();
			name.append("_StdDerr");
			arr_StdDerr->SetName(name.c_str());
			output->AddColumn(arr_StdDerr);
			arr_StdDerr->Delete();
		}
	}

	// init every field to 0
	counter.resize(nBin);
	output->SetNumberOfRows(nBin);
	for (int i = 0; i<nBin;i++)
	{
		counter.at(i) = 0;
		for (int j = 0; j<nArr; j++)
		{
			output->SetValue(i,j,0.0);
		}
	}
	if(this->DoStdDerr)
	{
		for (int i = 0; i<nBin;i++)
		{
			for (int j = nArr; j<2*nArr;j++)
			{
				output->SetValue(i,j,0.0);
			}
		}
	}



	// sum up the values
	for (int i = 0; i<input->GetNumberOfPoints(); i++)
	{
		int bin = range[0] + *filterArray->GetTuple(i);
		counter.at(bin)++;

		for (int j = 0; j<nArr; j++)
		{
			double oldval = output->GetValue(bin,j).ToDouble();
			double orgval = *pData->GetArray(j)->GetTuple(i);
			output->SetValue(bin,j, oldval + orgval);
		}
	}

	output->Update();

	// get the average
	for (int bin = 0; bin<nBin; bin++)
	{
		for (int j = 0; j<nArr; j++)
		{
			double oldval = output->GetValue(bin,j).ToDouble();
			output->SetValue(bin,j, oldval / (double)counter.at(bin));
		}
	}

	output->Update();

	// calculate stdderr
	if(this->DoStdDerr)
	{
		for (int i = 0; i<input->GetNumberOfPoints(); i++)
		{
			int bin = range[0] + *filterArray->GetTuple(i);
			for (int j = nArr; j<2*nArr; j++)
			{
				double oldval = output->GetValue(bin,j).ToDouble();
				double newval = *pData->GetArray(j-nArr)->GetTuple(i);
				double mean = output->GetValue(bin,j-nArr).ToDouble();
				output->SetValue(bin,j, oldval + (newval - mean)*(newval - mean));
			}
		}

		output->Update();

		for (int bin = 0; bin<nBin; bin++)
		{
			for (int j = nArr; j<2*nArr; j++)
			{
				double oldval = output->GetValue(bin,j).ToDouble();
				output->SetValue(bin,j, sqrt(oldval / (double)(counter.at(bin)-1)));
			}
		}

		output->Update();
	}


	timer->StopTimer();
	vtkErrorMacro(" binning took: " << timer->GetElapsedTime() << " s");

	
	return 1;
}
