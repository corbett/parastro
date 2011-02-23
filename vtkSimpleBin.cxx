/*=========================================================================

TODO:
	- this->bincount isnt updated when switching set #bins manually
	- implement better binning (workinprogress)


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
#include "vtkPointSet.h"

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
#include "vtkMath.h"

#include <vector>
#include <math.h>


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

	//this->DoStdDerr = false;
	this->IntBin = false;
	this->LogScale = false;
	this->BinCount = 10;
	this->Del0Row = true;
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
	info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
	info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
	info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkStructuredGrid");
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

	//vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);
	vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);
	vtkstd::string filterArrayName = filterArray->GetName();
	vtkPointData * pData = input->GetPointData();
	int nArr = pData->GetNumberOfArrays();
	
	// find id of filter array
	int filterArrayId=-1;
	for (int i = 0; i<nArr; ++i)
	{
		vtkstd::string name = pData->GetArray(i)->GetName();
		if (name.compare(filterArrayName)==0)
		{
			filterArrayId = i;
			break;
		}
	}
	if (filterArrayId==-1){vtkErrorMacro("Filterarray not found");return 0;}
	
	output->Initialize();

	vtkAbstractArray * arr;
	vtkDataArray * countArr;
	vtkDataArray * valueArr;

	// init counter column
	arr = vtkAbstractArray::CreateArray(VTK_UNSIGNED_INT);
	arr->SetName("NumberOfEntriesInBin");
	output->AddColumn(arr);
	//countArr->Delete();
	countArr = vtkDataArray::SafeDownCast(output->GetColumn(0));

	// init bin column
	arr = vtkAbstractArray::CreateArray(filterArray->GetDataType());
	arr->SetName(filterArrayName.c_str());
	output->AddColumn(arr);
	//valueArr->Delete();
	valueArr = vtkDataArray::SafeDownCast(output->GetColumn(1));

	// init the other columns (mean values)
	for (int i = 0; i<nArr; i++)
	{
		vtkAbstractArray * inArr = pData->GetArray(i);
		vtkstd::string arrayName = inArr->GetName();
		arr = vtkAbstractArray::CreateArray(inArr->GetDataType());
		arrayName.append("_mean");
		arr->SetName(arrayName.c_str());
		output->AddColumn(arr);
		//arr->Delete();
	}

	// get range
	double rangeV[2];
	filterArray->GetRange(rangeV);
	//define shortcuts
	double range = rangeV[1]-rangeV[0];
	double min = rangeV[0];
	double max = rangeV[1];
	if(this->LogScale)
	{
		if (range == 0 || min == 0 || max == 0)
			{vtkErrorMacro("Error, log10(0) in range, min, max");return 0;}
		min = log10(min);
		max = log10(max);
		range = max - min;
	}


	// calculate this->BinCount
	if(this->IntBin)
	//(simple, only works for int values..
	{
		++range;
		this->BinCount = range;
	}
	// else: bincount is set in gui
	int nBin = this->BinCount; //define shortcut


	// init every field to 0
	for (int i = 0; i<nBin;i++){output->InsertNextBlankRow(0.0);}
	output->Update();


	// set up binned value
	double val;
	for (int i = 0; i<nBin;i++)
	{
		if(this->IntBin)
		{
			val = min+range/(double)nBin*i;
		}
		else
		{
			val = min+range/(double)nBin*((double)i+0.5);
		}
		if(this->LogScale){val = pow(10,val);}
		valueArr->SetTuple1(i,val);
	}
	output->Update();


	// sum up the values
	int binnr;
	for (int i = 0; i<input->GetNumberOfPoints(); ++i)
	{
		//check for arrays with mass = 0, dont count them
		if(pData->GetArray("Mvir")->GetTuple1(i) == 0){continue;}
		
		double value = filterArray->GetTuple1(i);
		if(this->LogScale) {value = log10(value);}

		if(this->IntBin)
		{
			double binnrd = ( value - min )/range*nBin;
			binnr = binnrd;
		}
		else
		{
			double binnrd = ( value - min )/range*(nBin-1);
			binnr = vtkMath::Round(binnrd);
		}

		countArr->SetTuple1(binnr, countArr->GetTuple1(binnr)+1);

		//sum up the other arrays
		// reminder: +2 offset because of 2 additional arrays in output
		for (int j = 0; j<nArr; ++j)
		{
			double oldval = output->GetValue(binnr,j+2).ToDouble(); 
			double orgval = *pData->GetArray(j)->GetTuple(i);
			output->SetValue(binnr,j+2, oldval + orgval);
		}
	}

	// get the average
	vtkstd::vector<int> remrow; //saves the rows with later are deleted
	remrow.clear();

	for (int row = 0; row<nBin; ++row)
	{
		int num = output->GetValue(row,0).ToInt(); //get number of values in this bin
		if (num == 0)
		{
			if (this->Del0Row){remrow.push_back(row);}
			else
			{
				for (int j = 0; j<nArr; ++j)
				{
					output->SetValue(row,j+2, 0); // +2 offset because +2 arrays
				}
			}
		}
		else
		{
			for (int j = 0; j<nArr; ++j)
			{
				double val = output->GetValue(row,j+2).ToDouble() / (double)num;
				output->SetValue(row,j+2, val); // +2 offset because +2 arrays
			}
		}
	}

	if (this->Del0Row)
	{
		for (vtkstd::vector<int>::reverse_iterator rit = remrow.rbegin();
			rit!=remrow.rend();
			++rit)
		{
			output->RemoveRow(*rit);
		}
	}

	output->Update();

	timer->StopTimer();

	vtkstd::stringstream ss;
	ss<<"\n\nSimpleBin run successfully!\n\n";
	ss<<"   No of Input Points: ";
	ss<<input->GetNumberOfPoints()<<"\n";
	ss<<"   No of Input Arrays: ";
	ss<<nArr<<"\n\n";
	ss<<"   Integer Bins?     : ";
	ss<<this->IntBin<<"\n";
	ss<<"   Log Scale?        : ";
	ss<<this->LogScale<<"\n";
	ss<<"   delete NaN Bins?  : ";
	ss<<this->Del0Row<<"\n\n";
	ss<<"   No of Bins        : ";
	ss<<nBin<<"\n";
	ss<<"   deleted NaN Bins  : ";
	ss<<remrow.size()<<"\n\n";
	ss<<"   Remaining Bins    : ";
	ss<<nBin-remrow.size()<<"\n";
	ss<<"   Time taken        : ";
	ss<<timer->GetElapsedTime()<<" s\n";
	vtkErrorMacro(<<ss.str());

	return 1;
}