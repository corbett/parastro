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
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkPolyLine.h"
#include "vtkDataArray.h"

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
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
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
	vtkPolyData* output = vtkPolyData::GetData(outputVector);
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);
	vtkPointData * pData = input->GetPointData();

	double range [2];
	filterArray->GetRange(range);
	int nBin = range[1] - range[0] +1;

	for (int i = 0; i<pData->GetNumberOfArrays(); i++)
	{
		// erzeuge hier output arrays, mit anz tupel = anz bin
		// fuelle alle mit 0
		// evtl nach jeweils dazugehoeriges array mit standartabweichung
		
		vtkSmartPointer<vtkFloatArray> parray = vtkSmartPointer<vtkFloatArray>::New();
		parray->DeepCopy(pData->GetArray(i));
		parray->Initialize();
		parray->SetNumberOfTuples(nBin);
		parray->FillComponent(0,0);

		output->GetPointData()->AddArray(parray);
	}

	vtkSmartPointer<vtkPoints> pos = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i<nBin;i++)
	{
		pos->InsertNextPoint(i,i,i);
	}
	output->SetPoints(pos);

	std::vector<int> counter;

	int bin;

	for (int i = 0; i<input->GetNumberOfPoints(); i++)
	{
		bin = range[0] + filterArray->GetTuple1(i);
		//data.at(bin).num++;
		// addiere werte direkt in ouput
		for (int j = 0; i<pData->GetNumberOfArrays(); i++)
		{
			output->GetPointData()->GetArray(j)->SetTuple1(bin, 
				output->GetPointData()->GetArray(j)->GetTuple1(bin) +
				input->GetPointData()->GetArray(j)->GetTuple1(i));
		}

	}

	/*
	for (all bins)
	{
		for (all arrays)
		{
			ave_value = value / counter
			stdderr = ...
		}
	}
	*/

	timer->StopTimer();
	vtkErrorMacro(" binning took: " << timer->GetElapsedTime() << " s");
	
	return 1;
}
