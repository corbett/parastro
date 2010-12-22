/*=========================================================================

  Program:   
  Module:    vtkTrackFilter.cxx

  Copyright (c) Rafael Küng

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTrackFilter.h"

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

#include <vector>


vtkCxxRevisionMacro(vtkTrackFilter, "$Revision: 0.1 $");
vtkStandardNewMacro(vtkTrackFilter);

//----------------------------------------------------------------------------
vtkTrackFilter::vtkTrackFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->LowPoint = 0;
	this->HighPoint = 0;

	this->SetInputArrayToProcess(
		0,
		0,
		0,
		vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
		vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
vtkTrackFilter::~vtkTrackFilter()
{
}

//----------------------------------------------------------------------------
void vtkTrackFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkTrackFilter::FillInputPortInformation(int port, vtkInformation* info)
{
	this->Superclass::FillInputPortInformation(port, info);
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}
//----------------------------------------------------------------------------
int vtkTrackFilter::FillOutputPortInformation(
	int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	return 1;
}

//----------------------------------------------------------------------------
int vtkTrackFilter::RequestData(vtkInformation*,
								vtkInformationVector** inputVector,
								vtkInformationVector* outputVector)
{
	vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
	timer->StartTimer();

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector);
	this->input = vtkPolyData::GetData(inputVector[0]);

	vtkSmartPointer<vtkCellArray> tracks = this->input->GetLines();
	vtkErrorMacro(" anz lines: " << tracks->GetNumberOfCells());

	this->filterArray = this->GetInputArrayToProcess(0, inputVector);

	vtkIdType * pts;
	vtkIdType  npts;
	bool savetrack;
	vtkIdList * selectedPoints = vtkIdList::New();
	selectedPoints->Initialize();
	vtkSmartPointer<vtkCellArray> newTracks = vtkSmartPointer<vtkCellArray>::New();

	tracks->InitTraversal();
	for (int i = 0; i<tracks->GetNumberOfCells(); i++)
	{
		tracks->GetNextCell(npts,pts);
		savetrack = false;
		for (int j =0;j<npts;j++)
		{
			if (this->filterArray->GetTuple1(*(pts+j)) <= this->HighPoint &&
				this->filterArray->GetTuple1(*(pts+j)) >= this->LowPoint)
			{
				savetrack = true;
				break;
			}
			//vtkErrorMacro(" id: " << *(pts+j));
		}
		if (savetrack)
		{
			//vtkErrorMacro(" this track is seleced! tid: " << i);

			vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
		
			for (int j =0;j<npts;j++)
			{
				vtkIdType newId = selectedPoints->InsertNextId(*(pts+j));
				nextLine->GetPointIds()->InsertNextId(newId);
			}
			newTracks->InsertNextCell(nextLine);
		}
	}

	vtkPolyData* newDataSet = CopyPointsAndData(this->input,selectedPoints);
	output->Initialize();
	output->DeepCopy(newDataSet);
	output->SetLines(newTracks);
	newDataSet->Delete();

	timer->StopTimer();
	vtkErrorMacro(" track selection took: " << timer->GetElapsedTime() << " s");
	
	return 1;
}
