/*=========================================================================

  Program:   
  Module:    vtkTrackFilter2.cxx

  Copyright (c) Rafael Küng

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTrackFilter2.h"

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


vtkCxxRevisionMacro(vtkTrackFilter2, "$Revision: 0.1 $");
vtkStandardNewMacro(vtkTrackFilter2);

//----------------------------------------------------------------------------
vtkTrackFilter2::vtkTrackFilter2()
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
vtkTrackFilter2::~vtkTrackFilter2()
{
}

//----------------------------------------------------------------------------
void vtkTrackFilter2::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkTrackFilter2::FillInputPortInformation(int port, vtkInformation* info)
{
	this->Superclass::FillInputPortInformation(port, info);
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}
//----------------------------------------------------------------------------
int vtkTrackFilter2::FillOutputPortInformation(
	int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	return 1;
}

//----------------------------------------------------------------------------
int vtkTrackFilter2::RequestData(vtkInformation*,
								vtkInformationVector** inputVector,
								vtkInformationVector* outputVector)
{
	vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
	timer->StartTimer();

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector);
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);

	vtkSmartPointer<vtkCellArray> tracks = input->GetLines();
	vtkErrorMacro(" anz lines: " << tracks->GetNumberOfCells());

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);

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

		if (this->Mode==0)
		// one point on a track has to fullfill the condition
		{
			savetrack = false;
			for (int j =0;j<npts;j++)
			{
				if (filterArray->GetTuple1(*(pts+j)) <= this->HighPoint &&
					filterArray->GetTuple1(*(pts+j)) >= this->LowPoint)
				{
					savetrack = true;
					break;
				}
				//vtkErrorMacro(" id: " << *(pts+j));
			}
		}
		else if (this->Mode==1)
		// every point on a track has to fullfil the condition
		{
			savetrack = true;
			for (int j =0;j<npts;j++)
			{
				if (!(filterArray->GetTuple1(*(pts+j)) <= this->HighPoint &&
					filterArray->GetTuple1(*(pts+j)) >= this->LowPoint))
				{
					savetrack = false;
					break;
				}
				//vtkErrorMacro(" id: " << *(pts+j));
			}
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

	vtkPolyData* newDataSet = CopyPointsAndData(input,selectedPoints);
	output->Initialize();
	output->DeepCopy(newDataSet);
	output->SetLines(newTracks);
	newDataSet->Delete();

	timer->StopTimer();
	vtkErrorMacro(" track selection took: " << timer->GetElapsedTime() << " s");
	
	return 1;
}
