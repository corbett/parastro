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

	this->FilterHi = 0;
	this->FilterLow = 0;
	this->RestrictionHi = 0;
	this->RestrictionLow = 0;

	this->SetInputArrayToProcess(
		0,
		0,
		0,
		vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
		vtkDataSetAttributes::SCALARS);
	
	this->SetInputArrayToProcess(
		1,
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
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);

	vtkSmartPointer<vtkCellArray> tracks = input->GetLines();
	vtkErrorMacro(" anz lines: " << tracks->GetNumberOfCells());

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);
	vtkDataArray* restrictionArray = this->GetInputArrayToProcess(1, inputVector);

	vtkErrorMacro(" array1: "<<filterArray->GetName());
	vtkErrorMacro(" array2: "<<restrictionArray->GetName());

	vtkErrorMacro(" mode  : "<<this->Mode);

	vtkErrorMacro(" FiltB0: "<<this->FilterLow);
	vtkErrorMacro(" FiltB1: "<<this->FilterHi);
	vtkErrorMacro(" RestB0: "<<this->RestrictionLow);
	vtkErrorMacro(" RestB1: "<<this->RestrictionHi);







	bool savetrack;

	vtkSmartPointer<vtkIdList> selectedPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkCellArray> selectedTracks = vtkSmartPointer<vtkCellArray>::New();
	selectedPoints->Initialize();

	vtkIdType * pts;
	vtkIdType npts;

	tracks->InitTraversal();

	// select tracks where one point meets criteria
	if (this->Mode==0)
	{
		for (int i = 0; i<tracks->GetNumberOfCells(); ++i)
		{
			savetrack = false;
			tracks->GetNextCell(npts, pts);
			for (int j = 0; j < npts; ++j)
			{
				if (filterArray->GetTuple1(*(pts+j)) >= this->FilterLow &&
					filterArray->GetTuple1(*(pts+j)) <= this->FilterHi)
				{
					savetrack=true;
					break;
				}
			}
			
			if (savetrack)
			{
				vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
				for (int j =0;j<npts;j++)
				{
					vtkIdType newId = selectedPoints->InsertNextId(*(pts+j));
					nextLine->GetPointIds()->InsertNextId(newId);
				}
				selectedTracks->InsertNextCell(nextLine);
			}

		}
	}
	// select tracks where every point meets critria
	else if (this->Mode == 1)
	{
		for (int i = 0; i<tracks->GetNumberOfCells(); ++i)
		{
			savetrack = true;
			tracks->GetNextCell(npts, pts);
			for (int j = 0; j < npts; ++j)
			{
				if (!(filterArray->GetTuple1(*(pts+j)) >= this->FilterLow &&
					filterArray->GetTuple1(*(pts+j)) <= this->FilterHi))
				{
					savetrack=false;
					break;
				}
			}
			
			if (savetrack)
			{
				vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
				for (int j =0;j<npts;j++)
				{
					vtkIdType newId = selectedPoints->InsertNextId(*(pts+j));
					nextLine->GetPointIds()->InsertNextId(newId);
				}
				selectedTracks->InsertNextCell(nextLine);
			}

		}
	}

	//select tracks where any point within restriction meets criteria
	else if (this->Mode == 2)
	{
		for (int i = 0; i<tracks->GetNumberOfCells(); ++i)
		{
			tracks->GetNextCell(npts, pts);
			for (int j = 0; j < npts; ++j)
			{
				savetrack = false;
				if (restrictionArray->GetTuple1(*(pts+j)) >= this->RestrictionLow &&
					restrictionArray->GetTuple1(*(pts+j)) <= this->RestrictionHi)
				{
					if (filterArray->GetTuple1(*(pts+j)) >= this->FilterLow &&
						filterArray->GetTuple1(*(pts+j)) <= this->FilterHi)
					{
						savetrack = true;
						break;
					}

				}
			}

			if (savetrack)
			{
				vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
				for (int j =0;j<npts;j++)
				{
					vtkIdType newId = selectedPoints->InsertNextId(*(pts+j));
					nextLine->GetPointIds()->InsertNextId(newId);
				}
				selectedTracks->InsertNextCell(nextLine);
			}

		}		
	}
	// something strange went wrong...
	else
	{
		vtkErrorMacro("some nasty error occured...");
		return 0;
	}


	output->DeepCopy(CopyPointsAndData(input,selectedPoints));
	output->SetLines(selectedTracks);

/*	vtkPolyData* newDataSet = CopyPointsAndData(input,selectedPoints);
	output->Initialize();
	output->DeepCopy(newDataSet);
	output->SetLines(newTracks);
	newDataSet->Delete();*/

	timer->StopTimer();
	vtkErrorMacro(" track selection took: " << timer->GetElapsedTime() << " s");





	/*
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
	*/
	return 1;
}
