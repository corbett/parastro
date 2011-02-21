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

#include "vtkSMProperty.h"
#include "vtkSMProxy.h"
#include "vtkSMProxyManager.h"
#include "vtkSMPropertyHelper.h"


vtkCxxRevisionMacro(vtkTrackFilter2, "$Revision: 0.1 $");
vtkStandardNewMacro(vtkTrackFilter2);

//----------------------------------------------------------------------------
vtkTrackFilter2::vtkTrackFilter2()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	//this->ModeSelection.assign("");

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

/*
//----------------------------------------------------------------------------
int vtkTrackFilter2::RequestInformation(
		vtkInformation* vtkNotUsed(request),
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
	return 1;
}*/

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

	vtkDataArray* filterArray = this->GetInputArrayToProcess(0, inputVector);
	vtkDataArray* restrictionArray = this->GetInputArrayToProcess(1, inputVector);

	vtkErrorMacro("FilterArray: "<<filterArray->GetName());
	vtkErrorMacro("RestrictionArray: "<<restrictionArray->GetName());


	vtkErrorMacro("Filter_0: "<<this->FilterBounds[0]);
	vtkErrorMacro("Filter_1: "<<this->FilterBounds[1]);

	vtkErrorMacro("Restriction_0: "<<this->RestrictionBounds[0]);
	vtkErrorMacro("Restriction_1: "<<this->RestrictionBounds[1]);

	/*
	vtkErrorMacro("Filter_0: "<<this->Filter[0]);
	vtkErrorMacro("Filter_1: "<<this->Filter[1]);

	vtkErrorMacro("Restriction_0: "<<this->Restriction[0]);
	vtkErrorMacro("Restriction_1: "<<this->Restriction[1]);
	*/

	//vtkErrorMacro("ModeSelection: "<<this->ModeSelection.c_str());




	//---- actualt filtering code
/*
	vtkDataArray * filterArray = input->GetPointData()->GetArray(this->FilterArray);
	vtkDataArray * restrictionArray = input->GetPointData()->GetArray(this->RestrictionArray);



	vtkSmartPointer<vtkCellArray> tracks = input->GetLines();
	int mode = 0;//debug
	bool savetrack;

	vtkSmartPointer<vtkIdList> selectedPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkCellArray> selectedTracks = vtkSmartPointer<vtkCellArray>::New();
	selectedPoints->Initialize();

	vtkIdType * pts;
	vtkIdType npts;

	tracks->InitTraversal();

	// select tracks where one point meets criteria
	if (mode==0)
	{

		for (int i = 0; i<tracks->GetNumberOfCells(); ++i)
		{
			tracks->GetNextCell(npts, pts);
			for (int j = 0; j < npts; ++j)
			{
				if (filterArray->GetTuple1(*(pts+j)) >= this->Filter[0] &&
					filterArray->GetTuple1(*(pts+j)) <= this->Filter[1])
				{
					savetrack=true;
					break;
				}
			}
		}
	}
	// select tracks where every point meets critria
	else if (mode == 1)
	{
		
	}
	//select tracks where any point within restriction meets criteria
	else if (mode == 2)
	{
		for (int i = 0; i<tracks->GetNumberOfCells(); ++i)
		{
			tracks->GetNextCell(npts, pts);
			for (int j = 0; j < npts; ++j)
			{
				if (restrictionArray->GetTuple1(*(pts+j)) >= this->Restriction_0 &&
					restrictionArray->GetTuple1(*(pts+j)) <= this->Restriction_1)
				{
					if (!(filterArray->GetTuple1(*(pts+j)) >= this->Filter[0] &&
						filterArray->GetTuple1(*(pts+j)) <= this->Filter[1]))
					{
						savetrack = false;
						break;
					}

				}
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
*/
	timer->StopTimer();
	vtkErrorMacro(" track selection took: " << timer->GetElapsedTime() << " s");
	
	return 1;
}


/*
void vtkTrackFilter2::SetFilterArray(const char * sel)
{
	//vtkErrorMacro("FilterArray sub: "<<index);
	this->FilterArray = vtkstd::string(sel);
	this->Modified();
}

void vtkTrackFilter2::SetModeSelection(const char * sel)
{
	//vtkErrorMacro("SetModeSelection sub: "<<sel);
	this->ModeSelection = vtkstd::string(sel);
	this->Modified();
}

void vtkTrackFilter2::SetRestrictionArray(const char * sel)
{
	//vtkErrorMacro("RestrictionArray sub: "<<sel);
	this->RestrictionArray = vtkstd::string(sel);
	//vtkErrorMacro("RestrictionArray sub: " << this->RestrictionArray.c_str());
	this->Modified();
}
*/
/*
void vtkTrackFilter2::SetFilter(double f0, double f1)
{
	this->Filter[0] = f0;
	this->Filter[1] = f1;
	this->Modified();
}

void vtkTrackFilter2::SetRestriction(double f0, double f1)
{
	this->Restriction[0] = f0;
	this->Restriction[1] = f1;
	this->Modified();
}
*/