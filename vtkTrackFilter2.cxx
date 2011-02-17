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

	vtkErrorMacro("sdf"<<this->FilterHi);


	timer->StopTimer();
	vtkErrorMacro(" track selection took: " << timer->GetElapsedTime() << " s");
	
	return 1;
}
