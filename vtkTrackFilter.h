/*=========================================================================

  Program:   
  Module:    vtkTrackFilter.h

  Copyright (c) Rafael Küng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTrackFilter - Filter Tracks
// .SECTION Description
// vtkTrackFilter allows to select Lines, where a certain PointData Value meets specified condition

#ifndef __vtkTrackFilter_h
#define __vtkTrackFilter_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkstd/vector>

class VTK_EXPORT vtkTrackFilter :  public vtkPolyDataAlgorithm
{
public:
	static vtkTrackFilter* New();
	//vtkTypeMacro(vtkTrackFilter, vtkPolyDataAlgorithm);
	vtkTypeRevisionMacro(vtkTrackFilter,vtkPolyDataAlgorithm);

	void PrintSelf(ostream& os, vtkIndent indent);

	vtkSetMacro(HighPoint,double);
	vtkGetMacro(HighPoint,double);

	vtkSetMacro(LowPoint,double);
	vtkGetMacro(LowPoint,double);


protected:
	vtkTrackFilter();
	~vtkTrackFilter();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);


private:
	vtkTrackFilter(const vtkTrackFilter&);  // Not implemented.
	void operator=(const vtkTrackFilter&);  // Not implemented.

	double HighPoint;
	double LowPoint;
	vtkPolyData* input;
	vtkDataArray* filterArray;

};

#endif
