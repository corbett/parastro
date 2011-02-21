/*=========================================================================

  Program:   
  Module:    vtkTrackFilter2.h

  Copyright (c) Rafael Küng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTrackFilter2 - Filter Tracks
// .SECTION Description
// vtkTrackFilter2 allows to select Lines, where a certain PointData Value meets specified condition

#ifndef __vtkTrackFilter2_h
#define __vtkTrackFilter2_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkstd/vector>
#include "vtkStdString.h"


class VTK_EXPORT vtkTrackFilter2 :  public vtkPolyDataAlgorithm
{
public:
	static vtkTrackFilter2* New();
	//vtkTypeMacro(vtkTrackFilter2, vtkPolyDataAlgorithm);
	vtkTypeRevisionMacro(vtkTrackFilter2,vtkPolyDataAlgorithm);

	void PrintSelf(ostream& os, vtkIndent indent);

	vtkSetVectorMacro(FilterBounds, double, 2);
	vtkGetVectorMacro(FilterBounds, double, 2);

	vtkSetVectorMacro(RestrictionBounds, double, 2);
	vtkGetVectorMacro(RestrictionBounds, double, 2);


	//void SetFilter(double,double);
	//void SetRestriction(double, double);

protected:
	vtkTrackFilter2();
	~vtkTrackFilter2();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
	//virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);


private:
	vtkTrackFilter2(const vtkTrackFilter2&);  // Not implemented.
	void operator=(const vtkTrackFilter2&);  // Not implemented.


	//vtkStdString FilterArray;
	//vtkStdString ModeSelection;
	//vtkStdString RestrictionArray;

	//double Filter[2];
	//double Restriction[2];

	double FilterBounds[2];
	double RestrictionBounds[2];



};

#endif
