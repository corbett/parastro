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

class VTK_EXPORT vtkTrackFilter2 :  public vtkPolyDataAlgorithm
{
public:
	static vtkTrackFilter2* New();
	//vtkTypeMacro(vtkTrackFilter2, vtkPolyDataAlgorithm);
	vtkTypeRevisionMacro(vtkTrackFilter2,vtkPolyDataAlgorithm);

	void PrintSelf(ostream& os, vtkIndent indent);

	vtkSetMacro(Mode,int);
	vtkGetMacro(Mode,int);

	vtkSetMacro(HighPoint,double);
	vtkGetMacro(HighPoint,double);

	vtkSetMacro(LowPoint,double);
	vtkGetMacro(LowPoint,double);


protected:
	vtkTrackFilter2();
	~vtkTrackFilter2();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);


private:
	vtkTrackFilter2(const vtkTrackFilter2&);  // Not implemented.
	void operator=(const vtkTrackFilter2&);  // Not implemented.

	int Mode;
	double HighPoint;
	double LowPoint;
};

#endif
