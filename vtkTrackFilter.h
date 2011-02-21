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

	vtkSetMacro(Mode,int);
	vtkGetMacro(Mode,int);

/*	vtkSetVectorMacro(FilterBound, double, 2);
	vtkGetVectorMacro(FilterBound, double, 2);
	vtkSetVectorMacro(RestrictionBound, double, 2);
	vtkGetVectorMacro(RestrictionBound, double, 2);
*/
	
	vtkSetMacro(FilterHi,double);
	vtkGetMacro(FilterHi,double);
	vtkSetMacro(FilterLow,double);
	vtkGetMacro(FilterLow,double);

	vtkSetMacro(RestrictionHi,double);
	vtkGetMacro(RestrictionHi,double);
	vtkSetMacro(RestrictionLow,double);
	vtkGetMacro(RestrictionLow,double);


protected:
	vtkTrackFilter();
	~vtkTrackFilter();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);


private:
	vtkTrackFilter(const vtkTrackFilter&);  // Not implemented.
	void operator=(const vtkTrackFilter&);  // Not implemented.

	int Mode;
/*	double FilterBound[2];
	double RestrictionBound[2];
*/

	double FilterHi;
	double FilterLow;

	double RestrictionHi;
	double RestrictionLow;
};

#endif
