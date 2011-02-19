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

	void SetFilterArray(int);
	void SetModeSelection(int);
	void SetRestrictionArray(const char *);
	//vtkSetMacro(FilterArray,int);
	//vtkGetMacro(FilterArray,int);
	//vtkSetMacro(ModeSelection,int);
	//vtkGetMacro(ModeSelection,int);
	//vtkSetMacro(RestrictionArray,char);
	//vtkGetMacro(RestrictionArray,char*);


	void SetFilter(double,double);
/*	vtkSetMacro(Filter_0,double);
	vtkGetMacro(Filter_0,double);
	vtkSetMacro(Filter_1,double);
	vtkGetMacro(Filter_1,double);*/
	//vtkSetMacro(Filter,double);
	//vtkGetMacro(Filter,double);

	void SetRestriction(double, double);

protected:
	vtkTrackFilter2();
	~vtkTrackFilter2();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
	virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);


private:
	vtkTrackFilter2(const vtkTrackFilter2&);  // Not implemented.
	void operator=(const vtkTrackFilter2&);  // Not implemented.


	int FilterArray;
	const char * ModeSelection;
	const char * RestrictionArray;

	double Filter[2];
	double Restriction[2];



};

#endif
