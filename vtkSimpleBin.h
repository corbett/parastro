/*=========================================================================

  Program:   
  Module:    vtkSimpleBin.h

  Copyright (c) Rafael Küng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSimpleBin - Filter Tracks
// .SECTION Description
// vtkSimpleBin allows to select Lines, where a certain PointData Value meets specified condition

#ifndef __vtkSimpleBin_h
#define __vtkSimpleBin_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkstd/string>


class VTK_EXPORT vtkSimpleBin :  public vtkPolyDataAlgorithm
{
public:
	static vtkSimpleBin* New();
	//vtkTypeMacro(vtkSimpleBin, vtkPolyDataAlgorithm);
	vtkTypeRevisionMacro(vtkSimpleBin,vtkPolyDataAlgorithm);

	void PrintSelf(ostream& os, vtkIndent indent);

	//vtkSetMacro(DoStdDerr,bool);
	//vtkGetMacro(DoStdDerr,bool);
	vtkSetMacro(IntBin,bool);
	vtkGetMacro(IntBin,bool);
	vtkSetMacro(LogScale,bool);
	vtkGetMacro(LogScale,bool);
	vtkSetMacro(Del0Row,bool);
	vtkGetMacro(Del0Row,bool);
	vtkSetMacro(BinCount,int);
	vtkGetMacro(BinCount,int);


protected:
	vtkSimpleBin();
	~vtkSimpleBin();

	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);
	virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
	vtkSimpleBin(const vtkSimpleBin&);  // Not implemented.
	void operator=(const vtkSimpleBin&);  // Not implemented.

	//bool DoStdDerr;
	bool IntBin;
	bool LogScale;
	int BinCount;
	bool Del0Row;

};

#endif



