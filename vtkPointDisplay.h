/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkPointDisplay.h,v $

  Copyright (c)
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPointDisplay - Dummy Reader File, Shows point at origin
// .SECTION Description
// Dummy reader file, does nothing exept displaing a single point at 0,0,0
// Shows basic setup of a reader file


#ifndef __vtkPointDisplay_h
#define __vtkPointDisplay_h

#include "vtkPolyDataAlgorithm.h" // superclass

//#include "vtkSmartPointer.h"
//#include "tipsylib/ftipsy.hpp" // functions take tipsy particle objects
//#include <vtkstd/vector>
//#include "vtkFloatArray.h"

class VTK_EXPORT vtkPointDisplay : public vtkPolyDataAlgorithm
{
public:
	static vtkPointDisplay* New();
	vtkTypeRevisionMacro(vtkPointDisplay,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

protected:
	vtkPointDisplay();
	~vtkPointDisplay();
	char* FileName;

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

private:
  vtkPointDisplay(const vtkPointDisplay&);  // Not implemented.
  void operator=(const vtkPointDisplay&);  // Not implemented.
};

#endif