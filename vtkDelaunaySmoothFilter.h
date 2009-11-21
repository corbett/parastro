/*=========================================================================

		Program:   AstroViz plugin for ParaView
		Module:    $RCSfile: vtkDelaunaySmoothFilter.h,v $

		Copyright (c) Christine Corbett Moran
		All rights reserved.
   
	This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.


=========================================================================*/
// .NAME vtkDelaunaySmoothFilter 
// .SECTION Description
// vtkDelaunaySmoothFilter 
// Performs a delaunay tessellation, smoothing over neighbors of a delaunay
// cell or over volume of cell, as with density.
// .SECTION See Also
// vtkDelaunay3D

#ifndef __vtkDelaunaySmoothFilter_h
#define __vtkDelaunaySmoothFilter_h
#include "vtkPointSetAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkDelaunaySmoothFilter : public vtkPointSetAlgorithm
{
public:
  static vtkDelaunaySmoothFilter *New();
  vtkTypeRevisionMacro(vtkDelaunaySmoothFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

//BTX
protected:
  vtkDelaunaySmoothFilter();
  ~vtkDelaunaySmoothFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
   	vtkInformationVector**,
    vtkInformationVector*);


private:
  vtkDelaunaySmoothFilter(const vtkDelaunaySmoothFilter&);  // Not implemented.
  void operator=(const vtkDelaunaySmoothFilter&);  // Not implemented.

	// Description:
	// returns a string representing the name of the smoothed array
	vtkstd::string GetSmoothedArrayName(vtkstd::string baseName, int dataIndex);

//ETX
};
#endif
