/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkMassFunctionFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMassFunctionFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// Plots the mass function N(>M) as a scatter plot
// .SECTION See Also

#ifndef __vtkMassFunctionFilter_h
#define __vtkMassFunctionFilter_h

#include "vtkRectilinearGridAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkMassFunctionFilter : public vtkRectilinearGridAlgorithm
{
public:
  static vtkMassFunctionFilter *New();
  vtkTypeRevisionMacro(vtkMassFunctionFilter,vtkRectilinearGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkMassFunctionFilter();
  ~vtkMassFunctionFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);

private:
  vtkMassFunctionFilter(const vtkMassFunctionFilter&);  // Not implemented.
  void operator=(const vtkMassFunctionFilter&);  // Not implemented.
};

#endif
