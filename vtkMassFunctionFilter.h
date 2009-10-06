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

#include "vtkTableAlgorithm.h"

class VTK_EXPORT vtkMassFunctionFilter : public vtkTableAlgorithm
{
public:
  static vtkMassFunctionFilter* New();
  vtkTypeRevisionMacro(vtkMassFunctionFilter, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

//BTX
protected:
  vtkMassFunctionFilter();
  ~vtkMassFunctionFilter();

  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Description:
  // This is called within ProcessRequest when a request asks the algorithm
  // to do its work. This is the method you should override to do whatever the
  // algorithm is designed to do. This happens during the fourth pass in the
  // pipeline execution process.
  virtual int RequestData(vtkInformation*, 
                          vtkInformationVector**, 
                          vtkInformationVector*);

private:
  vtkMassFunctionFilter(const vtkMassFunctionFilter&); // Not implemented
  void operator=(const vtkMassFunctionFilter&); // Not implemented
//ETX
};

#endif


