/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNSmoothFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkNSmoothFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkNSmoothFilter 
// Build a Kd tree, then find the N nearest neighbors (as specified by NeighborNumber) 
// average over them to find smoothed variable value. For those variables which need a volume to be computed
// consider the volume as sphere around point with radius of the outermost neighbor point.
// .SECTION See Also
// vtkKdTree

#ifndef __vtkNSmoothFilter_h
#define __vtkNSmoothFilter_h

#include "vtkPointSetAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkNSmoothFilter : public vtkPointSetAlgorithm
{
public:
  static vtkNSmoothFilter *New();
  vtkTypeRevisionMacro(vtkNSmoothFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get/Set the number of neighbors to search
  vtkSetMacro(NeighborNumber, int);
  vtkGetMacro(NeighborNumber, int);

protected:
  vtkNSmoothFilter();
  ~vtkNSmoothFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);

  int NeighborNumber;

private:
  vtkNSmoothFilter(const vtkNSmoothFilter&);  // Not implemented.
  void operator=(const vtkNSmoothFilter&);  // Not implemented.
};

#endif
