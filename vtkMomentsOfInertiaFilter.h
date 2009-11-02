/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkMomentsOfInertiaFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMomentsOfInertiaFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkMomentsOfInertiaFilter 
// finds the moment of inertia tensor of a collection of particles, then
// displays graphically the principle moments of inertia
// .SECTION See Also
// vtkKdTree

#ifndef __vtkMomentsOfInertiaFilter_h
#define __vtkMomentsOfInertiaFilter_h

#include "vtkPointSetAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkMomentsOfInertiaFilter : public vtkPointSetAlgorithm
{
public:
  static vtkMomentsOfInertiaFilter *New();
  vtkTypeRevisionMacro(vtkMomentsOfInertiaFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
protected:
  vtkMomentsOfInertiaFilter();
  ~vtkMomentsOfInertiaFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
private:
  vtkMomentsOfInertiaFilter(const vtkMomentsOfInertiaFilter&);  // Not implemented.
  void operator=(const vtkMomentsOfInertiaFilter&);  // Not implemented.
};

#endif
