/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkCenterOfMassFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCenterOfMassFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkCenterOfMassFilter 
// finds the center of mass of a collection of particles. Either of all marked
// particles or of all particles
// .SECTION See Also
// vtkKdTree

#ifndef __vtkCenterOfMassFilter_h
#define __vtkCenterOfMassFilter_h
#include "vtkPointSetAlgorithm.h"
#include "vtkMultiProcessController.h"

enum Data 
{
	TOTAL_MASS,
	TOTAL_WEIGHTED_MASS
};

class VTK_GRAPHICS_EXPORT vtkCenterOfMassFilter : public vtkPointSetAlgorithm
{
public:
  static vtkCenterOfMassFilter *New();
  vtkTypeRevisionMacro(vtkCenterOfMassFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // By defualt this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

protected:
  vtkCenterOfMassFilter();
  ~vtkCenterOfMassFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  vtkMultiProcessController *Controller;
private:
  vtkCenterOfMassFilter(const vtkCenterOfMassFilter&);  // Not implemented.
  void operator=(const vtkCenterOfMassFilter&);  // Not implemented.
};

#endif
