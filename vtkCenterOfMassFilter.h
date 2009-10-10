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

class VTK_GRAPHICS_EXPORT vtkCenterOfMassFilter : public vtkPointSetAlgorithm
{
public:
  static vtkCenterOfMassFilter *New();
  vtkTypeRevisionMacro(vtkCenterOfMassFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Get/Set the density parameter
  vtkSetMacro(Overdensity, double);
  vtkGetMacro(Overdensity, double);

protected:
  vtkCenterOfMassFilter();
  ~vtkCenterOfMassFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
	double Overdensity;

private:
  vtkCenterOfMassFilter(const vtkCenterOfMassFilter&);  // Not implemented.
  void operator=(const vtkCenterOfMassFilter&);  // Not implemented.
	double* CalculateWeightedMass(double& mass,double* point);
};

#endif
