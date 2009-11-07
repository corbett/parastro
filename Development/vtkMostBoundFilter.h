/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkMostBoundFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMostBoundFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkMostBoundFilter 
// finds the center of mass of a collection of particles. Either of all marked
// particles or of all particles
// .SECTION See Also
// vtkKdTree

#ifndef __vtkMostBoundFilter_h
#define __vtkMostBoundFilter_h

#include "vtkPointSetAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkMostBoundFilter : public vtkPointSetAlgorithm
{
public:
  static vtkMostBoundFilter *New();
  vtkTypeRevisionMacro(vtkMostBoundFilter,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Get/Set the density parameter
  vtkSetMacro(Overdensity, double);
  vtkGetMacro(Overdensity, double);

protected:
  vtkMostBoundFilter();
  ~vtkMostBoundFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
	// Description:
	// Set in GUI, with defaults
	// Describes the softening of the simulation which can influence the 
	// root finding
	double Softening;
	// Description:
	// Set in GUI, with defaults
	// Overdensity
	double Overdensity;

private:
  vtkMostBoundFilter(const vtkMostBoundFilter&);  // Not implemented.
  void operator=(const vtkMostBoundFilter&);  // Not implemented.
	double* GetMostBoundParticle(vtkPointSet* input);
};

#endif
