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
  // Get/Set the number of neighbors to search
  vtkSetMacro(NeighborNumber, int);
  vtkGetMacro(NeighborNumber, int);

protected:
  vtkCenterOfMassFilter();
  ~vtkCenterOfMassFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  int NeighborNumber;

private:
  vtkCenterOfMassFilter(const vtkCenterOfMassFilter&);  // Not implemented.
  void operator=(const vtkCenterOfMassFilter&);  // Not implemented.
	// Description:
	// calculates the density by taking 4/3 pi r^3 to be the volume
	// where r=dist(pointOne,pointTwo), and diving the smoothed mass
	// which is the average mass in that volume by the volume
	float CalculateDensity(double pointOne[3],double pointTwo[3], float& smoothedMass);
};

#endif
