/*=========================================================================

		Program:   AstroViz plugin for ParaView
		Module:    $RCSfile: vtkNSmoothFilter.h,v $

		Copyright (c) Christine Corbett Moran
		All rights reserved.
   
	This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.


=========================================================================*/
// .NAME vtkNSmoothFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkNSmoothFilter 
// Build a Kd tree, then find the N nearest neighbors 
// (as specified by NeighborNumber) 
// average over them to find smoothed variable value. For those 
// variables which need a volume to be computed
// consider the volume as sphere around point with radius of the
// outermost neighbor point.
 // Outline of this filter:
// 1. Build Kd tree, or if running in parallel run D3 filter and use its KDT
// 2. Go through each point in output
// 		o calculate N nearest neighbors
//		o calculate smoothed quantities
// 		o add to their respective arrays.
// 3. Add the arrays of smoothed variables to the output
// .SECTION See Also
// vtkKdTree, vtkPKdTree

#ifndef __vtkNSmoothFilter_h
#define __vtkNSmoothFilter_h
#include "vtkPointSetAlgorithm.h"
class vtkMultiProcessController;

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

 	// Description:
	// By defualt this filter uses the global controller,
	// but this method can be used to set another instead.
	virtual void SetController(vtkMultiProcessController*);
	vtkGetObjectMacro(Controller, vtkMultiProcessController);

//BTX
protected:
  vtkNSmoothFilter();
  ~vtkNSmoothFilter();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(vtkInformation*,
   	vtkInformationVector**,
    vtkInformationVector*);
	vtkMultiProcessController* Controller;
  int NeighborNumber;

private:
  vtkNSmoothFilter(const vtkNSmoothFilter&);  // Not implemented.
  void operator=(const vtkNSmoothFilter&);  // Not implemented.
	// Description:
	// calculates the density by taking 4/3 pi r^3 to be the volume
	// where r=dist(pointOne,pointTwo), and diving the smoothed mass
	// which is the average mass in that volume by the volume

	double CalculateDensity(double pointOne[],double pointTwo[], double smoothedMass);
	// Description:
	// returns a string representing the name of the smoothed array
	vtkstd::string GetSmoothedArrayName(vtkstd::string baseName, int dataIndex);

//ETX
};
#endif
