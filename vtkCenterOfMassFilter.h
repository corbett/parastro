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
#include "vtkPointSetAlgorithm.h" // superclass
#include <vtkstd/string> // argument to function
class vtkMultiProcessController;
enum CenterOfMassMPIData 
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

	// Description:
	// Computes the center of mass of the vtkPointSet input
	// using the mass as the array with massArrayName. Functions in parallel
	// if a controller is set. Returns NULL if in parallel and process id != 0
	// (result isn't ready until process 0). So check for this.
	//BTX
	double* ComputeCenterOfMass(vtkPointSet* input, vtkstd::string massArrayName);
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
	// Description
	// helper function called by each process to increment totalMass and
	// totalWeighted mass based on the dataset of that process
	// compute COM is called with these variables as input
	// at the very last stage
	void UpdateCenterOfMassVariables(vtkPointSet* input,double& totalMass, 
		double totalWeightedMass[]);

	// Description:
	// helper function to compute center of mass of point set, must be called
	// at last stage, once update COM vars has been called on each process
	double* ComputeCenterOfMassFinal(vtkPointSet* input,double& totalMass, 
		double totalWeightedMass[]);
	// Description:
	// ComputeCOM helper function to calculate [m*x,m*y,m*z]
	double* ComputeWeightedMass(double& mass,double* point);
//ETX
};

#endif
