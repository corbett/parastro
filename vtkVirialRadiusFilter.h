/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkVirialRadiusFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVirialRadiusFilter
// Given an overdensity and a center, calculates and cuts off 
// the data set at the point where the density equals this overdensity.

#ifndef __vtkVirialRadiusFilter_h
#define __vtkVirialRadiusFilter_h
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkStringArray.h" // some class variables are vtkStringArrays

class vtkPointSet;
class vtkDataSet;
class vtkMultiProcessController;
//----------------------------------------------------------------------------
enum BinUpdateType
{
	ADD, 
	MULTIPLY,
	SET
};

enum ColumnType
{
	AVERAGE,
	TOTAL,
	CUMULATIVE
};


//----------------------------------------------------------------------------
class VTK_EXPORT vtkVirialRadiusFilter : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkVirialRadiusFilter* New();
  vtkTypeRevisionMacro(vtkVirialRadiusFilter, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Get/Set the softening parameter
  vtkSetMacro(Softening, double);
  vtkGetMacro(Softening, double);
  // Description:
  // Get/Set the density parameter
  vtkSetMacro(Delta, double);
  vtkGetMacro(Delta, double);
  // Description:
  // Get/Set the center
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);
 	// Description:
	// By defualt this filter uses the global controller,
	// but this method can be used to set another instead.
	virtual void SetController(vtkMultiProcessController*);
	vtkGetObjectMacro(Controller, vtkMultiProcessController);  
	
	// Description:
  // Specify the point locations used to probe input. Any geometry
  // can be used. New style. Equivalent to SetInputConnection(1, algOutput).
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);
	// Description:
	// overridden to only take in certain types of data
	virtual int FillInputPortInformation (int port, vtkInformation *info);

	// Override to specify different type of output
	virtual int FillOutputPortInformation(int vtkNotUsed(port), 
		vtkInformation* info);
	
//BTX
protected:
  vtkVirialRadiusFilter();
  ~vtkVirialRadiusFilter();
  // Description:
  // This is called within ProcessRequest when a request asks the algorithm
  // to do its work. This is the method you should override to do whatever the
  // algorithm is designed to do. This happens during the fourth pass in the
  // pipeline execution process.
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
	double Delta; 
  // Description:
	// Center around which to compute radial bins
	double Center[3];
	// Description:
	// Max distance from center point to the data set boundaries, or to
	// the virial radius if applicable
	double MaxR;
	// Description:
	// The MPI controller for this filter, if needed
	vtkMultiProcessController* Controller;
	// Description:
  // Calculates the center and the maximum distance from the center
	// based upon the user's input and the boundaries of the dataset.
	// Works in parallel if necessary
	void CalculateAndSetBounds(vtkPointSet* input, vtkDataSet* source);
private:
  vtkVirialRadiusFilter(const vtkVirialRadiusFilter&); // Not implemented
  void operator=(const vtkVirialRadiusFilter&); // Not implemented

//ETX
};

#endif


