/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkVirialRadiusFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVirialRadiusFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// Plots the mass function N(>M) as a scatter plot
// .SECTION See Also

#ifndef __vtkVirialRadiusFilter_h
#define __vtkVirialRadiusFilter_h
#include "vtkPointSetAlgorithm.h"
#include "vtkStringArray.h" // some class variables are vtkStringArrays

class vtkPointSet;
class vtkDataSet;
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
class VTK_EXPORT vtkVirialRadiusFilter : public vtkPointSetAlgorithm
{
public:
  static vtkVirialRadiusFilter* New();
  vtkTypeRevisionMacro(vtkVirialRadiusFilter, vtkPointSetAlgorithm);
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
  // Specify the point locations used to probe input. Any geometry
  // can be used. New style. Equivalent to SetInputConnection(1, algOutput).
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);
	// Description:
	// overridden to only take in certain types of data
	virtual int FillInputPortInformation (int port, vtkInformation *info);
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
  // Calculates the center and the maximum distance from the center
	// based upon the user's input and the boundaries of the dataset.
	void CalculateAndSetBounds(vtkPointSet* input, vtkDataSet* source);

private:
  vtkVirialRadiusFilter(const vtkVirialRadiusFilter&); // Not implemented
  void operator=(const vtkVirialRadiusFilter&); // Not implemented
//ETX
};

#endif


