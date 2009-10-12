/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkProfileFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkProfileFilter - shrink cells composing an arbitrary data set
// .SECTION Description
// Plots the mass function N(>M) as a scatter plot
// .SECTION See Also

#ifndef __vtkProfileFilter_h
#define __vtkProfileFilter_h

#include "vtkExtractHistogram.h"

class VTK_EXPORT vtkProfileFilter : public vtkExtractHistogram
{
public:
  static vtkProfileFilter* New();
  vtkTypeRevisionMacro(vtkProfileFilter, vtkExtractHistogram);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Get/Set the density parameter
  vtkSetMacro(Delta, double);
  vtkGetMacro(Delta, double);
  // Description:
  // Get/Set the number of bins
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);
  // Description:
  // Get/Set whether the bins are by radius
  vtkSetMacro(BinByRadius, int);
  vtkGetMacro(BinByRadius, int);

  // Description:
  // Get/Set whether the bins should be only from the center to the virial 
	// radius
	vtkSetMacro(CutOffAtVirialRadius,int);
	vtkGetMacro(CutOffAtVirialRadius,int);

  // Description:
  // Specify the point locations used to probe input. Any geometry
  // can be used. New style. Equivalent to SetInputConnection(1, algOutput).
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);
//BTX
protected:
  vtkProfileFilter();
  ~vtkProfileFilter();
  // Description:
  // This is called within ProcessRequest when a request asks the algorithm
  // to do its work. This is the method you should override to do whatever the
  // algorithm is designed to do. This happens during the fourth pass in the
  // pipeline execution process.
  virtual int RequestData(vtkInformation*, 
                          vtkInformationVector**, 
                          vtkInformationVector*);
	// Set in GUI, with defaults
	double Delta;
	double Center[3];
	int BinByRadius;
	int CutOffAtVirialRadius;
	void SetCenterFromGUI();
	
private:
  vtkProfileFilter(const vtkProfileFilter&); // Not implemented
  void operator=(const vtkProfileFilter&); // Not implemented

//ETX
};

#endif


