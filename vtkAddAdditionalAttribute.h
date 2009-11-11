/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkAddAdditionalAttribute.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkAddAdditionalAttribute - shrink cells composing an arbitrary data set
// .SECTION Description
// vtkAddAdditionalAttribute 
// finds the moment of inertia tensor of a collection of particles, then
// displays graphically the principle moments of inertia
// .SECTION See Also
// vtkKdTree

#ifndef __vtkAddAdditionalAttribute_h
#define __vtkAddAdditionalAttribute_h
#include "vtkPointSetAlgorithm.h"
#include <vtkstd/vector>

class vtkMultiProcessController;
class VTK_GRAPHICS_EXPORT vtkAddAdditionalAttribute : public vtkPointSetAlgorithm
{
public:
  static vtkAddAdditionalAttribute *New();
  vtkTypeRevisionMacro(vtkAddAdditionalAttribute,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // By defualt this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);
  // Description:
  // Set/Get the name of the file from which to get additional attributes
  vtkSetStringMacro(AttributeFile);
  vtkGetStringMacro(AttributeFile);
  // Description:
  // Set/Get the name of the additional attribute
  vtkSetStringMacro(AttributeName);
  vtkGetStringMacro(AttributeName);

//BTX
protected:
  vtkAddAdditionalAttribute();
  ~vtkAddAdditionalAttribute();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);
  vtkMultiProcessController *Controller;
	char* AttributeFile;
	char* AttributeName;

private:
  vtkAddAdditionalAttribute(const vtkAddAdditionalAttribute&);  // Not implemented.
  void operator=(const vtkAddAdditionalAttribute&);  // Not implemented.
	// Description:
	// Reads this file in as an additional attribute array.
 //  If a marked file was specified, it reads in 
	// only the particles at indices which were marked
	int ReadAdditionalAttributeFile(vtkstd::vector<int>& markedParticleIndices,
		vtkPointSet* output);
//ETX
};

#endif
