/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkAddAdditionalAttribute.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkAddAdditionalAttribute 
// .SECTION Description
// vtkAddAdditionalAttribute 
//Additional attribute files set are merged into existing input. 
// Works only in serial.
// .SECTION See Also
// vtkKdTree

#ifndef __vtkAddAdditionalAttribute_h
#define __vtkAddAdditionalAttribute_h
#include "vtkPointSetAlgorithm.h"

enum FileFormat 
{
  FORMAT_SKID_ASCII=0,
  FORMAT_HOP_DENSITY_BIN=1
};

class VTK_EXPORT vtkAddAdditionalAttribute : public vtkPointSetAlgorithm
{
public:
  static vtkAddAdditionalAttribute *New();
  vtkTypeRevisionMacro(vtkAddAdditionalAttribute,vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to get additional attributes
  vtkSetStringMacro(AttributeFile);
  vtkGetStringMacro(AttributeFile);
  // Description:
  // Set/Get the name of the additional attribute
  vtkSetStringMacro(AttributeName);
  vtkGetStringMacro(AttributeName);
  
  // Description:
  // Set/Get the name of the additional attribute
  vtkSetMacro(AttributeFileFormatType,int);
  vtkGetMacro(AttributeFileFormatType,int);


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
	char* AttributeFile;
	char* AttributeName;
  int AttributeFileFormatType;
private:
  vtkAddAdditionalAttribute(const vtkAddAdditionalAttribute&);  // Not implemented.
  void operator=(const vtkAddAdditionalAttribute&);  // Not implemented.
	// Description:
	// Reads this file in as an additional attribute array.
 //  If a marked file was specified, it reads in 
	// only the particles at indices which were marked
	int ReadAdditionalAttributeFile(
		vtkDataArray* globalIdArray, 
		vtkPointSet* output);
//ETX
};

#endif
