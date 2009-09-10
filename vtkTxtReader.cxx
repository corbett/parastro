/*=========================================================================
Modified from vtkSimplePointsReader -christine

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSimplePointsReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTxtReader.h"

#include "vtkCellArray.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

vtkCxxRevisionMacro(vtkTxtReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkTxtReader);

//----------------------------------------------------------------------------
vtkTxtReader::vtkTxtReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkTxtReader::~vtkTxtReader()
{
  this->SetFileName(0);
}

//----------------------------------------------------------------------------
void vtkTxtReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";

}

//----------------------------------------------------------------------------
int vtkTxtReader::RequestData(vtkInformation*,
                                       vtkInformationVector**,
                                       vtkInformationVector* outputVector)
{
  // Make sure we have a file to read.
  if(!this->FileName)
    {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }

  // Open the input file.
  ifstream fin(this->FileName);
  if(!fin)
    {
    vtkErrorMacro("Error opening file " << this->FileName);
    return 0;
    }

  // Allocate objects to hold points and vertex cells.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  double x[3];
  while(fin >> x[0] >> x[1] >> x[2])
    {
    vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
    }
  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts);

  return 1;
}
