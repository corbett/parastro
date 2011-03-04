/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkGraficReader.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGraficReader - Read points from a Grafic standard binary file
// .SECTION Description
// Read points from a Grafic standard binary file. Fully parallel. Has ability
// to read in additional attributes from an ascii file, and to only load in
// marked particles but both these functions are serial only.
#ifndef __vtkGraficReader_h
#define __vtkGraficReader_h

#include "vtkPolyDataAlgorithm.h" // superclass

#include "vtkSmartPointer.h"
#include "tipsylib/ftipsy.hpp" // functions take Grafic particle objects
#include <vtkstd/vector>

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkDoubleArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;
class vtkFloatArray;
class vtkIntArray;
class VTK_EXPORT vtkGraficReader : public vtkPolyDataAlgorithm
{
public:
  static vtkGraficReader* New();
  vtkTypeRevisionMacro(vtkGraficReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);
// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkGraficReader();
  ~vtkGraficReader();
	char* FileName;
	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

  int RequestData(vtkInformation*,vtkInformationVector**,
    vtkInformationVector*);

  vtkIdType                       ParticleIndex;
  vtkSmartPointer<vtkIdTypeArray> GlobalIds;
  vtkSmartPointer<vtkPoints>      Positions;
  vtkSmartPointer<vtkCellArray>   Vertices;

  vtkSmartPointer<vtkDoubleArray>   Velocity;
	
  vtkSmartPointer<vtkFloatArray>   Mass;
	vtkSmartPointer<vtkFloatArray>   EPS;
	vtkSmartPointer<vtkFloatArray>   Potential;
  vtkSmartPointer<vtkFloatArray>   RHO;
  vtkSmartPointer<vtkFloatArray>   Temperature;
  vtkSmartPointer<vtkFloatArray>   Metals;
  vtkSmartPointer<vtkFloatArray>   Tform;
	vtkSmartPointer<vtkIntArray>		 Type;
	
	
	
  //
  int           UpdatePiece;
  int           UpdateNumPieces;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

private:
  vtkGraficReader(const vtkGraficReader&);  // Not implemented.
  void operator=(const vtkGraficReader&);  // Not implemented.
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them 
	// in the output vector
	void AllocateAllGraficVariableArrays(vtkIdType numBodies,
																			vtkPolyData* output);
//ETX

};
#endif
