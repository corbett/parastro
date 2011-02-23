/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkRamsesReader.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRamsesReader - Read points from a Ramses standard binary file
// .SECTION Description
// Read points from a Ramses standard binary file. Fully parallel. Has ability
// to read in additional attributes from an ascii file, and to only load in
// marked particles but both these functions are serial only.
#ifndef __vtkRamsesReader_h
#define __vtkRamsesReader_h

#include "vtkPolyDataAlgorithm.h" // superclass

#include "vtkSmartPointer.h"
#include "tipsylib/ftipsy.hpp" // functions take Ramses particle objects
#include <vtkstd/vector>

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkDoubleArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

class VTK_EXPORT vtkRamsesReader : public vtkPolyDataAlgorithm
{
public:
  static vtkRamsesReader* New();
  vtkTypeRevisionMacro(vtkRamsesReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

  // Description:
  // An H5Part file may contain multiple arrays
  // a GUI (eg Paraview) can provide a mechanism for selecting which data arrays
  // are to be read from the file. The PointArray variables and members can
  // be used to query the names and number of arrays available
  // and set the status (on/off) for each array, thereby controlling which
  // should be read from the file. Paraview queries these point arrays after
  // the (update) information part of the pipeline has been updated, and before the
  // (update) data part is updated.
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAll();
  void        EnableAll();
  void        Disable(const char* name);
  void        Enable(const char* name);
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkRamsesReader();
  ~vtkRamsesReader();
	char* FileName;
	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

  int RequestData(vtkInformation*,vtkInformationVector**,
    vtkInformationVector*);

  vtkIdType                       ParticleIndex;
  vtkSmartPointer<vtkIdTypeArray> GlobalIds;
  vtkSmartPointer<vtkPoints>      Positions;
  vtkSmartPointer<vtkCellArray>   Vertices;

  vtkSmartPointer<vtkDoubleArray>   Potential;
  vtkSmartPointer<vtkDoubleArray>   Mass;
  vtkSmartPointer<vtkDoubleArray>   EPS;
  vtkSmartPointer<vtkDoubleArray>   RHO;
  vtkSmartPointer<vtkDoubleArray>   Hsmooth;
  vtkSmartPointer<vtkDoubleArray>   Temperature;
  vtkSmartPointer<vtkDoubleArray>   Metals;
  vtkSmartPointer<vtkDoubleArray>   Tform;
	vtkSmartPointer<vtkDoubleArray>		 Type;
	
	vtkSmartPointer<vtkDoubleArray>		 Age;
  vtkSmartPointer<vtkDoubleArray>   Velocity;

	
  //
  int           UpdatePiece;
  int           UpdateNumPieces;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

private:
  vtkRamsesReader(const vtkRamsesReader&);  // Not implemented.
  void operator=(const vtkRamsesReader&);  // Not implemented.
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them 
	// in the output vector
	void AllocateAllRamsesVariableArrays(vtkIdType numBodies,
																			vtkPolyData* output);
//ETX

};
#endif
