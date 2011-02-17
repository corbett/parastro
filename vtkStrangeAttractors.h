/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkStrangeAttractors.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkStrangeAttractors - Read points from a Tipsy standard binary file
// .SECTION Description
// Read points from a Tipsy standard binary file. Fully parallel. Has ability
// to read in additional attributes from an ascii file, and to only load in
// marked particles but both these functions are serial only.
#ifndef __vtkStrangeAttractors_h
#define __vtkStrangeAttractors_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "vtkSmartPointer.h"
#include "tipsylib/ftipsy.hpp" // functions take Tipsy particle objects
#include <vtkstd/vector>

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

class VTK_EXPORT vtkStrangeAttractors : public vtkPolyDataAlgorithm
{
public:
  static vtkStrangeAttractors* New();
  vtkTypeRevisionMacro(vtkStrangeAttractors,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	vtkSetMacro(PrimesOnly,int);
	vtkGetMacro(PrimesOnly,int);
	
	
// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkStrangeAttractors();
  ~vtkStrangeAttractors();
	char* FileName;
	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

  int RequestData(vtkInformation*,vtkInformationVector**,
    vtkInformationVector*);

  vtkIdType                       ParticleIndex;
  vtkSmartPointer<vtkIdTypeArray> GlobalIds;
  vtkSmartPointer<vtkPoints>      Positions;
  vtkSmartPointer<vtkCellArray>   Vertices;
  vtkSmartPointer<vtkIntArray>   Primes;
	vtkSmartPointer<vtkFloatArray>   Velocity;
	int PrimesOnly;

private:
  vtkStrangeAttractors(const vtkStrangeAttractors&);  // Not implemented.
  void operator=(const vtkStrangeAttractors&);  // Not implemented.
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them 
	// in the output vector
	void AllocateAllVariableArrays(vtkIdType numBodies,
																			vtkPolyData* output);
	
	// used to get the next point color
	void nextColor(struct point *point);
	
//ETX

};
#endif
