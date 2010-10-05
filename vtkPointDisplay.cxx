/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkPointDisplay.cxx,v $

  Dummy reader file, does nothing exept displaing a single point at 0,0,0
  Shows basic setup of a reader file

=========================================================================*/
#include "vtkPointDisplay.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
//#include "sqlite3.h"

vtkCxxRevisionMacro(vtkPointDisplay, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkPointDisplay);

//----------------------------------------------------------------------------
// Constructor
vtkPointDisplay::vtkPointDisplay()
{
	this->FileName          = 0;		// init filename
	this->SetNumberOfInputPorts(0);   // set no of input files (0 is just fine)
}

//----------------------------------------------------------------------------
// Deconstructor
vtkPointDisplay::~vtkPointDisplay()
{

}

//----------------------------------------------------------------------------
// not exactly sure what it does...
int vtkPointDisplay::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{

// Stuff for doing it in parallel, leave it for the moment...
	//	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	//	// means that the data set can be divided into an arbitrary number of pieces
	//	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);

	return 1;
}

//----------------------------------------------------------------------------

int vtkPointDisplay::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{

// this is for opening a file..

	//  //
	//	// Make sure we have a file to read.
	//  //
	//  if(!this->FileName)
	//	  {
	//    vtkErrorMacro("A FileName must be specified.");
	//    return 0;
	//    }

// Setup the output

	vtkPolyData* output = vtkPolyData::GetData(outputVector);
	vtkSmartPointer<vtkPoints>   newPoints = vtkSmartPointer<vtkPoints>::New();
	// set up the vertices and cells, for direct display in paraview, not needed right now
	//vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	output->SetPoints(newPoints);
	//output->SetVerts(vertices);

/* -----
	Here comes the loop to read in and store the points
   ----- */

// Create the to be displayed points
	double dbCenterOfMass[3] = {0,0,0};

// Store the points
	newPoints->SetNumberOfPoints(1);
	newPoints->SetPoint(0, dbCenterOfMass);
	
	// Cells are needed for direct display of the points in paraview, otherwise a glyphfilter is needed.
	// let it be for the moment, do it next time..

	//vtkIdType *cells = vertices->WritePointer(1, 2);
	//cells[0] = 1;
	//cells[1] = 0;

 	return 1;
}
