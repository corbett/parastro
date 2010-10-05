/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.cxx,v $
=========================================================================*/
#include "vtkSQLiteReader.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
//#include "sqlitelib/sqlite3.h"

vtkCxxRevisionMacro(vtkSQLiteReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkSQLiteReader);

//----------------------------------------------------------------------------
// Constructor
vtkSQLiteReader::vtkSQLiteReader()
{
	this->FileName          = 0;		// init filename
	this->SetNumberOfInputPorts(0);   // set no of input files (0 is just fine)
}

//----------------------------------------------------------------------------
// Deconstructor
vtkSQLiteReader::~vtkSQLiteReader()
{

}

//----------------------------------------------------------------------------
// not exactly sure what it does...
int vtkSQLiteReader::RequestInformation(
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

int vtkSQLiteReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{

//init vars
	sqlite3 *db; //handle for db

	// sets up output
	vtkPolyData* output = vtkPolyData::GetData(outputVector);

	//sets up the output of the points
	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	output->SetPoints(newPoints);

	// set up the vertices and cells, for direct display in paraview	, not needed right now
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	output->SetVerts(vertices);


// opening database

	if(!this->FileName) //check for existing filename
	{
		vtkErrorMacro("A FileName must be specified.");
		return 0;
	}
	
	if (sqlite3_open(this->FileName, &db)) // open db, returns SQLITE_OK if successful
	{
		vtkErrorMacro("Can't open database: " + *sqlite3_errmsg(db));
	}
	//vtkErrorMacro("opened successfully: " << this->FileName)

/* -----
	Here comes the loop to read in and store the points
   ----- */

// Create the to be displayed points
	double dbCenterOfMass[3] = {0,0,0};
	double dbCenterOfMass2[3] = {1,0,0};

// Store the points
	newPoints->SetNumberOfPoints(2);
	newPoints->SetPoint(0, dbCenterOfMass);
	newPoints->SetPoint(1, dbCenterOfMass2);
	
	// Cells are needed for direct display of the points in paraview, otherwise a glyphfilter is needed.
	
	vtkErrorMacro(vertices->Print)

	vtkIdType *cells = vertices->WritePointer(1, 2);
	cells[0] = 1;
	cells[1] = 2;
	cells[2] = 0;

 	return 1;
}
