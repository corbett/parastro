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

vtkCxxRevisionMacro(vtkSQLiteReader, "$Revision: 1.0.1 $");
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
	// set up the sql stuff
	sqlite3			*db;		// handle for db
	char			sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	int				sql_count;	// no of columns from the query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

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

// Prepare the query
	sql_query = ' ';

// Query the db
	sql_error = sqlite3_prepare_v2(db,
        "SELECT * FROM stat WHERE gid=1",
        1000, &res, &tail);

	if (sql_error != SQLITE_OK){
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	sql_count = sqlite3_column_count(res); //doesnt deliver the right no, but is >= than actual no
			// TODO: it it really?? only a assumption. need to check that!
	vtkErrorMacro("sql_count = " << sql_count);

// Create the points to be displayed 
	// alloc mem, is bit too much, but exact no isn't known yet, and souldn't make a big diff.
	newPoints->SetNumberOfPoints(sql_count);

	int i=0; //counts throu
	while (sqlite3_step(res) == SQLITE_ROW) {
		newPoints->SetPoint(i,
			sqlite3_column_double(res, 4),
			sqlite3_column_double(res, 5),
			sqlite3_column_double(res, 6));
		i++;
	}

	// set no of point to real value
	newPoints->SetNumberOfPoints(i);
	newPoints->Squeeze();
	
// Create the vertices (one point per vertex, for easy display)
	//vtkErrorMacro(vertices->Print);
	vtkIdType N = newPoints->GetNumberOfPoints();
	vtkErrorMacro("No of Points: " << N);
    vtkIdType *cells = vertices->WritePointer(N, N*2);
    for (vtkIdType i=0; i<N; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }

 	return 1;
}
