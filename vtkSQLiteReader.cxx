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
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	//int				sql_count;	// no of columns from the query
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
	sql_query = "SELECT * FROM stat WHERE gid=1";

// Query the db
	sql_error = sqlite3_prepare_v2(db,
        sql_query,
        1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

// --------------------------
	// version 1:
	/* uses insertpoint, possibly slow because everytime range check */
/*
	int i = 0;
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		newPoints->InsertPoint(i,
			sqlite3_column_double(res, 4),
			sqlite3_column_double(res, 5),
			sqlite3_column_double(res, 6));
		/*vtkErrorMacro("coords: "
			<< sqlite3_column_double(res, 4) << ";"
			<< sqlite3_column_double(res, 5) << ";"
			<< sqlite3_column_double(res, 6) << ";");
		i++;
	}

	newPoints->Squeeze();

*/
// --------------------------
	// version 2:
	/* saves data first in simple datastructre
	check for speed differences with version 1.. possibly this is
	easier to implent parallelly. already set sqlite3 to use multiply
	threads (check sqlite3.c: SQLITE_THREADSAFE 2), see doc here:
	(http://www.sqlite.org/threadsafe.html) and possibly this disc:
	(http://osdir.com/ml/sqlite-users/2010-03/msg00209.html) for help
	later on.
	*/

	//vtkErrorMacro("threadsafe: " << sqlite3_threadsafe()); // display threading mode

	// create temp data storage structre (simple list)
	//vtkErrorMacro("structdefs");
	struct tmp_pnt //elements of the list
	{
		double data[3]; //data array
		tmp_pnt* pnext; //pointer to next point in list
		//int count;
	};

	/* should first store the list itself, but i use 2 pointers and a counter directly instead.. could be deleted..
	struct tmp_list //the list itself
	{
		tmp_pnt * pfirst; //pointer to first element of list
		tmp_pnt * plast; //pointer to last element of list
		int count;
	};

	vtkErrorMacro("init mytmplist");
	tmp_list my_tmp_list = {NULL,NULL,0};*/

	//statt struct einfach direkt 2 pointer brauchen...
	tmp_pnt * pfirst; //pointer to first element of list
	tmp_pnt * plast; //pointer to last element of list
	int count = 0;

	// add a first dummy point
	tmp_pnt * firstpnt = new tmp_pnt;
	pfirst = firstpnt; 
	plast = firstpnt; 

	//vtkErrorMacro("add entry to list");

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		tmp_pnt * pnt = new tmp_pnt;
		//double (*data)[3] = new double;
		double data[3] = {sqlite3_column_double(res, 4),
			sqlite3_column_double(res, 5),
			sqlite3_column_double(res, 6)};
		pnt->data[0] = data[0];
		pnt->data[1] = data[1];
		pnt->data[2] = data[2];
		pnt->pnext = NULL;

		plast->pnext = pnt;
		plast = pnt;
		//vtkErrorMacro("read index: " << count << " cords: " << (pnt->data)[0]);
		count++;
	}
	//vtkErrorMacro("count tot: " << count);

	newPoints->SetNumberOfPoints(count);
	tmp_pnt * actpnt = pfirst->pnext;
	
	for (int n = 0; n < count; n++)
	{
		newPoints->SetPoint(n,(actpnt->data));
		
		// for debugging
		//double (* pnt1)[3] = &(actpnt->data);
		//vtkErrorMacro("index: " << n << " cords: " << **pnt1 <<" "<<  *(*(pnt1)+1) <<" "<< *(*(pnt1)+2));

		actpnt = actpnt->pnext;
	}

	newPoints->Squeeze();

// --------------------------

// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = newPoints->GetNumberOfPoints();
	//vtkErrorMacro("No of Points: " << N);
    vtkIdType *cells = vertices->WritePointer(N, N*2);
    for (vtkIdType i=0; i<N; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }

// cleaning up mem
	//TODO Clean up everything...
	sqlite3_close(db); //closing db handle (maybe do this earlier??)

 	return 1;
}
