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
#include "vtkFloatArray.h"
#include <vector>
#include <sstream>


vtkCxxRevisionMacro(vtkSQLiteReader, "$Revision: 1.0.1 $");
vtkStandardNewMacro(vtkSQLiteReader);




//----------------------------------------------------------------------------
// Constructor
vtkSQLiteReader::vtkSQLiteReader()
{
	this->FileName          = 0;		// init filename
	this->SetNumberOfInputPorts(0);   // set no of input files (0 is just fine)

	this->db	= NULL;
	this->dataIsRead = false;

}

//----------------------------------------------------------------------------
// Deconstructor
vtkSQLiteReader::~vtkSQLiteReader()
{

}

/*----------------------------------------------------------------------------
reads some head data in, is called directly when loading the file (database), before apply is clicked (presents the infos for the window)
*/
int vtkSQLiteReader::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{

// Stuff for doing it in parallel, leave it for the moment...
	//	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	//	// means that the data set can be divided into an arbitrary number of pieces
	//	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);


// opening database
	if(!this->openDB(this->FileName)){return 0;}


// check for right format

// read in database header
	if(!this->ReadHeader(outputVector)){return 0;}

	return 1;
}

/*----------------------------------------------------------------------------
is called after clicking apply, reads the actual, selected data
*/
int vtkSQLiteReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
	//DEBUG implement this into gui later, for coding it's hardcoded..
	int displayid = (this->DisplaySnapshot);
	
	vtkPolyData * out = vtkPolyData::GetData(outputVector);

	if(!this->dataIsRead){
		//Set up output
		std::vector<vtkSmartPointer<vtkPolyData>> data(this->numSnaps); //big fat nasty vetor containing all data..
		
		this->data = data;
			
		//read in all the data (do this only once)
		readSnapshots(&(this->data));
	}

	//read in snap infos

	//generate tracks (do this only once)
	//TODO

	// set witch data to display
	out->SetPoints(this->data.at(displayid)->GetPoints());
	out->SetVerts(this->data.at(displayid)->GetVerts());

	return 1;
}

/*----------------------------------------------------------------------------
opens the database
	arguments:
		char * filename: path to db
	returns:
		int	errorcode (1 =ok)
	sets:
		sqlite3*	db:		database handle
*/
int vtkSQLiteReader::openDB(char* filename)
{
	if(!filename) //check for existing filename
	{
		vtkErrorMacro("A FileName must be specified.");
		return 0;
	}
	
	if (sqlite3_open(this->FileName, &db)) // open db, returns SQLITE_OK if successful
	{
		vtkErrorMacro("Can't open database: " + *sqlite3_errmsg(db));
		return 0;
	}

	//vtkErrorMacro("opened successfully: " << this->FileName)
	return 1;
}

/*----------------------------------------------------------------------------
reads the snapshots in
	assumes:
		db set and openend
	arguments:
		pointer to vtk... array, where to store data
		int numSnap		Number of snapshots
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots(std::vector<vtkSmartPointer<vtkPolyData>> * output)
{
//init vars

	//std::vector<vtkSmartPointer<vtkPolyData>> output(numSnap);
	//std::vector<vtkSmartPointer<vtkPolyData>>::iterator pvect;
	//pvect = output->begin();
	output->clear();


	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	vtkSmartPointer<vtkPoints> Position; // stores the coordinates of the points
	vtkSmartPointer<vtkCellArray> Vertices;
	vtkSmartPointer<vtkFloatArray> Velocity;
	vtkSmartPointer<vtkIntArray> qid;
	vtkSmartPointer<vtkIntArray> npart;


	// create temp data storage structre (simple list)
	struct tmp_data //elements of the list
	{
		int qid;
		int npart; // no of particles in this halo
		double pos[3]; //position
		double velo[3]; //velocity
		tmp_data* pnext; //pointer to next point in list
	};

	tmp_data * pfirst; //pointer to first element of list
	tmp_data * plast; //pointer to last element of list
	int count;
	//int snap = 1;

	for (int snap = 1; snap<=this->numSnaps; snap++) //loop over all snapshots
	{

		//init datastructre
		tmp_data * dummydata = new tmp_data;
		pfirst = dummydata;
		plast = dummydata;
		count = 0;

		// Prepare the query
		sql_query = "SELECT * FROM stat WHERE snap_id=" + Int2Str(snap);

		//vtkErrorMacro("SQL query: " + sql_query);

		sql_error = sqlite3_prepare_v2(db,
			sql_query,
			1000, &res, &tail);

		if (sql_error != SQLITE_OK)
		{
			vtkErrorMacro("Error with sql query! Error: " << sql_error);
			return 0;
		}

		while (sqlite3_step(res) == SQLITE_ROW)
		{

			tmp_data * data = new tmp_data;

			data->qid = sqlite3_column_int(res, 1); //qid
			data->npart = sqlite3_column_int(res, 2); //npart
			// sqlite3_column_double(res, 3); //nvpart
			
			data->pos[0] = sqlite3_column_double(res, 4); //xc
			data->pos[1] = sqlite3_column_double(res, 5); //yc
			data->pos[2] = sqlite3_column_double(res, 6); //zc
			
			data->velo[0] = sqlite3_column_double(res, 7); //vxc
			data->velo[0] = sqlite3_column_double(res, 8); //vyc
			data->velo[0] = sqlite3_column_double(res, 9); //vzc
			// ... and other stuff

			// orga stuff
			data->pnext=NULL;
			plast->pnext = data;
			plast = data;
			count++;

			//debug
			//if ( snap==30 ) {vtkErrorMacro(" " << data->pos[0]);}
		}

		//store the data
		Position = vtkSmartPointer<vtkPoints>::New();
			Position->SetNumberOfPoints(count);
		Vertices = vtkSmartPointer<vtkCellArray>::New();
		Velocity = vtkSmartPointer<vtkFloatArray>::New();
			Velocity->SetNumberOfComponents(3);
			Velocity->SetNumberOfTuples(count);
			Velocity->SetName("Velocity");
		qid = vtkSmartPointer<vtkIntArray>::New();
			qid->SetNumberOfComponents(1);
			qid->SetNumberOfTuples(count);
			qid->SetName("qid");
		npart = vtkSmartPointer<vtkIntArray>::New();
			npart->SetNumberOfComponents(1);
			npart->SetNumberOfTuples(count);
			npart->SetName("npart");

		vtkIdType * cells = Vertices->WritePointer(count, count*2);

		tmp_data * actdata = pfirst->pnext;
		
		for (int n = 0; n<count;n++)
		{
			Position->SetPoint(n, actdata->pos);
			Velocity->InsertTuple(n, actdata->velo);
			qid->InsertValue(n,actdata->qid);
			npart->InsertValue(n,actdata->npart);

			cells[n*2] = 1;
			cells[n*2+1] = n;

			actdata = actdata->pnext;
		}
		
		output->push_back(vtkSmartPointer<vtkPolyData>::New());
		(*output)[snap-1]->SetPoints(Position);
		(*output)[snap-1]->SetVerts(Vertices);

		//output[snap]->GetPointData()->AddArray(Velocity);
		//output[snap]->GetPointData()->AddArray(qid);
		//output[snap]->GetPointData()->AddArray(npart);

		//vtkSmartPointer<vtkPolyData> tmp = vtkSmartPointer<vtkPolyData>::New();
		//tmp->SetPoints(Position);
		//tmp->SetVerts(Vertices);


		//*pvect = output->insert(pvect,tmp);
		//pvect++;

		//snap++;
	}

	this->dataIsRead = true;
	return 1;
}

/*----------------------------------------------------------------------------
Demonstration of reading in particles
*/
int vtkSQLiteReader::RequestDataDemo(vtkInformationVector* outputVector)
{

//init vars
	// set up the sql stuff
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




// Prepare the query
	sql_query = "SELECT * FROM stat WHERE gid=1";
	//sql_query = this->GetSqlQuery();

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

/*	--------------------------
	version 2:
	saves data first in simple datastructre
	check for speed differences with version 1.. possibly this is
	easier to implent parallelly. already set sqlite3 to use multiply
	threads (check sqlite3.c: SQLITE_THREADSAFE 2), see doc here:
	(http://www.sqlite.org/threadsafe.html) and possibly this disc:
	(http://osdir.com/ml/sqlite-users/2010-03/msg00209.html) for help
	later on.
	-------------------------- */

	//vtkErrorMacro("threadsafe: " << sqlite3_threadsafe()); // display threading mode

	// create temp data storage structre (simple list)
	struct tmp_pnt //elements of the list
	{
		double data[3]; //data array
		tmp_pnt* pnext; //pointer to next point in list
	};

	/* should first store the list itself, but i use 2 pointers and a counter directly instead.. could be deleted..
	struct tmp_list //the list itself
	{
		tmp_pnt * pfirst; //pointer to first element of list
		tmp_pnt * plast; //pointer to last element of list
		int count;
	};
	*/

	//statt struct einfach direkt 2 pointer brauchen...
	tmp_pnt * pfirst; //pointer to first element of list
	tmp_pnt * plast; //pointer to last element of list
	int count = 0;

	// add a first dummy point
	tmp_pnt * firstpnt = new tmp_pnt;
	pfirst = firstpnt; 
	plast = firstpnt; 

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


/*----------------------------------------------------------------------------
Converts a Integer to a string
*/
vtkStdString vtkSQLiteReader::Int2Str(int number)
{
	std::stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

/*----------------------------------------------------------------------------
Reads the header data from db
	assumes:
		opened database, db set
	arguments:
		vtkInformationVector Output information
	returns:
		int	errorcode (1 =ok)
	sets:
		numSnap		Number of snapshots

*/
int vtkSQLiteReader::ReadHeader(vtkInformationVector* outputVector)
{
	// set up the sql stuff
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	vtkStdString sql_query = "SELECT * FROM header";
	
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	if (sqlite3_step(res) == SQLITE_ROW) //there should be only one row..
	{
		this->numSnaps = sqlite3_column_int(res, 9);
	}

	return 1;
}

/*----------------------------------------------------------------------------
-- NOT USED --
Sets up an SQL Query
	assumes:
		opened database, db set
	arguments:
		vtkStdString sql_query the query command
		sqlite3_stmt *res pointer to result
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::SQLQuery(vtkStdString sql_query, sqlite3_stmt * res)
{
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	sql_error = sqlite3_prepare_v2(this->db, sql_query, 1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	return 1;
}