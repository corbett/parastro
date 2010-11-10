/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.cxx,v $
=========================================================================*/
#include "vtkSQLiteReader.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"
#include "vtkArrayData.h"
#include "vtkLookupTable.h"

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


// TODO check for right format
	/* implement this, so i can check for right database format in the beginning and omit all the errorhandling code with database queries.. should give some speed..*/

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

	///* OLD STUFF, USED WITH readSnapshot(..)
	//if(!this->dataIsRead){
	//	//Set up output
	//	std::vector<vtkSmartPointer<vtkPolyData>> data(this->numSnaps); //big fat nasty vetor containing all data..
	//	
	//	this->data = data;
	//		
	//	//read in all the data (do this only once)
	//	readSnapshots(&(this->data));
	//	
	//	//read in snap infos
	//	ReadSnapshotInfo();

	//	//read in the track info
	//	ReadTracks();
	//}
	//// --END OLD STUFF --- */
	//
	////new stuff here:

	//if(!this->dataIsRead){
	//	//read in all the data (do this only once)
	//	readSnapshots2();
	//	
	//	//read in snap infos
	//	ReadSnapshotInfo();

	//	//read in the track info
	//	ReadTracks();
	//}
	//
	////generate tracks (do this only once)
	////GenerateTracks();
	////TODO
	//
	////generate all the lines
	////vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	////CollectLines(&lines);
	//
	//// set witch data to display
	///* old stuff again..
	//out->SetPoints(this->data.at(displayid)->GetPoints());
	//out->SetVerts(this->data.at(displayid)->GetVerts());
	//out->SetLines(lines);
	//*/
	//out->SetPoints(this->data2.at(displayid).coord);
	//out->SetVerts(this->data2.at(displayid).cells);
	////out->SetLines(lines);

	////vtkSQLiteReader::GenerateOutput(out);


	// --- v3 here ---

	readSnapshots3();
	readSnapshotInfo3();
	readTracks3();
	generateColors();

	out->SetPoints(this->Position);
	out->SetVerts(this->Cells);
	out->SetLines(this->Tracks);
	out->GetPointData()->AddArray(this->Velocity);
	out->GetPointData()->AddArray(this->Qid);
	out->GetPointData()->AddArray(this->SnapId);
	out->GetPointData()->AddArray(this->RVir);
	out->GetPointData()->AddArray(this->TrackId);


	//out->GetPointData()->AddArray(this->colors);
	out->GetPointData()->SetScalars(this->colors);

	
	//vtkSmartPointer<vtkLookupTable> LuT = vtkSmartPointer<vtkLookupTable>::New();
	//LuT->SetTableRange(1,100);
	//double alpha [] = {0,1,0};
	//LuT->SetAlphaRange(alpha);
	//LuT->Build();

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
NOT USED ANYMORE, IS HERE FOR BACKUP PURPOSES
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

	//TODO maybee i've to put this in header, otherwise theyre destroyed
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
		
		//vtkSmartPointer<vtkPolyData> tmp = vtkSmartPointer<vtkPolyData>::New();

		
		output->push_back(vtkSmartPointer<vtkPolyData>::New());
		output->at(snap-1)->SetPoints(Position);
		output->at(snap-1)->SetVerts(Vertices);

		//(vtkSmartPointer<vtkPolyData>)output[snap]->GetPointData()->AddArray(Velocity);
		output->at(snap-1)->GetPointData()->AddArray(qid);
		output->at(snap-1)->GetPointData()->AddArray(npart);
		output->at(snap-1)->GetPointData()->AddArray(Velocity);

		//tmp->GetPointData()->AddArray(Velocity);
	
		//vtkSmartPointer<vtkPolyData> tmp = vtkSmartPointer<vtkPolyData>::New();
		//tmp->SetPoints(Position);
		//tmp->SetVerts(Vertices);



	}

	this->dataIsRead = true;
	return 1;
}
/*----------------------------------------------------------------------------
OUTDATED
reads the snapshots in
	assumes:
		db set and openend
	arguments:
		none
	sets:
		vector<> this->data where alllll the data is..
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots2()
{
	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	vtkSmartPointer<vtkPoints> tmp_points;
	vtkSmartPointer<vtkCellArray> tmp_cells;

	this->data2.clear();
	this->totNumPoints = 0;

	for (int snap = 1; snap<=this->numSnaps; snap++) //loop over all snapshots
	{
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

		int count = 0;
		
		velocity tmpv;
		std::vector<velocity> tmp_velo;
		tmp_velo.clear();

		tmp_points = vtkSmartPointer<vtkPoints>::New();

		while (sqlite3_step(res) == SQLITE_ROW)
		{
			// sqlite3_column_int(res, 1); //qid
			// sqlite3_column_int(res, 2); //npart
			// sqlite3_column_double(res, 3); //nvpart
			
			//maybee use here first a vector to store the item and alloc later and insert then fast
			tmp_points->InsertNextPoint(
				sqlite3_column_double(res, 4), //xc
				sqlite3_column_double(res, 5), //yc
				sqlite3_column_double(res, 6)); //zc
			
			tmpv.vx = sqlite3_column_double(res, 7); //vxc
			tmpv.vy = sqlite3_column_double(res, 8); //vyc
			tmpv.vz = sqlite3_column_double(res, 9); //vzc
			tmp_velo.push_back(tmpv);
			
			count++;
		}

		tmp_cells = vtkSmartPointer<vtkCellArray>::New();
		vtkIdType * cells = tmp_cells->WritePointer(count, count*2);

		for (int n = 0; n<count;n++)
		{
			cells[n*2] = 1;
			cells[n*2+1] = n;
		}

		snapshot actSnapshot;
		actSnapshot.coord = tmp_points;
		actSnapshot.cells = tmp_cells;
		actSnapshot.velo = tmp_velo;

		this->totNumPoints = this->totNumPoints + count;
		this->data2.push_back(actSnapshot);
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
/*----------------------------------------------------------------------------
Reads in the infos for the snapshots
	assumes:
		opened database, db set
	sets:
		this->snapinfo	information aout snapshots
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::ReadSnapshotInfo()
{
	this->snapinfoVector.clear();

	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM snapinfo";

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
		snapinfo tmp = {
			sqlite3_column_int(res, 0), //snap_id
			sqlite3_column_double(res, 2), //redshift
			sqlite3_column_double(res, 3), //time
			sqlite3_column_int(res, 5)}; //npart
			
		this->snapinfoVector.push_back(tmp);
	}
	return 1;
}
/*----------------------------------------------------------------------------
OUTDATED
Reads in the track information
	assumes:
		opened database, db set
	sets:
		std::vector<track> this->trackVector	information about tracks
		int this->numTracks		tot count of tracks
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::ReadTracks()
{
	//TODO all this...
	this->trackVector.clear();

	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	bool goOn = true;
	int tracksCount = 0;
	int pointCount;

	while (goOn)
	{
		sql_query = "SELECT * FROM tracks WHERE id=" + Int2Str(tracksCount+1);
		sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

		if (sql_error != SQLITE_OK)
		{
			vtkErrorMacro("Error with sql query! Error: " << sql_error);
			return 0;
		}
		//bool doIt = true;
		track tmpTrack;
		pointCount = 0;
		
		while (sqlite3_step(res) == SQLITE_ROW)
		{
			pointCount++;
			
			trackPoint tmpTrackPoint = {
				sqlite3_column_int(res, 1), //snap_id
				sqlite3_column_int(res, 2)}; //qid
					
			tmpTrack.point.push_back(tmpTrackPoint);
		}
		if (pointCount>0)
		{
			tmpTrack.noOfPoints = pointCount;
			this->trackVector.push_back(tmpTrack);
		}
		else
		{
			goOn = false;
		}
		
		tracksCount++;
	}

	this->numTracks = tracksCount-1;
	return 1;
}
/*----------------------------------------------------------------------------
OUTDATED, NOT USED ANYMORE
Generates the tracks
	assumes:
		opened database, db set
	sets:
		this->trackinfo	information about tracks
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::GenerateTracks()
{
/* DEMO for doing the lines
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	
	for (int i = 0;i<3;i++){
		vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
		nextLine->GetPointIds()->SetNumberOfIds(3);
		nextLine->GetPointIds()->SetId(0,i*3+0);
		nextLine->GetPointIds()->SetId(1,i*3+1);
		nextLine->GetPointIds()->SetId(2,i*3+2);

		lines->InsertNextCell(nextLine);
	}

	//and then outside
	out->SetLines(lines);
*/

	track * actTrack;
	vtkSmartPointer<vtkPolyLine> actLine;
	vtkSmartPointer<vtkDataArray> actQid;

	int qid;
	int snap_id;

	for (int trackno = 0; trackno < this->numTracks; trackno++)
	{
		actTrack = &(this->trackVector.at(trackno));
		actLine = vtkSmartPointer<vtkPolyLine>::New();
		actLine->GetPoints()->SetNumberOfPoints(numSnaps);
				
		//
		// Continue HERE
		// try to add simple lines first, probably this part here sucks..
		//
		//
		
		for(int i = 0; i<this->numSnaps;i++){
			
			qid = actTrack->point.at(i).qid;
			snap_id = (actTrack->point.at(i).snap_id)-1;
			
			// gets pointer to double array with data
			double * ptr = this->data2.at(snap_id).coord->GetData()->GetTuple3(qid);
			
			//sets the pointer to data
			actLine->GetPoints()->GetData()->SetTuple(i, ptr);
			
		}
		actTrack->line = actLine;
	}

	return 1;
}

/*----------------------------------------------------------------------------
OUTDATED, NOT USED ANYMORE
Collects all the tracks in trackVector and generates one vtkCellArray to display them
	assumes:
		trackVector filled up
	sets:
		this->trackinfo	information about tracks
	arguments:
		vtkSmartPointer<vtkCellArray> * Pointer to cell array to fill in the collected lines
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::CollectLines(vtkSmartPointer<vtkCellArray> * lines)
{

	for (int trackid = 0;trackid<this->numTracks;trackid++)
	{
		(*lines)->InsertNextCell(this->trackVector.at(trackid).line);
	}
	return 1;
}


/*----------------------------------------------------------------------------
Collects all the tracks in trackVector and generates one vtkCellArray to display them
Collects all Points
	assumes:
		data2 filled with points
		trackVector filled up
	sets:
		this->trackinfo	information about tracks
	arguments:
		vtkPolyData * out pointer to outputdata
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::GenerateOutput(vtkPolyData * out)
{

/*

	// START OVER, THIS SUCKS...

	// prepare datastructre

	// loop over all elements in data2 and put pointers to 

	
	vtkSmartPointer<vtkArrayData> arrayData0 = vtkSmartPointer<vtkArrayData>::New();
	arrayData0->AddArray(data2.at(0).coord->Data);
	arrayData0->Update();
	arrayData0->

	vtkSmartPointer<vtkPoints> tmp_points = vtkSmartPointer<vtkPoints>::New();

	tmp_points->SetNumberOfPoints(this->totNumPoints);


	vtkSmartPointer<vtkCellArray> tmp_cells;
	//vtkSmartPointer<vtkDataArray> tmp_pointdata = vtkSmartPointer<vtkDataArray>::New();
	vtkSmartPointer<vtkCellArray> tmp_lines;

	//tmp_pointdata->SetNumberOfComponents(3);

	// collect all points
	int totPoints = 0;
	int numPoints = 0;

	for (int i = 0; i<this->numSnaps;i++){
		numPoints = data2.at(i).coord->GetNumberOfPoints();
		//totPoints =+ numPoints;
		//out->GetPointData()->SetNumberOfTuples(totPoints);   //->Resize(totPoints);
		for (int j = 0; j<numPoints; j++){
			double tmpdata[3];
			data2.at(i).coord->GetPoint(j,tmpdata);
			//out->GetPoints()->GetData()->SetTuple(i,ptr);
			////tmp_pointdata->SetTuple(i,data2.at(i).coord->GetData()->GetTuple3(j));
			//tmp_points->InsertNextPoint(*ptr);
			out->GetPoints()->InsertNextPoint(tmpdata);
		}
	}

	//out->GetPoints()->SetData(tmp_pointdata);

	//out->SetPoints(tmp_points);
	//out->SetVerts(tmp_cells);
	//out->SetLines(tmp_lines);
	*/
	return 1;
}

/*----------------------------------------------------------------------------
reads the snapshots in, next try
	assumes:
		db set and openend
	arguments:
		none
	sets:
		
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots3()
{
// prepare the variables

	this->Position = vtkSmartPointer<vtkPoints>::New();
	this->Velocity = vtkSmartPointer<vtkFloatArray>::New();
	this->Cells = vtkSmartPointer<vtkCellArray>::New();

	this->Qid = vtkSmartPointer<vtkIdTypeArray>::New();
	this->SnapId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->RVir = vtkSmartPointer<vtkFloatArray>::New();

	this->Velocity->SetNumberOfComponents(3);
	this->Qid->SetNumberOfComponents(1);
	this->SnapId->SetNumberOfComponents(1);
	this->RVir->SetNumberOfComponents(1);

	this->Velocity->SetName("Velocity");
	this->Qid->SetName("QId");
	this->SnapId->SetName("SnapId");
	this->RVir->SetName("RVir");

	this->SnapInfo3.clear();

// sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

// Prepare the query
	sql_query = "SELECT * FROM stat";

// Query the db
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	int count = 0;
	int snap_id=1;
	int old_snap_id = 1;
	int qid = 0;

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		old_snap_id = snap_id;
		snap_id = sqlite3_column_double(res, 0);

		if (snap_id != old_snap_id)
		{
			SnapshotInfo3 tmp;
			tmp.Offset = count-qid;
			tmp.lenght = qid;
			this->SnapInfo3.push_back(tmp);
		}

		qid = sqlite3_column_double(res, 1);

		this->Position->InsertNextPoint(
			sqlite3_column_double(res, 4),
			sqlite3_column_double(res, 5),
			sqlite3_column_double(res, 6));
		this->Velocity->InsertNextTuple3(
			sqlite3_column_double(res, 7),
			sqlite3_column_double(res, 8),
			sqlite3_column_double(res, 9));
		this->Qid->InsertNextTuple1(qid);
		this->SnapId->InsertNextTuple1(snap_id);
		this->RVir->InsertNextTuple1(sqlite3_column_double(res, 11));

		count++;
	}

	SnapshotInfo3 tmp;
	tmp.Offset = count-qid;
	tmp.lenght = qid;
	this->SnapInfo3.push_back(tmp);

	this->nTracks3 = (int)this->SnapInfo3.size();

	this->nParticles3 = count;

	// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = this->Position->GetNumberOfPoints();
	//vtkErrorMacro("No of Points: " << N);
	vtkIdType *cells = this->Cells->WritePointer(N, N*2);
    for (vtkIdType i=0; i<N; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }

	this->dataIsRead = true;
	return 1;
}

int vtkSQLiteReader::readTracks3(){

	this->Tracks = vtkSmartPointer<vtkCellArray>::New();

	this->TrackId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->TrackId->SetName("TrackId");
	this->TrackId->SetNumberOfComponents(1);
	this->TrackId->SetNumberOfTuples(this->nParticles3);

	//init everything to 0 (0: halo belongs to no track)
	for (int i = 0; i<this->nParticles3;i++)
	{
		this->TrackId->InsertTuple1(i,0);
	}

	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	int snap_id;
	int qid;
	int offset;

	for (int TrackNo = 0; TrackNo < this->nTracks3; TrackNo++) //
	{
		sql_query = "SELECT * FROM tracks WHERE id=" + Int2Str(TrackNo+1);
		sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

		if (sql_error != SQLITE_OK)
		{
			vtkErrorMacro("Error with sql query! Error: " << sql_error);
			return 0;
		}
		
		vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
		
		while (sqlite3_step(res) == SQLITE_ROW)
		{
			snap_id = sqlite3_column_int(res, 1);
			qid = sqlite3_column_int(res, 2);

			//generate the lines
			if (qid<1) {continue;} //check for empty fields
			offset = this->SnapInfo3.at(snap_id-1).Offset + qid;

			nextLine->GetPointIds()->InsertNextId(offset-1);

			//fill the TrackId vector
			this->TrackId->InsertTuple1(offset-1,TrackNo+1);
		}

		this->Tracks->InsertNextCell(nextLine);
	}

	return 1;

}
/*----------------------------------------------------------------------------
Reads in the infos for the snapshots
	assumes:
		opened database, db set
	sets:
		this->snapinfo	information aout snapshots
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::readSnapshotInfo3()
{
	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM snapinfo";

	//vtkErrorMacro("SQL query: " + sql_query);

	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	int snap_id;
	
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		snap_id = sqlite3_column_int(res, 0);
		this->SnapInfo3.at(snap_id-1).redshift = sqlite3_column_double(res, 2);
		this->SnapInfo3.at(snap_id-1).time = sqlite3_column_double(res, 3);
		this->SnapInfo3.at(snap_id-1).npart = sqlite3_column_int(res, 5);
	}
	return 1;
}

/*----------------------------------------------------------------------------
Generates the colorvector for the data

at the moment, it just generates the color depending on snapid, but
this can later be made interactive..
	assumes:
		all data read in
	sets:
		this->color		color vector
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::generateColors()
{

	this->colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->colors->SetName("Colors");
	this->colors->SetNumberOfComponents(3);
	this->colors->SetNumberOfTuples(this->nParticles3);

	//int color, index;
	//int n;

	//int offset = 0;
	//for (int i = 0; i<this->numSnaps; i++)
	//{
	//	int l = this->SnapInfo3.at(i).lenght;
	//	for (n = 0; n<l;n++)
	//	{
	//		color = 255*i/(this->numSnaps-1);
	//		index = offset+n;

	//		this->colors->InsertTuple3(index,
	//			color,
	//			0,
	//			0);
	//	}
	//	offset =+ n;
	//}

	int snapid, trackid;
	int red, green, blue;

	for (int i = 0; i<this->nParticles3;i++)
	{
		snapid = this->SnapId->GetTuple1(i);
		trackid = this->TrackId->GetTuple1(i);

		red = 255 * (snapid-1) / (this->numSnaps-1);

		if (trackid == 2)
		{
			green =255;
		} else
		{
			green = 0;
		}

		blue = 0;

		this->colors->InsertTuple3(i,
				red,
				green,
				blue);
	}

	return 1;
}