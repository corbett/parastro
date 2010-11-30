/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.cxx,v $

  
	todo:
	- crashes when changing parameters... 1. run works  
	- check with big data
	- pointselection for display
	- optimizing of collision detection..
  
  
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
#include "vtkPolyDataMapper.h"
#include "vtkSMProxy.h"
#include "vtkSMProxyProperty.h"



#include <vector>
#include <list>
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

	//fasse freie vars in structs zusammen
	Gui.DisplayOnlySelectedData = &(this->DisplayOnlySelectedData); 
	Gui.DisplaySelected = &(this->DisplaySelected); 
	Gui.DisplaySelectedSnapshot = &(this->DisplaySelectedSnapshot); 
	Gui.DisplaySelectedSnapshotNr = &(this->DisplaySelectedSnapshotNr); 
	Gui.DisplaySelectedTrack = &(this->DisplaySelectedTrack); 
	Gui.DisplaySelectedTrackNr = &(this->DisplaySelectedTrackNr); 
	Gui.DisplayCalculated = &(this->DisplayCalculated);
	Gui.CalculationImpactParameter = &(this->CalculationImpactParameter);

	// use this only temp, later delete the free variables, only use struct
	dataInfo.nParticles = &(this->nParticles3);
	dataInfo.nSnapshots = &(this->nSnapshots);
	dataInfo.nTracks = &(this->nTracks3);
	dataInfo.nSelectedParticles = -1;
	dataInfo.nSelectedTracks = -1;
	dataInfo.selectedPoints = vtkSmartPointer<vtkIdTypeArray>::New();
	dataInfo.selectedSnapshots = vtkSmartPointer<vtkIdTypeArray>::New();
	dataInfo.selectedTracks = vtkSmartPointer<vtkIdTypeArray>::New();
	dataInfo.selectedPoints->SetNumberOfComponents(1);
	dataInfo.selectedSnapshots->SetNumberOfComponents(1);
	dataInfo.selectedTracks->SetNumberOfComponents(1);


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
	if(!this->ReadHeader()){return 0;}

	return 1;
}

/*----------------------------------------------------------------------------
is called after clicking apply, reads the actual, selected data
*/
int vtkSQLiteReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{

	vtkPolyData * out = vtkPolyData::GetData(outputVector);


	if (!this->dataIsRead)
		// only read if it's not been read bevore
	{
		readSnapshots();
		readSnapshotInfo();
		readTracks();

		this->dataIsRead = true;
	}
	
	if(*this->Gui.DisplaySelected)
	{
		if(*this->Gui.DisplaySelectedSnapshot){
			dataInfo.selectedSnapshots->InsertNextTuple1(*Gui.DisplaySelectedSnapshotNr);
		}
		if(*this->Gui.DisplaySelectedTrack){
			dataInfo.selectedTracks->InsertNextTuple1(*Gui.DisplaySelectedTrackNr);
		}
	}
	else if (*this->Gui.DisplayCalculated)
	{
		doCalculations(); //fills dataInfo.selected* fields
	}

	if (*this->Gui.DisplayOnlySelectedData)
	{
		selectPoints();
	}
	else 
	{
		DataSelection.Position = this->Position;
		DataSelection.Cells = this->Cells;
		DataSelection.Tracks = this->Tracks;
		DataSelection.Velocity = this->Velocity;
		DataSelection.GId = this->GId;
		DataSelection.SnapId = this->SnapId;
		DataSelection.RVir = this->RVir;
		DataSelection.TrackId = this->TrackId;
	}

	// update the colors anyways

	out->SetPoints(DataSelection.Position);
	out->SetVerts(DataSelection.Cells);
	out->SetLines(DataSelection.Tracks);
	out->GetPointData()->AddArray(DataSelection.Velocity);
	out->GetPointData()->AddArray(DataSelection.GId);
	out->GetPointData()->AddArray(DataSelection.SnapId);
	out->GetPointData()->AddArray(DataSelection.RVir);
	out->GetPointData()->AddArray(DataSelection.TrackId);
	generateColors();
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
int vtkSQLiteReader::ReadHeader()
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
		this->nSnapshots = sqlite3_column_int(res, 9);
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
reads the snapshots in (reads all the points, and the according data, generates the cells)
	assumes:
		db set and openend
	arguments:
		none
	sets:
		this->Position
		this->Velocity
		this->Cells
		this->Qid
		this->SnapId
		this->RVir
		this->SnapInfo3		stores information about the snapshots
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots()
{
// prepare the variables

	this->Position = vtkSmartPointer<vtkPoints>::New();
	this->Velocity = vtkSmartPointer<vtkFloatArray>::New();
	this->Cells = vtkSmartPointer<vtkCellArray>::New();

	this->GId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->SnapId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->RVir = vtkSmartPointer<vtkFloatArray>::New();

	this->Velocity->SetNumberOfComponents(3);
	this->GId->SetNumberOfComponents(1);
	this->SnapId->SetNumberOfComponents(1);
	this->RVir->SetNumberOfComponents(1);

	this->Velocity->SetName("Velocity");
	this->GId->SetName("GId");
	this->SnapId->SetName("SnapId");
	this->RVir->SetName("RVir");

	this->SnapInfo3.clear();

// sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

// Prepare the query
	sql_query = "SELECT * FROM stat ORDER BY snap_id";

// Query the db
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	int count = 0;
	int snap_id=0;
	int old_snap_id = 0;
	int GId = 0;

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		old_snap_id = snap_id;
		snap_id = sqlite3_column_double(res, 0);

		if (snap_id != old_snap_id)
		{
			SnapshotInfo tmp;
			tmp.Offset = count-GId-1;
			tmp.lenght = GId+1;
			this->SnapInfo3.push_back(tmp);
		}

		GId = sqlite3_column_double(res, 1);

		this->Position->InsertNextPoint(
			sqlite3_column_double(res, 4),
			sqlite3_column_double(res, 5),
			sqlite3_column_double(res, 6));
		this->Velocity->InsertNextTuple3(
			sqlite3_column_double(res, 7),
			sqlite3_column_double(res, 8),
			sqlite3_column_double(res, 9));
		this->GId->InsertNextTuple1(GId);
		this->SnapId->InsertNextTuple1(snap_id);
		this->RVir->InsertNextTuple1(sqlite3_column_double(res, 11));

		count++;
	}

	SnapshotInfo tmp;
	tmp.Offset = count-GId-1;
	tmp.lenght = GId+1;
	this->SnapInfo3.push_back(tmp);

	//this->nTracks3 = (int)10;
	//vtkErrorMacro("No of tracks: " << this->nTracks3);

	this->nSnapshots = 	(int)this->SnapInfo3.size();
	vtkErrorMacro("No of snapshots: " << this->nSnapshots);

	this->nParticles3 = count;
	vtkErrorMacro("No of particles: " << this->nParticles3);

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
/*----------------------------------------------------------------------------
reads the tracks, generates the lines
	assumes:
		db set and openend
	arguments:
		none
	sets:
		this->Tracks
		this->TrackId
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readTracks(){

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
	int GId;
	int offset;
	bool goOn = true;
	bool validtrack = false;
	int TrackNo = 0;
	this->nTracks3 = 0;

	while (goOn)
	//for (int TrackNo = 0; TrackNo < this->nTracks3; TrackNo++) //
	{
		sql_query = "SELECT * FROM tracks WHERE id=" + Int2Str(TrackNo) + " ORDER BY snap_id";
		sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

		if (sql_error != SQLITE_OK)
		{
			vtkErrorMacro("Error with sql query! Error: " << sql_error + sql_query);
			return 0;
		}
		
		vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
		
		while (sqlite3_step(res) == SQLITE_ROW)
		{
			validtrack = true;
			snap_id = sqlite3_column_int(res, 1);
			GId = sqlite3_column_int(res, 2);

			//generate the lines
			//if (GId<1) {continue;} //check for empty fields TODO ADAPT THIS
			

			offset = this->SnapInfo3.at(snap_id).Offset + GId;

			// check for exact 0 coordinates, then it's no real point
			// just a temp solution
			if (this->Position->GetData()->GetComponent(offset,0) == 0 && 
				this->Position->GetData()->GetComponent(offset,1) == 0 &&
				this->Position->GetData()->GetComponent(offset,2) == 0)
			{
				continue;
			}

			nextLine->GetPointIds()->InsertNextId(offset);

			//fill the TrackId vector
			this->TrackId->InsertTuple1(offset,TrackNo);
		}
		
		if (validtrack)
		{
			this->Tracks->InsertNextCell(nextLine);
			TrackNo++;
			validtrack = false;
		}
		else
		{
			if (this->nTracks3<TrackNo) {this->nTracks3 = TrackNo;}
			goOn = false;
		}
	}
	vtkErrorMacro("track count: " << this->nTracks3);
	return 1;

}
/*----------------------------------------------------------------------------
Reads in addidtional infos for the snapshots (redshift, time, npart)
	assumes:
		opened database, db set
	sets:
		this->snapinfo	information about snapshots
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::readSnapshotInfo()
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
	int i = 0;

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		i++;
		snap_id = sqlite3_column_int(res, 0);
		if (snap_id<1){
			snap_id = i;
		}
		snap_id = sqlite3_column_int(res, 0);
		//this->SnapInfo3.at(snap_id-1).redshift = sqlite3_column_double(res, 2);
		//this->SnapInfo3.at(snap_id-1).time = sqlite3_column_double(res, 3);
		//this->SnapInfo3.at(snap_id-1).npart = sqlite3_column_int(res, 5);
	}
	return 1;
}

/*----------------------------------------------------------------------------
Generates the colorvector for the data


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

	//this->opacity = vtkSmartPointer<vtkUnsignedCharArray>::New();
	//this->opacity->SetName("Opacity");
	//this->opacity->SetNumberOfComponents(1);
	//this->opacity->SetNumberOfTuples(this->nParticles3);

	int snapid, trackid;
	unsigned char red, green, blue, alpha;

	//vtkSmartPointer<vtkLookupTable> LuT = vtkSmartPointer<vtkLookupTable>::New();
	//LuT->SetTableRange(1,this->nSnapshots);
	//double value [] = {0.75,0.25};
	//LuT->SetValueRange(value);
	//LuT->Build();

	for (int i = 0; i<this->nParticles3;i++)
	{
		snapid = this->SnapId->GetTuple1(i);
		trackid = this->TrackId->GetTuple1(i);

		red = 190*(float)(snapid) / (float)(this->nSnapshots-1);

		if (trackid == this->DisplaySelectedTrackNr)
		{
			green = 255;
		} else
		{
			green = 0;
		}
		
		if (snapid == this->DisplaySelectedSnapshotNr)
		{
			red = 255;
			green = 255;
		} else
		{
		}

		blue = 0;
		alpha = 0;

		this->colors->InsertTuple3(i,
				red,
				green,
				blue);
		//this->opacity->InsertTuple1(i,alpha);
		
		
		// Above code would work, i'm trying to do this with a lookup table instead
		//maybe this will be easier.
		// LuT works so far, but next problem: how to get the colors displayed?
		// i probably need an actor or something (for opicaty and such...)
		
		//double rgb [4] = {0,0,0,0};
		//LuT->GetColor(snapid,rgb);
		//this->colors->InsertTuple(i,rgb);
	}

	return 1;
}


/*----------------------------------------------------------------------------
Selects a subset of points to display


	assumes:
		
	sets:
		
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::selectPoints()
{
std::vector<int> idlist;

	if(*dataInfo.nSnapshots >-1)
	{
		for (int i =0; i<*dataInfo.nSnapshots; i++){
			for (int j = 0 ; j<*dataInfo.nParticles;j++){
				if(this->SnapId->GetTuple1(j)==dataInfo.selectedSnapshots->GetTuple1(i)){
					idlist.push_back(j);
				}
			}
		}
	}
	return 1;
}

/*----------------------------------------------------------------------------
Does some calculations

	assumes:
		
	sets:
		
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
int vtkSQLiteReader::doCalculations(){

	int startId;
	int endId;
	
	std::list<int> keeptracks;

	double parameter = *Gui.CalculationImpactParameter;

	for (int snapId = 0; snapId<*dataInfo.nSnapshots;snapId++)
	{
		startId = SnapInfo3.at(snapId).Offset;
		endId = SnapInfo3.at(snapId).Offset+SnapInfo3.at(snapId).lenght;

		for (int gid1 = startId; gid1<endId; gid1++)
		{
			for (int gid2 = startId; gid2<endId; gid2++)
			{
				if (gid1 == gid2) {continue;}
				if (distance(gid1,gid2)<parameter){
					vtkErrorMacro("   dist: "<<distance(gid1,gid2));
					keeptracks.push_back(gid1);
					keeptracks.push_back(gid2);
				}
			}
		}
	}
	keeptracks.sort();
	keeptracks.unique();
	vtkErrorMacro("gefunden: anzahl pkt: "<<keeptracks.size());

	//insert the gids to the dataInfo, selected particles.

	return 1;
}

/*----------------------------------------------------------------------------
calculates the distance between the points with gid1 and 2

	assumes:
		
	sets:
		
	arguments:
		ids of two points
	returns:
		the distance
*/
double vtkSQLiteReader::distance(int gid1, int gid2){

	double x1[3];
	double x2[3];

	this->Position->GetPoint(gid1, x1);
	this->Position->GetPoint(gid2, x2);

	double r1 = (x1[0]-x2[0])*(x1[0]-x2[0]);
	double r2 = (x1[1]-x2[1])*(x1[1]-x2[1]);
	double r3 = (x1[2]-x2[2])*(x1[2]-x2[2]);
	double tmp = sqrt( r1 + r2 +r3);
	
	return tmp;
}