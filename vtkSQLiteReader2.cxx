/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader2.cxx,v $

  
	todo:
	- debug for wrong input of nsnapshot, resp. input from display snapshot nr is > nsnapshot
	- check for tracksdisplay
	- maybe more computationonal stuff
	- option for estimation good tolerance parameter
	bugfixing:
	- lost color bug
	- check resetting of variables, whats needed, what not..  
  
=========================================================================*/
#include "vtkTimerLog.h"
#include "vtkSQLiteReader2.h"
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
#include "vtkKdTree.h"
#include "vtkTable.h"
#include "vtkUnicodeString.h"

#include <vector>
#include <sstream>
#include <string>


vtkCxxRevisionMacro(vtkSQLiteReader2, "$Revision: 1.0.1 $");
vtkStandardNewMacro(vtkSQLiteReader2);


//----------------------------------------------------------------------------
// Constructor
vtkSQLiteReader2::vtkSQLiteReader2()
{
	this->FileName          = 0;		// init filename
	this->SetNumberOfInputPorts(0);   // set no of input files (0 is just fine)
	this->SetNumberOfOutputPorts(3);

	this->db	= NULL;

	
	AllData = vtkSmartPointer<vtkPolyData>::New();
	SelectedData = vtkSmartPointer<vtkPolyData>::New();
	EmptyData = vtkSmartPointer<vtkPolyData>::New();

	TrackData = vtkSmartPointer<vtkPolyData>::New();
	SnapshotData = vtkSmartPointer<vtkPolyData>::New();


	//collect free gui variables into struct, and init
	Gui.DisplayColliding = &this->DisplayColliding;
	Gui.DisplayMerging = &this->DisplayMerging;
	Gui.DisplayInverted = &this->DisplayInverted;

	Gui.LowerLimit = &this->LowerLimit;
	Gui.UpperLimit = &this->UpperLimit;

	*Gui.DisplayColliding = false;
	*Gui.DisplayMerging = false;
	*Gui.DisplayInverted = false;

	*Gui.LowerLimit = 0.001;
	*Gui.UpperLimit = 0.0;


	//init

	// init dataInfo (maybee just call this->reset!??)
	dataInfo.InitComplete = false;
	dataInfo.dataIsRead = false;
	dataInfo.nPoints = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.nTracks = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.Hubble = 0;

	//init data stuff

	/*
	for (int i = 0; i<3; i++)
	{
		Data * actData;
		if (i==0){actData = &this->allData;}
		if (i==1){actData = &this->selectedData;}
		if (i==2){actData = &this->emptyData;}

		actData->Position = vtkSmartPointer<vtkPoints>::New();
		actData->Velocity = vtkSmartPointer<vtkFloatArray>::New();
		actData->Cells = vtkSmartPointer<vtkCellArray>::New();
		actData->Tracks = vtkSmartPointer<vtkCellArray>::New();
		actData->TrackId = vtkSmartPointer<vtkIntArray>::New();
		actData->GId = vtkSmartPointer<vtkIntArray>::New();
		actData->SnapId = vtkSmartPointer<vtkIntArray>::New();
		actData->Mvir = vtkSmartPointer<vtkFloatArray>::New();
		actData->Rvir = vtkSmartPointer<vtkFloatArray>::New();
		actData->Vmax = vtkSmartPointer<vtkFloatArray>::New();
		actData->Rmax = vtkSmartPointer<vtkFloatArray>::New();
		actData->Redshift = vtkSmartPointer<vtkFloatArray>::New();
		actData->Time = vtkSmartPointer<vtkFloatArray>::New();
		actData->Cv = vtkSmartPointer<vtkFloatArray>::New();
		actData->CollisionTypePoint = vtkSmartPointer<vtkUnsignedCharArray>::New();
		actData->CollisionTypeTrack = vtkSmartPointer<vtkUnsignedCharArray>::New();
		actData->CvAverage = vtkSmartPointer<vtkFloatArray>::New();
		actData->RGc = vtkSmartPointer<vtkFloatArray>::New();

		actData->Velocity->SetNumberOfComponents(3);
		actData->TrackId->SetNumberOfComponents(1);
		actData->GId->SetNumberOfComponents(1);
		actData->SnapId->SetNumberOfComponents(1);
		actData->Mvir->SetNumberOfComponents(1);
		actData->Rvir->SetNumberOfComponents(1);
		actData->Vmax->SetNumberOfComponents(1);
		actData->Rmax->SetNumberOfComponents(1);
		actData->Redshift->SetNumberOfComponents(1);
		actData->Time->SetNumberOfComponents(1);
		actData->Cv->SetNumberOfComponents(1);
		actData->CollisionTypePoint->SetNumberOfComponents(1);
		actData->CollisionTypeTrack->SetNumberOfComponents(1);
		actData->CvAverage->SetNumberOfComponents(1);
		actData->RGc->SetNumberOfComponents(1);

		actData->Velocity->SetName("Velocity");
		actData->TrackId->SetName("TrackId");
		actData->GId->SetName("GId");
		actData->SnapId->SetName("SnapId");
		actData->Mvir->SetName("Mvir");
		actData->Rvir->SetName("Rvir");
		actData->Vmax->SetName("Vmax");
		actData->Rmax->SetName("Rmax");
		actData->Redshift->SetName("Redshift");
		actData->Time->SetName("Time");
		actData->Cv->SetName("Cv");
		actData->CollisionTypePoint->SetName("CollisionType of Point");
		actData->CollisionTypeTrack->SetName("CollisionType of Track");
		actData->CvAverage->SetName("Cv (SnapshotAverage)");
		actData->RGc->SetName("Distance from Galactic Center");
	}*/

	this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
	this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
	this->collisionCalc.CalcIsDone = false;

	//init calc struct
	/*
	this->calcInfo.calcDone = false;
	this->calcInfo.nCollisions = -1;
	this->calcInfo.tolerance = 0;
*/
	// inti calcesttol2
	/*this->calcEstTol2 = vtkSmartPointer<vtkFloatArray>::New();
	this->calcEstTol2->SetNumberOfComponents(4);
	this->calcEstTol2->SetName("CalculationOfTolerance");*/
}

//----------------------------------------------------------------------------
// Deconstructor
vtkSQLiteReader2::~vtkSQLiteReader2()
{

	sqlite3_close(this->db);
}

//----------------------------------------------------------------------------
// Print self
void vtkSQLiteReader2::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	//TODO print out more stuff...
}


//----------------------------------------------------------------------------
// Output port
int vtkSQLiteReader2::FillOutputPortInformation(int port,	vtkInformation* info)
{
	if(port == 0)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}
	
	else if(port == 1)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}
	
	else if(port == 2)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}

	return 1;
}

/*----------------------------------------------------------------------------
reads some head data in, is called directly when loading the file (database), before apply is clicked (presents the infos for the window)
this currently does the following:
- checks the database format (right tables, right columns, important data not empty)
- checks what data is available and allocs arrays
*/
int vtkSQLiteReader2::RequestInformation(
		vtkInformation* vtkNotUsed(request),
		vtkInformationVector** vtkNotUsed(inputVector),
		vtkInformationVector* outputVector)
{
	/* Stuff for doing it in parallel, leave it for the moment...
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
	*/
	if(!this->dataInfo.InitComplete)
	{
		vtkSmartPointer<vtkTimerLog> inittimer = vtkSmartPointer<vtkTimerLog>::New();
		inittimer->StartTimer();

		// opening database
		if(!this->openDB(this->FileName)){return 0;}

		//init data
		vtkSmartPointer<vtkPoints> points;
		points = vtkSmartPointer<vtkPoints>::New();
		this->AllData->SetPoints(points);
		points = vtkSmartPointer<vtkPoints>::New();
		this->TrackData->SetPoints(points);
		points = vtkSmartPointer<vtkPoints>::New();
		this->SnapshotData->SetPoints(points);


		// add hepler arrays
		this->dataInfo.pSnapIdArray = this->CreateIdTypeArray(this->AllData,"SnapId");
		this->dataInfo.pGIdArray = this->CreateIdTypeArray(this->AllData,"GId");
		this->dataInfo.pTrackIdArray = this->CreateIdTypeArray(this->AllData,"TrackId");
		this->dataInfo.dataArrayOffset = 3;

		// read in database header
		if(!this->ReadHeader()){return 0;}

		//add aditional arrays
		this->CreateArray(this->AllData, "Cv"); //halo concentration
		this->CreateArray(this->AllData, "Redshift");
		this->CreateArray(this->AllData, "a"); //scale factor = 1/(1+redshift)


		this->SelectedData->DeepCopy(this->AllData);
		this->EmptyData->DeepCopy(this->AllData);
		
		this->dataInfo.InitComplete = true;

		//init timers
		this->timing.tot = 0;
		this->timing.read = 0;
		this->timing.calc = 0;
		this->timing.refresh = 0;
		this->timing.generateSelect = 0;
		this->timing.executeSelect = 0;


		inittimer->StopTimer();
		//vtkErrorMacro("Initialisation took: " << inittimer->GetElapsedTime() << " s");
		this->timing.init = inittimer->GetElapsedTime();

	}
	return 1;
}

/*----------------------------------------------------------------------------
is called after clicking apply, reads the actual, selected data
*/
int vtkSQLiteReader2::RequestData(vtkInformation*,
								 vtkInformationVector**,
								 vtkInformationVector* outputVector)
{
	vtkSmartPointer<vtkTimerLog> tottimer = vtkSmartPointer<vtkTimerLog>::New();
	tottimer->StartTimer();

	vtkPolyData * out = vtkPolyData::GetData(outputVector->GetInformationObject(0));
	vtkPolyData * snapshotout = vtkPolyData::GetData(outputVector->GetInformationObject(1));
	vtkPolyData * trackout = vtkPolyData::GetData(outputVector->GetInformationObject(2));

	if (!this->dataInfo.dataIsRead)
	// only read if it's not been read before
	{
		vtkSmartPointer<vtkTimerLog> readtimer = vtkSmartPointer<vtkTimerLog>::New();
		readtimer->StartTimer();

		readSnapshotInfo();
		readSnapshots();
		readTracks();
		calculatePointData();

		this->dataInfo.dataIsRead = true;

		readtimer->StopTimer();
		//vtkErrorMacro(" reading data took: " << readtimer->GetElapsedTime() << " s");
		this->timing.read = readtimer->GetElapsedTime();
	}

	if (this->collisionCalc.MergingTolerance != *this->Gui.LowerLimit ||
		this->collisionCalc.CollisionTolerance != *this->Gui.UpperLimit)
	{
		this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
		this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
		this->collisionCalc.CalcIsDone = false;
	}

	if((!this->collisionCalc.CalcIsDone) && (*Gui.DisplayColliding || *Gui.DisplayMerging))
	{
		vtkSmartPointer<vtkTimerLog> calctimer = vtkSmartPointer<vtkTimerLog>::New();
		calctimer->StartTimer();

		this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
		this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
		findCollisions(&this->collisionCalc);
		this->collisionCalc.CalcIsDone = true;

		calctimer->StopTimer();
		//vtkErrorMacro(" calculations took: " << calctimer->GetElapsedTime() << " s");
		this->timing.calc = calctimer->GetElapsedTime();
	}

	// refresh the ouput
	vtkSmartPointer<vtkTimerLog> refreshtimer = vtkSmartPointer<vtkTimerLog>::New();
	refreshtimer->StartTimer();

	if(!*Gui.DisplayColliding && !*Gui.DisplayMerging && !*Gui.DisplayInverted)
	//nothing selected = display all points
	{
		out->DeepCopy(this->AllData);
	}
	else if (!*Gui.DisplayColliding && !*Gui.DisplayMerging && *Gui.DisplayInverted)
	//dont display anything..
	{
		out->DeepCopy(this->EmptyData);
	}
	else
	{
		vtkSmartPointer<vtkTimerLog> genStimer = vtkSmartPointer<vtkTimerLog>::New();
		genStimer->StartTimer();

		generateSelection(&this->collisionCalc, &this->selection);

		genStimer->StopTimer();
		//vtkErrorMacro(" generating selection: " << genStimer->GetElapsedTime() << " s");
		this->timing.generateSelect = genStimer->GetElapsedTime();


		vtkSmartPointer<vtkTimerLog> seldtimer = vtkSmartPointer<vtkTimerLog>::New();
		seldtimer->StartTimer();

		generateSelectedData(this->SelectedData, &this->selection);

		seldtimer->StopTimer();
		//vtkErrorMacro(" generating selected data: " << seldtimer->GetElapsedTime() << " s");
		this->timing.executeSelect = seldtimer->GetElapsedTime();


		out->DeepCopy(this->SelectedData);
	}
	
	refreshtimer->StopTimer();
	//vtkErrorMacro(" refreshing took: " << refreshtimer->GetElapsedTime() << " s");
	this->timing.refresh = refreshtimer->GetElapsedTime();

	// set additional data output
	snapshotout->DeepCopy(this->SnapshotData);
	trackout->DeepCopy(this->TrackData);

	tottimer->StopTimer();
	//vtkErrorMacro("Total elapsed time: " << tottimer->GetElapsedTime() << " s");
	this->timing.main = tottimer->GetElapsedTime();
	this->timing.tot = this->timing.main + this->timing.init;


		vtkstd::stringstream ss;
		ss<<"\n\nSQLite2Reader was successful!\n\n";
		ss<<"   Headerdata looked like\n";
		ss<<"      Points                  : ";
		ss<<this->initData.nPoints<<"\n";
		ss<<"      Tracks                  : ";
		ss<<this->initData.nTracks<<"\n";
		ss<<"      Snapshots               : ";
		ss<<this->initData.nSnapshots<<"\n\n";
		ss<<"   actually read in\n";
		ss<<"      Points                  : ";
		ss<<this->dataInfo.nPoints<<"\n";
		ss<<"      Tracks                  : ";
		ss<<this->dataInfo.nTracks<<"\n";
		ss<<"      Snapshots               : ";
		ss<<this->dataInfo.nSnapshots<<"\n\n";
	if(this->collisionCalc.CalcIsDone){
		ss<<"   calculated\n";
		ss<<"      colliding tracks        : ";
		ss<<this->collisionCalc.CollisionTrackIds->GetNumberOfIds()<<"\n";
		ss<<"      merging tracks          : ";
		ss<<this->collisionCalc.MergingTrackIds->GetNumberOfIds()<<"\n\n";
		ss<<"   selected\n";
		ss<<"      selected points         : ";
		ss<<this->selection.Points->GetNumberOfIds()<<"\n\n";
	}
		ss<<"   Timings\n";
		ss<<"        Init                  : ";
		ss<<this->timing.init<<" s\n";
		ss<<"        Runtime mainmodule    : ";
		ss<<this->timing.main<<" s\n\n";
		ss<<"          Reading             : ";
		ss<<this->timing.read<<" s\n";
		ss<<"          Calculations        : ";
		ss<<this->timing.calc<<" s\n";
		ss<<"          generating Selection: ";
		ss<<this->timing.generateSelect<<" s\n";
		ss<<"          executing Selection : ";
		ss<<this->timing.executeSelect<<" s\n";
		ss<<"          refreshing Output   : ";
		ss<<this->timing.refresh<<" s\n";
		ss<<"      TOTAL TIME TAKEN        : ";
		ss<<this->timing.tot<<" s\n";

	vtkErrorMacro(<<ss.str());


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
int vtkSQLiteReader2::openDB(char* filename)
{
	if(!filename) //check for existing filename
	{
		vtkErrorMacro("A FileName must be specified.");
		return 0;
	}
	
	if (sqlite3_open(this->FileName, &db) != SQLITE_OK) // open db, returns SQLITE_OK==0 if successful
	{
		vtkErrorMacro("Can't open database: " + *sqlite3_errmsg(db));
		return 0;
	}

	//vtkErrorMacro("opened successfully: " << this->FileName)
	return 1;
}


/*----------------------------------------------------------------------------
Converts a Integer to a string
*/
vtkStdString vtkSQLiteReader2::Int2Str(int number)
{
	std::stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

/*----------------------------------------------------------------------------
Reads the header data from db, sets up data structres, checks datastructres
	assumes:
		opened database, db set
	arguments:
		vtkInformationVector Output information
	returns:
		int	errorcode (1 =ok)
	sets:
		numSnap		Number of snapshots

*/
int vtkSQLiteReader2::ReadHeader()
{
	// set up local vars
	int counter;
	vtkstd::string name;

	// set up the sql stuff
	vtkStdString	sql_query;	// a querry
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query
	
	//------------------------------------------------------
	// read the header table
	sql_query = "SELECT * FROM header";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	if (sqlite3_step(res) == SQLITE_ROW) //there should be only one row..
	{
		/* assume header has always the same structre:
		0 name
		1 format
		2 version
		3 revision
		4 type
		5 Omega0
		6 OmegaLambda
		7 Hubble
		8 ndims
		9 nsnapshots
		10 nsteps

		read (some of) the data in, add more if needed (don't forget to add them in the header, DataInformation struct)
		*/

		// 5 Omega0
		name = sqlite3_column_name(res,5);
		if (name.compare("Omega0")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 5);
			if (value>0){this->dataInfo.Omega0 = value;}
			else {this->dataInfo.Omega0 = value;} // you could define a standart Omega0 to use here
		}

		// 6 OmegaLambda
		name = sqlite3_column_name(res,6);
		if (name.compare("OmegaLambda")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 6);
			if (value>0){this->dataInfo.OmegaLambda = value;}
			else {this->dataInfo.OmegaLambda = value;} // you could define a standart OmegaLambda to use here
		}

		// 7 Hubble
		name = sqlite3_column_name(res,7);
		if (name.compare("Hubble")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 7);
			if (value>0){this->dataInfo.Hubble = value;}
			else {this->dataInfo.Hubble = 70.4;}
		}

		// 9 nsnapshots
		/* not used, instead nSnapshots is gotten by reading snpinfo table
		int num = sqlite3_column_int(res, 9);
		if (num>0){this->dataInfo.nSnapshots = num;}
		else {this->dataInfo.nSnapshots = -1;}
		*/
	}
	else
	{
		vtkErrorMacro("error reading header table");
		return 0;
	}

	//------------------------------------------------------
	// read snapinfo
	sql_query = "SELECT * FROM snapinfo";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	counter = 0;
	if (sqlite3_step(res) == SQLITE_ROW)
	{
		// read in snapinfo table header
		this->dataInfo.SnapinfoDataColumns.clear();
		this->dataInfo.SnapinfoSnapidColumn = -1;

		for (int i = 0; i<sqlite3_data_count(res); i++)
		{
			name = sqlite3_column_name(res,i);
			const unsigned char * ptr = sqlite3_column_text(res, i);
			if (name.compare("snap_id")==0 && ptr != NULL)
			{
				this->dataInfo.SnapinfoSnapidColumn = i;
			}
			else if (ptr != NULL)
			{
				this->dataInfo.SnapinfoDataColumns.push_back(i);
				this->CreateArray(this->SnapshotData, name.c_str());
			}
		}
		
		//count the num of snapshots
		do {++counter;}	while (sqlite3_step(res) == SQLITE_ROW);

		if(counter!=0)
		{
			this->dataInfo.nSnapshots = counter;
			//this->InitAllArrays(this->SnapshotData,counter);
		}
		else
		{
			vtkErrorMacro("wrong db structre (table 'snapinfo': no data)");
		}
	}
	// something went wrong
	else
	{
		vtkErrorMacro("wrong db structre (table 'snapinfo': querry error)")
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);

	//------------------------------------------------------
	// read tracks
	sql_query = "SELECT * FROM tracks";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	int id, gid;
	int max=0, min=0, gidmax=0;
	if (sqlite3_step(res) == SQLITE_ROW)
	{
		vtkstd::string id1 = sqlite3_column_name(res,0);
		vtkstd::string id2 = sqlite3_column_name(res,1);
		vtkstd::string id3 = sqlite3_column_name(res,2);
		if (id1.compare("id")==0 &&
			id2.compare("snap_id")==0 &&
			id3.compare("gid")==0)
		{
			// if necessairy, implement here automatic column detection...
			// for now just assume id in col 0, snap_id in col 1, gid in col 2
			this->dataInfo.TracksTrackidColumn = 0;
			this->dataInfo.TracksSnapidColumn = 1;
			this->dataInfo.TracksGidColumn = 2;

			do
			{
				id = sqlite3_column_int(res,0);
				if(id>max){max=id;}
				if(id<min){min=id;}
				gid = sqlite3_column_int(res,2);
				if(gid>gidmax){gidmax=gid;}

			} while (sqlite3_step(res) == SQLITE_ROW);
			this->dataInfo.gidmax = gidmax;
		}
		// something wrong with db structre
		else
		{
			vtkErrorMacro("wrong db structre (table 'tracks': no id / snap_id / gid columns)");
			return 0;
		}
	}
	
	if (max!=0)
	{
		this->dataInfo.nTracks = max-min+1;
	}
	// something went wrong
	else
	{
		vtkErrorMacro("wrong db structre (table 'tracks': no data)");
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);

	//------------------------------------------------------
	// count the points
	sql_query = "SELECT * FROM stat";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	counter = 0;

	if (sqlite3_step(res) == SQLITE_ROW)
	{
		this->dataInfo.StatCordinateColumns.clear();
		this->dataInfo.StatCordinateColumns.resize(3);
		for (int i = 0; i<sqlite3_data_count(res); i++)
		{
			name = sqlite3_column_name(res,i);
			const unsigned char * ptr = sqlite3_column_text(res, i);

			if(name.compare("Xc")==0) {this->dataInfo.StatCordinateColumns.at(0) = i;}
			else if(name.compare("Yc")==0) {this->dataInfo.StatCordinateColumns.at(1) = i;}
			else if(name.compare("Zc")==0) {this->dataInfo.StatCordinateColumns.at(2) = i;}
			else if(name.compare("snap_id")==0) {this->dataInfo.StatSnapidColumn = i;}
			else if(name.compare("gid")==0) {this->dataInfo.StatGidColumn = i;}
			/* not used, instead just get simple colums from velo vector...
			else if(name.compare("VXc")==0) {this->dataInfo.StatCordinateColumns.at(0) = i;}
			else if(name.compare("VYc")==0) {this->dataInfo.StatCordinateColumns.at(1) = i;}
			else if(name.compare("VZc")==0) {this->dataInfo.StatCordinateColumns.at(2) = i;}
			*/
			else if (ptr != NULL)
			{
				this->dataInfo.StatDataColumns.push_back(i);
				this->CreateArray(this->AllData, name.c_str());
			}
		}

		do{counter++;}
		while (sqlite3_step(res) == SQLITE_ROW);

		this->dataInfo.nPoints = counter;
	} 
	else
	{
		vtkErrorMacro("wrong db structre (table 'stat': no data)");
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);


	//------------------------------------------------------
	// cleaning up
	/*vtkErrorMacro(" ready to read: " <<
		this->dataInfo.nPoints << " Points in "<<
		this->dataInfo.nSnapshots << " Snapshots, on "<<
		this->dataInfo.nTracks << " tracks.");*/
	this->initData.nPoints = this->dataInfo.nPoints;
	this->initData.nSnapshots = this->dataInfo.nSnapshots;
	this->initData.nTracks = this->dataInfo.nTracks;
	return 1;
}

/*----------------------------------------------------------------------------
reads the snapshots in (reads all the points, and the according data, generates the cells)
	assumes:
		db set and openend
	arguments:
		none
	sets:

	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader2::readSnapshots()
{
	// sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM stat";// ORDER BY snap_id";
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	//shotcuts
	int					xC = this->dataInfo.StatCordinateColumns.at(0);
	int					yC = this->dataInfo.StatCordinateColumns.at(1);
	int					zC = this->dataInfo.StatCordinateColumns.at(2);
	int					SnapidC = this->dataInfo.StatSnapidColumn;
	int					GIdC = this->dataInfo.StatGidColumn;
	vtkstd::vector<int>	*dC = &this->dataInfo.StatDataColumns;

	vtkIdTypeArray		*SnapidA = this->dataInfo.pSnapIdArray;// pointer to snapid array
	vtkIdTypeArray		*GIdA = this->dataInfo.pGIdArray;// pointer to GId array
	int					offsetA = this->dataInfo.dataArrayOffset; // wich is the first of the dataarrays.

	vtkstd::vector<vtkSmartPointer<vtkIdList>> * snapI = &this->SnapInfo2;
	vtkPointData		*pData = this->AllData->GetPointData();
	vtkPoints			*pnts = this->AllData->GetPoints();

	int					gidmax = this->dataInfo.gidmax;

	// prepare variables
	snapI->resize(this->dataInfo.nSnapshots);
	for (int i = 0; i<this->dataInfo.nSnapshots; ++i)
	{
		//snapI->at(i).PointIds.resize(gidmax+1,-1);
		snapI->at(i) = vtkSmartPointer<vtkIdList>::New();
		snapI->at(i)->SetNumberOfIds(gidmax+1);
		for (vtkIdType j = 0; j<gidmax+1;j++){snapI->at(i)->SetId(j,-1);}
	}
	pnts->SetNumberOfPoints(this->dataInfo.nPoints);
	this->InitAllArrays(this->AllData,this->dataInfo.nPoints);


	// loop variables
	int counter = 0;
	int i, snapid, gid;
	double x,y,z;


	while (sqlite3_step(res) == SQLITE_ROW)
	{
		x = sqlite3_column_double(res, xC);
		y = sqlite3_column_double(res, yC);
		z = sqlite3_column_double(res, zC);

		if (x==0.0 && y==0.0 && z==0.0){continue;}

		pnts->InsertPoint(counter,x,y,z);

		snapid = sqlite3_column_int(res, SnapidC);
		gid = sqlite3_column_int(res, GIdC);

		SnapidA->InsertTuple1(counter,snapid);
		GIdA->InsertTuple1(counter,gid);

		//snapI->at(snapid).PointIds.at(gid) = counter;
		snapI->at(snapid)->SetId(gid, counter);

		for (i=0;i<dC->size();++i)
		{
			pData->GetArray(i+offsetA)->InsertTuple1(counter,sqlite3_column_double(res, dC->at(i)));
		}
		++counter;
	}

	pnts->SetNumberOfPoints(counter);
	for (i=0; i<pData->GetNumberOfArrays(); ++i)
	{
		pData->GetArray(i)->SetNumberOfTuples(counter);
	}
	this->dataInfo.nPoints = counter;

	// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = counter;
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType *pverts = verts->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		pverts[i*2]   = 1;
		pverts[i*2+1] = i;
    }
	this->AllData->SetVerts(verts);

	/*
	int count = 0;
	int snapcount = 0;
	int snap_id=0;
	int old_snap_id = 0;
	int GId = 0;
	int offset = 0;
	double x,y,z;
	double cv, vmax, rmax;
	std::vector<int> pid;
	SnapshotInfo *actSnap;
	double redshift = this->SnapInfo.at(snap_id).redshift;
	double hubble = 70.4;//this->dataInfo.hubble;

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		old_snap_id = snap_id;
		snap_id = sqlite3_column_double(res, 0);

		if (snap_id != old_snap_id)
		{
			actSnap = &this->SnapInfo.at(old_snap_id);
			actSnap->Offset = offset;
			actSnap->lenght = snapcount;
			actSnap->PointId = pid;

			pid.clear();
			snapcount = 0;
			offset = count;
			redshift = this->SnapInfo.at(snap_id).redshift;
		}

		GId = sqlite3_column_double(res, 1);

		//filter out any points with coordinates = 0
		x = sqlite3_column_double(res, 4);
		y = sqlite3_column_double(res, 5);
		z = sqlite3_column_double(res, 6);

		if(x==0&&y==0&&z==0)
		{
			pid.push_back(-1);
			continue;
		}

		this->allData.Position->InsertNextPoint(x,y,z);
		this->allData.Velocity->InsertNextTuple3(
			sqlite3_column_double(res, 7),
			sqlite3_column_double(res, 8),
			sqlite3_column_double(res, 9));
		this->allData.GId->InsertNextTuple1(GId);
		this->allData.SnapId->InsertNextTuple1(snap_id);

		this->allData.Mvir->InsertNextTuple1(sqlite3_column_double(res, 10));
		this->allData.Rvir->InsertNextTuple1(sqlite3_column_double(res, 11));
		vmax = sqlite3_column_double(res, 12);
		this->allData.Vmax->InsertNextTuple1(vmax);
		rmax = sqlite3_column_double(res, 13);
		this->allData.Rmax->InsertNextTuple1(rmax);
		this->allData.Redshift->InsertNextTuple1(redshift);
		this->allData.Time->InsertNextTuple1(1.0/(1.0+redshift));

		//compute additional stuff
		cv = 2.0 * (vmax / (hubble * rmax)) * (vmax / (hubble * rmax));
		this->allData.Cv->InsertNextTuple1(cv);

		pid.push_back(snapcount);

		snapcount++;
		count++;
	}

	actSnap = &this->SnapInfo.at(snap_id);
	actSnap->Offset = offset;
	actSnap->lenght = snapcount;
	actSnap->PointId = pid;

	//this->nTracks3 = (int)10;
	//vtkErrorMacro("No of tracks: " << this->nTracks3);

	this->dataInfo.nSnapshots = (int)this->SnapInfo.size();
	vtkErrorMacro("No of snapshots: " << this->dataInfo.nSnapshots);

	this->dataInfo.nPoints = count;
	vtkErrorMacro("No of particles: " << this->dataInfo.nPoints);

	// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = this->allData.Position->GetNumberOfPoints();
	vtkIdType *cells = this->allData.Cells->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		cells[i*2]   = 1;
		cells[i*2+1] = i;
    }

	// init other field data
	this->allData.CollisionTypePoint->SetNumberOfTuples(this->dataInfo.nPoints);
	this->allData.CollisionTypePoint->FillComponent(0,0);

	this->allData.CollisionTypeTrack->SetNumberOfTuples(this->dataInfo.nPoints);
	this->allData.CollisionTypeTrack->FillComponent(0,0);

	*/
	this->dataInfo.dataIsRead = true;
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
int vtkSQLiteReader2::readTracks(){

	//set up sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	//prepare sql querry
	sql_query = "SELECT * FROM tracks";
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	//shotcuts
	int TrackidC = this->dataInfo.TracksTrackidColumn;
	int SnapidC = this->dataInfo.TracksSnapidColumn;
	int GidC = this->dataInfo.TracksGidColumn;

	vtkIdTypeArray * TracksA = this->dataInfo.pTrackIdArray;
	vtkstd::vector<vtkSmartPointer<vtkIdList>> * ti = &this->TrackInfo2;
	vtkstd::vector<vtkSmartPointer<vtkIdList>> * si = &this->SnapInfo2;

	//init
	vtkSmartPointer<vtkCellArray> Tracks = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyLine> nextLine;
	//vtkstd::vector<int> nPointsOnTrack;

	ti->clear();
	ti->resize(this->dataInfo.nTracks);
	//nPointsOnTrack.clear();
	//nPointsOnTrack.resize(this->dataInfo.nTracks,0);

	for (int i = 0; i<this->dataInfo.nTracks; ++i)
	{
		//ti->at(i).nPoints=0;
		//ti->at(i).PointIds.clear();
		//ti->at(i).PointIds.resize(this->dataInfo.nSnapshots,-1);
		ti->at(i) = vtkSmartPointer<vtkIdList>::New();
		ti->at(i)->SetNumberOfIds(this->dataInfo.nSnapshots);
		for (vtkIdType j = 0; j<this->dataInfo.nSnapshots;++j){ti->at(i)->SetId(j,-1);}
	}

	TracksA->SetNumberOfTuples(this->dataInfo.nPoints);
	TracksA->FillComponent(0,-1);

	// loop vars
	vtkIdType trackid, snapid, gid, newid;
	//vtkIdList * PointsOnTrack;
	//vtkIdList * Track;

	//read in the data
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		trackid = sqlite3_column_int(res, TrackidC);
		snapid = sqlite3_column_int(res, SnapidC);
		gid = sqlite3_column_int(res, GidC);

		//newid = si->at(snapid).PointIds.at(gid);
		newid = si->at(snapid)->GetId(gid);

		if(newid>-1)
		{
			//ti->at(trackid).PointIds.at(snapid) = newid;
			//ti->at(trackid).nPoints++;
			ti->at(trackid)->SetId(snapid,newid);
			//++nPointsOnTrack.at(trackid);

			TracksA->SetTuple1(newid,trackid);
		}
	}

	//delete the fields with no values (this is necessairy, because the data read in from tracks
	// isnt necessairy ordered.. otherwise this could be done much quicker by chaning the above loop
	// starting from an empty vector and just insert the values...
	for (int i = 0; i<ti->size(); ++i)
	{
		ti->at(i)->DeleteId(-1);
		
		//PointsOnTrack = &ti->at(i).PointIds;
		/*Track = &ti->at(i);

		

		

		vtkstd::vector<vtkIdType>::reverse_iterator rev_iter = PointsOnTrack->rbegin();
		// note: explanation of use of reverse iterators:
		// http://www.drdobbs.com/cpp/184401406;jsessionid=N3YFDK53XSCRPQE1GHPCKHWATMY32JVN
		// and sample code from here:
		// http://stackoverflow.com/questions/404258/how-do-i-erase-a-reverse-iterator-from-an-stl-data-structure

		while ( rev_iter != PointsOnTrack->rend())
		{
			// Find 3 and try to erase
			if (*rev_iter == -1)
			{
				//cout << "Erasing : " << *rev_iter;
				vtkstd::vector<vtkIdType>::iterator tempIter = PointsOnTrack->erase( --rev_iter.base());
				rev_iter = vtkstd::vector<vtkIdType>::reverse_iterator(tempIter);                        
			}
			else if (Track->nPoints == PointsOnTrack->size())
			{
				break;
			}
			else
			{
				++rev_iter;
			}
		}   */

		
		/* simple way with forward iterator

		vtkstd::vector<vtkIdType>::iterator iter = PointsOnTrack->begin();
		while( iter != PointsOnTrack->end() )
		{
			if (*iter == -1)
			{
				iter = PointsOnTrack->erase( iter );
			}
			else
			{
				++iter;
			}
		}
		*/
	}


	vtkPoints * tmpPoints = this->TrackData->GetPoints();
	tmpPoints->SetNumberOfPoints(this->dataInfo.nTracks);

	vtkSmartPointer<vtkIntArray> nPointsOnTrackDataArray = vtkSmartPointer<vtkIntArray>::New();
	nPointsOnTrackDataArray->SetNumberOfComponents(1);
	nPointsOnTrackDataArray->SetNumberOfTuples(this->dataInfo.nTracks);
	nPointsOnTrackDataArray->SetName("nPointsOnTrack");

	// you can add here more data if you want...

	for (int i = 0; i<ti->size(); ++i)
	{
		/*
		Track = ti->at(i);
		
		PointsOnTrack = &ti->at(i).PointIds;
		nextLine = vtkSmartPointer<vtkPolyLine>::New();

		for (int j=0; j<Track->nPoints;++j)
		{
			nextLine->GetPointIds()->InsertNextId(Track->PointIds.at(j));
		}
		Tracks->InsertNextCell(nextLine);
		*/
		Tracks->InsertNextCell(ti->at(i));

		//Fill Data Arrays here..
		tmpPoints->InsertPoint(i,0,0,i);
		//nPointsOnTrackDataArray->InsertTuple1(i,Track->nPoints);
		nPointsOnTrackDataArray->InsertTuple1(i,ti->at(i)->GetNumberOfIds());
	}

	this->AllData->SetLines(Tracks);

	this->TrackData->GetPointData()->AddArray(nPointsOnTrackDataArray);



	
/*
	this->allData.TrackId->SetNumberOfTuples(this->dataInfo.nPoints);

	//init everything to -1 (-1: halo belongs to no track)
	for (int i = 0; i<this->dataInfo.nPoints; i++)
	{
		this->allData.TrackId->InsertTuple1(i,-1);
	}

	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	int snap_id;
	int GId, pid;
	int offset;
	bool goOn = true;
	bool validtrack = false;
	int TrackNo = 0;
	this->dataInfo.nTracks = 0;


	while (goOn)
	{
		sql_query = "SELECT * FROM tracks WHERE id=" + Int2Str(TrackNo) + " ORDER BY snap_id";
		sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);

		if (sql_error != SQLITE_OK)
		{
			vtkErrorMacro("Error with sql query! Error: " << sql_error + sql_query);
			return 0;
		}
		
		vtkSmartPointer<vtkPolyLine> nextLine = vtkSmartPointer<vtkPolyLine>::New();
		Track tmpTrack; // for filling in the trackinfo struct
		tmpTrack.PointsIds.resize(this->dataInfo.nSnapshots, -1);
		tmpTrack.nPoints = 0;
		
		while (sqlite3_step(res) == SQLITE_ROW)
		{
			validtrack = true;
			snap_id = sqlite3_column_int(res, 1);
			GId = sqlite3_column_int(res, 2);

			//generate the lines
			//if (GId<1) {continue;} //check for empty fields TODO ADAPT THIS
			
			pid = this->SnapInfo.at(snap_id).PointId.at(GId);
			if(pid<0){continue;}

			offset = this->SnapInfo.at(snap_id).Offset + pid;

			// check for exact 0 coordinates, then it's no real point
			// just a temp solution
			// TODO fix this..
			//if (this->allData.Position->GetData()->GetComponent(offset,0) == 0 && 
			//	this->allData.Position->GetData()->GetComponent(offset,1) == 0 &&
			//	this->allData.Position->GetData()->GetComponent(offset,2) == 0)
			//{
			//	continue;
			//}

			nextLine->GetPointIds()->InsertNextId(offset);

			//fill the TrackId vector (associates points to tracks)
			this->allData.TrackId->InsertTuple1(offset,TrackNo);

			// fill in the data for trackinfo vector (for association track to point)
			tmpTrack.PointsIds.at(snap_id) = offset;
			tmpTrack.nPoints++;
		}
		
		if (validtrack)
		{
			this->allData.Tracks->InsertNextCell(nextLine);

			this->TracksInfo.push_back(tmpTrack);
			TrackNo++;
			validtrack = false;
		}
		else
		{
			if (this->dataInfo.nTracks < TrackNo) {this->dataInfo.nTracks = TrackNo;}
			goOn = false;
		}
	}
	
*/

	//vtkErrorMacro("track count: " << this->dataInfo.nTracks);
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

int vtkSQLiteReader2::readSnapshotInfo()
{
	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM snapinfo ORDER BY snap_id, redshift, time";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error!=SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}
	
	// init vars
	int snapId;
	int counter = 0;
	int SnapinfoSnapidColumn = this->dataInfo.SnapinfoSnapidColumn;
	bool SnapidExists = true;
	if (SnapinfoSnapidColumn == -1) {SnapidExists = false;}
	vtkPointData * pData = this->SnapshotData->GetPointData();

	this->InitAllArrays(this->SnapshotData,this->dataInfo.nSnapshots);

	// fill points
	vtkSmartPointer<vtkPoints> pnt = vtkSmartPointer<vtkPoints>::New();
	pnt->SetNumberOfPoints(this->dataInfo.nSnapshots);
	for (int i = 0; i<this->dataInfo.nSnapshots;i++)
	{
		pnt->InsertPoint(i,0,0,i);
	}
	this->SnapshotData->SetPoints(pnt);

	// main loop
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		if (SnapidExists){snapId = sqlite3_column_int(res, SnapinfoSnapidColumn);}
		else {snapId = counter;}

		for (int i = 0; i<this->dataInfo.SnapinfoDataColumns.size(); i++)
		{
			pData->GetArray(i)->InsertTuple1(snapId,sqlite3_column_double(res, this->dataInfo.SnapinfoDataColumns.at(i)));
		}
		++counter;
	}

	if (this->dataInfo.nSnapshots != counter) {vtkErrorMacro("reading snapshotInfo: something not good..");}
	//this->dataInfo.nSnapshots = counter;
	//this->ResizeAllArrays(this->SnapshotData,this->dataInfo.nSnapshots);

	return 1;
}


/*----------------------------------------------------------------------------
Generates additional data


	assumes:
		all data read in
	sets:
		this->color		color vector
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/
/*
int vtkSQLiteReader2::calculateAdditionalData()
{
	this->allData.CvAverage->SetNumberOfTuples(this->dataInfo.nPoints);
	this->allData.RGc->SetNumberOfTuples(this->dataInfo.nPoints);

	// calc average cv in one snapshot
	for (int i = 0; i<this->dataInfo.nSnapshots;i++)
	{
		int length = this->SnapInfo.at(i).lenght;
		int offset = this->SnapInfo.at(i).Offset;
		int count = 0;
		double sum = 0;
		double ave = 0;
		for (int j=0; j<length;j++)
		{
			sum =+ this->allData.Cv->GetTuple1(j+offset);
			count ++;
		}
		ave = sum / (double)count;

		this->SnapInfo.at(i).CvAverage = ave;
		for (int j = offset; j<offset+length;j++)
		{
			this->allData.CvAverage->InsertTuple1(j,ave);
		}


	}

	double dist;
	// calc dist from galactic center for every point
	for (int i = 0; i<this->dataInfo.nPoints; i++)
	{
		//dist = getDistanceToO(i);
		dist = sqrt(getDistance2(i,
			this->TracksInfo.at(this->dataInfo.nTracks-2).PointsIds.at(
			this->allData.SnapId->GetTuple1(i))));
		this->allData.RGc->InsertTuple1(i,dist);
	}
	return 1;
}

*/

/*----------------------------------------------------------------------------
tries to get an estimate for wich tolerance parameter yields to how many colliding tracks.
prints the endresult with vtkErrorMacro

	assumes:
		
	sets:
		int	nCollisions: Number of collisions
	arguments:
	returns:
		int errorcode
*/
int vtkSQLiteReader2::calcTolerance()
{
	/*
	// settings:
	double accFactor = 2; // to witch factor should the desired nCollision be reached
	int nTries = 5; //how many max tires to reach a goal
	double startFactor = 0.000001; // factor for the starting parameter

	bool done = false;
	double parameter;// = 0.0000001;
	int nCollisions;
	double acc;
	int counter;
	calcEstTol.clear();
	
	double minmaxx [2];
	this->allData.Position->GetData()->GetRange(minmaxx,0);
	double minmaxy [2];
	this->allData.Position->GetData()->GetRange(minmaxy,1);
	double minmaxz [2];
	this->allData.Position->GetData()->GetRange(minmaxz,2);

	double rangex = minmaxx[1]-minmaxx[0];
	double rangey = minmaxy[1]-minmaxy[0];
	double rangez = minmaxz[1]-minmaxz[0];

	double range = std::max(std::max(rangex,rangey),rangez);

	//estimate a good starting parameter based on the range of the points and their number
	parameter = range / pow(this->dataInfo.nTracks,0.333333) * startFactor;

	vtkErrorMacro("range: "<<range<<
		"; npoints: "<<this->dataInfo.nPoints<<
		"; starting parameter: "<<parameter);

	for (double goal = 1000; goal >= 10; goal = goal/10)
	//for (double goal = 10; goal >= 1; goal--)
	{
		counter = 0;
		vtkErrorMacro("goal: "<<goal);
		done = false;

		while(!done){
			counter++;
			//this->doCalculations(parameter,5);
			nCollisions = this->calcInfo.nCollisions;
			acc = (double)nCollisions / (double)goal;
			//vtkErrorMacro("  calc: goal: " << goal <<" parameter: " << parameter << " ncoll: "<< nCollisions<<" acc: "<<acc<<" counter: "<<counter );
			
			ResultOfEsimationOfTolerance tmpres;
			tmpres.goal = goal;
			tmpres.parameter = parameter;
			tmpres.nCollisions = nCollisions;
			calcEstTol.push_back(tmpres);

			this->calcEstTol2->InsertNextTuple4(parameter,nCollisions,goal,counter);

			if(counter == nTries)
			{
				done = true;
			}
			else if (acc == 0)
			{
				parameter = parameter * 10;
				continue;
			}
			else if (acc < (1.0 / accFactor))
			{
				parameter = parameter / (acc);
				continue;
			}
			else if (acc > accFactor)
			{
				parameter = parameter / (acc);
				continue;
			}
			vtkErrorMacro("RESULT: goal: " << goal <<" parameter: " << parameter << " ncoll: "<< nCollisions<<" counter: "<<counter );
			done = true;

		}
	}

	// filling up this array with nonsense data so it maches the others
	int n = this->calcEstTol2->GetNumberOfTuples(); //just making enough space
	this->calcEstTol2->Resize(this->dataInfo.nPoints);
	for (int i = n; i<this->dataInfo.nPoints;i++)
	{
		this->calcEstTol2->InsertTuple4(i,0,0,0,0);
	}
*/
	return 1;
}


/*----------------------------------------------------------------------------
calculates the distance (squared) between the points with gid1 and 2

	assumes:
		
	sets:
		
	arguments:
		ids of two points
	returns:
		the distance
*/
double vtkSQLiteReader2::getDistance2(int gid1, int gid2){

	double x1[3];
	double x2[3];

	this->allData.Position->GetPoint(gid1, x1);
	this->allData.Position->GetPoint(gid2, x2);

	return (x1[0]-x2[0])*(x1[0]-x2[0])+
		(x1[1]-x2[1])*(x1[1]-x2[1])+
		(x1[2]-x2[2])*(x1[2]-x2[2]);
}
/*----------------------------------------------------------------------------
calculates the distance  between a point and the origin

	assumes:
		
	sets:
		
	arguments:
		ids of two points
	returns:
		the distance
*/
/*double vtkSQLiteReader2::getDistanceToO(int gid1){

	double x1[3];
	double x2[] = {0.0,0.0,0.0};

	this->allData.Position->GetPoint(gid1, x1);

	return sqrt(x1[0]-x2[0])*(x1[0]-x2[0])+
		(x1[1]-x2[1])*(x1[1]-x2[1])+
		(x1[2]-x2[2])*(x1[2]-x2[2]);
}*/
/*----------------------------------------------------------------------------
generates the id map from alldata ids to only selected points ids
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
/*int vtkSQLiteReader2::generateIdMap()
{
	
	std::vector<int> * idMap1 = &this->dataInfo.idMap1;
	std::vector<int> * idMap2 = &this->dataInfo.idMap2;

	int counter = 0;
	int offset = 0;
	int length = 0;
	int pid = 0;
	int trackid = 0;

	idMap1->clear();
	idMap2->clear();
		
	idMap1->resize(this->dataInfo.nPoints,-1);
	idMap2->resize(this->dataInfo.nPoints,-1);

	
	for (int i = 0; i<this->dataInfo.nSelectedPoints; i++)
	{
		idMap1->at(this->dataInfo.selectedPoints.at(i)) = counter;
		idMap2->at(counter) = this->dataInfo.selectedPoints.at(i);
		counter++;
	}
	

	for (int i = 0; i<this->dataInfo.nSelectedSnapshots; i++)
	{
		if (i>=this->dataInfo.nSnapshots) {continue;} //overflow protection
		offset = this->SnapInfo.at(this->dataInfo.selectedSnapshots.at(i)).Offset;
		length = this->SnapInfo.at(this->dataInfo.selectedSnapshots.at(i)).lenght;

		for (int j = offset; j < offset+length; j++)
		{
			if(idMap1->at(j) < 0) //only if not already point assigned
			{
				idMap1->at(j) = counter;
				idMap2->at(counter) = j;
				counter++;
			}
		}
	}

	for (int i = 0; i < this->dataInfo.nSelectedTracks; i++)
	{
		trackid = this->dataInfo.selectedTracks.at(i);
		length =  this->TracksInfo.at(trackid).nPoints;

		for (int j = 0; j<length; j++)
		{
			pid = this->TracksInfo.at(trackid).PointsIds.at(j);

			if (pid>-1)// && idMap1->at(pid)>-1)
			{
				if (idMap1->at(pid)==-1) // check if point not already inserted
				{
					idMap1->at(pid) = counter;
					idMap2->at(counter) = pid;
					counter++;
				}
			}
		}

	}

	this->dataInfo.nAllSelectedPoints = counter;
	idMap2->resize(counter);

	
	return 1;
}*/

/*----------------------------------------------------------------------------
generates points in the output vector
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
/*int vtkSQLiteReader2::generatePoints(SelectionStruct* selection, Data* actData)
{
	std::vector<int> * idMap2 = &selection->vPointIdMapReverse;

	actData->Position->Reset();
	actData->Cells->Reset();
	
	actData->Position->SetNumberOfPoints(selection->nSelectedPoints);
	actData->Cells->SetNumberOfCells(selection->nSelectedPoints);
	actData->Velocity->SetNumberOfTuples(selection->nSelectedPoints);
	actData->TrackId->SetNumberOfTuples(selection->nSelectedPoints);
	actData->GId->SetNumberOfTuples(selection->nSelectedPoints);
	actData->SnapId->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Mvir->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Rvir->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Vmax->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Rmax->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Redshift->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Time->SetNumberOfTuples(selection->nSelectedPoints);
	actData->Cv->SetNumberOfTuples(selection->nSelectedPoints);
	actData->CvAverage->SetNumberOfTuples(selection->nSelectedPoints);
	actData->CollisionTypePoint->SetNumberOfTuples(selection->nSelectedPoints);
	actData->CollisionTypeTrack->SetNumberOfTuples(selection->nSelectedPoints);
	actData->RGc->SetNumberOfTuples(selection->nSelectedPoints);

	for (int i = 0; i<selection->nSelectedPoints; i++)
	{
		actData->Position->InsertPoint(i,
			this->allData.Position->GetPoint(idMap2->at(i)));
		actData->Velocity->InsertTuple(i,
			this->allData.Velocity->GetTuple(idMap2->at(i)));
		actData->TrackId->InsertTuple(i,
			this->allData.TrackId->GetTuple(idMap2->at(i)));
		actData->GId->InsertTuple(i,
			this->allData.GId->GetTuple(idMap2->at(i)));
		actData->SnapId->InsertTuple(i,
			this->allData.SnapId->GetTuple(idMap2->at(i)));
		actData->Mvir->InsertTuple(i,
			this->allData.Mvir->GetTuple(idMap2->at(i)));
		actData->Rvir->InsertTuple(i,
			this->allData.Rvir->GetTuple(idMap2->at(i)));
		actData->Vmax->InsertTuple(i,
			this->allData.Vmax->GetTuple(idMap2->at(i)));
		actData->Rmax->InsertTuple(i,
			this->allData.Rmax->GetTuple(idMap2->at(i)));
		actData->Redshift->InsertTuple(i,
			this->allData.Redshift->GetTuple(idMap2->at(i)));
		actData->Time->InsertTuple(i,
			this->allData.Time->GetTuple(idMap2->at(i)));
		actData->Cv->InsertTuple(i,
			this->allData.Cv->GetTuple(idMap2->at(i)));
		actData->CvAverage->InsertTuple(i,
			this->allData.CvAverage->GetTuple(idMap2->at(i)));
		actData->CollisionTypePoint->InsertTuple(i,
			this->allData.CollisionTypePoint->GetTuple(idMap2->at(i)));
		actData->CollisionTypeTrack->InsertTuple(i,
			this->allData.CollisionTypeTrack->GetTuple(idMap2->at(i)));
		actData->RGc->InsertTuple(i,
			this->allData.RGc->GetTuple(idMap2->at(i)));
	}

	// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = actData->Position->GetNumberOfPoints();
	vtkIdType *cells = actData->Cells->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		cells[i*2]   = 1;
		cells[i*2+1] = i;
    }

	return 1;
}*/



/*----------------------------------------------------------------------------
generates tracks in the output vector
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
/*int vtkSQLiteReader2::generateTracks(SelectionStruct* selection, Data* actData)
{
	int trackId, oldId, newId;
	vtkSmartPointer<vtkPolyLine> nextLine;

	actData->Tracks->Reset();

	//this->allData.Tracks->GetCell(

	for (int i = 0; i<selection->nSelectedTracks; i++)
	{
		nextLine = vtkSmartPointer<vtkPolyLine>::New();
		trackId = selection->vTrackIds.at(i);

		for (int j = 0; j<this->TracksInfo.at(trackId).nPoints; j++)
		{
			oldId = this->TracksInfo.at(trackId).PointsIds.at(j);
			newId = selection->vPointIdMap.at(oldId);
			nextLine->GetPointIds()->InsertNextId(newId);
		}
		actData->Tracks->InsertNextCell(nextLine);
	}

	return 1;
}*/

/*----------------------------------------------------------------------------
Finds the collision and merging Points

	assumes:
		
	sets:
		int	nCollisions: Number of collisions
	arguments:
		pCollCalc		pointer to struct with settings and to fill the data in
	returns:
		int errorcode
*/
int vtkSQLiteReader2::findCollisions(CollisionCalculation* pCollCalc){



	//double candidate[3];
	//double dist;
	//double CollTol = pCollCalc->CollisionTolerance;
	//double MergTol = pCollCalc->MergingTolerance;

	////int offset, length;
	//int PointId, SnapId, TrackId, useid;
	vtkIdType SnapId;

	//int nMerging = 0;
	//int nColliding = 0;

	vtkSmartPointer<vtkPoints> points;
	//vtkSmartPointer<vtkKdTree> kDTree;
	vtkSmartPointer<vtkIdList> IdsOfSnap = vtkSmartPointer<vtkIdList>::New();
	//vtkSmartPointer<vtkIdList> IdList;

	/* dont use this, use id list instead..
	// vectores to store the points/tracks / snaps of interesst (indexed by id)
	std::vector<vtkIdType> Points; //points of interest
	std::vector<vtkIdType> Tracks; // tracks of interest
	std::vector<vtkIdType> Snapshots; // 

	// stores info about tracks and points
	//TODO check whtas faster, vector or id list...
	//  -1 unchecked
	//	*0 checked, no collision
	//	1 collision on this point / track
	//	2 merging event (a smaller halo merges in thisone, track has such an event
	//	*3 is the minor of two mergin halos (points only)
	Points.resize(this->dataInfo.nPoints,-1);
	Tracks.resize(this->dataInfo.nTracks,-1);
	Snapshots.resize(this->dataInfo.nSnapshots,-1);
	*/
	vtkSmartPointer<vtkIdList> CollisionPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> CollisionTracks = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> MergingPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> MergingTracks = vtkSmartPointer<vtkIdList>::New();


	for (SnapId=0; SnapId < this->dataInfo.nSnapshots; SnapId++)
	{
		// filter out snapshots with 0 points and possibly other nasty stuff..
		IdsOfSnap->DeepCopy(this->SnapInfo2.at(SnapId));
		IdsOfSnap->DeleteId(-1);
		if (IdsOfSnap->GetNumberOfIds() <= 1){continue;} // dont do this for only 0 or 1 point

		// select points
		points = vtkSmartPointer<vtkPoints>::New();
		this->AllData->GetPoints()->GetPoints(IdsOfSnap,points);

		// build kd tree
		vtkSmartPointer<vtkKdTree> kDTree = vtkSmartPointer<vtkKdTree>::New();
		kDTree->BuildLocatorFromPoints(points);

		// check for best variant:
		//  - FindClosestNPoints
		//  - BuildMapForDuplicatePoints 

		vtkIdTypeArray * result;
		vtkIdType pointid;

		// FIND COLLIDING POINTS
		result = kDTree->BuildMapForDuplicatePoints(pCollCalc->CollisionTolerance);

		for (vtkIdType id = 0; id<result->GetNumberOfTuples(); ++id)
		{
			if (id != result->GetTuple1(id))
			{
				//SPEEDUP: maybe just add ids and check later for uniqueness
				CollisionPoints->InsertUniqueId(IdsOfSnap->GetId(id));
				CollisionPoints->InsertUniqueId(IdsOfSnap->GetId(result->GetTuple1(id)));
			}
		}
		
		/*vtkErrorMacro("Printing the result of kd tree:");
		for (int i = 0; i<result->GetNumberOfTuples(); ++i)
		{
			vtkErrorMacro(" " << i << ": " << result->GetTuple1(i));
		}*/

		result->Delete(); // free memory

		// FIND MERGING POINTS
		result = kDTree->BuildMapForDuplicatePoints(pCollCalc->MergingTolerance);

		for (vtkIdType id = 0; id<result->GetNumberOfTuples(); ++id)
		{
			if (id != result->GetTuple1(id))
			{
				pointid = IdsOfSnap->GetId(id);
				MergingPoints->InsertUniqueId(pointid);
				MergingTracks->InsertUniqueId(this->dataInfo.pTrackIdArray->GetTuple1(pointid));

				pointid = IdsOfSnap->GetId(result->GetTuple1(id));
				MergingPoints->InsertUniqueId(pointid);
				MergingTracks->InsertUniqueId(this->dataInfo.pTrackIdArray->GetTuple1(pointid));
			}
		}
		result->Delete(); // free memory

	}
		
	//DEBUG OUTPUT
	/*
	vtkErrorMacro("Printing the merging track ids (before cleanup):");
	for (int i = 0; i<MergingTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingTracks->GetId(i));}
	vtkErrorMacro("Printing the merging point ids (before cleanup):");
	for (int i = 0; i<MergingPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingPoints->GetId(i));}
	vtkErrorMacro("Printing the colliding point ids (before cleanup):");
	for (int i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionPoints->GetId(i));}
	*/



	//OUTPUT COLLIDING POINTS ONLY
	this->IdListComplement(CollisionPoints, MergingPoints);

	// find colliding, but not merging tracks
	vtkIdType trackid;
	for (vtkIdType i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{
		trackid = this->dataInfo.pTrackIdArray->GetTuple1(CollisionPoints->GetId(i));
		if (MergingTracks->IsId(trackid) == -1)
		{
			CollisionTracks->InsertUniqueId(trackid);
		}
	}

	// DEBUG OUTPUT
	/*
	vtkErrorMacro("Printing the merging track ids (after cleanup):");
	for (int i = 0; i<MergingTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingTracks->GetId(i));}
	vtkErrorMacro("Printing the Colliding track ids (after cleanup):");
	for (int i = 0; i<CollisionTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionTracks->GetId(i));}
	
	vtkErrorMacro("Printing the merging point ids (after cleanup):");
	for (int i = 0; i<MergingPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingPoints->GetId(i));}
	vtkErrorMacro("Printing the colliding point ids (after cleanup):");
	for (int i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionPoints->GetId(i));}
	*/
	
	/*
	vtkErrorMacro("No Merging Tracks: " << MergingTracks->GetNumberOfIds());
	vtkErrorMacro("No Collision Tracks: " << CollisionTracks->GetNumberOfIds());
	vtkErrorMacro("No Merging Points: " << MergingPoints->GetNumberOfIds());
	vtkErrorMacro("No Collision Points: " << CollisionPoints->GetNumberOfIds());
	*/


	
		/*vtkErrorMacro("the selected points ids are: 1 :");
		for (int i = 0; i<this->SnapInfo2.at(SnapId)->GetNumberOfIds(); ++i)
		{
			vtkErrorMacro(" I  " << i << ": id=" << this->SnapInfo2.at(SnapId)->GetId(i));
		}
		vtkErrorMacro("the selected points are::");
		for (int i = 0; i<points->GetNumberOfPoints(); ++i)
		{
			vtkErrorMacro(" II " << i << ": x=" << points->GetData()->GetComponent(i,0));
		}
		vtkErrorMacro("Printing the result of kd tree:");
		for (int i = 0; i<result->GetNumberOfTuples(); ++i)
		{
			vtkErrorMacro(" " << i << ": " << result->GetTuple1(i));
		}*/
		
		
	
	pCollCalc->CollisionPointIds = CollisionPoints;
	pCollCalc->CollisionTrackIds = CollisionTracks;
	pCollCalc->MergingPointIds = MergingPoints;
	pCollCalc->MergingTrackIds = MergingTracks;
	pCollCalc->CalcIsDone = true;

	/*
	vtkErrorMacro("Printing the colliding track ids: (findColl fnc)");
	for (int i = 0; i<CollisionTracks->GetNumberOfIds(); ++i)
	{
		vtkErrorMacro(" " << i << ": " << CollisionTracks->GetId(i));
	}
	*/


	return 1;
}
		
		//result->get

/*
		for (int PointidToCheck = 0;
			PointidToCheck < this->SnapInfo2.at(SnapId)->GetNumberOfIds();
			++PointidToCheck)
		{

			kDTree->FindClosestNPoints(2,candidate,IdList);

		}




		for (TrackId = 0; TrackId < this->dataInfo.nTracks; TrackId++)
		{
			// go on if theres already a collision on this track, but then you only get the
			// first point on a track wheres a collision..
			// if(Tracks.at(TrackId)==1){continue;}

			PointId = this->TracksInfo.at(TrackId).PointsIds.at(SnapId);
			if(PointId == -1){continue;}
			this->allData.Position->GetPoint(PointId, candidate);
			
			kDTree->FindClosestNPoints(2,candidate,IdList);
			
			if(IdList->GetNumberOfIds()<2){continue;}
			
			// check for identical points, select the other point
			// TODO make this better, check for identical coordinates too..
			if (this->allData.Mvir->GetTuple1(PointId) != this->allData.Mvir->GetTuple1(IdList->GetId(1)))
			{
				useid = 1;
			}
			else if (this->allData.Mvir->GetTuple1(PointId) != this->allData.Mvir->GetTuple1(IdList->GetId(0)))
			{
				useid = 0;
			}
			else // this shouln't happen
			{
				continue;
			}

			//dist = getDistance2(PointId,IdList->GetId(useid)+offset);
			dist = 1;
			
			//vtkErrorMacro(" for track: " << TrackId
			//	<< ", got point ids: " << IdList->GetId(0) << ", " << IdList->GetId(1)
			//	<< " with dist: " << dist);

			if(dist<uppTol)
			{
				//vtkErrorMacro("COLLISION: for track: " << TrackId
				//	<< " in snap: " << SnapId
				//	<< ", got point ids: " << IdList->GetId(1)
				//	<< " with dist: " << dist);
				
				if (dist<=lowTol)
				{
					nMerging++;

					Tracks.at(TrackId) = 2;
					Points.at(PointId) = 2;
					//Snapshots.at(SnapId) = 2;
				}
				else
				{
					nColliding++;
					if (Tracks.at(TrackId)!=2){Tracks.at(TrackId) = 1;}
					Points.at(PointId) = 1;
					//Snapshots.at(SnapId) = 1;
				}
			}
		}
	}

	CollisionResultStruct * result;
	int number;
	int state;
	for (int i=0; i<2; i++)
	{
		if(i==0){
			result = &pCollCalc->Colliding;
			number = nColliding;
			state = 1; //colliding tracks
		}
		else if (i==1){
			result = &pCollCalc->Merging;
			number = nMerging;
			state = 2; //merging tracks
		}

		result->nEvents = number;
		result->nPoints = 0;
		result->nTracks = 0;
		result->vPointIds.clear();
		result->vTrackIds.clear();

		for (int j=0; j < this->dataInfo.nTracks; j++)
		{
			if (Tracks.at(j)==state)
			{
				result->nTracks ++;
				result->vTrackIds.push_back(j);
				for (int l = 0; l<this->TracksInfo.at(j).nPoints; l++)
				{
					this->allData.CollisionTypeTrack->InsertTuple1(this->TracksInfo.at(j).PointsIds.at(l),state);
				}
			}
		}

		for (int j=0;j<this->dataInfo.nPoints;j++)
		{
			if (Points.at(j)==state)
			{
				result->nPoints ++;
				result->vPointIds.push_back(j);
				this->allData.CollisionTypePoint->InsertTuple1(j,state);
			}
		}
	}

*/


/* old version ---
	double candidate[3];
	double dist;
	double uppTol = pCollCalc->upperTolerance;
	double lowTol = pCollCalc->lowerTolerance;

	int offset, length;
	int PointId, SnapId, TrackId, useid;

	int nMerging = 0;
	int nColliding = 0;

	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkKdTree> kDTree;
	vtkSmartPointer<vtkIdList> IdList;

	// vectores to store the points/tracks / snaps of interesst (indexed by id)
	std::vector<int> Points; //points of interest
	std::vector<int> Tracks; // tracks of interest
	std::vector<int> Snapshots; // 

	// stores info about tracks and points
	//  -1 unchecked
	//	0 checked, no collision
	//	1 collision on this point / track
	//	2 merging event (a smaller halo merges in thisone, track has such an event
	//	3 is the minor of two mergin halos (points only)
	Points.resize(this->dataInfo.nPoints,-1);
	Tracks.resize(this->dataInfo.nTracks,-1);
	Snapshots.resize(this->dataInfo.nSnapshots,-1);

	for (SnapId=0; SnapId < this->dataInfo.nSnapshots; SnapId++)
	{
		// get all points from current snapshot
		offset = this->SnapInfo.at(SnapId).Offset;
		length = this->SnapInfo.at(SnapId).lenght;
		
		// copy data of interest
		points = vtkSmartPointer<vtkPoints>::New();
		points->SetNumberOfPoints(length);
		this->allData.Position->GetData()->GetTuples(offset, offset+length-1, points->GetData());

		//build kd tree
		kDTree = vtkSmartPointer<vtkKdTree>::New();
		kDTree->BuildLocatorFromPoints(points);

		// init idlist for result		
		IdList = vtkSmartPointer<vtkIdList>::New();

		for (TrackId = 0; TrackId < this->dataInfo.nTracks; TrackId++) //
		{
			// go on if theres already a collision on this track, but then you only get the
			// first point on a track wheres a collision..
			// if(Tracks.at(TrackId)==1){continue;} 

			PointId = this->TracksInfo.at(TrackId).PointsIds.at(SnapId);
			if(PointId == -1){continue;}
			this->allData.Position->GetPoint(PointId, candidate);
			
			kDTree->FindClosestNPoints(2,candidate,IdList);
			
			if(IdList->GetNumberOfIds()<2){continue;}
			
			// check for identical points, select the other point
			// TODO make this better, check for identical coordinates too..
			if (this->allData.Mvir->GetTuple1(PointId) != this->allData.Mvir->GetTuple1(IdList->GetId(1)))
			{
				useid = 1;
			}
			else if (this->allData.Mvir->GetTuple1(PointId) != this->allData.Mvir->GetTuple1(IdList->GetId(0)))
			{
				useid = 0;
			}
			else // this shouln't happen
			{
				continue;
			}

			dist = getDistance2(PointId,IdList->GetId(useid)+offset);
			
			//vtkErrorMacro(" for track: " << TrackId
			//	<< ", got point ids: " << IdList->GetId(0) << ", " << IdList->GetId(1)
			//	<< " with dist: " << dist);

			if(dist<uppTol)
			{
				//vtkErrorMacro("COLLISION: for track: " << TrackId
				//	<< " in snap: " << SnapId
				//	<< ", got point ids: " << IdList->GetId(1)
				//	<< " with dist: " << dist);
				
				if (dist<=lowTol)
				{
					nMerging++;

					Tracks.at(TrackId) = 2;
					Points.at(PointId) = 2;
					//Snapshots.at(SnapId) = 2;
				}
				else
				{
					nColliding++;
					if (Tracks.at(TrackId)!=2){Tracks.at(TrackId) = 1;}
					Points.at(PointId) = 1;
					//Snapshots.at(SnapId) = 1;
				}
			}
		}
	}

	CollisionResultStruct * result;
	int number;
	int state;
	for (int i=0; i<2; i++)
	{
		if(i==0){
			result = &pCollCalc->Colliding;
			number = nColliding;
			state = 1; //colliding tracks
		}
		else if (i==1){
			result = &pCollCalc->Merging;
			number = nMerging;
			state = 2; //merging tracks
		}

		result->nEvents = number;
		result->nPoints = 0;
		result->nTracks = 0;
		result->vPointIds.clear();
		result->vTrackIds.clear();

		for (int j=0; j < this->dataInfo.nTracks; j++)
		{
			if (Tracks.at(j)==state)
			{
				result->nTracks ++;
				result->vTrackIds.push_back(j);
				for (int l = 0; l<this->TracksInfo.at(j).nPoints; l++)
				{
					this->allData.CollisionTypeTrack->InsertTuple1(this->TracksInfo.at(j).PointsIds.at(l),state);
				}
			}
		}

		for (int j=0;j<this->dataInfo.nPoints;j++)
		{
			if (Points.at(j)==state)
			{
				result->nPoints ++;
				result->vPointIds.push_back(j);
				this->allData.CollisionTypePoint->InsertTuple1(j,state);
			}
		}
	}

	return 1;
	*/


/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader2::generateSelection(CollisionCalculation * collCalc, SelectionStruct* select)
{
	select->Points = vtkSmartPointer<vtkIdList>::New();
	select->Tracks = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> allTrackIds = vtkSmartPointer<vtkIdList>::New();
	
	// generate a list of all trackids
	//could be done faster, by assuming ids go from 0 to nTracks... just wanna play safe..
	/* cancel this, it's too slow...
	this->IdTypeArray2IdList(allTrackIds, this->dataInfo.pTrackIdArray);
	*/
	allTrackIds->SetNumberOfIds(this->dataInfo.nTracks);
	for (vtkIdType i = 0; i < this->dataInfo.nTracks; ++i)
	{
		allTrackIds->SetId(i,i);
	}
	
	/*
	vtkErrorMacro("Printing the selected track ids (be4 union):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i)<<"\n");}
	
	vtkErrorMacro("Printing the selected track ids (be4 union):");
	for (int i = 0; i<collCalc->CollisionTrackIds->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << collCalc->CollisionTrackIds->GetId(i)<<"\n");}
	*/

	if(*this->Gui.DisplayColliding)
	{
		select->Tracks = this->IdListUnionUnique(select->Tracks, collCalc->CollisionTrackIds);
	}

	/*vtkErrorMacro("Printing the selected track ids (after union):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i)<<"\n");}*/


	if(*this->Gui.DisplayMerging)
	{
		select->Tracks = this->IdListUnionUnique(select->Tracks, collCalc->MergingTrackIds);
	}

	if(*this->Gui.DisplayInverted)
	{
		select->Tracks = this->IdListComplement(select->Tracks, allTrackIds);
	}

	/*vtkErrorMacro("Printing the selected track ids (after invert):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i));}*/


	// collect the points, corresponding to the selected tracks
	for (vtkIdType i = 0 ; i < select->Tracks->GetNumberOfIds(); ++i)
	{
		vtkIdList *tmplist = this->TrackInfo2.at(select->Tracks->GetId(i));
		for (vtkIdType j = 0;
			j<tmplist->GetNumberOfIds();
			++j)
		{
			//SPEEDCHECK01
			//select->Points->InsertUniqueId(tmplist->GetId(j));
			select->Points->InsertNextId(tmplist->GetId(j));
		}
	}

	/*vtkErrorMacro("Printing the selected points:");
	for (int i = 0; i<select->Points->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Points->GetId(i));}*/


	return 1;
}

/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader2::generateSelectedData(
		vtkSmartPointer<vtkPolyData> output,
		SelectionStruct* selection)
{
	// generate points
	this->AllData->GetPoints()->GetPoints(
		selection->Points,
		output->GetPoints());
	
	// generate cells
	vtkIdType N = selection->Points->GetNumberOfIds();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType *pverts = verts->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		pverts[i*2]   = 1;
		pverts[i*2+1] = i;
    }
	output->SetVerts(verts);

	// generate tracks
	vtkSmartPointer<vtkCellArray> Tracks = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIdList> nextLine;
	//vtkIdType oldtrackid;
	//vtkIdList * oldTrack;
	
	int offset = 0;

	for (vtkIdType track = 0; track<selection->Tracks->GetNumberOfIds(); ++track)
	{
		//this is a speedhack, not sure how stable it is..
		vtkIdType nPnts = this->TrackInfo2.at(selection->Tracks->GetId(track))->GetNumberOfIds();
		Tracks->InsertNextCell(nPnts);
		for (vtkIdType i = offset; i<offset+nPnts; ++i)
		{
			Tracks->InsertCellPoint(i);
		}
		offset += nPnts;

		/* Here's the stable way to do it..
		oldTrack = this->TrackInfo2.at(selection->Tracks->GetId(track));
		nextLine = vtkSmartPointer<vtkIdList>::New();
		nextLine->SetNumberOfIds(oldTrack->GetNumberOfIds());

		for (vtkIdType pnt = 0;
			pnt < oldTrack->GetNumberOfIds();
			pnt++)
		{
			nextLine->SetId(pnt,
				selection->Points->IsId(	//gets newid
				oldTrack->GetId(pnt))); //gets oldid

			// if you wanna be really sure nothing goes wrong
			//if (newid>-1){nextLine->InsertId(pnt,newid);}
			//else {vtkErrorMacro("oldid not found");}
		}

		Tracks->InsertNextCell(nextLine);
		*/
	}

	output->SetLines(Tracks);


	// generate pointdata
	vtkPointData * PDall = this->AllData->GetPointData();
	vtkPointData * PDout = output->GetPointData();
	int nPnts = selection->Points->GetNumberOfIds();

	for (int i = 0; i<PDall->GetNumberOfArrays(); ++i)
	{
		PDout->GetArray(i)->SetNumberOfTuples(nPnts);
		PDall->GetArray(i)->GetTuples(selection->Points, PDout->GetArray(i));
	}

	return 1;
}


/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		p_map
		p_count
		t_map
		t_count
		result
		type		type of event (1 = collision, 2 = merging)
		
	returns:
		int errorcode
*/
/*int vtkSQLiteReader2::fillIdList(std::vector<int> * p_map, int * p_count,
								std::vector<int> * t_map, int * t_count, 
								CollisionResultStruct* result, int type)
{
	int i, j, pid, trackid, length;

	//fill in tracks
	for (i = 0; i < result->nTracks; i++)
	{
		trackid = result->vTrackIds.at(i);
		length =  this->TracksInfo.at(trackid).nPoints;

		if (t_map->at(trackid)>-1){continue;}
		
		t_map->at(trackid) = *t_count;
		*t_count = *t_count +1;

		for (j = 0; j<length; j++)
		{
			pid = this->TracksInfo.at(trackid).PointsIds.at(j);

			if (pid>-1 && p_map->at(pid)==-1) // if valid point on track AND check if point not already inserted
			{
				p_map->at(pid) = *p_count;
				*p_count = *p_count + 1;
			}
		}

	}
	return 1;
*/
	//fill in points
	/* VERY OLD STUFF HERE
	for (i = 0; i<result->nPoints; i++)
	{
		if (p_map->at(result->vPointIds.at(i)) ==-1)
		{
			trackid = this->allData.TrackId->GetTuple1(result->vPointIds.at(i))
			if(t_map->at(trackid) == type)
			{
				p_map->at(result->vPointIds.at(i)) = *p_count;
				*p_count = *p_count + 1;
			}
		}
	}
	*/





//}


/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkFloatArray> vtkSQLiteReader2::CreateArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkFloatArray> dataArray=vtkSmartPointer<vtkFloatArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkIntArray> vtkSQLiteReader2::CreateIntArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkIntArray> dataArray=vtkSmartPointer<vtkIntArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkIdTypeArray> vtkSQLiteReader2::CreateIdTypeArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkIdTypeArray> dataArray=vtkSmartPointer<vtkIdTypeArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Initialises all array
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/
int vtkSQLiteReader2::InitAllArrays(
		vtkDataSet *output,
		unsigned long numTuples)
{
	vtkSmartPointer<vtkDataArray> dataArray;

	for (int i = 0; i<output->GetPointData()->GetNumberOfArrays(); ++i)
	{
		dataArray = output->GetPointData()->GetArray(i);
		dataArray->SetNumberOfTuples(numTuples);
		for (int j=0; j < dataArray->GetNumberOfComponents(); ++j) {
			dataArray->FillComponent(j, 0.0);
		}
	}
	return 1;
}

/*----------------------------------------------------------------------------
Unions two IdLists, adds list 2 to list 1
( set notation: list1 u list2, list1 or list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader2::IdListUnion(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	int n1 = list1->GetNumberOfIds();
	int n2 = list1->GetNumberOfIds();
	list1->SetNumberOfIds(n1+n2);
	for (vtkIdType i = 0; i<n2; ++i)
	{
		list1->SetId(n1+i,list2->GetId(i));
	}
	return list1;
}

/*----------------------------------------------------------------------------
Unions two IdLists, adds list 2 to list 1
and makes sure, each id is unique
list1 for uniqueness!
( set notation: list1 u list2, list1 or list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader2::IdListUnionUnique(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
	//result->SetNumberOfIds(list1->GetNumberOfIds()+list2->GetNumberOfIds());
	
	for (vtkIdType i = 0; i<list1->GetNumberOfIds(); ++i)
	{
		//vtkErrorMacro(" union: l1: id-"<<i<<" - "<<list1->GetId(i));
		result->InsertUniqueId(list1->GetId(i));
	}
	for (vtkIdType i = 0; i<list2->GetNumberOfIds(); ++i)
	{
		//vtkErrorMacro(" union: l2: id-"<<i<<" - "<<list2->GetId(i));
		result->InsertUniqueId(list2->GetId(i));
	}
	result->Squeeze();
	//list1 = result;
	return result;

}
/*----------------------------------------------------------------------------
Intersects two IdLists, result in list1 (wrapper)
( set notation: list1 n list2, list1 and list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader2::IdListIntersect(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	list1->IntersectWith(*list2);
	return list1;
}
/*----------------------------------------------------------------------------
Gets the complement of two IdLists, result in list1
or: deletes all ids from list 2 in list 1
( set notation: list1 \ list2, list1 - list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader2::IdListComplement(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	for (vtkIdType i = 0; i<list2->GetNumberOfIds(); ++i)
	{
		list1->DeleteId(list2->GetId(i));
	}
	list1->Squeeze();
	return list1;
}

/*----------------------------------------------------------------------------
Converts a vtkidTypeArray to a vtkIdList
(adds only unique ids to the list)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

int vtkSQLiteReader2::IdTypeArray2IdList(vtkIdList* destination, vtkIdTypeArray* source)
{
	for (vtkIdType i = 0; i<source->GetNumberOfTuples(); ++i)
	{
		destination->InsertUniqueId(source->GetTuple1(i));
	}
	return 1;
}


/*----------------------------------------------------------------------------
Adds/calculates aditional point data
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
int vtkSQLiteReader2::calculatePointData()
{
	vtkDataArray * cvarr = this->AllData->GetPointData()->GetArray("Cv");
	vtkDataArray * rsarr = this->AllData->GetPointData()->GetArray("Redshift");
	vtkDataArray * aarr = this->AllData->GetPointData()->GetArray("a");

	cvarr->SetNumberOfTuples(this->dataInfo.nPoints);
	rsarr->SetNumberOfTuples(this->dataInfo.nPoints);
	aarr->SetNumberOfTuples(this->dataInfo.nPoints);
	
	//cvarr->FillComponent(0,0);

	vtkDataArray * Vmax = this->AllData->GetPointData()->GetArray("Vmax");
	vtkDataArray * Rmax = this->AllData->GetPointData()->GetArray("Rmax");

	vtkDataArray * ssdrsarr = this->SnapshotData->GetPointData()->GetArray("redshift");

	float cv, rs, a;
	int snapid;

	for(vtkIdType i = 0; i<this->dataInfo.nPoints; ++i)
	{
		if (Rmax->GetTuple1(i) == 0)
		{
			cv = 0;
		}
		else
		{
			cv = 2 * ( Vmax->GetTuple1(i) / (this->dataInfo.Hubble * Rmax->GetTuple1(i))) *
				( Vmax->GetTuple1(i) / (this->dataInfo.Hubble * Rmax->GetTuple1(i)));
		}
		cvarr->SetTuple1(i,cv);

		snapid = this->AllData->GetPointData()->GetArray(this->dataInfo.StatSnapidColumn)->GetTuple1(i);
		
		rs = ssdrsarr->GetTuple1(snapid);
		rsarr->SetTuple1(i,rs);
		
		a=1/(1+rs);
		aarr->SetTuple1(i,a);
	}

	return 1;
}