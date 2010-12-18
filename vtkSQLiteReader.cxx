/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.cxx,v $

  
	todo:
	- debug for wrong input of nsnapshot, resp. input from display snapshot nr is > nsnapshot
	- check for tracksdisplay
	- maybe more computationonal stuff
	- option for estimation good tolerance parameter
	bugfixing:
	- lost color bug
	- check resetting of variables, whats needed, what not..  
  
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
#include "vtkKdTree.h"
#include "vtkTable.h"


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

	/*
	
	Gui.DisplayOnlySelectedData = &(this->DisplayOnlySelectedData); 
	Gui.DisplaySelected = &(this->DisplaySelected); 
	Gui.DisplaySelectedSnapshot = &(this->DisplaySelectedSnapshot); 
	Gui.DisplaySelectedSnapshotNr = &(this->DisplaySelectedSnapshotNr); 
	Gui.DisplaySelectedTrack = &(this->DisplaySelectedTrack); 
	Gui.DisplaySelectedTrackNr = &(this->DisplaySelectedTrackNr); 
	Gui.DisplayCalculated = &(this->DisplayCalculated);
	Gui.CalculationImpactParameter = &(this->CalculationImpactParameter);
	Gui.DisplayEstimateTolerance = &(this->DisplayEstimateTolerance);

	Gui.calcHighlightCollisionPoints = &(this->calcHighlightCollisionPoints);
	Gui.calcHighlightTrack = &(this->calcHighlightTrack);
	Gui.calcHighlightSnapshot = &(this->calcHighlightSnapshot);
	Gui.calcHighlightPoint = &(this->calcHighlightPoint);
	Gui.calcHighlightTrackNr = &(this->calcHighlightTrackNr);
	Gui.calcHighlightSnapshotNr = &(this->calcHighlightSnapshotNr);
	Gui.calcHighlightPointNr = &(this->calcHighlightPointNr);

	*Gui.DisplayOnlySelectedData = false; 
	*Gui.DisplaySelected = false; 
	*Gui.DisplaySelectedSnapshot = false; 
	*Gui.DisplaySelectedSnapshotNr = -1; 
	*Gui.DisplaySelectedTrack = false; 
	*Gui.DisplaySelectedTrackNr = -1; 
	*Gui.DisplayCalculated = false;
	*Gui.CalculationImpactParameter = 0.0000001;
	*Gui.DisplayEstimateTolerance = false;
	*/

	// init dataInfo (maybee just call this->reset!??)
	dataInfo.dataIsRead = false;
	dataInfo.nPoints = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.nTracks = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.hubble = 0;

	/*dataInfo.nSelectedPoints = -1;
	dataInfo.nSelectedTracks = -1;
	dataInfo.nSelectedSnapshots = -1;
	dataInfo.selectedPoints.clear();
	dataInfo.selectedSnapshots.clear();
	dataInfo.selectedTracks.clear();*/

	//init data stuff

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

	}

	this->collisionCalc.lowerTolerance = *this->Gui.LowerLimit * *this->Gui.LowerLimit;
	this->collisionCalc.upperTolerance = *this->Gui.UpperLimit * *this->Gui.UpperLimit;
	this->collisionCalc.isDone = false;

	/*
			
			
			
			
	this->allData.Position = vtkSmartPointer<vtkPoints>::New();
	this->allData.Velocity = vtkSmartPointer<vtkFloatArray>::New();
	this->allData.Cells = vtkSmartPointer<vtkCellArray>::New();
	this->allData.Tracks = vtkSmartPointer<vtkCellArray>::New();
	this->allData.TrackId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->allData.GId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->allData.SnapId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->allData.RVir = vtkSmartPointer<vtkFloatArray>::New();
	this->allData.RVir = vtkSmartPointer<vtkFloatArray>::New();
	this->allData.RVir = vtkSmartPointer<vtkFloatArray>::New();
	this->allData.RVir = vtkSmartPointer<vtkFloatArray>::New();
	this->allData.Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();

	this->allData.Velocity->SetNumberOfComponents(3);
	this->allData.TrackId->SetNumberOfComponents(1);
	this->allData.GId->SetNumberOfComponents(1);
	this->allData.SnapId->SetNumberOfComponents(1);
	this->allData.RVir->SetNumberOfComponents(1);
	this->allData.RVir->SetNumberOfComponents(1);
	this->allData.RVir->SetNumberOfComponents(1);
	this->allData.RVir->SetNumberOfComponents(1);
	this->allData.Colors->SetNumberOfComponents(3);

	this->allData.Velocity->SetName("Velocity");
	this->allData.TrackId->SetName("TrackId");
	this->allData.GId->SetName("GId");
	this->allData.SnapId->SetName("SnapId");
	this->allData.RVir->SetName("RVir");
	this->allData.Colors->SetName("Colors");

	this->selectedData.Position = vtkSmartPointer<vtkPoints>::New();
	this->selectedData.Velocity = vtkSmartPointer<vtkFloatArray>::New();
	this->selectedData.Cells = vtkSmartPointer<vtkCellArray>::New();
	this->selectedData.Tracks = vtkSmartPointer<vtkCellArray>::New();
	this->selectedData.TrackId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->selectedData.GId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->selectedData.SnapId = vtkSmartPointer<vtkIdTypeArray>::New();
	this->selectedData.RVir = vtkSmartPointer<vtkFloatArray>::New();
	this->selectedData.Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();

	this->selectedData.Velocity->SetNumberOfComponents(3);
	this->selectedData.TrackId->SetNumberOfComponents(1);
	this->selectedData.GId->SetNumberOfComponents(1);
	this->selectedData.SnapId->SetNumberOfComponents(1);
	this->selectedData.RVir->SetNumberOfComponents(1);
	this->selectedData.Colors->SetNumberOfComponents(3);

	this->selectedData.Velocity->SetName("Velocity");
	this->selectedData.TrackId->SetName("TrackId");
	this->selectedData.GId->SetName("GId");
	this->selectedData.SnapId->SetName("SnapId");
	this->selectedData.RVir->SetName("RVir");
	this->selectedData.Colors->SetName("Colors");*/

	//init calc struct
	this->calcInfo.calcDone = false;
	this->calcInfo.nCollisions = -1;
	this->calcInfo.tolerance = 0;

	// inti calcesttol2
	this->calcEstTol2 = vtkSmartPointer<vtkFloatArray>::New();
	this->calcEstTol2->SetNumberOfComponents(4);
	this->calcEstTol2->SetName("CalculationOfTolerance");
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
	//vtkTable * table = vtkTable::GetData(outputVector);
	//table->Initialize();

	//vtkDebugMacro("test");

	if (!this->dataInfo.dataIsRead)
		// only read if it's not been read before
	{
		readSnapshotInfo();
		readSnapshots();
		readTracks();

		this->dataInfo.dataIsRead = true;
	}

	if (this->collisionCalc.lowerTolerance != *this->Gui.LowerLimit * *this->Gui.LowerLimit ||
		this->collisionCalc.upperTolerance != *this->Gui.UpperLimit * *this->Gui.UpperLimit)
	{
		this->collisionCalc.lowerTolerance = *this->Gui.LowerLimit * *this->Gui.LowerLimit;
		this->collisionCalc.upperTolerance = *this->Gui.UpperLimit * *this->Gui.UpperLimit;
		this->collisionCalc.isDone = false;
	}

	if(!this->collisionCalc.isDone && (*Gui.DisplayColliding || *Gui.DisplayMerging))
	{
		this->collisionCalc.lowerTolerance = *this->Gui.LowerLimit * *this->Gui.LowerLimit;
		this->collisionCalc.upperTolerance = *this->Gui.UpperLimit * *this->Gui.UpperLimit;
		findCollisions(&this->collisionCalc);
		this->collisionCalc.isDone = true;
	}

	Data * actData;

	if(!*Gui.DisplayColliding && !*Gui.DisplayMerging && !*Gui.DisplayInverted)
		//nothing selected = display all points
	{
		actData = &this->allData;
	}
	else if (!*Gui.DisplayColliding && !*Gui.DisplayMerging && *Gui.DisplayInverted)
		//dont display anything..
	{
		actData = &this->emptyData;
	}
	else
	{
		actData = &this->selectedData;
		generateSelection(&this->collisionCalc, &this->selection); // generates id maps
		generatePoints(&this->selection, actData);
		generateTracks(&this->selection, actData);
	}
	
	out->SetPoints(actData->Position);
	out->SetVerts(actData->Cells);
	out->SetLines(actData->Tracks);
	out->GetPointData()->AddArray(actData->Velocity);
	out->GetPointData()->AddArray(actData->GId);
	out->GetPointData()->AddArray(actData->SnapId);
	out->GetPointData()->AddArray(actData->TrackId);
	out->GetPointData()->AddArray(actData->Mvir);
	out->GetPointData()->AddArray(actData->Rvir);
	out->GetPointData()->AddArray(actData->Vmax);
	out->GetPointData()->AddArray(actData->Rmax);
	out->GetPointData()->AddArray(actData->Redshift);
	out->GetPointData()->AddArray(actData->Time);
	out->GetPointData()->AddArray(actData->Cv);
	out->GetPointData()->AddArray(actData->CollisionTypePoint);
	out->GetPointData()->AddArray(actData->CollisionTypeTrack);
	//out->GetPointData()->SetScalars(actData->Colors);


	//this->reset();
	
	/*
	if(*this->Gui.DisplaySelected)
	{
		this->dataInfo.selectedSnapshots.clear();
		this->dataInfo.nSelectedSnapshots = -1;
		this->dataInfo.selectedTracks.clear();
		this->dataInfo.nSelectedTracks = -1;

		if(*this->Gui.DisplaySelectedSnapshot)
		{
			this->dataInfo.selectedSnapshots.push_back(*this->Gui.DisplaySelectedSnapshotNr);
			this->dataInfo.nSelectedSnapshots = 1;
		}

		if(*this->Gui.DisplaySelectedTrack)
		{
			this->dataInfo.selectedTracks.push_back(*this->Gui.DisplaySelectedTrackNr);
			this->dataInfo.nSelectedTracks = 1;
		}
	}
	else if (*this->Gui.DisplayCalculated)
	{
		if (this->calcInfo.tolerance != *this->Gui.CalculationImpactParameter * *this->Gui.CalculationImpactParameter)
		{
			this->calcInfo.calcDone = false;
		}

		if(!this->calcInfo.calcDone)
		{
			doCalculations(*this->Gui.CalculationImpactParameter * *this->Gui.CalculationImpactParameter,3); //fills dataInfo.selected* fields
			this->calcInfo.calcDone = true;
			this->calcInfo.tolerance = *this->Gui.CalculationImpactParameter * *this->Gui.CalculationImpactParameter;
		}
		//else{vtkErrorMacro("Skipping recalculation...");}
		
		this->dataInfo.nSelectedPoints = this->calcInfo.nSelectedPoints;
		this->dataInfo.nSelectedTracks = this->calcInfo.nSelectedTracks;
		this->dataInfo.nSelectedSnapshots = -1;

		this->dataInfo.selectedPoints.clear();
		this->dataInfo.selectedPoints = this->calcInfo.selectedPoints;

		this->dataInfo.selectedTracks.clear();
		this->dataInfo.selectedTracks = this->calcInfo.selectedTracks;

		this->dataInfo.selectedSnapshots.clear();

	}
	else if (*this->Gui.DisplayEstimateTolerance)
	{

		//used this for approximation of parameter...
		this->calcTolerance();
		out->GetPointData()->AddArray(this->calcEstTol2);

	}*/

	/*
	Data * actData;

	if (*this->Gui.DisplayOnlySelectedData)
	{
		if (*this->Gui.DisplaySelected && !*this->Gui.DisplaySelectedSnapshot && !*this->Gui.DisplaySelectedTrack)
		{
			vtkErrorMacro("Please select at least one track or Snapshot to be displayed.\nI'm settig snapshot 0 to be displayed for you..");
			this->dataInfo.selectedSnapshots.clear();
			this->dataInfo.nSelectedSnapshots = 1;
			this->dataInfo.selectedSnapshots.push_back(0);
			this->dataInfo.selectedTracks.clear();
			this->dataInfo.nSelectedTracks = -1;
		}

		this->generateIdMap();
		this->generatePoints();
		this->generateTracks();
		this->generateColors();

		actData = &this->selectedData;

		//out->SetPoints(this->selectedData.Position);
		//out->SetVerts(this->selectedData.Cells);
		//out->SetLines(this->selectedData.Tracks);
		//out->GetPointData()->AddArray(this->selectedData.Velocity);
		//out->GetPointData()->AddArray(this->selectedData.GId);
		//out->GetPointData()->AddArray(this->selectedData.SnapId);
		//out->GetPointData()->AddArray(this->selectedData.Rvir);
		//out->GetPointData()->AddArray(this->selectedData.TrackId);
		//out->GetPointData()->SetScalars(this->selectedData.Colors);


	}
	else 
	{
		this->generateColors();

		actData = &this->allData;


		//out->SetPoints(this->allData.Position);
		//out->SetVerts(this->allData.Cells);
		//out->SetLines(this->allData.Tracks);
		//out->GetPointData()->AddArray(this->allData.Velocity);
		//out->GetPointData()->AddArray(this->allData.GId);
		//out->GetPointData()->AddArray(this->allData.SnapId);
		//out->GetPointData()->AddArray(this->allData.RVir);
		//out->GetPointData()->AddArray(this->allData.TrackId);
		//out->GetPointData()->SetScalars(this->allData.Colors);
	}
	*/

	// trying with luts, but ddoesnt work as expected...
	/*
	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetValueRange(0.5,1);
	lut->SetHueRange(0.9,1.0);
	lut->Build();
	double r [] = {0.0, 9.0};
	lut->SetTableRange(r);
	lut->Modified();
	//lut->SetNumberOfTableValues(9);
	//lut->Build();
	//lut->SetTableValue(2,1,0,0,1);
	//lut->SetNumberOfColors(24);
	//lut->SetHueRange(0,1.5);
	//lut->SetSaturationRange(1,1);
	//lut->SetValueRange(0,3.14);
	//lut->SetScale(24);
	//lut->Build();

	//double col[3];
	//lut->GetColor(2,col);
	//vtkErrorMacro(":" << col[1]);

	actData->Rvir->SetLookupTable(lut);

	vtkSmartPointer<vtkLookupTable> lut2 = vtkSmartPointer<vtkLookupTable>::New();
	//lut2->SetTableRange(0.5,1);
	//lut2->SetSaturationRange(1,1);
	//lut2->SetHueRange(0.5,0.75);
	//lut2->SetValueRange(0,3.14);
	//lut2->SetAlphaRange(1,1);
	//lut2->SetNumberOfTableValues(60);

	//lut2->Build();
	lut2->SetTableRange(0.0,3.0);
	double * range;
	range = lut2->GetRange();
	*(range) = 0.00;
	*(range+1) = 3.0;
	//lut2->SetNumberOfTableValues(4);
	//lut2->SetTableValue(0,1,0,0);
	//lut2->SetTableValue(1,1,1,0);
	//lut2->SetTableValue(2,1,1,1);
	//lut2->SetTableValue(3,0,0,1);
	lut2->Build();
	lut2->ForceBuild();

	actData->SnapId->SetLookupTable(lut);
	//out->GetCellData()->get


	//lut2->Print(

	//lut2->PrintSelf(cerr, vtkIndent(0)) ;
	*/



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
		this->dataInfo.hubble = sqlite3_column_int(res, 7);
		//this->dataInfo.nSnapshots = sqlite3_column_int(res, 9);
		// TODO here more data can be read in..
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
		this->SnapInfo		stores information about the snapshots
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots()
{
// prepare the variables

/*
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
*/

	//this->SnapInfo.clear();

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
			/*if (this->allData.Position->GetData()->GetComponent(offset,0) == 0 && 
				this->allData.Position->GetData()->GetComponent(offset,1) == 0 &&
				this->allData.Position->GetData()->GetComponent(offset,2) == 0)
			{
				continue;
			}*/

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
	vtkErrorMacro("track count: " << this->dataInfo.nTracks);
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
	sql_query = "SELECT * FROM snapinfo ORDER BY snap_id, redshift, time";

	//vtkErrorMacro("SQL query: " + sql_query);

	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);

	if (sql_error != SQLITE_OK)
	{
		vtkErrorMacro("Error with sql query! Error: " << sql_error);
		return 0;
	}

	int snap_id;
	int i = 0;
	SnapshotInfo tmpSnap;

	this->SnapInfo.clear();

	while (sqlite3_step(res) == SQLITE_ROW)
	{
		snap_id = sqlite3_column_int(res, 0);
		
		// actually there should be numbers in the snap_id column, but
		// in the converted vl2 data there aren't any, so just assume its
		// ordered... check this here:
		// TODO really check for empty field, because id can be 0 for first snapshot...
		if (snap_id<1)
		{
			snap_id = i;
		}
		tmpSnap.snapshotNr = sqlite3_column_double(res, 1);
		tmpSnap.redshift = sqlite3_column_double(res, 2);
		tmpSnap.time = sqlite3_column_double(res, 3);
		tmpSnap.npart = sqlite3_column_int(res, 5);
		
		this->SnapInfo.push_back(tmpSnap);
		i++;
	}
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
int vtkSQLiteReader::calculateAdditionalData()
{

	// calc average cv in one snapshot
	for (int i = 0; i<this->dataInfo.nSnapshots;i++)
	{
		int length = this->SnapInfo.at(i).lenght;
		int offset = this->SnapInfo.at(i).Offset;
		int count = 0;
		double sum = 0;
		double ave = 0;
		for (j=0; j<length;j++)
		{
			sum =+ this->allData.Cv->GetTuple1(j+offset);
			count ++;
		}
		ave = sum / (double)count;
		//TODO insert data into points


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
{/*
	unsigned char red, green, blue;
	int snapid, pid_old, pid_new;

	int nPoints;
	Data * data;

	if (*this->Gui.DisplayOnlySelectedData)
	{
		nPoints = this->dataInfo.nAllSelectedPoints;
		data = &this->selectedData;
	}
	else
	{
		nPoints = this->dataInfo.nPoints;
		data = &this->allData;
	}

	// init color vector
	data->Colors->SetNumberOfTuples(nPoints);

	//basic coloring of all points, according to snapid.
	for (int i = 0; i<nPoints; i++)
	{
		snapid = data->SnapId->GetTuple1(i);

		blue = 190*(float)(snapid) / (float)(this->dataInfo.nSnapshots-1);
		green = 0;
		red = 0;

		data->Colors->InsertTuple3(i, red, green, blue);
	}

	//highlightning
	if(! (*this->Gui.DisplayOnlySelectedData))
	{
		int length, pid, trackid, snapid, offset;
		
		// highlight tracks
		//COLORSETTING
		red = 255;
		green = 255;
		blue = 0;

		for (int j = 0; j<this->dataInfo.nSelectedTracks; j++)
		{
			trackid = this->dataInfo.selectedTracks.at(j);
			
			//overflowprotection
			if (trackid>=this->dataInfo.nTracks){trackid=this->dataInfo.nTracks-1;}

			length = this->TracksInfo.at(trackid).nPoints;
			for (int i = 0; i<length; i++)
			{
				pid = this->TracksInfo.at(trackid).PointsIds.at(i);
				data->Colors->InsertTuple3(pid,red,green,blue);
			}
		}

		// highlight snapshots
		//COLORSETTING
		red = 255;
		green = 0;
		blue = 0;

		for (int j = 0; j<this->dataInfo.nSelectedSnapshots; j++)
		{
			snapid = this->dataInfo.selectedSnapshots.at(j);

			//overflowprotection
			if (snapid>=this->dataInfo.nSnapshots){snapid=this->dataInfo.nSnapshots-1;}

			offset = this->SnapInfo.at(snapid).Offset;
			length = this->SnapInfo.at(snapid).lenght;
			
			for (int i = offset; i<offset+length; i++)
			{
				data->Colors->InsertTuple3(i,red,green,blue);
			}
		}

		// highlight points
		//COLORSETTING
		red = 255;
		green = 0;
		blue = 255;

		for (int j = 0; j<this->dataInfo.nSelectedPoints; j++)
		{
			pid = this->dataInfo.selectedPoints.at(j);
			data->Colors->InsertTuple3(pid,red,green,blue);
		}

	}

	if (*this->Gui.DisplayOnlySelectedData && ! (*this->Gui.DisplaySelected))
	{
		if(*this->Gui.calcHighlightCollisionPoints)
		{
			;
		}
		if (*this->Gui.calcHighlightTrack)
		{
			//COLORSETTING
			red = 0;
			green = 255;
			blue = 0;

			int trackid;

			if (*this->Gui.calcHighlightTrackNr < this->calcInfo.nSelectedTracks)
			{
				trackid = this->dataInfo.selectedTracks.at(*this->Gui.calcHighlightTrackNr);
			} 
			else
			{
				trackid = this->dataInfo.selectedTracks.at(this->calcInfo.nSelectedTracks - 1);
			}

			for (int i = 0; i<this->TracksInfo.at(trackid).nPoints; i++)
			{
				pid_old = this->TracksInfo.at(trackid).PointsIds.at(i);
				pid_new = this->dataInfo.idMap1.at(pid_old);
				data->Colors->InsertTuple3(pid_new,red,green,blue);
			}

		}
		if (*this->Gui.calcHighlightSnapshot)
		{
			;
		}
		if (*this->Gui.calcHighlightPoint)
		{
			;
		}
	}

*/
	return 1;
}

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
int vtkSQLiteReader::calcTolerance()
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
double vtkSQLiteReader::getDistance2(int gid1, int gid2){

	double x1[3];
	double x2[3];

	this->allData.Position->GetPoint(gid1, x1);
	this->allData.Position->GetPoint(gid2, x2);

	return (x1[0]-x2[0])*(x1[0]-x2[0])+
		(x1[1]-x2[1])*(x1[1]-x2[1])+
		(x1[2]-x2[2])*(x1[2]-x2[2]);
}

/*----------------------------------------------------------------------------
generates the id map from alldata ids to only selected points ids
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generateIdMap()
{
	/*
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

	*/
	return 1;
}

/*----------------------------------------------------------------------------
generates points in the output vector
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generatePoints(SelectionStruct* selection, Data* actData)
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
	actData->CollisionTypePoint->SetNumberOfTuples(selection->nSelectedPoints);
	actData->CollisionTypeTrack->SetNumberOfTuples(selection->nSelectedPoints);

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
		actData->CollisionTypePoint->InsertTuple(i,
			this->allData.CollisionTypePoint->GetTuple(idMap2->at(i)));
		actData->CollisionTypeTrack->InsertTuple(i,
			this->allData.CollisionTypeTrack->GetTuple(idMap2->at(i)));
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
}



/*----------------------------------------------------------------------------
generates tracks in the output vector
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generateTracks(SelectionStruct* selection, Data* actData)
{
	int trackId, oldId, newId;
	vtkSmartPointer<vtkPolyLine> nextLine;

	actData->Tracks->Reset();

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
}

/*----------------------------------------------------------------------------
resets all settings (not data)
	assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::reset()
{
	/*
	this->dataInfo.nSelectedPoints = -1;
	this->dataInfo.selectedPoints.clear();

	this->dataInfo.nSelectedSnapshots = -1;
	this->dataInfo.selectedSnapshots.clear();

	this->dataInfo.nSelectedTracks = -1;
	this->dataInfo.selectedTracks.clear();

	this->dataInfo.nAllSelectedPoints = -1;

	this->dataInfo.idMap1.clear();
	this->dataInfo.idMap2.clear();
	*/

	/*this->selectedData.Position->Delete();
	this->selectedData.Velocity->Delete();
	this->selectedData.Cells->Delete();
	this->selectedData.Tracks->Delete();
	this->selectedData.TrackId->Delete();
	this->selectedData.GId->Delete();
	this->selectedData.SnapId->Delete();
	this->selectedData.RVir->Delete();
	this->selectedData.Colors->Delete();
*/
	/*
	this->selectedData.Position = vtkSmartPointer<vtkPoints>::New();
	this->selectedData.Velocity = vtkSmartPointer<vtkFloatArray>::New();
	this->selectedData.Cells = vtkSmartPointer<vtkCellArray>::New();
	this->selectedData.Tracks = vtkSmartPointer<vtkCellArray>::New();
	this->selectedData.TrackId = vtkSmartPointer<vtkIntArray>::New();
	this->selectedData.GId = vtkSmartPointer<vtkIntArray>::New();
	this->selectedData.SnapId = vtkSmartPointer<vtkIntArray>::New();
	this->selectedData.Rvir = vtkSmartPointer<vtkFloatArray>::New();
	this->selectedData.Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();


	this->selectedData.Velocity->SetNumberOfComponents(3);
	this->selectedData.TrackId->SetNumberOfComponents(1);
	this->selectedData.GId->SetNumberOfComponents(1);
	this->selectedData.SnapId->SetNumberOfComponents(1);
	this->selectedData.Rvir->SetNumberOfComponents(1);
	this->selectedData.Colors->SetNumberOfComponents(3);

	this->selectedData.Velocity->SetName("Velocity");
	this->selectedData.TrackId->SetName("TrackId");
	this->selectedData.GId->SetName("GId");
	this->selectedData.SnapId->SetName("SnapId");
	this->selectedData.Rvir->SetName("RVir");
	this->selectedData.Colors->SetName("Colors");
*/
	return 1;
}


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
int vtkSQLiteReader::findCollisions(CollisionCalculationStruct* pCollCalc){

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
}


/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generateSelection(CollisionCalculationStruct * collCalc, SelectionStruct* select)
{
	std::vector<int> * idMap1 = &select->vPointIdMap;
	std::vector<int> * idMap2 = &select->vPointIdMapReverse;
	std::vector<int> * trackIdMap = &select->vTrackIdMap;

	int p_counter = 0;
	int t_counter = 0;
	int trackCount = 0;
	int offset = 0;
	int length = 0;
	int pid = 0;
	int trackid = 0;

	std::vector<int> * v1, * v2, * v3;
	std::vector<int>::iterator it1, it2;
	v1 = &collCalc->Colliding.vTrackIds;
	v2 = &collCalc->Merging.vTrackIds;
	v3 = &select->vTrackIds;


	idMap1->clear();
	idMap2->clear();
	trackIdMap->clear();
		
	idMap1->resize(this->dataInfo.nPoints,-1);
	idMap2->resize(this->dataInfo.nPoints,-1);
	trackIdMap->resize(this->dataInfo.nTracks,-1);


	//generate forwardmaps
	if(*this->Gui.DisplayColliding)
	{
		fillIdList(idMap1,&p_counter,trackIdMap,&t_counter,&collCalc->Colliding,1);
	}
	
	if(*this->Gui.DisplayMerging)
	{
		fillIdList(idMap1,&p_counter,trackIdMap,&t_counter,&collCalc->Merging,2);
	}
	
	if(*this->Gui.DisplayMerging && *this->Gui.DisplayColliding)
	{
		//generate track ids (merge two vectors, expect sorted!)
		v3->clear();

		it1 = v1->begin();
		it2 = v2->begin();
		bool goon = true;

		while (goon)
		{
			if (it1 == v1->end()){break;}
			if (it2 == v2->end()){it2 = it1; v2 = v1; break;}

			if(*it1 < *it2)	{v3->push_back(*it1); it1++; continue;}
			if(*it1 > *it2)	{v3->push_back(*it2); it2++; continue;}
			if(*it1 ==*it2)	{v3->push_back(*it1); it1++; it2++; continue;}
		}
		while (it2<v2->end()){v3->push_back(*it2); it2++;}
	}
	else if(*this->Gui.DisplayMerging)
	{
		v3->assign(v2->begin(),v2->end());
	}
	else if(*this->Gui.DisplayColliding)
	{
		v3->assign(v1->begin(),v1->end());
	}
	
	if(*this->Gui.DisplayInverted)
	{
		//invertIdList(idMap1,&counter);
		p_counter = 0;

		for (int i = 0; i<idMap1->size();i++)
		{
			if (idMap1->at(i) == -1)
			{
				idMap1->at(i) = p_counter;
				p_counter++;
			}
			else
			{
				idMap1->at(i) = -1;
			}
		}
	
		//invert tracks
		t_counter = 0;
		v3->clear();
		for (int i = 0; i<trackIdMap->size();i++)
		{
			if (trackIdMap->at(i) == -1)
			{
				trackIdMap->at(i) = t_counter;
				v3->push_back(i);
				t_counter++;
			}
			else
			{
				trackIdMap->at(i) = -1;
			}
		}

	}

	//generate reverse point map
	p_counter = 0;
	for(int i = 0; i<idMap1->size();i++)
	{
		pid=idMap1->at(i);
		if(pid>-1)
		{
			idMap2->at(pid) = i;
			p_counter++;
		}
	}
	idMap2->resize(p_counter);



	select->nSelectedPoints = p_counter;
	select->nSelectedTracks = t_counter;
	select->nSelectedSnapshots = -1;

	return 1;

/*
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
*/


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
int vtkSQLiteReader::fillIdList(std::vector<int> * p_map, int * p_count,
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


	//fill in points
	/*
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




	return 1;
}