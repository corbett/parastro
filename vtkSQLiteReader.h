/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.h,v $

  Copyright (c) Rafael Kueng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSQLiteReader - Read points from a SQLite database
// .SECTION Description
// Reads points from a SQLite DB and displays them on screen


#ifndef __vtkSQLiteReader_h
#define __vtkSQLiteReader_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "sqlitelib/sqlite3.h" // sqlite headerfile

#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkstd/vector>

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;


class VTK_EXPORT vtkSQLiteReader : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader* New();
	vtkTypeRevisionMacro(vtkSQLiteReader,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	vtkSetMacro(DisplayColliding,bool);
	vtkGetMacro(DisplayColliding,bool);

	vtkSetMacro(DisplayMerging,bool);
	vtkGetMacro(DisplayMerging,bool);

	vtkSetMacro(DisplayInverted,bool);
	vtkGetMacro(DisplayInverted,bool);

	vtkSetMacro(UpperLimit,double);
	vtkGetMacro(UpperLimit,double);

	vtkSetMacro(LowerLimit,double);
	vtkGetMacro(LowerLimit,double);


/*
	vtkSetMacro(DisplayOnlySelectedData,bool);
	vtkGetMacro(DisplayOnlySelectedData,bool);

	vtkSetMacro(DisplaySelected,bool);
	vtkGetMacro(DisplaySelected,bool);

		vtkSetMacro(DisplaySelectedSnapshot,bool);
		vtkGetMacro(DisplaySelectedSnapshot,bool);
		vtkSetMacro(DisplaySelectedSnapshotNr,int);
		vtkGetMacro(DisplaySelectedSnapshotNr,int);
		vtkSetMacro(DisplaySelectedTrack,bool);
		vtkGetMacro(DisplaySelectedTrack,bool);
		vtkSetMacro(DisplaySelectedTrackNr,int);
		vtkGetMacro(DisplaySelectedTrackNr,int);


	vtkSetMacro(DisplayCalculated,bool);
	vtkGetMacro(DisplayCalculated,bool);

		vtkSetMacro(CalculationImpactParameter,double);
		vtkGetMacro(CalculationImpactParameter,double);

		vtkSetMacro(calcHighlightCollisionPoints,bool);
		vtkGetMacro(calcHighlightCollisionPoints,bool);

		vtkSetMacro(calcHighlightTrack,bool);
		vtkGetMacro(calcHighlightTrack,bool);
		vtkSetMacro(calcHighlightSnapshot,bool);
		vtkGetMacro(calcHighlightSnapshot,bool);
		vtkSetMacro(calcHighlightPoint,bool);
		vtkGetMacro(calcHighlightPoint,bool);
		vtkSetMacro(calcHighlightTrackNr,int);
		vtkGetMacro(calcHighlightTrackNr,int);
		vtkSetMacro(calcHighlightSnapshotNr,int);
		vtkGetMacro(calcHighlightSnapshotNr,int);
		vtkSetMacro(calcHighlightPointNr,int);
		vtkGetMacro(calcHighlightPointNr,int);


	vtkSetMacro(DisplayEstimateTolerance,bool);
	vtkGetMacro(DisplayEstimateTolerance,bool);
	*/



//BTX
protected:

	//functions

	vtkSQLiteReader();
	~vtkSQLiteReader();

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);
	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

	//structs
	struct SnapshotInfo {
		vtkIdType Offset; //stores the id where this snapshot starts
		int lenght; // stores the amount of halos in this snapshot
		std::vector<int> PointId; //allocates gid to new id (PointId.at(gid) = id+offset) maybee save here some mem, by using smaller datatype..
		
		int snapshotNr; // snapshot Nr from simulation (this is not the id!)
		double redshift;
		double time;
		int npart;
		double CvAverage;
	};

	struct Track {
		int nPoints;
		//int StartSnap;
		//int EndSnap;
		std::vector<vtkIdType> PointsIds; //v of length nSnaps, entry of -1 if this track doesnt have a point in this snap
	};

	struct DataInformation {
		bool dataIsRead;
		int nPoints;
		int nTracks;
		int nSnapshots;
		double hubble;
		/*
		int nSelectedPoints; //-1 means all, single points selected to display
		int nSelectedTracks; // -1 means all
		int nSelectedSnapshots; 
		int nAllSelectedPoints; // # of all points (inkl tracks, snaps) to display
		std::vector<int>	selectedPoints; // holds globalids from points to display
		std::vector<int>	selectedSnapshots; //holds snapids of points from snapshot to display
		std::vector<int>	selectedTracks; //holds trackid from points to display
		std::vector<int>	idMap1; //mapps ids from read in data to displayed data
		std::vector<int>	idMap2; //mapps ids from displayed data to read in data
		*/
	};

	struct GuiStruct {
		bool * DisplayColliding;
		bool * DisplayMerging;
		bool * DisplayInverted;

		double * LowerLimit;
		double * UpperLimit;

		/*

		bool * DisplayOnlySelectedData;
		bool * DisplaySelected;
		bool * DisplaySelectedSnapshot;
		int * DisplaySelectedSnapshotNr;
		bool * DisplaySelectedTrack;
		int * DisplaySelectedTrackNr;
		bool * DisplayCalculated;
		double * CalculationImpactParameter;

		bool * calcHighlightCollisionPoints;
		bool * calcHighlightTrack;
		bool * calcHighlightSnapshot;
		bool * calcHighlightPoint;
		int * calcHighlightTrackNr;
		int * calcHighlightSnapshotNr;
		int * calcHighlightPointNr;

		bool * DisplayEstimateTolerance;
		*/
	};

	struct Data {
		vtkSmartPointer<vtkPoints>			Position;
		vtkSmartPointer<vtkFloatArray>		Velocity;
		vtkSmartPointer<vtkCellArray>		Cells;
		vtkSmartPointer<vtkCellArray>		Tracks;
		vtkSmartPointer<vtkIntArray>		TrackId;
		vtkSmartPointer<vtkIntArray>		GId;
		vtkSmartPointer<vtkIntArray>		SnapId;
		vtkSmartPointer<vtkFloatArray>		Mvir;
		vtkSmartPointer<vtkFloatArray>		Rvir;
		vtkSmartPointer<vtkFloatArray>		Vmax;
		vtkSmartPointer<vtkFloatArray>		Rmax;
		vtkSmartPointer<vtkFloatArray>		Redshift;
		vtkSmartPointer<vtkFloatArray>		Time;
		vtkSmartPointer<vtkFloatArray>		Cv;
		vtkSmartPointer<vtkFloatArray>		CvAverage;
		vtkSmartPointer<vtkFloatArray>		RGc;
		vtkSmartPointer<vtkUnsignedCharArray> CollisionTypePoint;
		vtkSmartPointer<vtkUnsignedCharArray> CollisionTypeTrack;
	};

	struct CalculationSettings {
		// OLD dont use this
		bool calcDone;
		double tolerance;
		int nCollisions;
		int nSelectedPoints; //-1 means all, single points selected to display
		int nSelectedTracks; // -1 means all
		int nAllSelectedPoints; // # of all points (inkl tracks, snaps) to display
		std::vector<int>	selectedPoints; // stores globalids from points to display
		std::vector<int>	selectedTracks; //stores trackid from points to display
	};

	struct CollisionResultStruct{
		int nEvents;
		int nPoints;
		std::vector<int> vPointIds;
		int nTracks;
		std::vector<int> vTrackIds;
	};

	struct CollisionCalculationStruct{
		bool isDone;
		double upperTolerance;
		double lowerTolerance;
		CollisionResultStruct Colliding;
		CollisionResultStruct Merging;
	};

	struct SelectionStruct{
		int nSelectedPoints;
		int nSelectedTracks;
		int nSelectedSnapshots;
		std::vector<int> vPointIds;
		std::vector<int> vTrackIds;
		std::vector<int> vSnapshotIds;
		std::vector<int> vPointIdMap;
		std::vector<int> vPointIdMapReverse; // equal to vPointIds!!
		std::vector<int> vTrackIdMap;

	};

	struct ResultOfEsimationOfTolerance {
		int goal;
		double parameter;
		int nCollisions;
	};

	//variables (not used!)
	/*
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkCellArray>		Tracks;
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkIdTypeArray>		GId;
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		RVir;
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	*/
	//vtkSmartPointer<vtkUnsignedCharArray> opacity;

	DataInformation dataInfo;
	GuiStruct Gui;
	Data allData, selectedData, emptyData;
	CalculationSettings calcInfo;
	SelectionStruct	selection;
	
	CollisionCalculationStruct collisionCalc;

	std::vector<ResultOfEsimationOfTolerance> calcEstTol;
	vtkSmartPointer<vtkFloatArray> calcEstTol2;
	
	std::vector<SnapshotInfo> SnapInfo;
	std::vector<Track> TracksInfo;

	//gui variables
	char* FileName;

	bool DisplayColliding;
	bool DisplayMerging;
	bool DisplayInverted;
	double LowerLimit;
	double UpperLimit;
	
	/*
	bool DisplayOnlySelectedData;

	bool DisplaySelected;
		bool DisplaySelectedSnapshot;
		int DisplaySelectedSnapshotNr;
		bool DisplaySelectedTrack;
		int DisplaySelectedTrackNr;

	bool DisplayCalculated;
		double CalculationImpactParameter;

		bool calcHighlightCollisionPoints;
		bool calcHighlightTrack;
		bool calcHighlightSnapshot;
		bool calcHighlightPoint;
		int calcHighlightTrackNr;
		int calcHighlightSnapshotNr;
		int calcHighlightPointNr;


	bool DisplayEstimateTolerance;
	*/


private:
	vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
	void operator=(const vtkSQLiteReader&);  // Not implemented.

	//variables
	sqlite3 * db;
	bool dataIsRead;

	//functions
	int vtkSQLiteReader::ReadHeader(); // reads the database header information

	int vtkSQLiteReader::readSnapshots(); // reads the snapshots
	int vtkSQLiteReader::readSnapshotInfo(); 
	int vtkSQLiteReader::readTracks();
	int vtkSQLiteReader::calculateAdditionalData();
	int vtkSQLiteReader::generateColors();

	int vtkSQLiteReader::findCollisions(CollisionCalculationStruct*);
	//int vtkSQLiteReader::doCalculations(double,int);
	int vtkSQLiteReader::calcTolerance();

	int vtkSQLiteReader::generateSelection(CollisionCalculationStruct*, SelectionStruct*);
	int vtkSQLiteReader::fillIdList(std::vector<int>*, int*,
		std::vector<int>*, int*,
		CollisionResultStruct*, int);
	int vtkSQLiteReader::generateIdMap();
	int vtkSQLiteReader::generatePoints(SelectionStruct*, Data*);
	int vtkSQLiteReader::generateTracks(SelectionStruct*, Data*);
	int vtkSQLiteReader::reset();

	// helper
	int openDB(char*);
	vtkStdString vtkSQLiteReader::Int2Str(int);
	double getDistance2(int, int);
	double getDistanceToO(int);

	// old stuff - not yet, or not anymore needed
	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);


//ETX
};

#endif