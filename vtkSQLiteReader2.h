/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader2.h,v $

  Copyright (c) Rafael Kueng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSQLiteReader2 - Read points from a SQLite database
// .SECTION Description
// Reads points from a SQLite DB and displays them on screen


#ifndef __vtkSQLiteReader2_h
#define __vtkSQLiteReader2_h

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


class VTK_EXPORT vtkSQLiteReader2 : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader2* New();
	vtkTypeRevisionMacro(vtkSQLiteReader2,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

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



//BTX
protected:

	//functions

	vtkSQLiteReader2();
	~vtkSQLiteReader2();

	virtual int FillOutputPortInformation(int vtkNotUsed(port),	vtkInformation* info);


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
	vtkSQLiteReader2(const vtkSQLiteReader2&);  // Not implemented.
	void operator=(const vtkSQLiteReader2&);  // Not implemented.

	//variables
	sqlite3 * db;
	bool dataIsRead;

	//functions
	int vtkSQLiteReader2::ReadHeader(); // reads the database header information

	int vtkSQLiteReader2::readSnapshots(); // reads the snapshots
	int vtkSQLiteReader2::readSnapshotInfo(); 
	int vtkSQLiteReader2::readTracks();
	int vtkSQLiteReader2::calculateAdditionalData();
	int vtkSQLiteReader2::generateColors();

	int vtkSQLiteReader2::findCollisions(CollisionCalculationStruct*);
	//int vtkSQLiteReader2::doCalculations(double,int);
	int vtkSQLiteReader2::calcTolerance();

	int vtkSQLiteReader2::generateSelection(CollisionCalculationStruct*, SelectionStruct*);
	int vtkSQLiteReader2::fillIdList(std::vector<int>*, int*,
		std::vector<int>*, int*,
		CollisionResultStruct*, int);
	int vtkSQLiteReader2::generateIdMap();
	int vtkSQLiteReader2::generatePoints(SelectionStruct*, Data*);
	int vtkSQLiteReader2::generateTracks(SelectionStruct*, Data*);
	int vtkSQLiteReader2::reset();

	// helper
	int openDB(char*);
	vtkStdString vtkSQLiteReader2::Int2Str(int);
	double getDistance2(int, int);
	double getDistanceToO(int);

	// old stuff - not yet, or not anymore needed
	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader2::SQLQuery(vtkStdString, sqlite3_stmt*);


//ETX
};

#endif