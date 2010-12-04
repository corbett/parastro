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
		double redshift;
		double time;
		int npart;
	};

	struct Track {
		int nPoints;
		//int StartSnap;
		//int EndSnap;
		std::vector<vtkIdType> PointsIds; //v of length nSnaps, -1 if this track doesnt have a point in this snap
	};

	struct DataInformation {
		int nPoints;
		int nTracks;
		int nSnapshots;
		int nSelectedPoints; //-1 means all, single points selected to display
		int nSelectedTracks; // -1 means all
		int nSelectedSnapshots; 
		int nAllSelectedPoints; // # of all points (inkl tracks, snaps) to display
		std::vector<int>	selectedPoints; // holds globalids from points to display
		std::vector<int>	selectedSnapshots; //holds snapids of points from snapshot to display
		std::vector<int>	selectedTracks; //holds trackid from points to display
		std::vector<int>	idMap1; //mapps ids from read in data to displayed data
		std::vector<int>	idMap2; //mapps ids from displayed data to read in data
	};

	struct GUISettings {
		bool * DisplayOnlySelectedData;
		bool * DisplaySelected;
		bool * DisplaySelectedSnapshot;
		int * DisplaySelectedSnapshotNr;
		bool * DisplaySelectedTrack;
		int * DisplaySelectedTrackNr;
		bool * DisplayCalculated;
		double * CalculationImpactParameter;
	};

	struct Data {
		vtkSmartPointer<vtkPoints>			Position;
		vtkSmartPointer<vtkFloatArray>		Velocity;
		vtkSmartPointer<vtkCellArray>		Cells;
		vtkSmartPointer<vtkCellArray>		Tracks;
		vtkSmartPointer<vtkIdTypeArray>		TrackId;
		vtkSmartPointer<vtkIdTypeArray>		GId;
		vtkSmartPointer<vtkIdTypeArray>		SnapId;
		vtkSmartPointer<vtkFloatArray>		RVir;
		vtkSmartPointer<vtkUnsignedCharArray> Colors;
		//vtkSmartPointer<vtkUnsignedCharArray> opacity;
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
	GUISettings Gui;
	Data allData, selectedData;
	
	std::vector<SnapshotInfo> SnapInfo;
	std::vector<Track> TracksInfo;

	//gui variables
	char* FileName;

	bool DisplayOnlySelectedData;

	bool DisplaySelected;
		bool DisplaySelectedSnapshot;
		int DisplaySelectedSnapshotNr;
		bool DisplaySelectedTrack;
		int DisplaySelectedTrackNr;

	bool DisplayCalculated;
		double CalculationImpactParameter;


private:
	vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
	void operator=(const vtkSQLiteReader&);  // Not implemented.


// used in v3
	//variables
	sqlite3 * db;
	bool dataIsRead;

	//functions
	int vtkSQLiteReader::ReadHeader(); // reads the database header information

	int vtkSQLiteReader::readSnapshots(); // reads the snapshots
	int vtkSQLiteReader::readSnapshotInfo(); 
	int vtkSQLiteReader::readTracks();
	int vtkSQLiteReader::selectPoints();
	int vtkSQLiteReader::generateColors();

	int vtkSQLiteReader::doCalculations(double,int);

	int vtkSQLiteReader::generateIdMap();
	int vtkSQLiteReader::generatePoints();


	// helper
	int openDB(char*);
	vtkStdString vtkSQLiteReader::Int2Str(int);
	double distance(int, int);

	// old stuff - not yet, or not anymore needed
	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);


//ETX
};

#endif