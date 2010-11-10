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

// probably there better ways than list...
//TODO check this...
//#include <list>
//using namespace std;

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

//BTX
//struct halo {
//	int id;
//	vtkPoints coordinates;
//};
//ETX

class VTK_EXPORT vtkSQLiteReader : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader* New();
	vtkTypeRevisionMacro(vtkSQLiteReader,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	vtkSetMacro(DisplaySnapshot,int);
	vtkGetMacro(DisplaySnapshot,int);

//BTX
protected:
	vtkSQLiteReader();
	~vtkSQLiteReader();
	char* FileName;

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

	

// functions

// structs
	struct velocity {
		double vx;
		double vy;
		double vz;
	};

	struct snapshot {
		vtkSmartPointer<vtkPoints> coord;
		vtkSmartPointer<vtkCellArray> cells;
		std::vector<velocity> velo;
		// std::vector<int> npart;
		// std::vector<int> Mvir;
		// std::vector<int> Rvir;
		// ...
	};

	struct snapinfo {
		int snap_id;
		double redshift;
		double time;
		int npart;
	};

	struct trackPoint {
		int snap_id;
		int qid;
	};

	struct track {
		std::vector<trackPoint> point;
		int noOfPoints;
		vtkSmartPointer<vtkPolyLine> line;
	};

// variables
	int numSnaps;
	int numTracks;
	int totNumPoints;
	int DisplaySnapshot;
	std::vector<vtkSmartPointer<vtkPolyData>> data;
	
	std::vector<snapshot> data2;
	std::vector<snapinfo> snapinfoVector;
	std::vector<track> trackVector;

// --- v3 ------
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkCellArray>		Tracks;
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkIdTypeArray>		Qid;
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		RVir;
	vtkSmartPointer<vtkUnsignedCharArray> colors;

	int nParticles3;
	int nTracks3;

	struct SnapshotInfo3 {
		vtkIdType Offset; //stores the id where this snapshot starts
		int lenght; // stores the amount of halos in this snapshot
		double redshift;
		double time;
		int npart;
	};

	std::vector<SnapshotInfo3> SnapInfo3;

/* not used
	// for halos
	vtkIdType							ParticleIndex;
	vtkSmartPointer<vtkIdTypeArray>		ParticleId;
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkFloatArray>		nParticles;

	vtkSmartPointer<vtkFloatArray>		mVir;
	vtkSmartPointer<vtkFloatArray>		rVir;
	vtkSmartPointer<vtkFloatArray>		RHO;

	// for snapshots
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		Redshift;
	vtkSmartPointer<vtkFloatArray>		Time;
	vtkSmartPointer<vtkDataArray>		Halos;

	// for tracks
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkPolyData>		Lines;
*/


// constants


private:
	vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
	void operator=(const vtkSQLiteReader&);  // Not implemented.

	//functions
	int openDB(char*);

	// used in v3
	int vtkSQLiteReader::readSnapshots3();
	int vtkSQLiteReader::readSnapshotInfo3();
	int vtkSQLiteReader::readTracks3();
	int vtkSQLiteReader::generateColors();

	// old stuff
	int vtkSQLiteReader::readSnapshots(
		std::vector<vtkSmartPointer<vtkPolyData>> *);
	int vtkSQLiteReader::readSnapshots2();
	int vtkSQLiteReader::ReadHeader(vtkInformationVector*);
	int vtkSQLiteReader::ReadTracks();
	int vtkSQLiteReader::GenerateTracks();
	int vtkSQLiteReader::CollectLines(vtkSmartPointer<vtkCellArray>*);
	int vtkSQLiteReader::GenerateOutput(vtkPolyData *);

	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader::ReadSnapshotInfo();

	// helpers
	vtkStdString vtkSQLiteReader::Int2Str(int);
	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);

	//variables
	sqlite3 * db;
	bool dataIsRead;

	//constants


//ETX
};

#endif