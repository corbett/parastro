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

	vtkSetMacro(HighlightSnapshot,int);
	vtkGetMacro(HighlightSnapshot,int);

	vtkSetMacro(HighlightTrack,int);
	vtkGetMacro(HighlightTrack,int);

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

	//variables
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkCellArray>		Tracks;
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkIdTypeArray>		GId;
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		RVir;
	vtkSmartPointer<vtkUnsignedCharArray> colors;

	int nParticles3;
	int nTracks3;
	int nSnapshots;

	std::vector<SnapshotInfo> SnapInfo3;

	//gui variables
	char* FileName;
	int HighlightSnapshot;
	int HighlightTrack;


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
	int vtkSQLiteReader::generateColors();

	int openDB(char*);
	vtkStdString vtkSQLiteReader::Int2Str(int);

// old stuff - not yet, or not anymore needed

	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);


//ETX
};

#endif