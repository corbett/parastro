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

//#include "vtkSmartPointer.h"
//#include "tipsylib/ftipsy.hpp" // functions take tipsy particle objects
//#include <vtkstd/vector>
//#include "vtkFloatArray.h"

class VTK_EXPORT vtkSQLiteReader : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader* New();
	vtkTypeRevisionMacro(vtkSQLiteReader,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	// set/get the db query
	//vtkSetStringMacro(SqlQuery);
	//vtkGetStringMacro(SqlQuery);

	char*	GetSqlQuery();
	void	SetSqlQuery(const char* name);

protected:
	vtkSQLiteReader();
	~vtkSQLiteReader();
	char* FileName;
	char* SqlQuery;

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

private:
  vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
  void operator=(const vtkSQLiteReader&);  // Not implemented.
};

#endif