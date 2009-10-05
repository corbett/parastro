/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib
 -christine
=========================================================================*/
// .NAME vtkGadgetIIReader - Read points from a Tipsy standard binary file
// .SECTION Description
// here is the desciprtion
#ifndef __vtkGadgetIIReader_h
#define __vtkGadgetIIReader_h
#include "vtkPolyDataAlgorithm.h" // needed as this class extends vtkPolyDataAlgorithm
#include "vtkSmartPointer.h" // needed for the functions to initialize arrays
#include "vtkPolyData.h" // needed as most helper functions modify output which is vtkPolyData

/* The actual class definition */
class VTK_IO_EXPORT vtkGadgetIIReader : public vtkPolyDataAlgorithm
{
public:
  static vtkGadgetIIReader* New();
  vtkTypeRevisionMacro(vtkGadgetIIReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  // Description:
  // Set/Get the name of the file from which to read the marked points.
  vtkSetStringMacro(MarkFileName);
  vtkGetStringMacro(MarkFileName);

protected:
  vtkGadgetIIReader();
  ~vtkGadgetIIReader();
  char* FileName;
	char* MarkFileName;
  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkGadgetIIReader(const vtkGadgetIIReader&);  // Not implemented.
  void operator=(const vtkGadgetIIReader&);  // Not implemented.
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them in the output vector
	void AllocateAllVariableArrays(vtkPolyData* output);
	int ReadSnapshot(FILE* gadgetInFile,vtkPolyData* output,int files);
};
#endif
