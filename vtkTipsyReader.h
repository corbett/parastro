// .NAME vtkTipsyReader - Read points from a Tipsy standard binary file
// .SECTION Description
// here is the desciprtion
#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h
#include "vtkPolyDataReader.h"
#include "tipsylib/ftipsy.hpp" // needed for functions which take Tipsy particles as arguments
#include "vtkSmartPointer.h" // needed for the functions to initialize arrays
#include "vtkPolyData.h" // needed as most helper functions modify output which is vtkPolyData
#include <queue> // needed for FIFO queue used to store marked particles
using std::vector;
class vtkCharArray;
class VTK_IO_EXPORT vtkTipsyReader : public vtkPolyDataReader
{
public:
  static vtkTipsyReader* New();
  vtkTypeRevisionMacro(vtkTipsyReader,vtkPolyDataReader);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to get additional attributes
  vtkSetStringMacro(AttributeFile);
  vtkGetStringMacro(AttributeFile);
  // Description:
  // Set/Get the name of the file from which to read the marked points.
  vtkSetStringMacro(MarkFileName);
  vtkGetStringMacro(MarkFileName);
  // Description:
  // Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);


  // Description:
  // Get/Set whether only the particles positions should be read in.
	vtkSetMacro(ReadPositionsOnly,int);
	vtkGetMacro(ReadPositionsOnly,int);

// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkTipsyReader();
  ~vtkTipsyReader();
	char* AttributeFile;
	char* MarkFileName;
	char* FileName;
	int ReadPositionsOnly;

  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
	/* Help functions for reading */
	// Description:
	// Reads the tipsy header. 
	TipsyHeader ReadTipsyHeader(ifTipsy& tipsyInfile);
	// Description:
	// Reads all particles from the tipsy file
	void ReadAllParticles(TipsyHeader& tipsyHeader,
		ifTipsy& tipsyInfile,vtkPolyData* output);
	// Description:
	// reads in a particle (either gas, dark or star as appropriate) 
	//from the tipsy in file of this class
	vtkIdType ReadParticle(ifTipsy& tipsyInFile,vtkPolyData* output);
	// Description:
	// reads variables common to all particles
	vtkIdType ReadBaseParticle(vtkPolyData* output, TipsyBaseParticle& b);
	// Description:
	// reads variables common to all gas particles
	vtkIdType ReadGasParticle(vtkPolyData* output, TipsyGasParticle& g);
	// Description:
	// reads variables common to all star particles
	vtkIdType ReadStarParticle(vtkPolyData* output, TipsyStarParticle& s);
	// Description:
	// reads variables common to all dark particles
	vtkIdType ReadDarkParticle(vtkPolyData* output, TipsyDarkParticle& d);
	// Description:
	// Reads only Marked particles from the tipsy file.
	// Must be called after function ReadMarkedParticleIndices.
	void ReadMarkedParticles(vector<int>& markedParticleIndices,
		TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile,vtkPolyData* output);
	// Description:
	// reads in an array of the indices of marked particles from a file, 
	// returns a queue of marked particles
	// which is empty if reading was unsucessful.
	vector<int> ReadMarkedParticleIndices(TipsyHeader& tipsyHeader,
		ifTipsy& tipsyInfile);
	// Description:
	// If an additional attribute file is specified, reads this file in as 
	// an additional attribute array. The points in the tipsy file must have
	// already been read in, and an additional attribute file must have
	// been input by the user. If a marked file was specified, it reads in 
	// only the particles at indices which were marked
	int ReadAdditionalAttributeFile(vector<int>& markedParticleIndices,
		TipsyHeader& tipsyHeader, vtkPolyData* output);
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them 
	// in the output vector
	void AllocateAllTipsyVariableArrays(int numBodies,vtkPolyData* output);
//ETX

};
#endif
