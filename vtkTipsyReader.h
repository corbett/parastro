/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib
 -christine
=========================================================================*/
// .NAME vtkTipsyReader - Read points from a Tipsy standard binary file
// .SECTION Description
// here is the desciprtion
#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h
#include "vtkPolyDataAlgorithm.h" // needed as this class extends vtkPolyDataAlgorithm
#include "ftipsy.hpp" // needed for functions which take Tipsy particles as arguments
#include "vtkSmartPointer.h" // needed for the functions to initialize arrays
#include "vtkFloatArray.h" // needed for the functions to initialize float arrays
#include "vtkIntArray.h" // needed for functions to initialize unsigned int arrays
#include <queue> //needed for FIFO queue used to store marked particles

using std::queue;


class VTK_IO_EXPORT vtkTipsyReader : public vtkPolyDataAlgorithm
{
public:
  static vtkTipsyReader* New();
  vtkTypeRevisionMacro(vtkTipsyReader,vtkPolyDataAlgorithm);
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
  vtkTipsyReader();
  ~vtkTipsyReader();
  char* FileName;
	char MarkFileName;
  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
	// The BTX, ETX comments bracket the portion of the code which should not be
	// attempted to wrap for use by python, specifically the code which uses
	// C++ templates as this code is unable to be wrapped.
	// private variables: points, scalars, and vectors
	int numDark;
	int numGas;
	int numStar;
	int numBodies;
	ifTipsy tipsyInFile;
	//BTX
	vtkSmartPointer<vtkPoints> Points;
	vtkSmartPointer<vtkCellArray> Verts; 
	vtkSmartPointer<vtkIntArray> ParticleTypes;
	vtkSmartPointer<vtkFloatArray> MassScalars;
	vtkSmartPointer<vtkFloatArray> PhiScalars;
	vtkSmartPointer<vtkFloatArray> EpsScalars;
	vtkSmartPointer<vtkFloatArray> VelocityVectors;
	vtkSmartPointer<vtkFloatArray> RhoScalars;   
	vtkSmartPointer<vtkFloatArray> TempScalars;       
	vtkSmartPointer<vtkFloatArray> HsmoothScalars;    
	vtkSmartPointer<vtkFloatArray> MetalsScalars;
	vtkSmartPointer<vtkFloatArray> TformScalars;
	queue<int> MarkedParticleIndices;
	// private functions: initialization and reading
	// Description:
  // create a vtkFloatArray with the  name arrayName, number of components 
  // numComponents and number of tuples numTuples
  // e.g. AllocateFloatArray("density",1,100) creates a array of 100 scalar densities
  // AllocateFloatArray("velocity",3,100) creates a array of 100 vector velocities
  vtkSmartPointer<vtkFloatArray> AllocateFloatArray(const char* arrayName,int numComponents,int numTuples);
	//ETX
	// Description:
	// reads in a particle (either gas, dark or star as appropriate) from the tipsy in file of this class
  vtkIdType ReadParticle(); 
	// Description:
	// reads variables common to all particles
  vtkIdType ReadBaseParticle(TipsyBaseParticle& b); 
	// Description:
	// reads in an array of the indices of marked particles from a file, returns 0 if unsucessful, 1 otherwise
	int ReadMarkedParticleIndices();
	// Description:
	// allocates all vtk arrays for Tipsy variables
	void AllocateAllTipsyVariableArrays();
	// Description:
	// stores the data read in a vtkPolyData output vector
	void StoreDataRead(vtkInformationVector* outputVector);
	// Description:
	// Reads the tipsy header. Must be called after the tipsy file is opened, but before any marked particle file is attempted to be open
	void ReadTipsyHeader();
	// Description:
	// Reads all particles from the tipsy file
	void ReadAllParticles();
	// Description:
	// Reads only Marked particles from the tipsy file. Must be called after function ReadMarkedParticleIndices.
	void ReadMarkedParticles();
};
#endif
