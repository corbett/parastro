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

protected:
  vtkTipsyReader();
  ~vtkTipsyReader();
  char* FileName;
  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
	// The BTX, ETX comments bracket the portion of the code which should not be
	// attempted to wrap for use by python, specifically the code which uses
	// C++ templates as this code is unable to be wrapped.
	//BTX
	// private variables: points, scalars, and vectors
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
	// private functions: initialization and reading
	// Description:
  // create a vtkFloatArray with the  name arrayName, number of components 
  // numComponents and number of tuples numTuples
  // e.g. AllocateFloatArray("density",1,100) creates a array of 100 scalar densities
  // AllocateFloatArray("velocity",3,100) creates a array of 100 vector velocities
  vtkSmartPointer<vtkFloatArray> AllocateFloatArray(const char* arrayName,int numComponents,int numTuples);
	//ETX
	// Description:
	// reads variables common to all particles
  vtkIdType ReadParticle(TipsyBaseParticle& baseParticle); 
	// Description:
	// reads variables common to gas particles
  void ReadGasParticle(TipsyGasParticle& gasParticle); 
	// Description:
	// reads variables common to dark particles
  void ReadDarkParticle(TipsyDarkParticle& darkParticle); 
	// Description:
	// reads variables common to star particles
  void ReadStarParticle(TipsyStarParticle& starParticle); 
};
#endif
