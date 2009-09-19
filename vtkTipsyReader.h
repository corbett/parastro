/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib
 -christine
=========================================================================*/
// .NAME vtkTipsyReader - Read points from a Tipsy standard binary file
// .SECTION Description
// here is the desciprtion
#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h
#include "vtkPolyDataAlgorithm.h"
#include "ftipsy.hpp"
#include "vtkFloatArray.h"
#include "vtkSmartPointer.h" // For collections

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
//The following tells the Python etc. wrappers to avoid trying to wrap the following code
//as they fail on C++ templates
//BTX
  vtkSmartPointer<vtkFloatArray> AllocateFloatArray(int numComponents, char* arrayName);
//ETX
  vtkIdType ReadParticle(TipsyBaseParticle& baseParticle); //reads variables common to all particles
  void ReadGasParticle(TipsyGasParticle& gasParticle); //reads variables common to gas particles
  void ReadDarkParticle(TipsyDarkParticle& darkParticle); //reads variables common to dark particles
  void ReadStarParticle(TipsyStarParticle& starParticle); //reads variables common to star particles
};

#endif
