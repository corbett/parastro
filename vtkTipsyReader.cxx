/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib
 -christine
=========================================================================*/
#include "ftipsy.hpp" 
#include "vtkTipsyReader.h"
#include "vtkCellArray.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

//Initializing 
//TODO: currently crashes if I try to use ifTipsy, this is because libtipsy.a is not included (likely one of the other files? )
ifTipsy in;          // The input file

TipsyHeader       h; // The header structure
TipsyGasParticle  g; // A gas particle
TipsyDarkParticle d; // A dark particle
TipsyStarParticle s; // A star particle
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkTipsyReader::~vtkTipsyReader()
{
  this->SetFileName(0);
}

//----------------------------------------------------------------------------
void vtkTipsyReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";

}

//----------------------------------------------------------------------------
int vtkTipsyReader::RequestData(vtkInformation*,
                                       vtkInformationVector**,
                                       vtkInformationVector* outputVector)
{
  // Make sure we have a file to read.
  if(!this->FileName)
    {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }

    // Open the tipsy standard file and abort if there is an error.
    in.open(this->FileName,"standard");
    if ( ! in.is_open() ) {
	    vtkErrorMacro("Error opening file " << this->FileName);
	    return 0;	
    }

  // Allocate objects to hold points and vertex cells.
/*
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  double x[3];
  while(fin >> x[0] >> x[1] >> x[2])
    {
    vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
    }
  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts);
*/
  return 1;
}
