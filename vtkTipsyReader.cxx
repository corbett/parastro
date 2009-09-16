/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib, 
this depends on a few header files as well as the Tipsylib library.

Currently only reads in standard format Tipsy files
@author corbett
=========================================================================*/
#include "ftipsy.hpp" 
#include "vtkTipsyReader.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
//Initializing 
ifTipsy in;          // The input file
TipsyHeader       h; // The header structure
TipsyGasParticle  g; // A gas particle
TipsyDarkParticle d; // A dark particle
TipsyStarParticle s; // A star particle
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
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
    if (!in.is_open()) {
	    vtkErrorMacro("Error opening file " << this->FileName);
	    return 0;	
    }

  // Allocate objects to hold points and vertex cells.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
//  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();//added this

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  double x[3];
  // Read the header from the input
  in >> h;
  // Read every particle and add their position to be displayed
  //TODO: this code can obviously be more general; it's repeating three times, fix this
  for( i=0; i<h.h_nSph;  i++ ) {
	in >> g;
	x[0]=g.pos[0];
	x[1]=g.pos[1];
	x[2]=g.pos[2];
	vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
//	scalars->InsertNextValue(g.mass);
  }
  for( i=0; i<h.h_nDark; i++ ) { 
	in >> d;
	x[0]=d.pos[0];
	x[1]=d.pos[1];
	x[2]=d.pos[2];
	vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
//	scalars->InsertNextValue(d.mass);
  }
  for( i=0; i<h.h_nStar; i++) {
	in >> s;
	x[0]=s.pos[0];
	x[1]=s.pos[1];
	x[2]=s.pos[2];
	vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
//	scalars->InsertNextValue(s.mass);
  }
  // Close the file.
  in.close();
  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts);
//  scalars->InsertNextValue(g.mass); //doesn't work
  return 1;
}
