/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib, 
this depends on a few header files as well as the Tipsylib library.

Currently only reads in standard format Tipsy files
@author corbett
=========================================================================*/
#include "vtkTipsyReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
//Initializing for reading  
ifTipsy in;          // The input file
TipsyHeader       h; // The header structure
TipsyGasParticle  g; // A gas particle
TipsyDarkParticle d; // A dark particle
TipsyStarParticle s; // A star particle
//points and scalars
vtkSmartPointer<vtkPoints> points;
vtkSmartPointer<vtkCellArray> verts; 
vtkSmartPointer<vtkDoubleArray> mass_scalars;
double x[3]; //for positions
double v[3]; //for velocities
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
 this->FileName = 0;
 this->SetNumberOfInputPorts(0);
  // Allocate objects to hold points and vertex cells.
 points = vtkSmartPointer<vtkPoints>::New();
 verts = vtkSmartPointer<vtkCellArray>::New();
 // Allocate Scalars
 mass_scalars = vtkSmartPointer<vtkDoubleArray>::New();//added this
 mass_scalars->SetNumberOfComponents(1);
 mass_scalars->SetName("mass");
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
void vtkTipsyReader::ReadBaseParticle(TipsyBaseParticle& baseParticle) {
	x[0]=baseParticle.pos[0];
	x[1]=baseParticle.pos[1];
	x[2]=baseParticle.pos[2];
	vtkIdType id = points->InsertNextPoint(x);
    verts->InsertNextCell(1, &id);
	mass_scalars->InsertNextValue(baseParticle.mass);
}

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

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  // Read the header from the input
  in >> h;
  
  // Read every particle and add their position to be displayed
  for( i=0; i<h.h_nSph;  i++ ) {
	in >> g;
	ReadBaseParticle(g);
  }
  for( i=0; i<h.h_nDark; i++ ) { 
	in >> d;
	ReadBaseParticle(d);
  }
  for( i=0; i<h.h_nStar; i++) {
	in >> s;
	ReadBaseParticle(s);
  }
  // Close the file.
  in.close();
  //set the number of points for each scalar array //TODO: this may not be necessary?
  mass_scalars->SetNumberOfTuples(points->GetNumberOfPoints());

  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts); 
  output->GetPointData()->SetScalars(mass_scalars);
  return 1;

//TODO: memory management, may need to call things like mass_scalars->Delete(); points->Delete(); verts->Delete(); etc.
}
