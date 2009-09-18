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
#include "vtkFloatArray.h"
#include "vtkSmartPointer.h"
//Initializing for reading  
ifTipsy in;          // The input file
TipsyHeader       h; // The header structure
TipsyGasParticle  g; // A gas particle
TipsyDarkParticle d; // A dark particle
TipsyStarParticle s; // A star particle
//points, scalars, and vectors
vtkSmartPointer<vtkPoints> points;
vtkSmartPointer<vtkCellArray> verts; 
vtkSmartPointer<vtkFloatArray> mass_scalars;
vtkSmartPointer<vtkFloatArray> phi_scalars;
/*
vtkSmartPointer<vtkFloatArray> rho_scalars;        
vtkSmartPointer<vtkFloatArray> temp_scalars;       
vtkSmartPointer<vtkFloatArray> hsmooth_scalars;    
vtkSmartPointer<vtkFloatArray> metals_scalars;
vtkSmartPointer<vtkFloatArray> eps_scalars;
vtkSmartPointer<vtkFloatArray> tform_scalars;
*/
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
 this->FileName = 0;
 this->SetNumberOfInputPorts(0);

///TODO: lots of copy and paste here, remove, but first need to figure out how exactly to deal with vtkSmartPointer template class in header file
  // Allocate objects to hold points and vertex cells.
 points = vtkSmartPointer<vtkPoints>::New();
 verts = vtkSmartPointer<vtkCellArray>::New();
 // Allocate Scalars
 //mass
 mass_scalars = vtkSmartPointer<vtkFloatArray>::New();
 mass_scalars->SetNumberOfComponents(1);
 mass_scalars->SetName("mass");
//potential
 phi_scalars = vtkSmartPointer<vtkFloatArray>::New();
 phi_scalars->SetNumberOfComponents(1);
 phi_scalars->SetName("potential");
/*
 //rho
 rho_scalars = vtkSmartPointer<vtkFloatArray>::New();
 rho_scalars->SetNumberOfComponents(1);
 rho_scalars->SetName("rho");
 //temp
 temp_scalars = vtkSmartPointer<vtkFloatArray>::New();
 temp_scalars->SetNumberOfComponents(1);
 temp_scalars->SetName("temperature");
 //smooth
 hsmooth_scalars = vtkSmartPointer<vtkFloatArray>::New();
 hsmooth_scalars->SetNumberOfComponents(1);
 hsmooth_scalars->SetName("hsmooth");
 //metals
 metals_scalars = vtkSmartPointer<vtkFloatArray>::New();
 metals_scalars->SetNumberOfComponents(1);
 metals_scalars->SetName("metals");
 //softening
 eps_scalars = vtkSmartPointer<vtkFloatArray>::New();
 eps_scalars->SetNumberOfComponents(1);
 eps_scalars->SetName("softening");
 //tform
 tform_scalars = vtkSmartPointer<vtkFloatArray>::New();
 tform_scalars->SetNumberOfComponents(1);
 tform_scalars->SetName("tform");
*/
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
vtkIdType vtkTipsyReader::ReadParticle(TipsyBaseParticle& baseParticle) 
{
  vtkIdType id = points->InsertNextPoint(baseParticle.pos[0],baseParticle.pos[1],baseParticle.pos[2]);
  verts->InsertNextCell(1, &id);
  mass_scalars->SetValue(id,baseParticle.mass);
  phi_scalars->SetValue(id,baseParticle.phi);
  return id;
}

void vtkTipsyReader::ReadGasParticle(TipsyGasParticle& gasParticle) 
{
	/*
  vtkIdType id = ReadParticle(gasParticle);
  rho_scalars->SetValue(id, gasParticle.rho);
  temp_scalars->SetValue(id, gasParticle.temp);
  hsmooth_scalars->SetValue(id, gasParticle.hsmooth);
  metals_scalars->SetValue(id, gasParticle.metals);
 */
}


void vtkTipsyReader::ReadDarkParticle(TipsyDarkParticle& darkParticle) 
{
/*
  vtkIdType id = ReadParticle(darkParticle);
  eps_scalars->SetValue(id, darkParticle.eps);
*/
}

void vtkTipsyReader::ReadStarParticle(TipsyStarParticle& starParticle) 
{
/*
  vtkIdType id = ReadParticle(starParticle);
  eps_scalars->SetValue(id, starParticle.eps);
  metals_scalars->SetValue(id, starParticle.metals);
  tform_scalars->SetValue(id, starParticle.metals);
*/
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
  // Read every particle and add their position to be displayed, as well as relevant scalars
  for( i=0; i<h.h_nSph;  i++ ) {
	in >> g;
	ReadGasParticle(g);
  }
  for( i=0; i<h.h_nDark; i++ ) { 
	in >> d;
	ReadDarkParticle(d);
  }
  for( i=0; i<h.h_nStar; i++) {
	in >> s;
	ReadStarParticle(s);
  }
  // Close the file.
  in.close();
  //set the number of points for each scalar array //TODO: this may not be necessary?
  mass_scalars->SetNumberOfTuples(points->GetNumberOfPoints());
  phi_scalars->SetNumberOfTuples(points->GetNumberOfPoints());
  
  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts); 
  output->GetPointData()->SetScalars(mass_scalars);
  output->GetPointData()->SetScalars(phi_scalars);
/*
  output->GetPointData()->SetScalars(rho_scalars);
  output->GetPointData()->SetScalars(temp_scalars);
  output->GetPointData()->SetScalars(hsmooth_scalars);
  output->GetPointData()->SetScalars(metals_scalars);
  output->GetPointData()->SetScalars(eps_scalars);
  output->GetPointData()->SetScalars(tform_scalars);
*/
  //Memory management
/*
  points->Delete();
  verts->Delete();
  mass_scalars->Delete();
  phi_scalars->Delete();
  rho_scalars->Delete();
  temp_scalars->Delete();
  hsmooth_scalars->Delete();
  metals_scalars->Delete();
  eps_scalars->Delete();
  tform_scalars->Delete();
*/
  return 1;
}
