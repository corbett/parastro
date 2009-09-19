/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib, 
this depends on a few header files as well as the Tipsylib library.

Currently only reads in standard format Tipsy files
@author corbett
=========================================================================*/
#include "math.h"
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
//for velocities
//float vel[3];
//points, scalars, and vectors
vtkSmartPointer<vtkPoints> points;
vtkSmartPointer<vtkCellArray> verts; 
vtkSmartPointer<vtkFloatArray> mass_scalars;
vtkSmartPointer<vtkFloatArray> phi_scalars;
vtkSmartPointer<vtkFloatArray> eps_scalars;
vtkSmartPointer<vtkFloatArray> velocity_vectors;
/*
TODO: ADD BACK IN WHEN READY
vtkSmartPointer<vtkFloatArray> rho_scalars;   
vtkSmartPointer<vtkFloatArray> temp_scalars;       
vtkSmartPointer<vtkFloatArray> hsmooth_scalars;    
vtkSmartPointer<vtkFloatArray> metals_scalars;
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
 mass_scalars = AllocateFloatArray(1,"mass");
//potential
 phi_scalars = AllocateFloatArray(1,"potential");
 //softening
 eps_scalars = AllocateFloatArray(1,"softening");
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

 //tform
 tform_scalars = vtkSmartPointer<vtkFloatArray>::New();
 tform_scalars->SetNumberOfComponents(1);
 tform_scalars->SetName("tform");
*/
 //Allocate Vectors
 //
 velocity_vectors = AllocateFloatArray(3,"velocity");
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
  vtkIdType id = points->InsertNextPoint(baseParticle.pos);
  verts->InsertNextCell(1, &id);
  velocity_vectors->InsertTupleValue(id,baseParticle.vel);
  mass_scalars->InsertValue(id,baseParticle.mass);
  phi_scalars->SetValue(id,baseParticle.phi);
  return id;
}

//allocates a float array for use as scalar (numComponents=1) or vector (numComponents > 1)
vtkSmartPointer<vtkFloatArray> vtkTipsyReader::AllocateFloatArray(int numComponents, char* arrayName)
{
  vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
  floatArray->SetNumberOfComponents(numComponents);
  floatArray>SetName(arrayName);  
  return floatArray;
}

//TODO: what to do with portions of the scalar arrays which are not set. will paraview segfault?, should these be set to some default value, if their values are not known
//why is something which should only be dark, reading in values such as metals?, for now only enabling reading dark particles

/*
//TODO: ADD BACK IN WHEN READY
void vtkTipsyReader::ReadGasParticle(TipsyGasParticle& gasParticle) 
{
  vtkIdType id = ReadParticle(gasParticle);
  rho_scalars->SetValue(id, gasParticle.rho);
  temp_scalars->SetValue(id, gasParticle.temp);
  hsmooth_scalars->SetValue(id, gasParticle.hsmooth);
  metals_scalars->SetValue(id, gasParticle.metals);
}
*/

void vtkTipsyReader::ReadDarkParticle(TipsyDarkParticle& darkParticle) 
{
  vtkIdType id = ReadParticle(darkParticle);
  eps_scalars->SetValue(id, darkParticle.eps);
}
/*
//TODO: ADD BACK IN WHEN READY
void vtkTipsyReader::ReadStarParticle(TipsyStarParticle& starParticle) 
{
  vtkIdType id = ReadParticle(starParticle);
  eps_scalars->SetValue(id, starParticle.eps);
  metals_scalars->SetValue(id, starParticle.metals);
  tform_scalars->SetValue(id, starParticle.metals);
}
*/

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
  //set the number of points for each scalar array; this is necessary if I want to use InsertValue for scalars by id
  int numberPoints = h.h_nDark;
  mass_scalars->SetNumberOfTuples(numberPoints);
  phi_scalars->SetNumberOfTuples(numberPoints);
  eps_scalars->SetNumberOfTuples(numberPoints);
  velocity_vectors->SetNumberOfTuples(numberPoints);
/*
 //TODO: ADD BACK IN WHEN READY
  rho_scalars->SetNumberOfTuples(numberPoints);
  temp_scalars->SetNumberOfTuples(numberPoints);
  hsmooth_scalars->SetNumberOfTuples(numberPoints);
  metals_scalars->SetNumberOfTuples(numberPoints);
  tform_scalars->SetNumberOfTuples(numberPoints);
*/

  // Read every particle and add their position to be displayed, as well as relevant scalars
  for( i=0; i<h.h_nDark; i++ ) 
  { 
	in >> d;
	ReadDarkParticle(d);
  }
/*
 //TODO: ADD BACK IN WHEN READY
  for( i=0; i<h.h_nSph;  i++ ) 
  {
	in >> g;
	ReadGasParticle(g);
  }
  for( i=0; i<h.h_nStar; i++) 
  {
	in >> s;
	ReadStarParticle(s);
  }
*/
  // Close the file.
  in.close();
  
  vtkDebugMacro("Read " << points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->SetVerts(verts); 
  output->GetPointData()->SetScalars(phi_scalars); //the default scalars to be displayed
  output->GetPointData()->AddArray(mass_scalars);
  output->GetPointData()->AddArray(eps_scalars);
  output->GetPointData()->SetVectors(velocity_vectors); //the default vectors to be displayed
/*
//TODO: ADD BACK IN WHEN READY
  output->GetPointData()->AddArray(rho_scalars);
  output->GetPointData()->AddArray(temp_scalars);
  output->GetPointData()->AddArray(hsmooth_scalars);
  output->GetPointData()->AddArray(metals_scalars);
  output->GetPointData()->AddArray(tform_scalars);
*/

  return 1;
}
