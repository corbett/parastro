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
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
 this->FileName = 0;
 this->SetNumberOfInputPorts(0);
	// Allocate objects to hold points and vertex cells.
 this->points = vtkSmartPointer<vtkPoints>::New();
 this->verts = vtkSmartPointer<vtkCellArray>::New();
 // Allocate scalars and vectors
 this->massScalars = AllocateFloatArray(1,"mass");
 this->phiScalars = AllocateFloatArray(1,"potential");
 this->epsScalars = AllocateFloatArray(1,"softening");
 this->velocityVectors = AllocateFloatArray(3,"velocity");
/*
 rhoScalars =  AllocateFloatArray(1,"rho")
 tempScalars =  AllocateFloatArray(1,"temp");
 hsmoothScalars =  AllocateFloatArray(1,"hsmooth");
 metalsScalars =  AllocateFloatArray(1,"metals");
 tformScalars =  AllocateFloatArray(1,"tform");
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
  vtkIdType id = this->points->InsertNextPoint(baseParticle.pos);
  this->verts->InsertNextCell(1, &id);
  this->velocityVectors->InsertTupleValue(id,baseParticle.vel);
  this->massScalars->InsertValue(id,baseParticle.mass);
  this->phiScalars->SetValue(id,baseParticle.phi);
  return id;
}

//allocates a float array for use as scalar (numComponents=1) or vector (numComponents > 1)
vtkSmartPointer<vtkFloatArray> vtkTipsyReader::AllocateFloatArray(int numComponents, const char* arrayName)
{
  vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
  	floatArray->SetNumberOfComponents(numComponents);
  	floatArray->SetName(arrayName);  
  return floatArray;
}

//TODO: what to do with portions of the scalar arrays which are not set. will paraview segfault?, should these be set to some default value, if their values are not known
//why is something which should only be dark, reading in values such as metals?, for now only enabling reading dark particles

/*
//TODO: ADD BACK IN WHEN READY
void vtkTipsyReader::ReadGasParticle(TipsyGasParticle& gasParticle) 
{
  vtkIdType id = ReadParticle(gasParticle);
  this->rhoScalars->SetValue(id, gasParticle.rho);
  this->tempScalars->SetValue(id, gasParticle.temp);
  this->hsmoothScalars->SetValue(id, gasParticle.hsmooth);
  this->metalsScalars->SetValue(id, gasParticle.metals);
}
*/

void vtkTipsyReader::ReadDarkParticle(TipsyDarkParticle& darkParticle) 
{
  vtkIdType id = ReadParticle(darkParticle);
  this->epsScalars->SetValue(id, darkParticle.eps);
}
/*
//TODO: ADD BACK IN WHEN READY
void vtkTipsyReader::ReadStarParticle(TipsyStarParticle& starParticle) 
{
  vtkIdType id = ReadParticle(starParticle);
  this->epsScalars->SetValue(id, starParticle.eps);
  this->metalsScalars->SetValue(id, starParticle.metals);
  this->tformScalars->SetValue(id, starParticle.metals);
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
    if (!in.is_open()) 
			{
	    vtkErrorMacro("Error opening file " << this->FileName);
	    return 0;	
    	}

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  // Read the header from the input
  in >> h;
  //set the number of points for each scalar array; this is necessary if I want to use InsertValue for scalars by id
  int numberPoints = h.h_nDark;
  massScalars->SetNumberOfTuples(numberPoints);
  phiScalars->SetNumberOfTuples(numberPoints);
  epsScalars->SetNumberOfTuples(numberPoints);
  velocityVectors->SetNumberOfTuples(numberPoints);
/*
 //TODO: ADD BACK IN WHEN READY
  rhoScalars->SetNumberOfTuples(numberPoints);
  tempScalars->SetNumberOfTuples(numberPoints);
  hsmoothScalars->SetNumberOfTuples(numberPoints);
  metalsScalars->SetNumberOfTuples(numberPoints);
  tformScalars->SetNumberOfTuples(numberPoints);
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
  output->GetPointData()->SetScalars(phiScalars); //the default scalars to be displayed
  output->GetPointData()->AddArray(massScalars);
  output->GetPointData()->AddArray(epsScalars);
  output->GetPointData()->SetVectors(velocityVectors); //the default vectors to be displayed
/*
//TODO: ADD BACK IN WHEN READY
  output->GetPointData()->AddArray(rhoScalars);
  output->GetPointData()->AddArray(tempScalars);
  output->GetPointData()->AddArray(hsmoothScalars);
  output->GetPointData()->AddArray(metalsScalars);
  output->GetPointData()->AddArray(tformScalars);
*/
  return 1;
}
