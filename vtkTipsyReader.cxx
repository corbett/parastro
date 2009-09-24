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
// Initializing for reading  
// the input file
ifTipsy in;
// used to read in header, gas, dark, and star particles respectively
TipsyHeader       h; 
TipsyGasParticle  g; 
TipsyDarkParticle d; 
TipsyStarParticle s; 
// Used to store which type a particle is in an int array. Later will separate 
// each type into a separate dataset
enum particle {STAR, DARK, GAS};
uint32_t i;
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
 this->FileName = 0;
 this->MarkFileName = 0;
 this->SetNumberOfInputPorts(1); // consumes a mark file
	// Allocate objects to hold points and vertex cells.
 this->Points = vtkSmartPointer<vtkPoints>::New();
 this->Verts = vtkSmartPointer<vtkCellArray>::New();
}

//----------------------------------------------------------------------------
vtkTipsyReader::~vtkTipsyReader()
{
  this->SetFileName(0);
  this->SetMarkFileName(0);
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
  vtkIdType id = this->Points->InsertNextPoint(baseParticle.pos);
  this->Verts->InsertNextCell(1, &id);
  this->VelocityVectors->InsertTupleValue(id,baseParticle.vel);
  this->MassScalars->InsertValue(id,baseParticle.mass);
  this->PhiScalars->SetValue(id,baseParticle.phi);
  return id;
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkFloatArray> vtkTipsyReader::AllocateFloatArray(const char* arrayName, int numComponents, int numTuples)
{
  vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
  	floatArray->SetNumberOfComponents(numComponents);
  	floatArray->SetName(arrayName);
		floatArray->SetNumberOfTuples(numTuples);
  return floatArray;
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadGasParticle(TipsyGasParticle& gasParticle) 
{
  vtkIdType id = ReadParticle(gasParticle);
	this->ParticleTypes->SetValue(id, GAS);
  this->RhoScalars->SetValue(id, gasParticle.rho);
  this->TempScalars->SetValue(id, gasParticle.temp);
  this->HsmoothScalars->SetValue(id, gasParticle.hsmooth);
  this->MetalsScalars->SetValue(id, gasParticle.metals);
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadDarkParticle(TipsyDarkParticle& darkParticle) 
{
  vtkIdType id = ReadParticle(darkParticle);
	this->ParticleTypes->SetValue(id, DARK);
  this->EpsScalars->SetValue(id, darkParticle.eps);
}

void vtkTipsyReader::ReadStarParticle(TipsyStarParticle& starParticle) 
{
  vtkIdType id = ReadParticle(starParticle);
	this->ParticleTypes->SetValue(id, STAR);
  this->EpsScalars->SetValue(id, starParticle.eps);
  this->MetalsScalars->SetValue(id, starParticle.metals);
  this->TformScalars->SetValue(id, starParticle.metals);
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
    if (!in.is_open()) 
			{
	    vtkErrorMacro("Error opening file " << this->FileName);
	    return 0;	
    	}

  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  // Read the header from the input
  in >> h;
  // Set the number of points for each scalar array; this is necessary if I want to use InsertValue for scalars by id
  int numTuples = h.h_nDark + h.h_nSph + h.h_nStar; //TODO: will need to be changed when particles other than dark particles are read
 	// Allocate scalars and vectors
	// Allocate object to hold particle types
 	this->ParticleTypes = vtkSmartPointer<vtkIntArray>::New();
		ParticleTypes->SetName("particle type");
		ParticleTypes->SetNumberOfTuples(numTuples);
 	this->MassScalars = AllocateFloatArray("mass",1,numTuples);
 	this->PhiScalars = AllocateFloatArray("potential",1,numTuples);
 	this->EpsScalars = AllocateFloatArray("softening",1,numTuples);
 	this->VelocityVectors = AllocateFloatArray("velocity",3,numTuples);
	this->RhoScalars =  AllocateFloatArray("rho",1,numTuples);
 	this->TempScalars =  AllocateFloatArray("temp",1,numTuples);
 	this->HsmoothScalars =  AllocateFloatArray("hsmooth",1,numTuples);
 	this->MetalsScalars =  AllocateFloatArray("metals",1,numTuples);
 	this->TformScalars =  AllocateFloatArray("tform",1,numTuples);
  // Read every particle and add their position to be displayed, as well as relevant scalars
  for( i=0; i<h.h_nDark; i++ ) 
  	{ 
			in >> d;
			ReadDarkParticle(d);
  	}
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
  // Close the file.
  in.close();
  
  vtkDebugMacro("Read " << this->Points->GetNumberOfPoints() << " points.");

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(this->Points);
  output->SetVerts(this->Verts); 
	// the default scalars to be displayed
  output->GetPointData()->SetScalars(this->PhiScalars);
	// the rest of the scalars
	output->GetPointData()->AddArray(this->ParticleTypes);
  output->GetPointData()->AddArray(this->MassScalars);
  output->GetPointData()->AddArray(this->EpsScalars);
  output->GetPointData()->AddArray(this->RhoScalars);
  output->GetPointData()->AddArray(this->TempScalars);
  output->GetPointData()->AddArray(this->HsmoothScalars);
  output->GetPointData()->AddArray(this->MetalsScalars);
  output->GetPointData()->AddArray(this->TformScalars);
	// the default vectors to be displayed
  output->GetPointData()->SetVectors(this->VelocityVectors); 
  return 1;
}
