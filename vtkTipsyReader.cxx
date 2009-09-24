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
ifTipsy tipsyIn;
// used to read in header, gas, dark, and star particles respectively
TipsyHeader       h; 
TipsyGasParticle  g; 
TipsyDarkParticle d; 
TipsyStarParticle s; 
//int, if 0 not marked, if 1 marked
int MarkedParticle;
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
 this->MarkFileName = 0; // this file is optional
 this->SetNumberOfInputPorts(0); 
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
     << (this->FileName ? this->FileName : "(none)") << "\n"
		 << indent << "MarkFileName: "
     << (this->MarkFileName ? this->MarkFileName : "(none)") << "\n";
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
void vtkTipsyReader::AllocateAllTipsyVariableArrays()
{
	this->ParticleTypes = vtkSmartPointer<vtkIntArray>::New();
		ParticleTypes->SetName("particle type");
		ParticleTypes->SetNumberOfTuples(this->numBodies);
 	this->MassScalars = AllocateFloatArray("mass",1,this->numBodies);
 	this->PhiScalars = AllocateFloatArray("potential",1,this->numBodies);
 	this->EpsScalars = AllocateFloatArray("softening",1,this->numBodies);
 	this->VelocityVectors = AllocateFloatArray("velocity",3,this->numBodies);
	this->RhoScalars =  AllocateFloatArray("rho",1,this->numBodies);
 	this->TempScalars =  AllocateFloatArray("temp",1,this->numBodies);
 	this->HsmoothScalars =  AllocateFloatArray("hsmooth",1,this->numBodies);
 	this->MetalsScalars =  AllocateFloatArray("metals",1,this->numBodies);
 	this->TformScalars =  AllocateFloatArray("tform",1,this->numBodies);
}

void vtkTipsyReader::StoreDataRead(vtkInformationVector* outputVector)
{
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

int vtkTipsyReader::ReadMarkedParticleIndices()
{
	ifstream fin(this->MarkFileName);
	if(!fin)
 		{
 		vtkErrorMacro("Error opening marked particle file: " << this->FileName << " please specify a valid mark file or none at all. For now reading all particles.");
		return 0;
 		}
	else
	{
	int mfIndex,mfBodies,mfGas,mfStar,mfDark;
	//first line
	if(fin >> mfBodies >> mfGas >> mfStar)
		{
	 	mfDark=mfBodies-mfGas-mfStar;
		if(mfBodies!=this->numBodies || mfDark!=this->numDark || mfGas!=this->numGas || mfStar !=this->numStar)
	 		{
	 		vtkErrorMacro("Error opening marked particle file, wrong format, number of particles do not match Tipsy file " << this->FileName << " please specify a valid mark file or none at all. For now reading all particles.");
	 		return 0;
	 		}
		else
			{
			//read in the file, note the marked particles
			//TODO: read here
			return 1;
			}
	 	}
	}
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
  tipsyIn.open(this->FileName,"standard");
  if (!tipsyIn.is_open()) 
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }
	//All the following is to handle marked particles
	//	Open the marked file is one is specified and abort if there is an error
  // Read points from the file.
  vtkDebugMacro("Reading points from file " << this->FileName);
  // Read the header from the input
  tipsyIn >> h;
  // Set the number of points for each scalar array; this is necessary if I want to use InsertValue for scalars by id
	this->numDark=h.h_nDark;
	this->numGas=h.h_nSph;
	this->numStar=h.h_nStar;
	this->numBodies=this->numDark+this->numGas+this->numStar;
	if(this->MarkFileName)
	{
		  vtkDebugMacro("Reading marked point indices from file " << this->MarkFileName);
			//TODO: this needs to actually read the marked indices into an array
			//todo: need to actually care if it returns 0
			ReadMarkedParticleIndices();
	}
 	// Allocate vtk scalars and vector arrays to hold particle data
	AllocateAllTipsyVariableArrays();
  // Read every particle and add their position to be displayed, as well as relevant scalars
	//TODO: this needs to actually consider the array of marked particles if it exists
  for( i=0; i<this->numDark; i++ ) 
  	{ 
			tipsyIn >> d;
			ReadDarkParticle(d);
  	}
  for( i=0; i<this->numGas; i++ ) 
  	{
			tipsyIn >> g;
			ReadGasParticle(g);
  	}
  for( i=0; i<this->numStar; i++) 
  	{
			tipsyIn >> s;
			ReadStarParticle(s);
  	}
  // Close the file.
  tipsyIn.close();
	//Storing the data in the output vector
	StoreDataRead(outputVector);  
  vtkDebugMacro("Read " << this->Points->GetNumberOfPoints() << " points.");

 	return 1;
}
