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
	//each array allocated in AllocateAllTipsyVariableArrays should be stored here
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
void vtkTipsyReader::ReadTipsyHeader()
{
	// Initializing for reading  
	TipsyHeader       h; 
	// Reading in the header
  this->tipsyInFile >> h;
  // Set the number of points for each scalar array; this is necessary if I want to use InsertValue for scalars by id
	this->numDark=h.h_nDark;
	this->numGas=h.h_nSph;
	this->numStar=h.h_nStar;
	this->numBodies=this->numDark+this->numGas+this->numStar;
}
//----------------------------------------------------------------------------
void vtkTipsyReader::ReadAllParticles()
{
	//allocate local variables for reading
	TipsyGasParticle  g; 
	TipsyDarkParticle d; 
	TipsyStarParticle s;
	// Allocate vtk scalars and vector arrays to hold particle data
	AllocateAllTipsyVariableArrays(); 	
	for( i=0; i<this->numDark; i++ ) 
  	{ 
		this->tipsyInFile >> d;
		ReadDarkParticle(d);
  	}
  for( i=0; i<this->numGas; i++ ) 
  	{
		this->tipsyInFile >> g;
		ReadGasParticle(g);
  	}
  for( i=0; i<this->numStar; i++) 
  	{
		this->tipsyInFile >> s;
		ReadStarParticle(s);
  	}
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadMarkedParticles()
{
	vtkErrorMacro("function to read only marked particles not yet implemented, for now reading all particles.");
	ReadAllParticles();
	/*
	//allocate local variables for reading
	TipsyGasParticle  g; 
	TipsyDarkParticle d; 
	TipsyStarParticle s;
	// Allocate vtk scalars and vector arrays to hold particle data
	AllocateAllTipsyVariableArrays(); 	
	//TODO:Implement
	*/
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
	ifstream markInFile(this->MarkFileName);
	if(!markInFile)
 		{
 		vtkErrorMacro("Error opening marked particle file: " << this->FileName << " please specify a valid mark file or none at all. For now reading all particles.");
		return 0;
 		}
	else
		{
		int mfIndex,mfBodies,mfGas,mfStar,mfDark;
		//first line
		if(markInFile >> mfBodies >> mfGas >> mfStar)
			{
	 		mfDark=mfBodies-mfGas-mfStar;
			if(mfBodies!=this->numBodies || mfDark!=this->numDark || mfGas!=this->numGas || mfStar !=this->numStar)
	 			{
	 			vtkErrorMacro("Error opening marked particle file, wrong format, number of particles do not match Tipsy file " << this->FileName << " please specify a valid mark file or none at all. For now reading all particles.");
	 			return 0;
	 			}
		else
			{
			while(markInFile >> mfIndex)
				{
				//read in the file, note the marked particles
				//insert the next marked particle index into the array keeping track of the total number of marked particles
				this->MarkedParticleIndices.push(mfIndex);
				this->totalMark++;
				}
			//done reading
			//TODO: set total bodies to total mark
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
	/*
	* Reads a file, optionally only the marked particles from the file, in the following order:
	* 1. Open Tipsy binary
	* 2. Read Tipsy header (tells us the number of particles of each type we are dealing with)
	* 3. Read mark file indices from marked particle file, if there is one
	* 4. Read either marked particles only or all particles
	*/
  // Make sure we have a file to read.
  if(!this->FileName)
    {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }
	// Open the tipsy standard file and abort if there is an error.
  this->tipsyInFile.open(this->FileName,"standard");
  if (!this->tipsyInFile.is_open()) 
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }
  // Read the header from the input
	ReadTipsyHeader();
	// Next considering whether to read in a mark file, and if so whether that reading was a success 
	int readMarkFile;
	if(this->MarkFileName) //TODO: even when empty this is returning true and trying to read.
		{
		vtkDebugMacro("Reading marked point indices from file " << this->MarkFileName);
		readMarkFile=ReadMarkedParticleIndices();
		}
  // Read every particle and add their position to be displayed, as well as relevant scalars
	if (readMarkFile)
		{
		vtkDebugMacro("Reading only the marked points in file " << this->MarkFileName << " from file " << this->FileName);
		ReadMarkedParticles(); 
		}
	else 
		{
		vtkDebugMacro("Reading all points from file " << this->FileName);
		ReadAllParticles();
		}
  // Close the tipsy file.
  this->tipsyInFile.close();
	//Storing the data in the output vector
	StoreDataRead(outputVector);  
	//Done reading
  vtkDebugMacro("Read " << this->Points->GetNumberOfPoints() << " points.");
 	return 1;
}
