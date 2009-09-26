/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib, 
this depends on a few header files as well as the Tipsylib library.

Currently only reads in standard format Tipsy files
@author corbett
=========================================================================*/
#include <math.h>
#include <assert.h>
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
vtkIdType vtkTipsyReader::ReadBaseParticle(TipsyBaseParticle& b) 
{
	vtkIdType id = this->Points->InsertNextPoint(b.pos);
  this->Verts->InsertNextCell(1, &id);
  this->VelocityVectors->InsertTupleValue(id,b.vel);
  this->MassScalars->InsertValue(id,b.mass);
  this->PhiScalars->SetValue(id,b.phi);
	return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadParticle() 
{
	//allocating variables for reading
	vtkIdType id;
	TipsyGasParticle  g;
  TipsyDarkParticle d;
  TipsyStarParticle s;
  switch(this->tipsyInFile.tellg().section()) 
		{
  	case tipsypos::gas:
			this->tipsyInFile >> g;
			id=ReadBaseParticle(g);
			this->ParticleTypes->SetValue(id, GAS);
		  this->RhoScalars->SetValue(id, g.rho);
		  this->TempScalars->SetValue(id, g.temp);
		  this->HsmoothScalars->SetValue(id, g.hsmooth);
		  this->MetalsScalars->SetValue(id, g.metals);
			break;
    case tipsypos::dark:
			this->tipsyInFile >> d;
			id=ReadBaseParticle(d);
			this->ParticleTypes->SetValue(id, DARK);
	  	this->EpsScalars->SetValue(id, d.eps);
			break;
    case tipsypos::star:
			this->tipsyInFile >> s;
			id=ReadBaseParticle(s);
			this->ParticleTypes->SetValue(id, STAR);
		  this->EpsScalars->SetValue(id, s.eps);
		  this->MetalsScalars->SetValue(id, s.metals);
		  this->TformScalars->SetValue(id, s.metals);
			break;
    default:
			assert(0);
			break;
    }
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
int vtkTipsyReader::ReadAllParticles()
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	AllocateAllTipsyVariableArrays();
	for( i=0; i<this->numDark+this->numGas+this->numStar; i++ ) 
  	{ 
		ReadParticle();
  	}
	return 1;
}

//----------------------------------------------------------------------------
int vtkTipsyReader::ReadMarkedParticles()
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	// As marked file was read, only allocates this->numBodies which now equals the number of marked particles
	AllocateAllTipsyVariableArrays();
	int nextMarkedParticleIndex;
  while(!this->MarkedParticleIndices.empty())
		{
			nextMarkedParticleIndex=this->MarkedParticleIndices.front();
			this->MarkedParticleIndices.pop();
			// navigating to the next marked particle
			if(nextMarkedParticleIndex < this->numGas)
				{
				//we are seeking a gas particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::gas,nextMarkedParticleIndex));	
				}
			else if (nextMarkedParticleIndex < this->numDark+this->numGas)
				{
				//we are seeking a dark particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::dark,nextMarkedParticleIndex));	
				}
			else if (nextMarkedParticleIndex < this->numDark+this->numGas+this->numStar)
				{
				//we are seeking a star particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::star,nextMarkedParticleIndex));	
				}
			else
				{
				vtkErrorMacro("A marked particle index is greater than the number of particles in the file, unable to read");
				return 0;
				}
			// reading in the particle
			ReadParticle();
		}
		return 1;
}

//----------------------------------------------------------------------------
int vtkTipsyReader::ReadMarkedParticleIndices()
{
	ifstream markInFile(this->MarkFileName);
	if(!markInFile)
 		{
 		vtkErrorMacro("Error opening marked particle file: " 
									<< this->MarkFileName 
									<< " please specify a valid mark file or none at all. For now reading all particles.");
		return 0;
 		}
	else
		{
		int mfIndex,mfBodies,mfGas,mfStar,mfDark;
		// first line of the mark file is of a different format: intNumBodies intNumGas intNumStars
		if(markInFile >> mfBodies >> mfGas >> mfStar)
			{
	 		mfDark=mfBodies-mfGas-mfStar;
			if(mfBodies!=this->numBodies || mfDark!=this->numDark || mfGas!=this->numGas || mfStar !=this->numStar)
	 			{
	 			vtkErrorMacro("Error opening marked particle file, wrong format, number of particles do not match Tipsy file: " 
											<< this->MarkFileName 
											<< " please specify a valid mark file or none at all. For now reading all particles.");
	 			return 0;
	 			}
			else
				{
				// Resetting the number of bodies as we will only count those marked
				this->numBodies=0; 
				// The rest of the file is is format: intMarkIndex\n
				// Read in the rest of the file file, noting the marked particles
				while(markInFile >> mfIndex)
					{
					// Insert the next marked particle index into the array keeping track of the total number of marked particles
					//subtracting one as arrays in marked particle file begin at 1, not 0
					this->MarkedParticleIndices.push(mfIndex-1);
					this->numBodies++;
					}
				// closing file
				markInFile.close();
				// read file successfully
				vtkDebugMacro("Read " << this->numBodies<< " marked point indices.");
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
 	//and if so whether that reading was a success 
	int readMarkFile,readParticles;
	// Next considering whether to read in a mark file,
	if(strcmp(this->MarkFileName,"")!=0)
		{
		vtkDebugMacro("Reading marked point indices from file:" << this->MarkFileName);
		readMarkFile=ReadMarkedParticleIndices();
		}
  // Read every particle and add their position to be displayed, as well as relevant scalars
	if (readMarkFile)
		{
		vtkDebugMacro("Reading only the marked points in file: " << this->MarkFileName << " from file " << this->FileName);
		readParticles=ReadMarkedParticles(); 
		}
	else 
		{
		vtkDebugMacro("Reading all points from file " << this->FileName);
		readParticles=ReadAllParticles();
		}
  // Close the tipsy file.
  this->tipsyInFile.close();
	// If reading was succesful, store the data
	if(readParticles)
		{
		// Storing the data in the output vector
		StoreDataRead(outputVector);  
  	vtkDebugMacro("Read " << this->Points->GetNumberOfPoints() << " points.");
 		return 1;
		}
	else
		{
		vtkErrorMacro("Was unable to read in the particles, please pick a valid tipsy binary, and if applicable a valid mark file");
		return 0;
		}
}
