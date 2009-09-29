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
#include "vtkFloatArray.h" 
#include "vtkIntArray.h" 
// Used to store which type a particle is in an int array. Later will separate 
// each type into a separate dataset
enum particle {STAR, DARK, GAS};
vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->FileName = 0;
  this->MarkFileName = 0; // this file is optional
  this->SetNumberOfInputPorts(0); 
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
vtkIdType vtkTipsyReader::ReadBaseParticle(vtkPolyData* output, TipsyBaseParticle& b) 
{
	//TODO: change this-> to use output vector instead
  vtkIdType id = this->Points->InsertNextPoint(b.pos);
  this-Verts->InsertNextCell(1, &id);
	///////////////////////////////////////////////////
	SetDataValue(output,"velocity",id,b.vel);
  SetDataValue(id,"mass",b.mass);
	SetDataValue(id,"potential",b.phi);
	return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadParticle(ifTipsy& tipsyInFile,vtkPolyData* output) 
{
  //allocating variables for reading
  vtkIdType id;
  TipsyGasParticle  g;
  TipsyDarkParticle d;
  TipsyStarParticle s;
  switch(tipsyInFile.tellg().section()) 
		{
   case tipsypos::gas:
     tipsyInFile >> g;
     id=ReadBaseParticle(g);
     SetDataValue(output,"particle",id,GAS);
     SetDataValue(output,"rho",id,g.rho);
     SetDataValue(output,"temperature",id,g.temp);
     SetDataValue(output,"hsmooth",id,g.hsmooth);
     SetDataValue(output,"metals",id,g.metals);
     break;
   case tipsypos::dark:
     tipsyInFile >> d;
     id=ReadBaseParticle(d);
     SetDataValue(output,"particle",id,DARK);
     SetDataValue(output,"eps",id,d.eps);
     break;
   case tipsypos::star:
     tipsyInFile >> s;
     id=ReadBaseParticle(s);
     SetDataValue(output,"particle",id,STAR);
     SetDataValue(output,"eps",id,s.eps);
     SetDataValue(output,"metals",id,s.metals);
     SetDataValue(output,"tform",id,s.tform);
     break;
   default:
     assert(0);
     break;
    }
  return id;
}

//----------------------------------------------------------------------------
void vtkTipsyReader::AllocateAllTipsyVariableArrays(TipsyHeader tipsyHeader,vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. Storing the points and cells in the output data object.
  output->SetPoints(vtkSmartPointer<vtkPoints>::New());
  output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	// the default scalars to be displayed
  output->GetPointData()->SetScalars(AllocateDataArray<vtkFloatArray>("potential",1,tipsyHeader.h_nBodies));
	// the rest of the scalars
	output->GetPointData()->AddArray(AllocateDataArray<vtkIntArray>("particle",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("mass",1,tipsyHeader.h_nBodies););
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("eps",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("rho",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("hsmooth",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("temperature",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("metals",1,tipsyHeader.h_nBodies));
  output->GetPointData()->AddArray(AllocateDataArray<vtkFloatArray>("tform",1,tipsyHeader.h_nBodies));
	// the default vectors to be displayed
  output->GetPointData()->SetVectors(AllocateDataArray<vtkFloatArray>("velocity",3,tipsyHeader.h_nBodies));
}

//----------------------------------------------------------------------------
template <class T> vtkSmartPointer<T> vtkTipsyReader::AllocateDataArray(const char* arrayName, int numComponents, int numTuples)
{
	vtkSmartPointer<T> dataArray=vtkSmartPointer<T>::New();
  dataArray->SetNumberOfComponents(numComponents);
  dataArray->SetName(arrayName);
	dataArray->SetNumberOfTuples(numTuples);
  return dataArray;
}

template <class T> void vtkTipsyReader::SetDataValue(vtkPolyData* output, const char* arrayName, vtkIdType id, T data)
{
	//TODO: implement
}

//----------------------------------------------------------------------------
TipsyHeader vtkTipsyReader::ReadTipsyHeader(ifTipsy tipsyInfile)
{
	// Initializing for reading  
	TipsyHeader       h; 
	// Reading in the header
  tipsyInFile >> h;
	return h;
}
//----------------------------------------------------------------------------
void vtkTipsyReader::ReadAllParticles(TipsyHeader tipsyHeader, ifTipsy tipsyInfile,vtkPolyData* output)
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	AllocateAllTipsyVariableArrays(tipsyHeader,output);
	for(int i=0; i<tipsyHeader.h_nBodies; i++) 
  	{ 
			ReadParticle(tipsyInfile,output);
  	}
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadMarkedParticles(TipsyHeader tipsyHeader,ifTipsy tipsyInfile,vtkPolyData* output)
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	// As marked file was read, only allocates numBodies which now equals the number of marked particles
	AllocateAllTipsyVariableArrays(tipsyHeader,output);
	int nextMarkedParticleIndex;
  while(!this->MarkedParticleIndices.empty())
		{
			//TODO: switch this from a global index
			nextMarkedParticleIndex=this->MarkedParticleIndices.front();
			this->MarkedParticleIndices.pop();
			// navigating to the next marked particle
			if(nextMarkedParticleIndex < tipsyHeader.h_nSph)
				{
				//we are seeking a gas particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::gas,nextMarkedParticleIndex));	
				}
			else if (nextMarkedParticleIndex < tipsyHeader.h_nDark+tipsyHeader.h_nSph)
				{
				//we are seeking a dark particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::dark,nextMarkedParticleIndex));	
				}
			else if (nextMarkedParticleIndex < this->tipsyHeader.h_nBodies)
				{
				//we are seeking a star particle
				this->tipsyInFile.seekg(tipsypos(tipsypos::star,nextMarkedParticleIndex));	
				}
			else
				{
				vtkErrorMacro("A marked particle index is greater than the number of particles in the file, unable to read")	
				break;
				}
			// reading in the particle
			ReadParticle(tipsyInfile,output);
		}
}

//----------------------------------------------------------------------------
int vtkTipsyReader::ReadMarkedParticleIndices(TipsyHeader tipsyHeader,ifTipsy tipsyInfile)
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
			int mfIndex,mfBodies,mfGas,mfStar,mfDark,numBodies;
		// first line of the mark file is of a different format: intNumBodies intNumGas intNumStars
		if(markInFile >> mfBodies >> mfGas >> mfStar)
			{
	 		mfDark=mfBodies-mfGas-mfStar;
			if(mfBodies!=tipsyHeader.h_nBodies || mfDark!=tipsyHeader.h_nDark || mfGas!=tipsyHeader.h_nGas || mfStar!=tipsyHeader.h_nStar)
	 			{
	 			vtkErrorMacro("Error opening marked particle file, wrong format, number of particles do not match Tipsy file: " 
											<< this->MarkFileName 
											<< " please specify a valid mark file or none at all. For now reading all particles.");
	 			return 0;
	 			}
			else
				{
				// Resetting the number of bodies as we will only count those marked
				numBodies=0; 
				// The rest of the file is is format: intMarkIndex\n
				// Read in the rest of the file file, noting the marked particles
				while(markInFile >> mfIndex)
					{
					// Insert the next marked particle index into the array keeping track of the total number of marked particles
					//subtracting one as arrays in marked particle file begin at 1, not 0
					//TODO: remove this global reference as well
					this->MarkedParticleIndices.push(mfIndex-1);
					numBodies++;
					}
				// closing file
				markInFile.close();
				// read file successfully
				vtkDebugMacro("Read " << numBodies<< " marked point indices.");
				return numBodies;
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
	ifTipsy tipsyInfile;
  tipsyInfile.open(this->FileName,"standard");
  if (!this->tipsyInfile.is_open()) 
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }
	//All helper functions will need access to this
	vtkPolyData* output = vtkPolyData::GetData(outputVector);
  // Read the header from the input
	TipsyHeader tipsyHeader=ReadTipsyHeader(tipsyInfile);
	// Next considering whether to read in a mark file, and if so whether that reading was a success 
	int numberMarkedParticles;
	if(strcmp(this->MarkFileName,"")!=0)
		{
		vtkDebugMacro("Reading marked point indices from file:" << this->MarkFileName);
		numberMarkedParticles=ReadMarkedParticleIndices(tipsyHeader,tipsyInfile);
		}
  // Read every particle and add their position to be displayed, as well as relevant scalars
	if (numberMarkedParticles)
		{
		vtkDebugMacro("Reading only the marked points in file: " << this->MarkFileName << " from file " << this->FileName);
		// now the number of bodies is equal to the number of marked particles
		tipsyHeader.h_nBodies=numberMarkedParticles;
		ReadMarkedParticles(tipsyHeader,tipsyInfile,output);
		}
	else 
		{
		// no marked particle file or there was an error reading the mark file
		vtkDebugMacro("Reading all points from file " << this->FileName);
		ReadAllParticles(tipsyHeader,tipsyInfile,output);
		}
  // Close the tipsy file.
	tipsyInfile.close();
	// Read Successfully
	//todo, change to access output vector
	//  vtkDebugMacro("Read " << this->Points->GetNumberOfPoints() << " points.");
 	return 1;
}
