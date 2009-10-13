/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib, 
this depends on a few header files as well as the Tipsylib library.

Only reads in standard format Tipsy files.
@author corbett
=========================================================================*/
#include <math.h>
#include <assert.h>
#include "vtkGadgetIIReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h" 
#include "vtkIntArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
vtkCxxRevisionMacro(vtkGadgetIIReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGadgetIIReader);
/* useful structures */
struct io_header_1
	{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
	} header1;

int NumPart, Ngas;

struct particle_data 
	{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  float  Rho, U, Temp, Ne;
	} *P;
int *Id;
double  Time, Redshift; 
//----------------------------------------------------------------------------
vtkGadgetIIReader::vtkGadgetIIReader()
{
  this->FileName = 0;
  this->MarkFileName = 0; // this file is optional
  this->SetNumberOfInputPorts(0); 
}

//----------------------------------------------------------------------------
vtkGadgetIIReader::~vtkGadgetIIReader()
{
  this->SetFileName(0);
  this->SetMarkFileName(0);
}

//----------------------------------------------------------------------------
void vtkGadgetIIReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n"
		 << indent << "MarkFileName: "
     << (this->MarkFileName ? this->MarkFileName : "(none)") << "\n";
}

int vtkGadgetIIReader::ReadSnapshot(FILE* gadgetInFile,\
													vtkPolyData* output,int files)
{
	/*this code is currently identital to a portion of the read_snapshot.c 
	code of Volker Springel */
	// TODO: implement
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
	#define SKIP fread(&dummy, sizeof(dummy), 1, gadgetInFile);
  for(i=0, pc=1; i<files; i++, pc=pc_new)
  {
  fread(&dummy, sizeof(dummy), 1, gadgetInFile);
  fread(&header1, sizeof(header1), 1, gadgetInFile);
  fread(&dummy, sizeof(dummy), 1, gadgetInFile);
	for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
 		{
  	NumPart+= header1.npart[k];
		}
	Ngas= header1.npart[0];

  for(k=0, ntot_withmasses=0; k<5; k++)
		{
		vtkWarningMacro("header mass " << header1.mass[k]);
		
		if(header1.mass[k]==0)
			{
			vtkWarningMacro("header npart " << header1.npart[k]);
			
	    ntot_withmasses+= header1.npart[k];
			}
		}
  SKIP;
// segfault is in the next section.

/*

  for(k=0,pc_new=pc;k<6;k++)
		{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, gadgetInFile);
	      pc_new++;
	    }
		}
  SKIP;

  SKIP;
  for(k=0,pc_new=pc;k<6;k++)
		{
	  for(n=0;n<header1.npart[k];n++)
	    {
	    fread(&P[pc_new].Vel[0], sizeof(float), 3, gadgetInFile);
	    pc_new++;
	    }
		}
  SKIP;

  SKIP;
  for(k=0,pc_new=pc;k<6;k++)
		{
	  for(n=0;n<header1.npart[k];n++)
	    {
	   	fread(&Id[pc_new], sizeof(int), 1, gadgetInFile);
	    pc_new++;
	    }
	}
  SKIP;

  if(ntot_withmasses>0)
		{
		SKIP;
		}
  for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	    P[pc_new].Type=k;
      if(header1.mass[k]==0)
				{
				fread(&P[pc_new].Mass, sizeof(float), 1, gadgetInFile);
				}
	    else
				{
				P[pc_new].Mass= header1.mass[k];
				}
	    pc_new++;
	    }
	}
  if(ntot_withmasses>0)
		{		
		SKIP;
		}
 	if(header1.npart[0]>0)
		{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	    fread(&P[pc_sph].U, sizeof(float), 1, gadgetInFile);
	    pc_sph++;
	    }
	  SKIP;
	
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
      fread(&P[pc_sph].Rho, sizeof(float), 1, gadgetInFile);
      pc_sph++;
	    }
	  SKIP;
	
	  if(header1.flag_cooling)
	    {
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
				{
		  	fread(&P[pc_sph].Ne, sizeof(float), 1, gadgetInFile);
		  	pc_sph++;
				}
	    SKIP;
	    }
	  else
			{
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
				P[pc_sph].Ne= 1.0;
				pc_sph++;
	      }
			}
		}
	}
  Time= header1.time;
  Redshift= header1.time;
	*/
}
}


//----------------------------------------------------------------------------
void vtkGadgetIIReader::AllocateAllVariableArrays(vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. Storing the points and cells in the output data object.
  output->SetPoints(vtkSmartPointer<vtkPoints>::New());
  output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
}
//----------------------------------------------------------------------------
int vtkGadgetIIReader::RequestData(vtkInformation*,
                                       vtkInformationVector**,
                                       vtkInformationVector* outputVector)
{
	// Make sure we have a file to read.
  if(!this->FileName)
	  {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }
	// Open the file and abort if there is an error.
	FILE *gadgetInFile;
	if(!(gadgetInFile=fopen(this->FileName,"r")))
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }
	//All helper functions will need access to this
	vtkPolyData* output = vtkPolyData::GetData(outputVector);
  // Read the file
	//TODO: fill in with the correct number of files
	ReadSnapshot(gadgetInFile,output,1);
  // Close the in file.
	fclose(gadgetInFile);
	// Read Successfully
	vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
								<< " points.");
 	return 1;
}
