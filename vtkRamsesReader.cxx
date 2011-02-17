/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkRamsesReader.cxx,v $
=========================================================================*/
#include "vtkRamsesReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h" 
#include "vtkIntArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDistributedDataFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkDataArraySelection.h"
#include <cmath>
#include <assert.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cfloat>
#include <iostream>
#include <functional>
#include <algorithm>
#include "assert.h"
#include "tipsylib/ftipsy.hpp"
#include "RAMSES_particle_data.hh"
#include "RAMSES_amr_data.hh"
#include "RAMSES_hydro_data.hh"

vtkCxxRevisionMacro(vtkRamsesReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkRamsesReader);



//----------------------------------------------------------------------------
//... define the RAMSES base cell type to be of cell_locally_essential type
//... - this type allows moving between refinement levels
typedef RAMSES::AMR::cell_locally_essential<> RAMSES_cell;

//----------------------------------------------------------------------------
//... define the tree type to be of the standard RAMSES::AMR:level type
//... with cells as defined above
typedef RAMSES::AMR::tree< RAMSES_cell, RAMSES::AMR::level< RAMSES_cell > > RAMSES_tree;

//----------------------------------------------------------------------------
//... associate hydro data with the tree type defined above
typedef RAMSES::HYDRO::data< RAMSES_tree > RAMSES_hydro_data;

//----------------------------------------------------------------------------
double dRandInRange(double min, double max) {
	double d = (double)rand()/RAND_MAX;
	return min + d*(max - min);
}


enum RamsesParticleTypes 
{
	RAMSES_DARK,
	RAMSES_STAR,
	RAMSES_SINK,
	RAMSES_GAS
};


//----------------------------------------------------------------------------
vtkSmartPointer<vtkFloatArray> AllocateRamsesDataArray(
  vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkFloatArray> dataArray=vtkSmartPointer<vtkFloatArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetNumberOfTuples(numTuples);
	dataArray->SetName(arrayName);
	// initializes everything to zero
	for (int i=0; i < numComponents; ++i) {
		dataArray->FillComponent(i, 0.0);
  }
  output->GetPointData()->AddArray(dataArray);
  return dataArray;
}

//----------------------------------------------------------------------------
vtkRamsesReader::vtkRamsesReader()
{
	srand((unsigned)time(0));
  this->FileName          = 0;
  this->UpdatePiece       = 0;
  this->UpdateNumPieces   = 0;
  this->SetNumberOfInputPorts(0); 
  //
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  //
  this->Positions     = NULL;
  this->Vertices      = NULL;
  this->GlobalIds     = NULL;
  this->ParticleIndex = 0;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
  this->Tform       = NULL;
  this->Velocity    = NULL;
}

//----------------------------------------------------------------------------
vtkRamsesReader::~vtkRamsesReader()
{
  this->SetFileName(0);
  this->PointDataArraySelection->Delete();
}

//----------------------------------------------------------------------------
void vtkRamsesReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";
}

		
//----------------------------------------------------------------------------
void vtkRamsesReader::AllocateAllRamsesVariableArrays(vtkIdType numBodies,
	vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. 
  this->Positions = vtkSmartPointer<vtkPoints>::New();
  this->Positions->SetDataTypeToFloat();
  this->Positions->SetNumberOfPoints(numBodies);
  //
  this->Vertices  = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType *cells = this->Vertices->WritePointer(numBodies, numBodies*2);
  for (vtkIdType i=0; i<numBodies; ++i) {
    cells[i*2]   = 1;
    cells[i*2+1] = i;
  }

  //
  this->GlobalIds = vtkSmartPointer<vtkIdTypeArray>::New();
  this->GlobalIds->SetName("global_id");
  this->GlobalIds->SetNumberOfTuples(numBodies);

 // Storing the points and cells in the output data object.
  output->SetPoints(this->Positions);
  output->SetVerts(this->Vertices); 

  // allocate velocity first as it uses the most memory and on my win32 machine 
  // this helps load really big data without alloc failures.
  if (this->GetPointArrayStatus("Velocity")) 
    this->Velocity = AllocateRamsesDataArray(output,"velocity",3,numBodies);
  else 
    this->Velocity = NULL;
  if (this->GetPointArrayStatus("Potential")) 
    this->Potential = AllocateRamsesDataArray(output,"potential",1,numBodies);
  else 
    this->Potential = NULL;
  if (this->GetPointArrayStatus("Mass"))
    this->Mass = AllocateRamsesDataArray(output,"mass",1,numBodies);
  else 
    this->Mass = NULL;
  if (this->GetPointArrayStatus("Eps")) 
    this->EPS = AllocateRamsesDataArray(output,"eps",1,numBodies);
  else 
    this->EPS = NULL;
  if (this->GetPointArrayStatus("Rho")) 
    this->RHO = AllocateRamsesDataArray(output,"rho",1,numBodies);
  else 
    this->RHO = NULL;
  if (this->GetPointArrayStatus("Hsmooth")) 
    this->Hsmooth = AllocateRamsesDataArray(output,"hsmooth",1,numBodies);
  else 
    this->Hsmooth = NULL;
  if (this->GetPointArrayStatus("Temperature"))
    this->Temperature = AllocateRamsesDataArray(output,"temperature",1,numBodies);
  else 
    this->Temperature = NULL;

  if (this->GetPointArrayStatus("Metals"))
    this->Metals = AllocateRamsesDataArray(output,"metals",1,numBodies);
  else 
    this->Metals = NULL;
  if (this->GetPointArrayStatus("Tform"))
    this->Tform = AllocateRamsesDataArray(output,"tform",1,numBodies);
  else 
    this->Tform = NULL;
}
//----------------------------------------------------------------------------
int vtkRamsesReader::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
		-1);

  this->PointDataArraySelection->AddArray("Potential");
  this->PointDataArraySelection->AddArray("Mass");
  this->PointDataArraySelection->AddArray("Eps");
  this->PointDataArraySelection->AddArray("Rho");
  this->PointDataArraySelection->AddArray("Hsmooth");
  this->PointDataArraySelection->AddArray("Temperature");
  this->PointDataArraySelection->AddArray("Metals");
  this->PointDataArraySelection->AddArray("Tform");
  this->PointDataArraySelection->AddArray("Velocity");

	return 1;
}
/*
* Reads a file, optionally only the marked particles from the file, 
* in the following order:
* 1. Open Ramses binary
* 2. Read Ramses header (tells us the number of particles of each type we are 
*    dealing with)
* NOTE: steps 3, 5 are currently not parallel
* 3. Read mark file indices from marked particle file, if there is one
* 4. Read either marked particles only or all particles
* 5. If an attribute file is additionally specified, reads this additional
* 	 attribute into a data array, reading only those marked if necessary.
*/
//----------------------------------------------------------------------------
int vtkRamsesReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
  //
	// Make sure we have a file to read.
  //
  if(!this->FileName)
	  {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }

	// TODO: Open the Ramses standard file and abort if there is an error.

	
	//  Open the snapshot info file
	std::string filename(this->FileName);
	RAMSES::snapshot rsnap(filename , RAMSES::version3);    
	vtkDebugMacro("simulation has " << rsnap.m_header.ncpu << " domains");
	
  // Get output information
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output polydata
  vtkPolyData *output = \
      vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkPolyData> RamsesReadInitialOutput = \
      vtkSmartPointer<vtkPolyData>::New();

  // get this->UpdatePiece information
  this->UpdatePiece = \
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	this->UpdateNumPieces = \
	    outInfo->Get(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  // reset counter before reading
  this->ParticleIndex = 0;

  // so reading all particles
	// Open particle data source for first cpu
	RAMSES::PART::data local_data(filename, 1 );
	// Retrieve the available variable names from the file
	std::vector< std::string > varnames;
	local_data.get_var_names( std::back_inserter(varnames));
	
	bool dark_only = true;
	for(unsigned i=0; i<varnames.size(); i++) {
		if(varnames[i]=="age" || varnames[i]=="metallicity") {
			dark_only=false;
		}
	}
	
	std::cout 
	<< "file contains the following variables:"
	<<"\n----------------------------------\n";
	std::copy( varnames.begin(), varnames.end(), std::ostream_iterator<std::string>(std::cout,"\n"));
	std::cout << "----------------------------------\n";
	
	
	std::vector<double> x, y, z, vx,vy,vz,  mass;
	std::vector<double> age, metals;
	std::vector<int> type;
	// TODO: use bigger numbers
	std::vector<int> ids;
	
	// Actually reading the data looping over all domains
	for( unsigned int icpu=1; icpu<=rsnap.m_header.ncpu; ++icpu ) {
		std::cout << "reading star/dark if app from domain " << icpu << std::endl;
		// Open particle data source
		RAMSES::PART::data local_data(rsnap, icpu);	
		try {
			local_data.get_var<double>("position_x",std::back_inserter(x));
			local_data.get_var<double>("position_y",std::back_inserter(y));
			local_data.get_var<double>("position_z",std::back_inserter(z));
			local_data.get_var<double>("velocity_x",std::back_inserter(vx));
			local_data.get_var<double>("velocity_y",std::back_inserter(vy));
			local_data.get_var<double>("velocity_z",std::back_inserter(vz));
			local_data.get_var<double>("mass",std::back_inserter(mass));
			local_data.get_var<int>("particle_ID",std::back_inserter(ids));
			if(!dark_only) {
				local_data.get_var<double>("age",std::back_inserter(age)); 
				local_data.get_var<double>("metallicity",std::back_inserter(metals)); 
			}
		} catch(...){
			std::cerr << "something bad happened in reading the variable.\n";
			throw;
		}	
		
	}
	
	// computing minimum_darkparticle_mass and separating dark, star (later gas)
	double min_darkparticle_mass=DBL_MAX;
	for(unsigned i=0;i< x.size();i++) {
		// IDs > 0: star or dark
		// IDs < 0: gas(used) or sink(thrownaway) 
		if(ids[i] > 0) {
			if(!dark_only && age[i]!=0 ){
				type.push_back(RAMSES_STAR);
			}
			else {
				type.push_back(RAMSES_DARK);
				min_darkparticle_mass = std::min(min_darkparticle_mass,mass[i]);
			}
		}
	}

	
	
	
	
	
	
  vtkDebugMacro("Reading all points from file " << this->FileName);
    // Read Successfully
  vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
		<< " points.");
  // release memory smartpointers - just to play safe.
  this->Vertices    = NULL;
  this->GlobalIds   = NULL;
  this->Positions   = NULL;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
  this->Tform       = NULL;
  this->Velocity    = NULL;
  //
 	return 1;
}
//----------------------------------------------------------------------------
// Below : Boiler plate code to handle selection of point arrays
//----------------------------------------------------------------------------
const char* vtkRamsesReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkRamsesReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkRamsesReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name)) {
    if (status) {
      this->PointDataArraySelection->EnableArray(name);
    }
    else {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}
//----------------------------------------------------------------------------
void vtkRamsesReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkRamsesReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkRamsesReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkRamsesReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkRamsesReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
