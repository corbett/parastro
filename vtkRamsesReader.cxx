/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkRamsesReader.cxx,v $
  Author: Christine Corbett Moran
=========================================================================*/
#include "vtkRamsesReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h" 
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


// NOTE: this is a simple port from inefficient Fortran code, needs to be optimized for memory and speed
// putting a baseline in place, however.

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
vtkSmartPointer<vtkDoubleArray> AllocateRamsesDataArray(
  vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkDoubleArray> dataArray=vtkSmartPointer<vtkDoubleArray>::New();
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

	
	if (this->GetPointArrayStatus("Age"))
    this->Age = AllocateRamsesDataArray(output,"age",1,numBodies);
  else 
    this->Age = NULL;
	
	if (this->GetPointArrayStatus("Type"))
    this->Type = AllocateRamsesDataArray(output,"type",1,numBodies);
  else 
    this->Type = NULL;
	
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
	this->PointDataArraySelection->AddArray("Age");
  this->PointDataArraySelection->AddArray("Tform");
	this->PointDataArraySelection->AddArray("Type");
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
	
	
	/*** GAS PARTICLE CONVERSION ***/
	// Here's where we want to extract gas particles. Perhaps take in a flag whether we should bother here, or not.
	
	double gas_mass_correction = 0.0;
	if(!dark_only) {
		
		// static variables 
		static int minlvl = 1, maxlvl = rsnap.m_header.levelmax;	
		
		// First ldoubleoop: 
		// compute total volume leaf cells
		// compute toal particle mass leaf cells
		//
		// From this we get:
		// averdens = tot_mass/tot_volume
		// partmass = tot_mass/dummy_parts
		// cumulative variables
		double total_mass = 0.0;
		double total_volume = 0.0;
		double rho,dx,dx3 = 0.0;
		std::vector< RAMSES::AMR::vec<double> > leaf_cell_pos;
		
		std::vector<double> leaf_cell_density;
		std::vector<double> leaf_cell_size;
		for( unsigned idomain=1; idomain<=rsnap.m_header.ncpu; ++idomain ) {
			RAMSES_tree local_tree( rsnap, idomain, maxlvl, minlvl );
			local_tree.read();		
			// Create a new hydro data object for same domain and ref. level
			// This HAS TO BE compatible with the corresponding AMR object
			RAMSES_hydro_data local_hydro_data( local_tree );
			local_hydro_data.read( "density" );
			// cumulative variables
			for(int ilevel = minlvl; ilevel < maxlvl; ++ilevel ) {
				RAMSES_tree::iterator grid_it = local_tree.begin(ilevel), temp_it;
				while( grid_it!=local_tree.end(ilevel) ){
					// Is current grid a local grid? ...
					if((unsigned)grid_it.get_domain() == idomain){
						// loop over its cells 
						for(int i=0; i<8; i++){
							// are they further refined or have we reached the max. refinement level?
							if(!grid_it.is_refined(i) || ilevel == maxlvl-1) {
								// determine its position ...
								RAMSES::AMR::vec<double> pos = local_tree.cell_pos<double>( grid_it, i );
								rho = local_hydro_data.cell_value( grid_it, i );  // I presume this is density: TODO: verify -corbett
								dx = pow(0.5,ilevel); 
								// It's a leaf
								leaf_cell_pos.push_back(pos);
								leaf_cell_density.push_back(rho);
								leaf_cell_size.push_back(dx);
								// adding to total volume and mass, the end goal of this
								dx3=pow(dx,3);
								total_volume += dx3;
								total_mass += rho*dx3;
							}
						}
					}
					++grid_it;
				}
			}
		}	
		float CORRECTIONFACTOR=8.0;
		total_mass=total_mass/CORRECTIONFACTOR;
		total_volume=total_volume/CORRECTIONFACTOR;
		// Second loop:
		// we convert to particle via above calculation, and randomly distribute them over a cell 
		// keeping track of the remainder for a final distribution
		//compute mdm prior to this
		// mdm=min dark matter particle mass*conversion factor
		//   particle_number_guess=int(partmass/mdm)+1
		double average_density = total_mass/total_volume;

		
		double particle_mass = this->ParticleMassGuess;
		if(particle_mass==0) {
			double conversion_factor = rsnap.m_header.omega_b / (rsnap.m_header.omega_m-rsnap.m_header.omega_b);
			unsigned particle_number_guess = floor(total_mass/(min_darkparticle_mass*conversion_factor))+1;
			particle_mass = total_mass/particle_number_guess;// ndm=mdm*omegab/(omegam-omegab), valid for cosmorun, otherwise m_sph, otherwise request from user
			std::cout << "Number of gas dummy particles = " << particle_number_guess << std::endl;
			std::cout << " conversion factor=" << conversion_factor << std::endl;
			std::cout << "Mdm =" << min_darkparticle_mass*conversion_factor<< std::endl;
			
		}		
		std::cout << "Total gas mass = " << total_mass << std::endl;
		std::cout << "Gas particle mass = " << particle_mass << std::endl;
		std::cout << "averdens =" << average_density << std::endl;
		std::cout <<"volume =" << total_volume << std::endl;
		
		double mass_leftover = 0.0;
		
		unsigned number_local_particles = 0;
		unsigned total_particles = 0;
		
		int gas_id=x.size();
		
		for(unsigned leaf_idx = 0; leaf_idx < leaf_cell_pos.size(); leaf_idx++) {
			// TODO TODO: where to put the conversion factor
			
			
			rho=leaf_cell_density[leaf_idx]/CORRECTIONFACTOR;
			dx=leaf_cell_size[leaf_idx];
			RAMSES::AMR::vec<double> pos=leaf_cell_pos[leaf_idx];
			
			dx3=pow(dx,3);
			number_local_particles = floor(rho * dx3/particle_mass); 
			
			// updating stats for taking care of leftovers
			total_particles+=number_local_particles;
			mass_leftover += (rho*dx3-number_local_particles * particle_mass);
			
			
			// Here we want to distribute number_local_particles over dx*3 of space
			double g_x,g_y,g_z=0;
			for(unsigned i=0; i < number_local_particles; i++) {
				gas_id+=1;				
				g_x=dRandInRange(pos.x-dx,pos.x+dx);
				g_y=dRandInRange(pos.y-dx,pos.y+dx);
				g_z=dRandInRange(pos.z-dx,pos.z+dx);
				//std::cout << "inserting new gas particle at " << g_x << "," << g_y << "," << g_z << std::endl;
				x.push_back(g_x);
				y.push_back(g_y);
				z.push_back(g_z);
				vx.push_back(0.0);
				vy.push_back(0.0);
				vz.push_back(0.0);
				mass.push_back(0.0);
				if(!dark_only){
					age.push_back(0.0);
					metals.push_back(0.0);
				}
				ids.push_back(gas_id);	
				type.push_back(RAMSES_GAS);
			}
		}
		
		
		
		std::cout << "total particles were created " << total_particles << std::endl;		
		
		// finally we want to distribute leftover_particles = floor(mass_leftover/particle_mass) over entire volume
		unsigned leftover_particles =floor(mass_leftover/particle_mass);
		double g_x,g_y,g_z=0;
		for(unsigned i=0; i < leftover_particles; i++) {
			gas_id+=1;
			g_x=dRandInRange(0,1);
			g_y=dRandInRange(0,1);
			g_z=dRandInRange(0,1);
			x.push_back(g_x);
			y.push_back(g_y);
			z.push_back(g_z);
			vx.push_back(0.0); 
			vy.push_back(0.0);
			vz.push_back(0.0);
			mass.push_back(0.0);
			if(!dark_only){
				age.push_back(0.0);
				metals.push_back(0.0);
			}
			ids.push_back(gas_id);
			type.push_back(RAMSES_GAS);
		}
		
		// finally we may have a very small amount of mass leftover which we will in principle want to distribute
		// over all particles 
		double final_mass_leftover=mass_leftover-leftover_particles*particle_mass;
		gas_mass_correction = final_mass_leftover/(-1*gas_id-1); // correct every gas particle by adding this much mass
		
		std::cout << " mass leftover " << mass_leftover << " vs. particle_mass " << particle_mass << " so finally there are some leftover particles to dist over entire volue: " << leftover_particles << std::endl;
		std::cout << "finally there is some leftover mass to distribute over all gas particles: " << final_mass_leftover << " with a mass correction of " << gas_mass_correction <<   std::endl;
		
	}		
	//--------------------------------
	// here's where the ParaView specific code comes in
	
	
	// Allocate the arrays
	this->AllocateAllRamsesVariableArrays(x.size(), output);
	// Loop through and add to PV arrays
	double pos[3];
	double vel[3];
	double rho[3] = {0.0,0.0,0.0};
	double pmass;
	
	
	// TODO: ignoring age and metals for now.
	for(unsigned i=0;i< x.size();i++) {
		pos[0]=x[i];
		pos[1]=y[i];
		pos[2]=z[i];
		vel[0]=vx[i];
		vel[1]=vy[i];
		vel[2]=vz[i];
		if(type[i]==RAMSES_GAS) {
			pmass+=gas_mass_correction;
		}

		// all particle types have this
		this->Positions->SetPoint(i, pos);
		if (this->Velocity)  this->Velocity->SetTuple(i, vel);
		if (this->Mass)      this->Mass->SetTuple1(i, mass[i]);
		if (this->Type)      this->Type->SetTuple1(i, type[i]);
		
		//
		if(!dark_only) {
			if (this->Metals)      this->Metals->SetTuple1(i, metals[i]);
			if (this->Age)      this->Age->SetTuple1(i, age[i]);
		}
		
		// These currently are unused, should be removed
		if (this->Potential) this->Potential->SetTuple1(i,0.0);		
	  if (this->RHO)         this->RHO->SetTuple(i, rho);
		if (this->Temperature) this->Temperature->SetTuple1(i, 0.0);
		if (this->Hsmooth)     this->Hsmooth->SetTuple1(i, 0.0);
		if (this->EPS)    this->EPS->SetTuple1(i, 0.0);
	}
	
	
	// Done, vis o'clock
	// can free memory allocated in vectors above
	
	 //x, y, z, vx,vy,vz,  mass, age, metals, type;
	x.clear();
	y.clear();
	z.clear();
	vx.clear();
	vy.clear();
	vz.clear();
	mass.clear();
	age.clear();
	metals.clear();
	type.clear();
	
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
	this->Age         = NULL;
  this->Tform       = NULL;
	this->Type        = NULL;
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
