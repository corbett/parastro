/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkGraficReader.cxx,v $
  Author: Christine Corbett Moran
=========================================================================*/
#include "vtkGraficReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h" 
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
#include "fio/fio.h"
#include "tipsylib/ftipsy.hpp"
#include "RAMSES_particle_data.hh"
#include "RAMSES_amr_data.hh"
#include "RAMSES_hydro_data.hh"

vtkCxxRevisionMacro(vtkGraficReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGraficReader);



//----------------------------------------------------------------------------
vtkSmartPointer<vtkDoubleArray> AllocateDoubleDataArray(
  vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkDoubleArray> dataArray=\
    vtkSmartPointer<vtkDoubleArray>::New();
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


vtkSmartPointer<vtkFloatArray> AllocateFloatDataArray(
																												vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkFloatArray> dataArray=\
	vtkSmartPointer<vtkFloatArray>::New();
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


vtkSmartPointer<vtkIntArray> AllocateIntDataArray(
																												vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkIntArray> dataArray=\
	vtkSmartPointer<vtkIntArray>::New();
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
vtkGraficReader::vtkGraficReader()
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
  this->Velocity    = NULL;
}

//----------------------------------------------------------------------------
vtkGraficReader::~vtkGraficReader()
{
  this->SetFileName(0);
}

//----------------------------------------------------------------------------
void vtkGraficReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";
}

		
//----------------------------------------------------------------------------
void vtkGraficReader::AllocateAllGraficVariableArrays(vtkIdType numBodies,
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
  this->Velocity = AllocateDoubleDataArray(output,"velocity",3,numBodies);
	
	this->Mass = AllocateFloatDataArray(output,"mass",1,numBodies);
	this->EPS = AllocateFloatDataArray(output,"eps",1,numBodies);
	this->RHO = AllocateFloatDataArray(output,"rho",1,numBodies);
	this->Potential = AllocateFloatDataArray(output,"potential",1,numBodies);
	this->Temperature = AllocateFloatDataArray(output,"temperature",1,numBodies);
	this->Metals = AllocateFloatDataArray(output,"metals",1,numBodies);
	this->Tform = AllocateFloatDataArray(output,"tform",1,numBodies);
	this->Type = AllocateIntDataArray(output,"type",1,numBodies);
	
}
//----------------------------------------------------------------------------
int vtkGraficReader::RequestInformation(
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
  this->PointDataArraySelection->AddArray("Velocity");

	return 1;
}
/**/
//----------------------------------------------------------------------------
int vtkGraficReader::RequestData(vtkInformation*,
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
  FIO grafic = \
  	  fioOpen(this->FileName, 0.01, 0.01);

	
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output polydata
  vtkPolyData *output = \
      vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkPolyData> GraficReadInitialOutput = \
      vtkSmartPointer<vtkPolyData>::New();

  // get this->UpdatePiece information
  this->UpdatePiece = \
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	this->UpdateNumPieces = \
	    outInfo->Get(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  // reset counter before reading
  this->ParticleIndex = 0;

	
	// not all header variables used in script, left in case of easy modifications
	// header variables
  double dExpansion,dEcosmo,dTimeOld,dUOld;
  uint64_t nTot,nGas,nDark,nStar;
  
  if(!fioGetAttr(grafic,"dTime",FIO_TYPE_DOUBLE,&dExpansion)) dExpansion=0.0;
  if (!fioGetAttr(grafic,"dEcosmo",FIO_TYPE_DOUBLE,&dEcosmo)) dEcosmo = 0.0;
  if (!fioGetAttr(grafic,"dTimeOld",FIO_TYPE_DOUBLE,&dTimeOld)) dTimeOld = 0.0;
  if (!fioGetAttr(grafic,"dUOld",FIO_TYPE_DOUBLE,&dUOld)) dUOld = 0.0;
  nTot = fioGetN(grafic,FIO_SPECIES_ALL);
  nGas  = fioGetN(grafic,FIO_SPECIES_SPH);
  nDark = fioGetN(grafic,FIO_SPECIES_DARK);
  nStar = fioGetN(grafic,FIO_SPECIES_STAR);
	
	// Allocate the arrays
	this->AllocateAllGraficVariableArrays(nTot, output);
	
  // particle variables
  uint64_t piOrder;
  double pdPos[3],pdVel[3];
  float pfMass,pfSoft,pfPot,pfRho,pfTemp,pfMetals,pfTform;
	// loop variable
	uint64_t i;
	// loop variable
	uint64_t idx;
	
	// read/write star, 
  for(i=0; i<nStar; i++) {
    fioSeek(grafic,i,FIO_SPECIES_STAR);    
    fioReadStar(grafic,
								&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot,&pfMetals,&pfTform);
		// TODO: read the rest of the variables!
		//fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//				FIO_SPECIES_STAR,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);
		idx=i;
		this->GlobalIds->SetTuple1(idx,piOrder);
		this->Type->SetTuple1(idx,FIO_SPECIES_STAR);
		this->Positions->SetPoint(idx, pdPos);
    this->Velocity->SetTuple(idx, pdVel);
		this->Mass->SetTuple1(idx,pfMass);
		this->EPS->SetTuple1(idx,pfSoft);
		this->Potential->SetTuple1(idx,pfPot);
		// only star has this		
		this->Metals->SetTuple1(idx,pfMetals);
		this->Tform->SetTuple1(idx,pfTform);
  }
	// read/write dark
  for(i=0; i<nDark; i++) {
    fioSeek(grafic,i,FIO_SPECIES_DARK);
    fioReadDark(grafic,
								&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot);
    //fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//				FIO_SPECIES_DARK,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);

		idx=i+nStar;
		// all particle types have this
		this->GlobalIds->SetTuple1(idx,piOrder);
		this->Type->SetTuple1(idx,FIO_SPECIES_DARK);
		this->Positions->SetPoint(idx, pdPos);		
    this->Velocity->SetTuple(idx, pdVel);
		this->Mass->SetTuple1(idx,pfMass);
		this->EPS->SetTuple1(idx,pfSoft);
		this->Potential->SetTuple1(idx,pfPot);
  }
	// read/write gas
  for(i=0; i<nGas; i++) {
    fioSeek(grafic,i,FIO_SPECIES_SPH);
    fioReadSph(grafic,
							 &piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot,&pfRho,&pfTemp,&pfMetals);
    //fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//				FIO_SPECIES_SPH,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);
		idx=i+nStar+nDark;
		// all particle types have this
		this->GlobalIds->SetTuple1(idx,piOrder);
		this->Type->SetTuple1(idx,FIO_SPECIES_SPH);
		this->Positions->SetPoint(idx, pdPos);		
    this->Velocity->SetTuple(idx, pdVel);
		this->Mass->SetTuple1(idx,pfMass);
		this->EPS->SetTuple1(idx,pfSoft);
		this->Potential->SetTuple1(idx,pfPot);
		// only gas has
		this->RHO->SetTuple1(idx,pfRho);
		this->Temperature->SetTuple1(idx,pfTemp);
		this->Metals->SetTuple1(idx,pfMetals);
  }
	
	
	
  vtkDebugMacro("Reading all points from file " << this->FileName);
    // Read Successfully
  vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
		<< " points.");
  // release memory smartpointers - just to play safe.
  this->Vertices    = NULL;
  this->GlobalIds   = NULL;
  this->Positions   = NULL;
  this->Velocity    = NULL;
  //
 	return 1;
}
