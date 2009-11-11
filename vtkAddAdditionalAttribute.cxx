/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkAddAdditionalAttribute.cxx,v $
=========================================================================*/
#include "vtkAddAdditionalAttribute.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkMultiProcessController.h"
#include "vtkCenterOfMassFilter.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkUnsignedCharArray.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"
#include "vtkMath.h"

vtkCxxRevisionMacro(vtkAddAdditionalAttribute, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkAddAdditionalAttribute);
vtkCxxSetObjectMacro(vtkAddAdditionalAttribute,Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkAddAdditionalAttribute::vtkAddAdditionalAttribute()
{
	this->AttributeFile = 0;
	this->AttributeName = 0; 
	this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkAddAdditionalAttribute::~vtkAddAdditionalAttribute()
{
 	this->SetController(0);
  this->SetAttributeFile(0);
  this->SetAttributeName(0);
}

//----------------------------------------------------------------------------
void vtkAddAdditionalAttribute::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
	os << indent << "AttributeFile: "
	 << (this->AttributeFile ? this->AttributeFile : "(none)") << "\n";

}

//----------------------------------------------------------------------------
int vtkAddAdditionalAttribute::FillInputPortInformation(int, 
	vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

int vtkAddAdditionalAttribute::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  return 1;
}

int vtkAddAdditionalAttribute::ReadAdditionalAttributeFile(
	vtkstd::vector<int>& markedParticleIndices, vtkPointSet* output)
{
	// open file
	ifstream attributeInFile(this->AttributeFile);
	if(!attributeInFile)
 		{
 		vtkErrorMacro("Error opening marked particle file: " 
			<< this->AttributeFile 
			<< " reading only attributes defined in binary.");
 		}
	else
		{
	if(strcmp(this->AttributeName,"")==0)	
		{
		// if default has been pummeled by user, we restore it
		this->AttributeName="additional attribute";
		}
		int numBodies;
		attributeInFile >> numBodies;
		if(numBodies==output->GetPoints()->GetNumberOfPoints())
			{
			// ready to read in
			int dataIndex=0;
			double attributeData;
			if(markedParticleIndices.empty())
				{
				// read additional attribute for all particles
				AllocateDataArray(output,this->AttributeName,1,
					output->GetPoints()->GetNumberOfPoints());
				while(attributeInFile >> attributeData)
					{
					// place attribute data in output
					SetDataValue(output,this->AttributeName,dataIndex,
						&attributeData);
					dataIndex++;
					}
				}
			else
				{
				// read additional attribute only for marked particles
				AllocateDataArray(output,this->AttributeName,1,
					markedParticleIndices.size());
				int nextMarkedParticleIndex=0;
				for(vtkstd::vector<int>::iterator it = markedParticleIndices.begin();
					it != markedParticleIndices.end(); ++it)		
					{
			 		nextMarkedParticleIndex=*it;
					while(attributeInFile >> attributeData)
						{
						if(nextMarkedParticleIndex == dataIndex)
							{
							// nextMarkedParticleIndex == current particle so store
							SetDataValue(output,this->AttributeName,
								dataIndex,&attributeData);
							dataIndex++;
							// break as we are done reading current marked particle's 
							// data and want to get the index of the next marked particle,
							// if there are more marked particles
							break;
							}
						else
							{
							// skipping, not marked
							dataIndex++;
							}
						}
					}
				}
			// closing file
			attributeInFile.close();
			// successful
			return 1;
			}
		else
			{
			vtkErrorMacro("Error opening marked particle file: " 
				<< this->AttributeFile 
				<< " reading only attributes defined in binary.");
			}
		}
	// unsuccessful
	return 0;
}

//----------------------------------------------------------------------------
int vtkAddAdditionalAttribute::RequestData(vtkInformation*,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get input and output data
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet* output = vtkPolyData::GetData(outputVector);
	output->Initialize();
	output->ShallowCopy(input);
	// Make sure we have a file to read.
  if(!this->AttributeFile)
	  {
    vtkErrorMacro("An attribute file must be specified.");
    return 0;
    }
	// For now not supporting marked particle indices. TODO: support
	vtkstd::vector<int> markedParticleIndices;
	ReadAdditionalAttributeFile(markedParticleIndices,output);	
	return 1;
}
