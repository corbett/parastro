/*=========================================================================

		Program:   AstroViz plugin for ParaView
		Module:    $RCSfile: vtkFriendsOfFriendsHaloFinder.h,v $

		Copyright (c) Christine Corbett Moran
		All rights reserved.
   
	This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.


=========================================================================*/
// .NAME vtkFriendsOfFriendsHaloFinder 
// .SECTION Description
// vtkFriendsOfFriendsHaloFinder 
// Finds groups of particles, defined to be haloes each within a specified
//  linking length of each other. Parallel by process, and halo ids are 
//  guaranteed to be globally unique, but does not merge haloes
//  accross processes. If run in parallel requires a unique ID list as input
// otherwise, this is unused.
// .SECTION See Also
// vtkKdTree, vtkKdTree, vtkPointSetAlgorithm.h

#ifndef __vtkFriendsOfFriendsHaloFinder_h
#define __vtkFriendsOfFriendsHaloFinder_h
#include "vtkPointSetAlgorithm.h"
class vtkPointSet;
class vtkKdTree;
class vtkMultiProcessController;
class vtkIdTypeArray;

enum FriendsOfFriendsMPIData 
{
	GHOST_POINTS_AND_LOCAL_HALO_IDS,
	GHOST_POINTS_AND_LOCAL_HALO_IDS_TO_GLOBAL,
};
class VTK_GRAPHICS_EXPORT vtkFriendsOfFriendsHaloFinder : public vtkPointSetAlgorithm
{
public:
  static vtkFriendsOfFriendsHaloFinder *New();
  vtkTypeRevisionMacro(vtkFriendsOfFriendsHaloFinder,
		vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get/Set the linking length, for which any particle within this radius
	// of another particle is linked into a halo
  vtkSetMacro(LinkingLength, double);
  vtkGetMacro(LinkingLength, double);

  // Description:
  // Get/Set the minimum number of particles to consider a halo
  vtkSetMacro(MinimumNumberOfParticles, int);
  vtkGetMacro(MinimumNumberOfParticles, int);

 	// Description:
	// By defualt this filter uses the global controller,
	// but this method can be used to set another instead.
	virtual void SetController(vtkMultiProcessController*);
	vtkGetObjectMacro(Controller, vtkMultiProcessController);
	
  // Description:
	// Groups all particles within a certain linking length of eachother into
	// a single halo. Considers particles to comprise a halo only if its group
	// has more than the requisite number of particles, as input by user. 
	// Output should contain the data set in which halos should be searched
	// before calling.
	vtkIdTypeArray* FindHaloes(vtkKdTree* pointTree, 
		vtkIdTypeArray* globalIdArray, vtkPointSet* input);

//BTX
protected:
  vtkFriendsOfFriendsHaloFinder();
  ~vtkFriendsOfFriendsHaloFinder();
  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  // Main implementation.
  virtual int RequestData(vtkInformation*,
   	vtkInformationVector**,
    vtkInformationVector*);
  double LinkingLength;
	int MinimumNumberOfParticles;
	vtkMultiProcessController* Controller;

	// Description:
	// Returns unique id, simply equal to index+1 if running in serial,
	// or equal to the id at index+1 in the globalIdArray if running in parallel
	// The plus one is so that we can tell whether a unique id has been assigned
	// yet in a map (has not been assigned if = 0, which would break down if 
	// assigned ids were not strictly greater than zero)
	vtkIdType GetUniqueId(int index, vtkIdTypeArray* globalIdArray);

private:
  vtkFriendsOfFriendsHaloFinder(const vtkFriendsOfFriendsHaloFinder&);  // Not implemented.
  void operator=(const vtkFriendsOfFriendsHaloFinder&);  // Not implemented.
//ETX
};
#endif
