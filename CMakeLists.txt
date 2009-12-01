cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
###############
# AstroViz plugin
# INCLUDES
#   o Tipsy binary reader
#   o color bars equivalent to thos available in tipsy
#   o smoothing filter
#   o center of mass filter 
#   o moments of inertia filter
#   o profile filter
#   o add additional attribute filter
#   o friends-of-friends halo finder filter
# Author: Christine Corbett Moran
###

#Defines the AstroViz plugin
ADD_PARAVIEW_PLUGIN(AstroVizPlugin "1.0"
	SERVER_MANAGER_SOURCES vtkTipsyReader.cxx 
	vtkNSmoothFilter.cxx vtkCenterOfMassFilter.cxx 
	vtkProfileFilter.cxx vtkMomentsOfInertiaFilter.cxx
	vtkVirialRadiusFilter.cxx vtkAddAdditionalAttribute.cxx 
	vtkFriendsOfFriendsHaloFinder.cxx vtkKdTreeDemo.cxx #<-- for server
	SERVER_MANAGER_XML  TipsyReaderSM.xml 
	NSmoothFilterSM.xml  CenterOfMassFilter.xml
	ProfileFilter.xml MomentsOfInertiaFilter.xml 
	VirialRadiusFilter.xml AddAdditionalAttribute.xml 
	FriendsOfFriendsHaloFinder.xml KdTreeDemo.xml  #<-- for server
	GUI_RESOURCES TipsyReader.qrc  #<-- for client
	GUI_RESOURCE_FILES AstroVizFilterMenu.xml) #<-- for client

# For helper functions often used, will later include these in a single
# VTK class.
ADD_LIBRARY(AstroVizHelpers AstroVizHelpersLib/AstroVizHelpers.cxx)
SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizPlugin AstroVizHelpers) 

# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
ADD_LIBRARY(TipsyHelpers SHARED tipsylib/adapter.cpp tipsylib/binner.cpp 
	tipsylib/ftipsy.cpp tipsylib/native.cpp 
	tipsylib/standard.cpp
	tipsylib/vtipsy.cpp)
TARGET_LINK_LIBRARIES(AstroVizPlugin TipsyHelpers)