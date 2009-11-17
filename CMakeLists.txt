cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
###############
# AstroViz plugin
# INCLUDES
#   o Tipsy binary reader
#   o smoothing filter
#   o color bars equivalent to thos available in tipsy
#   o center of mass filter 
#   o moments of inertia filter
#   o profile filter
#
# Author: Christine Corbett Moran
###

#Defines the AstroViz plugin
ADD_PARAVIEW_PLUGIN(AstroVizPlugin "1.0"
	SERVER_MANAGER_SOURCES vtkTipsyReader.cxx 
	vtkNSmoothFilter.cxx vtkCenterOfMassFilter.cxx 
	vtkProfileFilter.cxx vtkMomentsOfInertiaFilter.cxx
	vtkVirialRadiusFilter.cxx vtkAddAdditionalAttribute.cxx 
	vtkFriendsOfFriendsHaloFinder.cxx #<-- for server
	SERVER_MANAGER_XML  TipsyReaderSM.xml 
	NSmoothFilterSM.xml  CenterOfMassFilter.xml
	ProfileFilter.xml MomentsOfInertiaFilter.xml 
	VirialRadiusFilter.xml AddAdditionalAttribute.xml 
	FriendsOfFriendsHaloFinder.xml #<-- for server
	GUI_RESOURCES TipsyReader.qrc AstroVizFilterMenu.xml) #<-- for client

# For helper functions often used, will later include these in a single
# VTK class
#eventially I want this to be contained within a single class:
# vtkAstroVizHelpers.cxx, already working on this
ADD_LIBRARY(AstroVizHelpers astrovizhelpers/DataSetHelpers.cxx 
	astrovizhelpers/ProfileHelpers.cxx )
SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizPlugin AstroVizHelpers) 

# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib
ADD_LIBRARY(TipsyHelpers SHARED tipsylib/adapter.cpp tipsylib/binner.cpp 
	tipsylib/ftipsy.cpp tipsylib/native.cpp 
	tipsylib/standard.cpp
	tipsylib/vtipsy.cpp)
TARGET_LINK_LIBRARIES(AstroVizPlugin TipsyHelpers)

# Extras/Experiments
###########
# Toolbar
###
## We need to wrap for Qt stuff such as signals/slots etc. to work correctly.
#QT4_WRAP_CPP(MOC_SRCS SourceToolbarActions.h)
## This is a macro for adding QActionGroup subclasses automatically
# as toolbars.
#ADD_PARAVIEW_ACTION_GROUP(IFACES IFACE_SRCS 
#                        CLASS_NAME SourceToolbarActions
#                        GROUP_NAME "ToolBar/SourceToolbar")
## Now create a plugin for the toolbar. Here we pass IFACES and IFACE_SRCS
## which are filled up by the above macro with relevant entries
#ADD_PARAVIEW_PLUGIN(SourceToolbar "1.0"
#                  GUI_INTERFACES ${IFACES}
#                  SOURCES ${MOC_SRCS} ${IFACE_SRCS} 
#                  SourceToolbarActions.cxx)
#############

###############
#AstroViz Menu bar
###
#QT4_WRAP_CPP(MOC_SRCS_MENU SourceToolbarActions.h)
#ADD_PARAVIEW_ACTION_GROUP(IFACES_MENU IFACE_SRCS_MENU 
#                        CLASS_NAME SourceToolbarActions
#                        GROUP_NAME "MenuBar/AstroViz")
# Add the below to the ADD_PARAVIEW_PLUGIN to enable menu
#	GUI_INTERFACES ${IFACES_MENU} #<-- for AstroViz menu
#	SOURCES ${MOC_SRCS_MENU} ${IFACE_SRCS_MENU}  #<-- for AstroViz
#	SourceToolbarActions.cxx) #<-- for menu
##############

