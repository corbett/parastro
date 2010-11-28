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

INCLUDE_DIRECTORIES(AstroVizHelpersLib)

QT4_WRAP_CPP(MOC_SRCS pqSQLiteReaderPanel.h)
QT4_ADD_RESOURCES(RCS_SRCS pqSQLiteReaderPanel.qrc)

ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS CLASS_NAME pqSQLiteReaderPanel 
                          XML_NAME SQLiteReader XML_GROUP sources)

#Defines the AstroViz plugin
ADD_PARAVIEW_PLUGIN(AstroVizPlugin "1.0"
	SERVER_MANAGER_SOURCES 
	vtkTipsyReader.cxx
	vtkNSmoothFilter.cxx vtkCenterOfMassFilter.cxx 
	vtkProfileFilter.cxx vtkMomentsOfInertiaFilter.cxx
	vtkVirialRadiusFilter.cxx vtkAddAdditionalAttribute.cxx 
	vtkFriendsOfFriendsHaloFinder.cxx 
	#vtkPointDisplay.cxx
	vtkSQLiteReader.cxx #<-- for server
	SERVER_MANAGER_XML  TipsyReaderSM.xml
	NSmoothFilterSM.xml  CenterOfMassFilter.xml
	ProfileFilter.xml MomentsOfInertiaFilter.xml 
	VirialRadiusFilter.xml AddAdditionalAttribute.xml 
	FriendsOfFriendsHaloFinder.xml	
	#PointDisplaySM.xml
	SQLiteReaderSM.xml #<-- for server
	GUI_RESOURCES TipsyReader.qrc 
	PointDisplay.qrc #SQLiteReader.qrc #<-- for client
	GUI_RESOURCE_FILES AstroVizFilterMenu.xml
	GUI_INTERFACES ${IFACES} 
		SOURCES pqSQLiteReaderPanel.cxx ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS}) #<-- for client

	
# For helper functions often used, will later include these in a single
# VTK class.
ADD_LIBRARY(AstroVizHelpers AstroVizHelpersLib/AstroVizHelpers.cxx)
SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizPlugin AstroVizHelpers) 

# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
ADD_LIBRARY(TipsyHelpers STATIC tipsylib/adapter.cpp tipsylib/binner.cpp 
	tipsylib/ftipsy.cpp tipsylib/native.cpp 
	tipsylib/standard.cpp
	tipsylib/vtipsy.cpp)
SET_TARGET_PROPERTIES(TipsyHelpers PROPERTIES COMPILE_FLAGS "-fPIC")	
# add the winsock2 library for net lookup names
IF (WIN32)
  TARGET_LINK_LIBRARIES(TipsyHelpers ws2_32)  
ENDIF (WIN32)

# For SQLiteReader Plugin, SQLite lib
ADD_LIBRARY(SQLitelib STATIC sqlitelib/sqlite3.c)

TARGET_LINK_LIBRARIES(AstroVizPlugin TipsyHelpers SQLitelib)
