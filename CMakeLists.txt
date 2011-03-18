cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

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
# Author: Christine Corbett Moran, contributions by Rafael Kueng and John Biddiscombe
###

###############
# StrangeAttractors plugin
# o strange attractors random/binary visualization
# Author: Christine Corbett Moran
###



INCLUDE_DIRECTORIES(AstroVizHelpersLib)
#INCLUDE_DIRECTORIES(fio)

# for special gui
QT4_WRAP_CPP(MOC_SRCS pqTrackFilterPanel.h)

ADD_PARAVIEW_OBJECT_PANEL(
	IFACES IFACE_SRCS
	CLASS_NAME pqTrackFilterPanel
	XML_NAME TrackFilter
	XML_GROUP filters
	)

#Defines the AstroViz plugin
ADD_PARAVIEW_PLUGIN(AstroVizPlugin "1.1"
	SERVER_MANAGER_SOURCES 
		vtkTipsyReader.cxx
		vtkNSmoothFilter.cxx
		vtkCenterOfMassFilter.cxx 
		vtkProfileFilter.cxx
		vtkMomentsOfInertiaFilter.cxx
		vtkVirialRadiusFilter.cxx
		vtkAddAdditionalAttribute.cxx 
		vtkFriendsOfFriendsHaloFinder.cxx 
		vtkRamsesReader.cxx
		vtkGraficReader.cxx		
		vtkSQLiteReader2.cxx
		vtkTrackFilter.cxx
		vtkSimpleBin.cxx
	SERVER_MANAGER_XML
		TipsyReaderSM.xml
		NSmoothFilterSM.xml
		CenterOfMassFilter.xml
		ProfileFilter.xml
		MomentsOfInertiaFilter.xml 
		VirialRadiusFilter.xml
		AddAdditionalAttribute.xml 
		FriendsOfFriendsHaloFinder.xml	
		GraficReaderSM.xml
		SQLiteReader2SM.xml
		TrackFilterSM.xml
		SimpleBinSM.xml
		RamsesReaderSM.xml

	GUI_RESOURCES
		TipsyReader.qrc
		SQLiteReader2.qrc
		RamsesReader.qrc
		GraficReader.qrc		
	GUI_RESOURCE_FILES
		AstroVizFilterMenu.xml
		
	GUI_INTERFACES ${IFACES} 
		SOURCES pqTrackFilterPanel.cxx ${MOC_SRCS} ${IFACE_SRCS}
	)
#PARAVIEW_INCLUDE_WRAPPED_SOURCES(fio/fio.c fio/romberg.c)
#Defines the StrangeAttractors plugin

#ADD_PARAVIEW_PLUGIN(StrangeAttractorsPlugin "1.0"
#	SERVER_MANAGER_SOURCES 
#		vtkStrangeAttractors.cxx
#	SERVER_MANAGER_XML
#		StrangeAttractorsSM.xml
#	GUI_RESOURCES
#		StrangeAttractors.qrc 
#	GUI_INTERFACES ${IFACES} 
#	SOURCES ${MOC_SRCS} ${IFACE_SRCS}
#	)
	
# For helper functions often used, will later include these in a single
# VTK class.
ADD_LIBRARY(AstroVizHelpers AstroVizHelpersLib/AstroVizHelpers.cxx)

SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizPlugin AstroVizHelpers ) 

# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
ADD_LIBRARY(
	TipsyHelpers
	STATIC
		tipsylib/adapter.cpp
		tipsylib/binner.cpp 
		tipsylib/ftipsy.cpp tipsylib/native.cpp 
		tipsylib/standard.cpp
		tipsylib/vtipsy.cpp
	)
	
SET_TARGET_PROPERTIES(TipsyHelpers PROPERTIES COMPILE_FLAGS "-fPIC")	

# add the winsock2 library for net lookup names
IF (WIN32)
  TARGET_LINK_LIBRARIES(TipsyHelpers ws2_32)  
ENDIF (WIN32)


# For SQLiteReader2 Plugin, SQLite lib
ADD_LIBRARY(
	SQLitelib
	STATIC
		sqlitelib/sqlite3.c
	)
SET_TARGET_PROPERTIES(SQLitelib PROPERTIES COMPILE_FLAGS "-fPIC")
#For the Grafic reader plugin (and later) integrated tipsy HDF, tipsy binary, and grafic reader I use some of Doug Potter's/Joachim Stadel's fio
# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
# TODO: this is not working, hardcoding library path for now
ADD_LIBRARY(
	fio
	STATIC
	fio/fio.c
	fio/romberg.c
)
SET_TARGET_PROPERTIES(fio PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(fio PROPERTIES COMPILE_FLAGS "-fPIC")

TARGET_LINK_LIBRARIES(
	AstroVizPlugin
	fio
	TipsyHelpers
	SQLitelib
	#gmp
	)
#TARGET_LINK_LIBRARIES(AstroVizPlugin #/Users/corbett/Documents/Projects/pvaddons/ParaViz/ParaViz_src/fio/libFio.so)
