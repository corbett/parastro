cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
###############
# AstroViz plugin
# INCLUDES
#   o Tipsy binary reader
#   o AstroViz menu, menu options
#   o smoothing filter
#	  o color bars equivalent to thos available in tipsy
###
#AstroViz Menu bar
QT4_WRAP_CPP(MOC_SRCS_MENU SourceToolbarActions.h)
ADD_PARAVIEW_ACTION_GROUP(IFACES_MENU IFACE_SRCS_MENU 
                        CLASS_NAME SourceToolbarActions
                        GROUP_NAME "MenuBar/AstroViz")
#define the client-side portion of the AstroViz plugin
ADD_PARAVIEW_PLUGIN(AstroVizPlugin "1.0"
	SERVER_MANAGER_SOURCES vtkTipsyReader.cxx #<-- for server
	SERVER_MANAGER_XML  TipsyReaderSM.xml #<-- for server
	GUI_RESOURCES TipsyReader.qrc #<-- for Tipsy Binary reader
	GUI_INTERFACES ${IFACES_MENU} #<-- for AstroViz menu
	SOURCES ${MOC_SRCS_MENU} ${IFACE_SRCS_MENU}  #<-- for AstroViz
	SourceToolbarActions.cxx) #<-- for menu

# for helper functions often used-removing for now
ADD_LIBRARY(AstroVizHelpers astrovizhelpers/DataSetHelpers.cxx)
SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizPlugin AstroVizHelpers)

# for the Tipsy reader plugin, a build of Doug Potter's tipsylib. 
#ADD_LIBRARY(TipsyHelpers tipsylib/adapter.cpp tipsylib/binner.cpp 
#						tipsylib/ftipsy.cpp tipsylib/native.cpp 
#						tipsylib/standard.cpp tipsylib/zstream.cpp 
#						tipsylib/vtipsy.cpp)
#SET_TARGET_PROPERTIES(TipsyHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
# still need to work on cross-platform library build
# the above works on linux, but not on mac
TARGET_LINK_LIBRARIES(AstroVizPlugin /usr/lib/libtipsy.a)

######################
# Below are samples/playing around for possible later inclusion in AstroViz
######################
##########
# For reading text files of points
# Just a sample of the reader plugin architecture
###
##define the server-side portion of the reader plugin
##vtkTxtReader.cxx is just a copy of vtkDelimitedTextReader in the plugin ##directory (same as is used by CSV)
#ADD_PARAVIEW_PLUGIN(TxtReaderSMPlugin "1.0"
#	SERVER_MANAGER_SOURCES vtkTxtReader.cxx
#	SERVER_MANAGER_XML  TxtReaderSM.xml)
##define the client-side portion of the reader plugin
#ADD_PARAVIEW_PLUGIN(TxtReaderGUIPlugin "1.0"
#	GUI_RESOURCES TxtReader.qrc)
###########

###########
# For a toolbar
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

############
# For a graphics filter
# just a sample of the architecture for a filter plugin
####
ADD_PARAVIEW_PLUGIN(AstroVizFilters "1.0"
  SERVER_MANAGER_XML NSmoothFilterSM.xml  #<-- server manager xml
  SERVER_MANAGER_SOURCES vtkNSmoothFilter.cxx #<-- source files for the new classes
	GUI_RESOURCE_FILES AstroVizFilterMenu.xml #<-- adds a new menu of filters
)
TARGET_LINK_LIBRARIES(AstroVizFilters AstroVizHelpers)
