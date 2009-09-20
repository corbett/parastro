cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})


############
#Menu bar
###
# We need to wrap for Qt stuff such as signals/slots etc. to work correctly.
QT4_WRAP_CPP(MOC_SRCS_MENU SourceToolbarActions.h)
ADD_PARAVIEW_ACTION_GROUP(IFACES_MENU IFACE_SRCS_MENU 
                        CLASS_NAME SourceToolbarActions
                        GROUP_NAME "MenuBar/ParaViz")

#ADD_PARAVIEW_PLUGIN(MenuToolbar "1.0"
#	GUI_INTERFACES ${IFACES_MENU}
#	SOURCES ${MOC_SRCS_MENU} ${IFACE_SRCS_MENU} 
#	SourceToolbarActions.cxx)
###########

###############
#Tipsy binary reader
#requires an additional library, libtipsy 
#check out with
#%svn co https://svn.physik.uzh.ch/repos/itp/tipsylib
#build with
#%make
#may need to run
#%ranlib libtipsy.a
#install
#%sudo cp libtipsa.a /usr/lib
###
#define the server-side portion of the reader plugin
ADD_PARAVIEW_PLUGIN(TipsyReaderSMPlugin "1.0"
	SERVER_MANAGER_SOURCES vtkTipsyReader.cxx
	SERVER_MANAGER_XML  TipsyReaderSM.xml)
#can later do this more cleanly
TARGET_LINK_LIBRARIES(TipsyReaderSMPlugin /usr/lib/libtipsy.a)
#define the client-side portion of the reader plugin
#the bottom two were added for menu bar
ADD_PARAVIEW_PLUGIN(TipsyReaderGUIPlugin "1.0"
	GUI_RESOURCES TipsyReader.qrc
	GUI_INTERFACES ${IFACES_MENU}
	SOURCES ${MOC_SRCS_MENU} ${IFACE_SRCS_MENU} 
	SourceToolbarActions.cxx)
#requires an additional library, libtipsy 
#can later do this more cleanly
TARGET_LINK_LIBRARIES(TipsyReaderGUIPlugin /usr/lib/libtipsy.a)
#################

##########
#FOR READING TEXT FILES OF POINTS
#Just a sample of the reader plugin architecture
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
#For the toolbar
###
## We need to wrap for Qt stuff such as signals/slots etc. to work correctly.
#QT4_WRAP_CPP(MOC_SRCS SourceToolbarActions.h)
## This is a macro for adding QActionGroup subclasses automatically as toolbars.
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
#Menu bar
###
# We need to wrap for Qt stuff such as signals/slots etc. to work correctly.
#QT4_WRAP_CPP(MOC_SRCS_MENU SourceToolbarActions.h
#ADD_PARAVIEW_ACTION_GROUP(IFACES_MENU IFACE_SRCS_MENU 
#                        CLASS_NAME SourceToolbarActions
#                        GROUP_NAME "MenuBar/ParaViz")
#
#ADD_PARAVIEW_PLUGIN(MenuToolbar "1.0"
#	GUI_INTERFACES ${IFACES_MENU}
#	SOURCES ${MOC_SRCS_MENU} ${IFACE_SRCS_MENU} 
#	SourceToolbarActions.cxx)
###########