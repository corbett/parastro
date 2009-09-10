cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
##FOR READING TEXT FILES AS CSVs
#define the server-side portion of the reader plugin
#vtkTxtReader.cxx is just a copy of vtkDelimitedTextReader in the plugin directory (same as is used by CSV)
ADD_PARAVIEW_PLUGIN(TxtReaderSMPlugin "1.0"
	SERVER_MANAGER_SOURCES vtkTxtReader.cxx
	SERVER_MANAGER_XML  TxtReaderSM.xml)
#define the client-side portion of the reader plugin
ADD_PARAVIEW_PLUGIN(TxtReaderGUIPlugin "1.0"
	GUI_RESOURCES TxtReader.qrc)