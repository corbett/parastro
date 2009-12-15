#!/usr/bin/python
#########
# Make sure: PYTHONPATH includes the VTK/Wrapping/Python and the bin/ directories from your ParaView build:
#export PYTHONPATH=/Users/corbett/Documents/Projects/Work/Viz/ParaView/ParaView-3.6.1-Stable/ParaView-3.6.1-Stable_build/Utilities/VTKPythonWrapping/:/Users/corbett/Documents/Projects/Work/Viz/ParaView/ParaView-3.6.1-Stable/ParaView-3.6.1-Stable_build/bin/
from paraview.simple import *
from paraview import servermanager
###############
# Connect to Server
#
# Here we connect to the built in server.
# Can also connect to an external server e.g. with
# connection = servermanager.Connect('x01y01.zbox.physik.uzh.ch', 11111)
###############
if not servermanager.ActiveConnection:
 	connection = servermanager.Connect()

###############
# Load AstroViz plugin
#
# AstroViz data readers will be available under servermanager.sources
# AstroViz data filters will be available under servermanager.filters
###############
servermanager.LoadPlugin("/Users/corbett/Documents/Projects/Work/Viz/pvaddons/ParaViz/ParaViz1.3_build/libAstroVizPlugin.dylib")

###############
# Get Help 
#
# Lists available modules:
# dir(servermanager.sources.TipsyReader)
#
# Lists documentation:
# help(servermanager.sources.TipsyReader)
###############

###############
# Use the TipsyReader
#
# AaaFileName rather than simply FileName is an artifact of a ParaView bug, 
# it must be alphabetically first for GUI to function properly.
###############
b1_00300_d01000_std = servermanager.sources.TipsyReader( AaaFileName='/Users/corbett/Documents/Projects/Work/Viz/pvaddons/testdata/b1.00300.d0-1000.std' )

###############
# Create a Glyph
###############
glyph1=Glyph(Input=b1_00300_d01000_std)
glyph1.GlyphType="2D Glyph"
glyph1.GlyphType.GlyphType="Vertex"
glyph1.MaskPoints=0

###############
# Render what we have so far
###############
view1=servermanager.CreateRenderView()
rep1=servermanager.CreateRepresentation(glyph1,view1)
rep1.PointSize=1
view1.ResetCamera()
view1.StillRender()

###############
# Doing Some Analysis!
###############
# finding the center of mass
com1=servermanager.filters.CenterOfMass(Input=b1_00300_d01000_std)
com1.SelectInputArray=['POINTS','mass']
com1.UpdatePipeline()
view2=servermanager.CreateRenderView()
rep2=servermanager.CreateRepresentation(com1,view2)
rep2.PointSize=3
view2.ResetCamera()
view2.StillRender()
