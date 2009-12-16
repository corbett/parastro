#!/usr/bin/python
# Import libraries
from paraview.simple import *
from paraview import servermanager
# Connect to server 
if not servermanager.ActiveConnection:
	connection = servermanager.Connect('x01y01.zbox.physik.uzh.ch', 11111)
# Load AstroViz plugin
servermanager.LoadPlugin('libAstroVizPlugin.dylib')
# Use the TipsyReader
tipsyfile = servermanager.sources.TipsyReader(FileName='b1.00300.d0-1000.std')
# Finding the center of mass
com=servermanager.filters.CenterOfMass(Input=tipsyfile)
com.SelectInputArray=['POINTS','mass']
com.UpdatePipeline()
# Render the result
view=servermanager.CreateRenderView()
rep=servermanager.CreateRepresentation(com,view)
view.StillRender()


#########
# Make sure: PYTHONPATH includes the VTK/Wrapping/Python and the bin/ directories from your ParaView build:
#export PYTHONPATH=/Users/corbett/Documents/Projects/Work/Viz/ParaView/ParaView-3.6.1-Stable/ParaView-3.6.1-Stable_build/Utilities/VTKPythonWrapping/:/Users/corbett/Documents/Projects/Work/Viz/ParaView/ParaView-3.6.1-Stable/ParaView-3.6.1-Stable_build/bin/
###############
# Connect to Server
#
# Here we connect to the built in server.
# Can also connect to an external server e.g. with
# connection = servermanager.Connect('x01y01.zbox.physik.uzh.ch', 11111)
###############
if not servermanager.ActiveConnection:
	connection = servermanager.Connect('x01y01.zbox.physik.uzh.ch', 11111)

###############
#
# AstroViz data readers will be available under servermanager.sources
# AstroViz data filters will be available under servermanager.filters
###############

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

# Doing Some Analysis
