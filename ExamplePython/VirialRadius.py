try: paraview.simple
except: from paraview.simple import *

Glyph1 = GetActiveSource()
VirialRadius1 = VirialRadius( ProbeType="High Resolution Line Source" )

VirialRadius1.ProbeType.Point1 = [-0.010403109714388847, -0.061579469591379166, -0.11254894733428955]
VirialRadius1.ProbeType.Point2 = [-0.0036313242744654417, -0.055315900593996048, -0.090570949018001556]

VirialRadius1.SelectInputArray = ['POINTS', 'global id']
VirialRadius1.ProbeType = "High Resolution Line Source"

VirialRadius1.SelectInputArray = ['POINTS', 'mass']
VirialRadius1.ProbeType = "Fixed Radius Point Source"
VirialRadius1.Softening = 0.001
VirialRadius1.Delta = 0.13

my_representation1 = GetDisplayProperties(Glyph1)
DataRepresentation11 = Show()
DataRepresentation11.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation11.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation11.SelectionPointLabelJustification = 'Center'
DataRepresentation11.SelectionCellLabelJustification = 'Center'
DataRepresentation11.PointSize = 1.0
DataRepresentation11.ColorAttributeType = 'POINT_DATA'
DataRepresentation11.ColorArrayName = 'global id'
DataRepresentation11.SelectionLineWidth = 2.0
DataRepresentation11.Texture = []
DataRepresentation11.SelectionCellLabelFontSize = 24
DataRepresentation11.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation11.SelectionRepresentation = 'Wireframe'
DataRepresentation11.LookupTable = []

VirialRadius1.ProbeType.Center = [-0.015999002019237738, -0.055182199643817965, -0.10699443917338176]

my_representation1.Visibility = 0

Render()
