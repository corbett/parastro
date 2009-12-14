try: paraview.simple
except: from paraview.simple import *

Glyph1 = GetActiveSource()
NeighborSmooth2 = NeighborSmooth()

NeighborSmooth2.SelectInputArray = ['POINTS', 'global id']

NeighborSmooth2.NeighborNumber = 25
NeighborSmooth2.SelectInputArray = ['POINTS', 'mass']

my_representation1 = GetDisplayProperties(Glyph1)
DataRepresentation9 = Show()
DataRepresentation9.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation9.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation9.SelectionPointLabelJustification = 'Center'
DataRepresentation9.SelectionCellLabelJustification = 'Center'
DataRepresentation9.PointSize = 1.0
DataRepresentation9.ColorAttributeType = 'POINT_DATA'
DataRepresentation9.ColorArrayName = 'global id'
DataRepresentation9.SelectionLineWidth = 2.0
DataRepresentation9.Texture = []
DataRepresentation9.SelectionCellLabelFontSize = 24
DataRepresentation9.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation9.SelectionRepresentation = 'Wireframe'
DataRepresentation9.LookupTable = []

my_representation1.Visibility = 0

Render()
