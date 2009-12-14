try: paraview.simple
except: from paraview.simple import *

Glyph3 = GetActiveSource()
AddAdditionalAttribute1 = AddAdditionalAttribute()

AddAdditionalAttribute1.AttributeName = 'Density'
AddAdditionalAttribute1.AdditionalAttributeFile = '/Users/corbett/Documents/Projects/Work/Viz/pvaddons/testdata/b1.00300.d0-1000.den'

DataRepresentation4 = GetDisplayProperties(Glyph3)
DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation5.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation5.SelectionPointLabelJustification = 'Center'
DataRepresentation5.SelectionCellLabelJustification = 'Center'
DataRepresentation5.PointSize = 1.0
DataRepresentation5.ColorAttributeType = 'POINT_DATA'
DataRepresentation5.ColorArrayName = 'global id'
DataRepresentation5.SelectionLineWidth = 2.0
DataRepresentation5.Texture = []
DataRepresentation5.SelectionCellLabelFontSize = 24
DataRepresentation5.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation5.SelectionRepresentation = 'Wireframe'
DataRepresentation5.LookupTable = []

DataRepresentation4.Visibility = 0

Render()
