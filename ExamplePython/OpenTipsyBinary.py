try: paraview.simple
except: from paraview.simple import *

b1_00300_d01000_std = TipsyReader( AaaFileName='/Users/corbett/Documents/Projects/Work/Viz/pvaddons/testdata/b1.00300.d0-1000.std' )

b1_00300_d01000_std.MarkFileName = ''
b1_00300_d01000_std.DistributeDataOn = 0

DataRepresentation1 = Show()
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation1.SelectionPointLabelJustification = 'Center'
DataRepresentation1.SelectionCellLabelJustification = 'Center'
DataRepresentation1.SelectionLineWidth = 2.0
DataRepresentation1.ScalarOpacityUnitDistance = 0.023835305423260164
DataRepresentation1.SelectionCellLabelFontSize = 24
DataRepresentation1.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation1.SelectionRepresentation = 'Wireframe'

Render()
