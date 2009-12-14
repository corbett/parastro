try: paraview.simple
except: from paraview.simple import *

Glyph1 = GetActiveSource()
PrincipleMomentsofInertia1 = PrincipleMomentsofInertia()

PrincipleMomentsofInertia1.SelectInputArray = ['POINTS', 'global id']

PrincipleMomentsofInertia1.SelectInputArray = ['POINTS', 'mass']

DataRepresentation10 = Show()
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation10.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation10.SelectionPointLabelJustification = 'Center'
DataRepresentation10.SelectionCellLabelJustification = 'Center'
DataRepresentation10.PointSize = 1.0
DataRepresentation10.SelectionLineWidth = 2.0
DataRepresentation10.Texture = []
DataRepresentation10.SelectionCellLabelFontSize = 24
DataRepresentation10.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation10.SelectionRepresentation = 'Wireframe'

Render()
