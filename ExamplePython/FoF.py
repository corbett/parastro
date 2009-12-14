try: paraview.simple
except: from paraview.simple import *

Glyph1 = GetActiveSource()
FriendsOfFriendsHaloFinder1 = FriendsOfFriendsHaloFinder()

FriendsOfFriendsHaloFinder1.SelectInputArray = ['POINTS', 'global id']

FriendsOfFriendsHaloFinder1.LinkingLength = 0.001
FriendsOfFriendsHaloFinder1.MinimumNumberOfParticles = 40

my_representation1 = GetDisplayProperties(Glyph1)
DataRepresentation7 = Show()
DataRepresentation7.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation7.SelectionCellLabelColor = [0.0, 1.0, 0.0]
DataRepresentation7.SelectionPointLabelJustification = 'Center'
DataRepresentation7.SelectionCellLabelJustification = 'Center'
DataRepresentation7.PointSize = 1.0
DataRepresentation7.ColorAttributeType = 'POINT_DATA'
DataRepresentation7.ColorArrayName = 'global id'
DataRepresentation7.SelectionLineWidth = 2.0
DataRepresentation7.Texture = []
DataRepresentation7.SelectionCellLabelFontSize = 24
DataRepresentation7.SelectionColor = [0.048416876478217748, 0.63672846570534825, 1.0]
DataRepresentation7.SelectionRepresentation = 'Wireframe'
DataRepresentation7.LookupTable = []

my_representation1.Visibility = 0

Render()
