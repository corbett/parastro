try: paraview.simple
except: from paraview.simple import *

Glyph1 = FindSource("Glyph1")
my_representation1 = GetDisplayProperties(Glyph1)
NeighborSmooth1 = GetActiveSource()
DataRepresentation8 = GetDisplayProperties(NeighborSmooth1)
Delete(DataRepresentation8)
Delete(NeighborSmooth1)
my_representation1.Visibility = 1

Render()
