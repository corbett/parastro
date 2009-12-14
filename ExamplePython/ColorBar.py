try: paraview.simple
except: from paraview.simple import *

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='global id', Position2=[0.13, 0.5], Enabled=1, LabelFontSize=12, LookupTable=[], TitleFontSize=12, Position=[0.87, 0.25] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)
Render()
