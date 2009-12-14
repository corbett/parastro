try: paraview.simple
except: from paraview.simple import *

Glyph1 = GetActiveSource()
Profile1 = Profile( ProbeType="High Resolution Line Source" )

Profile1.ProbeType.Point1 = [-0.010403109714388847, -0.061579469591379166, -0.11254894733428955]
Profile1.ProbeType.Point2 = [-0.0036313242744654417, -0.055315900593996048, -0.090570949018001556]

Profile1.SelectInputArray = ['POINTS', 'global id']
Profile1.ProbeType = "High Resolution Line Source"

Profile1.SelectInputArray = ['POINTS', 'mass']
Profile1.BinNumber = 20

SpreadSheetView1 = CreateRenderView()

DataRepresentation12 = Show()
DataRepresentation12.FieldAssociation = 'Row Data'

Render()
