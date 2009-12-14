try: paraview.simple
except: from paraview.simple import *

Profile1 = GetActiveSource()
PlotData1 = PlotData()

XYPlotView1 = CreateXYPlotView()

DataRepresentation13 = Show()
DataRepresentation13.AttributeType = 'Row Data'
DataRepresentation13.SeriesVisibility = ['velocity_total (0)', '0', 'velocity_total (1)', '0', 'velocity_total (2)', '0', 'velocity_average (0)', '0', 'velocity_average (1)', '0', 'velocity_average (2)', '0', 'velocity_cumulative (0)', '0', 'velocity_cumulative (1)', '0', 'velocity_cumulative (2)', '0', 'GlyphVector_total (0)', '0', 'GlyphVector_total (1)', '0', 'GlyphVector_total (2)', '0', 'GlyphVector_average (0)', '0', 'GlyphVector_average (1)', '0', 'GlyphVector_average (2)', '0', 'GlyphVector_cumulative (0)', '0', 'GlyphVector_cumulative (1)', '0', 'GlyphVector_cumulative (2)', '0', 'angular momentum_average (0)', '0', 'angular momentum_average (1)', '0', 'angular momentum_average (2)', '0', 'radial velocity_average (0)', '0', 'radial velocity_average (1)', '0', 'radial velocity_average (2)', '0', 'tangential velocity_average (0)', '0', 'tangential velocity_average (1)', '0', 'tangential velocity_average (2)', '0', 'velocity dispersion_total (0)', '0', 'velocity dispersion_total (1)', '0', 'velocity dispersion_total (2)', '0', 'tangential velocity dispersion_total (0)', '0', 'tangential velocity dispersion_total (1)', '0', 'tangential velocity dispersion_total (2)', '0', 'radial velocity dispersion_total (0)', '0', 'radial velocity dispersion_total (1)', '0', 'radial velocity dispersion_total (2)', '0', '', '1']

XYPlotView1.ChartTitle = 'Density vs. Radius'
XYPlotView1.ShowLegend = 0
XYPlotView1.AxisLabelPrecision = [4, 3, 2, 2]
XYPlotView1.AxisTitle = ['Density', 'Radius', '', '']

DataRepresentation13.XArrayName = 'bin radius_total'
DataRepresentation13.SeriesColor = ['density_total', '1', '0.352972', '0.222614']
DataRepresentation13.SeriesLineStyle = ['density_total', '1']
DataRepresentation13.SeriesMarkerStyle = ['density_total', '1']
DataRepresentation13.SeriesVisibility = ['velocity_total (0)', '0', 'velocity_total (1)', '0', 'velocity_total (2)', '0', 'velocity_average (0)', '0', 'velocity_average (1)', '0', 'velocity_average (2)', '0', 'velocity_cumulative (0)', '0', 'velocity_cumulative (1)', '0', 'velocity_cumulative (2)', '0', 'GlyphVector_total (0)', '0', 'GlyphVector_total (1)', '0', 'GlyphVector_total (2)', '0', 'GlyphVector_average (0)', '0', 'GlyphVector_average (1)', '0', 'GlyphVector_average (2)', '0', 'GlyphVector_cumulative (0)', '0', 'GlyphVector_cumulative (1)', '0', 'GlyphVector_cumulative (2)', '0', 'angular momentum_average (0)', '0', 'angular momentum_average (1)', '0', 'angular momentum_average (2)', '0', 'radial velocity_average (0)', '0', 'radial velocity_average (1)', '0', 'radial velocity_average (2)', '0', 'tangential velocity_average (0)', '0', 'tangential velocity_average (1)', '0', 'tangential velocity_average (2)', '0', 'velocity dispersion_total (0)', '0', 'velocity dispersion_total (1)', '0', 'velocity dispersion_total (2)', '0', 'tangential velocity dispersion_total (0)', '0', 'tangential velocity dispersion_total (1)', '0', 'tangential velocity dispersion_total (2)', '0', 'radial velocity dispersion_total (0)', '0', 'radial velocity dispersion_total (1)', '0', 'radial velocity dispersion_total (2)', '0', '', '1', 'bin radius_total', '0', 'number in bin_total', '0', 'number in bin_cumulative', '0', 'global id_total', '0', 'global id_average', '0', 'global id_cumulative', '0', 'potential_total', '0', 'potential_average', '0', 'potential_cumulative', '0', 'mass_total', '0', 'mass_average', '0', 'mass_cumulative', '0', 'eps_total', '0', 'eps_average', '0', 'eps_cumulative', '0', 'rho_total', '0', 'rho_average', '0', 'rho_cumulative', '0', 'hsmooth_total', '0', 'hsmooth_average', '0', 'hsmooth_cumulative', '0', 'temperature_total', '0', 'temperature_average', '0', 'temperature_cumulative', '0', 'metals_total', '0', 'metals_average', '0', 'metals_cumulative', '0', 'tform_total', '0', 'tform_average', '0', 'tform_cumulative', '0', 'velocity_total (Magnitude)', '0', 'velocity_average (Magnitude)', '0', 'velocity_cumulative (Magnitude)', '0', 'GlyphVector_total (Magnitude)', '0', 'GlyphVector_average (Magnitude)', '0', 'GlyphVector_cumulative (Magnitude)', '0', 'angular momentum_average (Magnitude)', '0', 'radial velocity_average (Magnitude)', '0', 'tangential velocity_average (Magnitude)', '0', 'velocity squared_average', '0', 'radial velocity squared_average', '0', 'tangential velocity squared_average', '0', 'velocity dispersion_total (Magnitude)', '0', 'tangential velocity dispersion_total (Magnitude)', '0', 'radial velocity dispersion_total (Magnitude)', '0', 'circular velocity_total', '0', 'density_total', '1']
DataRepresentation13.UseIndexForXAxis = 0

Render()
