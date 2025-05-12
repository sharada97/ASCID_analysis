from ASCID_variable import *
from ASCID_fit_shortcut import *
from ASCID_analysis_extended import *
# from ASCID_analysis import *



files = glob.glob('%s/RAW/*.CSV'%read_folder)
data, Data_Time = read_data(files)
# energyPlot = data_preprocess(data)
# save_fit(energyPlot)
# print("Calibration parameter found")
# energyPlot = []

data = data_calibration(data, energy_source)
analysis = ASCID_analysis(data, Data_Time)
analysis.data_preprocess()
analysis.coincidence()
analysis.CoaddedEnergy()
# analysis.exclusive_search()
analysis.added_energy("All", plot = plottype)
analysis.added_energy("inner", plot = plottype)
analysis.multiplicity(plot = plottype, veto = "Veto")
analysis.multiplicity(plot = plottype, veto = "noVeto")
analysis.coincidence_1det(plot = plottype) 
analysis.coincidence_2det(plot = plottype)
analysis.save_pickle()
