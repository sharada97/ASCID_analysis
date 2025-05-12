from ASCID_variable import *
from ASCID_fit_shortcut import *
from ASCID_functions import *
from ASCID_analysis_extended import *




files = glob.glob('%s/RAW/*.CSV'%read_folder)
data, Data_Time = read_data(files)
energyPlot = data_preprocess(data)
save_fit(energyPlot)
print("Calibration parameter found")
energyPlot = []
data = data_calibration(data, energy_source)


for time_window in range(7,11,2):
    timeWindow = time_window*1e-7
    for energyThreshold_CsI in range(60, 140, 20):
        analysis = ASCID_analysis(data, Data_Time, energyThreshold_CsI, timeWindow)
        analysis.data_preprocess()
        analysis.coincidence()
        analysis.CoaddedEnergy()
        analysis.exclusive_search()
        analysis.added_energy("All", plot = plottype)
        analysis.added_energy("inner", plot = plottype)
        analysis.multiplicity(plot = plottype, veto = "Veto")
        analysis.multiplicity(plot = plottype, veto = "noVeto")
        time.sleep(60)
        print("Energy threshold done for %d keV"%(energyThreshold_CsI))
        print()
    print("Time coincidence done for %d ns"%(timeWindow*1e9))


