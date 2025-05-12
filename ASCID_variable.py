import glob
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from ASCID_functions import *

N_det_CsI = 26
N_det_panel = 11
N_det=N_det_CsI+N_det_panel
csiCh = np.arange(26)
innerCh=np.arange(16,26,1)
VetoCh = np.arange(41, 52, 1)


Det_Comb=N_det_CsI
energyThreshold_CsI = 100  #int(sys.argv[4])*10    #    #in keV
energyThreshold_Veto = 8 #in ADC, equivalent ~120 keV
energyThreshold = [energyThreshold_CsI, energyThreshold_Veto]
timeWindow = 5e-7   #int(sys.argv[3]) * 1e-7     #

SOURCE = ["Sc46", "Na22", "Background", "Cs137", "Ba133"]
Energy_Source = [889, 1275, 1460, 667, 356]
Vague_CalFactor = [1.8, 3.2, 4, 1.5, 1]
source_index = 4

source = SOURCE[source_index]  #"Na22", "Background"
energy_source = Energy_Source[source_index]  # in keV
Date="102324"
cut_data=900             #Select the required second of data

run_index = 5       # for the current set-up
folder_index = int(sys.argv[1])
index= int(sys.argv[2])

folder_all=glob.glob("../RUN5/%s/DAQ/CsI*Take4_*"%(source))                   #reading data folder
Total_folder = len(folder_all)
size_threshold = 0 #4e9  # Folders larger than 1GB
read_folder_all = filter_folders_by_size(folder_all, size_threshold)

read_folder_all.sort(key=lambda x: os.stat(x).st_ctime)
print(f"Total full data files yet:{len(read_folder_all)}")
# print(read_folder_all)

save_folder="../RUN5/%s/Figures/Take4"%source                                # choose plot save folder
# save_folder="../RUN5/Background/Figures/Take4"  

print(f"Folder name {read_folder_all[folder_index]}")

selected_folder = read_folder_all[folder_index]
folder_size = sum(
    os.path.getsize(os.path.join(dirpath, file))
    for dirpath, _, files in os.walk(selected_folder)
    for file in files
)
print(f"Actual memory of the selected folder: {folder_size / 1e9:.2f} G")





read_folder = read_folder_all[folder_index]
save_folder_each = os.path.join(save_folder, "Folder_%d"%folder_index)
try:
    os.makedirs(save_folder_each, exist_ok=False)
except:
    pass

sub_directory = "pulse_%d"%index
full_path = os.path.join(save_folder_each, sub_directory)
try:
    os.makedirs(full_path, exist_ok=False)
except:
    pass




date="%s"%(Date)
serial=np.array(["1st","2nd","3rd","4th","5th","6th", "7th", "8th"])




#histograms
energy_bin = 20
energy_end = 1000
plottype = True


def label(x,y,a,b):
    plt.xlabel(x,fontsize=18)
    plt.ylabel(y,fontsize=18)
    plt.rcParams['figure.figsize'] = [a,b]
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    return 
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],int(y[i]),fontsize=12)