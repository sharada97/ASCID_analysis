import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
import os
import pandas as pd
from scipy.optimize import curve_fit
import glob
import pickle
from matplotlib.colors import LogNorm

NDet = 13

def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i]*1e4,format(y[i], '.2e'),fontsize=12)
def label(x,y,a,b,pur="all"):
    plt.xlabel(x,fontsize=18)
    plt.ylabel(y,fontsize=18)
    if pur=="all":
        plt.rcParams['figure.figsize'] = [a,b]
    elif pur=="chi2":
        plt.title("Chi2 Analysis, S-%s, Channel %s"%(a,b), fontsize=20)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    return 

def find_files_partial_match_in_folders(partial_name, search_path):
    """
    Search for files with partially matching names in all subdirectories.

    :param partial_name: Part of the filename to search for.
    :param search_path: Root directory to start the search.
    :return: List of full paths to matching files.
    """
    file_paths = []
    for root, dirs, files in os.walk(search_path):  # Traverse all subdirectories
        for file in files:
            if partial_name in file:  # Check for partial match
                file_paths.append(os.path.join(root, file))  # Save the full path
    return file_paths


class DataMerge:
    def __init__(self, files):
        self.files = files

        self.count_pair = None
        self.multiplicity_coin0 = None
        self.Multiplicity = None


    def label(self, x,y,a,b):
        plt.xlabel(x,fontsize=18)
        plt.ylabel(y,fontsize=18)
        plt.rcParams['figure.figsize'] = [a,b]
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)
        return 

        
    def load_data(self):
        # all_files = len(self.files)
        # for i in range(all_files):
        file = self.files[0]
        f = open('%s'%(file), 'rb')
        dat=pickle.load(f)
        self.count_pair=dat['count_pair']
        self.multiplicity_coin0=dat['multiplicity_coin0']
        self.Multiplicity = dat['Multiplicity']
        f.close()

    def merging(self):
        all_files = len(self.files)
        for i in range(1, all_files):
            f = open('%s'%(self.files[i]), 'rb')
            dat=pickle.load(f)
            self.count_pair+= dat['count_pair']
            self.multiplicity_coin0+=dat['multiplicity_coin0']
            self.Multiplicity+= dat['Multiplicity']
            f.close()
        return self.Multiplicity, self.multiplicity_coin0, self.count_pair
            


save_folder = "../ASCID_efficiency/Figures/Pb_2inch/"
search_path = "%s"%save_folder
partial_name = "Output_Run"
save_folder_each = "%s/Result"%save_folder

file_array = find_files_partial_match_in_folders(partial_name, search_path)
print(f"Files found: {len(file_array)}")

# files = glob.glob("%s/Output_Run5*.pickle"%save_folder) 
files = file_array
data_merge = DataMerge(files)
data_merge.load_data()
Multiplicity, multiplicity_coin0, count_pair = data_merge.merging()




cut_data = 300*len(file_array)
count = Multiplicity
total = count.sum()
count=count/count.sum()*1000000
count=np.array(count, dtype="float")
for j in range(len(count)):
    count[j]=round(count[j],2)
data_count = {'1':count[1], '2':count[2], '3':count[3],
        '4':count[4], '5':count[5], '6':count[6],
             '7':count[7], '8':count[8], '9':count[9],
        '10':count[10], '11':count[11], '12':count[12]}
courses = list(data_count.keys())
values = list(data_count.values())
fig = plt.figure(figsize = (10, 5))
plt.bar(courses, values,
        width = 0.4, label="Total decay rate= %.1f decay/sec"%(total/(cut_data)))   #Change time here
values = count[1:NDet]/1e4
# addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
plt.title("Multiplicity Distribution - Na22", fontsize="15")
plt.savefig("%s/Multiplicity.jpg"%(save_folder_each))
plt.close()



count = multiplicity_coin0
total = count.sum()
count=count/total*1000000
count=np.array(count, dtype="float")
for j in range(len(count)):
    count[j]=round(count[j],2)
data_count = {'1':count[1], '2':count[2], '3':count[3],
        '4':count[4], '5':count[5], '6':count[6],
             '7':count[7], '8':count[8], '9':count[9],
        '10':count[10], '11':count[11], '12':count[12]}
courses = list(data_count.keys())
values = list(data_count.values())
fig = plt.figure(figsize = (10, 5))
plt.bar(courses, values,
        width = 0.4)   #Change time here
values = count[1:NDet]/1e4
addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
plt.title("Multiplicity Distribution - Na22, Ch-0: 511 keV triggered", fontsize="15")
plt.savefig("%s/Multiplicity_Coin0.jpg"%(save_folder_each))
plt.close()





count = count_pair
total = count.sum()
count=count/total*1000000
count=np.array(count, dtype="float")
for j in range(len(count)):
    count[j]=round(count[j],2)
data_count = {'1st':count[1], '2nd':count[2], '3rd':count[3],
        '4th':count[4], '5th':count[5], '6th':count[6],
             '7th':count[7], '8th':count[8], '9th':count[9],
        '10th':count[10], '11th':count[11], '12th':count[12]}
courses = list(data_count.keys())
values = list(data_count.values())
# fig = plt.figure(figsize = (10, 5))
plt.bar(courses, values,
        width = 0.4)   #Change time here
values = count[1:NDet]/1e4
addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("Index of the detector", "No. of events (percentage)", 8, 5)
plt.title("Detector wise number of 511 keV triggered trend, Ch-0: 511 keV triggered", fontsize="15")
plt.savefig("%s/Event_Detector_all.jpg"%(save_folder_each))
plt.close()



f = open('%s/Result_Run.pickle'%(save_folder_each), 'wb')
Parameter={'Multiplicity':Multiplicity,'multiplicity_coin0':multiplicity_coin0, 'count_pair':count_pair}
pickle.dump(Parameter, f)
f.close()


