import csv
import numpy as np
# import xlsxwriter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from natsort import natsorted
import sys
#from tabulate import tabulate
import os
from matplotlib.animation import FuncAnimation
# from xlsxwriter import Workbook
# get_ipython().run_line_magic('matplotlib', 'inline')
import glob
import pickle
import shutil
# importing required libraries
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

N_det = 13
Det_Comb=21

energyThreshold = 100
timeWindow=5e-7
         #Not required if you are doing the full run. 
#first cut_data event or 2nd, Ex: index=0 means collecting first 60 second of data and 
                                 #index=1 means collecting second 60 second of data and so on...     
# run_index=0      #int(sys.argv[1])    #Type the index of run mentioned in the data file, required while taking multiple dataset 


activity = 0.6
cut_data=300             #Select the required second of data
read_folder_all=glob.glob("../ASCID_efficiency/DAQ/CsI_Na22_Containment_test_2inchLead*")                   #reading data folder
Total_folder = len(read_folder_all)
read_folder_all.sort(key=lambda x: os.stat(x).st_ctime)
folder_index = int(sys.argv[1])
read_folder= read_folder_all[folder_index]            #glob.glob("../RUN4/Sc46/DAQ/*")
save_folder="..//ASCID_efficiency/Figures/Pb_2inch" 

print(f"Folder name {read_folder_all[folder_index]}")
index= int(sys.argv[2])

innerChannel=np.arange(16,26)


serial=np.array(["1st","2nd","3rd","4th","5th","6th", "7th", "8th"])
save_folder_each = os.path.join(save_folder, "Folder_%d"%folder_index)
try:
    os.makedirs(save_folder_each, exist_ok=False)
except:
    pass


files = glob.glob('%s/RAW/*.CSV'%read_folder)
files=natsorted(files)#[0:26]
# files = files
# print(files)
MinTime=np.zeros(N_det)
for nf,nfile in enumerate(files):
    # nfile1 = os.path.join(folder_path, nfile)
    file = open(nfile)
    csvreader = csv.reader(file, delimiter = ';')
    header = next(csvreader)
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            MinTime[nf] = float(row[2]) * 1e-12
            break
minTime=np.min(MinTime)


# In[5]:

data=[]
for j in range(3):                                     #3=channel, time, calib. energy
    b=[]
    data.append(b)
nEvent=np.zeros(len(files))
for nf,nfile in enumerate(files):
    # nfile1 = os.path.join(folder_path, nfile)
    file = open(nfile)
    csvreader = csv.reader(file, delimiter = ';')
    header = next(csvreader)
    line_count = 0
    for row in csvreader:
        if float(row[2])*1e-12 - minTime>=cut_data*index:
            data[0].append(int(row[1])-26)
            data[1].append(float(row[2])*1e-12 - minTime)
            data[2].append(float(row[3]))     # row[4] = Caliberated energy, row[3] = ADC energy sometime, else row[3] = Caliberated energy
            line_count += 1
            # if line_count%1000000==0:
            #     print("%d M data loaded"%int(line_count/1000000))
            if float(row[2])*1e-12 - minTime>cut_data*(index+1):
                data[0].pop()
                data[1].pop()
                data[2].pop()
                break
    nEvent[nf]=line_count-1
    print("file %d done"%nf)
nEvents=int(nEvent.sum())


Data_Time=round((data[1][nEvents-1]-data[1][0]),0)
print("Total Time for this data set is %d s"%Data_Time)



data=np.array(data).T
data=data[data[:, 1].argsort()]
data=np.array(data).T


mean = np.array([466, 791, 820, 962, 950, 776, 693, 658, 616, 600, 754, 503, 700])
Cal_keV = np.ones(N_det)*511




nEvents = len(data[0])
for j in range(nEvents):
    data[2][j]=data[2][j]/mean[int(data[0][j])]*Cal_keV[int(data[0][j])]


nEvents = len(data[0])
csiCh = np.arange(N_det)
timeWindow = 5e-7
Sig_index=-1
eventTime=0
Timed_Energy_event=[]
for j in range(nEvents):
    if  data[1][j]>eventTime+2*timeWindow and int(data[0][j]) in csiCh and data[2][j] > 100:   #Avoiding center detector overlap
        eventTime=data[1][j]
        eventChannel=int(data[0][j])
        Ev_indx=j
        Sig_index+=1
        EventEnergy=data[2][j]
        a=[]
        a.extend([Ev_indx])
        
        while Ev_indx>0 and abs(data[1][Ev_indx-1]-eventTime)<=timeWindow:
            Ev_indx-=1
            a.extend([Ev_indx])
        Ev_indx=j
        
        while Ev_indx<nEvents-1 and abs(data[1][Ev_indx+1]-eventTime)<=timeWindow:
            Ev_indx+=1    
            a.extend([Ev_indx])
        a = np.unique(np.sort(a)) # a contains the index within the defined coincidence window
        Timed_Energy_event.extend([a])
    if j %1000000 == 0:
        print("%d M done"%int(j/1000000))
# Timed_Energy_event = Timed_Energy_event
# Timed_Energy_event=[]
print("data coincidence done")



data = np.array(data)
sl=[]
count = 0
Multiplicity = np.zeros(30)

for k,j in enumerate(Timed_Energy_event):
    event_Channel = data[0][j]
    event_Energies = data[2][j]
    multiplicity = len(j)
    Multiplicity[multiplicity]+=1


def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i]*1e4,format(y[i], '.1e'),fontsize=12)


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
values = count[1:N_det]/1e4
addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
plt.title("Multiplicity Distribution - Na22", fontsize="15")
plt.savefig("%s/Multiplicity_%d.jpg"%(save_folder_each, index))
plt.close()



Coincident_energy = np.zeros((len(Timed_Energy_event), N_det), dtype = float)
for k,j in enumerate(Timed_Energy_event):
    event_Channel = data[0][j]
    event_Energies = data[2][j]
    # Coincidence_energy
    if 0 in event_Channel:
        for l,i in enumerate(event_Channel):
            Coincident_energy[k, int(i)] = event_Energies[l]

triggered = Coincident_energy[Coincident_energy[:, 0]!=0]


count_pair = np.zeros(N_det)
for i in range(1, N_det):
    filtered_rows = triggered[triggered[:, i] != 0]
    selected_columns = filtered_rows[:, [0, i]]
    for j,k in enumerate(selected_columns):
        if 511 - 2*34 <k[0]< 511 + 2*34 and 511 - 2*34 <k[1]< 511 + 2*34:
            count_pair[i] += 1


multiplicity_coin0 = np.zeros(N_det)
for k,j in enumerate(triggered):
    filtered_rows = j[np.nonzero(j)]
    if 511 - 2*34 <j[0]< 511 + 2*34:
        multiplicity_coin0[len(filtered_rows)-1] += 1



count = multiplicity_coin0
total = count.sum()
count=count/len(triggered)*1000000
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
values = count[1:N_det]/1e4
addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
plt.title("Multiplicity Distribution - Na22, Ch-0: 511 keV triggered", fontsize="15")
plt.savefig("%s/Multiplicity_Coin0_%d.jpg"%(save_folder_each, index))
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
values = count[1:N_det]/1e4
addlabels(courses, values)
plt.legend(fontsize=14)
plt.yscale("log")
label("Index of the detector", "No. of events (percentage)", 8, 5)
plt.title("Detector wise number of 511 keV triggered trend, Ch-0: 511 keV triggered", fontsize="15")
plt.savefig("%s/Event_Detector_%d.jpg"%(save_folder_each, index))
plt.close()




f = open('%s/Output_Run_%d_%d.pickle'%(save_folder_each, folder_index, index), 'wb')
Parameter={'Multiplicity':Multiplicity,'multiplicity_coin0':multiplicity_coin0, 'count_pair':count_pair}
pickle.dump(Parameter, f)
f.close()

