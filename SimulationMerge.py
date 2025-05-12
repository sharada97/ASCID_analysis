import csv
import numpy as np
# import xlsxwriter
import matplotlib.pyplot as plt
#from tabulate import tabulate
import os
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib.animation import FuncAnimation
import glob
import pickle
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
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],int(y[i]),fontsize=12)

N_det=26
Det_Comb=21
energyThreshold = 100
hist_endEnergy = 2500
gamma_source = "Sc46"   #Background, Na22
serial=np.array(["1st", "2nd", "3rd","4th","5th","6th"])
read_folder="/data7/CsI_Project/Simulation/Sc46/Figure/Smear/Each"
save_folder="/data7/CsI_Project/Simulation/Sc46/Figure/Smear"
index=0




##Inner and Outer channel combined spectra
n=[]
index=[5, 6, 7]
nDecay=0
for i in index:
    f = open('%s/Output_Run-folder%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    nEvents=dat['nEvent']
    nDecay=nDecay+nEvents
    n1=dat["n_EnergyInner"]
    bins=dat["bins_combined"]
    n.extend([n1])
    f.close()

print(nDecay)    




n_inner=np.zeros(len(n1))
for j in range(len(index)):
    n_inner+=n[j]
plt.rcParams['figure.figsize'] = [12,5]    
plt.step(bins[:-1], n_inner, 'k',where="mid",linewidth=1,label="%d events"%n_inner.sum())
plt.legend(fontsize=15)
label("Energy (keV)", "Counts", 12, 5)
plt.title("Inner detector coincidence Spectrum, combined data (%d Decays) "%nDecay, fontsize=14)
plt.savefig("%s/Energy-InnCoin.jpg"%(save_folder))
plt.close()



n=[]
for i in index:
    f = open('%s/Output_Run-folder%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    n1=dat["n_EnergyAll"]
    bins=dat["bins_combined"]
    n.extend([n1])
    f.close()
    
    
n_All=np.zeros(len(n1))
for j in range(len(index)):
    n_All+=n[j]
    
plt.step(bins[:-1], n_All, 'k',where="mid",linewidth=1,label="%d events"%n_All.sum())
plt.legend(fontsize=15)
label("Energy (keV)", "Counts", 12, 5)
plt.title("All detector coincidence Spectrum, combined data (%d Decays) "%nDecay, fontsize=14)
plt.savefig("%s/Energy-AllCoin.jpg"%(save_folder))
plt.close()
print("Coincidence energy spectrum done")
bins_combined = bins



#Multiplicity plot
count=np.zeros(N_det)   
for i in index:
    f = open('%s/Output_Run-folder%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    count_mult=dat["Multiplicity"]
#     print(count_mult)
    count=count+count_mult
    f.close()


data_count = {'1 Det':count[1], '2Det':count[2], '3 Det':count[3],
        '4 Det':count[4], '5 Det':count[5], '6 Det':count[6], '7 Det':count[7],
             '8 Det':count[8], '9 Det':count[9], '10 Det':count[10]}
courses = list(data_count.keys())
values = list(data_count.values())
  
fig = plt.figure(figsize = (12, 5))
 
plt.bar(courses, values,
        width = 0.4, label="Total decay events= %d decays"%(count.sum()))   #Change time here
addlabels(courses, count[1:])
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 12, 5) 
plt.title("Multiplicity Distribution-%d keV Threshold (%d decays- %s)"%( energyThreshold, nDecay, gamma_source), fontsize="15")
plt.savefig("%s/Multiplicity_%dkeV.jpg"%(save_folder,energyThreshold))
plt.close()




count=np.array(count, dtype=int)
nEvent=count.sum()
count=count/count.sum()*100
count=np.array(count, dtype="float")
for j in range(len(count)):
    count[j]=round(count[j],0)

data_count = {'1 Det':count[1], '2 Det':count[2], '3 Det':count[3],
        '4 Det':count[4], '5 Det':count[5], '6 Det':count[6]}
courses = list(data_count.keys())
values = list(data_count.values())
fig = plt.figure(figsize = (12, 5))
 
plt.bar(courses, values,
        width = 0.4, label="Total events= %d "%(nDecay))   #Change time here
addlabels(courses, count[1:])
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("Multiplicity Distribution in Percentage; Threshold-%dkeV (%d Decays- %s)"%(energyThreshold, nDecay, gamma_source), fontsize="15")
plt.savefig("%s/Multiplicity_Percent-%dkeV.jpg"%(save_folder,energyThreshold))    
plt.close()
print("Multiplicity done")



##1- Det coincidence

n_1det=[]
for i in index:
    f = open('%s/Output_Run-folder%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    n=dat["n_1det"]
    bins=dat["bins1_keV"]
    n_1det.extend([n])
    f.close()

n1=np.zeros(len(n))
for j in range(len(index)):
    n1+=n_1det[j]

plt.step(bins[:-1], n1, 'k',where="mid",linewidth=1,label="%d events"%n1.sum())
plt.legend(fontsize=15)
label("Energy (keV)", "Counts", 12, 5)
plt.title("1 detector: Spectrum, combined data (%d Decays)"%(nDecay), fontsize=14)
plt.savefig("%s/Energy-1Det.jpg"%(save_folder))
plt.close()
bins1_keV = bins




# 
count1_ch=np.zeros(26)
for i in index:
    f = open('%s/Output_Run-folder%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    count1_ch+=dat["count1det_ch"]
    f.close()
    
data_count = {'0':count1_ch[0], '1':count1_ch[1], '2':count1_ch[2],
        '3':count1_ch[3], '4':count1_ch[4], '5':count1_ch[5], '6':count1_ch[6],
             '7':count1_ch[7], '8':count1_ch[8], '10':count1_ch[9], '11':count1_ch[10],
             '11':count1_ch[11], '12':count1_ch[12], '13':count1_ch[13], '14':count1_ch[14],
             '15':count1_ch[15], '16':count1_ch[16], '17':count1_ch[17], '18':count1_ch[18],
             '19':count1_ch[19], '20':count1_ch[20], '21':count1_ch[21], '22':count1_ch[22],
             '23':count1_ch[23], '24':count1_ch[24], '25':count1_ch[25]}
courses = list(data_count.keys())
values = list(data_count.values())

fig = plt.figure(figsize = (10, 5))

plt.bar(courses, values,
        width = 0.4, label="Total events= %d "%int(count1_ch.sum()))   #Change time here
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 10, 5)
plt.title("1Det coincidence events, combined data (%d Decays)"%(nDecay), fontsize="15")
plt.savefig("%s/CountEvents-1det_5-5.jpg"%(save_folder))
plt.close()
print("1 Detector Coincidence done")


"""
# 2-det

Energy_2det=[[],[]]
for i in index:
    f = open('%s/Output_Run%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    Energy_2det1=dat["Energy_2det"]
    Energy_2det[0].extend(Energy_2det1[0])
    Energy_2det[1].extend(Energy_2det1[1])
    f.close()

   
N=300
h=plt.hist2d(Energy_2det[0], Energy_2det[1], bins=(np.arange(0,2000,2000/N), np.arange(0,2000,2000/N)), norm = LogNorm())
# h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("2 detector coincidence events, %d events, combined data (%d Decay) "%(len(Energy_2det[0]),nDecay), fontsize=16)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 10,5)
plt.savefig("%s/2Det-DistribuionPlot.jpg"%(save_folder))
plt.close()
print("2 Detector Coincidence done")

# 3 - det

Energy_3det_sl=[[],[]]
for i in index:
    f = open('%s/Output_Run%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    Energy_3det1=dat["Energy_3det_sl"]
    Energy_3det_sl[0].extend(Energy_3det1[:,0])
    Energy_3det_sl[1].extend(Energy_3det1[:,1])
    f.close()
    
N=300
h=plt.hist2d(Energy_3det_sl[0], Energy_3det_sl[1], bins=(np.arange(0,1200,2000/N), np.arange(0,1200,2000/N)), norm = LogNorm())
# h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("3 detector coincidence events, %d events, combined data (%d Decays)"%(len(Energy_3det_sl[0]), nDecay), fontsize=16)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 10,5)
plt.savefig("%s/3Det-DistribuionPlot.jpg"%(save_folder))
plt.close()
print("3 Detector Coincidence done")

# 4-det

## shared-511 keV

Energy_4det_sl511=[[],[]]
for i in index:
    f = open('%s/Output_Run%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    Energy_4det1=dat["Energy_4det_sl511"]
    Energy_4det_sl511[0].extend(Energy_4det1[:,0])
    Energy_4det_sl511[1].extend(Energy_4det1[:,1])
    f.close()

    
N=300
h=plt.hist2d(Energy_4det_sl511[0], Energy_4det_sl511[1], bins=(np.arange(0,800,1300/N), np.arange(0,800,1300/N)), norm = LogNorm())
# h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("4 detector, Other two hits - [1275, 511], %d events, combined data (%d Decays)"%(len(Energy_4det_sl511[0]), nDecay), fontsize=16)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 10,5)
plt.savefig("%s/4Det-DistribuionPlot_511.jpg"%(save_folder))
plt.close()



## shared-1275 keV
Energy_4det_sl1275=[[],[]]
for i in index:
    f = open('%s/Output_Run%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    Energy_4det1=dat["Energy_4det_sl1275"]
    Energy_4det_sl1275[0].extend(Energy_4det1[:,0])
    Energy_4det_sl1275[1].extend(Energy_4det1[:,1])
    f.close()

    
N=300
h=plt.hist2d(Energy_4det_sl1275[0], Energy_4det_sl1275[1], bins=(np.arange(0,1500,1300/N), np.arange(0,1500,1300/N)), norm = LogNorm())
# h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("4 detector, Other two hits - [511, 511], %d events, combined data (%d Decays)"%(len(Energy_4det_sl1275[0]), nDecay), fontsize=16)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 10,5)
plt.savefig("%s/4Det-DistribuionPlot_1275.jpg"%(save_folder))
plt.close()
print("4 Detector Coincidence done")


#5-det

Energy_5det_sl=[[],[]]
for i in index:
    f = open('%s/Output_Run%d.pickle'%(read_folder,i), 'rb')
    dat=pickle.load(f)
    Energy_5det1=dat["Energy_5det_sl"]
    Energy_5det_sl[0].extend(Energy_5det1[:,0])
    Energy_5det_sl[1].extend(Energy_5det1[:,1])
    f.close()

    
N=300
h=plt.hist2d(Energy_5det_sl[0], Energy_5det_sl[1], bins=(np.arange(0,1200,2000/N), np.arange(0,1200,2000/N)), norm = LogNorm())
plt.title("5 detector, other 3-[1275], %d events, combined data (%d Decays)"%(len(Energy_5det_sl[0]), nDecay), fontsize=16)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 10,5)
plt.savefig("%s/5Det-DistribuionPlot.jpg"%(save_folder))
plt.close()
print("5 Detector Coincidence done")
"""

f = open('%s/Result_Run.pickle'%(save_folder), 'wb')
Parameter={'nEvent': nDecay,  "Multiplicity":count_mult, 
           'n_EnergyInner':n_inner, 'n_EnergyAll':n_All, 'bins_combined':bins_combined,
            "n_1det":n1,"bins1_keV":bins1_keV, "count1det_ch":count1_ch}
pickle.dump(Parameter, f)
f.close()

