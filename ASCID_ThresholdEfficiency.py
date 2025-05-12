import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
import os
import csv
import pandas as pd
from scipy.optimize import curve_fit
import glob
import sys
import pickle
from matplotlib.colors import LogNorm
from ASCID_fit_shortcut import read_data, find_1peaks

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

def gaus_fit(ri,rf,fx,fy,nm="norm"):
    from scipy.optimize import curve_fit
    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    from scipy.stats import chisquare
    b=len(bins)-1
    hist_PI=[]
    for j in range(b):
        if ri<=bins[j]<=rf:                                                          #Bins range for gauss fit
            hist_PI.extend([j])
    x=bins[hist_PI]+(bins[1]-bins[0])/2
    y=n[hist_PI]
    l = len(x)
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))#the number of data
    popt,pcov = curve_fit(gaus,x,y, p0=[max(y), mean, sigma])
    m=x[1]-x[0]
    X=np.arange(x.min(),x.max(),m/10)
#     chi2=chisquare(n[hist_PI],f_exp=gaus(x,*popt))
#     b=chi2[0]/(len(hist_PI)-3)
    if nm=="norm":
        textstr = '\n'.join((
        r'$\mu=%.2f \pm %.2f$ ADC' % (popt[1], np.sqrt(pcov[1,1])),
        r'$\sigma=%.2f \pm %.2f$ ADC' % (popt[2], np.sqrt(pcov[2,2]))
        ))               #"r'$\chi^2/Dof=%.2f$' % (b, )"
    elif nm=="keV":
        textstr = '\n'.join((
            r'$\mu=%.2f \pm %.1f$ keV' % (popt[1], np.sqrt(pcov[1,1])),
            r'$\sigma=%.2f \pm %.1f$ keV' % (popt[2], np.sqrt(pcov[2,2]))
            )) 
    elif nm=="keV_br":
        textstr = '\n'.join((
            r'$\mu=%.2f \pm %.2f$ eV' % (popt[1]*1000, np.sqrt(pcov[1,1])*1000),
            r'$\sigma=%.2f \pm %.2f$ eV' % (popt[2]*1000, np.sqrt(pcov[2,2])*1000)
            ))               #"r'$\chi^2/Dof=%.2f$' % (b, )"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(fx, fy, textstr, fontsize=12,
             
             verticalalignment='top', bbox=props)
    return popt[0], popt[1],popt[2],plt.plot(X,gaus(X,*popt),'r')




N_det=26+11
Det_Comb=21

energyThreshold = 100
timeWindow=5e-7
Date="103024"
         #Not required if you are doing the full run. 
#first cut_data event or 2nd, Ex: index=0 means collecting first 60 second of data and 
                                 #index=1 means collecting second 60 second of data and so on...     
# run_index=0      #int(sys.argv[1])    #Type the index of run mentioned in the data file, required while taking multiple dataset 

low_Th = 0
activity = 1.47
cut_data= 900               #Select the required second of data
source = 'Ba133'
energy_source = 356
new_drive = '/mnt/e'
read_folder_all=glob.glob("%s/RUN5/Take5/Ba133/DAQ/CsI_Ba1*"%new_drive)                   #reading data folder
# print(read_folder_all)
Total_folder = len(read_folder_all)
read_folder_all.sort(key=lambda x: os.stat(x).st_ctime)

folder_index = int(sys.argv[1])
read_folder= read_folder_all[folder_index]            #glob.glob("../RUN4/Sc46/DAQ/*")
save_folder="%s/RUN5/Take5/Ba133/Figures"%new_drive

if low_Th == 1:
    save_folder="%s/RUN5/Take5/Ba133/Figures/low_Th"%new_drive

print(f"Folder name {read_folder_all[folder_index]}")
index= int(sys.argv[2])

innerChannel=np.arange(16,26)

date="%s"%(Date)
serial=np.array(["1st","2nd","3rd","4th","5th","6th", "7th", "8th"])
save_folder_each = os.path.join(save_folder, "Folder_%d"%folder_index)
try:
    os.makedirs(save_folder_each, exist_ok=False)
except:
    pass






files = glob.glob('%s/RAW/*.CSV'%read_folder)
files=natsorted(files)[0:26]

files = glob.glob('%s/RAW/*.CSV'%read_folder)
files=natsorted(files)[0:26]
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
            data[0].append(int(row[1]))
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

# data[1] = []

# print(int(data[0][0]))
# data0_list = list(data[0])  # Convert to list
nEvents = len(data[0])
energyPlot=[[] for _ in range(26)]
for i, value in enumerate(data[0]):
    if int(value) < 26:
        energyPlot[int(value)].extend([data[2][i]])
# data0_list = None


mean = np.array([30.93, 45.15,47.66, 38.26, 52.63,     46.61, 40.73, 75.99, 39.51, 56.96,      32.3, 82.92, 74.96, 26.51, 90.15,
                25.1, 47.35, 44.57, 43.99, 53.14,     41.35, 42.5, 49.62, 24.81, 76.86, 105.77])
Mean = np.zeros(26)




for k in range(26):
    bin_last=mean[k]*2
    bin_start = 4 if mean[k] < 40 else 10
    bin_w = 1 if mean[k] < 40 else 2
    n,bins=np.histogram(energyPlot[k], bins=np.arange(bin_start,bin_last,bin_w))
    bin_centers = bins[1:]

    window_length = int(mean[k]/6)   #length of the window this and that side of mean
    # Mean[k, :] = find_peaks(bin_centers, n, k, window_length, min_distance, full_path)
    Mean[k] = find_1peaks(bin_centers, n, k, window_length, save_folder_each)
    
    print(f"fit channel {k} done, mean value found {Mean[k], Mean[k]}, given {mean[k], mean[k]}")


# nEvents = len(data[0])
# for j in range(nEvents):
#     if data[0][j]<26:
#         data[2][j]=data[2][j]/Mean[int(data[0][j])]*energy_source
# print("Calibration done")


# Count_raw = np.zeros((200, 26))
# for j in  range(nEvents):
#     for ind in range(200):
#         if ind*10 <= data[2][j]< (ind+1)*10 and data[0][j]<30:
#             Count_raw[ind, data[0][j]]+= 1



# f = open('%s/Efficiency_count_%d_%d.pickle'%(save_folder_each, folder_index, index), 'wb')
# Parameter={'Count':Count_raw, "Energy_range":np.arange(0, 2000, 10) }
# pickle.dump(Parameter, f)
# f.close()


f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'wb')
Parameter={'Time': cut_data, 'Mean':Mean}
pickle.dump(Parameter, f)
f.close()










