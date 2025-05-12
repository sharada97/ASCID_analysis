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


# In[2]:


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


def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],int(y[i]),fontsize=12)




N_det=26+11
Det_Comb=21

energyThreshold = 100
timeWindow=5e-7
Date="101724"
         #Not required if you are doing the full run. 
#first cut_data event or 2nd, Ex: index=0 means collecting first 60 second of data and 
                                 #index=1 means collecting second 60 second of data and so on...     
# run_index=0      #int(sys.argv[1])    #Type the index of run mentioned in the data file, required while taking multiple dataset 




cut_data=900               #Select the required second of data

read_folder="../RUN5/Sc46/DAQ/101724_CsI_Sc46_2/"                          #glob.glob("../RUN4/Sc46/DAQ/*")
save_folder="../RUN5/Sc46/Figures/%s"%(Date)
folder_index = 0   #int(sys.argv[1])
index=int(sys.argv[1])



innerChannel=np.arange(16,26)

date="%s"%(Date)
serial=np.array(["1st","2nd","3rd","4th","5th","6th", "7th", "8th"])






# In[4]:

# print("Loading DATA %s"%read_folder[folder_index])

# read_fold = read_folder[folder_index]
files = glob.glob('%s/RAW/*.csv'%read_folder)
# import os

# # Path to the folder
# folder_path = '%s/RAW/'%read_folder

# # List all files in the folder
# files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]


# print(files)
files=natsorted(files)
print("\n".join(map(str, files)))
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
            if line_count%1000000==0:
                print("%d M data loaded"%int(line_count/1000000))
            if float(row[2])*1e-12 - minTime>cut_data*(index+1):
                data[0].pop()
                data[1].pop()
                data[2].pop()
                break
    nEvent[nf]=line_count-1
    print("file %d done"%nf)
nEvents=int(nEvent.sum())
print("Total time of data taking is %.3f seconds"%data[1][nEvents-1])
# print("Total time of data taking is %.3f seconds and total events is %d"%(data[1][nEvents-1], nEvents))

Data_Time=round((data[1][nEvents-1]-data[1][0]),0)
print("Total Time for this data set is %d s"%Data_Time)



data=np.array(data).T
data=data[data[:, 1].argsort()]
data=np.array(data).T



print("Calibration starts")

f = open('%s/DataInfo_Run.pickle'%(save_folder), 'rb')        #Change directory for saved pickle file of mean and sigma
dat = pickle.load(f)
# mean=dat['Mean']
# sigma=dat['Sigma']
bin_last=dat['bin_last']
bin_w=dat['bin_width']
gaus_x=dat['gaus_x']
gaus_y=dat['gaus_y']
gaus_diff=dat['gaus_diff']
f.close()

# multiplier = folder_index*3 + index 
multiplier = 1
# no=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# gaus_diff=np.array([[1.6,0,0,1.2,1,      0,0,0,1.2,0,     0,0,0,0,0.1,     0.2,0,0,0,0,     0,0,0,0,0,0],
#                     [1.6,0,0,1.2,1,      0,0,0,0.8,0,     0,0,0,3.5,0,     1,0,0,2,0,     1.2,0,0,0,0,0]])


#gaus_diff = 0
#Calibraion of the data
energyPlot=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
for i in range(nEvents):
    energyPlot[int(data[0][i])].extend([data[2][i]])
    
n_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]   
# gaus_x=gaus_x-gaus_diff*index
# gaus_y=gaus_y-gaus_diff*index

mean=np.zeros((26,2))
sigma=np.zeros((26,2))
for k in range(26):
    n,bins,patches=plt.hist(energyPlot[k], bins=np.arange(15,bin_last[k],bin_w[k]), histtype="step",label="%devents"%len(energyPlot[k]))
    label("Energy (ADC)", "Counts", 8,5)
    plt.title(r"Channel %d Energy Spectrum- %d min data (Sc46)"%(k,int(Data_Time/60)), fontsize=16)
    a1, mean[k,0],sigma[k,0], fig=gaus_fit(gaus_x[0][k]-gaus_diff[0][k]*multiplier,gaus_x[1][k]-gaus_diff[0][k]*multiplier,bin_last[k]/3,max(n)-max(n)/10)
    a2, mean[k,1],sigma[k,1], fig=gaus_fit(gaus_y[0][k]-gaus_diff[1][k]*multiplier,gaus_y[1][k]-gaus_diff[1][k]*multiplier,bin_last[k]/1.5,max(n)-max(n)/2)
    n_energy[k]=n
    plt.legend(fontsize=14)
    plt.savefig("%s/pulse/Energy_Ch%d_%d_%d.jpg"%(save_folder,k, folder_index, index))
    plt.close()
#     print("%d done"%k)

    




for j in range(nEvents):
    data[2][j]=data[2][j]/mean[int(data[0][j]),1]*1132
    
    

   
for k in range(26):
    mean[k,0]=mean[k,0]/mean[k,1]*1132
    sigma[k,0]=sigma[k,0]/mean[k,1]*1132
    sigma[k,1]=sigma[k,1]/mean[k,1]*1132
    mean[k,1]=mean[k,1]/mean[k,1]*1132

    
    
data_process = np.array(data).T
energyThreshold_Veto = 8
Veto_threshold = np.ones(11)*energyThreshold_Veto
Threshold_CsI=np.ones(26)*energyThreshold
Threshold_null = np.zeros(15)
Threshold = np.append(np.append(Threshold_CsI,Threshold_null), Veto_threshold)
print(Threshold)
np.shape(data_process)
data_process1 = []
for j in range(nEvents):
    if data_process[j, 2] > Threshold[int(data_process[j, 0])]:
        data_process1.extend(data_process[j])





## Multiplicity distribution starts

print("********Multiplicity Distribution starts*******")
# for energyThreshold in Thre:               #Threshold check


count_rest=0
count_main=0
eventTime=0
count=np.zeros(Det_Comb)
Sig_index=-1
Timed_Energy_event=[]
EventEnergyAll=[]
EventEnergyInner=[]
innerCh=np.arange(16,26,1)
for j in range(nEvents):
    if  data[1][j]>eventTime+2*timeWindow and data[2][j]> Threshold[int(data[0][j])] :   #Avoiding center detector overlap
        eventTime=data[1][j]
        eventChannel=int(data[0][j])
        Ev_indx=j
        Sig_index+=1
        EventEnergy=data[2][j]
        a=[]
        while Ev_indx>0 and abs(data[1][Ev_indx-1]-eventTime)<=timeWindow:
            Ev_indx-=1
            a.extend([Ev_indx])
        Ev_indx=j

        while Ev_indx<nEvents-1 and abs(data[1][Ev_indx+1]-eventTime)<=timeWindow:
            Ev_indx+=1    
            a.extend([Ev_indx])
        np.sort(a) 
#         Timed_event.extend([a])        

        EventEnergyall=EventEnergy
        if eventChannel in innerCh:
            EventEnergyinn=EventEnergy
        else:
            EventEnergyinn=0
        
        TE_event=[]
        TE_event.extend([j])              
        signal=0
        for k in range(len(a)):
            if data[2][a[k]] >= Threshold[int(data[0][a[k]])]:   #Considering all double detector, same one also
                EventEnergyall+=data[2][a[k]]
                if data[0][a[k]] in innerCh:
                    EventEnergyinn+=data[2][a[k]]
                signal+=1
                TE_event.extend([a[k]])
        np.sort(TE_event)
        Timed_Energy_event.extend([TE_event])  
        EventEnergyAll.extend([EventEnergyall])
        if EventEnergyinn!=0:
            EventEnergyInner.extend([EventEnergyinn])
#         if any(data[0][sl] ==3 and 1081<data[2][sl] <1281 for sl in TE_event):      #Change energy min max.
#             Timed_Energy_Ch3_event.extend([TE_event])    


        for sig in range(Det_Comb):                              
            if signal == sig:                    
                count[sig]+=1                     
#                 count_index[sig].extend([Sig_index]) 

    if j%1000000==0:
        print("%d M data loaded"%int(j/1000000))    

count_mult=count     
    
#Checking the inner detector combination hits
Detect=0
for j in EventEnergyAll:
    if j>1500:
        Detect+=1
print("Including all Total detection is %d out of %d"%(Detect, len(EventEnergyAll)))


Detectn=0
for j in EventEnergyInner:
    if j>1500:
        Detectn+=1
print("Including inner Total detection is %d out of %d"%(Detectn, len(EventEnergyInner)))

#All and Inner detector coincidence energy 
plt.rcParams['figure.figsize'] = [12,5]
n,bins,patches=plt.hist(EventEnergyAll, bins=np.arange(0,4000,20),label="%d events"%len(EventEnergyAll), histtype="step")
plt.legend(fontsize=15)
plt.xticks(np.arange(0, 4000, 250))
label("Energy (keV)", "Counts", 12, 5)
plt.title("All detector coincidence Spectrum  ", fontsize=14)
plt.savefig("%s/Energy-AllCoin.jpg"%(save_folder))
plt.close()
n_EnergyAll=n

    
n,bins,patches=plt.hist(EventEnergyInner, bins=np.arange(0,4000,20),label="%d events"%len(EventEnergyInner), histtype="step")
plt.legend(fontsize=15)
plt.xticks(np.arange(0, 4000, 250))
label("Energy (keV)", "Counts", 12, 5)
plt.title("Inner detector coincidence Spectrum  ", fontsize=14)
plt.savefig("%s/Energy-InnCoin.jpg"%(save_folder))
plt.close()
n_EnergyInner=n 
bins_combined=bins
        
print(n_EnergyInner.sum(), n_EnergyAll.sum())        
        

# count=np.array(count, dtype=int) 
# count_mult=count
# np.savetxt("%s/count5_%dkeV.txt"%(save_folder, energyThreshold), count, newline="\n")
        
        
        
data_time=Data_Time
data_count = {'1 Det':count[0], '2 Det':count[1], '3 Det':count[2],
        '4 Det':count[3], '5 Det':count[4], '6 Det':count[5], '7 Det':count[6],
             '8 Det':count[7], '9 Det':count[8], '10 Det':count[9]}
courses = list(data_count.keys())
values = list(data_count.values())
fig = plt.figure(figsize = (10, 5))
plt.bar(courses, values,
        width = 0.4, label="Total decay rate= %.1f decay/sec"%(count.sum()/(data_time)))   #Change time here
addlabels(courses, values)
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("Multiplicity Distribution-%d keV (%d s- Sc46)"%( energyThreshold, data_time), fontsize="15")
plt.savefig("%s/Multiplicity_%dkeV.jpg"%(save_folder,energyThreshold))
plt.close()




count=np.array(count, dtype=int)
count=count/count.sum()*100
count=np.array(count, dtype="float")
for j in range(len(count)):
    count[j]=round(count[j],2)
data_time=Data_Time
data_count = {'1 Det':count[0], '2 Det':count[1], '3 Det':count[2],
        '4 Det':count[3], '5 Det':count[4], '6 Det':count[5]}
courses = list(data_count.keys())
values = list(data_count.values())
fig = plt.figure(figsize = (10, 5))
 
plt.bar(courses, values,
        width = 0.4, label="Total events= %d "%(nEvents))   #Change time here
addlabels(courses, count)
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("Multiplicity Distribution in Percentage; Threshold-%dkeV (%d s- Sc46)"%(energyThreshold,Data_Time), fontsize="15")
plt.savefig("%s/Multiplicity_Percent-%dkeV.jpg"%(save_folder,energyThreshold))    
plt.close()
    
        

        
        
        
     
print("1 detector Coincidence starts") 
Energy_1det=[]
Channel_1det=[]
Energy_1det_ch=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
Energy_1det_inner=[]
count1_ch=np.zeros(26)
for l in Timed_Energy_event:
    if len(l)==1:
        E=data[2][l[0]]
        p=int(data[0][l[0]])
        Energy_1det.extend([E])
        Channel_1det.extend([p])
        Energy_1det_ch[p].extend([E])
        count1_ch[p]+=1
        if p in innerChannel:
            Energy_1det_inner.extend([E])




count1_ch=np.array(count1_ch, dtype=int) 


data_time=Data_Time
data_count = {'0':count1_ch[0], '1':count1_ch[1], '2':count1_ch[2],
        '3':count1_ch[3], '4':count1_ch[4], '5':count1_ch[5], '6':count1_ch[6],
             '7':count1_ch[7], '8':count1_ch[8], '9':count1_ch[9], '10':count1_ch[10],
             '11':count1_ch[11], '12':count1_ch[12], '13':count1_ch[13], '14':count1_ch[14],
             '15':count1_ch[15], '16':count1_ch[16], '17':count1_ch[17], '18':count1_ch[18],
             '19':count1_ch[19], '20':count1_ch[20], '21':count1_ch[21], '22':count1_ch[22],
             '23':count1_ch[23], '24':count1_ch[24], '25':count1_ch[25]}
courses = list(data_count.keys())
values = list(data_count.values())

fig = plt.figure(figsize = (10, 5))

plt.bar(courses, values,
        width = 0.4, label="Total events= %.1f "%(count1_ch.sum()))   #Change time here
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("1Det coincidence events", fontsize="15")
plt.savefig("%s/CountEvents-1det_5-5.jpg"%(save_folder))
plt.close()


def exp(x,a,b):
    return a+b*np.exp(-x/400)
n,bins,patches=plt.hist(Energy_1det, bins=np.arange(0,2500,20), histtype="step",label="%d events"%len(Energy_1det))
plt.legend(fontsize=15)
label("Energy (keV)", "Counts", 8, 5)
plt.title("1 detector: Spectrum  ", fontsize=14)
plt.savefig("%s/Energy-1Det.jpg"%(save_folder))
plt.close()
n_1det=n
bins1_keV=bins






n1det_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]  
for i in range(26):
    n,bins,patches=plt.hist(Energy_1det_ch[i], bins=np.arange(0,2500,20),label="%d events"%len(Energy_1det_ch[i]), histtype="step")
    plt.legend(fontsize=15)
    label("Energy (keV)", "Counts", 8, 5)
    n1det_energy[i]=n
    plt.title("1 detector: Energy Spectrum, Channel %d "%i, fontsize=16)
    plt.savefig("%s/Energy-1Det_Ch%d.jpg"%(save_folder,i))
    plt.close()

bins1det_ch=bins







    
    
  





print("2 Detector study starts")    
Energy_2det=[]
Channel_2det=[]
Energy_2det_inner=[]
count=0
for l in Timed_Energy_event:
    if len(l)==2:
        E1=data[2][l[0]]
        E2=data[2][l[1]]
        E=np.array([E1, E2])
        Energy_2det.extend([E])
        p=[int(data[0][l[0]]), int(data[0][l[1]])]
        Channel_2det.extend([p])
        E_inn=0
        for m in range(2):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_2det_inner.extend([E_inn])
        
Energy_2det=np.array(Energy_2det, dtype=object)
Energy_2det_add=Energy_2det[:,0]+Energy_2det[:,1]



N=300
h=plt.hist2d(Energy_2det[:,0], Energy_2det[:,1], bins=(np.arange(0,2000,2000/N), np.arange(0,2000,2000/N)), norm = LogNorm())
# h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("2 detector coincidence events distribution", fontsize=14)
label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 8,5)
plt.savefig("%s/2Det-DistribuionPlot.jpg"%(save_folder))
plt.close()




"""
print("3 Detector study starts") 
Energy_3det=[]
Channel_3det=[]
Energy_3det_inner=[]
for l in Timed_Energy_event:
    if len(l)==3:
        E1=data[2][l[0]]
        E2=data[2][l[1]]
        E3=data[2][l[2]]
        E=np.array([E1, E2, E3])
        Energy_3det.extend([E])
        p=[int(data[0][l[0]]), int(data[0][l[1]]),  int(data[0][l[2]])]
        Channel_3det.extend([p])
        E_inn=0
        for m in range(3):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_3det_inner.extend([E_inn])
Energy_3det=np.array(Energy_3det, dtype=object)
Energy_3det_add=Energy_3det[:,0]+Energy_3det[:,1]+Energy_3det[:,2]



print("4 Detector study starts") 
Energy_4det=[]
Channel_4det=[]
Energy_4det_inner=[]
for l in Timed_Energy_event:
    if len(l)==4:
        E1=data[2][l[0]]
        E2=data[2][l[1]]
        E3=data[2][l[2]]
        E4=data[2][l[3]]
        E=np.array([E1, E2, E3, E4])
        Energy_4det.extend([E])
        p=[int(data[0][l[0]]), int(data[0][l[1]]),  int(data[0][l[2]]) ,  int(data[0][l[3]])]
        Channel_4det.extend([p])
        E_inn=0
        for m in range(4):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_4det_inner.extend([E_inn])
Energy_4det=np.array(Energy_4det, dtype=object)
Energy_4det_add=Energy_4det[:,0]+Energy_4det[:,1]+Energy_4det[:,2]+Energy_4det[:,3]








 
        
        


        
        
print("5 Detector study starts") 
Energy_5det=[]
Channel_5det=[]
Energy_5det_inner=[]
for l in Timed_Energy_event:
    if len(l)==5:
        E1=data[2][l[0]]
        E2=data[2][l[1]]
        E3=data[2][l[2]]
        E4=data[2][l[3]]
        E5=data[2][l[4]]
        E=np.array([E1, E2, E3, E4, E5])
        Energy_5det.extend([E])
        p=[int(data[0][l[0]]), int(data[0][l[1]]),  int(data[0][l[2]]) ,  int(data[0][l[3]]),  int(data[0][l[4]])]
        Channel_5det.extend([p])
        E_inn=0
        for m in range(5):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_5det_inner.extend([E_inn])
Energy_5det=np.array(Energy_5det, dtype=object)
Energy_5det_add=Energy_5det[:,0]+Energy_5det[:,1]+Energy_5det[:,2]+Energy_5det[:,3]+Energy_5det[:,4]
"""






f = open('%s/Output_Run-%d_%d.pickle'%(save_folder,folder_index, index), 'wb')
# f = open('%s/Output_Run.pickle'%(save_folder), 'wb')
Parameter={'Date': Date, 'Time': Data_Time, 'Count': nEvent, 'Mean': mean, 'Sigma': sigma,
           'bin_last':bin_last,'bin_width':bin_w, "n_energy":n_energy, "count_mult":count_mult, 
           'n_EnergyInner':n_EnergyInner, 'n_EnergyAll':n_EnergyAll, 'bins_combined':bins_combined,
            "n_1det":n_1det, "n1det_Ch":n1det_energy, "bins1_keV":bins1_keV, "count1det_ch":count1_ch}
pickle.dump(Parameter, f)
f.close()





# f = open('%s/Energy_Output_Run-%d_%d.pickle'%(save_folder,folder_index, index), 'wb')
# # f = open('%s/Output_Run.pickle'%(save_folder), 'wb')
# Parameter={'Date': Date, 'Time': Data_Time, 'Count': nEvent, 'Mean': mean, 'Sigma': sigma,
#            "Energy_1det":Energy_1det, "Energy_2det":Energy_2det, "Energy_3det":Energy_3det,
#            "Energy_4det":Energy_4det, "Energy_5det":Energy_5det, 
#           "Channel_1det":Channel_1det, "Channel_2det":Channel_2det, "Channel_3det":Channel_3det,
#            "Channel_4det":Channel_4det, "Channel_5det":Channel_5det}
# pickle.dump(Parameter, f)
# f.close()