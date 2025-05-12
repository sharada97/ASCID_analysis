import os, sys, math
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from ROOT import gROOT, TFile, TTree, TGraph, TChain
import warnings
import glob
import ROOT
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
from collections import OrderedDict 



def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],int(y[i]),fontsize=15)

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



energy_threshold=100

# run_index=int(sys.argv[2])       #Specify file name inside folder



folder=int(sys.argv[1])    #Specify folder name
read_folder="/data7/CsI_Project/Output_new"
save_folder="Figure/Output_new"


# read_folder="../Simulation"
# save_folder="Figure/try"
run=np.arange(1)  #selcting n number of files
innerChannel=np.arange(16,26)




t=TChain("Deposited_energy1")
for p in run:
    t.Add("%s/output0_t%d.root"%(read_folder,p))   
N=t.GetEntries() 
print("Total entries is %d"%N)
nEvent=N


mean=np.zeros((26, 2))
sigma=np.zeros((26, 2))
n_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

Energy_Ch=np.zeros((N,26))
for i in range(1,27):
    t=TChain("Deposited_energy%d"%i)
    for p in run:
        t.Add("%s/output0_t%d.root"%(read_folder,p))   
    N=t.GetEntries() 
    print("Total entries for Channel %d is %d"%(i-1,N))
    Energy=[]
    for k in range(N):
        l=np.zeros(1)
        t.SetBranchAddress("Deposited_energy%d_event"%i,l)
        t.GetEntry(k)
        Energy.extend(l)
#         if k%5000000==0:
#             print("reading %d M"%int(k/1000000))
    Energy=np.array(Energy)
#     Energy=Energy[Energy!=0]
    Energy_Ch[:,i-1]=Energy
    n,bins,patches=plt.hist(Energy,bins=np.arange(15,2000,15), label="%d events"%len(Energy), histtype="step")
    n_energy[i-1]=n
    a,mean[i-1,0], sigma[i-1,0], fig=gaus_fit(365,630,0, max(n))
    a,mean[i-1,1], sigma[i-1,1], fig=gaus_fit(1100,1430,1200, max(n)/3)
    print(mean[i-1,0], mean[i-1,1])
    label("Energy (keV)", "Counts", 8, 5)
    plt.title("Channel %d Energy Spectrum"%(i-1), fontsize=17)
    plt.legend(fontsize=15)
    plt.savefig("%s/Energy_ch-%d.jpg"%(save_folder,i-1))
    plt.close()
    Energy=[]
    

  
Timed_SmEnergy_event=[]   #Energy group from event
Timed_SmEnergy_Channel=[]   #Channel group from event
SmEnergy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]  # individual detector hit
EventEnergy_inner=[]
EventEnergy_all=[]
count=np.zeros(26)    #26 detector multiplicity distribution
bins=np.arange(15,2400,15)
for n,j in enumerate(Energy_Ch):
    group_energy=[]
    group_channel=[]
    energyall=0
    energyinner=0
    for k in range(0,26):
        if j[k]>energy_threshold:            #100keV energy threshold
            mu=j[k]             
            smearE=mu
            SmEnergy[k].append(smearE)
            group_energy.append(smearE)
            energyall+=smearE
            if 16<k<26:
                energyinner+=smearE
            group_channel.append(k)
    count[int(len(group_energy))]+=1
    
    EventEnergy_all.append(energyall)
    if energyinner!=0:
        EventEnergy_inner.append(energyinner)
    
    Timed_SmEnergy_event.append(group_energy)
    Timed_SmEnergy_Channel.append(group_channel)
    if n%1000000==0:
        print("done %d M"%int(n/1000000))        
print("Grouping of the event done")
    
    

plt.rcParams['figure.figsize'] = [12,5]
n,bins,patches=plt.hist(EventEnergy_all, bins=np.arange(0,4000,20),label="%d events"%len(EventEnergy_all), histtype="step")
plt.legend(fontsize=15)
plt.xticks(np.arange(0, 4000, 250))
label("Energy (keV)", "Counts", 12, 5)
plt.title("All detector coincidence Spectrum  ", fontsize=14)
plt.savefig("%s/Energy-AllCoin.jpg"%(save_folder))
plt.close()
n_EnergyAll=n        
        
 

n,bins,patches=plt.hist(EventEnergy_inner, bins=np.arange(0,4000,20),label="%d events"%len(EventEnergy_inner), histtype="step")
plt.legend(fontsize=15)
plt.xticks(np.arange(0, 4000, 250))
label("Energy (keV)", "Counts", 12, 5)
plt.title("Inner detector coincidence Spectrum  ", fontsize=14)
plt.savefig("%s/Energy-InnCoin.jpg"%(save_folder))
plt.close()
n_EnergyInner=n 
bins_combined=bins


# Energy_Ch=[]
# n_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# for k in range(26):
#     n,bins,patches=plt.hist(SmEnergy[k],bins=np.arange(15,2400,20), label="%d events"%len(SmEnergy[k]), histtype="step")
#     n_energy[k]=n
#     label("Energy (keV)", "Counts", 8, 5)
#     plt.title("Channel %d Energy Spectrum"%k, fontsize=17)
#     plt.legend(fontsize=15)
#     plt.savefig("%s/Energy_ch-%d_Smear.jpg"%(save_folder,k))
#     plt.close()





count_Mult=count
count=np.array(count, dtype=int)
data_count = {'1 Det':count[1], '2 Det':count[2], '3 Det':count[3],
        '4 Det':count[4], '5 Det':count[5], '6 Det':count[6],'7 Det':count[7],
             '8 Det':count[8], '9 Det':count[9], '10 Det':count[10]}
courses = list(data_count.keys())
values = list(data_count.values())
  
fig = plt.figure(figsize = (10, 5))
 
plt.bar(courses, values,
        width = 0.4, label="Total decay events= %d"%N)   
addlabels(courses, count[1:])
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("Multiplicity Distribution - Simulation (Na22)- %d keV threshold"%energy_threshold, fontsize="15")
plt.savefig("%s/Multiplicity_Sim.jpg"%(save_folder))
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
fig = plt.figure(figsize = (10, 5))
 
plt.bar(courses, values,
        width = 0.4, label="Total events= %d "%(nEvent))   
addlabels(courses, count[1:])
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("Multiplicity Distribution in Percentage; Threshold-%dkeV"%(energy_threshold), fontsize="15")
plt.savefig("%s/Multiplicity_Percent-%dkeV.jpg"%(save_folder, energy_threshold))    
plt.close()



#1 -Det coincidence starts

print("1 -Det coincidence starts")
Energy_1det=[]
Energy_1det_ch=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
Energy_1det_inner=[]
count1_ch=np.zeros(26)
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==1:
        E=l[0]
        p=Timed_SmEnergy_Channel[n][0]
        Energy_1det.extend([E])
        Energy_1det_ch[p].extend([E])
        count1_ch[p]+=1
        if p in innerChannel:
            Energy_1det_inner.extend([E])


        
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
        width = 0.4, label="Total events= %d "%(count1_ch.sum()))   #Change time here
plt.legend(fontsize=14)
label("No. of coincident detectors", "No. of events", 8, 5)
plt.title("1 Det event event topology ( Na22)", fontsize="15")
plt.savefig("%s/CountEvents-1Coin_5-5-%dkeV.jpg"%(save_folder, energy_threshold))
plt.close()


n,bins,patches=plt.hist(Energy_1det, bins=np.arange(0,2500,20),label="%d events"%len(Energy_1det), histtype="step")
plt.legend(fontsize=15)
label("Energy (keV)", "Counts", 8, 5)
plt.title("1 detector: Combined Spectrum  ", fontsize=17)
plt.savefig("%s/Energy-1Det.jpg"%(save_folder))
plt.close()
n_1det=n
bins1_keV=bins


n1det_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]  
for k in range(26):
    n,bins,patches=plt.hist(Energy_1det_ch[k], bins=np.arange(0,2100,15),label="%d events"%len(Energy_1det_ch[k]), histtype="step")
    plt.legend(fontsize=15)
    n1det_energy[k]=n
    label("Energy (keV)", "Counts", 8, 5)
    plt.title("1 detector Ch %d: Spectrum  "%k, fontsize=14)
    plt.savefig("%s/Energy-1Det-ch%d.jpg"%(save_folder,k))
    plt.close()
bins1det_ch=bins    
    
    


    
#2 -Det coincidence starts

print("2 -Det coincidence starts")
Energy_2det=[]
Energy_2det_inner=[]
c=0
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==2:
        E1=l[0]
        E2=l[1]
        p=[Timed_SmEnergy_Channel[n][0], Timed_SmEnergy_Channel[n][1]]
        E=np.array([E1, E2])
        Energy_2det.extend([E])
        E_inn=0
        for m in range(2):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_2det_inner.extend([E_inn])
        c+=1
Energy_2det=np.array(Energy_2det)
Energy_2det_add=Energy_2det[:,0]+Energy_2det[:,1]

N=200
h=plt.hist2d(Energy_2det[:,0], Energy_2det[:,1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
plt.title("2 detector coincidence events distribution, %d events"%len(Energy_2det[:,1]), fontsize=17)
label("First detector Energy (keV)", "Second detector Energy (keV)", 8,5)
plt.savefig("%s/2Det-DistribuionPlot.jpg"%(save_folder))
plt.close()




#3 -Det coincidence starts

print("3 -Det coincidence starts")
Energy_3det=[]
Energy_3det_inner=[]
time_diff=[]
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==3:
        E1=l[0]
        E2=l[1]
        E3=l[2]
        E=np.array([E1, E2, E3])
        Energy_3det.extend([E])
        p=[Timed_SmEnergy_Channel[n][0], Timed_SmEnergy_Channel[n][1], Timed_SmEnergy_Channel[n][2]]
        E_inn=0
        for m in range(3):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_3det_inner.extend([E_inn])
            
Energy_3det=np.array(Energy_3det, dtype=object)
Energy_3det_add=Energy_3det[:,0]+Energy_3det[:,1]+Energy_3det[:,2]


# f = open('DataInfo_Run.pickle', 'rb')
# dat = pickle.load(f)
# mean=dat['mean']
# sigma=dat['sigma']
# f.close()
# for k in range(26):
#     mean[k,0]=mean[k,0]/mean[k,1]*1275
#     sigma[k,0]=sigma[k,0]/mean[k,1]*1275
#     sigma[k,1]=sigma[k,1]/mean[k,1]*1275
#     mean[k,1]=mean[k,1]/mean[k,1]*1275


    
    
Energy_3det_sl=[]
check=[]
check2=[]
co=0
sett=0
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==3:
        co+=1
        p=np.array([int(Timed_SmEnergy_Channel[n][0]), int(Timed_SmEnergy_Channel[n][1]), int(Timed_SmEnergy_Channel[n][2])])
        E=np.array([l[0], l[1], l[2]])
        k=0
        sett+=1
        for n1 in range(2):
            for n2 in range(n1+1,3):
                n3=list(set(np.arange(3))-set([n1,n2]))[0]
                mu1=mean[:,0][p[n1]]  #511   of first hit
                mu2=mean[:,0][p[n2]]  #511   of second hit
                mu3=mean[:,1][p[n3]] #1275   of third hit
                sig1=sigma[:,0][p[n1]]
                sig2=sigma[:,0][p[n2]]                
                sig3=sigma[:,1][p[n3]]
                
                if mu3-2*sig3<E[n3]<mu3+2*sig3:
                    k+=1
                    check.extend([sett])
                    if k==1:
                        Energy_3det_sl.extend([[E[n1], E[n2]]])
#                     if mu1-2*sig1<E[n1]<mu1+2*sig1 and mu2-2*sig2<E[n2]<mu2+2*sig2:
#                         check2.extend([sett])

Energy_3det_sl=np.array(Energy_3det_sl, dtype=object)
N=150
# print(len(Energy_3det_sl[:,0]))
h=plt.hist2d(Energy_3det_sl[:,0], Energy_3det_sl[:,1], bins=(np.arange(0,2100,2100/N), np.arange(0,2100,2100/N)), norm = LogNorm())
plt.title("3 detector: Hit Any-1275 keV, %d events"%(len(Energy_3det_sl[:,0])), fontsize=14)
label("1 hit detector Energy (keV)", "2 hit detector Energy (keV)", 8,5)
plt.savefig("%s/3Det-DistribuionPlot-Any-1275.jpg"%(save_folder))
plt.close()





#4 -Det coincidence starts

print("4 -Det coincidence starts")
Energy_4det_inner=[]
Energy_4det=[]
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==4:
        E1=l[0]
        E2=l[1]
        E3=l[2]
        E4=l[3]
        E=np.array([E1, E2, E3, E4])
        Energy_4det.extend([E])
        p=[Timed_SmEnergy_Channel[n][0], Timed_SmEnergy_Channel[n][1], Timed_SmEnergy_Channel[n][2], Timed_SmEnergy_Channel[n][3]]
        E_inn=0
        for m in range(4):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_4det_inner.extend([E_inn])
            
            
Energy_4det=np.array(Energy_4det, dtype=object)
Energy_4det_add=Energy_4det[:,0]+Energy_4det[:,1]+Energy_4det[:,2]+Energy_4det[:,3]



#shared 511 keV
print("Staring For 511 keV")
Energy_4det_sl=[]
check=[]
check2=[]
co=0
sett=0
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==4:
        co+=1
        E=np.array([l[0], l[1], l[2], l[3]])
        p=np.array([int(Timed_SmEnergy_Channel[n][0]), int(Timed_SmEnergy_Channel[n][1]), int(Timed_SmEnergy_Channel[n][2]), int(Timed_SmEnergy_Channel[n][3])])
        sett+=1
        k=0
        for n1 in range(3):
            for n2 in range(n1+1,4):
                n3=list(set(np.arange(4))-set([n1,n2]))
                mu1=[mean[p[n1],0],mean[p[n1],1]]
                mu2=[mean[p[n2],0],mean[p[n2],1]]
                sig1=[sigma[p[n1],0],sigma[p[n1],1]]
                sig2=[sigma[p[n2],0],sigma[p[n2],1]]
                
                if (mu1[0]-2*sig1[0]<E[n1]<mu1[0]+2*sig1[0] and mu2[1]-2*sig2[1]<E[n2]<mu2[1]+2*sig2[1])\
                or (mu1[1]-2*sig1[1]<E[n1]<mu1[1]+2*sig1[1] and mu2[0]-2*sig2[0]<E[n2]<mu2[0]+2*sig2[0]):  
                    k+=1
                    if k==1:
                        Energy_4det_sl.extend([[E[n3][0], E[n3][1]]])    #2-511, 1275, then spectrum of rest
#                     check.extend([sett])
#                     if 524-130<E[n3][0]+ E[n3][1]<524+130:
#                         check2.extend([sett])
#                         Energy_4det_2sl[k].extend([[E[n3][0], E[n3][1]]])  #above in specific region
        if co%1000000==0:
            print("well done %d M"%int(co/1000000))
Energy_4det_sl=np.array(Energy_4det_sl)
Energy_4det_sl511=Energy_4det_sl 
N=150
h=plt.hist2d(Energy_4det_sl[:,0], Energy_4det_sl[:,1], bins=(np.arange(0,2100,2100/N), np.arange(0,2100,2100/N)), norm = LogNorm())
plt.title("4 detector: Hit Any-511_1275 keV, %d events"%(len(Energy_4det_sl[:,0])), fontsize=14)
label("1 hit detector Energy (keV)", "2 hit detector Energy (keV)", 8,5)
plt.savefig("%s/4Det-DistribuionPlot-Any-511_1275.jpg"%(save_folder))
plt.close()





#shared 1275 keV
print("Staring For 1275 keV")
Energy_4det_sl=[]

check=[]
check2=[]
co=0
sett=0
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==4:
        co+=1
        E=np.array([l[0], l[1], l[2], l[3]])
        p=np.array([int(Timed_SmEnergy_Channel[n][0]), int(Timed_SmEnergy_Channel[n][1]), int(Timed_SmEnergy_Channel[n][2]), int(Timed_SmEnergy_Channel[n][3])])
        sett+=1
        k=0
        for n1 in range(3):
            for n2 in range(n1+1,4):
                n3=list(set(np.arange(4))-set([n1,n2]))
                mu1=[mean[p[n1],0],mean[p[n1],1]]
                mu2=[mean[p[n2],0],mean[p[n2],1]]
                sig1=[sigma[p[n1],0],sigma[p[n1],1]]
                sig2=[sigma[p[n2],0],sigma[p[n2],1]]
                
                if (mu1[0]-2*sig1[0]<E[n1]<mu1[0]+2*sig1[0] and mu2[0]-2*sig2[0]<E[n2]<mu2[0]+2*sig2[0]):
                    k+=1
                    if k==1:
                        Energy_4det_sl.extend([[E[n3][0], E[n3][1]]])    #2-511, 1275, then spectrum of rest
#                     check.extend([sett])
#                     if 1275-170<E[n3][0]+ E[n3][1]<1275+170:
#                         check2.extend([sett])
        if co%1000000==0:
            print("well done %d M"%int(co/1000000))
Energy_4det_sl=np.array(Energy_4det_sl)
Energy_4det_sl1275=Energy_4det_sl 



N=230
h=plt.hist2d(Energy_4det_sl[:,0], Energy_4det_sl[:,1], bins=(np.arange(0,2100,2100/N), np.arange(0,2100,2100/N)), norm = LogNorm())
plt.title("4 detector: Hit Any-511_511 keV, %d events"%(len(Energy_4det_sl[:,0])), fontsize=14)
label("1 hit detector Energy (keV)", "2 hit detector Energy (keV)", 8,5)
plt.savefig("%s/4Det-DistribuionPlot-Any-511_511.jpg"%(save_folder))
plt.close()










#5 -Det coincidence starts

print("5 -Det coincidence starts")
Energy_5det_inner=[]
Energy_5det=[]
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==5:
        E1=l[0]
        E2=l[1]
        E3=l[2]
        E4=l[3]
        E5=l[4]
        E=np.array([E1, E2, E3, E4, E5])
        Energy_5det.extend([E])
        p=[Timed_SmEnergy_Channel[n][0], Timed_SmEnergy_Channel[n][1], Timed_SmEnergy_Channel[n][2], Timed_SmEnergy_Channel[n][3],                   Timed_SmEnergy_Channel[n][4]]
        E_inn=0
        for m in range(5):
            if p[m] in innerChannel:
                E_inn+=E[m]
        if E_inn!=0:
            Energy_5det_inner.extend([E_inn])
Energy_5det=np.array(Energy_5det, dtype=object)
Energy_5det_add=Energy_5det[:,0]+Energy_5det[:,1]+Energy_5det[:,2]+Energy_5det[:,3]+Energy_5det[:,4]


Energy_5det_sl=[]

check=[]
check2=[]
co=0
sett=0
for n,l in enumerate(Timed_SmEnergy_event):
    if len(l)==5:
        co+=1
        E=np.array([l[0], l[1], l[2], l[3], l[4]])
        p=np.array([int(Timed_SmEnergy_Channel[n][0]), int(Timed_SmEnergy_Channel[n][1]), int(Timed_SmEnergy_Channel[n][2]), int(Timed_SmEnergy_Channel[n][3]),int(Timed_SmEnergy_Channel[n][4])])
        sett+=1
        k=0
        for n1 in range(4):
            for n2 in range(n1+1,5):
                n3=list(set(np.arange(5))-set([n1,n2]))
                mu1=[mean[p[n1],0],mean[p[n1],1]]
                mu2=[mean[p[n2],0],mean[p[n2],1]]
                sig1=[sigma[p[n1],0],sigma[p[n1],1]]
                sig2=[sigma[p[n2],0],sigma[p[n2],1]]
                
                
                if 1275-200<E[n3][0]+ E[n3][1]+ E[n3][2]<1275+200:
                    k+=1
                    check.extend([sett])
                    if k==1:
                        Energy_5det_sl.extend([np.array([E[n1], E[n2]])])
                    
                    #above in specific region
#                     if mu1[0]-2*sig1[0]<E[n1]<mu1[0]+2*sig1[0] and mu2[0]-2*sig2[0]<E[n2]<mu2[0]+2*sig2[0]:
#                         check2.extend([sett])
        if co%1000000==0:
            print("well done %d M"%int(co/1000000))

            
Energy_5det_sl=np.array(Energy_5det_sl)
N=250
h=plt.hist2d(Energy_5det_sl[:,0], Energy_5det_sl[:,1], bins=(np.arange(0,2100,2100/N), np.arange(0,2100,2100/N)), norm = LogNorm())
plt.title("5 detector: Hit Comb-1275 keV, %d events"%(len(Energy_5det_sl[:,0])), fontsize=14)
label("1 hit detector Energy (keV)", "2 hit detector Energy (keV)", 8,5)
plt.savefig("%s/5Det-DistribuionPlot-AnyComb-1275.jpg"%(save_folder))
plt.close()    
            

    
f = open('%s/OutputCoadded_Run.pickle'%(save_folder), 'wb')
Parameter={"n_energy":n_energy, "n1det_energy":n1det_energy, "bins1det_ch":bins1det_ch, 
            "Energy_1det_all":Energy_1det,
           "Energy_2det_all":Energy_2det_add,
           "Energy_3det_all":Energy_3det_add,
           "Energy_4det_all":Energy_4det_add,
           "Energy_5det_all":Energy_5det_add, 
          "Energy_1det_inner":Energy_1det_inner, "Energy_2det_inner":Energy_2det_inner, "Energy_3det_inner":Energy_3det_inner,
          "Energy_4det_inner":Energy_4det_inner, "Energy_5det_inner":Energy_5det_inner}
pickle.dump(Parameter, f)
f.close()


  
f = open('%s/Output_Run.pickle'%(save_folder), 'wb')
Parameter={'nEvent': nEvent, "count_mult":count_Mult, 
           'n_EnergyInner':n_EnergyInner, 'n_EnergyAll':n_EnergyAll, 'bins_combined':bins_combined,
           "n_1det":n_1det, "bins1_keV":bins1_keV, "count1det_ch":count1_ch, 
           "Energy_2det":Energy_2det, "Energy_3det_sl":Energy_3det_sl,
          "Energy_4det_sl511":Energy_4det_sl511, "Energy_4det_sl1275":Energy_4det_sl1275, "Energy_5det_sl":Energy_5det_sl}
pickle.dump(Parameter, f)
f.close()


print("***************Folder %d done*************"%folder)
    
    

