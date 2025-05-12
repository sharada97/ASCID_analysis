import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from natsort import natsorted
import sys
import os
import pandas as pd
import glob
import time
import pickle
import shutil
import matplotlib.pyplot as plt
from ASCID_variable import *
from matplotlib.colors import LogNorm

def label(self, x,y,a,b,pur="all"):
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

def read_data(files):
    """
    reads the .CSV files and provide a data array with channel, time and energy info.
    Parameters:
        files: list of .CSV files

    return:
        data: list, sorted with increase in detection time. data[0]: channel, data[1]: time (s), data[2]: energy
    """

    files=natsorted(files)
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

    return data, Data_Time


def data_calibration(data, energy_source):
    data=np.array(data).T
    data=data[data[:, 1].argsort()]
    data=np.array(data).T

    f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'rb')
    dat = pickle.load(f)
    mean=dat['Mean']
    f.close()
    
    nEvents = len(data[0])
    for j in range(nEvents):
        if data[0][j]<26:
            data[2][j]=data[2][j]/mean[int(data[0][j]),0]*energy_source   #change to 0 for Sc-46
    print("Calibration done")
    return data



class ASCID_analysis:
    # from ASCID_variable import *
    # Initialize the car with make, model, and year
    def __init__(self, data, Data_Time):
        self.data = data
        self.nEvent = len(data[0])
        data = None
        self.data_time = Data_Time
        self.energy_bin = energy_bin
        self.energy_end = energy_end
        self.energyThreshold_CsI = energyThreshold_CsI
        self.source = source
        self.timeWindow = timeWindow

        
        self.Bins = np.arange(self.energyThreshold_CsI, energy_end, energy_bin)
        self.Timed_Energy_event = None
        self.Multiplicity_raw = None
        self.Multiplicity = None
        self.EventEnergyAll = None
        self.EventEnergyInner = None
        self.EventEnergyAll_NVeto = None
        self.EventEnergyInner_NVeto = None
        self.antiveto_index = None
        self.count1_ch = None
        self.n_1det = None 
        self.n_EnergyAll = None
        self.n_EnergyAll_NVeto = None
        self.n_EnergyInner = None
        self.n_EnergyInner_NVeto = None
        self.energy_exclusive = None
        self.multipliciy_any = None
        
        
        """
        Find best range for gaussian fit by calculating minimum error over different ranges
        Parameters:
            x: list x-axis, energy bin
            y: list count values taken from the histogram.
            peak_position_1st: initial guess of the 1st peak
            window: width of the range from the peak position
            shift_width: shifting of the range to find the better range for the fit
            max_shift: maximum iteration you want to do to find the best fit
        
        return:
            plot: if plot = True is given
            params: [mean, sigma], [mean, sigma] of 2 peaks that fitted. 
        """


    def label(self, x,y,a,b,pur="all"):
        plt.xlabel(x,fontsize=18)
        plt.ylabel(y,fontsize=18)
        if pur=="all":
            plt.rcParams['figure.figsize'] = [a,b]
        elif pur=="chi2":
            plt.title("Chi2 Analysis, S-%s, Channel %s"%(a,b), fontsize=20)
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)
        return 

    def addlabels(self,x,y):
        for i in range(len(x)):
            plt.text(i,y[i],int(y[i]),fontsize=12)


    def data_preprocess(self):
        """
        cut out the data which are bwlow the given threshold value
        Parameters:
            data: array, data obtained from the .CSV files
        return:
            data: filtered it by putting energy threshold
        """
        
        # data_process = np.array(self.data).T
        # Veto_threshold = np.ones(N_det_panel)*energyThreshold[1]
        # Threshold_CsI=np.ones(N_det_CsI)*energyThreshold[0]
        # Threshold_null = np.zeros(15)
        # Threshold = np.append(np.append(Threshold_CsI,Threshold_null), Veto_threshold)
        # # print(Threshold)
        # # np.shape(data_process)
        # data_process1 = []
        # for j in range(self.nEvent):
        #     arr = [self.data[0][j], self.data[1][j], self.data[2][j]]
        #     if data_process[j, 2] > Threshold[int(data_process[j, 0])]:
        #         data_process1.extend([arr])
        # self.data = np.array(data_process1).T

        data_process = np.array(self.data).T

        # Threshold arrays
        Veto_threshold = np.ones(N_det_panel) * energyThreshold[1]
        Threshold_CsI = np.ones(N_det_CsI) * self.energyThreshold_CsI
        Threshold_null = np.zeros(15)
        Threshold = np.concatenate((Threshold_CsI, Threshold_null, Veto_threshold))
        
        # Boolean mask to filter rows that meet the threshold condition
        mask = data_process[:, 2] > Threshold[data_process[:, 0].astype(int)]
        
        # Filtered data without intermediate list creation
        self.data = data_process[mask, :3].T
        print("Data preprocessing done")
        # return data_process[mask, :3].T

    def coincidence(self):
        """
        Provides grouping of the events based on the coincidence window time
        parameter:
            taken from the class
        return: Timed_Energy_event
        """
        # from ASCID_variable import csiCh, innerCh, VetoCh, timeWindow
        nEvents = len(self.data[0])
        Sig_index=-1
        eventTime=0
        Timed_Energy_event=[]
        for j in range(nEvents):
            if  self.data[1][j]>eventTime+2*self.timeWindow and self.data[0][j] in csiCh:   #Avoiding center detector overlap
                eventTime=self.data[1][j]
                eventChannel=int(self.data[0][j])
                Ev_indx=j
                Sig_index+=1
                EventEnergy=self.data[2][j]
                a=[]
                a.extend([Ev_indx])
                
                while Ev_indx>0 and abs(self.data[1][Ev_indx-1]-eventTime)<=self.timeWindow:
                    Ev_indx-=1
                    a.extend([Ev_indx])
                Ev_indx=j
                
                while Ev_indx<nEvents-1 and abs(self.data[1][Ev_indx+1]-eventTime)<=self.timeWindow:
                    Ev_indx+=1    
                    a.extend([Ev_indx])
                a = np.unique(np.sort(a)) # a contains the index within the defined coincidence window
                Timed_Energy_event.extend([a])
        self.Timed_Energy_event = Timed_Energy_event
        Timed_Energy_event=[]
        print("data coincidence done")
        # return Timed_Energy_event


    def CoaddedEnergy(self):
        # from ASCID_variable import csiCh, innerCh, VetoCh
        """
        Provides co-added all 26 and inner 10 spectrum with and without the external veto coincidence. 
        parameter:
            taken from the class
        return:  
        """
        EventEnergyAll = []
        EventEnergyInner = []
        EventEnergyAll_NVeto = []
        antiveto_index = []
        EventEnergyInner_NVeto = []
        Multiplicity_count = []
        Multiplicity = np.zeros(len(csiCh))
        Multiplicity_raw = np.zeros(len(csiCh))
        Multiplicity_889 = np.zeros(len(csiCh))
        Multiplicity_1132 = np.zeros(len(csiCh))
        
        
        csiCh_set = set(csiCh)
        innerCh_set = set(innerCh)
        VetoCh_set = set(VetoCh)

        sl=[[],[]]
        min_val = [889-1.2*77, 1132-1.2*94, 889-1.2*77]
        max_val = [889+1.2*77, 1132+1.2*94, 889+1.2*77]
        min_val= np.array(min_val, dtype=int)
        max_val = np.array(max_val, dtype=int)

        energy_exclusive = [[], []]
        multipliciy_any = np.zeros((len(csiCh), 2))
        
        for k, j in enumerate(self.Timed_Energy_event):
            event_Channel = self.data[0][j]
            event_Energies = self.data[2][j]
            
            vetoMatch = VetoCh_set.intersection(event_Channel)
            non_veto_indices = np.isin(event_Channel, list(VetoCh), invert=True)
            multiplicity = np.sum(non_veto_indices)
            Multiplicity_raw[multiplicity - 1] += 1


            # In the loop
            event_in_csiCh = np.array([ch in csiCh_set for ch in event_Channel])
            event_in_innerCh = np.array([ch in innerCh_set for ch in event_Channel])

            
            EventEnergyall = np.sum(event_Energies[event_in_csiCh])
            EventEnergyinn = np.sum(event_Energies[event_in_innerCh])
            
            EventEnergyall_NVeto = 0
            EventEnergyinn_NVeto = 0

            
            
            if not vetoMatch:
                antiveto_index.append(k)
                multiplicity = len(j)
                Multiplicity[multiplicity - 1] += 1
                EventEnergyall_NVeto = np.sum(event_Energies[event_in_csiCh])
                EventEnergyinn_NVeto = np.sum(event_Energies[event_in_innerCh])
                # if 877-1.2*72<EventEnergyall_NVeto <877+1.2*72:
                #     Multiplicity_889[multiplicity - 1] += 1
                # elif 1109-1.2*94<EventEnergyall_NVeto <1109+1.2*94:
                #     Multiplicity_1132[multiplicity - 1] += 1

                # for element in event_Energies:
                #     if min_val[0] <= element <= max_val[0]:
                #         energy_rest = sum(event_Energies) - element
                #         energy_exclusive[0].extend([energy_rest])
                #         multipliciy_any[multiplicity, 0]+=1
                #     elif min_val[1] <= element <= max_val[1]:
                #         energy_rest = sum(event_Energies) - element
                #         energy_exclusive[1].extend([energy_rest])
                #         multipliciy_any[multiplicity, 1]+=1
            
            EventEnergyAll.append(EventEnergyall)
            EventEnergyInner.append(EventEnergyinn)
            EventEnergyAll_NVeto.append(EventEnergyall_NVeto)
            EventEnergyInner_NVeto.append(EventEnergyinn_NVeto)
            
            if k % 1000000 == 0:
                print(f"{k // 1000000} M events done")
                time.sleep(5)
         
        
        # print(EventEnergyall, EventEnergyInner)
        EventEnergyAll = np.array(EventEnergyAll)
        EventEnergyInner = np.array(EventEnergyInner)
        EventEnergyAll_NVeto = np.array(EventEnergyAll_NVeto)
        EventEnergyInner_NVeto = np.array(EventEnergyInner_NVeto)
    
        EventEnergyAll = EventEnergyAll[EventEnergyAll!=0]
        EventEnergyInner = EventEnergyInner[EventEnergyInner!=0]
        EventEnergyAll_NVeto = EventEnergyAll_NVeto[EventEnergyAll_NVeto!=0]
        EventEnergyInner_NVeto = EventEnergyInner_NVeto[EventEnergyInner_NVeto!=0]

        # self.energy_exclusive = energy_exclusive
        # self.multipliciy_any = multipliciy_any

        self.EventEnergyAll = EventEnergyAll
        self.EventEnergyInner = EventEnergyInner
        self.EventEnergyAll_NVeto = EventEnergyAll_NVeto
        self.EventEnergyInner_NVeto = EventEnergyInner_NVeto
        self.Multiplicity = Multiplicity
        self.Multiplicity_raw = Multiplicity_raw
        self.antiveto_index = antiveto_index

        energy_exclusive = None
        multipliciy_any = None
        EventEnergyAll = None
        EventEnergyInner = None
        EventEnergyAll_NVeto = None
        EventEnergyInner_NVeto = None
        Multiplicity_count = None
        Multiplicity = None
        Multiplicity_raw = None

        print("grouped energy done")

        # f = open('%s/Outputmultiplicity_%d_%d.pickle'%(save_folder_each, folder_index, index), 'wb')
        # Parameter={"Multiplicity":Multiplicity, "Multiplicity_raw":Multiplicity_raw,"Multiplicity_889":Multiplicity_889, "Multiplicity_1132":Multiplicity_1132}
        # pickle.dump(Parameter, f)
        # f.close()

    
    
    def added_energy(self, group, plot = False):
        
        """
        Provides plot for the co-added energy spectrums
        parameter:
            group: "All" or "inner". Provides plot for 26 channel combine spectrum if "All" or inner 10 channel combined spectrum if "inner"
            plot: True or False, saves .jpg if given True
        return: n_Energy, n_Energy_NVeto for All and Inner channels. 
        """
        
        #All and Inner detector coincidence energy 
        if group == "All":
            EventEnergy = self.EventEnergyAll
            EventEnergy_NVeto = self.EventEnergyAll_NVeto
        elif group == "inner":
            EventEnergy = self.EventEnergyInner
            EventEnergy_NVeto = self.EventEnergyInner_NVeto
        plt.rcParams['figure.figsize'] = [12,5]
        n1,bins,patches=plt.hist(EventEnergy, bins=self.Bins,label="%d events, all"%len(EventEnergy), histtype="step")
        n2,bins,patches=plt.hist(EventEnergy_NVeto, bins=self.Bins,label="%d events, Veto antocoincidence"%len(EventEnergy_NVeto), histtype="step")
        plt.legend(loc="upper left", fontsize=15)
        # plt.xticks(np.arange(0, 4000, 250))
        self.label("Energy (keV)", "Counts", 12, 5)
        plt.title("%s detector coincidence Spectrum"%(group), fontsize=14)
        plt.grid()
        if plot:
            plt.savefig("%s/Energy-%sCoin.jpg"%(save_folder_each, group))
        plt.close()
        n_Energy=n1
        n_Energy_NVeto=n2

        if group == "All":
            self.n_EnergyAll = n_Energy
            self.n_EnergyAll_NVeto = n_Energy_NVeto
        elif group == "inner":
            self.n_EnergyInner = n_Energy
            self.n_EnergyInner_NVeto = n_Energy_NVeto
        return EventEnergy, EventEnergy_NVeto
        

    def exclusive_search(self):
        min_val = [889-1.2*77, 1132-1.2*94, 889-1.2*77]
        max_val = [889+1.2*77, 1132+1.2*94, 889+1.2*77]
        min_val= np.array(min_val, dtype=int)
        max_val = np.array(max_val, dtype=int)
    
        any_gamma = ["889", "1132", "889"]
        n_exclusive = [[],[]]
        for i in range(2):
            n,bins, patches=plt.hist(self.energy_exclusive[i], bins=np.arange(0,self.energy_end, self.energy_bin), histtype="step", label="%d events"%len(self.energy_exclusive[i]))
            self.label("Energy (ADC)", "Counts", 10,5)
        
            count_889= np.zeros(3)
            for j in self.energy_exclusive[i]:
                if j == 0:
                    count_889[0]+=1
                elif j<min_val[i+1]:
                    count_889[1]+=1
                elif min_val[i+1]<j<max_val[i+1]:
                    count_889[2]+=1
            n_exclusive[i] = n
            plt.legend(fontsize=14)
            textstr = '\n'.join((
                    r'$N_0= %d$  = %.2f %%' % (count_889[0], count_889[0]/len(self.energy_exclusive[i])*100),
                    r'$N_{compton}= %d = %.2f %%' % (count_889[1], count_889[1]/len(self.energy_exclusive[i])*100),
                    r'$N_{%s}= %d = %.2f %%' % (any_gamma[i+1], count_889[2], count_889[2]/len(self.energy_exclusive[i])*100)
                    ))               #"r'$\chi^2/Dof=%.2f$' % (b, )"
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            plt.text(1200, max(n)*0.8, textstr, fontsize=12,
                 verticalalignment='top', bbox=props)
            plt.title("Rest energy to the detector at E_1 = %s keV"%(any_gamma[i]), fontsize = 12)
            plt.savefig("%s/EnergyRest_%s_%dns_%dkeV.png"%(save_folder_each, any_gamma[i], self.timeWindow*1e9, self.energyThreshold_CsI))
            plt.close()
        
        
            count12 = self.multipliciy_any[:, i]
            total = count12.sum()
            count12=count12/count12.sum()*100
            count12=np.array(count12, dtype="float")
            for j in range(len(count12)):
                count12[j]=round(count12[j],2)
            data_count = {'1 Det':count12[1], '2 Det':count12[2], '3 Det':count12[3],
                '4 Det':count12[4], '5 Det':count12[5], '6 Det':count12[6]}
            courses = list(data_count.keys())
            values = list(data_count.values())
            fig = plt.figure(figsize = (10, 5))
            plt.bar(courses, values,
                width = 0.4)   #Change time here
            addlabels(courses, values)
            # plt.legend(fontsize=14)
            self.label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
            plt.title("Multiplicity Distribution-%s keV (%s)"%(any_gamma[i], self.source), fontsize="15")
            plt.savefig("%s/MultiplicityAny_%s.png"%(save_folder_each, any_gamma[i]))
            plt.close()
        
        
        
        f = open('%s/OutputTest_Run5-folder%d_%dns_%dkeV.pickle'%(save_folder_each, folder_index, self.timeWindow*1e9, self.energyThreshold_CsI), 'wb')
        Parameter={'nEvent': self.nEvent, "energy_exclusive":n_exclusive, "Bins":np.arange(0,self.energy_end, self.energy_bin)}
        pickle.dump(Parameter, f)
        f.close()
        print("Exclusive search done")

    
    def multiplicity(self, plot = False, veto = "Veto"):
        """
        Provides plot for the multiplicity plot
        parameter:
            plot: True or False, saves .jpg if given True
            
        return: save plot if asked
        """
        if veto == "Veto":
            count = self.Multiplicity
        elif veto == "noVeto":
            count = self.Multiplicity_raw
        else:
            print("out of options, provide between Veto or noVeto")
            
        total = count.sum()
        count=count/count.sum()*100
        count=np.array(count, dtype="float")
        for j in range(len(count)):
            count[j]=round(count[j],2)
        data_count = {'1 Det':count[0], '2 Det':count[1], '3 Det':count[2],
                '4 Det':count[3], '5 Det':count[4], '6 Det':count[5]}
        courses = list(data_count.keys())
        values = list(data_count.values())
        fig = plt.figure(figsize = (10, 5))
        plt.bar(courses, values,
                width = 0.4, label="Total decay rate= %.1f decay/sec"%(total/(self.data_time)))   #Change time here
        self.addlabels(courses, values)
        plt.legend(fontsize=14)
        self.label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
        plt.title("Multiplicity Distribution (%s), %s"%(self.source, veto), fontsize="15")
        if plot:
            plt.savefig("%s/Multiplicity_%s.jpg"%(save_folder_each, veto))
        plt.close()


    def coincidence_1det(self, plot = False):
        """
        Provides plot for 1-det coincidence energy spectrum and channel wise no. of events distribution across the 26 channels
        parameter:
            plot: True or False, saves .jpg if given True
            
        return: save plot if asked
        """
    
        Energy_1det=[]
        Channel_1det=[]
        Energy_1det_ch=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        Energy_1det_inner=[]
        count1_ch=np.zeros(26)
        for l in self.Timed_Energy_event:
            if len(l)==1:
                E=self.data[2][l[0]]
                p=int(self.data[0][l[0]])
                Energy_1det.extend([E])
                Channel_1det.extend([p])
                Energy_1det_ch[p].extend([E])
                count1_ch[p]+=1
                if p in innerCh:
                    Energy_1det_inner.extend([E])
        
        
        
        
        count1_ch=np.array(count1_ch, dtype=int) 
        
    
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
        self.label("No. of coincident detectors", "No. of events", 8, 5)
        plt.title("1Det coincidence events", fontsize="15")
        if plot:
            plt.savefig("%s/CountEvents-1det.jpg"%(save_folder_each))
        plt.close()
        

        n,bins,patches=plt.hist(Energy_1det, bins=self.Bins, histtype="step",label="%d events"%len(Energy_1det))
        plt.legend(fontsize=15)
        self.label("Energy (keV)", "Counts", 8, 5)
        plt.title("1 detector: Spectrum  ", fontsize=14)
        if plot:
            plt.savefig("%s/Energy-1Det.jpg"%(save_folder_each))
        plt.close()
        n_1det=n
        n1det_energy=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]  
        for i in range(26):
            n,bins,patches=plt.hist(Energy_1det_ch[i], bins=np.arange(0,2800,20),label="%d events"%len(Energy_1det_ch[i]), histtype="step")
            plt.legend(fontsize=15)
            self.label("Energy (keV)", "Counts", 8, 5)
            n1det_energy[i]=n
            plt.title("1 detector: Energy Spectrum, Channel %d "%i, fontsize=16)
            plt.savefig("%s/Energy-1Det_Ch%d.jpg"%(save_folder_each,i))
            plt.close()
        
        self.count1_ch = count1_ch
        self.n_1det = n_1det  
    

    

    
    
    
    def coincidence_2det(self, plot = False):
        """
        Provides plot for the 2-det coincidence energy spectrum
        parameter:
            plot: True or False, saves .jpg if given True
            
        return: save plot if asked
        """
    
        Energy_2det=[]
        Channel_2det=[]
        Energy_2det_inner=[]
        count=0
        for l in self.Timed_Energy_event:
            if len(l)==2:
                E1=self.data[2][l[0]]
                E2=self.data[2][l[1]]
                E=np.array([E1, E2])
                Energy_2det.extend([E])
                p=[int(self.data[0][l[0]]), int(self.data[0][l[1]])]
                Channel_2det.extend([p])
                # E_inn=0
                # for m in range(2):
                #     if p[m] in innerCh:
                #         E_inn+=E[m]
                # if E_inn!=0:
                #     Energy_2det_inner.extend([E_inn])
                
        Energy_2det=np.array(Energy_2det, dtype=object)
        # Energy_2det_add=Energy_2det[:,0]+Energy_2det[:,1]
        
        
        
        N=300
        h=plt.hist2d(Energy_2det[:,0], Energy_2det[:,1], bins=(np.arange(0,self.energy_end,self.energy_end/N), np.arange(0,self.energy_end,self.energy_end/N)), norm = LogNorm())
        # h=plt.hist2d(Energy_2det_sl[0], Energy_2det_sl[1], bins=(np.arange(0,2500,2500/N), np.arange(0,2500,2500/N)), norm = LogNorm())
        plt.title("2 detector coincidence events distribution", fontsize=14)
        self.label("First hit detector Energy (keV)", "Second hit detector Energy (keV)", 8,5)
        if plot:
            plt.savefig("%s/2Det-DistribuionPlot.jpg"%(save_folder_each))
        plt.close()


    def save_pickle(self, energy = False):
        """
        save pickle file with the necessary information
        parameter:
  
        return: save .pickle file
        """
        f = open('%s/Output_Run%d-folder%d_index%d.pickle'%(save_folder_each,run_index, folder_index, index), 'wb')
        Parameter={'Run': run_index, 'Time': self.data_time, 'nEvent': self.nEvent,  "Multiplicity":self.Multiplicity, "Multiplicity_raw":self.Multiplicity_raw,
                   'n_EnergyInner':self.n_EnergyInner, 'n_EnergyAll':self.n_EnergyAll, "n_EnergyAll_NVeto":self.n_EnergyAll_NVeto, "Bins":self.Bins,     
                   "n_EnergyInner_NVeto":self.n_EnergyInner_NVeto, "n_1det":self.n_1det, "count1det_ch":self.count1_ch}
        pickle.dump(Parameter, f)
        f.close()
        if energy:
            f = open('%s/OutputEnergy_Run%d-folder%d_index%d.pickle'%(save_folder_each,run_index, folder_index, index), 'wb')
            Parameter={"antiveto_index":self.antiveto_index, "Timed_Energy_event":self.Timed_Energy_event, "data":self.data}
            pickle.dump(Parameter, f)
            f.close()

        

        
        






# files = glob.glob('%s/RAW/*.CSV'%read_folder)
# data = read_data(files)
# data = data_calibration(data, energy_source)

# analysis = ASCID_analysis(data)
# analysis.data_preprocess()
# analysis.coincidence()
# analysis.CoaddedEnergy()
# analysis.added_energy("All", plot = plottype)
# analysis.added_energy("inner", plot = plottype)
# analysis.multiplicity(plot = plottype, veto = "Veto")
# analysis.multiplicity(plot = plottype, veto = "noVeto")
# analysis.coincidence_1det(plot = plottype)
# analysis.coincidence_2det(plot = plottype)
# analysis.save_pickle()

            
        