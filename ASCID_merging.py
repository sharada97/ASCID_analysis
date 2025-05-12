from ASCID_variable import *
from ASCID_analysis import *
import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
import os
import pandas as pd
from scipy.optimize import curve_fit
import glob
import pickle
from matplotlib.colors import LogNorm


def gaus_fit(n, bins, ri,rf,fx,fy,nm="norm"):
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
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(fx, fy, textstr, fontsize=12,
             
             verticalalignment='top', bbox=props)
    return popt[0], popt[1],popt[2],plt.plot(X,gaus(X,*popt),'r')




class DataMerge:
    def __init__(self, files):
        self.files = files

        self.Data_time = None
        self.nEvent = None
        self.Multiplicity = None
        self.Multiplicity_raw = None
        self.n_EnergyInner = None
        self.n_EnergyAll = None
        self.n_EnergyAll_NVeto = None
        self.n_EnergyInner_NVeto = None
        self.Bins = None
        self.n_1det = None
        self.count1det_ch = None

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
        self.Data_time=dat['Time']
        self.nEvent=dat['nEvent']
        self.Multiplicity = dat['Multiplicity']
        self.Multiplicity_raw = dat['Multiplicity_raw']
        self.n_EnergyInner = dat['n_EnergyInner']
        self.n_EnergyAll = dat['n_EnergyAll']
        self.n_EnergyAll_NVeto = dat['n_EnergyAll_NVeto']
        self.Bins = dat['Bins']
        self.n_EnergyInner_NVeto = dat['n_EnergyInner_NVeto']
        self.n_1det = dat['n_1det']
        self.count1det_ch = dat['count1det_ch']
        f.close()

    def merging(self):
        all_files = len(self.files)
        for i in range(1, all_files):
            f = open('%s'%(self.files[i]), 'rb')
            dat=pickle.load(f)
            self.Data_time+= dat['Time']
            self.nEvent+=dat['nEvent']
            self.Multiplicity+= dat['Multiplicity']
            self.Multiplicity_raw+= dat['Multiplicity_raw']
            self.n_EnergyInner+= dat['n_EnergyInner']
            self.n_EnergyAll+= dat['n_EnergyAll']
            self.n_EnergyAll_NVeto+= dat['n_EnergyAll_NVeto']
            self.n_EnergyInner_NVeto += dat['n_EnergyInner_NVeto']
            self.n_1det += dat['n_1det']
            self.count1det_ch += dat['count1det_ch']
            f.close()
    

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
                width = 0.4, label="Total decay rate= %.1f decay/sec"%(total/(self.Data_time)))   #Change time here
        addlabels(courses, values)
        plt.legend(fontsize=14)
        self.label("No. of coincident detectors", "No. of events (percentage)", 8, 5)
        plt.title("Multiplicity Distribution (%s), %s,  %.2f-Hr data"%(source, veto, (self.Data_time/3600)), fontsize="15")
        if plot:
            plt.savefig("%s/Result/Multiplicity_%s.jpg"%(save_folder, veto))
        plt.close()

    
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
            EventEnergy = self.n_EnergyAll
            EventEnergy_NVeto = self.n_EnergyAll_NVeto
        elif group == "inner":
            EventEnergy = self.n_EnergyInner
            EventEnergy_NVeto = self.n_EnergyInner_NVeto
        plt.rcParams['figure.figsize'] = [12,5]
        plt.step(self.Bins[:-1], EventEnergy, 'b',where="mid",linewidth=1,label="%d events, All"%EventEnergy.sum())
        plt.step(self.Bins[:-1], EventEnergy_NVeto, 'orange',where="mid",linewidth=1,label="%d events, Veto AntiCoincidence"%EventEnergy_NVeto.sum())
        n = EventEnergy_NVeto
        bins = self.Bins[:-1]
        gaus_fit(n, bins, 800, 960, 100, 1e8,nm="keV")
        gaus_fit(n, bins, 1000, 1200, 1200, 1e8,nm="keV")
        plt.legend(loc="upper left", fontsize=15)
        self.label("Energy (keV)", "Counts", 12, 5)
        plt.title("%s detector coincidence Spectrum, %.2f-Hr data"%(group, self.Data_time/3600), fontsize=14)
        plt.grid()
        if plot:
            plt.savefig("%s/Result/Energy-%sCoin.jpg"%(save_folder, group))
        plt.close()

    def plot_1det(self, plot = False):
        plt.step(self.Bins[:-1], self.n_1det, 'k',where="mid",linewidth=1,label="%d events"%self.n_1det.sum())
        
        plt.legend(loc="upper left", fontsize=15)
        self.label("Energy (keV)", "Counts", 12, 5)
        plt.title("1 detector: Spectrum, %.2f-Hr data"%(self.Data_time/3600), fontsize=14)
        plt.grid()
        if plot:
            plt.savefig("%s/Energy-1Det.jpg"%(save_folder))
        plt.close()


        
        count1_ch=np.array(self.count1det_ch, dtype=int) 
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
        plt.title("1Det coincidence events, %.2f-Hr data"%(self.Data_time/3600), fontsize="15")
        if plot:
            plt.savefig("%s/Result/CountEvents-1det.jpg"%(save_folder))
        plt.close()
        
    def save_pickle(self):
        """
        save pickle file with the necessary information
        parameter:
  
        return: save .pickle file
        """
        f = open('%s/Result/Result_Run%d.pickle'%(save_folder,run_index), 'wb')
        Parameter={'Run': run_index, 'Time': self.Data_time, 'nEvent': self.nEvent,  "Multiplicity":self.Multiplicity, "Multiplicity_raw":self.Multiplicity_raw,
                   'n_EnergyInner':self.n_EnergyInner, 'n_EnergyAll':self.n_EnergyAll, "n_EnergyAll_NVeto":self.n_EnergyAll_NVeto, "Bins":self.Bins,     
                   "n_EnergyInner_NVeto":self.n_EnergyInner_NVeto, "n_1det":self.n_1det, "count1det_ch":self.count1det_ch}
        pickle.dump(Parameter, f)
        f.close()

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

# Example usage
save_folder = "../RUN5/Sc46/Figures/"
search_path = "%s"%save_folder
partial_name = "Output_Run5"

file_array = find_files_partial_match_in_folders(partial_name, search_path)
print(f"Files found: {len(file_array)}")

# files = glob.glob("%s/Output_Run5*.pickle"%save_folder) 
files = file_array
data_merge = DataMerge(files)
data_merge.load_data()
data_merge.merging()
data_merge.multiplicity(plot = True, veto = "Veto")
data_merge.multiplicity(plot = True, veto = "noVeto")
data_merge.added_energy("All", plot = True)
data_merge.added_energy("inner", plot = True)
data_merge.plot_1det(plot = True)
data_merge.save_pickle()


        

    
    