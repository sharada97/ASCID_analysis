import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from natsort import natsorted
import sys
import os
import pandas as pd
import glob
import pickle
import shutil
import matplotlib.pyplot as plt
from ASCID_variable import *






def read_data(files):

    files=natsorted(files)
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
    
    return data

def data_preprocess(data):
    nEvents = len(data[0])
    energyPlot=[[] for _ in range(26)]
    for i in range(nEvents):
        if data[0][i] < 26:
            energyPlot[int(data[0][i])].extend([data[2][i]])
    return energyPlot


class bestFit:
    # Initialize the car with make, model, and year
    def __init__(self, x, y, channel, peak_position_1st, mean, shift_width, n_iteration, max_rightshift = -5):
        self.x = x
        self.y = y
        self.channel = channel
        self.mean = mean
        self.peak_position_1st = peak_position_1st
        self.shift_width = shift_width
        self.max_rightshift = max_rightshift
        self.max_shift = n_iteration
        self.best_range = None 
        self.energyPlot = None
        
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

    

    
    def gaussian(self, x, amplitude, mean, sigma):
        """Gaussian function for peak fitting."""
        return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

    def Best_range (self, window1, window2, check_ratio, check_minfit):
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
            best_range: contains the start and end point of the best range found
            Params: contains the fit parameters of the best fit
        """
        
        
    
        sigm_mu = 0.07 # typical PMT sigma/mu ratio taken. which is 7 %
        sigma_initial = self.peak_position_1st* sigm_mu
        best_residual = float('inf')
        Params = []
        best_range = (self.peak_position_1st - window1, self.peak_position_1st + 2 * window2)  # Default initial range
        for shift_index in range(self.max_rightshift, self.max_shift):
            peak_position = self.peak_position_1st - self.shift_width*shift_index
            start = max(0, peak_position - window1)
            end = min(max(self.x) - 1, peak_position + window2)

            # print(start, end)
            # Ensure the range is valid
            if start >= end:
                continue  # Skip if the range is invalid
        
            closest_index_start = np.abs(self.x - start).argmin()
            closest_index_end = np.abs(self.x - end).argmin()

            x_window = self.x[closest_index_start:closest_index_end]
            y_window = self.y[closest_index_start:closest_index_end]


            # print(x_window, y_window)
            # Ensure there are enough data points for fitting
            if len(y_window) < 3:  # At least 3 points are needed to fit 3 parameters
                continue
        
            # print("okay")
            # Initial guess for Gaussian fit
            try:
                initial_guess = [max(y_window), x_window[np.argmax(y_window)], sigma_initial]
                
                # Fit and calculate residuals
                params, _ = curve_fit(self.gaussian, x_window, y_window, p0=initial_guess)
                residuals = y_window - self.gaussian(x_window, *params)
                residual_sum = np.sum(residuals**2)
                # print(f"Mean, max_mean{params[1],self.mean[self.channel, 0], self.max_rightshift*self.shift_width}")
                # Update if this range provides a better fit
                if residual_sum < best_residual and params[2]/params[1]<check_ratio and params[1]>check_minfit:
                    best_residual = residual_sum
                    best_range = [closest_index_start, closest_index_end]
                    Params = params
            except (ValueError, RuntimeError):
                continue  # Skip if fitting fails for this range
        # print(best_range)
        self.best_range =  best_range
        # print(f"best range is from {best_range}")
        
    
    def plot_fit(self, plot = False, folder=None, source = "Sc46"):
        """
        Plotting the spectrum by taking the best range from the last function.
        Parameters:
            x: x- value, energy bin
            y: counts
            best_range: best found range from the last function
            plot: False or True
            source = "Sc46" or "Na22" or "Background"
    
        returns:
            plot: if plot = True is given
            params: [amplitude, mean, sigma], [amplitude, mean, sigma] of 2 peaks that fitted. 
        """
        sigma_initial = self.peak_position_1st*0.07
        x_window = self.x[self.best_range[0]: self.best_range[1]]
        y_window = self.y[self.best_range[0]: self.best_range[1]]
        initial_guess = [max(y_window), x_window[np.argmax(y_window)], sigma_initial]
        params, _ = curve_fit(self.gaussian, x_window, y_window, p0=initial_guess)
        
    
        start = self.x[self.best_range[0]]*(1.27 if source == "Sc46" else 1)  #mu2/mu1 ratio
        end = self.x[self.best_range[1]]*(1.27 if source == "Sc46" else 1)
        # sigma_initial = Params[1]*0.07
        closest_index_start = np.abs(self.x - start).argmin()
        closest_index_end = np.abs(self.x - end).argmin() 
        x_window2 = self.x[closest_index_start:closest_index_end]
        y_window2 = self.y[closest_index_start:closest_index_end]
        initial_guess = [max(y_window2), x_window2[np.argmax(y_window2)], sigma_initial]
        try:
            params2, _ = curve_fit(self.gaussian, x_window2, y_window2, p0=initial_guess)
        except:
            print(f"2nd peak fitting not done, Channel {self.channel}")
            params2 = [None, None, None]

        
        if plot:
            plt.plot(self.x, self.y, label='Data')
            plt.plot(x_window, self.gaussian(x_window, *params), label=f'mean = {params[1]:.2f}, sigma = {params[2]:.2f}', linestyle='--')
            try:
                plt.plot(x_window2, self.gaussian(x_window2, *params2), label=f'mean = {params2[1]:.2f}, sigma = {params2[2]:.2f}', linestyle='--')
            except:
                pass
            plt.xlabel("Energy (keV)",fontsize=18)
            plt.title("Energy spectrum - Channel-%d"%self.channel)
            plt.legend(fontsize = 14)
            plt.ylabel("counts",fontsize=18)
            plt.rcParams['figure.figsize'] = [10,5]
            plt.grid()
            plt.savefig("%s/EnergyFit_Channel_%d.jpg"%(folder, self.channel))
            plt.close()
        return params[1], params[2], params2[1], params2[2]



    
def save_fit(energyPlot):
    f = open('%s/DataInfo_Run.pickle'%save_folder, 'rb')
    dat = pickle.load(f)
    mean=dat['mean']
    sigma=dat['sigma']
    bin_last = dat['bin_last']
    gaus_x = dat["gaus_x"]
    f.close()


    Mean = np.zeros((26, 2))
    Sigma = np.zeros((26, 2))

    for k in range(N_det_CsI):
        n,bins=np.histogram(energyPlot[k], bins=np.arange(12,bin_last[k],2))
        
        
        window1 = mean[k, 0] - gaus_x[0][k]
        window2 = gaus_x[1][k] - mean[k, 0]
        peak_position_1st = mean[k, 0] +window1/2
        check_ratio = 0.2  #maximum sigma/mu ratio can be 20 %, in the worst case
        check_minfit = bin_last[k]*0.3  # 1st peak should start at least 30% of the total bins
        shift_width = (1 if mean[k, 0]<50 else 2 if 50<mean[k, 0]<80 else 3 if 80<mean[k, 0]<100 else 4)
        n_iteration = 10
        
        x = bins[:-1]
        y = n

        fit = bestFit(x, y, k, peak_position_1st, mean, shift_width, n_iteration+5, max_rightshift = 0)
        fit.Best_range(window1, window2, check_ratio, check_minfit)
        Mean[k, 0], Sigma[k, 0], Mean[k, 1], Sigma[k, 1] = fit.plot_fit(plot = True, folder = full_path, source = "Sc46")
        if Sigma[k, 0]/Mean[k, 0] > 0.1:
            fit = bestFit(x, y, k, peak_position_1st, mean, shift_width, n_iteration, max_rightshift = -5)
            fit.Best_range(window1, window2, check_ratio, check_minfit)
            Mean[k, 0], Sigma[k, 0], Mean[k, 1], Sigma[k, 1] = fit.plot_fit(plot = True, folder = full_path, source = "Sc46")
            if Sigma[k, 1]/Mean[k, 1] > 0.1:
                print("wrong fit, Channel %d"%k)

        
        
        print(f"fit channel {k} done")
    f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'wb')
    Parameter={'Time': cut_data, 'Mean':Mean, 'Sigma': Sigma}
    pickle.dump(Parameter, f)
    f.close()

    if index == 0:
        f = open('%s/DataInfo_Run.pickle'%(save_folder), 'wb')
        Parameter={'mean':Mean, 'sigma': Sigma, 'bin_last':bin_last, 'gaus_x':gaus_x}
        pickle.dump(Parameter, f)
        f.close()
        

        


# files = glob.glob('%s/RAW/*.CSV'%read_folder)
# data = read_data(files)
# energyPlot = data_preprocess(data)
# save_fit(energyPlot)
