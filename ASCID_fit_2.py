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


def dual_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return (
        A1 * np.exp(-0.5 * ((x - mu1) / sigma1)**2) +
        A2 * np.exp(-0.5 * ((x - mu2) / sigma2)**2)
    )



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
    def __init__(self, x, y, channel, peak_position_1st, peak_distance, mean, shift_width, peak_value, n_iteration, max_rightshift = -5):
        self.x = x
        self.y = y
        self.channel = channel
        self.mean = mean
        self.peak_position_1st = peak_position_1st
        self.peak_distance = peak_distance
        self.peak_value = peak_value
        self.shift_width = shift_width
        self.max_rightshift = max_rightshift
        self.max_shift = n_iteration
        self.best_range = None 
        self.Params = None
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

    def Best_range_two_peaks(self, window1, window2, check_ratio, check_minfit):
        """
        Find the best range for dual-Gaussian fit by calculating the minimum error over different ranges.
        
        Parameters:
            window1: Width of the range on the left of the first peak.
            window2: Width of the range on the right of the second peak.
            check_ratio: Threshold for sigma/mu ratio to filter poor fits.
            check_minfit: Minimum mean value to consider the fit valid.
            peak_distance: Estimated distance between the two peaks.
            
        Returns:
            best_range: Tuple containing the start and end points of the best range found.
            Params: Contains the fit parameters of the best dual-Gaussian fit.
        """

        sigm_mu = 0.07  # Typical PMT sigma/mu ratio.
        sigma_initial = self.peak_position_1st * sigm_mu
        best_residual = float('inf')
        Params = []
        start, end = (self.peak_position_1st - window1, self.peak_position_1st + self.peak_distance + window2)  # Default initial range
        closest_index_start = np.abs(self.x - start).argmin()
        closest_index_end = np.abs(self.x - end).argmin()
        best_range = [closest_index_start, closest_index_end]
        
        # print(self.peak_position_1st, window1, window2, self.peak_distance)
        for shift_index in range(self.max_rightshift, self.max_shift):
            # Adjust range to include both peaks
            peak_position_1 = self.peak_position_1st - self.shift_width * shift_index
            peak_position_2 = peak_position_1 + self.peak_distance
            start = max(0, peak_position_1 - window1)
            end = min(max(self.x) - 1, peak_position_2 + window2)
    
            # Ensure the range is valid
            if start >= end:
                continue  # Skip invalid ranges
    
            closest_index_start = np.abs(self.x - start).argmin()
            closest_index_end = np.abs(self.x - end).argmin()
    
            x_window = self.x[closest_index_start:closest_index_end]
            y_window = self.y[closest_index_start:closest_index_end]
    
            # Ensure there are enough data points for fitting
            if len(y_window) < 6:  # At least 6 points for fitting 6 parameters (2 Gaussians)
                continue
    
            # Initial guess for dual-Gaussian fit
            try:
                initial_guess = [
                    max(y_window), peak_position_1, sigma_initial,  # Peak 1
                    max(y_window) * 0.8, peak_position_2, sigma_initial  # Peak 2
                ]
    
                # Define dual-Gaussian function
                def dual_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
                    return (
                        A1 * np.exp(-0.5 * ((x - mu1) / sigma1)**2) +
                        A2 * np.exp(-0.5 * ((x - mu2) / sigma2)**2)
                    )
    
                # Fit and calculate residuals
                params, _ = curve_fit(dual_gaussian, x_window, y_window, p0=initial_guess)
                residuals = y_window - dual_gaussian(x_window, *params)
                residual_sum = np.sum(residuals**2)

                # print(x_window[0],  params[1] - 4*self.shift_width) 
                # Validate and update best fit
                sigma1_ratio = params[2] / params[1]
                sigma2_ratio = params[5] / params[4]
                if (
                    residual_sum < best_residual and
                    sigma1_ratio < check_ratio and sigma2_ratio < check_ratio and
                    params[1] > check_minfit and params[4] > check_minfit  
                    # and params[0] > 0.9* self.peak_value[0] and params[3] > 0.9* self.peak_value[1]
                ):
                    best_residual = residual_sum
                    best_range = [closest_index_start, closest_index_end]
                    Params = params
    
            except (ValueError, RuntimeError):
                continue  # Skip if fitting fails for this range
    
        self.best_range = best_range
        self.Params = Params
        return best_range, Params

    
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
        # print(self.best_range)
        def dual_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
                    return (
                        A1 * np.exp(-0.5 * ((x - mu1) / sigma1)**2) +
                        A2 * np.exp(-0.5 * ((x - mu2) / sigma2)**2)
                    )
    
                # Fit and calculate residuals
        
        initial_guess = [*self.Params]
        x_window = self.x[self.best_range[0]: self.best_range[1]]
        y_window = self.y[self.best_range[0]: self.best_range[1]]

        # if self.channel<16:
        #     x_window = self.x[self.best_range[0]: self.best_range[1]]
        #     y_window = self.y[self.best_range[0]: self.best_range[1]]
        params, _ = curve_fit(dual_gaussian, x_window, y_window, p0=initial_guess)
        
        # params = self.Params
        if plot:
            plt.plot(self.x, self.y, label='Data')
            # plt.plot(x_window, y_window, label='fit_Data')
            plt.plot(
                    x_window, 
                    dual_gaussian(x_window, *params), 
                    label=f'mean1 = {params[1]:.2f}, sigma1 = {params[2]:.2f} \nmean2 = {params[4]:.2f}, sigma2 = {params[5]:.2f}', 
                    linestyle='--')
            plt.xlabel("Energy (keV)",fontsize=18)
            plt.title("Energy spectrum - Channel-%d"%self.channel)
            plt.legend(fontsize = 14)
            plt.ylabel("counts",fontsize=18)
            plt.rcParams['figure.figsize'] = [10,5]
            plt.grid()
            plt.savefig("%s/EnergyFit_Channel_%d.jpg"%(folder, self.channel))
            plt.close()
        # return params[1], params[2], params[4], params[5]



    
def save_fit(energyPlot):
    f = open('%s/DataInfo_Run.pickle'%save_folder, 'rb')
    dat = pickle.load(f)
    mean=dat['mean']
    amp=dat['amp']
    bin_last = dat['bin_last']
    good_range = dat["best_range"]
    f.close()


    Mean = np.zeros((26, 2))
    Sigma = np.zeros((26, 2))
    Amp = np.zeros((26, 2))
    Range = np.zeros((2, 26))

    for k in range(N_det_CsI):
        n,bins=np.histogram(energyPlot[k], bins=np.arange(12,bin_last[k],2))
        
        
        window1 = mean[k, 0] - good_range[0, k]
        window2 = good_range[1, k] - mean[k, 1]
        peak_position_1st = mean[k, 0] +window1/2
        peak_distance = mean[k, 1] - mean[k, 0]
        check_ratio = 0.15  #maximum sigma/mu ratio can be 20 %, in the worst case
        check_minfit = bin_last[k]*0.3  # 1st peak should start at least 30% of the total bins
        shift_width = (1 if mean[k, 0]<50 else 2 if 50<mean[k, 0]<80 else 3 if 80<mean[k, 0]<100 else 4)
        n_iteration = 7
        peak_value = amp[k, :]
        
        x = bins[:-1]
        y = n
        # print(peak_value)
        
        fit = bestFit(x, y, k, peak_position_1st, peak_distance, mean, shift_width, peak_value, n_iteration+5, max_rightshift = 0)
        BBest_range, PParams = fit.Best_range_two_peaks(window1, window2, check_ratio, check_minfit)
        fit.plot_fit(plot = True, folder = full_path, source = "Sc46")

        Amp[k, 0], Mean[k, 0], Sigma[k, 0], Amp[k, 1], Mean[k, 1], Sigma[k, 1] = PParams[0], PParams[1], PParams[2], PParams[3], PParams[4], PParams[5]
        Range[0, k], Range[1, k] = x[BBest_range[0]], x[BBest_range[1]]
        
        print(f"fit channel {k} done, best range {Range[:, k]}")
    f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'wb')
    Parameter={'Time': cut_data, 'Mean':Mean, 'Sigma': Sigma}
    pickle.dump(Parameter, f)
    f.close()

    if index == 0:
        f = open('%s/DataInfo_Run.pickle'%(save_folder), 'wb')
        Parameter={'mean':Mean, 'sigma': Sigma, "amp":Amp, 'bin_last':bin_last, 'best_range':Range}
        pickle.dump(Parameter, f)
        f.close()
        

        


# files = glob.glob('%s/RAW/*.CSV'%read_folder)
# data = read_data(files)
# energyPlot = data_preprocess(data)
# save_fit(energyPlot)
