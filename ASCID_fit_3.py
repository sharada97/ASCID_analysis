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
        if int(data[0][i]) < 26:
            energyPlot[int(data[0][i])].extend([data[2][i]])
    return energyPlot


class bestFit:
    # Initialize the car with make, model, and year
    def __init__(self, x, y, channel, peak_position_1st, peak_position_2nd, mean, shift_width, n_iteration, max_rightshift = -3):
        self.x = x
        self.y = y
        self.channel = channel
        self.mean = mean
        self.peak_position_1st = peak_position_1st
        self.peak_position_2nd = peak_position_2nd
        # self.peak_distance = peak_distance
        # self.peak_value = peak_value
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

    def fit_two_gaussians(self, window1, window2, check_ratio, check_minfit):
        """
        Fit two Gaussian peaks with custom ranges for each peak.
        Parameters:
            window1: List [a, b], defines left and right range around the 1st peak.
            window2: List [c, d], defines left and right range around the 2nd peak.
            check_ratio: Maximum acceptable ratio sigma/mean for Gaussian fits.
            check_minfit: Minimum acceptable mean value for the fit.
        Returns:
            params1, params2: Fitted parameters for the 1st and 2nd Gaussian peaks.
        """
        best_residual1, best_residual2 = float('inf'), float('inf')
        params1, params2 = None, None
        fit_range1 = np.inf*np.ones(10)
        # Fitting the first peak
        for shift_index in range(self.max_rightshift, self.max_shift, 1):
            peak_position = self.peak_position_1st - self.shift_width * shift_index
            start = max(0, peak_position - window1[0])  # window1[0] for left range
            end = min(max(self.x) - 1, peak_position + window1[1])  # window1[1] for right range
            
            closest_index_start = np.abs(self.x - start).argmin()
            closest_index_end = np.abs(self.x - end).argmin()
            
            x_window = self.x[closest_index_start:closest_index_end]
            y_window = self.y[closest_index_start:closest_index_end]
            
            
            if len(y_window) < 3:  # Ensure enough data points for fitting
                continue
            
            # Initial guess and fit
            try:
                initial_guess = [max(y_window), peak_position, 0.1 * peak_position]  # amplitude, mean, sigma
                params, _ = curve_fit(self.gaussian, x_window, y_window, p0=initial_guess)
                residuals = y_window - self.gaussian(x_window, *params)
                residual_sum = np.sum(residuals**2)
                
                if (
                    residual_sum < best_residual1
                    and params[2] / params[1] < check_ratio
                    and params[1] > check_minfit
                    and min(x_window) < params[1] - 2 * self.shift_width
                ):
                    best_residual1 = residual_sum
                    params1 = params
                    fit_range1 = x_window
            except (ValueError, RuntimeError):
                continue
        
        # Fitting the second peak
        fit_range2 = np.inf*np.ones(10)
        for shift_index in range(self.max_rightshift, self.max_shift, 1):
            peak_position = self.peak_position_2nd - self.shift_width * shift_index
            start = max(0, peak_position - window2[0])  # window1[0] for left range
            end = min(max(self.x) - 1, peak_position + window2[1])  # window1[1] for right range
            
            closest_index_start = np.abs(self.x - start).argmin()
            closest_index_end = np.abs(self.x - end).argmin()
            
            # print(peak_position, start, end, self.peak_position_2nd, shift_index, self.max_rightshift, self.max_shift)
            x_window = self.x[closest_index_start:closest_index_end]
            y_window = self.y[closest_index_start:closest_index_end]

            if  max(fit_range1) < min(x_window):
                continue
            
            if len(y_window) < 3:  # Ensure enough data points for fitting
                continue
            
            # Initial guess and fit
            try:
                initial_guess = [max(y_window), peak_position, 0.1 * peak_position]  # amplitude, mean, sigma
                params, _ = curve_fit(self.gaussian, x_window, y_window, p0=initial_guess)
                residuals = y_window - self.gaussian(x_window, *params)
                residual_sum = np.sum(residuals**2)
                
                if (residual_sum < best_residual2 and params[2] / params[1] < check_ratio and params[1] > check_minfit and min(x_window) < params[1] - 2 * self.shift_width):
                    best_residual2 = residual_sum
                    params2 = params
                    fit_range2 = x_window
            except (ValueError, RuntimeError):
                continue
    
        # Return the best fit parameters for both peaks
        return fit_range1, fit_range2, params1, params2



    def plot_gaussian_fits(self, fit_range1, fit_range2, params1, params2, plot = False, folder=None, source = "Sc46"):
        """
        Plot the data along with the Gaussian fits for both peaks.
        Parameters:
            params1: Parameters for the 1st Gaussian fit [amplitude, mean, sigma].
            params2: Parameters for the 2nd Gaussian fit [amplitude, mean, sigma].
        """
        # Define Gaussian function
        x_window = fit_range1
        x_window2 = fit_range2
        if plot:
            plt.plot(self.x, self.y, label='Data')
            plt.plot(x_window, self.gaussian(x_window, *params1), label=f'mean = {params1[1]:.2f}, sigma = {params1[2]:.2f}', linestyle='--')
            plt.plot(x_window2, self.gaussian(x_window2, *params2), label=f'mean = {params2[1]:.2f}, sigma = {params2[2]:.2f}', linestyle='--')
            plt.xlabel("Energy (keV)",fontsize=18)
            plt.title("Energy spectrum - Channel-%d"%self.channel)
            plt.legend(fontsize = 14)
            plt.ylabel("counts",fontsize=18)
            plt.rcParams['figure.figsize'] = [10,5]
            plt.grid()
            plt.savefig("%s/EnergyFit_Channel_%d.jpg"%(folder, self.channel))
            plt.close()
    
    # Example usage
    




    
def save_fit(energyPlot):
    f = open('%s/DataInfo_Run.pickle'%save_folder, 'rb')
    dat = pickle.load(f)
    mean=dat['mean']
    peak1 = dat['windows']['1']
    peak2 = dat['windows']['2']
    f.close()


    Mean = np.zeros((26, 2))
    Sigma = np.zeros((26, 2))

    Range1 = np.zeros((2, 26))
    Range2 = np.zeros((2, 26))

    for k in range(N_det_CsI):
        bin_last=mean[k, 0]*1.8
        n,bins=np.histogram(energyPlot[k], bins=np.arange(12,bin_last,2))
        
        
        window1 = [peak1[k, 0]*mean[k, 0], peak1[k, 1]*mean[k, 0]]
        window2 = [peak2[k, 0]*mean[k, 0], peak2[k, 1]*mean[k, 0]]
        peak_position_1st = mean[k, 0]  
        peak_position_2nd = mean[k, 1]

        # print(peak_position_1st, peak_position_2nd, window1, window2)
        
        check_ratio = 0.15  #maximum sigma/mu ratio can be 20 %, in the worst case
        check_minfit = bin_last*0.4  # 1st peak should start at least 30% of the total bins
        shift_width = int(mean[k, 0]/50)*0.75+1
        n_iteration = 14
                
        x = bins[:-1]
        y = n
        
        fit = bestFit(x, y, k, peak_position_1st, peak_position_2nd, mean, shift_width, n_iteration, max_rightshift = -7)
        best_range1, best_range2, params1, params2 = fit.fit_two_gaussians(window1, window2, check_ratio, check_minfit)
        fit.plot_gaussian_fits(best_range1, best_range2, params1, params2, plot = True, folder = full_path, source = "Sc46")

        Amp, Mean[k, 0], Sigma[k, 0], Amp2, Mean[k, 1], Sigma[k, 1] = params1[0], params1[1], params1[2], params2[0], params2[1], params2[2]
        # print(best_range1)
        Range1[0, k], Range1[1, k] = min(best_range1), max(best_range1)
        Range2[0, k], Range2[1, k] = min(best_range2), max(best_range2)
        
        print(f"fit channel {k} done, best range found {Range1[:, k], Range2[:, k]}, given {window1, window2}")
    f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'wb')
    Parameter={'Time': cut_data, 'Mean':Mean, 'Sigma': Sigma}
    pickle.dump(Parameter, f)
    f.close()

    if index == 0:
        f = open('%s/DataInfo_Run.pickle'%(save_folder), 'wb')
        Parameter={'mean':Mean, 'windows':{"1":peak1, "2":peak2}}
        pickle.dump(Parameter, f)
        f.close()
        

        


# files = glob.glob('%s/RAW/*.CSV'%read_folder)
# data = read_data(files)
# energyPlot = data_preprocess(data)
# save_fit(energyPlot)
