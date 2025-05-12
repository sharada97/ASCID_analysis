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
# from ASCID_variable import *



def read_data(files):
    N_det = len(files)
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


def find_peaks(x, y, channel, window_length, min_distance, folder):
    from scipy.signal import find_peaks, savgol_filter

    
    y_smooth = savgol_filter(y, window_length, polyorder=3)
    
    # Find the index where 40% of the data lies
    start_index = int(0.15 * len(x))
    
    # Subset the smoothed data to include only values after 40%
    x_subset = x[start_index:]
    y_subset_smooth = y_smooth[start_index:]
    
    # Find peaks in the smoothed subset of data (after 40%)
    # We will use the prominence criterion to select well-separated peaks
    min_prominence = 0.2  # Minimum prominence to filter out small fluctuations (adjust this)
    
    peaks, properties = find_peaks(y_subset_smooth, prominence=min_prominence, distance=min_distance)
    
    # Ensure we only pick the two most prominent peaks
    peak_prominences = properties['prominences']
    top_peaks_idx = np.argsort(peak_prominences)[-2:]  # Get indices of the two highest peaks
    
    # Get the x and y values of the top two highest peaks
    top_peaks = peaks[top_peaks_idx]
    x_peaks = x[start_index:][top_peaks]  # Adjust the peak indices to the original x-values
    y_peaks = y_smooth[start_index:][top_peaks]  # Adjust the peak indices to the smoothed y-values

    mean = np.sort(x_peaks)
    
    # Plot the original and smoothed data, and mark the detected peaks
    plt.plot(x, y, label='Original Data', alpha=0.5)
    plt.plot(x, y_smooth, label='Smoothed Data', color='orange', linewidth=2)
    plt.scatter(x_peaks, y_peaks, color='red', label='Detected Peaks', zorder=5)
    plt.xlabel("Energy (keV)",fontsize=18)
    plt.title("Energy spectrum - Channel-%d"%channel)
    plt.legend(fontsize = 14)
    plt.ylabel("counts",fontsize=18)
    plt.rcParams['figure.figsize'] = [10,5]
    plt.grid()
    plt.savefig("%s/EnergyFit_Channel_%d.jpg"%(folder, channel))
    plt.close()
    # Print the x and y values of the detected peaks
    #print(f"Detected top two prominent peaks at x-values: {x_peaks}")
    return mean


def find_1peaks(x, y, channel, window_length, folder):
    from scipy.signal import find_peaks
    from scipy.signal import find_peaks, savgol_filter
    y_smooth = savgol_filter(y, window_length, polyorder=3)
    
    # Find the minimum value of the spectrum
    min_value = np.min(x)
    max_value = np.max(x)
    
    # Define the threshold as 20% of the minimum value
    threshold = 0.2 * max_value + min_value 
    
    # Create a mask to select the data points where y is greater than 20% of the minimum value
    mask = x >= threshold
    
    # Subset the data starting from the point where the values are above the threshold
    x_subset = x[mask]
    y_subset = y_smooth[mask]
    
    # Find peaks in the subset of the data
    peaks, _ = find_peaks(y_subset, distance=10, height=0)  # Adjust distance to ensure separation
    
    # If no peaks are found, you can handle it accordingly
    if len(peaks) > 0:
        # Get the index of the single highest peak in the subset
        top_peak_idx = np.argmax(y_subset[peaks])
        peak_x = x_subset[peaks[top_peak_idx]]
        peak_y = y_subset[peaks[top_peak_idx]]
    else:
        peak_x = None
        peak_y = None
    
    # Print the results
    print(f"Peak found at x = {peak_x}, y = {peak_y}")
    
    textstr = r'$\mu=%d $ ADC' %(peak_x)                        #"r'$\chi^2/Dof=%.2f$' % (b, )"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(peak_x, peak_y/2, textstr, fontsize=12,
             verticalalignment='top', bbox=props)
    
    # Plot the original data, smoothed data, and the detected peak
    plt.plot(x, y, label='Original Data', alpha=0.5)
    plt.plot(x, y_smooth, label='Smoothed Data', color='orange', linewidth=2)
    if peak_x is not None:
        plt.scatter(peak_x, peak_y, color='red', label='Detected Peak', zorder=5)
    plt.legend()
    # plt.yscale("log")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Peak Detection, channel %d'%channel)
    # plt.show()
    plt.savefig("%s/Energy_Ch%d.jpg"%(folder,channel))
    plt.close()
    return peak_x


    
def save_fit(energyPlot):
    f = open('%s/DataInfo_Run.pickle'%save_folder, 'rb')
    dat = pickle.load(f)
    mean=dat['mean']
    f.close()


    Mean = np.zeros((26, 2))

    savetype = True
    for k in range(N_det_CsI):
        factor = Vague_CalFactor[source_index]
        bin_last=mean[k, 0]*factor
        bin_start = 4 if mean[k, 0] < 40 else 10
        bin_w = 1 if mean[k, 0] < 40 else 2
        n,bins=np.histogram(energyPlot[k], bins=np.arange(bin_start,bin_last,bin_w))
        bin_centers = bins[1:]

        min_distance = max(int((mean[k, 1] - mean[k, 0])/2) - 8, 4)
        window_length = int(mean[k, 0]/6)   #length of the window this and that side of mean
        # Mean[k, :] = find_peaks(bin_centers, n, k, window_length, min_distance, full_path)
        Mean[k, :] = find_1peaks(bin_centers, n, k, window_length, full_path)
        if Mean[k, 0] > mean[k, 0]*1.2 or Mean[k, 0] < mean[k, 0]*0.8:
            savetype = False
        
        print(f"fit channel {k} done, mean value found {Mean[k, 0], Mean[k, 1]}, given {mean[k, 0], mean[k, 1]}")
    f = open('%s/Calibration_parameter_folder%d_index%d.pickle'%(save_folder_each, folder_index, index), 'wb')
    Parameter={'Time': cut_data, 'Mean':Mean}
    pickle.dump(Parameter, f)
    f.close()

    if index == 0 and savetype == True:
        f = open('%s/DataInfo_Run.pickle'%(save_folder), 'wb')
        Parameter={'mean':Mean}
        pickle.dump(Parameter, f)
        f.close()
        

        


# files = glob.glob('%s/RAW/*.CSV'%read_folder)
# data = read_data(files)
# energyPlot = data_preprocess(data)
# save_fit(energyPlot)