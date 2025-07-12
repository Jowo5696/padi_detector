import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit

def drawHistograms(file, histogramArray, indexTranslation):
    for i in histogramArray:

        mm_strip = file[i]["raw"]["mm_strip"].array(library="np")
        apv_id = file[i]["raw"]["apv_id"].array(library="np")
        event_id = file[i]["data"]["apv_qmax"].array(library="np")
        mm_strip = np.concatenate(mm_strip)
        apv_id = np.concatenate(apv_id)
        event_id = np.concatenate(event_id)
        
        event_counts = Counter(event_id)


        #filter out unprecise data to eliminate local jitters
        valid_events = {eid for eid, count in event_counts.items() if count < 800}
        mask = np.isin(event_id, list(valid_events))

        mm_strip = mm_strip[mask]
        apv_id = apv_id[mask]

        #divide mm strip into x and y data with apv_id
        x_strip = mm_strip[(apv_id == 0) | (apv_id == 1)]
        y_strip = mm_strip[(apv_id == 2) | (apv_id == 3)]

        #swap out the left and right side to account for weird apv placement
        y_strip[y_strip >= 126] -= 250
        y_strip[y_strip >= 1] += 125
        y_strip[y_strip<=0] += 125 
        
        #convert the bin count into mm
        y_strip = y_strip / 2.5  +1

        #these are now arrays that list the individual counts, such that the size of this array is the amount of data collected.
        #I will transform this into two arrays, both 250 long (the amount of buckets)
        #such that one records the size of the bucket, and the other the bucket number.
        counts, edges = np.histogram(y_strip, bins=250, range=(1,100))
        centers = 0.5 * (edges[:-1] + edges[1:])
        
        #mirroring the Array so that the data becomes easier to analyze and the sensitivity difference between apvs is ignored.
        counts = mirrorArrayLeft(counts)





        #Now start a gauss fit by doing an initial guess
        #[amplitude, mean, stddev]
        p0 = [np.max(counts), centers[np.argmax(counts)], 10, 500]
        
        #fit with scipy
        popt, pcov = curve_fit(gauss, centers, counts, p0=p0)


        plt.bar(centers, counts, width=edges[1] - edges[0], alpha=0.6, label="Data")
        plt.plot(centers, gauss(centers, *popt), color='red', label="Gaussian Fit")
        
        plt.xlabel("x [mm]")
        plt.ylabel("Events")
        plt.text(60,np.max(counts) - np.max(counts)/10, f"μ={popt[1]:.2f}, σ={popt[2]:.2f}", fontsize=15)
        plt.text(60,np.max(counts) - 1.3*np.max(counts)/10, "θ=" + str(np.round(np.arctan(popt[2]/400), 5)), fontsize=15)
        plt.title(indexTranslation[i])
        plt.savefig("finished_plots/" + indexTranslation[i] + ".png")



        plt.show()
        


def openfile():
    file = [
            up.open("../../data/elsa/no_target/run82.root"),
            up.open("../../data/elsa/half_radiation_length/aluminium_5cm/distance_40cm/run87.root"),
           up.open("../../data/elsa/radiation_length/aluminium_9cm93/distance_40cm/run88.root"),
            up.open("../../data/elsa/double_radiation_length/aluminium_21cm/distance_40cm/run85.root"),
            up.open("../../data/elsa/triple_radiation_length/aluminium_27cm39/distance_40cm/run93.root"),

            up.open("../../data/elsa/half_radiation_length/copper_0cm7/distance_40cm/run91.root"),
            up.open("../../data/elsa/radiation_length/copper_1cm48/distance_40cm/run90.root"), 
            up.open("../../data/elsa/double_radiation_length/copper_2cm96/distance_40cm/run89.root"),
            up.open("../../data/elsa/triple_radiation_length/copper_4cm36/distance_40cm/run94.root")
            ]
    return file

def mirrorArrayLeft(Array):
    half = len(Array) // 2
    return np.concatenate([Array[:half], Array[half-1::-1]])

def gauss(x, A, mu, sigma, b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


noMaterialIndex = 0
AluminiumIndex = [1,2,3,4]
CopperIndex = [5,6,7,8]
#first index is half radiation length, last is triple.

indexTranslation = [
    "no target",
    "Aluminium, Half Radiation Length, 40cm Distance",
    "Aluminium, One Radiation Length, 40cm Distance",
    "Aluminium, Two Radiation Lengths, 40cm Distance",
    "Aluminium, Three Radiation Lengths, 40cm Distance",
    "Copper, Half Radiation Length, 40cm Distance",
    "Copper, One Radiation Length, 40cm Distance",
    "Copper, Two Radiation Lengths, 40cm Distance",
    "Copper, Three Radiation Lengths, 40cm Distance"
        ]
files = openfile()

#Use this array to show which histograms you want to analyze and plot
histogramArray = [0,1,2,3,4,5,6,7,8]
drawHistograms(files, histogramArray, indexTranslation)
