import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit

fixmu = 48

def drawHistograms(file, histogramArray, indexTranslation):
    
    #To later calculate the standard deviation from the difference of the no data run and the other runs.
    standardDevNoData = 0
    standardDevNoDataErr = 0
    for i in histogramArray:

        mm_strip = file[i]["raw"]["mm_strip"].array(library="np")
        apv_id = file[i]["raw"]["apv_id"].array(library="np")
        event_id = file[i]["data"]["apv_qmax"].array(library="np")
        
        events = file[i]["raw"]["apv_evt"].array(library="np")
        print(events)
        print(len(events))
        #plt.hist(events)
        #plt.show()


        mm_strip = np.concatenate(mm_strip)
        apv_id = np.concatenate(apv_id)
        event_id = np.concatenate(event_id)
        
        event_counts = Counter(event_id)


        #filter out unprecise data to eliminate local jitters
        valid_events = {eid for eid, count in event_counts.items() if count < 700}
        mask = np.isin(event_id, list(valid_events))

        mm_strip = mm_strip[mask]
        apv_id = apv_id[mask]

        #divide mm strip into x and y data with apv_id
        x_strip = mm_strip[(apv_id == 0) | (apv_id == 1)]
        y_strip = mm_strip[(apv_id == 2) | (apv_id == 3)]
        
        #doesn't work because the lengths aren't equal
        #print(len(x_strip))
        #print(len(y_strip))
        #if i == 0:
        #    twodHist(x_strip,y_strip)
        
        #if i == 0:
        #    plt.hist(y_strip, bins=250, range=(1,250))
        #    plt.xlabel("y-Detector")
        #    plt.ylabel("events")
        #    plt.show()



        #swap out the left and right side to account for weird apv placement
        #x_strip[x_strip >= 128] -= 256
        #x_strip[x_strip >= 1] += 127
        #y_strip[y_strip <=0] += 127

        y_strip[y_strip > 128] -= 256
        y_strip[y_strip >= 1] += 128
        y_strip[y_strip<=0] += 128

        if i == 0 and False:
            plt.hist(y_strip, bins=250, range=(1,250))
            plt.xlabel("y-Detector")
            plt.ylabel("events")
            plt.savefig("finished_plots/unfiltered_noMaterial.png")
            plt.show()
        
        #convert the bin count into mm
        y_strip = y_strip / 2.54
        x_strip = x_strip / 2.54
        
        #plt.hist(x_strip, bins = 250, range=(1,100))
        #plt.xlabel("x-Position")
        #plt.ylabel("Counts")

        #plt.show()

        #if i == 0:
        #    twodHist(x_strip, y_strip)

        #these are now arrays that list the individual counts, such that the size of this array is the amount of data collected.
        #I will transform this into two arrays, both 250 long (the amount of buckets)
        #such that one records the size of the bucket, and the other the bucket number.
        ycounts, yedges = np.histogram(y_strip, bins=254, range=(1,100))
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])
        
        #mirroring the Array so that the data becomes easier to analyze and the sensitivity difference between apvs is ignored.
        #yMirrorcounts = mirrorArrayLeft(ycounts)
        #for k in range(round(len(yMirrorcounts)/2), len(yMirrorcounts)):
        #    yMirrorcounts[k]=0
        yMirrorcounts = ycounts
        #yMirrorcounts = mirrorArrayRight(ycounts)

        #Now start a gauss fit by doing an initial guess
        #[amplitude, mean, stddev]
        #p0 = [np.max(yMirrorcounts), ycenters[np.argmax(yMirrorcounts)], 10]
        p0 = [np.max(yMirrorcounts), 10]

        
        #fit with scipy
        if False: 
            popt, pcov = curve_fit(gaussMu, ycenters, yMirrorcounts, p0=p0)
            perr = np.sqrt(np.diag(pcov))
            print(perr)
        else: 
            fitMask = (ycenters >= 40) & (ycenters <= 47)
            fitx = ycenters[fitMask]
            fity = ycounts[fitMask]
            
            popt, pcov = curve_fit(gaussMu, fitx, fity, p0=p0)
            perr = np.sqrt(np.diag(pcov))
            print(perr)

            
        if False:
            if i == 0:
                standardDevNoData = popt[2] #Saves noData standard deviation
        else:
            if i == 0:
                standardDevNoData = popt[1]
                standardDevNoDataErr = perr[1]


        plt.bar(ycenters, yMirrorcounts, width=yedges[1] - yedges[0], alpha=1, label="Data")
        plt.plot(ycenters, gaussMu(ycenters, *popt), color='red', label="Gaussian Fit")
        
        plt.xlabel("y-Position [mm]")
        plt.ylabel("Counts")
       
        if False:
            if i == 0:
                
                plt.text(60,np.max(yMirrorcounts) - np.max(yMirrorcounts)/10, rf"$\mu={popt[1]:.2f} \pm {perr[1]:.3f}$, $\sigma={popt[2]:.2f} \pm {perr[2]:.3f}$", fontsize=15)
                p = np.round(np.arctan(popt[2]/450), 5)
                thetaErr = np.sqrt(
                        (1/(1-(popt[1]/450)**2))**2 * perr[1]**2 +
                        (1/(1-(popt[1]/450)**2)) * 20**2
                        )
                plt.text(60,np.max(yMirrorcounts) - 3*np.max(yMirrorcounts)/10, rf"$\theta= {p:.5f}$", fontsize=15)
            else:
                popt[2] = np.sqrt( popt[2]**2 - standardDevNoData**2 )
                plt.text(60,np.max(yMirrorcounts) - np.max(yMirrorcounts)/10, rf"$\mu={popt[1]:.2f}$, $\sigma={popt[2]:.2f}$", fontsize=15)
                p = np.round(np.arctan(popt[2]/450),5)
                thetaErr = np.sqrt(
                        (1/(1-(popt[1]/450)**2))**2 * perr[1]**2 +
                        (1/(1-(popt[1]/450)**2)) * 20**2
                        )
                plt.text(60, np.max(yMirrorcounts) - 3*np.max(yMirrorcounts)/10, rf"$\theta= {p:.5f}$", fontsize=15)
        else:
            if i == 0:
                plt.text(60, np.max(yMirrorcounts) - np.max(yMirrorcounts)/10, rf"$\mu= ={fixmu:.2f}$, $\sigma= {popt[1]:.2f}\pm {perr[1]:.3f}$", fontsize=10)
                p = np.round(np.arctan(popt[1]/450), 5)
                thetaErr = np.sqrt(
                        ((1/450)*(1/(1-(popt[1]/450)**2)))**2 * perr[1]**2 +
                        ((popt[1]/(450**2))*(1/(1-(popt[1]/450)**2)))**2 * 20**2
                        )
                plt.text(60, np.max(yMirrorcounts) -3*np.max(yMirrorcounts)/10, rf"$\theta = {p:.5f} \pm {thetaErr:.5f} $" , fontsize=10)
            else:
                popt[1] = np.sqrt( popt[1]**2 - standardDevNoData**2 )
                perr[1] = np.sqrt(
                            ((2*popt[1])/(np.sqrt(
                                    popt[1]**2 - standardDevNoData**2
                                )))**2 * perr[1]**2 +
                            ((2*standardDevNoData**2)/(np.sqrt(
                                    popt[1]**2 - standardDevNoData**2
                                )))**2 * standardDevNoDataErr**2
                        )
                thetaErr = np.sqrt(
                        ((1/450)*(1/(1-(popt[1]/450)**2)))**2 * perr[1]**2 +
                        ((popt[1]/(450**2))*(1/(1-(popt[1]/450)**2)))**2 * 20**2
                        )
                p = np.round(np.arctan(popt[1]/450), 5)
                plt.text(60, np.max(yMirrorcounts) - np.max(yMirrorcounts)/10, rf"$\mu = {fixmu:.2f}$, $\sigma={popt[1]:.2f} + {perr[1]:.2f}$", fontsize=10)
                np.round(np.arctan(popt[1]/450), 5)                                 
                plt.text(60, np.max(yMirrorcounts) -3*np.max(yMirrorcounts)/10, rf"$\theta={p:.5f} \pm {thetaErr:.5f}$", fontsize=10)


        #plt.title(indexTranslation[i])
        plt.savefig("finished_plots/" + indexTranslation[i] + ".png")



        plt.show()
        
def twodHist(xHist, yHist):
    xbins, ybins = 250, 250
    xrange, yrange = (1,100), (1,100)

    H, xedges, yedges = np.histogram2d(
            xHist, yHist, bins=[xbins, ybins], 
            range=[xrange, yrange])


    plt.figure(figsize=(6, 5))
    plt.imshow(H.T, origin="lower", aspect="auto",
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           cmap="viridis")

    plt.colorbar(label="Counts")
    plt.xlabel("X Strip")
    plt.ylabel("Y Strip")
    plt.title("2D Histogram Heatmap")
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

def mirrorArrayRight(Array):
    half = len(Array) // 2
    return np.concatenate( [Array[-1:half-1:-1], Array[half:]] )

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def gaussMu(x,A,sigma):
    return A * np.exp(-(x - fixmu)**2 / (2*sigma**2))

def fit():
    yA = np.array([0.00904, 0.01419, 0.02185, 0.02594])
    xA = np.array([0.0034, 0.00497, 0.00744, 0.00857])

    yC = np.array([0.00866, 0.001322, 0.02223, 0.02371])
    xC = np.array([0.0031, 0.0047, 0.0069, 0.0085])
    
    p0 = [3,0]

    popt, pcov = curve_fit(linear, xA, yA, p0 = p0)
    perr = np.sqrt(np.diag(pcov))

    
    plt.plot(xA,linear(xA, *popt))
    #plt.title("Aluminium")
    plt.xlabel("theoretical value")
    plt.ylabel("experimental value")

    plt.text(xA[1], np.max(yA) - np.max(yA)/10, rf"$m={popt[0]:.2f} + {perr[0]:.2f}$", fontsize=10)
    plt.savefig("finished_plots/Aluminium.png")
    plt.show()

    popt, pcov = curve_fit(linear, xC, yC, p0=p0)
    perr = np.sqrt(np.diag(pcov))

    plt.plot(xC, linear(xC, *popt))
    #plt.title("Copper")
    plt.xlabel("theoretical value")
    plt.ylabel("experimental value")
    plt.text(xC[1], np.max(yC) - np.max(yC)/10, rf"$m={popt[0]:.2f} + {perr[0]:.2f}$", fontsize=10)
    plt.savefig("finished_plots/Copper.png")

    plt.show()



def linear(x, m, b):
    return m*x + b



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
#fit()



