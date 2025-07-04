import ROOT
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# Open the ROOT file
#file = ROOT.TFile.Open("../data/elsa/half_radiation_length/Aluminum_5cm/40cm_distance/run87.root")
#file = ROOT.TFile.Open("../data/elsa/no_target/run82.root")
file = ROOT.TFile.Open("../data/elsa/half_radiation_length/aluminium_5cm/distance_15cm/run83.root")

# Function to read all branches from a TTree
def read_tree(tree):
    print(f"\nReading tree: {tree.GetName()}")
    for entry in tree:
        for branch in tree.GetListOfBranches():
            name = branch.GetName()
            value = getattr(entry, name, None)

# Read the 'raw' tree
if file.Get("raw"):
    raw_tree = file.Get("raw")
    read_tree(raw_tree)

# Read the 'data' tree
if file.Get("data"):
    data_tree = file.Get("data")
    read_tree(data_tree)

# Check entry counts match
n_raw = raw_tree.GetEntries()
n_data = data_tree.GetEntries()
n = min(n_raw, n_data)  # Safest if they differ

print(f"Looping over {n} events...")

# Loop over entries in both trees
for i in range(n):
    raw_tree.GetEntry(i)
    data_tree.GetEntry(i)

# Create 1D histogram
hist_xpos = ROOT.TH1F("hist_xpos", "X Position (only ID=0 or 1);x position [strip];counts", 257, 0, 257)
#hist_ypos = ROOT.TH1F("hist_ypos", "Y Position (only ID=2 or 3);y position [strip];counts", 257, 0, 257)
hist_ypos = ROOT.TH1F("hist_ypos", "Y Position (only ID=2 or 3);y position [strip];counts", 100, 0, 101)
hist_charge = ROOT.TH1F("hist_charge", "title", 100, 0, 2500)
hist_apv = ROOT.TH1F("hist_apv", "APV ID", 4, 0, 4)
hist_apv_ch = ROOT.TH1F("hist_apv_ch", "APV CH", 128, 0, 128)
hist_cluster_count = ROOT.TH1F("hist_cluster_count", "Cluster count", 25, 1, 26)
hist_cluster_size = ROOT.TH1F("hist_cluster_size", "Cluster size", 25, 1, 26)

# Create 2D histogram: average X vs average Y
#hist_xy = ROOT.TH2F("hist_xy", "Average X vs Average Y;x position [mm];y position [mm]", 100, 0, 100, 100, 0, 100)
hist_xy = ROOT.TH2F("hist_xy", "Average X vs Average Y;x position [strip];y position [strip]", 256, 0, 257, 256, 0, 257)
hist_xy_no_avg = ROOT.TH2F("hist_xy_no_avg", "NOT Average X vs Average Y;x position [mm];y position [mm]", 100, 0, 100, 100, 0, 100)


# Loop through events
for entry_raw in raw_tree:

    # number of hits per triggered event (abt every 25ns)
    num_hits = len(entry_raw.apv_id)

    x_pos = 0
    sum_x = 0
    num_hits_x = 0
    y_pos = 0
    sum_y = 0
    num_hits_y = 0

    cluster_count = 0
    cluster_size = 0

    # a cluster is defined as a number of strips that fire next to each other
    for i in range(num_hits - 1):
        if (entry_raw.mm_strip[i] + 1) == entry_raw.mm_strip[i + 1]: 
            cluster_size += 1 # strips per cluster
            #continue
        else:
            cluster_count += 1 # number of clusters
            hist_cluster_size.Fill(cluster_size)
            cluster_size = 0

    for i in range(num_hits - 1):


        #charge_max = entry_data.apv_qmax[i]
        charge = entry_raw.apv_q[i] # type: {1,2,3,,,,}
        # len = 26 for every entry
        # print(len(charge)) # type seems to be vector
        # print(max(charge))

        # this is true
        # print(charge_max == max(charge))

        hist_charge.Fill(max(charge))
        
        # id of the apv which registers event
        apv_id = entry_raw.apv_id[i]
        apv_ch = entry_raw.apv_ch[i]
        strip = entry_raw.mm_strip[i]

        hist_apv.Fill(apv_id)
        hist_apv_ch.Fill(apv_ch)

        # {{{ if
        # to filter noise. this value can be adjusted to get a desired gaussian profile
        # highest charge doesn't seem to be the events we are looking for. see hist_charge plot. low charge are the most events
        #if True:
        if (
                #((strip + 1) == entry_raw.mm_strip[i + 1] and (strip + 2) == entry_raw.mm_strip[i + 2])
                #or ((strip - 1) == entry_raw.mm_strip[i - 1] and (strip + 1) == entry_raw.mm_strip[i + 1])
                #or ((strip - 1) == entry_raw.mm_strip[i - 1] and (strip - 2) == entry_raw.mm_strip[i - 2])
                #and ((strip + 3) != entry_raw.mm_strip[i + 3])
                cluster_count < 4 # if there are more than x clusters in total dont use the event
                #and max(charge) > 80 
                #and max(charge) < 1000
                #True
                ):
        #if ((strip + 1) != entry_raw.mm_strip[i+1] and (strip + 2) != entry_raw.mm_strip[i+2] and max(charge) > 80 and max(charge) < 1000):
    
            # sort where a particle hit
            # x
            if apv_id in [0,1]:
                x_pos = strip * 100. / 256. # 10cm long and 128 channels, 256 strips
                hist_xpos.Fill(x_pos) # counts events that happen at a certain x position
                sum_x += x_pos  
                num_hits_x += i
                # no weight
                #sum_x += x_pos
                #num_hits_x += 1

            # y
            if apv_id in [2,3]:
                y_pos = ((strip + 128) % 255) * 100. / 256.
                #y_pos = strip #* 100.0 / 256
                hist_ypos.Fill(y_pos)
                sum_y += y_pos * i
                num_hits_y += i
                #sum_y += y_pos
                #num_hits_y += 1
            #}}}

    hist_cluster_count.Fill(cluster_count)

    # average over all events
    # TODO weigh hits 
    # weigh with charge seems to make not much of a difference
    if num_hits_x > 0 and num_hits_y > 0:
        # average position which is hit, i.e. hit at 2 and hit at 10 will be hit at 6 = (2+10)/2
        # TODO actually don't average this??
        average_x = sum_x / num_hits_x
        average_y = sum_y / num_hits_y
        hist_xy.Fill(average_x, average_y)

# Draw the histogram
# TODO TeX font (computer modern serif)
# TODO fit gaussian curve to x and y

# fits
#fit_x = ROOT.TF1("fit_x", "gaus", 0, 100)
#hist_xpos.Fit("fit_x")
fit_y = ROOT.TF1("fit_y", "gaus", 0, 256)
hist_ypos.Fit("fit_y")

#fct = hist_ypos.GetFunction("gaus")
#sigma = fct.GetParameter(2)
#print("sigma: ",sigma)

c = ROOT.TCanvas("c", "Canvas", 800, 600)
hist_xpos.Draw()
c.SaveAs("x_position_histogram.png")  # Optional: save plot

hist_ypos.Draw()
c.SaveAs("y_position_histogram.png")  # Optional: save plot

hist_charge.Draw()
c.SaveAs("hist_charge.png")

hist_apv.Draw()
c.SaveAs("hist_apv.png")

hist_apv_ch.Draw()
c.SaveAs("hist_apv_ch.png")

hist_cluster_count.Draw()
c.SaveAs("hist_cluster_count.png")

hist_cluster_size.Draw()
c.SaveAs("hist_cluster_size.png")

# 2D hit map
hist_xy.Draw("COLZ")
c.SaveAs("xy_hitmap.png")

#hist_xy_no_avg.Draw("COLZ")
#c.SaveAs("xy_hitmap_no_avg.png")

# save bin centers, contents, and errors to file
with open("y_histogram_data.txt", "w") as f:
    f.write("# BinCenter  BinContent  BinError\n")
    for i in range(1, hist_ypos.GetNbinsX() + 1):
        bin_center = hist_ypos.GetBinCenter(i)
        bin_content = hist_ypos.GetBinContent(i)
        bin_error = hist_ypos.GetBinError(i)
        f.write(f"{bin_center}  {bin_content}  {bin_error}\n")

file.Close()
