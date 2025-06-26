import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile.Open("../data/run87.root")

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
h_xpos = ROOT.TH1F("h_xpos", "X Position (only ID=0 or 1);x position [mm];counts", 100, 0, 100)
h_ypos = ROOT.TH1F("h_ypos", "Y Position (only ID=2 or 3);y position [mm];counts", 100, 0, 100)

# Create 2D histogram: average X vs average Y
h_xy = ROOT.TH2F("h_xy", "Average X vs Average Y;x position [mm];y position [mm]", 100, 0, 100, 100, 0, 100)


# Loop through events
for entry in raw_tree:

    # number of hits per triggered event (abt every 25ns)
    n_hits = len(entry.apv_id)

    sum_x = 0
    nhits_x = 0
    sum_y = 0
    nhits_y = 0

    # sort where a particle hit
    for i in range(n_hits):

        charge = entry.apv_q[i]

        # id of the apv which registers event
        apv_id = entry.apv_id[i]
        strip = entry.mm_strip[i]

        # x
        if apv_id in [0,1]:
            x_pos = strip * 100. / 256.
            h_xpos.Fill(x_pos)
            sum_x += x_pos
            nhits_x += 1
        # y
        if apv_id in [2,3]:
            y_pos = ((strip + 128) % 256) * 100. / 256.
            #y_pos = strip * 100.0 / 256
            h_ypos.Fill(y_pos)
            sum_y += y_pos
            nhits_y += 1

    # average over all events
    if nhits_x > 0 and nhits_y > 0:
        average_x = sum_x / nhits_x
        average_y = sum_y / nhits_y
        h_xy.Fill(average_x, average_y)

# Draw the histogram
c = ROOT.TCanvas("c", "Canvas", 800, 600)
h_xpos.Draw()
c.SaveAs("x_position_histogram.png")  # Optional: save plot

h_ypos.Draw()
c.SaveAs("y_position_histogram.png")  # Optional: save plot

# 2D hit map
h_xy.Draw("COLZ")
c.SaveAs("xy_hitmap.png")

file.Close()
