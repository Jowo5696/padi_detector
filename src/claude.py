import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile.Open("../data/run87.root")

def fit_gaussian_core(histogram, fit_range_sigma=1.5):
    """
    Fit a Gaussian to the core of a histogram distribution
    
    Args:
        histogram: ROOT TH1 histogram
        fit_range_sigma: How many sigma around the peak to include in fit
    
    Returns:
        ROOT TF1 fit function or None if fit failed
    """
    if histogram.GetEntries() == 0:
        print("Warning: Empty histogram, cannot fit")
        return None
    
    # Find the peak bin
    peak_bin = histogram.GetMaximumBin()
    peak_value = histogram.GetBinCenter(peak_bin)
    peak_height = histogram.GetBinContent(peak_bin)
    
    # Estimate initial sigma from histogram RMS
    rms = histogram.GetRMS()
    if rms <= 0:
        print("Warning: Invalid RMS, using default sigma estimate")
        rms = (histogram.GetXaxis().GetXmax() - histogram.GetXaxis().GetXmin()) / 10
    
    # Define fit range around the peak
    fit_min = peak_value - fit_range_sigma * rms
    fit_max = peak_value + fit_range_sigma * rms
    
    # Ensure fit range is within histogram bounds
    hist_min = histogram.GetXaxis().GetXmin()
    hist_max = histogram.GetXaxis().GetXmax()
    fit_min = max(fit_min, hist_min)
    fit_max = min(fit_max, hist_max)
    
    print(f"  Peak at: {peak_value:.3f}")
    print(f"  Initial RMS: {rms:.3f}")
    print(f"  Fit range: [{fit_min:.3f}, {fit_max:.3f}]")
    
    # Create and configure the fit function
    fit_func = ROOT.TF1("gauss_core", "gaus", fit_min, fit_max)
    fit_func.SetParameters(peak_height, peak_value, rms)
    fit_func.SetParNames("Amplitude", "Mean", "Sigma")
    
    # Set reasonable parameter limits
    fit_func.SetParLimits(0, 0, peak_height * 2)  # Amplitude
    fit_func.SetParLimits(1, fit_min, fit_max)     # Mean
    fit_func.SetParLimits(2, rms * 0.1, rms * 3)  # Sigma
    
    # Perform the fit
    fit_result = histogram.Fit(fit_func, "QRS")  # Q=quiet, R=range, S=return fit result
    
    if fit_result and fit_result.Status() == 0:
        print(f"  Fit successful!")
        print(f"    Mean: {fit_func.GetParameter(1):.3f} ± {fit_func.GetParError(1):.3f}")
        print(f"    Sigma: {fit_func.GetParameter(2):.3f} ± {fit_func.GetParError(2):.3f}")
        print(f"    Chi2/NDF: {fit_func.GetChisquare()/fit_func.GetNDF():.2f}")
        
        # Set line style for visibility
        fit_func.SetLineColor(ROOT.kRed)
        fit_func.SetLineWidth(2)
        
        return fit_func
    else:
        print(f"  Fit failed with status: {fit_result.Status() if fit_result else 'None'}")
        return None

def explore_data_structure(tree, max_entries=5):
    """Explore the structure of the ROOT tree to understand data format"""
    print(f"\n=== Data Structure Analysis ===")
    print(f"Tree: {tree.GetName()}")
    print(f"Total entries: {tree.GetEntries()}")
    
    # List all branches
    print("\nAvailable branches:")
    branches = {}
    for branch in tree.GetListOfBranches():
        branch_name = branch.GetName()
        branches[branch_name] = branch
        print(f"  {branch_name}")
    
    # Examine first few entries
    print(f"\nExamining first {max_entries} entries:")
    for entry_idx in range(min(max_entries, tree.GetEntries())):
        tree.GetEntry(entry_idx)
        print(f"\nEntry {entry_idx}:")
        
        # Check apv_id structure
        if hasattr(tree, 'apv_id'):
            apv_ids = tree.apv_id
            print(f"  apv_id: type={type(apv_ids)}, len={len(apv_ids) if hasattr(apv_ids, '__len__') else 'N/A'}")
            if hasattr(apv_ids, '__len__') and len(apv_ids) > 0:
                print(f"    first few values: {list(apv_ids[:min(5, len(apv_ids))])}")
        
        # Check mm_strip structure
        if hasattr(tree, 'mm_strip'):
            strips = tree.mm_strip
            print(f"  mm_strip: type={type(strips)}, len={len(strips) if hasattr(strips, '__len__') else 'N/A'}")
            if hasattr(strips, '__len__') and len(strips) > 0:
                print(f"    first few values: {list(strips[:min(5, len(strips))])}")
        
        # Check all possible charge branches
        charge_branches = ['apv_q', 'apv_charge', 'charge', 'adc', 'amplitude']
        for charge_name in charge_branches:
            if hasattr(tree, charge_name):
                charge_data = getattr(tree, charge_name)
                print(f"  {charge_name}: type={type(charge_data)}, len={len(charge_data) if hasattr(charge_data, '__len__') else 'N/A'}")
                if hasattr(charge_data, '__len__') and len(charge_data) > 0:
                    try:
                        # Try to access individual elements
                        first_vals = []
                        for i in range(min(3, len(charge_data))):
                            val = charge_data[i]
                            first_vals.append(val)
                        print(f"    first few values: {first_vals}")
                        print(f"    element type: {type(charge_data[0])}")
                    except Exception as e:
                        print(f"    Error accessing elements: {e}")
                elif not hasattr(charge_data, '__len__'):
                    print(f"    scalar value: {charge_data}")
    
    return branches

def safe_get_charges(tree, max_sample=1000):
    """Safely extract charge values for analysis"""
    charges = []
    charge_branch = None
    
    # Find the charge branch
    charge_candidates = ['apv_q', 'apv_charge', 'charge', 'adc', 'amplitude']
    for candidate in charge_candidates:
        if hasattr(tree, candidate):
            charge_branch = candidate
            print(f"Using charge branch: {charge_branch}")
            break
    
    if charge_branch is None:
        print("No charge branch found - will use uniform weighting")
        return None, None
    
    print(f"Sampling {max_sample} entries to analyze charge distribution...")
    
    entry_count = 0
    for entry in tree:
        if entry_count >= max_sample:
            break
        entry_count += 1
        
        if not hasattr(entry, 'apv_id'):
            continue
            
        n_hits = len(entry.apv_id)
        charge_data = getattr(entry, charge_branch)
        
        # Handle different data structures
        try:
            if hasattr(charge_data, '__len__'):
                # Array-like data
                for i in range(min(n_hits, len(charge_data))):
                    val = charge_data[i]
                    if hasattr(val, '__len__'):
                        # If individual elements are arrays, take first element or sum
                        if len(val) > 0:
                            charges.append(float(val[0]))  # or sum(val) for total charge
                    else:
                        charges.append(float(val))
            else:
                # Scalar data
                charges.append(float(charge_data))
        except Exception as e:
            print(f"Warning: Error processing entry {entry_count}: {e}")
            continue
    
    if len(charges) == 0:
        print("No valid charge values found")
        return None, None
    
    charges = np.array(charges)
    print(f"Collected {len(charges)} charge values")
    print(f"Charge range: {np.min(charges):.2f} to {np.max(charges):.2f}")
    print(f"Charge mean: {np.mean(charges):.2f}, std: {np.std(charges):.2f}")
    
    return charges, charge_branch

def determine_threshold(charges, method='percentile'):
    """Determine noise threshold from charge distribution"""
    if charges is None or len(charges) == 0:
        return 0.0
    
    # Filter out negative charges first
    positive_charges = charges[charges > 0]
    if len(positive_charges) == 0:
        print("Warning: No positive charges found!")
        return 0.0
    
    print(f"Charge statistics (positive only):")
    print(f"  Total charges: {len(charges)}, Positive: {len(positive_charges)}")
    print(f"  Negative fraction: {(len(charges) - len(positive_charges))/len(charges):.1%}")
    
    if method == 'percentile':
        threshold = np.percentile(positive_charges, 75)  # Keep top 25%
    elif method == 'sigma':
        threshold = np.mean(positive_charges) + 1.0 * np.std(positive_charges)
    else:
        threshold = np.median(positive_charges)
    
    print(f"Noise threshold ({method}, positive charges only): {threshold:.2f}")
    return threshold

def process_event_safely(entry, charge_branch=None, threshold=0.0):
    """Safely process a single event, handling various data structures"""
    if not hasattr(entry, 'apv_id') or not hasattr(entry, 'mm_strip'):
        return None
    
    n_hits = len(entry.apv_id)
    if n_hits == 0:
        return None
    
    hits = {
        'apv_id': [],
        'mm_strip': [],
        'charge': []
    }
    
    # Get charge data if available
    charge_data = None
    if charge_branch and hasattr(entry, charge_branch):
        charge_data = getattr(entry, charge_branch)
    
    for i in range(n_hits):
        apv_id = entry.apv_id[i]
        strip = entry.mm_strip[i]
        
        # Get charge value - filter out negative charges
        charge = 1.0  # default
        if charge_data is not None:
            try:
                if hasattr(charge_data, '__len__') and i < len(charge_data):
                    val = charge_data[i]
                    if hasattr(val, '__len__'):
                        charge = float(val[0]) if len(val) > 0 else 1.0
                    else:
                        charge = float(val)
                elif not hasattr(charge_data, '__len__'):
                    charge = float(charge_data)
            except:
                charge = 1.0
        
        # Apply threshold and filter negative charges
        if charge > 0 and charge >= threshold:  # Positive charges only
            hits['apv_id'].append(apv_id)
            hits['mm_strip'].append(strip)
            hits['charge'].append(charge)
    
    return hits if len(hits['apv_id']) > 0 else None

def calculate_weighted_positions(hits):
    """Calculate charge-weighted positions"""
    sum_x, sum_y = 0.0, 0.0
    weight_x, weight_y = 0.0, 0.0
    
    for i in range(len(hits['apv_id'])):
        apv_id = hits['apv_id'][i]
        strip = hits['mm_strip'][i]
        charge = hits['charge'][i]
        
        if apv_id in [0, 1]:  # X strips
            x_pos = strip * 100.0 / 256.0
            sum_x += x_pos * charge
            weight_x += charge
        
        if apv_id in [2, 3]:  # Y strips
            y_pos = ((strip + 128) % 256) * 100.0 / 256.0
            sum_y += y_pos * charge
            weight_y += charge
    
    avg_x = sum_x / weight_x if weight_x > 0 else None
    avg_y = sum_y / weight_y if weight_y > 0 else None
    
    return avg_x, avg_y, weight_x, weight_y

# Main analysis
print("=== Micromegas Multiple Scattering Analysis ===")

# Get the tree
raw_tree = file.Get("raw")
if not raw_tree:
    print("Error: 'raw' tree not found!")
    exit()

# Explore data structure
branches = explore_data_structure(raw_tree)

# Analyze charge distribution
charges, charge_branch = safe_get_charges(raw_tree)
threshold = determine_threshold(charges) if charges is not None else 0.0

# Create histograms
h_xpos = ROOT.TH1F("h_xpos", "X Position;x position [mm];counts", 100, 0, 100)
h_ypos = ROOT.TH1F("h_ypos", "Y Position;y position [mm];counts", 100, 0, 100)
h_xy = ROOT.TH2F("h_xy", "Weighted Hit Map;x position [mm];y position [mm]", 100, 0, 100, 100, 0, 100)

if charges is not None:
    h_charge = ROOT.TH1F("h_charge", "Charge Distribution;charge;counts", 100, 
                        np.min(charges)*0.9, np.max(charges)*1.1)
    # Fill charge histogram
    for charge in charges:
        h_charge.Fill(charge)

# Process all events
print(f"\n=== Processing Events ===")
total_events = 0
valid_events = 0
xy_events = 0

for entry in raw_tree:
    total_events += 1
    if total_events % 1000 == 0:
        print(f"Processed {total_events} events...")
    
    hits = process_event_safely(entry, charge_branch, threshold)
    if hits is None:
        continue
    
    valid_events += 1
    
    # Fill individual hit positions
    for i in range(len(hits['apv_id'])):
        apv_id = hits['apv_id'][i]
        strip = hits['mm_strip'][i]
        
        if apv_id in [0, 1]:
            x_pos = strip * 100.0 / 256.0
            h_xpos.Fill(x_pos)
        
        if apv_id in [2, 3]:
            y_pos = ((strip + 128) % 256) * 100.0 / 256.0
            h_ypos.Fill(y_pos)
    
    # Calculate weighted average position for this event
    avg_x, avg_y, weight_x, weight_y = calculate_weighted_positions(hits)
    
    if avg_x is not None and avg_y is not None:
        # Weight by geometric mean of x and y total charges
        event_weight = np.sqrt(weight_x * weight_y)
        h_xy.Fill(avg_x, avg_y, event_weight)
        xy_events += 1

print(f"\n=== Results ===")
print(f"Total events: {total_events}")
print(f"Valid events (after filtering): {valid_events}")
print(f"Events with both X and Y: {xy_events}")
if total_events > 0:
    print(f"Efficiency: {valid_events/total_events:.1%}")

# Create plots
c = ROOT.TCanvas("c", "Micromegas Analysis", 1200, 900)
c.Divide(2, 2)

# X position
c.cd(1)
h_xpos.Draw()
print("\nFitting X position (core):")
fit_x_rough = fit_gaussian_core(h_xpos, fit_range_sigma=1.2)
c.Update()

# Y position  
c.cd(2)
h_ypos.Draw()
print("\nFitting Y position (core):")
fit_y_rough = fit_gaussian_core(h_ypos, fit_range_sigma=1.2)
c.Update()

# 2D map
c.cd(3)
h_xy.Draw("COLZ")
c.Update()

# Charge distribution (if available)
c.cd(4)
if charges is not None:
    h_charge.Draw()
    ROOT.gPad.SetLogy()
else:
    # Draw a text message if no charge data
    text = ROOT.TText(0.5, 0.5, "No charge data available")
    text.SetTextAlign(22)
    text.Draw()
c.Update()

c.SaveAs("micromegas_analysis.png")

# Individual high-quality plots with core Gaussian fits
c1 = ROOT.TCanvas("c1", "X Position", 800, 600)
h_xpos.Draw()
print("\nFitting X position (core only):")
fit_x = fit_gaussian_core(h_xpos, fit_range_sigma=1.2)  # Tighter core fit
c1.SaveAs("x_position.png")

c2 = ROOT.TCanvas("c2", "Y Position", 800, 600)
h_ypos.Draw()
print("\nFitting Y position (core only):")
fit_y = fit_gaussian_core(h_ypos, fit_range_sigma=1.2)
c2.SaveAs("y_position.png")

c3 = ROOT.TCanvas("c3", "2D Hit Map", 800, 600)
h_xy.Draw("COLZ")
c3.SaveAs("xy_hitmap.png")

# Print Gaussian fit results for the core fits
print(f"\n=== Multiple Scattering Analysis Results ===")
if fit_x and fit_x.Status() == 0:
    mean_x = fit_x.Parameter(1)
    sigma_x = fit_x.Parameter(2) 
    mean_err_x = fit_x.ParError(1)
    sigma_err_x = fit_x.ParError(2)
    chi2_x = fit_x.Chi2() / fit_x.Ndf() if fit_x.Ndf() > 0 else 0
    
    print(f"\nX Position Gaussian (core fit):")
    print(f"  Mean: {mean_x:.3f} ± {mean_err_x:.3f} mm")
    print(f"  Sigma: {sigma_x:.3f} ± {sigma_err_x:.3f} mm")
    print(f"  Chi2/NDF: {chi2_x:.2f}")
    print(f"  FWHM: {2.355 * sigma_x:.3f} mm")  # Full Width Half Maximum

if fit_y and fit_y.Status() == 0:
    mean_y = fit_y.Parameter(1)
    sigma_y = fit_y.Parameter(2)
    mean_err_y = fit_y.ParError(1)
    sigma_err_y = fit_y.ParError(2)
    chi2_y = fit_y.Chi2() / fit_y.Ndf() if fit_y.Ndf() > 0 else 0
    
    print(f"\nY Position Gaussian (core fit):")
    print(f"  Mean: {mean_y:.3f} ± {mean_err_y:.3f} mm")
    print(f"  Sigma: {sigma_y:.3f} ± {sigma_err_y:.3f} mm")
    print(f"  Chi2/NDF: {chi2_y:.2f}")
    print(f"  FWHM: {2.355 * sigma_y:.3f} mm")

# Calculate combined scattering angle (if both fits successful)
if (fit_x and fit_x.Status() == 0 and fit_y and fit_y.Status() == 0):
    # Assuming detector is at distance d from target (you'll need to specify this)
    # For now, let's assume a typical distance
    detector_distance = 400.0  # mm - adjust based on your setup
    
    # Convert position spread to angular spread (small angle approximation)
    theta_x_mrad = sigma_x / detector_distance * 1000  # mrad
    theta_y_mrad = sigma_y / detector_distance * 1000  # mrad
    theta_rms_mrad = np.sqrt(theta_x_mrad**2 + theta_y_mrad**2)
    
    print(f"\n=== Multiple Scattering Angles ===")
    print(f"(Assuming detector distance: {detector_distance} mm)")
    print(f"  θ_x (RMS): {theta_x_mrad:.2f} mrad")
    print(f"  θ_y (RMS): {theta_y_mrad:.2f} mrad") 
    print(f"  θ_total (RMS): {theta_rms_mrad:.2f} mrad")
    print(f"  θ_total (degrees): {theta_rms_mrad * 0.0573:.4f}°")  # mrad to degrees

file.Close()
print("\nAnalysis complete!")
