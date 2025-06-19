import uproot as ur
import numpy as np

# with function closes root file automatically after use
with ur.open("../data/run37.root") as file:

    # print(file.keys())
    # ['raw;1', 'pedestals;1', 'data;1', 'run_info;1', 'config;1', 'config/TEnv;1']

    tree_raw = file["raw"]

    tree_data = file["data"]

    # print(tree_raw.keys())
    # ['apv_evt', 'time_s', 'time_us', 'apv_fecNo', 'apv_id', 'apv_ch', 'mm_id', 'mm_readout', 'mm_strip', 'apv_q', 'apv_presamples']

    # print(tree_data.keys())
    # ['apv_qmax', 'apv_tbqmax']

    # each event has a number for it
    branch_raw = tree_raw["apv_evt"].array(library="np")

    print(type(branch_raw))

    print(branch_raw)
