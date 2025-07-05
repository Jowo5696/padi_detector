import uproot as ur
import numpy as np

# naming convention: 
# t = tree
# b = branch

# with function closes root file automatically after use
with ur.open("../data/run37.root") as file:

    # print(file.keys())
    # ['raw;1', 'pedestals;1', 'data;1', 'run_info;1', 'config;1', 'config/TEnv;1']

    info = file["run_info"]

    config = file["config"]

    t_raw = file["raw"]

    t_data = file["data"]

    print(config.keys())
    print(config["TEnv"])
    print(info.keys())
    print(info["comment"])

    # print(t_raw.keys())
    # ['apv_evt', 'time_s', 'time_us', 'apv_fecNo', 'apv_id', 'apv_ch', 'mm_id', 'mm_readout', 'mm_strip', 'apv_q', 'apv_presamples']

    # print(t_data.keys())
    # ['apv_qmax', 'apv_tbqmax']

    # each event has a number for it. not all numbers from 0 to 281. maybe empty events are left out
    b_apv_evt = t_raw["apv_evt"].array(library="np")

    b_apv_ch = t_raw["apv_ch"].array(library="np")
    b_mm_strip = t_raw["mm_strip"].array(library="np")
    b_apv_q = t_raw["apv_q"].array(library="np")
    b_apv_qmax = t_data["apv_qmax"].array(library="np")

    # 225 seems to be the number of arrays in the file
    for i in range(225):
        filename = f"b_apv_qmax_and_mm_strip_{i}"
        with open(filename, "w") as f:
            for j in range(len(b_apv_qmax[i])):
                f.write("{0} {1}\n".format(b_mm_strip[i][j], b_apv_qmax[i][j]))

    
