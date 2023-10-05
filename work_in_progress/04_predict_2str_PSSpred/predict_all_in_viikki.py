import sys
sys.path.append("../")
sys.path.append("../03_experimental_data_chemical_shifts/")
import numpy as np
import handle_BMBR as BMBR
import SRT_optimization as srt
import gc
import os
import time
gc.collect()

analysis_path="data/"

for folder in os.listdir(analysis_path):
    print(folder)
    start_time = time.time()
    os.system("/home/nenciric/Programmes/PSSpred_v4/PSSpred.pl "+analysis_path+folder+"/"+folder+".fasta")
    for file in os.listdir("."):
        try:
            if file not in ["try",".ipynb_checkpoints","04_predict_2str_PSSpred.ipynb","data","times.out","blast.out","predict_all_in_viikki.py"]:
                os.system("mv "+file+" "+analysis_path+folder)
            elif file=="blast.out":
                os.system("rm "+file)
        except:
            pass
    end_time = time.time()
    with open("times.out","a") as f:
        f.write(f"{folder}: {end_time-start_time:.2f} s \n")


