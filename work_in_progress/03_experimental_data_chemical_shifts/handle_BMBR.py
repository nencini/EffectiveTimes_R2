import os




def load_BMBR(BMBR_path:str):
    # load in experimental data
    experiments={}
    for root, dirs, files in os.walk(BMBR_path):
        if len(files)==1:
            file=root+"/"+files[0]
            with open(file,"r") as f:
                read_relaxations=False
                read_shifts=False
                counter_relaxations=0 # saves data when -1
                counter_shifts=0 # saves data when -1
                for idline,line in enumerate(f):
                    counter_relaxations+=1
                    counter_shifts+=1
                    if "_Entry.ID" in line:
                        ID=line.split()[1]
                        experiments[ID]={}
                        experiments[ID]["disordered"]=False
                        experiments[ID]["micelle"]=False
                        experiments[ID]["T2measur"]={}
                    if '_Heteronucl_T2_list.Sf_framecode' in line:
                        mesurment=line.split()[1]
                    if '_Heteronucl_T2_list.Sample_condition_list_label' in line:
                        condition=line.split()[1]
                        if condition not in experiments[ID]["T2measur"]:
                            experiments[ID]["T2measur"][condition]={}
                        experiments[ID]["T2measur"][condition][mesurment]={}
                    if "_Heteronucl_T2_list.Spectrometer_frequency_1H" in line:
                        field=line.split()[1]
                        try:
                            if float(field)>1100:
                                field=float(field)/1000000
                            experiments[ID]["T2measur"][condition][mesurment]['field']=float(field)
                        except:
                            print(ID,field)


                        
                        experiments[ID]["T2measur"][condition][mesurment]["results"]={}
                        experiments[ID]["T2measur"][condition][mesurment]["results"]["AA"]=[]
                        experiments[ID]["T2measur"][condition][mesurment]["results"]["atomID"]=[]
                        experiments[ID]["T2measur"][condition][mesurment]["results"]["atomIDreal"]=[]
                        experiments[ID]["T2measur"][condition][mesurment]["results"]["R2"]=[]
                        experiments[ID]["T2measur"][condition][mesurment]["results"]["R2error"]=[]
                    
                    if "_Heteronucl_T2_list.T2_val_units" in line:
                        units=line.split()[1]
                        experiments[ID]["T2measur"][condition][mesurment]["units"]=units
                    if "_T2.Heteronucl_T2_list_ID" in line:
                        counter_relaxations=-3
                    if counter_relaxations==-1:
                        read_relaxations=True
                    if read_relaxations:
                        if len(line.split())==21:
                            experiments[ID]["T2measur"][condition][mesurment]['expType']=line.split()[7]
                            if experiments[ID]["T2measur"][condition][mesurment]["units"]=="s-1" or experiments[ID]["T2measur"][condition][mesurment]["units"]=="Hz":
                                try:
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2"].append(float(line.split()[10]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2error"].append(float(line.split()[11]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["AA"].append(line.split()[6])
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomID"].append(int(line.split()[5]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomIDreal"].append(int(line.split()[16]))
                            
                                except:
                                    pass
                            elif experiments[ID]["T2measur"][condition][mesurment]["units"]=="s":
                                try:
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2"].append(1/float(line.split()[10]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2error"].append(1/float(line.split()[11]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["AA"].append(line.split()[6])
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomID"].append(int(line.split()[5]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomIDreal"].append(int(line.split()[16]))
                            
                                except:
                                    pass
                            elif experiments[ID]["T2measur"][condition][mesurment]["units"]=="ms":
                                try:
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2"].append(1/float(line.split()[10])*1000)
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2error"].append(1/float(line.split()[11])*1000)
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["AA"].append(line.split()[6])
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomID"].append(int(line.split()[5]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomIDreal"].append(int(line.split()[16]))
                            
                                except:
                                    pass
                            elif experiments[ID]["T2measur"][condition][mesurment]["units"]=="ms-1":
                                try:
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2"].append(float(line.split()[10])/1000)
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["R2error"].append(float(line.split()[11])/1000)
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["AA"].append(line.split()[6])
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomID"].append(int(line.split()[5]))
                                    experiments[ID]["T2measur"][condition][mesurment]["results"]["atomIDreal"].append(int(line.split()[16]))
                            
                                except:
                                    pass
                        else:
                            read_relaxations=False
                    

                    if ("disorder" in line or "Disorder" in line):
                        experiments[ID]["disordered"]=True

                    if ("micelle" in line or "Micelle" in line):
                        experiments[ID]["micelle"]=True
                    if "_Assembly.Molecular_mass" in line:

                        if line.split()[1]!=".":
                            experiments[ID]["weight"]=float(line.split()[1])
                        else:
                            experiments[ID]["weight"]=None
    return experiments



#used to be load_BMBR function , includes chemical chifts
def load_BMBR_v2(BMBR_path:str):
    # load in experimental data
    experiments={}
    for root, dirs, files in os.walk(BMBR_path):
        if len(files)==1:
            file=root+"/"+files[0]
            with open(file,"r") as f:
                read_relaxations=False
                read_shifts=False
                counter_relaxations=0 # saves data when -1
                counter_shifts=0 # saves data when -1
                for idline,line in enumerate(f):
                    counter_relaxations+=1
                    counter_shifts+=1
                    if "_Entry.ID" in line:
                        ID=line.split()[1]
                        experiments[ID]={}
                        experiments[ID]["disordered"]=False
                        experiments[ID]["micelle"]=False
                        experiments[ID]["shifts"]={}
                        experiments[ID]["fields"]={}
                        experiments[ID]["shifts"]["AA_shift"]=[]
                        experiments[ID]["shifts"]["atomID_shift"]=[]
                        shifts=["H","HA","C","CA","CB","N"]
                        for shift in shifts:
                            experiments[ID]["shifts"][shift]=[]
                    if "_Heteronucl_T2_list.Spectrometer_frequency_1H" in line:
                        field=line.split()[1]
                        try:
                            if float(field)>1100:
                                field=float(field)/1000000
                        except:
                            print(ID,field)


                        experiments[ID]["fields"][field]={}
                        experiments[ID]["fields"][field]["results"]={}
                        experiments[ID]["fields"][field]["results"]["AA"]=[]
                        experiments[ID]["fields"][field]["results"]["atomID"]=[]
                        experiments[ID]["fields"][field]["results"]["atomIDreal"]=[]
                        experiments[ID]["fields"][field]["results"]["R2"]=[]
                        experiments[ID]["fields"][field]["results"]["R2error"]=[]

                    if "_Heteronucl_T2_list.T2_val_units" in line:
                        units=line.split()[1]
                        experiments[ID]["fields"][field]["units"]=units
                    if "_T2.Heteronucl_T2_list_ID" in line:
                        counter_relaxations=-3
                    if counter_relaxations==-1:
                        read_relaxations=True
                    if read_relaxations:
                        if len(line.split())==21:
                            experiments[ID]["fields"][field]["results"]["AA"].append(line.split()[6])
                            experiments[ID]["fields"][field]["results"]["atomID"].append(int(line.split()[5]))
                            experiments[ID]["fields"][field]["results"]["atomIDreal"].append(int(line.split()[16]))
                            if experiments[ID]["fields"][field]["units"]=="s-1" or experiments[ID]["fields"][field]["units"]=="Hz":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(float(line.split()[10]))
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(float(line.split()[11]))
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="s":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(1/float(line.split()[10]))
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(1/float(line.split()[11]))
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="ms":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(1/float(line.split()[10])*1000)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(1/float(line.split()[11])*1000)
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="ms-1":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(float(line.split()[10])/1000)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(float(line.split()[11])/1000)
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                                    experiments[ID]["fields"][field]["results"]["R2error"].append(None)

                        else:
                            read_relaxations=False
                    if "_Atom_chem_shift.Assigned_chem_shift_list_ID" in line:
                        counter_shifts=-3
                    if counter_shifts==-1:
                        read_shifts=True
                    if read_shifts:
                        if "stop_" in line or len(line.split())==0:
                            read_shifts=False
                            for shift in shifts:
                                if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"]):
                                    experiments[ID]["shifts"][shift].append(None)

                        else:
                            if len(experiments[ID]["shifts"]["atomID_shift"])>0:
                                if int(line.split()[5])==experiments[ID]["shifts"]["atomID_shift"][-1]:
                                    if line.split()[3]==".":
                                        atom_position=8
                                        shift_position=11
                                    elif line.split()[3]=="1":
                                        atom_position=7
                                        shift_position=10
                                    for shift in shifts:
                                        if shift==line.split()[atom_position]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))
                                else:
                                    
                                    if line.split()[3]==".":
                                        atom_position=8
                                        shift_position=11
                                        AA_name_position=7
                                    elif line.split()[3]=="1":
                                        atom_position=7
                                        shift_position=10
                                        AA_name_position=6
                                    
                                    experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                    experiments[ID]["shifts"]["AA_shift"].append(line.split()[AA_name_position])
                                    
                                    
                                    for shift in shifts:
                                        if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"])-1:
                                            experiments[ID]["shifts"][shift].append(None)

                                        if shift==line.split()[atom_position]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))



                            else:
                                if line.split()[3]==".":
                                    atom_position=8
                                    shift_position=11
                                    AA_name_position=7
                                elif line.split()[3]=="1":
                                    atom_position=7
                                    shift_position=10
                                    AA_name_position=6
                            
                            
                                experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                experiments[ID]["shifts"]["AA_shift"].append(line.split()[AA_name_position])
                                
                                for index,shift in enumerate(shifts):
                                    if shift==line.split()[atom_position]:
                                        experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))
                        


                    if ("disorder" in line or "Disorder" in line):
                        experiments[ID]["disordered"]=True

                    if ("micelle" in line or "Micelle" in line):
                        experiments[ID]["micelle"]=True
                    if "_Assembly.Molecular_mass" in line:

                        if line.split()[1]!=".":
                            experiments[ID]["weight"]=float(line.split()[1])
                        else:
                            experiments[ID]["weight"]=None
    return experiments
    
    



#used to be load_BMBR function 
def load_BMBR_old(BMBR_path:str):
    # load in experimental data
    experiments={}
    for root, dirs, files in os.walk(BMBR_path):
        if len(files)==1:
            file=root+"/"+files[0]
            with open(file,"r") as f:
                read_relaxations=False
                read_shifts=False
                counter_relaxations=0 # saves data when -1
                counter_shifts=0 # saves data when -1
                for idline,line in enumerate(f):
                    counter_relaxations+=1
                    counter_shifts+=1
                    if "_Entry.ID" in line:
                        ID=line.split()[1]
                        experiments[ID]={}
                        experiments[ID]["disordered"]=False
                        experiments[ID]["micelle"]=False
                        experiments[ID]["shifts"]={}
                        experiments[ID]["shifts"]["AA_shift"]=[]
                        experiments[ID]["shifts"]["atomID_shift"]=[]
                        shifts=["H","HA","C","CA","CB","N"]
                        for shift in shifts:
                            experiments[ID]["shifts"][shift]=[]
                    if "_Heteronucl_T2_list.Spectrometer_frequency_1H" in line:
                        field=line.split()[1]
                        try:
                            if float(field)>1100:
                                field=float(field)/1000000
                        except:
                            print(ID,field)


                        experiments[ID]["fields"][field]={}
                        experiments[ID]["fields"][field]["results"]={}
                        experiments[ID]["fields"][field]["results"]["AA"]=[]
                        experiments[ID]["fields"][field]["results"]["atomID"]=[]
                        experiments[ID]["fields"][field]["results"]["R2"]=[]

                    if "_Heteronucl_T2_list.T2_val_units" in line:
                        units=line.split()[1]
                        experiments[ID]["fields"][field]["units"]=units
                    if "_T2.Heteronucl_T2_list_ID" in line:
                        counter_relaxations=-3
                    if counter_relaxations==-1:
                        read_relaxations=True
                    if read_relaxations:
                        if len(line.split())==21:
                            experiments[ID]["fields"][field]["results"]["AA"].append(line.split()[6])
                            experiments[ID]["fields"][field]["results"]["atomID"].append(int(line.split()[5]))
                            if experiments[ID]["fields"][field]["units"]=="s-1" or experiments[ID]["fields"][field]["units"]=="Hz":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(float(line.split()[10]))
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="s":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(1/float(line.split()[10]))
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="ms":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(1/float(line.split()[10])*1000)
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)
                            elif experiments[ID]["fields"][field]["units"]=="ms-1":
                                try:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(float(line.split()[10])/1000)
                                except:
                                    experiments[ID]["fields"][field]["results"]["R2"].append(None)

                        else:
                            read_relaxations=False
                    if "_Atom_chem_shift.Assigned_chem_shift_list_ID" in line:
                        counter_shifts=-3
                    if counter_shifts==-1:
                        read_shifts=True
                    if read_shifts:
                        if "stop_" in line or len(line.split())==0:
                            read_shifts=False
                            for shift in shifts:
                                if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"]):
                                    experiments[ID]["shifts"][shift].append(None)

                        else:
                            if len(experiments[ID]["shifts"]["atomID_shift"])>0:
                                if int(line.split()[5])==experiments[ID]["shifts"]["atomID_shift"][-1]:
                                    if line.split()[3]==".":
                                        atom_position=8
                                        shift_position=11
                                    elif line.split()[3]=="1":
                                        atom_position=7
                                        shift_position=10
                                    for shift in shifts:
                                        if shift==line.split()[atom_position]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))
                                else:
                                    
                                    if line.split()[3]==".":
                                        atom_position=8
                                        shift_position=11
                                        AA_name_position=7
                                    elif line.split()[3]=="1":
                                        atom_position=7
                                        shift_position=10
                                        AA_name_position=6
                                    
                                    experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                    experiments[ID]["shifts"]["AA_shift"].append(line.split()[AA_name_position])
                                    
                                    
                                    for shift in shifts:
                                        if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"])-1:
                                            experiments[ID]["shifts"][shift].append(None)

                                        if shift==line.split()[atom_position]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))



                            else:
                                if line.split()[3]==".":
                                    atom_position=8
                                    shift_position=11
                                    AA_name_position=7
                                elif line.split()[3]=="1":
                                    atom_position=7
                                    shift_position=10
                                    AA_name_position=6
                            
                            
                                experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                experiments[ID]["shifts"]["AA_shift"].append(line.split()[AA_name_position])
                                
                                for index,shift in enumerate(shifts):
                                    if shift==line.split()[atom_position]:
                                        experiments[ID]["shifts"][shift].append(float(line.split()[shift_position]))
                        


                    if ("disorder" in line or "Disorder" in line):
                        experiments[ID]["disordered"]=True

                    if ("micelle" in line or "Micelle" in line):
                        experiments[ID]["micelle"]=True
                    if "_Assembly.Molecular_mass" in line:

                        if line.split()[1]!=".":
                            experiments[ID]["weight"]=float(line.split()[1])
                        else:
                            experiments[ID]["weight"]=None
    return experiments
    
    
    
def load_PSSpred_structure(PPSpred_path:str):
    secondary_structure={}
    for root, dirs, files in os.walk(PPSpred_path):
        for dire in dirs:
            secondary_structure[dire[3:]]={}
            secondary_structure[dire[3:]]["residues"]=[]
            secondary_structure[dire[3:]]["codes"]=[]
            secondary_structure[dire[3:]]["disordered"]=[]
            secondary_structure[dire[3:]]["helix"]=[]
            secondary_structure[dire[3:]]["sheet"]=[]
            with open(root+dire+"/seq.dat.ss","r") as f:
                for line in f:
                    if len(line.split())==6:
                        secondary_structure[dire[3:]]["residues"].append(int(line.split()[0]))
                        secondary_structure[dire[3:]]["disordered"].append(float(line.split()[3]))
                        secondary_structure[dire[3:]]["helix"].append(float(line.split()[4]))
                        secondary_structure[dire[3:]]["sheet"].append(float(line.split()[5]))
                        if line.split()[2]=="C":
                            secondary_structure[dire[3:]]["codes"].append(0)
                        elif line.split()[2]=="E":
                             secondary_structure[dire[3:]]["codes"].append(-1)
                        elif line.split()[2]=="H":
                             secondary_structure[dire[3:]]["codes"].append(1)
    return secondary_structure
    
def get_indexes(experiments:dict, ID:str,shifts:dict):
    """extracts chemical shift index and scales the values a constant"""
    index={}
    #shifts=["H","HA","N","C"]
    ["C","CA","CB","N"]
    for shift in shifts:
        index[shift]=[]
    for i,res in enumerate(experiments[ID]["shifts"]["AA_shift"]):
        #print(res)
        for shift in shifts:
            if experiments[ID]["shifts"][shift][i]!=None and res in chemical_shifts[shift] and chemical_shifts[shift][res]!=None:
                #print(shift,res)
                index[shift].append((experiments[ID]["shifts"][shift][i]-chemical_shifts[shift][res])*shifts[shift])
            else:
                index[shift].append(None)
    return index    
    
    
chemical_shifts={
    "H":{
"ALA":8.11,
"ARG":8.17,
"ASN":8.33,
"ASP":8.39,
"CYS":7.81,
"CYS2":8.53,
"GLN":8.25,
"GLU":8.29,
"GLY":8.34,
"HIS":8.09,
"ILE":7.94,
"LEU":8.12,
"LYS":8.13,
"MET":8.37,
"PHE":7.95,
"PRO":None,
"SER":8.26,
"THR":8.22,
"TRP":7.59,
"TYR":7.9,
"VAL":7.88
},
"HA":
{
"ALA":4.25,
"ARG":4.33,
"ASN":4.60,
"ASP":4.64,
"CYS":4.63,
"CYS2":4.44,
"GLN":4.26,
"GLU":4.28,
"GLY":3.95,
"HIS":4.50,
"ILE":4.13,
"LEU":4.35,
"LYS":4.28,
"MET":4.55,
"PHE":4.62,
"PRO":4.41,
"SER":4.48,
"THR":4.33,
"TRP":4.54,
"TYR":4.55,
"VAL":4.14
},
"N":
{
"ALA":132.52,
"ARG":120.59,
"ASN":118.48,
"ASP":120.69,
"CYS":117.01,
"CYS2":118.62,
"GLN":119.73,
"GLU":120.87,
"GLY":109.94,
"HIS":118.87,
"ILE":121.07,
"LEU":121.53,
"LYS":121.44,
"MET":120.19,
"PHE":119.41,
"PRO":None,
"SER":115.94,
"THR":114.41,
"TRP":120.57,
"TYR":120.05,
"VAL":119.66
},  
"C":
{
"ALA":177.39,
"ARG":175.91,
"ASN":174.98,
"ASP":176.01,
"CYS":174.77,
"CYS2":175.85,
"GLN":175.88,
"GLU":176.11,
"GLY":174.30,
"HIS":174.88,
"ILE":175.46,
"LEU":176.61,
"LYS":176.15 ,
"MET":175.93,
"PHE":175.28,
"PRO":176.91,
"SER":174.33,
"THR":174.62,
"TRP":175.91,
"TYR":175.32,
"VAL":175.76
},

"CA":
{
"ALA":52.67,
"ARG":55.96,
"ASN":52.94,
"ASP":54.09,
"CYS":58.8,
"CYS2":57.68,
"GLN":55.94,
"GLU":56.39,
"GLY":45.34,
"HIS":55.78,
"ILE":60.64,
"LEU":54.85,
"LYS":56.40 ,
"MET":55.12,
"PHE":56.94,
"PRO":63.53,
"SER":58.35,
"THR":61.59,
"TRP":57.62,
"TYR":57.72,
"VAL":61.80
},

"CB":
{
"ALA":19.03,
"ARG":30.53,
"ASN":38.22,
"ASP":40.76,
"CYS":29.75,
"CYS2":38.38,
"GLN":28.67,
"GLU":30.02,
"GLY":None,
"HIS":29.62,
"ILE":38.26,
"LEU":41.87,
"LYS":32.57 ,
"MET":32.93,
"PHE":39.43,
"PRO":31.87,
"SER":63.88,
"THR":69.75,
"TRP":29.27,
"TYR":38.71,
"VAL":32.68
}
  
    
}


