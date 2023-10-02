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


                        experiments[ID][field]={}
                        experiments[ID][field]["results"]={}
                        experiments[ID][field]["results"]["AA"]=[]
                        experiments[ID][field]["results"]["atomID"]=[]
                        experiments[ID][field]["results"]["R2"]=[]

                    if "_Heteronucl_T2_list.T2_val_units" in line:
                        units=line.split()[1]
                        experiments[ID][field]["units"]=units
                    if "_T2.Heteronucl_T2_list_ID" in line:
                        counter_relaxations=-3
                    if counter_relaxations==-1:
                        read_relaxations=True
                    if read_relaxations:
                        if len(line.split())==21:
                            experiments[ID][field]["results"]["AA"].append(line.split()[6])
                            experiments[ID][field]["results"]["atomID"].append(int(line.split()[5]))
                            if experiments[ID][field]["units"]=="s-1" or experiments[ID][field]["units"]=="Hz":
                                try:
                                    experiments[ID][field]["results"]["R2"].append(float(line.split()[10]))
                                except:
                                    experiments[ID][field]["results"]["R2"].append(None)
                            elif experiments[ID][field]["units"]=="s":
                                try:
                                    experiments[ID][field]["results"]["R2"].append(1/float(line.split()[10]))
                                except:
                                    experiments[ID][field]["results"]["R2"].append(None)
                            elif experiments[ID][field]["units"]=="ms":
                                try:
                                    experiments[ID][field]["results"]["R2"].append(1/float(line.split()[10])*1000)
                                except:
                                    experiments[ID][field]["results"]["R2"].append(None)
                            elif experiments[ID][field]["units"]=="ms-1":
                                try:
                                    experiments[ID][field]["results"]["R2"].append(float(line.split()[10])/1000)
                                except:
                                    experiments[ID][field]["results"]["R2"].append(None)

                        else:
                            read_relaxations=False
                    if "_Atom_chem_shift.Assigned_chem_shift_list_ID" in line:
                        counter_shifts=-3
                    if counter_shifts==-1:
                        read_shifts=True
                    if read_shifts:

                        if len(line.split())==24:
                            if len(experiments[ID]["shifts"]["atomID_shift"])>0:
                                if int(line.split()[5])==experiments[ID]["shifts"]["atomID_shift"][-1]:
                                    for shift in shifts:
                                        if shift==line.split()[7]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[10]))
                                else:
                                    experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                    experiments[ID]["shifts"]["AA_shift"].append(line.split()[6])
                                    for shift in shifts:
                                        if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"])-1:
                                            experiments[ID]["shifts"][shift].append(None)

                                        if shift==line.split()[7]:
                                            experiments[ID]["shifts"][shift].append(float(line.split()[10]))



                            else:
                                experiments[ID]["shifts"]["atomID_shift"].append(int(line.split()[5]))
                                experiments[ID]["shifts"]["AA_shift"].append(line.split()[6])
                                for index,shift in enumerate(shifts):
                                    if shift==line.split()[7]:
                                        experiments[ID]["shifts"][shift].append(float(line.split()[10]))
                        else:
                            read_shifts=False
                            for shift in shifts:
                                if len(experiments[ID]["shifts"][shift])<len(experiments[ID]["shifts"]["atomID_shift"]):
                                    experiments[ID]["shifts"][shift].append(None)




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
