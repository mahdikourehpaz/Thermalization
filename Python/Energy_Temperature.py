address = "/home/mahdi/codes/F_H_Myself/Systematic_Data/HWBC/Test_Included/Data/Data_2/"
address+= "15_7/diff_v1/Vr=0.01/ETH/"
folder1 = glob.glob(address+"*13.0/")

for folder in folder1:
    if os.path.isfile(folder+"eigenenergies.dat") == False:
        continue
        print("skipped",folder)
    if os.path.isfile(folder+"Temperature.dat") == False:
        continue
        print("skipped",folder)
    if os.path.isfile(folder+"Energy_Temperature.dat") == False:
        continue
        print("skipped",folder)

    Temp_sub = np.loadtxt(folder+"Temperature.dat")
    E = np.loadtxt(folder+"eigenenergies.dat")
    E_in = np.zeros(Ensemble_Size)
    NNE = E.size


    for i_it in range(Ensemble_Size):
        WF_in = np.loadtxt(folder+"Wave_Function_Initial_"+str(i_it)+".dat")
        E_in[i_it] = np.sum( WF_in**2*E)
    TT = np.zeros((Ensemble_Size,2))
    TT[:,1] = Temp_sub
    TT[:,0] = E_in
    np.savetxt(folder+"Energy_Temperature.dat",TT)
