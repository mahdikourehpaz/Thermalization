address = "/home/mahdi/codes/F_H_Myself/Systematic_Data/HWBC/Test_Included/Data/Data_2/"
address+= "15_7/diff_v1/Vr=0.01/ETH/"
folder1 = glob.glob(address+"*13.0/")
first_col = 1


for folder in folder1:
    if os.path.isfile(folder+"eigenenergies.dat") == False:
        continue
        print("there is no eigenvalues file",folder)
    if os.path.isfile(folder+"Temperature.dat") == True:
        continue
        print(folder)
    E1 = np.loadtxt(folder+"eigenenergies.dat")

    inp = open(folder+"initial.in").read()
    initial = inp.split("\n")
    System_Size = int(initial[1].split()[0])
    inp = open(folder+"initial.in").read()
    initial = inp.split("\n")
    System_Size = int(initial[1].split()[0])
    find_text = [x for x in initial if 'which Site?' in x][0]
    index = initial.index(find_text)

    
    if initial[index + 1].split()[0] == "AL":
        Ensemble_Size = 1
    else:
        find_text = [x for x in initial if 'Ensembling over' in x][0]
        index = initial.index(find_text)
        Ensemble_Size = int(initial[index + 1])
        
    Temp = np.zeros((Ensemble_Size,2))
    f = open(folder+"Temperature.dat","w")

    if os.path.exists(folder+"Mean_Field_0.dat") is True:
        Mean_Field = 'yes'
        OnsiteP = np.loadtxt(folder+"Onsite_Random_Potential.dat")
    else:
        Mean_Field = 'No'
        OnsiteP = np.zeros(System_Size)

    for iii in range(Ensemble_Size):
        rr_real = np.loadtxt(folder+"Re_Natural_Orb_"+str(iii)+".dat")
        rr_imag = np.loadtxt(folder+"Im_Natural_Orb_"+str(iii)+".dat")
        l = np.loadtxt(folder+"Occ_Nat_Orb_"+str(iii)+".dat")
        t , N1 = l.shape
        N_sub = N1 - 1
        #########################################################################
        #########################################################################
        #########################################################################
        H_sub = np.zeros((N_sub,N_sub))
        if Mean_Field == 'yes':
            MF = np.loadtxt(folder+"Mean_Field_"+str(iii)+".dat")
            H_sub[0,0] = MF[0]
            for i in range(N_sub-1):
                H_sub[i+1,i+1] = MF[i+1]
        E_sub = np.zeros(N_sub)
        H_sub[0,0] = OnsiteP[0] + H_sub[0,0]
        for i in range(N_sub-1):
            H_sub[i+1,i+1] = OnsiteP[i+1] + H_sub[i+1,i+1]
            H_sub[i,i+1] = -1.000
            H_sub[i+1,i] = -1.000
        E_sub , H_sub = np.linalg.eigh(H_sub)
        #########################################################################
        #########################################################################
        #########################################################################

        N = N1-1
        l_av = np.zeros(N1-1)
        for i in range(t):
            for j in range(1,N1):
                l_av[j-1] = l_av[j-1] + l[i,j]
        l_av = l_av / t

        phit = np.zeros((N_sub,N_sub,t),dtype=np.complex)

        for i in range(N_sub):
            for j in range(N_sub):
                for tt in range(t):
                    phit[i,j,tt] = complex(rr_real[i+tt*N_sub,j],rr_imag[i+tt*N_sub,j])
        En_phi = np.zeros((t,N_sub))

        for it in range(t):
            for i in range(N_sub):
                for j in range(N_sub):
                     En_phi[it,i] = En_phi[it,i] + E_sub[j]*abs(np.dot(phit[:,i,it],H_sub[:,j]))**2
        
        En = np.zeros(N)
        for it in range(t):
            En[:] = En[:] + En_phi[it,:]
        En = En / t


        Temp_E = np.zeros((N,2))
        Temp_E[:,0] = En
        Temp_E[:,1] = l_av
        np.save(folder+"Gibbs_"+str(iii)+".dat",Temp_E)
        
        e = []
        l = []
        for i in range(g[:,0].size):
            if i > 3:
                l.append(np.log(l_av[i]))
                e.append(En[i])
        
        T = np.polyfit(e,l,1)
        f.write(str(-1.0/T)+"\n")









