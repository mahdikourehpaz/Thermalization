address = "/home/mahdi/codes/F_H_Myself/Systematic_Data/HWBC/Test_Included/Data/Data_2/"
address+= "15_7/diff_v1/Vr=0.01/ETH/"
folder1 = glob.glob(address+"*/")

for folder in folder1:
    if os.path.isfile(folder+"eigenenergies.dat") == False:
        continue
    if os.path.isfile(folder+"Rho_Energy_Basis_10.dat") == True:
        continue
    print(folder)
    Ensemble_Size = (np.loadtxt(folder+"Temperature.dat")).size
    for jjj in range(Ensemble_Size):
        l = np.loadtxt(folder+"Occ_Nat_Orb_"+str(jjj)+".dat")
        t , N = l.shape
        time = range(t)
        N = N-1

        rr_real = np.loadtxt(folder+"Re_Natural_Orb_"+str(jjj)+".dat")
        rr_imag = np.loadtxt(folder+"Im_Natural_Orb_"+str(jjj)+".dat")
        phit = np.zeros((N,N,t),dtype=np.complex)
        for i in range(N):
            for j in range(N):
                for tt in time:
                    phit[i,j,tt] = complex(rr_real[i+tt*N,j],rr_imag[i+tt*N,j])

        N_sub = N
        H_sub = np.zeros((N_sub,N_sub))
        E_sub = np.zeros(N_sub)
        OnsiteP = np.loadtxt(folder+"Onsite_Random_Potential.dat")

        if os.path.exists(folder+"Mean_Field_0.dat") is True:
            MF = np.loadtxt(folder+"Mean_Field_"+str(jjj)+".dat")
        else:
            MF = np.zeros(System_Size)

        H_sub[0] = OnsiteP[0] + MF[0]
        for i in range(N_sub-1):
                H_sub[i+1,i+1] = OnsiteP[i+1] + MF[i+1]
                H_sub[i,i+1] = -1.000
                H_sub[i+1,i] = -1.000
        E_sub , H_sub = np.linalg.eigh(H_sub)

        rho = np.zeros((N,N,t),dtype=np.complex)
        asis = np.zeros((N,N,t),dtype=np.complex)
        Assiss = np.zeros((N,N,t),dtype=np.complex)
        for it in time:
            for i in range(N):
                for j in range(N):
                    Assiss[i,j,it] = np.dot(phit[:,i,it],H_sub[:,j])
        for it in time:
            for j in range(N):
                for k in range(N):
                    for i in range(N):
                        asis[j,k,it] = asis[j,k,it]+l[it,i+1]*Assiss[i,j,it]*np.conj(Assiss[i,k,it])
        for n in range(N):
            for m in range(N):
                for j in range(N):
                    for k in range(N):
                        rho[n,m,:] = rho[n,m,:]+asis[k,j,:]*H_sub[n,j]*H_sub[m,k]

        np.save(folder+"Rho_Energy_Basis_"+str(jjj)+".dat",rho)

