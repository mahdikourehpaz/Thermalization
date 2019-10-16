import scipy.optimize as spopt
import sciyp.special as sp.special
import os

#########################################################
#########################################################
#########################################################
#########################################################
pi = np.pi
def Exp(x):
    return np.exp(x)
def Nc_pois(x):
    n = 1.0-Exp(-x)
    return n
def Nc_wig(x):
    n = 1.0 - Exp(-pi*x*x/4.0)
    return n
def Nc_Brody(x,b):
    al = sp.special.gamma((b+2.0)/(b+1.0))**(b+1.0)
    n = 1.0 - Exp(-al*x**(b+1.0))
    return n
def dis_wig(x):
    ans = 0.5*pi*x*np.exp(-pi*x*x/4.0)
    return ans
def dis_brody(x,q):
    b = (sp.special.gamma((2.0+q)/(1.0+q)))**(q+1.0)
    ans = b*(1+q)*(x**q) *np.exp(-b*x**(q+1))
    return ans
def dis_BR(x,q):
    b = sp.special.erfc(np.sqrt(pi)*q*x/2.0)
    ans = Exp(x*(q-1))* ((1.-q)**2 * b + (2.0*q*(1.-q)+pi/2. * q**3 * x ) * Exp(-pi*q*q*x*x/4.) )
    return ans
def Nc_BR(x,q):
    b = sp.special.erfc(np.sqrt(pi)*q*x/2.0)
    ans = 1. - Exp(x*(q-1))* (( q * Exp(-pi*q*q*x*x/4.) + (1.-q)*b  ))
    return ans
#########################################################
#########################################################
#########################################################
#########################################################



address = "/home/mahdi/codes/F_H_Myself/Systematic_Data/HWBC/Test_Included/Data/Data_2/"
address+= "15_7/diff_v1/Vr=0.01/ETH/"
folder1 = glob.glob(address+"*/")

for folder in folder1:
    if os.path.isfile(folder+"eigenenergies.dat") == False:
        continue
        print(folder)
        print("No 'eigenenergies.dat' here")
    if os.path.isfile(folder+"brody_parameter.dat") == True:
        continue
        print(folder)
        print("Done Before!!!")

    E1 = np.loadtxt(folder+"eigenenergies.dat")
    itteration = np.size(E1.shape)
    percentage = 10
    polyorder = 10



    NN = E1.size
    r = int(percentage*1.0*NN/100.0)
    E = E1[r:NN-r]
    N = E.size

    
    
    Nc = np.zeros(10*N)
    rho , _ = np.histogram(E,bins=Nc.size-1)
    for i in range(Nc.size-1):
        Nc[i+1] = Nc[i] + rho[i]
    EE_x = np.linspace(E[0],E[-1],Nc.size)
    eta_av = np.polyfit(EE_x,Nc,polyorder)

    E_uf = np.polyval(eta_av,E)
    E_uf = np.sort(E_uf)
    np.savetxt(folder + "Unfolded_Spectrum.dat",E_uf)
    s = np.zeros(E.size)
    for i in range(N-1):
        s[i] = E_uf[i+1] - E_uf[i]

    s = s / np.average(s)

    N_s = s.size*10
    Nc_s = np.zeros(N_s)
    h , _ = np.histogram(s,bins=Nc_s.size-1)
    x = np.linspace(np.min(s),np.max(s),Nc_s.size)
    for i in range(Nc_s.size-1):
        Nc_s[i+1] = Nc_s[i] + h[i] 
    Nc = Nc_s / np.max(Nc_s)
    
    b , err_b = spopt.curve_fit(Nc_Brody,x,Nc)
    br , err_br = spopt.curve_fit(Nc_BR,x,Nc)
    
    print("    ", b , br)
    hist , temp = np.histogram(s,bins=100,density=True)
    dtemp = temp[1] - temp[0]

    edge = np.arange(0,hist.size)*dtemp + dtemp/2.

    save1_file = np.zeros((edge.size,2))
    save1_file[:,0] = edge
    save1_file[:,1] = hist
    np.savetxt(folder+"NNLS_Distribution.dat",save1_file)
    save2_file = np.zeros((Nc.size,2))
    save2_file[:,0] = x
    save2_file[:,1] = Nc
    np.savetxt(folder+"NNLS_Cumulation.dat",save2_file)
    np.savetxt(folder+"brody_parameter.dat",[b , err_b])
    np.savetxt(folder+"berry_robnik.dat",[br , err_br])
    del Nc_s
    del x
    del s
    del Nc
    del rho
    del save1_file
    del save2_file





