import numpy as np


gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;



mu = 4 * np.pi * 10**(-7) #magnetic constant of vacuum permeability
h_planck = 1.055 * 10**(-34); #reduced Planck constant
rN = 0.101 * 10**(-9); # average cubic length of N-H bond
d = 1 * (mu * gammaN * gammaH * h_planck) / (4 * np.pi * rN**3); # dipolar coupling constant


# Define the Redfield R2 equations
def R2_function_max(w, x, magnetic_field):
    """-R2 Redfield equation for the maximum search"""
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    return  -2*np.sum(0.5 * (d**2 / 20) * (4 * (w *  x) 
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x))
                                + 6 * (w *  x / (1 + wh * wh * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 90 * (4 * (w  * x ) 
                                       + 3 *  (w  * x / (1 + (wn) *(wn) * x * x))))


def R2_function_min(w, x, magnetic_field):
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    return  2*np.sum(0.5 * (d**2 / 20) * (4 * (w *  x) 
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x))
                                + 6 * (w *  x / (1 + wh * wh * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 90 * (4 * (w *  x ) 
                                       + 3 *  (w  * x / (1 + (wn) *(wn) * x * x))))




def R1_function_max(w, x, magnetic_field):
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    return  -2*np.sum( (d**2 / 20) * (
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 15 * (   (w  * x / (1 + (wn) *(wn) * x * x))))


def R1_function_min(w, x, magnetic_field):
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    return  2*np.sum( (d**2 / 20) * (
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 15 * (   (w  * x / (1 + (wn) *(wn) * x * x))))



def NOE_function_min(w, x, magnetic_field):
    """Does not seem to work for more than 2 timescales"""
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    R1= 2*np.sum( (d**2 / 20) * (
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 15 * (   (w  * x / (1 + (wn) *(wn) * x * x))))
    
    return 1 +  2 * np.sum((d**2 / 20) *( 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))
                                    -  1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) ) * gammaH) / (gammaN * R1)

    NOE = 1 + (d**2 / 20) * (6 * JhPn - 1 * JhMn) * gammaH / (gammaN * R1);

    
def NOE_function_max(w, x, magnetic_field):
    """Does not seem to work for more than 2 timescales"""
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    R1= 2*np.sum( (d**2 / 20) * (
                                + 3 * (w *  x / (1 + wn * wn * x * x)) 
                                + 1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) 
                                + 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))) 
    + (wn * 160 * 10**(-6))**2 / 15 * (   (w  * x / (1 + (wn) *(wn) * x * x))))
    
    return -(1 + (d**2 / 20) * 2 * np.sum( 6 * (w *  x / (1 + (wh+wn) *(wh+wn) * x * x))
                                    -  1 * (w *  x / (1 + (wh-wn) *(wh-wn) * x * x)) ) * gammaH / (gammaN * R1))
                                    
                                    
                                    
                                    
# Define the constraint that the sum of w_i is 1
def constraint(w):
    return np.sum(w) - 1

# Define the constraint that the sum of w_i * x_i is equal to b
def constraint_sum(w, x, b):
    return np.sum(w * x) - b


# Define the constraint that w_i must be greater than 0
def constraint_positive(w):
    return w


# Function to get effective time from R2 values, using the Python np.poly1d function
def find_tau_from_R2_np_poly1d(magn_field_MHz,R2):
    """Function takes magnetic field, R2 value, returns Effective time
    
    Inputs:
         magn_field_MHz - magnetic field [MHz]
         R2             - spin relaxation R2 value [s-1]
         
    Output:
         eff_time       - Effective time calculated from
                          Redfield equation [s]
    """


    gammaH=267.513*10**6;  #r*s^(-1)*T^(-1)
    gammaN=-27.166*10**6;

    magnetic_field=magn_field_MHz*2*np.pi/gammaH*10**6

    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 

    mu = 4 * np.pi * 10**(-7) #magnetic constant of vacuum permeability
    h_planck = 1.055 * 10**(-34); #reduced Planck constant
    rN = 0.101 * 10**(-9); # average cubic length of N-H bond
    d = 1 * (mu * gammaN * gammaH * h_planck) / (4 * np.pi * rN**3); # dipolar coupling constant

    K1=d*d/20
    K2=(wn * 160 * 10**(-6))**2 / 45

    A=(wh-wn)**2
    B=wn**2
    C=wh**2
    D=(wh+wn)**2


    # coefficients calculated from Redfield R2 equation
    a=(4*K1+4*K2)*A*B*C*D
    b=-R2*A*B*C*D
    c=((10*K1+4*K2)*A*B*C+(10*K1+4*K2)*A*B*D+(7*K1+7*K2)*A*C*D+(5*K1+4*K2)*B*C*D)
    d=-R2*(A*B*C+A*B*D+A*C*D+B*C*D)
    e=(16*K1+4*K2)*A*B+(13*K1+7*K2)*A*C+(13*K1+7*K2)*A*D+(11*K1+4*K2)*B*C+(11*K1+4*K2)*B*D+(8*K1+7*K2)*C*D
    f=-R2*(A*B+A*C+A*D+B*C+B*D+C*D)
    g=(19*K1+7*K2)*A+(17*K1+4*K2)*B+(14*K1+7*K2)*C+(14*K1+7*K2)*D
    h=-R2*(A+B+C+D)
    j=20*K1+7*K2
    k=-R2
    
    



    p=np.poly1d([a,b,c,d,e,f,g,h,j,k])
    a=p.r
    eff_time=[]
    for i in a:
        if i.imag==0:
            eff_time.append(i.real)
    # this part is probably not needed,
    # the function should be monotonic and only 
    # one real solution should exist
    if len(eff_time)==1:
        eff_time=eff_time[0]
    
    return eff_time
    
    
    
    
def find_tau_from_R2_local_search(T2_threshold,magn_field_MHz,R2):
    
    #get magnetic field in Tesla
    gammaH=267.513*10**6;  #r*s^(-1)*T^(-1)
    magnetic_field=magn_field_MHz*2*np.pi/gammaH*10**6 # [T]
    
    # initiate the search
    T1s, T2s, NOEs, x = getSRT(-20,-3,10,magnetic_field)
    
    T2=1/R2
    
    minT2=list(np.absolute(T2s-T2))
    T2assign=T2s[minT2.index(min(np.absolute(minT2)))]

                        
    while np.absolute(T2assign-T2)>T2_threshold:
        if T2assign>T2:
            fast=np.log10(x[minT2.index(min(np.absolute(minT2)))])
            slow=np.log10(x[minT2.index(min(np.absolute(minT2)))+1])
            T1s, T2s, NOEs, x = getSRT(fast,slow,10,magnetic_field)
        else:
            fast=np.log10(x[minT2.index(min(np.absolute(minT2)))-1])
            slow=np.log10(x[minT2.index(min(np.absolute(minT2)))])
            T1s, T2s, NOEs, x = getSRT(fast,slow,10,magnetic_field)
            
        minT2=list(np.absolute(T2s-T2))
        T2assign=T2s[minT2.index(min(np.absolute(minT2)))]
        
    if T2assign>T2:
        T2sm=T2s[minT2.index(min(np.absolute(minT2)))+1]
        T2_borders_diff=T2assign-T2sm
        T2_diff_from_small=T2-T2sm
        factor=T2_diff_from_small/T2_borders_diff
        effT2=x[minT2.index(min(np.absolute(minT2)))+1]-factor*(x[minT2.index(min(np.absolute(minT2)))+1]-x[minT2.index(min(np.absolute(minT2)))])
    else:
        T2b=T2s[minT2.index(min(np.absolute(minT2)))-1]
        T2_borders_diff=T2b-T2assign
        T2_diff_from_small=T2assign-T2
        factor=T2_diff_from_small/T2_borders_diff
        effT2=x[minT2.index(min(np.absolute(minT2)))]-factor*(x[minT2.index(min(np.absolute(minT2)))]-x[minT2.index(min(np.absolute(minT2)))-1])

        
    return effT2
    
def getSRT(fast,slow,times,magnetic_field):

    x=np.logspace((fast),(slow),times) # x-akselin rajat

    T1s=[]
    T2s=[]
    NOEs=[]

    for b in x:
        T1,T2,NOE=get_relaxation_N(magnetic_field,[1],[b],0)
        T1s.append(T1)
        T2s.append(T2)
        NOEs.append(NOE)

    T1s=np.array(T1s)
    T2s=np.array(T2s)
    NOEs=np.array(NOEs)
    
    return T1s, T2s, NOEs, x
    
    
    
    
def get_relaxation_N(magnetic_field,Coeffs,Ctimes,OP):
    
    
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    
    #initiate spectral densities
    J0 = 0
    JhMn = 0
    JhPn = 0
    Jh = 0
    Jn = 0

    m = len(Ctimes)
    for i in range(0, m):
        w = 0
      
        J0 = J0 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])
        
        w = wh-wn;
        JhMn = JhMn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wn;
        Jn = Jn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])
        
        w = wh;
        Jh= Jh + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wn+wh;
        JhPn = JhPn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])


    mu = 4 * np.pi * 10**(-7) #magnetic constant of vacuum permeability
    h_planck = 1.055 * 10**(-34); #reduced Planck constant
    rN = 0.101 * 10**(-9); # average cubic length of N-H bond
    d = 1 * (mu * gammaN * gammaH * h_planck) / (4 * np.pi * rN**3); # dipolar coupling constant

    #units were corrected by S.Ollila and E.Mantzari, removed 2*pi from R1 and R2
    R1 = (d**2 / 20) * (1 * JhMn + 3 * Jn + 6 * JhPn) + Jn * (wn * 160 * 10**(-6))**2 / 15   ; 
    R2 = 0.5 * (d**2 / 20) * (4 * J0 + 3 * Jn + 1 * JhMn + 6 * Jh + 6 * JhPn) + (wn * 160 * 10**(-6))**2 / 90 * (4 * J0 + 3 * Jn);
    NOE = 1 + (d**2 / 20) * (6 * JhPn - 1 * JhMn) * gammaH / (gammaN * R1);


    #print("T1: {}, T2: {}, NOE: {}".format(1/R1, 1/R2, NOE))
    
    
           
    return 1/R1, 1/R2, NOE


def find_tau_from_R2_linear_approximation(magn_field_MHz,R2):
    magnetic_field=magn_field_MHz*2*np.pi/gammaH*10**6

    wn = gammaN * magnetic_field 

    
    K1=(d**2 / 20)
    K2=(wn * 160 * 10**(-6))**2 / 45

    eff_time=R2/(4*K1+4*K2)
    
    return eff_time
    
if __name__ == "__main__":
    magn_field_MHz=float(input("Magnetic field in MHz "))
    R2=float(input("R2 value in 1/s "))
    print(f'Effective time is {find_tau_from_R2_np_poly1d(magn_field_MHz,R2)*10**9:.3f} s')
