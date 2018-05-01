
import pandas as pds
import numpy as np
from scipy.interpolate import interp1d


__all__ = ['read_NRLMSIS','read_TIMED','Photo_Absorption','IonProducRate','Scaterring_Depth','Ionization_Rate']




def read_NRLMSIS(files, columns):
    df=pds.read_csv(files, skiprows=34, delim_whitespace=True,header=None,
                      index_col=0)
    df.index.name="Alt"
    df.columns=columns
    return df
    
    
def read_TIMED(files):
    df=pds.read_csv(files,index_col=1) # index= wavelength
    return df
    
    

#####################################################################
# Compute the Irradiance as a function of the optical depth
# if species=True, the output is the Intensity I for each species j
# if species=False, the output is the Intensity I of the total species
######################################################################
    
def Photo_Absorption(df,index,res,c_section, I_infinity,angle_deg, wavelength, j=None, species=False):

    Earth_Radii=6371 #km
    alt=index
    chsi=angle_deg*np.pi/180
    
    
    if species:
        
        I=[]
        
        for i in np.arange(len(alt)):
            z0=alt[i]
            z=np.array(alt[i:-1])
            dz=res*np.ones(len(z))# dz in cm
            ratio=(Earth_Radii+z0)/(Earth_Radii+z)
            density=df[df.columns[j]][z]
            function= density*1e6*(1-(ratio**2)*np.sin(chsi)**2)**(-0.5)
            Tau=c_section[j]*1e-4*np.trapz(function, x=z*1e3, dx=dz)
            I.append(I_infinity*np.exp(-Tau))
            i+=1
    else:
             
        Tau=np.zeros((len(c_section), len(alt)))
        
        print('============================')
        print('Wavelength=', str(wavelength) + ' Angstrom' )
        print('============================')
        
        for s in np.arange(len(c_section)):            

            print('species= ',df.columns[s] + ' // Absorption cross section=', str(c_section[s]) + ' cm^2' )
            #print('cross_section=', str(c_section[s]) + ' cm^2')

            for i in np.arange(len(alt)):
                z0=alt[i]
                z=np.array(alt[i:-1])
                dz=res*np.ones(len(z)) # dz in cm
                ratio=(Earth_Radii+z0)/(Earth_Radii+z)
        
                density=df[df.columns[s]][z]#m3
                function= 1e6*density*(1-(ratio**2)*np.sin(chsi)**2)**(-0.5)
                #function=density*1e6/(np.cos(chsi))
                Tau[s][i]= 1e-4*c_section[s]*np.trapz(function, x=z*1e3, dx=dz)
                i+=1
            s+=1
            

        Tau_tot=np.sum(Tau, axis=0)
        I=I_infinity*np.exp(-Tau_tot)
    
    return I


def IonProducRate(df,cross_ion,I_Ly,Wavelength):
    #I_Ly=W/m2/nm
    #Wavelength=A
    #cross_ion=cm2
    qs=df*cross_ion*I_Ly*1e-4*Wavelength*1e-1*(Wavelength*1e-10/(6.626*1e-34*3*1e8)) #cm3/s
    return qs
    
    
def Scaterring_Depth(alt, Tot_n):
    s=[]    
    for i in np.arange(len(alt)):
            z=np.array(alt[i:-1])
            rho=Tot_n[z]*1.66e-27*1e3
            s_val=np.trapz(rho, x=z)*1e5
            s.append(s_val)     
            i+=1
    return s
 
def Ionization_Rate(df,E, scattering_Depth, EnDissip_fct,Tot_n,delta_E,F):
    
    ind=np.where(EnDissip_fct.index >0)[0]
    f= interp1d(EnDissip_fct.index[ind], EnDissip_fct[EnDissip_fct.columns[0]][EnDissip_fct.index[ind]],
             fill_value = "extrapolate")
    rho=Tot_n*1.66e-27*1e3
    
    R_Ep=4.3e-7+5.36e-6*((E*1e-3)**1.67)
    R_ratio=np.array(scattering_Depth)/R_Ep
    R_ratio[R_ratio > 1] = 0
    Y=f(R_ratio)
    ind_Notzer0=np.where(R_ratio!=0)[0]
    Ionization_Rate= (E*F*Y[ind_Notzer0]*rho[ind_Notzer0])/(R_Ep*delta_E)

    return Ionization_Rate, df.index[ind_Notzer0]
