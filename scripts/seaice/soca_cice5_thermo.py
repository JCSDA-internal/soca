# -*- coding: utf-8 -*-
"""CICE5 Thermodynamic
"""

import numpy as np

rho_s = 330.0    #[kg/m3]
rho_i = 917.0    #[kg/m3]
c0   = 2106.0    #[J/kg/deg] specific heat of fresh ice at 0C
L0   = 3.34e5    #[J/kg] latent heat of fusion of fresh ice at 0C
cw   = 4218.0
mu   = 0.054     #[deg/ppt] liquidus ratio between the freezing temperature and salinity of brine.
K0   = 2.03      #[W/m/deg] conductivity of fresh ice
beta = 0.13      #[W/m/ppt] empirical constant
Ks   = 0.30      #[W/m/deg] Conductivity of snow 
hs_min = 1e-4    #[m] minimum snow depth (could be in cice namelist)
hi_min = 0.01    #[m] minimum sea-ice thickness (hard coded) 

def qs(Ts):
    """Enthalpy of snow given snow temperature
    
    Args:
        Ts   (float): Snow temperature

    Return: 
             (float): Enthalpy of snow     
    """
    return -rho_s*(-c0*Ts+L0)

def Ts(qs):
    """Inverse of qs(Ts): Temperature of snow given snow Enthalpy
    
    Args:
        qs   (float): Snow enthalpy

    Return: 
             (float): Snow temperature
    """
    return qs/(rho_s*c0)+L0/c0

def Tm(S):
    """Freezing temperature of seawater
    
    Args:
        S   (float): Salinity

    Return: 
             (float): Freezing temperature (linear liquidus approximation)     
    """
    return -mu*S

def qi(Ti,Si):
    """Enthalpy of sea-ice given sea-ice temperature and salinity
    
    Args:
        Ti   (float): Sea-Ice temperature
        Si   (float): Sea-Ice salinity

    Return: 
             (float): Enthalpy of Sea-Ice     
    """
    dum1=c0*(Tm(Si)-Ti)
    dum2=L0*(1-Tm(Si)/Ti)
    dum3=cw*Tm(Si)
    return -rho_i*(dum1+dum2-dum3)

def Ti(qi,Si):
    """Temperature of sea-ice given sea-ice enthalpy and salinity
    
    Args:
        qi   (float): Sea-Ice enthalpy
        Si   (float): Sea-Ice salinity

    Return: 
             (float): Temperature of Sea-Ice     
    """
    a=c0
    b=(cw-c0)*Tm(Si)-(qi/rho_i)-L0
    c=L0*Tm(Si)
    return -(b+np.sqrt(b**2-4*a*c))/(2*a)

def Ki(Ti,Si):
    """Thermal conductivity of sea-ice given sea-ice temperature and salinity (conduct='MU71')
    
    Args:
        Ti   (float): Sea-Ice temperature
        Si   (float): Sea-Ice salinity

    Return: 
             (float): Conductivity of Sea-Ice     
    """    
    return K0+beta*(Si/Ti)
