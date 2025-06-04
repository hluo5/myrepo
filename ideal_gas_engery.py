#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 12 2025

@author: Jie Deng & Haiyang Luo
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

# parameters
h_mass    = 1.66054e-27 # kg
ev2j      = 1.60218e-19
boltz_ev  = 8.61733e-5 #eV/K
boltz     = boltz_ev*ev2j
avogadro  = 6.022e23
planck    = 6.62607e-34 #si, j/hz
planck_ev = planck/boltz_ev

def _energy(atomic_mass, natoms,vol,temp):
    """
    Input 
    atomic mass : atomic mass #
    natoms     : # of atoms
    vol   : total volume in Ang^3
    temp  : temperature in K
    ------
    Return
    F : Helmholtz free energy in eV
    ------
    
    benchmark example from Dorner 2017, page 183
    atomic_mass  = 28.085
    natoms       = 64
    vol          = 20.49*1e-30*natoms
    temp         = 1687
    inverse_temp = 1/boltz/temp
    thermal_lambda = planck/(2*np.pi*h_mass*atomic_mass/inverse_temp)**.5
    F              = -1/inverse_temp*natoms*(np.log(vol/(thermal_lambda**3)/natoms) + 1)
    """
    vol = vol*1e-30
    inverse_temp = 1/boltz/temp
    thermal_lambda = planck/(2*np.pi*h_mass*atomic_mass/inverse_temp)**.5
    #F              = -1/inverse_temp*natoms*(np.log(vol/(thermal_lambda**3)/natoms) + 1) #Stirling approximation
    F              = -1/inverse_temp*(natoms*np.log(vol/thermal_lambda**3)-scipy.special.gammaln(natoms+1)) #exact equation
    return F/ev2j

def energy(atomic_mass, natoms,vol,temp):
    """
    wrapper function of _energy to handle >=1 element case
    
    benchmark example from Yuan & Steinle‐Neumann, 2020, Table S2
    ---------
    atomic_mass = 55.845
    natoms = 50
    temp = 4000
    vol = 424.19
    hf =energy(atomic_mass, natoms,vol,temp)
    print(hf)

    vol = 485.1
    natoms = np.array([15, 15, 45])
    atomic_mass = np.array([24.305, 28.085, 15.999])

    hf =energy(atomic_mass, natoms,vol,temp)
    print(hf)
    
    """
    try:
        nele     = len(atomic_mass)
        vol4each = vol*(np.ones(nele))#*(np.array(natoms)/sum(natoms))
        f4each   = np.zeros(nele)
        for i in range(nele):
            f4each[i] = _energy(atomic_mass[i], natoms[i], vol4each[i],temp) 
            f = sum(f4each)
    except:
        f = _energy(atomic_mass, natoms,vol,temp)
    return f

def TS(natoms,temp):
    nele     = len(natoms)
    lnfactorial4each   = np.zeros(nele)
    for i in range(nele):
            lnfactorial4each[i] = scipy.special.gammaln(natoms[i]+1)
    lnfactorial = sum(lnfactorial4each)
    total_atoms = sum(natoms)
    lnfactorial_total = scipy.special.gammaln(total_atoms+1)
    TS = boltz_ev*temp*(lnfactorial_total-lnfactorial)
    return TS

#H = 3.555555556
#ts=3000*boltz_ev*(scipy.special.gammaln(160+H+1)-scipy.special.gammaln(32+1)-scipy.special.gammaln(32+1)-scipy.special.gammaln(96+1)-scipy.special.gammaln(H+1))
#ts=ts/(32+H)
#print(ts)

#atomic_mass = 28.085
#natoms = 2
#temp = 1687
#vol = 10.53896**3 
#hf = energy(atomic_mass, natoms,vol,temp)/64
#hf_ = _energy(atomic_mass, natoms,vol,temp)
#print("Si:",hf_)
#print("Si system",hf)


# atomic_mass = 24.305
# natoms = 64
# temp = 2850
# vol = 1903.36
# hf =energy(atomic_mass, natoms,vol,temp)
# print(hf)

#atomic_mass = [55.845,30.974]
#natoms = [50,1]
#temp = 3000
#vol = 588.4934
#hf =energy(atomic_mass, natoms,vol,temp)
#print(hf)

#atomic_mass = [55.845,15.999,30.974]
#natoms = [50,0,1]
#temp = 3000
#vol = 588.4934
#hf =energy(atomic_mass, natoms,vol,temp)
#print(hf)

atomic_mass = [24.305, 28.085,15.999]
natoms = [32,32,96]
temp = 3000
vol = 1736.92749
hf =energy(atomic_mass, natoms,vol,temp)
print(hf)
TS = TS(natoms,temp)
print(TS)

##### WCA liquid papram ref mirzaeinia17 Table I### 

a1j = [-0.38722877,	-1.53187761, 2.36329135, -0.32895711]
a2j = [0.13526998, -0.63646675, 2.26372117, -1.22784362]
a3j = [0.03792191, -0.18568084, 0.09014619, 0.83258622]
a4j = [0, 0, 0, 0]
a5j = [0, 0, 0, 0]
b  = [2.12271065,	1.45331509,	1.34899324,	-0.29320675,-0.08480989]

aij = np.array([a1j,a2j, a3j, a4j, a5j])

# rho = 1
# T = 1
# z0, u = 0, 0, 0, 0
# for i in range(1,4):
#     for j in range(1,5):
#         ii = i - 1
#         jj = j - 1
#         z0 += 2*i*aij[ii,jj]*rho**(2*i)*(T**((2-j)/2))
# ii = np.array(range(0,5))
# Z = z0 + sum(b*(ii+1)*rho**(ii+1)) + 1

# for i in range(1,4):
#     for j in range(1,5):
#         ii = i - 1
#         jj = j - 1
#         u += -(2-j)/2*aij[ii,jj]*(rho**(2*i))*(T**((2-j)/2)) # compressibility factor




def _cal_int_helper(rho,T):
    """
    calculate the 

    Parameters
    ----------
    rho : float
        reduced density = N*sigma^3/V
    T : float
        reduced temperature = kT/eps

    Returns
    -------
    integrand : float
        (Z-1)/rho
    Ref
    -------
    eqn 11
    mirzaeinia17

    """
    integrand = 0
    for i in range(1,4):
        for j in range(1,5):
            ii = i - 1
            jj = j - 1
            integrand += 2*i*aij[ii,jj]*rho**(2*i-1)*(T**((2-j)/2))
    ii = np.array(range(0,5))
    integrand += sum(b*(ii+1)*rho**(ii+1-1))
    return integrand

def cal_int(rho,T):
    """
    

    Parameters
    ----------
    rho : float
        target reduced density = N*sigma^3/V
    T : float
        reduced temperature = kT/eps

    Returns
    -------
    fres: float
        residual helmholtz free energy
        LFS of eqn 17 of mirzaeinia17
     
    benchmark with Sun et al., 2018
    -------
    rho =  0.1*6/np.pi
    T   = 1.5
    rho_in = np.linspace(0.001, rho)
    int_out = np.zeros(rho_in.shape)
    for i in range(len(rho_in)):  
        int_out[i]=cal_int(rho_in[i],T)
        
    plt.plot(rho_in,int_out)

    sum(int_out)*(rho_in[1]- rho_in[0])

    """
    rho_in = np.linspace(1e-6, rho)
    int_out = np.zeros(rho_in.shape)
    for i in range(len(rho_in)):  
        int_out[i]=_cal_int_helper(rho_in[i],T)    
    
    fres = sum(int_out)*(rho_in[1]- rho_in[0])
    return fres

eta = 0.1 # packing fraction, encoded in modify vasp script finite_diff.F ! eta = 0.1, (eta*6*vol/NIONS/3.14159265359)**(1.0/3.0), ref. Sun et al
rho =  eta*6/np.pi
T   = 1.5 # reduced temperature, encoded in modify vasp script finite_diff.F 

print(cal_int(rho,T))
