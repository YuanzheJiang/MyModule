# To generate a thin disk SED

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

import astropy.units as u
import astropy.constants as c

import unit_conversion as uc

'''
# variable
a = 0.7 # normalized spin
eta = 0.104 # given by a, but i don't know how to calculate
mass = 5e7
nmass = mass/10**8
medd = 0.108 # give by paper, we should calculate it latter;LEdd = 1.5 * 10^38(MBH/M)
ledd = (1.5*10**38)*mass
lum = medd*ledd
nlum = lum/10**44
'''

# the flux
def l_nu(r,nu,rmin,nmass,nlum,eta):
    '''
    l_nu(r,nu,rmin,nmass,nlum,eta)
    r - independent variable, radius
    nu - frequency
    rmin - inner radius
    nmass - normalize BH mass, M/10^8 M0
    nlum - normalize luminosity, L/(10**4 erg s-1)
    eta - efficiency
    '''
    k1 = (1.2692431*10**(-19))*(nu**3)*(nmass**2)
    k2 = (2.7226414*10**(-16))*nu*(eta*nmass**2/nlum)**(1/4)
    f = (1-(rmin/r)**(1/2))**(1/4)
    k3 = k2*r**(3/4)/f
    return (k1*r)*np.exp(-k3)/(1-np.exp(-k3)) 

def radius(mass,a,medd):
    '''
    mass - BH mass, M0
    a - normalized spin paramter, such as a = 0.7
    medd - Eddington ratio
    '''
    # rg = gm/c^2
    rg = c.G*mass*u.solMass/((c.c)**2)
    rg = (rg.to(u.cm)).value

    # inner radius, or inner most stable circular radius
    # \chi = 2a/R_s
    z1 = 1 + (1-a**2)**(1/3) *((1+a)**(1/3) + (1-a)**(1/3))
    z2 = np.sqrt(3*a**2+z1**2)
    if a >= 0:
        zin = 3+z2-np.sqrt((3-z1)*(3+z1+2*z2))
    elif a < 0:
        zin = 3+z2+np.sqrt((3-z1)*(3+z1+2*z2))

    termin = c.G*mass*u.solMass/((c.c)**2)
    termin = termin.to(u.cm)
    rin = termin.value*zin

    # outter radius, given by self-gravity limit. May be set to \infy
    # for convenience, we set alpha =0.1 first
    # if 0<alpha<1, then alpha_0.1^(14/27) ranges from 0 to 3.3; not influence much yet
    alpha = 0.1
    nmass = mass/10**8
    rout = (3*10**16) * (alpha/0.1) * medd**(8/27) * nmass**(1/27)

    return [rin/rg,rout/rg]


def ssdisk_sed(nu,rmin,rmax,nmass,nlum,eta):
    '''
    l_nu(r,nu,rmin,nmass,nlum,eta)
    r - independent variable, radius
    nu - frequency
    rmin - inner radius
    nmass - normalize BH mass, M/10^8 M0
    nlum - normalize luminosity, L/(10**4 erg s-1)
    eta - efficiency
    '''   

        
    res = integrate.quadrature(l_nu,rmin,rmax,args=(nu,rmin,nmass,nlum,eta),rtol=5, maxiter=10000,miniter=1000)
    
    return res
