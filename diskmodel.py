import numpy as np
import astropy.units as u
import astropy.constants as c

def ssdisk(mass,medd,wav):
    # unit: lightday; maybe i will change it better someday
    # And there are more accurate correlations that are less significant 
    lightday = c.c*u.d
    wav = wav*u.AA
    mass = mass*u.M_sun
    X = 2.49 # X = 4.87
    eta = 0.1
    k = c.k_B
    h = c.h
    G = 6.6743*10**(-11)*u.m**3/(u.kg*u.s**2)
    sigma_sb = 5.6703744*10**(-8)*u.W/(u.m**2*u.K**4)
    Ledd = 4*np.pi*G*mass*c.c/(0.4*(u.cm)**2/u.g)
    r = (X*k*wav/(h*c.c))**(4/3) * ((4*G*mass/(8*np.pi*sigma_sb)) * (Ledd/(eta*c.c**2)) * medd)**(1/3)
    return (r.to(lightday)).value