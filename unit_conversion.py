import astropy.units as u
import astropy.constants as c
import numpy as np

def pc2cm(pc):
    pc = pc*u.pc
    cm = pc.to(u.cm)
    return [cm.value,np.log10(cm.value)]
def cm2pc(cm):
    cm = cm*u.pc
    pc = cm.to(u.pc)
    return [pc.value,np.log10(pc.value)]

def ryd2ev(ryd):
    ryd = ryd*u.Ry
    ev = ryd.to(u.eV)
    return [ev.value,np.log10(ev.value)]
def ev2ryd(ev):
    ev = ev*u.eV
    ryd = ev.to(u.Ry)
    return [ryd.value,np.log10(ryd.value)]

def nu2aa(nu):
    nu = nu*u.Hz
    aa = (c.c/nu).to(u.AA)
    return aa.value
def aa2nu(aa):
    aa = aa*u.AA
    nu = (c.c/aa).to(u.Hz)
    return nu.value

def energy2wav(eng,un1,un2):
    if un1 == 'ryd' or un1 == 'Ryd' or un1 == 'Rd':
        eng = eng*u.Ry
    elif un1 == 'eV' or un1 == 'ev':
        eng = eng*u.eV
    
    nu = eng/c.h
    nu = nu.to(u.Hz)
    res = 0
    if un2 == 'hz' or un2 == 'Hz':
        res = nu.value
    elif un2 == 'A' or un2 == 'a':
        aa = c.c/nu
        aa = aa.to(u.AA)
        res = aa.value
    elif un2 == '&':
        aa = c.c/nu
        aa = aa.to(u.AA)
        res = [nu.value,aa.value]
    
    return res

def wav2energy(wav,un1,un2):
    if un1 == 'hz' or un1 == 'Hz':
        nu = wav*u.Hz
    elif un1 == 'A' or un1 == 'a':
        wav = wav*u.AA
        nu = c.c/wav
    
    eng = c.h*nu
    if un2 == 'eV' or un2 == 'ev':
        eng = eng.to(u.eV)
        return eng.value
    elif un2 == 'ryd' or un2 == 'Ryd' or un2 == 'Rd':
        eng = eng.to(u.Ry)
        return eng.value