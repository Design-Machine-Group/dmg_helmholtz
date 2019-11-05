import numpy as np 
import scipy as sp
from math import sqrt


'''Helmholtz resonator array calculation based on Bart Van Der Aaa's thesis'''

def greens(r, k):
    """This is a greens function as described by [Bart2019]

    Parameters:
    r (float) - The distance between source and reciever (m)
    k (float) - Wavenumber ()

    Returns:
    (float) The greens function something.
    """
    return np.exp(-1j * k * r) / r

def hr_impedance(afreq, m, rn, s, sn):
    zhr = ((1j * afreq * m) + rn - ((1j / afreq) * s)) / sn
    return zhr

def neck_len_correction(sn, nr, br):
    linner = .48 * sqrt(sn) * (1 - 1.25 * (nr / br))
    louter = (8 / (3 * np.pi)) * nr
    return linner, louter

def body_mass_correction(rho, bl, sb, sn):
    mb = (1 / 3) * ((rho * bl) / sb) * (sn ** 2)
    return mb

def spring_stiffness(rho, c, sn, vol):
    s = (rho * (c ** 2) * (sn ** 2)) / vol
    return s

def neck_resistance(mu, gamma, thcon, cp, nl, nr, rho, afreq):
    #TODO neck resistance values still differ slightly from Barts chart
    mueff = mu * ((1 + ((gamma - 1) * sqrt(thcon / (mu * cp)))) ** 2)
    nwall = (nl / nr) * (sqrt(2 * mueff * rho * afreq) / (np.pi * nr ** 2))
    nends = 2 * (sqrt(2 * mu * rho * afreq) / (np.pi * nr ** 2))
    rn = nwall + nends
    return rn

def neck_resistance_paper(nl, nr, mu, rho, afreq):
    rn = (2 * (nl / nr)) * (sqrt(2 * mu * rho * afreq) / (np.pi * nr ** 2))
    return rn


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    flist = range(1, 2020, 20)
    # flist = [200]
    rnlist = []
    zhrlist = []
    for f in flist:

        # acoustic parameters - - -
        #f       = 250.             # frequency (Hz)
        c       = 344.3             # speed of sound (m/s)
        rho     = 1.205             # density of air (kg/m3c)
        mu      = 1.88e-5           # dynamic viscosity of air (Ns / m2)
        gamma   = 1.4               # ratio of specific heats (-)
        thcon   = .026              # thermal conductivity (W / mk)
        cp      = 1010.              # heat capacity at constant pressure (J / kgK)

        wlen = c / f                # wavelength (m)
        afreq = 2 * np.pi * f       # angular frequency (rad/s)
        k = (2. * np.pi) / wlen     #wavenumber (rad/m)

        # hr parameters - - -
        nl = .004        # neck length
        nr = .002       # neck radius
        bl = .016       # body length
        br = .0035        # body radius

        sn = np.pi * (nr ** 2) # neck opening surface (m2)
        sb = np.pi * (br ** 2) # body opening surface (m2)

        nl_= nl + (1.7 * nr)
        linner, louter = neck_len_correction(sn, nr, br)
        nl += louter + linner
        # nl = nl_
        # print(nl, nl_)



        mn = sn * rho * nl                              # neck mass (kg/m3)
        mb = body_mass_correction(rho, bl, sb, sn)      # body mass (kg/m3)
        m = mn + mb                                     # total mass (kg/m3)
        vol = (sb * bl)# (sn * nl) + (sb * bl)                     # volume (m3)

        rn = neck_resistance(mu, gamma, thcon, cp, nl, nr, rho, afreq)  # neck resitance (kg/s)
        #TODO neck resistance values still differ slightly from Barts chart
        # rn = neck_resistance_paper(nl, nr, mu, rho, afreq)
        rnlist.append(rn)

        s = spring_stiffness(rho, c, sn, vol)           # spring stiffness
        zhr = hr_impedance(afreq, m, rn, s, sn)

        
        za = (rho * c ** 2) / (1j * afreq * vol)
        tf = (za * sn) / zhr
        zhrlist.append(np.real(zhr))
        

    # plt.plot(flist, zhrlist)
    # plt.xlabel('f')
    # plt.ylabel('TFcalc')
    # plt.grid()
    # plt.show()


    plt.plot(flist, rnlist)
    plt.xlabel('f')
    plt.ylabel('Rn')
    plt.grid()
    plt.show()


    # Questions:
        # when is the corected mass used?
        # how is the volume calculated? is it corrected??






