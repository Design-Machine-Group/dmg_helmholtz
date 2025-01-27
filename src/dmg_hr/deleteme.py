import numpy as np 
import scipy as sp
from math import sqrt

from plotting import plot_grid

'''Helmholtz resonator array calculation based on Bart Van Der Aaa's thesis (2010)'''
#TODO: neck resistance values still differ slightly from Barts chart
#TODO: make sure the S matrix is in the right orientation (See bellow)
#TODO: this entire code needs to be checked, negative pressures???
#TODO: make a rhino check of distances and indices
#TODO: reproduce Bart's experiments, TF???, IL???
#TODO: check the np.solve. Probably why signs are weird???




def hr_resonant_freq(sn, nl, nr, br, bl):
    s = np.pi * nr ** 2 # neck cross sectional area
    lef = nl + (1.7 * nr)  # neck effective length
    v = bl * np.pi * br ** 2  # cavity volume
    f = (c / (2 * np.pi)) * np.sqrt( s / (lef * v))
    return f

def greens(r, k):
    """This is a greens function as described in Var der Aa (2010) eqs 2.12 in page 14.
    Parameters
    ----------
    r (float): The distance between source and reciever (m)
    k (float): Wavenumber ()
    Returns
    -------
      (float): The greens function something.
    """
    return np.exp(-1j * k * r) / r

def coupling_function(afreq, rho, rmn, s):
    """ This is the couppling function between two HRs (m and n) at a distance of rmn.
    Var der Aa (2010) eqs 2.24 in page 18.
    Parameters
    ----------
    afreq   (float): Angular frequency = 2 pi f (rad / s)
    rho     (float): Air density (kg/m3)
    rms     (float): The distance between the center of tne neck openings for HRs m and n (m)
    s       (float): The neck opening surface for resonator m (m2)
    Returns
    -------
    gmn (   float): 
    """

    gmn = ((1j * afreq * rho) / (2 * np.pi)) * greens(rmn) * s
    return gmn

def hr_impedance(afreq, m, rn, s, sn):
    """This function computes the mechanical impedance of a Helmholtz
    resonator. Taken from Var der Aa (2010) eqs 2.10 in page 12.
    Parameters
    ----------
    afreq   (float): Angular frequency = 2 pi f (rad/s)
    m       (float): Total mass neck + body (kg/m3)
    rn      (float): Neck resistance (kg/s)
    s       (float): Spring stiffness (N/m)
    sn      (float): Neck opening surface (m2)
    Returns
    -------
    zhr     (float): The mechanical impedance of the HR (Ns/m)
    """
    #TODO: compare the results of this function againts something. Cox?
    zhr = ((1j * afreq * m) + rn - ((1j / afreq) * s)) / sn
    return zhr

def neck_len_correction(sn, nr, br, nl):
    """This function applies corrections to the neck length as
    described by Var der Aa (2010) eqs 2.2 and 2.3 in page 9.
    Parameters
    ----------
    sn (float): Neck opening surface (m2)
    nr (float): Neck radius (m)
    br (float): Body radius (m)
    nl (float): Neck length (m)
    Returns
    -------
    nl (float): Corrected neck length (m)
    """
    linner = .48 * np.sqrt(sn) * (1 - 1.25 * (nr / br))
    louter = (8 / (3 * np.pi)) * nr
    nl += linner + louter
    return nl

def body_mass_correction(rho, bl, sb, sn):
    """This function corrects the body cavity mass as
    described by Var der Aa (2010) eqs 2.4 in page 9.
    Parameters
    ----------
    rho     (float): Air density (kg/m3)
    bl      (float): Body length (m)
    sb      (float): Body opening surface (m2)
    sn      (float): Neck opening surface (m2)
    Return
    ------
    mb      (float): Corrected body mass (kg/m3)
    """
    mb = (1 / 3) * ((rho * bl) / sb) * (sn ** 2)
    return mb

def spring_stiffness(rho, c, sn, vol):
    """This function computes the air spring stiffness
    according to Var der Aa (2010) eqs 2.5 in page 10.
    Parameters
    ----------
    rho     (float): Air density (kg/m3)
    c       (float): Speed of sound (m/s)
    sn      (float): Neck opening surface (m2)
    vol     (float): Total resonator volume (m3)
    Returns
    -------
    s       (float): Air spring stiffness (N/m)
    """
    s = (rho * (c ** 2) * (sn ** 2)) / vol
    return s

def neck_resistance(mu, gamma, thcon, cp, nl, nr, rho, afreq):
    """
    This function computes the neck resistance in an HR using the forulation
    presented in Var der Aa (2010) eqs 2.6 and 2.7 in page 10. 
    Parameters
    ----------
    mu      (float): Dynamical viscosity of air (Ns / m2)
    gamma   (float): Ratio of specific heats (-)
    thcon    (flat): Thermal conductivity (W / mK)
    cp      (float): Heat capacity at constant pressure (J / kgK)
    nl      (float): Neck length (m)
    nr      (float): Neck radius (m)
    rho     (float): Air density (kg/m3)
    afreq   (float): Angular frequency = 2 pi f (rad / s)
    Returns
    -------
    rn      (float): neck resistance (kg/s)
    """
    #TODO neck resistance values still differ slightly from Barts chart
    mueff = mu * ((1 + ((gamma - 1) * sqrt(thcon / (mu * cp)))) ** 2)
    nwall = (nl / nr) * (sqrt(2 * mueff * rho * afreq) / (np.pi * nr ** 2))
    nends = 2 * (sqrt(2 * mu * rho * afreq) / (np.pi * nr ** 2))
    rn = nwall + nends
    return rn

def neck_resistance_paper(nl, nr, mu, rho, afreq):
    """Probably wont keep this function...
    """
    rn = (sqrt(2 * rho * afreq * mu ) / (np.pi * nr ** 2)) * (2 + (nl / nr))  
    return rn

def calculate_corrections(sn, nr, br, nl, rho):
    """This function computes the neck and mass corrections necesary
    for the calculation of HR neck resistances and impedances, as
    described by Var der Aa (2010).
    Parameters
    ----------
    sn      (float): Neck opening surface (m2)
    nr      (float): Neck radius (m)
    br      (float): Body radius (m)
    nl      (float): Neck length (m)
    rho     (float): Air density (kg/m3)
    Returns
    -------
    nl  (float): Corrected neck length (m)
    mn  (float): neck mass (kg/m3)
    mb  (float): Corrected body mass (kg/m3)
    m   (float): Total mass (kg/m3)
    """

    nl = neck_len_correction(sn, nr, br, nl)
    mn = sn * rho * nl
    mb = body_mass_correction(rho, bl, sb, sn)
    m = mn + mb
    return nl, mn, mb, m

def coupling_functions(afreq, rho, D, S, k):
    """This function computes the coupling functions as 
    described in Var der Aa (2010) eqs 2.24 in page 18.
    Parameters
    ----------
    afreq   (float): Angular frequency = 2 pi f (rad / s)
    rho     (float): Air density (kg/m3)
    D       (float): Distance matrix of all resonator centers (m)
    S       (float): Neck opening surface Matrix (m)
    k       (float): Wave number ()
    Returns
    -------
    G (float): Coupling functions (?)
    """
    G = ((1j * afreq * rho) / ( 2 * np.pi)) * greens(D, k) * S
    return G

def coupling_pressures(zhr, rho, c, k, nr, G):
    """This function computes the coupling pressures of an HR array based
    on Var der Aa (2010) eqs 2.30 in page 22.
    Parameters
    zhr (float): HR mechanical impedance (Ns/m)
    rho (float): Air desnity (kg/m3)
    c (float): Speed of sound (m/s)
    k (float): Wave number (-)
    nr (float): Neck radius (m)
    G (float): coupling functions (?)
    Returns
    G (float): Coupling pressures (?)
    """
    diag = zhr - np.real(rho * c * (1 -  np.exp(-1j * k * nr)))
    np.fill_diagonal(G, diag)
    return G

def calculate_direct_qs(q0, sD, k):
    """This function computes the source-to-hr pressure as
    described in Var der Aa (2010) eqs 2.25 in page 18. 
    Parameters
    ----------
    q0  (float): Initial source stength (pressure) (Pa)
    sD  (float): Distance from source to HR centers (m)
    k   (float): wavenumber (-)
    Returns
    -------
    source-to-hr pressure
    """
    return q0 * 2 * greens(sD, k)

def calculate_coupled_hr_strenghts(afreq, rho, sn, v_):
    """
    """
    qhr = ((1j * rho * sn) / 2 * np.pi) * v_
    return qhr

def calculate_pressure(q0, src_rec_d, qhr, rD, k):
    """
    """
    # p = q0 * greens(src_rec_d, k) + np.sum(qhr * greens(rD, k))
    p = np.sum(qhr * greens(rD, k))
    return p


if __name__ == '__main__':
    from scipy.linalg import solve

    from spatial import make_othro_array
    from spatial import point_array_distance
    from spatial import make_spatial_matix
    from spatial import point_to_points_distance
    
    pressures = []
    frequencies = []
    for x in range(10):
        temp_pressures = []
        for z in range(1, 10):
            # scene parameters - - - - - - - - - - - - - - - -
            src_xyz = np.array([[10,0,2]])
            q0 = .1
            recs = np.array([[x*.1, 0, z*.1]])  # this should be only one for now. dont know what happens if more

            # acoustic parameters - - - - - - - - - - - - - -
            f       = 1000.             # frequency (Hz)
            c       = 344.3             # speed of sound (m/s)
            rho     = 1.205             # density of air (kg/m3c)
            mu      = 1.88e-5           # dynamic viscosity of air (Ns / m2)
            gamma   = 1.4               # ratio of specific heats (-)
            thcon   = .026              # thermal conductivity (W / mk)
            cp      = 1010.              # heat capacity at constant pressure (J / kgK)

            wlen = c / f                # wavelength (m)
            afreq = 2. * np.pi * f       # angular frequency (rad/s)
            k = (2. * np.pi) / wlen     # wavenumber (rad/m)

            # hr array parameters - - - - - - - - - - - - - -
            num_hr_x = 20
            num_hr_y = 20
            dx = .1
            dy = .1
            num_hr = num_hr_x * num_hr_y
            points = make_othro_array(num_hr_x, num_hr_y, dx, dy)
            D = point_array_distance(points)

            # hr parameters - - - - - - - - - - - - - - - - - 
            nl = .004 * np.ones(num_hr)       # neck length
            nr = .002 * np.ones(num_hr)       # neck radius
            bl = .016 * np.ones(num_hr)       # body length
            br = .0035 * np.ones(num_hr)      # body radius

            sn = np.pi * (nr ** 2) # neck opening surface (m2)
            sb = np.pi * (br ** 2) # body opening surface (m2)

            hr_rf = hr_resonant_freq(sn, nl, nr, br, bl)
            # print(hr_rf)

            # masses, volume - - - - - - - - - - - - - - - -
            nl, mn, mb, m = calculate_corrections(sn, nr, br, nl, rho)
            vol = (sb * bl)# (sn * nl) + (sb * bl)          # volume (m3) WICH VOLUME SHOULD BE USED????

            # reistance, stiffness and impedance - - - - - - 
            rn = neck_resistance(mu, gamma, thcon, cp, nl, nr, rho, afreq)  # neck resitance (kg/s)
            s = spring_stiffness(rho, c, sn, vol)           # spring stiffness
            zhr = hr_impedance(afreq, m, rn, s, sn)         # hr impedance 

            # coupling pressures - - - - - - - - - - - - - -
            S = make_spatial_matix(sn)  # not sure this matrix is in the right direction
            G = coupling_functions(afreq, rho, D, S, k)
            A = coupling_pressures(zhr, rho, c, k, nr, G)

            # direct to HRs - - - - - - - - - - - - - - - - -
            sD = point_to_points_distance(src_xyz, points)
            c_ = calculate_direct_qs(q0, sD, k).transpose()

            # velocities - - - - - - - - - - - - - - - - - - -
            v_ = solve(A, c_).flatten()

            # coupled pressures, pressures - - - - - - - - - -
            qhr = calculate_coupled_hr_strenghts(afreq, rho, sn, v_)
            rD = point_to_points_distance(recs, points)
            src_rec_d = point_to_points_distance(src_xyz, recs)
            p = calculate_pressure(q0, src_rec_d, qhr, rD, k)
            
            temp_pressures.append(np.real(p))
            # pressures.append(np.real(p))
            # frequencies.append(f)
        pressures.append(temp_pressures)
    
    # # plot pressures - - - - - - - - -
    # import matplotlib.pyplot as plt

    # plt.plot(frequencies, pressures)
    # plt.xlabel('Frequencies (Hz)')
    # plt.ylabel('Pressures (Pa)')
    # plt.show()

    plot_grid(np.array(pressures))