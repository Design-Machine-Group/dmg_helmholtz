import json
import numpy as np 
import scipy as sp
from math import sqrt
from scipy.linalg import solve
from scipy.linalg import lstsq
from spatial import make_othro_array
from spatial import point_array_distance
from spatial import make_spatial_matix
from spatial import point_to_points_distance

from plotting import plot_grid

'''Helmholtz resonator array calculation based on Bart Van Der Aaa's thesis (2010)'''
#TODO: neck resistance values still differ slightly from Barts chart
#TODO: make sure the S matrix is in the right orientation (See bellow)
#TODO: this entire code needs to be checked, negative pressures???
#TODO: make a rhino check of distances and indices
#TODO: reproduce Bart's experiments, TF???, IL???
#TODO: check the np.solve. Probably why signs are weird???
#TODO: Why is qhr so small compared to q0???
#TODO: Are the distances correct??? Check with Rhino
#TODO: Try non uniform HR geometries

class HelmholtzPanel(object):
    def __init__(self, f, num_hr_x, num_hr_y, dx, dy):
        # HR dimensions - - - 
        self.nl = None
        self.nr = None
        self.bl = None
        self.br = None
        self.sn = None
        self.sb = None
        self.bv = None

        # Panel parameters - - - -
        self.num_hr_x = num_hr_x
        self.num_hr_y = num_hr_y
        self.dx = dx
        self.dy = dy
        self.num_hr = self.num_hr_x * self.num_hr_y
        self.D = None
        self.SN = None
        self.sD = None

        # scene - - - - 
        self.src_xyz = None
        self.q0 = None
        self.recs = None

        # acoustic parameters - - - - - - - - - - - - - -
        self.f       = f             # frequency (Hz)
        self.c       = 344.3            # speed of sound (m/s)
        self.rho     = 1.205            # density of air (kg/m3c)
        self.mu      = 1.88e-5          # dynamic viscosity of air (Ns / m2)
        self.gamma   = 1.4              # ratio of specific heats (-)
        self.thcon   = .026             # thermal conductivity (W / mk)
        self.cp      = 1010.            # heat capacity at constant pressure (J / kgK)

        self.wlen = self.c / self.f                # wavelength (m)
        self.afreq = 2. * np.pi * self.f       # angular frequency (rad/s)
        self.k = (2. * np.pi) / self.wlen     # wavenumber (rad/m)

        # - - -
        self.nr = None
        self.p = None

    def compute_hr_distances(self):

        self.hr_pts = make_othro_array(self.num_hr_x, 
                                       self.num_hr_y, 
                                       self.dx, 
                                       self.dy, 
                                       sp=[0, 0, 0])
        self.D = point_array_distance(self.hr_pts)
    
    def compute_panel_recs_distances(self):
        self.rD = point_to_points_distance(self.recs, self.hr_pts)
        self.src_rec_d = point_to_points_distance(self.src_xyz, self.recs)

    def compute_panel_src_distances(self):
        self.sD = point_to_points_distance(self.src_xyz, self.hr_pts)

    def compute_si_rec_distance(self):
        s = self.src_xyz.tolist()[0]
        self.source_image = np.array([[s[0], s[1], -s[2]]])
        self.si_D = point_to_points_distance(self.source_image, self.recs)

    def compute_panel_areas(self):
        self.sn = np.pi * (self.nr ** 2) # neck opening surface (m2)
        self.sb = np.pi * (self.br ** 2) # body opening surface (m2)

    def hr_resonant_freq(self):
        lef = self.nl + (1.7 * self.nr)  # neck effective length
        self.bv = self.bl * np.pi * self.br ** 2  # cavity volume
        self.resonant_f = (self.c / (2 * np.pi)) * np.sqrt( self.sn / (lef * self.bv))

    def calculate_mass_corrections(self):
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

        self.neck_len_correction()
        self.mn = self.sn * self.rho * self.nl
        self.body_mass_correction()
        self.m = self.mn + self.mb
        self.vol = (self.sb * self.bl) # (sn * nl) + (sb * bl)     
        
    def neck_len_correction(self):
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
        linner = .48 * np.sqrt(self.sn) * (1 - 1.25 * (self.nr / self.br))
        louter = (8 / (3 * np.pi)) * self.nr
        self.nl += linner + louter

    def body_mass_correction(self):
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
        self.mb = (1 / 3) * ((self.rho * self.bl) / self.sb) * (self.sn ** 2)

    def neck_resistance(self):
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
        mueff = self.mu * ((1 + ((self.gamma - 1) * sqrt(self.thcon / (self.mu * self.cp)))) ** 2)
        nwall = (self.nl / self.nr) * (sqrt(2 * mueff * self.rho * self.afreq) / (np.pi * self.nr ** 2))
        nends = 2 * (sqrt(2 * self.mu * self.rho * self.afreq) / (np.pi * self.nr ** 2))
        self.rn = nwall + nends
                
    def spring_stiffness(self):
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
        self.s = (self.rho * (self.c ** 2) * (self.sn ** 2)) / self.vol

    def hr_impedance(self):
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
        self.zhr = ((1j * self.afreq * self.m) + self.rn - ((1j / self.afreq) * self.s)) / self.sn

    def compute_spatial_neck_surface(self):
        self.SN = make_spatial_matix(self.sn) # not sure this matrix is in the right direction

    def coupling_functions(self):
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
        self.G = ((1j * self.afreq * self.rho) / ( 2 * np.pi)) * self.greens(self.D, self.k) * self.SN

    def coupling_pressures(self):
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
        diag = self.zhr - np.real(self.rho * self.c * (1 -  np.exp(-1j * self.k * self.nr)))
        np.fill_diagonal(self.G, diag)
        
    def calculate_direct_qs(self):
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
        self.c_ =  self.q0 * 2 * self.greens(self.sD, self.k).transpose()

    def calculate_velocities(self):
        self.v_ = solve(self.G, self.c_).flatten()
        # self.v_ = lstsq(self.G, self.c_)[0].flatten()

    def calculate_coupled_hr_strenghts(self):
        """
        """
        self.qhr = ((1j * self.rho * self.sn) / 2 * np.pi) * self.v_
        
    def calculate_pressure(self):
        """
        """
        p_direct = self.q0 * self.greens(self.src_rec_d, self.k)  # direct component
        p_hr_panel  = self.qhr * self.greens(self.rD, self.k)  # without direct component
        p_si = self.q0 * self.greens(self.si_D, self.k) * .8  # source image contribution, hard coded coeff. 
        self.p = np.sum(p_direct) + np.sum(p_hr_panel) + np.sum(p_si)  # all
        # self.p = np.sum(p_hr_panel)  # panel only
        # self.p = np.sum(p_direct)  # direct only
        # self.p = np.sum(p_si)  # source image only
        
    def greens(self, r, k):
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

    def compute_panel_impedance(self):
        # start calling functions - - - 
        self.compute_hr_distances()
        self.compute_panel_areas()
        self.hr_resonant_freq()  # optional for the moment 
        self.calculate_mass_corrections()

        # reistance, stiffness and impedance - - - - - - 
        self.neck_resistance()  # neck resitance (kg/s)
        self.spring_stiffness()           # spring stiffness
        self.hr_impedance()         # hr impedance 

    def compute_panel_coupling(self):
        self.compute_spatial_neck_surface()
        self.coupling_functions()
        self.coupling_pressures()

    def compute_velocities(self):
        self.compute_panel_src_distances()
        self.calculate_direct_qs()
        self.calculate_velocities()

    def compute_rec_pressure(self):
        self.calculate_coupled_hr_strenghts()
        self.compute_panel_recs_distances()
        self.compute_si_rec_distance()
        self.calculate_pressure()


# exporting to json - - - - - - - - - - - - -
def data_to_json(filepath, panel, pressures, recs):
    data = {}
    data['src'] = panel.src_xyz.tolist()[0]
    data['hr_pts'] = panel.hr_pts.tolist()
    data['recs'] = recs
    data['pressures'] = pressures
    data['nr'] = panel.nr.tolist()
    data['nl'] = panel.nl.tolist()
    data['br'] = panel.br.tolist()
    data['bl'] = panel.bl.tolist()


    with open(filepath, 'w+') as fp:
        json.dump(data, fp)

if __name__ == '__main__':
    
    f = 1032.
    num_hr_x = 2
    num_hr_y = 2
    dx = .1
    dy = .1

    panel = HelmholtzPanel(f, num_hr_x, num_hr_y, dx, dy)

    # scene parameters - - - - - - - - - - - - - - - - -
    panel.src_xyz = np.array([[-2.5, .55, 2]])
    panel.q0 = .1

    # hr parameters - - - - - - - - - - - - - - - - - -
    panel.nl = .02 * np.ones(panel.num_hr)       # neck length
    panel.nr = .01 * np.ones(panel.num_hr)       # neck radius
    panel.bl = .06 * np.ones(panel.num_hr)       # body length
    panel.br = .035 * np.ones(panel.num_hr)      # body radius

    # panel.nr[25] = .05
    # panel.nl[25] = .01
    # panel.br[25] = .1
    # panel.bl[25] = .1

    # panel.nr[55] = .05
    # panel.nl[55] = .01
    # panel.br[55] = .1
    # panel.bl[55] = .1

    ## HR resonant frequencies - - - - - - - - - - - - -
    # panel.compute_panel_areas()
    # panel.hr_resonant_freq()
    # print(panel.resonant_f)

    ## start calling methods - - - - - - - - - - - - - - -
    panel.compute_panel_impedance()
    panel.compute_panel_coupling()
    panel.compute_velocities()

    n = 2
    # k = panel.k
    pressures = []
    r_pressures = []
    recs = []
    delta  = 1.01

    for z in range(n):
        tpress = []
        for x in range(n):
            panel.recs = np.array([[.5, x * delta, .03 + z * delta]])
            panel.compute_rec_pressure()
            tpress.append(np.real(panel.p))
            r_pressures.append(np.real(panel.p))
            recs.append(panel.recs.tolist()[0])
        pressures.append(tpress)

    for z in range(n):
        # tpress = []
        for x in range(n):
            panel.recs = np.array([[x * delta, .5, .03 + z * delta]])  
            # panel.recs = np.array([[.7, x * delta, .03 + z * delta]])
            panel.compute_rec_pressure()
            # tpress.append(np.real(panel.p))
            r_pressures.append(np.real(panel.p))
            recs.append(panel.recs.tolist()[0])
    #     # pressures.append(tpress)

    minmax = 4e-2
    plot_grid(np.array(pressures), vmin=None, vmax=None)
    # plot_grid(np.array(pressures), vmin=-minmax, vmax=minmax)
    filepath = '/Users/time/Documents/UW/04_code/dmg_helmholtz/data/panel_geometry.json'
    data_to_json(filepath, panel, r_pressures, recs)