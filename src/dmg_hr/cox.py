import numpy as np
import matplotlib.pyplot as plt

# TO BE CHECKED, contrast with ML, use is equal function


#Absorption of a perforated absorber
#Normal incidence

c = 340.                    # speed of sound 
rho = 1.21                  # density of air 
Z0 = c * rho
viscosity = 15e-6           # kinemetric viscosity of air

sigma = 20000               # Flow resistivity of mineral wool
l1 = 0.025                  # backing thickness air
l2 = 0.025                  # backing thickness porous absorber


f = list(range(100, 2550, 50))   # Frequency
f = np.array(f)
nf = len(f)
kair = 2 * np.pi * f / c
w = 2 * np.pi * f

#Impedance at top of air layer
z1 = -1j * Z0 * (1. / np.tan((kair * l1)))


#calculate impedance of porous material (Delany and Bazley)
#dimensionless quantity for Delany and Bazley
X = rho * f / sigma


#characteristic impedance
Zc = rho * c * (1 + 0.0571 * (X**-0.754)-1j*0.087*(X**-0.732))
k = (2*np.pi/c) * f *(1+0.0978*(X**-0.700)-1j*0.189*(X**-0.595)) #wavenumber


#Impedance at top of porous absorbent
z2 = (-1j*z1*Zc* (1. / (np.tan(k*l2)))+Zc**2 ) / (z1 -1j*Zc*( 1. / (np.tan((k*l2)))))


# Loop over different open areas
eta = [0.0625, 0.125, 0.25, 0.50, 1.00]
ne = len(eta)


for m in eta:
    a = 2.5e-3                                          # hole radius
    D = np.sqrt(np.pi/m)*a                              # Hole spacing
    delta = 1.6*(1-1.47*m**0.5+0.47*m**3/2)             # end correction       
    t = 6.3e-3                                          # plate thickness
    rm = (rho/m)*np.sqrt(8*viscosity*w)*(1+t/(2*a))     # surface resistance                               
    z3 = (1j/m)*(2*delta*a+t)*w*rho+z2+rm               # impedance of resonant absorber
    R = (z3-rho*c)/(z3+rho*c)                           # reflection factor
    alpha = 1-abs(R)**2                                 # absorption coefficient
    # plt.plot(f, np.real(z3))
    plt.plot(f, alpha)



plt.xlabel('f')
plt.ylabel('Absorption Coefficient')
plt.grid()
plt.show()