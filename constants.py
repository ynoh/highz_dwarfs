"""
constant in cgs
"""
import math


## physical constants
pi = math.pi
c = 2.99792e10 # cm
h_planck = 6.626e-27 #  erg s
k_B = 1.380658e-16 # erg K
e_charge = 4.8e-10 #  electron charge (esu)
m_e = 9.1e-28 #  electron mass(g)
m_H = 1.673534e-24 #  hyerogen mass (g)
m_p = 1.6726e-24 #  proton mass (g)
m_n = 1.675e-24 #  neutron mass (g)
G_newton = 6.67259e-8 #  gravitational const (eyne/cm^2/g^2)
alpha = 1./137. #  fine structure constant
sig_T = 6.65e-25 # ; thomson cross section(cm^2)
sig_sb = 5.67051e-5 #  stefan-boltzman cst
a_radcst = 7.56591e-15 # ; 4sg/c raeiation constant

## for the unit conversion [to cgs, rad]
eV = 1.6022e-12 #  erg
deg = pi/180. # rad
yrs = 365.*24.*3600. #  s
Joul = 1.e7 #  1 J =1.e7 erg 
cal = 4.184 #  lcal=4.184J
km = 1.e5
AU = 1.49598e13 # 1au= cm
pc = 3.086e18 # cm
Mpc = 1.e6*pc

## convenient energy scale
E_HIion = 13.6 #  eV hyerogen ionization energy
E_HeIion = 24.6 # eV HeI ionization energy
E_HeIIion = 54.4 # ev HeII ionization energy
T_HIion = 157807. # or E_Hion*eV/k_B ## taken from Hui & Gnedin
T_HeIion = 285335.
T_HeIIion = 631515.

## useful numbers for astronomy
M_sun = 1.9889e33 # g
L_sun = 3.846e33 #  erg/s
R_sun = 7.e10 #  cm
Tc_sun = 1.5e7 #  K
R_earth = 6.378e8 #  cm


## cosmologycal constants
Omg_m = .27
Omg_b = .046
sig_8 = .8
n_tilt = .96
h_hubble = .71
Hubble0 = h_hubble*100.*km/Mpc # s^-1
rho_cr = 3./(8.*pi*G_newton)*Hubble0**2 ## at z=0
rho_cr_practicalUnit = 3./(8.*pi*G_newton)*(100.*km/Mpc)**2*Mpc**3/M_sun # M_sun*h_hubble^2/Mpc^3
hubble_time = 1./Hubble0 # s
Y_He = .24
mu_cos = 0.6


## some radiation constant
sig_HI = 6.3e-18 # photoionization cross section at 13.6eV cm^-2
