"""
06/02/2012 Convert the mathematica code writen by Matt to python code

Edit the code to plot analytic curves with data

Editted based on Gnedin & Hui 97
"""

import os
import numpy as N
import matplotlib.pylab as plt
import constants as cst
from scipy.optimize import fsolve
from matplotlib.font_manager import FontProperties

## cosmology 
def E_z(z): ## ignore radiation at this point
    Omg_L = 1. - cst.Omg_m
    return N.sqrt(cst.Omg_m*(1. + z)**3. + Omg_L)

def mean_hydrogen_numden_z(z): ## cm^-3 : mean density at the given redshift
	#rho_cr = cst.rho_cr/cst.h_hubble**2*hubble**2 
	## cst.h_hubble is set by the hubble cst used in sim.
	return cst.rho_cr/cst.m_p*(1. - cst.Y_He)*cst.Omg_b*(1. + z)**3 

def mean_helium_numden_z(z): ## cm^-3 : mean density at the given redshift
	rho_cr = cst.rho_cr#/cst.h_hubble**2*hubble**2
	return rho_cr/cst.m_p*(cst.Y_He/4.)*cst.Omg_b*(1. + z)**3


def hydrogen_numden_to_baryon_den(n_H):
	#rho = n_H*cst.m_p ## not including He
	return n_H*cst.m_p/(1. - cst.Y_He)  ## to include He


def omega_m_z(zred):
	Omg_L = 1. - cst.Omg_m
	omgm = cst.Omg_m*(1. + zred)**3/(cst.Omg_m*(1. + zred)**3 + Omg_L)
	return omgm


def virial_overdensity_z(zred):
	"""
	from Bryan and Norman 1998
	"""
	omgmz = omega_m_z(zred)
	ovden_c = 18.*N.pi**2 + 82.*(omgmz - 1.) - 39.*(omgmz - 1.)**2
	return ovden_c



##### jeans #####

def sound_speed(temp, gamma):
	return N.sqrt(gamma*cst.k_B*temp/(cst.mu_cos*cst.m_p))

def calc_sound_speed_temp(c_s, gamma):
	"""
	get the temperature corresponding to the given sound speed
	"""
	return c_s**2*cst.mu_cos*cst.m_p/gamma/cst.k_B


def jeans_scale(temp, n_H, gamma=5./3.):
	"""
	k_J = (4 pi G rho)^0.5 / c_s
	L_J = 2 pi / k_J
	M_J propto rho(L_J/2.)**3
	return in M_sun 
	"""
	#rho = n_H*cst.m_p
	rho = hydrogen_numden_to_baryon_den(n_H)*cst.Omg_m/cst.Omg_b 
	c_s = sound_speed(temp, gamma)
	k_J = N.sqrt(4.*cst.pi*cst.G_newton*rho)/c_s
	return k_J


def jeans_mass(temp, n_H, gamma=5./3.):
	"""
	k_J = (4 pi G rho)^0.5 / c_s
	L_J = 2 pi / k_J
	M_J propto rho(L_J/2.)**3
	return in M_sun 
	"""
	rho = hydrogen_numden_to_baryon_den(n_H)*cst.Omg_m/cst.Omg_b 
	k_J = jeans_scale(temp, n_H, gamma=gamma)
	m_J = 4.*cst.pi/3.*rho*(2.*cst.pi/k_J/2.)**3
	#m_J = 4.*cst.pi/3.*rho*(2.*cst.pi/k_J)**3
	return m_J/cst.M_sun


def calc_jeans_mass_temp(n_H, m_J, gamma=5./3.):
	"""
	m_J is in the unit of M_sun
	get the temperature corresponding to the given jeans mass and the density
	"""
	m_J *= cst.M_sun
	rho = hydrogen_numden_to_baryon_den(n_H)*cst.Omg_m/cst.Omg_b 
	inv_ff_time = N.sqrt(4.*cst.pi*cst.G_newton*rho)
	c_s = (m_J/(4.*cst.pi/3.)/rho)**(1./3.)*2.*inv_ff_time/(2.*cst.pi)
	temp = calc_sound_speed_temp(c_s, gamma)
	return temp


def free_fall_time(n_H):
	rho_b = hydrogen_numden_to_baryon_den(n_H)
	rho_m = cst.Omg_m/cst.Omg_b*rho_b
	return (8.*cst.pi*cst.G_newton*rho_m/3.)**(-0.5)

##### self-shielding #####

def column_density_at_jeans_scale(temp, n_H, Gamma12):
	"""
	self shielding happens when tau = 1
	tau = N_H (over the Jeans scale) * photo ionization cross section
	"""
	f_g = cst.Omg_b/cst.Omg_m
	column_den_HI = 2.3e13*(n_H/1.e-5)**1.5*(temp/1.e4)**(-0.26)*Gamma12**(-1.)*N.sqrt(f_g/0.16)
	return column_den_HI 


def self_shielding(temp, n_H, Gamma12, sigma=cst.sig_HI):
	return column_density_at_jeans_scale(temp, n_H, Gamma12)*sigma - 1.
	

def calc_self_shielding_temp(n_H, Gamma12, sigma=cst.sig_HI):
	f_g = cst.Omg_b/cst.Omg_m
	xx = 2.3e13*(f_g/0.16)**0.5
	temp4 = (xx*sigma*(n_H/1.e-5)**1.5/Gamma12)**(1./0.26)
	#print (xx*cst.sig_HI*(1./1.e-5)**1.5/Gamma12)**(1./0.26)
	return temp4*1.e4


##### cooling and heating #####

## ratio of transition temperature to temperature
## it is defined as lambda in Hui & Gnedin &97)
get_Tion_to_T = lambda temp_TR, temp: 2.*temp_TR/temp


## neutral fraction:
## the given expression is only valid when e.g. x_H << 1, X_H = n_HI/n_H
neutral_fraction = lambda n_e, Gamma, recomb_coeff: recomb_coeff*n_e/Gamma
#neutral_frac_H = lambda temp, n_e, Gamma: recomb_rate_H(temp)*n_e/Gamma
#neutral_frac_HeI = lambda temp, n_e, Gamma: recomb_rate_HeI(temp)*n_e/Gamma
#neutral_frac_HeII = lambda : 1.


## recombination rate/coefficient 
#recomb_rate_HeI = lambda temp: 3.e-14*(2.*285335./temp)**0.654
#recomb_rate_H = lambda temp: 2.51e-13*(temp/2.e4)**(-0.7)
def recombination_coefficient(temp, species, case):

	func1 = lambda ratT, c0, c1, a0, a1, a2: c0*ratT**a0/(1.+(ratT/c1)**a1)**a2

	valsdict = dict(A = dict(HI = [cst.T_HIion, 1.269e-13, 0.522, 1.503, 0.470, 1.923],
	HeI = [cst.T_HeIion, 3.e-14, 0.654], 
	HeII = [cst.T_HeIIion, 2.*1.269e-13, 0.522, 1.503, 0.470, 1.923]),
	B = dict(HI = [cst.T_HIion, 2.753e-14, 2.740, 1.50 , 0.407, 2.242],
	HeI = [cst.T_HeIion, 1.26e-14, 0.750], 
	HeII = [cst.T_HeIIion, 2.*2.753e-14, 2.740, 1.50 , 0.407, 2.242]))

	vals = valsdict[case][species]
	ratT = get_Tion_to_T(vals[0], temp)
	#print "{:s}: {:}".format(species, vals)
	if species == "HeI":
		ri = vals[1]*ratT**vals[2]
	else:
		ri = func1(ratT, vals[1], vals[2], vals[3], vals[4], vals[5])

	return ri


def recombination_cooling_rate(temp, species, case):
	"""
	caseA: recombination cooling when optically thin (no local absorption)
	caseB: when optically thick. local absorption exists
	"""
	func1 = lambda ratT, c0, c1, a0, a1, a2: c0*ratT**a0/(1.+(ratT/c1)**a1)**a2

	valsdict = dict(A = dict(HI = [cst.T_HIion, 1.778e-29, 0.541, 1.965, 0.502, 2.697]),
	#"HeI": [cst.T_HeIion, ], 
	B = dict(HI = [cst.T_HeIIion, 3.435e-30, 2.250, 1.970, 0.376, 3.720]))

	vals = valsdict[case][species]
	ratT = get_Tion_to_T(vals[0], temp)
	#print "{:s}: {:}".format(species, vals)
	if species == "HeI":
		rc = cst.k_B*temp*recombination_coefficient(temp, case, species)
	else:
		rc = temp*func1(ratT, vals[1], vals[2], vals[3], vals[4], vals[5])

	return rc 



def collisional_ionization_coefficient(temp, species):
	"""
	currently, don't use He
	"""

	func1 = lambda ratT, c0, c1, a0, a1, a2: c0*ratT**a0/(1. + (ratT/c1)**a1)**a2

	valsdict = dict(HI = [cst.T_HIion, 21.11, 0.354, -1.089, 0.874, 1.101], HeI = [cst.T_HeIion], HeII = [cst.T_HeIIion])

	vals = valsdict[species]	
	ratT = get_Tion_to_T(vals[0], temp)
	#print "{:s}: {:}".format(species, vals)

	ci = temp**(-1.5)*N.exp(-ratT/2.)*func1(ratT, vals[1], vals[2], vals[3], vals[4], vals[5])
	
	return ci


"""
def collisional_cooling_rate(temp):
	return c0*temp**(-1.5)*N.exp(-157807./temp)*(2.*157807./temp)**a0/(1.+(2.*157807./temp/c1)**a1)**a2
"""

def collisional_cooling_rate(temp, species):
	tiondict = dict(HI = cst.T_HIion, HeI = cst.T_HeIion, HeII = cst.T_HeIIion)

	return cst.k_B*tiondict[species]*collisional_ionization_coefficient(temp, species)


def line_excitation_cooling_rate(temp, species):
	valsdict = dict(HI = [cst.T_HIion, 7.5e-19, 0.], HeII = [cst.T_HeIIion, 5.54e-17, -0.397])
	vals = valsdict[species]
	ratT = get_Tion_to_T(vals[0], temp)
	#print "{:s}: {:}".format(species, vals)
	ec = vals[1]*temp**vals[2]*N.exp(-0.75*ratT/2.)/(1. + (temp/1.e5)**0.5)
	return ec


def cooling_rate(temp, n_H, Gamma12, recomb_cooling_case):
	"""
	currently only hydrogen
	have to be changed if other species are included
	"""
	species = ["HI"]
	coolingR = 0.
	for sp in species:
		ri = recombination_coefficient(temp, sp, recomb_cooling_case)
		x_sp = neutral_fraction(n_H, Gamma12*1.e-12, ri)
		## I think n_H is multiplied only once because of the cancelation of one with n_H in heating
		rec = recombination_cooling_rate(temp, sp, recomb_cooling_case)*(1. - x_sp)**2*n_H 
		#print rec
		col = collisional_cooling_rate(temp, sp)*x_sp*(1. - x_sp)*n_H
		lin = line_excitation_cooling_rate(temp, sp)*x_sp*(1. - x_sp)*n_H
		#print col
		#coolingR += (rec + col)/(3./2.*cst.k_B*2.)
		coolingR += (rec + col + lin)/(3./2.*cst.k_B*2.)
	return coolingR


def heating_rate_simplest(temp, n_H, Gamma12, alpha):
	"""
	assumed n_e ~ n_H
	Gamma12 doesn't affect the result at all. Due to x_H, it canceled out.
	"""
	ri_H = recombination_coefficient(temp, "HI", recomb_cooling_case)
	x_H = neutral_fraction(n_H, Gamma12*1.e-12, ri_H)
	#x_H = neutral_frac_H(temp, n_H, Gamma12*1.e-12) 
	ionizationE = cst.E_HIion*cst.eV
	heatingR = x_H*Gamma12*1.e-12*ionizationE/(2. + alpha)/(3./2.*2.*cst.k_B)
	return heatingR


def heating_rate(temp, n_H, Gamma12Arr, alpha, recomb_cooling_case):
	"""
	assumed n_e ~ n_H
	include heating by HeI, HeII
	"""
	## photo heating rate function
	photo_heating_rate = lambda g, en, p1, p2, f1, f2: g*en*(p1/p2*f1/f2 - 1.)

	species = N.array(["HI", "HeI", "HeII"])
	valsdict = dict(HI = [cst.E_HIion, 2.8, 1.8, 4.], HeI = [cst.E_HeIion, 1.7, 0.7, 2.2], 
	HeII = [cst.E_HeIIion, 2.9, 1.9, 0.])

	GammaArr = 1.e-12*Gamma12Arr

	heatingR = 0.
	for sp, ga in zip(species[:len(GammaArr)], GammaArr):
		vals = valsdict[sp]
		f1 = 1.; f2 = 1.
		if sp != "HeII":
			ri = recombination_coefficient(temp, sp, recomb_cooling_case)
			x_sp = neutral_fraction(n_H, ga, ri)
			f1 -= vals[3]**(-(vals[2] + alpha)); f2 -= vals[3]**(-(vals[1] + alpha))	
		else:
			x_sp = 1.
		pr = photo_heating_rate(ga, vals[0]*cst.eV, vals[1], vals[2], f1, f2)
		heatingR += x_sp*pr	
	
	heatingR /= 3.*cst.k_B
	return heatingR
	


def func_coolingtime_eqto_fftime(temp, n_H, Gamma12Arr, alpha, recomb_cooling_case):
	inv_cooling_time = cooling_rate(temp, n_H, Gamma12Arr[0], recomb_cooling_case) - heating_rate(temp, n_H, Gamma12Arr, alpha, recomb_cooling_case)
	inv_fftime = 1./free_fall_time(n_H)
	return inv_cooling_time - inv_fftime


def coolingtime_eqto_fftime_initemp(n_H, Gamma12Arr, alpha=0, recomb_cooling_case="A"):
	"""
	get the initial temperature at which cooling time is equal to free fall time
	better to take inverse of free fall time because cooling_rate could be "0" in the given range of temperature but free fall time shouldn't be "0"
	"""
	ntemp = 100
	inv_fftime = 1./free_fall_time(n_H) 

	tempArr = N.logspace(2., 6., num=ntemp)
	inv_cooling_time = cooling_rate(tempArr, n_H, Gamma12Arr[0], recomb_cooling_case) - heating_rate(tempArr, n_H, Gamma12Arr, alpha, recomb_cooling_case)
	if inv_fftime > N.min(inv_cooling_time) and inv_fftime < N.max(inv_cooling_time):
		diff = N.abs(inv_cooling_time - inv_fftime)
		ii = N.where(diff == N.min(diff))[0]
		order = int(N.log10(N.min(tempArr[ii])))
		#print N.min(diff)
		while N.min(diff) > 10.**(order - 1):
			if len(ii) > 1:
				tempArr = N.logspace(tempArr[N.min(ii), N.max(ii)], num=ntemp)
			else:
				tempArr = N.logspace(tempArr[ii-1], tempArr[ii+1], num=ntemp)
			inv_cooling_time = cooling_rate(Gamma12Arr[0], n_H, tempArr, recomb_cooling_case) - heating_rate(Gamma12Arr, alpha, n_H, tempArr, recomb_cooling_case)
			diff = N.abs(inv_cooling_time - inv_fftime)
			ii = N.where(diff == N.min(diff))[0]
			order = int(N.log10(N.min(tempArr[ii])))
		return tempArr[ii]
	## need to fix for else case
	else:
		#print "reset temperature range"
		#print "inverse cooling time range: [%.2e, %.2e] and inverse free fall time: %.2e" % (N.min(inv_cooling_time), N.max(inv_cooling_time), inv_fftime)
		return 1.e5 ## anyhow set initial temperature for fsolve
			
def calc_coolingtime_eqto_fftime_temp(n_H, Gamma12Arr, alpha=0., recomb_cooling_case="A", temp0=1.e4):
	if type(n_H) != N.ndarray:
		n_H = N.array([n_H])

	tempArr = N.zeros_like(n_H)
	for i, n_H1 in enumerate(n_H):
		#print i
		temp0 = coolingtime_eqto_fftime_initemp(n_H1, Gamma12Arr)
		tempArr[i] = fsolve(func_coolingtime_eqto_fftime, temp0, args=(n_H1, Gamma12Arr, alpha, recomb_cooling_case))
	return tempArr
		


##### adibatic #####

def calc_adiabatic_heating_temp(n_H, temp_i, n_H_i=None, gamma=5./3):	
	"""
	T propto rho^(gamma-1)
	"""
	if type(n_H) == float:
		n_H = N.array([n_H])

	if n_H_i == None:
		n_H_i = n_H[0]

	temp = temp_i*(n_H/n_H_i)**(gamma-1)
	return temp
	

def coord_write_text(n_HArr, tempArr1, refval, ref):
	if ref == "temp":
		arr = tempArr1
	else:
		arr = n_HArr
	diff = N.abs(arr - refval)
	ii = N.where(diff == N.min(diff))[0]
	if len(ii) > 1:
		ii = ii[0]

	print ii, N.min(diff), n_HArr[ii], tempArr1[ii]

	#if ii >= len(diff)-1:
	#	return n_HArr[ii], ref
	#else:
	#	return n_HArr[ii+1], tempArr1[ii]

	return n_HArr[ii], tempArr1[ii]



def filtering_mass(zcol, zrei, temp_rei, isothermal=True):
	if isothermal:
		func1 = lambda zcol, zrei: 3./10.*(1. + 4.*((1. + zcol)/(1. + zrei))**2.5 - 5.*((1.+ zcol)/(1. + zrei))**2)
		temp0 = temp_rei
	else:
		func1 = lambda zcol, zrei: 1. + 2.*((1. + zcol)/(1. + zrei))**1.5 - 3.*(1.+ zcol)/(1. + zrei)
		temp0 = temp_rei*(1. + zcol)/(1. + zrei)


	n_H = mean_hydrogen_numden_z(zcol)
	k_J = jeans_scale(temp0, n_H, gamma=5./3.)
	k_F = k_J/N.sqrt(func1(zcol, zrei))
	print k_J, k_F

	rho = hydrogen_numden_to_baryon_den(n_H)*cst.Omg_m/cst.Omg_b 
	##m_F = 4.*cst.pi/3.*rho*(2.*cst.pi/k_F)**3
	##return 8.*m_F/cst.M_sun ## a factor of 8 is due to our definition of Jeans mass
	## since our definition of jeans mass is 8 times smaller than what Gnedin used 
	m_F = 4.*cst.pi/3.*rho*(2.*cst.pi/k_F)**3
	return m_F/cst.M_sun ## a factor of 8 is due to our definition of Jeans mass
	


def accretion_mass(zcol, Gamma12, delta_eq=60., mu=0.6, alpha=0., recomb_cooling_case='A', temp0=1.e4):
#def accretion_mass(zcol, Gamma12Arr, delta_eq=60., mu=0.6, alpha=0., recomb_cooling_case='A', temp0=1.e4):
	"""
	T_vir = T_eq(60.*mean density) (60.\sim Delta_vir/3.
	T_vir = \mu m_p /k_b/2 (Delta_vir *Omg_m/2^(1/3) *(1+z) (G M H_0)^(2/3)
	Okamoto et al. 2008 MNRAS eq. 4
	"""

	print "CALCULATING ACCRETION MASS"
	delta_v = virial_overdensity_z(zcol)
	if delta_eq == 0.:
		delta_eq = delta_v/3.
		print "calculating virial overdensity: delta_eq=delta_v/3.=%.2f" % delta_eq

	n_H = delta_eq*mean_hydrogen_numden_z(zcol)

	## getting T_vir = T_eq(Delta_vir/3.)
	##old code
	#temp_eq = calc_coolingtime_eqto_fftime_temp(n_H, Gamma12Arr, alpha=alpha, recomb_cooling_case=recomb_cooling_case, temp0=temp0)
	##new: using Matt's calculation
	abspath_current = os.path.abspath(".")
	ii = abspath_current.find("highz_dwarfs")
	abspath_current = abspath_current[:ii+12]
	path = os.path.join(abspath_current, "analytic", "curves")
	filename = "Teq_G%.1f_full.dat" % (Gamma12)
	n_HArr, temp_eqArr = N.loadtxt(os.path.join(path, filename), usecols=(0, 1), unpack=True)
	dff = N.abs(n_HArr/n_H - 1.)
	imin = N.where(dff == N.min(dff))[0]
	temp_eq = temp_eqArr[imin]
	print "mean num den times %.0f :%.2e which the num den from the file: %.2e wants to be close to. The corresponding equilibrium temp from the file: %.1f" % (delta_eq, n_H, n_HArr[imin], temp_eq)
	if (imin > 0 and imin < (len(temp_eqArr) - 1)):
		print "values around the chosen in the file:[%.2e, %.1f], [%.2e, %.1f]" % (n_HArr[imin-1], temp_eqArr[imin-1], n_HArr[imin+1], temp_eqArr[imin+1])
	elif imin == 0:
		print "values around the chosen in the file:[%.2e, %.1f]" % (n_HArr[imin+1], temp_eqArr[imin+1])
	else:
		print "values around the chosen in the file:[%.2e, %.1f]" % (n_HArr[imin-1], temp_eqArr[imin-1])

	xx = mu*cst.m_p/cst.k_B/2.
	#den = delta_v/2.*(cst.Omg_m/omega_m_z(zcol))
	den = delta_v/2.*cst.Omg_m
	print "1/Om_z=%.2f, Om_z=%.2f" % (1./omega_m_z(zcol), omega_m_z(zcol))

	mvir = temp_eq**3/(den*((1. + zcol)*xx)**3)
	mvir = mvir**(0.5)/(cst.G_newton*cst.Hubble0*cst.h_hubble)

	return mvir/cst.M_sun#*cst.h_hubble


def hoeft_characteristic_mass(zcol):
	"""
	Hoeft et al 2006 eq. 6. 
	Their overdensity definition is different from Bryan & Norman original definition
	"""
	func = lambda zred: 0.73*(1. + zred)**0.18*N.exp(-(0.25*zred)**2.1)
	tauz = func(zcol)
	overden0 = virial_overdensity_z(0.)/omega_m_z(0.)
	overdenz = virial_overdensity_z(zcol)/omega_m_z(zcol)

	mc = (tauz/(1. + zcol))**1.5*(overden0/overdenz)**0.5
	return mc*1.e10/cst.h_hubble


#### expected maximum Jeans mass
def turn_around_adiabat(zcol, n_HArr):
	#zturn = 1.68/1.06*(1. + zcol) - 1.
	zturn = 2.**(2./3.)*(1. + zcol) - 1. ## cooray&sheth p14
	#den_turn = 5.55*mean_hydrogen_numden_z(zcol)
	den_turn = 4.55*mean_hydrogen_numden_z(zturn)
	temp_turn = 1.e4
	tempArr = calc_adiabatic_heating_temp(n_HArr, temp_turn, n_H_i=den_turn, gamma=5./3)
	return tempArr
	
def estimate_max_jeans_mass(zcol, Gamma12=1.):
	import plot_Td_analytic as pTda
	n_HArr, tempCCarr = pTda.get_analytic_curve_cc([], Gamma12)
	tempADarr = turn_around_adiabat(zcol, n_HArr)
	if len(tempADarr.shape) != len(tempCCarr.shape):
		print "WARNING: tempADarr and tempCCarr shapes are not consistent -- one of them will be changed"
		if len(tempADarr.shape) > len(tempCCarr.shape):
			tempADarr = tempADarr.reshape(len(tempADarr))
		else:
			tempCCarr = tempCCarr.reshape(len(tempCCarr))
	diff = N.abs(tempADarr - tempCCarr)
	print tempADarr.shape, tempCCarr.shape, n_HArr.shape, diff.shape
	ieq = N.where(diff == N.min(diff))[0][0]
	#print diff[ieq-5:ieq+5]
	mJ = jeans_mass(tempADarr[ieq], n_HArr[ieq]) ## in solar mass
	return mJ, n_HArr[ieq], tempADarr[ieq]
