"""
07/30 group array getting from the file is changed
"""

import os, sys
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
from matplotlib.font_manager import FontProperties
import constants as cst
reload(cst)
import calc_Td_func as aclc
reload(aclc)


def plot_analytic_curve_altogether(ax, n_HArr, tempArrs, pltparam):

	tempJArr, tempAArr, tempCArr, tempSArr = tempArrs[0], tempArrs[1], tempArrs[2], tempArrs[3]

	xr, yr, write_text, mJ_init, ionizRL = pltparam[0], pltparam[1], pltparam[2], pltparam[3], pltparam[4]
	"""
	plot all together
	order in write_text is jeans, adiabatic, cooling, self-shielding
	mJ_ini contains the jeans mass used to get tempJArr
	""" 
	
	#colorArr = []
	#for icl in xrange(4):
	#	colorArr.append(plt.get_cmap("cool")(float(icl)/3))
	#colorArr[0] = "yellowgreen"
	colorArr = ['#99CC00', '#FF6699', '#8585FF', '#FFA347'] ## green-ish, pink, blue, orange, 
	lineArr = ['--', '-.', ':', '--']
	pltkwarg = dict(marker = "None", linewidth = 0.02)

	pltkwarg["color"] = colorArr[0]; pltkwarg["linestyle"] = lineArr[0]
	pltkwarg["dashes"] = [10, 10]
	if len(tempJArr) > 0:	
		if write_text[0]: txList = make_textList_in_power(mJ_ini)
		else: txList = []
		plot_analytic_curve_loop(ax, n_HArr, tempJArr, pltkwarg, txList, xr, yr, [0.05, 0.05])

	pltkwarg["color"] = colorArr[1]; pltkwarg["linestyle"] = lineArr[1]
	pltkwarg["dashes"] = [3, 3, 1, 3]
	if len(tempAArr) > 0:	
		plot_analytic_curve_loop(ax, n_HArr, tempAArr, pltkwarg, [], xr, yr, [0.05, 0.05]) 
	pltkwarg["color"] = colorArr[2]; pltkwarg["linestyle"] = lineArr[2]
	pltkwarg["dashes"] = [2, 1]
	if len(tempCArr) > 0:
		if write_text[2]: txList = ["%d" % ionizR for ionizR in ionizRList]
		else: txList = []
		plot_analytic_curve_loop(ax, n_HArr, tempCArr, pltkwarg, txList, xr, yr, [-0.05, -0.05]) 

	pltkwarg["color"] = colorArr[3]; pltkwarg["linestyle"] = lineArr[3]
	pltkwarg["dashes"] = [5, 3]
	if len(tempSArr) > 0:
		if write_text[3]: txList = ["%d" % ionizR for ionizR in ionizRList]
		else: txList = []
		plot_analytic_curve_loop(ax, n_HArr, tempSArr, pltkwarg, txList, xr, yr, [-0.05, -0.05]) 


def get_text_loc_axis_coord(xArr, yArr, xr, yr, delta):
	"""
	when calculating location, take absolute value for the possibility to locate the text on the bottom or left
	len(delta) should be 2. It sets the location of the text to be inside/outside of the figure 
	"""
	dff = N.abs(yArr - yr[1])
	ii = N.where(dff == N.min(dff))[0]
	if xArr[ii] < xr[1]:
		yloc = 1. + delta[1]
		xloc = N.abs(N.log10(xArr[ii]/xr[0]))/N.abs(N.log10(xr[1]/xr[0]))
	else:
		xloc = 1. + delta[0]
		dff = N.abs(xArr - xr[1])
		jj = N.where(dff == N.min(dff))[0]
		yloc = N.abs(N.log10(yArr[jj]/yr[0]))/N.abs(N.log10(yr[1]/yr[0]))
	return [xloc, yloc]


def plot_analytic_curve_loop(ax, n_HArr, tempArr, pltkwarg, txList, xr, yr, delta, fs=8):
	"""
	tempArr should be in 2D shape
	"""
	
	for icol in xrange(tempArr.shape[1]):
		ax.plot(n_HArr, tempArr[:, icol], **pltkwarg)

		if len(txList) > 0:
			txloc = get_text_loc_axis_coord(n_HArr, tempArr[:, icol], xr, yr, delta)
			ax.text(txloc[0], txloc[1], txList[icol], fontsize=fs, horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)

		

def plot_analytic_curve(ax, n_HArr, tempArr, pltkwarg, write_text, tx, txloc=[], fs=8):
	"""
	used in get_analytic_curve_** or it can be used independently
	write_text: den_coord, temp_coord, axis_coord, or False
	if write_text is den_coord, temp_coord, txloc is given by density or temperature, respectively, if axis coord, txloc should be 1x2 list
	"""
	ax.plot(n_HArr, tempArr, **pltkwarg)

	if write_text == "temp_coord":
		diff = N.abs(tempArr - txloc)
		ii = N.where(diff == N.min(diff))[0]
		ax.text(n_HArr[ii], tempArr[ii], r"$%s$" % tx, fontsize=fs)
	elif write_text == "den_coord":
		diff = N.abs(n_HArr - txloc)
		ii = N.where(diff == N.min(diff))[0]
		ax.text(n_HArr[ii], tempArr[ii], r"$%s$" % tx, fontsize=fs)
	elif write_text == "axis_coord":
			ax.text(txloc[0], txloc[1], r"$%s$" % tx, fontsize=fs, horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
	else:
		pass

def check_binsize(bins):
	binrat = bins[1:] / bins[:-1]
	ii = N.where(binrat < 2.)[0]
	if len(ii) > 0:
		ii += 1
		bins = N.delete(bins, ii)
	return bins


def make_textList_in_power(valArr):
	pwArr = N.round(N.log10(valArr))
	mulArr = 0.1*N.round(valArr/10.**(pwArr - 1.))
	txList = []
	for mul, pw in zip(mulArr, pwArr):
		if mul < 1.:
			mul *= 10.; pw -= 1
		txList.append(r"$%.1f\times 10^{%d}$" % (mul, pw))
	return txList


def get_analytic_curve_mJ(n_HArr, mJ_ini):
	## jeans mass temp
	npnts = len(n_HArr)

	if type(mJ_ini) == float: ## if it is only one number
		m_JArr = N.array([mJ_ini]) #N.array([mJ_ini/3., mJ_ini, mJ_ini*3., mJ_ini*10.])
	else:
		## if a list for the jeans mass is given
		if type(mJ_ini) != N.ndarray:
			mJ_ini = N.array(mJ_ini)
		m_JArr = N.sort(mJ_ini)

	temp_jeansArr = N.zeros((npnts, len(m_JArr)))
	for i, m_J in enumerate(m_JArr):
		temp_jeansArr[:, i] = aclc.calc_jeans_mass_temp(n_HArr, m_J, gamma=5./3.)
		
	#return m_JArr, temp_jeansArr
	return temp_jeansArr
	

def get_analytic_curve_ad(n_HArr, vals_i):

	den_ini, temp_ini = vals_i[0], vals_i[1]
	"""
	temp_ini can be a data array, a list or an array of the initial temperature or a number 
	"""
	npnts = len(n_HArr)
	if type(temp_ini) == N.ndarray:
		temp_iArr = temp_ini; den_iArr = den_ini
	elif type(temp_ini) == list: # list
		temp_iArr = N.array(temp_ini); den_iArr = N.array(den_ini)
	else:
		temp_iArr = N.array([temp_ini]); den_iArr = N.array([den_ini])

	#print temp_iArr, den_iArr
	temp_adArr = N.zeros((len(n_HArr), len(temp_iArr)))
	for i, temp_i in enumerate(temp_iArr):
		temp_adArr[:, i] = aclc.calc_adiabatic_heating_temp(n_HArr, temp_i, n_H_i=den_iArr[i])
	
	return temp_adArr


def get_analytic_curve_cc_old(n_HArr, Gamma12Arr):

	if type(Gamma12Arr) == float:
		Gamma12Arr = N.array([Gamma12Arr])
	elif type(Gamma12Arr) == list:
		Gamma12Arr = N.array(Gamma12Arr)

	## if Gamma12 is a number or a list of Gamma12_H = Gamma12_HeI
	if len(Gamma12Arr.shape) == 1:
		Gamma12sArr = N.zeros((len(Gamma12Arr), 3))
		for i, g0 in enumerate(Gamma12Arr):
			Gamma12sArr[i, :2].fill(g0)
	else:
		Gamma12sArr = Gamma12Arr.copy()
		
	temp_coolArr = N.zeros((len(n_HArr), len(Gamma12sArr)))
	for i, Gamma12s in enumerate(Gamma12sArr):
		temp_coolArr[:, i] = aclc.calc_coolingtime_eqto_fftime_temp(n_HArr, Gamma12s, alpha=0., recomb_cooling_case="A", temp0=1.e4)

	return temp_coolArr


def get_analytic_curve_cc(n_HArr, Gamma12Arr, inc_HeII=False):
	"""
	use the numbers Matt generated: n_HArr is not required
	"""
	abspath_current = os.path.abspath(".")
	ii = abspath_current.find("highz_dwarfs")
	abspath_current = abspath_current[:ii+13]
	path = os.path.join(abspath_current, "analytic", "curves")


	if type(Gamma12Arr) == float: Gamma12Arr = N.array([Gamma12Arr])
	elif type(Gamma12Arr) == list: Gamma12Arr = N.array(Gamma12Arr)
	else: pass

	for i, Gamma12 in enumerate(Gamma12Arr):
		filename = "Teq_G%.1f_full.dat" % (Gamma12)
		n_HArr1, temp_coolArr1 = N.loadtxt(os.path.join(path, filename), usecols=(0, 1), unpack=True)
		if i == 0:
			temp_coolArr = N.zeros((len(n_HArr1), len(Gamma12Arr)))
			n_HArr = n_HArr1.copy()
		temp_coolArr[:, i] = temp_coolArr1
		## check read n_HArr
		idff = N.where(n_HArr1 != n_HArr)[0]
		if len(idff) > 0:
			print "WARNING"
			print "read n_HArr is different from the previous one - don't use the previous n_HArr when plotting"

	if inc_HeII: ## only works now for G1.0
		filename = "Teq_G1.0_full_wHeII.dat" 
		n_HArr1, temp_coolArr1 = N.loadtxt(os.path.join(path, filename), usecols=(0, 1), unpack=True)
		temp_coolArr = N.column_stack((temp_coolArr, temp_coolArr1))
		Gamma12Arr = N.append(Gamma12Arr, 1.0)
		## check read n_HArr
		idff = N.where(n_HArr1 != n_HArr)[0]
		if len(idff) > 0:
			print "WARNING: including HeII"
			print "read n_HArr is different from the previous one - don't use the previous n_HArr when plotting"
	
	return n_HArr, temp_coolArr


def get_analytic_curve_ss(n_HArr, Gamma12Arr, sigma=None):
	## self-shielding
	## Gamma12 can be a number

	if type(Gamma12Arr) == float:
		Gamma12Arr = N.array([Gamma12Arr])
	elif type(Gamma12Arr) == list:
		Gamma12Arr = N.array(Gamma12Arr)

	temp_ssArr = N.zeros((len(n_HArr), len(Gamma12Arr)))
	for i, Gamma12 in enumerate(Gamma12Arr):
		if sigma != None:
			temp_ssArr[:, i] = aclc.calc_self_shielding_temp(n_HArr, Gamma12, sigma=sigma)
		else:
			temp_ssArr[:, i] = aclc.calc_self_shielding_temp(n_HArr, Gamma12)

	return temp_ssArr


def get_Tini(denArr, tempArr, n, denr0, tempr0):
	den0 = [N.min(denArr), N.max(denArr)]
	temp0 = [N.min(tempArr), N.max(tempArr)]
	if len(denr0.flatten()) > 2:
		denr = [N.min(denr0[:, 0]), N.max(denr0[:, 0])]
	else: denr = denr0
	if len(tempr0.flatten()) > 2:
		tempr = [N.min(tempr0[:, 0]), N.max(tempr0[:, 0])]
	else: tempr = tempr0

	for i in xrange(2):
		den0[i] = [denr[i], den0[i]][denr[i] < den0[i]]
		temp0[i] = [tempr[i], temp0[i]][tempr[i] < temp0[i]]

	den0 = N.log10(den0); temp0 = N.log10(temp0)
	den_ini = N.logspace(den0[0], den0[1], n)	
	temp_ini = N.logspace(temp0[0], temp0[1], n)	

	return den_ini, temp_ini
