"""
07/30 group array getting from the file is changed
"""

import os, sys
import numpy as N
from matplotlib.font_manager import FontProperties
sys.path.append(os.path.join("..", "analytic"))
import constants as cst
reload(cst)
import calc_Td_func as aclc
reload(aclc)


"""
def get_ptcls_accreted_at_zcol(denArr, zredArr_all, simI, deltaacc=200.):
def find_too_dense_ptcls(denArr, denlim=0.03):
def pick_ptcls_to_plot(nptcls, nptcls_toplot, method=0):
def find_maximum_curve(denArr, tempArr, templim=[1.e4, 2.e5], denlim=0.05, tempdiff=0.12):
def find_zred_accretion(denArr, zredArr_all, simI, denacc=200.):
def group_by_prop(pArr1_o, ngrp=10, fixed_interval=True, takelog=False, bins=[], pArr_ref=1., prec=1.e-2):
def group_by_pos(posArr1, ngrid=100, boxsize=10000.):
def get_high_density_ptcls(over_density, denArr_atzcol, zcol):
def calc_center_of_mass(posArr, grpidArr):
def get_ptcls_inside_virial_radius(posArr1, grpInfoArr, zcol, boxsize, virial_radius, delta=180., bkgd=True, unitconv=1.e3):
def get_virial_radius(mhalo, delta, bkgd=True):
def diff_periodic(arr, boxsize=10.e3):
"""

def convert_delta_to_nH(delta, hubble, z):
	## cst.rho_cr is calculated using cst.h_hubble so divide it by cs.h_hubble**2 and multiply the simulation's hubble**2 gives the correct rho_cr of the simulation
	
	rho_mean = aclc.mean_hydrogen_numden_z(zcol)
	return rho_mean*delta/cst.m_p


#### position handling
def diff_periodic(arr, boxsize=10.e3):
	"""
	arr is a position arr. n of columns is corresponding to the dimension to be dealed
	"""
	for icol in xrange(arr.shape[1]):
		ii = N.array([0])
		while len(ii) > 0:
			ii = N.where(arr[:, icol] < -0.5*boxsize)[0]
			arr[ii, icol] += boxsize
		ii = N.array([0])
		while len(ii) > 0:
			ii = N.where(arr[:, icol] >= 0.5*boxsize)[0]
			arr[ii, icol] -= boxsize
	return arr


def get_virial_radius(mhalo, delta, zred, bkgd=True):
	"""
	mhalo in the unit of M_sun/h, virial radius is returned in the unit of Mpc/h
	the above condition doesn't hold anymore
	"""
	if bkgd:
		rho = delta*cst.rho_cr_practicalUnit*cst.Omg_m
	else:
		rho = delta*cst.rho_cr_practicalUnit*(cst.Omg_m*(1. + zred)**3 + (1. - cst.Omg_m))/(1. + zred)**3
		#rho = delta*cst.rho_cr_practicalUnit

	rad = (3.*mhalo/rho/4./cst.pi)**(1./3.)
	return rad
	

def get_ptcls_inside_virial_radius(posArr1, grpInfoArr, zcol, boxsize, virial_radius, delta=180., bkgd=True, unitconv=1.e3):
#	if zcol 
	posdff = diff_periodic(posArr1 - grpInfoArr[:, 3:], boxsize=boxsize)
	dst = N.sqrt(N.sum(posdff**2, axis=1))

	if virial_radius: ## using a fixed radius
		ids_in_vir = N.where(dst <= virial_radius)[0]

	else:
		## grpmassArr was saved in the unit of M_sun, so multiply hubble parameter
		uniq_grps, uniq_grps_id = N.unique(grpInfoArr[:, 1], return_index=True)
		ids_in_vir = N.array([], dtype="int")
		for i, igrp in enumerate(uniq_grps):
			ii = N.where(grpInfoArr[:, 1] == igrp)[0]
			mgrp = grpInfoArr[ii[0], 2]
			vrad = get_virial_radius(mgrp, delta, bkgd=bkgd)*unitconv ## b/c simulation data is in kpc
			jj = N.where(dst[ii] <= vrad)[0]
			ids_in_vir = N.append(ids_in_vir, ii[jj])
			## to check
			if i % (len(uniq_grps)/10) == 0:
				print "m=%.2e, vr=%.1f, %d, %d" % (mgrp, vrad, len(ii), len(jj))
	
	return ids_in_vir



def calc_center_of_mass(posArr, grpidArr):
	"""
	pos information and the indicator of groups
	periodic condition is applied already to posArr
	"""
	print "calculating center of mass"
	uniq_grps = N.unique(grpidArr)
	com_gasArr = N.zeros((len(uniq_grps), 3))
	nptclsArr = N.zeros(len(uniq_grps), dtype="int")
	for i, igrp in enumerate(uniq_grps):
		ii = N.where(grpidArr == igrp)[0]
		com = N.sum(posArr[ii, :], axis=0)/len(ii)
		com_gasArr[i, :] = com
		nptclsArr[i] = len(ii)
	return com_gasArr, nptclsArr




################### to pick particles to plot #####################

def group_by_pos(posArr1, ngrid=100, boxsize=10000.):
	dgrid = boxsize/ngrid
	boxid = N.array(posArr1/dgrid, dtype="int")
	boxid = boxid[:, 0]*ngrid**2 + boxid[:, 1]*ngrid + boxid[:, 2]
	return boxid
	

def group_by_prop(pArr1_o, ngrp=10, fixed_interval=True, takelog=False, bins=[], pArr_ref=1., prec=1.e-2):
	"""
	binning by properties and 
	returning the group ids corresponding to the given properties
	if pArr_ref, which is only valid when bins are defined, is an array, bining by comparing two different properties
	i.e. (bins[i] * pArr_ref[j] < pArr1_o [j] < bins[i+1] * pArr_ref[j]
	"""
	if takelog:
		pArr1 = N.log10(pArr1_o)
		if len(bins) > 0: bins = N.log10(bins)
	else:
		pArr1 = pArr1_o.copy()

	if len(bins) > 0:
		## note that if p_id is negative, that is invalid index
		p_id = -1*N.ones_like(pArr1, dtype="int")
		for ibin in xrange(len(bins)-1):
			ii = N.where(N.all([pArr1 >= bins[ibin]*pArr_ref, pArr1 < bins[ibin+1]*pArr_ref], axis=0))[0]
			p_id[ii] = ibin
		print "among %d ptcls, %d are in the given bin range" % (len(pArr1), len(N.where(p_id != -1)[0])) 
	else:
		if fixed_interval:
			## ngrp can be varied
			p_range = [(1. - prec)*N.min(pArr1), (1. + prec)*N.max(pArr1)]
			dp = (p_range[1] - p_range[0])/ngrp
			dp = prec*N.round(dp/prec)
			p_id = N.array((pArr1 - p_range[0])/dp, dtype="int")
			bins = p_range[0] + dp*N.arange(ngrp + 1)
			if takelog:
				print "range: %.2e, %.2e, ngrp: %d, dp: %.2e" % (10.**bins[0], 10.**bins[-1], len(bins)-1, 10.**dp)
			else:
				print "range: %.2e, %.2e, ngrp: %d, dp: %.2e" % (bins[0], bins[-1], len(bins)-1, dp)
		else:
			## each bin contains the same number of particles
			isorted = N.argsort(pArr1)
			p_id = -1*N.ones_like(pArr1, dtype="int")
			bins = N.zeros(ngrp + 1)
			dn = len(pArr1)/ngrp ## n of elements in each bin
			for i in xrange(ngrp - 1):
				idgrp = isorted[i*dn:(i+1)*dn]
				p_id[idgrp] = i
				bins[i] = pArr1[isorted[i*dn]]
				print "min %.2e, max %.2e" % (N.min(pArr1[idgrp]), N.max(pArr1[idgrp]))
			## the last group may have more than dn if len(pArr1)/ngrp has a remaining
			idgrp = isorted[(ngrp-1)*dn:]
			p_id[idgrp] = ngrp-1
			bins[-2] = pArr1[isorted[(ngrp-1)*dn]]
			bins[-1] = N.max(pArr1)*1.000001
			if N.min(pArr1) < 0.:
				bins[0] = N.min(pArr1)*1.000001
			else:
				bins[0] = N.min(pArr1)*0.000009
	return p_id, bins


def get_ptclIDs_to_plot(x, y, param):
	ids, nmax, f_mul, discard = param[0], param[1], param[2], param[3]
	if len(ids) > nmax:
		dd = len(ids)/nmax; dd0 = int(dd*f_mul) ## dd0 is starting point
		iptcltoplot = N.arange(dd0, len(ids), dd)
		iptcltoplot = iptcltoplot[:nmax]

		if discard:
			print "before:", iptcltoplot
			for i, ip in enumerate(iptcltoplot):
				ipp = ip
				nzeros = x.shape[1]
				ndense = x.shape[1]
				while (nzeros > int(x.shape[1]*0.2) or ndense > int(x.shape[1]*0.7) or (ndense + nzeros) > int(x.shape[1]*0.9)) and (ipp - ip) < dd:
					zeros = N.where(N.any([x[ipp, :] == 0., y[ipp, :] == 0.], axis=0))[0]
					nzeros = len(zeros)
					denses = N.where(x[ipp, :] > 0.01)[0]
					ndense = len(denses)
					if (len(zeros) > int(x.shape[1]*0.2) or ndense > int(x.shape[1]*0.7)):
						print "discarding: %d b/c:" % ipp, x[ipp, :]
						ipp += 1
				iptcltoplot[i] = ipp 
			print "after:", iptcltoplot
	else:
		dd = 1; dd0 = 0
		iptcltoplot = N.arange(len(ids))
	return iptcltoplot


def get_non_zeros_only(denArr, tempArr):        
	i0 = N.where(N.all(N.column_stack((denArr, tempArr)) != 0., axis=1))[0]
	## test print
	print "%d will be kept among %d since density and temperature at any z's aren't 0" % (len(i0), nptcls)
	return i0



############################

def get_high_density_ptcls(over_density, denArr_atzcol, zcol):
	"""
	when choosing high density particles outside of get_data function
	"""
	#izcol = N.where(zredArr == zcol)[0]
	rho_mean = aclc.mean_hydrogen_numden_z(zcol)
	print "rho_mean: %.2e at z=%.2f" % (rho_mean, zcol)
	i0 = N.where(denArr_atzcol >= over_density*rho_mean)[0]
	print "%d are over the %.1f x mean density" % (len(i0), over_density)
	return i0



def find_too_dense_ptcls(denArr, zredArr, denlim=0.03, zlim=8.):
	"""
	too dense particles are numerical effect. discard them to find mJmax or to plot
	this function returns the row and the starting column of the particle. i.e. denArr[ip, N.min(jj):] should be neglected
	ilim is the limit of the redshift because at high-z, mean density is high and we don't want to include it.
	"""
	ilim = N.where(zredArr < zlim)[0][0]
	print "limit redshift is %.1f and the corresponding column is %d" % (zlim, ilim)
	i_denseptcls = N.where(N.any(denArr[:, ilim:] >= denlim, axis=1))[0]
	print "%d are too dense" % len(i_denseptcls)

	if len(i_denseptcls) > 0:
		jcols = N.zeros_like(i_denseptcls, dtype="int")
		for i, ip in enumerate(i_denseptcls):
			den = denArr[ip, ilim:]
			jj = N.where(den >= denlim)[0]
			jcols[i] = N.min(jj) + ilim
			#print "%d %d" % (ip, jcols[i])
			#print "{:}".format(denArr[ip, jcols[i]:])
		if len(i_denseptcls) < 20:
			print "dense particles are {:}".format(i_denseptcls)
			print "the column corresponding to redshifts when it exceeds the density limit are {:}".format(jcols)
			print "the densities are: {:}".format(denArr[i_denseptcls, jcols])
		return N.column_stack((i_denseptcls, jcols))	
	else:
		return []



def find_maximum_curve(denArr, tempArr, templim=[1.e4, 2.e5], denlim=0.05, tempdiff=0.12):
	"""
	find the maximum of the trajectory on density-temperature
	"""
	#tt = N.log10(tempArr[ivalid, :])
	#dd = N.log10(denArr[ivalid, :])

	nptcls = len(denArr)
	locslope = (tempArr[:, 1:] - tempArr[:, :-1])/(denArr[:, 1:] - denArr[:, :-1] + 1.e-10)
	slope_sign = N.sign(locslope)
	signchange = ((slope_sign[:, 1:] - slope_sign[:, :-1]) != 0).astype(int)

	imaxima = -1*N.ones(nptcls, dtype='int')

	for iptcl in xrange(nptcls):
		## prevent from using the non-valid case
		ics = N.where(signchange[iptcl, :])[0]
		ics += 1
		if len(ics) > 0:
			imaxcand = []
			for ic in ics:
				if locslope[iptcl, ic-1] > 0 and tempArr[iptcl, ic] >= templim[0] and tempArr[iptcl, ic] <= templim[1] and denArr[iptcl, ic] < denlim:
					imaxcand.append(ic)
			imaxcand = N.array(imaxcand)
			if len(imaxcand) > 1:
				#it = N.where(tempArr[iptcl, imaxcand] == N.max(tempArr[iptcl, imaxcand]))[0]
				it = N.where(N.abs(tempArr[iptcl, imaxcand]/N.max(tempArr[iptcl, imaxcand])-1.) < tempdiff)[0]
				if len(it) > 1:
					#print it, tempArr[iptcl, imaxcand[it]]
					mjeans_cand = aclc.jeans_mass(tempArr[iptcl, imaxcand[it]], denArr[iptcl, imaxcand[it]])
					imj = N.where(mjeans_cand == N.max(mjeans_cand))[0]
					imaxima[iptcl] = imaxcand[it[imj]]
				else:   
					imaxima[iptcl] = imaxcand[it]
				#imaxima[i] = min(imaxcand)
			elif len(imaxcand) == 1:
				imaxima[iptcl] = imaxcand[0]
			else:
				#print "given densities" 
				#print denArr[iptcl, :]
				#print "and temperatures are"
				#print tempArr[iptcl, :]
				"""
				ir = N.where(N.all((denArr[iptcl, :] > denArr[iptcl, 1], denArr[iptcl, :] < denlim), axis=0))[0]
				if len(ir) == 0:
					print "there is no gas particles which move to denser than the density at reionization"
					#print denArr[iptcl, :]
					imaxima[iptcl] = 1
				else:
					mjeans_cand = aclc.jeans_mass(tempArr[iptcl, ir], denArr[iptcl, ir])
					it = N.where(mjeans_cand == N.max(mjeans_cand))[0]
					imaxima[iptcl] = ir[it]
				"""
				mjeans_cand = aclc.jeans_mass(tempArr[iptcl, :], denArr[iptcl, :])
				it = N.where(mjeans_cand == N.max(mjeans_cand))[0]
				imaxima[iptcl] = it
				print "no peak candidate found for %d: assigned jeans mass: %.4e" % (iptcl, aclc.jeans_mass(tempArr[iptcl, imaxima[iptcl]], denArr[iptcl, imaxima[iptcl]]))

	return imaxima


def find_zred_accretion(denArr, zredArr_all, simI, deltaacc=200.):
    ## denArr should be just from the simulation, not coverted to n_H rho/rho_bg
    ## check
    if simI['convert_density']:
 		print 'if noneed to convert the overdensity to n_H, check simI and use the raw outputs not converting to n_H ------> converting to delta'
		rho_mean = aclc.mean_hydrogen_numden_z(simI['zcol'])
		deltaArr = denArr*cst.m_p/rho_mean

    zaccArr = -1*N.ones(len(denArr))

    for i, delta in enumerate(deltaArr):
        iz = N.where(delta >= deltaacc)[0]
        if len(iz) > 0:
            zaccArr[i] = zredArr_all[iz[0]]

    return zaccArr


def get_ptcls_accreted_at_zcol(denArr, zredArr_all, simI, deltaacc=200.):
    ## denArr should be just from the simulation, not coverted to n_H rho/rho_bg
    ## check
	
	zaccArr = find_zred_accretion(denArr, zredArr_all, simI, deltaacc=deltaacc)
	iacc = N.where(zaccArr == simI['zcol'])[0]
	print "%d are accreted at zcol among %d" % (len(iacc), len(denArr))
	return iacc



def find_fitting_params_cooling_curve(ioniz_param):
	from plot_Td_analytic import get_analytic_curve_cc 

	Gamma12 = ioniz_param[0]; inc_HeII = ioniz_param[1] ## incHeII = false is default

	n_HArr, tempCArr = get_analytic_curve_cc([], N.array([Gamma12]), inc_HeII=inc_HeII)
	tempCArr = tempCArr.reshape(len(tempCArr)) ## tempCArr returns (len(tempCArr), 1)
	x = N.log10(n_HArr)
	y = N.log10(tempCArr)
	ii = N.where(y == y[0])[0] ## to discard the values artificially set
	x = x[ii[-1]+1:]; y = y[ii[-1]+1:]
	coef = N.polyfit(x, y, 2)
	return coef


def get_max_jeans_mass_with_cooling(denArr, tempArr, ioniz_param):
	coef = find_fitting_params_cooling_curve(ioniz_param)
	print coef
	func = N.poly1d(coef)
	tempfit = func(N.log10(denArr + 1.e-10))
	ndata = len(denArr)
	mjmaxArr = N.zeros(ndata); imaxArr = N.zeros_like(mjmaxArr, dtype='int')

	for i in xrange(ndata):
		## fitting function works only for the case that density is larger than 1.e-5
		ikeep = N.where(N.all([tempArr[i, :] <= 10.**tempfit[i, :] + 5.e3, tempArr[i, :] > 1.e-10], axis=0))[0] ## 5.e3 is added to allow some errors 
		#ikeep = N.where(N.all([N.log10(tempArr[i, :] + 1.e-10) <= tempfit[i, :] + N.log10(5.e3), denArr[i, :] > 1.e-5], axis=0))[0] ## 5.e3 is added to allow some errors 
		#print N.column_stack((N.log10(tempArr[i, :] + 1.e-10), tempfit[i, :] + N.log10(5.e3)))
		if len(ikeep) < 1:
			print "nothing to keep for %dth row" % i
			print N.column_stack((tempArr[i, :], denArr[i, :], 10.**tempfit[i, :]))
			mjmaxArr[i] = 0.
		else:
			mj = aclc.jeans_mass(tempArr[i, ikeep] + 1.e-10, denArr[i, ikeep] + 1.e-10)
			mjmaxArr[i] = N.max(mj)
			imax = N.where(mj == mjmaxArr[i])[0]
			if len(imax) > 1:
				print "there are many maximums for {:d}th particle: {:}".format(i, mj)
				print N.column_stack((tempArr[i, :], denArr[i, :]))
			imaxArr[i] = ikeep[imax]
		## test
		"""
		if len(ikeep) < len(denArr[0, :]):
			mjtest = aclc.jeans_mass(tempArr[i, :] + 1.e-10, denArr[i, :] + 1.e-10)	
			mjtestmax = N.max(mjtest)
			imaxtest = N.where(mjtest == mjtestmax)[0]
			#print "%5d - %6d - %8.2e - %5d (%8.2e - %5d)" % (i, len(ikeep), mjmaxArr[i], imaxArr[i], mjtestmax, imaxtest)
			#print denArr[i, :], tempArr[i, :], tempfit[i, :]
		mjtest = aclc.jeans_mass(tempArr[i, :] + 1.e-10, denArr[i, :] + 1.e-10)	
		print "%s %d: %.2e (%d) %s" % ("*"*5, i, mjmaxArr[i], imaxArr[i], "*"*5)
		print mjtest
		"""
	return mjmaxArr, imaxArr


def get_max_jeans_mass_no_cooling(denArr, tempArr):
	mjArr = aclc.jeans_mass(tempArr + 1.e-10, denArr + 1.e-10)
	mjmaxArr = N.max(mjArr, axis=1); imaxArr = N.zeros_like(mjmaxArr, dtype='int')
	
	for i in len(mjArr):
		imaxArr[i] = N.where(mjArr[i, :] == mjmaxArr)[0]

	return mjmaxArr, imaxArr
		


def get_peak_jeans_mass():
	print "peak"	
