
import os, sys
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
from matplotlib.backends.backend_pdf import PdfPages


sys.path.append("..")
import plot_Td_func as pTdf
import read_outputs as RO
import plot_Td_analytic as pTda
import calc_Td_func as cTdf
import construct_dataArr as cdA


def bin_by_grpmass(massArr, bin_param, nptcl_each):
	## if I want to use specific binning, then set mrange_x as bins and nbins == 0
	## this can be used any times to make equal number particles in a given bin
	mrange_x = bin_param[0]; nbins = bin_param[1]

	if nbins == 0:
		bins = mrange_x; nbins = len(bins) - 1
	else:
		bins = N.logspace(N.log10(mrange_x[0]*0.99), N.log10(mrange_x[1]*1.01), nbins+1)

	print "in bin_by_grpmass{:}".format(bins)
	p_id, bins = pTdf.group_by_prop(massArr, takelog=False, bins=bins)

	nptcl_each_arr = N.zeros(nbins)
	for i in xrange(nbins):
		nptcl_each_arr[i] = len(N.where(p_id == i)[0])
	print "nptcl_each_arr{:}".format(nptcl_each_arr)

	print "%d ptcls will be used to find out which ptcls are in the given bins" % N.sum(nptcl_each_arr)

	## how many particles I want to have in each bin
	if nptcl_each == 0:
		jj = N.where(nptcl_each_arr > 100)[0]
		#jj = N.where(nptcl_each_arr != 0)[0]
		nptcl_each = N.min(nptcl_each_arr[jj])
	#else:
	#   jj = N.where(nptcl_each_arr > 100)[0]
	#   #jj = N.where(nptcl_each_arr != 0)[0] 
	#   nptcl_each = N.max([N.min(nptcl_each_arr[jj]), nptcl_each])
	print "each bin will contain %d ptcls" % nptcl_each

	i_p_id_sorted = N.argsort(p_id)

	## for the values not lying in the given bin range, p_ids are set to -1
	## I don't want to inclue them 
	ineg = N.where(p_id[i_p_id_sorted] == -1)[0]
	if len(ineg) > 0: nn = N.max(ineg) + 1
	else: nn = 0
	nids = 0
	idArr_plt = -1*N.ones(nptcl_each*nbins, dtype="int")
	print "check idArr_plt, length: %d = nptcl_each (%d) * nbins(%d)" % (len(idArr_plt), nptcl_each, nbins)
	for i in xrange(nbins):
		if nptcl_each_arr[i] > nptcl_each:
			dn = int(N.round(nptcl_each_arr[i]/(1.*nptcl_each)))
		else:
			dn = 1

		if dn == 1: 
			if nptcl_each_arr[i] == 0:
				ieach = N.array([])
			elif nptcl_each_arr[i] < nptcl_each:
				ieach = N.arange(nn, nn + nptcl_each_arr[i])
			else:
				ieach = N.arange(nn, nn + nptcl_each)
		else:
			ieach = N.arange(nn, nn + nptcl_each_arr[i], dn)
			## correct the case that len(ieach) doesn't match with nptcl_each
			dff = len(ieach) - nptcl_each
			if dff > 0:
				ieach = ieach[:nptcl_each]
			elif dff < 0:
				ieach1 = N.arange(nn+1, nn + nptcl_each_arr[i], dn)[:-dff]
				ieach = N.hstack((ieach, ieach1))
				ieach = N.sort(ieach)

		ieach = N.array(ieach, dtype='int'); nieach = len(ieach)
		#print "ieach:{:} (n_ieach={:d}), the id starting in idArr_plt: {:d}".format(ieach, len(ieach), nids)
		#idArr_plt[i*nptcl_each:(i+1)*nptcl_each] = idx = i_p_id_sorted[ieach]
		#print len(i_p_id_sorted), max(ieach), min(ieach), len(i_p_id_sorted[ieach]), len(idArr_plt[nids:nids+nieach])
		#print nids+nieach, nids, idArr_plt[nids:nids+nieach]
		#print "idArr_plt:", idArr_plt, "length", len(idArr_plt)
		idArr_plt[nids:nids+nieach] = i_p_id_sorted[ieach]
		#idArr_plt[nids:nids+len(ieach)] = idx = i_p_id_sorted[ieach]
		nn += nptcl_each_arr[i]; nids += len(ieach)
		#print "idArr_plt:", idArr_plt, "length", len(idArr_plt)
		print "%d are found among %d in [%.2e, %.2e]" % (len(ieach), nptcl_each_arr[i], bins[i], bins[i+1])
		#if len(ieach) > 0:
			#print "\ttrue min max is [%.2e, %.2e]" % (N.min(massArr[idx]), N.max(massArr[idx]))

	inone = N.where(idArr_plt == -1)[0]
	idArr_plt = N.delete(idArr_plt, inone)

	return idArr_plt





def choose_particles_by_jeans_mass(m_GArr, m_JArr1, propArr, coeff=[4., 4.], p_bins=N.array([1.e-5, 10., 100., 100000.]), equal_n_in_mbin = True, nbins_mass=10, nptcl_each=0):

	mrange_x = [0.999*N.min(m_GArr), 1.001*N.max(m_GArr)]
	mrange_y = [0.999*N.min(m_JArr1), 1.001*N.max(m_JArr1)]

	if equal_n_in_mbin:
		idArr_plt = bin_by_grpmass(m_GArr, (mrange_x, nbins_mass), nptcl_each)
		print "particles to plot: %d" % len(idArr_plt)
		m_xArr = m_GArr[idArr_plt]; m_yArr = m_JArr1[idArr_plt]
		propArr = propArr[idArr_plt]
	else:
		m_xArr = m_GArr; m_yArr = m_JArr1


	iptcls = N.zeros(len(propArr), dtype='int')

	## low density but seems to be accreted
	p_ids = N.where(N.all([propArr >= p_bins[0], propArr < p_bins[1]], axis=0))[0]
	ii = N.where(m_yArr[p_ids] < coeff[0]*m_xArr[p_ids])[0]
	iptcls[p_ids[ii]] = -1
	print "test print"
	for i in ii:
		print "%.2e %.2e %.1f" % (4.*m_xArr[p_ids[i]], m_yArr[p_ids[i]], propArr[p_ids[i]])

	## high density but doesn't seem to be accreted
	p_ids = N.where(N.all([propArr >= p_bins[-2], propArr < p_bins[-1]], axis=0))[0]
	ii = N.where(m_yArr[p_ids] > coeff[1]*m_xArr[p_ids])[0]
	iptcls[p_ids[ii]] = 1

	if equal_n_in_mbin:
		return idArr_plt, iptcls	
	else:
		return iptcls



def get_ids_mj_vs_mh(pids, i_indicator):
	## mj < mh but the density is low (i.e. not comllapsed)
	ilow = N.where(i_indicator == -1)[0] 
	ii = pids[ilow]

	ihigh = N.where(i_indicator == 1)[0]
	jj = pids[ihigh]
	return ii, jj


def remove_earlycollapsed_ptcls_coldIGM(simI, rparam):
	## previously, this function got rid of the particles which haven't had reached to 10^4K but now it gets rid of the particles which reached the collapsing condition either being inside r_vir or overdensity is larger than 200
	simI['simbkgd'] = 'coldIGM'
	print "simI in remove earlycollapsed particles coldIGM: {:}".format(simI)
	DC = RO.DataConstruc(**simI)
	if len(rparam) == 2:
		izcol, zarr, darr, tarr, garr = DC.get_entire_array2(rparam)
	else:
		izcol, zarr, darr, tarr, garr = DC.get_entire_array1()

	ivalid = N.where(N.all([darr[:, -1] > 200., N.all(darr[:, :-1] < 200., axis=1)], axis=0))[0]
	
	return ivalid, darr, DC.convert_density(), tarr, garr
