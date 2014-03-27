"""
find local M_J and plot max(M_J) vs. M_group
this is only possible when the data files have group masses in there

this version: in weirdo plot, no f_{n_H} plot
"""
import os, sys
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
from matplotlib.backends.backend_pdf import PdfPages


sys.path.append("..")
import plot_Td_func as pTdf
import plot_Td_analytic as pTda
import read_outputs as RO
import plot_mscale_func as pmf
reload(pmf)
import constants as cst
reload(cst)
import calc_Td_func as cTdf 
reload(cTdf)
from matplotlib.font_manager import FontProperties


## plot range
#mrange_x = [1.e8, 3.e11]; mrange_y = [1.e8, 3.e11]
#reset_mrange = True; equal_range = False ## equal range is valid only if reset_mrange is True

## want to plot M_h vs. n_H(z_col) > over_density2*mean(n_H(z_col))
#over_density2 = 200.

##  we want to see what is the behavior of the particles in left above in the plot of M_J^max vs. M_h. 
## mgrp_lim sets the left limit
## xx is the coefficient of M_h,i.e. M_J^max = xx*M_h. It sets the above boundary.
#mgrp_lim = 2.5e9; xx = 4.


## to plot mJ vs mh, after running the code, do followings:
"""
zredArr_all, denArr_all, tempArr_all, posArr_all, grpInfoArr_all, ivalid, izcol = get_data(simI)
nH_mean, collcond, data_abspath, subpath, masscArr, masscnameArr = get_some_values(**simI)
fignameh = "mJmax" # "mJmax"; z_mJ = zcol; "mJ%.1f" % z_mJ
fignamestr = [fignameh, collcond, subpath]
m_GArr = grpInfoArr_all[ivalid, 2]
m_JArr_all, m_JArr1 = get_jeans_mass_arr(fignameh, denArr_all[ivalid, :], tempArr_all[ivalid, :], zredArr_all)
propArr = denArr_all[ivalid, izcol]/nH_mean
#pr = [N.min(propArr), N.max(propArr)]
pr = N.log10(N.array([1., 400.]))
izcoll, idarr = choose_ptcls_by_collapsed_redshift(denArr_all[ivalid, :], zredArr_all, overdensity=200, inc_not_collapsed=True, **simI); ipltptcls=1 ## ipltptcls: collapsed any redshifts:0, not yet collpased until zcol: 1, collpased at zcol: 2
make_plot_mJ_mh(m_JArr1[idarr], m_GArr[idarr], nH_mean, fignamestr, masscArr=masscArr, masscnameArr=masscnameArr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, propArr=propArr[idarr], takelog_prop=True, mrange_x=mrange_x, mrange_y=mrange_y, ipltptcls=ipltptcls, pr=pr)
izcoll, idarr = choose_ptcls_by_collapsed_redshift(denArr_all[ivalid, :], zredArr_all, overdensity=200, inc_not_collapsed=False, **simI); ipltptcls=2 ## ipltptcls: collapsed any redshifts:0, not yet collpased until zcol: 1, collpased at zcol: 2
make_plot_mJ_mh(m_JArr1[idarr], m_GArr[idarr], nH_mean, fignamestr, masscArr=masscArr, masscnameArr=masscnameArr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, propArr=propArr[idarr], takelog_prop=True, mrange_x=mrange_x, mrange_y=mrange_y, ipltptcls=ipltptcls, pr=pr)
make_plot_mJ_mh(m_JArr1, m_GArr, nH_mean, fignamestr, masscArr=masscArr, masscnameArr=masscnameArr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, propArr=propArr, takelog_prop=True, mrange_x=mrange_x, mrange_y=mrange_y, ipltptcls=0, pr=pr)


fignameh = "mJmax"; # z_mJ = zcol; "mJ%.1f" % z_mJ
fignamestr = [fignameh, collcond, subpath]
m_JArr_all, m_JArr1 = get_jeans_mass_arr(fignameh, denArr_all[ivalid, :], tempArr_all[ivalid, :], zredArr_all)
make_plot_mJ_mh(m_JArr1, m_GArr, nH_mean, fignamestr, masscArr=masscArr, masscnameArr=masscnameArr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, propArr=propArr, takelog_prop=True, mrange_x=mrange_x, mrange_y=mrange_y)
"""

figext = "eps"
fs = dict(tc = 12, tt = 10, lb = 15, lg = 12, tx = 11, cb_tc = 14, cb_lb = 17)


def color_by_prop(propArr, nbins):
	colorArr = []
	for icl in xrange(nbins):
		colorArr.append(plt.get_cmap("hsv")(float(icl)/(nbins)))

	color_id_Arr = N.zeros_like(propArr)
	#propbins = N.logspace(N.log10(colorprops*0.99), N.log10(colorprops*1.01), ncolor+1)
	p_ids, p_bins = pTdf.group_by_prop(propArr, takelog=True, fixed_interval=True, bins=[])
	return p_ids, p_bins, colorArr




#### plot

#m_GArr = grpInfoArr_all[ivalid, 2]
#propArr = denArr_all[ivalid, izcol]

def get_mass_scales(param):
	zcol = param[0]; ionizR = param[1]
	mF = cTdf.filtering_mass(zcol, 9., 1.e4) 
	#macc = cTdf.accretion_mass(simI['zcol'], N.array([simI['ionizR'][0], simI['ionizR'][0], 0.]))  ## for old code
	
	macc = cTdf.accretion_mass(zcol, ionizR)  ## to use Matt's data 
	mhc = cTdf.hoeft_characteristic_mass(zcol)
	print "mF=%.2e, macc=%.2e, mhc=%.2e" % (mF, macc, mhc)
	masscArr = [mF, macc, mhc]
	masscnameArr = [r"$M_F$", r"$M_{acc}$", r"$M_c$"]
	#masscArr = [mF, mhc]
	#masscnameArr = [r"$M_F$", r"$M_c$"]
	return masscArr, masscnameArr


def choose_ptcls_by_collapsed_redshift(denArr_org, zredArr, overdensity=200, inc_not_collapsed=True, **simI):
	if simI['convert_density']:
		denArr = denArr_org*cst.m_p/cTdf.mean_hydrogen_numden_z(zredArr)
	else:
		denArr = denArr_org.copy()

	izcoll = -1*N.ones(len(denArr), dtype="int")
	n_zred = len(zredArr)
	for icol in xrange(n_zred):
		icolp = n_zred -1*icol - 1 
		ii = N.where(denArr[:, icolp] >= overdensity)[0]
		print icolp, len(ii)
		izcoll[ii] = icolp

	if inc_not_collapsed:
		idarr = N.where(N.any((izcoll == -1, izcoll == n_zred - 1), axis=0))[0]
	else:
		idarr = N.where(izcoll == n_zred - 1)[0]

	return izcoll, idarr


def get_colorbar(propArr, prange_param):
	takelog_prop = prange_param[0]; pr = prange_param[1]; nticks = prange_param[2] 

	cbar_labels = []
	if len(pr) == 0:
		pr = N.array([N.min(propArr), N.max(propArr)]) ## if takelog_prop is true, propArr contains probably already log values
	else:
		if takelog_prop:
			pr = N.log10(pr)

	pri = N.ceil(pr)
	#for i in xrange(len(pr)):
	#	pri[i] = [pri[i], pri[i] - 1.][pri[i] < 0.]

	cbar_ticks = N.linspace(pri[0], pri[1], nticks)
	print cbar_ticks

	if takelog_prop:
		for ctick in cbar_ticks:
			cbar_labels.append(r"$10^{%d}$" % ctick)
		for ctick in cbar_ticks:
			cbar_labels.append("%.2e" % ctick)

	print "color bar is going to be: {:}".format(cbar_ticks)

	return pr, cbar_ticks, cbar_labels


def get_plot_range(pltrangevals, xarr, yarr):

	xr = pltrangevals[0]; yr = pltrangevals[1]
	fr = [1., 1.] # fr = [0.999, 1.001]

	if len(xr) == 0:
		xr = [fr[0]*N.min(xarr), fr[1]*N.max(xarr)]
	if len(yr) == 0:
		yr = [fr[0]*N.min(yarr), fr[1]*N.max(yarr)]

	if pltrangevals[2] and len(xr) == 0 and len(yr) == 0:
		xr = [fr[0]*N.min([N.min(xarr), N.min(yarr)]), fr[1]*N.max([N.max(xarr), N.max(yarr)])]
		yr = xr

	xr = N.array(xr); yr = N.array(yr)
	print "plot ranges are got: xr={:}, yr={:}".format(xr, yr)
	return xr, yr

	
def make_plot_mJ_mh_range_by_den_zred_func_paper(data, figparam):
	"""
	for paper, only two density bins are plotted - delta < 10 and delta > 200
	"""

	fig, figdim, p_bins, ptext, mrange_x, mrange_y, pr, ylabel, ws, hs, masscArr, masscnameArr, mJmax_est, zcoltext, show_legend = figparam[0], figparam[1], figparam[2], figparam[3], figparam[4], figparam[5], figparam[6], figparam[7], figparam[8], figparam[9], figparam[10], figparam[11], figparam[12], figparam[13], figparam[14]
	m_xArr, m_yArr, propArr = data[0], data[1], data[2]
	nbins_prop = len(p_bins) - 1
	for ibin in xrange(2):
		ax = fig.add_subplot(figdim[0], figdim[1], figdim[2] + ibin + 1)

		p_ids = N.where(N.all([propArr >= p_bins[2*ibin], propArr < p_bins[2*ibin+1]], axis=0))[0]
		print "In log10(delta)=[%.2e, %.2e): %d" % (p_bins[ibin], p_bins[ibin+1], len(p_ids))
		#print "test print:"
		#print propArr[p_ids]
		im = ax.scatter(m_xArr[p_ids], m_yArr[p_ids], c=propArr[p_ids], alpha=0.7, edgecolor="none", s=7, vmin=pr[0], vmax=pr[1])

		## plot the mass given in the literature
		colorL = ['green', 'orange', 'cyan']; icl=0
		if len(masscArr) > 0:
			for massc, masscname in zip(masscArr, masscnameArr):
				print "%s=%.4e" % (masscname, massc)
				ax.axvline(massc, ls='--', color=colorL[icl], lw=2.5, label=masscname)
				print "%d: %d, %d" % (figdim[2] + ibin, (figdim[2] + ibin)/figdim[1], (figdim[2] + ibin) % figdim[1])
				##if figdim[2] + ibin == 0 :
				#if (figdim[2] + ibin)/figdim[1] == 1 and (figdim[2] + ibin) % figdim[1] == 0:
					#print "PRINTING TEXT: %d" % figdim[2]
					#ax.text(massc, 3.e11, "%s" % (masscname), fontsize=fs["tc"], weight='bold')
				icl += 1
		if mJmax_est != 0.:
			ax.axhline(mJmax_est, lw = 1.5, color='m', label=r"$M_J^{max,SC}$")
			print "Plotting estimated max Jeans mass"


		## plot guide line
		acstArr = N.array([1., 4.])
		#mrange = N.array([N.min([mrange_x[0], mrange_y[0]]), N.max([mrange_x[1], mrange_y[1]])])
		for iacst, acst in enumerate(acstArr):
			if iacst == 0:
				ax.plot(mrange_x, acst*mrange_x, color='k', ls=':', lw=2.5, label=r"$M_h \propto M_J^{max}$")
			else:
				ax.plot(mrange_x, acst*mrange_x, color='k', ls=':', lw=2.5)
			"""
			if acst*mrange_x[1] > mrange_y[1]:
				ax.text(0.9*mrange_y[1]/acst, mrange_y[1], "%.1f" % acst, fontsize=fs["tc"])
			else:
				ax.text(mrange_x[1], acst*mrange_x[1], "%.1f" % acst, fontsize=fs["tc"])
			"""

		#ax.set_aspect("equal")

		ax.set_xlim([0.6*mrange_x[0], 1.1*mrange_x[1]]); 
		ax.set_ylim([0.9*mrange_y[0], 1.1*mrange_y[1]]); 

		"""
		if figdim[2] == 0: ## density range
			print ptext
			if ibin == 0:
				ax.text(0.05, 0.9, "[ , %s)" % (ptext[ibin]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			elif ibin == nbins_prop - 1:
				ax.text(0.05, 0.9, "[%s, )" % (ptext[ibin-1]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			else:
				ax.text(0.05, 0.9, "[%s, %s)" % (ptext[ibin-1], ptext[ibin]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			#ax.text(mrange_x[0]*1.01, mrange_y[1]*0.8, "%d" % figdim[2])
		"""

		ax.set_xscale("log")
		ax.set_yscale("log")

		ax.text(0.05, 0.85, zcoltext, fontsize=fs['tx'], transform=ax.transAxes, horizontalalignment='left') 

		if ws == 0:
			if ibin == 0:
				ax.set_ylabel(ylabel, fontsize=fs['lb'])
				plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])
			else:
				plt.setp(ax.get_yticklabels(), visible=False)
		else:
			if ibin == 0:
				ax.set_ylabel(ylabel, fontsize=fs['lb'])
			plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])

		if hs == 0:
			if figdim[2] == figdim[1]*(figdim[0] - 1):
				ax.set_xlabel(r"$M_h(M_\odot)$", fontsize=fs['lb'])
				plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])
			else:
				plt.setp(ax.get_xticklabels(), visible=False)
		else:
			if figdim[2] == figdim[1]*(figdim[0] - 1):
				ax.set_xlabel(r"$M_h(M_\odot)$", fontsize=fs['lb'])
			plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])

		if show_legend and ibin == 0:
			lfont = legend_font_props = FontProperties(size=fs['lg'])
			#plt.legend(plot_axes, varprops_str, bbox_to_anchor=(-0.1, 1.), loc=3, borderaxespad=0., prop=lfont, ncol=5, frameon=False)
			plt.legend(bbox_to_anchor=(0., 1.), loc=3, borderaxespad=0., prop=lfont, ncol=5, frameon=False)

	return im, ax


def make_plot_mJ_mh_range_by_den_zred_func(data, figparam):

	fig, figdim, p_bins, ptext, mrange_x, mrange_y, pr, ylabel, ws, hs, masscArr, masscnameArr, mJmax_est, zcoltext, show_legend = figparam[0], figparam[1], figparam[2], figparam[3], figparam[4], figparam[5], figparam[6], figparam[7], figparam[8], figparam[9], figparam[10], figparam[11], figparam[12], figparam[13], figparam[14]
	m_xArr, m_yArr, propArr = data[0], data[1], data[2]
	nbins_prop = len(p_bins) - 1
	for ibin in xrange(nbins_prop):
		ax = fig.add_subplot(figdim[0], figdim[1], figdim[2] + ibin + 1)

		p_ids = N.where(N.all([propArr >= p_bins[ibin], propArr < p_bins[ibin+1]], axis=0))[0]
		print "In log10(delta)=[%.2e, %.2e): %d" % (p_bins[ibin], p_bins[ibin+1], len(p_ids))
		#print "test print:"
		#print propArr[p_ids]
		im = ax.scatter(m_xArr[p_ids], m_yArr[p_ids], c=propArr[p_ids], alpha=0.7, edgecolor="none", s=7, vmin=pr[0], vmax=pr[1])

		## plot the mass given in the literature
		colorL = ['green', 'orange', 'cyan']; icl=0
		if len(masscArr) > 0:
			for massc, masscname in zip(masscArr, masscnameArr):
				print "%s=%.4e" % (masscname, massc)
				ax.axvline(massc, ls='--', color=colorL[icl], lw=2.5)
				#if figdim[2] + ibin == 0 :
					#ax.text(massc, mrange_y[0]*2., "%s" % (masscname), fontsize=fs["tx"])
				icl += 1
				#ax.text(massc, mrange_y[0]*1.1, "%s=%.2e" % (masscname, massc), fontsize=fs["tc"])

		## plot guide line
		acstArr = N.array([1., 4.])
		#mrange = N.array([N.min([mrange_x[0], mrange_y[0]]), N.max([mrange_x[1], mrange_y[1]])])
		for acst in acstArr:
			ax.plot(mrange_x, acst*mrange_x, color='k', ls=':', lw=2.5)
			"""
			if acst*mrange_x[1] > mrange_y[1]:
				ax.text(0.9*mrange_y[1]/acst, mrange_y[1], "%.1f" % acst, fontsize=fs["tx"])
			else:
				ax.text(mrange_x[1], acst*mrange_x[1], "%.1f" % acst, fontsize=fs["tx"])
			"""

		#ax.set_aspect("equal")

		ax.set_xlim([0.6*mrange_x[0], 1.1*mrange_x[1]]); 
		ax.set_ylim([0.9*mrange_y[0], 1.1*mrange_y[1]]); 
		"""
		if figdim[2] == 0:
			print ptext
			if ibin == 0:
				ax.text(0.05, 0.9, "[ , %s)" % (ptext[ibin]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			elif ibin == nbins_prop - 1:
				ax.text(0.05, 0.9, "[%s, )" % (ptext[ibin-1]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			else:
				ax.text(0.05, 0.9, "[%s, %s)" % (ptext[ibin-1], ptext[ibin]), fontsize=fs["tx"], transform=ax.transAxes, horizontalalignment='left')
			#ax.text(mrange_x[0]*1.01, mrange_y[1]*0.8, "%d" % figdim[2])
		"""

		ax.set_xscale("log")
		ax.set_yscale("log")

		ax.text(0.05, 0.85, zcoltext, fontsize=fs['tx'], transform=ax.transAxes, horizontalalignment='left') 

		if ws == 0:
			if ibin == 0:
				ax.set_ylabel(ylabel, fontsize=fs['lb'])
				plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])
			else:
				plt.setp(ax.get_yticklabels(), visible=False)
		else:
			if ibin == 0:
				ax.set_ylabel(ylabel, fontsize=fs['lb'])
			plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])

		if hs == 0:
			if figdim[2] == figdim[1]*(figdim[0] - 1):
				ax.set_xlabel(r"$M_h$", fontsize=fs['lb'])
				plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])
			else:
				plt.setp(ax.get_xticklabels(), visible=False)
		else:
			if figdim[2] == figdim[1]*(figdim[0] - 1):
				ax.set_xlabel(r"$M_h$", fontsize=fs['lb'])
			plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])

	return im, ax




def make_plot_mJ_mh_range_by_den_zred(data, param, paper=True):
#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False
## subplt_bins: for the separate plot
	"""
	make separate plot for a few range of overdensity with different redshifts
	"""
	print param
	zcolarr, zcolarrtext, ionizR, figname, bin_by_grp, pltrangevals, takelog_prop, prange_param, subplt_bins, equal_n_subplt, saperate_star_ptcls, plot_mJEST = param[0], param[1], param[2], param[3], param[4], param[5], param[6][0], param[6], param[7][0], param[7][1], param[7][2], param[8]
	nzcol = len(zcolarr); narr = len(data)/nzcol

	if figname.find("max") != -1:
		ylabel = r"$M_J^{max}(M_\odot)$"
	elif figname.find("peak") != -1: 
		ylabel = r"$M_J^{peak}(M_\odot)$" 
	else:
		ylabel = r"$M_J(%s)(M_\odot)$" % fignameh[2:]

	## to make plot
	fig = plt.figure(1)
	plt.clf()
	#plt.suptitle(r"$\Gamma=%d$" % (ionizR))
	hs = [0., 0.2][len(pltrangevals[0]) == 0]
	ws = [0., 0.2][len(pltrangevals[1]) == 0]
	plt.subplots_adjust(hspace=hs, wspace=ws, top=0.9, bottom=0.15, left=0.15, right=0.9)
	if paper:
		nbins_prop = len(subplt_bins) - 2; figdim = [nzcol, nbins_prop, 0]
		ptext = ["%.0f" % subplt_bins[i] for i in xrange(1, nbins_prop)]
	else:
		nbins_prop = len(subplt_bins) - 1; figdim = [nzcol, nbins_prop, 0]
		ptext = ["%.0f" % subplt_bins[i] for i in xrange(1, nbins_prop)]

	if takelog_prop:
		subplt_bins = N.log10(subplt_bins) ## for the separate plot

	for izcol in xrange(nzcol):
		show_legend = 1 if izcol == 0 else 0
		zcoltext = "z=%.2f" % zcolarr[izcol] if len(zcolarrtext) == 0 else zcolarrtext[izcol]	
		xarr = data[izcol*narr]; yarr = data[izcol*narr+1]; parr = data[izcol*narr+2]
	
		if takelog_prop:
			parr = N.log10(parr + 1.e-10)

		## get the data array
		masscArr, masscnameArr = get_mass_scales((zcolarr[izcol], ionizR))

		## estimate max Jeans mass
		if plot_mJEST:
			mJEst, nHJest, tJest = cTdf.estimate_max_jeans_mass(zcolarr[izcol])
		else:
			mJEst = 0.

		## get color bar related properties
		clrrange, cbar_ticks, cbar_labels = get_colorbar(parr, prange_param)
		
		## plot range
		xr, yr = get_plot_range(pltrangevals, xarr, yarr)

		## plot particles by having equal number of particles in each of group mass bin
		## this only probably necessary when I take all the particles not using density or virial radius criterion
		if len(bin_by_grp) > 0: 
			## bin_by_grp[0] = nbins, bin_by_grp[1] = nptcls in each bin
			idArr_plt = pmf.bin_by_grpmass(xarr, (xr, bin_by_grp[0]), bin_by_grp[1])
			print "particles to plot: %d" % len(idArr_plt)
			xarr = xarr[idArr_plt]; yarr = yarr[idArr_plt]; parr = parr[idArr_plt]

		if equal_n_subplt:
			if not saperate_star_ptcls:
				idArr_plt2 = pmf.bin_by_grpmass(parr, (subplt_bins, 0), 0)
			else:
				idArr_plt2 = pmf.bin_by_grpmass(parr, (subplt_bins[:-1], 0), 0)
				istars = N.where(parr >= 4.)[0] ## 1.e4: random density for star particles -- make sure i'm using a correct density value
				idArr_plt2 = N.append(idArr_plt2, istars)
			xarr = xarr[idArr_plt2]; yarr = yarr[idArr_plt2]; parr = parr[idArr_plt2]

		## make plot
		#fig, figdim, p_bins, ptext, mrange_x, mrange_y, pr, ylabel, ws, hs, masscArr, masscnameArr, zcoltext
		if paper:
			#im, ax = make_plot_mJ_mh_range_by_den_zred_func_paper([xarr, yarr, parr], [fig, figdim, subplt_bins, ptext, xr, yr, clrrange, ylabel, ws, hs, masscArr, masscnameArr, mJEst, zcoltext, show_legend])
			## updated for the revision 02/2014 to remove reference mass scales
			im, ax = make_plot_mJ_mh_range_by_den_zred_func_paper([xarr, yarr, parr], [fig, figdim, subplt_bins, ptext, xr, yr, clrrange, ylabel, ws, hs, [], [], mJEst, zcoltext, show_legend])
		else:
			im, ax = make_plot_mJ_mh_range_by_den_zred_func([xarr, yarr, parr], [fig, figdim, subplt_bins, ptext, xr, yr, clrrange, ylabel, ws, hs, masscArr, masscnameArr, mJEst], zcoltext, show_legend)

		"""
		if len(zcolarrtext) == 0:
			ax.text(0.05, 0.85, "z=%.2f" % zcolarr[izcol], fontsize=fs['tx'], transform=ax.transAxes, horizontalalignment='left') 
		else:
			ax.text(0.05, 0.85, zcolarrtext[izcol], fontsize=fs['tx'], transform=ax.transAxes, horizontalalignment='left') 
		"""

		figdim[2] = (izcol + 1)*nbins_prop

	## color  bar
	cax = fig.add_axes([0.91, 0.15, 0.01, 0.75])
	cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks)
	cbar.ax.set_yticklabels(cbar_labels, fontsize=fs['cb_tc'])
	cbar.ax.set_ylabel(r"$\delta$", rotation=270, fontsize=fs['cb_lb'])
	cbar.ax.yaxis.set_label_coords(4., 0.45)
	print cbar_ticks, cbar_labels

	print figname
	plt.savefig(figname)


###################

def make_plot_mJ_mh_range_by_den_zred_old(zcolArr, pparam, rparam, simI):
#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False
	"""
	make separate plot for a few range of overdensity with different redshifts
	"""
	fignameh = pparam[0]; bin_by_grp = pparam[1]; pltrangevals = pparam[2]; pparam[3] = takelog_prop, color_pr = pparam[4]; p_bins = pparam[5] #=[1., 400.]=N.array([1.e-5, 10., 100., 100000.])

	fs = dict(tc = 9, tt = 10, lb = 9, lg = 6, tx = 6)
	if fignameh.find("max") != -1:
		ylabel = r"$M_J^{max}$"
	elif fignameh.find("peak") != -1: 
		ylabel = r"$M_J^{peak}$" 
	else:
		ylabel = r"$M_J(%s)$" % fignameh[2:]

	## to make plot
	fig = plt.figure(1)
	plt.clf()
	plt.suptitle(r"$\Gamma=%d$" % (simI["ionizR"][0]))
	plt.subplots_adjust(hspace=0.2, wspace=0., top=0.9, bottom=0.1)
	nbins_prop = len(p_bins) - 1; figdim = [len(zcolArr), nbins_prop, 0]

	if takelog_prop:
		p_bins = N.log10(p_bins) ## for the separate plot
		if len(color_pr) > 0:
			color_pr = N.log10(N.array(color_pr))

	for irow, zcol in enumerate(zcolArr):
		simI['zcol'] = zcol
		## get the data array
		DC = RO.DataConstruc(**simI)
		izcol, zarr, darr, tarr, garr = DC.get_entire_array2(rparam)
		#nH_mean, collcond, data_abspath, subpath, masscArr, masscnameArr = get_some_values(**simI)
		m_GArr = grpInfoArr_all[ivalid, 2]
		m_JArr_all, m_JArr1 = get_jeans_mass_arr(fignameh, denArr_all[ivalid, :], tempArr_all[ivalid, :], zredArr_all)

		## for color
		propArr = denArr_all[ivalid, izcol]/nH_mean
		if takelog_prop:
			propArr = N.log10(propArr) ## for the color
		if len(color_pr) == 0:
			color_pr = [N.min(propArr), N.max(propArr)] 

		## plot range
		if len(pltrangevals[0]) == 0 and len(pltrangevals[1]) == 0:
			if pltrangevals[2]:
				mrange_x = [N.min([N.min(m_GArr), N.min(m_JArr1)]), N.max([N.max(m_GArr), N.max(m_JArr1)])]
				mrange_x[0] = mrange_x[0]*0.999; mrange_x[1] = mrange_x[1]*1.001
				mrange_y = mrange_x
			else:
				mrange_x = [0.999*N.min(m_GArr), 1.001*N.max(m_GArr)]
				mrange_y = [0.999*N.min(m_JArr1), 1.001*N.max(m_JArr1)]
		else:
			mrange_x = pltrangevals[0].copy(); mrange_y = pltrangevals[1].copy()
		mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

		## plot particles by having equal number of particles in each of group mass bin
		## this only probably necessary when I take all the particles not using density or virial radius criterion
		if len(bin_by_grp) > 0: ## bin_by_grp[0] = nbins, bin_by_grp[1] = nptcls in each bin
			idArr_plt = bin_by_grpmass(m_GArr, mrange_x, bin_by_grp[0], bin_by_grp[1])
			print "particles to plot: %d" % len(idArr_plt)
			m_xArr = m_GArr[idArr_plt]; m_yArr = m_JArr1[idArr_plt]
			propArr = propArr[idArr_plt]
		else:
			m_xArr = m_GArr; m_yArr = m_JArr1

		## make plot
		im = make_plot_mJ_mh_range_by_den_zred_func(fig, figdim, m_xArr, m_yArr, propArr, p_bins, mrange_x, mrange_y, color_pr, ylabel, masscArr=masscArr, masscnameArr=masscnameArr)
		figdim[2] = (irow + 1)*nbins_prop

	## color  bar
	cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
	cbar_ticks = N.linspace(color_pr[0], color_pr[1], 5)
	cbar_labels = []
	for ctick in cbar_ticks:
		cbar_labels.append("%.1f" % 10.**ctick)
	cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks)
	#cbar = fig.colorbar(cax, ticks=[N.min(propArr), nH_mean, N.max(propArr)]) #cbar_ticks)
	cbar.ax.set_yticklabels(cbar_labels, fontsize=fs['lb'])
	#cbar.ax.set_yticklabels(["%.2e" % (10.**N.min(propArr)), "%.2e" % (10.**nH_mean), "%.2e" % 10.**N.max(propArr)], fontsize=fs['lb'])
	print cbar_ticks, cbar_labels

	istr = collcond.find("grpm")
	figname = "mG-%s_%s_ov%d_vr%d_eqN%d_clsubs.%s" % (fignameh, collcond[istr:], int(simI["over_density"]), simI["virial_radius"], int(bin_by_mgrp), figext)
	print figname
	plt.savefig(os.path.join(subpath, figname))


def get_mJ_mh_range_by_den_ids(data, param):
#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False
## subplt_bins: for the separate plot
	"""
	make separate plot for a few range of overdensity with different redshifts
	"""
	print param
	bin_by_grp = param[0]; pltrangevals = param[1]; takelog_prop = param[2][0]; prange_param = param[2]; subplt_bins = param[3][0]; equal_n_subplt = param[3][1]; saperate_star_ptcls = param[3][2] 

	nbins_prop = len(subplt_bins) - 1

	if takelog_prop:
		subplt_bins = N.log10(subplt_bins) ## for the separate plot

	xarr = data[0]; yarr = data[1]; parr = data[2]
	
	if takelog_prop:
		parr = N.log10(parr + 1.e-10)

	## get color bar related properties
	clrrange, cbar_ticks, cbar_labels = get_colorbar(parr, prange_param)
		
	## plot range
	xr, yr = get_plot_range(pltrangevals, xarr, yarr)

	## plot particles by having equal number of particles in each of group mass bin
	## this only probably necessary when I take all the particles not using density or virial radius criterion
	if len(bin_by_grp) > 0: 
		## bin_by_grp[0] = nbins, bin_by_grp[1] = nptcls in each bin
		idArr_plt = pmf.bin_by_grpmass(xarr, (xr, bin_by_grp[0]), bin_by_grp[1])
		print "particles to plot: %d" % len(idArr_plt)
		xarr = xarr[idArr_plt]; yarr = yarr[idArr_plt]; parr = parr[idArr_plt]

	if equal_n_subplt:
		if not saperate_star_ptcls:
			idArr_plt2 = pmf.bin_by_grpmass(parr, (subplt_bins, 0), 0)
		else:
			idArr_plt2 = pmf.bin_by_grpmass(parr, (subplt_bins[:-1], 0), 0)
			istars = N.where(parr >= 4.)[0] ## 1.e4: random density for star particles -- make sure i'm using a correct density value
			idArr_plt2 = N.append(idArr_plt2, istars)
	else:
		idArr_plt2 = idArr_plt

	return idArr_plt2


def make_plot_mJ_mh_range_by_den(m_JArr, m_GArr, propArr, nH_mean, fignamestr, masscArr=[], masscnameArr=[], bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, takelog_prop=True, mrange_x=[], mrange_y=[], ipltptcls=2, pr=[]):
#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False
	"""
	make separate plot for a few range of overdensity
	"""

	fs = dict(tc = 9, tt = 10, lb = 9, lg = 6, tx = 6)
	if fignamestr[0].find("max") != -1:
		ylabel = r"$M_J^{max}$"
	elif fignamestr[0].find("peak") != -1: 
		ylabel = r"$M_J^{peak}$" 
	else:
		ylabel = r"$M_J(%s)$" % fignamestr[0][2:]


	if reset_mrange:
		if equal_range:
			mrange_x = [N.min([N.min(m_GArr), N.min(m_JArr)]), N.max([N.max(m_GArr), N.max(m_JArr)])]
			mrange_x[0] = mrange_x[0]*0.999; mrange_x[1] = mrange_x[1]*1.001
			mrange_y = mrange_x
		else:
			mrange_x = [0.999*N.min(m_GArr), 1.001*N.max(m_GArr)]
			mrange_y = [0.999*N.min(m_JArr), 1.001*N.max(m_JArr)]

	mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

	## plot particles by having equal number of particles in each of group mass bin
	## this only probably necessary when I take all the particles not using density or virial radius criterion
	if bin_by_mgrp:
		idArr_plt = bin_by_grpmass(m_GArr, mrange_x, nbins, nptcl_each)
		print "particles to plot: %d" % len(idArr_plt)
		m_xArr = m_GArr[idArr_plt]; m_yArr = m_JArr[idArr_plt]
	else:
		m_xArr = m_GArr; m_yArr = m_JArr

	fig = plt.figure(1)
	plt.clf()
	mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

	""" simple histogram like color assignment
	id_clprop, clprop_bins, colorArr = color_by_prop(propArr, ncolor)
	for icl in xrange(ncolor):
		iptcl_cl = N.where(id_clprop == icl)[0]
		ax.plot(m_xArr[iptcl_cl], m_yArr[iptcl_cl], marker="o", color=colorArr[icl], linestyle='None', ms=2.5, mew=0., mec="b", alpha=0.7)
	"""
	cbar_labels = []
	if takelog_prop:
		propArr = N.log10(propArr)
		if len(pr) == 0:
			cbar_ticks = N.linspace(N.min(propArr), N.max(propArr), 5)
		else:
			cbar_ticks = N.linspace(pr[0], pr[1], 5)
		#cbar_ticks = [N.min(propArr), nH_mean, N.max(propArr)]
		for ctick in cbar_ticks:
			cbar_labels.append("%.1f" % 10.**ctick)
	else:
		if len(pr) == 0:
			cbar_ticks = [N.min(propArr), nH_mean, N.max(propArr)]
		else:
			cbar_ticks = N.linspace(pr[0], pr[1], 5)
		for ctick in cbar_ticks:
			cbar_labels.append("%.2e" % ctick)
	print N.min(propArr), N.max(propArr)

	##make plot
	p_bins = N.array([1.e-5, 10., 100., 100000.])
	if takelog_prop:
		p_bins = N.log10(p_bins)
	nbin = len(p_bins) - 1

	plt.suptitle(r"$z_{col}=%.1f,\ \Gamma=%d, n_H > %d\ \bar{n_H}(=%.2e)$" % (simI["zcol"], simI["ionizR"][0], int(simI["over_density"]), nH_mean))
	plt.subplots_adjust(hspace=0.3, wspace=0., top=0.7, bottom=0.3)

	for ibin in xrange(nbin):
		ax = fig.add_subplot(1, nbin, ibin+1)

		p_ids = N.where(N.all([propArr >= p_bins[ibin], propArr < p_bins[ibin+1]], axis=0))[0]
		print "In log10(delta)=[%.2e, %.2e): %d" % (p_bins[ibin], p_bins[ibin+1], len(p_ids))
		#print "test print:"
		#print propArr[p_ids]
		im = ax.scatter(m_xArr[p_ids], m_yArr[p_ids], c=propArr[p_ids], alpha=0.7, edgecolor="none", s=10, vmin=pr[0], vmax=pr[1])

		## plot the mass given in the literature
		if len(masscArr) > 0:
			for massc, masscname in zip(masscArr, masscnameArr):
				print "%s=%.4e" % (masscname, massc)
				ax.axvline(massc, ls = '--', color='green')
				ax.text(massc, mrange_y[0]*1.1, "%s" % (masscname), fontsize=fs["tc"])
				#ax.text(massc, mrange_y[0]*1.1, "%s=%.2e" % (masscname, massc), fontsize=fs["tc"])

		## plot guide line
		acstArr = N.logspace(-2., 2., 5, base=2)
		#mrange = N.array([N.min([mrange_x[0], mrange_y[0]]), N.max([mrange_x[1], mrange_y[1]])])
		for acst in acstArr:
			ax.plot(mrange_x, acst*mrange_x, "m:")
			if acst*mrange_x[1] > mrange_y[1]:
				ax.text(0.9*mrange_y[1]/acst, mrange_y[1], "%.1f" % acst, fontsize=fs["tc"])
			else:
				ax.text(mrange_x[1], acst*mrange_x[1], "%.1f" % acst, fontsize=fs["tc"])

		#ax.set_aspect("equal")

		ax.set_xlim(mrange_x); ax.set_ylim(mrange_y)
		ax.text(mrange_x[0], mrange_y[1]*1., "[%.0f-%.0f)" % (p_bins[ibin], p_bins[ibin+1]), fontsize=fs["tc"])

		ax.set_xscale("log")
		ax.set_yscale("log")

		if ibin == 0:
			ax.set_ylabel(ylabel)
			plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])
		else:
			plt.setp(ax.get_yticklabels(), visible=False)

		ax.set_xlabel(r"$M_h (M_\odot)$", fontsize=fs['lb'])
		plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])

	## color  bar

	cax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
	cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks)
	#cbar = fig.colorbar(cax, ticks=[N.min(propArr), nH_mean, N.max(propArr)]) #cbar_ticks)
	cbar.ax.set_yticklabels(cbar_labels, fontsize=fs['lb'])
	#cbar.ax.set_yticklabels(["%.2e" % (10.**N.min(propArr)), "%.2e" % (10.**nH_mean), "%.2e" % 10.**N.max(propArr)], fontsize=fs['lb'])
	print cbar_ticks, cbar_labels

	figname = "mG-%s_%s_ov%d_vr%d_eqN%d_%d_clsubs.%s" % (fignamestr[0], fignamestr[1], int(simI["over_density"]), simI["virial_radius"], int(bin_by_mgrp), ipltptcls, figext)
	plt.savefig(os.path.join(fignamestr[2], figname))



def make_plot_mJ_mh(m_JArr, m_GArr, nH_mean, fignamestr, masscArr=[], masscnameArr=[], bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False, propArr=[], takelog_prop=True, mrange_x=[], mrange_y=[], ipltptcls=2, pr=[]):
#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False

	fs = dict(tc = 9, tt = 10, lb = 9, lg = 6, tx = 6)
	if fignamestr[0].find("max") != -1:
		ylabel = r"$M_J^{max}$"
	elif fignamestr[0].find("peak") != -1: 
		ylabel = r"$M_J^{peak}$" 
	else:
		ylabel = r"$M_J(%s)$" % fignamestr[0][2:]


	if reset_mrange:
		if equal_range:
			mrange_x = [N.min([N.min(m_GArr), N.min(m_JArr)]), N.max([N.max(m_GArr), N.max(m_JArr)])]
			mrange_x[0] = mrange_x[0]*0.999; mrange_x[1] = mrange_x[1]*1.001
			mrange_y = mrange_x
		else:
			mrange_x = [0.999*N.min(m_GArr), 1.001*N.max(m_GArr)]
			mrange_y = [0.999*N.min(m_JArr), 1.001*N.max(m_JArr)]

	mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

	## plot particles by having equal number of particles in each of group mass bin
	## this only probably necessary when I take all the particles not using density or virial radius criterion
	if bin_by_mgrp:
		idArr_plt = bin_by_grpmass(m_GArr, mrange_x, nbins, nptcl_each)
		print "particles to plot: %d" % len(idArr_plt)
		m_xArr = m_GArr[idArr_plt]; m_yArr = m_JArr[idArr_plt]
	else:
		m_xArr = m_GArr; m_yArr = m_JArr

	fig = plt.figure(1)
	plt.clf()
	ax = fig.add_subplot(111)
	mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

	if len(propArr) > 0:
		propArr = propArr[idArr_plt]
		""" simple histogram like color assignment
		id_clprop, clprop_bins, colorArr = color_by_prop(propArr, ncolor)
		for icl in xrange(ncolor):
			iptcl_cl = N.where(id_clprop == icl)[0]
			ax.plot(m_xArr[iptcl_cl], m_yArr[iptcl_cl], marker="o", color=colorArr[icl], linestyle='None', ms=2.5, mew=0., mec="b", alpha=0.7)
		"""
		cbar_labels = []
		if takelog_prop:
			propArr = N.log10(propArr)
			if len(pr) == 0:
				cbar_ticks = N.linspace(N.min(propArr), N.max(propArr), 5)
			else:
				cbar_ticks = N.linspace(pr[0], pr[1], 5)
			#cbar_ticks = [N.min(propArr), nH_mean, N.max(propArr)]
			for ctick in cbar_ticks:
				cbar_labels.append("%.2e" % 10.**ctick)
		else:
			if len(pr) == 0:
				cbar_ticks = [N.min(propArr), nH_mean, N.max(propArr)]
			else:
				cbar_ticks = N.linspace(pr[0], pr[1], 5)
			for ctick in cbar_ticks:
				cbar_labels.append("%.2e" % ctick)
		cax = ax.scatter(m_xArr, m_yArr, c=propArr, alpha=0.7, edgecolor="none", s=10, vmin=pr[0], vmax=pr[1])
		#cbar = fig.colorbar(cax)
		print N.max(propArr)
		cbar = fig.colorbar(cax, ticks=cbar_ticks)
		#cbar = fig.colorbar(cax, ticks=[N.min(propArr), nH_mean, N.max(propArr)]) #cbar_ticks)
		cbar.ax.set_yticklabels(cbar_labels, fontsize=fs['lb'])
		#cbar.ax.set_yticklabels(["%.2e" % (10.**N.min(propArr)), "%.2e" % (10.**nH_mean), "%.2e" % 10.**N.max(propArr)], fontsize=fs['lb'])
		print cbar_ticks, cbar_labels
	else:
		ax.plot(m_xArr, m_yArr, marker="o", color="b", linestyle='None', ms=2.5, mew=0., mec="none", alpha=0.7)

	## plot the mass given in the literature
	if len(masscArr) > 0:
		for massc, masscname in zip(masscArr, masscnameArr):
			print "%s=%.4e" % (masscname, massc)
			ax.axvline(massc, ls = ':', color='magenta')
			ax.text(massc, mrange_y[0]*1.1, "%s" % (masscname), fontsize=fs["tc"])
			#ax.text(massc, mrange_y[0]*1.1, "%s=%.2e" % (masscname, massc), fontsize=fs["tc"])

	## plot guide line
	acstArr = N.logspace(-2., 2., 5, base=2)
	#mrange = N.array([N.min([mrange_x[0], mrange_y[0]]), N.max([mrange_x[1], mrange_y[1]])])
	for acst in acstArr:
		ax.plot(mrange_x, acst*mrange_x, "r-")
		if acst*mrange_x[1] > mrange_y[1]:
			ax.text(0.9*mrange_y[1]/acst, mrange_y[1], "%.1f" % acst, fontsize=fs["tc"])
		else:
			ax.text(mrange_x[1], acst*mrange_x[1], "%.1f" % acst, fontsize=fs["tc"])

	#ax.set_aspect("equal")

	ax.set_xlim(mrange_x); ax.set_ylim(mrange_y)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel(r"$M_h$", fontsize=fs['lb'])
	ax.set_ylabel(ylabel)

	ax.set_title(r"$z_{col}=%.1f,\ \Gamma=%d, n_H > %d\ \bar{n_H}(=%.2e)$" % (simI["zcol"], simI["ionizR"][0], int(simI["over_density"]), nH_mean))

	figname = "mG-%s_%s_ov%d_vr%d_eqN%d_%d.%s" % (fignamestr[0], fignamestr[1], int(simI["over_density"]), simI["virial_radius"], int(bin_by_mgrp), ipltptcls, figext)
	plt.savefig(os.path.join(fignamestr[2], figname))


##########

#denArr_all[ivalid, :], tempArr_all[ivalid, :]; 
## plot n_H - T for the particles below 1*M_h
def make_plot_conditional_nH_T(denArr, tempArr, m_GArr, m_JArr, fignamestr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, propArr=[], ncolor=5):

	if fignamestr[0].find("max") != -1:
		ylabel = r"$M_J^{max}$"
	elif fignamestr[0].find("peak") != -1: 
		ylabel = r"$M_J^{peak}$" 
	else:
		ylabel = r"$M_J(%s)$" % fignamestr[0][2:]

		mrange_x = [0.9*N.min(m_GArr), 1.1*N.max(m_GArr)]
		mrange_y = [0.9*N.min(m_JArr), 1.1*N.max(m_JArr)]

	ii = N.where(m_yArr < m_xArr)[0]
	print ii

	## make seperate plot

	if len(ii) > 0 and len(propArr) > 0:
		## choose the color based on the density
		id_clprop, clprop_bins, colorArr = color_by_prop(propArr[ii], ncolor)

		fig = plt.figure(2)
		figname = "n_H-T_%s_ov%d_vr%d_%sLTmG_eqN%s_%sclr.pdf" % (fignamestr[1], int(simI["over_density"]), simI["virial_radius"], fignamestr[0], int(bin_by_mass))
		figpp = PdfPages(figname)

		nrow = 2; ncol = 2; nplt = nrow*ncol
		pltkwarg = dict(linestyle="--", marker="o", ms=5., mew=0., alpha=0.7, mec="none")

		for icl, idx in enumerate(ii):
			iplt = icl % nplt

			if iplt == 0:
				plt.clf()

			ax = fig.add_subplot(2, 2, iplt + 1)
			pltkwarg["color"] = colorArr[id_clprop[icl]]
			ax.plot(xArr[idx, :], yArr[idx, :], **pltkwarg)
			ax.plot(xArr[idx, -1], yArr[idx, -1], "k+", ms=6.)
			ax.set_xscale("log"); ax.set_yscale("log")
			#ax.set_xlim([1.e-5, 1.e-1]); ax.set_ylim([8.e3, 1.2e5])
			ax.set_xlim(N.min(xArr[ii, :]), N.max(xArr[ii, :]))
			ax.set_ylim(N.min(yArr[ii, :]), N.max(yArr[ii, :]))
			ax.set_xlabel(r"$n_H$"); ax.set_ylabel(r"$T$")

			if iplt == 3:
				figpp.savefig()

		figpp.clsoe()

	## plot in one panel
	elif len(ii) > 0 and len(ii) <= 30:
		fig = plt.figure(2)
		plt.clf()
		ax = fig.add_subplot(111)
		colorArr = []
		for icl in xrange(len(ii)):
			colorArr.append(plt.get_cmap("hsv")(float(icl)/(len(ii))))

		pltkwarg = dict(linestyle="--", marker="o", ms=5., mew=0., alpha=0.7, mec="none")

		for icl, idx in enumerate(ii):
			pltkwarg["color"] = colorArr[icl]
			ax.plot(xArr[idx, :], yArr[idx, :], **pltkwarg)

		for icl, idx in enumerate(ii):
			ax.plot(xArr[idx, -1], yArr[idx, -1], "k+", ms=6.)

		ax.set_xscale("log")
		ax.set_yscale("log")
		#ax.set_xlim([1.e-5, 1.e-1]); ax.set_ylim([8.e3, 1.2e5])
		ax.set_xlim(N.min(denArr[ii, :]), N.max(denArr[ii, :]))
		ax.set_ylim(N.min(tempArr[ii, :]), N.max(tempArr[ii, :]))
		ax.set_xlabel(r"$n_H$"); ax.set_ylabel(r"$T$")
		figname = "n_H-T_%s_ov%d_vr%d_%sLTmG_eqN%s.%s" % (fignamestr[1], int(simI["over_density"]), simI["virial_radius"], fignamestr[0], int(bin_by_mass), figext)
		plt.savefig(figname)

		


#mgrpArr_all = grpInfoArr_all[:, 2]
### make plot of fraction of number of overdensity gas particles as a function of halo mass
def make_plot_mG_nH(mgrpArr_all, denArr_all, over_density2, fignamestr, ngrp=15):
	fig = plt.figure(3)
	plt.clf()
	ax = fig.add_subplot(111)

	idArr, bins = pTdf.group_by_prop(mgrpArr_all, ngrp=ngrp, takelog=True, fixed_interval=True)
	ngrp = len(bins) - 1
	bins = 10.**bins
		
	f_nHArr = N.zeros(len(bins)-1)

	rho_mean = cTdf.mean_hydrogen_numden_z(simI["zcol"])
	for igrp in xrange(ngrp):
		isame = N.where(idArr == igrp)[0]
		print "ptcls in [%.2e]: %d" % (bins[igrp], len(isame))
		if len(isame) > 0:
			ihigh = N.where(denArr_all[isame, izcol] >= over_density2*rho_mean)[0]
			f_nHArr[igrp] = 1.*len(ihigh)/(1.*len(isame))

	print "x({:d}): {:}".format(len(bins[:-1]+bins[1:]), (bins[:-1]+bins[1:])/2.)
	print "y({:d}): {:}".format(len(f_nHArr), f_nHArr)

	ax.plot((bins[:-1]+bins[1:])/2., f_nHArr)
	ax.set_xlabel(r"$M_h$") 
	ax.set_ylabel(r"$f_{n_H > %d \bar{n_H}}(%.1f)$" % (over_density2, simI["zcol"]))
	ax.set_xscale("log")
	ax.set_xlim(bins[0], bins[-1]); ax.set_ylim(0., 1.1)
	#plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])
	#plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])

	figname = "mgrp-fnH_%s_ov%d_vr%d.%s" % (fignamestr[1], int(over_density2), simI["virial_radius"], figext)
	plt.savefig(os.path.join(fignamestr[2], figname))



#bin_by_mgrp = True; nbins = 10; nptcl_each = 0; reset_mrange = True; equal_range = False
#### plot weirdos: both scatter plot and trajectory
def make_plot_weirdos(m_GArr, m_JArr, denArr_all, tempArr_all, grpInfoArr_all, fignamestr, bin_by_mgrp = True, nbins = 10, nptcl_each = 0, reset_mrange = True, equal_range = False):
	fs = dict(tc = 8, tt = 10, lb = 8, lg = 6, tx = 6)
	if fignamestr[0].find("max") != -1:
		ylabel = r"$M_J^{max}$"
	elif fignamestr[0].find("peak") != -1: 
		ylabel = r"$M_J^{peak}$" 
	else:
		ylabel = r"$M_J(%s)$" % fignamestr[0][2:]

	if reset_mrange:
		if equal_range:
			mrange_x = [N.min([N.min(m_GArr), N.min(m_JArr)]), N.max([N.max(m_GArr), N.max(m_JArr)])]
			mrange_x[0] = mrange_x[0]*0.9; mrange_x[1] = mrange_x[1]*1.1
			mrange_y = mrange_x
		else:
			mrange_x = [0.9*N.min(m_GArr), 1.1*N.max(m_GArr)]
			mrange_y = [0.9*N.min(m_JArr), 1.1*N.max(m_JArr)]

	mrange_x = N.array(mrange_x); mrange_y = N.array(mrange_y)

	if bin_by_mgrp:
		idArr_plt = bin_by_grpmass(grpInfoArr_all[ivalid, 2], mrange_x, nbins, nptcl_each)
		m_xArr = grpInfoArr_all[ivalid[idArr_plt], 2]
		m_yArr = m_yArr[idArr_plt]
	else:
		m_xArr = grpInfoArr_all[ivalid, 2]


	## find weirdos	
	jj = N.where(m_yArr > xx*m_xArr)[0]
	#jj = N.where(N.any([grpInfoArr[:, 2] < mgrp_lim, m_yArr > xx*grpInfoArr[:, 2]], axis=0))[0]
	print "%d of %d seem to be weirdos" % (len(jj), len(ivalid))

	if len(jj) > 0 and len(jj) <= 50:
		fs = dict(tc = 8, tt = 10, lb = 8, lg = 6, tx = 6)
		pltkwarg = dict(linestyle="none", marker="o", ms=3., alpha=0.6, mec="none")
		pltkwarg_a = dict(linestyle=":", color="y", marker="None")

		fig = plt.figure(5)
		plt.clf()
		plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.08, right=0.9)

		#### plot mJ vs m_G
		ax = fig.add_subplot(121)

		pltkwarg["color"] = 'b'; ax.plot(m_xArr, m_yArr, **pltkwarg)

		## write their ids
		if len(jj) <= 30:
			for j in jj:
				pltkwarg["color"] = 'r'; ax.plot(m_xArr[j], m_yArr[j], **pltkwarg)
				ax.text(m_xArr[j], m_yArr[j], "%d" % j, color="k", fontsize=fs["tx"])
		else:
			pltkwarg["color"] = 'r'; ax.plot(m_xArr[jj], m_yArr[jj], **pltkwarg)

		## plot guide line
		acstArr = N.logspace(-2., 5., 8, base=2)
		#mrange = N.array([N.min([mrange_x[0], mrange_y[0]]), N.max([mrange_x[1], mrange_y[1]])])
		for acst in acstArr:
			ax.plot(mrange_x, acst*mrange_x, "g:")
			if acst*mrange_x[1] > mrange_y[1]:
				ax.text(0.9*mrange_y[1]/acst, mrange_y[1], "%.1f" % acst, fontsize=fs["tc"])
			else:
				ax.text(mrange_x[1], acst*mrange_x[1], "%.1f" % acst, fontsize=fs["tc"])

		#ax.set_aspect("equal")

		ax.set_xlim(mrange_x); ax.set_ylim(mrange_y)
		ax.set_xscale("log"); ax.set_yscale("log")
		ax.set_xlabel(r"$M_h$"); ax.set_ylabel(ylabel)
		plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])
		plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])

		ax.set_title(r"$z_{col}=%.1f,\ \Gamma=%d,\ n_H > %d \bar{n_H}$" % (simI["zcol"], simI["ionizR"][0], int(simI["over_density"])))

		#### plot n_H vs T of ptcls plotted in red

		ax = fig.add_subplot(122)

		colorArr = []
		for icl in xrange(len(jj)):
			colorArr.append(plt.get_cmap("hsv")(float(icl)/(len(jj))))

		pltkwarg["linestyle"] = "--" 
		for icl, j in enumerate(ivalid[jj]):
			pltkwarg["color"] = colorArr[icl]; 
			idense = N.where(denArr_all[j, :] >= denlim)[0] 
			if len(idense) > 0:
				i_to = N.min(idense)
				ax.plot(denArr_all[j, :i_to], tempArr_all[j, :i_to], label="%d" % j, **pltkwarg)
				ax.text(denArr_all[j, i_to-1], tempArr_all[j, i_to-1], "%.1f" % zredArr[i_to-1], fontsize=fs["tx"])
			else:
				ax.plot(denArr_all[j, :], tempArr_all[j, :], label="%d" % j, **pltkwarg)
			if len(jj) <= 30:
				lfont = FontProperties(size=fs["lg"])
				plt.legend(bbox_to_anchor=(1.01, 1.), loc=2, frameon=False, prop=lfont, ncol=1, borderaxespad=0., numpoints=1)

		n_HArr = N.logspace(N.log10(1.e-6), N.log10(1.e-1), 100)

		pTda.plot_analytic_curve_altogether(ax, n_HArr, N.array([N.min(m_yArr[jj]), N.max(m_yArr[jj])]), N.median(tempArr_all[ivalid[jj], :], axis=0), N.median(denArr_all[ivalid[jj], :], axis=0), simI["ionizR"])

		ax.set_title(r"$%s > %dM_h$" % (ylabel[1:-1], int(xx)))
		ax.set_xlabel(r"$%s$" % "n_H"); ax.set_ylabel(r"$T$")
		ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		plt.setp(ax.get_xticklabels(), fontsize=fs["tc"])
		plt.setp(ax.get_yticklabels(), fontsize=fs["tc"])
		ax.set_xlim([5.e-6, 1.e-1]); ax.set_ylim([8.e3, 2.e5])
		ax.set_xscale("log"); ax.set_yscale("log")


		figname = "%sGT%dmG_%s_ov%d_vr%d.%s" % (fignamestr[0], int(xx), fignamestr[1], int(simI["over_density"]), int(simI["virial_radius"]), figext)
		plt.savefig(os.path.join(fignamestr[2], figname))

