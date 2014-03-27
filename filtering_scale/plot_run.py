import os, sys
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
from matplotlib.font_manager import FontProperties

sys.path.append("..")
import read_outputs as RO
import plot_Td_func as pTdf
import plot_mscale_func as pmf


simI = {
'boxsize': 10., ## simulation boxsize
'colorList': ['b', 'g', 'orange', 'r', 'm', 'purple', 'cyan'], ## to plot
'denlim': 0.05, ## a density limit set by GADGET to make gas form stars
'ids_ref': ["coldIGM", "ids_zc3.49_1.0rv"], ## to get the file name in outputs_idscond. ids are from density condition file or group mass condition file
'iminmax': [None, None], ## if only some ranges of sanpshots were read when picking up the particles from the snapshot, set imin and imax corresponding snap shot numbers
'ionizR': [1], ## 1 or 10 ionizR[0] is to set the ionization rate of the data, if want to plot additional analytic curve, ionizR[1:] is set
'npickptcl': 10000, ## the initial setting of the  number of particles to be picked in simulation. To print out all, it is set to -1
'npltdim': [1, 1], ## nrow, ncol to plot
'nsimptcls': 256.0, ## n of simulation particles (not cubed)
'over_density': 0.0, ## choose overdense particles
'ptcl_pick_cond': 'ids', ## 'grp', 'ids', 'den'
'simbkgd': 'coldIGM', ## gdm10n128: 'newbkgd', others: ''
'zcol': 5.07,
#'fnameE': 'n4-0'
'fnameE': ''
}


#zcolL = [5.07, 3.49, 1.39]; nzcol = len(zcolL)
zcolL = [5.93, 3.07, 1.55]; nzcol = len(zcolL)
#zcolL = [5.07] ; nzcol = len(zcolL)
data = (); datac = (); ivaldata = ()

for zcol in zcolL:
	simI['zcol'] = zcol
	simI['ids_ref'][1] = "ids_zc%.2f_n5-0_1.0rv" % zcol
	#simI['ids_ref'][1] = "ids_zc%.2f_1.0rv" % zcol
	simI['simbkgd'] = 'coldIGM'
	ivalid, oarr, darrc, tarrc, garrc = pmf.remove_earlycollapsed_ptcls_coldIGM(simI, ())
	mjmaxarrc, imaxarrc = pTdf.get_max_jeans_mass_with_cooling(darrc, tarrc, (1., False))

	simI['simbkgd'] = 'bkgd'
	simI['ionizR'][0] = 1
	DC = RO.DataConstruc(**simI)
	izcol, zarr, parr, tarr, garr = DC.get_entire_array1()
	mjmaxarr, imaxarr = pTdf.get_max_jeans_mass_with_cooling(DC.convert_density(), tarr, (1., False))

	## set star particles' densities
	## the order is important
	zs = N.where(parr[:, -1] == 0.)[0]
	parr[zs, -1] = 1.e4
	zs = N.where(parr[:, 1] == 0.)[0]
	parr[zs, -1] = 1.e5

	data += (garr[ivalid], mjmaxarr[ivalid], parr[ivalid, -1])
	ivaldata += (ivalid, )

	print "before taking valid data: %d after taking valid data: %d" % (len(garr), len(ivalid))


	#datac += (darrc, tarrc, garrc, mjmaxarrc, imaxarrc, DC.convert_density(), tarr, garr, mjmaxarr, imaxarr, parr)
	#print "# of particles collapsed after reionization: %d (out of %d)" % (len(ivalid), len(garr))


import plot_mJVSmG as pJG

paths = DC.get_paths_names()
figname = os.path.join(paths[2], "mG-mJmax_%s_ids_n5-0_1.0rv_eqN%d_clsubs_eqNs%d_onlyzcol_estdelta4.6_1.eps" % (simI["ids_ref"][0], True, False))

zcolLtext = ['z=6', 'z=3', 'z=1.5'] ## if it is a zero length list, %.2f real simulation collapsed redshift will be used
## plotting parameters are in order [zcolarr, ionizR, figname, bin_by_grp, pltrangevals, takelog_prop, prange_param, subplt_bins, equal_n_subplt, saperate_star_ptcls] plot_mJEst
param = (zcolL, zcolLtext, simI['ionizR'][0], figname, [8, 300], ([], [8.e8, 2.e11], False), (True, [.5, 5.e3], 5), ([1.e-10, 10., 100., 1.e4, 1.e6], False, False))
#param = (zcolL, zcolLtext, simI['ionizR'][0], figname, [8, 300], ([], [8.e8, 2.e11], False), (True, [.5, 5.e3], 5), ([1.e-10, 10., 100., 1.e4, 1.e6], False, False))
#param = (zcolL, zcolLtext, simI['ionizR'][0], figname, [8, 170], ([8.e7, 2.e11], [8.e8, 2.e11], False), (True, [.5, 5.e3], 5), ([1.e-10, 10., 100., 1.e4, 1.e6], False, False))

## for paper = True option
param = (zcolL, zcolLtext, simI['ionizR'][0], figname, [8, 300], ([1.e8, 2.e12], [2.e9, 1.e12], False), (True, [.5, 5.e3], 5), ([1.e-10, 10., 200., 1.e6], False, False), True)

## no argument paper=False, it automatically executes paper=True
pJG.make_plot_mJ_mh_range_by_den_zred(data, param)

## testing estimate max mJ
"""
import calc_Td_func as cTdf
import plot_Td_analytic as pTda
zcolL = [5.93, 3.07, 1.55]; nzcol = len(zcolL)
n_HArr, tempCCarr = pTda.get_analytic_curve_cc([], 1.)
tempADarr = N.zeros((len(n_HArr), len(zcolL)))
fig = plt.figure(1)
plt.clf()
#plt.subplots_adjust(wspace=0., left=0.1, right=0.95, top=0.8, bottom=0.4)
ax = fig.add_subplot(111)
for i, zcol in enumerate(zcolL):
	tempADarr[:, i] = cTdf.turn_around_adiabat(zcol, n_HArr)
	mJ, nHJ, tJ = cTdf.estimate_max_jeans_mass(zcol)
	ax.plot(n_HArr, tempADarr[:, i], marker='o', label='%.2f' % zcol)
	ax.plot(n_HArr, tempCCarr, color='k', marker='o')
	ax.plot(nHJ, tJ, marker='*')
	ax.text(nHJ, tJ, "%.2e" % mJ)
	print "%.2e, %.2e, %.2e" % (nHJ, tJ, mJ)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1.e-5, 1.e-2); ax.set_ylim(5.e3, 2.e5)
ax.legend(loc=2, frameon=False)
plt.savefig("estimate_maxMJ_test.pdf")
plt.close()
"""
"""
## to test the code
dataL = list(data)
for izcol, zcol in enumerate(zcolL):
	print "%d: %d" % (izcol, len(ivaldata[izcol]))
	pp = datac[izcol*11+10][ivaldata[izcol], :]
	zs = N.where(pp[:, -1] == 0.)[0]
	print len(zs)
	dataL[izcol*3+2][zs] = 1.e4
	zs = N.where(pp[:, 1] == 0.)[0]
	print len(zs)
	dataL[izcol*3+2][zs] = 1.e5

pJG.make_plot_mJ_mh_range_by_den_zred(tuple(dataL), param)
"""


"""
sys.path.append(os.path.join("..", "plot_Td"))
import plot_Td as pTd

figparam = [2, os.path.join(paths[2], "mjmaxTestLowden_%s_Td.pdf" % paths[1]), (3, 3), [1.e-4, 9.e-3], [1.e3, 5.e5], ["coldIGM", "bkgd"], [1.], True]


#for izcol, zcol in enumerate(zcolL):
izcol = 0
idsplt = pJG.get_mJ_mh_range_by_den_ids((data[izcol*3], data[izcol*3+1], data[izcol*3+2]), (param[3], param[4], param[5], param[6]))
den_indic = pmf.choose_particles_by_jeans_mass(data[izcol*3][idsplt], data[izcol*3+1][idsplt], data[izcol*3+2][idsplt], coeff=[4., 4.], p_bins=N.array([1.e-5, 100., 1.e6]), equal_n_in_mbin = False, nbins_mass = 10, nptcl_each = 0)
iL, iH = pmf.get_ids_mj_vs_mh(idsplt, den_indic)
print "Low: %d, High: %d" % (len(iL), len(iH))

darr1 = datac[izcol*11][ivaldata[izcol]]; tarr1 = datac[izcol*11+1][ivaldata[izcol]]; 
garr1= datac[izcol*11+2][ivaldata[izcol]]; mjmaxarr1 = datac[izcol*11+3][ivaldata[izcol]]; 
imaxarr1 = datac[izcol*11+4][ivaldata[izcol]]; 

darr2 = datac[izcol*11+5][ivaldata[izcol]]; tarr2 = datac[izcol*11+6][ivaldata[izcol]]; 
garr2= datac[izcol*11+7][ivaldata[izcol]]; mjmaxarr2 = datac[izcol*11+8][ivaldata[izcol]]; 
imaxarr2 = datac[izcol*11+9][ivaldata[izcol]]; 

pTd.plot_Td_compare((darr1[iL, :], tarr1[iL, :], garr1[iL, :], mjmaxarr1[iL], imaxarr1[iL], darr2[iL, :], tarr2[iL, :], garr2[iL, :], mjmaxarr2[iL, :], imaxarr2[iL]), figparamL) 

pTd.plot_Td_compare((darr1[iH, :], tarr1[iH, :], garr1[iH, :], mjmaxarr1[iH], imaxarr1[iH], darr2[iH, :], tarr2[iH, :], garr2[iH, :], mjmaxarr2[iH, :], imaxarr2[iH]), figparamH)


"""


"""
## plot very low jeans mass particles temperature and density

i0 = N.where(data[1] < 1.e8)[0]
i0_1 = ivalid[i0[::30]]
figparam = (2, os.path.join(paths[2], "lowmaxmJ_%s_Td.pdf" % paths[1]), (3, 3), [1.e-4, 9.e-3], [1.e1, 5.e5], ["coldIGM", "bkgd"], [1.], True)
pTd.plot_Td_compare((datac[0][i0_1], datac[1][i0_1], datac[2][i0_1], datac[3][i0_1], datac[4][i0_1], datac[5][i0_1], datac[6][i0_1], datac[7][i0_1], datac[8][i0_1], datac[9][i0_1]), figparam)


for i in xrange(3):
	#jj = N.argsort(data[i*3+2])
	#kk = N.where(data[i*3+2][jj] > 200.)[0]
	#print kk[0], len(jj)
	#print data[i*3+2][jj[kk[::10]]]
	zz = N.where(datac[i*11+5][:, 1] == 0)[0]
	print len(zz), len(datac[i*11+5])
	print datac[i*11+5][zz[::30], :3], datac[i*11+7][zz[::30]]


"""
