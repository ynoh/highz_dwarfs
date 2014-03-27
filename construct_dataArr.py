"""
there are common parts in most of the codes to construct data array

proparr is now getting from *_all

"""
import os, sys
import numpy as N
import plot_Td_func as pTd
sys.path.append(os.path.join("..", "analytic"))
import constants as cst
reload(cst)
import calc_Td_func as aclc
reload(aclc)

"""
simulationInfo = {'boxsize': 10.0,
 'colorList': ['b', 'g', 'orange', 'r', 'm', 'purple', 'cyan'],
 'denlim': 0.05,
 'discard_zeros': True,
 'grpMrange': [0.01, 20.0],
 'ids_based_on_grp': True,
 'iminmax': [None, None],
 'ionizR': [1],
 'npickptcl': 1000,
 'npltdim': [3, 3],
 'nsimptcls': 128.0,
 'over_density': 200.0,
 'plot_upto_zcol': True,
 'ptcl_pick_cond': 'grp',
 'simbkgd': 'newbkgd',
 'convert_density': True,
 'virial_radius': 0.0,
 'zcol': 3.5}
"""


class DataConstruc:

	def __init__(self, **simulationInfo):
		self.ptcl_pick_cond = simulationInfo['ptcl_pick_cond']
		self.boxsize = simulationInfo['boxsize']
		self.nsimptcls = simulationInfo['nsimptcls']
		self.simbkgd = simulationInfo['simbkgd']
		self.zcol = simulationInfo['zcol']
		self.grpMrange = simulationInfo['grpMrange']
		self.ids_ref = simulationInfo['ids_ref']
		self.iminmax = simulationInfo['iminmax']
		self.ionizR = simulationInfo['ionizR']
		self.npickptcl = simulationInfo['npickptcl']

		self.discard_zeros = simulationInfo['discard_zeros']
		self.convert_density = simulationInfo['convert_density']
		self.plot_upto_zcol = simulationInfo['plot_upto_zcol']
		self.npltdim = simulationInfo['npltdim']
		self.colorList = simulationInfo['colorList']
		self.denlim = simulationInfo['denlim']
		self.over_density = simulationInfo['over_density']
		self.virial_radius = simulationInfo['virial_radius']
		self.izcol = -1 ## random number is assigned to initialize
		self.collcond = pTd.get_collcond_name(self.ptcl_pick_cond, zcol=self.zcol, grpMrange=self.grpMrange, ids_ref=self.ids_ref, npick=self.npickptcl)

		self.data_abspath, self.subpath = pTd.get_data_paths(self.ptcl_pick_cond, self.boxsize, self.nsimptcls, self.simbkgd)
		self.zredArr = []; self.izcol = 0
		self.denArr = []; self.tempArr = []; self.posArr = []; self.grpInfoArr = []


	def get_paths_names(self):
		return self.collcond, self.data_abspath, self.subpath

	def get_entire_array(self):
		data_full_path = os.path.join(self.data_abspath, self.subpath)
		
		## redshift array
		self.zredArr_all = pTd.get_redshift_range(self.data_abspath, self.boxsize, self.nsimptcls, imin=self.iminmax[0], imax=self.iminmax[1], simbkgd=self.simbkgd)
		self.izcol = N.where(self.zredArr_all == self.zcol)[0]
		print "collapsed redshift: {:} and its snapshot number: {:}".format(self.zcol, self.izcol)

		## read data
		if self.plot_upto_zcol:
			self.zredArr_all = self.zredArr_all[:self.izcol+1]
			self.denArr_all, self.tempArr_all, self.posArr_all, self.grpInfoArr_all = pTd.get_data(self.zredArr_all, self.ptcl_pick_cond, data_full_path, self.collcond, ids_ref=self.ids_ref, density_convert = self.convert_density, discard_zeros = self.discard_zeros)
		else:
			self.denArr_all, self.tempArr_all, self.posArr_all, self.grpInfoArr_all = pTd.get_data(self.zredArr_all, self.ptcl_pick_cond, data_full_path, self.collcond, density_convert = self.convert_density, discard_zeros = self.discard_zeros)
		return self.izcol, self.zredArr_all, self.denArr_all, self.tempArr_all, self.posArr_all, self.grpInfoArr_all



	def process_data(self, zredArr_all=[], denArr_all=[], tempArr_all=[], posArr_all=[], grpInfoArr_all=[]):

		nptcls = len(denArr_all)

		if nptcls == 0:
			self.izcol, self.zredArr_all, self.denArr_all, self.tempArr_all, self.posArr_all, self.grpInfoArr_all = self.get_entire_array()
			nptcls = len(self.denArr_all)
		else:
			self.izcol = N.where(zredArr_all == self.zcol)[0] ## in case get_entire array wasn't called
			self.zredArr_all, self.denArr_all, self.tempArr_all, self.posArr_all, self.grpInfoArr_all = zredArr_all, denArr_all, tempArr_all, posArr_all, grpInfoArr_all

		## take over density particles only?
		if self.over_density > 0.:
			print "getting high density particles only: %.1f" % self.over_density
			idense = pTd.get_high_density_ptcls(self.over_density, denArr_all[:, self.izcol], self.zcol)
		else:
			print "no over density criterion"
			idense = N.arange(nptcls)

		## take the only particles inside the virial radius
		## if it is 0, virial_radius varies according to the given equaiton 
		if self.virial_radius >= 0.:
			print "getting particles inside of the virial radius %.1f" % self.virial_radius
			iinside = pTd.get_ptcls_inside_virial_radius(posArr_all[:, self.izcol*3:(self.izcol+1)*3], grpInfoArr_all, self.zcol, self.boxsize*1.e3, self.virial_radius)
		else:
			print "no virial radius condition"
			iinside = N.arange(nptcls)


		if len(idense) < nptcls and len(iinside) < nptcls:
			import collections
			iall = N.hstack((idense, iinside))
			iall_count = collections.Counter(iall)
			iall_intersect = [i for i in iall_count if iall_count[i] > 1]
		elif len(idense) == nptcls and len(iinside) == nptcls:
			iall_intersect = N.arange(nptcls)
		else:
			if len(idense) <= len(iinside):
				iall_intersect = idense
			else:
				iall_intersect = iinside
		
		#self.zredArr = zredArr_all[:self.izcol+1]
		#self.denArr = denArr_all[iall_intersect, :self.izcol+1]; self.tempArr = tempArr_all[iall_intersect, :self.izcol+1]
		#self.posArr = posArr_all[iall_intersect, :self.izcol+1]; self.grpInfoArr = grpInfoArr_all[iall_intersect, :self.izcol+1]

		return N.array(iall_intersect), self.izcol
		#return iall_intersect, self.izcol, zredArr_all, denArr_all, tempArr_all, posArr_all, grpInfoArr_all
		#return iall_intersect, self.zredArr, self.denArr, self.tempArr, self.posArr, self.grpInfoArr


	def choose_propArr(self, propnamew, ipropcol, bins, fixed_interval):
		"""
		to make groups. this function should be here
		"""
		## properties to be used for grouping
		pArr_ref = 1.
		if propnamew == "mgrp":
			propArr = self.grpInfoArr_all[:, 2]; propname = "M_G"
		elif propnamew.find("den")!= -1:
			if ipropcol > self.denArr_all.shape[1] or ipropcol < -1:
				ipropcol = self.izcol
			if self.convert_density:
				propArr = self.denArr_all[:, ipropcol]; propname = "n_H(%.1f)" % self.zredArr_all[ipropcol];
			else:
				propArr = self.denArr_all[:, ipropcol]; propname = "\delta(%.1f)" % self.zredArr_all[ipropcol];
		elif propnamew.find("mJ")!= -1:
			pArr_ref = self.grpInfoArr_all[:, 2]; propname = "X"
			if ipropcol >= 0:
				propArr = aclc.jeans_mass(tempArr_all[:, ipropcol], denArr_all[:, ipropcol])
			else:
				m_JArr = N.zeros_like(denArr_all)
				for i in xrange(denArr_all.shape[1]):
					m_JArr[:, i] = aclc.jeans_mass(tempArr_all[:, i], denArr_all[:, i])
				propArr = N.max(m_JArr, axis=1)
				## reassign max(m_JArr) if den > denlim
				denseptcls = pTd.find_too_dense_ptcls(denArr_all, denlim=denlim)
				for i in xrange(len(denseptcls)):
					irow = denseptcls[i, 0]; icol = denseptcls[i, 1]
					if icol > 0:
						propArr[irow] = N.max(m_JArr[irow, :ipropcol])
					else:
						print "mJmax-density check:", denArr_all[irow, :]
						propArr[irow] = m_JArr[irow, ipropcol]
		else:
			sys.exit("invalid propnamew: %s" % propnamew)

		if len(bins) > 0:
			csave_bins = "G" ##given bins
			ngrp = len(bins) - 1
		else:
			if fixed_interval:
				csave_bins = "F" ## fixed interval
			else:
				csave_bins = "S" ## sorted
		return propnamew, propname, propArr, csave_bins, pArr_ref
