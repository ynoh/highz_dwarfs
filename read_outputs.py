"""
1/24/2013
separate out reading outputs from plot_Td_func
combine construct data array here
get_data is changed

* npickptcls -> changed due to different conditions when sorting out the particles
"""

import os, sys
import numpy as N
import constants as cst
import calc_Td_func as aclc

def get_data_length(fname):
	f = open(fname, 'r')
	line = "#"; ncommentlines = 0
	while line[0] == '#':
		ncommentlines += 1
		line = f.readline()
	f.close()
	ncommentlines -= 1 ## due to one more loop of while

	fw, fr = os.popen2("wc -l %s" % fname)
	nline = int(fr.read().split()[0])
	print "in %s, %d lines in total and %d line(s) is/are for comments." % (fname, nline, ncommentlines)
	return nline - ncommentlines


def read_indv_file(fname, rparam, icolArr):
	"""
	[start, nptcls to pick, skip]
	istart: 0 - some number smaller than skip
	"""
	istart = rparam[0]; nptcls = rparam[1]; skip = rparam[2]

	print icolArr, nptcls, skip
	data = N.zeros((nptcls, len(icolArr)))
	f = open(fname, 'r')
	line = "#";
	while line[0] == '#':
		line = f.readline()

	#print "the first data is: {:}".format(line)
	lline = len(line)

	##to make istart points the first line of the data 
	f.seek(-1*lline, 1)
	for i in xrange(istart):
		f.readline()

	if skip == 1:
		for i, line in enumerate(f):
			data1 = N.array(line.split())
			data1 = data1[icolArr]
			data[i, :] = data1.astype(N.float)
		print "read %d lines" % i
	else:
		j = 0
		for i, line in enumerate(f):
			if i % skip == 0:
				data1 = N.array(line.split())
				data1 = data1[icolArr]
				if j < nptcls:
					data[j, :] = data1.astype(N.float)
					j += 1
				else:
					print "j=%d but possible jmax is %d. no more data will be saved" % (j, nptcls)
		print "read %d lines" % j
		data = data[:j, :]
			
	f.close()

	return data
	

def get_reading_interval(istart, nptcls, nptcls_tot):
	
	if nptcls > 0:
		skip = int(nptcls_tot/nptcls)
		nptcls = (int(nptcls_tot) - istart)/skip + 1
	else:
		skip = 1; nptcls = nptcls_tot
	
	print "get reading interval: istart=%d, nptcls to read =%d among %d, skip=%d" % (istart, nptcls, nptcls_tot, skip)
	return istart, nptcls, skip


def get_partial_data_grps(ptcl_pick_cond, nameparam, ids_ref, rparam, icolArr):
	"""
	by defining different columns, it can be group mass, position and so on	
	"""
	full_data_path = nameparam[0]; collcond = nameparam[1]
	
	if ptcl_pick_cond == "grp":
		fname = os.path.join(full_data_path, "%s.dat" % collcond)
	elif ptcl_pick_cond == "ids":
		id_simbkgd = ids_ref[0]; id_fname = ids_ref[1]
		refpath = os.path.join(os.path.dirname(full_data_path), id_simbkgd)
		fname = os.path.join(refpath, "%s.dat" % id_fname)
	else: pass
	print "To get the group information, reading %s" % (fname)

	nptcls_tot = get_data_length(fname)
	if len(rparam) > 2: 
		print "len(rparam) is longer than 2. given rparam is {:}".format(rparam)
		## when two reading parameters shold be combined i.e. when reading ids' n4-0 partially
		## the second parts in rparam should be the first condition
		## both istarts are meaningless. the first nptcls is used for the second one but no more	
		istart, nptcls, skip1 = get_reading_interval(rparam[2], rparam[3], nptcls_tot) 
		istart, nptcls, skip2 = get_reading_interval(rparam[0], rparam[1], nptcls)
		rparam1 = (rparam[2] + skip1*rparam[0], rparam[1], skip1*skip2)
	else:
		rparam1 = get_reading_interval(rparam[0], rparam[1], nptcls_tot)
		
	data = read_indv_file(fname, rparam1, icolArr)

	return data


def get_partial_data_ptcls_Td(zredArr, nameparam, rparam) :
	"""
	read Td
	"""
	full_data_path = nameparam[0]; collcond = nameparam[1]
	istart = rparam[0]; nptcls = rparam[1]

	## get how many particles in the file
	fname = os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zredArr[0], collcond))
	print fname
	nptcls_tot = get_data_length(fname)
	rparam1 = get_reading_interval(istart, nptcls, nptcls_tot)

	print "every %d particles will be read from %d in total %d" % (rparam1[2], rparam1[0], rparam1[1])
	nred = len(zredArr)

	## set up the arrays
	icolArr = N.array([1, 2])

	for i, zred in enumerate(zredArr):
		fname = os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zred, collcond))
		print "reading %s" % fname
		data = read_indv_file(fname, rparam1, icolArr) 
		if i == 0:
			denArr = N.zeros((len(data), nred)); tempArr = N.zeros_like(denArr); 
		denArr[:, i] = data[:, 0]; tempArr[:, i] = data[:, 1]
	

	return denArr, tempArr


def get_partial_data_ptcls_pos(zredArr, nameparam, rparam) :
	"""
	read positions
	"""
	full_data_path = nameparam[0]; collcond = nameparam[1]
	istart = rparam[0]; nptcls = rparam[1]

	## get how many particles in the file
	fname = os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zredArr[0], collcond))
	nptcls_tot = get_data_length(fname)
	rparam1 = get_reading_interval(istart, nptcls, nptcls_tot)

	print "every %d particles will be read from %d in total %d" % (rparam1[2], rparam1[0], rparam1[1])
	nred = len(zredArr)

	## set up the arrays
	icolArr = N.array([4, 5, 6])

	for i, zred in enumerate(zredArr):
		fname = os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zred, collcond))
		print "reading %s" % fname
		data = read_indv_file(fname, rparam1, icolArr)
		if i == 0:
			posArr = N.zeros((len(data), nred*3)); 
		posArr[:, i*3:(i+1)*3] = data

	return posArr


def get_data_grps(zredArr, nameparam, ids_ref, cols):
	#, high_density=False):
	"""
	read densities and temperatures in the given redshift range 
	from the files in read_data directory
	high_density can take numbers. That indicates how much overdensity we want to consider -> this condition is taken out 8/15
	the columns are increasing number, density, temperature, id, pos[3], grpmss, grppos[3]
	"""

	ptcl_pick_cond = nameparam[0]; full_data_path = nameparam[1]; collcond = nameparam[2]

	## when ptls picked by group condition, read the group information
	grpInfoArr = N.array([])

	if ptcl_pick_cond == "grp":
		filename = "%s.dat" % (collcond)
		print "To get the group information, reading %s/%s" % (full_data_path, filename)
		grpInfoArr = N.loadtxt(os.path.join(full_data_path, filename), usecols=cols)
	elif ptcl_pick_cond == "ids":
		id_simbkgd = ids_ref[0]; id_fname = ids_ref[1]
		refpath = os.path.join(os.path.dirname(full_data_path), id_simbkgd)
		print "To get the group information, reading %s/%s" % (refpath, id_fname)
		grpInfoArr = N.loadtxt(os.path.join(refpath, "%s.dat" % id_fname), usecols=cols)
	else: pass

	return grpInfoArr


def get_data_ptcls_Td(zredArr, nameparam):
	#, high_density=False):
	"""
	read densities and temperatures in the given redshift range 
	from the files in read_data directory
	"""

	ptcl_pick_cond = nameparam[0]; full_data_path = nameparam[1]; collcond = nameparam[2]

	## get the lines of each file (# of lines are the same in each file)
	nptcls = get_data_length(os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zredArr[-1], collcond)))
	nred = len(zredArr)
	print "n of particles: %d" % nptcls

	cols = (1, 2) #cols = (1, 2, 3)
	denArr = N.zeros((nptcls, nred))
	tempArr = N.zeros_like(denArr)

	for i, zred in enumerate(zredArr):
		filename = "Td_%.2f_%s.dat" % (zred, collcond)
		print "reading %s/%s" % (full_data_path, filename)
		data = N.loadtxt(os.path.join(full_data_path, filename), usecols=cols)

		denArr[:, i] = data[:, 0]
		tempArr[:, i] = data[:, 1]
		#print data[::100, 2]

	return denArr, tempArr


def get_data_ptcls_pos(zredArr, nameparam):
	#, high_density=False):
	"""
	read positions of partiles in the given redshift range 
	from the files in read_data directory
	"""
	ptcl_pick_cond = nameparam[0]; full_data_path = nameparam[1]; collcond = nameparam[2]

	## get the lines of each file (# of lines are the same in each file)
	nptcls = get_data_length(os.path.join(full_data_path, "Td_%.2f_%s.dat" % (zredArr[-1], collcond)))
	nred = len(zredArr)
	print "n of particles: %d" % nptcls

	cols = (4, 5, 6)
	posArr = N.zeros((nptcls, nred*3))

	for i, zred in enumerate(zredArr):
		filename = "Td_%.2f_%s.dat" % (zred, collcond)
		print "reading %s/%s" % (full_data_path, filename)
		data = N.loadtxt(os.path.join(full_data_path, filename), usecols=cols)
		posArr[:, 3*i:(i+1)*3] = data

	return posArr


def get_redshift_range(nameparam):
	"""
	find the redshift list corresponding to the snapshots in between imin and imax
	"""
	print nameparam
	path = nameparam[0]; boxsize = nameparam[1]; nsimptcls= nameparam[2]; simbkgd = nameparam[3];
	iminmax = nameparam[4];  

	identifier = "gdm%dn%d" % (int(boxsize), nsimptcls)

	snapshotNList = N.loadtxt(os.path.join(path, "zredList_%s_%s.dat" % (identifier, simbkgd)))

	if iminmax[0] == None:
		iz0 = 0
	else:
		iz0 = N.where(snapshotNList[:, 0] == iminmax[0])[0]
	if iminmax[1] == None:
		iz1 = len(snapshotNList)-1
	else:
		iz1 = N.where(snapshotNList[:, 0] == iminmax[1])[0]
	zredArr = snapshotNList[iz0:iz1+1, 1]
	return zredArr


def get_collcond_name(param):
	"""
	to get the last part of the data file name
	"""
	## currently to use files in old format
	ptcl_pick_cond = param[0]; zcol = param[1]; fname_end = param[2]; idref = param[3]


	if ptcl_pick_cond == "ids":
		return "%s_%s%s" % (idref[0], idref[1], fname_end)
	else:
		return "%s_zc%.2f_%s" % (ptcl_pick_cond, zcol, fname_end)


def get_mainpath():
	path = os.path.abspath(".")
	dirname = "highz_dwarfs"
	ii = path.find(dirname)
	path = path[:ii+len(dirname)]
	return path


def get_subpath(param):
	"""
	simbkgd: newbkgd, newbkgd10, newbkgd_zrei15, or just blanck
	"""
	ptcl_pick_cond = param[0]; boxsize = param[1]; nsimptcls= param[2]; simbkgd = param[3];

	data_dir1 = "outputs"
	data_dir2 = "gdm%dn%d" % (int(boxsize), nsimptcls)
	subpath = os.path.join(data_dir1, data_dir2, simbkgd)
	return subpath


class DataConstruc:
	def __init__(self, **simulationInfo):
		self.ptcl_pick_cond = simulationInfo['ptcl_pick_cond']
		self.boxsize = simulationInfo['boxsize']
		self.nsimptcls = simulationInfo['nsimptcls']
		self.simbkgd = simulationInfo['simbkgd']
		self.zcol = simulationInfo['zcol']
		self.ids_ref = simulationInfo['ids_ref']
		self.iminmax = simulationInfo['iminmax']
		self.ionizR = simulationInfo['ionizR']
		self.fnameE = simulationInfo['fnameE']
		#self.read_pos = simulationInfo['read_pos']

		self.mainpath = get_mainpath() 
		self.collcond = get_collcond_name((self.ptcl_pick_cond, self.zcol, self.fnameE, self.ids_ref))
		self.subpath = get_subpath((self.ptcl_pick_cond, self.boxsize, self.nsimptcls, self.simbkgd))
		#self.zredArr = []; self.izcol = -1  ## random number is assigned to initialize
		self.izcol, self.zredArr = self.get_zred_array()
		self.denArr = []; self.tempArr = []; self.posArr = []; self.grpInfoArr = []

		self.check_consistency()

	
	def check_consistency(self):
		if self.ptcl_pick_cond == "ids":
			ic = self.ids_ref[1].find("zc")
			zcol_ids = float(self.ids_ref[1][ic+2:ic+6])
			if zcol_ids != self.zcol:
				print self.ids_ref[1][ic+2:ic+4]
				print "zcol: %.2f and ids' zcol: %.2f" % (self.zcol, zcol_ids)
				print "they are NOT consistent"
				#raise
		else:
			pass

	def get_paths_names(self):
		return self.mainpath, self.collcond, self.subpath


	def get_zred_array(self):
		## now get the z array upto the collapsed redshift
		zredArr = get_redshift_range((os.path.join(self.mainpath, "read_data"), self.boxsize, self.nsimptcls, self.simbkgd, self.iminmax)) 
		izcol = N.where(zredArr == self.zcol)[0]
		zredArr = zredArr[:izcol+1]
		print self.zcol, izcol, zredArr
		return izcol, zredArr

	def get_entire_array1(self):
		"""
		old way, reading the files which doesn't contain too many particles
		"""
		data_full_path = os.path.join(self.mainpath, "read_data", self.subpath)
		
		## read density
		self.denArr, self.tempArr = get_data_ptcls_Td(self.zredArr, (self.ptcl_pick_cond, data_full_path, self.collcond)) 

		## read group mass
		if self.ptcl_pick_cond == "grp":
			self.grpmArr = get_data_grps(self.zredArr, (self.ptcl_pick_cond, data_full_path, self.collcond), self.ids_ref, N.array([3])) 
		elif self.ptcl_pick_cond == "ids":
			## if I don't want to read everything but read partially
			## ids case, there is only one file containing all the groups' gas particles
			if len(self.fnameE) == 0:
				self.grpmArr = get_data_grps(self.zredArr, (self.ptcl_pick_cond, data_full_path, self.collcond), self.ids_ref, N.array([3])) 
			elif self.fnameE.find("n-1") == -1:
				nameparam = (data_full_path, self.collcond)
				ii = self.fnameE.find("n"); jj = self.fnameE.find("-")
				rparam = (int(self.fnameE[jj+1:]), 10.**int(self.fnameE[ii+1:jj]))
				self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam, N.array([3]))
				self.grpmArr = self.grpmArr.reshape(len(self.grpmArr))
			else:
				self.grpmArr = get_data_grps(self.zredArr, (self.ptcl_pick_cond, data_full_path, self.collcond), self.ids_ref, N.array([3])) 

		return self.izcol, self.zredArr, self.denArr, self.tempArr, self.grpmArr


	def get_entire_array2(self, rparam):
		"""
		read file by selecting lines, rparam=(istart, nptcls)
		"""
		data_full_path = os.path.join(self.mainpath, "read_data", self.subpath)
		nameparam = (data_full_path, self.collcond)

		## read density
		self.denArr, self.tempArr = get_partial_data_ptcls_Td(self.zredArr, nameparam, rparam)

		## read group mass
		if self.ptcl_pick_cond == "ids":
			## files are changed ---- 2013/06/18
			if len(self.fnameE) == 0:
				self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam, N.array([3]))
			elif self.fnameE.find("n-1") == -1:
				## if I don't want to read everything but read partially
				## ids case, there is only one file containing all the groups' gas particles -- this is the old case (2013/06/18)
				## rparam0 and rparam should be combined
				## the data should be read first with rparam0 and then that data is re-read with rparam
				#ii = self.fnameE.find("n"); jj = self.fnameE.find("-")
				#ii = self.collcond.find("n"); jj = self.collcond.find("-")
				#print self.collcond
				#rparam0 = (rparam[0], rparam[1], int(self.collcond[jj+1:]), 10.**int(self.collcond[ii+1:jj])) 
				## e.g. n4-0 case: rparam[2] = 0 rparam[3] = 10000
				rparam0 = (rparam[0], rparam[1], rparam[2], rparam[3])
				self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam0, N.array([3])) 
		elif self.ptcl_pick_cond == "grp":
			self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam, N.array([3]))

		self.grpmArr = self.grpmArr.reshape(len(self.grpmArr))

		return self.izcol, self.zredArr, self.denArr, self.tempArr, self.grpmArr

	def get_array2(self, rparam, icols):
		"""
		read some data only
		"""
		data_full_path = os.path.join(self.mainpath, "read_data", self.subpath)
		nameparam = (data_full_path, self.collcond)

		if self.ptcl_pick_cond == "ids" and self.fnameE.find("n-1") == -1:
			## if I don't want to read everything but read partially
			## ids case, there is only one file containing all the groups' gas particles
			## rparam0 and rparam should be combined
			## the data should be read first with rparam0 and then that data is re-read with rparam
			ii = self.fnameE.find("n"); jj = self.fnameE.find("-")
			rparam0 = (rparam[0], rparam[1], int(self.fnameE[jj+1:]), 10.**int(self.fnameE[ii+1:jj]))
			self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam0, icols) 
		elif self.ptcl_pick_cond == "grp":
			self.grpmArr = get_partial_data_grps(self.ptcl_pick_cond, nameparam, self.ids_ref, rparam, icols)

		self.grpmArr = self.grpmArr.reshape(len(self.grpmArr))


	def convert_density(self):
		self.denArr = self.denArr*aclc.mean_hydrogen_numden_z(self.zredArr)
		return self.denArr


	def choose_propArr(self, propnamew, ipropcol, bins, fixed_interval):
		"""
		to make groups. this function should be here
		"""
		## properties to be used for grouping
		pArr_ref = 1.
		if propnamew == "mG":
			propArr = self.grpInfoArr_all[:, 2]; propname = "M_G"
		elif propnamew == "nH":
			propArr = self.denArr_all[:, ipropcol]; propname = "n_H(%.1f)" % self.zredArr_all[ipropcol]
		elif propnamew == "del":
			propArr = self.denArr_all[:, ipropcol]; propname = "\delta(%.1f)" % self.zredArr_all[ipropcol];
		elif propnamew.find("mJ")!= -1:
			pArr_ref = self.grpInfoArr_all[:, 2];  ##????
			if propnamew.find("peak") != -1:
				propArr = pTdf.get_peak_jeans_mass()
			elif propnamew.find("max") != -1:
				propArr = pTdf.get_max_jeans_mass()
			else:
				if ipropcol >= 0:
					propArr = aclc.jeans_mass(self.tempArr_all[:, ipropcol], self.denArr_all[:, ipropcol]); 
					propname = "M_J(%.1f)"
				else:
					print "ipropcol=%d is not valid" % ipropcol

				"""
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
				"""
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
