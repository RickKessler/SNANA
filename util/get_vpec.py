#!/usr/bin/env python
# D. Jones - 4/13/20
"""Peculiar velocity corrections from a 2M++ flow model.  

   Usage: get_vpec.py ra dec redshift

Beware that command line usage may be slow due to reading
large map; get_vpec.main(ra,dec,zcmb) should be faster.

Explanation adopted from Jones+18:
This code corrects for peculiar velocities using the nearby galaxy 
density field measured by the 2M++ catalog from 2MASS 
(Lavaux & Hudson 2011, arxiv:1105.6107). The uncorrelated 
uncertainty associated with each VPEC value is +_250 km/s 
(from Pantheon analysis in Scolnic+18).

The peculiar-velocity model is parameterized by 
        beta_I = Omega_M^0.55/b_I, 
where b_I describes the light-to-matter bias. 
(beta_I is unrelated to the SALT2 nuisance parameter)
Carrick et al. 2015 (arXiv:1504.04627) measured beta_I=0.43+-0.021,
which is propagated here as a systematic error on vpec.
However, the 250km/s error from Pantheon is not included here.

This algorithm and map, along with an older version of this code, 
has been used previously in the following SNIa-cosmology analyses:

https://ui.adsabs.harvard.edu/abs/2018ApJ...859..101S
   (Pantheon analysism, Scolnic et al., 2018)
https://ui.adsabs.harvard.edu/abs/2018ApJ...857...51J
   (Photometric SNIa-cosmology results, Jones 2018)
https://ui.adsabs.harvard.edu/abs/2019ApJ...881...19J
   (SNIa-cosmologuy from single telescope, Jones 2019)

"""


import argparse
import os
import numpy as np

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

# a few hard-coded variables that might need to be updated someday

_c=299792458/1000.0
_beta = 0.43
_beta_err = 0.021*5.0  # 5 sigma is industry standard
_v_helio = 371.
_l_h = 263.85 # galactic longitude
_b_h = 48.25 # galactic latitude

_v_LG = 318.
_l_LG = 105.7467
_b_LG = -5.9451
_l_0 = np.radians(264.14)
_b_0 = np.radians(48.26)

_vpecerr = 250.0  # only for printout

# define hard-wired default map file name
SNDATA_ROOT   = os.environ['SNDATA_ROOT']
vpec_mapfile_default = ('%s/models/VPEC/twomass++_velocity_LH11.npy' % SNDATA_ROOT)

# ====================================================

def sin_d(x):
	""" sine in degrees for convenience """
	return np.sin(x*np.pi/180.)

def cos_d(x):
	""" cosine in degrees, for convenience """
	return np.cos(x*np.pi/180.)

def dmdz(z):
	""" converts uncertainty in redshift to uncertainty in magnitude, for clarity"""
	return 5./(np.log(10)*z)

class VelocityCorrection(object):
	""" functions to correct for peculiar velocities using 2M++ velocity field (Carrick+ 2015) """

	def __init__(self, velfield):
		""" sets up constants and velocity field; includes:
		params of CMB dipole w.r.t. heliocentric frame (Bennett+ 2003)
		params of LG w.r.t. heliocentric frame (Tully+ 2008 table 3)
		all velocities in km/s; angles in degrees
		sets up 2M++ velocity field is in galactic cartesian coordinates,
		in local group frame and in Mpc/h
		this is linear with coefficient beta
		we perturb beta to measure the effect to compute the covariance matrix
		"""

		self.c = _c

		self.v_helio = _v_helio
		self.l_h = _l_h # galactic longitude
		self.b_h = _b_h # galactic latitude

		self.v_LG = _v_LG
		self.l_LG = _l_LG
		self.b_LG = _b_LG

		self.velocity_field = np.load(velfield)
		self.beta = _beta
		self.beta_err = _beta_err
		self.r_plus = 1 + self.beta_err/self.beta
		self.r_minus = 1 - self.beta_err/self.beta

	def convert_helio_to_LG(self, v, l, b):
		""" converts velocity from heliocentric frame to Local Group frame """

		return v - self.v_LG*(sin_d(b)*sin_d(self.b_LG)
							  + cos_d(b)*cos_d(self.b_LG)*cos_d(l-self.l_LG))
	
	def convert_to_helio(self, l, b, z):
		c_icrs = SkyCoord(ra=l*u.degree, dec=b*u.degree, frame = 'icrs')
		c_icrs = c_icrs.galactic
		b = c_icrs.b.degree
		l = c_icrs.l.degree

		b = np.radians(b)
		l = np.radians(l)

		v = (z*_c - _v_helio * (np.sin(b) * np.sin(_b_0)+np.cos(b) * np.cos(_b_0) *
								np.cos(l-_l_0)))/(_c)
		return v
																					
	
	def lookup_velocity(self, z_h, l, b):
		""" look up velocity field to find peculiar velocity vector
		input needs to be heliocentric position (and in degrees)
		output is vector in galactic cartesian coordinates, in CMB frame
		"""
		cz_LG = self.convert_helio_to_LG(self.c*z_h, l, b)
		x = cz_LG * cos_d(l) * cos_d(b)
		y = cz_LG * sin_d(l) * cos_d(b)
		z = cz_LG * sin_d(b)
		x, y, z = x/100, y/100, z/100 # to convert from km/s to Mpc/h
		i = int(round(128. + 256./400*x))
		j = int(round(128. + 256./400*y))
		k = int(round(128. + 256./400*z))
		try:
			vpec = self.velocity_field[:, i, j, k]
		except IndexError:
			vpec = np.array([0, 0, 0])#outside velocity field; approximate as zero
		return vpec

	def correct_redshift(self, z_h, vpec, l, b):
		""" convert helio redshift to cosmological redshift (zbar; in CMB frame)
		input needs to be vector in galactic cartesian coordinates, w.r.t. CMB frame
		components are: heliocentric motion, peculiar velocity (in radial direction)
		"""
		helio_corr = self.v_helio/self.c*((sin_d(b)*sin_d(self.b_h)
										   + cos_d(b)*cos_d(self.b_h)*cos_d(l-self.l_h)))
		pec_corr = vpec.dot(np.array([cos_d(l)*cos_d(b),
									  sin_d(l)*cos_d(b), sin_d(b)]))/self.c
		corr_term = 1 - helio_corr + pec_corr
		return (1+z_h)/corr_term - 1

	def covmat_pecvel(self, SNe):
		""" build covariance matrix """
		z, z_err = self.apply(SNe)
		nSNe = len(SNe)
		cpecvel = np.zeros((3*nSNe, 3*nSNe))

		for i in range(nSNe):
			for j in range(i, nSNe):
				cpecvel[3*i, 3*j] = dmdz(z[i])*dmdz(z[j]) * z_err[i] * z_err[j]
				cpecvel[3*j, 3*i] = cpecvel[3*i, 3*j]

		return cpecvel


# xxx RK mark delete vpec_mapfile_default = os.path.expandvars('$SNDATA_ROOT/models/VPEC/twomass++_velocity_LH11.npy')

if os.path.exists(vpec_mapfile_default): 
	_ini = VelocityCorrection(vpec_mapfile_default)
else:
	print('warning: default map file %s cannot be found' 
	      % vpec_mapfile_default)
	_ini = None

	
def main(ra,dec,z,vpec_mapfile=None):

	if vpec_mapfile is not None and os.path.expandvars(vpec_mapfile) != vpec_mapfile_default:
		ini = VelocityCorrection(vpec_mapfile)
	else: ini = _ini
	if _ini is None: raise RuntimeError('cannot find flow map file')
	
	sc = SkyCoord(ra,dec,unit=(u.deg, u.deg))
	gsc = sc.galactic

	z_hel    = ini.convert_to_helio(sc.ra.degree, sc.dec.degree, z)
	vpec     = ini.lookup_velocity(z_hel,gsc.l.degree,gsc.b.degree)
	pec_corr = ini.correct_redshift(z_hel,vpec,gsc.l.degree,gsc.b.degree)

	z_c = ini.correct_redshift(z_hel,vpec,gsc.l.degree,gsc.b.degree)

	ini.beta     = _beta
	ini.beta_err = _beta_err
	ini.r_plus   = 1 + ini.beta_err/ini.beta
	ini.r_minus  = 1 - ini.beta_err/ini.beta
	z_plus = ini.correct_redshift(z_hel, ini.r_plus*vpec,gsc.l.degree,gsc.b.degree)
	z_minus = ini.correct_redshift(z_hel, ini.r_minus*vpec,gsc.l.degree,gsc.b.degree)

	vpec0 = (pec_corr-z)*_c
	vpec1 = (z_plus-z)*_c  # includes error
	vsys  = abs(vpec0 - vpec1)
# xxx mark delete RK	return (pec_corr-z)*_c,(z_plus-z)*_c
	return vpec0,vsys

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("ra",type=float)
	parser.add_argument("dec",type=float)
	parser.add_argument("redshift",type=float)
# xxx RK mark delete  parser.add_argument("--vpec_mapfile",default="$SNDATA_ROOT/models/VPEC/twomass++_velocity_LH11.npy",type=str)
	parser.add_argument("--vpec_mapfile",default=vpec_mapfile_default,type=str)
	p = parser.parse_args()
	
	vpec,vsyserr = main(p.ra,p.dec,p.redshift,vpec_mapfile=p.vpec_mapfile)

	print('vpec = %.1f +- %.1f +- %.1f(5sigma_sys)  km/sec' 
	      % (vpec, _vpecerr, vsyserr) )
# xxx mark delete RK	print('VPEC VPEC_SYS')
	print(vpec,vsyserr)
