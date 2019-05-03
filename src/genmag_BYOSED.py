import numpy as np
import os
import optparse
import configparser
import sys

if not hasattr(sys, 'argv'):
	sys.argv  = ['']

required_keys = ['magsmear']

class genmag_BYOSED:
	def __init__(self,PATH_VERSION):
		self.PATH_VERSION = os.path.expandvars('$WFIRST_ROOT/BYOSED_dev/BYOSEDINPUT')
		phase,wave,flux = np.loadtxt('%s/Hsiao07.dat'%self.PATH_VERSION,unpack=True)
		fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
		self.flux = fluxarr
		self.phase = np.unique(phase)
		self.wave = np.unique(wave)
		self.wavelen = len(self.wave)

		self.paramfile = '%s/BYOSED.params'%self.PATH_VERSION
		if os.path.exists(self.paramfile):
			config = configparser.ConfigParser()
			config.read(self.paramfile)
		else: raise RuntimeError('param file %s not found!'%self.paramfile)
		parser = self.add_options(usage='',config=config)

		options,  args = parser.parse_args()
		for k in required_keys:
			if k not in options.__dict__.keys():
				raise RuntimeError('key %s not in parameter file'%k)

		self.options = options

	def add_options(self, parser=None, usage=None, config=None):
			
		if parser == None:
			parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
				
		# The basics
		parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
		parser.add_option('--clobber', default=False, action="store_true",
						  help='clobber output file')
		parser.add_option('--magsmear', default=config.get('main','magsmear'),
						  type="float",help='amount of Gaussian-random mag smearing (default=%default)')

		return parser
			
	def fetchSED_BYOSED(self,tobs,z,maxlam):
		
		iPhase = tobs-self.phase == np.min(np.abs(tobs-self.phase))
		fluxsmear = self.flux[iPhase,:]*10**(0.4*(np.random.normal(0,self.options.magsmear)))
		
		return self.wavelen,self.wave*(1+z),fluxsmear

	def fetchParNames_BYOSED(self):
		pass

