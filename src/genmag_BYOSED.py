#!/usr/bin/env python

import numpy as np	
import os
import optparse
import configparser
import sys
from scipy.interpolate import RectBivariateSpline
from ast import literal_eval
from scipy.stats import rv_continuous,norm as normal
from copy import copy

if not hasattr(sys, 'argv'):
		sys.argv  = ['']

required_keys = ['magsmear']

class genmag_BYOSED:
		def __init__(self,PATH_VERSION,OPTMASK):
			self.PATH_VERSION = os.path.dirname(PATH_VERSION)
			phase,wave,flux = np.loadtxt('%s/Hsiao07.dat'%self.PATH_VERSION,unpack=True)
			fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
			self.flux = fluxarr*10**(0.4*19.365)
			self.phase = np.unique(phase)
			self.wave = np.unique(wave)
			self.wavelen = len(self.wave)
			self.sedInterp=RectBivariateSpline(self.phase,self.wave,self.flux,kx=3,ky=3)
				
			self.paramfile = '%s/BYOSED.params'%self.PATH_VERSION
			self.warpfile='%s/sed_warps.dat'%self.PATH_VERSION
				
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

			self.warping_types=['stretch','color','dust']
			self.warp_bools={k:self.fetchParVals_BYOSED(k) for k in self.warping_types}
			self.warping_functions,self.warping_distributions,self.warping_params=self.fetchWarp_BYOSED()
			#self.applyWarping_Effects()
				
		def add_options(self, parser=None, usage=None, config=None):
						
				if parser == None:
						 parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
				# The basics
				parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
				parser.add_option('--clobber', default=False, action="store_true",help='clobber output file')
				parser.add_option('--magsmear', default=config.get('MAIN','magsmear'),
								  type="float",help='amount of Gaussian-random mag smearing (default=%default)')
				parser.add_option('--magoff', default=config.get('MAIN','magoff'),
								  type="float",help='mag offset (default=%default)')
				#parser.add_option('--stretch',action="store_true",
				#				  default=literal_eval(config.get('MAIN','stretch')),
				#				  help="include stretch function? (default=%default)")
				#parser.add_option('--color',action="store_true",
				#				  default=literal_eval(config.get('MAIN','color')),
				#				  help="include intrinsic color-law function? (default=%default)")
				#parser.add_option('--dust',action="store_true",
				#				  default=literal_eval(config.get('MAIN','dust')),
				#				  help="include dust function? (default=%default)")

				return parser
				
				
		
		def fetchWarp_BYOSED(self):
			if os.path.exists(self.warpfile):
				config=configparser.ConfigParser()
				config.read(self.warpfile)
			else: raise RuntimeError('warping file %s not found!'%self.warpfile)
				
									  
			#read in warp file
			func_dict=dict([])
			dist_dict=dict([])
			param_dict=dict([])

			for warp in self.warping_types:
				dist_data=literal_eval(config.get('distributions',warp))
				dist_dict[warp]=_skewed_normal(warp,dist_data)
				if self.warp_bools[warp]:
					func_data=np.array(literal_eval(config.get('functions',warp)))
					x0,x1,y=_read_griddata(func_data)
						
					func_dict[warp]=RectBivariateSpline(x0,x1,y,kx=3,ky=3)
					param_dict[warp]=dist_dict[warp]()[0]
				else:
					func_dict[warp]=RectBivariateSpline(self.phase,self.wave,
														np.ones((len(self.phase),len(self.wave))))
					param_dict[warp]=0.
						
			return(func_dict,dist_dict,param_dict)
		

		def updateWarping_Params(self):
				for warp in self.warping_types:
						if self.warp_bools[warp]:
								self.warping_params[warp]=self.warping_distributions[warp]()
						else:
								self.warping_params[warp]=0.
				return self
						
		def applyWarping_Effects(self):
				tempSED=None
				for warp in self.warping_types:
						if tempSED is None:
								tempSED=copy(self.flux)
						tempSED*=self.warping_params[warp]*\
								  self.warping_functions[warp](self.phase,self.wave)
				self.BYOFLUX_Warped=tempSED
				return self
				
		def fetchSED_NLAM(self):
				return self.wavelen

		def fetchSED_LAM(self):

				return list(self.wave)
		
		def fetchSED_BYOSED(self,trest,maxlam,external_id,new_event):
				if len(self.wave)>maxlam:
						raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))
				#iPhase = np.where(np.abs(trest-self.phase) == np.min(np.abs(trest-self.phase)))[0][0]

				fluxsmear = self.sedInterp(trest,self.wave)*10**(0.4*(np.random.normal(0,self.options.magsmear)))
				for warp in self.warping_types:
						if self.warp_bools[warp]:
								fluxsmear*=10**(self.warping_params[warp]*\
										  self.warping_functions[warp](trest,self.wave)[0,:])

				return list(fluxsmear) 
				
		def fetchParNames_BYOSED(self):
				parnames = []
				for k in self.options.__dict__.keys():
						if k != 'clobber' and k != 'verbose':
								parnames += [k]
				return parnames

		def fetchNParNames_BYOSED(self):
				parnames = []
				for k in self.options.__dict__.keys():
						if k != 'clobber' and k != 'verbose':
								parnames += [k]
				return len(parnames)
		
		def fetchParVals_BYOSED(self,parName):
				if parName in self.options.__dict__.keys():
						return self.options.__dict__[parName]
				else:
						return -1

def _read_griddata(gridded):	  
		x0=np.sort(np.unique([x[0] for x in gridded]))
		x1=np.sort(np.unique([x[1] for x in gridded]))
		
		y=np.array([x[2] for x in gridded]).reshape((len(x0),len(x1)))
		return x0,x1,y


class skewed_normal(rv_continuous):
		"Skewed Normal Distribution"
		def _pdf(self,x,mu,left_sigma,right_sigma):
				try:
						mu=list(mu)[0]
						left_sigma=list(left_sigma)[0]
						right_sigma=list(right_sigma)[0]
				except:
						pass
				
				left=normal(loc=mu,scale=left_sigma)
				right=normal(loc=mu,scale=right_sigma)
				return(np.piecewise(x,[x<mu,x>=mu],
									[lambda y : left.pdf(y)/np.max(left.pdf(y)),
									 lambda y : right.pdf(y)/np.max(right.pdf(y))]))
		def _argcheck(self,*args):
				return True


def _skewed_normal(name,dist_dat):
		dist = skewed_normal(name,a=dist_dat['lower_lim'],b=dist_dat['upper_lim'])
		return(lambda : dist.rvs(dist_dat['mu'],dist_dat['sigma_left'],dist_dat['sigma_right'],size=1))
		

def main():
		mySED=genmag_BYOSED('../../BYOSEDINPUT/')
		print(np.max(mySED.fetchSED_BYOSED(0,5000)))
		#import matplotlib.pyplot as plt
		#fig=plt.figure()
		#ax=fig.gca()
		#for i in range(5):
		#		 temp=mySED.fetchSED_BYOSED(0,5000)
		#		 ax.plot(mySED.wave,temp/np.max(temp))
		#		 mySED.updateWarping_Params()
		#plt.show()


if __name__=='__main__':
		main()
