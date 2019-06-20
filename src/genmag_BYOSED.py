import numpy as np  
import os,six,abc
import optparse
import configparser
import pandas
import sys
from scipy.interpolate import RectBivariateSpline,interp2d,interpn
from ast import literal_eval
from scipy.stats import rv_continuous,gaussian_kde,norm as normal
from copy import copy
import pickle

if not hasattr(sys, 'argv'):
		sys.argv  = ['']

required_keys = ['magsmear']


__mask_bit_locations__={'verbose':1,'dump':2}


class genmag_BYOSED:

		def __init__(self,PATH_VERSION,OPTMASK,ARGLIST,HOST_PARAM_NAMES):

			self.verbose = OPTMASK & (1 << __mask_bit_locations__['verbose']) > 0

				self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))

				self.NAMES_HOSTPAR = NAMES_HOSTPAR.split(',')
				self.PATH_VERSION = os.path.dirname(PATH_VERSION)

=======
			self.dump = OPTMASK & (1 << __mask_bit_locations__['dump'])>0
			self.sn_id=None

			self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))
			self.host_param_names=HOST_PARAM_NAMES

			self.paramfile = os.path.join(self.PATH_VERSION,'BYOSED.params')
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

			self.warp_effects=self.fetchParVals_BYOSED(config)
			
			self.sn_effects,self.host_effects=self.fetchWarp_BYOSED(config)


			phase,wave,flux = np.loadtxt(os.path.join(self.PATH_VERSION,self.options.sed_file),unpack=True)


			fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
			self.x0=10**(.4*19.365)
			self.flux = fluxarr*self.x0

			self.phase = np.unique(phase)
			self.wave = np.unique(wave)
			self.wavelen = len(self.wave)

			self.sedInterp=interp2d(self.phase,self.wave,self.flux.T,kind='linear',bounds_error=True)

	
			return
		
		def add_options(self, parser=None, usage=None, config=None):
						
				if parser == None:
						 parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
				# The basics
				parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
				parser.add_option('--clobber', default=False, action="store_true",help='clobber output file')
				parser.add_option('--magsmear', default=config.get('MAIN','MAGSMEAR'),
								  type="float",help='amount of Gaussian-random mag smearing (default=%default)')
				parser.add_option('--magoff', default=config.get('MAIN','MAGOFF'),
								  type="float",help='mag offset (default=%default)')
				parser.add_option('--sed_file',default=config.get('MAIN','SED_FILE'),
								  type='str',help='Name of sed file')
				
				
				return parser
				
				
		
		def fetchWarp_BYOSED(self,config):
				#read in warp effect data
				#if speed becomes a real issue, we could
				#combine the various additive functions
				#into a single function, which would have
				#to be repeated for every SN but would
				#save looping through all for various 
				#phases...not sure if that would actually
				#come out on top though.

				sn_dict=dict([])
				host_dict=dict([])

				
				for warp in self.warp_effects:
					warp_data={k.upper():np.array(config.get(warp,k).split()).astype(float) if k.upper() not in ['SN_FUNCTION','HOST_FUNCTION','DIST_FILE'] else config.get(warp,k) for k in config[warp]}
					distribution=_get_distribution(warp,warp_data,self.PATH_VERSION)
					
						
					if 'SN_FUNCTION' in warp_data:
						
								
						if 'SN_FUNCTION_SCALE' not in warp_data:
							scale_factor=1.
						else:
							scale_factor=warp_data['SN_FUNCTION_SCALE']
						try:
							sn_param_names,sn_function=_read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['SN_FUNCTION']))),scale_factor)
						except RuntimeError:
							raise RuntimeError("Do not recognize format of function for %s SN Function"%warp)
						warp_parameter=distribution()[0]

						sn_dict[warp]=warpModel(warp_function=sn_function,
												param_names=sn_param_names,
												parameters=np.array([0 if sn_param_names[i]!=warp.upper() else warp_parameter for i in range(len(sn_param_names))]),
												warp_parameter=warp_parameter,
												warp_distribution=distribution,
												name=warp)
						
					if 'HOST_FUNCTION' in warp_data:

						try:
							host_param_names,host_function=_read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['HOST_FUNCTION']))))
						except RuntimeError:
							raise RuntimeError("Do not recognize format of function for %s HOST Function"%warp)
						warp_parameter=distribution()[0]
						host_dict[warp]=warpModel(warp_function=host_function,
												param_names=host_param_names,
												parameters=np.array([0 if host_param_names[i]!=warp.upper() else warp_parameter for i in range(len(host_param_names))]),
												warp_parameter=warp_parameter,
												warp_distribution=distribution,
												name=warp)
								
				return(sn_dict,host_dict)
		

		#def updateWarping_Params(self):
		#       for warp in self.warp_effects:
		#           if warp in sn_effects.keys():
		#               self.sn_effects[warp].update(warp_parameter
		#       return self
						
				
		def fetchSED_NLAM(self):
				return self.wavelen

		def fetchSED_LAM(self):

				return list(self.wave)
		

		def fetchSED_BYOSED(self,trest,maxlam,external_id,new_event,hostpars):

				if len(self.wave)>maxlam:
						raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))
				#iPhase = np.where(np.abs(trest-self.phase) == np.min(np.abs(trest-self.phase)))[0][0]
				if self.sn_id is None:
						self.sn_id=external_id
				fluxsmear=self.sedInterp(trest,self.wave).flatten()
				orig_fluxsmear=copy(fluxsmear)
				if self.options.magsmear!=0.0:
						fluxsmear *= 10**(0.4*(np.random.normal(0,self.options.magsmear)))
				trest_arr=trest*np.ones(len(self.wave))
				for warp in [x for x in self.warp_effects if x!='COLOR']:
					if True:
						if external_id!=self.sn_id:
							if warp in self.sn_effects.keys():
								self.sn_effects[warp].updateWarp_Param()
								if warp in self.host_effects.keys():
									self.host_effects[warp].warp_parameter=self.sn_effects[warp].warp_parameter
								else:
									self.host_effects[warp].updateWarp_Param()
							self.sn_id=external_id
							
						
						# not sure about the multiplication by x0 here, depends on if SNANA is messing with the 
						# absolute magnitude somewhere else
						product=1.
						temp_warp_param = None
						if warp in self.sn_effects.keys():
							if self.verbose:
								print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].warp_parameter))
							product*=self.sn_effects[warp].flux(trest_arr,self.wave,HOST_PARAMS,self.host_param_names)

							if warp in self.sn_effects[warp]._param_names:
								temp_warp_param=1.
							else:
								temp_warp_param=self.sn_effects[warp].warp_parameter
						if warp in self.host_effects.keys():
							product*=self.host_effects[warp].flux(trest_arr,self.wave,HOST_PARAMS,self.host_param_names)
							if temp_warp_param is None:
								if warp in self.host_effects[warp]._param_names:
									temp_warp_param=1.
								else:
									temp_warp_param=self.host_effects[warp].warp_parameter

						fluxsmear+=temp_warp_param*product*self.x0
					#except:
						#import pdb; pdb.set_trace()

				if 'COLOR' in self.warp_effects:
						if external_id!=self.sn_id:
								self.sn_effects['COLOR'].updateWarping_Params()
								self.sn_id=external_id
						if self.verbose:
								print('Phase=%.1f, %s: %.2f'%(trest,'COLOR',self.sn_effects['COLOR'].warp_parameter))
						fluxsmear*=10**(-0.4*self.sn_effects['COLOR'].warp_parameter*\
								self.sn_effects['COLOR'].flux(trest_arr,self.wave,HOST_PARAMS).flatten())
							  
				return list(fluxsmear) 
				
		def fetchParNames_BYOSED(self):
				parnames = []
				for k in self.warping_params.keys():
					parnames += [k]
				return parnames

		def fetchNParNames_BYOSED(self):
				parnames = []
				for k in self.warping_params.keys():
					parnames += [k]
				return len(parnames)

		def fetchParVals_BYOSED_4SNANA(self,varname):
				return self.warping_params[varname]

		def fetchParVals_BYOSED(self,config):
				if 'FLAGS' in config.sections():
						return([k.upper() for k in list(config['FLAGS'].keys()) if config['FLAGS'][k]=='True'])
				else:
						return([x for x in config.sections() if x not in ['MAIN','FLAGS']])



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
				pdf=np.piecewise(x,[x<mu,x>=mu],
									[lambda y : left.pdf(y)/np.max(left.pdf(y)),
									 lambda y : right.pdf(y)/np.max(right.pdf(y))])
				return(pdf/np.sum(pdf))
		
		def _argcheck(self,*args):
				return True

class warpModel(object):
	"""Base class for anything with parameters.

	Derived classes must have properties ``_param_names`` (list of str)
	and ``_parameters`` (1-d numpy.ndarray).
	"""

	def __init__(self, warp_function,parameters,param_names,warp_parameter,warp_distribution,name):
		self.name = name
		self._parameters = parameters
		self._param_names = [x.upper() for x in param_names]
		self.warp_function=warp_function
		self.warp_parameter=warp_parameter
		self.warp_distribution=warp_distribution


	def updateWarp_Param(self):
		self.warp_parameter=self.warp_distribution()[0]
		if self.name in self._param_names:
			self.set(**{self.name:self.warp_parameter})

	def flux(self,phase,wave,host_params,host_param_names):
		phase_wave_dict={'PHASE':phase,'WAVELENGTH':wave}
		self.set(**{p:host_params[host_param_names.index(p)] for p in self._param_names if p in host_param_names})
		parameter_arrays=[np.ones(len(wave))*self._parameters[i] if self._param_names[i] not in ['PHASE','WAVELENGTH'] else phase_wave_dict[self._param_names[i]] for i in range(len(self._param_names))]

		return(self.warp_function(np.vstack(parameter_arrays).T).flatten())




	@property
	def param_names(self):
		"""List of parameter names."""
		return self._param_names

	@property
	def parameters(self):
		"""Parameter value array"""
		return self._parameters

	@parameters.setter
	def parameters(self, value):
		value = np.asarray(value)
		if value.shape != self._parameters.shape:
			raise ValueError("Incorrect number of parameters.")
		self._parameters[:] = value

	def set(self, **param_dict):
		"""Set parameters of the model by name."""
		self.update(param_dict)

	def update(self, param_dict):
		"""Set parameters of the model from a dictionary."""
		for key, value in param_dict.items():
			self[key] = value

	def __setitem__(self, key, value):
		"""Set a single parameter of the model by name."""
		try:
			i = self._param_names.index(key)
		except ValueError:
			raise KeyError("Unknown parameter: " + repr(key))
		self._parameters[i] = value

	def get(self, name):
		"""Get parameter of the model by name."""
		return self[name]

	def __getitem__(self, name):
		"""Get parameter of the model by name"""
		try:
			i = self._param_names.index(name)
		except ValueError:
			raise KeyError("Model has no parameter " + repr(name))
		return self._parameters[i]


	def __str__(self):
		parameter_lines = [self._headsummary(), 'parameters:']
		if len(self._param_names) > 0:
			m = max(map(len, self._param_names))
			extralines = ['  ' + k.ljust(m) + ' = ' + repr(v)
						  for k, v in zip(self._param_names, self._parameters)]
			parameter_lines.extend(extralines)
		return '\n'.join(parameter_lines)

	def __copy__(self):
		"""Like a normal shallow copy, but makes an actual copy of the
		parameter array."""
		new_model = self.__new__(self.__class__)
		for key, val in self.__dict__.items():
			new_model.__dict__[key] = val
		new_model._parameters = self._parameters.copy()
		return new_model



def _skewed_normal(name,dist_dat):
		a=dist_dat['DIST_PEAK']-3*dist_dat['DIST_SIGMA'][0]
		b=dist_dat['DIST_PEAK']+3*dist_dat['DIST_SIGMA'][1]
		dist = skewed_normal(name,a=a,b=b)
		sample=np.arange(a,b,.01)
		return(lambda : np.random.choice(sample,1,
										 p=dist._pdf(sample,dist_dat['DIST_PEAK'],dist_dat['DIST_SIGMA'][0],dist_dat['DIST_SIGMA'][1])))
		

def _param_from_dist(dist_file,path):
	dist=np.loadtxt(os.path.join(path,dist_file))
	a=np.min(dist)-abs(np.min(dist))
	b=np.max(dist)+abs(np.max(dist))
	sample=np.arange(a,b,.01)
	pdf=gaussian_kde(dist.T).pdf(np.arange(a,b,.01))
	return(lambda : np.random.choice(sample,1,p=pdf/np.sum(pdf)))

def _get_distribution(name,dist_dat,path):
	if 'DIST_FILE' in dist_dat.keys():
		return(_param_from_dist(dist_dat['DIST_FILE'],path))

	return(_skewed_normal(name,dist_dat))

def _meshgrid2(*arrs):
	arrs = tuple(arrs)  #edit
	lens = list(map(len, arrs))
	dim = len(arrs)

	sz = 1
	for s in lens:
		sz*=s

	ans = []    
	for i, arr in enumerate(arrs):
		slc = [1]*dim
		slc[i] = lens[i]
		arr2 = np.asarray(arr).reshape(slc)
		for j, sz in enumerate(lens):
			if j!=i:
				arr2 = arr2.repeat(sz, axis=j) 
		ans.append(arr2)

	return tuple(ans)


def _generate_ND_grids(func,filename=None,colnames=None,*arrs):
	g=_meshgrid2(*arrs)
	positions = np.vstack(list(map(np.ravel, g))).T
	res=func(*(positions[:,i] for i in range(positions.shape[1]))).reshape((positions.shape[0],1))
	gridded=np.hstack([positions,res])
	if filename is not None:
		if colnames is not None:
			header=' '.join(colnames)
		else:
			header=''
		np.savetxt(filename,gridded,fmt='%f',header=header)
	return(gridded)
	
	
def _read_ND_grids(filename,scale_factor=1.):   
	with open(filename,'r') as f:
		temp=f.readline()
		
		if temp[0]=='#':
			names=temp.strip('#').split()
			gridded=pandas.read_csv(filename,sep=' ',names=names,comment='#',header=None)
		else:
			gridded=pandas.read_csv(filename,sep=' ',comment='#',header=None)
	
	arrs=tuple(np.unique(gridded.values[:,i]) for i in range(len(gridded.columns)-1))
	
	dim=[len(x) for x in arrs]

	theta=np.array(gridded[gridded.columns[-1]]).reshape(dim)*scale_factor
	
	
	return([x.upper() for x in gridded.columns][:-1],lambda interp_array:interpn(arrs,theta,xi=interp_array,method='linear',bounds_error=False,fill_value=0))

	
	
	
	

def main():
		
		#func=lambda x,y,z:np.array(x+y+z)
		#_generate_ND_grids(func,'test_host.dat',['phase','mass','velocity','theta'],np.array([0,20]),np.array([2,3]),np.array([4,5]))
		#sys.exit()
		#p,test=_read_ND_grids('salt2_m0.dat')
		#print(test(np.array([[10,5000],[10,6000]])))
		import matplotlib.pyplot as plt
		#sys.exit()
		mySED=genmag_BYOSED('$WFIRST_ROOT/BYOSED_dev/BYOSEDINPUT/',2,[],['HOST_MASS','SFR','AGE','REDSHIFT'])
		#plt.plot(mySED.wave,mySED.sedInterp(0,mySED.wave)/mySED.x0)
		#f=mySED.sedInterp(0,mySED.wave).flatten()/mySED.x0
		#s=mySED.sn_effects['STRETCH'].flux(0*np.ones(len(mySED.wave)),mySED.wave,[],[])

		#plt.plot(mySED.wave,f+s)
		#plt.xlim(3400,10000)
		#plt.show()
		#sys.exit()
		mySED.sn_effects['VELOCITY'].set(VELOCITY=0)
		mySED.sn_effects['STRETCH'].warp_parameter=0
		plt.plot(mySED.wave,mySED.fetchSED_BYOSED(0,5000,3,3,[2.5,1,1,.5]),label=0)

		for p in range(3):
			mySED.sn_effects['VELOCITY'].updateWarp_Param()
			v=mySED.sn_effects['VELOCITY'].warp_parameter
			mySED.sn_effects['STRETCH'].updateWarp_Param()
			s=mySED.sn_effects['STRETCH'].warp_parameter
			plt.plot(mySED.wave,mySED.fetchSED_BYOSED(0,5000,3,3,[2.5,1,1,.5]),label=(v,s))
		plt.xlim(3400,10000)
		plt.legend()
		plt.show()



if __name__=='__main__':
		main()
