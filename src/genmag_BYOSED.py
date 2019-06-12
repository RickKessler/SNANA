import numpy as np	
import os,six
import optparse
import configparser
import sys
from scipy.interpolate import RectBivariateSpline,interp2d
from ast import literal_eval
from scipy.stats import rv_continuous,norm as normal
from copy import copy
import pickle

if not hasattr(sys, 'argv'):
		sys.argv  = ['']

required_keys = ['magsmear']

__mask_bit_locations__={'verbose':1,'dump':2}
__host_param_indices=['mass','sfr','age']

class genmag_BYOSED:
		def __init__(self,PATH_VERSION,OPTMASK,ARGLIST,NAMES_HOSTPAR):
			try:
				self.verbose = OPTMASK & (1 << __mask_bit_locations__['verbose']) > 0

				self.dump = OPTMASK & (1 << __mask_bit_locations__['dump'])>0
				self.sn_id=None

<<<<<<< HEAD
				self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))

=======
				self.NAMES_HOSTPAR = NAMES_HOSTPAR.split(',')
				self.PATH_VERSION = os.path.dirname(PATH_VERSION)
>>>>>>> 01492262c07b0234a23e0b7c77935d68274a4a89
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
				
				self.warping_functions,self.warping_distributions,self.warping_params=self.fetchWarp_BYOSED(config)


				phase,wave,flux = np.loadtxt(os.path.join(self.PATH_VERSION,self.options.sed_file),unpack=True)



				fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
				self.x0=10**(.4*19.365)
				self.flux = fluxarr*self.x0

				self.phase = np.unique(phase)
				self.wave = np.unique(wave)
				self.wavelen = len(self.wave)

				self.sedInterp=interp2d(self.phase,self.wave,self.flux.T,kind='linear',bounds_error=True)
			except:
				import pdb; pdb.set_trace()
	
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

				func_dict=dict([])
				dist_dict=dict([])
				param_dict=dict([])
				for warp in self.warp_effects:
						warp_data={k.upper():np.array(config.get(warp,k).split()).astype(float) if k.upper()!='WARP_FUNCTION' else config.get(warp,k) for k in config[warp]}
						dist_dict[warp]=_skewed_normal(warp,warp_data)
						
						
						func_data=np.array(warp_data['WARP_FUNCTION'])
								
						try:
								x0,x1,y=_read_griddata(os.path.expandvars(func_data))
						except:
								try:
										x0,x1,y=_sncosmo_read_griddata(os.path.expandvars(os.path.join(self.PATH_VERSION,str(func_data))))
								except RuntimeError:
										raise RuntimeError("Do not recognize format of function for %s"%warp)
								
								
						
						func_dict[warp]=RectBivariateSpline(x0,x1,y,kx=3,ky=3)
						param_dict[warp]=dist_dict[warp]()[0]
						
								
				return(func_dict,dist_dict,param_dict)
		

		def updateWarping_Params(self):
				for warp in self.warp_effects:
						self.warping_params[warp]=self.warping_distributions[warp]()[0]
				return self
						
		def applyWarping_Effects(self):
				tempSED=None
				for warp in self.warp_effects:
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
		
<<<<<<< HEAD
		def fetchSED_BYOSED(self,trest,maxlam,external_id,new_event,HOST_PARAMS):
=======
		def fetchSED_BYOSED(self,trest,maxlam,external_id,new_event,hostpars):
>>>>>>> 01492262c07b0234a23e0b7c77935d68274a4a89

				if len(self.wave)>maxlam:
						raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))
				#iPhase = np.where(np.abs(trest-self.phase) == np.min(np.abs(trest-self.phase)))[0][0]
				if self.sn_id is None:
						self.sn_id=external_id
				fluxsmear=self.sedInterp(trest,self.wave).flatten()
				if self.options.magsmear!=0.0:
						fluxsmear *= 10**(0.4*(np.random.normal(0,self.options.magsmear)))

				for warp in [x for x in self.warp_effects if x!='COLOR']:
						if True:
							#print(external_id,self.sn_id)
							if external_id!=self.sn_id:
								self.updateWarping_Params()
								self.sn_id=external_id
								
							if self.verbose:
								print('Phase=%.1f, %s: %.2f'%(trest,warp,self.warping_params[warp]))
							# not sure about the multiplication by x0 here, depends on if SNANA is messing with the 
							# absolute magnitude somewhere else
							fluxsmear+=self.warping_params[warp]*\
										self.warping_functions[warp](trest,self.wave).flatten()*self.x0
						#except:
							#import pdb; pdb.set_trace()

				if 'COLOR' in self.warp_effects:
						if external_id!=self.sn_id:
								self.updateWarping_Params()
								self.sn_id=external_id
						if self.verbose:
								print('Phase=%.1f, %s: %.2f'%(trest,'COLOR',self.warping_params['COLOR']))
						fluxsmear*=10**(-0.4*self.warping_params['COLOR']*\
								self.warping_functions['COLOR'](trest,self.wave).flatten())
							  
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
						return([k.upper() for k in list(config['FLAGS'].keys()) if config['FLAGS'][k]])
				else:
						return([x for x in config.sections() if x not in ['MAIN','FLAGS']])

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
				pdf=np.piecewise(x,[x<mu,x>=mu],
									[lambda y : left.pdf(y)/np.max(left.pdf(y)),
									 lambda y : right.pdf(y)/np.max(right.pdf(y))])
				return(pdf/np.sum(pdf))
		
		def _argcheck(self,*args):
				return True


def _skewed_normal(name,dist_dat):
		dist = skewed_normal(name,a=min(dist_dat['DIST_SIGMA']),b=max(dist_dat['DIST_SIGMA']))
		sample=np.arange(min(dist_dat['DIST_SIGMA']),max(dist_dat['DIST_SIGMA']),.01)
		return(lambda : np.random.choice(sample,1,
										 p=dist._pdf(sample,dist_dat['DIST_PEAK'],dist_dat['DIST_SIGMA'][0],dist_dat['DIST_SIGMA'][1])))
		
def _sncosmo_read_griddata(name_or_obj):
		if isinstance(name_or_obj, six.string_types):
				f = open(name_or_obj, 'r')
		else:
				f = name_or_obj
		x0 = []	   # x0 values.
		x1 = None  # x1 values for first x0 value, assume others are the same.
		y = []	   # 2-d array of internal values

		x0_current = None
		x1_current = []
		y1_current = []
		for line in f:
				stripped_line = _stripcomment(line)
				if len(stripped_line) == 0:
						continue
				x0_tmp, x1_tmp, y_tmp = map(float, stripped_line.split())
				if x0_current is None:
						x0_current = x0_tmp	 # Initialize first time

				# If there is a new x0 value, ingest the old one and reset values
				if x0_tmp != x0_current:
						x0.append(x0_current)
						if x1 is None:
								x1 = x1_current
						y.append(y1_current)

						x0_current = x0_tmp
						x1_current = []
						y1_current = []
				x1_current.append(x1_tmp)
				y1_current.append(y_tmp)

		# Ingest the last x0 value and y1 array
		x0.append(x0_current)
		y.append(y1_current)

		f.close()
		return np.array(x0), np.array(x1), np.array(y)

def _stripcomment(line, char='#'):
		pos = line.find(char)
		if pos == -1:
				return line
		else:
				return line[:pos]

def main():
		mySED=genmag_BYOSED('$WFIRST_ROOT/BYOSED_dev/BYOSEDINPUT/',2,[])
		temp=mySED.fetchSED_BYOSED(10,5000,3,0,[1,1,1])


if __name__=='__main__':
		main()
