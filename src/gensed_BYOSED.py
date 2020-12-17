import numpy as np  
import os,six,abc,math
import optparse
import configparser
import pandas
import sys
from scipy.interpolate import RectBivariateSpline,interp1d,interp2d,interpn
from ast import literal_eval
from scipy.stats import rv_continuous,gaussian_kde,norm as normal
from copy import copy,deepcopy
import pickle

if not hasattr(sys, 'argv'):
        sys.argv  = ['']

required_keys = []

__mask_bit_locations__={'verbose':1,'dump':2}

def print_err():
    print("""
               ______
             /    x  \\
            /   --------<  ABORT Python on Fatal Error.
        __ /  _______/
/^^^^^^^^^^^^^^/  __/
\________________/
                """)
    raise RuntimeError


class gensed_BYOSED:
        def __init__(self,PATH_VERSION,OPTMASK,ARGLIST,HOST_PARAM_NAMES):
            # TODO: write a print statement that warns if
            # HOST_PARAM_NAMES is a variable that the code
            # isn't going to do anything with
            try:
                self.verbose = OPTMASK & (1 << __mask_bit_locations__['verbose']) > 0

                self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))

                self.host_param_names = [x.upper() for x in HOST_PARAM_NAMES.split(',')]
                self.PATH_VERSION = os.path.dirname(PATH_VERSION)

                self.dump = OPTMASK & (1 << __mask_bit_locations__['dump'])>0
                self.sn_id=None

                self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))
                
                if os.path.exists(os.path.join(self.PATH_VERSION,'BYOSED.params')):
                    self.paramfile = os.path.join(self.PATH_VERSION,'BYOSED.params')
                elif os.path.exists(os.path.join(self.PATH_VERSION,'byosed.params')):
                    self.paramfile = os.path.join(self.PATH_VERSION,'byosed.params')
                else:
                    raise RuntimeError('param file %s not found!'%os.path.join(self.PATH_VERSION,'BYOSED.params'))

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

                self.warp_effects=self.fetchParNames_CONFIG(config)

                self.sn_effects,self.host_effects=self.fetchWarp_BYOSED(config)

                phase,wave,flux = np.loadtxt(os.path.join(self.PATH_VERSION,self.options.sed_file),unpack=True)


                fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
                self.norm =float(config['MAIN']['NORM']) if 'NORM' in config['MAIN'].keys() else -19.365
                self.alpha=float(config['MAIN']['ALPHA']) if 'ALPHA' in config['MAIN'].keys() else 0.14
                self.beta=float(config['MAIN']['BETA']) if 'BETA' in config['MAIN'].keys() else 3.1

                typ=config['MAIN']['SNTYPE'] if 'SNTYPE' in config['MAIN'].keys() else 'IA'
                if typ=='IA':
                    self.is_Ia=True
                else:
                    self.is_Ia=False

                self.x0=10**(-.4*self.norm)
                
                if self.options.magsmear!=0.0:
                    self.magsmear = np.random.normal(0,self.options.magsmear)
                else:
                    self.magsmear = 0.0
                self.magoff=self.options.magoff
                
                self.flux = fluxarr*self.x0*10**(-0.4*self.magoff)

                self.phase = np.unique(phase)
                self.wave = np.unique(wave)
                self.wavelen = len(self.wave)

                self.sedInterp=interp2d(self.phase,self.wave,self.flux.T,kind='linear',bounds_error=True)
                
                self.phase_data={}
                
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print('Python Error :',e)
                print('gensed_BYOSED.py, line number: %i'%exc_tb.tb_lineno)
                print_err()

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

            sn_dict=dict([])
            host_dict=dict([])


            for warp in self.warp_effects:
                warp_data={}
                for k in config[warp]:
                    try:
                        warp_data[k.upper()]=np.array(config.get(warp,k).split()).astype(float)
                    except:
                        warp_data[k.upper()]=config.get(warp,k)
                if 'SCALE_TYPE' not in warp_data.keys():
                    warp_data['SCALE_TYPE']='inner'
                elif warp_data['SCALE_TYPE'] not in ['inner','outer']:
                    raise RuntimeError("Do not recognize variable SCALE_TYPE, should be 'inner' or 'outer'")
                


                if 'SN_FUNCTION' in warp_data:
                    if 'DIST' in ' '.join([x for x in warp_data.keys() if 'SN' in x or 'SCALE' in x]):
                        distribution=_get_distribution(warp,{k:warp_data[k] for k in warp_data.keys() if 'SN' in k or 'SCALE' in k},self.PATH_VERSION,'SN')
                    else:
                        raise RuntimeError("Did not supply scale distribution information for SN effect %s."%warp)


                    
                    if 'SN_FUNCTION_SCALE' not in warp_data:
                        scale_factor=1.
                    else:
                        scale_factor=warp_data['SN_FUNCTION_SCALE']
                    try:
                        sn_param_names,sn_function=_read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['SN_FUNCTION']))),scale_factor)
                    except RuntimeError:
                        raise RuntimeError("Do not recognize format of function for %s SN Function"%warp)
                    if warp.upper() in sn_param_names and 'PARAM' not in distribution.keys():
                        raise RuntimeError("Must supply parameter distribution for SN effect %s"%warp)
                    sn_scale_parameter=distribution['SCALE'](-1)
                    warp_parameter=distribution['PARAM'](-1) if 'PARAM' in distribution.keys() else None
                    if warp.upper() in sn_param_names and warp_parameter is None:
                        raise RuntimeError("Woops, you are not providing a PARAM distribution for your %s effect."%warp.upper())

                    sn_dict[warp]=WarpModel(warp_function=sn_function,
                                            param_names=sn_param_names,
                                            parameters=np.array([0. if sn_param_names[i]!=warp.upper() else warp_parameter for i in range(len(sn_param_names))]),
                                            warp_parameter=warp_parameter,
                                            warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
                                            scale_parameter=sn_scale_parameter,
                                            scale_distribution=distribution['SCALE'],
                                            scale_type=warp_data['SCALE_TYPE'],
                                            name=warp)

                if 'HOST_FUNCTION' in warp_data:
                    if 'DIST' in ' '.join([x for x in warp_data.keys() if 'HOST' in x or 'SCALE' in x]):
                        distribution=_get_distribution(warp,{k:warp_data[k] for k in warp_data.keys() if 'HOST' in k or 'SCALE' in k},self.PATH_VERSION,'HOST')
                    else:
                        raise RuntimeError("Did not supply scale distribution information for HOST effect %s."%warp)
                    try:
                        host_param_names,host_function=_read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['HOST_FUNCTION']))))
                    except RuntimeError:
                        raise RuntimeError("Do not recognize format of function for %s HOST Function"%warp)

                    if warp.upper() in host_param_names and 'PARAM' not in distribution.keys() and warp.upper() not in self.host_param_names:
                        raise RuntimeError("Must supply parameter distribution for HOST effect %s"%warp)

                    host_scale_parameter=distribution['SCALE'](-1)
                    warp_parameter=distribution['PARAM'](-1) if 'PARAM' in distribution.keys() else None
                    if warp.upper() in host_param_names and warp_parameter is None and warp.upper() not in self.host_param_names:
                        raise RuntimeError("Woops, you are not providing a PARAM distribution for your %s effect."%warp.upper())
                    host_dict[warp]=WarpModel(warp_function=host_function,
                                              param_names=host_param_names,
                                              parameters=np.array([0. if host_param_names[i]!=warp.upper() else warp_parameter for i in range(len(host_param_names))]),

                                              warp_parameter=warp_parameter,
                                              warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
                                              scale_parameter=host_scale_parameter,
                                              scale_distribution=distribution['SCALE'],
                                              scale_type=warp_data['SCALE_TYPE'],
                                              name=warp)

            return(sn_dict,host_dict)
        
                
        def fetchSED_NLAM(self):
                return self.wavelen

        def fetchSED_LAM(self):
                return list(self.wave)
        

        def fetchSED_BYOSED(self,trest,maxlam=5000,external_id=1,new_event=1,hostpars='',mB_calc=False):
                    
            try:
                if len(self.wave)>maxlam:
                    raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))

                if self.sn_id is None or external_id==self.sn_id:
                    newSN=False
                else:
                    newSN=True
                    self.phase_data={}

                self.sn_id=external_id
                if not newSN and np.round(trest,6) in self.phase_data.keys():
                    return deepcopy(self.phase_data[np.round(trest,6)])
                    
                fluxsmear=self.sedInterp(trest,self.wave).flatten()
                orig_fluxsmear=copy(fluxsmear)

                if self.options.magsmear!=0.0 and (self.sn_id!=external_id or self.magsmear is None):
                    self.magsmear=np.random.normal(0,self.options.magsmear)
                else:
                    self.magsmear=0.0

                fluxsmear *= 10**(-0.4*(self.magsmear))
                trest_arr=trest*np.ones(len(self.wave))
                inner_product=np.zeros(len(self.wave))
            except Exception as e:
                print('Python Error :',e)
                print_err()
            
            
            outer_product=np.zeros(len(self.wave))

            for warp in [x for x in self.warp_effects]:
                try: 
                    
                    if newSN:
                        z=hostpars[self.host_param_names.index('REDSHIFT')] if\
                           'REDSHIFT' in self.host_param_names else None
                        if warp in self.sn_effects.keys():
                            self.sn_effects[warp].updateWarp_Param(z)
                            self.sn_effects[warp].updateScale_Param(z)
                            if warp in self.host_effects.keys():
                                self.host_effects[warp].updateWarp_Param(z)
                                self.host_effects[warp].scale_parameter=1.

                        else:
                            self.host_effects[warp].updateWarp_Param(z)
                            self.host_effects[warp].updateScale_Param(z)



                    product=np.ones(len(self.wave))
                    temp_scale_param = 0
                    temp_outer_product=np.ones(len(self.wave))
                    outer_scale_param=0


                    if warp in self.sn_effects.keys():
                        if self.verbose:
                            if self.sn_effects[warp].warp_parameter is not None:
                                print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].warp_parameter))
                            else:
                                print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].scale_parameter))

                        if self.sn_effects[warp].scale_type=='inner':

                            product*=self.sn_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
                            if temp_scale_param ==0:
                                temp_scale_param = self.sn_effects[warp].scale_parameter
                            else:
                                temp_scale_param*=self.sn_effects[warp].scale_parameter
                            
                        else:

                            temp_outer_product*=self.sn_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
                            if outer_scale_param ==0:
                                outer_scale_param = self.sn_effects[warp].scale_parameter
                            else:
                                outer_scale_param*=self.sn_effects[warp].scale_parameter
                            

                    if warp in self.host_effects.keys():
                        
                        if self.verbose:
                            if self.host_effects[warp].warp_parameter is not None:
                                print('Phase=%.1f, %s: %.2f'%(trest,warp,self.host_effects[warp].warp_parameter))
                            else:
                                print('Phase=%.1f, %s: %.2f'%(trest,warp,self.host_effects[warp].scale_parameter))
                        if self.host_effects[warp].scale_type=='inner':
                                if temp_scale_param==0:
                                        temp_scale_param=self.host_effects[warp].scale_parameter
                                else:
                                        temp_scale_param*=self.host_effects[warp].scale_parameter
                                product*=self.host_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
                                
                        else:
                                
                                if outer_scale_param==0:
                                        outer_scale_param=self.host_effects[warp].scale_parameter
                                else:
                                        outer_scale_param*=self.host_effects[warp].scale_parameter
                                temp_outer_product*=self.host_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
                        
                    inner_product+=product*temp_scale_param
                    outer_product+=temp_outer_product*outer_scale_param

                except Exception as e:
                    print('Python Error :',e)
                    print_err()
            
        

            fluxsmear*=((1+inner_product)*10**(-0.4*outer_product))

            if self.is_Ia:
                fluxsmear*=self.brightness_correct_Ia()

            self.phase_data[np.round(trest,6)]=deepcopy(list(fluxsmear))
            return list(fluxsmear)
            
        def brightness_correct_Ia(self):
            if 'COLOR' in self.sn_effects.keys():
                c=self.sn_effects['COLOR'].scale_parameter
            else:
                c=0
            if 'STRETCH' in self.sn_effects.keys():
                s=self.sn_effects['STRETCH'].scale_parameter
            else:
                s=0
            return(10**(-.4*(self.beta*c-self.alpha*s)))

        
        def fetchParNames_BYOSED(self):
                return list(np.append(self.warp_effects,['lum']))

        def fetchNParNames_BYOSED(self):
                return len(self.warp_effects)+1

        def fetchParVals_BYOSED_4SNANA(self,varname):
                if varname=='lum':
                    return self.x0*(10**(-0.4*self.magsmear))
                if varname in self.sn_effects.keys():
                    if self.sn_effects[varname].warp_parameter is not None:
                        return self.sn_effects[varname].warp_parameter
                    else:

                        return self.sn_effects[varname].scale_parameter
                else:
                    if self.host_effects[varname].warp_parameter is not None:
                        return self.host_effects[varname].warp_parameter
                    else:
                        return self.host_effects[varname].scale_parameter
                    

        def fetchParVals_BYOSED(self,varname):
                if varname in self.sn_effects.keys():
                    return self.sn_effects[varname].warp_parameter
                else:
                    return self.host_effects[varname].warp_parameter
        def fetchParNames_CONFIG(self,config):
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

class WarpModel(object):
    """Base class for anything with parameters.

    Derived classes must have properties ``_param_names`` (list of str)
    and ``_parameters`` (1-d numpy.ndarray).
    """

    def __init__(self, warp_function,parameters,param_names,warp_parameter,warp_distribution,
                 scale_parameter,scale_distribution,scale_type,name):
        self.name = name
        self._parameters = parameters
        self._param_names = [x.upper() for x in param_names]
        self.warp_function=warp_function
        self.warp_parameter=warp_parameter
        self.scale_parameter=scale_parameter
        self.warp_distribution=warp_distribution
        self.scale_distribution=scale_distribution
        self.scale_type=scale_type

    def updateWarp_Param(self,z=None):
        if self.warp_distribution is not None:
            self.warp_parameter=self.warp_distribution(z)

            if self.name in self._param_names:
                self.set(**{self.name:self.warp_parameter})
        

    def updateScale_Param(self,z=None):
        self.scale_parameter=self.scale_distribution(z)

        
    def flux(self,phase,wave,host_params,host_param_names):
        phase_wave_dict={'PHASE':phase,'WAVELENGTH':wave}
        self.set(**{p:host_params[host_param_names.index(p)] for p in self._param_names if p in host_param_names})
        
        parameter_arrays=[np.ones(len(wave))*self._parameters[i] if self._param_names[i] not in ['PHASE','WAVELENGTH'] 
                            else phase_wave_dict[self._param_names[i]] for i in range(len(self._param_names))]

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
        temp_ps=[]
        for i in range(len(self._param_names)):
            if self._param_names[i] in param_dict.keys():
                temp_ps.append(param_dict[self._param_names[i]])
            else:
                temp_ps.append(self._parameters[i])
        self._parameters[:]=temp_ps


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



def _skewed_normal(name,dist_dat,dist_type):
        if dist_type+'_DIST_LIMITS' in dist_dat:
            a,b=dist_dat[dist_type+'_DIST_LIMITS']
        else:
            a=dist_dat[dist_type+'_DIST_PEAK']-3*dist_dat[dist_type+'_DIST_SIGMA'][0]
            b=dist_dat[dist_type+'_DIST_PEAK']+3*dist_dat[dist_type+'_DIST_SIGMA'][1]
        if a==b:
                return(lambda :[a])
        dist = skewed_normal(name,a=a,b=b)
        sample=np.linspace(a,b,int(1e4))
        return(lambda : np.random.choice(sample,1,
                                         p=dist._pdf(sample,dist_dat[dist_type+'_DIST_PEAK'],dist_dat[dist_type+'_DIST_SIGMA'][0],dist_dat[dist_type+'_DIST_SIGMA'][1])))
        

def _param_from_dist(dist_file,path):
    dist=np.loadtxt(os.path.join(path,dist_file))
    a=np.min(dist)-abs(np.min(dist))
    b=np.max(dist)+abs(np.max(dist))
    sample=np.linspace(a,b,int(1e4))
    pdf=gaussian_kde(dist.T).pdf(np.linspace(a,b,int(1e4)))
    return(lambda : np.random.choice(sample,1,p=pdf/np.sum(pdf)))

def _get_zdepend(dist_file,path,typ):
    fv=0 if typ=='add' else 1
    z,factor=np.loadtxt(os.path.join(path,dist_file),unpack=True)
    return(interp1d(z,factor,fill_value=fv,bounds_error=False))


def _get_distribution(name,dist_dat,path,sn_or_host):
    dist_dict={}
    param=True
    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'PARAM' in x]):
        
        if sn_or_host+'_PARAM_DIST_FILE' in dist_dat.keys():
            param_func=_param_from_dist(dist_dat[sn_or_host+'_PARAM_DIST_FILE'],path)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'PARAM_DIST_FILE' for %s?"%sn_or_host)
        
        
    elif np.any(['PARAM' in x for x in dist_dat.keys()]):
        try:
            param_func=_skewed_normal(name,dist_dat,sn_or_host+'_PARAM')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s param distribution(s)."%sn_or_host)
    else:
        param=False

    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'SCALE' in x]):
        if 'SCALE_DIST_FILE' in dist_dat.keys():
            scale_func=_param_from_dist(dist_dat['SCALE_DIST_FILE'],path)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'SCALE_DIST_FILE'?")
        
    elif np.any(['SCALE' in x for x in dist_dat.keys()]):
        try:
            scale_func=_skewed_normal(name,dist_dat,'SCALE')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s scale distribution(s)."%sn_or_host)
    else:
        raise RuntimeError("Must supply scale distribution for every effect.")

    if np.any(['ZDEPEND' in x for x in dist_dat.keys() if 'PARAM' in x]):
        if sn_or_host+'_PARAM_ZDEPEND_FILE' in dist_dat.keys():
            
            if sn_or_host+'_PARAM_ZDEPEND_TYPE' in dist_dat.keys() and dist_dat[sn_or_host+'_PARAM_ZDEPEND_TYPE']=='MULTIPLY':
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_PARAM_ZDEPEND_FILE'],path,'multiply')
                dist_dict['PARAM']=lambda z: (param_func()[0])*zdepend(z)
            else:
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_PARAM_ZDEPEND_FILE'],path,'add')
                dist_dict['PARAM']=lambda z: param_func()[0]+zdepend(z)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'ZDEPEND' for %s?"%sn_or_host)
    elif param:
        dist_dict['PARAM']=lambda z:param_func()[0]
    if np.any(['ZDEPEND' in x for x in dist_dat.keys() if 'SCALE' in x]):
        if sn_or_host+'_SCALE_ZDEPEND_FILE' in dist_dat.keys():
            if sn_or_host+'_SCALE_ZDEPEND_TYPE' in dist_dat.keys() and\
               dist_dat[sn_or_host+'_SCALE_ZDEPEND_TYPE'].upper()=='MULTIPLY':
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_SCALE_ZDEPEND_FILE'],path,'multiply')
                dist_dict['SCALE']=lambda z: (scale_func()[0])*zdepend(z)
            else:
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_SCALE_ZDEPEND_FILE'],path,'add')
                dist_dict['SCALE']=lambda z: scale_func()[0]+zdepend(z)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'ZDEPEND' for %s scale distribution?"%sn_or_host)
    else:
        dist_dict['SCALE']=lambda z:scale_func()[0]
    return(dist_dict)

def _integration_grid(low, high, target_spacing):
    """Divide the range between `start` and `stop` into uniform bins
    with spacing less than or equal to `target_spacing` and return the
    bin midpoints and the actual spacing."""

    range_diff = high - low
    spacing = range_diff / int(math.ceil(range_diff / target_spacing))
    grid = np.arange(low + 0.5 * spacing, high, spacing)

    return grid, spacing

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
                    
        import matplotlib.pyplot as plt

        mySED=gensed_BYOSED('$WFIRST_USERS/jpierel/hyperion_wfirst/byosed/',2,[],'REDSHIFT,AGE,ZCMB,METALLICITY,HOSTMASS')
        mySED.sn_id=1
        print(np.sum(mySED.fetchSED_BYOSED(0,5000,0,0,[.1,1,1,.5,11])))
        print(np.sum(mySED.fetchSED_BYOSED(0,5000,0,1,[.1,1,1,.5,9])))

if __name__=='__main__':
        main()
