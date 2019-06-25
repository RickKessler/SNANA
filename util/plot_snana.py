#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os,glob,math,sys
from optparse import OptionParser
from scipy.interpolate import interp1d


__band_order__=np.append(['u','b','g','r','i','z','y','j','h','k'],
	 [x.upper() for x in ['u','b','g','r','i','z','y','j','h','k']])

def read_spec(cid,base_name):
	names=['wave','flux','fluxerr','tobs']
	id_to_obs=dict([])
	with open(base_name+".SPECLIST.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		else:
			id_to_obs[int(temp[varnames.index('CID')])]=float(temp[varnames.index('TOBS')])
	sn={k:[] for k in names}

	with open(base_name+".SPECPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		
		
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		elif len(temp)>0 and b'OBS:' in temp and\
			 str(temp[varnames.index('CID')].decode('utf-8'))in cid:
			sn['wave'].append((float(temp[varnames.index('LAMMAX')])+float(temp[varnames.index('LAMMIN')]))/2.)
			sn['flux'].append(float(temp[varnames.index('FLAM')]))
			sn['fluxerr'].append(float(temp[varnames.index('FLAMERR')]))
			sn['tobs'].append(id_to_obs[int(temp[varnames.index('CID')])])
	sn={k:np.array(sn[k]) for k in sn.keys()}
	return(sn)
def read_lc(cid,base_name):
	names=['time','flux','fluxerr','filter','chi2']
	peak=None
	sn={k:[] for k in names} 
	fit={k:[] for k in ['time','flux','filter']}
	with open(base_name+".LCPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		elif len(temp)>0 and b'OBS:' in temp and str(temp[varnames.index('CID')].decode('utf-8')) in cid:
			if int(temp[varnames.index('DATAFLAG')])==1:
				if peak is None:
					peak=float(temp[varnames.index('MJD')])-float(temp[varnames.index('Tobs')])
				sn['time'].append(float(temp[varnames.index('Tobs')]))
				sn['flux'].append(float(temp[varnames.index('FLUXCAL')]))
				sn['fluxerr'].append(float(temp[varnames.index('FLUXCAL_ERR')]))
				sn['filter'].append(str(temp[varnames.index('BAND')].decode('utf-8')))
				sn['chi2'].append(float(temp[varnames.index('CHI2')]))
			elif int(temp[varnames.index('DATAFLAG')])==0:
				fit['time'].append(float(temp[varnames.index('Tobs')]))
				fit['flux'].append(float(temp[varnames.index('FLUXCAL')]))
				fit['filter'].append(str(temp[varnames.index('BAND')].decode('utf-8')))
	
	sn={k:np.array(sn[k]) for k in sn.keys()}
	fit={k:np.array(fit[k]) for k in fit.keys()}
	if len(fit['filter'])>0:
		fits={k:interp1d(fit['time'][fit['filter']==k],
					 fit['flux'][fit['filter']==k]) for k in np.unique(fit['filter'])}
	else:
		fits=[]
	return(sn,fits,peak)

def plot_spec(cid,bin_size,base_name):
	sn=read_spec(cid,base_name)
	
	if len(sn['tobs'])==0:
		return
	if len(np.unique(sn['tobs']))>1:
		fig,ax=plt.subplots(nrows=len(np.unique(sn['tobs'])),ncols=1,figsize=(8,8),sharex=True)
		ax[0].set_title('SN%s'%cid[0],fontsize=16)
		for j in range(len(np.unique(sn['tobs']))):

			temp_sn=np.where(sn['tobs']==np.unique(sn['tobs'])[j])[0]
			if bin_size!=0:
				
				#bins=np.digitize(np.array(sn['wave']),np.arange(sn['wave'][0],sn['wave'][-1],bin_size))
				binned_wave=[]
				binned_flux=[]
				binned_fluxerr=[]
				bins=np.trunc(sn['wave'][temp_sn]/bin_size)
				for i in np.unique(bins):
					binned_wave=np.append(binned_wave,np.mean(sn['wave'][temp_sn][bins==i]))
					binned_flux=np.append(binned_flux,np.mean(sn['flux'][temp_sn][bins==i]))
					binned_fluxerr=np.append(binned_fluxerr,np.mean(sn['fluxerr'][temp_sn][bins==i]))
			else:
				binned_wave=sn['wave'][temp_sn]
				binned_flux=sn['flux'][temp_sn]
				binned_fluxerr=sn['fluxerr'][temp_sn]
				#sn=(sn.group_by(np.trunc(sn['wave']/bin_size))).groups.aggregate(np.mean)
				
			ax[j].plot(binned_wave,binned_flux,color='k',label='TOBS:%.2f'%np.unique(sn['tobs'])[j])
			ylim=ax[j].get_ylim()
			ax[j].fill_between(binned_wave,binned_flux-binned_fluxerr,binned_flux+binned_fluxerr,
							 color='r',alpha=.3,label=r'$1\sigma$ Error')
			ax[j].plot([binned_wave[0],binned_wave[-1]],[0,0],'k--',alpha=.5)
			ax[j].set_ylim(ylim)
			ax[j].legend(fontsize=16)
			
			ax[j].set_ylabel('Flux',fontsize=16)
		ax[j].set_xlabel('Observer Frame Wavelength ($\AA$)',fontsize=16)

	else:
		fig=plt.figure(figsize=(10,8))
		if bin_size!=0:
			
			#bins=np.digitize(np.array(sn['wave']),np.arange(sn['wave'][0],sn['wave'][-1],bin_size))
			binned_wave=[]
			binned_flux=[]
			binned_fluxerr=[]
			bins=np.trunc(sn['wave']/bin_size)
			for i in np.unique(bins):
				binned_wave=np.append(binned_wave,np.mean(sn['wave'][bins==i]))
				binned_flux=np.append(binned_flux,np.mean(sn['flux'][bins==i]))
				binned_fluxerr=np.append(binned_fluxerr,np.mean(sn['fluxerr'][bins==i]))
		else:
			binned_wave=sn['wave']
			binned_flux=sn['flux']
			binned_fluxerr=sn['fluxerr']
			#sn=(sn.group_by(np.trunc(sn['wave']/bin_size))).groups.aggregate(np.mean)
			
		plt.plot(binned_wave,binned_flux,color='k',label='TOBS:%.2f'%np.unique(sn['tobs'])[0])
		ylim=plt.ylim()
		plt.fill_between(binned_wave,binned_flux-binned_fluxerr,binned_flux+binned_fluxerr,
						 color='r',alpha=.3,label=r'$1\sigma$ Error')
		plt.plot([binned_wave[0],binned_wave[-1]],[0,0],'k--',alpha=.5)
		plt.ylim(ylim)
		plt.legend(fontsize=16)
		plt.xlabel('Observer Frame Wavelength ($\AA$)',fontsize=16)
		plt.ylabel('Flux',fontsize=16)
	
	plt.savefig('SNANA_SPEC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)

def plot_lc(cid,base_name):
	sn,fits,peak=read_lc(cid,base_name)
	if len(sn['time'])==0:
		return
	rows=int(math.ceil(len(np.unique(sn['filter']))))

	fig,ax=plt.subplots(nrows=rows,ncols=1,figsize=(8,8),sharex=True)
	ax[0].set_title('SN%s'%cid[0],fontsize=16)
	i=0
	for band in np.append([x for x in __band_order__ if x in np.unique(sn['filter'])],
						[x for x in np.unique(sn['filter']) if x not in __band_order__]):

		temp_sn={k:sn[k][np.where(sn['filter']==band)[0]] for k in sn.keys()}
		chi2=np.mean(temp_sn['chi2'])
		if chi2>0:
			lab=r'%s: $\chi^2$=%.1f'%(band,np.mean(temp_sn['chi2']))
			leg_size=12
		else:
			lab=band
			leg_size=16

		ax[i].errorbar(temp_sn['time'],temp_sn['flux'],yerr=temp_sn['fluxerr'],
						  fmt='.',markersize=8,color='k',
						  label=lab)
		if len(fits)>0:
			fit_time=np.arange(temp_sn['time'][0],temp_sn['time'][-1],1)
			ax[i].plot(fit_time,fits[band](fit_time),color='r',label='Best Fit',linewidth=3)
		ax[i].legend(fontsize=leg_size)
		ax[i].set_ylabel('Flux',fontsize=16)
		i+=1
	ax[i-1].set_xlabel('Time-%.2f (Rest Frame Days)'%peak,fontsize=16)
	#fig.text(0.5, 0.02, 'Time (Rest Frame Days)', ha='center',fontsize=16)
	#fig.text(0.04, .5, 'Flux', va='center', rotation='vertical',fontsize=16)
		
	plt.savefig('SNANA_LC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)

def plot_cmd(genversion,cid_list):
	if os.path.splitext(genversion)[1]=='.NML':
		plotter='salt2'
	else:
		plotter='normal'
	rand=str(np.random.randint(10000,100000))
	cmd="snana.exe NOFILE VERSION_PHOTOMETRY "+genversion+\
		" SNCCID_LIST "+cid_list+\
		" CUTWIN_CID 0 0 SNTABLE_LIST 'SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX 'OUT_TEMP_"+rand+\
		"' > OUT_TEMP_"+rand+".LOG"
	
	os.system(cmd)
	return(plotter,'OUT_TEMP_'+rand)

def main():
	parser = OptionParser()
	parser.add_option("--spec",action="store_true",dest="spec",default=False)
	parser.add_option("--lc",action="store_true",dest="lc",default=False)
	parser.add_option("--noclean",action="store_true",dest="noclean",default=False)
	parser.add_option("-i",action="store",type="string",dest="CID",default="None")
	parser.add_option("-b",action="store",type="float",dest='bin_size',default=0)
	parser.add_option("-v",action="store",type='string',dest='genversion',default=None)
	(options,args)=parser.parse_args()
	if options.CID=="None":
		raise RuntimeError("Need to define CID")
	if options.genversion is None:
		raise RuntimeError("Need to define genversion")
	
	plotter_choice,options.base_name=plot_cmd(options.genversion,options.CID)
	options.CID=options.CID.split(',')
	for cid in options.CID:
		if options.spec:
			plot_spec([cid],options.bin_size,options.base_name)
		elif options.lc:
			plot_lc([cid],options.base_name)
		else:
			plot_spec([cid],options.bin_size,options.base_name)
			plot_lc([cid],options.base_name)
	
	if not options.noclean:
		for x in glob.glob(options.base_name+'*'):
			os.remove(x)

if __name__=='__main__':
	main()
