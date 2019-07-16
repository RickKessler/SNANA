#!/usr/bin/env python
#
#June 2019 J. Pierel
#Plotter tool for SNANA LCs and Spectra

from __future__ import print_function
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os,glob,math,sys,textwrap
from optparse import OptionParser
from scipy.interpolate import interp1d


__band_order__=np.append(['u','b','g','r','i','z','y','j','h','k'],
	 [x.upper() for x in ['u','b','g','r','i','z','y','j','h','k']])

def read_spec(cid,base_name):
	names=['wave','flux','fluxerr','tobs','mjd']
	id_to_obs=dict([])
	mjds=[]
	with open(base_name+".SPECLIST.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		else:
			id_to_obs[int(temp[varnames.index('ID')])]=float(temp[varnames.index('TOBS')])
			mjds.append(float(temp[varnames.index('MJD')]))

	sn={k:[] for k in names}

	with open(base_name+".SPECPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	temp_id=None
	mjd_ind=0
	for line in dat:
		temp=line.split()
		
		
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		elif len(temp)>0 and b'OBS:' in temp and\
			 str(temp[varnames.index('CID')].decode('utf-8'))in cid:
			if temp_id is None:
				temp_id=int(temp[varnames.index('ID')])
			if temp_id!=int(temp[varnames.index('ID')]):	
				mjd_ind+=1
			temp_id=int(temp[varnames.index('ID')])
			sn['wave'].append((float(temp[varnames.index('LAMMAX')])+float(temp[varnames.index('LAMMIN')]))/2.)
			sn['flux'].append(float(temp[varnames.index('FLAM')]))
			sn['fluxerr'].append(float(temp[varnames.index('FLAMERR')]))
			sn['tobs'].append(id_to_obs[int(temp[varnames.index('ID')])])
			sn['mjd'].append(mjds[mjd_ind])
	sn={k:np.array(sn[k]) for k in sn.keys()}
	return(sn)
def read_lc(cid,base_name,plotter_choice):
	names=['time','flux','fluxerr','filter','chi2']
	peak=None
	sn={k:[] for k in names} 
	fit={k:[] for k in ['time','flux','filter']}
	with open(base_name+".LCPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	fitted=False
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
				fitted=True
				fit['time'].append(float(temp[varnames.index('Tobs')]))
				fit['flux'].append(float(temp[varnames.index('FLUXCAL')]))
				fit['filter'].append(str(temp[varnames.index('BAND')].decode('utf-8')))
	if fitted and plotter_choice=='salt2':
		with open(base_name+".FITRES.TEXT",'rb') as f:
			dat=f.readlines()
		for line in dat:
			temp=line.split()
			if len(temp)>0 and b'VARNAMES:' in temp:
				varnames=[str(x.decode('utf-8')) for x in temp]
			elif len(temp)>0 and b'SN:' in temp and str(temp[varnames.index('CID')].decode('utf-8')) in cid: 
				fit['params']={p:(float(temp[varnames.index(p)]),float(temp[varnames.index(p+'ERR')])) for p in ['x0','x1','c']}
				break
	sn={k:np.array(sn[k]) for k in sn.keys()}
	fit={k:np.array(fit[k]) if k !='params' else fit['params'] for k in fit.keys()}
	if len(fit['filter'])>0:
		fits={k:interp1d(fit['time'][fit['filter']==k],
					 fit['flux'][fit['filter']==k]) for k in np.unique(fit['filter'])}
		fits['params']=fit['params']
	else:
		fits=[]
	return(sn,fits,peak)

def plot_spec(cid,bin_size,base_name,noGrid):
	sn=read_spec(cid,base_name)

	if len(sn['tobs'])==0:
		return []
	if len(np.unique(sn['tobs']))>1:
		figs=[]
		m=0
		for nfig in range(int(math.ceil(len(np.unique(sn['tobs']))/4.))):
			fig,ax=plt.subplots(nrows=min(len(np.unique(sn['tobs'])),4),ncols=1,figsize=(8,8),sharex=True)
			ax[0].set_title('SN%s'%cid[0],fontsize=16)
			for j in range(min(len(np.unique(sn['tobs'])[m:]),4)):
				
				temp_sn=np.where(sn['tobs']==np.unique(sn['tobs'])[m])[0]
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
				
				if np.unique(sn['mjd'])[j]<0:
					spec_label='HOST'
				else:
					spec_label='SN:%.2f'%np.unique(sn['tobs'])[j]
				ax[j].plot(binned_wave,binned_flux,color='k',label=spec_label)
				ylim=ax[j].get_ylim()
				ax[j].fill_between(binned_wave,binned_flux-binned_fluxerr,binned_flux+binned_fluxerr,
							 color='r',alpha=.3,label=r'$1\sigma$ Error')
				ax[j].plot([binned_wave[0],binned_wave[-1]],[0,0],'k--',alpha=.5)
				ax[j].set_ylim(ylim)
				ax[j].legend(fontsize=12)

			
				ax[j].set_ylabel('Flux',fontsize=16)
				if not noGrid:
					ax[j].grid()
				m+=1
			ax[j].set_xlabel('Observer Frame Wavelength ($\AA$)',fontsize=16)
			
			figs.append(fig)
			plt.close()
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
		ind=np.argsort(binned_wave)	
		plt.plot(binned_wave[ind],binned_flux[ind],color='k',label='TOBS:%.2f'%np.unique(sn['tobs'])[0])
		ylim=plt.ylim()
		plt.fill_between(binned_wave[ind],binned_flux[ind]-binned_fluxerr[ind],binned_flux[ind]+binned_fluxerr[ind],
						 color='r',alpha=.3,label=r'$1\sigma$ Error')
		plt.plot([binned_wave[0],binned_wave[-1]],[0,0],'k--',alpha=.5)
		plt.ylim(ylim)
		plt.legend(fontsize=12)
		plt.xlabel('Observer Frame Wavelength ($\AA$)',fontsize=16)
		plt.ylabel('Flux',fontsize=16)
		plt.title('SN%s'%cid[0],fontsize=16)
		if not noGrid:
			plt.grid()
		figs=[fig]
		plt.close()
	#plt.savefig('SNANA_SPEC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)
	return(figs)

def plot_lc(cid,base_name,noGrid,plotter_choice):
	sn,fits,peak=read_lc(cid,base_name,plotter_choice)
	if len(sn['time'])==0:
		return []
	rows=int(math.ceil(len(np.unique(sn['filter']))))
	figs=[]
	all_bands=np.append([x for x in __band_order__ if x in np.unique(sn['filter'])],
                        [x for x in np.unique(sn['filter']) if x not in __band_order__])
	
	j=0
	minx=np.min(sn['time'])
	maxx=np.max(sn['time'])
	if minx<0:
		minx=min(minx*1.1,minx-5)
	else:
		minx=min(minx*.9,minx-5)
	if maxx<0:
		maxx=max(maxx*.9,maxx+5)
	else:
		maxx=max(maxx*1.1,maxx+5)
	xlims=(minx,maxx)
	sharedx=True
	for nfig in range(int(math.ceil(rows/4.))): 
		fig,ax=plt.subplots(nrows=min(len(all_bands),4),ncols=1,figsize=(8,8),sharex=sharedx)
		ax[0].set_title('SN%s'%cid[0],fontsize=16)
		fit_print=False
		for i in range(min(len(all_bands[j:]),4)):
			temp_sn={k:sn[k][np.where(sn['filter']==all_bands[j])[0]] for k in sn.keys()}
			chi2=np.mean(temp_sn['chi2'])
			if chi2>0:
				lab=r'%s: $\chi^2$=%.1f'%(all_bands[j],np.mean(temp_sn['chi2']))
				leg_size=10
			else:
				lab=all_bands[j]
				leg_size=12
			
			ax[i].errorbar(temp_sn['time'],temp_sn['flux'],yerr=temp_sn['fluxerr'],
						  fmt='.',markersize=8,color='k',
						  label=lab)
			if len(fits)>0:
				fit_time=np.arange(temp_sn['time'][0],temp_sn['time'][-1],1)
				ax[i].plot(fit_time,fits[all_bands[j]](fit_time),color='r',label='Best Fit',linewidth=3)
				if not fit_print:
					ax[i].annotate('\n'.join([r'$%s: %.2e\pm%.2e$'%(fit_key,fits['params'][fit_key][0],
						fits['params'][fit_key][1]) for fit_key in fits['params'].keys()]),xy=(.02,.65),xycoords='axes fraction',fontsize=6)
				fit_print=True
			ax[i].legend(fontsize=leg_size)
			ax[i].set_ylabel('Flux',fontsize=16)
			ax[i].set_ylim((-.1*np.max(temp_sn['flux']),1.1*np.max(temp_sn['flux'])))
			if not noGrid:
				ax[i].grid()
			j+=1
			#i+=1
		for k in range(i+1,min(len(all_bands),4)):
			fig.delaxes(ax[k])
		ax[i].tick_params(axis='x',labelbottom=True,bottom=True)
		ax[i].set_xlabel('MJD-%.2f'%peak,fontsize=16)
		ax[i].set_xlim(xlims)
		figs.append(fig)
		plt.close()
	#fig.text(0.5, 0.02, 'Time (Rest Frame Days)', ha='center',fontsize=16)
	#fig.text(0.04, .5, 'Flux', va='center', rotation='vertical',fontsize=16)
		
	#plt.savefig('SNANA_LC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)
	return(figs)

def plot_cmd(genversion,cid_list,nml):
	if nml is not None:
		if os.path.splitext(nml)[1]!='.NML':
			nml=os.path.splitext(nml)[0]+'.NML'
		plotter='salt2'
	else:
		plotter='normal'
	rand=str(np.random.randint(10000,100000))
	if nml is not None:
		cmd="snlc_fit.exe "+nml+" VERSION_PHOTOMETRY "+genversion+\
			" SNCCID_LIST "+cid_list+\
			" CUTWIN_CID 0 0 SNTABLE_LIST 'FITRES(text:key) SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX 'OUT_TEMP_"+rand+\
			"' > OUT_TEMP_"+rand+".LOG"
	else:
		cmd="snana.exe NOFILE VERSION_PHOTOMETRY "+genversion+\
			" SNCCID_LIST "+cid_list+\
			" CUTWIN_CID 0 0 SNTABLE_LIST 'SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX 'OUT_TEMP_"+rand+\
			"' > OUT_TEMP_"+rand+".LOG"
	
	os.system(cmd)
	with open('OUT_TEMP_'+rand+'.LOG','r+') as f:
		content=f.read()
		f.seek(0,0)
		f.write('SNANA COMMAND:\n\n'+textwrap.fill(cmd,80)+'\n'+content)
	if len(glob.glob('OUT_TEMP_'+rand+'*.TEXT'))==0:
		print("There was an error in retrieving your SN")
		sys.exit()
	return(plotter,'OUT_TEMP_'+rand)

def main():

	parser = OptionParser()
	parser.add_option("--spec",help='Plot only spectra',action="store_true",dest="spec",default=False)
	parser.add_option("--lc",help='Plot only LC',action="store_true",dest="lc",default=False)
	parser.add_option("--noclean",help='Leave intermediate files for debugging',action="store_true",dest="noclean",default=False)
	parser.add_option("-i",help='CID(s) as comma separated list',action="store",type="string",dest="CID",default="None")
	parser.add_option("-b",help='Bin size for spectral plotting',action="store",type="float",dest='bin_size',default=0)
	parser.add_option("-v",help='Version',action="store",type='string',dest='version',default=None)
	parser.add_option("-f",help='.NML filename',action="store",type='string',dest='nml_filename',default=None)
	parser.add_option("--silent",help="Do not print anything",action="store_true",dest="silent",default=False)
	parser.add_option("--nogrid",help="Do add a grid to the plots.",action="store_true",dest="noGrid",default=False)
	#parser.add_option("--help",action="store_true",dest='help',default=False)
	(options,args)=parser.parse_args()

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit()
	if options.CID=="None":
		raise RuntimeError("Need to define CID")
	if options.version is None:
		raise RuntimeError("Need to define genversion")
	
	plotter_choice,options.base_name=plot_cmd(options.version,options.CID,options.nml_filename)
	options.CID=options.CID.split(',')
	filename=options.version+'.pdf'
	num=0
	if os.path.exists(filename):
		filename=os.path.splitext(filename)[0]+'_'+str(num)+'.pdf'
	while os.path.exists(filename):
		num+=1
		filename=os.path.splitext(filename)[0][:-1]+str(num)+'.pdf'
	with PdfPages(filename) as pdf:
		for cid in options.CID:
			if not options.silent:
				print("Plotting SN %s"%cid)
			if options.spec:
				figs=plot_spec([cid],options.bin_size,options.base_name,options.noGrid)
				for f in figs:
					pdf.savefig(f)
			elif options.lc:
				figs=plot_lc([cid],options.base_name,options.noGrid,plotter_choice)
				for f in figs:
					pdf.savefig(f)
			else:
				figs=plot_spec([cid],options.bin_size,options.base_name,options.noGrid)
				for f in figs:
					pdf.savefig(f)
				figs=plot_lc([cid],options.base_name,options.noGrid,plotter_choice)
				for f in figs:
					pdf.savefig(f)
	
	if not options.noclean:
		for x in glob.glob(options.base_name+'*'):
			os.remove(x)
	if not options.silent:
		print('Done.')

if __name__=='__main__':
	main()
