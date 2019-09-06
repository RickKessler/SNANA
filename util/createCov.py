#!/usr/bin/env python
#
# Created May 23 2017 by D.Scolnic
# Static bins on Oct 10 2018 D.Brout
# Updated August 5 2019 by D.Brout
# Updated for SNANA Sep 5 2019 D.Brout

# Make a systematic covariance matrix from a folder of M0DIF files
#
#
# Usage:
#   createCov.py <inFile> [--covmatonly]
#
#def sysmat(base_output,fitopt='_',muopt='_',topdir='',do_each=0,sysfile='SYS.LIST',output_dir='COSMO',sysdefault=1,remove_extra=True,topfile=''):
# where <inFile> contains
#   COSMOMC_TEMPLATES: <path> #mandatory location of cosmomc template ini and sbatch files to be used
#   BASEOUTPUT:  <path>   # this is mandatory
#   TOPDIR:  <path>   # this must be filled in to know where to look for files
#   FITOPT: <option> # if you want to only include a specific FITOPT
#   MUOPT: <option> # if you want to only include a specific MUOPT
#  SYSFILE: <file> #if you want to specify the magnitude of each systematic
# SYSDEFAULT: <option> #Specify the default scale of each systematic - should be a number
# OUTPUTDIR: <path> #This is where output files are going
# REMOVE_EXTRA: <True> #Don't touch this for now
# TOPFILE: <path> # If you want to specify the base file for the cosmology (so not FITOPT000_MUOPT000)

# Outputs:
#    OUTPUTDIR/BASEOUTPUT.dataset
#    OUTPUTDIR/BASEOUTPUT_nosys.dataset
# OUTPUTDIR/lcparam_BASEOUTPUT.txt
# OUTPUTDIR/sys_BASE_OUTPUT.txt

import os
import sys
import numpy as np
import getpass
import time
import shutil
import string
import fnmatch

# globals

SNDATA_ROOT = os.environ['SNDATA_ROOT']
HOSTNAME    = os.environ['HOSTNAME']
NOW         = time.strftime("%c")
CWD         = os.getcwd()

def linef(file1,line1):
            co=0
            with open(file1, 'r') as inF:
                        for line in inF:
                                if line1 in line:
                                        return co
                                co=co+1

            print('Didnt find line in file, going to crash')
            stop
            
def parseLines(Lines,key,narg,vbose):
        # Lines is input array of lines in file
        # key is key to search
        # narg is number of args to return after key

        arg     = []
        rowList = Lines[np.char.startswith(Lines,key)]
        nrow    = len(rowList)
        print('parse', key, narg)
        if (( nrow == 1 )&(narg!=99)):
            if ( narg==1 ):
                arg = rowList[0].split()[1]
            else:
                arg = rowList[0].split()[1:narg+1]
        elif (( nrow > 1 )|(narg==99)):
          for row in rowList:
                    if (narg!=99):    
                                narg2=len(row.split())
                                arg.append(row.split()[1:narg2])
                    if (narg==99):
                                arg.append(row)
        if ( vbose > 0 ):
                print('\t ', key, arg)
                
        return(arg)

def dataset(output_dir,base_output,strex1,strex2,sys=1):
            # DILLON: I replaced getcwd() because I'm specifying full output
            g=open(output_dir+'/'+base_output+strex1+'.dataset','w+');   
            #g=open(os.getcwd()+'/'+output_dir+'/'+base_output+strex1+'.dataset','w+');
            g.write('name = JLA\n');
            #g.write('data_file = '+os.getcwd()+'/'+output_dir+'/lcparam_'+base_output+strex2+'.txt\n');           
            g.write('data_file = '+output_dir+'/lcparam_'+base_output+strex2+'.txt\n');
            #print 'data_file = '+os.getcwd()+'/'+output_dir+'/lcparam_'+base_output+strex2+'.txt\n'
            g.write('pecz = 0\n');  
            g.write('intrinsicdisp = 0\n'); 
            g.write('twoscriptmfit = F\n');
            g.write('scriptmcut = 10.0\n');
                  
            if (sys==1):
                        g.write('has_mag_covmat = T\n');
                        g.write('mag_covmat_file =  '+'/'+output_dir+'/sys_'+base_output+strex1+'.txt\n')                        
                        #g.write('mag_covmat_file =  '+os.getcwd()+'/'+output_dir+'/sys_'+base_output+strex1+'.txt\n')
            if (sys!=1): g.write('has_mag_covmat = F\n')
                    
            g.write('has_stretch_covmat = F\n');
            g.write('has_colour_covmat = F\n'); 
            g.write('has_mag_stretch_covmat = F\n');
            g.write('has_mag_colour_covmat = F\n'); 
            g.write('has_stretch_colour_covmat = F\n');
            g.close();
            return 2

def fullcosmo(base_output,file1,lc1,mat1,output_dir='COSMO'):
        from scipy.interpolate import interp2d
        headn=linef(file1,'zCMB')
        data1=np.genfromtxt(file1,skip_header=headn,names=True,comments='#')
        cid=np.genfromtxt(file1,skip_header=headn,usecols=(1),comments='#',dtype='str')[1:]
        z1 = data1['zHD'].astype(float)
        mu = data1['MU'].astype(float)
        muerr = data1['MUERR'].astype(float)
                        
        f1=open(output_dir+'/lcparam_'+base_output+'.txt','w') #this is the file for cosmomc
        f1.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec  \n')          #standard format
        for x in range(0,len(z1)):
                f1.write(cid[x]+' '+str(z1[x])+' '+str(z1[x])+' 0.0 '+str(mu[x]-19.35)+' '+str(muerr[x])+' 0 0 0 0 0 0 0 0 0 0 0\n')
        f1.close()
        g=open(output_dir+'/'+base_output+'.dataset','w'); h=open(output_dir+'/'+base_output+'_nosys.dataset','w');
        print('Shafer', lc1)
        ztemp1,whos = np.loadtxt(lc1, usecols=(1,2), unpack=True, dtype='str',skiprows=1)
        #stop
        #scount=sys1[0]
        #sys1=sys1.astype(float)
        ztemp1=ztemp1.astype(float)
        sys1 = np.loadtxt(mat1, unpack=True, dtype='str')
        scount=sys1[0]
        sys1=sys1.astype(float)
        bigmatmm=np.zeros((len(ztemp1), len(ztemp1)))+.000000
        co=1
        for x in range(0,len(ztemp1)):
                for y in range(0,len(ztemp1)):
                        bigmatmm[x,y]=sys1[co]
                        co=co+1
        
        gmm=open(output_dir+'/sys_'+base_output+'.txt','w')
        gmm.write(str(len(z1))+'\n')
        xvec=ztemp1; yvec=ztemp1
        f = interp2d(xvec, yvec, bigmatmm)
        #stop
        for x in range(0,len(z1)):
                linemm=''
                for y in range(0,len(z1)):
                        xx=np.argmin(np.absolute(z1[x]-ztemp1))
                        yy=np.argmin(np.absolute(z1[y]-ztemp1))
                        #big1=bigmatmm[xx-1:xx+1,yy-1:yy+1]
                        #xvec=ztemp1[xx-1:xx+1]; yvec=ztemp1[yy-1:yy+1]
                        #f = interp2d(xvec, yvec, big1)
                        #print 'mat',big1
                        #print 'xvec yvec', xvec, yvec
                        if ((z1[x]>1.312)&(z1[y]>1.312)):
                                    print('comparison', bigmatmm[xx,yy], f(z1[x],z1[y]))
                                    #stop
                                    #print z1[x], z1[y], ztemp1[xx], ztemp1[yy]
                        linemm=''
                        
                        #linemm=str("%.8f" % bigmatmm[xx,yy])
                        linemm=str("%.8f" % f(z1[x],z1[y]))
                        #linemm=str("%.8f" % bigmatmm[xx,yy])
                        gmm.write(linemm+'\n')
        gmm.close()
        
        dataset(output_dir,base_output,'','',sys=1)
        dataset(output_dir,base_output,'_nosys','',sys=0)
        #stop
        return 2

                                                                                                                                                    
                                                               
def avgmat(base_output,mat1,mat2,lc1,lc2,output_dir='COSMO'):
    import numpy as np
    import matplotlib.pyplot as plt
    import array
    import math
    import os
    from astropy import cosmology as cosmo
    from astropy.cosmology import FlatLambdaCDM
    import re
    cosmo2=FlatLambdaCDM(H0=70, Om0=0.3)
        
    list1, z1,mb1,mb1e = np.loadtxt(output_dir+'/lcparam_'+lc1+'.txt', usecols=(0, 1,4,5), unpack=True, dtype='string')
    z1 = z1.astype(float)
    mb1 = mb1.astype(float)
    mb1e = mb1e.astype(float)
    
    x=cosmo2.luminosity_distance(z1).value
    mu_syn1=5.0*(np.log10(x))+25.0-19.35
    mu1=mb1-mu_syn1
    
    list2, z2,mb2,mb2e = np.loadtxt(output_dir+'/lcparam_'+lc2+'.txt', usecols=(0, 1,4,5), unpack=True, dtype='string')
    
    z2 = z1 #using z1 so z lines up
    mb2 = mb2.astype(float)
    mb2e = mb2e.astype(float)
        
    x=cosmo2.luminosity_distance(z2).value
    mu_syn2=5.0*(np.log10(x))+25.0-19.35
    mu2=mb2-mu_syn2

    mua=(mu1+mu2)/2.0
    muae=(mb1e+mb2e)/2.0
    mua=mu_syn1+mua
    print(output_dir+'/lcparam_'+lc1+'.txt')
    print(output_dir+'/lcparam_'+lc2+'.txt')
    #stop
    #print z1
    #stop
    f1=open(output_dir+'/lcparam_'+base_output+'.txt','w') #this is the file for cosmomc    
    f1.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor \n')          #standard format
    for x in range(0,len(z1)):
        f1.write(str(list1[x])+' '+str(z1[x])+' '+str(z1[x])+' 0.0 '+str(mua[x])+' '+str(muae[x])+' 0 0 0 0 0 0 0 0 0 0 0 0\n')
    f1.close()
    print(output_dir+'/sys_'+mat1+'.txt')
    print(output_dir+'/sys_'+mat2+'.txt')
        
    sys1 = np.loadtxt(output_dir+'/sys_'+mat1+'.txt', unpack=True, dtype='string')
    sys2 = np.loadtxt(output_dir+'/sys_'+mat2+'.txt', unpack=True, dtype='string')
    scount=sys1[0]
    sys1=sys1.astype(float)
    sys2=sys2.astype(float)
    savg=(sys1+sys2)/2.0
    sys3=open(output_dir+'/sys_'+base_output+'.txt','w')
    sys3.write(str(scount)+'\n')
    for x in range(1,len(sys1)):
        sys3.write(str(savg[x])+'\n')
    sys3.close()
    dataset(output_dir,base_output,'','',sys=1);
    dataset(output_dir,base_output,'_nosys','',sys=0);
    print(output_dir+'/sys_'+base_output+'.txt')




def avgmat_Ngrid(base_output,mats,lcs,output_dir='COSMO'):
    import numpy as np
    import matplotlib.pyplot as plt
    import array
    import math
    import os
    from astropy import cosmology as cosmo
    from astropy.cosmology import FlatLambdaCDM
    import re
    cosmo2=FlatLambdaCDM(H0=70, Om0=0.3)


    lists, zs, mbs, mbes,xs,mu_syns,mus = [],[],[],[],[],[],[]
    
    for mat,lc in zip(mats,lcs):
                list1, z1,mb1,mb1e = np.loadtxt(output_dir+'/lcparam_'+lc+'.txt', usecols=(0, 1,4,5), unpack=True, dtype='string')
                z1 = z1.astype(float)
                mb1 = mb1.astype(float)
                mb1e = mb1e.astype(float)
                lists.append(list1)
                zs.append(z1)
                mbs.append(mb1)
                mbes.append(mb1e)

                x=cosmo2.luminosity_distance(zs[0]).value
                mu_syn1=5.0*(np.log10(x))+25.0-19.35
                mu1=mb1-mu_syn1

                xs.append(x)
                mu_syns.append(mu_syn1)
                mus.append(mu1)
                output_dir+'/lcparam_'+lc+'.txt'  
                
    #list2, z2,mb2,mb2e = np.loadtxt(output_dir+'/lcparam_'+lc2+'.txt', usecols=(0, 1,4,5), unpack=True, dtype='string')

    #z2 = z1 #using z1 so z lines up                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    #mb2 = mb2.astype(float)
    #mb2e = mb2e.astype(float)

    #x=cosmo2.luminosity_distance(z2).value
    #mu_syn2=5.0*(np.log10(x))+25.0-19.35
    #mu2=mb2-mu_syn2
    mua = np.mean(mus,axis=0)
    #print mua.shape
    muae = np.mean(mbes,axis=0)
    #print muae.shape
    #asdf
    mua = np.array(mu_syns[0]) + mua
    
    #mua=(mu1+mu2)/2.0
    #muae=(mb1e+mb2e)/2.0
    #mua=mu_syn1+mua
    #print output_dir+'/lcparam_'+lc1+'.txt'
    #print output_dir+'/lcparam_'+lc2+'.txt'
    #stop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    #print z1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    #stop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    f1=open(output_dir+'/lcparam_'+base_output+'.txt','w') #this is the file for cosmomc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    f1.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor \n')          #standard format                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    for x in range(0,len(zs[0])):
        f1.write(str(lists[0][x])+' '+str(zs[0][x])+' '+str(zs[0][x])+' 0.0 '+str(mua[x])+' '+str(muae[x])+' 0 0 0 0 0 0 0 0 0 0 0 0\n')
    f1.close()
    #print output_dir+'/sys_'+mat1+'.txt'
    #print output_dir+'/sys_'+mat2+'.txt'


    syss = []
    for mat,lc in zip(mats,lcs):
                
                sys1 = np.loadtxt(output_dir+'/sys_'+mat+'.txt', unpack=True, dtype='string')
                sys1=sys1.astype(float)
                syss.append(sys1)
    #print syss[0].shape
    savg = np.mean(syss,axis=0)
    #print savg.shape
    #asdf
    scount = syss[0][0]
    sys1 = syss[0]
    sys3=open(output_dir+'/sys_'+base_output+'.txt','w')
    sys3.write(str(scount)+'\n')
    for x in range(1,len(sys1)):
                sys3.write(str(savg[x])+'\n')
    sys3.close()   
    dataset(output_dir,base_output,'','',sys=1);
    dataset(output_dir,base_output,'_nosys','',sys=0);
    print(output_dir+'/sys_'+base_output+'.txt')





def sysmat(base_output,fitopt='_',muopt='_',topdir='',do_each=0,sysfile='NONE',output_dir='COSMO',sysdefault=1,remove_extra=True,covlines='',topfile='NONE',errscales='NONE',subdir='*'):
        
    import numpy as np
    import matplotlib.pyplot as plt
    import array
    import math
    import os
    from astropy import cosmology as cosmo
    from astropy.cosmology import FlatLambdaCDM
    import re
    if not output_dir: output_dir='COSMO'
    if not sysdefault: sysdefault=1
    if not remove_extra: remove_extra=True
            
    if (os.path.isdir(output_dir)==False): os.mkdir(output_dir)
    print(len(covlines))
    #stop
    if len(covlines)>1: sysnum=len(covlines)
    if covlines=='NONE': sysnum=0
    co=0
    sys_ratio=1
    print('subdir', subdir)
    print('topdir', topdir)
    print('fitopt', fitopt)
    print('muopt', muopt)
    print('base_output',base_output)
    print('ls '+topdir+'/'+subdir+'/*'+fitopt+'*'+muopt+'*'+'M0DIF')
    os.system('ls '+topdir+'/'+subdir+'/*'+fitopt+'*'+muopt+'*'+'M0DIF > '+base_output+'.list')
    print('subdir', subdir)
    print('topdir', topdir)
    #stop
    print('Created list of all M0DIF files called '+base_output+'.list')
    if os.path.isfile(base_output+'.list')==False:
        print('List file does not exist. No M0DIF files!!! This makes me sad!!! Im done here!!')
        return 0
    if (os.stat(base_output+'.list').st_size == 0):
        print('List exists but is empty. No M0DIF files!!! This makes me sad!!! Im done here!!')
        return 0

    if not os.path.exists(topdir+'/FITJOBS/FITJOBS_SUMMARY.LOG'):
        print(topdir+'/FITJOBS/FITJOBS_SUMMARY.LOG')        
        print('Log file not there. No M0DIF files!!! This makes me sad!!! Im done here!!')
        return 0

    if os.path.isfile(base_output+'.list'): file_lines=open(base_output+'.list','r').readlines()
    if os.path.isfile(topdir+'/FITJOBS/FITJOBS_SUMMARY.LOG'): log_lines=open(topdir+'/FITJOBS/FITJOBS_SUMMARY.LOG','r').readlines()    
    print(topdir+'/FITJOBS/FITJOBS_SUMMARY.LOG')

    filesize=len(file_lines) #read in number of M0DIF files
    print('Total of '+str(filesize)+' M0DIF files')

    MUOPT_var1=[]
    MUOPT_var2=[]
        
    FITOPT_var1=[]
    FITOPT_var2=[]

    SYSOPT_var1=[]
    SYSOPT_var2=[]
    SYSOPT_var3=[]
        
    INPDIR1=[]
    
    for xco in range(0,len(log_lines)):
        if 'MUOPT:' in log_lines[xco]: 
                mu_split=log_lines[xco].split()
                print(mu_split)
                MUOPT_var1=np.append(MUOPT_var1,'MUOPT'+mu_split[1])
                MUOPT_var2=np.append(MUOPT_var2,mu_split[2][1:-1])
                
        if 'INPDIR+:' in log_lines[xco]:
            mu_split=log_lines[xco].split()
            INPDIR1=np.append(INPDIR1,mu_split[1])

    
                        

    print(INPDIR1[0]+'/FITOPT.README')
    #stop
    #stop
    if os.path.isfile(INPDIR1[0]+'/FITOPT.README')==False:
         print('No FITOPT README in !!!'+INPDIR1[0]+'/FITOPT.README This makes me sad!!! Im done here!!')
         return 0
    if os.path.isfile(INPDIR1[0]+'/FITOPT.README'): fit_lines=open(INPDIR1[0]+'/FITOPT.README','r').readlines()

    for xco in range(0,len(fit_lines)):
        if 'FITOPT:' in fit_lines[xco]:
            mu_split=fit_lines[xco].split()
            FITOPT_var1=np.append(FITOPT_var1,'FITOPT'+mu_split[1])
            FITOPT_var2=np.append(FITOPT_var2,mu_split[2][1:-1])
    
                        
    if (((os.path.isfile(sysfile)&(sysfile!='NONE')&(errscales=='NONE'))|((sysfile=='NONE')&(errscales!='NONE')))):
      if ((os.path.isfile(sysfile)&(sysfile!='NONE')&(errscales=='NONE'))):
        if (((os.path.isfile(sysfile)==False)&(sysfile!='NONE'))):
           print('That '+ sysfile +' doesnt exist.  Grrrr.  Have to leave')
           stop
        sys_lines=open(sysfile,'r').readlines()
      if ((sysfile=='NONE')&(errscales!='NONE')):
        sys_lines=errscales
      print('syslines', sys_lines)
      #stop
      for xco in range(0,len(sys_lines)):
         if 'ERRSCALE:' in sys_lines[xco]:
             mu_split=sys_lines[xco].split()
             SYSOPT_var1=np.append(SYSOPT_var1,mu_split[1])
             SYSOPT_var2=np.append(SYSOPT_var2,mu_split[2])
             SYSOPT_var3=np.append(SYSOPT_var3,mu_split[3])
             
    if ((sysfile=='NONE')&(errscales=='NONE')):
             print('WARNING: All systematics have default scaling with no cuts.  This is really dangerous!')
             
             SYSOPT_var1=[]
    if ((sysfile!='NONE')&(errscales!='NONE')):
        print('You have a list of systematics in your inFile and in your included file.  That is one two many lists.  We have to stop')
        stop
        
    xco=0
    if (topfile!='NONE'): topfile=file_lines[xco][:-33]+topfile
    if ((topfile=='NONE')|(topfile=='')|(topfile=='None')): topfile=file_lines[xco][:-1]
    
    skipc=linef(topfile,'VARNAMES')
    if (topfile!=''):z1, mu1, mu1e   = np.loadtxt(topfile, usecols=(4,5,6), unpack=True, dtype='str', skiprows=skipc+1)
    if (topfile==''):z1, mu1, mu1e   = np.loadtxt(topfile, usecols=(4,5,6), unpack=True, dtype='str', skiprows=skipc+1)
    print('topfile', topfile)
    mu1=mu1.astype(float)
    mu1e=mu1e.astype(float)
    z1=z1.astype(float)
    #xxa=[mu1e<90]

    xxa=[mu1e<np.inf]#CHANGED BY DILLON HERE to get covmats all the same size for multiple sims
    z1=z1[xxa]
    mu1=mu1[xxa]
    mu1e=mu1e[xxa]
    cosmo2=FlatLambdaCDM(H0=70, Om0=0.3)
    x=cosmo2.luminosity_distance(z1).value
    mu_syn=5.0*(np.log10(x))+25.0-19.35
    mu_syn1=mu_syn+mu1
    
    f1=open(output_dir+'/lcparam_'+base_output+'.txt','w') #this is the file for cosmomc
    f1.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor \n')          #standard format
    for x in range(0,len(z1)):
        f1.write(str(x)+' '+str(z1[x])+' '+str(z1[x])+' 0.0 '+str(mu_syn1[x])+' '+str(mu1e[x])+' 0 0 0 0 0 0 0 0 0 0 0 0\n')
    f1.close()
    bigmatmm=np.zeros((len(z1), len(z1),sysnum+1))+.000000
                              
    logf=open(output_dir+'/'+base_output+'.log','w')    
    for xco in range(0,len(file_lines)):
        print(file_lines[xco].split('_')[-2], file_lines[xco].split('_')[-1][:-7])
        #SALT2mu_SNLS+SDSS+LOWZ+PS1_Scolnic2+HST/DS17/SALT2mu_FITOPT000_MUOPT000.M0DIF
        #stop
        xx1=(FITOPT_var1==file_lines[xco].split('_')[-2])
        xx2=(MUOPT_var1==file_lines[xco].split('_')[-1][:-7])
        skipc=linef(file_lines[xco][:-1],'VARNAMES')
        z2, mu2, mu2e   = np.loadtxt(file_lines[xco][:-1], usecols=(4,5,6), unpack=True, dtype='str', skiprows=skipc+1)
        print(file_lines[xco][:-1])
        mu2=mu2.astype(float)
        mu2e=mu2e.astype(float)
        z2=z2.astype(float)
        #xxa=[mu2e<900000]
        z2=z2[xxa]
        mu2=mu2[xxa]
        mu2e=mu2e[xxa]
                        
        cosmo2=FlatLambdaCDM(H0=70, Om0=0.3)
        x=cosmo2.luminosity_distance(z1).value            
        mu_syn2=5.0*(np.log10(x))+25.0-19.35
        print(len(z1), len(z2), len(mu1), len(mu_syn2), len(mu2))
        #35 32 35 35 32
        mu_syn2=mu_syn2+mu2
        xxb=((z1==0)|(z2==0))
        if len(z2[xxb])>0:
                    
                    mu_syn2[xxb]=mu_syn1[xxb]
        sys_ratio=float(sysdefault)
        print('sysopt', SYSOPT_var1)
        print(FITOPT_var2)
        #stop
        if len(SYSOPT_var1)>0:
                comatch=0
                for y1 in range(0,len(SYSOPT_var1)):
                   filtered1 = fnmatch.filter([FITOPT_var2[xx1][0]], SYSOPT_var1[y1])
                   filtered2 = fnmatch.filter([MUOPT_var2[xx2][0]], SYSOPT_var2[y1])
                   if ((len(filtered1)>0)&(len(filtered2)>0)):
                        print('sys', SYSOPT_var3)       
                        sys_ratio=float(SYSOPT_var3[y1])
                        #print sys_ratio
                        #stop
                        print('Have a systematic from '+str(SYSOPT_var1[y1])+str(SYSOPT_var2[y1])+' of '+str(SYSOPT_var3[y1]))
                        logf.write('Have a systematic from '+str(SYSOPT_var1[y1])+str(SYSOPT_var2[y1])+' of '+str(SYSOPT_var3[y1])+'\n')
                        #stop
                        if (comatch>0):
                           print('WARNING you have had multiple systematics match up!!! That is bad')
                        comatch=comatch+1
                         
                        #if ((np.amax(np.absolute(z1-z2)/z1)>0.1)&(sys_ratio>0)):
                                    
                        #            print z1-z2
                        #            print np.absolute(z1-z2)/z1
                        #            print 'There is a misalignment of z bins!!! We have to stop!'
                        #            print file_lines[xco][:-1]
                        #            print z1[0], z2[0], z1[0]
                        #            stop
                                                                                                                                                                                                           
        #if 'SALT2' in FITOPT_var2[xx1][0]:
                    #print sys_ratio
                    #stop
        distm=np.zeros((1))
        dm2=mu_syn1-mu_syn2
        dm2=np.multiply(dm2,sys_ratio)
        dm2t=np.matrix(dm2)
        dm2t=dm2t.T
        dmm=dm2t*np.matrix(dm2)
        
                    #stop
        x=0
        bigmatmm[:,:,x]=np.add(bigmatmm[:,:,x],np.multiply(dmm,1.0))
        print('bigmat', bigmatmm[1,1,0], dmm[1,1])
        if dmm[1,1]>.3:
                print(file_lines[xco])
                stop
        print('covlines all', covlines)
        print('sysnum', sysnum+1)
        #stop
        for x in range(1,sysnum+1):
                print(x)
                print('covlines', covlines[x-1])    
                syscheck1=covlines[x-1][1][1:-1].split(',')[0]; syscheck2=covlines[x-1][1][1:-1].split(',')[1]
                sys_flag1=False
                sys_flag2=False
                if (syscheck1[0]=='-'): sys_flag1=(syscheck1[1:] not in FITOPT_var2[xx1][0])
                if (syscheck1[0]=='+'): sys_flag1=(syscheck1[1:] in FITOPT_var2[xx1][0])
                if (syscheck1[0]=='='): sys_flag1=(syscheck1[1:]==FITOPT_var2[xx1][0])
                if (syscheck2[0]=='-'): sys_flag2=(syscheck2[1:] not in MUOPT_var2[xx2][0])
                if (syscheck2[0]=='+'): sys_flag2=(syscheck2[1:] in MUOPT_var2[xx2][0])
                if (syscheck2[0]=='='): sys_flag2=(syscheck2[1:]==MUOPT_var2[xx2][0])
                if (syscheck1[0]=='-'):
                            print(sys_flag1)
                            print(sys_flag2)
                            print(FITOPT_var2[xx1][0], MUOPT_var2[xx2][0], (sys_flag1)&(sys_flag2))
                            #stop
                if ((sys_flag1)&(sys_flag2)):
                        logf.write(FITOPT_var2[xx1][0]+' '+MUOPT_var2[xx2][0]+' '+syscheck1[0:]+' '+syscheck2[0:]+' '+str(x)+' '+str(sys_ratio)+' \n')
                        bigmatmm[:,:,x]=np.add(bigmatmm[:,:,x],np.multiply(dmm,1.0))
                                
        co=co+1
    for z in range(0,(sysnum+1)):
            gmm=open(output_dir+'/sys_'+base_output+'_'+str(z)+'.txt','w')
            gmm.write(str(len(z1))+'\n')
            for x in range(0,len(z1)):
                    linemm=''            
                    for y in range(0,len(z1)):
                            
                            linemm=''
                            linemm=str("%.8f" % bigmatmm[x,y,z])
                            gmm.write(linemm+'\n')
            gmm.close()
    
    dataset(output_dir,base_output,'_0','',sys=1);
    dataset(output_dir,base_output,'_nosys','',sys=0);
    print('sysnum', sysnum)
    #stop
    for z in range(1, sysnum+1):
            
            temp=dataset(output_dir,base_output,'_'+str(z),'',sys=1)
            
            
            #stop
    print('just did it')
    logf.close()
    return 2

class FILE_INFO:
        def __init__(self,filename):
                print('   Parse kcor input file: ', filename)
                
                filename_expandvars = os.path.expandvars(filename)
                # check local dir first; then check $SNDATA_ROOT/kcor
                if os.path.isfile(filename_expandvars):
                        fname_local = filename_expandvars
                else:
                        fname_local = TOPDIR_KCOR + '/' + filename_expandvars
                        
                # open file file and read all lines into Lines
                f = open(fname_local,"rt")
                Lines  = np.array(f.readlines())
                self.COSMOMC_TEMPLATES=parseLines(Lines,'COSMOMC_TEMPLATES:',1,1)
                self.BASEOUTPUT=parseLines(Lines,'BASEOUTPUT:',1,1)   
                self.TOPDIR=parseLines(Lines,'TOPDIR:',1,1)   
                self.FITOPT=parseLines(Lines,'FITOPT:',1,1)   
                self.MUOPT=parseLines(Lines,'MUOPT:',1,1)   
                self.SYSFILE=parseLines(Lines,'SYSFILE:',1,1)
                self.USEFILE=parseLines(Lines,'USE_SYSFILE:',1,1)                
                self.SYSDEFAULT=parseLines(Lines,'SYSDEFAULT:',1,1)   
                self.OUTPUTDIR=parseLines(Lines,'OUTPUTDIR:',1,1)
                self.ROOTDIR=parseLines(Lines,'ROOTDIR:',1,1)
                self.TOPFILE=parseLines(Lines,'TOPFILE:',1,1)
                self.COVOPT=parseLines(Lines,'COVOPT:',1,1)
                self.ERRSCALE=parseLines(Lines,'ERRSCALE:',99,1)
                self.MAKEINI=parseLines(Lines,'MAKEINI:',1,1)
                self.SUBDIR=parseLines(Lines,'SUBDIR:',1,1)
                                                
# =========================
# ======= MAIN ============
# =========================

import os.path
def makeini (file_root, outputdir,baseoutput,base,BASE_INI_AND_BATCH,extra=0,rootdir='',datasetnum=0):
            #dataset=outputdir+'/'+baseoutput+'.dataset'
            dataset='%s_%d.dataset'%(baseoutput,datasetnum)
            #dvin_nosn_ocmb_omol.ini
            print('we are making ini files!')
            svec=['snonly_omw','omw','wwa','omol','snonly_omol','sn_omw_validation']
            gvec=['','','bao_','obao_','ocmb_','nosn_ocmb_']
            if (extra==0):
                        ex1=0
                        ex2=0
            if (extra==1):
                        ex1=4
                        ex2=4                        
            for i in range(0,2+ex1):
                        for j in range(0,2+ex2):
                                    #print base+'_'+gvec[j]+svec[i]+'.ini'
                                    #print os.path.isfile(base+'_'+gvec[j]+svec[i]+'.ini')
                                    #stop
                                    #print(BASE_INI_AND_BATCH+'/'+base+'_'+gvec[j]+svec[i]+'.ini')
                                    if os.path.isfile(BASE_INI_AND_BATCH+'/'+base+'_'+gvec[j]+svec[i]+'.ini') :
                                                g=open(BASE_INI_AND_BATCH+'/'+base+'_'+gvec[j]+svec[i]+'.ini','r')
                                                h=open(outputdir+'/'+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'.ini','w')
                                                with open(BASE_INI_AND_BATCH+'/'+base+'_'+gvec[j]+svec[i]+'.ini','r') as f:
                                                            content = f.readlines()
                                                for x in content:
                                                            h.write(x)
                                                h.write('\nfile_root='+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'\n')
                                                h.write('jla_dataset='+dataset+'\n')
                                                if not os.path.exists(rootdir):
                                                            os.mkdir(rootdir)
                                                h.write('root_dir = '+rootdir+'\n')
                                                h.close()
                                                g.close()
                                                with open(BASE_INI_AND_BATCH+'/'+base+'_temp.sbatch','r') as f:
                                                            content = f.readlines()
                                                h=open(outputdir+'/'+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'.sbatch','w')
                                                for x in content:
                                                            h.write(x)
                                                #h.write('timeout 126000 mpirun /project/rkessler/SN/CosmoMC/v01 '+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'.ini\n')
                                                h.write('module unload openmpi\n')
                                                h.write('module load intelmpi/5.1+intel-16.0\n')
                                                h.write('timeout 126000 mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc '+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'.ini\n')
                                                h.write('if [ $? -eq 124 ]; then\n')
                                                h.write('    sleep 120\n')
                                                h.write('    sbatch '+file_root+'_'+gvec[j]+svec[i]+'_'+str(int(datasetnum))+'.sbatch\n')
                                                h.write('fi\n')
                                                h.close()
                                                                                    
                                                                                                                                                                        
                                                                                                
if __name__ == "__main__":
                                                                                                            
             # parse input argument(s)
  covmatonly = False
  if ( len(sys.argv) < 2 ):
          sys.exit("Must give INFILE argument\n-->ABORT")
  else:
            if sys.argv[1][0] == '-':
                        sys.exit("Must give INFILE argument\n-->ABORT")
            else:
                        INFILE = sys.argv[1]
                        print('Input file: ', INFILE)
            try:
                        if sys.argv[2] == '--covmatonly':
                                    covmatonly = True

  print('SNDATA_ROOT = ', SNDATA_ROOT)
        
  FileInfo = FILE_INFO(INFILE)
  print(FileInfo.ERRSCALE)
  
  print('COSMOMC_TEMPLATES',FileInfo.COSMOMC_TEMPLATES)
  print('BASE OUTPUT', FileInfo.BASEOUTPUT)
  print('TOPDIR', FileInfo.TOPDIR)      
  print('SYSFILE', FileInfo.SYSFILE)
  print('SYSDEFAULT', FileInfo.SYSDEFAULT)
  print('OUTPUTDIR', FileInfo.OUTPUTDIR)
  print('ROOTDIR', FileInfo.ROOTDIR)
  print('COVLINES', FileInfo.COVOPT)
  print('ERRSCALE', FileInfo.ERRSCALE)
  print('TOPFILE', FileInfo.TOPFILE)
 
  if not FileInfo.BASEOUTPUT: print('We need a BASE OUTPUT'); stop
  if not FileInfo.TOPDIR: print('We need a TOPDIR'); stop
  if not FileInfo.SYSFILE: FileInfo.SYSFILE='NONE'
  if not FileInfo.SYSDEFAULT: FileInfo.SYSDEFAULT='NONE'
  if not FileInfo.OUTPUTDIR: print('no OUTPUTDIR specified so making it COSMO/'); FileInfo.OUTPUTDIR='COSMO'
  if not FileInfo.ROOTDIR: print('no ROOTDIR. please provide a rootdir for your chains\n-->ABORT'); sys.exit()
  if not FileInfo.COVOPT:  FileInfo.COVOPT='NONE'
  if not FileInfo.ERRSCALE:  FileInfo.ERRSCALE='NONE'
  if not FileInfo.TOPFILE:  FileInfo.TOPFILE='NONE'
  if not FileInfo.SUBDIR:  FileInfo.SUBDIR='*'
  if not FileInfo.MAKEINI:  FileInfo.MAKEINI='NONE'
    
  #if not FileInfo.MAKEINI:  FileInfo.MAKEINI='NONE'
  
  #stop
  sysmat(FileInfo.BASEOUTPUT,topdir=FileInfo.TOPDIR,sysfile=FileInfo.SYSFILE,sysdefault=FileInfo.SYSDEFAULT, output_dir=FileInfo.OUTPUTDIR, covlines=FileInfo.COVOPT,errscales=FileInfo.ERRSCALE,topfile=FileInfo.TOPFILE,subdir=FileInfo.SUBDIR)
  #print FileInfo.MAKEINI
  #stop
  print(FileInfo.OUTPUTDIR)
  #DILLON: I'm editing here for giving full outputdir path not relative to cwd
  f = open('/'.join(FileInfo.OUTPUTDIR.split('/')[:-1])+'/covopt.dict','w')
  if (FileInfo.MAKEINI!='NONE'):
              for d in range(len(FileInfo.COVOPT)+1):
                          makeini(FileInfo.MAKEINI,FileInfo.OUTPUTDIR,FileInfo.BASEOUTPUT,'dvin',FileInfo.COSMOMC_TEMPLATES,extra=1,rootdir=FileInfo.ROOTDIR,datasetnum=d)
                          if d == 0:
                                      covwrite = 'ALLSYS'
                          else:
                                      covwrite = FileInfo.COVOPT[d-1][0].replace('[','').replace("'",'').replace(']','')
                          f.write('%d\t%s\n'%(d,covwrite))
  f.close()
  #print "\n Done parsing ", nkcor, " kcor-input files "
                                                                                                        
  # change filter char
  #change_filterChar(versionInfo,kcorInfo)
                                                                                                                                                                                         
