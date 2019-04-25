#!/usr/bin/env python
# D. Jones - 1/7/16

# map FITRES names to plot names
histvardict = {'SNRMAX1':'max SNR',
               'zHD':'$z_{CMB}$',
               'x1':'$X_1$',
               'c':'$C$',
               'PKMJDERR':'$\sigma_{pkMJD}$',
               'x1ERR':'$\sigma_{x_1}$',
               'cERR':'$\sigma_C$',
               'MURES':'Hubble Residual',
               'x1vzHD':'mean $X_1$',
               'cvzHD':'mean $C$',
               'x1ERRvzHD':'$\sigma_{X_1}$',
               'cERRvzHD':'$\sigma_C$',
               'mBvzHD':'$m_B$',
               'mBERRvzHD':'$\sigma_{m_B}$',
               'SNRMAX1vmB':'max SNR'}

class txtobj:
    def __init__(self,filename):
        import numpy as np
        fin = open(filename,'r')
        lines = fin.readlines()
        for l in lines:
            if l.startswith('VARNAMES:'):
                l = l.replace('\n','')
                coldefs = l.split()
                break
        reader = [x.split() for x in lines if x.startswith('SN:')]

        i = 0
        for column in zip(*reader):
            try:
                self.__dict__[coldefs[i]] = np.array(column[:]).astype(float)
            except:
                self.__dict__[coldefs[i]] = np.array(column[:])
            i += 1

class ovhist:
    def __init__(self):
        self.clobber = False
        self.verbose = False

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        # The basics
        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
        parser.add_option('--clobber', default=False, action="store_true",
                          help='overwrite output file if it exists')
        parser.add_option('--cutwin',default=[],
                          type='string',action='append',
                          help='parameter range for specified variable',nargs=3)
        parser.add_option('--defaultcuts',default=False,action='store_true',
                          help='make the default sample cuts')
        parser.add_option('--x1cellipse',default=False,action='store_true',
                          help='elliptical default cut, not box; use in conjunction w/ defaultcuts')
        parser.add_option('--alpha',type='float',default=0.147,
                          help='SALT2 alpha for computing distances (default=%default)')
        parser.add_option('--beta',type='float',default=3.13,
                          help='SALT2 beta (default=%default)')
        parser.add_option('--dataM',type='float',default=-19.3,
                          help='SALT2 M parameter for data only (default=%default)')
        parser.add_option('--simM',type='float',default=-19.3,
                          help='SALT2 M parameter for simulations only (default=%default)')
        parser.add_option('--sigint',type='float',default=0.115,
                          help='Intrinsic Dispersion for computing distance errors (default=%default)')

        parser.add_option('-o','--outfile',type='string',default='',
                          help='output figure file (default=ovplot_[varname].pdf)')
        parser.add_option('--interact', default=False, action="store_true",
                          help='open up the figure interactively')
        parser.add_option('--nbins',type='int',default=40,
                          help='number of histogram bins (default=%default)')
        parser.add_option('--bins',type='float',default=(None,None,None),
                          help="""3 arguments specifying number of bins and the min. and max. of the histogram range.  
                          For multi-plot mode, use the --cutwin options instead.""",nargs=3)
        parser.add_option('--ylog', default='', type="string",
                          help='log-scaled y-axis for specified variables, use "all" for everything ')
        parser.add_option('--scaleb4cuts', default=False, action="store_true",
                          help='if set, scale sims to data before making cuts')
        parser.add_option('--ylim', default=(None,None), type='float',
                          help='custom y-axis limits',nargs=2)
        parser.add_option('--journal', default=False, action="store_true",
                          help='make journal figure')
        parser.add_option('--nplots', default=(None,None), type='int',
                          help='number of x,y plots on page (for journal option)',nargs=2)

        
        return(parser)

    def main(self,datafile,simfile):
        data = txtobj(datafile)
        sim = txtobj(simfile)

        # getting distance modulus is slow, so don't do it unless necessary
        getMU = False
        if len(self.options.cutwin):
            for cutopt in self.options.cutwin:
                if 'MU' in cutopt[0]: getMU = True
        for h in self.options.histvar:
            if 'MU' in h: getMU = True
                
        if 'MU' in self.options.histvar or getMU:
            if not 'MU' in data.__dict__:
                data.MU,data.MUERR = salt2mu(x1=data.x1,x1err=data.x1ERR,c=data.c,cerr=data.cERR,mb=data.mB,mberr=data.mBERR,
                                             cov_x1_c=data.COV_x1_c,cov_x1_x0=data.COV_x1_x0,cov_c_x0=data.COV_c_x0,
                                             alpha=self.options.alpha,beta=self.options.beta,
                                             x0=data.x0,sigint=self.options.sigint,z=data.zHD,M=self.options.dataM)
                from astropy.cosmology import Planck13 as cosmo
                if not 'MURES' in data.__dict__:
                    data.MURES = data.MU - cosmo.distmod(data.zHD).value
            if not 'MU' in sim.__dict__:
                sim.MU,sim.MUERR = salt2mu(x1=sim.x1,x1err=sim.x1ERR,c=sim.c,cerr=sim.cERR,mb=sim.mB,mberr=sim.mBERR,
                                           cov_x1_c=sim.COV_x1_c,cov_x1_x0=sim.COV_x1_x0,cov_c_x0=sim.COV_c_x0,
                                           alpha=self.options.alpha,beta=self.options.beta,
                                           x0=sim.x0,sigint=self.options.sigint,z=sim.zHD,M=self.options.simM)
                from astropy.cosmology import Planck13 as cosmo
                if not 'MURES' in sim.__dict__:
                    sim.MURES = sim.MU - cosmo.distmod(sim.zHD).value

        if self.options.scaleb4cuts:
            cols_CC = np.where((sim.SIM_TYPE_INDEX != 1))[0]
            cols_Ia = np.where((sim.SIM_TYPE_INDEX == 1))[0]
            lenCC = float(len(cols_CC))
            lenIa = float(len(cols_Ia))

        sim = self.mkcuts(sim,fitresfile=simfile)
        data = self.mkcuts(data,fitresfile=datafile)

        if self.options.journal:
            mf = factors(len(self.options.histvar))
            if self.options.nplots[0]: ysubplot = self.options.nplots[1]; xsubplot = self.options.nplots[0]
            else:
                ysubplot = mf[len(mf)/2]
                xsubplot = len(self.options.histvar)/ysubplot
            plt.rcParams['figure.figsize'] = (xsubplot*7,ysubplot*7)
            if not self.options.outfile:
                self.options.outfile = 'ovplot_%s.png'%("_".join(self.options.histvar))
        else:
            plt.rcParams['figure.figsize'] = (8.5,11)
            from matplotlib.backends.backend_pdf import PdfPages
            if not self.options.outfile:
                self.options.outfile = 'ovplot_%s.pdf'%("_".join(self.options.histvar))
            if not os.path.exists(self.options.outfile) or self.options.clobber:
                pdf_pages = PdfPages(self.options.outfile)
            else:
                print('File %s exists!  Not clobbering...'%self.options.outfile)
                return(1)

        for histvar,i in zip(self.options.histvar,
                             np.arange(len(self.options.histvar))+1):
            if self.options.journal:
                ax = plt.subplot(ysubplot,xsubplot,i)
                import string
                ax.text(-0.1, 1.05, '%s)'%string.ascii_uppercase[i-1], transform=ax.transAxes, 
                         size=20, weight='bold')
            else:
                if i%3 == 1: fig = plt.figure()
                if i == 3: subnum = 3
                else: subnum = i%3
                if len(self.options.histvar) >= 3:
                    ax = plt.subplot(3,1,subnum)
                else:
                    ax = plt.subplot(len(self.options.histvar),1,subnum)
            if not self.options.journal:
                ax.set_xlabel(histvar,labelpad=0)
            else:
                try:
                    if '$' in histvardict[histvar]:
                        ax.set_xlabel(histvardict[histvar],fontsize=40)
                    else:
                        ax.set_xlabel(histvardict[histvar],fontsize=30)
                except KeyError:
                    ax.set_xlabel(histvar)
                ax.set_ylabel('$N_{SNe}$',labelpad=0,fontsize=30)
                if 'vzHD' in histvar: 
                    ax.set_ylabel(histvardict[histvar],fontsize=30)
                    ax.set_xlabel('$z_{CMB}$',fontsize=30)
                elif 'vmB' in histvar: 
                    ax.set_ylabel(histvardict[histvar],fontsize=30)
                    ax.set_xlabel('$m_{B}$',fontsize=30)

            if 'vzHD' in histvar:
                self.plt2var(data,sim,ax,histvar)
                continue
            if 'vmB' in histvar:
                self.plt2mB(data,sim,ax,histvar)
                continue
                
                
            self.options.histmin,self.options.histmax = None,None
            if len(self.options.cutwin):
                for cutopt in self.options.cutwin:
                    var,min,max = cutopt[0],cutopt[1],cutopt[2]; min,max = float(min),float(max)
                if var == histvar:
                    self.options.histmin = min; self.options.histmax = max
            if not self.options.histmin:
                self.options.histmin = np.min(np.append(sim.__dict__[histvar],data.__dict__[histvar]))
                self.options.histmax = np.max(np.append(sim.__dict__[histvar],data.__dict__[histvar]))


            cols_CC = np.where((sim.SIM_TYPE_INDEX != 1) & 
                               (sim.__dict__[histvar] >= self.options.histmin) &
                               (sim.__dict__[histvar] <= self.options.histmax))[0]
            cols_Ia = np.where((sim.SIM_TYPE_INDEX == 1) & 
                               (sim.__dict__[histvar] >= self.options.histmin) &
                               (sim.__dict__[histvar] <= self.options.histmax))[0]
            if not self.options.scaleb4cuts:
                lenCC = float(len(cols_CC))
                lenIa = float(len(cols_Ia))
            
            # bins command options
            if self.options.bins[0] != None: self.options.nbins = self.options.bins[0]
            if self.options.bins[1] != None: self.options.histmin = self.options.bins[1]
            if self.options.bins[2] != None: self.options.histmax = self.options.bins[2]
            print(histvar,self.options.histmin,self.options.histmax)

            histint = (self.options.histmax - self.options.histmin)/self.options.nbins
            histlen = float(len(np.where((data.__dict__[histvar] > self.options.histmin) &
                                         (data.__dict__[histvar] < self.options.histmax))[0]))
            n_nz = np.histogram(data.__dict__[histvar],bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins))
            
            errl,erru = poisson_interval(n_nz[0])
            ax.plot(n_nz[1][:-1]+(n_nz[1][1]-n_nz[1][0])/2.,n_nz[0],'o',color='k',lw=2,label='data')
            ax.errorbar(n_nz[1][:-1]+(n_nz[1][1]-n_nz[1][0])/2.,n_nz[0],yerr=[n_nz[0]-errl,erru-n_nz[0]],color='k',fmt=' ',lw=2)
            import copy
            n_nz_chi2 = copy.deepcopy(n_nz)
            n_nz = np.histogram(sim.__dict__[histvar],bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins))
            ax.plot((n_nz[1][:-1]+n_nz[1][1:])/2.,n_nz[0]/float(lenIa+lenCC)*histlen,
                    color='k',drawstyle='steps-mid',lw=4,label='All Sim. SNe',ls='--')
            chi2 = np.sum((n_nz[0]/float(lenIa+lenCC)*histlen-n_nz_chi2[0])**2./((erru-errl)/2.)**2.)/float(len(n_nz[0])-1)
            print('chi2 = %.3f for %s'%(chi2,histvar))

            n_nz = np.histogram(sim.__dict__[histvar][cols_CC],bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins))
            ax.plot((n_nz[1][:-1]+n_nz[1][1:])/2.,n_nz[0]/float(lenIa+lenCC)*histlen,
                    color='b',drawstyle='steps-mid',lw=4,label='Sim. CC SNe',ls='-.')
            n_nz = np.histogram(sim.__dict__[histvar][cols_Ia],bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins))
            ax.plot((n_nz[1][:-1]+n_nz[1][1:])/2.,n_nz[0]/float(lenIa+lenCC)*histlen,
                    color='r',drawstyle='steps-mid',lw=2,label='Sim. SNe Ia')

            if self.options.ylim[0] or self.options.ylim[1]: ax.set_ylim([self.options.ylim[0],self.options.ylim[1]])
            if self.options.ylog == 'all' or histvar in self.options.ylog:
                ax.set_yscale('log')
                if not self.options.ylim[0]: ax.set_ylim(bottom=0.5)

            print('Variable: %s'%histvar)
            print('NDATA: %i'%len(data.CID))
            print('MC Scale: %.1f'%(histlen/float(lenIa+lenCC)))
            if lenIa:
                print('N(CC Sim.)/N(Ia Sim.): %.3f'%(lenCC/float(lenIa)))
            else:
                print('N(CC Sim.)/N(Ia Sim.): inf')


            if not self.options.journal and i%3 == 1:
                box = ax.get_position()
                ax.set_position([box.x0, box.y0,# + box.height * 0.15,
                                 box.width, box.height * 0.85])
                # Put a legend below current axis
                ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.6),
                          fancybox=True, ncol=2,numpoints=1)
            if self.options.journal and i == 2:
                ax.legend(loc='upper center', bbox_to_anchor=(0.625, 1.0),
                          fancybox=True,numpoints=1,prop={'size':23})
        

            if self.options.interact:
                plt.show()
            
            if not self.options.journal:
                if i%3 == 1: self.plottitle(ax)
                if not i%3:
                    if not os.path.exists(self.options.outfile) or self.options.clobber:
                        pdf_pages.savefig(fig)

        if not self.options.journal:
            if i%3:
                pdf_pages.savefig(fig)
            pdf_pages.close()
        else:
            if not self.options.outfile:
                outfile = 'ovplot_%s.png'%("_".join(self.options.histvar))
            else: outfile = self.options.outfile
            if not os.path.exists(outfile) or self.options.clobber:
                plt.savefig(outfile)
            else:
                print('File %s exists!  Not clobbering...'%outfile)

    def plt2var(self,data,sim,ax,histvar):
        from scipy.stats import binned_statistic
        plotvar = histvar.split('vzHD')[0]
        histvar = 'zHD'

        self.options.histmin,self.options.histmax = None,None
        if len(self.options.cutwin):
            for cutopt in self.options.cutwin:
                var,min,max = cutopt[0],cutopt[1],cutopt[2]; min,max = float(min),float(max)
            if var == histvar:
                self.options.histmin = min; self.options.histmax = max
        if not self.options.histmin:
            self.options.histmin = np.min(np.append(sim.__dict__[histvar],data.__dict__[histvar]))
            self.options.histmax = np.max(np.append(sim.__dict__[histvar],data.__dict__[histvar]))


        cols_CC = np.where((sim.SIM_TYPE_INDEX != 1) & 
                           (sim.__dict__[histvar] >= self.options.histmin) &
                           (sim.__dict__[histvar] <= self.options.histmax))[0]
        cols_Ia = np.where((sim.SIM_TYPE_INDEX == 1) & 
                           (sim.__dict__[histvar] >= self.options.histmin) &
                           (sim.__dict__[histvar] <= self.options.histmax))[0]
        if not self.options.scaleb4cuts:
            lenCC = float(len(cols_CC))
            lenIa = float(len(cols_Ia))
            
        # bins command options
        if self.options.bins[0]: self.options.nbins = self.options.bins[0]
        if self.options.bins[1]: self.options.histmin = self.options.bins[1]
        if self.options.bins[2]: self.options.histmax = self.options.bins[2]

        histint = (self.options.histmax - self.options.histmin)/self.options.nbins
        histlen = float(len(np.where((data.__dict__[histvar] > self.options.histmin) &
                                     (data.__dict__[histvar] < self.options.histmax))[0]))

        bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins)
        if len(cols_CC):
            binned_cc = binned_statistic(sim.zHD[cols_CC],sim.__dict__[plotvar][cols_CC],bins=bins,statistic='median').statistic
            binned_ccerr = binned_statistic(sim.zHD[cols_CC],sim.__dict__[plotvar][cols_CC],bins=bins,statistic=errfnc).statistic
        binned_ia = binned_statistic(sim.zHD[cols_Ia],sim.__dict__[plotvar][cols_Ia],bins=bins,statistic='median').statistic
        binned_iaerr = binned_statistic(sim.zHD[cols_Ia],sim.__dict__[plotvar][cols_Ia],bins=bins,statistic=errfnc).statistic
        binned_sim = binned_statistic(sim.zHD,sim.__dict__[plotvar],bins=bins,statistic='median').statistic
        binned_simerr = binned_statistic(sim.zHD,sim.__dict__[plotvar],bins=bins,statistic=errfnc).statistic
        binned_data = binned_statistic(data.zHD,data.__dict__[plotvar],bins=bins,statistic='median').statistic
        binned_dataerr = binned_statistic(data.zHD,data.__dict__[plotvar],bins=bins,statistic=errfnc).statistic
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_sim,fmt='--',
                    label='sim spec sample',color='k',capsize=0,lw=4)
        if len(cols_CC):
            ax.errorbar((bins[1:]+bins[:-1])/2.,binned_cc,fmt='-.',
                        label='sim spec sample',color='b',capsize=0,lw=4)
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_ia,
                    label='sim spec sample',color='r',capsize=0,lw=2)
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_data,yerr=binned_dataerr,fmt='o',
                    label='sim spec sample',color='k',capsize=0,lw=2)
        return()

    def plt2mB(self,data,sim,ax,histvar):
        from scipy.stats import binned_statistic
        plotvar = histvar.split('vmB')[0]
        histvar = 'mB'

        self.options.histmin,self.options.histmax = None,None
        if len(self.options.cutwin):
            for cutopt in self.options.cutwin:
                var,min,max = cutopt[0],cutopt[1],cutopt[2]; min,max = float(min),float(max)
            if var == histvar:
                self.options.histmin = min; self.options.histmax = max
        if not self.options.histmin:
            self.options.histmin = np.min(np.append(sim.__dict__[histvar],data.__dict__[histvar]))
            self.options.histmax = np.max(np.append(sim.__dict__[histvar],data.__dict__[histvar]))


        cols_CC = np.where((sim.SIM_TYPE_INDEX != 1) & 
                           (sim.__dict__[histvar] >= self.options.histmin) &
                           (sim.__dict__[histvar] <= self.options.histmax))[0]
        cols_Ia = np.where((sim.SIM_TYPE_INDEX == 1) & 
                           (sim.__dict__[histvar] >= self.options.histmin) &
                           (sim.__dict__[histvar] <= self.options.histmax))[0]
        if not self.options.scaleb4cuts:
            lenCC = float(len(cols_CC))
            lenIa = float(len(cols_Ia))
            
        # bins command options
        if self.options.bins[0]: self.options.nbins = self.options.bins[0]
        if self.options.bins[1]: self.options.histmin = self.options.bins[1]
        if self.options.bins[2]: self.options.histmax = self.options.bins[2]

        histint = (self.options.histmax - self.options.histmin)/self.options.nbins
        histlen = float(len(np.where((data.__dict__[histvar] > self.options.histmin) &
                                     (data.__dict__[histvar] < self.options.histmax))[0]))

        bins=np.linspace(self.options.histmin,self.options.histmax,self.options.nbins)
        binned_cc = binned_statistic(sim.mB[cols_CC],sim.__dict__[plotvar][cols_CC],bins=bins,statistic='median').statistic
        binned_ccerr = binned_statistic(sim.mB[cols_CC],sim.__dict__[plotvar][cols_CC],bins=bins,statistic=errfnc).statistic
        binned_ia = binned_statistic(sim.mB[cols_Ia],sim.__dict__[plotvar][cols_Ia],bins=bins,statistic='median').statistic
        binned_iaerr = binned_statistic(sim.mB[cols_Ia],sim.__dict__[plotvar][cols_Ia],bins=bins,statistic=errfnc).statistic
        binned_sim = binned_statistic(sim.mB,sim.__dict__[plotvar],bins=bins,statistic='median').statistic
        binned_simerr = binned_statistic(sim.mB,sim.__dict__[plotvar],bins=bins,statistic=errfnc).statistic
        binned_data = binned_statistic(data.mB,data.__dict__[plotvar],bins=bins,statistic='median').statistic
        binned_dataerr = binned_statistic(data.mB,data.__dict__[plotvar],bins=bins,statistic=errfnc).statistic
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_sim,fmt='--',
                    label='sim spec sample',color='k',capsize=0,lw=4)
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_cc,fmt='-.',
                    label='sim spec sample',color='b',capsize=0,lw=4)
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_ia,
                    label='sim spec sample',color='r',capsize=0,lw=2)
        ax.errorbar((bins[1:]+bins[:-1])/2.,binned_data,yerr=binned_dataerr,fmt='o',
                    label='sim spec sample',color='k',capsize=0,lw=2)

        return()

        
    def plottitle(self,ax):
        from textwrap import wrap
        import datetime; dt = datetime.date.today()
        titlestr = 'Created %i/%i/%i;  '%(dt.month,dt.day,dt.year)
        if self.options.defaultcuts:
            if self.options.x1cellipse:
                titlestr += '0.3 < c < 0.3; 3.0 < x1 < 3.0; pkMJDERR < 2/(1+z); x1ERR < 1; FITPROB > 0.001; '
            else:
                titlestr += 'c/0.3^2 + x1/3.0^2 < 1; pkMJDERR < 2/(1+z); x1ERR < 1; FITPROB > 0.001; '
        if len(self.options.cutwin):
            for cutopt in self.options.cutwin:
                i,min,max = cutopt[0],cutopt[1],cutopt[2]; min,max = float(min),float(max)
                titlestr += '%s < %s < %s; '%(cutopt[1],cutopt[0],cutopt[2])
        ax.set_title("\n".join(wrap(titlestr[:-2],width=60)),fontsize=12)
        
    def mkcuts(self,fr,fitresfile=None):
        # uncertainties
        sf = -2.5/(fr.x0*np.log(10.0))
        cov_mb_c = fr.COV_c_x0*sf
        cov_mb_x1 = fr.COV_x1_x0*sf
    
        invvars = 1.0 / (fr.mBERR**2.+ self.options.alpha**2. * fr.x1ERR**2. + self.options.beta**2. * fr.cERR**2. + \
                             2.0 * self.options.alpha * (fr.COV_x1_x0*sf) - 2.0 * self.options.beta * (fr.COV_c_x0*sf) - \
                             2.0 * self.options.alpha*self.options.beta * (fr.COV_x1_c) )
        if self.options.defaultcuts:
            if self.options.x1cellipse:
                # I'm just going to assume cmax = abs(cmin) and same for x1
                cols = np.where((fr.x1**2./3.0**2. + fr.c**2./0.3**2. < 1) &
                                (fr.x1ERR < 1) & (fr.PKMJDERR < 2*(1+fr.zHD)) &
                                (fr.FITPROB >= 0.001) & (invvars > 0))
            else:
                cols = np.where((fr.x1 > -3.0) & (fr.x1 < 3.0) &
                                (fr.c > -0.3) & (fr.c < 0.3) &
                                (fr.x1ERR < 1) & (fr.PKMJDERR < 2*(1+fr.zHD)) &
                                (fr.FITPROB >= 0.001) & (invvars > 0))

            for k in fr.__dict__.keys():
                fr.__dict__[k] = fr.__dict__[k][cols]
                
        if len(self.options.cutwin):
            cols = np.arange(len(fr.CID))
            for cutopt in self.options.cutwin:
                i,min,max = cutopt[0],cutopt[1],cutopt[2]; min,max = float(min),float(max)
                if not i in fr.__dict__:
                    if i not in self.options.histvar:
                        print('Warning : key %s not in fitres file %s! Ignoring for this file...'%(i,fitresfile))
                    else:
                        raise RuntimeError('Error : key %s not in fitres file %s!'%(i,fitresfile))
                else:
                    cols = cols[np.where((fr.__dict__[i][cols] >= min) & (fr.__dict__[i][cols] <= max))]

            for k in fr.__dict__.keys():
                fr.__dict__[k] = fr.__dict__[k][cols]
            
        return(fr)

def poisson_interval(k, alpha=0.32): 
    """
    uses chisquared info to get the poisson interval. Uses scipy.stats 
    (imports in function). 
    (http://stackoverflow.com/questions/14813530/poisson-confidence-interval-with-numpy)
    """
    from scipy.stats import chi2
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    low[k == 0] = 0.0
    #if k == 0: 
    #    low = 0.0
    return low, high
    
def factors(n):
    from functools import reduce
    return np.array(reduce(list.__add__,([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))
    
def salt2mu(x1=None,x1err=None,
            c=None,cerr=None,
            mb=None,mberr=None,
            cov_x1_c=None,cov_x1_x0=None,cov_c_x0=None,
            alpha=None,beta=None,
            alphaerr=None,betaerr=None,
            M=19.3,x0=None,sigint=None,
            z=None,peczerr=0.0005):

    sf = -2.5/(x0*np.log(10.0))
    cov_mb_c = cov_c_x0*sf
    cov_mb_x1 = cov_x1_x0*sf
    mu_out = mb + x1*alpha - beta*c - M
    invvars = 1.0 / (mberr**2.+ alpha**2. * x1err**2. + beta**2. * cerr**2. + \
                         2.0 * alpha * (cov_x1_x0*sf) - 2.0 * beta * (cov_c_x0*sf) - \
                         2.0 * alpha*beta * (cov_x1_c) )
        
    zerr = peczerr*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
    muerr_out = np.sqrt(1/invvars + zerr**2. + 0.055**2.*z**2.)
    if sigint: muerr_out = np.sqrt(muerr_out**2. + sigint**2.)
    return(mu_out,muerr_out)
        
if __name__ == "__main__":
    usagestring="""
ovdatamc.py <DataFitresFile> <SimFitresFile>  <varName1:varName2:varName3....>  [--cutwin NN_ITYPE 1 1 --cutwin x1 -3 3]

Given a FITRES file for both data and an SNANA simulation, 
ovdatamc.py creates histograms that compare fit parameters or other 
variables between the two.  If distance moduli/residuals are not 
included in the fitres file, specify MU/MURES as the varName and 
they will be computed with standard SALT2 nuisance parameters.  To 
specify multiple variable names, use colons to separate them.

use -h/--help for full options list
"""

    import os
    import optparse

    hist = ovhist()

    # read in the options from the param file and the command line
    # some convoluted syntax here, making it so param file is not required
    parser = hist.add_options(usage=usagestring)
    options,  args = parser.parse_args()

    try:
        datafile,simfile,histvar = args[0],args[1],args[2]
    except:
        import sys
        print(usagestring)
        print('Error : incorrect or wrong number of arguments')
        sys.exit(1)
        
    import numpy as np
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    
    hist.options = options
    hist.verbose = options.verbose
    hist.clobber = options.clobber
    hist.options.histvar = histvar.split(':')

    hist.main(datafile,simfile)
