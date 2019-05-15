/*******************************************

  TO DO: 
   + May 2016: if input file is ROOT or HBOOK, then 
     output is also ROOT or HBOOK.
  
Created by J. Marriner.
Installed into snana v8_38, January 2010.

Program to take output from the SALT fitter dict files and
1.  Determine alpha and beta parameters
2.  Output a file of mu's for cosmological fits

==USAGE

./SALT2mu <parameter input file>
    file=file.fitres bins=10 zmin=0.02 zmax=1.02 u0=1 u1=1 u2=0 u3=0
    prefix=SALT2mu sigmB=0.0 sigx1=0.0 sigc=.1  p9=.73 

The 'file' argument is required, and it is usually the 
FITRES_DMPFILE argument in the &FITINP namelist that is
input to the snlc_fit.exe program. So if you have the
following namelist for the snlc_fit.exe program,
    &FITINP
       FITRES_DMPFILE = 'mySN.fitres'
    &END
then for the SALT2mu input   file='mySN.fitres'

Additional arguments below are optional:

file = name of fitres file to analyze

nmax=100                 ! fit first 100 events only
nmax=70(SDSS),200(PS1MD) ! fit 70 SDSS and 200 PS1MD
nmax=300,200(PS1MD)      ! fit 300 total, with 200 in PS1MD sub-sample

sigint_fix=0.11            ! fix sigint=0.11 for all data
sigint_fix=0.11,0.09,0.08  ! comma-sep list of sigint_fix for each IDSAMPLE


maxerr_abort_c=0      ! for USESIM_c=T,  increase this to epsilon>0
maxerr_abort_x1=0     ! for USESIM_x1=T, increase this to epsilon>0
maxerr_abort_x0=0

 - - - - -  biasCor options - - - - - 
simfile_biascor=name of simulated fitres file to compute mB bias vs. z
opt_biascor=<option>    ! grep MASK_BIASCOR  SALT2mu.c | grep define
opt_biascor=112    ! as in BBC paper

opt_biascor=368    ! idem, but use SigmaInt(m0,x1,c)=S3x3 from biasCor
                   ! Use S3x3 * SCALE, where SCALE is varied for chi2/dof=1.
                   ! Note that SCALE is sort'of the analog of sigint.

sigint_biascor=<sigmb_bias>   ! instead of auto-computed sigint_biascor
snrmin_sigint_biascor         ! SNRMIN to compute siginit_biascor 
prescale_biascor=<sample>,<preScale>

fieldGroup_biascor='shallow,medium,deep'             ! for biasCor & CCprior
fieldGroup_biascor='C3+X3,X1+E1+S1,C2,X2+E2+S2+C2'

surveygroup_biascor='CFA3+CSP,PS1MD'   ! combine CFA3+CSP into one biasCor
surveygroup_biascor='CFA3+CSP(zbin=.02),PS1MD' 
surveygroup_biascor='CFA3+CSP(zbin=.02:cbin=.04:x1bin=.4),PS1MD' 
surveygroup_biascor='CFA3+CSP(zbin=.02),SDSS(zbin=.04),PS1MD' 

  NOTE: if OPT_PHOTOZ column exists in the input FITRES tables, 
        then each biasCor group  is automatically split into 
        [GROUPNAME]-zSPEC and [GROUPNAME]-zPHOT groups. 
        OPT_PHOTOZ is from &SNCLINP input SNTABLE_LIST='FITRES NOZPHOT'

idsample_select=2+3                ! fit only IDSAMPLE = 2 and 3
surveylist_nobiascor='HST,LOWZ'    ! no biasCor for these surveys

To check sample stats for each surveyGroup and fieldGroup,
   grep IDSAMPLE  <stdout_file>

 - - - - -  CCprior options - - - - -  
simfile_ccprior=name of sim fitres file to compute CC prior and
                    flag to to a BEAMS-like fit
simfile_ccprior=same  --> use same file as simfile_bias
simfile_ccprior=H11   --> no sim; use CC MU-z function from Hlozek 2011
  BEWARE: must use 5D biasCor with simfile_ccprior option;
          1D biasCor + ccprior --> abort.

varname_pIa=name of fitres param containing Prob_Ia

maxprobcc_for_sigint --> compute sigInt from chi2=Ndof for 
                         this Ia-like subset

- - - - - -  binning - - - - - -
nzbin    = number of uniform z bins to use (beware: some bins may be empty)
nlogzbin = number of log-spaced z bins (Jun 21 2018)
zmin  = lower limit on redshift (e.g., 0.02)
zmax  = upper limit on redshift (e.g., 1.00)
powzbin: zbinSize(fit+output) \propto (1+z)^powzbin 
   [constrained by nzbin,zmin,zmax.  Default powzbin=0 --> uniform bins]
powzbin=2.0     --> binSize propto (1+z)^2 for all redshift range.
powzbin=2.0,0.4 --> binSize propto (1+z)^2 for nzbin/2 and z<0.4, then
                    constant binsize for z>0.4. Allows small bins at low-z,
                    without too-big bins at high-z.

min_per_zbin =  min number of SN in z-bin to keep z-bin (default=5)

x1min = lower limit on x1 (-6.0)
x1max = upper limit on x1 (+6.0)
cmin  = lower limit on color (-6.0)
cmax  = upper limit on color (+6.0)

sntype = list of types to select. Examples are
    sntype=120
    sntype=120,105
    sntype=120,105,106  ! comma-separated, no blank spaces
    sntype=-1           ! -> all types used

CUTWIN  <VARNAME>  <MIN> <MAX>
CUTWIN  FITPROB  .01 1.0
CUTWIN  SNRMAX3   8  999999
CUTWIN(NOABORT)  SIM_x1  -2 2    ! do not abort if SIM_x1 isn't there
CUTWIN(DATAONLY) LOGMASS  5 12   ! cut on data only (not on biasCor)

errmask_write=[allowed errcode mask to write output fitres]
errmask_write=-1  -> write everything
errmask_write=0   -> write only selected SN used in fit (default)
errmask_write=64  -> include SNe that fail CUTWIN 
errmask_write=-9  -> write comments only with fit results; no SN lines

sigmB= intrinsic magnitude error (mB)
sigx1= intrinsic stretch error (x1)
sigc= intrinsic color error (x3 a.k.a. c) 
xi01= reduced covariance between mB and x1
xi0c= reduced cov between mB and c
xi1c= reduced cov between x1 and c

ZPOLY_sigmB= [IDSURVEY] [a0] [a1] [a2] [a3]
    -> sigmB = a0 + a1*z + a2*(z^2) + a3*(z^3) for this IDSURVEY
ZPOLY_sigx1=  idem 
ZPOLY_sigc=   idem
ZPOLY_xi01=   idem
ZPOLY_xi0c=   idem
ZPOLY_xi1c=   idem

zVARNAME=zPHOT  ! use photo-z (zPHOT and ZPHOTERR); covariances not used
   or
varname_z=zPHOT

varname_gamma='HOST_LOGMASS' ! default name of variable to fit gamma=HR step
varname_gamma='sSFR'


p1 = alpha0 (0.13)
p2 = beta0  (3.2)
p3 = alpha1 (0.0)  # alpha = alpha0 + z*alpha1
p4 = beta1  (0.0)  # beta  = beta0  + z*beta1

p5 = gamma0 (0.0)        # mag step across hostMass
p6 = gamma1 (0.0)        # z-dependence, gamma = gamma0 + z*gamma1.
p7 = logmass_cen (10.0)  # SN magOff = gamma*(FermiFun-0.5)
                         #  where FermiFun = 1/{ 1 + exp[-(m-m_cen)/tau] }
p8 = logmass_tau (0.01)  #  and m=logmass, m_cen=logmass_cen, tau=logmass_tau

p9  = Omega_L (0.73)
p10 = Omega_k (0.0)
p11 = w (-1.0)
p12 = wa (0.0)
p13 = scalePCC (if using CC prior in BEAMS-like chi2)
p14 = sigint   (if using CC prior in BEAMS) DOES NOT WORK !!!
p15 = alphaHost (dalpha/dlogMhost)
p16 = betaHost  (dbeta/dlogMhost)

# if simfile_biascor=H11, use params 17-22
# see u-options below to pick subset of float params
p17 = H11mucc0
p18 = H11mucc1
p19 = H11mucc2
p20 = H11sigcc0
p21 = H11sigcc1  (added July 9 2018)
p22 = H11sigcc2  (added July 9 2018)

u1 = if true, vary parameter 1 (1=True)
u2 = if true, vary parameter 2 (1=True)
u3 = if true, vary parameter 3 (1=True)
.
.
.
u15=1 --> alphaHost = dalpha/dlogmass
u16=1 --> betaHost  = dbeta /dlogmass
u15=2 --> alphaHost = alpha shift about logmass_cen
u16=2 --> betaHost  = beta  shift about logmass_cen

# if H11 is set and none of the u[nn] are set, then all of them
# will be used as default H11 option
u17=1 --> float H11mucc0
u18=1 --> float H11mucc1
u19=1 --> float H11mucc2
u20=1 --> float H11sigcc0
u21=1 --> float H11sigcc1
u22=1 --> float H11sigcc2

uM0= 0 to fix M0 params to INPUTS.mag0
   = 1 to float fixed M0 in each bin
   = 2 to float M0 as knot with interpolation

fixpar_all=1 --> internally set all float logicals to zero, even if
                 they are set true in the input file
                  (i.e., uM0=0, u1=0, u2=0, etc ...)

blindflag=0  --> turn off blinding
blindflag=1  --> apply sinusoid function to MU-vs-z   
blindflag=2  --> DEFAULT: apply random shift ref OL,w (for data only)
blindflag=66 --> 2+64: apply blindflag=2 option to sim & real data

With default blindflag=2, user can set blinding parameters with
  blindpar9=.1,4522   --> OL -> OL + 0.1*cos(4522)  
  blindpar11=.1,207   --> w  ->  w + 0.1*cos(207) 
  [defaults are 0.06,23434 for OL and 0.2,8432 for w]

h0 = hubble parameter (for cosmological distances)
nommag0 = nominal magnitude 
uave = use average magnitude from cosmology (not nommag0)

zpecerr = error on vpec/c, to REPLACE original vpec/c, not add.
vpecerr = error on peculiar velocity (km/sec), replaces orig vpec.
           zpecerr = vpecerr/c
pecv = LEGACY key for zpecerr [should switch to zpecerr]
if zpecerr==0 then compute zpecerr from biasCor RMS(SIM_VPEC)

lensing_zpar --> add  "z*lensing_zpar" to sigma_int

fitflag_sigmb or sig1fit = 1 -->  find sigmB giving chi2/N = 1
fitflag_sigmb or sig1fit = 2 -->  idem, with extra fit adding 2log(sigma)
redchi2_tol   or sig1tol = tolerance (0.02) on chi2/dof-1

prescale_simdata=<preScale>  # pre scale for sim data
prescale_simcc=<preScale>    # pre-scale only the simulated CC

NSPLITRAN=[NRAN] = number of independent sub-samples to run SALT2mu.
                  A separate fitres file is created for each sub-sample.

iflag_duplicate=1  # 0=ignore, 1=abort, 2=merge

snid_mucovdump='5944'  # after each fit iteration, full muCOV dump 


Default output files (can change names with "prefix" argument)
  SALT2mu.log
  SALT2mu.fitres
  SALT2mu.aux


==HISTORY

  8 April 2009 Original version
  
  Jan 24, 2013 RK
    Fix bug in avemag0_calc(void) ; affected fits in which
    z-bins with too few SN are skipped.


  Feb 21 2013 RK - fix to compile with c++.  Inlcude prototype
                   for each minuit command.

  Apr 27, 2013: 
     in ppar, skip lines with a colon since these keys
     (e.g. ,"OUTDIR:") are used by a master script.

     In chi2 loop (fcn), skip SN that fail cut if errmask_write=0
     --> much faster when cutting lots of SN.

  Apr 30, 2013 RK - allow CUTWIN option to work on x1ERR and cERR;
                    see rawdata_malloc().

  Jun 3, 2013 RK - if z <= 1.0E-8 then skip calculations and set
                   MU-related variables to -999. Avoids NaN/inf
                   when z=0.

  Aug 27, 2013 RK
      - replace fitresCCID_init() with generic fitresVar_init().
        Corresponds with fitres-updates in sntools.c.
        rawdata.sn_name -> rawdata.name[isn], so no more need
        for ptrtok to split the string to extract ccid names.

 Jan 18 2014 RK - replace local NOVAR_ABORT with fitresVar_ABORT()
                  that is part of sntools.o.

 Feb 7 2014: RK - use new feature of fitresVar_init(..) by passing
                  multiple strings; e.g., "z Z zSPEC" Simplifies
                  logic of testing for different redshift names.

 Feb 9 2014: fix subtle z-bin coutning issue in setup_zbins().

 Apr 29 2014: MXSNTYPE -> 20 (was 10)

 May 7 2014: fix command-line override to work with 
                CUTWIN <VARNAME> <min> <max>

 May 27 2014: double fitsc() -> void fitsc() since nothing was returned.

 Sep 14 ,2014: replace one access to SNTOOLS_FITRES with call to
               SNTABLE_GET_VARINFO(ivar, varName, &ICAST).

 Oct 27 2014: use refactored table-read functions
              The reading and fitting works for ascii,hbook,root
              input file formats. However, the output ascii file
              works only if the input is also ascii ... will take
              more work to make the output work for any input
              format.

 Dec 08 2014: new input zVARNAME to specify zVARNAME=zPHOT

 Jan 07, 2015: set sig1flag=0 so that NSPLITRAN>1 works

 Mar 7 2015: refactor sig1 params and flags.
                INPUTS.sig1fit -> fitflag_sigmb
                INPUTS.sig1tol -> redchi2_tol
             New option with fitflag_sigmb=2 --> do additional fit-iter
             with chi2 += 2log(sigma).

 Dec 9  2015: new input key prescale_simcc
 Dec 21 2015: write muerr_nosigmB to output file (muerr without sigint)


 Jan 19 2016: begin major upgrade to incorporate 
      + muBias(Ia) vs. redshift from simulation
      + CC prior from simulation
      + modify CHI2(Ia) -> log(P_BEAMS)
          [similar to Hlokez 2011, but Gaussian CC prior 
           replaced with sim map]
        See new keys simfile_mubias and simfile_ccprior

      + remove global omega_[klm], wde and wa ;
        pass cosPar[4] so that cosmo params can be floated.
        Previously cosmo params were fixed and were not floated
        even if the corresponding u# = 1.

 Feb 17 2016: 
   + NCPAR -> 15 (was 12) and got BEAMS-like formalism sort-of working.
   + add blindflag user-input, and function BLINDED to control output
   + implement CCprior in fcn function; see simdata_CCprior.USE
     See new sim-input keys  
         simfile_bias=FITOPT000_SIM.FITRES
         simfile_ccprior=same
         varname_p1a=NN_FRAC_CELL

  April 8 2016: switch biasCor from 2D to 3D(z,x1,c)

  Apri 16 2016: 
    if MUDIF & MUDIFERR columns are in fitres file, then fit to
    MU = MU(REF) + MUDIF where MU(REF) is computed from -wref
    and -omref. Default wref=-1.0 and omref=0.30 .

    Fix bug writing VARANAMES to output fitres file.
    READTABLE_POINTERS got overwritten reading biasCor file.

 April 26 2016: switch biasCor from 3D to 5D (z,x1,c,alpha,beta).

 May 8 2016: 
    + add wrapper "update_covMatrix" around code to fix cov matrices
      that have bad eigenvalues. Allows using this fix on biasCor sim 
      as well as on data.
    + remove legacy fitres variable SIMZ (check only SIM_ZCMB)
 
           
 May 25 2016: muerr_nosiginit -> muerr_raw (no corrections)
 May 30 2016: read zHD[ERR] instead of z[ERR] (z=zHELIO)

 Jun 07 2016: new input typeIa_ccprior --> P_CC=0 for this SNTYPE.

 Jun 20 2016: use ENVreplace to allow ENV in any input file name.

 Jul 21 2016
   + do biasCor for each IDSAMPLE
   + prefix=NONE by default --> no fitres output
   + errmask_write=-9 --> do NOT write SN lines to output fitres file.

 Aug 7 2016:
   + flag biasCor SNIa with  SIM_NONIA_INDEX=0  instead of SIM_TYPE_INDEX=1.
     The latter (SIM_TYPE_INDEX) can have other values; e.g., =120 for SDSS.

 Aug 22 2016: call new utiloty remove_quote() to remove quotes from
              string inputs.  file='bla.txt' now works the same as
              file=bla.txt

 Aug 23 2016: 
   + major re-factor so that  CELLINFO_XXX -> CELLINFO_XX[IDSAMPLE],
     allowing different biasCor binning for each SAMPLE.
   + allow binning options in () for surveygroup_biascor key.

 Sep 05 2016: add input  host_logmass_split. NOT TESTED YET.

 Sep 9 2016: new input keys:
      zpecerr     vpecerr (km/sec)    lensing_zpar
      Legacy key 'pecv = zpecerr' is still allowed.

 Sep 18 2016:
   + fix zinterp to apply vs. sample
   + fix bug counting chi2/dof.

 Sep 22 2016: add input option   min_per_zbin (default=5)
 Sep 25 2016: fix a few stupid parsing bugs for lensing_zpar and blindflag
 Sep 26 2016: add input snrmin_sigint_biascor (default=60)

 Sep 29 2016: 
   + p5=gamma -> gamma0  and  p6=zeta -> gamma1 ;
     gamma = gamma0 + z*gamma1 = magDiff across host mass
   + p7 -> logmass_cen, p8 -> logmass_tau
   + CUTWIN(NOABORT) BLABLA min max
       --> do not abort if variable BLABLA does not exist.
   + remove input host_logmass_split; it is now fitPar p7 = logmass_cen

 Oct 8 2016:
   + MAXBIN_BIASCOR_1D -> 100,000 (was 50,000) to allow 4x4 axb grid
   + small fix in setup_BININFO_biasCor() to allow for TEXT truncation
     in SALT2 values.

 Nov 18 2016: set a few isnan traps, and set min(biasErr)=1.0E-12
             to avoid crash with WGT=1/ERR^2 = 1/0^2.

 Nov 22 2016: change error to MINOS error (avg of err+ and err-),
              instead of parab errors. With low-stat samples,
              errors seemed too small.

 Nov 26 2016: 
   Replace maxProbcc cut with weighted chi2_1a and weighted nfitsn1a
   so that sigma_int no longer depends on arbitrary Prob_cc cut.

 Dec 2 2016:
   + allow SPLITRAN to work with fitflag_simmb>0
 Dec 11 2016
   + refactor with new function exec_mnpout_mnerrs();
     needed to include SIGINT as part of SPLITRAN summary.

 jan 2 2017:
   to avoid Mac confict with SIGINT, refactor with name changes.
   
 Mar 1 2017:
   write blinding info to M0DIF file; see new func write_blindFlag_message

 Mar 12 2017: 
   + count NPASS vs. stage in biasMapSelect; print stats if no
     events pass cuts.

 April 17 2017: ABORT if command-line item is not recognized.
 April 18 2017: new input dumpflag_nobiascor=1
 April 28 2017: MXa & MXb -> 2 (was 5) 

 May 11 2017: fix few benign bugs found by valgrind.

 Jun 22 2017: 
   + fix NVAR bug in writing M0DIF file
   + skip all events with ERRMASK>0; fixes problem with
     too few events in last z-bin (ERRMASK=4)
   + MAXBIN_z->50 (was 30)
   + MAXBIN_BIASCOR_FITPAR->50 (was 40)
   + optional 2nd argument for powzbin
   + new parameter MAXPAR_MINUIT 80

 Jun 26 2017: fix some awful z-index bugs that showed up when
              there are empty z-bins.

 Jun 27 2017: REFACTOR z-indices
    + all arrays go from 0 to nbin_z-1. No more fortran-like arrays
      to accomodate mnparm starting at index=1.
      ipar=MXCOSPAR is the first M0(z)
    + pass all z-bins to minuit; fix M0 values for empty bins
       instead of suppressing them --> 
       avoid ugly sparse 'izuse' indices.

 Jul 5 2017: move read_SURVEY & get_IDSURVEY into sntools.c[h]

 July 16 2017: 
   + add nmax arg 
   + MINEVT_PER_ZBIN_DEFAULT -> 1 (was 5). See min_per_zbin arg to change.

 Aug 9 2017
   + new blindpar=2 (default) feature to blind cosmology offsets.
     See above description for new keys: 
         blindflag, blindpar9, blindpar11

 Aug 18 2017:
    if ( zpecerr==0 )then read SIM_VPEC from biasCor and compute vpecerr 
    from its RMS. If zpecerr>0 then do NOT read SIM_VPEC ... allows
    back-compatibility and user override.

 Sep 10 2017:
   +  MAXPAR_MINUIT-> 110 (was 80)
   +  fix dimension for BININFO_DEF; use MXz instead of MXpar
 
 Sep 13 2017:
   + fix to work if host mag-step (p5) is fixed to non-zero value,
     but not floated. See INPUTS.USE_GAMMA0

 Sep 25 2017:
   + for zinterp method, reject zDATA outside biasCor range,
     instead of aborting. Print reject stats for eachIDSAMPLE.

 Nov 22 2017: 
   + for fixed cosmo params, speed fitting by replacing  cosmodl 
     calculation with z-binned lookup. See prep_cosmodl_lookup().

  Nov 30 2017: 
    in prep_init(), add sanity checks related to USE_GAMMA0

  Dec 2 2017: works with gzipped FITRES files.

  Dec 21 2017: try to reduce init time with
    + makeSparseList_biasCor
    + replace pow(x,2.0) with x*x 
     
  Jan 8 2018: 
    + read optional VPEC, VPEC_ERR, ZHEL, ZHELERR and SIM_VPEC from 
      data files. If they don't exist, make sure to explicitly set
       values = 0.
    + zpecerr replaces vpec/c in data files by adding 
        [zpecerr^2 - (vpec/c)^2] to zerr^2.
      See new function zerr_adjust().
    + sigma_z -> vpec*(1+z)  [Eq A1 in Davis 2012]
    + for clarity, rename z -> zhd and zerr -> zhderr
    - Beware that D_L still uses incorrect (1+zHD) prefactor.

  Jan 25 2018:
    + remove obsolete vpec_biasCor()
    + new global flag SKIPZCUT_BIASCOR so that z-cut is suppressed
      when selecting SNR>60 events needed to compute sigma_int.

  Jan 26 2018: 
    + lots of refactor to prepare COVMAT[idsample,z,a,b]
    + remove obsolete input  COVMAT_SCALE

  Feb 23 2018:
    + add protection against LOGMASS_TAU=0. See IPAR_LOGMASS_TAU=8

  Mar 01 2018: write M0DIF to fitres output so that one can
               verify  MURES = MU - (MUMODEL+M0DIF)

  Mar 21 2018: 
    + read optional opt_photoz flag; if set, automatically
      split data & biasCor samples into zSPEC and zPHOT.

  Apr 2 2018:
    + add two optional fit parameters:
        p15 = dAlpha/dlog(Mhost)
        p16 = dBetaa/dlog(Mhost)

  May 1 2018:  few changes in zerr_adjust().

  Jun 10 2018
    + new inputs maxerr_abort_[x0,x1,c] = 0 by default.
      For fits with USESIM_[x0,x1,c]=T, the fit-errors are zero;
      setting maxerr_abort_[x] < 0 will allow fit to proceed.

  Jun 19 2018
    + abort if any pIa value is < 0 or > 1
    + new input snid_mucovdump='name' to dump full muCOV calculation.
    + fix some H11 bugs from previous refactor
    + new input option for log-spaced z bins:   nlogzbin=20

  July 5 2018
    + new input sigint_fix (see usage above) to specify sigma_int
      for each IDSAMPLE.

  July 9 2018
    + for H11 option, sigCC -> 2nd order poly(z) instead of constant.
    + u17-u22 options to pick subset of H11 parameters
    + fix awful bug in prob_CCprior_H11

  July 31 2018:
    + for CC prior, read either obsolete SIM_NONIA_INDEX, 
      or new SIM_TEMPLATE_INDEX.

  Oct 8 2018: 
    + add option (opt_biascor += 512) to compute sigint(biasCor)
      for each IDSAMPLE. SIGINT_ABGRID now has 3rd dimension
      [MXNUM_IDSAMPLE]

  Oct 9 2018:
    + for user sigint_fix option (sigint per IDSAMPLE), iterate with
      COV-scale to ensure chi2/dof=1.       

  Oct 23 2018:
    + new input   varname_gamma=HOST_LOGMASS  controls which variable
      is used to fit gamma = HR step.
    + new option   CUTWIN(DATAONLY) LOGMASS 5 15  will apply cut
      only to data and not to bias Cor sim, nor to CC-prior sim.

  Feb 5 2019:
   + minor refactor of ICC_FITCLASS_XXX: UNKNOWN -> OTHER, CC -> TOT

  Apr 29 2019:
   + for NSPLITRAN option, write summary to [prefix]_summary.out,
     and use KEY format.

  May 02 2019: remove TABLEFILE_CLOSE() calls since READ_EXEC closes.

  May 07 2019:
   + fix to work when all parameters are fixed;
     If uM0=0, then set M0 fit bound to be -30 +- 0.001 so that
     MINUIT still runs as usual. Might fix later so that fit is
     skipped when there are no fit params.
   + new fixpar_all option to use SALT2mu as MU-calculator

******************************************************/

#include <stdio.h>      
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_fit.h>  // Jun 13 2016

#include <sys/types.h>
#include <sys/stat.h>

#include "sntools.h" 
#include "sntools_output.h" 

//#define MAXPAR_MINUIT 80   // limit for MINUIT params (June 24 2017)
#define MAXPAR_MINUIT 110  // limit for MINUIT params (Sep 10 2017)

// Maximum number of bins
// Note that 5D biasCor array would take 30*20*20*5*5 = 300,000
// So instead we define a linear array of 50,000 to take less memory,
// especially if we need to raise some of the MAXBIN_ later.
#define MAXBIN_z              85  // 50->85 on 9/11/2017
#define MAXBIN_BIASCOR_FITPAR 50  // for biasCor map variables x1,c
#define MAXBIN_BIASCOR_ALPHA   2  // for biasCor map
#define MAXBIN_BIASCOR_BETA    2  // for biasCor map
#define MAXBIN_BIASCOR_1D  100000  // max on total 5D bins for {a,b,z,x1,c}
#define MINPERCELL_BIASCOR    10  // min per cell for RMS and mBslope corr

// shorter names for local declarations
#define MXa   MAXBIN_BIASCOR_ALPHA
#define MXb   MAXBIN_BIASCOR_BETA
#define MXz   MAXBIN_z
#define MXpar MAXBIN_BIASCOR_FITPAR 

#define MXFIELD_OVERLAP  20  // max number of overlapping fields for event

// define biasCor mask-options for user input opt_biascor
#define MASK_BIASCOR_1DZAVG   1     // bit0: interp MUBIAS vs. z (1D), WGT=1
#define MASK_BIASCOR_1DZWGT   2     // bit1: idem, but  WGT=1/muerr^2
#define MASK_BIASCOR_1DZ      4     // bit2: OR of above
#define MASK_BIASCOR_1D5DCUT  8     // bit3: 1DZWGT, but apply 5D cut (test) 
#define MASK_BIASCOR_5D      16     // bit4: 5D map of biasCor
#define MASK_BIASCOR_MUCOV   32     // bit5: apply MUCOVSCALE vs. z
#define MASK_BIASCOR_SAMPLE  64     // bit6: biasCor vs. IDSAMPLE
#define MASK_BIASCOR_2x2ab  128     // bit7: force 2x2 ab grid (min,max ab)
#define MASK_BIASCOR_COVINT 256     // bit8: biasCor-covint(3x3) x SCALE
#define MASK_BIASCOR_SIGINT_SAMPLE 512  // bit9: sigint(biasCor) vs. IDSAMPLE

#define MASK_BIASCOR_DEFAULT  MASK_BIASCOR_5D + MASK_BIASCOR_MUCOV + MASK_BIASCOR_SAMPLE

#define USEMASK_BIASCOR_COVFIT 1
#define USEMASK_BIASCOR_COVINT 2
#define USEMASK_BIASCOR_COVTOT 3

#define IFLAG_DUPLICATE_IGNORE 0
#define IFLAG_DUPLICATE_ABORT  1
#define IFLAG_DUPLICATE_AVG    2  // use weighted avg of SALT2 fit par.
#define MXSTORE_DUPLICATE     20  // always abort if more than this many

char PATH_SNDATA_ROOT[MXPATHLEN];

// hard wire bins to store biasCor
//                                    mB    x1     c
double  BIASCOR_MINVAL_LCFIT[3]  = {  5.0, -4.0, -0.30 } ;
double  BIASCOR_MAXVAL_LCFIT[3]  = { 30.0, +4.0, +0.50 } ;
double  BIASCOR_BINSIZE_LCFIT[3] = { 25.0,  0.5,  0.05 } ;
char    BIASCOR_NAME_LCFIT[3][4] = { "mB", "x1", "c" } ;
int     BIASCOR_MIN_PER_CELL     = 3  ; // at least this many per 5D cell
int     BIASCOR_MINSUM           = 10 ; // at least this many summed in 3x3x3
double  BIASCOR_SNRMIN_SIGINT    = 60. ; //compute biasCor sigInt for SNR>xxx


/* Number of "cosmological" and "SN" parameters  */
#define MXCOSPAR 24  // 22->24 (July 9 2018)


/* Maximum number parameters (must be >=MAXZBIN+MXCOSPAR) */
#define MAXPAR  MAXBIN_z+MXCOSPAR  // RK: compute instead of hard-wire at 30

/* Minimum SN in a bin for bin to be used in fitting */
#define MINEVT_PER_ZBIN_DEFAULT 1

/* Minimum eigenvalue for cov matrix */
#define EIGMIN_COV 0.0001

/* Nominal M0 */
#define M0_DEFAULT -30.0   

// options for uM0
#define M0FITFLAG_CONSTANT     0 // global constant M0 (no floated)
#define M0FITFLAG_ZBINS_FLAT   1 // float constant M0 in each zbin
#define M0FITFLAG_ZBINS_INTERP 2 // idem, but interpolate between z-bins

// max number of sntype to select
#define MXSNTYPE 20

// max number of random split jobs (RK - July 10, 2012)
#define MXSPLITRAN 100

// define indices for fit parameters
#define NLCPAR    3
#define INDEX_mB  0
#define INDEX_x1  1
#define INDEX_c   2

#define MXCUTWIN 20 // max number of CUTWIN definitions in input file.

#define MXCHAR_LINE 2000 // max character per line of fitres file
//#define MXCHAR_FILENAME 256 // max char len of file name

// define ERRMASK for each type of error
#define ERRMASK_x1         1    //  fails x1min-x1max cut
#define ERRMASK_c          2    //  fails cmin - cmax cut
#define ERRMASK_MINBIN     4    //  in z-bin with too few to fit
#define ERRMASK_z          8    //  fails zmin - zmax cut
#define ERRMASK_sntype    16    //  fails sntype cut
#define ERRMASK_HOST      32    //  fails required host info
#define ERRMASK_CUTWIN    64    //  fails a CUTWIN cut
#define ERRMASK_BADERR   128    //  a fit param error is <=0 or bad COV
#define ERRMASK_SPLITRAN 256    //  not in this random sample
#define ERRMASK_SIMPS    512    //  rejected SIM event from pre-scale

#define dotDashLine "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" 

#define MAXFITITER 20             // max number of fit iterations
#define FITFLAG_DONE         0    // no more fits
#define FITFLAG_CHI2         1    // repeat fit with nominal chi2
#define FITFLAG_2LOGSIGMA    2    // repeat fit with chi2 + 2log(sigma)

#define ICC_FITCLASS_OTHER    0
#define ICC_FITCLASS_IBC      1  // for CC prior
#define ICC_FITCLASS_II       2
#define ICC_FITCLASS_TOT      3  // all NON1A together
#define MXCC_FITCLASS         4
#define MAXMUBIN 20

#define NCOSPAR 4  // size of cosPar array (OL,Ok,w0,wa)

#define MXNUM_SAMPLE  25  // max number of SURVEY/FIELD samples
#define MXCHAR_SAMPLE 100 // max string length of sample name
#define USERFLAG_SURVEYGROUP_SAMPLE  1  // bookkeeping for biasCor IDSAMPLE
#define USERFLAG_FIELDGROUP_SAMPLE   2
#define AUTOFLAG_SURVEYGROUP_SAMPLE  3  // survey automatically added

// ---------------------
double LOGTEN  ;

int NSET_ERRMASK[10000]; // number set per errmask.

int MAXSN ;
int NJOB_SPLITRAN; // number of random split-jobs


int    NSIMDATA,NSIMCC ;  // used to implement prescale_sim[Data,CC]
double RNORM;
double dTWOPI;
double PIFAC; 

// ----------- TYPEDEF STRUCTURES -------------

typedef struct { double VAL[NLCPAR][NLCPAR]; } COV_DEF ;

typedef struct {
  // parameters used for bias correction
  double z, alpha, beta, FITPAR[NLCPAR];
  int idsample ;  //specifies SURVEY/FIELDGROUP
} BIASCORLIST_DEF ;

typedef struct {
  double VAL[NLCPAR];     // bias value on mB,x1,c
  double ERR[NLCPAR];     // error on above
  double RMS[NLCPAR];     // RMS spread of biasCor

  // bias w.r.t. mB is treated different than the binned corrections:
  // bias += (mB-mBoff) + fitPar*mBslope  [where fitPar=mB,x1,c]
  //  double mBoff[NLCPAR], mBslope[NLCPAR] ;

} FITPARBIAS_DEF ;

typedef struct {
  int    nbin;
  double lo[MXz] ;      // min in each bin
  double hi[MXz] ;      // max in each bin
  double avg[MXz];      // avg used for interpolation
  double binSize ;
  int    n_perbin[MXz] ;
  char   varName[20];
} BININFO_DEF ;

typedef struct  {      
  int ia_min, ia_max, ib_min, ib_max;
  double WGT[MAXBIN_BIASCOR_ALPHA][MAXBIN_BIASCOR_BETA] ;
} INTERPWGT_ALPHABETA ;
 

typedef struct {
  // things fixed at start
  BININFO_DEF  ZBIN ;
  BININFO_DEF  DMUBIN ;

  // things that vary in fit
  double    alpha, beta, M0 ;
  double    cosPar[NCOSPAR] ;

  // PDF along DMU, in each z bin and in each FITCLAS(Ib,II,CC)
  double DMUPDF[MXNUM_SAMPLE][MXCC_FITCLASS][MAXBIN_z][MAXMUBIN]; 
  int    NDMUPDF[MXNUM_SAMPLE][MXCC_FITCLASS][MAXBIN_z][MAXMUBIN]; 

  // MEAN and RMS of DMU in each z-bin, to approximate Gaussian
  double DMUAVG[MXNUM_SAMPLE][MXCC_FITCLASS][MAXBIN_z] ;
  double DMURMS[MXNUM_SAMPLE][MXCC_FITCLASS][MAXBIN_z] ;
  
  // CC params used in Hlozek 2011
  double H11_fitpar[10] ;  
} MUZMAP_DEF ;


// define CELL-INFO for multi-dimensional storage of biasCor.
// Actual biasCor stored elsewhere, since different kinds of
// info are stored.
typedef struct {
  // define 1D info separately for each dimension
  BININFO_DEF BININFO_z ;             // redshift bining
  BININFO_DEF BININFO_LCFIT[NLCPAR] ; // binning for mB, x1, c

  // collapse 5D space (a,b,z,x1,c)  into 1D space
  int     NCELL ;  // number of cells to store info 
  int     MAPCELL[MXa][MXb][MXz][MXpar][MXpar]; // J1D=0 to NCELL-1

  // define arrays which depend on SURVEYGROUP/FIELDGROUP sub-sample
  int     *NperCell ;           // Number of events in each cell, vs. J1D
  double  *AVG_z;               // wgt-avg z in each cell, vs J1D
  double  *AVG_LCFIT[NLCPAR];   // idem for mB,x1,c, vs. J1D
} CELLINFO_DEF ;



typedef struct {
  int  LEN_MALLOC ;   // number of lines in file used for malloc
  int  NROW       ;   // total number of rows read from simFile

  char   **name, **field ;       
  double *z, *c, *x1, *x0, *mB ;
  double *SIM_DLMAG ;  // read from simFile
  int    *SIM_TYPE_INDEX ;     // =1 for SNIa
  int    *SIM_NONIA_INDEX;     // more refined info for each CC
  int    *iz_index ;
  int    *icc_class ;          // CC class index to fit separate scales
  int    *idsurvey, *idsample, *opt_photoz ;

  double *CUTVAL[MXCUTWIN];
  int     DOFLAG_CUTWIN[MXCUTWIN]; // flag to apply each cutwin

  FITPARBIAS_DEF  ***FITPARBIAS_ALPHABETA ;
  
} SIMFILE_INFO_DEF ;


typedef struct   // sndata
{
  char   name[20];
  char   field[20];    // added July 2016
  int    opt_photoz;   // optional, >0 ==> photo-z
  int    snid ;       // integer SNID (for sims)
  double fitpar[3];   // mB, x1, c  [mB = -2.5*log(x0) ]
  double covmat_tot[3][3]; // fit + intrinsic scatter matrix
  double covmat_fit[3][3]; // covmat from LCfit (no intrinsic scatter)
  double zhd, zhderr, zhel, zhelerr, vpec, vpecerr ;
  double host_logmass; // Sep 2016

  double mumodel, M0, mu, muerr, muerr_raw, muerr_last, muerrsq_last ;
  double mures, pull, pIa, sigCC_last, sqsigCC_last ;
  double sim_mu, sim_mb, sims, simc, simz;

  int idsurvey, sntype, sim_nonIa_index; 
  int idsample; // computed from idsurvey & field

  // before fit, store bias at each alpha,beta bin
  FITPARBIAS_DEF 
  FITPARBIAS_ALPHABETA[MAXBIN_BIASCOR_ALPHA][MAXBIN_BIASCOR_BETA] ; 
  
  // before fit, store muCOVscale at each alpha,beta bin
  double MUCOVSCALE_ALPHABETA[MAXBIN_BIASCOR_ALPHA][MAXBIN_BIASCOR_BETA] ; 

  // determined from fit
  double fitParBias[NLCPAR];  // bias on mB,x1,c, interpolated over 5D
  double muBias, muBiasErr ;  // muBias applied to model 
  double muCOVscale ;         // muCOV scale applied to muErr
  double muBias_zinterp ;     // bias from interpolating vs. z

  // misc stuf
  int izbin;    // redshift index (user bins for fit & output)
  int errmask;  // error mask 
  int skipfit;  // flag to skip event during fitting
  int warn;     // 1 --> potential problem with COV matrix
} SNDATA_DEF ;


typedef struct {

  // WARNING: update copy_IDSAMPLE_biasCor for each new item

  // group of fields: e.g, '82N'
  char NAME_FIELDGROUP[MXCHAR_SAMPLE] ; 

  // group of surveys: e.g.,CSP+CFA
  char NAME_SURVEYGROUP[MXCHAR_SAMPLE]; 

  // string-option for each sample (user-defined binning)
  char STRINGOPT[MXCHAR_SAMPLE];

  int  OPT_PHOTOZ ;
  // xxx   char NAME_zGROUP[MXCHAR_SAMPLE];  // zSPEC or zPHOT or ''

  // Full name of each sample. e.g, 'SDSS(82N)'
  char NAME[MXCHAR_SAMPLE];  // full name of sample
  
  int  NDATA ;   // Number of data events per sample
  int  NBIASCOR ; // idem for biasCor
  int  NCCPRIOR ; // idem for CCprios sample

  // Dec 2017: store sparse list of rows for each IDSAMPLE with cuts
  int  NBIASCOR_CUTS;
  int *IROW_CUTS ;  // points to rawdata biasCor row

  double zMIN ; // zMIN among data in SAMPLE
  double zMAX ; // zMAX among data in SAMPLE

  // define CELLINFO [ranges and Nbin]
  double  RANGE_REDSHIFT[2]; 
  double  BINSIZE_REDSHIFT ;
  double  BINSIZE_FITPAR[NLCPAR] ;

  int DOFLAG_SELECT ;  // 1 --> select this sample for fit (default=all)
  int DOFLAG_BIASCOR ; // 1--> apply biasCor to this sample

  int IFLAG_ORIGIN ; // see USERFLAG_XXX

  // WARNING: update copy_IDSAMPLE_biasCor for each new item

} SAMPLE_INFO_DEF ;

// ------------ GLOBAL STRUCTURES ----------------

SNDATA_DEF   *data ;
CELLINFO_DEF *CELLINFO_BIASCOR ;  
CELLINFO_DEF *CELLINFO_MUCOVSCALE ; 


// misc global info about data, broken down by survey
struct {
  int    NEVT_PER_SURVEY[MXIDSURVEY];
  double zMIN[MXIDSURVEY];
  double zMAX[MXIDSURVEY];
  int    IFILETYPE ; // TEXT, ROOT or HBOOK

  double sigint_fix[MXIDSURVEY]; // see siginit_fix input (7.05.2018)
  double sqsigint_fix[MXIDSURVEY];
  
} SNDATA_INFO ;



// Mar 11, 2011 RK- define data struct for raw data (before cuts)
//                  Each element is allocated in rawdata_malloc().
// July 18 2016 - double -> int for integer quantities

struct rawdata {
  char   **name, **field ; 
  double *zhd,    *x0,    *x1,    *x3,    *zhel,    *vpec;
  double *zhderr, *x0err, *x1err, *x3err, *zhelerr, *vpecerr;
  double *x0x1_c, *x0x3_c, *x1x3_c;

  double *pIa ;  // use with CCprior
  double *simz, *sim_mu;
  double *simx0, *simx1, *simc; 

  int    *idsurvey, *sntype, *opt_photoz, *sim_nonIa_index; 
  double *dsurvey, *dtype;   // duplicate arrays needed for cutvar pointer

  double *CUTVAL[MXCUTWIN];  // used only to make cuts.
  int     DOFLAG_CUTWIN[MXCUTWIN]; // flag to apply each cutwin
  int     ICUTWIN_GAMMA;        // index for gamma-variable
} rawdata ;


int  ONE_SAMPLE_BIASCOR ;  // flag to lump all samples into one biasCor
int  NSAMPLE_BIASCOR ;
SAMPLE_INFO_DEF SAMPLE_BIASCOR[MXNUM_SAMPLE];

struct {
  // info translated directory from user input strings
  // surveygroup_biascor & fieldgroup_biascor

  // - - - - - - - - - - - - - - - - -
  // user-defined fieldGroups
  int  NFIELDGROUP_USR ;
  char FIELDGROUP_LIST[MXNUM_SAMPLE][MXCHAR_SAMPLE];
  char FIELDGROUP_OPTLIST[MXNUM_SAMPLE][MXCHAR_SAMPLE]; // e.., 'zbin=.04'

  // user-defined survey groups
  int  NSURVEYGROUP_USR ; // number of user-specified survey groups
  int  NSURVEYGROUP_TOT ; // total number, including auto-generated groups
  char SURVEYGROUP_LIST[MXNUM_SAMPLE][MXCHAR_SAMPLE]; // comma-separated list
  char SURVEYGROUP_OPTLIST[MXNUM_SAMPLE][MXCHAR_SAMPLE]; // e.g., 'zbin=.04'

  // IDSURVEY list per SURVEYGROUP
  int NSURVEY_PER_GROUP[MXNUM_SAMPLE];
  int IDSURVEYGROUP_LIST[MXNUM_SAMPLE][MXNUM_SAMPLE];

  // min & max for alpha & beta (Dec 2016)
  double ALPHA_MIN, ALPHA_MAX, BETA_MIN, BETA_MAX;
} INPUTS_SAMPLE_BIASCOR ;




// Mar 12 2017: track NEVT passing each stage of biasMapSelect;
// if all events fail, can pinpoint cut that fails.
#define NSTAGE_biasMapSelect 17
int  NPASS_biasMapSelect[NSTAGE_biasMapSelect] ;
int  SKIPZCUT_BIASCOR ;  // needed to compute sigma_int
char COMMENT_biasMapSelect[NSTAGE_biasMapSelect][40] = {
  "ALL",
  "IDSAMPLE is defined" , 
  "IDSAMPLE is selected" ,
  "true SNIa", 
  "CUTWIN", 
  "valid alpha grid", 
  "valid beta grid", 
  "x1min/x1max cut",
  "cmin/cmax cut", 
  "z in biasCor GRID",
  "x1 in biasCor GRID",
  "c in biasCor GRID",
  "sane COVx0x1", 
  "sane COVx0c", 
  "sane COVx1c", 
  "mBerr cut", 
  "sim preScale"
} ;


struct simdata_bias {
  int  LEN_MALLOC ;  // total number of lines in file = max for malloc
  int  NROW ;       // total number of rows read from simFile
  int  NUSE ;       // number after cuts
  
  // read from simFile. Use floats to save memory
  char **name, **field ;
  float  *z, *zerr ; 
  float  *vpec, *vpecerr ;
  float  *SNRMAX ;
  float  *FITVAL[NLCPAR], *FITERR[NLCPAR], *SIMVAL[NLCPAR];
  float  *FITVAL_IDEAL[NLCPAR], *x0_IDEAL ;
  float  *x0, *x0err ;
  float  *COV_x0x1, *COV_x0c, *COV_x1c ; // read from input table
  float  *COV_FITVAL[NLCPAR][NLCPAR] ;      // computed from above
  float  *SIM_ALPHA, *SIM_BETA;
  float  *SIM_DLMAG, *SIM_LENSDMU, *SIM_ZCMB, *SIM_VPEC ;
  int    *SIM_NONIA_INDEX ;   // =0 for SNIa, >0 for CC
  int    *idsurvey, *idsample, *opt_photoz ;
  float  *CUTVAL[MXCUTWIN];

  double *muerrsq;
  int     DOFLAG_CUTWIN[MXCUTWIN]; // flag to apply each cutwin

  // alpha,beta,z binning
  short  int *IA, *IB, *IZ, *iz ; // store alpha,beta index for each event
  BININFO_DEF BININFO_SIM_ALPHA ; 
  BININFO_DEF BININFO_SIM_BETA ;

  // bias in 5D cells, for each fitPar (mB,x1,c)
  // 5D cells are defined vs. J1D = 0 to CELLINFO_MUCOVSCALE.NCELL-1
  FITPARBIAS_DEF  **FITPARBIAS ;  // vs. SAMPLE and 5D cell

  // define bias-scale on muCOV
  double **MUCOVSCALE ; 

  // for biasCor sample; needed to determine bias on sigma_mu
  double SIGINT_ABGRID[MXNUM_SAMPLE][MXa][MXb] ; 
  double SIGINT_AVG;            // averaged over a,b bins
  
  // define intrinsic scatter matrix 
  int     NEVT_COVINT[MXNUM_SAMPLE][MXz][MXa][MXb] ;
  COV_DEF COVINT[MXNUM_SAMPLE][MXz][MXa][MXb] ;
  COV_DEF COVINT_AVG[MXNUM_SAMPLE];  // avg over a,b,z (for printing)
  BININFO_DEF zCOVINT;

  // wgted-avg redshift  corresponding to each M0 offset (for M0 outfile)
  double zM0[MAXBIN_z];       

  int DOCOR_1D;
  int DOCOR_5D;
  
} simdata_bias ;


 
struct simdata_ccprior {

  int USE ;      // use sim CC as prior
  int USEH11 ;   // use Hlozek CC prior if simfile_ccprior='H11'

  // quantities read from simFile_CCprior
  SIMFILE_INFO_DEF SIMFILE_INFO_ALL ;
  SIMFILE_INFO_DEF SIMFILE_INFO_CC ;

  // computed quantities
  MUZMAP_DEF  MUZMAP ;
  double PCC_biasCorScale[MXNUM_SAMPLE] ; // scale Prob_CC due to biasCor cut
  
} simdata_ccprior ;


#define ORDER_ZPOLY_COVMAT 3   // 3rd order polynom for intrinsic scatter cov.
#define MXSURVEY_ZPOLY_COVMAT 10 

struct INPUTS {
  char filename[MXCHAR_FILENAME]; // input fitres file with c,x1,x0 ...

  int    nmax_tot ; // Nmax to fit for all
  int    nmax[MXIDSURVEY];   // idem by survey
  char   nmaxString[100];

  double maxerr_abort_c, maxerr_abort_x1, maxerr_abort_x0;

  int    blindFlag  ; // suppress cosmo param printout for data
  double blindPar   ; // hard-wired to make cosine arg

  // optional sim-files with corrections/maps
  char simFile_biasCor[MXCHAR_FILENAME];     // to get SNIa mB bias vs. z
  int  opt_biasCor ;
  int  prescale_biasCor[2] ; // subset use [0] of [1]
  double sigint_biasCor ;     // force sigint instead of autocompute
  double snrmin_sigint_biasCor; // SNRMIN used to determine sigint
  double sigma_cell_biasCor; // to weight events in cell for biasCor value

  char fieldGroup_biasCor[120]; // sub-divide biasCor groups based on field
  int  use_fieldGroup_biasCor;  // logical flag for above

  char surveyGroup_biasCor[200]; // combine surveys into group
  int  use_surveyGroup_biasCor;

  char surveyList_noBiasCor[120]; // list of surveys fit, but skip biasCor
  char idsample_select[40];       // e.g., '0+3'

  // ----------
  char simFile_CCprior[MXCHAR_FILENAME];  // to get CC prior, dMU vs. z
  char varname_pIa[40];
  int  typeIa_ccprior ;       // PCC=0 for this sntype
  int  simtype_Ibc[2], simtype_II[2];   // SIM_TYPE ranges for CCprior
  double maxProbCC_for_sigint;  // max P_CC/ProbIa to sum chi2_1a

  int    fitflag_sigmb ;  // flag to fit for sigMb that gives chi2/dof=1
  double redchi2_tol ;    // tolerance in red chi2 (was sig1tol)
  double sigmb_step1 ;        // size of first sigMB step, OR ...
  double scale_covint_step1;  // size of first scale_covint step
  double covint_param_step1;  // one of the above

  char   sigint_fix[200];

  double prescale_simData ;   // prescale for simData (never real data)
  double prescale_simCC ;     // e.g., 2 --> use only 1/2 of sim CC

  double cmin,  cmax  ;
  double x1min, x1max ;
  double zmin,  zmax  ;

  int    NCUTWIN ;
  int    CUTWIN_RDFLAG[MXCUTWIN] ; // 1=> read, 0=> use existing var
  char   CUTWIN_NAME[MXCUTWIN][20];
  double CUTWIN_RANGE[MXCUTWIN][2];
  int    CUTWIN_ABORTFLAG[MXCUTWIN] ;  // 1=> abort if var does not exist
  int    CUTWIN_DATAONLY[MXCUTWIN] ;   // 1=> cut on data only (not biasCor)

  int  Nsntype ;
  int  sntype[MXSNTYPE]; // list of sntype(s) to select
  char sntypeString[100]; // comma-separated string from input file

  int     nzbin ;    // number of uniform-spaced redshift bins
  int     nlogzbin;  // number of log-spaced redshift bins
  double  powzbin ;  // fit & output zbin ~ (1+z)^powzbin
  double  znhalf ;   // z where half of nzbins are log and half are constant

  BININFO_DEF BININFO_z ; // Aug 20 2016
  int     min_per_zbin ;

  char varname_z[40]; // name of redshift variable (default = 'z Z zSPEC')

  char   PREFIX[100] ;
  double parval[MAXPAR];
  double parstep[MAXPAR];
  double parbndmin[MAXPAR]; // min par bound
  double parbndmax[MAXPAR]; // max par bound
  int    ipar[MAXPAR];      // boolean float flag for each param
  int    izpar[MAXPAR];   // iz index or -9 
  int    fixpar_all ;     // global flag to fix all parameters
	       
  double blindpar[MAXPAR][2]; // blind offset = [0] * cos([1])

  char varname_gamma[40]; // name of variable to fit gamma HR step;
                           // e.g, "HOST_LOGMASS" or "SSFR", etc...
  int USE_GAMMA0 ;   // true if p5 is floated or fixed to non-zero value
  int errmask_write; // mask of errors to write SN to output

  // define fixed intrinsic scatter matrix
  double sigmB, sigx1, sigc; // intrinsic sigma on mB, x1 ,c
  double xi01, xi0c, xi1c;   // Correlation coeffic

  double zpecerr ;     // replace obsolete pecv (Sep 9 2016)
  double lensing_zpar; // describes lensing constribution to sigma_int

  int    uave;       // flag to use avg mag from cosmology (not nommag0)
  int    uM0 ;       // flag to float the M0 (default=true)
  int    uzsim ;     // flat to cheat and use simulated zCMB (default=F)
  double nommag0 ;   // nominal SN magnitude
  double H0 ;

  int    FLOAT_COSPAR ;    // internal: TRUE if any COSPAR is floated
  double COSPAR[NCOSPAR];  // OL, Ok, w0, wa: input to dlcosmo

  int zpolyflag;   // z-dependent cov matrix JLM AUG 15 2012

  int  NDUMPLOG; // number of SN to dump to flog file

  int NSPLITRAN ;  // number of random subsets to split jobs (RK July 2012)

  int iflag_duplicate;
  
  int dumpflag_nobiasCor; // 1 --> dump info for data with no valid biasCor

  char SNID_MUCOVDUMP[40]; // dump MUERR info for this SNID (Jun 2018)

  
} INPUTS ;


// Nov 2017: define COSMODL storage for fast dl lookup
struct {
  int USE;
  int NBZ;
  double ZMIN, ZMAX, ZBIN, *z, *dl ;
} COSMODL_LOOKUP ;


// define a separate input structure for redshift polynomial params 
// for intrinsic scatter matrix
struct INPUTS_ZPOLY_COVMAT {
  int IDSURVEY_LIST[MXSURVEY_ZPOLY_COVMAT]; // sparse list of SURVEY IDs
  int NSURVEY; 
  double sigmB[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
  double sigx1[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
  double sigc[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
  double xi01[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
  double xi0c[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
  double xi1c[MXSURVEY_ZPOLY_COVMAT][ORDER_ZPOLY_COVMAT+1];
} INPUTS_ZPOLY_COVMAT ;



char FITPARNAMES_DEFAULT[MXCOSPAR][20] = {
  "blank" ,
  "alpha0      ", "beta0       ", "alpha1      ","beta1       ",
  "gamma0      ", "gamma1      ", "logmass_cen ","logmass_tau ",
  "Omega_L     ", "Omega_k     ", "w           ","wa          ",
  "scalePCC    ", "sigint      ", "alphaHost   ","betaHost    ",
  "H11mucc0    ", "H11mucc1    ", "H11mucc2    ",
  "H11sigcc0   ", "H11sigcc1   ", "H11sigcc2   ",
  "blank4      "
} ;


// define alpha & beta bound for fit and for biasCor extrapolation
double FITPARBOUND_ALPHA[2] = { 0.02, 0.30 } ;
double FITPARBOUND_BETA[2]  = { 1.00, 6.00 } ;

int IPAR_ALPHA0, IPAR_BETA0, IPAR_GAMMA0, IPAR_GAMMA1;
int IPAR_LOGMASS_CEN, IPAR_LOGMASS_TAU ;
int IPAR_scalePCC, IPAR_H11, NPAR_H11_TOT, NPAR_H11_USER ;
int IPAR_COVINT_PARAM ;  // sigint or SCALE_COVINT(biasCor)
int IPAR_OL, IPAR_Ok, IPAR_w0, IPAR_wa;


// define inputs to fit that are read or calculated
struct {
  int NSNTOT ;  // total number of SN read
  int NSNCUTS ; // number of SN readn & passing cuts
  int NSNERR  ; // number of SN flagged with error

  double COVINT_PARAM_FIX ;  //  sigint OR SCALE_COVINT
  double COVINT_PARAM_LAST; 
  char   COVINT_PARAM_NAME[20];

  double ZPOLY_COVMAT[3][3];   // z-dependent scatter matrix

  int NZBIN_TOT[MAXBIN_z]; // total number before cuts
  int NZBIN_FIT[MAXBIN_z]; // number after cuts

  int NFITPAR_ALL ; // fixed + floated
  int NFITPAR_FLOAT ;
  int NFITPAR_FLOAT_z   ;  // number of floated M0(z) bins

  int ISFLOAT_z[MAXBIN_z] ; // 1->floated, vs. z-index
  int ISFLOAT[MAXPAR];      // 1-> floated in fit; 0->fixed, vs. ipar index
  
} FITINP ; 


// define fit results
struct {
  int NSNFIT ;        // Number of SN used in fit; was nsnacc
  int NSNFIT_SPLITRAN[MXSPLITRAN]; // idem in SPLITRAN sub-samples
  int NDOF ;          // Ndof in fit
  int NCALL_FCN ;     // number of calls to FCN function
  int NSNFIT_TRUECC;   // for sim, number of true CC

  double HRMS ;        // rms residual around Hubble line
  double CHI2SUM_MIN;  // global min chi2
  double CHI2RED_ALL;  // global reduced chi2
  double CHI2SUM_1A;   // chi2 sum for Ia subset
  double CHI2RED_1A;   // reduced chi2 for Ia subset
  double NSNFIT_1A;    // number of fit SN with P(CC)/ProbTOT < epsilon
  double ALPHA, BETA;  

  double AVEMAG0 ; // average mag among z bins
  double SNMAG0 ;  // AVEMAG0, or user-input INPUTS.nommag0 

  double M0DIF[MAXBIN_z]; // M0-M0avg per z-bin
  double M0ERR[MAXBIN_z]; // error on above
  double zM0[MAXBIN_z];   // wgted-average z (not from fit)

  char   PARNAME[MAXPAR][40];
  double PARVAL[MXSPLITRAN+1][MAXPAR];
  double PARERR[MXSPLITRAN+1][MAXPAR];

  // keep track of chi2red and sigMB for each iteration
  int    NFIT_ITER ;    // Number of fit iterations (was sig1repeats)
  double CHI2RED_LIST[MAXFITITER]; // JLM AUG 15 2012
  double SIGINT_LIST[MAXFITITER];   // JLM AUG 15 2012

} FITRESULT ;

int  DOFIT_FLAG;    // non-zero --> do another fit iteration

int ISDATA;                // T if no SIM keys are found
int FOUNDKEY_SNTYPE = 0 ;  // T -> found sntype key in fitres file
int FOUNDKEY_SIM    = 0 ;  // T -> is simulation (formerly 'simok')
int FOUNDKEY_SIMx0  = 0 ;  // T -> is SALT2 (formerly 'simx0')
int FOUNDKEY_SIM_NONIA_INDEX = 0 ;
int FOUNDKEY_OPT_PHOTOZ = 0 ;

int NVAR_ORIG ; // NVAR from original ntuple
int NVAR_NEW  ; // NVAR added from SALTmu
int NVAR_TOT  ; // sum of above
char VARNAMES_NEW[MAXPAR*10][40] ;
char VARNAMES_ORIG[100][40] ;

time_t t_start, t_end_init, t_end_fit, t_read_biasCor[3] ;

// parameters for errmsg utility 
#define SEV_INFO   1  // severity flag => give info     
#define SEV_WARN   2  // severity flag => give warning  
#define SEV_ERROR  3  // severity flag => error         
#define SEV_FATAL  4  // severity flag => abort program 

#define BLINDMASK_MUz     1  // apply sinusoid function to MU-vs-z
#define BLINDMASK_FIXPAR  2  // blind fixed params (OL,w)
#define BLINDMASK_SIM    64  // blind applies to simulation

//Main function definitions
void apply_blindpar(void);
void check_data(void);
void check_duplicate_SNID(void);
void apply_nmax_cuts(void);
void merge_duplicates(int N, int *isnList);
void setup_zbins_fit(void);
void setup_BININFO_powz(void);
void setup_BININFO_logz(void);

void fcn(int* npar, double grad[], double* fval,
	 double xval[], int* iflag, void *);

void fitsc(int n,double *s, double *c);
void parse_parFile(char *parFile );
void override_parFile(int argc, char **argv);
void parse_ZPOLY_COVMAT(char *item);
void load_ZPOLY_COVMAT(int IDSURVEY, double Z ) ;
double sum_ZPOLY_COVMAT(double Z, double *polyPar) ;

void parse_CUTWIN(char *item);
int  apply_CUTWIN(int OPT, int *DOFLAG_CUTWIN, double *CUTVAL_LIST);
int  usesim_CUTWIN(char *varName) ;
int  set_DOFLAG_CUTWIN(int ivar, int icut, int isData );

void parse_sntype(char *item);
void parse_nmax(char *item);
void parse_prescale_biascor(char *item);
void parse_powzbin(char *item) ;
void parse_IDSAMPLE_SELECT(char *item);
void parse_sigint_fix(char *item);
void parse_blindpar(char *item) ;

void  prep_input(void);
void  prep_gamma_input(void) ;
void  prep_cosmodl_lookup(void);
int   ppar(char* string);
int   read_data(char *filename);
void  get_zString(char *str_z, char *str_zerr, char *cast) ;

void   prepare_IDSAMPLE_biasCor(void); 
void   set_FIELDGROUP_biasCor(void);
void   set_SURVEYGROUP_biasCor(void);
void   set_BINSIZE_SAMPLE_biasCor(int IDSAMPLE);

int    get_IDSAMPLE(int IDSURVEY, int OPT_PHOTOZ, 
		    char *FIELDGROUP, char *SURVEYGROUP, int DUMPFLAG );
void   sort_IDSAMPLE_biasCor(void) ;
void   copy_IDSAMPLE_biasCor(SAMPLE_INFO_DEF *S1, SAMPLE_INFO_DEF *S2) ;
void   dump_SAMPLE_INFO(char *what) ;
void   match_fieldGroup(char *SNID, char *FIELD, 
			char *FIELDGROUP, char *STRINGOPT ); 
void   match_surveyGroup(char *SNID, int IDSURVEY, 
			 char *SURVEYGROUP, char *STRINGOPT ) ;

void   prepare_biasCor(char *simFile);
void   prepare_biasCor_zinterp(void) ;
void   make_COV_FITVAL_biasCor(int ievt) ;

void   prepare_CCprior(char *simFile);

void   fcn_AlphaBetaWGT(double alpha, double beta, int DUMPFLAG,
			INTERPWGT_ALPHABETA *INTERPWGT, char *callFun );

void   fcn_AlphaBeta(double *xval, double z, double logmass,
		     double *alpha, double *beta);

void   read_simFile_CCprior(char *simFile, SIMFILE_INFO_DEF *SIMFILE_INFO);
void   store_simFile_CCprior(SIMFILE_INFO_DEF *INFO_ALL, 
			     SIMFILE_INFO_DEF *INFO_CC) ;
int    storeBias_CCprior(int n) ;
void   setup_zbins_CCprior (SIMFILE_INFO_DEF *INFO_CC, BININFO_DEF *ZBIN) ;
void   setup_MUZMAP_CCprior(int IDSAMPLE, SIMFILE_INFO_DEF *INFO_CC, 
			    MUZMAP_DEF *MUZMAP) ;
void   setup_DMUPDF_CCprior(int IDSAMPLE, SIMFILE_INFO_DEF *INFO_CC, 
			    MUZMAP_DEF *MUZMAP );
void  dump_DMUPDF_CCprior(int IDSAMPLE, int IZ, MUZMAP_DEF *MUZMAP) ;

double prob_CCprior_sim(int IDSAMPLE, MUZMAP_DEF *MUZMAP, 
		    double z, double dmu, int DUMPFLAG);
double prob_CCprior_H11(int n, double dmu, double *H11_fitpar, 
			double *sqsigCC, double *sigCC_chi2penalty);

void   load_FITPARBIAS_CCprior(int icc, FITPARBIAS_DEF
			       (*FITPARBIAS_ABGRID)[MAXBIN_BIASCOR_BETA]);

int IBINFUN(double VAL, BININFO_DEF *BIN, int OPT_ABORT, char *MSG_ABORT );
void copy_BININFO(BININFO_DEF *BIN1, BININFO_DEF *BIN2);

void  setup_CELLINFO_biasCor(int IDSAMPLE);
void  setup_BININFO_biasCor(int IDSAMPLE, int ipar_LCFIT, int MAXBIN, 
			    BININFO_DEF *BININFO);
void  get_BININFO_biasCor_alphabeta(char *varName, double *VAL_MIN, 
				    double *VAL_MAX, double *VAL_BIN ) ;
				    
void  makeMap_fitPar_biasCor(int ISAMPLE, int ipar_LCFIT);
void  makeMap_sigmu_biasCor(int ISAMPLE);   
void  vpec_biasCor(void);
void  init_sigInt_biasCor_legacy(int IDSAMPLE) ;

void  makeMap_binavg_biasCor(int ISAMPLE);
void  makeSparseList_biasCor(void);

void  calc_zM0_biasCor(void);
void  calc_zM0_data(void);
double zmu_solve(double mu, double *cosPar) ;
void   test_zmu_solve(void);

int   storeDataBias(int n, int dumpFlag ) ;
int   biasMapSelect(int i) ;
void  dumpStages_biasMapSelect(void);
void  read_simFile_biasCor(char *simFile);
void  set_MAPCELL_biasCor(int IDSAMPLE) ;
void  store_iaib_biasCor(void) ;

void  zero_FITPARBIAS(FITPARBIAS_DEF *FITPARBIAS) ;
int   get_fitParBias(char *CID, BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		     FITPARBIAS_DEF *FITPARBIAS);
int get_muCOVscale(char *cid,  BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		   double *muCOVscale ) ;

void   get_muBias(char *NAME, 
		  BIASCORLIST_DEF *BIASCORLIST,  
		  FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MAXBIN_BIASCOR_BETA],
		  double         (*MUCOVSCALE_ABGRID)[MAXBIN_BIASCOR_BETA],
		  INTERPWGT_ALPHABETA *INTERPWGT,  
		  double *FITPARBIAS_INTERP, 
		  double *muBias, double *muBiasErr, double *muCOVscale ) ;

double get_magoff_host(double *hostPar, int n);

void  set_defaults(void);
void  set_fitPar(int ipar, double val, double step, 
		double bndmin, double bndmax, int use) ;
void  set_ERRMASK(int isn_raw, int isn_data);
void  getComment_ERRMASK(int mask, char *comment) ;
void  countData_per_zbin(void) ;
int   reject_simData(int isn_raw);
void  write_fitres(char filnam[ ]);
void  write_fitres_misc(FILE *fout);
void  write_blindFlag_message(FILE *fout) ;
void  append_fitres(FILE *fp, char *CCID, int  indx);
void  get_CCIDindx(char *CCID, int *indx) ;
double  BLIND_OFFSET(int ipar);

int   ISBLIND_FIXPAR(int ipar);
void  printmsg_fitStart(FILE *fp); 
void  printmsg_repeatFit(char *msg) ;

void  outFile_driver(void);
void  write_M0(char *fileName);
void  write_MUERR_INCLUDE(FILE *fp) ;

int   SPLITRAN_ACCEPT(int isn, int snid);
void  SPLITRAN_errmask(void);
void  SPLITRAN_SUMMARY(void);
void  CPU_SUMMARY(void);

int keep_errcode(int errcode) ;

int     prepNextFit(void);

// JLM AUG 15 2012 for fitflag_sigmb
void    conflict_check(void);
// xxx double  get_sigint_calc(double rchi2resid, double orig_sigint);
double  next_covFitPar(double redchi2, double orig_parval, double parstep);
void    recalc_dataCov(); 

//Utility function definitions
double cosmodl_forFit(double z, double *cosPar);
double cosmodl(double z, double *cosPar);
double inc    (double z, double *cosPar);

void gengauss(double r[2]);
void ludcmp(double* a, const int n, const int ndim, int* indx, 
	    double* d, int* icon);
void lubksb(const double* a, const int n, const int ndim, 
       const int* indx, double* b);
double rombint(double f(double z, double *cosPar),
	       double a, double b, double *cosPar, double tol ) ;

double avemag0_calc(int opt_dump);
void   M0dif_calc(void) ;
double fcn_M0(int n, double *M0LIST );

void   printCOVMAT(int NPAR);
double fcn_muerrsq(char *name,double alpha,double beta,
		   double (*COV)[NLCPAR], double z, double zerr, int dumpFlag);
double fcn_muerrz(int OPT, double z, double zerr ) ;
double zerr_adjust(double z, double zerr, double vpecerr, char *name );
void   test_muerrz(void);


double muerrsq_biasCor(int ievt, int maskCov, int *istat_cov) ;

void   get_COVINT_model(int idsample, double (*COVINT)[NLCPAR] ); // returns COVINT
void   get_COVINT_biasCor(int IDSAMPLE, double z, double alpha, double beta, 
			  double (*COVINT)[NLCPAR] ); // returns COVINT
void   dump_COVINT_biasCor(int idsample, int iz, int ia, int ib) ;
void   write_COVINT_biasCor(void);
void   init_COVINT_biasCor(void);
void   zero_COV(double (*COV)[NLCPAR]);
void   scale_COV(double scale, double (*COV)[NLCPAR]);

double muresid_biasCor(int ievt);
int    J1D_biasCor(int ievt, char *msg);
void   J1D_invert_D(int IDSAMPLE, int J1D, double *a, double *b, 
		    double *z, double *x1, double *c);
void   J1D_invert_I(int IDSAMPLE, int J1D, int *ia, int *ib, 
		    int *iz, int *ix1, int *ic);
void   get_J1DNBR_LIST(int IDSAMPLE, int J1D, int *NJ1DNBR, int *J1DNBR_LIST) ;
double WGT_biasCor(int opt, int ievt, char *msg);

//Eigenvalue routine
void rs_(void *n, int *nm,double *err,double *w, bool *matz,
           double *z,double *fv1,double *fv2,int *ierr);  

//MINUIT function definitions
void mninit_(const int* ird,const int* iwr,const int* isav);
void mnparm_(const int* nun, const char chnam[], const double* stval,
	     const double* step, const double* bnd1, const double* bnd2,
	     int* ierflg, int nchnam);
void exec_mnparm(void); // local shell function to call mnparm_
void exec_mnpout_mnerrs(void); // loop over params and store VAL and ERR

void mnseti_(const char ctitle[], int nctitle);
/*
void mnexcm_(double fcn, const char chcom[], const double arglis[],
	     const int* narg, int* ierflg, void*, int nchcom);
*/
void mnpout_(const int* num, char chnam[], double* val, double* error,
	     double* bnd1, double* bnd2, int* ivarbl, int nchnam);


typedef void (mfcn)( int* npar, double grad[], double* fval,
	 double xval[], int* iflag, void*);
void mncomd_(mfcn *fcn, char *text, int *icondn, const int *opt, int len);


void mnstat_(double *chi2min, double *fedm, double *errdef,
	     int *npari, int *nparx, int *istat);

void mnerrs_(int *num, double *eplus, double *eminus, 
	     double *eparab, double *globcc);

void mnemat_(double(*emat)[MAXPAR], int *num) ;


int  rawdata_malloc(char *filname);
void rawdata_init(void);

void malloc_simdata_biasCor(int opt, int NROW);
void malloc_simdata_ccprior(int opt, int NROW, SIMFILE_INFO_DEF *SIMFILE_INFO);


// *********************************
int main(int argc,char* argv[ ])
{
  //History:
  // Dec 5 2016: call exec_mnparm, and cleanup local declarations

  // local declarations
  int inf=5, outf=6, savef=7;
  const int null=0;

  int icondn, i, j, len, num, npari, nparx, istat, ndof, idsample ;
  double chi2min, fedm, errdef, corr;
  char text[100], mcom[50];
  char fnam[] = "main" ;

  // ------------------ BEGIN MAIN -----------------

  RNORM  = RAND_MAX;
  dTWOPI = 8.0*atan(1.0);
  PIFAC  = 1.0/sqrt(dTWOPI);
  LOGTEN = log(10.0) ;

  t_start = time(NULL);

  set_defaults();

  if (argc < 2) {
    sprintf(c1err,"Must give param-input file as arguent");
    sprintf(c2err,"to SALT2mu program. See manual.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n\t Start program %s \n\n", argv[0]) ;
  fflush(stdout);

  // read SURVEY.DEF file and store each survey name vs. IDSURVEY
  read_SURVEYDEF(); 

  //parse input parameter file
  parse_parFile( argv[1] );

  // parse command-line args to override input file
  override_parFile(argc, argv);

  // check for conflicts between input variables
  conflict_check();

  // prepare input
  prep_input();


  //  test_zmu_solve();
  //  test_muerrz(); // xxxx

  // --------------------------------------
  //Read input data from SALT2 fit
  read_data(INPUTS.filename);


  // apply parameter blinding (after we know if DATA are real or sim)
  apply_blindpar();

  //Routine to check the input data.  
  //Simulated data can be manipulated as an option in this routine
  check_data();

  // abort if there are any duplicate SNIDs
  check_duplicate_SNID();

  //Set up z bins for data
  if ( INPUTS.nlogzbin > 0 ) 
    { setup_BININFO_logz(); }
  else
    { setup_BININFO_powz(); }


  // prepare mapindex for each IDSURVEY & FIELD --> for biasCor
  prepare_IDSAMPLE_biasCor();

  // read optional sim maps 
  prepare_biasCor(INPUTS.simFile_biasCor);

  prepare_CCprior(INPUTS.simFile_CCprior);

  // check for user-constraint on nmax (July 2017) AFTER prepare_biasCor 
  apply_nmax_cuts();

  /* Solve for minimum of chi-squared */

  t_end_init = time(NULL);

 DOFIT:
  NJOB_SPLITRAN++ ;
  DOFIT_FLAG = FITFLAG_CHI2 ; 
  printmsg_fitStart(stdout);
  
  SPLITRAN_errmask(); // check for random sub-samples

  setup_zbins_fit();

  FITRESULT.NCALL_FCN = 0 ;
  mninit_(&inf,&outf,&savef);

  strcpy(text,"SALT2mu"); len = strlen(text);  mnseti_(text,len);
  fflush(stdout);

  // execuate minuit mnparm_ commands
  exec_mnparm(); 

  FITRESULT.NFIT_ITER = 0 ;
  
  // check reasons to suppress MINUIT screen dump

  int NOMNPRI = 0;
  if ( ISDATA && INPUTS.blindFlag   ) { NOMNPRI = 1; } 
  if ( IGNOREFILE(INPUTS.PREFIX)    ) { NOMNPRI = 1; } 
  if ( INPUTS.errmask_write == -9   ) { NOMNPRI = 1; } 
  NOMNPRI = 1 ; 
  
  // Beginning of DOFIT loop
  while ( DOFIT_FLAG != FITFLAG_DONE  ) {

    if ( NOMNPRI ) 
      { strcpy(mcom,"SET PRI -1"); } // turn off MINUIT printing
    else
      { strcpy(mcom,"SET PRI 0"); }  // MIMUIT printing
    
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    //Miniut MINIMIZE (find minimum chi-squared)
    strcpy(mcom,"SIM 1000");
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    strcpy(mcom,"MINI");
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);
    //Minuit MINOS (compute errors)
    strcpy(mcom,"MINO");
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    //Final call to FCN at minimum of chi-squared
    strcpy(mcom,"CALL FCN 3");
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);   fflush(stdout);

    mnstat_(&chi2min, &fedm, &errdef,&npari,&nparx,&istat);
    ndof = FITRESULT.NSNFIT - npari; 
    FITRESULT.CHI2SUM_MIN = chi2min ;
    FITRESULT.NDOF        = ndof ;
    FITRESULT.CHI2RED_ALL = chi2min/(double)ndof;

    DOFIT_FLAG = prepNextFit();
    FITRESULT.NFIT_ITER++ ; 
  
  }   // End of fitflag_sigmb  loop


  // ===================================================
  //           SUMMARY and WRAP-UP
  // ===================================================

  fflush(stdout);

  printf("\n**********Fit summary**************\n");
  if (istat==0) 
    { printf("Fit status=%i (Fit invalid)\n",istat); }
  if (istat==1) 
    { printf("Fit status=%i (Bad fit. Diagonal errors only.)\n",istat); }
  if (istat==2) 
    { printf("Fit status=%i (Suspect fit.Errors forced positive.)\n",istat); }
  if (istat==3) 
    { printf("Fit status=%i (Good fit/Errors valid)\n",istat); }


  double redChi2 = FITRESULT.CHI2SUM_MIN/(double)FITRESULT.NDOF ;
  printf("-2lnL/dof = %.2f/%i = %6.3f (%i SN) \n",
	 FITRESULT.CHI2SUM_MIN, FITRESULT.NDOF, redChi2, 
	 FITRESULT.NSNFIT );

  fflush(stdout);

  exec_mnpout_mnerrs(); // Dec 12 2016

  //The exact value of the B-band magnitude shouldn't matter 
  //But take the average over bins (weighted by number of SN 
  FITRESULT.AVEMAG0 = avemag0_calc(1);  // call after PARVAL is loaded.

  // compute M0 - AVEMAG0, and apply optional blind offset
  M0dif_calc();

  // print reduced COV matrix
  printCOVMAT(npari);
  
  // ------------------------------------------------
  // check files to write
  outFile_driver();

  //---------
  if ( NJOB_SPLITRAN < INPUTS.NSPLITRAN ) { goto DOFIT ; }

  t_end_fit = time(NULL);

  SPLITRAN_SUMMARY();
  CPU_SUMMARY();

  printf("\n Done. \n"); fflush(stdout);
  
  return(0);

} // end of main


// ********************************************
void exec_mnparm(void) {

  // Created Dec 5 2016
  // Run loops over parameters and executed MINUIT's  mnparm_ function.
  // Also flag which parameters are/aren't fitted.
  // Note that INPUTS.ipar[i] is a global float-flag,
  // while FITINP.ISFLOAT[i] can turn off for z-bins with
  // too few events.
  //
  // Jun 27 2017: REFACTOR z bins
  //
  // May 08 2018: 
  //   + fix to work with uM0=0 by setting M0 bound to -30+_0.001 mag.

  int i, iz, iMN, len, ierflag, icondn, ISFLOAT ;
  int nzbin = INPUTS.nzbin ;
  double M0min, M0max;
  const int null=0;
  char text[100];
  char fnam[] = "exec_mnparm" ;

  // -------------- BEGIN --------------

  //Setup cosmology parameters for Minuit
  for (i=0; i<MXCOSPAR; i++ )    {
    iMN = i + 1;  // MINUIT fortral-line index
    len = strlen(FITPARNAMES_DEFAULT[i]);
    mnparm_(&iMN, FITPARNAMES_DEFAULT[i], 
	    &INPUTS.parval[i], &INPUTS.parstep[i],
	    &INPUTS.parbndmin[i], &INPUTS.parbndmax[i], &ierflag, len);
    sprintf(FITRESULT.PARNAME[i],"%s", FITPARNAMES_DEFAULT[i] );
  }


  //Setup M0(z) paramters for Minuit
  for (iz=0; iz<nzbin; iz++ )    {

    i   = iz + MXCOSPAR ;
    iMN = i + 1 ;

    M0min=-35.0;  M0max=-25.0;
    ISFLOAT = FITINP.ISFLOAT_z[iz] ;

    if ( INPUTS.uM0 == M0FITFLAG_CONSTANT ) 
      { M0min= M0_DEFAULT-0.001; M0max=M0_DEFAULT+.001; }
    // xxx mark delete    ISFLOAT = FITINP.ISFLOAT_z[iz] * INPUTS.uM0 ;
      
    //               val            step   min  max     boolean
    set_fitPar( i, INPUTS.nommag0, 0.1,   M0min,M0max,  ISFLOAT ); 
    INPUTS.izpar[i] = iz ; 
      
    sprintf(text,"m0_%2.2d",iz);
    len = strlen(text);
    
    mnparm_(&iMN,text,&INPUTS.parval[i],&INPUTS.parstep[i],
	    &INPUTS.parbndmin[i], &INPUTS.parbndmax[i], &ierflag,len);
    sprintf(FITRESULT.PARNAME[i],"%s", text );
  }
  
  FITINP.NFITPAR_ALL   = MXCOSPAR + nzbin ;
  FITINP.NFITPAR_FLOAT = 0 ;
  
  //Fix parameters as indicated by ipar array
  for (i=0;i<MAXPAR;++i)  {
    ISFLOAT           = INPUTS.ipar[i] ;
    FITINP.ISFLOAT[i] = ISFLOAT ;

    if ( ISFLOAT ) {
      FITINP.NFITPAR_FLOAT++ ;
    }
    else {
      iMN = i + 1;
      sprintf(text,"FIX %i",iMN);  len = strlen(text);
      mncomd_(fcn, text, &icondn, &null, len);
    }
  }

  fflush(stdout);

  // Jan 2018: after passing param names to MINIUT, remove blank
  //    spaces from COVINT_PARAM for better screen-print
  trim_blank_spaces(FITRESULT.PARNAME[IPAR_COVINT_PARAM]);

  return ;
  
} // end exec_mnparm


// ***********************************************
void exec_mnpout_mnerrs(void) {

  // Created Dec 12, 2016
  // Mostly moving code out of main.
  // Load FITRESULT.PARVAL and FITRESULT.PARERR 
  // for each fitPar and sigint.
  //
  // Jun 27 2017: REFACTOR z bins

  double PARVAL, PARERR, bnd1, bnd2, eplus,eminus,eparab, globcc;
  int    ipar, iMN, iv ;
  int    LEN_VARNAME = 10 ;
  char text[100], format[80], cPARVAL[40] ;
  //  char fnam[] = "exec_mnpout_mnerrs" ;

  // -------------- BEGIN ----------------

  printf("\nFinal parameter values and errors.\n");  fflush(stdout);

  for (ipar=0; ipar<FITINP.NFITPAR_ALL; ipar++ )  {
    iMN     = ipar + 1 ;
    text[0] = 0 ;
    mnpout_(&iMN, text, &PARVAL, &PARERR, &bnd1,&bnd2, &iv, LEN_VARNAME);
    sprintf(cPARVAL,"%7.4f", PARVAL);
    text[LEN_VARNAME] = 0;
          
    if (iv==0) {
      if (ISBLIND_FIXPAR(ipar) ) { sprintf(cPARVAL,"BLINDED"); }
      printf("par %2i      %10s %s         fixed\n",  ipar, text, cPARVAL);
    }
    else  {
      mnerrs_(&iMN, &eplus, &eminus, &eparab,&globcc);
      PARERR = 0.5*( fabs(eplus) + fabs(eminus) ) ;
	  	
      if ( BLIND_OFFSET(ipar) == 0.0   ) {
	strcpy(format,"par %2i (%2i) %10s %11.4e +/- %.4e "
	       "MINOS (+)%7.4e (-)%7.4e Global CC=%6.3f \n");
	printf(format,ipar,iv,text,PARVAL,eparab,eplus,eminus,globcc);
      }
      else {
	printf("par %2d (%2d) %10s   ***** BLINDED ***** \n",
	       ipar, iv, text);
      }
    }

    // fill global arrays for later
    FITRESULT.PARVAL[NJOB_SPLITRAN][ipar] = PARVAL ;
    FITRESULT.PARERR[NJOB_SPLITRAN][ipar] = PARERR ;

    fflush(stdout);
  } // end num loop

  
  // load sigInt separately in IPAR_SIGINT space 
  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_COVINT_PARAM] = FITINP.COVINT_PARAM_FIX;
  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_COVINT_PARAM] = 1.0E-8 ;
  
  return ;

} // end exec_mnpout_mnerrs


// ***********************************************
void setup_BININFO_logz(void) {

  int nzbin   = INPUTS.nlogzbin ;
  double logzmin = log10(INPUTS.zmin);
  double logzmax = log10(INPUTS.zmax);
  double logzbin = (logzmax-logzmin) / (double)nzbin ;

  int nz;
  double logz_lo, logz_hi, zlo, zhi;
  char fnam[] = "setup_BININFO_logz";

  // ----------- BEGIN ----------

  INPUTS.BININFO_z.nbin    = nzbin ;
  sprintf(INPUTS.BININFO_z.varName,"redshift");

  sprintf(BANNER,"Begin %s ", fnam);
  print_banner(BANNER);
  printf("\t nlogzbin = %d -> logz binsize = %.4f \n",
	 nzbin, logzbin );

  for (nz=0; nz < nzbin; nz++)  {    

    logz_lo = logzmin + logzbin * (double)(nz+0) ;
    logz_hi = logzmin + logzbin * (double)(nz+1);

    zlo = pow(10.0,logz_lo);
    zhi = pow(10.0,logz_hi);
    
    INPUTS.BININFO_z.lo[nz]  = zlo ;
    INPUTS.BININFO_z.hi[nz]  = zhi ;
    INPUTS.BININFO_z.avg[nz] = 0.5 * ( zlo + zhi ) ;
  }


  return ;

} // end  setup_BININFO_logz

// ***********************************************
void setup_BININFO_powz(void) {

  // Created Dec 5 2016, R.Kessler
  //
  // Sep 7 2017: for powzbin=1, set p1 = verySmallNumber to avoid
  //             divide-by-zero

  double zlo, zhi, zlast, Z1PMIN, Z1PMAX, Z1PBIN, Z1, ZNORM, p1 ;
  double zmin, zmax, zbin2, arg;
  double powzbin = INPUTS.powzbin ;
  int nzbin = INPUTS.nzbin ;
  int nz, nztmp, nztmp2 ;
  int LHALF   = (INPUTS.znhalf > 0.0) ;
  char fnam[] = "setup_BININFO_powz" ;

  // ----------- BEGIN ------------

  INPUTS.BININFO_z.nbin    = nzbin ;
  sprintf(INPUTS.BININFO_z.varName,"redshift");

  // - - - - - - - - - - - - - - - - - - - - - -
  //
  //           _ zmax
  //          |           dz  
  // ZNORM *  |      -------------  = 1
  //         _|zmin  (1+z)^powzbin 
  //
  //

  zmin = INPUTS.zmin ;
  zmax = INPUTS.zmax ;
  nztmp = nzbin ;
  if ( LHALF > 0.0 ) { 
    zmax   = INPUTS.znhalf ; 
    nztmp  = nzbin/2; 
    nztmp2 = nzbin - nztmp ;
    zbin2  = (INPUTS.zmax-INPUTS.znhalf) / (double)nztmp2 ;
  }

  // - - - - - -
  sprintf(BANNER,"Begin %s ", fnam);
  print_banner(BANNER);
  printf("\t zbinsize ~ (1+z)^%.2f for z < %.3f \n",
	 powzbin, zmax );
  if ( LHALF ) { 
    printf("\t zbinsize = %.3f for z > %.3f \n", zbin2, zmax );
  }
  // - - - - - -

  p1 = 1.0 - powzbin ;
  if ( p1 == 0.0 ) { p1 = 1.0E-12; } // avoid NaN below

  Z1 = 1.0 + zmin ; Z1PMIN = pow(Z1,p1);
  Z1 = 1.0 + zmax ; Z1PMAX = pow(Z1,p1);
  ZNORM  = p1/(Z1PMAX-Z1PMIN);
  Z1PBIN = 1.0/( ZNORM * (double)nztmp );
  INPUTS.BININFO_z.binSize = Z1PBIN ; // not clear if this is useful
  zlast = INPUTS.zmin ;

  for (nz=0; nz < nztmp; nz++)  {    

    zlo = zlast ;
    arg = p1*Z1PBIN + pow((1.0+zlast),p1) ;
    zhi = pow(arg,1.0/p1) - 1.0 ;
    zlast = zhi;
    
    INPUTS.BININFO_z.lo[nz]  = zlo ;
    INPUTS.BININFO_z.hi[nz]  = zhi ;
    INPUTS.BININFO_z.avg[nz] = 0.5 * ( zlo + zhi ) ;
  }

  // check for constant z-bins after log z-bins
  if ( LHALF > 0.0 ) {
    zmin   = zlast ;
    zmax   = INPUTS.zmax ;

    for (nz=nztmp; nz < nzbin; nz++)  {    
      zlo = zlast ;
      zhi = zlo + zbin2 ;
      zlast = zhi;
      INPUTS.BININFO_z.lo[nz]  = zlo ;
      INPUTS.BININFO_z.hi[nz]  = zhi ;
      INPUTS.BININFO_z.avg[nz] = 0.5 * ( zlo + zhi ) ;
    }
  }


  /* xxxx
  for (nz=0; nz < nzbin; nz++)  {    
    printf(" %s  iz=%2d : %.4f - %.4f \n",
	   fnam, nz, INPUTS.BININFO_z.lo[nz], INPUTS.BININFO_z.hi[nz] );
  }
  fflush(stdout); debugexit(fnam);
  xxxxxx */

  return ;

} // end setup_BININFO_powz


// ********************************************
void setup_zbins_fit(void) {

  // Created Dec 2016
  // Called after reading data and setting up biasCor.
  // Setup z-bins for fitting.
  //
  // Jun 27 2017: REFACTOR z bins

  double z;
  int nzbin = INPUTS.nzbin ;
  int n, nz, izbin, NZFLOAT ;
  char fnam[] = "setup_zbins_fit";

  // ----------- BEGIN --------
  
  // init counters
  for (nz=0; nz < nzbin; nz++)  {
    FITINP.NZBIN_TOT[nz] = 0;
    FITINP.NZBIN_FIT[nz] = 0;
  }
  
  //Put data into bins and count bin population
  for (n=0; n<FITINP.NSNCUTS; ++n) {
    z = data[n].zhd ;    
    izbin = IBINFUN(z, &INPUTS.BININFO_z, 0, "");
    data[n].izbin = izbin;
      
    if (izbin<0 || izbin >= nzbin) 
      { data[n].errmask |= ERRMASK_z ; }  // should be redundant

    
    FITINP.NZBIN_TOT[izbin]++ ;
    if ( data[n].errmask != 0 ) { continue ; }
    if ( data[n].skipfit != 0 ) { continue ; } 
    FITINP.NZBIN_FIT[izbin]++ ;

  } // end loop over NSNCUTS
  
  /*  Count the number of bins with enough SN data to fit */
  NZFLOAT=0;

  for (nz=0; nz < nzbin; ++nz) {

    if ( FITINP.NZBIN_FIT[nz] >= INPUTS.min_per_zbin ) { 
      FITINP.ISFLOAT_z[nz] = 1 ;  // logical flag
      NZFLOAT++ ; 
    }
    else {
      FITINP.ISFLOAT_z[nz] = 0 ;    
    }

    printf(" z=%8.5f - %8.5f  NZBIN(TOT,CUTS)=%6d,%6d   ISFLOAT=%i\n",
	   INPUTS.BININFO_z.lo[nz], INPUTS.BININFO_z.hi[nz],
	   FITINP.NZBIN_TOT[nz], FITINP.NZBIN_FIT[nz], 
	   FITINP.ISFLOAT_z[nz]);
  }

  FITINP.NFITPAR_FLOAT_z = NZFLOAT ;
  printf(" --> Use %d of %d z-bins in fit.\n", NZFLOAT, nzbin );

  // Flag SN in z-bins with fewer thann MINBIN;
  // i.e, Flag SN which are not in a zbin that is being used
  for (n=0; n<FITINP.NSNCUTS; ++n)  {
      izbin = data[n].izbin;
      
      if ( FITINP.ISFLOAT_z[izbin] == 0 )  { 
	data[n].errmask += ERRMASK_MINBIN ; 
	data[n].skipfit  = 1;  // Jun 23 2017
      }
  }


  /*  Check that there is at least one bin to fit */
  if ( NZFLOAT <= 0 )     {
    sprintf(c1err,"no z bin with at least %d SN", INPUTS.min_per_zbin);
    sprintf(c2err,"Check input fitres file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  fflush(stdout);
  return ;

} // end of setup_zbins_fit


// ****************************************
void  apply_blindpar(void) {

  // Created Aug  2017
  // Apply optional FIXPAR blinding to fixed parametes OL,w, etc.
  
  int  ipar ;
  double *blindpar ;
  char fnam[] = "apply_blindpar" ;

  // -------------- BEGIN ---------------
  
  if ( (INPUTS.blindFlag & BLINDMASK_FIXPAR)== 0 ) { return ; }
  
  printf("\n");
  for(ipar=0; ipar < MAXPAR; ipar++ ) {
    if ( ISBLIND_FIXPAR(ipar) ) {
      blindpar             = INPUTS.blindpar[ipar];
      INPUTS.parval[ipar] += (blindpar[0] * cos(blindpar[1]) ) ;
      printf("  BLIND FIXPAR: %s += %8.4f * cos(%f) \n",
	     FITPARNAMES_DEFAULT[ipar], blindpar[0], blindpar[1] );
      fflush(stdout);
    }
  }

  return ;

} // end apply_blindpar

// ****************************************
void check_data(void) {

  // Routine to check the data.  
  // Optionally, the data or the error matrices could be modified 
  // Nov 9 2018: replace fortran rs_ with C version of rs

  int dodiag = 0;
  double spr[3][3] = { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };
  double aveg[3][3]= { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };
  double sint[3][3]= { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };

  int nm=3;
  
  // xxx mark delete   int matz=1;
  bool matz = true ;
  int l, m, n, ierr;  
  double err[3][3], eigval[3], eigvec[3][3], fv1[3], fv2[3];
  double serr=1.0;
  double dm, ds, dc;
  int nsnacc=0;

  char fnam[] = "check_data";

  // ------------- BEGIN ------------------

  // Do manipulations only on fake data

  sprintf(BANNER,"Begin %s", fnam);
  print_banner(BANNER);

  //  printf("\n Begin check_data \n"); fflush(stdout);
 
  if (dodiag) {
    printf("   !!! Warning: Off diagonal elements of the "
	   "error input error matrices will be zeroed !!! ");
  }

  /*
  if (FOUNDKEY_SIM ) 
    { printf("\tt !!! Warning -> Fudging errors !!! \n"); }
  */

  for (n=0; n < FITINP.NSNCUTS; ++n) {

     //For simulated data only, compute error matrix of mb,s,c
     //and also the errors expected from measurement+prior
    if ( !data[n].errmask )  {
      ++nsnacc;
      if (FOUNDKEY_SIM) {       
	dm = data[n].fitpar[0] - data[n].sim_mb;
	ds = data[n].fitpar[1] - data[n].sims;
	dc = data[n].fitpar[2] - data[n].simc;
	spr[0][0] += dm*dm;
	spr[0][1] += dm*ds; 
	spr[1][1] += ds*ds; 
	spr[1][2] += ds*dc; 
	spr[2][2] += dc*dc; 
	spr[0][2] += dm*dc; 
      }  // FOUNDKEY_SIM


      aveg[0][0] += data[n].covmat_tot[0][0];
      aveg[0][1] += data[n].covmat_tot[0][1];
      aveg[0][2] += data[n].covmat_tot[0][2];
      aveg[1][1] += data[n].covmat_tot[1][1];
      aveg[1][2] += data[n].covmat_tot[1][2];
      aveg[2][2] += data[n].covmat_tot[2][2];
    }  // errmask = 0


    //Zero off diagonal terms
    //Scale elements [1][1] and [2][2] by 1.0
    if (dodiag)  {      
      data[n].covmat_tot[0][1] = 0.0;
      data[n].covmat_tot[0][2] = 0.0;
      data[n].covmat_tot[1][2] = 0.0;
      data[n].covmat_tot[1][0] = 0.0;
      data[n].covmat_tot[2][0] = 0.0;
      data[n].covmat_tot[2][1] = 0.0;
      data[n].covmat_tot[1][1] *= 1.0;
      data[n].covmat_tot[2][2] *= 1.0;
    }

     //Apply error scale factor
    for (l=0;l<3;++l) {
	for (m=0;m<3;++m)  {
	  data[n].covmat_tot[l][m] = serr*serr*data[n].covmat_tot[l][m];
	  err[l][m] = data[n].covmat_tot[l][m];
	}
    }
    

    // Find eigenvalues and eigenvectors, convention below
    // err[j][i]*eigvec[0][j] = eigval[0]*eigvec[0][i]
    // xxx rs_(&nm,&nm,&err[0][0],eigval,&matz,&eigvec[0][0],fv1,fv2,&ierr);
    ierr = rs(nm,&err[0][0], eigval, &matz, &eigvec[0][0] );


  } // end nsn loop over data


  // ============================================

  //For simulated data only, compute excess error
  if (FOUNDKEY_SIM)  {
    sint[0][0] = spr[0][0] - aveg[0][0];
    sint[0][1] = spr[0][1] - aveg[0][1];
    sint[0][2] = spr[0][2] - aveg[0][2];
    sint[1][1] = spr[1][1] - aveg[1][1];
    sint[1][2] = spr[1][2] - aveg[1][2];
    sint[2][2] = spr[2][2] - aveg[2][2];
    spr[0][0] = sqrt(spr[0][0]/nsnacc);
    spr[1][1] = sqrt(spr[1][1]/nsnacc);
    spr[2][2] = sqrt(spr[2][2]/nsnacc);
    spr[0][1] = spr[0][1]/(nsnacc*spr[0][0]*spr[1][1]);
    spr[0][2] = spr[0][2]/(nsnacc*spr[0][0]*spr[2][2]);
    spr[1][2] = spr[1][2]/(nsnacc*spr[1][1]*spr[2][2]);
    printf("\nAverage error (fit-sim).\n");
    printf("%f %f %f \n",spr[0][0],spr[0][1],spr[0][2]);
    printf("%f %f \n",spr[1][1],spr[1][2]);
    printf("%f \n",spr[2][2]);
  } // FOUNDKEY_SKM

  aveg[0][0] = sqrt(aveg[0][0]/nsnacc);
  aveg[1][1] = sqrt(aveg[1][1]/nsnacc);
  aveg[2][2] = sqrt(aveg[2][2]/nsnacc);
  aveg[0][1] = aveg[0][1]/(nsnacc*aveg[0][0]*aveg[1][1]);
  aveg[0][2] = aveg[0][2]/(nsnacc*aveg[0][0]*aveg[2][2]);
  aveg[1][2] = aveg[1][2]/(nsnacc*aveg[1][1]*aveg[2][2]);
  printf("\nAverage error (meas+prior) \n");
  printf("%f %f %f \n",aveg[0][0],aveg[0][1],aveg[0][2]);
  printf("%f %f \n",aveg[1][1],aveg[1][2]);
  printf("%f \n",aveg[2][2]);
  if (FOUNDKEY_SIM)   {
    if (sint[0][0]>0.0)sint[0][0] = sqrt(sint[0][0]/nsnacc);
    else sint[0][0] = 0.0;
    if (sint[1][1]>0.0)sint[1][1] = sqrt(sint[1][1]/nsnacc);
    else sint[1][1] = 0.0;
    if (sint[2][2]>0.0)sint[2][2] = sqrt(sint[2][2]/nsnacc);
    else sint[2][2] = 0.0;
    if (sint[0][0]>0.0 && sint[1][1]>0.0)
      {
	sint[0][1] = sint[0][1]/(nsnacc*sint[0][0]*sint[1][1]);
      }
    else sint[0][1] = 0.0;
    if (sint[0][0]>0.0 && sint[2][2]>0.0)
      {
	sint[0][2] = sint[0][2]/(nsnacc*sint[0][0]*sint[2][2]);
      }
    else sint[0][2] = 0.0;
    if (sint[1][1]>0.0 && sint[2][2]>0.0)
      {
	sint[1][2] = sint[1][2]/(nsnacc*sint[1][1]*sint[2][2]);
      }
    else sint[1][2] = 0.0;
    printf("\nExcess error. \n");
    printf("%f %f %f \n",sint[0][0],sint[0][1],sint[0][2]);
    printf("%f %f \n",sint[1][1],sint[1][2]);
    printf("%f \n",sint[2][2]);
  }


  fflush(stdout);
  return;

} // end of check_data

// ******************************************
void apply_nmax_cuts(void) {

  // Created July 17 2017
  // Check user nmax constraints and set skipfit
  // tro trim extra events. See nmax input and parse_nmax().

  int  NTOT, NPERSURVEY[MXIDSURVEY]; // before nmax cuts
  int  ntot, npersurvey[MXIDSURVEY]; // after
  int  id, isn, ntot_iter1, nrej_iter1 ;
  int  nexclude, NREMAIN ;
  char fnam[] = "apply_nmax_cuts" ;

  // --------------- BEGIN --------------

  if ( IGNOREFILE(INPUTS.nmaxString) ) { return ; }

  NTOT = ntot = ntot_iter1 = nrej_iter1 = nexclude = 0 ;
  for(id=0; id < MXIDSURVEY; id++ ) 
    {  NPERSURVEY[id] = npersurvey[id] = 0; }

  // ----- check survey-dependent nmax cuts ------------
  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {
    if ( data[isn].skipfit ) { continue ; }
    id = data[isn].idsurvey ;
    if ( NPERSURVEY[id] >= INPUTS.nmax[id] ) { data[isn].skipfit=1;}
    NPERSURVEY[id]++ ;    NTOT++ ;
    if ( data[isn].skipfit == 0 )
      { ntot_iter1++ ; }
    else
      { nrej_iter1++ ; }

    if ( INPUTS.nmax[id] < 999888000 ) { nexclude++ ; }
  }

  nexclude -= nrej_iter1; // number to exclude for next iteration
  

  if( INPUTS.nmax_tot > ntot_iter1 ) { goto SUMMARY ; }

  
  // compute how many events REMAIN after excluding IDSURVEYs
  // that have already been trimmed to requested size.
  NREMAIN = INPUTS.nmax_tot - nexclude ;

  /*
  printf(" xxx ntot_iter1=%d   nrej_iter[1,2]=%d,%d  nex=%d NREMAIN=%d \n",
	 ntot_iter1, nrej_iter1,nrej_iter2, nexclude, NREMAIN ); fflush(stdout);
  */
    
  // Now apply NTOT requirement
  ntot=0;
  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {
    if ( data[isn].skipfit ) { continue ; }
    id = data[isn].idsurvey ;

    // if this SURVEY is already cut, don't cut more
    if ( INPUTS.nmax[id] < 999888000 ) { continue ; }

    ntot++ ;    
    if ( ntot > NREMAIN ) { data[isn].skipfit=1; }
  }

  // ------------------------------------------------
  // now check stats again and write out before->after for each survey
  
  SUMMARY:

  ntot=0;
  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {
    if ( data[isn].skipfit ) { continue ; }
    id = data[isn].idsurvey ;    ntot++ ;    npersurvey[id]++ ;
  }
    
  print_banner("Data Cut Summary for nmax");

  printf("\t %-20s : %5d -> %5d \n", "ALL", NTOT, ntot );
  for(id = 0; id < MXIDSURVEY; id++ ) {
    if ( NPERSURVEY[id] == 0 ) { continue ; }

    printf("\t %-20s : %5d -> %5d \n", 
	   SURVEY_INFO.SURVEYDEF_LIST[id],
	   NPERSURVEY[id], npersurvey[id] );
  }

  return ;

} // end apply_nmax_cuts

// ******************************************
void check_duplicate_SNID(void) {

  // Sep 2016:
  // Use sorted redshift and its error to flag duplicates.
  // CHeck iflag_duplictes for what to do:
  // 0 -> nothing
  // 1 -> abort
  // 2 -> merge fitparams and cov into one LC
  //

  int  MXSTORE = MXSTORE_DUPLICATE ;
  int  isn, nsn, MEMD, MEMI, isort, *isortList, ORDER_SORT   ;
  double *zList ;
  char fnam[] = "check_duplicate_SNID" ;

  // ----------- BEGIN -----------

  sprintf(BANNER,"Begin %s", fnam);
  print_banner(BANNER);


  // for simulation there is only one duplicate option: abort
  if ( FOUNDKEY_SIM )  { INPUTS.iflag_duplicate = IFLAG_DUPLICATE_ABORT ; }

  nsn  = FITINP.NSNCUTS ;
  MEMD = nsn * sizeof(double) + 4 ;
  MEMI = nsn * sizeof(int) + 4 ;
  zList     = (double *)malloc(MEMD); // allocate integer ID list
  isortList = (int    *)malloc(MEMI); // allocate sorted list

  for ( isn=0; isn < nsn; isn++ )  { zList[isn] = data[isn].zhd ; }

  ORDER_SORT = + 1 ; // increasing order
  sortDouble( nsn, zList, ORDER_SORT, isortList ) ;


  int SAME_z, SAME_SNID, NTMP, idup, isort_last, LAST_DUPL=0 ;
  double z, zerr, z_last, zerr_last ;
  char *snid, *snid_last ;
  int  ISORT_DUPL[MXSTORE_DUPLICATE][MXSTORE_DUPLICATE];
  int  NDUPL_LIST[MXSTORE_DUPLICATE]; // how many duplicates per set
  int  NDUPL_SET ; // number of duplicate sets

  isort_last = -9 ;
  z_last = zerr_last = -9.0 ;
  snid_last = fnam; // anything to avoid compile warning

  NDUPL_SET = 0 ;
  for(idup=0; idup < MXSTORE; idup++ ) 
    {  NDUPL_LIST[idup] = 0 ; }

  for ( isn=0; isn < nsn; isn++ ) {
    isort = isortList[isn];
    z     = data[isort].zhd ;
    zerr  = data[isort].zhderr ;
    snid  = data[isort].name ;
    if ( isn==0 ) { goto SET_LAST; }

    SAME_z    = ( z==z_last && zerr==zerr_last );
    SAME_SNID = ( strcmp(snid,snid_last) == 0 );

    if ( SAME_z && SAME_SNID &&  NDUPL_SET < MXSTORE-1 ) {
      if ( LAST_DUPL == 0 ) { 
	NDUPL_SET++ ; 
	NTMP = NDUPL_LIST[NDUPL_SET-1] ;
	ISORT_DUPL[NDUPL_SET-1][NTMP] = isort_last ;
	NDUPL_LIST[NDUPL_SET-1]++ ;
      }
      NTMP = NDUPL_LIST[NDUPL_SET-1] ;
      if( NTMP < MXSTORE )  { ISORT_DUPL[NDUPL_SET-1][NTMP] = isort ; }
      NDUPL_LIST[NDUPL_SET-1]++ ;
      LAST_DUPL = 1 ;
    }
    else {
      LAST_DUPL = 0 ;
    }

  SET_LAST:
    z_last    = z;   zerr_last=zerr;   isort_last=isort ;
    snid_last = data[isort].name ;
  } // end loop over isn


  if ( NDUPL_SET==0 ) { printf("\t No duplicates found. \n"); goto END; }
  
  printf("   Found %d sets of duplicates: \n", NDUPL_SET);

  for(idup=0; idup < NDUPL_SET ; idup++ ) {

    isort = ISORT_DUPL[idup][0] ; // first duplicate has SNID and z
    NTMP = NDUPL_LIST[idup] ;
    snid = data[isort].name ;
    z    = data[isort].zhd ;

    printf("\t -> %2d with SNID=%8.8s at z=%.3f \n", 
	   NTMP, snid, z );	   
  }

  fflush(stdout);

  if ( NDUPL_SET >= MXSTORE ) {
    sprintf(c1err,"NDUPL_SET=%d exceeds bound, MXSTORE_DUPLICATE=%d",
	    NDUPL_SET, MXSTORE );
    sprintf(c2err,"Check duplicates");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

 
  int iflag = INPUTS.iflag_duplicate ;
  printf("  iflag_duplicate = %d --> ", iflag);
  if ( iflag == IFLAG_DUPLICATE_IGNORE ) {
    printf(" do nothing.\n");
  }
  else if ( iflag == IFLAG_DUPLICATE_ABORT ) {
    printf(" ABORT. \n");
    sprintf(c1err,"Duplicates not allowed.");
    sprintf(c2err,"Check input FITRES file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  else if ( iflag == IFLAG_DUPLICATE_AVG ) {
    printf(" merge duplicates.\n");
    for(idup=0; idup < NDUPL_SET ; idup++ ) 
      { merge_duplicates(NDUPL_LIST[idup], ISORT_DUPL[idup] ); }
  }
  else {
    sprintf(c1err,"Invalid iflag_duplicate=%d", INPUTS.iflag_duplicate );
    sprintf(c2err,"grep IFLAG_DUPLICATE  SALT2mu.c");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


 END:
  free(zList);  free(isortList);

  return;

} // end of check_duplicate_SNID


void merge_duplicates(int NDUPL, int *isnList) {

  // First entry in duplicate list has its SALT2 parameters
  // replaced with weighted average among duplicates.
  // COV -> [ sum 1/cov_orig ]^-1
  // 
  // 2nd, 3rd ... duplicates are flagged to SKIP in the fit.
  //
  // SNRMAX, PKMJD,PKMJDERR are NOT Modifed.
  //
  // Beware that output FITRES file does not have the
  // SALT2 parameters changed.


  int i, isn, ipar, ipar2, ISN_SAVE ;
  double fitpar[NLCPAR], fiterr[NLCPAR];
  double COVFIT_INV[NLCPAR][NLCPAR], COVINT[NLCPAR][NLCPAR];
  double wgt, sqerr, wgtsum[NLCPAR], covfit, covtot  ;
  char stringList_fitparOrig[NLCPAR][100], snid[20] ;
  //  char fnam[] = "merge_duplicates" ;

  // ----------- BEGIN ----------------
  
  ISN_SAVE = isnList[0] ;
  sprintf(snid,"%s", data[ISN_SAVE].name);

  printf("# ------ merge %d duplicates for %s ---------- \n", NDUPL, snid);
  fflush(stdout);

  for(ipar=0; ipar < NLCPAR; ipar++ ) {
    fitpar[ipar] = fiterr[ipar] = wgtsum[ipar] = 0.0 ;
    for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) 
      { COVFIT_INV[ipar][ipar2] = 0.0 ;  }

    stringList_fitparOrig[ipar][0] = 0 ;
  }


  for(i=0; i < NDUPL; i++ ) {
    isn = isnList[i] ;

    // keep first duplicate; remove the rest
    if ( i>0 ) { data[isn].skipfit=1; }

    for(ipar=0; ipar < NLCPAR; ipar++ ) {     // mB,x1,c
      sqerr = data[isn].covmat_tot[ipar][ipar] ;
      wgt   = 1.0/sqerr;
      wgtsum[ipar] += wgt ;
      fitpar[ipar] += wgt * data[isn].fitpar[ipar];

      sprintf(stringList_fitparOrig[ipar],"%s %6.3f",
	      stringList_fitparOrig[ipar], data[isn].fitpar[ipar] );

      for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) 
	{ COVFIT_INV[ipar][ipar2] += 1.0/data[isn].covmat_fit[ipar][ipar2] ; }
    } // loop over mB,x1,c
  }

  get_COVINT_model(-1,COVINT);

  for(ipar=0; ipar < NLCPAR; ipar++ ) { 
    fitpar[ipar] /= wgtsum[ipar] ;
    fiterr[ipar]  = sqrt(1.0/wgtsum[ipar]);
    data[ISN_SAVE].fitpar[ipar] = fitpar[ipar] ;
    for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) { 
      covfit = 1.0 / COVFIT_INV[ipar][ipar2] ; 
      covtot = covfit + COVINT[ipar][ipar2] ; 
      data[ISN_SAVE].covmat_fit[ipar][ipar2] = covfit ;
      data[ISN_SAVE].covmat_tot[ipar][ipar2] = covtot ;
    }
  }


  // print old & new SALT2 params
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    printf("  %s-%2.2s(%s) -> %6.3f \n",  
	   snid, BIASCOR_NAME_LCFIT[ipar], 
	   stringList_fitparOrig[ipar], fitpar[ipar] );
  }

  // print new COV
  printf("  %s-COV(no sigInt | +sigInt): \n", snid);
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    printf("\t");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ ) 
      { printf("%10.3le ", data[ISN_SAVE].covmat_fit[ipar][ipar2] ); }
    printf("  |  ");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ ) 
      { printf("%10.3le ", data[ISN_SAVE].covmat_tot[ipar][ipar2] ); }

    printf("\n");
  }

  fflush(stdout);

  return ;

} // end merge_duplicates

// ******************************************
int prepNextFit(void) {

  // Called after fit,  decide if there is another fit iteration, 
  // and if so determine sigint to use. 
  // Function returns 
  //   0 --> no more fits
  //   1 --> do another chi2 fit 
  //   2 --> another fit with chi2 + 2log(sigma) 
  //
  // July 5 2018: stop if input sigint_fix is set; see STOP_COVFIX
  
  double redchi2, covParam ;
  double step1 = INPUTS.covint_param_step1 ;
  int STOP_TOL, STOP_MXFIT, STOP_COV0, STOP_COVFIX, retCode ;
  int NFIT_ITER = FITRESULT.NFIT_ITER ;
  char msg[100];
  char fnam[] = "prepNextFit" ;

  // ----------------- BEGIN -------------

  retCode = FITFLAG_DONE ;
  FITINP.COVINT_PARAM_LAST = FITINP.COVINT_PARAM_FIX ;

  // check reasons to stop all fitting
  if ( INPUTS.fitflag_sigmb == 0 && simdata_ccprior.USE == 0 )  { 
    return(FITFLAG_DONE); 
  }   // fix sigint --> no more fits

  // ----------------------------------

  redchi2 = FITRESULT.CHI2RED_1A ; 

  // check reasons to stop fitting
  STOP_TOL    = ( fabs(redchi2-1.0) < INPUTS.redchi2_tol ) ;
  STOP_MXFIT  = ( NFIT_ITER == MAXFITITER-1 );
  STOP_COV0   = ( NFIT_ITER > 0 && FITINP.COVINT_PARAM_FIX == 0.0)  ;
  
  // for CC prior, require at least 2 iterations
  if ( simdata_ccprior.USE > 0 && NFIT_ITER == 0 ) { STOP_TOL = 0 ; }

  /*
  printf(" xxxx ------------------------------------- \n");
  printf(" xxxx NFIT_ITER=%d: DOFIT_FLAG = %d   chi2red=%.3f   tol=%.3f "
	 "SIGINT=%.3f \n", 
	 NFIT_ITER, DOFIT_FLAG, 
	 redchi2, INPUTS.redchi2_tol, FITINP.COVINT_PARAM_FIX );
	 fflush(stdout); 
  printf(" xxxx STOP(TOL,MXFIT,SIGMB0) = %d %d %d \n",
	 STOP_TOL, STOP_MXFIT, STOP_COV0 );
  fflush(stdout);
  */
  
  if ( STOP_TOL || STOP_MXFIT || STOP_COV0 ) {

    retCode = FITFLAG_DONE ; 
    printf("\t Final %s value = %0.3f  for chi2(Ia)/dof=%.4f\n",
	   FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	   FITINP.COVINT_PARAM_FIX, redchi2 );
  } 

  else if ( INPUTS.fitflag_sigmb == 0 && simdata_ccprior.USE > 0 ) {
    retCode = FITFLAG_DONE ; 
  } 
  else {
    // Try another covParam 
    // On 2nd iteration, use linear fit to estimate next covParam

    covParam = FITINP.COVINT_PARAM_FIX ;
    FITINP.COVINT_PARAM_FIX = next_covFitPar(redchi2,covParam,step1); 

    if ( FITINP.COVINT_PARAM_FIX < 0 ) {
      FITINP.COVINT_PARAM_FIX = 0.0 ;
    } 
    else {
      recalc_dataCov();
      check_data();
    }
    retCode = FITFLAG_CHI2 ;

    sprintf(msg,"last redchi2=%.3f --> next %s=%.3f)", 
	    redchi2, FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	    FITINP.COVINT_PARAM_FIX );
    printmsg_repeatFit(msg);
  }
  
  return(retCode) ;

} // end of prepNextFit


// ******************************************
void printmsg_repeatFit(char *msg) {
  printf("\n");
  printf("# %s\n", dotDashLine);
  printf("\t REPEAT SALT2mu fit: %s\n", msg);
  printf("# %s\n", dotDashLine);
  printf("\n");
  fflush(stdout);
} 

// ******************************************
void fcn(int *npar, double grad[], double *fval, double xval[],
	 int *iflag, void *not)
{
  //c flat=1 read input, flag 2=gradient, flag=3 is final value
  double M0, alpha, beta, alpha0, beta0, da_dz, db_dz, da_dm, db_dm ;
  double scalePCC, scalePCC_fitpar, nsnfit1a, afix, bfix ;
  int i, n, errmask, nsnfit, nsnfit_truecc, idsample ;
  int DOBIASCOR_1D, DOBIASCOR_5D, DUMPFLAG=0, dumpFlag_muerrsq=0 ;
  double chi2sum_tot, delta, sqdelta, hrms_sum;
  double muerrsq, muerrsq_last, muerrsq_raw, muerrsq_tmp, sqsigCC=0.001 ;
  double chi2, chi2_1a, chi2sum_1a, sigCC_chi2penalty=0.0 ;
  double mb, s, c, z, zerr, logmass;
  double dl, mumodel, muBias, muBiasErr, muCOVscale, magoff_host ;
  double muerr, muerr_raw, muerr_last, COVINT[NLCPAR][NLCPAR] ;
  int    ipar,  INTERPFLAG_ab;
  double omega_l, omega_k, wde, wa, cosPar[NCOSPAR] ;
  double *hostPar, gamma0, gamma1, logmass_cen, logmass_tau ;
  double ProbRatio_1a, ProbRatio_CC, betaTmp ;
  char   *name ;

  BIASCORLIST_DEF     BIASCORLIST ;
  INTERPWGT_ALPHABETA INTERPWGT ;
  int NSAMPLE = NSAMPLE_BIASCOR ;
  int IVAR_GAMMA = rawdata.ICUTWIN_GAMMA ;
  char fnam[]= "fcn";

  // --------------- BEGIN -------------

  FITRESULT.NCALL_FCN++ ;

  //Set input cosmology parameters

  alpha0       = xval[1];
  beta0        = xval[2];
  da_dz        = xval[3];
  db_dz        = xval[4];
  gamma0       = xval[5];
  gamma1       = xval[6];
  logmass_cen  = xval[7];
  logmass_tau  = xval[8];
  omega_l      = xval[9];
  omega_k      = xval[10];
  wde          = xval[11];
  wa           = xval[12];
  scalePCC_fitpar  = xval[13] ;

  da_dm        = xval[15];  // added Apr 2 2018
  db_dm        = xval[16];  // idem

  DOBIASCOR_1D = ( INPUTS.opt_biasCor & MASK_BIASCOR_1DZ );
  DOBIASCOR_5D = ( INPUTS.opt_biasCor & MASK_BIASCOR_5D  );
	
  hostPar = &xval[IPAR_GAMMA0];

  /*
  // xxxxxxxxx
  if ( (ncall_fcn > 1000 && ncall_fcn < 1065) || ncall_fcn<10 ) {
    printf(" xxx ncall=%4d: a=%f b=%f S_CC=%f  sigint=%f\n", 
	   ncall_fcn, xval[1], xval[2], xval[13], xval[14]); 
    fflush(stdout);
  }
  // xxxxxxxxx
  */


  // load cosPar array to pass to functions below (RK Jan 2016)
  cosPar[0] = omega_l ;
  cosPar[1] = omega_k ;
  cosPar[2] = wde ;
  cosPar[3] = wa ;

  // -------------------------------
  if ( simdata_ccprior.USE > 0 ) {   

    // Setup info needed to compute CC prior.
    // Do not include z-dependence for each event because it's too slow
    // to compute this inside the data loop.

    // get DMUPDF in each z bin
    
    simdata_ccprior.MUZMAP.alpha = alpha0 ;
    simdata_ccprior.MUZMAP.beta  = beta0 ;
    simdata_ccprior.MUZMAP.M0    = INPUTS.nommag0 ;
    for(i=0; i < NCOSPAR ; i++ ) 
      { simdata_ccprior.MUZMAP.cosPar[i] = cosPar[i] ; }      

    for(idsample=0; idsample < NSAMPLE; idsample++ ) {
      setup_DMUPDF_CCprior(idsample,
			   &simdata_ccprior.SIMFILE_INFO_CC,
			   &simdata_ccprior.MUZMAP );     
    } 

    if ( simdata_ccprior.USEH11 ) {
      for(i=0; i < NPAR_H11_TOT; i++ )
	{ simdata_ccprior.MUZMAP.H11_fitpar[i] = xval[IPAR_H11+i] ; }
    }

    ProbRatio_1a = ProbRatio_CC = 0.0 ;
  }

  // -------------------------------
  // For biasCor, get INTERP weights for alpha and beta grid.
  // Beware that this is not quite right for z-dependent alpha,beta,
  // but it's faster to compute ia,ib here outside the data loop.
  INTERPFLAG_ab = 0;
  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_5D) ||
       (INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT)  ) {
    INTERPFLAG_ab = 1; 
    // check for z-dependent alpha or beta
    if ( INPUTS.ipar[3] || INPUTS.ipar[4] ) { INTERPFLAG_ab = 2; } 

    // check for hostmass-dependent alpha or beta (April 2 2018)
    if ( INPUTS.ipar[15] || INPUTS.ipar[16] ) { INTERPFLAG_ab = 2; } 
  }

  if(INTERPFLAG_ab) { fcn_AlphaBetaWGT(alpha0, beta0, 0, &INTERPWGT, fnam ); }


  // -------------------------------
  chi2sum_tot = chi2sum_1a = 0.0;
  nsnfit      = nsnfit_truecc = 0 ;
  nsnfit1a = 0.0 ;
  hrms_sum  = 0.0 ;

  FITRESULT.NSNFIT = 0;
  FITRESULT.NSNFIT_TRUECC = 0;
  FITRESULT.HRMS   = 0.0;

  for (n=0; n<FITINP.NSNCUTS; ++n)  {

    errmask = data[n].errmask ;

    if ( data[n].skipfit ) { continue ; }

    data[n].mures   = -999. ;
    data[n].pull    = -999. ;
    data[n].mu      = -999. ;
    data[n].muerr   = -999. ;    
    data[n].muerr_raw = -999. ;  // no scale and no sigInt

    z       = data[n].zhd;    if(z<1.0E-8) {continue ;} // Jun 3 2013
    zerr    = data[n].zhderr ;
    mb      = data[n].fitpar[INDEX_mB];
    s       = data[n].fitpar[INDEX_x1];
    c       = data[n].fitpar[INDEX_c];
    //    logmass = data[n].host_logmass ;
    logmass = rawdata.CUTVAL[IVAR_GAMMA][n] ;  // Jan 20, 2019

    name     = data[n].name ;
    idsample = data[n].idsample ;

    fcn_AlphaBeta(xval, z, logmass, &alpha, &beta); // return alpha,beta

    // for z-dependent alpha,beta, interpolate each event
    if (INTERPFLAG_ab==2) { fcn_AlphaBetaWGT(alpha,beta,0,&INTERPWGT,fnam); }

    // get mag offset for this z-bin
    M0    = fcn_M0(n, &xval[MXCOSPAR] );

    // compute distance modulus from cosmology params
    if ( INPUTS.FLOAT_COSPAR ) {
      dl       = cosmodl_forFit(z, cosPar) ;
      mumodel  = 5.0*log10(dl) + 25.0 ;
    }
    else 
      { mumodel = data[n].mumodel ; }


    // compute error-squared on distance mod
    muerrsq  = fcn_muerrsq(name, alpha, beta, data[n].covmat_tot,
			   z, zerr, 0);
       
    // add optional user-input sigma_int (July 5 2018)    
    // xxx    muerrsq += SNDATA_INFO.sqsigint_fix[idsample]; 

    // --------------------------------
    // Compute bias from biasCor sample
    BIASCORLIST.z        = z;
    BIASCORLIST.alpha    = alpha ;
    BIASCORLIST.beta     = beta ;
    BIASCORLIST.idsample = idsample;

    for(ipar=0; ipar<NLCPAR; ipar++ )
      { BIASCORLIST.FITPAR[ipar] = data[n].fitpar[ipar] ; }

    muBias = muBiasErr = 0.0 ;  muCOVscale=1.0 ; 

    if ( DOBIASCOR_5D ) {
      get_muBias(name, &BIASCORLIST,           // (I) misc inputs
		 data[n].FITPARBIAS_ALPHABETA, // (I) bias at each a,b
		 data[n].MUCOVSCALE_ALPHABETA, // (I) muCOVscale at each a,b
		 &INTERPWGT,             // (I) wgt at each a,b grid point
		 data[n].fitParBias,     // (O) interp bias on mB,x1,c
		 &muBias,        // (O) interp bias on mu
		 &muBiasErr,     // (O) stat-error on above
		 &muCOVscale );  // (O) scale bias on muCOV     
    }
    else if ( DOBIASCOR_1D ) {
      muBias     = data[n].muBias_zinterp ; 
    }

      
    magoff_host = get_magoff_host(hostPar,n);

    // load muBias info into globals
    data[n].muBias      = muBias ;
    data[n].muBiasErr   = muBiasErr ;
    data[n].muCOVscale  = muCOVscale ;

    // zero out muBiasErr after storing it, since adding this
    // would contradict the muCOVscale correction.
    muBiasErr = 0.0 ; 
                   
    // apply bias corrections
    muerrsq  *= muCOVscale ;  // error scale 
    muerrsq  += ( muBiasErr*muBiasErr); 
	
    // ------------------------

    delta    = mb  + alpha*s - beta*c;
    //    delta   += s*s*zeta + s*c*eta + c*c*theta; // leftover from J.M. ??
    delta   -= muBias ;     // bias correction (Jun 2016)
    delta   += magoff_host;  // correct SN mag based on host (Sep 2016)

    data[n].mu = delta ;    // store corrected MU without M0

    delta   -= M0 ;       // subtract M0 from data --> add M0 to model
    delta   -= mumodel ;  // residual = mu - mu_model

    /* Note on mu errors
       The expression below is a reasonable approximation to the
       error on distance modulus mu.  It ignores the error in alpha 
       and beta, but that is generally a small effect except when one
       considers the correlated errors on many SN.  To estimate that
       effect, one can do cosmological fits varying alpha and beta 
       in a parametric way..
    */

    muerr = sqrt(muerrsq);

    // store info
    data[n].mures     = delta;
    data[n].pull      = delta/muerr;
    data[n].mumodel   = mumodel ;
    data[n].muerr     = muerr;
    data[n].M0        = M0 ;      // Mar 1 2018 

    sqdelta       = delta*delta ;
    muerrsq_last  = data[n].muerrsq_last ;
    muerr_last    = data[n].muerr_last ;

    if ( errmask == 0  && simdata_ccprior.USE == 0 ) {
      // original SALT2mu chi2 with only spec-confirmed SNIa 
      nsnfit++ ;        nsnfit1a  = (double)nsnfit ;
      hrms_sum       += sqdelta ;    

      chi2_1a       = sqdelta/muerrsq ;
      chi2          = sqdelta/muerrsq ;
      chi2sum_1a   += chi2_1a ;
      chi2sum_tot  += chi2 ;

      // check option to add log(sigma) term for 5D biasCor
      if ( INPUTS.fitflag_sigmb == 2 ) 
      	{ chi2sum_tot  += log(muerrsq/muerrsq_last); }

    } // end of errmask==0 loop


    if ( errmask == 0 && simdata_ccprior.USE > 0 ) {
      // BEAMS-like chi2 = -2ln [ PIa + PCC ]
      double Prob_1a, Prob_CC, Prob_SUM, dPdmu_1a, dPdmu_CC ;
      double PTOTRAW_1a, PTOTRAW_CC, PTOT_1a, PTOT_CC, PSUM ;
      double muerrsq_update, muerr_update ;
      DUMPFLAG = (n == -44);

      nsnfit++ ;
      scalePCC = scalePCC_fitpar ;

      PTOTRAW_1a  = data[n].pIa ;       // NN_PROB_Ia
      PTOTRAW_CC  = 1.0 - PTOTRAW_1a ;  // CC prob 
      PSUM        = PTOTRAW_1a + scalePCC*PTOTRAW_CC ;
      PTOT_CC     = scalePCC * PTOTRAW_CC/PSUM ;
      PTOT_1a     = PTOTRAW_1a/PSUM ;

      if ( data[n].sntype == INPUTS.typeIa_ccprior )
	{ PTOT_1a=1.0;  PTOT_CC=0.0 ;  } // spec-confirmed SNIa
      
      if ( INPUTS.fitflag_sigmb == 2 ) 
	{ muerrsq_update = muerrsq ; }  // current muerr for 5D biasCor
      else
	{ muerrsq_update = muerrsq_last ; } //  fixed muerr for 1D biasCor
      muerr_update   = sqrt(muerrsq_update);

      chi2_1a  = sqdelta/( muerrsq - muBiasErr*muBiasErr)  ;
      dPdmu_1a = ( exp(-0.5*chi2_1a) ) * (PIFAC/muerr_update) ; 
      Prob_1a  = PTOT_1a * dPdmu_1a ;

      // for CC LH, get raw MUERR without intrinsic scatter
      // This affects H11 option only.
      if ( simdata_ccprior.USEH11 ) {
	dPdmu_CC = prob_CCprior_H11(n, delta, &xval[IPAR_H11], 
				    &sqsigCC, &sigCC_chi2penalty );
      }
      else {
	dPdmu_CC = prob_CCprior_sim(idsample, &simdata_ccprior.MUZMAP, 
				    z, delta, DUMPFLAG );
      }

      Prob_CC  = PTOT_CC * dPdmu_CC ;

      Prob_SUM  = Prob_1a + Prob_CC ;

      // xxxxxxxxxxxxxxxxxxx
      if ( DUMPFLAG ) {
	printf(" xxx --------------------------- \n");
	printf(" xxx fcn dump for SN = '%s' \n", name);
	printf(" xxx dPdmu_CC=%le \n", dPdmu_CC);
	printf(" xxx PTOT_CC = scale*PTOTRAW/PSUM  = %le*%le/%le = %le\n",
	       scalePCC, PTOTRAW_CC, PSUM, PTOT_CC );
	printf(" xxx Prob(1a,CC) = %le, %le \n", Prob_1a, Prob_CC);
	debugexit(fnam);
      }
      // xxxxxxxxxxxxxxxxxxx

      
      // sum total chi2 that includes Ia + CC
      if ( Prob_SUM > 0.0 ) {
	ProbRatio_1a = Prob_1a / Prob_SUM ;
	ProbRatio_CC = Prob_CC / Prob_SUM ;
	nsnfit1a   +=  ProbRatio_1a ; 
	chi2sum_1a += (ProbRatio_1a * chi2_1a) ; 

	Prob_SUM    *= (0.15/PIFAC)  ;  
	chi2sum_tot += ( -2.0*log(Prob_SUM) ) ;
	//  chi2sum_tot += sigCC_chi2penalty ; // prevent sigCC<0 (7.17.2018)
      }
      else {
	chi2sum_tot += 1.0E8 ;
      }

    } // end CCprios.USE if block


    // check things on final pass
    if ( errmask == 0 && *iflag==3 ) {
	
      // Jan 26 2018: store raw muerr without intrinsic cov
      muerrsq_raw = fcn_muerrsq(name,alpha,beta, data[n].covmat_fit,
				z,zerr, 0);

      // check user dump with total error (Jun 19 2018)
      dumpFlag_muerrsq = ( strcmp(name,INPUTS.SNID_MUCOVDUMP) == 0 );
      if ( dumpFlag_muerrsq ) {
	muerrsq_tmp = fcn_muerrsq(name,alpha,beta, data[n].covmat_tot,
				  z, zerr, 1);
      }
      
      muerr_raw   = sqrt(muerrsq_raw);
      data[n].muerr_raw = muerr_raw ;
    
      if ( FOUNDKEY_SIM  ) {
	if ( data[n].sim_nonIa_index>0) { nsnfit_truecc++ ; }
      }

      // store reference errors for 1/sigma term
      if ( DOFIT_FLAG == FITFLAG_CHI2 ) {
	data[n].muerr_last   = muerr ; 
	data[n].muerrsq_last = muerrsq ;
	data[n].sigCC_last   = sqrt(sqsigCC) ;
	data[n].sqsigCC_last = sqsigCC ;
      }
    }  // end *iflag==3 
    
    // ----------------------------------

	  
  } // end loop over SN


  // ===============================================
  // ============= WRAP UP =========================
  // ===============================================

  FITRESULT.NSNFIT        = nsnfit ;
  FITRESULT.NSNFIT_TRUECC = nsnfit_truecc ;
  FITRESULT.NSNFIT_SPLITRAN[NJOB_SPLITRAN] = nsnfit ;
  
  if (*iflag==3)  { // done with fit
    FITRESULT.HRMS = sqrt(hrms_sum/nsnfit); 
    FITRESULT.CHI2SUM_1A = chi2sum_1a ;
    FITRESULT.CHI2RED_1A = chi2sum_1a/(double)(nsnfit1a-FITINP.NFITPAR_FLOAT) ;
    FITRESULT.NSNFIT_1A  = nsnfit1a ; 
    FITRESULT.ALPHA      = alpha ;
    FITRESULT.BETA       = beta ;
  }

  *fval = chi2sum_tot;

  return;

}    // end of fcn

// ==============================================
void  fcn_AlphaBeta(double *xval, double z, double logmass, 
		    double *alpha, double *beta) {

  // Apr 2 2018
  // Return alpha and beta for fcn fit function
  //
  // Aug 19 2018: set OPT_LOGMASS_SLOPE if either ipar<=1 (instead of ==1)
  //              Allows fixing dalpha/dmass & dbeta/dmass wihtout floating.
  //
  double a0          = xval[1] ; // alpha0
  double b0          = xval[2] ; // beta0
  double da_dz       = xval[3] ; // dalpha/dz
  double db_dz       = xval[4] ; // dbeta/dz
  double aHost       = xval[15]; // dalpha/dlog(Mhost)
  double bHost       = xval[16]; // dbeta/dlog(Mhost)
  double logmass_cen = xval[7];  // log(Msplit) for gamma0
  double dlogmass    = logmass - logmass_cen ;
  double alpha_local, beta_local ;
  int    OPT_LOGMASS_SLOPE=0, OPT_LOGMASS_SPLIT=0 ;

  alpha_local = a0  + z*da_dz ;
  beta_local  = b0  + z*db_dz ;  

  if ( INPUTS.ipar[15]<=1 || INPUTS.ipar[16]<=1 ) 
    { OPT_LOGMASS_SLOPE = 1; }

  if ( INPUTS.ipar[15]==2 || INPUTS.ipar[16]==2 ) 
    { OPT_LOGMASS_SPLIT = 1; }

  if ( OPT_LOGMASS_SLOPE ) {
    // alphaHost & betaHost are slopes w.r.t. logmass
    alpha_local += (aHost * dlogmass ) ; 
    beta_local  += (bHost * dlogmass ) ;
  }

  if ( OPT_LOGMASS_SPLIT ) {
    // aalphaHost and betaHost are split values around logmass_cen
    if ( dlogmass > 0.0 ) 
      { alpha_local += aHost/2;  beta_local += bHost/2.0; }
    else
      { alpha_local -= aHost/2;  beta_local -= bHost/2.0; }
  }

  *alpha = alpha_local ;
  *beta  = beta_local ;

  return ;

} // end fcn_AlphaBeta

// ================================
double fcn_M0(int n, double *M0LIST) {

  // return model M0 for this data index 'n'
  // and list of M0LIST in each z bin
  // Jun 27 2017: REFACTOR z bins
  // Jan 29 2019: if no iz1 bin, return(M0) instead of retrn(M0bin0)

  int LDMP=0;
  int iz, iz0, iz1, impar, NBINz ;
  double M0, zdata, zbin0, zbin1, zfrac, M0bin0, M0bin1 ;  
  char fnam[] = "fcn_M0";

  // ----------- BEGIN ----------

  M0 = M0_DEFAULT ;
 
  if ( INPUTS.uM0 == M0FITFLAG_CONSTANT ||
       INPUTS.uM0 == M0FITFLAG_ZBINS_FLAT ) {
    iz    = data[n].izbin;  
    M0    = M0LIST[iz] ;
  }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_INTERP ) {
    // linear interp
    if ( simdata_bias.NROW <= 0 ) {
      sprintf(c1err,"Cannot implement option uM0=%d", INPUTS.uM0 );
      sprintf(c2err,"without simFile_bias");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
    }
    
    NBINz   = INPUTS.BININFO_z.nbin ;
    zdata   = data[n].zhd ;
    iz0     = data[n].izbin;  
    zbin0   = simdata_bias.zM0[iz0]; // wgt z-avg in this z bin

    // get neighbor-bin z value to use for interp
    if ( iz0 == 0 )           // 1st z-bin
      { iz1 = iz0 + 1 ; }

    else if ( iz0 == NBINz-1 ) // last z-bin
      { iz1 = iz0 - 1 ; }

    else if ( zdata > zbin0 ) 
      { iz1 = iz0 + 1 ; }

    else if ( zdata < zbin0 ) 
      { iz1 = iz0 - 1 ; }

    else {
      iz1 = -9 ;
      sprintf(c1err,"Could not determine iz1" );
      sprintf(c2err,"iz0=%d zbin0=%f  NBINz=%d", iz0, zbin0, NBINz);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
    }

    // if interp z-bin does not float M0, then just return
    // constant M0 in this z-bin (Jun 2017)
    if ( FITINP.ISFLOAT_z[iz1] == 0  ) { return(M0); }

    zbin1 = simdata_bias.zM0[iz1];
    M0bin0 = M0LIST[iz0];
    M0bin1 = M0LIST[iz1];

    zfrac = ( zdata - zbin0 ) / ( zbin1 - zbin0) ;
    M0    = M0bin0 + (M0bin1-M0bin0) * zfrac ;

    LDMP = (iz0 >= 8888  ) ;
    if ( LDMP ) {    
      printf(" xxx -------------------------- \n");
      printf(" xxx iz0=%d  iz1=%d \n", iz0, iz1);
      printf(" xxx zdata=%.4f  zbin0=%.4f zbin1=%.4f  zfrac=%.3f\n",
	     zdata, zbin0, zbin1, zfrac );
      printf(" xxx M0bin0=%.3f  M0bin1=%.3f   M0=%.3f\n",
	     M0bin0, M0bin1, M0);
      fflush(stdout);
      //      debugexit(fnam);
    }

  }
  else {
    sprintf(c1err,"Invalid uM0=%d", INPUTS.uM0 );
    sprintf(c2err,"data index n=%d", n);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  return(M0);

} // end fcn_M0

// ==================================================
void fcn_AlphaBetaWGT(double alpha, double beta, int DUMPFLAG,
		      INTERPWGT_ALPHABETA *INTERPWGT, char *callFun ) {

  // Created April 26 2016 by R.Kessler
  //
  // Inputs:  alpha & beta
  //
  // Ouput struct *INTERP with WGT at each alpha,beta corner
  //     used for interpolation in get_muBias.
  //   
  //  This function is called from fcn before the loop over
  //  data in order to reduce the CPU overhead for data event.
  //
  // Sep 4 2017: init INTERPWGT->WGT[ia][ib] = 0.0 
  // Dec 21 2017: widen min,max bound to allow extrapolation

  int    IA, IB, ia, ib;
  int    NBINa     = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int    NBINb     = simdata_bias.BININFO_SIM_BETA.nbin ;
  double a_binSize = simdata_bias.BININFO_SIM_ALPHA.binSize ;
  double b_binSize = simdata_bias.BININFO_SIM_BETA.binSize ;


  double a_bound_min = simdata_bias.BININFO_SIM_ALPHA.avg[0];
  double a_bound_max = simdata_bias.BININFO_SIM_ALPHA.avg[NBINa-1] ;
  double b_bound_min = simdata_bias.BININFO_SIM_BETA.avg[0];
  double b_bound_max = simdata_bias.BININFO_SIM_BETA.avg[NBINb-1] ;

  double alpha_interp, beta_interp ;;
  char fnam[] = "fcn_AlphaBetaWGT" ;

  // ------------------ BEGIN ---------------

  IA = IBINFUN(alpha, &simdata_bias.BININFO_SIM_ALPHA, 2, "" );
  IB = IBINFUN(beta,  &simdata_bias.BININFO_SIM_BETA,  2, "" );

  // get local "interp" values of alpha&beta which do not
  // extend past the storage grid to avoid crazy extrapolations.
  alpha_interp = alpha;
  if ( alpha_interp < a_bound_min ) { alpha_interp = a_bound_min; }
  if ( alpha_interp > a_bound_max ) { alpha_interp = a_bound_max; }

  beta_interp = beta ;  
  if ( beta_interp < b_bound_min ) { beta_interp = b_bound_min; }
  if ( beta_interp > b_bound_max ) { beta_interp = b_bound_max; }

  int FIRST = 1;
  double a_grid, b_grid, Da, Db;
  double SUMWGT = 0.0 ;
  double WGT[MAXBIN_BIASCOR_ALPHA][MAXBIN_BIASCOR_BETA];

  // init WGT map
  for(ia=0; ia < MAXBIN_BIASCOR_ALPHA; ia++ ) {
    for(ib=0; ib < MAXBIN_BIASCOR_BETA; ib++ ) {
      WGT[ia][ib] = 0.0 ;
      INTERPWGT->WGT[ia][ib] = 0.0 ;
    }
  }

  for(ia=IA-1; ia <= IA+1; ia++ ) {
    if ( ia < 0      ) { continue ; }
    if ( ia >= NBINa ) { continue ; }
    for(ib=IB-1; ib <= IB+1; ib++ ) {
      if ( ib < 0      ) { continue ; }
      if ( ib >= NBINb ) { continue ; }

      a_grid = simdata_bias.BININFO_SIM_ALPHA.avg[ia] ;
      b_grid = simdata_bias.BININFO_SIM_BETA.avg[ib] ;
      Da     = fabs(alpha_interp - a_grid) ;
      Db     = fabs(beta_interp  - b_grid) ;

      if ( NBINa > 1 ) { Da /= a_binSize; } else { Da = 0.0 ; }
      if ( NBINb > 1 ) { Db /= b_binSize; } else { Db = 0.0 ; }

      if ( fabs(Da) > 0.99999 ) { continue ; }
      if ( fabs(Db) > 0.99999 ) { continue ; }
      
      WGT[ia][ib] = (1.0 - Da) * (1.0 - Db);
      SUMWGT     += WGT[ia][ib] ;

      if ( FIRST ) { INTERPWGT->ia_min = ia;  INTERPWGT->ib_min = ib; }
      INTERPWGT->ia_max = ia ;
      INTERPWGT->ib_max = ib ;

      FIRST = 0 ;

    } // end ib
  } // end ia



  // ------------------------------
  if ( SUMWGT < 1.0E-9 ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("\t Alpha=%.3f  Alpha_interp=%.3f \n", alpha, alpha_interp);
    printf("\t Beta=%.3f   Beta_interp=%.3f  \n", beta,  beta_interp);

    for(ia=0; ia < MAXBIN_BIASCOR_ALPHA; ia++ ) {
      for(ib=0; ib < MAXBIN_BIASCOR_BETA; ib++ ) {
	printf("\t WGT[ia=%d,ib=%d] = %f \n", ia,ib,INTERPWGT->WGT[ia][ib]);
      }
    }

    sprintf(c1err,"Invalid SUMWGT=%f  (called from %s)", SUMWGT, callFun);
    sprintf(c2err,"IA=%d   IB=%d", IA, IB);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  // divide WGT by SUMWGT so that \Sum WGT = 1
  int NEGWGT=0;
  for(ia=INTERPWGT->ia_min; ia<=INTERPWGT->ia_max; ia++ ) {
    for(ib=INTERPWGT->ib_min; ib<=INTERPWGT->ib_max; ib++ ) {
      INTERPWGT->WGT[ia][ib] = WGT[ia][ib]/SUMWGT ;
      if ( INTERPWGT->WGT[ia][ib] < 0.0 ) { NEGWGT = 1; }
    }
  }


  int LDMP = DUMPFLAG ;
  if ( LDMP || NEGWGT ) {
    int ia_min = INTERPWGT->ia_min ;
    int ia_max = INTERPWGT->ia_max ;
    int ib_min = INTERPWGT->ib_min ;
    int ib_max = INTERPWGT->ib_max ;
    printf("xxx -------------------------------- \n");
    printf("xxx Alpha/Alpha_interp = %f/%f \n", alpha, alpha_interp);
    printf("xxx Beta /Beta_interp  = %f/%f \n", beta,  beta_interp);
    printf("xxx SUMWGT = %le \n", SUMWGT);
    printf("xxx    ia=%d to %d  ib=%d to %d\n"
	   ,ia_min, ia_max, ib_min, ib_max );
    printf("xxx WGT[ia=%d, ib=%d and %d] =%9.6f  %9.6f \n",
	   ia_min, ib_min, ib_max,
	   INTERPWGT->WGT[ia_min][ib_min],
	   INTERPWGT->WGT[ia_min][ib_max] );
    printf("xxx WGT[ia=%d, ib=%d and %d] =%9.6f  %9.6f \n",
	   ia_max, ib_min, ib_max,
	   INTERPWGT->WGT[ia_max][ib_min],
	   INTERPWGT->WGT[ia_max][ib_max] );

    printf("xxx WGT = %f, %f, %f, %f \n",
	   INTERPWGT->WGT[0][0],
	   INTERPWGT->WGT[0][1],
	   INTERPWGT->WGT[1][0],
	   INTERPWGT->WGT[1][1] );
  }
  if ( NEGWGT ) {
    sprintf(c1err,"Invalid Negative weight (called from %s)", callFun);
    sprintf(c2err,"See dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  return ;
  
} // end fcn_AlphaBetaWGT


// ===========================================================
double fcn_muerrsq(char *name, double alpha, double beta, 
		   double (*COV)[NLCPAR], 
		   double z, double zerr, int dumpFlag) {

  // Created Feb 29 2016 by R.Kessler
  //
  // return mu-error for inputs as follows:
  //
  // *name      : name of SN (in case of error msg)
  // alpha,beta : SALT2 standardization parmas for stretch & color
  // COV        : fit cov matrix + intrinsic cov matrix
  // z,zerr     : redshift & its error
  //
  // Jun 27 2017: REFACTOR z bins
  // Sep 04 2017: use matrix-vector summing:
  //    MUERRSQ = (1,a,b) x COV x (1,a,b)
  //
  // Jun 19 2018: pass dumpFlag argument.
  
  double muerrsq, sqtmp, dmuz, dmuLens, VEC[NLCPAR];
  int    i,j ;
  char VECNAME[3][8] = { "ONE  ", "ALPHA", "-BETA" } ;
  
  char fnam[] = "fcn_muerrsq" ;

  // ---------------- BEGIN ------------

  if ( dumpFlag ) {
    printf(" xxx \n");
    printf(" xxx --- %s dump for SNID=%s ----- \n", fnam, name);
    printf(" xxx alpha=%le  beta=%le  z=%le   zerr=%le \n",
	   alpha, beta, z, zerr);
    fflush(stdout);
  }
  

  VEC[0] = +1.0;
  VEC[1] = +alpha ;
  VEC[2] = -beta; 

  muerrsq = 0.0 ;
  for(i=0; i<NLCPAR; i++ ) {
    for(j=0; j<NLCPAR; j++ ) {
      sqtmp    = VEC[j] * COV[i][j] * VEC[i] ;
      muerrsq += sqtmp ;

      if(dumpFlag) {
	printf(" xxx mucov += %13.6le (= %s * %s * %13.6le) \n",
	       sqtmp, VECNAME[j], VECNAME[i], COV[i][j] );
	fflush(stdout);
      }
    }
  }

  if(dumpFlag) { printf(" xxx mucov  = %le from COV \n", muerrsq);  }
  
  // add redshift error and peculiar velocity error  
  dmuz     = fcn_muerrz(1, z, zerr );
  muerrsq += (dmuz*dmuz) ;

  if(dumpFlag) {
    printf(" xxx mucov  = %le with dmu(z,vpec)=%le \n",
			muerrsq, dmuz);   fflush(stdout);
  }
    
  // add lensing error (Sep 9 2016)
  dmuLens  = INPUTS.lensing_zpar * z;
  muerrsq += (dmuLens * dmuLens);

  if(dumpFlag) {
    printf(" xxx mucov  = %le with dmuLens=%le \n",
	   muerrsq, dmuLens) ;
    printf(" xxx mucov(FINAL) = %le (does not include biasScale_muCOV)\n",
	   muerrsq);
    fflush(stdout);
  }

    
  if (muerrsq  <= 0.0 )  {
    sprintf(c1err,"non-positive muerrsq = %le", muerrsq);	
    sprintf(c2err,"for SN = %s  alpha=%f  beta=%f", 
	    name, alpha, beta );   
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return(muerrsq) ;

} // end fcn_muerrsq

// =============================================
double fcn_muerrz(int OPT, double z, double zerr) {

  // Created July 1 2017
  // Compute muerr from redshift (z) and its error (zerr)
  // and optional zpecerr.
  // 
  //
  // OPT=1 --> low-z approx for empty universe
  // OPT=2 --> exact calc.

  // Note that OPT=1 & 2 get closer as z increases.
  //
  // Jan 9 2018: remove zpecerr argment.

  double zerrsq, zpecsq, zerrtot ;
  double muerr = 0.0 ;
  double FAC   = 5.0/LOGTEN ;  
  double *cosPar = &INPUTS.parval[IPAR_OL] ;
  char fnam[]  = "fcn_muerrz" ;

  // --------------- BEGIN --------------

  if ( OPT == 1 ) {

    zerrtot = zerr ;
    muerr    = FAC*(zerrtot/z) * (1.0+z)/(1.0+z/2.0);
  }
  else if ( OPT == 2 )  {

    double dl, mu0, mu1, mu2 ;
    dl     = cosmodl(z,cosPar );
    mu0    = 5.0*log10(dl) + 25.0 ;

    dl     = cosmodl(z-zerr,cosPar);
    mu1    = 5.0*log10(dl) + 25.0 ;

    dl     = cosmodl(z+zerr,cosPar);
    mu2    = 5.0*log10(dl) + 25.0 ;

    muerr = (mu2-mu1)/2.0 ;
  }

  return(muerr);

} // end fcn_muerrz

// ================================================
void test_muerrz(void) {

  double z, muerr1, muerr2;
  double zerr = 0.001 ;
  double zpecerr = 0.0 ;
  char fnam[] = "test_muerrz" ;

  z = 0.01;

  while ( z < 2.5 ) { 
    muerr1 = fcn_muerrz(1, z, zerr);
    muerr2 = fcn_muerrz(2, z, zerr);
    printf(" xxx z=%.3f : muerr(OPT1,2) = %.4f, %.4f  DIF=%0.5f  Ratio=%.3f\n",
	   z, muerr1, muerr2, muerr1-muerr2, muerr1/muerr2 ) ; fflush(stdout);

    if ( z < 0.05 ) 
      { z += 0.002; }
    else if ( z < .5 ) 
      { z += 0.02; }
    else
      { z += 0.1; }
  }

  debugexit(fnam);

} // end test_muerrz


// ============================================
double zerr_adjust(double z, double zerr, double vpecerr, char *name) {

  // Created Jan 9 2018
  // If INPUTS.zpecerr == 0, returns zerr.
  // Otherwise, subtract vpecerr/c in quadrature, and add INPUTS.zpecerr.
  // Thus, zpecerr REPLACES vpec/c.
  //
  // Inputs:
  //   z       = CMB redshift, to apply (1+z) correction on vpecerr
  //   zerr    = zHDERR from input FITRES table
  //   vpecerr = pec-velocity error from input FITRES table
  //   name    = SN name, for error message only
  //
  // May 1 2018:
  //  + z -> CMB redshift
  //  + apply 1E-5 redshift tolerance on abort trap.
  //     --> allows for numerical truncation in FITRES file .
  //
  // Feb 18 2019: abort if arg of sqrt is negative

  double zerr_new = zerr;
  double zpecerr  = vpecerr/LIGHT_km;
  double z1       = (1+z) ;
  double sqz1     = z1*z1 ;
  double sqzerr1, sqzerr2, sqzerr3, sqarg;
  double zdif_tol = 1.0E-5 ;
  char fnam[] = "zerr_adjust" ;
  
  // ---------- BEGIN ------------

  if ( INPUTS.zpecerr < 1.0E-9 ) { return(zerr); }

  sqzerr1 = zerr * zerr;
  sqzerr2 = sqz1 * zpecerr*zpecerr ; // original 
  sqzerr3 = sqz1 * INPUTS.zpecerr * INPUTS.zpecerr ; // user-defined

  // allow tolerance on difference to allow for
  // precision truncation on FITRES file.
  if ( z1*zpecerr > zerr+zdif_tol ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("    z1*zpecerr = %f * %f = %f \n", z1, zpecerr, z1*zpecerr);
    printf("    zerr = %f \n", zerr);
    sprintf(c1err,"Cannot subtract zpecerr from zerr for SNID=%s", name );
    sprintf(c2err,"zerr=%f > (1+z)*zpecerr = (%f)*%f", 
	    zerr, 1+z, zpecerr);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
  }

  sqarg    = (sqzerr1 - sqzerr2 + sqzerr3 );

  if ( sqarg < 0.0 ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("   sqzerr[1,2,3] = %le, %le, %le \n",
	   sqzerr1, sqzerr2, sqzerr3 );
    printf("   user zpecerr = %f \n", INPUTS.zpecerr);
    sprintf(c1err,"sqrt(negatice number); Cannot adjust zHDERR");
    sprintf(c2err,"Check pre-abort dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  zerr_new = sqrt(sqarg );

  return(zerr_new) ;

} // end zerr_adjust 

// =================================================
void fitsc(int ns, double *sbar,double *cbar)
{
  /*
    This function computes the best value of x1 (sbar) 
    and c (cbar) for the SN indexed by ns.
    This function should be called after the solution
    for the fit parameters (alpha, beta, etc) has been
    obtained.  
    The fit factorizes into determination of the global parameters
    and the individual SN parameters, so no information is lost
    by separating the fit this way
    This fit does not use the non-linear terms (zeta, eta, and theta) and
    the results may not be valid if they are non-zero..
    Using the fit values (sbar, and cbar) results in a much smaller
    dispersion around the cosmology, but the SN errors are highly correlated
    and does not result in a reduction in the error on the distance modulus.
  */

  int l, m;
  double I[3];
  double ginv[3][3], gtemp[3][3];
  int indx[10];
  double d;
  int icon;
  double z;
  double alpha, beta;
  int ibin;
  double gi00, gi11, gi22, gi01, gi02, gi12;
  double R[2][2], Rinv[2][2], S[2], det;
  int impar;
  double s, c;
  double xi0, xi1, xi2;
  double chisq=0.0;
  double mag0;

  for (l=0;l<3;++l)
    {
      for (m=0;m<3;++m) gtemp[l][m] = data[ns].covmat_tot[l][m];
    }
  ludcmp(&gtemp[0][0],3,3,indx,&d,&icon);
  for (l=0;l<3;++l)
    {
      for (m=0;m<3;++m) I[m] = 0.0;
      I[l] = 1.0;
      lubksb(&gtemp[0][0],3,3,indx,I);
      for (m=0;m<3;++m) ginv[l][m] = I[m];
    }
  
  z = data[ns].zhd;

  if ( z < 1.0E-8 ) {
    *sbar = -999. ; 
    *cbar = -999. ;
    return ;
  }

  alpha = INPUTS.parval[1] + z*INPUTS.parval[3];
  beta  = INPUTS.parval[2] + z*INPUTS.parval[4];
  
  //Now get mag0 for bin
  if (data[ns].errmask < 4)
    {
      ibin = data[ns].izbin;
      mag0 = INPUTS.parval[ibin+MXCOSPAR];
    }
  else
    { mag0 = FITRESULT.SNMAG0; }

  
  gi00 = ginv[0][0];
  gi11 = ginv[1][1];
  gi22 = ginv[2][2];
  gi01 = ginv[0][1];
  gi02 = ginv[0][2];
  gi12 = ginv[1][2];
  R[0][0] = alpha*alpha*gi00 - 2.0*alpha*gi01 + gi11;
  R[0][1] =  -alpha*beta*gi00 - alpha*gi02 + beta*gi01 + gi12;
  R[1][0] = R[0][1];
  R[1][1] = beta*beta*gi00 + 2.0*beta*gi02 + gi22;
  
  det = R[0][0]*R[1][1] - R[0][1]*R[1][0];
  Rinv[0][0] = R[1][1]/det;
  Rinv[1][1] = R[0][0]/det;
  Rinv[0][1] = -R[1][0]/det;
  Rinv[1][0] = -R[0][1]/det;
  
	  
  xi0 = data[ns].fitpar[INDEX_mB] - mag0 - data[ns].mumodel; 
  xi1 = data[ns].fitpar[INDEX_x1];
  xi2 = data[ns].fitpar[INDEX_c];

  S[0] = (-alpha*gi00+gi01)*xi0 + (-alpha*gi01+gi11)*xi1 
    + (-alpha*gi02+gi12)*xi2;
  S[1] = (beta*gi00+gi02)*xi0 + (beta*gi01+gi12)*xi1
    + (beta*gi02+gi22)*xi2;
  
  
  s = Rinv[0][0]*S[0] + Rinv[0][1]*S[1];
  c = Rinv[1][0]*S[0] + Rinv[1][1]*S[1];
  
 //compute chi-squared
  xi0 = xi0 + alpha*s - beta*c;
  xi1 = data[ns].fitpar[INDEX_x1]-s;
  xi2 = data[ns].fitpar[INDEX_c ]-c;
  chisq = gi00*xi0*xi0 + 2.0*gi01*xi0*xi1 + 2.0*gi02*xi0*xi2;
  chisq += gi11*xi1*xi1 + 2.0*gi12*xi1*xi2 + gi22*xi2*xi2;
  /* Debug 
  printf(" s=%f sbar=%f c=%f cbar=%f chisq=%f\n",
	 data[ns].par[1],s,data[ns].par[2],c,chisq);
  */
  *sbar = s;
  *cbar = c;

  return ;  // (chisq);
   
} // end of fitsc


// ******************************************
void set_defaults(void) {

  int isurvey, order, ipar ;

  // ---------------- BEGIN -----------------
  
  set_EXIT_ERRCODE(EXIT_ERRCODE_SALT2mu);

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );

  sprintf(INPUTS.filename,   "NONE");  
  sprintf(INPUTS.PREFIX,     "NONE" );

  INPUTS.opt_biasCor           = 0 ;
  INPUTS.sigint_biasCor        = -9.0 ; 
  INPUTS.snrmin_sigint_biasCor = BIASCOR_SNRMIN_SIGINT ;
  INPUTS.prescale_biasCor[0]   = 0;
  INPUTS.prescale_biasCor[1]   = 1;
  INPUTS.sigma_cell_biasCor    = 50.0; // flat WGT distribution in each cell
  sprintf(INPUTS.simFile_biasCor,    "NONE");
  sprintf(INPUTS.fieldGroup_biasCor, "NONE" );
  INPUTS.use_fieldGroup_biasCor = 0 ;

  sprintf(INPUTS.surveyGroup_biasCor, "NONE" );
  INPUTS.use_surveyGroup_biasCor = 0 ;

  sprintf(INPUTS.surveyList_noBiasCor, "NONE" );
  INPUTS.idsample_select[0] = 0 ;

  // default is to blind cosmo params for data
  INPUTS.blindFlag = BLINDMASK_FIXPAR; 

  INPUTS.NSPLITRAN = 1; // default is all SN in one job

  INPUTS.iflag_duplicate = IFLAG_DUPLICATE_ABORT ;


  INPUTS.NCUTWIN = 0;

  INPUTS.uM0   = M0FITFLAG_ZBINS_FLAT;
  INPUTS.uzsim = 0 ; // option to set z=simz (i.e., to cheat)

  // stuff for CC prior
  sprintf(INPUTS.simFile_CCprior, "NONE");
  INPUTS.varname_pIa[0] = 0 ;
  INPUTS.simtype_II[0]  = 20 ;
  INPUTS.simtype_II[1]  = 29 ;
  INPUTS.simtype_Ibc[0] = 30 ; 
  INPUTS.simtype_Ibc[1] = 39 ; 
  INPUTS.maxProbCC_for_sigint = 0.2 ;
  INPUTS.typeIa_ccprior       = -9  ;
  INPUTS.zmin = 0.02 ;
  INPUTS.zmax = 1.02 ;
  
  // ---------------------
  // keep obsolete input parameters (4/24/2012 RK)
  INPUTS.x1min = -6.0;
  INPUTS.x1max = +6.0;
  INPUTS.cmin  = -6.0; 
  INPUTS.cmax  = +6.0;
  // ---------------------

  INPUTS.errmask_write = 0;  // default -> write events used in fit

  INPUTS.Nsntype = 0 ;
  sprintf(INPUTS.sntypeString,"NULL");

  INPUTS.nmaxString[0] = 0 ;
  INPUTS.nmax_tot = 999888777 ;
  for(isurvey=0; isurvey<MXIDSURVEY; isurvey++ ) 
    { INPUTS.nmax[isurvey] = 999888777 ; }

  INPUTS.nzbin    =  0 ;
  INPUTS.nlogzbin =  0 ;
  INPUTS.powzbin  =  0.0 ; 
  INPUTS.znhalf   = -9.0 ;
  INPUTS.varname_z[0] = 0 ; // use default z Z zSPEC

  INPUTS.min_per_zbin = MINEVT_PER_ZBIN_DEFAULT ;

  // global scatter matrix
  INPUTS.sigmB  = 0.00 ; // -> 0 on 4/23/2012 by RK
  INPUTS.sigx1  = 0.00 ;
  INPUTS.sigc   = 0.00 ;
  INPUTS.xi01   = 0.0  ;
  INPUTS.xi0c   = 0.0  ;
  INPUTS.xi1c   = 0.0  ;

  INPUTS.sigint_fix[0] = 0 ;
  
  // default abort is for any fitpar error <=0
  INPUTS.maxerr_abort_c=INPUTS.maxerr_abort_x1=INPUTS.maxerr_abort_x0=0.0;

  // redshift dependent scatter matrix
  INPUTS_ZPOLY_COVMAT.NSURVEY = 0;

  for ( isurvey=0; isurvey < MXSURVEY_ZPOLY_COVMAT; isurvey++ ) {
    INPUTS_ZPOLY_COVMAT.IDSURVEY_LIST[isurvey] = -9;
    for ( order=0; order <= ORDER_ZPOLY_COVMAT; order++ ) {
	INPUTS_ZPOLY_COVMAT.sigmB[isurvey][order] = 0.0 ;
	INPUTS_ZPOLY_COVMAT.sigx1[isurvey][order] = 0.0 ;
	INPUTS_ZPOLY_COVMAT.sigc[isurvey][order]  = 0.0 ;
	INPUTS_ZPOLY_COVMAT.xi01[isurvey][order]  = 0.0 ;
	INPUTS_ZPOLY_COVMAT.xi0c[isurvey][order]  = 0.0 ;
	INPUTS_ZPOLY_COVMAT.xi1c[isurvey][order]  = 0.0 ;
    }
  }


  //Defaults
  INPUTS.zpecerr = 0.;     // error on v/c
  INPUTS.lensing_zpar = 0.0 ;

  INPUTS.fitflag_sigmb       = 0;     // option to repeat fit until chi2/dof=1
  INPUTS.redchi2_tol         = 0.02;  // tolerance on chi2.dof
  INPUTS.sigmb_step1         = 0.05 ; // size of 1st step for sigMB, OR ...
  INPUTS.scale_covint_step1  = 0.04 ; // for scale_covint

  INPUTS.prescale_simData  = 1.0 ; // include all simData by default
  INPUTS.prescale_simCC    = 1.0 ; // include all simcc by default.

  DOFIT_FLAG = FITFLAG_CHI2 ;
  INPUTS.zpolyflag = 0;
  //Default distance modulus parameters
  INPUTS.H0      = 70.0;
  INPUTS.nommag0 = M0_DEFAULT ;
  INPUTS.uave    = 1;

  // Default parameters 

  double a0 = FITPARBOUND_ALPHA[0];
  double a1 = FITPARBOUND_ALPHA[1];
  double b0 = FITPARBOUND_BETA[0];
  double b1 = FITPARBOUND_BETA[1];

  //         ipar  val   step   bnd0  bnd1 float
  set_fitPar(  1,  0.13, 0.01,  a0,   a1,    1 ); // alpha
  set_fitPar(  2,  3.20, 0.30,  b0,   b1,    1 ); // beta
  set_fitPar(  3,  0.0,  0.02, -0.50, 0.50,  0 ); // dAlpha/dz
  set_fitPar(  4,  0.0,  0.10, -3.00, 3.00,  0 ); // dBeta/dz
  set_fitPar(  5,  0.0,  0.01, -0.50, 0.50,  0 ); // gamma0
  set_fitPar(  6,  0.0,  0.01, -2.00, 2.00,  0 ); // gamma1
  set_fitPar(  7, 10.0,  0.10,  9.00,11.00,  0 ); // logmass_cen
  set_fitPar(  8,  0.02, 0.01,  0.001,1.00,  0 ); // logmass_tau
  set_fitPar(  9,  0.73, 0.05,  0.00, 2.00,  0 ); // Omega_Lam
  set_fitPar( 10,  0.0,  0.05, -2.00, 5.00,  0 ); // Omega_k
  set_fitPar( 11, -1.0,  0.1,  -3.00, 1.00,  0 ); // w0
  set_fitPar( 12,  0.0,  0.5,  -8.00, 8.00,  0 ); // wa
  set_fitPar( 13,  1.0,  0.05, -0.10, 5.00,  0 ); // scale_PCC
  set_fitPar( 14,  0.1,  0.01,  0.02, 0.50,  0 ); // sigint
  set_fitPar( 15,  0.0,  0.02, -0.50, 0.50,  0 ); // dAlpha/dlogMhost
  set_fitPar( 16,  0.0,  0.10, -3.00, 3.00,  0 ); // dBeta/dlogMhost

  // set IPAR_XXX params to access FITRESULTS.
  IPAR_ALPHA0=1, IPAR_BETA0=2;
  IPAR_GAMMA0=5; IPAR_GAMMA1=6;  IPAR_LOGMASS_CEN=7; IPAR_LOGMASS_TAU=8;
  IPAR_scalePCC=13, IPAR_COVINT_PARAM=14 ;
  IPAR_OL=9; IPAR_Ok=10; IPAR_w0=11;  IPAR_wa=12;
  IPAR_H11=17; NPAR_H11_TOT=6; NPAR_H11_USER=0;

  // H11 params for CC prior (April 15 2016)
  //         ipar  val   step   bnd1  bnd2 float
  set_fitPar( 17, +1.2,  0.1,  -9.0,  9.0,  0 ); // H11mucc0
  set_fitPar( 18, +3.0,  0.2,  -9.0,  9.0,  0 ); // H11mucc1
  set_fitPar( 19, -2.0,  0.5,  -30.0,30.0,  0 ); // H11mucc2
  set_fitPar( 20,  0.3,  0.05,  0.01, 2.0,  0 ); // H11sig0
  set_fitPar( 21,  0.0,  0.05, -9.0,  9.0,  0 ); // H11sig1
  set_fitPar( 22,  0.0,  0.05, -9.0,  9.0,  0 ); // H11sig2

  for(ipar=0; ipar < MAXPAR; ipar++ )  {  INPUTS.izpar[ipar] = -99; }

  // ======== misc ======
  INPUTS.NDUMPLOG = 1000 ;
  INPUTS.SNID_MUCOVDUMP[0] = 0 ;
  
  // === set blind-par values to be used if blindflag=2 (Aug 2017)
  INPUTS.blindpar[IPAR_OL][0] = 0.06; 
  INPUTS.blindpar[IPAR_OL][1] = 23434. ;

  INPUTS.blindpar[IPAR_w0][0] = 0.20 ; 
  INPUTS.blindpar[IPAR_w0][1] = 8432. ;

  sprintf(INPUTS.varname_gamma,"HOST_LOGMASS");
  INPUTS.USE_GAMMA0  = 0 ;
  //  INPUTS.USE_LOGMASS = 0 ;

  return ;

}   // end of set_defaults


// ============================================
void set_fitPar(int ipar, double val, double step, 
		double bndmin, double bndmax, int use) {
  INPUTS.parval[ipar]    = val ;
  INPUTS.parstep[ipar]   = step ;
  INPUTS.parbndmin[ipar] = bndmin;  
  INPUTS.parbndmax[ipar] = bndmax ;
  INPUTS.ipar[ipar]      = use ;
  INPUTS.izpar[ipar]     = -9; // init

  INPUTS.blindpar[ipar][0] = 0.0 ;
  INPUTS.blindpar[ipar][1] = 0.0 ;

}

// ******************************************
int read_data(char *filnam)
{

  // Read data from input fitres-formatted file.
  // April 2012 - read optional CUTVAR variables for CUTWIN options
  // Jan   2016 - don't print each fudged error matrix unless LDMP=1
  // Jul   2016 - read integer (instead of double) for
  //              SNTYPE, IDSURYEY, SIM_NON1A_INDEX
  //            - Read FIELD if INPUTS.fieldSplit is specified
  //
  // Jan  8 2018: read VPEC and VPEC_ERR. 
  // Jan  9 2018: read ZHEL and ZHEL_ERR. 
  // Mar 21 2018: read optional opt_photoz
  // Jun 19 2018: abort if any pIa < 0 or pIa>1
  // Jun 26 2018: if SNID_MUCOVDUMP=snid, set dump option for update_covMatrix
  // Jan 20 2019: specify SIM_TEMPLATE_INDEX:I instead of SIM_TEMPLATE_INDEX
  //    
  double COVINT_PARAM_FIX_SAVE = FITINP.COVINT_PARAM_FIX ;
  double sf, covfit, covtot, SCALE, z, zhd, zcmb, zpec, dl;
  double zerr, vpecerr, zpecerr, pIa ;
  double avex1, avec, rmsx1, rmsc, cx1c, dc, dx1, x0tmp;
  double covMat_fit[NLCPAR][NLCPAR];
  double covMat_int[NLCPAR][NLCPAR];

  int idsurvey, icut, RDFLAG, ABORTFLAG, NCOVFIX, OPTMASK ;
  int vbose, noABORT, ivar, l, m, n, istat_cov, dumpFlag_covMatrix;
  int FOUND_VPEC=1, FOUND_ZHEL=1, FOUND_IDSURVEY=1, FOUND_OPT_PHOTOZ=1;
  int NERR_pIa=0;
  
  char vartmp[80], *cutname, str_z[40], str_zerr[40], *name;
  char fnam[] = "read_data" ;

  // -------------- BEGIN ---------------

   sprintf(BANNER,"%s", fnam);
   print_banner(BANNER);


   // -----------------------------------------------
   // prepare reading file 
   TABLEFILE_INIT();
   int IFILETYPE = TABLEFILE_OPEN(filnam,"read");
   NVAR_ORIG     = SNTABLE_READPREP(IFILETYPE,"FITRES");
   SNDATA_INFO.IFILETYPE = IFILETYPE ;

   // store VARNAMES_ORIG
   for(ivar=1; ivar <= NVAR_ORIG; ivar++ ) {
     sprintf(VARNAMES_ORIG[ivar-1], "%s",
	     READTABLE_POINTERS.VARNAME[ivar-1] );
   }

   MAXSN = rawdata_malloc(filnam);
   rawdata_init();

  
  FOUNDKEY_SNTYPE = 0 ;
  vbose = 3 ; // verbosity flag for fitresVar_init()
  noABORT = 1; // verbosity flag, but don't abort if var is missing

  // --------
  // cast CID as char to allow strings
  sprintf(vartmp, "CID:C*%d  CCID:C*%d", MXCHAR_CCID, MXCHAR_CCID );  
  SNTABLE_READPREP_VARDEF(vartmp, rawdata.name, MAXSN,vbose) ;

  if ( INPUTS.use_fieldGroup_biasCor ) {
    sprintf(vartmp, "FIELD:C*%d", MXCHAR_CCID);  
    SNTABLE_READPREP_VARDEF(vartmp, rawdata.field, MAXSN,vbose) ;
  }

  // allow different redshift keys
  str_z[0] = str_zerr[0] = 0 ;
  get_zString(str_z,str_zerr, "D");
  SNTABLE_READPREP_VARDEF(str_z,    rawdata.zhd,    MAXSN,vbose);
  SNTABLE_READPREP_VARDEF(str_zerr, rawdata.zhderr, MAXSN,vbose) ;

  // Jan 9 2018: read ZHEL and its error
  sprintf(vartmp,"zHEL");
  ivar = SNTABLE_READPREP_VARDEF(vartmp,rawdata.zhel,    MAXSN, noABORT);
  if ( ivar < 0 ) { FOUND_ZHEL=0; }
  sprintf(vartmp,"zHELERR");
  ivar = SNTABLE_READPREP_VARDEF(vartmp,rawdata.zhelerr, MAXSN, noABORT);
  if ( ivar < 0 ) { FOUND_ZHEL=0; }

  // Jan 8 2018: read VPEC and its error
  sprintf(vartmp,"VPEC");
  ivar = SNTABLE_READPREP_VARDEF(vartmp,rawdata.vpec,    MAXSN, noABORT);
  if ( ivar < 0 ) { FOUND_VPEC = 0; }
  sprintf(vartmp,"VPECERR VPEC_ERR") ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp,rawdata.vpecerr, MAXSN, noABORT);
  if ( ivar < 0 ) { FOUND_VPEC = 0; }

  sprintf(vartmp,"x0");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x0,    MAXSN,vbose);
  sprintf(vartmp,"x0ERR");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x0err, MAXSN,vbose) ;

  sprintf(vartmp,"x1");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x1,   MAXSN,vbose) ;
  sprintf(vartmp,"x1ERR");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x1err,MAXSN,vbose) ;

  sprintf(vartmp,"c");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x3,    MAXSN,vbose) ;
  sprintf(vartmp,"cERR");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x3err, MAXSN,vbose) ;


  //Link correlations
  sprintf(vartmp,"COVx0x1 COV_x1_x0");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x0x1_c,MAXSN,vbose);

  sprintf(vartmp,"COVx0c COV_c_x0");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x0x3_c,MAXSN,vbose) ;

  sprintf(vartmp,"COVx1c COV_x1_c");
  SNTABLE_READPREP_VARDEF(vartmp,rawdata.x1x3_c,MAXSN,vbose) ;


  // check option to read p1a for CCprior (Feb 2016)
  if ( simdata_ccprior.USE ) {

    sprintf(vartmp,"%s", INPUTS.varname_pIa);
    SNTABLE_READPREP_VARDEF(vartmp,rawdata.pIa,MAXSN,vbose) ;

    /* xxxxxxxxx mark delete Feb 5 2019 xxxxxxxxxxxxxx
    sprintf(vartmp,"%s", INPUTS.varname_pII);
    if ( strlen(vartmp) > 0 ) 
      { SNTABLE_READPREP_VARDEF(vartmp,rawdata.pII,MAXSN,vbose) ; }

    sprintf(vartmp,"%s", INPUTS.varname_pIbc);
    if ( strlen(vartmp) > 0 ) 
      { SNTABLE_READPREP_VARDEF(vartmp,rawdata.pIbc,MAXSN,vbose) ; }
    xxxxxxxxxxx */

  }
  fflush(stdout); 

  // read optional survey ID
  sprintf(vartmp,"IDSURVEY:I");
  if ( SNTABLE_READPREP_VARDEF(vartmp,rawdata.idsurvey,MAXSN,1) < 0 ) 
    { FOUND_IDSURVEY = 0; }

  // read optional PHOTOZ flag
  sprintf(vartmp,"OPT_PHOTOZ:I");
  if ( SNTABLE_READPREP_VARDEF(vartmp,rawdata.opt_photoz,MAXSN,1) >= 0 ) 
    { FOUNDKEY_OPT_PHOTOZ = 1 ; }


  //Look for SN type -- initially assume not found
  ivar = SNTABLE_READPREP_VARDEF("TYPE:I",rawdata.sntype,MAXSN,1) ;
  if ( ivar > 0 )  { FOUNDKEY_SNTYPE = 1; }
 

  //Simulated values
  ivar = SNTABLE_READPREP_VARDEF("SIM_ZCMB", rawdata.simz, MAXSN,1);
  if ( ivar > 0 )  { FOUNDKEY_SIM = 1; }

  ivar = SNTABLE_READPREP_VARDEF("SIM_DLMAG",rawdata.sim_mu,MAXSN,1) ;
  if ( ivar > 0 ) { FOUNDKEY_SIM = 1; }

  ivar = SNTABLE_READPREP_VARDEF("SIMx0 SIM_x0",   rawdata.simx0, MAXSN,1) ;
  if ( ivar > 0 ) { FOUNDKEY_SIMx0 = 1; }


  // Jul 31 2018: read either "NONIA" or "TEMPLATE"
  ivar = SNTABLE_READPREP_VARDEF("SIM_NONIA_INDEX:I  SIM_TEMPLATE_INDEX:I",
				 rawdata.sim_nonIa_index,MAXSN,1) ; 
  if ( ivar > 0 ) { FOUNDKEY_SIM_NONIA_INDEX = 1; }

  SNTABLE_READPREP_VARDEF("SIMx1 SIM_x1 S2x1",rawdata.simx1,MAXSN,1);
  SNTABLE_READPREP_VARDEF("SIMc  SIM_c  S2c", rawdata.simc, MAXSN,1);
 
  // April 2012 RK - read optional variables for CUTWIN
  
  for ( icut=0; icut < INPUTS.NCUTWIN; icut++ ) { 
    cutname   = INPUTS.CUTWIN_NAME[icut] ;
    RDFLAG    = INPUTS.CUTWIN_RDFLAG[icut] ;
    ABORTFLAG = INPUTS.CUTWIN_ABORTFLAG[icut] ;
    rawdata.DOFLAG_CUTWIN[icut] = 1;
    if ( RDFLAG ) { 
      ivar = SNTABLE_READPREP_VARDEF(cutname, rawdata.CUTVAL[icut], MAXSN, 1);
      rawdata.DOFLAG_CUTWIN[icut] = set_DOFLAG_CUTWIN(ivar,icut,1);      
      if ( strcmp(cutname,INPUTS.varname_gamma) == 0 ) 
	{ rawdata.ICUTWIN_GAMMA = icut ; }
    }
  }
  

  fflush(stdout);
  printf("\n");

  if (FOUNDKEY_SNTYPE) 
    { printf("\t SN TYPE is present.\n"); }
  else 
    { printf("\t SN TYPE is not present \n"); }

  if (FOUNDKEY_SIM ) 
    { printf("\t Simulated values are present.\n"); ISDATA=0; }
  else
    { printf("\t Simulated values are not present.\n"); ISDATA=1; }


  fflush(stdout);
  FITINP.NSNTOT = SNTABLE_READ_EXEC();

  // xxx obsolete  TABLEFILE_CLOSE(filnam);  // added Jan 19 2016

  if ( FITINP.NSNTOT > MAXSN) {
    sprintf(c1err,"NSNTOT = %d exceeds bound of MAXSN=%d",
	    FITINP.NSNTOT, MAXSN) ;
    sprintf(c2err,"increase MAXSN or check input file size.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // ------ check prescale_simcc option
  if ( INPUTS.prescale_simCC > 1.0 && FOUNDKEY_SIM_NONIA_INDEX == 0 ) {
    sprintf(c1err,"prescale_simcc=%f", INPUTS.prescale_simCC );
    sprintf(c2err,"but cant find requred SIM_TEMPLATE_INDEX variable");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
  }
	    // ---------
  // initialize ERRMASK counters
  for  ( n=0; n < 10000; n++ ) { NSET_ERRMASK[n] = 0; }


  // init IDSURVEY-dependent arrays that are filled below
  for(n=0; n < MXIDSURVEY; n++ ) { 
    SNDATA_INFO.NEVT_PER_SURVEY[n] = 0; 
    SNDATA_INFO.zMIN[n] = +999.0 ;
    SNDATA_INFO.zMAX[n] = -999.0 ;
    SNDATA_INFO.sigint_fix[n]   = 0.0 ;
    SNDATA_INFO.sqsigint_fix[n] = 0.0 ;
  }

  FITINP.NSNCUTS  = 0;
  FITINP.NSNERR   = 0 ;
  NCOVFIX = 0 ;
  int nsn = 0 ;
  int errmask ;
  // beware than nsn must equal n in loop below

  for (n = 0; n < FITINP.NSNTOT; ++n)
    {
      sprintf(data[nsn].name, "%s", rawdata.name[n] ); 
      sscanf( data[nsn].name, "%d", &data[nsn].snid ); // problem with string?
      idsurvey = (int)data[nsn].idsurvey ;

      if( INPUTS.use_fieldGroup_biasCor) 
	{ sprintf(data[nsn].field, "%s", rawdata.field[n] );  }
      else
	{ sprintf(data[nsn].field, "notRead" );  }

      data[nsn].errmask  = 0;
      data[nsn].warn     = 0;

      // July 2016: transfer a few integers to double arrays
      // to allow using in cutvar array.
      rawdata.dsurvey[n] = (double)rawdata.idsurvey[n];
      rawdata.dtype[n]   = (double)rawdata.sntype[n];
    
      set_ERRMASK(n,nsn);

      // check errmask reasons to skip this in fit (Mar 2016)
      data[nsn].skipfit = 0 ;
      errmask = data[nsn].errmask ;
      if ( errmask > 0 )  { data[nsn].skipfit=1; } // Jun 23 2017      
      
      data[nsn].idsurvey   = rawdata.idsurvey[n];

      if ( FOUNDKEY_OPT_PHOTOZ )
	{ data[nsn].opt_photoz = rawdata.opt_photoz[n]; }
      else
	{ data[nsn].opt_photoz = 0; }
 
	//simulated values
      if ( FOUNDKEY_SIM )
	{
	  data[nsn].simz   = rawdata.simz[n];    // ZCMB
	  data[nsn].sim_mu = rawdata.sim_mu[n];

	  if ( data[nsn].simz < 0 ) {
	    sprintf(c1err,"Invalid SIMz = %f ",   data[nsn].simz );
	    sprintf(c2err,"for SNID='%s'  n=%d ", data[nsn].name, nsn);
	    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
	  }

	
	  x0tmp = rawdata.simx0[n] ;
	  if ( x0tmp > 0.0 ) 
	    { data[nsn].sim_mb = -2.5*log10(x0tmp); }
	  else
	    { data[nsn].sim_mb = -9.0; } // not a SNIa

	  data[nsn].sims = rawdata.simx1[n];
	  data[nsn].simc = rawdata.simc[n];
	  data[nsn].sim_nonIa_index = rawdata.sim_nonIa_index[n];
	}
      else
	{
	  data[nsn].simz = 0.0;
	  data[nsn].sim_mu = 0.0;
	  data[nsn].sim_mb = 0.0;
	  data[nsn].sims = 0.0;
	  data[nsn].simc = 0.0;
	} // end of FOUNDKEY_SIM

      if (FOUNDKEY_SNTYPE) 
	{ data[nsn].sntype = rawdata.sntype[n]; }
      else 
	{ data[nsn].sntype = -1; }

      //redshift
      data[nsn].zhd    = rawdata.zhd[n];
      data[nsn].zhderr = rawdata.zhderr[n];

      // load zhel only if it was read in
      data[nsn].zhel = data[nsn].zhelerr = 0.0 ;
      if ( FOUND_ZHEL )   { 
	data[nsn].zhel    = rawdata.zhel[n]; 
	data[nsn].zhelerr = rawdata.zhelerr[n];
      }

      // load vpec only if it was read in
      data[nsn].vpec = data[nsn].vpecerr = 0.0 ;
      if ( FOUND_VPEC )   { 
	data[nsn].vpec    = rawdata.vpec[n]; 
	data[nsn].vpecerr = rawdata.vpecerr[n];
      }


      // Jan 8 2018: if user-input zpecerr > 0, then subtract out 
      //    original zpecerr and add user input
      // zhd = (1+zcmb)*(1+zpec) -1 = zcmb + zpec + zcmb*zpec
      // zcmb(1+zpec) = zhd - zpec
      
      zhd     =  data[nsn].zhd ;
      zpec    =  data[nsn].vpec/LIGHT_km ;
      zcmb    =  (zhd - zpec)/(1.0+zpec) ;
      zerr    =  data[nsn].zhderr;  vpecerr= data[nsn].vpecerr;
      name    =  data[nsn].name ;
      data[nsn].zhderr = zerr_adjust(zcmb, zerr, vpecerr, name);
      
      // if cosmo params are fixed, then store mumodel for each event.
      if ( INPUTS.FLOAT_COSPAR == 0 ) {
	z = data[nsn].zhd ;
	dl                = cosmodl_forFit(z,INPUTS.COSPAR);
	data[nsn].mumodel = 5.0*log10(dl) + 25.0 ;
      }

      // check z-cheat option (Jun 2016)
      if ( INPUTS.uzsim && FOUNDKEY_SIM ) { 
	data[nsn].zhd    = data[nsn].simz ;   // ZCMB
	data[nsn].zhderr = 1.0E-7 ;
      }

      if ( data[nsn].zhd < 0 ) {
	sprintf(c1err,"Invalid redshift z = %f ", data[nsn].zhd );
	sprintf(c2err,"for SNID='%s'  n=%d ", data[nsn].name, nsn);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
      }

      //input parameters x0=scale, x1=stretch, x3=color
      data[nsn].fitpar[INDEX_mB] = -2.5*log10(rawdata.x0[n]);
      data[nsn].fitpar[INDEX_x1] = rawdata.x1[n];
      data[nsn].fitpar[INDEX_c ] = rawdata.x3[n];

      if ( simdata_ccprior.USE )  { 
	data[nsn].pIa  = rawdata.pIa[n]; 

	/* xxxxxxx mark delete Feb 5 2019 xxxxxxx
	data[nsn].pII  = rawdata.pII[n]; 
	data[nsn].pIbc = rawdata.pIbc[n];
	xxxxxxxxxx */

	pIa = rawdata.pIa[n];
	if ( pIa < -1.0E-12 || pIa > 1.00000001 ) {
	  NERR_pIa++ ;
	  if ( NERR_pIa < 10 ) {
	    printf("ERROR: Invalid pIa = %le for SNID=%s \n", pIa, name);
	    fflush(stdout);
	  }
	}
      }


      //scale factor for errors after taking log
      // Mar 31, 2010 RK - include minus sign for sf
      sf = -2.5/(rawdata.x0[n]*LOGTEN) ;
      covMat_fit[INDEX_mB][INDEX_mB] = rawdata.x0err[n]*rawdata.x0err[n]*sf*sf;
      covMat_fit[INDEX_mB][INDEX_x1] = rawdata.x0x1_c[n]*sf;
      covMat_fit[INDEX_mB][INDEX_c]  = rawdata.x0x3_c[n]*sf;
      covMat_fit[INDEX_x1][INDEX_x1] = rawdata.x1err[n]*rawdata.x1err[n];
      covMat_fit[INDEX_x1][INDEX_c]  = rawdata.x1x3_c[n];
      covMat_fit[INDEX_c][INDEX_c]   = rawdata.x3err[n]*rawdata.x3err[n];
      covMat_fit[INDEX_x1][INDEX_mB] = covMat_fit[INDEX_mB][INDEX_x1];
      covMat_fit[INDEX_c][INDEX_mB]  = covMat_fit[INDEX_mB][INDEX_c];
      covMat_fit[INDEX_c][INDEX_x1]  = covMat_fit[INDEX_x1][INDEX_c];

      OPTMASK=0;
      if ( strcmp(name,INPUTS.SNID_MUCOVDUMP) == 0 ) { OPTMASK=4; } //dump

      update_covMatrix(data[n].name, OPTMASK, NLCPAR, covMat_fit, 
		       EIGMIN_COV, &istat_cov ); 
      if ( istat_cov != 0 ) { data[n].warn += 1;  NCOVFIX++; }

      // fill z-dependent scatter matrix, ZPOLY_COVMAT
      idsurvey = (int)data[nsn].idsurvey ;
      z        = data[nsn].zhd ;
      load_ZPOLY_COVMAT(idsurvey, z );

      // get user-defined intrinsic cov to get initial COVTOT below

      FITINP.COVINT_PARAM_FIX  = INPUTS.sigmB ;
      get_COVINT_model(-1,covMat_int); 
      FITINP.COVINT_PARAM_FIX = COVINT_PARAM_FIX_SAVE; 

      //Add approx intrinsic scatter to error matrix
      for (l=0; l<NLCPAR; ++l) {
        for (m=0; m<NLCPAR; ++m) {
	  covfit = covMat_fit[l][m] ;
	  covtot = covMat_int[l][m] + covfit ;
          data[nsn].covmat_fit[l][m] = covfit ;
	  data[nsn].covmat_tot[l][m] = covtot ;
        }
      }

      // update survey-dependent quantities needed for biasCor samples
      // and binning (Aug 23 2016)
      SNDATA_INFO.NEVT_PER_SURVEY[idsurvey]++ ;
      if (z<SNDATA_INFO.zMIN[idsurvey] ) { SNDATA_INFO.zMIN[idsurvey]=z;}
      if (z>SNDATA_INFO.zMAX[idsurvey] ) { SNDATA_INFO.zMAX[idsurvey]=z;}

      nsn++ ; // increment global counter

    } // end big loop over sn



  // --------------------
  if ( simdata_ccprior.USE && NERR_pIa > 0 ) {
    sprintf(c1err,"Found %d Invalid pIa values.", NERR_pIa);
    sprintf(c2err,"Required range is 0 <= pIa <= 1");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
	  
  FITINP.NSNCUTS = nsn ;

  // ----------------------------------
  // set initial muerr so that optional 1/sigma term 
  // can be fixed on first iteration
  int    isn;
  double muerrsq ;
  double alpha = INPUTS.parval[IPAR_ALPHA0] ;
  double beta  = INPUTS.parval[IPAR_BETA0] ;


  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {

    if ( data[isn].skipfit ) { 
      muerrsq = 9.0 ;
    }
    else {
      muerrsq = 
	fcn_muerrsq (data[isn].name, alpha, beta ,
		     data[isn].covmat_tot, data[isn].zhd, data[isn].zhderr,0);
    }
    data[isn].muerrsq_last = muerrsq ;
    data[isn].muerr_last   = sqrt(muerrsq);

    z = data[isn].zhd ;
    double *ptr_sigCC = &INPUTS.parval[IPAR_H11+3];
    double sigCC      = ptr_sigCC[0] + z*ptr_sigCC[1] + z*z*ptr_sigCC[2];
    data[isn].sigCC_last   = sigCC ;
    data[isn].sqsigCC_last = sigCC * sigCC ;
  }

  // ----------

  printf("\n");
  printf(" Prepared %i SN (%d with non-zero error flag). \n", 
	 FITINP.NSNCUTS, FITINP.NSNERR); 
  fflush(stdout);

  // dump number of errors per ERRMASK bit
  int NERR, mask ;
  char comment[200];
  for ( mask=0; mask<1000; mask++ ) {
    NERR = NSET_ERRMASK[mask];
    if ( NERR > 0 ) {
      getComment_ERRMASK(mask,comment) ;
      printf("\t NERR[MASK = %4.4d] = %3d  (%s) \n", 
	     mask, NERR, comment) ;
    }
  }
  printf("\t Error matrix fudged for %d SN\n", NCOVFIX);
  printf("\n");
  fflush(stdout);

  int nfit = (FITINP.NSNCUTS - FITINP.NSNERR) ;
  if ( nfit <= INPUTS.min_per_zbin ) {
    sprintf(c1err,"Only %d SN with non-zero error flag.", nfit);
    sprintf(c2err,"Check input file and cuts.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if (FITINP.NSNCUTS == 0)  { 
    sprintf(c1err,"zero SN pass cuts ???");
    sprintf(c2err,"Check input file and cuts.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  avex1 = 0.0;
  avec = 0.0;
  for (n=0; n< FITINP.NSNCUTS; ++n)    {
    avex1 += data[n].fitpar[INDEX_x1];
    avec  += data[n].fitpar[INDEX_c];
  }
  avex1 /= FITINP.NSNCUTS ;
  avec  /= FITINP.NSNCUTS ;

  rmsx1 = 0.0;
  rmsc = 0.0;
  cx1c = 0.0;
  for (n=0;n<FITINP.NSNCUTS; ++n)
    {
      dx1    = data[n].fitpar[INDEX_x1] - avex1;
      dc     = data[n].fitpar[INDEX_c] - avec;
      rmsx1 += dx1*dx1;
      rmsc  += dc*dc;
      cx1c  += dx1*dc;
    }
  rmsx1 = sqrt(rmsx1/FITINP.NSNCUTS);
  rmsc = sqrt(rmsc/FITINP.NSNCUTS);
  cx1c = cx1c/(FITINP.NSNCUTS*rmsx1*rmsc);

  printf("\t Mean x1=%.3f Mean c=%.4f \n", avex1,avec);
  printf("\t Rms  x1=%.3f Rms  c=%.4f Correlation=%6.3f \n",
	 rmsx1,rmsc,cx1c);

  fflush(stdout);

   return(0);

} // end of read_data

 

// ==================================================
void get_zString(char *str_z, char *str_zerr, char *cast) {

  // Created Jan 2016
  // Return redshift strings to pass to SNTABLE_READPREP_VARDEF .
  // These strings contain default names, or a specific name
  // specified in the input.
  //
  // May 30 2016: fix bug and use zHD,zHDERR instead of z,zERR= zHELIO.
  // Sep 05 2017: pass cast as input argument

  if ( strlen(INPUTS.varname_z) == 0 ) {
    sprintf(str_z,    "zHD:%s zSPEC:%s", cast, cast);
    sprintf(str_zerr, "zHDERR:%s zSPECERR:%s", cast, cast);
  }
  else {
    sprintf(str_z,    "%s:%s",    INPUTS.varname_z, cast);
    sprintf(str_zerr, "%sERR:%s", INPUTS.varname_z, cast);
  }

  return ;

} // end get_zString


// =========================================================
void prepare_IDSAMPLE_biasCor(void) {

  // Created July 18 2016
  // Set internal map index for each SURVEYGROUP & FIELDGROUP.
  // Called after reading data, but before reading biasCor and CCprior.
  // Used later to determine biasCor and CCprior for sub-samples.
  //
  // Example 1:  
  //    IDSURVEY list in data file = 54(CFA4), 1(SDSS), 10(DES)
  //  --> MAPINDEX_SURVEYFIELD = 0,1,2  for CFA4,SDSS,DES
  //
  // Example 2: 
  //   IDSURVEY list = 54(CFA4) and 1(SDSS)
  //   fieldGroup_biasCor = '82N,82S'
  // --> MAPINDEX_SURVEYFIELD = 0 [CFA4]
  //                          = 1 [SDSS-82N]
  //                          = 2 [SDSS-82S]
  //  
  // If fieldGroup_biasCor='A+B+C,X+Y+Z', then there are
  // two groups: first group as A and B and C, and the 
  // second groupd has X and Y and Z.
  // Also note that fieldGroup requires only a partial match.
  // Thus, suppose there are fields MEDIUM-01, MEDIUM-02, ... MEDIDUM-50,
  // DEEP-01 ... DEEP-10. ;
  //    fieldGroup_biasCor='MEDIUM,DEEP' 
  // is sufficient to split the biasCor sample into two groups.
  //
  // Apr 4 2017: few fixes for "ALL" option (set redshifts)
  //
  // Mar 20 2018: if OPT_PHOTOZ exists, split by zSPEC and zPHOT 

  int USE_FIELDGROUP  = INPUTS.use_fieldGroup_biasCor ;
  int isn, IDSURVEY, OPT_PHOTOZ, N, IDSAMPLE, i, NIDSURVEY[MXIDSURVEY] ;
  int  DUMPFLAG=0, NDMP = 0  ; 
  double z;
  char FIELD_TMP[MXCHAR_CCID],  FIELDGROUP[100],  FIELDDEF[MXCHAR_CCID];
  char SURVEYGROUP[100], SURVEYDEF[MXCHAR_CCID], zGROUP[20];
  char SNID[MXCHAR_CCID],  *ptrName, STRINGOPT[40]  ;
  char fnam[] = "prepare_IDSAMPLE_biasCor" ;

  // --------------- BEGIN -----------------

  NSAMPLE_BIASCOR      = 0 ;
  ONE_SAMPLE_BIASCOR   = 0 ;
  INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR  = 0 ;
  INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR = 0 ;
  INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT = 0 ;

  if ( INPUTS.opt_biasCor == 0 ) { return ; }

  sprintf(BANNER,"Begin %s", fnam);
  print_banner(BANNER);

      
  for(i=0; i < MXNUM_SAMPLE; i++ ) { 
    SAMPLE_BIASCOR[i].NDATA    = 0 ; 
    SAMPLE_BIASCOR[i].NBIASCOR = 0 ; 
    SAMPLE_BIASCOR[i].NCCPRIOR = 0 ; 
    SAMPLE_BIASCOR[i].DOFLAG_SELECT  = 1 ; // default -> select sample
    SAMPLE_BIASCOR[i].DOFLAG_BIASCOR = 1 ; // default -> do Bias cor
    SAMPLE_BIASCOR[i].IFLAG_ORIGIN   = 0 ;
    SAMPLE_BIASCOR[i].zMIN = +999. ;
    SAMPLE_BIASCOR[i].zMAX = -999. ;

    SAMPLE_BIASCOR[i].NAME_SURVEYGROUP[0] = 0 ;
    SAMPLE_BIASCOR[i].NAME_FIELDGROUP[0]  = 0 ;
    SAMPLE_BIASCOR[i].STRINGOPT[0]        = 0 ;
    // xxx    SAMPLE_BIASCOR[i].NAME_zGROUP[0]      = 0 ;
    SAMPLE_BIASCOR[i].OPT_PHOTOZ          = 0 ; // zSPEC is default
    SAMPLE_BIASCOR[i].NAME[0]        = 0 ;

  } // end i loop 


  for(i=0; i<MXIDSURVEY; i++ ) { 
    SURVEY_INFO.SURVEYFLAG[i] = 0 ; 
    NIDSURVEY[i]=0;  // to set DUMPFLAG
  }


  // --------------------------------------
  // check option to compute global biasCor for surveys lumped together
  if ( ( INPUTS.opt_biasCor & MASK_BIASCOR_SAMPLE ) == 0 )  {
    NSAMPLE_BIASCOR     = 1 ;
    ONE_SAMPLE_BIASCOR  = 1;
    IDSAMPLE            = 0 ;

    SAMPLE_BIASCOR[IDSAMPLE].NDATA = FITINP.NSNCUTS ;   
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME_FIELDGROUP,  "NONE" );
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME_SURVEYGROUP, "ALL"  );
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME,        "ALL" );
    SAMPLE_BIASCOR[IDSAMPLE].OPT_PHOTOZ   = 0 ; // zSPEC is default

    // set redshift info (Apr 4 2017)
    SAMPLE_BIASCOR[IDSAMPLE].zMIN = INPUTS.zmin ;
    SAMPLE_BIASCOR[IDSAMPLE].zMAX = INPUTS.zmax ;
    set_BINSIZE_SAMPLE_biasCor(IDSAMPLE);

    dump_SAMPLE_INFO("DATA");
    return ;  
  }


  // -----------------------------------
  set_FIELDGROUP_biasCor();  // optional
  set_SURVEYGROUP_biasCor(); // user-option, but always set internally

  // --------- start loop over data ---------------
  if ( NDMP ) 
    { printf(" xxx    NAME         FIELD        SURVEY  "
	     "IDSAMP  FIELDGROUP \n" );
    }

  sprintf(FIELD_TMP, "NONE"); 
  sprintf(FIELDGROUP,"NONE");  


  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {

    data[isn].idsample = -9 ;
    if ( data[isn].skipfit ) { continue ; }
    sprintf(SNID,"%s", data[isn].name );

    // strip off stuff into local var
    IDSURVEY   = data[isn].idsurvey ;
    OPT_PHOTOZ = data[isn].opt_photoz; 

    sprintf(SURVEYDEF,"%s", SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY] );
    sprintf(FIELDDEF, "%s", data[isn].field ); 
    NIDSURVEY[IDSURVEY]++ ;

    // check which FIELDGROUP group this event is in 
    if ( USE_FIELDGROUP )  { sprintf(FIELD_TMP,"%s",  FIELDDEF ); }

    // check if this IDSURVEY and FIELDGROUP have been defined.
    // Note that only one of the match_xxx functions can return STRINGOPT
    STRINGOPT[0] = 0 ;
    match_fieldGroup (SNID, FIELD_TMP,          // (I)
		      FIELDGROUP, STRINGOPT  ); // (O)

    if ( strlen(STRINGOPT) == 0 ) {
      match_surveyGroup(SNID, IDSURVEY,           // (I)
			SURVEYGROUP, STRINGOPT ); // (O)
    }	

    if ( OPT_PHOTOZ == 0 ) 
      { sprintf(zGROUP,"-zSPEC"); }
    else
      { sprintf(zGROUP,"-zPHOT"); }
    if ( FOUNDKEY_OPT_PHOTOZ == 0 ) { zGROUP[0] = 0 ; }

    DUMPFLAG = 0 ;
    IDSAMPLE = get_IDSAMPLE(IDSURVEY, OPT_PHOTOZ, 
			    FIELDGROUP, SURVEYGROUP ,DUMPFLAG);

    if ( IDSAMPLE < 0 ) {
      // store new SURVEY/FIELDGROUP entry
      N = NSAMPLE_BIASCOR ;  
      sprintf(SAMPLE_BIASCOR[N].NAME_FIELDGROUP,  "%s", FIELDGROUP  );
      sprintf(SAMPLE_BIASCOR[N].NAME_SURVEYGROUP, "%s", SURVEYGROUP ); 
      sprintf(SAMPLE_BIASCOR[N].STRINGOPT,        "%s", STRINGOPT   ); 
      SAMPLE_BIASCOR[N].OPT_PHOTOZ = OPT_PHOTOZ ;

      ptrName = SAMPLE_BIASCOR[N].NAME ;
      if ( IGNOREFILE(FIELDGROUP) ) { 
	sprintf(ptrName,"%s%s", SURVEYGROUP, zGROUP ); 
	SAMPLE_BIASCOR[N].IFLAG_ORIGIN = USERFLAG_SURVEYGROUP_SAMPLE ;
	// warning: this IFLAG_ORGIN does not distinguish USER and AUTO
      }
      else { 
	sprintf(ptrName,"%s%s(%s)", SURVEYGROUP, zGROUP, FIELDGROUP ); 
	SAMPLE_BIASCOR[N].IFLAG_ORIGIN = USERFLAG_FIELDGROUP_SAMPLE ;
      }

      NSAMPLE_BIASCOR++ ;

      N = NSAMPLE_BIASCOR ;
      if ( N >= MXNUM_SAMPLE ) {
	printf("\n PRE-ABORT DUMP: \n");
	dump_SAMPLE_INFO("DATA");
	sprintf(c1err,"NSAMPLE_BIASCOR=%d exceeds bound.", N);
	sprintf(c2err,"Check IDSURVEY and FIELDGROUP ");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }
   
    // execute mapindex function again to get final IDSAMPLE
    DUMPFLAG = NIDSURVEY[IDSURVEY] < 0 ; // for debug only
    IDSAMPLE = get_IDSAMPLE(IDSURVEY, OPT_PHOTOZ,
			    FIELDGROUP, SURVEYGROUP, DUMPFLAG );

    if ( IDSAMPLE < 0 ) {
      printf("\n PRE-ABORT DUMP in %s: \n", fnam);
      printf("   Current SURVEYGROUP = '%s' \n", SURVEYGROUP );
      printf("   Current FIELDGROUP  = '%s' \n", FIELDGROUP  );
      dump_SAMPLE_INFO("DATA");
      sprintf(c1err,"Final IDSAMPLE=%d for CID=%s(isn=%d)", 
	      IDSAMPLE, SNID, isn );
      sprintf(c2err,"SURVEY=%s(%d)  FIELD=%s", 
	      SURVEYDEF, IDSURVEY, FIELDDEF );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    data[isn].idsample = IDSAMPLE ; // store in global
    SAMPLE_BIASCOR[IDSAMPLE].NDATA++ ;

    // keep track of min/max redshift for each sample
    z = data[isn].zhd ;
    if( z < SAMPLE_BIASCOR[IDSAMPLE].zMIN ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMIN = z; }
    if( z > SAMPLE_BIASCOR[IDSAMPLE].zMAX ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMAX = z; }

    // check user option to skip biasCor for this sample
    if ( strstr(INPUTS.surveyList_noBiasCor,SURVEYDEF) != NULL ) 
      { SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR = 0 ; }

    // xxxxxxxxxxxxx
    if ( isn < NDMP ) {
      printf(" xxx '%8s'  '%8s'   %8.8s(%2d)   %2d    %s \n", 
	     SNID, FIELDDEF, SURVEYDEF,IDSURVEY,  IDSAMPLE, FIELDGROUP );
      fflush(stdout);
    }
    // xxxxxxxxxxx
    
  } // end isn loop

  // resort IDSAMPLE list to always have the same order
  // regardless of the data order.
  sort_IDSAMPLE_biasCor();

  // check for user-specified IDSAMPLEs to select ...AFTER sorting.
  parse_IDSAMPLE_SELECT(INPUTS.idsample_select); 

  // ------- prepare z-bins for each sample -------------
  // default for each sample is user-input z-bins (zmin,zmax,nzbin)
  int NSAMPLE = NSAMPLE_BIASCOR ;
  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE; IDSAMPLE++ ) 
    { set_BINSIZE_SAMPLE_biasCor(IDSAMPLE);  }
 

  // ---------------
  dump_SAMPLE_INFO("DATA");

  // check option to fix sigint for each IDSAMPLE
  parse_sigint_fix(INPUTS.sigint_fix);
   

  return ;

} // end prepare_IDSAMPLE_biasCor




// ======================================================
void  set_BINSIZE_SAMPLE_biasCor(int IDSAMPLE) {

  // Created Aug 2016
  // Decode STRINGOPT and set biasCor binSize for z,x1,c

  int ipar;
  char fnam[] = "set_BINSIZE_SAMPLE_biasCor" ;

  // ------------- BEGIN -----------
  
  // first set defaults
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[0] = INPUTS.zmin ;
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[1] = INPUTS.zmax ;
  SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT  = 
    (INPUTS.zmax - INPUTS.zmin) / (double)INPUTS.nzbin ;

  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_FITPAR[ipar]  
      = BIASCOR_BINSIZE_LCFIT[ipar] ;
  }


  // check STRINGOPT in () on surveygroup_biascor input.
  // STRINGOPT must be of the form:
  //   zbin=[value] 
  //   x1bin=[value] 
  //   cbin=[value] 
  // and multiple bin-specifiers can be separated by colons,
  //  STRINGOPT = 'zbin=[zbin]:x1bin=[x1bin]:cbin=[cbin]'
  //
  // Split STRING opt in two steps: first split by colon 
  // and then split each by equal sign.

#define MXBINSTRING 10
  char *ptr_binVar[MXBINSTRING], *ptr_binString[MXBINSTRING];
  char equal[] = "=", colon[] = ":"; 
  char *STRINGOPT = SAMPLE_BIASCOR[IDSAMPLE].STRINGOPT;
  int  NSPLIT_STRING, NSPLIT_TMP, isplit, LZBIN=0 ;
  int  MEMC = 40*sizeof(char);
  int  LDMP=0 ;
  double value ;     

  if ( LDMP ) {
    printf(" xxx ---------------------------------- \n");
    printf(" xxx %s : split STRINGOPT = '%s' for %s \n", 
	   fnam, STRINGOPT, SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }

  if ( strlen(STRINGOPT) == 0 ) { return ; }
    
  for(isplit=0; isplit < MXBINSTRING; isplit++ ) {
    ptr_binString[isplit] = (char*) malloc( MEMC ); // varName
    ptr_binVar[isplit]    = (char*) malloc( MEMC ); // varName
  }

  // split by colon
  splitString(STRINGOPT, colon, MXBINSTRING,            // inputs
	      &NSPLIT_STRING, ptr_binString );           // outputs
  
  // split each binString by equal sing
  for(isplit=0; isplit < NSPLIT_STRING; isplit++ ) {
    splitString(ptr_binString[isplit], equal, MXBINSTRING,   // inputs
		&NSPLIT_TMP, ptr_binVar );           // outputs
    if ( NSPLIT_TMP != 2 ) {
      sprintf(c1err,"Invalid NSPLIT_TMP=%d for '%s' (expect 2)", 
	      NSPLIT_TMP, ptr_binString[isplit] );
      sprintf(c2err,"check surveygroup_biascor key in input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;   
    }
    sscanf(ptr_binVar[1],"%le", &value); // strip off value

    // check which variable
    if ( strcmp(ptr_binVar[0],"zbin") == 0 ) 
      { SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT = value ; LZBIN=1; }
    else if ( strcmp(ptr_binVar[0],"x1bin") == 0 )
      { SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_FITPAR[INDEX_x1] = value ; }
    else if ( strcmp(ptr_binVar[0],"cbin") == 0 ) 
      { SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_FITPAR[INDEX_c]  = value ; }
    else {
      sprintf(c1err,"Invalid varName '%s' in '%s' ", 
	      ptr_binVar[0], STRINGOPT );
      sprintf(c2err,"check surveygroup_biascor key in input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;   
    }
    
    if ( LDMP ) 
      { printf(" xxx %s = %f \n", ptr_binVar[0], value); fflush(stdout); }


  } // end isplit

  // if zbin is specified, update the biasCor ZRANGE so that 
  // if zbin is small (for low-z), we do not have too many z-bins.
  if ( LZBIN ) {
    double zBIN = SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT ;
    double zMIN, zMAX, RANGE[2];
    int    NBz ;
    zMIN = SAMPLE_BIASCOR[IDSAMPLE].zMIN ; // true zMIN,zMAX
    zMAX = SAMPLE_BIASCOR[IDSAMPLE].zMAX ; // in data sample

    // round biasCor zmin to nearest 0.01
    RANGE[0] = .01 * (int)(100.*zMIN);
    NBz      = (int)((zMAX - RANGE[0])/zBIN) + 1 ;
    RANGE[1] = RANGE[0] + zBIN * (double)NBz ; 
    
    SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[0] = RANGE[0] ;
    SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[1] = RANGE[1] ;
  
    /*  
    printf(" zzz zRANGE(true)=%.3f,%.3f   "
	   "zRANGE(biasCor)=%.3f,%.3f  zBIN=%.2f  (NBZ=%d) \n",
	   zMIN, zMAX, RANGE[0], RANGE[1], zBIN, NBz );  */

  } // end LZBIN 

  // free local malloc
  for(isplit=0; isplit < MXBINSTRING; isplit++ ) {
    free( ptr_binString[isplit] );
    free( ptr_binVar[isplit] );
  }

  return ;

} // end set_BINSIZE_SAMPLE_biasCor
 

// =============================================
void set_FIELDGROUP_biasCor(void) {

  // extract FIELDGROUP sub strings from user-input keyed by
  // fieldgroup_biascor='xxx+yyy,zzz'
  
  int i, NGRP;
  int USE_FIELDGROUP  = INPUTS.use_fieldGroup_biasCor ;
  char *ptrFIELD[MXNUM_SAMPLE] ;
  //  char fnam[] = "set_FIELDGROUP_biasCor"; 

  // ---------- BEGIN ----------

  INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR = 0 ;

  if ( USE_FIELDGROUP == 0 ) { return ; }

  char comma[] = "," ;
  for(i=0; i < MXNUM_SAMPLE; i++ ) 
    { ptrFIELD[i] = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[i] ; }
  
  splitString(INPUTS.fieldGroup_biasCor, comma, MXNUM_SAMPLE, // inputs
	      &NGRP, ptrFIELD );   // outputs

  INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR = NGRP;
  for(i=0; i < NGRP; i++ ) {
    extractStringOpt(INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[i],
		     INPUTS_SAMPLE_BIASCOR.FIELDGROUP_OPTLIST[i] );
  }

  // To do:
  // Determine which survey(s) correspond to the fieldgroups, 
  // and mark them so that these surveys are excluded from SURVEYGROUP.
  //  SAMPLE_INFO.SURVEYFLAG[ID] = IFLAG_FIELDGROUP ;

  return ;

} // end set_FIELDGROUP_biasCor

// ==========================================
void  set_SURVEYGROUP_biasCor(void) {

  // Created Aug 22 2016
  // parse user-input surveygroup_biascor, and set SAMPLE_BIASCOR arrays.
  // If any surveys are missing from user input, append then as new
  // groups so that all surveys are included in a surveyGroup.
  //
  // May 14 2019: free ptrTmp[i]


  int  LDMP   = 0 ;
  int  USE_SURVEYGROUP = INPUTS.use_surveyGroup_biasCor ;
  int  i, i2 ;
  int  NGRP, ID ;
  char *ptrSURVEY[MXNUM_SAMPLE], *S ;
  char comma[] = "," ;
  char plus[]  = "+" ;
  char fnam[] = "set_SURVEYGROUP_biasCor" ;

  // ------------- BEGIN -------------

  // extract SURVEYGROUP sub strings (aug 18 2016)

  if ( USE_SURVEYGROUP ) {
    for(i=0; i < MXNUM_SAMPLE; i++ ) 
      { ptrSURVEY[i] = INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i] ; }

    splitString(INPUTS.surveyGroup_biasCor, comma, MXNUM_SAMPLE, // inputs
		&NGRP, ptrSURVEY ); // outputs
    INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR = NGRP ;
  }

  // convert each SURGROUP_LIST (string) into a list of IDSURVEYs

  // define temp string space for plus-separated strings.
  char *ptrTmp[MXNUM_SAMPLE] ;
  for(i=0; i < MXNUM_SAMPLE; i++ ) 
    { ptrTmp[i] = (char*) malloc ( MXCHAR_CCID * sizeof(char) ) ;  }


  NGRP = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR ;

  for(i=0; i < NGRP; i++ ) {

    // strip out OPTLIST option-string from () in SURVEYGROUP_LIST
    extractStringOpt(INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i],
		     INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_OPTLIST[i] ); 

    splitString(INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i], 
		plus, MXNUM_SAMPLE, // (I) 
		&INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[i], ptrTmp ); // (O)

    // store integer IDSURVEY for each plus-separated SURVEY 
    for(i2=0; i2 < INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[i]; i2++ ) {
      ID = get_IDSURVEY(ptrTmp[i2]);
      if( ID < 0 ) {
	sprintf(c1err,"Undefined SURVEY = '%s'", ptrTmp[i2] );
	sprintf(c2err,"check SURVEYGROUP = '%s' ", 
		INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      if ( LDMP ) {
	printf(" xxx %s: IGRP=%d i2=%d  IDSRVY=%3d  SURVEY=%s \n",
	       fnam, i, i2, ID, ptrTmp[i2] ); fflush(stdout);
      }

      INPUTS_SAMPLE_BIASCOR.IDSURVEYGROUP_LIST[i][i2] = ID ;
      SURVEY_INFO.SURVEYFLAG[ID] = USERFLAG_SURVEYGROUP_SAMPLE ;
    }

  } // end i loop over user-defined survey groups

  // -------------------------------------------
  // for each UN-USED survey, add a new SURVEYGROUP so that
  // every survey is in a group. Make sure to ignore survey 
  // that is already part of a FIELDGROUP.

  INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT = 
    INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR ;

  for(ID=0; ID < MXIDSURVEY; ID++ ) {
    S = SURVEY_INFO.SURVEYDEF_LIST[ID] ;

    if ( SNDATA_INFO.NEVT_PER_SURVEY[ID] == 0 ) 
      { continue;  }
    if ( SURVEY_INFO.SURVEYFLAG[ID] == USERFLAG_SURVEYGROUP_SAMPLE ) 
      { continue ; }
    if ( SURVEY_INFO.SURVEYFLAG[ID] == AUTOFLAG_SURVEYGROUP_SAMPLE ) 
      { continue ; }
     
    // add another survey group with only this one survey
    NGRP = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT ;
    if ( NGRP < MXNUM_SAMPLE ) {
      INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[NGRP] = 1 ;
      INPUTS_SAMPLE_BIASCOR.IDSURVEYGROUP_LIST[NGRP][0] = ID ;
      sprintf(INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[NGRP],"%s", S);

      SURVEY_INFO.SURVEYFLAG[ID] = AUTOFLAG_SURVEYGROUP_SAMPLE ; 
      if ( NGRP==0 ) 
	{ sprintf(INPUTS.surveyGroup_biasCor, "%s",S); }
      else  { 
	strcat(INPUTS.surveyGroup_biasCor,",") ;
	strcat(INPUTS.surveyGroup_biasCor,S) ;
      }
    }
    INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT++ ;

  } // end loop over IDSURVEY

  NGRP = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT ;
  if ( NGRP >= MXNUM_SAMPLE ) {
    sprintf(c1err,"NSURVEYGRPOUP=%d exceeds bound of MXNUM_SAMPLE=%d",
	    NGRP, MXNUM_SAMPLE );
    sprintf(c2err,"check samples." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  for(i=0; i < MXNUM_SAMPLE; i++ )     { free(ptrTmp[i]); }

  return ;

} // end set_SURVEYGROUP_biasCor


// ======================================================
void sort_IDSAMPLE_biasCor(void) {

  // Aug 26 2016
  // Sort IDSAMPLE list in the following order:
  // * user-defined surveygroup_biascor
  // * user-defined fieldgroup_biascor
  // * auto-defined surveys groups, in same order as in SURVEY.DEF
  //
  //  This sorting ensures that IDSAMPLE has the same meaning
  //  for any random simulation and any random subset of data.
  //  Reason is that order of data determines order of IDSAMPLE,
  //  and it is thus random.
  //

  int NSAMPLE = NSAMPLE_BIASCOR ;
  int ID, IDSURVEY, NSORT=0  ;
  int NCOPY[MXNUM_SAMPLE], IDMAP_SORT[MXNUM_SAMPLE];
  int igrp, NGRP_USR,  FLAG ; 
  char *s1, *s2, *NAME ;
  SAMPLE_INFO_DEF *SAMPLE_BIASCOR_TEMP ;
  int LDMP = 0 ;
  char fnam[] = "sort_IDSAMPLE_biasCor" ;
  // ---------------- BEGIN ---------------

  SAMPLE_BIASCOR_TEMP = 
    (SAMPLE_INFO_DEF*) malloc ( NSAMPLE * sizeof(SAMPLE_INFO_DEF));
  
  // copy pre-sorted struct into temp struct
  for(ID=0; ID < NSAMPLE; ID++ ) {
    copy_IDSAMPLE_biasCor(&SAMPLE_BIASCOR[ID], &SAMPLE_BIASCOR_TEMP[ID]);
    NCOPY[ID] = 0; 
  }


  if ( LDMP ) {
    printf("\n xxx ---------- Start DUMP for %s ---------------- \n", fnam);
    dump_SAMPLE_INFO("DATA");
  }

  // start with user-defined survey groups, in the order
  // given by the user
  NGRP_USR = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR ;
  for(igrp=0; igrp < NGRP_USR; igrp++ ) {
    s1 = INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[igrp];
    for(ID=0; ID < NSAMPLE; ID++ ) {
      s2   = SAMPLE_BIASCOR_TEMP[ID].NAME_SURVEYGROUP ;
      NAME = SAMPLE_BIASCOR_TEMP[ID].NAME ;
      if ( strcmp(s1,s2)==0 ) {
	copy_IDSAMPLE_biasCor(&SAMPLE_BIASCOR_TEMP[ID], 
			      &SAMPLE_BIASCOR[NSORT]);
	
	if ( LDMP ) {
	  printf(" 0. xxx copy USER-SURVEYGROUP ID=%d --> %d (%s)\n", 
		 ID, NSORT, NAME );
	}

	IDMAP_SORT[ID]=NSORT;  NSORT++ ;  NCOPY[ID]++ ;  
      }
    }
  }  // end i loop over user-define survey groups


  // check user-defined field groups
  NGRP_USR = INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR ;
  for(igrp=0; igrp < NGRP_USR; igrp++ ) {
    s1 = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[igrp];
    for(ID=0; ID < NSAMPLE; ID++ ) {
      s2   = SAMPLE_BIASCOR_TEMP[ID].NAME_FIELDGROUP ;
      NAME = SAMPLE_BIASCOR_TEMP[ID].NAME ;
      if ( strcmp(s1,s2)==0 ) {
	copy_IDSAMPLE_biasCor(&SAMPLE_BIASCOR_TEMP[ID], 
			      &SAMPLE_BIASCOR[NSORT]);

	if ( LDMP ) {
	  printf(" 1. xxx copy USER-FIELDGROUP  ID=%d --> %d (%s) \n", 
		 ID, NSORT, NAME);
	}
	IDMAP_SORT[ID]=NSORT;  NSORT++ ;  NCOPY[ID]++ ;
      }
    }
  }  // end i loop over user-define field groups


  // add auto-generated survey groups, 
  //  in order given in SURVEY.DEF file


  for(IDSURVEY=0; IDSURVEY < MXIDSURVEY; IDSURVEY++ ) {

    // check only those surveys which are auto-generated
    FLAG = SURVEY_INFO.SURVEYFLAG[IDSURVEY] ;
    if ( FLAG != AUTOFLAG_SURVEYGROUP_SAMPLE ) { continue ; }

    s1 = SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY];
    for(ID=0; ID < NSAMPLE; ID++ ) {

      // skip surveys that are part of FIELDGROUP
      FLAG = SAMPLE_BIASCOR_TEMP[ID].IFLAG_ORIGIN ;
      if ( FLAG == USERFLAG_FIELDGROUP_SAMPLE)	{ continue ; }

      s2   = SAMPLE_BIASCOR_TEMP[ID].NAME_SURVEYGROUP ;
      NAME = SAMPLE_BIASCOR_TEMP[ID].NAME ;
      if ( strcmp(s1,s2) == 0 ) {
	copy_IDSAMPLE_biasCor(&SAMPLE_BIASCOR_TEMP[ID], 
			      &SAMPLE_BIASCOR[NSORT]);
	if ( LDMP ) {
	  printf(" 2. xxx copy AUTO-SURVEYGROUP ID=%d --> %d (%s) \n", 
		 ID, NSORT, NAME ); fflush(stdout);
	}
	IDMAP_SORT[ID]=NSORT;  NSORT++ ;  NCOPY[ID]++ ;
      }
    } // end IDSAMPLE loop

  }   // end IDSURVEY


  // ------- check that each sample was copied once and only once.
  for(ID=0; ID < NSAMPLE; ID++ ) {
    if( NCOPY[ID] != 1 ) {
      sprintf(c1err,"IDSAMPLE=%d(%s) copied %d times.",
	      ID, SAMPLE_BIASCOR_TEMP[ID].NAME, NCOPY[ID] );
      sprintf(c2err,"but should be copied once (NSAMPLE=%d, NSORT=%d)",
	      NSAMPLE, NSORT );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }


  // ---------------------------------------------
  // remap all of the IDSAMPLE tags for the data
  int isn  ;
  for(isn=0; isn < FITINP.NSNCUTS; isn++ ) {
    ID = data[isn].idsample;
    if ( ID < 0 ) { continue ; }
    data[isn].idsample = IDMAP_SORT[ID]; // sorted IDSAMPLE
  }


  if ( LDMP  ) 
    { printf(" xxx ----- END DUMP for %s --------- \n\n", fnam); }

  free(SAMPLE_BIASCOR_TEMP);
  return ;

} // end sort_IDSAMPLE_biasCor


// ======================================================
void copy_IDSAMPLE_biasCor(SAMPLE_INFO_DEF *S1, SAMPLE_INFO_DEF *S2) {

  // Copy contents of struct S1 into struct S2.

  int  ipar ;
  //  char fnam[] = "copy_IDSAMPLE_biasCor" ;

  // -------------- BEGIN -------------
  
  sprintf(S2->NAME_FIELDGROUP,  "%s", S1->NAME_FIELDGROUP ) ;
  sprintf(S2->NAME_SURVEYGROUP, "%s", S1->NAME_SURVEYGROUP ) ;
  sprintf(S2->STRINGOPT,        "%s", S1->STRINGOPT ) ;
  sprintf(S2->NAME,             "%s", S1->NAME ) ;
  S2->OPT_PHOTOZ = S1->OPT_PHOTOZ ;

  S2->NDATA    = S1->NDATA ;
  S2->NBIASCOR = S1->NBIASCOR ;
  S2->NCCPRIOR = S1->NCCPRIOR ;

  S2->zMIN = S1->zMIN ;
  S2->zMAX = S1->zMAX ;

  S2->RANGE_REDSHIFT[0] = S1->RANGE_REDSHIFT[0] ;
  S2->RANGE_REDSHIFT[1] = S1->RANGE_REDSHIFT[1] ;
  S2->BINSIZE_REDSHIFT  = S1->BINSIZE_REDSHIFT ;

  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    { S2->BINSIZE_FITPAR[ipar] = S1->BINSIZE_FITPAR[ipar] ; }

  //  S2->DOFLAG_SELECT  = S1->DOFLAG_SELECT  ;
  S2->DOFLAG_BIASCOR = S1->DOFLAG_BIASCOR ;
  S2->IFLAG_ORIGIN   = S1->IFLAG_ORIGIN ;

  return ;

} // end copy_IDSAMPLE_biasCor


// ==========================================
void dump_SAMPLE_INFO(char *what) {

  // Dump SAMPLE_INFO structure.
  // Input *what = "DATA" or "BIASCOR"  or  "CCPRIOR"
  //   --> which NEVT to print.
  //
  // Mar 22 2018: 
  //  + for BIASCOR or CCPRIOR, abort if any sample has zero event.
  //  + allow some CCPRIOR samples to have NEVT=0 (e.g., spec subsets)
  //
  int  isdata, NOBIASCOR, NEVT, i, N = NSAMPLE_BIASCOR;
  int  NZERR=0, NSAMPLE_ZERO=0, ABORT_ON_NEVTZERO=0;
  char *NAME, Nstring[60], zString[60] ;
  char fnam[] = "dump_SAMPLE_INFO" ;

  // -------- BEGIN ----------

  printf("  SAMPLE_INFO DUMP for %s: \n", what )    ;
  isdata=0;

  for(i=0; i < N; i++ )  { 
    if ( strcmp(what,"DATA") == 0 ) 
      { NEVT = SAMPLE_BIASCOR[i].NDATA ; isdata=1; }
    else if ( strcmp(what,"BIASCOR") ==  0 ) 
      { NEVT = SAMPLE_BIASCOR[i].NBIASCOR ; ABORT_ON_NEVTZERO=1; }
    else
      { NEVT = SAMPLE_BIASCOR[i].NCCPRIOR ; ABORT_ON_NEVTZERO=0; }

    NAME      = SAMPLE_BIASCOR[i].NAME ;
    NOBIASCOR = ( SAMPLE_BIASCOR[i].DOFLAG_BIASCOR == 0 );
    if ( NEVT == 0 ) { NSAMPLE_ZERO++; }

    sprintf(Nstring,"N%-7.7s=%6d", what, NEVT);
    if ( NOBIASCOR && isdata==0 ) { sprintf(Nstring,"No %s", what); }

    sprintf(zString,"%.3f<z<%.3f", 
	    SAMPLE_BIASCOR[i].zMIN, SAMPLE_BIASCOR[i].zMAX );

    if( SAMPLE_BIASCOR[i].zMAX < SAMPLE_BIASCOR[i].zMIN ) { NZERR++; }
    // print one-line summary per SAMPLE
    printf("  IDSAMPLE=%2d --> %-20.20s  (%s, %s)\n",
	   i, NAME, Nstring, zString );
  }

  if ( NZERR > 0 ) {
    sprintf(c1err,"zMAX < zMIN for %d IDSAMPLEs", NZERR);
    sprintf(c2err,"See IDSAMPLE dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( ABORT_ON_NEVTZERO && NSAMPLE_ZERO ) {
    sprintf(c1err,"%d of %d SAMPLES have zero events.", NSAMPLE_ZERO,N);
    sprintf(c2err,"Check %s sample.", what);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  fflush(stdout);

  return ;

} // end dump_SAMPLE_INFO


void match_fieldGroup(char *SNID, char *FIELD, 
		       char *FIELDGROUP, char *STRINGOPT ) {

  // for input *FIELD, return *FIELDGROUP for which *FIELD matches.
  // *SNID is the SN name, and is used only for error messages.

  int  NMATCH_TOT, NMATCH_FIELD, igroup, NFIELD_OVP, ifield ;
  int  NFIELDGROUP = INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR ;
  char *FTMP1, *FTMP2, fieldList[MXFIELD_OVERLAP][MXCHAR_CCID] ;
  char fnam[] = "match_fieldGroup" ;

  // ----------- BEGIN ------------

  NMATCH_TOT = 0 ;

  sprintf(FIELDGROUP, "NONE" );
  STRINGOPT[0] = 0 ;
  if ( strcmp(FIELD,"NONE") == 0 ) { return ; }

  if ( ONE_SAMPLE_BIASCOR ) { return ; }

  // get list of all overlapping fields (separated by '+')
  if ( strstr(FIELD,"+") != NULL ) {
    char *ptrFlist[MXFIELD_OVERLAP];
    for(ifield=0; ifield<MXFIELD_OVERLAP; ifield++ ) 
      { ptrFlist[ifield] = fieldList[ifield] ;  } 
    splitString(FIELD, "+", MXFIELD_OVERLAP, &NFIELD_OVP, ptrFlist ); 
  }
  else {
    // just one field; no overlaps
    NFIELD_OVP = 1 ;
    sprintf(fieldList[0],"%s", FIELD);
  }

  // if(NFIELD_OVP>1) {printf(" xxx ---------------- \n"); fflush(stdout);}

  for(igroup=0; igroup < NFIELDGROUP; igroup++ ) {
    FTMP1 = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[igroup] ;

    NMATCH_FIELD=0;
    // check if any overlap field matches FIELDGROUP
    for(ifield = 0; ifield < NFIELD_OVP; ifield++ ) {
            
      FTMP2 = fieldList[ifield];
      if ( strstr(FTMP2,FTMP1)!=NULL || strstr(FTMP1,FTMP2)!=NULL ) { 

	// for field-overlap, first match has priority.
	if( NMATCH_TOT==0 ) { 
	  sprintf(FIELDGROUP,"%s", FTMP1);  
	  sprintf(STRINGOPT,"%s",
		  INPUTS_SAMPLE_BIASCOR.FIELDGROUP_OPTLIST[igroup]);
	}

	NMATCH_TOT++ ;
	NMATCH_FIELD++;   
      }

    } // end ifield
  } // end igroup
  

  // abort on ambiguous match
  if ( NMATCH_TOT > 1  && NFIELD_OVP==1 ) {
    printf("\n PRE-ABORT DUMP \n");
    printf("  fieldgroup_biasCor = '%s' \n", INPUTS.fieldGroup_biasCor);
    for(igroup=0; igroup < NFIELDGROUP ; igroup++ ) {
      FTMP1 = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[igroup] ;
      printf("  FIELDGROUP(%d) = '%s' \n", igroup, FTMP1);
    }
    sprintf(c1err,"Invalid NMATCH_TOT=%d for CID=%s",   NMATCH_TOT, SNID );
    sprintf(c2err,"FIELD='%s' and FIELDGROUP list", FIELD);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;

} // end match_fieldGroup


void match_surveyGroup(char *SNID, int IDSURVEY, 
		       char *SURVEYGROUP, char *STRINGOPT ) {

  // Created Aug 19 2016
  // for input IDSURVEY, return *SURVEYGROUP for which IDSURVEY matches.
  // *SNID is the SN name, and is used only for error messages.

  int  NMATCH_TOT, igroup, i2, IDSURVEY_TMP ;
  int  NSURVEYGROUP = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT ;
  char *STMP ;
  char fnam[] = "match_surveyGroup" ;

  // ----------- BEGIN ------------

  NMATCH_TOT = 0 ;
  sprintf(SURVEYGROUP, "NONE" );
  STRINGOPT[0] = 0 ;

  if ( ONE_SAMPLE_BIASCOR ) { return ; }

  for(igroup=0; igroup < NSURVEYGROUP; igroup++ ) {

    for(i2=0; i2 < INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[igroup]; i2++ ) {
      IDSURVEY_TMP = INPUTS_SAMPLE_BIASCOR.IDSURVEYGROUP_LIST[igroup][i2] ;
      if ( IDSURVEY == IDSURVEY_TMP ) {
	NMATCH_TOT++ ; 
	sprintf(SURVEYGROUP, "%s", 
		INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[igroup] ) ;
	sprintf(STRINGOPT,   "%s", 
		INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_OPTLIST[igroup] );
      }
    }
  } // end igroup

  
  // abort on ambiguous match
  if ( NMATCH_TOT != 1 ) {
    printf("\n PRE-ABORT DUMP \n");
    printf("  surveygroup_biasCor = '%s' \n", INPUTS.surveyGroup_biasCor);
    for(igroup=0; igroup < NSURVEYGROUP ; igroup++ ) {
      STMP = INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[igroup] ;
      printf("  SURVEYGROUP(%d) = '%s' \n", igroup, STMP);
    }
    sprintf(c1err,"Invalid NMATCH_TOT=%d for CID=%s",   NMATCH_TOT, SNID );
    sprintf(c2err,"IDSURVEY=%d(%s)", 
	    IDSURVEY, SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  return ;

} // end match_surveyGroup


int get_IDSAMPLE(int IDSURVEY, int OPT_PHOTOZ, 
		 char *FIELDGROUP, char *SURVEYGROUP,
		 int DUMPFLAG ) {

  // Created July 2016
  // Utility to return IDSAMPLE for inputs
  //  IDSURVEY   : survey ID
  //  OPT_PHOTOZ : 0=>zSPEC  >0 ==> zPHOT  
  //  FIELDGROUP 
  //  SURVEYGROUP
  //
  // Used to label IDSAMPLE for data, biasCor and CCprior.

  int N = NSAMPLE_BIASCOR ;
  int  MATCH_FIELDGROUP, MATCH_SURVEYGROUP, MATCH_PHOTOZ, MATCH_OR, i ;
  int  CHECK_FIELDGROUP  = ( IGNOREFILE(FIELDGROUP)  == 0 ) ;
  int  CHECK_SURVEYGROUP = ( IGNOREFILE(SURVEYGROUP) == 0 ) ;
  int  OPT_PHOTOZ_TMP ;
  char *FIELDGROUP_TMP, *SURVEYGROUP_TMP ;
  char *SURVEYDEF = SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY] ;
  char fnam[] = "get_IDSAMPLE" ;

  // ----------- BEGIN ------------

  // check option to lump all surveys & fields together
  if ( ONE_SAMPLE_BIASCOR ) { return 0; }

  /* xxx maybe add this back later if FIELDGROUP survey is identified
  // abort on invalid option to check both FIELDGROUP & SURVEYGROUP
  if ( CHECK_FIELDGROUP && CHECK_SURVEYGROUP ) { 
    printf("\n PRE-ABORT DUMP: \n");
    printf("  Input IDSURVEY    = %d (%s) \n", IDSURVEY, SURVEYDEF );
    printf("  Input FIELDGROUP  = %s \n", FIELDGROUP );
    printf("  Input SURVEYGROUP = %s \n", SURVEYGROUP );
    sprintf(c1err,"Invalid check on FIELDGROUP & SURVEYGROUP");
    sprintf(c2err,"Check surveygroup_biascor and fieldgroup_biascor keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  xxxx */

  if ( DUMPFLAG ) { 
    printf(" xxx ----------- %s DUMP ---------------- \n", fnam );
    printf(" xxx Inputs: IDSURVEY=%d  FIELDGROUP=%s  SURVEYGROUP=%s\n",
	   IDSURVEY, FIELDGROUP, SURVEYGROUP ); fflush(stdout);
  }

  // -----------------------------
  // check for survey/fieldGroup
  for(i=0; i < N; i++ ) {

    FIELDGROUP_TMP    = SAMPLE_BIASCOR[i].NAME_FIELDGROUP ;
    SURVEYGROUP_TMP   = SAMPLE_BIASCOR[i].NAME_SURVEYGROUP ;
    OPT_PHOTOZ_TMP    = SAMPLE_BIASCOR[i].OPT_PHOTOZ ;

    // check matches
    MATCH_FIELDGROUP  = 0;
    MATCH_SURVEYGROUP = 0 ; 

    if ( CHECK_FIELDGROUP ) 
      { MATCH_FIELDGROUP  = (strcmp(FIELDGROUP_TMP,FIELDGROUP) == 0); }
    else if ( CHECK_SURVEYGROUP )
      { MATCH_SURVEYGROUP = (strcmp(SURVEYGROUP_TMP,SURVEYGROUP) == 0); }
   
    MATCH_PHOTOZ = (OPT_PHOTOZ_TMP == OPT_PHOTOZ ) ;


    // abort if both match
    if ( MATCH_FIELDGROUP && MATCH_SURVEYGROUP ) { 
      printf("\n PRE-ABORT DUMP: \n");
      printf("  Input IDSURVEY    = %d (%s) \n",  IDSURVEY, SURVEYDEF );
      printf("  Input FIELDGROUP  = %s \n",  FIELDGROUP );
      printf("  Input SURVEYGROUP = %s \n",  SURVEYGROUP );
      printf("  Matches FIELDGROUP[%d]  = %s \n",  i, FIELDGROUP_TMP );
      printf("  Matches SURVEYGROUP[%d] = %s \n",  i, SURVEYGROUP_TMP );
      sprintf(c1err,"Invalid match to both FIELD & SURVEY group");
      sprintf(c2err,"Check surveygroup_biascor and fieldgroup_biascor keys");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( DUMPFLAG ) { 
      printf(" xxx INFO(%d): GROUP_TMP[F,S]=%s,%s  MATCH[F,S]=%d,%d\n", 
	     i, FIELDGROUP_TMP, SURVEYGROUP_TMP, 
	     MATCH_FIELDGROUP, MATCH_SURVEYGROUP ) ;
      fflush(stdout);
    }

    MATCH_OR = ( MATCH_FIELDGROUP || MATCH_SURVEYGROUP ) ;
    if ( MATCH_OR && MATCH_PHOTOZ ) { return i; }
  }

  return(-9);

} // end get_IDSAMPLE


// =========================================================
void prepare_biasCor(char *simFile) {

  // Created Jan 2016
  // If simFile is not null, then prepare mB bias correction vs. z
  // for SNIa.
  //  -> read simFile
  //  -> make map of mB-SIM_mB vs redshift
  //
  // Jun 2 2016: call sigInt_biasCor BEFORE makeMape_fitPar
  // Jan 16 2018: force correct value for INPUTS.fitflag_sigmb = 2 

  int INDX, IDSAMPLE, SKIP ;
  int  OPTMASK = INPUTS.opt_biasCor ;
  int  DOCOR_1DZAVG  = ( OPTMASK & MASK_BIASCOR_1DZAVG  );
  int  DOCOR_1DZWGT  = ( OPTMASK & MASK_BIASCOR_1DZWGT  );
  int  DOCOR_1D5DCUT = ( OPTMASK & MASK_BIASCOR_1D5DCUT );
  int  DOCOR_5D      = ( OPTMASK & MASK_BIASCOR_5D      );
  int  DOCOR_MUCOV   = ( OPTMASK & MASK_BIASCOR_MUCOV   );
  int  IDEAL         = ( OPTMASK & MASK_BIASCOR_COVINT ) ;
  char txt_biasCor[40];

  char fnam[] = "prepare_biasCor" ;

  // ------------- BEGIN -------------

  simdata_bias.DOCOR_1D = 0 ;
  simdata_bias.DOCOR_5D = DOCOR_5D ;
  
  simdata_bias.NROW = 0 ; // number of simulated entries
  simdata_bias.NUSE = 0 ;
  SKIPZCUT_BIASCOR  = 0 ; // Jan 25, 2018

  if ( IGNOREFILE(simFile) ) { return ; }


  if ( DOCOR_5D ) {  
    if ( DOCOR_MUCOV ) 
      { sprintf(txt_biasCor,"5D+MUCOV"); }
    else
      { sprintf(txt_biasCor,"5D"); }
  }
  else if ( DOCOR_1DZAVG )   
    {  sprintf(txt_biasCor,"mu-z (WGT=1)");  simdata_bias.DOCOR_1D=1; }
  else if ( DOCOR_1DZWGT ) 
    {  sprintf(txt_biasCor,"mu-z (WGT=1/COV)"); simdata_bias.DOCOR_1D=1; }
  else if ( DOCOR_1D5DCUT ) 
    {  sprintf(txt_biasCor,"mu-z (WGT=1/COV,5DCUT)"); simdata_bias.DOCOR_1D=1; }
  else {
    sprintf(c1err,"Invalid opt_biascor=%d", OPTMASK );
    sprintf(c2err,"grep MASK_BIASCOR  SALT2mu.c | grep define ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  
  sprintf(BANNER,"%s: read sim to interpolate %s biasCor", 
	  fnam, txt_biasCor);
  print_banner(BANNER);
  printf("\t\t (opt_biascor=%d, sigma_cell=%.2f)\n", 
	 OPTMASK, INPUTS.sigma_cell_biasCor);
  fflush(stdout);

  // -------------------------------------------
  int PS0 = INPUTS.prescale_biasCor[0];
  int PS1 = INPUTS.prescale_biasCor[1];
  if ( PS1 > 1 ) {
    printf("\t Apply pre-scale: select subset %d of %d \n",
	   PS0, PS1);  
  }


  // Jan 16 2018: force correct option for log(sigma) term in chi2
  if ( DOCOR_5D ) 
    { INPUTS.fitflag_sigmb = 2 ; } // include log(sigma) term
  else
    { INPUTS.fitflag_sigmb = 1 ; } // ignore log(sigma) term
  printf("\t Force fitflag_sigmb = %d \n", INPUTS.fitflag_sigmb);


  // check which cov parameter in fit
  if ( IDEAL )  { 
    printf("\t Compute COVINT(z) using IDEAL fit params.\n");
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ;
  }
    
  // --------------------------------
  read_simFile_biasCor(simFile);  

  if ( DOCOR_1DZAVG || DOCOR_1DZWGT ) { goto CHECK_1DCOR ; }
  
  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ ) 
    { setup_CELLINFO_biasCor(IDSAMPLE); }

  // store ia and ib for every biasCor event
  store_iaib_biasCor();

  // make sparse list for faster looping below (Dec 21 2017)
  makeSparseList_biasCor();

  // count number of biasCor events passing cuts (after setup_BININFO calls)
  int ievt,  NCUTS=0;
  int NROW = simdata_bias.NROW;
  for(ievt=0; ievt < NROW; ievt++ ) { if(biasMapSelect(ievt)){NCUTS++;} }
  printf("\t %d of %d biasCor events pass cuts. \n",  NCUTS, NROW);

  if ( NCUTS == 0 ) {
    dumpStages_biasMapSelect();
    sprintf(c1err,"No BiasCor events pass cuts." ); 
    sprintf(c2err,"Check cut stages above, and check BiasCor file." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  // determine sigInt for biasCor sample BEFORE makeMap since
  // sigInt is needed for 1/muerr^2 weight
  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ ) 
    {  init_sigInt_biasCor_legacy(IDSAMPLE);  }
  
  // determine intrinsic scatter matrix on grid of
  // idsample, redshift, alpha, beta
  init_COVINT_biasCor();

  // get wgted-average z in each user redshift bin for M0-vs-z output
  calc_zM0_biasCor();

  // ----------------------------------------------------------------
  // -------- START LOOP OVER SURVEY/FIELDGROUP SUB-SAMPLES ---------
  // ----------------------------------------------------------------

  char *NAME_SAMPLE;
  char borderLine[] = 
    "@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@" ;


  for( IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ ) {

    printf("\n %s\n", borderLine);   

    // print ISAMPLE info to screen
    NAME_SAMPLE  = SAMPLE_BIASCOR[IDSAMPLE].NAME ; 

    SKIP = 0 ;
    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_SELECT  == 0 ) { SKIP=1; }
    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { SKIP=1; }

    if ( SKIP ) {
      printf("\t SKIP BIASCOR PREP for %s\n", NAME_SAMPLE);
      continue ;
    }
    else  {
      printf("\t START BIASCOR PREP for SAMPLE = %s\n", NAME_SAMPLE);
    }
    fflush(stdout);

    // get wgted avg in each bin to use for interpolation 
    makeMap_binavg_biasCor(IDSAMPLE);

    // prepare 3D bias maps need to interpolate bias
    for(INDX=0; INDX<NLCPAR; INDX++ ) 
      { makeMap_fitPar_biasCor(IDSAMPLE,INDX); }
    printf("\n");  fflush(stdout);
    
    // make map of sigma_mu bias
    if ( DOCOR_MUCOV ) { makeMap_sigmu_biasCor(IDSAMPLE); }

    printf("\n\t END BIASCOR PREP for %s\n", NAME_SAMPLE);
    printf(" %s\n", borderLine);       fflush(stdout);
  }


  // -------------------------------------------------------------
  // -------- END LOOP OVER SURVEY/FIELDGROUP SUB-SAMPLES --------
  // -------------------------------------------------------------

  // compute and store bias(mB,x1,c) for each data event.
  // If data event lies in undefined biasBin(z,x1,c) then reject event.
  int n, istore ;
  int NSKIP_TOT=0, NUSE_TOT=0;
  int NSKIP[MXNUM_SAMPLE] ;
  int NUSE[MXNUM_SAMPLE];

  for(IDSAMPLE=0; IDSAMPLE<MXNUM_SAMPLE; IDSAMPLE++ ) {
    NSKIP[IDSAMPLE] = 0 ;
    NUSE[IDSAMPLE]  = 0 ;
  }

  printf("\n");
  printf(" For each DATA event before fit, store \n");
  printf("   * mb,x1,c-bias at each alpha,beta \n");
  if ( DOCOR_MUCOV) 
    { printf("   * muCOVscale   at each alpha,beta \n"); }

  for (n=0; n < FITINP.NSNCUTS; ++n) {
    if ( data[n].skipfit ) { continue ; }
    istore = storeDataBias(n,0);  

    IDSAMPLE = data[n].idsample ;
    NUSE[IDSAMPLE]++ ; NUSE_TOT++ ;
    if ( istore == 0 )  { 
      NSKIP_TOT++; NSKIP[IDSAMPLE]++ ; data[n].skipfit = 1; 
      if ( INPUTS.dumpflag_nobiasCor )
	{ storeDataBias(n,1); } // pass dump flag
    }    
  }

  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR; IDSAMPLE++ ) {
    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_SELECT == 0 ) { continue; }
    printf("  Rejected %d (of %d) DATA events with no 5D biasCor (%s). \n",
	   NSKIP[IDSAMPLE], NUSE[IDSAMPLE], SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }
  printf("\n");
  fflush(stdout);

  // ------------
  // free memory used to hold each simBias event;
  // we only need to keep the small map made in makeMap_fitPar_biasCor.

  if ( NSKIP_TOT == NUSE_TOT ) {
    sprintf(c1err,"Could not compute biasCor");
    sprintf(c2err,"Something is messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // check option for JLA-like correction
 CHECK_1DCOR:
  if ( simdata_bias.DOCOR_1D ) {
    INPUTS.opt_biasCor |= MASK_BIASCOR_1DZ ; // set OR-bit (AVG or WGT)
    prepare_biasCor_zinterp();  
  }

  malloc_simdata_biasCor(-1,0);
  return ;

} // end prepare_biasCor


// =============================================
void  read_simFile_biasCor(char *simFile) {

  // Created Jan 2016
  // Read simFile and load simdata_bias struct.
  // use SNTABLE_XXX utilities.
  //
  // May 11 2017: make sure IDSAMPLE is loaded before applying cuts.
  // Sep 04 2017: read IDEAL params
  // Sep 11 2017: read SIM_ZCMB
  // Jun 02 2018: fix bug by requiring true Ia in loop over events.
  // May 14,2019: if OPT_PHOTOZ is not read, init opt_photoz[i]=0

  int NROW, IFILETYPE, NVAR_ORIG, icut, LEN;
  int ivar, ivar2, ivar_opt_photoz ;
  int vbose = 1 ;  // verbose, but no abort on missing variable
  int ABORT = 3 ;
  int USE_FIELDGROUP  = INPUTS.use_fieldGroup_biasCor ;
  int IDEAL = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  int FOUND_VPEC=0, FOUND_SIM_VPEC=0;
  double z, zerr, zpecerr;
  char str_z[40], str_zerr[40], vartmp[60], STRINGOPT[40];
  char fnam[] = "read_simFile_biasCor" ;

  // --------------- BEGIN ---------------

  t_read_biasCor[0] = time(NULL); 

  LEN = SNTABLE_NEVT(simFile,TABLENAME_FITRES);  LEN += 10;
  simdata_bias.LEN_MALLOC = LEN ;

  printf("\t Found %d lines in %s. \n", LEN, simFile);
  malloc_simdata_biasCor(+1,LEN);

  IFILETYPE = TABLEFILE_OPEN(simFile,"read");
  NVAR_ORIG = SNTABLE_READPREP(IFILETYPE,"FITRES");

  sprintf(vartmp, "CID:C*%d  CCID:C*%d", MXCHAR_CCID, MXCHAR_CCID); 
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.name, LEN, vbose) ;

  if ( USE_FIELDGROUP ) {
    sprintf(vartmp, "FIELD:C*%d", MXCHAR_CCID ); 
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.field, LEN, vbose) ;
  }

  str_z[0] = str_zerr[0] = 0 ;
  get_zString(str_z,str_zerr,"F");
  SNTABLE_READPREP_VARDEF(str_z,   simdata_bias.z,    LEN ,vbose);
  SNTABLE_READPREP_VARDEF(str_zerr,simdata_bias.zerr, LEN ,vbose);

  sprintf(vartmp,"VPEC:F" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.vpec, LEN ,vbose);
  sprintf(vartmp,"VPEC_ERR:F VPECERR:F" );
  ivar2 = SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.vpecerr, LEN ,vbose);
  if ( ivar>0 && ivar2 > 0 ) { FOUND_VPEC=1; }

  sprintf(vartmp,"SIM_NONIA_INDEX:I SIM_TEMPLATE_INDEX:I" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_NONIA_INDEX, LEN ,vbose);

  sprintf(vartmp,"IDSURVEY:I" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.idsurvey, LEN ,vbose);

  // - - - - - 
  sprintf(vartmp,"OPT_PHOTOZ:I" );
  ivar_opt_photoz =
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.opt_photoz, LEN ,vbose);
  if ( FOUNDKEY_OPT_PHOTOZ && ivar_opt_photoz < 0 ) {
    sprintf(c1err,"OPT_PHOTOZ column in data, but not in biasCor.");
    sprintf(c2err,"OPT_PHOTOZ must be in both data and biasCor, or neither.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  } 
  if ( FOUNDKEY_OPT_PHOTOZ==0 && ivar_opt_photoz >= 0 ) {
    sprintf(c1err,"OPT_PHOTOZ column in biasCor, but not in data.");
    sprintf(c2err,"OPT_PHOTOZ must be in both data and biasCor, or neither.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  // - - - - - 

  sprintf(vartmp,"SNRMAX1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SNRMAX, LEN,vbose);

  // - - - - - mB - - - - - -
  sprintf(vartmp,"mB:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL[INDEX_mB], LEN,vbose);
  sprintf(vartmp,"mBERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITERR[INDEX_mB], LEN,vbose);
  sprintf(vartmp,"SIM_mb:F SIM_mB:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIMVAL[INDEX_mB], LEN,vbose);

  sprintf(vartmp,"x0:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.x0, LEN,vbose);
  sprintf(vartmp,"x0ERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.x0err, LEN,vbose);

  // - - - - - x1 - - - - - -
  sprintf(vartmp,"x1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL[INDEX_x1], LEN,vbose);
  sprintf(vartmp,"x1ERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITERR[INDEX_x1], LEN,vbose);
  sprintf(vartmp,"SIMx1:F SIM_x1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIMVAL[INDEX_x1], LEN,vbose);

  // - - - - - c - - - - - -
  sprintf(vartmp,"c:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL[INDEX_c], LEN,vbose);
  sprintf(vartmp,"cERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITERR[INDEX_c], LEN,vbose);
  sprintf(vartmp,"SIMc:F SIM_c:F" );
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIMVAL[INDEX_c], LEN,vbose);
 
  // - - - - covariances - -- - - - 
  sprintf(vartmp,"COVx0x1:F COV_x1_x0:F");
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.COV_x0x1, LEN,vbose);

  sprintf(vartmp,"COVx0c:F COV_c_x0:F");
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.COV_x0c, LEN,vbose);

  sprintf(vartmp,"COVx1c:F COV_x1_c:F");
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.COV_x1c, LEN,vbose);

  // - - - -  Generated alpha & beta - - - - -
  sprintf(vartmp,"SIMalpha:F SIM_alpha:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_ALPHA, LEN,vbose);

  sprintf(vartmp,"SIMbeta:F  SIM_beta:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_BETA,  LEN,vbose);

  // - - - - - - - - - - - - - - - - 
  // IDEAL fit params from fit where FLUX = TRUEFLUX(Sep 2017) 
  if ( IDEAL ) {
    sprintf(vartmp,"x0_IDEAL:F" ) ;
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.x0_IDEAL,  
			    LEN, ABORT);

    sprintf(vartmp,"mB_IDEAL:F" );
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL_IDEAL[INDEX_mB], 
			    LEN, ABORT);
    sprintf(vartmp,"x1_IDEAL:F" );
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL_IDEAL[INDEX_x1], 
			    LEN, ABORT);
    sprintf(vartmp,"c_IDEAL:F" );
    SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.FITVAL_IDEAL[INDEX_c], 
			    LEN, ABORT);
  }

  // - - - - - - - - - - - - - - - - - - - 

  // true z & MU
  sprintf(vartmp,"SIM_ZCMB:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_ZCMB,  LEN,vbose);

  sprintf(vartmp,"SIM_DLMAG:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_DLMAG,  LEN,vbose);

  // true VPEC (Aug 2017)
  sprintf(vartmp,"SIM_VPEC:F" ) ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp, simdata_bias.SIM_VPEC, LEN,vbose);
  if ( ivar > 0 ) { FOUND_SIM_VPEC = 1; }

  //read CUTWIN variables to apply same cuts as for data.
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
    sprintf(vartmp, "%s:F", INPUTS.CUTWIN_NAME[icut]);
    if ( usesim_CUTWIN(vartmp) == 0 ) { continue ; }
    ivar = SNTABLE_READPREP_VARDEF(vartmp,simdata_bias.CUTVAL[icut],LEN,vbose);
    simdata_bias.DOFLAG_CUTWIN[icut] = set_DOFLAG_CUTWIN(ivar,icut,0);
  }

  // ==========================================
  // ==========================================
  // read entire file and load arrays
  NROW = SNTABLE_READ_EXEC();
  simdata_bias.NROW = NROW ;
  // ==========================================
  // ==========================================


  // close file
  // xxx obsolete  TABLEFILE_CLOSE(simFile); 

  t_read_biasCor[1] = time(NULL); 

  // ----------------------------
  // prepare loop over bias cor to prepare lots of stuff:
  // + define mB = -2.5*log(x0), to be like the data.
  // + Set IDSAMPLE 
  // + set COV_FITVAL matrix
  // + set min/max for alpha and beta
  
  int i ;
  double mB, x0, mBoff_fit ; 
  x0 = (double)simdata_bias.x0[0] ;
  mB = (double)simdata_bias.FITVAL[INDEX_mB][0] ;
  mBoff_fit = -2.5*log10(x0) - mB ;

  printf("\t mB(sim) --> -2.5*log10(x0) = mB + %.4f \n", mBoff_fit);
  fflush(stdout);

  int   IDSURVEY, IDSAMPLE, OPT_PHOTOZ;
  double CUTVAL_LOCAL[MXCUTWIN], a, amin, amax, b, bmin, bmax;
  char *NAME, FIELD[MXCHAR_CCID], FIELDGROUP[100] ;
  char SURVEY[MXCHAR_CCID], SURVEYGROUP[100];
  sprintf(FIELD,"NONE");  sprintf(SURVEY,"NONE");

  amin = bmin =  99999.  ; // init min alpha,beta
  amax = bmax = -99999.  ; // init max alpha,beta

  for (i=0; i < simdata_bias.NROW ; i++ ) {

    if ( simdata_bias.SIM_NONIA_INDEX[i] != 0 ) { continue ; } 

    simdata_bias.FITVAL[INDEX_mB][i] += (float)mBoff_fit ;
    simdata_bias.SIMVAL[INDEX_mB][i] += (float)mBoff_fit ;
    
    // if OPT_PHOTOZ is not defined, set opt_photoz to zero since
    // it was not read..
    if ( ivar_opt_photoz < 0 ) { simdata_bias.opt_photoz[i]; }

    // convert SURVEY/FIELD into mapID
    IDSURVEY   = simdata_bias.idsurvey[i] ;
    OPT_PHOTOZ = simdata_bias.opt_photoz[i] ;
    NAME     = simdata_bias.name[i] ;
    if ( USE_FIELDGROUP ) 
      { sprintf(FIELD,"%s", simdata_bias.field[i] ); }

    match_fieldGroup (NAME, FIELD,    
		      FIELDGROUP, STRINGOPT  ); // (O)
    match_surveyGroup(NAME, IDSURVEY, 
		      SURVEYGROUP, STRINGOPT ); // (O)
    
    IDSAMPLE = get_IDSAMPLE(IDSURVEY,OPT_PHOTOZ,FIELDGROUP,SURVEYGROUP,0) ;
    simdata_bias.idsample[i] = IDSAMPLE ;

    // May 11 2017: apply cuts AFTER IDSAMPLE is loaded.
    // Apply user cuts for stats in dump_SAMPLE_INFO 
    // Do NOT use biasMapSelect() since it's not ready yet.
    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
      { CUTVAL_LOCAL[icut] = (double)simdata_bias.CUTVAL[icut][i] ; }
    if ( apply_CUTWIN(2,simdata_bias.DOFLAG_CUTWIN,CUTVAL_LOCAL) == 0 )  
      { continue  ; }

    SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR++ ;

    // update min,max alpha & beta
    a = (double)simdata_bias.SIM_ALPHA[i];
    b = (double)simdata_bias.SIM_BETA[i];
    if ( a < amin ) { amin=a; }
    if ( a > amax ) { amax=a; }
    if ( b < bmin ) { bmin=b; }
    if ( b > bmax ) { bmax=b; }

    // set covariance matrix
    make_COV_FITVAL_biasCor(i);

    // Jan 8 2018: if vpec variables are not there, set relevant 
    // values to zero to avoid crazy values.
    if ( FOUND_VPEC == 0 ) 
      { simdata_bias.vpec[i] = simdata_bias.vpecerr[i] = 0.0 ; }
    if ( FOUND_SIM_VPEC == 0 ) 
      { simdata_bias.SIM_VPEC[i] = 0.0 ; }

    // Jan 9 2018: back-compatibility code for versions before 10_56.
    // If vpecerr=0 in the biasCor table, then apply INPUTS.zpecerr.
    if ( simdata_bias.vpecerr[i] < 1.0E-9 ) {
      z       = simdata_bias.z[i] ;
      zerr    = simdata_bias.zerr[i] ;
      zpecerr = INPUTS.zpecerr * ( 1.0 + z ) ;
      simdata_bias.zerr[i] = sqrt(zerr*zerr + zpecerr*zpecerr);
    }

  } // end i-loop over events


  INPUTS_SAMPLE_BIASCOR.ALPHA_MIN = amin;
  INPUTS_SAMPLE_BIASCOR.ALPHA_MAX = amax ;
  INPUTS_SAMPLE_BIASCOR.BETA_MIN  = bmin;
  INPUTS_SAMPLE_BIASCOR.BETA_MAX  = bmax ;

  // print number of biasCor events vs. SURVEY/FIELD
  dump_SAMPLE_INFO("BIASCOR") ;

  t_read_biasCor[2] = time(NULL); 

  return ;
} // end read_simFile_biasCor


// ================================================================
void  make_COV_FITVAL_biasCor(int ievt) {

  // Created Sep 4 2017
  // Load covariance matrix for biasCor event 'ievt'
  // and store in COV_FITVAL[j0][j1]. This matrix
  // is based only on LC fits, and does not include
  // intrinsic scatter.

  
  double mBerr    = (double)simdata_bias.FITERR[INDEX_mB][ievt] ;
  double x1err    = (double)simdata_bias.FITERR[INDEX_x1][ievt] ;
  double cerr     = (double)simdata_bias.FITERR[INDEX_c][ievt] ;
  double COV_x0x1 = (double)simdata_bias.COV_x0x1[ievt] ;
  double COV_x0c  = (double)simdata_bias.COV_x0c[ievt] ;
  double COV_x1c  = (double)simdata_bias.COV_x1c[ievt] ;
  double x0       = (double)simdata_bias.x0[ievt] ;
  double x0err    = (double)simdata_bias.x0err[ievt] ;
  double sf       = -2.5/(x0*LOGTEN) ;  
  double COV00, COV01, COV02, COV12, COV11, COV22 ;
  char   *name    = simdata_bias.name[ievt];

  int IDEAL = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;

  int j0, j1;
  char fnam[] = "make_COV_FITVAL_biasCor" ;

  // ------------- begin ------------

  // build cov matrix 
  COV00 = (sf*x0err) * (sf*x0err) ;
  COV01 = sf * COV_x0x1 ;
  COV02 = sf * COV_x0c  ;
  COV12 = COV_x1c ;
  COV11 = x1err * x1err;
  COV22 = cerr * cerr ;

  // store COV in global struct
  simdata_bias.COV_FITVAL[INDEX_mB][INDEX_mB][ievt] = (float)COV00 ;
  simdata_bias.COV_FITVAL[INDEX_x1][INDEX_x1][ievt] = (float)COV11 ;  
  simdata_bias.COV_FITVAL[INDEX_c ][INDEX_c ][ievt] = (float)COV22 ;
  simdata_bias.COV_FITVAL[INDEX_mB][INDEX_x1][ievt] = (float)COV01 ;
  simdata_bias.COV_FITVAL[INDEX_mB][INDEX_c ][ievt] = (float)COV02 ;
  simdata_bias.COV_FITVAL[INDEX_x1][INDEX_c ][ievt] = (float)COV12 ;
  simdata_bias.COV_FITVAL[INDEX_x1][INDEX_mB][ievt] = (float)COV01 ;
  simdata_bias.COV_FITVAL[INDEX_c ][INDEX_mB][ievt] = (float)COV02 ;
  simdata_bias.COV_FITVAL[INDEX_c ][INDEX_x1][ievt] = (float)COV12 ;
  

  return ;

} // end make_COV_FITVAL_biasCor

// ================================================================
void  prepare_biasCor_zinterp(void) {

  // Jun 28 2016 
  // Special biasCor option to correct mb(z) using muBias(sample,z).
  // Evaluate mubias(sample,z) from biasCor sample.
  //
  // Sep 18 2016: fix to include idsample-dependence
  // Sep 17 2017: include lensing error
  // Sep 25 2017: if zDATA is outside biasCor range, reject data
  //              event instead of aborting.
  //

  int    NSAMPLE  = NSAMPLE_BIASCOR ;
  int    NROW     = simdata_bias.NROW ;
  double alpha    = INPUTS.parval[IPAR_ALPHA0] ;
  double beta     = INPUTS.parval[IPAR_BETA0] ;

  int    ievt, NCUTS, iz, izusr, idsample ;
  double mu_fit, mB_fit, x1_fit, c_fit ;
  double mu_sim, mB_sim, x1_sim, c_sim ;
  double mu_err, mB_err, x1_err, c_err ;
  double mu_true, WGT, z, zMIN[MXNUM_SAMPLE], zMAX[MXNUM_SAMPLE] ;
  double zlo[MXNUM_SAMPLE][MAXBIN_z], zhi[MXNUM_SAMPLE][MAXBIN_z];
  double zavg[MXNUM_SAMPLE][MAXBIN_z], zbinSize[MXNUM_SAMPLE], zbin;
  double USR_SUMz[MAXBIN_z], USR_SUMWGT[MAXBIN_z];
  double MUBIAS_FIT[MXNUM_SAMPLE][MAXBIN_z];
  double MUBIAS_SIM[MXNUM_SAMPLE][MAXBIN_z];
  double SUMWGT[MXNUM_SAMPLE][MAXBIN_z];
  double SUMWGT_HISNR, MUBIAS_SIM_HISNR ;
  double muerrsq, muerrLens, SNRMAX ;

  int ndata[MXNUM_SAMPLE], ndata_reject[MXNUM_SAMPLE];
  int NDATA_REJECT=0;
  char fnam[] = "prepare_biasCor_zinterp" ;
  
  // ---------------- BEGIN --------------

 
  simdata_bias.SIGINT_AVG = .11 ; // temp hard-wire
    
  printf("\n");
  printf("   Compute muBias vs. z with: \n" );
  printf("\t    alpha=%.3f(p1) beta=%.3f(p2) \n", alpha, beta);
  //  printf("\t    zbinSize = %.3f \n", zbinSize );
  fflush(stdout);

    
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    zMIN[idsample]     = SAMPLE_BIASCOR[idsample].RANGE_REDSHIFT[0] ;
    zMAX[idsample]     = SAMPLE_BIASCOR[idsample].RANGE_REDSHIFT[1] ;
    zbinSize[idsample] = SAMPLE_BIASCOR[idsample].BINSIZE_REDSHIFT ;
    ndata[idsample]    = ndata_reject[idsample] = 0 ;

    /*
    printf(" xxx Sample=%s -> zMIN=%.3f  zMAX=%.3f  zbin=%.3f \n",
	   SAMPLE_BIASCOR[idsample].NAME,
	   zMIN[idsample], zMAX[idsample], zbinSize[idsample] );
	   fflush(stdout);  */
  }
  
  for(iz=0; iz < MAXBIN_z; iz++ ) {

    for(idsample=0; idsample < NSAMPLE; idsample++ ) {
      MUBIAS_FIT[idsample][iz] = MUBIAS_SIM[idsample][iz] = 0.0;
      SUMWGT[idsample][iz]     = 0.0 ;

      zbin               = zbinSize[idsample] ;
      zlo[idsample][iz]  = zMIN[idsample]    + zbin * (double)iz ;
      zhi[idsample][iz]  = zlo[idsample][iz] + zbin ;
      zavg[idsample][iz] = 0.5 * ( zlo[idsample][iz] + zhi[idsample][iz] ) ;
    } 

    MUBIAS_SIM_HISNR = SUMWGT_HISNR = 0.0 ;
    USR_SUMz[iz] = USR_SUMWGT[iz] = 0.0 ;
  }

  // compute muBias in z-bins
  NCUTS = 0 ;

  for(ievt=0; ievt < NROW; ievt++ ) { 
    z         = (double)simdata_bias.z[ievt];
    idsample  = simdata_bias.idsample[ievt];

    if ( biasMapSelect(ievt) == 0 ) { continue ; }    
    if ( z < zMIN[idsample] ) { continue ; }
    if ( z > zMAX[idsample] ) { continue ; }

    NCUTS++ ; 

    mB_fit    = (double)simdata_bias.FITVAL[INDEX_mB][ievt] ;
    x1_fit    = (double)simdata_bias.FITVAL[INDEX_x1][ievt] ;
    c_fit     = (double)simdata_bias.FITVAL[INDEX_c ][ievt] ;

    mB_err    = (double)simdata_bias.FITERR[INDEX_mB][ievt] ;
    x1_err    = (double)simdata_bias.FITERR[INDEX_x1][ievt] ;
    c_err     = (double)simdata_bias.FITERR[INDEX_c ][ievt] ;
    
    mB_sim    = (double)simdata_bias.SIMVAL[INDEX_mB][ievt] ;
    x1_sim    = (double)simdata_bias.SIMVAL[INDEX_x1][ievt] ;
    c_sim     = (double)simdata_bias.SIMVAL[INDEX_c ][ievt] ;

    SNRMAX    = (double)simdata_bias.SNRMAX[ievt];

    mu_true = (double)simdata_bias.SIM_DLMAG[ievt]; 
    mu_fit  = mB_fit + alpha*x1_fit - beta*c_fit - M0_DEFAULT ;
    mu_sim  = mB_sim + alpha*x1_sim - beta*c_sim - M0_DEFAULT ;
    mu_err  = -9.0; // not used

    muerrLens = ( z * INPUTS.lensing_zpar ) ;

    iz = (int)( (z-zMIN[idsample]) / zbinSize[idsample] ) ;
    
    if ( iz >= MAXBIN_z ) {
      sprintf(c1err,"iz=%d exceeds bound for z=%.4f", iz, z);
      sprintf(c2err,"Check MAXBIN_z in SALT2mu.c");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // get approx error using only diagonal covariance

    muerrsq = 
      (INPUTS.sigmB * INPUTS.sigmB) +
      (mB_err * mB_err) +
      (x1_err * x1_err) * (alpha * alpha) +
      (c_err  * c_err ) * (beta  * beta ) +
      (muerrLens * muerrLens)       // 9.17.2017
      ;

    if ( INPUTS.opt_biasCor & MASK_BIASCOR_1DZAVG ) 
      { WGT = 1.0 ; }
    else 
      { WGT = 1.0/muerrsq ; }

    SUMWGT[idsample][iz]     += WGT ;
    MUBIAS_FIT[idsample][iz] += WGT*(mu_fit - mu_true) ;
    MUBIAS_SIM[idsample][iz] += WGT*(mu_sim - mu_true) ;

    if ( SNRMAX > INPUTS.snrmin_sigint_biasCor  ) {
      SUMWGT_HISNR     += WGT ;
      MUBIAS_SIM_HISNR += (mu_sim - mu_true) ;
    }

    // get izusr = iz bin for user z bins ... can be different
    // from the z-bins used to interpolate the muBias
    // Needed to write z-wtged average in the M0-vs-z file.
    // Note that these sums are over the entire biasCor sample,
    // regardless of 'idsample'.

    izusr = IBINFUN(z, &INPUTS.BININFO_z, 0, "");
    USR_SUMz[izusr]   += (WGT*z) ;
    USR_SUMWGT[izusr] +=  WGT ;

  } // end loop over biasCor events

  printf("\t %d of %d biasCor events pass cuts. \n",  NCUTS, NROW);
  simdata_bias.NUSE = NCUTS;

  if ( NCUTS < 10 ) {
    dumpStages_biasMapSelect();
    sprintf(c1err,"NCUTS=%d is too few.", NCUTS );
    sprintf(c2err,"Something is wrong.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // --------------------------------
  // get MU bias in each z-bin
  // Also store IZMIN, IZMAX and NBZ for interpolation below
  int izmin, izmax;
  int IZMIN[MXNUM_SAMPLE];
  int IZMAX[MXNUM_SAMPLE];
  int NBZ[MXNUM_SAMPLE];

  MUBIAS_SIM_HISNR /= SUMWGT_HISNR ;
  printf("\t muBias(sim,SNR>%.0f) = %6.3f \n",
	 INPUTS.snrmin_sigint_biasCor, MUBIAS_SIM_HISNR);

  for(idsample=0; idsample < NSAMPLE; idsample++ ) {

    if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { continue; }
    IZMIN[idsample] = +999;
    IZMAX[idsample] = -999;
    
    printf("     %s : \n", SAMPLE_BIASCOR[idsample].NAME );
    for(iz=0; iz < MAXBIN_z; iz++ ) {

      if ( SUMWGT[idsample][iz] <= 1.0E-9 ) { continue ; }

      if ( iz < IZMIN[idsample] ) { IZMIN[idsample] = iz; }
      if ( iz > IZMAX[idsample] ) { IZMAX[idsample] = iz; }
      
      MUBIAS_FIT[idsample][iz] /= SUMWGT[idsample][iz] ;
      MUBIAS_SIM[idsample][iz] /= SUMWGT[idsample][iz] ;

      // subtract muBias for hi-SNR so that muBias ~ 0 at low-z.
      MUBIAS_FIT[idsample][iz] -= MUBIAS_SIM_HISNR ;
      MUBIAS_SIM[idsample][iz] -= MUBIAS_SIM_HISNR ;

      printf("\t z=%.3f - %.3f : muBias = %6.3f  (iz=%d)\n",
	     zlo[idsample][iz], zhi[idsample][iz],
	     MUBIAS_FIT[idsample][iz], iz );
      
    }  // end iz
    fflush(stdout);

  } // end idsample

  for(idsample=0; idsample < NSAMPLE; idsample++ ) 
    {  NBZ[idsample] = IZMAX[idsample] - IZMIN[idsample] + 1; }


  // -----------------------------------
  // get z-wgted aveage in user z-bins 
  for(izusr=0; izusr < INPUTS.nzbin; izusr++ ) {
    if ( USR_SUMWGT[izusr] > 0.0 )
      { simdata_bias.zM0[izusr] = USR_SUMz[izusr]/USR_SUMWGT[izusr];}
    else
      { simdata_bias.zM0[izusr] = -9.0 ; }
  }

  // -----------------------------------------------
  printf("\n   Now compute muBias(z) for each data event. \n");
  fflush(stdout);

  double mB, muBias ;

  for(ievt=0; ievt < FITINP.NSNCUTS; ievt++ ) {

    idsample = data[ievt].idsample ;
    if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) 
      { data[ievt].skipfit=1; }

    data[ievt].muBias_zinterp = 0.0 ;
    if ( data[ievt].skipfit ) { continue ; }

    z        = data[ievt].zhd ;
    mB       = data[ievt].fitpar[INDEX_mB];
    izmin    = IZMIN[idsample];
    izmax    = IZMAX[idsample];
    
    ndata[idsample]++ ;
    if ( z < zlo[idsample][izmin] || z > zhi[idsample][izmax] ) {
      data[ievt].skipfit = 1; 
      ndata_reject[idsample]++ ;
      NDATA_REJECT++ ;
    }
    if ( data[ievt].skipfit ) { continue ; }

    if ( z < zavg[idsample][izmin] ) 
      { muBias = MUBIAS_FIT[idsample][izmin]; }
    else if ( z > zavg[idsample][izmax] ) 
      { muBias = MUBIAS_FIT[idsample][izmax]; }
    else {
      // interpolate
      int OPT_INTERP = 1; 
      muBias = interp_1DFUN(OPT_INTERP, z, NBZ[idsample], 
			    &zavg[idsample][izmin], 
			    &MUBIAS_FIT[idsample][izmin], fnam);
    } 

    data[ievt].muBias_zinterp = muBias ;

    // xxxxxxxxxxxxxxxxxxxx
    if ( ievt < -20 ) {
      printf(" xxx ievt=%3d : z=%.3f --> muBias=%6.3f \n",
	     ievt, z, muBias);
      fflush(stdout);
    }
    // xxxxxxxxxxxxxxxxxxxx

  } // end ievt loop over data


  // give warning if any data are rejected because z is outside biasCor range
  int n, nrej;  float ratio ;
  if ( NDATA_REJECT > 0  ){
    printf("\n");
    printf(" ):--):--):--):--):--):--):--):--):--):--):--):--):--):--\n") ;
    printf(" %s WARNING: %d data events with invalid z rejected:\n",
	   fnam, NDATA_REJECT );
    for(idsample=0; idsample < NSAMPLE; idsample++ ) {
      n = ndata[idsample] ;  nrej = ndata_reject[idsample] ;
      if ( nrej > 0 ) {
	ratio = 100.0 * (float)nrej / (float)n;
	printf("\t NREJECT(%s) = %d (%.2f) \n",
	       SAMPLE_BIASCOR[idsample].NAME, nrej, ratio );  
      }
    }
    printf(" ):--):--):--):--):--):--):--):--):--):--):--):--):--):--\n" ) ;
    fflush(stdout);
  }

    
  return ;

} // end prepare_biasCor_zinterp

// ================================================================
void set_MAPCELL_biasCor(int IDSAMPLE) {

  
  // set  CELLINFO_BIASCOR.MAPCELL[ia][ib][iz][ix1][ic]
  // to map 5D indices to 1D index J1D.

  int ia, ib, iz, ix1, ic, NCELL, ipar;
  char fnam[] = "set_MAPCELL_biasCor" ;

  // -------------- BEGIN ----------------------------
  int NBINz   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINx1  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;
  int NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;

  // ------------------------------------------
  // establish map between 1D array and 5D binning (only do once)
  
  NCELL = 0 ;
  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(iz=0; iz<NBINz; iz++ ) {
	for(ix1=0; ix1<NBINx1; ix1++ ) {
	  for(ic=0; ic<NBINc; ic++ ) {
	    if ( NCELL < MAXBIN_BIASCOR_1D ) 
	      { CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] 
		  = NCELL ; }
	    NCELL++ ;
	  }
	}
      }
    }
  }

  // check array bound.
  if ( NCELL >= MAXBIN_BIASCOR_1D ) {
    sprintf(c1err,"NCELL=%d > bound of MAXBIN_BIASCOR_1D=%d", 
	    NCELL, MAXBIN_BIASCOR_1D ) ;
    sprintf(c2err,"NBIN(a,b,z,x1,c) = %d %d %d %d %d",
	   NBINa, NBINb, NBINz, NBINx1, NBINc);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  CELLINFO_BIASCOR[IDSAMPLE].NCELL = NCELL ;

  printf("  %s : malloc for %d  5D-cells for biasCor \n", fnam, NCELL );
  fflush(stdout);

  // malloc other CELLINFO arrays
  int MEMD0   = NCELL   * sizeof(double);
  int MEMI0   = NCELL   * sizeof(int);

  CELLINFO_BIASCOR[IDSAMPLE].NperCell = (int   *) malloc(MEMI0);
  CELLINFO_BIASCOR[IDSAMPLE].AVG_z    = (double*) malloc(MEMD0);
  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    {  CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar] = (double*)malloc(MEMD0);}

  // -------------------------------------
  // malloc FITPARBIAS & MUCOVSCALE based on number of 5D cells
  // and number of SURVEY-FIELDs.

  int MEMBIAS = NCELL   * sizeof(FITPARBIAS_DEF ) ;
  int MEMCOV  = NCELL   * sizeof(double  ) ;

  simdata_bias.FITPARBIAS[IDSAMPLE] = (FITPARBIAS_DEF*) malloc ( MEMBIAS);
  simdata_bias.MUCOVSCALE[IDSAMPLE] = (double *       ) malloc ( MEMCOV );

  fflush(stdout);
  return ;

} // end set_MAPCELL_biasCor


// ================================================================
void makeMap_fitPar_biasCor(int IDSAMPLE, int ipar_LCFIT) {

  // Created Jan 2016 by R.Kessler
  // After reading simFile_bias, prepare bias vs. redshift & parval
  // so that fitting model can be corrected. Fill struct simdata_bias.
  // ipar_LCFIT is for mB, x1 or c.
  //
  // WARNING: computes bias with user z-bining, so beware of
  //          combined low-z plus hi-z samples.
  //
  //
  // - - - - - - - - - -

#define MXJ1DNBR 10000

  int NCELL = CELLINFO_BIASCOR[IDSAMPLE].NCELL ;
  int NBIASCOR_CUTS  = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;

  int MEMD, NEVT_USE, NEVT_SAMPLE,  LDMP, N ;
  int J1D, ia, ib, iz, ix1, ic, ievt, isp ;
  int J1DNBR_LIST[MXJ1DNBR], NJ1DNBR ;

  double fitval, simval, biasVal ;
  double WGT, VAL, ERR, RMS, SQRMS, XN, XNLIST, tmp1, tmp2 ;
  double *SUMBIAS, *SUMWGT, *sum, *sumsq ;
  char *PARNAME = BIASCOR_NAME_LCFIT[ipar_LCFIT] ;
   
  char fnam[] = "makeMap_fitPar_biasCor";

  // ----------------- BEGIN -----------------

  printf("  %s of %s-bias(z,x1,c,a,b)  \n", fnam, PARNAME ); 
  fflush(stdout);

  // malloc arrays to store info in each biasCor cell
  MEMD    = NCELL * sizeof(double) ;
  SUMBIAS = (double*) malloc(MEMD) ;
  SUMWGT  = (double*) malloc(MEMD) ;
  sum     = (double*) malloc(MEMD) ;
  sumsq   = (double*) malloc(MEMD) ;
  
  // ------------------------------------------
  for( J1D=0; J1D < NCELL ; J1D++ ) {
    // init global arrays  
    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT]   = 999. ;
    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT]   = 999. ;
    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].RMS[ipar_LCFIT]   = 999. ;

    CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D]  = 0 ;

    // init local arrays
    SUMBIAS[J1D] = SUMWGT[J1D] = 0.0 ;
    sum[J1D]     = sumsq[J1D]  = 0.0 ; 
  }


  // -----------------------------------------------
  // -------- LOOP OVER BIASCOR SIM ROWS -----------
  // -----------------------------------------------

  NEVT_SAMPLE = NEVT_USE = 0 ;

  for(isp=0; isp < NBIASCOR_CUTS; isp++ ) {

    ievt = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp];
    LDMP = ( ievt < -5 );
    NEVT_SAMPLE++ ;

    // get bias for ipar_LCFIT = mB,x1 or c
    fitval  = (double)simdata_bias.FITVAL[ipar_LCFIT][ievt] ;
    simval  = (double)simdata_bias.SIMVAL[ipar_LCFIT][ievt] ;
    biasVal = fitval - simval ;

    WGT = WGT_biasCor(2,ievt,fnam) ;  // WGT=1/muerr^2 for wgted average
    J1D = J1D_biasCor(ievt,fnam);     // 1D index

    SUMBIAS[J1D]  += (WGT * biasVal) ;
    SUMWGT[J1D]   += WGT ;

    // increment unweighted sums to get RMS and bias-error
    sum[J1D]   += biasVal ;
    sumsq[J1D] += (biasVal * biasVal) ;  

    // -----------------
    NEVT_USE++ ;
    CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D]++ ;
  }

  // convert sums in each IZ,LCFIT bin into mean bias,
  
  int NMAP_TOT = 0 ;
  int NMAP_USE = 0 ;
  double sumsq_nbr, sum_nbr ;
  int    Nsum_nbr, J1D_nbr, inbr ;
  int    NgetList = 0 ;

  for(J1D=0; J1D < NCELL; J1D++ ) {

    NMAP_TOT++ ;  
    N   = CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D] ;
    XN  = (double)N ;

    if ( N < 1 ) { continue ; }
	
    NMAP_USE++ ;
	
    WGT = SUMWGT[J1D] ;               // WGT
    VAL = SUMBIAS[J1D]/SUMWGT[J1D] ;  // wgted bias value

    // if too few events in cell, sum 3x3 nbr grid to get
    // better stats for RMS
    NJ1DNBR = 1;  J1DNBR_LIST[0] = J1D ;
    if ( N < MINPERCELL_BIASCOR ) { 
      get_J1DNBR_LIST(IDSAMPLE, J1D, &NJ1DNBR, J1DNBR_LIST) ; 
      if ( NJ1DNBR >= MXJ1DNBR ) {
	sprintf(c1err,"NJ1DNBR=%d exceeds bound of %d",
		NJ1DNBR, MXJ1DNBR);
	sprintf(c2err,"J1D=%d", J1D);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }
      NgetList++ ;
    }

    sumsq_nbr = sum_nbr = 0.0 ;  Nsum_nbr=0;
    for(inbr=0; inbr < NJ1DNBR; inbr++ ) {
      J1D_nbr  = J1DNBR_LIST[inbr] ;
      
      if ( J1D_nbr < 0 || J1D_nbr >= NCELL ) {
	sprintf(c1err,"Invalid J1D_nbr=%d for inbr=%d", J1D_nbr, inbr);
	sprintf(c2err,"J1D=%d", J1D);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }
      sumsq_nbr += sumsq[J1D_nbr] ;
      sum_nbr   += sum[J1D_nbr] ;
      Nsum_nbr  += CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D_nbr] ;
    }
    XNLIST = (double)Nsum_nbr ;

    tmp1   = sum_nbr/XNLIST;   tmp2=sumsq_nbr/XNLIST ;
    SQRMS  = tmp2 - tmp1*tmp1 ;

    RMS    = sqrt( fabs(SQRMS) ); // protect against -epsilon (Nov 18 2016)
    ERR    = RMS/sqrt(XN);   // note denom uses XN, not XNLIST
    if ( ERR == 0.0 ) { ERR = 1.0E-12; } // avoid nan = 1/ERR^2

    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT] = VAL ;
    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT] = ERR ;
    simdata_bias.FITPARBIAS[IDSAMPLE][J1D].RMS[ipar_LCFIT] = RMS ;
    
  } // end J1D loop


  // -----------------------------------------------
  // print grid-cell stats on last parameter
  // (since it's the same for each parameter)
  if ( ipar_LCFIT == INDEX_c ) {
    printf("  BiasCor computed for %d of %d grid-cells with >=1 events.\n",
	   NMAP_USE, NMAP_TOT ) ;
    printf("  BiasCor sample: %d of %d pass cuts for IDSAMPLE=%d.\n",
	   NEVT_USE, NEVT_SAMPLE, IDSAMPLE );

    if ( NEVT_USE == 0 ) {
      dumpStages_biasMapSelect();
      sprintf(c1err,"No BiasCor events passed for %s", 
	      SAMPLE_BIASCOR[IDSAMPLE].NAME );
      sprintf(c2err,"Check BiasCor file" );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }

    fflush(stdout);
  }

  if ( ipar_LCFIT == INDEX_mB ) 
    { simdata_bias.NUSE += NEVT_USE ; } // sum over IDSAMPLE
  
  // -----------------
  // debug dump
  LDMP = 0 ; 
  if ( LDMP ) {
    iz=7; ix1=5; ic=6; ia=0; ib=0;
    J1D = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] ;
    printf(" xxx --------------------------------------- \n");
    printf(" xxx %s-bias = %.3f +- %.3f for \n"
	   " xxx \t z[%.3f:%.3f] x1[%.3f:%.3f] c[%.3f:%.3f]  N=%d\n"
	   " xxx \t alpha=%.3f  beta=%.3f  \n"
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[ipar_LCFIT].varName
	   ,simdata_bias.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT] 
	   ,simdata_bias.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT] 
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.lo[iz]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.hi[iz]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].lo[ix1]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].hi[ix1]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].lo[ic]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].hi[ic]
	   ,CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D] 
	   ,simdata_bias.BININFO_SIM_ALPHA.avg[ia]
	   ,simdata_bias.BININFO_SIM_BETA.avg[ib]
	   );
    fflush(stdout);
    debugexit(fnam);
  }
  // ------------

  fflush(stdout);

  free(SUMBIAS); free(SUMWGT); free(sum); free(sumsq);

  return ;

} //end makeMap_fitPar_biasCor

// =============================
void get_J1DNBR_LIST(int IDSAMPLE, int J1D, int *NJ1DNBR, int *J1DNBR_LIST) {

  // Created May 2016
  // For input J1D index, return
  //  *NJ1DNBR = number of bin neighbors in space of z,x1,c.
  //  *J1DNBR_LIST = list of J1D indices for neighbors
  //
  // Usually NJ1DNBR = 3^3 = 27, but near boundaries it can be less.
  //

  int NBINz  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINx1 = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;

  int iz, ix1, ic ;
  int IA, IB, IZ, IX1, IC ;
  int NJ1D_local = 0 ;
  int J1D_local ;

  //  char fnam[] = "get_J1DNBR_LIST";

  // ----------------- BEGIN --------------

  // convert J1D into index for each dimension
  J1D_invert_I(IDSAMPLE, J1D, &IA, &IB, &IZ, &IX1, &IC );

  /*
  printf(" xxx J1D=%d -> IA,IB,IZ,IX1,IC = %d %d %d %d %d \n",
	 J1D, IA,IB,IZ,IX1,IC ); fflush(stdout); 
  */

  for(iz=IZ-1; iz <= IZ+1; iz++ ) {
    if ( iz <  0     ) { continue ; }
    if ( iz >= NBINz ) { continue ; }

    for(ix1=IX1-1; ix1<=IX1+1 ; ix1++ ) {
    if ( ix1 <  0      ) { continue ; }
    if ( ix1 >= NBINx1 ) { continue ; }

      for(ic=IC-1; ic<=IC+1; ic++ ) {
	
	if ( ic <  0     ) { continue ; }
	if ( ic >= NBINc ) { continue ; }

	J1D_local = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[IA][IB][iz][ix1][ic];
	J1DNBR_LIST[NJ1D_local] = J1D_local ;
	NJ1D_local++ ;
      }
    }
  }

  // load return arg.
  *NJ1DNBR  = NJ1D_local ;

  return ;

} // end get_J1DNBR_LIST

// ============================
double WGT_biasCor(int opt, int ievt, char *msg ) {

  // For biasCor event "ievt", return WGT to use for
  // averaging bias within each bin.
  // *msg is for error message only.
  //
  // opt=1 --> 1/muerr^2 only
  // opt=2 --> muerr^-2 x [wgt from cell location]
  //

  int  idsample = simdata_bias.idsample[ievt];
  int  istat_cov ;
  double WGT, muerrsq ;
  char fnam[] = "WGT_biasCor" ;
  
  // --------------- BEGIN -------------------

  muerrsq  = muerrsq_biasCor(ievt, USEMASK_BIASCOR_COVTOT, &istat_cov) ; 
  WGT      = 1.0/muerrsq ;

  if ( opt == 1 ) { return(WGT); }

  if ( INPUTS.sigma_cell_biasCor > 10.0 ) { return(WGT); }

  // If we get here, compute additional weight based on 3D 
  // separation from wgted-avg in cell.
  double z, x1, c, Dz, Dx1, Dc, SQD, ARG, WGT_CELL, SQSIGMA_CELL ;
  double binSize_z, binSize_x1, binSize_c ;
  int  J1D     = J1D_biasCor(ievt,fnam);   // 1D index
  binSize_z    = CELLINFO_BIASCOR[idsample].BININFO_z.binSize ;
  binSize_x1   = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1].binSize ;
  binSize_c    = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c].binSize ;
  SQSIGMA_CELL = INPUTS.sigma_cell_biasCor * INPUTS.sigma_cell_biasCor ;
  z   = (double)simdata_bias.z[ievt] ;
  x1  = (double)simdata_bias.FITVAL[INDEX_x1][ievt] ;
  c   = (double)simdata_bias.FITVAL[INDEX_c][ievt] ;
  Dz  = (z  - CELLINFO_BIASCOR[idsample].AVG_z[J1D])/binSize_z ;
  Dx1 = (x1 - CELLINFO_BIASCOR[idsample].AVG_LCFIT[INDEX_x1][J1D])/binSize_x1 ;
  Dc  = (c  - CELLINFO_BIASCOR[idsample].AVG_LCFIT[INDEX_c ][J1D])/binSize_c ;
  SQD = Dz*Dz + Dx1*Dx1 + Dc*Dc ;
  ARG = -0.5 * (SQD/SQSIGMA_CELL);
  WGT_CELL = exp(ARG);
  WGT     *= WGT_CELL ;

  return(WGT);

} // end WGT_biasCor

// ===========================
int J1D_biasCor(int ievt, char *msg ) {

  // Created May 12 2016
  // For biasCor event "ievt", return 1D index to biasCor map.
  // *msg is used only for error message.

  double a,b,z,x1,c ;
  char *name, MSGERR_IBINFUN[100] ;
  int ia, ib, iz, ix1, ic, J1D, idsample ;
  char fnam[] = "J1D_biasCor" ;

  // -------------- BEGIN ------------

  J1D = -9 ;

  
  // get variables defining GRID map
  a     = (double)simdata_bias.SIM_ALPHA[ievt] ;
  b     = (double)simdata_bias.SIM_BETA[ievt] ;
  z     = (double)simdata_bias.z[ievt];
  x1    = (double)simdata_bias.FITVAL[INDEX_x1][ievt] ;
  c     = (double)simdata_bias.FITVAL[INDEX_c][ievt] ;
  name     = simdata_bias.name[ievt];
  idsample = simdata_bias.idsample[ievt];

  // find integer bin for each map variable to store bias
  sprintf(MSGERR_IBINFUN,"%s : CID=%s  z=%.3f x1=%.3f c=%.3f", 
	  fnam, name, z, x1, c);

  ia  = IBINFUN(a,  &simdata_bias.BININFO_SIM_ALPHA,
		1, MSGERR_IBINFUN );
  ib  = IBINFUN(b,  &simdata_bias.BININFO_SIM_BETA,
       		  1, MSGERR_IBINFUN );
  iz  = IBINFUN(z,  &CELLINFO_BIASCOR[idsample].BININFO_z, 
		1, MSGERR_IBINFUN );
  ix1 = IBINFUN(x1, &CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1], 
		1, MSGERR_IBINFUN );
  ic  = IBINFUN(c,  &CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c], 
		1, MSGERR_IBINFUN );

  int NBINz  = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
  int NBINx1 = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc  = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c].nbin ;
  if ( iz  < 0 || iz  >= NBINz  ) { return J1D; }
  if ( ix1 < 0 || ix1 >= NBINx1 ) { return J1D; }
  if ( ic  < 0 || ic  >= NBINc  ) { return J1D; }
  
  J1D = CELLINFO_BIASCOR[idsample].MAPCELL[ia][ib][iz][ix1][ic];

  if ( J1D >= MAXBIN_BIASCOR_1D || J1D < 0 ) 
    {  J1D = -9 ;  }

  return(J1D);

} // end J1D_biasCor


// ======================================================
void  J1D_invert_D(int IDSAMPLE, int J1D, 
		   double *a, double *b, 
		   double *z, double *x1, double *c) {

  // May 2016: for input J1D index, return doubles a,b,z,x1,c

  int NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;
  int NBINz   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINx1  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;
  int ia, ib, iz, ix1, ic, J1TMP ;

  // --------------- BEGIN ----------------

  *a = *b = *z = *x1 = *c = -999.9 ;

  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(iz=0; iz<NBINz; iz++ ) {
	for(ix1=0; ix1<NBINx1; ix1++ ) {
	  for(ic=0; ic<NBINc; ic++ ) {

	    J1TMP = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] ;
	    if ( J1TMP == J1D ) {
	      *a  = simdata_bias.BININFO_SIM_ALPHA.avg[ia];
	      *b  = simdata_bias.BININFO_SIM_BETA.avg[ib];
	      *z  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.avg[iz];
	      *x1 = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].avg[ix1];
	      *c  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].avg[ic];
	      return ;
	    }

	  }
	}
      }
    }
  }


} // end J1D_invert_D

// ======================================================
void  J1D_invert_I(int IDSAMPLE, int J1D, int *ja, int *jb, 
		   int *jz, int *jx1, int *jc) {

  // May 2016: for input J1D index, return doubles a,b,z,x1,c

  int NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;
  int NBINz   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINx1  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;
  int ia, ib, iz, ix1, ic, J1TMP ;

  // --------------- BEGIN ----------------

  *ja = *jb = *jz = *jx1 = *jc = -9 ;

  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(iz=0; iz<NBINz; iz++ ) {
	for(ix1=0; ix1<NBINx1; ix1++ ) {
	  for(ic=0; ic<NBINc; ic++ ) {

	    J1TMP = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] ;
	    if ( J1TMP == J1D ) {
	      *ja  = ia ;
	      *jb  = ib ;
	      *jz  = iz ;
	      *jx1 = ix1 ;
	      *jc  = ic ;
	      return ;
	    }

	  }
	}
      }
    }
  }


} // end J1D_invert_I


// ======================================================
int biasMapSelect(int ievt) {

  // return TRUE if this biasCor event should be included -->
  // + true Type is SNIa
  // + z,x1,c are all within the bias grid range.
  // + passes user cuts
  //
  // Dec 9 2016: check option MASK_BIASCOR_2x2ab
  // Mar 12 2017: increment NPASS_biasMapSelect

  double z, x1, c, a,b, zmin, zmax  ;
  double CUTVAL_LOCAL[MXCUTWIN];
  int icut, idsample, istage ;
  //  char fnam[] = "biasMapSelect";
  
  // --------------- BEGIN ---------------

  // init stage counters on first event
  if ( ievt == 0 ) {
    for(istage=0; istage < NSTAGE_biasMapSelect; istage++ ) 
      { NPASS_biasMapSelect[istage] = 0 ; }
  }

  istage = 0 ;NPASS_biasMapSelect[istage]++ ; // ALL

  idsample = simdata_bias.idsample[ievt]; 
  z        = (double)simdata_bias.z[ievt];
  x1       = (double)simdata_bias.FITVAL[INDEX_x1][ievt] ;
  c        = (double)simdata_bias.FITVAL[INDEX_c][ievt] ;

  if ( idsample < 0 ) { return(0); }
  istage++ ; NPASS_biasMapSelect[istage]++ ; // IDSAMPLE is defined

  // return if this idample is not selected (Sep 1 2016)
  if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // IDSAMPLE is selected

  // require true SNIa for bias correction
  if ( simdata_bias.SIM_NONIA_INDEX[ievt] != 0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // true SNIa

  // apply same CUTWIN cuts as for data
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
    { CUTVAL_LOCAL[icut] = (double)simdata_bias.CUTVAL[icut][ievt] ; }
  if ( apply_CUTWIN(2,simdata_bias.DOFLAG_CUTWIN,CUTVAL_LOCAL) == 0 )  
    { return 0  ; }

  istage++ ; NPASS_biasMapSelect[istage]++ ;  // CUTWIN

  // Dec 9 2016: check option to pick only the boundary nodes alpha,beta 
  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_2x2ab) > 0 ) {
    double epsilon = 0.0001 ;
    double amin = INPUTS_SAMPLE_BIASCOR.ALPHA_MIN + epsilon ;
    double amax = INPUTS_SAMPLE_BIASCOR.ALPHA_MAX - epsilon ;
    double bmin = INPUTS_SAMPLE_BIASCOR.BETA_MIN  + epsilon ;
    double bmax = INPUTS_SAMPLE_BIASCOR.BETA_MAX  - epsilon ;
    int    KEEP_a, KEEP_b ;
    a = (double)simdata_bias.SIM_ALPHA[ievt];
    b = (double)simdata_bias.SIM_BETA[ievt];
    KEEP_a = ( (a > amax) || ( a < amin ) );
    KEEP_b = ( (b > bmax) || ( b < bmin ) );
    if ( KEEP_a==0  ) { return(0); }
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // keep alpha

    if ( KEEP_b==0 ) { return(0); }
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // keep beta
  }
  else {
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // keep alpha
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // keep beta
  }



  // ----------------------------
  // apply legacy window cuts
  if ( x1 < INPUTS.x1min ) { return 0 ; }
  if ( x1 > INPUTS.x1max ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // x1 cut

  if ( c  < INPUTS.cmin  ) { return 0 ; }
  if ( c  > INPUTS.cmax  ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // c cut

  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_5D) || 
       (INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT)    ) {
    // make sure values are within grid
    int NZ =  CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    if ( SKIPZCUT_BIASCOR==0 ) {
      zmin = CELLINFO_BIASCOR[idsample].BININFO_z.lo[0] ;
      zmax = CELLINFO_BIASCOR[idsample].BININFO_z.hi[NZ-1] ;
      if ( z  < zmin ) { return 0 ; }
      if ( z  > zmax ) { return 0 ; }
    }
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // z in GRID

    if ( x1 < BIASCOR_MINVAL_LCFIT[INDEX_x1]  ) { return 0; }
    if ( x1 > BIASCOR_MAXVAL_LCFIT[INDEX_x1]  ) { return 0; }
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // x1 in GRID

    if ( c  < BIASCOR_MINVAL_LCFIT[INDEX_c]   ) { return 0; }
    if ( c  > BIASCOR_MAXVAL_LCFIT[INDEX_c]   ) { return 0; }
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // c in GRID
  }
  else {
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // z  in GRID
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // x1 in GRID
    istage++ ; NPASS_biasMapSelect[istage]++ ;  // c  in GRID
  }


  // remove events with crazy COV
  if ( fabsf(simdata_bias.COV_x0x1[ievt])   > 8.0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // sane COVx0x1

  if ( fabsf(simdata_bias.COV_x0c[ievt])    > 8.0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // sane COVx0c

  if ( fabsf(simdata_bias.COV_x1c[ievt])    > 8.0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // sane COVx1c
  
  if ( simdata_bias.FITERR[INDEX_mB][ievt] > 1.0 ) { return 0 ; }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // mBerr cut

  // check pre-scale option
  if ( INPUTS.prescale_biasCor[1] > 1 ) {
    int   PS0  =  INPUTS.prescale_biasCor[0] ;
    int   PS1  =  INPUTS.prescale_biasCor[1] ;
    float XEVT =  (float)ievt ;
    float XPS1 =  (float)PS1 ;
    if ( fmodf( XEVT, XPS1 ) != PS0 ) { return 0 ; } 
  }
  istage++ ; NPASS_biasMapSelect[istage]++ ;  // sim prescale

  return(1);

} // end biasMapSelect

void  dumpStages_biasMapSelect(void) {

  // Created Mar 2017
  // If too few events pass biasMapSelect, called this function
  // to dump NPASS vs. stage from biasMapSelect.
  // 
  int istage;
  char fnam[] = "dumpStages_biasMapSelect" ;

  // ---------- BEGIN -----------
  
  printf("\n DUMP NPASS vs. Stage for biasMapSelect: \n");

  for(istage=0; istage < NSTAGE_biasMapSelect; istage++ ) {
    printf("\t %2d) NPASS = %5d   (%s) \n",
	   istage, 
	   NPASS_biasMapSelect[istage], 
	   COMMENT_biasMapSelect[istage]  ) ;
    fflush(stdout);
  }

  return ;

} // end dumpStages_biasMapSelect

// ======================================================
void  zero_FITPARBIAS(FITPARBIAS_DEF *FITPARBIAS) {
  int ipar;
  for(ipar=0; ipar < NLCPAR; ipar++ ) {
    FITPARBIAS->VAL[ipar] = 0.0 ;
    FITPARBIAS->ERR[ipar] = 0.0 ;
    FITPARBIAS->RMS[ipar] = 0.0 ;
  }
} // end zero_FITPARBIAS

// ======================================================
void makeMap_sigmu_biasCor(int IDSAMPLE) {

  // Created May 12, 2016
  // Make map of sigma_mu bias as a function of (z,x1,c,a,b).
  //
  // To get better stats, the bias is just a function of z
  // but the 5D binning is still used so that the
  // on-the-fly lookup is done the same way as other bias
  // corrections.
  //
  //
  // Beware that the binning for alpha,beta,z is the same
  // for biasCor and sigmu-scale ... but the color bins 
  // are different; courser bins for sigmu-scale.
  //

  int NBINa    = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int NBINb    = simdata_bias.BININFO_SIM_BETA.nbin ;  
  int NBIASCOR_CUTS = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;

  int    NBINz, NBINc ;
  double cmin,  cmax, cbin, c_lo, c_hi;

  int DUMPFLAG = 0 ;
  int    ia, ib, iz, ic, i1d, NCELL, isp ; 
  int    ievt, istat_cov, istat_bias, N, J1D, ipar ;
  double muErr, muErrsq, muDif, muDifsq, pull, tmp1, tmp2  ;
  double muBias, muBiasErr, muCOVscale, fitParBias[NLCPAR] ;
  double a,b, z, c ;
  double *SUM_MUERR, *SUM_SQMUERR;
  double *SUM_MUDIF, *SUM_SQMUDIF ;
  double *SQMUERR,   *SQMURMS ;
  double *SUM_PULL,  *SUM_SQPULL ;
  char  *name ;

  // MAXBIN_BIASCOR_1D

  BIASCORLIST_DEF     BIASCORLIST ;
  FITPARBIAS_DEF      FITPARBIAS[MXa][MXb] ;
  double              MUCOVSCALE[MXa][MXb] ;
  INTERPWGT_ALPHABETA INTERPWGT ;
 
  char fnam[]  = "makeMap_sigmu_biasCor" ;
  
  // ----------------- BEGIN -------------------

  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return; }

  printf(" %s: make map of muCOVscale(a,b,z,c) \n", fnam);
  fflush(stdout);

  // redshift bins are same as for biasCor
  copy_BININFO(&CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 
	       &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z);
  NBINz = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.nbin ;
  
  // setup color bins that are coarser than those for muBias.
  // Note that other bins (a,b,z) are same for muCOVscale and muBias.
  // Beware hard-wired values here.
  NBINc=3; cmin=-0.3; cmax=+0.3; cbin=0.2; 
  // NBINc=1; cmin=-0.3; cmax=+0.3; cbin=0.6;  // return to < Jun 30 2016


  NCELL    = NBINa * NBINb * NBINz * NBINc ;
  CELLINFO_MUCOVSCALE[IDSAMPLE].NCELL = NCELL ;

  // ---------------------------------------------
  // malloc global struct to store CELLINFO
  printf("\t Malloc CELL-INFO arrays with size=%d \n", NCELL);
  int MEMD     = NCELL   * sizeof(double);
  int MEMI     = NCELL   * sizeof(int);
    
  CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell =  (int    *) malloc(MEMI);
  CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z    =  (double *) malloc(MEMD);
  CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c] = (double *) malloc(MEMD);

  // but color bins are different ...
  // this is redundant for each IDSAMPLE ???
  CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin    = NBINc ;
  CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].binSize = cbin ;
  for(ic=0; ic < NBINc; ic++ ) {
    c_lo = cmin + cbin * (double)ic ;
    c_hi = c_lo + cbin ;    
    CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].lo[ic]  = c_lo ; ;
    CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].hi[ic]  = c_hi ; ;
    CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].avg[ic] 
      = 0.5*(c_lo+c_hi) ;
  }


  // malloc local 1D arrays to track local sums.
  SUM_MUERR    = (double*) malloc(MEMD);
  SUM_SQMUERR  = (double*) malloc(MEMD);
  SUM_MUDIF    = (double*) malloc(MEMD);
  SUM_SQMUDIF  = (double*) malloc(MEMD);
  SQMUERR      = (double*) malloc(MEMD);
  SQMURMS      = (double*) malloc(MEMD);
  SUM_PULL     = (double*) malloc(MEMD);
  SUM_SQPULL   = (double*) malloc(MEMD);

  int N1D=0;
  for(ia=0; ia< NBINa; ia++ ) {
    for(ib=0; ib< NBINb; ib++ ) {  
      MUCOVSCALE[ia][ib] = 1.0 ; // dummy arg for get_muBias below
      for(iz=0; iz < NBINz; iz++ ) {
	for(ic=0; ic < NBINc; ic++ ) {
	  SUM_MUERR[N1D] = SUM_SQMUERR[N1D] = 0.0 ;
	  SUM_MUDIF[N1D] = SUM_SQMUDIF[N1D] = 0.0;	
	  SUM_PULL[N1D]  = SUM_SQPULL[N1D]  = 0.0 ;

	  simdata_bias.MUCOVSCALE[IDSAMPLE][N1D] = 1.0 ;

	  CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[N1D]  = 0 ;
	  CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z[N1D]     = 0.0 ;
	  CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c][N1D] = 0.0 ;
	  CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][iz][0][ic] = N1D ;
	  N1D++ ;
	}
      }
    }
  }

  for(isp=0; isp < NBIASCOR_CUTS; isp++ ) {

    ievt = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;

    // check if there is valid biasCor for this event
    J1D = J1D_biasCor(ievt,fnam);
    if ( CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D] < BIASCOR_MIN_PER_CELL ) 
      { continue ; }

    for(ia=0; ia<MAXBIN_BIASCOR_ALPHA; ia++ ) {
      for(ib=0; ib<MAXBIN_BIASCOR_BETA; ib++ ) {
	zero_FITPARBIAS(&FITPARBIAS[ia][ib]);
      }
    }
    
    z    = (double)simdata_bias.z[ievt];
    a    = (double)simdata_bias.SIM_ALPHA[ievt];
    b    = (double)simdata_bias.SIM_BETA[ievt];
    c    = (double)simdata_bias.FITVAL[INDEX_c][ievt];
    ia   = simdata_bias.IA[ievt];
    ib   = simdata_bias.IB[ievt];
    name = simdata_bias.name[ievt];
    
    // allow color (c) to be outside map
    iz = IBINFUN(z, &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z, 
		 1, fnam );
    ic = IBINFUN(c, &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c], 
		 2, fnam );

    // ---------------------------------------------------
    // need bias corrected distance to compute pull
    for(ipar=0; ipar < NLCPAR; ipar++ ) 
      { BIASCORLIST.FITPAR[ipar] = (double)simdata_bias.FITVAL[ipar][ievt]; }
    BIASCORLIST.z     = z ;
    BIASCORLIST.alpha = a ;
    BIASCORLIST.beta  = b ;
    BIASCORLIST.idsample = IDSAMPLE ;

    istat_bias = 
      get_fitParBias(name, &BIASCORLIST, DUMPFLAG,
		     &FITPARBIAS[ia][ib] ); // <== returned

    // skip if bias cannot be computed, just like for data
    if ( istat_bias == 0 ) { continue ; }
    
    fcn_AlphaBetaWGT(a,b, DUMPFLAG, &INTERPWGT, fnam );
    get_muBias(name, &BIASCORLIST, FITPARBIAS, MUCOVSCALE, &INTERPWGT,
	       fitParBias, &muBias, &muBiasErr, &muCOVscale );  

    // ----------------------------
    muDif   =  muresid_biasCor(ievt);  // mu - muTrue
    muDif  -=  muBias ;    
    muDifsq =  muDif*muDif ;

    // compute error with intrinsic scatter
    muErrsq = muerrsq_biasCor(ievt, USEMASK_BIASCOR_COVTOT, &istat_cov) ; 
    muErr   = sqrt(muErrsq);
    
    pull = (muDif/muErr);

    // get 1d index
    i1d = CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][iz][0][ic] ;

    SUM_PULL[i1d]    += pull ;
    SUM_SQPULL[i1d]  += (pull*pull) ;
    SUM_SQMUDIF[i1d] +=  muDifsq ;
    SUM_MUDIF[i1d]   +=  muDif ;
    SUM_SQMUERR[i1d] +=  muErrsq ;
    SUM_MUERR[i1d]   +=  muErr ;

    // increment sums to get average in each cell   
    CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[i1d]++ ;
    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z[i1d]              += z ;
    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c][i1d] += c ;

  } // end ievt


  // -------------------------------------------------
  double XN, SQRMS, RMS ;
  
  for(i1d=0; i1d < NCELL; i1d++ ) {

    SQMURMS[i1d] = 0.0 ;
    SQMUERR[i1d] = 0.0 ;
    
    N  = CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[i1d] ;
    XN = (double)N ;

    if ( N < 5 ) { continue ; }

    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z[i1d]              /= XN ;
    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c][i1d] /= XN ;

    tmp1= SUM_MUDIF[i1d]/XN ;    tmp2= SUM_SQMUDIF[i1d]/XN ;
    SQRMS = tmp2 - tmp1*tmp1 ;
    SQMURMS[i1d] = SQRMS ;

    // average calculated muErrsq
    muErr        = SUM_MUERR[i1d]/XN ;
    muErrsq      = muErr*muErr  ;
    SQMUERR[i1d] = muErrsq ;
    
    // RMS of pull
    tmp1= SUM_PULL[i1d]/XN;  tmp2=SUM_SQPULL[i1d]/XN ;
    SQRMS = tmp2 - tmp1*tmp1 ;
    simdata_bias.MUCOVSCALE[IDSAMPLE][i1d] = SQRMS ;

  }  // i1d
  

  // -------------------------
  // print errBias info in z bins

  int LPRINT = 0 ;

  if ( LPRINT ) {

    double zlo, zhi ;  
    printf("\n");
    printf("                            "
	   "RMS(muDif)/RMS(Pull)/NSIM for \n");
    
    printf("  ia,ib  z-range :   "
	   "    ic=0               ic=1                ic=2 \n");
    
    printf("  -------------------------------------------------"
	   "------------------------\n");
    fflush(stdout);
    
    for(ia=0; ia< NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {      
	for(iz=0; iz < NBINz; iz++ ) {	
	  zlo = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.lo[iz];
	  zhi = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.hi[iz];
	  printf("  %d,%d  %.2f-%.2f : ",  ia,ib, zlo, zhi );
	  
	  for(ic=0; ic<NBINc; ic++ ) {  
	    i1d   = CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][iz][0][ic] ;
	    N     = CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[i1d] ;
	    muCOVscale = simdata_bias.MUCOVSCALE[IDSAMPLE][i1d] ;
	    RMS        = sqrt ( SQMURMS[i1d] );       
	    printf("%6.3f/%5.3f/%5d  ", RMS, sqrt(muCOVscale), N );	  
	  } // ic
	  
	  printf("\n");     fflush(stdout);
	} // iz
	printf("\n");     fflush(stdout);
      } // ib
    } // ia

  } // end LPRINT


  free(SUM_MUERR);   free(SUM_SQMUERR);
  free(SUM_MUDIF);   free(SUM_SQMUDIF) ;
  free(SQMUERR);     free(SQMURMS);
  free(SUM_PULL);    free(SUM_SQPULL);
  
  return;

} // end makeMap_sigmu_biasCor


// ======================================================
void  store_iaib_biasCor(void) {

  // store alpha, beta, zcmb (ia,ib, iz) for each biasCor event.

  double a, b, z, Z ;
  int i, IA, IB, IZ, iz, ID;  
  char fnam[] = "store_iaib_biasCor" ;

  // ------------ BEGIN ----------

  for(i=0; i < simdata_bias.NROW; i++ ) {

    if ( simdata_bias.SIM_NONIA_INDEX[i] != 0 ) { continue ; }
    a      = (double)simdata_bias.SIM_ALPHA[i] ;
    b      = (double)simdata_bias.SIM_BETA[i] ;  
    Z      = (double)simdata_bias.SIM_ZCMB[i] ;  // true ZCMB
    z      = (double)simdata_bias.z[i] ;         // measured zcmb
    ID     = simdata_bias.idsample[i];
    if ( ID < 0 ) { continue ; }

    IA     = IBINFUN(a, &simdata_bias.BININFO_SIM_ALPHA, 1, fnam ) ;
    IB     = IBINFUN(b, &simdata_bias.BININFO_SIM_BETA,  1, fnam ) ;
    simdata_bias.IA[i] = IA ;
    simdata_bias.IB[i] = IB ;

    IZ  = IBINFUN(Z, &CELLINFO_BIASCOR[ID].BININFO_z, 0, fnam ) ;
    iz  = IBINFUN(z, &CELLINFO_BIASCOR[ID].BININFO_z, 0, fnam ) ;
    simdata_bias.IZ[i] = IZ ;
    simdata_bias.iz[i] = iz ; 

  } // end loop over biasCor rows
  return ; 
} // end store_iaib_biasCor


// ==================================
double muresid_biasCor(int ievt ) {

  // Return muFit - muTrue for biasCor sample event "ievt"
  // Note that there is no bias correction here, just the
  // naive distance from Trip formula.

  double z, a, b, M0, mB, x1, c, zTrue, muFit, muTrue, muz, muDif ;
  double dlz, dlzTrue, dmu, cosPar[10] ;

  //  char fnam[] = "muresid_biasCor" ;

  // ----------------- BEGIN ----------------

  a      = (double)simdata_bias.SIM_ALPHA[ievt] ;
  b      = (double)simdata_bias.SIM_BETA[ievt] ;
  M0     = INPUTS.nommag0 ;
  z      = (double)simdata_bias.z[ievt] ;
  mB     = (double)simdata_bias.FITVAL[INDEX_mB][ievt] ; 
  x1     = (double)simdata_bias.FITVAL[INDEX_x1][ievt] ;
  c      = (double)simdata_bias.FITVAL[INDEX_c][ievt] ;

  // need true MU at observed redshift z. We don't have the biasCor
  // cosmology parameters, so we'll use a derivative to compute
  //     muTrue = SIM_DLMAG + (z-zTrue)*dmu/dz,
  // where dmu/dz is computed from the SALT2mu-reference cosmology.
  // Assumption here is that correction to SIM_DLMAG is very small
  // so that using approximate cosmology for dmu/dz won't matter.

  zTrue  = (double)simdata_bias.SIM_ZCMB[ievt];
  muTrue = (double)simdata_bias.SIM_DLMAG[ievt];

  // use measured z including zpec
  dlz      = cosmodl_forFit(z,    INPUTS.COSPAR); 

  // use true zCMB
  dlzTrue  = cosmodl_forFit(zTrue, INPUTS.COSPAR); 

  dmu    = 5.0*log10(dlz/dlzTrue) ;
  muz    = muTrue + dmu ;
  muFit  = mB + a*x1 - b*c - M0 ; // magOff doesn't matter for RMS
  muDif  = muFit - muz ;
  
  return(muDif);

} // end muresid_biasCor


// ===============================================================
double muerrsq_biasCor(int ievt, int maskCov, int *istat_cov) {

  // Created Sep 04 2017
  // Return squre of distance-modulus error for event "ievt"
  // of biasCor sample. 
  // 
  // Inputs:
  //   ievt   = event number for biasCor sample
  //   maskCov & 1 ==> include SALT2-fitted COV 
  //   maskCov & 2 ==> include intrinsic COV
  //
  // Output:
  //   istat_cov is returned non-zero if the cov matrix was tweaked.
  //   Function returns square of error on distance modulus.
  //

  int    IDSAMPLE = simdata_bias.idsample[ievt] ;
  double z        = (double)simdata_bias.z[ievt] ;
  double zerr     = (double)simdata_bias.zerr[ievt] ;
  int    iz       = IBINFUN(z,  &INPUTS.BININFO_z, 0, "" );
  double a        = (double)simdata_bias.SIM_ALPHA[ievt] ;
  double b        = (double)simdata_bias.SIM_BETA[ievt] ;
  int    ia       = simdata_bias.IA[ievt] ;
  int    ib       = simdata_bias.IB[ievt] ;
  char *name      = simdata_bias.name[ievt];

  int DOCOVFIT     = ( (maskCov & USEMASK_BIASCOR_COVFIT)>0 ) ;
  int DOCOVINT     = ( (maskCov & USEMASK_BIASCOR_COVINT)>0 ) ;
  int IDEAL_COVINT = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;

  double COVTOT[NLCPAR][NLCPAR] ;
  double COVINT[NLCPAR][NLCPAR] ;
  double COVTMP, muErrsq ;
  int  j0, j1 ;
  char fnam[] = "muerrsq_biasCor" ;
  
  // ------------- BEGIN -----------

  //  check option to include intrinsic scatter matrix
  if ( DOCOVINT ) {
    // fetch COVINT[3][3]
    get_COVINT_biasCor(IDSAMPLE,z,a,b, COVINT);
  }

  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0;j1<NLCPAR;j1++ ) { 

      COVTMP = 0.0 ;

      if ( DOCOVFIT )
	{ COVTMP += (double)simdata_bias.COV_FITVAL[j0][j1][ievt]; }

      if ( DOCOVINT )
	{ COVTMP += COVINT[j0][j1]; }

      COVTOT[j0][j1] = COVTMP ;

    } 
  }


  // fix matrix if it has bad eigenvalues (just like for data)
  update_covMatrix( name, 0, NLCPAR, COVTOT, EIGMIN_COV, istat_cov ); 
  
  // compute error on distance, including covariances
  muErrsq = fcn_muerrsq(name, a, b, COVTOT, z, zerr, 0 );

  return(muErrsq);

} // end muerrsq_biasCor


// ===============================================================
void get_COVINT_biasCor(int IDSAMPLE, double z, double alpha, double beta, 
			double (*COVINT)[NLCPAR] ) {

  // Created Sep 5 2017
  // For input values(IDSAMPLE,z,alpha,beta), return
  // 3x3 intrinsic scatter matrix for biasCor sample: 
  //  --> return COVINT[3][3]
  //
  // Note that this function can be used for both the
  // biasCor sample and the data.
  //
  // Intrinsic scatter matrix is stored at the a,b 
  // corners, so need to interpolate.
  //
  // Nov 16 2018: alpha,beta min -> 1.0E-5 (was 1.0E-4)

  int  Na       = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = simdata_bias.BININFO_SIM_BETA.nbin ;
  INTERPWGT_ALPHABETA INTERPWGT ;
  int j0, j1, iz, ia, ib ;
  int DUMPFLAG_abWGT = 0;
  int LDMP = 0 ;
  double WGT, COVTMP ;
  char fnam[] = "get_COVINT_biasCor" ;

  // -------------- BEGIN ----------------

  // check for valid arguments
  if ( IDSAMPLE < 0 || z<0.0 || alpha < 1.0E-5 || beta < 1.0E-5 ) {
    sprintf(c1err, "Invalid argument(s)");
    sprintf(c2err, "IDSAMPLE=%d  z=%.4f  a=%.4f  b=%.4f", 
	    IDSAMPLE, z, alpha, beta );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  zero_COV(COVINT);

  if ( LDMP ) {
    printf("\n xxx ------------------------ \n");
    printf(" xxx IDSAMP=%d  z=%.4f  a=%.4f  b=%.4f \n",
	   IDSAMPLE, z, alpha, beta);
    fflush(stdout);
  }

  // convert redshift to index
  iz  = IBINFUN(z, &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 0, "" );

  // get alpha,beta interp weights
  fcn_AlphaBetaWGT(alpha, beta, DUMPFLAG_abWGT, &INTERPWGT, fnam );

  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0; j1<NLCPAR; j1++ )  { 

      for(ia=0; ia<Na; ia++ ) {
	for(ib=0; ib<Nb; ib++ ) {
	  WGT    = INTERPWGT.WGT[ia][ib] ; 
	  COVTMP = simdata_bias.COVINT[IDSAMPLE][iz][ia][ib].VAL[j0][j1];
	  COVINT[j0][j1] += ( WGT * COVTMP );
	}
      }

      if ( LDMP ) {
	printf(" xxx COVINT[%d][%d] = %f \n",
	       j0, j1, COVINT[j0][j1] ); fflush(stdout);
      }

    }  // j1
  }  // j0

  if ( LDMP ) {  debugexit(fnam); } 

  return ;

} // end get_COVINT_biasCor
 

void  write_COVINT_biasCor(void) {

  // Open text file [prefix].COVINT and each COV matrix
  // to table row  in text file.

  int  NSAMPLE  = NSAMPLE_BIASCOR ;
  int  Na       = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = simdata_bias.BININFO_SIM_BETA.nbin ;
  int  idsample, Nz, ia ,ib, iz, NEVT, k0 ;
  int  NVAR, IROW=0 ;
  FILE *fp;
  double z, a, b, cov, rho, covDiag, sigma[NLCPAR] ;
  char outFile[200], cline[200] ;
  char fnam[] = "void   write_COVINT_biasCor" ;

  // ------------ BEGIN ----------------

  sprintf(outFile,"%s.COVINT", INPUTS.PREFIX );
  fp = fopen(outFile,"wt");

  fprintf(fp, "# REDCOV01 = COVINT(mB,x1)/(sigma_mB*sigma_x1) \n");
  fprintf(fp, "# REDCOV02 = COVINT(mB,c )/(sigma_mB*sigma_c ) \n");
  fprintf(fp, "# REDCOV12 = COVINT(x1,c )/(sigma_x1*sigma_c ) \n");
  fprintf(fp, "# \n");

  NVAR=12;
#ifdef TEXTFILE_NVAR
  fprintf(fp, "NVAR: %d\n", NVAR);
#endif
  fprintf(fp, "VARNAMES: "
	  "ROW NEVT IDSAMPLE z alpha beta "    // 6
	  "sigma_mB sigma_x1 sigma_c "         // 3
	  "REDCOV01  REDCOV02 REDCOV12 "       // 3
	  "\n" );
  
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    Nz       = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    for(iz=0; iz < Nz; iz++ ) {
      for(ia=0; ia < Na; ia++ ) {
	for(ib=0; ib < Nb; ib++ ) {
	  
	  IROW++ ;
	  NEVT = simdata_bias.NEVT_COVINT[idsample][iz][ia][ib]; 
	  z        = CELLINFO_BIASCOR[idsample].BININFO_z.avg[iz] ;
	  a    = simdata_bias.BININFO_SIM_ALPHA.avg[ia] ;
	  b    = simdata_bias.BININFO_SIM_ALPHA.avg[ib] ;
	  cline[0] = 0 ;

	  sprintf(cline,"ROW: "
		  "%3d %4d %d %.3f %4.2f %4.2f ",
		  IROW, NEVT, idsample, z, a, b );

	  for(k0=0; k0<NLCPAR; k0++ ) {
	    covDiag = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[k0][k0] ;
	    sigma[k0]  = sqrt(covDiag);
	    sprintf(cline,"%s %7.4f ", cline, sigma[k0] );
	  }
	  
	  // REDCOV01
	  cov = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[0][1] ;
	  rho = cov / (sigma[0]*sigma[1]);
	  sprintf(cline,"%s %7.4f ", cline, rho );

	  // REDCOV02
	  cov = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[0][2] ;
	  rho = cov / (sigma[0]*sigma[2]);
	  sprintf(cline,"%s %7.4f ", cline, rho );

	  // REDCOV12
	  cov = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[1][2] ;
	  rho = cov / (sigma[1]*sigma[2]);
	  sprintf(cline,"%s %7.4f ", cline, rho );

	  fprintf(fp,"%s\n", cline);

	}
      }
    }
  }

  fclose(fp);
  return ;

} // end void   write_COVINT_biasCor

void  dump_COVINT_biasCor(int idsample, int iz, int ia, int ib) {

  // Screen dump COVINT for the input indices.

  int k0, k1, NEVT ;
  double covDiag, cov, sigma0, sigma1, rho ;
  double zlo = CELLINFO_BIASCOR[idsample].BININFO_z.lo[iz] ;
  double zhi = CELLINFO_BIASCOR[idsample].BININFO_z.hi[iz] ;
  char *parName ;
  char fnam[] = "dump_COVINT_biasCor" ;

  // ------------- BEGIN --------------

  NEVT = simdata_bias.NEVT_COVINT[idsample][iz][ia][ib]; 

  printf("\n");
  printf("# ----------------------------------------- \n");
  printf("# %s: \n", fnam);
  printf("#   idsample=%d  iz=%d(%.3f-%.3f)  ia,ib=%d,%d   NEVT=%d\n", 
	 idsample, iz,zlo,zhi, ia,ib, NEVT);


  // simdata_bias.COVINT[idsample][iz][ia][ib].VAL[0][0] = COV ;

  for(k0=0; k0 < NLCPAR; k0++ ) {

    covDiag = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[k0][k0] ;
    sigma0  = sqrt(covDiag);
    parName = BIASCOR_NAME_LCFIT[k0] ;

    printf("\t sigma(%2s) = %7.4f | ", parName, sigma0 );
    for(k1=0; k1 < NLCPAR; k1++ ) {
      // print reduced COV here ...
      covDiag = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[k1][k1] ;
      sigma1  = sqrt(covDiag);
      cov     = simdata_bias.COVINT[idsample][iz][ia][ib].VAL[k0][k1] ;
      rho     = cov/(sigma0*sigma1);
      printf(" %8.4f ", rho);
    }
    printf("\n");
  }

  return ;

} // end dump_COVINT_biasCor


void init_COVINT_biasCor(void) {

  // Jan 2018
  // Compute and store intrinsic scatter matrix 
  //
  // * For traditional model (IDEAL_COVINT=0), 
  //   set COV[0][0] = sigma_int^2
  //
  // * For emprical method (IDEAL_COVINT=1), compute COV in bins of
  //   idsample,redshift,alpha,beta
  //

  int  DO_IDEAL_COVINT = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  int  DO_COV00_ONLY   = ( DO_IDEAL_COVINT ==0 ) ;
  int  NSAMPLE  = NSAMPLE_BIASCOR ;
  int  NROW_TOT = simdata_bias.NROW ;
  int  Na       = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = simdata_bias.BININFO_SIM_BETA.nbin ;
  int  Nz, idsample, iz, ia, ib, ievt, ipar, ipar2 ; 
  double sigInt, COV, z ;
  char fnam[] = "init_COVINT_biasCor" ;

  //  int DOTEST_SIGINT_ONLY = 1 ; // traditional sigma_int model (debug only)

  // ------------- BEGIN ---------------

  sprintf(BANNER,"%s:", fnam);
  print_banner(BANNER); 

  printf("\t Compute Intrinsic Matrix (COVINT) from BiasCor Sample.\n");
  printf("\t COVINT computed in bins of IDSAMPLE, Redshift, Alpha, Beta\n");
  fflush(stdout);

  // zero out entire COV matrix.
  // For traditional sigma_int model (IDEAL_COVINT=0), set COV[0][0]
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    zero_COV(simdata_bias.COVINT_AVG[idsample].VAL) ;
    Nz       = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    for(iz=0; iz < Nz; iz++ ) {
      for(ia=0; ia < Na; ia++ ) {
	for(ib=0; ib < Nb; ib++ ) {
	  zero_COV(simdata_bias.COVINT[idsample][iz][ia][ib].VAL) ;
	  simdata_bias.NEVT_COVINT[idsample][iz][ia][ib] = 0 ;
	  if ( DO_COV00_ONLY ) {
	    sigInt = simdata_bias.SIGINT_ABGRID[idsample][ia][ib];
	    COV    = (sigInt * sigInt) ;
	    simdata_bias.COVINT[idsample][iz][ia][ib].VAL[0][0] = COV ;
	  }

	} // end iz
      } // and ib
    } // ene ia
  }  // end idsample

  if ( DO_COV00_ONLY ) { return ; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  // If we get here, compute full COV matrix in bins of 
  // IDSAMPLE,z,a,b
  //   COV(x,y) = sum[(x-xtrue)*(x-ytrue) ] / N

  //  SKIPZCUT_BIASCOR  = 0 ; // remove z-cut from biasCor selection 

  int NBIASCOR_IDEAL=0, NBIASCOR_CUTS=0 ;
  double FITPAR_IDEAL[NLCPAR], FITPAR_TRUE[NLCPAR];
  double tmpVal, tmpVal2, x0_IDEAL, mB_IDEAL ;

  for(ievt=0; ievt < NROW_TOT; ievt++ ) {

    // apply selection
    if ( biasMapSelect(ievt) == 0  ) { continue ; }
    NBIASCOR_CUTS++ ;

    // check for valid IDEAL fit params 
    tmpVal = simdata_bias.FITVAL_IDEAL[INDEX_mB][ievt];
    if ( tmpVal > 900.0 ) { continue ; }
    
    // mB_IDEAL is corrupt, so hack fix based on x0
    x0_IDEAL = (double)simdata_bias.x0_IDEAL[ievt] ;
    mB_IDEAL = - 2.5*log10(x0_IDEAL);
    simdata_bias.FITVAL_IDEAL[INDEX_mB][ievt] = (float)mB_IDEAL;

    // reject mB outliers:
    tmpVal = 
      simdata_bias.FITVAL_IDEAL[INDEX_mB][ievt] -
      simdata_bias.SIMVAL[INDEX_mB][ievt] ;
    if ( fabs(tmpVal) > 1.0 ) { continue ; }

    idsample = simdata_bias.idsample[ievt];
    ia       = simdata_bias.IA[ievt] ; // true alpha index
    ib       = simdata_bias.IB[ievt] ; // true beta index
    iz       = simdata_bias.IZ[ievt] ; // true zcmb index
    z        = simdata_bias.SIM_ZCMB[ievt] ;

    simdata_bias.NEVT_COVINT[idsample][iz][ia][ib]++ ;

    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      for(ipar2=ipar; ipar2 < NLCPAR; ipar2++ ) {

	tmpVal = 
	  simdata_bias.FITVAL_IDEAL[ipar][ievt] -
	  simdata_bias.SIMVAL[ipar][ievt] ;

	tmpVal2 = 
	  simdata_bias.FITVAL_IDEAL[ipar2][ievt] -
	  simdata_bias.SIMVAL[ipar2][ievt] ;

	/* xxxxxxxxxxxx
	if ( ievt < 100 && ipar==0 && ipar2==0 ) {
	  printf(" xxx diff[%d][%d] = %f, %f \n",
		 ipar, ipar2, tmpVal, tmpVal2 );
	}
	xxxxxxxxxxx*/

	simdata_bias.COVINT[idsample][iz][ia][ib].VAL[ipar][ipar2] 
	  += (tmpVal * tmpVal2) ;

	simdata_bias.COVINT[idsample][iz][ia][ib].VAL[ipar2][ipar] =
	  simdata_bias.COVINT[idsample][iz][ia][ib].VAL[ipar][ipar2] ;
      }
    }

    NBIASCOR_IDEAL++ ;

  } // end ievt loop


  printf("\t %d of %d BiasCor events have IDEAL fit params. \n",
	 NBIASCOR_IDEAL, NBIASCOR_CUTS);
  fflush(stdout) ;

  // - - - - - - - - - - - - - - - - - 

  // divide each sum-term by N, and load symmetric part of matrix
  int N;
  double XNINV;

  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    Nz       = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    for(iz=0; iz < Nz; iz++ ) {
      for(ia=0; ia < Na; ia++ ) {
	for(ib=0; ib < Nb; ib++ ) {
	  
	  N = simdata_bias.NEVT_COVINT[idsample][iz][ia][ib]; 
	  if ( N > 0 ) {
	    XNINV = 1.0/(double)N;
	    scale_COV(XNINV,simdata_bias.COVINT[idsample][iz][ia][ib].VAL);
	  }

	  if ( iz < -4 ) { dump_COVINT_biasCor(idsample,iz,ia,ib); }
	}
      }
    }
  }

  write_COVINT_biasCor();


  return ;

} // end init_COVINT_biasCor


// ======================================================
void  init_sigInt_biasCor_legacy(int IDSAMPLE) {

  // Estimate sigInt for biasCor sample using very high SNR subset.
  // Jun 17 2016: compute sigInt in each alpha & beta bin. 
  // Jan 25,2018
  //   + set SKIPZFLAG_BIASCOR flag
  //   + abort if NUSE[ia][ib] < 10
  // Oct 08 2018: 
  //   + new input IDSAMPLE to allow option of sigint(IDSAMPLE)


  int  NROW_TOT = simdata_bias.NROW ;
  int  DO_SIGINT_SAMPLE = ( INPUTS.opt_biasCor & MASK_BIASCOR_SIGINT_SAMPLE ) ;
  int  NROW_malloc ;
  int  LDMP = 0 ;
  int  i, istat_cov, NCOVFIX, ia, ib ;

  double muErrsq, muErr, muDif, muOff, SNRMAX, sigInt, tmp1, tmp2 ;
  double zero = 0.0 ;
  double SUMDIF[MXa][MXb], SUMDIFSQ[MXa][MXb], SUMERRSQ[MXa][MXb] ;    
  double *MUDIF[MXa][MXb], *MUERRSQ[MXa][MXb];
  int    NUSE[MXa][MXb], NTMP ;

  int NBINa    = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int NBINb    = simdata_bias.BININFO_SIM_BETA.nbin ;

  char fnam[]   = "init_sigInt_biasCor_legacy" ;

  // ------------------- BEGIN -------------------

  printf("\n");


  // ---------------------------------------
  // check option for user to fix sigmb_biascor
  sigInt = INPUTS.sigint_biasCor ;
  if ( sigInt >= 0.0 ) {
    if ( IDSAMPLE == 0 ) 
      { printf(" sigInt -> %.3f from user input sigmb_biascor key\n", sigInt);}
    simdata_bias.SIGINT_AVG = sigInt ;
    for(ia=0; ia < NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {
	simdata_bias.SIGINT_ABGRID[IDSAMPLE][ia][ib] = sigInt ;
      }
    }
    return ;
  }

  // ---------------------------------------
  // check option to compute sigInt for each IDSAMPLE
  if ( DO_SIGINT_SAMPLE ) {
    printf(" %s for IDSAMPLE=%d (%s): \n",
	   fnam, IDSAMPLE,  SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }
  else if ( IDSAMPLE==0 ) {
    printf(" %s for all IDSAMPLEs combined: \n", fnam);
  }
  else {
    // IDSAMPLE>0 and one sigInt -> set sigInt = sigInt[IDSAMPLE=0]
    for(ia=0; ia < NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {
	simdata_bias.SIGINT_ABGRID[IDSAMPLE][ia][ib] = 
	  simdata_bias.SIGINT_ABGRID[0][ia][ib] ; 
      } // end ib
    } // end ia
        
  } 


  // ---------------------------------------
  // quick pass with SNR cut to estimate size for malloc  
  NROW_malloc = 0 ;
  for(i=0; i < NROW_TOT; i++ ) {
    SNRMAX = (double)simdata_bias.SNRMAX[i] ;
    if ( simdata_bias.SIM_NONIA_INDEX[i] != 0 ) { continue ; }
    if ( SNRMAX < INPUTS.snrmin_sigint_biasCor) { continue ; }
    NROW_malloc++ ;
  }

  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      MUDIF[ia][ib]   = (double*) malloc ( NROW_malloc * sizeof(double) );
      MUERRSQ[ia][ib] = (double*) malloc ( NROW_malloc * sizeof(double) );
      NUSE[ia][ib] = NCOVFIX = 0 ;
      SUMDIF[ia][ib]  = SUMDIFSQ[ia][ib] = SUMERRSQ[ia][ib] = 0.0 ;
    }
  }
  

  // skip z-cut to ensure events with SNR>60.
  // This allows, for example, doing BBC fits only at high-z.
  SKIPZCUT_BIASCOR  = 1 ; 

  for(i=0; i < NROW_TOT; i++ ) {

    // get variables defining GRID map
    SNRMAX = (double)simdata_bias.SNRMAX[i] ;
    if ( SNRMAX < INPUTS.snrmin_sigint_biasCor ) { continue ; }

    // apply selection
    if ( biasMapSelect(i) == 0  ) { continue ; }

    // check option for IDSAMPLE-dependent sigint (Oct 2018)
    if ( DO_SIGINT_SAMPLE && IDSAMPLE != simdata_bias.idsample[i] ) 
      { continue ; }

    // mu(mB,x1,c) = muTrue
    muDif  = muresid_biasCor(i); 

    // compute error with no intrinsic scatter (just use data COVFIT)
    muErrsq = muerrsq_biasCor(i, USEMASK_BIASCOR_COVFIT, &istat_cov) ;

    if ( istat_cov < 0 ) { NCOVFIX++ ; }    

    ia     = simdata_bias.IA[i] ; // alpha index
    ib     = simdata_bias.IB[i] ; // beta index
    
    SUMDIF[ia][ib]   +=  muDif ;
    SUMDIFSQ[ia][ib] += (muDif*muDif);
    SUMERRSQ[ia][ib] += muErrsq ;

    NTMP = NUSE[ia][ib];
    MUDIF[ia][ib][NTMP]   = muDif ;
    MUERRSQ[ia][ib][NTMP] = muErrsq ;
    NUSE[ia][ib]++ ;

  } // end loop over biasCor sample


  SKIPZCUT_BIASCOR  = 0 ;

  // -------------------------------------------------
  double XN, SQRMS, RMS ;
  double pull, sigInt_bin, sigTmp_lo, sigTmp_hi ;
  double sigTmp, sigInt_store[20], rmsPull_store[20] ;
  double sumdif, sumdifsq ;
  int    NBIN_SIGINT, DOPRINT ;
  int OPT_INTERP=1 ;
  double ONE = 1.0 ;


  simdata_bias.SIGINT_AVG = 0.0 ;

  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {    

      XN     = (double)NUSE[ia][ib] ;

      if ( XN < 10.0 ) {
	sprintf(c1err,"NUSE[ia=%d,ib=%d] = %d is too small.",
		ia, ib, NUSE[ia][ib] );
	sprintf(c2err,"Check biasCor events with SNR> %.1f",
		INPUTS.snrmin_sigint_biasCor ) ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }

      tmp1   = SUMDIF[ia][ib]/XN ;
      tmp2   = SUMDIFSQ[ia][ib]/XN ;
      SQRMS  = tmp2 - tmp1*tmp1 ;
      muOff  = SUMDIF[ia][ib]/XN ; 

      // compute approx sigInt using average distance error
      muErrsq = SUMERRSQ[ia][ib]/XN ;

      if ( muErrsq > SQRMS ) { 
	muErr = sqrt(muErrsq);        RMS = sqrt(SQRMS);
	sprintf(c1err,"Avg computed muErr > RMS(mures) for SNRMAX>%.0f",
		INPUTS.snrmin_sigint_biasCor );
	sprintf(c2err,"<muErr>=%f   RMS(mures)=%f   ia=%d ib=%d", 
		muErr, RMS, ia, ib );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }

      sigInt = sqrt(SQRMS - muErrsq) ;  // approx sigInt

      if ( LDMP ) {
	printf(" xxx ia,ib=%d,%d RMS(mu)=%.3f  muErr=%.3f  "
	       "muOff=%.3f   sigInt=%.3f\n", 
	       ia,ib, sqrt(SQRMS), sqrt(muErrsq), 
	       muOff, sigInt ); fflush(stdout);
      }

      // Loop over trial sigInt values . . . . 
      NBIN_SIGINT = 0 ;
      sigInt_bin = 0.01 ;
      sigTmp_lo  = sigInt - 4.*sigInt_bin - 1.0E-7 ;
      sigTmp_hi  = sigInt + 8.*sigInt_bin ;
      
      // start with largest sigInt and decrease so that RMS is increasing
      // for the interp function below
      for(sigTmp = sigTmp_hi; sigTmp >= sigTmp_lo; sigTmp -= sigInt_bin ) {
	sumdif = sumdifsq = 0.0 ;
	if ( sigTmp < 0.0 ) { continue ; }
	
	for(i=0; i < NUSE[ia][ib]; i++ ) {
	  muDif     = MUDIF[ia][ib][i] - muOff ;
	  muErrsq   = MUERRSQ[ia][ib][i] ;  
	  pull      = muDif/sqrt(muErrsq + sigTmp*sigTmp);
	  sumdif   += pull ;
	  sumdifsq += (pull*pull);
	}

	tmp1=sumdif/XN;    tmp2=sumdifsq/XN;
	SQRMS  = tmp2 - tmp1*tmp1 ;
	rmsPull_store[NBIN_SIGINT] = sqrt(SQRMS);
	sigInt_store[NBIN_SIGINT]  = sigTmp ;
	
	if ( LDMP ) {
	  printf("\t xxx ia,ib=%d,%d : "
		 "sigTmp(%d)=%.3f --> RMS(pull)=%.3f \n",
		 ia, ib, NBIN_SIGINT, sigInt_store[NBIN_SIGINT], 
		 rmsPull_store[NBIN_SIGINT])  ;    fflush(stdout);
	}
	
	NBIN_SIGINT++ ;
      }  // end sigTmp loop

      
      // interpolate sigInt vs. rmsPull at rmsPull=1

      sigInt = interp_1DFUN(OPT_INTERP, ONE, NBIN_SIGINT,
			    rmsPull_store, sigInt_store, fnam);

      // load globals
      simdata_bias.SIGINT_ABGRID[IDSAMPLE][ia][ib] = sigInt ;
      simdata_bias.SIGINT_AVG += sigInt ;

      if ( IDSAMPLE==0 || DO_SIGINT_SAMPLE) { DOPRINT=1; } else { DOPRINT=0; }
      if ( DOPRINT ) {
	printf("\t sigInt[ia=%d,ib=%d] = %.3f (%d events with SNR>%.0f) \n", 
	       ia,ib, sigInt, NUSE[ia][ib], INPUTS.snrmin_sigint_biasCor);
	fflush(stdout);
      }
      
    } // ib
  } // ia


  simdata_bias.SIGINT_AVG /= (double)( NBINa*NBINb) ;

  // free memory
  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      free( MUDIF[ia][ib] )   ;
      free( MUERRSQ[ia][ib] ) ;  
    }
  } 
  
  return ;

} // end init_sigInt_biasCor_legacy


// ===========================================
void zero_COV(double (*COV)[NLCPAR] ) {
  int j0, j1;  
  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0;j1<NLCPAR;j1++ ) {  COV[j0][j1] = 0.0 ; }
  }
  return ;
} // end zero_COV

void scale_COV(double scale, double (*COV)[NLCPAR] ) {
  int j0, j1;  
  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0;j1<NLCPAR;j1++ ) {  COV[j0][j1] *= scale ; }
  }
  return ;
} // end scale_COV



// ==============================================
void  makeSparseList_biasCor(void) {

  // Created Dec 21 2017
  // make sparse BIASCOR list for this IDSAMPLE with cuts,
  // so that later loops are faster.

  int NROW     = simdata_bias.NROW ;
  int irow, isp, idsample, MEMI, NBIASCOR ;
  char fnam[] = "makeSparseList_biasCor" ;

  // ---------- BEGIN ------------

  sprintf(BANNER," %s: make sparse row-list for each biasCor sample.",fnam);
  print_banner(BANNER);    

  // allocate memory for sparse irow list
  for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) {
    NBIASCOR = SAMPLE_BIASCOR[idsample].NBIASCOR ;
    MEMI     = NBIASCOR * sizeof(int); 
    SAMPLE_BIASCOR[idsample].IROW_CUTS = (int*)malloc(MEMI);
    SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS = 0 ;
  }



  for(irow=0; irow < NROW; irow++ ) {

    // apply selection cuts for making biasCor map
    if ( biasMapSelect(irow) == 0 ) { continue ; }

    idsample = simdata_bias.idsample[irow] ;
    if ( idsample < 0 ) { continue ; }

    isp = SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS ;
    SAMPLE_BIASCOR[idsample].IROW_CUTS[isp] = irow;
    SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS++ ;
  }

  return;

} // end makeSparseList_biasCor

// ===============================================================
void  makeMap_binavg_biasCor(int IDSAMPLE) {

  // in each 3D cell, compute wgted avg of z,x1,c to use
  // as interpolation nodes. The bin-center is not the
  // right quantity for interpolation.
  //
  // 
  int NCELL   = CELLINFO_BIASCOR[IDSAMPLE].NCELL;
  int NROW    = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;
  int MEMD    = NCELL * sizeof(double);
  int MEMI    = NCELL * sizeof(int);

  int irow, isp, J1D, ipar, iz, *NperCell ;
  double WGT, z, fitPar[NLCPAR];
  double *SUM_WGT_5D, *SUM_z_5D, *SUM_FITPAR_5D[NLCPAR];
  char fnam[] = "makeMap_binavg_biasCor" ;

  // --------------- BEGIN ---------------

  NperCell    = (int*)    malloc(MEMI) ;
  SUM_z_5D    = (double*) malloc(MEMD) ;
  SUM_WGT_5D  = (double*) malloc(MEMD) ;
  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    { SUM_FITPAR_5D[ipar] = (double*) malloc(MEMD) ; }

  for(J1D = 0; J1D < NCELL; J1D++ ) {
    NperCell[J1D] = 0 ;
    SUM_WGT_5D[J1D] = SUM_z_5D[J1D]    = 0.0 ;
    for(ipar=0; ipar < NLCPAR; ipar++ )  { SUM_FITPAR_5D[ipar][J1D] = 0.0 ; }  

    // init bin avg to nonsense values
    CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D] = 9999. ;
    for(ipar=0; ipar < NLCPAR; ipar++ ) 
      { CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar][J1D] = 9999. ;  }  
  }


  // loop over biasCor sample ... xyz
  for(isp=0; isp < NROW; isp++ ) {

    irow = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;

    WGT = WGT_biasCor(1,irow,fnam) ;  // WGT=1/muerr^2
    J1D = J1D_biasCor(irow,fnam);     // 1D index

    z   = (double)simdata_bias.z[irow] ;
    iz  = IBINFUN(z,  &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 0, "" );

    NperCell[J1D]++ ;
    SUM_WGT_5D[J1D] += WGT ;
    SUM_z_5D[J1D]   += (WGT * z) ;
    
    for(ipar=0; ipar < NLCPAR ; ipar++ ) {
      fitPar[ipar]  = (double)simdata_bias.FITVAL[ipar][irow] ;
      SUM_FITPAR_5D[ipar][J1D] += ( WGT * fitPar[ipar] ) ;
    }

  } // end loop over biasCor rows


  // ------------------
  // loop over cells and compute avg
  for(J1D = 0; J1D < NCELL; J1D++ ) {
    WGT = SUM_WGT_5D[J1D] ;
    if ( WGT < 1.0E-9 ) { continue ; }

    CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D]  = 
      SUM_z_5D[J1D] / WGT ;

    for(ipar=0; ipar < NLCPAR; ipar++ ) { 
      CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar][J1D] = 
	SUM_FITPAR_5D[ipar][J1D]/WGT;
    }
  }

  // ---------------------------------------
  // free temp memory
  free(NperCell); free(SUM_WGT_5D); free(SUM_z_5D) ;
  for(ipar=0; ipar < NLCPAR; ipar++ )  { free(SUM_FITPAR_5D[ipar]); }

  return ;

} // end makeMap_binavg_biasCor

// =====================================================
void calc_zM0_biasCor(void) {

  // Created July 2016
  // compute wgted-average z in each user-redshift bin,
  // from total biasCor sample.
  // These averages are used to get model M0.

  int NBINz   = INPUTS.BININFO_z.nbin ;
  int NROW    = simdata_bias.NROW ;

  int iz, i ;
  double SUM_WGT[MAXBIN_z], SUM_z[MAXBIN_z];
  double WGT, z ;
  char fnam[] = "calc_zM0_biasCor" ;

  // --------- BEGIN -------------

  for(iz=0; iz < NBINz; iz++ )  { SUM_WGT[iz] = SUM_z[iz] = 0.0 ; }

  // loop over biasCor sample
  for(i=0; i < NROW; i++ ) {

    // apply selection cuts for making biasCor map
    if ( biasMapSelect(i) == 0 ) { continue ; }

    WGT = WGT_biasCor(1,i,fnam) ;  // WGT=1/muerr^2

    z   = (double)simdata_bias.z[i] ;
    iz  = IBINFUN(z,  &INPUTS.BININFO_z, 0, "" );

    SUM_z[iz]    += (WGT * z) ;
    SUM_WGT[iz]  += WGT ;
  }
  
  // store avg redshift for writing to M0 outFile
  for(iz=0; iz < NBINz; iz++ ) {
    WGT = SUM_WGT[iz] ;
    if ( WGT > 0.0 ) 
      { simdata_bias.zM0[iz] = SUM_z[iz] / WGT ; }
    else
      { simdata_bias.zM0[iz] = INPUTS.BININFO_z.avg[iz] ; }
  }

  return ;

} // end calc_zM0_biasCor 

// =====================================================
void calc_zM0_data(void) {

  // Created June 2017
  // compute binned redshift from inverting <mumodel> in each bin,
  // where <mumodel> is wgted average MU.
  // Output is used to write M0-vs-<z> to output for data.
  //
  // Jun 27 2017: REFACTOR z bins
  // Dec 21 2017:
  //  + minor refactor to make more clear that <MUMODEL> (not <MU>)
  //    is inverted to get <z>. No change in function outputs.
  //

  int NBINz   = INPUTS.BININFO_z.nbin ;

  int iz, i ;
  double SUM_WGT[MAXBIN_z], SUM_z[MAXBIN_z];
  double SUM_mu[MAXBIN_z], SUM_mumodel[MAXBIN_z];
  double WGT, mu, mumodel, muerr, muerrsq, a, b, z, zerr ;
  char *name ;
  char fnam[] = "calc_zM0_data" ;
  int  LDMP=0 ;
  // --------- BEGIN -------------

  for(iz=0; iz < NBINz; iz++ ) { 
    SUM_WGT[iz] = SUM_z[iz] = SUM_mu[iz] = SUM_mumodel[iz] = 0.0 ; 
    FITRESULT.zM0[iz] = simdata_bias.zM0[iz] ;
  }


  for (i=0; i < FITINP.NSNCUTS; i++ ) {
    
    if ( data[i].skipfit ) { continue ; }

    name     = data[i].name ;    
    z        = data[i].zhd;    
    zerr     = data[i].zhderr ;
    mu       = data[i].mu - FITRESULT.SNMAG0; 
    mumodel  = data[i].mumodel ;

    muerr    = data[i].muerr ;
    muerrsq  = muerr * muerr ;

    if ( muerr < 1.0E-8 ) {
      sprintf(c1err,"muerr(%s) = %f  (i=%d) \n", name, muerr, i);
      sprintf(c2err,"z = %.3f +- %0.3f ", z, zerr);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    iz         = IBINFUN(z,  &INPUTS.BININFO_z, 0, "" );
    WGT        = 1.0 / muerrsq ;
    SUM_z[iz]        += (WGT * z) ;
    SUM_mu[iz]       += (WGT * mu);
    SUM_mumodel[iz]  += (WGT * mumodel);
    SUM_WGT[iz]      += WGT ;
  }


  // store avg redshift for writing to M0 outFile

  double M0DIF, M0ERR, muAvg, muM0, zM0, zAvg ;
  double cosPar[10];

  // load reference cosmology params
  cosPar[0] = INPUTS.parval[IPAR_OL] ;   // OL
  cosPar[1] = INPUTS.parval[IPAR_Ok] ;  // Ok
  cosPar[2] = INPUTS.parval[IPAR_w0] ;  // w0
  cosPar[3] = INPUTS.parval[IPAR_wa] ;  // wa


  for(iz=0; iz < NBINz; iz++ ) {
    WGT = SUM_WGT[iz] ;

    if ( WGT < 1.0E-9 ) { continue ; }

    muAvg = SUM_mumodel[iz]/WGT ;
    zM0   = zmu_solve(muAvg,cosPar); // new
    FITRESULT.zM0[iz] = zM0 ;

    if ( LDMP ) {
      zAvg = SUM_z[iz] / WGT ;  // wgted avg redshift (for comparison)
      printf(" xxx iz=%2d (%.3f-%.3f) : "
	     "<z>(data,biasCor)=%.3f,%.3f  zM0 = %.3f  <mu>=%.3f\n", 
	     iz, INPUTS.BININFO_z.lo[iz],INPUTS.BININFO_z.hi[iz],
	     zAvg, simdata_bias.zM0[iz], zM0, muAvg);
      fflush(stdout) ;
    }

  } // end iz loop

  return ;

} // end calc_zM0_data


// =======================================
void test_zmu_solve(void) {
  double zTmp, muTmp = 34.;
  double cosTmp[4] = { 0.7, 0.0, -1.0, 0.0 } ;
  zTmp = zmu_solve(muTmp, cosTmp) ;
  printf(" xxx zmu_solve returns z = %f \n", zTmp);
  debugexit("test_zmu_solve");
}

// =======================================
double zmu_solve(double mu, double *cosPar) {

  // Created Jun 26 2017
  // For input distance modulus (mu) and cosmology parameters (cospar),
  // return redshift.
  // cosPar = OL,ok,w0,wa

  int    NITER=0;
  double z, dmu, DMU, dl, mutmp ;
  double DMU_CONVERGE = 1.0E-4 ;
  char fnam[] = "zmu_solve" ;
  // -------------- BEGIN ---------------

  DMU = 9999.0 ;
  z   = 0.5 ;
  while ( DMU > DMU_CONVERGE ) {    
    dl     = cosmodl_forFit(z,cosPar);
    mutmp  = 5.0*log10(dl) + 25.0 ;
    dmu    = mutmp - mu ;
    DMU    = fabs(dmu);
    z     *= (1-dmu/10.0) ;

    NITER++ ;
    if ( NITER > 500 ) {
      sprintf(c1err,"Could not solve for z after NITER=%d", NITER);
      sprintf(c2err,"mu=%f  dmu=%f  ztmp=%f", mu, dmu, z);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  } // end dz

  return(z) ;

} // end zmu_solve

// ======================================================
int  storeDataBias(int n, int DUMPFLAG) {

  // Created Mar 22 2016
  // for data event 'n', store bias for ipar = mB,x1,c.
  // Returns 1 if bias can be determined;
  // returns 0 otherwise (to be rejected)
  //
  // Each bias is a 5D function,
  //    mb-bias = function of (a=alpha,b=beta, z,x1,c)
  //    x1-bias = function of (a=alpha,b=beta, z,x1,c)
  //    c-bias  = function of (a=alpha,b=beta, z,x1,c)
  //
  // The data values z,x1,c are known a-prior and thus we
  // do the 3D interp here (before fitting) to get the bias 
  // as a function of the subset: z,x1,c.  Alpha & Beta, 
  // however, are SALT2mu fit parameters that are floated
  // in the fit,  and thus we store the bias in bins of 
  // alpha+beta so that the fit-function (fcn) can interpolate 
  // for each fit-trial in the minimization.
  //
  // July 1 2016: also store muCOVscale[ia][ib]
  // Apr 18 2017: fix aweful index bug ia -> ib for beta


  int  NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;
  BIASCORLIST_DEF BIASCORLIST ;
  int    ia, ib, ipar, istat_bias, idsample ;
  char   cid[40];
  char   fnam[] = "storeDataBias" ;

  // ------------- BEGIN -------------

  sprintf(cid,"%s", data[n].name);  // for err messages

  if ( DUMPFLAG ) {
    printf("\n");
    printf(" XXX ======================================== \n");
    printf(" XXX =========== %s DUMP for CID=%s ================ \n", 
	   fnam, cid );
    printf(" XXX ======================================== \n");
    fflush(stdout);
  }

  idsample = data[n].idsample ;
  if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { return(0); }


  BIASCORLIST.idsample = data[n].idsample ;
  BIASCORLIST.z        = data[n].zhd ;
  for(ipar=0; ipar<NLCPAR; ipar++ ) 
    { BIASCORLIST.FITPAR[ipar] = data[n].fitpar[ipar]; }

  for(ia = 0; ia < NBINa; ia++ ) {
    for(ib = 0; ib < NBINb; ib++ ) {
   
      BIASCORLIST.alpha = simdata_bias.BININFO_SIM_ALPHA.avg[ia];
      BIASCORLIST.beta  = simdata_bias.BININFO_SIM_BETA.avg[ib];

      istat_bias = 
	get_fitParBias(cid, &BIASCORLIST, DUMPFLAG,             // in
		       &data[n].FITPARBIAS_ALPHABETA[ia][ib] ); //out

      if ( DUMPFLAG ) {
	printf(" xxx %s: a=%.2f b=%.2f (ia,ib=%d,%d) istat_bias=%d \n",
	       fnam, BIASCORLIST.alpha, BIASCORLIST.beta, ia,ib, istat_bias);
	fflush(stdout);
      }

      if ( istat_bias == 0 ) { return 0 ; }

      istat_bias = 
	get_muCOVscale(cid, &BIASCORLIST, DUMPFLAG,             // in
		       &data[n].MUCOVSCALE_ALPHABETA[ia][ib] ); //out

      if ( DUMPFLAG ) {
	printf(" xxx %s: a=%.2f b=%.2f (ia,ib=%d,%d) istat_muCOVscale=%d \n",
	       fnam, BIASCORLIST.alpha, BIASCORLIST.beta, ia,ib, istat_bias);
	fflush(stdout);
      }

      if ( istat_bias == 0 ) { return 0 ; }

    } // end ib
  }  // end ia

  return(1);

} // end storeDataBias

// ======================================================
int  storeBias_CCprior(int n) {

  // Created Jun 7 2016
  // for CCprior event 'n', store bias for ipar = mB,x1,c.
  // Returns 1 if bias can be determined;
  // returns 0 otherwise (to be rejected)
  //
  // This routine is analagous to storeDataBias for real data.
  //
  // Apr 20 2017: fix awful index bug ia->ib to get beta

  int  NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;

  BIASCORLIST_DEF BIASCORLIST ;
  int    ia, ib, istat_bias ;
  char   cid[40];
  //  char   fnam[] = "storeBias_CCprior" ;

  // ------------- BEGIN -------------

  // extract cid for err messages
  sprintf(cid,"%s",  simdata_ccprior.SIMFILE_INFO_ALL.name[n]);  

  // load local BIASCORLIST struct
  BIASCORLIST.z =
    simdata_ccprior.SIMFILE_INFO_ALL.z[n] ;
  BIASCORLIST.FITPAR[INDEX_mB] =
    simdata_ccprior.SIMFILE_INFO_ALL.mB[n] ;
  BIASCORLIST.FITPAR[INDEX_x1] =
    simdata_ccprior.SIMFILE_INFO_ALL.x1[n] ;
  BIASCORLIST.FITPAR[INDEX_c] =
    simdata_ccprior.SIMFILE_INFO_ALL.c[n] ;
  BIASCORLIST.idsample =
    simdata_ccprior.SIMFILE_INFO_ALL.idsample[n] ; 
  
  for(ia = 0; ia < NBINa; ia++ ) {
    for(ib = 0; ib < NBINb; ib++ ) {
      
      BIASCORLIST.alpha = simdata_bias.BININFO_SIM_ALPHA.avg[ia];
      BIASCORLIST.beta  = simdata_bias.BININFO_SIM_BETA.avg[ib];

      istat_bias =
	get_fitParBias(cid, &BIASCORLIST, 0,
		       &simdata_ccprior.SIMFILE_INFO_ALL.FITPARBIAS_ALPHABETA[n][ia][ib]);

      if ( istat_bias == 0 ) { return 0 ; }

    } // end ib
  }  // end ia
  
  return(1);

} // end storeBias_CCPrior

// ======================================================
int get_fitParBias(char *cid, 
		   BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		   FITPARBIAS_DEF  *FITPARBIAS) {

  // Created May 2016
  // For input BIASCORLIST = { mB,x1,c, z, alpha, beta },
  // return FITPARBIAS[VAL,ERR,RMS] for  mB, x1, c.
  // CID is the SN name used only for error message.
  //
  // Function returns 1 if bias is found; returns 0 if no bias found.
  //
  // Apr 18 2017: enhance dump output.
  //
  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z  = BIASCORLIST->z ;
  double mB = BIASCORLIST->FITPAR[INDEX_mB];
  double x1 = BIASCORLIST->FITPAR[INDEX_x1];
  double c  = BIASCORLIST->FITPAR[INDEX_c];
  double a  = BIASCORLIST->alpha ;
  double b  = BIASCORLIST->beta ;
  int IDSAMPLE = BIASCORLIST->idsample ;

  // get bin info into local variables
  int  NBINz   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int  NBINx1  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int  NBINc   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;

  double BINSIZE_z  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.binSize ;
  double BINSIZE_x1 = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].binSize ;
  double BINSIZE_c  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].binSize ;
    
  int BADBIAS = 0 ;

  int J1D, IZ, IX1, IC, OK, USEZ[MAXBIN_z];
  int IZMIN, IZMAX, ICMIN, ICMAX, IX1MIN, IX1MAX ;
  int j1d, ia, ib, iz, ix1, ic, ipar ;
  int NperCell, NSUM_Cell, NCELL_INTERP_TOT, NCELL_INTERP_USE ;

  double WGT, SUM_WGT, BINSIZE;
  double AVG_z, AVG_x1, AVG_c, avg_z, avg_x1, avg_c ;
  double SUM_VAL[NLCPAR], SUM_ERR[NLCPAR];
  double SUM_SQERRINV[NLCPAR], SUM_SQRMS[NLCPAR] ;
  double VAL, ERR, RMS, dif, Dc, Dz, Dx1 ;

  double SUM_MBOFF[NLCPAR], SUM_MBSLOPE[NLCPAR];
  double DEBUG_LIST_DIF[3][50];
  int    DEBUG_LIST_INDX[3][50];
  double DEBUG_LIST_WGT[50];

  char fnam[] = "get_fitParBias" ;

  // ------------- BEGIN ----------------

  int LDMP = DUMPFLAG ;  

 START:

  // init output structure
  zero_FITPARBIAS(FITPARBIAS);

  // check "noBiasCor" option for this sample 
  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return(1); }

  // init output to crazy values.
  for(ipar=0; ipar < NLCPAR; ipar++ ) { FITPARBIAS->VAL[ipar]  = 666.0 ; } 

  // strip off local indices
  ia  = IBINFUN(a,  &simdata_bias.BININFO_SIM_ALPHA, 0, "" );
  ib  = IBINFUN(b,  &simdata_bias.BININFO_SIM_BETA,  0, "" );
  IZ  = IBINFUN(z,  &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z,         0, "" );
  IX1 = IBINFUN(x1, &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1],0,"");
  IC  = IBINFUN(c,  &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c],0,"");
  J1D = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][IZ][IX1][IC] ;

  if ( LDMP ) {      
    printf("\n") ;
    printf(" xxx ---------------------------------------------------- \n") ;
    if ( BADBIAS ) { printf("\t !!!!! BAD BIAS DETECTED !!!!! \n"); }
    printf(" xxx %s DUMP for CID=%s \n", fnam, cid );
    printf(" xxx  input: z=%.4f  mb,x1,c = %.2f, %.3f, %.3f \n",
	   z, mB, x1, c ); 

    printf(" xxx \t    IZ,IX1,IC=%d,%d,%d   NBIN(z,x1,c)=%d,%d,%d \n",
	   IZ, IX1, IC,
	   CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin,
	   CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin,
	   CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin  );

    printf(" xxx  IDSAMPLE=%d -> %s(%s) \n", IDSAMPLE,
	   SAMPLE_BIASCOR[IDSAMPLE].NAME_SURVEYGROUP, 
	   SAMPLE_BIASCOR[IDSAMPLE].NAME_FIELDGROUP );
    fflush(stdout);
    fflush(stdout);
  }

  if ( IZ  < 0 ) { return 0 ; }
  if ( IX1 < 0 ) { return 0 ; }
  if ( IC  < 0 ) { return 0 ; }

  // reset counters and weight for this a,b cell
  SUM_WGT = 0.0 ;
  NCELL_INTERP_TOT = NCELL_INTERP_USE  = NSUM_Cell = 0 ;
  
  for(ipar=0; ipar < NLCPAR; ipar++ ) {
    SUM_VAL[ipar]      = SUM_ERR[ipar]   = 0.0 ;
    SUM_SQERRINV[ipar] = SUM_SQRMS[ipar] = 0.0  ;
    SUM_MBOFF[ipar]    = SUM_MBSLOPE[ipar] = 0.0 ;
  }

  for(iz=0; iz < MAXBIN_z; iz++ ) { USEZ[iz] = 0; }
  
  // ------------------------------------------
  // determine bins to interpolate such that current value
  // is between biasCor nodes (i.e., no extrapolation allowed)
  
  IZMIN=IZMAX=IZ;  ICMIN=ICMAX=IC ;    IX1MIN=IX1MAX=IX1 ;
  
  AVG_z = CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D] ;
  if ( z >= AVG_z ) { IZMAX++ ; } else { IZMIN--; }
  if (IZMIN<0){IZMIN=0;}  if(IZMAX>=NBINz){IZMAX = NBINz-1;}

  
  AVG_x1 = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][J1D] ;
  if ( x1 >= AVG_x1 ) { IX1MAX++; } else { IX1MIN--; }
  if (IX1MIN<0){IX1MIN=0;}  if(IX1MAX>=NBINx1){IX1MAX = NBINx1-1;}

  
  AVG_c = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c][J1D] ;
  if ( c >= AVG_c ) { ICMAX++ ;} else { ICMIN--; }
  if (ICMIN<0){ICMIN=0;}  if(ICMAX>=NBINc){ICMAX = NBINc-1;}


  // determine max possible binsize in each dimension,
  // to use below for interpolation. 

  //  BINSIZE_z = BINSIZE_x1 = BINSIZE_c = 0.0 ;
  
  for(iz = IZMIN; iz <= IZMAX; iz++ ) {    
    for(ix1 = IX1MIN; ix1 <= IX1MAX; ix1++ ) {
      for(ic = ICMIN; ic <= ICMAX; ic++ ) {

	if( iz == IZ && ix1==IX1 && ic==IC ) { continue ; }
	
	j1d = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] ;
	NperCell = CELLINFO_BIASCOR[IDSAMPLE].NperCell[j1d] ;
	if ( NperCell < BIASCOR_MIN_PER_CELL  ) { continue; }

	avg_z  = CELLINFO_BIASCOR[IDSAMPLE].AVG_z[j1d] ;
	avg_x1 = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][j1d] ;
	avg_c  = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c][j1d] ;
	if ( avg_z  > 9000.0 ) { continue; }
	if ( avg_x1 > 9000.0 ) { continue; }
	if ( avg_c  > 9000.0 ) { continue; }

	BINSIZE = fabs(AVG_z - avg_z);
	if( BINSIZE > BINSIZE_z && iz!=IZ) { BINSIZE_z = BINSIZE; }

	BINSIZE = fabs(AVG_x1 - avg_x1);
	if( BINSIZE > BINSIZE_x1 && ix1!=IX1) { BINSIZE_x1 = BINSIZE; }

	BINSIZE = fabs(AVG_c - avg_c);
	if( BINSIZE > BINSIZE_c && ic!=IC) { BINSIZE_c = BINSIZE; }	
      }
    }
  }

  if ( LDMP ) {
    printf(" xxx  AVGCELL(z,x1,c) = %f, %f, %f \n",
	   AVG_z, AVG_x1, AVG_c);
    printf(" xxx  BINSIZE(z,x1,c) = %f, %f, %f \n",
	   BINSIZE_z, BINSIZE_x1, BINSIZE_c ); fflush(stdout);

    printf(" xxx \n");
    printf(" xxx iz,ix,ic   <z>   <x1>   <c>       "
	   "dmB    dx1    dc  Ncell WGT\n");  
    fflush(stdout);
  }

  // return if any binSize is not defined
  if ( BINSIZE_z  == 0.0 || BINSIZE_z  > 9000. ) { return 0; }
  if ( BINSIZE_x1 == 0.0 || BINSIZE_x1 > 9000. ) { return 0; }
  if ( BINSIZE_c  == 0.0 || BINSIZE_c  > 9000. ) { return 0; }

  
  // ----------------------------------------------
  // -------- start 3D loop over cells -------------
  
  for(iz = IZMIN; iz <= IZMAX; iz++ ) {    
    for(ix1 = IX1MIN; ix1 <= IX1MAX; ix1++ ) {
      for(ic = ICMIN; ic <= ICMAX; ic++ ) {

	NCELL_INTERP_TOT++ ;
	
	j1d = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][iz][ix1][ic] ;

	if ( j1d >= CELLINFO_BIASCOR[IDSAMPLE].NCELL || j1d < 0 ) {
	  sprintf(c1err,"Invalid j1d=%d ", j1d);
	  sprintf(c2err,"ia=%d ib=%d iz=%d ix1=%d ic=%d (cid=%s)",
		  ia, ib, iz, ix1, ic, cid);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	// require something in a cell to use the bias estimate
	NperCell = CELLINFO_BIASCOR[IDSAMPLE].NperCell[j1d] ;
	if ( NperCell < BIASCOR_MIN_PER_CELL  ) { continue; }

	// get distance between current data value and wgted-avg in bin

	dif = z - CELLINFO_BIASCOR[IDSAMPLE].AVG_z[j1d] ;
	Dz  = fabs(dif/BINSIZE_z) ;

	dif = x1 - CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][j1d] ;
	Dx1 = fabs(dif/BINSIZE_x1);
	    
	dif = c - CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c][j1d] ;
	Dc  = fabs(dif/BINSIZE_c);

	/*
	// only keep closest bins to interpolate
	if ( Dz  > 1.1 ) { continue ; }
	if ( Dx1 > 1.1 ) { continue ; }
	if ( Dc  > 1.1 ) { continue ; }	
	*/
	
	WGT   = (1.0 - Dz) * (1.0 - Dx1) * (1.0 - Dc);

	// prepare DEBUG_LISTs in case of abort
	DEBUG_LIST_DIF[0][NCELL_INTERP_USE]  = Dz;
	DEBUG_LIST_DIF[1][NCELL_INTERP_USE]  = Dx1;
	DEBUG_LIST_DIF[2][NCELL_INTERP_USE]  = Dc; 
	DEBUG_LIST_INDX[0][NCELL_INTERP_USE] = iz ;
	DEBUG_LIST_INDX[1][NCELL_INTERP_USE] = ix1 ;
	DEBUG_LIST_INDX[2][NCELL_INTERP_USE] = ic ;
  	DEBUG_LIST_WGT[NCELL_INTERP_USE]     = WGT ;
	
	if ( WGT < 0.0 ) { WGT = 0.0 ; }
		
	SUM_WGT   += WGT ;
	NSUM_Cell += NperCell ; // sum Nevt in nbr cells
	NCELL_INTERP_USE++ ;
	USEZ[iz]++ ;
	
	// now loop over the data fit params (mB,x1,c) to update bias sums
	for(ipar=0; ipar < NLCPAR; ipar++ ) { 
	  VAL     = simdata_bias.FITPARBIAS[IDSAMPLE][j1d].VAL[ipar] ; 
	  ERR     = simdata_bias.FITPARBIAS[IDSAMPLE][j1d].ERR[ipar] ; 
	  RMS     = simdata_bias.FITPARBIAS[IDSAMPLE][j1d].RMS[ipar] ; 

	  SUM_VAL[ipar]      += (WGT * VAL );
	  SUM_SQRMS[ipar]    += (WGT * RMS*RMS );
	  SUM_SQERRINV[ipar] += (WGT / (ERR*ERR) );
	  SUM_ERR[ipar]      += (WGT * ERR ) ;	 
	} // end ipar

	if ( LDMP ) {
	  printf( " xxx %2d,%2d,%2d  %5.4f %5.2f %6.3f  ",
		  iz, ix1, ic,
		  CELLINFO_BIASCOR[IDSAMPLE].AVG_z[j1d],
		  CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][j1d],
		  CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c ][j1d] ) ;
	  printf("  %6.3f %5.2f %6.3f  %3d %.2f ",
		 simdata_bias.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_mB],
		 simdata_bias.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_x1],
		 simdata_bias.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_c],
		 NperCell, WGT );
	  printf("\n");  fflush(stdout);

	}
	
      } // end ic
    } // end ix1
   } // end iz


  if ( LDMP ) {
    printf("\n") ;
    printf(" xxx\t SUM_WGT = %le   NCELL_INTERP_USE=%d (CUT=3)\n", 
	   SUM_WGT, NCELL_INTERP_USE );
    fflush(stdout);
  }

   // require enough cells for interpolation (July 2016)
  //   if ( NCELL_INTERP_USE < NCELL_INTERP_TOT ) { return(0); }
   if ( NCELL_INTERP_USE < 3 ) { return(0); } 

   // require both z-bins to be used.
   int ISKIP = 0 ;
   if ( IZMIN >=  0 ) { ISKIP = (USEZ[IZMIN]==0) ; }
   if ( IZMAX >=  0 ) { ISKIP = (USEZ[IZMAX]==0) ; }

   if ( LDMP ) {
     printf(" xxx IZMIN=%d IZMAX=%d  USEZ=%d,%d  ISKIP=%d\n",
	    IZMIN, IZMAX, USEZ[IZMIN], USEZ[IZMAX], ISKIP );
     fflush(stdout);
   }
   if ( ISKIP ) { return(0); }

   
   OK = 0 ;
   int NTMP = NSUM_Cell ;
   if ( NTMP >= BIASCOR_MINSUM ) { OK = 1; }
   if ( NTMP >= 2 &&  z > INPUTS.zmax ) { OK = 1; }

   if ( LDMP ) {
     printf(" xxx NSUM_Cell=%d  (cut=%d) OK=%d \n",
	    NSUM_Cell, BIASCOR_MINSUM, OK); fflush(stdout);
   }

   if ( OK == 0 ) { return(0); }   

  // - - - - - - - - - - 
  // make sure that SUMWGT > 0 
  if ( SUM_WGT <= 1.0E-9 ) {
    int icell ;
    printf("\n PRE-ABORT DUMP: \n");

    printf("  IZMIN/MAX=%d/%d   IX1MIN/MAX=%d/%d   ICMIN/MAX=%d/%d\n",	   
	   IZMIN,IZMAX,  IX1MIN,IX1MAX,   ICMIN, ICMAX);
    printf("  BINSIZE(z,x1,c) = %.3f, %.3f, %.3f \n",
	   BINSIZE_z, BINSIZE_x1, BINSIZE_c );
    
    printf("\n");
    printf("   cell   iz ix1 ic    Dz     Dx1     Dc     WGT \n");
    
    for(icell=0; icell < NCELL_INTERP_USE; icell++ ) {
      printf("    %3d   %2d %2d %2d   %6.3f %6.3f %6.3f   %.3f\n", icell
	     ,DEBUG_LIST_INDX[0][icell]
	     ,DEBUG_LIST_INDX[1][icell]
	     ,DEBUG_LIST_INDX[2][icell]
	     ,DEBUG_LIST_DIF[0][icell]
	     ,DEBUG_LIST_DIF[1][icell]
	     ,DEBUG_LIST_DIF[2][icell]
	     ,DEBUG_LIST_WGT[icell] );
    }
    sprintf(c1err,"SUM_WGT=%f for  CID=%s", SUM_WGT, cid) ;
    sprintf(c2err,"a=%.3f b=%.3f  z=%.4f  mB=%.4f  x1=%.4f  c=%.4f", 
	    a, b, z, mB, x1, c ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }      
  
  // update bias for each fit param (ipar=mB,x1,c)
  for(ipar=0; ipar < NLCPAR; ipar++ ) {
    FITPARBIAS->VAL[ipar]  = SUM_VAL[ipar] / SUM_WGT ;
    FITPARBIAS->ERR[ipar]  = sqrt((1.0/SUM_SQERRINV[ipar]) / SUM_WGT);
    FITPARBIAS->RMS[ipar]  = sqrt(SUM_SQRMS[ipar]/SUM_WGT);
    if ( fabs( FITPARBIAS->VAL[ipar]) > 20.0 ) { BADBIAS=1; }
  } 
  
  if ( LDMP ) {
    printf("\n") ;
    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      printf(" xxx bias(%2.2s) = %6.3f +-  %6.3f  SQERRINV=%le\n"
	     ,BIASCOR_NAME_LCFIT[ipar]
	     ,FITPARBIAS->VAL[ipar]
	     ,FITPARBIAS->ERR[ipar] 
	     ,SUM_SQERRINV[ipar]
	     );
    }
    fflush(stdout);
    //    debugexit(fnam);
  }

  if ( BADBIAS ) { LDMP=1; goto START ; }

  return(1) ;

} // end get_fitParBias

// ======================================================
int get_muCOVscale(char *cid, 
		   BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		   double *muCOVscale ) {

  // Created July 1 2016: 
  // Analog of get_fitParBias(), but for scale on muCOV.
  //  
  // inputs:
  //   cid          = name of SN (for error message)
  //   BIASCORLIST  = list of input values
  //
  // Output:
  //   muCovScale
  //
  //   Function returns 1 on success; 0 on failure
  //

  int ia, ib, iz, ic, IZ, IC, j1d ;

  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z  = BIASCORLIST->z ;
  double c  = BIASCORLIST->FITPAR[INDEX_c];
  double a  = BIASCORLIST->alpha ;
  double b  = BIASCORLIST->beta ;
  int IDSAMPLE = BIASCORLIST->idsample ;  

  // get bin info into local variables
  int  NBINz   = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.nbin ;
  int  NBINc   = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;

  double BINSIZE_z  
    = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.binSize ;
  double BINSIZE_c  
    = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].binSize ;
  
  double dif, Dz, Dc, WGT, SUM_WGT, SUM_muCOVscale ;
  double muCOVscale_local  ;
  char fnam[] = "get_muCOVscale" ;

  // -------------- BEGIN --------------

  // init output
  muCOVscale_local = 1.0 ;
  *muCOVscale = muCOVscale_local ;

  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOV)==0   ) { return(1); }

  // check "noBiasCor" option (Aug 20 2016)
  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return(1); }

  // get indices
  ia    = IBINFUN(a,  &simdata_bias.BININFO_SIM_ALPHA, 0, "" );
  ib    = IBINFUN(b,  &simdata_bias.BININFO_SIM_BETA,  0, "" );
  IZ    = IBINFUN(z,  &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z,  0, "" );

  // for color bin, use OPT=2 to allow values outside BININFO 
  // range to be pulled into edge bin.
  IC = IBINFUN(c,&CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c],2,""); 
  
  SUM_WGT = SUM_muCOVscale = 0.0 ;

  for(iz = IZ-1; iz <= IZ+1; iz++ ) {

    if ( iz < 0 )      { continue ; }
    if ( iz >= NBINz ) { continue ; }
          
    for(ic = IC-1; ic <= IC+1; ic++ ) {
      if ( ic < 0      ) { continue ; } 
      if ( ic >= NBINc ) { continue ; }  

      j1d = CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][iz][0][ic] ;

      if ( j1d >= CELLINFO_MUCOVSCALE[IDSAMPLE].NCELL || j1d < 0 ) {
	  sprintf(c1err,"Invalid j1d=%d ", j1d);
	  sprintf(c2err,"ia=%d ib=%d iz=%d ic=%d",
		  ia, ib, iz, ic);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

      // get distance between current data value and bin-center
  
      dif = z - CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z[j1d];
      Dz  = fabs(dif/BINSIZE_z) ;
	      
      dif = c - CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c][j1d];
      Dc  = fabs(dif/BINSIZE_c);

      // only keep closest bins to interpolate
      if ( Dz  > 0.99 ) { continue ; }
      if ( Dc  > 0.99 ) { continue ; }
	
      WGT        = (1.0 - Dz) * (1.0 - Dc);
      SUM_WGT   += WGT ;

      SUM_muCOVscale += ( WGT * simdata_bias.MUCOVSCALE[IDSAMPLE][j1d] ) ;
    } // ic
   } // iz
  
  
  if( SUM_WGT > 0.01 ) {
    muCOVscale_local = SUM_muCOVscale / SUM_WGT ;
  }

  // load output arg.
  *muCOVscale = muCOVscale_local ;

  return(1);

} // end get_muCOVscale


// ======================================================
void setup_CELLINFO_biasCor(int IDSAMPLE) {

  // Created Aug 23 2016
  // Setup 5D biasCor cells for input IDSAMPLE

  int NSAMPLE = NSAMPLE_BIASCOR ;
  int INDX;
  char fnam[] = "setup_CELLINFO_biasCor" ;

  // ------------ BEGIN -----------

  // on first call, allocate CELLINFO structures
  if ( IDSAMPLE == 0 ) {

    printf("\n# ============== %s ================= \n", fnam);
    int MEM = NSAMPLE * sizeof(CELLINFO_DEF);
    CELLINFO_BIASCOR    = (CELLINFO_DEF*) malloc ( MEM );
    CELLINFO_MUCOVSCALE = (CELLINFO_DEF*) malloc ( MEM );

    // setup bining for SIMalpha,beta; note storage in different structure
    // since alpha,beta binning is fixed for all IDSAMPLEs
    printf("\n\t\t Global Alpha,Beta Binning \n"); fflush(stdout);
    setup_BININFO_biasCor(IDSAMPLE, 100*INDEX_x1, MAXBIN_BIASCOR_ALPHA,
			  &simdata_bias.BININFO_SIM_ALPHA );
    setup_BININFO_biasCor(IDSAMPLE, 100*INDEX_c, MAXBIN_BIASCOR_BETA,
			  &simdata_bias.BININFO_SIM_BETA );
    
    
    // allocate sample-dependent FITPARBIAS & MUCOVSCALE
    int MEMBIAS = NSAMPLE * sizeof(FITPARBIAS_DEF*) ;
    int MEMCOV  = NSAMPLE * sizeof(double *) ;
    simdata_bias.FITPARBIAS = (FITPARBIAS_DEF**) malloc ( MEMBIAS );
    simdata_bias.MUCOVSCALE = (double **       ) malloc ( MEMCOV );    
  }


  CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin  = 0 ; 
  CELLINFO_BIASCOR[IDSAMPLE].NCELL           = 0 ;

  if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR==0 ) { return ; }

  printf("\n\t Setup biasCor bins for SAMPLE = %s \n",
	 SAMPLE_BIASCOR[IDSAMPLE].NAME ); fflush(stdout);

  // setup binning for redshift
  setup_BININFO_biasCor(IDSAMPLE, -1, MAXBIN_z, 
			&CELLINFO_BIASCOR[IDSAMPLE].BININFO_z ); 

  // setup binning for mB,x1,c
  for(INDX=1; INDX<NLCPAR ; INDX++ ) { // skip mB
    setup_BININFO_biasCor(IDSAMPLE, INDX, MAXBIN_BIASCOR_FITPAR,
			  &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDX] );
  }
  fflush(stdout);

  // make map between 5D indices and 1D index
  set_MAPCELL_biasCor(IDSAMPLE);

  return ;

} // end setup_CELLINFO_biasCor

// ======================================================
void setup_BININFO_biasCor(int IDSAMPLE, int ipar_LCFIT, int MAXBIN, 
			   BININFO_DEF *BININFO) {

  // 
  // ipar_LCFIT = -1    --> redshift
  // ipar_LCFIT = 0,1,2 --> mB, x1, c   (fitted)
  // ipar_LCFIT = 100   --> SIM_ALPHA   (generated)
  // ipar_LCFIT = 200   --> SIM_BETA    (generated)
  //
  // Output: *BININFO
  // 
  // Oct 8 2016: adjust LO_MAX calc to allow for TEXT truncation
  //
  // Apr 30 2017: chnage "nbin>=MAXBIN" -> "nbin>MAXBIN" for abort trap.

  double val_lo, val_hi, VAL_MIN, VAL_MAX, VAL_BIN;
  int nbin ;
  char NAME[20];
  char fnam[] = "setup_BININFO_biasCor" ;

  // ------------- BEGIN ----------------

  if ( ipar_LCFIT >= 0 && ipar_LCFIT < 10 ) {
    // mb, x1, c
    VAL_MIN  = BIASCOR_MINVAL_LCFIT[ipar_LCFIT] ; 
    VAL_MAX  = BIASCOR_MAXVAL_LCFIT[ipar_LCFIT] ;
    VAL_BIN  = SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_FITPAR[ipar_LCFIT] ;
    sprintf(NAME,"%s", BIASCOR_NAME_LCFIT[ipar_LCFIT] );
  }
  else if ( ipar_LCFIT < 0 ) {  
    // redshift  
    VAL_MIN = SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[0];
    VAL_MAX = SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[1];
    VAL_BIN = SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT ;
    sprintf(NAME,"z");

    // add one more z bin on the high-z side to get better
    // interpolation in last z-bin
    VAL_MAX += VAL_BIN ;    
  }
  else if ( ipar_LCFIT == 100*INDEX_x1 ) {
    // get info for SIMalpha binning
    get_BININFO_biasCor_alphabeta("SIMalpha", &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIMalpha");
  }
  else if ( ipar_LCFIT == 100*INDEX_c ) {
    // get info for SIMbeta binning
    get_BININFO_biasCor_alphabeta("SIMbeta", &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIMbeta");
  }
  else {
    sprintf(c1err,"Invalid ipar_LCFIT=%d", ipar_LCFIT );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" );    
  }


  // -----------------------------
  nbin=0 ;
  double LO_MIN = VAL_MIN;
  double LO_MAX = VAL_MAX - 0.99*VAL_BIN ;
  for(val_lo = LO_MIN; val_lo <= LO_MAX; val_lo += VAL_BIN ) {   
    val_hi = val_lo + VAL_BIN ;

    if ( nbin < MAXBIN ) {
      BININFO->lo[nbin]  = val_lo ;
      BININFO->hi[nbin]  = val_hi ;
      BININFO->avg[nbin] = 0.5 * (val_lo + val_hi) ;
    }
    nbin++  ;
    if ( VAL_BIN == 0.0 ) { val_lo += 0.00001; }
  } // end val_lo loop

  BININFO->nbin = nbin ;
  BININFO->binSize = VAL_BIN ;
  sprintf(BININFO->varName, "%s", NAME);

  printf("  %s: define %2d %2s bins from %6.3f to %6.3f \n",
	 fnam, nbin, NAME, BININFO->lo[0], BININFO->hi[nbin-1] );


  if ( nbin > MAXBIN || nbin >= MXpar ) {
    sprintf(c1err,"nbin = %d exceeds MAXBIN=%d or MXpar=%d", 
	    nbin, MAXBIN, MXpar );
    sprintf(c2err,"Check %s binning.", NAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
  }

  fflush(stdout);

  return ;

} // end setup_BININFO_biasCor


// ==================================================
void  get_BININFO_biasCor_alphabeta(char *varName, 
				    double *VAL_MIN, double *VAL_MAX, 
				    double *VAL_BIN ) {

  // Return binning for SIMalpha or SIMbeta.
  // This is tricky because there are two ways to define
  // alpha,beta bins:
  //
  // 1) alpha,beta are defined on a grid
  //      --> use native grid from sim file
  // 2) alpha,beta are continuous
  //      --> force 2x2 grid
  //
  // Apr 28 2017: abort if NVAL > NBMAX

  double valmin_loc, valmax_loc, valbin_loc;
  double val, val_last, val1st, val2nd ;
  float  *ptrVal_f ;
  int    irow, isort, NVAL, NBMAX ;
  int    NROW = simdata_bias.NROW ;
  int    ISGRID;
  char fnam[] = "get_BININFO_biasCor_alphabeta" ;

  // ------------- BEGIN --------------

  *VAL_MIN = *VAL_MAX = *VAL_BIN = 0.0 ;

  valmin_loc = +1.0E8;
  valmax_loc = -1.0E8 ;
  valbin_loc =  0.0 ;

  if ( strcmp(varName,"SIMalpha") == 0 ) 
    { ptrVal_f = simdata_bias.SIM_ALPHA ; NBMAX=MXa; }
  else if ( strcmp(varName,"SIMbeta") == 0 )  
    { ptrVal_f = simdata_bias.SIM_BETA ;  NBMAX=MXb; }
  else {
    sprintf(c1err,"Invalid varName = '%s' ", varName);
    errmsg(SEV_FATAL, 0, fnam, c1err, "");  
  }

  int ORDER_SORT  = +1 ; // increasing order
  int *INDEX_SORT = (int*)malloc( NROW * sizeof(int) );
  //  sortDouble( NROW, ptrVal, ORDER_SORT, INDEX_SORT ) ;
  sortFloat( NROW, ptrVal_f, ORDER_SORT, INDEX_SORT ) ;

  val1st = val2nd = -99999. ;
  NVAL = 0;  val_last = -9999. ;
  for ( irow=0; irow < NROW; irow++ ) {
    isort = INDEX_SORT[irow];
    val   = (double)ptrVal_f[isort];

    if ( simdata_bias.SIM_NONIA_INDEX[isort] != 0 ) { continue ; }

    if ( val > val_last ) { 
      NVAL++ ; 
      // store 1st and send value to get binsize in case of grid
      if ( NVAL == 1 ) { val1st = val2nd = val ; }
      if ( NVAL == 2 ) { val2nd = val ; }
    }
    val_last = val ;
  }


  // add error check on number of a,b bins (Apr 28 2017)
  if ( NVAL > NBMAX ) {
    sprintf(c1err,"%d %s values exceeds bound of %d",
	    NVAL, varName, NBMAX);
    sprintf(c2err,"Check biasCor sample.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
  }


  // Dec 2016: check option to use only the boundardy bins
  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_2x2ab) > 0 ) {
    val2nd = val_last;  NVAL=2;
  }


  // ---------------------------
  if ( NVAL < 10 ) { 
    ISGRID = 1 ; 
    valbin_loc = val2nd - val1st ;
    valmin_loc = val1st   - 0.5*valbin_loc ;
    valmax_loc = val_last + 0.5*valbin_loc ;
    
  }
  else {
    ISGRID = 0 ;  // continuous
    isort = INDEX_SORT[0];        valmin_loc = (double)ptrVal_f[isort];
    isort = INDEX_SORT[NROW-1];   valmax_loc = (double)ptrVal_f[isort];

    // hard-wire 3 bins; maybe later add user parameter
    valbin_loc = (valmax_loc - valmin_loc)/3.0 ;  
  }


  /*
  printf("  %s : \n\t Found %d distinct %s values --> %s \n",
	 fnam, NVAL, varName, GRIDDEF[ISGRID] );
  printf("\t %s Range = %.4f to %.4f with binSize=%.4f \n",
	 varName, valmin_loc, valmax_loc, valbin_loc );
  fflush(stdout);
  */

  free(INDEX_SORT); // free memory
  
  // --------------------------
  // load function args
  *VAL_MIN = valmin_loc ;
  *VAL_MAX = valmax_loc ;
  *VAL_BIN = valbin_loc ;

  return ;

} // end get_BININFO_biasCor_alphabeta

// ================================================
void get_muBias(char *NAME,
		BIASCORLIST_DEF *BIASCORLIST, 
		FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MXb],
		double         (*MUCOVSCALE_ABGRID)[MXb],
		INTERPWGT_ALPHABETA *INTERPWGT,  
		double *fitParBias,
		double *muBias, double *muBiasErr, double *muCOVscale ) { 

  
  // Created Jan 2016
  // Inputs:
  //   NAME        = name of SN (for err msg only)
  //   BIASCORLIST = z,a,b,mB,x1,c 
  //   FITPARBIAS  = pre-computed bias for {mB,x1,c} on grid of {a,b}
  //   INTERPWGT   = wgt at each ia,ib
  //
  // Ouptuts:
  // *fitParBias  = bias on mB,x1,c, interpolated over alpha & beta
  //  muBias      = bias on distance
  //  muBiasErr   = error on above (based on biasCor sim stats)
  //  muCOVscale  = scale bias to apply to muErr 
  //

  double alpha = BIASCORLIST->alpha ;
  double beta  = BIASCORLIST->beta  ;
  double z     = BIASCORLIST->z  ;
  double mB    = BIASCORLIST->FITPAR[INDEX_mB] ;
  double x1    = BIASCORLIST->FITPAR[INDEX_x1] ;
  double c     = BIASCORLIST->FITPAR[INDEX_c] ;

  double muBias_local     = 0.0 ;
  double muBiasErr_local  = 0.0 ;
  double muCOVscale_local = 0.0 ;

  double VAL, ERR, SQERR ;
  double WGTab, WGTpar, WGTpar_SUM[NLCPAR];
  double biasVal[NLCPAR], biasErr[NLCPAR] ;
  double MUCOEF[NLCPAR];

  int  ia, ib, ipar ;
  char fnam[] = "get_muBias";

  // --------------- BEGIN ------------

  for(ipar=0; ipar < NLCPAR; ipar++ )  { 
    biasVal[ipar] = biasErr[ipar] = 0.0 ;
    WGTpar_SUM[ipar] = 0.0 ;
  }

  MUCOEF[INDEX_mB] = +1.0 ;
  MUCOEF[INDEX_x1] = +alpha ;
  MUCOEF[INDEX_c ] = -beta ;

  // interpolate bias on grid of alpha & beta.
  // The WGT at each corner was computed earlier.
  for(ia=INTERPWGT->ia_min; ia <= INTERPWGT->ia_max; ia++ ) {
    for(ib=INTERPWGT->ib_min; ib <= INTERPWGT->ib_max; ib++ ) {
      
      WGTab = INTERPWGT->WGT[ia][ib] ; // weight for this alpha,beta
      if ( WGTab < -1.0E-12 || WGTab > 1.0000001 ) {
	sprintf(c1err,"Undefined WGTab=%le for CID=%s", WGTab, NAME);
	sprintf(c2err,"ia=%d ib=%d  a=%.3f b=%.2f  z=%.3f",
		ia, ib, alpha, beta, z );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
      }
      
      for(ipar=0; ipar < NLCPAR; ipar++ ) {
	VAL = FITPARBIAS_ABGRID[ia][ib].VAL[ipar];
	ERR = FITPARBIAS_ABGRID[ia][ib].ERR[ipar];
	if ( VAL > 600.0 || isnan(ERR) ) {
	  sprintf(c1err,"Undefined %s-bias = %.3f +- %.3f for CID=%s", 
		  BIASCOR_NAME_LCFIT[ipar], VAL, ERR, NAME);
	  sprintf(c2err,"ia=%d ib=%d  a=%.3f b=%.2f",
		  ia, ib, alpha, beta );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
	}

	WGTpar = WGTab ;
	//	WGTpar = WGTab/ERR ;
	WGTpar_SUM[ipar] += WGTpar;
	biasVal[ipar] += ( WGTpar * VAL );
	biasErr[ipar] += ( WGTpar * ERR ); 
      } // end ipar

      // now the COV scale (Jul 1 2016)
      VAL = MUCOVSCALE_ABGRID[ia][ib];
      muCOVscale_local += ( WGTab * VAL ) ;

    } // end ib
  } // end ia


  // Note that sum(WGT) is already normalized to 1.000,
  // so no need to divide by SUMWGT here.

  muBias_local = SQERR = 0.0 ;
  
  for(ipar=0; ipar < NLCPAR; ipar++ )  { 

    biasVal[ipar] /= WGTpar_SUM[ipar] ;
    biasErr[ipar] /= WGTpar_SUM[ipar] ;

    fitParBias[ipar]   = biasVal[ipar] ;     //store output array
    muBias_local += ( MUCOEF[ipar] * biasVal[ipar] ) ;
    ERR           = ( MUCOEF[ipar] * biasErr[ipar] ) ;
    SQERR        += (ERR * ERR);
  }    

  // check for NaN on SQERR
  if ( isnan(SQERR) ) {
    printf("\n PRE-ABORT DUMP in %s : \n", fnam );
    for(ipar=0; ipar < NLCPAR; ipar++ )  { 
      printf("\t %2s = %.3f: bias = %.3f +- %.3f \n",
	     BIASCOR_NAME_LCFIT[ipar], BIASCORLIST->FITPAR[ipar],
	     biasVal[ipar], biasErr[ipar] );
    }    
    sprintf(c1err,"isnan(SQERR) for CID = %s", NAME);
    sprintf(c2err,"z=%.3f", z );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
  }
  
  muBiasErr_local = sqrt(SQERR);

  // ---------------------------------

  if ( fabs(muBias_local) > 5.0 ) {
    printf(" xxx ------------------------------- \n");
    printf(" xxx PRE-ABORT DUMP: \n");
    printf(" xxx alpha=%f  beta=%f \n", alpha, beta);

    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      printf(" xxx bias(%2s) = %7.3f +- %.3f  \n"
	     ,BIASCOR_NAME_LCFIT[ipar]
	     ,biasVal[ipar], biasErr[ipar]
	     );
    }

    sprintf(c1err,"Crazy muBias=%f for CID = %s",  muBias_local, NAME);
    sprintf(c2err,"z=%.3f  mB=%.3f  x1=%.3f  c=%.4f", z, mB, x1, c );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  *muBias     = muBias_local ;
  *muBiasErr  = muBiasErr_local ;
  *muCOVscale = muCOVscale_local ;

  return ;

} // end get_muBias


// ==================================================
double get_magoff_host(double *hostPar, int n) {

  // return SN mag offset based on host mass.
  // Inputs:
  //  hostPar   [see below]
  //  n     = index for data[n]
  //
  // Sep 13 2017: add gamma0==0 test to return magoff=0.
  //              Allows working with fixed gamma0 != 0
  //

  // strip off host properties into local variables
  // These are floated (or fixed) in the fit.
  double gamma0      = hostPar[0]; // magSPlit at z=0
  double gamma1      = hostPar[1]; // dgamma/dz
  double logmass_cen = hostPar[2]; // Prob=0.5 at this logmass
  double logmass_tau = hostPar[3]; // tau param in transition function

  double z, gamma, FermiFun, arg ;
  double magoff = 0.0 ;
  int    IVAR       = rawdata.ICUTWIN_GAMMA ;

  //  printf(" xxx IVAR = %d  n=%d\n", IVAR, n) ;
  double xval_gamma = rawdata.CUTVAL[IVAR][n] ; // value to correlate with HR step
  char   fnam[] = "get_magoff_host" ;

  // ------------ BEGIN ---------------


  if ( INPUTS.USE_GAMMA0 == 0 ) { return(magoff); }

  z        = data[n].zhd ;
  gamma    = gamma0 + z*gamma1 ;
  arg      = -( xval_gamma - logmass_cen ) / logmass_tau ;
  FermiFun = 1.0 / ( 1.0 + exp(arg) ) ;  // 0 to 1
  magoff   = gamma * ( FermiFun - 0.5 ) ;

  /*
  printf(" xxx xval=%.2f cen=%.1f tau=%.3f Fun=%.3f magoff=%.3f \n",
	 xval_gamma, logmass_cen, logmass_tau, FermiFun, magoff);
  debugexit(fnam); 
  */

  return(magoff);

} // end get_magoff_host

// ==================================================
void prepare_CCprior(char *simFile) {

  // Created Jan 2016
  // If simFile is not null, then prepare CC prior for BEAMS-like
  // usage in fit likelihood.
  //
  // Jun 20 2018: abort on 1D biasCor.
  
  char fnam[] = "prepare_CCprior" ;
  int idsample ;
  int NSAMPLE  ;

  // ------------- BEGIN -------------

  simdata_ccprior.USE = 0 ;

  if ( IGNOREFILE(simFile) ) { 
    INPUTS.ipar[IPAR_scalePCC] = 0 ; // make sure scalePCC is not floated
    return ; 
  }

  simdata_ccprior.USE = 1 ;

  
  if ( simdata_bias.DOCOR_1D ) {
    sprintf(c1err,"Cannot use 1D biasCor with CC likelihood.");
    sprintf(c2err,"Either remove CC term, or use 5D biasCor.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  if ( simdata_ccprior.USEH11 ) { 
    sprintf(BANNER,"%s: use CC mu-vs-z prior from Hlozek 2011", fnam);
    print_banner(BANNER);    
    return ;
  }

  sprintf(BANNER,"%s: read simulated CC mu-prior vs. z", fnam);
  print_banner(BANNER);


  // allocate memory and read simFile
  read_simFile_CCprior(simFile, &simdata_ccprior.SIMFILE_INFO_ALL);

  // transfer list of true CC events to another INFO_CC structure,
  // so that large-original INFO_ALL struct can be deleted from memory.
  store_simFile_CCprior(&simdata_ccprior.SIMFILE_INFO_ALL,
			&simdata_ccprior.SIMFILE_INFO_CC);

  // free memory 
  malloc_simdata_ccprior(-1, 0, &simdata_ccprior.SIMFILE_INFO_ALL);

  // set z-bin index for each event to speed things up in fitting
  setup_zbins_CCprior(&simdata_ccprior.SIMFILE_INFO_CC,
		      &simdata_ccprior.MUZMAP.ZBIN );

  NSAMPLE = NSAMPLE_BIASCOR ;
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    setup_MUZMAP_CCprior(idsample, &simdata_ccprior.SIMFILE_INFO_CC,
			 &simdata_ccprior.MUZMAP );
  }

  
  // ----------------------------------------------
  printf("\n Finished preparing CC prior. \n");
  fflush(stdout);

  return ;

} // end prepare_CCprior


void  read_simFile_CCprior(char *simFile, SIMFILE_INFO_DEF *SIMFILE_INFO) {

  // Created Jan 2016
  // Read simFile and load SIMFILE_INFO struct.

  char fnam[] = "read_simFile_CCprior" ;
  char str_z[40], str_zerr[40], vartmp[60], STRINGOPT[40];
  int  vbose=1;
  int  IFILETYPE, NVAR_ORIG, LEN, NROW, icut, ivar ;
  int  USE_FIELDGROUP  = INPUTS.use_fieldGroup_biasCor ;

  // -------------------- BEGIN --------------------

  LEN = SNTABLE_NEVT(simFile,TABLENAME_FITRES); LEN+=10;
  SIMFILE_INFO->LEN_MALLOC = LEN ; ;

  printf("\t Found %d lines in %s. \n", LEN, simFile) ;
  malloc_simdata_ccprior(+1, LEN, SIMFILE_INFO);

  IFILETYPE = TABLEFILE_OPEN(simFile,"read");
  NVAR_ORIG = SNTABLE_READPREP(IFILETYPE,"FITRES");
  
  sprintf(vartmp, "CID:C*%d  CCID:C*%d", MXCHAR_CCID, MXCHAR_CCID); 
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->name, LEN ,vbose) ;

  if( USE_FIELDGROUP ) {
    sprintf(vartmp, "FIELD:C*%d", MXCHAR_CCID); 
    SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->field, LEN ,vbose) ;
  }

  get_zString(str_z,str_zerr,"D");
  SNTABLE_READPREP_VARDEF(str_z, SIMFILE_INFO->z, LEN ,vbose);
  
  sprintf(vartmp,"c:D" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->c, LEN ,vbose);

  sprintf(vartmp,"x1:D" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->x1, LEN ,vbose);

  sprintf(vartmp,"x0:D" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->x0, LEN ,vbose);

  sprintf(vartmp,"SIM_DLMAG:D" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->SIM_DLMAG,LEN,vbose);
 
  sprintf(vartmp,"SIM_TYPE_INDEX:I" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->SIM_TYPE_INDEX,LEN,vbose);

  sprintf(vartmp,"SIM_NONIA_INDEX:I SIM_TEMPLATE_INDEX:I" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->SIM_NONIA_INDEX,LEN,vbose);

  sprintf(vartmp,"IDSURVEY:I" );
  SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->idsurvey,LEN,vbose);

  if ( FOUNDKEY_OPT_PHOTOZ ) {
    sprintf(vartmp,"OPT_PHOTOZ:I" );
    SNTABLE_READPREP_VARDEF(vartmp, SIMFILE_INFO->opt_photoz,LEN,vbose);
  }

  //read CUTWIN variables to apply same cuts as for data.
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
    sprintf(vartmp, "%s", INPUTS.CUTWIN_NAME[icut]);
    if ( usesim_CUTWIN(vartmp) == 0 ) { continue ; }
    ivar=SNTABLE_READPREP_VARDEF(vartmp,SIMFILE_INFO->CUTVAL[icut],LEN,vbose);
    SIMFILE_INFO->DOFLAG_CUTWIN[icut] = set_DOFLAG_CUTWIN(ivar,icut,0);
  }
  
 
  // read entire file and load arrays
  NROW = SNTABLE_READ_EXEC();
  SIMFILE_INFO->NROW = NROW ;

  // close file
  // xxx obsolete   TABLEFILE_CLOSE(simFile); 

  // ---------------------------------------
  // set IDSAMPLE for each CCprior event
  int   IDSURVEY, IDSAMPLE, OPT_PHOTOZ=0, i ;
  double CUTVAL_LOCAL[MXCUTWIN] ;
  char *NAME, FIELD[MXCHAR_CCID], FIELDGROUP[100] ;
  char  SURVEY[MXCHAR_CCID], SURVEYGROUP[100];
  int   NOBIASCOR =  ( (INPUTS.opt_biasCor & MASK_BIASCOR_SAMPLE)==0 );

  if ( NOBIASCOR ) {
    NSAMPLE_BIASCOR = 1;
    sprintf(SAMPLE_BIASCOR[0].NAME_FIELDGROUP,   "NONE");
    sprintf(SAMPLE_BIASCOR[0].NAME,              "ALL");    
  }

  sprintf(FIELD,"NONE");  sprintf(SURVEY,"NONE");

  for (i=0; i < NROW ; i++ ) {    
    
    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
      { CUTVAL_LOCAL[icut] = SIMFILE_INFO->CUTVAL[icut][i] ; }
    if ( apply_CUTWIN(2, SIMFILE_INFO->DOFLAG_CUTWIN, CUTVAL_LOCAL) == 0 )  
      { continue  ; }

    IDSURVEY   = SIMFILE_INFO->idsurvey[i] ;
    if ( FOUNDKEY_OPT_PHOTOZ ) 
      { OPT_PHOTOZ = SIMFILE_INFO->opt_photoz[i] ; }

    NAME     = SIMFILE_INFO->name[i] ;
    if ( USE_FIELDGROUP ) 
      { sprintf(FIELD,"%s", SIMFILE_INFO->field[i] ); }
    
    if ( NOBIASCOR ) {      
      IDSAMPLE = 0 ;   FIELDGROUP[0] = 0 ;
    }
    else {      
      match_fieldGroup (NAME, FIELD,   
			FIELDGROUP, STRINGOPT);  // (O)
      match_surveyGroup(NAME, IDSURVEY, 
			SURVEYGROUP, STRINGOPT); // (O)
      IDSAMPLE = get_IDSAMPLE(IDSURVEY,OPT_PHOTOZ,FIELDGROUP,SURVEYGROUP,0) ;
    }

    SIMFILE_INFO->idsample[i] = IDSAMPLE ;
    SAMPLE_BIASCOR[IDSAMPLE].NCCPRIOR++ ;    
  }
  dump_SAMPLE_INFO("CCPRIOR");

  return ;

} // end  read_simFile_CCprior


void  store_simFile_CCprior(SIMFILE_INFO_DEF *INFO_ALL, 
			    SIMFILE_INFO_DEF *INFO_CC) {

  // Created Jan 2016
  // 
  // Read INFO_ALL to count number of true CC events.
  // Then allocate INFO_CC memory.
  // Finally, copy CC-subset from INFO_ALL into INFO_CC array
  //

  int  NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;
  int  NSAMPLE = NSAMPLE_BIASCOR ;

  int  isn, icc, icut, NCC_CUTS_ALL, NCC_BIASCOR_ALL, SIM_TYPE ;
  int  NCC_CUTS_SAMPLE[MXNUM_SAMPLE], NCC_BIASCOR_SAMPLE[MXNUM_SAMPLE];
  int  ia, ib, ipar, istat_bias, idsample ;
  int *KEEP ;
  double z, x0, mB, x1, c, CUTVAL_LOCAL[MXCUTWIN] ;
  char fnam[] = "store_simFile_CCprior" ;

  // ------------- BEGIN -------------

  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    NCC_CUTS_SAMPLE[idsample] = NCC_BIASCOR_SAMPLE[idsample] = 0;
  }

  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
    { INFO_CC->DOFLAG_CUTWIN[icut] = INFO_ALL->DOFLAG_CUTWIN[icut]; }


  NCC_CUTS_ALL = NCC_BIASCOR_ALL = 0 ;
  KEEP = (int*) malloc ( INFO_ALL->NROW * sizeof(int) );

  // count number of true CC with cuts

  for(isn=0; isn < INFO_ALL->NROW; isn++ ) {
    KEEP[isn] = 0;

    SIM_TYPE = INFO_ALL->SIM_TYPE_INDEX[isn] ;
    if ( SIM_TYPE == 1 ) { continue ; } // exclude true SNIa

    z  = INFO_ALL->z[isn];
    c  = INFO_ALL->c[isn];
    x1 = INFO_ALL->x1[isn];
    idsample = INFO_ALL->idsample[isn];

    if ( z  < INPUTS.zmin  ) { continue ; }
    if ( z  > INPUTS.zmax  ) { continue ; }
    if ( x1 < INPUTS.x1min ) { continue ; }
    if ( x1 > INPUTS.x1max ) { continue ; }
    if ( c  < INPUTS.cmin  ) { continue ; }
    if ( c  > INPUTS.cmax  ) { continue ; }

    // apply same CUTWIN cuts as for data
    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
      { CUTVAL_LOCAL[icut] = INFO_ALL->CUTVAL[icut][isn] ; }
    if ( apply_CUTWIN(2, INFO_ALL->DOFLAG_CUTWIN, CUTVAL_LOCAL) == 0 )  
      { continue ; }

    NCC_CUTS_ALL++ ;        
    NCC_CUTS_SAMPLE[idsample]++ ; 

    // make sure bias correction is valid
    istat_bias = 1;
    if ( simdata_bias.NROW > 0 )  { istat_bias = storeBias_CCprior(isn);  }
    if ( istat_bias == 0 ) { continue ; }

    NCC_BIASCOR_ALL++ ;
    NCC_BIASCOR_SAMPLE[idsample]++ ;
    KEEP[isn] = 1;
  }

  // ------------------------------------------------------

  printf("\n");
  printf("   Found %d true CC after cuts. \n",
	 NCC_CUTS_ALL ) ;
  printf("   Found %d true CC after cuts with valid biasCor. \n",
	 NCC_BIASCOR_ALL ) ;

  // transfer CC events passing cuts to a much smaller array;
  // then delete the big array.

  malloc_simdata_ccprior(+1, NCC_BIASCOR_ALL, INFO_CC);

  icc = 0 ;
  for(isn=0; isn < INFO_ALL->NROW; isn++ ) {
    if ( KEEP[isn] == 0 ) { continue ; }
    
    x0 = INFO_ALL->x0[isn] ;
    mB = -2.5*log10(x0);
    INFO_CC->z[icc]   = INFO_ALL->z[isn];
    INFO_CC->x0[icc]  = x0 ;
    INFO_CC->mB[icc]  = mB ;
    INFO_CC->x1[icc]  = INFO_ALL->x1[isn];
    INFO_CC->c[icc]   = INFO_ALL->c[isn];
    INFO_CC->SIM_TYPE_INDEX[icc]  = INFO_ALL->SIM_TYPE_INDEX[isn];
    INFO_CC->SIM_NONIA_INDEX[icc] = INFO_ALL->SIM_NONIA_INDEX[isn];

    INFO_CC->idsurvey[icc] = INFO_ALL->idsurvey[isn];
    INFO_CC->idsample[icc] = INFO_ALL->idsample[isn];

    for(ia = 0; ia < NBINa; ia++ ) {
      for(ib = 0; ib < NBINb; ib++ ) {
	for(ipar=0; ipar < NLCPAR; ipar++ ) { 
	  INFO_CC->FITPARBIAS_ALPHABETA[icc][ia][ib].VAL[ipar] =
	  INFO_ALL->FITPARBIAS_ALPHABETA[isn][ia][ib].VAL[ipar] ;
	  INFO_CC->FITPARBIAS_ALPHABETA[icc][ia][ib].ERR[ipar] =
	    INFO_ALL->FITPARBIAS_ALPHABETA[isn][ia][ib].ERR[ipar] ;
	  INFO_CC->FITPARBIAS_ALPHABETA[icc][ia][ib].RMS[ipar] =
	    INFO_ALL->FITPARBIAS_ALPHABETA[isn][ia][ib].RMS[ipar] ;
	} //end ipar
      } // end ib
    }  // end ia
    
    icc++ ;  
  }
  INFO_CC->NROW = NCC_BIASCOR_ALL ;

  if ( icc != NCC_BIASCOR_ALL ) {
    sprintf(c1err,"icc=%d but NCC=%d", icc, NCC_BIASCOR_ALL );
    sprintf(c2err,"but they should be equal");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

 
  free(KEEP);

  return ;
  
} // end store_simFile_CCprior

void setup_zbins_CCprior(SIMFILE_INFO_DEF *INFO_CC, BININFO_DEF *ZBIN) {

  // Setup z-grid for CCprior, and label each data point with
  // iz index to save time in fitting. Note that this z binning
  // can be different from the user z-binning.
  // Mar 22 2017: Use z0,z1 to define z-range
  // Jun 04 2018: fix bug computing DZBIN_CCPRIOR

  char MSGERR[100];
  double z, zlo, zhi, z0, z1, dnbz, zrange ;
  int  SIM_TYPE, nbz;
  int  n_perbin_Ibc[MAXBIN_z], n_perbin_II[MAXBIN_z], n_perbin_other[MAXBIN_z];
  
  double DZBIN_CCPRIOR = 0.10 ;

  char fnam[] = "setup_zbins_CCprior" ;

  // ------------- BEGIN -------------

  // find exact DZBIN that fits in user z range
  zrange = INPUTS.zmax - INPUTS.zmin ;
  dnbz   = zrange/DZBIN_CCPRIOR ;
  nbz    = (int)dnbz;
  if ( dnbz > (double)nbz ) { nbz++ ; }
  DZBIN_CCPRIOR = zrange/(double)nbz ;


  // make uniform z-grid, although a non-uniform grid is allowed.
  nbz = 0;
  z0 = INPUTS.zmin;  z1=INPUTS.zmax-DZBIN_CCPRIOR+1.0E-5;
  sprintf( ZBIN->varName,"z(CCprior)" );

  for ( z=z0; z < z1; z += DZBIN_CCPRIOR ) {
    zlo = z;
    zhi = z + DZBIN_CCPRIOR ;
    ZBIN->lo[nbz] = zlo ;
    ZBIN->hi[nbz] = zhi ;
    ZBIN->avg[nbz] = 0.5*(zlo+zhi);
    ZBIN->n_perbin[nbz] = 0 ;
    n_perbin_Ibc[nbz] = 0 ; // just for local table
    n_perbin_II[nbz]  = 0 ;
    n_perbin_other[nbz]  = 0 ;
    nbz++ ;  ZBIN->nbin = nbz;
  }
  

  int icc, iz, IZ, IS_Ibc, IS_II;
  for(icc=0; icc < INFO_CC->NROW; icc++ ) {

    z = INFO_CC->z[icc];
    SIM_TYPE = INFO_CC->SIM_TYPE_INDEX[icc];

    sprintf(MSGERR,"%s: z=%f", fnam, z);
    IZ = IBINFUN(z, ZBIN, 1, MSGERR );

    // printf(" xxxx icc=%3d  z=%f IZ=%d \n", icc, z, IZ); fflush(stdout);

    INFO_CC->iz_index[icc] = IZ ;
    ZBIN->n_perbin[IZ]++ ;

    IS_Ibc =  ( SIM_TYPE >= INPUTS.simtype_Ibc[0] && 
		SIM_TYPE <= INPUTS.simtype_Ibc[1] ) ;
    IS_II  =  ( SIM_TYPE >= INPUTS.simtype_II[0] && 
		SIM_TYPE <= INPUTS.simtype_II[1] ) ;

    if ( IS_Ibc ) { 
      n_perbin_Ibc[IZ]++ ; 
      INFO_CC->icc_class[icc] = ICC_FITCLASS_IBC ;
    }
    else if( IS_II ) { 
      n_perbin_II[IZ]++ ; 
      INFO_CC->icc_class[icc] = ICC_FITCLASS_II ;
    }
    else {
      n_perbin_other[IZ]++ ; 
      INFO_CC->icc_class[icc] = ICC_FITCLASS_OTHER ;
    }

  }


  // print summary info for each z-bin
  printf("\n\t   iz  zmin   zmax     NCC = N(Ibc) + N(II) + Other\n" );
  for(iz=0; iz < ZBIN->nbin; iz++ ) {
    printf("\t  %3d  %.3f  %.3f  %5d = %4d  + %4d  + %4d \n",
	   iz, ZBIN->lo[iz], ZBIN->hi[iz], 
	   ZBIN->n_perbin[iz] ,
	   n_perbin_Ibc[iz],
	   n_perbin_II[iz],
	   n_perbin_other[iz]   );
  }

  fflush(stdout);

  return ;

} // end setup_zbins_CCprior


void setup_MUZMAP_CCprior(int IDSAMPLE, 
			  SIMFILE_INFO_DEF *INFO_CC, MUZMAP_DEF *MUZMAP) {

  // setup dmu vs. z map binning in each z bins.
  // Do test with input alpha,beta,M0.
  // Only MUZMAP is used; INFO_CC is just used to test the PDF function.
  //
  // Jun 26 2017: fix index bug for cospar.

  // Feb 2016: hard-code DMU bins, but maybe later allow 
  //     user-defined binning
  //
  // Jun 19 2018:
  //  + if no alpha/beta binning (e..g, 1D biasCor) then set
  //     MUZMAP->alpha,beta = user-init values.
  //
#define DMUMIN_CCPRIOR  -1.5  // add buffer bin to allow for interp
#define DMUMAX_CCPRIOR  +4.0  // add buffer bin to allow for interp
#define DMUBIN_CCPRIOR   0.5  

  int NBINa, NBINb, ia, ib;
  double  SUMa=0.0, SUMb=0.0 ;
    
  int    nbin_dmu ;
  double dmu, dmulo, dmuhi ;
  char fnam[] = "setup_MUZMAP_CCprior" ;

  // ------------ BEGIN -------------

    
  nbin_dmu = 0;
  for ( dmu  = DMUMIN_CCPRIOR; dmu < DMUMAX_CCPRIOR-DMUBIN_CCPRIOR+.0001; 
	dmu += DMUBIN_CCPRIOR ) {
    dmulo = dmu;
    dmuhi = dmu + DMUBIN_CCPRIOR ;
    MUZMAP->DMUBIN.lo[nbin_dmu] = dmulo ;
    MUZMAP->DMUBIN.hi[nbin_dmu] = dmuhi ;
    MUZMAP->DMUBIN.avg[nbin_dmu] = 0.5 * ( dmulo + dmuhi );
    nbin_dmu++ ;  MUZMAP->DMUBIN.nbin = nbin_dmu ;
  }
  

  // print DMU distribution for first redshift bin
  MUZMAP->alpha = INPUTS.parval[IPAR_ALPHA0] ;
  MUZMAP->beta  = INPUTS.parval[IPAR_BETA0] ;
  MUZMAP->M0    = INPUTS.nommag0 ;


  MUZMAP->cosPar[0] = INPUTS.parval[IPAR_OL] ;   // OL
  MUZMAP->cosPar[1] = INPUTS.parval[IPAR_Ok] ;   // Ok
  MUZMAP->cosPar[2] = INPUTS.parval[IPAR_w0] ;  // w0
  MUZMAP->cosPar[3] = INPUTS.parval[IPAR_wa] ;  // wa
    
  // if there is a biasCor map, set alpha,beta to the average
  // since that should be a better estimate in case user input
  // starts alpha,beta way off.
  if ( simdata_bias.NROW > 0 ) {
    NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
    NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;

    for(ia=0; ia<NBINa; ia++ )
      { SUMa += simdata_bias.BININFO_SIM_ALPHA.avg[ia]; }
    for(ib=0; ib<NBINb; ib++ )
      { SUMb += simdata_bias.BININFO_SIM_BETA.avg[ib]; }

    if ( NBINa>0  && NBINb>0 ) {
      MUZMAP->alpha = SUMa / (double)NBINa ;
      MUZMAP->beta  = SUMb / (double)NBINb ;
    }
    else {
      // user input is only choice (Jun 19 2018)
      MUZMAP->alpha = INPUTS.parval[IPAR_ALPHA0] ;
      MUZMAP->beta  = INPUTS.parval[IPAR_BETA0] ;
    }
  }

    
  
  // get DMUPDF in each z bin
  setup_DMUPDF_CCprior(IDSAMPLE, INFO_CC, MUZMAP );  
  
  // -------------------------------
  // print sumary of MU-bins
  int IZTEST = 2;
  dump_DMUPDF_CCprior(IDSAMPLE, IZTEST, MUZMAP) ;
    
  return ;

}  // end setup_MUZMAP_CCprior


void  dump_DMUPDF_CCprior(int IDSAMPLE, int IZ, MUZMAP_DEF *MUZMAP) {

  int imu, NCC ;
  double z    = MUZMAP->ZBIN.avg[IZ] ;
  //  char fnam[] = "dump_DMUPDF_CCprior" ;

  // -------------- BEGIN ---------------

  printf("\n  %s : ", SAMPLE_BIASCOR[IDSAMPLE].NAME  );
  printf("  a=%.3f b=%.3f \n",
	 MUZMAP->alpha, MUZMAP->beta );

  printf("   <DMU>=%.3f  RMS(DMU)=%.3f \n"
	 , MUZMAP->DMUAVG[IDSAMPLE][ICC_FITCLASS_TOT][IZ] 
	 , MUZMAP->DMURMS[IDSAMPLE][ICC_FITCLASS_TOT][IZ]
	 );
  /*
  printf("  OL,Ok,w0,wa = %.3f, %.3f, %.3f, %.3f \n",	 
	 MUZMAP->cosPar[0], MUZMAP->cosPar[1],
	 MUZMAP->cosPar[2], MUZMAP->cosPar[3]  );  */
  
  printf("\t                          PDF(z=%.2f) for \n", z );
  printf("\t  imu  dMU(min) dMU(max)     Ibc    II    OTHER    Total\n" );

  for(imu=0; imu < MUZMAP->DMUBIN.nbin; imu++ ) {
    NCC = MUZMAP->NDMUPDF[IDSAMPLE][ICC_FITCLASS_TOT][IZ][imu] ;    
    printf("\t  %3d  %6.3f  %6.3f       %.3f  %.3f  %.3f   %.3f (%d) \n"
	   ,imu, MUZMAP->DMUBIN.lo[imu], MUZMAP->DMUBIN.hi[imu]
	   ,MUZMAP->DMUPDF[IDSAMPLE][ICC_FITCLASS_IBC][IZ][imu]
	   ,MUZMAP->DMUPDF[IDSAMPLE][ICC_FITCLASS_II][IZ][imu] 
	   ,MUZMAP->DMUPDF[IDSAMPLE][ICC_FITCLASS_OTHER][IZ][imu] 
	   ,MUZMAP->DMUPDF[IDSAMPLE][ICC_FITCLASS_TOT][IZ][imu]
	   ,NCC
	   ) ;
  }  
  fflush(stdout);

  return ;
  
} // end dump_DMUPDF_CCprior


void  setup_DMUPDF_CCprior(int IDSAMPLE, 
			   SIMFILE_INFO_DEF *INFO_CC, MUZMAP_DEF *MUZMAP ) {

  // Created Feb 2016
  // for input redshift (z) and DMU map (*MUZMAP),
  // return the normalized PDF of the DMU distribution (*DMUPDF)
  // Compute PDF at center of two neighboring z bins 
  // and interpolate.
  //
  // This function is called for each fcn call, so no print statements !
  
  int imu, NMUBIN, NZBIN, ia, ib, iz, icc, iclass, idsample, NCCUSE ;
  int NCC[MXCC_FITCLASS][MAXBIN_z][MAXMUBIN]; 
  int NCC_SUM[MXCC_FITCLASS][MAXBIN_z];  // integrals over dmu
  int idata, IDATA;
  double SUM_DMU[MXCC_FITCLASS][MAXBIN_z];
  double SUMSQ_DMU[MXCC_FITCLASS][MAXBIN_z];
  double z, c, x1, mB, dl, mu, mumodel, dmu;
  double a, b, M0, zdif, zdifmin ;
  double XMU, XCC, XCC_SUM, DMUBIN ;
  char *name ;

  // muBias declarations
  int USE_BIASCOR  = ( simdata_bias.NROW > 0 ) ;
  BIASCORLIST_DEF     BIASCORLIST ;
  INTERPWGT_ALPHABETA INTERPWGT ;
  FITPARBIAS_DEF FITPARBIAS_TMP[MXa][MXb] ; 
  double         MUCOVSCALE_TMP[MXa][MXb] ; 
  double fitParBias[NLCPAR], muBias, muBiasErr, muCOVscale ;

  char fnam[] = "setup_DMUPDF_CCprior" ;
    
  // ----------------- BEGIN ---------------
  
  if ( simdata_ccprior.USEH11 ) { return ; }
    
  // number of DMU bins to make PDF
  NMUBIN = MUZMAP->DMUBIN.nbin ;
  NZBIN  = MUZMAP->ZBIN.nbin ;
  XMU    = (double)NMUBIN ;
  DMUBIN = MUZMAP->DMUBIN.hi[0] - MUZMAP->DMUBIN.lo[0]; 
  
  // init default PDF to flat in case there are not enough CC sim stats
  for(iclass=0; iclass < MXCC_FITCLASS; iclass++ ) {
    for(iz=0; iz < NZBIN; iz++ ) {
      NCC_SUM[iclass][iz]   = 0 ;
      SUM_DMU[iclass][iz]   = 0.0 ;
      SUMSQ_DMU[iclass][iz] = 0.0 ;
      MUZMAP->DMUAVG[IDSAMPLE][iclass][iz] = 0.0 ;
      MUZMAP->DMURMS[IDSAMPLE][iclass][iz] = 9999.0 ; // --> flat distribution
	
      for(imu=0; imu < NMUBIN; imu++ ) { 
	NCC[iclass][iz][imu] = 0 ;
	MUZMAP->NDMUPDF[IDSAMPLE][iclass][iz][imu] = 0 ;
	MUZMAP->DMUPDF[IDSAMPLE][iclass][iz][imu]  = 1.0/XMU ;
	MUZMAP->DMUPDF[IDSAMPLE][iclass][iz][imu] /= DMUBIN ; // Jul 2016
      }
    }
  }


    
  for(ia=0; ia < MXa; ia++ ) {
    for(ib=0; ib < MXb; ib++ ) { MUCOVSCALE_TMP[ia][ib] = 1.0 ; }
  }

  //  printf(" xxx %s: NROW = %d \n", fnam, INFO_CC->NROW);

  a  = MUZMAP->alpha ;
  b  = MUZMAP->beta ;
  M0 = MUZMAP->M0 ;

  if ( simdata_bias.DOCOR_5D ) 
    { fcn_AlphaBetaWGT(a,b,0, &INTERPWGT,fnam ); } // returns INTERPWGT struct

    
  NCCUSE = 0;
  for(icc=0; icc < INFO_CC->NROW; icc++ ) {  

    // require correct IDSAMPLE to continue
    idsample = INFO_CC->idsample[icc] ;
    if( idsample != IDSAMPLE ) { continue ; }

    // keep only CC events in current IZ bin
    iz       = INFO_CC->iz_index[icc] ;   // redshift bin
    iclass   = INFO_CC->icc_class[icc] ;  // Ibc, II
    z    = INFO_CC->z[icc] ;
    c    = INFO_CC->c[icc] ;
    x1   = INFO_CC->x1[icc] ;
    mB   = INFO_CC->mB[icc] ;
    name = INFO_CC->name[icc] ;
    
    mu = mB + a*x1 - b*c - M0 ;
    
    // apply mu-bias Correction if simfile_bias is given
    muBias = muBiasErr = 0.0 ;  muCOVscale = 1.0 ;
    if ( USE_BIASCOR ) {    
      BIASCORLIST.z                = z;
      BIASCORLIST.FITPAR[INDEX_mB] = mB ;
      BIASCORLIST.FITPAR[INDEX_x1] = x1 ;
      BIASCORLIST.FITPAR[INDEX_c ] = c  ;
      BIASCORLIST.alpha            = a ;
      BIASCORLIST.beta             = b ;
      BIASCORLIST.idsample         = idsample ;

      // transfer FITPARBIAS_ALPHABETA to FITPARBIAS_TMP that has
      // the right 2D structure definition for get_muBias
      load_FITPARBIAS_CCprior(icc,FITPARBIAS_TMP);

      if ( simdata_bias.DOCOR_5D) {
	get_muBias(name, &BIASCORLIST,     // (I) misc inputs
		   FITPARBIAS_TMP,         // (I) bias vs. ia,ib
		   MUCOVSCALE_TMP,         // (I) muCOVscale vs. ia,ib
		   &INTERPWGT,             // (I) wgt at each a,b grid point
		   fitParBias,     // (O) interp bias on mB,x1,c
		   &muBias,        // (O) interp bias on mu
		   &muBiasErr,     // (O) stat-error on above
		   &muCOVscale );  // (O) scale bias on muCOV (not used below)
      }
      else if ( simdata_bias.DOCOR_1D ) {
	debugexit("CCPRIOR does NOT WORK WITH BBC-1D");  // Jun 19 2018
	// muBias   =  ?? 	
      }
	
    } // end biasCor if-block
    
    // compute mucos for each SIM CC event
    dl       = cosmodl_forFit(z, MUZMAP->cosPar) ;
    mumodel  = 5.0*log10(dl) + 25.0 ;
    mumodel += muBias ;     // add bias
    dmu      = mu - mumodel ;
    imu      = IBINFUN(dmu, &MUZMAP->DMUBIN, 0, "" );

    // xxxxxxxxxxxxxxxxxxxxxxx
    if ( icc == -5 ) {
      printf(" xxx -----------  icc=%d  ---------- \n", icc);
      printf(" xxx z=%.3f  c=%.3f  x1=%.3f  mB=%.3f  \n",
	      z, c, x1, mB );
      printf(" xxx alpha=%.3f  beta=%.3f  M0=%.3f \n", 
	     MUZMAP->alpha,  MUZMAP->beta,  MUZMAP->M0 );
      printf(" xxx iz=%d  iclass=%d \n", iz, iclass);
      printf(" xxx dmu = %.3f - %.3f = %.3f  (imu=%d)\n",
	     mu, mumodel, dmu, imu );
    }
    // xxxxxxxxxxxxxxxxxxxxxxx

    if ( imu >= 0 ) {
      NCC_SUM[iclass][iz]++ ;
      NCC[iclass][iz][imu]++ ;

      // always increment the class with CC=total
      NCC_SUM[ICC_FITCLASS_TOT][iz]++ ;
      NCC[ICC_FITCLASS_TOT][iz][imu]++ ;
      NCCUSE++ ;

      // increment sums for AVG & RMS ( July 2016)
      SUM_DMU[iclass][iz]   += dmu ;
      SUMSQ_DMU[iclass][iz] += (dmu*dmu) ;
      SUM_DMU[ICC_FITCLASS_TOT][iz]   += dmu ;
      SUMSQ_DMU[ICC_FITCLASS_TOT][iz] += (dmu*dmu) ;
    }
  
  } // end icc loop


  
  // Normalize each DMU distribution to get PDF in each class,iz bin.
  // If the stats are too low then leave flat/default PDF.
  double SQRMS, S1TMP, S2TMP;   int NTMP ;
  
  for(iclass=0; iclass < MXCC_FITCLASS; iclass++ ) {
    for(iz=0; iz < NZBIN; iz++ ) {

      if ( NCC_SUM[iclass][iz] > 5 ) {
	NTMP    = NCC_SUM[iclass][iz] ;
	XCC_SUM = (double)NTMP;

	// mean and rms
	S1TMP = SUM_DMU[iclass][iz];
	S2TMP = SUMSQ_DMU[iclass][iz];

	MUZMAP->DMUAVG[IDSAMPLE][iclass][iz] = S1TMP/XCC_SUM ;
	MUZMAP->DMURMS[IDSAMPLE][iclass][iz] = RMSfromSUMS(NTMP,S1TMP,S2TMP);
	
	//  SQRMS  = sumsq_nbr/XNLIST - pow( (sum_nbr/XNLIST),2.0) ;

      // prob in each bin
	for(imu=0; imu < NMUBIN; imu++ ) { 
	  XCC = (double)NCC[iclass][iz][imu];
	  MUZMAP->DMUPDF[IDSAMPLE][iclass][iz][imu] = XCC / XCC_SUM ;
	  MUZMAP->NDMUPDF[IDSAMPLE][iclass][iz][imu]= NCC[iclass][iz][imu];
	  
	  // divide by DMU bin size to get normalized integral
	  // assuming that each DMU bin is linearly interpolated.
	  MUZMAP->DMUPDF[IDSAMPLE][iclass][iz][imu] /= DMUBIN ; 
	}
      }

    }   // and iz loop
  }  // end iclass loop

  
  return ;
  
} // end setup_DMUPDF_CCprior


// ===================================================
void load_FITPARBIAS_CCprior(int icc, FITPARBIAS_DEF
			     (*FITPARBIAS_ABGRID)[MAXBIN_BIASCOR_BETA] )
{

  // Created Jun 7, 2016
  // load FITPARBIAS struct from global function of icc= CC index.
  // The global structure does not have the right form to pass
  // to get_muBias, so load FITPARBIAS that can be passed to get_muBias.

  int  NBINa   = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = simdata_bias.BININFO_SIM_BETA.nbin ;
  int  ia, ib, ipar ;
  //  char fnam[] = "load_FITPARBIAS_CCprior" ;

  // ----------- BEGIN ------------

  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(ipar=0; ipar<NLCPAR; ipar++ ) {

	FITPARBIAS_ABGRID[ia][ib].VAL[ipar] =
	  simdata_ccprior.SIMFILE_INFO_CC.FITPARBIAS_ALPHABETA[icc][ia][ib].VAL[ipar] ;

	FITPARBIAS_ABGRID[ia][ib].ERR[ipar] =
	  simdata_ccprior.SIMFILE_INFO_CC.FITPARBIAS_ALPHABETA[icc][ia][ib].ERR[ipar] ;

	FITPARBIAS_ABGRID[ia][ib].RMS[ipar] =
	  simdata_ccprior.SIMFILE_INFO_CC.FITPARBIAS_ALPHABETA[icc][ia][ib].RMS[ipar] ;
      }
    }
  }
      
} // end load_FITPARBIAS_CCprior


// ===================================================
double prob_CCprior_H11(int n, double dmu, 
			double *H11_fitpar, 
			double *sqsigCC_out, double *sigCC_chi2penalty ) {

  // Created Jun 19 2018
  // Compute CC prob using Hlozek 2011 formalism.
  // Inputs
  //   n   = index of data event.
  //   dmu = mu-residual
  //   H11_fitpar = H11 fit params
  //
  // Ouput
  //   sqsigCC = square of muerr
  //   sigCC_chi2penalty = penalty chi2 to prevent sigcc<0
  //
  // July 9 2018: 
  //  * sigCC -> 2nd order poly(z) instead of constant
  //  * fix awful bug setting alpha, beta; use INPUTS.PARVAL, not INPUTS.ipar
  //
  // July 17 2018: add sigCC_chi2penalty  output

  double *zPoly_mu  = &H11_fitpar[0];
  double *zPoly_sig = &H11_fitpar[3];
  char  *name  = data[n].name ;
  double z     = data[n].zhd ;
  double zerr  = data[n].zhderr ;
  double zsq   = z*z ;
  double prob  = 0.0 ;
  double arg, sqsigCC, sigCC, sigintCC, CCpoly, dmuCC ;
  double alpha, beta, muerrsq_raw ;
  char fnam[] = "prob_CCprior_H11" ;

  // ----------- BEGIN ------------

  CCpoly   = zPoly_mu[0]  + (zPoly_mu[1]*z)  + (zPoly_mu[2]*zsq) ;
  sigintCC = zPoly_sig[0] + (zPoly_sig[1]*z) + (zPoly_sig[2]*zsq) ;

  // compute muerrsq without intrinsic scatter
  alpha = INPUTS.parval[IPAR_ALPHA0];  // fix a,b as in Jone 2016
  beta  = INPUTS.parval[IPAR_BETA0] ;

  muerrsq_raw = fcn_muerrsq(name, alpha, beta, data[n].covmat_fit,
			    z,zerr, 0);

  sqsigCC = muerrsq_raw + (sigintCC*sigintCC) ;
  dmuCC   = (dmu-CCpoly) ;  // i.e., mucalc -> mucalc+CCpoly
  arg     = 0.5*(dmuCC*dmuCC)/sqsigCC ;
  
  //  sigCC   = data[n].sigCC_last ;
  sigCC   = sqrt(sqsigCC);
  prob    = (PIFAC/sigCC) * exp(-arg) ;

  *sqsigCC_out = sqsigCC ;

  // compute penalty chi2
  *sigCC_chi2penalty = 0.0 ;
  /* xxxx
  arg = sigCC/0.01 ;
  *sigCC_chi2penalty = 1.0 - (1.0/( exp(-arg) + 1.0 )) ;
  *sigCC_chi2penalty *= 1000.0 ;
  xxxx */

  return(prob) ;

} // end prob_CCprior_H11

// ===================================================
double prob_CCprior_sim(int IDSAMPLE, MUZMAP_DEF *MUZMAP, 
		    double z, double dmu, int DUMPFLAG) {

  // Called by fcn function in minuit to return CC prob 
  // for input MUZMAP (includes cosPar) and Hubble resid dmu.
  // Jun 18 2018: pass DUMPFLAG argument for debugging
  
#define FUNDMU_INTERP 1  // sometimes gives crazy errors on alpha,beta
#define FUNDMU_GAUSS  2  // no crazy errors.
  
  int FUNDMU = FUNDMU_GAUSS ;
  int OPT_INTERP = 1; // 1=linear
  int IZ, NBMU ;
  double dmumin, dmumax, dmu_local, prob;
  char fnam[] = "prob_CCprior_sim" ;

  // ----------------- BEGIN ------------
  prob = 0.0 ;
  
  IZ   = IBINFUN(z, &MUZMAP->ZBIN, 0, "");

  if ( FUNDMU == FUNDMU_GAUSS ) {
    // make Guassian approximation with <DMU> and RMS(DMU).
    double AVG, RMS, resid, arg ;
    AVG    = MUZMAP->DMUAVG[IDSAMPLE][ICC_FITCLASS_TOT][IZ];
    RMS    = MUZMAP->DMURMS[IDSAMPLE][ICC_FITCLASS_TOT][IZ];
    resid  = (dmu-AVG)/RMS ; arg=0.5*resid*resid ;
    prob   = (PIFAC/RMS) * exp(-arg);
  }
  
  // ------------------------------------------
  // if we get here, interpolate in dmu bins

  dmu_local = -9.0 ;
  if ( FUNDMU == FUNDMU_INTERP ) {
    NBMU   = MUZMAP->DMUBIN.nbin ;  
    dmumin = MUZMAP->DMUBIN.avg[0];
    dmumax = MUZMAP->DMUBIN.avg[NBMU-1];
    
    // make sure dmu is contained withing grid
    if ( dmu < dmumin ) 
      { dmu_local = dmumin + 0.0001 ; }
    else if ( dmu > dmumax ) 
      { dmu_local = dmumax - 0.0001 ; }
    else
      { dmu_local = dmu ; }
    
    prob = interp_1DFUN(OPT_INTERP, dmu_local, NBMU, 
			&MUZMAP->DMUBIN.avg[0],
			&MUZMAP->DMUPDF[IDSAMPLE][ICC_FITCLASS_TOT][IZ][0],
			fnam) ;
  }



  // xxxxxxxxxxxxxxxxxxx
  if ( IZ == -2 && dmu_local > 1.0 ) {
    printf(" xxx IDSAMPLE=%d  dmu_local=%6.3f  prob=%.3f \n",
	   IDSAMPLE, dmu_local, prob); fflush(stdout);
  }
  // xxxxxxxxxxxxxxxxxxx

  
  // -------------------------------
  // if we get here, return zero prob.
  return(prob);

} // end prob_CCprior_sim


// ===========================================
int IBINFUN(double VAL, BININFO_DEF *BIN, int OPT, char *MSG_ABORT ) {

  // Generic function to  return integer bin for VAL
  // and bins defined by struct *BIN.
  // If VAL is outside range, return -9
  // 
  // If OPT == 1  abort if VAL is outside range
  // If OPT == 2  if VAL is outside, return edge bin (do not abort)
  // If OPT == 7  dump
  
  int ibin, IBIN, nbin;
  int LDMP = ( OPT == 7 ) ;
  double lo, hi;  
  char fnam[] = "IBINFUN" ;
  // --------- BEGIN -----------

  IBIN = -9; ibin=0;
  nbin = BIN->nbin ;
  
  if ( LDMP ) {
    printf(" 1. xxx %s DUMP: nbin=%d  VAL(%s)=%f \n",
	   fnam, nbin, BIN->varName, VAL);
    fflush(stdout);
  }
  

  while ( IBIN < 0 && ibin < nbin) {
    lo = BIN->lo[ibin] ;
    hi = BIN->hi[ibin] ;  if ( ibin==nbin-1 ) { hi += 1.0E-10; }

    if ( LDMP ) {
      printf(" 2. xxx %s DUMP: ibin=%d  lo=%.3f  hi=%.3f \n",
	     fnam, ibin, lo, hi);   fflush(stdout);
    }
    
    if ( VAL >= lo  &&  VAL <= hi )
      { IBIN = ibin ;  return(IBIN); }
    
    ibin++ ;
  }

  if ( LDMP ) {
    printf(" 3. xxx %s DUMP: IBIN,ibin=%d,%d \n",
	   fnam, IBIN, ibin );   fflush(stdout);
  }
    
  if ( OPT == 2 ) {
    if ( VAL < BIN->avg[0]      ) { IBIN=0; }
    if ( VAL > BIN->avg[nbin-1] ) { IBIN=nbin-1; }
  }

  if ( OPT == 1 ) {
    sprintf(c1err,"Cannot find IBIN for %s=%f "
	    "(range: %.3f to %.3f)", 
	    BIN->varName, VAL, BIN->lo[0], BIN->hi[BIN->nbin-1] );
    errmsg(SEV_FATAL, 0, fnam, c1err, MSG_ABORT);  
  }
  
  // xxxx  if ( VAL > BIN-hi[nbin-1] ) { IBIN = 999999 ; }

  return(IBIN);

} // end IBINFUN

void copy_BININFO(BININFO_DEF *BIN1, BININFO_DEF *BIN2) {

  int i, nbin ;
  //  char fnam[] = "copy_BININFO" ;

  // copy contents of BIN1 into BIN2
  BIN2->nbin     = BIN1->nbin ;
  BIN2->binSize  = BIN1->binSize ;
  sprintf(BIN2->varName, "%s", BIN1->varName);

  nbin = BIN1->nbin ;
  
  for(i=0; i < nbin; i++ ) {
    BIN2->lo[i]       = BIN1->lo[i];
    BIN2->hi[i]       = BIN1->hi[i];
    BIN2->avg[i]      = BIN1->avg[i];
    BIN2->n_perbin[i] = 0 ;
  }

} // end copy_BININFO


// ==================================================
void set_ERRMASK(int isn_raw, int isn_data) {

  // Created May 30, 2012 by R.Kessler
  // 
  // Check cuts and errors;
  // load data[isn_data].errmask
  // 
  // Inputs:
  //  isn_raw : index to rawdata struct
  //  isn_data : index to data struct.
  //
  // Sep 05 2016: remove reference to hgtype and check host_logmass.
  //
  // May 11 2017: fix bug accessing host_logmass.

  int  LCUT, icut, ITYPE, FOUND, LDMP ;
  double cutvar_local[MXCUTWIN], z, x1, c ;
  char   *SNname ;
  //  char fnam[]=  "set_ERRMASK";

  // ---------- BEGIN ---------


  // apply user CUTWIN options (4/24/2012)
  LCUT = 1;
  LDMP = 0;

  SNname = data[isn_data].name ;

  for ( icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
    { cutvar_local[icut] = rawdata.CUTVAL[icut][isn_raw] ; }
  LCUT = apply_CUTWIN(1, rawdata.DOFLAG_CUTWIN, cutvar_local);


  if ( LCUT == 0 ) { 
    data[isn_data].errmask += ERRMASK_CUTWIN ; 
    NSET_ERRMASK[ERRMASK_CUTWIN]++ ;
  }
  

  // -----------------------------------------
  // apply legacy cuts (RK - March 31, 2010)
  // Flag events -- not used in fit but processed to output file      
  // These cuts are redundant with above LCUT to allow older inputs.
  z = rawdata.zhd[isn_raw] ;

  if ( (z<INPUTS.zmin) || (z>INPUTS.zmax))   { 
    data[isn_data].errmask += ERRMASK_z; 
    NSET_ERRMASK[ERRMASK_z]++ ;
  }

  x1 = rawdata.x1[isn_raw] ;
  if ( (x1<INPUTS.x1min) || (x1>INPUTS.x1max))  { 
    data[isn_data].errmask += ERRMASK_x1; 
    NSET_ERRMASK[ERRMASK_x1]++ ;
  }

  c = rawdata.x3[isn_raw] ;
  if ( (c<INPUTS.cmin) || (c>INPUTS.cmax) ) { 
    data[isn_data].errmask += ERRMASK_c; 
    NSET_ERRMASK[ERRMASK_c]++ ;
  }
  
  
  // check simcc pre-scale (Dec 2015)
  if ( reject_simData(isn_raw) ) {
    data[isn_data].errmask += ERRMASK_SIMPS; 
    NSET_ERRMASK[ERRMASK_SIMPS]++ ;
  }

  //June 20 -- Move code & change to keep SN but flag not to use in fit
  if ( rawdata.x0err[isn_raw] <= INPUTS.maxerr_abort_x0 || 
       rawdata.x1err[isn_raw] <= INPUTS.maxerr_abort_x1 ||  
       rawdata.x3err[isn_raw] <= INPUTS.maxerr_abort_c     ) {
    printf(" ERROR-WARNING: SN %s (n=%d): one or more x0,x1,c errors <=0 \n",
	   data[isn_data].name, isn_raw );
    data[isn_data].errmask += ERRMASK_BADERR ;
    NSET_ERRMASK[ERRMASK_BADERR]++ ;
    goto SKIP_COVTEST ;
  }

  // May 2016  flag crazy COV and mBerr 
  int BADCOV   =  0;
  double sf    = -2.5/(rawdata.x0[isn_raw]*LOGTEN) ;
  double mBerr = fabs(sf) * rawdata.x0err[isn_raw] ;
  if ( fabs(rawdata.x0x1_c[isn_raw]) > 8.0 ) { BADCOV=1; } // COV_x0x1
  if ( fabs(rawdata.x0x3_c[isn_raw]) > 8.0 ) { BADCOV=1; } // COV_x0c
  if ( fabs(rawdata.x1x3_c[isn_raw]) > 8.0 ) { BADCOV=1; } // COV_x1c
  if ( mBerr > 1.0 ) { BADCOV=1; }
  
  if ( BADCOV ) {
    printf(" ERROR-WARNING: SN %s (n=%d): one or more crazy COV. \n",
	   data[isn_data].name, isn_raw );
    data[isn_data].errmask += ERRMASK_BADERR ;
    NSET_ERRMASK[ERRMASK_BADERR]++ ;
  }

 SKIP_COVTEST:

  ITYPE                 = (int)rawdata.sntype[isn_raw] ;
  data[isn_data].sntype = ITYPE ;
  if ( FOUNDKEY_SNTYPE && INPUTS.Nsntype > 0 ) {
    FOUND = 0 ;
    for ( icut=0; icut < INPUTS.Nsntype; icut++ ) {
      if ( ITYPE == INPUTS.sntype[icut] ) { FOUND = 1; }
    }

    if ( FOUND == 0 ) { 
      data[isn_data].errmask += ERRMASK_sntype ; 
      NSET_ERRMASK[ERRMASK_sntype]++ ;
    }
  }
  
  int errmask ;
  errmask = data[isn_data].errmask ;
  if ( errmask != 0 ) { FITINP.NSNERR++ ; }

  return ;

} // end of set_ERRMASK


// ===================================================
void  getComment_ERRMASK(int mask, char *comment) {

  // for input error 'mask', return comment.

  // --------------- BEGIN ---------------

  comment[0] = 0 ;

  if ( (mask & ERRMASK_x1) > 0 ) 
    { strcat(comment,"x1 "); }

  if ( (mask & ERRMASK_c) > 0 ) 
    { strcat(comment,"c "); }

  if ( (mask & ERRMASK_MINBIN) > 0 ) 
    { strcat(comment,"zbin "); }

  if ( (mask & ERRMASK_z) > 0 ) 
    { strcat(comment,"z "); }

  if ( (mask & ERRMASK_sntype) > 0 ) 
    { strcat(comment,"sntype "); }

  if ( (mask & ERRMASK_HOST) > 0 ) 
    { strcat(comment,"HOST "); }

  if ( (mask & ERRMASK_CUTWIN) > 0 ) 
    { strcat(comment,"CUTWIN "); }

  if ( (mask & ERRMASK_BADERR) > 0 ) 
    { strcat(comment,"BADERR "); }

  if ( (mask & ERRMASK_SPLITRAN) > 0 ) 
    { strcat(comment,"not-in-random-sample "); }

  if ( (mask & ERRMASK_SIMPS) > 0 ) 
    { strcat(comment,"sim-PreScale "); }

  return ;

} // end getComment_ERRMASK


// =============================================
int reject_simData(int isn_raw) {

  // Created Apr 14 2017 by R.Kessler
  // For simulated event check prescale_simData and prescale_simCC
  // to see if this event is rejected.
  //   * prescale_simData pre-scales all sim events randomly
  //   * prescale_simCC   pre-scales only the true CC subset
  // Note that data are always accepted.
  // Function returns 1 to reject, 0 to keep.

  int REJECT = 0 ;
  float XN, XNPS ;
  char fnam[] = "reject_simData" ;

  // ------------- BEGIN -----------------

  // return accept for data
  if ( FOUNDKEY_SIM == 0 ) { return(REJECT); }

  // increment NSIMDATA counter
  NSIMDATA++ ;
  XN    = (float)NSIMDATA ;
  XNPS  = (float)INPUTS.prescale_simData ;
  if ( fmodf( XN, XNPS ) != 0 )  { REJECT = 1 ; }

  // increment NSIMCC counter
  if ( rawdata.sim_nonIa_index[isn_raw] > 0 ) {
    NSIMCC++ ;
    XN    = (float)NSIMCC ;
    XNPS  = (float)INPUTS.prescale_simCC ;
    if ( fmodf( XN, XNPS ) != 0 ) { REJECT = 1 ; }
    if ( INPUTS.prescale_simCC > 9999.0 ) { REJECT=1 ; }
  }
  return(REJECT);

} // end reject_simData


// ====================================================
int  rawdata_malloc(char *filename) {

  // Created Jul 13, 2011 by R.Kessler
  // Read number of rows in file and
  // allocate memory for rawdata and data structures.
  // Return NROW.
  //
  // Apr 24, 2012: malloc rawdata.CUTVAL
  // Apr 11, 2013: Z -> z and ZERR -> zERR
  // Apr 30, 2013: add x1ERR and cERR to list of allowed cut-vars.
  // Jan 19, 2016: use nrow_read() function to get number of rows
  // Feb 16, 2016: malloc pIa
  // Jul 18, 2016: double -> int for a few items
  // Sep 05, 2016: add host_logmass
  // Jan 08, 2018: add vpec and vpecerr, and zhel[err]
  // Mar 21, 2018: add opt_photoz

  int NROW, MEMTOT, MEMD, MEMI, isn ;
  char fnam[] = "rawdata_malloc" ;

  // ---------- BEGIN ------------

  NROW = SNTABLE_NEVT(filename,TABLENAME_FITRES);
  NROW += 10;

  printf("\t Found %d rows \n", NROW);

  MEMTOT = 0;
  MEMD   = NROW * sizeof(double);
  MEMI   = NROW * sizeof(int);
  

  rawdata.name =  (char**)malloc(NROW * sizeof(char*));
  for(isn=0; isn < NROW; isn++ )  { 
    rawdata.name[isn] =  (char*)malloc(MXCHAR_CCID * sizeof(char));  
    MEMTOT += MXCHAR_CCID ;
  }

  if ( INPUTS.use_fieldGroup_biasCor ) { 
    rawdata.field =  (char**)malloc(NROW * sizeof(char*));
    for(isn=0; isn < NROW; isn++ )  { 
      rawdata.field[isn] =  (char*)malloc(MXCHAR_CCID * sizeof(char));  
      MEMTOT += MXCHAR_CCID ;
    }    
  }

  MEMTOT += MEMD ; rawdata.zhd       = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.zhel      = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.vpec      = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x0        = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x1        = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x3        = (double*)malloc(MEMD);

  MEMTOT += MEMD ; rawdata.zhderr    = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.zhelerr   = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.vpecerr   = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x0err     = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x1err     = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x3err     = (double*)malloc(MEMD);

  MEMTOT += MEMD ; rawdata.x0x1_c    = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x0x3_c    = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.x1x3_c    = (double*)malloc(MEMD);

  if ( simdata_ccprior.USE )  { 
    MEMTOT += MEMD ; rawdata.pIa   = (double*)malloc(MEMD); 

    /* xxxxxxxx mark delete Feb 5 2019 xxxxxxxxx
    MEMTOT += MEMD ; rawdata.pII   = (double*)malloc(MEMD); 
    MEMTOT += MEMD ; rawdata.pIbc  = (double*)malloc(MEMD); 
    xxxxxxxxxxx */
  }

  MEMTOT += MEMD ; rawdata.simz      = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.sim_mu    = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.simx0     = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.simx1     = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.simc      = (double*)malloc(MEMD);


  MEMTOT += MEMI ; rawdata.idsurvey   = (int*)malloc(MEMI);
  MEMTOT += MEMI ; rawdata.opt_photoz = (int*)malloc(MEMI);
  MEMTOT += MEMI ; rawdata.sntype     = (int*)malloc(MEMI);
  MEMTOT += MEMI ; rawdata.sim_nonIa_index = (int*)malloc(MEMI);

  MEMTOT += MEMD ; rawdata.dsurvey  = (double*)malloc(MEMD);
  MEMTOT += MEMD ; rawdata.dtype    = (double*)malloc(MEMD);

  // tack on the cut-variables. If cut-variable already exists above,
  // then just set pointer to rawdata.xxx above; otherwise allocate
  // more memory and set RDFLAG. The main point here is to use only
  // one array per variable.

  int icut;
  char *cutname ;

  for ( icut=0; icut < INPUTS.NCUTWIN ; icut++ )  { 
    INPUTS.CUTWIN_RDFLAG[icut] = 0;
    cutname = INPUTS.CUTWIN_NAME[icut] ;

    if ( strcmp(cutname,"x0") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.x0 ; }

    else if ( strcmp(cutname,"x1") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.x1 ; }
    else if ( strcmp(cutname,"x1ERR") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.x1err ; }

    else if ( strcmp(cutname,"c") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.x3 ; }
    else if ( strcmp(cutname,"cERR") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.x3err ; }

    else if ( strcmp(cutname,"z") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.zhd ; }
    else if ( strcmp(cutname,"Z") == 0 )    // allow old name (Apr 2013)
      { rawdata.CUTVAL[icut] = rawdata.zhd ; }

    else if ( strcmp(cutname,"zERR") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.zhderr ; }
    else if ( strcmp(cutname,"ZERR") == 0 )     // allow old name
      { rawdata.CUTVAL[icut] = rawdata.zhderr ; }

    else if ( strcmp(cutname,"IDSURVEY") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.dsurvey ; }

    else if ( strcmp(cutname,"TYPE") == 0 ) 
      { rawdata.CUTVAL[icut] = rawdata.dtype ; }

    else {
      INPUTS.CUTWIN_RDFLAG[icut] = 1 ;
      MEMTOT += MEMD ; rawdata.CUTVAL[icut] = (double*)malloc(MEMD);  
    }
  } // icut

  rawdata.ICUTWIN_GAMMA = -9 ;  // index for gamma-variable

  printf("\t Allocated %6.3f MB for rawdata structure. \n",
	 (double)MEMTOT/1.0E6 );

  // malloc data structure containing rawdata after cuts.
  MEMTOT = NROW*sizeof(SNDATA_DEF) ;
  data = (SNDATA_DEF*)malloc( MEMTOT );

  printf("\t Allocated %6.3f MB for data structure. \n",
	 (double)MEMTOT/1.0E6 );

  fflush(stdout);

  return(NROW);

} // endof rawdata_malloc


// ====================================================
void malloc_simdata_biasCor(int opt, int NROW) {

  // Created Jan 2016 by R.Kessler
  //
  // Inputs:
  //   opt> 0 --> allocate & return number of rows in simFile
  //   opt< 0 --> free
  // 
  //   NROW = size of array to malloc

  int MEMTOT, MEMD, MEMF, MEMI, MEMI2, isn, icut, ipar, ipar2 ;
  int IDEAL_COVINT = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  char fnam[] = "malloc_simdata_biasCor" ;

  // ---------- BEGIN ------------

  if( opt > 0 ) {
    
    MEMTOT = 0;
    MEMD   = NROW * sizeof(double);
    MEMF   = NROW * sizeof(float);
    MEMI   = NROW * sizeof(int);
    MEMI2  = NROW * sizeof(short int);
    
    simdata_bias.name =  (char**)malloc(NROW * sizeof(char*));
    for(isn=0; isn < NROW; isn++ )  {  
      simdata_bias.name[isn] =  (char*)malloc(MXCHAR_CCID * sizeof(char));  
      MEMTOT += MXCHAR_CCID; 
    }
    
    if ( INPUTS.use_fieldGroup_biasCor ) {
      simdata_bias.field =  (char**)malloc(NROW * sizeof(char*));
      for(isn=0; isn < NROW; isn++ )  {  
	simdata_bias.field[isn] = (char*)malloc(MXCHAR_CCID * sizeof(char));  
	MEMTOT += MXCHAR_CCID; 
      }
    }

    MEMTOT += MEMF ; simdata_bias.z              = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.zerr           = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.vpec           = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.vpecerr        = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.SNRMAX         = (float*)malloc(MEMF);
    MEMTOT += MEMI ; simdata_bias.SIM_NONIA_INDEX = (int*)malloc(MEMI);
    MEMTOT += MEMI ; simdata_bias.idsurvey       = (int*)malloc(MEMI);
    MEMTOT += MEMI ; simdata_bias.opt_photoz     = (int*)malloc(MEMI);
    MEMTOT += MEMI ; simdata_bias.idsample       = (int*)malloc(MEMI);
		                          
    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      MEMTOT += MEMF ; simdata_bias.FITVAL[ipar]  = (float*)malloc(MEMF);
      MEMTOT += MEMF ; simdata_bias.FITERR[ipar]  = (float*)malloc(MEMF);
      MEMTOT += MEMF ; simdata_bias.SIMVAL[ipar]  = (float*)malloc(MEMF);

      for(ipar2=0; ipar2<NLCPAR; ipar2++ ) {
	MEMTOT += MEMF ; simdata_bias.COV_FITVAL[ipar][ipar2]
			   = (float*)malloc(MEMF);
      }
      
    }
   
    MEMTOT += MEMF ; simdata_bias.COV_x0x1 = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.COV_x0c  = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.COV_x1c  = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.x0       = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.x0err    = (float*)malloc(MEMF);

    MEMTOT += MEMF ; simdata_bias.SIM_ALPHA  = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.SIM_BETA   = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.SIM_DLMAG  = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.SIM_ZCMB   = (float*)malloc(MEMF);
    MEMTOT += MEMF ; simdata_bias.SIM_VPEC   = (float*)malloc(MEMF);

    MEMTOT += MEMI2; simdata_bias.iz  = (short int*)malloc(MEMI2);
    MEMTOT += MEMI2; simdata_bias.IA  = (short int*)malloc(MEMI2);
    MEMTOT += MEMI2; simdata_bias.IB  = (short int*)malloc(MEMI2);
    MEMTOT += MEMI2; simdata_bias.IZ  = (short int*)malloc(MEMI2);

    if ( IDEAL_COVINT ) {
      for(ipar=0; ipar < NLCPAR; ipar++ ) {
	MEMTOT += MEMF ; simdata_bias.FITVAL_IDEAL[ipar]  
			   = (float*)malloc(MEMF);
      }
      MEMTOT += MEMF ; simdata_bias.x0_IDEAL = (float*)malloc(MEMF);      
    }

    // tack on variables defined by CUTWIN.
    // Initialize each CUTVAL in case some are not read
    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {    
      MEMTOT += MEMF ; simdata_bias.CUTVAL[icut] = (float*)malloc(MEMF);
      for(isn=0; isn < NROW; isn++ ) 
	{ simdata_bias.CUTVAL[icut][isn] = -9.0; } 
    }

    printf("\t Allocate %6.3f MB for simdata_bias structure. \n",
	   (double)MEMTOT/1.0E6 );   
    fflush(stdout);
    
  }
  else {
    printf("\t Free simdata_bias memory. \n");
    fflush(stdout);
    
    free(simdata_bias.z) ;
    free(simdata_bias.SIM_NONIA_INDEX) ;
    free(simdata_bias.SNRMAX) ;
    free(simdata_bias.idsurvey) ;
    free(simdata_bias.idsample);
    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      free(simdata_bias.FITVAL[ipar]) ;
      free(simdata_bias.FITERR[ipar]) ;
      free(simdata_bias.SIMVAL[ipar]) ;
      for(ipar2=0; ipar2<NLCPAR; ipar2++ ) 
	{ free(simdata_bias.COV_FITVAL[ipar][ipar2]);	}

    }
  
    free(simdata_bias.COV_x0x1);
    free(simdata_bias.COV_x0c);
    free(simdata_bias.COV_x1c);
    free(simdata_bias.x0);
    free(simdata_bias.x0err);

    free(simdata_bias.SIM_ALPHA);
    free(simdata_bias.SIM_BETA);
    free(simdata_bias.SIM_DLMAG);
    free(simdata_bias.SIM_ZCMB);
    free(simdata_bias.SIM_VPEC);
    free(simdata_bias.IA);
    free(simdata_bias.IB);
    
    if ( IDEAL_COVINT ) {
      free(simdata_bias.x0_IDEAL ) ;
      for(ipar=0; ipar < NLCPAR; ipar++ ) {
	free(simdata_bias.FITVAL_IDEAL[ipar]) ;
      }
    }

    for(isn=0; isn < simdata_bias.LEN_MALLOC; isn++ ) 
      { free(simdata_bias.name[isn]);}
    free(simdata_bias.name);
    
    if ( INPUTS.use_fieldGroup_biasCor ) {
      for(isn=0; isn < simdata_bias.LEN_MALLOC; isn++ ) 
	{ free(simdata_bias.field[isn]);}
      free(simdata_bias.field);
    }

    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
      { free(simdata_bias.CUTVAL[icut]);  }    
  }

  return ;

} // endof malloc_simdata_biasCor



// ====================================================
void malloc_simdata_ccprior(int opt, int NROW, 
			    SIMFILE_INFO_DEF *SIMFILE_INFO) {

  // Created Jan 2016 by R.Kessler
  // opt> 0 --> allocate & return number of rows in simFile
  // opt< 0 --> free

  int  USE_BIASCOR  =   ( simdata_bias.NROW > 0 ) ;
  int  MEMTOT, MEMD, MEMI, MEMB, MEMa, MEMb, isn, icut, ia, ib ;
  int  NBINa     = simdata_bias.BININFO_SIM_ALPHA.nbin ;
  int  NBINb     = simdata_bias.BININFO_SIM_BETA.nbin ;
  //  char fnam[] = "malloc_simdata_ccprior" ;

  // ---------- BEGIN ------------

  
  if( opt > 0 ) {

    SIMFILE_INFO->LEN_MALLOC = NROW ;
    
    MEMTOT = 0;
    MEMD   = NROW * sizeof(double);
    MEMI   = NROW * sizeof(int);
    
    SIMFILE_INFO->name =  (char**)malloc(NROW * sizeof(char*));
    for(isn=0; isn < NROW; isn++ )  {  
      SIMFILE_INFO->name[isn] =  (char*)malloc(MXCHAR_CCID * sizeof(char));  
      MEMTOT += MXCHAR_CCID; 
    }
    
    if( INPUTS.use_fieldGroup_biasCor ) {
      SIMFILE_INFO->field =  (char**)malloc(NROW * sizeof(char*));
      for(isn=0; isn < NROW; isn++ )  {  
	SIMFILE_INFO->field[isn] = (char*)malloc(MXCHAR_CCID * sizeof(char));
	MEMTOT += MXCHAR_CCID; 
      }
    }

    MEMTOT += MEMD ; SIMFILE_INFO->z         = (double*)malloc(MEMD);
    MEMTOT += MEMD ; SIMFILE_INFO->mB        = (double*)malloc(MEMD);
    MEMTOT += MEMD ; SIMFILE_INFO->c         = (double*)malloc(MEMD);
    MEMTOT += MEMD ; SIMFILE_INFO->x1        = (double*)malloc(MEMD);
    MEMTOT += MEMD ; SIMFILE_INFO->x0        = (double*)malloc(MEMD);
    MEMTOT += MEMI ; SIMFILE_INFO->SIM_DLMAG       = (double*)malloc(MEMD);
    MEMTOT += MEMI ; SIMFILE_INFO->SIM_TYPE_INDEX  = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->SIM_NONIA_INDEX = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->iz_index        = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->icc_class       = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->idsurvey        = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->opt_photoz      = (int*)malloc(MEMI);
    MEMTOT += MEMI ; SIMFILE_INFO->idsample        = (int*)malloc(MEMI);


    if ( USE_BIASCOR ) {      
      MEMB   = NROW * sizeof(FITPARBIAS_DEF***);
      SIMFILE_INFO->FITPARBIAS_ALPHABETA =
	    (FITPARBIAS_DEF***)malloc(MEMB);
      MEMTOT += MEMB ;

      for(isn=0; isn < NROW; isn++ ) {
	MEMa = MAXBIN_BIASCOR_ALPHA * sizeof(FITPARBIAS_DEF*) ;
	MEMb = MAXBIN_BIASCOR_BETA  * sizeof(FITPARBIAS_DEF ) ;

	SIMFILE_INFO->FITPARBIAS_ALPHABETA[isn] =
	  (FITPARBIAS_DEF**)malloc(MEMa);
	MEMTOT += MEMa ;
	for(ia=0; ia < MAXBIN_BIASCOR_ALPHA; ia++ ) {
	  SIMFILE_INFO->FITPARBIAS_ALPHABETA[isn][ia] =
	    (FITPARBIAS_DEF*)malloc(MEMb);
	  MEMTOT += MEMb ;
	}  //ia
      } // isn

    } // end USE_BIASCOR
    
    // tack on variables defined by CUTWIN.
    // Initialize each CUTVAL in case some are not read
    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {    
      MEMTOT += MEMD ; SIMFILE_INFO->CUTVAL[icut] = (double*)malloc(MEMD);
      for(isn=0; isn < NROW; isn++ ) 
	{ SIMFILE_INFO->CUTVAL[icut][isn] = -9.0; } 
    }

    printf("\t Allocate %6.3f MB for simdata_ccprior structure. \n",
	   (double)MEMTOT/1.0E6 ); 
    fflush(stdout);
  }
  else {
    printf("\t Free simdata_ccprior memory. \n");
    fflush(stdout);
        
    free(SIMFILE_INFO->z) ;
    free(SIMFILE_INFO->mB) ;
    free(SIMFILE_INFO->c) ;
    free(SIMFILE_INFO->x1) ;
    free(SIMFILE_INFO->SIM_TYPE_INDEX) ;
    free(SIMFILE_INFO->SIM_NONIA_INDEX) ;
    free(SIMFILE_INFO->iz_index) ;
    free(SIMFILE_INFO->icc_class) ;
    free(SIMFILE_INFO->idsurvey ) ;
    free(SIMFILE_INFO->opt_photoz ) ;
    free(SIMFILE_INFO->idsample );

    for(isn=0; isn < SIMFILE_INFO->LEN_MALLOC; isn++ ) 
      { free(SIMFILE_INFO->name[isn]);}
    free(SIMFILE_INFO->name);
    
    if( INPUTS.use_fieldGroup_biasCor ) {
      for(isn=0; isn < SIMFILE_INFO->LEN_MALLOC; isn++ ) 
	{ free(SIMFILE_INFO->field[isn]);}
      free(SIMFILE_INFO->field);
    }

    for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) 
      { free(SIMFILE_INFO->CUTVAL[icut]);  }

    if ( USE_BIASCOR ) {
      for(ia=0; ia <NBINa; ia++ ) {
	for(ib=0; ib<NBINb; ib++ ) {
	  free(SIMFILE_INFO->FITPARBIAS_ALPHABETA[ia][ib]);
	} // ib
      } // ia
      free(SIMFILE_INFO->FITPARBIAS_ALPHABETA);
    }

  }

  return ;

} // endof malloc_simdata_ccprior



// ==============================
void  rawdata_init(void) {
  int n ;

  // Created Jul 13, 2011 by R.Kessler
  // Init rawdata arrays.

  for (n=0;n < MAXSN; ++n) { sprintf(rawdata.name[n],"NULL");  }

  // Apr 2011 (RK) - init values to obvious garbage
  for (n=0;n < MAXSN; ++n) { 
    rawdata.idsurvey[n]  =  -9 ;
    rawdata.sntype[n]    =  -9.0 ;
    rawdata.simx0[n]     =   1.0 ; // => so that mB = 0
    rawdata.simx1[n]     =  -9.0 ; 
    rawdata.simc[n]      =  -9.0 ; 
    rawdata.simz[n]      =  -9.0 ; 
    rawdata.sim_mu[n]    =  -9.0 ; 
    rawdata.sim_nonIa_index[n] = 0 ;
  }

} // end of rawdata_ini


// ********************************************
void parse_parFile(char *parFile ) {

  // April 19, 2010 RK
  // parse file with SALT2mu parameter options
  // (Note: this is NOT the fitres file with SALT2 fit parameters)
  //
  // Apr 12, 2011 RK - read 120 chars intead of 80
  // Apr 16, 2016 RK - read 400 chars instead of 120 to handle long
  //                   file names after a "FITOPT:" key.
  //

  FILE *fdef;
  char *sptr;
  char fnam[] = "parse_parFile" ;

  // ------------------ BEGIN --------------

  // allow for some null options to skip reading file
  if ( strcmp(parFile,"NULL")  == 0 ) return;
  if ( strcmp(parFile,"null")  == 0 ) return;
  if ( strcmp(parFile,"BLANK") == 0 ) return;

  // Look for a file of default options
  fdef = fopen(parFile,"rt");
  if (!fdef) {
    printf("\n FATAL ERROR: could not open SALT2mu parameter-input file:\n");
    printf(" %s \n", parFile);
    printf("\t ***** ABORT ***** \n");
    exit(1); 
  }

  // Just return if file does not exist
  uniqueOverlap("INIT","SALT2mu-input file");

  printf("Reading SALT2mu parameter-input file '%s' : \n" ,  parFile);

  char line[400];
  while (fgets(line,400,fdef)) {

    // skip blank lines and lines starting with comment char (RK 4/23/2012)
    if ( strlen(line) < 3 )       { continue ; }  
    if ( commentchar(line) == 1 ) { continue ; } // see sntools.c

    sptr = strtok(line,"\n");

    // skip keys with colon that are used by a master 
    // SALT2mu_fit.pl script (RK Apr 2013)
    // Aug 24 2016: make exception for key containing group_biasCor
    //    to allow colon-delimeter in input string
    // NOTE: should re-facto to check list of keys allowing colon
    if ( strstr(sptr,"group_biascor") == NULL &&
	 strstr(sptr,":") != NULL ) { continue ; }

    ppar(sptr); // pass entire line
  }

  printf("\n");
  fclose(fdef);

  // --------------------------------------
  // Jun 20 2016: allow ENV for any input file
  ENVreplace(INPUTS.filename,fnam,1);
  ENVreplace(INPUTS.simFile_biasCor,fnam,1);
  ENVreplace(INPUTS.simFile_CCprior,fnam,1);

  return ;

} // end  of parse_parFile


// **************************************
void override_parFile(int argc, char **argv) {

  // April 23, 2012 by R.Kessler
  // Moved from MAIN
  //
  // Jan 15 2018: bug fix to work with CUTWIN; see i+=3

  int  ntmp, i, found ;
  char item[256], line[256];
  char fnam[] = "override_parFile";

  // ---------- BEGIN ------------

  ntmp = 0;
  printf("\n Parse command-line options: \n");
  uniqueOverlap("INIT","SALT2mu command-line override");

  for (i=2; i < argc; ++i) {
    strncpy(item,argv[i],255);
    if (strlen(item) > 200) {
      printf("\n FATAL ERROR in %s: \n", fnam );
      printf("\t Argument %i exceeds 200 characters \n", i  );
      printf("\t check string '%s' \n", item);
      printf("\t ***** ABORT ***** \n");
      fflush(stdout);
      exit(2);
    }
    ntmp++;

    //  if ( strcmp(item,"CUTWIN") == 0 ) {
    if ( !strncmp(item,"CUTWIN",6) ) {  // allow CUTWIN(option)
      // glue together 4 contiguous words into one string
      sprintf(line,"%s %s %s %s", argv[i],argv[i+1],argv[i+2],argv[i+3] ) ;
      found = ppar(line);
      i += 3 ; // Jan 2018 bug fix
    }
    else {
      found = ppar(item);
    }

    if ( found == 0 ) {
      sprintf(c1err,"Invalid command-line override:");
      sprintf(c2err,"%s", item);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  } // end i loop 

  if ( ntmp > 0 ) 
    { printf("\n"); }
  else
    { printf(" None found. Use all options from parameter-input file.\n"); }

} // end of override_parFile


// ********************************************
int ppar(char* item) {

  // Dec 08, 2014 - read new zVARNAME
  // Aug 22, 2016 - call remove_quote
  // Apr 17, 2017 - return(1) if item is found; 0 otherwise
  // Jun 18, 2018 - fix INPUTS.ipar for H11 option
  // Apr 10, 2019 - replace !strcmp with uniqueOverlap function
  //                which aborts on duplicate key.
  
  char *s ;
  int j, ipar, len ;  char tmpString[80], key[40];
  char fnam[] = "ppar" ;

  // --------- BEGIN ----------

  printf(" Parse '%s' \n",item);  fflush(stdout);

  //  if (!strncmp(item,"prefix=",7)) 
  if ( uniqueOverlap(item,"prefix=") )
    { sscanf(&item[7],"%s",INPUTS.PREFIX); return(1); }
  
  //  if (!strncmp(item,"errmask_write=",14)) 
  if ( uniqueOverlap(item,"errmask_write=") )
    { sscanf(&item[14],"%i", &INPUTS.errmask_write ); return(1); }


  //  if (!strncmp(item,"file=",5)) {
  if ( uniqueOverlap(item,"file=") ) {
    s=INPUTS.filename;  
    sscanf(&item[5],"%s",s); remove_quote(s); 
    return(1);
  }

  // -------------------------
  // allow two different keys to define biasCor file name
  //  if ( !strncmp(item,"simfile_bias=",13) ) { 
  if ( uniqueOverlap(item,"simfile_bias=") ) {
    s=INPUTS.simFile_biasCor ;  
    sscanf(&item[13],"%s",s); remove_quote(s); 
    return(1);
  }
  //  if ( !strncmp(item,"simfile_biascor=",16) ) { 
  if ( uniqueOverlap(item,"simfile_biascor=") ) {
    s=INPUTS.simFile_biasCor; 
    sscanf(&item[16],"%s",s); remove_quote(s); 
    return(1);
  }
  // - - - - - - 

  //  if (!strncmp(item,"prescale_biascor=",17)) 
  if ( uniqueOverlap(item,"prescale_biascor=") ) 
    { parse_prescale_biascor(&item[17]); return(1); }
 
  //  if (!strncmp(item,"opt_biascor=",12))
  if ( uniqueOverlap(item,"opt_biascor=")  )
    { sscanf(&item[12],"%d", &INPUTS.opt_biasCor);  return(1); }

  //  if (!strncmp(item,"sigmb_biascor=",14))   
  if ( uniqueOverlap(item,"sigmb_biascor=") )  // legacy name (9/2016)
    { sscanf(&item[14],"%le", &INPUTS.sigint_biasCor); return(1); }

  //  if (!strncmp(item,"sigint_biascor=",15)) 
  if ( uniqueOverlap(item,"sigint_biascor=")  )
    { sscanf(&item[15],"%le", &INPUTS.sigint_biasCor); return(1); }

  //  if (!strncmp(item,"snrmin_sigint_biascor=",22)) 
  if ( uniqueOverlap(item,"snrmin_sigint_biascor=")  )
    { sscanf(&item[22],"%le", &INPUTS.snrmin_sigint_biasCor); return(1); }

  // allow two different keys for fieldgroup_biascor
  //  if (!strncmp(item,"fieldsplit_biascor=",19)) {  
  if ( uniqueOverlap(item,"fieldsplit_biascor=")  ) { // legacy name
    s=INPUTS.fieldGroup_biasCor ;
    sscanf(&item[19],"%s",s); remove_quote(s);
    return(1);
  }


  //  if (!strncmp(item,"fieldgroup_biascor=",19)) { 
  if ( uniqueOverlap(item,"fieldgroup_biascor=")  ) {
    s=INPUTS.fieldGroup_biasCor ;  
    sscanf(&item[19],"%s",s); remove_quote(s);
    return(1);
  }

  //  if (!strncmp(item,"surveygroup_biascor=",20)) { 
  if ( uniqueOverlap(item,"surveygroup_biascor=")  ) {
    s=INPUTS.surveyGroup_biasCor ; 
    sscanf(&item[20],"%s",s); remove_quote(s);
    return(1);
  }

  //  if (!strncmp(item,"surveylist_nobiascor=",21)) { 
  if ( uniqueOverlap(item,"surveylist_nobiascor=")  ) {
    s=INPUTS.surveyList_noBiasCor ;
    sscanf(&item[21],"%s",s); remove_quote(s); 
    return(1);
  }

  //  if (!strncmp(item,"sigma_cell_biascor=",19)) 
  if ( uniqueOverlap(item,"sigma_cell_biascor=") ) 
    { sscanf(&item[19],"%le", &INPUTS.sigma_cell_biasCor); return(1); }
  
  // -------- CC prior ------------
  //  if (!strncmp(item,"simfile_ccprior=",16))  { 
  if ( uniqueOverlap(item,"simfile_ccprior=")  ) {
    s=INPUTS.simFile_CCprior ;
    sscanf(&item[16],"%s",s); remove_quote(s);
    if ( IGNOREFILE(s) ) { simdata_ccprior.USE=0; return(1); }
    simdata_ccprior.USE = 1;
    if ( strcmp(s,"H11") == 0 ) { simdata_ccprior.USEH11 = 1; }

    return(1);
  }

  //  if (!strncmp(item,"maxprobcc_for_sigint=",21)) 
  if ( uniqueOverlap(item,"maxprobcc_for_sigint=") )
    { sscanf(&item[21],"%le", &INPUTS.maxProbCC_for_sigint); return(1); } 

  //  if (!strncmp(item,"varname_pIa=",12)) { 
  if ( uniqueOverlap(item,"varname_pIa=")  ) {
    s=INPUTS.varname_pIa ;
    sscanf(&item[12],"%s",s); remove_quote(s);
    return(1);
  }

  //  if (!strncmp(item,"typeIa_ccprior=",15))  
  if ( uniqueOverlap(item,"typeIa_ccprior=") ) 
    { sscanf(&item[15],"%d", &INPUTS.typeIa_ccprior); return(1); } 

  if ( uniqueOverlap(item,"bins="))  // obsolete, but still allowed (use nzbin)
    { sscanf(&item[5],"%i",&INPUTS.nzbin); return(1); }
  if ( uniqueOverlap(item,"nzbin=")) 
    { sscanf(&item[6],"%i",&INPUTS.nzbin); return(1); }
  if ( uniqueOverlap(item,"nlogzbin=")) 
    { sscanf(&item[9],"%i",&INPUTS.nlogzbin); return(1); }
  if ( uniqueOverlap(item,"powzbin=")) 
    { parse_powzbin(&item[8]); return(1); }  

  if ( uniqueOverlap(item,"blindpar")) 
    { parse_blindpar(item); return(1); }

  if ( uniqueOverlap(item,"min_per_zbin=")) 
    { sscanf(&item[13],"%i", &INPUTS.min_per_zbin); return(1); }

  if ( uniqueOverlap(item,"zVARNAME=")) {
    s=INPUTS.varname_z;    sscanf(&item[9],"%s",s); remove_quote(s); 
    return(1);
  }

  if ( uniqueOverlap(item,"varname_gamma=")) {
    s=INPUTS.varname_gamma;  sscanf(&item[14],"%s",s); remove_quote(s); 
    return(1);
  }
  if ( uniqueOverlap(item,"varname_z=")) { 
    s=INPUTS.varname_z ;  sscanf(&item[10],"%s",s); remove_quote(s); 
    return(1);
  }

  if ( uniqueOverlap(item,"zmin=")) 
    { sscanf(&item[5],"%lf",&INPUTS.zmin) ; return(1); }
  if ( uniqueOverlap(item,"zmax=")) 
    { sscanf(&item[5],"%lf",&INPUTS.zmax) ; return(1); }

  if ( uniqueOverlap(item,"x1min=")) 
    { sscanf(&item[6],"%lf",&INPUTS.x1min); return(1); }
  if ( uniqueOverlap(item,"x1max="))  
    { sscanf(&item[6],"%lf",&INPUTS.x1max); return(1); }

  if ( uniqueOverlap(item,"maxerr_abort_x0="))  
    { sscanf(&item[16],"%lf",&INPUTS.maxerr_abort_x0); return(1); }
  if ( uniqueOverlap(item,"maxerr_abort_x1="))  
    { sscanf(&item[16],"%lf",&INPUTS.maxerr_abort_x1); return(1); }
  if ( uniqueOverlap(item,"maxerr_abort_c="))  
    { sscanf(&item[15],"%lf",&INPUTS.maxerr_abort_c); return(1); }
  
  if ( uniqueOverlap(item,"cmin="))  
    { sscanf(&item[5],"%lf",&INPUTS.cmin); return(1); }
  if ( uniqueOverlap(item,"cmax="))  
    { sscanf(&item[5],"%lf",&INPUTS.cmax); return(1); }

  if ( uniqueOverlap(item,"sig1="))  
    { sscanf(&item[5],"%lf",&INPUTS.sigmB);return(1); }
  if ( uniqueOverlap(item,"sig2="))  
    { sscanf(&item[5],"%lf",&INPUTS.sigx1); return(1); }
  if ( uniqueOverlap(item,"sig3="))  
    { sscanf(&item[5],"%lf",&INPUTS.sigc); return(1); }
  if ( uniqueOverlap(item,"xi01="))  
    { sscanf(&item[5],"%lf",&INPUTS.xi01); return(1); }
  if ( uniqueOverlap(item,"xi0c="))  
    { sscanf(&item[5],"%lf",&INPUTS.xi0c); return(1); }
  if ( uniqueOverlap(item,"xi1c="))  
    { sscanf(&item[5],"%lf",&INPUTS.xi1c); return(1); }

  if ( uniqueOverlap(item,"sigint_fix="))  { 
    s = INPUTS.sigint_fix ;
    sscanf(&item[11],"%s", INPUTS.sigint_fix);  remove_quote(s); 
    return(1); 
  }

  if ( uniqueOverlap(item,"sigmB=")) 
    { sscanf(&item[6],"%lf",&INPUTS.sigmB); return(1); }
  if ( uniqueOverlap(item,"sigmb=")) 
    { sscanf(&item[6],"%lf",&INPUTS.sigmB); return(1); }

  if ( uniqueOverlap(item,"sigx1=")) 
    { sscanf(&item[5],"%lf",&INPUTS.sigx1); return(1); }
  if ( uniqueOverlap(item,"sigc="))  
    { sscanf(&item[5],"%lf",&INPUTS.sigc); return(1); }

  if ( !strncmp(item,"ZPOLY_",6) )  // multiple ZPOLY keys allowed
    { INPUTS.zpolyflag = 1;  parse_ZPOLY_COVMAT(item); return(1); } 

  if ( !strncmp(item,"CUTWIN",6) )  // multiple CUTWIN keys allowed
    { parse_CUTWIN(item); return(1); } 

  if ( uniqueOverlap(item,"idsample_select=") ) 
    {  sscanf(&item[16], "%s", INPUTS.idsample_select );  return(1); } 

  if ( uniqueOverlap(item,"sntype=") )  
    { parse_sntype(&item[7]); return(1); }

  if ( uniqueOverlap(item,"nmax=") ) {
    s = INPUTS.nmaxString ;
    sscanf(&item[5],"%s",s);  remove_quote(s); return(1); } 

  if ( uniqueOverlap(item,"host_logmass_split="))  {
    sprintf(c1err,"Input host_logmass_split is obsolete !");
    sprintf(c2err,"Use p7 and u7 instead.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // read initial values for parameters in fit
  for(ipar=1; ipar < MXCOSPAR; ipar++ ) {
    sprintf(key,"p%d=", ipar); len=strlen(key);
    if ( uniqueOverlap(item,key)) 
      { sscanf(&item[len],"%lf",&INPUTS.parval[ipar]); return(1); }
  }

  // read alternate names for p%d params
  if ( uniqueOverlap(item,"alpha0=")) 
    { sscanf(&item[7],"%lf",&INPUTS.parval[1]); return(1); }
  if ( uniqueOverlap(item,"beta0=")) 
    { sscanf(&item[6],"%lf",&INPUTS.parval[2]); return(1); }
  if ( uniqueOverlap(item,"alpha1=")) 
    { sscanf(&item[7],"%lf",&INPUTS.parval[3]); return(1); }
  if ( uniqueOverlap(item,"beta1=")) 
    { sscanf(&item[6],"%lf",&INPUTS.parval[4]); return(1); }
  if ( uniqueOverlap(item,"OL=")) 
    { sscanf(&item[3],"%lf",&INPUTS.parval[9]); return(1); }
  if ( uniqueOverlap(item,"Ok=")) 
    { sscanf(&item[3],"%lf",&INPUTS.parval[10]); return(1); }
  if ( uniqueOverlap(item,"w="))  
    { sscanf(&item[2],"%lf",&INPUTS.parval[11]); return(1); }
  if ( uniqueOverlap(item,"wa=")) 
    { sscanf(&item[3],"%lf",&INPUTS.parval[12]); return(1); }
  if ( uniqueOverlap(item,"scalePCC=")) 
    { sscanf(&item[9],"%lf",&INPUTS.parval[13]); return(1); }
  if ( uniqueOverlap(item,"sigint=")) 
    { sscanf(&item[7],"%lf",&INPUTS.parval[14]); return(1); }
  if ( uniqueOverlap(item,"alphaHost=")) 
    { sscanf(&item[10],"%lf",&INPUTS.parval[15]); return(1); }
  if ( uniqueOverlap(item,"betaHost=")) 
    { sscanf(&item[9],"%lf",&INPUTS.parval[16]); return(1); }


  /* xxx mark delete xxxxxx
  // check for H11 fitpar options
  for(j=0; j < NPAR_H11_TOT; j++ ) {
    ipar = IPAR_H11 + j;
    sprintf(tmpString,"p%2.2d=", ipar);  len=strlen(tmpString);
    if (!strncmp(item,tmpString,len))  { 
      sscanf(&item[len],"%lf",&INPUTS.parval[ipar]); return(1); 
    }
  }
  xxxxxxxx */

  // ---

  // read initial step sizes for parameters in fit
  for(ipar=1; ipar < MXCOSPAR; ipar++ ) {
    sprintf(key,"s%d=", ipar); len=strlen(key);
    if ( uniqueOverlap(item,key)) 
      { sscanf(&item[len],"%lf",&INPUTS.parstep[ipar]); return(1); }
  }


  // read integer flag to fix or float each param
  for(ipar=1; ipar < MXCOSPAR; ipar++ ) {
    sprintf(key,"u%d=", ipar); len=strlen(key);
    if ( uniqueOverlap(item,key)) 
      { sscanf(&item[len],"%i",&INPUTS.ipar[ipar]); return(1); }
  }

  // read param to fix all params (i.e., use SALT2mu as distance calculator)
  if ( uniqueOverlap(item,"fixpar_all=") ) 
      { sscanf(&item[11],"%i",&INPUTS.fixpar_all ); return(1); }
  
  if ( uniqueOverlap(item,"uM0="))   
    { sscanf(&item[4],"%i",&INPUTS.uM0); return(1); }
  if ( uniqueOverlap(item,"uzsim=")) 
    { sscanf(&item[6],"%i",&INPUTS.uzsim); return(1); }


  /* xxxxxxxxxxxx mark delete Apr 10 2019 xxxxxxxxxxx
  // check for H11 fitpar options
  for(j=0; j < NPAR_H11_TOT; j++ ) {
    ipar = IPAR_H11 + j;
    sprintf(tmpString,"u%2.2d=", ipar);  len=strlen(tmpString);
    if (!strncmp(item,tmpString,len))  { 
      sscanf(&item[len],"%i",&INPUTS.ipar[ipar]); 
      NPAR_H11_USER++ ;
      return(1); 
    }
  }
  xxxxxxxxx end mark xxxxxxxxx */


  if ( uniqueOverlap(item,"blindflag=")) 
    { sscanf(&item[10],"%d",&INPUTS.blindFlag); return(1); }

  if ( uniqueOverlap(item,"h0=")) 
    { sscanf(&item[3],"%lf",&INPUTS.H0); return(1); }
  if ( uniqueOverlap(item,"mag0=")) 
    { sscanf(&item[5],"%lf",&INPUTS.nommag0); return(1); }
  if ( uniqueOverlap(item,"uave="))  
    { sscanf(&item[5],"%i", &INPUTS.uave); return(1); }

  // ------------
  // check peculiar velocity error
  if ( uniqueOverlap(item,"zpecerr="))
    { sscanf(&item[8],"%lf",&INPUTS.zpecerr);  return(1); }

  if ( uniqueOverlap(item,"vpecerr="))  { 
    sscanf(&item[8],"%lf",&INPUTS.zpecerr); 
    INPUTS.zpecerr /= LIGHT_km ;
    return(1); 
  }

  if ( uniqueOverlap(item,"pecv="))   // LEGACY key
    { sscanf(&item[5],"%lf",&INPUTS.zpecerr);  return(1); }

  // lensing term (Sep 2016)
  if ( uniqueOverlap(item,"lensing_zpar=") )
    { sscanf(&item[13],"%lf",&INPUTS.lensing_zpar);  return(1); }

  // ------------
  if ( uniqueOverlap(item,"sig1fit="))   // legacy key from JLM
    { sscanf(&item[8],"%i", &INPUTS.fitflag_sigmb); return(1); }
  if ( uniqueOverlap(item,"fitflag_sigmb="))  // new key from RK (Mar 2015)
    { sscanf(&item[14],"%i", &INPUTS.fitflag_sigmb); return(1); }

  if ( uniqueOverlap(item,"sig1tol=")) // legacy key from JLM
    { sscanf(&item[8],"%lf",&INPUTS.redchi2_tol); return(1); }
  if ( uniqueOverlap(item,"redchi2_tol="))  // new key from RK (Mar 2015)
    { sscanf(&item[12],"%lf",&INPUTS.redchi2_tol); return(1); }

  if ( uniqueOverlap(item,"prescale_simdata="))  
    { sscanf(&item[17],"%lf",&INPUTS.prescale_simData); return(1); }
  if ( uniqueOverlap(item,"prescale_simcc="))  // new key from RK (Dec 2015)
    { sscanf(&item[15],"%lf",&INPUTS.prescale_simCC); return(1); }


  // misc.  
  if ( uniqueOverlap(item,"NDUMPLOG=")) 
    { sscanf(&item[9],"%d", &INPUTS.NDUMPLOG ); return(1); }

  if ( uniqueOverlap(item,"NSPLITRAN=")) 
    { sscanf(&item[10],"%d", &INPUTS.NSPLITRAN); return(1); }

  if ( uniqueOverlap(item,"iflag_duplicate=")) 
    { sscanf(&item[16],"%d", &INPUTS.iflag_duplicate ); return(1); }

  if ( uniqueOverlap(item,"dumpflag_nobiascor=")) 
    { sscanf(&item[19],"%d", &INPUTS.dumpflag_nobiasCor ); return(1); }

  if ( uniqueOverlap(item,"snid_mucovdump=")) 
    { sscanf(&item[15],"%s", INPUTS.SNID_MUCOVDUMP); return(1); }
    
  return(0);
  
} // end ppar


// **************************************************
void parse_nmax(char *item) {

  // Created July 17 2017
  // Parse comma-separated list of NMAX per survey.
  // Examples:
  //  nmax=200                ! limit on total, regardless of survey
  //  nmax=80(SDSS),200(DES)  ! survey-specific nmax
  //  nmax=300,200(DES)       ! nmax=300 for total, and 200 max for DES subset
  //

#define MXARG_nmax 20
  int  i, NARG, nmax, ID ;
  char stringArg[MXARG_nmax][40];
  char *ptrArg[MXARG_nmax];
  char comma[] = "," ;
  char survey[60], tmpString[40] ;
  char fnam[] = "parse_nmax" ;

  // ------------- BEGIN ---------------

  if ( IGNOREFILE(item) ) { return ; }

  for(i=0; i < MXARG_nmax; i++ ) {  ptrArg[i] = stringArg[i]; }

  splitString(item, comma, MXARG_nmax,    // inputs
	      &NARG, ptrArg );            // outputs


  // parse nmax for each survey
  for(i=0; i < NARG; i++ ) {
    sprintf(tmpString, "%s", ptrArg[i] );
    extractStringOpt(tmpString,survey); // modify tmpString & return survey
    sscanf(tmpString, "%d", &nmax );

    if ( strlen(survey) == 0 ) 
      { INPUTS.nmax_tot = nmax; ID = 0; }
    else {
      ID = get_IDSURVEY(survey);
      if ( ID < 0  ){
	sprintf(c1err,"Invalid survey = '%s'", survey) ;
	sprintf(c2err,"Check nmax string in the input file");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
      }
      INPUTS.nmax[ID] = nmax;
    }

    /*
    printf(" xxx '%s' --> nmax=%d and survey='%s'(%d) \n",
    ptrArg[i], nmax, survey, ID ); */
  }
  
  return ;

} // end parse_nmax

// **************************************************
void parse_sntype (char *item) {

  // May 23, 2012
  // parse comma-separated list of types and load
  // INPUTS.Nsntype and INPUTS.sntype[i].
  // Input file example :
  //   sntype= 120,105,104  ! no blank spaces

  char *ptrtok, local_item[100] ;
  int  ITYPE, N;
  //  char fnam[] = "parse_sntype" ;

  // --------------- BEGIN ---------------

  sprintf(INPUTS.sntypeString, "%s", item);
  sprintf(local_item,          "%s", item);

  INPUTS.Nsntype = 0 ;

  N = 0;
  ptrtok = strtok(local_item,",");
  while ( ptrtok != NULL  ) {
    sscanf(ptrtok,"%d", &ITYPE);
    INPUTS.sntype[N]   = ITYPE ;

    if ( ITYPE == -1 ) { return ; } // flag to process all types.

    N++ ;  
    ptrtok = strtok(NULL, ",");
  }


  INPUTS.Nsntype = N ;
  
} // end of parse_sntype

// **************************************************
void parse_powzbin(char *item) {

  // parse powzbin arguments: either 1 arg, or 2 comma-separated args
  // First arg is powzbin such that bin size is proportional to 
  //  (1+z)^powzbin.
  // Second arg (optional) is znhalf, the redshift where half the
  // z-bins obey powzbin, and the remaining half have constant bins.
  //

  int  i, NARG, MXARG=3;
  char stringArg[2][20];
  char *ptrArg[2] = { stringArg[0], stringArg[1] } ;
  char comma[] = "," ;
  char fnam[] = "parse_powzbin" ;

  // ------------- BEGIN ---------------
  
  splitString(item, comma, MXARG,    // inputs
	      &NARG, ptrArg );       // outputs

  sscanf(ptrArg[0], "%le", &INPUTS.powzbin); 
  if ( NARG == 2 ) 
    { sscanf(ptrArg[1], "%le", &INPUTS.znhalf ); }

  /* xxxxxxxxxxxx
  printf(" xxx %s: powzbin = %f \n", fnam, INPUTS.powzbin );
  printf(" xxx %s: znhalf  = %f \n", fnam, INPUTS.znhalf );
  fflush(stdout);
  debugexit(fnam); // xxxx REMOVE
  xxxxxxx */

  return ;

} // end parse_powzbin


// **************************************************
void parse_blindpar(char *item) {

  // Created Aug 2017
  // Parse 2 parameters separate by comma.
  // E.g., blindpar9=.15,343.2 --> blind offset is 0.15*cos(343.2)

  int  ipar, i, NARG, MXARG=3;
  char stringArg[2][20], item_local[60] ;
  char *ptrArg[2] = { stringArg[0], stringArg[1] } ;
  char comma[] = "," ;
  char fnam[] = "parse_blindpar" ;

  // ------------- BEGIN ---------------

  if (!strncmp(item,"blindpar9=",10)) 
    { sprintf(item_local, "%s", &item[10] ); ipar=9; }
  else if (!strncmp(item,"blindpar11=",11)) 
    { sprintf(item_local, "%s", &item[11] ); ipar=11; }
  else {
    sprintf(c1err,"Invalid input '%s'", item);
    sprintf(c2err,"Check valid blindpar keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  splitString(item_local, comma, MXARG,    // inputs
	      &NARG, ptrArg );            // outputs

  sscanf(ptrArg[0], "%le", &INPUTS.blindpar[ipar][0] ); 
  sscanf(ptrArg[1], "%le", &INPUTS.blindpar[ipar][1] ); 

  /* xxx
  printf(" xxx %s: blindpar[%d] = %f, %f \n", fnam, ipar, 
	 INPUTS.blindpar[ipar][0], INPUTS.blindpar[ipar][1] );
	 debugexit(fnam); 
  xxxx */

  return ;

} // end parse_blindpar


// **************************************************
void parse_prescale_biascor (char *item) {

  // June 25, 2016
  // parse comma-separated prescale of the form
  //
  //   prescale_simbias=3,4
  //
  //   --> prescale by 4 and use 3rd subset.
  //   --> INPUTS.prescale_simbias[0] = 3
  //   --> INPUTS.prescale_simbias[1] = 4
  //
  //  Note that PS0 can have values 0,1,2 ... PS1-1
  //

  char *ptrtok, local_item[100] ;
  int  N, PS;
  char fnam[] = "parse_prescale_biascor" ;

  // --------------- BEGIN ---------------

  sprintf(local_item,          "%s", item);

  N = 0;
  ptrtok = strtok(local_item,",");
  while ( ptrtok != NULL  ) {
    sscanf(ptrtok,"%d", &PS);
    INPUTS.prescale_biasCor[N]   = PS ;
    N++ ;  
    ptrtok = strtok(NULL, ",");
  }
  
  int PS0 = INPUTS.prescale_biasCor[0] ;
  int PS1 = INPUTS.prescale_biasCor[1] ;

  if ( PS0 < 0 || PS0 >= PS1 ) {
    sprintf(c1err,"Invalid prescale_biascor=%d,%d", PS0, PS1);
    sprintf(c2err,"Check SALT2mu input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

} // end of parse_prescale_biascor

// **************************************************
void parse_IDSAMPLE_SELECT(char *item) {

  // Aug 31 2016
  // parse plus-separated list of IDSAMPLEs to fit.
  // Example:
  //   item = '0+2' --> 
  //

  int  NTMP, i, ID ; 
  char itemLocal[60], *ptrID[MXNUM_SAMPLE], strID[MXNUM_SAMPLE][4] ;
  char plus[] = "+" ;
  char fnam[] = "parse_IDSAMPLE_SELECT" ;

  // --------- BEGIN -----------

  if ( strlen(item) == 0 ) { return; }

  sprintf(itemLocal, "%s", item );
  for(i=0; i < MXNUM_SAMPLE; i++ ) { ptrID[i] = strID[i]; }

  if ( strstr(itemLocal,",") != NULL ) {
    sprintf(c1err,"Illegal comma delimeter in idsample_select key");
    sprintf(c2err,"Try plus (+) delimeter instead.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  splitString(itemLocal, plus, MXNUM_SAMPLE,      // inputs
	      &NTMP, ptrID );                   // outputs

  // reset all DOFLAGs to zero
  for(ID=0; ID < MXNUM_SAMPLE; ID++ )  { 
    SAMPLE_BIASCOR[ID].DOFLAG_SELECT  = 0; 
    SAMPLE_BIASCOR[ID].DOFLAG_BIASCOR = 0; 
  }

  for(i=0; i < NTMP; i++ ) {
    sscanf(ptrID[i], "%d", &ID); 
    SAMPLE_BIASCOR[ID].DOFLAG_SELECT  = 1; 
    SAMPLE_BIASCOR[ID].DOFLAG_BIASCOR = 1; 
  }

  printf(" Will select %d IDSAMPLEs : %s \n",  NTMP,  itemLocal);
  fflush(stdout);

  return ;

} // end parse_IDSAMPLE_SELECT


// **************************************************
void parse_sigint_fix(char *item) {

  // Created July 5 2018
  // parse sigint_fix argument and load SNDATA_INFO.sigint_fix[IDAMPLE]
  //
  
  // Example
  //   item  = '0.106' 
  //       --> sigint_fix = 0.106 for all IDSAMPLEs.
  //
  //   item  = '0.106,0.988,0.904 
  //      --> 3 sigint_fix values for 3 IDSAMPLEs
  //


  int  NSAMPLE = NSAMPLE_BIASCOR ;
  int  idsample, Nsigint, i ;
  double sigint;
  char *name, comma[] = "," ;
  char fnam[] = "parse_sigint_fix";
  char itemLocal[200], *ptrSIG[MXNUM_SAMPLE], strSIG[MXNUM_SAMPLE][8] ;

  // --------- BEGIN -----------

  if ( strlen(item) == 0 ) { return ;}

  // check for comma
  if ( strstr(item,comma) == NULL ) {
    // no comma --> fix same sigint for all IDSAMPLE
    sscanf(item, "%le", &sigint) ;
    for(idsample=0; idsample<NSAMPLE; idsample++ ) 
      { SNDATA_INFO.sigint_fix[idsample] = sigint; }
  }
  else {
    // strip sigint for each IDSAMPLE
    for(i=0; i < MXNUM_SAMPLE; i++ ) { ptrSIG[i] = strSIG[i]; }
    sprintf(itemLocal,"%s", item);
    splitString(itemLocal, comma, MXNUM_SAMPLE,      // inputs
		&Nsigint, ptrSIG );                    // outputs
    if ( Nsigint != NSAMPLE ) {
      sprintf(c1err,"Nsiginit=%d != N_IDSAMPLE=%d", Nsigint, NSAMPLE);
      sprintf(c2err,"sigint_fix=%s must include each IDSAMPLE",	 item) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }


    for(idsample=0; idsample < Nsigint; idsample++ ) {
      sscanf(strSIG[idsample], "%le", &sigint) ;
      SNDATA_INFO.sigint_fix[idsample] = sigint; 
    }

    
  }


  // - - - - -
  // turn off other sigint inputs
  INPUTS.sigmB = 0.0 ;  
  
  // ----------
  // print sigint_fix for each IDSAMPLE, and compute sqsigint_fix
  
  printf("\n");
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    name   = SAMPLE_BIASCOR[idsample].NAME ;
    sigint = SNDATA_INFO.sigint_fix[idsample];
    SNDATA_INFO.sqsigint_fix[idsample] = sigint*sigint ;

    printf("  sigint_fix=%6.4f for IDSAMPLE=%d (%s)\n",
	   sigint, idsample, name);
  }
  fflush(stdout);
    
  
  return ;

} // end parse_sigint_fix

// **************************************************
void parse_CUTWIN(char *item) {

  // Created 4/24/2012 by R.Kessler
  // parse line beginning with 
  //   CUTWIN <VARNAME>  <MIN> <MAX>
  // and fill
  //   INPUTS.NCUTWIN
  //   INPUTS.CUTWIN_NAME
  //   INPUTS.CUTWIN_RANGE
  //   INPUTS.CUTWIN_ABORTFLAG   ! 9.29.2016

  char  *ptrtok, local_item[200], stringOpt[60], KEY[60] ;
  int   ICUT, i ;
  char  fnam[] = "parse_CUTWIN" ;

  // ---------- BEGIN ------------

  INPUTS.NCUTWIN++ ;
  ICUT = INPUTS.NCUTWIN-1 ;
  INPUTS.CUTWIN_ABORTFLAG[ICUT] = 1;  // default is to abort on missing var
  INPUTS.CUTWIN_DATAONLY[ICUT]  = 0;  // default cut on data and biasCor

  sprintf(local_item,"%s", item);
  ptrtok = strtok(local_item," ");

  for ( i=0; i < 4 ; i++ ) {

    if ( i == 0 ) {
      // check for option in CUTWIN(option)
      sscanf ( ptrtok, "%s", KEY ); 
      extractStringOpt(KEY, stringOpt); // return stringOpt

      if ( strcmp(stringOpt,"NOABORT") == 0 ) 
	{ INPUTS.CUTWIN_ABORTFLAG[ICUT] = 0; } // allow missing variable 

      if ( strcmp(stringOpt,"DATAONLY") == 0 ) 
	{ INPUTS.CUTWIN_DATAONLY[ICUT] = 1; } // cut on data only
    }

    if ( i == 1 ) 
      { sscanf ( ptrtok, "%s", INPUTS.CUTWIN_NAME[ICUT] ); } 
        
    if ( i == 2 ) 
      { sscanf ( ptrtok, "%le", &INPUTS.CUTWIN_RANGE[ICUT][0] ); } 

    if ( i == 3 ) 
      { sscanf ( ptrtok, "%le", &INPUTS.CUTWIN_RANGE[ICUT][1] ); } 

    ptrtok = strtok(NULL, " ");
  }

  printf("\t Will apply CUTWIN on %12s from %10.4f to %10.4f  (ABORTFLAG=%d)\n"
	 ,INPUTS.CUTWIN_NAME[ICUT]
	 ,INPUTS.CUTWIN_RANGE[ICUT][0]
	 ,INPUTS.CUTWIN_RANGE[ICUT][1]
	 ,INPUTS.CUTWIN_ABORTFLAG[ICUT]
	 ) ;

} // end of parse_CUTWIN


// **************************************************
int apply_CUTWIN(int OPT, int *DOFLAG_CUTWIN, double *CUTVAL_LIST) {

  // Created Jan 2016
  // Returns TRUE if all CUTVAL_LIST values pass CUTWIN cuts.
  // OPT=1 for data    --> check all
  // OPT=2 for simFile --> check all except for SIM_TYPE_INDEX
  //
  // 9/29/2016: pass new array DOFLAG_CUTWIN to determine which
  //            cuts to apply or ignore
  //

  int LDMP = 0 ; // (OPT==666);
  int icut, LCUT ;
  char *cnam;
  char fnam[] = "apply_CUTWIN" ;

  // ------------- BEGIN -----------

  if ( LDMP ) { printf(" xxx --------------------------- \n"); }

  LCUT = 1 ;    
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
   
    if ( DOFLAG_CUTWIN[icut] == 0 ) { continue; }

    if ( LDMP ) {
      printf(" xxx cut on %s = %f  (cutwin=%.3f to %.3f, OPT=%d)\n",
	     INPUTS.CUTWIN_NAME[icut], CUTVAL_LIST[icut],
	     INPUTS.CUTWIN_RANGE[icut][0], 
	     INPUTS.CUTWIN_RANGE[icut][1], OPT ); 
    }

    // check SIM-option to skip SIM_TYPE_INDEX
    cnam   = INPUTS.CUTWIN_NAME[icut] ;
    if( OPT == 2 && usesim_CUTWIN(cnam)==0 ) { continue ; }

    if ( CUTVAL_LIST[icut] < INPUTS.CUTWIN_RANGE[icut][0] ) 
      { LCUT = 0 ; goto DONE ; }

    if ( CUTVAL_LIST[icut] > INPUTS.CUTWIN_RANGE[icut][1] ) 
      { LCUT = 0 ; goto DONE ; }
  }


 DONE:

  if (LDMP ) 
    { printf(" xxx ---> LCUT=%d\n", LCUT); fflush(stdout); debugexit(fnam); } 

  return(LCUT);


} // end apply_CUTWIN

// **************************************************
int set_DOFLAG_CUTWIN(int ivar, int icut, int isData) {

  // return 1 if ivar>=0.
  // If ivar<0 and abortflag is set for icut, then abort.
  // 
  // Oct 23 2018: 
  //  + new input arg isData=1 for data, zero for sim biasCor 
  //  + check DATAONLY flag.

  int  NOVAR     = ( ivar < 0 );
  int  ABORTFLAG = INPUTS.CUTWIN_ABORTFLAG[icut];
  int  DATAONLY  = INPUTS.CUTWIN_DATAONLY[icut];
  int  ISTAT;
  char *VARNAME  = INPUTS.CUTWIN_NAME[icut];
  char  fnam[] = "set_DOFLAG_CUTWIN" ;

  // ------------- BEGIN -------------
  
  if ( DATAONLY && isData==0 ) { return(0) ; } // Oct 23 2018

  if ( NOVAR && ABORTFLAG ) {
    sprintf(c1err,"Invalid CUTWIN on '%s'", VARNAME);
    sprintf(c2err,"Check CUTWIN keys in input file" ); 
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NOVAR )
    { ISTAT = 0 ; }
  else
    { ISTAT = 1 ; }

  return(ISTAT);

} // end set_DOFLAG_CUTWIN


// **************************************************
int usesim_CUTWIN(char *varName) {

  // Created Jan 2016
  // Return TRUE if input *varName should be used for
  // simFiles when checking the CUTWIN list.
  // Note that this function does not apply a cut,
  // but determines if a cut should be applied.

  int use ;
  use = 1; 
  if (strcmp(varName,"SIM_TYPE_INDEX" ) == 0 ) { use = 0 ; }  
  if (strcmp(varName,"SIM_NONIA_INDEX") == 0 ) { use = 0 ; }  
  if (strcmp(varName,"SIM_NON1A_INDEX") == 0 ) { use = 0 ; }  // 7.31.2018
  if (strcmp(varName,"SIM_TEMPLATE_INDEX") == 0 ) { use = 0 ; }  // 7.31.2018
  return(use);
}

// **************************************************
void  SPLITRAN_errmask(void) {

  // Created July 2012 by R.Kessler
  // set data[n].errmask for SPLITRAN option.

  int n, snid ;

  if ( INPUTS.NSPLITRAN <= 1 ) { return ; }

  for (n=0; n<FITINP.NSNCUTS; ++n)  {

    snid = data[n].snid ;

    // subtract this errbit if it was set earlier
    data[n].errmask -= (ERRMASK_SPLITRAN & data[n].errmask) ;

    // check random-subset option
    if ( SPLITRAN_ACCEPT(n,snid) == 0 ) { 
      data[n].errmask += ERRMASK_SPLITRAN ; 
    }
  }

} // end of SPLITRAN_errmask

// **************************************************
int SPLITRAN_ACCEPT(int isn, int snid) {

  // Created july 2012 by R.Kessler
  // return 1 if this isn is accepted; 0 otherwise
  // Use integer SNID so that the same SN are used
  // for analyses with slightly different cuts.
  //
  
  int jj;

  if ( INPUTS.NSPLITRAN <= 1 ) { return 1; }

  //jj = (isn % INPUTS.NSPLITRAN) + 1 ;
  jj = (snid % INPUTS.NSPLITRAN) + 1 ;

  if ( jj == NJOB_SPLITRAN ) 
    { return 1; }
  else
    { return 0; }

} // end of SPLITRAN_ACCEPT


// **************************************************
void SPLITRAN_SUMMARY(void) {

  // Created July 10, 2012 by R.Kessler
  // if SPLIT-RAN option is used, compute MEAN and RMS
  // for each parameter, and also compare with avgerage
  // fit-error.
  //
  // Dec 2 2016: compute & print rms of sample size

  
  int ipar, isplit, NSPLIT, NPAR, NERR_VALID ;
  double 
    VAL, ERR, XN, XNTMP, SQRMS, NSNAVG, NSNRMS
    ,SUMVAL1[MAXPAR], SUMERR1[MAXPAR]
    ,SUMVAL2[MAXPAR], SUMERR2[MAXPAR]
    ,AVG_VAL[MAXPAR], RMS_VAL[MAXPAR]
    ,AVG_ERR[MAXPAR], RMS_ERR[MAXPAR]
    ,SUMN, SQSUMN
    ;
  
  char fnam[] = "SPLITRAN_SUMMARY" ;

  // -------------- BEGIN -------------

  NSPLIT = INPUTS.NSPLITRAN ;
  NPAR   = FITINP.NFITPAR_ALL ;

  if ( NSPLIT <= 1 ) { return ; }

  XN     = (double)NSPLIT ;

  // get average sample size & RMS
  SUMN = SQSUMN = 0.0 ;
  for ( isplit=1; isplit <= NSPLIT; isplit++ ) {
    XNTMP = (double)FITRESULT.NSNFIT_SPLITRAN[isplit] ;
    SUMN   += XNTMP ;
    SQSUMN += (XNTMP*XNTMP) ;
  } 
  NSNAVG = SUMN/XN ;
  NSNRMS = RMSfromSUMS(NSPLIT, SUMN, SQSUMN) ;

  for ( ipar=1; ipar <= NPAR ; ipar++ ) {
    SUMVAL1[ipar] = 0.0 ;
    SUMERR1[ipar] = 0.0 ;
    SUMVAL2[ipar] = 0.0 ;
    SUMERR2[ipar] = 0.0 ;

    AVG_VAL[ipar] =  RMS_ERR[ipar] = 0.0 ;
    AVG_ERR[ipar] = -9.0 ;

    NERR_VALID = 0 ;
    
    for ( isplit=1; isplit <= NSPLIT; isplit++ ) {
      VAL = FITRESULT.PARVAL[isplit][ipar] ;
      ERR = FITRESULT.PARERR[isplit][ipar] ;
      if ( ERR > 0.0 ) {
	NERR_VALID++ ; 
	SUMVAL1[ipar] += VAL;
	SUMERR1[ipar] += ERR;
	SUMVAL2[ipar] += (VAL*VAL) ;
	SUMERR2[ipar] += (ERR*ERR) ;
      }
    }

    if ( NERR_VALID == 0 ) { continue ; }
    XN = (double)NERR_VALID ;
    
    AVG_VAL[ipar] = SUMVAL1[ipar]/XN ;
    AVG_ERR[ipar] = SUMERR1[ipar]/XN ;
    RMS_VAL[ipar] = RMSfromSUMS(NERR_VALID, SUMVAL1[ipar], SUMVAL2[ipar] );
    RMS_ERR[ipar] = RMSfromSUMS(NERR_VALID, SUMERR1[ipar], SUMERR2[ipar] );

  }  // end of ipar loop


  // ---- print results -----
  sleep(2);
  fflush(stdout);

  // write summary to outFile. 
  char OUTFILE[MXPATHLEN] ;
  FILE *fp;
  sprintf(OUTFILE,"%s_summary.out", INPUTS.PREFIX );
  fp = fopen(OUTFILE,"wt");
  if ( !fp )  {
    sprintf(c1err,"Could not open SPLITRAN summary file:");
    sprintf(c2err,"%s", OUTFILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n SPLITRAN SUMMARY is in %s \n", OUTFILE);

  //.xyz
  fprintf(fp, "# SPLITRAN SUMMARY for %d jobs (Avg sample size: %d +- %d) \n",
	 NJOB_SPLITRAN, (int)NSNAVG, (int)NSNRMS );
  fprintf(fp," \n");
  fprintf(fp,"VARNAMES: ROW PARNAME AVG_VAL AVG_ERR RMS_VAL RMS_ERR \n" );

  for ( ipar=1; ipar <= NPAR ; ipar++ ) {
    if ( AVG_ERR[ipar] <= 0.0 ) { continue ; }
    fprintf(fp, "ROW: %2d  %-10s  %8.3f  %8.3f  %8.3f %8.3f \n"
	   ,ipar
	   ,FITRESULT.PARNAME[ipar]
	   ,AVG_VAL[ipar]
	   ,AVG_ERR[ipar]
	   ,RMS_VAL[ipar]
	   ,RMS_ERR[ipar]
	   );
  }
  fclose(fp);


  /* xxxxxxx mark delete Apr 29 2019 xxxxxxxx
  printf("\n SPLITRAN SUMMARY for %d jobs (Avg sample size: %d +- %d) \n",
	 NJOB_SPLITRAN, (int)NSNAVG, (int)NSNRMS );
  printf(" \n");
  printf("                   Average   Average    RMS      RMS   \n" );
  printf(" Param             Value     Error      Value    Error \n" );
  printf(" ------------------------------------------------------ \n");

  for ( ipar=1; ipar <= NPAR ; ipar++ ) {
    if ( AVG_ERR[ipar] <= 0.0 ) { continue ; }
    printf(" %-14s  %8.3f  %8.3f  %8.3f %8.3f \n"
	   ,FITRESULT.PARNAME[ipar]
	   ,AVG_VAL[ipar]
	   ,AVG_ERR[ipar]
	   ,RMS_VAL[ipar]
	   ,RMS_ERR[ipar]
	   );
  }
  printf(" ------------------------------------------------------ \n");
  fflush(stdout);
  xxxxxxxxxxxxx end mark xxxxxxx */

}  // end of SPLITRAN_SUMMARY


// ======================================
void  CPU_SUMMARY(void) {

  // Created Nov 22 2017
  printf("\n PROCESS-TIME SUMMARY: \n");

  printf("   Ini-time:  %.2f minutes \n",
	 (t_end_init-t_start)/60.0 ) ;

  if ( INPUTS.opt_biasCor > 0 ) {
    printf("\t BiasCor read time %.2f minutes \n",
	   (t_read_biasCor[1]- t_read_biasCor[0])/60.0 );
    //	   (t_read_biasCor[2]- t_read_biasCor[1])/60.0 );
  }

  printf("   Fit-time:  %.2f minutes \n",
	 (t_end_fit-t_end_init)/60.0 ) ;

  printf("\n");
  fflush(stdout);

} // end CPU_SUMMARY

// **************************************************
void parse_ZPOLY_COVMAT(char *item) {

  // Created 4.23.2012 by R.Kessler
  // Parse string of the form
  //   ZPOLY_[covterm] = IDSURVEY a0 a1 a2 a3 a4
  // and fill INPUTS_ZPOLY_COVMAT.XXX array

  int order, i, IDSURVEY, LKEY, ISP_SURVEY, NTMP ;
  double tmp[50] ;
  double *ptr_ZPOLY ;

  char  local_item[200], KEY[60], key[60], *ptrtok ;
  char  fnam[] = "parse_ZPOLY_COVMAT" ;

  // ----------- BEGIN ----------


  sprintf(local_item,"%s", item);
  ptrtok = strtok(local_item," ");

  for ( i=0; i <= ORDER_ZPOLY_COVMAT+2; i++ ) {

    if ( i == 0 ) 
      { sscanf ( ptrtok, "%s", KEY ); } // ZPOLY_XXX
    else if ( i == 1 ) 
      { sscanf ( ptrtok, "%d", &IDSURVEY ); } 
    else {
      order = i-2;
      sscanf ( ptrtok, "%le", &tmp[order] ); 
    } 
    
    ptrtok = strtok(NULL, " ");
  }

  // translate absolute IDSURYVEY into sparse ISP_SURVEY=0,1,2 ...

  NTMP = INPUTS_ZPOLY_COVMAT.NSURVEY ;
  ISP_SURVEY = -9;
  for ( i=0; i < NTMP; i++ ) {
    if ( IDSURVEY == INPUTS_ZPOLY_COVMAT.IDSURVEY_LIST[i] ) {
      ISP_SURVEY = i;
    }
  }

  // if we have not found IDSURVEY above, then increment sparse index
  if ( ISP_SURVEY < 0 ) { 
    INPUTS_ZPOLY_COVMAT.NSURVEY++ ;
    ISP_SURVEY = INPUTS_ZPOLY_COVMAT.NSURVEY-1 ;
    INPUTS_ZPOLY_COVMAT.IDSURVEY_LIST[ISP_SURVEY] = IDSURVEY ;
  }

  // strip off part of KEY that is after ZPOLY but before '='
  // i.e, if KEY  = 'ZPOLY_sigmB=' then key = 'sigmB'
  LKEY = strlen(KEY);
  key[0] = 0 ;
  ptr_ZPOLY = &tmp[0] ; // anything to avoid compile error

  for ( i=6; i < LKEY-1; i++ ) 
    { sprintf(key,"%s%c", key, KEY[i] );  }

  if ( strcmp(key,"sigmB") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.sigmB[ISP_SURVEY] ; }
  else if ( strcmp(key,"sigx1") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.sigx1[ISP_SURVEY] ; }
  else if ( strcmp(key,"sigc") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.sigc[ISP_SURVEY] ; }
  else if ( strcmp(key,"xi01") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.xi01[ISP_SURVEY] ; }
  else if ( strcmp(key,"xi0c") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.xi0c[ISP_SURVEY] ; }
  else if ( strcmp(key,"xi1c") == 0 ) 
    { ptr_ZPOLY = INPUTS_ZPOLY_COVMAT.xi1c[ISP_SURVEY] ; }
  else {
    sprintf(c1err,"Could not parse '%s' (key=%s)", KEY, key );
    sprintf(c2err,"Check ZPOLY keys in SALT2mu input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // load polynomial parameters
  for ( i=0; i <= ORDER_ZPOLY_COVMAT; i++ ) 
    {  *(ptr_ZPOLY+i) = tmp[i] ; }

  /*
  printf("\t xxx KEY = %s(%s)   |  IDSURVEY=%d  ISP_SURVEY=%d\n", 
	 KEY, key, IDSURVEY, ISP_SURVEY);
  printf("\t xxx a0-a3 = %f,%f,%f,%f \n",
	 tmp[0], tmp[1], tmp[2], tmp[3] );
  */

} // end of parse_ZPOLY


// *********************************************
void  load_ZPOLY_COVMAT(int IDSURVEY, double Z ) {

  // Created 4/23/2012 by R.Kessler
  // Load ZPOLY_COVMAT scatter matrix for this IDSURYEY and redshift Z.
  // Bail if ZPOLY params are not stored for this IDSURVEY.

  double sigmB, sigx1, sigc, xi01, xi0c, xi1c ;
  int i, l,m, ISP_SURVEY;
  //  char fnam[] = "load_ZPOLY_COVMAT" ;

  // ------------ BEGIN -----------

  // init matrix to zero
   for (l=0; l < 3; l++) {
     for (m=0; m<3; m++) { 
       FITINP.ZPOLY_COVMAT[l][m] = 0.0 ;
     }
   }

  // first check of this IDSURVEY is defined in the 
  // INPUTS_ZPOLY_COVMAT struct
  ISP_SURVEY = -9;
  for ( i=0; i < INPUTS_ZPOLY_COVMAT.NSURVEY; i++ ) {
    if ( IDSURVEY == INPUTS_ZPOLY_COVMAT.IDSURVEY_LIST[i] ) 
      {  ISP_SURVEY = i; }    
  }


  if ( ISP_SURVEY < 0 ) { return ; }

  sigmB = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.sigmB[ISP_SURVEY] );
  sigx1 = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.sigx1[ISP_SURVEY] );
  sigc  = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.sigc[ISP_SURVEY] );
  xi01  = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.xi01[ISP_SURVEY] );
  xi0c  = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.xi0c[ISP_SURVEY] );
  xi1c  = sum_ZPOLY_COVMAT(Z, INPUTS_ZPOLY_COVMAT.xi1c[ISP_SURVEY] );  
 
  FITINP.ZPOLY_COVMAT[INDEX_mB][INDEX_mB] = sigmB*sigmB;
  FITINP.ZPOLY_COVMAT[INDEX_x1][INDEX_x1] = sigx1*sigx1;
  FITINP.ZPOLY_COVMAT[INDEX_c][INDEX_c]   = sigc*sigc;
  FITINP.ZPOLY_COVMAT[INDEX_mB][INDEX_x1] = xi01*(sigmB*sigx1);
  FITINP.ZPOLY_COVMAT[INDEX_mB][INDEX_c]  = xi0c*(sigmB*sigc);
  FITINP.ZPOLY_COVMAT[INDEX_x1][INDEX_mB] = 
    FITINP.ZPOLY_COVMAT[INDEX_mB][INDEX_x1];
  FITINP.ZPOLY_COVMAT[INDEX_x1][INDEX_c]  = xi1c*(sigx1*sigc);
  FITINP.ZPOLY_COVMAT[INDEX_c][INDEX_mB]  = 
    FITINP.ZPOLY_COVMAT[INDEX_mB][INDEX_c] ;
  FITINP.ZPOLY_COVMAT[INDEX_c][INDEX_x1]  = 
    FITINP.ZPOLY_COVMAT[INDEX_x1][INDEX_c] ;

} // end of load_ZPOLY_COVMAT


double sum_ZPOLY_COVMAT(double Z, double *polyPar) {
  int i;
  double xi, sum, par ;
  //  char fnam[] = "sum_ZPOLY_COVMAT" ;
  // ----------- BEGIN -----------------

  sum = 0.0 ;
  for ( i=0; i <= ORDER_ZPOLY_COVMAT; i++ ) {
    xi   = (double)i;
    par  = polyPar[i];
    sum += (par * pow(Z,xi) );
  }
  return sum ;
} // end of sum_ZPOLY_COVMAT


// *********************************************
void prep_input(void) {

  // May 9 2019: check INPUTS.fixpar_all

  int i,  NFITPAR, NTMP=0 ;
  double *blindpar ;
  char usage[10];
  char fnam[] = "prep_input";

  // ------------ BEGIN -----------

  // check option to fix all parameters and use SALT2mu 
  // as a distance calculator
  if ( INPUTS.fixpar_all ) {
    printf("\n FIXPAR_ALL option -> "
	   "Compute MU from initial values (no floated params)\n\n" );
    fflush(stdout);
	   
    for (i=0; i < MAXPAR; i++)  { INPUTS.ipar[i]=0; }
    INPUTS.uM0           = 0 ;
    INPUTS.fitflag_sigmb = 0 ;
    INPUTS.nzbin         = 1 ;
    INPUTS.nlogzbin      = 0 ;
  }

  // make sure that either nzbin or nlogzbin is specified, but not both
  if ( INPUTS.nzbin    > 0 ) { NTMP++; }
  if ( INPUTS.nlogzbin > 0 ) { NTMP++; }
  if ( NTMP != 1 ) {
    sprintf(c1err,"Error specifying nzbin=%d and nlogzbin=%d",
	    INPUTS.nzbin, INPUTS.nlogzbin );
    sprintf(c2err,"One of these must be specified (not both).");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( INPUTS.nlogzbin > 0 ) { INPUTS.nzbin = INPUTS.nlogzbin ; }

  NFITPAR = INPUTS.nzbin + MXCOSPAR; // estimate
  if ( NFITPAR >= MAXPAR_MINUIT ) {
    sprintf(c1err,"NFITPAR=%d (%d z bins + %d params)",
	    NFITPAR, INPUTS.nzbin, MXCOSPAR );
    sprintf(c2err,"But MAXPAR_MINUIT = %d", MAXPAR_MINUIT );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  NSIMCC = NSIMDATA = 0 ;
  NJOB_SPLITRAN = 0;

  if ( strlen(INPUTS.sigint_fix) > 0 ) { INPUTS.sigmB = 0.0 ; }
  
  FITINP.COVINT_PARAM_FIX   = INPUTS.sigmB ; // Mar 2016
  FITINP.COVINT_PARAM_LAST  = INPUTS.sigmB ; 
  INPUTS.covint_param_step1 = INPUTS.sigmb_step1 ; // default COVINT param

  // Oct 9 2018: for sigInt(IDSAMPLE), vary global COV scale instead
  //             of varying sigInt.
  if ( strlen(INPUTS.sigint_fix) > 0 ) {
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ; 
  }

    
  printf("zmin =%f zmax =%f bins=%i \n",
	 INPUTS.zmin, INPUTS.zmax, INPUTS.nzbin);

  if ( INPUTS.nzbin >= MAXBIN_z ) {
    sprintf(c1err,"nzbin=%d exceeds bound of MAXBIN_z=%d", 
	    INPUTS.nzbin, MAXBIN_z);
    sprintf(c2err,"Check 'bins=' key in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  if (INPUTS.Nsntype <= 0 ) 
    { printf("SN selection type=ALL\n"); }
  else 
    { printf("SN selection sntype= %s\n", INPUTS.sntypeString); }


  printf("x1min=%f x1max=%f \n", INPUTS.x1min, INPUTS.x1max);
  printf("cmin =%f cmax =%f \n", INPUTS.cmin, INPUTS.cmax);
  printf("H0=%f \n", INPUTS.H0);
  printf("Nominal M0=%f \n", INPUTS.nommag0);
  printf("zpecerr = %.4f \n", INPUTS.zpecerr );
  printf("dsigint/dz(lensing) =%4f \n", INPUTS.lensing_zpar );

  if (INPUTS.uave) 
    { printf("Use average M0.\n"); }
  else 
    { printf("Use nominal M0.\n"); }
  if ( INPUTS.fitflag_sigmb )  { 
    printf("Will fit for sig1 such that chi2/N ~ 1.0\n");
    printf("Will use sig1tol=%f . \n", INPUTS.redchi2_tol);
    printf("Will allow max sigmb iterations=%d .\n", MAXFITITER); 
  }
  

  printf("uM0=%d --> ", INPUTS.uM0 ) ;
  if ( INPUTS.uM0 == M0FITFLAG_CONSTANT ) 
    { printf("M0 fixed. \n"); }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_FLAT ) 
    { printf("M0 floated in each z bin. \n"); }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_INTERP ) 
    { printf("M0 floated in each z bin,and interpolated. \n"); }
  else
    { printf("M0 option not defined ???\n"); }

  if ( INPUTS.uzsim ) { printf("REDSHIFT CHEAT: z -> simz \n");  }

  fflush(stdout);

  if ( INPUTS.prescale_simData > 1.0 ) {
    printf("PRE-SCALE SIMDATA by %.1f \n", INPUTS.prescale_simData);
  }
  if ( INPUTS.prescale_simCC > 1.0 ) {
    printf("PRE-SCALE SIMCC by %.1f \n", INPUTS.prescale_simCC);
  }

  char *f_biasCor = INPUTS.simFile_biasCor;
  char *f_CCprior = INPUTS.simFile_CCprior;
  int  ISFILE_BIASCOR = ( IGNOREFILE(f_biasCor)==0 );

  // check for default biasCor option
  if ( ISFILE_BIASCOR ) {
    if ( INPUTS.opt_biasCor == 0 ) 
      { INPUTS.opt_biasCor = MASK_BIASCOR_DEFAULT; }

    if ( IGNOREFILE(INPUTS.fieldGroup_biasCor)==0 )
      {  INPUTS.use_fieldGroup_biasCor = 1; }

    if ( IGNOREFILE(INPUTS.surveyGroup_biasCor)==0 )
      {  INPUTS.use_surveyGroup_biasCor = 1; }

    if ( INPUTS.use_fieldGroup_biasCor || INPUTS.use_surveyGroup_biasCor ) 
      { INPUTS.opt_biasCor |= MASK_BIASCOR_SAMPLE ;  }
  }

  // check for incompatible biasCor args
  if ( INPUTS.opt_biasCor > 0 && ISFILE_BIASCOR==0 ) {
    sprintf(c1err,"opt_biascor=%d but simfile_biascor=%s .",
	    INPUTS.opt_biasCor, f_biasCor);
    sprintf(c2err,"Provide simfile_biascor, or set opt_biascor=0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  // check for 'same' arg in either simfile_bias or simfile_ccprior.
  if ( strcmp(f_biasCor,"same") == 0 ) 
    { sprintf(f_biasCor,"%s", f_CCprior); }

  if ( strcmp(f_CCprior,"same") == 0 ) 
    { sprintf(f_CCprior,"%s", f_biasCor); }

  // if CCprior is not set then make sure to fix scalePCC and sigint
  if ( simdata_ccprior.USE == 0 ) {
    INPUTS.ipar[IPAR_scalePCC] = 0 ; // do NOT float scalePCC
  }
  else {
    INPUTS.ipar[IPAR_COVINT_PARAM] = 0 ;  // do not float sigint
    INPUTS.parval[IPAR_COVINT_PARAM] = INPUTS.sigmB; 
  }

  // if scalePCC=0 and is fixed, turn off CC prior
  if( INPUTS.ipar[IPAR_scalePCC] == 0  && 
      fabs(INPUTS.parval[IPAR_scalePCC]) < 0.000001 )  { 
    simdata_ccprior.USE = 0 ; 
    sprintf(INPUTS.simFile_CCprior,"NONE");
  }



  // if there is no user-defined selection of fit params,
  // then float everything. Otherwise stick with user choices.
  if ( simdata_ccprior.USEH11 ) {
    int ipar, IPAR;
    INPUTS.ipar[IPAR_scalePCC] = 1 ; // Jun 18 2018
    if ( NPAR_H11_USER == 0 ) {
      // use all of the default H11 parameters
      for(ipar=0; ipar < NPAR_H11_TOT; ipar++ ) 
	{ INPUTS.ipar[IPAR_H11+ipar]=1; }
    }
    else {
      // user params from input file only; 
      // make sure un-used parameter values are zero
      for(ipar=0; ipar < NPAR_H11_TOT; ipar++ ) {
	IPAR = IPAR_H11 + ipar ;
	if ( INPUTS.ipar[IPAR] == 0 ) { INPUTS.parval[IPAR]=0.0; }	
      }
    }
  }


  // ----------------------------------------

  for (i=0; i < MXCOSPAR;++i)    {

    if ( INPUTS.parstep[i] == 0.0 )   { INPUTS.ipar[i] = 0; } 

    if (INPUTS.ipar[i]==0) 
      { strcpy(usage,"fixed "); }
    else 
      { strcpy(usage,"varied"); }
   
    printf("%s  %s  starting value=%8.3f  (step=%.3f)\n",
	   FITPARNAMES_DEFAULT[i], usage,
	   INPUTS.parval[i], INPUTS.parstep[i] );         
  }


  if (INPUTS.zmax <= INPUTS.zmin) {
    printf("Invalid z range zmin=%f zmax=%f \n",INPUTS.zmin,INPUTS.zmax);
    exit(3);
  }

  // prepare stuff related to gamma = HR step
  prep_gamma_input();


  // sanity checks on fitting for GAMMA0 (mag step across host logmass)
  if ( INPUTS.USE_GAMMA0 ) {
    int    FLOAT_LOGMASS = INPUTS.ipar[IPAR_LOGMASS_CEN] ;
    double LOGMASS       = INPUTS.parval[IPAR_LOGMASS_CEN] ;
    double TAU           = INPUTS.parval[IPAR_LOGMASS_TAU] ;

    if ( TAU < 0.001 ) {
      sprintf(c1err,"Invalid LOGMASS_TAU = %.4f for gammma0 fit", TAU);  
      sprintf(c2err,"logmasss_tau(p8) should be at least .01");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  INPUTS.COSPAR[0] = INPUTS.parval[IPAR_OL] ;
  INPUTS.COSPAR[1] = INPUTS.parval[IPAR_Ok] ;
  INPUTS.COSPAR[2] = INPUTS.parval[IPAR_w0] ;
  INPUTS.COSPAR[3] = INPUTS.parval[IPAR_wa] ;
  INPUTS.FLOAT_COSPAR=0;
  if ( INPUTS.ipar[IPAR_OL] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_Ok] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_w0] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_wa] ) { INPUTS.FLOAT_COSPAR=1; }

  parse_nmax(INPUTS.nmaxString);

  // check for fast lookup option if all cosmo params are fixed.
  // arg, not needed.  prep_cosmodl_lookup();

  printf("\n"); fflush(stdout);

  return ;

} // end of prep_input

// **********************************************
void  prep_gamma_input(void) {

  // Created Oct 2018
  // prepare stuff related to gamma0 & gamma1 for HR step.

  int  LDMP = 1;
  char *varname_gamma = INPUTS.varname_gamma ;
  char fnam[] = "prep_gamma_input" ;

  // ---------- BEGIN -------------

  // set USE_GAMMA0 flag if 1) gamma is floated, or 2) initial gamma !=0

  if ( INPUTS.ipar[IPAR_GAMMA0] ) 
    { INPUTS.USE_GAMMA0 = 1; }
  if ( fabs(INPUTS.parval[IPAR_GAMMA0]) > 1.0E-8 ) 
    { INPUTS.USE_GAMMA0 = 1; }

  if ( INPUTS.USE_GAMMA0 == 0 ) { return ; }

  if ( LDMP ) {  printf("\n %s: set USE_GAMMA0 flag \n", fnam); }

  // check if varname_gamma is already on CUTWIN list;
  // if not, then add it to ensure that varname_gamma is read
  // from data file.

  int icut, FOUND=0;
  int NCUTWIN = INPUTS.NCUTWIN ;
  char *tmpName;
  for(icut=0; icut < NCUTWIN; icut++ ) {
    tmpName = INPUTS.CUTWIN_NAME[icut] ;
    if ( strcmp(varname_gamma,tmpName) == 0 ) { FOUND=1; }
  }

  if ( LDMP ) {
    if ( FOUND ) {
      printf("\t %s already on CUTWIN list. \n", varname_gamma );
    }
    else {
      printf("\t Append %s to CUTWIN list. \n", varname_gamma );
    }
  }

  if ( FOUND == 0 ) {
    INPUTS.NCUTWIN++ ;
    icut = INPUTS.NCUTWIN - 1;
    sprintf(INPUTS.CUTWIN_NAME[icut],"%s", varname_gamma);       
    INPUTS.CUTWIN_RANGE[icut][0] = -9.0E12 ;
    INPUTS.CUTWIN_RANGE[icut][1] = +9.0E12 ;
  }

  // load name of fitpar
  sprintf(FITPARNAMES_DEFAULT[IPAR_LOGMASS_CEN],"%s_cen", varname_gamma);
  sprintf(FITPARNAMES_DEFAULT[IPAR_LOGMASS_TAU],"%s_tau", varname_gamma);

  //  debugexit(fnam); // xxxxxx

  return ;

} // end prep_gamma_input

// **********************************************
void  prep_cosmodl_lookup(void) {

  // Nov 22 2017
  // If all cosmo params are fixed, then make binned
  // dl vs, z lookup for faster computation.

  int NBZ, iz ;
  double ZMIN, ZMAX, ZBIN, z, di, dl ;
  char fnam[] = "prep_cosmodl_lookup" ;

  // ------------ BEGIN --------------

  COSMODL_LOOKUP.USE = 0 ;
  if ( INPUTS.FLOAT_COSPAR ) { return ; }

  ZMIN = INPUTS.zmin - 0.002 ;
  ZMAX = INPUTS.zmax + 0.100 ;

  NBZ  = (int)( 1000.0*(ZMAX - ZMIN) ) ;
  ZBIN = (ZMAX - ZMIN)/ (double)NBZ ;

  printf(" Create cosmdl-vs-z lookup: %d z-bins from %.4f to %.4f \n",
	 NBZ, ZMIN, ZMAX);

  int MEMD = NBZ * sizeof(double) ;
  COSMODL_LOOKUP.NBZ  = NBZ;  
  COSMODL_LOOKUP.ZMIN = ZMIN ;
  COSMODL_LOOKUP.ZMAX = ZMAX ;
  COSMODL_LOOKUP.ZBIN = ZBIN ;
  COSMODL_LOOKUP.z   = (double*) malloc(MEMD);
  COSMODL_LOOKUP.dl  = (double*) malloc(MEMD);

  for(iz=0; iz < NBZ; iz++ ) {
    di = (double)iz + 0.5 ;
    z  = ZMIN + (ZBIN * di) ;
    dl = cosmodl_forFit(z,INPUTS.COSPAR);
    COSMODL_LOOKUP.z[iz]  = z;
    COSMODL_LOOKUP.dl[iz] = dl;
  }

  COSMODL_LOOKUP.USE = 1 ;
  //  debugexit(fnam);
  return ;

} // end prep_cosmodl_lookup

// *******************************
void  printmsg_fitStart(FILE *fp) {

  char fnam[] = "printmsg_fitStart" ;

  if ( INPUTS.NSPLITRAN > MXSPLITRAN ) {
    sprintf(c1err,"NSPLITRAN = %d exceeds bound.", INPUTS.NSPLITRAN);
    sprintf(c2err,"Check MXSPLITRAN = %d", MXSPLITRAN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  fprintf(fp,"\n %s\n\n", dotDashLine );

  if ( INPUTS.NSPLITRAN == 1 ) {
    fprintf(fp,"\t BEGIN FIT \n" );
  }
  else {
    fprintf(fp,"\t BEGIN FIT with Random subset %d of %d\n",
	    NJOB_SPLITRAN, INPUTS.NSPLITRAN);
  }

  fprintf(fp,"\n %s\n\n", dotDashLine );
  fflush(stdout);

} // end of printmsg_fitStart


// =================================================
void get_COVINT_model(int IDSAMPLE, double (*COVINT)[NLCPAR] ) {

  // Created Jan 26 2018
  // Return intrinsic scatter matrix based on user model.
  //
  // If IDSAMPLE < 0 --> flag to set COV00 = sigMB**2 
  //
  // Oct 9 2018: 
  //  check USE_SIGINT_FIX option, and pass IDSAMPLE as input arg.


  int USE_SIGINT_FIX    = ( strlen(INPUTS.sigint_fix) > 0 );
  int l,m;
  double COVINT_PARAM  = FITINP.COVINT_PARAM_FIX; 
  double SQSIGINT ;

  // ----------- BEGIN -------------

  for(l=0; l < NLCPAR; l++ ) {
    for(m=0; m < NLCPAR; m++ ) {
      COVINT[l][m] = FITINP.ZPOLY_COVMAT[l][m] ;
    }
  }


  // add conventional sigma_int^2, which is varied
  // with each fit iteration.

  if ( IDSAMPLE < 0 ) {
    // just need something approxiate
    SQSIGINT = INPUTS.sigmB * INPUTS.sigmB;
  }
  else if ( USE_SIGINT_FIX ) {
    // user inputs * iteratve scale to have chi2/dof=1
    SQSIGINT = SNDATA_INFO.sqsigint_fix[IDSAMPLE] * COVINT_PARAM ;
  }
  else {
    // interative sigint to have chi2/dof=1
    SQSIGINT    = (COVINT_PARAM * COVINT_PARAM ) ;
  }

  COVINT[0][0] += SQSIGINT ;

  return ;

}  // end get_COVINT_model


// ******************************************
void recalc_dataCov(){

  // Jan 26 2018 - Major refactor 
  //
  // Called between fit iterations,
  // Recompute full data covariance matrix (covtot) for each
  // event.  Covfit (from SALT2 LC fit) is awlays the same, 
  // but covint changes in the recalc.
  //
  // SCALE_COVINT:
  // * For traditional sigma_int, sigma_int is varied while
  //   SCALE_COVINT is fixed at 1.0.
  // 
  // * For COVINT_BIASCOR approach, COVINT_BIASCOR depends
  //   slightly on alpha & beta, and SCALE_COVINT is varied 
  //   to get chi2/dof ~ 1.0.
  //

  int n, l, m, idsample ;
  double z, covtot, covint, covfit, COVMAT_INT[3][3];
  double SCALE_COVINT ;

  int USE_IDEAL_COVINT  = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  int USE_SIGINT_FIX    = ( strlen(INPUTS.sigint_fix) > 0 );
  double a = FITRESULT.ALPHA ;
  double b = FITRESULT.BETA ;
  char fnam[] = "recalc_dataCov";

  // ------------ BEGIN --------------

  if ( USE_IDEAL_COVINT == 0  && USE_SIGINT_FIX==0 ) { 
    SCALE_COVINT = 1.0;     // traditional sigma_int model
  }
  else {    
    // Covint from biasCor or from sigint_fix(IDAMPLE)
    SCALE_COVINT = FITINP.COVINT_PARAM_FIX; 
  }

  // - - - - - - - - - - 
  for (n=0; n<FITINP.NSNCUTS; ++n)  {

    if ( data[n].skipfit ) { continue ; }
    idsample = data[n].idsample ;
    z        = data[n].zhd ;

    // fetch COVINT depending on method
    if ( USE_IDEAL_COVINT ) {
      get_COVINT_biasCor(idsample,z,a,b,  COVMAT_INT); 
    }
    else {
      get_COVINT_model(idsample,COVMAT_INT); 
    }

    //Add intrinsic scatter to error matrix
    for (l=0; l<3; ++l) {
      for (m=0; m<3; ++m) {
	covfit = data[n].covmat_fit[l][m];
	covint = COVMAT_INT[l][m] ;
	covtot = covfit + covint*SCALE_COVINT ;
	data[n].covmat_tot[l][m] = covtot ;
      }
    }
  }
}  // end recalc_dataCov


// *******************************************************
double next_covFitPar(double redchi2, double parval_orig, double parval_step) {

  // Created Jan 26 2018
  // Returns next covFitPar; either sigmB or covScale.

  double parval_next;
  double step, slope, num, denom;
  int NFIT_ITER = FITRESULT.NFIT_ITER ;


  char fnam[] = "next_covFitPar" ;

  // ------------- BEGIN ------------


  parval_next = parval_orig; // init

  FITRESULT.CHI2RED_LIST[NFIT_ITER] = redchi2 ;
  FITRESULT.SIGINT_LIST[NFIT_ITER]  = parval_orig ;

  if ( NFIT_ITER == 0 ){
    // decide if sigint needs to be larger or smaller
    if (redchi2 > 1.0) 
      { step = +parval_step ; } 
    else 
      { step = -parval_step ; }
    // calculate new sigint
    parval_next = parval_orig + step;

  } else { 

    num   = 
      FITRESULT.CHI2RED_LIST[NFIT_ITER] - 
      FITRESULT.CHI2RED_LIST[NFIT_ITER-1] ;

    denom = 
      FITRESULT.SIGINT_LIST[NFIT_ITER] - 
      FITRESULT.SIGINT_LIST[NFIT_ITER-1] ;

    slope = num/denom ;
    parval_next = parval_orig - (redchi2-1.0)/slope ;
  }

  return(parval_next);

} // end next_covFitPar



// ******************************************
void conflict_check() {

  char var1[200], var2[200];
  int exitnow;

  exitnow = 0;

  if ( INPUTS.zpolyflag == 1 && INPUTS.fitflag_sigmb > 0 ) {
    sprintf(var1, "zpolyflag");
    sprintf(var2, "fitflag_sigmb");
    exitnow = exitnow+1;
  }

  /*
  if (INPUTS.NSPLITRAN > 1 && INPUTS.fitflag_sigmb > 0 ) {
    sprintf(var1, "NSPLITRAN");
    sprintf(var2, "fitflag_sigmb");
    exitnow = exitnow+1;
  }
  */
  
  if (exitnow > 0){
    printf("\n FATAL ERROR:  \n");
    printf("\t Identified %d input variable conflicts. \n",exitnow);
    printf("\t Last conflict is between %s and %s. \n", var1, var2);
    printf("\t See SALT2mu.c subroutine conflict_check for more info.\n");
    printf("\t Check inputs  and try again. \n");
    printf("\t ***** ABORT ***** \n ");
    exit(1);
  }

  return ;
} // end conflict_check


// ******************************************
void outFile_driver(void) {

  // Created Nov 30 2017
  // [move code out of main]

  char *prefix = INPUTS.PREFIX ;
  char tmpFile1[200], tmpFile2[200];

  // --------------- BEGIN -------------

  if ( strlen(prefix) > 0 && IGNOREFILE(prefix) == 0 ) {

    if ( INPUTS.NSPLITRAN == 1 )  { 
      sprintf(tmpFile1,"%s.fitres", prefix ); 
      sprintf(tmpFile2,"%s.M0DIF",  prefix ); 
    }
    else  { 
      sprintf(tmpFile1,"%s-SPLIT%3.3d.fitres", prefix, NJOB_SPLITRAN);
      sprintf(tmpFile2,"%s-SPLIT%3.3d.M0DIF",  prefix, NJOB_SPLITRAN);
    }
    
    write_fitres(tmpFile1);  // write result for each SN
    write_M0(tmpFile2);      // write M0 vs. redshift
  } 
  else {
    printf("\n PREFIX not specified --> no fitres output.\n");
    fflush(stdout);
  }

  return ;

} // end outFile_driver

// ******************************************
void  write_M0(char *fileName) {

  // write M0 vs. z to fitres-formatted file.
  //
  //  VARNAME:  ROW  z  M0 M0err
  //
  // Jun 21 2017: 
  //  + NVAR->6 (not 4) to fix aweful bug
  // Jun 26 2017:
  //  + use get_M0_data() to get zM0 values.
  //
  // Jun 27 2017: REFACTOR z bins
  // Mar 26 2018: write MUREF column (NVAR->8)
  //

  int iz, izuse, NVAR, irow, NFIT ;
  double z, zMIN, zMAX, VAL, ERR, dl, MUREF;
  FILE *fp;
  char fnam[] = "write_M0" ;

  // ---------- BEGIN -----------

  if ( INPUTS.errmask_write == -9 ) { return ; } // July 2016

  //  if ( simdata_ccprior.USE == 0 ) { return ; }

  calc_zM0_data(); // fill FITRESULT.zM0[iz]

  fp = fopen(fileName,"wt");

  if ( !fp )  {
    sprintf(c1err,"Could not open M0-outFile");
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n Write MUDIF vs. redshift to %s\n" , fileName); 

  double OL = FITRESULT.PARVAL[NJOB_SPLITRAN][9] ;
  double w0 = FITRESULT.PARVAL[NJOB_SPLITRAN][11] ;
  fprintf(fp,"# Reference cosmology params for MUDIF: \n");
  fprintf(fp,"#   Omega_DE(ref):  %.3f \n", OL );
  fprintf(fp,"#   w_DE(ref):      %.3f \n", w0 );

  fprintf(fp,"#   alpha(fit):   %.4f +- %0.4f \n",
	  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_ALPHA0],
	  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_ALPHA0]  );
  fprintf(fp,"#   beta(fit):    %.4f +- %0.4f \n",
	  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_BETA0],
	  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_BETA0]  );

  fprintf(fp,"# \n");

  // write blindFlag info (Mar 1 2017)
  if ( INPUTS.blindFlag > 0 && ISDATA ) { write_blindFlag_message(fp); }

  write_MUERR_INCLUDE(fp);

  NVAR=8 ;
#ifdef TEXTFILE_NVAR
  fprintf(fp,"NVAR:  %d \n", NVAR);
#endif
  fprintf(fp,"VARNAMES: ROW  zMIN     zMAX     z        "
	  "MUDIF  MUDIFERR   MUREF  NFIT \n");
  irow = 0 ;

  for(iz=0; iz < INPUTS.nzbin; iz++ ) {

    irow++ ;    
    z   = INPUTS.BININFO_z.avg[iz] ;

    if ( INPUTS.opt_biasCor > 0 ) { z = FITRESULT.zM0[iz]; }

    zMIN   = INPUTS.BININFO_z.lo[iz] ;
    zMAX   = INPUTS.BININFO_z.hi[iz] ;
	
    VAL   = FITRESULT.M0DIF[iz];
    ERR   = FITRESULT.M0ERR[iz];

    dl    = cosmodl_forFit(z, INPUTS.COSPAR) ;
    MUREF = 5.0*log10(dl) + 25.0 ;
    //    MUREF = 33.33 ;
    NFIT = FITINP.NZBIN_FIT[iz] ;

    if ( VAL == 0.0 && ERR == 0.0 )
      { ERR = 999.0 ; }
    if ( FITINP.ISFLOAT_z[iz] == 0 )
      { VAL=0.0 ; ERR=999.0;  }
    
    fprintf(fp, "ROW:     "
	    "%2d  %7.4f %7.4f %7.4f  "
	    "%9.4f %9.4f  %.4f %4d\n",
	    irow, zMIN, zMAX, z, 
	    VAL, ERR, MUREF, NFIT );
    fflush(fp);
  }

  fclose(fp);

  return ;

} // end write_M0



// ******************************************
void write_fitres(char* filnam) {

  // Write outout in fitres format.
  // Combine original fitres variables with those
  // calculated here in SALT2mu.
  // Jun 2012: increase line-string lenght from 400 to MXCHAR_LINE
  // Jul 2013: write averge M0 to fitres file (avemag0)
  // Apr 2016: for SIM, write true CC/tot ratio of fitted events.
  // Jul 2016: if errmask_write=-9, do NOT write the SN lines
  // Jun 27 2017: REFACTOR z bins
  // Mar 01 2018: add M0 to output

  double VAL, ERR, PULL ;
  FILE  *fout, *fin;

  int n, ivar, indx, NCUT, icut, errcode, NWR, ISFLOAT, iz, GZIPFLAG ;
  int idsample ;
  char 
     line[MXCHAR_LINE]
    ,tmpLine[MXCHAR_LINE] 
    ,tmpName[60], ztxt[60]
    ,*ptrtok, KEY[40], CCID[40], *ptrCR
    ,mode_rt[] = "rt"
    ;

  //  char fnam[] = "write_fitres" ;
  // ------------------ BEGIN ----------------


  // first define the new fitres variables to add to the
  // original list  
  NVAR_NEW = 0;
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MU");
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MUMODEL");  // June 2016
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MUERR");
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MUERR_RAW"); // May 2016
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MURES");
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"MUPULL");
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"M0DIF");     // Mar 1 2018
  NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"ERRCODE");

  if  ( FOUNDKEY_SIM   ) {
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"SIM_MU"); // for back-compat.
    if ( FOUNDKEY_SIMx0 ) 
      { NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"SIM_MB"); }
  }


  if ( simdata_bias.NROW > 0 ) {
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasCor_mu");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasCorErr_mu");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasCor_mB");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasCor_x1");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasCor_c");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"biasScale_muCOV");
    NVAR_NEW++ ;  sprintf(VARNAMES_NEW[NVAR_NEW],"IDSAMPLE");
  }

  // - - - - - - - - - -
  NVAR_TOT = NVAR_ORIG + NVAR_NEW ;

  fout = fopen(filnam,"wt");
  fin  = open_TEXTgz(INPUTS.filename, "rt", &GZIPFLAG); // check for .gz file

  printf("\n Open output file with  errmask_write=%d : \n", 
	 INPUTS.errmask_write );
  printf("\t %s \n", filnam );


  fprintf(fout,"# MU-RESIDUAL NOTE: MURES = MU-(MUMODEL+M0DIF) \n\n");

  write_MUERR_INCLUDE(fout);

  if ( INPUTS.errmask_write != -9 ) { 
    // write standard header keys
#ifdef TEXTFILE_NVAR
    fprintf(fout,"NVAR: %i \n", NVAR_TOT ); 
#endif
    fprintf(fout,"VARNAMES: ");
    
    for ( ivar=1; ivar <= NVAR_ORIG; ivar++ ) {      
      fprintf(fout,"%s ", VARNAMES_ORIG[ivar-1] );
    }
    for ( ivar=1; ivar <= NVAR_NEW; ivar++ )   {  
      fprintf(fout,"%s ", VARNAMES_NEW[ivar] );  
    }
  }
  else {
    // no SN lines, hence no header keys
  }


  NCUT = INPUTS.NCUTWIN ;
  if ( NCUT > 0 ) {
    fprintf(fout,"\n\n# CUTWIN Selection: \n");
    for ( icut=0; icut<NCUT; icut++ ) {
      fprintf(fout, "#\t %10.4f <= %12s <= %10.4f \n"
	     ,INPUTS.CUTWIN_RANGE[icut][0]
	     ,INPUTS.CUTWIN_NAME[icut]
	     ,INPUTS.CUTWIN_RANGE[icut][1]
	     ) ;
    }
  }

  if ( INPUTS.NSPLITRAN > 1 ) {
    fprintf(fout,"#\t SPLIT-RAN job %d of %d \n",
	    NJOB_SPLITRAN, INPUTS.NSPLITRAN  );
  }


  fprintf(fout,"\n\n# SALT2mu fit results: \n");

  write_fitres_misc(fout);

  // ---------------- now write fitted params ----------------


  if ( strlen(INPUTS.sigint_fix) == 0 ) {
    // fitted sigint
    fprintf(fout,"#  %s   = %0.4f  (%i iterations)\n",
	    FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	    FITINP.COVINT_PARAM_FIX, FITRESULT.NFIT_ITER );    
  }
  else {
    // print user specified sigint for each IDSAMPLE
    char *NAME;
    double sigint;
    for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) {
      NAME   = SAMPLE_BIASCOR[idsample].NAME ;
      sigint = SNDATA_INFO.sigint_fix[idsample]* FITINP.COVINT_PARAM_FIX; 
      fprintf(fout,"#  sigint(%s) = %6.4f   [sigint(user) * %.3f]\n",
	      NAME, sigint, FITINP.COVINT_PARAM_FIX );
    }
    
  }


  

  double chi2sum_m0 = 0.0 ;
  int    NBIN_m0 = 0 ;
  int    ISM0;

  for ( n=0; n < FITINP.NFITPAR_ALL ; n++ ) {

    ISFLOAT = FITINP.ISFLOAT[n] ;
    ISM0    = n >= MXCOSPAR ; // it's z-binned M0

    // skip fixed cosmo params; keep fixed M0 params to 
    // print clear message about unused z-bins
    if ( ISM0 == 0  && ISFLOAT==0 ) { continue ; } 
      
    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;
    sprintf(tmpName,"%s", FITRESULT.PARNAME[n]);
    ztxt[0] = 0 ;
    if ( ERR < 0.0 ) { continue ; }

    if ( n >= MXCOSPAR ) { 
      iz    = INPUTS.izpar[n] ;
      // one of the m0 params in a z bin
      VAL = FITRESULT.M0DIF[iz]; 
      sprintf(tmpName,"%s-<M0avg>", FITRESULT.PARNAME[n] );

      sprintf(ztxt,"(%5.3f < z < %5.3f, N=%d)", 
	      INPUTS.BININFO_z.lo[iz], INPUTS.BININFO_z.hi[iz],
	      FITINP.NZBIN_FIT[iz] ) ;
    }
    else {
      VAL += BLIND_OFFSET(n); // offset cosmo params besides M0
    }

    if ( ISFLOAT ) {
      fprintf(fout,"#  %-14s = %10.5f +- %8.5f  %s \n",
	      tmpName, VAL, ERR, ztxt );
    }
    else {
      fprintf(fout,"#  %-14s  too few events->IGNORED  %s \n",
	      tmpName, ztxt );
    }

    if ( ISM0 && ERR > 0.0 && ISFLOAT ) 
      { NBIN_m0++ ; PULL = VAL/ERR ;  chi2sum_m0 += (PULL*PULL); }

  } // end loop over fit params

  // --------------
  fprintf(fout,"#  %-14s = %12.5f \n", "M0avg", 
	  FITRESULT.AVEMAG0 ); // Jul 8 2013

  fprintf(fout,"#  m0-M0avg chi2/dof = %.1f / %d \n",
	  chi2sum_m0, NBIN_m0-1);

  // ----------
  if ( INPUTS.blindFlag > 0 && ISDATA ) { write_blindFlag_message(fout); }

  if ( INPUTS.uave)
    { FITRESULT.SNMAG0 = FITRESULT.AVEMAG0; }
  else 
    { FITRESULT.SNMAG0 = INPUTS.nommag0; }


  fprintf(fout," \n");
  
  indx = 0;
  NWR  = 0;
  fflush(fout);

  // check option to NOT write each SN to have smaller file
  // with only the fit results
  if ( INPUTS.errmask_write == -9 ) { 
    fclose(fout);
    printf(" Nothing written out  (%d/%d used in fit) \n", 
	   FITRESULT.NSNFIT , FITINP.NSNTOT );  fflush(stdout);    
    return ;
  }

  int NDUM=0;
  while ( fgets (line, MXCHAR_LINE, fin) !=NULL  ) {

    if ( line[0] == ' '   ) { continue ; }
    if ( strlen(line) < 3 ) { continue ; }

    // break tmpline into blank-separated strings
    sprintf(tmpLine,"%s", line);
    ptrtok = strtok(tmpLine," ");
    sscanf ( ptrtok, "%s", KEY );
    ptrtok = strtok(NULL, " ");
    sscanf ( ptrtok, "%s", CCID );
    //    ptrtok = strtok(NULL, " ");

    //    NDUM++ ; if(NDUM<30) { printf(" xxx read line='%s' \n", line); }

    if ( strcmp(KEY,"SN:") != 0 ) { continue ; }

    // remove <CR> from end of line (requested by Rahul)
    ptrCR = strchr(line,'\n');
    if (ptrCR) {*ptrCR = ' ';}

    // get index for data[] array
    get_CCIDindx(CCID, &indx) ;

    if ( data[indx].skipfit ) { continue ; } // May 2016

    // check which error codes to allow (RK, May 2012)
    errcode = data[indx].errmask ;
    if ( keep_errcode(errcode) == 1 ) { 

      // Print SN line from input fiters file
      fprintf(fout, "%s ", line);

      // append variables computed by SALT2mu
      append_fitres(fout, CCID, indx);      
      fprintf(fout,"\n");  fflush(fout);
      NWR++ ;
    }

  }  // fgets

  fclose(fout);
  printf(" Wrote %d SN  (%d/%d used in fit) \n", 
	 NWR, FITRESULT.NSNFIT , FITINP.NSNTOT );
  fflush(stdout);

  return ;

} // end of write_fitres


// ===============================================
void write_MUERR_INCLUDE(FILE *fp) {

  // Created Jun 25 2017
  // write vpec and lensing info in output file 
  // so that cosmology fitting program knows if/what
  // is already included in MUERR. Should help avoid
  // double-counting MUERR contributions.

  int USE=0;
  double tmpErr;
  char fnam[] = "write_MUERR_INCLUDE";

  // --------------- BEGIN ----------------

  fprintf(fp,"\n");

  fprintf(fp, "MUERR_INCLUDE: zERR \n" );  USE=1;

  tmpErr = INPUTS.zpecerr ;
  if ( tmpErr > 0.0 ) {
    fprintf(fp, "MUERR_INCLUDE: zPECERR=%.4f \n", tmpErr );
    USE=1;
  }

  tmpErr = INPUTS.lensing_zpar;
  if ( tmpErr > 0.0 ) {
    fprintf(fp, "MUERR_INCLUDE: SIGMA_LENS=%4f*z \n", tmpErr );
    USE=1;
  }

  if ( USE )  { fprintf(fp,"\n"); fflush(fp); }

  return ;

} // end write_MUERR_INCLUDE

// ===============================================
void write_blindFlag_message(FILE *fout) {
  int  blindFlag = INPUTS.blindFlag ;
  int  ipar;
  double *blindpar, parval_orig ;
  // ----------- BEGIN --------------
  fprintf(fout,"#\n");
  fprintf(fout,"#  *** WARNING: RESULTS ARE BLINDED *** \n");

  if ( (blindFlag & BLINDMASK_MUz)>0 ) {
    fprintf(fout,"#  m0_## values have blind offset = cos( z * %.3f ) \n",
	    INPUTS.blindPar ); 
  }
  else {
    for(ipar=0; ipar < MAXPAR; ipar++ ) {
      blindpar = INPUTS.blindpar[ipar];
      if ( blindpar[0] != 0.0 ) {
	parval_orig = INPUTS.parval[ipar] - blindpar[0]*cos(blindpar[1]);
	fprintf(fout,"#     %s = %8.4f + %5.2f*cos(%f) \n",
		FITRESULT.PARNAME[ipar], parval_orig,
		blindpar[0], blindpar[1] ); fflush(stdout);
      }
    }
  }

  fprintf(fout,"#\n");
  return ;

} // end write_blindFlag_message


// ===============================================
void write_fitres_misc(FILE *fout) {

  // write misc info to fitres file depending on whether
  // CC prior is used.
  // Apr 2016: write CC/tot ratio for simulation

  double chi2min, chi2red ;
  int    NDOF, NSN ;
  char   textCC[100];

  // ------------- BEGIN -----------

  if ( simdata_ccprior.USE == 0 ) {
    chi2min = FITRESULT.CHI2SUM_1A; 
    chi2red = FITRESULT.CHI2RED_1A; 
    NDOF    = FITRESULT.NDOF;
    NSN     = FITRESULT.NSNFIT ;
    textCC[0] = 0 ;
  }
  else {
    chi2min = FITRESULT.CHI2SUM_1A ;
    chi2red = FITRESULT.CHI2RED_1A ; 
    NDOF    = (int)FITRESULT.NSNFIT_1A - FITINP.NFITPAR_FLOAT ;
    NSN     = FITRESULT.NSNFIT ;
    sprintf(textCC,"(for ProbCC/ProbSum < %.3f)", 
	    INPUTS.maxProbCC_for_sigint );
  }

  char string_truecc[60] = "" ;
  double frac ;
  if ( FOUNDKEY_SIM ) {
    frac = (double)FITRESULT.NSNFIT_TRUECC/(double)FITRESULT.NSNFIT ;
    sprintf(string_truecc,"[ true CC/(Ia+CC) = %.4f) ]", frac);
  }

  int MASK  = INPUTS.opt_biasCor ;
  int NUMD = 5 ;    // default number of biasCor dimensions
  if ( MASK > 0 ) {
    if ( MASK & MASK_BIASCOR_1DZ ) { NUMD=1; } 
    fprintf(fout,"#  NSIM(%dD-BIASCOR) = %d   "
	    "(N_alpha x N_beta = %d x %d) \n"
	    ,NUMD
	    ,simdata_bias.NUSE
	    ,simdata_bias.BININFO_SIM_ALPHA.nbin
	    ,simdata_bias.BININFO_SIM_BETA.nbin
	    );

    fprintf(fout,"#  sigint(%dD-BIASCOR) = %.3f \n",
	    NUMD, simdata_bias.SIGINT_AVG );
     
  }
    
  fprintf(fout,"#  NSNFIT   = %d    %s \n", NSN, string_truecc );
    
  fprintf(fout,"#  -2log(L) = %.2f \n", FITRESULT.CHI2SUM_MIN );

  fprintf(fout,"#  chi2(Ia)/dof = %.2f/%i = %.3f  \n",
	  chi2min, NDOF, chi2red );
      

  fflush(fout);
  return ;

} // end write_fitres_misc

// =====================================
int ISBLIND_FIXPAR(int ipar) {

  // Created Aug 2017
  // Returns 1 if this 'ipar' parameter is blinded
  
  int blindFlag = INPUTS.blindFlag ;
  int NOTBLIND=0, BLIND=1;
  int DATAFLAG;
  // --------------- BEGIN --------------

  // first check if blindFlag is set for FIXPAR option
  if ( (blindFlag & BLINDMASK_FIXPAR)==0 ) 
    { return(NOTBLIND); } 


  // continue for DATA or if special flag is set to blind sim.
  DATAFLAG = ( ISDATA || ( blindFlag & BLINDMASK_SIM)>0 );
  if ( DATAFLAG == 0 ) { return(NOTBLIND); }

  if ( INPUTS.blindpar[ipar][0] != 0.0 ) 
    { return(BLIND); } 
  else 
    { return(NOTBLIND); } 

} // end ISBLIND_FIXPAR


// =====================================
double BLIND_OFFSET(int ipar) {

  // Created Feb 2016
  // Input ipar is the MINUIT M0off parameter 
  // Returns blind offset for M0Off parameter.
  // This is called after fitting to blind the fitted M0
  // before printing.

  double zero = 0.0 ;
  // ---------------------
  // never blind for sims or if blindflag is turned off
  if ( FOUNDKEY_SIM )           
    { return zero ; }

  // bail of MU-vs-z blinding is NOT set
  if ( (INPUTS.blindFlag & BLINDMASK_MUz)== 0 ) 
    { return zero ; }

  // do NOT blind for nuissance parameters.
  if ( ipar < MXCOSPAR )   { return zero ; }

  // blind only for the magOff params
  double z, arg, off ;
  int iz ;
  iz  = INPUTS.izpar[ipar];
  z   = INPUTS.BININFO_z.lo[iz] ;  
  INPUTS.blindPar = 10.0 ;
  arg = z * INPUTS.blindPar ;
  off = cos(arg) ;
  return(off);

} // end BLIND_OFFSET


// ==========================================
void get_CCIDindx(char *CCID, int *indx) {


  // May 23, 2012
  // return integer *indx corresponding to CCID

  int n;
  char fnam[] = "get_CCIDindx";

  // -------------- BEGIN -----------------

  n = *indx;
  while ( strcmp(data[n].name,CCID) != 0  ) {

    n++ ;  

    // abort if  we cannot find the indx
    if( n > FITINP.NSNCUTS ) {
      sprintf(c1err,"Cannot find SALT2mu index for CCID='%s'", CCID);
      sprintf(c2err,"after checking all %d SN", FITINP.NSNCUTS ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }
  *indx = n ; // load output arg

} // end of get_CCIDindx



void append_fitres(FILE *fp, char *CCID, int  indx ) {

  // Oct 17, 2011 R.Kessler
  // append line in fitres file with SALT2mu info
  // corresponding to this CID.
  // *indx is the current data[] index to quickly
  // lookup the next index for the current CID. 
  //
  // Mar 1 2018: add M0
  
  double 
    z, zerr, mu, muerr, muerr2, mumodel, mures, pull, M0DIF
    ,sim_mb, sim_mu, mb, mberr, s, serr, sbar,  c, cerr, cbar 
    ;

  int n, errcode, NWR ;

  char fnam[20] = "append_fitres" ;

  // ------------ BEGIN --------------

  n =  indx;

  z       = data[n].zhd;
  zerr    = data[n].zhderr;
  mumodel = data[n].mumodel;
  mu      = data[n].mu - FITRESULT.SNMAG0; 
  muerr   = data[n].muerr;
  muerr2  = data[n].muerr_raw;
  mures   = data[n].mures;
  pull    = data[n].pull;
  M0DIF   = data[n].M0 - FITRESULT.AVEMAG0 ;
  errcode = data[n].errmask  ;

  sim_mb  = data[n].sim_mb ;
  sim_mu  = data[n].sim_mu ;

  mb     = data[n].fitpar[0];    mberr  = sqrt(data[n].covmat_tot[0][0]);
  s      = data[n].fitpar[1];    serr   = sqrt(data[n].covmat_tot[1][1]);
  c      = data[n].fitpar[2];    cerr   = sqrt(data[n].covmat_tot[2][2]);
  
  if (pull > 99.999) { pull=99.999; }
  //Get fitted values for s and c
  //These are the best values that put the SN on the Hubble diagram
  fitsc(n,&sbar,&cbar);
  
  NWR=0;

  NWR++ ; fprintf(fp, "%7.4f ", mu    );
  NWR++ ; fprintf(fp, "%7.4f ", mumodel); // calculated distance from model
  NWR++ ; fprintf(fp, "%7.4f ", muerr ); 
  NWR++ ; fprintf(fp, "%7.4f ", muerr2); // no sigInit, no corrections
  NWR++ ; fprintf(fp, "%7.4f ", mures );
  NWR++ ; fprintf(fp, "%6.3f ", pull  );
  NWR++ ; fprintf(fp, "%7.4f ", M0DIF );
  NWR++ ; fprintf(fp, "%d ",    errcode );
  if ( FOUNDKEY_SIM ) {
    NWR++ ; fprintf(fp, "%7.4f ",  sim_mu ); // same as SIMDLMAG
    if ( FOUNDKEY_SIMx0 ) 
      { NWR++ ; fprintf(fp, "%7.4f ",  sim_mb ); }
  }


  if ( simdata_bias.NROW > 0 ) {
    NWR++ ; fprintf(fp, "%6.3f ", data[n].muBias );
    NWR++ ; fprintf(fp, "%6.3f ", data[n].muBiasErr );
    NWR++ ; fprintf(fp, "%6.3f ", data[n].fitParBias[INDEX_mB] );
    NWR++ ; fprintf(fp, "%6.3f ", data[n].fitParBias[INDEX_x1] );
    NWR++ ; fprintf(fp, "%6.3f ", data[n].fitParBias[INDEX_c] );
    NWR++ ; fprintf(fp, "%6.3f ", data[n].muCOVscale ) ;
    NWR++ ; fprintf(fp, "%d "   , data[n].idsample ) ;
  }


  fflush(fp);

  if ( NWR != NVAR_NEW ) {
    sprintf(c1err,"Expected to write NVAR_NEW=%d SALT2mu variables",
	    NVAR_NEW);
    sprintf(c2err,"but wrote only %d variables.", NWR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


} // end of append_fitres


// **************************************
double avemag0_calc(int opt_dump) {

  // Created Aug 22, 2012 
  // code moved from main and fixed bug.
  //
  // Jan 24, 2013: 
  //   fix aweful bug found by Y.Lo;  iz -> izuse to compute mag0.
  //   Fix affects fit when iz bins with too few events are skipped.
  //
  // Feb 2016: opt_dump>0 --> print result to stdout
  // Jun 27 2017: REFACTOR z bins

  int nzfit, iz, isplit, IOFF_MAG0 ;
  double norm, mag0, d_nzfit, ave ;
  
  // -------------- BEGIN -----------------

  norm    = 0.0;
  ave     = 0.0;
  isplit  = NJOB_SPLITRAN ;
  IOFF_MAG0 = MXCOSPAR ;

  for (iz=0; iz < INPUTS.nzbin; iz++ )     {

    if ( FITINP.ISFLOAT_z[iz] == 0 ) { continue ; }

    nzfit    = FITINP.NZBIN_FIT[iz] ;  // number of SN fitted in this z bin
    d_nzfit  = (double)nzfit ;
    mag0     = FITRESULT.PARVAL[isplit][IOFF_MAG0+iz] ;

    norm    += d_nzfit ;   
    ave     += d_nzfit * mag0 ;

    /*
    printf(" xxx iz=%d  nzbin=%d  norm=%f mag0=%f  \n", 
	   iz, nzfit , norm, mag0 );
    */
    
  }

  ave /= norm; // global
  
  if ( opt_dump ) {
    printf("Average mag0 offset = %f  (wgted by NSN per z bin)\n", ave);
    fflush(stdout);
  }

  return(ave);

}  // end of avemag0_calc


// *********************************************
void  M0dif_calc(void) {

  // compute M0dif[iz] = M0(iz) - AVEMAG0,
  // and apply optional blind offset.
  // Note that iz-index is sparse index from 0 to NZUSE-1
  //
  // Jun 27 2017: REFACTOR z bins

  int iz, n;
  double VAL, ERR;

  for ( n=MXCOSPAR ; n < FITINP.NFITPAR_ALL ; n++ ) {

    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;

    // one of the m0 params in a z bin
    VAL -= FITRESULT.AVEMAG0 ; 
    VAL += BLIND_OFFSET(n);

    iz  = INPUTS.izpar[n] ; 
    FITRESULT.M0DIF[iz] = VAL ;
    FITRESULT.M0ERR[iz] = ERR ;
  }

} // end M0dif_calc

// *********************************************
void  printCOVMAT(int NPAR) {

  // Created Nov 30 2017
  // [move code out of main]

  int num, i, j ;
  double emat[MAXPAR][MAXPAR], terr[MAXPAR], corr ;

  // --------- BEGIN ----------
  printf("\nFull reduced COV matrix\n"); fflush(stdout);
  num = MAXPAR;
  mnemat_(emat,&num);
  for ( i=0 ; i < NPAR ; ++i )
    {
      terr[i] = sqrt(emat[i][i]);
      //      printf("%2i %.4e ",i+1,terr[i]); 
      printf("%2i ",i+1); 
      for (j=0;j<i;++j)	{
	corr = emat[i][j]/(terr[i]*terr[j]);
	printf("%6.3f ",corr);
      }
      printf("\n");  fflush(stdout);
    }

  return ;

} // end printCOVMAT

// *********************************************
int keep_errcode(int errcode) {

  // Returns 1 if this errcode is allowed and the
  // SN is written to the output fitres file.
  // Returns 0 if this error code is not allowed;
  // the SN is not written out.
  //
  // User input 'errmask_write' is a mask of
  // errocode bits that are allowed; if any other
  // bit is set then return 0.
  //

  int ovp;
  //  char fnam[] = "keep_errcode" ;

  // ---------------- BEGIN -------------

  // -1 is the same as allowing all bits.
  if (INPUTS.errmask_write == -1 ) { return 1 ; }

  ovp = ( INPUTS.errmask_write & errcode );
  
  if ( (errcode - ovp) > 0 ) 
    { return 0 ; }
  else
    { return 1 ; }

} // end of keep_errcode


// ==============================================
double cosmodl_forFit(double z, double *cosPar) {

  // Created Jan 2016.
  //
  // This is a shell function to avoid nan when Omega_L > 1 with
  // fixed curvature.  When OL > 0.97 do an extrapolation based
  // on the slope ddl/dOL near OL ~ 1.

  int ipar, i;
  int IPAR_OL = 0 ;
  double dl, DL[2], OL, slp, cosPar_local[10] ;
  double OL_extrap[2] = { 0.97, 0.99 } ;
  //  char fnam[] = "cosmodl_forFit" ;

  // -------------- BEGIN -------------

  OL = cosPar[0];

  if ( OL > OL_extrap[0] ) {
    for(ipar=0 ; ipar< NCOSPAR ; ipar++ ) 
      { cosPar_local[ipar] = cosPar[ipar] ; }

    for(i=0; i <2; i++ ) {
      cosPar_local[IPAR_OL] = OL_extrap[i] ;
      DL[i] = cosmodl(z,cosPar_local);
    }

    slp = (DL[1] - DL[0]) / (OL_extrap[1] - OL_extrap[0]) ;
    dl  = DL[0] + ( OL - OL_extrap[0] )*slp ;
  }
  else {
    dl = cosmodl(z,cosPar);
  }

  /*
  if ( OL > 2.8 ) {
    printf(" zzz z=%.4f  OL = %.4f  dl=%f (dl0,dl1=%f,%f) \n", 
	   z, OL, dl, dl0, dl1 ); fflush(stdout);
    fflush(stdout);
  }
  */


  return(dl);

} // end cosmodl_forFit

double cosmodl(double z, double *cosPar)
{
  const double  cvel = LIGHT_km; // 2.99792458e5;
  const double  tol  = 1.e-6;
  double dflat, distance, H0inv, dl ;
  double omega_m, omega_l, omega_k, OK, wde, wa ;
  double dz, zbin, dl0, dl1, z0, z1 ;
  int iz;
  char fnam[] = "cosmodl" ;

  // ------------- BEGIN --------------

  omega_l = cosPar[0]; 
  omega_k = cosPar[1];
  wde     = cosPar[2];
  wa      = cosPar[3];

  if(fabs(omega_k)<tol) { omega_k = 0.0; }
  omega_m = 1.0 - omega_l - omega_k;
  OK      = fabs(omega_k);

  //comoving distance to redshift

  dflat = rombint(inc, 0.0, z, cosPar, tol);

  H0inv = 1.0/INPUTS.H0 ;

  if( omega_k == 0.0  ) 
    { distance = cvel*H0inv*dflat; }
  else if(omega_k<0.0) 
    { distance = cvel*H0inv*(1.0/sqrt(OK))*sin(sqrt(OK)*dflat); }
  else 
    { distance = cvel*H0inv*(1.0/sqrt(OK))*sinh(sqrt(OK)*dflat); }

  return( (1.0+z)*distance );

} // end cosmodl



double rombint(double f(double z, double *cosPar),
	       double a, double b, double *cosPar, double tol) {

  /*
    Rombint returns the integral from a to b of using Romberg integration.
                                                                               
    The method converges provided that f(x) is continuous in (a,b).
    f must be double precision and must be declared external in the 
    calling routine.  Input "tol" indicates the desired relative 
    accuracy in the integral.
  */
  
#define MAXITER 23
#define MAXJ 5
  double g[MAXJ+2];
  double h,gmax,error,g0,fourj,g1, ztmp;
  int nint,i,k,j,jmax;

  g0 = 0.0 ; // init outut value (Jan 2017)
  h = 0.5*(b-a);
  gmax = h*( f(a,cosPar) + f(b,cosPar) );
  g[0] = gmax;
  nint = 1;
  error = 1.0e20;
  for (i=1;i<=MAXITER;++i) {
    if (fabs(error)<tol) break;
    //Calculate next trapezoidal rule approximation to integral.
    g0=0.0;
    for (k=0;k<nint;++k) { 
      ztmp = a+(2*k+1)*h ;
      g0 += f(ztmp,cosPar) ; 
    }
    g0 = 0.5*g[0] + h*g0;
    h  = 0.5*h;
    nint = nint + nint;
    if (i<MAXJ) { jmax=i; }
    else        { jmax=MAXJ; }
    fourj=1.0;
    
    for (j=0;j<jmax;++j)      {
      //  ! Use Richardson extrapolation.
      fourj=4.0*fourj;
      g1=g0+(g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }
    
    if (fabs(g0)>tol) 
      { error = 1.0 - gmax/g0; }
    else 
      { error = gmax; }
    
    gmax     = g0;
    g[jmax]  = g0;
  }
  
  if (i==MAXITER && fabs(error)>tol) {
    printf("Rombint failed to converge; integral = %e +- %e \n",
	   g0, error);
  }
  
  return(g0);
} // end rombint


double inc(double z, double *cosPar) {

  // Nov 22 2017: speed things up a little with fabs checks on zzpow and wa

  double hubble, rhode, omega_m, omega_k, omega_l, wde, wa ;
  double zz, zzpow;
  
  zz = 1.0 + z;
  
  omega_l = cosPar[0]; 
  omega_k = cosPar[1];
  wde     = cosPar[2];
  wa      = cosPar[3];

  if ( fabs(omega_k) < 1.0E-6 ) { omega_k = 0.0 ; }
  omega_m = 1.0 - omega_l - omega_k;

  zzpow = 3.0 * ( 1.0 + wde + wa ) ;

  rhode = omega_l; 
  if ( fabs(zzpow) > 1.0E-9 ) { rhode *= pow(zz,zzpow); }
  if ( fabs(wa)    > 1.0E-6 ) { rhode *= exp(-3.0*(wa*z/zz)); }
  
  hubble = sqrt( (omega_m*(zz*zz*zz)) + rhode + (omega_k*(zz*zz)) );
  
  return(1.0/hubble);
					  
} // end inc

void gengauss(double r[2])
{
  double radius, phi;
  radius = sqrt(-2.0*log(rand()/RNORM));
  phi = dTWOPI*rand()/RNORM;
  //  printf("radius %f phi %f \n",radius,phi);
  r[0] = radius*cos(phi);
  r[1] = radius*sin(phi);
  return;
}
void ludcmp(double* a, const int n, const int ndim, int* indx, 
       double* d, int* icon)
{
  /* System generated locals */
  int i1, i2, i3;

  /* Local variables */
  int imax = 0;
  int  i, j, k;
  double aamax, dum, sum;
  double vv[100] = { 0.0 };
    
  /* Function Body */
  *icon = 1;
  *d = 1.0;
  i1 = i2 = -9 ; // init to avoid compiler warnings (Jan 2017)

  for (i = 0; i < n; ++i) 
  {
    aamax = 0.0;
  
    for (j = 0; j < n; ++j) 
    {
      if (fabs(a[i + j * ndim]) > aamax) 
      {
	aamax = fabs(a[i + j * ndim]);
      }
    }

    /* MH       IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.' */
    if (aamax == 0.0) 
    {
      printf("LU decomposition (ludcmp.c) : Singular matrix, but continue. %d %d %d %d %d\n",i,j,i1,i2,ndim);
      *icon = -10;
      return;
      
    }

    vv[i] = 1.0 / aamax;
  }

  /* OK to HERE */
  i1 = n;
  for (j = 0; j < n; ++j) 
  {
    if (j > 0) 
    {
      for (i = 0; i < j; ++i) 
      {
	sum = a[i + j * ndim];
	if (i > 0) 
	{
	  i3 = i - 1;
	  for (k = 0; k < i; ++k) 
	  {
	    sum -= a[i + k * ndim] * a[k + j * ndim];
	  }
	  a[i + j * ndim] = sum;
	}
      }
    }

    aamax = 0.0;

    for (i = j; i < n; ++i) 
    {
      sum = a[i + j * ndim];
      if (j > 0) 
      {
     
	for (k = 0; k < j; ++k) 
	{
	  sum -= a[i + k * ndim] * a[k + j * ndim];
	}
	a[i + j * ndim] = sum;
      }

      dum = vv[i] * fabs(sum);
      if (dum >= aamax) 
      {
	imax = i;
	aamax = dum;
      }
    }

    if (j != imax) 
    {
    
      for (k = 0; k < n; ++k) 
      {
	dum = a[imax + k * ndim];
	a[imax + k * ndim] = a[j + k * ndim];
	a[j + k * ndim] = dum;
      }
      *d = - *d;
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if (j != (n-1)) 
    {
      /*          IF(A(J,J).EQ.0.)A(J,J)=TINY */
      if (fabs(a[j + j * ndim]) <= 1.e-20) 
      {
	if (a[j + j * ndim] > 0.0) 
	{
	  a[j + j * ndim] = 1.0e-20;
	} 
	else 
	{
	  a[j + j * ndim] = -1.0e-20;
	}
      }

      dum = 1.0 / a[j + j * ndim];
     
      for (i = j + 1; i < n; ++i) 
      {
	a[i + j * ndim] *= dum;
      }
    }
  }

  /*      IF(A(N,N).EQ.0.)A(N,N)=TINY */
  if (fabs(a[n-1 + (n-1) * ndim]) <= 1.0e-20) 
  {
    if (a[n-1 + (n-1) * ndim] > 0.0) 
    {
      a[n-1 + (n-1) * ndim] = 1.0e-20;
    }
    else 
    {
      a[n-1 + (n-1) * ndim] = -1.0e-20;
    }
  }
  return;
} /* ludcmp */

void lubksb(const double* a, const int n, const int ndim, 
	    const int* indx, double* b)
{
  /* System generated locals */
  int i, j, i1, ii, ll;
  double sum;

  
  /* Function Body */
  ii = -1;
  i1 = n;
  for (i = 0; i < i1; ++i) 
  {
    ll = indx[i];
    sum = b[ll];
    b[ll] = b[i];
    if (ii != -1) 
    {
      for (j = ii; j < i; ++j) 
      {
	sum -= a[i + j * ndim] * b[j];
      }
    } 
    else 
    {
      if (sum != 0.0) 
      {
	ii = i;
      }
    }
    b[i] = sum;
  }

  for (i = n-1; i >= 0; --i) 
  {
    sum = b[i];
    if (i < n-1) 
    {
   
      for (j = i + 1; j < n; ++j) 
      {
	sum -= a[i + j * ndim] * b[j];
      }
    }
    b[i] = sum / a[i + i * ndim];
  }
  return;
} /* lubksb */

