/*******************************************
Created by J. Marriner.
Installed into snana v8_38, January 2010. 

Program to take output from the SALT fitter dict files and
1.  Determine alpha and beta parameters
2.  Output a file of bias-corrected distances for cosmological fits

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

file = comma-sep list of fitres file names to analyze
     
nmax=100                 ! fit first 100 events only
nmax=70(SDSS),200(PS1MD) ! fit 70 SDSS and 200 PS1MD
nmax=300,200(PS1MD)      ! fit 300 total, with 200 in PS1MD sub-sample

sigint_fix=0.11            ! fix sigint=0.11 for all data
sigint_fix=0.11,0.09,0.08  ! comma-sep list of sigint_fix for each IDSAMPLE

maxerr_abort_c=0      ! for USESIM_c=T,  increase this to epsilon>0
maxerr_abort_x1=0     ! for USESIM_x1=T, increase this to epsilon>0
maxerr_abort_x0=0

 - - - - -  biasCor options - - - - - 
simfile_biascor=name    ! sim fitres file to compute bias mape
simfile_biascor=name1,name2,etc   ! idem with comma-sep list

opt_biascor=<option>    ! grep MASK_BIASCOR  SALT2mu.c | grep define
opt_biascor=112    ! as in BBC paper

opt_biascor=368    ! idem, but use SigmaInt(m0,x1,c)=S3x3 from biasCor
                   ! Use S3x3 * SCALE, where SCALE is varied for chi2/dof=1.
                   ! Note that SCALE is sort'of the analog of sigint.

opt_biascor=240   ! +=128 to correct mu instead of correcting mB,x1,c

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
interp_biascor_logmass=1           ! allows turning OFF logmass interpolation

dumpflag_nobiascor=1;   ! dump info for data with no valid biasCor (up to 10 events)
frac_warn_nobiascor=0.02  ! print warning in output fitres file if nobiascor
                          ! cut-loss exceeds this fraction (applies to each IDSAMPLE)

To check sample stats for each surveyGroup and fieldGroup,
   grep IDSAMPLE  <stdout_file>

 - - - - -  CCprior options - - - - -  
simfile_ccprior=name    ! sim fitres file to compute CC prior and
                        !  flag to to a BEAMS-like fit
simfile_ccprior=same  --> use same file as simfile_bias
simfile_ccprior=name1,name2,etc  ! comma-sep list
simfile_ccprior=H11   --> no sim; use CC MU-z function from Hlozek 2011
  BEWARE: must use 5D biasCor with simfile_ccprior option;
          1D biasCor + ccprior --> abort.

varname_pIa=name of fitres param containing Prob_Ia
force_pIa=forced value of Prob_Ia for every event
force_pIa='perfect' --> pIa = 1 or 0 for true Ia or CC (sim only)

maxprobcc_for_sigint --> compute sigInt from chi2=Ndof for 
                         this Ia-like subset

# to force P(SNIa)=1 and P(CC)=0 for spectroscopic-confirmed subset,
type_list_probcc0=1,2,11 
     (list of integer TYPE values in data header)
idsurvey_list_probcc=5,50,51,53 
       or
idsurvey_list_probcc=CSP,JRK07,KAIT,CFA3
       or
idsurvey_list_probcc=CSP,JRK07,51,53
   (list of survey names and/or integers from SURVEY.DEF file)

    grep Force <stdout>  # to verify parsing

# to allow for missing data columns in some input files, need to append 
# a value of -9 in the FITRES output to avoid mis-aligned output columns.
# Default will append varname_pIa and anything with PROB_ prefix to allow 
# for multiple classification probs in photometric sample, while missing
# in the spec-sample.
append_varname_missing = 'PROB_*'         ! default: wildcard for PROB_* 
append_varname_missing = 'PROB_NOTSOGOOD' ! only this one varname
append_varname_missing = 'PROB_*,PIA_*'   ! wildcard for PROB_* or PIA_*
# note that varname_pIa is automatically included, even 
# if not specified in this key. 

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

zbinuser=.01,0.012,0.1,0.2,0.3,0.4   # user-defined z-bins 
 
min_per_zbin =  min number of SN in z-bin to keep z-bin (default=1)

nzbin_ccprior = number of redshift bins for CC prior

x1min = lower limit on x1 (-6.0)
x1max = upper limit on x1 (+6.0)
cmin  = lower limit on color (-6.0)
cmax  = upper limit on color (+6.0)

logmass_min = min cut on logmass (-20)
logmass_max = max  cut on logmass (+20)
nbin_logmass = number of logmass bins for BBC7D

sntype = list of types to select. Examples are
    sntype=120
    sntype=120,105
    sntype=120,105,106  ! comma-separated, no blank spaces
    sntype=-1           ! -> all types used

CUTWIN  <VARNAME>  <MIN> <MAX>
CUTWIN  FITPROB  .01 1.0
CUTWIN  SNRMAX3   8  999999
CUTWIN(NOABORT)  SIM_x1  -2 2    ! do not abort if SIM_x1 isn't there
CUTWIN(DATAONLY)    LOGMASS  5 12   ! cut on data only (not on biasCor)
CUTWIN(BIASCORONLY) BLA_LOGMASS  5 12 ! cut on biasCor (not on data)

#select field(s) for data and biasCor with
fieldlist=X1,X2   # X1 and X2
fieldlist=X3      # X3 only
fieldlist=X       # any field with X in name

chi2max  = chi2-outlier cut applied before fit. Beware that initial
           and final chi2 can differ, so allow slop in the cut.
           = -2log10(ProbIa_BEAMS + ProbCC_BEAMS); see Eq 6 of BBC paper

cutmask_write=[allowed errcode mask to write output fitres]
cutmask_write=-1  -> write everything
cutmask_write=0   -> write only selected SN used in fit (default)
cutmask_write=64  -> include SNe that fail CUTWIN 
cutmask_write=-9  -> write comments only with fit results; no SN lines

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
p13 = scalePCC  if u13=1 (scale P_CC in BEAMS-like chi2)
p13 = scalePIa  if u13=2 (= A in Eq 4 of https://arxiv.org/abs/1111.5328)
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

u1 = 1  if true, vary parameter 1 (1=True)
u2 = 2  if true, vary parameter 2 (1=True)
u3 = 1  if true, vary parameter 3 (1=True)
.
.
u13 = 1 or 2 (see p13 comments above)
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
   = 1 to float fixed M0 in each bin (default)
   = 2 to float M0 as knot with interpolation

fixpar_all=1 --> internally set all float logicals to false, even if
                 they are set true in the input file.
                  (i.e., uM0=0, u1=0, u2=0, etc ...)
                This option turns BBC into a distance calculator.

blindflag=0  --> turn off blinding
blindflag=1  --> add cos(10*z)  to MUDIF(z)
blindflag=2  --> DEFAULT: apply random shift to ref OL,w0 (for data only)
blindflag=66 --> 2+64: apply blindflag=2 option to sim & real data

With default blindflag=1, user can set blinding parameters with
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
JOBID_SPLITRAN=[JOBID] 
   --> do only this splitran job, JOBID=1,2 ... NSPLITRAN
   --> write summary file if JOBID > NSPLITRAN

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

     In chi2 loop (fcn), skip SN that fail cut if cutmask_write=0
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
   + cutmask_write=-9 --> do NOT write SN lines to output fitres file.

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

 May 20 2019: 
   + refactor biasCor to read comma-separated list of simFile_biasCor 
   + same for CCprior.

 May 29 2019:
   + new input JOBID_SPLITRAN to process only 1 of the NSPLITRAN jobs.
     The summary file is created if JOBID_SPLITRAN > NSPLITRAN
   + fix bugs re-initializing NSPLITRAN jobs. See SPLITRAN_prep_input(),
     and bug-fix in setup_zbins_fit().

 Jun 3 2019:
   + major refactor to apply cuts more uniformly for data,biasCor,CCprior.
     See new functions set_CUTMASK() and setbit_CUTMASK().
     Found and fixed a few inconsistencies:
       CCprior : was missing cuts on bad LC fit errors and bad COV.
       BiasCor : for sigma_int calc, global z-cut was removed.

 Jun 22 2019: fix bug using varname_pIa as CUTWIN
 Jul 08 2019: remove USECODE_LEGACY pre-proc flag.
 Jul 18 2019: remove legacy code and USE_REFACOR 
 Jul 19 2019: new chi2max input; see function applyCut_chi2()
 Jul 26 2019: new input zbinuser

 August 2019: implement 7D biasCor based on gammaDM and LOGMASS
 
 Sep 4 2019: finally remove INPUTS.LEGACY_BUGS[READFILE]

 Sep 25-26: 
    + write contamination info to fitres output. See new functions
      _contam_CCprior

 Sep 29 2019: u13=2 --> scalePCC is switched to scalePIa as in Eq 4 of H11.
                (https://arxiv.org/abs/1111.5328)

 Oct 13 2019:
   +  new option force_pIa=xx.yy
   +  in write_MUERR_INCLUDE(), include hash for comment line so that
      it's easier to parse with python.
   +  if MUDIFERR=0, write MUDIFERR_ZERO = 666. Also write
      WARNING lines in FITRES and M0DIF files (for MUDIFERR=0)
   + fix bug implementing cutmask_write
   + if cutmask_write!=0, add WARNING message to FITRES file.

 Oct 18 2019:
   + add user input keys for logmass_min, logmass_max, nbin_logmass
     Used only for 7D correction.
 Oct 22 2019: bix blinding options.
 Oct 24 2019: 
    + scalePCC bound -> 15 (was 5)
    + clean up stdout messaging in prepare_biasCor(1D,5D,6D,7D)

 Nov 14 2019: MAXBIN_BIASCOR_1D -> 500k (was 200k)
 Dec 11 2019: 
     + fix few bugs so that 1D5DCUT option works
       (i.e., apply 1D biasCor, but require 5D biasCor exists)
     + abort if length(varname_gamma) > MXCHAR_VARNAME =  60

 Jan 6 2020:
   +  use print_preAbort_banner(fnam)
   +  add ABORT logic in sort_IDSAMPLE.

 Jan 09 2020:
   + fix write_fitres to loop over all input data files, not just first one.

 Jan 17. 2020:
   + in makeMap_sigmu_biasCor, fix gamma dimension.
   + fix gamma sign errors. Does not affact BBC5D, but affects BBC7D.

 Jan 23 2020: 
   + add new input nzbin_ccprior 
   + new option force_pIa=perfect (see above)

 Feb 5 2020: ignore gammaDM dimension if >5 bins to allow for a
             physically motivated distribution

 Feb 22 2020:
   + begin refactoring with option to bias-correct MU instead of
     mB,x1,c. See MASK_BIASCOR_MU = 128, and see ILCPAR_MIN, ICLPAR_MAX.
     Invoke mu-bias correction with opt_biascor += 128.

  Feb 27 2020: MXSPLITRAN -> 1000 (was 100)

  Mar 31 2020: abort if a,b,g binning is not the same in each IDSURVEY
             see check_abg_minmax_biasCor().

  Apr 14 2020: 
    + in write_NWARN, add warning about excessive loss from biasCor cut.

  Apr 15 2020: if MINOS errors are zero, report parabolic error.

  May 10 2020: 
   + use new/refactored write_fitres that allows for different
     columns in the input data-fitres files (file= key). 
     Also see new append_varname_missing key to specify which
     columns to append for files that are missing a column.
     Initial motivation is to combined photometric samples
     that have PROB_[classifier] columns and spec samples that
     do not have these columns. 

  May 20 2020: new input fieldlist=<comma sep list>

  May 29 2020:
    + for ALL surveys lumped together (no 64 bit), fix bug where
      IDSAMPLE was not set and was thus -9; set all IDSAMPLE=0.
    + ABORT if surveygroup_biascor includes survey which does
      not exist in the data.

  Jun 8 2020: add abort traps for crazy IDSURVEY

  Jun 11 
   + write SPLITRAN_read_wfit() to check for wfit output in NSPLITRAN summary
   + fix compile warnings with Wall flag
   + allow iflag_duplicate=0 (IGNORE duplicates) for sims 
 
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

// ==============================================
// define data types to track selection cuts

#define BBC_VERSION 2

#define EVENT_TYPE_DATA     1
#define EVENT_TYPE_BIASCOR  2
#define EVENT_TYPE_CCPRIOR  3
#define MXEVENT_TYPE        4
char STRING_EVENT_TYPE[MXEVENT_TYPE][12] = 
  { "NULL", "DATA", "BIASCOR", "CCPRIOR" } ;

#define MAXPAR_MINUIT 110  // limit for MINUIT params (Sep 10 2017)

#define MXFILE_DATA     20  // max number of data files to read
#define MXFILE_BIASCOR  20  // max number of biasCor files to read
#define MXFILE_CCPRIOR  20  // max number of CCprior files to read

// Maximum number of bins
// Note that 5D biasCor array would take 30*20*20*5*5 = 300,000
// So instead we define a linear array of 50,000 to take less memory,
// especially if we need to raise some of the MAXBIN_ later.
#define MAXBIN_z              85  // 50->85 on 9/11/2017
#define MAXBIN_LOGMASS        25  // Aug 2019
#define MAXBIN_BIASCOR_FITPAR 50  // for biasCor map variables x1,c
#define MAXBIN_BIASCOR_ALPHA   2  // for biasCor map
#define MAXBIN_BIASCOR_BETA    2  // for biasCor map
#define MAXBIN_BIASCOR_GAMMADM 2  // for biasCor map
#define MAXBIN_BIASCOR_1D  500000  // max on total 5D bins for {a,b,z,x1,c}
#define MINPERCELL_BIASCOR    10  // min per cell for RMS and mBslope corr

// shorter names for local declarations
#define MXz   MAXBIN_z
#define MXm   MAXBIN_LOGMASS
#define MXa   MAXBIN_BIASCOR_ALPHA
#define MXb   MAXBIN_BIASCOR_BETA
#define MXg   MAXBIN_BIASCOR_GAMMADM
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
//#define MASK_BIASCOR_2x2ab  128     // bit7: force 2x2 ab grid (min,max ab)
#define MASK_BIASCOR_MU     128     // bit7: bias on MU instead of mB,x1,c
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

#define MUDIFERR_EMPTY 999.0   // if no events in bin, set error to 999
#define MUDIFERR_ZERO  666.0   // if MUDIFFERR=0, set to 666
int NWARN_MUDIFERR_ZERO ;
int NWARN_MUDIFERR_EMPTY ;


char PATH_SNDATA_ROOT[MXPATHLEN];

// hard wire bins to store biasCor
//                                    mB    x1     c
double  BIASCOR_MINVAL_LCFIT[3]  = {  5.0, -4.0, -0.30 } ;
double  BIASCOR_MAXVAL_LCFIT[3]  = { 30.0, +4.0, +0.50 } ;
double  BIASCOR_BINSIZE_LCFIT[3] = { 25.0,  0.5,  0.05 } ;
char    BIASCOR_NAME_LCFIT[4][4] = { "mB", "x1", "c", "mu" } ;
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
#define MXSPLITRAN 1000

// define indices for fit parameters
#define NLCPAR    3  // mB,x1,c
#define INDEX_mB  0
#define INDEX_x1  1
#define INDEX_c   2
#define INDEX_mu  3  // Feb 2020

#define INDEX_z   5  // for internal flags,  not for arrays

#define MXCUTWIN 20 // max number of CUTWIN definitions in input file.

#define MXCHAR_LINE 2000 // max character per line of fitres file


// define CUTBIT for each type of cut (lsb=0)
// (bit0->1, bit4->16, bit6->64, bit8->256, bit10->1024,  bit12 --> 4096
#define CUTBIT_z          0    //  zmin - zmax cut
#define CUTBIT_x1         1    //  x1min-x1max cut
#define CUTBIT_c          2    //  cmin - cmax cut
#define CUTBIT_logmass    3    //  logmass cut
#define CUTBIT_sntype     4    //  sntype cut 
#define CUTBIT_HOST       5    //  required host info
#define CUTBIT_CUTWIN     6    //  user-CUTWIN 
#define CUTBIT_BADERR     7    //  a fit param error is <=0
#define CUTBIT_BADCOV     8    //  bad COV for a fit param
#define CUTBIT_MINBIN     9    //  z-bin with too few to fit 
#define CUTBIT_SPLITRAN  10    //  random subsample for NSPLITRAN
#define CUTBIT_SIMPS     11    //  SIM event from pre-scale 
#define CUTBIT_BIASCOR   12    //  valid biasCor (data only)
#define CUTBIT_zBIASCOR  13    //  BIASCOR z cut
#define CUTBIT_x1BIASCOR 14    //  BIASCOR x1 cut
#define CUTBIT_cBIASCOR  15    //  BIASCOR c cut
#define CUTBIT_TRUESNIa  16    //  true SNIa 
#define CUTBIT_TRUESNCC  17    //  true !SNIa (i.e., SNCC,SNIax, etc ...)
#define CUTBIT_NMAXCUT   18    //  NMAX-sample cut
#define CUTBIT_DUPL      19    //  duplicate 
#define CUTBIT_IDSAMPLE  20    //  IDSAMPLE
#define CUTBIT_CHI2      21    //  chi2 outlier (data only)
#define CUTBIT_CID       22    //  for specifying a list of cids to process
#define CUTBIT_FIELD     23    //  see fieldlist= input
#define MXCUTBIT         25  

#define dotDashLine "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" 
#define dashLine    "- - - - - - - - - - - - - - - - - - - - - - - - - - - - "

#define MAXFITITER 20             // max number of fit iterations
#define FITFLAG_DONE         0    // no more fits
#define FITFLAG_CHI2         1    // repeat fit with nominal chi2
#define FITFLAG_2LOGSIGMA    2    // repeat fit with chi2 + 2log(sigma)

#define MAXMUBIN 20

#define NCOSPAR 4  // size of cosPar array (OL,Ok,w0,wa)

#define MXNUM_SAMPLE  25  // max number of SURVEY/FIELD samples (max IDSAMPLE)
#define MXCHAR_SAMPLE 100 // max string length of sample name
#define USERFLAG_SURVEYGROUP_SAMPLE  1  // bookkeeping for biasCor IDSAMPLE
#define USERFLAG_FIELDGROUP_SAMPLE   2
#define AUTOFLAG_SURVEYGROUP_SAMPLE  3  // survey automatically added
#define USERFLAG_IGNORE_SAMPLE       6  // ignore this sample

// ---------------------
double LOGTEN  ;

// define CUTBIT and CUTMASK variables to determine events to
// reject, and to track statistics.
int  NSTORE_CUTBIT[MXEVENT_TYPE][MXCUTBIT] ;
int  CUTMASK_LIST[MXCUTBIT];          // ERRMASK for each bit
char CUTSTRING_LIST[MXCUTBIT][40];    // text definition per cut
int  *CUTMASK_POINTER[MXEVENT_TYPE];
int  *NALL_CUTMASK_POINTER[MXEVENT_TYPE];
int  *NPASS_CUTMASK_POINTER[MXEVENT_TYPE];
int  *NREJECT_CUTMASK_POINTER[MXEVENT_TYPE];
int  NDATA_BIASCORCUT[2][MXNUM_SAMPLE]; // track biascor reject for data

int MAXSN ;
int NJOB_SPLITRAN; // number of random split-jobs

int    NSIMDATA,NSIMCC ;  // used to implement prescale_sim[Data,CC]
double RNORM;
double dTWOPI;
double PIFAC; 

// ----------- TYPEDEF STRUCTURES -------------

typedef struct {
  double VAL[NLCPAR+1] ;     // bias value on mB,x1,c ... and mu
  double ERR[NLCPAR+1] ;     // error on above
  double RMS[NLCPAR+1] ;     // RMS spread of biasCor
} FITPARBIAS_DEF ;



// Jun 4 2019: define general data structure for DATA,BIASCOR,CCPRIOR

typedef struct {
  int  EVENT_TYPE ;      // DATA, BIASCOR, or CCPRIOR
  char INPUT_FILE[MXFILE_DATA][MXCHAR_FILENAME]; // input file name(s)
  int  NVAR[MXFILE_DATA];   // number of variables read per file
  char **VARNAMES_LIST[MXFILE_DATA];  // varnames per [file][ivar]

  // global counters
  int  LEN_MALLOC ;    // total number of lines in file = max for malloc
  int  NSN_ALL ;       // total number of SN rows read from file
  int  NSN_PASSCUTS ;  // number passing cuts
  int  NSN_REJECT;     // number rejected by cuts
  int  NSN_CUTBIT[MXCUTBIT]; // number cut per CUTBIT

  int  EVENT_RANGE[MXFILE_DATA][2]; // first/last event for each data file

  // variables read directly from table file.
  // use float to save memory, particularly for biasCor.
  char   **name, **field ; 
  float  *fitpar[NLCPAR+1], *fitpar_err[NLCPAR+1], *x0, *x0err ;
  float  *COV_x0x1, *COV_x0c, *COV_x1c ;
  float  *zhd,    *zhel,    *vpec ;
  float  *zhderr, *zhelerr, *vpecerr ;
  float  *logmass, *pIa, *snrmax ;
  short int  *IDSURVEY, *SNTYPE, *OPT_PHOTOZ ; 
  float  *fitpar_ideal[NLCPAR+1], *x0_ideal;

  // SIM_xxx read from table file
  short int *SIM_NONIA_INDEX ;
  float *SIM_ALPHA, *SIM_BETA, *SIM_GAMMADM;
  float *SIM_X0, *SIM_FITPAR[NLCPAR+1]; 
  float *SIM_ZCMB, *SIM_VPEC, *SIM_MU ;

  // scalar flags & counters computed from above
  bool   IS_DATA, IS_SIM ;
  int    DOFLAG_CUTWIN[MXCUTWIN]; // flag to apply each cutwin
  int    ICUTWIN_GAMMA;           // CUTVAL index for gamma-variable
  int    IVAR_VPEC, IVAR_SIM_VPEC, IVAR_OPT_PHOTOZ ; // logicals
  int    IVAR_SNTYPE, IVAR_SIM_GAMMADM;
  int    IVAR_pIa[MXFILE_DATA]; // track this variable for each file 

  int    NSN_PER_SURVEY[MXIDSURVEY] ;
  float  zMIN_PER_SURVEY[MXIDSURVEY] ;
  float  zMAX_PER_SURVEY[MXIDSURVEY] ;

  // quantities determined from table var
  float  *CUTVAL[MXCUTWIN];  // used only to make cuts.
  short int *IZBIN, *IDSAMPLE ;
  int       *CUTMASK;
  float     *mumodel;

  // covariance matrix(mB,x1,c)
  bool   *warnCov;
  float  ***covmat_tot; // fit + intrinsic scatter matrix
  float  ***covmat_fit; // covmat from LCfit (no intrinsic scatter)
  int     NCOVFIX;

} TABLEVAR_DEF ;


typedef struct { double VAL[NLCPAR][NLCPAR]; } COV_DEF ;

typedef struct {
  // parameters used for bias correction
  double z, logmass, alpha, beta, gammadm, FITPAR[NLCPAR+1];
  int idsample ;  //specifies SURVEY/FIELDGROUP
} BIASCORLIST_DEF ;


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
  int ia_min, ia_max, ib_min, ib_max, ig_min, ig_max ;
  double WGT[MXa][MXb][MXg];
} INTERPWGT_AlphaBetaGammaDM ;
 


typedef struct {
  // things fixed at start
  BININFO_DEF  ZBIN ;
  BININFO_DEF  DMUBIN ;

  // things that vary in fit
  double    alpha, beta, gammadm, M0 ;
  double    cosPar[NCOSPAR] ;

  // PDF along DMU, in each z bin and in each FITCLAS(Ib,II,CC)
  double DMUPDF[MXNUM_SAMPLE][MXz][MAXMUBIN]; 
  int    NDMUPDF[MXNUM_SAMPLE][MXz][MAXMUBIN]; 

  // MEAN and RMS of DMU in each z-bin, to approximate Gaussian
  double DMUAVG[MXNUM_SAMPLE][MXz] ;
  double DMURMS[MXNUM_SAMPLE][MXz] ;
  
  // CC params used in Hlozek 2011
  double H11_fitpar[10] ;  
} MUZMAP_DEF ;


// define CELL-INFO for multi-dimensional storage of biasCor.
// Actual biasCor stored elsewhere, since different kinds of
// info are stored.
typedef struct {
  // define 1D info separately for each dimension
  BININFO_DEF BININFO_z ;             // redshift bining
  BININFO_DEF BININFO_m ;             // logmass bining 
  BININFO_DEF BININFO_LCFIT[NLCPAR] ; // binning for mB, x1, c

  // collapse 5D space (a,b,z,x1,c)  into 1D space
  int     NCELL ;  // number of cells to store info 
  int     MAPCELL[MXa][MXb][MXg][MXz][MXm][MXpar][MXpar]; // J1D=0 to NCELL-1
  //  int  MAPCELL[MXa][MXb][MXg][MXz][MXm][MXpar][MXpar]; 

  // define arrays which depend on SURVEYGROUP/FIELDGROUP sub-sample
  int     *NperCell ;           // Number of events in each cell, vs. J1D
  double  *AVG_z;               // wgt-avg z in each cell, vs J1D
  double  *AVG_m;
  double  *AVG_LCFIT[NLCPAR];   // idem for mB,x1,c, vs. J1D
} CELLINFO_DEF ;




typedef struct {

  // WARNING: update copy_IDSAMPLE_biasCor for each new item

  // group of fields: e.g, '82N'
  char NAME_FIELDGROUP[MXCHAR_SAMPLE] ; 

  // group of surveys: e.g.,CSP+CFA
  char NAME_SURVEYGROUP[MXCHAR_SAMPLE]; 

  // string-option for each sample (user-defined binning)
  char STRINGOPT[MXCHAR_SAMPLE];

  int  OPT_PHOTOZ ;

  // Full name of each sample. e.g, 'SDSS(82N)'
  char NAME[MXCHAR_SAMPLE];  // full name of sample
  int  IDSURVEY;             // int ID of survey for this sample

  int  NSN[MXEVENT_TYPE]; // per DATA, BIASCOR, CCPRIOR

  // Dec 2017: store sparse list of rows for each IDSAMPLE with cuts
  int  NBIASCOR_CUTS;
  int *IROW_CUTS ;  // points to rawdata biasCor row

  double zMIN_DATA ; // zMIN among data in SAMPLE
  double zMAX_DATA ; // zMAX among data in SAMPLE

  // define CELLINFO [ranges and Nbin]
  double  BINSIZE_REDSHIFT, RANGE_REDSHIFT[2]; 
  double  BINSIZE_LOGMASS,  RANGE_LOGMASS[2]; 
  double  BINSIZE_FITPAR[NLCPAR] ;

  int DOFLAG_SELECT ;  // 1 --> select this sample for fit (default=all)
  int DOFLAG_BIASCOR ; // 1--> apply biasCor to this sample

  int IFLAG_ORIGIN ; // see USERFLAG_XXX

  // WARNING: update copy_IDSAMPLE_biasCor for each new item

} SAMPLE_INFO_DEF ;


struct {
  TABLEVAR_DEF TABLEVAR;

  double *mumodel, *M0, *mu, *muerr, *muerr_raw,  *mures, *mupull;
  double *muerr_last, *muerrsq_last, *sigCC_last, *sqsigCC_last ;
  double *muCOVscale, *muBias, *muBiasErr, *muBias_zinterp; 
  double *chi2, *probcc_beams;
  double **fitParBias;

  // before fit, store bias[isn][ialpha][ibeta][igammadm]
  FITPARBIAS_DEF ****FITPARBIAS_ALPHABETA ; 
  
  // before fit, store muCOVscale at each alpha,beta,gamma bin
  double ****MUCOVSCALE_ALPHABETA ; 

  float MEMORY;  // Mbytes
} INFO_DATA;


struct {
  int ILCPAR_MIN, ILCPAR_MAX ; // either mB,x1,c  OR just mu
  TABLEVAR_DEF TABLEVAR ;

  int  NDIM;     // number of dimensions for biasCor e.g., 1, 5, 7
  char STRING_PARLIST[20]; // e.g., z,x1,c,a,b

  // alpha,beta,z binning
  int8_t  *IA, *IB, *IG; // store alpha,beta,gamma index for each event
  int8_t  *IZ, *iz ;     // store redshift index for each event

  BININFO_DEF BININFO_SIM_ALPHA ; 
  BININFO_DEF BININFO_SIM_BETA ;
  BININFO_DEF BININFO_SIM_GAMMADM ;
  double      GAMMADM_OFFSET ;

  // bias in 5D cells, for each fitPar (mB,x1,c)
  // 5D cells are defined vs. J1D = 0 to CELLINFO_MUCOVSCALE.NCELL-1
  FITPARBIAS_DEF  **FITPARBIAS ;  // vs. SAMPLE and 5D cell

  // for biasCor sample; needed to determine bias on sigma_mu
  double SIGINT_ABGRID[MXNUM_SAMPLE][MXa][MXb][MXg] ; 
  double SIGINT_AVG;            // averaged over a,b bins
  
  // define intrinsic scatter matrix 
  int     NEVT_COVINT[MXNUM_SAMPLE][MXz][MXa][MXb][MXg] ;
  COV_DEF COVINT[MXNUM_SAMPLE][MXz][MXa][MXb][MXg] ;
  COV_DEF COVINT_AVG[MXNUM_SAMPLE];  // avg over a,b,z (for printing)
  BININFO_DEF zCOVINT;

  // define bias-scale on muCOV (use float to save memory)
  float **MUCOVSCALE ; 

  // wgted-avg redshift  corresponding to each M0 offset (for M0 outfile)
  double zM0[MXz];       

  float MEMORY;  // Mbytes

} INFO_BIASCOR;


struct {
  TABLEVAR_DEF TABLEVAR;       // entire file
  TABLEVAR_DEF TABLEVAR_CUTS;  // subset passing cuts (to free above memory)

  // flags transferred from INPUTS
  int USE ;      // use sim CC as prior
  int USEH11 ;   // use Hlozek CC prior if simfile_ccprior='H11'

  FITPARBIAS_DEF  ****FITPARBIAS_ALPHABETA ;
  FITPARBIAS_DEF  ****FITPARBIAS_ALPHABETA_CUTS ;

  // computed quantities
  MUZMAP_DEF  MUZMAP ;
  double PCC_biasCorScale[MXNUM_SAMPLE] ; // scale Prob_CC due to biasCor cut

  float MEMORY;  // Mbytes
} INFO_CCPRIOR;


// ------------ GLOBAL STRUCTURES ----------------

CELLINFO_DEF *CELLINFO_BIASCOR ;  
CELLINFO_DEF *CELLINFO_MUCOVSCALE ; 


// misc global info about data, broken down by survey
struct {
  int    NEVT_PER_SURVEY[MXIDSURVEY];
  double zMIN[MXIDSURVEY];
  double zMAX[MXIDSURVEY];

  double sigint_fix[MXIDSURVEY]; // see siginit_fix input (7.05.2018)
  double sqsigint_fix[MXIDSURVEY];
  
} SNDATA_INFO ;



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


#define MXPROBCC_ZERO 40
struct {
  bool USE ;
  char str_type_list[100];        // e.g., 119,120
  char str_idsurvey_list[200];   // e.g., CFA3,CFA4,CSP

  int  ntype, nidsurvey;
  int  type_list[MXPROBCC_ZERO];      // integer types
  int  idsurvey_list[MXPROBCC_ZERO];  // integer IDSURVEYs
} INPUTS_PROBCC_ZERO;


#define MXVARNAME_MISSING 10
struct {
  int    ndef ; 
  char   *varname_list[MXVARNAME_MISSING];
  bool   wildcard[MXVARNAME_MISSING];
} INPUTS_VARNAME_MISSING; // extracted from input append_varname_missing='xxx'

#define ORDER_ZPOLY_COVMAT 3   // 3rd order polynom for intrinsic scatter cov.
#define MXSURVEY_ZPOLY_COVMAT 10 

struct INPUTS {
  int  nfile_data;
  char **dataFile;  

  bool   cat_only;    // cat fitres files and do nothing else
  char   catfile_out[MXCHAR_FILENAME] ;

  int    nmax_tot ;   // Nmax to fit for all
  int    nmax[MXIDSURVEY];   // idem by survey
  char   nmaxString[100];

  double maxerr_abort_c, maxerr_abort_x1, maxerr_abort_x0;

  int    blindFlag  ; // suppress cosmo param printout for data
  double blind_cosinem0  ;       // M0DF -> M0DIF + cos(z*blind_cosinem0)
  double blind_cosinePar[MAXPAR][2]; // blind offset = [0] * cos([1])
  char   blindString[MAXPAR][60];

  // optional sim-files with corrections/maps
  int  nfile_biasCor;        // number of biascor files
  char **simFile_biasCor ;   // sim info for BBC method
  int  opt_biasCor ;
  int  prescale_biasCor[2] ; // subset use [0] of [1]
  double sigint_biasCor ;     // force sigint instead of autocompute
  double snrmin_sigint_biasCor; // SNRMIN used to determine sigint
  double sigma_cell_biasCor; // to weight events in cell for biasCor value

  char fieldGroup_biasCor[200]; // sub-divide biasCor groups based on field
  int  use_fieldGroup_biasCor;  // logical flag for above

  char surveyGroup_biasCor[400]; // combine surveys into group
  int  use_surveyGroup_biasCor;

  char surveyList_noBiasCor[200]; // list of surveys fit, but skip biasCor
  char idsample_select[40];       // e.g., '0+3'

  int interp_biascor_logmass;
  // ----------
  int  nfile_CCprior;
  char **simFile_CCprior;    // to get CC prior, dMU vs. z
  char   varname_pIa[100];

  char   append_varname_missing[100]; // force missing varname(s) with -9

  double force_pIa;
  bool   perfect_pIa;       // internally set if force_pIa='perfect'
  int  typeIa_ccprior ;       // PCC=0 for this sntype
  double maxProbCC_for_sigint;  // max P_CC/ProbIa to sum chi2_1a

  int    fitflag_sigmb ;  // flag to fit for sigMb that gives chi2/dof=1
  double redchi2_tol ;    // tolerance in red chi2 (was sig1tol)
  double sigmb_step1 ;        // size of first sigMB step, OR ...
  double scale_covint_step1;  // size of first scale_covint step
  double covint_param_step1;  // one of the above

  char   sigint_fix[200];

  double prescale_simData ;   // prescale for simData (never real data)
  double prescale_simCC ;     // e.g., 2 --> use only 1/2 of sim CC

  // - - - - - -  cuts - - - - - 
  double cmin,  cmax  ;
  double x1min, x1max ;
  double zmin,  zmax  ;
  double logmass_min, logmass_max ;
  int    nbin_logmass;
  double chi2max;     // chi2-outlier cut (Jul 2019)

  int    NCUTWIN ;
  char   CUTWIN_NAME[MXCUTWIN][MXCHAR_VARNAME];
  double CUTWIN_RANGE[MXCUTWIN][2];

  int    NFIELD ;
  char   *FIELDLIST[MXFIELD_OVERLAP] ;

  bool   LCUTWIN_RDFLAG[MXCUTWIN] ; // 1=> read, 0=> use existing var
  bool   LCUTWIN_ABORTFLAG[MXCUTWIN] ;  // 1=> abort if var does not exist
  bool   LCUTWIN_DATAONLY[MXCUTWIN] ;   // 1=> cut on real or sim data 
  bool   LCUTWIN_BIASCORONLY[MXCUTWIN]; // 1=> cut on biasCor

  int  Nsntype ;
  int  sntype[MXSNTYPE]; // list of sntype(s) to select
  char sntypeString[100]; // comma-separated string from input file

  // CID select for data only (data can be real or sim) djb
  char   cidFile_data[MXCHAR_FILENAME];
  char  **cidList_data; //2d array to be allocated later
  int   ncidList_data; //number of cids provided in listfile
  

  // - - - - - - redshift bins - - - - - - 
  int     nzbin ;    // number of uniform-spaced redshift bins
  int     nlogzbin;  // number of log-spaced redshift bins
  double  powzbin ;  // fit & output zbin ~ (1+z)^powzbin
  double  znhalf ;   // z where half of nzbins are log and half are constant
  char    zbinuser[MXPATHLEN]; // e.g., 0.01,0.04,0.01,0.3,0.7

  int      nzbin_ccprior; // number of z bins for CC prior (default=4)
  BININFO_DEF BININFO_z ; // Aug 20 2016
  int     min_per_zbin ;

  char varname_z[100]; // name of redshift variable (default = 'z Z zSPEC')

  // - - - - - user parameter bounds - - - - - 
  double parval[MAXPAR];
  double parstep[MAXPAR];
  double parbndmin[MAXPAR]; // min par bound
  double parbndmax[MAXPAR]; // max par bound
  int    ipar[MAXPAR];      // boolean float flag for each param
  int    izpar[MAXPAR];   // iz index or -9 
  int    fixpar_all ;     // global flag to fix all parameters	      

  char varname_gamma[MXCHAR_VARNAME]; // name of variable to fit gamma HR step;
                           // e.g, "HOST_LOGMASS" or "SSFR", etc...
  int USE_GAMMA0 ;   // true if p5 is floated or fixed to non-zero value
  int cutmask_write; // mask of errors to write SN to output

  // define fixed intrinsic scatter matrix
  double sigmB, sigx1, sigc; // intrinsic sigma on mB, x1 ,c
  double xi01, xi0c, xi1c;   // Correlation coeffic

  double zpecerr ;     // replace obsolete pecv (Sep 9 2016)
  double lensing_zpar; // describes lensing constribution to sigma_int

  int    uave;       // flag to use avg mag from cosmology (not nommag0)
  int    uM0 ;       // flag to float the M0 (default=true)
  int    uzsim ;     // flag to cheat and use simulated zCMB (default=F)
  double nommag0 ;   // nominal SN magnitude
  double H0 ;

  int    FLOAT_COSPAR ;    // internal: TRUE if any COSPAR is floated
  double COSPAR[NCOSPAR];  // OL, Ok, w0, wa: input to dlcosmo

  int zpolyflag;   // z-dependent cov matrix JLM AUG 15 2012

  int  NDUMPLOG; // number of SN to dump to flog file

  int NSPLITRAN ;       // number of random subsets to split jobs
  int JOBID_SPLITRAN ;  // do only this JOBID among NSPLITRAN

  int iflag_duplicate;

  char   PREFIX[100] ; // out file names = [PREFIX].extension
  
  int dumpflag_nobiasCor; // 1 --> dump info for data with no valid biasCor
  double frac_warn_nobiasCor; // give warning of nobiasCor frac exceeds this

  char SNID_MUCOVDUMP[MXCHAR_VARNAME]; // dump MUERR info for this SNID

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
  "Omega_L     ", "Omega_k     ", "w0          ","wa          ",
  "scalePCC    ", "sigint      ", "alphaHost   ","betaHost    ",
  "H11mucc0    ", "H11mucc1    ", "H11mucc2    ",
  "H11sigcc0   ", "H11sigcc1   ", "H11sigcc2   ",
  "blank4      "
} ;


// define alpha & beta bound for fit and for biasCor extrapolation
double FITPARBOUND_ALPHA[2] = {  0.02, 0.30 } ;
double FITPARBOUND_BETA[2]  = {  1.00, 6.00 } ;
double FITPARBOUND_GAMMA[2] = { -0.50, 0.50 } ;

int IPAR_ALPHA0, IPAR_BETA0, IPAR_GAMMA0, IPAR_GAMMA1;
int IPAR_LOGMASS_CEN, IPAR_LOGMASS_TAU ;
int IPAR_scalePCC, IPAR_H11, NPAR_H11_TOT, NPAR_H11_USER ;
int IPAR_COVINT_PARAM ;  // sigint or SCALE_COVINT(biasCor)
int IPAR_OL, IPAR_Ok, IPAR_w0, IPAR_wa;


// define inputs to fit that are read or calculated
struct {
  int    NSNTOT ;            // total number of SN read
  int    NPASS, NREJECT ;     // temporary during refactor
  double COVINT_PARAM_FIX ;  //  sigint OR SCALE_COVINT
  double COVINT_PARAM_LAST; 
  char   COVINT_PARAM_NAME[20];

  double ZPOLY_COVMAT[3][3];   // z-dependent scatter matrix

  int NZBIN_TOT[MXz]; // total number before cuts
  int NZBIN_FIT[MXz]; // number after cuts

  int NFITPAR_ALL ; // fixed + floated
  int NFITPAR_FLOAT ;
  int NFITPAR_FLOAT_z   ;  // number of floated M0(z) bins

  int ISFLOAT_z[MXz] ; // 1->floated, vs. z-index
  int ISFLOAT[MAXPAR];      // 1-> floated in fit; 0->fixed, vs. ipar index
  
  int IPARMAP_MN[MAXPAR];  // map minuit ipar to user ipar
} FITINP ; 


typedef struct {
  char        COMMENT[MXCHAR_VARNAME];
  BININFO_DEF BININFO;

  // start with SUMPROB summed over full sample
  double      SUMPROB_TOT_IA;    
  double      SUMPROB_TOT_CC;    
  double      SUMPROB_TOT_RATIO; 

  // now in bins. MXz is max for whatever variables ... 
  // doesn't have to be redshift bins.
  double      sumProb_Ia[MXz];   
  double      sumProb_cc[MXz];   
  double      sumProb_ratio[MXz];

  // true values for sim-data only, full sample
  int    NTRUE_TOT_IA;
  int    NTRUE_TOT_CC;
  double TRUE_TOT_RATIO;

  // and again in bins
  int    ntrue_Ia[MXz];
  int    ntrue_cc[MXz];
  double true_ratio[MXz];
} CONTAM_INFO_DEF ;

CONTAM_INFO_DEF CONTAM_MURES_BINS;
CONTAM_INFO_DEF CONTAM_REDSHIFT_BINS;



// define fit results
struct {
  int NSNFIT ;        // Number of SN used in fit; was nsnacc
  int NSNFIT_SPLITRAN[MXSPLITRAN]; // idem in SPLITRAN sub-samples
  int NDOF ;          // Ndof in fit
  int NCALL_FCN ;     // number of calls to FCN function

  double CHI2SUM_MIN;  // global min chi2
  double CHI2RED_ALL;  // global reduced chi2
  double CHI2SUM_1A;   // chi2 sum for Ia subset
  double CHI2RED_1A;   // reduced chi2 for Ia subset
  double ALPHA, BETA, GAMMA;  

  // maybe replace these with CONTAMIN_INFO ??
  double NSNFIT_1A;    // Sum of PROB(Ia)
  double NSNFIT_CC;    // Sum of PROB(CC)
  int NSNFIT_TRUECC;   // for sim, number of true CC

  double AVEMAG0 ; // average mag among z bins
  double SNMAG0 ;  // AVEMAG0, or user-input INPUTS.nommag0 

  double M0DIF[MXz]; // M0-M0avg per z-bin
  double M0ERR[MXz]; // error on above
  double zM0[MXz];   // wgted-average z (not from fit)

  char   PARNAME[MAXPAR][MXCHAR_VARNAME];
  double PARVAL[MXSPLITRAN+1][MAXPAR];
  double PARERR[MXSPLITRAN+1][MAXPAR];

  // keep track of chi2red and sigMB for each iteration
  int    NFIT_ITER ;    // Number of fit iterations (was sig1repeats)
  double CHI2RED_LIST[MAXFITITER]; // JLM AUG 15 2012
  double SIGINT_LIST[MAXFITITER];   // JLM AUG 15 2012

} FITRESULT ;


// May 2020: for multiple input data files, columns may not match.
// This OUTPUT_VARNAMES struct is a map from final ivar
// to original ivar in each file.
struct {
  int  NVAR_TOT ; 
  int  IVARMAP[MXFILE_DATA][MXVAR_TABLE];  // = original ivar from each file
  char *LIST[MXVAR_TABLE]; 

  int  NVAR_DROP;
  char DROPLIST[MXCHAR_FILENAME]; // list of dropped colums (for comment)

  bool PERFECT_COLUMN_MATCH ;  // true if all columns match
} OUTPUT_VARNAMES ;


int  DOFIT_FLAG;    // non-zero --> do another fit iteration

int ISDATA_REAL;           // T if no SIM keys are found -> real data
int FOUNDKEY_SNTYPE = 0 ;  // T -> found sntype key in fitres file
int FOUNDKEY_SIM    = 0 ;  // T -> is simulation (formerly 'simok')
int FOUNDKEY_SIMx0  = 0 ;  // T -> is SALT2 (formerly 'simx0')
int FOUNDKEY_SIM_NONIA_INDEX = 0 ;
int FOUNDKEY_OPT_PHOTOZ = 0 ;

int NVAR_ORIG ;    // NVAR from original ntuple
int NVAR_APPEND ;  // NVAR appended from SALTmu
char VARNAMES_APPEND[MAXPAR*10][MXCHAR_VARNAME] ;

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
void check_duplicate_SNID(void);
void applyCut_nmax(void);
void applyCut_chi2(void);
void merge_duplicates(int N, int *isnList);
void setup_zbins_fit(void);
void setup_BININFO_redshift(void); // setup BBC redshift bins
void setup_BININFO_powz(void);
void setup_BININFO_logz(void);
void setup_BININFO_userz(void);

void fcn(int* npar, double grad[], double* fval,
	 double xval[], int* iflag, void *);

void parse_parFile(char *parFile );
void override_parFile(int argc, char **argv);
void parse_datafile(char *item);
void parse_simfile_biasCor(char *item);
void parse_simfile_CCprior(char *item);
void parse_ZPOLY_COVMAT(char *item);
void load_ZPOLY_COVMAT(int IDSURVEY, double Z ) ;
double sum_ZPOLY_COVMAT(double Z, double *polyPar) ;

void parse_CUTWIN(char *item);
void parse_FIELDLIST(char *item);
int  apply_CUTWIN(int EVENT_TYPE, int *DOFLAG_CUTWIN, double *CUTVAL_LIST);
int  usesim_CUTWIN(char *varName) ;
int  set_DOFLAG_CUTWIN(int ivar, int icut, int isData );

void parse_sntype(char *item);
void parse_cidFile_data(char *item); //djb
void parse_prescale_biascor(char *item);
void parse_powzbin(char *item) ;
void parse_IDSAMPLE_SELECT(char *item);
void parse_sigint_fix(char *item);
void parse_blindpar(char *item) ;

void  prep_input_driver(void);
void  prep_input_nmax(char *item);
void  prep_input_gamma(void) ;
void  prep_input_probcc0(void);
void  prep_input_load_COSPAR(void);
void  prep_input_varname_missing(void);

int   force_probcc0(int itype, int idsurvey);
void  prep_cosmodl_lookup(void);
int   ppar(char* string);
void  read_data(void);
void  store_input_varnames(int ifile, TABLEVAR_DEF *TABLEVAR) ;
void  store_output_varnames(void);
bool  exist_varname(int ifile,char *varName, TABLEVAR_DEF *TABLEVAR);
void  get_zString(char *str_z, char *str_zerr, char *cast) ;
void  SNTABLE_READPREP_TABLEVAR(int ifile, int ISTART, int LEN, 
				TABLEVAR_DEF *TABLEVAR);
void  SNTABLE_CLOSE_TEXT(void) ;
void  compute_more_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR) ;
void  compute_more_INFO_DATA(void);
void  prepare_IDSAMPLE_biasCor(void); 
void  set_FIELDGROUP_biasCor(void);
void  set_SURVEYGROUP_biasCor(void);
void  set_BINSIZE_SAMPLE_biasCor(int IDSAMPLE);

int    get_IDSAMPLE(int IDSURVEY, int OPT_PHOTOZ, 
		    char *FIELDGROUP, char *SURVEYGROUP, int DUMPFLAG );
void   sort_IDSAMPLE_biasCor(void) ;
void   copy_IDSAMPLE_biasCor(SAMPLE_INFO_DEF *S1, SAMPLE_INFO_DEF *S2) ;
void   dump_SAMPLE_INFO(int EVENT_TYPE) ;
void   match_fieldGroup(char *SNID, char *FIELD, 
			char *FIELDGROUP, char *STRINGOPT ); 
void   match_surveyGroup(char *SNID, int IDSURVEY, 
			 char *SURVEYGROUP, char *STRINGOPT ) ;
int    match_fitParName(char *parName);

void   prepare_biasCor(void);
void   prepare_biasCor_zinterp(void) ;
void   init_COVINT_biasCor(void);
void   prepare_CCprior(void);

void   get_INTERPWGT_abg(double alpha,double beta,double gammadm, int DUMPFLAG,
			INTERPWGT_AlphaBetaGammaDM *INTERPWGT, char *callFun );

void   fcnFetch_AlphaBetaGamma(double *xval, double z, double logmass,
			       double *alpha, double *beta, double *gammadm );

void   read_simFile_CCprior(void);
void   store_INFO_CCPRIOR_CUTS(void);
int    storeBias_CCprior(int n) ;

void   setup_zbins_CCprior (TABLEVAR_DEF *TABLEVAR, BININFO_DEF *ZBIN) ;

void   setup_MUZMAP_CCprior(int IDSAMPLE, TABLEVAR_DEF *TABLEVAR,
			    MUZMAP_DEF *MUZMAP) ;

void   setup_DMUPDF_CCprior(int IDSAMPLE, TABLEVAR_DEF *TABLEVAR,
			    MUZMAP_DEF *MUZMAP );

void setup_contam_CCprior(void);
void zero_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO) ;
void sum_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO, double Prob_Ia,
			double xhisto, int SIM_NONIA_INDEX) ;
void print_contam_CCprior(FILE *fp, CONTAM_INFO_DEF *CONTAM_INFO);

void  dump_DMUPDF_CCprior(int IDSAMPLE, int IZ, MUZMAP_DEF *MUZMAP) ;

double prob_CCprior_sim(int IDSAMPLE, MUZMAP_DEF *MUZMAP, 
		    double z, double dmu, int DUMPFLAG);
double prob_CCprior_H11(int n, double dmu, double *H11_fitpar, 
			double *sqsigCC, double *sigCC_chi2penalty);

void   load_FITPARBIAS_CCprior(int icc, FITPARBIAS_DEF
			       (*FITPARBIAS_ABGRID)[MXb][MXg]);

int IBINFUN(double VAL, BININFO_DEF *BIN, int OPT_ABORT, char *MSG_ABORT );
void copy_BININFO(BININFO_DEF *BIN1, BININFO_DEF *BIN2);

void  setup_CELLINFO_biasCor(int IDSAMPLE);
void  setup_BININFO_biasCor(int IDSAMPLE, int ipar_LCFIT, int MAXBIN, 
			    BININFO_DEF *BININFO);
void  get_BININFO_biasCor_alphabeta(char *varName, double *VAL_MIN, 
				    double *VAL_MAX, double *VAL_BIN ) ;
void check_abg_minmax_biasCor(char *varName, double *valmin_list,
			      double *valmax_list) ;				    
void  makeMap_fitPar_biasCor(int ISAMPLE, int ipar_LCFIT);
void  makeMap_sigmu_biasCor(int ISAMPLE);   
void  vpec_biasCor(void);
void  init_sigInt_biasCor_SNRCUT(int IDSAMPLE) ;

void  makeMap_binavg_biasCor(int ISAMPLE);
void  makeSparseList_biasCor(void);

void  calc_zM0_biasCor(void);
void  calc_zM0_data(void);
double zmu_solve(double mu, double *cosPar) ;
void   test_zmu_solve(void);

int   storeDataBias(int n, int dumpFlag ) ;
int   biasMapSelect(int i) ;
void  read_simFile_biasCor(void);
void  set_MAPCELL_biasCor(int IDSAMPLE) ;
void  store_iaib_biasCor(void) ;

void  zero_FITPARBIAS(FITPARBIAS_DEF *FITPARBIAS) ;
int   get_fitParBias(char *CID, BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		     FITPARBIAS_DEF *FITPARBIAS);
int get_muCOVscale(char *cid,  BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		   double *muCOVscale ) ;

void   get_muBias(char *NAME, 
		  BIASCORLIST_DEF *BIASCORLIST,  
		  FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MXb][MXg],
		  double         (*MUCOVSCALE_ABGRID)[MXb][MXg],
		  INTERPWGT_AlphaBetaGammaDM *INTERPWGT,  
		  double *FITPARBIAS_INTERP, 
		  double *muBias, double *muBiasErr, double *muCOVscale ) ;

double get_gammadm_host(double z, double logmass, double *hostPar);

void  set_defaults(void);
void  set_fitPar(int ipar, double val, double step, 
		double bndmin, double bndmax, int use) ;

void  init_CUTMASK(void);
void  set_CUTMASK(int isn, TABLEVAR_DEF *TABLEVAR );
void  setbit_CUTMASK(int isn, int bitnum, TABLEVAR_DEF *TABLEVAR );


void  countData_per_zbin(void) ;
int   prescale_reject_simData(int SIM_NONIA_INDEX);
int   prescale_reject_biasCor(int isn);
int   outside_biasCor_grid(int isn);

int selectCID_data(char *cid); //djb

void  write_fitres_driver(char *fileName);
void  write_fitres_misc(FILE *fout);
void  write_version_info(FILE *fp) ;
void  define_varnames_append(void) ;
int   write_fitres_line(int indx, int ifile, char *line, FILE *fout) ;
void  write_fitres_line_append(FILE *fp, int indx);
void  write_cat_info(FILE *fout) ;
void  prep_blindVal_strings(void);
void  write_blindFlag_message(FILE *fout) ;
void  get_CCIDindx(char *CCID, int *indx) ;
double  BLIND_OFFSET(int ipar);

int   ISBLIND_FIXPAR(int ipar);
void  printmsg_fitStart(FILE *fp); 
void  printmsg_repeatFit(char *msg) ;
void  print_eventStats(int event_type);

void  outFile_driver(void);
void  write_M0(char *fileName);
void  write_MUERR_INCLUDE(FILE *fp) ;
void  write_NWARN(FILE *fp, int FLAG) ;

int   SPLITRAN_ACCEPT(int isn, int snid);
void  SPLITRAN_cutmask(void);
void  SPLITRAN_SUMMARY(void); 
void  SPLITRAN_write_fitpar(char *fileName);
void  SPLITRAN_read_fitpar(int isplit);
int   SPLITRAN_read_wfit(int isplit);
void  SPLITRAN_prep_input(void);

void  CPU_SUMMARY(void);

int keep_cutmask(int errcode) ;

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

void   printCOVMAT(FILE *fp, int NPAR, int NPARz_write);
double fcn_muerrsq(char *name,double alpha,double beta, double gamma,
		   double (*COV)[NLCPAR], double z, double zerr, int dumpFlag);
double fcn_muerrz(int OPT, double z, double zerr ) ;
void   fcn_ccprior_muzmap(double *xval, int USE_CCPRIOR_H11, 
			  MUZMAP_DEF *MUZMAP);
double zerr_adjust(double z, double zerr, double vpecerr, char *name );
void   test_muerrz(void);


double muerrsq_biasCor(int ievt, int maskCov, int *istat_cov) ;

void   get_COVINT_model(int idsample, double (*COVINT)[NLCPAR] ); // returns COVINT
void   get_COVINT_biasCor(int IDSAMPLE, double z, double alpha, double beta, 
			  double gammadm, 
			  double (*COVINT)[NLCPAR] ); // returns COVINT
void   dump_COVINT_biasCor(int idsample, int iz, int ia, int ib, int ig) ;
void   write_COVINT_biasCor(void);
void   zero_COV(double (*COV)[NLCPAR]);
void   scale_COV(double scale, double (*COV)[NLCPAR]);

double muresid_biasCor(int ievt);
int    J1D_biasCor(int ievt, char *msg);
void   J1D_invert_D(int IDSAMPLE, int J1D, double *a, double *b, double *g,
		    double *z, double *m, double *x1, double *c);
void   J1D_invert_I(int IDSAMPLE, int J1D, int *ia, int *ib, int *ig,
		    int *iz, int *im, int *ix1, int *ic);
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


void  malloc_INFO_DATA(int opt, int LEN_MALLOC);
void  malloc_INFO_BIASCOR(int opt, int LEN_MALLOC);
void  malloc_INFO_CCPRIOR(int opt, int LEN_MALLOC, int LEN_MALLOC_CUTS);
float malloc_TABLEVAR(int opt, int LEN_MALLOC, TABLEVAR_DEF *TABLEVAR);
int   malloc_TABLEVAR_CUTVAL(int LEN_MALLOC, int icut, 
			     TABLEVAR_DEF *TABLEVAR ) ;
float malloc_FITPARBIAS_ALPHABETA(int opt, int LEN_MALLOC, 
				  FITPARBIAS_DEF *****FITPARBIAS );

float malloc_double2D(int opt, int LEN1, int LEN2, double ***array2D );
float malloc_double3D(int opt, int LEN1, int LEN2, int LEN3, 
		      double ****array3D );
float malloc_double4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
		      double *****array4D );


// ======================================================================
// ==== GLOBALS AND FUNCTIONS TO USE SALT2mu as CHI2-INFO FUNCTION ======
// ======================================================================

void SNFUNPAR_CHI2INFO_INIT(void); // one time init (binning, malloc ...)
void SNFUNPAR_CHI2INFO_LOAD_OUTPUT(void); 
void SNFUNPAR_CHI2INFO_WRITE(void); // write output
int  SNFUNPAR_CHI2INFO_LOAD_BININFO(char *varName, double xmin, double xmax, 
				    double binSize,  BININFO_DEF *BININFO,
				    FILE *fp);

#define ISNPAR_x1 0
#define ISNPAR_c  1
#define ISNPAR_m  2
// VARNAMES: ic ix1 im NEVT MURES_SQSUM MURES_SUM
struct {
  bool  USE;
  char  OUTFILE[MXCHAR_FILENAME] ;
  FILE *FP_OUT ;
  char  LINE_VARNAMES[200];
  BININFO_DEF binInfo[3]; // c,x1,m
  double alpha, alphaErr, beta, betaErr;
  int    ***MAP3D_to_1D;
  int    NBIN_TOT ; // total number of storage bins
  int    *NEVT;  // Nevt per in NBIN_TOT bins
  double *MURES_SQSUM, *MURES_SUM; // to reconstruct MU-bias and MU-RMS
} SNFUNPAR_CHI2INFO_OUTPUT ;

#define MXSIM_SNFUNPAR 10
typedef struct {
  char PARNAME[40] ; // e.g., SIM_RV, SIM_c, etc ...
  char FUNNAME[40] ; // e.g., EXP or Gauss
  GENPOLY_DEF PARLIST[MXSIM_SNFUNPAR]; // each param can be polyfun of z
  bool ISGAUSS, ISEXP;
} SIM_SNFUNPAR_DEF ;

#define MXSIMPAR_REWGT 10
struct {
  bool USE;
  int N_BOUNDFUN;
  int N_REWGTFUN;
  SIM_SNFUNPAR_DEF  BOUNDFUN[MXSIMPAR_REWGT];
  SIM_SNFUNPAR_DEF  REWGTFUN[MXSIMPAR_REWGT];
} SIM_SNFUNPAR_STORE ;





// *********************************
int main(int argc,char* argv[ ])
{
  //History:
  // Dec 5 2016: call exec_mnparm, and cleanup local declarations

  // local declarations
  int inf=5, outf=6, savef=7;
  const int null=0;

  int icondn, len, npari, nparx, istat, ndof ;
  double chi2min, fedm, errdef ;
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
  prep_input_driver();

  //  test_zmu_solve();
  //  test_muerrz(); // xxxx

  // --------------------------------------
  //Read input data from SALT2 fit
  TABLEFILE_INIT();

  read_data(); 
  compute_more_INFO_DATA();  

  if( INPUTS.cat_only ) 
    { write_fitres_driver(INPUTS.catfile_out);  return(0); }

  // abort if there are any duplicate SNIDs
  check_duplicate_SNID();

  // setup BBC redshift bins.
  setup_BININFO_redshift();

  // prepare mapindex for each IDSURVEY & FIELD --> for biasCor
  prepare_IDSAMPLE_biasCor();

  // check option for SPLITRAN summary
  if ( INPUTS.JOBID_SPLITRAN > INPUTS.NSPLITRAN ) { goto DOFIT; }

  // read optional sim maps 
  prepare_biasCor();

  prepare_CCprior();

  // check for user-constraint on nmax (July 2017) AFTER prepare_biasCor 
  applyCut_nmax();

  // check option to turn SALT2mu into chi2-info function 
  SNFUNPAR_CHI2INFO_INIT();

  t_end_init = time(NULL);

  /* Solve for minimum of chi-squared */
 DOFIT:

  if ( INPUTS.JOBID_SPLITRAN > 0 ) 
    { NJOB_SPLITRAN = INPUTS.JOBID_SPLITRAN; } // do only this one SPLIT job
  else
    { NJOB_SPLITRAN++ ; }  // keep going to do them all

  DOFIT_FLAG = FITFLAG_CHI2 ; 
  printmsg_fitStart(stdout);
  
  if ( INPUTS.NSPLITRAN > 1 ) { 
    SPLITRAN_cutmask();     // check for random sub-samples
    SPLITRAN_prep_input();  // re-initialize a few things (May 2019)
  }

  setup_zbins_fit();    // set z-bins for data and biasCor

  FITRESULT.NCALL_FCN = 0 ;
  mninit_(&inf,&outf,&savef);

  strcpy(text,"SALT2mu"); len = strlen(text);  
  mnseti_(text,len);    fflush(stdout);

  // execuate minuit mnparm_ commands
  exec_mnparm(); 

  // check option to fetch summary of all previous SPLITRAN jobs
  if ( INPUTS.JOBID_SPLITRAN > INPUTS.NSPLITRAN ) 
    {  SPLITRAN_SUMMARY(); return(0); }

  // use FCN call and make chi2-outlier cut (Jul 19 2019)
  applyCut_chi2();

  FITRESULT.NFIT_ITER = 0 ;

  // print stats for data after ALL cuts are applied
  print_eventStats(EVENT_TYPE_DATA);
  
  // check reasons to suppress MINUIT screen dump
  int MNPRINT = 0 ;
  //  if ( ISDATA_REAL && INPUTS.blindFlag   ) { MNPRINT = 1; } 
  //  if ( IGNOREFILE(INPUTS.PREFIX)    ) { MNPRINT = 1; } 
  //  if ( INPUTS.cutmask_write == -9   ) { MNPRINT = 1; } 
  
  // Beginning of DOFIT loop
  while ( DOFIT_FLAG != FITFLAG_DONE  ) {

    if ( MNPRINT ) 
      { strcpy(mcom,"SET PRINTOUT 0"); }  // MIMUIT printing on
    else
      { strcpy(mcom,"SET PRINTOUT -1"); } // turn off MINUIT printing
    
    len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    //Miniut MINIMIZE (find minimum chi-squared)
    strcpy(mcom,"SIM 1000");   len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    strcpy(mcom,"MINI");   len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);
    //Minuit MINOS (compute errors)
    strcpy(mcom,"MINO");   len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(stdout);

    //Final call to FCN at minimum of chi-squared
    strcpy(mcom,"CALL FCN 3");  len = strlen(mcom);
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
  printCOVMAT(stdout,npari,999);
  
  // ------------------------------------------------
  // check files to write
  outFile_driver();


  //---------
  if ( NJOB_SPLITRAN < INPUTS.NSPLITRAN  &&  INPUTS.JOBID_SPLITRAN<0 ) 
    { goto DOFIT ; }

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
  //
  // Oct 22 2019: beware that mnparm_ prints blinded params.

  int i, iz, iMN, iMN_tmp, len, ierflag, icondn, ISFLOAT ;
  int nzbin = INPUTS.nzbin ;
  double M0min, M0max;
  const int null=0;
  char text[100];
  //  char fnam[] = "exec_mnparm" ;

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
      iMN_tmp = FITINP.NFITPAR_FLOAT-1;
      FITINP.IPARMAP_MN[iMN_tmp]  = i ; // map minuit ipar to user ipar
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
  char text[100], format[80], cPARVAL[MXCHAR_VARNAME] ;
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
      mnerrs_(&iMN, &eplus, &eminus, &eparab, &globcc);

      if ( fabs(eplus) > 1.0E-6 && fabs(eminus) > 1.0E-6 ) 
	{ PARERR = 0.5*( fabs(eplus) + fabs(eminus) ) ; }
      else
	{ PARERR = eparab; } // Apr 15 2020 : better than nothing
	  	
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
void setup_BININFO_redshift(void) {

  int LEN;
  char fnam[] = "setup_BININFO_redshift";

  // ------------ BEGIN ------------

  //Set up BBC z bins for data

  LEN = strlen(INPUTS.zbinuser);
  if ( LEN > 0 )  { 
    if ( LEN > MXPATHLEN - 10 ) {
      sprintf(c1err,"len(zbinuser) = %d is too long", LEN );
      sprintf(c2err,"Reduce size or increase array bound.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    setup_BININFO_userz(); 
  }

  else if ( INPUTS.nlogzbin > 0 ) 
    { setup_BININFO_logz(); }

  else if ( INPUTS.nzbin > 0 )
    { setup_BININFO_powz(); }

  else {
    sprintf(c1err,"Found no BBC z-bin option.");
    sprintf(c2err,"Check nzbin, nlogzbin, and zbinuser keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  return ;

} // end   setup_BININFO_redshift


// ***********************************************
void setup_BININFO_userz(void) {

  // Created Jul 26 2019
  // Parse zbinuser string and load z-bin structure

  int nzbin=0, Nsplit, iz ;
  int MEMC = 20*sizeof(char);
  double zlo, zhi, zmin=-9.0, zmax=-9.0 ;
  char *ptr_z[MXz];
  char comma[] = ",";
  char fnam[] = "setup_BININFO_userz" ;

  // --------------- BEGIN -------------

  for(iz=0; iz < MXz; iz++ ) { ptr_z[iz] = (char*)malloc(MEMC); }

  splitString(INPUTS.zbinuser, comma, MXz,    // inputs
	      &Nsplit, ptr_z );                    // outputs
  nzbin  = Nsplit-1 ;

  if ( Nsplit <= 1 || Nsplit >= MXz ) {
    sprintf(c1err,"Invalid Nsplit=%d for", Nsplit);
    sprintf(c2err,"zbinuser=%s\n", INPUTS.zbinuser);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  INPUTS.BININFO_z.nbin  = nzbin ;
  sprintf(INPUTS.BININFO_z.varName,"redshift");

  for(iz=0; iz < nzbin; iz++ ) {
    sscanf(ptr_z[iz+0], "%le", &zlo);
    sscanf(ptr_z[iz+1], "%le", &zhi);
    INPUTS.BININFO_z.lo[iz]  = zlo ;
    INPUTS.BININFO_z.hi[iz]  = zhi ;
    INPUTS.BININFO_z.avg[iz] = 0.5 * ( zlo + zhi ) ;

    if ( iz == 0       ) { zmin = zlo; }
    if ( iz == nzbin-1 ) { zmax = zhi; }

  }

  INPUTS.nzbin = nzbin;
  INPUTS.zmin  = zmin; 
  INPUTS.zmax  = zmax;
  return ;
  
} // setup_BININFO_userz

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
  double zmin, zmax, zbin2=0.0, arg;
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

  int nzbin = INPUTS.nzbin ;
  int NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;
  int n, nz, izbin, NZFLOAT, CUTMASK ;
  double z;
  char *name;
  char fnam[] = "setup_zbins_fit";

  // ----------- BEGIN --------
  
  // init counters
  for (nz=0; nz < nzbin; nz++)  {
    FITINP.NZBIN_TOT[nz] = 0;
    FITINP.NZBIN_FIT[nz] = 0;
  }
  
  //Put data into bins and count bin population
  for (n=0; n < NSN_DATA; ++n) {

    CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[n];
    if ( CUTMASK ) { continue; }

    z     = INFO_DATA.TABLEVAR.zhd[n] ;    
    name  = INFO_DATA.TABLEVAR.name[n] ;
    izbin = IBINFUN(z, &INPUTS.BININFO_z, 0, fnam);
    INFO_DATA.TABLEVAR.IZBIN[n] = izbin;

      
    if (izbin<0 || izbin >= nzbin)  
      {  setbit_CUTMASK(n, CUTBIT_z, &INFO_DATA.TABLEVAR);  }

    FITINP.NZBIN_TOT[izbin]++ ;
    CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[n];
    if ( CUTMASK == 0 ) { FITINP.NZBIN_FIT[izbin]++ ; }

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
  // i.e, Flag SN which are not in a valid zbin
  for (n=0; n< NSN_DATA; ++n)  {
    CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[n];
    if ( CUTMASK ) { continue; }
    izbin = INFO_DATA.TABLEVAR.IZBIN[n];

      if ( FITINP.ISFLOAT_z[izbin] == 0 )    
	{ setbit_CUTMASK(n, CUTBIT_MINBIN, &INFO_DATA.TABLEVAR );  }
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
  int  LDMP = 0 ;
  double *blindpar ;
  char fnam[] = "apply_blindpar" ;

  // -------------- BEGIN ---------------

  if ( INPUTS.cat_only ) { return; }
  if ( (INPUTS.blindFlag & BLINDMASK_FIXPAR) == 0 ) { return ; }
  
  if ( LDMP ) {
    printf(" xxx %s: blindFlag=%d ISDATAREAL=%d blindpar(LAM)=%f\n", 
	   fnam, INPUTS.blindFlag, ISDATA_REAL, 
	   INPUTS.blind_cosinePar[IPAR_OL][0] ); fflush(stdout);
  }

  for(ipar=0; ipar < MAXPAR; ipar++ ) {

    if ( ISBLIND_FIXPAR(ipar) ) {
      blindpar             = INPUTS.blind_cosinePar[ipar];
      INPUTS.parval[ipar] += (blindpar[0] * cos(blindpar[1]) ) ;
      printf("  BLIND FIXPAR: %s += %8.4f * cos(%f) \n",
	     FITPARNAMES_DEFAULT[ipar], blindpar[0], blindpar[1] );
      fflush(stdout);
    }
  }

  // re-load INPUTS.COSPAR (bugfix, Oct 22 2019)
  if ( LDMP ) {
    printf(" xxx %s: 1. COSPAR(OL,w0) = %f, %f \n",
	   fnam, INPUTS.COSPAR[0], INPUTS.COSPAR[2]); fflush(stdout);
  }

  // set INPUTS.COSPAR with blinded cosmo params
  prep_input_load_COSPAR();

  if ( LDMP ) {
    printf(" xxx %s: 2. COSPAR(OL,w0) = %f, %f \n",
	   fnam, INPUTS.COSPAR[0], INPUTS.COSPAR[2]); fflush(stdout);
  }

  return ;

} // end apply_blindpar


// ******************************************
void applyCut_nmax(void) {

  // Created July 17 2017
  // Check user nmax constraints and set skipfit
  // tro trim extra events. See nmax input and prep_input_nmax().

  int  NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  NTOT, NPERSURVEY[MXIDSURVEY]; // before nmax cuts
  int  ntot, npersurvey[MXIDSURVEY]; // after
  int  cutmask, id, isn, ntot_iter1, nrej_iter1 ;
  int  nexclude, NREMAIN ;
  //  char fnam[] = "applyCut_nmax" ;

  // --------------- BEGIN --------------

  if ( IGNOREFILE(INPUTS.nmaxString) ) { return ; }

  NTOT = ntot = ntot_iter1 = nrej_iter1 = nexclude = 0 ;
  for(id=0; id < MXIDSURVEY; id++ ) 
    {  NPERSURVEY[id] = npersurvey[id] = 0; }

  // ----- check survey-dependent nmax cuts ------------
  for(isn=0; isn < NSN_DATA ; isn++ ) {

    cutmask = (int)INFO_DATA.TABLEVAR.CUTMASK[isn];
    id      = (int)INFO_DATA.TABLEVAR.IDSURVEY[isn];
    if ( cutmask ) { continue ; }

    if ( NPERSURVEY[id] >= INPUTS.nmax[id] ) 
      { setbit_CUTMASK(isn, CUTBIT_NMAXCUT, &INFO_DATA.TABLEVAR); }

    NPERSURVEY[id]++ ;    NTOT++ ;
    cutmask = (int)INFO_DATA.TABLEVAR.CUTMASK[isn];
    if ( cutmask == 0 )
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
  for(isn=0; isn < NSN_DATA; isn++ ) {

    cutmask = (int)INFO_DATA.TABLEVAR.CUTMASK[isn];
    id      = (int)INFO_DATA.TABLEVAR.IDSURVEY[isn];
    if ( cutmask ) { continue ; }

    // if this SURVEY is already cut, don't cut more
    if ( INPUTS.nmax[id] < 999888000 ) { continue ; }

    ntot++ ;    
    if ( ntot > NREMAIN ) 
      {  setbit_CUTMASK(isn, CUTBIT_NMAXCUT, &INFO_DATA.TABLEVAR); }
  }

  // ------------------------------------------------
  // now check stats again and write out before->after for each survey
  
  SUMMARY:

  ntot=0;
  for(isn=0; isn < NSN_DATA; isn++ ) {
    cutmask = (int)INFO_DATA.TABLEVAR.CUTMASK[isn];
    id      = (int)INFO_DATA.TABLEVAR.IDSURVEY[isn];
    if ( cutmask ) { continue ; }
    ntot++ ;    npersurvey[id]++ ;
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

} // end applyCut_nmax


// ******************************************
void applyCut_chi2(void) {

  // Created Jul 19 2019
  // Call FCN (fit function) to evaluate chi2 with
  // initial params, then applhy chi2-outlier cut to data.
  // Note that this is BEAMS chi2 when there is CC contam,
  // and thus this is NOT the same as cutting on HR resids.

  int NSN_DATA   = INFO_DATA.TABLEVAR.NSN_ALL ;
  double chi2max = INPUTS.chi2max;
  int len, icondn, n, cutmask;
  const int null=0;
  double chi2;
  char mcom[60], *name ;
  //  char fnam[] = "applyCut_chi2" ;

  // ----------- BEGIN ------------

  if ( chi2max > 0.99E9 ) { return; }

  strcpy(mcom,"CALL FCN 1");  len = strlen(mcom);
  mncomd_(fcn,mcom,&icondn,&null,len);   fflush(stdout);

  for (n=0; n< NSN_DATA; ++n)  {
    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n] ; 
    if ( cutmask ) { continue; }

    chi2 = INFO_DATA.chi2[n];
    name = INFO_DATA.TABLEVAR.name[n];

    /*
    if ( chi2 > 8.0 ) {
      printf(" xxx %s: chi2(%8s) = %.2f \n", fnam, name, chi2);
      fflush(stdout);
    }
    */

    if ( chi2 > chi2max )  { 
      printf("\t Chi2(%s) = %.2f -> reject \n", name, chi2);
      setbit_CUTMASK(n, CUTBIT_CHI2, &INFO_DATA.TABLEVAR);       
    }
  }

  fflush(stdout);


  
  return ;

} // end applyCut_chi2


// ******************************************
void check_duplicate_SNID(void) {

  // Sep 2016:
  // Use sorted redshift and its error to flag duplicates.
  // CHeck iflag_duplictes for what to do:
  // 0 -> nothing
  // 1 -> abort
  // 2 -> merge fitparams and cov into one LC
  //

  int  iflag   = INPUTS.iflag_duplicate;
  int  MXSTORE = MXSTORE_DUPLICATE ;
  int  isn, nsn, MEMD, MEMI, unsort, *unsortList, ORDER_SORT   ;
  bool IS_SIM;
  double *zList ;
  char fnam[] = "check_duplicate_SNID" ;

  // ----------- BEGIN -----------

  sprintf(BANNER,"Begin %s", fnam);
  print_banner(BANNER);

  IS_SIM  = INFO_DATA.TABLEVAR.IS_SIM ;
  nsn     = INFO_DATA.TABLEVAR.NSN_ALL ;

  /* xxxxx mark delete Jun 11 2020
  // for simulation there is only one duplicate option: abort
  if ( IS_SIM == true )  { INPUTS.iflag_duplicate = IFLAG_DUPLICATE_ABORT; }
  xxxxxxxx  */


  // Jun 11 2020: for sims, don't bother tracking duplicates
  if ( IS_SIM &&  iflag == IFLAG_DUPLICATE_IGNORE ) { return; }

  MEMD = nsn * sizeof(double) + 4 ;
  MEMI = nsn * sizeof(int)    + 4 ;
  zList      = (double *)malloc(MEMD); // allocate integer ID list
  unsortList = (int    *)malloc(MEMI); // allocate sorted list

  for(isn=0; isn<nsn; isn++)  
    { zList[isn] = (double)INFO_DATA.TABLEVAR.zhd[isn]; }
  

  ORDER_SORT = + 1 ; // increasing order
  sortDouble( nsn, zList, ORDER_SORT, unsortList ) ;

  int SAME_z, SAME_SNID, NTMP, idup, unsort_last, LAST_DUPL=0 ;
  double z, zerr, z_last, zerr_last ;
  char *snid, *snid_last ;
  int  UNSORT_DUPL[MXSTORE_DUPLICATE][MXSTORE_DUPLICATE];
  int  NDUPL_LIST[MXSTORE_DUPLICATE]; // how many duplicates per set
  int  NDUPL_SET ; // number of duplicate sets

  unsort_last = -9 ;
  z_last = zerr_last = -9.0 ;
  snid_last = fnam; // anything to avoid compile warning

  NDUPL_SET = 0 ;
  for(idup=0; idup < MXSTORE; idup++ ) 
    {  NDUPL_LIST[idup] = 0 ; }

  for ( isn=0; isn < nsn; isn++ ) {
    unsort  = unsortList[isn];
    z      = INFO_DATA.TABLEVAR.zhd[unsort];
    zerr   = INFO_DATA.TABLEVAR.zhderr[unsort];
    snid   = INFO_DATA.TABLEVAR.name[unsort];

    if ( isn==0 ) { goto SET_LAST; }

    SAME_z    = ( z==z_last && zerr==zerr_last );
    SAME_SNID = ( strcmp(snid,snid_last) == 0 );

    if ( SAME_z && SAME_SNID &&  NDUPL_SET < MXSTORE-1 ) {
      if ( LAST_DUPL == 0 ) { 
	NDUPL_SET++ ; 
	NTMP = NDUPL_LIST[NDUPL_SET-1] ;
	UNSORT_DUPL[NDUPL_SET-1][NTMP] = unsort_last ;
	NDUPL_LIST[NDUPL_SET-1]++ ;
      }
      NTMP = NDUPL_LIST[NDUPL_SET-1] ;
      if( NTMP < MXSTORE )  { UNSORT_DUPL[NDUPL_SET-1][NTMP] = unsort ; }
      NDUPL_LIST[NDUPL_SET-1]++ ;
      LAST_DUPL = 1 ;
    }
    else {
      LAST_DUPL = 0 ;
    }

  SET_LAST:
    z_last    = z;   zerr_last=zerr;   unsort_last=unsort ;
    snid_last = INFO_DATA.TABLEVAR.name[unsort]; 

  } // end loop over isn


  if ( NDUPL_SET==0 ) { printf("\t No duplicates found. \n"); goto DONE; }
  
  printf("   Found %d sets of duplicates: \n", NDUPL_SET);

  for(idup=0; idup < NDUPL_SET ; idup++ ) {

    unsort = UNSORT_DUPL[idup][0] ; // first duplicate has SNID and z
    NTMP   = NDUPL_LIST[idup] ;
    snid   = INFO_DATA.TABLEVAR.name[unsort]; 
    z      = INFO_DATA.TABLEVAR.zhd[unsort];      
    printf("\t -> %2d with SNID=%8.8s at z=%.5f \n", 
	   NTMP, snid, z );	   
  }

  fflush(stdout);

  if ( NDUPL_SET >= MXSTORE ) {
    sprintf(c1err,"NDUPL_SET=%d exceeds bound, MXSTORE_DUPLICATE=%d",
	    NDUPL_SET, MXSTORE );
    sprintf(c2err,"Check duplicates");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

 
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
      { merge_duplicates(NDUPL_LIST[idup], UNSORT_DUPL[idup] ); }
  }
  else {
    sprintf(c1err,"Invalid iflag_duplicate=%d", INPUTS.iflag_duplicate );
    sprintf(c2err,"grep IFLAG_DUPLICATE  SALT2mu.c");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


 DONE:
  free(zList);  free(unsortList);

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
  //
  //   *** WARNING: need to refactor ****

  int i, isn, ipar, ipar2, ISN_SAVE ;
  double fitpar[NLCPAR], fiterr[NLCPAR];
  double COVFIT_INV[NLCPAR][NLCPAR], COVINT[NLCPAR][NLCPAR];
  double wgt, sqerr, wgtsum[NLCPAR], covfit, covtot, fitpar_tmp  ;
  char stringList_fitparOrig[NLCPAR][100], *name ;
  //  char fnam[] = "merge_duplicates" ;

  // ----------- BEGIN ----------------
  
  ISN_SAVE = isnList[0] ;
  name = INFO_DATA.TABLEVAR.name[ISN_SAVE] ;

  printf("# ------ merge %d duplicates for %s ---------- \n", NDUPL, name);
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
    if ( i>0 ) { setbit_CUTMASK(isn, CUTBIT_DUPL, &INFO_DATA.TABLEVAR); }

    for(ipar=0; ipar < NLCPAR; ipar++ ) {     // mB,x1,c
      sqerr = INFO_DATA.TABLEVAR.covmat_tot[isn][ipar][ipar] ;
      wgt    = 1.0/sqerr;
      fitpar_tmp = INFO_DATA.TABLEVAR.fitpar[ipar][isn];
      wgtsum[ipar] += wgt ;
      fitpar[ipar] += wgt * fitpar_tmp;

      sprintf(stringList_fitparOrig[ipar],"%s %6.3f",
	      stringList_fitparOrig[ipar], fitpar_tmp );

      for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) {
	covfit = INFO_DATA.TABLEVAR.covmat_fit[isn][ipar][ipar2];
	COVFIT_INV[ipar][ipar2] += 1.0/covfit ;
      }
    } // loop over mB,x1,c
  }

  get_COVINT_model(-1,COVINT);

  for(ipar=0; ipar < NLCPAR; ipar++ ) { 
    fitpar[ipar] /= wgtsum[ipar] ;
    fiterr[ipar]  = sqrt(1.0/wgtsum[ipar]);
    INFO_DATA.TABLEVAR.fitpar[ipar][ISN_SAVE] = fitpar[ipar] ;
    for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) { 
      covfit = 1.0 / COVFIT_INV[ipar][ipar2] ; 
      covtot = covfit + COVINT[ipar][ipar2] ; 
      INFO_DATA.TABLEVAR.covmat_fit[ISN_SAVE][ipar][ipar2] = covfit ;
      INFO_DATA.TABLEVAR.covmat_tot[ISN_SAVE][ipar][ipar2] = covtot ;
    }
  }


  // print old & new SALT2 params
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    printf("  %s-%2.2s(%s) -> %6.3f \n",  
	   name, BIASCOR_NAME_LCFIT[ipar], 
	   stringList_fitparOrig[ipar], fitpar[ipar] );
  }

  // print new COV
  printf("  %s-COV(no sigInt | +sigInt): \n", name );
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    printf("\t");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ )  { 
      covfit = INFO_DATA.TABLEVAR.covmat_fit[ISN_SAVE][ipar][ipar2] ;
      printf("%10.3le ", covfit);
    }
    printf("  |  ");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ ) {
      covtot = INFO_DATA.TABLEVAR.covmat_tot[ISN_SAVE][ipar][ipar2] ;
      printf("%10.3le ", covtot );
    }

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
  // Sep 7 2019: STOP if INPUTS.fixpar_all is set.

  double redchi2, covParam ;
  double step1 = INPUTS.covint_param_step1 ;
  int STOP_TOL, STOP_MXFIT, STOP_COV0, retCode, USE_CCPRIOR ;
  int NFIT_ITER = FITRESULT.NFIT_ITER ;
  char msg[100];
  //  char fnam[] = "prepNextFit" ;

  // ----------------- BEGIN -------------

  retCode = FITFLAG_DONE ;
  FITINP.COVINT_PARAM_LAST = FITINP.COVINT_PARAM_FIX ;
  USE_CCPRIOR = INFO_CCPRIOR.USE; 

  // check reasons to stop all fitting

  if ( INPUTS.fitflag_sigmb == 0 && !USE_CCPRIOR  )  { 
    return(FITFLAG_DONE); 
  }   // fix sigint --> no more fits

  // ----------------------------------

  redchi2 = FITRESULT.CHI2RED_1A ; 

  // check reasons to stop fitting
  STOP_TOL    = ( fabs(redchi2-1.0) < INPUTS.redchi2_tol ) ;
  STOP_MXFIT  = ( NFIT_ITER == MAXFITITER-1 || INPUTS.fixpar_all ) ;
  STOP_COV0   = ( NFIT_ITER > 0 && FITINP.COVINT_PARAM_FIX == 0.0)  ;
  
  // for CC prior, require at least 2 iterations
  if ( USE_CCPRIOR > 0 && NFIT_ITER == 0 ) { STOP_TOL = 0 ; }

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

  else if ( INPUTS.fitflag_sigmb == 0 && USE_CCPRIOR > 0 ) {
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
  // function to be minimized by MINUIT
  //c flat=1 read input, flag 2=gradient, flag=3 is final value
  double M0, alpha, beta, gamma, alpha0, beta0, gamma0;
  double da_dz, db_dz, dg_dz, da_dm, db_dm ;
  double scalePCC, scalePIa, scalePROB_fitpar, nsnfitIa, nsnfitcc;
  int    NSN_DATA, n, nsnfit, nsnfit_truecc, idsample, cutmask ;
  int    NDIM_BIASCOR; 
  int    DUMPFLAG=0, dumpFlag_muerrsq=0 ;
  double chi2sum_tot, mures, sqmures;
  double muerrsq, muerrsq_last, muerrsq_raw, muerrsq_tmp, sqsigCC=0.001 ;
  double chi2evt_Ia, chi2sum_Ia, chi2evt, sigCC_chi2penalty=0.0 ;
  double mu, mb, x1, c, z, zerr, logmass, dl, mumodel, mumodel_store;
  double muBias, muBiasErr, muBias_zinterp, muCOVscale, gammaDM ;
  double muerr, muerr_raw, muerr_last;
  int    ipar,  ipar2, ia, ib, ig, INTERPFLAG_abg ;
  int    sntype=0, idsurvey=0, SIM_NONIA_INDEX, IS_SIM ;
  double omega_l, omega_k, wde, wa, cosPar[NCOSPAR] ;
  double *hostPar, logmass_cen, logmass_tau ;
  double ProbRatio_Ia, ProbRatio_CC ;
  double covmat_tot[NLCPAR][NLCPAR], covmat_fit[NLCPAR][NLCPAR] ;
  char   *name ;
  MUZMAP_DEF  *CCPRIOR_MUZMAP ;
  int  USE_CCPRIOR=0, USE_CCPRIOR_H11=0 ;
  BIASCORLIST_DEF     BIASCORLIST ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
  FITPARBIAS_DEF FITPARBIAS_ALPHABETA[MXa][MXb][MXg]; // bias at each a,b
  double   MUCOVSCALE_ALPHABETA[MXa][MXb][MXg]; // (I) muCOVscale at each a,b
  double   *fitParBias;

  double Prob_Ia, Prob_CC, Prob_SUM, dPdmu_Ia, dPdmu_CC ;
  double PTOTRAW_Ia=0.0, PTOTRAW_CC, PTOT_Ia, PTOT_CC, PSUM ;
  double muerrsq_update, muerr_update ;

  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;  

  char fnam[]= "fcn";

  // --------------- BEGIN -------------

  FITRESULT.NCALL_FCN++ ;

  // Nov 2019, bail on inf or nan.
  for(ipar=1; ipar<=5; ipar++ ) {
    if ( isnan(xval[ipar]) ) { *fval = 1.0E14; return; }
    if ( isinf(xval[ipar]) ) { *fval = 1.0E14; return; }
  }

  //Set input cosmology parameters

  alpha0       = xval[1] ;
  beta0        = xval[2] ;
  da_dz        = xval[3] ;
  db_dz        = xval[4] ;
  gamma0       = xval[5] ;
  dg_dz        = xval[6] ;
  logmass_cen  = xval[7] ; 
  logmass_tau  = xval[8] ;
  omega_l      = xval[9] ;
  omega_k      = xval[10] ;
  wde          = xval[11] ;
  wa           = xval[12] ;
  scalePROB_fitpar  = xval[13] ;
  hostPar      = &xval[IPAR_GAMMA0];
  da_dm        = xval[15];  // added Apr 2 2018
  db_dm        = xval[16];  // idem

  NDIM_BIASCOR = INFO_BIASCOR.NDIM;


  // load cosPar array to pass to functions below (RK Jan 2016)
  cosPar[0] = omega_l ;
  cosPar[1] = omega_k ;
  cosPar[2] = wde ;
  cosPar[3] = wa ;

  /*
  // xxxxxxxxx
  if ( (ncall_fcn > 1000 && ncall_fcn < 1065) || ncall_fcn<10 ) {
    printf(" xxx ncall=%4d: a=%f b=%f S_CC=%f  sigint=%f\n", 
	   ncall_fcn, xval[1], xval[2], xval[13], xval[14]); 
    fflush(stdout);
  }
  // xxxxxxxxx
  */

  // -------------------------------

  NSN_DATA         = INFO_DATA.TABLEVAR.NSN_ALL ;
  USE_CCPRIOR      = INFO_CCPRIOR.USE; 
  USE_CCPRIOR_H11  = INFO_CCPRIOR.USEH11; 
  CCPRIOR_MUZMAP   = &INFO_CCPRIOR.MUZMAP;
  

  if ( USE_CCPRIOR  ) {   
    // load CCPRIOR_MUZMAP
    fcn_ccprior_muzmap(xval, USE_CCPRIOR_H11, CCPRIOR_MUZMAP);   
    ProbRatio_Ia = ProbRatio_CC = 0.0 ;
  }


  // -------------------------------
  // For biasCor, get INTERP weights for alpha and beta grid.
  // Beware that this is not quite right for z-dependent alpha,beta,
  // but it's faster to compute ia,ib here outside the data loop.
  INTERPFLAG_abg = 0;
  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_5D) ||
       (INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT)  ) {
    INTERPFLAG_abg = 1; 
    // check for z-dependent alpha or beta
    if ( INPUTS.ipar[3] || INPUTS.ipar[4] ) { INTERPFLAG_abg = 2; } 
    // check for hostmass-dependent alpha or beta (April 2 2018)
    if ( INPUTS.ipar[15] || INPUTS.ipar[16] ) { INTERPFLAG_abg = 2; } 
  }



  // -------------------------------
  chi2sum_tot = chi2sum_Ia = 0.0;
  nsnfit      = nsnfit_truecc = 0 ;
  nsnfitIa = nsnfitcc = 0.0 ;
  zero_contam_CCprior(&CONTAM_MURES_BINS);
  zero_contam_CCprior(&CONTAM_REDSHIFT_BINS);

  FITRESULT.NSNFIT = 0;
  FITRESULT.NSNFIT_TRUECC = 0;

  for (n=0; n< NSN_DATA; ++n)  {

    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n] ; 

    if ( cutmask ) { continue; }

    // - - - - -

    INFO_DATA.mures[n]     = -999. ;
    INFO_DATA.mupull[n]    = -999. ;
    INFO_DATA.mu[n]        = -999. ;
    INFO_DATA.muerr[n]     = -999. ;    
    INFO_DATA.muerr_raw[n] = -999. ;  // no scale and no sigInt
    name     = INFO_DATA.TABLEVAR.name[n] ;
    idsample = (int)INFO_DATA.TABLEVAR.IDSAMPLE[n] ;
    z        = (double)INFO_DATA.TABLEVAR.zhd[n] ;     
    logmass  = (double)INFO_DATA.TABLEVAR.logmass[n];
    zerr     = (double)INFO_DATA.TABLEVAR.zhderr[n] ;
    mb       = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_mB][n] ;
    x1       = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_x1][n] ;
    c        = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_c][n] ;
    mumodel_store   = (double)INFO_DATA.TABLEVAR.mumodel[n] ; 

    SIM_NONIA_INDEX = (int)INFO_DATA.TABLEVAR.SIM_NONIA_INDEX[n];
    IS_SIM          = (INFO_DATA.TABLEVAR.IS_SIM == true);
    
    if ( USE_CCPRIOR ) { 
      PTOTRAW_Ia  = (double)INFO_DATA.TABLEVAR.pIa[n] ; 
      sntype      = (int)INFO_DATA.TABLEVAR.SNTYPE[n] ; 
      idsurvey    = (int)INFO_DATA.TABLEVAR.IDSURVEY[n];
    }
    muBias_zinterp = INFO_DATA.muBias_zinterp[n] ; 
    muerrsq_last   = INFO_DATA.muerrsq_last[n] ;
    muerr_last     = INFO_DATA.muerr_last[n] ;

    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) 
	{ covmat_tot[ipar][ipar2] = 
	    INFO_DATA.TABLEVAR.covmat_tot[n][ipar][ipar2] ; 
	}
    }
            
    for(ia=0; ia<MXa; ia++ ) {
      for(ib=0; ib<MXb; ib++ ) {
	for(ig=0; ig<MXg; ig++ ) {
	  MUCOVSCALE_ALPHABETA[ia][ib][ig] = 
	    INFO_DATA.MUCOVSCALE_ALPHABETA[n][ia][ib][ig] ; 
	  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {
	    FITPARBIAS_ALPHABETA[ia][ib][ig].VAL[ipar] = 
	      INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig].VAL[ipar] ; 
	    FITPARBIAS_ALPHABETA[ia][ib][ig].ERR[ipar] = 
	      INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig].ERR[ipar] ; 
	  }
	}
      }
    }
    fitParBias = INFO_DATA.fitParBias[n] ; 
      
    if ( z < 1.0E-8) { continue ; } // Jun 3 2013 (obsolete?)

    // fetch alpha,beta,gamma (include z-dependence)
    fcnFetch_AlphaBetaGamma(xval, z, logmass, &alpha, &beta, &gamma); 

    gammaDM = get_gammadm_host(z, logmass, hostPar );

    //    DUMPFLAG = ( strcmp(name,"93018")==0 ) ; // xxx REMOVE 44
    if ( INTERPFLAG_abg ) {
      get_INTERPWGT_abg(alpha, beta, gammaDM, DUMPFLAG, &INTERPWGT, name);
    }

    //    DUMPFLAG = 0 ;

    // get mag offset for this z-bin
    M0    = fcn_M0(n, &xval[MXCOSPAR] );

    // compute distance modulus from cosmology params
    if ( INPUTS.FLOAT_COSPAR ) {
      dl       = cosmodl_forFit(z, cosPar) ;
      mumodel  = 5.0*log10(dl) + 25.0 ;
    }
    else 
      { mumodel = mumodel_store ; }


    // compute error-squared on distance mod
    muerrsq  = fcn_muerrsq(name, alpha, beta, gamma, covmat_tot,
			   z, zerr, 0);

    // --------------------------------
    // Compute bias from biasCor sample
    BIASCORLIST.z           = z;
    BIASCORLIST.logmass     = logmass;
    BIASCORLIST.alpha       = alpha ;
    BIASCORLIST.beta        = beta ;
    BIASCORLIST.gammadm     = gammaDM ;
    BIASCORLIST.idsample    = idsample;
    BIASCORLIST.FITPAR[INDEX_mB] = mb ;
    BIASCORLIST.FITPAR[INDEX_x1] = x1 ;
    BIASCORLIST.FITPAR[INDEX_c]  = c ;
    muBias = muBiasErr = 0.0 ;  muCOVscale=1.0 ; 

    if ( NDIM_BIASCOR >= 5 ) {
      get_muBias(name, &BIASCORLIST,      // (I) misc inputs
		 FITPARBIAS_ALPHABETA,    // (I) bias at each a,b,g
		 MUCOVSCALE_ALPHABETA,    // (I) muCOVscale at each a,b
		 &INTERPWGT,              // (I) wgt at each a,b,g grid point
		 fitParBias,     // (O) interp bias on mB,x1,c
		 &muBias,        // (O) interp bias on mu
		 &muBiasErr,     // (O) stat-error on above
		 &muCOVscale );  // (O) scale bias on muCOV     
    }
    else if ( NDIM_BIASCOR == 1 ) {
      muBias  = muBias_zinterp ; 
    }
       
    if ( n == -1 ) { 
      printf(" xxx %s-%3.3d %s  a,b,g=%.4f,%.4f,%.4f  "
	     "muBias=%7.4f  (fpb0=%.4f,%.4f) \n",
	     fnam, FITRESULT.NCALL_FCN, name, alpha,beta,gamma, 
	     muBias,
	     FITPARBIAS_ALPHABETA[0][0][0].VAL[ILCPAR_MAX],
	     FITPARBIAS_ALPHABETA[0][1][0].VAL[ILCPAR_MAX]
	     );
      fflush(stdout);
    }

    // load muBias info into globals

    INFO_DATA.muBias[n]     = muBias ;
    INFO_DATA.muBiasErr[n]  = muBiasErr ;
    INFO_DATA.muCOVscale[n] = muCOVscale ;

    // zero out muBiasErr after storing it, since adding this
    // would contradict the muCOVscale correction.
    muBiasErr = 0.0 ; 
                   
    // apply bias corrections
    muerrsq  *= muCOVscale ;  // error scale 
    muerrsq  += ( muBiasErr*muBiasErr); 
    muerr     = sqrt(muerrsq);	
    // ------------------------

    mu       = mb  + alpha*x1 - beta*c - gammaDM ; 

    mu      -= muBias ;      // bias correction 
    mures    = mu - M0 - mumodel ;
    sqmures  = mures*mures ;

    // store info
    INFO_DATA.mu[n]       = mu ;
    INFO_DATA.mures[n]    = mures ;
    INFO_DATA.mupull[n]   = mures/muerr ;
    INFO_DATA.mumodel[n]  = mumodel ;
    INFO_DATA.muerr[n]    = muerr ;
    INFO_DATA.M0[n]       = M0 ;

    if ( !USE_CCPRIOR  ) {
      // original SALT2mu chi2 with only spec-confirmed SNIa 
      nsnfit++ ;        nsnfitIa  = (double)nsnfit ;
      chi2evt_Ia    = sqmures/muerrsq ;
      chi2sum_Ia   += chi2evt_Ia ;
      chi2evt       = chi2evt_Ia ;

      // check option to add log(sigma) term for 5D biasCor
      if ( INPUTS.fitflag_sigmb == 2 ) 
      	{ chi2evt  += log(muerrsq/muerrsq_last); } 
    } 


    if ( USE_CCPRIOR  ) {
      // BEAMS-like chi2 = -2ln [ PIa + PCC ]
      DUMPFLAG = (n == -44);
      nsnfit++ ;

      if ( INPUTS.ipar[IPAR_scalePCC] <= 1 ) {
	scalePCC    = scalePROB_fitpar ;
	PTOTRAW_CC  = 1.0 - PTOTRAW_Ia ;  // CC prob 
	PSUM        = PTOTRAW_Ia + scalePCC*PTOTRAW_CC ;
	PTOT_CC     = scalePCC * PTOTRAW_CC/PSUM ;
	PTOT_Ia     = PTOTRAW_Ia/PSUM ;
      }
      else {
	// u13=2 --> Use H11, Eq4  Bayes factor (Sep 30 2019)
	scalePIa = scalePROB_fitpar;
	PTOT_Ia  = (PTOTRAW_Ia * scalePIa)/(1.0 + (scalePIa-1.0)*PTOTRAW_Ia);
	PTOT_CC  = 1.0 - PTOT_Ia ;
      }

      /* xxxx mark delete May 8 2020 xxxxxxxxxxx
      if ( force_probcc0(sntype,idsurvey) ) 
	{ PTOT_Ia = 1.0;  PTOT_CC = 0.0 ;  } // spec-confirmed SNIa
      xxxxxxx end mark xxxxx */

      if ( INPUTS.fitflag_sigmb == 2 ) 
	{ muerrsq_update = muerrsq ; }  // current muerr for 5D biasCor
      else
	{ muerrsq_update = muerrsq_last ; } //  fixed muerr for 1D biasCor
      muerr_update   = sqrt(muerrsq_update);

      chi2evt_Ia = sqmures/( muerrsq - muBiasErr*muBiasErr)  ;
      dPdmu_Ia   = ( exp(-0.5*chi2evt_Ia) ) * (PIFAC/muerr_update) ; 
      Prob_Ia    = PTOT_Ia * dPdmu_Ia ;

      if ( USE_CCPRIOR_H11 ) {
	dPdmu_CC = prob_CCprior_H11(n, mures, &xval[IPAR_H11], 
				    &sqsigCC, &sigCC_chi2penalty );
      }
      else { // CC prob from sim
	dPdmu_CC = prob_CCprior_sim(idsample, CCPRIOR_MUZMAP, 
				    z, mures, DUMPFLAG );
      }
      Prob_CC   = PTOT_CC * dPdmu_CC ;

      Prob_SUM  = Prob_Ia + Prob_CC ; // BEAMS prob in Eq. 6 of BBC paper

      // xxxxxxxxxxxxxxxxxxx
      if ( DUMPFLAG ) {
	printf(" xxx --------------------------- \n");
	printf(" xxx fcn dump for SN = '%s' \n", name);
	printf(" xxx dPdmu_CC=%le \n", dPdmu_CC);
	printf(" xxx PTOT_CC = scale*PTOTRAW/PSUM  = %le*%le/%le = %le\n",
	       scalePCC, PTOTRAW_CC, PSUM, PTOT_CC );
	printf(" xxx Prob(Ia,CC) = %le, %le \n", Prob_Ia, Prob_CC);
	debugexit(fnam);
      }
      // xxxxxxxxxxxxxxxxxxx

      
      // sum total chi2 that includes Ia + CC
      if ( Prob_SUM > 0.0 ) {
	ProbRatio_Ia = Prob_Ia / Prob_SUM ;
	ProbRatio_CC = Prob_CC / Prob_SUM ;
	nsnfitIa   +=  ProbRatio_Ia ; 
	nsnfitcc   +=  ProbRatio_CC ; 
	chi2sum_Ia += (ProbRatio_Ia * chi2evt_Ia) ; 
	
	if ( *iflag == 3 ) {
	  INFO_DATA.probcc_beams[n] = ProbRatio_CC;
	  sum_contam_CCprior(&CONTAM_MURES_BINS, ProbRatio_Ia, mures,
			     SIM_NONIA_INDEX); 
	  sum_contam_CCprior(&CONTAM_REDSHIFT_BINS, ProbRatio_Ia, z,
			     SIM_NONIA_INDEX); 
	}

	Prob_SUM    *= (0.15/PIFAC)  ;  
	chi2evt      = -2.0*log(Prob_SUM);
	//  chi2evt += sigCC_chi2penalty ; // prevent sigCC<0 (7.17.2018)
      }
      else {
	chi2evt      = 1.0E8 ;
      }

    } // end CCprios.USE if block

    chi2sum_tot      += chi2evt;
    INFO_DATA.chi2[n] = chi2evt; // store each chi2 to allow for outlier cut

    // check things on final pass
    if (  *iflag==3 ) {	

      // Jan 26 2018: store raw muerr without intrinsic cov
      for(ipar=0; ipar < NLCPAR; ipar++ ) {
	for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) {	  
	  covmat_fit[ipar][ipar2] =
	    (double)INFO_DATA.TABLEVAR.covmat_fit[n][ipar][ipar2];
	}
      }
      muerrsq_raw = fcn_muerrsq(name,alpha,beta,gamma,covmat_fit, z,zerr, 0);
      muerr_raw   = sqrt(muerrsq_raw);
      INFO_DATA.muerr_raw[n] = muerr_raw;   

      // check user dump with total error (Jun 19 2018)
      dumpFlag_muerrsq = ( strcmp(name,INPUTS.SNID_MUCOVDUMP) == 0 );
      if ( dumpFlag_muerrsq ) {
	muerrsq_tmp = fcn_muerrsq(name,alpha,beta,gamma,covmat_fit,z,zerr,1) ;
	muerrsq_tmp = fcn_muerrsq(name,alpha,beta,gamma,covmat_tot,z,zerr,1) ;
      }
      
      if ( IS_SIM && SIM_NONIA_INDEX > 0 ) { nsnfit_truecc++ ; } 

      // store reference errors for 1/sigma term
      if ( DOFIT_FLAG == FITFLAG_CHI2 ) {
	INFO_DATA.muerr_last[n]   = muerr;
	INFO_DATA.muerrsq_last[n] = muerrsq;
	INFO_DATA.sigCC_last[n]   = sqrt(sqsigCC) ;
	INFO_DATA.sqsigCC_last[n] = sqsigCC ;
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
    FITRESULT.CHI2SUM_1A = chi2sum_Ia ;
    FITRESULT.CHI2RED_1A = chi2sum_Ia/(double)(nsnfitIa-FITINP.NFITPAR_FLOAT) ;
    FITRESULT.NSNFIT_1A  = nsnfitIa ; 
    FITRESULT.NSNFIT_CC  = nsnfitcc ;  // 9.24.2019
    FITRESULT.ALPHA      = alpha ;
    FITRESULT.BETA       = beta ;
    FITRESULT.GAMMA      = gamma ;
  }

  *fval = chi2sum_tot;

  return;

}    // end of fcn

// ==============================================
void  fcnFetch_AlphaBetaGamma(double *xval, double z, double logmass, 
			      double *alpha, double *beta, double *gamma) {

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
  double g0          = xval[5] ; // gamma
  double dg_dz       = xval[6] ;
  double aHost       = xval[15]; // dalpha/dlog(Mhost)
  double bHost       = xval[16]; // dbeta/dlog(Mhost)
  double logmass_cen = xval[7];  // log(Msplit) for gamma0
  double dlogmass    = logmass - logmass_cen ;
  double alpha_local, beta_local, gamma_local ;
  int    OPT_LOGMASS_SLOPE=0, OPT_LOGMASS_SPLIT=0 ;

  alpha_local = a0  + z*da_dz ;
  beta_local  = b0  + z*db_dz ;  
  gamma_local = g0  + z*dg_dz ;

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
  *gamma = gamma_local ;

  return ;

} // end fcnFetch_AlphaBetaGamma

// ================================
double fcn_M0(int n, double *M0LIST) {

  // return model M0 for this data index 'n'
  // and list of M0LIST in each z bin
  // Jun 27 2017: REFACTOR z bins
  // Jan 29 2019: if no iz1 bin, return(M0) instead of retrn(M0bin0)

  int LDMP=0;
  int iz, iz0, iz1, NBINz, NSN_BIASCOR ;
  double M0, zdata, zbin0, zbin1, zfrac, M0bin0, M0bin1, *ptr_zM0 ;  
  char fnam[] = "fcn_M0";

  // ----------- BEGIN ----------

  M0 = M0_DEFAULT ;
  iz0     = INFO_DATA.TABLEVAR.IZBIN[n];
  zdata   = INFO_DATA.TABLEVAR.zhd[n];
  ptr_zM0 = INFO_BIASCOR.zM0 ;
  NSN_BIASCOR = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  

  if ( INPUTS.uM0 == M0FITFLAG_CONSTANT ||
       INPUTS.uM0 == M0FITFLAG_ZBINS_FLAT ) {
    iz    = iz0;
    M0    = M0LIST[iz] ;
  }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_INTERP ) {
    // linear interp
    if ( NSN_BIASCOR <= 0 ) {
      sprintf(c1err,"Cannot implement option uM0=%d", INPUTS.uM0 );
      sprintf(c2err,"without simFile_bias");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
    }
    
    NBINz   = INPUTS.BININFO_z.nbin ;
    zbin0   = ptr_zM0[iz0]; // wgt z-avg in this z bin

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

    zbin1  = ptr_zM0[iz1];
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
void get_INTERPWGT_abg(double alpha, double beta, double gammadm, int DUMPFLAG,
		      INTERPWGT_AlphaBetaGammaDM *INTERPWGT, char *callFun ) {

  // Created April 26 2016 by R.Kessler
  //
  // Inputs:  alpha, beta, gammadm, DUMPFLAG
  //   (beware that gammadm is the mag offset, NOT gamma)
  //
  // Ouput struct *INTERPWGT with WGT at each alpha,beta corner
  //     used for interpolation in get_muBias.
  //   
  //  This function is called from fcn before the loop over
  //  data in order to reduce the CPU overhead for data event.
  //
  // Sep 4 2017: init INTERPWGT->WGT[ia][ib] = 0.0 
  // Dec 21 2017: widen min,max bound to allow extrapolation

  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM ;
  int    IA, IB, IG, ia, ib, ig, NBINa, NBINb, NBINg ;
  double a_binSize, b_binSize, g_binSize ;
  double a_bound_min, a_bound_max, b_bound_min, b_bound_max ;
  double g_bound_min, g_bound_max ;
  double a_interp, b_interp, g_interp ;
  char fnam[] = "get_INTERPWGT_abg" ;

  // ------------------ BEGIN ---------------

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;

  NBINa     = (*BININFO_SIM_ALPHA).nbin ;
  NBINb     = (*BININFO_SIM_BETA).nbin ;
  NBINg     = (*BININFO_SIM_GAMMADM).nbin ;
  a_binSize = (*BININFO_SIM_ALPHA).binSize ;
  b_binSize = (*BININFO_SIM_BETA).binSize ;
  g_binSize = (*BININFO_SIM_GAMMADM).binSize ;

  a_bound_min = (*BININFO_SIM_ALPHA).avg[0];
  a_bound_max = (*BININFO_SIM_ALPHA).avg[NBINa-1] ;
  b_bound_min = (*BININFO_SIM_BETA).avg[0];
  b_bound_max = (*BININFO_SIM_BETA).avg[NBINb-1] ;
  g_bound_min = (*BININFO_SIM_GAMMADM).avg[0];
  g_bound_max = (*BININFO_SIM_GAMMADM).avg[NBINg-1] ;

  IA = IBINFUN(alpha,   BININFO_SIM_ALPHA,   2, fnam );
  IB = IBINFUN(beta,    BININFO_SIM_BETA,    2, fnam );
  IG = IBINFUN(gammadm, BININFO_SIM_GAMMADM, 2, fnam );

  // get local "interp" values of alpha&beta which do not
  // extend past the storage grid to avoid crazy extrapolations.
  a_interp = alpha ;
  if ( a_interp < a_bound_min ) { a_interp = a_bound_min; }
  if ( a_interp > a_bound_max ) { a_interp = a_bound_max; }

  b_interp = beta ;  
  if ( b_interp < b_bound_min ) { b_interp = b_bound_min; }
  if ( b_interp > b_bound_max ) { b_interp = b_bound_max; }

  g_interp = gammadm ;  
  if ( g_interp < g_bound_min ) { g_interp = g_bound_min; }
  if ( g_interp > g_bound_max ) { g_interp = g_bound_max; }

  int FIRST = 1;
  double a_grid, b_grid, g_grid, Da, Db, Dg ;
  double DUMPLIST_Da[MXa][MXb][MXg], DUMPLIST_Db[MXa][MXb][MXg];
  double DUMPLIST_Dg[MXa][MXb][MXg];
  double SUMWGT = 0.0 ;
  double WGT[MXa][MXb][MXg];

  /* xxxxxx
  printf(" xxx  ------------ DUMP --------------- \n", fnam);
  printf(" xxx NBIN[a,b,g] = %d %d %d    IA,IB,IG=%d,%d,%d\n", 
	 NBINa, NBINb, NBINg, IA,IB,IG );
  printf(" xxx [a,b,g]_interp  = %f %f %f \n", 
	 a_interp, b_interp, g_interp);
  printf(" xxx [a,b,g]_min     = %f %f %f \n", 
	 a_bound_min, b_bound_min, g_bound_min);
  printf(" xxx [a,b,g]_max     = %f %f %f \n", 
	 a_bound_max, b_bound_max, g_bound_max);
  printf(" xxx [a,b,g]_binSize = %f %f %f \n", 
	 a_binSize, b_binSize, g_binSize );
  printf(" xxx [a,b,g] input: %f %f %f \n",
	 alpha, beta, gammadm );
  fflush(stdout);
  xxxxxx  */


  // init WGT map
  for(ia=0; ia < MXa; ia++ ) {
    for(ib=0; ib < MXb; ib++ ) {
      for(ig=0; ig < MXg; ig++ ) {
	WGT[ia][ib][ig] = 0.0 ;
	INTERPWGT->WGT[ia][ib][ig] = 0.0 ;
      }
    }
  }

  for(ia=IA-1; ia <= IA+1; ia++ ) {
    if ( ia < 0      ) { continue ; }
    if ( ia >= NBINa ) { continue ; }
    
    for(ib=IB-1; ib <= IB+1; ib++ ) {
      if ( ib < 0      ) { continue ; }
      if ( ib >= NBINb ) { continue ; }

      for(ig=IG-1; ig <= IG+1; ig++ ) {
	if ( ig < 0      ) { continue ; }
	if ( ig >= NBINg ) { continue ; }

	a_grid = (*BININFO_SIM_ALPHA).avg[ia] ;
	b_grid = (*BININFO_SIM_BETA).avg[ib] ;
	g_grid = (*BININFO_SIM_GAMMADM).avg[ig] ;

	Da     = fabs(a_interp  - a_grid) ;
	Db     = fabs(b_interp  - b_grid) ;
	Dg     = fabs(g_interp  - g_grid) ;

	if ( NBINa > 1 ) { Da /= a_binSize; } else { Da = 0.0 ; }
	if ( NBINb > 1 ) { Db /= b_binSize; } else { Db = 0.0 ; }
	if ( NBINg > 1 ) { Dg /= g_binSize; } else { Dg = 0.0 ; }

	if ( DUMPFLAG ) {
	  DUMPLIST_Da[ia][ib][ig] = Da;
	  DUMPLIST_Db[ia][ib][ig] = Db;
	  DUMPLIST_Dg[ia][ib][ig] = Dg;
	}

	if ( fabs(Da) > 0.99999 ) { continue ; }
	if ( fabs(Db) > 0.99999 ) { continue ; }
	if ( fabs(Dg) > 0.99999 ) { continue ; }
      
	WGT[ia][ib][ig] = (1.0 - Da) * (1.0 - Db) * (1.0 - Dg) ;
	SUMWGT     += WGT[ia][ib][ig] ;
	
	if ( FIRST ) { 
	  INTERPWGT->ia_min = ia ;
	  INTERPWGT->ib_min = ib ; 
	  INTERPWGT->ig_min = ig ;  
	}
	INTERPWGT->ia_max = ia ;
	INTERPWGT->ib_max = ib ;
	INTERPWGT->ig_max = ig ;
	
	FIRST = 0 ;

      } // end ig
    } // end ib
  } // end ia


  // ------------------------------
  if ( SUMWGT < 1.0E-9 ) {
    print_preAbort_banner(fnam);
    printf("   Alpha  =%.3f   Alpha_interp=%.3f  (bnd: %.3f to %.3f)\n",  
	   alpha, a_interp, a_bound_min, a_bound_max);
    printf("   Beta   =%.3f   Beta_interp=%.3f   (bnd: %.3f to %.3f)\n",
	   beta,    b_interp, b_bound_min, b_bound_max );
    printf("   Gammadm=%.3f   Gammadm_interp=%.3f (bnd: %.3f to %.3f)\n", 
	   gammadm, g_interp, g_bound_min, g_bound_max );

    for(ia=0; ia < MXa; ia++ ) {
      for(ib=0; ib < MXb; ib++ ) {
	for(ig=0; ig < MXg; ig++ ) {
	printf("\t WGT[ia,ib,ig=%d,%d,%d] = %f \n", 
	       ia,ib,ig, WGT[ia][ib][ig] );
	}
      }
    }

    sprintf(c1err,"Invalid SUMWGT=%f  (called from %s)", SUMWGT, callFun);
    sprintf(c2err,"IA=%d  IB=%d  IG=%d", IA, IB, IG);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  // divide WGT by SUMWGT so that \Sum WGT = 1
  int NEGWGT=0;
  for(ia=INTERPWGT->ia_min; ia<=INTERPWGT->ia_max; ia++ ) {
    for(ib=INTERPWGT->ib_min; ib<=INTERPWGT->ib_max; ib++ ) {
      for(ig=INTERPWGT->ig_min; ig<=INTERPWGT->ig_max; ig++ ) {
	INTERPWGT->WGT[ia][ib][ig] = WGT[ia][ib][ig]/SUMWGT ;
	if ( INTERPWGT->WGT[ia][ib][ig] < 0.0 ) { NEGWGT = 1; }
      }
    }
  }

  int LDMP = DUMPFLAG ;
  if ( LDMP || NEGWGT ) {
    int ia_min = INTERPWGT->ia_min ;    int ia_max = INTERPWGT->ia_max ;
    int ib_min = INTERPWGT->ib_min ;    int ib_max = INTERPWGT->ib_max ;
    int ig_min = INTERPWGT->ig_min ;    int ig_max = INTERPWGT->ig_max ;
    printf("xxx ------------------ [%s] ------------------------------ \n", 
	   callFun);
    printf("xxx alpha/alpha_interp = %f/%f \n", alpha,   a_interp);
    printf("xxx beta /beta_interp  = %f/%f \n", beta,    b_interp);
    printf("xxx gDM/gDM_interp     = %f/%f \n", gammadm, g_interp);
    printf("xxx SUMWGT = %le \n", SUMWGT);
    printf("xxx ia=%d to %d  ib=%d to %d  ig=%d to %d\n",
	   ia_min,ia_max, ib_min,ib_max, ig_min,ig_max );
    printf("xxx binsize(a,b,g) = %.4f, %.4f, %.4f \n",
	   a_binSize, b_binSize, g_binSize);
    
    for(ia=ia_min; ia<=ia_max; ia++ ) {
      for(ib=ib_min; ib<=ib_max; ib++ ) {
	for(ig=ig_min; ig<=ig_max; ig++ ) {
	  INTERPWGT->WGT[ia][ib][ig] = WGT[ia][ib][ig]/SUMWGT ;
	  printf("xxx ia,ib,ig=%d,%d,%d: "
		 "WGT = %.4f  Da,Db,Dg=%.3f,%.3f,%.3f\n",
		 ia,ib,ig, INTERPWGT->WGT[ia][ib][ig],
		 DUMPLIST_Da[ia][ib][ig], DUMPLIST_Db[ia][ib][ig],
		 DUMPLIST_Dg[ia][ib][ig] );
	}
      }
    }
    
  } // end LDMP

  if ( NEGWGT ) {
    sprintf(c1err,"Invalid Negative weight (called from %s)", callFun);
    sprintf(c2err,"See dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  return ;
  
} // end get_INTERPWGT_abg


// ===========================================================
double fcn_muerrsq(char *name, double alpha, double beta, double gamma,
		   double (*COV)[NLCPAR], 
		   double z, double zerr, int dumpFlag) {

  // Created Feb 29 2016 by R.Kessler
  //
  // return mu-error for inputs as follows:
  //
  // *name      : name of SN (in case of error msg)
  // alpha,beta : SALT2 standardization parmas for stretch & color
  // gamma      : SALT2 standard param for logmass (NOT USED)
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
  // 
  //
  // OPT=1 --> low-z approx for empty universe
  // OPT=2 --> exact calc.

  // Note that OPT=1 & 2 get closer as z increases.
  //
  // Jan 9 2018: remove zpecerr argment.

  double zerrtot, muerr = 0.0 ;
  double FAC   = 5.0/LOGTEN ;  
  double *cosPar = &INPUTS.parval[IPAR_OL] ;
  //  char fnam[]  = "fcn_muerrz" ;

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
void fcn_ccprior_muzmap(double *xval, int USE_CCPRIOR_H11, 
			MUZMAP_DEF *MUZMAP ) {

  // Created Jun 2019     [code moved from fcn]
  //
  // Setup info needed to compute CC prior.
  // Do not include z-dependence for each event because it's
  // too slow to compute this inside the data loop.
  // Get DMUPDF in each z bin
  //
  // Inputs:
  //   xval[]          : array of parameters passed to fcn
  //   USE_CCPRIOR_H11 : H11 option for BEAMS
  //
  // Output:
  //   MUZMAP structure

  int NSAMPLE = NSAMPLE_BIASCOR ;
  int i, idsample ;
  double cosPar[10];
  //  char fnam[] = "fcn_ccprior_muzmap";

  // ----------- BEGIN ------------

  MUZMAP->alpha = xval[IPAR_ALPHA0];
  MUZMAP->beta  = xval[IPAR_BETA0];
  MUZMAP->M0    = INPUTS.nommag0 ;

  cosPar[0] = xval[IPAR_OL] ;
  cosPar[1] = xval[IPAR_Ok] ;
  cosPar[2] = xval[IPAR_w0] ;
  cosPar[3] = xval[IPAR_wa] ;
  for(i=0; i < NCOSPAR ; i++ ) { MUZMAP->cosPar[i] = cosPar[i] ; } 
  
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
      setup_DMUPDF_CCprior(idsample, 
			   &INFO_CCPRIOR.TABLEVAR_CUTS, MUZMAP );
  }

  // - - - - -

  if ( USE_CCPRIOR_H11 ) {
    for(i=0; i<NPAR_H11_TOT; i++)  { MUZMAP->H11_fitpar[i]=xval[IPAR_H11+i];}
  }

  return;

} // end fcn_ccprior_muzmap

// ================================================
void test_muerrz(void) {

  double z, muerr1, muerr2;
  double zerr = 0.001 ;
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
    print_preAbort_banner(fnam);
    printf("    z1*zpecerr = %f * %f = %f \n", z1, zpecerr, z1*zpecerr);
    printf("    zerr = %f \n", zerr);
    sprintf(c1err,"Cannot subtract zpecerr from zerr for SNID=%s", name );
    sprintf(c2err,"zerr=%f > (1+z)*zpecerr = (%f)*%f", 
	    zerr, 1+z, zpecerr);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
  }

  sqarg    = (sqzerr1 - sqzerr2 + sqzerr3 );

  if ( sqarg < 0.0 ) {
    print_preAbort_banner(fnam);
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


// ******************************************
void set_defaults(void) {

  int isurvey, order, ipar ;

  // ---------------- BEGIN -----------------
  
  set_EXIT_ERRCODE(EXIT_ERRCODE_SALT2mu);

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );

  INPUTS.cat_only   = false ;
  INPUTS.catfile_out[0] = 0 ;

  INPUTS.nfile_data = 0;
  sprintf(INPUTS.PREFIX,     "NONE" );

  INPUTS.opt_biasCor           = 0 ;
  INPUTS.sigint_biasCor        = -9.0 ; 
  INPUTS.snrmin_sigint_biasCor = BIASCOR_SNRMIN_SIGINT ;
  INPUTS.prescale_biasCor[0]   = 0;
  INPUTS.prescale_biasCor[1]   = 1;
  INPUTS.sigma_cell_biasCor    = 50.0; // flat WGT distribution in each cell
  INPUTS.nfile_biasCor         = 0 ;
  sprintf(INPUTS.fieldGroup_biasCor, "NONE" );
  INPUTS.use_fieldGroup_biasCor = 0 ;

  INPUTS.dumpflag_nobiasCor  = 0 ;
  INPUTS.frac_warn_nobiasCor = 0.02 ;

  sprintf(INPUTS.surveyGroup_biasCor, "NONE" );
  INPUTS.use_surveyGroup_biasCor = 0 ;

  sprintf(INPUTS.surveyList_noBiasCor, "NONE" );
  INPUTS.idsample_select[0] = 0 ;

  INPUTS.interp_biascor_logmass=1; // default is to do biasCor interp

  // default is to blind cosmo params for data
  INPUTS.blindFlag = BLINDMASK_FIXPAR; 

  INPUTS.NSPLITRAN      = 1; // default is all SN in one job
  INPUTS.JOBID_SPLITRAN = -9;

  INPUTS.iflag_duplicate = IFLAG_DUPLICATE_ABORT ;


  INPUTS.NCUTWIN = 0 ;
  INPUTS.NFIELD  = 0 ;

  INPUTS.uM0   = M0FITFLAG_ZBINS_FLAT;
  INPUTS.uzsim = 0 ; // option to set z=simz (i.e., to cheat)

  // stuff for CC prior
  INPUTS.nfile_CCprior  = 0 ;
  INPUTS.varname_pIa[0] = 0 ;
  sprintf(INPUTS.append_varname_missing,"PROB*");
  INPUTS.force_pIa      = -9.0;
  INPUTS.perfect_pIa    = false ;
  INPUTS.maxProbCC_for_sigint = 0.2 ;
  INPUTS.typeIa_ccprior       = -9  ;

  INPUTS_PROBCC_ZERO.USE       = false ;
  INPUTS_PROBCC_ZERO.ntype     = 0 ;
  INPUTS_PROBCC_ZERO.nidsurvey = 0 ;
  INPUTS_PROBCC_ZERO.str_type_list[0]     = 0;
  INPUTS_PROBCC_ZERO.str_idsurvey_list[0] = 0;

  INPUTS.zmin = 0.02 ;
  INPUTS.zmax = 1.02 ;
  

  // ---------------------
  // keep obsolete input parameters (4/24/2012 RK)
  INPUTS.x1min = -6.0;
  INPUTS.x1max = +6.0;
  INPUTS.cmin  = -6.0; 
  INPUTS.cmax  = +6.0;

  INPUTS.logmass_min  = -20.0 ;
  INPUTS.logmass_max  = +20.0 ;
  INPUTS.nbin_logmass =  1 ; 

  INPUTS.chi2max = 1.0E9 ;

  // ---------------------

  INPUTS.cutmask_write = 0;  // default -> write events used in fit

  INPUTS.Nsntype = 0 ;
  sprintf(INPUTS.sntypeString,"NULL");

  INPUTS.nmaxString[0] = 0 ;
  INPUTS.cidFile_data[0] = 0;
  INPUTS.ncidList_data = 0; //djb
  INPUTS.nmax_tot = 999888777 ;
  for(isurvey=0; isurvey<MXIDSURVEY; isurvey++ ) 
    { INPUTS.nmax[isurvey] = 999888777 ; }

  INPUTS.nzbin    =  0 ;
  INPUTS.nlogzbin =  0 ;
  INPUTS.powzbin  =  0.0 ; 
  INPUTS.znhalf   = -9.0 ;
  INPUTS.varname_z[0]   = 0 ; // use default z Z zSPEC
  INPUTS.zbinuser[0]    = 0 ; // comma-sep list of bin edges

  INPUTS.min_per_zbin = MINEVT_PER_ZBIN_DEFAULT ;

  INPUTS.nzbin_ccprior = 0 ; // 0-> use default z-bin size of 0.1

  // global scatter matrix
  INPUTS.sigmB  = 0.00 ; // -> 0 on 4/23/2012 by RK
  INPUTS.sigx1  = 0.00 ;
  INPUTS.sigc   = 0.00 ;
  INPUTS.xi01   = 0.0  ;
  INPUTS.xi0c   = 0.0  ;
  INPUTS.xi1c   = 0.0  ;

  INPUTS.sigint_fix[0] = 0 ;
  
  // default abort is for any fitpar error <=0
  INPUTS.maxerr_abort_c  = 0.0 ;
  INPUTS.maxerr_abort_x1 = 0.0 ;
  INPUTS.maxerr_abort_x0 = 0.0 ;

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
  double g0 = FITPARBOUND_GAMMA[0];
  double g1 = FITPARBOUND_GAMMA[1];

  //         ipar  val   step   bnd0  bnd1 float
  set_fitPar(  1,  0.13, 0.01,  a0,   a1,    1 ); // alpha
  set_fitPar(  2,  3.20, 0.30,  b0,   b1,    1 ); // beta
  set_fitPar(  3,  0.0,  0.02, -0.50, 0.50,  0 ); // dAlpha/dz
  set_fitPar(  4,  0.0,  0.10, -3.00, 3.00,  0 ); // dBeta/dz
  set_fitPar(  5,  0.0,  0.01,  g0,   g1,    0 ); // gamma0
  set_fitPar(  6,  0.0,  0.01, -2.00, 2.00,  0 ); // gamma1
  set_fitPar(  7, 10.0,  0.10,  9.00,11.00,  0 ); // logmass_cen
  set_fitPar(  8,  0.02, 0.01,  0.001,1.00,  0 ); // logmass_tau
  set_fitPar(  9,  0.73, 0.05,  0.00, 2.00,  0 ); // Omega_Lam
  set_fitPar( 10,  0.0,  0.05, -2.00, 5.00,  0 ); // Omega_k
  set_fitPar( 11, -1.0,  0.1,  -3.00, 1.00,  0 ); // w0
  set_fitPar( 12,  0.0,  0.5,  -8.00, 8.00,  0 ); // wa
  set_fitPar( 13,  1.0,  0.05, -0.10,15.00,  0 ); // scale_PCC
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
  ISDATA_REAL = 1;
  INPUTS.blind_cosinePar[IPAR_OL][0] = 0.06; 
  INPUTS.blind_cosinePar[IPAR_OL][1] = 23434. ;

  INPUTS.blind_cosinePar[IPAR_w0][0] = 0.20 ; 
  INPUTS.blind_cosinePar[IPAR_w0][1] = 8430. ;

  sprintf(INPUTS.varname_gamma,"HOST_LOGMASS");
  INPUTS.USE_GAMMA0  = 0 ;

  init_CUTMASK();

  SNFUNPAR_CHI2INFO_OUTPUT.USE        = false ;
  SNFUNPAR_CHI2INFO_OUTPUT.OUTFILE[0] = 0;

  return ;

}   // end of set_defaults

// =====================================
void init_CUTMASK(void) {
  int bit,evtype;

  for (bit=0; bit < MXCUTBIT; bit++ ) {
    CUTMASK_LIST[bit] = (1 << bit);
    for(evtype=0; evtype < MXEVENT_TYPE; evtype++ ) { 
      NSTORE_CUTBIT[evtype][bit] = 0 ; 
    }
  }

  for(bit=0; bit < MXNUM_SAMPLE; bit++ ) 
    { NDATA_BIASCORCUT[0][bit] = NDATA_BIASCORCUT[1][bit] = 0 ; }

  sprintf(CUTSTRING_LIST[CUTBIT_z],         "z"  );
  sprintf(CUTSTRING_LIST[CUTBIT_x1],        "x1" );
  sprintf(CUTSTRING_LIST[CUTBIT_c],         "c"  );
  sprintf(CUTSTRING_LIST[CUTBIT_logmass],   "logmass"  );
  sprintf(CUTSTRING_LIST[CUTBIT_MINBIN],    "min per zbin");
  sprintf(CUTSTRING_LIST[CUTBIT_CUTWIN],    "CUTWIN");
  sprintf(CUTSTRING_LIST[CUTBIT_FIELD],     "CUTFIELD");
  sprintf(CUTSTRING_LIST[CUTBIT_sntype],    "sntype");
  sprintf(CUTSTRING_LIST[CUTBIT_HOST],      "HOST");
  sprintf(CUTSTRING_LIST[CUTBIT_BADERR],    "BADERR among SALT2 fitPar");
  sprintf(CUTSTRING_LIST[CUTBIT_BADCOV],    "BADCOV among SALT2 fitPar");
  sprintf(CUTSTRING_LIST[CUTBIT_SPLITRAN],  "not in SPLITRAN-subset");
  sprintf(CUTSTRING_LIST[CUTBIT_SIMPS],     "sim Prescale" ) ;
  sprintf(CUTSTRING_LIST[CUTBIT_BIASCOR],   "valid BIASCOR" );
  sprintf(CUTSTRING_LIST[CUTBIT_zBIASCOR],  "z-range BIASCOR" );
  sprintf(CUTSTRING_LIST[CUTBIT_x1BIASCOR], "x1-range BIASCOR" );
  sprintf(CUTSTRING_LIST[CUTBIT_cBIASCOR],  "c-range BIASCOR" );
  sprintf(CUTSTRING_LIST[CUTBIT_TRUESNIa],  "TRUESNIa" );
  sprintf(CUTSTRING_LIST[CUTBIT_TRUESNCC],  "TRUESNCC" );
  sprintf(CUTSTRING_LIST[CUTBIT_DUPL],      "DUPLICATE" );
  sprintf(CUTSTRING_LIST[CUTBIT_NMAXCUT],   "NMAXCUT" );
  sprintf(CUTSTRING_LIST[CUTBIT_IDSAMPLE],  "IDSAMPLE" );
  sprintf(CUTSTRING_LIST[CUTBIT_CHI2],      "CHI2" );
  sprintf(CUTSTRING_LIST[CUTBIT_CID],       "CID"  );

} // end init_CUTMASK

// ============================================
void set_fitPar(int ipar, double val, double step, 
		double bndmin, double bndmax, int use) {
  INPUTS.parval[ipar]    = val ;
  INPUTS.parstep[ipar]   = step ;
  INPUTS.parbndmin[ipar] = bndmin;  
  INPUTS.parbndmax[ipar] = bndmax ;
  INPUTS.ipar[ipar]      = use ;
  INPUTS.izpar[ipar]     = -9; // init

  INPUTS.blind_cosinePar[ipar][0] = 0.0 ;
  INPUTS.blind_cosinePar[ipar][1] = 0.0 ;
  sprintf(INPUTS.blindString[ipar], "UNDEFINED");

}


// ******************************************
void read_data(void) {

  // Created Jun 4 2019: Refactored data-read.

  int  NFILE = INPUTS.nfile_data;
  //  int  EVENT_TYPE = EVENT_TYPE_DATA;
  int  NEVT_TOT, NEVT[MXFILE_DATA], ifile, IFILETYPE, LEN_MALLOC ;
  int  ISTART, NROW, isn ;
  char *dataFile ;
  char fnam[] = "read_data" ;

  // -------------- BEGIN --------------

   sprintf(BANNER,"%s", fnam);
   print_banner(BANNER);   

   // read each file quickly to get total size of arrays to malloc
   NEVT_TOT = 0 ;
   for(ifile=0; ifile < NFILE; ifile++ ) {
     dataFile     = INPUTS.dataFile[ifile];
     NEVT[ifile]  = SNTABLE_NEVT(dataFile,TABLENAME_FITRES);  
     NEVT_TOT    += NEVT[ifile];
     printf("\t Found %d events in %s. \n", NEVT[ifile], dataFile);
     fflush(stdout);
   }  

  // malloc arrays for all sim files
  LEN_MALLOC = NEVT_TOT + 10 ;
  INFO_DATA.TABLEVAR.LEN_MALLOC = LEN_MALLOC ;
  malloc_INFO_DATA(+1,LEN_MALLOC);

  // loop again over each data file: read and append INFO_DATA arrays
  for(ifile=0; ifile < NFILE; ifile++ ) {
    dataFile    = INPUTS.dataFile[ifile];
    IFILETYPE   = TABLEFILE_OPEN(dataFile,"read");
    NVAR_ORIG   = SNTABLE_READPREP(IFILETYPE,"FITRES");
    ISTART      = INFO_DATA.TABLEVAR.NSN_ALL ;
    SNTABLE_READPREP_TABLEVAR(ifile, ISTART, NEVT[ifile], &INFO_DATA.TABLEVAR);

    if ( INPUTS.cat_only ) 
      { NROW = NEVT[ifile];  SNTABLE_CLOSE_TEXT(); }
    else
      { NROW = SNTABLE_READ_EXEC(); }    // read entire file; load arrays
    INFO_DATA.TABLEVAR.NSN_ALL += NROW ;

    INFO_DATA.TABLEVAR.EVENT_RANGE[ifile][0] = ISTART ;
    INFO_DATA.TABLEVAR.EVENT_RANGE[ifile][1] = ISTART + NROW - 1;
    sprintf(INFO_DATA.TABLEVAR.INPUT_FILE[ifile],"%s", dataFile);
    store_input_varnames(ifile, &INFO_DATA.TABLEVAR) ;
  }

  // apply parameter blinding (after we know if DATA are real or sim)
  apply_blindpar();

  store_output_varnames(); // May 2020

  int NPASS=0;
  for(isn=0; isn < INFO_DATA.TABLEVAR.NSN_ALL; isn++ ) { 
    compute_more_TABLEVAR(isn, &INFO_DATA.TABLEVAR ); 
    if ( INFO_DATA.TABLEVAR.CUTMASK[isn] == 0 ) { NPASS++ ; }
  }
  if ( NPASS == 0 ) {
    sprintf(c1err,"All DATA events fail cuts");
    sprintf(c2err,"Check cut-windows in SALT2mu input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
  }


  return;

} // end read_data 


// ****************************************
void malloc_INFO_DATA(int opt, int LEN_MALLOC ) {

  // Created June 5 2019
  // Malloc or free INFO_DATA strucure.
  // opt > 0 --> malloc
  // opt < 0 --> free

  int EVENT_TYPE = EVENT_TYPE_DATA;
  int MEMD, MEMI, MEMD1, MEMTOT=0 ;
  float f_MEMORY = 0.0 ;
  char fnam[] = "malloc_INFO_DATA";

  // ------------- BEGIN --------------

  if ( opt > 0 ) {

    // start with generic malloc to read any FITRES file
    f_MEMORY = malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_DATA.TABLEVAR);

    CUTMASK_POINTER[EVENT_TYPE]         = &INFO_DATA.TABLEVAR.CUTMASK[0];
    NALL_CUTMASK_POINTER[EVENT_TYPE]    = &INFO_DATA.TABLEVAR.NSN_ALL;
    NPASS_CUTMASK_POINTER[EVENT_TYPE]   = &INFO_DATA.TABLEVAR.NSN_PASSCUTS;
    NREJECT_CUTMASK_POINTER[EVENT_TYPE] = &INFO_DATA.TABLEVAR.NSN_REJECT;

    // DATA-specific mallocs ...
    MEMD  = LEN_MALLOC * sizeof(double);
    MEMI  = LEN_MALLOC * sizeof(int);
    MEMD1 = LEN_MALLOC * sizeof(double*);

    INFO_DATA.mumodel        = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.M0             = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mu             = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerr          = (double*) malloc(MEMD); MEMTOT+=MEMD; 
    INFO_DATA.muerr_raw      = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mures          = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mupull         = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerr_last     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerrsq_last   = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.sigCC_last     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.sqsigCC_last   = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muCOVscale     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muBias         = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muBiasErr      = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muBias_zinterp = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.chi2           = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.probcc_beams   = (double*) malloc(MEMD); MEMTOT+=MEMD;
    f_MEMORY += (float)(MEMTOT/1.0E6) ;

    // fitParBias[isn][ipar] to allow passing mB,x1,c via array
    f_MEMORY += malloc_double2D(opt, LEN_MALLOC, NLCPAR+1, 
				&INFO_DATA.fitParBias ); // <== returned

    f_MEMORY += malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
					    &INFO_DATA.FITPARBIAS_ALPHABETA );

    // MUCOVSCALE_AB[isn][ia][ib][ig]
    f_MEMORY += malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg,
				&INFO_DATA.MUCOVSCALE_ALPHABETA ); //<==return
    

    INFO_DATA.MEMORY = f_MEMORY;
    printf("\t %s: %6.3f MB \n", fnam, f_MEMORY); fflush(stdout);
  }
  else {
    malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_DATA.TABLEVAR);
    free(INFO_DATA.mumodel);
    free(INFO_DATA.M0);
    free(INFO_DATA.mu);
    free(INFO_DATA.muerr);
    free(INFO_DATA.muerr_raw);
    free(INFO_DATA.mures);
    free(INFO_DATA.mupull);
    free(INFO_DATA.muerr_last);
    free(INFO_DATA.muerrsq_last);
    free(INFO_DATA.sigCC_last);
    free(INFO_DATA.sqsigCC_last);
    free(INFO_DATA.muBias);
    free(INFO_DATA.muBiasErr);
    free(INFO_DATA.muBias_zinterp);

    malloc_double2D(opt, LEN_MALLOC, NLCPAR+1, &INFO_DATA.fitParBias ); 
    malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
				&INFO_DATA.FITPARBIAS_ALPHABETA );
    malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg, 
		    &INFO_DATA.MUCOVSCALE_ALPHABETA ); 

  }


  return ;

} // end malloc_INFO_DATA

// ****************************************
void malloc_INFO_BIASCOR(int opt, int LEN_MALLOC ) {

  // Created June 5 2019
  // Malloc or free INFO_BIASCOR structure
  // opt > 0 --> malloc
  // opt < 0 --> free

  int  NSAMPLE = NSAMPLE_BIASCOR;
  int  EVENT_TYPE = EVENT_TYPE_BIASCOR ;
  int  MEMB = LEN_MALLOC * sizeof(int8_t);     // one byte
  int  MEMBIAS = NSAMPLE * sizeof(FITPARBIAS_DEF*);
  int  MEMCOV  = NSAMPLE * sizeof(float*);

  int  MEMADD=0 ;
  float f_MEMORY = 0.0 ;
  char fnam[] = "malloc_INFO_BIASCOR";

  // ------------- BEGIN --------------

  if ( opt > 0 ) {
    
    // start with generic malloc to read any FITRES file
    f_MEMORY = malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_BIASCOR.TABLEVAR);

    // set pointers to track event stats vs. EVENT_TYPE
    CUTMASK_POINTER[EVENT_TYPE]         = &INFO_BIASCOR.TABLEVAR.CUTMASK[0];
    NALL_CUTMASK_POINTER[EVENT_TYPE]    = &INFO_BIASCOR.TABLEVAR.NSN_ALL;
    NPASS_CUTMASK_POINTER[EVENT_TYPE]   = &INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS;
    NREJECT_CUTMASK_POINTER[EVENT_TYPE] = &INFO_BIASCOR.TABLEVAR.NSN_REJECT ;

    // BIASCOR-specific mallocs ...

    MEMADD += MEMB; INFO_BIASCOR.iz  = (int8_t*)malloc(MEMB);
    MEMADD += MEMB; INFO_BIASCOR.IZ  = (int8_t*)malloc(MEMB);
    MEMADD += MEMB; INFO_BIASCOR.IA  = (int8_t*)malloc(MEMB);
    MEMADD += MEMB; INFO_BIASCOR.IB  = (int8_t*)malloc(MEMB);
    MEMADD += MEMB; INFO_BIASCOR.IG  = (int8_t*)malloc(MEMB);

    // allocate IDSAMPLE-dimention here. The NCELL dimension
    // is allocate later in set_MAPCELL_biascor (when NCELL is known)
    INFO_BIASCOR.FITPARBIAS = (FITPARBIAS_DEF**) malloc ( MEMBIAS );
    INFO_BIASCOR.MUCOVSCALE = (float **        ) malloc ( MEMCOV  );    

    // print memory consumption to stdout
    f_MEMORY += ((float)MEMADD) / 1.0E6 ;
    INFO_BIASCOR.MEMORY = f_MEMORY;
    printf("\t %s: %6.3f MB \n", fnam, f_MEMORY); fflush(stdout);


  }
  else {
    // free memory
    malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_BIASCOR.TABLEVAR);
    free(INFO_BIASCOR.iz); free(INFO_BIASCOR.IZ);
    free(INFO_BIASCOR.IA); free(INFO_BIASCOR.IB);
  }


  return ;
  

} // end malloc_INFO_BIASCOR

// ****************************************
void malloc_INFO_CCPRIOR(int opt, int LEN_MALLOC, int LEN_MALLOC_CUTS) {

  // Created June 5 2019
  // Malloc or free INFO_CCPRIOR structure.
  // opt > 0 --> malloc
  // opt < 0 --> free

  //  int EVENT_TYPE   = EVENT_TYPE_CCPRIOR ;
  int USE_BIASCOR  = ( INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS > 0 ) ;
  int MEMTOT=0 ;
  float f_MEMORY = 0.0 ;
  char fnam[] = "malloc_INFO_CCPRIOR";

  // ------------- BEGIN --------------

  if ( opt > 0 ) {
    
    // start with generic malloc to read any FITRES file
    if ( LEN_MALLOC ) { 
      f_MEMORY = malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_CCPRIOR.TABLEVAR); 
      if ( USE_BIASCOR ) {
	f_MEMORY += 
	  malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
				      &INFO_CCPRIOR.FITPARBIAS_ALPHABETA );
      }
    }

    if ( LEN_MALLOC_CUTS ) {
      f_MEMORY = 
	malloc_TABLEVAR(opt, LEN_MALLOC_CUTS, &INFO_CCPRIOR.TABLEVAR_CUTS); 
   
      if ( USE_BIASCOR ) {
	f_MEMORY += 
	  malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC_CUTS,
				      &INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS);
      }
    }

    f_MEMORY += (float)(MEMTOT)/1.0E6 ;
    INFO_CCPRIOR.MEMORY = f_MEMORY;
    printf("\t %s: %6.3f MB \n", fnam, f_MEMORY); fflush(stdout);

  }
  else {
    if ( LEN_MALLOC ) {
      malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_CCPRIOR.TABLEVAR);
      if ( USE_BIASCOR ) {
	malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
				    &INFO_CCPRIOR.FITPARBIAS_ALPHABETA );
      }
    }

  }


  return ;
  

} // end malloc_INFO_CCPRIOR


// ****************************************
float malloc_TABLEVAR(int opt, int LEN_MALLOC, TABLEVAR_DEF *TABLEVAR) {

  // Created June 5 2019
  // Malloc or free TABLEVAR elements (works for DATA,BIASCOR,CCPRIOR).
  // opt > 0 --> malloc TABLEVAR arrays with LEN_MALLOC elements
  // opt < 0 --> free TABLEVAR
  //
  // Functions returns memory allocated in Mega-bytes

  int EVENT_TYPE = TABLEVAR->EVENT_TYPE;
  int IS_DATA    = (EVENT_TYPE == EVENT_TYPE_DATA);
  //int IS_BIASCOR = (EVENT_TYPE == EVENT_TYPE_BIASCOR);
  //  int IS_CCPRIOR = (EVENT_TYPE == EVENT_TYPE_CCPRIOR);


  int MEMF    = LEN_MALLOC  * sizeof(float);
  int MEMF2   = LEN_MALLOC  * sizeof(float**);
  int MEMI    = LEN_MALLOC  * sizeof(int);
  int MEMS    = LEN_MALLOC  * sizeof(short int);
  int MEMB    = LEN_MALLOC  * sizeof(bool);
  int MEMC    = LEN_MALLOC  * sizeof(char*);
  int MEMC2   = MXCHAR_CCID * sizeof(char);

  bool DOBIAS_MU = ( INPUTS.opt_biasCor & MASK_BIASCOR_MU     ) ;
  bool IDEAL     = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;

  int  NLCPAR_LOCAL = NLCPAR ; // beware that ILCPAR_MIN/MAX isn't set yet
  if ( DOBIAS_MU ) { NLCPAR_LOCAL++ ; }

  bool USE_FIELD = (INPUTS.use_fieldGroup_biasCor || INPUTS.NFIELD>0);
  int long long MEMTOT=0;
  float f_MEMTOT;
  int  i, isn, MEMF_TMP1, MEMF_TMP ;
  //  char fnam[] = "malloc_TABLEVAR" ;

  // ------------- BEGIN --------------

  if ( opt > 0 ) {   
    
    TABLEVAR->name =  (char**)malloc(MEMC);
    for(i=0; i<LEN_MALLOC; i++ ) 
      { TABLEVAR->name[i] = (char*)malloc(MEMC2); MEMTOT+=MEMC2; }

    if ( USE_FIELD ) {
      TABLEVAR->field =  (char**)malloc(MEMC);
      for(i=0; i<LEN_MALLOC; i++ ) 
	{ TABLEVAR->field[i] = (char*)malloc(MEMC2); MEMTOT+=MEMC2; }  
    }

    for (i=0; i < NLCPAR_LOCAL; i++ ) {
      TABLEVAR->fitpar[i]     = (float *) malloc(MEMF); MEMTOT+=MEMF;
      TABLEVAR->fitpar_err[i] = (float *) malloc(MEMF); MEMTOT+=MEMF;
    }

    // allocate 3D array: covmat[isn][ipar0][ipar1]
    // Note that this array cannot be read or passed as ISN-array,
    // but it can be passed as 2D array to functions.
    TABLEVAR->covmat_fit = (float ***) malloc(MEMF2); MEMTOT+=MEMF2;
    TABLEVAR->covmat_tot = (float ***) malloc(MEMF2); MEMTOT+=MEMF2;    
    MEMF_TMP1 = NLCPAR * sizeof(float*);
    MEMF_TMP  = NLCPAR * sizeof(float);
    for (isn=0; isn < LEN_MALLOC; isn++ ) {
      TABLEVAR->covmat_fit[isn] = 
	(float **) malloc(MEMF_TMP1); MEMTOT += MEMF_TMP1 ;
      TABLEVAR->covmat_tot[isn] = 
	(float **) malloc(MEMF_TMP1); MEMTOT += MEMF_TMP1 ;
      for (i=0; i < NLCPAR; i++ ) {
	TABLEVAR->covmat_fit[isn][i] = 
	  (float *) malloc(MEMF_TMP); MEMTOT+=MEMF_TMP ;
	TABLEVAR->covmat_tot[isn][i] = 
	  (float *) malloc(MEMF_TMP); MEMTOT+=MEMF_TMP ;
      }
    }

    if ( IDEAL ) { 
      TABLEVAR->x0_ideal = (float *) malloc(MEMF); MEMTOT+=MEMF;
      for (i=0; i < NLCPAR_LOCAL ; i++ ) 
	{ TABLEVAR->fitpar_ideal[i] = (float *) malloc(MEMF); MEMTOT+=MEMF; }
    }

    TABLEVAR->warnCov       = (bool  *) malloc(MEMB); MEMTOT+=MEMB;
    TABLEVAR->x0            = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->x0err         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x0x1      = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x0c       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x1c       = (float *) malloc(MEMF); MEMTOT+=MEMF;

    TABLEVAR->zhd           = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhderr        = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhel          = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhelerr       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->vpec          = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->vpecerr       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->logmass       = (float *) malloc(MEMF); MEMTOT+=MEMF; 
    TABLEVAR->snrmax        = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->pIa           = (float *) malloc(MEMF); MEMTOT+=MEMF; 

    TABLEVAR->IDSURVEY      = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->IDSAMPLE      = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->SNTYPE        = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->OPT_PHOTOZ    = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->IZBIN         = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->CUTMASK       = (int *) malloc(MEMI); MEMTOT+=MEMI;

    TABLEVAR->ICUTWIN_GAMMA = -9 ;
    for(i=0; i < INPUTS.NCUTWIN; i++ ) 
      { MEMTOT += malloc_TABLEVAR_CUTVAL(LEN_MALLOC,i, TABLEVAR ); }
  
    TABLEVAR->SIM_NONIA_INDEX  = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->SIM_ZCMB         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_VPEC         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_MU           = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_ALPHA        = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_BETA         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_GAMMADM      = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_X0           = (float *) malloc(MEMF); MEMTOT+=MEMF;

    for (i=0; i < NLCPAR_LOCAL ; i++ ) 
      { TABLEVAR->SIM_FITPAR[i] = (float *) malloc(MEMF); MEMTOT+=MEMF; } 

    if ( IS_DATA ) {
      TABLEVAR->mumodel       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  }
  else {
    // free memory

    for(i=0; i<LEN_MALLOC; i++ ) { free(TABLEVAR->name[i]); }
    free(TABLEVAR->name);

    if ( USE_FIELD ) {
      for(i=0; i<LEN_MALLOC; i++ ) { free(TABLEVAR->field[i]) ; }
      free(TABLEVAR->field);
    }

    for (i=0; i < NLCPAR_LOCAL; i++ ) {
      free(TABLEVAR->fitpar[i]);
      free(TABLEVAR->fitpar_err[i]);
    }

    
    if ( IDEAL ) {
      for (i=0; i < NLCPAR_LOCAL ; i++ ) 
	{ free(TABLEVAR->fitpar_ideal[i]); }
      free(TABLEVAR->x0_ideal);
    }

    free(TABLEVAR->x0);
    free(TABLEVAR->x0err);
    free(TABLEVAR->zhd);
    free(TABLEVAR->zhderr);
    free(TABLEVAR->zhel);
    free(TABLEVAR->zhelerr);
    free(TABLEVAR->vpec);
    free(TABLEVAR->vpecerr);
    free(TABLEVAR->snrmax);


    free(TABLEVAR->COV_x0x1);
    free(TABLEVAR->COV_x0c);
    free(TABLEVAR->COV_x1c);
    free(TABLEVAR->warnCov);
    
    free(TABLEVAR->IDSURVEY);
    free(TABLEVAR->IDSAMPLE);
    free(TABLEVAR->SNTYPE);
    free(TABLEVAR->OPT_PHOTOZ);
    free(TABLEVAR->IZBIN);
    free(TABLEVAR->CUTMASK);

    for(i=0; i < INPUTS.NCUTWIN; i++ ) { 
      if ( INPUTS.LCUTWIN_RDFLAG[i] ) { free(TABLEVAR->CUTVAL[i]); }
    }

    free(TABLEVAR->SIM_NONIA_INDEX);
    free(TABLEVAR->SIM_ZCMB);
    free(TABLEVAR->SIM_VPEC);
    free(TABLEVAR->SIM_MU);
    free(TABLEVAR->SIM_ALPHA);
    free(TABLEVAR->SIM_BETA);
    free(TABLEVAR->SIM_GAMMADM);
    free(TABLEVAR->SIM_X0);

    for (i=0; i < NLCPAR_LOCAL ; i++ )
      { free(TABLEVAR->SIM_FITPAR[i]); }

    free(TABLEVAR->pIa);
    if ( IS_DATA ) {
      free(TABLEVAR->mumodel);
    }

    return(0.0);
  }

  // - - - - -
  return(0.0) ;

} // end malloc_TABLEVAR


// ***************************************************
int malloc_TABLEVAR_CUTVAL(int LEN_MALLOC, int icut,
			     TABLEVAR_DEF *TABLEVAR ) {

  int MEMTOT = 0;
  int MEMF   = LEN_MALLOC * sizeof(float);
  char *CUTNAME = INPUTS.CUTWIN_NAME[icut] ;
  //  char fnam[] = "malloc_TABLEVAR_CUTVAL" ;

  // --------------- BEGIN -----------------

  INPUTS.LCUTWIN_RDFLAG[icut] = false ;

  if ( strcmp(CUTNAME,"x0") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->x0; }

  else if ( strcmp(CUTNAME,"x1") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->fitpar[INDEX_x1]; }
  else if ( strcmp(CUTNAME,"x1ERR") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->fitpar_err[INDEX_x1]; }

  else if ( strcmp(CUTNAME,"c") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->fitpar[INDEX_c]; }
  else if ( strcmp(CUTNAME,"cERR") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->fitpar_err[INDEX_c]; }

  else if ( strcmp(CUTNAME,"z") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->zhd; }
  else if ( strcmp(CUTNAME,"zhd") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->zhd; }
  else if ( strcmp(CUTNAME,"zERR") == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->zhderr; }

  else if ( strcmp(CUTNAME,INPUTS.varname_pIa) == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->pIa; }

  // IDSURVEY and SNTYPE are int ??

  else {
    INPUTS.LCUTWIN_RDFLAG[icut] = true ;
    TABLEVAR->CUTVAL[icut] = (float*)malloc(MEMF);  MEMTOT += MEMF ; 
    if ( strcmp(CUTNAME,INPUTS.varname_gamma) == 0 ) 
      { TABLEVAR->ICUTWIN_GAMMA = icut;       }
  }


  return(MEMTOT) ;

} // end malloc_TABLEVAR_CUTVAL

// ***************************************************
float malloc_FITPARBIAS_ALPHABETA(int opt, int LEN_MALLOC,
				  FITPARBIAS_DEF *****FITPARBIAS ) {

  // Created Jun 6 2019
  // Malloc 3D FITPARBIAS[LEN_MALLOC][MXa][MXb][MXg]
  // opt > 0 --> malloc
  // opt < 0 --> free
  //
  // Function returns memory allocated, in Mbytes.
  // 
  int ia, ib, isn;
  int MEML  = sizeof(FITPARBIAS_DEF***) * LEN_MALLOC ; 
  int MEMa  = sizeof(FITPARBIAS_DEF**)  * MXa ;
  int MEMb  = sizeof(FITPARBIAS_DEF*)   * MXb ;
  int MEMg  = sizeof(FITPARBIAS_DEF )   * MXg ;
  int MEMTOT = 0 ;
  float f_MEMTOT;
  // ----------- BEGIN ----------

  if ( opt > 0 ) {

    *FITPARBIAS = (FITPARBIAS_DEF****) malloc(MEML) ; MEMTOT+=MEML;
    for(isn=0; isn<LEN_MALLOC; isn++ ) {
      (*FITPARBIAS)[isn] = (FITPARBIAS_DEF***) malloc(MEMa) ; MEMTOT+=MEMa;
      for(ia=0; ia<MXa; ia++ ) {
	(*FITPARBIAS)[isn][ia] = (FITPARBIAS_DEF**)malloc(MEMb); MEMTOT+=MEMb;
	for(ib=0; ib<MXb; ib++ ) {
	  (*FITPARBIAS)[isn][ia][ib] = (FITPARBIAS_DEF*) malloc(MEMg);
	  MEMTOT+=MEMg;
	}
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(isn=0; isn<LEN_MALLOC; isn++ ) {
      for(ia=0; ia<MXa; ia++ ) {	
	for(ib=0; ib<MXb; ib++ ) {
	  free((*FITPARBIAS)[isn][ia][ib]);  
	}
	free((*FITPARBIAS)[isn][ia]);  
      }
      free( (*FITPARBIAS)[isn]);
    }
    free( (*FITPARBIAS) );
  }

  return(0.0) ;

} // end malloc_FITPARBIAS_ALPHABETA


// ***************************************
float malloc_double2D(int opt, int LEN1, int LEN2, double ***array2D ) {
  // Created Jun 11 2019
  // Malloc array2D[LEN1][LEN2]  (intended for LEN1=NSN, LEN2=NCLPAR)
  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1 ;
  int MEM1 = LEN1 * sizeof(double*); 
  int MEM2 = LEN2 * sizeof(double);
  // ----------- BEGIN -------------

  if ( opt > 0 ) {

    *array2D = (double**) malloc(MEM1) ; MEMTOT += MEM1;
    for(i1=0; i1< LEN1; i1++ ) {
      (*array2D)[i1] = (double*) malloc(MEM2) ; MEMTOT += MEM2;
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1 < LEN1; i1++ ) { free((*array2D)[i1]); }
    free(array2D[i1]) ;    
  }


  return(f_MEMTOT);
}

float malloc_double3D(int opt, int LEN1, int LEN2, int LEN3, 
		      double ****array3D ) {
  // Created Jun 11 2019
  // Malloc array3D[LEN1][LEN2][LEN3] 
  //   (intended for LEN1=NSN, LEN2=MXa, LEN3=MXb)

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2 ;
  int MEM1 = LEN1 * sizeof(double**); 
  int MEM2 = LEN2 * sizeof(double*);
  int MEM3 = LEN3 * sizeof(double);
  // ----------- BEGIN -------------

  if ( opt > 0 ) {

    *array3D = (double***) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array3D)[i1] = (double**) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array3D)[i1][i2] = (double*) malloc(MEM3); MEMTOT+=MEM3;
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {  free( (*array3D)[i1][i2]); }
      free( (*array3D)[i1]) ;
    }
    free(array3D);
  }


  return(f_MEMTOT);

} // end malloc_double3D


float malloc_double4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
		      double *****array4D ) {
  // Created July 2019
  // Malloc array3D[LEN1][LEN2][LEN3][LEN4] 
  //   (intended for LEN1=NSN, LEN2=MXa, LEN3=MXb, LEN4=MXg)

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2, i3 ;
  int MEM1 = LEN1 * sizeof(double***); 
  int MEM2 = LEN2 * sizeof(double**);
  int MEM3 = LEN3 * sizeof(double*);
  int MEM4 = LEN4 * sizeof(double);
  // ----------- BEGIN -------------

  if ( opt > 0 ) {

    *array4D = (double****) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array4D)[i1] = (double***) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array4D)[i1][i2] = (double**) malloc(MEM3); MEMTOT+=MEM3;
	for(i3=0; i3<LEN3; i3++ ) {
	  (*array4D)[i1][i2][i3] = (double*) malloc(MEM4); MEMTOT+=MEM4;
	}
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {
	for(i3=0; i3<LEN3; i3++ ) {
	  free( (*array4D)[i1][i2][i3] ); 
	}
	free( (*array4D)[i1][i2]); 
      }
      free( (*array4D)[i1]) ;
    }
    free(array4D);
  }


  return(f_MEMTOT);

}   // end malloc_double3D

// ***********************************************
void SNTABLE_READPREP_TABLEVAR(int IFILE, int ISTART, int LEN, 
			       TABLEVAR_DEF *TABLEVAR) {

  // Created Jun 4 2019
  // prepare reading TABLEVAR arrays.
  // EVENT_TYPE specifies DATA,BIASCOR, or CCPRIOR,
  // but most of this function is independent of EVENT_TYPE.
  //
  // Inputs:
  //   IFILE      : file index 
  //   EVENT_TYPE : specifies DATA, BIASCOR, or CCPRIOR
  //   ISTART     : start index to store data
  //   LEN        : size of array to read
  //   TABLEVAR   : structure with arrays to initialize
  //
  //  May  8 2020: add IFILE argument.
  //  May 20 2020: check NFIELD

  int EVENT_TYPE = TABLEVAR->EVENT_TYPE;
  int IS_DATA    = ( EVENT_TYPE == EVENT_TYPE_DATA);
  int IS_BIASCOR = ( EVENT_TYPE == EVENT_TYPE_BIASCOR);
  //  int IS_CCPRIOR = ( EVENT_TYPE == EVENT_TYPE_CCPRIOR);

  int VBOSE = 1 ;  // verbose, but no abort on missing variable
  int FIRSTFILE = ( ISTART == 0 ) ;
  int USE_FIELD = ( INPUTS.use_fieldGroup_biasCor>0 || INPUTS.NFIELD>0);
  int IDEAL           = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;

  int  icut, ivar, ivar2, irow, id ;
  bool RDFLAG ;
  char vartmp[MXCHAR_VARNAME], *cutname, str_z[MXCHAR_VARNAME]; 
  char str_zerr[MXCHAR_VARNAME]; 

  char fnam[] = "SNTABLE_READPREP_TABLEVAR" ;

  // ----------- BEGIN -------------


  // init flags on first file
  if ( FIRSTFILE ) {
    TABLEVAR->NSN_ALL           = 0;
    TABLEVAR->NSN_REJECT        = 0;
    TABLEVAR->NSN_PASSCUTS      = 0;
    TABLEVAR->NSN_REJECT        = 0;
    TABLEVAR->IS_SIM            = false ;
    TABLEVAR->IS_DATA           = false ;
    TABLEVAR->IVAR_VPEC         = -9 ;
    TABLEVAR->IVAR_SIM_VPEC     = -9 ;
    TABLEVAR->IVAR_OPT_PHOTOZ   = -9 ;
    TABLEVAR->IVAR_SNTYPE       = -9 ;
    TABLEVAR->IVAR_SIM_GAMMADM  = -9 ;

    for(icut=0; icut < MXCUTBIT; icut++ ) 
      { TABLEVAR->NSN_CUTBIT[icut] = 0 ; }    

    for(icut=0; icut < INPUTS.NCUTWIN; icut++ )  
      { TABLEVAR->DOFLAG_CUTWIN[icut] = 0; }

    for(id=0; id<MXIDSURVEY; id++ ) {
      TABLEVAR->NSN_PER_SURVEY[id]  = 0 ; 
      TABLEVAR->zMIN_PER_SURVEY[id] = +999.0 ; 
      TABLEVAR->zMAX_PER_SURVEY[id] = -999.0 ; 
    }    

  }

  // initialize a few arrays
  for(irow=ISTART; irow < ISTART+LEN; irow++ ) {
    if ( USE_FIELD ) { TABLEVAR->field[irow][0] = 0 ; }
    TABLEVAR->IDSURVEY[irow]   = -9 ;
    TABLEVAR->IDSAMPLE[irow]   = -9 ;
    TABLEVAR->CUTMASK[irow]    =  0 ;
    TABLEVAR->OPT_PHOTOZ[irow] =  0 ;
    TABLEVAR->SNTYPE[irow]     = -9 ;

    TABLEVAR->vpec[irow]       =  0.0 ;
    TABLEVAR->vpecerr[irow]    =  0.0 ;
    TABLEVAR->zhd[irow]        = -9.0 ;
    TABLEVAR->zhderr[irow]     = -9.0 ;
    TABLEVAR->zhel[irow]       = -9.0 ;
    TABLEVAR->zhelerr[irow]    = -9.0 ;
    TABLEVAR->snrmax[irow]     =  0.0 ;
    TABLEVAR->warnCov[irow]    =  false ;

    TABLEVAR->SIM_NONIA_INDEX[irow]  = -9 ;
    TABLEVAR->SIM_X0[irow]           = -9.0 ;
    TABLEVAR->SIM_FITPAR[0][irow]    = -9.0 ;
    TABLEVAR->SIM_FITPAR[1][irow]    = -9.0 ;
    TABLEVAR->SIM_FITPAR[2][irow]    = -9.0 ;
    TABLEVAR->SIM_ALPHA[irow]        = -9.0 ;
    TABLEVAR->SIM_BETA[irow]         = -9.0 ;
    TABLEVAR->SIM_GAMMADM[irow]      =  0.0 ; // allow missing gammadm
    TABLEVAR->SIM_MU[irow]           = -9.0 ;
    TABLEVAR->SIM_ZCMB[irow]         = -9.0 ;
    TABLEVAR->SIM_VPEC[irow]         =  0.0 ;
  }

  // - - - - - - - - prep strings - - - - - - - 
  sprintf(vartmp, "CID:C*%d  CCID:C*%d", MXCHAR_CCID, MXCHAR_CCID); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->name[ISTART], 
			  LEN, VBOSE) ;

  if ( USE_FIELD ) {
    sprintf(vartmp, "FIELD:C*%d", MXCHAR_CCID ); 
    ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->field[ISTART], 
				   LEN, VBOSE) ;
    if ( ivar < 0 ) {
      sprintf(c1err,"Required FIELD column missing");
      sprintf(c2err,"Check CUTFIELD or fieldGroup_biascor keys");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
    }
  }

  sprintf(vartmp,"IDSURVEY:S" ); // S -> short int
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->IDSURVEY[ISTART], 
			  LEN, VBOSE );

  sprintf(vartmp,"OPT_PHOTOZ:S" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->OPT_PHOTOZ[ISTART],
				 LEN, VBOSE);
  TABLEVAR->IVAR_OPT_PHOTOZ = ivar;

  sprintf(vartmp,"TYPE:S" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SNTYPE[ISTART],
				 LEN, VBOSE);
  TABLEVAR->IVAR_SNTYPE = ivar;


  // readshift and peculiar velocity
  str_z[0] = str_zerr[0] = 0 ;
  get_zString(str_z,str_zerr,"F");
  SNTABLE_READPREP_VARDEF(str_z, &TABLEVAR->zhd[ISTART], 
			  LEN, VBOSE);
  SNTABLE_READPREP_VARDEF(str_zerr, &TABLEVAR->zhderr[ISTART], 
			  LEN, VBOSE);

  sprintf(vartmp,"zHEL:F"); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zhel[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"zHELERR:F"); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zhelerr[ISTART], 
			  LEN, VBOSE);

  sprintf(vartmp,"VPEC:F" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->vpec[ISTART], 
				 LEN, VBOSE);
  sprintf(vartmp,"VPEC_ERR:F VPECERR:F" );
  ivar2 = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->vpecerr[ISTART], 
				  LEN, VBOSE);
  if ( ivar>0 && ivar2 > 0 ) { TABLEVAR->IVAR_VPEC = ivar; }

  sprintf(vartmp,"SNRMAX:F  SNRMAX1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->snrmax[ISTART],  
			  LEN, VBOSE );


  // - - - - - - - - - - LC fit params - - - - - - - - - - 
  
  sprintf(vartmp,"x0:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->x0[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"x0ERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->x0err[ISTART], 
			  LEN, VBOSE);

  sprintf(vartmp,"x1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar[INDEX_x1][ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"x1ERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar_err[INDEX_x1][ISTART], 
			  LEN, VBOSE);

  sprintf(vartmp,"c:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar[INDEX_c][ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"cERR:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar_err[INDEX_c][ISTART], 
			  LEN, VBOSE);


  sprintf(vartmp,"COVx0x1:F COV_x1_x0:F");
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->COV_x0x1[ISTART],  
			  LEN, VBOSE);

  sprintf(vartmp,"COVx0c:F COV_c_x0:F");
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->COV_x0c[ISTART],
			  LEN, VBOSE);

  sprintf(vartmp,"COVx1c:F COV_x1_c:F");
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->COV_x1c[ISTART], 
			  LEN, VBOSE );

  char *varname_pIa = INPUTS.varname_pIa ;
  if ( strlen(varname_pIa) > 0 ) {
    sprintf(vartmp,"%s:F", varname_pIa);
    ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->pIa[ISTART],
				   LEN, VBOSE );
    TABLEVAR->IVAR_pIa[IFILE] = ivar; // map valid ivar with each file
  }

  // - - - - - - -
  //read CUTWIN variables 
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
    cutname = INPUTS.CUTWIN_NAME[icut]; 
    RDFLAG  = INPUTS.LCUTWIN_RDFLAG[icut] ;
    sprintf(vartmp, "%s:F", cutname );
    if ( !usesim_CUTWIN(vartmp)  ) { continue ; }

    if ( RDFLAG ) {
      ivar = SNTABLE_READPREP_VARDEF(vartmp,&TABLEVAR->CUTVAL[icut][ISTART], 
				     LEN, VBOSE );
    }
    else {
      ivar = IVAR_READTABLE_POINTER(cutname); // May 8 2020 
    }
    TABLEVAR->DOFLAG_CUTWIN[icut] = set_DOFLAG_CUTWIN(ivar,icut,IS_DATA);
  }


  // - - - - - - - - SIM_XXX - - - - - - - - - -
  sprintf(vartmp,"SIM_NONIA_INDEX:S SIM_TEMPLATE_INDEX:S" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_NONIA_INDEX[ISTART], 
				 LEN, VBOSE );
  FOUNDKEY_SIM=0;

  if ( ivar >=0 ) 
    { TABLEVAR->IS_SIM = true ;   FOUNDKEY_SIM=1; }
  else
    { TABLEVAR->IS_DATA = true;   return ; }

  // here and below is for simulated data 

  // note that IS_DATA refers to datafile= argument, and can be
  // real data or simulated data. ISDATA_REAL is false for sim data.
  if ( IS_DATA ) { 
    ISDATA_REAL = 0 ;   // not real data 

    // if the 64 blind-sim bit isn't set by user, set blindFlag=0
    if ( (INPUTS.blindFlag & BLINDMASK_SIM)==0 ) { INPUTS.blindFlag=0; }
  }
  
  sprintf(vartmp,"SIMx0:F SIM_x0:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_X0[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"SIMmB:F SIM_mB:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_FITPAR[INDEX_mB][ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"SIMx1:F SIM_x1:F SIM_SALT2x1:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_FITPAR[INDEX_x1][ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"SIMc:F SIM_c:F SIM_SALT2c:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_FITPAR[INDEX_c][ISTART], 
			  LEN, VBOSE);

  // - - - - - c - - - - - -

  sprintf(vartmp,"SIMalpha:F SIM_alpha:F" ) ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_ALPHA[ISTART],
				 LEN, VBOSE );
  sprintf(vartmp,"SIMbeta:F  SIM_beta:F" ) ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_BETA[ISTART],
				 LEN, VBOSE );

  
  
  sprintf(vartmp,"SIM_gammaDM:F" ) ; 
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_GAMMADM[ISTART],
				 LEN, VBOSE );
  TABLEVAR->IVAR_SIM_GAMMADM = ivar;
  

  // - - - - - - - - - - - - - - - - - - - 

  // true z & MU
  sprintf(vartmp,"SIM_ZCMB:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_ZCMB[ISTART],
			  LEN, VBOSE);

  sprintf(vartmp,"SIM_DLMAG:F" ) ;
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_MU[ISTART],
			  LEN, VBOSE);

  // true VPEC (Aug 2017)
  sprintf(vartmp,"SIM_VPEC:F" ) ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_VPEC[ISTART], 
				 LEN, VBOSE);
  if ( ivar > 0 ) { TABLEVAR->IVAR_SIM_VPEC = ivar ; }

  if ( IS_BIASCOR && IDEAL ) {
  sprintf(vartmp,"x0_ideal:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->x0_ideal[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"x1_ideal:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar_ideal[INDEX_x1][ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"c_ideal:F" );
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->fitpar_ideal[INDEX_c][ISTART], 
			  LEN, VBOSE);
    
  }

  return ;

} // end  SNTABLE_READPREP_TABLEVAR


// ***********************************************
void compute_more_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR ) {

  // Created Jun 2019
  // After reading table and loadig TABLEVAR arrays, 
  // this function defines 'more' TABLEVAR quantities;
  //  * mB, mBerr
  //  * mumodel
  //  * cov matrices
  //  * fix cov matries
  //  * set_CUTMASK
  //
  // Inputs:
  //   EVENT_TYPE : points to DATA, BIASCOR, or CCPRIOR
  //   ISN        : SN index for TABLEVAR arrays
  //   TABLEVAR   : structure of arrays
  //
  //  Beware that we work here with float to save memory,
  //  particularly for large BIASCOR samples.
  //
  // Aug 22 2019: set logmass to table value, or to p7
  //
  // Feb 24 2020: load TABLEVAR->fitpar[INDEX_mu][ISN]
  //

  int EVENT_TYPE   = TABLEVAR->EVENT_TYPE;
  bool IS_DATA     = (EVENT_TYPE == EVENT_TYPE_DATA);
  bool IS_BIASCOR  = (EVENT_TYPE == EVENT_TYPE_BIASCOR);
  //  int IS_CCPRIOR  = (EVENT_TYPE == EVENT_TYPE_CCPRIOR);
  bool IDEAL       = (INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  bool FIRST_EVENT = (ISN == 0) ;
  //  bool LAST_EVENT  = (ISN == TABLEVAR->NSN_ALL-1 );
  char *STRTYPE   = STRING_EVENT_TYPE[EVENT_TYPE];  

  int IVAR_GAMMA   = INFO_DATA.TABLEVAR.ICUTWIN_GAMMA ;

  bool  DO_BIASCOR_SAMPLE = (INPUTS.opt_biasCor & MASK_BIASCOR_SAMPLE);
  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  bool  USE_FIELDGROUP    = (INPUTS.use_fieldGroup_biasCor > 0) ;

  int  IDSURVEY   = (int)TABLEVAR->IDSURVEY[ISN];
  int  SNTYPE     = (int)TABLEVAR->SNTYPE[ISN] ; 
  int  OPT_PHOTOZ = TABLEVAR->OPT_PHOTOZ[ISN];
  char *name      = TABLEVAR->name[ISN];
  float x0        = TABLEVAR->x0[ISN];
  float x0err     = TABLEVAR->x0err[ISN];
  float x1err     = TABLEVAR->fitpar_err[INDEX_x1][ISN];
  float cerr      = TABLEVAR->fitpar_err[INDEX_c][ISN];
  float zhd       = TABLEVAR->zhd[ISN];
  float zhderr    = TABLEVAR->zhderr[ISN];
  float vpec      = TABLEVAR->vpec[ISN];
  float vpecerr   = TABLEVAR->vpecerr[ISN];
  float COV_x0x1  = TABLEVAR->COV_x0x1[ISN];
  float COV_x0c   = TABLEVAR->COV_x0c[ISN];
  float COV_x1c   = TABLEVAR->COV_x1c[ISN];
  float SIM_X0    = TABLEVAR->SIM_X0[ISN];
  int   SIM_NONIA_INDEX = TABLEVAR->SIM_NONIA_INDEX[ISN];
  
  float mB, mBerr, mB_orig, mB_off, x1, c, sf;
  float zpec, zcmb, zMIN, zMAX, logmass ;
  double covmat8_fit[NLCPAR][NLCPAR], covmat8_int[NLCPAR][NLCPAR];
  double covmat8_tot[NLCPAR][NLCPAR], covtmp8 ;
  int   OPTMASK_COV, ISTAT_COV, i, i2, IDSAMPLE;
  char  *field, SURVEYGROUP[100], FIELDGROUP[100], STRINGOPT[40];
  char NONE[] = "NONE";
  char fnam[] =  "compute_more_TABLEVAR";

  // ----------------- BEGIN ------------------

  if ( INPUTS.cat_only ) { return; }

  if ( FIRST_EVENT ) {
    TABLEVAR->NCOVFIX  =  0 ;
  }

  // convert x0 and error to mB[err]
  x0 = TABLEVAR->x0[ISN];  x0err=TABLEVAR->x0err[ISN]; 
  mB = mBerr = -9.0; sf = 0.0 ;
  if ( x0 > 0.0 ) { 
    mB = -2.5*log10f(x0); 
    sf = -2.5/(x0*LOGTEN) ;
    mBerr = fabs(x0err * sf);
  }
  TABLEVAR->fitpar[INDEX_mB][ISN]     = mB ;  // this mB has no offset
  TABLEVAR->fitpar_err[INDEX_mB][ISN] = mBerr ;

  if ( IDSURVEY < 0 || IDSURVEY > MXIDSURVEY ) {
    sprintf(c1err,"Invalid IDSURVEY=%d (unsigned short) for %s SNID=%s", 
	    IDSURVEY, STRTYPE, name );
    sprintf(c2err,"Check $SNDATA_ROOT/SURVEY.DEF" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // increment IDSURVEY-dependent stuff
  TABLEVAR->NSN_PER_SURVEY[IDSURVEY]++ ;
  zMIN = TABLEVAR->zMIN_PER_SURVEY[IDSURVEY] ;
  zMAX = TABLEVAR->zMAX_PER_SURVEY[IDSURVEY] ;
  if ( zhd > 0.0 ) {
    if ( zhd < zMIN ) { TABLEVAR->zMIN_PER_SURVEY[IDSURVEY] = zhd; }
    if ( zhd > zMAX ) { TABLEVAR->zMAX_PER_SURVEY[IDSURVEY] = zhd; }
  }

  // - - - - covariance matrix - - - - - - - -  
  covmat8_fit[INDEX_mB][INDEX_mB] = (double)(x0err*x0err*sf*sf) ;
  covmat8_fit[INDEX_mB][INDEX_x1] = (double)(COV_x0x1 * sf);
  covmat8_fit[INDEX_mB][INDEX_c]  = (double)(COV_x0c  * sf) ;
  covmat8_fit[INDEX_x1][INDEX_x1] = (double)(x1err * x1err) ;
  covmat8_fit[INDEX_x1][INDEX_c]  = (double)(COV_x1c) ;
  covmat8_fit[INDEX_c][INDEX_c]   = (double)(cerr * cerr) ;

  // symmetric part of off-diagonals
  covmat8_fit[INDEX_x1][INDEX_mB] = (double)(COV_x0x1 * sf) ;
  covmat8_fit[INDEX_c][INDEX_mB]  = (double)(COV_x0c  * sf) ;
  covmat8_fit[INDEX_c][INDEX_x1]  = (double)(COV_x1c);

  OPTMASK_COV=0;
  if ( strcmp(name,INPUTS.SNID_MUCOVDUMP)==0 ) { OPTMASK_COV=4; } //dump COV
  update_covMatrix(name, OPTMASK_COV, NLCPAR, covmat8_fit, 
		   EIGMIN_COV, &ISTAT_COV ); 
  if ( ISTAT_COV != 0 ) { TABLEVAR->warnCov[ISN]=true;  TABLEVAR->NCOVFIX++; }

  // transfer local covmat8 to TABLEVAR float-storage
  for(i=0; i<NLCPAR; i++ ) {
    for(i2=0; i2<NLCPAR; i2++ )
      {  TABLEVAR->covmat_fit[ISN][i][i2] = (float)covmat8_fit[i][i2]; }
  }

  // estimate covmat_tot = covmat_fit + covmat_intrinsic
  get_COVINT_model(-1,covmat8_int); // just adds COV00 = sigmb^2
  for (i=0; i<NLCPAR; i++) {
    for (i2=0; i2<NLCPAR; i2++) {
      covtmp8 = covmat8_fit[i][i2] + covmat8_int[i][i2];
      covmat8_tot[i][i2]               = (float)covtmp8 ; // for local use
      TABLEVAR->covmat_tot[ISN][i][i2] = (float)covtmp8 ;
    }
  }

  // Aug 22 2019: logmass
  if ( INPUTS.USE_GAMMA0 && IVAR_GAMMA >= 0 )
    { logmass = (double)TABLEVAR->CUTVAL[IVAR_GAMMA][ISN]; }
  else
    { logmass = INPUTS.parval[IPAR_LOGMASS_CEN]; }

  TABLEVAR->logmass[ISN] = logmass;


  // - - - - - - - - - - - - - 
  // make sure that zHD is positive. Skip test for user-input zVARNAME.
  if ( zhd < 0 && strlen(INPUTS.varname_z)==0 ) { 
    sprintf(c1err,"Invalid %s redshift = %f ", STRTYPE, zhd );
    sprintf(c2err,"for ISN=%d, SNID='%s' ", ISN, name );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
  }

  // - - - - - - - 
  // Stuff for data only,
  if ( IS_DATA && zhd > 0 ) {

    // if user-input zpecerr > 0, subtract out original zpecerr and 
    // add user input zpecerr:
    // zhd = (1+zcmb)*(1+zpec) -1 = zcmb + zpec + zcmb*zpec
    // zcmb(1+zpec) = zhd - zpec     
    zpec    =  vpec/LIGHT_km ;
    zcmb    =  (zhd - zpec)/(1.0+zpec) ;
    TABLEVAR->zhderr[ISN] = zerr_adjust(zcmb, zhderr, vpecerr, name);

    // if cosmo params are fixed, then store mumodel for each event.
    if ( INPUTS.FLOAT_COSPAR == 0 ) {
      double z8 = (double)zhd;
      double dl = cosmodl_forFit(z8,INPUTS.COSPAR);
      TABLEVAR->mumodel[ISN] = (float)(5.0*log10(dl) + 25.0);
    }

    // check z-cheat option 
    if ( INPUTS.uzsim && TABLEVAR->IS_SIM == true ) { 
      TABLEVAR->zhd[ISN]     = TABLEVAR->SIM_ZCMB[ISN];
      TABLEVAR->zhderr[ISN]  = 1.0E-7;
    }


    if ( INPUTS.force_pIa >= 0.0 ) 
      { TABLEVAR->pIa[ISN] = INPUTS.force_pIa; } 

    if ( INPUTS.perfect_pIa )  {
      if ( !TABLEVAR->IS_SIM ) {
	sprintf(c1err,"Cannot force_pIa=perfect for real data.");
	sprintf(c2err,"This option works only for sim data.") ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      if ( SIM_NONIA_INDEX == 0 )
	{ TABLEVAR->pIa[ISN] = 1.0; } // it's a true SNIa
      else
	{ TABLEVAR->pIa[ISN] = 0.0; } // it's not SNIa
    } // end perfect_pIa


    // check option for z-dependent intrinsic COVMAT
    load_ZPOLY_COVMAT(IDSURVEY,zhd);

  } // end IS_DATA


  // check option to force pIa = 1 for spec confirmed SNIa
  if ( force_probcc0(SNTYPE,IDSURVEY) ) { TABLEVAR->pIa[ISN] = 1.0 ;  } 


  // - - - - - - - - - - - - - - - - - - - - - - 
  // IDSAMPLE is not ready for data yet,
  // but get it for BIASCOR and/or CCPRIOR

  if ( !IS_DATA ) {
    if ( USE_FIELDGROUP ) 
      { field = TABLEVAR->field[ISN]; }
    else
      { field = NONE ; }

    if ( DO_BIASCOR_SAMPLE) {
      match_fieldGroup (name, field,   
			FIELDGROUP, STRINGOPT);  // (O)
      match_surveyGroup(name, IDSURVEY, 
			SURVEYGROUP, STRINGOPT); // (O)
      IDSAMPLE = get_IDSAMPLE(IDSURVEY,OPT_PHOTOZ,FIELDGROUP,SURVEYGROUP,0) ;
    }
    else {
      IDSAMPLE = 0 ;   FIELDGROUP[0] = 0 ;
    }
    TABLEVAR->IDSAMPLE[ISN] = IDSAMPLE;

    // misc stuff for BIASCOR & CCPRIOR
    // increment number of events per IDSAMPLE (for biasCor)
    SAMPLE_BIASCOR[IDSAMPLE].NSN[EVENT_TYPE]++ ;

    // load sim_mb = -2.5log10(sim_x0)
    mB = -9.0;
    if ( SIM_X0 > 0.0 ) { mB = -2.5*log10(SIM_X0); }
    mB_orig = TABLEVAR->SIM_FITPAR[INDEX_mB][ISN];
    mB_off  = mB_orig - mB;
    TABLEVAR->SIM_FITPAR[INDEX_mB][ISN] = mB;  // no offset

    // legacy check on user zpecerr 
    if ( INPUTS.zpecerr > 1.0E-9 ) {
      double zpecerr = INPUTS.zpecerr * ( 1.0 + zhd ) ;
      TABLEVAR->zhderr[ISN]= (float)sqrt(zhderr*zhderr + zpecerr*zpecerr);
    }

  }  // end not-DATA

  // - - - - -
  if ( !IS_DATA  && DO_BIASCOR_MU ) { 
    // load mu for opton bias-correct MU instead of correcting mB,x1,c 
    // Note that true Alpha,Beta,GammaDM are used for mu_obs.
    // Beware that M0_DEFAULT may be fragile.
    double Alpha, Beta, GammaDM, mu_obs, mu_true; 
    mB      = TABLEVAR->fitpar[INDEX_mB][ISN] ;
    x1      = TABLEVAR->fitpar[INDEX_x1][ISN] ;
    c       = TABLEVAR->fitpar[INDEX_c][ISN] ;

    /* xxx NO 4/10/2020 xxxxx
    Alpha   = TABLEVAR->SIM_ALPHA[ISN] ;
    Beta    = TABLEVAR->SIM_BETA[ISN] ;
    xxxxxxx */

    Alpha   = INPUTS.parval[IPAR_ALPHA0]; // 4/10/2020: temp hack;later should
    Beta    = INPUTS.parval[IPAR_BETA0];  // measure A,B from slopes of HR vs x1,c
    GammaDM = TABLEVAR->SIM_GAMMADM[ISN] ;

    mu_true = TABLEVAR->SIM_MU[ISN] ;
    mu_obs  = mB - M0_DEFAULT + Alpha*x1 - Beta*c - GammaDM ;
    TABLEVAR->fitpar[INDEX_mu][ISN]      = mu_obs ; 
    TABLEVAR->SIM_FITPAR[INDEX_mu][ISN]  = mu_true ;
  }

  if ( IS_BIASCOR && IDEAL ) {
    float x0_ideal = TABLEVAR->x0_ideal[ISN];
    TABLEVAR->fitpar_ideal[INDEX_mB][ISN] = -2.5*log10f(x0_ideal); 
  }


  // - - - -  -
  // finally, set cutmask in TABLEVAR
  set_CUTMASK(ISN, TABLEVAR);


  // on last event, set NPASS ... but beware for data because additional
  // DATA cuts are applied later. print_eventStat determines final NPASS.
  int NALL   = *NALL_CUTMASK_POINTER[EVENT_TYPE];  
  int NREJ   = *NREJECT_CUTMASK_POINTER[EVENT_TYPE];  
  int *NPASS = NPASS_CUTMASK_POINTER[EVENT_TYPE];
  if ( ISN == NALL-1 ) {     *NPASS = NALL-NREJ ;   }
  

  return ;

} // end compute_more_TABLEVAR

// ==================================================
void compute_more_INFO_DATA(void) {

  // Created Jun 2019
  // set initial muerr so that optional log(1/sigma) term in fcn-chi2
  // can be fixed on first iteration

  int NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL ;
  double alpha      = INPUTS.parval[IPAR_ALPHA0] ;
  double beta       = INPUTS.parval[IPAR_BETA0] ;
  double gamma      = INPUTS.parval[IPAR_GAMMA0] ;
  double *ptr_sigCC = &INPUTS.parval[IPAR_H11+3];
  double muerrsq, sigCC, zhd, zhderr, cov, covmat_tot[NLCPAR][NLCPAR] ;
  int    isn, CUTMASK, i, i2 ;
  char *name;
  //  char fnam[] = "compute_more_INFO_DATA";

  // ----------- BEGIN ------------

  if ( INPUTS.cat_only) { return; }

  for(isn=0; isn < NSN_DATA; isn++ ) {

    muerrsq = 9.0;
    CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[isn];
    name    = INFO_DATA.TABLEVAR.name[isn];
    if ( CUTMASK ) 
      { muerrsq = 9.0; sigCC = 9.0; }
    else {
      zhd    = (double)INFO_DATA.TABLEVAR.zhd[isn];
      zhderr = (double)INFO_DATA.TABLEVAR.zhderr[isn];
      for(i=0; i < NLCPAR; i++ ) {
	for(i2=0; i2 < NLCPAR; i2++ ) {
	  cov = (double)INFO_DATA.TABLEVAR.covmat_tot[isn][i][i2];
	  covmat_tot[i][i2] = cov;
	}
      }
      muerrsq = fcn_muerrsq(name,alpha,beta,gamma, covmat_tot,zhd,zhderr, 0);
      sigCC   = ptr_sigCC[0] + zhd*ptr_sigCC[1] + zhd*zhd*ptr_sigCC[2];
    }

    INFO_DATA.muerrsq_last[isn]  = muerrsq ;
    INFO_DATA.muerr_last[isn]    = sqrt(muerrsq);
    INFO_DATA.sigCC_last[isn]    = sigCC ;
    INFO_DATA.sqsigCC_last[isn]  = sigCC * sigCC ;
  }

  return;

} // end compute_more_INFO_DATA


// ==================================================
void  store_input_varnames(int ifile, TABLEVAR_DEF *TABLEVAR) {

  // Created May 8 2020 by R.Kessler
  // Check for missing columns and decide whether or not to abort.
  // Main goal is to allow missing classifier PROB columns for 
  // spec-confirmed samples like LOWZ. However, can also be used
  // if a subset of samples have extra columns from APPEND_FITRES.

  char *INPFILE     = TABLEVAR->INPUT_FILE[ifile];
  int   FIRST_EVENT = TABLEVAR->EVENT_RANGE[ifile][0];
  int   LAST_EVENT  = TABLEVAR->EVENT_RANGE[ifile][1];
  int   IVAR_pIa    = TABLEVAR->IVAR_pIa[ifile];
  char *varname_pIa = INPUTS.varname_pIa ;
  bool USE_PROBCC_ZERO = (INPUTS_PROBCC_ZERO.nidsurvey > 0 );

  int evt, sntype, idsurvey, NFORCE, NOTFORCE ;
  char fnam[] = "store_input_varnames" ;

  // ------------ BEGIN --------------

  /*
  printf(" xxx %s hello, IVAR_pIa=%d  N_PCC_ZERO=%d\n",
	 fnam, IVAR_pIa, INPUTS_PROBCC_ZERO.nidsurvey );
  */

  if ( strlen(varname_pIa) > 0 ) {

    // variable for probIa is defined.
    // default is that no events have forced pIa
    NFORCE = 0;
    NOTFORCE = (LAST_EVENT - FIRST_EVENT + 1);

    // if force_probcc is set, explicitly set pIa=1 in this file
    // regardless of whether this column was read from the input file.
    // Thus the varname_pIa column need not exist in a file as long 
    // as all events in the file have forced pIa=1.
    if ( USE_PROBCC_ZERO ) {
      NFORCE = NOTFORCE = 0;
      for(evt=FIRST_EVENT; evt < LAST_EVENT; evt++ ) {
	sntype      = (int)TABLEVAR->SNTYPE[evt] ; 
	idsurvey    = (int)TABLEVAR->IDSURVEY[evt];
	if ( force_probcc0(sntype,idsurvey) ) 
	  { NFORCE++; }
	else
	  { NOTFORCE++ ; }
      }
      printf("\t Force pIa=1.0 for %d events from %s\n",
	     NFORCE, INPFILE); fflush(stdout);
    }

    // abort if varname_pIa column is missing and there are events
    // which don't have forced pIa=1.
    if (  IVAR_pIa < 0 && NOTFORCE > 0 ) {
      print_preAbort_banner(fnam);
      printf("\t input file: %s \n", INPFILE);
      printf("\t has %d events with forced pIa=1\n", NFORCE);
      printf("\t and %d unforced-pIa events.\n", NOTFORCE);
      sprintf(c1err,"Missing column varname_pIa='%s' .",  varname_pIa);
      sprintf(c2err,"Either add this column, or check type_list_probcc0");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
    }

    
  } // end varname_pIa


  // - - - - - - - - - - -
  // store varnames for each input file; used later to
  // prepare/align output table columns.
  int NVAR =  READTABLE_POINTERS.NVAR_TOT ;
  int MEMC1 =  NVAR * sizeof(char*) ;
  int MEMC0 =  MXCHAR_VARNAME * sizeof(char);
  int ivar ;

  TABLEVAR->VARNAMES_LIST[ifile] = (char**) malloc(MEMC1);
  TABLEVAR->NVAR[ifile]          = NVAR;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    TABLEVAR->VARNAMES_LIST[ifile][ivar] = (char*) malloc(MEMC0);
    sprintf(TABLEVAR->VARNAMES_LIST[ifile][ivar], "%s", 
	    READTABLE_POINTERS.VARNAME[ivar] );
  }

  return ;

} // end  store_input_varnames


// ==================================================
void store_output_varnames(void) {

  // Created May 2020
  // Analyze columns for each input data file and prepare
  // final varnames list for output fitres file, and also
  // a map from final column (ivar) location back to origina
  // column (ivar) in each input data file.
  // 
  // Here is the logic for keeping output variables:
  //  + keep all variables appearing in every data file
  //  + for variables in subset of data files, keep if
  //    on the append_missing_varnames list, or varname_pIa;
  //    otherwise drop from output.
  //
  // Notes: 
  // + SALT2mu/BBC info is NOT included here; this function works 
  //    only with the input VARNAMES.
  // + works on data only since only data are written out.
  // + Computed info stored in OUTPUT_VARNAMES struct.
  
  int  NFILE = INPUTS.nfile_data;
  //  int  MEMC  = MXCHAR_VARNAME * sizeof(char);
  int  ifile, ifile2, NVAR, NVAR2, MXVAR, NVAR_TOT; 
  int  ivar_match, ivar, ivar2, NMATCH, i ;
  char *varName, *varName2 ;
  bool wildcard, MATCH ;
  int LDMP = 0 ;
  char fnam[] = "store_output_varnames" ;
  // ------------- BEGIN ------------
  
  // xxx  OUTPUT_VARNAMES.PERFECT_COLUMN_MATCH = true;

  // alloate VARNAMES memory for x2 number of variables in first file ...
  // should be enough
  MXVAR = 2*INFO_DATA.TABLEVAR.NVAR[0];
  for(ivar=0; ivar < MXVAR; ivar++ ) {
    OUTPUT_VARNAMES.LIST[ivar] = (char*) malloc(MXVAR*sizeof(char) ) ;
    for(ifile=0; ifile < NFILE; ifile++ ) 
      { OUTPUT_VARNAMES.IVARMAP[ifile][ivar] = -888 ; }
  }

	    

  // check trivial case with just one input file
  if ( NFILE == 1 ) {
    ifile = 0 ;
    NVAR = INFO_DATA.TABLEVAR.NVAR[ifile];
    OUTPUT_VARNAMES.NVAR_TOT = NVAR;
    for(ivar=0; ivar < NVAR; ivar++ ) { 
      varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile][ivar];
      sprintf(OUTPUT_VARNAMES.LIST[ivar],"%s", varName);
      OUTPUT_VARNAMES.IVARMAP[ifile][ivar] = ivar ; 
    }
    return ;
  }

  // - - - - - - 
  // here we have multiple input data files, so check them all.
  NVAR_TOT = 0 ;

  // start by collecting columns that appear in every data file.
  // Loop over variables in first file and check if var exists
  // in the other files
  NVAR = INFO_DATA.TABLEVAR.NVAR[0];
  for(ivar=0; ivar < NVAR; ivar++ ) { 
    varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[0][ivar];
    NMATCH  = 0 ; 
    for(ifile=0; ifile < NFILE; ifile++ ) {
      NVAR2 = INFO_DATA.TABLEVAR.NVAR[ifile] ;
      ivar2 = ivar_matchList(varName, NVAR2, 
			     INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile]);
      OUTPUT_VARNAMES.IVARMAP[ifile][NVAR_TOT] = ivar2 ; 
      if ( ivar2 >= 0 ) { NMATCH++ ; }
      /*
      printf(" xxx %s: ivar0=%2d(%s) ->  ivar[ifile=%d]=%d \n",
	     fnam, ivar, varName,  ifile, ivar2); fflush(stdout);
      */
    }

    // if this varName is in all files, add it to final list
    if( NMATCH == NFILE ) {
      if ( LDMP ) {
	printf(" xxx %s: all files VAR[%2d] = %s \n",
	       fnam, NVAR_TOT, varName); fflush(stdout);
      }
      sprintf(OUTPUT_VARNAMES.LIST[NVAR_TOT],"%s", varName);
      NVAR_TOT++;
    }
  }  // end ivar loop
  

  // - - - - - -
  // next, check for variables in a subset of files
  for(i=0; i < INPUTS_VARNAME_MISSING.ndef; i++ ) {
    wildcard = INPUTS_VARNAME_MISSING.wildcard[i] ;
    varName  = INPUTS_VARNAME_MISSING.varname_list[i] ;
    
    if ( LDMP ) {
      printf(" xxx ----------------------------------- \n");
      printf(" xxx check for varName = %s  (wildcard=%d) \n",
	     varName, wildcard );
    }
    // skip if already included
    for(ifile=0; ifile < NFILE; ifile++ ) {
      NVAR = INFO_DATA.TABLEVAR.NVAR[ifile];
      for(ivar2=0; ivar2 < NVAR; ivar2++ ) {
	varName2 = INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile][ivar2] ;
	ivar_match = ivar_matchList(varName2, NVAR_TOT, 
				    OUTPUT_VARNAMES.LIST );
	if ( ivar_match >= 0 ) { continue; }
	if ( wildcard ) 
	  { MATCH = (strstr(varName2,varName) != NULL ); }
	else
	  { MATCH = (strcmp(varName,varName2) == 0 ); }

	if ( MATCH ) {
	  sprintf(OUTPUT_VARNAMES.LIST[NVAR_TOT],"%s", varName2);
	  for(ifile2=0; ifile2 < NFILE; ifile2++ ) {
	    NVAR2 = INFO_DATA.TABLEVAR.NVAR[ifile2];
	    OUTPUT_VARNAMES.IVARMAP[ifile2][NVAR_TOT] = 
	      ivar_matchList(varName2, NVAR2, 
			     INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile2] );
	  }
	  if ( LDMP ) {
	    printf(" xxx %s: append VAR[%2d] = %s \n", 
		   fnam, NVAR_TOT, varName2  );
	  }
	  NVAR_TOT++ ;
	}
	/*
	printf(" xxx ifile=%d  ivar2=%2d  check varName2 = '%s'  MATCH=%d \n", 
	       ifile, ivar, varName2, MATCH );
	*/
      }
    }
			  
  } // end i loop over append_missing
  
  OUTPUT_VARNAMES.NVAR_TOT = NVAR_TOT ;

  // build list of dropped variables, for comment use only
  if ( LDMP ) {
    printf(" xxx ----------------------------------- \n");
    printf(" xxx check for dropped table columns: \n");
  }

  int NDROP=0 ;
  OUTPUT_VARNAMES.DROPLIST[0] = 0;
  for(ifile=0; ifile < NFILE; ifile++ ) {
    NVAR = INFO_DATA.TABLEVAR.NVAR[ifile];
    for(ivar=0; ivar < NVAR; ivar++ ) {
      varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile][ivar] ;
      ivar_match = ivar_matchList(varName, NVAR_TOT, 
				  OUTPUT_VARNAMES.LIST );
      if ( ivar_match < 0  ){
	NDROP++ ;
	strcat(OUTPUT_VARNAMES.DROPLIST,varName);
	strcat(OUTPUT_VARNAMES.DROPLIST," ");	
	if ( LDMP) {
	  printf(" xxx %s: DROP varName='%s'  match=%d (ifile=%d,ivar=%d)\n",
		 fnam,  varName, ivar_match, ifile, ivar );
	}
      }
    }
  }
  OUTPUT_VARNAMES.NVAR_DROP = NDROP ;

  if ( LDMP )     { printf(" xxx %s DONE. \n", fnam); }

  //  debugexit(fnam);

  return ;

} // end store_output_varnames

// ============
bool exist_varname(int ifile, char *varName, TABLEVAR_DEF *TABLEVAR) {
  int NVAR = TABLEVAR->NVAR[ifile];
  int ivar;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    if ( strcmp(varName,TABLEVAR->VARNAMES_LIST[ifile][ivar]) == 0 )
      { return true; }
  }
  return false;
} // end exist_varname

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
  // Examine data to determine IDSAMPLE for each data event.
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
  //
  // Jun 6 2019: refactor to use INFO_DATA struct.
  //
  // May 29 2020: 
  //   + fix bug for ALL surveys lumped together; set all IDSAMPLE=0
  //

  int USE_FIELDGROUP  = INPUTS.use_fieldGroup_biasCor ;
  int isn, IDSURVEY, OPT_PHOTOZ, N, IDSAMPLE, i, NIDSURVEY[MXIDSURVEY] ;
  int  DUMPFLAG=0, NDMP = 0, NSN_DATA, CUTMASK  ; 
  double z;
  char FIELD_TMP[MXCHAR_CCID],  FIELDGROUP[100],  *FIELDDEF=NULL;
  char SURVEYGROUP[100], SURVEYDEF[MXCHAR_CCID], zGROUP[20];
  char *NAME_SN,  *NAME_SAMPLE, STRINGOPT[40]  ;
  char fnam[] = "prepare_IDSAMPLE_biasCor" ;

  // --------------- BEGIN -----------------

  NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL; 

  NSAMPLE_BIASCOR      = 0 ;
  ONE_SAMPLE_BIASCOR   = 0 ;
  INPUTS_SAMPLE_BIASCOR.NFIELDGROUP_USR  = 0 ;
  INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR = 0 ;
  INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_TOT = 0 ;

  if ( INPUTS.opt_biasCor == 0 ) { return ; }

  sprintf(BANNER,"Begin %s", fnam);
  print_banner(BANNER);

      
  for(i=0; i < MXNUM_SAMPLE; i++ ) { 
    SAMPLE_BIASCOR[i].NSN[EVENT_TYPE_DATA]     = 0 ; 
    SAMPLE_BIASCOR[i].NSN[EVENT_TYPE_BIASCOR]  = 0 ; 
    SAMPLE_BIASCOR[i].NSN[EVENT_TYPE_CCPRIOR]  = 0 ; 
    SAMPLE_BIASCOR[i].DOFLAG_SELECT  = 1 ; // default -> select sample
    SAMPLE_BIASCOR[i].DOFLAG_BIASCOR = 1 ; // default -> do Bias cor
    SAMPLE_BIASCOR[i].IFLAG_ORIGIN   = 0 ;
    SAMPLE_BIASCOR[i].zMIN_DATA = +999. ;
    SAMPLE_BIASCOR[i].zMAX_DATA = -999. ;
    SAMPLE_BIASCOR[i].NAME_SURVEYGROUP[0] = 0 ;
    SAMPLE_BIASCOR[i].NAME_FIELDGROUP[0]  = 0 ;
    SAMPLE_BIASCOR[i].STRINGOPT[0]        = 0 ;
    SAMPLE_BIASCOR[i].OPT_PHOTOZ          = 0 ; // zSPEC is default
    SAMPLE_BIASCOR[i].NAME[0]             = 0 ;
    SAMPLE_BIASCOR[i].IDSURVEY            = -9 ;
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

    SAMPLE_BIASCOR[IDSAMPLE].NSN[EVENT_TYPE_DATA] = NSN_DATA;
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME_FIELDGROUP,  "NONE" );
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME_SURVEYGROUP, "ALL"  );
    sprintf(SAMPLE_BIASCOR[IDSAMPLE].NAME,        "ALL" );
    SAMPLE_BIASCOR[IDSAMPLE].OPT_PHOTOZ   = 0 ; // zSPEC is default

    // set redshift info (Apr 4 2017)
    SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA = INPUTS.zmin ;
    SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA = INPUTS.zmax ;
    set_BINSIZE_SAMPLE_biasCor(IDSAMPLE);

    for(isn=0; isn < NSN_DATA; isn++ )
      { INFO_DATA.TABLEVAR.IDSAMPLE[isn] = IDSAMPLE ;} // May 29 2020

    dump_SAMPLE_INFO(EVENT_TYPE_DATA);
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


  for(isn=0; isn < NSN_DATA; isn++ ) {

    CUTMASK    = INFO_DATA.TABLEVAR.CUTMASK[isn];
    IDSURVEY   = INFO_DATA.TABLEVAR.IDSURVEY[isn];
    OPT_PHOTOZ = INFO_DATA.TABLEVAR.OPT_PHOTOZ[isn];
    NAME_SN    = INFO_DATA.TABLEVAR.name[isn];
    if(USE_FIELDGROUP) { FIELDDEF = INFO_DATA.TABLEVAR.field[isn]; }
    z          = INFO_DATA.TABLEVAR.zhd[isn];

    if ( CUTMASK ) { continue ; }

    if ( IDSURVEY < 0 || IDSURVEY > MXIDSURVEY ) {
      sprintf(c1err,"Invalid IDSURVEY=%d for SNID=%s", IDSURVEY, NAME_SN);
      sprintf(c2err,"Check $SNDATA_ROOT/SURVEY.DEF" );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // strip off stuff into local var

    sprintf(SURVEYDEF,"%s", SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY] );
    NIDSURVEY[IDSURVEY]++ ;

    // check which FIELDGROUP group this event is in 
    if ( USE_FIELDGROUP )  { sprintf(FIELD_TMP,"%s",  FIELDDEF ); }


    // check if this IDSURVEY and FIELDGROUP have been defined.
    // Note that only one of the match_xxx functions can return STRINGOPT
    STRINGOPT[0] = 0 ;
    match_fieldGroup (NAME_SN, FIELD_TMP,       // (I)
		      FIELDGROUP, STRINGOPT  ); // (O)

    if ( strlen(STRINGOPT) == 0 ) {
      match_surveyGroup(NAME_SN, IDSURVEY,        // (I)
			SURVEYGROUP, STRINGOPT ); // (O)
    }	

    if ( OPT_PHOTOZ == 0 ) 
      { sprintf(zGROUP,"-zSPEC"); }
    else
      { sprintf(zGROUP,"-zPHOT"); }
    if ( FOUNDKEY_OPT_PHOTOZ == 0 ) { zGROUP[0] = 0 ; }

    DUMPFLAG = 0 ;
    IDSAMPLE = get_IDSAMPLE(IDSURVEY, OPT_PHOTOZ, 
			    FIELDGROUP, SURVEYGROUP, DUMPFLAG);

    if ( IDSAMPLE < 0 ) {
      // store new SURVEY/FIELDGROUP entry
      N = NSAMPLE_BIASCOR ;  
      sprintf(SAMPLE_BIASCOR[N].NAME_FIELDGROUP,  "%s", FIELDGROUP  );
      sprintf(SAMPLE_BIASCOR[N].NAME_SURVEYGROUP, "%s", SURVEYGROUP ); 
      sprintf(SAMPLE_BIASCOR[N].STRINGOPT,        "%s", STRINGOPT   ); 
      SAMPLE_BIASCOR[N].OPT_PHOTOZ = OPT_PHOTOZ ;
      SAMPLE_BIASCOR[N].IDSURVEY = IDSURVEY ; // Jun 2020

      NAME_SAMPLE = SAMPLE_BIASCOR[N].NAME ;
      if ( IGNOREFILE(FIELDGROUP) ) { 
	sprintf(NAME_SAMPLE,"%s%s", SURVEYGROUP, zGROUP ); 
	SAMPLE_BIASCOR[N].IFLAG_ORIGIN = USERFLAG_SURVEYGROUP_SAMPLE ;
	// warning: this IFLAG_ORGIN does not distinguish USER and AUTO
      }
      else { 
	sprintf(NAME_SAMPLE,"%s%s(%s)", SURVEYGROUP, zGROUP, FIELDGROUP ); 
	SAMPLE_BIASCOR[N].IFLAG_ORIGIN = USERFLAG_FIELDGROUP_SAMPLE ;
      }

      NSAMPLE_BIASCOR++ ;

      N = NSAMPLE_BIASCOR ;
      if ( N >= MXNUM_SAMPLE ) {
	print_preAbort_banner(fnam);
	dump_SAMPLE_INFO(EVENT_TYPE_DATA);
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
      print_preAbort_banner(fnam);
      printf("   Current SURVEYGROUP = '%s' \n", SURVEYGROUP );
      printf("   Current FIELDGROUP  = '%s' \n", FIELDGROUP  );
      dump_SAMPLE_INFO(EVENT_TYPE_DATA);
      sprintf(c1err,"Final IDSAMPLE=%d for CID=%s(isn=%d)", 
	      IDSAMPLE, NAME_SN, isn );
      sprintf(c2err,"SURVEY=%s(%d)  FIELD=%s", 
	      SURVEYDEF, IDSURVEY, FIELDDEF );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    INFO_DATA.TABLEVAR.IDSAMPLE[isn] = IDSAMPLE ; 

    SAMPLE_BIASCOR[IDSAMPLE].NSN[EVENT_TYPE_DATA]++ ;

    // keep track of min/max redshift for each sample

    if( z < SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA = z; }
    if( z > SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA = z; }

    // check user option to skip biasCor for this sample
    if ( strstr(INPUTS.surveyList_noBiasCor,SURVEYDEF) != NULL ) 
      { SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR = 0 ; }

    // xxxxxxxxxxxxx
    if ( isn < NDMP ) {
      printf(" xxx '%8s'  '%8s'   %8.8s(%2d)   %2d    %s \n", 
	     NAME_SN, FIELDDEF, SURVEYDEF,IDSURVEY,  IDSAMPLE, FIELDGROUP );
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
  dump_SAMPLE_INFO(EVENT_TYPE_DATA);

  // check option to fix sigint for each IDSAMPLE
  parse_sigint_fix(INPUTS.sigint_fix);
   

  return ;

} // end prepare_IDSAMPLE_biasCor




// ======================================================
void  set_BINSIZE_SAMPLE_biasCor(int IDSAMPLE) {

  // Created Aug 2016
  // Decode STRINGOPT and set biasCor binSize for z,x1,c

  int ipar, nb;
  double dif;
  char fnam[] = "set_BINSIZE_SAMPLE_biasCor" ;

  // ------------- BEGIN -----------
  
  // first set defaults
  dif = INPUTS.zmax - INPUTS.zmin; nb=INPUTS.nzbin;
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[0] = INPUTS.zmin ;
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[1] = INPUTS.zmax ;
  SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT  = dif / (double)nb ;

  // logmass 
  dif = INPUTS.logmass_max - INPUTS.logmass_min; nb=INPUTS.nbin_logmass ;
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_LOGMASS[0] =  INPUTS.logmass_min ;
  SAMPLE_BIASCOR[IDSAMPLE].RANGE_LOGMASS[1] =  INPUTS.logmass_max ;
  SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_LOGMASS  =  dif / (double)nb ;

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
    zMIN = SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA ; // true zMIN,zMAX
    zMAX = SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA ; // in data sample

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
  // May 29 2020: 
  //  + abort if NO data exists for survey in user-input surveygroup_biascor
  //

  int  LDMP   = 0 ;
  int  USE_SURVEYGROUP = INPUTS.use_surveyGroup_biasCor ;
  int  i, i2, NGRP, ID, NEVT ;
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

      // May 29 2020: abort if there is no data in this requested sample
      NEVT = INFO_DATA.TABLEVAR.NSN_PER_SURVEY[ID];
      if ( NEVT == 0 ) {
        sprintf(c1err,"No data for requested SURVEY = %s(%d)", ptrTmp[i2],ID);
        sprintf(c2err,"Check input key surveygroup_biascor");
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

    NEVT = INFO_DATA.TABLEVAR.NSN_PER_SURVEY[ID];

    if ( NEVT == 0 )
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
  int ID, ID2, IDSURVEY, NSORT=0  ;
  int NCOPY[MXNUM_SAMPLE], IDMAP_SORT[MXNUM_SAMPLE];
  int igrp, NGRP_USR,  FLAG ; 
  char *s0, *s1, *s2, *NAME ;
  SAMPLE_INFO_DEF *SAMPLE_BIASCOR_TEMP ;
  int LDMP = 0 ;
  char fnam[] = "sort_IDSAMPLE_biasCor" ;
  // ---------------- BEGIN ---------------

  SAMPLE_BIASCOR_TEMP = 
    (SAMPLE_INFO_DEF*) malloc ( NSAMPLE * sizeof(SAMPLE_INFO_DEF));
  

  // Jan 2020: abort if some SURVEY events are part of a FIELDGROUP,
  //           and some events are not.
  //   e..g, if DES and DES(C3+X3) are defined, ABORT.
  bool ISGRP0, ISGRP2, SMATCH ;
  char *f0, *f2;
  for(ID=0; ID < NSAMPLE; ID++ ) {
    for(ID2=0; ID2<NSAMPLE; ID2++ ) {
      if ( ID == ID2 ) { continue ; } 
      f0     = SAMPLE_BIASCOR[ID].NAME_FIELDGROUP ;
      f2     = SAMPLE_BIASCOR[ID2].NAME_FIELDGROUP ;
      s0     = SAMPLE_BIASCOR[ID].NAME_SURVEYGROUP ;  
      s2     = SAMPLE_BIASCOR[ID2].NAME_SURVEYGROUP ;  
      SMATCH = strcmp(s0,s2) == 0 ;
      ISGRP0 = !IGNOREFILE(f0);      ISGRP2 = !IGNOREFILE(f2);

      if ( !ISGRP0 && ISGRP2 && SMATCH ) {
	sprintf(c1err,"Found %s events NOT in FIELDGROUP", s0);
        sprintf(c2err,"Define all fields in fieldgroup_biascor key.");
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

      }

      /*
      printf(" xxx ------------------------------ \n");
      printf(" xxx %s: %d,%d\n", fnam, ID, ID2 );
      printf(" xxx    SURVEYGROUP = '%s' and %s' \n",  s0, s2);
      printf(" xxx    FIELDGROUP  = %s(%d) and %s(%d) \n", 
	     f0,ISGRP0, f2,ISGRP2 );
      printf(" xxx    IFLAG_ORIGIN -> %d \n", 
      SAMPLE_BIASCOR[ID].IFLAG_ORIGIN ); */

    }
  }


  // copy pre-sorted struct into temp struct
  ID2 = 0;
  for(ID=0; ID < NSAMPLE; ID++ ) {
    if ( SAMPLE_BIASCOR[ID].IFLAG_ORIGIN == USERFLAG_IGNORE_SAMPLE ) 
      { continue; }       
    copy_IDSAMPLE_biasCor(&SAMPLE_BIASCOR[ID], &SAMPLE_BIASCOR_TEMP[ID2]);
    NCOPY[ID2] = 0; 
    ID2++ ;
  }
  NSAMPLE_BIASCOR = NSAMPLE = ID2;


  if ( LDMP ) {
    printf("\n xxx ---------- Start DUMP for %s ---------------- \n", fnam);
    dump_SAMPLE_INFO(EVENT_TYPE_DATA);
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

  for(isn=0; isn < INFO_DATA.TABLEVAR.NSN_ALL ; isn++ ) {
    ID = INFO_DATA.TABLEVAR.IDSAMPLE[isn] ;
    if ( ID < 0 ) { continue ; }
    INFO_DATA.TABLEVAR.IDSAMPLE[isn] = IDMAP_SORT[ID]; //sorted sample
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


  S2->NSN[EVENT_TYPE_DATA]    = S1->NSN[EVENT_TYPE_DATA] ;
  S2->NSN[EVENT_TYPE_BIASCOR] = S1->NSN[EVENT_TYPE_BIASCOR] ;
  S2->NSN[EVENT_TYPE_CCPRIOR] = S1->NSN[EVENT_TYPE_CCPRIOR] ;

  S2->zMIN_DATA = S1->zMIN_DATA ;
  S2->zMAX_DATA = S1->zMAX_DATA ;

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
void dump_SAMPLE_INFO(int EVENT_TYPE) {

  // Dump SAMPLE_INFO structure.
  // Input EVENT_TYPE --> "DATA" or "BIASCOR"  or  "CCPRIOR"
  //
  // Mar 22 2018: 
  //  + for BIASCOR or CCPRIOR, abort if any sample has zero event.
  //  + allow some CCPRIOR samples to have NEVT=0 (e.g., spec subsets)
  //
  // Jun 6 2019: change argument from *what to integer event type

  char *STRTYPE    = STRING_EVENT_TYPE[EVENT_TYPE];
  int  IS_DATA     = (EVENT_TYPE == EVENT_TYPE_DATA);
  int  IS_BIASCOR  = (EVENT_TYPE == EVENT_TYPE_BIASCOR);
  int  N           = NSAMPLE_BIASCOR;

  int  NOBIASCOR, NEVT, i;
  int  NZERR=0, NSAMPLE_ZERO=0, ABORT_ON_NEVTZERO=0;
  char *NAME, Nstring[60], zString[60] ;
  char fnam[] = "dump_SAMPLE_INFO" ;

  // -------- BEGIN ----------

  printf("  SAMPLE_INFO DUMP for %s: \n", STRTYPE )    ;

  if ( IS_BIASCOR ) { ABORT_ON_NEVTZERO=1; }
  for(i=0; i < N; i++ )  { 

    NEVT = SAMPLE_BIASCOR[i].NSN[EVENT_TYPE] ; 

    NAME      = SAMPLE_BIASCOR[i].NAME ;
    NOBIASCOR = ( SAMPLE_BIASCOR[i].DOFLAG_BIASCOR == 0 );
    if ( NEVT == 0 ) { NSAMPLE_ZERO++; }

    sprintf(Nstring,"N%-7.7s=%6d", STRTYPE, NEVT);
    if ( NOBIASCOR && !IS_DATA ) { sprintf(Nstring,"No %s", STRTYPE); }

    sprintf(zString,"%.3f<z<%.3f", 
	    SAMPLE_BIASCOR[i].zMIN_DATA, SAMPLE_BIASCOR[i].zMAX_DATA );

    if( SAMPLE_BIASCOR[i].zMAX_DATA < SAMPLE_BIASCOR[i].zMIN_DATA ) { NZERR++; }
    // print one-line summary per SAMPLE
    printf("  IDSAMPLE=%2d --> %-20.20s  (%s, %s)\n",
	   i, NAME, Nstring, zString );

    //.xyz printf(" xxx %s: IDSURVEY = %d \n", fnam, SAMPLE_BIASCOR[i].IDSURVEY);
  }

  if ( NZERR > 0 ) {
    sprintf(c1err,"zMAX < zMIN for %d IDSAMPLEs", NZERR);
    sprintf(c2err,"See IDSAMPLE dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( ABORT_ON_NEVTZERO && NSAMPLE_ZERO ) {
    sprintf(c1err,"%d of %d SAMPLES have zero events.", NSAMPLE_ZERO,N);
    sprintf(c2err,"Check %s sample.", STRTYPE );
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
    print_preAbort_banner(fnam);
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
    print_preAbort_banner(fnam);
    printf("  NSURVEYGROUP = %d \n", NSURVEYGROUP );
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
    print_preAbort_banner(fnam);
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
      print_preAbort_banner(fnam);
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
void prepare_biasCor(void) {

  // Created Jan 2016
  // If simFile is not null, then prepare mB bias correction vs. z
  // for SNIa.
  //  -> read simFile
  //  -> make map of mB-SIM_mB vs redshift
  //
  // Jun 2 2016: call sigInt_biasCor BEFORE makeMape_fitPar
  // Jan 16 2018: force correct value for INPUTS.fitflag_sigmb = 2 
  //
  // Oct 24 2019: 
  //  + refactor to define NDIM_BIASCOR and global INFO_BIASCOR.NDIM.
  //    Print more sensible messages with appropriate 1D,5D,6D,7D.
  //
  // Dec 11 2019: fix to work with DOCOR_1D5DCUT
  //
  // Feb 5 2020: if > 5 SIM_gammaDM bins, do NOT add another dimension

  int INDX, IDSAMPLE, SKIP, NSN_DATA, CUTMASK, ievt ;
  int  NBINm = INPUTS.nbin_logmass;
  int  NBINg = 0; // is set below 
  int  OPTMASK        = INPUTS.opt_biasCor ;
  int  EVENT_TYPE     = EVENT_TYPE_BIASCOR;
  int  nfile_biasCor  = INPUTS.nfile_biasCor ;
  bool  DOCOR_1DZAVG  = ( OPTMASK & MASK_BIASCOR_1DZAVG  );
  bool  DOCOR_1DZWGT  = ( OPTMASK & MASK_BIASCOR_1DZWGT  );
  bool  DOCOR_1D5DCUT = ( OPTMASK & MASK_BIASCOR_1D5DCUT );
  bool  DOCOR_MUCOV   = ( OPTMASK & MASK_BIASCOR_MUCOV   );
  bool  DOCOR_5D      = ( OPTMASK & MASK_BIASCOR_5D      );
  bool  DOCOR_1D      = ( DOCOR_1DZAVG || DOCOR_1DZWGT || DOCOR_1D5DCUT);
  bool  IDEAL         = ( OPTMASK & MASK_BIASCOR_COVINT ) ;
  bool  DOCOR_MU      = ( OPTMASK & MASK_BIASCOR_MU ) ;
  char txt_biasCor[40]  ;
  
  bool USEDIM_GAMMADM, USEDIM_LOGMASS;
  int NDIM_BIASCOR=0, ILCPAR_MIN, ILCPAR_MAX ;
  char fnam[] = "prepare_biasCor" ;

  // ------------- BEGIN -------------
  
  INFO_BIASCOR.NDIM = 0;
  INFO_BIASCOR.TABLEVAR.NSN_ALL       = 0 ;
  INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS  = 0 ;
  INFO_BIASCOR.TABLEVAR.NSN_REJECT    = 0 ;
  INFO_BIASCOR.GAMMADM_OFFSET         = 0.0 ;
  NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;

  if ( DOCOR_MU ) 
    { INFO_BIASCOR.ILCPAR_MIN = INDEX_mu; INFO_BIASCOR.ILCPAR_MAX = INDEX_mu;}
  else
    { INFO_BIASCOR.ILCPAR_MIN = 0; INFO_BIASCOR.ILCPAR_MAX = NLCPAR-1; }
  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;

  if ( nfile_biasCor == 0 ) { return ; }

  print_banner(fnam);

  // read biasCor file
  read_simFile_biasCor();

  // setup 5D bins 
  if ( DOCOR_5D || DOCOR_1D5DCUT ) {
    for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ ) 
      { setup_CELLINFO_biasCor(IDSAMPLE); }
  }


  if ( DOCOR_5D ) {  
    NDIM_BIASCOR = 5;  
    sprintf(INFO_BIASCOR.STRING_PARLIST,"z,x1,c,a,b");

    // if there are 2 gammaDM bins, automatically add this dimension
    // to biasCor; if more than 5 bins, assume it's a continuous 
    // distribution and do NOT add dimension to biasCor.
    // If >= 2 logMass bins, add this dimension to biasCor.

    NBINg = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin;
    USEDIM_GAMMADM = (NBINg > 1 && NBINg < 5);
    USEDIM_LOGMASS = (NBINm > 1 ) ;

    if ( USEDIM_GAMMADM && !USEDIM_LOGMASS ) {
      NDIM_BIASCOR = 6; 
      sprintf(INFO_BIASCOR.STRING_PARLIST,"z,x1,c,a,b,g");  
    }

    if ( !USEDIM_GAMMADM && USEDIM_LOGMASS ) {
      NDIM_BIASCOR = 6; 
      sprintf(INFO_BIASCOR.STRING_PARLIST,"z,m,x1,c,a,b");  
    }

    if ( USEDIM_GAMMADM && USEDIM_LOGMASS ) {
      NDIM_BIASCOR = 7; 
      sprintf(INFO_BIASCOR.STRING_PARLIST,"z,m,x1,c,a,b,g");  
    }

    if ( DOCOR_MUCOV ) 
      { sprintf(txt_biasCor,"%dD+MUCOV", NDIM_BIASCOR); }
    else
      { sprintf(txt_biasCor,"%dD", NDIM_BIASCOR); }
  }
  else {
    NDIM_BIASCOR = 1; 
    sprintf(INFO_BIASCOR.STRING_PARLIST,"z");
    if ( DOCOR_1DZAVG )   
      { sprintf(txt_biasCor,"1D (WGT=1)") ; }
    else if ( DOCOR_1DZWGT ) 
      { sprintf(txt_biasCor,"1D (WGT=1/COV)"); }
    else if ( DOCOR_1D5DCUT ) 
      { sprintf(txt_biasCor,"1D (WGT=1/COV,5DCUT)"); }
    else {
      sprintf(c1err,"Invalid opt_biascor=%d", OPTMASK );
      sprintf(c2err,"grep MASK_BIASCOR  SALT2mu.c | grep define ");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  if ( DOCOR_1D && DOCOR_5D ) {
    sprintf(c1err,"Cannot do 1D and 5D biasCor.");
    sprintf(c2err,"opt_biascor logic is messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  printf("\n\t Use simulaton to interpolate %s biasCor\n", 
	 txt_biasCor);
  printf("\t\t (opt_biascor=%d, sigma_cell=%.2f)\n", 
	 OPTMASK, INPUTS.sigma_cell_biasCor);
  fflush(stdout);

  // -------------------------------------------
  int PS0 = INPUTS.prescale_biasCor[0];
  int PS1 = INPUTS.prescale_biasCor[1];
  if ( PS1 > 1 ) {
    printf("\t Apply pre-scale: select subset %d of %d \n", PS0, PS1);  
  }

  // Jan 16 2018: force correct option for log(sigma) term in chi2
  if ( NDIM_BIASCOR >=5  ) 
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
  

  // count number of biasCor events passing cuts (after setup_BININFO calls)

  for(ievt=0; ievt < INFO_BIASCOR.TABLEVAR.NSN_ALL; ievt++ ) 
    { compute_more_TABLEVAR(ievt, &INFO_BIASCOR.TABLEVAR ); }
  if ( NDIM_BIASCOR >=5 ) { store_iaib_biasCor(); }
  

  print_eventStats(EVENT_TYPE);

  //  if ( NDIM_BIASCOR == 1 ) { goto CHECK_1DCOR ; }
  if ( NDIM_BIASCOR == 1 && !DOCOR_1D5DCUT ) { goto CHECK_1DCOR ; }

  // make sparse list for faster looping below (Dec 21 2017)
  makeSparseList_biasCor();

  // determine sigInt for biasCor sample BEFORE makeMap since
  // sigInt is needed for 1/muerr^2 weight
  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ ) 
    {  init_sigInt_biasCor_SNRCUT(IDSAMPLE);  } 

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
    for(INDX = ILCPAR_MIN; INDX <= ILCPAR_MAX ; INDX++ ) 
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
  printf("   * mb,x1,c-bias at each alpha,beta,gammaDM \n");
  if ( DOCOR_MUCOV) 
    { printf("   * muCOVscale   at each alpha,beta,gammaDM \n"); }


  int DUMPFLAG = 0 ;
  for (n=0; n < NSN_DATA; ++n) {

    CUTMASK  = INFO_DATA.TABLEVAR.CUTMASK[n]; 
    IDSAMPLE = INFO_DATA.TABLEVAR.IDSAMPLE[n]; 
    if ( CUTMASK ) { continue ; }

    DUMPFLAG = 0 ; // (NUSE_TOT == 1 ) ; // xxx REMOVE
    istore = storeDataBias(n,DUMPFLAG);
    
    NUSE[IDSAMPLE]++ ; NUSE_TOT++ ;
    if ( istore == 0 )  { 
      NSKIP_TOT++; NSKIP[IDSAMPLE]++ ;  
      setbit_CUTMASK(n, CUTBIT_BIASCOR, &INFO_DATA.TABLEVAR); 
      if( INPUTS.dumpflag_nobiasCor && NSKIP_TOT<10 ) 
	{ storeDataBias(n,1); } 
      
    }    
  }


  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR; IDSAMPLE++ ) {

    // store NSKIP and NUSE for later to print warning about excess loss
    NDATA_BIASCORCUT[0][IDSAMPLE] = NUSE[IDSAMPLE];
    NDATA_BIASCORCUT[1][IDSAMPLE] = NSKIP[IDSAMPLE];

    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_SELECT == 0 ) { continue; }
    printf("  Rejected %d (of %d) DATA events with no %dD biasCor (%s).\n",
	   NSKIP[IDSAMPLE], NUSE[IDSAMPLE], 
	   NDIM_BIASCOR, SAMPLE_BIASCOR[IDSAMPLE].NAME );
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
  
  if ( NDIM_BIASCOR == 1 ) {
    INPUTS.opt_biasCor |= MASK_BIASCOR_1DZ ; // set OR-bit (AVG or WGT)
    prepare_biasCor_zinterp();  
  }

  INFO_BIASCOR.NDIM = NDIM_BIASCOR ;

  return ;

} // end prepare_biasCor


// =============================================
void  read_simFile_biasCor(void) {

  // Created June 2019
  // Read simFile and load simdata_bias struct.


  int NFILE      = INPUTS.nfile_biasCor;
  int NROW, ISTART, IFILETYPE, ifile, NVAR_ORIG, LEN_MALLOC ;   
  int NEVT[MXFILE_BIASCOR], NEVT_TOT;
  char *simFile ;
  //  char fnam[] = "read_simFile_biasCor" ;

  // --------------- BEGIN ---------------

  t_read_biasCor[0] = time(NULL); 

  // read each file quickly to get total size of arrays to malloc
  NEVT_TOT = 0 ;
  for(ifile=0; ifile < NFILE; ifile++ ) {
    simFile     = INPUTS.simFile_biasCor[ifile];
    NEVT[ifile] = SNTABLE_NEVT(simFile,TABLENAME_FITRES);  
    NEVT_TOT   += NEVT[ifile];
    printf("\t Found %d events in %s. \n", NEVT[ifile], simFile);
    fflush(stdout);
  }

  // malloc arrays for all sim files
  LEN_MALLOC = NEVT_TOT + 10 ;
  INFO_BIASCOR.TABLEVAR.LEN_MALLOC = LEN_MALLOC ;
  malloc_INFO_BIASCOR(+1,LEN_MALLOC);
  
  // loop again over each data file: read and append arrays
  for(ifile=0; ifile < NFILE; ifile++ ) {
    simFile     = INPUTS.simFile_biasCor[ifile];
    IFILETYPE   = TABLEFILE_OPEN(simFile,"read");
    NVAR_ORIG   = SNTABLE_READPREP(IFILETYPE,"FITRES");
    ISTART      = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
    SNTABLE_READPREP_TABLEVAR(ifile, ISTART, NEVT[ifile], 
			      &INFO_BIASCOR.TABLEVAR);
    NROW = SNTABLE_READ_EXEC(); // read entire file; load arrays
    INFO_BIASCOR.TABLEVAR.NSN_ALL += NROW ;

    INFO_BIASCOR.TABLEVAR.EVENT_RANGE[ifile][0] = ISTART ;
    INFO_BIASCOR.TABLEVAR.EVENT_RANGE[ifile][1] = ISTART + NROW - 1;
    sprintf(INFO_BIASCOR.TABLEVAR.INPUT_FILE[ifile],"%s", simFile);
    store_input_varnames(ifile, &INFO_BIASCOR.TABLEVAR) ;
  }


  t_read_biasCor[1] = time(NULL); 
  

  return ;
} // end read_simFile_biasCor



// ================================================================
void  prepare_biasCor_zinterp(void) {

  // Jun 28 2016 
  // Special biasCor option to correct mb(z) using muBias(sample,z).
  // Evaluate mubias(sample,z) from biasCor sample.
  //

  int    NSAMPLE  = NSAMPLE_BIASCOR ;
  double alpha    = INPUTS.parval[IPAR_ALPHA0] ;
  double beta     = INPUTS.parval[IPAR_BETA0] ;

  int    NSN_DATA, NSN_BIASCOR, ievt, NCUTS, iz, izusr, idsample, cutmask ;
  double mu_fit, mB_fit, x1_fit, c_fit ;
  double mu_sim, mB_sim, x1_sim, c_sim ;
  double mu_err, mB_err, x1_err, c_err ;
  double mu_true, WGT, z, zMIN[MXNUM_SAMPLE], zMAX[MXNUM_SAMPLE] ;
  double zlo[MXNUM_SAMPLE][MXz], zhi[MXNUM_SAMPLE][MXz];
  double zavg[MXNUM_SAMPLE][MAXBIN_z], zbinSize[MXNUM_SAMPLE], zbin;
  double USR_SUMz[MXz], USR_SUMWGT[MXz];
  double MUBIAS_FIT[MXNUM_SAMPLE][MXz];
  double MUBIAS_SIM[MXNUM_SAMPLE][MXz];
  double SUMWGT[MXNUM_SAMPLE][MXz];
  double SUMWGT_HISNR, MUBIAS_SIM_HISNR ;
  double muerrsq, muerrLens, SNRMAX, *ptr_zM0 ;

  int ndata[MXNUM_SAMPLE], ndata_reject[MXNUM_SAMPLE];
  int NDATA_REJECT=0;
  char fnam[] = "prepare_biasCor_zinterp" ;
  
  // ---------------- BEGIN --------------

  INFO_BIASCOR.SIGINT_AVG = .11 ;  
  NSN_DATA     = INFO_DATA.TABLEVAR.NSN_ALL ;
  NSN_BIASCOR  = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptr_zM0      = INFO_BIASCOR.zM0;


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
	   fflush(stdout); 
    */
  }
  
  for(iz=0; iz < MXz; iz++ ) {

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

  for(ievt=0; ievt < NSN_BIASCOR; ievt++ ) { 

    z         = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt];
    idsample  = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt];
    cutmask   = (int)INFO_BIASCOR.TABLEVAR.CUTMASK[ievt] ;

    if ( cutmask )            { continue; }
    if ( z < zMIN[idsample] ) { continue ; }
    if ( z > zMAX[idsample] ) { continue ; }

    NCUTS++ ; 

    mB_fit    = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_mB][ievt] ;
    x1_fit    = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_x1][ievt] ;
    c_fit     = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c ][ievt] ;
    mB_err    = (double)INFO_BIASCOR.TABLEVAR.fitpar_err[INDEX_mB][ievt] ;
    x1_err    = (double)INFO_BIASCOR.TABLEVAR.fitpar_err[INDEX_x1][ievt] ;
    c_err     = (double)INFO_BIASCOR.TABLEVAR.fitpar_err[INDEX_c ][ievt] ;
    mB_sim    = (double)INFO_BIASCOR.TABLEVAR.SIM_FITPAR[INDEX_mB][ievt] ;
    x1_sim    = (double)INFO_BIASCOR.TABLEVAR.SIM_FITPAR[INDEX_x1][ievt] ;
    c_sim     = (double)INFO_BIASCOR.TABLEVAR.SIM_FITPAR[INDEX_c ][ievt] ;
    SNRMAX    = (double)INFO_BIASCOR.TABLEVAR.snrmax[ievt];
    mu_true   = (double)INFO_BIASCOR.TABLEVAR.SIM_MU[ievt]; 

    mu_fit  = mB_fit + alpha*x1_fit - beta*c_fit - M0_DEFAULT ;
    mu_sim  = mB_sim + alpha*x1_sim - beta*c_sim - M0_DEFAULT ;
    mu_err  = -9.0; // not used

    muerrLens = ( z * INPUTS.lensing_zpar ) ;

    iz = (int)( (z-zMIN[idsample]) / zbinSize[idsample] ) ;
    
    if ( iz >= MXz ) {
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

    izusr = IBINFUN(z, &INPUTS.BININFO_z, 0, fnam);
    USR_SUMz[izusr]   += (WGT*z) ;
    USR_SUMWGT[izusr] +=  WGT ;

  } // end loop over biasCor events

  printf("\t %d of %d biasCor events pass cuts. \n",  NCUTS, NSN_BIASCOR);
  INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS = NCUTS; 

  if ( NCUTS < 10 ) {
    print_eventStats(EVENT_TYPE_BIASCOR);
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
    for(iz=0; iz < MXz; iz++ ) {

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
      { ptr_zM0[izusr] = USR_SUMz[izusr]/USR_SUMWGT[izusr];}
    else
      { ptr_zM0[izusr] = -9.0 ; }
  }

  // -----------------------------------------------
  printf("\n   Now compute muBias(z) for each data event. \n");
  fflush(stdout);

  double mB, muBias ;

  for(ievt=0; ievt < NSN_DATA; ievt++ ) {

    INFO_DATA.muBias_zinterp[ievt]  = 0.0;

    idsample = INFO_DATA.TABLEVAR.IDSAMPLE[ievt];
    cutmask = INFO_DATA.TABLEVAR.CUTMASK[ievt];
    if ( cutmask ) { continue ; }

    if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) 
      { setbit_CUTMASK(ievt, CUTBIT_BIASCOR, &INFO_DATA.TABLEVAR); }

    z        = INFO_DATA.TABLEVAR.zhd[ievt] ;
    mB       = INFO_DATA.TABLEVAR.fitpar[INDEX_mB][ievt];
    izmin    = IZMIN[idsample];
    izmax    = IZMAX[idsample];
    
    if ( z < zlo[idsample][izmin] || z > zhi[idsample][izmax] ) {
      setbit_CUTMASK(ievt, CUTBIT_BIASCOR, &INFO_DATA.TABLEVAR); 
      ndata_reject[idsample]++ ;
      NDATA_REJECT++ ;
    }

    cutmask = INFO_DATA.TABLEVAR.CUTMASK[ievt];
    if ( cutmask ) { continue ; }

    ndata[idsample]++ ;

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


    INFO_DATA.muBias_zinterp[ievt] = muBias ;

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

  
  // set  CELLINFO_BIASCOR.MAPCELL[ia][ib][ig][iz][im][ix1][ic]
  // to map 5D indices to 1D index J1D.
  // July 31 2019: add logmass dimension [im]

  int ID = IDSAMPLE;
  int ia, ib, ig, iz, im, ix1, ic, NCELL, ipar, NDIM ;
  char fnam[] = "set_MAPCELL_biasCor" ;

  // -------------- BEGIN ----------------------------
  int NBINz   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINm   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_m.nbin ;
  int NBINx1  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc   = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;
  int NBINa   = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int NBINb   = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int NBINg   = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;

  // ------------------------------------------
  // establish map between 1D array and 7D binning (only do once)
  
  NCELL = 0 ;
  
  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(ig=0; ig<NBINg; ig++ ) {
	for(iz=0; iz<NBINz; iz++ ) {
	  for(im=0; im<NBINm; im++ ) {
	    for(ix1=0; ix1<NBINx1; ix1++ ) {
	      for(ic=0; ic<NBINc; ic++ ) {
		if ( NCELL < MAXBIN_BIASCOR_1D )  {
		  CELLINFO_BIASCOR[ID].MAPCELL[ia][ib][ig][iz][im][ix1][ic] 
		    = NCELL ; 
		}
		NCELL++ ;
	      }
	    }
	  }	 
	}
      }
    }
  }

 

  // check array bound.
  if ( NCELL >= MAXBIN_BIASCOR_1D ) {
    sprintf(c1err,"NCELL=%d > bound of MAXBIN_BIASCOR_1D=%d", 
	    NCELL, MAXBIN_BIASCOR_1D ) ;
    sprintf(c2err,"NBIN(a,b,g,  z,m,x1,c) = %d,%d,%d   %d,%d,%d,%d",
	    NBINa, NBINb, NBINg, NBINz, NBINm, NBINx1, NBINc);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  CELLINFO_BIASCOR[IDSAMPLE].NCELL = NCELL ;

  NDIM = 5;
  if ( NBINg > 1 ) { NDIM++ ; }
  if ( NBINm > 1 ) { NDIM++ ; }
  // xxx mark del  if  ( NBINg <= 1 ) { NDIM=5; } else { NDIM=7; }

  printf("  %s : malloc for %d  %dD-cells for biasCor \n", fnam, NCELL, NDIM );
  fflush(stdout);

  // malloc other CELLINFO arrays
  int MEMD0   = NCELL   * sizeof(double);
  int MEMI0   = NCELL   * sizeof(int);

  CELLINFO_BIASCOR[IDSAMPLE].NperCell = (int   *) malloc(MEMI0);
  CELLINFO_BIASCOR[IDSAMPLE].AVG_z    = (double*) malloc(MEMD0);
  CELLINFO_BIASCOR[IDSAMPLE].AVG_m    = (double*) malloc(MEMD0);
  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    {  CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar] = (double*)malloc(MEMD0);}

  // -------------------------------------
  // malloc FITPARBIAS & MUCOVSCALE based on number of 7D cells
  // and number of SURVEY-FIELDs.

  int MEMBIAS, MEMCOV;
  MEMBIAS = NCELL   * sizeof(FITPARBIAS_DEF) ;
  MEMCOV  = NCELL   * sizeof(float) ;
  INFO_BIASCOR.FITPARBIAS[IDSAMPLE] = (FITPARBIAS_DEF*) malloc ( MEMBIAS);
  INFO_BIASCOR.MUCOVSCALE[IDSAMPLE] = (float *        ) malloc ( MEMCOV );

  fflush(stdout);
  return ;

} // end set_MAPCELL_biasCor


// ================================================================
void makeMap_fitPar_biasCor(int IDSAMPLE, int ipar_LCFIT) {

  // Created Jan 2016 by R.Kessler
  // After reading simFile_bias, prepare bias vs. redshift & parval
  // so that fitting model can be corrected. Fill struct simdata_bias.
  // ipar_LCFIT is for mB, x1 or c.
  // Note that mB includes gammadm correction from host.
  //
  // WARNING: computes bias with user z-bining, so beware of
  //          combined low-z plus hi-z samples.
  //
  // Aug 26 2019: account for gammadm
  // Feb 24 2020: update for ipar_LCFIT = index_mu
  //
  // - - - - - - - - - -

#define MXJ1DNBR 10000

  int NCELL          = CELLINFO_BIASCOR[IDSAMPLE].NCELL ;
  int NBIASCOR_CUTS  = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;
  char *PARNAME      = BIASCOR_NAME_LCFIT[ipar_LCFIT] ;

  int MEMD, NEVT_USE, NEVT_SAMPLE,  LDMP, N ;
  int J1D, ia, ib, ig, iz, im, ix1, ic, ievt, isp ;
  int J1DNBR_LIST[MXJ1DNBR], NJ1DNBR ;

  double fit_val, sim_val, biasVal, sim_gammadm ;
  double WGT, VAL, ERR, RMS, SQRMS, XN, XNLIST, tmp1, tmp2 ;
  double *SUMBIAS, *SUMWGT, *sum, *sumsq ;

  float *ptr_fitpar, *ptr_simpar, *ptr_gammadm;
  char fnam[] = "makeMap_fitPar_biasCor";

  // ----------------- BEGIN -----------------

  printf("  %s of %s-bias(%s)  \n", 
	 fnam, PARNAME, INFO_BIASCOR.STRING_PARLIST);
  fflush(stdout);

  // malloc arrays to store info in each biasCor cell
  MEMD    = NCELL * sizeof(double) ;
  SUMBIAS = (double*) malloc(MEMD) ;
  SUMWGT  = (double*) malloc(MEMD) ;
  sum     = (double*) malloc(MEMD) ;
  sumsq   = (double*) malloc(MEMD) ;

  ptr_fitpar  = INFO_BIASCOR.TABLEVAR.fitpar[ipar_LCFIT] ;
  ptr_simpar  = INFO_BIASCOR.TABLEVAR.SIM_FITPAR[ipar_LCFIT] ;
  ptr_gammadm = INFO_BIASCOR.TABLEVAR.SIM_GAMMADM ;

  // ------------------------------------------
  for( J1D=0; J1D < NCELL ; J1D++ ) {
    // init global arrays  
    
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT] = 999. ;
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT] = 999. ;
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].RMS[ipar_LCFIT] = 999. ;    
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
    fit_val  = (double)ptr_fitpar[ievt];
    sim_val  = (double)ptr_simpar[ievt];

    // Aug 26 2019: apply gammadm correction to simval
    if ( ipar_LCFIT == INDEX_mB ) {
      sim_gammadm = ptr_gammadm[ievt];  // from biasCor table
      sim_val    += sim_gammadm ;       // true mB is SIM_mB + true gammadm
    }

    biasVal = fit_val - sim_val ; 

    WGT = WGT_biasCor(2,ievt,fnam) ;  // WGT=1/muerr^2 for wgted average
    J1D = J1D_biasCor(ievt,fnam);     // 1D index

    SUMBIAS[J1D]  += (WGT * biasVal) ;
    SUMWGT[J1D]   += WGT ;

    // increment unweighted sums to get RMS and bias-error
    sum[J1D]   += (biasVal) ;
    sumsq[J1D] += (biasVal * biasVal) ;  

    // -----------------
    NEVT_USE++ ;
    CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D]++ ;
  }

  // convert sums in each IZ,LCFIT bin into mean bias,
  
  int    NMAP_TOT = 0 ;
  int    NMAP_USE = 0 ;
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
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT] = VAL ;
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT] = ERR ;
    INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].RMS[ipar_LCFIT] = RMS ;    

  } // end J1D loop


  // -----------------------------------------------
  // print grid-cell stats on last parameter
  // (since it's the same for each parameter)
  if ( ipar_LCFIT == INFO_BIASCOR.ILCPAR_MAX ) {
    printf("  BiasCor computed for %d of %d grid-cells with >=1 events.\n",
	   NMAP_USE, NMAP_TOT ) ;
    printf("  BiasCor sample: %d of %d pass cuts for IDSAMPLE=%d.\n",
	   NEVT_USE, NEVT_SAMPLE, IDSAMPLE );

    if ( NEVT_USE == 0 ) {
      print_eventStats(EVENT_TYPE_BIASCOR);
      sprintf(c1err,"No BiasCor events passed for %s", 
	      SAMPLE_BIASCOR[IDSAMPLE].NAME );
      sprintf(c2err,"Check BiasCor file" );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }

    fflush(stdout);
  }
  

  // -----------------
  // debug dump
  LDMP = 0 ; 
  if ( LDMP ) {
    iz=7; im=0; ix1=5; ic=6; ia=0; ib=0;
    J1D = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][ig][iz][im][ix1][ic] ;
    VAL = INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].VAL[ipar_LCFIT];
    ERR = INFO_BIASCOR.FITPARBIAS[IDSAMPLE][J1D].ERR[ipar_LCFIT];


    printf(" xxx --------------------------------------- \n");
    printf(" xxx %s-bias = %.3f +- %.3f for \n"
	   " xxx \t z[%.3f:%.3f] x1[%.3f:%.3f] c[%.3f:%.3f]  N=%d\n"
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[ipar_LCFIT].varName
	   ,VAL, ERR
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.lo[iz]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.hi[iz]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].lo[ix1]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].hi[ix1]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].lo[ic]
	   ,CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].hi[ic]
	   ,CELLINFO_BIASCOR[IDSAMPLE].NperCell[J1D] 	   );
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

  int ID     = IDSAMPLE;
  int NBINz  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin ;
  int NBINm  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_m.nbin ;
  int NBINx1 = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc  = CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;

  int iz, im, ix1, ic ;
  int IA, IB, IG, IZ, IM, IX1, IC ;
  int NJ1D_local = 0 ;
  int J1D_local ;

  //  char fnam[] = "get_J1DNBR_LIST";

  // ----------------- BEGIN --------------

  // convert J1D into index for each dimension
  J1D_invert_I(IDSAMPLE, J1D, &IA, &IB, &IG, &IZ, &IM, &IX1, &IC );

  /*
  printf(" xxx J1D=%d -> IA,IB,IZ,IX1,IC = %d %d %d %d %d \n",
	 J1D, IA,IB,IZ,IX1,IC ); fflush(stdout); 
  */

 
  for(iz=IZ-1; iz <= IZ+1; iz++ ) {
    if ( iz <  0     ) { continue ; }
    if ( iz >= NBINz ) { continue ; }

    for(im=IM-1; im <= IM+1; im++ ) {
      if ( im <  0     ) { continue ; }
      if ( im >= NBINm ) { continue ; }
      
      for(ix1=IX1-1; ix1<=IX1+1 ; ix1++ ) {
	if ( ix1 <  0      ) { continue ; }
	if ( ix1 >= NBINx1 ) { continue ; }
	
	for(ic=IC-1; ic<=IC+1; ic++ ) {
	  if ( ic <  0     ) { continue ; }
	  if ( ic >= NBINc ) { continue ; }
	  J1D_local=CELLINFO_BIASCOR[ID].MAPCELL[IA][IB][IG][iz][im][ix1][ic];
	  J1DNBR_LIST[NJ1D_local] = J1D_local ;
	  NJ1D_local++ ;
	}
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
  // Aug 22 2019: include logmass dependence

  int  istat_cov, J1D, idsample ;
  double WGT, muerrsq ;
  char fnam[] = "WGT_biasCor" ;
  
  // --------------- BEGIN -------------------

  muerrsq  = muerrsq_biasCor(ievt, USEMASK_BIASCOR_COVTOT, &istat_cov) ; 
  WGT      = 1.0/muerrsq ;

  if ( opt == 1 ) { return(WGT); }

  if ( INPUTS.sigma_cell_biasCor > 10.0 ) { return(WGT); } // default

  // If we get here, compute additional weight based on 3D/4D
  // separation from wgted-avg in cell.
  double z, m, x1, c, Dz, Dm, Dx1, Dc, SQD, ARG, WGT_CELL, SQSIGMA_CELL ;
  double binSize_z, binSize_m, binSize_x1, binSize_c ;

  idsample = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt];
  z   = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt] ;
  x1  = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_x1][ievt] ;
  c   = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][ievt] ;   
  

  J1D          = J1D_biasCor(ievt,fnam);   // 1D index
  binSize_z    = CELLINFO_BIASCOR[idsample].BININFO_z.binSize ;
  binSize_x1   = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1].binSize ;
  binSize_c    = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c].binSize ;

  
  Dz  = (z  - CELLINFO_BIASCOR[idsample].AVG_z[J1D])/binSize_z ;
  Dx1 = (x1 - CELLINFO_BIASCOR[idsample].AVG_LCFIT[INDEX_x1][J1D])/binSize_x1 ;
  Dc  = (c  - CELLINFO_BIASCOR[idsample].AVG_LCFIT[INDEX_c ][J1D])/binSize_c ;
  SQD = Dz*Dz + Dx1*Dx1 + Dc*Dc ;

  // check for logmass contribution
  int NBINm  = CELLINFO_BIASCOR[idsample].BININFO_m.nbin ;
  if ( NBINm > 1 ) {
    m          = (double)INFO_BIASCOR.TABLEVAR.logmass[ievt] ;
    binSize_m  = CELLINFO_BIASCOR[idsample].BININFO_m.binSize ;
    Dm         = (m  - CELLINFO_BIASCOR[idsample].AVG_m[J1D])/binSize_m ;
    SQD       += (Dm*Dm);
  }


  SQSIGMA_CELL = INPUTS.sigma_cell_biasCor * INPUTS.sigma_cell_biasCor ;
  ARG      = -0.5 * (SQD/SQSIGMA_CELL);
  WGT_CELL = exp(ARG) ;
  WGT     *= WGT_CELL ;

  return(WGT);

} // end WGT_biasCor

// ===========================
int J1D_biasCor(int ievt, char *msg ) {

  // Created May 12 2016
  // For biasCor event "ievt", return 1D index to biasCor map.
  // *msg is used only for error message.

  double a, b, g, z, m, x1, c ;
  char *name, MSGERR[100] ;
  int ia, ib, ig, iz, im, ix1, ic, J1D, idsample ;
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM;
  BININFO_DEF *BININFO_z, *BININFO_m, *BININFO_x1, *BININFO_c;

  char fnam[] = "J1D_biasCor" ;

  // -------------- BEGIN ------------

  J1D = -9 ;

  a        = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt] ;
  b        = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt] ;
  g        = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt] ;
  z        = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt];
  m        = (double)INFO_BIASCOR.TABLEVAR.logmass[ievt];
  x1       = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_x1][ievt] ;
  c        = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][ievt] ;
  name     = INFO_BIASCOR.TABLEVAR.name[ievt];
  idsample = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt];
  BININFO_SIM_ALPHA    = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA     = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM  = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;
  
  BININFO_z  = &CELLINFO_BIASCOR[idsample].BININFO_z ;
  BININFO_m  = &CELLINFO_BIASCOR[idsample].BININFO_m ;
  BININFO_x1 = &CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1];
  BININFO_c  = &CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c];

  // find integer bin for each map variable to store bias
  sprintf(MSGERR, "%s : CID=%s  z=%.3f m=%.2f x1=%.3f c=%.3f", 
	  fnam, name, z, m, x1, c);

  ia  = IBINFUN(a,  BININFO_SIM_ALPHA,    1, MSGERR ) ;
  ib  = IBINFUN(b,  BININFO_SIM_BETA,     1, MSGERR ) ;
  ig  = IBINFUN(g,  BININFO_SIM_GAMMADM,  1, MSGERR ) ;
  iz  = IBINFUN(z,  BININFO_z,            1, MSGERR ) ;
  im  = IBINFUN(m,  BININFO_m,            1, MSGERR ) ;
  ix1 = IBINFUN(x1, BININFO_x1,           1, MSGERR ) ;
  ic  = IBINFUN(c,  BININFO_c,            1, MSGERR ) ;

  int NBINz  = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
  int NBINm  = CELLINFO_BIASCOR[idsample].BININFO_m.nbin ;
  int NBINx1 = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_x1].nbin ;
  int NBINc  = CELLINFO_BIASCOR[idsample].BININFO_LCFIT[INDEX_c].nbin ;
  if ( iz  < 0 || iz  >= NBINz  ) { return J1D; }
  if ( im  < 0 || im  >= NBINm  ) { return J1D; }
  if ( ix1 < 0 || ix1 >= NBINx1 ) { return J1D; }
  if ( ic  < 0 || ic  >= NBINc  ) { return J1D; }
  
  J1D = CELLINFO_BIASCOR[idsample].MAPCELL[ia][ib][ig][iz][im][ix1][ic] ;

  if ( J1D >= MAXBIN_BIASCOR_1D || J1D < 0 ) 
    {  J1D = -9 ;  }

  return(J1D);

} // end J1D_biasCor


// ======================================================
void  J1D_invert_D(int IDSAMPLE, int J1D, 
		   double *a, double *b, double *g,
		   double *z, double *m, double *x1, double *c) {

  // May 2016: for input J1D index, return doubles a,b,z,x1,c

  int ID = IDSAMPLE;
  int NBINa, NBINb, NBINg, NBINz, NBINm, NBINx1, NBINc;
  int ia, ib, ig, iz, im, ix1, ic, J1TMP ;

  BININFO_DEF *BININFO_z, *BININFO_m, *BININFO_x1, *BININFO_c;
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM;

  // --------------- BEGIN ----------------

  *a = *b = *g = *z = *m = *x1 = *c = -999.9 ;

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;

  BININFO_z         = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z;
  BININFO_m         = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_m;
  BININFO_x1        = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1];
  BININFO_c	    = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c];

  NBINa   = (*BININFO_SIM_ALPHA).nbin ;
  NBINb   = (*BININFO_SIM_BETA).nbin ;
  NBINg   = (*BININFO_SIM_GAMMADM).nbin ;
  NBINz   = (*BININFO_z).nbin ;
  NBINm   = (*BININFO_m).nbin ;
  NBINx1  = (*BININFO_x1).nbin ;
  NBINc   = (*BININFO_c).nbin ;
	  
 
  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(ig=0; ig<NBINg; ig++ ) {
	for(iz=0; iz<NBINz; iz++ ) {
	  for(im=0; im<NBINm; im++ ) {
	    for(ix1=0; ix1<NBINx1; ix1++ ) {
	      for(ic=0; ic<NBINc; ic++ ) {
		
		J1TMP=CELLINFO_BIASCOR[ID].MAPCELL[ia][ib][ig][iz][im][ix1][ic];
		if ( J1TMP == J1D ) {
		  *a  = (*BININFO_SIM_ALPHA).avg[ia];
		  *b  = (*BININFO_SIM_BETA).avg[ib];
		  *g  = (*BININFO_SIM_GAMMADM).avg[ib];
		  *z  = (*BININFO_z).avg[iz];
		  *m  = (*BININFO_m).avg[im];
		  *x1 = (*BININFO_x1).avg[ix1];
		  *c  = (*BININFO_c).avg[ic];
		  return ;
		}	      
	      }
	    }
	  }
	}
      }
    }
  }
  
} // end J1D_invert_D

// ======================================================
void  J1D_invert_I(int IDSAMPLE, int J1D, int *ja, int *jb, int *jg,
		   int *jz, int *jm, int *jx1, int *jc) {

  // May 2016: for input J1D index, return doubles a,b,z,x1,c

  int ID = IDSAMPLE;
  int NBINa, NBINb, NBINg, NBINz, NBINm, NBINx1, NBINc;
  int ia, ib, ig, iz, im, ix1, ic, J1TMP ;
  BININFO_DEF *BININFO_z, *BININFO_m, *BININFO_x1, *BININFO_c;
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM;

  // --------------- BEGIN ----------------

  *ja = *jb = *jg = *jz = *jm = *jx1 = *jc = -9 ;

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;

  BININFO_z         = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z;
  BININFO_m         = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_m;
  BININFO_x1        = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_x1];
  BININFO_c	    = &CELLINFO_BIASCOR[IDSAMPLE].BININFO_LCFIT[INDEX_c];

  NBINa   = (*BININFO_SIM_ALPHA).nbin ;
  NBINb   = (*BININFO_SIM_BETA).nbin ;
  NBINg   = (*BININFO_SIM_GAMMADM).nbin ;
  NBINz   = (*BININFO_z).nbin ;
  NBINm   = (*BININFO_m).nbin ;
  NBINx1  = (*BININFO_x1).nbin ;
  NBINc   = (*BININFO_c).nbin ;

  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(ig=0; ig<NBINg; ig++ ) {
	for(iz=0; iz<NBINz; iz++ ) {
	  for(im=0; im<NBINm; im++ ) {
	    for(ix1=0; ix1<NBINx1; ix1++ ) {
	      for(ic=0; ic<NBINc; ic++ ) {
		
		J1TMP=CELLINFO_BIASCOR[ID].MAPCELL[ia][ib][ig][iz][im][ix1][ic];
		if ( J1TMP == J1D ) {
		  *ja  = ia ;
		  *jb  = ib ;
		  *jg  = ig ;
		  *jz  = iz ;
		  *jm  = im ;
		  *jx1 = ix1 ;
		  *jc  = ic ;
		  return ;
		}	      	  
	      }    
	    }
	  }
	}
      }
    }
  }


} // end J1D_invert_I


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
  // Jan 17 2020: fix gamma-dimension that tripped valgrind errors.
  //

  //  int ID = IDSAMPLE;
  int NBIASCOR_CUTS = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;
  int    NBINa, NBINb, NBINg, NBINz, NBINc ;
  double cmin,  cmax, cbin, c_lo, c_hi;

  int DUMPFLAG = 0 ;
  int    ia, ib, ig, iz, ic, i1d, NCELL, isp ; 
  int    ievt, istat_cov, istat_bias, N, J1D, ipar ;
  double muErr, muErrsq, muDif, muDifsq, pull, tmp1, tmp2  ;
  double muBias, muBiasErr, muCOVscale, fitParBias[NLCPAR+1] ;
  double a, b, gDM, z, m, c ;
  double *SUM_MUERR, *SUM_SQMUERR;
  double *SUM_MUDIF, *SUM_SQMUDIF ;
  double *SQMUERR,   *SQMURMS ;
  double *SUM_PULL,  *SUM_SQPULL ;
  char   *name ;

  float    *ptr_MUCOVSCALE;
  BIASCORLIST_DEF     BIASCORLIST ;
  FITPARBIAS_DEF      FITPARBIAS[MXa][MXb][MXg] ;
  double              MUCOVSCALE[MXa][MXb][MXg] ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
 
  char fnam[]  = "makeMap_sigmu_biasCor" ;
  
  // ----------------- BEGIN -------------------

  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return; }

  printf(" %s: make map of muCOVscale(a,b,g,z,c) \n", fnam);
  fflush(stdout);

  NBINa    = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  NBINb    = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;  
  NBINg    = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;  
  ptr_MUCOVSCALE = INFO_BIASCOR.MUCOVSCALE[IDSAMPLE];   

  // redshift bins are same as for biasCor
  copy_BININFO(&CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 
	       &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z);
  NBINz = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.nbin ;


  // setup color bins that are coarser than those for muBias.
  // Note that other bins (a,b,z) are same for muCOVscale and muBias.
  // Beware hard-wired values here.
  NBINc=3; cmin=-0.3; cmax=+0.3; cbin=0.2; 
  // NBINc=1; cmin=-0.3; cmax=+0.3; cbin=0.6;  // return to < Jun 30 2016

  // xxx mark delete  NCELL    = NBINa * NBINb * NBINz * NBINc ;
  NCELL    = NBINa * NBINb * NBINg * NBINz * NBINc ;
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
  ig = 0 ;

  int N1D=0;
  for(ia=0; ia< NBINa; ia++ ) {
    for(ib=0; ib< NBINb; ib++ ) {  
      for(ig=0; ig< NBINg; ig++ ) {  
	MUCOVSCALE[ia][ib][ig] = 1.0 ; // dummy arg for get_muBias below
	for(iz=0; iz < NBINz; iz++ ) {
	  for(ic=0; ic < NBINc; ic++ ) {
	    SUM_MUERR[N1D] = SUM_SQMUERR[N1D] = 0.0 ;
	    SUM_MUDIF[N1D] = SUM_SQMUDIF[N1D] = 0.0 ;	
	    SUM_PULL[N1D]  = SUM_SQPULL[N1D]  = 0.0 ;
	    ptr_MUCOVSCALE[N1D] = 1.0 ;
	    CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[N1D]  = 0 ;
	    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_z[N1D]     = 0.0 ;
	    CELLINFO_MUCOVSCALE[IDSAMPLE].AVG_LCFIT[INDEX_c][N1D] = 0.0 ;
	    CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][ig][iz][0][0][ic]=N1D;
	    N1D++ ;
	  }	  
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

    for(ia=0; ia<MXa; ia++ ) {
      for(ib=0; ib<MXb; ib++ ) {
	for(ig=0; ig<MXg; ig++ ) {
	  zero_FITPARBIAS(&FITPARBIAS[ia][ib][ig] ); 
	}
      }
    }
    
    z    = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt];
    m    = (double)INFO_BIASCOR.TABLEVAR.logmass[ievt];
    a    = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt];
    b    = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt];
    gDM  = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt];
    c    = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][ievt];
    ia   = (int)INFO_BIASCOR.IA[ievt];
    ib   = (int)INFO_BIASCOR.IB[ievt];
    ig   = (int)INFO_BIASCOR.IG[ievt];
    name = INFO_BIASCOR.TABLEVAR.name[ievt];
    for(ipar=0; ipar < NLCPAR; ipar++ ) 
      { BIASCORLIST.FITPAR[ipar] = 
	  (double)INFO_BIASCOR.TABLEVAR.fitpar[ipar][ievt]; 
      }
    BIASCORLIST.FITPAR[INDEX_mu] = 0.0; // mu slot not used

    // allow color (c) to be outside map
    iz = IBINFUN(z, &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z, 
		 1, fnam );

    ic = IBINFUN(c, &CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c], 
		 2, fnam );

    // ---------------------------------------------------
    // need bias corrected distance to compute pull

    BIASCORLIST.z           = z ;
    BIASCORLIST.logmass     = m ;
    BIASCORLIST.alpha       = a ;
    BIASCORLIST.beta        = b ;
    BIASCORLIST.gammadm     = gDM ;
    BIASCORLIST.idsample    = IDSAMPLE ;
    
    //    DUMPFLAG = ( isp < 5 ); // xxxx REMOVE
    istat_bias = 
      get_fitParBias(name, &BIASCORLIST, DUMPFLAG,
		     &FITPARBIAS[ia][ib][ig] ); // <== returned

    //  DUMPFLAG = 0 ; // xxx REMOVE

    // skip if bias cannot be computed, just like for data
    if ( istat_bias == 0 ) { continue ; }

    get_INTERPWGT_abg(a,b,gDM, DUMPFLAG, &INTERPWGT, fnam );
    get_muBias(name, &BIASCORLIST, FITPARBIAS, MUCOVSCALE, &INTERPWGT,
	       fitParBias, &muBias, &muBiasErr, &muCOVscale );  

    // ----------------------------
    muDif   =  muresid_biasCor(ievt);  // mu - muTrue
    muDif  -=  muBias ;  
    muDifsq =  muDif*muDif ;

    // compute error with intrinsic scatter
    muErrsq = muerrsq_biasCor(ievt, USEMASK_BIASCOR_COVTOT, &istat_cov) ; 

    if ( muErrsq <= 1.0E-14 || muErrsq > 100.0 || isnan(muErrsq) ) {
      print_preAbort_banner(fnam);
      printf("\t z=%f  a=%f  b=%f  gDM=%f\n",
	     z, a, b, gDM);
      printf("\t ia,ib.ig = %d, %d, %d \n", ia, ib, ig);

      sprintf(c1err,"Invalid muErrsq=%f for ievt=%d", muErrsq, ievt);
      sprintf(c2err,"Something is messed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }

    muErr   = sqrt(muErrsq) ;    
    pull    = (muDif/muErr) ;

    // get 1d index
    i1d = CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][ig][iz][0][0][ic] ;

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
    ptr_MUCOVSCALE[i1d] = (float)SQRMS ;

  }  // i1d
  

  // -------------------------
  // print errBias info in z bins

  int LPRINT = 0 ;

  if ( LPRINT ) {

    double zlo, zhi ;  
    printf("\n");
    printf("                            "
	   "RMS(muDif)/RMS(Pull)/NSIM for \n");
    
    printf("  ia,ib,ig  z-range :   "
	   "    ic=0               ic=1                ic=2 \n");
    
    printf("  -------------------------------------------------"
	   "------------------------\n");
    fflush(stdout);
    
    for(ia=0; ia< NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {      
	for(ig=0; ig < NBINg; ig++ ) {      
	  for(iz=0; iz < NBINz; iz++ ) {	
	    zlo = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.lo[iz];
	    zhi = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.hi[iz];
	    printf("  %d,%d,%d  %.2f-%.2f : ",  ia,ib,ig, zlo, zhi );
	    
	    for(ic=0; ic<NBINc; ic++ ) {  
	      i1d=CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][ig][iz][0][0][ic];
	      N     = CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[i1d] ;
	      muCOVscale = (double)ptr_MUCOVSCALE[i1d] ;
	      RMS        = sqrt ( SQMURMS[i1d] );       
	      printf("%6.3f/%5.3f/%5d  ", RMS, sqrt(muCOVscale), N );	  
	    } // ic
	    
	    printf("\n");     fflush(stdout);
	  } // iz
	  printf("\n");     fflush(stdout);
	} // ig
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

  // store alpha, beta, gamma, zcmb (ia,ib,ig iz) for each biasCor event.
  // 
  // Jul 16 2019: include gammaDM index (ig).

  int NSN = INFO_BIASCOR.TABLEVAR.NSN_ALL; 
  double a, b, g, z, Z ;
  int i, IA, IB, IG, IZ, iz, ID;
  char fnam[] = "store_iaib_biasCor" ;

  // ------------ BEGIN ----------

  for(i=0; i < NSN; i++ ) {

    if ( INFO_BIASCOR.TABLEVAR.SIM_NONIA_INDEX[i] != 0 ) { continue ; }

    ID  = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[i];
    if ( ID < 0 ) { continue ; }

    a      = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[i] ;
    b      = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[i] ;  
    g      = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[i] ;  
    Z      = (double)INFO_BIASCOR.TABLEVAR.SIM_ZCMB[i] ;  // true ZCMB
    z      = (double)INFO_BIASCOR.TABLEVAR.zhd[i] ;       // measured zcmb

    IA  = IBINFUN(a, &INFO_BIASCOR.BININFO_SIM_ALPHA,   1, fnam ) ;
    IB  = IBINFUN(b, &INFO_BIASCOR.BININFO_SIM_BETA,    1, fnam ) ;
    IG  = IBINFUN(g, &INFO_BIASCOR.BININFO_SIM_GAMMADM, 1, fnam ) ;
    IZ  = IBINFUN(Z, &CELLINFO_BIASCOR[ID].BININFO_z, 0, fnam ) ;
    iz  = IBINFUN(z, &CELLINFO_BIASCOR[ID].BININFO_z, 0, fnam ) ;

    INFO_BIASCOR.IA[i] = (int8_t)IA ;
    INFO_BIASCOR.IB[i] = (int8_t)IB ;
    INFO_BIASCOR.IG[i] = (int8_t)IG ;
    INFO_BIASCOR.IZ[i] = (int8_t)IZ ;
    INFO_BIASCOR.iz[i] = (int8_t)iz ;     

  } // end loop over biasCor rows

  return ; 

} // end store_iaib_biasCor



// ==================================
double muresid_biasCor(int ievt ) {

  // Return muFit - muTrue for biasCor sample event "ievt"
  // Note that there is no bias correction here, just the
  // naive distance from Trip formula.
  //
  // Aug 22 2019: incoude dmHost term based on initial hostPar values.
  // Apr 02 2020: for BIASCOR_MU option, a,b = user input p1,p2  

  double z, a, b, g, M0, mB, x1, c, logmass, dmHost, hostPar[10];
  double zTrue, muFit, muTrue, muz, muDif ;
  double dlz, dlzTrue, dmu  ;
  //  char fnam[] = "muresid_biasCor" ;

  // ----------------- BEGIN ----------------

  bool DOBIAS_MU = ( INPUTS.opt_biasCor & MASK_BIASCOR_MU     ) ;

  M0     = INPUTS.nommag0 ;

  if (DOBIAS_MU) {
    a = INPUTS.parval[IPAR_ALPHA0];
    b = INPUTS.parval[IPAR_BETA0];
  } 
  else{
    a        = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt] ;
    b        = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt] ;
  }

  g        = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt] ;
  z        = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt] ;
  logmass  = (double)INFO_BIASCOR.TABLEVAR.logmass[ievt];
  mB       = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_mB][ievt] ; 
  x1       = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_x1][ievt] ;
  c        = (double)INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][ievt] ;
  zTrue    = (double)INFO_BIASCOR.TABLEVAR.SIM_ZCMB[ievt];
  muTrue   = (double)INFO_BIASCOR.TABLEVAR.SIM_MU[ievt];
  
  dmHost = 0.0 ;
  if ( INPUTS.USE_GAMMA0 ) {
    hostPar[0] = INPUTS.parval[IPAR_GAMMA0];   // user-initial values
    hostPar[1] = INPUTS.parval[IPAR_GAMMA0+1];
    hostPar[2] = INPUTS.parval[IPAR_LOGMASS_CEN];
    hostPar[3] = INPUTS.parval[IPAR_LOGMASS_TAU];
    dmHost     = get_gammadm_host(z,logmass,hostPar);
  }

  // need true MU at observed redshift z. We don't have the biasCor
  // cosmology parameters, so we'll use a derivative to compute
  //     muTrue = SIM_DLMAG + (z-zTrue)*dmu/dz,
  // where dmu/dz is computed from the SALT2mu-reference cosmology.
  // Assumption here is that correction to SIM_DLMAG is very small
  // so that using approximate cosmology for dmu/dz won't matter.


  // use measured z including zpec
  dlz      = cosmodl_forFit(z, INPUTS.COSPAR); 

  // use true zCMB
  dlzTrue  = cosmodl_forFit(zTrue, INPUTS.COSPAR); 

  dmu    = 5.0*log10(dlz/dlzTrue) ;
  muz    = muTrue + dmu ;
  muFit  = mB + a*x1 - b*c + dmHost - M0 ; 
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

  int DOCOVFIT     = ( (maskCov & USEMASK_BIASCOR_COVFIT)>0 ) ;
  int DOCOVINT     = ( (maskCov & USEMASK_BIASCOR_COVINT)>0 ) ;

  int IDSAMPLE, iz, ia, ib, ig, OPTMASK;
  double z, zerr, a, b, gDM;
  double COVTOT[NLCPAR][NLCPAR] ;
  double COVINT[NLCPAR][NLCPAR] ;
  double COVTMP, muErrsq ;
  int  j0, j1 ;
  char *name ;
  //  char fnam[] = "muerrsq_biasCor" ;
  
  // ------------- BEGIN -----------


  IDSAMPLE = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt] ;
  z        = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt] ;
  zerr     = (double)INFO_BIASCOR.TABLEVAR.zhderr[ievt] ;
  a        = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt] ;
  b        = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt] ;
  gDM      = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt] ;
  ia       = (int)INFO_BIASCOR.IA[ievt] ;
  ib       = (int)INFO_BIASCOR.IB[ievt] ;
  ig       = (int)INFO_BIASCOR.IG[ievt] ;
  name     = INFO_BIASCOR.TABLEVAR.name[ievt];
  iz       = IBINFUN(z,  &INPUTS.BININFO_z, 0, "" );
  if ( DOCOVINT ) 
    { get_COVINT_biasCor(IDSAMPLE,z,a,b,gDM, COVINT); } // return COVINT
   

  //  check option to include intrinsic scatter matrix

  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0;j1<NLCPAR;j1++ ) { 

      COVTMP = 0.0 ;

      if ( DOCOVFIT ) 
	{ COVTMP += (double)INFO_BIASCOR.TABLEVAR.covmat_fit[ievt][j0][j1]; }
      
      if ( DOCOVINT )	
	{ COVTMP += COVINT[j0][j1]; }

      COVTOT[j0][j1] = COVTMP ;

    } 
  }


  // fix matrix if it has bad eigenvalues (just like for data)
  OPTMASK=0;
  update_covMatrix( name, OPTMASK, NLCPAR, COVTOT, EIGMIN_COV, istat_cov ); 
  
  // compute error on distance, including covariances
  muErrsq = fcn_muerrsq(name, a, b, gDM, COVTOT, z, zerr, 0 );

  return(muErrsq);

} // end muerrsq_biasCor


// ===============================================================
void get_COVINT_biasCor(int IDSAMPLE, double z, 
			double alpha, double beta, double gammadm,
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

  int  Na       = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  Ng       = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
  int j0, j1, iz, ia, ib, ig ;
  int LDMP = 0 ;
  int DUMPFLAG_abWGT = LDMP ;
  double WGT, COVTMP ;
  char fnam[] = "get_COVINT_biasCor" ;

  // -------------- BEGIN ----------------

  // check for valid arguments
  if ( IDSAMPLE < 0 || z<0.0 || alpha < 1.0E-5 || beta < 1.0E-5 ) {
    sprintf(c1err, "Invalid argument(s)");
    sprintf(c2err, "IDSAMPLE=%d  z=%.4f  a=%.4f  b=%.4f  g=%.3f", 
	    IDSAMPLE, z, alpha, beta, gammadm );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  zero_COV(COVINT);

  if ( LDMP ) {
    printf("\n xxx ------------------------ \n");
    printf(" xxx IDSAMP=%d  z=%.4f  a=%.4f  b=%.4f  g=%.3f\n",
	   IDSAMPLE, z, alpha, beta, gammadm );
    fflush(stdout);
  }

  // convert redshift to index
  iz  = IBINFUN(z, &CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 0, "" );

  // get alpha,beta interp weights
  get_INTERPWGT_abg(alpha, beta, gammadm, DUMPFLAG_abWGT, 
		    &INTERPWGT, fnam );

  for(j0=0; j0<NLCPAR; j0++ )  { 
    for(j1=0; j1<NLCPAR; j1++ )  { 

      for(ia=0; ia<Na; ia++ ) {
	for(ib=0; ib<Nb; ib++ ) {
	  for(ig=0; ig<Ng; ig++ ) {
	    WGT    = INTERPWGT.WGT[ia][ib][ig] ; 
	    COVTMP = INFO_BIASCOR.COVINT[IDSAMPLE][iz][ia][ib][ig].VAL[j0][j1];
	    COVINT[j0][j1] += ( WGT * COVTMP );
	  } // end ig
	} // end ib
      } // end ia

      if ( LDMP ) {
	printf(" xxx COVINT[%d][%d] = %f \n",
	       j0, j1, COVINT[j0][j1] ); fflush(stdout);
      }

    }  // j1
  }  // j0

  if ( LDMP ) {  debugexit(fnam); } 

  return ;

} // end get_COVINT_biasCor
 

// =======================================
void  write_COVINT_biasCor(void) {

  // Open text file [prefix].COVINT and each COV matrix
  // to table row  in text file.

  int  NSAMPLE  = NSAMPLE_BIASCOR ;
  int  Na       = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  Ng       = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  int  idsample, Nz, ia ,ib, ig, iz, NEVT, k0 ;
  int  NVAR, IROW=0 ;
  FILE *fp;
  double z, a, b, g, cov, rho, covDiag, sigma[NLCPAR] ;
  char outFile[200], cline[200] ;
  //  char fnam[] = "write_COVINT_biasCor" ;

  // ------------ BEGIN ----------------

  sprintf(outFile,"%s.COVINT", INPUTS.PREFIX );
  fp = fopen(outFile,"wt");

  fprintf(fp, "# REDCOV01 = COVINT(mB,x1)/(sigma_mB*sigma_x1) \n");
  fprintf(fp, "# REDCOV02 = COVINT(mB,c )/(sigma_mB*sigma_c ) \n");
  fprintf(fp, "# REDCOV12 = COVINT(x1,c )/(sigma_x1*sigma_c ) \n");
  fprintf(fp, "# \n");

  NVAR=13 ;
#ifdef TEXTFILE_NVAR
  fprintf(fp, "NVAR: %d\n", NVAR);
#endif
  fprintf(fp, "VARNAMES: "
	  "ROW NEVT IDSAMPLE z alpha beta gammaDM"    // 7
	  "sigma_mB sigma_x1 sigma_c "            // 3
	  "REDCOV01  REDCOV02 REDCOV12 "          // 3
	  "\n" );
  
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    Nz       = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    for(iz=0; iz < Nz; iz++ ) {
      for(ia=0; ia < Na; ia++ ) {
	for(ib=0; ib < Nb; ib++ ) {
	  for(ig=0; ig < Ng; ig++ ) {
	    IROW++ ; 
	    NEVT = INFO_BIASCOR.NEVT_COVINT[idsample][iz][ia][ib][ig]; 
	    z    = CELLINFO_BIASCOR[idsample].BININFO_z.avg[iz] ;
	    a    = INFO_BIASCOR.BININFO_SIM_ALPHA.avg[ia] ;
	    b    = INFO_BIASCOR.BININFO_SIM_ALPHA.avg[ib] ;
	    g    = INFO_BIASCOR.BININFO_SIM_GAMMADM.avg[ig] ;
	    cline[0] = 0 ;

	    sprintf(cline,"ROW: "
		    "%3d %4d %d %.3f %4.2f %4.2f %4.2f",
		    IROW, NEVT, idsample, z, a, b, g );

	    for(k0=0; k0<NLCPAR; k0++ ) {
	      covDiag=INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[k0][k0];
	      sigma[k0] = sqrt(covDiag);
	      sprintf(cline,"%s %7.4f ", cline, sigma[k0] );
	    }
	  
	    // REDCOV01
	    cov = INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[0][1] ;
	    rho = cov / (sigma[0]*sigma[1]);
	    sprintf(cline,"%s %7.4f ", cline, rho );

	    // REDCOV02
	    cov = INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[0][2] ;
	    rho = cov / (sigma[0]*sigma[2]);
	    sprintf(cline,"%s %7.4f ", cline, rho );

	    // REDCOV12
	    cov = INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[1][2] ;
	    rho = cov / (sigma[1]*sigma[2]);
	    sprintf(cline,"%s %7.4f ", cline, rho );
	    
	    fprintf(fp,"%s\n", cline);
	    
	  } 
	}
      }
    }
  }

  fclose(fp);
  return ;

} // end void   write_COVINT_biasCor


// ===============================================
void  dump_COVINT_biasCor(int idsample, int iz, int ia, int ib, int ig) {

  // Screen dump COVINT for the input indices.

  int k0, k1, NEVT ;
  double covDiag, cov, sigma0, sigma1, rho ;
  double zlo = CELLINFO_BIASCOR[idsample].BININFO_z.lo[iz] ;
  double zhi = CELLINFO_BIASCOR[idsample].BININFO_z.hi[iz] ;
  char *parName ;
  COV_DEF *COV_LOCAL;
  char fnam[] = "dump_COVINT_biasCor" ;

  // ------------- BEGIN --------------

  //INFO_BIASCOR.COVINT[idsample][iz][ia][ib].VAL[k0][k0] ; }

  NEVT = INFO_BIASCOR.NEVT_COVINT[idsample][iz][ia][ib][ig] ; 
  COV_LOCAL = &INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig] ;
  
  printf("\n");
  printf("# ----------------------------------------- \n");
  printf("# %s: \n", fnam);
  printf("#   idsample=%d  iz=%d(%.3f-%.3f)  ia,ib,ig=%d,%d,%d   NEVT=%d\n", 
	 idsample, iz,zlo,zhi, ia,ib,ig, NEVT);


  for(k0=0; k0 < NLCPAR; k0++ ) {   
    covDiag = (*COV_LOCAL).VAL[k0][k0] ;     
    sigma0  = sqrt(covDiag);
    parName = BIASCOR_NAME_LCFIT[k0] ;

    printf("\t sigma(%2s) = %7.4f | ", parName, sigma0 );
    for(k1=0; k1 < NLCPAR; k1++ ) {
      // print reduced COV here ...
      covDiag = (*COV_LOCAL).VAL[k1][k1] ;
      cov     = (*COV_LOCAL).VAL[k0][k1] ;
      sigma1  = sqrt(covDiag);
      rho     = cov/(sigma0*sigma1);
      printf(" %8.4f ", rho);
    }
    printf("\n");
  }

  return ;

} // end dump_COVINT_biasCor


// =====================================
void init_COVINT_biasCor(void) {

  // Jan 2018
  // Compute and store intrinsic scatter matrix for biasCor.
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
  int  NROW_TOT = INFO_BIASCOR.TABLEVAR.NSN_ALL ; 
  int  Na       = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  Nb       = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  Ng       = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  int  Nz, idsample, iz, ia, ib, ig, ievt, ipar, ipar2 ; 
  double sigInt, COV, z ;
  char fnam[] = "init_COVINT_biasCor" ;

  //  int DOTEST_SIGINT_ONLY = 1 ; // traditional sigma_int model (debug only)

  // ------------- BEGIN ---------------

  sprintf(BANNER,"%s:", fnam);
  print_banner(BANNER); 

  printf("\t Compute Intrinsic Matrix (COVINT) from BiasCor Sample.\n");
  printf("\t COVINT -> bins of IDSAMPLE, Redshift, Alpha, Beta, GammaDM\n");
  fflush(stdout);

  // zero out entire COV matrix.
  // For traditional sigma_int model (IDEAL_COVINT=0), set COV[0][0]
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    zero_COV(INFO_BIASCOR.COVINT_AVG[idsample].VAL) ;
    Nz       = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
    for(iz=0; iz < Nz; iz++ ) {
      for(ia=0; ia < Na; ia++ ) {
	for(ib=0; ib < Nb; ib++ ) {
	  for(ig=0; ig < Ng; ig++ ) {
	    zero_COV(INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL) ;
	    INFO_BIASCOR.NEVT_COVINT[idsample][iz][ia][ib][ig] = 0 ;
	    if ( DO_COV00_ONLY ) {
	      sigInt = INFO_BIASCOR.SIGINT_ABGRID[idsample][ia][ib][ig]; 
	      COV    = (sigInt * sigInt) ;
	      INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[0][0] = COV ;
	    }

	  } // end ig
	} // end ib
      } // and ia
    } // ene iz
  }  // end idsample


  if ( DO_COV00_ONLY ) { return ; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  // If we get here, compute full COV matrix in bins of 
  // IDSAMPLE,z,a,b
  //   COV(x,y) = sum[(x-xtrue)*(y-ytrue) ] / N


  int NBIASCOR_IDEAL=0, NBIASCOR_CUTS=0 ;
  double tmpVal, tmpVal2, x0_IDEAL, mB_IDEAL ;

  for(ievt=0; ievt < NROW_TOT; ievt++ ) {

    // apply selection
    if ( INFO_BIASCOR.TABLEVAR.CUTMASK[ievt] ) { continue; }
    NBIASCOR_CUTS++ ;

    // check for valid IDEAL fit params 
    tmpVal = INFO_BIASCOR.TABLEVAR.fitpar_ideal[INDEX_mB][ievt];
    if ( tmpVal > 900.0 ) { continue ; }
    
    // mB_IDEAL is corrupt, so hack fix based on x0
    x0_IDEAL = (double)INFO_BIASCOR.TABLEVAR.x0_ideal[ievt] ;
    mB_IDEAL = -2.5*log10(x0_IDEAL);
    INFO_BIASCOR.TABLEVAR.fitpar_ideal[INDEX_mB][ievt] = (float)mB_IDEAL;

    // reject mB outliers:
    tmpVal = 
      INFO_BIASCOR.TABLEVAR.fitpar_ideal[INDEX_mB][ievt] -
      INFO_BIASCOR.TABLEVAR.SIM_FITPAR[INDEX_mB][ievt] ;
    if ( fabs(tmpVal) > 1.0 ) { continue ; }


    idsample = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt];
    ia       = (int)INFO_BIASCOR.IA[ievt] ; // true alpha index
    ib       = (int)INFO_BIASCOR.IB[ievt] ; // true beta index
    ig       = (int)INFO_BIASCOR.IG[ievt] ; // true gamma DM
    iz       = (int)INFO_BIASCOR.IZ[ievt] ; // true zcmb index
    z        = (int)INFO_BIASCOR.TABLEVAR.SIM_ZCMB[ievt] ;

    INFO_BIASCOR.NEVT_COVINT[idsample][iz][ia][ib][ig]++ ;

    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      for(ipar2=ipar; ipar2 < NLCPAR; ipar2++ ) {

	tmpVal = 
	  INFO_BIASCOR.TABLEVAR.fitpar_ideal[ipar][ievt] -
	  INFO_BIASCOR.TABLEVAR.SIM_FITPAR[ipar][ievt] ;

	tmpVal2 = 
	  INFO_BIASCOR.TABLEVAR.fitpar_ideal[ipar2][ievt] -
	  INFO_BIASCOR.TABLEVAR.SIM_FITPAR[ipar2][ievt] ;

	/* xxxxxxxxxxxx
	if ( ievt < 100 && ipar==0 && ipar2==0 ) {
	  printf(" xxx diff[%d][%d] = %f, %f \n",
		 ipar, ipar2, tmpVal, tmpVal2 );
	}
	xxxxxxxxxxx*/

	INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[ipar][ipar2] 
	  += (tmpVal * tmpVal2) ;

	INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[ipar2][ipar] =
	  INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL[ipar][ipar2] ;
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
	  for(ig=0; ig < Ng; ig++ ) {
	  
	    N = INFO_BIASCOR.NEVT_COVINT[idsample][iz][ia][ib][ig]; 
	    if ( N > 0 ) {
	      XNINV = 1.0/(double)N;
	      scale_COV(XNINV,INFO_BIASCOR.COVINT[idsample][iz][ia][ib][ig].VAL);
	    }
	    
	    
	    if ( iz < -4 ) { dump_COVINT_biasCor(idsample,iz,ia,ib,ig); }
	  } // end ig
	}
      }
    }
  }

  
  write_COVINT_biasCor();


  return ;

} // end init_COVINT_biasCor


// ======================================================
void  init_sigInt_biasCor_SNRCUT(int IDSAMPLE) {

  // Estimate sigInt for biasCor sample using very high SNR subset.
  // Jun 17 2016: compute sigInt in each alpha & beta bin. 
  // Jan 25,2018
  //   + set SKIPZFLAG_BIASCOR flag
  //   + abort if NUSE[ia][ib] < 10
  // Oct 08 2018: 
  //   + new input IDSAMPLE to allow option of sigint(IDSAMPLE)
  //
  // Jun 2 2019:
  //   MUDIF and MUERSQ are now global instead of local in 
  //   case they can take too much memory for local.
  //
  // Jul 16 2019: add gamma (ig) dimension
  // Dec 11 2019: 
  //   + abort if muErrsq < 0 or if NUSE exceeds bound.
  //   + return immediately of 1D+5DCUT is set.
  //
  int  DO_SIGINT_SAMPLE = ( INPUTS.opt_biasCor & MASK_BIASCOR_SIGINT_SAMPLE ) ;
  int  DOCOR_1D5DCUT    = ( INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT );
  int  NROW_TOT, NROW_malloc, i, istat_cov, NCOVFIX, ia, ib, ig, MEMD, cutmask ;
  int  LDMP = 0 ;

  double muErrsq, muErr, muDif, muOff, SNRMAX, sigInt, tmp1, tmp2 ;
  double SUMDIF[MXa][MXb][MXg],SUMDIFSQ[MXa][MXb][MXg],SUMERRSQ[MXa][MXb][MXg];
  int    NUSE[MXa][MXb][MXg];
  int    NSNRCUT, NTMP, NBINa, NBINb, NBINg, *ptr_CUTMASK ;
  short int *ptr_IDSAMPLE;

  double  *MUDIF[MXa][MXb][MXg];
  double  *MUERRSQ[MXa][MXb][MXg];

  float  *ptr_SNRMAX;
  double SIGINT_AVG, SIGINT_ABGRID[MXa][MXb][MXg] ; 

  char *NAME;
  char fnam[]   = "init_sigInt_biasCor_SNRCUT" ;

  // ------------------- BEGIN -------------------

  if ( DOCOR_1D5DCUT ) { return ; }

  printf("\n");

  NROW_TOT     = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  NBINa        = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  NBINb        = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  NBINg        = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  ptr_CUTMASK  = INFO_BIASCOR.TABLEVAR.CUTMASK;
  ptr_IDSAMPLE = INFO_BIASCOR.TABLEVAR.IDSAMPLE ;
  ptr_SNRMAX   = INFO_BIASCOR.TABLEVAR.snrmax ;

  // ---------------------------------------
  // check option for user to fix sigmb_biascor
  sigInt = INPUTS.sigint_biasCor ;
  if ( sigInt >= 0.0 ) {
    if ( IDSAMPLE == 0 ) 
      { printf(" sigInt -> %.3f from user input sigmb_biascor key\n", sigInt);}
    fflush(stdout);
    SIGINT_AVG = sigInt ;
    for(ia=0; ia < NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {
	for(ig=0; ig < NBINg; ig++ ) {
	  SIGINT_ABGRID[ia][ib][ig] = sigInt ;
	}  
      }
    }
    goto  LOAD_GLOBALS ;
  }


  // ---------------------------------------
  // check option to compute sigInt for each IDSAMPLE
  if ( DO_SIGINT_SAMPLE ) {
    printf(" %s for IDSAMPLE=%d (%s): \n",
	   fnam, IDSAMPLE,  SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }
  else if ( IDSAMPLE == 0 ) {
    printf(" %s for all IDSAMPLEs combined: \n", fnam);
  }
  else {
    // do nothing
  } 

  fflush(stdout) ;

  // ---------------------------------------
  // quick pass with SNR cut to estimate size for malloc  
  NROW_malloc = 0 ;
  for(i=0; i < NROW_TOT; i++ ) {

    // apply selection without BIASCOR-z cut (but keep global z cut)
    cutmask  = ptr_CUTMASK[i] ;
    cutmask -= (cutmask & CUTMASK_LIST[CUTBIT_zBIASCOR]);
    if ( cutmask ) { continue; }

    SNRMAX = (double)(ptr_SNRMAX[i]) ;
    if ( SNRMAX < INPUTS.snrmin_sigint_biasCor) { continue ; }
    NROW_malloc++ ;
  }


  NROW_malloc += 100; // safety margin
  MEMD = NROW_malloc * sizeof(double) ;
  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      for(ig=0; ig < NBINg; ig++ ) {
	MUDIF[ia][ib][ig]   = (double*) malloc ( MEMD );
	MUERRSQ[ia][ib][ig] = (double*) malloc ( MEMD ) ;
	NUSE[ia][ib][ig]    = NCOVFIX = 0 ;
	SUMDIF[ia][ib][ig]  = SUMDIFSQ[ia][ib][ig] = 0.0 ;
	SUMERRSQ[ia][ib][ig] = 0.0 ;
      }
    }
  }

  // skip z-cut to ensure events with SNR>60.
  // This allows, for example, doing BBC fits only at high-z.

  NSNRCUT=0;

  for(i=0; i < NROW_TOT; i++ ) {

    // get variables defining GRID map
    SNRMAX = (double)(ptr_SNRMAX[i]) ;
    if ( SNRMAX < INPUTS.snrmin_sigint_biasCor ) { continue ; }

    // apply selection without BIASCOR-z cut (but keep global z cut)
    cutmask  = ptr_CUTMASK[i] ;
    cutmask -= (cutmask & CUTMASK_LIST[CUTBIT_zBIASCOR]);

    if ( cutmask ) { continue; }

    // check option for IDSAMPLE-dependent sigint (Oct 2018)
    if ( DO_SIGINT_SAMPLE && (IDSAMPLE != ptr_IDSAMPLE[i]) ) 
      { continue ; }
    

    NSNRCUT++ ;
    NAME   = INFO_BIASCOR.TABLEVAR.name[i] ;
    muDif  = muresid_biasCor(i); 

    // compute error with no intrinsic scatter (just use data COVFIT)
    muErrsq = muerrsq_biasCor(i, USEMASK_BIASCOR_COVFIT, &istat_cov) ;

    if ( muErrsq < 0.0 ) {
      sprintf(c1err,"Invalid muErrsq[%d] = %f < 0  SNID=%s",
	      i, muErrsq, NAME );
      sprintf(c2err,"muDif=%le ", muDif);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);   
    }

    if ( isnan(muErrsq) || isnan(muDif) ) {
      sprintf(c1err,"isnan trap for SNID = %s (irow=%d)", NAME, i);
      sprintf(c2err,"muErrsq=%le, muDif=%le ", muErrsq, muDif);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);   
    }


    if ( istat_cov < 0 ) { NCOVFIX++ ; }    

    ia = (int)INFO_BIASCOR.IA[i];
    ib = (int)INFO_BIASCOR.IB[i];
    ig = (int)INFO_BIASCOR.IG[i]; 

    SUMDIF[ia][ib][ig]   +=  muDif ;
    SUMDIFSQ[ia][ib][ig] += (muDif*muDif);
    SUMERRSQ[ia][ib][ig] +=  muErrsq ;

    NTMP = NUSE[ia][ib][ig];
    MUDIF[ia][ib][ig][NTMP]   = muDif ;
    MUERRSQ[ia][ib][ig][NTMP] = muErrsq ;
    NUSE[ia][ib][ig]++ ;

    // xxxxxxxx

    if ( NUSE[ia][ib][ig] >= NROW_malloc ) {
      sprintf(c1err,"NUSE[ia,ib,ig=%d,%d,%d] = %d exceeds malloc size",
	      ia, ib, ig, NUSE[ia][ib][ig] );
      sprintf(c2err,"NROW_malloc = %d", NROW_malloc);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }
    // xxxxxxxxxx

  } // end loop over biasCor sample


  // -------------------------------------------------
  double XN, SQRMS, RMS, pull, sigInt_bin, sigTmp_lo, sigTmp_hi ;
  double sigTmp, sqsigTmp, sigInt_store[20], rmsPull_store[20] ;
  double sumdif, sumdifsq ;
  int    NBIN_SIGINT, DOPRINT,  OPT_INTERP=1 ;
  double ONE = 1.0 ;

  SIGINT_AVG = 0.0 ;

  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {    
      for(ig=0; ig < NBINg; ig++ ) {

	XN     = (double)NUSE[ia][ib][ig] ;
	if ( XN < 10.0 ) {
	  sprintf(c1err,"NUSE[ia,ib,ig=%d,%d,%d] = %d is too small.",
		  ia, ib, ig, NUSE[ia][ib][ig] );
	  sprintf(c2err,"Check biasCor events with SNR> %.1f",
		  INPUTS.snrmin_sigint_biasCor ) ;
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
	}
	
	tmp1   = SUMDIF[ia][ib][ig]/XN ;
	tmp2   = SUMDIFSQ[ia][ib][ig]/XN ;
	SQRMS  = tmp2 - tmp1*tmp1 ;
	muOff  = SUMDIF[ia][ib][ig]/XN ; 

	// compute approx sigInt using average distance error
	muErrsq = SUMERRSQ[ia][ib][ig]/XN ;
	
	if ( muErrsq > SQRMS ) { 
	  muErr = sqrt(muErrsq);        RMS = sqrt(SQRMS);
	  sprintf(c1err,"Avg computed muErr > RMS(mures) for SNRMAX>%.0f",
		  INPUTS.snrmin_sigint_biasCor );
	  sprintf(c2err,"<muErr>=%f   RMS(mures)=%f   ia,ib,ig = %d,%d,%d", 
		  muErr, RMS, ia, ib, ig );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
	}

	sigInt = sqrt(SQRMS - muErrsq) ;  // approx sigInt
	
	if ( LDMP ) {
	  printf(" xxx ia,ib,ig=%d,%d,%d  N = %d \n", ia,ib,ig, (int)XN );
	  printf(" xxx --> RMS(mu)=%.4f  muErr=%.4f  muOff=%.4f  sigInt=%.4f\n",
		 sqrt(SQRMS), sqrt(muErrsq),  muOff, sigInt );        
	  fflush(stdout);
	}

	// Loop over trial sigInt values . . . . 
	NBIN_SIGINT = 0 ;
	sigInt_bin = 0.01 ;
	sigTmp_lo  = sigInt - (4.*sigInt_bin) - 1.0E-7 ;
	sigTmp_hi  = sigInt + (8.*sigInt_bin) ;
	
	// start with largest sigInt and decrease so that RMS is increasing
	// for the interp function below
	for(sigTmp = sigTmp_hi; sigTmp >= sigTmp_lo; sigTmp -= sigInt_bin ) {
	  sumdif = sumdifsq = 0.0 ;
	  if ( sigTmp < 0.0 ) { continue ; }
	  sqsigTmp = sigTmp * sigTmp ;
	  
	  for(i=0; i < NUSE[ia][ib][ig]; i++ ) {
	    muDif     = MUDIF[ia][ib][ig][i] - muOff ;
	    muErrsq   = MUERRSQ[ia][ib][ig][i] ;  
	    pull      = muDif/sqrt(muErrsq + sqsigTmp);
	    if ( isnan(pull) ) {
	      sprintf(c1err,"crazy pull = %f  ia,ib,ig=%d,%d,%d",
		      pull, ia, ib, ig );
	      sprintf(c2err,"i=%d  muDif=%f muErrsq = %f", i, muDif, muErrsq);
	      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	    }
	    sumdif   += pull ;
	    sumdifsq += (pull*pull);
	  }
	  
	  tmp1=sumdif/XN;    tmp2=sumdifsq/XN;
	  SQRMS  = tmp2 - tmp1*tmp1 ;
	  rmsPull_store[NBIN_SIGINT] = sqrt(SQRMS);
	  sigInt_store[NBIN_SIGINT]  = sigTmp ;
	  
       
	  if ( LDMP ) {
	    printf("\t xxx ia,ib,ig=%d,%d,%d : "
		   "sigTmp(%d)=%.3f --> RMS(pull)=%.3f \n",
		   ia, ib, ig, NBIN_SIGINT, sigInt_store[NBIN_SIGINT], 
		   rmsPull_store[NBIN_SIGINT])  ;    
	    printf("\t xxx XN=%f, sumdif=%le, sumdifsq=%le \n",
		   XN, sumdif, sumdifsq );
	    fflush(stdout);
	  }
	  
	  NBIN_SIGINT++ ;
	}  // end sigTmp loop
	
      
	// interpolate sigInt vs. rmsPull at rmsPull=1
	sigInt = interp_1DFUN(OPT_INTERP, ONE, NBIN_SIGINT,
			      rmsPull_store, sigInt_store, fnam);
	
	// load SIGINT value
	SIGINT_ABGRID[ia][ib][ig] = sigInt ;
	SIGINT_AVG += sigInt ;

      if ( IDSAMPLE==0 || DO_SIGINT_SAMPLE) { DOPRINT=1; } else { DOPRINT=0; }
      if ( DOPRINT ) {
	printf("\t sigInt[ia,ib,ig=%d,%d,%d] = %.4f "
	       "(%d events with SNR>%.0f) \n", 
	       ia, ib, ig, sigInt, NUSE[ia][ib][ig], 
	       INPUTS.snrmin_sigint_biasCor);
	fflush(stdout);
      }
      
      } // ig
    } // ib
  } // ia

 
  SIGINT_AVG /= (double)( NBINa*NBINb*NBINg) ;

  
  // free memory
  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      for(ig=0; ig < NBINg; ig++ ) {
	free(MUDIF[ia][ib][ig] )   ;
	free(MUERRSQ[ia][ib][ig] ) ;  
      }
    }
  } 
 

 LOAD_GLOBALS:

  INFO_BIASCOR.SIGINT_AVG = SIGINT_AVG;
  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      for(ig=0; ig < NBINg; ig++ ) {
	INFO_BIASCOR.SIGINT_ABGRID[IDSAMPLE][ia][ib][ig] = 
	  SIGINT_ABGRID[ia][ib][ig] ; 
      }
    }
  }
  

  return ;

} // end init_sigInt_biasCor_SNRCUT


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

  short int *ptr_IDSAMPLE, idsample ;
  int       *ptr_CUTMASK  ;
  int NROW, irow, isp, cutmask, MEMI, NBIASCOR ;
  char fnam[] = "makeSparseList_biasCor" ;

  // ---------- BEGIN ------------

  NROW         = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptr_CUTMASK  = INFO_BIASCOR.TABLEVAR.CUTMASK ; // set pointer
  ptr_IDSAMPLE = INFO_BIASCOR.TABLEVAR.IDSAMPLE ;   

  sprintf(BANNER," %s: make sparse row-list for each biasCor sample.",fnam);
  print_banner(BANNER);    

  // allocate memory for sparse irow list
  for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) {
    NBIASCOR = SAMPLE_BIASCOR[idsample].NSN[EVENT_TYPE_BIASCOR] ;
    MEMI     = NBIASCOR * sizeof(int); 
    SAMPLE_BIASCOR[idsample].IROW_CUTS = (int*)malloc(MEMI);
    SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS = 0 ;
  }


  for(irow=0; irow < NROW; irow++ ) {
    cutmask  = ptr_CUTMASK[irow];
    idsample = ptr_IDSAMPLE[irow] ;

    // apply selection cuts for making biasCor map
    if ( cutmask ) { continue; }

    isp = SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS ;
    if ( isp < 0 || isp > NROW ) {
      sprintf(c1err,"Invalid isp=%d (must be 0 to %d)", isp, NROW);
      sprintf(c2err,"Something is messed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
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

  int irow, isp, J1D, ipar,  *NperCell ;
  double WGT, z, m, fitPar[NLCPAR];
  double *SUM_WGT_5D, *SUM_z_5D,  *SUM_m_5D, *SUM_FITPAR_5D[NLCPAR];
  float *ptr_z, *ptr_m, *ptr_fitPar[NLCPAR];
  char fnam[] = "makeMap_binavg_biasCor" ;

  // --------------- BEGIN ---------------

  ptr_z  = INFO_BIASCOR.TABLEVAR.zhd ;
  ptr_m  = INFO_BIASCOR.TABLEVAR.logmass ; 

  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    { ptr_fitPar[ipar]  = INFO_BIASCOR.TABLEVAR.fitpar[ipar] ; }  

  NperCell    = (int*)    malloc(MEMI) ;
  SUM_z_5D    = (double*) malloc(MEMD) ;
  SUM_m_5D    = (double*) malloc(MEMD) ;
  SUM_WGT_5D  = (double*) malloc(MEMD) ;

  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    { SUM_FITPAR_5D[ipar] = (double*) malloc(MEMD) ; }

  for(J1D = 0; J1D < NCELL; J1D++ ) {
    NperCell[J1D] = 0 ;
    SUM_WGT_5D[J1D] = SUM_z_5D[J1D] = SUM_m_5D[J1D] = 0.0 ;
    for(ipar=0; ipar < NLCPAR; ipar++ )  { SUM_FITPAR_5D[ipar][J1D] = 0.0 ; }  

    // init bin avg to nonsense values
    CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D] = 9999. ;
    CELLINFO_BIASCOR[IDSAMPLE].AVG_m[J1D] = 9999. ;
    for(ipar=0; ipar < NLCPAR; ipar++ ) 
      { CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar][J1D] = 9999. ;  }  
  }


  // loop over biasCor sample ... xyz
  for(isp=0; isp < NROW; isp++ ) {

    irow = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;

    WGT = WGT_biasCor(1,irow,fnam) ;  // WGT=1/muerr^2
    J1D = J1D_biasCor(irow,fnam);     // 1D index
    NperCell[J1D]++ ;
    SUM_WGT_5D[J1D] += WGT ;

    z   = (double)ptr_z[irow] ;
    m   = (double)ptr_m[irow] ;
    SUM_z_5D[J1D]   += (WGT * z) ;
    SUM_m_5D[J1D]   += (WGT * m) ;

    for(ipar=0; ipar < NLCPAR ; ipar++ ) {
      fitPar[ipar]  = (double)ptr_fitPar[ipar][irow] ;
      SUM_FITPAR_5D[ipar][J1D] += ( WGT * fitPar[ipar] ) ;
    }

  } // end loop over biasCor rows


  // ------------------
  // loop over cells and compute avg
  for(J1D = 0; J1D < NCELL; J1D++ ) {
    WGT = SUM_WGT_5D[J1D] ;
    if ( WGT < 1.0E-9 ) { continue ; }
    CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D]  =  SUM_z_5D[J1D] / WGT ;
    CELLINFO_BIASCOR[IDSAMPLE].AVG_m[J1D]  =  SUM_m_5D[J1D] / WGT ;
    for(ipar=0; ipar < NLCPAR; ipar++ ) { 
      CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[ipar][J1D] = 
	SUM_FITPAR_5D[ipar][J1D]/WGT;
    }
  }


  // ---------------------------------------
  // free temp memory
  free(NperCell); free(SUM_WGT_5D); free(SUM_z_5D) ; free(SUM_m_5D) ;
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

  int NROW, iz, i, CUTMASK, *ptr_CUTMASK ;
  double SUM_WGT[MXz], SUM_z[MXz];
  double WGT, z, zM0, *ptr_zM0 ;
  float *ptr_z;
  char fnam[] = "calc_zM0_biasCor" ;

  // --------- BEGIN -------------

  NROW        = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptr_CUTMASK = INFO_BIASCOR.TABLEVAR.CUTMASK;
  ptr_z       = INFO_BIASCOR.TABLEVAR.zhd ;
  ptr_zM0     = INFO_BIASCOR.zM0 ; 

  for(iz=0; iz < NBINz; iz++ )  { SUM_WGT[iz] = SUM_z[iz] = 0.0 ; }

  // loop over biasCor sample
  for(i=0; i < NROW; i++ ) {

    // apply selection cuts for making biasCor map
    CUTMASK = ptr_CUTMASK[i] ;
    if ( CUTMASK ) { continue; }

    WGT = WGT_biasCor(1,i,fnam) ;  // WGT=1/muerr^2

    z   = (double)(ptr_z[i]) ;
    iz  = IBINFUN(z,  &INPUTS.BININFO_z, 0, "" );

    SUM_z[iz]    += (WGT * z) ;
    SUM_WGT[iz]  += WGT ;
  }
  
  // store avg redshift for writing to M0 outFile
  for(iz=0; iz < NBINz; iz++ ) {
    WGT = SUM_WGT[iz] ;
    if ( WGT > 0.0 ) 
      { zM0 = SUM_z[iz] / WGT ; }
    else
      { zM0 = INPUTS.BININFO_z.avg[iz] ; }

    ptr_zM0[iz] = zM0;
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

  int iz, i, CUTMASK, NSN_DATA ;
  double SUM_WGT[MXz], SUM_z[MXz];
  double SUM_mu[MXz], SUM_mumodel[MXz];
  double WGT, mu, mumodel, muerr, muerrsq, z, zerr ;
  double *ptr_zM0;
  char *name ;
  char fnam[] = "calc_zM0_data" ;
  int  LDMP   = 0 ;
  // --------- BEGIN -------------
  
  NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ; 
  ptr_zM0  = INFO_BIASCOR.zM0 ;

  for(iz=0; iz < NBINz; iz++ ) { 
    SUM_WGT[iz] = SUM_z[iz] = SUM_mu[iz] = SUM_mumodel[iz] = 0.0 ; 
    FITRESULT.zM0[iz] = ptr_zM0[iz] ;
  }


  // - - - - - 
  for (i=0; i < NSN_DATA; i++ ) {

    CUTMASK  = INFO_DATA.TABLEVAR.CUTMASK[i] ;
    name     = INFO_DATA.TABLEVAR.name[i] ;    
    z        = INFO_DATA.TABLEVAR.zhd[i] ;    
    zerr     = INFO_DATA.TABLEVAR.zhderr[i] ;
    mu       = INFO_DATA.mu[i] - FITRESULT.SNMAG0; 
    mumodel  = INFO_DATA.mumodel[i] ;
    muerr    = INFO_DATA.muerr[i] ;
    muerrsq  = muerr * muerr ;

    if ( CUTMASK ) { continue ; }

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

  double muAvg, zM0, zAvg, cosPar[10];

  // load reference cosmology params
  cosPar[0] = INPUTS.parval[IPAR_OL] ;  // OL
  cosPar[1] = INPUTS.parval[IPAR_Ok] ;  // Ok
  cosPar[2] = INPUTS.parval[IPAR_w0] ;  // w0
  cosPar[3] = INPUTS.parval[IPAR_wa] ;  // wa

  if ( LDMP ) {
    printf(" xxx -------------------------------------- \n");
    printf(" xxx %s: ref OL,Ok,w0,wa = %.3f, %.3f, %.3f, %.3f \n",
	   fnam, cosPar[0], cosPar[1], cosPar[2], cosPar[3] );
    fflush(stdout);
  }

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
	     zAvg, ptr_zM0[iz], zM0, muAvg);
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
  // Each bias is a 6D function,
  //    mb-bias = function of (a=alpha,b=beta,g=gammadm z,x1,c)
  //    x1-bias = function of    "
  //    c-bias  = function of    "
  //
  // The data values z,x1,c are known a-prior and thus we
  // do the 3D interp here (before fitting) to get the bias 
  // as a function of the subset: z,x1,c.  Alpha,Beta,GammaDM
  // however, are SALT2mu fit parameters that are floated
  // in the fit,  and thus we store the bias in bins of 
  // alpha+beta+gammadm so that the fit-function (fcn) can 
  // interpolate for each fit-trial in the minimization.
  //
  // July 1 2016: also store muCOVscale[ia][ib]
  // Apr 18 2017: fix aweful index bug ia -> ib for beta

  BIASCORLIST_DEF BIASCORLIST ;
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM;
  int    NBINa, NBINb, NBINg, ia, ib, ig, ipar, istat_bias, idsample ;
  double z, m ;
  char   *name ;
  char   fnam[] = "storeDataBias" ;

  // ------------- BEGIN -------------

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;
  NBINa    = (*BININFO_SIM_ALPHA).nbin ;
  NBINb    = (*BININFO_SIM_BETA).nbin ;
  NBINg    = (*BININFO_SIM_GAMMADM).nbin ;
  name     = INFO_DATA.TABLEVAR.name[n];
  idsample = (int)INFO_DATA.TABLEVAR.IDSAMPLE[n];
  z        = (double)INFO_DATA.TABLEVAR.zhd[n];
  m        = (double)INFO_DATA.TABLEVAR.logmass[n];


  if ( DUMPFLAG ) {
    printf("\n");
    printf(" xxx ======================================== \n");
    printf(" xxx =========== %s DUMP for CID=%s ================ \n", 
	   fnam, name );
    printf(" xxx ======================================== \n");
    fflush(stdout);
  }

  if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { return(0); }

  BIASCORLIST.idsample = idsample;
  BIASCORLIST.z        = z ;
  BIASCORLIST.logmass  = m ;

  for(ipar=0; ipar<NLCPAR; ipar++ ) 
    { BIASCORLIST.FITPAR[ipar] = INFO_DATA.TABLEVAR.fitpar[ipar][n]; }

  for(ia = 0; ia < NBINa; ia++ ) {
    for(ib = 0; ib < NBINb; ib++ ) {
      for(ig = 0; ig < NBINg; ig++ ) {

	BIASCORLIST.alpha    = (*BININFO_SIM_ALPHA).avg[ia];
	BIASCORLIST.beta     = (*BININFO_SIM_BETA).avg[ib];
	BIASCORLIST.gammadm  = (*BININFO_SIM_GAMMADM).avg[ig];

	istat_bias = 
	  get_fitParBias(name, &BIASCORLIST, DUMPFLAG,             // in
			 &INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig]);//out

	if ( DUMPFLAG ) {
	  printf(" xxx %s: a=%.2f b=%.2f gDM=%5.2f (ia,ib,ig=%d,%d,%d) "
		 "istat_bias=%d \n",
		 fnam, BIASCORLIST.alpha,BIASCORLIST.beta,BIASCORLIST.gammadm,
		 ia,ib,ig, istat_bias);
	  fflush(stdout);
	}
	
	if ( istat_bias == 0 ) { return 0 ; }

	istat_bias = 
	  get_muCOVscale(name, &BIASCORLIST, DUMPFLAG,             // in
			 &INFO_DATA.MUCOVSCALE_ALPHABETA[n][ia][ib][ig]); //out

	if ( DUMPFLAG ) {
	  printf(" xxx %s: a=%.2f b=%.2f gDM=%.2f (ia,ib,ig=%d,%d,%d) "
		 "istat_muCOVscale=%d \n",
		 fnam, BIASCORLIST.alpha,BIASCORLIST.beta,BIASCORLIST.gammadm,
		 ia, ib, ig, istat_bias);
	  fflush(stdout);
      }
	
      if ( istat_bias == 0 ) { return 0 ; }

      } // end ig
    } // end ib
  }  // end ia

  return(1);

} // end storeDataBias



// ======================================================
int  storeBias_CCprior(int n) {

  // Created Jun 7 2019
  // for CCprior event 'n', store bias for ipar = mB,x1,c.
  // Returns 1 if bias can be determined;
  // returns 0 otherwise (to be rejected)
  //
  // This routine is analagous to storeDataBias for real data.
  //
  // Jan 17 2020: fix index bug setting BIASCORLIST.gammadm.

  int  NBINa   = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  NBINg   = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  char *name   = INFO_CCPRIOR.TABLEVAR.name[n];

  BIASCORLIST_DEF BIASCORLIST ;

  int    DUMPFLAG = ( n < -4 );
  int    ia, ib, ig, istat_bias ;
  char   fnam[] = "storeBias_CCprior" ;

  // ------------- BEGIN -------------

  // load local BIASCORLIST struct
  BIASCORLIST.z = 
    INFO_CCPRIOR.TABLEVAR.zhd[n];

  BIASCORLIST.logmass = 
    INFO_CCPRIOR.TABLEVAR.logmass[n];

  BIASCORLIST.FITPAR[INDEX_mB] = 
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_mB][n];

  BIASCORLIST.FITPAR[INDEX_x1] =
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_x1][n];

  BIASCORLIST.FITPAR[INDEX_c] =
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_c][n];

  BIASCORLIST.idsample =
    INFO_CCPRIOR.TABLEVAR.IDSAMPLE[n];
  
  if ( DUMPFLAG ) {
    printf(" xxx --------------------------------- \n");
    printf(" xxx %s: mB,x1,c = %.3f, %.3f, %.3f \n", fnam,
	   BIASCORLIST.FITPAR[INDEX_mB],
	   BIASCORLIST.FITPAR[INDEX_x1], 
	   BIASCORLIST.FITPAR[INDEX_c] );
    printf(" xxx %s: zhd=%.3f  idsample=%d \n",
	   fnam, BIASCORLIST.z, BIASCORLIST.idsample );
    fflush(stdout);
  }

  for(ia = 0; ia < NBINa; ia++ ) {
    for(ib = 0; ib < NBINb; ib++ ) {
      for(ig = 0; ig < NBINg; ig++ ) {
      
      BIASCORLIST.alpha    = INFO_BIASCOR.BININFO_SIM_ALPHA.avg[ia];
      BIASCORLIST.beta     = INFO_BIASCOR.BININFO_SIM_BETA.avg[ib];
      BIASCORLIST.gammadm  = INFO_BIASCOR.BININFO_SIM_GAMMADM.avg[ig];

      istat_bias =
	get_fitParBias(name, &BIASCORLIST, DUMPFLAG,
		       &INFO_CCPRIOR.FITPARBIAS_ALPHABETA[n][ia][ib][ig] );

      /* 
      if ( DUMPFLAG ) {
	printf(" xxx %s: a,b=%.2f,%.2f  BIAS(mB,x1,c)=%.3f,%.3f,%.3f \n",
	       fnam, BIASCORLIST.alpha, BIASCORLIST.beta,
	       FITPARBIAS_TMP.VAL[0], FITPARBIAS_TMP.VAL[1],
	       FITPARBIAS_TMP.VAL[2] );    fflush(stdout);
      }
      */

      if ( istat_bias == 0 ) { return 0 ; }

      } // end ig
    } // end ib
  }  // end ia
  
  return(1);

} // end storeBias_CCPrior


// ======================================================
int get_fitParBias(char *cid, 
		   BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		   FITPARBIAS_DEF  *FITPARBIAS) {

  // Created May 2016
  // For input BIASCORLIST = { mB,x1,c, z, alpha, beta, gammadm },
  // return FITPARBIAS[VAL,ERR,RMS] for  mB, x1, c.
  // CID is the SN name used only for error message.
  //
  // Function returns 1 if bias is found; returns 0 if no bias found.
  //
  // Apr 18 2017: enhance dump output.
  //
  // Nov 18 2019: check option to interpolate biasCor vs. logMass.
  //
  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z   = BIASCORLIST->z ;
  double m   = BIASCORLIST->logmass ;
  double mB  = BIASCORLIST->FITPAR[INDEX_mB];
  double x1  = BIASCORLIST->FITPAR[INDEX_x1];
  double c   = BIASCORLIST->FITPAR[INDEX_c];
  double a   = BIASCORLIST->alpha ;
  double b   = BIASCORLIST->beta ;
  double gDM = BIASCORLIST->gammadm ;
  int IDSAMPLE = BIASCORLIST->idsample ;
  int ID = IDSAMPLE;

  bool DO_BIASCOR_MU  = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;

  // get bin info into local variables
  int  NBINz   = CELLINFO_BIASCOR[ID].BININFO_z.nbin ;
  int  NBINm   = CELLINFO_BIASCOR[ID].BININFO_m.nbin ;
  int  NBINx1  = CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_x1].nbin ;
  int  NBINc   = CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_c].nbin ;

  double BINSIZE_z  = CELLINFO_BIASCOR[ID].BININFO_z.binSize ;
  double BINSIZE_m  = CELLINFO_BIASCOR[ID].BININFO_m.binSize ;
  double BINSIZE_x1 = CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_x1].binSize ;
  double BINSIZE_c  = CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_c].binSize ;
    
  int BADBIAS = 0 ;

  int J1D, IZ, IM, IX1, IC, OK, USEZ[MXz];
  int IZMIN, IZMAX, IMMIN, IMMAX, ICMIN, ICMAX, IX1MIN, IX1MAX ;
  int j1d, ia, ib, ig, iz, im, ix1, ic, ipar ;
  int NperCell, NSUM_Cell, NCELL_INTERP_TOT, NCELL_INTERP_USE ;

  double WGT, SUM_WGT, BINSIZE;
  double AVG_z, AVG_m, AVG_x1, AVG_c, avg_z, avg_m, avg_x1, avg_c ;
  double SUM_VAL[NLCPAR+1], SUM_ERR[NLCPAR+1];
  double SUM_SQERRINV[NLCPAR], SUM_SQRMS[NLCPAR+1] ;
  double VAL, ERR, RMS, dif, Dc, Dz, Dm, Dx1 ;

  double SUM_MBOFF[NLCPAR+1], SUM_MBSLOPE[NLCPAR+1];
  double DEBUG_LIST_DIF[4][50];
  int    DEBUG_LIST_INDX[4][50];
  double DEBUG_LIST_WGT[50];
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM ;

  char fnam[] = "get_fitParBias" ;

  // ------------- BEGIN ----------------

  int LDMP = DUMPFLAG ;  

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;

 START:

  // init output structure
  zero_FITPARBIAS(FITPARBIAS);

  // check "noBiasCor" option for this sample 
  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return(1); }

  // init output to crazy values.
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ )
    { FITPARBIAS->VAL[ipar]  = 666.0 ; } 

  // strip off local indices
  ia  = IBINFUN(a,   BININFO_SIM_ALPHA,    0, fnam );
  ib  = IBINFUN(b,   BININFO_SIM_BETA,     0, fnam );
  ig  = IBINFUN(gDM, BININFO_SIM_GAMMADM,  0, fnam );
  IZ  = IBINFUN(z,  &CELLINFO_BIASCOR[ID].BININFO_z,         0, fnam );
  IM  = IBINFUN(m,  &CELLINFO_BIASCOR[ID].BININFO_m,         0, fnam );
  IX1 = IBINFUN(x1, &CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_x1],0, fnam);
  IC  = IBINFUN(c,  &CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_c], 0, fnam);
  J1D = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][ig][IZ][IM][IX1][IC] ;


  if ( IZ  < 0 ) { return 0 ; }
  if ( IM  < 0 ) { return 0 ; }
  if ( IX1 < 0 ) { return 0 ; }
  if ( IC  < 0 ) { return 0 ; }

  // reset counters and weight for this a,b cell
  SUM_WGT = 0.0 ;
  NCELL_INTERP_TOT = NCELL_INTERP_USE  = NSUM_Cell = 0 ;
  
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {	    
    SUM_VAL[ipar]      = SUM_ERR[ipar]   = 0.0 ;
    SUM_SQERRINV[ipar] = SUM_SQRMS[ipar] = 0.0  ;
    SUM_MBOFF[ipar]    = SUM_MBSLOPE[ipar] = 0.0 ;
  }

  for(iz=0; iz < MXz; iz++ ) { USEZ[iz] = 0; }
  
  // ------------------------------------------
  // determine bins to interpolate such that current value
  // is between biasCor nodes (i.e., no extrapolation allowed)
  
  IZMIN=IZMAX=IZ;  ICMIN=ICMAX=IC ;    IX1MIN=IX1MAX=IX1 ;
  IMMIN=IMMAX=IM;

  AVG_z = CELLINFO_BIASCOR[IDSAMPLE].AVG_z[J1D] ;
  if ( z >= AVG_z ) { IZMAX++ ; } else { IZMIN--; }
  if (IZMIN<0){IZMIN=0;}  if(IZMAX>=NBINz){IZMAX = NBINz-1;}

  AVG_m = CELLINFO_BIASCOR[IDSAMPLE].AVG_m[J1D] ;
  if ( m >= AVG_m ) { IMMAX++ ; } else { IMMIN--; }
  if (IMMIN<0){IMMIN=0;}  if(IMMAX>=NBINm){IMMAX = NBINm-1;}
  if ( !INPUTS.interp_biascor_logmass ) { IMMIN = IMMAX = IM; }
  
  AVG_x1 = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][J1D] ;
  if ( x1 >= AVG_x1 ) { IX1MAX++; } else { IX1MIN--; }
  if (IX1MIN<0){IX1MIN=0;}  if(IX1MAX>=NBINx1){IX1MAX = NBINx1-1;}

  AVG_c = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c][J1D] ;
  if ( c >= AVG_c ) { ICMAX++ ;} else { ICMIN--; }
  if (ICMIN<0){ICMIN=0;}  if(ICMAX>=NBINc){ICMAX = NBINc-1;}


  if ( LDMP ) {  
    printf("\n") ;
    printf(" xxx ---------------------------------------------------- \n") ;
    if ( BADBIAS ) { printf("\t !!!!! BAD BIAS DETECTED !!!!! \n"); }
    printf(" xxx %s DUMP for CID=%s \n", fnam, cid );
    printf(" xxx  input: z=%.4f m=%.2f  mb,x1,c = %.2f, %.3f, %.3f \n",
	   z, m, mB, x1, c ); 

    printf(" xxx \t IZ,IM,IX1,IC=%d,%d,%d,%d   "
	   "NBIN(z,m,x1,c)=%d,%d,%d,%d \n",
	   IZ, IM, IX1, IC,     NBINz, NBINm, NBINx1, NBINc );

    printf(" xxx \t min/max = %d/%d(IZ) %d/%d(IM) %d/%d(IX1) %d/%d(IC) \n",
	   IZMIN,IZMAX, IMMIN,IMMAX, IX1MIN,IX1MAX, ICMIN,ICMAX);

    printf(" xxx  IDSAMPLE=%d -> %s(%s) \n", IDSAMPLE,
	   SAMPLE_BIASCOR[ID].NAME_SURVEYGROUP, 
	   SAMPLE_BIASCOR[ID].NAME_FIELDGROUP );
    fflush(stdout);
  }



  // determine max possible binsize in each dimension,
  // to use below for interpolation. 
  
  for(iz = IZMIN; iz <= IZMAX; iz++ ) {    
    for(im = IMMIN; im <= IMMAX; im++ ) {    
      for(ix1 = IX1MIN; ix1 <= IX1MAX; ix1++ ) {
	for(ic = ICMIN; ic <= ICMAX; ic++ ) {
	
	  if( iz == IZ && im==IM && ix1==IX1 && ic==IC ) { continue ; }
	
	  j1d = CELLINFO_BIASCOR[ID].MAPCELL[ia][ib][ig][iz][im][ix1][ic];
	  NperCell = CELLINFO_BIASCOR[ID].NperCell[j1d] ;
	  if ( NperCell < BIASCOR_MIN_PER_CELL  ) { continue; }

	  avg_z  = CELLINFO_BIASCOR[IDSAMPLE].AVG_z[j1d] ;
	  avg_m  = CELLINFO_BIASCOR[IDSAMPLE].AVG_m[j1d] ;
	  avg_x1 = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][j1d] ;
	  avg_c  = CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c][j1d] ;

	  if ( avg_z  > 9000.0 ) { continue; }
	  if ( avg_m  > 9000.0 ) { continue; }
	  if ( avg_x1 > 9000.0 ) { continue; }
	  if ( avg_c  > 9000.0 ) { continue; }
	  
	  BINSIZE = fabs(AVG_z - avg_z);
	  if( BINSIZE > BINSIZE_z && iz!=IZ) { BINSIZE_z = BINSIZE; }

	  BINSIZE = fabs(AVG_m - avg_m);
	  if( BINSIZE > BINSIZE_m && im!=IM) { BINSIZE_m = BINSIZE; }
	  
	  BINSIZE = fabs(AVG_x1 - avg_x1);
	  if( BINSIZE > BINSIZE_x1 && ix1!=IX1) { BINSIZE_x1 = BINSIZE; }
	  
	  BINSIZE = fabs(AVG_c - avg_c);
	  if( BINSIZE > BINSIZE_c && ic!=IC) { BINSIZE_c = BINSIZE; }	
	}      
      }
    }
  }


  if ( LDMP ) {
    printf(" xxx  AVGCELL(z,m,x1,c) = %7.3f, %7.3f, %7.3f, %7.3f \n",
	   AVG_z, AVG_m, AVG_x1, AVG_c);
    printf(" xxx  BINSIZE(z,m,x1,c) = %7.3f, %7.3f, %7.3f, %7.3f \n",
	   BINSIZE_z, BINSIZE_m, BINSIZE_x1, BINSIZE_c ); fflush(stdout);

    printf(" xxx \n");
    printf(" xxx iz,im,ix,ic    <z>   <m>   <x1>   <c>      "
	   "dmB    dx1    dc  Ncell WGT\n");  
    fflush(stdout);
  }

  // return if any binSize is not defined
  if ( BINSIZE_z  == 0.0 || BINSIZE_z  > 9000. ) { return 0; }
  if ( BINSIZE_m  == 0.0 || BINSIZE_m  > 9000. ) { return 0; }
  if ( BINSIZE_x1 == 0.0 || BINSIZE_x1 > 9000. ) { return 0; }
  if ( BINSIZE_c  == 0.0 || BINSIZE_c  > 9000. ) { return 0; }

  
  // ----------------------------------------------
  // -------- start 4D loop over cells -------------
  
  for(iz = IZMIN; iz <= IZMAX; iz++ ) {    
    for(im = IMMIN; im <= IMMAX; im++ ) {    
      for(ix1 = IX1MIN; ix1 <= IX1MAX; ix1++ ) {
	for(ic = ICMIN; ic <= ICMAX; ic++ ) {
	  
	  NCELL_INTERP_TOT++ ;
	
	  j1d = CELLINFO_BIASCOR[ID].MAPCELL[ia][ib][ig][iz][im][ix1][ic] ;
	  
	  if ( j1d >= CELLINFO_BIASCOR[ID].NCELL || j1d < 0 ) {
	    sprintf(c1err,"Invalid j1d=%d ", j1d);
	    sprintf(c2err,"ia=%d ib=%d ig=%d iz=%d im=%d ix1=%d ic=%d (cid=%s)",
		    ia, ib, ig, iz, im, ix1, ic, cid);
	    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	  }
	
	  // require something in a cell to use the bias estimate
	  NperCell = CELLINFO_BIASCOR[ID].NperCell[j1d] ;
	  if ( NperCell < BIASCOR_MIN_PER_CELL  ) { continue; }
	  
	  // get distance between current data value and wgted-avg in bin
	  
	  dif = z - CELLINFO_BIASCOR[ID].AVG_z[j1d] ;
	  Dz  = fabs(dif/BINSIZE_z) ;

	  if ( INPUTS.interp_biascor_logmass ) {
	    dif = m - CELLINFO_BIASCOR[ID].AVG_m[j1d] ;
	    Dm  = fabs(dif/BINSIZE_m) ;
	  }
	  else
	    { Dm = 0.0 ; }
	  
	  dif = x1 - CELLINFO_BIASCOR[ID].AVG_LCFIT[INDEX_x1][j1d] ;
	  Dx1 = fabs(dif/BINSIZE_x1);
	  
	  dif = c - CELLINFO_BIASCOR[ID].AVG_LCFIT[INDEX_c][j1d] ;
	  Dc  = fabs(dif/BINSIZE_c);
	  
	  WGT = (1.0 - Dz) * (1.0 - Dx1) * (1.0 - Dc) * ( 1.0 - Dm); 
	  
	  // prepare DEBUG_LISTs in case of abort
	  DEBUG_LIST_DIF[0][NCELL_INTERP_USE]  = Dz;
	  DEBUG_LIST_DIF[1][NCELL_INTERP_USE]  = Dm;
	  DEBUG_LIST_DIF[2][NCELL_INTERP_USE]  = Dx1;
	  DEBUG_LIST_DIF[3][NCELL_INTERP_USE]  = Dc; 
	  DEBUG_LIST_INDX[0][NCELL_INTERP_USE] = iz ;
	  DEBUG_LIST_INDX[1][NCELL_INTERP_USE] = im ;
	  DEBUG_LIST_INDX[2][NCELL_INTERP_USE] = ix1 ;
	  DEBUG_LIST_INDX[3][NCELL_INTERP_USE] = ic ;
	  DEBUG_LIST_WGT[NCELL_INTERP_USE]     = WGT ;
	  
	  if ( WGT < 0.0 ) { WGT = 0.0 ; }
	
	  SUM_WGT   += WGT ;
	  NSUM_Cell += NperCell ; // sum Nevt in nbr cells
	  NCELL_INTERP_USE++ ;
	  USEZ[iz]++ ;
	
	  // now loop over the data fit params (mB,x1,c) to update bias sums
	  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {	    
	    VAL  = INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].VAL[ipar] ; 
	    ERR  = INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].ERR[ipar] ; 
	    RMS  = INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].RMS[ipar] ; 
	    
	    SUM_VAL[ipar]      += (WGT * VAL );
	    SUM_SQRMS[ipar]    += (WGT * RMS*RMS );
	    SUM_SQERRINV[ipar] += (WGT / (ERR*ERR) );
	    SUM_ERR[ipar]      += (WGT * ERR ) ;	 
	  } // end ipar
	  
	  if ( LDMP ) {
	    printf( " xxx %2d,%2d,%2d,%2d  %5.4f %5.2f %5.2f %6.3f  ",
		    iz, im, ix1, ic,  
		    CELLINFO_BIASCOR[IDSAMPLE].AVG_z[j1d],
		    CELLINFO_BIASCOR[IDSAMPLE].AVG_m[j1d],
		    CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_x1][j1d],
		    CELLINFO_BIASCOR[IDSAMPLE].AVG_LCFIT[INDEX_c ][j1d] ) ;
	  
	    printf("  %6.3f %5.2f %6.3f  %3d %.2f ",
		   INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_mB],
		   INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_x1],
		   INFO_BIASCOR.FITPARBIAS[IDSAMPLE][j1d].VAL[INDEX_c],
		   NperCell, WGT );	  
	    printf("\n");  fflush(stdout);
	  }
	  
	} // end ic
      } // end ix1
    }  // end im
  } // end iz
  
  
  if ( LDMP ) {
    printf("\n") ;
    printf(" xxx\t SUM_WGT = %le   NCELL_INTERP_USE=%d of %d (CUT=3)\n", 
	   SUM_WGT, NCELL_INTERP_USE, NCELL_INTERP_TOT );
    fflush(stdout);
  }

   // require enough cells for interpolation (July 2016)
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
    print_preAbort_banner(fnam);
    printf("  IZMIN/MAX=%d/%d   IMMIN/MAX=%d,%d  "
	   "IX1MIN/MAX=%d/%d   ICMIN/MAX=%d/%d\n",	   
	   IZMIN,IZMAX,  IMMIN, IMMAX, IX1MIN,IX1MAX,   ICMIN, ICMAX);
    printf("  BINSIZE(z,m,x1,c) = %.3f, %.3f, %.3f, %.3f \n",
	   BINSIZE_z, BINSIZE_m, BINSIZE_x1, BINSIZE_c );
    
    printf("\n");
    printf("   cell   iz im ix1 ic    Dz     Dm     Dx1     Dc     WGT \n");
    
    for(icell=0; icell < NCELL_INTERP_USE; icell++ ) {
      printf("    %3d   %2d %2d %2d %2d   %6.3f %6.3f %6.3f %6.3f   %.3f\n"
	     ,icell
	     ,DEBUG_LIST_INDX[0][icell]
	     ,DEBUG_LIST_INDX[1][icell]
	     ,DEBUG_LIST_INDX[2][icell]
	     ,DEBUG_LIST_INDX[3][icell]
	     ,DEBUG_LIST_DIF[0][icell]
	     ,DEBUG_LIST_DIF[1][icell]
	     ,DEBUG_LIST_DIF[2][icell]
	     ,DEBUG_LIST_DIF[3][icell]
	     ,DEBUG_LIST_WGT[icell] );
    }
    sprintf(c1err,"SUM_WGT=%f for  CID=%s", SUM_WGT, cid) ;
    sprintf(c2err,"a=%.3f b=%.2f gDM=%.3f  "
	    "z=%.4f  m=%.2f mB=%.4f x1=%.4f c=%.4f", 
	    a, b, gDM, z, m, mB, x1, c ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }      
  
  // update bias for each fit param (ipar=mB,x1,c)
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {
    FITPARBIAS->VAL[ipar]  = SUM_VAL[ipar] / SUM_WGT ;
    FITPARBIAS->ERR[ipar]  = sqrt((1.0/SUM_SQERRINV[ipar]) / SUM_WGT);
    FITPARBIAS->RMS[ipar]  = sqrt(SUM_SQRMS[ipar]/SUM_WGT);
    if ( fabs( FITPARBIAS->VAL[ipar]) > 20.0 ) { BADBIAS=1; }
  } 
  
  if ( LDMP ) {
    printf("\n") ;
    for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {
      printf(" xxx bias(%2.2s) = %7.4f +-  %7.4f   SQERRINV=%le\n"
	     ,BIASCOR_NAME_LCFIT[ipar]
	     ,FITPARBIAS->VAL[ipar], FITPARBIAS->ERR[ipar] 
	     ,SUM_SQERRINV[ipar]  );
    }

    if ( DO_BIASCOR_MU == false ) {
      double bias_mB = FITPARBIAS->VAL[INDEX_mB];
      double bias_x1 = FITPARBIAS->VAL[INDEX_x1];
      double bias_c  = FITPARBIAS->VAL[INDEX_c] ;
      double muBias  = bias_mB + a*bias_x1 - b*bias_c ;
      if ( LDMP ) 
	{ printf(" xxx bias(mu) = %7.4f \n", muBias ); }
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

  int ia, ib, ig, iz, ic, IZ, IC, j1d ;

  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z  = BIASCORLIST->z ;
  double c  = BIASCORLIST->FITPAR[INDEX_c];
  double a  = BIASCORLIST->alpha ;
  double b  = BIASCORLIST->beta ;
  double g  = BIASCORLIST->gammadm ;
  int IDSAMPLE = BIASCORLIST->idsample ;  

  // get bin info into local variables
  int  NBINz   = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.nbin ;
  int  NBINc   = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].nbin ;
  
  double BINSIZE_z  
    = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_z.binSize ;
  double BINSIZE_c  
    = CELLINFO_MUCOVSCALE[IDSAMPLE].BININFO_LCFIT[INDEX_c].binSize ;
  
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM ;
  float *ptr_MUCOVSCALE ;

  double dif, Dz, Dc, WGT, SUM_WGT, SUM_muCOVscale, muCOVscale_biascor ;
  double muCOVscale_local  ;
  char fnam[] = "get_muCOVscale" ;

  // -------------- BEGIN --------------

  BININFO_SIM_ALPHA    = &INFO_BIASCOR.BININFO_SIM_ALPHA;
  BININFO_SIM_BETA     = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM  = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;
  ptr_MUCOVSCALE       = INFO_BIASCOR.MUCOVSCALE[IDSAMPLE] ;

  // init output
  muCOVscale_local = 1.0 ;
  *muCOVscale = muCOVscale_local ;

  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOV)==0   ) { return(1); }

  // check "noBiasCor" option (Aug 20 2016)
  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return(1); }

  // get indices
  ia    = IBINFUN(a,  BININFO_SIM_ALPHA,   0, "" );
  ib    = IBINFUN(b,  BININFO_SIM_BETA,    0, "" );
  ig    = IBINFUN(g,  BININFO_SIM_GAMMADM, 0, "" );
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

      j1d = CELLINFO_MUCOVSCALE[IDSAMPLE].MAPCELL[ia][ib][ig][iz][0][0][ic] ;

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

      muCOVscale_biascor = (double)ptr_MUCOVSCALE[j1d];
      SUM_muCOVscale += ( WGT * muCOVscale_biascor );
    } // ic
   } // iz
  
  
  if( SUM_WGT > 0.01 ) {
    muCOVscale_local = SUM_muCOVscale / SUM_WGT ;
  }

  if ( DUMPFLAG) {
    printf(" xxx %s: muCOVscale = %f/%f = %f \n", 
	   fnam, SUM_muCOVscale, SUM_WGT, muCOVscale_local);
    fflush(stdout);
  }

  // load output arg.
  *muCOVscale = muCOVscale_local ;

  return(1);

} // end get_muCOVscale


// ======================================================
void setup_CELLINFO_biasCor(int IDSAMPLE) {

  // Created Aug 23 2016
  // Setup 5D biasCor cells for input IDSAMPLE

  int NSAMPLE   = NSAMPLE_BIASCOR ;
  int MEMCELL   = NSAMPLE * sizeof(CELLINFO_DEF);
  int INDX;
  char fnam[] = "setup_CELLINFO_biasCor" ;

  // ------------ BEGIN -----------

  // on first call, allocate CELLINFO structures
  if ( IDSAMPLE == 0 ) {

    printf("\n# ============== %s ================= \n", fnam);

    CELLINFO_BIASCOR    = (CELLINFO_DEF*) malloc ( MEMCELL );
    CELLINFO_MUCOVSCALE = (CELLINFO_DEF*) malloc ( MEMCELL );

    // setup bining for SIMalpha,beta; note storage in different structure
    // since alpha,beta binning is fixed for all IDSAMPLEs
    printf("\n\t\t Global Alpha,Beta Binning \n"); fflush(stdout);

    setup_BININFO_biasCor(IDSAMPLE, 100, MXa,  
			  &INFO_BIASCOR.BININFO_SIM_ALPHA );     
    setup_BININFO_biasCor(IDSAMPLE, 200, MXb,
			  &INFO_BIASCOR.BININFO_SIM_BETA );   
    setup_BININFO_biasCor(IDSAMPLE, 300, MXg,			    
			  &INFO_BIASCOR.BININFO_SIM_GAMMADM ); 
  }


  CELLINFO_BIASCOR[IDSAMPLE].BININFO_z.nbin  = 0 ; 
  CELLINFO_BIASCOR[IDSAMPLE].NCELL           = 0 ;

  if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR==0 ) { return ; }

  printf("\n\t Setup biasCor bins for SAMPLE = %s \n",
	 SAMPLE_BIASCOR[IDSAMPLE].NAME ); fflush(stdout);

  // setup binning for redshift
  setup_BININFO_biasCor(IDSAMPLE, -1, MXz, 
			&CELLINFO_BIASCOR[IDSAMPLE].BININFO_z ); 

  // setup binning for logmass (Aug 2019)
  setup_BININFO_biasCor(IDSAMPLE, -2, MXm, 
			&CELLINFO_BIASCOR[IDSAMPLE].BININFO_m ); 

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
  // ipar_LCFIT = -2    --> logmass
  // ipar_LCFIT = 0,1,2 --> mB, x1, c    (fitted)
  // ipar_LCFIT = 100   --> SIM_ALPHA    (generated)
  // ipar_LCFIT = 200   --> SIM_BETA     (generated)
  // ipar_LCFIT = 300   --> SIM_GAMMADM  (generated)
  // Output: *BININFO
  // 
  // Oct 8 2016: adjust LO_MAX calc to allow for TEXT truncation
  //
  // Apr 30 2017: change "nbin>=MAXBIN" -> "nbin>MAXBIN" for abort trap.
  // Jul 16 2019: process SIM_gammaDM
  // Mar 31 2020: allow 1 bin for beta as well as gamma; see OK1BIN

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
  else if ( ipar_LCFIT == -1 ) {  
    // redshift
    VAL_MIN = SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[0];
    VAL_MAX = SAMPLE_BIASCOR[IDSAMPLE].RANGE_REDSHIFT[1];
    VAL_BIN = SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_REDSHIFT ;
    sprintf(NAME,"z");

    // add one more z bin on the high-z side to get better
    // interpolation in last z-bin
    VAL_MAX += VAL_BIN ;    
  }
  else if ( ipar_LCFIT == -2 ) {  
    // logmass
    VAL_MIN = SAMPLE_BIASCOR[IDSAMPLE].RANGE_LOGMASS[0];
    VAL_MAX = SAMPLE_BIASCOR[IDSAMPLE].RANGE_LOGMASS[1];
    VAL_BIN = SAMPLE_BIASCOR[IDSAMPLE].BINSIZE_LOGMASS ;
    sprintf(NAME,"m");
  }
  else if ( ipar_LCFIT == 100*INDEX_x1 ) {
    // get info for SIMalpha binning
    get_BININFO_biasCor_alphabeta("SIM_alpha", 
				  &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIM_alpha");
  }
  else if ( ipar_LCFIT == 100*INDEX_c ) {
    // get info for SIM_beta binning
    get_BININFO_biasCor_alphabeta("SIM_beta", 
				  &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIM_beta");
  }
  else if ( ipar_LCFIT == 300 ) {
    // get info for SIM_gammaDM binning (Jul 2019)
    get_BININFO_biasCor_alphabeta("SIM_gammaDM", 
				  &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIM_gammaDM");
    INFO_BIASCOR.GAMMADM_OFFSET = (VAL_MIN + VAL_MAX)/2.0;

    /*
    printf(" xxx %s: VAL_MIN/MAX=%.3f/%.3f VAL_BIN=%.3f  "
	   "OFFSET = %.3f \n", 
	   fnam, VAL_MIN, VAL_MAX, VAL_BIN, INFO_BIASCOR.GAMMADM_OFFSET);
    */
  }
  else {
    sprintf(c1err,"Invalid ipar_LCFIT=%d", ipar_LCFIT );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" );    
  }


  // allow 1 bin for beta and gamma; but not for other params
  bool ONEBIN = (VAL_MAX == VAL_MIN) ;
  bool OK1BIN = ( ONEBIN && (ipar_LCFIT >= 100) ) ;
  if ( VAL_MAX <= VAL_MIN  && !OK1BIN ) {
    sprintf(c1err,"%s VAL_MAX=%f  <=  VAL_MIN=%f", 
	    NAME, VAL_MAX, VAL_MIN);
    sprintf(c2err,"VAL_BIN = %f \n", VAL_BIN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
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

  BININFO->nbin    = nbin ;
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
  // Jul 16 2019: allow for gammaDM
  // Mar 31 2020: store min/max per IDSURVEY (since IDSAMPLE not yet known)

  double valmin_loc, valmax_loc, valbin_loc;
  double val, val_last, val1st, val2nd ;
  double valmin_sample[MXIDSURVEY], valmax_sample[MXIDSURVEY];
  float  *ptrVal_f = NULL;
  short int *ptrVal_index, *ptr_IDSURVEY ;
  bool   IS_GRID, IS_GAMMADM = false ;
  int    irow, unsort, NVAL, NBMAX ;
  int    NROW, IDSURVEY, NONIA_INDEX ;
  char fnam[] = "get_BININFO_biasCor_alphabeta" ;

  // ------------- BEGIN --------------

  *VAL_MIN = *VAL_MAX = *VAL_BIN = 0.0 ;

  valmin_loc = +1.0E8 ;
  valmax_loc = -1.0E8 ;
  valbin_loc =  0.0 ;

  NROW         = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptrVal_index = INFO_BIASCOR.TABLEVAR.SIM_NONIA_INDEX;
  ptr_IDSURVEY = INFO_BIASCOR.TABLEVAR.IDSURVEY ;   
  if ( strcmp(varName,"SIM_alpha") == 0 ) 
    { ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_ALPHA ; NBMAX=MXa; }
  else if ( strcmp(varName,"SIM_beta") == 0 )  
    { ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_BETA ;  NBMAX=MXb; }
  else if ( strcmp(varName,"SIM_gammaDM") == 0 ) { 
    ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_GAMMADM ; NBMAX=MXg; 
    IS_GAMMADM = true ;
  }
  else {
    sprintf(c1err,"Invalid varName = '%s' ", varName);
    errmsg(SEV_FATAL, 0, fnam, c1err, "Must be SIMalpha or SIMbeta");  
  }
  

  for(IDSURVEY=0; IDSURVEY < MXIDSURVEY; IDSURVEY++ ) {
    valmin_sample[IDSURVEY] = +1000.0 + IDSURVEY ;
    valmax_sample[IDSURVEY] = -1000.0 - IDSURVEY ;
  }


  int  ORDER_SORT   = +1 ; // increasing order
  int *INDEX_UNSORT = (int*)malloc( NROW * sizeof(int) );
  sortFloat( NROW, ptrVal_f, ORDER_SORT, INDEX_UNSORT ) ;

  val1st = val2nd = -99999. ;
  NVAL = 0;  val_last = -99999. ;
  for ( irow=0; irow < NROW; irow++ ) {
    unsort      = INDEX_UNSORT[irow];
    val         = (double)ptrVal_f[unsort];
    NONIA_INDEX = ptrVal_index[unsort];
    IDSURVEY    = (int)ptr_IDSURVEY[unsort];

    if ( NONIA_INDEX != 0 ) { continue ; }

    if ( val > val_last ) { 
      NVAL++ ; 
      // store 1st and send value to get binsize in case of grid
      if ( NVAL == 1 ) { val1st = val2nd = val ; }
      if ( NVAL == 2 ) { val2nd = val ; }
    }
    val_last = val ;

    if ( val < valmin_sample[IDSURVEY] ) { valmin_sample[IDSURVEY]=val; }
    if ( val > valmax_sample[IDSURVEY] ) { valmax_sample[IDSURVEY]=val; }
  }

  // Feb 2020: if lots of gammaDM bins, allow physical distribution
  //  and do NOT add this dimenstion to BBC biasCor
  if ( IS_GAMMADM && NVAL > 5 ) { NVAL = 1; val2nd = val_last; }

  // add error check on number of a,b bins (Apr 28 2017)
  if ( NVAL > NBMAX ) {
    sprintf(c1err,"%d %s values exceeds bound of %d",
	    NVAL, varName, NBMAX);
    sprintf(c2err,"Check biasCor sample.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
  }


  // ---------------------------
  if ( NVAL < 10 ) { 
    IS_GRID = true ; 
    valbin_loc = val2nd - val1st ;
    valmin_loc = val1st   - 0.5*valbin_loc ;
    valmax_loc = val_last + 0.5*valbin_loc ;
    
    // check that GRID min & max is the same for each IDSAMPLE
    check_abg_minmax_biasCor(varName,valmin_sample,valmax_sample);

  }
  else {
    IS_GRID = false ;  // continuous
    unsort = INDEX_UNSORT[0];        valmin_loc = (double)ptrVal_f[unsort];
    unsort = INDEX_UNSORT[NROW-1];   valmax_loc = (double)ptrVal_f[unsort];

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

  free(INDEX_UNSORT); // free memory
  
  // --------------------------
  // load function args
  *VAL_MIN = valmin_loc ;
  *VAL_MAX = valmax_loc ;
  *VAL_BIN = valbin_loc ;

  return ;

} // end get_BININFO_biasCor_alphabeta


// ================================================
void check_abg_minmax_biasCor(char *varName, double *valmin_list,
			      double *valmax_list) {

  // Created Mar 31 2020
  // If min or max is different for any IDSAMPLE, abort.
  // This enforces requirement that alpha, beta, and gamma
  // grid must be the same for each biasCor sub-sample.
  //

  int idsurvey ;
  int NERR=0;
  double valmin, valmax, valmin_ref, valmax_ref;
  char fnam[] = "check_abg_minmax_biasCor" ;
  // ----------- BEGIN ------------

  valmin_ref = -9.0 ;
  valmax_ref = -9.0 ;

  for(idsurvey=0; idsurvey < MXIDSURVEY; idsurvey++ ) {
    valmin = valmin_list[idsurvey];
    valmax = valmax_list[idsurvey];

    if ( valmax < -900.0 ) { continue; }
    if ( valmin_ref == -9.0 ) { valmin_ref = valmin; }
    if ( valmax_ref == -9.0 ) { valmax_ref = valmax; }

    if ( valmin != valmin_ref ) { NERR++ ; }
    if ( valmax != valmax_ref ) { NERR++ ; } 

  }

  if ( NERR > 0 ) {
    print_preAbort_banner(fnam);        
    printf("   IDSURVEY   SURVEY    min(%s)  max(%s) \n", varName, varName);   
    for(idsurvey=0; idsurvey < MXIDSURVEY; idsurvey++ ) {
      char *SURVEY = SURVEY_INFO.SURVEYDEF_LIST[idsurvey] ;
      valmin = valmin_list[idsurvey];
      valmax = valmax_list[idsurvey];
      if ( valmax < -900.0 ) { continue; }
      printf("     %2d %12s     %8.4f       %8.4f \n",  
	     idsurvey, SURVEY, valmin, valmax);
      fflush(stdout);
    }

    sprintf(c1err,"%d min/max errors in alpha,beta,gamma grid",	NERR);
    sprintf(c2err,"Check min & max vs. IDSAMPLE above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
  }

  return ;
} // end check_abg_minmax_biasCor

// ================================================
void get_muBias(char *NAME,
		BIASCORLIST_DEF *BIASCORLIST, 
		FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MXb][MXg],
		double         (*MUCOVSCALE_ABGRID)[MXb][MXg],
		INTERPWGT_AlphaBetaGammaDM *INTERPWGT,  
		double *fitParBias,
		double *muBias, double *muBiasErr, double *muCOVscale ) { 

  
  // Created Jan 2016
  // Inputs:
  //   NAME        = name of SN (for err msg only)
  //   BIASCORLIST = z,logmass,a,b,g,mB,x1,c 
  //   FITPARBIAS  = pre-computed bias for {mB,x1,c} on grid of {a,b}
  //   INTERPWGT   = wgt at each ia,ib
  //
  // Ouptuts:
  // *fitParBias  = bias on mB,x1,c, interpolated over alpha & beta
  //  muBias      = bias on distance
  //  muBiasErr   = error on above (based on biasCor sim stats)
  //  muCOVscale  = scale bias to apply to muErr 
  //

  double alpha    = BIASCORLIST->alpha ;
  double beta     = BIASCORLIST->beta  ;
  double gammadm  = BIASCORLIST->gammadm  ;
  //  double logmass  = BIASCORLIST->logmass  ;
  double z        = BIASCORLIST->z  ;
  double mB       = BIASCORLIST->FITPAR[INDEX_mB] ;
  double x1       = BIASCORLIST->FITPAR[INDEX_x1] ;
  double c        = BIASCORLIST->FITPAR[INDEX_c] ;

  double muBias_local     = 0.0 ;
  double muBiasErr_local  = 0.0 ;
  double muCOVscale_local = 0.0 ;

  double VAL, ERR, SQERR ;
  double WGTabg, WGTpar, WGTpar_SUM[NLCPAR+1];
  double biasVal[NLCPAR+1], biasErr[NLCPAR+1] ;
  double MUCOEF[NLCPAR+1];

  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;

  int  ia, ib, ig, ipar ;
  char fnam[] = "get_muBias";

  // --------------- BEGIN ------------

  for(ipar=0; ipar < NLCPAR+1; ipar++ )  { 
    biasVal[ipar] = biasErr[ipar] = 0.0 ;
    WGTpar_SUM[ipar] = 0.0 ;
  }

  MUCOEF[INDEX_mB] = +1.0 ;
  MUCOEF[INDEX_x1] = +alpha ;
  MUCOEF[INDEX_c ] = -beta ;
  MUCOEF[INDEX_mu] = +1.0 ;
 

  // interpolate bias on grid of alpha,beta,gammadm
  // The WGT at each corner was computed earlier.
  for(ia=INTERPWGT->ia_min; ia <= INTERPWGT->ia_max; ia++ ) {
    for(ib=INTERPWGT->ib_min; ib <= INTERPWGT->ib_max; ib++ ) {
      for(ig=INTERPWGT->ig_min; ig <= INTERPWGT->ig_max; ig++ ) {
      
	WGTabg = INTERPWGT->WGT[ia][ib][ig] ; // weight for this a,b,g
	if ( WGTabg < -1.0E-12 || WGTabg > 1.0000001 ) {
	  sprintf(c1err,"Undefined WGTabg=%le for CID=%s", WGTabg, NAME);
	  sprintf(c2err,"ia,ib,ig=%d,%d,%d  a,b,g=%.3f,%.2f,%.3f  z=%.3f",
		  ia, ib, ig, alpha, beta, gammadm, z );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
	}
      
	for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX ; ipar++ ) {
	  VAL = FITPARBIAS_ABGRID[ia][ib][ig].VAL[ipar];
	  ERR = FITPARBIAS_ABGRID[ia][ib][ig].ERR[ipar];
	  if ( VAL > 600.0 || isnan(ERR) ) {
	    sprintf(c1err,"Undefined %s-bias = %.3f +- %.3f for CID=%s", 
		    BIASCOR_NAME_LCFIT[ipar], VAL, ERR, NAME);
	    sprintf(c2err,"ia,ib,ig=%d,%d,%d  a,b,g=%.3f,%.2f,%.3f",
		    ia, ib, ig, alpha, beta, gammadm );
	    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
	  }
	  
	  WGTpar = WGTabg ;
	  WGTpar_SUM[ipar] += WGTpar;
	  biasVal[ipar] += ( WGTpar * VAL );
	  biasErr[ipar] += ( WGTpar * ERR ); 
	} // end ipar
	
	// now the COV scale (Jul 1 2016)
	VAL = MUCOVSCALE_ABGRID[ia][ib][ig] ;
	muCOVscale_local += ( WGTabg * VAL ) ;
	
      } // end ig
    } // end ib
  } // end ia


  // Note that sum(WGT) is already normalized to 1.000,
  // so no need to divide by SUMWGT here.

  muBias_local = SQERR = 0.0 ;
    
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX ; ipar++ ) {
    biasVal[ipar] /= WGTpar_SUM[ipar] ;
    biasErr[ipar] /= WGTpar_SUM[ipar] ;

    fitParBias[ipar]   = biasVal[ipar] ;     //store output array
    muBias_local += ( MUCOEF[ipar] * biasVal[ipar] ) ;
    ERR           = ( MUCOEF[ipar] * biasErr[ipar] ) ;
    SQERR        += (ERR * ERR);
  } 

  // check for NaN on SQERR
  if ( isnan(SQERR) ) {
    print_preAbort_banner(fnam);
    for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX ; ipar++ ) {
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

    print_preAbort_banner(fnam);
    printf(" xxx alpha=%f  beta=%f \n", alpha, beta);

    for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX ; ipar++ ) {
      printf(" xxx bias(%2s) = %7.3f +- %.3f  \n"
	     ,BIASCOR_NAME_LCFIT[ipar],biasVal[ipar], biasErr[ipar]  );
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
double get_gammadm_host(double z, double logmass, double *hostPar) {

  // return SN mag offset based on host mass;
  // Brighter in massive hosts.
  //ur
  // Inputs:
  //  z        = redshift
  //  logmass  = logmass, or whatever host property is used
  //  hostPar   [see below]
  //
  // strip off host properties into local variables
  // These are floated (or fixed) in the fit.
  //
  // Jan 18 2020:
  //  Holy crapola, just noticed bug that returned mag is
  //  fainter on bright hosts, not brighter. This bug was 
  //  canclled in fcn (chi2 function) that had wrong sign
  //  for gamma in Tripp equation.

  double gamma0      = hostPar[0]; // magSPlit at z=0
  double gamma1      = hostPar[1]; // dgamma/dz
  double logmass_cen = hostPar[2]; // Prob=0.5 at this logmass
  double logmass_tau = hostPar[3]; // tau param in transition function

  double gamma, FermiFun, arg ;
  double magoff = 0.0 ;

  //  char   fnam[] = "get_gammadm_host" ;

  // ------------ BEGIN ---------------

  if ( !INPUTS.USE_GAMMA0 ) { return(magoff); }

  gamma      = gamma0 + z*gamma1 ;
  arg        = -( logmass - logmass_cen ) / logmass_tau ;
  FermiFun   = 1.0/(1.0 + exp(arg)) ; 
  // xxx REMOVE BUG  magoff     = gamma * (  FermiFun - 0.5 ) ;
  magoff     = gamma * ( 0.5 - FermiFun ) ;
  magoff -= INFO_BIASCOR.GAMMADM_OFFSET ;

  /*
  printf(" xxx xval=%.2f cen=%.1f tau=%.3f Fun=%.3f magoff=%.3f \n",
	 XVAL_GAMMA, logmass_cen, logmass_tau, FermiFun, magoff);
  debugexit(fnam); 
  */

  return(magoff);

} // end get_gammadm_host


// ==================================================
void prepare_CCprior(void) {

  // Created Jan 2016
  // If simFile is not null, then prepare CC prior for BEAMS-like
  // usage in fit likelihood.
  //
  // Jun 20 2018: abort on 1D biasCor.
  
  int  EVENT_TYPE   = EVENT_TYPE_CCPRIOR ;
  int  NSAMPLE      = NSAMPLE_BIASCOR ;
  int  NDIM_BIASCOR = INFO_BIASCOR.NDIM ;
  int idsample, USE_CCPRIOR_H11 ;
  char fnam[] = "prepare_CCprior" ;

  // ------------- BEGIN -------------

  USE_CCPRIOR_H11 = INFO_CCPRIOR.USEH11;

  if ( INPUTS.nfile_CCprior == 0 ) { return; }
  
  if ( NDIM_BIASCOR == 1 ) {
    sprintf(c1err,"Cannot use 1D biasCor with CC likelihood.");
    sprintf(c2err,"Either remove CC term, or use 5D biasCor.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  // setup MURES bins and z bins for storing binned contamination info. 
  setup_contam_CCprior();

  if ( USE_CCPRIOR_H11 ) { 
    sprintf(BANNER,"%s: use CC mu-vs-z prior from Hlozek 2011", fnam);
    print_banner(BANNER);    
    return ;
  }

  sprintf(BANNER,"%s: read simulated CC mu-prior vs. z", fnam);
  print_banner(BANNER);


  // allocate memory and read simFile
  read_simFile_CCprior(); 
  print_eventStats(EVENT_TYPE); // call this before store_INFO_CCPRIOR
  store_INFO_CCPRIOR_CUTS();    // transfer INFO to smaller arrays

  // set z-bin index for each event to speed things up in fitting
  setup_zbins_CCprior(&INFO_CCPRIOR.TABLEVAR_CUTS, 
		      &INFO_CCPRIOR.MUZMAP.ZBIN );
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    setup_MUZMAP_CCprior(idsample, &INFO_CCPRIOR.TABLEVAR_CUTS,
			 &INFO_CCPRIOR.MUZMAP );
  }


  //  debugexit(fnam); // xxx REMOVE

  // ----------------------------------------------
  printf("\n Finished preparing CC prior. \n");
  fflush(stdout);

  return ;

} // end prepare_CCprior



// =========================================
void  read_simFile_CCprior(void) {

  // Created June 2019
  // Read simFile(s) and load INFO_CCPRIOR struct.
  //

  int  EVENT_TYPE = EVENT_TYPE_CCPRIOR ;
  int  NFILE      = INPUTS.nfile_CCprior;
  int  NEVT[MXFILE_CCPRIOR], NEVT_TOT;
  int  IFILETYPE, NVAR_ORIG, LEN_MALLOC, NROW, ifile, ISTART, isn ;
  char *simFile ;
  //  char fnam[] = "read_simFile_CCprior" ;

  // -------------------- BEGIN --------------------

  // do quick read of each file get NEVT count for each file
  NEVT_TOT = 0 ;
  for(ifile=0; ifile < NFILE; ifile++ ) {
    simFile = INPUTS.simFile_CCprior[ifile];
    NEVT[ifile] = SNTABLE_NEVT(simFile,TABLENAME_FITRES); 
    NEVT_TOT += NEVT[ifile];
    printf("\t Found %d events in %s. \n", NEVT[ifile], simFile) ;
  }


  // malloc enough memory to read all files
  LEN_MALLOC = NEVT_TOT + 10;
  INFO_CCPRIOR.TABLEVAR.LEN_MALLOC = LEN_MALLOC ; ;
  malloc_INFO_CCPRIOR(+1, LEN_MALLOC, 0);

  CUTMASK_POINTER[EVENT_TYPE]       = &INFO_CCPRIOR.TABLEVAR.CUTMASK[0];
  NALL_CUTMASK_POINTER[EVENT_TYPE]  = &INFO_CCPRIOR.TABLEVAR.NSN_ALL;
  NPASS_CUTMASK_POINTER[EVENT_TYPE] = &INFO_CCPRIOR.TABLEVAR.NSN_PASSCUTS ;
  NREJECT_CUTMASK_POINTER[EVENT_TYPE] = &INFO_CCPRIOR.TABLEVAR.NSN_REJECT ;

  // loop again over each file and read them

  for(ifile=0; ifile < NFILE; ifile++ ) {
    simFile   = INPUTS.simFile_CCprior[ifile];
    IFILETYPE = TABLEFILE_OPEN(simFile,"read");
    NVAR_ORIG = SNTABLE_READPREP(IFILETYPE,"FITRES");
  
    ISTART = INFO_CCPRIOR.TABLEVAR.NSN_ALL;

    SNTABLE_READPREP_TABLEVAR(ifile, ISTART, NEVT[ifile], 
			      &INFO_CCPRIOR.TABLEVAR);

    NROW = SNTABLE_READ_EXEC();
    INFO_CCPRIOR.TABLEVAR.NSN_ALL += NROW ;

    INFO_CCPRIOR.TABLEVAR.EVENT_RANGE[ifile][0] = ISTART ;
    INFO_CCPRIOR.TABLEVAR.EVENT_RANGE[ifile][1] = ISTART + NROW - 1;
    sprintf(INFO_CCPRIOR.TABLEVAR.INPUT_FILE[ifile],"%s", simFile);
    store_input_varnames(ifile, &INFO_CCPRIOR.TABLEVAR) ;
  }

  for(isn=0; isn < INFO_CCPRIOR.TABLEVAR.NSN_ALL; isn++ )  { 
    compute_more_TABLEVAR(isn, &INFO_CCPRIOR.TABLEVAR ); 
  }


  return ;

} // end  read_simFile_CCprior



// ===================================================
void store_INFO_CCPRIOR_CUTS(void) {


  // Created Jan 2016
  // 
  // Apply cuts to CCPRIOR sample, then copy CC-subset
  // to separate structures so that large (ALL) arrays 
  // can be de-alloated.
  //

  int  NSN_ALL      = INFO_CCPRIOR.TABLEVAR.NSN_ALL ;
  int  NSN_PASS     = INFO_CCPRIOR.TABLEVAR.NSN_PASSCUTS ;
  int  NBINa        = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  NBINb        = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  NBINg        = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;

  //  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  //  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;

  int  isn, icc, ia, ib, ig, ipar, cutmask;
  char *name ;
  char fnam[] = "store_INFO_CCPRIOR_CUTS" ;

  // ------------- BEGIN -------------

  printf("  %s: %d(ALL) -> %d(CUTS) \n", fnam, NSN_ALL, NSN_PASS);
  fflush(stdout);

  // allocate memory to hold events with cuts
  malloc_INFO_CCPRIOR(+1, 0, NSN_PASS);

  INFO_CCPRIOR.TABLEVAR_CUTS.NSN_ALL      = NSN_PASS ;
  INFO_CCPRIOR.TABLEVAR_CUTS.NSN_PASSCUTS = NSN_PASS ;
  
  // count number of true CC with cuts
  icc=0;
  for(isn=0; isn < NSN_ALL; isn++ ) {

    cutmask  = INFO_CCPRIOR.TABLEVAR.CUTMASK[isn];
    if ( cutmask ) { continue; }

    name = INFO_CCPRIOR.TABLEVAR.name[icc];
    sprintf(INFO_CCPRIOR.TABLEVAR_CUTS.name[icc],"%s", name);

    INFO_CCPRIOR.TABLEVAR_CUTS.zhd[icc] = 
      INFO_CCPRIOR.TABLEVAR.zhd[isn];
    INFO_CCPRIOR.TABLEVAR_CUTS.x0[icc]  = 
      INFO_CCPRIOR.TABLEVAR.x0[isn];

    for(ipar=0; ipar < NLCPAR; ipar++ ) { 
      INFO_CCPRIOR.TABLEVAR_CUTS.fitpar[ipar][icc] = 
	INFO_CCPRIOR.TABLEVAR.fitpar[ipar][isn];
    }

    INFO_CCPRIOR.TABLEVAR_CUTS.SIM_NONIA_INDEX[icc]  = 
      INFO_CCPRIOR.TABLEVAR.SIM_NONIA_INDEX[isn];

    INFO_CCPRIOR.TABLEVAR_CUTS.IDSURVEY[icc] = 
      INFO_CCPRIOR.TABLEVAR.IDSURVEY[isn];

    INFO_CCPRIOR.TABLEVAR_CUTS.IDSAMPLE[icc] = 
      INFO_CCPRIOR.TABLEVAR.IDSAMPLE[isn];


    for(ia = 0; ia < NBINa; ia++ ) {
      for(ib = 0; ib < NBINb; ib++ ) {
	for(ig = 0; ig < NBINg; ig++ ) {
	  for(ipar=0; ipar < NLCPAR; ipar++ ) { 
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].VAL[ipar] =
	      INFO_CCPRIOR.FITPARBIAS_ALPHABETA[isn][ia][ib][ig].VAL[ipar] ;
	    
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].ERR[ipar] =
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA[isn][ia][ib][ig].ERR[ipar] ;

	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].RMS[ipar] =
	      INFO_CCPRIOR.FITPARBIAS_ALPHABETA[isn][ia][ib][ig].RMS[ipar] ;

	  } //end ipar
	} // end ig
      } // end ib
    }  // end ia
    
    icc++ ;
  }    // end isn loop over all events
  
  
  // free ALL memory
  printf("  %s: free memory for ALL \n", fnam); fflush(stdout);
  malloc_INFO_CCPRIOR(-1,INFO_CCPRIOR.TABLEVAR.LEN_MALLOC,0);


  return ;
  
} // end store_INFO_CCPRIOR_CUTS


// =================================================
void setup_zbins_CCprior(TABLEVAR_DEF *TABLEVAR, BININFO_DEF *ZBIN) {

  // Setup z-grid for CCprior, and label each data point with
  // iz index to save time in fitting. Note that this z binning
  // can be different from the user z-binning.
  //
  // Jan 24 2020: check optional user input from nzbin_ccprior

  int    nzbin_ccprior = INPUTS.nzbin_ccprior ;
  double zrange = INPUTS.zmax - INPUTS.zmin ;
  double DZBIN_CCPRIOR = 0.10 ;

  double z, zlo, zhi, z0, z1, d_nbz ;
  int  nbz, icc, iz ;
  char MSGERR[100];
  char fnam[] = "setup_zbins_CCprior" ;

  // ------------- BEGIN -------------

  if ( nzbin_ccprior > 0 ) {
    nbz = nzbin_ccprior;  // user input nzbin
  }
  else {
    // use hard wired bin size of DZBIN_CCPRIOR
    d_nbz   = zrange/DZBIN_CCPRIOR ;
    nbz     = (int)d_nbz;
    if ( d_nbz > (double)nbz ) { nbz++ ; }
  }
  DZBIN_CCPRIOR = zrange/(double)nbz ;


  // make uniform z-grid, although a non-uniform grid is allowed.
  nbz = 0;
  z0 = INPUTS.zmin;  
  z1 = INPUTS.zmax - DZBIN_CCPRIOR + 1.0E-5;
  sprintf( ZBIN->varName,"z(CCprior)" );

  for ( z=z0; z < z1; z += DZBIN_CCPRIOR ) {
    zlo = z;
    zhi = z + DZBIN_CCPRIOR ;
    ZBIN->lo[nbz]  = zlo ;
    ZBIN->hi[nbz]  = zhi ;
    ZBIN->avg[nbz] = 0.5*(zlo+zhi);
    ZBIN->n_perbin[nbz] = 0 ;
    nbz++ ;  ZBIN->nbin = nbz;
  }
    
  for(icc=0; icc < TABLEVAR->NSN_ALL; icc++ ) {
    z = TABLEVAR->zhd[icc];
    sprintf(MSGERR,"%s: z=%f", fnam, z);
    iz = IBINFUN(z, ZBIN, 1, MSGERR );
    TABLEVAR->IZBIN[icc] = iz ;
    ZBIN->n_perbin[iz]++ ;
  }


  // print summary info for each z-bin
  printf("\n");
  printf("  %s: %d z-bins from %.3f to %.3f\n",
	 fnam, ZBIN->nbin, INPUTS.zmin, INPUTS.zmax);

  printf("\t   iz  zmin   zmax     N(NONIA) \n" );
  for(iz=0; iz < ZBIN->nbin; iz++ ) {
    printf("\t  %3d  %.3f  %.3f  %5d  \n",
	   iz, ZBIN->lo[iz], ZBIN->hi[iz], ZBIN->n_perbin[iz] );
  }

  fflush(stdout);

  return ;

} // end setup_zbins_CCprior


// ================================================
void setup_MUZMAP_CCprior(int IDSAMPLE, TABLEVAR_DEF *TABLEVAR,
			  MUZMAP_DEF *MUZMAP) {

  // setup dmu vs. z map binning in each z bins.
  // Do test with input alpha,beta,M0.
  // Only MUZMAP is used; INFO_CC is just used to test the PDF function.
  //

#define DMUMIN_CCPRIOR  -1.5  // add buffer bin to allow for interp
#define DMUMAX_CCPRIOR  +4.0  // add buffer bin to allow for interp
#define DMUBIN_CCPRIOR   0.5  

  int NBINa, NBINb, ia, ib, nbin_dmu ;
  double SUMa=0.0, SUMb=0.0, dmu, dmulo, dmuhi ;

  int NSN_BIASCOR = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  BININFO_DEF *BININFO_SIM_ALPHA = &INFO_BIASCOR.BININFO_SIM_ALPHA;
  BININFO_DEF *BININFO_SIM_BETA  = &INFO_BIASCOR.BININFO_SIM_BETA ;

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
  
  if ( IDSAMPLE == 0 ) {
    printf("\n");
    printf("  %s: %d dMU-bins from %.3f to %.3f\n",
	   fnam, MUZMAP->DMUBIN.nbin, DMUMIN_CCPRIOR, DMUMAX_CCPRIOR);
    fflush(stdout);
  }

  // print DMU distribution for first redshift bin
  MUZMAP->alpha = INPUTS.parval[IPAR_ALPHA0] ;
  MUZMAP->beta  = INPUTS.parval[IPAR_BETA0] ;
  MUZMAP->M0    = INPUTS.nommag0 ;

  MUZMAP->cosPar[0] = INPUTS.parval[IPAR_OL] ;   // OL
  MUZMAP->cosPar[1] = INPUTS.parval[IPAR_Ok] ;   // Ok
  MUZMAP->cosPar[2] = INPUTS.parval[IPAR_w0] ;   // w0
  MUZMAP->cosPar[3] = INPUTS.parval[IPAR_wa] ;   // wa
    
  // if there is a biasCor map, set alpha,beta to the average
  // since that should be a better estimate in case user input
  // starts alpha,beta way off.
  if ( NSN_BIASCOR > 0 ) {
    NBINa   = (*BININFO_SIM_ALPHA).nbin ;
    NBINb   = (*BININFO_SIM_BETA).nbin ;

    for(ia=0; ia<NBINa; ia++ )
      { SUMa += (*BININFO_SIM_ALPHA).avg[ia]; }
    for(ib=0; ib<NBINb; ib++ )
      { SUMb += (*BININFO_SIM_BETA).avg[ib]; }

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
  setup_DMUPDF_CCprior(IDSAMPLE, TABLEVAR, MUZMAP );  
  
  // -------------------------------
  // print sumary of MU-bins
  int IZTEST = 2;
  dump_DMUPDF_CCprior(IDSAMPLE, IZTEST, MUZMAP) ;
    
  return ;

}  // end setup_MUZMAP_CCprior


// =============================================
void setup_DMUPDF_CCprior(int IDSAMPLE, TABLEVAR_DEF *TABLEVAR,
			  MUZMAP_DEF *MUZMAP ) {

  // Created Feb 2016
  // for input redshift (z) and DMU map (*MUZMAP),
  // return the normalized PDF of the DMU distribution (*DMUPDF)
  // Compute PDF at center of two neighboring z bins 
  // and interpolate.
  //
  // This function is called for each fcn call, so no print statements !
  
  int imu, NMUBIN, NZBIN, ia, ib, ig, iz, icc, idsample ;
  int NCC[MXz][MAXMUBIN], NCC_SUM[MXz]; 
  double SUM_DMU[MXz], SUMSQ_DMU[MXz];
  double z, c, x1, mB, dl, mu, mumodel, dmu, a, b, gDM, M0 ;
  double XMU, XCC, XCC_SUM, DMUBIN ;
  char *name ;

  // muBias declarations
  bool USE_BIASCOR  = ( INFO_BIASCOR.TABLEVAR.NSN_ALL > 0 ) ;
  int  NDIM_BIASCOR = INFO_BIASCOR.NDIM ;
  BIASCORLIST_DEF     BIASCORLIST ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
  FITPARBIAS_DEF   FITPARBIAS_TMP[MXa][MXb][MXg] ; 
  double           MUCOVSCALE_TMP[MXa][MXb][MXg] ; 
  double fitParBias[NLCPAR], muBias, muBiasErr, muCOVscale ;

  char fnam[] = "setup_DMUPDF_CCprior" ;
    
  // ----------------- BEGIN ---------------
  
  if ( INFO_CCPRIOR.USEH11 ) { return ; }
    
  // number of DMU bins to make PDF
  NMUBIN = MUZMAP->DMUBIN.nbin ;
  NZBIN  = MUZMAP->ZBIN.nbin ;
  XMU    = (double)NMUBIN ;
  DMUBIN = MUZMAP->DMUBIN.hi[0] - MUZMAP->DMUBIN.lo[0]; 
  
  // init default PDF to flat in case there are not enough CC sim stats

  for(iz=0; iz < NZBIN; iz++ ) {
    NCC_SUM[iz]   = 0 ;
    SUM_DMU[iz]   = 0.0 ;
    SUMSQ_DMU[iz] = 0.0 ;
    MUZMAP->DMUAVG[IDSAMPLE][iz] = 0.0 ;
    MUZMAP->DMURMS[IDSAMPLE][iz] = 9999.0 ; // --> flat distribution
    
    for(imu=0; imu < NMUBIN; imu++ ) { 
      NCC[iz][imu] = 0 ;
      MUZMAP->NDMUPDF[IDSAMPLE][iz][imu] = 0 ;
      MUZMAP->DMUPDF[IDSAMPLE][iz][imu]  = 1.0/XMU ;
      MUZMAP->DMUPDF[IDSAMPLE][iz][imu] /= DMUBIN ; // Jul 2016
    }
  }
  
    
  for(ia=0; ia < MXa; ia++ ) {
    for(ib=0; ib < MXb; ib++ ) { 
      for(ig=0; ig < MXg; ig++ ) {
	MUCOVSCALE_TMP[ia][ib][ig] = 1.0 ; 
      }
    }
  }

  //  printf(" xxx %s: NROW = %d \n", fnam, INFO_CC->NROW);

  a    = MUZMAP->alpha ;
  b    = MUZMAP->beta ;
  gDM  = MUZMAP->gammadm;
  M0   = MUZMAP->M0 ;

  if ( NDIM_BIASCOR >= 5 ) 
    { get_INTERPWGT_abg(a,b,gDM, 0, &INTERPWGT,fnam ); } // returns INTERPWGT
   
  for(icc=0; icc < TABLEVAR->NSN_ALL ; icc++ ) {  

    // require correct IDSAMPLE to continue
    idsample = TABLEVAR->IDSAMPLE[icc] ;
    if( idsample != IDSAMPLE ) { continue ; }

    // keep only CC events in current IZ bin
    iz   = TABLEVAR->IZBIN[icc] ;   // redshift bin
    z    = TABLEVAR->zhd[icc] ;
    c    = TABLEVAR->fitpar[INDEX_c][icc] ;
    x1   = TABLEVAR->fitpar[INDEX_x1][icc] ;
    mB   = TABLEVAR->fitpar[INDEX_mB][icc] ;
    name = TABLEVAR->name[icc] ;
    
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
      BIASCORLIST.gammadm          = gDM ;
      BIASCORLIST.idsample         = idsample ;

      // transfer FITPARBIAS_ALPHABETA to FITPARBIAS_TMP that has
      // the right 2D structure definition for get_muBias
      load_FITPARBIAS_CCprior(icc,FITPARBIAS_TMP);

      if ( NDIM_BIASCOR >= 5 ) {
	get_muBias(name, &BIASCORLIST,     // (I) misc inputs
		   FITPARBIAS_TMP,         // (I) bias vs. ia,ib
		   MUCOVSCALE_TMP,         // (I) muCOVscale vs. ia,ib
		   &INTERPWGT,             // (I) wgt at each a,b grid point
		   fitParBias,     // (O) interp bias on mB,x1,c
		   &muBias,        // (O) interp bias on mu
		   &muBiasErr,     // (O) stat-error on above
		   &muCOVscale );  // (O) scale bias on muCOV (not used below)
      }
      else if ( NDIM_BIASCOR == 1 ) {
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
    if ( icc < -5 ) {
      printf(" xxx -----------  icc=%d  ---------- \n", icc);
      printf(" xxx z=%.3f  c=%.3f  x1=%.3f  mB=%.3f  \n",
	      z, c, x1, mB );
      printf(" xxx alpha=%.3f  beta=%.3f  M0=%.3f \n", 
	     MUZMAP->alpha,  MUZMAP->beta,  MUZMAP->M0 );
      printf(" xxx dmu = %.3f - %.3f = %.3f  (imu=%d)\n",
	     mu, mumodel, dmu, imu );
    }
    // xxxxxxxxxxxxxxxxxxxxxxx

    if ( imu >= 0 ) {
      NCC_SUM[iz]++ ;
      NCC[iz][imu]++ ;
      // increment sums for AVG & RMS ( July 2016)
      SUM_DMU[iz]   += dmu ;
      SUMSQ_DMU[iz] += (dmu*dmu) ;      
    }
  
  } // end icc loop


  
  // Normalize each DMU distribution to get PDF in each class,iz bin.
  // If the stats are too low then leave flat/default PDF.
  double S1TMP, S2TMP;   int NTMP ;  

  for(iz=0; iz < NZBIN; iz++ ) {

    if ( NCC_SUM[iz] > 5 ) {
      NTMP    = NCC_SUM[iz] ;
      XCC_SUM = (double)NTMP;

      // mean and rms
      S1TMP = SUM_DMU[iz];
      S2TMP = SUMSQ_DMU[iz];

      MUZMAP->DMUAVG[IDSAMPLE][iz] = S1TMP/XCC_SUM ;
      MUZMAP->DMURMS[IDSAMPLE][iz] = RMSfromSUMS(NTMP,S1TMP,S2TMP);
	
      // prob in each bin
      for(imu=0; imu < NMUBIN; imu++ ) { 
	XCC = (double)NCC[iz][imu];
	MUZMAP->DMUPDF[IDSAMPLE][iz][imu] = XCC / XCC_SUM ;
	MUZMAP->NDMUPDF[IDSAMPLE][iz][imu]= NCC[iz][imu];
	
	// divide by DMU bin size to get normalized integral
	// assuming that each DMU bin is linearly interpolated.
	MUZMAP->DMUPDF[IDSAMPLE][iz][imu] /= DMUBIN ; 
      }
    }

  }   // and iz loop
  
  return ;
  
} // end setup_DMUPDF_CCprior

// ============================================
void  dump_DMUPDF_CCprior(int IDSAMPLE, int IZ, MUZMAP_DEF *MUZMAP) {

  int imu, NCC ;
  double PDF;
  double z  = MUZMAP->ZBIN.avg[IZ] ;
  //  char fnam[] = "dump_DMUPDF_CCprior" ;

  // -------------- BEGIN ---------------

  printf("\n  %s : ", SAMPLE_BIASCOR[IDSAMPLE].NAME  );
  printf("  a=%.3f b=%.3f \n",
	 MUZMAP->alpha, MUZMAP->beta );

  printf("   <DMU>=%.3f  RMS(DMU)=%.3f \n"
	 , MUZMAP->DMUAVG[IDSAMPLE][IZ] 
	 , MUZMAP->DMURMS[IDSAMPLE][IZ]
	 );
  /*
  printf("  OL,Ok,w0,wa = %.3f, %.3f, %.3f, %.3f \n",	 
	 MUZMAP->cosPar[0], MUZMAP->cosPar[1],
	 MUZMAP->cosPar[2], MUZMAP->cosPar[3]  );  */
  

  printf("\t  imu  dMU(min) dMU(max)     PDF(z=%.2f)\n", z );

  for(imu=0; imu < MUZMAP->DMUBIN.nbin; imu++ ) {
    NCC = MUZMAP->NDMUPDF[IDSAMPLE][IZ][imu] ; 
    PDF = MUZMAP->DMUPDF[IDSAMPLE][IZ][imu] ;
    printf("\t  %3d  %6.3f  %6.3f       %.3f  (%d) \n"
	   ,imu, MUZMAP->DMUBIN.lo[imu], MUZMAP->DMUBIN.hi[imu]
	   ,PDF, NCC);
  }  
  fflush(stdout);

  return ;
  
} // end dump_DMUPDF_CCprior

// =======================
void  setup_contam_CCprior(void) {

  // Sep 26 2019
  // setup MURES bins and z bins for storing binned contamination info. 
  // Global CONTAM_INFO structures are filled.
  // This information is used to monitor contamination after the fit,
  // but is not used in the fit.

  int nb, i;
  double lo[MXz], hi[MXz];
  double tmp_lo, tmp_hi, tmp_avg, tmp_binsize;
  //  char fnam[] = "setup_contam_CCprior";

  // -------------- BEGIN ----------------

  // store with three hard-coded MURES bins

  sprintf(CONTAM_MURES_BINS.BININFO.varName,"MURES");
  nb=0;
  lo[nb] = -4.0;  hi[nb]=-0.5; nb++ ;
  lo[nb] = -0.5;  hi[nb]=+0.5; nb++ ;
  lo[nb] = +0.5;  hi[nb]=+4.0; nb++ ;

  for(i=0; i < nb; i++  ) {
    tmp_lo  = lo[i]; tmp_hi=hi[i];  tmp_binsize=tmp_hi-tmp_lo;
    tmp_avg = 0.5*(tmp_hi+tmp_lo);
    CONTAM_MURES_BINS.BININFO.lo[i]      = tmp_lo;
    CONTAM_MURES_BINS.BININFO.hi[i]      = tmp_hi;
    CONTAM_MURES_BINS.BININFO.avg[i]     = tmp_avg;
  }
  CONTAM_MURES_BINS.BININFO.nbin    = nb;
  CONTAM_MURES_BINS.BININFO.binSize = -9.0;  // N/A

  // next store redshift bins
  sprintf(CONTAM_REDSHIFT_BINS.BININFO.varName,"REDSHIFT");
  nb=4;
  double zbin = (INPUTS.zmax - INPUTS.zmin) / (double)nb;
  for(i=0; i < nb; i++  ) {
    tmp_lo  = INPUTS.zmin + zbin*(double)i;
    tmp_hi  = tmp_lo + zbin;
    tmp_avg = 0.5*(tmp_lo + tmp_hi);

    CONTAM_REDSHIFT_BINS.BININFO.lo[i]      = tmp_lo;
    CONTAM_REDSHIFT_BINS.BININFO.hi[i]      = tmp_hi;
    CONTAM_REDSHIFT_BINS.BININFO.avg[i]     = tmp_avg;
  }
  CONTAM_REDSHIFT_BINS.BININFO.nbin    = nb;
  CONTAM_REDSHIFT_BINS.BININFO.binSize = zbin;

  return ;

} // end setup_contam_CCprior

void zero_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO) {

  int i;
  int nbin = CONTAM_INFO->BININFO.nbin ;

  CONTAM_INFO->SUMPROB_TOT_IA    = 0.0 ;
  CONTAM_INFO->SUMPROB_TOT_CC    = 0.0 ;
  CONTAM_INFO->SUMPROB_TOT_RATIO = 0.0 ;

  CONTAM_INFO->NTRUE_TOT_IA    = 0;
  CONTAM_INFO->NTRUE_TOT_CC    = 0;
  CONTAM_INFO->TRUE_TOT_RATIO  = 0.0 ;

  for(i=0; i < nbin; i++ ) {
    CONTAM_INFO->sumProb_Ia[i] = 0.0 ;
    CONTAM_INFO->sumProb_cc[i] = 0.0 ;
    CONTAM_INFO->ntrue_Ia[i]   = 0 ;
    CONTAM_INFO->ntrue_cc[i]   = 0 ;
    CONTAM_INFO->true_ratio[i] = 0.0 ;
  }

} // end zero_contam_CCprior
	 
void sum_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO, double Prob_Ia,
			double xhisto, int SIM_NONIA_INDEX) {

  // Sep 26 2019:
  // Increment histogram of CC contamination vs. xhisto (mures,z, etc ...)

  int    ibin;
  bool   IS_SIM  =  INFO_DATA.TABLEVAR.IS_SIM ;
  double Prob_CC = 1.0 - Prob_Ia ;
  double sum_Ia, sum_cc, ratio ;
  char fnam[] = "sum_contam_CCprior";

  // -------------- BEGIN ----------------

  // un-binned sumprob 
  CONTAM_INFO->SUMPROB_TOT_IA   += Prob_Ia ;
  CONTAM_INFO->SUMPROB_TOT_CC   += Prob_CC ;
  sum_Ia = CONTAM_INFO->SUMPROB_TOT_IA;
  sum_cc = CONTAM_INFO->SUMPROB_TOT_CC;
  ratio  = sum_cc/(sum_Ia+sum_cc);
  CONTAM_INFO->SUMPROB_TOT_RATIO = ratio;

  
  // repeat for binning vs. xhisto
  // xxx mark delete  ibin = IBINFUN(xhisto, &CONTAM_INFO->BININFO, 1, fnam );
  ibin = IBINFUN(xhisto, &CONTAM_INFO->BININFO, 2, fnam );
  CONTAM_INFO->sumProb_Ia[ibin] += Prob_Ia ;
  CONTAM_INFO->sumProb_cc[ibin] += Prob_CC ;
  sum_Ia = CONTAM_INFO->sumProb_Ia[ibin];
  sum_cc = CONTAM_INFO->sumProb_cc[ibin];
  ratio  = sum_cc/(sum_Ia+sum_cc); 
  CONTAM_INFO->sumProb_ratio[ibin] = ratio;

  
  // - - - - - - SUM TRUTH for SIM DATA - - - - - 
  if ( IS_SIM ) {    
    if ( SIM_NONIA_INDEX == 0 ) {
      CONTAM_INFO->NTRUE_TOT_IA++ ;
      CONTAM_INFO->ntrue_Ia[ibin]++ ;
    }
    else {
      CONTAM_INFO->NTRUE_TOT_CC++ ;
      CONTAM_INFO->ntrue_cc[ibin]++ ;
    }


    sum_Ia = (double)CONTAM_INFO->NTRUE_TOT_IA;
    sum_cc = (double)CONTAM_INFO->NTRUE_TOT_CC;
    ratio  = sum_cc/(sum_Ia+sum_cc); 
    CONTAM_INFO->TRUE_TOT_RATIO = ratio;

    sum_Ia = (double)CONTAM_INFO->ntrue_Ia[ibin];
    sum_cc = (double)CONTAM_INFO->ntrue_cc[ibin];
    ratio  = sum_cc/(sum_Ia+sum_cc); 
    CONTAM_INFO->true_ratio[ibin] = ratio;

  } // end of IS_SIM
  

  return ;

} // end sum_contam_CCprior

void print_contam_CCprior(FILE *fp, CONTAM_INFO_DEF *CONTAM_INFO) {

  // print contamination total, and table vs. variable in BININFO.
  //
  // E.g.
  //   varName Range    CC/TOT(SUMPROB)     CC/TOT(NTRUE)
  //    ALL  
  //   -4.0 to -0.5
  //   -0.5 to +0.5
  //
  int  nbin     = CONTAM_INFO->BININFO.nbin;
  char *varName = CONTAM_INFO->BININFO.varName;
  int   IS_SIM  =  (INFO_DATA.TABLEVAR.IS_SIM == true) ;
  int  i;
  double lo, hi, xnIa, xncc, ratio, true_ratio;
  int    ntrue_cc, ntrue_Ia;
  char cRange[40], str_contam_data[80], str_contam_true[80];
  //  char fnam[] = "print_contam_CCprior";

  // -------------- BEGIN --------------
  
  fprintf(fp,"\n# CC Contamination vs. %s \n", varName);

  fprintf(fp,"#  %8s Range     CC/TOT(SUMPROB)         CC/TOT(TRUE) \n", 
	  varName);
  fprintf(fp,"# %s \n", dashLine);
  for(i=-1; i < nbin; i++ ) {  // -1 ==> all 
    str_contam_data[0]  =  str_contam_true[0] = 0 ;

    if ( i < 0 ) {
      // total sum
      sprintf(cRange, "TOTAL");
      xncc  = CONTAM_INFO->SUMPROB_TOT_CC ;
      xnIa  = CONTAM_INFO->SUMPROB_TOT_IA ;
      ratio = CONTAM_INFO->SUMPROB_TOT_RATIO ;
      ntrue_cc   = CONTAM_INFO->NTRUE_TOT_CC;
      ntrue_Ia   = CONTAM_INFO->NTRUE_TOT_IA;
      true_ratio = CONTAM_INFO->TRUE_TOT_RATIO ;
    } else {
      // specific bin of varName
      lo = CONTAM_INFO->BININFO.lo[i] ;
      hi = CONTAM_INFO->BININFO.hi[i] ;
      sprintf(cRange,"%5.2f to %5.2f", lo, hi); 
      xncc  = CONTAM_INFO->sumProb_cc[i] ;
      xnIa  = CONTAM_INFO->sumProb_Ia[i] ;
      ratio = CONTAM_INFO->sumProb_ratio[i] ;

      ntrue_cc   = CONTAM_INFO->ntrue_cc[i];
      ntrue_Ia   = CONTAM_INFO->ntrue_Ia[i];
      true_ratio = CONTAM_INFO->true_ratio[i] ;
    }

    sprintf(str_contam_data,"%4.1f/%7.1f = %.3f", 
	    xncc, xnIa+xncc, ratio);

    if ( IS_SIM ) {
      sprintf(str_contam_true,"%4d/%7d = %.3f", 
	      ntrue_cc, ntrue_Ia+ntrue_cc, true_ratio);
    }

    fprintf(fp,"#  %14s   %s    %s\n", 
	    cRange, str_contam_data, str_contam_true);

  } // end loop over bins

  fprintf(fp,"# %s \n\n", dashLine);
  fflush(fp);

} // end print_contam_CCprior

// ===================================================
void load_FITPARBIAS_CCprior(int icc, FITPARBIAS_DEF
			     (*FITPARBIAS_ABGRID)[MXb][MXg] )
{

  // Created Jun 7, 2016
  // load FITPARBIAS struct from global function of icc= CC index.
  // The global structure does not have the right form to pass
  // to get_muBias, so load FITPARBIAS that can be passed to get_muBias.

  int  NBINa, NBINb, NBINg, ia, ib, ig, ipar ;

  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;  

  //  char fnam[] = "load_FITPARBIAS_CCprior" ;

  // -------- BEGIN ------------

  NBINa   = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  NBINb   = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  NBINg   = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  

  for(ia=0; ia<NBINa; ia++ ) {
    for(ib=0; ib<NBINb; ib++ ) {
      for(ig=0; ig<NBINg; ig++ ) {
	for(ipar=ILCPAR_MIN; ipar<= ILCPAR_MAX; ipar++ ) {

	  FITPARBIAS_ABGRID[ia][ib][ig].VAL[ipar] =
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].VAL[ipar] ;
	  FITPARBIAS_ABGRID[ia][ib][ig].ERR[ipar] =
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].ERR[ipar] ;
	  FITPARBIAS_ABGRID[ia][ib][ig].RMS[ipar] =
	    INFO_CCPRIOR.FITPARBIAS_ALPHABETA_CUTS[icc][ia][ib][ig].RMS[ipar] ;
	
	} //ipar
      } //ig
    } // ib
  } //ia
      
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
  char  *name;
  double z,zerr, zsq ;
  double prob  = 0.0 ;
  double arg, sqsigCC, sigCC, sigintCC, CCpoly, dmuCC ;
  double alpha, beta, gamma, muerrsq_raw, covmat_fit[NLCPAR][NLCPAR] ;
  int i,i2;
  //  char fnam[] = "prob_CCprior_H11" ;

  // ----------- BEGIN ------------

  name  = INFO_DATA.TABLEVAR.name[n];
  z     = INFO_DATA.TABLEVAR.zhd[n];
  zerr  = INFO_DATA.TABLEVAR.zhderr[n];
  for(i=0; i < NLCPAR; i++ ) {
    for(i2=0; i2 < NLCPAR; i2++ )
      { covmat_fit[i][i2] = (double)INFO_DATA.TABLEVAR.covmat_fit[n][i][i2];}
  }
  
  zsq      = z*z ;
  CCpoly   = zPoly_mu[0]  + (zPoly_mu[1]*z)  + (zPoly_mu[2]*zsq) ;
  sigintCC = zPoly_sig[0] + (zPoly_sig[1]*z) + (zPoly_sig[2]*zsq) ;

  // compute muerrsq without intrinsic scatter
  alpha  = INPUTS.parval[IPAR_ALPHA0];  // fix a,b as in Jone 2016
  beta   = INPUTS.parval[IPAR_BETA0] ;
  gamma  = INPUTS.parval[IPAR_GAMMA0] ;

  muerrsq_raw = fcn_muerrsq(name, alpha, beta, gamma, covmat_fit,
			    z,zerr, 0);

  sqsigCC = muerrsq_raw + (sigintCC*sigintCC) ;
  dmuCC   = (dmu-CCpoly) ;  // i.e., mucalc -> mucalc+CCpoly
  arg     = 0.5*(dmuCC*dmuCC)/sqsigCC ;
  
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
  
  IZ   = IBINFUN(z, &MUZMAP->ZBIN, 0, fnam);

  if ( FUNDMU == FUNDMU_GAUSS ) {
    // make Guassian approximation with <DMU> and RMS(DMU).
    double AVG, RMS, resid, arg ;
    AVG    = MUZMAP->DMUAVG[IDSAMPLE][IZ];
    RMS    = MUZMAP->DMURMS[IDSAMPLE][IZ];
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
			&MUZMAP->DMUPDF[IDSAMPLE][IZ][0],
			fnam) ;
  }



  // xxxxxxxxxxxxxxxxxxx
  if ( IZ == -2 && dmu_local > 1.0 ) {
    printf(" xxx IDSAMPLE=%d  dmu_local=%6.3f  prob=%.3f \n",
	   IDSAMPLE, dmu_local, prob); fflush(stdout);
  }
  // xxxxxxxxxxxxxxxxxxx
 
  // -------------------------------
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


// =========================================
void print_eventStats(int event_type) {

  // Created Jun 1 2019
  // for input "event_type" (DATA,BIASCOR, or CCPRIOR),
  // Print number of events: total, pass, reject.
  // Note that NPASS is computed and stored here/
  // 
  // Tricky part is counting how many events are rejected by
  // a single cut ... and to to it efficiently.

  char *STRTYPE     = STRING_EVENT_TYPE[event_type];
  int  NSN_TOT, *CUTMASK_PTR ;
  int  NSN_REJ, *NSN_PASS, *NBIT, bit, isn, cutmask, NCUT ;
  int  NCUT_SOLO[MXCUTBIT];
  char fnam[] = "print_eventStats" ;

  // ---------- BEGIN ------------

  NSN_TOT      = *NALL_CUTMASK_POINTER[event_type];
  NSN_REJ      = *NREJECT_CUTMASK_POINTER[event_type];
  NSN_PASS     =  NPASS_CUTMASK_POINTER[event_type]; // pointer filled below
  CUTMASK_PTR  =  CUTMASK_POINTER[event_type];
  NBIT         =  &NSTORE_CUTBIT[event_type][0];

  // loop over all events and for each cutbit, count how many
  // events were cut by ONLY this one cut (i.e.. passed all other cuts)
  for(bit=0; bit < MXCUTBIT; bit++ ) { NCUT_SOLO[bit]=0; }

  for(isn=0; isn < NSN_TOT; isn++ ) {
    cutmask = CUTMASK_PTR[isn];
    if ( cutmask ==0 ) { continue; }
    // check if one and only 1 bit is set
    if ( (cutmask & (cutmask-1) ) == 0 ) {
      for(bit=0; bit < MXCUTBIT; bit++ ) {
	if ( cutmask & CUTMASK_LIST[bit] ) 
	  { NCUT_SOLO[bit]++; bit=MXCUTBIT+1; }
      }
    }
  }

  // ----------------------------
  *NSN_PASS = NSN_TOT - NSN_REJ ;

  printf("\n#%s\n", dashLine);
  printf(" %s STAT SUMMARY: %d(TOTAL) = %d(ACCEPT) + %d(REJECT) \n",
	 STRTYPE, NSN_TOT, *NSN_PASS, NSN_REJ);

  for(bit=0; bit < MXCUTBIT; bit++ ) {
    NCUT = NBIT[bit] ;
    if ( NCUT > 0 ) {
      printf(" %s NCUT[%2.2d] = %6d(ALL) %6d(onlyCut)   [%s] \n",
	     STRTYPE, bit, NCUT, NCUT_SOLO[bit], CUTSTRING_LIST[bit] );
    }
  }

  if ( *NSN_PASS == 0 ) {
    sprintf(c1err,"All %s events fail cuts.", STRTYPE);
    sprintf(c2err,"Check NCUT vs. cut listed above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("#%s\n\n", dashLine);
  fflush(stdout);

  return ;

} // end  print_eventStats


// ==================================================
void set_CUTMASK(int isn, TABLEVAR_DEF *TABLEVAR ) {

  // Created June, 2019 by R.Kessler
  // Set 'cutmask'  for this 'isn' and 'event_type'
  // Also check CUTFIELD (May 2020)
  //
  // Inputs: 
  //   isn        -> SN index
  //   event_type -> data, biasCor or CCprior
  //
  // Jun 25 2019: fux bug setting CUTBIT_IDSAMPLE for IS_BIASCOR
  //
  int event_type = TABLEVAR->EVENT_TYPE;
  int IS_DATA    = ( event_type == EVENT_TYPE_DATA );
  int IS_BIASCOR = ( event_type == EVENT_TYPE_BIASCOR );
  int IS_CCPRIOR = ( event_type == EVENT_TYPE_CCPRIOR );
  int NFIELD = INPUTS.NFIELD ;

  //  int  LDMP = 0;
  int  LCUTWIN, DOFLAG_CUTWIN[MXCUTWIN], icut, outside ;
  int  CUTMASK, REJECT, ACCEPT ;
  int  sntype, FOUND, SIM_NONIA_INDEX, idsample ;
  bool BADERR=false, BADCOV=false ;
  double cutvar_local[MXCUTWIN];
  double z, x0, x1, c, logmass, x0err, x1err, cerr  ;
  double COV_x0x1, COV_x0c, COV_x1c,  mBerr ;
  char   *name ;
  char fnam[]=  "set_CUTMASK";

  // ---------- BEGIN ---------

  // strip off local variables

  for ( icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
    DOFLAG_CUTWIN[icut] = (int)TABLEVAR->DOFLAG_CUTWIN[icut];
    cutvar_local[icut]  = (double)TABLEVAR->CUTVAL[icut][isn] ; 
  }

    
  name      =  TABLEVAR->name[isn];
  idsample  =  (int)TABLEVAR->IDSAMPLE[isn] ;
  sntype    =  (int)TABLEVAR->SNTYPE[isn] ;
  
  z         =  (double)TABLEVAR->zhd[isn];
  x0        =  (double)TABLEVAR->x0[isn];  
  x1        =  (double)TABLEVAR->fitpar[INDEX_x1][isn] ;
  c         =  (double)TABLEVAR->fitpar[INDEX_c ][isn] ;  
  logmass   =  (double)TABLEVAR->logmass[isn];
  x0err     =  (double)TABLEVAR->x0err[isn] ;
  mBerr     =  (double)TABLEVAR->fitpar_err[INDEX_mB][isn] ;
  x1err     =  (double)TABLEVAR->fitpar_err[INDEX_x1][isn] ;
  cerr      =  (double)TABLEVAR->fitpar_err[INDEX_c ][isn] ;
    
  COV_x0x1  = (double)TABLEVAR->COV_x0x1[isn];
  COV_x0c   = (double)TABLEVAR->COV_x0c[isn];
  COV_x1c   = (double)TABLEVAR->COV_x1c[isn];

  SIM_NONIA_INDEX = (int)TABLEVAR->SIM_NONIA_INDEX[isn] ;


  // =======================================
  // check CUTWIN options
  LCUTWIN = apply_CUTWIN(event_type, DOFLAG_CUTWIN, cutvar_local);
  if ( LCUTWIN == 0 ) { setbit_CUTMASK(isn, CUTBIT_CUTWIN, TABLEVAR); }

  // -----------------
  // check CUTFIELD (May 2020)
  if ( NFIELD > 0 ) {
    bool MATCH = false;
    char *tmpField, *field = TABLEVAR->field[isn];  
    for(icut=0; icut < NFIELD; icut++ ) {
      tmpField = INPUTS.FIELDLIST[icut] ;
      if ( strstr(field,tmpField) != NULL )   {  MATCH=true; }
    }
    if ( !MATCH ) { setbit_CUTMASK(isn, CUTBIT_FIELD, TABLEVAR); }
  }

  // -----------------------------------------
  // apply legacy cuts 
  // These cuts are redundant with above LCUTWIN to allow older inputs.

  if ( (z < INPUTS.zmin) || ( z > INPUTS.zmax)) 
    { setbit_CUTMASK(isn, CUTBIT_z, TABLEVAR);  }

  if ( (x1 < INPUTS.x1min) || (x1 > INPUTS.x1max)) 
    { setbit_CUTMASK(isn, CUTBIT_x1, TABLEVAR); }

  if ( (c < INPUTS.cmin) || ( c > INPUTS.cmax) ) 
    { setbit_CUTMASK(isn, CUTBIT_c, TABLEVAR); }

  if ( (logmass<INPUTS.logmass_min) || (logmass>INPUTS.logmass_max) ) 
    { setbit_CUTMASK(isn, CUTBIT_logmass, TABLEVAR);  }

  if ( x0err <= INPUTS.maxerr_abort_x0 ) { BADERR = true ; }
  if ( x1err <= INPUTS.maxerr_abort_x1 ) { BADERR = true ; }
  if (  cerr <= INPUTS.maxerr_abort_c  ) { BADERR = true ; }
  if ( BADERR ) {
    setbit_CUTMASK(isn, CUTBIT_BADERR, TABLEVAR);
    if ( IS_DATA ) {
      printf(" %s WARNING: "
	     "one or more x0,x1,c errors <=0 (SNID=%s)\n", fnam, name );
    }
  }


  // May 2016  flag crazy COV and mBerr 

  if ( fabs(COV_x0x1) > 8.0 ) { BADCOV = true ; } // COV_x0x1
  if ( fabs(COV_x0c ) > 8.0 ) { BADCOV = true ; } // COV_x0c
  if ( fabs(COV_x1c ) > 8.0 ) { BADCOV = true ; } // COV_x1c  
  if ( mBerr > 1.0  )         { BADCOV = true ; }
  
  if ( BADCOV ) {
    setbit_CUTMASK(isn, CUTBIT_BADCOV, TABLEVAR );
    if ( IS_DATA ) {
      printf(" %s WARNING: one or more crazy COV (SNID=%s). \n",
	     fnam, name);
    }
  }


  if ( TABLEVAR->IVAR_SNTYPE > 0 && INPUTS.Nsntype > 0 ) {
    FOUND = 0 ;
    for ( icut=0; icut < INPUTS.Nsntype; icut++ ) {
      if ( sntype == INPUTS.sntype[icut] ) { FOUND = 1; }
    }
    if ( FOUND==0 ) { setbit_CUTMASK(isn, CUTBIT_sntype, TABLEVAR); }
  }

   
  int CUT_IDSAMPLE=0;
  if ( idsample >= 0 ) {
    if (  SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { CUT_IDSAMPLE = 1;}
  }
  if ( CUT_IDSAMPLE ) { setbit_CUTMASK(isn, CUTBIT_IDSAMPLE, TABLEVAR ); }
  
  // - - - - - - - - - - - - - - - - - - - - - - -  -
  // check cuts specific to event_type

  if ( IS_DATA ) {
    // check pre-scale if data is from a simulation
    REJECT = prescale_reject_simData(SIM_NONIA_INDEX);
    if ( REJECT )
      { setbit_CUTMASK(isn, CUTBIT_SIMPS, TABLEVAR); }


    // djb we only want to apply the cid list for the data                 
    ACCEPT = selectCID_data(name);
    if ( !ACCEPT )
      { setbit_CUTMASK(isn, CUTBIT_CID, TABLEVAR); } //the mask is in tablevar

  }
  else if ( IS_BIASCOR ) { 
    
    if ( SIM_NONIA_INDEX != 0 ) 
      { setbit_CUTMASK(isn, CUTBIT_TRUESNIa, TABLEVAR); }
        
    if ( idsample < 0 )
      { setbit_CUTMASK(isn, CUTBIT_IDSAMPLE, TABLEVAR); }

    if ( prescale_reject_biasCor(isn) ) 
      { setbit_CUTMASK(isn, CUTBIT_SIMPS, TABLEVAR); }

    
    // check biasCor grid-cuts only if cutmask is not already set
    // (to avoid confusion with global cuts on z,x1,c)

    CUTMASK = TABLEVAR->CUTMASK[isn];
    if ( CUTMASK == 0 ) {
      outside = outside_biasCor_grid(isn);
      if ( outside == INDEX_z )  
	{ setbit_CUTMASK(isn, CUTBIT_zBIASCOR, TABLEVAR); }
      else if ( outside == INDEX_x1 ) 
	{ setbit_CUTMASK(isn, CUTBIT_x1BIASCOR, TABLEVAR); }
      else if ( outside == INDEX_c ) 
	{ setbit_CUTMASK(isn, CUTBIT_cBIASCOR, TABLEVAR); }
    }
    
  }
  else if ( IS_CCPRIOR ) {
   
    if ( SIM_NONIA_INDEX == 0 ) // require true SNCC
      { setbit_CUTMASK(isn, CUTBIT_TRUESNCC, TABLEVAR); }
    
    if ( idsample < 0 )
      { setbit_CUTMASK(isn, CUTBIT_IDSAMPLE, TABLEVAR); }

    if ( !storeBias_CCprior(isn) ) 
      { setbit_CUTMASK(isn, CUTBIT_BIASCOR, TABLEVAR); }    
  }
  
  return ;

} // end of set_CUTMASK 

// ============================================
void setbit_CUTMASK(int isn, int bitnum, TABLEVAR_DEF *TABLEVAR ) {

  // Inputs:
  //   bitnum     : bit to set
  //   TABLEVAR   : structure with arrays and cutmask
  //

  int  EVENT_TYPE   = TABLEVAR->EVENT_TYPE ;
  int  *CUTMASK     = &TABLEVAR->CUTMASK[isn];
  int  CUTMASK_SET  = CUTMASK_LIST[bitnum] ;
  //  char fnam[] = "setbit_CUTMASK" ;

  // ------------- BEGIN -------------

  if ( ( *CUTMASK & CUTMASK_SET) == 0 ) {

    if ( *CUTMASK == 0 ) { TABLEVAR->NSN_REJECT++ ;  }

    *CUTMASK |= CUTMASK_LIST[bitnum];
    TABLEVAR->NSN_CUTBIT[bitnum]++ ;
    NSTORE_CUTBIT[EVENT_TYPE][bitnum]++ ;
  }

  return ;

} // end setbit_CUTMASK 

int selectCID_data(char *cid){
  // Created Sep 5 2019 by D.Brout
  // for file= data. determines if cid is in cidlist_data

  //  char fnam[] = "selectCID_data";

  int ACCEPT = 1;
  int i;
  char *tmpCID;

  // ------- BEGIN -------------

  if (strlen(INPUTS.cidFile_data) == 0) { return ACCEPT; }

  for (i=0;i<INPUTS.ncidList_data;i++) {
    tmpCID = INPUTS.cidList_data[i];
    if (strcmp(cid,tmpCID)==0) {
      return ACCEPT;
    }
  }
  ACCEPT = 0;

  return ACCEPT;
} // END selectCID_data

// =============================================
int prescale_reject_simData(int SIM_NONIA_INDEX) {

  // Created Apr 14 2017 by R.Kessler
  // For simulated data, check prescale_simData and prescale_simCC
  // to see if this event is rejected.
  //   * prescale_simData pre-scales all sim events randomly
  //   * prescale_simCC   pre-scales only the true CC subset
  // Note that real data are always accepted.
  // Function returns 1 to reject, 0 to keep.
  //
  // Do not use this function for biasCor or CCprior.

  int REJECT = 0 ;
  float XN, XNPS;
  //  char fnam[] = "prescale_reject_simData" ;

  // ------------- BEGIN -----------------

  // return accept for real data
  if ( FOUNDKEY_SIM == 0 ) { return(REJECT); }

  // increment NSIMDATA counter.
  // Note that prescale_simData can be non-integer
  NSIMDATA++ ;
  XN    = (float)NSIMDATA ;
  XNPS  = (float)INPUTS.prescale_simData ;
  if ( fmodf( XN, XNPS ) != 0 )  { REJECT = 1 ; }

  // increment NSIMCC counter
  if ( SIM_NONIA_INDEX > 0 ) {
    NSIMCC++ ;
    XN    = (float)NSIMCC ;
    XNPS  = (float)INPUTS.prescale_simCC ;
    if ( fmodf( XN, XNPS ) != 0 ) { REJECT = 1 ; }
    if ( INPUTS.prescale_simCC > 9999.0 ) { REJECT=1 ; }
  }
  return(REJECT);

} // end prescale_reject_simData

// =============================================
int prescale_reject_biasCor(int isn) {

  // Created Jun 2019 
  //   [coded moved from obsolete biasMapSelect]
  // Returns 1 to reject; returns 0 to keep.
  //
  int  PS0, PS1;
  double Xisn, XPS1;
  //  char fnam[] = "prescale_reject_biasCor" ;

  // --------------- BEGIN ------------

  if ( INPUTS.prescale_biasCor[1] > 1 ) {
    PS0  =  INPUTS.prescale_biasCor[0] ;
    PS1  =  INPUTS.prescale_biasCor[1] ;
    Xisn =  (double)isn ;
    XPS1 =  (double)PS1 ;
    if ( fmod(Xisn,XPS1) != PS0 ) { return(1); } 
  }


  return(0) ;

}  // end prescale_reject_biasCor

// ====================================================
int outside_biasCor_grid(int isn) {

  // Created Jun 2019
  // For 5D biasCor grid, 
  //  --> Returns 1 (reject) if z,x1, or c is outside the valid 
  //      grid range of the biasCor.
  // --> Returns 0 --> keep it.

  int  APPLY = 0 ;
  int  idsample, NZ;
  double z, x1, c, zmin, zmax ;
  char *name;
  char fnam[] = "outside_biasCor_grid";

  // -------------- BEGIN -------------

  if (INPUTS.opt_biasCor & MASK_BIASCOR_5D)       { APPLY = 1; }
  if (INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT)  { APPLY = 1; }
  if ( APPLY == 0 ) { return(0); }

  z        = INFO_BIASCOR.TABLEVAR.zhd[isn];
  x1       = INFO_BIASCOR.TABLEVAR.fitpar[INDEX_x1][isn];
  c        = INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][isn];
  idsample = INFO_BIASCOR.TABLEVAR.IDSAMPLE[isn];
  name     = INFO_BIASCOR.TABLEVAR.name[isn];
  

  if ( idsample < 0 ) {
    sprintf(c1err,"Invalid idsample = %d ", idsample);
    sprintf(c2err,"isn=%d (SNID=%s)", isn, name);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // make sure values are within grid
  NZ   = CELLINFO_BIASCOR[idsample].BININFO_z.nbin ;
  zmin = CELLINFO_BIASCOR[idsample].BININFO_z.lo[0] ;
  zmax = CELLINFO_BIASCOR[idsample].BININFO_z.hi[NZ-1] ;
  if ( z  < zmin ) { return INDEX_z ; }
  if ( z  > zmax ) { return INDEX_z ; }

  if ( x1 < BIASCOR_MINVAL_LCFIT[INDEX_x1]  ) { return INDEX_x1; }
  if ( x1 > BIASCOR_MAXVAL_LCFIT[INDEX_x1]  ) { return INDEX_x1; }

  if ( c  < BIASCOR_MINVAL_LCFIT[INDEX_c]   ) { return INDEX_c; }
  if ( c  > BIASCOR_MAXVAL_LCFIT[INDEX_c]   ) { return INDEX_c; }

  return(0) ;

} // end outside_biasCor_grid





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
  bool SKIP, EXCEPTION;
  char *sptr;
  //  char fnam[] = "parse_parFile" ;

  // ------------------ BEGIN --------------

  // check for special mode to cat data files and do NOTHING else.
  if ( strcmp(parFile,"cat_only") == 0 ) 
    { INPUTS.cat_only = true; return ;  }

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

  char line[MXCHAR_LINE];
  while (fgets(line,MXCHAR_LINE,fdef)) {

    // skip blank lines and lines starting with comment char (RK 4/23/2012)
    if ( strlen(line) < 3 )       { continue ; }  
    if ( commentchar(line) == 1 ) { continue ; } // see sntools.c

    sptr = strtok(line,"\n");

    // skip keys with colon that are used by a master 
    // SALT2mu_fit.pl script (RK Apr 2013)
    // Aug 24 2016: make exception for key containing group_biasCor
    //    to allow colon-delimeter in input string
    // NOTE: should re-facto to check list of keys allowing colon
    SKIP      = ( strstr(sptr,":")             != NULL );
    EXCEPTION = ( strstr(sptr,"group_biascor") != NULL );
    if ( SKIP && !EXCEPTION ) { continue ; }

    /* xxx mark delete May 19 2020 xxxxxxxx
    if ( strstr(sptr,"group_biascor") == NULL &&
	 strstr(sptr,":") != NULL ) { continue ; }
    xxxxxxxx end mark xxxxxxx */

    ppar(sptr); // pass entire line
  }

  printf("\n");
  fclose(fdef);

  return ;

} // end  of parse_parFile


// **************************************
void override_parFile(int argc, char **argv) {

  // April 23, 2012 by R.Kessler
  // Moved from MAIN
  //
  // Read command-line arguments and override input file.
  //
  // Jan 15 2018: bug fix to work with CUTWIN; see i+=3
  // May 15 2019: increase item & line size to allow for
  //              comma-separated simfile_biascor
  //
  // Aug 6 2019: item[MXCHAR_LINE] -> *item, and remove obsolete
  //             abort on more than 200 characters.

  int  ntmp, i, found ;
  char *item, tmpLine[256];
  char fnam[] = "override_parFile";

  // ---------- BEGIN ------------

  ntmp = 0;
  printf("\n Parse command-line options: \n");
  uniqueOverlap("INIT","SALT2mu command-line override");

  for (i=2; i < argc; ++i) {

    item = argv[i];
    ntmp++;

    //  if ( strcmp(item,"CUTWIN") == 0 ) {
    if ( !strncmp(item,"CUTWIN",6) ) {  // allow CUTWIN(option)
      // glue together 4 contiguous words into one string
      sprintf(tmpLine,"%s %s %s %s", argv[i],argv[i+1],argv[i+2],argv[i+3] ) ;
      found = ppar(tmpLine);
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
  // Parses command line or input files

  // Dec 08, 2014 - read new zVARNAME
  // Aug 22, 2016 - call remove_quote
  // Apr 17, 2017 - return(1) if item is found; 0 otherwise
  // Jun 18, 2018 - fix INPUTS.ipar for H11 option
  // Apr 10, 2019 - replace !strcmp with uniqueOverlap function
  //                which aborts on duplicate key.
  
  int  ipar, len ;  
  char key[MXCHAR_VARNAME], *s, tmpString[60];
  char fnam[] = "ppar" ;

  // --------- BEGIN ----------

  printf(" Parse '%s' \n",item);  fflush(stdout);

  if ( uniqueOverlap(item,"prefix=") )
    { sscanf(&item[7],"%s",INPUTS.PREFIX); return(1); }

  if ( uniqueOverlap(item,"catfile_out=") ) {
    s = INPUTS.catfile_out ;
    sscanf(&item[12],"%s",s); remove_quote(s); return(1);
  }

  if ( uniqueOverlap(item,"snparfile_out=") ) {
    s = SNFUNPAR_CHI2INFO_OUTPUT.OUTFILE ;
    SNFUNPAR_CHI2INFO_OUTPUT.USE = true ;
    sscanf(&item[14],"%s",s); remove_quote(s); return(1);
  }
  
  if ( uniqueOverlap(item,"cutmask_write=") )
    { sscanf(&item[14],"%i", &INPUTS.cutmask_write ); return(1); }
  if ( uniqueOverlap(item,"errmask_write=") ) // allow legacy name
    { sscanf(&item[14],"%i", &INPUTS.cutmask_write ); return(1); }

  if ( uniqueOverlap(item,"file=") )
    { parse_datafile(&item[5]);  return(1);  }
  if ( uniqueOverlap(item,"datafile=") ) 
    { parse_datafile(&item[9]);  return(1);  }


  // -------------------------
  // allow two different keys to define biasCor file name
  if ( uniqueOverlap(item,"simfile_bias=") ) 
    { parse_simfile_biasCor(&item[13]);  return(1);  }
  if ( uniqueOverlap(item,"simfile_biascor=") ) 
    { parse_simfile_biasCor(&item[16]);  return(1); }

  // - - - - - - 

  if ( uniqueOverlap(item,"prescale_biascor=") ) 
    { parse_prescale_biascor(&item[17]); return(1); }
 
  if ( uniqueOverlap(item,"opt_biascor=")  )
    { sscanf(&item[12],"%d", &INPUTS.opt_biasCor);  return(1); }

  if ( uniqueOverlap(item,"sigmb_biascor=") )  // legacy name (9/2016)
    { sscanf(&item[14],"%le", &INPUTS.sigint_biasCor); return(1); }

  if ( uniqueOverlap(item,"sigint_biascor=")  )
    { sscanf(&item[15],"%le", &INPUTS.sigint_biasCor); return(1); }

  if ( uniqueOverlap(item,"snrmin_sigint_biascor=")  )
    { sscanf(&item[22],"%le", &INPUTS.snrmin_sigint_biasCor); return(1); }

  // allow two different keys for fieldgroup_biascor
  if ( uniqueOverlap(item,"fieldsplit_biascor=")  ) { // legacy name
    s=INPUTS.fieldGroup_biasCor ;
    sscanf(&item[19],"%s",s); remove_quote(s);
    return(1);
  }
  if ( uniqueOverlap(item,"fieldgroup_biascor=")  ) {
    s=INPUTS.fieldGroup_biasCor ;  
    sscanf(&item[19],"%s",s); remove_quote(s);
    return(1);
  }


  if ( uniqueOverlap(item,"surveygroup_biascor=")  ) {
    s=INPUTS.surveyGroup_biasCor ; 
    sscanf(&item[20],"%s",s); remove_quote(s);
    return(1);
  }

  if ( uniqueOverlap(item,"surveylist_nobiascor=")  ) {
    s=INPUTS.surveyList_noBiasCor ;
    sscanf(&item[21],"%s",s); remove_quote(s); 
    return(1);
  }

  if ( uniqueOverlap(item,"interp_biascor_logmass=")  )
    { sscanf(&item[23],"%d", &INPUTS.interp_biascor_logmass);  return(1); }

  if ( uniqueOverlap(item,"sigma_cell_biascor=") ) 
    { sscanf(&item[19],"%le", &INPUTS.sigma_cell_biasCor); return(1); }
  
  // -------- CC prior ------------
  if ( uniqueOverlap(item,"simfile_ccprior=")  ) 
    { parse_simfile_CCprior(&item[16]); return(1);  }

  if ( uniqueOverlap(item,"nzbin_ccprior=")) 
    { sscanf(&item[14],"%i",&INPUTS.nzbin_ccprior); return(1); }

  if ( uniqueOverlap(item,"maxprobcc_for_sigint=") )
    { sscanf(&item[21],"%le", &INPUTS.maxProbCC_for_sigint); return(1); } 

  if ( uniqueOverlap(item,"varname_pIa=")  ) {
    s = INPUTS.varname_pIa ;
    sscanf(&item[12],"%s",s); remove_quote(s);  return(1);
  }

  if ( uniqueOverlap(item,"append_varname_missing=")  ) {
    s = INPUTS.append_varname_missing ;    s[0]=0; // clobber default
    sscanf(&item[23],"%s",s); remove_quote(s);  return(1);
  }

  if ( uniqueOverlap(item,"force_pIa=")  ) { 
    sscanf(&item[10],"%s", tmpString); 
    if ( strcmp(tmpString,"perfect")==0 || strcmp(tmpString,"PERFECT")==0 )
      { INPUTS.perfect_pIa = true ;  INPUTS.ipar[IPAR_scalePCC]=0; }
    else
      { sscanf(tmpString, "%le", &INPUTS.force_pIa); }
    return(1);  
  }

  if ( uniqueOverlap(item,"typeIa_ccprior=") ) 
    { sscanf(&item[15],"%d", &INPUTS.typeIa_ccprior); return(1); } 

  if ( uniqueOverlap(item,"type_list_probcc0=") ) 
    { sscanf(&item[18],"%s",INPUTS_PROBCC_ZERO.str_type_list); return(1); } 
  if ( uniqueOverlap(item,"idsurvey_list_probcc0=") ) 
    { sscanf(&item[22],"%s",INPUTS_PROBCC_ZERO.str_idsurvey_list); return(1);} 

  if ( uniqueOverlap(item,"bins="))  // obsolete, but still allowed (use nzbin)
    { sscanf(&item[5],"%i",&INPUTS.nzbin); return(1); }
  if ( uniqueOverlap(item,"nzbin=")) 
    { sscanf(&item[6],"%i",&INPUTS.nzbin); return(1); }
  if ( uniqueOverlap(item,"nlogzbin=")) 
    { sscanf(&item[9],"%i",&INPUTS.nlogzbin); return(1); }
  if ( uniqueOverlap(item,"powzbin=")) 
    { parse_powzbin(&item[8]); return(1); }  
  if ( uniqueOverlap(item,"zbinuser=")) 
    { sscanf(&item[9],"%s", INPUTS.zbinuser); return(1); }

  if ( uniqueOverlap(item,"blindpar")) 
    { parse_blindpar(item); return(1); }

  if ( uniqueOverlap(item,"min_per_zbin=")) 
    { sscanf(&item[13],"%i", &INPUTS.min_per_zbin); return(1); }

  if ( uniqueOverlap(item,"zVARNAME=")) {
    s=INPUTS.varname_z;    sscanf(&item[9],"%s",s); remove_quote(s); 
    return(1);
  }


  // xxx mark delete   sprintf(key,"varname_gamma=");
  if ( uniqueOverlap(item,"varname_gamma=") ) {
    s=INPUTS.varname_gamma;  sscanf(&item[14],"%s",s); remove_quote(s); 
    return(1);
  }

  if ( uniqueOverlap(item,"varname_z=")) { 
    s=INPUTS.varname_z ;  sscanf(&item[10],"%s",s); remove_quote(s); 
    return(1);
  }

  if ( uniqueOverlap(item,"zmin=")) 
    { sscanf(&item[5],"%le",&INPUTS.zmin) ; return(1); }
  if ( uniqueOverlap(item,"zmax=")) 
    { sscanf(&item[5],"%le",&INPUTS.zmax) ; return(1); }

  if ( uniqueOverlap(item,"x1min=")) 
    { sscanf(&item[6],"%le",&INPUTS.x1min); return(1); }
  if ( uniqueOverlap(item,"x1max="))  
    { sscanf(&item[6],"%le",&INPUTS.x1max); return(1); }

  if ( uniqueOverlap(item,"cmin="))  
    { sscanf(&item[5],"%le",&INPUTS.cmin); return(1); }
  if ( uniqueOverlap(item,"cmax="))  
    { sscanf(&item[5],"%le",&INPUTS.cmax); return(1); }


  if ( uniqueOverlap(item,"logmass_min="))  
    { sscanf(&item[12],"%le",&INPUTS.logmass_min); return(1); }
  if ( uniqueOverlap(item,"logmass_max="))  
    { sscanf(&item[12],"%le",&INPUTS.logmass_max); return(1); }
  if ( uniqueOverlap(item,"nbin_logmass="))  
    { sscanf(&item[13],"%d",&INPUTS.nbin_logmass); return(1); }

  if ( uniqueOverlap(item,"chi2max=")) 
    { sscanf(&item[8],"%le",&INPUTS.chi2max); return(1); }

  if ( uniqueOverlap(item,"maxerr_abort_x0="))  
    { sscanf(&item[16],"%le",&INPUTS.maxerr_abort_x0); return(1); }
  if ( uniqueOverlap(item,"maxerr_abort_x1="))  
    { sscanf(&item[16],"%le",&INPUTS.maxerr_abort_x1); return(1); }
  if ( uniqueOverlap(item,"maxerr_abort_c="))  
    { sscanf(&item[15],"%le",&INPUTS.maxerr_abort_c); return(1); }
  

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

  if ( uniqueOverlap(item,"fieldlist=") ) 
    { parse_FIELDLIST(&item[10]); return(1); } 

  if ( uniqueOverlap(item,"idsample_select=") ) 
    {  sscanf(&item[16], "%s", INPUTS.idsample_select );  return(1); } 

  if ( uniqueOverlap(item,"cid_select_file=") ) // djb
    {  sscanf(&item[16], "%s", INPUTS.cidFile_data );  return(1); }


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


  /* xxx mark delete (params read below) xxxxxx
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
  if ( uniqueOverlap(item,"JOBID_SPLITRAN=")) 
    { sscanf(&item[15],"%d", &INPUTS.JOBID_SPLITRAN); return(1); }

  if ( uniqueOverlap(item,"iflag_duplicate=")) 
    { sscanf(&item[16],"%d", &INPUTS.iflag_duplicate ); return(1); }

  if ( uniqueOverlap(item,"dumpflag_nobiascor=")) 
    { sscanf(&item[19],"%d", &INPUTS.dumpflag_nobiasCor ); return(1); }

  if ( uniqueOverlap(item,"frac_warn_nobiascor=")) 
    { sscanf(&item[20],"%le", &INPUTS.frac_warn_nobiasCor ); return(1); }

  if ( uniqueOverlap(item,"snid_mucovdump=")) 
    { sscanf(&item[15],"%s", INPUTS.SNID_MUCOVDUMP); return(1); }

  return(0);
  
} // end ppar


// **************************************************
void parse_datafile(char *item) {

  // Created June 4 2019
  // parse comma-separate list of data files names.

  int ifile;
  int MEMC = MXCHAR_FILENAME*sizeof(char);
  char comma[] = ",";
  //  char fnam[]  = "parse_dataFile" ;

  // ------------------ BEGIN -----------------

  // first allocate memory for file names 
  INPUTS.dataFile = (char**)malloc( MXFILE_DATA*sizeof(char*));
  for(ifile=0; ifile < MXFILE_DATA; ifile++ ) 
    { INPUTS.dataFile[ifile] = (char*)malloc(MEMC); }

  // split item string
  splitString(item, comma, MXFILE_DATA,    // inputs
	      &INPUTS.nfile_data, INPUTS.dataFile ); // outputs 
  
  char *f0 = INPUTS.dataFile[0];
  if ( IGNOREFILE(f0) ) { INPUTS.nfile_data = 0 ; }

  return;

} // parse_datafile

// **************************************************
void parse_simfile_biasCor(char *item) {

  // Created May 15 2019
  // parse comma-separate list of biascor files names.

  int ifile;
  int MEMC = MXCHAR_FILENAME*sizeof(char);
  char comma[] = ",";
  //  char fnam[]  = "parse_simfile_biasCor" ;

  // ------------------ BEGIN -----------------

  // first allocate memory for file names
  INPUTS.simFile_biasCor = (char**)malloc( MXFILE_BIASCOR*sizeof(char*));
  for(ifile=0; ifile < MXFILE_BIASCOR; ifile++ ) 
    { INPUTS.simFile_biasCor[ifile] = (char*)malloc(MEMC); }

  // split item string
  splitString(item, comma, MXFILE_BIASCOR,    // inputs
	      &INPUTS.nfile_biasCor, INPUTS.simFile_biasCor ); // outputs 
  
  char *f0 = INPUTS.simFile_biasCor[0];
  if ( IGNOREFILE(f0) ) { INPUTS.nfile_biasCor = 0 ; }

  return;

} // parse_simfile_biasCor


// **************************************************
void parse_simfile_CCprior(char *item) {

  // Created May 16 2019
  // parse comma-separate list of CCprior files names.

  int ifile;
  int MEMC           = MXCHAR_FILENAME*sizeof(char);
  char comma[] = ",";
  //  char fnam[]  = "parse_simfile_CCprior" ;

  // ------------------ BEGIN -----------------

  // first allocate memory for file names
  INPUTS.simFile_CCprior = (char**)malloc( MXFILE_CCPRIOR*sizeof(char*));
  for(ifile=0; ifile < MXFILE_CCPRIOR; ifile++ ) 
    { INPUTS.simFile_CCprior[ifile] = (char*)malloc(MEMC); }

  // split item string
  splitString(item, comma, MXFILE_BIASCOR,    // inputs
	      &INPUTS.nfile_CCprior,
	      INPUTS.simFile_CCprior ); // outputs  

  char *f0 = INPUTS.simFile_CCprior[0];
  if ( IGNOREFILE(f0) ) { 
    INPUTS.nfile_CCprior          = 0 ;
    INPUTS.ipar[IPAR_scalePCC]    = 0 ; // make sure scalePCC is not floated
    return; 
  }

  // - - - - - - - - - - - 
  INFO_CCPRIOR.USE =  1;
  if ( strcmp(f0,"H11") == 0 ) { INFO_CCPRIOR.USEH11 =  1; }

  return ;

} // parse_simfile_CCprior

// **************************************************                                                                                  
void parse_cidFile_data(char *filename) {
  // Created Sep 5 2019 - Dillon Brout        

  int NCID    = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,filename);
  int MEMC    = NCID        * sizeof(char*);
  int MEMC2   = MXCHAR_CCID * sizeof(char);

  char *cid ;
  int i;
  int LDMP = 0 ;
  char fnam[]="parse_cidFile_data";

  // ------ BEGIN --------------


  // for an array of strings this is a 2d array
  INPUTS.cidList_data =  (char**)malloc(MEMC);

  for(i=0; i<NCID; i++ ) { 
    INPUTS.cidList_data[i] = (char*)malloc(MEMC2); 
    cid = INPUTS.cidList_data[i];
    get_PARSE_WORD(0,i,cid); // fill the list

    if ( strstr(cid,",") != NULL || 
	 strstr(cid,":") != NULL || 
	 strstr(cid,"=") != NULL ) {
      sprintf(c1err,"Invalid cid string = '%s'",cid);
      sprintf(c2err,"Check cid_select_file %s",filename);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    if (LDMP) {
      printf ("xxx %s: select cid = %s \n", fnam, INPUTS.cidList_data[i]);
    }
  }

  INPUTS.ncidList_data = NCID;
 
  return ;

} // END of parse_cidFile_data()



// **************************************************
void prep_input_nmax(char *item) {

  // Created July 17 2017
  // Parse comma-separated list of NMAX per survey.
  // Examples:
  //  nmax=200                ! limit on total, regardless of survey
  //  nmax=80(SDSS),200(DES)  ! survey-specific nmax
  //  nmax=300,200(DES)       ! nmax=300 for total, and 200 max for DES subset
  //

#define MXARG_nmax 20
  int  i, NARG, nmax, ID ;
  char stringArg[MXARG_nmax][MXCHAR_VARNAME];
  char *ptrArg[MXARG_nmax];
  char comma[] = "," ;
  char survey[60], tmpString[MXCHAR_VARNAME] ;
  char fnam[] = "prep_input_nmax" ;

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

} // end prep_input_nmax

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

  int  NARG, MXARG=3;
  char stringArg[2][20];
  char *ptrArg[2] = { stringArg[0], stringArg[1] } ;
  char comma[] = "," ;
  //  char fnam[] = "parse_powzbin" ;

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
  xxxxxxx */

  return ;

} // end parse_powzbin


// **************************************************
void parse_blindpar(char *item) {

  // Created Aug 2017
  // Parse 2 parameters separate by comma.
  // E.g., blindpar9=.15,343.2 --> blind offset is 0.15*cos(343.2)

  int  ipar=-9, NARG, MXARG=3;
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

  sscanf(ptrArg[0], "%le", &INPUTS.blind_cosinePar[ipar][0] ); 
  sscanf(ptrArg[1], "%le", &INPUTS.blind_cosinePar[ipar][1] ); 

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
  //  char  fnam[] = "parse_CUTWIN" ;

  // ---------- BEGIN ------------

  INPUTS.NCUTWIN++ ;
  ICUT = INPUTS.NCUTWIN-1 ;
  INPUTS.LCUTWIN_ABORTFLAG[ICUT]   = true ;   //  abort on missing var
  INPUTS.LCUTWIN_DATAONLY[ICUT]    = false ;  //  cut on data 
  INPUTS.LCUTWIN_BIASCORONLY[ICUT] = false ;  //  cut on sim data and biasCor

  sprintf(local_item,"%s", item);
  ptrtok = strtok(local_item," ");

  for ( i=0; i < 4 ; i++ ) {

    if ( i == 0 ) {
      // check for option in CUTWIN(option)
      sscanf ( ptrtok, "%s", KEY ); 
      extractStringOpt(KEY, stringOpt); // return stringOpt

      if ( strcmp(stringOpt,"NOABORT") == 0 ) 
	{ INPUTS.LCUTWIN_ABORTFLAG[ICUT] = false; } // allow missing variable 

      if ( strcmp(stringOpt,"DATAONLY") == 0 ) 
	{ INPUTS.LCUTWIN_DATAONLY[ICUT] = true ; } // cut on data only

      if ( strcmp(stringOpt,"BIASCORONLY") == 0 ) 
	{ INPUTS.LCUTWIN_BIASCORONLY[ICUT] = true ; } // cut on sim & biascor
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
	 ,INPUTS.LCUTWIN_ABORTFLAG[ICUT]
	 ) ;

} // end of parse_CUTWIN


// **************************************************
int apply_CUTWIN(int EVENT_TYPE, int *DOFLAG_CUTWIN, double *CUTVAL_LIST) {

  // Created Jan 2016
  // Returns TRUE if all CUTVAL_LIST values pass CUTWIN cuts.
  // Input EVENT_TYPE specifies DATA or BIASCOR
  //

  //  int IS_DATA    = ( EVENT_TYPE == EVENT_TYPE_DATA );
  int IS_BIASCOR = ( EVENT_TYPE == EVENT_TYPE_BIASCOR );
  int IS_CCPRIOR = ( EVENT_TYPE == EVENT_TYPE_CCPRIOR );

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
      printf(" xxx cut on %s = %f  (cutwin=%.3f to %.3f, EVENT_TYPE=%d)\n",
	     INPUTS.CUTWIN_NAME[icut], CUTVAL_LIST[icut],
	     INPUTS.CUTWIN_RANGE[icut][0], 
	     INPUTS.CUTWIN_RANGE[icut][1], EVENT_TYPE ); 
    }

    // check SIM-option to skip SIM_TYPE_INDEX
    cnam   = INPUTS.CUTWIN_NAME[icut] ;
    if( (IS_BIASCOR || IS_CCPRIOR) && usesim_CUTWIN(cnam)==0 ) 
      { continue ; }

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
  // isData=1 for datafile argument (real or sim dat);
  // isData=0 for biasCor sample.
  //
  // Oct 23 2018: 
  //  + new input arg isData=1 for data, zero for sim biasCor 
  //  + check DATAONLY flag.
  //
  // May 18 2020:
  //   + check BIASCORONLY flag.
  //

  bool  NOVAR       = ( ivar < 0 );
  bool  ABORTFLAG   = INPUTS.LCUTWIN_ABORTFLAG[icut];
  bool  DATAONLY    = INPUTS.LCUTWIN_DATAONLY[icut];
  bool  BIASCORONLY = INPUTS.LCUTWIN_BIASCORONLY[icut];
  int   ISTAT ;
  char *VARNAME   = INPUTS.CUTWIN_NAME[icut];
  char  fnam[] = "set_DOFLAG_CUTWIN" ;

  // ------------- BEGIN -------------
  
  if ( DATAONLY    && !isData ) { return(0) ; } // Oct 23 2018
  if ( BIASCORONLY &&  isData ) { return(0) ; } // May 18 2020

  if ( NOVAR && ABORTFLAG ) {
    sprintf(c1err,"Invalid CUTWIN on '%s' (ivar=%d, icut=%d, isData=%d)", 
	    VARNAME, ivar, icut, isData );
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
void parse_FIELDLIST(char *item) {

  // Created May 2020
  // break comma-separated list and load NFIELD values

  int  i ;
  int LDMP = 0 ;
  char comma[] = "," ;
  char fnam[] = "parse_FIELDLIST" ;

  // ------------ BEGIN ------------
  
  for(i=0; i < MXFIELD_OVERLAP; i++ ) 
    { INPUTS.FIELDLIST[i] = (char*) malloc(20*sizeof(char) ); }

  splitString(item, comma, MXFIELD_OVERLAP,               // inputs
	      &INPUTS.NFIELD, INPUTS.FIELDLIST ); // outputs
  
  if ( LDMP ) {
    for(i=0; i < INPUTS.NFIELD; i++ ) {
      printf(" xxx %s: select FIELD = '%s' \n", 
	     fnam, INPUTS.FIELDLIST[i] ); fflush(stdout);
    }
  }

  //  debugexit(fnam); // xxxx REMOVE

  return ;

} // end parse_FIELDLIST

// **************************************************
void SPLITRAN_prep_input(void) {

  // for NSPLITRAN>1, call this function before each new fit.

  //  char fnam[] = "SPLITRAN_prep_input" ;

  // ------------ BEGIN -----------

  FITINP.COVINT_PARAM_FIX   = INPUTS.sigmB ; 
  FITINP.COVINT_PARAM_LAST  = INPUTS.sigmB ; 
  INPUTS.covint_param_step1 = INPUTS.sigmb_step1 ; // default COVINT param

  if ( strlen(INPUTS.sigint_fix) > 0 ) {
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ; 
  }

  recalc_dataCov(); 

  return;

} // end SPLITRAN_prep_input


// **************************************************
void  SPLITRAN_cutmask(void) {

  // Created July 2012 by R.Kessler
  // set data[n].cutmask for SPLITRAN option.

  int NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;
  int n, snid, bit, CUTMASK ;
  char *name ;

  // ------------------ BEGIN ---------------

  if ( INPUTS.NSPLITRAN <= 1 ) { return ; }

  for (n=0; n < NSN_DATA; ++n)  {

    CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[n];

    // subtract SPLITRAN errbit if it was set earlier
    bit = CUTBIT_SPLITRAN;
    CUTMASK -= (CUTMASK & CUTMASK_LIST[bit]) ;

    // subtract ERRMASK_MINBIN if it was set earlier (May 2019)
    bit = CUTBIT_MINBIN;
    CUTMASK -= (CUTMASK & CUTMASK_LIST[bit] );

    INFO_DATA.TABLEVAR.CUTMASK[n] = CUTMASK; 
 
    // check random-subset option. Assume SN name is integer (i.e., sim)
    // to convert string *name to integer snid.
    // Not clear what happens for data names that are string.
    name = INFO_DATA.TABLEVAR.name[n];
    sscanf( name, "%d", &snid ); // problem with string?
    if ( SPLITRAN_ACCEPT(n,snid) == 0 ) 
      { setbit_CUTMASK(n, CUTBIT_SPLITRAN, &INFO_DATA.TABLEVAR );    }
  }

  return ;

} // end of SPLITRAN_cutmask


// **************************************************
int SPLITRAN_ACCEPT(int isn, int snid) {

  // Created july 2012 by R.Kessler
  // return 1 if this isn is accepted; 0 otherwise
  // Use integer SNID so that the same SN are used
  // for analyses with slightly different cuts.
  //
  
  int jj;

  if ( INPUTS.NSPLITRAN <= 1 ) { return 1; }

  jj = (snid % INPUTS.NSPLITRAN) + 1 ;

  if ( jj == NJOB_SPLITRAN ) 
    { return 1; }
  else if ( INPUTS.JOBID_SPLITRAN > INPUTS.NSPLITRAN ) 
    { return 1; } // for SPLITRAN summary only
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

  
  int ipar, isplit, NSPLIT, NPAR, NERR_VALID, JOBID_SPLIT ; 
  int FOUND_wfit = 1;
  double 
    VAL, ERR, XN, XNTMP, NSNAVG, NSNRMS
    ,SUMVAL1[MAXPAR], SUMERR1[MAXPAR]
    ,SUMVAL2[MAXPAR], SUMERR2[MAXPAR]
    ,AVG_VAL[MAXPAR], RMS_VAL[MAXPAR]
    ,AVG_ERR[MAXPAR], RMS_ERR[MAXPAR]
    ,SUMN, SQSUMN
    ;
  
  char PARNAME[MXCHAR_VARNAME];
  char fnam[] = "SPLITRAN_SUMMARY" ;

  // -------------- BEGIN -------------

  NSPLIT       = INPUTS.NSPLITRAN ;
  JOBID_SPLIT  = INPUTS.JOBID_SPLITRAN;
  NPAR         = FITINP.NFITPAR_ALL ;

  if ( NSPLIT <= 1 ) { return ; }

  // check for individual JOBID instead of entire set.
  if ( JOBID_SPLIT > 0 ) {
    if ( JOBID_SPLIT > NSPLIT ) {
      store_PARSE_WORDS(-1,"");
      for(isplit=1; isplit<=NSPLIT; isplit++ )
	{ SPLITRAN_read_fitpar(isplit); }

      // check for optional wfit results too (Jun 11 2020)
      // Note separate isplit loop to make sure that all jobs
      // have finished. If any wfit file is missing, STOP trying.
      for(isplit=1; isplit<=NSPLIT; isplit++ ) {
	if ( FOUND_wfit ) 
	  { FOUND_wfit = SPLITRAN_read_wfit(isplit);  }
      }

    }
    else
      { return; }
  }



  // get average sample size & RMS
  SUMN = SQSUMN = 0.0 ;
  for ( isplit=1; isplit <= NSPLIT; isplit++ ) {
    XNTMP = (double)FITRESULT.NSNFIT_SPLITRAN[isplit] ;
    SUMN   += XNTMP ;
    SQSUMN += (XNTMP*XNTMP) ;
  } 
  XN     = (double)NSPLIT ;
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

  fprintf(fp, "# SPLITRAN SUMMARY for %d jobs (Avg sample size: %d +- %d) \n",
	  INPUTS.NSPLITRAN, (int)NSNAVG, (int)NSNRMS );
  fprintf(fp," \n");
  fprintf(fp,"NSPLITRAN: %d \n", INPUTS.NSPLITRAN);
  fprintf(fp,"VARNAMES: ROW PARNAME AVG_VAL AVG_ERR RMS_VAL RMS_ERR \n" );

  for ( ipar=1; ipar <= NPAR ; ipar++ ) {
    if ( AVG_ERR[ipar] <= 0.0 ) { continue ; }
    sprintf(PARNAME, "%s", FITRESULT.PARNAME[ipar]);
    trim_blank_spaces(PARNAME);
    fprintf(fp, "ROW: %2d  %-15s  %8.3f  %8.3f  %8.3f %8.3f \n"
	   ,ipar
	   ,PARNAME
	   ,AVG_VAL[ipar]
	   ,AVG_ERR[ipar]
	   ,RMS_VAL[ipar]
	   ,RMS_ERR[ipar]
	   );
  }
  fclose(fp);

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
void prep_input_driver(void) {

  // May 9 2019: check INPUTS.fixpar_all

  int i,  NFITPAR, ifile, NTMP=0, USE_CCPRIOR, USE_CCPRIOR_H11 ;
  char usage[10];
  char *varname_pIa = INPUTS.varname_pIa;
  char fnam[] = "prep_input_driver";

  // ------------ BEGIN -----------

  if ( INPUTS.cat_only ) { prep_input_varname_missing(); return; }

  USE_CCPRIOR      = INFO_CCPRIOR.USE; 
  USE_CCPRIOR_H11  = INFO_CCPRIOR.USEH11; 

  ENVreplace("init",fnam,1); 

  // substitute ENV for filenames
  for(ifile=0; ifile < INPUTS.nfile_data; ifile++ ) 
    { ENVreplace(INPUTS.dataFile[ifile],fnam,1); }

  for(ifile=0; ifile < INPUTS.nfile_biasCor; ifile++ ) 
    { ENVreplace(INPUTS.simFile_biasCor[ifile],fnam,1); }

  for(ifile=0; ifile < INPUTS.nfile_CCprior; ifile++ ) 
    { ENVreplace(INPUTS.simFile_CCprior[ifile],fnam,1); }

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

  // make sure that only one method is specified for number of z bins.
  if ( INPUTS.nzbin    > 0 ) { NTMP++; }
  if ( INPUTS.nlogzbin > 0 ) { NTMP++; }
  if ( strlen(INPUTS.zbinuser) > 0 ) { NTMP++ ; }
  if ( NTMP != 1 ) {
    print_preAbort_banner(fnam);
    printf("\t nzbin=%d      \n",  INPUTS.nzbin);
    printf("\t nlogzbin=%d   \n",  INPUTS.nlogzbin);
    printf("\t zbinuser='%s' \n",  INPUTS.zbinuser);
    sprintf(c1err,"%d inputs specify number of redshift bins.", NTMP);
    sprintf(c2err,"Only one of these must be specified.");
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

  if ( INPUTS.nzbin >= MXz ) {
    sprintf(c1err,"nzbin=%d exceeds bound of MAXBIN_z=%d", 
	    INPUTS.nzbin, MXz);
    sprintf(c2err,"Check 'bins=' key in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  if (INPUTS.Nsntype <= 0 ) 
    { printf("SN selection type=ALL\n"); }
  else 
    { printf("SN selection sntype= %s\n", INPUTS.sntypeString); }


  printf("x1min/x1max = %.3f / %.3f \n", INPUTS.x1min, INPUTS.x1max);
  printf("cmin/cmax   = %.3f / %.3f \n", 
	 INPUTS.cmin, INPUTS.cmax);
  printf("logmass_min/max = %.3f/%.3f \n", 
	 INPUTS.logmass_min, INPUTS.logmass_max);
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

  int  ISFILE_BIASCOR = ( INPUTS.nfile_biasCor > 0 );
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
    sprintf(c1err,"opt_biascor=%d but simfile_biascor is not defined .",
	    INPUTS.opt_biasCor);
    sprintf(c2err,"Provide simfile_biascor, or set opt_biascor=0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  


  // if CCprior is not set then make sure to fix scalePCC and sigint
  if ( !USE_CCPRIOR  ) {
    if ( strlen(varname_pIa) > 0 ) {
      sprintf(c1err,"Illegal varname_pIa=%s without CC prior.", varname_pIa);
      sprintf(c2err,"Must set simfile_ccprior with varname_pIa");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
    }
    INPUTS.ipar[IPAR_scalePCC] = 0 ; // do NOT float scalePCC   
  }
  else {
    INPUTS.ipar[IPAR_COVINT_PARAM] = 0 ;  // do not float sigint
    INPUTS.parval[IPAR_COVINT_PARAM] = INPUTS.sigmB; 
    prep_input_probcc0();

    if ( INPUTS.ipar[IPAR_scalePCC] == 2 ) 
      { sprintf(FITPARNAMES_DEFAULT[IPAR_scalePCC], "scalePIa"); }
  } 

  // if scalePCC=0 and is fixed, turn off CC prior
  if( INPUTS.ipar[IPAR_scalePCC] == 0  && 
      fabs(INPUTS.parval[IPAR_scalePCC]) < 0.000001 )  { 

    if ( INPUTS.nfile_CCprior > 0 ) {
      printf("\t WARNING: scalePCC(u13,p13)=0,0 --> turn off CC prior\n");
      fflush(stdout);
    }

    INPUTS.nfile_CCprior = 0 ;
    INFO_CCPRIOR.USE = 0; 

  }


  if ( INPUTS.force_pIa >= 0.0 ) 
    { printf(" force_pIa = %.3f \n", INPUTS.force_pIa ); }

  if ( INPUTS.perfect_pIa ) 
    { printf(" force_pIa = 1 or 0 for true SNIa or SNCC \n"); }


  // if there is no user-defined selection of fit params,
  // then float everything. Otherwise stick with user choices.
  if ( USE_CCPRIOR_H11 ) {
    int ipar, IPAR;
    if ( !INPUTS.ipar[IPAR_scalePCC] )
      { INPUTS.ipar[IPAR_scalePCC] = 1; }

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
  prep_input_gamma();

  // sanity checks on fitting for GAMMA0 (mag step across host logmass)
  if ( INPUTS.USE_GAMMA0 ) {
    double TAU           = INPUTS.parval[IPAR_LOGMASS_TAU] ;

    if ( TAU < 0.001 ) {
      sprintf(c1err,"Invalid LOGMASS_TAU = %.4f for gammma0 fit", TAU);  
      sprintf(c2err,"logmasss_tau(p8) should be at least .01");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  prep_input_load_COSPAR();

  INPUTS.FLOAT_COSPAR=0;
  if ( INPUTS.ipar[IPAR_OL] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_Ok] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_w0] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_wa] ) { INPUTS.FLOAT_COSPAR=1; }

  prep_input_nmax(INPUTS.nmaxString);

  parse_cidFile_data(INPUTS.cidFile_data); // djb

  prep_input_varname_missing();

  // check for fast lookup option if all cosmo params are fixed.
  // arg, not needed.  prep_cosmodl_lookup();

  printf("\n"); fflush(stdout);

  INFO_DATA.TABLEVAR.EVENT_TYPE    = EVENT_TYPE_DATA ;
  INFO_BIASCOR.TABLEVAR.EVENT_TYPE = EVENT_TYPE_BIASCOR ;
  INFO_CCPRIOR.TABLEVAR.EVENT_TYPE = EVENT_TYPE_CCPRIOR ;

  return ;

} // end of prep_input_driver

// =====================================
void prep_input_load_COSPAR(void) {
  INPUTS.COSPAR[0] = INPUTS.parval[IPAR_OL] ;
  INPUTS.COSPAR[1] = INPUTS.parval[IPAR_Ok] ;
  INPUTS.COSPAR[2] = INPUTS.parval[IPAR_w0] ;
  INPUTS.COSPAR[3] = INPUTS.parval[IPAR_wa] ;
}

// **********************************************
void prep_input_probcc0(void) {

  // Created Sep 20 2019 by R.Kessler
  // If using CC prior, process input strings of the form
  //
  //   type_list_probcc0=3,33  
  //        (list of int TYPE values in data header)
  //
  //   idsurvey_list_probcc=DES,52,53 
  //       (list of string and/or ID from SURVEY.DEF file)
  //
  //  The type and idsurvey lists correspond to spec-confirmed
  //  SNIa, and thus PROBCC is forced explicitly to zero.
  //
  char *str_type_list     = INPUTS_PROBCC_ZERO.str_type_list ;
  char *str_idsurvey_list = INPUTS_PROBCC_ZERO.str_idsurvey_list ;
  int  LEN_type_list      = strlen(str_type_list);
  int  LEN_idsurvey_list  = strlen(str_idsurvey_list);
  int  DO_PROBCC0, i, itype, id, nval, NERR=0 ;
  int  NUSE, NUSE_IDSURVEY[MXIDSURVEY];
  char *str_values[MXPROBCC_ZERO], *surveyName ;
  char comma[] = "," ;
  char fnam[]  = "prep_input_probcc0" ;

  // ---------------- BEGIN ------------------
  
  DO_PROBCC0 = ( LEN_type_list > 0 || LEN_idsurvey_list > 0 ) ;
  if ( DO_PROBCC0 ) {
    INPUTS_PROBCC_ZERO.USE = true;
    for(i=0; i < MXPROBCC_ZERO; i++ ) 
      { str_values[i] = (char*)malloc( 80 * sizeof(char) ); }
  }
  else {
    return ;
  }


  // check TYPE from data header
  if ( LEN_type_list > 0 ) {
    splitString(str_type_list, comma, MXPROBCC_ZERO,    // inputs
		&nval, str_values ) ;                    // outputs    
    INPUTS_PROBCC_ZERO.ntype = nval ;
    for(i=0; i < nval; i++ ) {
      sscanf(str_values[i],"%d", &itype);
      INPUTS_PROBCC_ZERO.type_list[i] = itype ;
      printf("\t Force PROB_BBC(CC) = 0 for TYPE = %d\n", itype);
      fflush(stdout);
    }
  }
 

  // check survey ID from $SNDATA_ROOT/SURVEY.DEF
  if ( LEN_idsurvey_list > 0 ) {
    splitString(str_idsurvey_list, comma, MXPROBCC_ZERO,    // inputs
		&nval, str_values ) ;                    // outputs    
    INPUTS_PROBCC_ZERO.nidsurvey = nval ;

    for(id=0; id < MXIDSURVEY; id++ )  { NUSE_IDSURVEY[id] = 0; }

    for(i=0; i < nval; i++ ) {
      INPUTS_PROBCC_ZERO.idsurvey_list[i] = -9 ;

      id = NOINT  ;
      sscanf(str_values[i],"%d", &id);
      if ( id == NOINT ) { 
	// it's a survey name string, so convert back to ID
	id = get_IDSURVEY(str_values[i]);
      }

      if ( id > 0 ) {
	INPUTS_PROBCC_ZERO.idsurvey_list[i] = id ;
	surveyName = SURVEY_INFO.SURVEYDEF_LIST[id];
	NUSE_IDSURVEY[id]++ ;
	if ( IGNOREFILE(surveyName) ) 
	  { sprintf(surveyName,"UNKNOWN ID->ERROR"); NERR++; }
      }
      else {
	surveyName = str_values[i] ;
	strcat(surveyName," is UNKNOWN->ERROR");
	NERR++ ;
      }

      printf("\t Force PROB_BBC(CC) = 0 for IDSURVEY = %3d (%s)\n", 
	     id, surveyName );

      NUSE=NUSE_IDSURVEY[id];
      if( NUSE > 1 ) {
	printf("\t\t ERROR: IDSURVEY=%d used %d times\n", id, NUSE);
	NERR++ ;
      }

      fflush(stdout);
    }
  }

  
  // free string memory
  for(i=0; i < MXPROBCC_ZERO; i++ ) 
    { free(str_values[i]); }

  // abort on error
  if ( NERR > 0 ) {
    sprintf(c1err,"Found %d errors above.", NERR );
    sprintf(c2err,"check inputs "
	    "type_list_probcc0 and idsurvey_list_probcc0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  //  debugexit(fnam);

  return;

} // end prep_input_probcc0

int force_probcc0(int itype, int id) {

  // return force = 1 (true) if this type or id should be
  // forced to have probcc=0 (probIa=1).
  //
  // Inputs:
  //   itype = TYPE in data header
  //   id    = IDSURVEY (from SURVEY.DEF file)

  int  ntype     = INPUTS_PROBCC_ZERO.ntype ; 
  int  nidsurvey = INPUTS_PROBCC_ZERO.nidsurvey ; 
  int  force=0, i;
  //  char fnam[] = "force_probcc0";

  // ------------- BEGIN ---------------

  // check legacy option ...
  if ( itype == INPUTS.typeIa_ccprior ) { return(1); }

  // now check newer option.
  if ( !INPUTS_PROBCC_ZERO.USE ) { return(force); }  

  for(i=0; i < ntype; i++ )  { 
    if ( itype == INPUTS_PROBCC_ZERO.type_list[i] ) 
      { return(1) ; }
  }

  for(i=0; i < nidsurvey; i++ )  { 
    if ( id == INPUTS_PROBCC_ZERO.idsurvey_list[i] ) 
      { return(1) ; }
  }

  return(force);

} // end force_probcc0

// **********************************************
void  prep_input_gamma(void) {

  // Created Oct 2018
  // prepare stuff related to gamma0 & gamma1 for HR step.

  int  LDMP = 1;
  char *varname_gamma = INPUTS.varname_gamma ;
  char fnam[] = "prep_input_gamma" ;

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

} // end prep_input_gamma


// **********************************************
void  prep_input_varname_missing(void) {

  // Created May 8 2020
  // prep comma-sep list of varnames to allow missing in
  // one or more data files. Load INPUTS_VARNAME_MISSING struct.
  // Wildcards can also be used.
  // Example:
  //   append_varname_missing='PROB*,SNR*,FAKEID'
  //
  // will allow variables missing that contain PROB,
  // contain SNR, or have the exact variable FAKEID.
  //

  char *varname_missing = INPUTS.append_varname_missing ;
  char *varname_pIa     = INPUTS.varname_pIa ;
  char tmpName[100];
  char comma[] = ",", *ptrTmp ;
  int  MXVAR = MXVARNAME_MISSING ;
  int  MEMC  = 60*sizeof(char);
  int  ndef, i, LEN ; 
  bool wildcard;
  //  char fnam[] = "prep_input_varname_missing" ;

  // ----------- BEGIN ----------

  INPUTS_VARNAME_MISSING.ndef = 0 ;

  // if varname_pIa is specified, always add it to list
  if ( strlen(varname_pIa) > 0 ) {
    sprintf(tmpName, "%s", varname_missing);
    if ( strlen(varname_missing) == 0 ) 
      { sprintf(varname_missing,"%s", varname_pIa); return; }
    else
      { sprintf(varname_missing,"%s,%s", varname_pIa, tmpName); }
  }


  // - - - - 
  for(i=0; i < MXVAR; i++ ) 
    { INPUTS_VARNAME_MISSING.varname_list[i] = (char*)malloc(MEMC); }

  splitString(varname_missing, comma, MXVAR,         // inputs
	      &INPUTS_VARNAME_MISSING.ndef,             //output
	      INPUTS_VARNAME_MISSING.varname_list ) ;   //output

  ndef = INPUTS_VARNAME_MISSING.ndef ;
  if ( ndef == 0 ) { return; }

  for(i=0; i < ndef; i++ ) {
    ptrTmp = INPUTS_VARNAME_MISSING.varname_list[i];  LEN=strlen(ptrTmp);
    wildcard = false;
    if ( ptrTmp[LEN-1] == '*' ) { wildcard = true; ptrTmp[LEN-1] = 0; }
    INPUTS_VARNAME_MISSING.wildcard[i] = wildcard ;

    if ( !wildcard ) 
      { printf("\t Append missing column with varname = '%s' \n", ptrTmp ); }
    else
      { printf("\t Append missing column(s) with substring = '%s'\n",ptrTmp);}

  }

  return ;
} // end   prep_input_varname_missing


// **********************************************
void  prep_cosmodl_lookup(void) {

  // Nov 22 2017
  // If all cosmo params are fixed, then make binned
  // dl vs, z lookup for faster computation.

  int NBZ, iz ;
  double ZMIN, ZMAX, ZBIN, z, di, dl ;
  //  char fnam[] = "prep_cosmodl_lookup" ;

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
    // just need something approximate for initialization
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
  // Jun 2 2019: make sure idsample >= 0 for argument to get_COVINT

  int  NSN_DATA, n, l, m, idsample, cutmask ;
  double z, covtot, covint, covfit, COVMAT_INT[3][3];
  double SCALE_COVINT ;

  int USE_IDEAL_COVINT  = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  int USE_SIGINT_FIX    = ( strlen(INPUTS.sigint_fix) > 0 );
  double a = FITRESULT.ALPHA ;
  double b = FITRESULT.BETA ;
  double g = 0.0 ; 
  //  char fnam[] = "recalc_dataCov";

  // ------------ BEGIN --------------

  if ( USE_IDEAL_COVINT == 0  && USE_SIGINT_FIX==0 ) { 
    SCALE_COVINT = 1.0;     // traditional sigma_int model
  }
  else {    
    // Covint from biasCor or from sigint_fix(IDAMPLE)
    SCALE_COVINT = FITINP.COVINT_PARAM_FIX; 
  }

  NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ; 

  // - - - - - - - - - - 
  for (n=0; n< NSN_DATA; ++n)  {


    idsample = INFO_DATA.TABLEVAR.IDSAMPLE[n];
    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n];
    z        = INFO_DATA.TABLEVAR.zhd[n];

    if ( cutmask ) { continue ; }
    if ( idsample < 0 ) { idsample = 0; }  // Jun 2 2019

    // fetch COVINT depending on method
    if ( USE_IDEAL_COVINT ) {
      get_COVINT_biasCor(idsample,z,a,b,g,  COVMAT_INT); 
    }
    else {
      // traditional
      get_COVINT_model(idsample,COVMAT_INT); 
    }

    //Add intrinsic scatter to error matrix

    for (l=0; l<3; ++l) {
      for (m=0; m<3; ++m) {
	covfit = (double)INFO_DATA.TABLEVAR.covmat_fit[n][l][m];
	covint = COVMAT_INT[l][m] ;
	covtot = covfit + covint*SCALE_COVINT ;
	INFO_DATA.TABLEVAR.covmat_tot[n][l][m] = (float)covtot ;
      }
    }
    

  } // end of n-loop over data

  return ;

}  // end recalc_dataCov


// *******************************************************
double next_covFitPar(double redchi2, double parval_orig, double parval_step) {

  // Created Jan 26 2018
  // Returns next covFitPar; either sigmB or covScale.

  double parval_next;
  double step, slope=0.0, num=0.0, denom=0.0 ;
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

  if ( parval_next > 100. ) {
    print_preAbort_banner(fnam);
    if ( NFIT_ITER > 0 ) {
      printf("\t slope = %f/%f = %f \n", num, denom, slope);
      printf("\t SIGINT_LIST[iter=%d,%d] = %f, %f \n",
	     NFIT_ITER, NFIT_ITER-1,
	     FITRESULT.SIGINT_LIST[NFIT_ITER],
	     FITRESULT.SIGINT_LIST[NFIT_ITER-1] );
      printf("\t CHI2RED_LIST[iter=%d,%d] = %f, %f \n",
	     NFIT_ITER, NFIT_ITER-1,
	     FITRESULT.CHI2RED_LIST[NFIT_ITER],
	     FITRESULT.CHI2RED_LIST[NFIT_ITER-1] );
    }
    char *PARNAME = FITRESULT.PARNAME[IPAR_COVINT_PARAM];
    sprintf(c1err,"Crazy parval_next(%s) = %f  at NFIT_ITER=%d",
	    PARNAME, parval_next, NFIT_ITER);
    sprintf(c2err,"redchi2=%.2f, parval[orig,step]=%.3f,%.3f ", 
	    redchi2, parval_orig, parval_step);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
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
  // May 29 2019: call SPLITRAN_write_fitpar

  int  JOBID     = INPUTS.JOBID_SPLITRAN ;
  int  NSPLITRAN = INPUTS.NSPLITRAN ;
  char *prefix = INPUTS.PREFIX ;
  char tmpFile1[200], tmpFile2[200], tmpFile3[200];


  // --------------- BEGIN -------------

  if ( SNFUNPAR_CHI2INFO_OUTPUT.USE ) {
    SNFUNPAR_CHI2INFO_LOAD_OUTPUT();
    SNFUNPAR_CHI2INFO_WRITE();
    return ;
  }

  if ( strlen(prefix) > 0 && !IGNOREFILE(prefix)  ) {

    if ( INPUTS.NSPLITRAN == 1 )  { 
      sprintf(tmpFile1,"%s.fitres", prefix ); 
      sprintf(tmpFile2,"%s.M0DIF",  prefix ); 
    }
    else  { 
      sprintf(tmpFile1,"%s-SPLIT%3.3d.fitres", prefix, NJOB_SPLITRAN);
      sprintf(tmpFile2,"%s-SPLIT%3.3d.M0DIF",  prefix, NJOB_SPLITRAN);
      sprintf(tmpFile3,"%s-SPLIT%3.3d.fitpar", prefix, NJOB_SPLITRAN);
    }
    
    prep_blindVal_strings();
    write_fitres_driver(tmpFile1);  // write result for each SN
    write_M0(tmpFile2);      // write M0 vs. redshift


    // for single JOBID_SPLITRAN, write fitpar file so that they
    // can be scooped up later to make summary.
    if ( JOBID >=1 && JOBID <= NSPLITRAN ) 
      { SPLITRAN_write_fitpar(tmpFile3); }

  } 
  else {
    printf("\n PREFIX not specified --> no fitres output.\n");
    fflush(stdout);
  }

  return ;

} // end outFile_driver

// ******************************************
void write_version_info(FILE *fp) {

  fprintf(fp,"# SNANA_VERSION: %s \n", SNANA_VERSION_CURRENT);
  fprintf(fp,"# BBC_VERSION:   %d \n", BBC_VERSION);
  fprintf(fp,"\n");
  fflush(fp);

} // end write_version_info


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
  // Oct 13 2019: check bad bins to write MUDIFERR_ZERO[EMPTY]

  int iz, NVAR, irow, NFIT ;
  double z, zMIN, zMAX, VAL, ERR, dl, MUREF;
  char *tmpName, strval_OL[80], strval_w0[80];
  FILE *fp;
  char fnam[] = "write_M0" ;

  // ---------- BEGIN -----------

  if ( INPUTS.cutmask_write == -9 ) { return ; } // July 2016

  calc_zM0_data(); // fill FITRESULT.zM0[iz]

  fp = fopen(fileName,"wt");

  if ( !fp )  {
    sprintf(c1err,"Could not open M0-outFile");
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n Write MUDIF vs. redshift to %s\n" , fileName); 

  // - - - - 
  write_version_info(fp);

  fprintf(fp,"# Reference cosmology params for MUDIF: \n");

  if ( (INPUTS.blindFlag & BLINDMASK_FIXPAR)>0 ) {
    sprintf(strval_OL, "%s", INPUTS.blindString[IPAR_OL] );
    sprintf(strval_w0, "%s", INPUTS.blindString[IPAR_w0] );
  }
  else {
    sprintf(strval_OL, "%.3f", FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_OL]) ;
    sprintf(strval_w0, "%.3f", FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_w0]) ;
  }


  tmpName = FITRESULT.PARNAME[IPAR_OL];
  fprintf(fp,"#   REF %10s:  %s \n", tmpName, strval_OL );

  tmpName = FITRESULT.PARNAME[IPAR_w0];
  fprintf(fp,"#   REF %10s:  %s \n", tmpName, strval_w0 );

  tmpName = FITRESULT.PARNAME[IPAR_ALPHA0];
  fprintf(fp,"#   FIT %10s:  %.4f +- %0.4f \n", tmpName,
	  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_ALPHA0],
	  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_ALPHA0]  );

  tmpName = FITRESULT.PARNAME[IPAR_BETA0];
  fprintf(fp,"#   FIT %10s:  %.4f +- %0.4f \n", tmpName,
	  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_BETA0],
	  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_BETA0]  );

  fprintf(fp,"# \n");
  fflush(fp);

  // write blindFlag info (Mar 1 2017)
  if ( INPUTS.blindFlag > 0 && ISDATA_REAL ) 
    { write_blindFlag_message(fp); }

  write_NWARN(fp,0);
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
    NFIT = FITINP.NZBIN_FIT[iz] ;

    /* xxx mark delete Oct 13 2019 xxxxxxxxxxx
    if ( VAL == 0.0 && ERR == 0.0 )
      { ERR = MUDIFERR_EMPTY ; }

    if ( ERR == 0.0 )
      { ERR = MUDIFERR_ZERO ; }

    if ( FITINP.ISFLOAT_z[iz] == 0 )
      { VAL=0.0 ; ERR=MUDIFERR_EMPTY;  }
    xxxxxxx  */

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
void SPLITRAN_write_fitpar(char *fileName) {

  // May 29 2019
  // Write fit params to machine-parsable file so that they 
  // can be read back later for summary file.

  FILE *fout;
  int n, ISFLOAT, ISM0, iz ;
  double VAL, ERR ;
  char tmpName[60] ;
  char KEY[]  = "FITPAR:" ;
  char fnam[] = "SPLITRAN_write_fitpar";

  // ------------- BEGIN -------------

  printf(" %s: open %s \n", fnam, fileName);
  fout = fopen(fileName,"wt");

  // write NSNFIT
  sprintf(tmpName,"NSNFIT" );
  VAL = (double)FITRESULT.NSNFIT;  ERR=0.0;
  fprintf(fout,"%s  %-20s  %8.2f %8.2f \n",  KEY, tmpName, VAL, ERR );

  // write sigint
  sprintf(tmpName,"%s", FITRESULT.PARNAME[IPAR_COVINT_PARAM] );
  VAL = FITINP.COVINT_PARAM_FIX ;    ERR = 0.0 ;
  fprintf(fout,"%s  %-20s  %8.5f %8.5f \n",  KEY, tmpName, VAL, ERR );

  for ( n=0; n < FITINP.NFITPAR_ALL ; n++ ) {

    ISFLOAT = FITINP.ISFLOAT[n] ;
    ISM0    = n >= MXCOSPAR ; // it's z-binned M0

    // skip fixed cosmo params; keep fixed M0 params to 
    // print clear message about unused z-bins
    if ( ISM0 == 0  && ISFLOAT==0 ) { continue ; } 
      
    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;
    sprintf(tmpName, "%s", FITRESULT.PARNAME[n]);
    if ( ERR < 0.0 ) { continue ; }

    if ( n >= MXCOSPAR ) { 
      iz  = INPUTS.izpar[n] ;
      VAL = FITRESULT.M0DIF[iz]; 
      sprintf(tmpName,"%s-<M0avg>", FITRESULT.PARNAME[n] );
    }
    else {
      VAL += BLIND_OFFSET(n); // offset cosmo params besides M0
    }

    if ( !ISFLOAT ) { VAL = -9.0 ; ERR = -9.0; }

    fprintf(fout,"%s  %-20s  %8.5f %8.5f \n",
	    KEY, tmpName, VAL, ERR );
    

  } // end loop over fit params
  
  fclose(fout);
  return ;

} // end SPLITRAN_write_fitpar

// ******************************************
void SPLITRAN_read_fitpar(int isplit) {

  // May 29 2019
  // Read all of the fitpar files as if this were all in one job.

  FILE *fp;
  int  iwd, NWD, ipar ;
  int  NTRY_OPEN=0;
  double VAL, ERR;
  char tmpFile[200], LINE[100], WORD[6][MXCHAR_VARNAME], *PARNAME;
  char *prefix = INPUTS.PREFIX ;
  char fnam[] = "SPLITRAN_read_fitpar";

  // --------------- BEGIN --------------
  
  sprintf(tmpFile,"%s-SPLIT%3.3d.fitpar", prefix, isplit);
  printf(" Read %s \n", tmpFile); fflush(stdout);

 TRY_OPEN:
  fp = fopen(tmpFile,"rt");
  if( !fp ) {
    printf("\t fitpar file not created yet ... try again in 10 sec.\n"); 
    fflush(stdout);
    sleep(10);
    NTRY_OPEN++;
    if ( NTRY_OPEN < 10 ) {
      goto TRY_OPEN ;
    }      
    else {
      sprintf(c1err,"Could open SPLITRAN file for reading:");
      sprintf(c2err,"%s", tmpFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);           
    }
  }

  
  while(fgets(LINE,100,fp) != NULL ) {
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);    
    if ( NWD < 4 ) { continue; }
    iwd=0; get_PARSE_WORD(0,iwd,WORD[iwd]);
    if ( strcmp(WORD[iwd],"FITPAR:") != 0 ) { continue; }

    iwd=1; get_PARSE_WORD(0,iwd,WORD[1]);  // parName
    iwd=2; get_PARSE_WORD(0,iwd,WORD[2]);  // VAL
    iwd=3; get_PARSE_WORD(0,iwd,WORD[3]);  // ERR     
    PARNAME = WORD[1];
    sscanf(WORD[2],"%le", &VAL);
    sscanf(WORD[3],"%le", &ERR);   
    
    if ( strcmp(PARNAME,"NSNFIT") == 0 ) {
      FITRESULT.NSNFIT_SPLITRAN[isplit] = (int)VAL;
    }
    else if ( strcmp(PARNAME,"sigint") == 0 ) {
      FITRESULT.PARVAL[isplit][IPAR_COVINT_PARAM] = VAL ;
      FITRESULT.PARERR[isplit][IPAR_COVINT_PARAM] = 1.0E-8 ;
    }
    else {
      ipar = match_fitParName(PARNAME);
      FITRESULT.PARVAL[isplit][ipar] = VAL ;
      FITRESULT.PARERR[isplit][ipar] = ERR ;      
    }

    
  } // end while over lines in file

  fclose(fp);

  return;

} // end SPLITRAN_read_fitpar


// ******************************************
int SPLITRAN_read_wfit(int isplit) {

  // Created Jun 11 2020
  // Read optional wfit output (with .COSPAR extension).
  // Note that wfit is run by SALT2mu_fit.pl batch script,
  // and NOT run here by SALT2mu. This utility compiles
  // the wFit stats along with the SALT2mu stats.

  int  NWD, iwd, NLINE=0, EXIST = 0 ;
  FILE *fp;
  double w, werr;
  char *prefix = INPUTS.PREFIX ;
  char tmpFile[200], LINE[100], WORD[20] ;
  //  char fnam[] = "SPLITRAN_read_wfit" ;

  // ------------ BEGIN ---------------

  sprintf(tmpFile,"wfit_%s-SPLIT%3.3d.COSPAR", prefix, isplit);
  printf(" Read %s \n", tmpFile); fflush(stdout);

  fp = fopen(tmpFile,"rt");
  if( !fp ) {
    printf("\t no wfit output files. \n");
    return(EXIST);
  }
  
  // set w parameter as floated so that it gets printed;
  // also modify parName to indicate wfit origin.
  FITINP.ISFLOAT[IPAR_w0] = 1 ;  
  sprintf(FITRESULT.PARNAME[IPAR_w0],"w(wfit)" );

  EXIST = 1;

  while(fgets(LINE,100,fp) != NULL ) {
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);    
    if ( NLINE == 1 ) {
      iwd=0; get_PARSE_WORD(0,iwd,WORD);  sscanf(WORD, "%le", &w);
      iwd=1; get_PARSE_WORD(0,iwd,WORD);  sscanf(WORD, "%le", &werr);
      FITRESULT.PARVAL[isplit][IPAR_w0] = w ;
      FITRESULT.PARERR[isplit][IPAR_w0] = werr ; 
    }
    NLINE++ ;
  }
  
  fclose(fp);

  return(EXIST);

} // end SPLITRAN_read_wfit

// ************************************
int match_fitParName(char *parName) {
  int  ipar = -9;
  char *PARNAME;
  char fnam[] = "match_fitParName";

  // ------------- BEGIN -------------

  if ( strstr(parName,"m0_") != NULL ) { parName[5] = 0 ; }

  for(ipar=1; ipar < FITINP.NFITPAR_ALL; ipar++ ) {
    PARNAME = FITRESULT.PARNAME[ipar] ;
    trim_blank_spaces(PARNAME);
    if( strcmp(parName,PARNAME)==0 ) { return(ipar); }
  }

  sprintf(c1err,"Cannot find ipar index for parName='%s'", parName);
  sprintf(c2err,"Check valid parameter names.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

  return(-9);

} // end match_fitParName

// ******************************************
void write_fitres_driver(char* fileName) {

  // Write outout in fitres format.
  // Combine original fitres variables with those
  // calculated here in SALT2mu.
  // Jun 2012: increase line-string length from 400 to MXCHAR_LINE
  // Jul 2013: write averge M0 to fitres file (avemag0)
  // Apr 2016: for SIM, write true CC/tot ratio of fitted events.
  // Jul 2016: if cutmask_write=-9, do NOT write the SN lines
  // Jun 27 2017: REFACTOR z bins
  // Mar 01 2018: add M0 to output
  // Jun 10 2019: call printCOVINT 
  // Jan 09 2019: fix to loop over each datafile instead of only the first.
  // May 09 2020: 
  //   + refactor to write columns defined in OUTPUT_VARNAMES struct.
  //   + call define_varnames_append()
  //   + check cat_only option that skips fit and just catenates

  //  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  //  bool  IS_SIM            = INFO_DATA.TABLEVAR.IS_SIM ;
  int   IWD_KEY  = 0;
  int   IWD_CCID = 1;

  double VAL, ERR, PULL ;
  FILE  *fout, *finp;

  int n, ivar, indx, NCUT, icut, cutmask, NWR, ISFLOAT, iz, GZIPFLAG ;
  int idsample, NSN_DATA, NSN_BIASCOR, ifile ;
  char  line[MXCHAR_LINE], tmpName[60];
  char  ztxt[60], KEY[MXCHAR_VARNAME], CCID[40];

  char fnam[] = "write_fitres_driver" ;

  // ------------------ BEGIN ----------------

  NSN_DATA    =  INFO_DATA.TABLEVAR.NSN_ALL;
  NSN_BIASCOR =  INFO_BIASCOR.TABLEVAR.NSN_ALL;

  // define the new fitres variables to add to the original list  
  define_varnames_append();  // sets NVAR_APPEND and VARNAMES_APPEND

  // - - - - - - - - - -

  printf("\n Open output file with  cutmask_write=%d : \n", 
	 INPUTS.cutmask_write );
  printf("\t %s \n", fileName );  fflush(stdout);

  
  // - - - - - - - -  -
  fout = fopen(fileName,"wt");
  if (!fout ) {
    if ( INPUTS.cat_only ) {
      sprintf(c1err,"Could not open catfile_out='%s'", INPUTS.catfile_out);
      sprintf(c2err,"Check catfile_out key." );
    }
    else {
      sprintf(c1err,"Could not open output fitres file");
      sprintf(c2err,"'%s' ", fileName );
    }
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( INPUTS.cat_only ) {
    write_cat_info(fout); // write cat info to output file header
  }
  else {
    write_version_info(fout);
    if ( INPUTS.blindFlag > 0 && ISDATA_REAL ) 
      { write_blindFlag_message(fout); }      
    fprintf(fout,"# MU-RESIDUAL NOTE: MURES = MU-(MUMODEL+M0DIF) \n");
    write_NWARN(fout,1);
    write_MUERR_INCLUDE(fout);
  }

  if ( INPUTS.cutmask_write != -9 ) { 
    // write standard header keys
    fprintf(fout,"VARNAMES: ");
    
    for ( ivar=0; ivar < OUTPUT_VARNAMES.NVAR_TOT; ivar++ ) { 
      fprintf(fout,"%s ", OUTPUT_VARNAMES.LIST[ivar] );
    }
    

    // tack on SALT2mu/BBC variables (e.g., MU, MURES, etc ...)
    for ( ivar=0; ivar < NVAR_APPEND; ivar++ )   {  
      fprintf(fout,"%s ", VARNAMES_APPEND[ivar] );  
    }
    fprintf(fout, "\n\n");
  }
  else {
    // no SN lines, hence no header keys
  }

  if ( INPUTS.cat_only ) { goto WRITE_TABLE_ROWS; }

  NCUT = INPUTS.NCUTWIN ;
  if ( NCUT > 0 ) {
    fprintf(fout,"# CUTWIN Selection: \n");
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
    ISM0    = (n >= MXCOSPAR) ; // it's z-binned M0

    // skip fixed cosmo params; keep floated & fixed M0 params 
    // to print clear message about unused z-bins
    if ( ISM0 == 0  && ISFLOAT==0 ) { continue ; } 
      
    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;
    sprintf(tmpName,"%s", FITRESULT.PARNAME[n]);
    ztxt[0] = 0 ;
    if ( ERR < 0.0 ) { continue ; }

    if ( ISM0 ) { 
      iz    = INPUTS.izpar[n] ;
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

  if ( INPUTS.uave)
    { FITRESULT.SNMAG0 = FITRESULT.AVEMAG0; }
  else 
    { FITRESULT.SNMAG0 = INPUTS.nommag0; }


  fprintf(fout," \n");
  int NZwrite = 4; // include this many zM0 bins in COV dump
  printCOVMAT(fout, FITINP.NFITPAR_FLOAT, NZwrite);

  indx = 0;
  NWR  = 0;
  fflush(fout);

  // print contamination tables if CC prior is used
  if ( INFO_CCPRIOR.USE ) {
    print_contam_CCprior(fout, &CONTAM_MURES_BINS);
    print_contam_CCprior(fout, &CONTAM_REDSHIFT_BINS);
  }

  // check option to NOT write each SN to have smaller file
  // with only the fit results
  if ( INPUTS.cutmask_write == -9 ) { 
    fclose(fout);
    printf(" Nothing written out  (%d/%d used in fit) \n", 
	   FITRESULT.NSNFIT , NSN_DATA );  fflush(stdout);    
    return ;
  }

  // - - - - - - - -
  // re-read each data file

 WRITE_TABLE_ROWS:


  for(ifile=0; ifile < INPUTS.nfile_data; ifile++ ) {
    finp  = open_TEXTgz(INPUTS.dataFile[ifile], "rt", &GZIPFLAG); 

    while ( fgets (line, MXCHAR_LINE, finp) !=NULL  ) {

      if ( line[0] == ' '   ) { continue ; }
      if ( strlen(line) < 3 ) { continue ; }
      if ( commentchar(line) ) { continue; }

      store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,line);

      // skip if first wd of line is not a valid row key
      get_PARSE_WORD(0, IWD_KEY, KEY);
      if ( strcmp(KEY,"SN:") != 0 ) { continue ; }

      // check cutmask for writing events
      if ( !INPUTS.cat_only ) {
	get_PARSE_WORD(0, IWD_CCID, CCID);
	get_CCIDindx(CCID, &indx) ;      
	cutmask = INFO_DATA.TABLEVAR.CUTMASK[indx]; 
	if ( !keep_cutmask(cutmask)  ) { continue; }
      }


      NWR += write_fitres_line(indx,ifile,line,fout);

    }  // end reading line with fgets

    fclose(finp);   
  } // end ifile

  fclose(fout);

  if ( INPUTS.cat_only ) {
    printf(" Wrote %d SN to cat table.\n", NWR); 
  }
  else {
    printf(" Wrote %d SN  (%d/%d used in fit) \n", 
	   NWR, FITRESULT.NSNFIT , NSN_DATA );
  }
  fflush(stdout);

  return ;

} // end of write_fitres_driver

// ===============================================
int write_fitres_line(int indx, int ifile, char *line, FILE *fout) {

  // for input 'ifile' and 'line', write to fout.
  // Note that store_PARSE_WORDS has already been called,
  // so here use get_PARSE_WORDS to retrieve.
  //
  //  

  int NVAR_TOT = OUTPUT_VARNAMES.NVAR_TOT ;
  int ISTAT = 0 ;
  int  ivar_tot, ivar_file, ivar_word ;
  char word[MXCHAR_VARNAME], line_out[MXCHAR_LINE];
  char blank[] = " ";
  char fnam[] = "write_fitres_line" ;
  int  LDMP = 0 ;
  // ----------- BEGIN -----------

  sprintf(line_out,"SN: ");
  for(ivar_tot=0; ivar_tot < NVAR_TOT; ivar_tot++ ) {
    ivar_file = OUTPUT_VARNAMES.IVARMAP[ifile][ivar_tot];
    ivar_word = ivar_file + 1; // add 1 to skip SN: key
    if ( ivar_word >=0 ) 
      { get_PARSE_WORD(0, ivar_word, word); }
    else 
      { sprintf(word,"-9.0"); }

    if ( LDMP ) {
      printf(" xxx %s: ivar[tot,file,word]=%2d,%2d,%2d  word='%s' \n",
	     fnam, ivar_tot,ivar_file,ivar_word, word); fflush(stdout);
    }

    strcat(line_out,word);
    strcat(line_out,blank);
  }

  fprintf(fout,"%s", line_out);
  if ( !INPUTS.cat_only) { write_fitres_line_append(fout, indx); }

  fprintf(fout,"\n");

  ISTAT = 1;

  if ( LDMP ) { debugexit(fnam); }

  return(ISTAT);

} // end write_fitres_line


// ===============================================
void define_varnames_append(void) {

  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  int   NSN_BIASCOR       =  INFO_BIASCOR.TABLEVAR.NSN_ALL;
  char  tmpName[MXCHAR_VARNAME];
  //  char fnam[] = "define_varnames_append";

  // ----------- BEGIN -----------

  NVAR_APPEND = 0 ;
  if ( INPUTS.cat_only) { return; }

  sprintf(VARNAMES_APPEND[NVAR_APPEND],"CUTMASK");     NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MU");          NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUMODEL");     NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR");       NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR_RAW");   NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MURES");       NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUPULL");      NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"M0DIF");       NVAR_APPEND++ ;  

  if ( INFO_CCPRIOR.USE ) 
    { sprintf(tmpName,"CHI2_BEAMS"); }
  else                    
    { sprintf(tmpName,"CHI2"); }
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"%s",tmpName);   NVAR_APPEND++ ;  

  if ( INFO_CCPRIOR.USE ) 
    { sprintf(VARNAMES_APPEND[NVAR_APPEND],"PROBCC_BEAMS"); NVAR_APPEND++;  }


  if ( NSN_BIASCOR  > 0 ) {
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_mu");       NVAR_APPEND++ ;  
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCorErr_mu");    NVAR_APPEND++ ;  

    if ( DO_BIASCOR_MU == false ) {
      sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_mB");      NVAR_APPEND++ ;  
      sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_x1");      NVAR_APPEND++ ;  
      sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_c");       NVAR_APPEND++ ;  
    }
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasScale_muCOV");   NVAR_APPEND++ ;  
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"IDSAMPLE");          NVAR_APPEND++ ;  
  }

  return;
} // end define_varnames_append

// ===============================================
void write_NWARN(FILE *fp, int FLAG) {

  // FLAG=0 --> MUDIF file (1 row per z bin)
  // FLAG=1 --> FITRES file (1 row per SN)

  // ------------- BEGIN -------------

  if ( FLAG == 1 && INPUTS.cutmask_write != 0 ) {
    fprintf(fp,"# WARNING(MINOR): cutmask_write=%d "
	    "(grep CUTBIT SALT2mu.c | grep define) \n", INPUTS.cutmask_write );
    fflush(fp);
  }

  if ( NWARN_MUDIFERR_EMPTY > 0 ) {
    fprintf(fp,"# WARNING(MINOR): %d z bins excluded --> MUDIFFERR = %.0f\n", 
	    NWARN_MUDIFERR_EMPTY, MUDIFERR_EMPTY );
    fflush(fp);
  }


  if ( NWARN_MUDIFERR_ZERO > 0 ) {
    fprintf(fp,"# WARNING(SEVERE): %d z bins have MUDIFFERR = 0 --> %.0f\n", 
	    NWARN_MUDIFERR_ZERO, MUDIFERR_ZERO );
    fflush(fp);
  }

  // Apr 14 2020: check excessive loss from nobiascor cut
  int idsample, NREJ, NTOT; 
  double frac;
  char *NAME ;
  for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) {
    NTOT = NDATA_BIASCORCUT[0][idsample];
    NREJ = NDATA_BIASCORCUT[1][idsample];
    NAME = SAMPLE_BIASCOR[idsample].NAME ;
    frac = 0.0 ;
    if ( NTOT > 0 ) { frac = (double)NREJ / (double)NTOT; }

    if ( frac > INPUTS.frac_warn_nobiasCor ) {
      double frac_percent = 100.0*frac;
      fprintf(fp,"# WARNING(SEVERE): %d of %d events (%.1f %%) "
	      "have no biasCor for %s\n",
	      NREJ, NTOT, frac_percent, NAME );	   
    }
  }

  return ;

} // end write_NWARN

// ===============================================
void write_MUERR_INCLUDE(FILE *fp) {

  // Created Jun 25 2017
  // write vpec and lensing info in output file 
  // so that cosmology fitting program knows if/what
  // is already included in MUERR. Should help avoid
  // double-counting MUERR contributions.
  //
  // Oct 13 2019: write hash to make comment field

  int USE=0;
  double tmpErr;
  //  char fnam[] = "write_MUERR_INCLUDE";

  // --------------- BEGIN ----------------

  fprintf(fp, "# MUERR_INCLUDE: zERR \n" );  USE=1;

  tmpErr = INPUTS.zpecerr ;
  if ( tmpErr > 0.0 ) {
    fprintf(fp, "# MUERR_INCLUDE: zPECERR=%.4f \n", tmpErr );
    USE=1;
  }

  tmpErr = INPUTS.lensing_zpar;
  if ( tmpErr > 0.0 ) {
    fprintf(fp, "# MUERR_INCLUDE: SIGMA_LENS=%4f*z \n", tmpErr );
    USE=1;
  }

  if ( USE )  { fprintf(fp,"\n"); fflush(fp); }

  return ;

} // end write_MUERR_INCLUDE

// =====================================
void prep_blindVal_strings(void) {

  // Created Oct 22 2019
  // For each blinded param, create string of the form
  //   A + Bcos(C)
  //
  // so that the same blind-val can be written in
  // several different places.

  int ipar;
  double *blindpar, parval_orig; ;
  char *s;
  //  char fnam[] = "prep_blindVal_strings" ;

  // --------- BEGIN ---------

  if ( (INPUTS.blindFlag & BLINDMASK_FIXPAR) == 0 ) { return; }

  for(ipar=0; ipar < MAXPAR; ipar++ ) {
    blindpar = INPUTS.blind_cosinePar[ipar];
    if ( blindpar[0] == 0.0 ) { continue; }
    s           = INPUTS.blindString[ipar] ;
    parval_orig = INPUTS.parval[ipar] - blindpar[0]*cos(blindpar[1]);
    sprintf(s, "%8.4f + %6.4f*cos(%f)", 
	    parval_orig, blindpar[0], blindpar[1] );     
  }

  return;

} // end prep_blindVal_strings

// ===============================================
void write_blindFlag_message(FILE *fout) {
  int  blindFlag = INPUTS.blindFlag ;
  int  ipar;
  double *blindpar ;
  // ----------- BEGIN --------------

  if ( !blindFlag ) { return; }

  fprintf(fout,"#\n");
  fprintf(fout,"#  *** WARNING: RESULTS ARE BLINDED *** \n");

  if ( (blindFlag & BLINDMASK_MUz)>0 ) {
    fprintf(fout,"#  m0_## values have blind offset = cos( z * %.3f ) \n",
	    INPUTS.blind_cosinem0 ); 
  }
  else {
    for(ipar=0; ipar < MAXPAR; ipar++ ) {
      blindpar = INPUTS.blind_cosinePar[ipar];
      if ( blindpar[0] != 0.0 ) {
	fprintf(fout,"#     %s = %s \n",
		FITRESULT.PARNAME[ipar], INPUTS.blindString[ipar] );
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
  //
  // Sep 24 2019: write SUM_PROBIa and SUM_PROBCC
  //
  
  bool   IS_SIM = INFO_DATA.TABLEVAR.IS_SIM ;
  double chi2min, chi2red, SIGINT_AVG ;
  int    NDOF, NSNFIT, USE_CCPRIOR, NSN_BIASCOR ;
  BININFO_DEF  *BININFO_SIM_ALPHA, *BININFO_SIM_BETA ;

  // ------------- BEGIN -----------

  NSN_BIASCOR = INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS ;
  USE_CCPRIOR = INFO_CCPRIOR.USE ; 
  SIGINT_AVG  = INFO_BIASCOR.SIGINT_AVG;
  BININFO_SIM_ALPHA = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA  = &INFO_BIASCOR.BININFO_SIM_BETA ;
  

  if ( !USE_CCPRIOR  ) {
    chi2min = FITRESULT.CHI2SUM_1A; 
    chi2red = FITRESULT.CHI2RED_1A; 
    NDOF    = FITRESULT.NDOF;
    NSNFIT  = FITRESULT.NSNFIT ;
  }
  else {
    chi2min = FITRESULT.CHI2SUM_1A ;
    chi2red = FITRESULT.CHI2RED_1A ; 
    NDOF    = (int)FITRESULT.NSNFIT_1A - FITINP.NFITPAR_FLOAT ;
    NSNFIT  = FITRESULT.NSNFIT ;
  }

  int MASK  = INPUTS.opt_biasCor ;
  int NUMD =  INFO_BIASCOR.NDIM ;    // default number of biasCor dimensions
  if ( MASK > 0 ) {
    // xxx mark delete    if ( MASK & MASK_BIASCOR_1DZ ) { NUMD=1; } 
    fprintf(fout,"#  NSIM(%dD-BIASCOR)   = %d   "
	    "(N_alpha x N_beta = %d x %d) \n"
	    ,NUMD, NSN_BIASCOR
	    ,(*BININFO_SIM_ALPHA).nbin, (*BININFO_SIM_BETA).nbin    );

    fprintf(fout,"#  sigint(%dD-BIASCOR) = %.3f \n",
	    NUMD, SIGINT_AVG );
     
  }
    
  fprintf(fout,"#  NSNFIT        = %d \n", NSNFIT );

  if ( USE_CCPRIOR ) { // write out contamination info
    double xn1a = FITRESULT.NSNFIT_1A;
    double xncc = FITRESULT.NSNFIT_CC;
    double xncc_true = (double)FITRESULT.NSNFIT_TRUECC ;
    double xnsn      = (double)FITRESULT.NSNFIT;
    double contam = xncc/(xn1a+xncc);
    double contam_true = xncc_true/xnsn;
    fprintf(fout,"#  SUM_PROBIa    = %.2f \n", xn1a);
    fprintf(fout,"#  SUM_PROBCC    = %.2f \n", xncc);
    fprintf(fout,"#  CONTAM_DATA   = %.4f    "
	    "# SUM_PROB ratio\n", contam);
    if ( IS_SIM ) {
      fprintf(fout,"#  CONTAM_TRUE   = %.4f    "  
	      "# true NCC/(NIa+NCC) for sim-data\n", contam_true );
    }
    
    fflush(stdout);
  }

  fprintf(fout,"#  -2log(L)     = %.2f \n", FITRESULT.CHI2SUM_MIN );

  fprintf(fout,"#  chi2(Ia)/dof = %.2f/%i = %.3f  \n",
	  chi2min, NDOF, chi2red );
      

  fflush(fout);
  return ;

} // end write_fitres_misc

// ================================
void write_cat_info(FILE *fout) {

  int NFILE = INPUTS.nfile_data; 
  int ifile;
  // write cat info to fout file header.

  fprintf(fout,"# SNANA_VERSION: %s \n", SNANA_VERSION_CURRENT);

  fprintf(fout,"# Catenated data files: \n");
  for(ifile=0; ifile < NFILE; ifile++ ) {
    fprintf(fout,"#   + %s \n", INPUTS.dataFile[ifile] );
  }
  
  fprintf(fout,"# Appended columns: %s \n", 
	  INPUTS.append_varname_missing);
  fprintf(fout,"# Dropped columns: %s \n", 
	  OUTPUT_VARNAMES.DROPLIST );

  fprintf(fout,"#\n");
  fflush(fout);

  return;
} // end write_cat_info

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
  DATAFLAG = ( ISDATA_REAL || ( blindFlag & BLINDMASK_SIM)>0 );
  if ( DATAFLAG == 0 ) { return(NOTBLIND); }

  if ( INPUTS.blind_cosinePar[ipar][0] != 0.0 ) 
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
  // xxx mark delete  if ( FOUNDKEY_SIM )     { return zero ; }

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
  INPUTS.blind_cosinem0 = 10.0 ;
  arg = z * INPUTS.blind_cosinem0 ;
  off = cos(arg) ;
  return(off);

} // end BLIND_OFFSET


// ==========================================
void get_CCIDindx(char *CCID, int *indx) {

  // May 23, 2012
  // return integer *indx corresponding to CCID

  int n, NSN_DATA ;
  char fnam[] = "get_CCIDindx";

  // -------------- BEGIN -----------------

  n = *indx;

  while ( strcmp(INFO_DATA.TABLEVAR.name[n],CCID) != 0  ) {
    n++ ;  
    NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;
    if( n > NSN_DATA ) {
      sprintf(c1err,"Cannot find SALT2mu index for CCID='%s'", CCID);
      sprintf(c2err,"after checking all %d SN", NSN_DATA);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }


  *indx = n ; // load output arg

} // end of get_CCIDindx


// ============================================
void write_fitres_line_append(FILE *fp, int indx ) {

  // Oct 17, 2011 R.Kessler
  // append line in fitres file with SALT2mu info corresponding to 
  // this CID.
  // indx is the data index for INFO_DATA.TABLEVAR arrays.
  //
  // Mar 1 2018: add M0
  // Sep 26 2019: write probcc_beams
  // Apr 18 2020: write extra digit of precision for bias(mb,c,mu)
  // May 13 2020: write to char line, then single fprintf for entire line.

  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  // bool  IS_SIM            = INFO_DATA.TABLEVAR.IS_SIM ;
  double z, zerr, mu, muerr, muerr2, mumodel, mures, pull, M0DIF ;
  double sim_mb, sim_mu  ;
  double muBias=0.0, muBiasErr=0.0,  muCOVscale=0.0, chi2=0.0 ;
  double fitParBias[NLCPAR] = { 0.0, 0.0, 0.0 } ;
  int    n, cutmask, NWR, NSN_BIASCOR, idsample ;
  char line[400], word[40] ;	 
  char fnam[] = "write_fitres_line_append" ;

  // ------------ BEGIN --------------

  n =  indx;

  NSN_BIASCOR =  INFO_BIASCOR.TABLEVAR.NSN_ALL;
  z        = INFO_DATA.TABLEVAR.zhd[n] ;  
  zerr     = INFO_DATA.TABLEVAR.zhderr[n] ;
  mumodel  = INFO_DATA.mumodel[n];
  mu       = INFO_DATA.mu[n] - FITRESULT.SNMAG0; 
  muerr    = INFO_DATA.muerr[n];
  muerr2   = INFO_DATA.muerr_raw[n] ;
  mures    = INFO_DATA.mures[n] ;
  pull     = INFO_DATA.mupull[n] ;
  M0DIF    = INFO_DATA.M0[n] - FITRESULT.AVEMAG0 ;
  chi2     = INFO_DATA.chi2[n] ;
  cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n]  ;
  idsample = INFO_DATA.TABLEVAR.IDSAMPLE[n]  ;
  sim_mb   = INFO_DATA.TABLEVAR.SIM_FITPAR[INDEX_mB][n] ;
  sim_mu   = INFO_DATA.TABLEVAR.SIM_MU[n] ;

  if ( NSN_BIASCOR > 0 ) { 
    muBias               = INFO_DATA.muBias[n] ;
    muBiasErr            = INFO_DATA.muBiasErr[n] ;
    if ( DO_BIASCOR_MU == false ) {
      fitParBias[INDEX_mB] = INFO_DATA.fitParBias[n][INDEX_mB] ;
      fitParBias[INDEX_x1] = INFO_DATA.fitParBias[n][INDEX_x1] ;
      fitParBias[INDEX_c]  = INFO_DATA.fitParBias[n][INDEX_c] ;    
    }
    muCOVscale           = INFO_DATA.muCOVscale[n]  ;
  }
  
  if (pull > 99.999) { pull=99.999; }
  
  NWR=0;  line[0] = 0 ;
  sprintf(word, "%d ",    cutmask);   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", mu     );   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", mumodel);   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", muerr  );   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", muerr2 );   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", mures  );   NWR++ ; strcat(line,word);
  sprintf(word, "%6.3f ", pull   );   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", M0DIF  );   NWR++ ; strcat(line,word);
  sprintf(word, "%.2f ",  chi2   );   NWR++ ; strcat(line,word);

  if ( INFO_CCPRIOR.USE ) { 
    sprintf(word,"%.3e ", INFO_DATA.probcc_beams[n]);  NWR++; 
    strcat(line,word);
  }

  if ( NSN_BIASCOR > 0 ) {
    sprintf(word, "%6.4f ", muBias );       NWR++ ; strcat(line,word);
    sprintf(word, "%6.3f ", muBiasErr );    NWR++ ; strcat(line,word);
    if ( DO_BIASCOR_MU == false ) {
      sprintf(word,"%6.4f ", fitParBias[INDEX_mB]); NWR++; strcat(line,word);
      sprintf(word,"%6.3f ", fitParBias[INDEX_x1]); NWR++; strcat(line,word);
      sprintf(word,"%6.4f ", fitParBias[INDEX_c]);  NWR++; strcat(line,word);
    }
    sprintf(word, "%6.3f ", muCOVscale ) ;    NWR++ ; strcat(line,word);
    sprintf(word, "%d "   , idsample ) ;      NWR++ ; strcat(line,word);
  }

  
  fprintf(fp,"%s", line);
  fflush(fp);

  if ( NWR != NVAR_APPEND ) {
    sprintf(c1err,"Expected to write NVAR_APPEND=%d SALT2mu variables",
	    NVAR_APPEND);
    sprintf(c2err,"but wrote only %d variables.", NWR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;

} // end of write_fitres_line_append


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
  // Oct 13 2019: if MUDIFERR=0, set to 999 or 666 and increment NWARN

  int iz, n;
  double VAL, ERR;

  // --------------- BEGIN --------------
  NWARN_MUDIFERR_ZERO   = 0 ;
  NWARN_MUDIFERR_EMPTY  = 0 ;

  for ( n=MXCOSPAR ; n < FITINP.NFITPAR_ALL ; n++ ) {

    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;

    // one of the m0 params in a z bin
    VAL -= FITRESULT.AVEMAG0 ; 
    VAL += BLIND_OFFSET(n);

    iz  = INPUTS.izpar[n] ; 

    if ( FITINP.ISFLOAT_z[iz] == 0 )
      { VAL=0.0 ; ERR=MUDIFERR_EMPTY;  NWARN_MUDIFERR_EMPTY++; }
    else if ( ERR < 1.0E-12 )
      { ERR = MUDIFERR_ZERO ; NWARN_MUDIFERR_ZERO++; }

    FITRESULT.M0DIF[iz] = VAL ;
    FITRESULT.M0ERR[iz] = ERR ;

  }

} // end M0dif_calc

// *********************************************
void printCOVMAT(FILE *fp, int NPAR_FLOAT, int NPARz_write) {

  // Created Nov 30 2017
  //   [move code out of main]
  //
  // Inputs:
  //   fp          : file pointer to write to
  //   NPAR_FLOAT  : total number of floated params
  //   NPARz_write : number of zM0 params to include
  //
  // Jun 2019: pass fp to allow writing to file or stdout.
  // Dec 9 2019: avoid truncating FITRESULTS.PARNAME to 10 chars.

  int num, iMN, jMN, i, j, NZwr ;
  double emat[MAXPAR][MAXPAR], terr[MAXPAR], corr ;
  char tmpName[MXCHAR_VARNAME], msg[100];

  // --------- BEGIN ----------

  if ( NPARz_write > MXCOSPAR ) 
    { sprintf(msg, "Reduced COV matrix:"); }
  else
    { sprintf(msg, "Reduced COV matrix, includes %d m0_ bins:", 
	      NPARz_write); }

  fprintf(fp, "\n# %s\n", msg); fflush(fp);
  num = MAXPAR;
  mnemat_(emat,&num);
  NZwr = 0 ;

  for ( iMN=0 ; iMN < NPAR_FLOAT ; iMN++ )
    {
      i  = FITINP.IPARMAP_MN[iMN]  ;
      if ( i - MXCOSPAR >= NPARz_write ) { continue; }
      sprintf(tmpName, "%s", FITRESULT.PARNAME[i] );
      tmpName[10] = 0;    // force truncation of name
      terr[iMN] = sqrt(emat[iMN][iMN]);

      //      fprintf(fp,"%2i ",iMN+1); 
      fprintf(fp,"# %-10.10s ", tmpName ); 

      for (jMN=0; jMN <= iMN; ++jMN)	{
	j  = FITINP.IPARMAP_MN[jMN]  ;
	if ( j - MXCOSPAR >= NPARz_write ) { continue; }
	corr = emat[iMN][jMN]/(terr[iMN]*terr[jMN]);
	fprintf(fp, "%6.3f ",corr);
      }
      fprintf(fp,"\n");  fflush(fp);
    }

  fprintf(fp,"#\n\n");  fflush(fp);

  return ;

} // end printCOVMAT

// *********************************************
int keep_cutmask(int cutmask) {

  // Returns 1 if this cutmask is allowed and the
  // SN is written to the output fitres file.
  // Returns 0 if this error code is not allowed;
  // the SN is not written out.
  //
  // User input 'cutmask_write' is a mask of
  // errocode bits that are allowed; if any other
  // bit is set then return 0.
  //

  int ovp;
  int cutmask_write = INPUTS.cutmask_write ;
  //  char fnam[] = "keep_cutmask" ;

  // ---------------- BEGIN -------------

  // -1 is the same as allowing all bits.
  if (cutmask_write == -1 )                { return 1 ; }
  if (cutmask_write == 0 && cutmask == 0 ) { return 1 ; }

  // check specific cut bits to keep
  ovp = ( INPUTS.cutmask_write & cutmask );
  
  if ( (cutmask - ovp) > 0 ) 
    { return 0 ; }
  else
    { return 1 ; }

} // end of keep_cutmask


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
  double dflat, distance, H0inv ;
  double omega_m, omega_l, omega_k, OK, wde, wa ;
  //  char fnam[] = "cosmodl" ;

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


// =======================================================
//  SNPAR_CHI2INFO functions to use SALT2mu as 
//  function to return chi2-info for higher-level
//  MCMC fitter for SN population parameters.
// =======================================================

void  SNFUNPAR_CHI2INFO_INIT(void) {

  int   nbc, nbx1, nbm, ix1, ic, im;
  BININFO_DEF *BININFO ;
  char fnam[] = "SNFUNPAR_CHI2INFO_INIT" ;

  // ------------- BEGIN ----------

  if ( !SNFUNPAR_CHI2INFO_OUTPUT.USE ) { return; }

  print_banner(fnam);

  // - - - - - - - - - - - - - - - - - - - 
  // open output file
  char *outFile = SNFUNPAR_CHI2INFO_OUTPUT.OUTFILE ;
  FILE *fp;
  fp = fopen(outFile, "wt");
  if ( !fp ) {
    sprintf(c1err,"Could not open SNPAR outfile='%s'", outFile);
    c2err[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  SNFUNPAR_CHI2INFO_OUTPUT.FP_OUT = fp;
  printf("\t Open outFile: %s \n", outFile);
  fflush(stdout);


  // - - - - - - - - - - - - - - - - - - - 
  // load BININFO for c,x1,m
  // Hard-code for now, maybe later take user input

  BININFO = &SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_x1];
  nbx1 = SNFUNPAR_CHI2INFO_LOAD_BININFO("x1", -3.0, 3.0, 2.0, BININFO, fp);

  BININFO = &SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_c];
  nbc = SNFUNPAR_CHI2INFO_LOAD_BININFO("c", -0.3, 0.5, 0.2, BININFO, fp);
  
  BININFO = &SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_m];
  nbm = SNFUNPAR_CHI2INFO_LOAD_BININFO("m", 5.0, 15.0, 5.0, BININFO, fp);

  printf("\t Setup CHI2INFO bins: nbx1=%d nbc=%d nbm=%d \n",
	 nbx1, nbc, nbm );

  // - - - - - - - - - - - - - - - - - - - 
  // prepare 3D -> 1D map
  int MEMx1 = nbx1 * sizeof(int**) ;
  int MEMc  = nbc  * sizeof(int*) ;
  int MEMm  = nbm  * sizeof(int) ;
  int NBIN_TOT = 0;
  SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D = (int***) malloc(MEMx1);
  for(ix1=0; ix1 < nbx1; ix1++ ) {
    SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D[ix1]= (int**) malloc(MEMc);
    for(ic=0; ic < nbc; ic++ ) {
      SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D[ix1][ic]= (int*) malloc(MEMm);
      for(im=0; im < nbm; im++ ) {
	SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D[ix1][ic][im] = NBIN_TOT ;
	NBIN_TOT++ ;
      }
    }
  }


  fprintf(fp,"NBIN_TOT: %d \n", NBIN_TOT);

  int MEMI = NBIN_TOT * sizeof(int) ;
  int MEMD = NBIN_TOT * sizeof(double) ;
  SNFUNPAR_CHI2INFO_OUTPUT.NBIN_TOT = NBIN_TOT ;
  SNFUNPAR_CHI2INFO_OUTPUT.NEVT = (int*) malloc(MEMI);
  SNFUNPAR_CHI2INFO_OUTPUT.MURES_SUM   = (double*) malloc(MEMD);
  SNFUNPAR_CHI2INFO_OUTPUT.MURES_SQSUM = (double*) malloc(MEMD);

  // store VARNAMES line
  sprintf(SNFUNPAR_CHI2INFO_OUTPUT.LINE_VARNAMES,
	  "VARNAMES: ROW       ix1 ic im   NEVT  MURES_SUM   MURES_SQSUM");

  //  debugexit(fnam);
  return ;

} // end SNFUNPAR_CHI2INFO_INIT


// ===================================================
int SNFUNPAR_CHI2INFO_LOAD_BININFO(char *varName, 
				   double xmin, double xmax, double binSize,
				   BININFO_DEF *BININFO, FILE *fp ) {
  int nbin;
  double x;
  //  char fnam[] = "SNFUNPAR_CHI2INFO_LOAD_BININFO" ;
  // ------------ BEGIN -----------
  nbin = 0 ;
  sprintf(BININFO->varName,"%s", varName );
  for(x = xmin ; x < xmax; x += binSize ) {
    BININFO->lo[nbin]  = x ;
    BININFO->hi[nbin]  = x + binSize ;
    BININFO->avg[nbin] = x + 0.5*binSize ;
    nbin++ ;
  }
  BININFO->nbin = nbin ;
  BININFO->binSize = binSize ;

  fprintf(fp,"BININFO: %2s %2d %6.3f %6.3f\n",   
	  varName, nbin, xmin, xmax);

  return(nbin) ;
} // end SNFUNPAR_CHI2INFO_LOAD_BININFO


// ===========================================
void SNFUNPAR_CHI2INFO_LOAD_OUTPUT(void) {

  // called after each fit, load output struct.

  int NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL ;
  int NBIN_TOT      = SNFUNPAR_CHI2INFO_OUTPUT.NBIN_TOT ;
  int i, isn, cutmask, ix1,ic,im, isnpar ;
  int i3d[3], J1D ;
  double xval, mures ;
  char *CCID;
  BININFO_DEF *BININFO ;
  char fnam[] = "SNFUNPAR_CHI2INFO_LOAD_OUTPUT";

  // ---------- BEGIN ----------
  for(i=0; i < NBIN_TOT; i++ ) {
    SNFUNPAR_CHI2INFO_OUTPUT.NEVT[i]        = 0 ;
    SNFUNPAR_CHI2INFO_OUTPUT.MURES_SUM[i]   = 0.0 ;
    SNFUNPAR_CHI2INFO_OUTPUT.MURES_SQSUM[i] = 0.0 ;
  }

  // loop over data

  for(isn=0; isn < NSN_DATA; isn++ ) {

    CCID    = INFO_DATA.TABLEVAR.name[isn]; 
    cutmask = INFO_DATA.TABLEVAR.CUTMASK[isn]; 
    if ( !keep_cutmask(cutmask)  ) { continue; }
    ix1 = ic = im = -9;

    for(isnpar = 0; isnpar < 3; isnpar++ ) {

      if ( isnpar == ISNPAR_x1 ) 
	{ xval  = INFO_DATA.TABLEVAR.fitpar[INDEX_x1][isn]; }
      else if ( isnpar == ISNPAR_c ) 
	{ xval  = INFO_DATA.TABLEVAR.fitpar[INDEX_c][isn]; }
      else if ( isnpar == ISNPAR_m ) 
	{ xval  = INFO_DATA.TABLEVAR.logmass[isn]; }
      
      BININFO = &SNFUNPAR_CHI2INFO_OUTPUT.binInfo[isnpar];
      i3d[isnpar] = IBINFUN(xval, BININFO, 1, fnam);

      if ( i3d[isnpar] < 0 ) {
	sprintf(c1err,"Invalid i3d[%d] = %d for xval=%f", 
		isnpar, i3d[isnpar], xval );
	sprintf(c2err,"Check CCID = '%s'", CCID);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    ix1 = i3d[ISNPAR_x1];  ic=i3d[ISNPAR_c]; im=i3d[ISNPAR_m];
    J1D = SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D[ix1][ic][im];
    SNFUNPAR_CHI2INFO_OUTPUT.NEVT[J1D]++ ;

    mures    = INFO_DATA.mures[isn] ;
    SNFUNPAR_CHI2INFO_OUTPUT.MURES_SUM[J1D]   +=  mures;
    SNFUNPAR_CHI2INFO_OUTPUT.MURES_SQSUM[J1D] +=  (mures*mures);
  } // end isn 


  return ;

} // end  SNFUNPAR_CHI2INFO_LOAD_OUTPUT

// ===========================================
void SNFUNPAR_CHI2INFO_WRITE(void) {

  //  int NBIN_TOT = SNFUNPAR_CHI2INFO_OUTPUT.NBIN_TOT ;
  int nbx1     = SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_x1].nbin;
  int nbc      = SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_c].nbin;
  int nbm      = SNFUNPAR_CHI2INFO_OUTPUT.binInfo[ISNPAR_m].nbin ;
  FILE *FP     = SNFUNPAR_CHI2INFO_OUTPUT.FP_OUT ;

  int ITER=1;
  char NAME[40], tmpName[40];
  int ix1, ic, im, IBIN1D, NEVT, n, ISFLOAT, ISM0 ;
  double SUM, SQSUM, VAL, ERR ;
  //  char fnam[] = "SNFUNPAR_CHI2INFO_WRITE" ;

  // ----------- BEGIN -------------

  fprintf(FP,"\n# ================================================ \n");
  // write fitted nuisance params 
  for ( n=0; n < FITINP.NFITPAR_ALL ; n++ ) {

    ISFLOAT = FITINP.ISFLOAT[n] ;
    ISM0    = (n >= MXCOSPAR) ; // it's z-binned M0

    if ( ISFLOAT && !ISM0 ) {
      VAL = FITRESULT.PARVAL[1][n] ;
      ERR = FITRESULT.PARERR[1][n] ;
      sprintf(tmpName,"%s", FITRESULT.PARNAME[n]);
      fprintf(FP,"FITPAR:  %-14s = %10.5f +- %8.5f \n",
	      tmpName, VAL, ERR );
    }
  } // end loop over SALT2mu fit params

  // - - - - - - 

  fprintf(FP,"%s\n", SNFUNPAR_CHI2INFO_OUTPUT.LINE_VARNAMES);

  for(ix1=0; ix1 < nbx1; ix1++ ) {
    for(ic=0; ic < nbc; ic++ ) {
      for(im=0; im < nbm; im++ ) {
	IBIN1D = SNFUNPAR_CHI2INFO_OUTPUT.MAP3D_to_1D[ix1][ic][im];
	NEVT   = SNFUNPAR_CHI2INFO_OUTPUT.NEVT[IBIN1D];
	SUM    = SNFUNPAR_CHI2INFO_OUTPUT.MURES_SUM[IBIN1D];
	SQSUM  = SNFUNPAR_CHI2INFO_OUTPUT.MURES_SQSUM[IBIN1D];

	sprintf(NAME,"ITER%5.5d-%4.4d", ITER, IBIN1D);
	fprintf(FP,"ROW: %s %2d %2d %2d  %4d %14.6le %14.6le \n", 
		NAME, ix1,ic,im,  NEVT, SUM, SQSUM); 
      }
    }    
  }

  return ;

} // end SNFUNPAR_CHI2INFO_WRITE
