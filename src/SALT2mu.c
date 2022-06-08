/*******************************************
Created by J. Marriner.
Installed into snana v8_38, January 2010. 

Program to take output from the SALT fitter dict files and
1.  Determine alpha and beta parameters
2.  Output a file of bias-corrected distances for cosmological fits

==USAGE 

 To dump a list of all valid keys and then quit,
    ./SALT2mu.exe KEY_DUMP

./SALT2mu.exe  <parameter input file>
    file=file.fitres bins=10 zmin=0.02 zmax=1.02 u0=1 u1=1 u2=0 u3=0
    prefix=SALT2mu sigmB=0.0 sigx1=0.0 sigc=.1  p9=.73 

The 'datafile=' argument (or file=) is required, and it is usually the 
FITRES table output from snlc_fit.exe of submit_batch_jobs.sh.

Additional arguments below are optional:

file=<comma-sep list of fitres file names to analyze>
datafile=<same as file=>

datafile_override=over1.dat,over2.dat,etc... ! comma-sep list of data-overrides
     ! enabled for zHEL, VPEC, VPEC_ERR, HOST_LOGMASS
     ! if VPEC [VPEC_ERR] is changed, so is zhd [zhderr]

minos=1  ! default; minos errors
minos=0  ! switch to MIGRAD (should be faster)
     
nmax=100                 ! fit first 100 events only
nmax=70(SDSS),200(PS1MD) ! fit 70 SDSS and 200 PS1MD
nmax=300,200(PS1MD)      ! fit 300 total, with 200 in PS1MD sub-sample

cid_select_file=<file with CID 'accept-only' list  
      (FITRES or unkeyed list format)>
cid_reject_file=<file with CID reject list 
      (FITRES or unkeyed format)>

izbin_from_cid_file=1 ! use izbin in cid_selecr_file

uzsim=1                  ! cheat and use true zCMB for redshift

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

opt_biascor=1136  ! +1024 -> if no valid biasCor, keep event with bias=0
                  !   (default is to discard events with no biasCor)

sigint_biascor=<sigmb_bias>   ! instead of auto-computed sigint_biascor
snrmin_sigint_biascor         ! SNRMIN to compute siginit_biascor 
prescale_biascor=<subset>,<preScale>  ! select <subset> from <prescale>

fieldGroup_biascor='shallow,medium,deep'             ! for biasCor & CCprior
fieldGroup_biascor='C3+X3,X1+E1+S1,C2,X2+E2+S2+C2'

surveygroup_biascor='CFA3+CSP,PS1MD'   ! combine CFA3+CSP into one biasCor
surveygroup_biascor='CFA3+CSP(zbin=.02),PS1MD' 
surveygroup_biascor='CFA3+CSP(zbin=.02:cbin=.04:x1bin=.4),PS1MD' 
surveygroup_biascor='CFA3+CSP(zbin=.02),SDSS(zbin=.04),PS1MD' 
surveygroup_biascor_abortflag=1  ! 0->allow survey(s) that are not in data

  NOTE: if OPT_PHOTOZ column exists in the input FITRES tables, 
        then each biasCor group  is automatically split into 
        [GROUPNAME]-zSPEC and [GROUPNAME]-zPHOT groups. 
        OPT_PHOTOZ is from &SNCLINP input SNTABLE_LIST='FITRES NOZPHOT'
zspec_errmax_idsample=0.002  ! default=0
   ! IS_SPECZ = OPT_PHOTOZ==0 || zhelerr < zspec_errmax_idsample
   ! Thus if all events have OPT_PHOTOZ=2, user input
   ! zspec_errmax_idsample  defines zSpec IDSAMPLE

idsample_select=2+3                ! fit only IDSAMPLE = 2 and 3
surveylist_nobiascor='HST,LOWZ'    ! no biasCor for these surveys
interp_biascor_logmass=1           ! allows turning OFF logmass interpolation

select_trueIa=1          ! select only true SNIa, disable CC prior
force_realdata=1         ! treat sim like real data

ndump_nobiascor=20       ! dump for first 20 data events with no biasCor
dumpflag_nobiascor=20;   ! idem; legacy input variable

frac_warn_nobiascor=0.02  ! print warning in output fitres file if nobiascor
                  ! cut-loss exceeds this fraction (applies to each IDSAMPLE)

cidlist_debug_biascor  ! comma-sep list to dump biasCor info

To check sample stats for each surveyGroup and fieldGroup,
   grep IDSAMPLE  <stdout_file>

 - - - - -  CCprior options - - - - -  
simfile_ccprior=name    ! sim fitres file to compute CC prior and
                        !  flag to to a BEAMS-like fit
simfile_ccprior=same  --> use same file(s) as simfile_bias
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
idsurvey_list_probcc0=5,50,51,53 
       or
idsurvey_list_probcc0=CSP,JRK07,KAIT,CFA3
       or
idsurvey_list_probcc0=CSP,JRK07,51,53
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
nzbin    = number of z bins to use (beware: some bins may be empty)
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
CUTWIN(NOABORT)  SIM_x1  -2 2       ! do not abort if SIM_x1 isn't there
CUTWIN(DATAONLY)    LOGMASS  5 12   ! cut on data only (not on biasCor)
CUTWIN(BIASCORONLY) BLA_LOGMASS  5 12 ! cut on biasCor (not on data)
CUTWIN varname_pIa  0.9 1.0   ! substitute argument of varname_pIa, and 
                              ! easy to change varname_pIa on command line

CUTWIN(FITWGT0) varname_pIa  0.9 1.0   ! MUERR->888 instead of cut

CUTWIN NONE  ! command-line override to disable all cuts;
             ! e.g., useful with cid_select_file

#select field(s) for data and biasCor with
fieldlist=X1,X2   # X1 and X2
fieldlist=X3      # X3 only
fieldlist=X       # any field with X in name

chi2max=16             # default cut to all events
chi2max(FITWGT0)=16    # no cut; instead set fit wgt=0 with large MUERR
chi2max(DES,PS1)=12    # apply cut only to DES & PS1
chi2max(CSP)=10        # apply cut to CSP
     = chi2-outlier cut applied before fit, using initial values. Beware 
       that initial and final chi2 can differ, so allow slop in the cut.
  Cut applied to -2log10(ProbIa_BEAMS + ProbCC_BEAMS); see Eq 6 of BBC paper.
  Note that multiple chi2max inputs are allowed. In the avove example, 
  chi2max=12 is applied to DES & PS1 events; chi2max=10 is applied to CSP; 
  chi2max=16 is for all other surveys (e.g., SDSS, other lowz, etc...)


cutmask_write=[allowed errcode mask to write output fitres]
cutmask_write=-1  -> write everything
cutmask_write=0   -> write only selected SN used in fit (default)
cutmask_write=64  -> include SNe that fail CUTWIN 
cutmask_write=-9  -> write comments only with fit results; no SN lines

write_yaml=1 -> write yaml info output for batch script
write_csv=1  -> write M0DIF vs. z in csv format for CosmoMC input

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

u1 = 1  --> if true, vary parameter 1 (1=True)
u2 = 1  --> if true, vary parameter 2 (1=True)
u3 = 1  --> if true, vary parameter 3 (1=True)
.
.
u13 = 1 or 2 (see p13 comments above)
u15=1 --> alphaHost = dalpha/dlogmass
u16=1 --> betaHost  = dbeta /dlogmass
u15=2 --> alphaHost = alpha shift about logmass_cen
u16=2 --> betaHost  = beta  shift about logmass_cen

u1 = 3 --> float alpha, but force 1 biasCor bin
u2 = 3 --> float beta,  but force 1 biasCor bin
u5 = 3 --> float gamma, but force 1 biasCor bin

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

zwin_vpec_check=0.01,0.05 [default] -> compute RMS(HR) for 0.01<z<0.05 
     using zHD and again with vpec sign-flip; 
     abort if sign-flip RMS(HR) is smaller than no flip.
zwin_vpec_check=0,0 -> disable check on sign-flip.


lensing_zpar --> add  "z*lensing_zpar" to sigma_int

fitflag_sigmb or sig1fit = 1 -->  find sigmB giving chi2/N = 1
fitflag_sigmb or sig1fit = 2 -->  idem, with extra fit adding 2log(sigma)
redchi2_tol   or sig1tol = tolerance (0.02) on chi2/dof-1

prescale_simdata=<preScale>  # pre scale for sim data
prescale_simcc=<preScale>    # pre-scale only the simulated CC
prescale_simIa=<preScale>    # pre-scale only the simulated Ia
nthread=<n>                  # use pthread for multiple cores on same node

NSPLITRAN=[NRAN] = number of independent sub-samples to run SALT2mu.
                  A separate fitres file is created for each sub-sample.
JOBID_SPLITRAN=[JOBID] 
   --> do only this splitran job, JOBID=1,2 ... NSPLITRAN
   --> write summary file if JOBID > NSPLITRAN

iflag_duplicate=1  # 0=ignore, 1=abort, 2=merge

snid_mucovdump='5944'  # after each fit iteration, full muCOV dump 

debug_mucovscale=44  # print info for j1d=44 (mucovscale cell), and also
                    # write biasCor-fitres file with mucovScale info

Default output files (can change names with "prefix" argument)
  SALT2mu.log
  SALT2mu.fitres
  SALT2mu.aux

utility to catenate FITRES files with different columns:
If PROB is in a subset of files, can force keeping PROB column
with append_varname_missing,
  SALT2mu.exe cat_only \
              datafile=<commaSepList> \
              append_varname_missing='PROB*' \  # keep PROB columuns
              catfile_out={cat_file_out} 

 restore_mucovscale_bug=1 -> restore bug of keeping cells where muCOVscale
                             is not defined.

==============HISTORY========================

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
              instead of parab errors. With low-stat samples,              errors seemed too small.

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
   + for NSPLITRAN option, write summary to [prefix]_SPLITRAN_summary.out,
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
 
 Jun 17 2020:
   New input surveygroup_biascor_abortflag=0 will NOT abort if  
   surveygroup_biascor includes survey(s) not in the data.

 Jun 24 2020: 
    minor refactor for NSPLITRAN to remove SPLIT[nnn] from file names
    since the files are written to [nnnn]/ sub-directories.
    NSPLITRAN probably doesn't work interactively; only via batch.
     
 Jun 25 2020:
   + when using MU in biasCor (opt_biadcor += 128), fix mu_true
     to be based on measured redshift instead of true redshift.

  Jul 1 2020: abort if NEVT(biasCor)=0 for any IDSAMPLE

  Jul 2 2020
    + add SIM_MUz[ISN] array to better track calc MU at observed zHD.

  Jul 3 2020
    + new option u1=3, u2=3, u3=3 -> force one biasCor bin for 
      continuous distribution.
    + SUBPROCESS refactoring so that SALT2mu can be called by python driver.

 Jul 17 2020: new option to dump all keys and then quit
     SALT2mu.exe KEY_DUMP

  Jul 24 2020: skip reading YAML keys for re-written batch script.

  Jul 29 2020
    + new input option minos=0 (switch to MIGRAD errors)
    + add CPU timer per fit, and write CPU to output fitres file

 Aug 12 2020:
    + new input flag write_yaml=1 to produce YAML file for
      batch-submit script. 
 
 Sep 2 2020:
   + refactor print_contam_CCprior to be computed after the fit,
     instead of during, so that it doesn't need threading logic.

   + new nthread=<n> argument breaks up chi2 calc into threads
     using pthread_create. Default nthread=1 does not use pthread.

 Sep 23 2020: finish subprocess OUTPUT_TABLEs

 Sep 24 2020: 
    + new key cid_reject_file
    + cid_reject_file and cid_select_file can be either FITRES format
      or list without any keys or commas.

 Sep 29 2020: few tweaks to get_fitParBias

 Oct 12 2020: fix rare numerical artifact for COV in merge_duplicates
 Oct 14 2020
    + allow comma-sep list for cid_select_file and cid_reject_file
    + refactor parsing comma-sep lists using new sntools utility
       parse_commaSep_list(...)

 Oct 28 2020:
   +  new input zwin_vpec_check
   +  modify set_DOFLAG_CUTWIN to allow CUTWIN on non-existing 
      PROB* variable for BiasCor sample.
   +  allow "CUTWIN varname_pIa xxx yyy", and make substitution internally

 Nov 12 2020:
   +  dumpflag_nobiascor is now number of events to dump.
      Add ndump_nobiascor key to do same thing.
   +  begin datafile_override feature.
 
 Nov 20 2020: few fixes to avoid conflicts with cat_only option.

 Dec 01 2020: 
    + call diagnostic function M0dif_rebin_check
    + add IZBIN & M0ERR to output FITRES file -->
         enable re-binning for makeCOV
    + write Nz x Nz cov matrix to [prefix].cov

 Dec 4 2020; 
    allow passing output fitres file as input (recycling);
    the original appended variables are excluded.

 Dec 6, 2020: write muerr_renorm to output such that sum of WGT
              in each z-bin matches 1/M0DIFERR^2.

 Dec 10 2020: remove legacy SPLITAN-summary functions; no longer
              needed with submit_batch_jobs script. Interactive
              NSPLITRAN probably won't work any more.

 Dec 11 2020: 
   + read zhel and zhelerr
   + cosmodl: split z arg -> zhel and zcmb.

 Dec 21 2020:   
   + opt_biaskcor += 1024 -> keep event if no valid biasCor; set biasCor=0.

 Jan 27 2021
   + fix few issues related to header_override
   + MXSTORE_DUPLICATE -> 200 (was 20)

 May 02 2021: new input zspec_maxerr_idsample
 May 12 2021: move read_data_override call before set_CUTMASK call.
 May 24 2021: disable cuts with "CUTWIN NONE"
 May 25 2021: new debug_malloc=1 input
 Jun 02 2021: replace all errmsg with new errlog that includes
               FP_STDOUT argument.
 Jun 07 2021:
    + write full command
    + new arg SUBPROCESS_STDOUT_CLOBBER=0 to turn off default clobber
    + fix major subprocess bug; reset cut bits for REWGT and PRESCALE
    + debug_flag=68 --> muCOVscale is NOT applied to vpec part of MUERR
    + use hash table for cid_select_file & cid_reject_file.

 Jun 15 2021:
    + release change of NOT applying MUCOVSCALE to vpec part of MUERR.
      debug_flag=-68 to go back.

 Jun 16 2021: 
    + remove LEGACY table functions for SUBPROCESS
    + new input args force_realdata and select_trueIa

 Jun 29 2021: 
   + if interp_biascor_logmass=0, implement for muCOVscale
        as well as biasCor.
   + if no biasCor, disable malloc for lots of mubias arrays;
      needed to save memory for large samples used by SUBPROCESS.

 Aug 02 2021:
   + Dillon added sigint in bins determined from biascor (motivated by BS20).
       Enabled opt_biascor+=4096
       Beware of notable code refactors.

 Sep 14 2021:
   + fix subtle bug by rejecting cells where MUCOVSCALE is not defined.
   + CUTWIN varname matching previous varname replaces it, enabling
     looser cuts on command line.

 Sep 27 2021:  
   + require muCOVadd>0 to use it (function MNCHI2FUN)
   + define MINPERCELL_MUCOVSCALE 5

 Nov 24 2021: new input key prescale_simIa
 Dec 28 2021: few tweaks so that prescale works with HOSTLIB

 Jan 22 2022
    + fix SUBPROCESS bug reading ref sim-input file
    + integrate REFAC_SUBPROC_STD to be default (no more debug_flag=930) 
    + increase DROPLIST array size to avoid overwrite bug

 Feb 26 2022 RK
    + tweak zspec_errmax_idsample to operate on zhelerr rather than zhderr
      because the latter includes VPECERR.

 Mar 02 2022 RK 
    + read zPRIOR[ERR] from data files if zspec_errmax_idsample>0
      so that zSPEC vs zPHOT is determined from originak zSPEC err
      and not from fitted zSPEC err.
    + new IS_SPECZ_TABLEVAR function to determine if event is specz or photoz.
    + improve debug_mucovscale to write two files: info vs. bin and 
        info vs. biasCor event.

 Mar 04 2022 :  see new DUST_FLAG

 Mar 07 2022 M.Vincenzi and R.Kessler
    + set default INPUTS.cmin/cmax/x1min/x1max to nominal SALT2 cuts,
      and add input nbinc_mucovscale with default=3 bins.
      
 Mar 17 2022: few fixes for sigint_fix 

 Mar 25 2022:
    + fix bug and read input sigint_step1
    + new input dchi2red_dsigint to specify slope so that fewer
      fit-iterations are needed.

 Apr 18 2022: new input izbin_from_cid_file=1 (same as debug_flag=401)
 Apr 22 2022: default minos=0 (was 1) and sigint_step1=0.01 (was .05)
               --> faster fitting

 ******************************************************/

#include "sntools.h" 
#include "sntools_cosmology.h" 
#include "sntools_output.h" 
#include <gsl/gsl_fit.h>  // Jun 13 2016
#include <sys/types.h>
#include <sys/stat.h>

#define USE_THREAD   // Sep 2020 : used in SUBPROCESS mode

#ifdef USE_THREAD
#include <pthread.h>
#endif

// ==============================================
// define data types to track selection cuts

//#define BBC_VERSION  2
//#define BBC_VERSION  3   // Jul 3 2020: add SUBPROCESS functions
#define BBC_VERSION  4     // Sep 2020: add pthread option
#define MXTHREAD    20

#define EVENT_TYPE_DATA     1
#define EVENT_TYPE_BIASCOR  2
#define EVENT_TYPE_CCPRIOR  3
#define MXEVENT_TYPE        4
char STRING_EVENT_TYPE[MXEVENT_TYPE][12] = 
  { "NULL", "DATA", "BIASCOR", "CCPRIOR" } ;

char STRING_MINUIT_ERROR[2][8] = { "MIGRAD", "MINOS" };

#define MAXPAR_MINUIT 110  // limit for MINUIT params (Sep 10 2017)

#define MXFILE_DATA     20  // max number of data files to read
#define MXFILE_BIASCOR  20  // max number of biasCor files to read
#define MXFILE_CCPRIOR  20  // max number of CCprior files to read

#define FLAG_EXEC_STOP   1
#define FLAG_EXEC_REPEAT 2

#define MXVAR_OVERRIDE 10

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
#define MINPERCELL_BIASCOR    10  // min per cell for biasCor
#define MINPERCELL_MUCOVSCALE  5  // min per cell for MUCOVSCALE

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
#define MASK_BIASCOR_MUCOVSCALE   32     // bit5: apply MUCOVSCALE vs. z
#define MASK_BIASCOR_SAMPLE  64     // bit6: biasCor vs. IDSAMPLE
#define MASK_BIASCOR_MU     128     // bit7: bias on MU instead of mB,x1,c
#define MASK_BIASCOR_COVINT 256     // bit8: biasCor-covint(3x3) x SCALE
#define MASK_BIASCOR_SIGINT_SAMPLE 512  // bit9: sigint(biasCor) vs. IDSAMPLE
#define MASK_BIASCOR_noCUT  1024    // do NOT reject event if no biasCor
#define MASK_BIASCOR_DEFAULT  MASK_BIASCOR_5D + MASK_BIASCOR_MUCOVSCALE + MASK_BIASCOR_SAMPLE
#define MASK_BIASCOR_MAD    2048 // bit11: median instead of rms in COV scale
#define MASK_BIASCOR_MUCOVADD 4096 // bit12: use error floor for error tuning


#define USEMASK_BIASCOR_COVFIT 1
#define USEMASK_BIASCOR_COVINT 2
#define USEMASK_BIASCOR_COVTOT 3
#define USEMASK_BIASCOR_ZMUERR 4 // include zMUERR(VPECERR) 

#define IFLAG_DUPLICATE_IGNORE 0
#define IFLAG_DUPLICATE_ABORT  1
#define IFLAG_DUPLICATE_AVG    2  // use weighted avg of SALT2 fit par.
#define MXSTORE_DUPLICATE    800  // abort if more than this many

#define MUERR_FITWGT0  8888.8  // MUERR-> this value in fit for FITWGT0 option
#define STRING_FITWGT0 "FITWGT0"

#define MUDIFERR_EMPTY 999.0   // if no events in bin, set error to 999
#define MUDIFERR_ZERO  666.0   // if MUDIFFERR=0, set to 666
int NWARN_MUDIFERR_ZERO ;
int NWARN_MUDIFERR_EMPTY ;


FILE *FP_STDOUT ;  // direct stdout to screen (stdout) or log file
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

#define DOFLAG_CUTWIN_IGNORE  0 // do not apply CUTWIN
#define DOFLAG_CUTWIN_APPLY   1 // apply CUTWIN
#define DOFLAG_CUTWIN_FITWGT0 2 // do not apply CUTWIN; deweight instead 

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

// define FITRES var names for misc variables
#define VARNAME_VPEC       "VPEC"
#define VARNAME_VPECERR    "VPECERR" 
#define VARNAME_VPECERR2   "VPEC_ERR"
#define VARNAME_zHD        "zHD"
#define VARNAME_zHDERR     "zHDERR"
#define VARNAME_zHEL       "zHEL"
#define VARNAME_zHELERR    "zHELERR"
#define VARNAME_zCMB       "zCMB"
#define VARNAME_LOGMASS    "HOST_LOGMASS"
#define VARNAME_LOGSFR     "HOST_LOGSFR"
#define VARNAME_LOGsSFR    "HOST_LOGsSFR"
#define VARNAME_COLOR      "HOST_COLOR"

#define VARNAME_IZBIN    "IZBIN"

// ---------------------
double LOGTEN  ;

// define CUTBIT and CUTMASK variables to determine events to
// reject, and to track statistics.
int  NSTORE_CUTBIT[MXEVENT_TYPE][MXCUTBIT] ;
int  CUTMASK_LIST[MXCUTBIT];          // ERRMASK for each bit
char CUTSTRING_LIST[MXCUTBIT][40];    // text definition per cut
int  *CUTMASK_POINTER[MXEVENT_TYPE];
short int  *IDSAMPLE_POINTER[MXEVENT_TYPE];
int  *NALL_CUTMASK_POINTER[MXEVENT_TYPE];
int  *NPASS_CUTMASK_POINTER[MXEVENT_TYPE];
int  *NREJECT_CUTMASK_POINTER[MXEVENT_TYPE];
int  NDATA_BIASCORCUT[2][MXNUM_SAMPLE]; // track biascor reject for data
int  NPASS_CUTMASK_BYSAMPLE[MXEVENT_TYPE][MXNUM_SAMPLE]; // 9.18.2021

int MAXSN ;
int NJOB_SPLITRAN; // number of random split-jobs

int    NSIMDATA, NSIMIa, NSIMCC ;  // to implement prescale_sim[Data,Ia,CC]
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
  float  *zhd,     *zcmb,    *zhel,    *zprior,    *vpec ;
  float  *zhderr,  *zcmberr, *zhelerr, *zpriorerr, *vpecerr, *zmuerr  ;
  float  *host_logmass, *host_logsfr, *host_logssfr, *host_color ;
  float  *pIa, *snrmax ;
  short int  *IDSURVEY, *SNTYPE, *OPT_PHOTOZ ; 
  bool   *IS_PHOTOZ;
  float  *fitpar_ideal[NLCPAR+1], *x0_ideal, *peakmjd ;

  // SIM_xxx read from table file
  short int *SIM_NONIA_INDEX ;
  float *SIM_ALPHA, *SIM_BETA, *SIM_GAMMADM;
  float *SIM_X0, *SIM_FITPAR[NLCPAR+1]; 
  float *SIM_ZCMB, *SIM_VPEC, *SIM_MU;
  float *SIM_MUz ; // calculated SIM_MU at observed redshift

  // scalar flags & counters computed from above
  bool   IS_DATA, IS_SIM ;
  int    DOFLAG_CUTWIN[MXCUTWIN]; // flag to apply each cutwin
  int    ICUTWIN_GAMMA;           // CUTWIN index for gamma-variable
  int    ICUTWIN_VARNAME_PIA ;      // CUTWIN index for varname_pIa
  int    IVAR_VPEC, IVAR_SIM_VPEC, IVAR_OPT_PHOTOZ ; // logicals
  int    IVAR_ZPRIOR;
  int    IVAR_SNTYPE, IVAR_SIM_GAMMADM;
  int    IVAR_pIa[MXFILE_DATA]; // track this variable for each file 

  int    NSN_PER_SURVEY[MXIDSURVEY] ;
  float  zMIN_PER_SURVEY[MXIDSURVEY] ;
  float  zMAX_PER_SURVEY[MXIDSURVEY] ;

  // quantities determined from table var
  float  *CUTVAL[MXCUTWIN];  // used only to make cuts.
  short int *IZBIN, *IDSAMPLE;
  int       *CUTMASK, *IMUCOV ;
  float     *mumodel;

  // covariance matrix(mB,x1,c)
  bool   *warnCov;
  float  ***covmat_fit ; // covmat from LCfit (no intrinsic scatter)
  float  ***covmat_tot ; // fit + intrinsic scatter matrix
  int     NCOVFIX;
  
} TABLEVAR_DEF ;


typedef struct { double VAL[NLCPAR][NLCPAR]; } COV_DEF ;

typedef struct {
  // parameters used for bias correction
  double z, host_logmass, alpha, beta, gammadm, FITPAR[NLCPAR+1];
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

  // define arrays which depend on SURVEYGROUP/FIELDGROUP sub-sample
  bool    *USE;                 // T -> enough events in cell to use
  int     *NperCell ;           // Number of events in each cell, vs. J1D
  double  *AVG_z;               // wgt-avg z in each cell, vs J1D
  double  *AVG_m;
  double  *AVG_LCFIT[NLCPAR];   // idem for mB,x1,c, vs. J1D

  double   **ABSPULL; //used only for MUCOVSCALE MAD option
  double   **MURES; //used only for MUCOVSCALE MAD option
  double   **MUCOV; //used only for MUCOVSCALE MAD option

} CELLINFO_DEF ;




typedef struct {

  // WARNING: update copy_IDSAMPLE_biasCor for each new item

  // group of fields: e.g, '82N'
  char NAME_FIELDGROUP[MXCHAR_SAMPLE] ; 

  // group of surveys: e.g.,CSP+CFA
  char NAME_SURVEYGROUP[MXCHAR_SAMPLE]; 

  // string-option for each sample (user-defined binning)
  char STRINGOPT[MXCHAR_SAMPLE];

  bool IS_PHOTOZ ;

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

  double *mumodel, *M0, *mu, *muerr, *muerr_renorm, *muerr_raw, *muerr_vpec;
  double *mures, *mupull;
  double *muerr_last, *muerrsq_last, *sigCC_last, *sqsigCC_last ;
  double *muCOVscale, *muCOVadd, *muBias, *muBiasErr, *muBias_zinterp; 
  double *chi2, *probcc_beams;
  double **fitParBias;

  bool *set_fitwgt0; // flag to set fit wgt=0 with MUERR -> large value

  // before fit, store bias[isn][ialpha][ibeta][igammadm]
  FITPARBIAS_DEF ****FITPARBIAS_ALPHABETA ; 
 
  // before fit, store muCOVscale and muCOVadd at each alpha,beta,gamma bin
  double ****MUCOVSCALE_ALPHABETA ; 
  double ****MUCOVADD_ALPHABETA ;
  short int ****I1D_MUCOVSCALE ;

  float MEMORY;  // Mbytes

  // Nov 2020: add data-override info
  int  NVAR_OVERRIDE;
  char **VARNAMES_OVERRIDE ;           // list of override varnames
  float *PTRVAL_OVERRIDE[MXVAR_OVERRIDE]; // pointers to override values
  int  *IVAR_OUTPUT_INVMAP;  // map ivar_out to ivar_over

  bool USE_IZBIN_from_CIDFILE ;
  int *IZBIN_from_CIDFILE; // Apri 2022
  int  NCHANGE_IZBIN;  // Number of events with IZBIN change

} INFO_DATA;


struct {

  int ILCPAR_MIN, ILCPAR_MAX ; // either mB,x1,c  OR just mu
  TABLEVAR_DEF TABLEVAR ;

  int  NDIM;     // number of dimensions for biasCor e.g., 1, 5, 7
  char STRING_PARLIST[20]; // e.g., z,x1,c,a,b
  bool DUST_FLAG ;  // True for dust-based model using nbin_beta=1

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
  float **MUCOVADD ;

  // wgted-avg redshift  corresponding to each M0 offset (for M0 outfile)
  double zM0[MXz];
  
  float MEMORY;  // Mbytes

  // diagnostic arrays if "write biascor" option is set
  float  *muerr_raw, *muerr, *mu, *mures ;
  int    *J1D_COV;

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
CELLINFO_DEF *CELLINFO_MUCOVADD ;


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
} INPUTS_PROBCC_ZERO ;


#define MXVARNAME_MISSING 10
struct {
  int    ndef ; 
  char   *varname_list[MXVARNAME_MISSING];
  bool   wildcard[MXVARNAME_MISSING];
} INPUTS_VARNAME_MISSING; // extracted from input append_varname_missing='xxx'

#define ORDER_ZPOLY_COVMAT 3   // 3rd order polynom for intrinsic scatter cov.
#define MXSURVEY_ZPOLY_COVMAT 10 

struct INPUTS {

  bool   KEYNAME_DUMPFLAG; // flag to dump all key names, then quit

  int  nfile_data;
  char **dataFile;  

  int nfile_data_override ; 
  char **dataFile_override ; // e.g, change all VPEC, VPEC_ERR
  
  bool   cat_only;    // cat data fitres files and do nothing else
  int    cat_prescale; // scale to reduce file size for debugging
  char   cat_file_out[MXCHAR_FILENAME] ;
  
  int    write_yaml;  // used by submit_batch_jobs.py
  int    write_csv ;  // M0DIF formatted for CosmoMC

  int    minos;

  int    nmax_tot ;   // Nmax to fit for all
  int    nmax[MXIDSURVEY];   // idem by survey
  char   nmaxString[200];

  double maxerr_abort_c, maxerr_abort_x1, maxerr_abort_x0;

  int    blindFlag  ; // suppress cosmo param printout for data
  double blind_cosinem0  ;       // M0DF -> M0DIF + cos(z*blind_cosinem0)
  double blind_cosinePar[MAXPAR][2]; // blind offset = [0] * cos([1])
  char   blindString[MAXPAR][60];

  // optional sim-files with corrections/maps
  int  nfile_biasCor;        // number of biascor files
  char **simFile_biasCor ;   // comma-sep list of biasCor inputs
  char *simFile_biasCor_arg;  // store input string in case CCprior needs it
  int  opt_biasCor ;
  int  prescale_biasCor[3] ; // subset use [0] of [1]; [2]=write flag

  double sigint_biasCor ;     // force sigint instead of autocompute
  double snrmin_sigint_biasCor; // SNRMIN used to determine sigint
  double sigma_cell_biasCor; // to weight events in cell for biasCor value

  char fieldGroup_biasCor[200]; // sub-divide biasCor groups based on field
  int  use_fieldGroup_biasCor;  // logical flag for above

  char surveyGroup_biasCor[400]; // combine surveys into group
  int  use_surveyGroup_biasCor;   // internall set flag
  int  surveyGroup_biasCor_abortFlag; // 1-> abort on missing survey

  char surveyList_noBiasCor[200]; // list of surveys fit, but skip biasCor
  char idsample_select[40];       // e.g., '0+3'
  int    select_trueIa;      // T -> select only true SNIa, disable CC prior
  int    force_realdata ;    // T -> treat SIM like real data
  double zspec_errmax_idsample; // used to create [SAMPLE]-zSPEC IDSAMPLE

  int interp_biascor_logmass;
  // ----------
  int  nfile_CCprior;
  char **simFile_CCprior;    // to get CC prior, dMU vs. z
  char   varname_pIa[100];
  bool   sameFile_flag_CCprior; // True -> same file(s) as biasCor

  char   append_varname_missing[100]; // force missing varname(s) with -9

  double force_pIa;
  bool   perfect_pIa;       // internally set if force_pIa='perfect'
  int  typeIa_ccprior ;       // PCC=0 for this sntype
  double maxProbCC_for_sigint;  // max P_CC/ProbIa to sum chi2_1a

  int    fitflag_sigmb ;  // flag to fit for sigMb that gives chi2/dof=1
  double redchi2_tol ;    // tolerance in red chi2 (was sig1tol)

  double dchi2red_dsigint;    // option to input slope instead of computing it
  double sigint_step1 ;        // size of first sigint step, OR ...
  double scale_covint_step1;  // size of first scale_covint step
  double covint_param_step1;  // one of the above

  char   sigint_fix[200];

  double prescale_simData ;   // prescale for simData (never real data)
  double prescale_simCC ;     // e.g., 2 --> use only 1/2 of sim CC
  double prescale_simIa ;     // prescale Ia only 

  // - - - - - -  cuts - - - - - 
  double cmin,  cmax  ;
  double x1min, x1max ;
  double zmin,  zmax  ;
  double logmass_min, logmass_max ;
  int    nbin_logmass;

  double chi2max ;         // global HR chi2-outlier cut (uses PROB_BEAMS)
  double *chi2max_list;    // list vs. IDSURVEY
  int    iflag_chi2max;    // 1->cut, 2->fitwgt0; 4->global, 8-> vs. IDSURVEY

  int    NCUTWIN ;
  char   CUTWIN_NAME[MXCUTWIN][MXCHAR_VARNAME];
  double CUTWIN_RANGE[MXCUTWIN][2];

  int    NFIELD ;
  char   *FIELDLIST[MXFIELD_OVERLAP] ;

  bool   LCUTWIN_RDFLAG[MXCUTWIN] ; // T=> read, 0=> use existing var
  bool   LCUTWIN_ABORTFLAG[MXCUTWIN] ;  // T=> abort if var does not exist
  bool   LCUTWIN_DATAONLY[MXCUTWIN] ;   // T=> cut on real or sim data 
  bool   LCUTWIN_BIASCORONLY[MXCUTWIN]; // T=> cut on biasCor
  bool   LCUTWIN_FITWGT0[MXCUTWIN];     // T=> MUERR->888 instead of cut
  bool   LCUTWIN_DISABLE;          // only if "CUTWIN NONE" command

  int  Nsntype ;
  int  sntype[MXSNTYPE]; // list of sntype(s) to select
  char sntypeString[100]; // comma-separated string from input file

  // CID select or reject for data only (data can be real or sim)
  int   ncidFile_data;  // number of cid-select files in comma-sep list
  char  **cidFile_data ; // list of cidFiles
  int   izbin_from_cidFile;

  int   ncidList_data;  //number of cids provided in listfile(s)
  int   acceptFlag_cidFile_data ; // +1 to accept, -1 to reject
  bool  match_on_cid_idsurvey, match_on_cid_only ; 
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
  double zwin_vpec_check[2]; // check vpec sign in this z-range

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
  
  int ndump_nobiasCor; // dump info for data with no valid biasCor
  double frac_warn_nobiasCor; // give warning of nobiasCor frac exceeds this

  char SNID_MUCOVDUMP[MXCHAR_VARNAME]; // dump MUERR info for this SNID

  int nthread ; // number of threads (default = 0 -> no threads)

  int restore_sigz ; // 1-> restore original sigma_z(measure) x dmu/dz
  int restore_mucovscale_bug ; // Sep 14 2021 allow restoring bug
  int restore_mucovadd_bug ; // +=1 to restore wrong beta for BS21 , +=2 for bug in covadd logic, March 14 2022   
  int debug_flag;    // for internal testing/refactoring
  int debug_malloc;  // >0 -> print every malloc/free (to catch memory leaks)
  int debug_mucovscale; //write mucovscale info for every biascor event
  int nbinc_mucovscale; //number of colour bins to determine muCOVSCALE and muCOVADD
  char cidlist_debug_biascor[100];

  bool LEGACY_NBINc;
  bool LEGACY_IZBIN;
  // set internal LEGACY and REFAC flags for development

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
  double DCHI2RED_DSIGINT;   // Mar 25 2022
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
  int IPARMAPINV_MN[MAXPAR];

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



// Aug 31 2020: define typedef for threads
typedef struct {
  int id_thread, nthread;
  int isn_min, isn_max;

  double  xval_fcn[MAXPAR] ;
  int     npar_fcn, iflag_fcn ;

  double chi2sum_tot, chi2sum_Ia ;
  double nsnfitIa, nsnfitcc ;   // note double for sum of BBC Probs
  int    nsnfit, nsnfit_truecc ;

} thread_chi2sums_def ;


// define fit results
struct {
  int NSNFIT ;        // Number of SN used in fit; was nsnacc
  int NSNFIT_SPLITRAN[MXSPLITRAN]; // idem in SPLITRAN sub-samples
  int NDOF ;          // Ndof in fit
  int NCALL_FCN ;     // number of calls to FCN function

  int   MNSTAT ;       // store istat returned from mnstat call (7/2020)
  double CHI2SUM_MIN;  // global min chi2
  double CHI2RED_ALL;  // global reduced chi2
  double CHI2SUM_1A;   // chi2 sum for Ia subset
  double CHI2RED_1A;   // reduced chi2 for Ia subset
  double ALPHA, BETA, GAMMA;  

  // maybe replace these with CONTAMIN_INFO ??
  double NSNFIT_1A;    // Sum of PROB(Ia)
  double NSNFIT_CC;    // Sum of PROB(CC)
  int NSNFIT_TRUECC;   // for sim, number of true CC

  double AVEMAG0 ; // average M0 among z bins
  double SNMAG0 ;  // AVEMAG0, or user-input INPUTS.nommag0 

  double M0DIF[MXz]; // M0-M0avg per z-bin
  double M0ERR[MXz]; // error on above
  double zM0[MXz];       // wgted-average z (not from fit)
  double MUREF_M0[MXz];  // wgted-avg MUREF per iz bin

  char   PARNAME[MAXPAR][MXCHAR_VARNAME];


  double PARVAL[MXSPLITRAN+1][MAXPAR];
  double PARERR[MXSPLITRAN+1][MAXPAR];

  double COVMAT[MAXPAR][MAXPAR];

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
  char *LIST[MXVAR_TABLE];   // output array of varnames

  int  NVAR_DROP;
  char DROPLIST[MXCHAR_VARLIST]; // list of dropped colums (for comment)

  bool PERFECT_COLUMN_MATCH ;  // true if all columns match
} OUTPUT_VARNAMES ;


int  DOFIT_FLAG;    // non-zero --> do another fit iteration

int ISDATA_REAL;           // T if no SIM keys are found -> real data
int FOUNDKEY_SNTYPE = 0 ;  // T -> found sntype key in fitres file
int FOUNDKEY_SIM    = 0 ;  // T -> is simulation (formerly 'simok')
int FOUNDKEY_SIMx0  = 0 ;  // T -> is SALT2 (formerly 'simx0')
int FOUNDKEY_SIM_NONIA_INDEX = 0 ;

int NVAR_APPEND ;  // NVAR appended from SALTmu
char VARNAMES_APPEND[MAXPAR*10][MXCHAR_VARNAME] ;

time_t t_start, t_end_init, t_start_fit, t_end_fit, t_read_biasCor[3] ;

// parameters for errmsg utility 
#define SEV_INFO   1  // severity flag => give info     
#define SEV_WARN   2  // severity flag => give warning  
#define SEV_ERROR  3  // severity flag => error         
#define SEV_FATAL  4  // severity flag => abort program 

#define BLINDMASK_MUz     1  // apply sinusoid function to MU-vs-z
#define BLINDMASK_FIXPAR  2  // blind fixed params (OL,w)
#define BLINDMASK_SIM    64  // blind applies to simulation

//Main function definitions

void SALT2mu_DRIVER_INIT(int argc, char **argv);
void SALT2mu_DRIVER_EXEC(void);
int  SALT2mu_DRIVER_SUMMARY(void);

void apply_blindpar(void);
void check_duplicate_SNID(void);
void check_redshifts(void) ;
void check_vpec_sign(void);
void check_zhel(void) ;
void applyCut_nmax(void);
void applyCut_chi2max(void);
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

void parse_cat_only(char *string_cat_only); 

void parse_simfile_CCprior(char *item);
void parse_ZPOLY_COVMAT(char *item);
void load_ZPOLY_COVMAT(int IDSURVEY, double Z ) ;
double sum_ZPOLY_COVMAT(double Z, double *polyPar) ;

void parse_CUTWIN(char *item);
void parse_FIELDLIST(char *item);
int  reject_CUTWIN(int EVENT_TYPE, int *DOFLAG_CUTWIN, double *CUTVAL_LIST);
int  usesim_CUTWIN(char *varName) ;
int  set_DOFLAG_CUTWIN(int ivar, int icut, int isData );
void copy_CUTWIN(int icut0,int icut1);
void parse_sntype(char *item);
void parse_cidFile_data(int OPT, char *item); 
void parse_prescale_biascor(char *item, int wrflag);
void parse_powzbin(char *item) ;
void parse_IDSAMPLE_SELECT(char *item);
void parse_sigint_fix(char *item);
void parse_blindpar(char *item) ;
void parse_chi2max(char *item);

void  prep_input_driver(void);
void  prep_input_trueIa(void);
void  prep_input_nmax(char *item);
void  prep_input_gamma(void) ;
void  prep_input_probcc0(void);
void  prep_input_load_COSPAR(void);
void  prep_input_varname_missing(void);
void  prep_input_repeat(void); 
void  prep_debug_flag(void);

int   force_probcc0(int itype, int idsurvey);
void  prep_cosmodl_lookup(void);
int   ppar(char* string);
void  read_data(void);
void  read_data_override(void);
double zhd_data_override(int isn, double vpec_over ) ;
double zhderr_data_override(int isn, double vpecerr_over ) ;
void  write_word_override(int ivar_tot, int indx, char *word) ;

void  store_input_varnames(int ifile, TABLEVAR_DEF *TABLEVAR) ;
void  store_output_varnames(void);
bool  exist_varname(int ifile,char *varName, TABLEVAR_DEF *TABLEVAR);
void  get_zString(char *str_z, char *str_zerr, char *cast) ;
void  SNTABLE_READPREP_TABLEVAR(int ifile, int ISTART, int LEN, 
				TABLEVAR_DEF *TABLEVAR);
int   SNTABLE_READPREP_HOST(char *VARNAME, int ISTART, int LEN, 
			    TABLEVAR_DEF *TABLEVAR );

void  SNTABLE_CLOSE_TEXT(void) ;
void  compute_more_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR) ;
bool  IS_SPECZ_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR) ;
void  compute_CUTMASK(int ISN, TABLEVAR_DEF *TABLEVAR );
void  compute_more_INFO_DATA(void);
void  prepare_IDSAMPLE_biasCor(void); 
void  set_FIELDGROUP_biasCor(void);
void  set_SURVEYGROUP_biasCor(void);
void  set_BINSIZE_SAMPLE_biasCor(int IDSAMPLE);

int    get_IDSAMPLE(int IDSURVEY, bool IS_PHOTOZ, 
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

void print_contam_CCprior(FILE *fp);
void print_table_CONTAM_INFO(FILE *fp,  CONTAM_INFO_DEF *CONTAM_INFO);
void setup_contam_CCprior(char *which, CONTAM_INFO_DEF *CONTAM_INFO) ;
void zero_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO) ;
void sum_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO, double Prob_Ia,
			double xhisto, int SIM_NONIA_INDEX) ;

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
void  get_BININFO_biasCor_abg(char *varName, double *VAL_MIN, 
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
		     char *callFun, FITPARBIAS_DEF *FITPARBIAS);
int get_muCOVcorr(char *cid,  BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		  double *muCOVscale, double *muCOVadd, int *i1d_cen ) ;
void dump_muCOVcorr(int n);

void   get_muBias(char *NAME, 
		  BIASCORLIST_DEF *BIASCORLIST,  
		  FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MXb][MXg],
		  double         (*MUCOVSCALE_ABGRID)[MXb][MXg],
		  double         (*MUCOVADD_ABGRID)[MXb][MXg],
		  INTERPWGT_AlphaBetaGammaDM *INTERPWGT,  
		  double *FITPARBIAS_INTERP, 
		  double *muBias, double *muBiasErr, double *muCOVscale, double *muCOVadd ) ;

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

int selectCID_data(char *cid, int IDSURVEY, int *IZBIN); 

void  write_fitres_driver(char *fileName);
void  write_fitres_misc(FILE *fout);
void  write_version_info(FILE *fp) ;
void  write_yaml_info(char *fileName);
void  write_debug_mucovcorr(int IDSAMPLE, double *muDif_list, double *muErr_list);

void  define_varnames_append(void) ;
int   write_fitres_line(int indx, int ifile, char *rowkey, 
			char *line, FILE *fout) ;
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
void  write_M0_fitres(char *fileName);
void  write_M0_csv(char *fileName);  
void  write_M0_cov(char *fileName) ;
void  write_MUERR_INCLUDE(FILE *fp) ;
void  write_NWARN(FILE *fp, int FLAG) ;

int   SPLITRAN_ACCEPT(int isn, int snid);
void  SPLITRAN_cutmask(void);

void  CPU_SUMMARY(void);

int keep_cutmask(int errcode) ;

int     prepNextFit(void);

void    conflict_check(void);
double  next_covFitPar(double redchi2, double orig_parval, double parstep);
void    recalc_dataCov(void); 

//Utility function definitions
double cosmodl_forFit(double zhel, double zcmb, double *cosPar);
double cosmodl(double zhel, double zcmb, double *cosPar);
double inc    (double zcmb, double *cosPar);

void ludcmp(double* a, const int n, const int ndim, int* indx, 
	    double* d, int* icon);
void lubksb(const double* a, const int n, const int ndim, 
       const int* indx, double* b);
double rombint(double f(double z, double *cosPar),
	       double a, double b, double *cosPar, double tol ) ;

double avemag0_calc(int opt_dump);
void   M0dif_calc(void) ;
double fcn_M0(int n, double *M0LIST );

void   muerr_renorm(void);
void   printCOVMAT(FILE *fp, int NPAR, int NPARz_write);
double fcn_muerrsq(char *name,double alpha,double beta, double gamma,
		   double (*COV)[NLCPAR], double z, double zerr, int optmask);
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

int validRowKey_TEXT(char *string) ; // see sntools_output_text.c

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

void mnpout_(const int* num, char chnam[], double* val, double* error,
	     double* bnd1, double* bnd2, int* ivarbl, int nchnam);


void *MNCHI2FUN(void *thread);

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
int   malloc_TABLEVAR_HOST(int LEN_MALLOC, TABLEVAR_DEF *TABLEVAR, 
			   char *VARNAME);
int   malloc_TABLEVAR_CUTVAL(int LEN_MALLOC, int icut, 
			     TABLEVAR_DEF *TABLEVAR ) ;
float malloc_FITPARBIAS_ALPHABETA(int opt, int LEN_MALLOC, 
				  FITPARBIAS_DEF *****FITPARBIAS );
float malloc_MUCOV(int opt, int IDSAMPLE, CELLINFO_DEF *cellinfo);

// ======================================================================
// ==== GLOBALS AND FUNCTIONS TO USE SALT2mu as SUBPROCESS ==============
// ======================================================================

#define USE_SUBPROCESS

#ifdef USE_SUBPROCESS
void SUBPROCESS_HELP(void);
void SUBPROCESS_INIT(void); // one time init driver (binning, malloc ...)
void SUBPROCESS_MALLOC_INPUTS(void);
void SUBPROCESS_PREP_NEXTITER(void); // prepare for next iteration
void SUBPROCESS_READPREP_TABLEVAR(int IFILE, int ISTART, int LEN, 
				  TABLEVAR_DEF *TABLEVAR); 
void SUBPROCESS_SIM_REWGT(int ITER_EXPECT);
double SUBPROCESS_PROB_SIMREF(int ITER, int imap, double XVAL);
void SUBPROCESS_SIM_PRESCALE(void);
int  SUBPROCESS_IVAR_TABLE(char *varName_GENPDF);
void SUBPROCESS_INIT_DUMP(void);
void SUBPROCESS_INIT_RANFLAT(int iter);
void SUBPROCESS_STORE_EBV(void);
void SUBPROCESS_READ_SIMREF_INPUTS(void);
void SUBPROCESS_OUTPUT_TABLE_PREP(int itable);
void SUBPROCESS_OUTPUT_LOAD(void);
void SUBPROCESS_OUTPUT_TABLE_LOAD(int isn, int itable);
void SUBPROCESS_OUTPUT_TABLE_RESET(int itable) ;
void SUBPROCESS_OUTPUT_WRITE(void); 
void SUBPROCESS_OUTPUT_TABLE_WRITE(int itable);
void SUBPROCESS_COMPUTE_STD(int ITABLE) ;

void SUBPROCESS_EXIT(void);
void SUBPROCESS_REMIND_STDOUT(void) ;

void SUBPROCESS_STORE_BININFO(int itable, int ivar, char *string);
void SUBPROCESS_MAP1D_BININFO(int itable);
void SUBPROCESS_OUTPUT_TABLE_HEADER(int itable);

#include "sntools_genPDF.h" 
#include "sntools_genPDF.c"

#define KEYNAME_SUBPROCESS_STDOUT          "SALT2mu_SUBPROCESS:"
#define KEYNAME_SUBPROCESS_ITERATION_BEGIN "ITERATION_BEGIN:"
#define KEYNAME_SUBPROCESS_ITERATION_END   "ITERATION_END:"
#define MXTABLE_SUBPROCESS        6  // max number of output tables
#define MXVAR_TABLE_SUBPROCESS    3  // max number of dimensions per table
#define SUBPROCESS_OPTMASK_WRFITRES 1 // write fitres file each iteration
#define SUBPROCESS_OPTMASK_WRM0DIF  2 // write M0DIF file for each iteration
#define SUBPROCESS_OPTMASK_RANSEED  4 // Use different set of randoms for each reweight event

#define VARNAME_SIM_AV   "SIM_AV"
#define VARNAME_SIM_RV   "SIM_RV"
#define VARNAME_SIM_EBV  "SIM_EBV"

typedef struct {

  // info stored before fit
  int NVAR;
  int IVAR_FITRES[MXVAR_TABLE_SUBPROCESS];
  BININFO_DEF BININFO[MXVAR_TABLE_SUBPROCESS];
  float  *PTRVAL[MXVAR_TABLE_SUBPROCESS] ;

  char VARNAMES_HEADER[200]; // varnames list for table header

  int  NBINTOT; // total number of multi-dimensional bins
  int *INDEX_BININFO[MXVAR_TABLE_SUBPROCESS] ; // for 0 to NBINTOT-1

  // info computed for each iteration.
  // Each array is malloced with NBINTOT storage
  int    *NEVT;   
  double *MURES_SQSUM, *MURES_SUM; // to reconstruct MU-bias and MU-RMS
  double *MURES_SUM_WGT, *SUM_WGT; // Brodie June 14 2021
  double **MURES_LIST;        // mures list per table bin; needed for median
  double *MURES_STD, *MURES_STD_ROBUST ;

} SUBPROCESS_TABLE_DEF ;


struct {
  bool  USE;

  // INPUT_xxx are read directly from command line
  char  *INPUT_FILES ; // comma-sep list of INPFILE,OUTFILE,STDOUT_FILE
  char  *INPUT_VARNAMES_GENPDF_STRING;
  char  *INPUT_CID_REWGT_DUMP ;
  char **INPUT_OUTPUT_TABLE ; 
  char  *INPUT_SIMREF_FILE ; 
  int    INPUT_OPTMASK;    // e.g., write fitres file for debug
  int    N_OUTPUT_TABLE ;
  int    NEVT_SIM_PRESCALE ;   // tune sim prescale to fit this many
  int    INPUT_ISEED;         // random seed
  int    STDOUT_CLOBBER; // default=T ==> rewind FP_STDOUT each iter
  
  // variables below are computed/extracted from INPUT_xxx
  char  *INPFILE ; // read PDF map from here
  char  *OUTFILE ; // write info back to python driver
  char  *STDOUT_FILE ; // direct stdout here (used only for visual debug)
  FILE  *FP_INP, *FP_OUT ;
  char   VARNAMES_GENPDF[MXVAR_GENPDF][40];
  int   NVAR_GENPDF;
  int   IVAR_TABLE_GENPDF[MXMAP_GENPDF][MXVAR_GENPDF]; // map GENPDF <-> TABLE

  int  IVAR_EBV;  // flags EBV = AV/RV [since SIM_EBV is not in FITRES file]
  bool *DUMPFLAG_REWGT;

  // store fitres columns used in GENPDF maps
  float *TABLEVAR[MXVAR_GENPDF]; // use float to save memory
  
  // store random number for each event.
  float *RANFLAT_REWGT ;
  float *RANFLAT_PRESCALE ;

  // variables filled during each subprocess iteration
  int  ITER ;
  bool *KEEP_AFTER_REWGT;

  // define info for each SALT2mu-output table.
  SUBPROCESS_TABLE_DEF OUTPUT_TABLE[MXTABLE_SUBPROCESS] ;

  // For bounding function info
  bool ISFLAT_SIM ; //True -> all sim distributions are flat; else read bounding functions
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2c ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2x1 ;
  GENGAUSS_ASYM_DEF GENGAUSS_RV ;
  GEN_EXP_HALFGAUSS_DEF EXP_HALFGAUSS_EBV ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2BETA ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2ALPHA ;
  double MAXPROB_RATIO ; 
} SUBPROCESS ;


#endif

// *******************************************
//              BEGIN CODE
// *******************************************

int main(int argc,char* argv[ ]) {

  int FLAG, N_EXEC=0;
  char fnam[] = "main";
  // ------------------ BEGIN MAIN -----------------

  SALT2mu_DRIVER_INIT(argc,argv);
  
 DRIVER_EXEC:
  N_EXEC++ ;

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) { SUBPROCESS_PREP_NEXTITER(); }
#endif

  SALT2mu_DRIVER_EXEC();

  FLAG = SALT2mu_DRIVER_SUMMARY();

  if ( FLAG == FLAG_EXEC_REPEAT ) { goto DRIVER_EXEC; } // e.g., NSPLITRAN

  fprintf(FP_STDOUT, "\n Done. \n"); fflush(FP_STDOUT);
  
  return(0) ;

} // end of main


// ********************************************
void SALT2mu_DRIVER_INIT(int argc, char **argv) {

  // Created July 2 2020 by R.Kessler
  // Part of refactor to prepare for higher-level python scripts
  // calling SALT2mu.

  char fnam[] = "SALT2mu_DRIVER_INIT";

  // ------------ BEGIN -----------

  FP_STDOUT = stdout;
  PIFAC  = 1.0/sqrt(TWOPI);
  LOGTEN = log(10.0) ;

  t_start = time(NULL);

  set_defaults();

  if (argc < 2) {
    sprintf(c1err,"Must give param-input file as arguent");
    sprintf(c2err,"to SALT2mu program. See manual.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  fprintf(FP_STDOUT, "\n\t Start program %s \n\n", argv[0]) ;
  fflush(FP_STDOUT);
  print_full_command(FP_STDOUT,argc,argv);

  // read SURVEY.DEF file and store each survey name vs. IDSURVEY
  read_SURVEYDEF(); 

  //parse input parameter file
  parse_parFile( argv[1] );

  // parse command-line args to override input file
  override_parFile(argc, argv);

  // check for conflicts between input variables
  conflict_check();

  TABLEFILE_INIT();  // call before prep_input_driver , 9.28.2020

  // prepare input
  prep_input_driver();

  //  test_zmu_solve();
  //  test_muerrz(); // xxxx

  // --------------------------------------
  //Read input data from SALT2 fit
  read_data(); 

  compute_more_INFO_DATA();  

  if( INPUTS.cat_only ) 
    { write_fitres_driver(INPUTS.cat_file_out);  exit(0); }

  // check for duplicate SNIDs and take action based on iflag_duplicate
  check_duplicate_SNID();

  // misc redshift checks
  check_redshifts();

  // setup BBC redshift bins.
  setup_BININFO_redshift();

  // prepare mapindex for each IDSURVEY & FIELD --> for biasCor
  prepare_IDSAMPLE_biasCor();

  // check option for SPLITRAN summary
  if ( INPUTS.JOBID_SPLITRAN > INPUTS.NSPLITRAN ) { return ; }

  // read optional sim biasCor maps 
  prepare_biasCor();

  // read optional CC sim for CC prior
  prepare_CCprior();

  // check for user-constraint on nmax (July 2017) AFTER prepare_biasCor 
  applyCut_nmax();

  //  recalc_dataCov(); // xxx REMOVE ??

  // check option to turn SALT2mu into a subprocess
#ifdef USE_SUBPROCESS
  SUBPROCESS_INIT();
#endif

  t_end_init = time(NULL);

  return ;

} // end SALT2mu_DRIVER_INIT

// ********************************************
void SALT2mu_DRIVER_EXEC(void) {

  // Created July 2 2020 by R.Kessler
  // Execute MINUIT-based fit.
  // Part of refactor to prepare for higher-level python scripts
  // calling SALT2mu.
  //
  // Remove MINUIT printing by moving "SET PRI -1" command to 
  // be right after MNINIT (instead of further down)
  //
  const int null=0 ;
  int inf = 5, outf = 6, savef = 7;
  int icondn, len, npari, nparx, istat, ndof ;
  double chi2min, fedm, errdef ;
  char text[100], mcom[50];
  char fnam[] = "SALT2mu_DRIVER_EXEC" ;

  // -------------- BEGIN ---------------
  t_start_fit = time(NULL);

  if ( INPUTS.JOBID_SPLITRAN > 0 ) 
    { NJOB_SPLITRAN = INPUTS.JOBID_SPLITRAN; } // do only this one SPLIT job
  else
    { NJOB_SPLITRAN++ ; }  // keep going to do them all

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) { NJOB_SPLITRAN=1; }
#endif

  DOFIT_FLAG = FITFLAG_CHI2 ; 
  printmsg_fitStart(FP_STDOUT);
  
  if ( INPUTS.NSPLITRAN > 1 ) { 
    SPLITRAN_cutmask();     // check for random sub-samples
    prep_input_repeat();  // re-initialize a few things (May 2019)
  }

  FITRESULT.NCALL_FCN = 0 ;
  mninit_(&inf,&outf,&savef);

  strcpy(mcom,"SET PRI -1");     len = strlen(mcom);
  mncomd_(fcn, mcom, &icondn, &null, len);  fflush(FP_STDOUT);

  strcpy(text,"SALT2mu"); len = strlen(text);  
  mnseti_(text,len);    fflush(FP_STDOUT);

  setup_zbins_fit();    // set z-bins

  // execuate minuit mnparm_ commands
  exec_mnparm(); 

  // use FCN call and make chi2-outlier cut (Jul 19 2019)
  applyCut_chi2max();

  FITRESULT.NFIT_ITER = 0 ;

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) { SUBPROCESS_SIM_PRESCALE(); } // Jun 2021
#endif

  // print stats for data after ALL cuts are applied
  print_eventStats(EVENT_TYPE_DATA);
  
  // Beginning of DOFIT loop
  while ( DOFIT_FLAG != FITFLAG_DONE  ) {

    //Miniut MINIMIZE (find minimum chi-squared)
    strcpy(mcom,"SIM 1000");   len = strlen(mcom);
    mncomd_(fcn, mcom, &icondn, &null, len);  fflush(FP_STDOUT);

    strcpy(mcom,"MINI");   len = strlen(mcom);
    mncomd_(fcn,mcom,&icondn,&null,len);  fflush(FP_STDOUT); 

    //Minuit MINOS (compute errors)
    strcpy(mcom,STRING_MINUIT_ERROR[INPUTS.minos]);

    len = strlen(mcom); 
    mncomd_(fcn, mcom, &icondn, &null, len);  fflush(FP_STDOUT);

    //Final call to FCN at minimum of chi-squared
    strcpy(mcom,"CALL FCN 3");  len = strlen(mcom);
    mncomd_(fcn, mcom, &icondn, &null, len);   fflush(FP_STDOUT);

    mnstat_(&chi2min, &fedm, &errdef, &npari, &nparx, &istat);
    ndof = FITRESULT.NSNFIT - npari; 
    FITRESULT.MNSTAT      = istat;
    FITRESULT.CHI2SUM_MIN = chi2min ;
    FITRESULT.NDOF        = ndof ;
    FITRESULT.CHI2RED_ALL = chi2min/(double)ndof;


    DOFIT_FLAG = prepNextFit();
    FITRESULT.NFIT_ITER++ ; 

    fflush(FP_STDOUT);  
  }   // End of fitflag_sigmb  loop

  // - - - - -
  // May 26 2021: free genPDF maps
#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE && NMAP_GENPDF>0 ) { free_memory_genPDF(); }
#endif

  t_end_fit = time(NULL);

} // end SALT2mu_DRIVER_EXEC


// ***********************************************
int SALT2mu_DRIVER_SUMMARY(void) {

  // ===================================================
  //           SUMMARY and WRAP-UP
  // ===================================================

  int MNSTAT = FITRESULT.MNSTAT ;

  char COMMENT_MNSTAT[4][40] = {
    "Fit Invalid", 
    "Bad fit; Diagonal errors only", 
    "Suspect fit; errors forced positive",
    "Good fit; errors valid"
  };

  char fnam[] = "SALT2mu_DRIVER_SUMMARY" ;

  // ------------ BEGIN ----------

  fprintf(FP_STDOUT, "\n**********Fit summary**************\n");

  fprintf(FP_STDOUT, "MNFIT status=%i (%s)\n", 
	  MNSTAT, COMMENT_MNSTAT[MNSTAT] );

  double redChi2 = FITRESULT.CHI2SUM_MIN/(double)FITRESULT.NDOF ;
  fprintf(FP_STDOUT, "-2lnL/dof = %.2f/%i = %6.3f (%i SN) \n",
	 FITRESULT.CHI2SUM_MIN, FITRESULT.NDOF, redChi2, 
	 FITRESULT.NSNFIT );

  fflush(FP_STDOUT);

  exec_mnpout_mnerrs(); // Dec 12 2016

  //The exact value of M0 shouldn't matter,
  //But take the average over bins (weighted by number of SN)
  FITRESULT.AVEMAG0 = avemag0_calc(1);  // call after PARVAL is loaded.

  // compute M0 - AVEMAG0, and apply optional blind offset
  M0dif_calc();

  // print reduced COV matrix
  printCOVMAT(FP_STDOUT, FITINP.NFITPAR_FLOAT, 999);

  // renormalize individual MUERR values so that weighted
  // variance in each z-bin matches fitted M0DIFERR
  muerr_renorm();

  // ------------------------------------------------
  // check files to write
  outFile_driver();

  //---------
  if ( NJOB_SPLITRAN < INPUTS.NSPLITRAN  &&  INPUTS.JOBID_SPLITRAN<0 ) 
    { return(FLAG_EXEC_REPEAT); }

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) {
    printf("%s CHI2_MIN = %.2f   <M0> = %.4f  NFIT_ITER=%d\n",
	   KEYNAME_SUBPROCESS_STDOUT, FITRESULT.CHI2SUM_MIN,
	   FITRESULT.AVEMAG0, FITRESULT.NFIT_ITER );
    fflush(stdout);
    if ( ISDATA_REAL ) 
      { SUBPROCESS_EXIT(); return(FLAG_EXEC_STOP); }
    else
      { return(FLAG_EXEC_REPEAT); }
 
  }
#endif

  CPU_SUMMARY();
  
  
  return(FLAG_EXEC_STOP);

} // end SALT2mu_DRIVER_SUMMARY

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

  int i, iz, iMN, iMN_tmp, len, ierflag=0, icondn, ISFLOAT ;
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
      FITINP.IPARMAPINV_MN[i]     = iMN_tmp ; // inverse map
    } 
    else {
      iMN = i + 1;
      sprintf(text,"FIX %i",iMN);  len = strlen(text);
      mncomd_(fcn, text, &icondn, &null, len);
      FITINP.IPARMAPINV_MN[i]     = -9 ;
    }
  }

  fflush(FP_STDOUT);

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
  // Jul 29 2020: incluce MINOS/MIGRAD string in a few places
  // Dec 02 2020: call mnemat_ and fill COVMAT 

  int    minos = INPUTS.minos ;
  double PARVAL, PARERR, bnd1, bnd2, eplus,eminus,eparab, globcc;
  int    ipar, iMN, iv ;
  int    LEN_VARNAME = 10 ;
  char text[100], format[80], cPARVAL[MXCHAR_VARNAME] ;
  //  char fnam[] = "exec_mnpout_mnerrs" ;

  // -------------- BEGIN ----------------

  fprintf(FP_STDOUT, "\nFinal parameter values and %s errors.\n",
	  STRING_MINUIT_ERROR[minos] );  
  fflush(FP_STDOUT);

  for (ipar=0; ipar<FITINP.NFITPAR_ALL; ipar++ )  {
    iMN     = ipar + 1 ;
    text[0] = 0 ;
    mnpout_(&iMN, text, &PARVAL, &PARERR, &bnd1,&bnd2, &iv, LEN_VARNAME);
    sprintf(cPARVAL,"%7.4f", PARVAL);
    text[LEN_VARNAME] = 0;
          
    if (iv==0) {
      if (ISBLIND_FIXPAR(ipar) ) { sprintf(cPARVAL,"BLINDED"); }
      fprintf(FP_STDOUT, "par %2i      %10s %s         fixed\n", 
	      ipar, text, cPARVAL);
    }
    else  {
      mnerrs_(&iMN, &eplus, &eminus, &eparab, &globcc);

      if ( fabs(eplus) > 1.0E-6 && fabs(eminus) > 1.0E-6 ) 
	{ PARERR = 0.5*( fabs(eplus) + fabs(eminus) ) ; }
      else
	{ PARERR = eparab; } // Apr 15 2020 : better than nothing
	  	
      if ( BLIND_OFFSET(ipar) == 0.0   ) {
	char string_asymerr[40] = "";
	if ( minos ) { 
	  sprintf(string_asymerr,"%s (+)%7.4e (-)%7.4e ",
		  "MINOS",  eplus, eminus );
	}
	
	strcpy(format,"par %2i (%2i) %10s %11.4e +/- %.4e "
	       "%s Global CC=%6.3f \n");
	fprintf(FP_STDOUT,  format,
		ipar, iv, text, PARVAL,eparab, string_asymerr, globcc);
      }
      else {
	fprintf(FP_STDOUT, "par %2d (%2d) %10s   ***** BLINDED ***** \n",
	       ipar, iv, text);
      }
    }

    // fill global arrays for later
    FITRESULT.PARVAL[NJOB_SPLITRAN][ipar] = PARVAL ;
    FITRESULT.PARERR[NJOB_SPLITRAN][ipar] = PARERR ;

    fflush(FP_STDOUT);
  } // end num loop

  
  // load sigInt separately in IPAR_SIGINT space 
  FITRESULT.PARVAL[NJOB_SPLITRAN][IPAR_COVINT_PARAM] = FITINP.COVINT_PARAM_FIX;
  FITRESULT.PARERR[NJOB_SPLITRAN][IPAR_COVINT_PARAM] = 1.0E-8 ;
  
  // load full cov matrix 
  int num = MAXPAR;
  mnemat_(FITRESULT.COVMAT,&num);

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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }
  
  return ;

} // end   setup_BININFO_redshift


// ***********************************************
void setup_BININFO_userz(void) {

  // Created Jul 26 2019
  // Parse zbinuser string and load z-bin structure

  int nzbin=0, Nsplit, iz ;
  int MEMC = 20*sizeof(char);
  int debug_malloc = INPUTS.debug_malloc ;
  double zlo, zhi, zmin=-9.0, zmax=-9.0 ;
  char *ptr_z[MXz];
  char fnam[] = "setup_BININFO_userz" ;

  // --------------- BEGIN -------------

  print_debug_malloc(+1*debug_malloc,fnam);
  for(iz=0; iz < MXz; iz++ ) { ptr_z[iz] = (char*)malloc(MEMC); }

  splitString(INPUTS.zbinuser, COMMA, MXz,    // inputs
	      &Nsplit, ptr_z );                    // outputs
  nzbin  = Nsplit-1 ;

  if ( Nsplit <= 1 || Nsplit >= MXz ) {
    sprintf(c1err,"Invalid Nsplit=%d for", Nsplit);
    sprintf(c2err,"zbinuser=%s\n", INPUTS.zbinuser);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  fprint_banner(FP_STDOUT,BANNER);
  fprintf(FP_STDOUT, "\t nlogzbin = %d -> logz binsize = %.4f \n",
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
  fprint_banner(FP_STDOUT,BANNER);
  fprintf(FP_STDOUT, "\t zbinsize ~ (1+z)^%.2f for z < %.3f \n",
	 powzbin, zmax );
  if ( LHALF ) { 
    fprintf(FP_STDOUT, "\t zbinsize = %.3f for z > %.3f \n", zbin2, zmax );
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
  int NSN_CUTS = 0 ;
  int n, nz, izbin, iztmp, NZFLOAT, CUTMASK ;
  double z;
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

    NSN_CUTS++ ;

    z     = INFO_DATA.TABLEVAR.zhd[n] ;
    izbin = IBINFUN(z, &INPUTS.BININFO_z, 0, fnam);

    if ( INFO_DATA.USE_IZBIN_from_CIDFILE ) { 
      iztmp = INFO_DATA.TABLEVAR.IZBIN[n]; 
      if ( iztmp >= 0 ) { 
	if ( iztmp != izbin ) { INFO_DATA.NCHANGE_IZBIN++ ; }
	izbin = iztmp; 
      }
    } 

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

    fprintf(FP_STDOUT, " z=%8.5f - %8.5f  NZBIN(TOT,CUTS)=%6d,%6d   "
	    "ISFLOAT=%i\n",
	   INPUTS.BININFO_z.lo[nz], INPUTS.BININFO_z.hi[nz],
	   FITINP.NZBIN_TOT[nz], FITINP.NZBIN_FIT[nz], 
	   FITINP.ISFLOAT_z[nz]);
  }

  FITINP.NFITPAR_FLOAT_z = NZFLOAT ;
  fprintf(FP_STDOUT," --> Use %d of %d z-bins in fit.\n", NZFLOAT, nzbin );

  if ( INFO_DATA.NCHANGE_IZBIN > 0 ) {
    fprintf(FP_STDOUT,"   ALERT: %d of %d events change IZBIN to "
	    "match cid_select_file\n",
	    INFO_DATA.NCHANGE_IZBIN, NSN_CUTS );   
  }

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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  fprintf(FP_STDOUT,"\n");
  fflush(FP_STDOUT);
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
    fprintf(FP_STDOUT," xxx %s: blindFlag=%d ISDATAREAL=%d blindpar(LAM)=%f\n", 
	   fnam, INPUTS.blindFlag, ISDATA_REAL, 
	   INPUTS.blind_cosinePar[IPAR_OL][0] ); fflush(FP_STDOUT);
  }

  for(ipar=0; ipar < MAXPAR; ipar++ ) {

    if ( ISBLIND_FIXPAR(ipar) ) {
      blindpar             = INPUTS.blind_cosinePar[ipar];
      INPUTS.parval[ipar] += (blindpar[0] * cos(blindpar[1]) ) ;
      fprintf(FP_STDOUT, "  BLIND FIXPAR: %s += %8.4f * cos(%f) \n",
	     FITPARNAMES_DEFAULT[ipar], blindpar[0], blindpar[1] );
      fflush(FP_STDOUT);
    }
  }

  // re-load INPUTS.COSPAR (bugfix, Oct 22 2019)
  if ( LDMP ) {
    fprintf(FP_STDOUT, " xxx %s: 1. COSPAR(OL,w0) = %f, %f \n",
	   fnam, INPUTS.COSPAR[0], INPUTS.COSPAR[2]); fflush(stdout);
  }

  // set INPUTS.COSPAR with blinded cosmo params
  prep_input_load_COSPAR();

  if ( LDMP ) {
    fprintf(FP_STDOUT, " xxx %s: 2. COSPAR(OL,w0) = %f, %f \n",
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
    
  fprint_banner(FP_STDOUT, "Data Cut Summary for nmax");

  fprintf(FP_STDOUT, "\t %-20s : %5d -> %5d \n", "ALL", NTOT, ntot );
  for(id = 0; id < MXIDSURVEY; id++ ) {
    if ( NPERSURVEY[id] == 0 ) { continue ; }

    fprintf(FP_STDOUT, "\t %-20s : %5d -> %5d \n", 
	   SURVEY_INFO.SURVEYDEF_LIST[id],
	   NPERSURVEY[id], npersurvey[id] );
  }

  return ;

} // end applyCut_nmax


// ******************************************
void applyCut_chi2max(void) {

  // Created Jul 19 2019
  // Call FCN (fit function) to evaluate chi2 with
  // initial params from input file (i.e., p1,p2,sigint ...) 
  // then apply chi2-outlier cut to data.
  // Note that this is BEAMS chi2 when there is CC contam,
  // and thus this is NOT the same as cutting on HR resids.
  //
  // Jan 22 2021: refactor to check for FITWGT0 option
  //
  // Jul 08 2021: if setbit_CUTMASK is called, call setup_zbins_fit()
  //              to check for z-bin dropouts and update NFIT(z)
  //
  int  NSN_DATA       = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  iflag_chi2max  = INPUTS.iflag_chi2max ;

  int  IFLAG_APPLY    = DOFLAG_CUTWIN_APPLY ;
  int  IFLAG_FITWGT0  = DOFLAG_CUTWIN_FITWGT0 ;
  int  IFLAG_GLOBAL   = 4 ;
  int  IFLAG_SURVEY   = 8 ;

  bool DOCUT_APPLY   = (iflag_chi2max & IFLAG_APPLY)   > 0 ;
  bool DOCUT_FITWGT0 = (iflag_chi2max & IFLAG_FITWGT0) > 0 ;
  bool DOCUT_GLOBAL  = (iflag_chi2max & IFLAG_GLOBAL ) > 0 ;

  double chi2max ;
  int len, icondn, n, cutmask, idsurvey, NREJ=0 ;
  const int null=0 ;
  double chi2;
  bool FAILCUT;
  char mcom[60], *name ;
  //  char fnam[] = "applyCut_chi2max" ;

  // ----------- BEGIN ------------

  if ( iflag_chi2max == 0 ) { return; }
  // xxxx  if ( chi2max > 0.99E9 ) { return; }

  strcpy(mcom,"CALL FCN 1");  len = strlen(mcom);
  mncomd_(fcn, mcom, &icondn, &null, len); 
  fflush(FP_STDOUT);

  for (n=0; n< NSN_DATA; ++n)  {
    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n] ; 
    if ( cutmask ) { continue; }

    chi2     = INFO_DATA.chi2[n];
    name     = INFO_DATA.TABLEVAR.name[n];
    idsurvey = INFO_DATA.TABLEVAR.IDSURVEY[n];

    if ( DOCUT_GLOBAL ) 
      { chi2max = INPUTS.chi2max ; }
    else 
      { chi2max = INPUTS.chi2max_list[idsurvey]; }

    /*
    if ( chi2 > 8.0 ) {
      printf(" xxx %s: chi2(%8s) = %.2f \n", fnam, name, chi2);
      fflush(stdout);
    }
    */

    FAILCUT = ( chi2 > chi2max );
    if ( FAILCUT ) {
      char str_chi2[40];
      sprintf(str_chi2, "Chi2(%s) = %.2f", name, chi2);

      if ( DOCUT_APPLY )  { 
	fprintf(FP_STDOUT, "\t %s -> reject \n", str_chi2);
	setbit_CUTMASK(n, CUTBIT_CHI2, &INFO_DATA.TABLEVAR);
	NREJ++ ;
      }
      else {
	// do NOT cut; instead, set fit wgt = 0 via MUERR -> huge number
	fprintf(FP_STDOUT, "\t %s -> fit wgt = 0 \n", str_chi2);
	INFO_DATA.set_fitwgt0[n] = true;
      }
    } // end FAILCUT
    
  } // end n loop over SN


  if ( NREJ > 0 ) { 
    printf("\n Setup z-bins again after chi2max cut: \n");
    fflush(stdout);
    setup_zbins_fit(); 
  }  

  fflush(FP_STDOUT);
 
  return ;

} // end applyCut_chi2max


// *******************************
void check_redshifts(void) {
  check_vpec_sign();
}  // check_redshifts

// ******************************************
void check_vpec_sign(void) {

  // Created Oct 28 2020
  // For z in zwin_vpec_check, measure Hubble rms twice:
  // with current zHD, and with zHD recomputed with vpec sign flip.
  // if the sign flip has smaller RMS, abort with error message.
  //
  // Misc task: if zhel < 0 (does not exist), set zhel = zcmb
  // to allow processing very old FITRES files.
  //
  // Feb 10 2021: fix nasty delcaration bug: rms[0] -> rms[2]

  double *zwin           = INPUTS.zwin_vpec_check ;
  double alpha           = INPUTS.parval[IPAR_ALPHA0] ;
  double beta            = INPUTS.parval[IPAR_BETA0] ;
  int    NSN_ALL         = INFO_DATA.TABLEVAR.NSN_ALL ;

  int isn, i, cutmask, NSN_SUM=0;
  double SUM_MURES[2], SUM_SQMURES[2], mean[2], rms[2], sgn_flip ;
  double zHD, zCMB, zHD_tmp, vpec, zpec;
  double mB, x1, c, mumodel, mures, dl ;
  char fnam[] = "check_vpec_sign" ;

  // ------- BEGIN -------

  if ( zwin[1] < 0.0001 ) { return; }

  for(i=0; i < 2; i++ ) 
    { SUM_MURES[i] = SUM_SQMURES[i] = 0.0 ;  }

  for(isn=0; isn < NSN_ALL; isn++ ) {
    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[isn] ; 
    zHD      = (double)INFO_DATA.TABLEVAR.zhd[isn] ;

    if ( cutmask ) { continue; }
    if ( zHD > zwin[1] ) { continue; }
    if ( zHD < zwin[0] ) { continue; }

    zCMB = (double)INFO_DATA.TABLEVAR.zcmb[isn] ;
    mB   = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_mB][isn] ;
    x1   = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_x1][isn] ;
    c    = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_c][isn] ;
    vpec = (double)INFO_DATA.TABLEVAR.vpec[isn] ; 
    zpec = fabs(vpec / LIGHT_km) ;
    if ( vpec == 0.0 ) { continue; }

    if (zHD < zCMB)
      { sgn_flip = +1.0 ; }
    else
      { sgn_flip = -1.0 ; }

    NSN_SUM++ ;  
    for(i=0; i < 2; i++ ) {

      if ( i ==0 ) 
	{ zHD_tmp = zHD;  } // nominal
      else 
	{ zHD_tmp += (sgn_flip*2.0*zpec); } // flip vpec sign

      dl = cosmodl_forFit(zHD_tmp, zHD_tmp, INPUTS.COSPAR);
      mumodel        = 5.0*log10(dl) + 25.0 ;
      mures          = mB  + alpha*x1 - beta*c - M0_DEFAULT - mumodel;
      SUM_MURES[i]   += mures ;
      SUM_SQMURES[i] += (mures*mures) ;
    }

  } // end isn loop

  // - - - -

  if ( NSN_SUM < 70 )   { return; }

  for(i=0; i < 2; i++ ) {
    mean[i] = SUM_MURES[i] / (double)NSN_SUM ;
    rms[i]  = STD_from_SUMS(NSN_SUM, SUM_MURES[i], SUM_SQMURES[i]);
  }

  printf("\n %s with RMS(MURES) using N(%.3f<z<%.3f) = %d events:\n", 
	 fnam, zwin[0], zwin[1], NSN_SUM);
  printf("\t RMS(MURES,nominal)   = %.4f \n", rms[0] );
  printf("\t RMS(MURES,flip-vpec) = %.4f \n", rms[1] );
  fflush(stdout);

  if ( rms[1] < rms[0] ) {
    sprintf(c1err,"RMS(MURES) is smaller with vpec sign-flip.");
    sprintf(c2err,"See RMS(MURES) values above.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  //  debugexit(fnam);

  return ;

} // end check_vpec_sign


// ******************************************
void check_zhel(void) {
  
  // Created Dec 11 2020
  // If zhel does not exist, compute it from zHD and set zhelerr = zhderr.
  // This allows process old FITRES files that do not have zhel.

  float z0      = INFO_DATA.TABLEVAR.zhel[0]  ;
  bool SKIP     = (z0 > 0.0);
  if (SKIP) { return; }

  int  NSN_ALL  = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  isn;
  double zhd, zhderr, zhel, zhelerr;
  char fnam[] = "check_zhel" ;

  // ------------- BEGIN ------------

  fprint_banner(FP_STDOUT,fnam);
  fprintf(FP_STDOUT,"\t zhel data column does not exist -> "
	  "set zhel = zhd" );

  for(isn=0; isn < NSN_ALL; isn++ ) {
    zhd     = (double)INFO_DATA.TABLEVAR.zhd[isn] ;
    zhel    = zhd; // WARNING: compute this later
    INFO_DATA.TABLEVAR.zhel[isn] = (float)zhel;    
  }


} // end check_zhel

// ******************************************
void check_duplicate_SNID(void) {

  // Sep 2016:
  // Use sorted redshift and its error to flag duplicates.
  // CHeck iflag_duplictes for what to do:
  // 0 -> nothing
  // 1 -> abort
  // 2 -> merge fitparams and cov into one LC
  //
  // Sep 24 2021: malloc local UNSORT_DUPL and NDUPL_LIST
  //              Introduce 2nd dimension MXSET_DUPLICATE

#define MXSET_DUPLICATE 10 // abort if more than this many per dupl set

  int  iflag   = INPUTS.iflag_duplicate;
  int  MXSTORE = MXSTORE_DUPLICATE ;
  int debug_malloc = INPUTS.debug_malloc ;

  int  isn, isn2, nsn, MEMD, MEMI, MEMB;
  int  unsort, unsort2, *unsortList, ORDER_SORT   ;
  bool IS_SIM, LDMP ;
  double *zList ;
  bool   *IS_DUPL;
  char fnam[] = "check_duplicate_SNID" ;

  // ----------- BEGIN -----------

  sprintf(BANNER,"Begin %s", fnam);
  fprint_banner(FP_STDOUT,BANNER);

  IS_SIM  = INFO_DATA.TABLEVAR.IS_SIM ;
  nsn     = INFO_DATA.TABLEVAR.NSN_ALL ;

  // Jun 11 2020: for sims, don't bother tracking duplicates
  if ( IS_SIM ) { return; }

  print_debug_malloc(+1*debug_malloc,fnam);
  MEMD = (nsn+1) * sizeof(double)  ;
  MEMI = (nsn+1) * sizeof(int)     ;
  MEMB = (nsn+1) * sizeof(bool)    ;
  zList      = (double *)malloc(MEMD); // allocate redshift list
  unsortList = (int    *)malloc(MEMI); // allocate sorted list
  IS_DUPL    = (bool   *)malloc(MEMB);

  for(isn=0; isn<nsn; isn++)  {
    zList[isn] = (double)INFO_DATA.TABLEVAR.zhd[isn]; 
    IS_DUPL[isn] = false;
  }
  

  ORDER_SORT = + 1 ; // increasing order
  sortDouble( nsn, zList, ORDER_SORT, unsortList ) ;

  bool SAME, SAME_z, SAME_SNID, FOUND_DUPL ;
  int   NTMP, idup ;
  double z, z2 ;
  char *snid, *snid2 ;
  int  **UNSORT_DUPL;
  int  *NDUPL_LIST; // how many duplicates per set
  int  NDUPL_SET ; // number of duplicate sets
  int  NDUPL_TOT ; // includes those beyone storage capacity
  int  NDUPL_SN  ; // total number of SN controbuting to duplicates

  MEMI = MXSTORE_DUPLICATE * sizeof(int);
  NDUPL_LIST  = (int*)  malloc(MEMI);
  UNSORT_DUPL = (int**) malloc( MXSTORE * sizeof(int*) );

  for(idup=0; idup < MXSTORE; idup++ )  {  
    NDUPL_LIST[idup] = 0 ; 
    UNSORT_DUPL[idup] = (int*) malloc(MXSET_DUPLICATE*sizeof(int)) ;
  }
  NDUPL_SET = NDUPL_TOT = NDUPL_SN = 0 ;

  // - - - - 

  for ( isn=0; isn < nsn-1; isn++ ) {
    unsort  = unsortList[isn];
    z      = INFO_DATA.TABLEVAR.zhd[unsort];
    snid   = INFO_DATA.TABLEVAR.name[unsort];    

    if ( IS_DUPL[isn] ) { continue; }
    isn2 = isn  ;      z2 = z;   FOUND_DUPL=false;

    while ( z == z2 && isn2 < nsn-1 ) {
      isn2++ ;
      unsort2    = unsortList[isn2];
      z2         = INFO_DATA.TABLEVAR.zhd[unsort2];
      snid2      = INFO_DATA.TABLEVAR.name[unsort2];
      SAME_SNID  = ( strcmp(snid,snid2) == 0 );

      if ( IS_DUPL[isn2] ) { continue; }
      if ( !SAME_SNID    ) { continue; }

      // we have a duplicate
      NDUPL_TOT++ ; 
      IS_DUPL[isn] = IS_DUPL[isn2] = true; 

      /* xxx
      if ( isn < 10 ) {
	printf(" xxx %s: compare %s/%s  %d/%d  FOUND=%d\n", 
	       fnam, snid, snid2, isn,isn2, FOUND_DUPL ); fflush(stdout);
      }
      xxxxxx*/

      if ( NDUPL_SET >= MXSTORE ) { continue; }

      // store this duplicate
      if ( !FOUND_DUPL ) {
	// first duplicate for isn; store isn info
	FOUND_DUPL = true;
	NDUPL_SET++ ;  // increment number of duplicate sets
	NDUPL_SN++ ;
	NTMP = NDUPL_LIST[NDUPL_SET-1] ;
	UNSORT_DUPL[NDUPL_SET-1][NTMP] = unsort ;
	NDUPL_LIST[NDUPL_SET-1]++ ;
      }

      // always store isn2 info
      NTMP = NDUPL_LIST[NDUPL_SET-1] ;
      if( NTMP < MXSET_DUPLICATE ) 
	{ UNSORT_DUPL[NDUPL_SET-1][NTMP] = unsort2 ; }
      NDUPL_LIST[NDUPL_SET-1]++ ;
      NDUPL_SN++ ;
      
    }  // end z==z2 (isn2 loop)
  } // end loop over isn

  // - - - - -
  if ( NDUPL_SET == 0 ) 
    { fprintf(FP_STDOUT, "\t No duplicates found. \n"); goto DONE; }
  
  fprintf(FP_STDOUT, "   Found %d sets of duplicates from %d input SN: \n", 
	  NDUPL_SET, NDUPL_SN );

  for(idup=0; idup < NDUPL_SET ; idup++ ) {

    unsort = UNSORT_DUPL[idup][0] ; // first duplicate has SNID and z
    NTMP   = NDUPL_LIST[idup] ;
    snid   = INFO_DATA.TABLEVAR.name[unsort]; 
    z      = INFO_DATA.TABLEVAR.zhd[unsort];      
    fprintf(FP_STDOUT, "\t -> DUPL-%3.3d: %2d with SNID=%14.14s at z=%.5f \n",
	    idup, NTMP, snid, z );	
  }

  fflush(stdout);
  //  debugexit(fnam); // xxx REMOVE

  // - - - - - - -
  fflush(FP_STDOUT);

  if ( NDUPL_SET >= MXSTORE  && iflag>0 ) {
    sprintf(c1err,"NDUPL=%d exceeds bound, MXSTORE_DUPLICATE=%d",
	    NDUPL_TOT, MXSTORE );
    sprintf(c2err,"Check duplicates");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

 
  fprintf(FP_STDOUT, "  iflag_duplicate = %d --> ", iflag);
  if ( iflag == IFLAG_DUPLICATE_IGNORE ) {
    fprintf(FP_STDOUT, " do nothing.\n");
  }
  else if ( iflag == IFLAG_DUPLICATE_ABORT ) {
    sprintf(c1err,"Duplicates not allowed.");
    sprintf(c2err,"Check input FITRES file.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }
  else if ( iflag == IFLAG_DUPLICATE_AVG ) {
    fprintf(FP_STDOUT, " merge duplicates.\n");
    for(idup=0; idup < NDUPL_SET ; idup++ ) 
      { merge_duplicates(NDUPL_LIST[idup], UNSORT_DUPL[idup] ); }
  }
  else {
    sprintf(c1err,"Invalid iflag_duplicate=%d", INPUTS.iflag_duplicate );
    sprintf(c2err,"grep IFLAG_DUPLICATE  SALT2mu.c");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }


 DONE:
  print_debug_malloc(-1*debug_malloc,fnam);
  free(zList);  free(unsortList); free(IS_DUPL);
  
  free(NDUPL_LIST);
  for(idup=0; idup < MXSTORE; idup++ )  { free(UNSORT_DUPL[idup]); }
  free(UNSORT_DUPL);

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
  // Oct 12 2020:
  //   found numerical problem where two COV terms have opposite
  //   signs with nearly equal abs value ... summing gives ~zero,
  //   and then taking inverse results in insanely HUGE cov.
  //   If SIGN_CHANGE = True, take average COV to avoid artifact.
  //
  //   *** WARNING: need to refactor ??? ****

  int i, isn, ipar, ipar2, ISN_SAVE ;
  double fitpar[NLCPAR], fiterr[NLCPAR];
  double COVFIT_INV[NLCPAR][NLCPAR], COVINT[NLCPAR][NLCPAR];
  double COVFIT_SUM[NLCPAR][NLCPAR] ;
  double wgt, sqerr, wgtsum[NLCPAR], covfit, covfit0, covtot, fitpar_tmp  ;
  bool   SIGN_CHANGE[NLCPAR][NLCPAR], sign_change ;
  char stringList_fitparOrig[NLCPAR][100], *name ;
  char string_tmp[100];
  //  char fnam[] = "merge_duplicates" ;

  // ----------- BEGIN ----------------

  ISN_SAVE = isnList[0] ;
  name = INFO_DATA.TABLEVAR.name[ISN_SAVE] ;

  fprintf(FP_STDOUT, 
	  "# ------ merge %d duplicates for %s ---------- \n", NDUPL, name);
  fflush(FP_STDOUT);

  for(ipar=0; ipar < NLCPAR; ipar++ ) {
    fitpar[ipar] = fiterr[ipar] = wgtsum[ipar] = 0.0 ;
    for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) {
      COVFIT_INV[ipar][ipar2] = COVFIT_SUM[ipar][ipar2] = 0.0 ;  
      SIGN_CHANGE[ipar][ipar2] = false ;
    }

    stringList_fitparOrig[ipar][0] = 0 ;
  }

  int  isn0 = isnList[0] ;

  for(i=0; i < NDUPL; i++ ) {
    isn  = isnList[i] ;
    // keep first duplicate; remove the rest
    if ( i > 0 ) { setbit_CUTMASK(isn, CUTBIT_DUPL, &INFO_DATA.TABLEVAR); }

    for(ipar=0; ipar < NLCPAR; ipar++ ) {     // mB,x1,c
      sqerr  = INFO_DATA.TABLEVAR.covmat_tot[isn][ipar][ipar] ;
      wgt    = 1.0/sqerr ;
      fitpar_tmp = INFO_DATA.TABLEVAR.fitpar[ipar][isn];
      wgtsum[ipar] += wgt ;
      fitpar[ipar] += wgt * fitpar_tmp ;

      sprintf(string_tmp, "%s", stringList_fitparOrig[ipar]);
      sprintf(stringList_fitparOrig[ipar],"%s %6.3f",
	      string_tmp, fitpar_tmp );

      for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) {
	covfit  = INFO_DATA.TABLEVAR.covmat_fit[isn][ipar][ipar2];
	covfit0 = INFO_DATA.TABLEVAR.covmat_fit[isn0][ipar][ipar2];
	sign_change = ( covfit*covfit0 < 0.0 ) ;
	if ( sign_change ) { SIGN_CHANGE[ipar][ipar2] = true; }
	COVFIT_INV[ipar][ipar2] += 1.0/covfit ;
	COVFIT_SUM[ipar][ipar2] += covfit ; // backup in case of sign change
      }
    } // loop over mB,x1,c
  }

  get_COVINT_model(-1,COVINT);

  for(ipar=0; ipar < NLCPAR; ipar++ ) { 
    fitpar[ipar] /= wgtsum[ipar] ;
    fiterr[ipar]  = sqrt(1.0/wgtsum[ipar]);
    INFO_DATA.TABLEVAR.fitpar[ipar][ISN_SAVE] = fitpar[ipar] ;
    for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) { 
      if ( SIGN_CHANGE[ipar][ipar2] ) 
	{ covfit = COVFIT_SUM[ipar][ipar2] / (double)NDUPL ; }
      else
	{ covfit = 1.0 / COVFIT_INV[ipar][ipar2] ; }
      covtot = covfit + COVINT[ipar][ipar2] ; 
      INFO_DATA.TABLEVAR.covmat_fit[ISN_SAVE][ipar][ipar2] = covfit ;
      INFO_DATA.TABLEVAR.covmat_tot[ISN_SAVE][ipar][ipar2] = covtot ;
    }
  }


  // print old & new SALT2 params
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    fprintf(FP_STDOUT, "  %s-%2.2s(%s) -> %6.3f \n",  
	   name, BIASCOR_NAME_LCFIT[ipar], 
	   stringList_fitparOrig[ipar], fitpar[ipar] );
  }

  // print new COV
  fprintf(FP_STDOUT, "  %s-COV(no sigInt | +sigInt): \n", name );
  for(ipar=0; ipar<NLCPAR; ipar++ ) {
    printf("\t");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ )  { 
      covfit = INFO_DATA.TABLEVAR.covmat_fit[ISN_SAVE][ipar][ipar2] ;
      fprintf(FP_STDOUT, "%10.3le ", covfit);
    }
    fprintf(FP_STDOUT, "  |  ");
    for(ipar2=0; ipar2<NLCPAR; ipar2++ ) {
      covtot = INFO_DATA.TABLEVAR.covmat_tot[ISN_SAVE][ipar][ipar2] ;
      fprintf(FP_STDOUT, "%10.3le ", covtot );
    }

    fprintf(FP_STDOUT, "\n");
  }

  fflush(FP_STDOUT);

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
  // Sep 3 2021: REPEAT fit msg incldues alpha,beta,redchi2.
  // Feb 22 2022: fix bug to recalc_dataCov if sigint=0

  double redchi2, covParam ;
  double step1 = INPUTS.covint_param_step1 ;
  double COVINT_PARAM_MIN = 0.0 ;
  int STOP_TOL, STOP_MXFIT, STOP_COV0, STOP_COVFIX, retCode, USE_CCPRIOR ;
  int NFIT_ITER = FITRESULT.NFIT_ITER ;
  char msg[100];
  char fnam[] = "prepNextFit" ;

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
  STOP_COV0   = ( NFIT_ITER > 0 && 
		  FITINP.COVINT_PARAM_FIX <= COVINT_PARAM_MIN ) ;
  
  STOP_COVFIX = ( strlen(INPUTS.sigint_fix) > 0);

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

  
  if ( STOP_TOL || STOP_MXFIT || STOP_COV0 || STOP_COVFIX ) {

    retCode = FITFLAG_DONE ; 
    fprintf(FP_STDOUT, "\t Final %s value = %0.3f  for chi2(Ia)/dof=%.4f\n",
	   FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	   FITINP.COVINT_PARAM_FIX, redchi2 );
  } 

  else if ( INPUTS.fitflag_sigmb == 0 && USE_CCPRIOR > 0 ) {
    retCode = FITFLAG_DONE ; 
  } 
  else {
    // Try another covParam 
    // On 2nd iteration, use linear approx and dchi2red/dsigint 
    // to estimate next covParam

    covParam = FITINP.COVINT_PARAM_FIX ;
    FITINP.COVINT_PARAM_FIX = next_covFitPar(redchi2,covParam,step1); 
    if ( FITINP.COVINT_PARAM_FIX < COVINT_PARAM_MIN ) 
      {  FITINP.COVINT_PARAM_FIX = COVINT_PARAM_MIN ; } 
    recalc_dataCov();

    retCode = FITFLAG_CHI2 ;

    /* 
    sprintf(msg,"last redchi2=%.3f --> next %s=%.3f (NCALL_FCN=%d)", 
	    redchi2, FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	    FITINP.COVINT_PARAM_FIX,  FITRESULT.NCALL_FCN );
     */
    sprintf(msg,"last a,b,redchi2=%.3f,%.2f,%.3f --> next %s=%.3f", 
	    FITRESULT.ALPHA, FITRESULT.BETA, redchi2, 
	    FITRESULT.PARNAME[IPAR_COVINT_PARAM],
	    FITINP.COVINT_PARAM_FIX );
    printmsg_repeatFit(msg);
  }
  
  return(retCode) ;

} // end of prepNextFit



// ******************************************
void printmsg_repeatFit(char *msg) {
  fprintf(FP_STDOUT, "\n");
  fprintf(FP_STDOUT, "# %s\n", dotDashLine);
  fprintf(FP_STDOUT, "    REPEAT fit: %s\n", msg);
  fprintf(FP_STDOUT, "# %s\n", dotDashLine);
  fprintf(FP_STDOUT, "\n");
  fflush(FP_STDOUT);
} 


#ifdef USE_THREAD
void fcn(int *npar, double grad[], double *fval, double xval[],
	 int *iflag, void *not_used) {

  int  NSN_DATA    = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  nthread     = INPUTS.nthread ;
  int  NFITPAR_ALL = FITINP.NFITPAR_ALL ; // Ncospar + Nzbin
  int  ipar, t, rc, NERR, NSN_per_thread, isn_min, isn_max ;

  pthread_t thread[MXTHREAD];
  thread_chi2sums_def  thread_chi2sums[MXTHREAD];
  char fnam[] = "fnam";

  // ----------- BEGIN ----------------

  FITRESULT.NCALL_FCN++ ;

  // bail on inf or nan.
  for(ipar=1; ipar<=5; ipar++ ) {
    if ( isnan(xval[ipar]) ) { *fval = 1.0E14; return; }
    if ( isinf(xval[ipar]) ) { *fval = 1.0E14; return; }
  }

  if ( nthread == 1 ) 
    { NSN_per_thread = NSN_DATA; }
  else
    { NSN_per_thread = (int)( (float)NSN_DATA/(float)nthread )  + 1 ; } 

  // - - - - - - - - - - - - - - - - - - -
  for ( t = 0; t < nthread; t++ ) {

    isn_min = t    * NSN_per_thread;
    isn_max = (t+1)* NSN_per_thread;
    if ( isn_max > NSN_DATA ) { isn_max = NSN_DATA; }

    thread_chi2sums[t].nthread   = nthread;
    thread_chi2sums[t].id_thread = t ;
    thread_chi2sums[t].isn_min   = isn_min ;
    thread_chi2sums[t].isn_max   = isn_max ;
    

    // load fcn args to typedef struct
    thread_chi2sums[t].npar_fcn  = *npar ;
    thread_chi2sums[t].iflag_fcn = *iflag ;
    for(ipar=0; ipar < NFITPAR_ALL ; ipar++ ) 
      { thread_chi2sums[t].xval_fcn[ipar] = xval[ipar];   }

    if ( nthread == 1 )
      {  MNCHI2FUN(&thread_chi2sums);  } 
#ifdef USE_THREAD
    else  { 
      rc = pthread_create(&thread[t], NULL, MNCHI2FUN, 
			  &thread_chi2sums[t] ) ; 
    }
#endif
  }  // end t loop over threads

  // - - - - - - - 
#ifdef USE_THREAD
  // for threads, wait for them all to finish
  if ( nthread > 1 ) {
    NERR = 0 ;
    for ( t = 0; t < nthread; t++ ) { 
      rc = pthread_join(thread[t], NULL); 
      if ( rc != 0 ) {
	NERR++; 
	printf(" ERROR: thread return errcode=%d for t=%d\n", rc,t); }
    }

    if ( NERR > 0 ) {
      sprintf(c1err,"%d thread return code errors", NERR);
      sprintf(c2err,"NCALL_FCN=%d", FITRESULT.NCALL_FCN );
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
    }  
  } // end ntrhread>1
#endif

  // ===============================================
  // ============= WRAP UP =========================
  // ===============================================

  // sum each thread
  int nsnfit = 0, nsnfit_truecc=0;
  double chi2sum_Ia=0.0, chi2sum_tot=0.0, nsnfitIa=0.0, nsnfitcc=0.0 ;
  //  double alpha, beta, gamma, logmass  ;
  for ( t = 0; t < nthread; t++ ) { 
    nsnfit        += thread_chi2sums[t].nsnfit ;
    nsnfit_truecc += thread_chi2sums[t].nsnfit_truecc ;
    nsnfitIa      += thread_chi2sums[t].nsnfitIa ;
    nsnfitcc      += thread_chi2sums[t].nsnfitcc ;
    chi2sum_Ia    += thread_chi2sums[t].chi2sum_Ia ;
    chi2sum_tot   += thread_chi2sums[t].chi2sum_tot ;
  }


  // load globals
  FITRESULT.NSNFIT        = nsnfit ;
  FITRESULT.NSNFIT_TRUECC = nsnfit_truecc ;
  FITRESULT.NSNFIT_SPLITRAN[NJOB_SPLITRAN] = nsnfit ;
    
  if ( *iflag == 3 )  {   // done with fit
    double xdof = nsnfitIa - (double)FITINP.NFITPAR_FLOAT ;
    FITRESULT.CHI2SUM_1A = chi2sum_Ia ;
    FITRESULT.CHI2RED_1A = chi2sum_Ia/xdof ;
    FITRESULT.NSNFIT_1A  = nsnfitIa ; 
    FITRESULT.NSNFIT_CC  = nsnfitcc ; 

    // a,b,g stored for re-computing COV between fit iterations
    FITRESULT.ALPHA      = xval[IPAR_ALPHA0];
    FITRESULT.BETA       = xval[IPAR_BETA0];
    FITRESULT.GAMMA      = xval[IPAR_GAMMA0];
  }
  
  *fval = chi2sum_tot;

  return ;
    
} // end fcn for pthread

// =================================================================
void *MNCHI2FUN(void *thread) {

  // Created Aug 31 2020
  // Essentially this is fcn, but called from wrapper than has
  // pthread option to distribute chi2 loop over many cores.
  // The data loop starts at id_thread, and skips in steps of
  // nthread. For default nthread=1 (no thread), it is a normal 
  // loop over all events from 0 to NSN_DATA-1 
  //
  // Apr 8 2021: subtract muerr_vpec from muerr_raw
  // Sep 24 2021: abort on muerrsq < 0
  // Sep 27 2021: require muCOVadd>0 to implement; fixes rare muerrsq<0 problem.

  thread_chi2sums_def *thread_chi2sums = (thread_chi2sums_def *)thread;
  //  int  npar      = thread_chi2sums->npar_fcn ;
  int  iflag     = thread_chi2sums->iflag_fcn ;
  double *xval   = thread_chi2sums->xval_fcn ;
  int  nthread   = thread_chi2sums->nthread;
  int  id_thread = thread_chi2sums->id_thread ;
  int  isn_min   = thread_chi2sums->isn_min ;
  int  isn_max   = thread_chi2sums->isn_max ;
  char fnam[]    = "MNCHI2FUN" ;
  char *name ;

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD  ) > 0;
  bool APPLY_COVADD ;

  int NDIM_BIASCOR, INTERPFLAG_abg;
  double logmass, omega_l, omega_k, wde, wa;
  double alpha, beta, gamma, scalePROB_fitpar, *hostPar ;
  double cosPar[NCOSPAR], ProbRatio_Ia, ProbRatio_CC ;
  double chi2sum_tot, chi2sum_Ia, sqmures, mures, mu, muBias, muBiasErr ;
  double z, zmuerr, mb, x1, c, mumodel_store, muCOVscale, muCOVadd ;
  double PTOTRAW_Ia=0.0, PTOTRAW_CC=0.0, muBias_zinterp;
  double muerr, muerrsq, muerrsq_last, muerr_last ;
  double PSUM,Prob_SUM,PTOT_Ia,PTOT_CC, dPdmu_Ia, dPdmu_CC, Prob_Ia, Prob_CC ;
  double sqsigCC=0.001, sigCC_chi2penalty=0.0;
  double covmat_tot[NLCPAR][NLCPAR], covmat_fit[NLCPAR][NLCPAR] ;
  double gammaDM, M0, dl, mumodel ;
  double muerr_raw, muerrsq_raw, muerr_vpec, muerrsq_vpec;
  double muerrsq_tmp, muerr_update, muerrsq_update ; 
  double chi2evt, chi2evt_Ia, scalePIa, scalePCC, nsnfitIa=0.0, nsnfitcc=0.0;
  int    n, nsnfit, nsnfit_truecc, ipar, ipar2 ;
  int    cutmask, idsample, SIM_NONIA_INDEX, IS_SIM ; 
  int    ia, ib, ig, optmask_muerrsq ;
  int    dumpFlag_muerrsq=0, DUMPFLAG=0 ;
  int    USE_CCPRIOR=0, USE_CCPRIOR_H11=0 ;
  bool   set_fitwgt0 = false ;
  MUZMAP_DEF  *CCPRIOR_MUZMAP ;

  int  ILCPAR_MIN = INFO_BIASCOR.ILCPAR_MIN ;
  int  ILCPAR_MAX = INFO_BIASCOR.ILCPAR_MAX ;  

  BIASCORLIST_DEF     BIASCORLIST ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
  FITPARBIAS_DEF FITPARBIAS_ALPHABETA[MXa][MXb][MXg]; // bias at each a,b
  double   MUCOVSCALE_ALPHABETA[MXa][MXb][MXg]; // (I) muCOVscale at each a,b
  double   MUCOVADD_ALPHABETA[MXa][MXb][MXg]; // (I) muCOVadd at each a,b
  double   *fitParBias;

  // -------------- BEGIN ------------

  //Set input cosmology parameters
  //  alpha0       = xval[IPAR_ALPHA0] ;
  //  beta0        = xval[IPAR_BETA0] ;
  //  da_dz        = xval[3] ;
  //  db_dz        = xval[4] ;
  //  gamma0       = xval[IPAR_GAMMA0] ;
  //  dg_dz        = xval[6] ;
  //  logmass_cen  = xval[7] ; 
  //  logmass_tau  = xval[8] ;
  omega_l      = xval[IPAR_OL] ;
  omega_k      = xval[IPAR_Ok] ;
  wde          = xval[IPAR_w0] ;
  wa           = xval[IPAR_wa] ;
  scalePROB_fitpar  = xval[IPAR_scalePCC] ;
  hostPar      = &xval[IPAR_GAMMA0];
  //  da_dm        = xval[15];  // added Apr 2 2018
  //  db_dm        = xval[16];  // idem

  NDIM_BIASCOR = INFO_BIASCOR.NDIM;

  // load cosPar array to pass to functions below
  cosPar[0] = omega_l ;
  cosPar[1] = omega_k ;
  cosPar[2] = wde ;
  cosPar[3] = wa ;

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

  chi2sum_tot = chi2sum_Ia    = 0.0;
  nsnfit      = nsnfit_truecc = 0 ;
  nsnfitIa    = nsnfitcc      = 0.0 ;

  // - - - - - - - - - - - - - - - - -
  for ( n = isn_min; n < isn_max; n++ ) {

    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[n] ; 
    if ( cutmask ) { continue; }

    // - - - - -

    INFO_DATA.mures[n]     = -999. ;
    INFO_DATA.mupull[n]    = -999. ;
    INFO_DATA.mu[n]        = -999. ;
    INFO_DATA.muerr[n]     = -999. ;    
    INFO_DATA.muerr_raw[n] = -999. ;  // no scale and no sigInt
    INFO_DATA.muerr_vpec[n] = -999. ;  // muerr from vpec only
    name     = INFO_DATA.TABLEVAR.name[n] ;
    idsample = (int)INFO_DATA.TABLEVAR.IDSAMPLE[n] ;
    z        = (double)INFO_DATA.TABLEVAR.zhd[n] ;     
    logmass  = (double)INFO_DATA.TABLEVAR.host_logmass[n];
    zmuerr   = (double)INFO_DATA.TABLEVAR.zmuerr[n] ; // for muerr calc
    mb       = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_mB][n] ;
    x1       = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_x1][n] ;
    c        = (double)INFO_DATA.TABLEVAR.fitpar[INDEX_c][n] ;
    mumodel_store   = (double)INFO_DATA.TABLEVAR.mumodel[n] ; 

    SIM_NONIA_INDEX = (int)INFO_DATA.TABLEVAR.SIM_NONIA_INDEX[n];
    IS_SIM          = (INFO_DATA.TABLEVAR.IS_SIM == true);
    
    if ( USE_CCPRIOR ) { 
      PTOTRAW_Ia  = (double)INFO_DATA.TABLEVAR.pIa[n] ; 
    }

    muerrsq_last   = INFO_DATA.muerrsq_last[n] ;
    muerr_last     = INFO_DATA.muerr_last[n] ;

    for(ipar=0; ipar < NLCPAR; ipar++ ) {
      for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) 
	{ covmat_tot[ipar][ipar2] = 
	    INFO_DATA.TABLEVAR.covmat_tot[n][ipar][ipar2] ; 
	}
    }

    if ( NDIM_BIASCOR > 0 ) {
      muBias_zinterp = INFO_DATA.muBias_zinterp[n] ; 
      set_fitwgt0    = INFO_DATA.set_fitwgt0[n];
            
      for(ia=0; ia<MXa; ia++ ) {
	for(ib=0; ib<MXb; ib++ ) {
	  for(ig=0; ig<MXg; ig++ ) {
	    
	    MUCOVSCALE_ALPHABETA[ia][ib][ig] = 
	      INFO_DATA.MUCOVSCALE_ALPHABETA[n][ia][ib][ig] ; 

	    if ( DO_COVADD ) {
	      MUCOVADD_ALPHABETA[ia][ib][ig] = 
		INFO_DATA.MUCOVADD_ALPHABETA[n][ia][ib][ig] ; 
	    }

	    if (set_fitwgt0) { MUCOVSCALE_ALPHABETA[ia][ib][ig]=1.0; }
	    if (set_fitwgt0) { MUCOVADD_ALPHABETA[ia][ib][ig]=1.0; }
	    
	    for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {
	      FITPARBIAS_ALPHABETA[ia][ib][ig].VAL[ipar] = 
		INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig].VAL[ipar] ; 
	      FITPARBIAS_ALPHABETA[ia][ib][ig].ERR[ipar] = 
		INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig].ERR[ipar] ; 
	    } // end ipar
	  }  // end ig
	}    // end ib 
      }     // end ia
      fitParBias = INFO_DATA.fitParBias[n] ; 
    } // end NDIM_BIASCOR if block

    if ( z < 1.0E-8) { continue ; } // Jun 3 2013 (obsolete?)

    // fetch alpha,beta,gamma (include z-dependence)
    fcnFetch_AlphaBetaGamma(xval, z, logmass, &alpha, &beta, &gamma); 

    gammaDM = get_gammadm_host(z, logmass, hostPar );

    DUMPFLAG = 0 ; // ( strcmp(name,"93018")==0 ) ; 
    if ( INTERPFLAG_abg ) {
      get_INTERPWGT_abg(alpha, beta, gammaDM, DUMPFLAG, &INTERPWGT, name);
    }

    DUMPFLAG = 0 ;

    // get mag offset for this z-bin
    M0    = fcn_M0(n, &xval[MXCOSPAR] );

    // compute distance modulus from cosmology params
    if ( INPUTS.FLOAT_COSPAR ) {
      // not tested, so beware !!!
      dl       = cosmodl_forFit(z, z, cosPar) ;
      mumodel  = 5.0*log10(dl) + 25.0 ;
    }
    else {
      // intended use for BBC
      mumodel = mumodel_store ; 
    }


    optmask_muerrsq = 0;
    if ( set_fitwgt0 ) { optmask_muerrsq += DOFLAG_CUTWIN_FITWGT0; }

    // compute error-squared on distance mod
    muerrsq  = fcn_muerrsq(name, alpha, beta, gamma, covmat_tot,
			   z, zmuerr, optmask_muerrsq );

    // --------------------------------
    // Compute bias from biasCor sample
    BIASCORLIST.z            = z;
    BIASCORLIST.host_logmass = logmass;
    BIASCORLIST.alpha        = alpha ;
    BIASCORLIST.beta         = beta ;
    BIASCORLIST.gammadm      = gammaDM ;
    BIASCORLIST.idsample     = idsample;
    BIASCORLIST.FITPAR[INDEX_mB] = mb ;
    BIASCORLIST.FITPAR[INDEX_x1] = x1 ;
    BIASCORLIST.FITPAR[INDEX_c]  = c ;
    muBias = muBiasErr = 0.0 ;  
    muCOVscale =1.0 ;
    muCOVadd = 0.0;

    if ( NDIM_BIASCOR >= 5 ) {
      get_muBias(name, &BIASCORLIST,      // (I) misc inputs
		 FITPARBIAS_ALPHABETA,    // (I) bias at each a,b,g
		 MUCOVSCALE_ALPHABETA,    // (I) muCOVscale at each a,b
		 MUCOVADD_ALPHABETA,    // (I) muCOVscale at each a,b
		 &INTERPWGT,              // (I) wgt at each a,b,g grid point
		 fitParBias,     // (O) interp bias on mB,x1,c
		 &muBias,        // (O) interp bias on mu
		 &muBiasErr,     // (O) stat-error on above
		 &muCOVscale,   // (O) scale bias on muCOV  
		 &muCOVadd );   // (O) add bias on muCOV     
    }
    else if ( NDIM_BIASCOR == 1 ) {
      muBias  = muBias_zinterp ; 
    }
       
    if ( n == -95 ) { 
      fprintf(FP_STDOUT," xxx %s-%3.3d %s  a,b,g=%.4f,%.4f,%.4f  "
	     "muBias=%7.4f  (fpb0=%.4f,%.4f) \n",
	     fnam, FITRESULT.NCALL_FCN, name, alpha,beta,gamma, 
	     muBias,
	     FITPARBIAS_ALPHABETA[0][0][0].VAL[ILCPAR_MAX],
	     FITPARBIAS_ALPHABETA[0][1][0].VAL[ILCPAR_MAX]
	     );
      fflush(FP_STDOUT);
    }

    // load muBias info into globals

    if ( NDIM_BIASCOR ) {
      INFO_DATA.muBias[n]     = muBias ;
      INFO_DATA.muBiasErr[n]  = muBiasErr ;
      INFO_DATA.muCOVscale[n] = muCOVscale ;
      if ( DO_COVADD ) {
	INFO_DATA.muCOVadd[n] = muCOVadd ;
      }
    }

    // zero out muBiasErr after storing it, since adding this
    // would contradict the muCOVscale correction.
    muBiasErr = 0.0 ; 

    APPLY_COVADD = ( DO_COVADD && muCOVadd > 0.0  );
                 
    bool restore_mucovadd_bug =(INPUTS.restore_mucovadd_bug&2)>0;
    if (restore_mucovadd_bug){
      APPLY_COVADD = ( DO_COVADD && muCOVscale > 1.0  );}
    
    if ( APPLY_COVADD ) {
      // Aug 2 2021: Dillon's sigint in bins. note that global sigint = 0
      muerrsq += muCOVadd; 
      
      if ( muerrsq < 0.0 ) {
	double muerrsq_orig = muerrsq - muCOVadd;
	print_preAbort_banner(fnam);
	printf("   SNID=%s  IDSAMPLE=%d \n", name, idsample );
	printf("   z=%.5f, mB=%.3f  x1=%.4f  c=%.4f  logmass=%.3f \n",
	       z, mb, x1, c, logmass);
	printf("   alpha=%.4f  beta=%.4f  gDM=%.4f\n", 
	       alpha, beta, gammaDM);
	dump_muCOVcorr(n);
	
	sprintf(c1err,"Insane muerrsq = %.5f for snid=%s (muerrsq_orig=%.5f)", 
		muerrsq, muerrsq_orig, name);
	sprintf(c2err,"muCOV[scale,add] = %.5f,%.5f ", 
		muCOVscale, muCOVadd );
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
      }

    } else {
      // Scale as in original BBC
      muerr_vpec    = fcn_muerrz(1, z, zmuerr); 
      muerrsq_vpec  = muerr_vpec * muerr_vpec ;
      muerrsq       = (muerrsq-muerrsq_vpec) * muCOVscale + muerrsq_vpec ;
    }
    
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
      nsnfit++ ;      nsnfitIa  = (double)nsnfit ;
      chi2evt_Ia    = sqmures/muerrsq ;
      chi2sum_Ia   += chi2evt_Ia ;
      chi2evt       = chi2evt_Ia ;

      // check option to add log(sigma) term for 5D biasCor
      if ( INPUTS.fitflag_sigmb == 2 ) 
      	{ chi2evt  += log(muerrsq/muerrsq_last); } 

      if ( iflag == 3 ) { INFO_DATA.probcc_beams[n] = 0.0 ; } // Dec 2020
    } 
    
    if ( USE_CCPRIOR  ) {
      // BEAMS-like chi2 = -2ln [ PIa + PCC ]
      DUMPFLAG = 0 ; // (n == 45 && FITRESULT.NCALL_FCN < 3700);
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
	fprintf(FP_STDOUT, " xxx ---------------------------------- \n");
	fprintf(FP_STDOUT, " xxx fcn dump for SN = '%s'  (NCALL=%d) \n", 
		name, FITRESULT.NCALL_FCN );
	fprintf(FP_STDOUT, " xxx dPdmu_CC=%le \n", dPdmu_CC);
	fprintf(FP_STDOUT, " xxx PTOT_CC = scale*PTOTRAW/PSUM  = "
		"%le*%le/%le = %le\n",
	       scalePCC, PTOTRAW_CC, PSUM, PTOT_CC );
	fprintf(FP_STDOUT, " xxx Prob(Ia,CC) = %le, %le \n", Prob_Ia, Prob_CC);
	fflush(FP_STDOUT);
	//	if ( FITRESULT.NCALL_FCN > 3000 )  { debugexit(fnam); }
      }
      // xxxxxxxxxxxxxxxxxxx

      
      // sum total prob that includes Ia + CC
      if ( Prob_SUM > 0.0 ) {
	ProbRatio_Ia = Prob_Ia / Prob_SUM ;
	ProbRatio_CC = Prob_CC / Prob_SUM ;
	nsnfitIa    +=  ProbRatio_Ia ; 
	nsnfitcc    +=  ProbRatio_CC ; 
	chi2sum_Ia  += (ProbRatio_Ia * chi2evt_Ia) ; 
	
	if ( iflag == 3 ) { INFO_DATA.probcc_beams[n] = ProbRatio_CC;}

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
    if (  iflag==3 ) {	

      // Jan 26 2018: store raw muerr without intrinsic cov
      for(ipar=0; ipar < NLCPAR; ipar++ ) {
	for(ipar2=0; ipar2 < NLCPAR; ipar2++ ) {	  
	  covmat_fit[ipar][ipar2] =
	    (double)INFO_DATA.TABLEVAR.covmat_fit[n][ipar][ipar2];
	}
      }

      muerrsq_raw = fcn_muerrsq(name,alpha,beta,gamma,covmat_fit, z,zmuerr,
				optmask_muerrsq );

      // compute diagnostic muerr_raw(from LC fit, no vpec, no sigint)
      // and muerr_vpec (from vpec only). Neither include COVscale.
      muerr_vpec      = fcn_muerrz(1, z, zmuerr); 
      muerrsq_vpec    = muerr_vpec * muerr_vpec;
      muerr_raw       = sqrt(muerrsq_raw - muerrsq_vpec);
      INFO_DATA.muerr_raw[n]  = muerr_raw ;
      INFO_DATA.muerr_vpec[n] = muerr_vpec ;

      // check user dump with total error (Jun 19 2018)
      dumpFlag_muerrsq = ( strcmp(name,INPUTS.SNID_MUCOVDUMP) == 0 );
      if ( dumpFlag_muerrsq ) {
	int optmask_dump = optmask_muerrsq + 1 ;
	muerrsq_tmp = fcn_muerrsq(name,alpha,beta,gamma,covmat_fit,z,zmuerr,
				  optmask_dump );
	muerrsq_tmp = fcn_muerrsq(name,alpha,beta,gamma,covmat_tot,z,zmuerr,
				  optmask_dump );
      }
      
      if ( IS_SIM && SIM_NONIA_INDEX > 0 ) { nsnfit_truecc++ ; } 

      // store reference errors for 1/sigma term
      if ( DOFIT_FLAG == FITFLAG_CHI2 ) {
	INFO_DATA.muerr_last[n]   = muerr;
	INFO_DATA.muerrsq_last[n] = muerrsq;
	INFO_DATA.sigCC_last[n]   = sqrt(sqsigCC) ;
	INFO_DATA.sqsigCC_last[n] = sqsigCC ;
      }
    }  // end iflag==3 
    
    // ----------------------------------
	  
  } // end loop over SN

  // - - - - - - - - - - - - - - - - - - - - 
  // load sums in output typedef  
 
  thread_chi2sums->nsnfit        = nsnfit ;
  thread_chi2sums->nsnfit_truecc = nsnfit_truecc ;
  thread_chi2sums->nsnfitIa      = nsnfitIa ;
  thread_chi2sums->nsnfitcc      = nsnfitcc ;

  thread_chi2sums->chi2sum_tot   = chi2sum_tot ;
  thread_chi2sums->chi2sum_Ia    = chi2sum_Ia  ;

  // check CPU-load balance on first FCN call
  if ( FITRESULT.NCALL_FCN == 1 && nthread > 1 ) {
    printf("\t %s-%3.3d: id_thread = %d of %d  nsnfit=%d \n", 
	   fnam, FITRESULT.NCALL_FCN, id_thread, nthread, nsnfit );
    fflush(stdout);
  }

  return(void *) 0 ;

} // end MNCHI2FUN

#endif 


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

  M0      = M0_DEFAULT ;
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
    }

    // if interp z-bin does not float M0, then just return
    // constant M0 in this z-bin (Jun 2017)
    if ( FITINP.ISFLOAT_z[iz1] == 0  ) { return(M0); }

    zbin1  = ptr_zM0[iz1];
    M0bin0 = M0LIST[iz0];
    M0bin1 = M0LIST[iz1];

    zfrac = ( zdata - zbin0 ) / ( zbin1 - zbin0) ;
    M0    = M0bin0 + (M0bin1-M0bin0) * zfrac ;

    LDMP = (n == -95 ); // xxx REMOVE
    if ( LDMP ) {    
      fprintf(FP_STDOUT, " xxx -------------------------- \n");
      fprintf(FP_STDOUT, " xxx iz0=%d  iz1=%d \n", iz0, iz1);
      fprintf(FP_STDOUT, " xxx zdata=%.4f  zbin0=%.4f zbin1=%.4f  zfrac=%.3f\n",
	     zdata, zbin0, zbin1, zfrac );
      fprintf(FP_STDOUT, " xxx M0bin0=%.3f  M0bin1=%.3f   M0=%.3f\n",
	     M0bin0, M0bin1, M0);
      fflush(FP_STDOUT);
      //      debugexit(fnam);
    }

  }
  else {
    sprintf(c1err,"Invalid uM0=%d", INPUTS.uM0 );
    sprintf(c2err,"data index n=%d", n);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
  }

  return ;
  
} // end get_INTERPWGT_abg


// ===========================================================
double fcn_muerrsq(char *name, double alpha, double beta, double gamma,
		   double (*COV)[NLCPAR], 
		   double z, double zerr, int optmask ) {

  // Created Feb 29 2016 by R.Kessler
  //
  // return mu-error for inputs as follows:
  //
  // *name      : name of SN (in case of error msg)
  // alpha,beta : SALT2 standardization parmas for stretch & color
  // gamma      : SALT2 standard param for logmass (NOT USED)
  // COV        : fit cov matrix + intrinsic cov matrix
  // z,zerr     : redshift & vpec error
  // 
  // optmask += 1 -> dump flag
  // optmask += 2 -> return MUERR_FITWGT0^2
  //
  // Jun 27 2017: REFACTOR z bins
  // Sep 04 2017: use matrix-vector summing:
  //    MUERRSQ = (1,a,b) x COV x (1,a,b)
  //
  // Jun 19 2018: pass dumpFlag argument.
  // Jan 22 2021: change dumpFlag arg to optmask arg; check MUERR_WGT option

  bool   flag_dump    = (optmask & 1) > 0 ;
  bool   flag_fitwgt0 = (optmask & 2) > 0 ;

  double muerrsq, sqtmp, muerr_z, muerrsq_z, dmuLens;
  double VEC[NLCPAR], COVMU[NLCPAR][NLCPAR];
  int    i,j ;
  char VECNAME[3][8] = { "ONE  ", "ALPHA", "-BETA" } ;
  
  char fnam[] = "fcn_muerrsq" ;

  // ---------------- BEGIN ------------

  if ( flag_dump ) {
    printf(" xxx \n");
    printf(" xxx --- %s dump for SNID=%s ----- \n", fnam, name);
    printf(" xxx alpha=%le  beta=%le  z=%le   zerr=%le \n",
	   alpha, beta, z, zerr);
    fflush(stdout);
  }
  
  
  if ( flag_fitwgt0 ) { return(MUERR_FITWGT0*MUERR_FITWGT0) ; }

  VEC[0] = +1.0;
  VEC[1] = +alpha ;
  VEC[2] = -beta; 

  muerrsq = 0.0 ;
  for(i=0; i<NLCPAR; i++ ) {
    for(j=0; j<NLCPAR; j++ ) {
      sqtmp       = VEC[j] * COV[i][j] * VEC[i] ;
      COVMU[i][j] = sqtmp;  // for diagnostic only
      muerrsq    += sqtmp ;

      if( flag_dump ) {
	printf(" xxx mucov += %13.6le (= %s * %s * %13.6le) \n",
	       sqtmp, VECNAME[j], VECNAME[i], COV[i][j] );
	fflush(stdout);
      }
    }
  }

  if( flag_dump ) { printf(" xxx mucov  = %le from COV \n", muerrsq);  }
  
  // add peculiar velocity error
  if ( zerr > 0.0 ) {
    muerr_z     = fcn_muerrz(1, z, zerr );
    muerrsq_z   = muerr_z*muerr_z;
    muerrsq    += muerrsq_z ;
  }
  else {
    muerr_z = muerrsq_z = 0.0 ;
  }

  if( flag_dump ) {
    printf(" xxx mucov  = %le with dmu(z,vpec)=%le \n",
			muerrsq, muerr_z);   fflush(stdout);
  }
    
  // add lensing error (Sep 9 2016)
  dmuLens  = INPUTS.lensing_zpar * z;
  muerrsq += (dmuLens * dmuLens);

  if( flag_dump ) {
    printf(" xxx mucov  = %le with dmuLens=%le \n",
	   muerrsq, dmuLens) ;
    printf(" xxx mucov(FINAL) = %le (does not include biasScale_muCOV)\n",
	   muerrsq);
    fflush(stdout);
  }

    
  if ( muerrsq < muerrsq_z )  {
    print_preAbort_banner(fnam);

    printf("   COV(SALT2) = \n" );
    for(i=0; i<NLCPAR; i++ ) {      
      for(j=0; j<NLCPAR; j++ )
	{ printf(" %12.6f ", COV[i][j] ); }
      printf("\n");
    }

    printf("   COV(MU) = \n" );
    for(i=0; i<NLCPAR; i++ ) {      
      for(j=0; j<NLCPAR; j++ )
	{ printf(" %12.6f ", COVMU[i][j] ); }
      printf("\n");
    }

    printf("   muerr_z  = %f (zerr=%f)\n", muerr_z, zerr);
    printf("   muerrsq(noVpec) = %le \n", muerrsq - muerrsq_z);
    printf("   dmuLens  = %f \n", dmuLens);
    printf("   alpha = %f  beta=%f \n", alpha, beta);
    sprintf(c1err,"muerrsq = %le is less than muerrsq_z = %le", 
	    muerrsq, muerrsq_z );	
    sprintf(c2err,"for SN = %s ", name);   
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

  double zerrtot, zlo, zhi, muerr = 0.0 ;
  double FAC   = 5.0/LOGTEN ;  
  double *cosPar = &INPUTS.parval[IPAR_OL] ;
  //  char fnam[]  = "fcn_muerrz" ;

  // --------------- BEGIN --------------

  if ( OPT == 1 ) {
    zerrtot = zerr ;
    muerr    = FAC*(zerrtot/z) * (1.0+z)/(1.0+z/2.0);
  }
  else if ( OPT == 2 )  {

    double dl, mu1, mu2 ;
    zlo    = z-zerr;  zhi = z+zerr;

    dl     = cosmodl(zlo,zlo,cosPar);
    mu1    = 5.0*log10(dl) + 25.0 ;

    dl     = cosmodl(zhi,zhi,cosPar);
    mu2    = 5.0*log10(dl) + 25.0 ;
    muerr  = (mu2-mu1)/2.0 ;
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

  double zerr_new      = zerr;
  double zpecerr_orig  = vpecerr/LIGHT_km; // from data file
  double zpecerr_user  = INPUTS.zpecerr;   // override with this user value
  double z1            = (1+z) ;
  double sqz1          = z1*z1 ;
  double sqzerr1, sqzerr2, sqzerr3, sqarg;
  double zdif_tol = 1.0E-5 ;
  char fnam[] = "zerr_adjust" ;
  
  // ---------- BEGIN ------------

  if ( zpecerr_user < 1.0E-9 ) { return(zerr); }

  sqzerr1 = zerr * zerr;
  sqzerr2 = sqz1 * zpecerr_orig * zpecerr_orig ; // original 
  sqzerr3 = sqz1 * zpecerr_user * zpecerr_user ; // user-defined

  // allow tolerance on difference to allow for
  // precision truncation on FITRES file.
  if ( z1*zpecerr_orig > zerr+zdif_tol ) {
    print_preAbort_banner(fnam);
    printf("    z1*zpecerr = %f * %f = %f \n", 
	   z1, zpecerr_orig, z1*zpecerr_orig);
    printf("    zerr = %f \n", zerr);
    sprintf(c1err,"Cannot subtract zpecerr from zerr for SNID=%s", name );
    sprintf(c2err,"zerr=%f > (1+z)*zpecerr = (%f)*%f", 
	    zerr, 1+z, zpecerr_orig);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
  }

  sqarg    = (sqzerr1 - sqzerr2 + sqzerr3 );

  if ( sqarg < 0.0 ) {
    print_preAbort_banner(fnam);
    printf("   sqzerr[1,2,3] = %le, %le, %le \n",
	   sqzerr1, sqzerr2, sqzerr3 );
    printf("   zpecerr_user = %f \n", zpecerr_user);
    sprintf(c1err,"sqrt(negative number); Cannot adjust zHDERR for CID=%s", name);
    sprintf(c2err,"Check pre-abort dump above.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }
  zerr_new = sqrt(sqarg );

  return(zerr_new) ;

} // end zerr_adjust 


// ******************************************
void set_defaults(void) {

  int isurvey, order, ipar, N ;
  char fnam[] = "set_defaults";

  // ---------------- BEGIN -----------------
  
  set_EXIT_ERRCODE(EXIT_ERRCODE_SALT2mu);
  N = store_PARSE_WORDS(-1,""); // May 26 2021                                               
  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );

  INPUTS.cat_only   = false ;
  INPUTS.cat_prescale = 1;
  INPUTS.cat_file_out[0] = 0 ;
  INPUTS.write_yaml = 0 ;
  INPUTS.write_csv  = 0 ;

  INPUTS.minos      = 0 ; // disable default minos, Apr 22 2022
  INPUTS.nfile_data = 0 ;
  INPUTS.nfile_data_override = 0 ;
  sprintf(INPUTS.PREFIX,     "NONE" );
  INPUTS.KEYNAME_DUMPFLAG      = false ;

  INPUTS.opt_biasCor           = 0 ;
  INPUTS.sigint_biasCor        = -9.0 ; 
  INPUTS.snrmin_sigint_biasCor = BIASCOR_SNRMIN_SIGINT ;

  INPUTS.prescale_biasCor[0]   = 0;  // subset number
  INPUTS.prescale_biasCor[1]   = 1;  // integer pre-scale

  INPUTS.sigma_cell_biasCor    = 50.0; // flat WGT distribution in each cell
  INPUTS.nfile_biasCor         = 0 ;
  sprintf(INPUTS.fieldGroup_biasCor, "NONE" );
  INPUTS.use_fieldGroup_biasCor = 0 ;

  INPUTS.ndump_nobiasCor  = 0 ;
  INPUTS.frac_warn_nobiasCor = 0.02 ;

  INPUTS.surveyGroup_biasCor_abortFlag = 1 ;
  sprintf(INPUTS.surveyGroup_biasCor, "NONE" );
  INPUTS.use_surveyGroup_biasCor = 0 ;

  sprintf(INPUTS.surveyList_noBiasCor, "NONE" );
  INPUTS.idsample_select[0] = 0 ;
  INPUTS.zspec_errmax_idsample = 0.0 ;

  INPUTS.select_trueIa  = 0;
  INPUTS.force_realdata = 0 ;
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
  INPUTS.sameFile_flag_CCprior = false ;

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

  INPUTS.zmin = -0.01 ;
  INPUTS.zmax =  9.99 ;

  // ---------------------
  // keep obsolete input parameters (4/24/2012 RK)
  // Mar 2022: update to use nominal SALT2 cuts
  INPUTS.x1min = -3.0;
  INPUTS.x1max = +3.0;
  INPUTS.cmin  = -0.3; 
  INPUTS.cmax  = +0.3;

  INPUTS.LEGACY_NBINc = false ;
  if ( INPUTS.LEGACY_NBINc ) {
    printf("\n xxx %s: BEWARE LEGACY_NBINc = True !! \n\n", fnam );
    INPUTS.x1min = -6.0;
    INPUTS.x1max = +6.0;
    INPUTS.cmin  = -0.6; 
    INPUTS.cmax  = +0.6;    
  }

  INPUTS.logmass_min  = -20.0 ;
  INPUTS.logmass_max  = +20.0 ;
  INPUTS.nbin_logmass =  1 ; 

  INPUTS.chi2max       = 1.0E9 ;
  INPUTS.iflag_chi2max = 0;       // Dec 2020

  // ---------------------

  INPUTS.cutmask_write = 0;  // default -> write events used in fit

  INPUTS.Nsntype = 0 ;
  sprintf(INPUTS.sntypeString,"NULL");

  INPUTS.nmaxString[0]   = 0 ;

  INPUTS.ncidFile_data   = 0;
  INPUTS.ncidList_data   = 0;
  INPUTS.acceptFlag_cidFile_data = 0 ;
  INPUTS.izbin_from_cidFile = 0;

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

  INPUTS.zwin_vpec_check[0] = 0.01 ;  // Oct 28 2020
  INPUTS.zwin_vpec_check[1] = 0.05 ; 

  INPUTS.fitflag_sigmb       = 0;     // option to repeat fit until chi2/dof=1
  INPUTS.redchi2_tol         = 0.02;  // tolerance on chi2.dof
  INPUTS.sigint_step1        = 0.01 ; // size of 1st sigint step 
  INPUTS.dchi2red_dsigint    = 0.0 ;
  INPUTS.scale_covint_step1  = 0.04 ; // for scale_covint

  INPUTS.prescale_simData  = 1.0 ; // include all simData by default
  INPUTS.prescale_simCC    = 1.0 ; // include all simcc by default.
  INPUTS.prescale_simIa    = 1.0 ; // include all simIa by default.

  DOFIT_FLAG = FITFLAG_CHI2 ;
  INPUTS.zpolyflag = 0;
  //Default distance modulus parameters
  INPUTS.H0      = 70.0;
  INPUTS.nommag0 = M0_DEFAULT ;
  INPUTS.uave    = 1;

  // Default parameters 

  // set IPAR_XXX params to access FITRESULTS.
  IPAR_ALPHA0=1, IPAR_BETA0=2;
  IPAR_GAMMA0=5; IPAR_GAMMA1=6;  IPAR_LOGMASS_CEN=7; IPAR_LOGMASS_TAU=8;
  IPAR_scalePCC=13, IPAR_COVINT_PARAM=14 ;
  IPAR_OL=9; IPAR_Ok=10; IPAR_w0=11;  IPAR_wa=12;
  IPAR_H11=17; NPAR_H11_TOT=6; NPAR_H11_USER=0;

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
  INPUTS.debug_flag        = 0 ;
  INPUTS.debug_malloc      = 0 ;
  INPUTS.debug_mucovscale  = -9 ; // negative to avoid i1d dump
  INPUTS.nbinc_mucovscale  = 3 ;
  INPUTS.restore_sigz      = 0 ; // 0->new, 1->old(legacy)
  INPUTS.restore_mucovscale_bug = 0 ;
  INPUTS.restore_mucovadd_bug = 0 ;
  INPUTS.nthread           = 1 ; // 1 -> no thread

  INPUTS.cidlist_debug_biascor[0] = 0 ;

  // === set blind-par values to be used if blindflag=2 (Aug 2017)
  ISDATA_REAL = 1 ;
  INPUTS.blind_cosinePar[IPAR_OL][0] = 0.06; 
  INPUTS.blind_cosinePar[IPAR_OL][1] = 23434. ;

  INPUTS.blind_cosinePar[IPAR_w0][0] = 0.20 ; 
  INPUTS.blind_cosinePar[IPAR_w0][1] = 8430. ;

  sprintf(INPUTS.varname_gamma,VARNAME_LOGMASS);
  INPUTS.USE_GAMMA0  = 0 ;

  INPUTS.LCUTWIN_DISABLE = false;
  init_CUTMASK();

#ifdef USE_SUBPROCESS
  SUBPROCESS.USE         = false ;
  SUBPROCESS.ITER        = -9;
  SUBPROCESS.INPUT_ISEED = 12345;
  SUBPROCESS.STDOUT_CLOBBER  = 1; // default is to clobber each stdout
  SUBPROCESS.NEVT_SIM_PRESCALE     = -9; 
#endif

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

  for(bit=0; bit < MXNUM_SAMPLE; bit++ ) {
    NDATA_BIASCORCUT[0][bit] = NDATA_BIASCORCUT[1][bit] = 0 ; 
  }


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
  sprintf(CUTSTRING_LIST[CUTBIT_SPLITRAN],  "SPLITRAN-subset"); 
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
  // 
  // May 12 2021: call read_data_override here before set_CUTMASK is called.

  int  NFILE = INPUTS.nfile_data;
  //  int  EVENT_TYPE = EVENT_TYPE_DATA;
  int  NEVT_TOT, NEVT[MXFILE_DATA], ifile, NVAR_ORIG, IFILETYPE, LEN_MALLOC ;
  int  ISTART, NROW, isn ;
  char *dataFile ;
  char fnam[] = "read_data" ;

  // -------------- BEGIN --------------

   sprintf(BANNER,"%s", fnam);
   fprint_banner(FP_STDOUT,BANNER);   

   // read each file quickly to get total size of arrays to malloc
   NEVT_TOT = 0 ;
   for(ifile=0; ifile < NFILE; ifile++ ) {
     dataFile     = INPUTS.dataFile[ifile];
     NEVT[ifile]  = SNTABLE_NEVT(dataFile,TABLENAME_FITRES);  
     NEVT_TOT    += NEVT[ifile];
     fprintf(FP_STDOUT, "\t Found %d events in %s. \n", NEVT[ifile], dataFile);
     fflush(FP_STDOUT);
   }  
   

   // malloc arrays to read fitres data file
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

  // compute more table variables
  for(isn=0; isn < INFO_DATA.TABLEVAR.NSN_ALL; isn++ ) { 
    compute_more_TABLEVAR(isn, &INFO_DATA.TABLEVAR ); 
  }

  // Nov 2020: e.g., replace VPEC, HOST_LOGMASS, etc ..
  read_data_override(); 

  // apply cuts after data override
  int NPASS=0;
  for(isn=0; isn < INFO_DATA.TABLEVAR.NSN_ALL; isn++ ) { 
    compute_CUTMASK(isn, &INFO_DATA.TABLEVAR ); 
    if ( INFO_DATA.TABLEVAR.CUTMASK[isn] == 0 ) { NPASS++ ; }
  }

  if ( NPASS == 0 ) {
    print_preAbort_banner(fnam);
    print_eventStats(EVENT_TYPE_DATA);
    sprintf(c1err,"All DATA events fail cuts");
    sprintf(c2err,"Check cut-windows in SALT2mu input file.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
  }


  return;

} // end read_data 


// ****************************************            
float malloc_MUCOV(int opt, int IDSAMPLE, CELLINFO_DEF *CELLINFO ) {

  // Created July 28 2021 by Dillon Brout
  // Compute CELLINFO bins and malloc arrays.
  // Used later for adjusting muErr based on RMS from sim.
  // Input opt not used,.
  // Sep 14 2021: malloc USE element

  int debug_malloc = INPUTS.debug_malloc ;
  bool DO_MAD      = (INPUTS.opt_biasCor & MASK_BIASCOR_MAD) > 0;
  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  float f_MEMORY = 0.0;
  int NCELL;

  int NBINa,NBINb,NBINg,NBINm,NBINz,NBINc;
  double cmin,cmax,cbin,c_lo,c_hi,c_avg;
  double mmin,mmax,mbin,m_lo,m_hi,m_avg;
  int ic,im,i1d;

  int NPERCELL_REALLOC=2000;
  char fnam[] = "malloc_MUCOV";

  // ------------- BEGIN --------------

  NBINa    = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  NBINb    = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;  
  NBINg    = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;  

  //  ptr_MUCOVSCALE = INFO_BIASCOR.MUCOVSCALE[IDSAMPLE];   

  // redshift bins are same as for biasCor
  copy_BININFO(&CELLINFO_BIASCOR[IDSAMPLE].BININFO_z, 
	       &CELLINFO->BININFO_z);
  NBINz = CELLINFO->BININFO_z.nbin ;

  // setup color bins that are coarser than those for muBias.
  // Note that other bins (a,b,z) are same for muCOVscale and muBias.
  NBINc=INPUTS.nbinc_mucovscale; // user input instead of hard-wired (March 2022 MVincenzi)
  cmin=INPUTS.cmin; cmax=INPUTS.cmax;
  cbin=(cmax-cmin)/(double)NBINc;

  if ( INPUTS.LEGACY_NBINc ) { NBINc=3; cmin = -0.3; cmax=0.3; cbin=0.2; }

  // Jun 11 2021: include logmass bins.
  NBINm = CELLINFO_BIASCOR[IDSAMPLE].BININFO_m.nbin ;
  if ( NBINm == 1 ) 
    { NBINm = 1; mmin=-15.0; mmax=15.0; mbin=30.0 ; } // legacy behavior
  else  { 
    // start with 2 logmass bins for devel; later decide on either
    // biasCor logmass bins, or hard-wired bins such as for color
    NBINm = 2; mmin=5.0; mmax=15.0; mbin=5.0 ;  // refactor
  }

  // - - - - 
  NCELL    = NBINa * NBINb * NBINg * NBINz * NBINc * NBINm ;
  CELLINFO->NCELL = NCELL ; 

  // ---------------------------------------------
  // malloc global struct to store CELLINFO
  fprintf(FP_STDOUT, "\t Malloc CELL-INFO arrays with size=%d "
	  "(IDSAMPLE=%d)\n", NCELL, IDSAMPLE);

  int MEMD     = NCELL   * sizeof(double);
  int MEMI     = NCELL   * sizeof(int);
  int MEMB     = NCELL   * sizeof(bool);

  print_debug_malloc(+1*debug_malloc,fnam);

  // allocate redshift bins
  CELLINFO->USE      =  (bool   *) malloc(MEMB);
  CELLINFO->NperCell =  (int    *) malloc(MEMI);
  CELLINFO->AVG_z    =  (double *) malloc(MEMD);
  CELLINFO->AVG_m    =  (double *) malloc(MEMD);
  CELLINFO->AVG_LCFIT[INDEX_c] = (double *) malloc(MEMD);

  // set color bins; and note that color bins are different ...
  CELLINFO->BININFO_LCFIT[INDEX_c].nbin    = NBINc ;
  CELLINFO->BININFO_LCFIT[INDEX_c].binSize = cbin ;
  for(ic=0; ic < NBINc; ic++ ) {
    c_lo  = cmin + cbin * (double)ic ;
    c_hi  = c_lo + cbin ;
    c_avg = 0.5*(c_lo + c_hi) ;
    CELLINFO->BININFO_LCFIT[INDEX_c].lo[ic]  = c_lo ;
    CELLINFO->BININFO_LCFIT[INDEX_c].hi[ic]  = c_hi ;
    CELLINFO->BININFO_LCFIT[INDEX_c].avg[ic] = c_avg ;
  }

  
  // June 11 2021: set logmass bins
  CELLINFO->BININFO_m.nbin    = NBINm ;
  CELLINFO->BININFO_m.binSize = mbin ;
  for(im=0; im < NBINm; im++ ) {
    m_lo  = mmin + mbin * (double)im ;
    m_hi  = m_lo + mbin ;
    m_avg = 0.5*(m_lo + m_hi) ;
    CELLINFO->BININFO_m.lo[im]  = m_lo ;
    CELLINFO->BININFO_m.hi[im]  = m_hi ;
    CELLINFO->BININFO_m.avg[im] = m_avg ;
  }

  if (DO_MAD) {
    if (DO_COVSCALE || DO_COVADD) {
	CELLINFO->ABSPULL =  (double **) malloc(sizeof(double*)*NCELL);
	for (i1d=0; i1d<NCELL; i1d++){
	  CELLINFO->ABSPULL[i1d] = (double *) malloc(sizeof(double)*NPERCELL_REALLOC);
	}
    }
    if (DO_COVADD) {
      CELLINFO->MURES =  (double **) malloc(sizeof(double*)*NCELL);
      CELLINFO->MUCOV =  (double **) malloc(sizeof(double*)*NCELL);
      for (i1d=0; i1d<NCELL; i1d++){
	CELLINFO->MURES[i1d] = (double *) malloc(sizeof(double)*NPERCELL_REALLOC);
	CELLINFO->MUCOV[i1d] = (double *) malloc(sizeof(double)*NPERCELL_REALLOC);
      }
    }
  }



  return f_MEMORY;
} // end malloc_MUCOV


// ****************************************
void malloc_INFO_DATA(int opt, int LEN_MALLOC ) {

  // Created June 5 2019
  // Malloc or free INFO_DATA strucure.
  // opt > 0 --> malloc
  // opt < 0 --> free
  //
  // Jan 22 2021: malloc set_fitwgt0
  // Jun 29 2021: check nfile_biasCor for biascor-dependent mallocs

  int nfile_biasCor = INPUTS.nfile_biasCor ;
  int EVENT_TYPE    = EVENT_TYPE_DATA;
  int debug_malloc = INPUTS.debug_malloc ;
  int MEMD, MEMI, MEMB, i_mem, N_MEM=0;
  long long MEMTOT=0 ;
  float f_MEMORY = 0.0 ;
  float f_MEM[10];
  char  COMMENT_MEM[10][80];
  char fnam[] = "malloc_INFO_DATA";

  // ------------- BEGIN --------------

  print_debug_malloc(opt*debug_malloc,fnam);

  if ( opt > 0 ) {

    // start with generic malloc to read any FITRES file
    f_MEM[N_MEM] = malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_DATA.TABLEVAR);
    sprintf(COMMENT_MEM[N_MEM], "TABLEVAR");	  N_MEM++;

    CUTMASK_POINTER[EVENT_TYPE]         = &INFO_DATA.TABLEVAR.CUTMASK[0];
    IDSAMPLE_POINTER[EVENT_TYPE]        = &INFO_DATA.TABLEVAR.IDSAMPLE[0];
    NALL_CUTMASK_POINTER[EVENT_TYPE]    = &INFO_DATA.TABLEVAR.NSN_ALL;
    NPASS_CUTMASK_POINTER[EVENT_TYPE]   = &INFO_DATA.TABLEVAR.NSN_PASSCUTS;
    NREJECT_CUTMASK_POINTER[EVENT_TYPE] = &INFO_DATA.TABLEVAR.NSN_REJECT;

    // DATA-specific mallocs ...
    MEMD  = LEN_MALLOC * sizeof(double);
    MEMI  = LEN_MALLOC * sizeof(int);
    MEMB  = LEN_MALLOC * sizeof(bool);

    INFO_DATA.mumodel        = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.M0             = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mu             = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerr          = (double*) malloc(MEMD); MEMTOT+=MEMD; 
    INFO_DATA.muerr_raw      = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerr_vpec     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mures          = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.mupull         = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerr_last     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.muerrsq_last   = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.sigCC_last     = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.sqsigCC_last   = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.chi2           = (double*) malloc(MEMD); MEMTOT+=MEMD;
    INFO_DATA.probcc_beams   = (double*) malloc(MEMD); MEMTOT+=MEMD;

    if ( nfile_biasCor > 0 ) {
      INFO_DATA.muCOVscale     = (double*) malloc(MEMD); MEMTOT+=MEMD;
      INFO_DATA.muCOVadd       = (double*) malloc(MEMD); MEMTOT+=MEMD;
      INFO_DATA.muBias         = (double*) malloc(MEMD); MEMTOT+=MEMD;
      INFO_DATA.muBiasErr      = (double*) malloc(MEMD); MEMTOT+=MEMD;
      INFO_DATA.muBias_zinterp = (double*) malloc(MEMD); MEMTOT+=MEMD;
      INFO_DATA.set_fitwgt0    = (bool  *) malloc(MEMB); MEMTOT+=MEMB;
    }

    f_MEM[N_MEM] = (float)(MEMTOT/1.0E6);
    sprintf(COMMENT_MEM[N_MEM], "INFO_DATA");  N_MEM++;

    if ( nfile_biasCor > 0 ) {
      // fitParBias[isn][ipar] to allow passing mB,x1,c via array
      f_MEM[N_MEM] = malloc_double2D(opt, LEN_MALLOC, NLCPAR+1, 
				     &INFO_DATA.fitParBias ); // <== returned
      sprintf(COMMENT_MEM[N_MEM], "INFO_DATA.fitparBias");  N_MEM++;

      f_MEM[N_MEM] = 
	malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
				    &INFO_DATA.FITPARBIAS_ALPHABETA );
      sprintf(COMMENT_MEM[N_MEM], "INFO_DATA.FITPARBIAS_AB");  N_MEM++;

      // MUCOVSCALE_AB[isn][ia][ib][ig]
      f_MEM[N_MEM] = 
	malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg,
			&INFO_DATA.MUCOVSCALE_ALPHABETA ); //<==return
      sprintf(COMMENT_MEM[N_MEM], "INFO_DATA.MUCOVSCALE");  N_MEM++; 

      f_MEM[N_MEM] = 
	malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg,
			&INFO_DATA.MUCOVADD_ALPHABETA ); //<==return
      sprintf(COMMENT_MEM[N_MEM], "INFO_DATA.MUCOVADD");  N_MEM++; 

      f_MEM[N_MEM] = 
	malloc_shortint4D(opt, LEN_MALLOC, MXa, MXb, MXg,
			  &INFO_DATA.I1D_MUCOVSCALE ); //<==return
      sprintf(COMMENT_MEM[N_MEM], "INFO_DATA.I1D_MUCOVSCALE");  N_MEM++; 

    }

    for(i_mem=0; i_mem < N_MEM; i_mem++ ) {
      f_MEMORY += f_MEM[i_mem];
      fprintf(FP_STDOUT, "\t %s: allocate %7.3f MB for %s\n", 
	      fnam, f_MEM[i_mem], COMMENT_MEM[i_mem] );
    }

    INFO_DATA.MEMORY = f_MEMORY; 
    fprintf(FP_STDOUT, "\t %s:   TOTAL  %7.3f MB \n\n", fnam, f_MEMORY); 
    fflush(FP_STDOUT);
  }
  else {
    malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_DATA.TABLEVAR);
    free(INFO_DATA.mumodel);
    free(INFO_DATA.M0);
    free(INFO_DATA.mu);
    free(INFO_DATA.muerr);
    free(INFO_DATA.muerr_raw);
    free(INFO_DATA.muerr_vpec);
    free(INFO_DATA.mures);
    free(INFO_DATA.mupull);
    free(INFO_DATA.muerr_last);
    free(INFO_DATA.muerrsq_last);
    free(INFO_DATA.sigCC_last);
    free(INFO_DATA.sqsigCC_last);
    free(INFO_DATA.muBias);
    free(INFO_DATA.muBiasErr);
    free(INFO_DATA.muBias_zinterp);
    free(INFO_DATA.chi2);
    free(INFO_DATA.probcc_beams);
    free(INFO_DATA.set_fitwgt0);

    malloc_double2D(opt, LEN_MALLOC, NLCPAR+1, &INFO_DATA.fitParBias ); 
    malloc_FITPARBIAS_ALPHABETA(opt, LEN_MALLOC,
				&INFO_DATA.FITPARBIAS_ALPHABETA );
    malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg, 
		    &INFO_DATA.MUCOVSCALE_ALPHABETA ); 
    malloc_double4D(opt, LEN_MALLOC, MXa, MXb, MXg, 
		    &INFO_DATA.MUCOVADD_ALPHABETA ); 

  }


  return ;

} // end malloc_INFO_DATA


// ****************************************
void read_data_override(void) {

  // Created Nov 13 2020
  // Read optional data_override file, and replace values
  // in data arrays. For VPEC and VPEC_ERR, also modify
  // zHD and zHDERR, respectively. Similarly, for zHEL
  // override, update zCMB and zHD.
  //
  // May 12 2021: for HOST_LOGMASS, also set CUTVAL for applying cuts.
  // Nov 18 2021: 
  //   + allow VPECERR or VPEC_ERR
  //   + avoid double-counting zHD override if zHEL and VPEC are both
  //      on override list. Same for zHDERR with ZHELERR and VPECERR.
  //

  int IVAR_OVER_VPEC = -9, IVAR_OVER_VPECERR = -9 ;
  int IVAR_OVER_zHEL = -9, IVAR_OVER_zHELERR = -9 ;
  int IVAR_OVER_zHD  = -9, IVAR_OVER_zHDERR  = -9 ;
  int IVAR_OVER_zCMB = -9, IVAR_OVER_LOGMASS = -9 ;
  int NSN_CHANGE[MXVAR_OVERRIDE];
  int nfile_over   = INPUTS.nfile_data_override;
  int IVAR_GAMMA   = INFO_DATA.TABLEVAR.ICUTWIN_GAMMA ;
  int debug_malloc = INPUTS.debug_malloc ;
  int ifile_data, ifile_over, NROW;
  int ivar_data, NVAR_DATA, ivar_over, NVAR_OVER, OPTMASK, ntmp;
  char *ptrFile, *varName, *VARNAMES_STRING_DATA, *VARNAMES_STRING_OVER ;
  char fnam[] = "read_data_override" ;

  // ------------- BEGIN --------------

  INFO_DATA.NVAR_OVERRIDE = 0;
  if ( nfile_over == 0 ) { return; }

  print_banner(fnam);  printf("\n");

  // prepare comma sep list of all varNames in first data files
  // (later may need to check all data files?)
  print_debug_malloc(+1*debug_malloc,fnam);
  VARNAMES_STRING_DATA = (char*) malloc(MXCHAR_VARLIST*sizeof(char));
  VARNAMES_STRING_DATA[0] = 0;
  ifile_data = 0 ; // primary file index
  NVAR_DATA = INFO_DATA.TABLEVAR.NVAR[ifile_data];
  for(ivar_data=0; ivar_data < NVAR_DATA; ivar_data++ ) { 
      varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile_data][ivar_data];
      catVarList_with_comma(VARNAMES_STRING_DATA,varName);
  }


  // read in all of the override files and store variables
  // that appear in the data
  OPTMASK = 4; // append each file
  for(ifile_over=0; ifile_over < nfile_over; ifile_over++ ) {
    ptrFile = INPUTS.dataFile_override[ifile_over];
    ENVreplace(ptrFile,fnam,1);
    NROW = SNTABLE_AUTOSTORE_INIT(ptrFile,"OVERRIDE", 
				  VARNAMES_STRING_DATA, OPTMASK);
  }

  // check each data varName and store those which exist in 
  // an override file

  VARNAMES_STRING_OVER    = (char*) malloc(MXCHAR_VARLIST*sizeof(char));
  VARNAMES_STRING_OVER[0] = NVAR_OVER = 0 ;
  for(ivar_data=0; ivar_data < NVAR_DATA; ivar_data++ ) { 
    varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[0][ivar_data];
    if ( EXIST_VARNAME_AUTOSTORE(varName) )  { 
      catVarList_with_comma(VARNAMES_STRING_OVER,varName); 

      if ( strcmp(varName,VARNAME_VPEC) == 0     ) 
	{ IVAR_OVER_VPEC = NVAR_OVER ; }

      if ( strcmp(varName,VARNAME_VPECERR) == 0 ) 
	{ IVAR_OVER_VPECERR = NVAR_OVER ; }
      if ( strcmp(varName,VARNAME_VPECERR2) == 0 ) // alternate key 11.18.2021
	{ IVAR_OVER_VPECERR = NVAR_OVER ; }

      if ( strcmp(varName,VARNAME_zHEL) == 0     ) 
	{ IVAR_OVER_zHEL = NVAR_OVER ; }

      if ( strcmp(varName,VARNAME_zHELERR) == 0     ) 
	{ IVAR_OVER_zHELERR = NVAR_OVER ; }

      if ( strcmp(varName,VARNAME_zCMB) == 0     ) 
	{ IVAR_OVER_zCMB = NVAR_OVER ; }

      if ( strcmp(varName,VARNAME_LOGMASS) == 0     ) 
	{ IVAR_OVER_LOGMASS = NVAR_OVER ; }

      NVAR_OVER++ ;
    }
  }


  // - - - - - - - -
  // check for conflicts

  if ( NVAR_OVER == 0 ) {
    sprintf(c1err,"Found zero override variables.");
    sprintf(c2err,"Check VARNAMES in datafile_override");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }
  

  if ( IVAR_OVER_zHD >= 0 && IVAR_OVER_zHEL >= 0 ) {
    sprintf(c1err,"Cannot override both zHD and zHEL; pick one");
    sprintf(c2err,"and SALT2mu auto-computes override of the other.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( IVAR_OVER_zHD >= 0 && IVAR_OVER_VPEC >= 0 ) {
    sprintf(c1err,"Cannot override both zHD and VPEC; pick one");
    sprintf(c2err,"and SALT2mu auto-computes override of the other.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( IVAR_OVER_zCMB >= 0 ) {
    sprintf(c1err,"Cannot override zCMB.");
    sprintf(c2err,"Try overrid for zHEL or zHD");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( IVAR_OVER_zHELERR >= 0 ) {
    sprintf(c1err,"Sorry, zHELERR override not implemented yet.");
    sprintf(c2err,"Please post Github issue if this is important.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  // - - - - - - 
  // if VPEC or zHEL is on override list, add recalc zHD to override list.

  if ( IVAR_OVER_VPEC >= 0 || IVAR_OVER_zHEL >= 0 ) {     
    IVAR_OVER_zHD = NVAR_OVER ;
    NVAR_OVER++; catVarList_with_comma(VARNAMES_STRING_OVER,VARNAME_zHD); 
  }

  // if zHELERR or VPEC_ERR is on override list, add recalculated zHDERR
  if ( IVAR_OVER_VPECERR >= 0 || IVAR_OVER_zHELERR >= 0 )  {
    IVAR_OVER_zHDERR = NVAR_OVER ; 
    NVAR_OVER++;   catVarList_with_comma(VARNAMES_STRING_OVER,VARNAME_zHDERR); 
  }


  INFO_DATA.NVAR_OVERRIDE = NVAR_OVER;

  if ( NVAR_OVER >= MXVAR_OVERRIDE )  {
    sprintf(c1err,"NVAR_OVERRIDE=%d exceeds bound.", NVAR_OVER);
    sprintf(c2err,"Check MXVAR_OVERRIDE=%d parameter", MXVAR_OVERRIDE);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  /*
  printf(" xxx %s: IVAR_OVER_VPEC[ERR] = %d, %d \n", 
	 fnam, IVAR_OVER_VPEC, IVAR_OVER_VPECERR);
  printf(" xxx %s: IVAR_OVER_ZHD[ERR] = %d, %d \n", 
	 fnam, IVAR_OVER_ZHD, IVAR_OVER_ZHDERR);
  */

  // - - - - - -
  // convert comma sep list to array
  parse_commaSepList("OVERRIDE", VARNAMES_STRING_OVER, 
		     MXVAR_OVERRIDE, MXCHAR_VARNAME, &ntmp, 
		     &INFO_DATA.VARNAMES_OVERRIDE);

  // - - - - - - -
  // assign pointer to float array for each override variable
  for(ivar_over=0; ivar_over < NVAR_OVER; ivar_over++ ) {
    varName   = INFO_DATA.VARNAMES_OVERRIDE[ivar_over] ;
    NSN_CHANGE[ivar_over] = 0 ;

    if ( strcmp(varName,VARNAME_VPEC) == 0 )
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.vpec ;  }   

    // Nov 18 2021: allow VPECERR or VPEC_ERR
    else if ( strcmp(varName,VARNAME_VPECERR) == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.vpecerr; }  
    else if ( strcmp(varName,VARNAME_VPECERR2) == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.vpecerr; }  

    else if ( strcmp(varName,"HOST_LOGMASS") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = 
	  INFO_DATA.TABLEVAR.host_logmass ;  }

    else if ( strcmp(varName,"zHD") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zhd ;  }
    else if ( strcmp(varName,"zHDERR") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zhderr ;  }

    else if ( strcmp(varName,"zCMB") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zcmb ;  }
    else if ( strcmp(varName,"zCMBERR") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zcmberr ; }

    else if ( strcmp(varName,"zHEL") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zhel ;  }
    else if ( strcmp(varName,"zHELERR") == 0 ) 
      { INFO_DATA.PTRVAL_OVERRIDE[ivar_over] = INFO_DATA.TABLEVAR.zhelerr ;  }

    else {
      sprintf(c1err,"Unable to implement override for %s", varName);
      sprintf(c2err,"Might need to update function %s", fnam );
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
  }  // end ivar_over loop

  // - - - - - - - - -  
  // loop over each data event and each varname to override, 
  // and replace value.

  int NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;
  int istat, isn;
  bool  override_zhd, override_zhderr;
  float *fval ;   double dval;    char *name, *cval ;
  double zhd_over, zhderr_over, dl, zhel_over, zhd_orig, zhel_orig ;

  for(isn=0; isn < NSN_DATA; isn++ ) { 

    zhd_orig   = (double)INFO_DATA.TABLEVAR.zhd[isn];
    zhel_orig  = (double)INFO_DATA.TABLEVAR.zhel[isn];
    override_zhd = override_zhderr = false;

    for(ivar_over=0; ivar_over < NVAR_OVER; ivar_over++ ) {

      if ( ivar_over == IVAR_OVER_zHD    ) { continue ; }
      if ( ivar_over == IVAR_OVER_zHDERR ) { continue ; }

      name    = INFO_DATA.TABLEVAR.name[isn];
      varName = INFO_DATA.VARNAMES_OVERRIDE[ivar_over];
      SNTABLE_AUTOSTORE_READ(name, varName, &istat, &dval, cval ); 
      
      // xxx mark delete override_zhd = false;

      // xxxxxxx
      if ( istat == -99990 ) {
	printf(" xxx %s: isn=%4d ivar_over=%d  istat=%d, dval=%f \n",
	       fnam, isn, ivar_over, istat, dval); fflush(stdout);
      }// xxxx
     
      if ( istat == 0 ) {

	// check computed zhd changes based on changes to vpec[err]
	if ( ivar_over == IVAR_OVER_VPEC ) {
	  double vpec_over = dval;
	  zhd_over = zhd_data_override(isn,vpec_over); 
	  INFO_DATA.PTRVAL_OVERRIDE[IVAR_OVER_zHD][isn] = zhd_over;
	  if(!override_zhd) { NSN_CHANGE[IVAR_OVER_zHD]++ ; }
	  override_zhd = true ;
	}
	else if ( ivar_over == IVAR_OVER_VPECERR ) {	 
	  double vpecerr_over = dval;
	  zhderr_over = zhderr_data_override(isn,vpecerr_over); 
	  INFO_DATA.PTRVAL_OVERRIDE[IVAR_OVER_zHDERR][isn] = zhderr_over ;
	  INFO_DATA.TABLEVAR.zmuerr[isn] = vpecerr_over/LIGHT_km ;// Erik Peterson 11/19
	  NSN_CHANGE[IVAR_OVER_zHDERR]++ ; 
          override_zhderr = true ;
	}
	else if ( ivar_over == IVAR_OVER_zHEL ) {
	  zhel_over = dval ;
	  zhd_over  = zhd_orig + (zhel_over-zhel_orig); 
	  INFO_DATA.PTRVAL_OVERRIDE[IVAR_OVER_zHD][isn] = zhd_over;
	  if(!override_zhd) { NSN_CHANGE[IVAR_OVER_zHD]++ ; }
	  override_zhd = true ;
	}
	else if ( ivar_over == IVAR_OVER_LOGMASS && IVAR_GAMMA >= 0 ) {
	  double logmass = dval;
	  INFO_DATA.TABLEVAR.CUTVAL[IVAR_GAMMA][isn] = logmass ;
	}

	// apply override AFTER checking zhd[err] overrides
	INFO_DATA.PTRVAL_OVERRIDE[ivar_over][isn] = dval;
	NSN_CHANGE[ivar_over]++ ;

	// if zhd is modified, update zhel & mumodel.
	if ( override_zhd ) {
	  zhel_over = zhel_orig + (zhd_over - zhd_orig);
	  dl = cosmodl_forFit(zhel_over,zhd_over,INPUTS.COSPAR); 
	  INFO_DATA.TABLEVAR.mumodel[isn] = (float)(5.0*log10(dl) + 25.0);
	}
       
      } // end istat if-block

    } // end ivar_over
  } // end isn 


  // - - - - - - 
  // Finally, map ivar_out to ivar_over so that output function can
  // quickly identify which output variables get replaced.
  int ivar_out, NVAR_TOT = OUTPUT_VARNAMES.NVAR_TOT ;  
  INFO_DATA.IVAR_OUTPUT_INVMAP = (int*) malloc(NVAR_TOT*sizeof(int));
  for(ivar_out=0; ivar_out < NVAR_TOT; ivar_out++ ) {
    varName   = OUTPUT_VARNAMES.LIST[ivar_out] ; 
    ivar_over = ivar_matchList(varName, NVAR_OVER, 
			       INFO_DATA.VARNAMES_OVERRIDE );   
    INFO_DATA.IVAR_OUTPUT_INVMAP[ivar_out] = ivar_over ;
  }

  // - - - - -
  // print summary of changes
  for(ivar_over=0; ivar_over < NVAR_OVER; ivar_over++ ) {
    varName   = INFO_DATA.VARNAMES_OVERRIDE[ivar_over] ;
    printf("   Override %-20.20s for %d of %d events \n",
	   varName, NSN_CHANGE[ivar_over], NSN_DATA ); 
    fflush(stdout);
  }


  print_debug_malloc(-1*debug_malloc,fnam);
  free(VARNAMES_STRING_DATA);
  free(VARNAMES_STRING_OVER);

  //  debugexit(fnam); // xxx remove
  return ;

} // end  read_data_override


// *********************************
double zhd_data_override(int isn, double vpec_over ) {

  // Created Nov 13 2020
  // For input isn index and vpec(override), return new redshift 
  // zhd_over to override old value.
  //
  // zHD_orig + 1  = (1+zCMB)/(1+zpec_orig) where zCMB is computed from zHEL.
  //  -->
  // zHD_over + 1 = (zHD_orig+1)*(1+zpec_orig)/(1+zpec_over)

  double zhd_orig  = (double)INFO_DATA.TABLEVAR.zhd[isn];
  double vpec_orig = (double)INFO_DATA.TABLEVAR.vpec[isn];

  double zpec_orig = vpec_orig / LIGHT_km;
  double zpec_over = vpec_over / LIGHT_km;
  double zhd_over;
  char fnam[] = "zhd_data_override" ;

  // -------- BEGIN ---------
  
  // debug with approx formula; later replace with exact formula
  zhd_over = (zhd_orig+1.0) * ( 1.0+zpec_orig) / (1+zpec_over) - 1.0 ;
  return zhd_over ;

} // end zhd_data_override

 
double zhderr_data_override(int isn, double vpecerr_over ) {

  // Nov 15 2020
  // Return zhderr with orginal vpecerr_orig replaced with vpecerr_over
  // From snana.car:
  //     ZERR1 = SNLC_ZCMB_ERR
  //	 ZERR2 = ZPECERR * (1.0 + SNLC_ZCMB) ! Eq A1, Davis 2012
  //     SNLC_ZHD_ERR = sqrt(ZERR1*ZERR1 + ZERR2*ZERR2)
  //
 
  double zhderr_orig  = (double)INFO_DATA.TABLEVAR.zhderr[isn];
  double vpecerr_orig = (double)INFO_DATA.TABLEVAR.vpecerr[isn];
  double zcmb         = (double)INFO_DATA.TABLEVAR.zcmb[isn];
  double zcmberr      = (double)INFO_DATA.TABLEVAR.zcmberr[isn];
  double zpecerr_orig = vpecerr_orig / LIGHT_km ;
  double zpecerr_over = vpecerr_over / LIGHT_km ;

  double ZERR1 = zcmberr;
  double ZERR2 = zpecerr_over * ( 1.0 + zcmb );

  double zhderr_over, zhderrsq ;
  char fnam[] = "zhderr_data_override" ;

  //---------- BEGIN ------------

  zhderrsq = ZERR1*ZERR1 + ZERR2*ZERR2 ;
  zhderr_over = sqrt(zhderrsq);
  return zhderr_over ;

} // end zhderr_data_override


// ****************************************
void malloc_INFO_BIASCOR(int opt, int LEN_MALLOC ) {

  // Created June 5 2019
  // Malloc or free INFO_BIASCOR structure
  // opt > 0 --> malloc
  // opt < 0 --> free

  int  NSAMPLE    = NSAMPLE_BIASCOR;
  int  EVENT_TYPE = EVENT_TYPE_BIASCOR ;
  int  MEMB       = LEN_MALLOC * sizeof(int8_t);     // one byte
  int  MEMI       = NSAMPLE * sizeof(int);
  int  MEMF       = NSAMPLE * sizeof(float);
  int  MEMBIAS    = NSAMPLE * sizeof(FITPARBIAS_DEF*);
  int  MEMCOV     = NSAMPLE * sizeof(float*);
  int debug_malloc = INPUTS.debug_malloc ;

  int  MEMADD=0 ;
  float f_MEMORY = 0.0 ;
  char fnam[] = "malloc_INFO_BIASCOR";

  // ------------- BEGIN --------------

  print_debug_malloc(opt*debug_malloc,fnam);

  if ( opt > 0 ) {
    
    // start with generic malloc to read any FITRES file
    f_MEMORY = malloc_TABLEVAR(opt, LEN_MALLOC, &INFO_BIASCOR.TABLEVAR);

    // set pointers to track event stats vs. EVENT_TYPE
    CUTMASK_POINTER[EVENT_TYPE]         = &INFO_BIASCOR.TABLEVAR.CUTMASK[0];
    IDSAMPLE_POINTER[EVENT_TYPE]        = &INFO_BIASCOR.TABLEVAR.IDSAMPLE[0];
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
    INFO_BIASCOR.MUCOVADD   = (float **        ) malloc ( MEMCOV  );

    // print memory consumption to stdout
    f_MEMORY += ((float)MEMADD) / 1.0E6 ;
    INFO_BIASCOR.MEMORY = f_MEMORY;
    fprintf(FP_STDOUT, "\t %s: %6.3f MB \n", fnam, f_MEMORY); 
    fflush(FP_STDOUT);

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
  int debug_malloc = INPUTS.debug_malloc ;
  int MEMTOT=0 ;
  float f_MEMORY = 0.0 ;
  char fnam[] = "malloc_INFO_CCPRIOR";

  // ------------- BEGIN --------------

  print_debug_malloc(opt*debug_malloc,fnam);

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
    fprintf(FP_STDOUT, "\t %s: %6.3f MB \n", fnam, f_MEMORY); 
    fflush(FP_STDOUT);

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
  //
  // Jun 30 2021: fix bug summing memory from malloc_float3D; see f_MEM

  int EVENT_TYPE = TABLEVAR->EVENT_TYPE;
  int IS_DATA    = (EVENT_TYPE == EVENT_TYPE_DATA);
  //int IS_BIASCOR = (EVENT_TYPE == EVENT_TYPE_BIASCOR);
  //  int IS_CCPRIOR = (EVENT_TYPE == EVENT_TYPE_CCPRIOR);

  int MEMF    = LEN_MALLOC  * sizeof(float);
  int MEMI    = LEN_MALLOC  * sizeof(int);
  int MEMS    = LEN_MALLOC  * sizeof(short int);
  int MEMB    = LEN_MALLOC  * sizeof(bool);
  int MEMC    = LEN_MALLOC  * sizeof(char*);
  int MEMC2   = MXCHAR_CCID * sizeof(char);

  bool DOBIAS_MU = ( INPUTS.opt_biasCor & MASK_BIASCOR_MU     ) ;
  bool IDEAL     = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  bool CHECK_DUPLICATE = ( INPUTS.iflag_duplicate > 0 );

  int  NLCPAR_LOCAL = NLCPAR ; // beware that ILCPAR_MIN/MAX isn't set yet
  if ( DOBIAS_MU ) { NLCPAR_LOCAL++ ; }

  bool USE_FIELD = (INPUTS.use_fieldGroup_biasCor || INPUTS.NFIELD>0);
  long long MEMTOT=0;
  float f_MEMTOT, f_MEM ;
  int  i, isn, MEMF_TMP2, MEMF_TMP1, MEMF_TMP ;
  int debug_malloc = INPUTS.debug_malloc ;
  char fnam[] = "malloc_TABLEVAR" ;

  // ------------- BEGIN --------------

  print_debug_malloc(opt*debug_malloc,fnam);

  if ( opt > 0 ) {   
    
    TABLEVAR->name =  (char**)malloc(MEMC);
    for(i=0; i<LEN_MALLOC; i++ )  { 
      TABLEVAR->name[i] = (char*)malloc(MEMC2); 
      MEMTOT += MEMC2; 
      TABLEVAR->name[i][0] = 0 ;
    }

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

    f_MEM = 0.0;
    f_MEM += malloc_float3D(+1, LEN_MALLOC, NLCPAR, NLCPAR,
			    &TABLEVAR->covmat_fit );
    f_MEM += malloc_float3D(+1, LEN_MALLOC, NLCPAR, NLCPAR,
			    &TABLEVAR->covmat_tot );
    MEMTOT += (int)(f_MEM * 1.0E6); // convert MB back to bytes (Jun 2021)

    if ( IDEAL ) { 
      TABLEVAR->x0_ideal = (float *) malloc(MEMF); MEMTOT+=MEMF;
      for (i=0; i < NLCPAR_LOCAL ; i++ ) 
	{ TABLEVAR->fitpar_ideal[i] = (float *) malloc(MEMF); MEMTOT+=MEMF; }
    }

    if ( IS_DATA ) 
      { TABLEVAR->peakmjd = (float *) malloc(MEMF); MEMTOT += MEMF;  }

    TABLEVAR->warnCov       = (bool  *) malloc(MEMB); MEMTOT+=MEMB;
    TABLEVAR->x0            = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->x0err         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x0x1      = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x0c       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->COV_x1c       = (float *) malloc(MEMF); MEMTOT+=MEMF;

    TABLEVAR->zhd           = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhderr        = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zcmb          = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zcmberr       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhel          = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zhelerr       = (float *) malloc(MEMF); MEMTOT+=MEMF;

    if ( INPUTS.zspec_errmax_idsample > 0.0 ) {
      TABLEVAR->zprior      = (float *) malloc(MEMF); MEMTOT+=MEMF;
      TABLEVAR->zpriorerr   = (float *) malloc(MEMF); MEMTOT+=MEMF;
    }

    TABLEVAR->vpec          = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->vpecerr       = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->zmuerr        = (float *) malloc(MEMF); MEMTOT+=MEMF; // 6/2020
    TABLEVAR->snrmax        = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->pIa           = (float *) malloc(MEMF); MEMTOT+=MEMF; 

    MEMTOT += malloc_TABLEVAR_HOST(LEN_MALLOC,TABLEVAR,VARNAME_LOGMASS);

    TABLEVAR->IDSURVEY      = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->IDSAMPLE      = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->SNTYPE        = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->OPT_PHOTOZ    = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->IS_PHOTOZ     = (bool      *) malloc(MEMB); MEMTOT+=MEMB;
    TABLEVAR->IZBIN         = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->CUTMASK       = (int *)       malloc(MEMI); MEMTOT+=MEMI;

    TABLEVAR->ICUTWIN_GAMMA       = -9 ;
    TABLEVAR->ICUTWIN_VARNAME_PIA = -9 ;
    for(i=0; i < INPUTS.NCUTWIN; i++ ) 
      { MEMTOT += malloc_TABLEVAR_CUTVAL(LEN_MALLOC, i, TABLEVAR ); }
  
    TABLEVAR->SIM_NONIA_INDEX  = (short int *) malloc(MEMS); MEMTOT+=MEMS;
    TABLEVAR->SIM_ZCMB         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_VPEC         = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_MU           = (float *) malloc(MEMF); MEMTOT+=MEMF;
    TABLEVAR->SIM_MUz          = (float *) malloc(MEMF); MEMTOT+=MEMF;
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
    free(TABLEVAR->zhd);      free(TABLEVAR->zhderr);
    free(TABLEVAR->zcmb);     free(TABLEVAR->zcmberr);
    free(TABLEVAR->zhel);     free(TABLEVAR->zhelerr);
    free(TABLEVAR->vpec);     free(TABLEVAR->vpecerr);
    free(TABLEVAR->zmuerr);
    free(TABLEVAR->snrmax);

    if ( INPUTS.zspec_errmax_idsample > 0.0 ) 
      { free(TABLEVAR->zprior); free(TABLEVAR->zpriorerr); }

    free(TABLEVAR->COV_x0x1);
    free(TABLEVAR->COV_x0c);
    free(TABLEVAR->COV_x1c);
    free(TABLEVAR->warnCov);
    
    free(TABLEVAR->IDSURVEY);
    free(TABLEVAR->IDSAMPLE);
    free(TABLEVAR->SNTYPE);
    free(TABLEVAR->OPT_PHOTOZ);
    free(TABLEVAR->IS_PHOTOZ);
    free(TABLEVAR->IZBIN);
    free(TABLEVAR->CUTMASK);

    for(i=0; i < INPUTS.NCUTWIN; i++ ) { 
      if ( INPUTS.LCUTWIN_RDFLAG[i] ) { free(TABLEVAR->CUTVAL[i]); }
    }

    free(TABLEVAR->SIM_NONIA_INDEX);
    free(TABLEVAR->SIM_ZCMB);
    free(TABLEVAR->SIM_VPEC);
    free(TABLEVAR->SIM_MU);
    free(TABLEVAR->SIM_MUz);
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
int malloc_TABLEVAR_HOST(int LEN_MALLOC, TABLEVAR_DEF *TABLEVAR, 
			   char *VARNAME) {

  int MEMF    = LEN_MALLOC  * sizeof(float);
  int MEMTOT  = 0 ;
  if ( strcmp(VARNAME_LOGMASS,VARNAME)== 0 )  {
    TABLEVAR->host_logmass = (float*) malloc(MEMF);
    MEMTOT += MEMF;
  }
  else if ( strcmp(VARNAME_LOGSFR,VARNAME)== 0 ) {
    TABLEVAR->host_logsfr = (float*) malloc(MEMF);
    MEMTOT += MEMF;
  }
  else if ( strcmp(VARNAME_LOGsSFR,VARNAME)== 0 ) {
    TABLEVAR->host_logssfr = (float*) malloc(MEMF);
    MEMTOT += MEMF;
  }
  else if ( strcmp(VARNAME_COLOR,VARNAME)== 0 ) {
    TABLEVAR->host_color = (float*) malloc(MEMF);
    MEMTOT += MEMF;
  }

  return MEMTOT ;

} // end malloc_TABLEVAR_HOST

// *************************************************
int SNTABLE_READPREP_HOST(char *VARNAME, int ISTART, int LEN, 
			  TABLEVAR_DEF *TABLEVAR) {

  // Created Mar 7 2022
  // Call SNTABLE_READPREP_VARDEF for appropriate TABLEVAR pointer
  // to host variable based on input VARNAME.
  //
  // Inputs:
  //   *VARNAME   name of host variable in FITRES file
  //   ISTART     start index for TABLEVAR array
  //   LEN        length of TABLEVAR array
  //   TABLEVAR   structure with host_xxx arrays.

  int  ivar = -9, VBOSE = 1;
  char varCast[80];
  char fnam[] = "SNTABLE_READPREP_HOST";
  
  // ----------- BEGIN ------------

  sprintf(varCast,"%s:F", VARNAME);

  if ( strcmp(VARNAME,VARNAME_LOGMASS) == 0 ) {
      ivar = SNTABLE_READPREP_VARDEF(varCast, &TABLEVAR->host_logmass[ISTART], 
				     LEN, VBOSE);
  }
  else if ( strcmp(VARNAME,VARNAME_LOGSFR) == 0 ) {
      ivar = SNTABLE_READPREP_VARDEF(varCast, &TABLEVAR->host_logsfr[ISTART], 
				     LEN, VBOSE);
  }
  else if ( strcmp(VARNAME,VARNAME_LOGsSFR) == 0 ) {
      ivar = SNTABLE_READPREP_VARDEF(varCast, &TABLEVAR->host_logssfr[ISTART], 
				     LEN, VBOSE);
  }
  else if ( strcmp(VARNAME,VARNAME_COLOR) == 0 ) {
      ivar = SNTABLE_READPREP_VARDEF(varCast, &TABLEVAR->host_color[ISTART], 
				     LEN, VBOSE);
  }

  return ivar ;
} // end SNTABLE_READPREP_HOST


// ***************************************************
int malloc_TABLEVAR_CUTVAL(int LEN_MALLOC, int icut,
			     TABLEVAR_DEF *TABLEVAR ) {

  // Oct 29 2020: make substitution if CUTNAME = "varname_pIa"

  int MEMTOT = 0;
  int MEMF   = LEN_MALLOC * sizeof(float);
  int debug_malloc = INPUTS.debug_malloc ;
  char *CUTNAME = INPUTS.CUTWIN_NAME[icut] ;
  char *varname_pIa = INPUTS.varname_pIa ;
  char fnam[] = "malloc_TABLEVAR_CUTVAL" ;

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

  else if ( strcmp(CUTNAME,varname_pIa) == 0 )
    { TABLEVAR->CUTVAL[icut] = TABLEVAR->pIa; }

  else if ( strcmp(CUTNAME,"varname_pIa") == 0 )  {  // Oct 29 2020
    TABLEVAR->CUTVAL[icut] = TABLEVAR->pIa; 
    TABLEVAR->ICUTWIN_VARNAME_PIA = icut;  
    sprintf(CUTNAME, "%s", varname_pIa) ;
  }

  // IDSURVEY and SNTYPE are int ??

  else {
    print_debug_malloc(+1*debug_malloc,fnam);
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
  long long MEMTOT = 0 ;
  float f_MEMTOT;
  int debug_malloc = INPUTS.debug_malloc ;
  char fnam[] = "malloc_FITPARBIAS_ALPHABETA" ;
  // ----------- BEGIN ----------

  print_debug_malloc(opt*debug_malloc,fnam);

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
  //  Jul 06 2020: check SUBPROCESS GENPDF-variables
  //  Dec 11 2020: read zhel and zhelerr
  //  Jan 28 2021: read peakmjd for duplicate check (data only)
  //  Jun 16 2021: read explicit logmass rather than thru CUTWIN
  //  Mar 02 2022: read zPRIOR[ERR]

  int EVENT_TYPE = TABLEVAR->EVENT_TYPE;
  int IS_DATA    = ( EVENT_TYPE == EVENT_TYPE_DATA);
  int IS_BIASCOR = ( EVENT_TYPE == EVENT_TYPE_BIASCOR);
  //  int IS_CCPRIOR = ( EVENT_TYPE == EVENT_TYPE_CCPRIOR);

  int VBOSE = 1 ;  // verbose, but no abort on missing variable
  bool FIRSTFILE = ( ISTART == 0 ) ;
  bool USE_FIELD = ( INPUTS.use_fieldGroup_biasCor>0 || INPUTS.NFIELD>0);
  bool IDEAL          = ( INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  bool CHECK_DUPL     = ( INPUTS.iflag_duplicate > 0 ) ;
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
    TABLEVAR->IVAR_ZPRIOR       = -9 ;

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
    TABLEVAR->IS_PHOTOZ[irow]  =  false ;
    TABLEVAR->SNTYPE[irow]     = -9 ;

    TABLEVAR->vpec[irow]       =  0.0 ;
    TABLEVAR->vpecerr[irow]    =  0.0 ;
    TABLEVAR->zhd[irow]        = -9.0 ;
    TABLEVAR->zhderr[irow]     = -9.0 ;
    TABLEVAR->zmuerr[irow]     = -9.0 ;
    TABLEVAR->zcmb[irow]       = -9.0 ;
    TABLEVAR->zcmberr[irow]    = -9.0 ;
    TABLEVAR->zhel[irow]       = -9.0 ;
    TABLEVAR->zhelerr[irow]    = -9.0 ;
    TABLEVAR->host_logmass[irow] = -9.0 ;
    TABLEVAR->snrmax[irow]     =  0.0 ;
    TABLEVAR->warnCov[irow]    =  false ;

    if ( INPUTS.zspec_errmax_idsample > 0.0 ) 
      { TABLEVAR->zprior[irow]  = TABLEVAR->zpriorerr[irow]  = -9.0 ; }

    TABLEVAR->SIM_NONIA_INDEX[irow]  = -9 ;
    TABLEVAR->SIM_X0[irow]           = -9.0 ;
    TABLEVAR->SIM_FITPAR[0][irow]    = -9.0 ;
    TABLEVAR->SIM_FITPAR[1][irow]    = -9.0 ;
    TABLEVAR->SIM_FITPAR[2][irow]    = -9.0 ;
    TABLEVAR->SIM_ALPHA[irow]        = -9.0 ;
    TABLEVAR->SIM_BETA[irow]         = -9.0 ;
    TABLEVAR->SIM_GAMMADM[irow]      =  0.0 ; // allow missing gammadm
    TABLEVAR->SIM_MU[irow]           = -9.0 ;
    TABLEVAR->SIM_MUz[irow]          = -9.0 ;
    TABLEVAR->SIM_ZCMB[irow]         = -9.0 ;
    TABLEVAR->SIM_VPEC[irow]         =  0.0 ;
  }

  // - - - - - - - - prep strings - - - - - - - 
  sprintf(vartmp, "CID:C*%d  CCID:C*%d  GALID:C*%d", 
	  MXCHAR_CCID, MXCHAR_CCID, MXCHAR_CCID);
 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->name[ISTART], 
			  LEN, VBOSE) ;

  if ( USE_FIELD ) {
    sprintf(vartmp, "FIELD:C*%d", MXCHAR_CCID ); 
    ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->field[ISTART], 
				   LEN, VBOSE) ;
    if ( ivar < 0 ) {
      sprintf(c1err,"Required FIELD column missing");
      sprintf(c2err,"Check CUTFIELD or fieldGroup_biascor keys");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
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

  sprintf(vartmp,"zCMB:F"); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zcmb[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"zCMBERR:F"); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zcmberr[ISTART], 
			  LEN, VBOSE);

  sprintf(vartmp,"zHEL:F");   // Dec 11 2020
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zhel[ISTART], 
			  LEN, VBOSE);
  sprintf(vartmp,"zHELERR:F"); 
  SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zhelerr[ISTART], 
			  LEN, VBOSE);

  if ( INPUTS.zspec_errmax_idsample > 0.0 ) {
    sprintf(vartmp,"zPRIOR:F"); 
    ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zprior[ISTART], 
				   LEN, VBOSE);
    sprintf(vartmp,"zPRIORERR:F"); 
    ivar2 = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->zpriorerr[ISTART], 
				    LEN, VBOSE);
    if ( ivar > 0 && ivar2>0 ) { TABLEVAR->IVAR_ZPRIOR = ivar; }
  }

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

  char *varname_gamma = INPUTS.varname_gamma ;
  if ( strlen(varname_gamma) > 0 ) {
    sprintf(vartmp,"%s:F", varname_gamma);
    ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->host_logmass[ISTART],
                                   LEN, VBOSE );
  }

  if ( IS_DATA ) { // Jan 28 2021
      sprintf(vartmp,"PKMJD:F" ) ;
      ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->peakmjd[ISTART], 
				     LEN, VBOSE);
  }

  // - - - - - - -
  //read CUTWIN variables 
  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
    cutname = INPUTS.CUTWIN_NAME[icut]; 
    RDFLAG  = INPUTS.LCUTWIN_RDFLAG[icut] ;
    sprintf(vartmp, "%s:F", cutname );
    if ( !usesim_CUTWIN(vartmp)  ) { continue ; }

    if ( RDFLAG ) {
      ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->CUTVAL[icut][ISTART], 
				     LEN, VBOSE );
    }
    else {
      ivar = IVAR_READTABLE_POINTER(cutname) ;   // May 8 2020 
    }
    TABLEVAR->DOFLAG_CUTWIN[icut] = 
      set_DOFLAG_CUTWIN(ivar, icut, IS_DATA) ;

  } // end icut


  // - - - - - - - - SIM_XXX - - - - - - - - - -
  sprintf(vartmp,"SIM_NONIA_INDEX:S SIM_TEMPLATE_INDEX:S" );
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_NONIA_INDEX[ISTART], 
				 LEN, VBOSE );
  if ( IS_DATA && INPUTS.force_realdata ) { ivar = -9; } // Jun 17 2021

  FOUNDKEY_SIM=0;

  if ( ivar >=0 ) 
    { TABLEVAR->IS_SIM = true ;   FOUNDKEY_SIM=1; }
  else
    { TABLEVAR->IS_DATA = true;   goto CHECK_SUBPROCESS ; }


  // --------------------------------------
  // here and below is for simulated data 
  // ---------------------------------------

  // note that IS_DATA refers to datafile= argument, and can be
  // real data or simulated data. ISDATA_REAL is false for sim data.
  if ( IS_DATA ) { 
    ISDATA_REAL = 0 ;   // not real data -> sim data
    // if the 64 blind-sim bit isn't set by xuser, set blindFlag=0
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

  sprintf(vartmp,"SIM_alpha:F SIMalpha:F" ) ;
  ivar = SNTABLE_READPREP_VARDEF(vartmp, &TABLEVAR->SIM_ALPHA[ISTART],
				 LEN, VBOSE );
  sprintf(vartmp,"SIM_beta:F  SIMbeta:F" ) ;
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


 CHECK_SUBPROCESS:
#ifdef USE_SUBPROCESS  
  if ( SUBPROCESS.USE ) { 
    SUBPROCESS_READPREP_TABLEVAR(IFILE, ISTART, LEN, TABLEVAR); 
  }
#endif

  return ;

} // end  SNTABLE_READPREP_TABLEVAR


// ***********************************************
void compute_more_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR ) {

  // Created Jun 2019
  // After reading table and loading TABLEVAR arrays, 
  // this function defines 'more' TABLEVAR quantities;
  //  * mB, mBerr
  //  * mumodel
  //  * cov matrices
  //  * fix cov matries
  //  * set_CUTMASK
  //
  // Inputs:
  //   ISN        : SN index for TABLEVAR arrays
  //   TABLEVAR   : structure of arrays
  //
  //  Beware that TABLEVAR are float to save memory,
  //  particularly for large BIASCOR samples. But local variables
  //  and calculations are done with double.
  //
  // Aug 22 2019: set logmass to table value, or to p7
  //
  // Feb 24 2020: load TABLEVAR->fitpar[INDEX_mu][ISN]
  // Jun 27 2020: local float -> double
  // Dec 11 2020: if zhel<0, set zhel = zhd. Pass zhel to cosmodl().

  int EVENT_TYPE   = TABLEVAR->EVENT_TYPE;
  bool IS_DATA     = (EVENT_TYPE == EVENT_TYPE_DATA);
  bool IS_BIASCOR  = (EVENT_TYPE == EVENT_TYPE_BIASCOR);
  bool IS_SIM      = TABLEVAR->IS_SIM ; // data are sim

  //  int IS_CCPRIOR  = (EVENT_TYPE == EVENT_TYPE_CCPRIOR);
  bool IDEAL       = (INPUTS.opt_biasCor & MASK_BIASCOR_COVINT ) ;
  bool FIRST_EVENT = (ISN == 0) ;
  //  bool LAST_EVENT  = (ISN == TABLEVAR->NSN_ALL-1 );
  char *STRTYPE   = STRING_EVENT_TYPE[EVENT_TYPE];  

  int IVAR_GAMMA   = INFO_DATA.TABLEVAR.ICUTWIN_GAMMA ;

  bool  DO_BIASCOR_SAMPLE = (INPUTS.opt_biasCor & MASK_BIASCOR_SAMPLE);
  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  bool  USE_FIELDGROUP    = (INPUTS.use_fieldGroup_biasCor > 0) ;

  double vpecerr_orig   = (double)TABLEVAR->vpecerr[ISN];
  if ( INPUTS.zpecerr > 1.0E-9 ) 
    { TABLEVAR->vpecerr[ISN] = INPUTS.zpecerr * LIGHT_km; } // 6.30.2020

  int  IDSURVEY    = (int)TABLEVAR->IDSURVEY[ISN];
  int  SNTYPE      = (int)TABLEVAR->SNTYPE[ISN] ; 
  int  OPT_PHOTOZ  = (int)TABLEVAR->OPT_PHOTOZ[ISN];
  char *name       =      TABLEVAR->name[ISN];
  double x0        = (double)TABLEVAR->x0[ISN];
  double x0err     = (double)TABLEVAR->x0err[ISN];
  double x1err     = (double)TABLEVAR->fitpar_err[INDEX_x1][ISN];
  double cerr      = (double)TABLEVAR->fitpar_err[INDEX_c][ISN];
  double zhd       = (double)TABLEVAR->zhd[ISN];
  double zhderr    = (double)TABLEVAR->zhderr[ISN];
  double vpec      = (double)TABLEVAR->vpec[ISN];
  double zhel      = (double)TABLEVAR->zhel[ISN];   
  double vpecerr   = (double)TABLEVAR->vpecerr[ISN];
  double COV_x0x1  = (double)TABLEVAR->COV_x0x1[ISN];
  double COV_x0c   = (double)TABLEVAR->COV_x0c[ISN];
  double COV_x1c   = (double)TABLEVAR->COV_x1c[ISN];

  double SIM_X0, SIM_ZCMB, SIM_MU, SIM_MUz;
  int    SIM_NONIA_INDEX;

  double mB, mBerr, mB_orig, x1, c, sf, x0_ideal, mB_ideal ;
  double zpec, zpecerr, zcmb, zMIN, zMAX, zhderr_tmp, logmass;
  double dl_hd, dl_sim, dl, dl_ratio, dmu ;
  double covmat8_fit[NLCPAR][NLCPAR], covmat8_int[NLCPAR][NLCPAR];
  double covmat8_tot[NLCPAR][NLCPAR], covtmp8 ;
  int   OPTMASK_COV, ISTAT_COV, i, i2, IDSAMPLE;
  bool  IS_SPECZ, IS_PHOTOZ;
  char  *field, SURVEYGROUP[100], FIELDGROUP[100], STRINGOPT[40];
  char NONE[] = "NONE";
  char fnam[] =  "compute_more_TABLEVAR";

  // ----------------- BEGIN ------------------

  if ( INPUTS.cat_only ) { return; }

  if ( FIRST_EVENT ) { TABLEVAR->NCOVFIX  =  0 ;  }

  // convert x0 and error to mB[err]
  x0 = TABLEVAR->x0[ISN];  x0err=TABLEVAR->x0err[ISN]; 
  mB = mBerr = -9.0; sf = 0.0 ;
  if ( x0 > 0.0 ) { 
    mB = -2.5*log10f(x0); 
    sf = -2.5/(x0*LOGTEN) ;
    mBerr = fabs(x0err * sf);
  }
  TABLEVAR->fitpar[INDEX_mB][ISN]     = (float)mB ;  // this mB has no offset
  TABLEVAR->fitpar_err[INDEX_mB][ISN] = (float)mBerr ;

  if ( IDSURVEY < 0 || IDSURVEY > MXIDSURVEY ) {
    sprintf(c1err,"Invalid IDSURVEY=%d (unsigned short) for %s SNID=%s", 
	    IDSURVEY, STRTYPE, name );
    sprintf(c2err,"Check $SNDATA_ROOT/SURVEY.DEF" );
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }


  if ( !IS_DATA) {  // biascor or CC prior
    SIM_X0    = (double)TABLEVAR->SIM_X0[ISN];
    SIM_ZCMB  = (double)TABLEVAR->SIM_ZCMB[ISN];
    SIM_MU    = (double)TABLEVAR->SIM_MU[ISN] ; 
    SIM_MUz   = -9.0 ;   // SIM_MU at zHD, computed below
    SIM_NONIA_INDEX = TABLEVAR->SIM_NONIA_INDEX[ISN];
  }

  // - - - - - - - - - - - - - - - - - -   
  // 6.29.2020: load zmuerr used for muerr += dmu/dz x zerr
  if ( INPUTS.restore_sigz ) 
    { TABLEVAR->zmuerr[ISN] = zhderr; } // legacy: full zerr
  else
    { TABLEVAR->zmuerr[ISN] = vpecerr/LIGHT_km; } // only vpec component

  // increment IDSURVEY-dependent stuff
  TABLEVAR->NSN_PER_SURVEY[IDSURVEY]++ ;
  zMIN = TABLEVAR->zMIN_PER_SURVEY[IDSURVEY] ;
  zMAX = TABLEVAR->zMAX_PER_SURVEY[IDSURVEY] ;
  if ( zhd > 0.0 ) {
    if ( zhd < zMIN ) { TABLEVAR->zMIN_PER_SURVEY[IDSURVEY] = (float)zhd; }
    if ( zhd > zMAX ) { TABLEVAR->zMAX_PER_SURVEY[IDSURVEY] = (float)zhd; }
  }

  // - - - - covariance matrix - - - - - - - -  
  covmat8_fit[INDEX_mB][INDEX_mB] = (x0err*x0err*sf*sf) ;
  covmat8_fit[INDEX_mB][INDEX_x1] = (COV_x0x1 * sf);
  covmat8_fit[INDEX_mB][INDEX_c]  = (COV_x0c  * sf) ;
  covmat8_fit[INDEX_x1][INDEX_x1] = (x1err * x1err) ;
  covmat8_fit[INDEX_x1][INDEX_c]  = (COV_x1c) ;
  covmat8_fit[INDEX_c][INDEX_c]   = (cerr * cerr) ;

  // symmetric part of off-diagonals
  covmat8_fit[INDEX_x1][INDEX_mB] = (COV_x0x1 * sf) ;
  covmat8_fit[INDEX_c][INDEX_mB]  = (COV_x0c  * sf) ;
  covmat8_fit[INDEX_c][INDEX_x1]  = (COV_x1c);

  OPTMASK_COV=0;
  if ( strcmp(name,INPUTS.SNID_MUCOVDUMP)==0 ) { OPTMASK_COV=4; } //dump COV
  update_covMatrix(name, OPTMASK_COV, NLCPAR, covmat8_fit, 
		   EIGMIN_COV, &ISTAT_COV ); 
  if ( ISTAT_COV != 0 ) { TABLEVAR->warnCov[ISN]=true;  TABLEVAR->NCOVFIX++; }

  // transfer local covmat8 to TABLEVAR float-storage

  for(i=0; i < NLCPAR; i++ ) {
    for(i2=0; i2 < NLCPAR; i2++ ) {
      TABLEVAR->covmat_fit[ISN][i][i2] = (float)covmat8_fit[i][i2];
    }
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

  // - - - - - - - - - - - - - 
  // make sure that zHD is positive. Skip test for user-input zVARNAME.
  if ( zhd < 0 && strlen(INPUTS.varname_z)==0 ) { 
    sprintf(c1err,"Invalid %s redshift = %f ", STRTYPE, zhd );
    sprintf(c2err,"for ISN=%d, SNID='%s' ", ISN, name );
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
  }

  // - - - - - - - 
  // Stuff for data only,
  if ( IS_DATA && zhd > 0 ) {

    // if user-input zpecerr > 0, subtract out original zpecerr and 
    // add user input zpecerr:
    // zhd = (1+zcmb)*(1+zpec) -1 = zcmb + zpec + zcmb*zpec
    // zcmb(1+zpec) = zhd - zpec     
    zpec       =  vpec/LIGHT_km ;
    zcmb       =  (zhd - zpec)/(1.0+zpec) ;
    zhderr_tmp = (float)zerr_adjust(zcmb, zhderr, vpecerr_orig, name);
    TABLEVAR->zhderr[ISN] = zhderr_tmp ;
    
    // xxxxxxxxxxxxx
    if ( ISN < -10 ) {
      double dmuz  = fcn_muerrz(1, zhd, zhderr );
      printf(" xxx %s: ISN=%2d z = %.3f +_ (%.3f->%.3f)  dmuz=%.3f \n",
	     fnam, ISN, zhd, zhderr, TABLEVAR->zhderr[ISN], dmuz);
    }
    // xxxxxxxxxxxxx

    if ( zhel < 0.0 ) {  // Dec 11 2020
      zhel = zhd;     // approx fix since we don't have RA,DEC to transform 
      TABLEVAR->zhel[ISN] = (float)zhel ;
    }

    // if cosmo params are fixed, then store mumodel for each event.
    if ( INPUTS.FLOAT_COSPAR == 0 ) {
      dl = cosmodl_forFit(zhd,zhd,INPUTS.COSPAR);  // legacy
      // dl = cosmodl_forFit(zhel,zhd,INPUTS.COSPAR); // fixed
      TABLEVAR->mumodel[ISN] = (float)(5.0*log10(dl) + 25.0);
    }

    // check z-cheat option 
    if ( INPUTS.uzsim && IS_SIM == true ) { 
      TABLEVAR->zhd[ISN]     = TABLEVAR->SIM_ZCMB[ISN];
      TABLEVAR->zhderr[ISN]  = 1.0E-7;
    }


    if ( INPUTS.force_pIa >= 0.0 ) 
      { TABLEVAR->pIa[ISN] = INPUTS.force_pIa; } 

    if ( INPUTS.perfect_pIa )  {
      if ( !IS_SIM ) {
	sprintf(c1err,"Cannot force_pIa=perfect for real data.");
	sprintf(c2err,"This option works only for sim data.") ;
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

  IS_SPECZ  = IS_SPECZ_TABLEVAR(ISN,TABLEVAR); 
  IS_PHOTOZ = !IS_SPECZ ;
  TABLEVAR->IS_PHOTOZ[ISN] = IS_PHOTOZ;

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
      IDSAMPLE = get_IDSAMPLE(IDSURVEY,IS_PHOTOZ,FIELDGROUP,SURVEYGROUP,0) ;
    }
    else {
      IDSAMPLE = 0 ;   FIELDGROUP[0] = 0 ;
    }
    TABLEVAR->IDSAMPLE[ISN] = IDSAMPLE;

    // misc stuff for BIASCOR & CCPRIOR
    // increment number of events per IDSAMPLE (for biasCor)
    if ( IDSAMPLE >= 0 ) 
      { SAMPLE_BIASCOR[IDSAMPLE].NSN[EVENT_TYPE]++ ; }

    // load sim_mb = -2.5log10(sim_x0)
    mB = -9.0;
    if ( SIM_X0 > 0.0 ) { mB = -2.5*log10(SIM_X0); }
    mB_orig = (double)TABLEVAR->SIM_FITPAR[INDEX_mB][ISN];
    // xxx    mB_off  = mB_orig - mB;
    TABLEVAR->SIM_FITPAR[INDEX_mB][ISN] = (float)mB;  // no offset

    // legacy check on user zpecerr 
    if ( INPUTS.zpecerr > 1.0E-9 ) {
      zpecerr = INPUTS.zpecerr * ( 1.0 + zhd ) ;
      TABLEVAR->zhderr[ISN]= (float)sqrt(zhderr*zhderr + zpecerr*zpecerr);
    }

  }  // end not-DATA

  if ( IS_DATA && ISN < -10 ) {
    double zmuerr = TABLEVAR->zmuerr[ISN] ;
    double dmuz  = fcn_muerrz(1, zhd, zmuerr );
    printf(" xxx %s: DATA ISN=%d(%s)  z=%5.3f  zmuerr=%.4f  (vpecerr=%.1f)\n",
           fnam, ISN, name, zhd, zmuerr, vpecerr );
    fflush(stdout);
  }

  // Compute SIM_MUz from zHD (SIM_MU is from true SIM_zCMB)
  TABLEVAR->SIM_MUz[ISN] = -9.0 ;
  if ( !IS_DATA ) {
    // beware that INPUTS.COSPAR are not always the same as SIM-COSPAR,
    // so be careful computing SIM_MU at zHD
    dl_hd    = cosmodl(zhd,zhd,INPUTS.COSPAR);
    dl_sim   = cosmodl(SIM_ZCMB,SIM_ZCMB,INPUTS.COSPAR);
    dl_ratio = dl_hd/dl_sim ;
    dmu = 5.0*log10(dl_ratio);
    SIM_MUz = SIM_MU + dmu ;
    TABLEVAR->SIM_MUz[ISN] = (float)SIM_MUz ;
  }
 
  // - - - - -
  if ( !IS_DATA  && DO_BIASCOR_MU ) { 
    int NBINb = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;  
    if ( SIM_NONIA_INDEX == 0 && NBINb == 1 )
      { INFO_BIASCOR.DUST_FLAG=true; }

    bool restore_mucovadd_bug =(INPUTS.restore_mucovadd_bug&1)>0;
    if (restore_mucovadd_bug){ INFO_BIASCOR.DUST_FLAG=false; }

    // Mainly for BS21 model:
    // Prepare option to bias-correct MU instead of correcting mB,x1,c 
    // Note that true Alpha,Beta,GammaDM are used for mu_obs.
    // Beware that M0_DEFAULT may be fragile.
    double Alpha, Beta, GammaDM, mu_obs, mu_true; 
    mB      = (double)TABLEVAR->fitpar[INDEX_mB][ISN] ;
    x1      = (double)TABLEVAR->fitpar[INDEX_x1][ISN] ;
    c       = (double)TABLEVAR->fitpar[INDEX_c][ISN] ;

    // 4/10/2020: temp hack;later should measure A,B from slopes of HR vs x1,c
    // Maybe need option to use true A,B from sims ?
    Alpha   = INPUTS.parval[IPAR_ALPHA0]; 
    Beta    = INPUTS.parval[IPAR_BETA0]; 
    GammaDM = TABLEVAR->SIM_GAMMADM[ISN] ;

    // Mar 3 2022: hack on hack; set sim_beta to user p2 ~ fitted beta,
    //   and override the intrinsic beta. Needed in muerrsq_biasCor.
    //   Still need to define a true Tripp-Beta for BS21 model
    //   This fix helps fix the anomalously low chi2.
    if ( INFO_BIASCOR.DUST_FLAG ) { 
      // intended for for BS21 model
      Alpha   = INPUTS.parval[IPAR_ALPHA0]; 
      Beta    = INPUTS.parval[IPAR_BETA0]; 
      GammaDM = TABLEVAR->SIM_GAMMADM[ISN] ;
      INFO_BIASCOR.TABLEVAR.SIM_BETA[ISN]  = Beta ; // overwrite !!!
    }
    else {
      // non-BS21 models with single color law
      Alpha   = TABLEVAR->SIM_ALPHA[ISN] ;
      Beta    = TABLEVAR->SIM_BETA[ISN] ;
      GammaDM = TABLEVAR->SIM_GAMMADM[ISN] ;
    }

    /*
    if ( ISN < -20 ) {
      printf(" xxx %s: ISN=%2d zhd=%.3f SIM_ZCMB=%.3f dmu=%.3f \n",
	     fnam, ISN, zhd, SIM_ZCMB, dmu );
    }
    */

    mu_obs  = mB - M0_DEFAULT + Alpha*x1 - Beta*c - GammaDM ;
    TABLEVAR->fitpar[INDEX_mu][ISN]      = (float)mu_obs ; 
    TABLEVAR->SIM_FITPAR[INDEX_mu][ISN]  = (float)SIM_MUz ;
  }

  if ( IS_BIASCOR && IDEAL ) {
    x0_ideal = (double)TABLEVAR->x0_ideal[ISN];
    mB_ideal = -2.5*log10f(x0_ideal); 
    TABLEVAR->fitpar_ideal[INDEX_mB][ISN] = (float)mB_ideal;
  }

  return ;

} // end compute_more_TABLEVAR


// ***********************************************
bool IS_SPECZ_TABLEVAR(int ISN, TABLEVAR_DEF *TABLEVAR) {

  // Created Mar 2 2022 by R.Kessler
  // Function returns True if this event (ISN index) has a
  // spectroscopic redshift base on either
  //  + is not a photo-z fit (OPT_PHOTOZ=0) or
  //  + has accurate zpriorerr < INPUTS.zspec_errmax_idsample
  //

  int  OPT_PHOTOZ  = (int)TABLEVAR->OPT_PHOTOZ[ISN];
  bool IS_SPECZ = false;
  double zerr;
  char fnam[] = "IS_SPECZ_TABLEVAR" ;

  // ------------ BEGIN ------------

  if ( TABLEVAR->IVAR_ZPRIOR > 0 ) {
    // use original zhelio uncertainty for redshift prior in LC fit
    zerr = (double)TABLEVAR->zpriorerr[ISN];   
  }
  else {
    // try using zhelio uncertainty, but it might be inflated
    // from photo-z fit
    zerr = (double)TABLEVAR->zhelerr[ISN];   
  }

  IS_SPECZ  = ( OPT_PHOTOZ == 0 || zerr < INPUTS.zspec_errmax_idsample);

  return IS_SPECZ;

} // end IS_SPECZ_TABLEVAR

// ***********************************************
void compute_CUTMASK(int ISN, TABLEVAR_DEF *TABLEVAR ) {

  // Created May 12 2021
  // Call set_CUTMASK and set some counters.
  // Was part of compute_more_TABLEVAR, but seperated here
  // so that read_data_override can take effect before cuts.

  int EVENT_TYPE   = TABLEVAR->EVENT_TYPE;
  if ( INPUTS.cat_only ) { return; }

  set_CUTMASK(ISN, TABLEVAR);

  return;

} // end compute_CUTMASK

// ==================================================
void compute_more_INFO_DATA(void) {

  // Created Jun 2019
  // set initial muerr so that optional log(1/sigma) term in fcn-chi2
  // can be fixed on first iteration
  // 
  // Jan 2021: init set_fitwgt0[isn] = false

  int NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL ;
  int nfile_biasCor = INPUTS.nfile_biasCor ;
  double *ptr_sigCC = &INPUTS.parval[IPAR_H11+3];
  double muerrsq, sigCC, zhd, zmuerr, cov, covmat_tot[NLCPAR][NLCPAR] ;
  int    isn, CUTMASK, i, i2 ;
  char *name;
  char fnam[] = "compute_more_INFO_DATA";

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
      zmuerr = (double)INFO_DATA.TABLEVAR.zmuerr[isn];
      for(i=0; i < NLCPAR; i++ ) {
	for(i2=0; i2 < NLCPAR; i2++ ) {
	  cov = (double)INFO_DATA.TABLEVAR.covmat_tot[isn][i][i2];
	  covmat_tot[i][i2] = cov;
	}
      }

      // strip INPUTS.parval here to avoid crash with cat_only option
      double alpha      = INPUTS.parval[IPAR_ALPHA0] ;
      double beta       = INPUTS.parval[IPAR_BETA0] ;
      double gamma      = INPUTS.parval[IPAR_GAMMA0] ;
      muerrsq = fcn_muerrsq(name,alpha,beta,gamma, covmat_tot,zhd,zmuerr, 0);
      sigCC   = ptr_sigCC[0] + zhd*ptr_sigCC[1] + zhd*zhd*ptr_sigCC[2];
    }

    INFO_DATA.muerrsq_last[isn]  = muerrsq ;
    INFO_DATA.muerr_last[isn]    = sqrt(muerrsq);
    INFO_DATA.sigCC_last[isn]    = sigCC ;
    INFO_DATA.sqsigCC_last[isn]  = sigCC * sigCC ;

    if ( nfile_biasCor > 0 ) {
      INFO_DATA.set_fitwgt0[isn]   = false; // Jan 2021
    }
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
  //
  // Oct 20 2020:
  //   Clarify error message for undefined pIa; give explicit idsurvey and names.

  char *INPFILE     = TABLEVAR->INPUT_FILE[ifile];
  int   FIRST_EVENT = TABLEVAR->EVENT_RANGE[ifile][0];
  int   LAST_EVENT  = TABLEVAR->EVENT_RANGE[ifile][1];
  int   IVAR_pIa    = TABLEVAR->IVAR_pIa[ifile];
  char *varname_pIa = INPUTS.varname_pIa ;
  bool USE_PROBCC_ZERO = (INPUTS_PROBCC_ZERO.nidsurvey > 0 );
  int debug_malloc = INPUTS.debug_malloc ;

  char *survey, survey_missing_list[200];
  int evt, sntype, idsurvey, NFORCE, NOTFORCE ;
  int NIDSURVEY_NOTFORCE[MXIDSURVEY];
  char fnam[] = "store_input_varnames" ;

  // ------------ BEGIN --------------

  /*
  printf(" xxx %s hello, IVAR_pIa=%d  N_PCC_ZERO=%d\n",
	 fnam, IVAR_pIa, INPUTS_PROBCC_ZERO.nidsurvey );
  */

  survey_missing_list[0] = 0 ;
  for(idsurvey=0; idsurvey < MXIDSURVEY; idsurvey++ ) 
    { NIDSURVEY_NOTFORCE[idsurvey] = 0 ; }

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
	  { NOTFORCE++ ; NIDSURVEY_NOTFORCE[idsurvey]++; }
      }
      fprintf(FP_STDOUT,"\t Force pIa=1.0 for %d events from %s\n",
	     NFORCE, INPFILE); fflush(FP_STDOUT);
    }

    // abort if varname_pIa column is missing and there are events
    // which don't have forced pIa=1.
    if (  IVAR_pIa < 0 && NOTFORCE > 0 ) {
      print_preAbort_banner(fnam);
      printf("   input file: %s \n", INPFILE);
      printf("   has %d events with forced pIa=1\n", NFORCE);
      printf("   and %d undefined-pIa events.\n", NOTFORCE);
      printf("   Breakdown : \n");

      for(idsurvey=0; idsurvey < MXIDSURVEY; idsurvey++ ) {
	int NID = NIDSURVEY_NOTFORCE[idsurvey];
	if ( NID > 0 ) {
	  survey = SURVEY_INFO.SURVEYDEF_LIST[idsurvey] ;
	  printf("\t %d undefined-pIa have IDSURVEY=%d (%s) \n",
		 NID, idsurvey, survey );
	  catVarList_with_comma(survey_missing_list,survey);
	}
      }

      sprintf(c1err,"Missing column varname_pIa='%s' .",  varname_pIa);
      sprintf(c2err,"Try adding %s to check type_list_probcc0",
	      survey_missing_list) ;
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
    }

    
  } // end varname_pIa


  // - - - - - - - - - - -
  // store varnames for each input file; used later to
  // prepare/align output table columns.
  int NVAR =  READTABLE_POINTERS.NVAR_TOT ;
  int MEMC1 =  NVAR * sizeof(char*) ;
  int MEMC0 =  MXCHAR_VARNAME * sizeof(char);
  int ivar ;

  print_debug_malloc(+1*debug_malloc,fnam);
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
  //
  // Dec 4 2020: allow recylcing output fitres to input by chopping
  //             out variables from CUTMASK to end of list.

  int debug_malloc = INPUTS.debug_malloc ;
  int  NFILE = INPUTS.nfile_data;
  //  int  MEMC  = MXCHAR_VARNAME * sizeof(char);
  int  ifile, ifile2, NVAR, NVAR2, MXVAR, NVAR_TOT; 
  int  ivar_match, ivar, ivar2, NMATCH, i, NVAR_KEEP ;
  char *varName, *varName2 ;
  bool wildcard, MATCH, CHOP_VAR ;
  int LDMP = 0 ;

  char FIRST_VARNAME_APPEND[] = "CUTMASK" ;
  char fnam[] = "store_output_varnames" ;
  // ------------- BEGIN ------------
  
  // xxx  OUTPUT_VARNAMES.PERFECT_COLUMN_MATCH = true;

  // alloate VARNAMES memory for x2 number of variables in first file ...
  // should be enough
  MXVAR = 2*INFO_DATA.TABLEVAR.NVAR[0];
  print_debug_malloc(+1*debug_malloc,fnam);
  for(ivar=0; ivar < MXVAR; ivar++ ) {
    OUTPUT_VARNAMES.LIST[ivar] = (char*) malloc(MXVAR*sizeof(char) ) ;
    for(ifile=0; ifile < NFILE; ifile++ ) 
      { OUTPUT_VARNAMES.IVARMAP[ifile][ivar] = -888 ; }
  }

 
  // check trivial case with just one input file
  if ( NFILE == 1 ) {
    ifile = 0 ;    NVAR_KEEP=0;    CHOP_VAR = false;  
    NVAR = INFO_DATA.TABLEVAR.NVAR[ifile];
    for(ivar=0; ivar < NVAR; ivar++ ) { 
      varName = INFO_DATA.TABLEVAR.VARNAMES_LIST[ifile][ivar];
      if (strcmp(varName,FIRST_VARNAME_APPEND) == 0 ) { CHOP_VAR=true;}
      if ( CHOP_VAR ) { break; }
      sprintf(OUTPUT_VARNAMES.LIST[ivar],"%s", varName);
      OUTPUT_VARNAMES.IVARMAP[ifile][ivar] = ivar ; 
      NVAR_KEEP++ ;
    }

    OUTPUT_VARNAMES.NVAR_TOT = NVAR_KEEP;
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
    if (strcmp(varName,FIRST_VARNAME_APPEND) == 0 ) { continue; }

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
	if (strcmp(varName2,FIRST_VARNAME_APPEND) == 0 ) { continue; }

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
      if ( ivar_match < 0  ) {
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
  int IVAR_OPT_PHOTOZ = INFO_DATA.TABLEVAR.IVAR_OPT_PHOTOZ ;

  int isn, IDSURVEY, OPT_PHOTOZ, N, IDSAMPLE, i, NIDSURVEY[MXIDSURVEY] ;
  int  DUMPFLAG=0, NDMP = 0, NSN_DATA, CUTMASK  ; 
  bool IS_SPECZ, IS_PHOTOZ ;
  double zhd ;
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

  SURVEYGROUP[0] = 0; // initalize blank string

  if ( INPUTS.opt_biasCor == 0 ) { return ; }

  sprintf(BANNER,"Begin %s", fnam);
  fprint_banner(FP_STDOUT,BANNER);

      
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
    // xxx    SAMPLE_BIASCOR[i].OPT_PHOTOZ          = 0 ; // zSPEC is default
    SAMPLE_BIASCOR[i].IS_PHOTOZ           = false ; // zSPEC is default
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
    SAMPLE_BIASCOR[IDSAMPLE].IS_PHOTOZ   = false ; // zSPEC is default

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
    zhd        = INFO_DATA.TABLEVAR.zhd[isn];

    if(USE_FIELDGROUP) { FIELDDEF = INFO_DATA.TABLEVAR.field[isn]; }

    if ( CUTMASK ) { continue ; }

    if ( IDSURVEY < 0 || IDSURVEY > MXIDSURVEY ) {
      sprintf(c1err,"Invalid IDSURVEY=%d for SNID=%s", IDSURVEY, NAME_SN);
      sprintf(c2err,"Check $SNDATA_ROOT/SURVEY.DEF" );
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

    IS_SPECZ  = IS_SPECZ_TABLEVAR(isn,&INFO_DATA.TABLEVAR);
    IS_PHOTOZ = !IS_SPECZ ;

    if ( IS_SPECZ ) 
      { sprintf(zGROUP,"-zSPEC"); }
    else
      { sprintf(zGROUP,"-zPHOT"); }

    if ( IVAR_OPT_PHOTOZ < 0 ) { zGROUP[0] = 0 ; }

    DUMPFLAG = 0 ;
    IDSAMPLE = get_IDSAMPLE(IDSURVEY, IS_PHOTOZ, 
			    FIELDGROUP, SURVEYGROUP, DUMPFLAG);

    if ( IDSAMPLE < 0 ) {
      // store new SURVEY/FIELDGROUP entry
      N = NSAMPLE_BIASCOR ;  
      sprintf(SAMPLE_BIASCOR[N].NAME_FIELDGROUP,  "%s", FIELDGROUP  );
      sprintf(SAMPLE_BIASCOR[N].NAME_SURVEYGROUP, "%s", SURVEYGROUP ); 
      sprintf(SAMPLE_BIASCOR[N].STRINGOPT,        "%s", STRINGOPT   ); 
      SAMPLE_BIASCOR[N].IS_PHOTOZ = IS_PHOTOZ ;
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
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
      }
    }
   
    // execute mapindex function again to get final IDSAMPLE
    DUMPFLAG = NIDSURVEY[IDSURVEY] < 0 ; // for debug only
    IDSAMPLE = get_IDSAMPLE(IDSURVEY, IS_PHOTOZ,
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }

    INFO_DATA.TABLEVAR.IDSAMPLE[isn] = IDSAMPLE ; 

    SAMPLE_BIASCOR[IDSAMPLE].NSN[EVENT_TYPE_DATA]++ ;

    // keep track of min/max redshift for each sample

    if( zhd < SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMIN_DATA = zhd; }
    if( zhd > SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA ) 
      { SAMPLE_BIASCOR[IDSAMPLE].zMAX_DATA = zhd; }

    // check user option to skip biasCor for this sample
    if ( strstr(INPUTS.surveyList_noBiasCor,SURVEYDEF) != NULL ) 
      { SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR = 0 ; }
    
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
  int debug_malloc = INPUTS.debug_malloc ;
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
    
  print_debug_malloc(+1*debug_malloc,fnam);
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  print_debug_malloc(-1*debug_malloc,fnam);
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

  for(i=0; i < MXNUM_SAMPLE; i++ ) 
    { ptrFIELD[i] = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[i] ; }
  
  splitString(INPUTS.fieldGroup_biasCor, COMMA, MXNUM_SAMPLE, // inputs
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
  int debug_malloc = INPUTS.debug_malloc ;
  int  i, i2, NGRP, ID, NEVT ;
  char *ptrSURVEY[MXNUM_SAMPLE], *S ;
  char fnam[] = "set_SURVEYGROUP_biasCor" ;

  // ------------- BEGIN -------------

  // extract SURVEYGROUP sub strings (aug 18 2016)

  if ( USE_SURVEYGROUP ) {
    for(i=0; i < MXNUM_SAMPLE; i++ ) 
      { ptrSURVEY[i] = INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i] ; }

    splitString(INPUTS.surveyGroup_biasCor, COMMA, MXNUM_SAMPLE, // inputs
		&NGRP, ptrSURVEY ); // outputs
    INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR = NGRP ;
  }

  // convert each SURGROUP_LIST (string) into a list of IDSURVEYs

  // define temp string space for plus-separated strings.
  print_debug_malloc(+1*debug_malloc,fnam);
  char *ptrTmp[MXNUM_SAMPLE] ;
  for(i=0; i < MXNUM_SAMPLE; i++ ) 
    { ptrTmp[i] = (char*) malloc ( MXCHAR_CCID * sizeof(char) ) ;  }


  NGRP = INPUTS_SAMPLE_BIASCOR.NSURVEYGROUP_USR ;

  for(i=0; i < NGRP; i++ ) {

    // strip out OPTLIST option-string from () in SURVEYGROUP_LIST
    extractStringOpt(INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i],
		     INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_OPTLIST[i] ); 

    splitString(INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i], 
		PLUS, MXNUM_SAMPLE, // (I) 
		&INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[i], ptrTmp ); // (O)

    // store integer IDSURVEY for each plus-separated SURVEY 
    for(i2=0; i2 < INPUTS_SAMPLE_BIASCOR.NSURVEY_PER_GROUP[i]; i2++ ) {
      ID = get_IDSURVEY(ptrTmp[i2]);
      if( ID < 0 ) {
	sprintf(c1err,"Undefined SURVEY = '%s'", ptrTmp[i2] );
	sprintf(c2err,"check SURVEYGROUP = '%s' ", 
		INPUTS_SAMPLE_BIASCOR.SURVEYGROUP_LIST[i] );
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
      }

      // May 29 2020: abort if there is no data in this requested sample
      NEVT = INFO_DATA.TABLEVAR.NSN_PER_SURVEY[ID];
      if ( NEVT == 0 && INPUTS.surveyGroup_biasCor_abortFlag ) {
        sprintf(c1err,"No data for requested SURVEY = %s(%d)", ptrTmp[i2],ID);
        sprintf(c2err,"Check input key surveygroup_biascor");
        errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  print_debug_malloc(-1*debug_malloc,fnam);
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
  // Feb 9 2022: fix bug to allow [SURVEY] fieldgroups and also allow
  //             [SURVEY] in surveygroup_biascor. See new DO_COPY bool.

  int NSAMPLE = NSAMPLE_BIASCOR ;
  int debug_malloc = INPUTS.debug_malloc ;
  int ID, ID2, IDSURVEY, NSORT=0  ;
  int NCOPY[MXNUM_SAMPLE], IDMAP_SORT[MXNUM_SAMPLE];
  int igrp, NGRP_USR,  FLAG ; 
  char *s0, *s1, *s2, *f0, *f2, *NAME ;
  SAMPLE_INFO_DEF *SAMPLE_BIASCOR_TEMP ;
  bool ISGRP0, ISGRP2, SMATCH, NOGRP, DO_COPY ;
  bool LOGIC_FIX_FEB09 = true;
  int LDMP = 0 ;
  char fnam[] = "sort_IDSAMPLE_biasCor" ;
  // ---------------- BEGIN ---------------

  print_debug_malloc(+1*debug_malloc,fnam);
  SAMPLE_BIASCOR_TEMP = 
    (SAMPLE_INFO_DEF*) malloc ( NSAMPLE * sizeof(SAMPLE_INFO_DEF));
  

  // Jan 2020: abort if some SURVEY events are part of a FIELDGROUP,
  //           and some events are not.
  //   e..g, if DES and DES(C3+X3) are defined, ABORT.
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
        errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);

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
      s2      = SAMPLE_BIASCOR_TEMP[ID].NAME_SURVEYGROUP ;
      NAME    = SAMPLE_BIASCOR_TEMP[ID].NAME ;
      f0      = SAMPLE_BIASCOR_TEMP[ID].NAME_FIELDGROUP;  
      SMATCH  = strcmp(s1,s2)==0;
      NOGRP   = IGNOREFILE(f0);
      DO_COPY = SMATCH && NOGRP;
      if ( !LOGIC_FIX_FEB09 ) { DO_COPY = SMATCH; } // only for emergency

      if ( DO_COPY ) {
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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


  print_debug_malloc(-1*debug_malloc,fnam);
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
// xxx  S2->OPT_PHOTOZ = S1->OPT_PHOTOZ ;
  S2->IS_PHOTOZ = S1->IS_PHOTOZ ;


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
  
  fprintf(FP_STDOUT,"  SAMPLE_INFO DUMP for %s: \n", STRTYPE )    ;

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
    fprintf(FP_STDOUT,"  IDSAMPLE=%2d --> %-20.20s  (%s, %s)\n",
	   i, NAME, Nstring, zString );

    // printf(" xxx %s: IDSURVEY = %d \n", fnam, SAMPLE_BIASCOR[i].IDSURVEY);
  }

  if ( NZERR > 0 ) {
    sprintf(c1err,"zMAX < zMIN for %d IDSAMPLEs", NZERR);
    sprintf(c2err,"See IDSAMPLE dump above.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( ABORT_ON_NEVTZERO && NSAMPLE_ZERO ) {
    sprintf(c1err,"%d of %d SAMPLES have zero events.", NSAMPLE_ZERO,N);
    sprintf(c2err,"Check %s sample.", STRTYPE );
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
  }

  fflush(FP_STDOUT);

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
    fprintf(FP_STDOUT, "  fieldgroup_biasCor = '%s' \n", 
	    INPUTS.fieldGroup_biasCor);
    for(igroup=0; igroup < NFIELDGROUP ; igroup++ ) {
      FTMP1 = INPUTS_SAMPLE_BIASCOR.FIELDGROUP_LIST[igroup] ;
      fprintf(FP_STDOUT,"  FIELDGROUP(%d) = '%s' \n", igroup, FTMP1);
    }
    sprintf(c1err,"Invalid NMATCH_TOT=%d for CID=%s",   NMATCH_TOT, SNID );
    sprintf(c2err,"FIELD='%s' and FIELDGROUP list", FIELD);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  return ;

} // end match_surveyGroup


int get_IDSAMPLE(int IDSURVEY, bool IS_PHOTOZ, 
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
  bool  IS_PHOTOZ_TMP ;
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    IS_PHOTOZ_TMP     = SAMPLE_BIASCOR[i].IS_PHOTOZ ;

    // check matches
    MATCH_FIELDGROUP  = 0;
    MATCH_SURVEYGROUP = 0 ; 

    if ( CHECK_FIELDGROUP ) 
      { MATCH_FIELDGROUP  = (strcmp(FIELDGROUP_TMP,FIELDGROUP) == 0); }
    else if ( CHECK_SURVEYGROUP )
      { MATCH_SURVEYGROUP = (strcmp(SURVEYGROUP_TMP,SURVEYGROUP) == 0); }
   
    MATCH_PHOTOZ = (IS_PHOTOZ_TMP == IS_PHOTOZ ) ;


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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  // Jun 25 2020: if INPUTS.fitflag_sigmb=0, leave it at zero
  // Dec 21 2020: check option to NOT require valid biasCor 

  int INDX, IDSAMPLE, SKIP, NSN_DATA, CUTMASK, ievt ;
  int  NBINm = INPUTS.nbin_logmass;
  int  NBINg = 0; // is set below 
  int  OPTMASK        = INPUTS.opt_biasCor ;

  int  EVENT_TYPE     = EVENT_TYPE_BIASCOR;
  int  nfile_biasCor  = INPUTS.nfile_biasCor ;
  bool  DOCOR_1DZAVG  = ( OPTMASK & MASK_BIASCOR_1DZAVG  );
  bool  DOCOR_1DZWGT  = ( OPTMASK & MASK_BIASCOR_1DZWGT  );
  bool  DOCOR_1D5DCUT = ( OPTMASK & MASK_BIASCOR_1D5DCUT );
  bool  DOCOR_MUCOVSCALE   = ( OPTMASK & MASK_BIASCOR_MUCOVSCALE);
  bool  DOCOR_MUCOVADD   = ( OPTMASK & MASK_BIASCOR_MUCOVADD);
  bool  DOCOR_5D      = ( OPTMASK & MASK_BIASCOR_5D      );
  bool  DOCOR_1D      = ( DOCOR_1DZAVG || DOCOR_1DZWGT || DOCOR_1D5DCUT);
  bool  IDEAL         = ( OPTMASK & MASK_BIASCOR_COVINT ) ;
  bool  DOCOR_MU      = ( OPTMASK & MASK_BIASCOR_MU ) ;
  bool  REQUIRE_VALID_BIASCOR = (OPTMASK & MASK_BIASCOR_noCUT) == 0 ;

  char txt_biasCor[40], *name  ;
  
  bool USEDIM_GAMMADM, USEDIM_LOGMASS;
  int NDIM_BIASCOR=0, ILCPAR_MIN, ILCPAR_MAX ;
  char fnam[] = "prepare_biasCor" ;

  // ------------- BEGIN -------------
  
  INFO_BIASCOR.NDIM = 0;
  INFO_BIASCOR.TABLEVAR.NSN_ALL       = 0 ;
  INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS  = 0 ;
  INFO_BIASCOR.TABLEVAR.NSN_REJECT    = 0 ;
  INFO_BIASCOR.GAMMADM_OFFSET         = 0.0 ;
  INFO_BIASCOR.DUST_FLAG          = false ;

  NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ;

  if (DOCOR_MUCOVADD) { // link MAD with MUCOVADD; doesnt work without MAD
    INPUTS.opt_biasCor |= MASK_BIASCOR_MAD;
  }

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

    if ( DOCOR_MUCOVSCALE ) 
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
  }

  if ( DOCOR_1D && DOCOR_5D ) {
    sprintf(c1err,"Cannot do 1D and 5D biasCor.");
    sprintf(c2err,"opt_biascor logic is messed up.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }
  
  fprintf(FP_STDOUT, "\n\t Use simulaton to interpolate %s biasCor\n", 
	 txt_biasCor);
  fprintf(FP_STDOUT, "\t\t (opt_biascor=%d, sigma_cell=%.2f)\n", 
	 OPTMASK, INPUTS.sigma_cell_biasCor);

  if ( !REQUIRE_VALID_BIASCOR ) {
    fprintf(FP_STDOUT,"\t Keep events without valid biasCor. \n");
  }

  fflush(FP_STDOUT);

  // -------------------------------------------
  int PS0 = INPUTS.prescale_biasCor[0];
  int PS1 = INPUTS.prescale_biasCor[1];
  if ( PS1 > 1 ) {
    fprintf(FP_STDOUT, "\t Apply pre-scale: select subset %d of %d \n", 
	    PS0, PS1);  
  }

  // Jan 16 2018: force correct option for log(sigma) term in chi2
  if ( INPUTS.fitflag_sigmb > 0 ) {
    if ( NDIM_BIASCOR >=5  ) 
      { INPUTS.fitflag_sigmb = 2 ; } // include log(sigma) term
    else
      { INPUTS.fitflag_sigmb = 1 ; } // ignore log(sigma) term
    fprintf(FP_STDOUT, "\t Force fitflag_sigmb = %d \n", INPUTS.fitflag_sigmb);
  }

  // check which cov parameter in fit
  if ( IDEAL )  { 
    fprintf(FP_STDOUT, "\t Compute COVINT(z) using IDEAL fit params.\n");
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ;
  }
  

  // count number of biasCor events passing cuts (after setup_BININFO calls)
  for(ievt=0; ievt < INFO_BIASCOR.TABLEVAR.NSN_ALL; ievt++ )  { 
    compute_more_TABLEVAR(ievt, &INFO_BIASCOR.TABLEVAR );
    compute_CUTMASK(ievt, &INFO_BIASCOR.TABLEVAR );
  }

  if ( NDIM_BIASCOR >=5 ) { store_iaib_biasCor(); }
  

  print_eventStats(EVENT_TYPE);

  //  if ( NDIM_BIASCOR == 1 ) { goto CHECK_1DCOR ; }
  if ( NDIM_BIASCOR == 1 && !DOCOR_1D5DCUT ) { goto CHECK_1DCOR ; }

  // make sparse list for faster looping below (Dec 21 2017)
  makeSparseList_biasCor();

  // determine sigInt for biasCor sample BEFORE makeMap since
  // sigInt is needed for 1/muerr^2 weight
  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR ; IDSAMPLE++ )  {  
    init_sigInt_biasCor_SNRCUT(IDSAMPLE);
  } 

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

    fprintf(FP_STDOUT, "\n %s\n", borderLine);   

    // print ISAMPLE info to screen
    NAME_SAMPLE  = SAMPLE_BIASCOR[IDSAMPLE].NAME ; 

    SKIP = 0 ;
    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_SELECT  == 0 ) { SKIP=1; }
    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { SKIP=1; }

    if ( SKIP ) {
      fprintf(FP_STDOUT, "\t SKIP BIASCOR PREP for %s\n", NAME_SAMPLE);
      continue ;
    }
    else  {
      fprintf(FP_STDOUT, "\t START BIASCOR PREP for SAMPLE = %s\n", 
	      NAME_SAMPLE);
    }
    fflush(FP_STDOUT);

    // get wgted avg in each bin to use for interpolation 
    makeMap_binavg_biasCor(IDSAMPLE);

    // prepare 3D bias maps need to interpolate bias
    for(INDX = ILCPAR_MIN; INDX <= ILCPAR_MAX ; INDX++ ) 
      { makeMap_fitPar_biasCor(IDSAMPLE,INDX); }

    fprintf(FP_STDOUT,"\n");  fflush(FP_STDOUT);   

    // make map of sigma_mu bias
    if ( DOCOR_MUCOVSCALE ) { makeMap_sigmu_biasCor(IDSAMPLE); }

    fprintf(FP_STDOUT, "\n\t END BIASCOR PREP for %s\n", NAME_SAMPLE);
    fprintf(FP_STDOUT, " %s\n", borderLine);       
    fflush(FP_STDOUT);
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

  fprintf(FP_STDOUT, "\n");
  fprintf(FP_STDOUT, " For each DATA event before fit, store \n");
  fprintf(FP_STDOUT, "   * mb,x1,c-bias at each alpha,beta,gammaDM \n");
  if ( DOCOR_MUCOVSCALE) 
    { fprintf(FP_STDOUT, "   * muCOVscale   at each alpha,beta,gammaDM \n"); }


  int DUMPFLAG = 0 ;
  int ndump_nobiasCor = INPUTS.ndump_nobiasCor;

  /* xxxx mark delete Apr 20 2022 xxxxxxxxxx
  // user input ndump_nobiasCor is number of noBiasCor events to dump;
  // however, dump at least 10 if user uses this input as a logic flag.
  int ndump_nobiasCor = INPUTS.ndump_nobiasCor;
  if ( ndump_nobiasCor>0 && ndump_nobiasCor < 10 ) { ndump_nobiasCor=10; }
  xxxxxxx end mark xxxxxx */

  for (n=0; n < NSN_DATA; ++n) {

    CUTMASK  = INFO_DATA.TABLEVAR.CUTMASK[n]; 
    IDSAMPLE = INFO_DATA.TABLEVAR.IDSAMPLE[n]; 
    if ( CUTMASK ) { continue ; }

    name = INFO_DATA.TABLEVAR.name[n];
    DUMPFLAG = ( strstr(INPUTS.cidlist_debug_biascor,name) != NULL );

    istore = storeDataBias(n,DUMPFLAG);
    
    NUSE[IDSAMPLE]++ ; NUSE_TOT++ ;
    if ( istore == 0 && REQUIRE_VALID_BIASCOR )  { 
      NSKIP_TOT++; NSKIP[IDSAMPLE]++ ;  
      setbit_CUTMASK(n, CUTBIT_BIASCOR, &INFO_DATA.TABLEVAR); 
      if( NSKIP_TOT < ndump_nobiasCor ) 
	{ storeDataBias(n,1); } 
    }    
  }


  for(IDSAMPLE=0; IDSAMPLE < NSAMPLE_BIASCOR; IDSAMPLE++ ) {

    // store NSKIP and NUSE for later to print warning about excess loss
    NDATA_BIASCORCUT[0][IDSAMPLE] = NUSE[IDSAMPLE];
    NDATA_BIASCORCUT[1][IDSAMPLE] = NSKIP[IDSAMPLE];

    if ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_SELECT == 0 ) { continue; }
    fprintf(FP_STDOUT,
	    "  Rejected %d (of %d) DATA events with no %dD biasCor (%s).\n",
	   NSKIP[IDSAMPLE], NUSE[IDSAMPLE], 
	   NDIM_BIASCOR, SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }
  fprintf(FP_STDOUT, "\n");
  fflush(FP_STDOUT);


  // ------------
  // free memory used to hold each simBias event;
  // we only need to keep the small map made in makeMap_fitPar_biasCor.

  if ( NSKIP_TOT == NUSE_TOT ) {
    sprintf(c1err,"Could not compute biasCor");
    sprintf(c2err,"Something is messed up.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  int NROW, ISTART, IFILETYPE, ifile, LEN_MALLOC ;   
  int NEVT[MXFILE_BIASCOR], NEVT_TOT, NVAR_ORIG ;
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
    fprintf(FP_STDOUT, "\t Found %d events in %s. \n", NEVT[ifile], simFile);
    fflush(FP_STDOUT);
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
  // Special biasCor option to correct vs. redshift only.
  // Evaluate mubias(sample,z) from biasCor sample.
  //
  // Jun 25 2020: fix determination of mu_true based on DO_BIASCOR_MU
  // Jul 02 2020: fix mu_true = SIM_MUz, regardless of DO_BIASCOR_MU

  int    NSAMPLE  = NSAMPLE_BIASCOR ;
  double alpha    = INPUTS.parval[IPAR_ALPHA0] ;
  double beta     = INPUTS.parval[IPAR_BETA0] ;
  //  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );

  int    NSN_DATA, NSN_BIASCOR, ievt, NCUTS, iz, izusr, idsample, cutmask ;
  double mu_fit, mB_fit, x1_fit, c_fit ;
  double mu_sim, mB_sim, x1_sim, c_sim ;
  double mB_err, x1_err, c_err ;
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


  fprintf(FP_STDOUT, "\n");
  fprintf(FP_STDOUT, "   Compute muBias vs. z with: \n" );
  fprintf(FP_STDOUT, "\t    alpha=%.3f(p1) beta=%.3f(p2) \n", alpha, beta);
  //  fprintf(FP_STDOUT,"\t    zbinSize = %.3f \n", zbinSize );
  fflush(FP_STDOUT);

    
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

    mu_true   = (double)INFO_BIASCOR.TABLEVAR.SIM_MUz[ievt];  // at zHD
    mu_fit    = mB_fit + alpha*x1_fit - beta*c_fit - M0_DEFAULT ;
    mu_sim    = mB_sim + alpha*x1_sim - beta*c_sim - M0_DEFAULT ;

    muerrLens = ( z * INPUTS.lensing_zpar ) ;

    iz = (int)( (z-zMIN[idsample]) / zbinSize[idsample] ) ;
    
    if ( iz >= MXz ) {
      sprintf(c1err,"iz=%d exceeds bound for z=%.4f", iz, z);
      sprintf(c2err,"Check MAXBIN_z in SALT2mu.c");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

  fprintf(FP_STDOUT, "\t %d of %d biasCor events pass cuts. \n",  
	 NCUTS, NSN_BIASCOR);
  INFO_BIASCOR.TABLEVAR.NSN_PASSCUTS = NCUTS; 

  if ( NCUTS < 10 ) {
    print_eventStats(EVENT_TYPE_BIASCOR);
    sprintf(c1err,"NCUTS=%d is too few.", NCUTS );
    sprintf(c2err,"Something is wrong.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  // --------------------------------
  // get MU bias in each z-bin
  // Also store IZMIN, IZMAX and NBZ for interpolation below
  int izmin, izmax;
  int IZMIN[MXNUM_SAMPLE];
  int IZMAX[MXNUM_SAMPLE];
  int NBZ[MXNUM_SAMPLE];

  MUBIAS_SIM_HISNR /= SUMWGT_HISNR ;
  fprintf(FP_STDOUT, "\t muBias(sim,SNR>%.0f) = %6.3f \n",
	 INPUTS.snrmin_sigint_biasCor, MUBIAS_SIM_HISNR);

  for(idsample=0; idsample < NSAMPLE; idsample++ ) {

    if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { continue; }
    IZMIN[idsample] = +999;
    IZMAX[idsample] = -999;
    
    fprintf(FP_STDOUT, "     %s : \n", SAMPLE_BIASCOR[idsample].NAME );
    for(iz=0; iz < MXz; iz++ ) {

      if ( SUMWGT[idsample][iz] <= 1.0E-9 ) { continue ; }

      if ( iz < IZMIN[idsample] ) { IZMIN[idsample] = iz; }
      if ( iz > IZMAX[idsample] ) { IZMAX[idsample] = iz; }
      
      MUBIAS_FIT[idsample][iz] /= SUMWGT[idsample][iz] ;
      MUBIAS_SIM[idsample][iz] /= SUMWGT[idsample][iz] ;

      // subtract muBias for hi-SNR so that muBias ~ 0 at low-z.
      MUBIAS_FIT[idsample][iz] -= MUBIAS_SIM_HISNR ;
      MUBIAS_SIM[idsample][iz] -= MUBIAS_SIM_HISNR ;

      fprintf(FP_STDOUT, "\t z=%.3f - %.3f : muBias = %6.3f  (iz=%d)\n",
	     zlo[idsample][iz], zhi[idsample][iz],
	     MUBIAS_FIT[idsample][iz], iz );
      
    }  // end iz
    fflush(FP_STDOUT);

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
  fprintf(FP_STDOUT, "\n   Now compute muBias(z) for each data event. \n");
  fflush(FP_STDOUT);

  double muBias ;

  for(ievt=0; ievt < NSN_DATA; ievt++ ) {

    INFO_DATA.muBias_zinterp[ievt]  = 0.0;

    idsample = INFO_DATA.TABLEVAR.IDSAMPLE[ievt];
    cutmask = INFO_DATA.TABLEVAR.CUTMASK[ievt];
    if ( cutmask ) { continue ; }

    if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) 
      { setbit_CUTMASK(ievt, CUTBIT_BIASCOR, &INFO_DATA.TABLEVAR); }

    z        = INFO_DATA.TABLEVAR.zhd[ievt] ;
    //    mB       = INFO_DATA.TABLEVAR.fitpar[INDEX_mB][ievt];
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
    fprintf(FP_STDOUT,"\n");
    fprintf(FP_STDOUT,
	    " ):--):--):--):--):--):--):--):--):--):--):--):--):--):--\n") ;
    fprintf(FP_STDOUT,
	    " %s WARNING: %d data events with invalid z rejected:\n",
	   fnam, NDATA_REJECT );
    for(idsample=0; idsample < NSAMPLE; idsample++ ) {
      n = ndata[idsample] ;  nrej = ndata_reject[idsample] ;
      if ( nrej > 0 ) {
	ratio = 100.0 * (float)nrej / (float)n;
	fprintf(FP_STDOUT, "\t NREJECT(%s) = %d (%.2f) \n",
	       SAMPLE_BIASCOR[idsample].NAME, nrej, ratio );  
      }
    }
    fprintf(FP_STDOUT,
	    " ):--):--):--):--):--):--):--):--):--):--):--):--):--):--\n");
    fflush(FP_STDOUT);
  }

    
  return ;

} // end prepare_biasCor_zinterp



// ================================================================
void set_MAPCELL_biasCor(int IDSAMPLE) {

  
  // set  CELLINFO_BIASCOR.MAPCELL[ia][ib][ig][iz][im][ix1][ic]
  // to map 5D indices to 1D index J1D.
  // July 31 2019: add logmass dimension [im]

  int ID = IDSAMPLE;
  int debug_malloc = INPUTS.debug_malloc ;
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  CELLINFO_BIASCOR[IDSAMPLE].NCELL = NCELL ;

  NDIM = 5;
  if ( NBINg > 1 ) { NDIM++ ; }
  if ( NBINm > 1 ) { NDIM++ ; }

  fprintf(FP_STDOUT,"  %s : malloc for %d  %dD-cells for biasCor \n", 
	 fnam, NCELL, NDIM );
  fflush(FP_STDOUT);

  // malloc other CELLINFO arrays
  int MEMD0   = NCELL  * sizeof(double);
  int MEMI0   = NCELL  * sizeof(int);

  print_debug_malloc(+1*debug_malloc,fnam);
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
  INFO_BIASCOR.MUCOVADD[IDSAMPLE] = (float *        ) malloc ( MEMCOV );

  fflush(FP_STDOUT);
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
  int debug_malloc = INPUTS.debug_malloc ;

  int MEMD, NEVT_USE, NEVT_SAMPLE,  LDMP, N ;
  int J1D, ia, ib, ig, iz, im, ix1, ic, ievt, isp ;
  int J1DNBR_LIST[MXJ1DNBR], NJ1DNBR ;

  double fit_val, sim_val, biasVal, sim_gammadm ;
  double WGT, VAL, ERR, RMS, SQRMS, XN, XNLIST, tmp1, tmp2 ;
  double *SUMBIAS, *SUMWGT, *sum, *sumsq ;

  float *ptr_fitpar, *ptr_simpar, *ptr_gammadm;
  char fnam[] = "makeMap_fitPar_biasCor";

  // ----------------- BEGIN -----------------

  fprintf(FP_STDOUT, "  %s of %s-bias(%s)  \n", 
	 fnam, PARNAME, INFO_BIASCOR.STRING_PARLIST);
  fflush(FP_STDOUT);

  // malloc arrays to store info in each biasCor cell
  print_debug_malloc(+1*debug_malloc,fnam);
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
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
      }
      NgetList++ ;
    }

    sumsq_nbr = sum_nbr = 0.0 ;  Nsum_nbr=0;
    for(inbr=0; inbr < NJ1DNBR; inbr++ ) {
      J1D_nbr  = J1DNBR_LIST[inbr] ;
      
      if ( J1D_nbr < 0 || J1D_nbr >= NCELL ) {
	sprintf(c1err,"Invalid J1D_nbr=%d for inbr=%d", J1D_nbr, inbr);
	sprintf(c2err,"J1D=%d", J1D);
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
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
    fprintf(FP_STDOUT,
	   "  BiasCor computed for %d of %d grid-cells with >=1 events.\n",
	   NMAP_USE, NMAP_TOT ) ;
    fprintf(FP_STDOUT,
	   "  BiasCor sample: %d of %d pass cuts for IDSAMPLE=%d.\n",
	   NEVT_USE, NEVT_SAMPLE, IDSAMPLE );

    if ( NEVT_USE == 0 ) {
      print_eventStats(EVENT_TYPE_BIASCOR);
      sprintf(c1err,"No BiasCor events passed for %s", 
	      SAMPLE_BIASCOR[IDSAMPLE].NAME );
      sprintf(c2err,"Check BiasCor file" );
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
    }

    fflush(FP_STDOUT);
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

  fflush(FP_STDOUT);

  print_debug_malloc(-1*debug_malloc,fnam);
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
  // Sep 27 2021: check INPUTS.interp_biascor_logmass

  int  USEMASK = USEMASK_BIASCOR_COVTOT + USEMASK_BIASCOR_ZMUERR;
  int  istat_cov, J1D, idsample ;
  double WGT, muerrsq ;
  char fnam[] = "WGT_biasCor" ;
  
  // --------------- BEGIN -------------------

  muerrsq  = muerrsq_biasCor(ievt, USEMASK, &istat_cov) ; 
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
  if ( NBINm > 1 && INPUTS.interp_biascor_logmass ) {
    m          = (double)INFO_BIASCOR.TABLEVAR.host_logmass[ievt] ;
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
  m        = (double)INFO_BIASCOR.TABLEVAR.host_logmass[ievt];
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

  // May 2016: for input J1D index, return int indices for a,b,z,x1,c

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
  // Make map of MUCOVSCALE, where MUCOVSCALE = RMS(MURES)/AVG(MUERR)
  // as a function of (z,x1,c,a,b).
  //
  // To get better stats, the bias is just a function of z & c,
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
  // Jun 29 2021: fix bug setting im index for logmass.
  // Sep 14 2021: little cleanup/refac 
  // Sep 16 2021: add dump utils; see i1d_dump_mucovscale and OPTMASK
  // Jun 05 2022: write SALT2 fit params in abort msg for crazy muErr

  int NBIASCOR_CUTS    = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;
  int NBIASCOR_ALL     = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  bool DO_MAD          = (INPUTS.opt_biasCor & MASK_BIASCOR_MAD) > 0;
  int debug_malloc     = INPUTS.debug_malloc ;
  int debug_mucovscale = INPUTS.debug_mucovscale ;

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  int    NBINa, NBINb, NBINg, NBINz, NBINm, NBINc, NperCell ;
  int    OPTMASK, DUMPFLAG = 0 ;
  int    ia, ib, ig, iz, im, ic, i1d, NCELL, isp ; 
  int    ievt, istat_cov, istat_bias, N, J1D, ipar, USEMASK ;
  double muErr, muErrsq, muErrsq_raw, muDif, muDifsq, pull, tmp1, tmp2  ;
  double muBias, muBiasErr, muCOVscale, muCOVadd, fitParBias[NLCPAR+1] ;
  double a, b, gDM, z, m, c ;
  double *SUM_MUERR, *SUM_SQMUERR;
  double *SUM_MUDIF, *SUM_SQMUDIF ;
  double *SQMUERR,   *SQMUSTD ;
  double *SUM_PULL,  *SUM_SQPULL ;
  double *SIG_PULL_MAD; //1.48*MedianAbsDev
  double *SIG_PULL_STD;
  char   *name ;

  int NperCell_min = MINPERCELL_MUCOVSCALE;

  // Declare lists for debug_mucovscale
  double *muErr_list, *muErr_raw_list, *muDif_list;
  int NPERCELL_REALLOC=2000;
  int N_REALLOC=0;

  double    UNDEFINED = 9999.0 ;
  float    *ptr_MUCOVSCALE;
  float    *ptr_MUCOVADD;
  BIASCORLIST_DEF     BIASCORLIST ;
  FITPARBIAS_DEF      FITPARBIAS[MXa][MXb][MXg] ;
  double              MUCOVSCALE[MXa][MXb][MXg] ;
  double              MUCOVADD[MXa][MXb][MXg] ;
  INTERPWGT_AlphaBetaGammaDM INTERPWGT ;
 
  CELLINFO_DEF *CELL_BIASCOR    = &CELLINFO_BIASCOR[IDSAMPLE];
  CELLINFO_DEF *CELL_MUCOVSCALE = &CELLINFO_MUCOVSCALE[IDSAMPLE];
  CELLINFO_DEF *CELL_MUCOVADD   = &CELLINFO_MUCOVADD[IDSAMPLE];

  char fnam[]  = "makeMap_sigmu_biasCor" ;
  
  // ----------------- BEGIN -------------------

  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return; }

  if ( debug_mucovscale > 0 ) {
    int memd   = sizeof(double) * NBIASCOR_ALL;
    int memi   = sizeof(int   ) * NBIASCOR_ALL;
    if ( IDSAMPLE == 0 ) 
      { INFO_BIASCOR.TABLEVAR.IMUCOV = (int   *) malloc(memi); }
    muErr_list                = (double*) malloc(memd);
    muErr_raw_list            = (double*) malloc(memd);
    muDif_list                = (double*) malloc(memd);
  }

  fprintf(FP_STDOUT, " %s: make map of muCOVscale(a,b,g,z,c) \n", fnam);
  fflush(FP_STDOUT);

  malloc_MUCOV(+1,IDSAMPLE, CELL_MUCOVSCALE );
  malloc_MUCOV(+1,IDSAMPLE, CELL_MUCOVADD   );
  NCELL = CELL_MUCOVSCALE->NCELL;

  int MEMD     = NCELL   * sizeof(double);
  int MEMI     = NCELL   * sizeof(int);

  ptr_MUCOVSCALE = INFO_BIASCOR.MUCOVSCALE[IDSAMPLE];         
  ptr_MUCOVADD   = INFO_BIASCOR.MUCOVADD[IDSAMPLE];         

  NBINa    = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  NBINb    = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;  
  NBINg    = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;  
  NBINz    = CELL_MUCOVSCALE->BININFO_z.nbin ;
  NBINm    = CELL_MUCOVSCALE->BININFO_m.nbin ;  
  NBINc    = CELL_MUCOVSCALE->BININFO_LCFIT[INDEX_c].nbin;
 
  // malloc local 1D arrays to track local sums.
  print_debug_malloc(+1*debug_malloc,fnam);
  SUM_MUERR      = (double*) malloc(MEMD);
  SUM_SQMUERR    = (double*) malloc(MEMD);
  SUM_MUDIF      = (double*) malloc(MEMD);
  SUM_SQMUDIF    = (double*) malloc(MEMD);
  SQMUERR        = (double*) malloc(MEMD);
  SQMUSTD        = (double*) malloc(MEMD);
  SUM_PULL       = (double*) malloc(MEMD);
  SUM_SQPULL     = (double*) malloc(MEMD);
  SIG_PULL_MAD   = (double*) malloc(MEMD);
  SIG_PULL_STD   = (double*) malloc(MEMD);
  ig = 0 ;


  int N1D=0;
  for(ia=0; ia< NBINa; ia++ ) {
    for(ib=0; ib< NBINb; ib++ ) {  
      for(ig=0; ig< NBINg; ig++ ) {  
	MUCOVSCALE[ia][ib][ig] = 1.0 ; // dummy arg for get_muBias below
	MUCOVADD[ia][ib][ig] = 1.0 ; // dummy arg for get_muBias below
	for(iz=0; iz < NBINz; iz++ ) {
	  for(im=0; im < NBINm; im++ ) {
	    for(ic=0; ic < NBINc; ic++ ) {
	      SUM_MUERR[N1D] = SUM_SQMUERR[N1D] = 0.0 ;
	      SUM_MUDIF[N1D] = SUM_SQMUDIF[N1D] = 0.0 ;	
	      SUM_PULL[N1D]  = SUM_SQPULL[N1D]  = 0.0 ;
              SIG_PULL_MAD[N1D] = UNDEFINED ;
	      SIG_PULL_STD[N1D] = UNDEFINED ;

	      ptr_MUCOVSCALE[N1D] = 1.0 ;
	      CELL_MUCOVSCALE->NperCell[N1D]  = 0 ;
	      CELL_MUCOVSCALE->AVG_z[N1D]     = 0.0 ;
	      CELL_MUCOVSCALE->AVG_m[N1D]     = 0.0 ;
	      CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][N1D] = 0.0 ;
	      CELL_MUCOVSCALE->MAPCELL[ia][ib][ig][iz][im][0][ic]=N1D;

	      ptr_MUCOVADD[N1D] = 1.0e-12 ;
	      CELL_MUCOVADD->NperCell[N1D]  = 0 ;
	      CELL_MUCOVADD->AVG_z[N1D]     = 0.0 ;
	      CELL_MUCOVADD->AVG_m[N1D]     = 0.0 ;
	      CELL_MUCOVADD->AVG_LCFIT[INDEX_c][N1D] = 0.0 ;
	      CELL_MUCOVADD->MAPCELL[ia][ib][ig][iz][im][0][ic]=N1D;

	      N1D++ ;
	    }	  
	  }
	}
      }
    }
  }

  for(isp=0; isp < NBIASCOR_CUTS; isp++ ) {

    ievt = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;

    if ( debug_mucovscale > 0 ) { INFO_BIASCOR.TABLEVAR.IMUCOV[ievt] = -9;  }

    // check if there is valid biasCor for this event
    J1D = J1D_biasCor(ievt,fnam);
    if ( CELL_BIASCOR->NperCell[J1D] < BIASCOR_MIN_PER_CELL ) 
      { continue ; } 

    for(ia=0; ia<MXa; ia++ ) {
      for(ib=0; ib<MXb; ib++ ) {
	for(ig=0; ig<MXg; ig++ ) {
	  zero_FITPARBIAS(&FITPARBIAS[ia][ib][ig] ); 
	}
      }
    }
    
    z    = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt];
    m    = (double)INFO_BIASCOR.TABLEVAR.host_logmass[ievt];
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

    // allow color (c) and logmass to be outside map
    iz = IBINFUN(z, &CELL_MUCOVSCALE->BININFO_z, 
		 1, fnam );

    im = IBINFUN(m, &CELL_MUCOVSCALE->BININFO_m, 
		 2, fnam );

    ic = IBINFUN(c, &CELL_MUCOVSCALE->BININFO_LCFIT[INDEX_c], 
		 2, fnam );

    // ---------------------------------------------------
    // need bias corrected distance to compute pull

    BIASCORLIST.z            = z ;
    BIASCORLIST.host_logmass = m ;
    BIASCORLIST.alpha        = a ;
    BIASCORLIST.beta         = b ;
    BIASCORLIST.gammadm      = gDM ;
    BIASCORLIST.idsample     = IDSAMPLE ;

    istat_bias = 
      get_fitParBias(name, &BIASCORLIST, DUMPFLAG, fnam, 
		     &FITPARBIAS[ia][ib][ig] ); // <== returned

    //  DUMPFLAG = 0 ; // xxx REMOVE

    // skip if bias cannot be computed, just like for data
    if ( istat_bias <= 0 ) { continue ; }

    get_INTERPWGT_abg(a,b,gDM, DUMPFLAG, &INTERPWGT, fnam );
    get_muBias(name, &BIASCORLIST, FITPARBIAS,MUCOVSCALE,MUCOVADD, &INTERPWGT,
	       fitParBias, &muBias, &muBiasErr, &muCOVscale, &muCOVadd );  

    // ----------------------------
    muDif   =  muresid_biasCor(ievt);  // mu - muTrue
    muDif  -=  muBias ;  
    muDifsq =  muDif*muDif ;

    // compute error with intrinsic scatter 
    USEMASK = USEMASK_BIASCOR_COVTOT;
    muErrsq = muerrsq_biasCor(ievt, USEMASK, &istat_cov) ; 

    if ( muErrsq <= 1.0E-14 || muErrsq > 100.0 || isnan(muErrsq) ) {
      print_preAbort_banner(fnam);
      printf("\t z=%f  a=%f  b=%f  gDM=%f\n",
	     z, a, b, gDM);
      printf("\t ia,ib,ig = %d, %d, %d \n", ia, ib, ig);
      printf("\t istat_cov = %d \n", istat_cov);

      for(ipar=0; ipar < NLCPAR; ipar++ ) { 
	char *name = BIASCOR_NAME_LCFIT[ipar];
	float val  = INFO_BIASCOR.TABLEVAR.fitpar[ipar][ievt]; 
	float err  = INFO_BIASCOR.TABLEVAR.fitpar_err[ipar][ievt]; 
	printf("\t %3s = %f +_ %f \n", name, val, err); 
	fflush(stdout);
      }

      sprintf(c1err,"Invalid muErrsq=%f for ievt=%d (SNID=%s)", 
	      muErrsq, ievt, name );
      sprintf(c2err,"Something is messed up.");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
    }

    muErr   = sqrt(muErrsq) ;    
    pull    = (muDif/muErr) ;

    // get 1d index
    i1d = CELL_MUCOVSCALE->MAPCELL[ia][ib][ig][iz][im][0][ic] ;

    if ( debug_mucovscale > 0 ) {
      INFO_BIASCOR.TABLEVAR.IMUCOV[ievt] = i1d;
      muDif_list[ievt] = muDif;
      muErr_list[ievt] = muErr;

      // store raw muErr (RK, 9.14.2021)
      USEMASK     = USEMASK_BIASCOR_COVFIT;
      muErrsq_raw = muerrsq_biasCor(ievt, USEMASK, &istat_cov) ; 
      muErr_raw_list[ievt] = sqrt(muErrsq_raw);
    }

    SUM_PULL[i1d]    += pull ;
    SUM_SQPULL[i1d]  += (pull*pull) ;
    SUM_SQMUDIF[i1d] +=  muDifsq ;
    SUM_MUDIF[i1d]   +=  muDif ;
    SUM_SQMUERR[i1d] +=  muErrsq ;
    SUM_MUERR[i1d]   +=  muErr ;

    if (DO_MAD) {
      NperCell = CELL_MUCOVSCALE->NperCell[i1d];
      int  newmem = (NperCell+1+NPERCELL_REALLOC) * sizeof(double);
      bool DO_REALLOC = (NperCell+1)%NPERCELL_REALLOC == 0 && NperCell > 0;
      CELL_MUCOVSCALE->ABSPULL[i1d][NperCell] = fabs(pull);
      if ( DO_COVADD ) {
	CELL_MUCOVADD->MURES[i1d][NperCell] = muDif;
	CELL_MUCOVADD->MUCOV[i1d][NperCell] = muErrsq;
      }
      if ( DO_REALLOC ){
	CELL_MUCOVSCALE->ABSPULL[i1d] = 
	  (double *)realloc(CELL_MUCOVSCALE->ABSPULL[i1d], newmem);
	if (DO_COVADD) {
	  CELL_MUCOVADD->MURES[i1d] = 
	    (double *)realloc(CELL_MUCOVADD->MURES[i1d], newmem);
	  CELL_MUCOVADD->MUCOV[i1d] = 
	    (double *)realloc(CELL_MUCOVADD->MUCOV[i1d], newmem);
        }
	N_REALLOC++;
      } // end DO_REALLOC
    } // end USE_MAD_MUCOVSCALE


    // increment sums to get average in each cell   
    CELL_MUCOVSCALE->NperCell[i1d]++ ;
    CELL_MUCOVSCALE->AVG_z[i1d]              += z ;
    CELL_MUCOVSCALE->AVG_m[i1d]              += m ;
    CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d] += c ;

  } // end ievt

  //printf("xxx N_REALLOC=%d\n",N_REALLOC);

  // -------------------------------------------------
  double XN, SQSTD, AVG=-9.0, STD=-9.0, MAD=-9.0 ;
  
  for(i1d=0; i1d < NCELL; i1d++ ) {

    SQMUSTD[i1d] = 0.0 ;
    SQMUERR[i1d] = 0.0 ;
    
    N  = CELL_MUCOVSCALE->NperCell[i1d] ;
    XN = (double)N ;

    if ( N < NperCell_min ) {       
      CELL_MUCOVSCALE->USE[i1d]             = false;

      if ( INPUTS.restore_mucovscale_bug ) { 
	CELL_MUCOVSCALE->USE[i1d]=true;
      }
      else {
	CELL_MUCOVSCALE->AVG_z[i1d]           = UNDEFINED ;
	CELL_MUCOVSCALE->AVG_m[i1d]           = UNDEFINED ;
	CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d] = UNDEFINED ;
      }
      continue ; 
    }
    
    CELL_MUCOVSCALE->USE[i1d]                 = true;
    CELL_MUCOVSCALE->AVG_z[i1d]              /= XN ;
    CELL_MUCOVSCALE->AVG_m[i1d]              /= XN ;
    CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d] /= XN ;

    tmp1= SUM_MUDIF[i1d]/XN ;    tmp2= SUM_SQMUDIF[i1d]/XN ;
    SQSTD = tmp2 - tmp1*tmp1 ;
    SQMUSTD[i1d] = SQSTD ;

    // average calculated muErrsq
    muErr        = SUM_MUERR[i1d]/XN ;
    muErrsq      = muErr*muErr  ;
    SQMUERR[i1d] = muErrsq ;
    
    // RMS of pull
    tmp1= SUM_PULL[i1d]/XN;  tmp2=SUM_SQPULL[i1d]/XN ;
    SQSTD = tmp2 - tmp1*tmp1 ;

    SIG_PULL_STD[i1d] = sqrt(SQSTD);

    if ( DO_MAD ) {
      // beware that MAD is meaningfill only for |PULL|; ignore AVG and STD
      arrayStat( N, CELL_MUCOVSCALE->ABSPULL[i1d], &AVG, &STD, &MAD);
      SIG_PULL_MAD[i1d]   = 1.48 * MAD;
      ptr_MUCOVSCALE[i1d] = (float)(SIG_PULL_MAD[i1d]*SIG_PULL_MAD[i1d]) ;
    } 
    else {
      ptr_MUCOVSCALE[i1d] = (float)(SIG_PULL_STD[i1d]*SIG_PULL_STD[i1d]) ;
    }

    if ( DO_COVADD ) {
      double sigInt;
      char   callfun[100];
      sprintf(callfun,"%s(j1d=%d,z=%.2f,c=%.2f,m=%.2f,IDSAMPLE=%d)",
	      fnam, i1d,
	      CELL_MUCOVSCALE->AVG_z[i1d],
	      CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d],
	      CELL_MUCOVSCALE->AVG_m[i1d],
	      IDSAMPLE);
      OPTMASK = 1;  // 1 --> do NOT abort if sigInt < 0

      if ( i1d == INPUTS.debug_mucovscale ) {
	OPTMASK += 64;
	printf(" xxx %s: ========================================= \n", fnam);
	printf(" xxx %s: ptr_MUCOVSCALE[%d] = %f   NperCell=%d\n", 
	       fnam, i1d, ptr_MUCOVSCALE[i1d], N );
	printf(" xxx %s: SIG_PULL_[STD,MAD][%d] = %f, %f \n",
	       fnam, i1d, SIG_PULL_STD[i1d], SIG_PULL_MAD[i1d] );
      }
      sigInt =
	sigint_muresid_list(N, CELL_MUCOVADD->MURES[i1d],
			    CELL_MUCOVADD->MUCOV[i1d], OPTMASK, callfun );
      if (sigInt == 0.) { sigInt = 1.0e-12;}
      ptr_MUCOVADD[i1d] = sigInt*fabs(sigInt); // preserve the sign 
    }
  }  // end i1d loop
  
 
  // -------------------------
  // print errBias info in z bins

  int LPRINT = 1 ;

  if ( LPRINT ) {

    double zlo, zhi ;  
    //    printf("NUMBER OF MASS BINS %d\n",NBINm);
    printf("\n");
    if ( DO_COVADD ){
      printf("                            "
	     "RMS(muDif)/RMS(Pull)/SIGINT/NSIM for \n");
    }
    else {
      printf("                            "
	     "RMS(muDif)/RMS(Pull)/NSIM for \n");
    }

    printf("  ia,ib,ig,im  z-range :   "
	   "    ic=0               ic=1                ic=2 \n");
    
    printf("  -------------------------------------------------"
	   "------------------------\n");
    fflush(stdout);
    
    for(ia=0; ia< NBINa; ia++ ) {
      for(ib=0; ib < NBINb; ib++ ) {      
	for(ig=0; ig < NBINg; ig++ ) {      
	  for(iz=0; iz < NBINz; iz++ ) {	
	    zlo = CELL_MUCOVSCALE->BININFO_z.lo[iz];
	    zhi = CELL_MUCOVSCALE->BININFO_z.hi[iz];
	    
	    for(im=0; im < NBINm; im++ ) {	
	      printf("  %d,%d,%d,%d  %.2f-%.2f : ",  ia,ib,ig,im, zlo, zhi );
	      for(ic=0; ic<3; ic++ ) {  
		i1d=CELL_MUCOVSCALE->MAPCELL[ia][ib][ig][iz][im][0][ic];
		N     = CELL_MUCOVSCALE->NperCell[i1d] ;
		muCOVscale = (double)ptr_MUCOVSCALE[i1d] ;
		STD        = sqrt ( SQMUSTD[i1d] );       
		if ( DO_COVADD ) {
		  double covint = (double)ptr_MUCOVADD[i1d];
		  double sigint = (covint/fabs(covint))*sqrt(fabs(covint));
		  printf("%6.3f/%5.3f/%5.3f/%5d ", 
			 STD, sqrt(muCOVscale), sigint, N );	  
		}
		else {
		  printf("%6.3f/%5.3f/%5d ",
			 STD, sqrt(muCOVscale), N );	  
		}
	      } // ic
	      printf("\n");     fflush(stdout);
	    } // im
	    printf("\n");     fflush(stdout);
	  } // iz
	  printf("\n");     fflush(stdout);
	} // ig
      } // ib
    } // ia
    
  } // end LPRINT

  // - - - - - - - - - - - - - - - - - - - 

  if ( debug_mucovscale > 0 ) {   
    write_debug_mucovcorr(IDSAMPLE, muDif_list, muErr_list);

    // xxxxxxxx legacy file Mar 3 2022 RK xxxxxxx
    char outfile[200], line[200], *name; 
    sprintf(outfile,"%s_IDSAMPLE%d_LEGACY.dat", INPUTS.PREFIX, IDSAMPLE); 
    printf("DEBUG: Create diagnostic file %s\n", outfile);
    FILE *fp = fopen(outfile,"wt");
    fprintf(fp,"VARNAMES: CID BIN Ncell "
	    "zMEAN cMEAN mMEAN "
	    "MUCOVSCALE_STD MUCOVSCALE_MAD "
	    "RMS_MUDIF MUDIF "
	    "MUERR MUERR_RAW\n");

    for(isp=0; isp < NBIASCOR_CUTS; isp++ ) {
      ievt = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;
      name = INFO_BIASCOR.TABLEVAR.name[ievt];

      // check if there is valid biasCor for this event
      J1D = J1D_biasCor(ievt,fnam);
      NperCell = CELL_BIASCOR->NperCell[J1D];
      i1d      = INFO_BIASCOR.TABLEVAR.IMUCOV[ievt] ;
      if ( NperCell < BIASCOR_MIN_PER_CELL )  { continue ; }
      if ( i1d < 0 )                          { continue ; }
      if ( !CELL_MUCOVSCALE->USE[i1d] )       { continue; }

      STD = sqrt ( SQMUSTD[i1d] );
      sprintf(line,"SN: "
	      "%8s %d %d "         // name bin Ncell
	      "%.3f %.3f %.3f "    // zMEAN cMEAN mMEAN
	      "%.3f %.3f "         // MUCOVSCALE_RMS MUCOVSCALE_MAD
	      "%.3f %.3f "         // RMS_MUDIF  MUDIF
	      "%.3f %.3f "         // MUERR MUERR_RAW
      	      ,name, i1d, CELL_MUCOVSCALE->NperCell[i1d]
	      ,CELL_MUCOVSCALE->AVG_z[i1d]
	      ,CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d]
              ,CELL_MUCOVSCALE->AVG_m[i1d]
	      ,SIG_PULL_STD[i1d], SIG_PULL_MAD[i1d]
	      ,STD, muDif_list[ievt]
	      ,muErr_list[ievt], muErr_raw_list[ievt]
      	      );
      fprintf(fp, "%s\n", line);
    }  // end isp loop over sparse events

    fclose(fp);
    // xxxxxxxxxxxxxx 

    free(muErr_list); free(muErr_raw_list); free(muDif_list);
  }  // end debug_mucovscale

  // - - - - - 
  print_debug_malloc(-1*debug_malloc,fnam);
  free(SUM_MUERR);   free(SUM_SQMUERR);
  free(SUM_MUDIF);   free(SUM_SQMUDIF) ;
  free(SQMUERR);     free(SQMUSTD);
  free(SUM_PULL);    free(SUM_SQPULL);
 
  if ( DO_MAD ) {
    for (i1d=0; i1d<NCELL; i1d++) { free(CELL_MUCOVSCALE->ABSPULL[i1d]); }
    free(CELL_MUCOVSCALE->ABSPULL);
  }
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
  // Aug 22 2019: include dmHost term based on initial hostPar values.
  // Apr 02 2020: for BIASCOR_MU option, a,b = user input p1,p2  

  double z, a, b, M0, mB, x1, c, logmass, dmHost, hostPar[10];
  double zTrue, muFit, muTrue, muz, muDif ;
  double dlz, dlzTrue, dmu  ;
  char fnam[] = "muresid_biasCor" ;

  // ----------------- BEGIN ----------------

  bool DOBIAS_MU = ( INPUTS.opt_biasCor & MASK_BIASCOR_MU     ) ;

  M0     = INPUTS.nommag0 ;

  if ( DOBIAS_MU ) {
    // sim_alpha[beta] may not exist or make sense, so set
    // a,b to user-input values
    //   fragile alert !!! .xyz
    a  = INPUTS.parval[IPAR_ALPHA0];
    b  = INPUTS.parval[IPAR_BETA0];
  }  
  else{
    a  = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt] ;
    b  = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt] ;
  }

  //  g        = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt] ;
  z        = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt] ;
  logmass  = (double)INFO_BIASCOR.TABLEVAR.host_logmass[ievt];
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

  // need true MU at observed redshift z for biasCor COSPAR. 
  // We don't have access to the biasCor COSPAR, so add 
  // mu-shift (dmu) to SIM_MU(biasCor) where
  //   dmu = mu(z,COSPAR_BBC) - mu(zTrue,COSPAR_BBC) 
  //
  // and COSPAR_BBC are the BBC-input cosmology params

  // get d_l for measured redshift 
  dlz      = cosmodl_forFit(z, z, INPUTS.COSPAR); 

  // get d_l for true redshift
  dlzTrue  = cosmodl_forFit(zTrue, zTrue, INPUTS.COSPAR); 

  dmu    = 5.0*log10(dlz/dlzTrue) ;
  muz    = muTrue + dmu ;  // mu at measured z and biasCor COSPAR
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
  //   maskCov & 4 ==> include zMUERR from VPECERR  (Jun 8 2021)
  //
  // Output:
  //   istat_cov is returned non-zero if the cov matrix was tweaked.
  //   Function returns square of error on distance modulus.
  //
  //  Jun 8 2021: check USEMASK_BIASCOR_ZMUERR option
  //

  bool DOCOVFIT     = ( (maskCov & USEMASK_BIASCOR_COVFIT)>0 ) ;
  bool DOCOVINT     = ( (maskCov & USEMASK_BIASCOR_COVINT)>0 ) ;
  bool DOZMUERR     = ( (maskCov & USEMASK_BIASCOR_ZMUERR)>0 ) ;

  int IDSAMPLE, OPTMASK;
  double zhd, zmuerr=0.0, a, b, gDM;
  double COVTOT[NLCPAR][NLCPAR] ;
  double COVINT[NLCPAR][NLCPAR] ;
  double COVTMP, muErrsq ;
  int  j0, j1 ;
  char *name ;
  char fnam[] = "muerrsq_biasCor" ;
  
  // ------------- BEGIN -----------

  IDSAMPLE = (int)INFO_BIASCOR.TABLEVAR.IDSAMPLE[ievt] ;
  zhd      = (double)INFO_BIASCOR.TABLEVAR.zhd[ievt] ;
  a        = (double)INFO_BIASCOR.TABLEVAR.SIM_ALPHA[ievt] ;
  b        = (double)INFO_BIASCOR.TABLEVAR.SIM_BETA[ievt] ;
  gDM      = (double)INFO_BIASCOR.TABLEVAR.SIM_GAMMADM[ievt] ;
  //  ia       = (int)INFO_BIASCOR.IA[ievt] ;
  //  ib       = (int)INFO_BIASCOR.IB[ievt] ;
  //  ig       = (int)INFO_BIASCOR.IG[ievt] ;
  name     = INFO_BIASCOR.TABLEVAR.name[ievt];
  //  iz       = IBINFUN(zhd,  &INPUTS.BININFO_z, 0, "" );
  if ( DOCOVINT ) 
    { get_COVINT_biasCor(IDSAMPLE,zhd,a,b,gDM, COVINT); } // return COVINT
   
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
  if ( DOZMUERR )  
    { zmuerr = (double)INFO_BIASCOR.TABLEVAR.zmuerr[ievt] ; }

  muErrsq = fcn_muerrsq(name, a, b, gDM, COVTOT, zhd, zmuerr, 0 );

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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  int  IROW=0 ;
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

#ifdef TEXTFILE_NVAR
  int NVAR=13 ;
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
  
  fprintf(FP_STDOUT,"\n");
  fprintf(FP_STDOUT,"# ----------------------------------------- \n");
  fprintf(FP_STDOUT,"# %s: \n", fnam);
  fprintf(FP_STDOUT,
	 "#   idsample=%d  iz=%d(%.3f-%.3f)  ia,ib,ig=%d,%d,%d   NEVT=%d\n", 
	 idsample, iz,zlo,zhi, ia,ib,ig, NEVT);


  for(k0=0; k0 < NLCPAR; k0++ ) {   
    covDiag = (*COV_LOCAL).VAL[k0][k0] ;     
    sigma0  = sqrt(covDiag);
    parName = BIASCOR_NAME_LCFIT[k0] ;

    fprintf(FP_STDOUT, "\t sigma(%2s) = %7.4f | ", parName, sigma0 );
    for(k1=0; k1 < NLCPAR; k1++ ) {
      // print reduced COV here ...
      covDiag = (*COV_LOCAL).VAL[k1][k1] ;
      cov     = (*COV_LOCAL).VAL[k0][k1] ;
      sigma1  = sqrt(covDiag);
      rho     = cov/(sigma0*sigma1);
      fprintf(FP_STDOUT," %8.4f ", rho);
    }
    fprintf(FP_STDOUT,"\n");
  }
  fflush(FP_STDOUT);

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
  double sigInt, COV ;
  char fnam[] = "init_COVINT_biasCor" ;

  //  int DOTEST_SIGINT_ONLY = 1 ; // traditional sigma_int model (debug only)

  // ------------- BEGIN ---------------

  sprintf(BANNER,"%s:", fnam);
  fprint_banner(FP_STDOUT,BANNER); 

  fprintf(FP_STDOUT,
	 "\t Compute Intrinsic Matrix (COVINT) from BiasCor Sample.\n");
  fprintf(FP_STDOUT,
	 "\t COVINT -> bins of IDSAMPLE, Redshift, Alpha, Beta, GammaDM\n");
  fflush(FP_STDOUT);

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

    // use true or meaured redshift ??? [7.01.2020]
    iz       = (int)INFO_BIASCOR.IZ[ievt] ; // true zcmb index
    //    z        = (int)INFO_BIASCOR.TABLEVAR.SIM_ZCMB[ievt] ;

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


  fprintf(FP_STDOUT, "\t %d of %d BiasCor events have IDEAL fit params. \n",
	 NBIASCOR_IDEAL, NBIASCOR_CUTS);
  fflush(FP_STDOUT) ;

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
  // July 25 2021: refactor to use sigint_muresid_list utility.

  int  DO_SIGINT_SAMPLE = ( INPUTS.opt_biasCor & MASK_BIASCOR_SIGINT_SAMPLE ) ;
  int  DOCOR_1D5DCUT    = ( INPUTS.opt_biasCor & MASK_BIASCOR_1D5DCUT );
  int  debug_malloc     = INPUTS.debug_malloc ;
  int  MINEVT_SIGINT_COMPUTE = 50; // abort if fewer events in ia,ib,ig bin

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  int  NROW_TOT, NROW_malloc, istat_cov, NCOVFIX, MEMD, cutmask ;
  int  i, ia, ib, ig ;
  int  LDMP = 0 ;

  double muErrsq, muErr, muDif, SNRMAX, sigInt ;
  int    NSNRCUT, NTMP, NBINa, NBINb, NBINg, USEMASK, DOPRINT ;
  int        *ptr_CUTMASK ;
  short int  *ptr_IDSAMPLE ;
  float      *ptr_SNRMAX;

  int      NUSE[MXa][MXb][MXg];
  double  *MUDIF[MXa][MXb][MXg];
  double  *MUERRSQ[MXa][MXb][MXg];

  double SIGINT_AVG, SIGINT_ABGRID[MXa][MXb][MXg] ; 

  char *NAME, callFun[60];
  char fnam[]   = "init_sigInt_biasCor_SNRCUT" ;

  // ------------------- BEGIN -------------------

  if ( DOCOR_1D5DCUT ) { return ; }

  fprintf(FP_STDOUT,"\n");

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
  if ( sigInt >= 0.0 || DO_COVADD ) {
    if ( DO_COVADD ) { // Aug 2 2021 Dillon
      sigInt = 0.0 ;
      fprintf(FP_STDOUT,
	      " sigInt -> 0: for COV_ADD option IDSAMPLE=%d\n",IDSAMPLE);
    }
    else {
      fprintf(FP_STDOUT,
	      " sigInt -> %.3f from user input sigmb_biascor key IDSAMPLE=%d\n", 
	      sigInt, IDSAMPLE);
    }
    fflush(FP_STDOUT);
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
    fprintf(FP_STDOUT," %s for IDSAMPLE=%d (%s): \n",
	   fnam, IDSAMPLE,  SAMPLE_BIASCOR[IDSAMPLE].NAME );
  }
  else if ( IDSAMPLE == 0 ) {
    fprintf(FP_STDOUT, " %s for all IDSAMPLEs combined: \n", fnam);
  }
  else {
    // do nothing
  } 

  fflush(FP_STDOUT) ;

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
  print_debug_malloc(+1*debug_malloc,fnam);
  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      for(ig=0; ig < NBINg; ig++ ) {
	MUDIF[ia][ib][ig]    = (double*) malloc ( MEMD );
	MUERRSQ[ia][ib][ig]  = (double*) malloc ( MEMD ) ;
	NUSE[ia][ib][ig]     = NCOVFIX = 0 ;
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
    USEMASK = USEMASK_BIASCOR_COVFIT + USEMASK_BIASCOR_ZMUERR;
    muErrsq = muerrsq_biasCor(i, USEMASK, &istat_cov) ;

    if ( muErrsq < 0.0 ) {
      sprintf(c1err,"Invalid muErrsq[%d] = %f < 0  SNID=%s",
	      i, muErrsq, NAME );
      sprintf(c2err,"muDif=%le ", muDif);
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
    }

    if ( isnan(muErrsq) || isnan(muDif) ) {
      sprintf(c1err,"isnan trap for SNID = %s (irow=%d)", NAME, i);
      sprintf(c2err,"muErrsq=%le, muDif=%le ", muErrsq, muDif);
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
    }


    if ( istat_cov < 0 ) { NCOVFIX++ ; }    

    ia = (int)INFO_BIASCOR.IA[i];
    ib = (int)INFO_BIASCOR.IB[i];
    ig = (int)INFO_BIASCOR.IG[i]; 

    NTMP = NUSE[ia][ib][ig];
    MUDIF[ia][ib][ig][NTMP]   = muDif ;
    MUERRSQ[ia][ib][ig][NTMP] = muErrsq ;
    NUSE[ia][ib][ig]++ ;

    if ( NUSE[ia][ib][ig] >= NROW_malloc ) {
      sprintf(c1err,"NUSE[ia,ib,ig=%d,%d,%d] = %d exceeds malloc size",
	      ia, ib, ig, NUSE[ia][ib][ig] );
      sprintf(c2err,"NROW_malloc = %d", NROW_malloc);
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
    }
  
  } // end loop over biasCor sample


  // -------------------------------------------------
  DOPRINT = ( IDSAMPLE==0 || DO_SIGINT_SAMPLE ) ;
  SIGINT_AVG = 0.0 ;

  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {    
      for(ig=0; ig < NBINg; ig++ ) {
	
	NTMP     = NUSE[ia][ib][ig] ;
	if ( NTMP < MINEVT_SIGINT_COMPUTE ) {
	  sprintf(c1err,"NUSE[ia,ib,ig=%d,%d,%d] = %d < %d",
		  ia, ib, ig, NTMP, MINEVT_SIGINT_COMPUTE );
	  sprintf(c2err,"Check biasCor events with SNR> %.1f",
		  INPUTS.snrmin_sigint_biasCor ) ;
	  errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);     
	}
	
	sprintf(callFun,"%s(ia,ib,ig=%d,%d,%d)", fnam, ia, ib, ig);
	sigInt =
	  sigint_muresid_list(NUSE[ia][ib][ig], MUDIF[ia][ib][ig],
			      MUERRSQ[ia][ib][ig], 0, callFun );
	
	// load SIGINT value
	SIGINT_ABGRID[ia][ib][ig] = sigInt ;
	SIGINT_AVG               += sigInt ;

	if ( DOPRINT ) {
	  fprintf(FP_STDOUT,"    sigInt[ia,ib,ig=%d,%d,%d] = %.4f "
		  "(%d evt with IDSAMPLE=%d & SNR>%.0f) \n", 
		  ia, ib, ig, sigInt, NUSE[ia][ib][ig], 
		  IDSAMPLE, INPUTS.snrmin_sigint_biasCor );
	  fflush(FP_STDOUT);
	}
      
      } // ig
    } // ib
  } // ia

 
  SIGINT_AVG /= (double)( NBINa*NBINb*NBINg) ;

  
  // free memory
  print_debug_malloc(-1*debug_malloc,fnam);
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

  int debug_malloc = INPUTS.debug_malloc ;
  short int *ptr_IDSAMPLE, idsample ;
  int       *ptr_CUTMASK  ;
  int NROW, irow, isp, cutmask, MEMI, NBIASCOR ;
  char fnam[] = "makeSparseList_biasCor" ;

  // ---------- BEGIN ------------

  NROW         = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptr_CUTMASK  = INFO_BIASCOR.TABLEVAR.CUTMASK ; // set pointer
  ptr_IDSAMPLE = INFO_BIASCOR.TABLEVAR.IDSAMPLE ;   

  sprintf(BANNER," %s: make sparse row-list for each biasCor sample.",fnam);
  fprint_banner(FP_STDOUT,BANNER);    

  // allocate memory for sparse irow list
  print_debug_malloc(+1*debug_malloc,fnam);
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
    SAMPLE_BIASCOR[idsample].IROW_CUTS[isp] = irow;
    SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS++ ;
  }

  // - - - - -
  for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) {
    printf("\t NBIASCOR_CUTS = %7d for IDSAMPLE=%d \n",
           SAMPLE_BIASCOR[idsample].NBIASCOR_CUTS, idsample);
  }
  fflush(stdout);


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
  int debug_malloc = INPUTS.debug_malloc ;

  int irow, isp, J1D, ipar,  *NperCell ;
  double WGT, z, m, fitPar[NLCPAR];
  double *SUM_WGT_5D, *SUM_z_5D,  *SUM_m_5D, *SUM_FITPAR_5D[NLCPAR];
  float *ptr_z, *ptr_m, *ptr_fitPar[NLCPAR];
  char fnam[] = "makeMap_binavg_biasCor" ;

  // --------------- BEGIN ---------------

  ptr_z  = INFO_BIASCOR.TABLEVAR.zhd ;
  ptr_m  = INFO_BIASCOR.TABLEVAR.host_logmass ; 

  for(ipar=0; ipar < NLCPAR; ipar++ ) 
    { ptr_fitPar[ipar]  = INFO_BIASCOR.TABLEVAR.fitpar[ipar] ; }  

  print_debug_malloc(+1*debug_malloc, fnam);
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


  // loop over biasCor sample ...
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
  print_debug_malloc(-1*debug_malloc, fnam);
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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    zM0   = zmu_solve(muAvg,cosPar); 
    FITRESULT.zM0[iz]      = zM0 ;
    FITRESULT.MUREF_M0[iz] = muAvg;

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
    dl     = cosmodl_forFit(z,z,cosPar);
    mutmp  = 5.0*log10(dl) + 25.0 ;
    dmu    = mutmp - mu ;
    DMU    = fabs(dmu);
    z     *= (1-dmu/10.0) ;

    NITER++ ;
    if ( NITER > 500 ) {
      sprintf(c1err,"Could not solve for z after NITER=%d", NITER);
      sprintf(c2err,"mu=%f  dmu=%f  ztmp=%f", mu, dmu, z);
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  // Sep 14 2021: check istat for get_muCOVcor

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD  ) > 0;

  BIASCORLIST_DEF BIASCORLIST ;
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM;
  int    NBINa, NBINb, NBINg, ia, ib, ig, ipar, istat_bias, idsample, i1d_cen ;
  int    ISTAT = 1;
  double z, m, muCOVscale, muCOVadd ;
  char   *name ;
  char   fnam[] = "storeDataBias" ;

  // ------------- BEGIN -------------

  BININFO_SIM_ALPHA   = &INFO_BIASCOR.BININFO_SIM_ALPHA ;
  BININFO_SIM_BETA    = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;
  NBINa    = (*BININFO_SIM_ALPHA).nbin ;
  NBINb    = (*BININFO_SIM_BETA).nbin ;
  NBINg    = (*BININFO_SIM_GAMMADM).nbin ;

  name        = INFO_DATA.TABLEVAR.name[n];
  idsample    = (int)INFO_DATA.TABLEVAR.IDSAMPLE[n];
  z           = (double)INFO_DATA.TABLEVAR.zhd[n];
  m           = (double)INFO_DATA.TABLEVAR.host_logmass[n];

  if ( DUMPFLAG ) {
    printf("\n");
    printf(" xxx ======================================== \n");
    printf(" xxx =========== %s DUMP for CID=%s ================ \n", 
	   fnam, name );
    printf(" xxx ======================================== \n");
    fflush(stdout);
  }

  if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { return(0); }

  BIASCORLIST.idsample      = idsample;
  BIASCORLIST.z             = z ;
  BIASCORLIST.host_logmass  = m ;

  for(ipar=0; ipar<NLCPAR; ipar++ ) 
    { BIASCORLIST.FITPAR[ipar] = INFO_DATA.TABLEVAR.fitpar[ipar][n]; }

  for(ia = 0; ia < NBINa; ia++ ) {
    for(ib = 0; ib < NBINb; ib++ ) {
      for(ig = 0; ig < NBINg; ig++ ) {

	BIASCORLIST.alpha    = (*BININFO_SIM_ALPHA).avg[ia];
	BIASCORLIST.beta     = (*BININFO_SIM_BETA).avg[ib];
	BIASCORLIST.gammadm  = (*BININFO_SIM_GAMMADM).avg[ig];

	istat_bias = 
	  get_fitParBias(name, &BIASCORLIST, DUMPFLAG, fnam,            // in
			 &INFO_DATA.FITPARBIAS_ALPHABETA[n][ia][ib][ig]);//out
	if ( istat_bias <= 0 ) { ISTAT = 0 ; }
	if ( DUMPFLAG ) {
	  printf(" xxx %s: a=%.2f b=%.2f gDM=%5.2f (ia,ib,ig=%d,%d,%d) "
		 "istat_bias=%d \n",
		 fnam, BIASCORLIST.alpha,BIASCORLIST.beta,BIASCORLIST.gammadm,
		 ia,ib,ig, istat_bias);
	  fflush(stdout);
	}
	

	istat_bias = 
	  get_muCOVcorr(name, &BIASCORLIST, DUMPFLAG,    // in
			&muCOVscale, &muCOVadd, &i1d_cen );        // out
	if ( istat_bias <= 0 ) { ISTAT = 0 ; }        //  Sep 2021

	INFO_DATA.MUCOVSCALE_ALPHABETA[n][ia][ib][ig] = muCOVscale ;
	INFO_DATA.I1D_MUCOVSCALE[n][ia][ib][ig]       = i1d_cen ;

	if ( DO_COVADD ) {
	  INFO_DATA.MUCOVADD_ALPHABETA[n][ia][ib][ig] = muCOVadd ;
	}

	if ( DUMPFLAG ) {
	  printf(" xxx %s: a=%.2f b=%.2f gDM=%.2f (ia,ib,ig=%d,%d,%d) "
		 "istat_muCOVscale=%d \n",
		 fnam, BIASCORLIST.alpha,BIASCORLIST.beta,BIASCORLIST.gammadm,
		 ia, ib, ig, istat_bias);
	  fflush(stdout);
      }
	
	if ( istat_bias == 0 ) { ISTAT = 0 ; }

      } // end ig
    } // end ib
  }  // end ia



  if ( INPUTS.debug_mucovscale > 0 && n < 40 )   { dump_muCOVcorr(n); }


  return(ISTAT);

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
  // Nov 29 2021: return if idsample < 0

  int  NBINa   = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin ;
  int  NBINb   = INFO_BIASCOR.BININFO_SIM_BETA.nbin ;
  int  NBINg   = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin ;
  char *name   = INFO_CCPRIOR.TABLEVAR.name[n];

  BIASCORLIST_DEF BIASCORLIST ;

  int    DUMPFLAG = 0; // (strcmp(name,"184000") == 0) ;
  int    ia, ib, ig, istat_bias, idsample ;
  char   fnam[] = "storeBias_CCprior" ;

  // ------------- BEGIN -------------

  // load local BIASCORLIST struct
  BIASCORLIST.z = 
    INFO_CCPRIOR.TABLEVAR.zhd[n];

  BIASCORLIST.host_logmass = 
    INFO_CCPRIOR.TABLEVAR.host_logmass[n];

  BIASCORLIST.FITPAR[INDEX_mB] = 
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_mB][n];

  BIASCORLIST.FITPAR[INDEX_x1] =
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_x1][n];

  BIASCORLIST.FITPAR[INDEX_c] =
    INFO_CCPRIOR.TABLEVAR.fitpar[INDEX_c][n];

  BIASCORLIST.idsample =
    INFO_CCPRIOR.TABLEVAR.IDSAMPLE[n];

  idsample = BIASCORLIST.idsample ;
  if ( SAMPLE_BIASCOR[idsample].DOFLAG_SELECT == 0 ) { return(0); }

  if ( idsample < 0 ) {
    sprintf(c1err,"Undefined IDSAMPLE for CID=%s", name);
    sprintf(c2err," ");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

 
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
	get_fitParBias(name, &BIASCORLIST, DUMPFLAG, fnam, 
		       &INFO_CCPRIOR.FITPARBIAS_ALPHABETA[n][ia][ib][ig] );

      /* 
      if ( DUMPFLAG ) {
	printf(" xxx %s: a,b=%.2f,%.2f  BIAS(mB,x1,c)=%.3f,%.3f,%.3f \n",
	       fnam, BIASCORLIST.alpha, BIASCORLIST.beta,
	       FITPARBIAS_TMP.VAL[0], FITPARBIAS_TMP.VAL[1],
	       FITPARBIAS_TMP.VAL[2] );    fflush(stdout);
      }
      */

      if ( istat_bias < 0 ) { return 0 ; }

      } // end ig
    } // end ib
  }  // end ia
  
  return(1);

} // end storeBias_CCPrior


// ======================================================
int get_fitParBias(char *cid, 
		   BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG, char *callFun,
		   FITPARBIAS_DEF  *FITPARBIAS) {

  // Created May 2016
  // For input BIASCORLIST = { mB,x1,c, z, alpha, beta, gammadm },
  // return FITPARBIAS[VAL,ERR,RMS] for  mB, x1, c.
  // CID is the SN name used only for error message.
  //
  // Function returns 1 if bias is found; 
  // returns negative error code if no bias found.
  //
  // Apr 18 2017: enhance dump output.
  //
  // Nov 18 2019: check option to interpolate biasCor vs. logMass.
  //
  // Sep 29 2020:
  //   + pass callFun for error message
  //   + requier central cell is used; see USE_CENTER_CELL
  //
  // Dec 21 2020: return negative number on error to flag error type
  //              Still return 1 for success.
  // 
  // Jun 2022: improve dump output
  //
  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z   = BIASCORLIST->z ;
  double m   = BIASCORLIST->host_logmass ;
  double mB  = BIASCORLIST->FITPAR[INDEX_mB];
  double x1  = BIASCORLIST->FITPAR[INDEX_x1];
  double c   = BIASCORLIST->FITPAR[INDEX_c];
  double a   = BIASCORLIST->alpha ;
  double b   = BIASCORLIST->beta ;
  double gDM = BIASCORLIST->gammadm ;
  int IDSAMPLE = BIASCORLIST->idsample ;
  int ID = IDSAMPLE;

  bool DO_BIASCOR_MU  = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  bool REQUIRE_VALID_BIASCOR = (INPUTS.opt_biasCor & MASK_BIASCOR_noCUT)==0;

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
  int BAD_BINSIZE=0 ;
  bool USE_CENTER_CELL;

  double WGT, SUM_WGT, BINSIZE;
  double AVG_z, AVG_m, AVG_x1, AVG_c, avg_z, avg_m, avg_x1, avg_c ;
  double SUM_VAL[NLCPAR+1], SUM_ERR[NLCPAR+1];
  double SUM_SQERRINV[NLCPAR], SUM_SQRMS[NLCPAR+1] ;
  double VAL, ERR, RMS, dif, Dc, Dz, Dm, Dx1, VAL_DEFAULT ;

  double DEBUG_LIST_DIF[4][50];
  int    DEBUG_LIST_INDX[4][50], DEBUG_LIST_NPERCELL[50];
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

  // init output
  if (REQUIRE_VALID_BIASCOR) {VAL_DEFAULT=666.0;}  else{VAL_DEFAULT=0.0;}
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ )
    { FITPARBIAS->VAL[ipar]  = VAL_DEFAULT ; }

  // strip off local indices
  ia  = IBINFUN(a,   BININFO_SIM_ALPHA,    0, fnam );
  ib  = IBINFUN(b,   BININFO_SIM_BETA,     0, fnam );
  ig  = IBINFUN(gDM, BININFO_SIM_GAMMADM,  0, fnam );
  IZ  = IBINFUN(z,  &CELLINFO_BIASCOR[ID].BININFO_z,         0, fnam );
  IM  = IBINFUN(m,  &CELLINFO_BIASCOR[ID].BININFO_m,         0, fnam );
  IX1 = IBINFUN(x1, &CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_x1],0, fnam);
  IC  = IBINFUN(c,  &CELLINFO_BIASCOR[ID].BININFO_LCFIT[INDEX_c], 0, fnam);

  if ( IZ  < 0 ) { return -1 ; }
  if ( IM  < 0 ) { return -2 ; }
  if ( IX1 < 0 ) { return -3 ; }
  if ( IC  < 0 ) { return -4 ; }

  J1D = CELLINFO_BIASCOR[IDSAMPLE].MAPCELL[ia][ib][ig][IZ][IM][IX1][IC] ;

  // reset counters and weight for this a,b cell
  SUM_WGT = 0.0 ;
  NCELL_INTERP_TOT = NCELL_INTERP_USE  = NSUM_Cell = 0 ;
  
  for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX; ipar++ ) {	    
    SUM_VAL[ipar]      = SUM_ERR[ipar]   = 0.0 ;
    SUM_SQERRINV[ipar] = SUM_SQRMS[ipar] = 0.0  ;
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
    printf(" xxx %s DUMP for CID=%s  BIASCOR(alpha,beta) = %.2f, %.2f\n", 
	   fnam, cid, a, b );
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
  if ( BINSIZE_z  == 0.0 || BINSIZE_z  > 9000. ) { BAD_BINSIZE = -5; }
  if ( BINSIZE_m  == 0.0 || BINSIZE_m  > 9000. ) { BAD_BINSIZE = -6; }
  if ( BINSIZE_x1 == 0.0 || BINSIZE_x1 > 9000. ) { BAD_BINSIZE = -7; }
  if ( BINSIZE_c  == 0.0 || BINSIZE_c  > 9000. ) { BAD_BINSIZE = -8; }

  if ( BAD_BINSIZE < 0 ) {
    if ( LDMP ) { printf(" xxx\t BAD BINSIZE --> FAIL\n"); fflush(stdout);  }
    return BAD_BINSIZE ;
  }

  
  // ----------------------------------------------
  // -------- start 4D loop over cells -------------
  
  USE_CENTER_CELL = false;

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
	    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
	  }
	
	  // require something in a cell to use the bias estimate
	  NperCell = CELLINFO_BIASCOR[ID].NperCell[j1d] ;
	  if ( NperCell < BIASCOR_MIN_PER_CELL  ) { continue; }
	  
	  if ( iz==IZ && im==IM && ix1==IX1 && ic==IC ) 
	    { USE_CENTER_CELL = true ; }

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
	  DEBUG_LIST_NPERCELL[NCELL_INTERP_USE] = NperCell ;

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
  
  
  int  NCELL_INTERP_MIN = 3;
  bool PASSCUT_NCELL_INTERP = ( NCELL_INTERP_USE >= NCELL_INTERP_MIN );

  if ( LDMP ) {
    char passfail[] = "PASS";
    if ( !PASSCUT_NCELL_INTERP ) { sprintf(passfail,"FAIL"); }
    
    printf("\n") ;
    printf(" xxx\t SUM_WGT = %le  \n", 	SUM_WGT );
    printf(" xxx\t NCELL_INTERP_USE=%d of %d --> %s CUT>=%d\n", 
	   SUM_WGT, NCELL_INTERP_USE, NCELL_INTERP_TOT, 
	   passfail, NCELL_INTERP_MIN );
    fflush(stdout);
  }

   // require enough cells for interpolation (July 2016)
   if ( !PASSCUT_NCELL_INTERP ) { return(-9) ; } 
   if ( !USE_CENTER_CELL      ) { return(-10); } // 9.29.2020

   // require both z-bins to be used.
   int ISKIP = 0 ;
   if ( IZMIN >=  0 ) { ISKIP = (USEZ[IZMIN]==0) ; }
   if ( IZMAX >=  0 ) { ISKIP = (USEZ[IZMAX]==0) ; }

   if ( LDMP ) {
     printf(" xxx IZMIN=%d IZMAX=%d  USEZ=%d,%d  ISKIP=%d\n",
	    IZMIN, IZMAX, USEZ[IZMIN], USEZ[IZMAX], ISKIP );
     fflush(stdout);
   }
   if ( ISKIP ) { return(-11); }

   
   OK = 0 ;
   int NTMP = NSUM_Cell ;
   if ( NTMP >= BIASCOR_MINSUM ) { OK = 1; }
   if ( NTMP >= 2 &&  z > INPUTS.zmax ) { OK = 1; }

   if ( LDMP ) {
     printf(" xxx NSUM_Cell=%d  (cut=%d) OK=%d \n",
	    NSUM_Cell, BIASCOR_MINSUM, OK); fflush(stdout);
   }

   if ( OK == 0 ) { return(-12); }   

  // - - - - - - - - - - 
  // make sure that SUMWGT > 0 
  if ( SUM_WGT <= 1.0E-9 ) {
    int icell ;
    print_preAbort_banner(fnam);
    printf("  Called by function: %s \n", callFun );
    printf("  IZMIN/MAX=%d/%d   IMMIN/MAX=%d,%d  "
	   "IX1MIN/MAX=%d/%d   ICMIN/MAX=%d/%d\n",	   
	   IZMIN,IZMAX,  IMMIN, IMMAX, IX1MIN,IX1MAX,   ICMIN, ICMAX);
    printf("  VALUE(z,m,x1,c)   = %.4f, %.4f, %.4f, %.4f \n",
	   z, m, x1, c );
    printf("  BINSIZE(z,m,x1,c) = %.4f, %.4f, %.4f, %.4f \n",
	   BINSIZE_z, BINSIZE_m, BINSIZE_x1, BINSIZE_c );
    
    printf("\n");
    printf("   cell   iz im ix1 ic    Dz     Dm     Dx1     Dc  "
	   "Ncell  WGT \n");
    
    for(icell=0; icell < NCELL_INTERP_USE; icell++ ) {
      printf("    %3d   %2d %2d %2d %2d   %6.3f %6.3f %6.3f %6.3f  %2d  %.3le\n"
	     ,icell
	     ,DEBUG_LIST_INDX[0][icell]
	     ,DEBUG_LIST_INDX[1][icell]
	     ,DEBUG_LIST_INDX[2][icell]
	     ,DEBUG_LIST_INDX[3][icell]
	     ,DEBUG_LIST_DIF[0][icell]
	     ,DEBUG_LIST_DIF[1][icell]
	     ,DEBUG_LIST_DIF[2][icell]
	     ,DEBUG_LIST_DIF[3][icell]
	     ,DEBUG_LIST_NPERCELL[icell]
	     ,DEBUG_LIST_WGT[icell] );
    }
    sprintf(c1err,"SUM_WGT=%le for  CID=%s", SUM_WGT, cid) ;
    sprintf(c2err,"a=%.3f b=%.2f gDM=%.3f  "
	    "z=%.4f  m=%.2f mB=%.4f x1=%.4f c=%.4f", 
	    a, b, gDM, z, m, mB, x1, c ) ;
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
int get_muCOVcorr(char *cid, 
		   BIASCORLIST_DEF *BIASCORLIST, int DUMPFLAG,
		  double *muCOVscale, double *muCOVadd, int *i1d_cen ) {

  // Created July 1 2016: 
  // Analog of get_fitParBias(), but for scale on muCOV.
  // 
  // inputs:
  //   cid          = name of SN (for error message)
  //   BIASCORLIST  = list of input values
  //
  // Output:
  //   muCovScale   # scale for covariance
  //   muCOVadd     # optional additive term for cov
  //   i1d_cen      #  central i1d bin index
  //
  //   Function returns 1 on success; 0 on failure
  //
  // Jun 11 2021: include logmass dependence.
  // Jun 29 2021: check INPUTS.interp_biascor_logmass
  // Sep 14 2021: check CELL_MUCOVSCALE->USE[j1d] and return 0 for
  //              events that have no valid MUSCALE -> bug fix !
  //
  // Sep 27 2021: return i1d_cen for diagnostic/debug/abort-message.

  int ia, ib, ig, iz, im, ic, IZ, IM, IC, j1d, IMMIN, IMMAX ;

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  // -----------------------------------------
  // strip BIASCORLIST inputs into local variables
  double z  = BIASCORLIST->z ;
  double m  = BIASCORLIST->host_logmass ;
  double c  = BIASCORLIST->FITPAR[INDEX_c];
  double a  = BIASCORLIST->alpha ;
  double b  = BIASCORLIST->beta ;
  double g  = BIASCORLIST->gammadm ;
  int IDSAMPLE = BIASCORLIST->idsample ;  

  // get bin info into local variables
  CELLINFO_DEF *CELL_MUCOVSCALE = &CELLINFO_MUCOVSCALE[IDSAMPLE];
  CELLINFO_DEF *CELL_MUCOVADD   = &CELLINFO_MUCOVADD[IDSAMPLE];
  int  NBINz   = CELL_MUCOVSCALE->BININFO_z.nbin ;
  int  NBINm   = CELL_MUCOVSCALE->BININFO_m.nbin ;
  int  NBINc   = CELL_MUCOVSCALE->BININFO_LCFIT[INDEX_c].nbin ;
 
  double BINSIZE_z = CELL_MUCOVSCALE->BININFO_z.binSize ;
  double BINSIZE_m = CELL_MUCOVSCALE->BININFO_m.binSize ;
  double BINSIZE_c = CELL_MUCOVSCALE->BININFO_LCFIT[INDEX_c].binSize ;
  
  BININFO_DEF *BININFO_SIM_ALPHA, *BININFO_SIM_BETA, *BININFO_SIM_GAMMADM ;
  float *ptr_MUCOVSCALE ;
  float *ptr_MUCOVADD ;

  double dif, Dz, Dm, Dc, WGT, SUM_WGT, SUM_muCOVscale, SUM_muCOVadd; 
  double muCOVscale_biascor, muCOVadd_biascor ;
  double muCOVscale_local, muCOVadd_local  ;
  bool USEBIN_CENTER ;
  char fnam[] = "get_muCOVcorr" ;

  // -------------- BEGIN --------------

  BININFO_SIM_ALPHA    = &INFO_BIASCOR.BININFO_SIM_ALPHA;
  BININFO_SIM_BETA     = &INFO_BIASCOR.BININFO_SIM_BETA ;
  BININFO_SIM_GAMMADM  = &INFO_BIASCOR.BININFO_SIM_GAMMADM ;
  ptr_MUCOVSCALE       = INFO_BIASCOR.MUCOVSCALE[IDSAMPLE] ;
  if ( DO_COVADD ) { ptr_MUCOVADD = INFO_BIASCOR.MUCOVADD[IDSAMPLE]; }

  // init output
  muCOVscale_local = 1.0 ;
  *muCOVscale = muCOVscale_local ;

  muCOVadd_local = 0.0 ;
  *muCOVadd = muCOVadd_local ;

  *i1d_cen = -9;

  if ( (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE)==0   ) { return(1); }

  // check "noBiasCor" option (Aug 20 2016)
  if  ( SAMPLE_BIASCOR[IDSAMPLE].DOFLAG_BIASCOR == 0 ) { return(1); }

  // get indices
  ia  = IBINFUN(a,  BININFO_SIM_ALPHA,   0, "" );
  ib  = IBINFUN(b,  BININFO_SIM_BETA,    0, "" );
  ig  = IBINFUN(g,  BININFO_SIM_GAMMADM, 0, "" );
  IZ  = IBINFUN(z,  &CELL_MUCOVSCALE->BININFO_z,  0, "" );

  // for color & logmass bin, use OPT=2 to allow values outside BININFO 
  // range to be pulled into edge bin.
  IM  = IBINFUN(m, &CELL_MUCOVSCALE->BININFO_m,              2, "");
  IC  = IBINFUN(c, &CELL_MUCOVSCALE->BININFO_LCFIT[INDEX_c], 2, ""); 
  
  SUM_WGT = SUM_muCOVscale = SUM_muCOVadd = 0.0 ;

  if ( INPUTS.interp_biascor_logmass ) 
    { IMMIN = IM-1; IMMAX = IM+1;}
  else
    { IMMIN = IMMAX = IM; }  // do NOT interp logmass dimension

  USEBIN_CENTER = false;
  if ( INPUTS.restore_mucovscale_bug ) { USEBIN_CENTER = true; }

  // - - - - - - 

  for(iz = IZ-1; iz <= IZ+1; iz++ ) {
    if ( iz < 0 )      { continue ; }
    if ( iz >= NBINz ) { continue ; }

    for(im = IMMIN; im <= IMMAX; im++ ) {
      if ( im < 0 )      { continue ; }
      if ( im >= NBINm ) { continue ; }
          
      for(ic = IC-1; ic <= IC+1; ic++ ) {
	if ( ic < 0      ) { continue ; } 
	if ( ic >= NBINc ) { continue ; }  

	j1d = CELL_MUCOVSCALE->MAPCELL[ia][ib][ig][iz][im][0][ic] ;

	if ( !CELL_MUCOVSCALE->USE[j1d] ) { continue; } // 9 2021

	if ( j1d >= CELL_MUCOVSCALE->NCELL || j1d < 0 ) {
	  sprintf(c1err,"Invalid j1d=%d ", j1d);
	  sprintf(c2err,"ia=%d ib=%d iz=%d im=%d ic=%d",
		  ia, ib, iz, im, ic);
	  errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
	}
       
	if ( iz==IZ && ic==IC && im==IM ) 
	  { USEBIN_CENTER = true; *i1d_cen = j1d; }

	// get distance between current data value and bin-center  
	dif = z - CELL_MUCOVSCALE->AVG_z[j1d];
	Dz  = fabs(dif/BINSIZE_z) ;

	dif = m - CELL_MUCOVSCALE->AVG_m[j1d];	
	if ( INPUTS.interp_biascor_logmass ) 
	  { Dm  = fabs(dif/BINSIZE_m) ;  }
	else 
	  { Dm  = 0.0 ; } // Jun 2021	
	
	dif = c - CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][j1d];
	Dc  = fabs(dif/BINSIZE_c);
	
	// only keep closest bins to interpolate
	if ( Dz  > 0.99 ) { continue ; }
	if ( Dc  > 0.99 ) { continue ; }
	if ( Dm  > 0.99 ) { continue ; }
	
	WGT        = (1.0 - Dz) * (1.0 - Dm) * (1.0 - Dc) ;
	SUM_WGT   += WGT ;
	
	muCOVscale_biascor = (double)ptr_MUCOVSCALE[j1d];
	SUM_muCOVscale    += ( WGT * muCOVscale_biascor );
	
	if ( DO_COVADD ) {
	  muCOVadd_biascor = (double)ptr_MUCOVADD[j1d];
	  SUM_muCOVadd += ( WGT * muCOVadd_biascor );
	}

	if ( DUMPFLAG) {
	  printf(" xxx %s: j1d=%d, ic=%d im=%d iz=%d "
		 "muCOVscale/add = %.3f/%.5f WGT=%f\n", 
		 fnam, j1d, ic,im,iz, 
		 muCOVscale_biascor,muCOVadd_biascor, WGT);
	  
	  fflush(stdout);
	}

      } // ic
    }  // im
  } // iz
  
  // - - - - -
  if( SUM_WGT > 0.01 && USEBIN_CENTER ) {
    muCOVscale_local = SUM_muCOVscale / SUM_WGT ;
    if ( DO_COVADD ) {
      muCOVadd_local = SUM_muCOVadd / SUM_WGT ;
    }
  }
  else if ( !INPUTS.restore_mucovscale_bug ) {
    return(0); // 9.2021
  }

  // - - - - - - - -
  if ( DUMPFLAG) {
    printf(" xxx %s: IC=%d IM=%d IZ=%d muCOVscale = %f/%f = %f \n", 
	   fnam,IC,IM,IZ, SUM_muCOVscale, SUM_WGT, muCOVscale_local);
    if ( DO_COVADD ) {
      printf(" xxx %s: IC=%d IM=%d IZ=%d muCOVadd = %f/%f = %f \n", 
	     fnam,IC,IM,IZ, SUM_muCOVadd, SUM_WGT, muCOVadd_local);
    }
    fflush(stdout);
  }

  // load output arg.
  *muCOVscale = muCOVscale_local ;
  *muCOVadd   = muCOVadd_local ;

  return(1);

} // end get_muCOVcorr



// ======================================================
void write_debug_mucovcorr(int IDSAMPLE, double *muDif_list, double *muErr_list) {

  // if debug_mucovscale is set, write two files:
  // Info per I1D mucov bin
  // Info per biasCor event.

  CELLINFO_DEF *CELL_MUCOVSCALE = &CELLINFO_MUCOVSCALE[IDSAMPLE];
  CELLINFO_DEF *CELL_MUCOVADD   = &CELLINFO_MUCOVADD[IDSAMPLE];

  int NCELL  = CELLINFO_BIASCOR[IDSAMPLE].NCELL ;
  FILE *fp;
  int  i1d, NperCell;
  bool USE;
  char outfile[200], line[200], *name; 
  char fnam[] = "write_debug_mucovcorr";

  // --------- BEGIN ----------

  // start with muCOV[add,scale] per bin

  sprintf(outfile,"%s_IDSAMPLE%d_MUCOVCORR.DAT", INPUTS.PREFIX, IDSAMPLE); 
  printf("DEBUG: Create diagnostic file %s\n", outfile);
  fp = fopen(outfile,"wt");
  fprintf(fp,"VARNAMES: ROW BIN_MUCOV Ncell  zMEAN cMEAN mMEAN "
	  "MUCOVSCALE MUCOVADD\n");

  for(i1d=0; i1d < NCELL ; i1d++ ) {

      // check if there is valid biasCor for this event
    NperCell = CELLINFO_MUCOVSCALE[IDSAMPLE].NperCell[i1d];
    USE      = CELLINFO_MUCOVSCALE[IDSAMPLE].USE[i1d] ;

    if ( NperCell < BIASCOR_MIN_PER_CELL )  { continue ; }
    if ( !USE ) { continue; }

    sprintf(line,"SN: "
	    "%4d %4d %4d  "       // ROW bin Ncell
	    "%5.3f %6.3f %6.3f "    // zMEAN cMEAN mMEAN
	    "%7.3f %7.4f"              // MUCOVSCALE MUCOVADD
	    ,i1d, i1d
	    ,CELL_MUCOVSCALE->NperCell[i1d]
	    ,CELL_MUCOVSCALE->AVG_z[i1d]
	    ,CELL_MUCOVSCALE->AVG_LCFIT[INDEX_c][i1d]
	    ,CELL_MUCOVSCALE->AVG_m[i1d]
	    ,INFO_BIASCOR.MUCOVSCALE[IDSAMPLE][i1d]
	    ,INFO_BIASCOR.MUCOVADD[IDSAMPLE][i1d]
	    );

    fprintf(fp, "%s\n", line);
  }  // end i1d loop over cells

  fclose(fp);

  // - - - - - - - - - -  -
  // 2nd file is each biasCor event for this IDSAMPLE

  sprintf(outfile,"%s_IDSAMPLE%d_BIASCOR.DAT", INPUTS.PREFIX, IDSAMPLE); 
  printf("DEBUG: Create diagnostic file %s\n", outfile);
  fp = fopen(outfile,"wt");
  fprintf(fp,"VARNAMES: CID  BIN_MUCOV  zHD c "
	    "MUDIF  MUERR  MUCOVSCALE MUCOVADD\n");

  CELLINFO_DEF *CELL_BIASCOR  = &CELLINFO_BIASCOR[IDSAMPLE];
  int NBIASCOR_CUTS    = SAMPLE_BIASCOR[IDSAMPLE].NBIASCOR_CUTS ;
  int NBIASCOR_ALL     = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  int isp, ievt, J1D;

  for(isp=0; isp < NBIASCOR_CUTS; isp++ ) {
      ievt = SAMPLE_BIASCOR[IDSAMPLE].IROW_CUTS[isp] ;

      // check if there is valid biasCor for this event
      J1D = J1D_biasCor(ievt,fnam);
      NperCell = CELL_BIASCOR->NperCell[J1D];
      if ( NperCell < BIASCOR_MIN_PER_CELL )  { continue ; }

      // check for valid muCov correction
      i1d      = INFO_BIASCOR.TABLEVAR.IMUCOV[ievt] ;
      if ( i1d < 0 )                          { continue ; }
      if ( !CELL_MUCOVSCALE->USE[i1d] )       { continue; }

      name = INFO_BIASCOR.TABLEVAR.name[ievt];
      sprintf(line,"SN: "
	      "%8s %4d  "         // name bin 
	      "%5.3f %6.3f "      // zHD, c
	      "%6.3f %.3f "         // MUDIF MUERR 
	      "%7.3f %7.4f"         // MUCOVSCALE MUCOVADD	      
      	      ,INFO_BIASCOR.TABLEVAR.name[ievt]
	      ,INFO_BIASCOR.TABLEVAR.IMUCOV[ievt]
	      ,INFO_BIASCOR.TABLEVAR.zhd[ievt]
	      ,INFO_BIASCOR.TABLEVAR.fitpar[INDEX_c][ievt]
	      ,muDif_list[ievt]
	      ,muErr_list[ievt]
	      ,INFO_BIASCOR.MUCOVSCALE[IDSAMPLE][i1d]
	      ,INFO_BIASCOR.MUCOVADD[IDSAMPLE][i1d]
      	      );
      fprintf(fp, "%s\n", line);
    }  // end isp loop over sparse events

    fclose(fp);

  return;

} // end write_debug_mucovcorr


// ======================================================
void dump_muCOVcorr(int n) {

  int NBINa = INFO_BIASCOR.BININFO_SIM_ALPHA.nbin;
  int NBINb = INFO_BIASCOR.BININFO_SIM_BETA.nbin;
  int NBINg = INFO_BIASCOR.BININFO_SIM_GAMMADM.nbin;
  char *name = INFO_DATA.TABLEVAR.name[n];

  int ia, ib, ig, i1d;
  double MUCOVSCALE, MUCOVADD;
  char fnam[] = "dump_muCOVcorr" ;

  // -----------BEGIN ------------

  for(ia=0; ia < NBINa; ia++ ) {
    for(ib=0; ib < NBINb; ib++ ) {
      for(ig=0; ig < NBINg; ig++ ) {

        MUCOVSCALE = INFO_DATA.MUCOVSCALE_ALPHABETA[n][ia][ib][ig];
	MUCOVADD   = INFO_DATA.MUCOVADD_ALPHABETA[n][ia][ib][ig];
	i1d        = INFO_DATA.I1D_MUCOVSCALE[n][ia][ib][ig];
	printf(" %s: CID=%s ia,ib,ig=%d,%d,%d  "
	       "MUCOV[i1d,SCALE,ADD]= %d, %.3f, %.4f\n",
	       fnam, name, ia,ib,ig, i1d, MUCOVSCALE, MUCOVADD );	
      }
    }
  }

  fflush(stdout);

  return;

} // end dump_muCOVcorr

// ======================================================
void setup_CELLINFO_biasCor(int IDSAMPLE) {

  // Created Aug 23 2016
  // Setup 5D biasCor cells for input IDSAMPLE

  int NSAMPLE   = NSAMPLE_BIASCOR ;
  int MEMCELL   = NSAMPLE * sizeof(CELLINFO_DEF);
  int INDX;
  int debug_malloc = INPUTS.debug_malloc ;
  char fnam[] = "setup_CELLINFO_biasCor" ;

  // ------------ BEGIN -----------

  // on first call, allocate CELLINFO structures
  if ( IDSAMPLE == 0 ) {

    fprintf(FP_STDOUT, "\n# ============== %s ================= \n", fnam);

    print_debug_malloc(+1*debug_malloc,fnam);
    CELLINFO_BIASCOR    = (CELLINFO_DEF*) malloc ( MEMCELL );
    CELLINFO_MUCOVSCALE = (CELLINFO_DEF*) malloc ( MEMCELL );
    CELLINFO_MUCOVADD   = (CELLINFO_DEF*) malloc ( MEMCELL );

    // setup bining for SIMalpha,beta; note storage in different structure
    // since alpha,beta,gammaDM binning is fixed for all IDSAMPLEs
    fprintf(FP_STDOUT, "\n\t\t Global Alpha,Beta Binning \n"); 
    fflush(FP_STDOUT);

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

  fprintf(FP_STDOUT, "\n\t Setup biasCor bins for SAMPLE = %s \n",
	 SAMPLE_BIASCOR[IDSAMPLE].NAME ); fflush(FP_STDOUT);

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
  fflush(FP_STDOUT);

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
    get_BININFO_biasCor_abg("SIM_alpha", 
			    &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIM_alpha");
  }
  else if ( ipar_LCFIT == 100*INDEX_c ) {
    // get info for SIM_beta binning
    get_BININFO_biasCor_abg("SIM_beta", 
			    &VAL_MIN, &VAL_MAX, &VAL_BIN );
    sprintf(NAME,"SIM_beta");
  }
  else if ( ipar_LCFIT == 300 ) {
    // get info for SIM_gammaDM binning (Jul 2019)
    get_BININFO_biasCor_abg("SIM_gammaDM", 
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, ""); 
  }


  // allow 1 bin for beta and gamma; but not for other params
  bool ONEBIN = (VAL_MAX == VAL_MIN) ;
  bool OK1BIN = ( ONEBIN && (ipar_LCFIT >= 100) ) ;
  if ( VAL_MAX <= VAL_MIN  && !OK1BIN ) {
    sprintf(c1err,"%s VAL_MAX=%f  <=  VAL_MIN=%f", 
	    NAME, VAL_MAX, VAL_MIN);
    sprintf(c2err,"VAL_BIN = %f \n", VAL_BIN);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

  fprintf(FP_STDOUT, "  %s: define %2d %2s bins from %6.3f to %6.3f \n",
	 fnam, nbin, NAME, BININFO->lo[0], BININFO->hi[nbin-1] );

  if ( nbin > MAXBIN || nbin >= MXpar ) {
    sprintf(c1err,"nbin = %d exceeds MAXBIN=%d or MXpar=%d", 
	    nbin, MAXBIN, MXpar );
    sprintf(c2err,"Check %s binning.", NAME);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);       
  }

  fflush(FP_STDOUT);

  return ;

} // end setup_BININFO_biasCor


// ==================================================
void  get_BININFO_biasCor_abg(char *varName, 
			      double *VAL_MIN, double *VAL_MAX, 
			      double *VAL_BIN ) {

  // Return binning for SIM_alpha or SIM_beta or SIM_gammaDM.
  // This is tricky because there are two ways to define
  // alpha,beta bins:
  //
  // 1) alpha,beta,gammaDM are defined on a grid
  //      --> use native grid from sim file
  // 2) alpha,beta,gammaDM are continuous
  //      --> force 2x2 grid
  //
  // Apr 28 2017: abort if NVAL > NBMAX
  // Jul 16 2019: allow for gammaDM
  // Mar 31 2020: store min/max per IDSURVEY (since IDSAMPLE not yet known)
  // Jul 03 2020: if ipar[IPAR_ABG] == 3, force one bin
  // Feb 09 2022: fix to work for fixed param; e.g. ipar[IPAR_ABG]=0
  //
  double valmin_loc, valmax_loc, valbin_loc;
  double val, val_last, val1st, val2nd ;
  double valmin_sample[MXIDSURVEY], valmax_sample[MXIDSURVEY];
  float  *ptrVal_f = NULL;
  short int *ptrVal_index, *ptr_IDSURVEY ;
  bool   FORCE_ONEBIN = false ; 
  int    irow, unsort, NVAL, NBMAX, ipar, IPAR_ABG=-9 ;
  int    NROW, IDSURVEY, NONIA_INDEX ;
  int debug_malloc = INPUTS.debug_malloc ;
  char fnam[] = "get_BININFO_biasCor_abg" ;

  // ------------- BEGIN --------------

  *VAL_MIN = *VAL_MAX = *VAL_BIN = 0.0 ;

  valmin_loc = +1.0E8 ;
  valmax_loc = -1.0E8 ;
  valbin_loc =  0.0 ;

  NROW         = INFO_BIASCOR.TABLEVAR.NSN_ALL ;
  ptrVal_index = INFO_BIASCOR.TABLEVAR.SIM_NONIA_INDEX;
  ptr_IDSURVEY = INFO_BIASCOR.TABLEVAR.IDSURVEY ;   
  if ( strcmp(varName,"SIM_alpha") == 0 ) { 
    ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_ALPHA ; NBMAX=MXa; 
    IPAR_ABG = IPAR_ALPHA0;
  }
  else if ( strcmp(varName,"SIM_beta") == 0 ) { 
    ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_BETA ;  NBMAX=MXb; 
    IPAR_ABG = IPAR_BETA0;
  }
  else if ( strcmp(varName,"SIM_gammaDM") == 0 ) { 
    ptrVal_f = INFO_BIASCOR.TABLEVAR.SIM_GAMMADM ; NBMAX=MXg; 
    IPAR_ABG = IPAR_GAMMA0;
  }
  else {
    sprintf(c1err,"Invalid varName = '%s' ", varName);
    sprintf(c2err,"Must be SIMalpha or SIMbeta");  
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);    
  }
  

  for(IDSURVEY=0; IDSURVEY < MXIDSURVEY; IDSURVEY++ ) {
    valmin_sample[IDSURVEY] = +1000.0 + IDSURVEY ;
    valmax_sample[IDSURVEY] = -1000.0 - IDSURVEY ;
  }


  print_debug_malloc(+1*debug_malloc, fnam);
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

  // if there is continuous gamma distribution, force 1 bing
  if ( IPAR_ABG==IPAR_GAMMA0 && NVAL > 5 ) { FORCE_ONEBIN = true ; }

  // if user input u[#] = 0 or 3, then force 1 bin (any parameter)
  ipar = INPUTS.ipar[IPAR_ABG];
  if ( ipar==0 || ipar == 3 )  { FORCE_ONEBIN = true ; }

  // implement FORCE_ONEBIN
  if ( FORCE_ONEBIN ) { NVAL = 1; val2nd = val_last; }

  // if floated param, add error check on number bins
  if ( NVAL > NBMAX && ipar > 0 ) {
    sprintf(c1err,"%d %s values exceeds bound of %d",
	    NVAL, varName, NBMAX);
    sprintf(c2err,"Check biasCor sample.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
  }

  // ---------------------------
  if ( NVAL < 3 ) { 
    valbin_loc = val2nd - val1st ;
    valmin_loc = val1st   - 0.5*valbin_loc ;
    valmax_loc = val_last + 0.5*valbin_loc ;
      
    if ( FORCE_ONEBIN ) { valbin_loc = valmax_loc - valmin_loc; }

    // check that GRID min & max is the same for each IDSAMPLE
    if ( !FORCE_ONEBIN ) 
      { check_abg_minmax_biasCor(varName,valmin_sample,valmax_sample); }

  }
  else {
    sprintf(c1err,"Invalid NVAL=%d for %s", NVAL, varName);
    sprintf(c2err,"val(1st,2nd) = %f, %f", val1st, val2nd);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
  }

  /* xxxx
  printf(" xxx ------------------------------------------ \n");
  printf(" xxx %s: Found %d distinct %s values \n",
	 fnam, NVAL, varName  );
  printf(" xxx %s val[min,max,bin]_loc = %.4f to %.4f, bin=%.4f\n",
	 varName, valmin_loc, valmax_loc, valbin_loc );
  printf(" xxx va1st=%.3f  val2nd=%.3f  val_last=%.3f\n", 
	 val1st, val2nd, val_last );
  printf(" xxx %s: IPAR_ABG=%d  FORCE_ONEBIN = %d \n", 
	 fnam, IPAR_ABG, FORCE_ONEBIN);
  fflush(stdout);
  xxxxxx */

  print_debug_malloc(-1*debug_malloc,fnam);
  free(INDEX_UNSORT); // free memory
  
  // --------------------------
  // load function args
  *VAL_MIN = valmin_loc ;
  *VAL_MAX = valmax_loc ;
  *VAL_BIN = valbin_loc ;

  return ;

} // end get_BININFO_biasCor_abg


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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
  }

  return ;
} // end check_abg_minmax_biasCor

// ================================================
void get_muBias(char *NAME,
		BIASCORLIST_DEF *BIASCORLIST, 
		FITPARBIAS_DEF (*FITPARBIAS_ABGRID)[MXb][MXg],
		double         (*MUCOVSCALE_ABGRID)[MXb][MXg],
		double         (*MUCOVADD_ABGRID)[MXb][MXg],
		INTERPWGT_AlphaBetaGammaDM *INTERPWGT,  
		double *fitParBias,
		double *muBias, double *muBiasErr, 
		double *muCOVscale, double *muCOVadd ) { 

  
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
  //  muCOVadd    = floor to apply to muErr Aug 2 2021 Dillon

  double alpha    = BIASCORLIST->alpha ;
  double beta     = BIASCORLIST->beta  ;
  double gammadm  = BIASCORLIST->gammadm  ;
  double z        = BIASCORLIST->z  ;
  double mB       = BIASCORLIST->FITPAR[INDEX_mB] ;
  double x1       = BIASCORLIST->FITPAR[INDEX_x1] ;
  double c        = BIASCORLIST->FITPAR[INDEX_c] ;

  double muBias_local     = 0.0 ;
  double muBiasErr_local  = 0.0 ;
  double muCOVscale_local = 0.0 ;
  double muCOVadd_local = 0.0 ;

  bool DO_COVSCALE = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD  ) > 0;

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
	  errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
	}
      
	for(ipar = ILCPAR_MIN; ipar <= ILCPAR_MAX ; ipar++ ) {
	  VAL = FITPARBIAS_ABGRID[ia][ib][ig].VAL[ipar];
	  ERR = FITPARBIAS_ABGRID[ia][ib][ig].ERR[ipar];
	  if ( VAL > 600.0 || isnan(ERR) ) {
	    sprintf(c1err,"Undefined %s-bias = %.3f +- %.3f for CID=%s", 
		    BIASCOR_NAME_LCFIT[ipar], VAL, ERR, NAME);
	    sprintf(c2err,"ia,ib,ig=%d,%d,%d  a,b,g=%.3f,%.2f,%.3f",
		    ia, ib, ig, alpha, beta, gammadm );
	    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
	  }
	  
	  WGTpar = WGTabg ;
	  WGTpar_SUM[ipar] += WGTpar;
	  biasVal[ipar] += ( WGTpar * VAL );
	  biasErr[ipar] += ( WGTpar * ERR ); 
	} // end ipar
	
	// now the COV scale (Jul 1 2016)
	VAL = MUCOVSCALE_ABGRID[ia][ib][ig] ;
	muCOVscale_local += ( WGTabg * VAL ) ;

	if ( DO_COVADD ) {
	  VAL = MUCOVADD_ABGRID[ia][ib][ig] ;
	  muCOVadd_local += ( WGTabg * VAL ) ;
	}	
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
  }

  *muBias     = muBias_local ;
  *muBiasErr  = muBiasErr_local ;
  *muCOVscale = muCOVscale_local ;
  *muCOVadd   = muCOVadd_local ;

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
  // Sep 28 2020: check option to use "same" file(s) as for biasCor
  // Feb 25 2021: for H11, set NPASS_CUTMASK_POINTER so that writing
  //              YAML file later doesn't crash.

  int  EVENT_TYPE   = EVENT_TYPE_CCPRIOR ;
  int  NSAMPLE      = NSAMPLE_BIASCOR ;
  int  NDIM_BIASCOR = INFO_BIASCOR.NDIM ;
  int idsample, USE_CCPRIOR_H11 ;
  char fnam[] = "prepare_CCprior" ;

  // ------------- BEGIN -------------

  // check for H11 polynomial prior that does use sim files
  USE_CCPRIOR_H11 = INFO_CCPRIOR.USEH11;

  if ( INPUTS.nfile_CCprior == 0 ) { return; }
  
  if ( NDIM_BIASCOR == 1 ) {
    sprintf(c1err,"Cannot use 1D biasCor with CC likelihood.");
    sprintf(c2err,"Either remove CC term, or use 5D biasCor.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
  }


  if ( USE_CCPRIOR_H11 ) { 
    sprintf(BANNER,"%s: use CC mu-vs-z prior from Hlozek 2011", fnam);
    fprint_banner(FP_STDOUT,BANNER);    
    NPASS_CUTMASK_POINTER[EVENT_TYPE] = &INFO_CCPRIOR.TABLEVAR.NSN_PASSCUTS ;
    INFO_CCPRIOR.TABLEVAR.NSN_PASSCUTS = 0;
    return ;
  }

  sprintf(BANNER,"%s: read simulated CC mu-prior vs. z", fnam);
  fprint_banner(FP_STDOUT,BANNER);


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

  // ----------------------------------------------
  fprintf(FP_STDOUT, "\n Finished preparing CC prior. \n");
  fflush(FP_STDOUT);

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
  int  LEN_MALLOC, NROW, ifile, ISTART, isn, IFILETYPE, NVAR_ORIG ;
  char *simFile ;
  //  char fnam[] = "read_simFile_CCprior" ;

  // -------------------- BEGIN --------------------

  // do quick read of each file get NEVT count for each file
  NEVT_TOT = 0 ;
  for(ifile=0; ifile < NFILE; ifile++ ) {
    simFile = INPUTS.simFile_CCprior[ifile];
    NEVT[ifile] = SNTABLE_NEVT(simFile,TABLENAME_FITRES); 
    NEVT_TOT += NEVT[ifile];
    fprintf(FP_STDOUT, "\t Found %d events in %s. \n", NEVT[ifile], simFile) ;
  }


  // malloc enough memory to read all files
  LEN_MALLOC = NEVT_TOT + 10;
  INFO_CCPRIOR.TABLEVAR.LEN_MALLOC = LEN_MALLOC ; ;
  malloc_INFO_CCPRIOR(+1, LEN_MALLOC, 0);

  CUTMASK_POINTER[EVENT_TYPE]       = &INFO_CCPRIOR.TABLEVAR.CUTMASK[0];
  IDSAMPLE_POINTER[EVENT_TYPE]      = &INFO_CCPRIOR.TABLEVAR.IDSAMPLE[0];
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
    compute_CUTMASK(isn, &INFO_CCPRIOR.TABLEVAR ); 
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

  fprintf(FP_STDOUT, "  %s: %d(ALL) -> %d(CUTS) \n", fnam, NSN_ALL, NSN_PASS);
  fflush(FP_STDOUT);

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
  fprintf(FP_STDOUT, "  %s: free memory for ALL \n", fnam); 
  fflush(FP_STDOUT);
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
  fprintf(FP_STDOUT, "\n");
  fprintf(FP_STDOUT, "  %s: %d z-bins from %.3f to %.3f\n",
	 fnam, ZBIN->nbin, INPUTS.zmin, INPUTS.zmax);

  fprintf(FP_STDOUT, "\t   iz  zmin   zmax     N(NONIA) \n" );
  for(iz=0; iz < ZBIN->nbin; iz++ ) {
    fprintf(FP_STDOUT, "\t  %3d  %.3f  %.3f  %5d  \n",
	   iz, ZBIN->lo[iz], ZBIN->hi[iz], ZBIN->n_perbin[iz] );
  }

  fflush(FP_STDOUT);

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
    fprintf(FP_STDOUT, "\n");
    fprintf(FP_STDOUT, "  %s: %d dMU-bins from %.3f to %.3f\n",
	   fnam, MUZMAP->DMUBIN.nbin, DMUMIN_CCPRIOR, DMUMAX_CCPRIOR);
    fflush(FP_STDOUT);
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
  double           MUCOVADD_TMP[MXa][MXb][MXg] ; 
  double fitParBias[NLCPAR], muBias, muBiasErr, muCOVscale, muCOVadd ;

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
		   MUCOVADD_TMP,         // (I) muCOVscale vs. ia,ib
		   &INTERPWGT,             // (I) wgt at each a,b grid point
		   fitParBias,     // (O) interp bias on mB,x1,c
		   &muBias,        // (O) interp bias on mu
		   &muBiasErr,     // (O) stat-error on above
		   &muCOVscale,  // (O) scale bias on muCOV (not used below) 
		   &muCOVadd );  // (O) add bias on muCOV (not used below)
      }
      else if ( NDIM_BIASCOR == 1 ) {
	debugexit("CCPRIOR does NOT WORK WITH BBC-1D");  // Jun 19 2018
	// muBias   =  ?? 	
      }
	
    } // end biasCor if-block
    

    // compute mucos for each SIM CC event
    dl       = cosmodl_forFit(z, z, MUZMAP->cosPar) ;
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
      MUZMAP->DMURMS[IDSAMPLE][iz] = STD_from_SUMS(NTMP,S1TMP,S2TMP);
	
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

  fprintf(FP_STDOUT, "\n  %s : ", SAMPLE_BIASCOR[IDSAMPLE].NAME  );
  fprintf(FP_STDOUT, "  a=%.3f b=%.3f \n",
	 MUZMAP->alpha, MUZMAP->beta );

  fprintf(FP_STDOUT, "   <DMU>=%.3f  RMS(DMU)=%.3f \n"
	 , MUZMAP->DMUAVG[IDSAMPLE][IZ] 
	 , MUZMAP->DMURMS[IDSAMPLE][IZ]
	 );
  /*
  printf("  OL,Ok,w0,wa = %.3f, %.3f, %.3f, %.3f \n",	 
	 MUZMAP->cosPar[0], MUZMAP->cosPar[1],
	 MUZMAP->cosPar[2], MUZMAP->cosPar[3]  );  */
  

  fprintf(FP_STDOUT, "\t  imu  dMU(min) dMU(max)     PDF(z=%.2f)\n", z );

  for(imu=0; imu < MUZMAP->DMUBIN.nbin; imu++ ) {
    NCC = MUZMAP->NDMUPDF[IDSAMPLE][IZ][imu] ; 
    PDF = MUZMAP->DMUPDF[IDSAMPLE][IZ][imu] ;
    fprintf(FP_STDOUT, "\t  %3d  %6.3f  %6.3f       %.3f  (%d) \n"
	   ,imu, MUZMAP->DMUBIN.lo[imu], MUZMAP->DMUBIN.hi[imu]
	   ,PDF, NCC);
  }  
  fflush(FP_STDOUT);

  return ;
  
} // end dump_DMUPDF_CCprior


// =====================================================
void  setup_contam_CCprior(char *which, CONTAM_INFO_DEF *CONTAM_INFO) {

  // Sep 2020
  // setup bins for storing binned contamination info vs. 
  // *which = 'MURES' or *which = 'REDSHIFT'.
  // This information is used to monitor contamination after the fit,
  // but is not used in the fit.

  int nb, i;
  double lo[MXz], hi[MXz];
  double tmp_lo, tmp_hi, tmp_avg, tmp_binsize;
  char fnam[] = "setup_contam_CCprior";

  // -------------- BEGIN ----------------

  zero_contam_CCprior(CONTAM_INFO) ;

  sprintf(CONTAM_INFO->BININFO.varName, "%s", which);  
  nb=0;

  if ( strcmp(which,"MURES") == 0 ) {
    lo[nb] = -4.0;  hi[nb] = -0.5;  nb++ ; 
    lo[nb] = -0.5;  hi[nb] = +0.5;  nb++ ;
    lo[nb] = +0.5;  hi[nb] = +4.0;  nb++ ;
  }
  else if ( strcmp(which,"REDSHIFT") ==0 ) {
    nb=4;
    double zbin = (INPUTS.zmax - INPUTS.zmin) / (double)nb;
    for(i=0; i < nb; i++  ) {
      lo[i]  = INPUTS.zmin + zbin*(double)i;
      hi[i]  = lo[i] + zbin;
    }
  }
  else {
    sprintf(c1err,"Invalid which = '%s'", which);
    sprintf(c2err,"Valid which options: MURES and REDSHIFT");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
  }
  
  for(i=0; i < nb; i++  ) {
    tmp_lo  = lo[i]; tmp_hi=hi[i];  
    tmp_binsize = tmp_hi - tmp_lo;
    tmp_avg = 0.5*(tmp_hi+tmp_lo);
    CONTAM_INFO->BININFO.lo[i]      = tmp_lo;
    CONTAM_INFO->BININFO.hi[i]      = tmp_hi;
    CONTAM_INFO->BININFO.avg[i]     = tmp_avg;
  }
  tmp_binsize = lo[1] - lo[0];
  CONTAM_INFO->BININFO.nbin    = nb;
  CONTAM_INFO->BININFO.binSize = tmp_binsize;

  return ;

} // end setup_contam_CCprior

void zero_contam_CCprior(CONTAM_INFO_DEF *CONTAM_INFO) {

  int i;
  int nbin = MXz;

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


void print_contam_CCprior(FILE *fp) {

  // Created Sep 2020
  // Driver function to compute CC contamination in bins of 
  // MURES and REDSHIFT, and print table to file pointer fp.
  // This info is for diagnostics and was not used in the fit.
  //
  // E.g.
  //   varName Range    CC/TOT(SUMPROB)     CC/TOT(NTRUE)
  //    ALL  
  //   -4.0 to -0.5
  //   -0.5 to +0.5
  //

  CONTAM_INFO_DEF CONTAM_MURES_BINS;
  CONTAM_INFO_DEF CONTAM_REDSHIFT_BINS;

  int   NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ; 
  int   i, SIM_NONIA_INDEX;
  //  char fnam[] = "print_contam_CCprior";

  // -------------- BEGIN --------------

  setup_contam_CCprior("MURES",    &CONTAM_MURES_BINS);
  setup_contam_CCprior("REDSHIFT", &CONTAM_REDSHIFT_BINS);

  int cutmask;  double zhd, mures, probcc, probIa ;
  for ( i = 0; i < NSN_DATA; i++ ) {
    cutmask  = INFO_DATA.TABLEVAR.CUTMASK[i] ;
    if ( cutmask ) { continue; }

    zhd      = INFO_DATA.TABLEVAR.zhd[i] ;
    mures    = INFO_DATA.mures[i];
    probcc   = INFO_DATA.probcc_beams[i] ;
    probIa   = 1.0 - probcc ;
    SIM_NONIA_INDEX = INFO_DATA.TABLEVAR.SIM_NONIA_INDEX[i];

    sum_contam_CCprior(&CONTAM_MURES_BINS,    probIa, mures,
		       SIM_NONIA_INDEX); 
    sum_contam_CCprior(&CONTAM_REDSHIFT_BINS, probIa, zhd,
		       SIM_NONIA_INDEX); 
  }

  print_table_CONTAM_INFO(fp, &CONTAM_MURES_BINS);
  print_table_CONTAM_INFO(fp, &CONTAM_REDSHIFT_BINS);

} // end print_contam_CCprior 

void print_table_CONTAM_INFO(FILE *fp,  CONTAM_INFO_DEF *CONTAM_INFO) {


  bool IS_SIM   = INFO_DATA.TABLEVAR.IS_SIM ;
  int  nbin     = CONTAM_INFO->BININFO.nbin ;
  double lo, hi, xnIa, xncc, ratio=-9.0, true_ratio = -9.0 ;
  int    ntrue_Ia, ntrue_cc, i ;
  char *varName = CONTAM_INFO->BININFO.varName;
  char cRange[40], str_contam_data[80], str_contam_true[80];

  // - - - - - -
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

  return ;

} // print_table_CONTAM_INFO

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
  double z, zmuerr, zsq ;
  double prob  = 0.0 ;
  double arg, sqsigCC, sigCC, sigintCC, CCpoly, dmuCC ;
  double alpha, beta, gamma, muerrsq_raw, covmat_fit[NLCPAR][NLCPAR] ;
  int i,i2;
  char fnam[] = "prob_CCprior_H11" ;

  // ----------- BEGIN ------------

  name    = INFO_DATA.TABLEVAR.name[n];
  z       = INFO_DATA.TABLEVAR.zhd[n];
  zmuerr  = INFO_DATA.TABLEVAR.zmuerr[n];
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
			    z,zmuerr, 0);

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
  // OPT == 0  -> return =9 if VAL is outside range
  // OPT == 1  -> abort if VAL is outside range
  // OPT == 2  -> if VAL is outside, return edge bin (do not abort)
  // OPT == 7  -> dump
  
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, MSG_ABORT);   
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
  // Note that NPASS is computed and stored here.
  // 
  // Tricky part is counting how many events are rejected by
  // a single cut ... and to to it efficiently.
  //
  // Jul 1 2020: abort if NEVT(biasCor)=0 for any IDSAMPLE
  // Feb 25 2021: abort on missing biasCor only if it is required.
  // Jun 08 2021: refactor to compute NSN_PASS and NSN_REJ here ...
  //              needed for SUBPROCESS iterations.
  //
  // Sep 18 2021: increment NPASS_CUTMASK_BYSAMPLE for YAML print out later.

  int  NSN_TOT      = *NALL_CUTMASK_POINTER[event_type];
  char *STRTYPE     = STRING_EVENT_TYPE[event_type];
  int  NSN_REJ=0, NSN_PASS=0, *NBIT, bit, isn, cutmask, idsample, NCUT=0 ;
  int  *CUTMASK_PTR, NCUT_SOLO[MXCUTBIT];
  short int *IDSAMPLE_PTR ;
  char fnam[] = "print_eventStats" ;

  // ---------- BEGIN ------------

  CUTMASK_PTR   =  CUTMASK_POINTER[event_type];
  IDSAMPLE_PTR  =  IDSAMPLE_POINTER[event_type];
  NBIT          =  &NSTORE_CUTBIT[event_type][0];

  for(idsample=0; idsample < NSAMPLE_BIASCOR; idsample++ ) 
    { NPASS_CUTMASK_BYSAMPLE[event_type][idsample] = 0; }

  // loop over all events and for each cutbit, count how many
  // events were cut by ONLY this one cut (i.e.. passed all other cuts)
  for(bit=0; bit < MXCUTBIT; bit++ ) { NCUT_SOLO[bit]=0; }

  for(isn=0; isn < NSN_TOT; isn++ ) {
    cutmask  = CUTMASK_PTR[isn];
    idsample = IDSAMPLE_PTR[isn];
    if ( cutmask ==0 ) {       
      NPASS_CUTMASK_BYSAMPLE[event_type][idsample]++ ;
      NSN_PASS++ ; 
      continue; 
    }
    // check if one and only 1 bit is set
    if ( (cutmask & (cutmask-1) ) == 0 ) {
      for(bit=0; bit < MXCUTBIT; bit++ ) {
	if ( cutmask & CUTMASK_LIST[bit] ) 
	  { NCUT_SOLO[bit]++; bit=MXCUTBIT+1; }
      }
    }
  }

  // ----------------------------
  NSN_REJ = NSN_TOT - NSN_PASS;
  *NREJECT_CUTMASK_POINTER[event_type] = NSN_REJ;
  *NPASS_CUTMASK_POINTER[event_type]   = NSN_PASS ;


  fprintf(FP_STDOUT, 
	  "\n#%s\n", dashLine);
  fprintf(FP_STDOUT, 
	  " %s STAT SUMMARY: %d(TOTAL) = %d(ACCEPT) + %d(REJECT) \n",
	 STRTYPE, NSN_TOT, NSN_PASS, NSN_REJ);

  for(bit=0; bit < MXCUTBIT; bit++ ) {
    NCUT = NBIT[bit] ;
    if ( NCUT > 0 ) {
      fprintf(FP_STDOUT, " %s NCUT[2**%2.2d=%5d] = "
	      "%6d(ALL) %6d(onlyCut)   [%s] \n",
	      STRTYPE, bit, CUTMASK_LIST[bit], 
	      NCUT, NCUT_SOLO[bit], CUTSTRING_LIST[bit] );
    }
  }

  if ( NSN_PASS == 0 ) {
    sprintf(c1err,"All %s events fail cuts.", STRTYPE);
    sprintf(c2err,"Check NCUT vs. cut listed above.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  fprintf(FP_STDOUT, "#%s\n\n", dashLine);
  fflush(FP_STDOUT);


  // July 1 2020:
  // Abort if BIASCOR sample is missing for any IDSAMPLE.
  int NMISSING = 0 ;
  if ( event_type == EVENT_TYPE_BIASCOR ) {
    int  NSAMPLE = NSAMPLE_BIASCOR ;
    int idsample, NSN ;      char *NAME;
    bool IS_PHOTOZ ;
    for(idsample=0 ; idsample < NSAMPLE; idsample++ ) {
      if ( SAMPLE_BIASCOR[idsample].DOFLAG_BIASCOR==0 ) { continue ; } 
      NSN        = SAMPLE_BIASCOR[idsample].NSN[event_type] ;
      IS_PHOTOZ  = SAMPLE_BIASCOR[idsample].IS_PHOTOZ ;
      NAME       = SAMPLE_BIASCOR[idsample].NAME ;
      if ( NSN == 0 ) {
	NMISSING++ ;
	printf(" WARNING: No BIASCOR events for %s [IDSAMPLE=%d]\n", 
	       NAME, idsample);
	printf("\t IS_PHOTOZ(data) = %d for %s \n", IS_PHOTOZ, NAME );
	fflush(stdout);
      }
    }

    if ( NMISSING > 0 ) {
      sprintf(c1err,"%d missing biasCor samples; see WARNINGS above.", 
	      NMISSING);
      sprintf(c2err,"Check data and biasCor samples for mis-match.") ; 
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
    }
  }

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
  // May 28 2021: fix COV cut bug; cut on cov(mB,X) instead of cov(x0,X)
  //
  int  event_type = TABLEVAR->EVENT_TYPE;
  bool IS_DATA    = ( event_type == EVENT_TYPE_DATA );
  bool IS_BIASCOR = ( event_type == EVENT_TYPE_BIASCOR );
  bool IS_CCPRIOR = ( event_type == EVENT_TYPE_CCPRIOR );
  int NFIELD = INPUTS.NFIELD ;

  bool  LCUTWIN_DISABLE   = IS_DATA && INPUTS.LCUTWIN_DISABLE ;

  //  int  LDMP = 0;
  int  DOFLAG_CUTWIN[MXCUTWIN], icut, outside ;
  int  CUTMASK, REJECT, ACCEPT ;
  int  sntype, FOUND, SIM_NONIA_INDEX, idsample, idsurvey, IZBIN ;
  bool BADERR=false, BADCOV=false ;
  double cutvar_local[MXCUTWIN];
  double z, x1, c, logmass, x0err, x1err, cerr  ;
  double COV_mBx1, COV_mBc, COV_x1c,  mBerr ;
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
  idsurvey  =  (int)TABLEVAR->IDSURVEY[isn] ;
  sntype    =  (int)TABLEVAR->SNTYPE[isn] ;
  
  z         =  (double)TABLEVAR->zhd[isn];
  x1        =  (double)TABLEVAR->fitpar[INDEX_x1][isn] ;
  c         =  (double)TABLEVAR->fitpar[INDEX_c ][isn] ;  
  logmass   =  (double)TABLEVAR->host_logmass[isn];

  // undo temporary bug that wrote logmass values of -9999 (Feb 2022)
  if ( logmass < -9990.0 ) 
    { TABLEVAR->host_logmass[isn] = logmass = -9.0; }

  x0err     =  (double)TABLEVAR->x0err[isn] ;
  mBerr     =  (double)TABLEVAR->fitpar_err[INDEX_mB][isn] ;
  x1err     =  (double)TABLEVAR->fitpar_err[INDEX_x1][isn] ;
  cerr      =  (double)TABLEVAR->fitpar_err[INDEX_c ][isn] ;
    
  COV_mBx1  = (double)TABLEVAR->covmat_fit[isn][INDEX_mB][INDEX_x1];     
  COV_mBc   = (double)TABLEVAR->covmat_fit[isn][INDEX_mB][INDEX_c];
  COV_x1c   = (double)TABLEVAR->covmat_fit[isn][INDEX_x1][INDEX_c];

  SIM_NONIA_INDEX = (int)TABLEVAR->SIM_NONIA_INDEX[isn] ;


  // =======================================
  // check CUTWIN options; reject_CUTWIN returns 
  // 0->accepted, 1->rejected, 2->dweight with large MUERR

  REJECT = reject_CUTWIN(event_type, DOFLAG_CUTWIN, cutvar_local);
  if ( REJECT == DOFLAG_CUTWIN_APPLY ) 
    { setbit_CUTMASK(isn, CUTBIT_CUTWIN, TABLEVAR); }

  else if ( REJECT == DOFLAG_CUTWIN_FITWGT0 )  {
    // Jan 2021: set flag to set fit wgt ~ 0 by setting MUERR -> large
    INFO_DATA.set_fitwgt0[isn]   = true; 
  }

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
  // (might have to always apply z-cut, even if DISABLE flag is set??)

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
      fprintf(FP_STDOUT, " %s WARNING: "
	      "one or more x0,x1,c errors <=0 (SNID=%s)\n", fnam, name );
    }
  }

  
  // May 2016  flag crazy COV and mBerr 
  if ( fabs(COV_mBx1) > 8.0 ) { BADCOV = true ; } 
  if ( fabs(COV_mBc ) > 8.0 ) { BADCOV = true ; } 
  if ( fabs(COV_x1c ) > 8.0 ) { BADCOV = true ; } 
  if ( mBerr > 1.0  )         { BADCOV = true ; }
  
  if ( BADCOV ) {
    setbit_CUTMASK(isn, CUTBIT_BADCOV, TABLEVAR );
    if ( IS_DATA ) {
      fprintf(FP_STDOUT," %s WARNING: one or more crazy COV (SNID=%s). \n",
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
    ACCEPT = selectCID_data(name, idsurvey, &IZBIN);
    if ( !ACCEPT )
      { setbit_CUTMASK(isn, CUTBIT_CID, TABLEVAR); } //the mask is in tablevar

    if ( INFO_DATA.USE_IZBIN_from_CIDFILE ) 
      { INFO_DATA.TABLEVAR.IZBIN[isn] = IZBIN; }

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
  bool IS_DATA      = ( EVENT_TYPE == EVENT_TYPE_DATA );
  int  *CUTMASK     = &TABLEVAR->CUTMASK[isn];
  int  CUTMASK_SET  = CUTMASK_LIST[bitnum] ;
  bool LCUTWIN_DISABLE  = IS_DATA && INPUTS.LCUTWIN_DISABLE ;
  bool APPLY ;
  //  char fnam[] = "setbit_CUTMASK" ;

  // ------------- BEGIN -------------

  if ( LCUTWIN_DISABLE ) {
    // if cuts are disabled for data, only allow cuts on
    // CID (read from list), biasCor, and few others
    APPLY = 
      (bitnum == CUTBIT_CID     ) || 
      (bitnum == CUTBIT_BIASCOR ) ||
      (bitnum == CUTBIT_z       ) ||
      (bitnum == CUTBIT_BADCOV  ) 
      ;
    if ( !APPLY ) { return; }
  }

  // - - - - -
  if ( ( *CUTMASK & CUTMASK_SET) == 0 ) {

    if ( *CUTMASK == 0 ) { TABLEVAR->NSN_REJECT++ ;  }

    *CUTMASK |= CUTMASK_LIST[bitnum] ;
    TABLEVAR->NSN_CUTBIT[bitnum]++ ;
    NSTORE_CUTBIT[EVENT_TYPE][bitnum]++ ;
  }

  return ;

} // end setbit_CUTMASK 


// ===========================================
int selectCID_data(char *cid, int IDSURVEY, int *IZBIN){

  // Created Sep 5 2019 by D.Brout
  // for file= data. determines if cid is in cidlist_data
  //
  // Sep 2020 RK - Refactor to accept or reject based on user input.
  // Jun 2021 RK - use match_cidlist_exec util based on hash table.
  // Jun 2021 DB - added boolean logic to match on IDSURVERY
  // Apr 2022 RK - return IZBIN if it is part of cid_select table

  
  int ncidList   = INPUTS.ncidList_data ;
  int acceptFlag = INPUTS.acceptFlag_cidFile_data ;
  bool match_on_cid_idsurvey = INPUTS.match_on_cid_idsurvey;
  bool match_on_cid_only = INPUTS.match_on_cid_only;
  int ACCEPT = 1, REJECT = 0, i, isn_match ;
  bool MATCH ;
  char *tmpCID, STRINGID[60];
  char fnam[] = "selectCID_data";

  // ------- BEGIN -------------

  *IZBIN = -9 ;

  if ( ncidList == 0 ) { return ACCEPT ; }

  if (match_on_cid_only) {
      sprintf(STRINGID,"%s",cid);
    }
  else if (match_on_cid_idsurvey) {
      sprintf(STRINGID,"%s_%d",cid,IDSURVEY);
    }
  else {
    sprintf(c1err,"Boolean logic failed, unable to match.");
    sprintf(c2err,"match_on_cid_only=false and match_on_cid_idsurvey=false");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  isn_match = match_cidlist_exec(STRINGID);
  MATCH     = (isn_match >= 0);
  if ( INFO_DATA.USE_IZBIN_from_CIDFILE && MATCH ) {
    *IZBIN = match_cidlist_parval(isn_match, VARNAME_IZBIN, 1);

    if ( INPUTS.LEGACY_IZBIN ) {
      int IZBIN_old = INFO_DATA.IZBIN_from_CIDFILE[isn_match];
      if ( IZBIN_old != *IZBIN ) {
	sprintf(c1err,"IZBIN(old,new) = %d, %d for CID=%s",
		IZBIN_old, *IZBIN, cid);
	sprintf(c2err,"something wrong with new IZBIN");
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
      }
    } // end LEGACY_IZBIN

  }

  /* xxx
  if ( strstr(STRINGID,"9159") != NULL ) {
    printf(" xxx %s: STRINGID='%s'  MATCH=%d \n", fnam, STRINGID, MATCH);
  }
  xxx */

  // - - - -
  if ( MATCH ) {
    if ( acceptFlag > 0 )  { return ACCEPT; }
    else     	           { return REJECT; }
  }
  else {
    // if there are no cid matches to file:
    //   REJECT if mode is to accept only CIDs in file
    //   ACCEPT if mode is to reject only CIDs in file
    if ( acceptFlag > 0 ) { return REJECT ; }
    else                  { return ACCEPT ; }
  }

} // END selectCID_data

// =============================================
int prescale_reject_simData(int SIM_NONIA_INDEX) {

  // Created Apr 14 2017 by R.Kessler
  // For simulated data, check prescales to see if this event is rejected:
  //   * prescale_simData pre-scales all sim events randomly
  //   * prescale_simCC   pre-scales only the true CC subset
  //   * prescale_simIa   pre-scales only the true Ia subset
  // Note that real data are always accepted.
  // Function returns 1 to reject, 0 to keep.
  // Prescale can be float, e.g., 2.6.
  //
  // Do not use this function for biasCor or CCprior.
  //
  // Nov 25 2021: check prescale_simIa

  int REJECT = 0 ;
  float XN, XNPS;
  char fnam[] = "prescale_reject_simData" ;

  // ------------- BEGIN -----------------

  // return accept for real data
  if ( FOUNDKEY_SIM == 0 ) { return(REJECT); }

  // Increment NSIMDATA counter.
  // Note that prescale_simData can be non-integer
  NSIMDATA++ ;
  XN    = (float)NSIMDATA ;
  XNPS  = (float)INPUTS.prescale_simData ;
  if ( fmodf( XN, XNPS ) != 0 )  { REJECT = 1 ; }

  // increment separate NSIMCC counter
  if ( SIM_NONIA_INDEX > 0 ) {
    NSIMCC++ ;
    XN    = (float)NSIMCC ;
    XNPS  = (float)INPUTS.prescale_simCC ;
    if ( fmodf( XN, XNPS ) != 0 ) { REJECT = 1 ; }
    if ( INPUTS.prescale_simCC > 9999.0 ) { REJECT=1 ; }
  }

  if ( SIM_NONIA_INDEX == 0 ) {
    NSIMIa++ ;
    XN    = (float)NSIMIa ;
    XNPS  = (float)INPUTS.prescale_simIa ;
    if ( fmodf( XN, XNPS ) != 0 ) { REJECT = 1 ; }
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  // Apr 27 2021: remove SKIP on colon to allow colon in comment field

  FILE *fdef;
  bool YAML;  
  char *sptr ;
  char fnam[] = "parse_parFile" ;

  // ------------------ BEGIN --------------

  if ( strstr(parFile,"cat_only") != NULL ) 
    { parse_cat_only(parFile); return; }



#ifdef USE_SUBPROCESS
  if ( strcmp(parFile,"SUBPROCESS_HELP") == 0 )  { SUBPROCESS_HELP(); }
#endif

  // allow for some null options to skip reading file
  if ( strcmp(parFile,"NULL")  == 0 ) return;
  if ( strcmp(parFile,"null")  == 0 ) return;
  if ( strcmp(parFile,"BLANK") == 0 ) return;
  
  if ( strcmp(parFile,STRINGMATCH_KEY_DUMP) == 0 ) {
    // prepare to dump all valid keys and quit.
    uniqueOverlap(STRINGMATCH_KEY_DUMP,"SALT2mu-input file"); 
    sprintf(parFile,"%s/SURVEY.DEF", PATH_SNDATA_ROOT);
    INPUTS.KEYNAME_DUMPFLAG = true ;
  }
  else {
    // normal init to read each key, and to abort on duplicate keys
    uniqueOverlap(STRINGMATCH_INIT,"SALT2mu-input file"); 
  }

  // Look for a file of default options
  fdef = fopen(parFile,"rt");
  if (!fdef) {
    sprintf(c1err,"could not open SALT2mu parameter-input file:");
    sprintf(c2err," %s ", parFile);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }


  fprintf(FP_STDOUT,
	  "Reading SALT2mu parameter-input file '%s' : \n",  parFile);

  char line[MXCHAR_LINE];
  YAML = false;

  while (fgets(line,MXCHAR_LINE,fdef)) {

    // skip blank lines and lines starting with comment char (RK 4/23/2012)
    if ( strlen(line) < 3 )       { continue ; }  
    if ( !YAML && commentchar(line) == 1 ) { continue ; } // see sntools.c

    sptr = strtok(line,"\n");
    
    // July 24 2020: skip YAML CONFIG 
    if ( strcmp(sptr,"CONFIG:")    == 0    ) { YAML = true  ; }
    if ( strstr(sptr,"END_YAML")   != NULL ) { YAML = false ; }
    if ( strstr(sptr,"END_CONFIG") != NULL ) { YAML = false ; }
    if ( strstr(sptr,"CONFIG_END") != NULL ) { YAML = false ; }
    if ( YAML ) { continue; }

    ppar(sptr); // pass entire line
  }

  fprintf(FP_STDOUT,"\n");
  fclose(fdef);

  if ( INPUTS.KEYNAME_DUMPFLAG ) { happyend(); }

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
  fprintf(FP_STDOUT, "\n Parse command-line options: \n");
  uniqueOverlap("INIT","SALT2mu command-line override");
  
  for (i=2; i < argc; ++i) {

    item = argv[i];
    ntmp++;


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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }

  } // end i loop 

  if ( ntmp > 0 ) 
    { fprintf(FP_STDOUT, "\n"); }
  else 
    { fprintf(FP_STDOUT,
	      " None found. Use all options from parameter-input file.\n"); 
    }

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
  //
  // Oct 14 2020: refactor to use parse_commaSepList utility
  //
  int  ipar, len, ikey, ntmp ;  
  char key[MXCHAR_VARNAME], *s, tmpString[60];
  char fnam[] = "ppar" ;

  // --------- BEGIN ----------
  
  if ( !INPUTS.KEYNAME_DUMPFLAG ) 
    { fprintf(FP_STDOUT, " Parse '%s' \n",item);  fflush(FP_STDOUT); }

  if ( uniqueOverlap(item,"prefix=") )
    { sscanf(&item[7],"%s",INPUTS.PREFIX); return(1); }

  // - - - 
  if ( uniqueOverlap(item,"catfile_out=") ) {  // legacy key name
    s = INPUTS.cat_file_out ;
    sscanf(&item[12],"%s",s); remove_quote(s); return(1);
  }
  if ( uniqueOverlap(item,"cat_file_out=") ) {
    s = INPUTS.cat_file_out ;
    sscanf(&item[13],"%s",s); remove_quote(s); return(1);
  }

  // Aug 2020: check for flag to write YAML formatted output; 
  // intended for submit_batch_jobs
  if ( uniqueOverlap(item,"write_yaml=") ) {
    sscanf(&item[11],"%d", &INPUTS.write_yaml);  return(1);
  }

  if ( uniqueOverlap(item,"write_csv=") ) {
    sscanf(&item[10],"%d", &INPUTS.write_csv);  return(1);
  }

  // - - - -

#ifdef USE_SUBPROCESS
  if ( uniqueOverlap(item,"SUBPROCESS_HELP") ) {
    SUBPROCESS_HELP();
  }
  if ( uniqueOverlap(item,"SUBPROCESS_FILES=") ) {
    SUBPROCESS_MALLOC_INPUTS();
    s = SUBPROCESS.INPUT_FILES ; // comma-sep list of INPFILE,OUTFILE
    SUBPROCESS.USE = true ;
    sscanf(&item[17],"%s",s); remove_quote(s); return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_VARNAMES_GENPDF=") ) {
    s = SUBPROCESS.INPUT_VARNAMES_GENPDF_STRING ; // comma-sep list of varnames
    sscanf(&item[27],"%s",s); remove_quote(s); return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_SIMREF_FILE=") ) {
    s = SUBPROCESS.INPUT_SIMREF_FILE ; // Name/location of input file used to generate sim reference 
    sscanf(&item[23],"%s",s); remove_quote(s); return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_CID_REWGT_DUMP=") ) {
    s = SUBPROCESS.INPUT_CID_REWGT_DUMP ; // comma-sep list of SNIDs
    sscanf(&item[26],"%s",s); remove_quote(s); return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_ISEED=") ) {
    sscanf(&item[17], "%d", &SUBPROCESS.INPUT_ISEED ); return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_STDOUT_CLOBBER=") ) {
    sscanf(&item[26], "%d", &SUBPROCESS.STDOUT_CLOBBER );  return(1);
  }
  if ( uniqueOverlap(item,"SUBPROCESS_NEVT_SIM_PRESCALE=") ) {
    sscanf(&item[29], "%d", &SUBPROCESS.NEVT_SIM_PRESCALE ); return(1);
  }
  if (  !strncmp(item,"SUBPROCESS_OUTPUT_TABLE=",24) ) {
    int N = SUBPROCESS.N_OUTPUT_TABLE ;
    sscanf(&item[24], "%s", SUBPROCESS.INPUT_OUTPUT_TABLE[N] ); 
    s = SUBPROCESS.INPUT_OUTPUT_TABLE[N]; remove_quote(s);
    //printf(" xxx %s: N_OUTPUT_TABLE = %d  VARDEF = '%s' \n", 
    //	   fnam, N, SUBPROCESS.INPUT_OUTPUT_TABLE[N] ); 
    SUBPROCESS.N_OUTPUT_TABLE++ ;
    return(1) ;
  }
  if ( uniqueOverlap(item,"SUBPROCESS_OPTMASK=") ) {
    sscanf(&item[19], "%d", &SUBPROCESS.INPUT_OPTMASK ); return(1);
  }

#endif

  if ( uniqueOverlap(item,"cutmask_write=") )
    { sscanf(&item[14],"%i", &INPUTS.cutmask_write ); return(1); }
  if ( uniqueOverlap(item,"errmask_write=") ) // allow legacy name
    { sscanf(&item[14],"%i", &INPUTS.cutmask_write ); return(1); }

  if ( uniqueOverlap(item,"minos=") ) 
    { sscanf(&item[6],"%i", &INPUTS.minos ); return(1); }


  // - - - - - -
  // allow two different keys for data file name
  char keyList_data[2][12] = { "file=", "datafile=" } ;
  for ( ikey=0; ikey < 2; ikey++ ) {
    len = strlen( keyList_data[ikey] );
    if ( uniqueOverlap(item,keyList_data[ikey]) ) {
      parse_commaSepList("DATAFILE", &item[len], MXFILE_DATA, MXCHAR_FILENAME, 
			 &INPUTS.nfile_data, &INPUTS.dataFile );
      return(1);
    }
  }

  if ( uniqueOverlap(item,"datafile_override=") ) {
    parse_commaSepList("DATAFILE_OVERRIDE", &item[18], 
		       MXFILE_BIASCOR, MXCHAR_FILENAME, 
		       &INPUTS.nfile_data_override, &INPUTS.dataFile_override);
    if ( IGNOREFILE(INPUTS.dataFile_override[0]) ) 
      { INPUTS.nfile_data_override = 0; }

    return(1);
  }

  // - - - - - - - 
  // allow two different keys to define biasCor file name
  char keyList_biasCor[2][20] = { "simfile_bias=", "simfile_biascor=" } ;
  for ( ikey=0; ikey < 2; ikey++ ) {
    len = strlen( keyList_biasCor[ikey] );
    if ( uniqueOverlap(item,keyList_biasCor[ikey]) ) {

      // save biasCor arg in case simfile_ccprior = 'same'
      print_debug_malloc(+2,fnam);
      INPUTS.simFile_biasCor_arg = (char*) malloc(strlen(item)*sizeof(char));
      sprintf(INPUTS.simFile_biasCor_arg, "%s", &item[len]);

      parse_commaSepList("SIMFILE_BIASCOR", &item[len], 
			 MXFILE_BIASCOR, MXCHAR_FILENAME, 
			 &INPUTS.nfile_biasCor, &INPUTS.simFile_biasCor );
      return(1);
    }
  }

  // - - - - - - 
  if ( uniqueOverlap(item,"prescale_biascor=") ) 
    { parse_prescale_biascor(&item[17],0); return(1); }
 
  if ( uniqueOverlap(item,"opt_biascor=")  )
    { sscanf(&item[12],"%d", &INPUTS.opt_biasCor);   return(1); }

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
  if ( uniqueOverlap(item,"fieldgroup_biascor=") ||
       uniqueOverlap(item,"fieldGroup_biascor=")     ) {
    s=INPUTS.fieldGroup_biasCor ;  
    sscanf(&item[19],"%s",s); remove_quote(s);
    return(1);
  }

  if ( uniqueOverlap(item,"surveygroup_biascor_abortflag=")  ) {
    sscanf(&item[30],"%d", &INPUTS.surveyGroup_biasCor_abortFlag); 
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
  if ( uniqueOverlap(item,"simfile_ccprior=") )   { 
    parse_simfile_CCprior(&item[16]); 
    return(1);  
  }

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

  //  if ( uniqueOverlap(item,"chi2max=")) 
  //  { sscanf(&item[80,"%le",&INPUTS.chi2max); return(1); }
  if ( !strncmp(item,"chi2max",7) )  // multiple chi2max keys allowed
    { parse_chi2max(item); return(1); }

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

  if ( uniqueOverlap(item,"sigint_step1=")) 
    { sscanf(&item[13],"%lf",&INPUTS.sigint_step1); return(1); }
  if ( uniqueOverlap(item,"dchi2red_dsigint=")) 
    { sscanf(&item[17],"%lf",&INPUTS.dchi2red_dsigint); return(1); }

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
    { sscanf(&item[16], "%s", INPUTS.idsample_select );  return(1); } 

  if ( uniqueOverlap(item,"select_trueIa=") ) 
    { sscanf(&item[14], "%d", &INPUTS.select_trueIa );  return(1); } 

  if ( uniqueOverlap(item,"force_realdata=") ) 
    { sscanf(&item[15], "%d", &INPUTS.force_realdata );  return(1); } 

  if ( uniqueOverlap(item,"zspec_errmax_idsample=") ) 
    { sscanf(&item[22], "%le", &INPUTS.zspec_errmax_idsample );  return(1); } 

  if ( uniqueOverlap(item,"cid_select_file=") )  {  
    parse_commaSepList("CID_SELECT_FILE", &item[16], 6, MXCHAR_FILENAME, 
		       &INPUTS.ncidFile_data, &INPUTS.cidFile_data );
    INPUTS.acceptFlag_cidFile_data = +1;   return(1); 
  }

  if ( uniqueOverlap(item,"izbin_from_cid_file=") ) 
    { sscanf(&item[19], "%d", &INPUTS.izbin_from_cidFile );  return(1); } 

  if ( uniqueOverlap(item,"cid_reject_file=") ) {
    parse_commaSepList("CID_REJECT_FILE", &item[16], 6, MXCHAR_FILENAME, 
		       &INPUTS.ncidFile_data, &INPUTS.cidFile_data );
    INPUTS.acceptFlag_cidFile_data = -1; return(1); 
  }

  if ( uniqueOverlap(item,"sntype=") )  
    { parse_sntype(&item[7]); return(1); }

  if ( uniqueOverlap(item,"nmax=") ) {
    s = INPUTS.nmaxString ;
    sscanf(&item[5],"%s",s);  remove_quote(s); return(1); } 

  if ( uniqueOverlap(item,"host_logmass_split="))  {
    sprintf(c1err,"Input host_logmass_split is obsolete !");
    sprintf(c2err,"Use p7 and u7 instead.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    { sscanf(&item[8],"%lf",&INPUTS.zpecerr);  return(1);   }

  if ( uniqueOverlap(item,"vpecerr="))  { 
    sscanf(&item[8],"%lf",&INPUTS.zpecerr); 
    INPUTS.zpecerr /= LIGHT_km ;  return(1); 
  }


  sprintf(key,"zwin_vpec_check=");
  if ( uniqueOverlap(item,key)) { 
    sscanf(&item[16], "%s", tmpString) ;
    char **str_zwin;
    parse_commaSepList(key, tmpString, 2, 10, &ntmp, &str_zwin);
    sscanf(str_zwin[0], "%le", &INPUTS.zwin_vpec_check[0]);
    sscanf(str_zwin[1], "%le", &INPUTS.zwin_vpec_check[1]);
    return(1); 
  }

  if ( uniqueOverlap(item,"pecv=")) {  // LEGACY key
    legacyKey_abort(fnam, "pecv", "zpecerr" ); 
  }

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
  if ( uniqueOverlap(item,"redchi2_tol="))  
    { sscanf(&item[12],"%lf",&INPUTS.redchi2_tol); return(1); }

  if ( uniqueOverlap(item,"prescale_simdata="))  
    { sscanf(&item[17],"%lf",&INPUTS.prescale_simData); return(1); }
  if ( uniqueOverlap(item,"simdata_prescale="))  // allow mental flip, Jan 2021
    { sscanf(&item[17],"%lf",&INPUTS.prescale_simData); return(1); }

  if ( uniqueOverlap(item,"prescale_simcc=")) 
    { sscanf(&item[15],"%lf",&INPUTS.prescale_simCC); return(1); }

  if ( uniqueOverlap(item,"prescale_simIa="))  
    { sscanf(&item[15],"%lf",&INPUTS.prescale_simIa); return(1); }
  if ( uniqueOverlap(item,"prescale_sim1a="))  // allow 1a or Ia
    { sscanf(&item[15],"%lf",&INPUTS.prescale_simIa); return(1); }


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
    { sscanf(&item[19],"%d", &INPUTS.ndump_nobiasCor ); return(1); }
  if ( uniqueOverlap(item,"ndump_nobiascor=")) 
    { sscanf(&item[16],"%d", &INPUTS.ndump_nobiasCor ); return(1); }

  if ( uniqueOverlap(item,"frac_warn_nobiascor=")) 
    { sscanf(&item[20],"%le", &INPUTS.frac_warn_nobiasCor ); return(1); }

  if ( uniqueOverlap(item,"snid_mucovdump=")) 
    { sscanf(&item[15],"%s", INPUTS.SNID_MUCOVDUMP); return(1); }

  if ( uniqueOverlap(item,"restore_sigz=")) 
    { sscanf(&item[13],"%d", &INPUTS.restore_sigz); return(1); }

  if ( uniqueOverlap(item,"restore_mucovscale_bug=")) 
    { sscanf(&item[23],"%d", &INPUTS.restore_mucovscale_bug); return(1); }
  if ( uniqueOverlap(item,"restore_mucovadd_bug="))
    { sscanf(&item[21],"%d", &INPUTS.restore_mucovadd_bug); return(1); }

  if ( uniqueOverlap(item,"debug_flag=")) { 
    sscanf(&item[11],"%d", &INPUTS.debug_flag); 
    if ( INPUTS.debug_flag==401 ) { INPUTS.izbin_from_cidFile=1; } 
    return(1); 
  }

  if ( uniqueOverlap(item,"debug_malloc=")) { 
    sscanf(&item[13],"%d", &INPUTS.debug_malloc); 
    INPUTS.debug_malloc *= 2;  // must be >1 to take effect
    return(1); 
  }

  if ( uniqueOverlap(item,"debug_mucovscale=")) 
    { sscanf(&item[17],"%d", &INPUTS.debug_mucovscale); return(1); }
  if ( uniqueOverlap(item,"nbinc_mucovscale="))
    { sscanf(&item[17],"%d", &INPUTS.nbinc_mucovscale); return(1); }
  if ( uniqueOverlap(item,"cidlist_debug_biascor=")) 
    { sscanf(&item[22],"%s", INPUTS.cidlist_debug_biascor); return(1); }

  if ( uniqueOverlap(item,"nthread=")) 
    { sscanf(&item[8],"%d", &INPUTS.nthread); return(1); }

  return(0);
  
} // end ppar

// **************************************************
void parse_cat_only(char *string_cat_only) {
  
  char fnam[]="parse_cat_only";

  // check for cat_only or cat_only/5 (where 5 = prescale)
  int LDMP=0;
  int OPT=0;
  

  INPUTS.cat_only = true ;


  // for exact cat_only match, exit with default PS=1

  // xxx  if ( strcmp(string_cat_only,"cat_only") == 0 ) { return; }
  if ( strstr(string_cat_only,"/") == NULL ) { return; }

  // check for prescale; e.g., cat_only/4.
  // Try python-like dictionary function parse_string_prescales

  STRING_DICT_DEF DICT_PRESCALE ;
  init_string_dict(&DICT_PRESCALE, "cat_only", 2);
 
  parse_string_prescales(string_cat_only, &DICT_PRESCALE);

  double ps = get_string_dict(OPT, "cat_only", &DICT_PRESCALE);

  INPUTS.cat_prescale = (int)ps;

  if ( LDMP ) {
    //printf("xxx %s: ps=%d",fnam,(int)ps);
    // print ps
    debugexit(fnam);
  }

} // end parse_cat_only

// **************************************************
void parse_simfile_CCprior(char *item) {

  // Created May 16 2019
  // Refactored Oct 2020 to use parse_commaSepList utility.
  // Parse comma-separate list of CCprior files names.
  //
  // Use item_local instead of item, so that we don't overwrite
  // input item and clobber other arguments.
  // 
  // Apr 11 2022: fix memory-overwrite bug related to MEMC value

  int ifile;
  int debug_malloc = INPUTS.debug_malloc ;
  int MEMC;
  char *item_local;
  char fnam[]  = "parse_simfile_CCprior" ;

  // ------------------ BEGIN -----------------

  print_debug_malloc(+1*debug_malloc,fnam);

  // 9.28.2020:check "same" option
  if ( strcmp(item,"same") == 0 ) {
    MEMC        = sizeof(char) * (strlen(INPUTS.simFile_biasCor_arg)+10);
    item_local  = (char*) malloc(MEMC);
    sprintf(item_local, "%s",  INPUTS.simFile_biasCor_arg); 
    INPUTS.sameFile_flag_CCprior = true; 
  }
  else {
    MEMC       = sizeof(char) * (strlen(item)+10);
    item_local = (char*) malloc(MEMC);
    sprintf(item_local, "%s", item); 
  }

  parse_commaSepList("SIMFILE_CCPRIOR", item_local, 
		     MXFILE_CCPRIOR, MXCHAR_FILENAME, 
		     &INPUTS.nfile_CCprior, &INPUTS.simFile_CCprior );

  char *f0 = INPUTS.simFile_CCprior[0] ;
  if ( IGNOREFILE(f0) ) { 
    INPUTS.ipar[IPAR_scalePCC]  = 0 ; // make sure scalePCC is not floated
  }
  else {
    INFO_CCPRIOR.USE = 1;
    if ( strcmp(f0,"H11") == 0 ) { INFO_CCPRIOR.USEH11 =  1; }
  }

  print_debug_malloc(-1*debug_malloc,fnam);
  free(item_local);

  return ;

} // parse_simfile_CCprior


// **************************************************     
void parse_cidFile_data(int OPT, char *fileName) {

  // Refactored April 2022 to read IZBIN.
  // Read inpt fileName for list of CIDs to accept or reject 
  // based on
  //
  //    OPT > 0 -> list to accept
  //    OPT < 0 -> list to reject
  //
  // If IDSURVEY column exists, match by CID_SURVEY
  // If IZBIN column exists, store it for use with event syncing (Apr 2022)

  int  ncidList_data = INPUTS.ncidList_data  ;
  int  ncid, isn, ISNOFF=0, IZBIN, ISTAT, IVAR_IZBIN, IFILE, ifile ;
  int  OPTMASK_MATCH=1; //match CID_IDSURVEY
  double DVAL;
  char id_name[20], CCID[MXCHAR_CCID], CVAL[12], VARLIST_STORE[40]="" ;
  char fnam[] = "parse_cidFile_data" ;

  int use_izbin = INPUTS.izbin_from_cidFile ;
  INPUTS.LEGACY_IZBIN = false;

  // ------------- BEGIN ------------

  // NOTE: OPTMASK_MATCH->0 in this function if IDSURVEY column doesnt exist
  if ( use_izbin ) { sprintf(VARLIST_STORE,"%s", VARNAME_IZBIN);  }

  ncid = match_cidlist_init(fileName, &OPTMASK_MATCH, VARLIST_STORE); 

  // set logical if IZBIN was found.
  if ( use_izbin ) {
    IVAR_IZBIN = IVAR_VARNAME_AUTOSTORE(VARNAME_IZBIN);
    if ( IVAR_IZBIN >=0 ) { INFO_DATA.USE_IZBIN_from_CIDFILE = true; }
  }

  // xxxx prepare to delete LEGACY_IZBIN ...
  if ( INPUTS.LEGACY_IZBIN && use_izbin ) {
    IVAR_IZBIN = IVAR_VARNAME_AUTOSTORE(VARNAME_IZBIN);
    if ( IVAR_IZBIN >= 0 ) { 
      INFO_DATA.IZBIN_from_CIDFILE = (int*) malloc( ncid * sizeof(int) ) ;
      INFO_DATA.USE_IZBIN_from_CIDFILE = true;
      IFILE = NFILE_AUTOSTORE-1;
      for(ifile=0; ifile < NFILE_AUTOSTORE-1; ifile++ ) 
	{ ISNOFF += SNTABLE_AUTOSTORE[ifile].NROW; }
      for ( isn=0; isn < ncid; isn++ ) {
	DVAL = SNTABLE_AUTOSTORE[IFILE].DVAL[IVAR_IZBIN][isn];
	IZBIN = (int)DVAL ;
	INFO_DATA.IZBIN_from_CIDFILE[ISNOFF+isn] = IZBIN ;
      }    
    } // end IVAR_IZBIN>0
  } // end use_izbin
  // xxxx end prep delete  xxxxxxxxxxxxxxx


  // - - - - - - - - 
  // D.Brout Jun 2021
  if ( (OPTMASK_MATCH & 1) == 0) {
    INPUTS.match_on_cid_idsurvey = false;
    INPUTS.match_on_cid_only     = true;
    sprintf(id_name,"CID");
  }
  else {
    INPUTS.match_on_cid_idsurvey = true;
    INPUTS.match_on_cid_only     = false;
    sprintf(id_name,"CID_IDSURVEY");
  }

  INPUTS.ncidList_data += ncid ;
 
  if ( OPT > 0 ) {
    printf("  %s: Accept only %d %s in %s\n", 
	   fnam, ncid, id_name, fileName);
  }
  else {
    printf("  %s: Reject %d  %s in %s\n", 
	   fnam, ncid, id_name, fileName);
  }
  fflush(stdout);

  return ;

} // end parse_cidFile_data

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
  char survey[60], tmpString[MXCHAR_VARNAME] ;
  char fnam[] = "prep_input_nmax" ;

  // ------------- BEGIN ---------------

  if ( IGNOREFILE(item) ) { return ; }

  for(i=0; i < MXARG_nmax; i++ ) {  ptrArg[i] = stringArg[i]; }

  splitString(item, COMMA, MXARG_nmax,    // inputs
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
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
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
  //  char fnam[] = "parse_powzbin" ;

  // ------------- BEGIN ---------------

  INPUTS.powzbin =  0.0 ;  
  INPUTS.znhalf  = -9.0 ;

  splitString(item, COMMA, MXARG,    // inputs
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
void parse_chi2max(char *item) {

  // Created Dec 8 2020
  // Parse strings of the form
  //
  //  chi2max=16             # apply cut to all events
  //  chi2max(DES,PS1)=12    # apply cut only to DES & PS1
  //  chi2max(CSP)=10        # apply cut to CSP
  //
  //  chi2max(CSP,FITWGT0)=10   # set fit wgt = 0 for CSP, but do NOT cut
  //
  //  Note that multiple chi2max inputs are allowed.
  //  In the above example, chi2max=12 is applied to DES & PS1
  //  events; chi2max=10 is applied to CSP; chi2max=16 is applied 
  //  to all other events (e.g., SDSS, low-z that are not CSP, etc...).
  //
  // Jan 2021: check for FITWGT0 option (not a survey name)
  // Mar 4 2021: few bug fixes

  int  IFLAG_GLOBAL = 4;
  int  IFLAG_SURVEY = 8;
  int  debug_malloc = INPUTS.debug_malloc ;

  int  lenkey = strlen("chi2max=");
  int  n_survey, idsurvey, i ;
  double chi2max;
  char string[80], stringOpt[80], **survey_list, *survey ;  
  char fnam[] = "parse_chi2max" ;
  bool SET_FITWGT0, NO_STRINGOPT;
  bool SURVEY_CUTS = ( (INPUTS.iflag_chi2max & IFLAG_SURVEY ) > 0 );
  bool GLOBAL_CUTS = ( (INPUTS.iflag_chi2max & IFLAG_GLOBAL ) > 0 );
  int  LDMP = 0 ;

  // -------------- BEGIN -------------

  // check for argument in ()
  sprintf(string,"%s", item);
  extractStringOpt(string, stringOpt); // return stringOpt
  sscanf( &string[lenkey], "%le", &chi2max );

  if ( LDMP ) {
    printf(" xxx %s: item=%s string=%s opt='%s'   chi2max=%.1f \n",
	   fnam, item, string, stringOpt, chi2max);
    fflush(stdout);
  }

  NO_STRINGOPT = ( strlen(stringOpt) == 0 ) ;
  SET_FITWGT0  = ( strcmp(stringOpt,STRING_FITWGT0) == 0 );

  if ( SET_FITWGT0 ) 
    { INPUTS.iflag_chi2max |= DOFLAG_CUTWIN_FITWGT0 ; }
  else
    { INPUTS.iflag_chi2max |= DOFLAG_CUTWIN_APPLY ; }

  // check trivial case with no argument -> global cut
  if ( NO_STRINGOPT || SET_FITWGT0 ) {
    INPUTS.iflag_chi2max |= IFLAG_GLOBAL ; 
    INPUTS.chi2max = chi2max ;
    if ( SURVEY_CUTS  ) { 
      sprintf(c1err,"Cannot define chi2max after chi2max(SURVEYLIST)");
      sprintf(c2err,"Remove chi2max or define it BEFORE chi2max(SURVEY)");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
    return ;
  }

  // - - - - -

  // allocate SURVEY-dependent chi2max and init all surveys to 
  // global INPUTS.chi2max value
  if ( !SURVEY_CUTS ) {
    print_debug_malloc(+1*debug_malloc,fnam);
    INPUTS.iflag_chi2max |= IFLAG_SURVEY ;
    int MEMD = MXIDSURVEY * sizeof(double) ;
    INPUTS.chi2max_list = (double*)malloc(MEMD);
    for(i=0; i < MXIDSURVEY; i++ ) 
      { INPUTS.chi2max_list[i] = INPUTS.chi2max; }

    // remove global cut flag if we get here (Mar 4 2021)
    if ( GLOBAL_CUTS ) { INPUTS.iflag_chi2max -= IFLAG_GLOBAL; }
  }

  // parse comma-sep list of survey names;
  parse_commaSepList("SURVEYLIST(chi2max)", stringOpt, 
		     MXIDSURVEY, MXCHAR_SAMPLE, &n_survey, &survey_list);

  for(i=0; i < n_survey; i++ ) {
    survey   = survey_list[i] ;
    if ( strcmp(survey,STRING_FITWGT0) == 0 ) {
      sprintf(c1err,"Cannot mix %s option with SURVEY", STRING_FITWGT0);
      sprintf(c2err,"Either remove %s or remove survey names", 
	      STRING_FITWGT0 );
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
    idsurvey = get_IDSURVEY(survey);
    INPUTS.chi2max_list[idsurvey] = chi2max ; 
    printf("\t %s: store cut chi2max=%.1f for SURVEY=%s (%d)\n",
	   fnam, chi2max, survey, idsurvey ); fflush(stdout);
   
  }

  return ;

} // end parse_chi2max

// ************************************************
void parse_blindpar(char *item) {

  // Created Aug 2017
  // Parse 2 parameters separate by comma.
  // E.g., blindpar9=.15,343.2 --> blind offset is 0.15*cos(343.2)

  int  ipar=-9, NARG, MXARG=3;
  char stringArg[2][20], item_local[60] ;
  char *ptrArg[2] = { stringArg[0], stringArg[1] } ;
  char fnam[] = "parse_blindpar" ;

  // ------------- BEGIN ---------------

  if (!strncmp(item,"blindpar9=",10)) 
    { sprintf(item_local, "%s", &item[10] ); ipar=9; }
  else if (!strncmp(item,"blindpar11=",11)) 
    { sprintf(item_local, "%s", &item[11] ); ipar=11; }
  else {
    sprintf(c1err,"Invalid input '%s'", item);
    sprintf(c2err,"Check valid blindpar keys");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  splitString(item_local, COMMA, MXARG,    // inputs
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
void parse_prescale_biascor (char *item, int wrflag) {

  // June 25, 2016
  // parse input *item that is a comma-separated prescale of the form
  //
  //   prescale_simbias=3,4
  //
  //   --> prescale by 4 and use 3rd subset.
  //   --> INPUTS.prescale_simbias[0] = 3
  //   --> INPUTS.prescale_simbias[1] = 4
  //
  //  Note that PS0 can have values 0,1,2 ... PS1-1
  //
  // wrflag = 1 -> set flag to write pre-scaled biasCor to out file.

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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  return ;

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
  char fnam[] = "parse_IDSAMPLE_SELECT" ;

  // --------- BEGIN -----------

  if ( strlen(item) == 0 ) { return; }

  sprintf(itemLocal, "%s", item );
  for(i=0; i < MXNUM_SAMPLE; i++ ) { ptrID[i] = strID[i]; }

  if ( strstr(itemLocal,",") != NULL ) {
    sprintf(c1err,"Illegal comma delimeter in idsample_select key");
    sprintf(c2err,"Try plus (+) delimeter instead.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  splitString(itemLocal, PLUS, MXNUM_SAMPLE,      // inputs
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

  fprintf(FP_STDOUT, " Will select %d IDSAMPLEs : %s \n",  NTMP,  itemLocal);
  fflush(FP_STDOUT);

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
  // Mar 17 2022; sigmB -> sigint_fix (not 0)

  int  NSAMPLE = NSAMPLE_BIASCOR ;
  int  idsample, Nsigint, i ;
  double sigint;
  char *name ;
  char fnam[] = "parse_sigint_fix";
  char itemLocal[200], *ptrSIG[MXNUM_SAMPLE], strSIG[MXNUM_SAMPLE][8] ;

  // --------- BEGIN -----------

  if ( strlen(item) == 0 ) { return ;}

  // check for comma
  if ( strstr(item,COMMA) == NULL ) {
    // no comma --> fix same sigint for all IDSAMPLE
    sscanf(item, "%le", &sigint) ;
    for(idsample=0; idsample<NSAMPLE; idsample++ ) 
      { SNDATA_INFO.sigint_fix[idsample] = sigint; }
    INPUTS.sigmB = sigint; 
  }
  else {
    // strip sigint for each IDSAMPLE
    for(i=0; i < MXNUM_SAMPLE; i++ ) { ptrSIG[i] = strSIG[i]; }
    sprintf(itemLocal,"%s", item);
    splitString(itemLocal, COMMA, MXNUM_SAMPLE,      // inputs
		&Nsigint, ptrSIG );                    // outputs
    if ( Nsigint != NSAMPLE ) {
      sprintf(c1err,"Nsiginit=%d != N_IDSAMPLE=%d", Nsigint, NSAMPLE);
      sprintf(c2err,"sigint_fix=%s must include each IDSAMPLE",	 item) ;
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }

    for(idsample=0; idsample < Nsigint; idsample++ ) {
      sscanf(strSIG[idsample], "%le", &sigint) ;
      SNDATA_INFO.sigint_fix[idsample] = sigint; 
    }
    INPUTS.sigmB = SNDATA_INFO.sigint_fix[0]; 
    
  }


  // - - - - -
  // turn off other sigint inputs
  // xxx mark delete  INPUTS.sigmB = 0.0 ;  
  
  // ----------
  // print sigint_fix for each IDSAMPLE, and compute sqsigint_fix
  
  fprintf(FP_STDOUT, "\n");
  for(idsample=0; idsample < NSAMPLE; idsample++ ) {
    name   = SAMPLE_BIASCOR[idsample].NAME ;
    sigint = SNDATA_INFO.sigint_fix[idsample];
    SNDATA_INFO.sqsigint_fix[idsample] = sigint*sigint ;

    fprintf(FP_STDOUT, "  sigint_fix=%6.4f for IDSAMPLE=%d (%s)\n",
	   sigint, idsample, name);
  }
  fflush(FP_STDOUT);
    
  
  return ;

} // end parse_sigint_fix

// **************************************************
void parse_CUTWIN(char *line_CUTWIN) {

  // Created 4/24/2012 by R.Kessler
  // parse *line_CUTWIN with four space-separated items:
  //   CUTWIN <VARNAME>  <MIN> <MAX>
  //     or
  //   CUTWIN([stringOpt]) <VARNAME>  <MIN> <MAX>
  //
  // where stringOpt can be: NOABORT DATAONLY BIASCORONLY FITWGT0
  //
  // Note that CUTWIN can inlcude multiple options, e.g., 
  //    CUTWIN(OPT1,OPT2,..ETC)
  //
  // Fill following globals:
  //   INPUTS.NCUTWIN
  //   INPUTS.CUTWIN_NAME
  //   INPUTS.CUTWIN_RANGE
  //
  //   INPUTS.LCUTWIN_DATAONLY
  //   INPUTS.LCUTWIN_BIASCORONLY
  //   INPUTS.LCUTWIN_ABORTFLAG   ! Sep 2016
  //   INPUTS.LCUTWIN_FITWGT0     ! Jan 2021
  //
  //
  // Oct 8 2020: check for bad input
  // Jan 19 2021: 
  //   + check FITWGT0 option
  //   + enable comma-sep list of options
  //   + abort on invalid option
  //

  int   NITEM = 4;
  char  item_list[4][60], line_local[200], string[60], KEY[60] ;
  char  *item, *ptrtok, *cutwinOpt, **cutwinOpt_list ;

  int   ICUT, i, opt, nread, NOPT ;
  char  fnam[] = "parse_CUTWIN" ;

  // ---------- BEGIN ------------

  ICUT = INPUTS.NCUTWIN ;

  INPUTS.CUTWIN_NAME[ICUT][0]      = 0 ;
  INPUTS.CUTWIN_RANGE[ICUT][0]     = -99999.0 ;
  INPUTS.CUTWIN_RANGE[ICUT][1]     = -99999.0 ;

  INPUTS.LCUTWIN_ABORTFLAG[ICUT]   = true ;   //  abort on missing var
  INPUTS.LCUTWIN_DATAONLY[ICUT]    = false ;  //  cut on data 
  INPUTS.LCUTWIN_BIASCORONLY[ICUT] = false ;  //  cut on sim data and biasCor
  INPUTS.LCUTWIN_FITWGT0[ICUT]     = false ;  //  MUERR->888 instead of cut

  // - - - - - -  -
  // strip each line_CUTWIN item into item_list, and check for missing items
  sprintf(line_local,"%s", line_CUTWIN);
  ptrtok = strtok(line_local," ");
  for ( i=0; i < 4 ; i++ ) {
    sprintf(item_list[i], "%s", ptrtok);

    //    printf(" xxx %s: i=%d item='%s' \n", fnam, i, item_list[i]);
    // check disable option with "CUTWIN NONE" 
    if ( i==1 && strcmp(item_list[i],"NONE") == 0 )  {  
      INPUTS.LCUTWIN_DISABLE = true ; 
      printf("\n\t DISABLE CUTS EXCEPT CIDLIST,BIASCOR,z,BADCOV\n"); 
      fflush(stdout); 
      return;
    }

    ptrtok = strtok(NULL, " ");
    // Oct 8 2020: check for missing CUTWIN element
    if ( i < 3 && ptrtok == NULL ) {
      sprintf(c1err,"Problem reading CUTWIN element i=%d", i+1 );
      sprintf(c2err,"for line_CUTWIN = '%s' ", line_CUTWIN);
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
    }
  }

  INPUTS.NCUTWIN++ ;

  // - - - - - 
  for ( i=0; i < NITEM ; i++ ) {

    item = item_list[i];

    if ( i == 0 ) {
      // check for option in CUTWIN(string)
      sscanf ( item, "%s", KEY ); 
      extractStringOpt(KEY, string); // return string

      // return list of args separated by commas.
      // E.g., if string = 'DATAONLY,FITWGT0' then
      // stringOpt_list = 'DATAONLY', 'FITWGT0'
      parse_commaSepList("CUTWIN_OPTION", string,
			 NITEM, 40,   // Max number of options and strlen 
			 &NOPT, &cutwinOpt_list); // <== returned
      
      //      printf(" xxx %s: string='%s' -> NOPT=%d \n",
      //	     fnam, string, NOPT); fflush(stdout);

      for ( opt=0; opt < NOPT; opt++ ) {
	cutwinOpt = cutwinOpt_list[opt];

	if ( strcmp(cutwinOpt,"NOABORT") == 0 ) 
	  { INPUTS.LCUTWIN_ABORTFLAG[ICUT] = false; } // allow missing var 

	else if ( strcmp(cutwinOpt,"DATAONLY") == 0 ) 
	  { INPUTS.LCUTWIN_DATAONLY[ICUT] = true ; } // cut on data only
	
	else if ( strcmp(cutwinOpt,"BIASCORONLY") == 0 ) 
	  { INPUTS.LCUTWIN_BIASCORONLY[ICUT] = true ; } // cut on sim & biascor

	else if ( strcmp(cutwinOpt,STRING_FITWGT0) == 0 ) 
	  { INPUTS.LCUTWIN_FITWGT0[ICUT] = true ; } 
      
	else {
	  sprintf(c1err,"Invalid CUTWIN option: '%s'", cutwinOpt);
	  sprintf(c2err,"Valid options: NOABORT DATAONLY BIASCORONLY FITWGT0");
	  errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
	}

      } // end opt loop

    } // end i= loop over CUTWIN(option) item


    if ( i == 1 ) { 
      nread = sscanf ( item, "%s", INPUTS.CUTWIN_NAME[ICUT] ); 
      if ( nread != 1 ) { abort_bad_input(KEY, ptrtok, i, fnam); }
    } 
        
    if ( i == 2 ) {
      nread = sscanf (item, "%le", &INPUTS.CUTWIN_RANGE[ICUT][0] ); 
      if ( nread != 1 ) { abort_bad_input(KEY, ptrtok, i, fnam); }
    } 

    if ( i == 3 ) {
      nread = sscanf (item, "%le", &INPUTS.CUTWIN_RANGE[ICUT][1] ); 
      if ( nread != 1 ) { abort_bad_input(KEY, ptrtok, i, fnam); }
    } 

  } // end i loop


  // - - - - -
  // 9.15.2021: allow command line override of CUTWIN with looser
  //   cut by repacing previous CUTWIN with same variable name.
  int icut;
  char *name, *NAME = INPUTS.CUTWIN_NAME[ICUT];
  for(icut=0; icut < ICUT; icut++ ) {
    name = INPUTS.CUTWIN_NAME[icut] ;
    if ( strcmp(NAME,name) == 0 ) {
      copy_CUTWIN(ICUT,icut);
      INPUTS.NCUTWIN-- ;
    }
  }

  // - - - - - -
  char cMUERR[20] = "" ;

  if ( INPUTS.LCUTWIN_FITWGT0[ICUT] )
    { sprintf(cMUERR,"MUERR->%.1f", MUERR_FITWGT0); }

  fprintf(FP_STDOUT, 
	 "\t Apply CUTWIN on %12s from %10.4f to %10.4f "
	  " (ABORTFLAG=%d,%s)\n"
	 ,INPUTS.CUTWIN_NAME[ICUT]
	 ,INPUTS.CUTWIN_RANGE[ICUT][0]
	 ,INPUTS.CUTWIN_RANGE[ICUT][1]
	 ,INPUTS.LCUTWIN_ABORTFLAG[ICUT]
	 ,cMUERR
	  ) ;

  return ;
} // end of parse_CUTWIN

// **************************************************
void copy_CUTWIN(int icut0,int icut1) {

  // Created Sep 15 2021
  // Copy CUTWIN contents from icut0 to icut1.
  // Used to enable command line CUTWIN overrides with relaxed cut
  // compared to CUTWIN in the input file.

  char fnam[] = "copy_CUTWIN" ;

  // ---------- BEGIN ---------

  sprintf(INPUTS.CUTWIN_NAME[icut1], "%s", INPUTS.CUTWIN_NAME[icut0] );
  INPUTS.CUTWIN_RANGE[icut1][0] = INPUTS.CUTWIN_RANGE[icut0][0];
  INPUTS.CUTWIN_RANGE[icut1][1] = INPUTS.CUTWIN_RANGE[icut0][1];

  INPUTS.LCUTWIN_ABORTFLAG[icut1]   = INPUTS.LCUTWIN_ABORTFLAG[icut0] ;
  INPUTS.LCUTWIN_DATAONLY[icut1]    = INPUTS.LCUTWIN_DATAONLY[icut1] ;
  INPUTS.LCUTWIN_BIASCORONLY[icut1] = INPUTS.LCUTWIN_BIASCORONLY[icut1] ;
  INPUTS.LCUTWIN_FITWGT0[icut1]     = INPUTS.LCUTWIN_FITWGT0[icut1] ;

  return;
} // end copy_CUTWIN

// **************************************************
int reject_CUTWIN(int EVENT_TYPE, int *DOFLAG_CUTWIN, double *CUTVAL_LIST) {

  // Created Jan 2016 [major refactor Jan 2021]
  //
  // Returns  0 if all CUTVAL_LIST values pass CUTWIN cuts.
  // Returns  1 if any CUTVAL cut fail CUTWIN
  // Returns  2 if CUTVAL that fails has FITWGT0 option (data only)
  //
  // Input EVENT_TYPE specifies DATA or BIASCOR
  //
  // Input array *DOFLAG_CUTWIN are instructions :
  //  DOFLAG_CUTINW[icut] = 0 -> do not apply cut (ignore)
  //  DOFLAG_CUTINW[icut] = 1 -> apply explicit cut
  //  DOFLAG_CUTINW[icut] = 2 -> no cut, but deweight with MUERR->large val
  //
  // CUTVAL_LIST is array of values to apply CUTWIN
  //
  // Jan 22 2021: 
  //   rename apply_CUTWIN -> reject_CUTWIN, and return reject flag
  //   with similar meaning as DOFLAG_CUTWIN.
  
  int IS_DATA    = ( EVENT_TYPE == EVENT_TYPE_DATA );
  int IS_BIASCOR = ( EVENT_TYPE == EVENT_TYPE_BIASCOR );
  int IS_CCPRIOR = ( EVENT_TYPE == EVENT_TYPE_CCPRIOR );

  int LDMP = 0 ; // (OPT==666);
  int icut, reject, DOFLAG ;
  double CUTVAL, *CUTWIN ;
  char *NAME;
  char fnam[] = "reject_CUTWIN" ;

  
  // ------------- BEGIN -----------

  if ( LDMP ) { printf(" xxx --------------------------- \n"); }

  reject = 0 ;    // init to pass cuts

  for(icut=0; icut < INPUTS.NCUTWIN; icut++ ) {
   
    DOFLAG = DOFLAG_CUTWIN[icut] ;
    if ( DOFLAG == DOFLAG_CUTWIN_IGNORE ) { continue; }

    // if not data (i.e.,  biasCor & CCprior), DOFLAG must be
    // "APPLY". The FITWGT0 feature works only for data.
    if ( !IS_DATA ) { DOFLAG = DOFLAG_CUTWIN_APPLY; }

    CUTVAL = CUTVAL_LIST[icut];
    CUTWIN = &INPUTS.CUTWIN_RANGE[icut][0];
    NAME   = INPUTS.CUTWIN_NAME[icut] ;

    if ( LDMP ) {
      printf(" xxx cut on %s = %f  (cutwin=%.3f to %.3f, EVENT_TYPE=%d)\n",
	     NAME, CUTVAL, CUTWIN[0], CUTWIN[1], EVENT_TYPE ); 
    }

    // check SIM-option to skip SIM_TYPE_INDEX
    if( (IS_BIASCOR || IS_CCPRIOR) && !usesim_CUTWIN(NAME)  ) 
      { continue ; }

    if ( CUTVAL < CUTWIN[0] || CUTVAL > CUTWIN[1] ) { 
      // be careful here; once reject is set to DOFLAG_CUTWIN_APPLY,
      // the reject flag cannot be changed (e.g, to DOFLAG_CUTWIN_FITWGT0)
      if ( reject != DOFLAG_CUTWIN_APPLY )  {  reject = DOFLAG ; }
    }
  }

  // - - - - 
  if (LDMP ) { 
    printf(" xxx ---> reject = %d\n", reject ); 
    fflush(stdout); debugexit(fnam);
  } 

  return(reject);


} // end reject_CUTWIN

// **************************************************
int set_DOFLAG_CUTWIN(int ivar, int icut, int isData) {

  // Return 1 if ivar >= 0 -> apply explicit cut.
  // Return 2 to set MUERR = large (implicit cut by deweight) [Jan 2021]
  // Return 0 to ignore this cut.
  //
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
  // Oct 28 2020: new ISVAR_PROB logic
  // Jan 21 2021: minor refac using DOFLAG_CUTWIN_XXX params
  // May 22 2021: check L_DISABLE

  bool  L_VALID_VAR   = ( ivar >= 0 );
  bool  L_NOVAR       = !L_VALID_VAR ;
  bool  L_ABORTFLAG   = INPUTS.LCUTWIN_ABORTFLAG[icut];
  bool  L_DATAONLY    = INPUTS.LCUTWIN_DATAONLY[icut];
  bool  L_BIASCORONLY = INPUTS.LCUTWIN_BIASCORONLY[icut];
  bool  L_FITWGT0     = INPUTS.LCUTWIN_FITWGT0[icut];
  bool  L_DISABLE     = INPUTS.LCUTWIN_DISABLE ;
  char *VARNAME     = INPUTS.CUTWIN_NAME[icut];
  bool  ISVAR_PROB  = (strstr(VARNAME,"PROB_") != NULL ); // Oct 2020
  int   DOFLAG ;
  char  DATATYPE[12]; // DATA or BIASCOR
  char  fnam[] = "set_DOFLAG_CUTWIN" ;

  // ------------- BEGIN -------------
  
  if ( L_DATAONLY    && !isData ) { return(DOFLAG_CUTWIN_IGNORE); }
  if ( L_BIASCORONLY &&  isData ) { return(DOFLAG_CUTWIN_IGNORE); }

  // Oct 28 2020: Apply cut to biasCor sample if varname doesn't exist
  //    and starts with PROB. This assumes that idsurvey_list_probcc0
  //    sets pIa=1 ... if not, all these events will be rejected. 
  //    This logic is not needed for data because the data-catenate 
  //    process ensures existing PROB columns.
  if ( !isData && ivar < 0 && ISVAR_PROB ) { return(DOFLAG_CUTWIN_APPLY); }

  if ( L_NOVAR && L_ABORTFLAG ) {
    if ( isData ) { sprintf(DATATYPE,"DATA"); }
    else          { sprintf(DATATYPE,"BIASCOR"); }

    sprintf(c1err,"Invalid CUTWIN on '%s' for %s (ivar=%d, icut=%d)", 
	    VARNAME, DATATYPE, ivar, icut );
    sprintf(c2err,"Check CUTWIN keys in input file" ); 
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( L_VALID_VAR ) { 
    DOFLAG = DOFLAG_CUTWIN_APPLY ;    
    if (L_FITWGT0)  { DOFLAG = DOFLAG_CUTWIN_FITWGT0; }
  }
  else
    { DOFLAG = DOFLAG_CUTWIN_IGNORE ; }

  return(DOFLAG);

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
  int debug_malloc = INPUTS.debug_malloc ;
  int LDMP = 0 ;
  char fnam[] = "parse_FIELDLIST" ;

  // ------------ BEGIN ------------
  
  print_debug_malloc(+1*debug_malloc,fnam);
  for(i=0; i < MXFIELD_OVERLAP; i++ ) 
    { INPUTS.FIELDLIST[i] = (char*) malloc(20*sizeof(char) ); }

  splitString(item, COMMA, MXFIELD_OVERLAP,               // inputs
	      &INPUTS.NFIELD, INPUTS.FIELDLIST ); // outputs
  
  if ( LDMP ) {
    for(i=0; i < INPUTS.NFIELD; i++ ) {
      printf(" xxx %s: select FIELD = '%s' \n", 
	     fnam, INPUTS.FIELDLIST[i] ); fflush(stdout);
    }
  }

  return ;

} // end parse_FIELDLIST

// **************************************************
void prep_input_repeat(void) {

  // Call this prep function before repeating a fit.
  // Examples:
  //   + NSPLITRAN>1
  //   + SUBPROCESS called from python fitter
  //

  char fnam[] = "prep_input_repeat" ;

  // ------------ BEGIN -----------

  FITINP.COVINT_PARAM_FIX   = INPUTS.sigmB ; 
  FITINP.COVINT_PARAM_LAST  = INPUTS.sigmB ; 
  INPUTS.covint_param_step1 = INPUTS.sigint_step1 ; // default COVINT param

  if ( strlen(INPUTS.sigint_fix) > 0 ) {
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ; 
  }

  recalc_dataCov(); 

  // - - - - - - -

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) {
    //printf("   Reset a few CUTBITS \n" );
    int isn, CUTMASK ;
    int NSN_DATA = INFO_DATA.TABLEVAR.NSN_ALL ; 
    for (isn=0; isn < NSN_DATA; isn++)  {

      CUTMASK = INFO_DATA.TABLEVAR.CUTMASK[isn];
      
      // reset... cut bits that get re-applied in SUBPROCESS
      CUTMASK -= ( CUTMASK & CUTMASK_LIST[CUTBIT_SPLITRAN] ) ;
      //      CUTMASK -= ( CUTMASK & CUTMASK_LIST[CUTBIT_MINBIN]   ) ;
      if ( SUBPROCESS.NEVT_SIM_PRESCALE > 0 )
	{ CUTMASK -= (CUTMASK & CUTMASK_LIST[CUTBIT_SIMPS] ); }
      
      INFO_DATA.TABLEVAR.CUTMASK[isn] = CUTMASK ;
      
    } // end isn
  }
#endif

  return;

} // end prep_input_repeat


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


// ======================================
void  CPU_SUMMARY(void) {

  // Created Nov 22 2017
  fprintf(FP_STDOUT, "\n PROCESS-TIME SUMMARY: \n");

  fprintf(FP_STDOUT, "   Ini-time:  %.2f minutes \n",
	 (t_end_init-t_start)/60.0 ) ;

  if ( INPUTS.opt_biasCor > 0 ) {
    fprintf(FP_STDOUT, "\t BiasCor read time %.2f minutes \n",
	   (t_read_biasCor[1]- t_read_biasCor[0])/60.0 );
  }

  fprintf(FP_STDOUT, "   Fit-time:  %.2f minutes \n",
	 (t_end_fit-t_end_init)/60.0 ) ;

  fprintf(FP_STDOUT, "\n");
  fflush(FP_STDOUT);

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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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

  int ICUTWIN_VARNAME_PIA = INFO_DATA.TABLEVAR.ICUTWIN_VARNAME_PIA;
  char *varname_pIa       = INPUTS.varname_pIa ;

  int i, icut;
  int  NFITPAR, ifile, NTMP=0, USE_CCPRIOR, USE_CCPRIOR_H11, OPT ;
  int  OPTMASK;
  char usage[10], *tmpName ;
  char fnam[] = "prep_input_driver";

  // ------------ BEGIN -----------

  prep_debug_flag();

  if ( INPUTS.cat_only ) { prep_input_varname_missing(); return; }

  // check option to select only true SNIA and ignore CC prior
  if ( INPUTS.select_trueIa ) { prep_input_trueIa(); }

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
    fprintf(FP_STDOUT, "\n FIXPAR_ALL option -> "
	   "Compute MU from initial values (no floated params)\n\n" );
    fflush(FP_STDOUT);
	   
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( INPUTS.nlogzbin > 0 ) { INPUTS.nzbin = INPUTS.nlogzbin ; }

  NFITPAR = INPUTS.nzbin + MXCOSPAR; // estimate
  if ( NFITPAR >= MAXPAR_MINUIT ) {
    sprintf(c1err,"NFITPAR=%d (%d z bins + %d params)",
	    NFITPAR, INPUTS.nzbin, MXCOSPAR );
    sprintf(c2err,"But MAXPAR_MINUIT = %d", MAXPAR_MINUIT );
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  NSIMIa = NSIMCC = NSIMDATA = 0 ;
  NJOB_SPLITRAN = 0;
  
  FITINP.COVINT_PARAM_FIX   = INPUTS.sigmB ;
  FITINP.COVINT_PARAM_LAST  = INPUTS.sigmB ; 
  FITINP.DCHI2RED_DSIGINT   = INPUTS.dchi2red_dsigint; // Mar 25 2022
  INPUTS.covint_param_step1 = INPUTS.sigint_step1 ; // default COVINT param

  // Oct 9 2018: for sigInt(IDSAMPLE), vary global COV scale instead
  //             of varying sigInt.
  if ( strlen(INPUTS.sigint_fix) > 0 ) {
    sprintf(FITPARNAMES_DEFAULT[IPAR_COVINT_PARAM], "scale_covint"); 
    FITINP.COVINT_PARAM_FIX    = 1.0 ; 
    FITINP.COVINT_PARAM_LAST   = 1.0 ;
    INPUTS.covint_param_step1  = INPUTS.scale_covint_step1 ; 
  }
    
  fprintf(FP_STDOUT, "zmin =%f zmax =%f bins=%i \n",
	 INPUTS.zmin, INPUTS.zmax, INPUTS.nzbin);

  if ( INPUTS.nzbin >= MXz ) {
    sprintf(c1err,"nzbin=%d exceeds bound of MAXBIN_z=%d", 
	    INPUTS.nzbin, MXz);
    sprintf(c2err,"Check 'bins=' key in input file.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);   
  }

  if (INPUTS.Nsntype <= 0 ) 
    { fprintf(FP_STDOUT, "SN selection type=ALL\n"); }
  else 
    { fprintf(FP_STDOUT, "SN selection sntype= %s\n", INPUTS.sntypeString); }


  fprintf(FP_STDOUT, "x1min/x1max = %.3f / %.3f \n", 
	  INPUTS.x1min, INPUTS.x1max);
  fprintf(FP_STDOUT, "cmin/cmax   = %.3f / %.3f \n", 
	 INPUTS.cmin, INPUTS.cmax);
  fprintf(FP_STDOUT, "logmass_min/max = %.3f/%.3f \n", 
	 INPUTS.logmass_min, INPUTS.logmass_max);
  fprintf(FP_STDOUT, "H0=%f \n", INPUTS.H0);
  fprintf(FP_STDOUT, "Nominal M0=%f \n", INPUTS.nommag0);
  fprintf(FP_STDOUT, "zpecerr(override)   = %.5f \n", INPUTS.zpecerr );
  fprintf(FP_STDOUT, "dsigint/dz(lensing) = %.4f \n", INPUTS.lensing_zpar );

  if (INPUTS.uave) 
    { fprintf(FP_STDOUT, "Use average M0.\n"); }
  else 
    { fprintf(FP_STDOUT, "Use nominal M0.\n"); }
  if ( INPUTS.fitflag_sigmb )  { 
    fprintf(FP_STDOUT, "Will fit for sig1 such that chi2/N ~ 1.0\n");
    fprintf(FP_STDOUT, "Will use sig1tol=%f . \n", INPUTS.redchi2_tol);
    fprintf(FP_STDOUT, "Will allow max sigmb iterations=%d .\n", MAXFITITER); 
  }
  

  fprintf(FP_STDOUT, "uM0=%d --> ", INPUTS.uM0 ) ;
  if ( INPUTS.uM0 == M0FITFLAG_CONSTANT ) 
    { fprintf(FP_STDOUT, "M0 fixed. \n"); }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_FLAT ) 
    { fprintf(FP_STDOUT, "M0 floated in each z bin. \n"); }
  else if ( INPUTS.uM0 == M0FITFLAG_ZBINS_INTERP ) 
    { fprintf(FP_STDOUT, "M0 floated in each z bin,and interpolated. \n"); }
  else
    { fprintf(FP_STDOUT, "M0 option not defined ???\n"); }

  if ( INPUTS.uzsim ) { fprintf(FP_STDOUT, "REDSHIFT CHEAT: z -> simz \n");  }

  fflush(FP_STDOUT);

  if ( INPUTS.prescale_simData > 1.0 ) {
    fprintf(FP_STDOUT, "PRE-SCALE SIMDATA by %.1f \n", 
	    INPUTS.prescale_simData);
  }
  if ( INPUTS.prescale_simIa > 1.0 ) {
    fprintf(FP_STDOUT, "PRE-SCALE SIMIa by %.1f \n", 
	    INPUTS.prescale_simIa);
  }
  if ( INPUTS.prescale_simCC > 1.0 ) {
    fprintf(FP_STDOUT, "PRE-SCALE SIMCC by %.1f \n", 
	    INPUTS.prescale_simCC);
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }
  


  // if CCprior is not set then make sure to fix scalePCC and sigint
  if ( !USE_CCPRIOR  ) {
    if ( strlen(varname_pIa) > 0 ) {
      sprintf(c1err,"Illegal varname_pIa=%s without CC prior.", varname_pIa);
      sprintf(c2err,"Must set simfile_ccprior with varname_pIa");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
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
      fprintf(FP_STDOUT, "\t WARNING: scalePCC(u13,p13)=0,0 --> "
	     "turn off CC prior\n");
      fflush(FP_STDOUT);
    }

    INPUTS.nfile_CCprior = 0 ;
    INFO_CCPRIOR.USE = 0; 

  }


  if ( INPUTS.force_pIa >= 0.0 ) 
    { fprintf(FP_STDOUT, " force_pIa = %.3f \n", INPUTS.force_pIa ); }

  if ( INPUTS.perfect_pIa ) 
    { fprintf(FP_STDOUT," force_pIa = 1 or 0 for true SNIa or SNCC \n"); }


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
   
    fprintf(FP_STDOUT, "%s  %s  starting value=%8.3f  (step=%.3f)\n",
	   FITPARNAMES_DEFAULT[i], usage,
	   INPUTS.parval[i], INPUTS.parstep[i] );         
  }


  if (INPUTS.zmax <= INPUTS.zmin) {
    sprintf(c1err, "Invalid z range zmin=%f zmax=%f",
	    INPUTS.zmin,INPUTS.zmax);
    sprintf(c2err,"Check inputs zmin and zmax");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 	
  }

  // prepare stuff related to gamma = HR step
  prep_input_gamma();

  // sanity checks on fitting for GAMMA0 (mag step across host logmass)
  if ( INPUTS.USE_GAMMA0 ) {
    double TAU   = INPUTS.parval[IPAR_LOGMASS_TAU] ;
    if ( TAU < 0.001 ) {
      sprintf(c1err,"Invalid LOGMASS_TAU = %.4f for gammma0 fit", TAU);  
      sprintf(c2err,"logmasss_tau(p8) should be at least .01");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
  }


  // - - - - -
  prep_input_load_COSPAR();

  INPUTS.FLOAT_COSPAR=0;
  if ( INPUTS.ipar[IPAR_OL] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_Ok] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_w0] ) { INPUTS.FLOAT_COSPAR=1; }
  if ( INPUTS.ipar[IPAR_wa] ) { INPUTS.FLOAT_COSPAR=1; }

  prep_input_nmax(INPUTS.nmaxString);

  INFO_DATA.USE_IZBIN_from_CIDFILE = false;
  INFO_DATA.NCHANGE_IZBIN = 0;
  if ( INPUTS.ncidFile_data > 0 ) {
    printf("\n");
    OPT = INPUTS.acceptFlag_cidFile_data;
    OPTMASK = 0;
    if ( INPUTS.izbin_from_cidFile ) { OPTMASK += 64; }

    // init hash table 
    match_cidlist_init(BLANK_STRING, &OPTMASK, BLANK_STRING); 
    for(ifile=0; ifile < INPUTS.ncidFile_data; ifile++ )  { 
      parse_cidFile_data(OPT, INPUTS.cidFile_data[ifile] ); 
    }
    printf("\n"); fflush(stdout);
  }


  prep_input_varname_missing();
  
  fprintf(FP_STDOUT, "\n"); fflush(FP_STDOUT);

  INFO_DATA.TABLEVAR.EVENT_TYPE    = EVENT_TYPE_DATA ;
  INFO_BIASCOR.TABLEVAR.EVENT_TYPE = EVENT_TYPE_BIASCOR ;
  INFO_CCPRIOR.TABLEVAR.EVENT_TYPE = EVENT_TYPE_CCPRIOR ;

  // Sep 2020: thread checks
  int nthread = INPUTS.nthread ;
  if ( nthread > 1 ) {
    if ( nthread >= MXTHREAD ) {
      sprintf(c1err, "nthread=%d exceeds bound of MXTHREAD=%d", 
	      nthread, MXTHREAD );
      sprintf(c2err, "Either reduce nthread or increase MXTHREAD");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
  }


  return ;

} // end of prep_input_driver

// =========================
void  prep_input_trueIa(void) {
  // Created Jun 17 2021
  // Disable true CC and CC prior; fit only SNIa 

  int  NFILE = INPUTS.nfile_data;
  int  ifile;
  char fnam[] = "prep_input_trueIa" ;
  // ---------- BEGIN -------------

  sprintf(BANNER,"%s: select only true SNIa; PROBcc -> 0", fnam);
  fprint_banner(FP_STDOUT,BANNER);

  INPUTS.prescale_simCC = 999999.;
  INPUTS.nfile_CCprior  = 0 ;
  INFO_CCPRIOR.USE = INFO_CCPRIOR.USEH11 = 0;

  INPUTS_PROBCC_ZERO.USE       = false ;
  INPUTS_PROBCC_ZERO.ntype     = 0 ;
  INPUTS_PROBCC_ZERO.nidsurvey = 0 ;  

  INPUTS.varname_pIa[0]        = 0;
  for(ifile=0; ifile < NFILE; ifile++ ) 
    { INFO_DATA.TABLEVAR.IVAR_pIa[ifile] = -9;  }

  return;

} // end prep_input_trueIa

// =========================
void prep_debug_flag(void) {

  // -----------------------------------
  // check INPUTS.debug_flag and set internal LEGACY/REFAC options.
  // These options are intended to be temporary to allow switching
  // back and forth between LEGACY/REFAC codes.
  // General convention: 
  //   debug_flag > 0 -> test unreleased refactor
  //   debug_flag < 0 -> switch back to legacy code 

  if ( INPUTS.restore_mucovscale_bug ) {
    printf("\n RESTORE mucovscale bug \n");
    fflush(stdout);
  }
  if ( INPUTS.restore_mucovadd_bug ) {
    printf("\n RESTORE mucovadd bug (set to %d)\n", INPUTS.restore_mucovadd_bug);
    fflush(stdout);
  }
  
  if ( INPUTS.debug_flag!=0) {
    printf("\n debug flag set to %d\n", INPUTS.debug_flag);
    fflush(stdout);
  }
  fflush(FP_STDOUT);

}  // end prep_debug_flag

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
  //   idsurvey_list_probcc0=DES,52,53 
  //       (list of string and/or ID from SURVEY.DEF file)
  //
  //  The type and idsurvey lists correspond to spec-confirmed
  //  SNIa, and thus PROBCC is forced explicitly to zero.
  //
  char *str_type_list     = INPUTS_PROBCC_ZERO.str_type_list ;
  char *str_idsurvey_list = INPUTS_PROBCC_ZERO.str_idsurvey_list ;
  int debug_malloc = INPUTS.debug_malloc ;
  int  LEN_type_list      = strlen(str_type_list);
  int  LEN_idsurvey_list  = strlen(str_idsurvey_list);
  int  DO_PROBCC0, i, itype, id, nval, NERR=0 ;
  int  NUSE, NUSE_IDSURVEY[MXIDSURVEY];
  char *str_values[MXPROBCC_ZERO], *surveyName ;
  char fnam[]  = "prep_input_probcc0" ;

  // ---------------- BEGIN ------------------
  
  DO_PROBCC0 = ( LEN_type_list > 0 || LEN_idsurvey_list > 0 ) ;
  if ( DO_PROBCC0 ) {
    print_debug_malloc(+1*debug_malloc,fnam);
    INPUTS_PROBCC_ZERO.USE = true;
    for(i=0; i < MXPROBCC_ZERO; i++ ) 
      { str_values[i] = (char*)malloc( 80 * sizeof(char) ); }
  }
  else {
    return ;
  }


  // check TYPE from data header
  if ( LEN_type_list > 0 ) {
    splitString(str_type_list, COMMA, MXPROBCC_ZERO,    // inputs
		&nval, str_values ) ;                    // outputs    
    INPUTS_PROBCC_ZERO.ntype = nval ;
    for(i=0; i < nval; i++ ) {
      sscanf(str_values[i],"%d", &itype);
      INPUTS_PROBCC_ZERO.type_list[i] = itype ;
      fprintf(FP_STDOUT, "\t Force PROB_BBC(CC) = 0 for TYPE = %d\n", itype);
      fflush(FP_STDOUT);
    }
  }
 

  // check survey ID from $SNDATA_ROOT/SURVEY.DEF
  if ( LEN_idsurvey_list > 0 ) {
    splitString(str_idsurvey_list, COMMA, MXPROBCC_ZERO,    // inputs
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

      fprintf(FP_STDOUT, "\t Force PROB_BBC(CC) = 0 for IDSURVEY = %3d (%s)\n", 
	     id, surveyName );

      NUSE=NUSE_IDSURVEY[id];
      if( NUSE > 1 ) {
	fprintf(FP_STDOUT,"\t\t ERROR: IDSURVEY=%d used %d times\n", id, NUSE);
	NERR++ ;
      }

      fflush(FP_STDOUT);
    }
  }

  
  // free string memory
  print_debug_malloc(-1*debug_malloc,fnam);
  for(i=0; i < MXPROBCC_ZERO; i++ ) 
    { free(str_values[i]); }

  // abort on error
  if ( NERR > 0 ) {
    sprintf(c1err,"Found %d errors above.", NERR );
    sprintf(c2err,"check inputs "
	    "type_list_probcc0 and idsurvey_list_probcc0");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  // Oct 29 2020: also replace varname_pIa with actual name
  int icut, FOUND_GAMMA=0;
  int NCUTWIN = INPUTS.NCUTWIN ;
  char *tmpName;
  for(icut=0; icut < NCUTWIN; icut++ ) {
    tmpName = INPUTS.CUTWIN_NAME[icut] ;
    if ( strcmp(tmpName,varname_gamma) == 0 ) { FOUND_GAMMA=1; }
  }
  if ( LDMP ) {
    if ( FOUND_GAMMA ) 
      { printf("\t %s already on CUTWIN list. \n", varname_gamma ); }
    else 
      { printf("\t Append %s to CUTWIN list. \n", varname_gamma ); }
  }

  if ( FOUND_GAMMA == 0 ) {
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
  char tmpName[100], *ptrTmp ;
  int  MXVAR = MXVARNAME_MISSING ;
  int  MEMC  = 60*sizeof(char);
  int  ndef, i, LEN ; 
  bool wildcard;
  char fnam[] = "prep_input_varname_missing" ;

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
  print_debug_malloc(+1,fnam);
  for(i=0; i < MXVAR; i++ ) 
    { INPUTS_VARNAME_MISSING.varname_list[i] = (char*)malloc(MEMC); }

  splitString(varname_missing, COMMA, MXVAR,         // inputs
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
      { fprintf(FP_STDOUT,
	       "\t Append missing column with varname = '%s' \n", ptrTmp ); 
      }
    else
      { fprintf(FP_STDOUT,
		"\t Append missing column(s) with substring = '%s'\n",ptrTmp);
      }

  }

  return ;
} // end   prep_input_varname_missing


// **********************************************
void  prep_cosmodl_lookup(void) {

  // Nov 22 2017
  // If all cosmo params are fixed, then make binned
  // dl vs, z lookup for faster computation.

  int debug_malloc = INPUTS.debug_malloc ;
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

  fprintf(FP_STDOUT, 
	  " Create cosmdl-vs-z lookup: %d z-bins from %.4f to %.4f \n",
	 NBZ, ZMIN, ZMAX);

  int MEMD = NBZ * sizeof(double) ;
  COSMODL_LOOKUP.NBZ  = NBZ;  
  COSMODL_LOOKUP.ZMIN = ZMIN ;
  COSMODL_LOOKUP.ZMAX = ZMAX ;
  COSMODL_LOOKUP.ZBIN = ZBIN ;

  print_debug_malloc(+1*debug_malloc,fnam);
  COSMODL_LOOKUP.z   = (double*) malloc(MEMD);
  COSMODL_LOOKUP.dl  = (double*) malloc(MEMD);

  for(iz=0; iz < NBZ; iz++ ) {
    di = (double)iz + 0.5 ;
    z  = ZMIN + (ZBIN * di) ;
    dl = cosmodl_forFit(z,z,INPUTS.COSPAR);
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }


  fprintf(fp,"\n %s\n\n", dotDashLine );

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) {
    fprintf(fp,"\t BEGIN FIT for SUBPROCESS ITER=%d\n", SUBPROCESS.ITER );
    fprintf(fp,"\n %s\n\n", dotDashLine );
    fflush(fp); return;
  }
#endif

  if ( INPUTS.NSPLITRAN == 1 ) {
    fprintf(fp,"\t BEGIN FIT \n" );
  }
  else {
    fprintf(fp,"\t BEGIN FIT with Random subset %d of %d\n",
	    NJOB_SPLITRAN, INPUTS.NSPLITRAN);
  }

  fprintf(fp,"\n %s\n\n", dotDashLine );
  fflush(fp);
  return;

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
void recalc_dataCov(void) {

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
  double a   = FITRESULT.ALPHA ; // ??? replace
  double b   = FITRESULT.BETA ;
  double g   = 0.0 ; 
  char fnam[] = "recalc_dataCov";

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
  // Mar 25 2022: check user input INPUTS.dchi2red_dsigint

  double parval_next;
  double step, slope=0.0, numer=0.0, denom=0.0 ;
  double slope_user = INPUTS.dchi2red_dsigint;
  int NFIT_ITER = FITRESULT.NFIT_ITER ;
  char fnam[] = "next_covFitPar" ;

  // ------------- BEGIN ------------

  parval_next = parval_orig; // init

  FITRESULT.CHI2RED_LIST[NFIT_ITER] = redchi2 ;
  FITRESULT.SIGINT_LIST[NFIT_ITER]  = parval_orig ;

  if ( slope_user != 0.0  ) {
    parval_next = parval_orig - (redchi2-1.0)/slope_user ;
    slope = slope_user ;
  }
  else if ( NFIT_ITER == 0 ){
    // decide if sigint needs to be larger or smaller
    if (redchi2 > 1.0) 
      { step = +parval_step ; } 
    else 
      { step = -parval_step ; }
    // calculate new sigint
    parval_next = parval_orig + step;

  } 
  else { 
    numer = 
      FITRESULT.CHI2RED_LIST[NFIT_ITER] - 
      FITRESULT.CHI2RED_LIST[NFIT_ITER-1] ;
    denom = 
      FITRESULT.SIGINT_LIST[NFIT_ITER] - 
      FITRESULT.SIGINT_LIST[NFIT_ITER-1] ;
    slope = numer/denom ;
    parval_next = parval_orig - (redchi2-1.0)/slope ;
  }

  // - - - - - - - - - - - 
  if ( parval_next > 100. ) {
    print_preAbort_banner(fnam);
    if ( NFIT_ITER > 0 ) {
      printf("\t slope = %f/%f = %f \n", numer, denom, slope);
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
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  FITINP.DCHI2RED_DSIGINT = slope ;
  return(parval_next);

} // end next_covFitPar



// ******************************************
void conflict_check() {

  // abort on conflict between incompatible inputs.

  int  NERR = 0, i ;
  char varName[10][2][20]; // [NERR][IVAR][varName]
  char fnam[] = "conflict_check" ;

  // ---------- BEGIN -----------

  if ( INPUTS.zpolyflag == 1 && INPUTS.fitflag_sigmb > 0 ) {
    sprintf(varName[NERR][0], "zpolyflag");
    sprintf(varName[NERR][1], "fitflag_sigmb");
    NERR++ ;
  }

#ifdef USE_SUBPROCESS   
  if (INPUTS.NSPLITRAN > 1 && SUBPROCESS.USE  ) {
    sprintf(varName[NERR][0], "NSPLITRAN");
    sprintf(varName[NERR][1], "SUBPROCESS");
    NERR++ ;
  }
#endif  
  
  if ( NERR > 0 ) {      
    print_preAbort_banner(fnam);
    for(i=0; i < NERR; i++ ) {
      printf("\t ERROR: conflict between inputs %s and %s\n", 
	     varName[i][0], varName[i][1] );
      fflush(stdout);
    }
    sprintf(c1err,"Found %d input conflicts (see above).", NERR );
    sprintf(c2err,"Check inputs.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( !INPUTS.cat_only ) {
    double zmin = INPUTS.zmin;
    double zmax = INPUTS.zmax;
    bool   define_zbinuser  = strlen(INPUTS.zbinuser) > 0 ;
    bool   define_zminzmax  = ( zmin > 0.0 && zmax < 9.0 );
    if ( !define_zminzmax && !define_zbinuser ) {
      sprintf(c1err,"Input file must specify zmin and zmax,");
      sprintf(c2err,"or specify zbinuser.");
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
  }

  return ;
} // end conflict_check


// ******************************************
void outFile_driver(void) {

  // Created Nov 30 2017
  // [move code out of main]
  // May 29 2019: call SPLITRAN_write_fitpar_legacy
  // Jun 24 2020: remove SPLIT[nnn] from NSPLITRAN file names.
  // Dec 04 2020: check option to write_M0_csv().

  int  JOBID     = INPUTS.JOBID_SPLITRAN ;
  int  NSPLITRAN = INPUTS.NSPLITRAN ;
  char *prefix   = INPUTS.PREFIX ;

  char tmpFile1[200], tmpFile2[200], tmpFile3[200], tmpFile[200];
  char fnam[] = "outFile_driver" ; 

  // --------------- BEGIN -------------

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) {
    SUBPROCESS_OUTPUT_LOAD();
    SUBPROCESS_OUTPUT_WRITE();

    if ( INPUTS.write_yaml ) {
      sprintf(tmpFile,"%s.YAML", prefix );
      write_yaml_info(tmpFile);
    }

    int OPTMASK   = SUBPROCESS.INPUT_OPTMASK ;
    bool WRFITRES = ( (OPTMASK & SUBPROCESS_OPTMASK_WRFITRES) > 0 );
    bool WRM0DIF  = ( (OPTMASK & SUBPROCESS_OPTMASK_WRM0DIF) > 0 );
    if ( WRFITRES ) {
      sprintf(tmpFile,"%s.FITRES", prefix );
      write_fitres_driver(tmpFile);  // write result for each SN
    }
    if ( WRM0DIF ) {
      sprintf(tmpFile,"%s.M0DIF", prefix );
      write_M0_fitres(tmpFile);  // write result for each SN 
    }

    return ;
  }
#endif
 

  if ( strlen(prefix) > 0 && !IGNOREFILE(prefix)  ) {

    sprintf(tmpFile1,"%s.fitres", prefix ); 
    sprintf(tmpFile2,"%s.M0DIF",  prefix ); 
    sprintf(tmpFile3,"%s.COV",    prefix );  // Dec 2 2020

    // Aug 12 2020 temp HACK: 
    // if YAML output is specified, it's from the new
    // submit_batch_jobs script. Write output fitres file with
    // .FITRES instead of .fitres to avoid the silly file-move
    // in batch script. Should fix this permanenetly, but need to 
    // check/fix SALT2mu_fit.pl.
    if ( INPUTS.write_yaml ) { sprintf(tmpFile1,"%s.FITRES", prefix ); }
    
    prep_blindVal_strings();
    write_fitres_driver(tmpFile1);  // write result for each SN
    write_M0_fitres(tmpFile2);      // write M0 vs. redshift
    write_M0_cov(tmpFile3);         // write cov_stat matrix for CosmoMC

    if ( INPUTS.write_yaml ) {
      sprintf(tmpFile,"%s.YAML", prefix );
      write_yaml_info(tmpFile);
    }
    if ( INPUTS.write_csv ) {
      sprintf(tmpFile,"%s.CSV", prefix );
      write_M0_csv(tmpFile);
    }

  } 
  else {  
    fprintf(FP_STDOUT, "\n PREFIX not specified --> no fitres output.\n");
    fflush(FP_STDOUT);
  }

  return ;

} // end outFile_driver

// ******************************************
void write_version_info(FILE *fp) {

  fprintf(fp,"# ISDATA_REAL:   %d \n", ISDATA_REAL);
  fprintf(fp,"# SNANA_VERSION: %s \n", SNANA_VERSION_CURRENT);
  fprintf(fp,"# BBC_VERSION:   %d \n", BBC_VERSION);
  fprintf(fp,"\n");
  fflush(fp);

} // end write_version_info

// ******************************************
void write_yaml_info(char *fileName) {

  // Aug 12, 2020
  // Write summary info to YAML file; to be used by batch-sumit script.
  // Write value and error for every entry; if no error, write error=0.
  // 
  // Jun 7 2021: write subprocess iteration
  // Sep 18 2021: write stats_bySAMPLE
  // Oct 06 2021: write ISDATA_REAL

  int  NDATA_REJECT_BIASCOR = NSTORE_CUTBIT[EVENT_TYPE_DATA][CUTBIT_BIASCOR] ;
  int  NDATA_PASS  = *NPASS_CUTMASK_POINTER[EVENT_TYPE_DATA]; 
  bool USE_DATA    = true ;
  bool USE_BIASCOR = (INPUTS.nfile_biasCor > 0) ;
  bool USE_CCPRIOR = (INPUTS.nfile_CCprior > 0) ;
  bool USE_EVENT_TYPE[MXEVENT_TYPE] = 
    { false, USE_DATA, USE_BIASCOR, USE_CCPRIOR };

  double t_cpu = (t_end_fit-t_start)/60.0 ;

  FILE *fp;
  int  NSN_PASS, NSN_PASS_LIST[10], evtype, idsample ;
  bool USE; 
  char KEY[60], ctmp[100], *str ;
  char fnam[] = "write_yaml_info" ;

  // -------------- BEGIN ----------------

  print_banner(fnam);

  fp = fopen(fileName,"wt");
  if ( !fp )  {
    sprintf(c1err,"Could not open YAML summary file:");
    sprintf(c2err,"%s", fileName);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  for(evtype=1; evtype < MXEVENT_TYPE; evtype++ ) {
    if ( USE_EVENT_TYPE[evtype] )
      { NSN_PASS  =  *NPASS_CUTMASK_POINTER[evtype]; }  // pass cuts
    else
      { NSN_PASS = 0 ; }
    str       =   STRING_EVENT_TYPE[evtype] ;
    NSN_PASS_LIST[evtype] = NSN_PASS; // used below for 'ALL' bySAMPLE
    sprintf(KEY,"NEVT_%s:", str);
    fprintf(fp,"%-22.22s %d \n", KEY, NSN_PASS );
  }

  sprintf(KEY,"NEVT_REJECT_BIASCOR:");
  fprintf(fp,"%-22.22s %d\n", KEY, NDATA_REJECT_BIASCOR );

  sprintf(KEY,"ABORT_IF_ZERO:");
  fprintf(fp,"%-22.22s %d\n", KEY, NDATA_PASS );

  sprintf(KEY,"CPU_MINUTES:");
  fprintf(fp,"%-22.22s %.2f\n", KEY, t_cpu);

  sprintf(KEY,"ISDATA_REAL:");
  fprintf(fp,"%-22.22s %d\n", KEY, ISDATA_REAL );

  // write NEVT_[WHAT]_bySAMPLE (e.g., LOWZ, SDSS, PS1, DES)
  if ( NSAMPLE_BIASCOR > 0 ) {
    fprintf(fp,"\n") ;  

    char STR_SAMPLE_LIST[200], STR_NEVT_LIST[MXEVENT_TYPE][100];
    STR_SAMPLE_LIST[0] = 0;
    for(evtype=1; evtype < MXEVENT_TYPE; evtype++ )
      { STR_NEVT_LIST[evtype][0] = 0 ; }

    // idsample = -1 is for 'ALL'
    for(idsample=-1; idsample < NSAMPLE_BIASCOR; idsample++ ) {

      if ( idsample < 0 ) 
	{ sprintf(ctmp,"ALL"); }
      else
	{ sprintf(ctmp," %s", SAMPLE_BIASCOR[idsample].NAME); }
      catVarList_with_comma(STR_SAMPLE_LIST, ctmp);

      for(evtype=1; evtype < MXEVENT_TYPE; evtype++ ) {
	if(idsample < 0 ) 
	  { NSN_PASS = NSN_PASS_LIST[evtype]; } // ALL samples combined
	else
	  { NSN_PASS = NPASS_CUTMASK_BYSAMPLE[evtype][idsample] ; }
	sprintf(ctmp," %d", NSN_PASS);
	catVarList_with_comma(STR_NEVT_LIST[evtype], ctmp);	
      }      
    }

    sprintf(KEY,"SAMPLE_LIST:");
    fprintf(fp,"%-22.22s %s\n", KEY, STR_SAMPLE_LIST);
    for(evtype=1; evtype < MXEVENT_TYPE; evtype++ ) {
      str       =   STRING_EVENT_TYPE[evtype] ;
      sprintf(KEY,"NEVT_%s_bySAMPLE:", str);
      fprintf(fp,"%-22.22s %s\n", KEY, STR_NEVT_LIST[evtype] );
    }
    fflush(fp);
  } // end NSAMPLE_BIASCOR>0

  

  if ( SUBPROCESS.ITER >= 0 ) {
    sprintf(KEY,"SUBPROCESS_ITERATION:");
    fprintf(fp,"%-22.22s %d\n", KEY, SUBPROCESS.ITER);
  }

  fprintf(fp,"\n") ;

  fprintf(fp,"BBCFIT_RESULTS:\n") ;

  int SIG_NSNFIT = (int)(sqrt(FITRESULT.NSNFIT)+0.5);
  fprintf(fp,"  - NSNFIT:       %5d     %d\n",   
	  FITRESULT.NSNFIT, SIG_NSNFIT ) ;

  fprintf(fp,"  - SIGINT:       %.5f    0.0\n", 
	  FITINP.COVINT_PARAM_FIX ) ;

  //fprintf(fp,"\n") ;
  
  bool ISFLOAT, ISM0;
  int n;   double VAL,ERR;  char tmpName[40];
  for ( n=0; n < FITINP.NFITPAR_ALL ; n++ ) {

    ISFLOAT = FITINP.ISFLOAT[n] ;
    ISM0    = n >= MXCOSPAR ; // it's z-binned M0
    if ( !ISFLOAT ) { continue ; }
    if ( ISM0     ) { continue ; }

    VAL = FITRESULT.PARVAL[NJOB_SPLITRAN][n] ;
    ERR = FITRESULT.PARERR[NJOB_SPLITRAN][n] ;
    sprintf(tmpName, "%s", FITRESULT.PARNAME[n]);
    trim_blank_spaces(tmpName);       strcat(tmpName,":") ;
    // if ( ERR < 0.0 ) { continue ; }
    
    fprintf(fp,"  - %-12.12s  %.5f  %.5f \n", tmpName, VAL, ERR ) ;
  }


  fclose(fp);

  return;

} // end write_yaml_info

// ******************************************
void  write_M0_fitres(char *fileName) {

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
  // Dec 02 2020: write redshift with %.5f instad of %.4f
  // Apr 22 2022: write %.5f format for distances/error instead of %.4f
  //               (for high-precision consistency tests)
  //
  int iz, irow, NFIT ;
  double z, zMIN, zMAX, VAL, ERR, dl, MUREF;
  char *tmpName, strval_OL[80], strval_w0[80];
  FILE *fp;
  char fnam[] = "write_M0_fitres" ;

  // ---------- BEGIN -----------

  if ( INPUTS.cutmask_write == -9 ) { return ; } // July 2016

  calc_zM0_data(); // fill FITRESULT.zM0[iz]

  fp = fopen(fileName,"wt");

  if ( !fp )  {
    sprintf(c1err,"Could not open M0-outFile");
    sprintf(c2err,"%s", fileName);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  fprintf(FP_STDOUT, "\n Write MUDIF vs. redshift to %s\n" , fileName); 

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

#ifdef TEXTFILE_NVAR
  int NVAR=8 ;
  fprintf(fp,"NVAR:  %d \n", NVAR);
#endif
  fprintf(fp,"VARNAMES: ROW  zMIN     zMAX     z        "
	  "MUDIF  MUDIFERR   MUREF  NFIT \n");
  irow = 0 ;

  for(iz=0; iz < INPUTS.nzbin; iz++ ) {

    irow++ ;    

    if ( INPUTS.opt_biasCor > 0 ) { 
      z     = FITRESULT.zM0[iz]; 
      MUREF = FITRESULT.MUREF_M0[iz]; 
    }  
    else {
      z = INPUTS.BININFO_z.avg[iz] ;
      dl    = cosmodl_forFit(z, z, INPUTS.COSPAR) ;
      MUREF = 5.0*log10(dl) + 25.0 ;
    }

    zMIN   = INPUTS.BININFO_z.lo[iz] ;
    zMAX   = INPUTS.BININFO_z.hi[iz] ;

    VAL   = FITRESULT.M0DIF[iz];
    ERR   = FITRESULT.M0ERR[iz];	

    NFIT = FITINP.NZBIN_FIT[iz] ;

    fprintf(fp, "ROW:     "
	    "%2d  %7.5f %7.5f %7.5f  "
	    "%9.5f %9.5f  %.5f %4d\n",
	    irow, zMIN, zMAX, z, 
	    VAL, ERR, MUREF, NFIT );
    fflush(fp);
  }

  fclose(fp);

  return ;

} // end write_M0_fitres


// ******************************************
void  write_M0_csv(char *fileName) {

  // Created Dec 4 2020
  // write M0 vs. z to csv file ; formatted for input to CosmoMC
  //

  int iz ;
  double z, M0DIF, M0ERR, dl, MUREF, MU, MUERR, zerr=0.0 ;
  FILE *fp ;
  char NAME[40];
  char fnam[] = "write_M0_fitres" ;

  // ---------- BEGIN -----------

  if ( INPUTS.cutmask_write == -9 ) { return ; } // July 2016

  fp = fopen(fileName,"wt");

  if ( !fp )  {
    sprintf(c1err,"Could not open M0-outFile");
    sprintf(c2err,"%s", fileName);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  fprintf(FP_STDOUT, "\n Write input HD for CosmoMC to %s\n" , fileName); 

  // - - - - 
  // write_version_info(fp);

  
  fprintf(fp,"# name   zcmb zhel zerr   mu muerr \n");

  for(iz=0; iz < INPUTS.nzbin; iz++ ) {

    sprintf(NAME,"BIN%2.2d", iz);

    z   = INPUTS.BININFO_z.avg[iz] ;
    if ( INPUTS.opt_biasCor > 0 ) { z = FITRESULT.zM0[iz]; }

    M0DIF   = FITRESULT.M0DIF[iz];
    M0ERR   = FITRESULT.M0ERR[iz];
    dl    = cosmodl_forFit(z, z, INPUTS.COSPAR) ;
    MUREF = 5.0*log10(dl) + 25.0 ;
    
    MU    = MUREF + M0DIF;
    MUERR = M0ERR ;

    fprintf(fp, "%s   %.5f %.5f %.1f   %.4f %.4f \n",
	    NAME, z, z, zerr, MU, MUERR);
    fflush(fp);
  }

  fclose(fp);

  return ;

} // end write_M0_csv


// ******************************************
void write_M0_cov(char *fileName) {

  // Created Dec 2, 2020
  // write M0 cov_stat (Nz x Nz) from fit to text file.
  // Write format for CosmoMC to read:
  //

  int NFITPAR_ALL =  FITINP.NFITPAR_ALL ;
  int iz0, iz1, ipar0, ipar1, iMN0, iMN1 ;
  double COV ;
  FILE *fp;
  char fnam[] = "write_M0_cov" ;

  // ---------- BEGIN -----------

  if ( INPUTS.cutmask_write == -9 ) { return ; } 

  fp = fopen(fileName,"wt");

  if ( !fp )  {
    sprintf(c1err,"Could not open cov-outFile");
    sprintf(c2err,"%s", fileName);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  fprintf(FP_STDOUT, "\n Write COV_stat matrix to %s\n" , fileName); 

  // - - - - 
  //  write_version_info(fp);

  fprintf(fp,"%d\n", INPUTS.nzbin );

  for ( ipar0=MXCOSPAR ; ipar0 < NFITPAR_ALL ; ipar0++ ) {

    iz0  = INPUTS.izpar[ipar0] ; 
    if ( iz0 >= INPUTS.nzbin ) { continue; }

    fprintf(fp, "# ---------- begin row %d ------------ \n", iz0 );
    for ( ipar1=MXCOSPAR ; ipar1 < NFITPAR_ALL ; ipar1++ ) {

      iz1  = INPUTS.izpar[ipar1] ; 
      if ( iz1 >= INPUTS.nzbin ) { continue; }

      iMN0 = FITINP.IPARMAPINV_MN[ipar0];
      iMN1 = FITINP.IPARMAPINV_MN[ipar1];
      if ( iMN0 >=0 && iMN1 >= 0 ) 
	{ COV = FITRESULT.COVMAT[iMN0][iMN1]; }
      else  { 
	COV = 0.0 ;
	if ( ipar0==ipar1 ) { COV = 1.0E5 ; } // huge for diagonal
      }

      fprintf(fp, "%le\n", COV);      
    } // end ipar1
    fflush(fp);
  }  // end ipar0


  fclose(fp);

  return ;

} // end write_M0_cov


// ******************************************
void write_fitres_driver(char* fileName) {

  // Write outout in fitres format.
  // Combine original fitres variables with those
  // calculated here in SALT2mu.
  //
  // Mar 01 2018: add M0 to output
  // Jun 10 2019: call printCOVINT 
  // Jan 09 2019: fix to loop over each datafile instead of only the first.
  // May 09 2020: 
  //   + refactor to write columns defined in OUTPUT_VARNAMES struct.
  //   + call define_varnames_append()
  //   + check cat_only option that skips fit and just catenates

  //  bool  DO_BIASCOR_MU     = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  //  bool  IS_SIM            = INFO_DATA.TABLEVAR.IS_SIM ;

  bool  cat_only         = INPUTS.cat_only ;
  int   IWD_KEY  = 0;
  int   IWD_CCID = 1;

  double VAL, ERR, PULL ;
  FILE  *fout, *finp;

  int n, ivar, indx, NCUT, icut, cutmask, NWR, NLINE, ISFLOAT, iz, GZIPFLAG ;
  int idsample, NSN_DATA, nfile, ifile ;
  bool VALID_ROWKEY;
  char  line[MXCHAR_LINE], tmpName[60], *ptrFile ;
  char  ztxt[60], KEY[MXCHAR_VARNAME], CCID[40];

  char fnam[] = "write_fitres_driver" ;

  // ------------------ BEGIN ----------------

  NSN_DATA    =  INFO_DATA.TABLEVAR.NSN_ALL;

  // define the new fitres variables to add to the original list  
  define_varnames_append();  // sets NVAR_APPEND and VARNAMES_APPEND

  // - - - - - - - - - -

  fprintf(FP_STDOUT, "\n Open output file with  cutmask_write=%d : \n", 
	 INPUTS.cutmask_write );
  fprintf(FP_STDOUT, "\t %s \n", fileName );  
  fflush(FP_STDOUT);

  
  // - - - - - - - -  -
  fout = fopen(fileName,"wt");
  if (!fout ) {
    if ( cat_only ) {
      sprintf(c1err,"Could not open catfile_out='%s'", INPUTS.cat_file_out);
      sprintf(c2err,"Check cat_file_out key." );
    }
    else {
      sprintf(c1err,"Could not open output fitres file");
      sprintf(c2err,"'%s' ", fileName );
    }
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
  }

  if ( cat_only ) {
    write_cat_info(fout); // write cat info to output file header
  }
  else {
    write_version_info(fout);
    fprintf(fout,"# %s\n", STRING_MINUIT_ERROR[INPUTS.minos]);
    fprintf(fout,"# NCALL_FCN: %d \n", FITRESULT.NCALL_FCN );
    fprintf(fout,"# CPU: %.2f minutes\n",
	    (t_end_fit-t_start_fit)/60.0  );
    if ( INPUTS.blindFlag > 0 && ISDATA_REAL ) 
      { write_blindFlag_message(fout); }      
    fprintf(fout,"# MU-RESIDUAL NOTE: MURES = MU-(MUMODEL+M0DIF) \n");
    write_NWARN(fout,1);
    write_MUERR_INCLUDE(fout);
  }

  // - - - - -

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

  if ( cat_only ) { goto WRITE_TABLE_ROWS; }

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

  indx = NWR  = NLINE = 0;
  fflush(fout);

  // print contamination tables if CC prior is used
  if ( INFO_CCPRIOR.USE ) { print_contam_CCprior(fout);  }

  // check option to NOT write each SN to have smaller file
  // with only the fit results
  if ( INPUTS.cutmask_write == -9 ) { 
    fclose(fout);
    fprintf(FP_STDOUT, " Nothing written out  (%d/%d used in fit) \n", 
	   FITRESULT.NSNFIT , NSN_DATA );  fflush(FP_STDOUT);    
    return ;
  }

  // - - - - - - - -
  // re-read each data file

 WRITE_TABLE_ROWS:

  nfile = INPUTS.nfile_data;
  for(ifile=0; ifile < nfile; ifile++ ) {

    ptrFile = INPUTS.dataFile[ifile];
    finp  = open_TEXTgz(ptrFile, "rt", &GZIPFLAG); 

    while ( fgets (line, MXCHAR_LINE, finp) !=NULL  ) {

      if ( line[0] == ' '    ) { continue ; }
      if ( strlen(line) < 3  ) { continue ; }
      if ( commentchar(line) ) { continue; }

      store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,line);

      // skip if first wd of line is not a valid row key
      get_PARSE_WORD(0, IWD_KEY, KEY);
      VALID_ROWKEY = 
	( strcmp(KEY,"SN:" )    == 0 ) ||
	( strcmp(KEY,"ROW:")    == 0 ) ||
	( strcmp(KEY,"GAL:")    == 0 ) ;
      if ( !VALID_ROWKEY ) { continue ; }
      // xxx mark      if ( strcmp(KEY,"SN:") != 0 ) { continue ; }

      if ( cat_only ) {
	// check prescale
	NLINE++ ;
	if ( NLINE % INPUTS.cat_prescale != 0 ) { continue ; }
      }
      else {
	// check cutmask for writing events
	get_PARSE_WORD(0, IWD_CCID, CCID);
	get_CCIDindx(CCID, &indx) ;      
	cutmask = INFO_DATA.TABLEVAR.CUTMASK[indx]; 
	if ( !keep_cutmask(cutmask)  ) { continue; }
      }

      NWR += write_fitres_line(indx,ifile, KEY,line,fout);

    }  // end reading line with fgets

    fclose(finp);   
  } // end ifile

  fclose(fout);

  if ( cat_only ) {
    fprintf(FP_STDOUT, " Wrote %d SN to cat table.\n", NWR); 
  }
  else {
    fprintf(FP_STDOUT, " Wrote %d SN  (%d/%d used in fit) \n", 
	   NWR, FITRESULT.NSNFIT , NSN_DATA );
  }
  fflush(FP_STDOUT);

  return ;

} // end of write_fitres_driver

// ===============================================
int write_fitres_line(int indx, int ifile, char *rowkey, 
		      char *line, FILE *fout) {

  // for input 'ifile' and 'line', write to fout.
  // Note that store_PARSE_WORDS has already been called,
  // so here use get_PARSE_WORDS to retrieve.
  //
  // Each *line is split into words in case different files have
  // different columns.
  //  
  // Nov 14 2020: check for datafile_overrides.
  // Dec 28 2021: pass rowkey to allow for prescaling HOSTLIB

  int NVAR_TOT = OUTPUT_VARNAMES.NVAR_TOT ;  
  int ISTAT = 0 ;
  int  ivar_tot, ivar_file, ivar_word ;
  char word[MXCHAR_VARNAME], line_out[MXCHAR_LINE], *varName ;
  char blank[] = " ";
  char fnam[] = "write_fitres_line" ;
  int  LDMP = 0 ;
  // ----------- BEGIN -----------


  // xxx mark delete  sprintf(line_out,"SN: ");
  sprintf(line_out,"%s ", rowkey);

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

    // Nov 15 2020: check for data override 
    write_word_override(ivar_tot,indx,word); 

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

// ===================
void  write_word_override(int ivar_tot, int indx, char *word) {

  // Nov 15 2020
  // Check to overwrite word for ivar_tot and SN indx.
  //
  // Jan 27 2021: write any zXXX variable with %.5f format.

  int NVAR_OVERRIDE  = INFO_DATA.NVAR_OVERRIDE ;
  int ivar_over ;
  bool IS_z;
  char *varName ;
  char fnam[] = "write_word_override";

  // ---------- BEGIN -------------

  if ( NVAR_OVERRIDE == 0 ) { return; }

  ivar_over = INFO_DATA.IVAR_OUTPUT_INVMAP[ivar_tot];
  if ( ivar_over < 0 ) { return; }

  // this variable gets an override value fval:
  float fval = INFO_DATA.PTRVAL_OVERRIDE[ivar_over][indx];

  // for formatting, check of variable contains zHD
  varName   = OUTPUT_VARNAMES.LIST[ivar_tot]; 
  IS_z      = ( varName[0] == 'z' );

  // replace word string
  if ( IS_z ) 
    { sprintf(word,"%.5f", fval); }
  else
    { sprintf(word,"%.3f", fval); }

  return;

} // end write_word_override

// ===============================================
void define_varnames_append(void) {

  // Nov 12 2020: add MUERR_VPEC
  // Dec 02 2020: add IZBIN & M0DIFERR

  bool DO_BIASCOR_MU  = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  bool DO_COVSCALE    = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD      = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  int   NSN_BIASCOR       =  INFO_BIASCOR.TABLEVAR.NSN_ALL;
  char  tmpName[MXCHAR_VARNAME];
  //  char fnam[] = "define_varnames_append";

  // ----------- BEGIN -----------

  NVAR_APPEND = 0 ;
  if ( INPUTS.cat_only) { return; }

  sprintf(VARNAMES_APPEND[NVAR_APPEND],"CUTMASK");      NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MU");           NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUMODEL");      NVAR_APPEND++ ;  

  // distance error from BBC
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR");        NVAR_APPEND++ ;  

  // see ??
  if ( !SUBPROCESS.USE )
    { sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR_RENORM"); NVAR_APPEND++ ; }

  // contribution from LC fit only (no sigInt, no VPEC, no scale)
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR_RAW");    NVAR_APPEND++ ;  

  // contribution from VPEC only, no scale
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUERR_VPEC");   NVAR_APPEND++ ;  

  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MURES");        NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"MUPULL");       NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"M0DIF");        NVAR_APPEND++ ;  
  sprintf(VARNAMES_APPEND[NVAR_APPEND],"M0DIFERR");     NVAR_APPEND++ ;  

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
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_muCOVSCALE");   NVAR_APPEND++ ;  
    if ( DO_COVADD ) {
      sprintf(VARNAMES_APPEND[NVAR_APPEND],"biasCor_muCOVADD");   NVAR_APPEND++ ;
    }
    sprintf(VARNAMES_APPEND[NVAR_APPEND],"IDSAMPLE");          NVAR_APPEND++ ;  
    sprintf(VARNAMES_APPEND[NVAR_APPEND],VARNAME_IZBIN);    NVAR_APPEND++ ;  

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
  // Oct 06 2021: remove obsolete line with zERR

  int USE=0;
  double tmpErr;
  //  char fnam[] = "write_MUERR_INCLUDE";

  // --------------- BEGIN ----------------

  // xxx mark delete   fprintf(fp, "# MUERR_INCLUDE: zERR \n" );  USE=1;

  tmpErr = INPUTS.zpecerr ;
  if ( tmpErr > 0.0 ) {
    fprintf(fp, "# MUERR_INCLUDE: zPECERR=%.5f \n", tmpErr );
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
    
    fflush(FP_STDOUT);
  }

  fprintf(fout,"#  -2log(L)     = %.2f \n", FITRESULT.CHI2SUM_MIN );

  fprintf(fout,"#  dchi2red/dsigint = %.3f\n",
	  FITINP.DCHI2RED_DSIGINT);

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

  if (INPUTS.cat_prescale > 1) {
    fprintf(fout,"# PRESCALE = %d\n",
	    INPUTS.cat_prescale);
  }

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
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
  // Nov 12 2020: write muerr_vpec for Dan.
  // Dec 02 2020: write izbin
  // Apr 22 2022: mu/mumodel/muerr -> %8.5f format instead of %7.4f
  //              Might be needed later for high-precision tests.
  //
  bool DO_BIASCOR_MU = (INPUTS.opt_biasCor & MASK_BIASCOR_MU );
  bool DO_COVSCALE   = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVSCALE) > 0;
  bool DO_COVADD     = (INPUTS.opt_biasCor & MASK_BIASCOR_MUCOVADD) > 0;

  double mu, muerr, muerr_renorm, muerr_raw, muerr_vpec, mumodel, mures, pull;
  double M0DIF, M0ERR ;
  double muBias=0.0, muBiasErr=0.0,  muCOVscale=0.0, chi2=0.0, muCOVadd=0.0 ;
  double fitParBias[NLCPAR] = { 0.0, 0.0, 0.0 } ;
  int    n, cutmask, NWR, NSN_BIASCOR, idsample, izbin;
  char line[400], word[40] ;	 
  char fnam[] = "write_fitres_line_append" ;

  // ------------ BEGIN --------------

  n =  indx;

  NSN_BIASCOR =  INFO_BIASCOR.TABLEVAR.NSN_ALL;

  cutmask    = INFO_DATA.TABLEVAR.CUTMASK[n]  ;
  idsample   = INFO_DATA.TABLEVAR.IDSAMPLE[n]  ;
  izbin      = INFO_DATA.TABLEVAR.IZBIN[n]  ;

  //  z        = INFO_DATA.TABLEVAR.zhd[n] ;  
  //  zerr     = INFO_DATA.TABLEVAR.zhderr[n] ;
  mumodel    = INFO_DATA.mumodel[n];
  mu            = INFO_DATA.mu[n] - FITRESULT.SNMAG0; 
  muerr         = INFO_DATA.muerr[n];
  muerr_raw     = INFO_DATA.muerr_raw[n] ; // from LC fit only
  muerr_vpec    = INFO_DATA.muerr_vpec[n] ;
  mures         = INFO_DATA.mures[n] ;
  pull          = INFO_DATA.mupull[n] ;
  M0DIF         = INFO_DATA.M0[n] - FITRESULT.AVEMAG0 ;
  M0ERR         = FITRESULT.M0ERR[izbin];
  chi2          = INFO_DATA.chi2[n] ;

  //  sim_mb   = INFO_DATA.TABLEVAR.SIM_FITPAR[INDEX_mB][n] ;
  //  sim_mu   = INFO_DATA.TABLEVAR.SIM_MU[n] ;

  if ( NSN_BIASCOR > 0 ) { 
    muBias               = INFO_DATA.muBias[n] ;
    muBiasErr            = INFO_DATA.muBiasErr[n] ;
    if ( DO_BIASCOR_MU == false ) {
      fitParBias[INDEX_mB] = INFO_DATA.fitParBias[n][INDEX_mB] ;
      fitParBias[INDEX_x1] = INFO_DATA.fitParBias[n][INDEX_x1] ;
      fitParBias[INDEX_c]  = INFO_DATA.fitParBias[n][INDEX_c] ;    
    }
    muCOVscale  = INFO_DATA.muCOVscale[n]  ;
    if ( DO_COVADD ) {
      muCOVadd    = INFO_DATA.muCOVadd[n] ;
      if ( muCOVscale > 1.0 && muCOVadd > 0.0 ) { 
	muCOVscale = 1.0; 
      } else {
	muCOVadd = 0.0;
      }
    }
  }
  
  if (pull > 99.999) { pull=99.999; }
  
  NWR=0;  line[0] = 0 ;
  sprintf(word, "%d ",    cutmask);       NWR++ ; strcat(line,word);
  sprintf(word, "%8.5f ", mu     );       NWR++ ; strcat(line,word);
  sprintf(word, "%8.5f ", mumodel);       NWR++ ; strcat(line,word);
  sprintf(word, "%8.5f ", muerr  );       NWR++ ; strcat(line,word);

  if ( !SUBPROCESS.USE ) {
    // Jun 21 2021: beware that muerr_renorm is not malloced or filled
    //     for SUBPROCESS   
    muerr_renorm  = INFO_DATA.muerr_renorm[n];
    sprintf(word, "%7.4f ", muerr_renorm ); NWR++ ; strcat(line,word);
  }
  sprintf(word, "%7.4f ", muerr_raw );    NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", muerr_vpec );   NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", mures  );       NWR++ ; strcat(line,word);
  sprintf(word, "%6.3f ", pull   );       NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", M0DIF  );       NWR++ ; strcat(line,word);
  sprintf(word, "%7.4f ", M0ERR  );       NWR++ ; strcat(line,word);
  sprintf(word, "%.2f ",  chi2   );       NWR++ ; strcat(line,word);

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
    if ( DO_COVADD ) {
      sprintf(word, "%6.4f ", muCOVadd ) ;    NWR++ ; strcat(line,word);
    }
    sprintf(word, "%d "   , idsample ) ;      NWR++ ; strcat(line,word);
    sprintf(word, "%2d "  , izbin ) ;         NWR++ ; strcat(line,word);

  }

  
  fprintf(fp,"%s", line);
  fflush(fp);

  if ( NWR != NVAR_APPEND ) {
    sprintf(c1err,"Expected to write NVAR_APPEND=%d SALT2mu variables",
	    NVAR_APPEND);
    sprintf(c2err,"but wrote only %d variables.", NWR);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
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
    fprintf(FP_STDOUT,
	   "Average mag0 offset = %f  (wgted by NSN per z bin)\n", ave);
    fflush(FP_STDOUT);
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
void muerr_renorm(void) {

  // Created Dec  2020
  //
  // Compute *muerr_renorm for each event such that
  // in each BBC redshift bin,
  //
  //    sum [1/muerr_renorm^2] = 1/M0DIFERR^2
  //
  // Use constant scale within each redshift bin until somebody
  // figures out a better recipe.
  // 
  // As a diagnostic crosscheck:
  // For each z-bin, use unbinned values to compute weighted 
  // M0DIF-mean and M0DIF-error; flag discrepancies > 0.001 mag.
  //
  // May 25 2021: write to FP_STDOUT
  // Nov 22 2021: include missing factor of 1/sqrt(pia)

  int NSN_DATA   = INFO_DATA.TABLEVAR.NSN_ALL ;  
  int MEMD       = NSN_DATA * sizeof(double);
  int debug_malloc = INPUTS.debug_malloc ;
  int iz, isn, cutmask ;
  double SUM_WGT[MXz], SUM_MURES[MXz], MURES_check[MXz], M0ERR_check[MXz];
  double RATIO_MUERR[MXz], DIF_MURES[MXz], DIF_VARIANCE[MXz];
  double mumodel, mu, muerr, mures, WGT, ratio, dif, pcc, pia ;
  double tol_warn = 0.001;
  double pia_min  = 1.0E-6;
  char star_mures[2] ;
  char fnam[] = "muerr_renorm" ;

  // --------- BEGIN -----------

#ifdef USE_SUBPROCESS
  if ( SUBPROCESS.USE ) { return; }
#endif

  fprintf(FP_STDOUT,
	  "\n  %s: compute MUERR_RENORM to preserve M0DIF wgt per z bin\n",
	 fnam );
  fflush(FP_STDOUT);

  print_debug_malloc(+1*debug_malloc,fnam);
  INFO_DATA.muerr_renorm = (double*) malloc(MEMD);

  for(iz=0; iz < MXz; iz++ )  { SUM_WGT[iz] = SUM_MURES[iz] = 0.0 ; }

  for(isn=0; isn < NSN_DATA; isn++ ) {

    cutmask    = INFO_DATA.TABLEVAR.CUTMASK[isn]  ;
    if ( cutmask ) { continue; }

    iz         = INFO_DATA.TABLEVAR.IZBIN[isn] ;
    mumodel    = INFO_DATA.mumodel[isn];
    mu         = INFO_DATA.mu[isn] - FITRESULT.SNMAG0; 
    muerr      = INFO_DATA.muerr[isn];
    mures      = INFO_DATA.mures[isn] ;
    pcc        = INFO_DATA.probcc_beams[isn];
    pia        = 1.0 - pcc;

    WGT             = pia / (muerr*muerr) ;
    SUM_WGT[iz]    += WGT;
    SUM_MURES[iz]  += (mures * WGT) ;
    
  } // end isn


  // - - - - -
  // compute diagnostic for each z bin
  for(iz=0; iz < INPUTS.nzbin ; iz++ ) {
    RATIO_MUERR[iz]  = 0.0;
    DIF_MURES[iz]    = 0.0;

    if ( SUM_WGT[iz] == 0.0 ) { continue ; }

    MURES_check[iz] = SUM_MURES[iz] / SUM_WGT[iz] ;
    M0ERR_check[iz] = 1.0 / sqrt(SUM_WGT[iz]) ;

    RATIO_MUERR[iz] = M0ERR_check[iz] / FITRESULT.M0ERR[iz] ;
    DIF_MURES[iz]   = MURES_check[iz] ; 
  }

  // - - - - - 
  // loop over data again and compute muerr_renorm
  for(isn=0; isn < NSN_DATA; isn++ ) {
    cutmask    = INFO_DATA.TABLEVAR.CUTMASK[isn]  ;
    if ( cutmask ) { continue; }
    muerr      = INFO_DATA.muerr[isn];
    iz         = INFO_DATA.TABLEVAR.IZBIN[isn] ;
    pcc        = INFO_DATA.probcc_beams[isn];
    pia        = 1.0 - pcc;
    if ( pia < pia_min ) { pia = pia_min; } // avoid divide-by-zero below

    ratio      = RATIO_MUERR[iz];
    INFO_DATA.muerr_renorm[isn] = (muerr/ratio)/sqrt(pia) ;
  }


  // print stuff
  int NERR = 0 ;
  for(iz=0; iz < INPUTS.nzbin ; iz++ ) {
    if ( SUM_WGT[iz] == 0.0 ) { continue ; }

    ratio = RATIO_MUERR[iz];
    dif   = DIF_MURES[iz];

    sprintf(star_mures," ");
    if ( fabs(dif) > tol_warn ) { NERR++; sprintf(star_mures,"*"); }

    fprintf(FP_STDOUT,
	     "    <z_%2.2d>=%.3f  muerr *= %.3f  [MURES check = %7.4f%s]\n", 
	   iz, INPUTS.BININFO_z.avg[iz], 1.0/ratio, dif, star_mures );
    fflush(FP_STDOUT);

    /* xxx
    printf("     %2d  %6.4f   %7.4f/%7.4f%s  \n",
	   iz, INPUTS.BININFO_z.avg[iz],
	   FITRESULT.M0DIF[iz], M0DIF_check[iz], star_avg );
    */
  }


  if ( NERR > 0 ) {
    fprintf(FP_STDOUT,
	    " WARNING: %d of %d z bins fail M0DIF check with tol=%f \n",
	   NERR, INPUTS.nzbin, tol_warn );
    fflush(FP_STDOUT);
  }



  return ;

} // end muerr_renorm



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
  // Dec 2 2020: move mnemat_ call to exec_mnpout_mnerrs.

  int num, iMN, jMN, i, j ;
  double cov, cov_diag, terr[MAXPAR], corr ;
  char tmpName[MXCHAR_VARNAME], msg[100];

  // --------- BEGIN ----------

  if ( NPARz_write > MXCOSPAR ) 
    { sprintf(msg, "Reduced COV matrix:"); }
  else
    { sprintf(msg, "Reduced COV matrix, includes %d m0_ bins:", 
	      NPARz_write); }

  fprintf(fp, "\n# %s\n", msg); fflush(fp);

  for ( iMN=0 ; iMN < NPAR_FLOAT ; iMN++ )
    {
      i  = FITINP.IPARMAP_MN[iMN]  ;
      if ( i - MXCOSPAR >= NPARz_write ) { continue; }
      sprintf(tmpName, "%s", FITRESULT.PARNAME[i] );
      tmpName[10] = 0;    // force truncation of name
      cov_diag    = FITRESULT.COVMAT[iMN][iMN];
      terr[iMN] = sqrt(cov_diag);

      //      fprintf(fp,"%2i ",iMN+1); 
      fprintf(fp,"# %-10.10s ", tmpName ); 

      for (jMN=0; jMN <= iMN; ++jMN)	{
	j  = FITINP.IPARMAP_MN[jMN]  ;
	if ( j - MXCOSPAR >= NPARz_write ) { continue; }
	cov  = FITRESULT.COVMAT[iMN][jMN];
	corr = cov/(terr[iMN]*terr[jMN]);
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
double cosmodl_forFit(double zhel, double zcmb, double *cosPar) {

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
      DL[i] = cosmodl(zhel,zcmb,cosPar_local);
    }

    slp = (DL[1] - DL[0]) / (OL_extrap[1] - OL_extrap[0]) ;
    dl  = DL[0] + ( OL - OL_extrap[0] )*slp ;
  }
  else {
    dl = cosmodl(zhel,zcmb,cosPar);
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

double cosmodl(double zhel, double zhd, double *cosPar)
{
  // Dec 11 2020: 
  // pass both zhel and zhd, where zhd has both cmb and vpec corrections.
 
  const double  cvel = LIGHT_km; // 2.99792458e5;
  const double  tol  = 1.e-6;
  double dflat, distance, H0inv ;
  double omega_k, OK, dl  ;
  //  char fnam[] = "cosmodl" ;

  // ------------- BEGIN --------------

  //  omega_l = cosPar[0];  // not used
  omega_k = cosPar[1];
  //  wde     = cosPar[2]; // not used
  //  wa      = cosPar[3]; // not used

  if(fabs(omega_k)<tol) { omega_k = 0.0; }
  OK      = fabs(omega_k);

  //comoving distance to redshift
  dflat = rombint(inc, 0.0, zhd, cosPar, tol);

  H0inv = 1.0/INPUTS.H0 ;

  if( omega_k == 0.0  ) 
    { distance = cvel*H0inv*dflat; }
  else if(omega_k<0.0) 
    { distance = cvel*H0inv*(1.0/sqrt(OK))*sin(sqrt(OK)*dflat); }
  else 
    { distance = cvel*H0inv*(1.0/sqrt(OK))*sinh(sqrt(OK)*dflat); }

  dl = (1.0+zhel)*distance ;
  return( dl );

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
    fprintf(FP_STDOUT, "Rombint failed to converge; integral = %e +- %e \n",
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
      fprintf(FP_STDOUT, 
	     "LU decomposition (ludcmp.c) : Singular matrix, but continue. "
	     "%d %d %d %d %d\n",i,j,i1,i2,ndim);
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


#ifdef USE_SUBPROCESS

// =======================================================
//  SUBPROCESS functions to use SALT2mu as 
//  function to return chi2-info for higher-level
//  MCMC fitter for SN population parameters.
// =======================================================


void SUBPROCESS_MALLOC_INPUTS(void) {
  // malloc SUBPROCESS.INPUT_xxx arrays; called just before reading
  // SUBPROCESS_FILES argument
  int i;
  int debug_malloc = INPUTS.debug_malloc ;
  char fnam[] = "SUBPROCESS_MALLOC_INPUTS" ;

  print_debug_malloc(+1*debug_malloc,fnam);

  SUBPROCESS.INPUT_FILES = 
    (char*) malloc( MXCHAR_FILENAME*3*sizeof(char) );

  SUBPROCESS.INPUT_CID_REWGT_DUMP =
    (char*) malloc( 2*MXCHAR_VARNAME*MXVAR_GENPDF*sizeof(char) );

  SUBPROCESS.INPUT_VARNAMES_GENPDF_STRING = 
    (char*) malloc( 2*MXCHAR_VARNAME*MXVAR_GENPDF*sizeof(char) );

  SUBPROCESS.INPUT_SIMREF_FILE =
    (char*) malloc( MXCHAR_FILENAME*sizeof(char) );

  SUBPROCESS.N_OUTPUT_TABLE = 0 ;
  SUBPROCESS.INPUT_OUTPUT_TABLE =  (char**) malloc( 10*sizeof(char*) );
  for(i=0; i < 10; i++ ) {
    SUBPROCESS.INPUT_OUTPUT_TABLE[i] =  (char*) malloc( 100*sizeof(char) );
    SUBPROCESS.INPUT_OUTPUT_TABLE[i][0] =  0;
  }

  SUBPROCESS.INPUT_FILES[0]          = 0;
  SUBPROCESS.INPUT_CID_REWGT_DUMP[0] = 0 ;
  SUBPROCESS.INPUT_VARNAMES_GENPDF_STRING[0] = 0;
  SUBPROCESS.INPUT_SIMREF_FILE[0] = 0;

    return ;
} // end SUBPROCESS_MALLOC_INPUTS


// ==================================
void  SUBPROCESS_HELP(void) {

  printf("\n");
  printf("\t   ********** SUBPROCESS HELP MENU ********** \n");
  printf("\n");
  printf("First argument is name of SALT2mu input file.\n"
	 "\n"
	 "SUBPROCESS_FILES=inpFile,outFile,stdoutFile \n"
	 "\t inpFile = name of file with input GENPDF maps (to read)\n"
	 "\t outFile = name of file with SALT2mu output (to write)\n"
	 "\t stdoutFile = name of log file containint stdout from SALT2mu\n"
	 "\n"
	 "SUBPROCESS_VARNAMES_GENPDF="
	 "<comma sep list of FITRES columns used for PDF maps>\n"
	 "\t e.g., SIM_x1,HOST_LOGMASS,SIM_c,SIM_RV,HOST_LOGMASS \n"
	 "\t Duplicate column names ok; they are internally handled.\n"
	 "\n"
	 "SUBPROCESS_OUTPUT_TABLE="
	 "<multiD output table: NBIN,min,max per dimension>\n"
	 "\t e.g., 'c(12,-0.3:0.3)*HOST_LOGMASS(2,0:20)'\n"
	 "\n"
	 "SUBPROCESS_CID_REWGT_DUMP=1,2,4,8654,9874    (optional) \n"
	 "\t list of CIDs to dump reweight info at each iteration.\n"
	 "\t CID < 10 -> isn index (e.g. CID=2 -> dump 2nd event)\n"
	 "\t CID > 10 -> dump this exact CID\n"
	 "\n" 
	 "Example of full SUBPROCESS command:\n"
	 "SALT2mu.exe SALT2mu_SIMDATA.input \\\n"
	 "   SUBPROCESS_FILES="
	 "SUBPROC_MAPS.DAT,SUBPROC_OUT.DAT,SUBPROC_LOG.STDOUT \\\n"
	 "   SUBPROCESS_VARNAMES_GENPDF="
	 "SIM_x1,HOST_LOGMASS,SIM_c,SIM_RV,HOST_LOGMASS \\\n"
	 "   SUBPROCESS_OUTPUT_TABLE="
	 "'c(12,-0.3:0.3)*RV(4,1:5)*HOST_LOGMASS(2,0:20)' \\\n"
	 "   SUBPROCESS_SNID_REWGT_DUMP=1,2,5177316\n"
	 );
	
  exit(0);
} // end SUBPROCESS_HELP

// ==================================
void  SUBPROCESS_INIT(void) {

  // Created July 2020
  // Parse SUBPROCESS.INPUT_FILES to get file names for:
  //   input PDF map file (written by python driver)
  //   output info        (passed from SALT2mu to python driver)
  //   log file           (std out created by SALT2mu)
  // These files communicate with python program.

  int NSN_DATA   = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  MEMC      =  MXCHAR_FILENAME * sizeof(char) ;
  int debug_malloc = INPUTS.debug_malloc ;
  int  NSPLIT, itable ;
  char *tmpFiles[3];
  char fnam[] = "SUBPROCESS_INIT" ;

  // ------------- BEGIN ----------

  if ( !SUBPROCESS.USE ) { return; }

  printf("\n\n%s  Begin %s\n", KEYNAME_SUBPROCESS_STDOUT, fnam );

  print_debug_malloc(+1*debug_malloc,fnam);

  // break comma-sep FILES into INPFILE and OUTFILE
  SUBPROCESS.INPFILE     = (char*) malloc(MEMC);
  SUBPROCESS.OUTFILE     = (char*) malloc(MEMC);
  SUBPROCESS.STDOUT_FILE = (char*) malloc(MEMC);
  tmpFiles[0]   = SUBPROCESS.INPFILE ;
  tmpFiles[1]   = SUBPROCESS.OUTFILE ;
  tmpFiles[2]   = SUBPROCESS.STDOUT_FILE ;
  splitString(SUBPROCESS.INPUT_FILES, ",", 3, &NSPLIT, tmpFiles);
  
  // open INPFILE in read mode, but only for sim data.
  // skip for real data since there is nothing to rewgt.
  if ( !ISDATA_REAL ) {
    SUBPROCESS.FP_INP = fopen(SUBPROCESS.INPFILE, "rt");
    if ( !SUBPROCESS.FP_INP ) {
      sprintf(c1err,"Could not open input GENPDF file to read:" );
      sprintf(c2err," '%s' ", SUBPROCESS.INPFILE);
      SUBPROCESS_REMIND_STDOUT();
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err); 
    }
    else {
      printf("%s  Opened input  file (GENPDF map): %s\n", 
	     KEYNAME_SUBPROCESS_STDOUT, SUBPROCESS.INPFILE );
      fflush(stdout);
    }
  }

  // open OUTFILE in write mode
  SUBPROCESS.FP_OUT = fopen(SUBPROCESS.OUTFILE, "wt");
  if ( !SUBPROCESS.FP_OUT ) {
    sprintf(c1err,"Could not open output file to write:" );
    sprintf(c2err," '%s' ", SUBPROCESS.OUTFILE) ;
    SUBPROCESS_REMIND_STDOUT();
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }
  else {
    printf("%s  Opened output file (fit info): %s\n", 
	   KEYNAME_SUBPROCESS_STDOUT, SUBPROCESS.OUTFILE );
    fflush(stdout);
  }

  // open SALT2mu-LOGFILE in write mode
  FP_STDOUT = fopen(SUBPROCESS.STDOUT_FILE, "wt");
  if ( !FP_STDOUT ) {
    sprintf(c1err,"Could not open STDOUT file to write:" );
    sprintf(c2err," '%s' ", SUBPROCESS.STDOUT_FILE) ;
    SUBPROCESS_REMIND_STDOUT();
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }
  else {
    printf("%s  Opened STDOUT file (stdout): %s\n", 
	   KEYNAME_SUBPROCESS_STDOUT, SUBPROCESS.STDOUT_FILE );
    fflush(stdout);
  }

  // prepare optional dumps
  SUBPROCESS_INIT_DUMP();

  printf("\n");

  for(itable=0; itable < SUBPROCESS.N_OUTPUT_TABLE; itable++ )
    { SUBPROCESS_OUTPUT_TABLE_PREP(itable) ; }
    // debugexit(fnam);
  
  // prep flat random for each event
  SUBPROCESS_INIT_RANFLAT(-1);

  // malloc logicals to keep/reject based on random re-wgt
  SUBPROCESS.KEEP_AFTER_REWGT = (bool*) malloc( NSN_DATA* sizeof(bool) );

  if ( !ISDATA_REAL ) { 
    SUBPROCESS_STORE_EBV(); 
    SUBPROCESS_READ_SIMREF_INPUTS();
  } 

  printf("%s  Finished %s\n", KEYNAME_SUBPROCESS_STDOUT, fnam );
  fflush(stdout);

  //  debugexit(fnam);
  return ;

} // end SUBPROCESS_INIT


// ==========================================================
void SUBPROCESS_READPREP_TABLEVAR(int IFILE, int ISTART, int LEN, 
				  TABLEVAR_DEF *TABLEVAR) {

  // Jul 2020
  // Called during read_data() [before SUBPROCESS_INIT] to 
  // read GENPDF varnames from FITRES file. This function
  //  + parses command-line input SUBPROCESS_VARNAMES_GENPDF to get list
  //    of GENPDF varNames and load SUBPROCESS.VARNAMES_GENPDF[ivar].
  //  + malloc arrays for GENPDF VARNAMES
  //  + call SNTABLE_READPREP_VARDEF to prep table read.
  // 
  // GENPDF varnames are read & stored here regardless of whether
  // they have already been read for SALT2mu fit. Reading is done
  // only for sim-data, so doesn't take much extra memory if a few
  // duplicate colummns are read.
  //
  // Apr 18 2022: set better abort trap for NVAR_ALL

  int  EVENT_TYPE       = TABLEVAR->EVENT_TYPE ;
  char *VARNAMES_STRING = SUBPROCESS.INPUT_VARNAMES_GENPDF_STRING ; 
  int  LEN_MALLOC       = TABLEVAR->LEN_MALLOC ;
  int  debug_malloc     = INPUTS.debug_malloc ;
  int  MEMF             = LEN_MALLOC*sizeof(float) ;
  int  N_TABLE          = SUBPROCESS.N_OUTPUT_TABLE ;
  
  char *ptrVarAll[MXVAR_GENPDF], *varName, varCast[60] ;
  char *VARLIST_READ = (char*) malloc(100*sizeof(char));
  char VARNAME[60], *OUTPUT_TABLE;
  int  VBOSE  = 3;         // print each var; abort on missing var
  int  ivar, ivar2, IVAR_TABLE, NVAR_GENPDF, NVAR_ALL, i, itab, MEM=0 ;
  bool SKIP, MATCH ;
  char fnam[] = "SUBPROCESS_READPREP_TABLEVAR" ;

  // ---------- BEGIN -----------

  // start by adding optional HOST_XXX TABLEVAR columns for both data and sim.
  // logmass is already read by default, so ignore logmass here.
#define  NVAR_HOST 3
  char VARNAMES_HOST[NVAR_HOST][40] = 
    { VARNAME_LOGSFR, VARNAME_LOGsSFR, VARNAME_COLOR } ;

  for(i=0; i < NVAR_HOST; i++ ) {
    sprintf(VARNAME,"%s", VARNAMES_HOST[i] );
    MATCH = false;

    // check if VARNAME appears in any output table
    for(itab=0; itab < N_TABLE; itab++ ) {
      OUTPUT_TABLE = SUBPROCESS.INPUT_OUTPUT_TABLE[itab];
      if ( strstr(OUTPUT_TABLE,VARNAME) != NULL )  { MATCH = true; }
    }

    // if varname appears, read it
    if ( MATCH ) {
      MEM = malloc_TABLEVAR_HOST(LEN_MALLOC,TABLEVAR,VARNAME);
      ivar = SNTABLE_READPREP_HOST(VARNAME, ISTART, LEN, TABLEVAR);
      if ( ivar < 0 ) {
	sprintf(c1err,"Output table includes VARNAME=%s", VARNAME);
	sprintf(c2err,"but %s is not in %s FITRES file", 
		VARNAME, STRING_EVENT_TYPE[EVENT_TYPE] );
	SUBPROCESS_REMIND_STDOUT();
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
      }
    } // end MATCH

  } // end ivar loop over HOST_xxx columns

  if ( TABLEVAR->IS_DATA ) { return; } // return on REAL data



  // - - - - - - - - - 
  // Below is for SIM data only
  if ( EVENT_TYPE != EVENT_TYPE_DATA ) { return; }

  SUBPROCESS.NVAR_GENPDF =  0;
  SUBPROCESS.IVAR_EBV    = -9;
  if ( strlen(VARNAMES_STRING) == 0 ) { return; }

  if ( ISDATA_REAL ) {
    sprintf(c1err,"Woah! Cannot process real data here.");
    sprintf(c2err,"Only SIM data allowed here.");
    SUBPROCESS_REMIND_STDOUT();
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  fprintf(FP_STDOUT, "%s: \n", fnam); fflush(FP_STDOUT);

  print_debug_malloc(+1*debug_malloc,fnam);
  for(ivar=0; ivar < MXVAR_GENPDF; ivar++ ) 
    { ptrVarAll[ivar] = (char*)malloc(MXCHAR_VARNAME*sizeof(char) ); }

  splitString(VARNAMES_STRING, COMMA, MXVAR_GENPDF,    // inputs
	      &NVAR_ALL, ptrVarAll );                  // outputs

  // - - - -- 
  // if SIM_EBV is on list, add AV and RV. Don't worry about
  //  duplicates since duplicates are checked below.
  int ivar_ebv_tmp    = -9 ;
  SUBPROCESS.IVAR_EBV = -9 ;
  for(ivar=0; ivar < NVAR_ALL; ivar++ ) {
    varName = ptrVarAll[ivar] ; 
    if ( strcmp(varName,VARNAME_SIM_EBV) == 0 ) {
      ivar_ebv_tmp = ivar;
    }
  }

  if ( ivar_ebv_tmp >= 0 ) {
    int NVAR_ADD = 2;
    if ( NVAR_ALL > MXVAR_GENPDF-NVAR_ADD ) {

      print_preAbort_banner(fnam);
      for(ivar=0; ivar < NVAR_ALL; ivar++ ) {
	printf("\t Host var[%2d] = %s\n", ivar, ptrVarAll[ivar]);
	fflush(stdout);
      }
      sprintf(c1err,"Cannot add %s and %s because", 
	      VARNAME_SIM_AV, VARNAME_SIM_RV);
      sprintf(c2err,"NVAR_ALL+%d=%d and MXVAR_GENPDF=%d",
	      NVAR_ADD, NVAR_ALL+NVAR_ADD, MXVAR_GENPDF);
      SUBPROCESS_REMIND_STDOUT();
      errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);      
    }
    sprintf(ptrVarAll[NVAR_ALL], "%s", VARNAME_SIM_AV); NVAR_ALL++ ;
    sprintf(ptrVarAll[NVAR_ALL], "%s", VARNAME_SIM_RV); NVAR_ALL++ ;
  }
  

  // Store each FITRES column used by GENPDF maps.
  // Don't bother checking if already read for SALT2mu, but avoid
  // reading duplicates passed by SUBPROCESS_VARNAMES_GENPDF arg.
  // NVAR_ALL is number of ALL variables in SUBPROCESS_VARNAMES_GENPDF;
  // NVAR_GENPDF is number of unique variables after removing
  // duplicates.

  NVAR_GENPDF = VARLIST_READ[0] = 0 ;
  for(ivar=0; ivar < NVAR_ALL; ivar++ ) {
    varName = ptrVarAll[ivar] ;

    sprintf(varCast, "%s:F", varName);
    
    // skip  if duplicate
    SKIP = false;
    if ( ivar > 0 ) {
      for(ivar2=0; ivar2 < ivar; ivar2++ ) {
	if (strcmp(varName,ptrVarAll[ivar2])== 0 ) { SKIP = true; }
      }
    }
    if ( SKIP ) { continue; }

    // malloc on IFILE==0
    if ( IFILE == 0 ) {
      SUBPROCESS.TABLEVAR[NVAR_GENPDF] = (float*)malloc(MEMF); 
      sprintf(SUBPROCESS.VARNAMES_GENPDF[NVAR_GENPDF],"%s", varName);
      catVarList_with_comma(VARLIST_READ,varName);
    }
    
    if ( ivar_ebv_tmp == ivar )   
      { SUBPROCESS.IVAR_EBV = NVAR_GENPDF;  NVAR_GENPDF++; continue; }

    IVAR_TABLE = 
      SNTABLE_READPREP_VARDEF(varCast, 
			      &SUBPROCESS.TABLEVAR[NVAR_GENPDF][ISTART],
			      LEN, VBOSE );   
    NVAR_GENPDF++ ;
  }

  SUBPROCESS.NVAR_GENPDF = NVAR_GENPDF;
  
  // summary for IFILE==0
  if ( IFILE == 0 ) {
    printf("%s loaded %s\n", KEYNAME_SUBPROCESS_STDOUT, VARLIST_READ );
    fflush(stdout);
  }

  return;

} // end SUBPROCESS_READPREP_TABLEVAR

// ========================================
void  SUBPROCESS_INIT_RANFLAT(int iter) {

  // iter < 0 means we need to init malloc before fitting
  // iter > 0 is during fit
  // Init fixed random number [0-1] for each event;
  // used later to select weighted sub-sample of sim-data.

  int NSN = INFO_DATA.TABLEVAR.NSN_ALL ;
  int MEMF = NSN * sizeof(float) ;
  int debug_malloc = INPUTS.debug_malloc ;
  int isn ;
  int OPTMASK = SUBPROCESS.INPUT_OPTMASK ; 
  bool INITSTEP = (iter < 0); 
  double r1, r2;
  char fnam[] = "SUBPROCESS_INIT_RANFLAT" ;

  // ------------ BEGIN -------------

  if (INITSTEP) {
    printf("%s  init randoms with ISEED=%d \n",
	   KEYNAME_SUBPROCESS_STDOUT, SUBPROCESS.INPUT_ISEED );

    init_random_seed(SUBPROCESS.INPUT_ISEED,1);

    print_debug_malloc(+1*debug_malloc,fnam);
    SUBPROCESS.RANFLAT_REWGT    = (float*) malloc ( MEMF );
    SUBPROCESS.RANFLAT_PRESCALE = (float*) malloc ( MEMF );

  }

  bool DORAN = (INITSTEP || (OPTMASK & SUBPROCESS_OPTMASK_RANSEED) > 0) ;
  if (!DORAN) { return; } 

  for(isn=0; isn < NSN; isn++ ) {
    r1 = unix_getRan_Flat1(1);
    r2 = unix_getRan_Flat1(1);

    SUBPROCESS.RANFLAT_REWGT[isn] = (float)r1; 
    SUBPROCESS.RANFLAT_PRESCALE[isn] = (float)r2; 
    if ( isn < 4 )   { 
      printf("%s\t RANFLAT[%d] = %8.5f %8.5f\n", 
	     KEYNAME_SUBPROCESS_STDOUT,isn, r1, r2); 
    }
  }

  return ;

} // end  SUBPROCESS_INIT_RANFLAT

// =======================================
void SUBPROCESS_READ_SIMREF_INPUTS(void) { 

  // read reference sim-inputs to get bounding function
  // for population parameters.
  //
  // Jan 20 2022: increase size of ptr_ITEMLIST to avoid memory corruption

  char *input_simref_file             = SUBPROCESS.INPUT_SIMREF_FILE ;
  GENGAUSS_ASYM_DEF *GENGAUSS_SALT2x1 = &SUBPROCESS.GENGAUSS_SALT2x1;
  GENGAUSS_ASYM_DEF *GENGAUSS_SALT2c  = &SUBPROCESS.GENGAUSS_SALT2c;
  GENGAUSS_ASYM_DEF *GENGAUSS_RV      = &SUBPROCESS.GENGAUSS_RV;
  GEN_EXP_HALFGAUSS_DEF *EXP_EBV      = &SUBPROCESS.EXP_HALFGAUSS_EBV ; 
  GENGAUSS_ASYM_DEF *GENGAUSS_SALT2BETA  = &SUBPROCESS.GENGAUSS_SALT2BETA;
  GENGAUSS_ASYM_DEF *GENGAUSS_SALT2ALPHA  = &SUBPROCESS.GENGAUSS_SALT2ALPHA;

  FILE *finp ; 
  int GZIPFLAG, NITEM, i, NWORD ;
  bool is_salt2, is_rv, is_ebv ; 
  char c_get[MXCHAR_FILENAME], **ptr_ITEMLIST;
  char LINE[MXPATHLEN], TMPLINE[MXPATHLEN] ; 
  char varlist[40] = "";
  char fnam[] = "SUBPROCESS_READ_SIMREF_INPUTS" ; 

  // ---------- BEGIN -------------

  // init optional profiles for REMREF bounding function
  init_GENGAUSS_ASYM(GENGAUSS_SALT2c,  0.0 );
  init_GENGAUSS_ASYM(GENGAUSS_SALT2x1, 0.0 );
  init_GENGAUSS_ASYM(GENGAUSS_RV,      0.0 );
  init_GEN_EXP_HALFGAUSS(EXP_EBV, 0.0) ;
  init_GENGAUSS_ASYM(GENGAUSS_SALT2BETA, 0.0);
  init_GENGAUSS_ASYM(GENGAUSS_SALT2ALPHA, 0.0);

  if (IGNOREFILE(input_simref_file)) {
    SUBPROCESS.ISFLAT_SIM = true ;
    return ;
  } 

  print_banner(fnam);

  SUBPROCESS.ISFLAT_SIM = false ; 

  finp  = open_TEXTgz(input_simref_file, "rt", &GZIPFLAG);
  if (!finp) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Could not open input simref file:");
    sprintf(c2err,"%s", input_simref_file);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }
  else {
    printf("\t Read bounding functions from SIMREF input file:\n"
	    "\t %s\n", input_simref_file );
    fflush(stdout);
  }

  ptr_ITEMLIST = (char**)malloc( 30*sizeof(char*));
  for(i=0; i<30; i++) 
    { ptr_ITEMLIST[i] = (char*)malloc(MXPATHLEN*sizeof(char)); }


  while (fscanf(finp, "%s", c_get) != EOF ) {
    is_salt2  = ( strstr(c_get,"SALT2")  != NULL ) ; 
    is_rv      = ( strstr(c_get,"RV")      != NULL ) ;
    is_ebv     = ( strstr(c_get,"EBV")     != NULL ) ;

    // will need to add EBV later
    if ( is_salt2 || is_rv || is_ebv ) {  // SALT2 or RV is in c_get

      fgets(LINE,MXPATHLEN,finp); 
      sprintf(TMPLINE,"%s %s", c_get, LINE);
      splitString(TMPLINE, " ", MXPATHLEN,          // inputs
                  &NITEM, ptr_ITEMLIST );      // outputs

      if ( is_salt2 ) {
	NWORD = parse_input_GENGAUSS("SALT2c", ptr_ITEMLIST, KEYSOURCE_FILE,
				     GENGAUSS_SALT2c);
	NWORD = parse_input_GENGAUSS("SALT2x1", ptr_ITEMLIST, KEYSOURCE_FILE,
				     GENGAUSS_SALT2x1);
        NWORD = parse_input_GENGAUSS("SALT2BETA", ptr_ITEMLIST, KEYSOURCE_FILE,
                                     GENGAUSS_SALT2BETA);
        NWORD = parse_input_GENGAUSS("SALT2ALPHA", ptr_ITEMLIST, KEYSOURCE_FILE,
                                     GENGAUSS_SALT2ALPHA);
      }

      if ( is_rv ) {
	NWORD = parse_input_GENGAUSS("RV", ptr_ITEMLIST, KEYSOURCE_FILE,
				     GENGAUSS_RV);
      }

      if ( is_ebv ) {
	NWORD = parse_input_EXP_HALFGAUSS("EBV", ptr_ITEMLIST, KEYSOURCE_FILE,
					  EXP_EBV);
        NWORD = parse_input_EXP_HALFGAUSS("EBV_HOST", ptr_ITEMLIST, KEYSOURCE_FILE,
                                          EXP_EBV);

      }

    } // end if is_xxx 
  } // end while

  fclose(finp) ; 

  if ( GENGAUSS_SALT2x1->USE ) { catVarList_with_comma(varlist,"SALT2x1"); }
  if ( GENGAUSS_SALT2c->USE  ) { catVarList_with_comma(varlist,"SALT2c");  }
  if ( GENGAUSS_RV->USE      ) { catVarList_with_comma(varlist,"RV");      }
  if ( EXP_EBV->USE          ) { catVarList_with_comma(varlist,"EBV") ;    }
  if ( GENGAUSS_SALT2BETA->USE ) { catVarList_with_comma(varlist,"SALT2BETA"); }
  if ( GENGAUSS_SALT2ALPHA->USE ) { catVarList_with_comma(varlist,"SALT2ALPHA"); }

  printf("\t Stored bounding functions for '%s' \n", varlist);
  fflush(stdout);


  int LDMP = 0;
  if ( LDMP ) {
    dump_GENGAUSS_ASYM(GENGAUSS_SALT2c);
    dump_GENGAUSS_ASYM(GENGAUSS_SALT2x1);
    dump_GENGAUSS_ASYM(GENGAUSS_RV);
    // xxx maybe will need EBV dump
    dump_GENGAUSS_ASYM(GENGAUSS_SALT2BETA);
    dump_GENGAUSS_ASYM(GENGAUSS_SALT2ALPHA);
  }

  //debugexit(fnam) ; 
  return ; 
} //END SUBPROCESS_READ_SIMREF_INPUTS


// ========================================
void SUBPROCESS_STORE_EBV(void) {

  // Jun 2021
  // if SIM_EBV is requested in SUBPROCESS_VARNAMES_GENPDF,
  // compute EBV array as if it were in the FITRES file.

  int NVAR_GENPDF = SUBPROCESS.NVAR_GENPDF ;
  int IVAR_EBV    = SUBPROCESS.IVAR_EBV ;
  int NSN         = INFO_DATA.TABLEVAR.NSN_ALL ;
  int isn ;
  char  *varName ;
  float AV, RV;
  char fnam[] = "SUBPROCESS_STORE_EBV" ;

  // ---------- BEGIN -----------

  if ( IVAR_EBV < 0 ) { return; }

  printf("%s: compute EBV columm\n", KEYNAME_SUBPROCESS_STDOUT);

  // find AV & RV
  int ivar, ivar_av=-9, ivar_rv=-9;
  for(ivar=0; ivar < NVAR_GENPDF; ivar++ ) {
    varName = SUBPROCESS.VARNAMES_GENPDF[ivar] ;
    if ( strcmp(varName,VARNAME_SIM_AV) == 0 ) { ivar_av = ivar; }
    if ( strcmp(varName,VARNAME_SIM_RV) == 0 ) { ivar_rv = ivar; }
  }

  if ( ivar_av < 0 || ivar_rv < 0 ) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Invalid ivar_[av,rv] = %d, %d", ivar_av, ivar_rv);
    sprintf(c2err,"Cannot compute EBV column");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  for(isn=0; isn < NSN; isn++ ) {
    AV = SUBPROCESS.TABLEVAR[ivar_av][isn] ;
    RV = SUBPROCESS.TABLEVAR[ivar_rv][isn] ;
    SUBPROCESS.TABLEVAR[IVAR_EBV][isn] = AV/RV;
  }

  return ;

} // end SUBPROCESS_STORE_EBV

// ========================================
void SUBPROCESS_PREP_NEXTITER(void) {

  // Created July 2020
  // For real data, do nothing.
  // For sim, prepare for next iteration.

  int  ITER_EXPECT = -9 ;
  char fnam[] = "SUBPROCESS_PREP_NEXTITER";

  // --------- BEGIN -------------

  if ( ISDATA_REAL ) { return ; }

  // request expected iteration number to be entered as std input
  // or q to quit

  printf("\n%s Enter expected ITERATION number (-1 to quit) => \n",
	 KEYNAME_SUBPROCESS_STDOUT );   fflush(stdout);
  scanf( "%d", &ITER_EXPECT); // read response
  if ( ITER_EXPECT < 0 ) { SUBPROCESS_EXIT(); }
  SUBPROCESS.ITER = ITER_EXPECT ; 

  prep_input_repeat();

  // rewind all SUBPROCESS files
  rewind(SUBPROCESS.FP_INP);   
  rewind(SUBPROCESS.FP_OUT);   
  if ( SUBPROCESS.STDOUT_CLOBBER ) { rewind(FP_STDOUT); }

  // - - - - - -
  
  SUBPROCESS_SIM_REWGT(ITER_EXPECT);

  //  debugexit(fnam);

  return ;

} // end void SUBPROCESS_PREP_NEXTITER


// ========================================
void SUBPROCESS_SIM_REWGT(int ITER_EXPECT) {

  // Created July 2020
  // For real data, do nothing.
  // For sim, read PDF population map(s), same map as for sim-input
  // GENPDF_FILE, and reweight sim data assuming that sim was 
  // generated with a flat distribution in each variable.
  // 
  // July 13 2021
  // Updated to include bounding function option

  int  OPTMASK  = OPTMASK_GENPDF_EXTERNAL_FP ;
  int  ITER     = SUBPROCESS.ITER ;
  FILE *FP_INP  = SUBPROCESS.FP_INP ;
  char *INPFILE = SUBPROCESS.INPFILE ;

  bool FOUND_ITER_BEGIN = false;
  char c_get[60];
  int  ISTAT_READ=-9, ITER_FOUND = -9  ;
  char fnam[] = "SUBPROCESS_SIM_REWGT" ;

  // -------- BEGIN -----------

  // read input file until we reach iteration key
  while ( !FOUND_ITER_BEGIN && ISTAT_READ != EOF ) {
    ISTAT_READ = fscanf(FP_INP, "%s", c_get) ;
    if ( strcmp(c_get,KEYNAME_SUBPROCESS_ITERATION_BEGIN) == 0 ) {
      FOUND_ITER_BEGIN = true ;
      fscanf(FP_INP, "%d", &ITER_FOUND);
    }
  }

  // make sure that ITERATION key was found, and that it
  // matches ITER_EXPECT entered on command line.

  if ( ITER_FOUND < 0 ) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Could not find required '%s' key in", 
	    KEYNAME_SUBPROCESS_ITERATION_BEGIN);
    sprintf(c2err,"file %s", INPFILE );
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  if ( ITER_EXPECT != ITER_FOUND ) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Found ITERATION=%d in PDF file %s",
	    ITER_FOUND, INPFILE);
    sprintf(c2err,"But expected ITERATION=%d passed via std input", 
	    ITER_EXPECT);
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  printf("%s Read PDF map(s) for ITERATION=%d\n", 
	 KEYNAME_SUBPROCESS_STDOUT, ITER );

  // re-init uniqueOverlap in case a key is parsed again
  uniqueOverlap(STRINGMATCH_INIT,"SUBPROCESS"); 

  init_genPDF(OPTMASK, FP_INP, INPFILE, "");

  // over-write CUTBIT_SPLITRAN 
  sprintf(CUTSTRING_LIST[CUTBIT_SPLITRAN],  "GENPDF rewgt");

  // - - - - -
  // prepare index map between IVAR(MAP) and IVAR(TABLE)
  int imap, NVAR_GENPDF, ivar, IVAR_TABLE ;
  char *varName_GENPDF ;
  for(imap=0; imap < NMAP_GENPDF; imap++ ) {
    NVAR_GENPDF = GENPDF[imap].GRIDMAP.NDIM; 

    for(ivar=0; ivar < NVAR_GENPDF; ivar++ ) {
      varName_GENPDF = GENPDF[imap].VARNAMES[ivar];
      IVAR_TABLE = SUBPROCESS_IVAR_TABLE(varName_GENPDF);
      if ( IVAR_TABLE < 0 ) {
	sprintf(c1err,"Could not find IVAR_TABLE for GENPDF var = '%s'", 
		varName_GENPDF);
	sprintf(c2err,"Check command-line arg SUBPROCESS_VARNAMES_GENPDF");
	errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
      }

      /*
      int IDMAP = IDGRIDMAP_GENPDF + imap;
      printf(" xxx %s: imap=%d ivar=%d IDMAP=%d --> IVAR_TABLE=%d \n",
	     fnam, imap, ivar, IDMAP, IVAR_TABLE);
      */
      SUBPROCESS.IVAR_TABLE_GENPDF[imap][ivar] = IVAR_TABLE ;

    } // end ivar-GENPDF loop
  }  // end imap loop

  // -------------------------------------

  // loop over sim data and set make to keep/reject based on GENPDF map
  int isn, istat, SIM_NONIA_INDEX, NKEEP_ORIG=0, NKEEP_REWGT=0 ;
  bool LDMP, KEEP ;
  int NSN = INFO_DATA.TABLEVAR.NSN_ALL;
  char *name;
  double XVAL, XVAL_for_GENPDF[MXVAR_GENPDF], PROB, PROB_TOT, RANFLAT, PROB_SIMREF, PROB_RATIO ;
  SUBPROCESS.MAXPROB_RATIO = 0.;

  SUBPROCESS_INIT_RANFLAT(ITER_EXPECT); 

  for ( isn=0 ; isn < NSN; isn++ ) {
    PROB_SIMREF     = 1.0;
    PROB_TOT        = 1.0;
    name            = INFO_DATA.TABLEVAR.name[isn];
    SIM_NONIA_INDEX = INFO_DATA.TABLEVAR.SIM_NONIA_INDEX[isn];
    //    CUTMASK         = INFO_BIASCOR.TABLEVAR.CUTMASK[isn] ;
    SUBPROCESS.KEEP_AFTER_REWGT[isn] = KEEP = false;
    LDMP            = SUBPROCESS.DUMPFLAG_REWGT[isn];

    if ( SIM_NONIA_INDEX > 0 ) { continue; } // reject of not true SNIa
    //    if ( CUTMASK != 0        ) { continue; }

    NKEEP_ORIG++ ;

    if ( LDMP ) {
      printf(" xxx \n");
      printf(" xxx %s -------------------------------------- \n", fnam ); 
      printf(" xxx %s DUMP for isn=%d, SNID=%s \n", fnam, isn, name ); 
      fflush(stdout);
    }

    for(imap=0; imap < NMAP_GENPDF; imap++ ) {
      NVAR_GENPDF = GENPDF[imap].GRIDMAP.NDIM ;

      for(ivar=0; ivar < NVAR_GENPDF; ivar++ ) {
	IVAR_TABLE = SUBPROCESS.IVAR_TABLE_GENPDF[imap][ivar] ;
	XVAL       = SUBPROCESS.TABLEVAR[IVAR_TABLE][isn] ;
	XVAL_for_GENPDF[ivar] = XVAL ;
      } // end ivar loop

      istat = interp_GRIDMAP(&GENPDF[imap].GRIDMAP, XVAL_for_GENPDF, &PROB);
      // check for bounding function here
      if (! SUBPROCESS.ISFLAT_SIM) {
	PROB_SIMREF =  SUBPROCESS_PROB_SIMREF(ITER_FOUND, imap, XVAL_for_GENPDF[0]) ;
      }

      PROB_RATIO = (PROB/ PROB_SIMREF) ;
      PROB_TOT *= PROB_RATIO ;

      if (PROB_RATIO > SUBPROCESS.MAXPROB_RATIO) {
	SUBPROCESS.MAXPROB_RATIO = PROB_RATIO ; 
      }
      // printf("xxx PROB = %le, PROB_SIMREF=%le, PROB_RATIO=%le \n", PROB, PROB_SIMREF, PROB_RATIO) ;  

      if ( LDMP ) {
	XVAL = XVAL_for_GENPDF[0]; 
	printf(" xxx %s: PROB(%s=%.3f) = %f \n", 
		 fnam, GENPDF[imap].VARNAMES[0], XVAL, PROB); 
	fflush(stdout);
      }
      
    } // end GENPDF map loop

    // apply random number against PROB_TOT to keep or reject this event.
    RANFLAT = (double)SUBPROCESS.RANFLAT_REWGT[isn];

    if ( PROB_TOT >= RANFLAT ) { KEEP = true; NKEEP_REWGT++ ; } 
    if ( LDMP ) {
      printf(" xxx %s: PROB_TOT=%7.4f  FLATRAN=%7.4f  KEEP=%d\n", 
	     fnam, PROB_TOT, RANFLAT, KEEP );
      fflush(stdout);
    }

    SUBPROCESS.KEEP_AFTER_REWGT[isn] = KEEP;

    if ( !KEEP ) 
      { setbit_CUTMASK(isn, CUTBIT_SPLITRAN, &INFO_DATA.TABLEVAR ); }

  } // end isn loop

  // - - - - -
  printf("%s  Keep %d of %d events after GENPDF reweight\n", 
	 KEYNAME_SUBPROCESS_STDOUT, NKEEP_REWGT, NKEEP_ORIG );
  fflush(stdout);

  return ;

} // end SUBPROCESS_SIM_REWGT

double SUBPROCESS_PROB_SIMREF(int ITER, int imap, double XVAL) {

  // Created July 2021
  // Return bound function prob from reference sim.

  char fnam[] = "SUBPROCESS_PROB_SIMREF" ; 
  double PROB_SIMREF = 1.0 ; 
  int NVAR = SUBPROCESS.NVAR_GENPDF ;
  char *VARNAME = GENPDF[imap].VARNAMES[0] ; 
  bool  SET_INDEX = true;

  // begin
  if ( SET_INDEX ) { 

    if (strcmp(VARNAME, "SIM_RV") == 0 ) 
      { SUBPROCESS.GENGAUSS_RV.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_c") == 0 ) 
      { SUBPROCESS.GENGAUSS_SALT2c.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_x1") == 0 ) 
      { SUBPROCESS.GENGAUSS_SALT2x1.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_beta") == 0 ) 
      { SUBPROCESS.GENGAUSS_SALT2BETA.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_alpha") == 0 ) 
      { SUBPROCESS.GENGAUSS_SALT2ALPHA.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_EBV") == 0 ) 
      { SUBPROCESS.EXP_HALFGAUSS_EBV.INDEX = imap ; }
    if (strcmp(VARNAME, "SIM_EBV_HOST") == 0 ) 
      { SUBPROCESS.EXP_HALFGAUSS_EBV.INDEX = imap ; }

  } // end SET_INDEX

  // - - - - - - - - 
  if (SUBPROCESS.GENGAUSS_RV.INDEX == imap) {
    PROB_SIMREF = funVal_GENGAUSS_ASYM(XVAL ,&SUBPROCESS.GENGAUSS_RV) ; 
  }
 
  else if (SUBPROCESS.GENGAUSS_SALT2c.INDEX == imap) {
    PROB_SIMREF = funVal_GENGAUSS_ASYM(XVAL,&SUBPROCESS.GENGAUSS_SALT2c) ; 
  }
  
  else if (SUBPROCESS.GENGAUSS_SALT2x1.INDEX == imap) {
    PROB_SIMREF = funVal_GENGAUSS_ASYM(XVAL,&SUBPROCESS.GENGAUSS_SALT2x1) ;
  }

  else if (SUBPROCESS.GENGAUSS_SALT2BETA.INDEX == imap) {
    PROB_SIMREF = funVal_GENGAUSS_ASYM(XVAL,&SUBPROCESS.GENGAUSS_SALT2BETA);
  }
  else if (SUBPROCESS.GENGAUSS_SALT2ALPHA.INDEX == imap) {
    PROB_SIMREF = funVal_GENGAUSS_ASYM(XVAL,&SUBPROCESS.GENGAUSS_SALT2ALPHA);
  }

  else if (SUBPROCESS.EXP_HALFGAUSS_EBV.INDEX == imap) {
    PROB_SIMREF = funVal_GEN_EXP_HALFGAUSS(XVAL,&SUBPROCESS.EXP_HALFGAUSS_EBV);
  }
  
  else {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Did not find bounding functions for '%s'", VARNAME) ;
    sprintf(c2err,"Bounding func must be specified for all variables or none!") ;
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  if ( PROB_SIMREF < 0.) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"Invalid PROB_SIMREF=%f for VARNAME=%s", PROB_SIMREF, VARNAME) ;
    sprintf(c2err,"ITER = %d, imap=%d, XVAL=%f", ITER, imap, XVAL) ;
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);  
  }

  if (ITER < -1) {
    printf("xxx ITER = %d, imap = %d, PROB_SIMREF = %f, XVAL = %f, VARNAME = %s \n", ITER, imap, PROB_SIMREF, XVAL, VARNAME); 
  }

  return PROB_SIMREF; 
} // end SUBPROCESS_PROB_SIMREF


// ================================
void SUBPROCESS_SIM_PRESCALE(void) {

  // Created Jun 7 2021
  // if SUBPROCESS.NEVT_SIM_PRESCALE > 0, figure out pre-scale
  // to get this NEVT, and apply it. This allows control code
  // to determine the sim subset size for BBC fit.
  
  bool IS_SIM     = INFO_DATA.TABLEVAR.IS_SIM ;
  int  NEVT_SIM   = SUBPROCESS.NEVT_SIM_PRESCALE ;
  int  NSN_TOT    = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  isn;
  int  LDMP = 1;
  double RANFLAT ;
  char fnam[] = "SUBPROCESS_SIM_PRESCALE" ;

  // ---------- BEGIN -----------

  if ( !IS_SIM || NEVT_SIM < 10 ) { return; }

  print_eventStats(EVENT_TYPE_DATA);

  int  NSN_PASS   = *NPASS_CUTMASK_POINTER[EVENT_TYPE_DATA]; 
  if ( NEVT_SIM > NSN_PASS ) { return; }

  // compute prescale (nearest integer) to get close to NEVT_SIM
  float PS = (float)NSN_PASS / (float)NEVT_SIM ;
  float PSinv = 1.0/PS;

  fprintf(FP_STDOUT,"%s: force prescale=%.1f to get NEVT_SIM ~ %d \n",
	  fnam, PS, NEVT_SIM) ;
  fflush(FP_STDOUT);

  for(isn=0; isn < NSN_TOT; isn++) {
    RANFLAT = SUBPROCESS.RANFLAT_PRESCALE[isn];
    if ( RANFLAT > PSinv ) 
      { setbit_CUTMASK(isn, CUTBIT_SIMPS, &INFO_DATA.TABLEVAR ); }
  }

  return;

} // end SUBPROCESS_SIM_PRESCALE

// =======================================
int SUBPROCESS_IVAR_TABLE(char *varName) {
  int ivar;
  int IVAR_TABLE = -9;
  char *varTmp;
  for(ivar=0; ivar < SUBPROCESS.NVAR_GENPDF; ivar++ ) {
    varTmp = SUBPROCESS.VARNAMES_GENPDF[ivar];
    if ( strcmp(varName,varTmp)==0 ) { IVAR_TABLE = ivar; }
  }


  return(IVAR_TABLE);
} // end SUBPROCESS_IVAR_TABLE

// =======================================
void SUBPROCESS_INIT_DUMP(void) {

  // parse command-line input SUBPROCESS_CID_REWGT_DUMP
  // and load DUMPFLAG_REWGT make for each isn index.
  
  int NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL ;
  int MXSPLIT=20, NSPLIT=0, isn, i, SNID ;
  int debug_malloc = INPUTS.debug_malloc ;
  bool MATCH, PICK_isn;
  char *ptrSNID[20], *name ;
  char *string = SUBPROCESS.INPUT_CID_REWGT_DUMP ;
  char fnam[] = "SUBPROCESS_INIT_DUMP" ;

  // ------------- BEGIN ---------------

  print_debug_malloc(+1*debug_malloc,fnam);
  SUBPROCESS.DUMPFLAG_REWGT = (bool*) malloc( NSN_DATA* sizeof(bool) );

  if ( strlen(string) > 0 ) { 
    for(i=0; i < MXSPLIT; i++ ) 
      { ptrSNID[i] = (char*) malloc( 20*sizeof(char) ); }
    
    splitString(string, COMMA, MXSPLIT,    // inputs
		&NSPLIT, ptrSNID );        // outputs
  }

  for(isn=0; isn < NSN_DATA; isn++ ) {
    SUBPROCESS.DUMPFLAG_REWGT[isn] = false;
    name = INFO_DATA.TABLEVAR.name[isn];

    // check user-input list for SNIDs to dump
    for(i=0; i < NSPLIT; i++ ) {
      sscanf(ptrSNID[i], "%d", &SNID);
      MATCH = strcmp(name,ptrSNID[i]) == 0 ;
      PICK_isn = (SNID == isn && isn < 10);
      if ( MATCH || PICK_isn ) { SUBPROCESS.DUMPFLAG_REWGT[isn] = true; }
    }
  }  

  return;
} // end SUBPROCESS_INIT_DUMP

// =======================================
void SUBPROCESS_OUTPUT_TABLE_PREP(int itable) {

  // Sep 17 2020
  // prep output tables.
  // Jun 1 2021: split variables by % or * because % wreaks havoc with python.

  // strup string from user input; e.g, 'x1(10,-4:4)%c(6,-0.3:0.3)'
  char *TABLE_STRING = SUBPROCESS.INPUT_OUTPUT_TABLE[itable]; 
  int  MXVAR         = MXVAR_TABLE_SUBPROCESS ;
  int  NVAR, ivar ;
  int debug_malloc = INPUTS.debug_malloc ;
  char *ptrVarDef[MXVAR], DELIM[2]; 
  char BININFO_STRING[40];
  char fnam[] = "SUBPROCESS_OUTPUT_TABLE_PREP" ;

  // ----------- BEGIN -----------

  print_debug_malloc(+1*debug_malloc,fnam);

  SUBPROCESS.OUTPUT_TABLE[itable].NVAR     = 0 ;

  // first split by "%" or "*" to get each variable/dimension
  for(ivar=0; ivar < MXVAR; ivar++ ) 
    { ptrVarDef[ivar] = (char*) malloc( 60*sizeof(char) ); }
  
  if ( strstr(TABLE_STRING,PERCENT) != NULL ) 
    { sprintf(DELIM,"%s", PERCENT); }
  else
    { sprintf(DELIM,"%s", STAR); }
  splitString(TABLE_STRING, DELIM, MXVAR,       // inputs
	      &NVAR, ptrVarDef );               // outputs

  // - - - - 
  // Now we have [VARNAME]([nbin],[min]:[max])
  // so extract varname and bin-info from ()
  for(ivar=0; ivar < NVAR; ivar++ ) {
    SUBPROCESS_STORE_BININFO(itable, ivar, ptrVarDef[ivar] ) ;
  }
  
  // convert N-D tables into 1D arrays for each access later.
  SUBPROCESS_MAP1D_BININFO(itable);

  // construct VARNAMES list for table header
  SUBPROCESS_OUTPUT_TABLE_HEADER(itable);

  // - - - - - 

  print_debug_malloc(-1*debug_malloc,fnam);
  for(ivar=0; ivar < MXVAR; ivar++ ) 
    { free(ptrVarDef[ivar]);  }


} // end SUBPROCESS_OUTPUT_TABLE_PREP


// ==============================================
void SUBPROCESS_STORE_BININFO(int ITABLE, int IVAR, char *VARDEF_STRING ) {

  // Input VARDEF_STRING is of the form
  //    'x1(10,-4:4)'
  //
  // Parse VARDEF_STRING and load info into global 
  //     SUBPROCESS.OUTPUT_TABLE[ITABLE][IVAR]
  // with all info related to this VARDEF.
  //
  // Mar 7 2022: 
  //  fix to check non-standard BBC host variables (e.g. HOST_COLOR)

  int debug_malloc = INPUTS.debug_malloc ;
  bool LDMP = false ;
  int    NSPLIT, nbin, i, IVAR_FITRES ;
  double xmin, xmax, lo, hi, binSize ;
  char VARNAME[40], stringOpt[40], *ptrSplit[2], *ptrRange[2] ;
  char fnam[] = "SUBPROCESS_STORE_BININFO" ;

  // ------------ BEGIN ------------

  // extract varname out (), and bin info inside ()
  sprintf(VARNAME,"%s", VARDEF_STRING);
  extractStringOpt(VARNAME, stringOpt); // return stringOpt

  print_debug_malloc(+1*debug_malloc,fnam);
  for(i=0;  i < 2; i++ ) {
    ptrSplit[i]  = (char*)malloc(40*sizeof(char) ) ;
    ptrRange[i]  = (char*)malloc(40*sizeof(char) ) ;
  }

  // split by comma to get nbin and xmin:xmax
  splitString(stringOpt, COMMA, 2,       // inputs
	      &NSPLIT, ptrSplit );        // outputs

  sscanf(ptrSplit[0], "%d", &nbin);

  if ( nbin >= MXz ) {
    SUBPROCESS_REMIND_STDOUT();
    sprintf(c1err,"NBIN(%s) = %d exceeds bound of %d", VARNAME, nbin, MXz) ;
    sprintf(c2err,"Check SUBPROCESS_OUTPUT_TABLE args") ;
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  // split xmin:xmax by colon
  splitString(ptrSplit[1], COLON, 2,       // inputs
	      &NSPLIT, ptrRange );        // outputs

  sscanf(ptrRange[0], "%le", &xmin) ;
  sscanf(ptrRange[1], "%le", &xmax) ;
  binSize = (xmax - xmin)/(double)nbin;
  
  // - - - - 
  printf("%s  TABLE-%d  VARNAME = %s  NBIN=%d  RANGE=%.3f to %.3f \n",
	 KEYNAME_SUBPROCESS_STDOUT, ITABLE, VARNAME, nbin, xmin, xmax );
  fflush(stdout);

 
  // - - - - - - - - - - - 
  // load global info  
  SUBPROCESS.OUTPUT_TABLE[ITABLE].NVAR++ ;
  SUBPROCESS.OUTPUT_TABLE[ITABLE].IVAR_FITRES[IVAR] = -9; //IVAR_FITRES ;

  BININFO_DEF *BININFO = &SUBPROCESS.OUTPUT_TABLE[ITABLE].BININFO[IVAR];
  sprintf(BININFO->varName, "%s", VARNAME);
  BININFO->nbin    = nbin ;
  BININFO->binSize = binSize ;

  for(i=0; i < nbin; i++ ) {
    lo = xmin + binSize * (double)i ;
    hi = lo + binSize;
    BININFO->lo[i]  = lo ;
    BININFO->hi[i]  = hi ;
    BININFO->avg[i] = 0.5*(lo+hi) ;
    BININFO->n_perbin[i] = 0;
  }


  // - - - - - - - - - - - - - - -
  // Assign float pointer to INFO_DATA.TABLEVAR array
  // For now it's hard-wired, but later should match based
  // on column name and column index.

  float *PTRVAL ; 
  int IVAR_TABLE = SUBPROCESS_IVAR_TABLE(VARNAME);
  int ivar, NVAR_VALID=0 ;
  bool MATCH     = false;

  // define space separated list of valid varnames for output table
  char VARNAME_VALID_LIST[] = "x1 c zHD "				\
    "HOST_"HOSTGAL_PROPERTY_BASENAME_LOGMASS " " HOSTGAL_PROPERTY_BASENAME_LOGMASS " "\
    "HOST_"HOSTGAL_PROPERTY_BASENAME_LOGSFR  " " HOSTGAL_PROPERTY_BASENAME_LOGSFR  " "\
    "HOST_"HOSTGAL_PROPERTY_BASENAME_LOGsSFR " " HOSTGAL_PROPERTY_BASENAME_LOGsSFR " "\
    "HOST_"HOSTGAL_PROPERTY_BASENAME_COLOR   " " HOSTGAL_PROPERTY_BASENAME_COLOR   " "\
      ; 

  // define list of pointers corresponding to VARNAME_VALID_LIST
  float *PTRVAL_VALID_LIST[] = {
    INFO_DATA.TABLEVAR.fitpar[INDEX_x1],
    INFO_DATA.TABLEVAR.fitpar[INDEX_c],
    INFO_DATA.TABLEVAR.zhd,
    INFO_DATA.TABLEVAR.host_logmass, INFO_DATA.TABLEVAR.host_logmass,
    INFO_DATA.TABLEVAR.host_logsfr,  INFO_DATA.TABLEVAR.host_logsfr,
    INFO_DATA.TABLEVAR.host_logssfr, INFO_DATA.TABLEVAR.host_logssfr,
    INFO_DATA.TABLEVAR.host_color,   INFO_DATA.TABLEVAR.host_color
  } ;

  if ( IVAR_TABLE >= 0 ) {
    // Mar 7 2022: load value from supplemental table 
    PTRVAL = SUBPROCESS.TABLEVAR[IVAR_TABLE];
    MATCH  = true;
  }
  else {
    //    debugexit(VARNAME_VALID_LIST);
    char VARNAME_VALID_TMP[60];
    NVAR_VALID = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, VARNAME_VALID_LIST);
    for(ivar=0; ivar < NVAR_VALID; ivar++ ) {
      get_PARSE_WORD(0, ivar, VARNAME_VALID_TMP);
      if ( strcmp(VARNAME,VARNAME_VALID_TMP) == 0  ) 
	{ PTRVAL = PTRVAL_VALID_LIST[ivar];  MATCH=true; }
    } 
  }

  // - - - - - -
  if ( !MATCH ) {
    print_preAbort_banner(fnam);    
    SUBPROCESS_REMIND_STDOUT();
    fprintf(FP_STDOUT,"   Valid varnames for output table:\n  %s\n", VARNAME_VALID_LIST);
    sprintf(c1err,"Unknown output table var = '%s'", VARNAME);
    sprintf(c2err,"Check valid varnames above.");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }

  /* xxx mark delete xxxxxx
  // below, check for standard SALT2mu table variables
  else if ( strcmp(VARNAME,"x1") == 0  ) 
    { PTRVAL = INFO_DATA.TABLEVAR.fitpar[INDEX_x1];  }
  else if ( strcmp(VARNAME,"c") == 0  ) 
    { PTRVAL = INFO_DATA.TABLEVAR.fitpar[INDEX_c];   }

  else if ( strcmp(VARNAME,"zhd") == 0     || 
	    strcmp(VARNAME,"zHD") == 0  ) 
    { PTRVAL = INFO_DATA.TABLEVAR.zhd ; }

  else if ( strcmp(VARNAME,"HOST_LOGMASS") == 0  ||
	    strcmp(VARNAME,"LOGMASS") == 0   ) 
    { PTRVAL = INFO_DATA.TABLEVAR.host_logmass;  }
  else if ( strcmp(VARNAME,"HOST_LOGSFR") == 0  ||
	    strcmp(VARNAME,"LOGSFR") == 0   ) 
    { PTRVAL = INFO_DATA.TABLEVAR.host_logsfr;  }
  else if ( strcmp(VARNAME,"HOST_LOGsSFR") == 0  ||
	    strcmp(VARNAME,"LOGsSFR") == 0   ) 
    { PTRVAL = INFO_DATA.TABLEVAR.host_logssfr;  }
  else if ( strcmp(VARNAME,"HOST_COLOR") == 0  ||
	    strcmp(VARNAME,"COLOR") == 0   ) 
    { PTRVAL = INFO_DATA.TABLEVAR.host_color;  }

  else {
    sprintf(c1err,"Unknown output table var = '%s'", VARNAME);
    sprintf(c2err,"Check SUBPROCESS_OUTPUT_TABLE args");
    errlog(FP_STDOUT, SEV_FATAL, fnam, c1err, c2err);
  }
  xxxxxx end mark xxxxx */

  SUBPROCESS.OUTPUT_TABLE[ITABLE].PTRVAL[IVAR] = PTRVAL ;

  // free local memory
  print_debug_malloc(-1*debug_malloc,fnam);
  for(i=0; i < 2; i++ ) { free(ptrSplit[i]);  free(ptrRange[i]); }

  return ;

} // end SUBPROCESS_STORE_BININFO

// ===============================
void SUBPROCESS_MAP1D_BININFO(int ITABLE) {

  // Convert multi-D tables into 1D arrays for easy access.
  // 
  int debug_malloc = INPUTS.debug_malloc ;
  int NVAR = SUBPROCESS.OUTPUT_TABLE[ITABLE].NVAR ;
  int i, ivar, NBINTOT=1, nbin, nbin_per_var[MXVAR_TABLE_SUBPROCESS];
  int MEMI, MEMD, MEMDD, IB1D, ib0, ib1, ib2;
  int ib_per_var[MXVAR_TABLE_SUBPROCESS];
  int  IDMAP = 10 + ITABLE;
  int  *INDEX_BININFO[MXVAR_TABLE_SUBPROCESS];
  char fnam[] = "SUBPROCESS_MAP1D_BININFO"; 

  // ------------ BEGIN ------------

  // store array with nbin per ivar
  for(ivar=0; ivar < MXVAR_TABLE_SUBPROCESS; ivar++ )  { 
    nbin = 1;
    if(ivar<NVAR) 
      { nbin = SUBPROCESS.OUTPUT_TABLE[ITABLE].BININFO[ivar].nbin ; }
    nbin_per_var[ivar] = nbin; 
    NBINTOT *= nbin;
  }

  // utility for N-dim -> 1-Dim map
  init_1DINDEX(IDMAP, NVAR, nbin_per_var);

  // allocate NBINTOT memory for each variable/dimension
  SUBPROCESS.OUTPUT_TABLE[ITABLE].NBINTOT = NBINTOT;
  MEMI = NBINTOT * sizeof(int);
  MEMD = NBINTOT * sizeof(double);
  MEMDD = NBINTOT * sizeof(double*);

  print_debug_malloc(+1*debug_malloc,fnam);
  for(ivar=0; ivar < NVAR; ivar++ ) {

    SUBPROCESS.OUTPUT_TABLE[ITABLE].INDEX_BININFO[ivar] = (int*)malloc(MEMI);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].NEVT          = (int   *)malloc(MEMI);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SQSUM   = (double*)malloc(MEMD);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SUM     = (double*)malloc(MEMD);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SUM_WGT = (double*)malloc(MEMD);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].SUM_WGT       = (double*)malloc(MEMD);

    
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_STD        = (double*)malloc(MEMD);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_STD_ROBUST = (double*)malloc(MEMD);
    SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_LIST = (double**)malloc(MEMDD);

    INDEX_BININFO[ivar] = SUBPROCESS.OUTPUT_TABLE[ITABLE].INDEX_BININFO[ivar];
    for(i=0; i < NBINTOT; i++ )  { INDEX_BININFO[ivar][i] = -9 ; }
  }

  // Clumsy: hard-wire 3D -> 1D index map, even if NVAR<3.
  // Note that there are more elegant methods for arbitrary dimensions.
  for (ib0=0; ib0 < nbin_per_var[0]; ib0++ ) {
    for (ib1=0; ib1 < nbin_per_var[1]; ib1++ ) {
      for (ib2=0; ib2 < nbin_per_var[2]; ib2++ ) {
	ib_per_var[0] = ib0;
	ib_per_var[1] = ib1;
	ib_per_var[2] = ib2;
	IB1D = get_1DINDEX(IDMAP, NVAR, ib_per_var);
	
	//	  printf(" 5. xxx %s ib[0,1,2] = %d, %d, %d IB1D = %d \n", 
	//	 fnam, ib0, ib1, ib2, IB1D ); fflush(stdout);
	
	for(ivar=0; ivar < NVAR; ivar++ )
	  { INDEX_BININFO[ivar][IB1D] = ib_per_var[ivar]; }
	
      } // end ib2
    } // end ib1
  }  // end ib0
  

  return;

} // end SUBPROCESS_MAP1D_BININFO


// ====================================
void SUBPROCESS_OUTPUT_TABLE_HEADER(int ITABLE) { 
  //Brodie Jun 14 2021 Added WGT columns 
  int NVAR = SUBPROCESS.OUTPUT_TABLE[ITABLE].NVAR ;
  int ivar;
  char VARNAMES[200], varName[40] ;
  char VARNAMES_FIX[100]; // = "NEVT MURES_SUM MURES_SQSUM  SUM_WGT MURES_SUM_WGT" ;
  BININFO_DEF *BININFO;

  // ----------- BEGIN ---------

  sprintf(VARNAMES_FIX,"NEVT MURES_SUM STD STD_ROBUST");

  sprintf(VARNAMES,"ROW ");

  for(ivar=0; ivar < NVAR; ivar++ ) {
    BININFO = &SUBPROCESS.OUTPUT_TABLE[ITABLE].BININFO[ivar];
    sprintf(varName,"ibin_%s ", BININFO->varName);
    strcat(VARNAMES,varName);
  }
  strcat(VARNAMES,VARNAMES_FIX);

  sprintf(SUBPROCESS.OUTPUT_TABLE[ITABLE].VARNAMES_HEADER,"%s", VARNAMES);

  return ;
} // SUBPROCESS_TABLE_HEADER


// ===========================================
void SUBPROCESS_OUTPUT_LOAD(void) {

  // driver function to load output tables after subprocess iteration
  int  NSN_DATA      = INFO_DATA.TABLEVAR.NSN_ALL ;
  int  N_TABLE       = SUBPROCESS.N_OUTPUT_TABLE;
  int  isn, ITABLE, cutmask; 
  char *TABLE_NAME ;
  char fnam[] = "SUBPROCESS_OUTPUT_LOAD" ;

  // ---------- BEGIN ----------

  for(ITABLE=0; ITABLE < N_TABLE; ITABLE++ ) {

    SUBPROCESS_OUTPUT_TABLE_RESET(ITABLE);

    for(isn=0; isn < NSN_DATA; isn++ ) {
      cutmask = INFO_DATA.TABLEVAR.CUTMASK[isn]; 
      if ( !keep_cutmask(cutmask)  ) { continue; }

      SUBPROCESS_OUTPUT_TABLE_LOAD(isn,ITABLE);
    }  // end isn

    // RK 9.29.2021 compute and store STD and STD_ROBUST
    SUBPROCESS_COMPUTE_STD(ITABLE); 

  } // end ITABLE
  
  return ;

} // end SUBPROCESS_OUTPUT_LOAD


// ========================================
void SUBPROCESS_COMPUTE_STD(int ITABLE) {

  // Created Sep 29 2021
  // for ITABLE, compute STD and STD_ROBUST for each 1d bin,
  // and store it in SUBPROCESS.OUTPUT_TABLE[ITABLE].

  SUBPROCESS_TABLE_DEF *OUTPUT_TABLE = &SUBPROCESS.OUTPUT_TABLE[ITABLE] ;

  int    NBINTOT, ibin, NEVT, i;
  double SUM, SQSUM, STD, STD_ROBUST, *MURES_LIST;
  double MURES_AVG, *ABSMURES_LIST, DUM_AVG, DUM_STD, MEDIAN ;
  char fnam[] = "SUBPROCESS_COMPUTE_STD" ;

  // -------------- BEGIN ------------

  NBINTOT  = OUTPUT_TABLE->NBINTOT ;
  for(ibin=0; ibin < NBINTOT; ibin++ ) {
    NEVT  = OUTPUT_TABLE->NEVT[ibin];
    SUM   = OUTPUT_TABLE->MURES_SUM[ibin];
    SQSUM = OUTPUT_TABLE->MURES_SQSUM[ibin];
    MURES_LIST = OUTPUT_TABLE->MURES_LIST[ibin];

    if (NEVT <= 1 ) {
      OUTPUT_TABLE->MURES_STD[ibin]        = 0.0 ;
      OUTPUT_TABLE->MURES_STD_ROBUST[ibin] = 0.0 ;
      continue ;
    }
    MURES_AVG     = SUM/(double)NEVT;
    ABSMURES_LIST = (double*) malloc( NEVT*sizeof(double) );
    for(i=0; i < NEVT; i++ )
      { ABSMURES_LIST[i] = fabs(MURES_LIST[i]-MURES_AVG); }
	
    arrayStat(NEVT, ABSMURES_LIST, &DUM_AVG, &DUM_STD, &MEDIAN);

    STD        = STD_from_SUMS(NEVT, SUM, SQSUM) ;
    STD_ROBUST = 1.48 * MEDIAN ;
	
    OUTPUT_TABLE->MURES_STD[ibin]        = STD;
    OUTPUT_TABLE->MURES_STD_ROBUST[ibin] = STD_ROBUST ;
	
    free(ABSMURES_LIST);
  } // end ibin loop

  return ;

}  // end SUBPROCESS_COMPUTE_STD

// =====================================================
void  SUBPROCESS_OUTPUT_TABLE_RESET(int ITABLE) {

  // For each 1D bin in ITABLE, zero NEVT and MURES sums.
  int NBINTOT = SUBPROCESS.OUTPUT_TABLE[ITABLE].NBINTOT ;
  SUBPROCESS_TABLE_DEF *OUTPUT_TABLE = &SUBPROCESS.OUTPUT_TABLE[ITABLE] ;

  int  ibin1d;
  // ----------- BEGIN ---------
  for ( ibin1d=0; ibin1d < NBINTOT; ibin1d++ ) {
    OUTPUT_TABLE->NEVT[ibin1d]           = 0;
    OUTPUT_TABLE->MURES_SQSUM[ibin1d]    = 0.0 ;
    OUTPUT_TABLE->MURES_SUM[ibin1d]      = 0.0 ;
    OUTPUT_TABLE->SUM_WGT[ibin1d]        = 0.0 ;
    OUTPUT_TABLE->MURES_SUM_WGT[ibin1d]  = 0.0 ;

    OUTPUT_TABLE->MURES_STD[ibin1d]         = 0.0 ;
    OUTPUT_TABLE->MURES_STD_ROBUST[ibin1d]  = 0.0 ;
    if ( ITABLE > 0 ) {
      free(OUTPUT_TABLE->MURES_LIST[ibin1d]); // fragile?
    }
    
  }
  return;
} // end  SUBPROCESS_OUTPUT_TABLE_RESET

// ==============================================================
void SUBPROCESS_OUTPUT_TABLE_LOAD(int ISN, int ITABLE) {

  // increment table info for event index ISN and table index ITABLE.
  // 
  // Sep 30 2021: OPT_BININFO -> 0 so that events outside bin range
  //              are excluded.
  //

  SUBPROCESS_TABLE_DEF *OUTPUT_TABLE = &SUBPROCESS.OUTPUT_TABLE[ITABLE] ;

  int  NVAR         = OUTPUT_TABLE->NVAR ;
  int  NBINTOT      = OUTPUT_TABLE->NBINTOT ;
  char   *name      = INFO_DATA.TABLEVAR.name[ISN] ;
  double mures      = INFO_DATA.mures[ISN] ;
  double muerr      = INFO_DATA.muerr[ISN] ;
  double WGT        = 1.0/(muerr*muerr);
  int   ibin_per_var[MXVAR_TABLE_SUBPROCESS];
  int   NEVT, IBIN1D, IVAR, OPT_BININFO=0;
  float FVAL;        double DVAL ;
  BININFO_DEF *BININFO ;
  bool  VALID = true;
  int   LPRINT_MALLOC = 0 ;
  char fnam[] = "SUBPROCESS_OUTPUT_TABLE_LOAD" ;

  // ---------- BEGIN -----------

  for(IVAR=0; IVAR < NVAR; IVAR++ ) {
    BININFO = &OUTPUT_TABLE->BININFO[IVAR];

    // get data value for this variable and ISN event number  
    FVAL = OUTPUT_TABLE->PTRVAL[IVAR][ISN] ;
    DVAL = (double)FVAL ;   

    // convert data value to table index
    ibin_per_var[IVAR] = IBINFUN(DVAL, BININFO, OPT_BININFO, fnam );	
    
    if ( ibin_per_var[IVAR] < 0 ) { VALID = false; }
  } // end ivar      

  // bail if any event is outside bin range 

  
  if ( !VALID ) { 
    if ( ISDATA_REAL ) { 
      fprintf(FP_STDOUT,"%s: REJECT SNID=%s -> outside table bin range\n",
	      fnam, name);
    }
    return; 
  }

  // - - - - 
  // convert multiple table indices to global 1D index for table
  IBIN1D = get_1DINDEX(10+ITABLE, NVAR, ibin_per_var);


  int MEMD, LEN_REALLOC = 100; // increase to 1000 after valgrind debug 
  NEVT       = OUTPUT_TABLE->NEVT[IBIN1D];
  if ( NEVT == 0 ) {
    // malloc before any events are stored 
    MEMD = LEN_REALLOC * sizeof(double);
    OUTPUT_TABLE->MURES_LIST[IBIN1D] = (double*) malloc(MEMD); 
    if ( LPRINT_MALLOC ) {
      printf(" xxx %s: malloc MURES_LIST with MEMD=%d ITABLE=%d ibin1d=%d\n",
	     fnam, MEMD, ITABLE, IBIN1D); fflush(stdout);
    }
  }
  else if ( (NEVT % LEN_REALLOC) == 0 )  {
    // realloc
    MEMD       = (NEVT+LEN_REALLOC) * sizeof(double);
    OUTPUT_TABLE->MURES_LIST[IBIN1D] = 
      (double *)realloc(OUTPUT_TABLE->MURES_LIST[IBIN1D],MEMD);
    if ( LPRINT_MALLOC ) {
      printf(" xxx %s: realloc MURES_LIST with MEMD=%d \n",
	     fnam, MEMD); fflush(stdout);
    }
  }


  // increment table contents
  OUTPUT_TABLE->NEVT[IBIN1D]++ ;
  OUTPUT_TABLE->MURES_SUM[IBIN1D]    += mures ;
  OUTPUT_TABLE->MURES_SQSUM[IBIN1D]  += (mures*mures) ;

  //Now the weighted sums 
  SUBPROCESS.OUTPUT_TABLE[ITABLE].SUM_WGT[IBIN1D]        += WGT ;
  SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SUM_WGT[IBIN1D]  += (mures*WGT) ;

  OUTPUT_TABLE->MURES_LIST[IBIN1D][NEVT] = mures; 


  return;

} // end SUBPROCESS_OUTPUT_TABLE_LOAD


// ===========================================
void SUBPROCESS_OUTPUT_WRITE(void) {

  // write SALT2mu output:
  //   + fit params
  //   + tables

  FILE *FP_OUT = SUBPROCESS.FP_OUT ;
  int  ITER    = SUBPROCESS.ITER ;
  int  N_TABLE = SUBPROCESS.N_OUTPUT_TABLE ;

  char tmpName[40];
  int  ISFLOAT, ISM0, itable, n ;
  double VAL, ERR;
  char fnam[] = "SUBPROCESS_OUTPUT_WRITE" ;

  // ----------- BEGIN -------------

  printf("%s write SALT2mu output\n",  KEYNAME_SUBPROCESS_STDOUT );
  fflush(stdout);

  fprintf(FP_OUT,"# Created by SALT2mu SUBPROCESS\n"); // Aug 30 2021

  fprintf(FP_OUT,"# ITERATION: %d\n#\n", ITER);
  fflush(FP_OUT);

  // CPU summary  (July 29 2020)
  double t_min = (t_end_fit-t_start_fit)/60.0;
  double t_per_event = (t_end_fit-t_start_fit)/(double)FITRESULT.NSNFIT;
  fprintf(FP_OUT, "# CPU:           %.2f minutes  \n", t_min );
  fprintf(FP_OUT, "# CPU_PER_EVENT: %.1f msec/event  \n", t_per_event*1000.);
  //  fprintf(FP_OUT, "#\n");
  fflush(FP_OUT);


  fprintf(FP_OUT,"# NSNFIT: %d \n", FITRESULT.NSNFIT);
  fflush(FP_OUT);

  // always write fitted nuisance params 
  for ( n=0; n < FITINP.NFITPAR_ALL ; n++ ) {

    ISFLOAT = FITINP.ISFLOAT[n] ;
    ISM0    = (n >= MXCOSPAR) ; // it's z-binned M0

    if ( ISFLOAT && !ISM0 ) {
      VAL = FITRESULT.PARVAL[1][n] ;
      ERR = FITRESULT.PARERR[1][n] ;
      sprintf(tmpName,"%s", FITRESULT.PARNAME[n]);
      fprintf(FP_OUT, "# FITPAR:  %-14s = %10.5f +- %8.5f \n",
	      tmpName, VAL, ERR );
    }
  } // end loop over SALT2mu fit params

  // July 6 2021
  double sigint = FITINP.COVINT_PARAM_FIX ;
  sprintf(tmpName, "%s", FITRESULT.PARNAME[IPAR_COVINT_PARAM]) ;
  ERR = 0.0 ;
  fprintf(FP_OUT, "# FITPAR:  %-14s = %10.5f +- %8.5f \n",
	  tmpName, sigint, ERR );

  sprintf(tmpName, "MAXPROB_RATIO") ;
  fprintf(FP_OUT, "# FITPAR:  %-14s = %10.5f  "
	  "#Beware, value > 1 violates bounding function \n",
          tmpName, SUBPROCESS.MAXPROB_RATIO); 

  fflush(FP_OUT);

  // - - - - - - 

  for(itable=0; itable < N_TABLE; itable++ )
    { SUBPROCESS_OUTPUT_TABLE_WRITE(itable); }

  return ;

} // end SUBPROCESS_OUTPUT_WRITE

// ===================
void SUBPROCESS_OUTPUT_TABLE_WRITE(int ITABLE) {

  // Write table contents for ITABLE.
  // Write to global output file pointer SUBPROCESS.FP_OUT.

  FILE *FP_OUT       = SUBPROCESS.FP_OUT ;
  char *TABLE_NAME   = SUBPROCESS.INPUT_OUTPUT_TABLE[ITABLE];
  int   NVAR         = SUBPROCESS.OUTPUT_TABLE[ITABLE].NVAR ;
  int   NBINTOT      = SUBPROCESS.OUTPUT_TABLE[ITABLE].NBINTOT ;
  char *VARNAMES     = SUBPROCESS.OUTPUT_TABLE[ITABLE].VARNAMES_HEADER ;

  int  ivar, ibin1d, IBIN1D, NEVT, NEVT_SUM=0 ;
  double MURES_SUM, MURES_SQSUM, SUM_WGT, MURES_SUM_WGT, STD, STD_ROBUST ;
  char cLINE[200], cVAL0[100], cVAL1[100];
  char fnam[]  = "SUBPROCESS_OUTPUT_TABLE_WRITE" ;

  // ----------- BEGIN ------------

  fprintf(FP_OUT,"\n");
  fprintf(FP_OUT,"TABLE_NAME: %s\n", TABLE_NAME);
  fprintf(FP_OUT,"VARNAMES: %s\n", VARNAMES);
  fflush(FP_OUT);

  for(IBIN1D=0; IBIN1D < NBINTOT; IBIN1D++ ) {
    cLINE[0] = 0 ;

    MURES_SUM   = SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SUM[IBIN1D];
    MURES_SQSUM = SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SQSUM[IBIN1D];
    NEVT        = SUBPROCESS.OUTPUT_TABLE[ITABLE].NEVT[IBIN1D];
    NEVT_SUM   += NEVT; // diagnostic

    MURES_SUM_WGT   = SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_SUM_WGT[IBIN1D];
    SUM_WGT         = SUBPROCESS.OUTPUT_TABLE[ITABLE].SUM_WGT[IBIN1D];   

    STD = SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_STD[IBIN1D]; 
    STD_ROBUST = SUBPROCESS.OUTPUT_TABLE[ITABLE].MURES_STD_ROBUST[IBIN1D];

    for(ivar=0; ivar < NVAR; ivar++ ) {
      // get 1D bin for this variable
      ibin1d = SUBPROCESS.OUTPUT_TABLE[ITABLE].INDEX_BININFO[ivar][IBIN1D];
      sprintf(cVAL0,"%3d ", ibin1d);
      strcat(cLINE,cVAL0);
    } // end ivar


    sprintf(cVAL0," %5d  %12.4le  %12.4le %12.4le", 
	    NEVT, MURES_SUM, STD, STD_ROBUST);
    strcat(cLINE,cVAL0);

    /* xxxxxxxxxxxxx mark delete Jan 20 2022 xxx
    else {
      sprintf(cVAL0," %5d  %12.4le  %12.4le", NEVT, MURES_SUM, MURES_SQSUM);
      sprintf(cVAL1," %12.4le  %12.4le", SUM_WGT, MURES_SUM_WGT);
      strcat(cLINE,cVAL0);
      strcat(cLINE,cVAL1);
    }
    xxxxxxxxxx */

    fprintf(FP_OUT,"ROW: %4.4d %s\n", IBIN1D, cLINE);

  } // end IBIN1D
  
  fprintf(FP_OUT,"# NEVT_SUM: %d    # diagnostic\n", NEVT_SUM);
  fflush(FP_OUT);
 
  return ;

} //  end SUBPROCESS_OUTPUT_TABLE_WRITE


// =======================================
void SUBPROCESS_REMIND_STDOUT(void) {
  printf("\n");
  printf("%s STDOUT/FATAL reminder: %s\n", 
	 KEYNAME_SUBPROCESS_STDOUT, SUBPROCESS.STDOUT_FILE);
  fflush(stdout);
}

// ===============================
void SUBPROCESS_EXIT(void) {

  SUBPROCESS_REMIND_STDOUT();
  printf("%s Graceful Program Exit. Bye.\n", KEYNAME_SUBPROCESS_STDOUT);
  fflush(stdout);
  exit(0);
}

#endif
