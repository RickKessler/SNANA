! F77 -> F90 translation created 2025-11-29 with command:
!   ../util/convert_snana_f77_to_f90.py psnid.car

! Include base snana code to read/write data and apply cuts
#include "snana.F90" 

! -------------------------------------------------
! Created Oct 2012 by R.Kessler for integration into SNANA.
! PSNID = Photometric SN Identification.
!  [pulled from psnic.cra to have same structure as snana and snlc_fit]
! 
! This SNANA-interface program allows arbitrary psnid
! methods (psnid_XXX.c) analogous to the way that snlc_fit
! allows arbitrary SNIa models (genmag_XXX.c)
! 
! To add more namelist variables:
!  - add new variable(s) inside  KEEP,PSNIDINP.
!    Make sure that new variables are in a common block
!    and in the &PSNIDINP namelist block.
!  - Increment REAL*8 INPUT_ARRAY(NVAR) in PSNID_LOAD_INPUT,
!      regardless of cast
!  - Increment INPUT_STRING if needed.
!      Syntax is  "KEY1=  VALUE1  KEY2=  VALUE2 etc ..."
!  - In psnid_tools.h  define inputs in struct PSNID_INPUTS.
!  - In psnid_tools.c  update function PSNID_USER_INPUT to
!     fill new input into PSNID_INPUTS struct.
!                     OR
!   Add keys inside the namelist file, but outside the namelists.
!   In this case, you parse it yourself from PSNID_INPUTS.NMLFILE
! 
! 
! 
!   HISTORY
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Aug 28 2017: add params for MODEL1, MODEL2, MODEL3, MODEL4
!              (in addition to II,Ibc,PEC1A).
! 
! Sep 15 2017: add new &PSNIDINP inputs for
!      DMU_NON1A_MIN, DMU_NON1A_MAX,
!      TMAX_START(3), TMAX_STOP(3), TMAX_STEP(3)
! 
! Jun 19 2019:
!   + add calls to empty BEST2 functions in psnid_BEST2.c.
!     Functions will be filled in over the summer by Masao.
! 
! Nov 28 2025: rename PSNIDANA -> PSNIDVAR to avoid conflict with SUBROUTINE PSNIDANA
!              in F90 conversion
! ---------------------------------------------------

! ###############################
! ###############################



! =====================================================================
  MODULE PSNIDPAR
    USE SNPAR
    IMPLICIT NONE

! define integer parameter for each  method
! (translated from user-namelist string METHOD_NAME)

    INTEGER, PARAMETER ::       & 
         METHOD_BEST       = 1  &  ! Bayesian Evidence from SN Templates (Sako)
        ,METHOD_BEsT2      = 2  &  ! faster code, more options (Sako)
        ,MXITER_PSNID      = 3     ! Sep 2017
            

  END MODULE PSNIDPAR

! =====================================================================
  MODULE PSNIDINP_NML
    USE SNPAR
    USE PSNIDPAR
    IMPLICIT NONE

    INTEGER  & 
         NOBSMIN           &  ! I: min NOBS to process
        ,OPT_ZPRIOR        &  ! I: 0,1,2 -> FLAT, ZSPEC, Zphot(host)
        ,OPT_RATEPRIOR     &  ! I: use rate prior; default=0 --> no rate prior
        ,OPT_SIMCHEAT      &  ! I: 0 --> don't allow same PSNID template as in sim
        ,MCMC_NSTEP        &  ! I: number of MCMC steps (set to <=0 to turn off)
        ,NCOLOR, NDMU      &  ! I: number of color and delta-mu bins in grid search
        ,NREJECT_OUTLIER  ! I: max number of outliers points to reject

    CHARACTER  & 
         METHOD_NAME*60                 &  ! I: pick method name/acronym
        ,FILTLIST_FIT*(MXFILT_ALL)      &  ! I: list of filters to fit
        ,FILTLIST_PEAKMAG_STORE*(MXFILT_ALL)   &  ! list of filters to store PKMAG
        ,PRIVATE_TEMPLATES_PATH*(MXCHAR_PATH)  &  ! I: replace $SNDATA_ROOT/models/psnid
        ,TEMPLATES_SNIA*(MXCHAR_FILENAME)      &  ! I: dir or FITS file with SNIa templates
        ,TEMPLATES_NONIA*(MXCHAR_FILENAME)     &  ! I: idem for nonIa
        ,TEMPLATES_NONIA_IGNORE*(MXCHAR_FILENAME)  &  ! I: list of NONIA templates to ignore
        ,TEMPLATES_NONIA_LIST*(MXCHAR_FILENAME)    &  ! I: select this list only
        ,MODELNAME_MAGERR*60                      ! I: name of mag-error model

    REAL*8  & 
         AV_TAU, AV_SMEAR, AV_PRIOR_STR    &  ! place-holder color-prior parameter
        ,CUTWIN_ZERR(2)          &  ! I: redshift error (required for ZSPEC prior)
        ,WGT_ZPRIOR              &  ! I: weight of redshift prior (default=1)
        ,COLOR_MIN, COLOR_MAX    &  ! I: min and max color parameter in grid search
        ,DMU_MIN, DMU_MAX        &  ! I: min and max delta-mu par in grid search
        ,DMU_NON1A_MIN           &  ! I: idem, but only for NON1A
        ,DMU_NON1A_MAX           &  ! I: idem
        ,TEMPLERR_SCALE_SNIA     &  ! I: template error scale for SNIA
        ,TEMPLERR_SCALE_NONIA    &  ! I: template error scale for NONIA
        ,PBAYES_IA_CUT           &  ! I: cut on IA P_ayes
        ,PBAYES_IBC_CUT          &  ! I: cut on IBC P_Bayes
        ,PBAYES_II_CUT           &  ! I: idem for SII
        ,PBAYES_PEC1A_CUT        &  ! I: idem for Peculiar SNIa
        ,PBAYES_MODEL1_CUT       &  ! I: idem for MODEL1 (e.g., KN)
        ,PBAYES_MODEL2_CUT       &  ! I:
        ,PBAYES_MODEL3_CUT       &  ! I:
        ,PBAYES_MODEL4_CUT       &  ! I:
        ,FITPROB_IA_CUT          &  ! I: min fitprob to type a Ia
        ,FITPROB_IBC_CUT         &  ! I: idem for IBC
        ,FITPROB_II_CUT          &  ! I: idem for SNII
        ,FITPROB_PEC1A_CUT       &  ! I: idem for Pec. SNIa
        ,FITPROB_MODEL1_CUT      &  ! I: idem for MODEL1
        ,FITPROB_MODEL2_CUT      &  ! I:
        ,FITPROB_MODEL3_CUT      &  ! I:
        ,FITPROB_MODEL4_CUT      &  ! I: idem for MODEL4
        ,ZRATEPRIOR_SNIA(3)      &  ! I: a1*(1+z)**a2, then flat for z>a3
        ,ZRATEPRIOR_NONIA(3)     &  ! I: a1*(1+z)**a2, then flat for z>a3
        ,CHISQMIN_OUTLIER        &  ! I: min chi2 for outlier rejection
        ,MJDFIT_RANGE(2)         &  ! I: restrict fit to this MJD range
        ,TMAX_START(MXITER_PSNID)    &  ! I: TMAX start for each iteration
        ,TMAX_STOP(MXITER_PSNID)     &  ! I: TMAX end for each iteration
        ,TMAX_STEP(MXITER_PSNID)    ! I: TMAX step for each iteration



! define namelist to read from input file.

    NAMELIST / PSNIDINP /  & 
         METHOD_NAME, NOBSMIN  & 
        ,FILTLIST_FIT, FILTLIST_PEAKMAG_STORE  & 
        ,PRIVATE_TEMPLATES_PATH, TEMPLATES_SNIA, TEMPLATES_NONIA  & 
        ,TEMPLATES_NONIA_IGNORE, TEMPLATES_NONIA_LIST  & 
        ,AV_TAU, AV_SMEAR, AV_PRIOR_STR, WGT_ZPRIOR, CUTWIN_ZERR  & 
        ,OPT_ZPRIOR, OPT_RATEPRIOR, OPT_SIMCHEAT  & 
        ,MCMC_NSTEP, NCOLOR, NDMU  & 
        ,COLOR_MIN, COLOR_MAX, DMU_MIN, DMU_MAX  & 
        ,DMU_NON1A_MIN, DMU_NON1A_MAX  & 
        ,TEMPLERR_SCALE_SNIA, TEMPLERR_SCALE_NONIA  & 
        ,PBAYES_IA_CUT, PBAYES_IBC_CUT, PBAYES_II_CUT, PBAYES_PEC1A_CUT  & 
        ,PBAYES_MODEL1_CUT, PBAYES_MODEL2_CUT  & 
        ,PBAYES_MODEL3_CUT, PBAYES_MODEL4_CUT  & 
        ,FITPROB_IA_CUT, FITPROB_IBC_CUT, FITPROB_II_CUT  & 
        ,FITPROB_PEC1A_CUT  & 
        ,FITPROB_MODEL1_CUT, FITPROB_MODEL2_CUT  & 
        ,FITPROB_MODEL3_CUT, FITPROB_MODEL4_CUT  & 
        ,ZRATEPRIOR_SNIA, ZRATEPRIOR_NONIA  & 
        ,MODELNAME_MAGERR  & 
        ,CHISQMIN_OUTLIER, NREJECT_OUTLIER, MJDFIT_RANGE  & 
        ,TMAX_START, TMAX_STOP, TMAX_STEP

  END MODULE PSNIDINP_NML

! =====================================================================
  MODULE PSNIDVAR
    USE SNPAR
    USE PSNIDPAR
    IMPLICIT NONE

! variables computed from PSNIDINP
    INTEGER  & 
         METHOD_PSNID  & 
        ,NFILT_PSNID              &  ! size of FILTLIST_FIT
        ,IFILTMAP_PSNID(100)      &  ! array of IFILT_OBS
        ,IFILTINV_PSNID(100)     ! inverse of above mmap

    LOGICAL  & 
         USEFILT_PSNID(100) ! TRUE if this IFILT_OBS is fit


! --------

    CHARACTER CCID_forC*(MXCHAR_CCID)


  END MODULE PSNIDVAR

! =====================================================================
  MODULE PSNIDMON
    USE SNPAR
    USE PSNIDPAR
    IMPLICIT NONE

    INTEGER, PARAMETER :: MXTYPE_PSNID = 400

    INTEGER  & 
         ITYPE_IA            &  ! index of IA
        ,NTYPE_PSNID(0:MXTYPE_PSNID)   &  ! number of TYPE per TYPE
        ,NSIM_GENIA          &  ! Number of generate IA
        ,NSIM_PSNIDIA       ! Number of correctly typed IA in sim



  END MODULE PSNIDMON

! =====================================================================
  MODULE PSNIDCOM
    USE SNPAR

    USE PSNIDPAR
    USE PSNIDINP_NML
    USE PSNIDVAR
    USE PSNIDMON


  END MODULE PSNIDCOM


! ###############################
! ###############################



! ===========================================
! =====================================================================
! =====================================================================
! =====================================================================
    SUBROUTINE PSNIDINI(IERR)
! 
! User init-routine called before SN are read in.
! 
! Nov 1 2013 RK - fix init sequence by calling PSNID_LOAD_INPUT()
!                 right after PSNID_INIT_VAR().
! 
! Aug 27 2017 RK - call PSNID_READ_TEMPLATE only if string len > 0
!                  (e.g., allows using only CC or SNIa templates)
! 
! -------------------------------------

    USE SNPAR
    USE PSNIDCOM
    USE SNCUTS
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IERR ! (O) 0 => OK

! local var

    INTEGER  & 
         LEN_FILT, TYPEINDX, i  & 
        ,L1, L2, LEN_PATH

    CHARACTER  & 
         FILTLIST_forC*(MXFILT_ALL)  & 
        ,FILE1_forC*(MXCHAR_FILENAME)  & 
        ,FILE2_forC*(MXCHAR_FILENAME)  & 
        ,STRING1_forC*(MXCHAR_FILENAME)  & 
        ,STRING2_forC*(MXCHAR_FILENAME)  & 
        ,C1ERR*80, C2ERR*80


! functions
    EXTERNAL  & 
         PSNID_INIT_VAR  & 
        ,PSNID_INIT_FILTERS  & 
        ,PSNID_READ_TEMPLATES  & 
        ,PSNID_INIT_TEMPLATES  & 
        ,PSNID_INIT_XTMW  & 
        ,PSNID_BEST_INIT

! ------------------- BEGIN ---------------------

    IERR = 0  ! set status to OK

    print*,'  PSNIDINI: Call PSNID-init Functions  '

! read user namelist
    CALL RDPSNIDNML()

    call ENVreplace(TEMPLATES_SNIA)
    call ENVreplace(TEMPLATES_NONIA)
    call ENVreplace(PRIVATE_TEMPLATES_PATH )

! -----------------------------------------------
!  prepare strings to pass as arguments to C functions

    LEN_FILT      = INDEX(FILTLIST_FIT,' ') - 1
    FILTLIST_forC = FILTLIST_FIT(1:LEN_FILT) // char(0)

! ------------------------------------------------
! start with generic inits for any method.

    CALL PSNID_INIT_VAR()  ! init variables

! prepare  list of namelist inputs and pass to C code
    CALL PSNID_LOAD_INPUT()

! init filter list
    CALL PSNID_INIT_FILTERS(FILTLIST_forC,LEN_FILT)

! init SN templates
    L1 = INDEX(TEMPLATES_SNIA, ' ') - 1
    L2 = INDEX(TEMPLATES_NONIA,' ') - 1
    FILE1_forC = TEMPLATES_SNIA(1:L1)  // char(0)
    FILE2_forC = TEMPLATES_NONIA(1:L2) // char(0)

! check for PRIVATE path option

    LEN_PATH = INDEX(PRIVATE_TEMPLATES_PATH,' ') - 1
    if ( LEN_PATH .GT. 1 ) then
      FILE1_forC = PRIVATE_TEMPLATES_PATH(1:LEN_PATH)  & 
                       // '/' // FILE1_forC
      FILE2_forC = PRIVATE_TEMPLATES_PATH(1:LEN_PATH)  & 
                       // '/' // FILE2_forC
    endif


    IF ( L1 > 0 ) THEN
       TYPEINDX = 1  ! SNIa
       CALL PSNID_READ_TEMPLATES(TYPEINDX, FILE1_forC, L1)
    ENDIF

    IF ( L2 > 0 ) THEN
       TYPEINDX = 2  ! NONIA
      CALL PSNID_READ_TEMPLATES(TYPEINDX, FILE2_forC, L2)
    ENDIF

! check and init templates after they have all been read.
! Pass list of NONIA templates to ignore
    L1 = INDEX(TEMPLATES_NONIA_LIST, ' ', BACK = .TRUE.)
    STRING1_forC = TEMPLATES_NONIA_LIST(1:L1) // char(0)

    L2 = INDEX(TEMPLATES_NONIA_IGNORE, ' ', BACK = .TRUE.)
    STRING2_forC = TEMPLATES_NONIA_IGNORE(1:L2) // char(0)

    CALL PSNID_INIT_TEMPLATES(STRING1_forC, STRING2_forC, L1, L2)

! ---------------------------------------------------
! init MilkyWay Galactic extinction (after reading templates with LAMAVG)
    CALL PSNID_INIT_XTMW(OPT_MWCOLORLAW)

! -----------------------------------------------
! method-specific inits.

    METHOD_PSNID = -9

    IF ( METHOD_NAME .EQ. 'BEST' ) THEN

! set integer flag for this method
      METHOD_PSNID = METHOD_BEST
      CALL PSNID_BEST_INIT()

! ------------
    ELSE IF ( METHOD_NAME .EQ. 'BEST2' ) THEN  ! begin June 19 2019

      METHOD_PSNID = METHOD_BEST2
      CALL PSNID_BEST2_INIT()

!      ELSE IF ( METHOD_NAME .EQ. 'SOFT' ) THEN
! ccccccccccc PLACE-HOLDER
! ------------

    ELSE
      C1ERR = 'Unknown METHOD_NAME = ' //  METHOD_NAME(1:40)
      C2ERR = 'Check  &PSNIDINP  namelist.'
      CALL MADABORT('PSNIDINI', C1ERR,  C2ERR )
    ENDIF

    print*,'   '
    CALL FLUSH(6)


! init monitor array
    DO i = 0, MXTYPE_PSNID
       NTYPE_PSNID(i) = 0
    ENDDO
    NSIM_GENIA   = 0
    NSIM_PSNIDIA = 0

    RETURN
  END SUBROUTINE PSNIDINI

! =====================================
    SUBROUTINE PSNIDINI2(IERR)
! 
! Dec 4, 2012:
! Second init, called after first SN is read so that
! we have more info such as the SIM_XXX parameter values.
! 


    USE SNPAR
    USE CTRLCOM
    USE PSNIDCOM

    IMPLICIT NONE

    INTEGER IERR ! (O) 0=> ok

! local args

    CHARACTER  & 
        C1ERR*80, C2ERR*80, FMT*20, FMT_forC*20, FNAM*10

    INTEGER LEN_FMT, OPT

! variables for legacy options
    CHARACTER DMPFILE_forC*(MXCHAR_FILENAME)
    INTEGER LEN_DMPFILE

! external functions
    EXTERNAL  PSNID_BEST_INIT_SNTABLE

! ----------- BEGIN ---------
    IERR = 0
    FNAM = 'PSNIDINI2'

    OPT      = OPT_TABLE(ITABLE_FITRES)
    FMT      = TEXTFORMAT_TABLE(ITABLE_FITRES)
    LEN_FMT  = INDEX(FMT,' ') - 1
    FMT_forC = FMT(1:LEN_FMT) // char(0)


! method-specific inits.

    IF ( METHOD_NAME .EQ. 'BEST' ) THEN

      CALL PSNID_BEST_INIT_SNTABLE(OPT,FMT_forC,LSIM_SNANA,LEN_FMT)

! ------------
    ELSE IF ( METHOD_NAME .EQ. 'BEST2' ) THEN

      CALL PSNID_BEST2_INIT_SNTABLE(OPT,FMT_forC,LSIM_SNANA,LEN_FMT)



!      ELSE IF ( METHOD_NAME .EQ. 'SOFT' ) THEN
! ccccccccccc PLACE-HOLDER
! ------------

    ELSE
      C1ERR = 'Unknown METHOD_NAME = ' //  METHOD_NAME(1:40)
      C2ERR = 'Check  &PSNIDINP  namelist.'
      CALL MADABORT(FNAM, C1ERR,  C2ERR )
    ENDIF

    RETURN
  END SUBROUTINE PSNIDINI2
! ===========================================
    SUBROUTINE PSNIDANA(IERR)
! 
! User analysis routine called after each Supernova has been
! read and analyzed by SNANA.
! See SNLC_XXX varibales in $SNANA_DIR/src/snana.car
! 
! Apr  2, 2013: add XTMW8 argument to PSNID_BEST_DOFIT
! Aug 23, 2015: add MJD cut to select epochs for fit
! Mar 10, 2017: pass SIM_NON1a to PSNID_BEST_DOFIT
! Jan 08, 2020: pass FLUXSIM
! -------------------------------------

    USE SNDATCOM
    USE PSNIDCOM
    USE FILTCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! declare subroutine args

    INTEGER ISN  ! (I) sparse SN index
    INTEGER IERR ! (O) 0 => OK

! local variables

    INTEGER  & 
         NOBS, EPMIN, EPMAX, ep, IMJD  & 
        ,IFILTLIST_PSNID(MXEPOCH)       &  ! sparse index, 1-NFILT
        ,IFILT, IFILT_OBS, LENCCID

    REAL*8  & 
         MJD8(MXEPOCH)  & 
        ,FLUXDATA8(MXEPOCH)  & 
        ,FLUXERR8(MXEPOCH)  & 
        ,FLUXSIM8(MXEPOCH)  & 
        ,XTMW8(MXEPOCH)  & 
        ,Z8(3), ZERR8(3)  & 
        ,MWEBV8, MWEBVERR8


    CHARACTER CCID*(MXCHAR_CCID), CFILT*2

! functions
    INTEGER   PSNID_BEST_DOFIT, PSNID_BEST2_DOFIT

    EXTERNAL  & 
         PSNID_BEST_DOFIT, PSNID_BEST2_DOFIT  & 
        ,PSNID_BEST_UPDATE_OUTPUT, PSNID_BEST2_UPDATE_OUTPUT

! ------------------- BEGIN ---------------------

    IERR = 0  ! return status is OK

    LENCCID    = ISNLC_LENCCID
    CCID       = SNLC_CCID(1:LENCCID)
    CCID_forC  = SNLC_CCID(1:LENCCID) // char(0)

    Z8(1)    = SNLC_REDSHIFT
    ZERR8(1) = SNLC_REDSHIFT_ERR
    Z8(2)    = SNHOST_ZPHOT(1)
    ZERR8(2) = SNHOST_ZPHOT_ERR(1)
    Z8(3)    = -9.0000
    ZERR8(3) = -9.0000

    IF ( USE_MWCOR ) THEN  ! fluxes already corrected for MWEBV
      MWEBV8    = 0.0
      MWEBVERR8 = 0.0
    ELSE                   ! apply MWEBV to fit-model
      MWEBV8    = SNLC_MWEBV
      MWEBVERR8 = SNLC_MWEBV_ERR
    ENDIF

! strip light curve into local double-precision array

    NOBS = 0
    DO 101 IMJD  = 1, ISNLC_NEWMJD_STORE
       EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,IMJD)
       EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,IMJD)
    DO 102 ep = EPMIN, EPMAX


! MJD cut (Aug 23 2015)
       IF ( SNLC8_MJD(ep) < MJDFIT_RANGE(1) ) GOTO 102
       IF ( SNLC8_MJD(ep) > MJDFIT_RANGE(2) ) GOTO 102

       IFILT_OBS   = ISNLC_IFILT_OBS(ep)
       IFILT       = IFILTDEF_INVMAP_SURVEY(IFILT_OBS)

       NOBS            = NOBS + 1
       MJD8(NOBS)      = SNLC8_MJD(ep)
       FLUXDATA8(NOBS) = DBLE( SNLC_FLUXCAL(ep) )
       FLUXERR8(NOBS)  = DBLE ( SNLC_FLUXCAL_ERRTOT(ep) )
       FLUXSIM8(NOBS)  = DBLE ( SIM_EPFLUXCAL(ep) )  ! Jan 2020
!          XTMW8(NOBS)     = SNLC_XTMW_FLUXFRAC(ifilt)

       IFILTLIST_PSNID(NOBS) = IFILTINV_PSNID(IFILT_OBS)

102   CONTINUE
101   CONTINUE

    IF ( NOBS .LT. NOBSMIN ) RETURN

        N_SNLC_FIT = N_SNLC_FIT + 1

    IF ( METHOD_PSNID .EQ. METHOD_BEST ) THEN

       ERRFLAG_FIT = PSNID_BEST_DOFIT(CCID_forC, NOBS,  & 
                           IFILTLIST_PSNID,  & 
                           MJD8, FLUXDATA8, FLUXERR8,  & 
                           FLUXSIM8,       &  ! Jan 2020
                           Z8, ZERR8,  & 
                           MWEBV8, MWEBVERR8, SIM_TEMPLATE_INDEX,  & 
                           LENCCID ) ;

    ELSE IF ( METHOD_PSNID .EQ. METHOD_BEST2 ) THEN

       ERRFLAG_FIT = PSNID_BEST2_DOFIT(CCID_forC, NOBS,  & 
                           IFILTLIST_PSNID,  & 
                           MJD8, FLUXDATA8, FLUXERR8,  & 
                           Z8, ZERR8,  & 
                           MWEBV8, MWEBVERR8, SIM_TEMPLATE_INDEX,  & 
                           LENCCID ) ;

    ELSE
      write(C1ERR,666) METHOD_PSNID, METHOD_NAME(1:20)
666     format('Invalid METHOD_PSNID = ', I5,'  for METHOD_NAME=',A)
      C2ERR = 'METHOD_NAME satisfied PSNIDINI,  but fails here ??'
      CALL MADABORT('PSNIDANA', C1ERR,  C2ERR )
    ENDIF

! increment number with useful fit for stat summary at end
    IF ( ERRFLAG_FIT .EQ. 0 ) THEN
       N_SNLC_FITCUTS = N_SNLC_FITCUTS + 1

       IF ( METHOD_PSNID .EQ. METHOD_BEST ) THEN
          CALL PSNID_BEST_UPDATE_OUTPUT(CCID_forC,LENCCID)

       ELSE IF ( METHOD_PSNID .EQ. METHOD_BEST2 ) THEN
          CALL PSNID_BEST2_UPDATE_OUTPUT(CCID_forC,LENCCID)

       ELSE

       ENDIF

       CALL PSNID_SNLCPLOT()

    ENDIF


! monitor Number of SN per fitted TYPE
    CALL PSNID_MONANA();

    RETURN
  END SUBROUTINE PSNIDANA

! ===========================================
    SUBROUTINE PSNIDEND(IERR)
! 
! User end-routine after all analysis/fits are done.
! [close files, summarize statistics, global analysis, etc ...]
! -------------------------------------

    USE SNDATCOM

    IMPLICIT NONE

! xxxxxxxx mark del Nov 28 2025 xxxxxx
! xx+CDE,SNANAFIT.
! xx+CDE,SNFITCOM.
! xxxxxxxxxxxxxxxxx

    INTEGER IERR ! (O) 0 => OK
! local var

    EXTERNAL PSNID_BEST_SUMMARY

! ------------------- BEGIN ---------------------
    IERR = 0

    CALL PSNID_BEST_SUMMARY()

! ---
    CALL PSNID_MONEND()

    RETURN
  END SUBROUTINE PSNIDEND


! =============================
    SUBROUTINE RDPSNIDNML()

! Read PSNIDINP namelist


    USE SNDATCOM
    USE SNLCINP_NML
    USE PSNIDCOM

    IMPLICIT NONE


    INTEGER   IFILT, IFILT_OBS, i, LEN_ROOT, LEN_PATH
    CHARACTER CFILT*2, TMPDIR*200

! functions
    INTEGER FILTINDX

! ------------ BEGIN -----------

! init namelist parameters before reading
! (don't be lazy ... INIT THEM ALL !!!)

    METHOD_NAME = 'UNKNOWN'

    NOBSMIN       = 1
    FILTLIST_FIT            = ' '
    FILTLIST_PEAKMAG_STORE  = ' '

    TEMPLATES_SNIA  = ''
    TEMPLATES_NONIA = ''
    PRIVATE_TEMPLATES_PATH = ''
    TEMPLATES_NONIA_IGNORE = ''
    TEMPLATES_NONIA_LIST   = ''

    MODELNAME_MAGERR = 'S13'  ! default mag-error model

! -----------
    AV_TAU       = -9.9
    AV_SMEAR     = -9.9    ! AV_TAU,AV_SMEAR < 0 -> flat AV prior
    AV_PRIOR_STR =  1.0

    WGT_ZPRIOR      =  1.0    ! external z-prior weight

    COLOR_MIN    =  -1.0
    COLOR_MAX    =   1.0
    NCOLOR       =  21

    DMU_MIN      =   0.0
    DMU_MAX      =   0.0
    NDMU         =   1

    DMU_NON1A_MIN = 0.0
    DMU_NON1A_MAX = 0.0

    OPT_ZPRIOR    = 0  ! default is flat prior
    CUTWIN_ZERR(1) = -9.9
    CUTWIN_ZERR(2) = +9.9

    OPT_RATEPRIOR = 0             ! default is no rate prior, as in S11
    ZRATEPRIOR_SNIA(1)  = 2.6E-5  ! SNLS rate result
    ZRATEPRIOR_SNIA(2)  = 2.2     ! (1+z)**2.2
    ZRATEPRIOR_SNIA(3)  = 1.0     ! flat rate for z>1
    ZRATEPRIOR_NONIA(1) = 6.8E-5  ! Bazin et al 2009 (SNLS)
    ZRATEPRIOR_NONIA(2) = 3.6     ! follows Galaxy formation
    ZRATEPRIOR_NONIA(3) = 2.0     ! flat rate for z>2

! - - - - - -  -
! Mar 10 2017 RK -
!   default is to cheat and use PSNID template if it is the
!   same as the simulated template for each event.
!   Set to 0 to turn off cheating.
    OPT_SIMCHEAT=1
! - - - - - -  -

    MCMC_NSTEP = 0  ! default is no MCMC

    TEMPLERR_SCALE_SNIA  = 1.0
    TEMPLERR_SCALE_NONIA = 1.0

    PBAYES_IA_CUT     = 0.0 ! default min
    PBAYES_IBC_CUT    = 0.6 ! default min (0.9 -> 0.6 Nov 2025)
    PBAYES_II_CUT     = 0.6 ! default min (0.9 -> 0.6 Nov 2025)
    PBAYES_PEC1A_CUT  = 0.9
    PBAYES_MODEL1_CUT = 0.0
    PBAYES_MODEL2_CUT = 0.0
    PBAYES_MODEL3_CUT = 0.0
    PBAYES_MODEL4_CUT = 0.0

    FITPROB_IA_CUT     = 0.01
    FITPROB_IBC_CUT    = 0.01
    FITPROB_II_CUT     = 0.01
    FITPROB_PEC1A_CUT  = 0.01
    FITPROB_MODEL1_CUT = 0.0
    FITPROB_MODEL2_CUT = 0.0
    FITPROB_MODEL3_CUT = 0.0
    FITPROB_MODEL4_CUT = 0.0

    CHISQMIN_OUTLIER = 1.0E9   ! default to large min chi2
    NREJECT_OUTLIER  = 0       ! default is to reject nothing

    MJDFIT_RANGE(1)  = 0.
    MJDFIT_RANGE(2)  = 9999999.

! TMAX grid search value taken from psnid_best_setup_searchgrid()
    TMAX_START(1) = -40.0
    TMAX_START(2) = -15.0
    TMAX_START(3) =  -5.0

    TMAX_STOP(1) =  60.0
    TMAX_STOP(2) =  15.0
    TMAX_STOP(3) =   5.0

    TMAX_STEP(1) =  10.0
    TMAX_STEP(2) =   3.0
    TMAX_STEP(3) =   1.0

! -------------------------------------------------
! read the namelist file for the PSNIDINP namelist

    OPEN (  & 
          UNIT   = LUNNML  & 
         ,file   = nmlfile   &  ! nml file with  SNLCINP
         ,status = 'OLD'  & 
         ,ERR    = 900  & 
             )

    READ(LUNNML, NML = PSNIDINP, ERR=901, END=902 )
    CLOSE ( UNIT = LUNNML )

    print*,' '

! check for command-line overrides
    CALL PSNIDINP_OVERRIDE

! dump namelist params to std out

    IF ( .not. REDUCE_STDOUT_BATCH ) THEN
      WRITE ( 6 , NML = PSNIDINP )
      CALL FLUSH(6)
    ENDIF

! -----------------------------------------------------
! compute number of filters to fit and load IFILTMAP_PSNID
    DO IFILT = 1, MXFILT_OBS
      IFILTMAP_PSNID(IFILT) = -9
      IFILTINV_PSNID(IFILT) = -9
      USEFILT_PSNID(IFILT)  = .FALSE.
    ENDDO
    NFILT_PSNID = INDEX(FILTLIST_FIT,' ' ) - 1
    IF ( NFILT_PSNID .LE. 0 .OR.  & 
           NFILT_PSNID .GT. MXFILT_OBS ) then
      write(c1err,661) NFILT_PSNID
661     format('Invalid NFILT_PSNID = ', I3  )
      c2err = 'Check &PSNIDINP namelist variable FILTLIST_FIT'
      CALL MADABORT('RDPSNIDNML', C1ERR,  nmlfile)
    ENDIF

    DO IFILT = 1, NFILT_PSNID
      cfilt     = FILTLIST_FIT(IFILT:IFILT)
      IFILT_OBS = FILTINDX(cfilt)
      IFILTMAP_PSNID(IFILT)     = IFILT_OBS
      IFILTINV_PSNID(IFILT_OBS) = IFILT - 1
      USEFILT_PSNID(IFILT_OBS)  = .TRUE.
    ENDDO

! ----------------------------------------------
! Feb 2015: allow $SNDATA_ROOT to be used in PRIVATE_TEMPLATES_PATH

    if ( PRIVATE_TEMPLATES_PATH(1:12) .EQ. '$SNDATA_ROOT' ) then
       LEN_ROOT = INDEX(SNDATA_ROOT,' ') - 1
       LEN_PATH = INDEX(PRIVATE_TEMPLATES_PATH,' ') - 1
       PRIVATE_TEMPLATES_PATH =  & 
             SNDATA_ROOT(1:LEN_ROOT) //  & 
             PRIVATE_TEMPLATES_PATH(13:LEN_PATH)
    endif
! -----------------------------------------------------

    RETURN

! ---------------------------------------------------
900   C1ERR = 'Could not open namelist file '
    CALL MADABORT('RDPSNIDNML', C1ERR,  nmlfile)

901   C1ERR = 'could not read &PSNIDINP namelist from '
    CALL MADABORT('RDPSNIDNML', C1ERR,  nmlfile)

902   C1ERR = 'Could not find &PSNIDINP namelist in'
    CALL MADABORT('RDPSNIDNML', C1ERR,  nmlfile)

    RETURN
  END SUBROUTINE RDPSNIDNML

! =========================
    SUBROUTINE PSNIDINP_OVERRIDE()

! Check for command-line overrides of namelist variables
! Each new nml variable must be included below.
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE PSNIDPAR
    USE PSNIDINP_NML

    IMPLICIT NONE

! local var
    INTEGER  i, LL, ilast, iuse, iter
    CHARACTER ARG*(MXCHAR_ARG), ARGLIST(MXKEY_ARGS)*(MXCHAR_ARG)
    LOGICAL   MATCH_NMLKEY

! ----------------- BEGIN -------------

    i = 2
    ilast = 2

    DO WHILE ( i .LT. NLINE_ARGS )

       CALL GETARG(i,ARG)
       LL = INDEX ( ARG, ' ' ) - 1

     if ( MATCH_NMLKEY('METHOD_NAME', 1,i,ARGLIST) ) then
          METHOD_NAME = ARGLIST(1)

      else if ( MATCH_NMLKEY('FILTLIST_FIT', 1,i,ARGLIST) ) then
          FILTLIST_FIT = ARGLIST(1)

      else if ( MATCH_NMLKEY('FILTLIST_PEAKMAG_STORE',  & 
                   1, i, ARGLIST) ) then
          FILTLIST_PEAKMAG_STORE = ARGLIST(1)

      else if ( MATCH_NMLKEY('TEMPLATES_SNIa TEMPLATES_SNIA',  & 
                   1, i, ARGLIST) ) then
          TEMPLATES_SNIa = ARGLIST(1)

      else if ( MATCH_NMLKEY('TEMPLATES_NONIa TEMPLATES_NONIA',  & 
                   1, i, ARGLIST) ) then
          TEMPLATES_NONIa = ARGLIST(1)

      else if(MATCH_NMLKEY('PRIVATE_TEMPLATES_PATH',1,i,ARGLIST))then
          PRIVATE_TEMPLATES_PATH = ARGLIST(1)

      else if(MATCH_NMLKEY('MODELNAME_MAGERR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) MODELNAME_MAGERR

      else if(MATCH_NMLKEY('TEMPLATES_NONIA_LIST', 1,i,ARGLIST))then
          TEMPLATES_NONIA_LIST = ARGLIST(1)

      else if(MATCH_NMLKEY('TEMPLATES_NONIA_IGNORE', 1,i,ARGLIST))then
          TEMPLATES_NONIA_IGNORE = ARGLIST(1)

      else if(MATCH_NMLKEY('AV_TAU', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) AV_TAU

      else if(MATCH_NMLKEY('AV_SMEAR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) AV_SMEAR

      else if(MATCH_NMLKEY('AV_SMEAR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) AV_SMEAR

      else if(MATCH_NMLKEY('AV_PRIOR_STR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) AV_PRIOR_STR

      else if(MATCH_NMLKEY('WGT_ZPRIOR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) WGT_ZPRIOR

      else if(MATCH_NMLKEY('OPT_ZPRIOR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) OPT_ZPRIOR

      else if(MATCH_NMLKEY('CUTWIN_ZERR', 2,i,ARGLIST))then
          READ(ARGLIST(1),*) CUTWIN_ZERR(1)
          READ(ARGLIST(2),*) CUTWIN_ZERR(2)

      else if(MATCH_NMLKEY('OPT_RATEPRIOR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) OPT_RATEPRIOR

      else if(MATCH_NMLKEY('ZRATEPRIOR_SNIA', 3,i,ARGLIST))then
          READ(ARGLIST(1),*) ZRATEPRIOR_SNIA(1)
          READ(ARGLIST(2),*) ZRATEPRIOR_SNIA(2)
          READ(ARGLIST(3),*) ZRATEPRIOR_SNIA(3)

      else if(MATCH_NMLKEY('ZRATEPRIOR_NONIA', 3,i,ARGLIST))then
          READ(ARGLIST(1),*) ZRATEPRIOR_NONIA(1)
          READ(ARGLIST(2),*) ZRATEPRIOR_NONIA(2)
          READ(ARGLIST(3),*) ZRATEPRIOR_NONIA(3)

      else if(MATCH_NMLKEY('OPT_SIMCHEAT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) OPT_SIMCHEAT

      else if(MATCH_NMLKEY('MCMC_NSTEP', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) MCMC_NSTEP

      else if(MATCH_NMLKEY('COLOR_MIN', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) COLOR_MIN

      else if(MATCH_NMLKEY('COLOR_MAX', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) COLOR_MAX

      else if(MATCH_NMLKEY('NCOLOR', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) NCOLOR

      else if(MATCH_NMLKEY('DMU_MIN', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) DMU_MIN

      else if(MATCH_NMLKEY('DMU_MAX', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) DMU_MAX

      else if(MATCH_NMLKEY('DMU_NON1A_MIN', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) DMU_NON1A_MIN

      else if(MATCH_NMLKEY('DMU_NON1A_MAX', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) DMU_NON1A_MAX

      else if(MATCH_NMLKEY('NDMU', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) NDMU

      else if(MATCH_NMLKEY('TEMPLERR_SCALE_SNIA TEMPLERR_SCALE_SNIa',  & 
                 1, i, ARGLIST) )then
          READ(ARGLIST(1),*) TEMPLERR_SCALE_SNIA

      else if(MATCH_NMLKEY(  & 
          'TEMPLERR_SCALE_NONIA TEMPLERR_SCALE_NONIa',1,i,ARGLIST))then
           READ(ARGLIST(1),*) TEMPLERR_SCALE_NONIA

      else if(MATCH_NMLKEY('TMAX_START',MXITER_PSNID,i,ARGLIST))then
          do iter = 1, MXITER_PSNID
             READ(ARGLIST(iter),*) TMAX_START(iter)
          enddo

      else if(MATCH_NMLKEY('TMAX_STOP', MXITER_PSNID,i,ARGLIST) )then
          do iter = 1, MXITER_PSNID
             READ(ARGLIST(iter),*) TMAX_STOP(iter)
          enddo

      else if(MATCH_NMLKEY('TMAX_STEP', MXITER_PSNID,i,ARGLIST) )then
          do iter = 1, MXITER_PSNID
             READ(ARGLIST(iter),*) TMAX_STEP(iter)
          enddo

! - - - - - -

      else if(MATCH_NMLKEY('PBAYES_IA_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_IA_CUT

      else if(MATCH_NMLKEY('PBAYES_IBC_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_IBC_CUT

      else if(MATCH_NMLKEY('PBAYES_II_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_II_CUT

      else if(MATCH_NMLKEY('PBAYES_PEC1A_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_PEC1A_CUT

      else if(MATCH_NMLKEY('PBAYES_MODEL1_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_MODEL1_CUT
      else if(MATCH_NMLKEY('PBAYES_MODEL2_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_MODEL2_CUT
      else if(MATCH_NMLKEY('PBAYES_MODEL3_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_MODEL3_CUT
      else if(MATCH_NMLKEY('PBAYES_MODEL4_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) PBAYES_MODEL4_CUT

      else if(MATCH_NMLKEY('FITPROB_IA_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_IA_CUT
      else if(MATCH_NMLKEY('FITPROB_IBC_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_IBC_CUT
      else if(MATCH_NMLKEY('FITPROB_II_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_II_CUT

      else if(MATCH_NMLKEY('FITPROB_PEC1A_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_PEC1A_CUT

      else if(MATCH_NMLKEY('FITPROB_MODEL1_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_MODEL1_CUT
      else if(MATCH_NMLKEY('FITPROB_MODEL2_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_MODEL2_CUT
      else if(MATCH_NMLKEY('FITPROB_MODEL3_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_MODEL3_CUT
      else if(MATCH_NMLKEY('FITPROB_MODEL4_CUT', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) FITPROB_MODEL4_CUT

      else if(MATCH_NMLKEY('CHISQMIN_OUTLIER', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) CHISQMIN_OUTLIER

      else if(MATCH_NMLKEY('MJDFIT_RANGE', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) MJDFIT_RANGE(1)
          READ(ARGLIST(2),*) MJDFIT_RANGE(2)

      else if(MATCH_NMLKEY('NREJECT_OUTLIER', 1,i,ARGLIST))then
          READ(ARGLIST(1),*) NREJECT_OUTLIER

! xxx add more here ....

       endif

! set logical flag for used LINE ARGS

       IF ( i .GT. ilast ) THEN
         DO iuse = ilast, i
            USE_LINE_ARGS(iuse) = .TRUE.
         ENDDO
       ENDIF

       i = i + 1
       ilast = i

    ENDDO

    RETURN
  END SUBROUTINE PSNIDINP_OVERRIDE

! ==================================
    SUBROUTINE PSNID_LOAD_INPUT
! 
! Transfer namelist arguments to linear arrays,
! then pass them to C function.
! 
! Feb 25 2020: pass DEBUG_FLAG to C code.


    USE SNDATCOM
    USE SNLCINP_NML
    USE PSNIDCOM

    IMPLICIT NONE


    INTEGER NVAR, L1, LSTR, LNN, NFILE, i, IVAR, iter
    REAL*8    INPUT_ARRAY(200)
    CHARACTER INPUT_STRING*2000

    EXTERNAL PSNID_USER_INPUT

! ---------------- BEGIN -----------

    NVAR = 0

! start with &SNLCINP values
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = H0_REF(1)
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = OMAT_REF(1)
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = OLAM_REF(1)
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = W0_REF(1)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(DEBUG_FLAG)

! next the  INPUT ARRAY from &PSNIDINP

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = AV_TAU

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = AV_SMEAR

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = AV_PRIOR_STR

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = WGT_ZPRIOR

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(OPT_ZPRIOR)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = CUTWIN_ZERR(1)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = CUTWIN_ZERR(2)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(MCMC_NSTEP)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = COLOR_MIN

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = COLOR_MAX

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(NCOLOR)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DMU_MIN
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DMU_MAX

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DMU_NON1A_MIN
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DMU_NON1A_MAX

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(NDMU)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = TEMPLERR_SCALE_SNIA

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = TEMPLERR_SCALE_NONIA

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_IA_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_IBC_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_II_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_PEC1A_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_MODEL1_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_MODEL2_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_MODEL3_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = PBAYES_MODEL4_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_IA_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_IBC_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_II_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_PEC1A_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_MODEL1_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_MODEL2_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_MODEL3_CUT
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = FITPROB_MODEL4_CUT

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = CHISQMIN_OUTLIER

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(NREJECT_OUTLIER)

    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(OPT_SIMCHEAT)

    DO iter=1, MXITER_PSNID
      NVAR = NVAR + 1
      INPUT_ARRAY(NVAR) = TMAX_START(iter)
      NVAR = NVAR + 1
      INPUT_ARRAY(NVAR) = TMAX_STOP(iter)
      NVAR = NVAR + 1
      INPUT_ARRAY(NVAR) = TMAX_STEP(iter)
    ENDDO

! -----------
! rate-prior info (Sep 6 2013 - RK)
    NVAR = NVAR + 1
    INPUT_ARRAY(NVAR) = DBLE(OPT_RATEPRIOR)
    DO i = 1, 3
       NVAR = NVAR + 1
       INPUT_ARRAY(NVAR) = ZRATEPRIOR_SNIA(i)
    ENDDO
    DO i = 1, 3
       NVAR = NVAR + 1
       INPUT_ARRAY(NVAR) = ZRATEPRIOR_NONIA(i)
    ENDDO

! ---------

! make INPUT_STRING
    L1 = 0

    LSTR = INDEX(NMLFILE,' ') - 1
    write(INPUT_STRING,20) '', 'NMLFILE:', NMLFILE(1:LSTR)
    L1 = L1 + LEN('NMLFILE: ')  + LSTR + 2

    LSTR = INDEX(MODELNAME_MAGERR,' ') - 1
    INPUT_STRING = INPUT_STRING(1:L1)  & 
              // 'MODELNAME_MAGERR: ' // MODELNAME_MAGERR(1:LSTR)
    L1 = L1 + LEN('MODELNAME_MAGERR: ')  + LSTR + 2


    LSTR = INDEX(FILTLIST_PEAKMAG_STORE,' ') - 1
    IF ( LSTR > 0 ) THEN
      INPUT_STRING = INPUT_STRING(1:L1)  & 
              // 'FILTLIST_PEAKMAG_STORE: '  & 
              // FILTLIST_PEAKMAG_STORE(1:LSTR)
      L1 = L1 + LEN('FILTLIST_PEAKMAG_STORE: ')  + LSTR + 2
    ENDIF

20    format(A,1x, A, 1x, A )

! ------- LOAD ARRAY AND STRING -----------------

    CALL PSNID_USER_INPUT(NVAR, INPUT_ARRAY,  & 
                 INPUT_STRING(1:L1) // char(0), L1)

    RETURN
  END SUBROUTINE PSNID_LOAD_INPUT


! =================================
    SUBROUTINE PSNID_SNLCPLOT()

! Created Feb 2013 by R.Kessler
! Determine how many best-fit functions,
! create subdir for each best-fit, and then
! call method-specific SNLCPLOT_PSNID_[METHOD]
! to make the plots
! 
! Aug 28 2017: use function DOPLOT_PSNID to determine which IPLOT to do


    USE SNDATCOM
    USE PSNIDCOM

    IMPLICIT NONE

! xxxxxxxx mark del Nov 28 2025 xxxxxx
! xx+CDE,SNFITPAR.
! xxxxxxxxxxxxxx

    INTEGER MXFUNPLOT, NFUNPLOT, iplot, LENCCID

! functions
    INTEGER  PSNID_GET_NFIT, DOPLOT_PSNID
    LOGICAL  DOPLOT_SNLC

! external functions (declared external Aug 17 2014)
    EXTERNAL  & 
         MAKEDIR_OUTPUT, CDTOPDIR_OUTPUT  & 
        ,SNLCPLOT_PSNID_BEST,  MONPLOT_PSNID_BEST, DOPLOT_PSNID  & 
        ,SNLCPLOT_PSNID_BEST2, MONPLOT_PSNID_BEST2

! ---------------- BEGIN -------------

    IF ( .NOT. DOPLOT_SNLC() ) RETURN

! determine max number of fit-functions to generate
    MXFUNPLOT = PSNID_GET_NFIT()

! ---------------------------------------

! figure out how many plots to make
    NFUNPLOT = 0
    DO iplot = 0, MXFUNPLOT-1
       IF ( DOPLOT_PSNID(iplot) > 0 ) NFUNPLOT=NFUNPLOT+1
    ENDDO
    CALL SNLCPAK_NFIT(NFUNPLOT)      ! store

! - - - - - - - - -
    DO 3000 iplot = 0, MXFUNPLOT-1

       IF ( DOPLOT_PSNID(iplot) .EQ. 0 ) GOTO 3000

       CALL MAKEDIR_OUTPUT(CCID_forC, SNLC_CID, ISNLC_LENCCID)

       IF ( METHOD_PSNID .EQ. METHOD_BEST ) THEN
          CALL SNLCPLOT_PSNID_BEST(iplot)  ! tables for light curve plots
          CALL MONPLOT_PSNID_BEST(iplot)   ! misc plots
       ELSE IF ( METHOD_PSNID .EQ. METHOD_BEST2 ) THEN
          CALL SNLCPLOT_PSNID_BEST2(iplot)  ! tables for light curve plots
          CALL MONPLOT_PSNID_BEST2(iplot)   ! misc plots
       ELSE

       ENDIF

       CALL CDTOPDIR_OUTPUT(STDOUT_UPDATE)
 3000 CONTINUE

    MADE_LCPLOT = .TRUE.  ! global flag

    RETURN
  END SUBROUTINE PSNID_SNLCPLOT


! =======================
    INTEGER FUNCTION PSNID_GET_NFIT()

! Return number of separate PSNID fits to be
! performed for each SN.
! 

    USE SNPAR
    USE PSNIDCOM

    IMPLICIT NONE

    INTEGER   NFIT
    CHARACTER C1ERR*80, C2ERR*80

! functions
    INTEGER  PSNID_BEST_GET_NFIT
    EXTERNAL PSNID_BEST_GET_NFIT

! ---------- BEGIN -----------

! set NFIT using default method so that valid NFIT is
! returned before initialization. This will be a problem
! if a future psnid option has a different NFIT.
    NFIT = PSNID_BEST_GET_NFIT()

! determine how many fit-function to generate
    IF ( METHOD_PSNID .EQ. METHOD_BEST ) THEN
       NFIT = PSNID_BEST_GET_NFIT()
    ELSE IF ( METHOD_PSNID .EQ. METHOD_BEST2 ) THEN
       NFIT = PSNID_BEST_GET_NFIT()  ! same as for BEST method
    ELSE

    ENDIF

    IF ( NFIT .EQ. 0 ) THEN
       C1ERR = 'NFIT = 0 ??'
       C2ERR = 'Check why there are no fits'
       CALL MADABORT('PSNID_GET_NFIT', C1ERR, C2ERR)
    ENDIF

    PSNID_GET_NFIT = NFIT

    RETURN
  END FUNCTION PSNID_GET_NFIT

! ==============================
    SUBROUTINE PSNID_MONANA()

! Basic monitoring of PSNID; results are reported in PSNID_MONEND
! at end of job.


    USE SNDATCOM
    USE PSNIDCOM

    IMPLICIT NONE

    INTEGER ITYPE

! function
    INTEGER  PSNID_BEST_GET_BESTTYPE, PSNID_BEST2_GET_BESTTYPE
    EXTERNAL PSNID_BEST_GET_BESTTYPE, PSNID_BEST2_GET_BESTTYPE

! ------------------- BEGIN ------------------

    ITYPE = -9
    IF ( METHOD_PSNID .EQ. METHOD_BEST ) THEN
      ITYPE_IA = 0
      ITYPE    = PSNID_BEST_GET_BESTTYPE()
    ELSE IF ( METHOD_PSNID .EQ. METHOD_BEST2 ) THEN
      ITYPE_IA = 0   ! ???
      ITYPE    = PSNID_BEST2_GET_BESTTYPE()

    ENDIF

    if ( ITYPE < 0 ) ITYPE = MXTYPE_PSNID  ! unknown

    NTYPE_PSNID(ITYPE) = NTYPE_PSNID(ITYPE) + 1

! monitor purity and effic for simulation

    IF ( LSIM_SNANA .and. SIM_TEMPLATE_INDEX .EQ. 0 ) THEN
!          generated true IA
         NSIM_GENIA = NSIM_GENIA + 1

        IF ( ITYPE .EQ. ITYPE_IA ) then
            NSIM_PSNIDIA = NSIM_PSNIDIA + 1  ! correctly typed IA
        endif
    ENDIF

    RETURN
  END SUBROUTINE PSNID_MONANA


! ==============================
    SUBROUTINE PSNID_MONEND()
! 
! Nov 30, 2012
! print TYPE-stat summary for each returned TYPE from psnid fitting.
! For simulation, include one-line summary of FOM = EFF x Purity.
! 
! Mar 31 2013: EFF->0 of there are no simulated SNIa (fix NaN bug)
! 

    USE SNDATCOM
    USE PSNIDCOM

    IMPLICIT NONE

    INTEGER ITYPE, NTOT, N, NIA_PSNID, LEN
    REAL  FRAC, EFF, PUR, FOM

    CHARACTER CTYPE*12, STRTYPE*12

! ------------------- BEGIN ------------------

    NTOT = 0
    DO itype = 0, MXTYPE_PSNID
       NTOT = NTOT + NTYPE_PSNID(itype)
    ENDDO

    IF ( NTOT .EQ. 0 ) RETURN

    CALL PRBANNER('PSNID Typing Summary')

! Now print Number and fraction for each type

    DO 200 itype = 0, MXTYPE_PSNID
       N = NTYPE_PSNID(itype)
       IF ( N .LE. 0 ) GOTO 200

       IF ( itype .EQ. MXTYPE_PSNID ) then
         STRTYPE = 'UNKNOWN'
         LEN     = INDEX(STRTYPE,' ') - 1
       ELSE
         CALL PSNID_BEST_GET_TYPENAME(ITYPE,CTYPE,LEN)
         LEN = INDEX(CTYPE,char(0)) - 1
         STRTYPE = CTYPE(1:LEN)
       ENDIF

       FRAC = FLOAT(N)/FLOAT(NTOT)
       WRITE(6,20) STRTYPE(1:8), itype, N, 100.*FRAC
20       format(T8,'N[PSNID-Type=',A,'(',I3.3,')] = ',  & 
              I6, 3x, '(', F5.1,'%)'  )

200   CONTINUE

! for simulation, include FoM(Ia) = Eff x Purity

    IF ( LSIM_SNANA ) THEN
      NIA_PSNID = NTYPE_PSNID(ITYPE_IA)

      IF ( NSIM_GENIA > 0 ) THEN
         EFF = float(NSIM_PSNIDIA)/float(NSIM_GENIA)
      ELSE
         EFF = 0.0
      ENDIF

      IF ( NIA_PSNID > 0 ) then  ! Number typed as IA
         PUR = float(NSIM_PSNIDIA)/float(NIA_PSNID)
      else
         PUR = 0.0
      endif

      print*,' '
      print*,'   Simulation Summary: '

! print raw stats to allow easy summing over log files.
      write(6,31) ' SIMGEN-Ia(ALL) ',   NSIM_GENIA
      write(6,31) ' PSNID-Ia(TRUE) ',   NSIM_PSNIDIA
      write(6,31) ' PSNID-Ia(ALL)  ',   NIA_PSNID
 31     format(T8,'N[', A, '] = ', I8 )

      FOM = EFF * PUR

      write(6,40) EFF, PUR, FOM
40      format(T8,'Simulated FOM(Ia) = EFF x PURITY = ',  & 
                F5.3,' x ', F5.3, ' = ', F5.3 )

      print*,' '
      call flush(6)

    ENDIF

    RETURN
  END SUBROUTINE PSNID_MONEND

