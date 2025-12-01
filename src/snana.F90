! F77 -> F90 translation created 2025-11-29 with command:
!   ../util/convert_snana_f77_to_f90.py snana.car





! =====================================================================
  MODULE SNPAR
    IMPLICIT NONE


! parameters used by snana code

    INTEGER, PARAMETER  ::           &  ! BEGIN SNANA PARAMS
         ONE           = 1           & 
        ,MXVERS        = 50          &  ! max number of versions to read
        ,MXSURVEY      = 100         &  ! max number of SURVEY_NAMEs to store
        ,MXITER        = 12          &  ! max # fit iterations
        ,MXSNLC        = 10000000     &  ! max number of SNe (fits format)
        ,MXIGNORE_LIST = 500         &  ! max size of IGNORE_LIST
        ,MXCID       = 499999999   &  ! max CID = 500 million -1 (9 digits)
        ,MXCID_CHECK =  99999999   &  ! MXCID to check duplicates
        ,MXEPOCH     = 30000       &  ! max number of filter-epochs per SN
        ,MXEP_MODELGRID = 200      &  ! max model grid for SNANA+SIM_MAGOBS table
        ,MXIDFIELD   = 200         &  ! max FIELD ID in SURVEY.DEF
        ,MXFIELD_OVP = 12          &  ! max number of overlapping fields
        ,MXSEASON    = 100         &  ! max number of seasons
        ,MXSNHOST    = 2           &  ! max number of host matches to read/write
        ,MXZPHOT_Q   = 20          &  ! max number of zphot quantiles (May 2022)
        ,MXVAR_PRIVATE = 40        &  ! max number of private variables
        ,MXCUT_PRIVATE = 10         &  ! max number of cuts on private var
        ,MXVAR_TERSE   = 30        &  ! max number of text-data  columms
        ,MXBIT_PHOTFLAG = 30       &  ! max number of PHOTFLAG bits
        ,MXLISTNML   = 52     &  ! max list size for some NML lists
        ,MXFILT_OBS  = 80     &  ! max number of used observer filters
        ,MXFILT_ALL  = 100    &  ! max number of all possible filter defs
        ,MXFILT_REST = 40     &  ! max number of rest-frame filter types (was 12)
        ,MXFILT_CALIB = MXFILT_REST  &  ! max number of filters used in calib table
        ,MXFILT_SNRMAX =  8        &  ! no more than this many SNRMAX(filt) cuts
        ,MNTYPE       =   1        &  ! min sn "type"
        ,MXTYPE       = 1000       &  ! max sn "type"
        ,LUNNML    = 22   &  ! LUN for input namelist
        ,LUNDAT    = 23   &  ! LUN for data
        ,LUNTMP    = 24  & 
        ,LUNOUT    = 25  & 
        ,LUNFIT    = 26  & 
        ,LUNLIST   = 27  & 
        ,LUNIGNORE = 28  & 
        ,LUNDMP    = 29  & 
        ,LUNRES1   = 31   &  ! for DMP_FITRES
        ,LUNRES2   = 32  & 
        ,LUNINTERP = 33   &  ! for FLUX,MAGs interpolated at SN,MJD
        ,LUNSALT2  = 34   &  ! for SALT2 dictFile
        ,LUNPKMJD  = 35  & 
        ,LUNCID    = 36   &  ! reserved for reading SNCID_LIST_FILE
        ,LUNSPEC   = 37   &  ! to write reformatted SPEC
        ,LUNLIST2  = 47  & 
        ,LUNIGNORE2= 48  & 
        ,ISTAGE_INIT    = 10        &  ! init has finished
        ,ISTAGE_RDSN    = 20        &  ! SN have been read
        ,ISTAGE_CUTS    = 30        &  ! cuts have been applied
        ,ISTAGE_TEST    = 90        &  ! set for testing
        ,INTERP_LINEAR  = 1     &  ! option flag to use linear DFINT
        ,INTERP_SMOOTH  = 2     &  ! option flag to use smoothing
        ,INTERP_ZSMOOTH = 3     &  ! linear DFINT, but smooth along z
        ,MXERRTYPE      = 10  & 
        ,ERRTYPE_MINOS  = 1  & 
        ,ERRTYPE_PARAB  = 2  & 
        ,ERRTYPE_MARG   = 3   &  ! marginalized error
        ,ERRTYPE_BAD    = 6  & 
        ,FCNFLAG_LAST     = 3     &  ! last MINUIT call
        ,FCNFLAG_USER     = 90    &  ! pass this IFLAG for user-FCNSNLC calls
        ,FCNFLAG_FAST     = 91    &  ! go as fast as possible (no LAST if-block)
        ,FCNFLAG_PRIOR_ONLY = 92  & 
        ,FCNFLAG_SIGMA_ONLY = 93  & 
        ,FCNFLAG_USESIM     = 99    &  ! pass this IFLAG to use SIM params
        ,FCNFLAG_MAX        = 100  & 
        ,OPT_INTEGPDF_QUITCHI2 = 2   &  ! abort FCNSNLC if chi2 > quitchi2
        ,OPT_INTEGPDF_FULL     = 1   &  ! do full FCNSNLC evaluation
        ,OPT_SNXT_CCM89   = 1   &  ! exact SN etinction using INIT_XTHOST
        ,OPT_SNXT_SJPAR   = 2   &  ! SN extinction with Jha's parameters
        ,OPT_KCORERR_SMOOTH = 1  &  ! use smooth half-Gaussian
        ,OPT_KCORERR_SJ     = 2    &  ! use Saurabh's Kcor error
        ,OPT_KCORERR_SJ5    = 5   &  ! x5 Saurabh's Kcor error
        ,OPTMASK_SNDATA_GLOBAL =  1  & 
        ,OPTMASK_SNDATA_HEAD   =  2  & 
        ,OPTMASK_SNDATA_OBS    =  4  & 
        ,OPTMASK_SNDATA_SPEC   =  8  & 
        ,OPTMASK_SNDATA_ALL    = 2+4+8   &  ! HEAD + OBS + SPEC
        ,OPTMASK_SNDATA_DONE   = 16      &  ! done reading data version
        ,OPTMASK_SNDATA_DUMP   = 32      &  ! flag to dump variable
        ,OPTMASK_SNDATA_REQUIRE= 64      &  ! flag to require variable
        ,OPT_MWCOLORLAW_DEFAULT = 99     &  ! Fitzpatrick 99
        ,OPT_MWEBV_DEFAULT      =  1     &  ! whatever is in the data file
        ,MXLINE_ARGS      = 100   &  ! max number of command line args
        ,MXKEY_ARGS       =   5   &  ! max number of args per key
        ,MXEPOCH_IGNORE   = 1000  & 
        ,NPAR_ANYLC       = 10    &  ! for MNFIT_PKMJD
        ,MXPAR_SIMSED     = 100   &  ! max number of SIMSED parameters
        ,MXPAR_LCLIB      = 40    &  ! should be same as in genmag_LCLIB.h
        ,MXLAMBIN_SNSED   = 6000  &  ! May 2024: raised from 4k -> 6k
        ,MXCUTBIT         = 64    &  ! max number of cut bits
        ,OPT_FILTUPD_EACHSN   = 1   &  ! default filter-updates
        ,OPT_FILTUPD_MAP      = 2   &  ! default filter-updates
        ,OPT_FILTUPD_SAMEFILT = 3   &  ! test: use same filter each SN
        ,OPT_FILTREST         = 1  & 
        ,OPT_FILTOBS          = 2  & 
        ,ISTAT_READAGAIN      = 7  & 
        ,ISTAT_SKIP           = -1  & 
        ,ITABLE_SNANA=1, ITABLE_FITRES=2, ITABLE_OUTLIER=3  & 
        ,ITABLE_SNLCPAK=4, ITABLE_SPECPAK=5  & 
        ,ITABLE_MODELSPEC=6, ITABLE_MARZ=7, ITABLE_DMPFCN=8  & 
        ,MXTABLE = 10  & 
        ,IDTABLE_SNANA     = 7100   &  ! anaysis variables, every event
        ,IDTABLE_FITRES    = 7788   &  ! SNANA table + fit results, passing cuts
        ,IDTABLE_OUTLIER   = 7800   &  ! outlier fluxes (Mar 2021)
        ,IDTABLE_MODELSPEC = 8000   &  ! SALT2 model spectra from LC fit
        ,IDTABLE_MARZ      = 8100   &  ! MARZ table for spectra
        ,IDTABLE_DMPFCN    = 8200   &  ! for DMPFCN (jan 2025)
        ,IDTABLE_MCMC      = 7711   &  ! obsolete ?
        ,OPT_PARSTORE_TEXTTABLE = 1   &  ! tag subset for TEXT table.
        ,MSKOPT_PARSE_WORDS_FILE    = 1    &  ! parse file
        ,MSKOPT_PARSE_WORDS_STRING  = 2    &  ! for store_PARSE_WORDS: string
        ,MSKOPT_PARSE_WORDS_IGNORECOMMA   = 4   &  ! ignore comma for string
        ,MSKOPT_PARSE_WORDS_IGNORECOMMENT = 8   &  ! ignore comment lines
        ,MSKOPT_PARSE_WORDS_FIRSTLINE     = 16  &  ! read only 1st line of file
        ,MODEL_STRETCH    = 1  & 
        ,MODEL_STRETCH2   = 2  & 
        ,MODEL_MLCS2k2    = 3  & 
        ,MODEL_SNOOPY     = 4  & 
        ,MODEL_SALT2      = 6  & 
        ,MODEL_SALT3      = 9  &  ! May 31 2019
        ,MODEL_BAYESN     = 13  &  ! Oct 7 2022: GN
        ,MODEL_SIMSED     = 7  & 
        ,MODEL_BYOSED     = 8  & 
        ,MODEL_SNEMO      = 9  & 
        ,MODEL_NON1A      = 10  & 
        ,MODEL_LCLIB      = 12  & 
        ,MODEL_FIXMAG     = 20  & 
        ,MXMODEL_INDEX    = 20  & 
        ,MXCHAR_CCID       = 20    &  ! max len of CCID string (i.e, SN name)
        ,MXCHAR_VERSION    = 72    &  ! max len of VERSION_PHOTOMETRY
        ,MXCHAR_SURVEY     = 40    &  ! max len of SURVEY_NAME
        ,MXCHAR_PATH       = 160   &  ! max len of path
        ,MXCHAR_FILENAME   = 300   &  ! max len of filename with full path
        ,MXCHAR_MODELNAME  = 72    &  ! max len of model name
        ,MXCHAR_FIELDNAME  = 20    &  ! max len of field name
        ,MXCHAR_PARNAME    = 20    &  ! max len of parameter name
        ,MXCHAR_CUTNAME    = 300   &  ! to define cut names
        ,MXCHAR_FILEWORD   =  60   &  ! size of FILEWORD_LIST
        ,MXCHAR_ARG        = 800   &  ! max lenth of command line arg
        ,MXFILE_LIST       =  10   &  ! max number of files in input list
        ,FLAG_DOCANA_START  =  1   &  ! flags DOCUMENTATION key
        ,FLAG_DOCANA_END    =  2   &  ! flags DOCUMENTATION_END
        ,FLAG_DOCANA_ERROR  = -1   &  ! flags missing DOCANA keys
        ,SNLCPAK_EPFLAG_FLUXDATA    = 1     &  ! epoch-dependent
        ,SNLCPAK_EPFLAG_FLUXMODEL   = 11  & 
        ,SNLCPAK_EPFLAG_REJECT      = 2  & 
        ,SNLCPAK_EPFLAG_CHI2        = 3  & 
        ,SNLCPAK_EPFLAG_FITFUN      = 4  & 
        ,SNLCPAK_EPFLAG_FLUXSIM     = 5     &  ! epoch-dependent
        ,SNLCPAK_EPFLAG_FLUXREST    = 6  & 
        ,SNLCPAK_EPFLAG_KCOR        = 7  & 
        ,SNLCPAK_EPFLAG_AVWARP      = 8  & 
        ,SNLCPAK_EPFLAG_SIMFLUXREST = 9  & 
        ,SNLCPAK_EPFLAG_ERRCALC     = 10  & 
        ,SNLCPAK_BANDFLAG_NDOF    = 100  & 
        ,SNLCPAK_BANDFLAG_PKFLUX  = 101   &  ! filter-dependent
        ,SNLCPAK_BANDFLAG_PKMJD   = 102  & 
        ,SNLCPAK_BANDFLAG_CHI2    = 103  & 
        ,IFLAG_INI=1, IFLAG_ANA=2, IFLAG_END=3  & 
        ,MASK_FLUXCOR_SNANA    = 1  & 
        ,MASK_FLUXERRCOR_SNANA = 2  & 
        ,EXIT_ERRCODE_SNANA = 21  & 
        ,EXIT_ERRCODE_SNFIT = 22  & 
        ,EXIT_ERRCODE_PSNID = 23  & 
        ,ERRFLAG_MNFIT_INITPAR = 91  & 
        ,ERRFLAG_MNFIT_FIXPAR  = 92  & 
        ,ERRFLAG_MNFIT_NAN     = 93  & 
        ,MXMASK_zSOURCE        = 128    ! max mask value for MASK_zSOURCE
            

   REAL*8, PARAMETER ::           &
         XTMW_FRACERR   = 0.16    &  ! default error on MW Xtinc is 16% of XTMW
        ,ZERO8          = 0.0     & 
        ,ONE8           = 1.0     & 
        ,TEN8           = 10.0    & 
        ,LOGTEN         = 2.302585092994   & 
        ,PI             = 3.1415926535898  & 
        ,PDFMIN         = 1.0E-5    &  ! used to speed up PDF integration
        ,PDFMAX_EDGE    = 0.03      &  ! max allowed PDF value at edges
        ,PDFMIN_GOOD    = 1.0E-4    &  ! PDF > PDFMIN_GOOD counts as good point
        ,MJDOFF         = 0.0       &
        ,RV_MWCOLORLAW_DEFAULT  = 3.1    &  ! A_V/E(B-V)
        ,CUTVAL_OPEN  = 1.0E12              ! cutwin value to accept everything.

   REAL, PARAMETER ::           &
        NULLVAL          = -99999.  &      
       ,MAG_SATURATE     = -7.0     &  ! for sim only
       ,LEGACY_INIT_VAL  = 1.0E8    

    CHARACTER, PARAMETER ::        &
       SNTABLE_LIST_DEFAULT*60           = 'SNANA  FITRES  LCPLOT'  &
      ,METHOD_SPLINE_QUANTILES_DEFAULT*8 = 'CUBIC'  ! default; may change in LC fit

! - - - - - - - - - - - - - - 
! physical constants

    REAL*8, PARAMETER ::                      & 
         PARSEC               = 3.085678E13   &  ! 1 parsec (km)
        ,CLIGHT               = 2.998E5       &  ! c (km/sec)
        ,H0_DEFAULT           = 70.0     &  ! standard value
        ,W0_DEFAULT           = -1.0     & 
        ,WA_DEFAULT           =  0.0     & 
        ,OMAT_DEFAULT         =  0.315   &  ! Planck 2018 (updated Mar 2020)
        ,OLAM_DEFAULT         =  0.685   &  ! idem
        ,ORAD_DEFAULT         =  1.2E-5  & 
        ,H0ERR_DEFAULT        =  7.0     & 
        ,W0ERR_DEFAULT        =  0.1     & 
        ,WAERR_DEFAULT        =  0.0     & 
        ,OMATERR_DEFAULT      =  0.03    & 
        ,OLAMERR_DEFAULT      =  0.03    & 
        ,Zat10pc              = 2.34E-9      ! magic redshift at 10 pc
          

  END MODULE SNPAR

! =====================================================================
  MODULE SNFILECOM
    USE SNPAR
    IMPLICIT NONE

    CHARACTER  & 
         SNDATA_ROOT*(MXCHAR_PATH)  & 
        ,SNANA_DIR*(MXCHAR_PATH)  & 
        ,HOST_MACHINE*(MXCHAR_PATH)      &  ! name of host machine
        ,SNDATA_PATH*(MXCHAR_PATH)       &  ! subdir with data or sim files
        ,SNLIST_FILE*(MXCHAR_FILENAME)   &  ! input list of SNDATA Files
        ,SNREADME_FILE(MXVERS)*(MXCHAR_FILENAME)    &  ! name of EVERY README file
        ,SNDATA_FILE_CURRENT*(MXCHAR_FILENAME)      &  ! current file being read
        ,GLOBAL_BANNER*120  & 
        ,SNDATA_PREFIX*(MXCHAR_FILENAME)   &  ! $SNDATA_ROOT/lcmerge/$VERSION
        ,C1ERR*88, C2ERR*88    ! generic error strings

    LOGICAL LFLAG_RDHEAD_ONLY


  END MODULE SNFILECOM

! =====================================================================
  MODULE CTRLCOM
    USE SNPAR
    IMPLICIT NONE

! control variables and counters.

    INTEGER  & 
         ISTAGE_SNANA        &  ! current stage of processing
        ,NACCEPT_CUT(MXCUTBIT)  &  ! Number of SN that pass each cut
        ,NACCEPT_CID            &  ! # SN with valid CID
        ,NACCEPT_TYPE           &  ! # SN with valid TYPE
        ,NACCEPT_Z              &  ! idem with valid redshift
        ,NACCEPT_ZERR           &  ! idem with valid redshift error
        ,NCALL_SNANA_DRIVER  & 
        ,NCALL_FCNFLAG(FCNFLAG_MAX)  & 
        ,NPASSCUT_INCREMENT(-1:MXTYPE,100)   &  ! 100 > NCUTBIT_SNLC
        ,NPASSCUT_FIT(-1:MXTYPE)  & 
        ,NSTORE_MAGCOR          &  ! number of stored MAGCOR values
        ,NUSE_MAGCOR            &  ! number of used MAGCOR values
        ,SIGN_MAGCOR            &  ! add or subtract
        ,FORCEMASK_FLUXCOR    &  ! mask to force fluxCor, even if already applied
        ,EXIT_ERRCODE        ! used for abort

    INTEGER*8  & 
         JTIME_START  & 
        ,JTIME_LOOPSTART         &  ! time at start of fits with TIME()
        ,JTIME_LOOPEND          ! time and end

    LOGICAL  & 
         DO_FIT  & 
        ,DO_GETINFO         &  ! T=>GETINFO flag to print and quit
        ,DO_FLUXERRCALC     &  ! T => compute error from PSF,SKY   ZPT
        ,LSIM_SNANA         &  ! simulated with SNANA
        ,LSIM_MAGOBS    &  ! data-like, but with SIM_MAGOBS (e.g. fakes on images)
        ,ISJOB_SNANA           &  ! =T for snana.exe only.
        ,ISJOB_SNFIT           &  ! =T for snlc_fit or psnid
        ,ISJOB_PSNID  & 
        ,ISJOB_SIM     &  ! either SNANA sim or fakes (LSIM_SNANA or LSIM_MAGOBS)
        ,ISJOB_BATCH   &  ! true if JOBSPLIT passed on command line (Aug 2020)
        ,REDUCE_STDOUT_BATCH    &  ! T => suppress stdout for JOBSPLIT(1)>1
        ,REFORMAT               &  ! T => a reformat option has been set
        ,REFORMAT_SNANA         &  ! T => reformat SNANA format
        ,REFORMAT_SAVE_BADEPOCHS  &  ! T => include reject epochs (default=F)
        ,REFORMAT_BAND_NAME     &  ! T => [band] -> [SURVEY]-[band]
        ,REFORMAT_SPECTRA_INCLUDE  &  ! include spectra in data reformat (default=F)
        ,REFORMAT_SPECTRA_ONLY     &  ! reformat only spectra (OPT_REFORMAT_SPECTRA>0)
        ,REFORMAT_PRIVATE       &  ! T-> include private variables (default=T)
        ,REFORMAT_SIMTRUTH      &  ! T -> include SIM truth (default=T)
        ,STDOUT_UPDATE          &  ! T => update event to screen
        ,DOFUDGE_HOSTNOISE      &  ! T => FUDGE_HOSTNOISE_FILE is set
        ,DOFUDGE_NONLIN         &  ! T => NONLINEARITY_FILE is set
        ,DOFUDGE_FLUXERRMODEL   &  ! T => FLUXERRMODEL_FILE
        ,DOzSHIFT               &  ! T => apply a zSHIFT systematic
        ,ISCORRECT_SIGN_VPEC    &  ! T => data has correct VPEC sign (10/2020)
        ,DOFIX_WRONG_SIGN_VPEC  &  ! T => flip VPEC sign
        ,UNIT_PSF_NEA           &  ! T => read NEA from data; else it's PSF
        ,FOUND_ATMOS            &  ! T=> found atmos variables (dRA,dDEC,AIRMASS)
        ,USE_SNCID_FILE          &  ! bit0 of OPT_SNCID_LIST
        ,USE_INIVAL_SNCID_FILE   &  ! bit1 of OPT_SNCID_LIST
        ,USE_PRIOR_SNCID_FILE   ! bit2 of OPT_SNCID_LIST

    CHARACTER  & 
          SNANA_VERSION*60         &  ! e.g., v11_04h-[commit-id]
         ,SNANA_VERSION_DATA*60    &  ! SNANA version used create FITS data
         ,REFORMAT_VERSION*(MXCHAR_VERSION)  & 
         ,DATATYPE*12  ! e..g, 'DATA', 'FAKE', 'SNANA_SIM'

! global survey info
    CHARACTER  & 
         SURVEY_NAME*(MXCHAR_SURVEY)  & 
        ,SURVEY_NAME_LIST(MXSURVEY)*(MXCHAR_SURVEY)  &  ! in SURVEY.DEF file
        ,SUBSURVEY_NAME*(MXCHAR_SURVEY)   &  ! e.g.,  SURVEY:  BLA(SUBSURVEY)
        ,SUBSURVEY_NAME_LIST*(MXCHAR_FILENAME)  &  ! comma-sep list
        ,SURVEY_FILTERS*(MXFILT_ALL)   &  ! read from data file
        ,SURVEY_FIELDNAME(MXIDFIELD)*(MXCHAR_FIELDNAME) ! from SURVEY.DEF

    INTEGER  & 
         NFIELD_SURVEY              &  ! number of survey fields in SURVEY.DEF
        ,SURVEY_IDFIELD(MXIDFIELD) ! integer ID for each field

    REAL  & 
         ZP_FLUXCAL  ! defines calibrated flux

    INTEGER*4  & 
         N_VERSION              &  ! number of photometry version to read
        ,N_SNLC_READ(MXVERS)    &  ! total numbrer of LC per version
        ,N_SNLC_PROC            &  ! number of SNLC processed with SNANA_DRIVER
        ,N_SNLC_CUTS            &  ! Number of SN after cuts (bookkeeping only)
        ,N_SNLC_SPEC            &  ! idem, but with 1 or more spectra
        ,N_SNLC_FIT             &  ! Number of fitted SN
        ,N_SNLC_FITCUTS         &  ! Number of fitted SN after fit cuts
        ,N_SNLC_COVFIX          &  ! Number of SN with fixed COV to be invertible
        ,N_SNHOST_ZSPEC         &  ! for YAML out file
        ,N_SNHOST_ZPHOT         &  ! idem
        ,N_MASK_zSOURCE_LC_CUTS(0:MXMASK_zSOURCE)      &  ! Nevt per bit mask (after SNANA/LC cuts)
        ,N_MASK_zSOURCE_LCFIT_CUTS(0:MXMASK_zSOURCE)   &  ! Nevt per bit mask (after LCFIT   cuts)
        ,N_DUPLICATE_CID        &  ! Number of duplicate CIDs
        ,N_DUPLICATE_MJD        &  ! Number of duplicate MJD+BAND (Jun 2017)
        ,NSTORE_DUPLICATE_MJD   &  ! Number stored
        ,N_SNFILE               &  ! # of SN files to read per version
        ,N_SNFILE_LAST          &  ! idem, as of last version
        ,ABSO_OFFSET            &  ! used to compute ABSO_INDEX
        ,ABSO_INDEX             &  ! absolute index (row or file number)
        ,NEPOCH_TOT             &  ! total number of epochs read from files
        ,NEPOCH_USE             &  ! total number of used epochs
        ,NEPOCH_CUT             &  ! total number of epochs passing cuts
        ,NEPOCH_BADPHOT         &  ! # epochs with bad PHOTFLAG (per event)
        ,NEPOCH_BADPHOT_SUM     &  ! # epochs with bad PHOTFLAG (summed)
        ,IDSURVEY               &  ! survey ID from SURVEY.DEF
        ,IDSUBSURVEY            &  ! =IDSURVEY unless subSurvey is different
        ,IDSURVEY_LIST(MXSURVEY)  &  ! corresponds to SURVEY_NAME_LIST
        ,NSURVEY_LIST          ! size of SURVEY_NAME_LIST   IDSURVEY_LIST

    LOGICAL*1  & 
          EXIST_CALIB_FILE  & 
         ,EXIST_FILT(MXFILT_OBS)   &  ! T => at least one point per filt
         ,FOUND_SURVEY  & 
         ,FORMAT_TEXT     &  ! ascii/txt for input data
         ,FORMAT_FITS    ! snfitsio for input data.

    INTEGER N_SNLC_PLOT
    LOGICAL MADE_LCPLOT  ! SAVE:  T if LC plot was made




! logical *1 stuff


! SNTABLE control variables (May 2014)
    LOGICAL  & 
         USE_TABLEFILE_ROOT    &  ! write table in ROOT format
        ,USE_TABLEFILE_TEXT    &  ! write table in TEXT format
        ,USE_TABLEFILE_MARZ    &  ! write table in MARZ/FITS format
        ,USE_TABLEFILE         &  ! T if either any of the above are set
        ,WRTABLEFILE_IAUC      &  ! T-> add IAUC column(for writing tables)
        ,WRTABLEFILE_SIMVAR    &  ! T-> include SIM_XXX vars for simulation
        ,WRTABLEFILE_ZPHOT     &  ! T-> include ZPHOT info for PHOTOZ fit
        ,WRTABLEFILE_HOST_TEXT   &  ! T-> add HOST info in SNANA TEXT table
        ,WRTABLEFILE_HOST2_TEXT  &  ! T-> add HOST and HOST2 info in SNANA TEXT table
        ,WRTABLEFILE_ERRCALC_TEXT! T-> add ERRCALC info in LCPLOT TEXT table

    INTEGER  & 
         OPT_TABLE(MXTABLE)        &  ! option(s) for each table
        ,CUTMASK_SNANA_TABLE      ! select CUTFLAG_SNANA for SNANA table

    REAL  & 
         PRESCALE_TABLE(MXTABLE)  ! table pre-scale (to reduce size)

    CHARACTER  & 
         TEXTFORMAT_TABLE(MXTABLE)*8

! arrays for TABLE_FILTER_REMAP (Feb 2017)
    INTEGER NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE(MXFILT_ALL)
    CHARACTER FILTLIST_REMAP_TABLE*(MXFILT_ALL)

! Sep 23 2017: allow user codes to add private variables to SNTABLE
    INTEGER   NTABLEVAR_USER  ! e.g., set in snana_private.cra
    CHARACTER TABLEVARNAME_USER(MXVAR_PRIVATE)*40
    REAL      TABLEVALUE_USER(MXVAR_PRIVATE)


! May 6, 2008: define lists for epochs to ignore

    INTEGER NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
    CHARACTER  & 
         EPOCH_IGNORE_CCID(MXEPOCH_IGNORE)*(MXCHAR_CCID)  & 
        ,EPOCH_IGNORE_FILT(MXEPOCH_IGNORE)*4  & 
        ,EPOCH_IGNORE_LASTFILE*(MXCHAR_FILENAME)

    REAL*8  & 
         EPOCH_IGNORE_MJD(MXEPOCH_IGNORE)




    REAL*8  DUPLICATE_MJDLIST(200)

  END MODULE CTRLCOM

! =====================================================================
  MODULE OUTLIERCOM
    USE SNPAR
    IMPLICIT NONE

    INTEGER, PARAMETER ::  MXFILT_OUTLIER = 100 ! avoid includig FILTCOM

    CHARACTER OUTLIER_TABLE_NAME*12

    REAL  & 
         NSIGCUT_OUTLIER     &  ! e.g., SNTABLE_LIST='OUTLIER(nsig:3.5)'
        ,FTRUECUT_OUTLIER    &  ! e.g., SNTABLE_LIST='OUTLIER(Ftrue:0.0)'
        ,SBMAGCUT_OUTLIER   ! e.g., SNTABLE_LIST='OUTLIER(sbmag:25.0)

    REAL*8 OUTLIER_MJD
    REAL*4  & 
         OUTLIER_FLUXCAL_DATA, OUTLIER_FLUXCAL_ERR_DATA  & 
        ,OUTLIER_FLUXCAL_ERR_CALC  & 
        ,OUTLIER_FLUXCAL_TRUE  & 
        ,OUTLIER_FLUXCAL_FIT, OUTLIER_FLUXCAL_ERR_FIT  & 
        ,OUTLIER_FLUXCAL_ERRTOT_FIT   &  ! data + fit errors in quadrature
        ,OUTLIER_AREAFRAC  & 
        ,OUTLIER_ZP, OUTLIER_NEA, OUTLIER_PSF, OUTLIER_SKYSIG  & 
        ,OUTLIER_NSIG, OUTLIER_LOGSNR, OUTLIER_SBMAG  & 
        ,OUTLIER_z, OUTLIER_TREST  ! for LC FIT only

    INTEGER  & 
          NEP_OUTLIER_TOT, NEP_OUTLIER_PEREVT  & 
         ,NEP_OUTLIER_CHECK(MXFILT_OUTLIER)   &  ! total epochs checked per band
         ,NEP_OUTLIER_FOUND(MXFILT_OUTLIER)   &  ! Noutlier per band
         ,OUTLIER_IFILTOBS, OUTLIER_DETNUM, OUTLIER_PHOTFLAG  & 
         ,OUTLIER_IXPIX, OUTLIER_IYPIX

    CHARACTER OUTLIER_BAND*4, OUTLIER_FIELD*20




! efficiency for PHOTFLAG bits

    INTEGER  & 
         NEP_OUTLIER_VS_PHOTBIT(2,MXBIT_PHOTFLAG)     &  ! 1=TOT, 2=bit set
        ,NEP_ALL_VS_PHOTBIT(2,MXBIT_PHOTFLAG)        ! 1=TOT, 2=bit set


  END MODULE OUTLIERCOM

! =====================================================================
  MODULE USRTAGCM
    USE SNPAR
    IMPLICIT NONE

    INTEGER, PARAMETER :: MXUSERTAG = 1000  ! max number of user tags

    INTEGER  & 
         N_USERTAGS            &  ! number of tags in USERTAGS_FILE
        ,USERTAG_VALUELIST(MXUSERTAG)   &  ! list of user tags
        ,USERTAG                   ! USERTAG value  for current SN

    BYTE  & 
         USERTAG_USED(MXUSERTAG)   ! mark when each USERTAG is used

    CHARACTER  & 
         USERTAG_CCIDLIST(MXUSERTAG)*(MXCHAR_CCID)  ! list of CCID



  END MODULE USRTAGCM

! =====================================================================
  MODULE REQEPCOM
    USE SNPAR
    IMPLICIT NONE
    CHARACTER FILTLIST_REQEP*(MXFILT_ALL)
    INTEGER NFILT_REQEP, IFILTLIST_REQEP(MXFILT_ALL)
    REAL  TRANGE_REQEP(3), SNRMIN_REQEP
    LOGICAL ISFRAME_REST_REQEP, ISFRAME_OBS_REQEP


    REAL NDAYS_ABOVE_SNRMIN_REQEP


  END MODULE REQEPCOM

! =====================================================================
  MODULE EARLYCOM
    USE SNPAR
    IMPLICIT NONE

    REAL, PARAMETER :: DT_SAMENIGHT = 0.400 ! same night within 0.4 day

! variables for selecting the early part of a light curve
    INTEGER    NOBS_EARLYLC, NNIGHT_EARLYLC
    INTEGER    NPHOTMASK_START_EARLYLC
    REAL       NSNR_START_EARLYLC
    REAL       MJDLAST_EARLYLC, MJDLAST_SELECT


! variables parsed from EARLYLC_STRING:
    INTEGER  MAXOBS_EARLYLC, MAXNIGHT_EARLYLC, PHOTMASK_EARLYLC,  & 
               NDAYADD_EARLYLC, PHOTMASK_START_EARLYLC
    REAL  SNRMIN_EARLYLC, PHOTPROBMIN_EARLYLC, SNR_START_EARLYLC
    CHARACTER  FILTERS_EARLYLC*(MXFILT_ALL)


  END MODULE EARLYCOM

! =====================================================================
  MODULE INTERPCM
    USE SNPAR
    IMPLICIT NONE

! arrays for interpolation

    INTEGER, PARAMETER :: MXINTERP = 1000 ! ->1000 Nov 19 2019 (was 100)

    INTEGER*4  & 
         N_INTERP_MJDLIST        ! # SN MJDs to interpolate

    REAL*8  & 
         INTERP8_MJDLIST(MXINTERP)

    CHARACTER  & 
         INTERP_CCIDLIST(MXINTERP)*(MXCHAR_CCID)

    INTEGER N_INTERP_MJDLIST_DONE
    LOGICAL INTERP_MJDLIST_DONE(MXINTERP)



  END MODULE INTERPCM

! =====================================================================
  MODULE SNLCCOM
    USE SNPAR
    IMPLICIT NONE

! Main common block variables for light curves and
! analysis-related variables.

    INTEGER  & 
         SNLC_CID                   &  ! integer CID
        ,INDEX_CID_MATCH            &  ! ISN index from SNCID_LIST_FILE
        ,ISNLC_SNRECON_USE(MXEPOCH)   &  ! 1 -> epoch used in SNRECON
        ,ISNLC_PHOTFLAG(MXEPOCH)      &  ! photomety flags
        ,ISNLC_zFLAG                  &  ! zspec-quality flag (May 2020), data only
        ,ISNLC_zSOURCE                &  ! mask indicating the source of redshift
        ,ISNLC_LENCCID      &  ! char-len of CCID
        ,ISNLC_LENIAUC      &  ! char-len of IAUC name
        ,ISNLC_LENNAME      &  ! char-len of transient name
        ,ISNLC_VERSION      &  ! photometry version index vs. ISN
        ,ISNLC_TYPE         &  ! type (120=confirmed Ia, etc ...)
        ,ISNLC_IFILE        &  ! ifile index, text format only
        ,ISNLC_IFILT_OBS(MXEPOCH)  &  ! filt-index for each epoch/SN
        ,ISNLC_NEWMJD_HEAD    &  ! Number of NEWMJDs in header
        ,ISNLC_NEWMJD_FOUND   &  ! Number of NEWMJDs found in file
        ,ISNLC_NEWMJD_STORE   &  ! Number of NEWMJDs stored
        ,ISNLC_NEWMJD_CUTS    &  ! Number of NEWMJDs after cuts
        ,ISNLC_EPOCH_RANGE_NEWMJD(2,MXEPOCH)  & 
        ,ISNLC_NFILT_NEWMJD(MXEPOCH)    &  ! # observed filters per NEWMJD
        ,ISNLC_NFILT_SNRMAX       &  ! number of filters passing SNR cut
        ,ISNLC_NFILT_SNRMAX2      &  ! idem for 2nd SNRMAX cut
        ,ISNLC_NFILT_TRESTMIN    &  ! Nfilt passing TRESTMIN cut
        ,ISNLC_NFILT_TRESTMAX    &  ! idem for TRESTMAX
        ,ISNLC_NFILT_TREST2      &  ! Nfilt passing CUTWIN_TREST2 cut
        ,ISNLC_NEPOCH_FOUND      &  ! actual number of epochs found
        ,ISNLC_NEPOCH_STORE      &  ! number of epochs stored in memory
        ,ISNLC_NEPOCH_USE        &  ! used after PHOTMASK cuts
        ,ISNLC_NEPOCH_PHOTPROB   &  ! NEPOCH with PHOTPROB >= 0
        ,ISNLC_NEPOCH_FILT(MXFILT_OBS)   &  ! NEPOCH vs. filter
        ,ISNLC_NOBS_DETECT       &  ! NOBS with detection
        ,ISNLC_NOBS_PREDETECT    &  ! NOBS in CUTWIN_TOBS_PREDETECT window
        ,ISNLC_FAKE              &  ! => real data, else it's a fake
        ,ISNLC_DETNUM(MXEPOCH)   &  ! read from header (May 2017)
        ,ISNLC_IMGNUM(MXEPOCH)   &  ! Oct 2021: IMAGE NUNBER (e.g., EXPNUM, VISIT_ID)
        ,ISNLC_IDFIELD(MXEPOCH)  &  ! integer field id
        ,ISNLC_NFIELD_OVP        &  ! number of fields (>=2 for overlap)
        ,ISNLC_CUTFLAG_REQEP         &  ! idem
        ,ISNLC_NOBS_TREST            &  ! mask for SNCUT_NOBS_TREST (Nov 2025)
        ,ISNLC_CUTFLAG_PRIVATE       &  ! idem for private var cuts
        ,ISNLC_CUTFLAG_SIMVAR        &  ! idem for SIMVAR cuts
        ,ISNLC_WRMASK_FLUXCOR_SNANA    &  ! write if fudges applied to data
        ,ISNLC_RDMASK_FLUXCOR_SNANA   ! read if SNANA fudges applied to data


    REAL*8  & 
         SNLC8_RA         &  ! RA   vs. SN
        ,SNLC8_DEC        &  ! DEC  vs. SN
        ,SNLC8_MJD(MXEPOCH)   &  ! MJD for each epoch
        ,SNLC8_MJDMIN         &  !   min MJD among all measurements
        ,SNLC8_MJDMAX         &  !   max MJD ...
        ,SNLC8_MJD_TRIGGER    &  ! MJD of trigger bit (see PHOTFLAG_TRIGGER)
        ,SNLC8_MJD_DETECT_FIRST  & 
        ,SNLC8_MJD_DETECT_LAST  & 
        ,SNLC8_FLUXCAL(MXEPOCH)  &  ! to call functions requiring double
        ,SNLC8_FLUXCAL_ERRTOT(MXEPOCH)

    REAL  & 
         SNLC_ZHELIO         &  ! redshift used for final analysis
        ,SNLC_ZHELIO_ERR     &  ! error on above
        ,SNLC_ZCMB           &  ! redshift used for final analysis
        ,SNLC_ZCMB_ERR       &  ! error on above
        ,SNLC_ZHD            &  ! redshift used for Hubble diagram
        ,SNLC_ZHD_ERR        &  ! error on above
        ,SNLC_ZSN            &  ! redshift of SN only (ignoring host-z)
        ,SNLC_ZSN_ERR        &  ! error on above
        ,SNLC_REDSHIFT       &  ! redshift used for light curve fit
        ,SNLC_REDSHIFT_ERR   &  ! error on above
        ,SNLC_VPEC           &  ! pec. velocity
        ,SNLC_VPEC_ERR       &  ! error on above
        ,SNLC_LENSDMU        &  ! measured LENSDMU  from WL (Jan 2025)
        ,SNLC_LENSDMU_ERR    &  ! error on measurement
        ,SNLC_ZPEC           & 
        ,SNLC_ZPEC_ERR       & 
        ,SNLC_Trestmin       &  ! earliest epoch, rest frame days since peak
        ,SNLC_Trestmax       &  ! latest   epoch, rest frame days since peak
        ,SNLC_TrestRange     & 
        ,SNLC_Tobsmin        &  ! earliest epoch, obs frame days since peak
        ,SNLC_Tobsmax        &  ! latest   epoch, obs frame days since peak
        ,SNLC_TGAPMAX        &  ! max gap within TREST-range
        ,SNLC_T0GAPMAX       &  ! max gap near peak
        ,SNLC_SNRMAX_FILT(0:MXFILT_OBS)    &  ! max S/N per filter/SN
        ,SNLC_SNRMAX_SORT(MXFILT_OBS)      &  ! 1st, 2nd ... SNRMAX by filt
        ,SNLC_SNRSUM                       &  ! quadrature sum of SNR 
        ,SNLC_FLUXCALMAX(MXFILT_OBS)       &  ! max flux per filter/SN
        ,SNLC_FLUXCALMAX_ERR(MXFILT_OBS)   &  ! uncertainty on above
        ,SNLC_SNANAFIT_PEAKMJD             &  ! SNANA-estimate of PEAKKMJD
        ,SNLC_SNANAFIT_PEAKMJD_FITPAR(MXFILT_OBS,NPAR_ANYLC)  & 
        ,SNLC_SNANAFIT_PEAKMJD_FITERR(MXFILT_OBS,NPAR_ANYLC)  & 
        ,SNLC_PHOTPROB(MXEPOCH)    &  ! generic 'fit probability' per epoch
        ,SNLC_PHOTPROB_MIN         &  ! min photprob for PHOTPROB>0
        ,SNLC_TOBS(MXEPOCH)        &  ! MJD-SET_PEAKMJD
        ,SNLC_TREST(MXEPOCH)       &  !  MJD-SET_PEAKMJD)/(1+z)
        ,SNLC_GAIN(MXEPOCH)       &  ! e/AUD
        ,SNLC_RDNOISE(MXEPOCH)    &  ! read noise per pix, e-
        ,SNLC_PIXSIZE             &  ! pixel size
        ,SNLC_NXPIX               &  ! total number of X-pixels (Aug 7 2014)
        ,SNLC_NYPIX               &  ! total number of Y-pixels
        ,SNLC_XPIX(MXEPOCH)       &  ! pixel location
        ,SNLC_YPIX(MXEPOCH)       &  ! pixel location
        ,SNLC_AREAFRAC(MXEPOCH)   &  ! area-frac contained by XPIX,YPIX
        ,SNLC_AREAFRAC_AVG        &  ! average over epochs (May 2020)
        ,SNLC_dRA(MXEPOCH)        &  ! RA(obs) - RA_AVG(band)
        ,SNLC_dDEC(MXEPOCH)       & 
        ,SNLC_AIRMASS(MXEPOCH)    & 
        ,SNLC_MWEBV               &  ! Milky Way Galactic E(B-V)
        ,SNLC_MWEBV_ERR           &  ! error on above
        ,SNLC_SKYSIG(MXEPOCH)     &  ! sigma on above
        ,SNLC_SKYSIG_T(MXEPOCH)   &  ! sigma on template run
        ,SNLC_PSF_SIG1(MXEPOCH)   &  ! sigma, pixels
        ,SNLC_PSF_SIG2(MXEPOCH)   & 
        ,SNLC_PSF_RATIO(MXEPOCH)  & 
        ,SNLC_PSF_FWHM_ARCSEC(MXEPOCH)    & 
        ,SNLC_PSF_NEA(MXEPOCH)            &  ! Noise-equivalent area (pixels)
        ,SNLC_SNR(MXEPOCH)                &  ! SNR
        ,SNLC_FLUXCAL_OFF(MXFILT_OBS)     &  ! add SN light from template
        ,SNLC_FLUXCAL_ERRCALC(MXEPOCH)    &  ! calc a-la simulation
        ,SNLC_FLUXCAL_HOSTERRCALC(MXEPOCH)   &  ! idem for host error
        ,SNLC_FLUXCAL(MXEPOCH)               & 
        ,SNLC_FLUXCAL_ERRTOT(MXEPOCH)        & 
        ,SNLC_FLUXCAL_ERRTEST(MXEPOCH)    &  ! ERRCALC/ERRTRUE (Nov 2019
        ,SNLC_FLUXCAL_FIT(MXEPOCH)        &  ! used for outlier table
        ,SNLC_FLUXCAL_ERR_FIT(MXEPOCH)    &  ! idem
        ,SNLC_FLUXCAL_ERRTOT_FIT(MXEPOCH) &  ! sqrt[cov(data)+cov(fit)]
        ,SNLC_MAG(MXEPOCH)                & 
        ,SNLC_MAG_ERRPLUS(MXEPOCH)        & 
        ,SNLC_MAG_ERRMINUS(MXEPOCH)       & 
        ,SNLC_ZEROPT(MXEPOCH)             & 
        ,SNLC_ZEROPT_ERR(MXEPOCH)         & 
        ,SNLC_ZEROPT_forCUT(MXEPOCH)      &  ! depends on CUTWIN_ZPADU or CUTWIN_ZPNPE
        ,SNLC_TEXPOSE(MXEPOCH)            & 
        ,SNLC_DLMAG                       &  ! 5*log10(10pc/DL)
        ,SNLC_SKYFLUXCAL(MXEPOCH)         &  ! calculated sky fluxcal/pixel
        ,SNLC_MWXT_MAG(MXFILT_OBS)        &  ! mag stellar extinct
        ,SNLC_MWXT_FLUXFRAC(MXFILT_OBS)   &  ! same for flux (< 1)
        ,SNLC_MWXT_MAGERR(MXFILT_OBS)     &  ! Galactic mag err per filter
        ,SNLC_SEARCH_PEAKMJD              &  ! external PEAKMJD
        ,SNLC_DTOBS(MXEPOCH)              &  ! time since last obs
        ,SNLC_DTOBS_SAMEFILT(MXEPOCH)     &  ! idem, but same filter
        ,SNLC_TLIVE_DETECT     ! MJD(last detection) - MJD(1st detection)


    INTEGER  & 
         NSEASON_TOT      &  ! total number of seasons
        ,NSEASON_ACTIVE  ! total number of active seasons

    REAL*4  & 
         MULTISEASON_CHI2RED(MXSEASON)  &  ! borrow another 'MX' param
        ,MULTISEASON_AVGFLUX(MXSEASON)  & 
        ,MULTISEASON_MJDMIN(MXSEASON)  & 
        ,MULTISEASON_MJDMAX(MXSEASON)

    CHARACTER  & 
         SNLC_CCID*(MXCHAR_CCID)            &  ! public integer or string name
        ,SNLC_NAME_IAUC*(MXCHAR_CCID)       &  ! public IAUC name
        ,SNLC_NAME_TRANSIENT*(MXCHAR_CCID)  &  ! internal transient name (July 2024)
        ,SNLC_CASTCID*8              &  ! = 'INT' or 'CHAR' to indicate cast
        ,SNLC_FIELD(MXEPOCH)*(MXCHAR_FIELDNAME)     &  ! char field name
        ,SNLC_FIELD_OVPLIST(MXFIELD_OVP)*(MXCHAR_FIELDNAME)   &  ! sparse list of overlap fields
        ,SNLC_FIELDLIST*(MXCHAR_FIELDNAME) ! e.g., '82S', '82N+82S'




! - - - - -  ADDCOL stuff - - - - - - -
! ADDCOL arrays are loaded with original filter indices,
! or with REMAPed filter indices. Allows using unique
! pointer for calls to SNTABLE_ADDCOL.

    CHARACTER  & 
          ADDCOL_FILTERS*(MXFILT_ALL) ! SURVEY_FILTERS

    REAL  & 
          ADDCOL_SNHOST_MAGOBS(MXFILT_ALL,MXSNHOST)  & 
         ,ADDCOL_SNHOST_SBFLUXCAL(MXFILT_ALL)  & 
         ,ADDCOL_SNHOST_SBMAG(MXFILT_ALL)  & 
         ,ADDCOL_FLUXCALMAX(MXFILT_ALL)  & 
         ,ADDCOL_FLUXCALMAX_ERR(MXFILT_ALL)  & 
         ,ADDCOL_CHI2_FITPKMJD(MXFILT_ALL)  & 
         ,ADDCOL_SNRMAX(MXFILT_ALL)  & 
         ,ADDCOL_XTMW(MXFILT_ALL)  & 
         ,ADDCOL_PROB_TRUEFLUX(MXFILT_ALL)

    INTEGER  & 
          ADDCOL_NDOF_FITPKMJD(MXFILT_ALL)  & 
         ,ADDCOL_NDOF_TRUEFLUX(MXFILT_ALL)



! - - - -- MULTI-SEASON STUFF - - - -

  END MODULE SNLCCOM

! =====================================================================
  MODULE SPECCOM
    USE SNPAR
    IMPLICIT NONE
!  Apr 2019
!  Read spectra to pass into tables for plotting.
!  Not intended for analysis.

    INTEGER, PARAMETER ::         & 
         MXVAR_SPECTRUM = 10      &  ! max number of variables after SPEC: key
        ,MXSPECTRUM     = 200     &  ! max number of spectra to read per event
        ,MXLAM_SPECTRUM = 20000   &  ! max number of wave bins per spectrum
        ,MXVAR_SPEC_TABLE = 20       ! max number of columns in summary spec table
              

    INTEGER NVAR_SPEC_TABLE  & 
         ,IVAR_SPEC_SUBTYPE, IVAR_SPEC_MJD, IVAR_SPEC_TOBS  & 
         ,IVAR_SPEC_TREST, IVAR_SPEC_NBIN

    CHARACTER  & 
         SPEC_TABLE_FILE*40  & 
        ,VARNAMES_SPEC_TABLE(MXVAR_SPEC_TABLE)*40   &  ! for summary table, not individual spectrum files
        ,COMMENTS_SPEC_TABLE(MXVAR_SPEC_TABLE)*80  ! comment string for each VARNAME
    LOGICAL USECOL_SPEC_TABLE(MXVAR_SPEC_TABLE)  ! logical flag

    LOGICAL WRSPEC_FLAM, WRSPEC_FLUX, WRSPEC_KEY, WRSPEC_SNID_SAGE

    INTEGER  & 
          NVAR_SPECTRUM            &  ! number of spectrum variables to read
         ,NSPECTRUM                     &  ! Number of spectra
         ,ID_SPECTRUM(MXSPECTRUM)       &  ! ID for each spectrum
         ,NLAMBIN_SPECTRUM(MXSPECTRUM)  &  ! Number of wave bins per spectrum
         ,NLAMBIN_READ                 ! keep track of wave index

    REAL*8  & 
          MJD_SPECTRUM(MXSPECTRUM)  & 
         ,TOBS_SPECTRUM(MXSPECTRUM)  & 
         ,TEXPOSE_SPECTRUM(MXSPECTRUM)  & 
         ,LAMMIN_SPECTRUM(MXLAM_SPECTRUM)  & 
         ,LAMMAX_SPECTRUM(MXLAM_SPECTRUM)  & 
         ,FLAM_SPECTRUM(MXLAM_SPECTRUM)  & 
         ,FLAMERR_SPECTRUM(MXLAM_SPECTRUM)  & 
         ,SIM_FLAM_SPECTRUM(MXLAM_SPECTRUM)




  END MODULE SPECCOM

! =====================================================================
  MODULE SNSIMCOM
    USE SNPAR
    IMPLICIT NONE

! simulation parameters (if FAKE=2)

    REAL*8  SIM8_RA, SIM8_DECL

    REAL  & 
         SIM_REDSHIFT_HELIO  & 
        ,SIM_REDSHIFT_CMB  & 
        ,SIM_REDSHIFT_HOST  & 
        ,SIM_REDSHIFT_HOST_MATCH  & 
        ,SIM_VPEC  & 
        ,SIM_DLMAG  & 
        ,SIM_LENSDMU  & 
        ,SIM_MUSHIFT  & 
        ,SIM_MWEBV  & 
        ,SIM_MWRV  & 
        ,SIM_WGT_POPULATION      &  ! Dec 2023: for biasCor sims only
        ,SIM_SALT2x0  & 
        ,SIM_SALT2mb  & 
        ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV  & 
        ,SIM_SHAPEPAR, SIM_SHAPELAW  & 
        ,SIM_PEAKMJD, SIM_MJD_EXPLODE  & 
        ,SIM_EXPOSURE_TIME(MXFILT_OBS)   &  ! relative exposure time
        ,SIM_PEAKMAG(MXFILT_OBS)  & 
        ,SIM_EPMAGOBS(MXEPOCH)       &  ! true epoch mag at each filter/epoch
        ,SIM_EPFLUXCAL(MXEPOCH)      &  ! true fluxcal at each filter/epoch
        ,SIM_EPSNRMON(MXEPOCH)       &  ! optional SNR at MAGMONITOR
        ,SIM_EPMAGREST(MXEPOCH)      &  ! true rest-frame mag
        ,SIM_EPFLUXCAL_HOSTERR(MXEPOCH)  &  ! true error from host noise
        ,SIM_EPCHI2FLUX(MXEPOCH)     &  ! data-sim chi2 per epoch
        ,SIM_EPPULL(MXEPOCH)         &  ! pull per epoch (F-Ftrue)/sig)
        ,SIM_EPdRA(MXEPOCH)  & 
        ,SIM_EPdDEC(MXEPOCH)  & 
        ,SIM_EPdMAG(MXEPOCH)  & 
        ,SIMSED_PARVAL(MXPAR_SIMSED)  & 
        ,PySEDMODEL_PARVAL(MXPAR_SIMSED)     &  ! Dec 10 2018
        ,LCLIB_PARVAL(MXPAR_LCLIB)  & 
        ,SIM_HOSTLIB_PARVAL(MXPAR_SIMSED)  &  ! HOSTLIB params
        ,SIM_MAGSMEAR_COH  & 
        ,SIM_SALT2gammaDM              &  ! usually this is gamma or gamma/2
        ,SIM_TEMPLATEMAG(MXFILT_ALL)   &  ! image-sub template, not LCLIB template
        ,SIM_LCWIDTH(MXFILT_ALL)       &  ! computed from SIM_EPFLUXCAL
        ,SIM_MODELGRID_TOBS(MXEP_MODELGRID)     &  ! for SNANA+SIM_MAGOBS table
        ,SIM_MODELGRID_MAGOBS(MXEP_MODELGRID,MXFILT_OBS/2) ! idem

    INTEGER  & 
         SIM_MODEL_INDEX       &  ! model index (MLCS,SALT2,NON1a ...)
        ,SIM_GENTYPE           &  ! generated SNTYPE; SIM_TYPE_INDEX
        ,SIM_TEMPLATE_INDEX    &  ! template index for NON1ASED, SIMSED, LCLIB ...
        ,SIM_SEARCHEFF_MASK    &  ! bits 1,2 => found by software,humans
        ,SIM_LIBID             &  ! LIBID for each event
        ,SIM_NGEN_LIBID        &  ! NGEN for this LIBID (usually 1)
        ,SIM_NOBS_UNDEFINED    &  ! NOBS where model is undefined
        ,SIM_NSUBSAMPLE_MARK   &  ! Number of marked sub-samples
        ,SIM_SUBSAMPLE_INDEX   &  ! sub-sample index
        ,SIM_REDSHIFT_FLAG     &  ! points to source of redshift
        ,SIMOPT_MWCOLORLAW     &  ! option for MW color law
        ,SIMOPT_MWEBV          &  ! option to modify MWEBV_SFD map
        ,NPAR_SIMSED  & 
        ,NPAR_PySEDMODEL       &  ! BYOSED, SNEMO ...
        ,NPAR_LCLIB  & 
        ,NPAR_SIM_HOSTLIB  & 
        ,SIM_EPFILTREST(MXEPOCH)  & 
        ,SIMLIB_MSKOPT         &  ! SIMLIB option mask (Dec 2015)
        ,SIM_BIASCOR_MASK      &  ! non-zero -> it's a biasCor sim
        ,NEP_SIM_MODELGRID  & 
        ,NEP_SIM_MAGOBS          ! number of sim epochs with MAGOBS < 99


    INTEGER*8  SIM_HOSTLIB_GALID
    REAL*8    DSIM_HOSTLIB_GALID  ! for table only

    CHARACTER  & 
         SIMNAME_MODEL*(MXCHAR_MODELNAME)      &  ! SALT2, mlcs2k2, NON1A, etc ...
        ,SIMNAME_TYPE*20       &  ! Ia, Ib, II, IIN, etc ...
        ,SIMNAME_SHAPEPAR*40   &  ! DELTA, x1, stretch ...
        ,SIMNAME_SHAPELAW*40   &  ! alpha ...
        ,SIMNAME_COLORPAR*40   &  ! AV, c ...
        ,SIMNAME_COLORLAW*40   &  ! RV,  beta ...
        ,SIMLIB_FILENAME*200   &  ! SIMLIB file used to generate sim
        ,SIMSED_KEYWORD(MXPAR_SIMSED)*80  &  ! full keyname in header
        ,SIMSED_PARNAME(MXPAR_SIMSED)*20  &  ! parname for fitres and ntuple
        ,PySEDMODEL_NAME*20               &  ! e.g., 'BYOSED' , 'SNEMO'
        ,PySEDMODEL_KEYWORD(MXPAR_SIMSED)*80  &  ! full keyname in header
        ,PySEDMODEL_PARNAME(MXPAR_SIMSED)*20  &  ! parname for fitres and ntuple
        ,LCLIB_KEYWORD(MXPAR_LCLIB)*80  &  ! full keyname in header
        ,LCLIB_PARNAME(MXPAR_LCLIB)*20  &  ! parname for fitres and ntuple
        ,SIM_HOSTLIB_KEYWORD(MXPAR_SIMSED)*60  &  ! full keyname
        ,SIM_HOSTLIB_PARNAME(MXPAR_SIMSED)*60  &  ! parname
        ,SIMNAME_SNRMON*40

    LOGICAL  & 
        LSIM_TRUE_SNIa    &  ! true SNIa (SALT3 or BAYES or MLCS ...)
       ,LSIM_TRUE_NONIa  ! true NONIa (SNCC, AGN, TDE ...)







  END MODULE SNSIMCOM

! =====================================================================
  MODULE PRIVCOM
    USE SNPAR
    IMPLICIT NONE

! Nov 2012: PRIVATE_VAR variables in data files
! Aug 2014: MXVAR_PRIVATE -> 28 (was 20)

    REAL*8, PARAMETER :: PRIVATE_NULL = -999.0909 

    INTEGER  & 
         NVAR_PRIVATE         &  ! Number of private  variables
        ,NCUT_PRIVATE

! for  the following header key
!    PRIVATE(BLA): 44.4
! we have
!   VARANME = 'BLA'
!   KEYWORD = 'PRIVATE(BLA)'
!   VALUE   = 44.4

    CHARACTER  & 
         PRIVATE_VARNAME(MXVAR_PRIVATE)*(MXCHAR_FILEWORD)  & 
        ,PRIVATE_KEYWORD(MXVAR_PRIVATE)*(MXCHAR_FILEWORD)

    REAL*8  & 
          PRIVATE_VALUE(MXVAR_PRIVATE)     &  ! value for each var
         ,PRIVATE_CUTWIN(2,MXVAR_PRIVATE) ! cut-window on each var

    INTEGER USE_PRIVATE_CUTWIN(MXVAR_PRIVATE)  ! 1=cutwin, -1=veto




  END MODULE PRIVCOM

! =====================================================================
  MODULE SNHOSTCOM
    USE SNPAR
    IMPLICIT NONE

! Host galaxy parameters.

! logical flags for output tables
    LOGICAL  & 
         EXIST_SNHOST_ANGSEP  & 
        ,EXIST_SNHOST_DDLR  & 
        ,EXIST_SNHOST_CONFUSION  & 
        ,EXIST_SNHOST_ZPHOT  & 
        ,EXIST_SNHOST_LOGMASS  & 
        ,EXIST_SNHOST_LOGSFR  & 
        ,EXIST_SNHOST_LOGsSFR  & 
        ,EXIST_SNHOST_COLOR  & 
        ,EXIST_SNHOST_MAGOBS  & 
        ,EXIST_SNHOST_SB

    INTEGER*8  SNHOST_OBJID(MXSNHOST)    ! int id
    REAL*8    DSNHOST_OBJID(MXSNHOST)   ! for tables only
    CHARACTER VARNAME_ZPHOT_Q(MXSNHOST,MXZPHOT_Q)*40

    REAL  & 
         SNHOST_ANGSEP(MXSNHOST)       &  ! SN-host sep, arcsec
        ,SNHOST_DDLR(MXSNHOST)         &  ! SNSEP/DLR
        ,SNHOST_CONFUSION              &  ! HC analog from Gupta 2016
        ,SNHOST_ZPHOT(MXSNHOST), SNHOST_ZPHOT_ERR(MXSNHOST)  & 
        ,SNHOST_ZPHOT_Q(MXSNHOST,MXZPHOT_Q)  & 
        ,SNHOST_ZPHOT_PERCENTILE(MXZPHOT_Q)  & 
        ,SNHOST_QZPHOT_MEAN(MXSNHOST), SNHOST_QZPHOT_STD(MXSNHOST)  & 
        ,SNHOST_ZSPEC(MXSNHOST), SNHOST_ZSPEC_ERR(MXSNHOST)  & 
        ,SNHOST_LOGMASS(MXSNHOST)  & 
        ,SNHOST_LOGMASS_ERR(MXSNHOST)  & 
        ,SNHOST_LOGSFR(MXSNHOST)  & 
        ,SNHOST_LOGSFR_ERR(MXSNHOST)  & 
        ,SNHOST_LOGsSFR(MXSNHOST)  & 
        ,SNHOST_LOGsSFR_ERR(MXSNHOST)  & 
        ,SNHOST_COLOR(MXSNHOST)  & 
        ,SNHOST_COLOR_ERR(MXSNHOST)  & 
        ,SNHOST_SBFLUXCAL(MXFILT_ALL)   &  ! surface brightness FLUXCAL /asec^2
        ,SNHOST_MAGOBS(MXFILT_ALL,(MXSNHOST))        &  ! observer-frame mags
        ,SNHOST_MAGOBS_ERR(MXFILT_ALL,(MXSNHOST))   ! error on above

    REAL*8  & 
          SNHOST8_RA(MXSNHOST)  & 
         ,SNHOST8_DEC(MXSNHOST)

! quantites which do not depend on which host
    INTEGER*4  & 
          SNHOST_NMATCH         &  ! number of host matches, e.g., d_DLR<4
         ,SNHOST_NMATCH2        &  ! number of host matches, e.g., d_DLR<7
         ,SNHOST_NZPHOT_Q       &  ! May 2022: number of zphot quantiles
         ,SNHOST_FLAG(MXSNHOST)    ! May 21 2021: indicate problems with host

    REAL  & 
         SNHOST_SBFLUXCAL_ERR(MXFILT_ALL)  & 
        ,SNHOST_SBMAG(MXFILT_ALL)        ! surface brightness mag/asec^2



  END MODULE SNHOSTCOM

! =====================================================================
  MODULE SIMLIBCOM
    USE SNPAR
    IMPLICIT NONE

! Feb 2017: common block for writing SIMLIB file to SIMLIB_OUT

    INTEGER, PARAMETER :: MXLIST_SKY = 20
    INTEGER  NLIST_SKY
    REAL*8   SKYLAM_LIST(MXLIST_SKY), SKYMAG_LIST(MXLIST_SKY)
    REAL*8   PSF_FWHM_GUESS, PIXSIZE_GUESS, ADD_SKYSIG_PIX

    LOGICAL ISGROUND, FOUND_TEMPLATE_INFO ! T=> ground based survey
    REAL  & 
          SIMLIB_ZPERR(MXFILT_ALL)  & 
         ,SIMLIB_TEMPLATE_SKYSIG(MXFILT_ALL)  & 
         ,SIMLIB_TEMPLATE_ZPT(MXFILT_ALL)





! --- hard-wire ground based array of SKYMAG vs. lambda
! Apr 6 2020:
!  Add J,H,K from Table 3 of arxiv:0809.4988, and increase NLIST->10.
!  This CALAR ALTO paper doesn't have a table of mean filter wavelengths,
!  so we take <SKYLAM> from CSP JHK bands.
! Jul 7 2025:
!  L' from https://www.gemini.edu/instrumentation/niri/capability

    INTEGER, PARAMETER :: NLIST_SKY_GROUND = 11
    REAL*8 SKYLAM_GROUND_LIST(NLIST_SKY_GROUND)
    REAL*8 SKYMAG_GROUND_LIST(NLIST_SKY_GROUND)

! ugrizY mag/asec^2 are from LSST DEEP-drill cadence;
! v is made up to extend bluer.
    DATA SKYLAM_GROUND_LIST /  & 
           2700.0, 3714.0, 4790.0, 6220.0,   &  ! v,u,g,r
           7544.0, 8679.0, 10095,            &  ! i,z,Y
           12500.0, 18900.0, 21500.0,        &  ! J,H,K (Apr 2020)
           37700.0 /                        ! L' (Jul 2025)
    DATA SKYMAG_GROUND_LIST /  & 
           23.8, 22.7,  21.0,  20.4,     &  ! v,u,g,r
           19.4, 18.1,  17.9,            &  ! i,z,Y
           16.0, 14.0,  12.4,            &  ! J,H,K (Apr 2020)
           3.5 /                        ! L' (Jul 2025)

! --- hard-write space based array pf SKYMAG vs. lambda
    INTEGER, PARAMETER :: NLIST_SKY_SPACE = 13
    REAL*8 SKYLAM_SPACE_LIST(NLIST_SKY_SPACE)
    REAL*8 SKYMAG_SPACE_LIST(NLIST_SKY_SPACE)

! from WFIRST sims
    DATA SKYLAM_SPACE_LIST /  & 
              1000.0,  4330.0,  6258.0,  8052.0,           &  ! x,B,R,I
              8745.0, 10653.0, 12976.0, 15848.0, 20000.,   &  ! Z,Y,J,H,FUDGE
              27700.0, 35600.0, 44400.0, 50000.0        / ! NIRCAM LW
    DATA SKYMAG_SPACE_LIST /  & 
              24.0,    23.8,     23.6,   23.5,  & 
              23.5,    23.7,     23.8,   23.9, 24.2,  & 
              23.1,    23.3,     22.3,   22.0	        /

  END MODULE SIMLIBCOM

! =====================================================================
  MODULE TRUECHI2COM
    USE SNPAR
    IMPLICIT NONE

! for fakes and sim, compute sum[F-Ftrue]^2/sigma^2 per band
! The zeroth element is all bands together.

    REAL  & 
         CHI2_TRUEFLUX(0:MXFILT_ALL)  & 
        ,PROB_TRUEFLUX(0:MXFILT_ALL)

    INTEGER NDOF_TRUEFLUX(0:MXFILT_ALL)


  END MODULE TRUECHI2COM

! =====================================================================
  MODULE SNCUTS
    USE SNPAR
    IMPLICIT NONE

! Define cut-window for each variables
! and cut-mask for each epoch.
! User must fill cutwin_XXX arrays in main routine.


    INTEGER, PARAMETER ::           & 
         CUTBIT_CID            = 1  & 
        ,CUTBIT_SNTYPE         = 2  & 
        ,CUTBIT_REDSHIFT       = 3  & 
        ,CUTBIT_REDSHIFT_ERR   = 4  & 
        ,CUTBIT_RA             = 5  & 
        ,CUTBIT_DEC            = 6  & 
        ,CUTBIT_HOSTSEP        = 7  & 
        ,CUTBIT_TrestMIN       = 8  & 
        ,CUTBIT_TrestMAX       = 9  & 
        ,CUTBIT_TrestRange     = 10   &  ! Dec 2017 (TrestMax-TrestMin)
        ,CUTBIT_Trest2         = 11  & 
        ,CUTBIT_Tgapmax        = 12  & 
        ,CUTBIT_T0gapMAX       = 13  & 
        ,CUTBIT_TobsMIN        = 14  & 
        ,CUTBIT_TobsMAX        = 15  & 
        ,CUTBIT_PEAKMJD        = 16  & 
        ,CUTBIT_NOBS_PREDETECT = 17   &  ! Mar 2025
        ,CUTBIT_NSEASON_ACTIVE = 18   &  ! May 2019
        ,CUTBIT_NEPOCH         = 19   &  ! total # of measurements
        ,CUTBIT_SEARCH         = 20   &  ! found by search (SIM only)
        ,CUTBIT_NFILT_SNRMAX   = 21  & 
        ,CUTBIT_NFILT_SNRMAX2  = 22  & 
        ,CUTBIT_NFILT_TRESTMIN = 23  & 
        ,CUTBIT_NFILT_TRESTMAX = 24  & 
        ,CUTBIT_NFILT_TREST2   = 25  & 
        ,CUTBIT_SNRMAX         = 26  & 
        ,CUTBIT_SNRMAX2        = 27  & 
        ,CUTBIT_SNRSUM         = 28   &  ! Jul 2024
        ,CUTBIT_NFIELD         = 29  & 
        ,CUTBIT_MWEBV          = 30  & 
        ,CUTBIT_REQEP          = 31   &  ! required epochs (9/17/2017)
        ,CUTBIT_PRIVATE        = 32   &  ! all private-var cuts
        ,CUTBIT_SIMVAR         = 33   &  !
        ,CUTBIT_MASK_NOBS_TREST= 34   &  ! Nov 24 2025 require Nobs in multiple Trest ranges
        ,CUTBIT_OFFSET_SBFLUX  = 35  & 
        ,CUTBIT_OFFSET_SNRMAX  = 35 + MXFILT_SNRMAX  & 
        ,CUTBIT_MJD_MARKER     = CUTBIT_OFFSET_SNRMAX + MXFILT_SNRMAX  & 
        ,NCUTBIT_SNLC          = CUTBIT_MJD_MARKER  & 
! SN-dependent cuts above
!  MJD-dependent cuts below
        ,CUTBIT_PSF          = CUTBIT_MJD_MARKER + 1   &  ! PSF cut
        ,CUTBIT_ZP           = CUTBIT_MJD_MARKER + 2   &  ! ZP  cut (Feb 2018)
        ,CUTBIT_ZPERR        = CUTBIT_MJD_MARKER + 3   &  ! ZPERR (Feb 2020)
        ,CUTBIT_PHOTPROB     = CUTBIT_MJD_MARKER + 4  & 
        ,CUTBIT_Trest        = CUTBIT_MJD_MARKER + 5  & 
        ,CUTBIT_Tobs         = CUTBIT_MJD_MARKER + 6  & 
        ,CUTBIT_MJD          = CUTBIT_MJD_MARKER + 7  & 
        ,CUTBIT_TREST_TRUEFLUX2= CUTBIT_MJD_MARKER + 8  & 
        ,CUTBIT_ERRTEST      = CUTBIT_MJD_MARKER + 9  & 
        ,CUTBIT_SIM_PULL     = CUTBIT_MJD_MARKER + 10  &  ! Mar 2021
        ,CUTBIT_OFFSET       = CUTBIT_SIM_PULL  & 
        ,NCUTBIT             = CUTBIT_OFFSET 
          

! Define cut-masks with 64-BIT integers

    INTEGER*8  & 
         CUTMASK8_SN              &  ! cutmask for each SN
        ,CUTMASK8_MJD(MXEPOCH)    &  ! cutmask vs. MJD, isn
        ,CUTMASK8_SN_ALL          &  ! all bits for SN cuts
        ,CUTMASK8_MJD_ALL        ! all bits for MJD cuts

    BYTE MJDMASK(MXEPOCH)  ! logical for each epoch,SN

    LOGICAL  & 
         LSNCUTS        &  ! T=> passes cuts 1 to CUTBIT_MJD_MARKER
        ,PASS_PRIVCUTS  &  ! T=> pass cuts on private var in data header
        ,PASS_SIMCUTS  ! T=> pass cuts on SIMVAR

    INTEGER  & 
         NCCID_LIST       &  ! size of SNCCID_LIST
        ,NCID_LIST        &  ! size of SNCID_LIST
        ,NCCID_IGNORE     &  ! size of SNCCID_IGNORE
        ,NCID_IGNORE      &  ! size of NCID_IGNORE
        ,CUTFLAG_SNANA    &  ! bits 0,1 -> LSNCUTS,LSNFITOK=T (for ntuple)
        ,ERRFLAG_FIT     ! error flag from fit (0=OK)

! stuff cut vars
    CHARACTER CUTVAR_NAME(NCUTBIT)*28

! define namelist cuts

    INTEGER  & 
         SNTYPE_LIST(MXLISTNML)     &  ! I: user list of types to select
        ,SNTYPE_IGNORE(MXLISTNML)   &  ! I: types to ignore
        ,SNCID_LIST(MXLISTNML)      &  ! I: user-list of integer CIDs
        ,MXEVT_PROCESS              &  ! I: quit after processing this many
        ,MXEVT_CUTS                 &  ! I: quit after Nevt passing cuts
        ,DETNUM_LIST(MXLISTNML)     &  ! I: user-list of Detector/CCD NUMBERS (Dec 2017)
        ,SIM_TEMPLATE_INDEX_LIST(MXLISTNML)  & 
        ,SNCID_IGNORE(MXLISTNML)    &  ! I: list of CIDs to ignore
        ,PHOTFLAG_MSKREJ(5)         &  ! I: PHOTFLAG mask-bits to reject, 1,2 => logical OR,AND
        ,PHOTFLAG_BITLIST_REJECT(MXBIT_PHOTFLAG)   &  ! I: bits to reject
        ,IDFIELD_LIST(MXLISTNML)    &  !
        ,NIDFIELD_LIST  & 
        ,NSNTYPE_LIST            &  ! size of SNTYPE_LIST
        ,NDETNUM_LIST            &  ! size of DETNUM list
        ,NSIM_TEMPLATE_INDEX_LIST  & 
        ,NSNTYPE_IGNORE  & 
        ,NFILT_SNRMAX            &  ! number of CUTWIN_SNRMAX_FILT cuts
        ,IFILT_SNRMAX(MXFILT_SNRMAX)  &  ! store filt index for SNRMAX cuts
        ,NFILT_HOST_SBFLUX       &  ! Nfilt to cut on HOST_SBFLUX
        ,IFILT_HOST_SBFLUX(MXFILT_SNRMAX)

    LOGICAL DOALL_SNTEL, DOALL_SNFIELD

    character  & 
          SNFIELD_LIST(MXLISTNML)*60   &  ! I:  list of fields to use (or 'ALL')
         ,SNCCID_LIST(MXLISTNML)*(MXCHAR_CCID)    &  ! I: list of CIDs to process
         ,SNCCID_IGNORE(MXLISTNML)*(MXCHAR_CCID)  &  ! I: list of CIDs to ignore
         ,SNCID_IGNORE_FILE*(MXCHAR_FILENAME)     &  ! I: file with CIDs to ignore
         ,SNCCID_IGNORE_ALL(MXIGNORE_LIST)*(MXCHAR_CCID)

!  -----------  SN-dependent cut-window  ! I: --------------------------------
    REAL  & 
         snlc_cutvar(NCUTBIT_SNLC)  & 
        ,snep_cutvar(CUTBIT_MJD_MARKER:NCUTBIT,MXEPOCH)  & 
        ,cutwin_var(2,NCUTBIT)  & 
        ,cutwin_cid(2)            &  ! candidate id
        ,cutwin_sntype(2)         &  !
        ,cutwin_redshift(2)       &  ! I: redshift cut
        ,cutwin_redshift_err(2)   &  ! I: cut on redshift uncertainty
        ,cutwin_ra(2)             &  ! I: cut on RA
        ,cutwin_dec(2)            &  ! I: cut on DEC
        ,cutwin_hostsep(2)        &  ! I: cut on host-SN sep, arcsec
        ,cutwin_sbflux_filt(2,MXFILT_SNRMAX)  & 
        ,cutwin_Nepoch(2)         &  ! I: cut on Nepoch
        ,cutwin_snrmax_filt(2,MXFILT_SNRMAX)    &  ! filled from SNCUT_SNRMAX
        ,cutwin_snrmax(2)         &  ! I: global SNRMAX cut
        ,cutwin_snrmax2(2)        &  ! I: 2nd global SNRMAX cut (different band)
        ,cutwin_snrsum(2)         &  ! I: cut on sqrt[  sum(SNR_i^2) ]
        ,cutwin_nfield(2)         &  ! I: number of fields (usually 1)
        ,cutwin_mwebv(2)          &  ! I: Galactic extinction
        ,cutwin_nseason_active(2)       &  ! I: number of seasons with detection
        ,cutwin_cutflag_reqep(2)        &  ! for internal use only
        ,cutwin_cutflag_private(2)      &  ! for internal use only
        ,cutwin_cutflag_simvar(2)       &  ! for internal use only
        ,cutwin_mask_nobs_trest(2)      &  ! for internal use only
        ,cutwin_nfilt_snrmax(2)     &  ! I: Nfilt passing global SNRMAX cut
        ,cutwin_nfilt_snrmax2(2)    &  ! I: Nfilt passing 2nd-best SNRMAX
        ,cutwin_nfilt_trestmin(2)   &  ! I: Nfilt passing Trestmin cut
        ,cutwin_nfilt_trestmax(2)   &  ! I: Nfilt passing Trestmax cut
        ,cutwin_nfilt_trest2(2)     &  ! I: Nfilt passinng Trest2 cut
        ,cutwin_nband_thresh(2)     &  ! I: number of bands above $band_THRESH
        ,cutwin_Trestmin(2)         &  ! I: cut on earliest epoch, rest frame days
        ,cutwin_Trestmax(2)         &  ! I: cut on  latest epoch, rest frame day
        ,cutwin_TrestRange(2)       &  ! I: cut on Trestmax - Trestmin
        ,cutwin_Trest2(2)           &  ! I: cut on CUTWIN_NFILT_TREST2
        ,cutwin_Tgapmax(2)          &  ! I: max Trest-gap within cutwin_TREST(2)
        ,cutwin_T0gapmax(2)         &  ! I: max Trest-gap near peak
        ,cutwin_Tobsmin(2)          &  ! I: cut on min Tobs
        ,cutwin_Tobsmax(2)          &  ! I: cut on max Tobs
        ,cutwin_peakmjd(2)          &  ! I: cut on search peakMJD
        ,cutwin_tobs_predetect(2)   &  ! I: pre-detection Tobs window  e.g. -10, 0 w.r.t. MJD_DETECT_FIRST
        ,cutwin_nobs_predetect(2)  ! I: Nobs required in pre-detect time window


!  Below are Epoch-dependent cuts
!  -----------  epoch-dependent cut windows  ! I: --------------------------------
    REAL  & 
         cutwin_psf(2)            &  ! I: PSF cut, FWHM, ARCSEC
        ,cutwin_zp(2)             &  ! I: ZEROPT cut, ADU or NPE
        ,cutwin_zpADU(2)          &  ! I: ZEROPT cut, ADU
        ,cutwin_zpNPE(2)          &  ! I: ZEROPT cut, NPE
        ,cutwin_ZPERR(2)          &  ! I: ZPERR cut
        ,cutwin_photprob(2)       &  ! I: cut on PHOTPROP (e.g., real/bogus)
        ,cutwin_errtest(2)        &  ! I: cut on ERRCALC/ERRTRUE
        ,cutwin_sim_pull(2)         &  ! I: cut on (F-Ftrue)/sigmaF
        ,cutwin_Trest_trueflux(2)   &  ! I: Trest range requiring valid true flux
        ,cutwin_Trest_trueflux2(2)  &  ! I: cut on logical flag
        ,cutwin_Trest(2)          &  ! I: cut on all epochs, rest frame days
        ,cutwin_Tobs(2)           &  ! I: cut on all epochs, obs-frame days
        ,cutwin_MJD(2)            &  ! I: cut on MJD range
        ,cutwin_mjd_exclude(2)    &  ! I: MJD window to exclude
        ,cutwin_SEARCHEFF_MASK(2)        &  ! I: cut on SIM_SEARCHEFF_MASK
        ,cutwin_snrmin_filt(2,MXFILT_OBS)  &  ! filled from EPCUT_SNRMIN
        ,cutwin_restlam(2)                 &  ! I: cut on <LAMREST>, no cutBit
        ,cutwin_lamrest(2)                 &  ! I: same as above
        ,cutwin_lamobs(2)                  &  ! I: cut on <LAMOBS>, no cutBit
        ,cutwin_snr_nodetect(2)    ! I: SNR cut for non-detections


! define character strings to specify cuts and filters
! 'SNCUT' specifies cut on each SN
! 'EPCUT' specifies selection cut on each epoch

  ! I:# --------------------------------
    CHARACTER*(MXCHAR_CUTNAME)  & 
        SNCUT_SNRMAX   &  ! I: max SNR required in each passband
                    ! I: example: 'u 10.  g 5.0  r 5.0  i 5.0  z -10.'
! 
       ,EPCUT_SNRMIN   &  ! I: min SNR accepted for each epoch/filter
                    ! I: example: 'u 20000.  g -4.  r -4.  i -4.  z 20000.'
! 
       ,EPCUT_SKYSIG   &  ! I: max SKY noise accepted each epoch/filter
                    ! I: example: 'u 20.  g 50.  r 80.  i 120.  z 200.'
! 
       ,EPCUT_PSFSIG   &  ! I: max PSF (arcsec) accepted each epoch/filter
                    ! I: example: 'u 1.8  g 1.5  r 1.6  i 1.8  z 2.0'
! 
       ,SNCUT_HOST_SBFLUX  &  ! I: require min surface brightness,
                        ! I: example 'r 1000.' -> SBFLUX_r > 1000
! 
       ,SNCUT_NOBS_TREST   &  ! I: multiple NOBS cuts in Trest ranges
                        ! I: e.g, '1 -10 0  2 0 30' -> >= 1 obs in -10 to 0 days & >=2 obs in 0 to 30 days
    ! I:# --------------------------------
! ----------
! systematic tests for calibration:
! 'U 01 B -0.01' => shift U & B mags of primary
! 
        ,MAGOBS_SHIFT_PRIMARY   &  ! I: shift primary mags  (for syst test)
        ,MAGOBS_SHIFT_ZP        &  ! I: shift zero points (e.g.,AB off for SDSS)
        ,MAGREST_SHIFT_PRIMARY  &  ! I: idem for rest-frame mags
        ,MAGREST_SHIFT_ZP       &  ! I: idem for rest-frame ZP
        ,FILTER_LAMSHIFT        &  ! I: e.g., 'r -2.4  i 6.2'  ! in Angstroms
! -------
! Fluxcal fudges for systematic tests (fudge photometry offsets and errors)
! Note that error is added in quadrature; ERROR<0 is subtracted in quadrature.
        ,FUDGE_FLUXCAL_OFFSET   &  ! I: shift all FLUXCAL
        ,FUDGE_FLUXCAL_ERROR    &  ! I: fudge net FLUXCAL error in each band
        ,FUDGE_FLUXCAL_ERRPIX   &  ! I: per-pixel error --> FLUXCAL error propto PSF
        ,FUDGE_MAG_ERROR        &  ! I: add stat error per band
        ,FUDGE_MAG_COVERR       &  ! I: add covariant error per band
        ,FUDGE_FLUXERR_SCALE    &  ! I: scale error in each band
        ,SIM_FUDGE_MAG_ERROR   ! I: add stat error per band, sim only

    REAL  & 
         MAGOBS_SHIFT_PRIMARY_FILT(MXFILT_ALL)  &  ! shift mag of primary ref
        ,MAGOBS_SHIFT_ZP_USER(MXFILT_ALL)       &  ! user ZP shift
        ,MAGOBS_SHIFT_ZP_FILT(MXFILT_ALL)       &  ! user ZP shift + system ZP
        ,FLUXSCALE_ZP_FILT(MXFILT_ALL)          &  ! corresponding flux scale
        ,MAGREST_SHIFT_PRIMARY_FILT(MXFILT_ALL)  &  ! shift mag of primary ref
        ,MAGREST_SHIFT_ZP_USER(MXFILT_ALL)       &  ! shift zero points
        ,MAGREST_SHIFT_ZP_FILT(MXFILT_ALL)       &  ! shift zero points
        ,FILTER_LAMSHIFT_FILT(MXFILT_ALL)        &  ! shift filter trans
! 
        ,MAGOBS_SHIFT_PRIMARY_PARAMS(3)        &  ! poly-fun of lambda;
        ,MAGOBS_SHIFT_ZP_PARAMS(3)             &  ! A0 + A1*LAM + A2*LAM^2
        ,FUDGE_FLUXCAL_OFFSET_FILT(MXFILT_ALL)  & 
        ,FUDGE_FLUXCAL_ERROR_FILT(MXFILT_ALL)  & 
        ,FUDGE_FLUXCAL_ERRPIX_FILT(MXFILT_ALL)  & 
        ,FUDGE_FLUXERR_SCALE_FILT(MXFILT_ALL)  & 
        ,FUDGE_MAG_ERROR_FILT(MXFILT_ALL)  & 
        ,FUDGE_MAG_COVERR_FILT(MXFILT_ALL)  & 
        ,SIM_FUDGE_MAG_ERROR_FILT(MXFILT_ALL)  & 
        ,MWEBV_SCALE         &  ! I: scale MW extinc for syst test
        ,MWEBV_SHIFT         &  ! I: shift MW extinc
        ,MWEBV_FORCE        ! I: force specific MWEBV value (Feb 2020)

    REAL  & 
         RV_MWCOLORLAW   ! I: RV for Galactic extinction

    INTEGER  & 
         OPT_MWCOLORLAW     &  ! I: MW color law opt (89, 94, 99 ...)
        ,OPT_MWEBV         ! I: option to modify SFD maps
    REAL  & 
         PARLIST_MWCOLORLAW(10) ! I: params for MWCOLORLAW calc

    LOGICAL  & 
         USE_MWCOR          ! I:  T=> correct data flux for MW extinc;
                          ! I:  F=> leave data, correct fit model for MW.
    REAL  & 
         REDSHIFT_FINAL_SHIFT      &  ! I: artificial z-shift
        ,HOSTGAL_PHOTOZ_SHIFT      &  ! I: obsolete
        ,HOSTGAL_SPECZ_SHIFT       &  ! I: obsolete
        ,HOSTGAL_ZPHOT_SHIFT       &  ! I: apply only when zFINAL = HOST_ZPHOT
        ,HOSTGAL_ZSPEC_SHIFT       &  ! I: apply only when zFINAL = HOST_ZPHOT
        ,FLUXERRCALC_ZPTERR



    EQUIVALENCE  & 
          ( cutwin_var(1,cutbit_cid),       cutwin_cid )  & 
         ,( cutwin_var(1,cutbit_sntype),    cutwin_sntype )  & 
         ,( cutwin_var(1,cutbit_redshift),  cutwin_redshift )  & 
         ,( cutwin_var(1,cutbit_redshift_err),cutwin_redshift_err )  & 
         ,( cutwin_var(1,cutbit_ra),        cutwin_ra )  & 
         ,( cutwin_var(1,cutbit_dec),       cutwin_dec )  & 
         ,( cutwin_var(1,cutbit_hostsep),   cutwin_hostsep )  & 
         ,( cutwin_var(1,cutbit_Nepoch),    cutwin_Nepoch )  & 
         ,( cutwin_var(1,cutbit_psf),       cutwin_psf )  & 
         ,( cutwin_var(1,cutbit_zp),        cutwin_zp  )  & 
         ,( cutwin_var(1,cutbit_zperr),     cutwin_zperr  )  & 
         ,( cutwin_var(1,cutbit_photprob),  cutwin_photprob  )  & 
         ,( cutwin_var(1,cutbit_errtest),   cutwin_errtest   )  & 
         ,( cutwin_var(1,cutbit_sim_pull),  cutwin_sim_pull   )  & 
         ,( cutwin_var(1,cutbit_trest_trueflux2),cutwin_trest_trueflux2)  & 
         ,( cutwin_var(1,cutbit_Nfilt_snrmax), cutwin_nfilt_snrmax )  & 
         ,( cutwin_var(1,cutbit_Nfilt_snrmax2),cutwin_nfilt_snrmax2 )  & 
         ,( cutwin_var(1,cutbit_Nfilt_trestmin),cutwin_nfilt_trestmin)  & 
         ,( cutwin_var(1,cutbit_Nfilt_trestmax),cutwin_nfilt_trestmax)  & 
         ,( cutwin_var(1,cutbit_Nfilt_trest2),  cutwin_nfilt_trest2)  & 
         ,( cutwin_var(1,cutbit_Trestmin),  cutwin_Trestmin )  & 
         ,( cutwin_var(1,cutbit_Trestmax),  cutwin_Trestmax )  & 
         ,( cutwin_var(1,cutbit_TrestRange),cutwin_TrestRange )  & 
         ,( cutwin_var(1,cutbit_Trest2),    cutwin_Trest2 )  & 
         ,( cutwin_var(1,cutbit_Tgapmax),   cutwin_Tgapmax )  & 
         ,( cutwin_var(1,cutbit_T0gapmax),  cutwin_T0gapmax )  & 
         ,( cutwin_var(1,cutbit_Tobsmin),   cutwin_Tobsmin )  & 
         ,( cutwin_var(1,cutbit_Tobsmax),   cutwin_Tobsmax )  & 
         ,( cutwin_var(1,cutbit_Trest),     cutwin_Trest )  & 
         ,( cutwin_var(1,cutbit_Trest_trueflux2),cutwin_Trest_trueflux2)  & 
         ,( cutwin_var(1,cutbit_Tobs),      cutwin_Tobs )  & 
         ,( cutwin_var(1,cutbit_MJD),       cutwin_MJD   )  & 
         ,( cutwin_var(1,cutbit_peakmjd),   cutwin_peakmjd )  & 
         ,( cutwin_var(1,cutbit_nobs_predetect),cutwin_nobs_predetect)  & 
         ,( cutwin_var(1,cutbit_nseason_active),cutwin_nseason_active)  & 
         ,( cutwin_var(1,cutbit_search),        cutwin_searcheff_mask)  & 
         ,( cutwin_var(1,cutbit_snrmax),        cutwin_snrmax  )  & 
         ,( cutwin_var(1,cutbit_snrmax2),       cutwin_snrmax2  )  & 
         ,( cutwin_var(1,cutbit_snrsum),        cutwin_snrsum  )  & 
         ,( cutwin_var(1,cutbit_nfield),        cutwin_nfield  )  & 
         ,( cutwin_var(1,cutbit_mwebv),         cutwin_mwebv  )  & 
         ,( cutwin_var(1,cutbit_reqep),         cutwin_cutflag_reqep)  & 
         ,( cutwin_var(1,cutbit_private),       cutwin_cutflag_private)  & 
         ,( cutwin_var(1,cutbit_simvar),         cutwin_cutflag_simvar)  & 
         ,( cutwin_var(1,cutbit_mask_nobs_trest),cutwin_mask_nobs_trest)  & 
         ,( cutwin_var(1,cutbit_offset_sbflux+1),cutwin_sbflux_filt)  & 
         ,( cutwin_var(1,cutbit_offset_snrmax+1),cutwin_snrmax_filt)  & 
         ,( cutwin_restlam(1), cutwin_lamrest(1) )


  END MODULE SNCUTS

! =====================================================================
  MODULE SNLCINP_NML
    USE SNPAR
    USE SNCUTS
    IMPLICIT NONE

! define general input to program:
! 
    CHARACTER*(MXCHAR_FILENAME)   &  ! user input files
         NMLFILE           &  ! name of input namelist file
        ,ROOTFILE_OUT      &  ! I: output filename for root for tables
        ,TEXTFILE_PREFIX   &  ! I: prefix for text file tables
        ,MARZFILE_OUT      &  ! I: output FITS file for input to MARZ (spectra)
        ,SNTABLE_LIST      &  ! I: list of SNTABLEs to create; e.g., 'SNANA FITRES'
        ,SNTABLE_APPEND_VARNAME   &  ! I: add this varname to output tables
        ,SNTABLE_FILTER_REMAP  &  ! I: remap filters in tables
        ,CALIB_FILE            &  ! I: fits-formatted KCOR/calib tables
        ,KCOR_FILE             &  ! I: legacy key for CALIB_FILE
        ,OVERRIDE_RESTLAM_BOUNDARY*80  &  ! I: for kcors; e.g., '5100(gV) 6200(rV)'
        ,EPOCH_IGNORE_FILE     &  ! I: user file with epochs to ignore
        ,OUT_EPOCH_IGNORE_FILE  &  ! I: write this file with discarded epochs
        ,SNMJD_LIST_FILE       &  ! I: list of "CID MJD" to process
        ,SNMJD_OUT_FILE        &  ! I: MJD_LIST out file
        ,MNFIT_PKMJD_LOGFILE   &  ! I: separate log for PKMJD fits
        ,USERTAGS_FILE         &  ! I: optional int user tag per SN
        ,VPEC_FILE             &  ! obsolete? pec. velocity corrections
        ,SNCID_LIST_FILE       &  ! I: list of CID's to process
        ,SIMLIB_OUT            &  ! legacy name
        ,SIMLIB_OUTFILE        &  ! I: write simlib entry for each event
        ,SIMLIB_ZPERR_LIST     &  ! I: e.g., 'abc .01 def .02 ghi .005'
        ,NONLINEARITY_FILE     &  ! I: non-linearity map
! 
        ,FUDGE_HOSTNOISE_FILE   &  ! legacy: inflate error vs. hostSB (Mar 2015)
        ,FLUXERRMODEL_FILE      &  ! I: DATA ONLY: err-fudge maps
        ,SIM_FLUXERRMODEL_FILE  &  ! I: idem, SIM only
        ,MAGCOR_FILE            &  ! I: DATA ONLY: mag-cor for each CID-MJD-band
        ,SIM_MAGCOR_FILE       ! I: idem, SIM only

    CHARACTER*(MXFILE_LIST*MXCHAR_FILENAME)  & 
          HEADER_OVERRIDE_FILE        &  ! I: comma-sep list of files with CID VAR
         ,SIM_HEADER_OVERRIDE_FILE   ! I: same, but for sims

    CHARACTER   &  ! versions
         VERSION_PHOTOMETRY(MXVERS)*(MXCHAR_VERSION)    &  ! I: SN versions to read
        ,VERSION_PHOTOMETRY_WILDCARD*(MXCHAR_FILENAME)  &  ! I: get list from glob
        ,VERSION_REFORMAT_FITS*(MXCHAR_VERSION)    &  ! I: reformat TEXT -> FITS
        ,VERSION_REFORMAT_TEXT*(MXCHAR_VERSION)   ! I: reformat FITS -> TEXT

    CHARACTER  &  ! paths
         PRIVATE_DATA_PATH*(MXCHAR_PATH)     &  ! I: private data subdir
        ,FILTER_UPDATE_PATH*(MXCHAR_PATH)   ! I: SN-dependent filter response

    CHARACTER   &  ! misc
         REFORMAT_KEYS*(MXCHAR_FILENAME)   &  !   global reformat info
        ,NONSURVEY_FILTERS*(MXFILT_ALL)    &  ! I: non-survey filters to add
        ,SNRMAX_FILTERS*(MXFILT_ALL)       &  ! I: list of filters for SNRMAX cuts
        ,FILTER_REPLACE*(MXFILT_ALL)       &  ! I: e.g., 'UGRIZ -> ugriz'
        ,FILTLIST_LAMSHIFT*(MXFILT_ALL)    &  ! I: list of lam-shifted filters
        ,PRIVATE_CUTWIN_STRING(MXCUT_PRIVATE)*(MXCHAR_CUTNAME)  &  ! I: cut on privat variables
        ,PRIVATE_VARNAME_READLIST*200    &  ! I: list of private vars to read (default=ALL)
        ,SIMVAR_CUTWIN_STRING*(MXCHAR_CUTNAME)  &  ! I: cuts on SIM_XXX
        ,EARLYLC_STRING*(MXCHAR_CUTNAME)        &  ! I: see manual
        ,REQUIRE_EPOCHS_STRING*100   &  ! I: e.g., 'riz 10 7 20' uses CUTWIN_SNRMAX
        ,DUMP_STRING*100            ! 'funName CID-list'

    INTEGER  & 
         NFIT_ITERATION        &  ! I: number of fit iterations
        ,MINUIT_PRINT_LEVEL    &  ! -1 -> none
        ,INTERP_OPT            &  ! I:  interp option (see INTERP_XXX params)
        ,NLINE_ARGS  & 
        ,OPT_YAML              &  ! I: 1=> write YAML out even if not batch job
        ,OPT_REFORMAT_TEXT     &  ! I: 1=> re-write SNDATA files in text format; see manual options
        ,OPT_REFORMAT_FITS     &  ! I: 1=> re-write in FITSformat (data only); see manual options
        ,OPT_REFORMAT_SALT2    &  ! I: 1=> re-write SNDATA files in SALT2 format
        ,OPT_REFORMAT_SPECTRA  &  ! I: 1=> 3-col text: lam, Flam, Flamerr (ignore photometry)
        ,OPTSIM_LCWIDTH        &  ! I: 1=> option to compute LC width
        ,DEBUG_FLAG            &  ! I: for internal debug
        ,JOBSPLIT(2)           &  ! I: for batch; process (1)-range of (2)=TOTAL
        ,JOBSPLIT_EXTERNAL(2)  &  ! passed by submit_batch for text format
        ,MXLC_FIT              &  ! I: stop after this many fits passing all cuts
        ,PHOTFLAG_DETECT       &  ! I: used to count NEPOCH_DETECT and TLIVE_DETECT
        ,PHOTFLAG_TRIGGER      &  ! I: determine MJD(trigger) for survey
        ,FLUXERRMODEL_OPTMASK  &  ! I: see manual
        ,OPT_SIMLIB_OUT        &  ! I: bit mask of SIMLIB_OUT options
        ,SIMLIB_OUT_TMINFIX    &  ! I: choose PEAKMJD so that min(MJD-PEAKMJD)=TMINFIX
        ,REQUIRE_DOCANA  & 
        ,OPT_SNCID_LIST       &  ! I: 1=force all and ignore cuts;
                           ! I: 2=set INIVAL=FITPAR
                           ! I: 4=use each FITPAR and ERROR as prior
        ,OPT_VPEC_COR        ! I: 1=apply vpec cor (default)

    LOGICAL  & 
         LSIM_SEARCH_SPEC    &  ! I: T => require simulated SPEC-tag
        ,LSIM_SEARCH_ZHOST   &  ! I: T => require simulated zHOST
        ,LTEST_INTERP        &  ! I: T => calls TEST_KCOR
        ,LTEST_KCOR          &  ! I: T => calls TEST_INTERP
        ,LTEST_MAG           &  ! I: T => test GET_MAGLC  and GET_MWXT
        ,LTEST_U3BAND        &  ! I: T => require exactly 3 bands that include U
        ,USE_LINE_ARGS(MXLINE_ARGS)  & 
        ,USE_SNHOST_ZPHOT     &  ! I: T=> replace SNLC_REDSHIFT -> SNHOST_ZPHOT
        ,USE_HOSTGAL_PHOTOZ   &  ! I: idem, but matches keyName in data files
! 
        ,ABORT_ON_NOEPOCHS    &  ! I: T=> abort if there are no epochs to fit
        ,ABORT_ON_BADAVWARP   &  ! I: T=> abort if AVwarp cannot be determined
        ,ABORT_ON_BADZ        &  ! I: T=> abort on bad z in GET_KCOR8
        ,ABORT_ON_BADKCOR     &  ! I: T=> affects only the init stage (RDKCOR)
        ,ABORT_ON_BADSURVEY   &  ! I: T=> abort if SURVEY changes
        ,ABORT_ON_BADFILTER   &  ! I: T=> abort if fit band is not in SURVEY filters
        ,ABORT_ON_MARGPDF0    &  ! I: T=> abort if marginalized pdf=0 everywhere
        ,ABORT_ON_TRESTCUT    &  ! I: T=> abort if any Trest cut is set (for photoz)
        ,ABORT_ON_DUPLCID     &  ! I: T=> abort on duplicate CID
        ,ABORT_ON_DUPLMJD     &  ! I: T=> abort on repeat MJD+band (Jun 2017)
        ,ABORT_ON_NOPKMJD     &  ! I: T=> abort if no PKMJDINI (see OPT_SETPKMJD)
        ,USE_MINOS            &  ! I: T=> use MINOS instead of MINIMIZE
        ,LDMP_SNFAIL          &  ! I: T => dump reason for each SN failure
        ,LDMP_SNANA_VERSION   &  ! I: T => dump SNANA version and SNANA_DIR
        ,LDMP_AVWARP          &  ! I: dump GET_AVWARP8 (debug only)
        ,LDMP_KCORFUN         &  ! I: dump KCORFUN8 (debug only)
        ,LDMP_SATURATE        &  ! I: dump saturated observations
        ,USESIM_SNIA          &  ! I: default True -> process simulated SNIa
        ,USESIM_NONIA         &  ! I: default True -> process simulated nonIa
        ,USESIM_TRUEFLUX      &  ! I: SNLC_FLUXCAL -> SIM_FLUXCAL
        ,USESIM_REDSHIFT      &  ! I: replace REDSHIFT with SIM_REDSHIFT
        ,LPROB_TRUEFLUX       &  ! I: T=> compute F-Ftrue chi2 and PROB
        ,RESTORE_WRONG_VPEC   &  ! I: restore incorrect VPEC sign correction
        ,RESTORE_OVERRIDE_ZBUG  &  ! I: restore z=zCMB insteead of zHEL
        ,RESTORE_MWEBV_ERR_BUG  &  ! I: bug used MWEBV/6 if MWEBV_ERR < 0.001
        ,RESTORE_DES5YR         &  ! I: restore features for DES-SN5YR/V24
        ,FORCE_STDOUT_BATCH    ! I: force all stdout in batch mode

! LC plot info

    INTEGER  & 
        MXLC_PLOT       &  ! I: max number of plots to make with SNLCPAK
       ,NCCID_PLOT     ! I: size of SNCCID_PLOT array
    CHARACTER  & 
        SNCCID_PLOT(MXLISTNML)*(MXCHAR_CCID) ! I: string-list of CCIDs to plot
    REAL  & 
         MJDPERIOD_PLOT   &  ! I: fold LC onto this period (periodic transients)
        ,MJDSHIFT_PLOT   ! I: shift MJD for LC plot output

    REAL  & 
         DTOBS_MODEL_PLOT   &  ! I: binning (days) for overlaid best-fit model
        ,zTOL_HELIO2CMB    ! I: give warning if  zCMB and zHEL are off

! Oct 2014: variables to check for multi-season transient activity

    INTEGER  & 
         MULTISEASON_OPTMASK      ! option(s) for GET_MULTISEASON
    REAL  & 
         MULTISEASON_TGAP               &  ! new season if time-gap > this
        ,MULTISEASON_NREJECT_OUTLIER    &  ! num outliers to reject per season
        ,MULTISEASON_CHI2RED_ACTIVE    ! min chi2red for active season

    REAL  & 
         SIM_PRESCALE         &  ! I: pre-scale applied only to SIM
        ,VPEC_ERR_OVERRIDE    &  ! I: override VPEC_ERR in data header (km/sec)
        ,SNTABLE_APPEND_VALUE    &  ! I: value for SNTABLE_APPEND_VARNAME
        ,QUANTILE_ZERRMIN    ! I: use quantile prior for zerr> this; else use Guass prior

! variables to control estimate of PEAKMJD
    INTEGER  & 
         OPT_SETPKMJD            ! I: bit-mask options  to determined PKMJD
    REAL  & 
         SNRCUT_SETPKMJD          &  ! I: define SNR for Fmax in clump method
        ,MJDWIN_SETPKMJD          &  ! I: MJD window for Fmax in clump method
        ,SHIFT_SETPKMJD          ! I: shift initial PKMJD for systematic test

! define reference cosmological parameters (double precision!)
! Aug 7, 2009: each parameter has dimension 2 for value & error

    REAL*8  & 
         H0_REF(2)     &  ! I: reference H0  (70 => 70 km/s/MPc)
        ,OLAM_REF(2)   &  ! I: OMEGA_LAMBDA   uncertainty (for x0 prior)
        ,OMAT_REF(2)   &  ! I: OMEGA_MATTER   uncertainty (for x0 prior)
        ,ORAD_REF(2)   &  ! I: OMEGA_RAD   uncertainty (for x0 prior)
        ,W0_REF(2)     &  ! I: DE eq stat w = p/rho   uncertainty (for x0 prior)
        ,WA_REF(2)    ! I: DE wa [a = 1/(1+z)]   uncertainty

    INTEGER   NSIMVAR_CUTWIN
    REAL      SIMVAR_CUTWIN(2,MXCUT_PRIVATE)
    CHARACTER SIMVAR_CUTWIN_LIST(MXCUT_PRIVATE)*40



    NAMELIST / SNLCINP /  & 
            VERSION_PHOTOMETRY, VERSION_PHOTOMETRY_WILDCARD  & 
          , VERSION_REFORMAT_FITS, VERSION_REFORMAT_TEXT  & 
          , PRIVATE_DATA_PATH, FILTER_UPDATE_PATH  & 
          , NONSURVEY_FILTERS, SNRMAX_FILTERS, VPEC_ERR_OVERRIDE  & 
          , FILTER_REPLACE, FILTLIST_LAMSHIFT  & 
          , JOBSPLIT, JOBSPLIT_EXTERNAL, SIM_PRESCALE, MXLC_FIT  & 
          , OPT_YAML  & 
          , OPTSIM_LCWIDTH, OPT_REFORMAT_SPECTRA, OPT_REFORMAT_TEXT  & 
          , OPT_REFORMAT_SALT2, REFORMAT_KEYS, OPT_REFORMAT_FITS  & 
          , SNMJD_LIST_FILE, SNMJD_OUT_FILE, MNFIT_PKMJD_LOGFILE  & 
          , rootfile_out, textfile_prefix  & 
          , SNTABLE_LIST, SNTABLE_FILTER_REMAP, MARZFILE_OUT  & 
          , SNTABLE_APPEND_VARNAME, SNTABLE_APPEND_VALUE  & 
          , CALIB_FILE, KCOR_FILE, OVERRIDE_RESTLAM_BOUNDARY  & 
          , USERTAGS_FILE  & 
          , VPEC_FILE, HEADER_OVERRIDE_FILE, SIM_HEADER_OVERRIDE_FILE  & 
          , EPOCH_IGNORE_FILE, OUT_EPOCH_IGNORE_FILE, NONLINEARITY_FILE  & 
          , SNCID_LIST_FILE, OPT_SNCID_LIST, OPT_VPEC_COR  & 
          , SIMLIB_OUT, SIMLIB_OUTFILE, SIMLIB_ZPERR_LIST  & 
          , OPT_SIMLIB_OUT, SIMLIB_OUT_TMINFIX  & 
          , NFIT_ITERATION, MINUIT_PRINT_LEVEL, INTERP_OPT, USE_MINOS  & 
          , OPT_SETPKMJD, QUANTILE_ZERRMIN  & 
          , SNRCUT_SETPKMJD, MJDWIN_SETPKMJD, SHIFT_SETPKMJD, DEBUG_FLAG  & 
          , LSIM_SEARCH_SPEC, LSIM_SEARCH_ZHOST  & 
          , LDMP_SNFAIL, LDMP_SNANA_VERSION, LDMP_SATURATE  & 
          , LTEST_KCOR, LTEST_INTERP, LTEST_MAG, LTEST_U3BAND  & 
          , USESIM_SNIA, USESIM_NONIA, USESIM_TRUEFLUX, USESIM_REDSHIFT  & 
          , USE_SNHOST_ZPHOT, USE_HOSTGAL_PHOTOZ  & 
          , RESTORE_WRONG_VPEC, RESTORE_OVERRIDE_ZBUG  & 
          , RESTORE_MWEBV_ERR_BUG, RESTORE_DES5YR  & 
          , REQUIRE_DOCANA, FORCE_STDOUT_BATCH  & 
          , LPROB_TRUEFLUX  & 
          , ABORT_ON_NOEPOCHS, ABORT_ON_BADAVWARP, ABORT_ON_NOPKMJD  & 
          , ABORT_ON_BADZ, ABORT_ON_BADKCOR, ABORT_ON_BADSURVEY  & 
          , ABORT_ON_MARGPDF0, ABORT_ON_TRESTCUT  & 
          , ABORT_ON_DUPLCID, ABORT_ON_DUPLMJD, ABORT_ON_BADFILTER  & 
          , H0_REF, OLAM_REF, OMAT_REF, W0_REF, WA_REF  & 
          , USE_MWCOR  & 
          , MXLC_PLOT, NCCID_PLOT, SNCCID_PLOT  & 
          , MJDPERIOD_PLOT, MJDSHIFT_PLOT  & 
          , DTOBS_MODEL_PLOT  & 
          , MULTISEASON_OPTMASK, MULTISEASON_TGAP, zTOL_HELIO2CMB  & 
          , MULTISEASON_NREJECT_OUTLIER, MULTISEASON_CHI2RED_ACTIVE  & 
          , PHOTFLAG_DETECT, PHOTFLAG_TRIGGER  & 
! 
! Below are cut-variables defined in SNCUTS
! 
         ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST  & 
         ,MXEVT_PROCESS, MXEVT_CUTS  & 
         ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE  & 
         ,SNTYPE_IGNORE, DETNUM_LIST, SIM_TEMPLATE_INDEX_LIST  & 
         ,SNFIELD_LIST  & 
         ,PHOTFLAG_MSKREJ, PHOTFLAG_BITLIST_REJECT  & 
         ,cutwin_cid, cutwin_sntype  & 
         ,cutwin_redshift, cutwin_redshift_err  & 
         ,cutwin_ra, cutwin_dec  & 
         ,cutwin_hostsep,   cutwin_Nepoch  & 
         ,cutwin_snrmax,    cutwin_snrmax2, cutwin_snrsum, cutwin_nfield  & 
         ,cutwin_mwebv,     cutwin_nseason_active  & 
         ,cutwin_Trestmin,   cutwin_Trestmax  & 
         ,cutwin_TrestRange, cutwin_Trest, cutwin_Tobs  & 
         ,cutwin_Tgapmax,  cutwin_T0gapmax  & 
         ,cutwin_Tobsmin,  cutwin_Tobsmax  & 
         ,cutwin_peakmjd,  cutwin_nobs_predetect, cutwin_mjd  & 
! 
         ,cutwin_psf, cutwin_zp,  CUTWIN_ZPADU, CUTWIN_ZPNPE  & 
         ,cutwin_zperr, cutwin_photprob,  cutwin_errtest  & 
         ,cutwin_sim_pull  & 
         ,cutwin_trest_trueflux, cutwin_nband_thresh  & 
         ,cutwin_nfilt_snrmax,   cutwin_nfilt_snrmax2  & 
         ,cutwin_nfilt_trestmin, cutwin_nfilt_trestmax  & 
         ,cutwin_trest2, cutwin_nfilt_trest2  & 
         ,cutwin_restlam, cutwin_lamrest, cutwin_lamobs  & 
         ,cutwin_snr_nodetect  & 
         ,PRIVATE_CUTWIN_STRING, PRIVATE_VARNAME_READLIST  & 
         ,SIMVAR_CUTWIN_STRING  & 
         ,EARLYLC_STRING, REQUIRE_EPOCHS_STRING, DUMP_STRING  & 
! 
         ,SNCUT_SNRMAX, SNCUT_HOST_SBFLUX, SNCUT_NOBS_TREST  & 
         ,EPCUT_SNRMIN, EPCUT_SKYSIG, EPCUT_PSFSIG  & 
         ,MAGOBS_SHIFT_PRIMARY,  MAGOBS_SHIFT_ZP  & 
         ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_ZP, FILTER_LAMSHIFT  & 
         ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS  & 
         ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX  & 
         ,FUDGE_FLUXERR_SCALE  & 
         ,FUDGE_MAG_ERROR, FUDGE_MAG_COVERR, SIM_FUDGE_MAG_ERROR  & 
         ,FLUXERRMODEL_FILE,SIM_FLUXERRMODEL_FILE,FLUXERRMODEL_OPTMASK  & 
         ,MAGCOR_FILE, SIM_MAGCOR_FILE, FUDGE_HOSTNOISE_FILE  & 
         ,RV_MWCOLORLAW, OPT_MWCOLORLAW, PARLIST_MWCOLORLAW, OPT_MWEBV  & 
         ,MWEBV_SCALE, MWEBV_SHIFT, MWEBV_FORCE  & 
         ,HOSTGAL_ZPHOT_SHIFT, HOSTGAL_ZSPEC_SHIFT  & 
         ,HOSTGAL_PHOTOZ_SHIFT, HOSTGAL_SPECZ_SHIFT    &  ! obsolete
         ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR  & 
         ,CUTWIN_MJD_EXCLUDE, CUTWIN_TOBS_PREDETECT


! ---------- end of SNLCINP ---------


  END MODULE SNLCINP_NML

! =====================================================================
  MODULE FITSCOM
    USE SNPAR
    IMPLICIT NONE

! define external C functions to read from fits files

! max length of catenated strings from all epochs;
! e.g., FIELD
! Aug 31 2023: 40k -> 80k
! Nov 16 2023: 80k -> 160k

    INTEGER, PARAMETER :: MXLEN_EPSTRING = 160000

    INTEGER  & 
         RD_SNFITSIO_PREP  & 
        ,RD_SNFITSIO_STR   &  ! read string
        ,RD_SNFITSIO_INT   &  ! read int
        ,RD_SNFITSIO_SHT   &  ! read short int
        ,RD_SNFITSIO_FLT   &  ! read float
        ,RD_SNFITSIO_DBL  ! read double

    EXTERNAL  & 
         RD_SNFITSIO_PREP  & 
        ,RD_SNFITSIO_STR  & 
        ,RD_SNFITSIO_INT  & 
        ,RD_SNFITSIO_SHT  & 
        ,RD_SNFITSIO_FLT  & 
        ,RD_SNFITSIO_DBL  & 
        ,SET_RDMASK_SNFITSIO ! void

! Define index-lookup arrays used for faster reading of fits files
    INTEGER, PARAMETER :: MXFITSPAR = 400

    INTEGER  & 
         INDXFITS_HEAD(MXFITSPAR)     &  ! used in RDSNFITS_HEAD
        ,INDXFITS_EPCUTS(MXFITSPAR)   &  ! used in RDSNFITS_EPCUTS
        ,INDXFITS_PHOT(MXFITSPAR)     &  ! used in RDSNFITS_PHOT
        ,IPARFITS_HEAD  & 
        ,SNFITSIO_CODE_IVERSION     ! FITSIO version used to write data

! defin epoch mask for which

    INTEGER  & 
         EPOCH_RDMASK_SNFITSIO(MXEPOCH)  ! epochs to keep/reject

! define long fits-string returned by RD_SNFITSIO_STR
    CHARACTER STRFITS*(MXLEN_EPSTRING)



  END MODULE FITSCOM

! =====================================================================
  MODULE FILTCOM
    USE SNPAR
    IMPLICIT NONE

! filter bandpasses.
! Jan 16 2017: MXLAMBIN_FILT -> 4000 (was 3000) for JWST test
! May 31 2024 sync MXLAMBIN_PRIM = MXLAMBIN_SNSED
! 
    INTEGER, PARAMETER ::                 & 
         MXLAMBIN_FILT = 4000             &  ! max number of lambda bins per filter
        ,MXLAMBIN_PRIM = MXLAMBIN_SNSED      ! max lambda bins for primary

! Define filter chars A-Z, a-z and 0-9, but do not change order of
! original FILTDEF characters so that we can still use older Kcor files.
! 
!    CHARACTER FILTDEF_STRING*100
    CHARACTER,  PARAMETER :: FILTDEF_STRING*100 =  &  
        "ugrizYJHK UBVRIXy0123456789 abcdef " //   &  !  32
        "ACDEFGLMNOPQSTWZ " //                     &  ! +16
        "hjklmnopqstvwx "                             ! +14 = 62


    CHARACTER  & 
         FILTOBS_NAME(MXFILT_ALL)*40    &  ! full filternames
        ,FILTREST_NAME(MXFILT_ALL)*40   &  ! full filternames
        ,SURVEY_FILT_NAME(MXFILT_ALL)*60     &  ! survey(s) for this band
        ,PRIMARY_NAME*40                &  ! name of primary
        ,NONSURVEY_FILTERS_ADD*80      ! NONSURVEY_FILTERS that were added

    LOGICAL  & 
         LFILTDEF_OBS(MXFILT_ALL)   &  ! filter-trans define in KCOR file
        ,LFILTDEF_REST(MXFILT_ALL)  &  ! filter-trans define in KCOR file
        ,LFILTDEF_NONSURVEY(MXFILT_ALL)  &  ! part of NONSURVEY_FILTERS
        ,LFILTDEF_SNRMAX(MXFILT_ALL)     &  ! use for SNRMAX cut
        ,FREEZE_SURVEY_FILTERS      &  ! T=> do not re-read on 2nd version
        ,EXIST_BXFILT_OBS           &  ! T=> observer BX exists
        ,EXIST_BXFILT_REST         ! T=> idem for rest-frame

! define indices for legacy filters
    INTEGER  & 
         IFILT_SDSS_u, IFILT_SDSS_g, IFILT_SDSS_r  & 
        ,IFILT_SDSS_i, IFILT_SDSS_z  & 
        ,IFILT_BESS_U, IFILT_BESS_B, IFILT_BESS_V  & 
        ,IFILT_BESS_R, IFILT_BESS_I, IFILT_BESS_BX  & 
        ,IFILT_Y, IFILT_J, IFILT_H, IFILT_K

! define filter set for survey (observer) and for rest-frame.
! Each set goes 1 - NFILTDEF_[SURVEY,REST]

    INTEGER  & 
         NFILTDEF_SURVEY           &  ! no. survey (obs) filters
        ,NFILTDEF_READ             &  ! no. filters to read (excludes BX)
        ,NFILTDEF_OBS         &  ! actual NFILTDEF from kcor file; should = NFILTDEF_SURVEY
        ,IFILTDEF_MAP_SURVEY(MXFILT_OBS)      &  ! IFILT_OBS vs. sparse filter index
        ,IFILTDEF_INVMAP_SURVEY(MXFILT_ALL)  &  ! sparse index vs. IFILT_OBS
! 
        ,NFILTDEF_REST                  &  ! no. rest-filts defined in Kcor file
        ,IFILTDEF_MAP_REST(MXFILT_OBS)  &  ! IFILT_REST vs. sparse filter index
        ,IFILTDEF_INVMAP_REST(MXFILT_ALL)  &  ! sparse index vs. IFILT_REST
! 
        ,NFILTDEF_IGNORE_REST  & 
        ,IFILTDEF_IGNORE_REST(MXFILT_OBS)  ! nearest rest-filters to ignore

! variables used in fitter; these variables are over-written
! in FITPAR_PREP for each SN

    INTEGER  & 
         IFILT_REST_MAP(MXFILT_OBS)   &  ! idem for rest-frame filter
        ,NFILT_OBS_USEFIT                ! number of filters used in fit

    INTEGER*8  & 
        IFILT_OBS_EVAL_MASK(2,MXFILT_ALL) ! set for each filter

    character FILTLIST_FIT_USE*64  ! filter-list USED each fit

! define filter properties using MXFILT_ALL
! Nov 12, 2010: split FILT_XXX into FILTOBS_XXX and FILTREST_XXX
! 
    INTEGER  & 
         NLAMBIN_FILTOBS(MXFILT_ALL)    &  ! number of bins per transmission curve
        ,NLAMBIN_FILTREST(MXFILT_ALL)   &  ! number of bins per transmission curve
        ,NLAMBIN_PRIMARY               ! lambda bins for primary spec


    REAL  & 
         FILTOBS_TRANS(MXLAMBIN_FILT,MXFILT_ALL)    &  ! filter transmissions
        ,FILTOBS_TRANSMAX(MXFILT_ALL)               &  ! max trans
        ,FILTOBS_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)   &  ! corresponding lambda's (A)
        ,FILTOBS_LAMAVG(MXFILT_ALL)            &  ! effective central lambda
        ,FILTOBS_LAMRMS(MXFILT_ALL)            &  ! lambda-RMS
        ,FILTOBS_LAMRANGE(2,MXFILT_ALL)        &  ! min,max range in obs-frames
        ,FILTOBS_MAG_PRIMARY(MXFILT_ALL)       &  ! primary mag vs. ifilt_obs
        ,FILTOBS_ZPOFF_PRIMARY(MXFILT_ALL)     &  ! mag(native) - mag(synth)
        ,FILTOBS_ZPOFF_SNPHOT(MXFILT_ALL)      &  ! apply these ZPOFF to SNphot
! xxx mark      &  ,FILTOBS_LAMSHIFT(MXFILT_ALL)         ! lambda shift per filter
! 
        ,FILTREST_TRANS(MXLAMBIN_FILT,MXFILT_ALL)    &  ! filter transmissions
        ,FILTREST_TRANSMAX(MXFILT_ALL)               &  ! max trans
        ,FILTREST_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)   &  !  lambda's (A)
        ,FILTREST_LAMAVG(MXFILT_ALL)            &  ! effective central lambda
        ,FILTREST_LAMRMS(MXFILT_ALL)            &  ! lambda-RMS
        ,FILTREST_LAMRANGE(2,MXFILT_ALL)        &  ! min,max range in rest-frames
        ,FILTREST_MAG_PRIMARY(MXFILT_ALL)       &  ! primary mag vs. ifilt
        ,FILTREST_ZPOFF_PRIMARY(MXFILT_ALL)     &  ! mag(native) - mag(synth)
        ,FILTREST_LAMSHIFT(MXFILT_ALL)          &  ! lambda shift per filter
! 
        ,PRIMARY_FLUX(MXLAMBIN_PRIM)        &  ! primary spec for obs filters
        ,PRIMARY_LAM(MXLAMBIN_PRIM)

! for remap
    INTEGER NFILT_REPLACE, IFILTOBS_REPLACE(MXFILT_ALL)



  END MODULE FILTCOM

! =====================================================================
  MODULE FILTUPDCM
    USE SNPAR
    USE FILTCOM
    IMPLICIT NONE

! Nov 2010
! information related to updating filter transmissions for each SN
! 

    INTEGER  & 
        OPT_FILTER_UPDATE    &  ! global flag for filter update
       ,DOFLAG_FILTER_UPDATE(MXFILT_ALL)  ! track which filters to update

    CHARACTER  & 
        FILTER_UPDATE_TOPDIR*(MXCHAR_FILENAME)  !full path to filt-trans files


! ------------------------------------

! define updated filt-tran for SN and REF;
! the FILTOBS_TRANS arrays change for each SN
    REAL  & 
        FILTOBS_TRANSSN_UPD(MXLAMBIN_FILT,MXFILT_ALL)  & 
       ,FILTOBS_TRANSREF_UPD(MXLAMBIN_FILT,MXFILT_ALL)  & 
       ,FILTOBS_ZPOFF_UPD(MXFILT_ALL)

! filter subdir is
!    [PREFIX_UPD_FILTDIR][SNID][SUFFIX_UPD_FILTDIR]
! and filter filenames are
!    [PREFIX_TRANSSN][filt][SUFFIX_UPD_TRANSSN]

    character  & 
         PREFIX_UPD_TRANSSN*80          &  ! Filename prefix for SN trans
        ,PREFIX_UPD_TRANSREF*80  & 
        ,PREFIX_UPD_FILTDIR*80  & 
        ,SUFFIX_UPD_TRANSSN*80  & 
        ,SUFFIX_UPD_TRANSREF*80  & 
        ,SUFFIX_UPD_FILTDIR*80  & 
        ,FILTER_UPDATE_DIR*(MXCHAR_PATH)  ! directory with filter updates

    LOGICAL FILTINFO_UPD_SN, FILTINFO_UPD_REF

! define common for filter-set per SN


! Define common for filter set mapping.
! There are two maps:
!   MAP1: subdir name vs. index
!   MAP2: indx vs. SNID

    INTEGER, PARAMETER :: MXMAP1_FILTER_UPDATE = 20 
    INTEGER, PARAMETER :: MXSNLC_FILTUPD       = 20000  ! replaces MXSNLC_FILES

! MAP1:
    INTEGER  & 
         NMAP1_FILTER_UPDATE  & 
        ,MAP1_FILTER_UPDATE_INDX(MXMAP1_FILTER_UPDATE)
    CHARACTER  & 
         MAP1_FILTER_UPDATE_SUBDIR(MXMAP1_FILTER_UPDATE)*100

! MAP2:
! WARNING: works only up to MXSNLC_FILTUPD << MXSNLC

    INTEGER  & 
         NMAP2_FILTER_UPDATE  & 
        ,MAP2_FILTER_UPDATE_PTRMAP1(MXSNLC_FILTUPD)  ! points to MAP1 index
    CHARACTER  & 
         MAP2_FILTER_UPDATE_CCID(MXSNLC_FILTUPD)*(MXCHAR_CCID)  ! CCID list



  END MODULE FILTUPDCM

! =====================================================================
  MODULE SNANAFIT
    USE SNPAR
    IMPLICIT NONE

! store fit results for use in snana
! define arrays used in fit.

    INTEGER, PARAMETER ::      & 
         MXFITPAR    =  12     &  ! max number of fit parameters
        ,MXFITSTORE  =  800       ! max number of storage params

    INTEGER  & 
         NFITPAR(0:MXFILT_OBS)    &  ! Nfitpar for 0=all, or 'ifilt'
        ,NFITPAR_MN          &  ! Minuit NFITPAR includes NFIXPAR
        ,NFIXPAR             &  ! Number of fixed parameters (i.e, INISTP=0)
        ,NPLOTPAR             &  ! Number of extra FITVAL parameters to plot
        ,NEPOCH_ALL(0:MXFILT_OBS)     &  ! fit+rejected epochs per filter
        ,NEPOCH_FIT(0:MXFILT_OBS)     &  ! number of epochs fit per filter
        ,NEPOCH_REJECT(0:MXFILT_OBS)  &  ! number of epochs rejected by fit
        ,NEPOCH_REJECT2           &  ! total rejected ep within TREST window
        ,ERRTYPE(MXFITPAR)   &  ! MINOS or PARAB
        ,PARPTR_CHI2  & 
        ,PARPTR_IFILT_REST(MXFILT_ALL)  & 
        ,PARPTR_FITRESTMAG

    DOUBLE PRECISION  & 
         INIVAL(MXFITPAR)  & 
        ,INISTP(MXFITPAR)  & 
        ,INIBND(2,MXFITPAR)  & 
! 
        ,FITVAL(MXFITPAR,0:MXITER)  & 
        ,FITERR(MXFITPAR,0:MXITER)  & 
        ,FITERR_PLUS(MXFITPAR,MXITER)  & 
        ,FITERR_MINUS(MXFITPAR,MXITER)  & 
        ,FITERR_RATIO(MXFITPAR,MXITER)   &  ! abs(EPLUS/Eminus)
! 
        ,FITCHI2_MIN    &  ! min chi2 after minimization
        ,FITCHI2_QUIT   &  ! quit FCN calc if chi2 exceeds this value
        ,SIMCHI2_CHEAT  &  ! chi2 using exact sim values
        ,FITERRMAT_SPARSE(MXFITPAR,MXFITPAR)  & 
        ,FITERRMAT(MXFITPAR,MXFITPAR)  & 
        ,FITCORMAT(MXFITPAR,MXFITPAR)  & 
! 
        ,PDFVAL(MXFITPAR)  & 
        ,PDFERR(MXFITPAR)  & 
        ,PDFPROB2(MXFITSTORE)  & 
        ,PDFERRMAT(MXFITPAR,MXFITPAR)  & 
        ,PDFCORMAT(MXFITPAR,MXFITPAR)

! declare fit-storage arrays for all SN

    REAL  & 
         FITVAL_STORE(MXFITSTORE)  & 
        ,FITERR_STORE(MXFITSTORE)  & 
        ,PDFVAL_STORE(MXFITSTORE)   &  ! from average over PDF integral
        ,PDFERR_STORE(MXFITSTORE)   &  !
        ,PDFPROB2_STORE(MXFITSTORE)   &  !
        ,LCVAL_STORE(MXFITSTORE)    &  !  either FITVAL or PDFVAL
        ,LCERR_STORE(MXFITSTORE)    &  !  either FITERR or PDFERR
        ,LCFRACERRDIF_STORE(MXFITSTORE)   &  ! frac change on last 2 iter
! 
        ,LCCHI2_STORE(4)         &  !  1,2,3,4=> TOTAL,data,prior,sigma
        ,LCPROBCHI2_STORE(4)     &  !  idem for prob
        ,FITCHI2_STORE(4)        &  !  idem for fit minimum
        ,FITPROBCHI2_STORE(4)  & 
! 
!  define a few intermediate-fit variables to identify cuts to more
!  quickly remove events.
        ,FITPROB_ITER1           &  ! total FITPROB after ITER=1
        ,FITCHI2RED_INI          &  ! first guess at fit chi2red
        ,FITCHI2RED_INI2         &  ! fit chi2red after rough adjustment
! 
        ,INIVAL_STORE(MXFITSTORE)   &  ! store all inital values
        ,SIMVAL_STORE(MXFITPAR)      &  ! sim values
! 
        ,ERRMAX_BAD(MXFITPAR)      &  ! err < ERRMAX_BAD labelled bad
        ,MAG_XTMW_REF(MXFILT_ALL)  &  ! ref MWXT mags subtracted from data
! 
! ------ define FITANA_CUTWIN_XXX that depend on spec-vs-photoz fit
        ,FITANA_CUTWIN_TREST(2)  & 
        ,FITANA_CUTWIN_TRESTMIN(2)  & 
        ,FITANA_CUTWIN_TRESTMAX(2)  & 
        ,FITANA_CUTWIN_TREST2(2)  & 
        ,FITANA_CUTWIN_TGAPMAX(2)  & 
        ,FITANA_CUTWIN_T0GAPMAX(2)  & 
        ,FITANA_CUTWIN_SNRSUM(2)      &  ! Jul 2024
! ------ SNLC_XXX -> XXX_FIT re-evaluated after final fit
        ,TRESTMIN_FIT, TRESTMAX_FIT, TRESTRANGE_FIT  & 
        ,TOBSMIN_FIT,  TOBSMAX_FIT  & 
        ,TGAPMAX_FIT,  T0GAPMAX_FIT  & 
        ,SNRSUM_FIT  & 
        ,SNRMAX_FILT_FIT(MXFILT_OBS)       &  ! sparse ordering
        ,SNRMAX_SORT_FIT(MXFILT_OBS)

    INTEGER  & 
         ERRTYPE_STORE(MXFITPAR)   &  ! MINOS or PARABOLIC
        ,NERRTYPE(MXERRTYPE)       &  ! # params per error type
        ,NDOF_STORE(4)             &  ! 1,2,3 => TOTAL, DATA, PRIOR, SIGMA
        ,NDOF_PRIOR  & 
        ,NCALL_INTEGPDF  & 
        ,TIME_INTEGPDF  & 
        ,TIMESUM_INTEGPDF  & 
        ,TIMEAVG_INTEGPDF  & 
! ---- ISNLC_XXX  -> XXX_FIT  re-evaluated after fitting
        ,NFILT_TRESTMIN_FIT  & 
        ,NFILT_TRESTMAX_FIT  & 
        ,NFILT_TREST2_FIT  & 
        ,NFILT_SNRMAX_FIT     &  ! added Sep 30 2012
        ,NFILT_SNRMAX2_FIT    &  ! added Sep 20 2012
        ,MNSTAT_COV          ! ISTAT returned from MNSTAT: see minuit manual

    CHARACTER  & 
         PARNAME_STORE(MXFITSTORE)*(MXCHAR_PARNAME)  ! stored parameter names
    INTEGER  & 
         PAROPT_STORE(MXFITSTORE)  ! storage options for output

    LOGICAL  & 
         FLOATPAR(MXFITPAR)    &  ! T => ipar is floated in fit
        ,USEPDF_MARG           &  ! T => use margin. pdf-avg instead of fit-values
        ,LREPEAT_ITER          &  ! internal flag for repeated iteration
        ,LREPEAT_MINOS        ! repeat entire fit with MINOS (May 2018)




! ----------------------------
  END MODULE SNANAFIT

! =====================================================================
  MODULE PKMJDCOM
    USE SNPAR
    IMPLICIT NONE

! define bits for user option-mask OPT_SETPKMJD  (lsb = 0 )

    INTEGER, PARAMETER ::             & 
        OPTBIT_SETPKMJD_ANYFUN  = 0   &  ! (=1) use exact function from Bazin 09
       ,OPTBIT_SETPKMJD_POLYCOR = 1   &  ! (=2) include POLY cor in Bazin function
       ,OPTBIT_SETPKMJD_NOABORT = 2   &  ! (=4) do not abort if can't find PKMJD
       ,OPTBIT_SETPKMJD_FLUXMAX = 3   &  ! (=8) PEAKMJD=MJD(maxFlux) like JG
       ,OPTBIT_SETPKMJD_FLUXMAX2= 4   &  ! (=16) PKMJD for Fmax clump
       ,OPTBIT_SETPKMJD_FLUXMAX3= 5   &  ! (=32) idem, but log10(SNR) weight
       ,OPTBIT_SETPKMJD_TRIGGER = 6   &  ! (=64) PKMJD = MJD_TRIGGER
       ,OPTBIT_SETPKMJD_SAVEPAR = 9   &  ! (=512) save fit params to SNANA table
       ,OPTBIT_SETPKMJD_DUMP    =10   &  !(=1024) extra screen dump per SN/filter
       ,OPTBIT_SETPKMJD_SIM     =11      !(=2048) PKMJDINI = SIM_PEAKMJD
        

! Define fit params

    INTEGER, PARAMETER ::  & 
         IPAR_ISN    = 1   &  ! fixed param
        ,IPAR_FILT   = 2   &  ! fixed param
        ,IPAR_MJDMIN = 3   &  ! fixed param (May 2019)
        ,IPAR_MJDMAX = 4   &  ! fixed param (May 2019)
        ,IPAR_T0     = 5  & 
        ,IPAR_TRISE  = 6  & 
        ,IPAR_TFALL  = 7  & 
        ,IPAR_A0     = 8  & 
        ,IPAR_A1     = 9  & 
        ,IPAR_A2     = 10 

     REAL, PARAMETER ::  SNRMIN_forFLUXMAX  = 3.0   ! min SNR to consider for fluxmax

     CHARACTER  PKPARNAME(NPAR_ANYLC)*16
     DATA  PKPARNAME /  & 
           'ISN', 'IFILTOBS', 'MJDMIN', 'MJDMAX',  & 
           'T0', 'TRISE', 'TFALL','A0', 'A1', 'A2' /

! -----
    REAL*8 FITERRMAT_PKMJD(NPAR_ANYLC,NPAR_ANYLC,MXFILT_ALL)


    REAL*4  & 
         PKMJD_FIT(MXFILT_ALL)  & 
        ,PKMJD_ERR(MXFILT_ALL)  & 
        ,PKMJD_ERRMIN, PKMJD_ERRWGT   &  ! min and weighted error
        ,PKFLUX_FIT(MXFILT_ALL)  & 
        ,PKFLUX_ERR(MXFILT_ALL)  & 
        ,PKFLUX_ERRMIN, PKFLUX_ERRWGT   &  ! min and weighted error
        ,CHI2_FITPKMJD(MXFILT_ALL)

    INTEGER NFIT_PKMJD, NDOF_FITPKMJD(MXFILT_ALL)


  END MODULE PKMJDCOM

! =====================================================================
  MODULE WRS2COM
    USE SNPAR
    IMPLICIT NONE

! Nov 03, 2011:
! common block for translating SNANA format into SALT2 format.

    INTEGER  & 
         LEN_SURVEY, LEN_INST, LEN_MAGSYS,  LEN_PREFIX  & 
        ,NREPLACE, IMAP_REPLACE(MXFILT_ALL)  & 
        ,NEWKEY

    CHARACTER  & 
         NAMEof_INSTRUMENT*(MXCHAR_FILEWORD)  & 
        ,NAMEof_MAGSYS*(MXCHAR_FILEWORD)  & 
        ,NAMEof_SURVEY*(MXCHAR_FILEWORD)  & 
        ,NAMEof_PREFIX*(MXCHAR_FILEWORD)  & 
        ,NAMEof_REPLACE(2)*(MXCHAR_FILEWORD)  & 
        ,NEWKEY_NAME(20)*(MXCHAR_FILEWORD)  & 
        ,NEWKEY_ARG(20)*(MXCHAR_FILEWORD)



  END MODULE WRS2COM

! =====================================================================
  MODULE SNDATCOM

    USE SNPAR
    USE SNFILECOM
    USE CTRLCOM
    USE SNLCCOM
    USE SNCUTS
    USE SNHOSTCOM
    USE SNSIMCOM


  END MODULE SNDATCOM







! =============================================
! =====================================================================
! =====================================================================
! =====================================================================
    PROGRAM MAIN
! 
! Main snana program.
! 
    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT
    USE FILTCOM

    IMPLICIT NONE

    INTEGER   IERR, IVERS, NFIT_PER_SN, JDIFF

! funtions

#if defined(SNANA)
    INTEGER SNANA_GET_NLCPLOT
#endif

#if defined(PSNID)
    INTEGER PSNID_GET_NFIT
#endif

    EXTERNAL PRINT_CPUTIME

! ------------------- BEGIN -------------------

    JTIME_START = TIME()

    CALL PARSE_SNANA_ARGS()

    ISTAGE_SNANA = 0

! init some variables
    CALL INIT_SNVAR()

    CALL WARN_OLDINPUTS("init"//char(0), 0);

! read namelist
    CALL RDSNNML(IERR)
      IF ( IERR .NE. 0 ) GOTO 666

    ISTAGE_SNANA = ISTAGE_INIT

! check for dump flags (call to C code)
    CALL INIT_SNANA_DUMP(DUMP_STRING//char(0) )

! check if format is FITS or ASCII (set logicals  FORMAT_FITS[TEXT])
    CALL CHECK_FORMAT()

! ------------------------------------------------
!            Misc. inits for options

! check for CID list from &SNLCINP and/or separate file
    CALL INIT_SNCID_LISTS()

! check for SN,MJD list from file ... used for interpolation
    CALL RDFILE_SNMJDLIST()

! read optional tags
    CALL RDFILE_USERTAGS()

! init multi-season analysis
    CALL MULTISEASON(IFLAG_INI)  ! 1 = init flag

! check reformat options
    CALL INIT_REFORMAT(1)

! -----------------------------------------------
! prepare reading data; determine FITS or TEXT; read global info ...
    CALL INIT_READ_DATA()

! check reformat options
    CALL INIT_REFORMAT(2)

! continue with initialization that may depend on SURVEY and FILTERS.

! process namelist strings after we know
! the filters and survey
    CALL PROCSNNML(IERR)
      IF ( IERR .NE. 0 ) GOTO 666

    CALL INIT_FUDGE_FLUXCAL(IERR) ! fudge fluxes and/or errors
      IF ( IERR .NE. 0 ) GOTO 666

! intialize cutmask after reading data
    CALL INIT_CUTMASK(IERR)
      IF ( IERR .NE. 0 ) GOTO 666

! init early-LC selection
    CALL PARSE_EARLYLC_STRING()

! init REQUIRE_EPOCHS (Sep 2017)
    CALL PARSE_REQUIRE_EPOCHS_STRING()

! init cut-names now that we have the filter names
    CALL INIT_CUTNAMES(IERR)
      IF ( IERR .NE. 0 ) GOTO 666

    CALL READ_CALIB_WRAPPER(IERR)
      IF ( IERR .NE. 0 ) GOTO 666

! option to create SIMLIB from data (call after RDKCOR, Feb 2017)
    CALL MAKE_SIMLIB_FILE(1)

! check for valid MW GALextinct option  and print message
    CALL INIT_GALextinct()

    ISTAGE_SNANA = ISTAGE_RDSN

! parse SNTABLE_LIST string to know what tables & LCPLOTs to make
    CALL INIT_SNTABLE_OPTIONS()

! --------------------------------------------
! init output file; set OPTFIT based on program.
#if defined(SNANA)
    NFIT_PER_SN = SNANA_GET_NLCPLOT()
    ISJOB_SNANA = .TRUE.
    CALL MNDUMMY()   ! just to require minuit library.
#elif defined(SNFIT)
    NFIT_PER_SN = 1  ! only 1 fit per SN
    ISJOB_SNFIT = .TRUE.
    CALL MNDUMMY()   ! just to require minuit library.
#elif defined(PSNID)
    NFIT_PER_SN = PSNID_GET_NFIT()
    ISJOB_SNFIT = .TRUE.
    ISJOB_PSNID = .TRUE.
    IF ( OPT_SETPKMJD == 0 ) THEN  ! if default, then ...
      OPT_SETPKMJD = -1  ! DON'T estimate PKMJD before psnid (4/22/2014)
    ENDIF
#endif

    CALL INIT_OUTFILES(NFIT_PER_SN)

! Mar 2013, create subdir for monitor-init (CDTOPDIR below)
    CALL MAKEDIR_OUTPUT("MONINIT"//char(0), -1, 7 )

! --------------------------------------------
! plot a few Hubble diagrams for references
    IF ( USE_TABLEFILE )  CALL MON_HUBBLEREF ( IERR )

#if defined(SNFIT)
! initialize fitter before SN are selected
    CALL PRBANNER ( "CALL FITPAR_INI" )
    CALL FITPAR_INI ( IERR )
       IF ( IERR .NE. 0 ) GOTO 666
#endif

    CALL CDTOPDIR_OUTPUT(STDOUT_UPDATE)  ! climb out of MONINIT subdir (Mar 2013)

#if defined(PSNID)
    CALL PRBANNER ( "CALL PSNIDINI" )
    CALL PSNIDINI(IERR)
#endif

! -------------------------------------------------
! Let user initialize their stuff before SN are selected
! Aug 14 2014: move after CDTOPDIR call.

       IF ( IERR .NE. 0 ) GOTO 666

! check command-line args after all inits
    CALL CHECK_LINE_ARGS()

! pak global survey info after all initialization (Mar 24 2013)
    IF ( USE_TABLEFILE .and. NFIT_PER_SN>0 ) CALL SNLCPAK_SURVEY()


! ######################################################
! 
! 
! ######################################################


! ------------------------------
    JTIME_LOOPSTART = TIME()
    JDIFF = JTIME_LOOPSTART - JTIME_START

    CALL PRINT_CPUTIME(JTIME_START, "CPUTIME_INIT"//char(0),  & 
                  "second"//char(0), 0, 20,20)

    call flush(6)
! ------------------------------

    DO ivers = 1, N_VERSION
       CALL PROCESS_DATA_VERSION(IVERS)
    ENDDO

! ------------------------------
    JTIME_LOOPEND = TIME()
! ------------------------------

     CALL SNANA_END()

! ###################
!  graceful end here
    CALL EXIT(0)
! ###################


! abort-end here

666   CONTINUE
    c1err = 'Fatal error in MAIN'
    CALL MADABORT("MAIN", c1err, "Check what gaver abort.")

    CALL EXIT(0)

    END  ! end MAIN program

! ==========================================


! ==========================================


! ==========================================
    SUBROUTINE INIT_READ_DATA()

! Created Feb 2021
! Initialization driver for reading data files.
! Called once during global init stage.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IVERS, INIT_NUM

! ------------ BEGIN -----------

! one-time init for both FITS and TEXT (since we don't know format yet)
    INIT_NUM = 0

    INIT_NUM = INIT_NUM+1
    CALL RD_SNFITSIO_INIT(INIT_NUM)

    INIT_NUM = INIT_NUM+1
    CALL RD_SNTEXTIO_INIT(INIT_NUM)

! set flag to read header of 1st light curve to get global info
! (SURVEY, FILTERS ..) needed for initialization

    IVERS = 1

    CALL GETINFO_PHOTOMETRY(IVERS)

    CALL EXEC_READ_DATA(IVERS,OPTMASK_SNDATA_GLOBAL)

    RETURN
  END SUBROUTINE INIT_READ_DATA

! ===============================
    SUBROUTINE INIT_READ_OVERRIDE()

! Created May 2023
! Wrapper to call C-function RD_OVERRIDE_INIT
! 
! Oct 2023: pass REQ_DOC=1 to require DOCANA.
! Oct 2025: require ISDATA=T to read HEADER_OVERRIDE_FILE, and require
!           LSIM_SNANA=T to read SIM_HEADER_OVERRIDE_FILE
! 


    USE SNDATCOM
    USE SNLCINP_NML
    INTEGER LENF_DATA, LENF_SIM, LEN_PRIV, REQ_DOC
    LOGICAL ISDATA
    CHARACTER STR_TMP*(10*MXCHAR_FILENAME)

! ------------ BEGIN -----------

    CALL ENVreplace(HEADER_OVERRIDE_FILE)
    CALL ENVreplace(SIM_HEADER_OVERRIDE_FILE)

    LENF_DATA = INDEX(HEADER_OVERRIDE_FILE,' ') - 1
    LENF_SIM  = INDEX(SIM_HEADER_OVERRIDE_FILE,' ') - 1
    ISDATA    = .NOT. LSIM_SNANA

    IF ( LENF_DATA > 0 .AND. LENF_SIM == 0 .AND. LSIM_SNANA ) THEN
        print*,' '
        print*,' !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!='
        print*,'    WARNING: HEADER_OVERRIDE_FILE ignored for sim'
        print*,' !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!='
        print*,' '
    ENDIF

    REQ_DOC = 1  ! require DOCANA

    IF ( LENF_DATA > 0 .AND. ISDATA ) THEN
     STR_TMP = HEADER_OVERRIDE_FILE(1:LENF_DATA)//char(0)
     CALL RD_OVERRIDE_INIT(STR_TMP, REQ_DOC, LENF_DATA)
    ENDIF
    IF ( LENF_SIM > 0 .AND. LSIM_SNANA ) THEN
      STR_TMP = SIM_HEADER_OVERRIDE_FILE(1:LENF_SIM)//char(0)
      CALL RD_OVERRIDE_INIT(STR_TMP, REQ_DOC, LENF_SIM)
    ENDIF

    LEN_PRIV = INDEX(PRIVATE_VARNAME_READLIST,' ') -1
    IF ( LEN_PRIV > 0 ) THEN
       STR_TMP = PRIVATE_VARNAME_READLIST(1:LEN_PRIV)//char(0)
       CALL RD_PRIVATE_INIT(STR_TMP, REQ_DOC, LEN_PRIV) ! Sep 2023
    ENDIF

    RETURN
  END SUBROUTINE INIT_READ_OVERRIDE

! ==============================================
    SUBROUTINE PREP_VERSION_SUBFOLDER()

! Created Apr 7 2021
! If VERSION_PHOTOMETRY includes a sub-folder, e.g.,
!     'JLA2014/JLA2014_CSP'
! then reset VERSION_PHOTOMTRY = 'JLA2014_CSP' and update
! PRIVATE_DATA_PATH to be
!      SNDATA_ROOT/JLA2014          ! if input path is ''
!      [PRIVATE_DATA_PATH]/JLA2014  ! if input path is not ''
! 
! Constraints
!   - multiple versions must be under same folder
!   - more than one slash in VERSION_PHOTOMETRY -> abort
!   - $ in VERSION_PHOTOMETRY -> abort (no ENVs)
! 
! -------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER iver, jslash, LENV, LENP
    CHARACTER  & 
         VERSION_ORIG*(MXCHAR_VERSION)  & 
        ,VERSION*(MXCHAR_VERSION)  & 
        ,FOLDER*(MXCHAR_VERSION)

! ------------- BEGIN ------------

    DO iver = 1, N_VERSION
       VERSION_ORIG = VERSION_PHOTOMETRY(iver)
       jslash  = INDEX(VERSION_ORIG,'/')
       LENV    = INDEX(VERSION_ORIG,' ') - 1
       IF ( jslash > 0 ) THEN
          folder  = VERSION_ORIG(1:jslash-1)
          VERSION = VERSION_ORIG(jslash+1:LENV)
           VERSION_PHOTOMETRY(iver) = VERSION
           IF ( PRIVATE_DATA_PATH == ' ' ) THEN
             LENP = INDEX(SNDATA_ROOT,' ') - 1
             PRIVATE_DATA_PATH = SNDATA_ROOT(1:LENP)  & 
                    // '/lcmerge/' // folder
           ELSE
             LENP = INDEX(PRIVATE_DATA_PATH,' ') - 1
             PRIVATE_DATA_PATH = PRIVATE_DATA_PATH(1:LENP)  & 
                    // '/' // folder
           ENDIF

! xxxxxxxxx
!            print*,' xxx -------------------------------------- '
!            print*,' xxx VERSION_ORIG = ', VERSION_ORIG(1:40)
!            print*,' xxx --> VERSION  = ', VERSION(1:40)
!            print*,' xxx --> PATH = ', PRIVATE_DATA_PATH(1:80)
! xxxxxxxxx

       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE PREP_VERSION_SUBFOLDER


! ===========================================
    SUBROUTINE READ_CALIB_WRAPPER(IERR)
! 
! Ceated Nov 2022
! Wrapper subroutine to check if kcor/calib file is defined,
! and to call read-kcor function with appropriate args.
! [uses refactored C code]

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IERR  ! (O) return non-zero on error
    LOGICAL IGNOREFILE_fortran, USE_CALIB

    INTEGER LENF, OPT_FRAME, IFILT, IFILTDEF_OBS, IFILTDEF_REST
    CHARACTER*(MXCHAR_FILENAME) cCALIB_FILE, cFILTERS
    LOGICAL EXIST_KCOR_FILE

    REAL*8  & 
         D_MAGOBS_SHIFT_PRIM(MXFILT_ALL)  & 
        ,D_MAGREST_SHIFT_PRIM(MXFILT_ALL)  & 
        ,LAMAVG, LAMRMS, LAMMIN, LAMMAX, ZPOFF, MAG_PRIM

! external C codes
    REAL*8  GET_CALIB_ZPOFF_FILE, GET_CALIB_PRIMARY_MAG

    EXTERNAL  & 
          FLOAT2DOUBLE  & 
         ,READ_CALIB_DRIVER  & 
         ,GET_CALIB_FILTINDEX_MAP  & 
         ,GET_CALIB_FILTLAM_STATS  & 
         ,GET_CALIB_ZPOFF_FILE  & 
         ,GET_CALIB_PRIMARY_MAG  & 
         ,GET_CALIB_NFILTDEF  & 
         ,EXIST_CALIB_BXFILT  & 
         ,GET_KCOR_ZRANGE

! ---------- BEGIN ---------

    IERR = 0


! Feb 20 2025: check for legacy KCOR_FILE input
    EXIST_KCOR_FILE  = .NOT. IGNOREFILE_fortran(KCOR_FILE)
    IF ( EXIST_KCOR_FILE ) CALIB_FILE = KCOR_FILE

! ----------
    EXIST_CALIB_FILE = .NOT. IGNOREFILE_fortran(CALIB_FILE)

    if ( .NOT. EXIST_CALIB_FILE ) RETURN

    USE_CALIB = .false.
#if defined(SNFIT)
    USE_CALIB = .true.
#endif


! prep arguments ...

    CALL ENVreplace(CALIB_FILE)  ! Mar 2016

    LENF = INDEX(CALIB_FILE,' ') - 1
    cCALIB_FILE = CALIB_FILE(1:LENF) // char(0)
    cFILTERS    = SURVEY_FILTERS(1:NFILTDEF_SURVEY) // char(0)

! convert primary mag-shifts from sngl to double precision
    CALL FLOAT2DOUBLE(MXFILT_ALL, MAGREST_SHIFT_PRIMARY_FILT,  & 
                 D_MAGREST_SHIFT_PRIM )  ! <== returned
    CALL FLOAT2DOUBLE(MXFILT_ALL, MAGOBS_SHIFT_PRIMARY_FILT,  & 
                 D_MAGOBS_SHIFT_PRIM )  ! <== returned


! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! read kcor file and store information (C code)

    CALL READ_CALIB_DRIVER(cCALIB_FILE, cFILTERS, USE_CALIB,  & 
                 D_MAGREST_SHIFT_PRIM, D_MAGOBS_SHIFT_PRIM)

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! check for option to update filter trans for each SN
    CALL FILTER_UPDATE_INIT(IERR)

! - - - - - - -  - -
! load filter information from C struct to fortran common blocks

! OBS frame map is already loaded; get rest-frame map
    OPT_FRAME = OPT_FILTREST - 1
    CALL GET_CALIB_FILTINDEX_MAP(OPT_FRAME, NFILTDEF_REST,  & 
              IFILTDEF_MAP_REST, IFILTDEF_INVMAP_REST)

! fetch wavelength stats per band: AVG, RMS, MIN, MAX
    OPT_FRAME = OPT_FILTOBS - 1
    DO ifilt = 1, NFILTDEF_SURVEY
       ifiltdef_obs = IFILTDEF_MAP_SURVEY(ifilt)
       LFILTDEF_OBS(ifiltdef_obs) = .TRUE.

       CALL GET_CALIB_FILTLAM_STATS(OPT_FRAME, IFILTDEF_OBS,  & 
             LAMAVG, LAMRMS, LAMMIN, LAMMAX)
       FILTOBS_LAMAVG(ifiltdef_obs)     = sngl(LAMAVG)
       FILTOBS_LAMRMS(ifiltdef_obs)     = sngl(LAMRMS)
       FILTOBS_LAMRANGE(1,ifiltdef_obs) = sngl(LAMMIN)
       FILTOBS_LAMRANGE(2,ifiltdef_obs) = sngl(LAMMAX)

!     store primary mags
       MAG_PRIM = GET_CALIB_PRIMARY_MAG(OPT_FRAME,IFILTDEF_OBS)
       FILTOBS_MAG_PRIMARY(IFILTDEF_OBS) = MAG_PRIM

!       fetch ZPOFF from ZPOFF.DAT file in filter dir
       ZPOFF = GET_CALIB_ZPOFF_FILE(OPT_FRAME,ifiltdef_obs)
       FILTOBS_ZPOFF_SNPHOT(ifiltdef_obs) = sngl(ZPOFF)
    ENDDO

    OPT_FRAME = OPT_FILTREST - 1
    DO ifilt = 1, NFILTDEF_REST
       ifiltdef_rest = IFILTDEF_MAP_REST(ifilt)
       LFILTDEF_REST(ifiltdef_rest) = .TRUE.
       CALL GET_CALIB_FILTLAM_STATS(OPT_FRAME, IFILTDEF_REST,  & 
             LAMAVG, LAMRMS, LAMMIN, LAMMAX)
       FILTREST_LAMAVG(ifiltdef_rest)     = sngl(LAMAVG)
       FILTREST_LAMRMS(ifiltdef_rest)     = sngl(LAMRMS)
       FILTREST_LAMRANGE(1,ifiltdef_rest) = sngl(LAMMIN)
       FILTREST_LAMRANGE(2,ifiltdef_rest) = sngl(LAMMAX)

!     store primary mags
       MAG_PRIM = GET_CALIB_PRIMARY_MAG(OPT_FRAME,IFILTDEF_REST)
       FILTREST_MAG_PRIMARY(IFILTDEF_REST) = MAG_PRIM
    ENDDO

! misc tasks to store info

    CALL SET_ZPOFF()

! store number of obs-frame and rest-frame filters
    CALL GET_CALIB_NFILTDEF(NFILTDEF_OBS,NFILTDEF_REST)

! check if BX filter is defined for OBS and REST fram
    CALL EXIST_CALIB_BXFILT(EXIST_BXFILT_REST,EXIST_BXFILT_OBS)
!       print*,' xxx EXIST(BX) = ', EXIST_BXFILT_REST,EXIST_BXFILT_OBS

    RETURN
  END SUBROUTINE READ_CALIB_WRAPPER


! ==========================================
    SUBROUTINE EXEC_READ_DATA(IVERS,OPTMASK)

! Created Feb 2021
! Driver to read & analyze data for input version.
! Use logicals FORMAT_FITS and FORMAT_TEXT for format-specific
! function calls.
! This subroutine replaces legacy RDVERSION_FITS & RDVERSION_TEXT
! 
! Oct 17 2023: read spectra only if using a REFORMAT option.
! Sep 03 2024: set LRDFLAG_SPEC=T for SIMLIB_OUTFILE
! ------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER  & 
         IVERS      &  ! (I) ! version index to read
        ,OPTMASK   ! (I) ! indicates GLOBAL or ALL

! local var

    INTEGER   NSN_VERS, LEN_VERS, LEN_PATH, OPTRD
    INTEGER   IJOB, NJOBTOT, ISN, ISTAT
    INTEGER*8 JTIME_EVENTSTART
    LOGICAL   LRDFLAG_GLOBAL, LRDFLAG_ALL, LRDFLAG_SPEC
    LOGICAL   REFORMAT_LOCAL, WR_SIMLIB_OUTFILE
    CHARACTER cVERSION*(MXCHAR_VERSION), cPATH*(MXCHAR_PATH)
    CHARACTER FNAM*16

! functions
    LOGICAL  LDONE_CIDLIST, IGNOREFILE_fortran
    INTEGER  RD_SNFITSIO_PREP, RD_SNTEXTIO_PREP
    EXTERNAL RD_SNFITSIO_PREP, RD_SNTEXTIO_PREP, RD_SNFITSIO_EVENT
    EXTERNAL PRINT_CPUTIME

! --------------- BEGIN -----------

    FNAM = 'EXEC_READ_DATA'
    LRDFLAG_GLOBAL = ( OPTMASK .EQ. OPTMASK_SNDATA_GLOBAL)
    LRDFLAG_ALL    = ( OPTMASK .EQ. OPTMASK_SNDATA_ALL   )

    REFORMAT_LOCAL    = (REFORMAT_SNANA .or. REFORMAT_SPECTRA_ONLY)
    WR_SIMLIB_OUTFILE = (.not. IGNOREFILE_fortran(SIMLIB_OUTFILE) )

    LRDFLAG_SPEC   = (REFORMAT_LOCAL .or. USE_TABLEFILE_MARZ .or.  & 
                        WR_SIMLIB_OUTFILE)  & 
                       .and. (DEBUG_FLAG .NE. -333)  ! -333 suppresses reading SPEC

    ISNLC_VERSION    = IVERS

! get version and private data path for this version
    cVERSION = VERSION_PHOTOMETRY(ivers)
    LEN_VERS = INDEX(cVERSION,' ') - 1
    cVERSION = VERSION_PHOTOMETRY(ivers)(1:LEN_VERS) // char(0)

    LEN_PATH = INDEX(PRIVATE_DATA_PATH,' ') - 1
    cPATH    = PRIVATE_DATA_PATH(1:LEN_PATH) // char(0)

! read global info and store in SNDATA struct in sndata.h
    OPTRD = 0   ;
    IF ( FORMAT_FITS ) THEN
       IF ( LRDFLAG_GLOBAL )          OPTRD = OPTRD + 2
       IF ( .NOT. LRDFLAG_SPEC      ) OPTRD = OPTRD + 128  ! Jul 9 2025
       IF ( .NOT. REFORMAT_SIMTRUTH ) OPTRD = OPTRD + 256  ! Mar 2022
	 IF ( DO_GETINFO )              OPTRD = 64   ! Oct 2025

       NSN_VERS  = RD_SNFITSIO_PREP(OPTRD, cPATH, cVERSION,  & 
                          LEN_PATH, LEN_VERS )
    ELSE
       IF ( DEBUG_FLAG == 1024 ) OPTRD = OPTRD + 1024 ! generic DUMP
       if ( LRDFLAG_GLOBAL .or. IVERS > 1 ) then
          NSN_VERS = RD_SNTEXTIO_PREP(OPTRD, cPATH, cVERSION,  & 
                    LEN_PATH, LEN_VERS)
       endif
    ENDIF

    IF ( NSN_VERS > MXSNLC-1 ) THEN
      write(C1err,161) NSN_VERS, MXSNLC
161     format('NSN_VERS=',I8,' exceeds bound MXSNLC=',I8 )
      C2err = 'Check MXSNLC parameter in snana.F90'
      CALL MADABORT(FNAM, c1err, c2err )
    ENDIF

    N_SNLC_READ(IVERS) = NSN_VERS

! read SNDATA struct and tranfer global info to fortran variables
    if ( LRDFLAG_GLOBAL .OR. IVERS > 1 ) then
! xxx mark Oct 2025         CALL INIT_READ_OVERRIDE()      ! May 2023
       CALL RDGLOBAL_DRIVER()
       CALL INIT_SURVEY(NSN_VERS)     ! init a few things
       IF ( LRDFLAG_GLOBAL ) RETURN
    endif

! - - - - - - - -
! if we get here, read & process events

! strip off info for SPLITTING jobs among multiple CPUs.
    NJOBTOT = JOBSPLIT(2)    ! total number of split jobs
    IJOB    = JOBSPLIT(1)    ! do IJOB of NJOBTOT
    IF ( NJOBTOT .GT. 1 ) THEN
      write(6,20) JOBSPLIT
20      format(T5,'Process SPLIT-JOB ',I3,' of ', I3 )
    ENDIF

! ------------------------------------------------------
! LOOP OVER EVENTS

    DO 100 isn = IJOB, NSN_VERS, NJOBTOT ! every NJOBTOT'th SN

       IF ( N_SNLC_FITCUTS >= MXLC_FIT ) GOTO 100

       CALL INIT_SNLC()

       ! absolute index independent of cuts or SPLIT jobs;
       ! used as integer index if CID is a string (see CIDASSIGN)
       ABSO_INDEX = ABSO_OFFSET + isn

! read header for event and load SNDATA struct.
! Note that ISN is a fortran index starting at 1
      IF ( FORMAT_FITS ) THEN
         CALL RD_SNFITSIO_EVENT(OPTMASK_SNDATA_HEAD,ISN) ! C func
      ELSE IF ( FORMAT_TEXT ) THEN
         CALL RD_SNTEXTIO_EVENT(OPTMASK_SNDATA_HEAD,ISN) ! C func
      ENDIF

! transfer from C struct to global fortran variables
      CALL RDHEAD_DRIVER(istat)

      if ( ISTAT < 0 ) GOTO 100    ! failed cut on header var/CID

! Read observations for event and load SNDATA C-struct.
! Read optional SPECTRA and load GENSPEC C-struct.
! All reading is done here to avoid STORE_PARSE_WORDS conflicts below.
      IF ( FORMAT_FITS ) THEN
         CALL RD_SNFITSIO_EVENT(OPTMASK_SNDATA_OBS,ISN)
         IF ( LRDFLAG_SPEC ) THEN
            CALL RD_SNFITSIO_EVENT(OPTMASK_SNDATA_SPEC,ISN)
         ENDIF

      ELSE IF ( FORMAT_TEXT ) THEN
         CALL RD_SNTEXTIO_EVENT(OPTMASK_SNDATA_OBS, ISN)
         IF ( LRDFLAG_SPEC ) THEN
            CALL RD_SNTEXTIO_EVENT(OPTMASK_SNDATA_SPEC,ISN)
         ENDIF
      ENDIF

! transfer variables from C struct to fortran global variables

      CALL RDOBS_DRIVER()

      IF ( LRDFLAG_SPEC ) CALL RDSPEC_DRIVER(0)

! ---------------------------------------------------------------

       IF ( ISTAT .EQ. 0 ) THEN
          N_SNLC_PROC = N_SNLC_PROC + 1 ! number of processed LC
          call PRINT_RDSN()      ! print one-line summary
          CALL SNANA_DRIVER( ISN, N_SNLC_PROC, IVERS)
       ENDIF


! stop reading if/when all CIDs are processed
      IF ( LDONE_CIDLIST() ) GOTO 500

 100  CONTINUE          ! end ISN loop

 500  CONTINUE

! Mar 2025: indicate DONE for this version so that last FITS files can be closed
    IF ( FORMAT_FITS ) THEN
       CALL RD_SNFITSIO_EVENT(OPTMASK_SNDATA_DONE, 0)
    ENDIF

    ABSO_OFFSET = ABSO_OFFSET + N_SNLC_READ(IVERS)

    RETURN
  END SUBROUTINE EXEC_READ_DATA

! ==========================================
    SUBROUTINE GETINFO_PHOTOMETRY(ivers)

! Shell to use getInfo_PHOTOMETRY_VERION() function to
! get path to data, full name of list file and
! full name of readme file. Also construct
! SNDATA_PREFIX
! 


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IVERS  ! (I) version index

    INTEGER LEN_VERS, istat, L1, L2, L3
    CHARACTER*(2*MXCHAR_FILENAME)  & 
         VERSION, PATH, LIST_FILE, README_FILE
    CHARACTER FNAM*20, FIRST_WORD*60

    INTEGER   GETINFO_PHOTOMETRY_VERSION
    EXTERNAL  GETINFO_PHOTOMETRY_VERSION, CHECK_FILE_DOCANA

! ----------------- BEGIN --------------
    FNAM = "GETINFO_PHOTOMETRY"

    VERSION  = VERSION_PHOTOMETRY(ivers)
    LEN_VERS = INDEX(VERSION,' ') - 1
    VERSION  = VERSION(1:LEN_VERS)  // char(0)

    L1   =  INDEX(PRIVATE_DATA_PATH,' ') - 1
    PATH =  PRIVATE_DATA_PATH(1:L1) // char(0)

    ISTAT = GETINFO_PHOTOMETRY_VERSION (  & 
                     VERSION,                         &  ! input
                     PATH, LIST_FILE, README_FILE,    &  ! returned
                     LEN_VERS, L1, L2, L3 )

    L1 = INDEX(PATH,char(0)) - 1
    IF ( L1 > MXCHAR_PATH ) then
       CALL PRINT_PREABORT_BANNER(FNAM//char(0),40)
       print*,' SNDATA_PATH = ', PATH(1:L1)
       write(c1err,61) 'SNDATA_PATH', L1, MXCHAR_PATH
       write(c2err,62)
       CALL MADABORT(FNAM, C1err, C2err)
    ENDIF
    SNDATA_PATH = PATH(1:L1)           ! fill global

61    format('LEN(',A,')=',I3,  & 
            ' is too long (MXCHAR_PATH=',I3,')')
161   format('LEN(',A,')=',I3,  & 
            ' is too long (MXCHAR_FILENAME=',I3,')')
62    format('Check $SNDATA_ROOT');

    L2 = INDEX(LIST_FILE,char(0)) - 1
    IF ( L2 > MXCHAR_FILENAME ) then
       CALL PRINT_PREABORT_BANNER(FNAM//char(0),40)
       print*,' LIST_FILE = ', LIST_FILE(1:L2)

       write(c1err,161) 'LIST_FILE', L2, MXCHAR_FILENAME
       write(c2err,62)
       CALL MADABORT(FNAM, C1err, C2err)
    ENDIF

    SNLIST_FILE = LIST_FILE(1:L2)      ! fill global

    L3 = INDEX(README_FILE,char(0)) - 1
    IF ( L3 > MXCHAR_FILENAME ) then
       CALL PRINT_PREABORT_BANNER(FNAM//char(0),40)
       print*,' README_FILE = ', README_FILE(1:L3)
       write(c1err,161) 'README_FILE', L3, MXCHAR_FILENAME
       write(c2err,62)
       CALL MADABORT(FNAM, C1err, C2err)
    ENDIF
    SNREADME_FILE(IVERS) = README_FILE(1:L3)  ! fill global

    CALL CHECK_FILE_DOCANA(REQUIRE_DOCANA,README_FILE,MXCHAR_FILENAME)

! - - - - - - - -
! construct prefix = path/[version]
    SNDATA_PREFIX = SNDATA_PATH(1:L1) //  & 
             '/'  // VERSION(1:LEN_VERS)

    RETURN
  END SUBROUTINE GETINFO_PHOTOMETRY


! ==================================================
    SUBROUTINE DUMP_README(IVERS)

! Created Jun 2011
! Dump README file to stdout for this version.
! 
! Use global SNDATA_PREFIX to construct the name
! of the README file for this version.
! This code was moved from RDSNDATA so that it can be
! called for both TEXT and FITS format.
! 
! Feb 26 2015: SNFILE_README already filled, so don't fill it here.
! Feb 14 2020: fix first abort message, use FNAM.
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IVERS, IERR  ! (I/O)

    INTEGER LEN_README, IERROPEN, L1, LL

    CHARACTER  & 
         CMD*(MXCHAR_FILENAME+20)  & 
        ,FILENAME*(MXCHAR_FILENAME)  & 
        ,VERSION*(MXCHAR_VERSION)  & 
        ,FNAM*12

    EXTERNAL READ_YAML_VALS

! --------------- BEGIN --------------

    IERR = 0
    FNAM = 'DUMP_README'

    L1 = INDEX( SNDATA_PREFIX, ' ' ) - 1

    FILENAME = SNREADME_FILE(IVERS)

    OPEN(UNIT = LUNDMP, FILE = FILENAME,  & 
           IOSTAT = IERROPEN, STATUS='OLD')
    LEN_README = INDEX ( FILENAME, ' ' ) - 1
    CLOSE ( UNIT = LUNDMP )
    IF ( IERROPEN .NE. 0 ) THEN
      CALL PRINT_PREABORT_BANNER(FNAM(1:10)//char(0),10)
      print*,' FILENAME = ', FILENAME
      C1err = 'Could not find README file '
      C2err = 'For FILENAME above.'
      CALL MADABORT(FNAM, c1err, c2err )
    ENDIF


    LEN_README = INDEX ( FILENAME, ' ' ) - 1
    write(cmd,50) FILENAME(1:LEN_README)
50    format('cat ', A )

    VERSION  = VERSION_PHOTOMETRY(ivers)
    LL       = index(VERSION,' ') - 1

    GLOBAL_BANNER = " DUMP README FILE for VERSION "  & 
             // VERSION(1:LL)
    CALL PRBANNER ( global_banner(1:80) )

    print*, cmd
      CALL FLUSH(6)
      CALL SYSTEM ( cmd )
      CALL FLUSH(6)
    print*,' '


    RETURN
  END SUBROUTINE DUMP_README


! ===========================================
    SUBROUTINE RDGLOBAL_DRIVER()

! Created Feb 2021
! Refactor to use SNDATA struct; same for FITS and TEXT format.
! Transfer global variables from SNDATA struct to fortran:
! SURVEY, FILTERS, etc ...
! 
! Apr 24 2021: read SIM_BIASCOR_MASK
! Oct 08 2021: read SIM_MODEL_INDEX
! Dec 02 2022: fix subtle (harmless?) bug setting ISCORRECT_SIGN_VPEC
! Jul 23 2025: read ZP_FLUXCAL
! Oct 14 2025: call INIT_READ_OVERRIDE after LSIM_SNANA is set (instead of before RDGLOBAL)



! local var

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE FITSCOM
    USE SPECCOM

    IMPLICIT NONE

    REAL*8    DARRAY(MXEPOCH)
    INTEGER   OPT, LEN
    CHARACTER STRING*100
    LOGICAL   ISCORRECT_BUG, ISCORRECT_FIX

    LOGICAL   correct_sign_vpec_data
    EXTERNAL  correct_sign_vpec_data
! ------------------ BEGIN -------------

    OPT  = OPTMASK_SNDATA_GLOBAL
!     &     + OPTMASK_SNDATA_DUMP

    global_banner =  & 
         "RDGLOBAL: FETCH GLOBAL SURVEY INFO from HEADER"
    CALL PRBANNER(global_banner(1:60) )

    SNANA_VERSION_DATA = ''
    CALL FETCH_SNDATA_WRAPPER("SNANA_VERSION",  & 
           ONE, SNANA_VERSION_DATA, DARRAY, OPT)
    LEN = INDEX(SNANA_VERSION_DATA,' ') - 1

    ISCORRECT_SIGN_VPEC =  & 
            correct_sign_vpec_data(SNANA_VERSION_DATA(1:LEN)//char(0))

    CALL FETCH_SNDATA_WRAPPER("SURVEY", ONE,  & 
              SURVEY_NAME, DARRAY, OPT)

    CALL FETCH_SNDATA_WRAPPER("ZP_FLUXCAL", ONE,  & 
              STRING, DARRAY, OPT)
    ZP_FLUXCAL = SNGL(DARRAY(1))

    IF ( .NOT. FREEZE_SURVEY_FILTERS ) THEN
       CALL FETCH_SNDATA_WRAPPER("FILTERS",  & 
              ONE, SURVEY_FILTERS, DARRAY, OPT)
       NFILTDEF_READ = INDEX(SURVEY_FILTERS,' ')
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("DATATYPE",  & 
           ONE, DATATYPE, DARRAY, OPT)

    ! set global logical based on DATATYPE
    IF ( DATATYPE(1:4) .EQ. 'DATA' ) THEN
       ! do nothing
    ELSE IF ( DATATYPE(1:9)  .EQ. 'SIM_SNANA' ) THEN
       LSIM_SNANA = .TRUE.   ! SNANA sim
    ELSE IF ( DATATYPE(1:10) .EQ. 'SIM_MAGOBS' ) THEN
       LSIM_MAGOBS = .TRUE.  ! e..g, fakes overlaid on images
    ELSE
      C1err = 'Unrecognized DATATYPE = ' // DATATYPE
      C2err = 'Check DATATYPE key in FITS header.'
      CALL MADABORT("RDGLOBAL", c1err, c2err )
    ENDIF

! Mar 2022: check option to re-write sims to look like real data
    IF ( LSIM_SNANA .and. .not. REFORMAT_SIMTRUTH ) then
        LSIM_SNANA = .FALSE.
        DATATYPE   = 'DATA'
    ENDIF

! check HEADER_OVERRIDE here after LSIM_SNANA is set
    CALL INIT_READ_OVERRIDE()

    ISJOB_SIM = (LSIM_SNANA .or. LSIM_MAGOBS)

    CALL FETCH_SNDATA_WRAPPER("NXPIX", ONE, STRING, DARRAY, OPT)
    SNLC_NXPIX = int(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("NYPIX", ONE, STRING, DARRAY, OPT)
    SNLC_NYPIX = int(DARRAY(1))

    CALL RDGLOBAL_PRIVATE(OPT)

    CALL RDGLOBAL_ZPHOT_Q(OPT) 

    IF ( LSIM_SNANA ) THEN
       CALL FETCH_SNDATA_WRAPPER("SIMLIB_FILE", ONE, SIMLIB_FILENAME, DARRAY, OPT)

       CALL FETCH_SNDATA_WRAPPER("SIMLIB_MSKOPT",  ONE, STRING, DARRAY, OPT)
       SIMLIB_MSKOPT = int(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("SIMOPT_MWCOLORLAW",  ONE, STRING, DARRAY, OPT)
       SIMOPT_MWCOLORLAW = int(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("SIMOPT_MWEBV",  ONE, STRING, DARRAY, OPT)
       SIMOPT_MWEBV = int(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("SIM_MWRV", ONE, STRING, DARRAY, OPT)
       SIM_MWRV = SNGL(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("SIM_VARNAME_SNRMON", ONE, SIMNAME_SNRMON, DARRAY, OPT)

       CALL FETCH_SNDATA_WRAPPER("SIM_BIASCOR_MASK", ONE, STRING, DARRAY, OPT)
       SIM_BIASCOR_MASK = int(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("SIM_MODEL_INDEX",  ONE, STRING, DARRAY, OPT)
       SIM_MODEL_INDEX = INT(DARRAY(1))

       CALL RDGLOBAL_SIMSED(OPT)

       CALL FETCH_SNDATA_WRAPPER("PySEDMODEL",  & 
             ONE, STRING, DARRAY, OPT)
       CALL RDGLOBAL_PySEDMODEL(OPT,STRING)

       CALL RDGLOBAL_LCLIB(OPT)

       CALL RDGLOBAL_SIM_HOSTLIB(OPT)

    ENDIF  ! end LSIM_SNANA block


    RETURN
  END SUBROUTINE RDGLOBAL_DRIVER


! =======================================
    SUBROUTINE FETCH_SNDATA_WRAPPER(KEY, NARG, STRINGVAL, DVAL, OPT)

! Created Feb 7 2021
! Fortran wrapper to call C function fetch_SNDATA
! Note that NARG must  be 1 for string; NARG can be >1 for DVALs
! Return STRINGVAL or DVAL


    USE SNPAR
    USE FITSCOM

    IMPLICIT NONE

! function args
    CHARACTER KEY*(*)       ! (I) key name
    INTEGER   NARG          ! (I) number of args to return
    CHARACTER STRINGVAL*(*) ! (O) string value for key
    REAL*8    DVAL(*)       ! (O) double value(s) for key
    INTEGER   OPT           ! (I) options: 1-> dump
! local args
    INTEGER   LEN_KEY, LEN_STR, LEN_WHAT, COPYFLAG
    LOGICAL   L_GLOBAL, L_HEAD, L_OBS, L_DUMP, L_REQUIRE
    CHARACTER cSTRING*(MXLEN_EPSTRING), FNAM*22
    CHARACTER cKEY*60, c1err*80, c2err*80, WHAT*8
    EXTERNAL  COPY_SNDATA_GLOBAL, COPY_SNDATA_HEAD
    EXTERNAL  COPY_SNDATA_OBS

! ------------ BEGIN ----------

    FNAM = 'FETCH_SNDATA_WRAPPER'
    L_GLOBAL  = ( IAND(OPT,OPTMASK_SNDATA_GLOBAL)   > 0 )
    L_HEAD    = ( IAND(OPT,OPTMASK_SNDATA_HEAD)     > 0 )
    L_OBS     = ( IAND(OPT,OPTMASK_SNDATA_OBS )     > 0 )
    L_DUMP    = ( IAND(OPT,OPTMASK_SNDATA_DUMP)     > 0 )
    L_REQUIRE = ( IAND(OPT,OPTMASK_SNDATA_REQUIRE ) > 0 )

    LEN_KEY = INDEX(KEY//' ',' ') - 1
    cKEY    = KEY(1:LEN_KEY) // char(0)

    DVAL(1) = -9.0 ; STRINGVAL = 'NOTSET'

    COPYFLAG = -1  ! flag to copy from SNDATA struct to cSTRING or DVAL

    IF ( L_GLOBAL ) THEN
       CALL copy_SNDATA_GLOBAL(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
       WHAT = "GLOBAL"
    ELSE IF ( L_HEAD ) THEN
       CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
       WHAT = "HEAD"
    ELSE IF ( L_OBS ) THEN
       CALL copy_SNDATA_OBS(copyFlag, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
       WHAT = "OBS"
    ELSE
       write(C1ERR, '("Invalid OPT=",I6)' ) OPT
       C2ERR = 'KEY = ' // KEY
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF
!       print*,' xxx FETCH_WRAP: ',cKEY, cSTRING


    LEN_STR   = INDEX(cSTRING,char(0)) - 1
    STRINGVAL = cSTRING(1:LEN_STR)
    LEN_WHAT  = INDEX(WHAT,' ') - 1

    if ( LEN_STR < 0 ) then
       LEN_STR   = INDEX(cSTRING//' ', ' ') - 1
       write(C1ERR,161) LEN_STR, MXLEN_EPSTRING
161      format('Returned strlen=', I6,' ~=  MXLEN_EPSTRING=', I6)
	 write(C2ERR,162) KEY, NARG
162      format('KEY=',A, 4x, 'NARG=', I6)
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    if ( L_DUMP ) then
       write(6,66) WHAT(1:LEN_WHAT), KEY(1:LEN_KEY),  & 
              STRINGVAL(1:LEN_STR), DVAL(1)
 66      format(' xxx FETCH_SNDATA_',A,': ',  & 
              A,' = |', A, '|', 5x, G10.4 )
    endif

    RETURN
  END SUBROUTINE FETCH_SNDATA_WRAPPER

! =======================================
    SUBROUTINE FETCH_GENSPEC_WRAPPER(KEY, ISPEC, DVAL, OPT)

! Created Feb 18 2021
! Fortran wrapper to call C function copy_GENSPEC.
! Returns DVAL array. Calling function must be careful
! to provide adequate DVAL array size.
! 


    USE SNPAR
    USE FITSCOM

    IMPLICIT NONE

! function args
    CHARACTER KEY*(*)       ! (I) key name
    INTEGER   ISPEC         ! (I) spectrum index
    REAL*8    DVAL(*)       ! (O) double value(s) for key
    INTEGER   OPT           ! (I) options: 1-> dump
! local args
    INTEGER   LEN_KEY, COPYFLAG, cISPEC
    CHARACTER cKEY*60
    LOGICAL   L_DUMP
    EXTERNAL  COPY_GENSPEC

! ------------ BEGIN ----------

    LEN_KEY = INDEX(KEY//' ',' ') - 1
    cKEY    = KEY(1:LEN_KEY) // char(0)
    cISPEC  = ISPEC - 1  ! C index for GENSPEC

    DVAL(1) = -9.0 ! init output arg
    COPYFLAG = -1  ! flag to copy from GENSPEC struct to DVAL

    CALL copy_GENSPEC(COPYFLAG, cKEY, cISPEC, DVAL)

    L_DUMP   = ( IAND(OPT,OPTMASK_SNDATA_DUMP) > 0 )
    if ( L_DUMP ) then
       write(6,66) KEY(1:LEN_KEY), ISPEC, DVAL(1),DVAL(2),DVAL(3)
 66      format(' xxx FETCH_GENSPEC: ', A8, '(',I2,') = ', 3G12.3 )
    endif

    RETURN
  END SUBROUTINE FETCH_GENSPEC_WRAPPER


! =================================
    SUBROUTINE RDGLOBAL_PRIVATE(OPT)

! Refactored Feb 2021 to use SNDATA struct in sndata.h.
!  (previous name was RDPRIVATE_FITS)
! Read/parse names of PRIVATE keys from header if NVAR_PRIVATE > 0
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM

    IMPLICIT NONE

    INTEGER   OPT  ! (I)  1 -> dump flag for FETCH_SNDATA_WRAPPER

    INTEGER   IVAR, LEN
    REAL*8    DARRAY(MXVAR_PRIVATE)
    CHARACTER DUMSTRING*10, cnum*2, KEYNAME*20, KEYWORD*60
    CHARACTER VARNAME(20)*60
    LOGICAL   MATCH_VARNAME

! -------------- BEGIN ---------------

    CALL FETCH_SNDATA_WRAPPER("NVAR_PRIVATE",  & 
           ONE, DUMSTRING, DARRAY, OPT)
    NVAR_PRIVATE = int(DARRAY(1))

!       print*,' xxx OPT, NVAR_PRIVATE = ', OPT, NVAR_PRIVATE

    IF ( NVAR_PRIVATE .LE. 0 ) RETURN

    IF( NVAR_PRIVATE > MXVAR_PRIVATE ) THEN
       write(C1ERR,61) NVAR_PRIVATE, MXVAR_PRIVATE
 61      format('NVAR_PRIVATE=',I4,' exceeds MXVAR_PRIVATE=',I4 )
       C2ERR = 'Check XXX_HEAD.FITS file'
       CALL MADABORT("RDGLOBAL_PRIVATE", C1ERR, C2ERR)
    ENDIF
! ----------------------------------------------

    DO 100 ivar = 1, NVAR_PRIVATE

      if ( ivar < 10 ) then
         write(cnum, '(I1)') ivar
      else
         write(cnum, '(I2)') ivar
      endif

      KEYNAME = "PRIVATE" // CNUM
      CALL FETCH_SNDATA_WRAPPER(KEYNAME,  & 
           ONE, KEYWORD, DARRAY, OPTMASK_SNDATA_GLOBAL )

      PRIVATE_KEYWORD(ivar) = KEYWORD    ! load common block
      CALL PARSE_PRIVATE_KEYWORD(KEYWORD,ivar)

100   CONTINUE

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE RDGLOBAL_PRIVATE


! =================================
    SUBROUTINE RDGLOBAL_ZPHOT_Q(OPT)

! Created May 12 2022
! Read/parse names host galaxy photo-z quantiles that determine
! column names. E.g., percentile 20 means HOSTGAL_ZPHOT_Q020 exists.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM

    IMPLICIT NONE

    INTEGER   OPT  ! (I)  1 -> dump flag for FETCH_SNDATA_WRAPPER

    INTEGER   IVAR, q, PCT
    REAL*8    DARRAY(MXZPHOT_Q)
    CHARACTER DUMSTRING*10, cnum*2, KEYNAME*20, KEYWORD*80

! -------------- BEGIN ---------------

    CALL FETCH_SNDATA_WRAPPER("NZPHOT_Q",  & 
           ONE, DUMSTRING, DARRAY, OPT)
    SNHOST_NZPHOT_Q = int(DARRAY(1))

    IF ( SNHOST_NZPHOT_Q .LE. 0 ) RETURN

    IF( SNHOST_NZPHOT_Q > MXZPHOT_Q ) THEN
       write(C1ERR,61) SNHOST_NZPHOT_Q, MXZPHOT_Q
 61      format('NZPHOT_Q=',I4,' exceeds MXZPHOT_Q=',I4 )
       C2ERR = 'Check XXX_HEAD.FITS file'
       CALL MADABORT("RDGLOBAL_ZPHOT_Q", C1ERR, C2ERR)
    ENDIF
! ----------------------------------------------

    DO 100 ivar = 1, SNHOST_NZPHOT_Q

      q = ivar - 1 ! C index 0 to N-1
      write(cnum, '(I2.2)') q
      KEYNAME = "PERCENTILE_ZPHOT_Q" // CNUM
      CALL FETCH_SNDATA_WRAPPER(KEYNAME,  & 
           ONE, KEYWORD, DARRAY, OPTMASK_SNDATA_GLOBAL )

      PCT = INT(DARRAY(1))
      SNHOST_ZPHOT_PERCENTILE(ivar) = SNGL(DARRAY(1))

! load varname for each HOSTGAL match; e.g., HOSTGAL_ZPHOT_Q030
      write(VARNAME_ZPHOT_Q(1,ivar),102) 'HOSTGAL',  PCT
      write(VARNAME_ZPHOT_Q(2,ivar),102) 'HOSTGAL2', PCT
102     format(A,'_ZPHOT_Q', I3.3)

100   CONTINUE

    write(6,40) SNHOST_NZPHOT_Q
40    format(T5,'Found NZPHOT_Q = ', I3, ' quantiles for HOST-zPHOT')
    CALL FLUSH(6)

    RETURN
  END SUBROUTINE RDGLOBAL_ZPHOT_Q

! =======================================
    SUBROUTINE RDGLOBAL_SIMSED(OPT)

! Refactored Feb 2021 (previous name: RDFITSHEAD_SIMSED_LEGACY)
! Read/parse SIMSED keys from header if NPAR_SIMSED exists.
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE FITSCOM

    IMPLICIT NONE

    INTEGER OPT   ! (I) 1 -> dump


    INTEGER   ipar, ipar_read, LENWD, IPAR_OFF, OPT_LOCAL
    REAL*8    DARRAY(MXPAR_SIMSED)
    CHARACTER DUM*20, KEYNAME*40, KEYWORD*80, c2*2, PARNAME*60

! -------------- BEGIN ---------------

    OPT_LOCAL = OPT

    CALL FETCH_SNDATA_WRAPPER("SIMSED_NPAR",ONE,DUM,DARRAY,OPT_LOCAL)
    NPAR_SIMSED = int(DARRAY(1))

    IF ( NPAR_SIMSED .LE. 0 ) RETURN

! ----------------------------------------------
! SIMSED before Dec 2018 had SIMSED_PAR01  to SIMSED_PAR[NN].
! After refactor, it is SIMSED_PAR00 to N-1.
! For back compatibility, check for SIMSED_PAR00/

    CALL FETCH_SNDATA_WRAPPER("SIMSED_PAR00",  & 
           ONE, KEYNAME, DARRAY, OPT_LOCAL)

    IF ( INDEX(KEYNAME,' ') > 2 ) THEN
       IPAR_OFF = 0  ! 0 to N-1 for SIMSED index after Dec 2018
    ELSE
       IPAR_OFF = 1  ! 1-N before Dec 2018 (e.g., PLASTICC)
    ENDIF

    DO 100 ipar = 1, NPAR_SIMSED
      ipar_read = ipar + IPAR_OFF - 1   ! C-like index
      write(c2, '(I2.2)' ) ipar_read
      KEYNAME = "SIMSED_PAR" // C2 // char(0)
      CALL FETCH_SNDATA_WRAPPER(KEYNAME,ONE,KEYWORD,DARRAY,OPT_LOCAL)

      LENWD = INDEX(KEYWORD,' ') - 1
      CALL PARSE_PARENTHESES(KEYWORD(1:LENWD), PARNAME)
      SIMSED_KEYWORD(IPAR) = KEYWORD(1:LENWD) ! key(varname)
      SIMSED_PARNAME(IPAR) = PARNAME(1:20)    ! extracted varname

100   CONTINUE


    CALL FLUSH(6)

    RETURN
  END SUBROUTINE RDGLOBAL_SIMSED

! =======================================
    SUBROUTINE RDGLOBAL_PySEDMODEL(OPT,MODEL_NAME)

! Created Dec 10 2018
! Refactored Feb 2021 to use SNDATA struct.
! 
! Read/parse PySEDMODEL keys from header if NPAR_[BYOSED,SNEMO] exists.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER   OPT            ! (I) 1 -> dump
    CHARACTER MODEL_NAME*(*) ! (I) e.g., 'BYOSED', 'SNEMO'

    INTEGER   ipar, LENM, LENVAL, NPAR_LOCAL
    CHARACTER KEYNAME*80, KEYVAL*80, c2*2, PARNAME*60
    LOGICAL   LDMP, LVALID
    REAL*8    DARRAY(4)

! -------------- BEGIN ---------------

    LENM   = INDEX(MODEL_NAME,' ') - 1
    LVALID = .FALSE.
    IF ( MODEL_NAME(1:6) .EQ. 'BYOSED' ) LVALID = .TRUE.
    IF ( MODEL_NAME(1:5) .EQ. 'SNEMO'  ) LVALID = .TRUE.

!      print*,'    xxx MODEL_NAME, LVALID = ',
!     &      MODEL_NAME(1:LENM), LVALID

    if ( .not. LVALID ) RETURN

    KEYNAME   = MODEL_NAME(1:LENM) // "_NPAR"
    CALL FETCH_SNDATA_WRAPPER(KEYNAME, ONE, KEYVAL, DARRAY, OPT)

    LDMP = .FALSE.
    IF ( LDMP ) THEN
       print*,' xxx -------------------------------------- '
       print*,' xxx **** DUMP RDFITSHEAD_PySEDMODEL ***** '
       print*,' xxx MODEL_NAME = ', MODEL_NAME
       print*,' xxx KEYNAME    = ', KEYNAME(1:40)
       print*,' xxx KEYVAL     = ', KEYVAL(1:40)
    ENDIF

    NPAR_LOCAL = int(DARRAY(1))
    IF ( NPAR_LOCAL .LE. 0 ) RETURN

! ----------------------------------------------
! load global NPAR and MODEL_NAME
    NPAR_PySEDMODEL = NPAR_LOCAL
    PySEDMODEL_NAME = MODEL_NAME(1:LENM)

    DO 100 ipar = 1, NPAR_PySEDMODEL
      write(c2, '(I2.2)') ipar-1
      KEYNAME = MODEL_NAME(1:LENM) // "_PAR" // C2

      CALL FETCH_SNDATA_WRAPPER(KEYNAME, ONE, KEYVAL, DARRAY, OPT)

      LENVAL  = INDEX(KEYVAL,' ') - 1

      CALL PARSE_PARENTHESES(KEYVAL(1:LENVAL), PARNAME)
      PySEDMODEL_KEYWORD(IPAR) = KEYVAL(1:LENVAL) ! key(varname)
      PySEDMODEL_PARNAME(IPAR) = PARNAME(1:20)  ! extracted varname

100   CONTINUE

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE RDGLOBAL_PySEDMODEL

! =======================================
    SUBROUTINE RDGLOBAL_LCLIB(OPT)

! Feb 2018
! Read/parse LCLIB keys from header if NPAR_LCLIB exists.
! 


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER  OPT  ! (I)

    INTEGER   ipar, LENVAL
    REAL*8    DARRAY(4)
    CHARACTER KEYNAME*80, KEYVAL*80, c2*2, PARNAME*60

! -------------- BEGIN ---------------

    CALL FETCH_SNDATA_WRAPPER("LCLIB_NPAR",  & 
             ONE, KEYVAL, DARRAY, OPT)

    NPAR_LCLIB = INT(DARRAY(1))
    IF ( NPAR_LCLIB .LE. 0 ) RETURN

! ----------------------------------------------

    DO 100 ipar = 1, NPAR_LCLIB
      write(c2, '(I2.2)') ipar-1
      KEYNAME = "LCLIB_PAR" // C2
      CALL FETCH_SNDATA_WRAPPER(KEYNAME,  & 
             ONE, KEYVAL, DARRAY, OPT)

      LENVAL  = INDEX(KEYVAL,' ') - 1
      CALL PARSE_PARENTHESES(KEYVAL(1:LENVAL), PARNAME)
      LCLIB_KEYWORD(IPAR) = KEYVAL(1:LENVAL) ! key(varname)
      LCLIB_PARNAME(IPAR) = PARNAME(1:20)  ! extracted varname

100   CONTINUE
    CALL FLUSH(6)
    RETURN
  END SUBROUTINE RDGLOBAL_LCLIB

! =======================================
    SUBROUTINE RDGLOBAL_SIM_HOSTLIB(OPT)

! Created Feb 2021 [refactored from RDFITSHEAD_SIM_HOSTLIB]
! Read/parse SIM_HOSTLIB keys from header if NPAR_SIM_HOSTLIB exists.
! These are extra HOSTLIB params used in wgtmaps, etc ...



    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT  ! (I)

    INTEGER   ipar, LENVAL
    REAL*8    DARRAY(4)
    CHARACTER KEYNAME*80, KEYVAL*80, c2*2, PARNAME*60

! -------------- BEGIN ---------------

    CALL FETCH_SNDATA_WRAPPER("SIM_HOSTLIB_NPAR",  & 
             ONE, KEYVAL, DARRAY, OPT)

    NPAR_SIM_HOSTLIB = INT(DARRAY(1))
    IF ( NPAR_SIM_HOSTLIB .LE. 0 ) RETURN
! ----------------------------------------------

    DO 100 ipar = 1, NPAR_SIM_HOSTLIB
      write(c2, '(I2.2)' ) ipar-1  ! C like index
      KEYNAME = "SIM_HOSTLIB_PAR" // C2

      CALL FETCH_SNDATA_WRAPPER(KEYNAME,  & 
             ONE, KEYVAL, DARRAY, OPT)

      LENVAL = INDEX(KEYVAL,' ') - 1
      CALL PARSE_PARENTHESES(KEYVAL(1:LENVAL), PARNAME)
      SIM_HOSTLIB_KEYWORD(IPAR) = KEYVAL(1:LENVAL) ! key(varname)
      SIM_HOSTLIB_PARNAME(IPAR) = PARNAME ! extracted varname

100   CONTINUE
    CALL FLUSH(6)
    RETURN
  END SUBROUTINE RDGLOBAL_SIM_HOSTLIB

! ========================================
    SUBROUTINE RDGLOBAL_SIM_COSPAR()
! 
! Created Mar 2023
! Read simulated cosmology parameters from yaml-README
! (instead of from data header)

    USE SNDATCOM

    IMPLICIT NONE

    INTEGER LEN
    REAL*8 SIM_COSPAR(20)
    CHARACTER FNAM*24, KEYLIST*100, FILENAME*(MXCHAR_FILENAME)

! ------------- BEGIN --------------

    FNAM    = 'RDGLOBAL_SIM_COSPAR' // char(0)
    KEYLIST = 'OMEGA_MATTER,OMEGA_LAMBDA,w0_LAMBDA,wa_LAMBDA,MUSHIFT'  & 
                 //char(0)

    LEN = INDEX(SNREADME_FILE(1),' ') - 1
    FILENAME = SNREADME_FILE(1)(1:LEN) // char(0)

    CALL READ_YAML_VALS(FILENAME, KEYLIST, FNAM,  & 
                            SIM_COSPAR) ! return SIM_COSPAR

! store cosmo params in globals

! xxxxxxx mark delete May 19 2025 xxxxxxx
!      SIM_OM = SIM_COSPAR(1)
!      SIM_OL = SIM_COSPAR(2)
!      SIM_w0 = SIM_COSPAR(3)
!      SIM_wa = SIM_COSPAR(4)
!      SIM_MUSHIFT = SIM_COSPAR(5)
! xxxxxxxxxxxxx end mark xxxxxxx

    RETURN
  END SUBROUTINE RDGLOBAL_SIM_COSPAR

! ===========================================
    SUBROUTINE RDHEAD_DRIVER(ISTAT)

! Created Feb 2021 [refactored from PARSE_HEAD]
! Read header info and apply a few selection cuts to rapidly
! select small subsets:
! 
!   CID
!   SNTYPE
!   REDSHIFT_ERR ??
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER ISTAT  ! (O) < 0 -> reject event

    INTEGER OPT, igal, CID
    REAL*8  DARRAY(10), TMPCUT
    CHARACTER STRING*100, FNAM*14
    LOGICAL   USECID

! functions
    INTEGER  GET_IDSURVEY
    LOGICAL  PASS_SNTYPE
!  xxx mark      EXTERNAL RD_SNFITSIO_EVENT, RD_SNTEXTIO_EVENT

! -------------- BEGIN ------------

    ISTAT = 0
    FNAM  = 'RDHEAD_DRIVER'

! - - - - -
! fetch SNDATA values and load fortran globals

    OPT = OPTMASK_SNDATA_HEAD
!     &    + OPTMASK_SNDATA_DUMP

    DARRAY(1) = -999.0
    STRING    = ''

    CALL FETCH_SNDATA_WRAPPER("SUBSURVEY",  & 
           ONE, STRING, DARRAY, OPT)
    IF ( INDEX(STRING,' ') > 2 ) THEN
       SUBSURVEY_NAME = STRING(1:MXCHAR_SURVEY)
       IDSUBSURVEY    = GET_IDSURVEY(SUBSURVEY_NAME)
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SNID",  & 
           ONE, SNLC_CCID, DARRAY, OPT)
    ISNLC_LENCCID = INDEX(SNLC_CCID,' ') - 1
    call checkString_CCID(SNLC_CCID)   ! abort if illegal char in name

    CALL FETCH_SNDATA_WRAPPER("NAME_IAUC",  & 
           ONE, SNLC_NAME_IAUC, DARRAY, OPT)
    ISNLC_LENIAUC = INDEX(SNLC_NAME_IAUC,' ') - 1
    call checkString_CCID(SNLC_NAME_IAUC) ! abort if illegal char in name
    IF ( SNLC_NAME_IAUC(1:4) == 'NULL' ) SNLC_NAME_IAUC = 'NONE'  ! null confuses plot_table.py

! xxxxxx mark delete 9.12.2025 xxxxxxxxxxxxxx
!      IF ( WRTABLEFILE_IAUC .and. ISNLC_LENIAUC > 0 ) THEN
!         SNLC_CCID      = SNLC_NAME_IAUC
!         ISNLC_LENCCID  = ISNLC_LENIAUC
!      ENDIF
! xxxxxxxxxxxxxxxxx

    CALL FETCH_SNDATA_WRAPPER("NAME_TRANSIENT",   &  ! July 2024
           ONE, SNLC_NAME_TRANSIENT, DARRAY, OPT)
    ISNLC_LENNAME = INDEX(SNLC_NAME_TRANSIENT,' ') - 1
    call checkString_CCID(SNLC_NAME_TRANSIENT)   ! abort if illegal char in name

! - - - - - - - - - - - - - - - - - - - - - - - - - -
! Convert CCID to integer CID, and check if this CCID is selected.
    CALL PARSE_CID( SNLC_CCID, SNLC_NAME_IAUC,     &  ! inputs
                        CID, USECID)                ! returns args
    SNLC_CID  = CID   ! load integer CID
    if ( .NOT. USECID ) then
       ISTAT = ISTAT_SKIP;  RETURN
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - -


    CALL FETCH_SNDATA_WRAPPER("FAKE",  & 
           ONE, STRING, DARRAY, OPT)
    ISNLC_FAKE = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("RA",  ONE, STRING, DARRAY,  OPT)
    SNLC8_RA = DARRAY(1)

    CALL FETCH_SNDATA_WRAPPER("DEC", ONE, STRING, DARRAY, OPT)
    SNLC8_DEC = DARRAY(1)

    CALL FETCH_SNDATA_WRAPPER("PIXSIZE", ONE, STRING, DARRAY, OPT)
    SNLC_PIXSIZE = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("NXPIX", ONE, STRING, DARRAY, OPT)
    SNLC_NXPIX = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("NYPIX", ONE, STRING, DARRAY, OPT)
    SNLC_NYPIX = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SNTYPE", ONE, STRING, DARRAY, OPT)
    ISNLC_TYPE = INT(DARRAY(1))
    IF ( ISNLC_TYPE == -9 ) ISNLC_TYPE = 0 ! Jun 17 2022
    IF ( .NOT. PASS_SNTYPE(ISNLC_TYPE) ) THEN
       ISTAT = ISTAT_SKIP;  RETURN
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("MWEBV", ONE, STRING, DARRAY, OPT)
    SNLC_MWEBV = SNGL(DARRAY(1))
    SNLC_MWEBV_ERR = SNLC_MWEBV * sngl(XTMW_FRACERR)  ! default error

    CALL FETCH_SNDATA_WRAPPER("MWEBV_ERR", ONE, STRING, DARRAY, OPT)

    TMPCUT = 1.0E-12
    IF ( RESTORE_MWEBV_ERR_BUG ) TMPCUT = 0.001
    IF ( DARRAY(1) .GE. TMPCUT ) THEN
       SNLC_MWEBV_ERR = SNGL(DARRAY(1))  ! override default error
    ENDIF

! - - - - -
! read MJD-related variables (PEAK, trigger, first & last detection)
    CALL FETCH_SNDATA_WRAPPER("PEAKMJD", ONE, STRING, DARRAY, OPT)
    SNLC_SEARCH_PEAKMJD = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("MJD_TRIGGER",  & 
              ONE, STRING, DARRAY, OPT)
    SNLC8_MJD_TRIGGER = DARRAY(1)

    CALL FETCH_SNDATA_WRAPPER("MJD_DETECT_FIRST",  & 
              ONE, STRING, DARRAY, OPT)
    SNLC8_MJD_DETECT_FIRST = DARRAY(1)

    CALL FETCH_SNDATA_WRAPPER("MJD_DETECT_LAST",  & 
              ONE, STRING, DARRAY, OPT)
    SNLC8_MJD_DETECT_LAST = DARRAY(1)

! - - - - - -

    CALL FETCH_SNDATA_WRAPPER("NOBS", ONE, STRING, DARRAY, OPT)
    ISNLC_NEWMJD_FOUND = INT(DARRAY(1))

    CALL RDHEAD_REDSHIFT(OPT)

    DO igal = 1, MXSNHOST
       CALL RDHEAD_HOSTGAL(OPT,igal)
    ENDDO

    CALL RDHEAD_PRIVATE(OPT)

    IF ( LSIM_SNANA ) THEN
       CALL RDHEAD_SIM_SNANA(OPT)

!   reject event if user doesn't want to process true SNIa or NONIa
       if ( .NOT. USESIM_SNIA  .and. LSIM_TRUE_SNIa) THEN
	     ISTAT = ISTAT_SKIP ; RETURN
	 endif
       if ( .NOT. USESIM_NONIa .and. LSIM_TRUE_NONIa)  THEN
	    ISTAT = ISTAT_SKIP; RETURN
       endif

    ENDIF

    RETURN
  END SUBROUTINE RDHEAD_DRIVER

! ============================================
    SUBROUTINE RDHEAD_PRIVATE(OPT)

! Created Feb 2021 to read PRIVATE variables (pass from SNDATA struct)


    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM

    IMPLICIT NONE

    INTEGER OPT    ! OPT is arg for WRAPPER

    INTEGER ivar, NVAR
    REAL*8 DARRAY(4)
    CHARACTER STRING*12, KEY*60

! ------------ BEGIN ---------

    NVAR = NVAR_PRIVATE
    DO ivar = 1, NVAR
       KEY  = PRIVATE_KEYWORD(IVAR)
       CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
       PRIVATE_VALUE(IVAR) = DARRAY(1)
    ENDDO

    RETURN
  END SUBROUTINE RDHEAD_PRIVATE

! ============================================
    SUBROUTINE RDHEAD_HOSTGAL(OPT,igal)

! Created Feb 9 2021 [refactored from PARSE_HOSTGAL]
! Transfer HOSTGAL info from SNDATA struct (sndata.h)
! to common block variable.s
! 
! Mar 25 2021: fix bug setting EXIST_SNHOST_SB
! May 11 2021: read HOSTGAL_MAG (forgotten in I/O refactor)
! May 21 2021: read HOSTGAL_FLAG
! May 11 2022: read ZPHOT_Q


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER OPT, igal  ! OPT is arg for WRAPPER, igal is host index

    INTEGER   LENPRE, ifilt, q
    CHARACTER PREFIX*20, STRING*40, KEY*60, KEY_PREFIX*60
    REAL*8    DARRAY(MXFILT_OBS), SB, MAG, zq
    LOGICAL LTMP, LZQ

! -------------- BEGIN ----------

    CALL SET_HOSTGAL_PREFIX(IGAL,PREFIX,LENPRE)

! Start with items that only appear once;
! i.e., do not depend on host-match.

    IF ( IGAL == 1 ) THEN

       CALL FETCH_SNDATA_WRAPPER("HOSTGAL_NMATCH",  & 
              ONE, STRING, DARRAY, OPT)
       SNHOST_NMATCH = INT(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("HOSTGAL_NMATCH2",  & 
              ONE, STRING, DARRAY, OPT)
       SNHOST_NMATCH2 = INT(DARRAY(1))

       KEY_PREFIX = "HOSTGAL_SB_FLUXCAL"
       CALL RDHEAD_FILTERLOOP(KEY_PREFIX, SNHOST_SBFLUXCAL, OPT)
       DO ifilt = 1, NFILTDEF_READ
          SB = SNHOST_SBFLUXCAL(ifilt)
          if ( SB > -998.0 ) EXIST_SNHOST_SB = .TRUE.
       ENDDO

       KEY_PREFIX = "HOSTGAL_MAG"
       CALL RDHEAD_FILTERLOOP(KEY_PREFIX, SNHOST_MAGOBS, OPT)
       DO ifilt = 1, NFILTDEF_READ
          MAG   = SNHOST_MAGOBS(ifilt,1)
          if ( MAG > -998.0 ) EXIST_SNHOST_MAGOBS = .TRUE.
       ENDDO

       KEY_PREFIX = "HOSTGAL_SB_FLUXCAL"
       CALL RDHEAD_FILTERLOOP(KEY_PREFIX, SNHOST_SBFLUXCAL, OPT)
       DO ifilt = 1, NFILTDEF_READ
          SB = SNHOST_SBFLUXCAL(ifilt)
          if ( SB > -998.0 ) EXIST_SNHOST_SB = .TRUE.
       ENDDO

    ENDIF  ! end IGAL=1

! - - - - - - - -

    KEY    = PREFIX(1:LENPRE)//'_OBJID'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_OBJID(igal)  = INT8(DARRAY(1))
    DSNHOST_OBJID(igal) = DARRAY(1)

    KEY    = PREFIX(1:LENPRE)//'_FLAG'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_FLAG(igal)  = INT(DARRAY(1))

! - - -
    KEY    = PREFIX(1:LENPRE)//'_PHOTOZ'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_ZPHOT(igal)     = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_PHOTOZ_ERR'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_ZPHOT_ERR(igal) = SNGL(DARRAY(1))

! read ZPHOT_Q (May 2022)
    if ( SNHOST_NZPHOT_Q > 0 ) then
      DO q = 1, SNHOST_NZPHOT_Q
        KEY = VARNAME_ZPHOT_Q(igal,q)
        CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
        SNHOST_ZPHOT_Q(igal,q) = SNGL(DARRAY(1))
      ENDDO

      SNHOST_QZPHOT_MEAN(igal) = -9.0
      SNHOST_QZPHOT_STD(igal)  = -9.0

!          print*,' xxx ------------------------------------------- '
!          print*,' xxx CID ,igal = ', SNLC_CCID, igal
!          print*,' xxx SNHOST_ZPHOT_Q = ', SNHOST_ZPHOT_Q(igal,1:12)
!          call flush(6)
    endif

    IF ( SNHOST_ZPHOT(igal) > 0.0 ) EXIST_SNHOST_ZPHOT = .TRUE.
! - - - -

    KEY    = PREFIX(1:LENPRE)//'_SPECZ'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_ZSPEC(igal)     = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_SPECZ_ERR'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_ZSPEC_ERR(igal)     = SNGL(DARRAY(1))

! - - -

    KEY  = PREFIX(1:LENPRE)//'_RA'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST8_RA(igal) = DARRAY(1)

    KEY  = PREFIX(1:LENPRE)//'_DEC'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST8_DEC(igal) = DARRAY(1)

    KEY  = PREFIX(1:LENPRE)//'_SNSEP'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_ANGSEP(igal) = SNGL(DARRAY(1))
    IF ( SNHOST_ANGSEP(igal) > 0.0 ) EXIST_SNHOST_ANGSEP = .TRUE.

    KEY  = PREFIX(1:LENPRE)//'_DDLR'
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_DDLR(igal) = SNGL(DARRAY(1))
    IF ( SNHOST_DDLR(igal) > 0.0 )  EXIST_SNHOST_DDLR = .TRUE.

! - - - - - - - - - - - -  -

    KEY    = PREFIX(1:LENPRE)//'_LOGMASS' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGMASS(igal)   = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_LOGMASS_ERR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGMASS_ERR(igal)   = SNGL(DARRAY(1))

! xxx   if(SNHOST_LOGMASS(igal) > -9000.) EXIST_SNHOST_LOGMASS = .TRUE.
    EXIST_SNHOST_LOGMASS = .TRUE.

! - - - -

    KEY    = PREFIX(1:LENPRE)//'_LOGSFR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGSFR(igal)     = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_LOGSFR_ERR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGSFR_ERR(igal)     = SNGL(DARRAY(1))

! xxx   IF(SNHOST_LOGSFR(igal) > -9000.0) EXIST_SNHOST_LOGSFR = .TRUE.
    EXIST_SNHOST_LOGSFR = .TRUE.
! - - - -

    KEY    = PREFIX(1:LENPRE)//'_LOGsSFR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGsSFR(igal)     = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_LOGsSFR_ERR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_LOGsSFR_ERR(igal)     = SNGL(DARRAY(1))

! xxx   IF(SNHOST_LOGsSFR(igal) > -9000.0) EXIST_SNHOST_LOGsSFR = .TRUE.
    EXIST_SNHOST_LOGsSFR = .TRUE.  ! Jun 13 2024

! - - - -

    KEY    = PREFIX(1:LENPRE)//'_COLOR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_COLOR(igal)     = SNGL(DARRAY(1))

    KEY    = PREFIX(1:LENPRE)//'_COLOR_ERR' ;
    CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
    SNHOST_COLOR_ERR(igal)     = SNGL(DARRAY(1))

! xxx   IF(SNHOST_COLOR(igal) > -9000.0) EXIST_SNHOST_COLOR = .TRUE.
    EXIST_SNHOST_COLOR = .TRUE.   ! Jun 13 2024
! - - -
    KEY_PREFIX = PREFIX(1:LENPRE) // "_MAG"
    CALL RDHEAD_FILTERLOOP(KEY_PREFIX, SNHOST_MAGOBS(ifilt,igal), OPT)

    RETURN
  END SUBROUTINE RDHEAD_HOSTGAL

! ===========================================
    SUBROUTINE RDHEAD_FILTERLOOP(KEY_PREFIX, FARRAY, OPT)

! Created Feb 2021
! Check keys KEY_PREFIX_[band] and return float FARRAY(ifilt).
! Mar 23 2021: fix bug using LENPRE


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    CHARACTER KEY_PREFIX*(*)  ! (I) check keys KEY_PREFIX_[band]
    REAL      FARRAY(*)       ! (O) value vs. ifilt
    INTEGER   OPT             ! (I) argument for FETCH_SNDATA_WRAPPER

    INTEGER   IFILT, IFILT_OBS, LENPRE
    CHARACTER STRING*10, BAND*2, KEY*60
    REAL*8    DARRAY(4)

! --------------- BEGIN ---------

    LENPRE = INDEX(KEY_PREFIX//' ',' ') - 1

    DO ifilt     = 1, NFILTDEF_READ
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       band      = FILTDEF_STRING(ifilt_obs:ifilt_obs)
       KEY       = KEY_PREFIX(1:LENPRE) // '_' // band
       CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
       FARRAY(ifilt) = SNGL(DARRAY(1))

!         write(6,66) KEY, FARRAY(ifilt)     ! xxx
! 6       format(' xxx ', A20,' = ', F12.5 ) ! xxx

    ENDDO
    RETURN
  END SUBROUTINE RDHEAD_FILTERLOOP

! ===========================================
    SUBROUTINE RDHEAD_REDSHIFT(OPT)

! Created Feb 2021 [refactored from PARSE_REDSHIFT]
! Read/store redshift & vpec info
! Note tat REDSHIFT_FINAL is used for LC fitting.
! 
! Jun 5 2021: no longer apply REDSHIFT_FINAL_SHIFT here; see SET_zSHIFT().
! Aug 19 2024: check OPT_VPEC_COR option to disable VPEC
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT    ! (I) argument for FETCH_SNDATA_WRAPPER

    CHARACTER STRING*40
    REAL*8 DARRAY(10)

! ----------- BEGIN ------------

    CALL FETCH_SNDATA_WRAPPER("REDSHIFT_HELIO",  & 
             ONE, STRING, DARRAY, OPT)
    SNLC_zHELIO = SNGL(DARRAY(1))
    SNLC_REDSHIFT = SNLC_zHELIO  ! used for LC fitting

    CALL FETCH_SNDATA_WRAPPER("REDSHIFT_HELIO_ERR",  & 
             ONE, STRING, DARRAY, OPT)
    SNLC_zHELIO_ERR = SNGL(DARRAY(1))
    SNLC_REDSHIFT_ERR = SNLC_zHELIO_ERR

! - - - - -

    CALL FETCH_SNDATA_WRAPPER("REDSHIFT_FINAL",  & 
             ONE, STRING, DARRAY, OPT)
    SNLC_zCMB = SNGL(DARRAY(1))  !


    CALL FETCH_SNDATA_WRAPPER("REDSHIFT_FINAL_ERR",  & 
             ONE, STRING, DARRAY, OPT)
    SNLC_zCMB_ERR = SNGL(DARRAY(1))


! if no REDSHIFT_FINAL, try alternative REDSHIFT_CMB
! (beware that VPEC correction should NOT be included here)
    IF ( SNLC_zCMB < 0.0 ) THEN
       CALL FETCH_SNDATA_WRAPPER("REDSHIFT_CMB",  & 
               ONE, STRING, DARRAY, OPT)
       SNLC_zCMB = SNGL(DARRAY(1))

       CALL FETCH_SNDATA_WRAPPER("REDSHIFT_CMB_ERR",  & 
               ONE, STRING, DARRAY, OPT)
       SNLC_zCMB_ERR = SNGL(DARRAY(1))
    ENDIF

! - - - - - - - - - -
! assign light curve redshift if not already assigned by Z_HELIO
! i.e., light curve analysis should use Z_HELIO if it's defined;
! otherwise use Z_CMB

    IF ( SNLC_ZHELIO < 0.0 ) then
       SNLC_REDSHIFT      = SNLC_zCMB
       SNLC_REDSHIFT_ERR  = SNLC_zCMB_ERR
    ENDIF

! ---------------------------
! LENSDMU (Feb 2025)

    CALL FETCH_SNDATA_WRAPPER("LENSDMU", ONE, STRING, DARRAY, OPT)
    SNLC_LENSDMU = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("LENSDMU_ERR", ONE, STRING, DARRAY, OPT)
    SNLC_LENSDMU_ERR = SNGL(DARRAY(1))

! -------------------------
! VPEC

    CALL FETCH_SNDATA_WRAPPER("VPEC", ONE, STRING, DARRAY, OPT)
    SNLC_VPEC = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("VPEC_ERR", ONE, STRING, DARRAY, OPT)
    SNLC_VPEC_ERR = SNGL(DARRAY(1))

! Oct 26 2020: check option to fix VPEC sign convention in older FITS data
    if ( DOFIX_WRONG_SIGN_VPEC ) then
       SNLC_VPEC = -1.0 * SNLC_VPEC
    endif

    IF (  abs(VPEC_ERR_OVERRIDE) > .001 ) THEN
       SNLC_VPEC_ERR = VPEC_ERR_OVERRIDE  ! Jan 11 2018
    ENDIF

    if ( OPT_VPEC_COR == 0 ) then  ! Aug 19 2024
        SNLC_VPEC     = 0.0
	  SNLC_VPEC_ERR = 0.0
    endif

! redshift quality flag (real data only)
    CALL FETCH_SNDATA_WRAPPER("REDSHIFT_QUALITYFLAG",  & 
          ONE, STRING, DARRAY, OPT)
    ISNLC_zFLAG = INT(DARRAY(1))

! mask redshift source (march 14 2024)
    CALL FETCH_SNDATA_WRAPPER("MASK_REDSHIFT_SOURCE",  & 
          ONE, STRING, DARRAY, OPT)
    ISNLC_zSOURCE = INT(DARRAY(1))

    RETURN
  END SUBROUTINE RDHEAD_REDSHIFT

! =============================
    SUBROUTINE RDHEAD_SIM_SNANA(OPT)


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT    ! (I) argument for FETCH_SNDATA_WRAPPER

    INTEGER ipar
    CHARACTER STRING*100, KEY*80
    REAL*8 DARRAY(10)

! functions
    INTEGER ISMODEL_SNIa

! ----------- BEGIN ----------

    STRING  = 'DUMMY'

    CALL FETCH_SNDATA_WRAPPER("SIM_MODEL_NAME",   &  ! e.g., SALT2, BYOSED
          ONE, STRING, DARRAY, OPT)
    SIMNAME_MODEL = STRING(1:MXCHAR_MODELNAME)

    CALL FETCH_SNDATA_WRAPPER("SIM_TYPE_NAME",   &  ! e.g., Ia, Ibc
          ONE, STRING, DARRAY, OPT)
    SIMNAME_TYPE = STRING(1:12)

    CALL FETCH_SNDATA_WRAPPER("SIM_MODEL_INDEX",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_MODEL_INDEX = INT(DARRAY(1))

      ! set TRUE_SNIa and NONIa  logicals (Apr 2024)
    LSIM_TRUE_SNIa = .FALSE.
    IF ( ISMODEL_SNIa(SIM_MODEL_INDEX) > 0 ) LSIM_TRUE_SNIa=.TRUE.
    LSIM_TRUE_NONIa = .NOT. LSIM_TRUE_SNIa

    CALL FETCH_SNDATA_WRAPPER("SIM_TYPE_INDEX",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_GENTYPE = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_TEMPLATE_INDEX",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_TEMPLATE_INDEX = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_SUBSAMPLE_INDEX",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_SUBSAMPLE_INDEX = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_LIBID",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_LIBID = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_NGEN_LIBID",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_NGEN_LIBID = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_NOBS_UNDEFINED",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_NOBS_UNDEFINED = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_SEARCHEFF_MASK",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_SEARCHEFF_MASK = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_REDSHIFT_HELIO",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_REDSHIFT_HELIO = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_REDSHIFT_CMB",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_REDSHIFT_CMB = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_REDSHIFT_HOST",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_REDSHIFT_HOST = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_REDSHIFT_HOST_MATCH",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_REDSHIFT_HOST_MATCH = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_REDSHIFT_FLAG",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_REDSHIFT_FLAG = INT(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_VPEC",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_VPEC = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_HOSTLIB_GALID",  & 
          ONE, STRING, DARRAY, OPT)
    DSIM_HOSTLIB_GALID = DARRAY(1)
    SIM_HOSTLIB_GALID  = DSIM_HOSTLIB_GALID

    do ipar = 1, NPAR_SIM_HOSTLIB
       CALL FETCH_SNDATA_WRAPPER(SIM_HOSTLIB_KEYWORD(ipar),  & 
              ONE, STRING, DARRAY, OPT)
       SIM_HOSTLIB_PARVAL(ipar) = SNGL(DARRAY(1))
    enddo

    CALL FETCH_SNDATA_WRAPPER("SIM_DLMU",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_DLMAG = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_LENSDMU",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_LENSDMU = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_MUSHIFT",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_MUSHIFT = SNGL(DARRAY(1))


    CALL FETCH_SNDATA_WRAPPER("SIM_MWEBV",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_MWEBV = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_PEAKMJD",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_PEAKMJD = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_MJD_EXPLODE",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_MJD_EXPLODE = SNGL(DARRAY(1))

! - - - - model params - - - - -

    CALL FETCH_SNDATA_WRAPPER("SIM_MAGSMEAR_COH",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_MAGSMEAR_COH = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_WGT_POPULATION",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_WGT_POPULATION = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2x0",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_SALT2x0 = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2mB",  & 
          ONE, STRING, DARRAY, OPT)
    SIM_SALT2mB = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2c", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_COLORPAR     = SNGL(DARRAY(1))
       SIMNAME_COLORPAR = 'SIM_c'
       SIMNAME_COLORLAW = 'SIM_beta'
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2alpha",ONE,STRING,DARRAY,OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPELAW     = SNGL(DARRAY(1))
       SIMNAME_SHAPELAW = 'SIM_alpha'
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2beta",ONE,STRING,DARRAY,OPT)
    if ( DARRAY(1) > -8.0 )  SIM_COLORLAW = SNGL(DARRAY(1))

    CALL FETCH_SNDATA_WRAPPER("SIM_AV", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_AV = SNGL(DARRAY(1))
! copy SIM_AV as color param if NOT SALT2 model; otherwise leave
! SIM_AV as separate from SALT2c.
      if ( SIM_MODEL_INDEX .NE. MODEL_SALT2 ) then
         SIM_COLORPAR     =  SIM_AV
         SIMNAME_COLORPAR = 'SIM_AV'  ! Jul 2016: SIMAV->SIM_AV
         SIMNAME_COLORLAW = 'SIM_RV'  ! Jul 2016: SIMRV->SIM_RV
      endif

    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_RV", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_RV = SNGL(DARRAY(1))
       if ( SIM_MODEL_INDEX .NE. MODEL_SALT2 ) then
          SIM_COLORLAW = SIM_RV
       endif
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2gammaDM",  & 
           ONE, STRING, DARRAY, OPT)
    SIM_SALT2gammaDM = SNGL(DARRAY(1))

! - - -  shape params - - - -

    CALL FETCH_SNDATA_WRAPPER("SIM_STRETCH", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPEPAR     = SNGL(DARRAY(1))
       SIMNAME_SHAPEPAR = 'SIM_STRETCH'
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_DELTA", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPEPAR     = SNGL(DARRAY(1))
       SIMNAME_SHAPEPAR = 'SIM_DELTA'
    ENDIF

!   Jun 2023 - check for THETA for BayeSN model
    CALL FETCH_SNDATA_WRAPPER("SIM_THETA", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPEPAR     = SNGL(DARRAY(1))
       SIMNAME_SHAPEPAR = 'SIM_THETA'
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_DM15", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPEPAR     = SNGL(DARRAY(1))
       SIMNAME_SHAPEPAR = 'SIM_DM15'
    ENDIF

    CALL FETCH_SNDATA_WRAPPER("SIM_SALT2x1", ONE, STRING, DARRAY, OPT)
    IF ( DARRAY(1) > -8.0 ) THEN
       SIM_SHAPEPAR     = SNGL(DARRAY(1))
       SIMNAME_SHAPEPAR = 'SIM_x1'
    ENDIF

! - - - -  filter dependent keys  - - - - -

    CALL RDHEAD_FILTERLOOP("SIM_PEAKMAG",  SIM_PEAKMAG, OPT)
    CALL RDHEAD_FILTERLOOP("SIM_EXPOSURE", SIM_EXPOSURE_TIME, OPT)
    CALL RDHEAD_FILTERLOOP("SIM_TEMPLATEMAG", SIM_TEMPLATEMAG, OPT) ! 9.20.2021


! - - -  model params: SIMSED, LCLIB, etc   - - - -

    DO ipar = 1, NPAR_SIMSED
       KEY  = SIMSED_KEYWORD(ipar)
       CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
       SIMSED_PARVAL(ipar) = SNGL(DARRAY(1))
    ENDDO


    DO ipar = 1, NPAR_PySEDMODEL
       KEY  = PySEDMODEL_KEYWORD(ipar)
       CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
       PySEDMODEL_PARVAL(ipar) = SNGL(DARRAY(1))
    ENDDO

    DO ipar = 1, NPAR_LCLIB
       KEY = LCLIB_KEYWORD(ipar)
       CALL FETCH_SNDATA_WRAPPER(KEY, ONE, STRING, DARRAY, OPT)
       LCLIB_PARVAL(ipar) = SNGL(DARRAY(1))
    ENDDO

    RETURN
  END SUBROUTINE RDHEAD_SIM_SNANA


! ===========================================
    SUBROUTINE RDOBS_DRIVER()

! Created Feb 11 2021
! Read event (FITS or TEXT format) and transfer epoch-dependent
! information from SNDATA C-struct to fortran variables.
! 
! Jun 7 2021: abort on undefined filter


! local var

    USE SNDATCOM
    USE SNLCINP_NML
    USE FITSCOM
    USE FILTCOM

    IMPLICIT NONE

    REAL*8, PARAMETER :: MJD_SAFETEY = 10.0  ! extra margin for MJD window

    INTEGER OPT, NOBS_STORE, o, LENTMP, IFILT, IFILT_OBS
    REAL*8 DARRAY(MXEPOCH), MJD_WINDOW(2), MJD
    REAL   z, z1, TOBS_MIN, TOBS_MAX
    CHARACTER STRING*12, FNAM*14, CFILT*2

    INTEGER  SELECT_MJD_SNDATA
    EXTERNAL SELECT_MJD_SNDATA

! ---------- BEGIN --------

    FNAM = 'RDOBS_DRIVER'

! set MJD window for copy SNDATA obs
    MJD_WINDOW(1) = 0.0;   MJD_WINDOW(2) = 1.0E6  ! default

! check to refine window around approx peakmjd
    if ( SNLC_SEARCH_PEAKMJD > 40000.0 ) then
       z  = SNLC_REDSHIFT
       z1 = 1.0 + z ;
       if ( z < 0.0 ) z1 = 1.0 + CUTWIN_REDSHIFT(2)
       TOBS_MIN = CUTWIN_TREST(1)*z1  ! note this is usually negative
       TOBS_MAX = CUTWIN_TREST(2)*z1  ! usually positive
       MJD_WINDOW(1) = SNLC_SEARCH_PEAKMJD + TOBS_MIN - MJD_SAFETEY
       MJD_WINDOW(2) = SNLC_SEARCH_PEAKMJD + TOBS_MAX + MJD_SAFETEY
    endif

    NOBS_STORE = SELECT_MJD_SNDATA(MJD_WINDOW)  ! C function

    if ( NOBS_STORE > MXEPOCH ) then
       write(C1ERR,660) NOBS_STORE, MXEPOCH
660      format('NOBS_STORE = ', I6, ' exceeds MXEPOCH=', I6)
       C2ERR = 'Need to increase MXEPOCH in snana.F90'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    ISNLC_NEWMJD_STORE = NOBS_STORE  ! store in global
    ISNLC_NEPOCH_STORE = NOBS_STORE

! - - - - - - - - - - -
! fetch SNDATA values and load fortran globals

    OPT = OPTMASK_SNDATA_OBS
!     &    + OPTMASK_SNDATA_DUMP

    DARRAY(1) = -999.0 ;      STRING = ''

    CALL FETCH_SNDATA_WRAPPER("MJD",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC8_MJD(o) = DARRAY(o);  ENDDO

    CALL FETCH_SNDATA_WRAPPER("BAND",  & 
            NOBS_STORE, STRFITS, DARRAY, OPT)
    CALL UNPACK_SNFITSIO_STR(NOBS_STORE, "FLT", STRFITS)

    CALL FETCH_SNDATA_WRAPPER("DETNUM",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; ISNLC_DETNUM(o)=int(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("IMGNUM",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; ISNLC_IMGNUM(o)=int(DARRAY(o)); ENDDO


    CALL FETCH_SNDATA_WRAPPER("FIELD",  & 
            NOBS_STORE, STRFITS, DARRAY, OPT)
    CALL UNPACK_SNFITSIO_STR(NOBS_STORE, "FIELD", STRFITS)

    CALL FETCH_SNDATA_WRAPPER("PHOTFLAG",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; ISNLC_PHOTFLAG(o)=int(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("PHOTPROB",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_PHOTPROB(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("FLUXCAL",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_FLUXCAL(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("FLUXCALERR",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_FLUXCAL_ERRTOT(o)=SNGL(DARRAY(o)); ENDDO

! - - - -
    UNIT_PSF_NEA = .FALSE.
    CALL FETCH_SNDATA_WRAPPER("PSF_NEA",  & 
           NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE
       SNLC_PSF_NEA(o)=SNGL(DARRAY(o))
       if ( SNLC_PSF_NEA(o) > 0.0 ) UNIT_PSF_NEA = .true.
    ENDDO

    IF ( .NOT. UNIT_PSF_NEA ) THEN
       CALL FETCH_SNDATA_WRAPPER("PSF_SIG1",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SNLC_PSF_SIG1(o)=SNGL(DARRAY(o)); ENDDO

       CALL FETCH_SNDATA_WRAPPER("PSF_SIG2",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SNLC_PSF_SIG2(o)=SNGL(DARRAY(o)); ENDDO

       CALL FETCH_SNDATA_WRAPPER("PSF_RATIO",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SNLC_PSF_RATIO(o)=SNGL(DARRAY(o)); ENDDO
    ENDIF
! - - - - -

    CALL FETCH_SNDATA_WRAPPER("SKY_SIG",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_SKYSIG(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("SKY_SIG_T",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_SKYSIG_T(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("ZEROPT",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_ZEROPT(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("ZEROPT_ERR",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_ZEROPT_ERR(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("TEXPOSE",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_TEXPOSE(o)=SNGL(DARRAY(o)); ENDDO

    CALL FETCH_SNDATA_WRAPPER("GAIN",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SNLC_GAIN(o)=SNGL(DARRAY(o)); ENDDO

    IF ( SNLC_NXPIX > 0.0 ) THEN
       CALL FETCH_SNDATA_WRAPPER("XPIX",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SNLC_XPIX(o)=SNGL(DARRAY(o)); ENDDO

       CALL FETCH_SNDATA_WRAPPER("YPIX",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SNLC_YPIX(o)=SNGL(DARRAY(o)); ENDDO
    ENDIF

! read optional variables for Atmos/DCR (Jul 2023)
    FOUND_ATMOS = .FALSE.
    CALL FETCH_SNDATA_WRAPPER("AIRMASS",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE
        SNLC_AIRMASS(o)=SNGL(DARRAY(o))
        if ( SNLC_AIRMASS(o) > .9 ) FOUND_ATMOS = .true.
    ENDDO

    IF ( FOUND_ATMOS ) THEN
      CALL FETCH_SNDATA_WRAPPER("dRA",  & 
              NOBS_STORE, STRING, DARRAY, OPT)
      DO o=1,NOBS_STORE; SNLC_dRA(o)=SNGL(DARRAY(o)); ENDDO

      CALL FETCH_SNDATA_WRAPPER("dDEC",  & 
                NOBS_STORE, STRING, DARRAY, OPT)
      DO o=1,NOBS_STORE; SNLC_dDEC(o)=SNGL(DARRAY(o)); ENDDO
    ENDIF


! ----------------------------------------------
! ----------- SIM QUANTITIES BELOW -------------
! ----------------------------------------------

    if ( .NOT. ISJOB_SIM ) GOTO 800

! read SIM_MAGOBS for SNANA sim or FAKES ...

    CALL FETCH_SNDATA_WRAPPER("SIM_MAGOBS",  & 
            NOBS_STORE, STRING, DARRAY, OPT)
    DO o=1,NOBS_STORE; SIM_EPMAGOBS(o)=SNGL(DARRAY(o)); ENDDO

! the rest is for SNANA sim only ...

    if ( .NOT. LSIM_SNANA ) GOTO 800

    LENTMP = INDEX(SIMNAME_SNRMON,' ') - 1
    IF ( LENTMP > 2 ) THEN
       CALL FETCH_SNDATA_WRAPPER(SIMNAME_SNRMON,  & 
              NOBS_STORE, STRING, DARRAY, OPT)
       DO o=1,NOBS_STORE; SIM_EPSNRMON(o)=SNGL(DARRAY(o)); ENDDO
    ENDIF

! ---------------------------------
 800  CONTINUE

! compute a few misc variables related to input data;
! these are essentially change of units.

    DO o = 1, NOBS_STORE
      CALL SET_EPVAR_MISC(o)   ! set SNLC_ZP (ADU or Npe) and PSF_FWHM

      ISNLC_EPOCH_RANGE_NEWMJD(1,o) = o
      ISNLC_EPOCH_RANGE_NEWMJD(2,o) = o

      IFILT_OBS = ISNLC_IFILT_OBS(o)
      IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

      IF ( IFILT < 0 ) THEN
         MJD = SNLC8_MJD(o)
         cfilt  = filtdef_string(ifilt_obs:ifilt_obs)
         write(c1err,61) IFILT_OBS, CFILT, SNLC_CCID(1:ISNLC_LENCCID)
61         format('Invalid IFILT_OBS=',I3, '(',A1,') for CID=',A)
         write(c2err,62) MJD
62         format('Check FILTERS key and MJD=', F9.3)
         CALL MADABORT(FNAM, c1err,c2err)
      ENDIF

      ISNLC_NEPOCH_FILT(ifilt) =  & 
             ISNLC_NEPOCH_FILT(ifilt) + 1

!        print*,' xxx o, ifilt, ifilt_obs = ',
!     &        o, ifilt, ifilt_obs
    ENDDO

! Since SET_PEAKMJD has not been called [yet],
! apply all epoch cuts EXCEPT for Trest.

    CALL SELECT_EPOCH_DRIVER(1) ! 1 -> all obs/epoch cuts EXCEPT for Trest

    RETURN
  END SUBROUTINE RDOBS_DRIVER


! ===========================================
    SUBROUTINE RDSPEC_DRIVER(ISPEC)

! 
! Created Feb 17 2021
! Read spectra (FITS or TEXT format) and transfer
! information from GENSPEC C-struct to fortran variables.
! 
! ISPEC = 0 -> read header info
! ISPEC > 0 -> read wave-dependent arrays for ISPEC


    USE SNDATCOM
    USE SPECCOM

    IMPLICIT NONE

    INTEGER ISPEC  ! (I)

    INTEGER    OPT, NBLAM, I
    REAL*8     DVAL(10)
    CHARACTER  FNAM*14

! ------------ BEGIN ----------

    FNAM = 'RDSPEC_DRIVER'
    OPT = 0
!     &    + OPTMASK_SNDATA_DUMP

    IF ( ISPEC == 0 ) THEN
       CALL FETCH_GENSPEC_WRAPPER("NSPECTRA", -9, DVAL, OPT)
       NSPECTRUM = INT(DVAL(1))

       IF ( NSPECTRUM .LE. 0 ) RETURN

       IF ( NSPECTRUM > MXSPECTRUM .or. NSPECTRUM < 0 ) THEN
          write(c1err, '(A,I3)' ) 'Invalid NSPECTRUM = ', NSPECTRUM
          write(c2err, '(A,I3)' ) 'MXSPECTRUM = ', MXSPECTRUM
          CALL MADABORT(FNAM, c1err,c2err)
       ENDIF


! load up all the spec-header info now; the big LAM-dependent
! arrays are loaded later with ISPEC > 0.

       DO I = 1, NSPECTRUM
          CALL FETCH_GENSPEC_WRAPPER("NBLAM", I, DVAL, OPT)
          NLAMBIN_SPECTRUM(I) = INT(DVAL(1))
          NBLAM = NLAMBIN_SPECTRUM(I)

          if ( NBLAM > MXLAM_SPECTRUM .or. NBLAM < 0 ) then
             write(c1err,'(A,I8)') 'Invalid NLAMBIN = ', NBLAM
             write(c2err,'(A,I8)') 'MXLAM_SPECTRUM=',MXLAM_SPECTRUM
             CALL MADABORT(FNAM, c1err,c2err)
          endif

          CALL FETCH_GENSPEC_WRAPPER("ID", I, DVAL, OPT)
          ID_SPECTRUM(I)    = INT(DVAL(1))

          CALL FETCH_GENSPEC_WRAPPER("MJD", I, DVAL, OPT)
          MJD_SPECTRUM(I)    = DVAL(1)
          TOBS_SPECTRUM(I)   = DVAL(1) - SNLC_SEARCH_PEAKMJD

          CALL FETCH_GENSPEC_WRAPPER("TEXPOSE", I,  & 
                 TEXPOSE_SPECTRUM(I), OPT)
       ENDDO

       RETURN
    ENDIF

! - - - - - - -  -
! retreive SPECTRUM for ISPEC.
! wave-dependent arrays
    CALL FETCH_GENSPEC_WRAPPER("LAMMIN", ISPEC,  & 
           LAMMIN_SPECTRUM, OPT)
    CALL FETCH_GENSPEC_WRAPPER("LAMMAX", ISPEC,  & 
           LAMMAX_SPECTRUM, OPT)

    CALL FETCH_GENSPEC_WRAPPER("FLAM", ISPEC,  & 
           FLAM_SPECTRUM, OPT)
    CALL FETCH_GENSPEC_WRAPPER("FLAMERR", ISPEC,  & 
           FLAMERR_SPECTRUM, OPT)

    IF (  LSIM_SNANA ) THEN
       CALL FETCH_GENSPEC_WRAPPER("SIM_FLAM", ISPEC,  & 
              SIM_FLAM_SPECTRUM, OPT)
    ENDIF

    RETURN
  END SUBROUTINE RDSPEC_DRIVER


! ==============================================
    LOGICAL FUNCTION LDONE_CIDLIST()

! Jun 2011:
! Returns TRUE if all CID in the list have been processed;
! allows parsing function to stop early.
! 
! May 18, 2012: return FALSE if interp-option is set
! Jun 25, 2019: check MXEVT_PROCESS
! Oct 12, 2021: check MXEVT_CUTS
! -----------

    USE SNDATCOM
    USE SNLCINP_NML
    USE INTERPCM

    IMPLICIT NONE

    LOGICAL LNCID

! ------------- BEGIN --------------

    LDONE_CIDLIST = .FALSE.  ! default

    if ( N_SNLC_PROC .GE. MXEVT_PROCESS ) then
       LDONE_CIDLIST = .TRUE.
       RETURN
    endif

    if ( N_SNLC_CUTS .GE. MXEVT_CUTS ) then
       LDONE_CIDLIST = .TRUE.
       RETURN
    endif

! if CUTWIN_CID is set then return FALSE
    if ( cutwin_cid(1) .GT. 0 ) RETURN
    if ( cutwin_cid(2) .GT. 0 ) RETURN

    if ( SNCID_LIST_FILE .NE. ' ' ) RETURN ! Nov 2013

! return FALSE if interp-option is set but not finished

    if ( N_INTERP_MJDLIST .GT. 0 .and.  & 
           N_INTERP_MJDLIST_DONE .LT. N_INTERP_MJDLIST ) THEN
         LDONE_CIDLIST = .FALSE.
         RETURN
    endif

    LNCID = NACCEPT_CID .EQ. (NCID_LIST + NCCID_LIST)
    IF ( LNCID ) THEN
      LDONE_CIDLIST = .TRUE.
    ENDIF

    RETURN
  END FUNCTION LDONE_CIDLIST


! =======================
    SUBROUTINE PARSE_PRIVATE_KEYWORD(CWORD, IVAR_PRIVATE)

! For input CWORD = 'PRIVATE($VARNAME)', strip out the $VARNAME
! and load PRIVATE_VARNAME(IVAR_PRIVATE) = $VARNAME
! 
! Apr 27 2021: check veto option with VARNAME != VALUE
! 

! local args
    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM

    IMPLICIT NONE

! input args
    CHARACTER CWORD*60     ! (I) text string to parse
    INTEGER IVAR_PRIVATE   ! (I) PRIVATE parameter index

! local args
    INTEGER J1, J2, LK, IVAR, LEN_NAME, NMATCH, icut
    INTEGER NWD, LWD1, LWD2, LWD3, MSKOPT
    CHARACTER VARNAME*60, STR*(MXCHAR_CUTNAME), CUTSTRING*30
    CHARACTER CWD1*60, CWD2*60, CWD3*60, FNAM*30

    LOGICAL LCUTWIN
    REAL*8 CUTWIN(2)

! function
    INTEGER NMATCH_PRIVATE

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! --------------- BEGIN --------------

    FNAM = 'PARSE_PRIVATE_KEYWORD'

! store keyword WITHOUT the ':' on the end.
    LK = INDEX( CWORD , ':' ) - 1
    IF ( LK .LT. 0 ) LK = INDEX(CWORD,' ' ) - 1
    PRIVATE_KEYWORD(IVAR_PRIVATE) = CWORD(1:LK)

! now extract the parameter name from inside the ().
! For example, if KEYWORD = 'PRIVATE(BLA):'
! then the VARNAME is 'BLA'.

    J1 = INDEX( CWORD , '(' )
    J2 = INDEX( CWORD , ')' )

    IF ( J1 .EQ. 0 .or. J2 .EQ. 0 ) THEN
       c1err = 'Invalid PRIVATE_KEYWORD = ' // CWORD(1:LK)
       c2err = 'Expecting to find key-string  PRIVATE(PARNAME)'
       CALL MADABORT("PARSE_PRIVATE_KEYWORD", c1err,c2err)
    ENDIF


    VARNAME = CWORD(j1+1:j2-1)
    LEN_NAME = INDEX(VARNAME // ' ', ' ' ) - 1

! load parameter name into global variable

    PRIVATE_VARNAME(IVAR_PRIVATE) = VARNAME

! ----------------------------------------------------
! make sure that this NAME has not already been used

    NMATCH = NMATCH_PRIVATE(VARNAME,IVAR)

    IF ( NMATCH > 1 ) then
       c1err = 'found duplicate PRIVATE variable: '  & 
                      // VARNAME(1:LEN_NAME)
       c2err = 'Check PRIVATE(XXX) keys'
       CALL MADABORT("PARSE_PRIVATE_KEYWORD", c1err,c2err)
    ENDIF

! -------- check for cut on private var (Nov 4 2014) ---------
! Allow two syntaxes:
!  VARNAME CUTMIN CUTMAX    or
!  VARNAME != VETO_VALUE

    LCUTWIN = .FALSE.
    CUTSTRING = ''
    DO 400 icut = 1, MXCUT_PRIVATE

       STR = PRIVATE_CUTWIN_STRING(icut)  ! 'VARNAME CUTMIN CUTMAX'
       if ( STR .EQ. '' ) GOTO 400

!        skip if already found
       CUTWIN(2) = PRIVATE_CUTWIN(2,IVAR_PRIVATE)
       if ( CUTWIN(2) < .99*CUTVAL_OPEN ) GOTO 400

! break up string into pieces
       MSKOPT = MSKOPT_PARSE_WORDS_STRING
       NWD = STORE_PARSE_WORDS(MSKOPT,STR//' '//char(0),  & 
              FNAM//char(0), 60, 30)
       call get_PARSE_WORD_fortran(1,CWD1,LWD1)

       if ( VARNAME .EQ. CWD1 ) then
          NCUT_PRIVATE = NCUT_PRIVATE + 1
          LCUTWIN = .TRUE.

          call get_PARSE_WORD_fortran(2,CWD2,LWD2)
          call get_PARSE_WORD_fortran(3,CWD3,LWD3)

          if ( CWD2(1:2) .EQ. '!=' ) then
             read(CWD3,*) CUTWIN(1)
             CUTWIN(2) = CUTWIN(1)
             USE_PRIVATE_CUTWIN(IVAR_PRIVATE) = -1 ! veto flag
             write(CUTSTRING,446) CUTWIN(1)
          else
             read(CWD2,*) CUTWIN(1)
             read(CWD3,*) CUTWIN(2)
             USE_PRIVATE_CUTWIN(IVAR_PRIVATE) = 1 ! normal cutwin
             write(CUTSTRING,444) CUTWIN
          endif

          PRIVATE_CUTWIN(1,IVAR_PRIVATE) = CUTWIN(1)
          PRIVATE_CUTWIN(2,IVAR_PRIVATE) = CUTWIN(2)

 444        format('CUT: ',G9.3,' to ', G9.3)
 446        format('VETO: ', G9.3)
       endif
 400  CONTINUE

! ----------------------------
! print one-line summary

    IF ( .NOT. DO_GETINFO ) THEN
      write(6,500) CWORD(1:LK), CUTSTRING
 500    format(T5,'Found ', A, T50,A)
      call flush(6)
    ENDIF

    RETURN
  END SUBROUTINE PARSE_PRIVATE_KEYWORD

! ===============================
    INTEGER FUNCTION NMATCH_PRIVATE(NAME,IVAR_PRIVATE)

! 
! Return number of NAME matches ; check both
! PRIVATE_VARNAME and PRIVATE_KEYWORD.
! Also returns IVAR_PRIVATE that is valid if NMATCH=1


    USE SNPAR
    USE PRIVCOM

    IMPLICIT NONE

    CHARACTER NAME*(*)
    INTEGER   IVAR_PRIVATE

    INTEGER  & 
         IVAR  & 
        ,LEN_TMP, LEN_NAME  & 
        ,NMATCH

    CHARACTER TMPNAME*60

! --------------- BEGIN --------------

    IVAR_PRIVATE   = -9
    NMATCH_PRIVATE =  0
    NMATCH         =  0
    LEN_NAME       = INDEX(NAME//' ' , ' ') - 1

    DO 30 IVAR = 1, NVAR_PRIVATE
      TMPNAME  = PRIVATE_VARNAME(IVAR)
      LEN_TMP  = INDEX(TMPNAME,' ') - 1
      IF ( LEN_TMP .NE. LEN_NAME ) goto 30

      if ( NAME(1:LEN_NAME) .EQ. TMPNAME(1:LEN_NAME) ) then
         IVAR_PRIVATE = IVAR
         NMATCH = NMATCH + 1
      endif
30    CONTINUE

! check the full name 'PRIVATE($VARNAME)'

    DO 32 IVAR = 1, NVAR_PRIVATE
      TMPNAME  = PRIVATE_KEYWORD(IVAR)
      LEN_TMP  = INDEX(TMPNAME,' ') - 1
      IF ( LEN_TMP .NE. LEN_NAME ) goto 32

      if ( NAME(1:LEN_NAME) .EQ. TMPNAME(1:LEN_NAME) ) then
         IVAR_PRIVATE = IVAR
         NMATCH = NMATCH + 1
      endif

32    CONTINUE


    NMATCH_PRIVATE = NMATCH

    RETURN
  END FUNCTION NMATCH_PRIVATE


! END_DATA_IO

! ==========================================

    SUBROUTINE PROCESS_DATA_VERSION(IVERS)

! Created Feb 16 2021 as part of I/O refactor.
! Called from MAIN to read & ananlyze events.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IVERS  ! (I) version index to process

! --------------- BEGIN -------------

    CALL GETINFO_PHOTOMETRY(ivers)

    CALL READ_EPIGNORE_FILE(ivers)

    CALL EXEC_READ_DATA(ivers,OPTMASK_SNDATA_ALL)

    RETURN
  END SUBROUTINE PROCESS_DATA_VERSION

! ==========================================
    SUBROUTINE PARSE_SNANA_ARGS()
! 
! Nov 08, 2012
! Parse command-line arguments and store in LINE_ARGS array.
! (code moved from main)
! 
! May 2019: check option to call DUMP_SNANA_VERSION
! Dec 2020: abort on ARG string length too long.
! Aug 2020: of no args, call function to print help menu
! Dec 2022: get SNANA_VERSION here to print it after full command


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, LL, LL_SNANA, istat
    CHARACTER ARG1*(MXCHAR_ARG), ARG2*(MXCHAR_ARG)
    CHARACTER FNAM*20, COMMAND*2000, CODE_FILENAME*40
    CHARACTER CMD_HELP*100, CMD_SNANA_HELP*100

    EXTERNAL GET_SNANA_VERSION

! ------------------- BEGIN ---------------

    FNAM = 'PARSE_SNANA_ARGS'

    NLINE_ARGS = IARGC()  ! intrinsic fortran function

! get SNANA version from C function
    call GET_SNANA_VERSION(SNANA_VERSION,40)
    LL_SNANA = index(SNANA_VERSION, char(0) ) - 1
    SNANA_VERSION(LL_SNANA+1:LL_SNANA+2) = ' '

! - - - - - - - - - - - - - - - - - -
! Aug 2022: if there are no arguments, print help menu
    IF ( NLINE_ARGS == 0 ) THEN

!  always print help for &SNLCINP base code inputs
       CMD_SNANA_HELP = 'help_inputs_fortran.py ' // 'snana.F90'
       CALL SYSTEM(CMD_SNANA_HELP)

#if defined(SNANA)
        CODE_FILENAME = ''
#elif defined(SNFIT)
        CODE_FILENAME = 'snlc_fit.F90'
#elif defined(PSNID)
        CODE_FILENAME = 'psnid.F90'
#endif

	 if ( CODE_FILENAME .NE. '' )  then
          CMD_HELP = 'help_inputs_fortran.py ' // CODE_FILENAME
          CALL SYSTEM(CMD_HELP)
       endif
       STOP
    ENDIF
! - - - - - - - - - - - - - - - - - -

    IF ( NLINE_ARGS .GT. MXLINE_ARGS ) THEN
       write(c1err,661) NLINE_ARGS
       write(c2err,662) MXLINE_ARGS
661      format(I3,' command-line arguments exceeds bound of')
662      format('MXLINE_ARGS = ', I4)
       CALL  MADABORT("PARSE_SNANA_ARGS", C1ERR, C2ERR)
    ENDIF


    CALL GET_COMMAND(COMMAND,LL, istat)

! abort on bad string length
    DO i = 1, NLINE_ARGS
      CALL GETARG(i, ARG1)

      if ( ARG1(MXCHAR_ARG:MXCHAR_ARG) .NE. ' ' ) then
         CALL PRINT_PREABORT_BANNER(FNAM//char(0),20)
         print*,'  ARG = ', ARG1
         c1err = 'String length too long for command-line ARG ' //  & 
                   ' (see above)'
         write(c2err, '(A,I6)' ) 'MXCHAR_ARG = ', MXCHAR_ARG
         CALL MADABORT(FNAM, c1err, c2err)
      endif

      USE_LINE_ARGS(i) = .FALSE.
    ENDDO

    CALL GETARG(1, ARG1)
    CALL GETARG(2, ARG2)

! check option to dump SNANA version and SNANA_DIR
! (for SNANA_CodeTests.py)
!    snana.exe --snana_version

    IF ( ARG1 .EQ. '--SNANA_VERSION' .or.  & 
           ARG1 .EQ. '--snana_version'  ) THEN
       CALL DUMP_SNANA_VERSION()
    ELSE
       print*,' Full command:  ', COMMAND(1:LL)
       print*,' '
       print*,' SNANA_VERSION: ', SNANA_VERSION(1:LL_SNANA)
       print*,' '

    ENDIF

! snana.exe GETINFO <VERSION_PHOTOMETRY>
    DO_GETINFO = .FALSE.
    IF ( ARG1 .EQ. 'GETINFO' ) THEN
       CALL DUMP_INFO_VERSION_PHOTOMETRY(ARG2)
    ENDIF

    RETURN
  END SUBROUTINE PARSE_SNANA_ARGS


! ==================================
    SUBROUTINE DUMP_SNANA_VERSION()

! Created May 2019
! Dump SNANA version and SNANA_DIR, then quit.

    INTEGER LL
    CHARACTER SNANA_VERSION*40, SNANA_DIR*200

    call GET_SNANA_VERSION(SNANA_VERSION,LL)
    LL = index(SNANA_VERSION, char(0) ) - 1
    print*,'SNANA_VERSION: ', SNANA_VERSION(1:LL)

    CALL GETENV ( 'SNANA_DIR', SNANA_DIR )
    LL = INDEX(SNANA_DIR,' ') - 1
    print*,'SNANA_DIR: ', SNANA_DIR(1:LL)

    CALL EXIT(0)

    RETURN
  END SUBROUTINE DUMP_SNANA_VERSION

! ======================================
    SUBROUTINE DUMP_INFO_VERSION_PHOTOMETRY(VERSION_ARG)

! Created Aug 29 2019
! dump path and format for VERSION_PHOTOMETRY; then quit.
! 
! Apr 2021: pass VERSION_ARG and dump more info
! Jan 10 2022: write NVAR_PRIVATE


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE PRIVCOM

    IMPLICIT NONE

    CHARACTER VERSION_ARG*(*)  ! VERSION or PATH/VERSION

! local var
    CHARACTER VERSION*(MXCHAR_FILENAME), COMMENT*40
    INTEGER IERR, LEN, JSLASH
! -------------- BEGIN ------------

    DO_GETINFO = .TRUE.

! if VERSION_ARG has full path, split it into PRIVATE_DATA_PATH and
! VERSION_PHOTOMETRY

    LEN     = INDEX(VERSION_ARG,' ') -1
    JSLASH  = INDEX(VERSION_ARG, '/', BACK = .True. )
    IF ( JSLASH > 1 ) THEN
       VERSION_PHOTOMETRY(1) = VERSION_ARG(JSLASH+1:LEN)
       PRIVATE_DATA_PATH     = VERSION_ARG(1:JSLASH-1)
    ELSE
       VERSION_PHOTOMETRY(1) = VERSION_ARG(1:LEN)
    ENDIF

    CALL INIT_SNVAR()
    CALL CHECK_FORMAT()
    CALL INIT_READ_DATA()

! - - - - - - - - - - - -
    print*,' '
    print*,' '

    VERSION = VERSION_PHOTOMETRY(1)
    LEN     = INDEX(VERSION,' ')-1
    write(6,20) 'VERSION:', VERSION(1:LEN), ''

    LEN = INDEX(SNDATA_PATH,' ') - 1
    write(6,20) 'SNDATA_PATH:', SNDATA_PATH(1:LEN), ' '

    LEN = INDEX(SNLIST_FILE,' ') - 1
    write(6,20) 'LIST_FILE:  ', SNLIST_FILE(1:LEN), ' '

    LEN = INDEX(SNREADME_FILE(1),' ') - 1
    write(6,20) 'README_FILE:', SNREADME_FILE(1)(1:LEN), ' '

 20   format(A, 2x, A,  3x, A)
 21   format(A, 2x, I8, 3x, A)
! - - - - - - - - - - - -
    if ( FORMAT_TEXT ) then
       write(6,20) 'FORMAT: ', 'TEXT',  & 
               '# one TEXT file per SN'
    else
       write(6,20) 'FORMAT: ', 'FITS',  & 
               '# FITS-HEADER and FITS-PHOT file'
    endif

    LEN = INDEX(SURVEY_NAME,' ') - 1
    write(6,20) 'SURVEY: ', SURVEY_NAME(1:LEN), ' '

    write(6,20) 'FILTERS:', SURVEY_FILTERS(1:NFILTDEF_SURVEY), ' '

! xxx mark delete Oct 16 2025 xxxxxxx
!      write(6,21) 'NEVENT: ',  N_SNLC_READ(1), ' '
! xxxxxxxxxx

    COMMENT = ''
    IF ( DATATYPE == 'SIM_MAGOBS' ) THEN
       COMMENT='# Fakes overlaid on images'
    ENDIF
    write(6,20) 'DATATYPE: ',  DATATYPE, COMMENT

    write(6,21) 'NVAR_PRIVATE: ',  NVAR_PRIVATE, ' '

    CALL FLUSH(6)
    CALL EXIT(0)

    RETURN
  END SUBROUTINE DUMP_INFO_VERSION_PHOTOMETRY

! ==================================
    SUBROUTINE RDSNNML(IERR)
    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

! 
! Feb 27 2018: check CUTWIN_ZP
! Apr 07 2022; use function VALIDFILENAME
! 
! ---------------

    INTEGER IERR   ! 0=>OK,  else error

! local ars

    integer i, ifilt, LL, ipar, NZPCUT, bit
    character nmlfile_default*40, FNAM*8

! function
    INTEGER PARSE_NML_STRINGLIST
    LOGICAL VALIDFILENAME

! -------------- BEGIN -----------
    IERR = 0

    CALL PRBANNER ( " READ SNLCINP NAMELIST. " )

    FNAM = 'RDSNNML'

! ----------------------------------------
! set namelist defaults

    VERSION_PHOTOMETRY_WILDCARD = 'NULL'
    DO i = 1, MXVERS
      VERSION_PHOTOMETRY(i) = 'NULL'
    ENDDO
    VERSION_REFORMAT_FITS = 'NULL'
    VERSION_REFORMAT_TEXT = 'NULL'

    PRIVATE_DATA_PATH  = ' '
    FILTER_UPDATE_PATH = ' '
    OPT_YAML              = 0
    OPT_REFORMAT_SALT2    = 0
    OPT_REFORMAT_TEXT     = 0
    OPT_REFORMAT_FITS     = 0
    OPT_REFORMAT_SPECTRA  = 0
    REFORMAT_SAVE_BADEPOCHS  = .FALSE.
    REFORMAT_BAND_NAME       = .FALSE.
    REFORMAT_SPECTRA_INCLUDE = .FALSE.
    REFORMAT_SPECTRA_ONLY    = .FALSE.
    REFORMAT_PRIVATE         = .TRUE.
    REFORMAT_SIMTRUTH        = .TRUE.
    REFORMAT_KEYS     = 'NULL'

    REDUCE_STDOUT_BATCH = .FALSE.
    FORCE_STDOUT_BATCH  = .FALSE.

    NONSURVEY_FILTERS = ''
    SNRMAX_FILTERS    = ''  ! blank -> use all filters
    FILTER_REPLACE    = ''
    FILTLIST_LAMSHIFT = ''

    VPEC_ERR_OVERRIDE = 0.0

    JOBSPLIT(1) = 1
    JOBSPLIT(2) = 1

    JOBSPLIT_EXTERNAL(1) = 1
    JOBSPLIT_EXTERNAL(2) = 1

    SIM_PRESCALE   = 1.0
    OPTSIM_LCWIDTH = 0

    MNFIT_PKMJD_LOGFILE = 'MNFIT_PKMJD.LOG'

    SNMJD_LIST_FILE = ' '
    SNMJD_OUT_FILE  = ' '

    CUTWIN_MJD_EXCLUDE(1) = 0.0
    CUTWIN_MJD_EXCLUDE(2) = 0.0

    CUTWIN_TOBS_PREDETECT(1) = -CUTVAL_OPEN
    CUTWIN_TOBS_PREDETECT(2) = +CUTVAL_OPEN

    CUTWIN_SNR_NODETECT(1) = -CUTVAL_OPEN
    CUTWIN_SNR_NODETECT(2) = +CUTVAL_OPEN

    SNCUT_HOST_SBFLUX = ''
    SNCUT_NOBS_TREST  = ''  ! Nov 2025
    SNCUT_SNRMAX = ''
    EPCUT_SNRMIN = ''
    EPCUT_SKYSIG = ''
    EPCUT_PSFSIG = ''

    CUTWIN_NFIELD(1) = 1.0
    CUTWIN_NFIELD(2) = 999.0

    CUTWIN_NFILT_SNRMAX(1) = 0.0
    CUTWIN_NFILT_SNRMAX(2) = 100.

    CUTWIN_NFILT_SNRMAX2(1) = 0.0
    CUTWIN_NFILT_SNRMAX2(2) = 100.

    CUTWIN_NFILT_TRESTMIN(1) =   1.
    CUTWIN_NFILT_TRESTMIN(2) = 100.
    CUTWIN_NFILT_TRESTMAX(1) =   1.
    CUTWIN_NFILT_TRESTMAX(2) = 100.
    CUTWIN_NFILT_TREST2(1)   =   1.
    CUTWIN_NFILT_TREST2(2)   = 100.

    CUTWIN_LAMREST(1)      =  2000.  ! Angstroms
    CUTWIN_LAMREST(2)      = 25000.
    CUTWIN_LAMOBS(1)       =  2000.
    CUTWIN_LAMOBS(2)       = 70000.  ! increase from 25k -> 70k, for JWST (May 2024)

    CUTWIN_NBAND_THRESH(1) = -9.0   !  obsolete as of Nov 15 2019
    CUTWIN_NBAND_THRESH(2) = -9.0

    do ifilt = 1, MXFILT_OBS
       CUTWIN_SNRMIN_FILT(1,ifilt) = -CUTVAL_OPEN
       CUTWIN_SNRMIN_FILT(2,ifilt) = +CUTVAL_OPEN
    enddo

    CUTWIN_TREST_TRUEFLUX(1)  = 999.0 ! user input
    CUTWIN_TREST_TRUEFLUX(2)  = 999.0
    CUTWIN_TREST_TRUEFLUX2(1) = 0.9  ! logical cutwin
    CUTWIN_TREST_TRUEFLUX2(2) = 1.1

! Setting ignore file to 'NONE' => do not read any file.
! Setting INGORE file to 'DEFAULT' is a flag to read
! epochs-to-ignore from $SNDATA_ROOT/lcmerge/[version].IGNORE

    EPOCH_IGNORE_FILE = 'DEFAULT'
    OUT_EPOCH_IGNORE_FILE = ''

    CALIB_FILE    = 'NULL'
    KCOR_FILE     = 'NULL'
    OVERRIDE_RESTLAM_BOUNDARY = ''

    rootfile_out      = ' '
    textfile_prefix   = ' '
    MARZFILE_OUT      = ' '  ! May 2020

    SNCID_LIST_FILE   = ' '
    OPT_SNCID_LIST    = 0  ! default 0 => AND with other cuts
    OPT_VPEC_COR      = 1  ! default is to correct zHD with vpec

    MXEVT_PROCESS      = 999888777
    MXEVT_CUTS         = 999888777
    SIMLIB_OUT     = ''
    SIMLIB_OUTFILE = ''
    SIMLIB_ZPERR_LIST = ''
    OPT_SIMLIB_OUT    = 1
    SIMLIB_OUT_TMINFIX = 99999 ! large value -> not used

    NONLINEARITY_FILE = ''

    DO i = 1, MXLISTNML
       SNTYPE_LIST(i)    = 0
       SNTYPE_IGNORE(i)  = -9  ! note that type=0(no type) can be ignored
       DETNUM_LIST(i)    = 0
       SNCID_LIST(i)     = 0
       SNCID_IGNORE(i)   = 0
       SIM_TEMPLATE_INDEX_LIST(i) = 0

       SNCCID_LIST(i)    = ''
       SNCCID_IGNORE(i)  = ''
       SNCCID_PLOT(i)    = ''

       SNFIELD_LIST(i)     = ''
       IDFIELD_LIST(i)     = -9
    ENDDO

! set default field to ALL
    SNFIELD_LIST(1)   = 'ALL'

    NSNTYPE_LIST = 0
    NDETNUM_LIST = 0
    NSNTYPE_IGNORE = 0

    NIDFIELD_LIST = 0
    NCCID_LIST = 0
    NCID_LIST  = 0

    SNCID_IGNORE_FILE = ''

! set photometry-flag bits to reject;
! default is to reject nothing (accept everything)

    PHOTFLAG_MSKREJ(1) = -1  ! negative -> not set
    PHOTFLAG_MSKREJ(2) = 0
    PHOTFLAG_MSKREJ(3) = 0
    PHOTFLAG_MSKREJ(4) = 0
    PHOTFLAG_MSKREJ(5) = 0

    do i = 1, MXBIT_PHOTFLAG
       PHOTFLAG_BITLIST_REJECT(i) = -9 ! beware that 0 is valid bit
    enddo

    SNCID_IGNORE(1)    = 0

    MAGOBS_SHIFT_PRIMARY  = ''
    MAGOBS_SHIFT_ZP       = ''
    MAGREST_SHIFT_PRIMARY = ''
    MAGREST_SHIFT_ZP      = ''
    FILTER_LAMSHIFT       = ''

    DO i = 1, 3
      MAGOBS_SHIFT_PRIMARY_PARAMS(i) = 0.0
      MAGOBS_SHIFT_ZP_PARAMS(i)      = 0.0
    ENDDO

    FUDGE_FLUXCAL_OFFSET = ''
    FUDGE_FLUXCAL_ERROR  = ''
    FUDGE_FLUXCAL_ERRPIX = ''
    FUDGE_FLUXERR_SCALE  = ''
    FUDGE_MAG_ERROR      = ''
    FUDGE_MAG_COVERR     = ''
    FUDGE_HOSTNOISE_FILE = ''
    SIM_FUDGE_MAG_ERROR  = ''
    FLUXERRMODEL_FILE    = ''
    SIM_FLUXERRMODEL_FILE = ''
    FLUXERRMODEL_OPTMASK = 0
    MAGCOR_FILE          = ''
    SIM_MAGCOR_FILE      = ''

    MWEBV_SCALE    = 1.0
    MWEBV_SHIFT    = 0.0
    MWEBV_FORCE    = -9.0
    OPT_MWCOLORLAW = -999   ! use default if not specified in &SNLCINP
    OPT_MWEBV      = -9   ! idem
    RV_MWCOLORLAW  = -9.0 ! idem

    DO i = 1, 10
       PARLIST_MWCOLORLAW(i) = 0.0
    ENDDO

    HOSTGAL_ZPHOT_SHIFT  = 0.0
    HOSTGAL_ZSPEC_SHIFT  = 0.0
    HOSTGAL_PHOTOZ_SHIFT = LEGACY_INIT_VAL
    HOSTGAL_SPECZ_SHIFT  = LEGACY_INIT_VAL
    REDSHIFT_FINAL_SHIFT = 0.0
    FLUXERRCALC_ZPTERR  = -9.9

    NFIT_ITERATION     = 0
    MINUIT_PRINT_LEVEL = -1  ! default is no MINUIT printing
    INTERP_OPT     = INTERP_LINEAR
    USE_MINOS      = .FALSE.  ! change from T to F, Jan 27 2017

    MXLC_FIT         = 999888777
    MXLC_PLOT        = 5    ! 100->5  (Apr 19 2022)
    NCCID_PLOT       = 0
    DTOBS_MODEL_PLOT = 2.0  ! 2 day binning
    MJDPERIOD_PLOT   = 0.0
    MJDSHIFT_PLOT    = 0.0

    OPT_SETPKMJD      = 0
    SNRCUT_SETPKMJD   =  5.0  ! only for FLUXMAX2 option (Fmax-clump)
    MJDWIN_SETPKMJD   = 60.0  ! only for FLUXMAX2 option (Fmax-clump)
    SHIFT_SETPKMJD    =  0.0  ! for systematic test only

    DEBUG_FLAG          =  0
    LDMP_SNFAIL         = .FALSE.
    LDMP_SNANA_VERSION  = .FALSE.
    LDMP_SATURATE       = .FALSE.

    LSIM_SEARCH_SPEC    = .FALSE.
    LSIM_SEARCH_ZHOST   = .FALSE.

    USESIM_SNIA      = .TRUE.
    USESIM_NONIA     = .TRUE.
    USESIM_TRUEFLUX  = .FALSE.
    USESIM_REDSHIFT  = .FALSE.

    LPROB_TRUEFLUX     = .FALSE.
    RESTORE_WRONG_VPEC = .FALSE.    ! Nov 3 2020
    RESTORE_OVERRIDE_ZBUG = .FALSE. ! Dec 12 2021
    RESTORE_MWEBV_ERR_BUG = .FALSE. ! Jul 2022
    RESTORE_DES5YR        = .FALSE. ! May 28 2025

    REQUIRE_DOCANA     =  0       ! use integer to match sim usage

    USE_HOSTGAL_PHOTOZ = .FALSE.
    USE_SNHOST_ZPHOT   = .FALSE.
    USE_MWCOR          = .FALSE.  ! apply MWCOR to fit-model

! by default do NOT check for multiseason activity (Oct 2014)
    MULTISEASON_TGAP               = 1.0E9
    MULTISEASON_NREJECT_OUTLIER    = 0
    MULTISEASON_CHI2RED_ACTIVE     = 5.0

    PHOTFLAG_DETECT   = 0
    PHOTFLAG_TRIGGER  = 0

    CUTWIN_NSEASON_ACTIVE(1) = -99.
    CUTWIN_NSEASON_ACTIVE(2) = +99.

    CUTWIN_CUTFLAG_REQEP(1) = 0.99
    CUTWIN_CUTFLAG_REQEP(2) = 1.01

! hard-wire cut on private cut-flag (Nov 4 2014)
    CUTWIN_CUTFLAG_PRIVATE(1) = 0.99
    CUTWIN_CUTFLAG_PRIVATE(2) = 1.01

    CUTWIN_CUTFLAG_SIMVAR(1) = 0.99
    CUTWIN_CUTFLAG_SIMVAR(2) = 1.01

! set reference cosmo params to defaults.
    H0_REF(1)   = H0_DEFAULT
    H0_REF(2)   = H0ERR_DEFAULT
    OMAT_REF(1) = OMAT_DEFAULT
    OMAT_REF(2) = OMATERR_DEFAULT
    OLAM_REF(1) = OLAM_DEFAULT
    OLAM_REF(2) = OLAMERR_DEFAULT
    ORAD_REF(1) = ORAD_DEFAULT
    ORAD_REF(2) = 0.0
    W0_REF(1)   = W0_DEFAULT
    W0_REF(2)   = W0ERR_DEFAULT
    WA_REF(1)   = WA_DEFAULT
    WA_REF(2)   = WAERR_DEFAULT

    zTOL_HELIO2CMB   = 0.02   ! tolerance on dzHELIO/zHELIO
    QUANTILE_ZERRMIN = 0.00   ! 0 -> always use quantiles

! set LTEST logicals to false.

    LTEST_KCOR    = .FALSE.
    LTEST_INTERP  = .FALSE.
    LTEST_U3BAND  = .FALSE.
    LTEST_MAG     = .FALSE.

    ABORT_ON_NOEPOCHS  = .TRUE.
    ABORT_ON_BADAVWARP = .TRUE.
    ABORT_ON_BADZ      = .TRUE.
    ABORT_ON_BADKCOR   = .TRUE.
    ABORT_ON_BADSURVEY = .TRUE.
    ABORT_ON_BADFILTER = .TRUE.
    ABORT_ON_MARGPDF0  = .TRUE.
    ABORT_ON_TRESTCUT  = .TRUE. ! for photo-z fits only
    ABORT_ON_DUPLCID   = .TRUE.
    ABORT_ON_DUPLMJD   = .FALSE.  ! Jun 2017
    ABORT_ON_NOPKMJD   = .TRUE.

    SNTABLE_LIST          = 'DEFAULT'  ! Sep 08 2014
    SNTABLE_FILTER_REMAP  = ''
    WRTABLEFILE_SIMVAR    = .TRUE.     ! Jul 17 2016
    WRTABLEFILE_ZPHOT     = .TRUE.     ! Mar 19 2018
    WRTABLEFILE_HOST_TEXT  = .FALSE.    ! May 23 2019
    WRTABLEFILE_HOST2_TEXT = .FALSE.    ! Sep 2023
    WRTABLEFILE_ERRCALC_TEXT = .FALSE. ! Oct 30 2019

    SNTABLE_APPEND_VARNAME = ''
    SNTABLE_APPEND_VALUE   = -9.0

! ----------------------------------------
! check command line arg for namelist file;

     CALL GETARG(1, nmlfile )
     if ( nmlFile(1:6) .EQ. 'NOFILE' ) GOTO 400  ! Jan 5 2016


! ------------------------
! read the nml

300   CONTINUE

    LL = INDEX(nmlfile,' ' ) - 1
    print*,'   Read namelist file: ', nmlfile(1:LL)

    OPEN (  & 
          UNIT   = LUNNML  & 
         ,file   = nmlfile  & 
         ,status = 'OLD'  & 
         ,ERR    = 900  & 
             )

    READ(LUNNML, NML=SNLCINP, ERR=901, END=902 )

    CLOSE ( UNIT = LUNNML )

400   CONTINUE

! check for namelist over-rides from command line.
    CALL SNLCINP_OVERRIDE()  ! refact with MATCH_NMLKEY (Enable 8.19.2020)

! check for cut on simulated search eff;
! if set, then require SEARCHEFF_MASK = 3 =>
! bits 1 & 2 are set.

    IF ( LSIM_SEARCH_SPEC ) THEN
      cutwin_searcheff_mask(1) = 2.9
      cutwin_searcheff_mask(2) = 3.1
    ENDIF

    IF ( LSIM_SEARCH_ZHOST ) THEN  ! Jun 19 2018
      cutwin_searcheff_mask(1) = 4.9
      cutwin_searcheff_mask(2) = 5.1
    ENDIF

! Mar 31 2021: check PHOTFLAG_BITLIST (LSB=0) and compute PHOTFLAG_MSKREJ
    IF ( PHOTFLAG_BITLIST_REJECT(1) .GE. 0 ) THEN
       PHOTFLAG_MSKREJ(1) = 0
       i = 1
       DO WHILE ( PHOTFLAG_BITLIST_REJECT(i) .GE. 0 )
          bit = PHOTFLAG_BITLIST_REJECT(i)
          PHOTFLAG_MSKREJ(1) = IBSET(PHOTFLAG_MSKREJ(1),BIT)
          i = i + 1
       END DO
    ENDIF


#if defined(SNANA)
    NFIT_ITERATION = 0
#endif

    IF ( .not. REDUCE_STDOUT_BATCH ) THEN
      WRITE ( 6 , NML = SNLCINP )
      CALL FLUSH(6)
    ENDIF

#if defined(SNFIT)
    IF ( NFIT_ITERATION .LE. 0 ) THEN
       write(c1err,821) NFIT_ITERATION
821      format('NFIT_ITERATION = ', I3 )
       C2ERR = 'NFIT_ITERATION must be > 0 (see &SNLCINP namelist).'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF
#endif

! -----------------------------------------
! Compute useful variables from namelist

    NIDFIELD_LIST = PARSE_NML_STRINGLIST(SNFIELD_LIST,60)
    DOALL_SNFIELD = ( SNFIELD_LIST(1) .EQ. 'ALL' )

    i = 1 ; NSNTYPE_LIST=0
    DO WHILE ( SNTYPE_LIST(i) .GE. MNTYPE )
        NSNTYPE_LIST = NSNTYPE_LIST + 1
        i = i + 1
    END DO
    if ( NSNTYPE_LIST .GT. MXLISTNML ) then
        write(c1err,678) 'SNTYPE_LIST', NSNTYPE_LIST, MXLISTNML
678       format(A,2x,I3,' elements exceeds bound.')
        CALL  MADABORT(FNAM, C1ERR, ' ')
    endif

    i = 1 ; NDETNUM_LIST=0
    DO WHILE ( DETNUM_LIST(i) > 0 )
        NDETNUM_LIST = NDETNUM_LIST + 1
        i = i + 1
    END DO

    i = 1 ; NSIM_TEMPLATE_INDEX_LIST=0
    DO WHILE ( SIM_TEMPLATE_INDEX_LIST(i) > 0 )
        NSIM_TEMPLATE_INDEX_LIST = NSIM_TEMPLATE_INDEX_LIST+1
        i = i + 1
    END DO

    i = 1 ; NSNTYPE_IGNORE=0
    DO WHILE ( SNTYPE_IGNORE(i) .GE. 0 )
        NSNTYPE_IGNORE = NSNTYPE_IGNORE + 1
        i = i + 1
    END DO
    if ( NSNTYPE_IGNORE .GT. MXLISTNML ) then
        write(c1err,678) 'SNTYPE_IGNORE',NSNTYPE_IGNORE,MXLISTNML
        CALL  MADABORT(FNAM, C1ERR, ' ')
    endif

! ------------------------
! count number of photometry versions

    CALL PREP_VERSION_PHOTOMETRY()

! Feb 27 2018: check ZP cut
    NZPCUT = 0
    if ( CUTWIN_ZPADU(1) > 0. ) then
       CUTWIN_ZP(1) = CUTWIN_ZPADU(1)
       CUTWIN_ZP(2) = CUTWIN_ZPADU(2)
       NZPCUT = NZPCUT+1
    endif
    if ( CUTWIN_ZPNPE(1) > 0. ) then
       CUTWIN_ZP(1) = CUTWIN_ZPNPE(1)
       CUTWIN_ZP(2) = CUTWIN_ZPNPE(2)
       NZPCUT = NZPCUT+1
    endif
    IF ( NZPCUT > 1 ) THEN
       C1ERR = 'Cannot specify both CUTWIN_ZPADU and CUTWIN_ZPNPE'
       C2ERR = 'Pick one of these CUTWIN_ZP options.'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

! -----------------------------

    DO_FIT = NFIT_ITERATION .GT. 0

! init some INIXXX parameter for fit
    IF ( DO_FIT ) THEN
       do ipar = 1, MXFITPAR
          INIBND(1,ipar) =  -1.0E4  ! avoid crazy MINUIT excursions
          INIBND(2,ipar) =  +1.0E4
          INIVAL(ipar)   =  -9.0
          INISTP(ipar)   =   0.0
       enddo
    ENDIF

    IF ( MWEBV_SCALE .NE. 1.0 .or. MWEBV_SHIFT .NE. 0.0 ) THEN
       GLOBAL_BANNER = ''
       write(GLOBAL_BANNER,541) MWEBV_SCALE, MWEBV_SHIFT
541      format('Galactic MWEBV -> MWEBV *', F6.3,'  +  ', F6.3 )
      CALL PRBANNER ( GLOBAL_BANNER )
    ENDIF

    IF ( MWEBV_FORCE > -0.01 ) THEN
       GLOBAL_BANNER = ''
       write(GLOBAL_BANNER,542) MWEBV_FORCE
542      format('Force Galactic MWEBV -> ', F6.3 )
      CALL PRBANNER ( GLOBAL_BANNER )
    ENDIF

    IF ( CUTWIN_PHOTPROB(1) > -0.001 .and.  & 
           PHOTFLAG_DETECT==0) THEN
       C1ERR = 'Cannot define CUTWIN_PHOTPROB without also'
       C2ERR = 'defining PHOTFLAG_DETECT'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

! ---------------------------
    IF ( VALIDFILENAME(VERSION_REFORMAT_FITS) ) THEN
      OPT_REFORMAT_FITS = IBSET(OPT_REFORMAT_FITS,0) ! set 1st bit
    ENDIF
    IF ( VALIDFILENAME(VERSION_REFORMAT_TEXT) ) THEN    ! Feb 2021
      OPT_REFORMAT_TEXT = IBSET(OPT_REFORMAT_TEXT,0) ! set 1st bit
    ENDIF

! - - - - -
    IF ( SIMLIB_OUT .NE. '' ) SIMLIB_OUTFILE = SIMLIB_OUT

! ---------------------------
! if ENV is in path, make substitution
    CALL ENVreplace(PRIVATE_DATA_PATH)

! check 'pwd' option, which is meant to run split_and_fit on
! different clusters using private data in your `pwd' .
! If interactive, replace invalid pwd/ with a valid './' path.

    IF ( PRIVATE_DATA_PATH .EQ. 'pwd' ) THEN
       PRIVATE_DATA_PATH = '.'
    ENDIF

! if PRIVATE_DATA_PATH points to public $SNDATA_ROOT/lcmerge,
! disalbe PRIVATE_DATA_PATH.
    LL = INDEX(SNDATA_ROOT,' ' ) - 1
    IF ( PRIVATE_DATA_PATH.EQ.SNDATA_ROOT(1:LL)//'/lcmerge')THEN
        PRIVATE_DATA_PATH = ' '
    ENDIF

! Jan 2025; set global flags for reading SNIDs from separate file
    USE_SNCID_FILE        = BTEST(OPT_SNCID_LIST,0) ! use CIDs from file
    USE_INIVAL_SNCID_FILE = BTEST(OPT_SNCID_LIST,1) ! set each INIVAL = FITPAR from file
    USE_PRIOR_SNCID_FILE  = BTEST(OPT_SNCID_LIST,2) ! use each FITPAR & ERR as prior


! Nov 15 2019: abort on obsolete key
    IF ( CUTWIN_NBAND_THRESH(1) > -1.0 .or.  & 
           CUTWIN_NBAND_THRESH(2) > -1.0 ) then

       C1ERR = 'Remove obsolete CUTWIN_NBAND_THRESH from &SNLCINP.'
       C2ERR = 'Sorry for the inconvenience.'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

! #################################################

             RETURN

! #################################################

900   c1err = 'Could not open namelist file: '
    c2err = nmlfile(1:80)
    CALL  MADABORT(FNAM, C1ERR, C2ERR)
    GOTO 990

901   C1ERR = 'could not read &SNLCINP namelist'
    c2err = "Check input-namelist file."
    CALL  MADABORT(FNAM, C1ERR, C2ERR )
    GOTO 990

902   C1ERR = 'could not find &SNLCINP namelist'
    c2err = "Check input-namelist file."
    CALL  MADABORT(FNAM, C1ERR, C2ERR )
    GOTO 990

990   IERR = -1
    RETURN
  END SUBROUTINE RDSNNML

! =======================================
    SUBROUTINE PREP_VERSION_PHOTOMETRY()

! Created Mar 2022 by R.Kessler
! Check wildcard option to get photometry folders using glob.
! Abort if too many versions.
! 
! Aug 2022: set abort trap on VERSION string length

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, index_C, LENV, LENF, JSLASH,  NVERS, langFlag
    CHARACTER  & 
         FNAM*24  & 
        ,cSTRING*(MXCHAR_FILENAME), cFILE*(MXCHAR_FILENAME)  & 
        ,VERSION*(2*MXCHAR_VERSION)
    LOGICAL  IGNOREFILE_fortran, IGNORE
    INTEGER  STORE_GLOB_FILE_LIST
    EXTERNAL STORE_GLOB_FILE_LIST
    EXTERNAL GET_GLOB_FILE, RESET_GLOB_FILE_LIST

! ---------- BEGIN -------------

    FNAM = 'PREP_VERSION_PHOTOMETRY'

! check for wild card; e.g.,

    CALL ENVreplace(VERSION_PHOTOMETRY_WILDCARD)
    LENV    = INDEX(VERSION_PHOTOMETRY_WILDCARD,' ') - 1

    IGNORE  = IGNOREFILE_fortran(VERSION_PHOTOMETRY_WILDCARD)

    IF ( .NOT. IGNORE ) THEN
       print*,' '
       print*,' Use glob to find data folders for'
       print*,'     ', VERSION_PHOTOMETRY_WILDCARD(1:LENV)
       call flush(6)

       cSTRING = VERSION_PHOTOMETRY_WILDCARD(1:LENV) // char(0)
       NVERS = STORE_GLOB_FILE_LIST(cSTRING,LENV)

       if ( NVERS == 0 ) then
	    C1ERR = 'Found zero data versions matching wildcard'
	    C2ERR = 'Check &SNLCINP  VERSION_PHOTOMETRY_WILDCARD'
          CALL MADABORT(FNAM, C1ERR, C2ERR)
 	 endif

       langFlag = 1 ! fortran
       DO i = 1, NVERS
          index_C = i - 1  ! C index
          CALL GET_GLOB_FILE(langFlag, index_C, cFILE, LENF)

!  strip off base name after last slash
          JSLASH = INDEX(cFILE,'/', BACK=.TRUE.)
          LENF   = INDEX(cFILE,' ') - 1
          VERSION_PHOTOMETRY(i) = cFILE(JSLASH+1:LENF)

          LENV = INDEX(VERSION_PHOTOMETRY(i),' ') - 1
          print*,' Found ', VERSION_PHOTOMETRY(i)(1:LENV)
       ENDDO
       call flush(6)
       CALL RESET_GLOB_FILE_LIST()
!           STOP ! xxx
    ENDIF

! - - - - -
! Count number of versions, and check string length

    i = 1
    DO WHILE ( VERSION_PHOTOMETRY(i) .NE. 'NULL' )
      N_VERSION = i
      VERSION = VERSION_PHOTOMETRY(i)
      LENV = INDEX(VERSION,' ') - 1
      IF ( LENV .LE. 0 .or. LENV .GE. MXCHAR_VERSION ) THEN
        CALL PRINT_PREABORT_BANNER(FNAM//char(0),40)
        print*,' VERSION = ', VERSION
        print*,' LEN(VERSION) = ', LENV
        print*,' Beware that MXCHAR_VERSION = ', MXCHAR_VERSION
        C1ERR = 'Invalid string len for VERSION'
        write(c2err,17) MXCHAR_VERSION
17        format('VERSION is either blank, or exceeds ',I3,' chars')
        CALL MADABORT(FNAM, C1err, C2err)
      ENDIF


      if ( N_VERSION .GT. MXVERS ) then
        write(C1ERR,61) N_VERSION, MXVERS
61        format('N_VERSION=',I3,' exceeds MXVERS=',I3)
        C2ERR = 'Reduce number of data versions or increase MXVERS'
        CALL  MADABORT(FNAM, C1ERR, C2ERR )
      endif
      i = i + 1
    END DO

    RETURN
  END SUBROUTINE PREP_VERSION_PHOTOMETRY

! ======================================
    SUBROUTINE CHKARRBOUND(i,MIN,MAX, VARNAME, COMMENT, FUNNAME)

! Nov 2013
! One-line array bound check.
! Abort if index i is outside range MIN to MAX.
! Char input args are used for abort message.


    INTEGER i, MIN, MAX  ! (I) index, MIN,MAX of array
    CHARACTER  & 
         VARNAME*(*)   &  ! (I) name of variable
        ,COMMENT*(*)   &  ! (I) comment string
        ,FUNNAME*(*)  ! (I) name of calling function

! local var

    CHARACTER MSG1*80, MSG2*80

! ------------------- BEGIN ------------------

    IF ( i .GE. MIN .and. i .LE. MAX ) RETURN

    write(MSG1,61) VARNAME, i, MIN, MAX, FUNNAME
 61   format(A,'-index ',I4,' outside valid range ',I2,'-',I5, 3x,  & 
                '(FUN=',A,')'  )

    MSG2 = COMMENT

    CALL  MADABORT("CHKARRBOUND", MSG1, MSG2)

    RETURN
  END SUBROUTINE CHKARRBOUND

! ===================================
    SUBROUTINE NMLCHECK_SNCCID(NMLVARNAME, CCID)

! Apr 26 2013
! Abort if CCID string contains a comma or blank space between
! strings.
! For example, abort on strings such as
!  '04D46,05D446' ,  '04D46 05D446'
! 
! NMLVARNAME is the name of the namelist variable,
! for error message only.
! 
! ---------------- BEGIN --------------


    USE SNPAR

    IMPLICIT NONE

    CHARACTER NMLVARNAME*(*)  ! (I) name of namelist variable
    CHARACTER CCID*(*)        ! (I) SN CID name

! local variables

    INTEGER   JCOMMA, NWD, LEN, MSKOPT
    CHARACTER MSG1*72, MSG2*72, FNAM*30

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------- BEGIN -----------------

    FNAM = 'NMLCHECK_SNCCID'
    JCOMMA   = INDEX(CCID,',')

    LEN      = INDEX(CCID,' ' ) - 1
    MSKOPT   = MSKOPT_PARSE_WORDS_STRING
    NWD      = STORE_PARSE_WORDS(MSKOPT,CCID(1:LEN)//char(0),  & 
                     FNAM//char(0), LEN, 30)

    IF ( JCOMMA > 0 .or. NWD > 1 ) THEN
       MSG1 = 'Invalid nml variable '  & 
                 // NMLVARNAME // ' = ' // CCID
       MSG2 = 'No blanks or commas allowed.'
       CALL  MADABORT("NMLCHECK_SNCCID", MSG1, MSG2)
    ENDIF

    RETURN
  END SUBROUTINE NMLCHECK_SNCCID


! ========================================
    SUBROUTINE  INIT_OUTFILES(NFIT_PER_SN)


! Created Feb 2012
! 
! Call appropriate init function for  root, etc ...
! Note that multiple output files can be initialized
! corresponding to different formats.
! If no files are opened, print warning message but do not abort.
! 
! Note that NFIT_PER_SN is the max number of expected fits per SN.
! 
! Apr 26 2014: replace HFILE_OPEN and ROOTFILE_OPEN with generic
!              wrapper TABLEFILE_OPEN.
! 
! May  3, 2014: call SNLCPAK_INIT only if NFIT_PER_SN > 0
! May 12, 2014: overhaul the ascii format; see TEXTFILE_PREFIX.
!               rename SNFILE_OUT_INIT -> INIT_OUTFILES
! 
! Sep 07, 2014: pass TEXTFORMAT to SNLCPAK_INIT()
! 
! Dec 10 2017: if NFILT_REMAP_TABLE>0, pass REMAP list of filters
!              to SNLCPAK_INIT
! 
! Dec 19 2019: pass SIMFLAG to SNLCPAK_INIT.
! 
! Apr 03 2020: add call to  ENVreplace(TEXTFILE_PREFIX)
! 
! May 02 2020: begin MARZFILE_OUT (to store WFIRST host spectra for MARZ)
! 
! Mar 19 2021: use IGNOREFILE to check usage; e.g., allows
!                ROOTFILE_OUT NONE
! 
! Sep 21 2023: fix mistaken 1:40 truncation of VER_PHOT_forC
! 
! ------------------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER NFIT_PER_SN ! (I) max number of fits to plot per SN

    INTEGER LENFILE,  IFILETYPE, SIMFLAG

    CHARACTER  & 
         FNAM*16  & 
        ,COPT*20  & 
        ,CFILE*(MXCHAR_FILENAME)  & 
        ,TEXTFMT*20

    INTEGER LENV1, LENV2, LENSRVY, LENFILT, LENFMT, LENC
    CHARACTER  & 
         SURVEY_forC*(MXCHAR_SURVEY)  & 
        ,VER_PHOT_forC*(MXCHAR_VERSION)  & 
        ,VER_SNANA_forC*40  & 
        ,FILTERS_forC*(MXFILT_OBS)  & 
        ,TEXTFMT_forC*20  & 
        ,COMMENT_forC*(MXCHAR_VERSION)

! function
    INTEGER  IGNOREFILE, TABLEFILE_OPEN
    EXTERNAL IGNOREFILE, TABLEFILE_OPEN

! --------------- BEGIN -------------

    FNAM = 'INIT_OUTFILES'

    CALL TABLEFILE_INIT();  ! 4/26/2014

! -----------------------------------------
! set strings needed for the SNLCPAK_INIT functio

    LENSRVY      = INDEX(SURVEY_NAME,' ') - 1
    SURVEY_forC  = SURVEY_NAME(1:LENSRVY) // char(0)

    LENV1         = INDEX(VERSION_PHOTOMETRY(1),' ') - 1
    VER_PHOT_forC = VERSION_PHOTOMETRY(1)(1:LENV1) // char(0)

    LENV2         = INDEX(SNANA_VERSION,' ') - 1
    VER_SNANA_forC = SNANA_VERSION(1:LENV2) // char(0)

    LENFILT      = INDEX(SURVEY_FILTERS,' ') - 1
    FILTERS_forC = SURVEY_FILTERS(1:LENFILT) // char(0)
    if ( NFILT_REMAP_TABLE > 0 ) then  ! Dec 2017
       LENFILT      = INDEX(FILTLIST_REMAP_TABLE,' ') - 1
       FILTERS_forC = FILTLIST_REMAP_TABLE(1:LENFILT) // char(0)
    endif

    TEXTFMT      = TEXTFORMAT_TABLE(ITABLE_SNLCPAK)
    LENFMT       = INDEX(TEXTFMT,' ') - 1
    TEXTFMT_forC = TEXTFMT(1:LENFMT) // char(0)

    IF ( OPT_TABLE(ITABLE_SNLCPAK) > 0 ) THEN
       SIMFLAG = 0
       IF ( ISJOB_SIM ) SIMFLAG=1
       CALL SNLCPAK_INIT(SURVEY_forC, VER_PHOT_forC, VER_SNANA_forC,  & 
              FILTERS_forC, SIMFLAG, NFIT_PER_SN, TEXTFMT_forC,  & 
              LENSRVY, LENV1, LENV2, LENFILT, LENFMT )
    ENDIF

    IF ( OPT_TABLE(ITABLE_SPECPAK) > 0 ) THEN
       TEXTFMT      = TEXTFORMAT_TABLE(ITABLE_SPECPAK)
       LENFMT       = INDEX(TEXTFMT,' ') - 1
       TEXTFMT_forC = TEXTFMT(1:LENFMT) // char(0)
       CALL SPECPAK_INIT(SURVEY_forC, VER_PHOT_forC, TEXTFMT_forC,  & 
                LENSRVY, LENV1, LENFMT )
    ENDIF

! -------------------------------------------

    USE_TABLEFILE        = .FALSE.
    USE_TABLEFILE_ROOT   = .FALSE.
    USE_TABLEFILE_TEXT   = .FALSE.     ! May 2 2020 (forgot init earlier)
    USE_TABLEFILE_MARZ   = .FALSE.     ! May 2 2020 (new table)

! MARZ (host spectra only)

    LENFILE   = INDEX(MARZFILE_OUT,' ') - 1
    CFILE     = MARZFILE_OUT(1:LENFILE) // char(0)
    IF ( IGNOREFILE(cFILE,LENFILE) == 0 ) THEN
       COPT          = "new" // char(0)
       CALL ENVreplace(MARZFILE_OUT)
       IFILETYPE = TABLEFILE_OPEN(CFILE, COPT, LENFILE, 20)

       USE_TABLEFILE_MARZ = .TRUE.
       USE_TABLEFILE      = .TRUE.

       IF ( IFILETYPE < 0 ) THEN
          C1ERR = 'Could not open MARZ fits-file'
          C2ERR = MARZFILE_OUT(1:LENFILE)
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       ENDIF
    ENDIF


    LENFILE   = INDEX(ROOTFILE_OUT,' ') - 1
    CFILE     = ROOTFILE_OUT(1:LENFILE) // char(0)
    IF ( IGNOREFILE(cFILE,LENFILE) == 0 ) THEN
       COPT          = "new root" // char(0)
       USE_TABLEFILE_ROOT = .TRUE.
       USE_TABLEFILE      = .TRUE.
       IFILETYPE = TABLEFILE_OPEN(CFILE, COPT, LENFILE, 20)

       IF ( IFILETYPE < 0 ) THEN
          C1ERR = 'Could not open ROOT table-file'
          C2ERR = ROOTFILE_OUT(1:LENFILE)
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       ENDIF
    ENDIF


    LENFILE   = INDEX(TEXTFILE_PREFIX,' ') - 1
    CFILE     = TEXTFILE_PREFIX(1:LENFILE) // char(0)
    IF ( IGNOREFILE(cFILE,LENFILE) == 0 ) THEN
       CALL ENVreplace(TEXTFILE_PREFIX)
       COPT          = "new text" // char(0)
       USE_TABLEFILE_TEXT = .TRUE.
       USE_TABLEFILE      = .TRUE.
       IFILETYPE = TABLEFILE_OPEN(CFILE, COPT, LENFILE, 20)

       IF ( IFILETYPE < 0 ) THEN
          C1ERR = 'Could not open TEXT table-files with PREFIX='
          C2ERR = TEXTFILE_PREFIX(1:LENFILE)
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       ENDIF
    ENDIF

! ----------------------------------------
! Oct 23, 2014: store global comments

    LENC = MXCHAR_VERSION

    COMMENT_forC = 'VERSION_PHOTOMETRY:  '//VER_PHOT_forC(1:LENC)
    CALL STORE_TABLEFILE_COMMENT(COMMENT_forC, LENC)

! ----------------------------------------
!     print warning message if there is no table file
    IF ( .not. USE_TABLEFILE ) THEN
       print*,' '
       write(6,60)
       write(6,61)  & 
           'no table-file specified for tables and light curves'
       write(6,60)
       print*,' '
 60      format(T3, 'WARNING: ', 30('<>') )
 61      format(T3, 'WARNING: ', A)
    ENDIF

    RETURN
  END SUBROUTINE  INIT_OUTFILES


! ================================================
    SUBROUTINE INIT_SNTABLE_OPTIONS()
! 
! Created May 2014 by R.Kessler
! 
! Analyze user input string SNTABLE_LIST (&SNLCINP) and set the
! following arrays vs. sparse ITABLE index:
! 
!   OPT_TABLE         ! 0,1 --> on,off;  >1 --> other option
!   TEXTFORMAT_TABLE  ! key, csv, col, none
! 
! This input string controls which tables are produced
! (SNANA, FITRES, SNLCPAK) and what text formats if
! TEXTFILE_PREFIX is given.
! 
! 
! Apr 26 2017: check pre-scale option: e.g.  'SNANA(ps:100)'
! Mar 19 2018: check nozphot option
! Jan 29 2019: serach for last non-null char in SNTABLE_LIST to get LENTB
! Jan 30 2019: MXTYPE -> 1000 (was 400)
! May 09 2019:
!  + ITABLE_SPECTRA renamed to ITABLE_MODELSPEC (for SALT2 only)
!       (to avoid confusion with ITABLE_SPECPAK
!  + fix bug parsing options with comma; see new MSKOPT argument
!    to store_PARSE_WORDS()
! 
! May 23 2019: add text:host    option ; see WRTABLEFILE_HOST_TEXT
! Oct 30 2019: add text:errcalc option ; see WRTABLEFILE_ERRCALC_TEXT
! Jul 13 2020: no default SNTABLEs if writing rejected epochs
! Mar 17 2021: begin OUTLIER table
! Apr 25 2021: OUTLIER table is automatically disabled for large BIASCOR sim
! Apr 19 2022: disable LCPLOT table for sims in batch mode
!                (to avoid excessive output)
! 
! -------------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER  & 
         itab, iwd, LENS, LENTB, LENOPT, LEN1,LEN2, NOPT, iopt,  & 
         PS, NWD, i, MSKOPT, j_colon

    CHARACTER USERSTRING*60, locase_STRING*60, TBNAME*60
    CHARACTER COPT_LIST(10)*20, COPT*40, FNAM*22, C_ARG*100
    CHARACTER cTMPOPT*60
    LOGICAL   MATCH, USEOPT, DO_SNTABLE, DISABLE

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ---------------- BEGIN ----------------

    FNAM = 'INIT_SNTABLE_OPTIONS'

! if reformating option is set, turn off all table output
    IF ( OPT_REFORMAT_SALT2 > 0 ) SNTABLE_LIST = ' '
    IF ( OPT_REFORMAT_FITS  > 0 ) SNTABLE_LIST = ' '
    IF ( OPT_REFORMAT_TEXT  > 0 ) SNTABLE_LIST = ' '

! init options to nothing created.

    DO itab = 1, MXTABLE
       OPT_TABLE(itab)          = 0
       PRESCALE_TABLE(itab)     = 1.0
    ENDDO

! Set default text format to output text-table for FITRES only:
! with default SNTABLE_LIST = 'SNANA  FITRES  LCPLOT' ,
! setting TEXTFILE_PREFIX  results in an ascii table for the
! FITRES table only ... just as before. To get ascii table for
! SNANA or LC, user must explicitly request it in the &SNLCINP
! namelist, SNTABLE_LIST = 'SNANA(text:[format])' and SNANA->LCPLOT.

! xxx mark delete      TEXTFORMAT_TABLE(ITABLE_SNANA)     = 'none'
    TEXTFORMAT_TABLE(ITABLE_SNANA)     = 'key'  ! Mar 2025
    TEXTFORMAT_TABLE(ITABLE_FITRES)    = 'key'
    TEXTFORMAT_TABLE(ITABLE_OUTLIER)   = 'key'
    TEXTFORMAT_TABLE(ITABLE_SNLCPAK)   = 'key'  ! xxx 'none'
    TEXTFORMAT_TABLE(ITABLE_SPECPAK)   = 'none' ! Apr 2019
    TEXTFORMAT_TABLE(ITABLE_MODELSPEC) = 'none'
    TEXTFORMAT_TABLE(ITABLE_MARZ)      = 'none'


#if defined(SNANA)
!  ISJOB_SNANA isn't set yet, so use preproc flag
    TEXTFORMAT_TABLE(ITABLE_SNANA) = 'key' ! Mar 2021
#endif

    CUTMASK_SNANA_TABLE = 63  ! default is to add everything in snana table

! -------------------------
    CALL PRBANNER(FNAM)
! ----------------------------------------------
! If using default, replace 'DEFAULT' with actual default value.

    IF ( SNTABLE_LIST .EQ. 'DEFAULT' ) THEN
       SNTABLE_LIST = SNTABLE_LIST_DEFAULT

!  No default SNTABLEs if writing rejected epochs
       IF ( OUT_EPOCH_IGNORE_FILE .NE. ''  ) SNTABLE_LIST = ' '

!  No default SNTABLEs if there is no NML file (Feb 2021)
       IF ( nmlFile(1:6) .EQ. 'NOFILE' ) SNTABLE_LIST = ' '
    ENDIF

! by default the SNID is written out as is
    WRTABLEFILE_IAUC = .FALSE.

! by default, SIM_XXX (truth) values are included for simulation
    WRTABLEFILE_SIMVAR = .TRUE.

! by default, ZPHOT[ERR] and ZPHOT-covariances are included
    WRTABLEFILE_ZPHOT  = .TRUE.

! for TEXT format, exclude HOSTGAL info by default
    WRTABLEFILE_HOST_TEXT  = .FALSE.
    WRTABLEFILE_HOST2_TEXT = .FALSE.

! --------------------------------------
! bail of there is nothing to do.
    DO_SNTABLE = .TRUE.
    IF ( SNTABLE_LIST      .EQ. ' '    ) DO_SNTABLE = .FALSE.
    IF ( SNTABLE_LIST(1:4) .EQ. 'NONE' ) DO_SNTABLE = .FALSE.

    IF ( .NOT. DO_SNTABLE    ) RETURN

! ---------------------------------------
! parse SNTABLE_LIST string

    LENTB = 0
    DO i = 1, MXCHAR_FILENAME
      if ( SNTABLE_LIST(i:i) .NE. ' ' ) LENTB=i  ! last non-null char
    ENDDO
    C_ARG = SNTABLE_LIST(1:LENTB) // char(0)

    MSKOPT = MSKOPT_PARSE_WORDS_STRING+MSKOPT_PARSE_WORDS_IGNORECOMMA
    NWD   = STORE_PARSE_WORDS(MSKOPT, C_ARG,  & 
                  FNAM//char(0), MXCHAR_FILENAME, 30)

    DO 100 iwd = 1, NWD

       CALL get_PARSE_WORD_fortran(iwd, USERSTRING, LENS);

! convert everything to lower case for parsing so that
! the input is case-insensitive
       CALL loCase(USERSTRING, locase_STRING)

! split string into table name and optional option inside ().
! STRING is of the form "tbName(copt_list)".

       CALL PARSE_SNTABLE_LIST(locase_STRING, TBNAME, NOPT, COPT_LIST)

! check MARZ outpout option (May 2020)
       if ( MARZFILE_OUT .NE. ' ' ) OPT_TABLE(ITABLE_SPECPAK)=2

! - - - - - - - - - - - - - - - - - -
! check table-names

       if ( TBname(1:5) .EQ. 'snana' ) then
          ITAB = ITABLE_SNANA
          OPT_TABLE(ITAB) = 1
          if ( TBname(1:11) .EQ. 'snana+epoch' ) then
             OPT_TABLE(ITAB) = 2     ! include epoch info (Mar 2015)
          endif
          if ( TBname(1:16) .EQ. 'snana+sim_magobs'  & 
                 .and. LSIM_SNANA) then
             OPT_TABLE(ITAB) = 4     ! include SIM_MAGOBS vs. MJD (1/2019)
          endif

       else if ( TBname(1:6) .EQ. 'fitres' ) then
          ITAB = ITABLE_FITRES
          OPT_TABLE(ITAB) = 1
          if ( TBname(1:15) .EQ. 'fitres+residual' ) then
             OPT_TABLE(ITAB) = 2     ! include epoch fit resids
          endif

       else if ( TBname(1:7) .EQ. 'outlier' ) then
!           OUTLIER table is automatically disabled for large BIASCOR sim
          if ( SIM_BIASCOR_MASK == 0 ) then
            ITAB = ITABLE_OUTLIER
            OPT_TABLE(ITAB) = 1
            NSIGCUT_OUTLIER   = 3.0   ! default outlier def is 3sigma
            FTRUECUT_OUTLIER  = 0.01  ! default requires Ftrue > 0
            SBMAGCUT_OUTLIER  = 999.  ! default is no SBMAG cut (Aug 2021)
          endif
       else if ( TBname(1:6) .EQ. 'lcplot' ) then
          ITAB = ITABLE_SNLCPAK
          DISABLE = ( ISJOB_BATCH .and. ISJOB_SIM )
	    if ( TBname(1:12) .EQ. 'lcplot_batch' ) DISABLE = .FALSE.  ! Oct 27 2025; special key to allow
          IF ( DISABLE ) THEN
             NOPT = 0
             print*,  & 
                 '  WARNING: disable LCPLOT table for sim-data in batch'
             call flush(6)
          ELSE
             OPT_TABLE(ITAB) = 1
          ENDIF
       else if ( TBname(1:8) .EQ. 'specplot' ) then
          ITAB = ITABLE_SPECPAK  ! spectra in the data files
          OPT_TABLE(ITAB) = 1

       else if ( TBname(1:9) .EQ. 'modelspec' ) then
          ITAB = ITABLE_MODELSPEC  ! SALT2 model spec
          OPT_TABLE(ITAB) = 1
       else if ( TBname(1:7) .EQ. 'spectra' ) then
          ITAB = ITABLE_MODELSPEC  ! legacy key for SALT2 model spec
          OPT_TABLE(ITAB) = 1

       else if ( TBname(1:8) .EQ. 'nosimvar' ) then ! option, not a table
          WRTABLEFILE_SIMVAR = .FALSE.  ! July 2016

       else if ( TBname(1:7) .EQ. 'nozphot' ) then ! option, not a table
          WRTABLEFILE_ZPHOT = .FALSE.  ! Mar 2018

       else
          C1ERR = 'Unrecognized table name : ' // TBname(1:20)
          C2ERR = 'Check manual'
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       endif

! check for options inside ()
       DO 200 iopt = 1, NOPT
          COPT     = COPT_LIST(iopt)
          LENOPT   = index(COPT,' ') - 1
          USEOPT   = .FALSE.
          j_colon  = INDEX(copt,':')

          MATCH= ( copt(1:5) .EQ. 'text:' )
          if ( MATCH ) then
             cTMPOPT = copt(6:LENOPT) // ' '
             if( cTMPOPT(1:5) .EQ. 'host2' ) then
                WRTABLEFILE_HOST2_TEXT = .TRUE.  ! for SNANA table
             else if( cTMPOPT(1:4) .EQ. 'host' ) then
                WRTABLEFILE_HOST_TEXT = .TRUE.  ! for SNANA table
             else if ( cTMPOPT(1:7) .EQ. 'errcalc' ) then
                WRTABLEFILE_ERRCALC_TEXT = .TRUE.  ! for LCPLOT table
             else
                TEXTFORMAT_TABLE(itab) = cTMPOPT(1:8) ! format
             endif
             USEOPT = .TRUE.
          endif

          MATCH = (copt(1:5) .EQ. 'mxlc:' )
          if ( MATCH  .and. ITAB .EQ. ITABLE_SNLCPAK) then
             read(copt(6:LENOPT),*) MXLC_PLOT
             USEOPT = .TRUE.
          endif

          MATCH = ( copt(1:8) .EQ. 'cutmask:' )
          if ( MATCH  .and. ITAB .EQ. ITABLE_SNANA ) then
             read(copt(9:LENOPT),*) CUTMASK_SNANA_TABLE
             USEOPT = .TRUE.
          endif

          MATCH = ( copt(1:4) .EQ. 'iauc' )
          if ( MATCH ) then
             WRTABLEFILE_IAUC = .TRUE.
             USEOPT           = .TRUE.
          endif

          MATCH = ( copt(1:3) .EQ. 'ps:' )
          if ( MATCH ) then
             read(copt(4:LENOPT),*) PS
             PRESCALE_TABLE(itab) = PS
             USEOPT           = .TRUE.
          endif

          MATCH = copt(1:5) .EQ. 'nsig:'
          IF ( MATCH ) THEN
             read(copt(6:LENOPT),*) NSIGCUT_OUTLIER
             USEOPT = .TRUE.
          ENDIF
          MATCH = copt(1:6) .EQ. 'ftrue:'
          IF ( MATCH ) THEN
             read(copt(7:LENOPT),*) FTRUECUT_OUTLIER
             USEOPT = .TRUE.
          ENDIF
          MATCH = copt(1:6) .EQ. 'sbmag:'
          IF ( MATCH ) THEN
             read(copt(7:LENOPT),*) SBMAGCUT_OUTLIER
             USEOPT = .TRUE.
          ENDIF

          if ( .NOT. USEOPT ) then
             LENTB = index(TBNAME,' ') - 1
             C1ERR = 'Invalid option: ' // COPT(1:LENOPT)
             C2ERR = 'for Table = ' // TBNAME(1:LENTB)
             CALL  MADABORT(FNAM, C1ERR, C2ERR)
          endif

 200     CONTINUE

 100  CONTINUE  ! end loop over table names

! ------------
! if LCPLOT is not specified, set MXLC_PLOT to zero.
    IF  ( OPT_TABLE(ITABLE_SNLCPAK) == 0 .and.  & 
            OPT_TABLE(ITABLE_SPECPAK) == 0  )    MXLC_PLOT = 0

! -------------
! summarize

    ITAB = ITABLE_SNANA
    write(6,40) 'SNANA',  OPT_TABLE(ITAB), TEXTFORMAT_TABLE(ITAB)

    ITAB = ITABLE_FITRES
    write(6,40) 'FITRES', OPT_TABLE(ITAB), TEXTFORMAT_TABLE(ITAB)

    ITAB = ITABLE_SNLCPAK
    write(6,41) 'SNLCPAK', MXLC_PLOT, TEXTFORMAT_TABLE(ITAB)
    ITAB = ITABLE_SPECPAK
    write(6,41) 'SPECPAK', MXLC_PLOT, TEXTFORMAT_TABLE(ITAB)

    ITAB = ITABLE_MODELSPEC
    write(6,40) 'MODELSPEC', OPT_TABLE(ITAB),TEXTFORMAT_TABLE(ITAB)

    ITAB = ITABLE_MARZ
    write(6,40) 'MARZ', OPT_TABLE(ITAB),TEXTFORMAT_TABLE(ITAB)

 40   format(T5,A9,'-table: USE=',I1,       3x, 'TEXTFORMAT=',A4 )
 41   format(T5,A9,'-table: MXLC_PLOT=',I8, 3x, 'TEXTFORMAT=',A4 )

    write(6,51) 'SNANA', CUTMASK_SNANA_TABLE
51    format(T5,A9,'-table: select CUTFLAG_SNANA & ', I3 )

    if ( WRTABLEFILE_HOST_TEXT ) then
       print*,'   Include best-match HOSTGAL info in TEXT tables.'
    endif
    if ( WRTABLEFILE_HOST2_TEXT ) then
       print*,'   Include  HOSTGAL & HOSTGAL2 info in TEXT tables.'
    endif

    if ( WRTABLEFILE_IAUC ) then
       print*,'    All tables: append column with IAUC name  '
    endif

    call flush(6)

! --------

! check option to remap filters for table output.
! e..g, multiple V-bands that are similar in different
! subsurveys can be combined into a single V band.
! See manual explanation of SNTABLE_FILTER_REMAP

    IF ( SNTABLE_FILTER_REMAP .EQ. '' ) THEN
      LEN1 = 0
    ELSE
      LEN1 = INDEX(SNTABLE_FILTER_REMAP,'->', BACK=.TRUE.) + 3
    ENDIF

    LEN2 = INDEX(SURVEY_FILTERS,' ') - 1

    CALL  FILTER_REMAP_INIT(  & 
           SNTABLE_FILTER_REMAP(1:LEN1)//char(0),    &  ! (I) user-input string
           SURVEY_FILTERS(1:LEN2)//char(0),          &  ! (I) valid filters
           NFILT_REMAP_TABLE,            &  ! (O) number of remapped bands
           IFILTLIST_REMAP_TABLE,        &  ! (O) index list of remapped bands
           FILTLIST_REMAP_TABLE,         &  ! (O) string list of ...
           LEN1, LEN2, MXFILT_ALL )

    FILTLIST_REMAP_TABLE =  & 
        FILTLIST_REMAP_TABLE(1:NFILT_REMAP_TABLE) // ' '

! option to compute width variables
    CALL GET_SIM_LCWIDTH(0)

    RETURN
  END SUBROUTINE INIT_SNTABLE_OPTIONS

! =================================
    SUBROUTINE GET_SIM_LCWIDTH(MODE)

! Feb 2019
! For each filter, compute and store SIM_LCWIDTH(ifilt)
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER MODE  ! (I) 0=init, 1=get width

    INTEGER, PARAMETER :: MXEP_LCWIDTH = 200

    REAL*8  SIM_MAGOBS(MXEP_LCWIDTH), SIM_MJD(MXEP_LCWIDTH)
    REAL*8  WIDTH
    INTEGER  & 
        IFILT, IFILT_OBS, EP, EPMIN, EPMAX, NEWMJD, NOBS, ERRFLAG

    CHARACTER FNAM*16
! functions
    REAL*8   get_lightCurveWidth
    EXTERNAL INIT_LIGHTCURVEWIDTH, get_lightCurveWidth

! -------- BEGIN ----------

    FNAM = 'GET_SIM_LCWIDTH'
    IF ( OPTSIM_LCWIDTH == 0 ) RETURN
    IF ( .NOT. LSIM_SNANA    ) RETURN

    IF (MODE == 0 ) THEN
!  init only
       CALL INIT_LIGHTCURVEWIDTH()
       RETURN
    ENDIF

! loop over each filter and strip off MJD and true MAG into
! separate array to pass to get_lightCurveWidth function.

    DO 100 ifilt = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       NOBS = 0
       DO 150  NEWMJD = 1, ISNLC_NEWMJD_STORE
          EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
          EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
          DO 151 EP = EPMIN, EPMAX
             if ( ISNLC_IFILT_OBS(ep) .NE. ifilt_obs ) goto 151
             NOBS = NOBS + 1
             SIM_MJD(NOBS)    = SNLC8_MJD(ep)
             SIM_MAGOBS(NOBS) = SIM_EPMAGOBS(ep)
151        CONTINUE
150      CONTINUE

       WIDTH = get_lightCurveWidth(OPTSIM_LCWIDTH,  & 
                     NOBS,SIM_MJD, SIM_MAGOBS, ERRFLAG, FNAM,16)
       SIM_LCWIDTH(ifilt) = SNGL(WIDTH)

!         print*,' xxx WIDTH: ', ifilt, ifilt_obs, WIDTH
 100  CONTINUE ! end IFILT loop

    RETURN
  END SUBROUTINE GET_SIM_LCWIDTH

! =============================================
    SUBROUTINE PROB_TRUEFLUX_CALC(ep)

! Created Apr 6 2021
! For fakes or SIM, compute chi2 = \Sum[F-Ftrue]/sigma_F^2 per band.
! Compute PROB(chi2,NDOF) for each band.
! 
!  ep = 0  -> init sums to zero
!  ep > 0  -> increment chi2
!  ep = -1 -> compute prob(chi2,Ndof)
! --------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE TRUECHI2COM

    IMPLICIT NONE

    INTEGER EP  ! (I)

    INTEGER IFILT_OBS, IFILT, NDOF
    REAL CHI2
    REAL*8 CHI8, PROB8

    REAL*8 PROB_CHI2NDOF
    EXTERNAL PROB_CHI2NDOF
! ------------- BEGIN -------------

    IF ( .NOT. LPROB_TRUEFLUX ) RETURN
    IF ( .NOT. ISJOB_SIM      ) RETURN ! fakes or SNANA sim

    IF ( ep == 0 ) THEN
       DO IFILT = 0, NFILTDEF_SURVEY
          NDOF_TRUEFLUX(ifilt) =   -1
          CHI2_TRUEFLUX(ifilt) =  0.0
          PROB_TRUEFLUX(ifilt) = -9.0
       ENDDO
       RETURN
    ENDIF

! - - - - -
    IF ( EP > 0 ) THEN
       IF ( SIM_EPFLUXCAL(ep) > 0.1 ) THEN
       ! increment chi2
          IFILT_OBS = ISNLC_IFILT_OBS(ep)
          IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
          CHI2      = SIM_EPCHI2FLUX(ep)
          CHI2_TRUEFLUX(ifilt) = CHI2_TRUEFLUX(ifilt) + CHI2
          NDOF_TRUEFLUX(ifilt) = NDOF_TRUEFLUX(ifilt) + 1

          CHI2_TRUEFLUX(0) = CHI2_TRUEFLUX(0) + CHI2
          NDOF_TRUEFLUX(0) = NDOF_TRUEFLUX(0) + 1
       ENDIF
    ELSE
       ! compute FITPROB per band
       DO ifilt    = 1, NFILTDEF_SURVEY
          CHI8     = DBLE( CHI2_TRUEFLUX(ifilt) )
          NDOF     = NDOF_TRUEFLUX(ifilt)
          if ( NDOF > 2 ) then
             PROB8    = PROB_CHI2NDOF(CHI8,NDOF)
             PROB_TRUEFLUX(ifilt) = SNGL(PROB8)
          endif

! global with all filters
          CHI8     = DBLE( CHI2_TRUEFLUX(0) )
          NDOF     = NDOF_TRUEFLUX(0)
          if ( NDOF > 2 ) then
             PROB8    = PROB_CHI2NDOF(CHI8,NDOF)
             PROB_TRUEFLUX(0) = SNGL(PROB8)
          endif

      ENDDO
    ENDIF

    RETURN
  END SUBROUTINE PROB_TRUEFLUX_CALC


! =============================================
    SUBROUTINE PARSE_SNTABLE_LIST(USERSTRING, TABLENAME, NOPT, COPT)

! parse USERSTRING from SNTABLE_LIST;
! rdturn name of table and string of options.


    CHARACTER  & 
         USERSTRING*(*)   &  ! (I) user-input string from SNTABLE_LIST
        ,TABLENAME*20     &  ! (O) name of table
        ,COPT(10)*20     ! (O) list of options

    INTEGER NOPT       ! (O) numbef of COPT options

    LOGICAL DONEXT

! local args

    INTEGER jp1, jp2, j, i, LENOPT
    CHARACTER COPT_LIST*60, CTMP*60 ! comma-separated list

! ----------------- BEGIN ----------------

! init output args
    TABLENAME = 'UNKNOWN'
    NOPT = 0
    do i = 1, 10
      COPT(i) = ''
    enddo

! ---------------------

    jp1  = index(USERSTRING,'(')
    jp2  = index(USERSTRING,')')
    if ( jp1 > 0 ) then
       TABLENAME  = USERSTRING(1:jp1-1)
       COPT_LIST  = USERSTRING(jp1+1:jp2-1)
    else
       TABLENAME  = USERSTRING
       RETURN
    endif

! if we get here then parse comma-separated list of options
! in COPT_LIST

    LENOPT = index(COPT_LIST,' ' ) - 1
    DONEXT = .TRUE.
    CTMP   = COPT_LIST

    DO 100 WHILE ( DONEXT )
       NOPT  = NOPT + 1
       j     = index(CTMP,',')
       if ( j .LE. 0  ) j = LENOPT + 1
       COPT(NOPT) = CTMP(1:j-1)
       DONEXT     = j < LENOPT
       CTMP       = CTMP(j+1:LENOPT)

       if ( NOPT > 10 ) then
          CALL  MADABORT("PARSE_SNTABLE_LIST",  & 
                     'NOPT > 10 ???          ',  & 
                     'Something is wrong     ' )
       endif
 100  CONTINUE

    RETURN
  END SUBROUTINE PARSE_SNTABLE_LIST


! ================================================
    SUBROUTINE SNLCINP_OVERRIDE()
! 
! Created July 2020 by R.Kessler [refactored]
! Check command-line args for user-option to override
! namelist values. Check for two styles of override:
! 
! Legacy style is space-sparated such as    CUTWIN_SNRMAX 5 999
! New style has equal sign and no spaces,   CUTWIN_SNRMAX=5,999
! 
! Beware that some arguments are comma-sep lists, so need
! to be careful about when to split by comma, and when
! not so split by comma.
! 
! Oct 22 2020:
!   + double string length for LINE_ARGS and ARGLIST to allow for
!     long comma-sep strings.
! 
! Nov 23 2020: read OPT_YAML
! 
! Dec 08 2020:
!   + fix bug setting USE_LINE_ARGS(iuse) to trap bad args.
!     last few months we were vulnerable to bad args that were
!     ignored instead of causing abort.
! 
! Aug 14 2024
!   call new function FOUND_MATCH_NMLKEY that has logic to figure out
!   valid key with space sep or '=' delimeter. The old refactor code
!   had never been properly integrated until now. Advantage is that
!   FITNML_OVERRIDE can use the same function.
! 
! Aug 25 2025: STRING_LIST*60 -> STRING_LIST*(MXCHAR_PATH) where MXCHAR_PATH=160;
!              Initial motivation is to pass bigger list for SNCID_LIST.
! 
! ----------

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! local var

    INTEGER iArg, LL, ilast, iuse, NVERLOC
    CHARACTER STRING_LIST*(MXCHAR_PATH)
    CHARACTER ARG*(MXCHAR_ARG), ARGLIST(MXKEY_ARGS)*(MXCHAR_ARG)
    LOGICAL LSNIGNORE, DO_STOP, FOUND_MATCH_LEGACY, FOUND_MATCH_REFAC
    LOGICAL FOUND_MATCH, LDMP

! functions
    LOGICAL MATCH_NMLKEY, FOUND_MATCH_NMLKEY
    INTEGER PARSE_INTLIST

! ------------- BEGIN -------------

    LDMP = .FALSE.
    iArg = 2 ;   ilast = 2 ;      NVERLOC = 0

    DO 4444 WHILE ( iArg .LE. NLINE_ARGS )

      CALL GETARG(iArg, ARG)
      LL = INDEX (ARG, ' ' ) - 1

      print*,'     PROCESS COMMAND LINE ARG: ', ARG(1:LL)

      if ( LDMP ) then
        print*,' xxx ===================================== '
        print*,' xxx iArg,ARG = ', iArg, ARG(1:40)
      endif

      CALL FLUSH(6)
      ARGLIST(1) = ''

       if ( MATCH_NMLKEY('VERSION_PHOTOMETRY',1,iArg,ARGLIST) ) then
          CALL PARSE_COMMASEP_LIST('VERSION_PHOTOMETRY',  & 
                                       ARGLIST(1))

       else if(MATCH_NMLKEY('VERSION_PHOTOMETRY_WILDCARD',  & 
                      1, iArg, ARGLIST) ) then
          VERSION_PHOTOMETRY_WILDCARD = ARGLIST(1)(1:MXCHAR_VERSION)
       else if(MATCH_NMLKEY('VERSION_REFORMAT_FITS',  & 
                      1, iArg, ARGLIST) ) then
          VERSION_REFORMAT_FITS = ARGLIST(1)(1:MXCHAR_VERSION)
       else if(MATCH_NMLKEY('VERSION_REFORMAT_TEXT',  & 
                      1, iArg, ARGLIST) ) then
          VERSION_REFORMAT_TEXT = ARGLIST(1)(1:MXCHAR_VERSION)

       else if ( MATCH_NMLKEY('OPT_YAML',  & 
                   1, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) OPT_YAML

       else if ( MATCH_NMLKEY('OPT_REFORMAT_FITS',  & 
                   1, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) OPT_REFORMAT_FITS

       else if ( MATCH_NMLKEY('OPT_REFORMAT_TEXT',  & 
                    1, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) OPT_REFORMAT_TEXT

       else if ( MATCH_NMLKEY('OPT_REFORMAT_SPECTRA',  & 
                   1, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) OPT_REFORMAT_SPECTRA

       else if ( MATCH_NMLKEY('OPT_REFORMAT_SALT2',  & 
                    1, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) OPT_REFORMAT_SALT2

       else if ( MATCH_NMLKEY('REFORMAT_KEYS',  & 
                    1, iArg, ARGLIST) ) then
           REFORMAT_KEYS = ARGLIST(1)(1:100)

       else if ( MATCH_NMLKEY('PRIVATE_DATA_PATH',  & 
                     1, iArg, ARGLIST) ) then
          PRIVATE_DATA_PATH = ARGLIST(1)(1:MXCHAR_PATH)

       else if ( MATCH_NMLKEY('PRIVATE_CUTWIN_STRING',  & 
                  1, iArg, ARGLIST) ) then
          PRIVATE_CUTWIN_STRING = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('PRIVATE_VARNAME_READLIST',  & 
                  1, iArg, ARGLIST) ) then
          PRIVATE_VARNAME_READLIST = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('EARLYLC_STRING',  & 
                   1, iArg, ARGLIST) ) then
         EARLYLC_STRING = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('REQUIRE_EPOCHS_STRING',  & 
                   1, iArg, ARGLIST) ) then
         REQUIRE_EPOCHS_STRING = ARGLIST(1)(1:100)

       else if ( MATCH_NMLKEY('DUMP_STRING',  & 
                   1, iArg, ARGLIST) ) then
         DUMP_STRING = ARGLIST(1)(1:100)

       else if ( MATCH_NMLKEY('FILTER_UPDATE_PATH',  & 
                   1,iArg, ARGLIST) ) then
         FILTER_UPDATE_PATH = ARGLIST(1)(1:MXCHAR_PATH)

       else if ( MATCH_NMLKEY('NONSURVEY_FILTERS',  & 
                   1, iArg, ARGLIST) ) then
         NONSURVEY_FILTERS = ARGLIST(1)(1:MXFILT_ALL)

       else if ( MATCH_NMLKEY('FILTER_REPLACE',  & 
                   1, iArg, ARGLIST) ) then
         FILTER_REPLACE = ARGLIST(1)(1:MXFILT_ALL)

       else if ( MATCH_NMLKEY('FILTLIST_LAMSHIFT',  & 
                   1, iArg, ARGLIST) ) then
         FILTLIST_LAMSHIFT = ARGLIST(1)(1:MXFILT_ALL)

       else if ( MATCH_NMLKEY('USERTAGS_FILE',  & 
                   1, iArg, ARGLIST) ) then
         USERTAGS_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('VPEC_FILE',  & 
                   1, iArg, ARGLIST) ) then
         VPEC_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('CALIB_FILE',  & 
                   1, iArg, ARGLIST) ) then
         CALIB_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)
         KCOR_FILE = ''

       else if ( MATCH_NMLKEY('KCOR_FILE',   &  ! legacy key for CALIB_FILE
                   1, iArg, ARGLIST) ) then
         CALIB_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)
         KCOR_FILE = ''

       else if ( MATCH_NMLKEY('OVERRIDE_RESTLAM_BOUNDARY',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) OVERRIDE_RESTLAM_BOUNDARY

       else if ( MATCH_NMLKEY('HEADER_OVERRIDE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         HEADER_OVERRIDE_FILE = ARGLIST(1)
       else if ( MATCH_NMLKEY('SIM_HEADER_OVERRIDE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SIM_HEADER_OVERRIDE_FILE = ARGLIST(1)

! LSNIGNORE

       else if ( MATCH_NMLKEY('INTERP_OPT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) INTERP_OPT

       else if ( MATCH_NMLKEY('ROOTFILE_OUT',  & 
                   1, iArg, ARGLIST) ) then
         ROOTFILE_OUT = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('TEXTFILE_PREFIX',  & 
                   1, iArg, ARGLIST) ) then
         TEXTFILE_PREFIX = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('MARZFILE_OUT',  & 
                   1, iArg, ARGLIST) ) then
         MARZFILE_OUT = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SNTABLE_LIST',  & 
                   1, iArg, ARGLIST) ) then
         SNTABLE_LIST = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SNTABLE_FILTER_REMAP',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SNTABLE_FILTER_REMAP

       else if ( MATCH_NMLKEY('SNTABLE_APPEND',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SNTABLE_APPEND_VARNAME
         READ(ARGLIST(2),*) SNTABLE_APPEND_VALUE

       else if ( MATCH_NMLKEY('OPTSIM_LCWIDTH',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) OPTSIM_LCWIDTH

       else if ( MATCH_NMLKEY('SNCID_IGNORE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SNCID_IGNORE_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SNCID_IGNORE',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SNCID_IGNORE(1)

       else if ( MATCH_NMLKEY('SNMJD_OUT_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SNMJD_LIST_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SNMJD_OUT_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SNMJD_LIST_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('EPOCH_IGNORE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         EPOCH_IGNORE_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('OUT_EPOCH_IGNORE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         OUT_EPOCH_IGNORE_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('NFIT_ITERATION',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) NFIT_ITERATION

       else if ( MATCH_NMLKEY('USE_MINOS',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USE_MINOS

       else if ( MATCH_NMLKEY('MINUIT_PRINT_LEVEL',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MINUIT_PRINT_LEVEL

       else if ( MATCH_NMLKEY('MXLC_PLOT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MXLC_PLOT

       else if ( MATCH_NMLKEY('MJDPERIOD_PLOT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MJDPERIOD_PLOT

       else if ( MATCH_NMLKEY('MJDSHIFT_PLOT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MJDSHIFT_PLOT

       else if ( MATCH_NMLKEY('DTOBS_MODEL_PLOT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) DTOBS_MODEL_PLOT

       else if ( MATCH_NMLKEY('PHOTFLAG_MSKREJ',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) PHOTFLAG_MSKREJ(1)
         READ(ARGLIST(2),*) PHOTFLAG_MSKREJ(2)

       else if ( MATCH_NMLKEY('PHOTFLAG_BITLIST_REJECT',  & 
                   1, iArg, ARGLIST) ) then
         if ( ARGLIST(1) .NE. 'NONE' ) then
            CALL PARSE_COMMASEP_LIST('PHOTFLAG_BITLIST',ARGLIST(1))
         endif

       else if ( MATCH_NMLKEY('PHOTFLAG_DETECT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) PHOTFLAG_DETECT

       else if ( MATCH_NMLKEY('PHOTFLAG_TRIGGER',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) PHOTFLAG_TRIGGER
! - - - - - -
       else if ( MATCH_NMLKEY('MULTISEASON_OPTMASK',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MULTISEASON_OPTMASK

       else if ( MATCH_NMLKEY('MULTISEASON_TGAP',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MULTISEASON_TGAP

       else if ( MATCH_NMLKEY('MULTISEASON_NREJECT_OUTLIER',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MULTISEASON_NREJECT_OUTLIER

       else if ( MATCH_NMLKEY('MULTISEASON_CHI2RED_ACTIVE',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MULTISEASON_CHI2RED_ACTIVE

       else if ( MATCH_NMLKEY('CUTWIN_NSEASON_ACTIVE',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) CUTWIN_NSEASON_ACTIVE(1)
         READ(ARGLIST(2),*) CUTWIN_NSEASON_ACTIVE(2)

! - - - - -
       else if ( MATCH_NMLKEY('LSIM_SEARCH_SPEC',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LSIM_SEARCH_SPEC

       else if ( MATCH_NMLKEY('LSIM_SEARCH_ZHOST',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LSIM_SEARCH_ZHOST
       else if ( MATCH_NMLKEY('LSIM_SEARCH_zHOST',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LSIM_SEARCH_ZHOST

       else if ( MATCH_NMLKEY('DEBUG_FLAG',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) DEBUG_FLAG

       else if ( MATCH_NMLKEY('LPROB_TRUEFLUX',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LPROB_TRUEFLUX

       else if ( MATCH_NMLKEY('LDMP_SNFAIL',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LDMP_SNFAIL

       else if ( MATCH_NMLKEY('LDMP_SATURATE',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LDMP_SATURATE

       else if ( MATCH_NMLKEY('USESIM_SNIA',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USESIM_SNIA
       else if ( MATCH_NMLKEY('USESIM_NONIA',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USESIM_NONIA

       else if ( MATCH_NMLKEY('USE_SNHOST_ZPHOT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USE_SNHOST_ZPHOT

       else if ( MATCH_NMLKEY('USE_HOSTGAL_PHOTOZ',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USE_HOSTGAL_PHOTOZ

       else if ( MATCH_NMLKEY('USESIM_TRUEFLUX',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USESIM_TRUEFLUX

       else if ( MATCH_NMLKEY('USESIM_REDSHIFT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USESIM_REDSHIFT

! - - - - -
       else if ( MATCH_NMLKEY('RESTORE_MWEBV_ERR_BUG',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) RESTORE_MWEBV_ERR_BUG

       else if ( MATCH_NMLKEY('RESTORE_WRONG_VPEC',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) RESTORE_WRONG_VPEC

       else if ( MATCH_NMLKEY('RESTORE_OVERRIDE_ZBUG',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) RESTORE_OVERRIDE_ZBUG

       else if ( MATCH_NMLKEY('RESTORE_DES5YR',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) RESTORE_DES5YR

! - - - - -
       else if ( MATCH_NMLKEY('REQUIRE_DOCANA',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) REQUIRE_DOCANA

       else if ( MATCH_NMLKEY('USE_MWCOR',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) USE_MWCOR

! - - - -
       else if ( MATCH_NMLKEY('LTEST_U3BAND',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LTEST_U3BAND
       else if ( MATCH_NMLKEY('LTEST_MAG',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LTEST_MAG
       else if ( MATCH_NMLKEY('LTEST_INTERP',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LTEST_INTERP
       else if ( MATCH_NMLKEY('LTEST_KCOR',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) LTEST_KCOR

! - - - -
       else if ( MATCH_NMLKEY('OPT_SETPKMJD',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) OPT_SETPKMJD

       else if ( MATCH_NMLKEY('MJDWIN_SETPKMJD',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MJDWIN_SETPKMJD

       else if ( MATCH_NMLKEY('SNRCUT_SETPKMJD',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SNRCUT_SETPKMJD

       else if ( MATCH_NMLKEY('SHIFT_SETPKMJD',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SHIFT_SETPKMJD

       else if ( MATCH_NMLKEY('SNCCID_LIST',  & 
                   1, iArg, ARGLIST) ) then
         if ( ARGLIST(1) .NE. 'NONE' ) then
            CALL PARSE_COMMASEP_LIST('SNCCID_LIST', ARGLIST(1))
! xxx mark              CUTWIN_CID(1) = 0;  CUTWIN_CID(2)=0
         endif

       else if ( MATCH_NMLKEY('SNCCID_IGNORE',   &  ! July 2023
                   1, iArg, ARGLIST) ) then
         if ( ARGLIST(1) .NE. 'NONE' ) then
            CALL PARSE_COMMASEP_LIST('SNCCID_IGNORE', ARGLIST(1))
         endif

       else if ( MATCH_NMLKEY('SNCID_LIST_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SNCID_LIST_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)
! xxx mark           CUTWIN_CID(1) = 0;  CUTWIN_CID(2)=0

       else if ( MATCH_NMLKEY('OPT_SNCID_LIST',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) OPT_SNCID_LIST


       else if ( MATCH_NMLKEY('MXEVT_PROCESS',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MXEVT_PROCESS

       else if ( MATCH_NMLKEY('MXEVT_CUTS',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MXEVT_CUTS

       else if ( MATCH_NMLKEY('MXLC_FIT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) MXLC_FIT

       else if ( MATCH_NMLKEY('SIM_TEMPLATE_INDEX_LIST',  & 
                   1, iArg, ARGLIST) ) then
         STRING_LIST = ARGLIST(1)(1:MXCHAR_PATH)
         NSIM_TEMPLATE_INDEX_LIST =   &  ! select SIM_TEMPLATE INDICE
              PARSE_INTLIST(STRING_LIST,SIM_TEMPLATE_INDEX_LIST)

       else if ( MATCH_NMLKEY('SIMLIB_OUT',  & 
                   1, iArg, ARGLIST) ) then
         SIMLIB_OUT = ARGLIST(1)(1:MXCHAR_FILENAME)
       else if ( MATCH_NMLKEY('SIMLIB_OUTFILE',  &  ! same as SIMLIB_OUT
                   1, iArg, ARGLIST) ) then
         SIMLIB_OUTFILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SIMLIB_ZPERR_LIST',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SIMLIB_ZPERR_LIST

       else if ( MATCH_NMLKEY('OPT_SIMLIB_OUT',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) OPT_SIMLIB_OUT

       else if ( MATCH_NMLKEY('SIMLIB_OUT_TMINFIX',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SIMLIB_OUT_TMINFIX

       else if ( MATCH_NMLKEY('NONLINEARITY_FILE',  & 
                   1, iArg, ARGLIST) ) then
         NONLINEARITY_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SNFIELD_LIST',  & 
                   1, iArg, ARGLIST) ) then
         SNFIELD_LIST(1) = ARGLIST(1)(1:MXCHAR_PATH)
         SNFIELD_LIST(2) = ''

       else if ( MATCH_NMLKEY('FORCE_STDOUT_BATCH',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) FORCE_STDOUT_BATCH

       else if ( MATCH_NMLKEY('JOBSPLIT',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) JOBSPLIT(1)
         READ(ARGLIST(2),*) JOBSPLIT(2)
         ISJOB_BATCH = .TRUE.
         REDUCE_STDOUT_BATCH = ( JOBSPLIT(1) .GT. 1 )

       else if ( MATCH_NMLKEY('JOBSPLIT_EXTERNAL',  & 
                   2, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) JOBSPLIT_EXTERNAL(1)
         READ(ARGLIST(2),*) JOBSPLIT_EXTERNAL(2)

       else if ( MATCH_NMLKEY('SIM_PRESCALE',  & 
                   1, iArg, ARGLIST) ) then
         READ(ARGLIST(1),*) SIM_PRESCALE

      else if ( MATCH_NMLKEY('CUTWIN_CID',  & 
                  2, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) CUTWIN_CID(1)
           READ(ARGLIST(2),*) CUTWIN_CID(2)

      else if ( MATCH_NMLKEY('CUTWIN_SNTYPE',  & 
                  2, iArg, ARGLIST) ) then
           READ(ARGLIST(1),*) CUTWIN_SNTYPE(1)
           READ(ARGLIST(2),*) CUTWIN_SNTYPE(2)

! - - - -
       else if ( MATCH_NMLKEY('SNTYPE_LIST',  & 
                   1, iArg, ARGLIST) ) then
         STRING_LIST = ARGLIST(1)(1:MXCHAR_PATH)
         NSNTYPE_LIST = PARSE_INTLIST(STRING_LIST,SNTYPE_LIST)

       else if ( MATCH_NMLKEY('DETNUM_LIST',  & 
                   1, iArg, ARGLIST) ) then
         STRING_LIST = ARGLIST(1)(1:MXCHAR_PATH)
         NDETNUM_LIST = PARSE_INTLIST(STRING_LIST,DETNUM_LIST)

       else if ( MATCH_NMLKEY('SNTYPE_IGNORE',  & 
                   1, iArg, ARGLIST) ) then
         STRING_LIST    = ARGLIST(1)(1:MXCHAR_PATH)
         NSNTYPE_IGNORE = PARSE_INTLIST(STRING_LIST,SNTYPE_IGNORE)

       else if ( MATCH_NMLKEY('SNCID_LIST',   &  ! Aug 18 2020
                   1, iArg, ARGLIST) ) then
         STRING_LIST = ARGLIST(1)(1:MXCHAR_PATH)
         NCID_LIST   = PARSE_INTLIST(STRING_LIST,SNCID_LIST)

! xxxx mark delete Jun 26 2025 xxx
!           if ( SNCID_LIST(1) > 0 ) then
!               CUTWIN_CID(1) = 0;  CUTWIN_CID(2)=0
!           endif
! xxx end mark xxx
! - - - -
      else if ( MATCH_NMLKEY('CUTWIN_NEPOCH',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NEPOCH(1)
        READ(ARGLIST(2),*) CUTWIN_NEPOCH(2)

      else if ( MATCH_NMLKEY('CUTWIN_TRESTMIN',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TRESTMIN(1)
        READ(ARGLIST(2),*) CUTWIN_TRESTMIN(2)

      else if ( MATCH_NMLKEY('CUTWIN_TRESTMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TRESTMAX(1)
        READ(ARGLIST(2),*) CUTWIN_TRESTMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_TRESTRANGE',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TRESTRANGE(1)
        READ(ARGLIST(2),*) CUTWIN_TRESTRANGE(2)

      else if ( MATCH_NMLKEY('CUTWIN_TREST2',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TREST2(1)
        READ(ARGLIST(2),*) CUTWIN_TREST2(2)

      else if ( MATCH_NMLKEY('CUTWIN_TGAPMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TGAPMAX(1)
        READ(ARGLIST(2),*) CUTWIN_TGAPMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_T0GAPMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_T0GAPMAX(1)
        READ(ARGLIST(2),*) CUTWIN_T0GAPMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_TREST',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TREST(1)
        READ(ARGLIST(2),*) CUTWIN_TREST(2)

      else if ( MATCH_NMLKEY('CUTWIN_TOBS',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TOBS(1)
        READ(ARGLIST(2),*) CUTWIN_TOBS(2)

      else if ( MATCH_NMLKEY('CUTWIN_TOBSMIN',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TOBSMIN(1)
        READ(ARGLIST(2),*) CUTWIN_TOBSMIN(2)

      else if ( MATCH_NMLKEY('CUTWIN_TOBSMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TOBSMAX(1)
        READ(ARGLIST(2),*) CUTWIN_TOBSMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_RA',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_RA(1)
        READ(ARGLIST(2),*) CUTWIN_RA(2)

      else if ( MATCH_NMLKEY('CUTWIN_DEC',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_DEC(1)
        READ(ARGLIST(2),*) CUTWIN_DEC(2)

      else if ( MATCH_NMLKEY('CUTWIN_PEAKMJD',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_PEAKMJD(1)
        READ(ARGLIST(2),*) CUTWIN_PEAKMJD(2)

      else if ( MATCH_NMLKEY('CUTWIN_MJD',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_MJD(1)
        READ(ARGLIST(2),*) CUTWIN_MJD(2)

      else if ( MATCH_NMLKEY('CUTWIN_MJD_EXCLUDE',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_MJD_EXCLUDE(1)
        READ(ARGLIST(2),*) CUTWIN_MJD_EXCLUDE(2)

! - - -
      else if ( MATCH_NMLKEY('CUTWIN_NOBS_PREDETECT',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NOBS_PREDETECT(1)
        READ(ARGLIST(2),*) CUTWIN_NOBS_PREDETECT(2)

      else if ( MATCH_NMLKEY('CUTWIN_TOBS_PREDETECT',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TOBS_PREDETECT(1)
        READ(ARGLIST(2),*) CUTWIN_TOBS_PREDETECT(2)

! - - -
      else if ( MATCH_NMLKEY('CUTWIN_REDSHIFT',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_REDSHIFT(1)
        READ(ARGLIST(2),*) CUTWIN_REDSHIFT(2)

      else if ( MATCH_NMLKEY('CUTWIN_REDSHIFT_ERR',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_REDSHIFT_ERR(1)
        READ(ARGLIST(2),*) CUTWIN_REDSHIFT_ERR(2)

      else if ( MATCH_NMLKEY('CUTWIN_HOSTSEP',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_HOSTSEP(1)
        READ(ARGLIST(2),*) CUTWIN_HOSTSEP(2)

      else if ( MATCH_NMLKEY('CUTWIN_NFILT_TRESTMIN',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFILT_TRESTMIN(1)
        READ(ARGLIST(2),*) CUTWIN_NFILT_TRESTMIN(2)

      else if ( MATCH_NMLKEY('CUTWIN_NFILT_TRESTMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFILT_TRESTMAX(1)
        READ(ARGLIST(2),*) CUTWIN_NFILT_TRESTMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_NFILT_TREST2',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFILT_TREST2(1)
        READ(ARGLIST(2),*) CUTWIN_NFILT_TREST2(2)

      else if ( MATCH_NMLKEY('SNCUT_SNRMAX',  & 
                  1, iArg, ARGLIST) ) then
        SNCUT_SNRMAX = ARGLIST(1)(1:MXCHAR_CUTNAME)

      else if ( MATCH_NMLKEY('SNCUT_HOST_SBFLUX',  & 
                  1, iArg, ARGLIST) ) then
        SNCUT_HOST_SBFLUX = ARGLIST(1)(1:MXCHAR_CUTNAME)

      else if ( MATCH_NMLKEY('SNCUT_NOBS_TREST',  & 
                  1, iArg, ARGLIST) ) then
        SNCUT_NOBS_TREST = ARGLIST(1)(1:MXCHAR_CUTNAME)

      else if ( MATCH_NMLKEY('SIMVAR_CUTWIN_STRING',  & 
                  1, iArg, ARGLIST) ) then
        SIMVAR_CUTWIN_STRING = ARGLIST(1)(1:MXCHAR_CUTNAME)

      else if ( MATCH_NMLKEY('CUTWIN_NFILT_SNRMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFILT_SNRMAX(1)
        READ(ARGLIST(2),*) CUTWIN_NFILT_SNRMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_NFILT_SNRMAX2',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFILT_SNRMAX2(1)
        READ(ARGLIST(2),*) CUTWIN_NFILT_SNRMAX2(2)

      else if ( MATCH_NMLKEY('CUTWIN_NFIELD',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_NFIELD(1)
        READ(ARGLIST(2),*) CUTWIN_NFIELD(2)

      else if ( MATCH_NMLKEY('CUTWIN_MWEBV',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_MWEBV(1)
        READ(ARGLIST(2),*) CUTWIN_MWEBV(2)

      else if ( MATCH_NMLKEY('CUTWIN_RESTLAM',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_RESTLAM(1)
        READ(ARGLIST(2),*) CUTWIN_RESTLAM(2)
      else if ( MATCH_NMLKEY('CUTWIN_LAMREST',  &  ! same as CUVTIN_RESTLAM
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_LAMREST(1)
        READ(ARGLIST(2),*) CUTWIN_LAMREST(2)

      else if ( MATCH_NMLKEY('CUTWIN_LAMOBS',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_LAMOBS(1)
        READ(ARGLIST(2),*) CUTWIN_LAMOBS(2)

      else if ( MATCH_NMLKEY('EPCUT_SNRMIN',  & 
                  1, iArg, ARGLIST) ) then
        EPCUT_SNRMIN = ARGLIST(1)(1:MXCHAR_CUTNAME)

      else if ( MATCH_NMLKEY('SNRMAX_FILTERS',  & 
                  1, iArg, ARGLIST) ) then
        SNRMAX_FILTERS = ARGLIST(1)(1:MXFILT_ALL)

      else if ( MATCH_NMLKEY('CUTWIN_SNRMAX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_SNRMAX(1)
        READ(ARGLIST(2),*) CUTWIN_SNRMAX(2)

      else if ( MATCH_NMLKEY('CUTWIN_SNRMAX2',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_SNRMAX2(1)
        READ(ARGLIST(2),*) CUTWIN_SNRMAX2(2)

      else if ( MATCH_NMLKEY('CUTWIN_SNRSUM',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_SNRSUM(1)
        READ(ARGLIST(2),*) CUTWIN_SNRSUM(2)

      else if ( MATCH_NMLKEY('CUTWIN_SNR_NODETECT',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_SNR_NODETECT(1)
        READ(ARGLIST(2),*) CUTWIN_SNR_NODETECT(2)

      else if ( MATCH_NMLKEY('CUTWIN_PSF',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_PSF(1)
        READ(ARGLIST(2),*) CUTWIN_PSF(2)

      else if ( MATCH_NMLKEY('CUTWIN_ZPERR',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_ZPERR(1)
        READ(ARGLIST(2),*) CUTWIN_ZPERR(2)

      else if ( MATCH_NMLKEY('CUTWIN_PHOTPROB',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_PHOTPROB(1)
        READ(ARGLIST(2),*) CUTWIN_PHOTPROB(2)

      else if ( MATCH_NMLKEY('CUTWIN_ZPADU',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_ZPADU(1)
        READ(ARGLIST(2),*) CUTWIN_ZPADU(2)

      else if ( MATCH_NMLKEY('CUTWIN_ZPNPE',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_ZPNPE(1)
        READ(ARGLIST(2),*) CUTWIN_ZPNPE(2)

      else if ( MATCH_NMLKEY('CUTWIN_ERRTEST',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_ERRTEST(1)
        READ(ARGLIST(2),*) CUTWIN_ERRTEST(2)

      else if ( MATCH_NMLKEY('CUTWIN_SIM_PULL',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_SIM_PULL(1)
        READ(ARGLIST(2),*) CUTWIN_SIM_PULL(2)

      else if ( MATCH_NMLKEY('CUTWIN_TREST_TRUEFLUX',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) CUTWIN_TREST_TRUEFLUX(1)
        READ(ARGLIST(2),*) CUTWIN_TREST_TRUEFLUX(2)

! - - - -
       else if ( MATCH_NMLKEY('MAGOBS_SHIFT_PRIMARY',  & 
                   1, iArg, ARGLIST) ) then
         MAGOBS_SHIFT_PRIMARY = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('MAGOBS_SHIFT_ZP',  & 
                   1, iArg, ARGLIST) ) then
         MAGOBS_SHIFT_ZP = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('MAGREST_SHIFT_PRIMARY',  & 
                   1, iArg, ARGLIST) ) then
         MAGREST_SHIFT_PRIMARY = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('MAGREST_SHIFT_ZP',  & 
                   1, iArg, ARGLIST) ) then
         MAGREST_SHIFT_ZP = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FILTER_LAMSHIFT',  & 
                   1, iArg, ARGLIST) ) then
         FILTER_LAMSHIFT = ARGLIST(1)(1:MXCHAR_CUTNAME)
! - - - - - -
      else if ( MATCH_NMLKEY('MAGOBS_SHIFT_PRIMARY_PARAMS',  & 
                  3, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) MAGOBS_SHIFT_PRIMARY_PARAMS(1)
        READ(ARGLIST(2),*) MAGOBS_SHIFT_PRIMARY_PARAMS(2)
        READ(ARGLIST(3),*) MAGOBS_SHIFT_PRIMARY_PARAMS(3)

      else if ( MATCH_NMLKEY('MAGOBS_SHIFT_ZP_PARAMS',  & 
                  3, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) MAGOBS_SHIFT_ZP_PARAMS(1)
        READ(ARGLIST(2),*) MAGOBS_SHIFT_ZP_PARAMS(2)
        READ(ARGLIST(3),*) MAGOBS_SHIFT_ZP_PARAMS(3)
! - - - - - - -
       else if ( MATCH_NMLKEY('FUDGE_FLUXCAL_OFFSET',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_FLUXCAL_OFFSET = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_FLUXCAL_ERROR',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_FLUXCAL_ERROR = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_FLUXCAL_ERRPIX',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_FLUXCAL_ERRPIX = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_FLUXERR_SCALE',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_FLUXERR_SCALE = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_MAG_ERROR',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_MAG_ERROR = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('SIM_FUDGE_MAG_ERROR',  & 
                   1, iArg, ARGLIST) ) then
         SIM_FUDGE_MAG_ERROR = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_MAG_COVERR',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_MAG_COVERR = ARGLIST(1)(1:MXCHAR_CUTNAME)

       else if ( MATCH_NMLKEY('FUDGE_HOSTNOISE_FILE',  & 
                   1, iArg, ARGLIST) ) then
         FUDGE_HOSTNOISE_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('FLUXERRMODEL_FILE',  & 
                   1, iArg, ARGLIST) ) then
         FLUXERRMODEL_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)
       else if ( MATCH_NMLKEY('SIM_FLUXERRMODEL_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SIM_FLUXERRMODEL_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

      else if ( MATCH_NMLKEY('FLUXERRMODEL_OPTMASK',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) FLUXERRMODEL_OPTMASK

       else if ( MATCH_NMLKEY('MAGCOR_FILE',  & 
                   1, iArg, ARGLIST) ) then
         MAGCOR_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

       else if ( MATCH_NMLKEY('SIM_MAGCOR_FILE',  & 
                   1, iArg, ARGLIST) ) then
         SIM_MAGCOR_FILE = ARGLIST(1)(1:MXCHAR_FILENAME)

! - - - - - - -
      else if ( MATCH_NMLKEY('RV_MWCOLORLAW',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) RV_MWCOLORLAW

      else if ( MATCH_NMLKEY('OPT_MWCOLORLAW',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) OPT_MWCOLORLAW

      else if ( MATCH_NMLKEY('PARLIST_MWCOLORLAW',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) PARLIST_MWCOLORLAW(1)  ! just read 1 param for now (Oct 2024)

      else if ( MATCH_NMLKEY('OPT_MWEBV',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) OPT_MWEBV

      else if ( MATCH_NMLKEY('MWEBV_SCALE',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) MWEBV_SCALE

      else if ( MATCH_NMLKEY('MWEBV_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) MWEBV_SHIFT

      else if ( MATCH_NMLKEY('MWEBV_FORCE',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) MWEBV_FORCE

      else if ( MATCH_NMLKEY('REDSHIFT_FINAL_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) REDSHIFT_FINAL_SHIFT

! - - - -
      else if ( MATCH_NMLKEY('HOSTGAL_ZSPEC_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) HOSTGAL_ZSPEC_SHIFT

      else if ( MATCH_NMLKEY('HOSTGAL_ZPHOT_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) HOSTGAL_ZPHOT_SHIFT

! Dec 2024: obsolete/legacy keys
      else if ( MATCH_NMLKEY('HOSTGAL_SPECZ_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) HOSTGAL_ZSPEC_SHIFT

      else if ( MATCH_NMLKEY('HOSTGAL_PHOTOZ_SHIFT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) HOSTGAL_ZPHOT_SHIFT

! - - - -

      else if ( MATCH_NMLKEY('VPEC_ERR_OVERRIDE',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) VPEC_ERR_OVERRIDE

      else if ( MATCH_NMLKEY('FLUXERRCALC_ZPTERR',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) FLUXERRCALC_ZPTERR

      else if ( MATCH_NMLKEY('H0_REF',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) H0_REF

      else if ( MATCH_NMLKEY('OLAM_REF',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) OLAM_REF(1)
        READ(ARGLIST(2),*) OLAM_REF(2)

      else if ( MATCH_NMLKEY('OMAT_REF',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) OMAT_REF(1)
        READ(ARGLIST(2),*) OMAT_REF(2)

      else if ( MATCH_NMLKEY('W0_REF',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) W0_REF(1)
        READ(ARGLIST(2),*) W0_REF(2)

      else if ( MATCH_NMLKEY('WA_REF',  & 
                  2, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) WA_REF(1)
        READ(ARGLIST(2),*) WA_REF(2)

      else if ( MATCH_NMLKEY('ZTOL_HELIO2CMB',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ZTOL_HELIO2CMB
      else if ( MATCH_NMLKEY('zZTOL_HELIO2CMB',  &  ! same, with 'z'
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ZTOL_HELIO2CMB
! - - - -

      else if ( MATCH_NMLKEY('ABORT_ON_BADSURVEY',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_BADSURVEY

      else if ( MATCH_NMLKEY('ABORT_ON_BADFILTER',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_BADFILTER

      else if ( MATCH_NMLKEY('ABORT_ON_NOEPOCHS',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_NOEPOCHS

      else if ( MATCH_NMLKEY('ABORT_ON_NOPKMJD',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_NOPKMJD

      else if ( MATCH_NMLKEY('ABORT_ON_BADAVWARP',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_BADAVWARP

      else if ( MATCH_NMLKEY('ABORT_ON_BADKCOR',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_BADKCOR

      else if ( MATCH_NMLKEY('ABORT_ON_MARGPDF0',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_MARGPDF0

      else if ( MATCH_NMLKEY('ABORT_ON_TRESTCUT',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_TRESTCUT

      else if ( MATCH_NMLKEY('ABORT_ON_DUPLCID',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_DUPLCID

      else if ( MATCH_NMLKEY('ABORT_ON_DUPLMJD',  & 
                  1, iArg, ARGLIST) ) then
        READ(ARGLIST(1),*) ABORT_ON_DUPLMJD

! - - -
      endif

       FOUND_MATCH = FOUND_MATCH_NMLKEY(iArg,ilast,ARG,ARGLIST(1))
       IF ( FOUND_MATCH ) THEN
         DO iuse = ilast, iArg
            USE_LINE_ARGS(iuse) = .TRUE.
         ENDDO
       ENDIF
       iArg = iArg + 1 ;       ilast = iArg

4444  CONTINUE

    DO_STOP = .false.
    IF ( DO_STOP ) then
       print*,' '
       print*,' xxx ------------------------------------- '
       print*,' xxx CUTWIN_SNRMAX = ', CUTWIN_SNRMAX
       print*,' xxx DEBUG STOP '
       STOP
    ENDIF

    RETURN
  END SUBROUTINE SNLCINP_OVERRIDE


! ================================================
    LOGICAL FUNCTION MATCH_NMLKEY(KEYNAME,NARG,iArg,argList)

! Created July 25 2020
! 
! Example 1:
!   KEYNAME           = 'CUTWIN_SNRMAX'
!   LINE_ARGS(iArg)   = 'CUTWIN_SNRMAX'
!   LINE_ARGS(iArg+1) = '5'
!   LINE_ARGS(iArg+2) = '9999'
!   NARG    = 2
!   iArg    = iArg + 2     ! <== returned
!   ARGLIST = '5', '9999'  ! <== returned
! 
! Example 2:
!   KEYNAME           = 'CUTWIN_SNRMAX'
!   LINE_ARGS(iArg)   = 'CUTWIN_SNRMAX=5,9999'
!   NARG    = 2
!   iArg    = iArg         ! <== returned unchanged
!   ARGLIST = '5', '9999'  ! <== returned
! 
! ----------------------
! 
    USE SNPAR
    USE SNCUTS
! CDE,CTRLCOM.
    USE SNLCINP_NML

    IMPLICIT NONE

! function args
! MXCHAR_FILENAME
    CHARACTER KEYNAME*(*)      ! (I) keyname to check
    INTEGER   NARG             ! (I) number of arguments to check
    INTEGER   iArg             ! (I) current index for LINE_ARGS(i)
    CHARACTER ARGLIST(*)*(*)   ! (O) 'NARG' string args for KEYNAME

! local args
    INTEGER LENKEY, LENARG1_EQ, LENARG1_TOT, i
    INTEGER jcomma, jlast
    LOGICAL KEYPLUSARG_ONESTRING, LDMP
    CHARACTER ARG*(MXCHAR_ARG), ARGSTRING*(MXCHAR_ARG)

! ------------- BEGIN ------------

    MATCH_NMLKEY = .false.
    DO i = 1, NARG
       ARGLIST(i) = 'NULL'
    ENDDO
    LDMP = .false.

    LENKEY = INDEX(KEYNAME//' ', ' ') - 1

    CALL GETARG(iArg,ARG)
    LENARG1_EQ  = INDEX(ARG,'=') - 1
    LENARG1_TOT = INDEX(ARG,' ') - 1

    IF ( LENARG1_EQ < 0 ) THEN
       LENARG1_EQ = LENARG1_TOT
       KEYPLUSARG_ONESTRING = .false. ! e.g., CUTWIN_SNRMAX 5 999
    ELSE
       KEYPLUSARG_ONESTRING = .true. ! e.g., CUTWIN_SNRMAX=5,999
    ENDIF

    if ( LDMP ) then
      print*,' xxx --------- MATCH_NMLKEY DUMP ------------- '
      print*,' xxx  LINE_ARG  = ', ARG(1:MXCHAR_PATH)
      print*,' xxx  KEYNAME   = ', KEYNAME(1:LENKEY)
      print*,' xxx  iArg,NARG = ', iArg, NARG
      print*,' xxx  LENARG1   = ', LENARG1_EQ, LENARG1_TOT
      print*,' xxx  LENKEY    = ', LENKEY
      print*,' xxx  KEYPLUSARG_ONESTRING=',KEYPLUSARG_ONESTRING
    endif

    if ( LENARG1_EQ .NE. LENKEY ) RETURN

    MATCH_NMLKEY = ( ARG(1:LENARG1_EQ) .EQ. KEYNAME(1:LENKEY) )

    IF ( .NOT. MATCH_NMLKEY ) RETURN

    IF ( KEYPLUSARG_ONESTRING ) THEN

! xxxx         ARGSTRING = LINE_ARGS(iArg)(LENARG1_EQ+2:LENARG1_TOT)
       ARGSTRING = ARG(LENARG1_EQ+2:LENARG1_TOT)

! break comma-sep ARGSTRING into sub strings
       jlast = 1
       DO i = 1, NARG
         if ( i .EQ. NARG ) then
             ARGLIST(i) = ARGSTRING(jlast:LENARG1_TOT)
          else
             jcomma     = INDEX(ARGSTRING,',') - 1
             ARGLIST(i) = ARGSTRING(jlast:jcomma)
             jlast      = jcomma + 2
          endif
       ENDDO
    ELSE
       DO i = 1, NARG
         iArg = iArg + 1
         CALL GETARG(iArg, ARG)
         ARGLIST(i) = ARG   ! xxx LINE_ARGS(iArg)
       ENDDO
    ENDIF

    if ( LDMP ) then
      do i = 1, NARG
         write(6,38) i, arglist(i)(1:40)
38         format(T10,' xxx out argList(',I2,') = ', A )
      enddo
    endif

    RETURN
  END FUNCTION MATCH_NMLKEY

! ===================================
    LOGICAL FUNCTION FOUND_MATCH_NMLKEY(iArg,ilast,KEY,FIRSTARG)

! Created Aug 2024
! Handle logic to identify valid key+arg using space-separated
! delimeter (e.g., CUTWIN_TREST -20 80) or equal sign with
! no spaces (e.g., CUTWIN_TREST=-20,80).
! In the latter case, passed args should be
! 
!   KEY=CUTWIN_TREST=-20,80 and FIRSTARG=-20
! 
! Move code output SNLCINP_OVERRIDE so that it can be
! used for other NML overrides.


! arguments
    INTEGER iArg, ilast
    CHARACTER KEY*(*), FIRSTARG*(*)

! local args
    LOGICAL FOUND_MATCH_SPSEP, FOUND_MATCH_EQ
    LOGICAL LDMP, VALID_FIRSTARG, HAS_EQUAL

! ----------- BEGIN -----------
    FOUND_MATCH_NMLKEY = .FALSE.
    LDMP = .FALSE.

    VALID_FIRSTARG = .TRUE.
    if ( FIRSTARG .EQ. ''     ) VALID_FIRSTARG = .FALSE.
    if ( FIRSTARG .EQ. 'NULL' ) VALID_FIRSTARG = .FALSE.

    HAS_EQUAL = ( INDEX(KEY,'=') > 0 )
    FOUND_MATCH_SPSEP =  & 
          ( iArg .GT. ilast )   ! space seperated args
    FOUND_MATCH_EQ  =          &  ! key=arg
          ( iArg .EQ. ilast .and. HAS_EQUAL .and. VALID_FIRSTARG )

    FOUND_MATCH_NMLKEY = FOUND_MATCH_SPSEP .or. FOUND_MATCH_EQ

    if ( LDMP ) then
        print*,' xxx - - - - - - - - '
        print*,' xxx  ilast, iArg = ', ilast, iArg
        print*,' xxx  FIRSTARG    = ', FIRSTARG(1:40)
	  print*,' xxx  HAS_EQUAL = ', HAS_EQUAL
        print*,' xxx  FOUND_MATCH[SPSEP,EQ] = ',  & 
                FOUND_MATCH_SPSEP, FOUND_MATCH_EQ
    endif

    RETURN
  END FUNCTION FOUND_MATCH_NMLKEY




! ==================================
    SUBROUTINE PROCSNNML(IERR)
! 
! Created May 17, 2008 by R.Kessler
! 
! Translate namelist strings that depend on filters
! into cutwin_*filt(2,filt) variables.
! 
! Must call after first event is read so that filters are defined.
! Note that RDKCOR is not called yet.
! 
! 
! Aug 31, 2010:
!  - set new variables  NFILT_SNRMAX & IFILT_SNRMAX(i=1,NFILT_SNRMAX)
!  - CUTWIN_SNRMAX_FILT now depends on sparse ifilt_SNRMAX index
! 
! Nov 15, 2011: new variables FUDGE_FLUXCAL_[OFFSET,ERROR]
! 
! Mar 06, 2012: parse new variable SNCUT_HOST_SBFLUX
! Apr 02, 2014: parse SNRMAX_FILTERS, set LFILTDEF_SNRMAX
! Dec 16, 2018: parse SIMVAR_CUTWIN_STRING
! -----------------------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IERR  ! (I) return status

    INTEGER IAFILT(MXFILT_OBS)
    REAL    XAFILT(MXFILT_OBS)

    INTEGER ifiltdef, NF, i, ifilt, ifiltinv, MASK_REQ, L
    CHARACTER cfilt*2
    REAL TMP_ZP, TMP_PRIM
    LOGICAL LDMP, LTMP

! function
    INTEGER FILTINDX

    INTEGER  INIT_SNCUT_NOBS_TREST
    EXTERNAL INIT_SNCUT_NOBS_TREST
! ---------------- BEGIN -------------

    CALL PRBANNER('PROCSNNML')

    IERR = 0
    CALL PARSE_FILTSTRING(1,SNCUT_SNRMAX, NF, iafilt, xafilt )
    NFILT_SNRMAX = NF
    if ( NFILT_SNRMAX .GT. MXFILT_SNRMAX ) then
       write(c1err,661) NFILT_SNRMAX, MXFILT_SNRMAX
661      format(I2,' SNRMAX cuts exceeds bound of MXFILT_SNRMAX=',I2)
       c2err = 'Check &SNLCINP variable SNCUT_SNRMAX'
       CALL  MADABORT("PROCSNNML", C1ERR, C2ERR)
    endif
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       ifilt    = IFILTDEF_INVMAP_SURVEY(ifiltdef)
       cutwin_snrmax_filt(1,i) = xafilt(i) ! Aug 31, 2010
       ifilt_SNRMAX(i)         = ifiltdef  ! Aug 31, 2010
    ENDDO

! -------

    CALL PARSE_FILTSTRING(1,SNCUT_HOST_SBFLUX, NF, iafilt, xafilt )
    NFILT_HOST_SBFLUX = NF
    if ( NF .GT. MXFILT_SNRMAX ) then
       write(c1err,662) NFILT_HOST_SBFLUX, MXFILT_SNRMAX
662      format(I2,' HOST_SBFLUX cuts exceeds bound of ',I2)
       c2err = 'Check &SNLCINP variable SNCUT_SNRMAX'
       CALL  MADABORT("PROCSNNML", C1ERR, C2ERR)
    endif

    DO i = 1, NF
       ifiltdef = IAFILT(i)
       ifilt    = IFILTDEF_INVMAP_SURVEY(ifiltdef)
       cutwin_sbflux_filt(1,i) = xafilt(i)
       ifilt_HOST_SBFLUX(i)    = ifiltdef
    ENDDO

! ----------

    CALL PARSE_FILTSTRING(1,EPCUT_SNRMIN, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       ifilt    = IFILTDEF_INVMAP_SURVEY(ifiltdef)
       cutwin_snrmin_filt(1,ifilt) = xafilt(i)
    ENDDO

! ------
    CALL PARSE_FILTSTRING(0, MAGOBS_SHIFT_PRIMARY, NF,iafilt,xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       magobs_shift_primary_filt(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,MAGOBS_SHIFT_ZP, NF, iafilt, xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       magobs_shift_zp_filt(ifiltdef) = xafilt(i)
       magobs_shift_zp_user(ifiltdef) = xafilt(i)
    ENDDO

! ---
    CALL PARSE_FILTSTRING(1,FUDGE_FLUXCAL_OFFSET, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_FLUXCAL_OFFSET_FILT(ifiltdef) = xafilt(i)
    ENDDO

! - - - -
    CALL PARSE_FILTSTRING(1,FUDGE_FLUXCAL_ERROR, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_FLUXCAL_ERROR_FILT(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,FUDGE_FLUXCAL_ERRPIX, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_FLUXCAL_ERRPIX_FILT(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,FUDGE_FLUXERR_SCALE, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_FLUXERR_SCALE_FILT(ifiltdef) = xafilt(i)
    ENDDO

! - - -

    CALL PARSE_FILTSTRING(1,FUDGE_MAG_ERROR, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_MAG_ERROR_FILT(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,SIM_FUDGE_MAG_ERROR,NF,iafilt,xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       SIM_FUDGE_MAG_ERROR_FILT(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,FUDGE_MAG_COVERR, NF, iafilt, xafilt )
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FUDGE_MAG_COVERR_FILT(ifiltdef) = xafilt(i)
    ENDDO

! --- REST FRAME
    CALL PARSE_FILTSTRING(0, MAGREST_SHIFT_PRIMARY, NF,iafilt,xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       magrest_shift_primary_filt(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,MAGREST_SHIFT_ZP, NF, iafilt, xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       magrest_shift_zp_filt(ifiltdef) = xafilt(i)
       magrest_shift_zp_user(ifiltdef) = xafilt(i)
    ENDDO

    CALL PARSE_FILTSTRING(1,FILTER_LAMSHIFT, NF, iafilt, xafilt)
    DO i = 1, NF
       ifiltdef = IAFILT(i)
       FILTER_LAMSHIFT_FILT(ifiltdef) = xafilt(i)
    ENDDO

! ------
    LDMP = .FALSE.

    do ifiltdef  = 1, MXFILT_ALL
       cfilt     = filtdef_string(ifiltdef:ifiltdef)
       tmp_zp    = MAGOBS_SHIFT_ZP_FILT(ifiltdef)
       tmp_prim  = MAGOBS_SHIFT_PRIMARY_FILT(ifiltdef)

       LTMP = ( tmp_zp .ne. 0.0 .or. tmp_prim .NE. 0.0 )
       if ( LDMP .and. LTMP ) then
         print*, 'DBUG: ', cfilt, ifiltdef  & 
            , ' : MAGSHIFT_[ZP,PRIMARY] = ', tmp_zp, tmp_prim
       endif
    enddo

! -------
    NF = INDEX(SNRMAX_FILTERS,' ') - 1

    IF ( NF > 0 ) THEN

      do ifilt = 1, MXFILT_ALL
         LFILTDEF_SNRMAX(ifilt) = .FALSE.
      enddo
      do ifilt = 1, NF
         cfilt = SNRMAX_FILTERS(ifilt:ifilt)
         IFILTDEF = FILTINDX(cfilt)
         ifiltinv = IFILTDEF_INVMAP_SURVEY(ifiltdef)
         if ( ifiltinv < 1 ) then
            C1ERR = 'Invalid filter = ' // Cfilt
            C2err = 'Check &SNLCINP variable SNRMAX_FILTERS = '  & 
                         // SNRMAX_FILTERS(1:NF)
            CALL  MADABORT("PROCSNNML", C1ERR, C2ERR)
         endif
         LFILTDEF_SNRMAX(ifiltdef) = .TRUE.
      enddo

    ENDIF


! ------------------------------------------------------------
! Dec 2018: check for cuts on SIMVAR (e.g., SIM_c, SIM_ZCMB)

    CALL PARSE_SIMVAR_CUTS()

! Jan 12 2021: check option to fix sign VPEC for older data
    DOFIX_WRONG_SIGN_VPEC = (.not. RESTORE_WRONG_VPEC ) .and.  & 
                              (.not. ISCORRECT_SIGN_VPEC)
    if ( DOFIX_WRONG_SIGN_VPEC ) then
       print*,' '
       print*,'    *** ALERT: VPEC sign will be flipped *** '
       print*,' '
    endif

    if ( RESTORE_WRONG_VPEC ) then
       print*,' '
       print*,'    *** ALERT: use wrong VPEC formula *** '
       print*,' '
    endif

       CALL FLUSH(6)

! set flag for systematic zshift (Jun 2021)
    CALL CHECK_LEGACY_INPUT(HOSTGAL_SPECZ_SHIFT,HOSTGAL_ZSPEC_SHIFT)
    CALL CHECK_LEGACY_INPUT(HOSTGAL_PHOTOZ_SHIFT,HOSTGAL_ZPHOT_SHIFT)

    DOzSHIFT =  & 
          ( REDSHIFT_FINAL_SHIFT .NE. 0.0 ) .or.  & 
          ( HOSTGAL_ZSPEC_SHIFT  .NE. 0.0 ) .or.  & 
          ( HOSTGAL_ZPHOT_SHIFT  .NE. 0.0 )

! -------------------
! check for Nobs cut in multiple Trest ranges

    L = LEN_TRIM(SNCUT_NOBS_TREST) + 1
    MASK_REQ=INIT_SNCUT_NOBS_TREST(SNCUT_NOBS_TREST(1:L)//char(0),L)

! set cut window to require MASK in which all trest ranges are satisfied;
! this feature also enables checking which Trest range(s) fail.

    cutwin_mask_nobs_trest(1) = float(MASK_REQ) - 0.5
    cutwin_mask_nobs_trest(2) = float(MASK_REQ) + 0.5

    RETURN
  END SUBROUTINE PROCSNNML

! ========================================
    SUBROUTINE CHECK_LEGACY_INPUT(VAL_LEGACY, VAL_NOMINAL)

! Created Dec 29 2024
! If the legacy value has been set, transfer it to the VAL_NOMINAL


    USE SNPAR

    IMPLICIT NONE

    REAL VAL_LEGACY   ! Input
    REAL VAL_NOMINAL  ! Output

    IF(VAL_LEGACY < .99*LEGACY_INIT_VAL ) THEN
       VAL_NOMINAL = VAL_LEGACY
    ENDIF

    RETURN
  END SUBROUTINE CHECK_LEGACY_INPUT

! ==================================
    SUBROUTINE RDSNIGNORE
! 
! Nov 10, 2008
! Read file SNCID_IGNORE_FILE and load logical array LSNCID_IGNORE
! This is designed for the thousands of SDSS-II candidates,
! so that large numbers of candidates can be ignored from a
! list file.
! 
! Sep 2 2019: allow up to MXIGNORE_LIST=500 (instead of 52)
! Sep 10 2021:
!   + fix to work properly with/without DOCANA
!   + ignore comments
!   + print first/last CID for visual crosscheck
! 
! --------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER CCID*(MXCHAR_CCID), cFILE*(MXCHAR_FILENAME), FNAM*12
    CHARACTER WD1*40
    INTEGER MSKOPT, i, LENFILE, LENCID, NWD, FLAG_DOCANA
    INTEGER NSTORE, NCCID_IGNORE_ORIG
    LOGICAL FOUND_DOCANA, FOUND_DOCANA_END, LREQUIRE_DOCANA

    INTEGER  STORE_PARSE_WORDS, ISTAT_FILE_DOCANA
    EXTERNAL STORE_PARSE_WORDS, REACT_MISSING_DOCANA
! ------------- BEGIN ---------

    IF ( SNCID_IGNORE_FILE .EQ. ' ' ) RETURN

    FNAM = 'RDSNIGNORE'
    CALL ENVreplace(SNCID_IGNORE_FILE)
    LENFILE   = INDEX(SNCID_IGNORE_FILE,' ') - 1
    cFILE     = SNCID_IGNORE_FILE(1:LENFILE)//char(0)
    MSKOPT    = MSKOPT_PARSE_WORDS_FILE  & 
                + MSKOPT_PARSE_WORDS_IGNORECOMMENT
    NWD = STORE_PARSE_WORDS(MSKOPT, cFILE, FNAM//char(0), LENFILE, 12)

    LENCID = INDEX(CCID,' ')-1

! check if 1st word is DOCANA key
    CALL get_PARSE_WORD_fortran(1, WD1, LENCID)
    FOUND_DOCANA     = ( WD1 .EQ. 'DOCUMENTATION:' )
    FOUND_DOCANA_END = (.FALSE.)
    LREQUIRE_DOCANA  = (REQUIRE_DOCANA > 0 )

    if ( .NOT. FOUND_DOCANA ) then
       CALL REACT_MISSING_DOCANA(LREQUIRE_DOCANA,cFILE,LENFILE)
    endif

    NSTORE = 0
    NCCID_IGNORE_ORIG = NCCID_IGNORE !might include CIDs from elsewhere

! - - - - - -
    DO 200 i = 1, NWD
       CALL get_PARSE_WORD_fortran(i, CCID, LENCID)

       IF ( FOUND_DOCANA .and. .NOT. FOUND_DOCANA_END) THEN
          FLAG_DOCANA = ISTAT_FILE_DOCANA(i,CCID)
          IF ( FLAG_DOCANA==FLAG_DOCANA_END ) FOUND_DOCANA_END=.TRUE.
          goto 200
       ENDIF

       NSTORE       = NSTORE + 1
       NCCID_IGNORE = NCCID_IGNORE + 1
       SNCCID_IGNORE_ALL(NCCID_IGNORE) = CCID

       IF ( NCCID_IGNORE > MXIGNORE_LIST ) then
         write(c1err,661 ) NCCID_IGNORE
         write(c2err,662 ) MXIGNORE_LIST
         CALL  MADABORT(FNAM, C1ERR, C2ERR)
       endif
661        format('NCCID_IGNORE=',I4,' exceeds bound.')
662        format('Check limit MXIGNORE_LIST = ', I4 )

200   CONTINUE

    CALL PRBANNER(FNAM)
    print*,' Loaded ', NCCID_IGNORE,' IGNORE-candidates from '
    print*, '  ', SNCID_IGNORE_FILE(1:LENFILE)
    print*,'  First IGNORED CID: ',  & 
           SNCCID_IGNORE_ALL(NCCID_IGNORE_ORIG+1)
    if ( NSTORE > 1 ) then
      print*,'  Last  IGNORED CID: ',  & 
           SNCCID_IGNORE_ALL(NCCID_IGNORE)
    endif
    print*,' '

    RETURN
  END SUBROUTINE RDSNIGNORE


! ====================================
    INTEGER FUNCTION ISTAT_FILE_DOCANA(iwd,string)

! use this utility for files where DOCANA keys break the
! file structure (e.g., SNCID_IGNORE_FILE). This utility
! easily enables skipping the DOCANA stuff.
! 
!   Returns  1 if STRING = 'DOCUMENTATION:'
!   Returns  2 if STRING = 'DOCUMENTATION_END:'
!   Returns -1 if first word is not 'DOCUMENTATION:'
!        (this flags calling function to abort)
! 
! - - - - -


    USE SNPAR

    IMPLICIT NONE

    INTEGER   IWD          ! (I) word number in file
    CHARACTER STRING*(*)   ! (I)  string to check for DOCANA keys
! ---------- BEGIN -----------

    ISTAT_FILE_DOCANA = 0
    if ( STRING .EQ. 'DOCUMENTATION:'    ) then
       ISTAT_FILE_DOCANA = FLAG_DOCANA_START
     endif
     if ( STRING .EQ. 'DOCUMENTATION_END:'    ) then
        ISTAT_FILE_DOCANA = FLAG_DOCANA_END
     endif

    if (IWD==1 .and. ISTAT_FILE_DOCANA.NE.FLAG_DOCANA_START) then
       ISTAT_FILE_DOCANA = FLAG_DOCANA_ERROR
    endif

    RETURN
  END FUNCTION ISTAT_FILE_DOCANA

! ====================================
    SUBROUTINE PRBANNER ( banner )
    USE SNPAR
    USE CTRLCOM
    character banner*(*)

    if ( .NOT. REDUCE_STDOUT_BATCH ) print*,' '
    print*,  & 
         ' ######################################################## '
    print*,'   ', BANNER

    print*,  & 
         ' ######################################################## '

    if ( .NOT. REDUCE_STDOUT_BATCH ) print*,' '
    CALL FLUSH(6)

    RETURN
  END SUBROUTINE PRBANNER 

! ========================================
    SUBROUTINE PRINT_SNSTATS
! 
! Print SN stats : reaad and passing cuts.
! If no SN pass cuts, list cut(s) that always fail.
! 
! Get <z> and RMS from hid 11
! 
! Apr 26, 2017: I6 -> I7 in format 48.


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

    integer icut, i,  JDIFF
    LOGICAL LTEST, LFAIL

    REAL Tsn, Ttot, f_dupl
    EXTERNAL PRINT_CPUTIME

! ------------ BEGIN --------------

    global_banner = 'Job Summary '
    CALL PRBANNER ( GLOBAL_BANNER )

! ----------------------------------------------
! Start with Npasscut vs. ICUT, and for each type.
    CALL PRINT_NPASSCUTS()

! Aug 7 2020: for batch job, write stats to separate YAML file for easy parsing
    CALL PRINT_JOBSPLIT_OUT()

! ------------------------
    write(6,48) N_SNLC_CUTS, N_SNLC_PROC
48    format(/,T5,'Finished processing ',I7,  & 
         ' SN after snana  cuts',  & 
          3x,'(',I7,' before cuts)'  )
    CALL FLUSH(6)

! Oct 21, 2010: skip other "NEVER FOUND" message if CID cut always fails
    IF ( NACCEPT_CID .EQ. 0 ) GOTO 600

! check for valid SNTYPES (Jul 2012)
    LFAIL = NSNTYPE_LIST .GT. 0 .and. NACCEPT_TYPE .LE. 0
    IF ( LFAIL ) THEN
       print*,' '
       print*,' WARNING-SNTYPE: ',  & 
             'NEVER FOUND A VALID SNTYPE ?!?!?! '
       print*,' WARNING-SNTYPE: Check nml variable SNTYPE_LIST = ',  & 
              SNTYPE_LIST(1), SNTYPE_LIST(2)
       print*,' '
    ENDIF


#if defined(CCDLATER)
! Need to define NACCEPT_DETNUM, and increment it.
! check for valid DETNUM
    LFAIL = NDETNUM_LIST .GT. 0 .and. NACCEPT_DETNUM .LE. 0
    IF ( LFAIL ) THEN
       print*,' '
       print*,' WARNING-DETNUM: ',  & 
             'NEVER FOUND A VALID DETNUM ?!?!?! '
       print*,' WARNING-DETNUM: Check nml variable DETNUM_LIST = ',  & 
              DETNUM_LIST(1), DETNUM_LIST(2)
       print*,' '
    ENDIF
#endif

! print fit-cut stats
    if ( N_SNLC_FIT .GT. 0 ) then
      write(6,49) N_SNLC_FITCUTS
49      format(T5,'Finished processing ',I6,  & 
         ' SN after fitter cuts' )
    endif

! if no supernova pass cuts, look for cuts that always fail.

       print*,' '

600   CONTINUE

    IF ( NACCEPT_CID .EQ. 0 ) THEN
       write(6,66)  CUTBIT_CID, 'CID'
       print*,'     check CUTWIN_CID and SNCID_LIST in &SNLCINP'

    ELSE IF  ( N_SNLC_CUTS .EQ. 0 ) THEN
      DO icut = 1, CUTBIT_MJD_MARKER
        LTEST = BTEST ( CUTMASK8_SN_ALL , icut-1 )
        if ( LTEST .and. NACCEPT_CUT(icut) .EQ. 0 ) then
           i = INDEX(cutvar_name(icut), ':' )
           write(6,66) icut, cutvar_name(icut)(1:i-1)
        endif
      ENDDO
    ENDIF
66    FORMAT(T5,'** WARNING ** EVERY SN FAILS CUTBIT ', I3,  & 
             ' :  ', A )

! -----------------------------------

!       write(6,170) NEPOCH_BADPHOT_SUM, NEPOCH_USE
    write(6,170) NEPOCH_BADPHOT_SUM, NEPOCH_TOT
170   format(T5,'Bad PHOTOMETRY FLAG for ',I8,' of ', I9,' epochs')

    print*,' '
! -------------------------------------------------
! print out processing time, and process time per SN.


    CALL PRINT_CPUTIME(JTIME_LOOPSTART,  & 
              "CPUTIME_PROCESS_ALL"//char(0),  & 
              "minute"//char(0), 0,    20,20)

    CALL PRINT_CPUTIME(JTIME_LOOPSTART,  & 
              "CPUTIME_PROCESS_RATE"//char(0),  & 
              "second"//char(0), N_SNLC_CUTS,    20,20)

! --------------------------
! print DUPLICATE-CID WARNING

    IF ( N_DUPLICATE_CID > 0 ) THEN
       f_dupl = float(N_DUPLICATE_CID)/ float(N_SNLC_PROC)
       write(6,60)
       write(6,60)
       write(6,606) N_DUPLICATE_CID, 100.*f_dupl
       write(6,60)
       write(6,60)
 606     format(/, T5,'SEVERE WARNING: rejected ',I4,' Duplicate CIDs',  & 
                3x,'(',F6.2,'%)', / )
    ENDIF

60      format(T2, 32('@-') )

! print DUPLICATE-MJD WARNING
    IF ( N_DUPLICATE_MJD > 0 ) THEN
       f_dupl = float(N_DUPLICATE_MJD)/ float(NEPOCH_USE)
       write(6,60)
       write(6,60)
       write(6,616) N_DUPLICATE_MJD, 100.*f_dupl
       write(6,617) NSTORE_DUPLICATE_MJD
       write(6,60)
       write(6,60)
       call flush(6)
 616     format(/,T5,'SEVERE WARNING: FOUND ',I4,' Duplicate MJD+BAND',  & 
                3x,'(',F6.3,'%)' )
 617     format(T20,I4,' Duplicate MJDs are unique', /)
    ENDIF

    IF ( OPT_TABLE(ITABLE_OUTLIER) > 0 ) THEN
       CALL PRINT_OUTLIER_SUMMARY()  ! Mar 2021
    ENDIF

    RETURN
  END SUBROUTINE PRINT_SNSTATS


! =======================
    SUBROUTINE PRINT_NPASSCUTS()

! 
! Jan 3, 2013
! 
! Print Number of SN passing each incremental cut, and for each type.
! Table has up to 6 types, if more then make new table so that
! table does not exceed 80 columns.
! 
! Feb 22, 2013: fix to work for TYPE=0
! 
! For ITYPE=-1, write string 'ALL' to be more clear about meaning.
! 
! Jun 19 2019: incluce ROWKEY = "CUTSTAT:" to make log file parsable.
! 

    USE SNPAR
    USE SNCUTS
    USE CTRLCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    LOGICAL USETYPE(-1:MXTYPE), LPRIN

    INTEGER  & 
         ITYPE, NTYPE, ITYPE_LIST(MXTYPE)  & 
        ,I1, I2, i, icut  & 
        ,NTYPE_PER_TABLE  & 
        ,N, NLAST, jcolon

    CHARACTER CUTNAME*28, ROWKEY*8

! ------------- BEGIN -------------

! first get list of used TYPES
    NTYPE = 0
    DO itype = -1, MXTYPE
      USETYPE(itype) = .FALSE.
      if (  NPASSCUT_INCREMENT(itype,1) .GT. 0 )  then
         USETYPE(itype) = .TRUE.
         NTYPE =  NTYPE + 1
         ITYPE_LIST(NTYPE) = ITYPE
      endif
    ENDDO

! print table header

    ROWKEY = "CUTSTAT:"

    print*,'                         ',  & 
            'Number of SN passing incremental cut for'

    NTYPE_PER_TABLE = 6
    DO  777 I1 = 1, NTYPE, NTYPE_PER_TABLE

      write(6,20) ROWKEY, 'CUT-NAME      ITYPE='
20      format(A, 2x, A20, $)

      I2   = MIN(I1+NTYPE_PER_TABLE-1,NTYPE)
      DO i = I1, I2
         itype = ITYPE_LIST(i)

         if ( itype .EQ. -1 ) then
            write(6,27) 'ALL'
         else
            write(6,28) ITYPE
         endif
 27        format(A8,$)
 28        format(I8,$)
      ENDDO
      print*,' '
      write(6,80)
80      format(T2,  72('-') )


      DO 888 ICUT = 1, NCUTBIT_SNLC

         if(.NOT.LSIM_SNANA .and. icut.EQ.CUTBIT_SEARCH ) goto 888

! print line if any entry has a change in NPASSCUT
         LPRIN = .FALSE.
         IF ( ICUT .EQ. 1 ) THEN
            LPRIN = .TRUE.
         ELSE
            DO i = I1, I2
              itype = ITYPE_LIST(i)
              NLAST = NPASSCUT_INCREMENT(itype,icut-1)
              N     = NPASSCUT_INCREMENT(itype,icut)
              if ( N .LT. NLAST ) LPRIN = .TRUE.
            ENDDO
         ENDIF
         IF ( LPRIN ) THEN
            jcolon    = INDEX(cutvar_name(icut), ':' )
            CUTNAME   = CUTVAR_NAME(icut)(1:jcolon-1) ! remove ':'
            write(6,20) ROWKEY, CUTNAME
            DO i = I1, I2
              itype = ITYPE_LIST(i)
              write(6,28) NPASSCUT_INCREMENT(itype,icut)
            ENDDO
            print*,' '
         ENDIF
888     CONTINUE

! for fit, write number passing fit cuts.

      IF ( NFIT_ITERATION > 0 ) THEN
          CUTNAME = 'FIT+CUTS'
          write(6,20) ROWKEY, CUTNAME
          DO i = I1, I2
             itype = ITYPE_LIST(i)
             write(6,28) NPASSCUT_FIT(itype)
          ENDDO
      ENDIF
            print*,' '
            print*,' '
777   CONTINUE

    write(6,80)
    RETURN
  END SUBROUTINE PRINT_NPASSCUTS

! =========================
    SUBROUTINE PRINT_CPU_REMAIN(ISN)
! 
! Created Apr 4 2024
! Print fraction of total events processed, elapsed time,
! and estimate time remaining:
! 
!   Processed ~12% of events in mm minutes -> mm minutes remaining.
! 
! Used only for BATCH mode.
! ----------------



    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER ISN  ! (I) integer index (1-NSN)

    INTEGER NSN_TOT, JTIME, JTDIF, PERCENT
    REAL FRAC, T_ELAPSE, T_REMAIN
! ------------ BEGIN ---------

    if ( .NOT. STDOUT_UPDATE ) RETURN

! bail for subset option for which this update makes no sense
    if ( .NOT. ISJOB_BATCH     ) RETURN
    if ( NCID_LIST > 0         ) RETURN
    if ( MXEVT_PROCESS < 10000 ) RETURN

    NSN_TOT = N_SNLC_READ(1)
    FRAC = FLOAT(ISN) / FLOAT(NSN_TOT)
    PERCENT = INT(FRAC*100.0 + 0.5)

    JTIME = TIME()
    JTDIF = JTIME - JTIME_LOOPSTART
    T_ELAPSE = FLOAT(JTDIF)/60.0
    T_REMAIN = T_ELAPSE * ( 1.0/FRAC - 1.0 )  ! prediction

    if ( T_ELAPSE < 1.0 ) RETURN  ! avoid printing crazy predictions

    write(6,20) PERCENT, T_ELAPSE, T_REMAIN
20    format(T6,'Processed ', I3,'% of events in ',  & 
           F6.1, ' minutes -> ', F6.1,' minutes remaining.')
    call flush(6)

    RETURN
  END SUBROUTINE PRINT_CPU_REMAIN

! ==================================
    SUBROUTINE PRINT_JOBSPLIT_OUT()

! Created Aug 7 2020
! If snana job is submitted by batch job (i.e., JOBSPLIT is set),
! create supplmental output file
!    [VERSION].YAML
! and write stats for batch process to read.
! This YAML file avoids submit_batch script having to parse
! potentially large log files. Also note that ABORT_IF_ZERO
! is a standard key telling batch script to ABORT if its
! value is zero; this avoids the batch script having to know
! program-specific keys for NEVT_what?
! 
! Oct 12 2020: check OPT_YAML
! Jul 08 2021: write N_SNHOST_ZSPEC[ZPHOT]


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT
    USE FILTCOM

    IMPLICIT NONE

    INTEGER LEN, JDIFF, AIZ
    REAL    T_CPU
    CHARACTER OUTFILE*(MXCHAR_FILENAME)

! ------------- BEGIN -------------

! if user does NOT request YAML file, then check default
! to create YAML file only if this is a batch job.

    IF ( OPT_YAML <= 0 ) THEN
      IF ( .NOT. ISJOB_BATCH       ) RETURN
      IF ( TEXTFILE_PREFIX .EQ. '' ) RETURN
    ENDIF

! write [REPFIX].YAML file...

    LEN     = INDEX(TEXTFILE_PREFIX,' ' ) - 1
    OUTFILE = TEXTFILE_PREFIX(1:LEN) // '.YAML'

! if TEXTFILE_PREFIX is not defined, then use VERSION_PHOTOMETRY
    if ( LEN == 0 ) then
       LEN = INDEX(VERSION_PHOTOMETRY(1),' ') - 1
       OUTFILE = VERSION_PHOTOMETRY(1)(1:LEN) // '.YAML'
    endif

    write(6,40) OUTFILE(1:LEN+5)
40    format(' Write stats to YAML output: ', A )

    OPEN(   UNIT   = LUNDAT  & 
            , FILE   = OUTFILE  & 
            , STATUS = 'UNKNOWN'  & 
                 )

    JDIFF = JTIME_LOOPEND - JTIME_LOOPSTART
    T_CPU = float(JDIFF)/60.

    write(LUNDAT,19) 'SURVEY:          ', SURVEY_NAME
    write(LUNDAT,19) 'FILTERS:         ',  & 
             SURVEY_FILTERS(1:NFILTDEF_SURVEY)
    write(LUNDAT,20) 'IDSURVEY:        ', IDSURVEY
    write(LUNDAT,20) 'NEVT_TOT:        ', N_SNLC_PROC
    write(LUNDAT,20) 'NEVT_LC_CUTS:    ', N_SNLC_CUTS

    write(LUNDAT,20) 'NEVT_HOST_ZSPEC: ', N_SNHOST_ZSPEC
    write(LUNDAT,20) 'NEVT_HOST_ZPHOT: ', N_SNHOST_ZPHOT
    write(LUNDAT,20) 'NEVT_SPECTRA: ',    N_SNLC_SPEC ! nevt with spectra

    write(LUNDAT,20) 'NEVT_LCFIT_CUTS: ', N_SNLC_FITCUTS

! write NEVT per redshift source
    CALL PRINT_JOBSPLIT_zSRC(LUNDAT, "LC_CUTS")
    CALL PRINT_JOBSPLIT_zSRC(LUNDAT, "LCFIT_CUTS")

    write(LUNDAT,21) 'CPU_MINUTES:     ', T_CPU
19    format(A, A)
20    format(A, I8)
21    format(A, F8.2)


! check options for ABORT_IF_ZERO (AIZ)

#if defined(SNANA)
    AIZ = N_SNLC_CUTS
#elif defined(SNFIT)
    AIZ = N_SNLC_FITCUTS
#elif defined(PSNID)
    AIZ = N_SNLC_FITCUTS
#endif
    IF ( MXEVT_CUTS < 10 ) AIZ = N_SNLC_CUTS


    write(LUNDAT,20) 'ABORT_IF_ZERO:   ', AIZ

    CLOSE ( UNIT = LUNDAT )

    RETURN
  END SUBROUTINE PRINT_JOBSPLIT_OUT

! ==================================
    SUBROUTINE PRINT_JOBSPLIT_zSRC(LUN, WHICH)

! Created Nov 26 2024
! write NEVT per redshift source to YAML file;
! intended to help figure out what redshift sources are
! in the data or sim.


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT
    USE FILTCOM

    IMPLICIT NONE

    INTEGER LUN  ! (I)
    CHARACTER WHICH*(*)  ! (I)  LC_CUTS or LCFIT_CUTS

    INTEGER NEVT, NEVT2, MASK, MASK2, NEVT_or
    INTEGER N_MASK_LOCAL(0:MXMASK_zSOURCE)
    CHARACTER KEY*28

! ----------  BEGIN ------------

    IF ( WHICH(1:7) .EQ. 'LC_CUTS' ) THEN
       KEY = 'MASK_zSOURCE_LC_CUTS:'
	 DO mask = 0, MXMASK_zSOURCE
	     N_MASK_LOCAL(mask) = N_MASK_zSOURCE_LC_CUTS(mask)
	 ENDDO
    ELSE
       KEY = 'MASK_zSOURCE_LCFIT_CUTS:'
	 DO mask = 0, MXMASK_zSOURCE
	     N_MASK_LOCAL(mask) = N_MASK_zSOURCE_LCFIT_CUTS(mask)
	 ENDDO
    ENDIF

    write(LUN,28) ' '
    write(LUN,28) KEY // '   # redshift sources for ' // WHICH
    write(LUN,38) '   += 1 -> zSPEC_HOST'
    write(LUN,38) '   += 2 -> zSPEC_SN'
    write(LUN,38) '   += 4 -> zPHOT_HOST(mean+stdev or quantiles)'
    write(LUN,38) '   += 8 -> zPHOT_HOST(quantiles)'
    write(LUN,38)  & 
         '    MASK:     NEVT(MASK=val)    NEVT(MASK val)'
    DO 200 MASK = 0, MXMASK_zSOURCE
       NEVT    = N_MASK_LOCAL(MASK)
       NEVT_or = 0
       DO 202 MASK2 = 0, MXMASK_zSOURCE
	    NEVT2 = N_MASK_LOCAL(MASK2)
	    if ( NEVT2 == 0 ) GOTO 202
	    if ( IAND(MASK,MASK2) == MASK ) then
	        NEVT_or = NEVT_or + NEVT2
          endif
202	 CONTINUE

	 if ( NEVT > 0 ) then
          write(LUNDAT,212) MASK, NEVT, NEVT_or
212	    format('  - ', I6,':', I10, 8x, I10 )
       endif
200   CONTINUE
    write(LUNDAT,28) ' '

28    format(A)
38    format('# ', A)

    RETURN
  END SUBROUTINE PRINT_JOBSPLIT_zSRC


! ======================================
    SUBROUTINE SNANA_DRIVER(ISN_ALL, ISN_PROC, IVERS)

! Created May 18, 2012 by R.Kessler
! 
! Analysis driver for SN with
!   ISN_ALL  = absolute index for all events, regardless of cuts
!   ISN_PROC = sparse index of processed events (after cuts)
!   IVERS    = index of VERSION_PHOTOMETRY (added Oct 2025)
! 
! Move code from MAIN to here so that each SN is fit
! right after it is read.
! 
! Jun 21 2018:
!   if ERRFLAG_FIT > 0 on any iteraton, REJECT fit immediately
!   See REJECT_FIT logical.
! 
! Dec 19 2024: print CPU time per event in stdout update
! ----------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER ISN_ALL, ISN_PROC, IVERS  ! (I)

#if defined(SNFIT)
    INTEGER NCALL_MNFIT, ITER
    LOGICAL LAST_ITER, REJECT_FIT
#endif

    INTEGER IERR, i
    REAL*8  PS8
    REAL t_start, t_end  ! Dec 2024
    LOGICAL REJECT_PRESCALE, USE_MINOS_LOCAL
    CHARACTER FNAM*14

! ----------------- BEGIN -------------

    FNAM = 'SNANA_DRIVER'

    if ( STDOUT_UPDATE ) CALL CPU_TIME(t_start)

    if ( ISN_PROC > MXSNLC-1 ) THEN
      write(c1err,666) N_SNLC_PROC
666     format('ISN_PROC=',I8,' exceeds MXSNLC array bound.')
      c2err = 'Try splitting job into multiple jobs with CID-ranges'
      CALL MADABORT(FNAM, c1err, c2err)
    endif

    if ( LSIM_SNANA ) THEN
       PS8 = DBLE(SIM_PRESCALE)
       if ( REJECT_PRESCALE(N_SNLC_PROC,PS8) ) RETURN
    endif

! - - - - - - - - - - -
    NCALL_SNANA_DRIVER  = NCALL_SNANA_DRIVER  + 1

! one-time init after first event is read.

    IF ( NCALL_SNANA_DRIVER .EQ. 1 )  THEN

      CALL MON_SNANA(IFLAG_INI)     ! book 1d distributions

#if defined(SNFIT)
      CALL FITPAR_INI2(IERR)
#endif

#if defined(PSNID)
      CALL PSNIDINI2(IERR)
#endif


        if ( IERR .NE. 0 ) then
          c1err = 'Problem with 2nd init when NCALL_SNANA_DRIVER=1'
          c2err = 'CID = ' // SNLC_CCID
          CALL  MADABORT(FNAM, C1ERR, C2ERR )
        endif
    ENDIF


#if defined(SNFIT)
! check for SN-dependent filter response
    CALL FILTER_UPDATE_DRIVER()
#endif

! compute quantities, apply AB offsets, set cut-mask
    CALL SNRECON()       ! compute variables, apply AB offsets

    CALL U3BAND()        ! special U-band test

    CALL SET_CUTMASK()   ! set cutmask

#if defined(SNFIT)
! check option to include IDEAL fit, after SNRECON (Aug 29 2017)
    CALL SIMFIT_IDEAL_PREP(0)
#endif

! set default to use minimized fitpar unless user calls MARG_DRIVER
    USEPDF_MARG = .FALSE.

! check dump-option for SN that fail cuts.

    IF ( .not. LSNCUTS ) THEN
       if ( STDOUT_UPDATE ) write(6,41) SNLC_CCID
41       format(T20,'SNANA cuts REJECT CID = ',A8, 2x,'=> SKIP ')
	 CALL FLUSH(6)
       IF ( LDMP_SNFAIL ) CALL DMP_SNFAIL()
       CALL MON_SNANA(IFLAG_ANA)       ! update tables for rejected SN
       RETURN  ! skip this SN
    ENDIF

    IF ( STDOUT_UPDATE ) THEN
      write(6,42) isn_proc  & 
              , SNLC_CCID  & 
              , ISNLC_NEWMJD_STORE  & 
              , ISNLC_NEWMJD_CUTS
    ENDIF

    CALL EXEC_REFORMAT(IVERS)

#if defined(SNFIT)

! do the fit

     ITER          = 0
     LREPEAT_ITER  = .FALSE.
     LREPEAT_MINOS = .FALSE.
     REJECT_FIT    = .FALSE.
     NCALL_MNFIT   = 0

     N_SNLC_FIT = N_SNLC_FIT + 1
     DO i = 1, FCNFLAG_MAX
       NCALL_FCNFLAG(i) = 0
     ENDDO

     DO 410 WHILE( iter < NFIT_ITERATION .and. .NOT.REJECT_FIT)

        ITER = ITER + 1
        LAST_ITER = ( ITER .EQ. NFIT_ITERATION )

        IF ( ITER==1 ) THEN
            USE_MINOS_LOCAL = .FALSE.
        ELSE
            USE_MINOS_LOCAL = (USE_MINOS .or. LREPEAT_MINOS)
        ENDIF

        CALL FITPAR_PREP ( iter, IERR )  ! init fit params
        if ( IERR .NE. 0 ) then
          ERRFLAG_FIT = IERR
 	    REJECT_FIT = .TRUE.
          write(6,40) ITER, SNLC_CCID(1:ISNLC_LENCCID)
40          format(T5,'FITPAR_PREP PROBLEM => SKIP FIT-ITER=',  & 
                    I2,'  for  CID=',A )
          GOTO 410          ! bad init => skip fit-iteration
        endif

        if ( STDOUT_UPDATE ) THEN
          WRITE(GLOBAL_BANNER,44) ISN_PROC,  & 
               SNLC_CCID(1:ISNLC_LENCCID), NFITPAR_MN, ITER
          CALL PRBANNER(GLOBAL_BANNER)
44          format(T5,'CALL MNFIT_DRIVER for ISN_PROC=',I5,  & 
                 3x,'CID=',A, 3x, 'NFITPAR=',I3, 3x,'iter=',I1 )
        endif

! call main driver for fitter.
! Note that all returned output from MNFIT_DRIVER goes to common blocks.

        CALL MNFIT_DRIVER (  & 
             SNLC_CCID, NFITPAR_MN                 &  ! (I)
            ,INIVAL, INISTP, INIBND                &  ! (I)
            ,PARNAME_STORE                         &  ! (I)
            ,USE_MINOS_LOCAL                       &  ! (I)
            ,MINUIT_PRINT_LEVEL                    &  ! (I)
            ,FITVAL(1,iter)                        &  ! (O)
            ,FITERR_PLUS(1,iter)                   &  ! (O)
            ,FITERR_MINUS(1,iter)                  &  ! (O)
            ,FITCHI2_MIN                           &  ! (O)
            ,NFIXPAR, ERRTYPE                      &  ! (O)
            ,MNSTAT_COV                            &  ! (O)
            ,IERR                                  &  ! (O)
                  )

! bail on error
       IF ( IERR > 90 ) THEN  ! May 2024
	   ERRFLAG_FIT = IERR
	   REJECT_FIT = .TRUE.
         GOTO 410
       ENDIF

! avoid infinite loop with MNFIT calls ...
       NCALL_MNFIT = NCALL_MNFIT + 1
       IF  ( NCALL_MNFIT .GT. 2*MXITER ) THEN
          write(c1err,641) 2*MXITER, SNLC_CCID
          write(c2err,642) NFIT_ITERATION
641         format('MNFIT_DRIVER called > ',I3,' times for CID=',A)
642         format('even though NFIT_ITERATION = ', I2 )
          CALL  MADABORT(FNAM, C1ERR, C2ERR )
       ENDIF

!  store fit results
       CALL MNFIT_STOREPAR(iter,IERR)

! analyze results: returns ERRFLAG_FIT
!  =  0 => OK; keep SN
!  >  0 => discard fit
!  = -1 => decrement ITER and try again (e.g., for photoZ spike)

        CALL FITPAR_ANA ( isn_proc, iter, ERRFLAG_FIT )

! check ERRFLAG_FIT
        CALL SNANA_ERRFLAG_FIT(iter,ERRFLAG_FIT)
        if ( ERRFLAG_FIT > 0 ) REJECT_FIT = .TRUE.

410     CONTINUE  ! end of iteration loop

      if ( STDOUT_UPDATE ) THEN
        CALL CPU_TIME(t_end)
	  write(6,450) SNLC_CCID(1:ISNLC_LENCCID), t_end-t_start
450       format(T5,'Finished fitting CID = ', A12,  & 
                 ' in ', F6.3,' seconds')
        CALL PRINT_CPU_REMAIN(ISN_ALL)     ! Apr 4 2024
        call flush(6)
      endif
#endif


#if defined(PSNID)
    CALL PSNIDANA(IERR)
    if ( STDOUT_UPDATE ) THEN
        CALL PRINT_CPU_REMAIN(ISN_ALL)     ! Nov 2025
    endif
#endif


42    format(T8,'Analyze ISN_PROC=',I8, 2x,  & 
          'for ',A14,'  NMJD(stored,cuts)=', I4, '->', I4 )
    CALL FLUSH(6)


#if defined(SNANA)
! pack the meta data
    IF ( OPT_TABLE(ITABLE_SNLCPAK) > 0  ) THEN
      CALL SNLCPLOT()       ! prepare light curves for plotting
    ENDIF
#endif

    IF ( OPT_TABLE(ITABLE_SPECPAK) > 0  ) THEN
      CALL SPECPLOT()       ! prepare spectra for plotting (Apr 2019)
    ENDIF


! Beware that SNANA table gets filled after FITRES table.
    CALL MON_SNANA(IFLAG_ANA)     ! monitor-driver

    CALL MAKE_SIMLIB_FILE(2)   ! update SIMLIB Feb 2016


    RETURN
  END SUBROUTINE SNANA_DRIVER

! =====================================================
#if defined(SNFIT)
    SUBROUTINE SNANA_ERRFLAG_FIT(ITER,ERRFLAG)

! Created May 24 2018
! Minor refactor: move some code out of SNANA_DRIVER for cleanup.


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER ITER, ERRFLAG ! (I)

    LOGICAL LAST_ITER, ALREADY_MINOS

! ------------ BEGIN ----------

    LAST_ITER = ( ITER .EQ. NFIT_ITERATION )

    if ( LAST_ITER ) then
       if ( ERRFLAG == 0 ) then
          N_SNLC_FITCUTS = N_SNLC_FITCUTS + 1
          CUTFLAG_SNANA  = IBSET(CUTFLAG_SNANA,1) ! for SNANA ntuple
          NPASSCUT_FIT(ISNLC_TYPE) = NPASSCUT_FIT(ISNLC_TYPE)+1
          NPASSCUT_FIT(-1)         = NPASSCUT_FIT(-1)+1 ! all types

          N_MASK_zSOURCE_LCFIT_CUTS(ISNLC_zSOURCE) =  & 
            N_MASK_zSOURCE_LCFIT_CUTS(ISNLC_zSOURCE) + 1  ! Dec 2024

       endif

! check option to redo entire fit with MINUIT's MINOS option;
! this is a better fit, particularly for errors, but it is
! much slower. So use this only when a serious fit problem is found.
       ALREADY_MINOS = LREPEAT_MINOS .or. USE_MINOS
       if ( ERRFLAG .EQ. -2 .and. (.NOT. ALREADY_MINOS) ) then
          ITER = 0
          LREPEAT_MINOS = .TRUE.
       endif

    endif

! check option to redo fit-iteration
    if ( ERRFLAG .EQ. -1 ) then
       ITER         = ITER - 1
       LREPEAT_ITER = .TRUE.
    else
       LREPEAT_ITER = .FALSE.
    endif

    RETURN
  END SUBROUTINE SNANA_ERRFLAG_FIT
#endif

! ======================================
    SUBROUTINE SNANA_END

! Created MAy 2012 by R.Kessler
! 
! Shell to call functions after all SN have been processed.
! 
! Nov 10 2018: move call to PRINT_SNSTATS() to end of SNANA_END
! 
    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER IERR, LENF
    EXTERNAL write_epoch_list_summary

! ----------------- BEGIN -------------

! -----------------------------

#if defined(SNFIT)
    CALL PRBANNER ( "CALL FITPAR_END" )
    CALL FITPAR_END(IERR)
#endif

#if defined(PSNID)
    CALL PRBANNER ( "CALL PSNIDEND" )
    CALL PSNIDEND(IERR)
#endif


    IF ( USE_TABLEFILE_ROOT ) THEN
       LENF = index(ROOTFILE_OUT,' ') - 1
       CALL TABLEFILE_CLOSE(ROOTFILE_OUT(1:LENF)//char(0), LENF)
    ENDIF

    IF ( USE_TABLEFILE_TEXT ) THEN  ! April 2019

! Apr 18 2022: if no events were processed then there are no text files.
!      In this case, create blank text-table file for each table
       IF ( N_SNLC_PROC == 0 ) CALL SNTABLE_CREATE_TEXT_EMPTY()
       CALL TABLEFILE_CLOSE("NULL"//char(0),10)
    ENDIF

    IF ( USE_TABLEFILE_MARZ ) THEN  ! April 2019
       LENF = index(MARZFILE_OUT,' ') - 1
       CALL TABLEFILE_CLOSE(MARZFILE_OUT(1:LENF)//char(0), LENF)
    ENDIF

    CALL MAKE_SIMLIB_FILE(3)  ! end optional SIMLIB FILE

    CALL WARN_OLDINPUTS("list"//char(0), 0);

    CALL END_MAGCOR()  ! check for warnings

! - - - - - -
    CALL PRINT_SNSTATS()

    IF ( OUT_EPOCH_IGNORE_FILE .NE. '' ) then
       CALL write_epoch_list_summary()
    ENDIF

    IF ( REFORMAT_SNANA )  THEN
       CALL END_REFORMAT()
    ENDIF

    call flush(6)
    print*,' '
    print*,'   ENDING PROGRAM GRACEFULLY. '

    RETURN
  END SUBROUTINE SNANA_END


! ======================================
    SUBROUTINE SNTABLE_CREATE_TEXT_EMPTY()
! 
! Created Apr 18 2022 by R.Kessler
! Goofy routine called at end of job if no events have been read;
! create empty TEXT table so that submit_batch_jobs finds all of
! the expected TEXT-table files. The underlying problem is that
! TEXT table file is not created until after 1st event is read ...
! so if there is no 1st event then there is no text table file.
! For ROOT, the file is created regardless of 1st event
! because these table formats hold multiple tables in one file.
! --------

    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER, PARAMETER :: NTABLE_CHECK = 2

    INTEGER  & 
         ITABLE_LIST(NTABLE_CHECK)  & 
        ,IDTABLE_LIST(NTABLE_CHECK)  & 
        ,itab, ITABLE, IDTABLE, LENNAM, LENFMT
    CHARACTER  & 
          TBNAME_LIST(NTABLE_CHECK)*40  & 
         ,TBNAME*40, TXTFMT*20  & 
         ,cTBNAME*40, cTXTFMT*20

    EXTERNAL SNTABLE_CREATE_TEXT

! --------- BEGIN ----------

    write(6,20)
20    format(/,T3,'No events read -> create empty text-table file:')
    call flush(6)

! hard-wire list of table info to check

    ITABLE_LIST(1)  = ITABLE_SNANA
    ITABLE_LIST(2)  = ITABLE_FITRES

    IDTABLE_LIST(1) = IDTABLE_SNANA
    IDTABLE_LIST(2) = IDTABLE_FITRES

    TBNAME_LIST(1)  = 'SNANA'
    TBNAME_LIST(2)  = 'FITRES'

! - - - - - -
    DO 100 itab =1, NTABLE_CHECK
      ITABLE  = ITABLE_LIST(itab)
      IDTABLE = IDTABLE_LIST(itab)
      TBNAME  = TBNAME_LIST(itab)
      TXTFMT  = 'key'

      IF ( OPT_TABLE(ITABLE) > 0 ) THEN
        LENNAM  = INDEX(TBNAME,' ') - 1
        LENFMT  = INDEX(TXTFMT,' ') -1
        cTBNAME = TBNAME(1:LENNAM) // char(0)
        cTXTFMT = TXTFMT(1:LENFMT) // char(0)
        CALL SNTABLE_CREATE_TEXT(IDTABLE, cTBNAME, cTXTFMT, 40,20)
      ENDIF
100   CONTINUE

    RETURN
  END SUBROUTINE SNTABLE_CREATE_TEXT_EMPTY


! ===========================================
    LOGICAL FUNCTION  REJECT_PRESCALE(N,PRESCALE)

! Created Jan 2015
! Apply &SNLCINP SIM_PRESCALE (if >1) to randomly pre-scale
! the simulation. Use reproducible random numbers based
! on input integer N (incremental index) ,
!   RAN = mod(N,PI)/PI
! 
! Apr 26 2017: pass PRESCALE as argument.


    USE SNPAR

    IMPLICIT NONE

    INTEGER N           ! (I) use for pseudo-random pre-scale
    REAL*8  PRESCALE    ! (I) prescale

    DOUBLE PRECISION XN, RAN, RANMAX

! ------------- BEGIN --------------

    REJECT_PRESCALE = .FALSE.

    IF ( PRESCALE < 1.000001 ) RETURN

    XN     = DBLE(N)
    RAN    = mod(XN,PI)/PI
    RANMAX = 1.0/PRESCALE
    IF ( RAN > RANMAX ) REJECT_PRESCALE = .TRUE.

    RETURN
  END FUNCTION  REJECT_PRESCALE

! =============================
    SUBROUTINE INIT_REFORMAT(OPT)
! 
! Create Feb 7 2021;
! Initialize option to reformat data
! Input OPT=1 -> Auto-set bits
! Input OPT=2 -> call INIT_REFORMAT_[FITS,TEXT,SPECTRA]
! 
! Note that REFORMAT_FITS[TEXT] operates on all of the data: header, photometry,spectra;
! however, REFORMAT_SPECTRA operates only on the spectra and ignores header & photometry
! 
! OPT_REFORMAT_FITS[TEXT] bits:
!  lsb   mask  comment
!   0     1    write fits;
!   1     2    include rejected epochs
!   2     4     SPARE
!   3     8     SPARE
!   4    16     SPARE
!   5    32    rename band -> [SURVEY]-[band]
!   6    64    exclude private variables
!   7   128    include spectra
!   8   256    exclude SIM-truth so that output looks like real data
!   9   512     SPARE
! 
! 
! OPT_REFORMAT_SPECTRA bits:
!  lsb  mask
!   0    1      3-column text format: lam, Flam, Flam_err
!   1    2      3-column text format: lam, Flux, Flux_err
!   2    4      SNANA keyed format for plot_table ; include SIM_FLAM for sims
!   3    8      format for SNID-SAVE (Oct 2025)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!       HISTORY
! July 2023: check OPT_REFORMAT_SPECTRA
! Dec  2024: new OPT_REFORMAT_SPECTRA = 4 for SNANA keyed format
! Apr  2025: new OPT += 32 to rename band -> [SURVEY]-band
!            Minor refactor using IAND and MASKs instead of BTEST
! Oct  2025: add new option for SNID-SAGE format
! 
! -------------------------
    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT

    INTEGER MASK_REFMT_INCLUDE_REJECT  /  2 /   ! include rejected epochs
    INTEGER MASK_REFMT_BAND_NAME       / 32 /   ! band -> [SURVEY]-[band]
    INTEGER MASK_REFMT_EXCLUDE_PRIVATE / 64 /   ! exclude private variables
    INTEGER MASK_REFMT_INCLUDE_SPECTRA / 128 /  ! include spectra
    INTEGER MASK_REFMT_EXCLUDE_SIM     / 256 /  ! exclude SIM truth to look like real data

    INTEGER LEN_VER, MASK_REFMT
    LOGICAL LEXIST, VALIDFILENAME, LMASK
    CHARACTER  & 
         FNAM*14  & 
        ,LIST_FILE*(MXCHAR_FILENAME), IGNORE_FILE*(MXCHAR_FILENAME)  & 
        ,CMD*(MXCHAR_FILENAME)
! ------------- BEGIN -------------

    FNAM = 'INIT_REFORMAT'

    IF ( OPT == 2 ) GOTO 500

    REFORMAT       = .FALSE.  ! any reformat option
! - - - - - - - - -
! convert &SNLCINP strings and OPT_REFORMAT_TEXT[FITS] into logicals

    IF ( VALIDFILENAME(VERSION_REFORMAT_FITS) ) THEN
      OPT_REFORMAT_FITS = IBSET(OPT_REFORMAT_FITS,0) ! set 1st bit
    ENDIF
    IF ( VALIDFILENAME(VERSION_REFORMAT_TEXT) ) THEN
      OPT_REFORMAT_TEXT = IBSET(OPT_REFORMAT_TEXT,0) ! set 1st bit
    ENDIF

! check option to save bad epochs for OPT_REFORMAT
    MASK_REFMT = MASK_REFMT_INCLUDE_REJECT
    LMASK = IAND(OPT_REFORMAT_FITS,MASK_REFMT) > 0 .or.  & 
              IAND(OPT_REFORMAT_TEXT,MASK_REFMT) > 0
    REFORMAT_SAVE_BADEPOCHS = LMASK

! check option to alter band names
    MASK_REFMT = MASK_REFMT_BAND_NAME
    LMASK = IAND(OPT_REFORMAT_FITS,MASK_REFMT) > 0 .or.  & 
              IAND(OPT_REFORMAT_TEXT,MASK_REFMT) > 0
    REFORMAT_BAND_NAME = LMASK

! check option to suppress PRIVATE variables: MASK+=64 or bit 6 (lsb=0)
    MASK_REFMT = MASK_REFMT_EXCLUDE_PRIVATE
    LMASK = IAND(OPT_REFORMAT_FITS,MASK_REFMT) > 0 .or.  & 
              IAND(OPT_REFORMAT_TEXT,MASK_REFMT) > 0
    REFORMAT_PRIVATE = .NOT. LMASK

! check option to reformat spectra: note that for TEXT -> FITS,
! spectra may not appear in 1st text file, so user must explicitly
! request reformatting spectra. MASK=128 -> bit 7

    MASK_REFMT = MASK_REFMT_INCLUDE_SPECTRA
    LMASK = IAND(OPT_REFORMAT_FITS,MASK_REFMT) > 0 .or.  & 
              IAND(OPT_REFORMAT_TEXT,MASK_REFMT) > 0
    REFORMAT_SPECTRA_INCLUDE = LMASK

! check option to exclude SIM truth; MASK += 256
    MASK_REFMT = MASK_REFMT_EXCLUDE_SIM
    LMASK = IAND(OPT_REFORMAT_FITS,MASK_REFMT) > 0 .or.  & 
              IAND(OPT_REFORMAT_TEXT,MASK_REFMT) > 0
    REFORMAT_SIMTRUTH = .NOT. LMASK

! - - - - - - - - - - - - - - - - - -
! check for reformat option using SNANA format

    REFORMAT_SNANA = (OPT_REFORMAT_FITS>0 .or. OPT_REFORMAT_TEXT>0)

    IF ( REFORMAT_SNANA ) THEN

       IF ( OPT_REFORMAT_FITS > 0 ) THEN
          REFORMAT_VERSION = VERSION_REFORMAT_FITS
       ELSE IF ( OPT_REFORMAT_TEXT > 0 ) THEN
          REFORMAT_VERSION = VERSION_REFORMAT_TEXT
       ENDIF

! check if reformat dir already exists
       LEN_VER = INDEX(REFORMAT_VERSION,' ')-1
       INQUIRE(FILE = REFORMAT_VERSION, EXIST = LEXIST)
       IF ( LEXIST ) THEN
          c1err = 'REFORMAT DIR already exists: '  & 
                 // REFORMAT_VERSION(1:LEN_VER)
          c2err = 'Remove it, or change name'
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       ENDIF

       ! create reformat directory
       CMD = 'mkdir ' // REFORMAT_VERSION(1:LEN_VER)
       CALL SYSTEM(CMD)

! open aux files and leave them open
       CALL OPEN_REFORMAT_FILE(LUNLIST2,   'LIST',   LIST_FILE)

       CALL OPEN_REFORMAT_FILE(LUNIGNORE2, 'IGNORE', IGNORE_FILE)
       write(LUNIGNORE2,420)     ! write header
 420     format(T10,'CID      MJD     FILTER')
       NEPOCH_IGNORE_WRFITS = 0
    ENDIF

    REFORMAT = OPT_REFORMAT_FITS>0 .or. OPT_REFORMAT_TEXT>0  & 
            .or. OPT_REFORMAT_SALT2>0

    REFORMAT_SPECTRA_ONLY = (OPT_REFORMAT_SPECTRA > 0 )

! - - - - -
! print info

    IF ( REFORMAT ) THEN
       print*,' '
       print*,'   INIT_REFORMAT options: '
       print*,'     OPT_REFORMAT_[FITS,TEXT] = ',  & 
                    OPT_REFORMAT_FITS, OPT_REFORMAT_TEXT
       print*,'     REFORMAT_BADEPOCHS       = ',  & 
                       REFORMAT_SAVE_BADEPOCHS
       print*,'     REFORMAT_BAND_NAME       = ',  & 
                       REFORMAT_BAND_NAME
       print*,'     REFORMAT_SPECTRA_INCLUDE = ',  & 
                        REFORMAT_SPECTRA_INCLUDE
       print*,'     REFORMAT_PRIVATE         = ',  & 
                        REFORMAT_PRIVATE
       print*,'     REFORMAT_SIMTRUTH        = ',  & 
                        REFORMAT_SIMTRUTH
    ENDIF

    IF ( REFORMAT_SPECTRA_ONLY ) THEN
       print*,' '
       print*,'   INIT_REFORMAT options: '
       print*,'     OPT_REFORMAT_SPECTRA  = ', OPT_REFORMAT_SPECTRA
    ENDIF

       call flush(6)

    RETURN
! - - - - - - - - - -

500   CONTINUE

    IF ( OPT_REFORMAT_FITS  >  0 ) THEN
       CALL INIT_REFORMAT_FITS()
    ENDIF

    IF ( OPT_REFORMAT_TEXT  > 0 ) THEN
       CALL INIT_REFORMAT_TEXT()
    ENDIF

    IF ( OPT_REFORMAT_SALT2 .EQ. 2 ) THEN
       CALL INIT_REFORMAT_SALT2()
    ENDIF

    IF ( OPT_REFORMAT_SPECTRA > 0 ) THEN
       CALL INIT_REFORMAT_SPECTRA_ONLY()
    ENDIF

    RETURN
  END SUBROUTINE INIT_REFORMAT

! ============================================
    SUBROUTINE INIT_REFORMAT_FITS

! Prepare for translating text file format into FITS format.
! 
    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    CHARACTER  & 
         cPATH*(MXCHAR_PATH)  & 
        ,cVERSION*(MXCHAR_VERSION)  & 
        ,cPREFIX*(MXCHAR_VERSION)  & 
        ,cHEADFILE*(MXCHAR_FILENAME)  & 
        ,LIST_FILE*(MXCHAR_FILENAME)  & 
        ,IGNORE_FILE*(MXCHAR_FILENAME)

    INTEGER  & 
         WRITE_FLAG, LL, NMARK  & 
        ,LEN_DIR  & 
        ,LEN_SURVEY  & 
        ,LEN_VERS  & 
        ,LEN_HEAD

! ------------------ BEGIN ----------------

    LEN_SURVEY = INDEX(SURVEY_NAME,' ') - 1

    LEN_VERS = INDEX(VERSION_REFORMAT_FITS,' ') - 1
    cVERSION = VERSION_REFORMAT_FITS(1:LEN_VERS) // char(0)
    cPREFIX  = cVERSION

    LEN_DIR    = INDEX(REFORMAT_VERSION,' ') - 1
    cPATH      = REFORMAT_VERSION(1:LEN_DIR) // char(0)

! to see flags in C code, "grep WRITE_MASK sndata.h"
    IF ( LSIM_MAGOBS ) THEN
       WRITE_FLAG = 8   ! Jan 23 2018
    ELSE IF ( LSIM_SNANA ) THEN
       WRITE_FLAG = 4   ! Jan 23 2018
       if ( .NOT. REFORMAT_SIMTRUTH ) WRITE_FLAG = 2 ! Mar 2022
    ELSE
       WRITE_FLAG = 2
    ENDIF

    if ( REFORMAT_SPECTRA_INCLUDE ) WRITE_FLAG = WRITE_FLAG + 128
    if ( FOUND_ATMOS              ) WRITE_FLAG = WRITE_FLAG + 256

! -----------------------
    CALL WR_SNFITSIO_INIT(cPATH, cVERSION, cPREFIX, WRITE_FLAG  & 
                 , NMARK, cHEADFILE   &  ! output name
                 , LEN_DIR, LEN_VERS, LEN_VERS, LEN_HEAD)

! --- write name of HEAD file to already-opend LIST file

    LL = INDEX(cHEADFILE, char(0) ) - 1
    WRITE(LUNLIST2,50) cHEADFILE(1:LL)
50    FORMAT(A)
    CLOSE(LUNLIST2)

    RETURN
  END SUBROUTINE INIT_REFORMAT_FITS

! ===============================
    SUBROUTINE INIT_REFORMAT_TEXT()
!       print*,'  INI_WRTEXT: placeholder -> do nothing ... yet'
    RETURN
  END SUBROUTINE INIT_REFORMAT_TEXT

! ====================================
    SUBROUTINE INIT_REFORMAT_SALT2
! -----------------------------------------------
! Nov 2011 R.Kessler
! 
! Initialization for translating SNANA formatted data
! into SALT2-formatted data. The input options are fed
! vis &SNLCINP namelist string
! 
!  &SNLCINP
!     ...
!   REFORMAT_KEYS =
!   '@INSTRUMENT <instr> @MAGSYS <magsys> @PREFIX <prefix> @REPLACE <f1> <f2>
!     ...
!  &END
! 
! where @INSTRUMENT and @MAGSYS are required and the others are optional.
! Definition of above keys
! 
! @INSTRUMENT = name of telescope or survey defined by SALT2
! @MAGSYS     = mag system defined by SALT2 (i.e, AB, VEGA ...)
! @PREFIX     = file-name prefix (default prefix is name of survey)
! @REPLACE    = replace filter list <f1> with <f2>. For example,
!               if f1 = UGRIZ and f2 = ugriz, then U -> u, G -> g,
!               etc in the output SALT2 files.
! 
! 
! ------------------------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE WRS2COM

    IMPLICIT NONE

    INTEGER iwd, NWD, ifilt, NTMP, IFILT_OBS, IFILT_TMP(2)
    INTEGER LL, LL0, LL1, LL2, MSKOPT
    character  & 
         cwd*(MXCHAR_FILEWORD-2)  & 
        ,cwd1*(MXCHAR_FILEWORD-2)  & 
        ,cwd2*(MXCHAR_FILEWORD-2)  & 
        ,ctmp*(MXCHAR_FILEWORD-2)  & 
        ,cfilt(2)*2, cfilt1*1  & 
        ,cutvar_file*(MXCHAR_FILENAME)  & 
        ,list_file*(MXCHAR_FILENAME)  & 
        ,FNAM*30

! function
    INTEGER FILTINDX

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS
! -------------------- BEGIN -------------------

    FNAM = 'INIT_REFORMAT_SALT2'
    CALL PRBANNER("INI_WRSALT2: Translate SNANA -> SALT2 format")

    LL = INDEX(REFORMAT_KEYS,' ') -1
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    NWD = STORE_PARSE_WORDS(MSKOPT, REFORMAT_KEYS(1:LL)//char(0),  & 
               FNAM//char(0), LL, 30)

    NAMEof_SURVEY        = 'NULL'
    NAMEof_INSTRUMENT    = 'NULL'
    NAMEof_MAGSYS        = 'NULL'
    NAMEof_PREFIX        =  SURVEY_NAME
    NAMEof_REPLACE(1)    = 'NULL'
    NAMEof_REPLACE(2)    = 'NULL'

    NREPLACE = 0
    NEWKEY = 0

    DO ifilt = 1, MXFILT_ALL
         IMAP_REPLACE(IFILT) = -9
    ENDDO

    DO iwd = 1, NWD

      CALL get_PARSE_WORD_fortran(iwd+0, cwd,  LL0 )
      CALL get_PARSE_WORD_fortran(iwd+1, cwd1, LL1 )

      if ( cwd(1:LL0) .EQ. '@SURVEY' ) then
        NAMEof_SURVEY = cwd1

      else if ( cwd(1:LL0) .EQ. '@INSTRUMENT' ) then
        NAMEof_INSTRUMENT = cwd1

      else if ( cwd(1:LL0) .EQ. '@MAGSYS' ) then
        NAMEof_MAGSYS = cwd1

      else if ( cwd(1:LL0) .EQ. '@PREFIX' ) then
        NAMEof_PREFIX = cwd1

      else if ( cwd(1:LL0) .EQ. '@REPLACE' ) then
         CALL get_PARSE_WORD_fortran(iwd+2, cwd2, LL2 )
        NAMEof_REPLACE(1) = cwd1
        NAMEof_REPLACE(2) = cwd2
        NREPLACE    = INDEX(NAMEof_REPLACE(1), ' ' ) - 1
      else if ( cwd(1:1) .EQ. '@' ) then
        NEWKEY = NEWKEY + 1
        NEWKEY_NAME(NEWKEY)(1:LL0-1) = cwd(2:LL0)
        NEWKEY_ARG(NEWKEY)           = cwd1
      endif

    ENDDO

    LEN_SURVEY  = INDEX(NAMEof_SURVEY,    ' ' ) - 1
    LEN_INST    = INDEX(NAMEof_INSTRUMENT,' ' ) - 1
    LEN_MAGSYS  = INDEX(NAMEof_MAGSYS,    ' ' ) - 1
    LEN_PREFIX  = INDEX(NAMEof_PREFIX,    ' ' ) - 1

! make sure that required keys are specified.

    if ( NAMEof_INSTRUMENT .EQ. 'NULL' ) then
      c1err = 'MUST specify @INSTRUMENT <instrument> in '
      c2err = '&SNLCINP namelist string REFORMAT_KEYS'
      CALL  MADABORT("WRSALT2_2", C1ERR, C2ERR)
    endif
    if ( NAMEof_MAGSYS .EQ. 'NULL' ) then
      c1err = 'MUST specify @MAGSYS <magsys> in '
      c2err = '&SNLCINP namelist string REFORMAT_KEYS'
      CALL  MADABORT("WRSALT2_2", C1ERR, C2ERR)
    endif

! ----------
! Check for filter-name substitutions;
! i.e., UGRIZ -> ugriz for the SDSS

    IF ( NREPLACE .GT. 0  ) THEN

! abort if filter-strings have different length
      NTMP   = INDEX(NAMEof_REPLACE(2), ' ' ) - 1
      IF ( NREPLACE .NE. NTMP ) THEN
        c1err = 'Filter @REPALCE strings have different length'
        c2err = 'Cannot replace ' //  & 
                    NAMEof_REPLACE(1)(1:NREPLACE) // ' with ' //  & 
                    NAMEof_REPLACE(2)(1:NTMP)
        CALL  MADABORT("WRSALT2_2", C1ERR, C2ERR)
      ENDIF

! create map between original filter and subst. filter

      DO 55 ifilt = 1, NREPLACE
         cfilt(1) = NAMEof_REPLACE(1)(ifilt:ifilt)
         cfilt(2) = NAMEof_REPLACE(2)(ifilt:ifilt)
         IFILT_TMP(1) = FILTINDX ( cfilt(1) )
         IFILT_TMP(2) = FILTINDX ( cfilt(2) )
         IMAP_REPLACE(IFILT_TMP(1)) = IFILT_TMP(2)

         write(6,56) cfilt, IFILT_TMP
56         format(t10,'Prepare @REPLACE map for ' , A2,' -> ', A2,  & 
             '(', I2,' -> ', I2, ')' )

55      CONTINUE
    ENDIF

! ----------------------------
!  open list-file and leave it open.

    LIST_FILE = NAMEof_PREFIX(1:LEN_PREFIX) // '.LIST'
    print*,'   Open list-file: ', LIST_FILE(1:LEN_PREFIX+12)

    OPEN( UNIT   = LUNSALT2  & 
          , FILE   = LIST_FILE  & 
          , STATUS = 'UNKNOWN'  & 
                 )

! ----------------------------
! write cut-def file, then close it.

    CUTVAR_FILE = NAMEof_PREFIX(1:LEN_PREFIX) // '_CUTVAR.LOG'
    print*,'   Write cut-definitions to ' //  & 
          CUTVAR_FILE(1:LEN_PREFIX+12)

    OPEN( UNIT   = LUNTMP  & 
          , FILE   = CUTVAR_FILE  & 
          , STATUS = 'UNKNOWN'  & 
                 )

    write(LUNTMP,111) 'Z_HELIO   ',  & 
             'heliocentric redshift'
    write(LUNTMP,111) 'MWEBV     ',  & 
             'Galactic E(B-V)'
    write(LUNTMP,111) 'MWEBV_ERR     ',  & 
             'error on Galactic E(B-V)'

      write(ctmp,401) int(cutwin_trest(1)), int(cutwin_trest(2))
401     format(I3,'<Trest<' , I3,' days')
      write(LUNTMP,111) 'NOBS      ',  & 
             'Nobs total with any S/N and ' // ctmp(1:20)
      write(LUNTMP,111) 'TRESTMIN  ',  & 
             'min Trest(days) relative to peak (any S/N)'
      write(LUNTMP,111) 'TRESTMAX  ',  & 
             'max Trest(days) relative to peak (any S/N)'

      write(LUNTMP,111) 'T0GAPMAX  ',  & 
             'max rest-frame gap (days) that overlaps peak epoch'

      write(LUNTMP,111) 'SNRMAX    ',  & 
             'max S/N among all observations'
      write(LUNTMP,111) 'SNRMAX2   ',  & 
             'max S/N excluding filter with SNRMAX'

      write(LUNTMP,111) 'SNRMAX3   ',  & 
             'max S/N excluding filters with SNRMAX   SNRMAX2'


      write(LUNTMP,600) ' '
      DO IFILT     = 1, NFILTDEF_SURVEY
         ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
         cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
         ctmp      = 'SNRMAX_' // cfilt1
         LL        = INDEX(ctmp,' ') - 1
         write(LUNTMP,111) ctmp(1:10),  & 
             'max S/N for indicated filter'
      ENDDO

    IF ( OPT_SETPKMJD > 0 ) THEN
      write(LUNTMP,600) ' '
      DO IFILT     = 1, NFILTDEF_SURVEY
         ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
         cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
         ctmp      = 'ERRT0_' // cfilt1
         LL        = INDEX(ctmp,' ') - 1
         write(LUNTMP,111) ctmp(1:10),  & 
             'fitted T0 error (days) for indicated filter'
      ENDDO

      write(LUNTMP,600) ' '
      DO IFILT     = 1, NFILTDEF_SURVEY
         ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
         cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
         ctmp      = 'ERRF0_' // cfilt1
         LL        = INDEX(ctmp,' ') - 1
         write(LUNTMP,111) ctmp(1:10),  & 
             'ERR(F0)/F0 for indicated filter'
      ENDDO

    ENDIF

111     format('@', A, ' : ', A)
600     format(A)
      CLOSE(UNIT=LUNTMP)

    print*,' -----------------------------------------------'
    print*,' '
    RETURN
  END SUBROUTINE INIT_REFORMAT_SALT2

! ====================================
    SUBROUTINE INIT_REFORMAT_SPECTRA_ONLY()

! Open SPEC.LIST summary file and write table header
! Beware that using OPT_REFORMAT_SPECTRA option on multiple data folders
! results in clobbering SPEC.LIST file.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

    INTEGER NVAR, ivar, LL, LENV
    CHARACTER VARNAMES_LIST*100, FNAM*28, VARNAME*20, COMMENT*80
    CHARACTER LINE*60

! --------------- BEGIN ------------

    FNAM = 'INIT_REFORMAT_SPECTRA_ONLY'
    CALL PRBANNER(FNAM)

    LINE = '========================================'  & 
            // '================== '
    WRSPEC_FLAM = BTEST(OPT_REFORMAT_SPECTRA,0)
    WRSPEC_FLUX = BTEST(OPT_REFORMAT_SPECTRA,1)
    WRSPEC_KEY  = BTEST(OPT_REFORMAT_SPECTRA,2)
    WRSPEC_SNID_SAGE = BTEST(OPT_REFORMAT_SPECTRA,3)

    write(6,120)  & 
         WRSPEC_FLAM, WRSPEC_FLUX, WRSPEC_KEY, WRSPEC_SNID_SAGE
120   format(T5, 'WRSPEC[FLAM,FLUX,KEY,SNID_SAGE] = ', 4L3)


    IF ( WRSPEC_SNID_SAGE ) THEN  ! set FLAM bit in case user doesn't
       OPT_REFORMAT_SPECTRA = IBSET(OPT_REFORMAT_SPECTRA,0)
	 WRSPEC_FLAM = .TRUE.
    ENDIF

! define all possible columns for spec table

    NVAR = 0

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'CID'
    COMMENTS_SPEC_TABLE(NVAR) = 'Candidate ID (or SNID) of event'

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'SPEC_FILE'
    COMMENTS_SPEC_TABLE(NVAR) =  & 
           'name of table file with WAVE FLAM FLAMERR'

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'TYPE'
    COMMENTS_SPEC_TABLE(NVAR) = 'string type, e.g. Ia, Ib, II ... '

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'SUBTYPE'
    COMMENTS_SPEC_TABLE(NVAR) = 'string sub type '
    IVAR_SPEC_SUBTYPE = NVAR

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'REDSHIFT'
    COMMENTS_SPEC_TABLE(NVAR) = 'redshift (-9 for unknown) '

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'MJD'
    COMMENTS_SPEC_TABLE(NVAR) = 'MJD of spectrum'
    IVAR_SPEC_MJD = NVAR

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'TOBS'
    COMMENTS_SPEC_TABLE(NVAR) = 'MJD - PEAKMJD (days) '
    IVAR_SPEC_TOBS = NVAR

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'TREST'
    COMMENTS_SPEC_TABLE(NVAR) = 'TOBS/(1+z)'
    IVAR_SPEC_TREST = NVAR

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'NBIN_WAVE'
    COMMENTS_SPEC_TABLE(NVAR) = 'Number of WAVE bins'
    IVAR_SPEC_NBIN = NVAR

    NVAR = NVAR + 1
    VARNAMES_SPEC_TABLE(NVAR) = 'ISDATA_REAL'
    COMMENTS_SPEC_TABLE(NVAR) = '1 -> real data, 0 -> sim data'

    if ( NVAR > MXVAR_SPEC_TABLE ) then
      write(c1err,61) NVAR, MXVAR_SPEC_TABLE
61      format('NVAR(SPEC_TABLE) = ',I3,' exceeds bound of ', I3)
      c2err = 'Consider increasing MXVAR_SPEC_TABLE in snana.F90'
      CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    NVAR_SPEC_TABLE = NVAR
    DO ivar = 1, NVAR
       USECOL_SPEC_TABLE(ivar) = .TRUE.
    ENDDO

! Now remove columns based on output format option
    IF ( WRSPEC_SNID_SAGE ) THEN
       USECOL_SPEC_TABLE(IVAR_SPEC_MJD)   = .FALSE.
       USECOL_SPEC_TABLE(IVAR_SPEC_TOBS)  = .FALSE.
       USECOL_SPEC_TABLE(IVAR_SPEC_NBIN)  = .FALSE.
    ELSE
       USECOL_SPEC_TABLE(IVAR_SPEC_SUBTYPE) = .FALSE.
    ENDIF

! 

! ---------------------------------------------------------
! open table list file to give one per spectrum summary

    SPEC_TABLE_FILE = 'SPEC.LIST'
    print*,' Open ', SPEC_TABLE_FILE

    OPEN (UNIT=LUNTMP, FILE=SPEC_TABLE_FILE, STATUS='UNKNOWN' )
    WRITE(LUNTMP,700) LINE

! construct list of varnames based on USECOL_SPEC_TABLE list of logicals
! Also write comment to SPEC_TABLE_FILE for each used variable
    VARNAMES_LIST = ''
    LL = 0
    DO ivar = 1, NVAR
       if (  USECOL_SPEC_TABLE(ivar) ) then
           VARNAME = VARNAMES_SPEC_TABLE(ivar)
           LENV    = INDEX(VARNAME,' ')
           VARNAMES_LIST = VARNAMES_LIST(1:LL)//VARNAME(1:LENV) // ' '
           LL = LL + LENV + 1

           COMMENT = COMMENTS_SPEC_TABLE(ivar)
           WRITE(LUNTMP,700) VARNAME(1:13) // ' : ' // COMMENT(1:40)
       endif
    ENDDO

700   FORMAT('# ', A)

    WRITE(LUNTMP,700) LINE
    WRITE(LUNTMP,701) VARNAMES_LIST  ! avoid blank space
701   FORMAT(A)

    RETURN
  END SUBROUTINE INIT_REFORMAT_SPECTRA_ONLY


! ================================================
    SUBROUTINE OPEN_REFORMAT_FILE(LUN, SUFFIX, REFORMAT_FILE)
! Open text file name REFORMAT_VERSION/REFORMAT_VERSION.SUFFIX
! use imput LUN

    USE SNDATCOM

    IMPLICIT NONE

    INTEGER   LUN          ! (I) open with this logical unit
    CHARACTER SUFFIX*(*)   ! (I) open with this suffix
    CHARACTER REFORMAT_FILE*(*)  ! (O) name of file opened

    INTEGER   LEN_VER

! ------------ BEGIN -----------
    LEN_VER = INDEX(REFORMAT_VERSION,' ') - 1
    REFORMAT_FILE =  & 
         REFORMAT_VERSION(1:LEN_VER) // '/' //  & 
         REFORMAT_VERSION(1:LEN_VER) // '.' // SUFFIX

    OPEN (LUN, FILE = REFORMAT_FILE, status='UNKNOWN' )

    RETURN
  END SUBROUTINE OPEN_REFORMAT_FILE

! =============================
    SUBROUTINE EXEC_REFORMAT(IVERS)
! 
! Create Feb 7 2021;
! Update reformat output:
!  + append FITS files for FITS output
!  + write new TEXT file for TEXT output
! 
! Jan 20 2022: move copy_SNDATA_MISC() to work for either TEXT for FITS.
! Jul 18 2023: check OPT_REFORMAT_SPECTRA
! Apr 19 2025: check new option for BAND_NAME
! Oct 22 2025: pass IVERS index for VERSION_PHOTOMETRY
! -------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IVERS  ! (I)

    INTEGER, PARAMETER :: MXWR_TEXT = 400

    LOGICAL    NOT_TOO_MANY
    CHARACTER  CCID*(MXCHAR_CCID), TEXTFILE*100, TEXTFILE_FULL*100
    INTEGER    LEN_VER, LEN_TXT, LEN_S
    EXTERNAL   WR_SNTEXTIO_DATAFILE

! ---------- BEGIN ----------

    if ( .NOT. ( REFORMAT .or. REFORMAT_SPECTRA_ONLY)  ) RETURN

    IF ( REFORMAT_SNANA ) THEN
       CALL WRITE_REFORMAT_IGNORE()  ! update ignore file
    ENDIF

    CALL copy_SNDATA_MISC() ! update SNDATA structure

    IF ( OPT_REFORMAT_FITS > 0 )  THEN
       CALL WR_SNFITSIO_UPDATE()
    ENDIF

    NOT_TOO_MANY = (N_SNLC_CUTS .LE. MXWR_TEXT)
    IF ( OPT_REFORMAT_TEXT > 0  .and. NOT_TOO_MANY ) THEN
!        prepare name of text file
       CCID    = SNLC_CCID(1:ISNLC_LENCCID)
       LEN_VER = INDEX(REFORMAT_VERSION,' ') - 1
       LEN_S   = INDEX(SURVEY_NAME,' ') - 1
       textFile =  & 
!  xxx     &       REFORMAT_VERSION(1:LEN_VER) // '_' //
             SURVEY_NAME(1:LEN_S) // '_' //  & 
             SNLC_CCID(1:ISNLC_LENCCID) // '.DAT'
       LEN_TXT = INDEX(textfile,' ') - 1
       TEXTFILE_FULL = REFORMAT_VERSION(1:LEN_VER) // '/'  //  & 
             textfile(1:len_txt) // char(0)

       CALL WR_SNTEXTIO_DATAFILE(TEXTFILE_FULL, 100)
       write(LUNLIST2, '(A)' ) TEXTFILE(1:len_txt)

       if ( N_SNLC_CUTS == MXWR_TEXT ) then
          print*,' '
          print*,'  WARNING: wrote max number of TEXT files: ',  & 
                      N_SNLC_CUTS
          print*,' '
       endif
    ENDIF

! check Julien's original SALT2 format (not SNANA format)
    IF ( OPT_REFORMAT_SALT2 .EQ. 2 ) THEN
       CALL WRSALT2_2()
    ENDIF

    IF ( REFORMAT_SPECTRA_ONLY ) THEN
      CALL WR_SPECTRA(IVERS)
    ENDIF

    RETURN
  END SUBROUTINE EXEC_REFORMAT

! =====================================
    SUBROUTINE copy_SNDATA_MISC()

! Created Mar 14 2021
! Update the following SNDATA struct variables for FITS or TEXT format:
! 
!   + PEAKMJD                ( if OPT_SETPKMJD > 0)
!   + MWEBV[_ERR]            ( if OPT_MWEBV    > 0)
!   + REDSHIFT_CMB[_ERR]
!   + REDSHIFT_HELIO[_ERR]
!   + VPEC[_ERR]
!   + NPRIVATE               set to 0 OPT_REFORMAT_PRIVATE=F (Jan 2022)
!   + HOSTGAL_LOGMASS[_ERR]  ??? TO DO (how to handle mutliple hosts?)
!   + HOSTGAL_sSFR[_ERR]     ??? TO DO (idem)
! 
!   + FLUXCAL      # if MAGCOR_INFILE is set  (Dec 2021)
!   + FLUXCALERR   # if FLUXMODELERR_FILE is set (Dec 2021)
! 
!   + MJD_TRIGGER        # if PHOTFLAG_TRIGGER > 0
!   + MJD_DETECT_FIRST   # if PHOTFLAG_DETECT > 0
!   + MJD_DETECT_LAST    # idem
! 
!   + MWMAGCOR_[band]  # Nov 2023
!   + BAND name        # Apr 2025
! 
!   + photo-z quantiles  # Oct 14 2025
! -------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM
    USE FILTCOM
    USE FITSCOM

    IMPLICIT NONE

    INTEGER COPYFLAG, NARG, NOBS, LEN_KEY, LEN_STR, o
    INTEGER IFILT, IFILT_OBS, LEN_SURVEY, NQ, q, igal
    REAL*8 DVAL(MXEPOCH)
    CHARACTER cKEY*20, cSTRING*20, cfilt*2
    EXTERNAL COPY_SNDATA_HEAD

! -------------- BEGIN ----------

    NARG     =  1
    NOBS     =  ISNLC_NEWMJD_FOUND
    COPYFLAG = +1
    LEN_KEY  = 20
    LEN_STR  = 20
    cSTRING  = "DUMMY" // char(0)

    IF ( OPT_SETPKMJD > 0 ) THEN
      cKEY     = "PEAKMJD" // char(0)
      DVAL(1)     = DBLE(SNLC_SEARCH_PEAKMJD)
      CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
    ENDIF

    IF ( OPT_MWEBV > 0 ) THEN
      cKEY     = "MWEBV" // char(0)
      DVAL(1)  = DBLE(SNLC_MWEBV)
      CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

      cKEY     = "MWEBV_ERR" // char(0)
      DVAL(1)  = DBLE(SNLC_MWEBV_ERR)
      CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

! write MWXT mag as diagnostic (not used)
      do ifilt     = 1, NFILTDEF_SURVEY
         ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
         cfilt     = filtdef_string(ifilt_obs:ifilt_obs)
         cKEY      = 'MWXT_MAG_' // cfilt(1:1) // char(0)
         DVAL(1)   = SNLC_MWXT_MAG(ifilt)
         CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
                 cSTRING, DVAL, LEN_KEY, LEN_STR)
      end do

    ENDIF

! always update redshift and vpec info in case of header-override
! or other calculated update.

! zCMB
    cKEY     = "REDSHIFT_CMB" // char(0)
    DVAL(1)  = DBLE(SNLC_zCMB)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

    cKEY     = "REDSHIFT_CMB_ERR" // char(0)
    DVAL(1)  = DBLE(SNLC_zCMB_ERR)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

! zHELIO
    cKEY     = "REDSHIFT_HELIO" // char(0)
    DVAL(1)  = DBLE(SNLC_zHELIO)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

    cKEY     = "REDSHIFT_HELIO_ERR" // char(0)
    DVAL(1)  = DBLE(SNLC_zHELIO_ERR)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

! VPEC
    cKEY     = "VPEC" // char(0)
    DVAL(1)  = DBLE(SNLC_VPEC)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

    cKEY     = "VPEC_ERR" // char(0)
    DVAL(1)  = DBLE(SNLC_VPEC_ERR)
    CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)

! Oct 2025: photo-z quantiles
    NQ = SNHOST_NZPHOT_Q
    IF ( NQ > 0 ) THEN
      IGAL = 1
      do q = 1, NQ
         cKEY    = VARNAME_ZPHOT_Q(IGAL,q) // char(0)
         DVAL(1) = SNHOST_ZPHOT_Q(IGAL,q)
         CALL copy_SNDATA_HEAD(COPYFLAG, cKEY, NARG,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
      enddo
    ENDIF

! check option to exclude PRIVATE variables from output.
! Do NOT modify NVAR_PRIVATE.
    IF ( .NOT. REFORMAT_PRIVATE .and. NVAR_PRIVATE > 0 ) then
      cKEY     = "NVAR_PRIVATE"  // char(0)
      DVAL(1)  = -1.0   ! -1 means disable to avoid abort on NVAR change
      CALL copy_SNDATA_GLOBAL(COPYFLAG, cKEY, NARG,  & 
                cSTRING, DVAL, LEN_KEY, LEN_STR)
    ENDIF

! ----
! update flux[err] if modified
! ------
    IF ( DOFUDGE_FLUXERRMODEL .or. FUDGE_MAG_ERROR.NE.'') THEN
       cKEY  = "FLUXCALERR" // char(0)
       DO o=1,NOBS; DVAL(o)=SNLC_FLUXCAL_ERRTOT(o) ; END DO
       CALL copy_SNDATA_OBS(copyFlag, cKEY, NOBS,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
    ENDIF

    IF ( NSTORE_MAGCOR > 0 ) THEN
       cKEY  = "FLUXCAL" // char(0)
       DO o=1,NOBS; DVAL(o)=SNLC_FLUXCAL(o) ; END DO
       CALL copy_SNDATA_OBS(copyFlag, cKEY, NOBS,  & 
              cSTRING, DVAL, LEN_KEY, LEN_STR)
    ENDIF

! Apr 19 2025: check modifying band name
    IF ( REFORMAT_BAND_NAME ) THEN
       cKEY = "BAND" // char(0)
       LEN_SURVEY = INDEX(SURVEY_NAME,' ') - 1
       STRFITS = ''
       LEN_STR = 0
       DO o=1,NOBS
          IFILT_OBS  = ISNLC_IFILT_OBS(o)
          cfilt      = filtdef_string(ifilt_obs:ifilt_obs)
          cSTRING    = SURVEY_NAME(1:LEN_SURVEY) // '-' // cfilt(1:1)
          STRFITS    = STRFITS(1:LEN_STR) // ' ' // cSTRING
          LEN_STR    = LEN_STR + LEN_SURVEY + 3
       ENDDO

       CALL copy_SNDATA_OBS(copyFlag, cKEY, NOBS,  & 
              STRFITS(1:LEN_STR)//char(0), DVAL, LEN_KEY, LEN_STR)
       LEN_STR = 20
    ENDIF


! Sep 2023: update MJD_TRIGGER and MJD_DETECT_XYZ variables
!           if they have been computed from PHOTFLAG .

    if ( PHOTFLAG_TRIGGER > 0 ) then
       cKEY = "MJD_TRIGGER" // char(0)
       DVAL(1) = SNLC8_MJD_TRIGGER
       CALL copy_SNDATA_HEAD(copyFlag, cKEY, NARG,  & 
               cSTRING, DVAL, LEN_KEY, LEN_STR)
    endif

    if ( PHOTFLAG_DETECT > 0 ) then
       cKEY = "MJD_DETECT_FIRST" // char(0)
       DVAL(1) = SNLC8_MJD_DETECT_FIRST
       CALL copy_SNDATA_HEAD(copyFlag, cKEY, NARG,  & 
               cSTRING, DVAL, LEN_KEY, LEN_STR)

       cKEY = "MJD_DETECT_LAST" // char(0)
       DVAL(1) = SNLC8_MJD_DETECT_LAST
       CALL copy_SNDATA_HEAD(copyFlag, cKEY, NARG,  & 
               cSTRING, DVAL, LEN_KEY, LEN_STR)
    endif

! - - - - -

    RETURN
  END SUBROUTINE copy_SNDATA_MISC

! =========================================
    SUBROUTINE WRITE_REFORMAT_IGNORE()

! Open and write IGNORE file


    USE SNDATCOM

    IMPLICIT NONE

    INTEGER i, LEN1_CCID, LEN2_CCID
    CHARACTER CCID*(MXCHAR_CCID)
! ----------- BEGIN ----------

    LEN1_CCID = INDEX(SNLC_CCID,' ') - 1
    DO 400 i = 1, NEPOCH_IGNORE
       CCID = EPOCH_IGNORE_CCID(i)
       LEN2_CCID = INDEX(CCID,' ') - 1
       IF( LEN1_CCID .NE. LEN2_CCID ) GOTO 400
       IF ( CCID(1:LEN1_CCID) .EQ. SNLC_CCID(1:LEN1_CCID) ) THEN
          write(LUNIGNORE2,410)  & 
            CCID(1:LEN1_CCID), EPOCH_IGNORE_MJD(i), EPOCH_IGNORE_FILT(i)
 410        format('IGNORE:  ', A, 3x, F9.3, 3x, A1)
          call flush(LUNIGNORE2)
          NEPOCH_IGNORE_WRFITS = NEPOCH_IGNORE_WRFITS + 1
       ENDIF
 400  CONTINUE

    RETURN
  END SUBROUTINE WRITE_REFORMAT_IGNORE


! ===============================
    SUBROUTINE WR_SPECTRA(IVERS)

! Created July 2023
! For each spectrum, re-write using format based on OPT_REFORMAT_SPECTRA.
! File name is
!   [version]-snid[snid]-mjd[mjd].spec
! 
! Iniitial code writes only one format, but later this could be expanded
! to write different formats.
! 
! Mar 20 2025: update output to include SIM_FLAM for all options
!              (instead of only KEY option)
! 
! Oct 22 2025: pass IVERS index for VERSION_PHOTOMETRY
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

    INTEGER IVERS ! (I)

! local args

    INTEGER LENV, LENC, LENROW, LENF, LENT, LENTT, LENR
    INTEGER ispec, NLAMBIN, ilam, ISDATA_REAL
    REAL*8  MJD, LAM, LAMBIN, F, FERR, FSIM, TOBS, TREST, z
    CHARACTER SPEC_FILE*(MXCHAR_FILENAME), MSG*80, MSGTMP*80
    CHARACTER VERS*(MXCHAR_FILENAME), TYPE*20,SUBTYPE*20

! ------------- BEGIN -----------

! open table list file to give summary info for each spectrum
    OPEN (UNIT=LUNTMP, FILE=SPEC_TABLE_FILE, STATUS='UNKNOWN',  & 
             ACCESS='APPEND' )

    VERS = VERSION_PHOTOMETRY(IVERS)
    LENV = INDEX(VERS,' ') - 1
    LENC = ISNLC_LENCCID

    ISDATA_REAL = 1
    if ( LSIM_SNANA ) ISDATA_REAL = 0

    z       = SNLC_REDSHIFT

    DO 400 ispec = 1, NSPECTRUM

       CALL RDSPEC_DRIVER(ispec)  ! load LAM and FLAM arrays

       MJD     = MJD_SPECTRUM(ispec)
	 TOBS    = TOBS_SPECTRUM(ispec)
	 TREST   = TOBS/(1.0+z)
       NLAMBIN = NLAMBIN_SPECTRUM(ispec)

       if ( MJD > 0.0 ) THEN
           WRITE(SPEC_FILE,40) VERS(1:LENV), SNLC_CCID(1:LENC), MJD
40          format(A,'_',A,'_', F9.3, '.SPEC' )
       ELSE
           WRITE(SPEC_FILE,41)  & 
             VERS(1:LENV), SNLC_CCID(1:LENC), "HOST"
41          format(A,'_',A,'_', A, '.SPEC' )
       ENDIF
       LENF = INDEX(SPEC_FILE,' ') - 1


       TYPE    = 'NULL'
       SUBTYPE = 'NULL'
       IF ( LSIM_SNANA ) THEN
          CALL GET_SIMTYPE_SPEC_FILE(TYPE,SUBTYPE)
       ENDIF
       LENT  = INDEX(TYPE,' ') - 1  !
       LENTT = INDEX(SUBTYPE,' ') - 1

       if ( WRSPEC_SNID_SAGE ) then
         WRITE(LUNTMP,709)  & 
              SNLC_CCID(1:LENC), SPEC_FILE(1:LENF),  & 
              TYPE(1:LENT), SUBTYPE(1:LENTT), z,  & 
              TREST, ISDATA_REAL
709        FORMAT(A, 2x, A, 2x,  & 
                  A, 2x, A, 2x, E9.3, 2x,  & 
                  F8.1, I3)
       else
!             write few more columns for generic SNANA SPEC table
         WRITE(LUNTMP,710)  & 
              SNLC_CCID(1:LENC), SPEC_FILE(1:LENF),  & 
              TYPE(1:LENT), z,  & 
              MJD, TOBS, TREST,  & 
              NLAMBIN, ISDATA_REAL
710        FORMAT(A, 2x, A, 2x,  & 
                  A, 2x, E9.3, 2x,  & 
                  F10.4, 2F8.1,  & 
                  I5, I3 )
       endif

       OPEN(UNIT= LUNSPEC, FILE=SPEC_FILE, STATUS='UNKNOWN' )

! write comments at top of SPECTRUM file

       MSG = 'VERSION_PHOTOMETRY: ' // VERS(1:LENV)
	 write(LUNSPEC,350) MSG

       MSG = 'CID:      ' // SNLC_CCID(1:LENC)  ! Oct 29 2025
	 write(LUNSPEC,350) MSG

       write(MSG,301) 'MJD:      ', MJD
	 write(LUNSPEC,350) MSG

       write(MSG,301) 'TOBS:     ', TOBS
	 write(LUNSPEC,350) MSG

       write(MSG,301) 'REDSHIFT: ', z
	 write(LUNSPEC,350) MSG

301      format(A, 3x, F9.3)

       write(MSG,302) 'NLAMBIN:  ', NLAMBIN
302      format(A, 3x, I5)
	 write(LUNSPEC,350) MSG

       write(LUNSPEC,350) ' '

       if ( WRSPEC_FLAM ) then
     	      write(MSG,360) 'WAVE FLAM FLAMERR'
       else if ( WRSPEC_FLUX ) then
            write(MSG,360)	'WAVE FLUX FLUXERR'
	 else if ( WRSPEC_KEY ) then
            write(MSG,351)	'ROW WAVE FLAM FLAMERR'
       endif

       if ( LSIM_SNANA ) MSG = MSG(1:35) // '  SIM_FLAM  '
       write(LUNSPEC,360) MSG(1:50)

350      format('# ', A)
351      format('VARNAMES: ', A)
360      format(A)

! write spectrum
       DO 401 ilam = 1, NLAMBIN
         LAM = 0.5* ( LAMMIN_SPECTRUM(ilam) + LAMMAX_SPECTRUM(ilam) )
         LAMBIN = LAMMAX_SPECTRUM(ilam) - LAMMIN_SPECTRUM(ilam)
         F      = FLAM_SPECTRUM(ilam)
         FERR   = FLAMERR_SPECTRUM(ilam)
         if ( WRSPEC_FLUX ) then
            F    = F * LAMBIN
            FERR = FERR * LAMBIN
         endif

         if ( WRSPEC_FLAM .or. WRSPEC_FLUX ) THEN
            WRITE(MSG,440) '', LAM, F, FERR
440           format(A, F8.2, 3x, 2E14.5 )
            LENROW = 42
         else if ( WRSPEC_KEY ) THEN
            WRITE(MSG,442) 'ROW: ', ilam, LAM, F, FERR
442           format(A, I5, 2x, F8.2, 3x, 2E14.5 )
            LENROW = 60
         endif

	   if ( LSIM_SNANA) then
	      FSIM = SIM_FLAM_SPECTRUM(ilam)
	      MSGTMP = MSG
	      write(MSG,446) MSGTMP(1:LENROW), FSIM
446           format(A,E14.5)
         endif

	   write(LUNSPEC,360) MSG

401      CONTINUE

       CLOSE(LUNSPEC)  ! close spec file
400   CONTINUE

    CLOSE ( UNIT = LUNTMP )

    RETURN
  END SUBROUTINE WR_SPECTRA

! ===========================
    SUBROUTINE  GET_SIMTYPE_SPEC_FILE(SIM_TYPE, SIM_SUBTYPE)

! Created Oct 23 2025
! Return TYPE and SUBTYPE based on SIMNAMEP_TYPE;
! This really should be read from a map.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

    CHARACTER SIM_TYPE*20, SIM_SUBTYPE*20

! ---------- BEGIN -----------

    SIM_TYPE = SIMNAME_TYPE
    if ( INDEX(SIMNAME_TYPE,'SALT')   > 0 ) SIM_TYPE = 'Ia'
    if ( INDEX(SIMNAME_TYPE,'BAYESN') > 0 ) SIM_TYPE = 'Ia'

! get subtypes for SNID-SAGE only using hard-wired strings sent
! by Fiore on Oct 23 2025

    if ( WRSPEC_SNID_SAGE ) then

       if ( SIM_TYPE(1:2) .EQ. 'Ia' ) SIM_SUBTYPE = 'Ia-norm'

       if ( SIMNAME_TYPE .EQ. 'IIP' ) then
           SIM_TYPE = 'II' ; SIM_SUBTYPE = 'IIP'
       endif

       if ( SIMNAME_TYPE .EQ. 'IIL' ) then
           SIM_TYPE = 'II' ; SIM_SUBTYPE = 'IIL'
       endif

       if ( SIMNAME_TYPE .EQ. 'Ib' ) then
           SIM_TYPE = 'Ib' ; SIM_SUBTYPE = 'Ib-norm'
       endif

       if ( SIMNAME_TYPE .EQ. 'Ic' ) then
           SIM_TYPE = 'Ic' ; SIM_SUBTYPE = 'Ic-norm'
       endif

       if ( SIMNAME_TYPE .EQ. 'iax' ) then
           SIM_TYPE = 'Ia' ; SIM_SUBTYPE = 'Ia-02cx'
       endif
    endif
    RETURN
  END SUBROUTINE  GET_SIMTYPE_SPEC_FILE

! ===============================
    SUBROUTINE END_REFORMAT

! 
! Close FITS files and create README file.
! 
! Jul 6 2021: don't write anything else to README so that
!             original DOCUMENTATION block appears at the top.
!             (maybe later, call function to append yaml format)
! 
! Mar 7 2022: do explicit gzip
! 
! -------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER LEN_VERS, iver, OPTMASK
    CHARACTER  & 
         README_FILE*(MXCHAR_FILENAME)  & 
        ,CMD*(2*MXCHAR_FILENAME)

! ----------------- BEGIN ----------------

! call C function to close out the FITS files.
    IF ( OPT_REFORMAT_FITS > 0 ) THEN
       OPTMASK = 0  ! don't gzip because gzip logic doesn't work
       CALL wr_snfitsio_end(OPTMASK)

       LEN_VERS  = INDEX(VERSION_REFORMAT_FITS,' ' ) - 1
       CMD = 'cd ' // VERSION_REFORMAT_FITS(1:LEN_VERS) //  & 
              ' ; gzip *.FITS'
       print*,'   gzip output FITS files in ',  & 
              VERSION_REFORMAT_FITS(1:LEN_VERS)
       call flush(6)
       CALL SYSTEM(CMD)
    ENDIF


! create README blank readme file
    CALL OPEN_REFORMAT_FILE(LUNTMP, 'README', README_FILE)
    CLOSE(LUNTMP)  ! close README file

! --- catenate original README contents for each version --------

    DO iver = 1, N_VERSION

! now just catenate the final readme with the original readme file.

       LEN_VERS  = INDEX(SNREADME_FILE(iver),' ' ) - 1
       CMD = 'cat ' // SNREADME_FILE(iver)(1:LEN_VERS)  & 
               // ' >> ' // README_FILE

       CALL SYSTEM(CMD)

    ENDDO

! ---- update IGNORE file with comment on total number of IGNORE epochs
    write(LUNIGNORE2,430) NEPOCH_IGNORE_WRFITS
 430  format(/,'Total number of IGNORE epochs: ', I5)
    close(LUNIGNORE2)

    IF ( OPT_REFORMAT_FITS > 0 ) then
       print*,' Done writing reformatted FITS version: ',  & 
              VERSION_REFORMAT_FITS
    ELSE
       print*,' Done writing reformatted TEXT version: ',  & 
              VERSION_REFORMAT_TEXT
    ENDIF

    RETURN
  END SUBROUTINE END_REFORMAT


! =========================================
    SUBROUTINE WRSALT2_2()
! 
! Created July 30, 2010 by R.Kessler
! Re-write SN data files into new SALT2 format with
! one file per SN.
! Output used by J. Guy's SALT2 'snfit' and 'pcafit' programs,
! version 2.3.0 and higher.
! Mainly used for SDSS data, and for simulations.
! 
! 
! See &SNLCINP namelist variables
!  OPT_REFORMAT_SALT2 = 2
! 
! so that arbitray keys can be specified.
! 
! 
! Jan 25, 2013: use REAL*8 MJD8 instead of REAL*4 MJD to avoid
!               round-off error.
! 
! May 20 2013: include RA and DEC in SALT2 files.
! 
! Jul 8 2013: write name of ascii SNDATA file ONLY for FORMAT_TEXT=T
! 
! Jul 17, 2013: write @FIELD  <list of fields>
! 
! Feb 25, 2014: remove logic with ITERSE_FLUXZPT
! -----------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE WRS2COM
    USE PKMJDCOM

    IMPLICIT NONE

    INTEGER  & 
         EPMIN, EPMAX, EP, NEWMJD  & 
        ,IFILT, IFILT_OBS, IFILT_REPLACE  & 
        ,LENCID, LS, IERR, CID, ikey, L1, L2, LL, i

    REAL  & 
         PEAKMJD, FLUX, FLUXERR  & 
        ,ZP, ZPOFF, arg, Fscale, TMP, ERR


    REAL*8  MJD8

    LOGICAL LWR_FILT(MXFILT_OBS)

    CHARACTER  & 
         cfilt1*1,  ccid*(MXCHAR_CCID)  & 
        ,salt2_file*(MXCHAR_FILENAME)  & 
        ,SURVEY_NAME_LOCAL*(MXCHAR_FILEWORD)  & 
        ,SIMFILE*(MXCHAR_FILENAME)  & 
        ,ctmp*60

! ------------ BEGIN -------------

    CCID   = SNLC_CCID

! if CCID is an integer, write in I6.6 (or I8.8) format
    read( ccid, 25 , iostat = IERR ) CID
25    format(I9)

    IF ( IERR .EQ. 0  ) THEN
       CALL CIDSTRING(CID,CCID,LENCID)
    ENDIF

!  ----------------------

    LS  = INDEX(SURVEY_NAME,' ' ) - 1
    LENCID = INDEX(CCID,' ' ) - 1

    SALT2_FILE = NAMEof_PREFIX(1:LEN_PREFIX) // '_' //  & 
                   CCID(1:LENCID) // '.DAT'

     LL = INDEX(SALT2_FILE,' ') - 1
     WRITE(LUNSALT2,700) SALT2_FILE(1:LL)
700    FORMAT(A)

    IF ( LSIM_SNANA ) THEN
        PEAKMJD = SIM_PEAKMJD
    ELSE
        PEAKMJD = SNLC_SEARCH_PEAKMJD
    ENDIF

    OPEN(   UNIT   = LUNDAT  & 
            , FILE   = SALT2_FILE  & 
            , STATUS = 'UNKNOWN'  & 
                 )

    IF ( NAMEof_SURVEY(1:4) .EQ.  'NULL' ) then
      SURVEY_NAME_LOCAL = SURVEY_NAME  ! SNANA survey name
    ELSE
      SURVEY_NAME_LOCAL = NAMEof_SURVEY ! user-specified survey name
    ENDIF


    ctmp = 'Translated into SALT2 format by SNANA ' // SNANA_VERSION
    write(LUNDAT,120) ctmp
    write(LUNDAT,120) ' '

    write(LUNDAT,120) 'Required variables'
    write(LUNDAT,11) 'SURVEY',   SURVEY_NAME_LOCAL
    write(LUNDAT,11) 'SN',       CCID(1:LENCID)

    write(LUNDAT,126) 'RA',    SNLC8_RA
    write(LUNDAT,126) 'DEC',   SNLC8_DEC

! ---- write field(s) - include all overlapping fields for this SN -----
    write(LUNDAT,1220)
    DO  i = 1, ISNLC_NFIELD_OVP
       LL = INDEX( SNLC_FIELD_OVPLIST(i), ' ') - 1
       write(LUNDAT,1230) SNLC_FIELD_OVPLIST(i)(1:LL)
    ENDDO
    write(LUNDAT,1240)

 1220 FORMAT('@FIELD  ', $)
 1230 FORMAT(A, '  ', $)
 1240 FORMAT(' ')

! ------------------------------------

    write(LUNDAT,12) 'Z_HELIO',  SNLC_ZHELIO
    write(LUNDAT,12) 'Z_CMB  ',  SNLC_ZCMB

    write(LUNDAT,12) 'MWEBV  ',  SNLC_MWEBV
    write(LUNDAT,12) 'MWEBV_ERR ',SNLC_MWEBV_ERR
    write(LUNDAT,13) 'DayMax  ', PEAKMJD

! --------------------------------------
! add new keys, if any are given (Dec 2010)

    DO 22 ikey = 1, NEWKEY
       L1 = INDEX(NEWKEY_NAME(ikey),' ') - 1
       L2 = INDEX(NEWKEY_ARG(ikey),' ') - 1
       write(LUNDAT,11)  & 
            NEWKEY_NAME(ikey)(1:L1), NEWKEY_ARG(ikey)(1:L2)
22    CONTINUE


! Set LWR_FILT logical array for filters to write out.
! Use SNRMAX_FILT to determine if a filter is used.

    DO IFILT     = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       TMP       = SNLC_SNRMAX_FILT(ifilt)
       if ( TMP .NE. -9.0 ) then
          LWR_FILT(IFILT_OBS) = .TRUE.
       else
          LWR_FILT(IFILT_OBS) = .FALSE.
       endif
    ENDDO

! --------------------------------------
! add extra cut variables (Sep 29 2010)

    write(LUNDAT,120) ' '
    write(LUNDAT,120) 'Analysis variables'
    write(LUNDAT,30) 'SNTYPE   ',  ISNLC_TYPE
    write(LUNDAT,30) 'NOBS     ',  ISNLC_NEPOCH_USE
    write(LUNDAT,14) 'TRESTMIN ',  SNLC_TRESTMIN
    write(LUNDAT,14) 'TRESTMAX ',  SNLC_TRESTMAX
    write(LUNDAT,14) 'T0GAPMAX ',  SNLC_T0GAPMAX
    write(LUNDAT,14) 'SNRMAX   ',  SNLC_SNRMAX_SORT(1)
    write(LUNDAT,14) 'SNRMAX2  ',  SNLC_SNRMAX_SORT(2)
    write(LUNDAT,14) 'SNRMAX3  ',  SNLC_SNRMAX_SORT(3)

! write SNRMAX for each filter.
    DO IFILT     = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
       TMP       = SNLC_SNRMAX_FILT(ifilt)
       if ( LWR_FILT(ifilt_obs) ) then
          write(LUNDAT,114) 'SNRMAX_', cfilt1, TMP
       endif

    ENDDO

114   format('@', A, A, 2x, F8.3)

! to quanititative FoM, write error on flux and T0 for each band
    IF ( OPT_SETPKMJD > 0 ) THEN

      write(LUNDAT,120) ' '
      write(LUNDAT,120) 'Fitted errors => exp(dT/T1)/[1+exp(dT/T2)]'
      write(LUNDAT,120) 'ERRT0 => error on peak MJD'
      write(LUNDAT,120) 'ERRF0 => ERROR(F0)/F0 where F0 = peak flux'


      DO IFILT    = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
        if ( LWR_FILT(ifilt_obs) ) then
          ERR       = PKMJD_ERR(ifilt_obs)
          write(LUNDAT,114) 'ERRT0_', cfilt1, ERR
        endif
      ENDDO
        write(LUNDAT,12) 'ERRT0_MIN', PKMJD_ERRMIN
        write(LUNDAT,12) 'ERRT0_WGT', PKMJD_ERRWGT

      DO IFILT    = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
        if ( LWR_FILT(ifilt_obs) ) then
          ERR       = PKFLUX_ERR(ifilt_obs)/PKFLUX_FIT(ifilt_obs)
          write(LUNDAT,114) 'ERRF0_', cfilt1, ERR
        endif
      ENDDO
        write(LUNDAT,12) 'ERRF0_MIN', PKFLUX_ERRMIN
        write(LUNDAT,12) 'ERRF0_WGT', PKFLUX_ERRWGT

    ENDIF  ! end of OPT_SETPKMJD


    IF ( LSIM_SNANA ) THEN
      write(LUNDAT,120) ' '
      write(LUNDAT,120) 'SIMULATION Parameters'

! start with the exact name of the SNANA file that was translated
      IF ( FORMAT_TEXT ) THEN
        L1 = INDEX(VERSION_PHOTOMETRY(1),' ' ) - 1
        L2 = INDEX(SNDATA_FILE_CURRENT,' ' ) - 1
        SIMFILE = '$SNDATA_ROOT/SIM/' //  & 
              VERSION_PHOTOMETRY(1)(1:L1) // '/' //  & 
              SNDATA_FILE_CURRENT(1:L2)
        LL = INDEX(SIMFILE,' ') - 1
        WRITE(LUNDAT,11) 'SIM_SNANAFILE' , SIMFILE(1:LL)
      ENDIF

      write(LUNDAT,30)  'SIM_SNTYPE   ', SIM_GENTYPE
      write(LUNDAT,30)  'SIM_TEMPLATE_INDEX  ', SIM_TEMPLATE_INDEX
      write(LUNDAT,93)  'SIM_PEAKMJD  ', SIM_PEAKMJD
      write(LUNDAT,12)  'SIM_MWEBV    ', SIM_MWEBV
      write(LUNDAT,12)  'SIM_REDSHIFT ', SIM_REDSHIFT_CMB
      write(LUNDAT,12)  'SIM_MU       ', SIM_DLMAG

      IF ( SIM_TEMPLATE_INDEX .EQ. 0 ) THEN
        write(LUNDAT,12)  'SIM_x0 ',    SIM_SALT2x0
        write(LUNDAT,12)  'SIM_mb ',    SIM_SALT2mb
        write(LUNDAT,12)  'SIM_x1 ',    SIM_SHAPEPAR
        write(LUNDAT,12)  'SIM_c  ',    SIM_COLORPAR
        write(LUNDAT,12)  'SIM_alpha ', SIM_SHAPELAW
        write(LUNDAT,12)  'SIM_beta  ', SIM_COLORLAW
      ENDIF

      DO IFILT    = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        cfilt1    = filtdef_string(ifilt_obs:ifilt_obs)
        TMP       = SIM_PEAKMAG(ifilt)
        write(LUNDAT,114) 'SIM_PEAKMAG_', cfilt1, TMP
      ENDDO

    ENDIF  ! end of LSIM_SNANA

    write(LUNDAT,120) ' '
    write(LUNDAT,120) '--------------------------------------'
    write(LUNDAT,120) ' '

! ---------------------------------------
! define data columns
    write(LUNDAT,20) 'Date'
    write(LUNDAT,20) 'Flux'
    write(LUNDAT,20) 'Fluxerr'
    write(LUNDAT,20) 'ZP'
    write(LUNDAT,20) 'Filter'
    write(LUNDAT,20) 'MagSys'
    write(LUNDAT,20) 'end'

11    format('@', A, 2x, A)
12    format('@', A, 2x, F8.5)
13    format('@', A, 2x, F9.3, 2x, '1')
14    format('@', A, 2x, F7.2)
20    format('#', A, ' : ' )
30    format('@', A, 2x, I5)
93    format('@', A, 2x, F9.3)
126   format('@', A, 2x, F12.6)

120   format('# ! ', A )

! write observations.

    DO 301 NEWMJD = 1, ISNLC_NEWMJD_STORE
        EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
        EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

        DO 302 EP = EPMIN, EPMAX

          IFILT_OBS     = ISNLC_IFILT_OBS(ep)
          IFILT_REPLACE = IMAP_REPLACE(IFILT_OBS)

!  skip obs that are cut in SNRECON (Jan 2, 2012)
          if ( .not. LWR_FILT(ifilt_obs) ) GOTO 302


          IF ( IFILT_REPLACE .GT. 0 ) THEN
             cfilt1 = filtdef_string(ifilt_replace:ifilt_replace)
          ELSE
             cfilt1 = filtdef_string(ifilt_obs:ifilt_obs)
          ENDIF

          ZPOFF   = MAGOBS_SHIFT_ZP_FILT(ifilt_obs)
          MJD8    = SNLC8_MJD(EP)

          ZP  = SNGL(ZP_FLUXCAL) + ZPOFF
          arg = 0.4*ZPOFF

          Fscale  = 10.0**ARG

! get fluxes in native system (i.e., undo AB offsets)
          Flux    = Fscale * SNLC_FLUXCAL(EP)
          Fluxerr = Fscale * SNLC_FLUXCAL_ERRTOT(EP)
          WRITE(LUNDAT,330) MJD8, FLUX, FLUXERR, ZP  & 
                , NAMEof_INSTRUMENT(1:LEN_INST)  & 
                , cfilt1  & 
                , NAMEof_MAGSYS(1:LEN_MAGSYS)

302      CONTINUE  ! EP
301   CONTINUE  ! NEWMJD

330   FORMAT(F9.3, 2x, 2G14.6, 1x, F7.3, 3x,A,'::',A,1x,A )

    CLOSE ( LUNDAT )

    RETURN
  END SUBROUTINE WRSALT2_2


! ============================
    SUBROUTINE MADABORT(funname,c1err,c2err)


    USE SNPAR
    USE CTRLCOM

    IMPLICIT NONE

    character funname*(*)  ! name of function calling MADABORT
    character c1err*(*), c2err*(*) ! two error messages

    INTEGER LENF, LEN1, LEN2
    INTEGER, PARAMETER :: ISEV_FATAL = 4, IPROMPT=0 

    EXTERNAL ERRMSG
    INTEGER  LENSTR
! ----------------- BEGIN --------------
    LENF = LENSTR(funname//' ')
    LEN1 = LENSTR(c1err)
    LEN2 = LENSTR(c2err)
    call errmsg(ISEV_FATAL, IPROMPT,  & 
           funname(1:LENF)//char(0),  & 
           c1err(1:LEN1)//char(0),  & 
           c2err(1:LEN2)//char(0),  & 
           LENF, LEN1, LEN2 )

    RETURN
  END SUBROUTINE MADABORT

! =============================
    INTEGER FUNCTION LENSTR(string)
! Created Dec 17 2019
! STRING = 'Hi, what is this.' --> LENSTR = 17
! 

    CHARACTER STRING*(*)
    INTEGER LENTOT, LEN, i
! ---------- BEGIN -------------
    LENTOT = index(string,' ', back = .TRUE. ) - 1
    DO i = 1, LENTOT
      IF ( STRING(i:i) .NE. ' ' ) LEN=i
    ENDDO
    LENSTR = LEN
    RETURN
  END FUNCTION LENSTR


! ======================================
    SUBROUTINE CHECK_LINE_ARGS()
! 
! Feb 7, 2007
! Check that all NLINE_ARG command-line args
! were used. If any were not used, then give
! error message and abort.
! 
! Aug 2020: use MADABORT for consistent error messaging
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, LL, NBAD
    CHARACTER FNAM*16, ARG*(MXCHAR_ARG)
! ------------------ BEGIN --------------

    NBAD = 0
    FNAM = 'CHECK_LINE_ARGS'

    DO i = 2, NLINE_ARGS

       if ( .NOT. USE_LINE_ARGS(i) ) then
          NBAD = NBAD + 1

          if ( NBAD .EQ. 1 ) then
            CALL PRINT_PREABORT_BANNER(FNAM//char(0),20)
          endif
          CALL GETARG(i,ARG)
          LL = INDEX(ARG, ' ' ) - 1
          print*,' ERROR: Detected un-used command-line argument: ',  & 
                       '"' , ARG(1:LL), '"'
       endif
    ENDDO

    IF ( NBAD .GT. 0 ) THEN
       write(C1ERR,61) NBAD
61       format('Detected ', I3,' unknown command line arguments')
       C2ERR = 'Check args ERROR messages above'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

    RETURN
  END SUBROUTINE CHECK_LINE_ARGS


! ================================================
    SUBROUTINE INIT_SNVAR()
! 
! Init some SN variables.
! Called once at beginning of program.
! See INIT_SNLC to initilaize each SN.
! 
! Feb 19 2014: abort if SNDATA_ROOT is not defined.
! Jan 16 2019: init EXIT_ERRCODE
! ---------------------------------------------------


! include everything except the memory hogging SNLCCOM

    USE SNPAR
    USE SNFILECOM
    USE CTRLCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE SNHOSTCOM
! +CDE,PARSECOM.
    USE SNSIMCOM
    USE SNANAFIT
    USE FILTCOM
    USE PKMJDCOM
    USE PRIVCOM
    USE EARLYCOM
    USE SPECCOM

    IMPLICIT NONE

    INTEGER icut, i, ibit, itype, ipar, ifilt, LM,istat
    INTEGER LL, LL_SNANA, NWD, LL_HOST
    CHARACTER ENVtemp*400, MACH*40, FNAM*12

! function
    INTEGER FILTINDX, EXEC_CIDMASK, STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS, SET_EXIT_ERRCODE

! ---------------- BEGIN ------------

    FNAM = 'INIT_SNVAR'

#if defined(SNANA)
    EXIT_ERRCODE = EXIT_ERRCODE_SNANA
#elif defined(SNFIT)
    EXIT_ERRCODE = EXIT_ERRCODE_SNFIT
#elif defined(PSNID)
    EXIT_ERRCODE = EXIT_ERRCODE_PSNID
#endif
    CALL SET_EXIT_ERRCODE(EXIT_ERRCODE) ! for C function aborts

    IDSURVEY = -9 ;  IDSUBSURVEY=-9
    NCALL_SNANA_DRIVER  = 0

    CALL PRBANNER ( " INIT_SNVAR: Init variables." )

    FORMAT_TEXT    = .FALSE.
    FORMAT_FITS    = .FALSE.

    ISJOB_SNANA = .FALSE.
    ISJOB_SNFIT = .FALSE.
    ISJOB_PSNID = .FALSE.
    ISJOB_BATCH = .FALSE.

    NEPOCH_TOT = 0
    NEPOCH_USE = 0
    NEPOCH_CUT = 0
    NEPOCH_BADPHOT     = 0
    NEPOCH_BADPHOT_SUM = 0

    INTERP_OPT  = INTERP_LINEAR  ! used by simulation

! define filter indices for Legacy filters
    IFILT_SDSS_u = FILTINDX("u ")
    IFILT_SDSS_g = FILTINDX("g ")
    IFILT_SDSS_r = FILTINDX("r ")
    IFILT_SDSS_i = FILTINDX("i ")
    IFILT_SDSS_z = FILTINDX("z ")
    IFILT_BESS_U = FILTINDX("U ")
    IFILT_BESS_B = FILTINDX("B ")
    IFILT_BESS_V = FILTINDX("V ")
    IFILT_BESS_R = FILTINDX("R ")
    IFILT_BESS_I = FILTINDX("I ")
    IFILT_BESS_BX = FILTINDX("X ")
    IFILT_Y = FILTINDX("Y ")
    IFILT_J = FILTINDX("J ")
    IFILT_H = FILTINDX("H ")
    IFILT_K = FILTINDX("K ")

    SURVEY_NAME    = ''
    SUBSURVEY_NAME = ''
    SUBSURVEY_NAME_LIST=''

    SURVEY_FILTERS = ''
    FREEZE_SURVEY_FILTERS = .FALSE.
    EXIST_BXFILT_OBS      = .FALSE.
    EXIST_BXFILT_REST     = .FALSE.

    NFILTDEF_SURVEY = 0
    NFILTDEF_READ   = 0
    NFILTDEF_REST   = 0

    NFILTDEF_IGNORE_REST = 0

    DO ifilt = 1, MXFILT_OBS
        IFILTDEF_MAP_SURVEY(ifilt) = -9
        IFILTDEF_MAP_REST(ifilt)   = -9
        IFILTDEF_IGNORE_REST(ifilt) = -9
    ENDDO

    DO ifilt = 1, MXFILT_ALL

        IFILTDEF_INVMAP_SURVEY(ifilt)  = -9
        IFILTDEF_INVMAP_REST(ifilt)    = -9

        MAGOBS_SHIFT_PRIMARY_FILT(ifilt)  = 0.0
        MAGOBS_SHIFT_ZP_FILT(ifilt)       = 0.0
        MAGOBS_SHIFT_ZP_USER(ifilt)       = 0.0
        MAGREST_SHIFT_PRIMARY_FILT(ifilt) = 0.0
        MAGREST_SHIFT_ZP_FILT(ifilt)      = 0.0
        MAGREST_SHIFT_ZP_USER(ifilt)      = 0.0
        FILTOBS_ZPOFF_SNPHOT(ifilt)       = 0.0

        FILTOBS_TRANSMAX(ifilt)         =  0.0
        FILTOBS_LAMAVG(ifilt)           = -9.0
        FILTOBS_LAMRMS(ifilt)           = -9.0
        FILTOBS_LAMRANGE(1,ifilt)       = -9.0
        FILTOBS_LAMRANGE(2,ifilt)       = -9.0
! xxx mark           FILTOBS_LAMSHIFT(ifilt)         =  0.0
        NLAMBIN_FILTOBS(ifilt)          =  0

        FILTREST_TRANSMAX(ifilt)         =  0.0
        FILTREST_LAMAVG(ifilt)           = -9.0
        NLAMBIN_FILTREST(ifilt)          =  0
        FILTREST_LAMRANGE(1,ifilt)       = -9.0
        FILTREST_LAMRANGE(2,ifilt)       = -9.0

        FILTOBS_NAME(ifilt)       = ''
        SURVEY_FILT_NAME(ifilt)   = 'NONE'
        LFILTDEF_OBS(ifilt)       = .FALSE.
        LFILTDEF_REST(ifilt)      = .FALSE.
        LFILTDEF_NONSURVEY(ifilt) = .FALSE.
        LFILTDEF_SNRMAX(ifilt)    = .TRUE.  ! use all filters by default

        FUDGE_FLUXCAL_ERROR_FILT(ifilt)   =  0.0
        FUDGE_FLUXCAL_ERRPIX_FILT(ifilt)  =  0.0
        FUDGE_FLUXCAL_OFFSET_FILT(ifilt)  =  0.0
        FUDGE_FLUXERR_SCALE_FILT(ifilt)   =  1.0 ! Feb 2020
        FUDGE_MAG_ERROR_FILT(ifilt)       =  0.0
        FUDGE_MAG_COVERR_FILT(ifilt)      =  0.0
        SIM_FUDGE_MAG_ERROR_FILT(ifilt)   =  0.0

    END DO

! extract SNDATA_ROOT and SNANA_DIR  env variables

    CALL GETENV ( 'HOST',     ENVtemp )
    LL = INDEX ( ENVtemp , ' ' ) - 1
    HOST_MACHINE = ENVtemp(1:LL)

    CALL GETENV ( 'MACHTYPE', ENVtemp )
    LM = INDEX ( ENVtemp , ' ' ) - 1
    MACH = ENVtemp(1:LM)

    print*,'  HOST MACHINE = ', HOST_MACHINE(1:LL),  & 
            '   (',MACH(1:LM),')'

    CALL GETENV ( 'SNDATA_ROOT', ENVtemp )
    LL = INDEX ( ENVtemp , ' ' ) - 1

    IF ( LL > MXCHAR_PATH ) then
       write(C1ERR,61) 'SNDATA_ROOT', LL
       write(C2ERR,62) MXCHAR_PATH
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF
61       FORMAT('LEN($', A, ') = ',I3,' is too long')
62       format('MXCHAR_PATH = ', I3)
    if ( LL .EQ. 0 ) then
       C1ERR = 'SNDATA_ROOT is not defined.'
       C2ERR = 'Check your snana setup.'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    SNDATA_ROOT = ENVtemp(1:LL)
    print*,'  SNDATA_ROOT = ', SNDATA_ROOT(1:LL)

! - - - - -
    CALL GETENV ( 'SNANA_DIR', ENVtemp )
    LL = INDEX ( ENVtemp , ' ' ) - 1
    IF ( LL > MXCHAR_PATH ) then
       write(C1ERR,61) 'SNANA_DIR', LL
       write(C2ERR,62) MXCHAR_PATH
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

    SNANA_DIR = ENVtemp(1:LL)
    print*,'  SNANA_DIR = ', SNANA_DIR(1:LL)


    N_SNFILE      = 0
    N_SNFILE_LAST = 0

    N_SNLC_PROC     = 0
    N_SNLC_CUTS     = 0
    N_SNLC_SPEC     = 0
    N_SNLC_FIT      = 0
    N_SNLC_FITCUTS  = 0
    N_SNLC_COVFIX   = 0
    N_SNLC_PLOT     = 0
    N_SNHOST_ZSPEC  = 0
    N_SNHOST_ZPHOT  = 0
    N_DUPLICATE_CID = 0
    N_DUPLICATE_MJD = 0
    NSTORE_DUPLICATE_MJD = 0

    do i = 0, MXMASK_zSOURCE
      N_MASK_zSOURCE_LC_CUTS(i)    = 0
      N_MASK_zSOURCE_LCFIT_CUTS(i) = 0
    enddo

    ABSO_OFFSET = 0
    ABSO_INDEX  = 0

    DO itype = -1, MXTYPE
      NPASSCUT_FIT(itype) = 0
    DO ibit  = 1, NCUTBIT_SNLC+1
      NPASSCUT_INCREMENT(itype,ibit) = 0
    ENDDO
    ENDDO

    LSIM_SNANA   = .FALSE.  ! snana sim
    LSIM_MAGOBS  = .FALSE.
    ISJOB_SIM    = .FALSE.

    NACCEPT_Z    = 0
    NACCEPT_ZERR = 0
    NACCEPT_CID  = 0
    NACCEPT_TYPE = 0
    DO icut = 1, NCUTBIT
      NACCEPT_CUT(icut) = 0
    ENDDO

! init HIDKCOR array to 0

    ISTAT = EXEC_CIDMASK(0,MXCID_CHECK)  ! init CID bit-mask

    DO IPAR = 1, MXFITPAR
       ERRMAX_BAD(ipar) = 1.0E-4  ! err < ERRMAX_BAD is labelled bad
    ENDDO

    SIMNAME_MODEL        = 'SIMNULL0'
    SIMNAME_SHAPEPAR     = 'SIM_SHAPE'  ! changed from SIMSHAPE (4/11/2017)
    SIMNAME_SHAPELAW     = 'SIMNULL1'  ! each SIMNULL has unique name
    SIMNAME_COLORPAR     = 'SIMNULL2'
    SIMNAME_COLORLAW     = 'SIMNULL3'
    SIMNAME_SNRMON       = 'NONE'

    SIMOPT_MWCOLORLAW = -9
    SIMOPT_MWEBV      = -9

    DO ifilt = 1, MXFILT_OBS
       SIM_EXPOSURE_TIME(ifilt)  = 1.0
    ENDDO

    SIMLIB_FILENAME = ""
    SIMLIB_MSKOPT   = 0

    DO ibit = 1, NCUTBIT
       CUTWIN_VAR(1,ibit) = -CUTVAL_OPEN
       CUTWIN_VAR(2,ibit) = +CUTVAL_OPEN
    ENDDO

    NEPOCH_IGNORE =  0
    EPOCH_IGNORE_LASTFILE = ''

!       ZPOFF_FILE    = 'NULL'  ! NULL => use offsets in hid 292
    USERTAGS_FILE = ''
    VPEC_FILE            = ''
    HEADER_OVERRIDE_FILE = ''
    SIM_HEADER_OVERRIDE_FILE = ''


    NPAR_SIMSED      = 0
    NPAR_PySEDMODEL  = 0
    NPAR_LCLIB       = 0

    NPAR_SIM_HOSTLIB = 0

    C1ERR = ' '
    C2ERR = ' '

    DO_FLUXERRCALC = .FALSE.

    NFIT_PKMJD = 0

    NVAR_PRIVATE = 0
    DO i = 1, MXVAR_PRIVATE
      PRIVATE_KEYWORD(i)  = 'UNKNOWN'
      PRIVATE_VARNAME(i)  = 'UNKNOWN'
      PRIVATE_VALUE(i)    = PRIVATE_NULL
      PRIVATE_CUTWIN(1,i) = -CUTVAL_OPEN
      PRIVATE_CUTWIN(2,i) = +CUTVAL_OPEN
    ENDDO

    NCUT_PRIVATE = 0
    DO i = 1, MXCUT_PRIVATE
       PRIVATE_CUTWIN_STRING(i) = ''
       USE_PRIVATE_CUTWIN(i)    = 0
    ENDDO
    PRIVATE_VARNAME_READLIST = ''

    SIMVAR_CUTWIN_STRING  = ''
    NSIMVAR_CUTWIN        = 0

    EARLYLC_STRING        = ''
    REQUIRE_EPOCHS_STRING =  ''
    DUMP_STRING           = ''

! init HOST logicals
    EXIST_SNHOST_ZPHOT    = .FALSE.
    EXIST_SNHOST_LOGMASS  = .FALSE.
    EXIST_SNHOST_LOGSFR   = .FALSE.
    EXIST_SNHOST_LOGsSFR  = .FALSE.
    EXIST_SNHOST_COLOR    = .FALSE.
    EXIST_SNHOST_MAGOBS   = .FALSE.
    EXIST_SNHOST_SB       = .FALSE.
    EXIST_SNHOST_ANGSEP   = .FALSE.
    EXIST_SNHOST_DDLR     = .FALSE.
    EXIST_SNHOST_CONFUSION= .FALSE.

! misc.
    NONSURVEY_FILTERS = ''
    FILTER_REPLACE    = ''
    FILTLIST_LAMSHIFT = ''

    DO i = 1, MXFITSTORE
       PARNAME_STORE(i) = ''
       PAROPT_STORE(i)  = 0
    ENDDO

    CALL FLUSH(6)

    NWD = STORE_PARSE_WORDS(-1, ""//char(0), FNAM//char(0), 2, 12);


    RETURN
  END SUBROUTINE INIT_SNVAR


! =======================================
    SUBROUTINE INIT_SNCID_LISTS()

! Created Jun 9 2021
! Check CID lists from:
!   + integer list in NML (use and ignore)
!   + char list in NML    (use and ignore)
!   + list in separate file
! 
! For string lists that can be very long (from list file),
! use match_cidlist_init[exec] utility in sntools.c[h];
! this match utility uses hash table for faster access.
! 
! Jun 26 2025: set CUTWIN_CID=0 for any input list including SNCID_LIST_FILE
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, LENF, LENV, LENC, NCCID_TMP, NCCID_FILE, OPTMASK
    CHARACTER FNAM*20, cFILE*(MXCHAR_FILENAME)
    CHARACTER CCID*(MXCHAR_CCID), cCCID*(MXCHAR_CCID)
    CHARACTER cBLANK*4, cVARLIST_STORE*60
    LOGICAL IS_CCID_LIST, IS_CIDINT_LIST, IS_FILE, IS_LIST

    LOGICAL  VALIDFILENAME
    INTEGER  MATCH_CIDLIST_INIT, PARSE_NML_STRINGLIST
    EXTERNAL MATCH_CIDLIST_INIT

! ----------- BEGIN ------------

    IF ( SNCID_LIST_FILE(1:20) .EQ. 'HEADER_OVERRIDE_FILE' ) THEN
        SNCID_LIST_FILE = HEADER_OVERRIDE_FILE
    ENDIF

    cBLANK = ""//char(0)

    CALL ENVreplace(SNCID_LIST_FILE)  ! Mar 30 2025

    LENF    = INDEX(SNCID_LIST_FILE,' ') - 1
    cFILE   = SNCID_LIST_FILE(1:LENF) // char(0)

! check input sources of CID list:
    IS_FILE        = VALIDFILENAME(SNCID_LIST_FILE)  ! from file
    IS_CCID_LIST   = ( SNCCID_LIST(1) .NE. ' ' )     ! user input list of strings
    IS_CIDINT_LIST = SNCID_LIST(1) > 0               ! user input list of integers

    IS_LIST = IS_FILE .or. IS_CCID_LIST .or. IS_CIDINT_LIST
    NCID_LIST = 0; NCCID_LIST = 0; NCCID_FILE = 0

! - - - - - - -

    FNAM    = 'INIT_SNCID_LISTS'
    CALL PRBANNER(FNAM)

! disable CUTWIN_CID and prescale for user-specified list,
! but allow these features for list read from file.
    IF ( IS_LIST ) THEN
       CUTWIN_CID(1) = 0;  CUTWIN_CID(2) = 0
	 SIM_PRESCALE = 1  ! disable prescale Dec 11 2024
    ENDIF

!      print*,' xxx IS_CCID/IS_CIDINT/FILE = ',
!     &    IS_CCID_LIST, IS_CIDINT_LIST, IS_FILE
!      print*,' xxx CUTWIN_CID = ', CUTWIN_CID

! - - - - - - -
! check integer lists the old way
    IF ( IS_CIDINT_LIST ) THEN
       i = 1
       DO WHILE ( SNCID_LIST(i) .GE. 1 )
          i = i + 1
       END DO
       NCID_LIST = i - 1
    ENDIF

    IF ( IS_CCID_LIST ) THEN
       OPTMASK = 1  ! match CID only

!  check SNCCID_LIST in &SNLCINP ...
       i = 1; CCID = SNCCID_LIST(i)
       DO WHILE( CCID .NE. ' ' )
          CALL NMLCHECK_SNCCID("SNCCID_LIST",CCID)
          LENC = INDEX(CCID,' ') - 1
          cCCID = CCID(1:LENC) // char(0)
          NCCID_TMP = MATCH_CIDLIST_INIT(cCCID,OPTMASK,cBLANK,LENC,4)
          i = i + 1 ; CCID = SNCCID_LIST(i)
       END DO
       NCCID_LIST = i - 1
    ENDIF

! check CIDs from separate file (table format, or unformatted)
    IF ( IS_FILE ) THEN
       if ( NFIT_ITERATION > 0 .and. OPT_SNCID_LIST>1 ) then
           ! store SALT2 fit pars for INIVAL in snlc_fit evt sync
          cVARLIST_STORE = 'BANDLIST,zHD,PKMJD,x0,x1,c' // ',' //  & 
       	              'zHDERR,PKMJDERR,x0ERR,x1ERR,cERR' // char(0)
       else
          cVARLIST_STORE = cBLANK
       endif
       OPTMASK = 1+8 ! send reset for AUTOSTORE
       NCCID_FILE =  & 
            MATCH_CIDLIST_INIT(cFILE,OPTMASK,cVARLIST_STORE, LENF, 60)
    ENDIF

    write(6,31) NCID_LIST,  'integer SNCID_LIST'
    write(6,31) NCCID_LIST, 'char SNCCID_LIST'
    write(6,31) NCCID_FILE, 'file SNCID_LIST_FILE'

 31   format(T5,'Store ',I6,  ' string  CIDs from ', A)


    CALL FLUSH(6)

! - - - - - -
! Check ignore lists, and use brute-force matching since
! these lists should be small.

    i = 1 ; NCID_IGNORE=0
    DO WHILE ( SNCID_IGNORE(i) .GE. 1 )
       i = i + 1
    END DO
    NCID_IGNORE = i - 1

    i = 1
    DO WHILE ( SNCCID_IGNORE(i) .NE. ' ' )
       CALL NMLCHECK_SNCCID("SNCCID_IGNORE",SNCCID_IGNORE(i))
       i = i + 1
       SNCCID_IGNORE_ALL(i) = SNCCID_IGNORE(i)
    END DO
    NCCID_IGNORE = i - 1

! - - - -
    CALL  RDSNIGNORE          ! read IGNORE-list from user file
    NCCID_PLOT = PARSE_NML_STRINGLIST(SNCCID_PLOT,MXCHAR_CCID)

    if ( NCID_IGNORE > 0 .or. NCCID_IGNORE>0 ) then
       write(6,61) NCID_IGNORE,  'SNCID_IGNORE'
       write(6,61) NCCID_IGNORE, 'SNCCID_IGNORE'
 61      format(T5,'IGNORE ',I4,  ' CIDs from ', A)
    endif

    CALL FLUSH(6)


    RETURN
  END SUBROUTINE INIT_SNCID_LISTS


! =======================================
    SUBROUTINE RDFILE_SNMJDLIST
! 
! Created May 2008 by R.Kessler
! 
! Read list of "SNID MJD" to process SNID list,
! so that fitter can interpolate flux at MJD.
! Initial use is for Ryan Foley's analysis with
! KECK spectra for SDSS-discovered SN.
! Second usage is comparing CSP-vs-SDSS photometry (J.Mosher)
! 
! May 28, 2008:
!   add screen-dump to see what's going on.
! 
! Jan 2012:
!   check 'ALL' option to interpolate ALL SN at a MJD
!   After SN are read, a call to UPDSNMJDLIST checks the
!   ALL option.
! 
! May 18, 2012: complete re-write

! ---------


    USE SNDATCOM
    USE SNLCINP_NML
    USE INTERPCM

    IMPLICIT NONE

    INTEGER   iwd, NWD, LEN, N, j, LL0, LL1
    CHARACTER CCID*(MXCHAR_CCID), NAME_forC*(MXCHAR_FILENAME)
    CHARACTER CTMP0*40, CTMP1*40, FNAM*18
    REAL*8  MJD8
    LOGICAL LOK

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ---------- BEGIN ----------

    N_INTERP_MJDLIST    = 0


    IF ( SNMJD_LIST_FILE .EQ. ' ' ) RETURN

    FNAM = 'RDFILE_SNMJDLIST'
    N_INTERP_MJDLIST_DONE = 0
! init arrays
    N_INTERP_MJDLIST       = 0
    DO j = 1, MXINTERP
       INTERP8_MJDLIST(j)  = 0.0
       INTERP_CCIDLIST(j)  = ''
       INTERP_MJDLIST_DONE(j) = .FALSE.
    ENDDO

    CALL PRBANNER(FNAM)

    LEN = INDEX(SNMJD_LIST_FILE,' ') - 1
    print*,' '
    print*,'  Read flux-interpolation instructions from '
    print*,'  ', SNMJD_LIST_FILE(1:LEN)


    NAME_ForC = SNMJD_LIST_FILE(1:LEN)//char(0)
    NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE, NAME_forC,  & 
               FNAM//char(0), LEN, 18 )

    DO 100 iwd = 1, NWD-1, 2

       CALL get_PARSE_WORD_fortran(iwd+0, CTMP0, LL0 )
       CALL get_PARSE_WORD_fortran(iwd+1, CTMP1, LL1 )

       READ ( CTMP0, *) CCID
       READ ( CTMP1, *) MJD8

       N_INTERP_MJDLIST  = N_INTERP_MJDLIST + 1
       N = N_INTERP_MJDLIST

       INTERP_CCIDLIST(N) = CCID
       INTERP8_MJDLIST(N) = MJD8

       LEN = INDEX(CCID,' ') - 1
       write(6,70) CCID, MJD8
70       format(T6,'Will interpolate flux for SN ', A12,  & 
              ' at MJD=', F10.3)

! --- error checking

       LOK = ( MJD8 .EQ. 0.0) .or.  & 
               ( MJD8 .GT. 40000. .and. MJD8 .LT. 70000. )
       if ( .NOT. LOK ) then
          write(c1err,61) MJD8
61          format('Read Invalid MJD = ', F12.3,' from ' )
          c2err = SNMJD_LIST_FILE(1:80)
          CALL  MADABORT(FNAM, C1ERR, C2ERR)
       endif

! abort if we exceed array bound
       IF ( N_INTERP_MJDLIST .GT. MXINTERP ) THEN
          CALL PRINT_PREABORT_BANNER("RDFILE_SNMJDLIST"//char(0),40)
          print*,'   SNMJD_LIST_FILE = ', SNMJD_LIST_FILE
          write(c1err,161) CCID, MJD8
161         format('N_INTERP_MJDLIST exceeds bound at CID=',A8, 2x,  & 
                'MJD=',F9.3)
          write(c2err,162) MXINTERP
162         format('Check MXINTERP=',I5,' and SNMJD_LIST_FILE above.')

         CALL  MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

100   CONTINUE  ! iwd

    RETURN
  END SUBROUTINE RDFILE_SNMJDLIST

! =============================
    SUBROUTINE GET_INTERP_MJDLIST(CCID, NMJD, MJDLIST, ADDFLAG)
! 
! Return NMJD and MJDLIST to interpolate for this CID.
! ADDFLAG = T if this CCID is specifically listed;
! ADDFLAG = F if this CCID is added because of the ALL option.
! 

    USE SNPAR
    USE INTERPCM

    IMPLICIT NONE

! subroutine args
    CHARACTER CCID*(*)  ! (I) SN character id
    INTEGER NMJD        ! (O) number of MJDs to interpolate
    LOGICAL ADDFLAG     ! (O) T => add to process list
    REAL*8  & 
         MJDLIST(MXINTERP) ! (O) list of MJDs

! local var

    LOGICAL LMATCH, LALL
    INTEGER i, L0, L1
    character CCID_tmp*(MXCHAR_CCID)

! --------------- BEGIN --------------

    NMJD    = 0
    ADDFLAG = .FALSE.

    L0   = INDEX(CCID,' ') - 1

    DO 100 i = 1, N_INTERP_MJDLIST
       CCID_tmp = INTERP_CCIDLIST(i)
       L1 = INDEX(CCID_tmp,' ') - 1

       LMATCH = ( L0 .EQ. L1 ) .and.  & 
                  ( CCID(1:L0) .EQ. CCID_tmp(1:L1) )
       LALL    = CCID_tmp(1:3) .EQ. 'ALL'

       IF ( LMATCH ) ADDFLAG = .TRUE.  ! set return arg

!         print*,' xxx compare ', CCID,'  to ', CCID_tmp,
!     &     '  L0,L1=',L0,L1

       if ( LALL .or. LMATCH )then
         NMJD = NMJD + 1
         MJDLIST(NMJD) = INTERP8_MJDLIST(i)

! keep track of which ones are DONE to help know when to quit reading
         if ( .not. INTERP_MJDLIST_DONE(i) ) then
            N_INTERP_MJDLIST_DONE  = N_INTERP_MJDLIST_DONE + 1
            INTERP_MJDLIST_DONE(i) = .TRUE.
         endif

       endif
100   CONTINUE

    RETURN
  END SUBROUTINE GET_INTERP_MJDLIST

! =========================================
    SUBROUTINE ENVreplace(STRING)

! If input STRING starts with '$' then evaluate ENV
! with GETENV and return STRING with substitution.
! 
! There is a C version of this function in sntools.c,
! but passing back the return string to fortran is not
! so easy, so safer to write this one in fortran too.


    USE SNPAR

    IMPLICIT NONE

    CHARACTER STRING*(*)  ! I/O

    INTEGER jslash, LTOT, LENV, LENS
    CHARACTER ENVname*200, ENV*400, SUFFIX*400
    CHARACTER C1ERR*72, C2ERR*72
    CHARACTER FNAM*12

! -------------- BEGIN -------------

    IF ( STRING(1:1) .NE. '$' ) RETURN

    FNAM = 'ENVreplace'

    LTOT    = INDEX(STRING,' ' ) - 1
    IF (LTOT <= 0 ) THEN
      CALL PRINT_PREABORT_BANNER(FNAM(1:10)//char(0), 40)
      print*,' STRING=', STRING
      write(c1err,31) LTOT
      c2err = 'STRING could be too long or zero length'
31      format( 'LEN(STRING)=', I3, ' ???' )
      CALL MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

    JSLASH  = INDEX(STRING,'/') - 1
    IF ( JSLASH <= 0 ) JSLASH = LTOT
    SUFFIX  = STRING(JSLASH+1:LTOT)
    LENS    = INDEX(SUFFIX,' ') - 1

    ENVname = STRING(2:JSLASH)
    CALL GETENV(ENVname, ENV) ! return ENV
    LENV    = INDEX(ENV,' ')-1

    if ( LENV > MXCHAR_PATH ) then
       write(C1ERR,61) ENVname(1:JSLASH-1), LENV
61       format('LEN($', A, ') = ',I3,' is too long.')
       write(C2ERR,62) MXCHAR_PATH
62       format('MXCHAR_PATH = ', I3 )
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    if ( LENV+LENS > MXCHAR_FILENAME ) then
       write(C1ERR,161) STRING(1:LTOT), LENV+LENS
161       format('Decoded LEN($', A, ') = ',I3,' is too long.')
       write(C2ERR,162) MXCHAR_FILENAME
162       format('MXCHAR_FILENAME = ', I3 )
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    if ( LENV .EQ. 0 ) then
       C1ERR = 'ENV = $' // ENVname(1:JSLASH) // ' is not defined.'
       C2ERR = 'Check your input namelists.'
       CALL MADABORT(FNAM, C1ERR, C2ERR)
    endif

    STRING = ENV(1:LENV) // SUFFIX

    RETURN
  END SUBROUTINE ENVreplace

! ==============================
    LOGICAL FUNCTION IGNOREFILE_fortran(fileName)

! Created Apri 2022
! Wrapper to call IGNOREFILE C-utility to return True if
! file = 'NULL' or 'NONE' or 'BLANK' or ''

    USE SNPAR
    CHARACTER FILENAME*(*), cFILE*(MXCHAR_FILENAME)
    INTEGER LENF, IGNORE

    INTEGER  IGNOREFILE
    EXTERNAL IGNOREFILE
! --------- BEGIN -------------

    LENF = INDEX(FILENAME//' ',' ') - 1
    cFILE = FILENAME(1:LENF) // char(0)
    IGNORE = IGNOREFILE(cFILE,LENF)

! xxxxxxxx
!      print*,' xxx ------------------------------- '
!      print*,' xxx FILENAME = ', FILENAME(1:LENF)
!      print*,' xxx LENF = ', LENF
!      print*,' xxx IGNORE = ', IGNORE
! xxxxxxxxxx

    IGNOREFILE_fortran = (IGNORE .NE. 0 )

    RETURN
  END FUNCTION IGNOREFILE_fortran

    LOGICAL FUNCTION VALIDFILENAME(fileName)
    CHARACTER FILENAME*(*)
    LOGICAL   IGNOREFILE_fortran
    VALIDFILENAME = ( .NOT. IGNOREFILE_fortran(FILENAME) )
    RETURN
  END FUNCTION VALIDFILENAME

! =======================================
    SUBROUTINE get_PARSE_WORD_fortran(iwd,WORD,LEN)
! 
! Jan 2018
! Wrapper to call C-function get_PARSE_WORD to retrieve
! iwd word stored.
! 
! Jan 2019: use MXCHAR_FILEWORD for char array abound.
! Jul 2020: extend WORD_C size to 4*MXCHAR_FILENAME
! 

    USE SNPAR

    IMPLICIT NONE

! subroutine args
    INTEGER   IWD           ! (I) word index to fetch
    CHARACTER WORD*(*)      ! (O) stored word
    INTEGER   LEN           ! (O) length of word

! local args
    INTEGER  iwd_C
    CHARACTER WORD_C*(4*MXCHAR_FILENAME)
    EXTERNAL GET_PARSE_WORD

! ------------- BEGIN -----------
    iwd_C = iwd-1
    CALL GET_PARSE_WORD(ONE, iwd_C, WORD_C,  & 
              "get_PARSE_WORD_fortran"//char(0), LEN, 20)
    LEN  = INDEX(WORD_C, ' ' ) - 1
    WORD = WORD_C(1:LEN)
    RETURN
  END SUBROUTINE get_PARSE_WORD_fortran


! ======================
    SUBROUTINE UPCASE(strin, strupper)
    character strin*(*), strupper*(*)
    integer  j, LL

    strupper = ''
    LL = INDEX(strin,' ') - 1
    do j = 1, LL
     if(strin(j:j) >= "a" .and. strin(j:j) <= "z") then
       strupper(j:j) = achar(iachar(strin(j:j)) - 32)
     else
       strupper(j:j) = strin(j:j)
     end if
    end do
  END SUBROUTINE UPCASE

! ======================
    SUBROUTINE LOCASE(strin, strlower)
! return strlower.
    character strin*(*), strlower*(*)
    integer  j, LL

    strlower = ''
    LL = INDEX(strin,' ') - 1
    do j = 1, LL
     if(strin(j:j) >= "A" .and. strin(j:j) <= "Z") then
       strlower(j:j) = achar(iachar(strin(j:j)) + 32)
     else
       strlower(j:j) = strin(j:j)
     end if
    end do
  END SUBROUTINE LOCASE

! ===========================================
    SUBROUTINE RDNINT( cwords, NINT, IARRAY )

! ---------------------------------------------------
! Read NINT consecutive integers from cwords array
! Store results in IARRAY.
! --------------------------------------------------
    INTEGER   NINT
    character cwords(NINT)*60  ! (I) list of char words
    INTEGER   IARRAY(NINT)     ! (O) array of parsed integers
! local args
    INTEGER i
! ------------ BEGIN -----------
    DO i = 1, NINT
        read ( cwords(i), * ) IARRAY(i)
    ENDDO
    RETURN
  END SUBROUTINE RDNINT

! ===========================================
    SUBROUTINE RDNFLOAT ( CWORDS, NFLOAT, FARRAY )

! ---------------------------------------------------
! Read NFLOAT consecutive floats from CWORDS array
! Store results in FARRAY.
! 
! Jul 2008: ABORT on 'nan'.
! Sep 2010: increase cwords bytes from 40 to 60
! 
! --------------------------------------------------
    INTEGER NFLOAT
    character cwords(NFLOAT)*60  ! (I) list of char words
    REAL      FARRAY(NFLOAT)     ! (O) array of parsed floats
! local args
    INTEGER i
    character c1err*60, c2err*60
! ------------ BEGIN -----------
    DO i = 1, NFLOAT
        if ( cwords(i) .EQ. 'nan' ) then
           print*,' '
           c1err = "Found 'nan' instead of float. "
           c2err = "Use grep to find the nan in this version. "
           CALL MADABORT("RDNFLOAT", c1err, c2err )
        endif
        read ( cwords(i), * ) FARRAY(i)

    ENDDO
    RETURN
  END SUBROUTINE RDNFLOAT 

! ===========================================
    SUBROUTINE RDNDBLE ( CWORDS, NDBLE, DARRAY )

! ---------------------------------------------------
! Created Jun 20, 2011
! Read NDBLE consecutive doubles from CWORDS array
! Store results in DARRAY.
! 
! --------------------------------------------------
    INTEGER NDBLE
    character cwords(NDBLE)*60  ! (I) list of char words
    REAL*8    DARRAY(NDBLE)     ! (O) array of parsed floats
! local args
    INTEGER i
    character c1err*60, c2err*60
! ------------ BEGIN -----------
    DO i = 1, NDBLE
        if ( cwords(i) .EQ. 'nan' ) then
           print*,' '
           c1err = "Found 'nan' instead of dble. "
           c2err = "Use grep to find the nan in this version. "
           CALL MADABORT("RDNDBLE", c1err, c2err )
        endif
        read ( cwords(i), * ) DARRAY(i)
!           read ( cwords(i), "(D20.10)" ) DARRAY(i)
    ENDDO
    RETURN
  END SUBROUTINE RDNDBLE 

! ==============================================
    CHARACTER*80 FUNCTION RDSTRING( OPT, FILE )

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! Read single-word STRING from FILE
! OPT=0 => read from file
! OPT=1 => read from $SNDATA_ROOT/file

    INTEGER OPT
    CHARACTER FILE*(*)

! local var

    INTEGER LL
    CHARACTER LOCFILE*(MXCHAR_FILENAME), string*80

! -----------

    IF ( OPT .EQ. 0 ) THEN
       LOCFILE = FILE
    ELSE
       LL = INDEX(SNDATA_ROOT,' ') - 1
       LOCFILE = SNDATA_ROOT(1:LL) // '/' // FILE
    ENDIF

    OPEN (UNIT=LUNTMP, FILE = LOCFILE, STATUS='OLD', ERR=666 )
    READ(LUNTMP,*) string
    CLOSE ( UNIT = LUNTMP )

    RDSTRING = STRING
    RETURN

! --------------------
666   CONTINUE
    LL = INDEX(FILE,' ') - 1
    C1ERR = 'Cannot open file: ' // file(1:LL)
    CALL MADABORT("RDSTRING", C1ERR, "")

    RETURN
  END FUNCTION RDSTRING



! ================================================
    SUBROUTINE INIT_SURVEY(NSN_VERS)
! 
! Created Jun 2011
! [re-named from RDSURVEY, and read-part moved to RDSURVEY_TEXT]
! 
! Process a few things from the global survey variables
! read from RDSURVEY_TEXT[FITS]
! 
! May 26, 2012: remove BX check; moved to LANDOLT_PREP.
! 
! Jun 15, 2013: check FILTER_REPLACE option
! 
! Aug 08, 2013: read and fill SURVEY_FIELD(i) for each survey field;
!               needed to sort overlapping fields so that
!               E2+E1 gets written as E1+E2
! 
! Dec 16 2014: fix string logic checking SURVEY; see L1.
! 
! Dec 09 2020: fix index bug setting SURVEY_TMP
! -------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER NSN_VERS  ! (I) number of data versions

    INTEGER  & 
         iwd, NWD, IFILT_OBS, LL, L1, L2, LNAM,  i, IDTMP  & 
        ,IDFIELD, NFILT_NONSURVEY, NFILT_ADD, iadd, jtmp  & 
        ,LWD0, LWD1, LWD2, LWD3

    character  & 
         FNAM*12, cwd*60, cwd_next*60, cwd_next2*60, cwd_next3*60  & 
        ,cfilt1*2, SURVEY_TMP*(MXCHAR_SURVEY)  & 
        ,SURVEYFILE*(MXCHAR_FILENAME)  & 
        ,NAME_forC*(MXCHAR_FILENAME)  & 
        ,upper*60, ctel*60, cfield*60  & 
        ,SURVEY_FILTERS_ORIG*(MXFILT_ALL)

    LOGICAL LTMP, LKEY, LSRVY

! functions
    INTEGER FILTINDX
    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------------- BEGIN --------------

    FNAM            = 'INIT_SURVEY'
    IDSURVEY        = 0
    NFIELD_SURVEY   = 0
    NSURVEY_LIST    = 0

    LNAM  = INDEX ( SURVEY_NAME,    ' ' ) - 1
    L2    = 0

    IF ( FREEZE_SURVEY_FILTERS ) THEN
!        leave SURVEY_FILTERS and NFILTDEF_XXX alone
    ELSE
      L2  = INDEX ( SURVEY_FILTERS, ' ' ) - 1
      NFILTDEF_SURVEY = L2
      NFILTDEF_READ   = L2
    ENDIF

! ----------------------------------------
! Apr 2013: check user-option to add non-survey filters.
!           These extra filters must be defined in the kcor/calib file,
!           otherwise RDKCOR will abort.
! Jun 1 2013:
!  Only add NONSURVEY filter if it's NOT already defined. For example,
!  consider SDSS ugrizUGRIZ, but the simulations are done only with ugriz.
!  Thus we can define NONSURVEY_FILTERS ='UGRIZ' which are added to the
!  simulation list to avoid a KCOR-abort, but for data these NONSURVEY_FILTERS
!  are ignored since they are already defined.

    NFILT_NONSURVEY = INDEX(NONSURVEY_FILTERS,' ') - 1
    SURVEY_FILTERS_ORIG = SURVEY_FILTERS
    NFILT_ADD = 0
    NONSURVEY_FILTERS_ADD = ''

    DO 601 iadd = 1, NFILT_NONSURVEY
      cfilt1 = NONSURVEY_FILTERS(iadd:iadd)
      jtmp   = INDEX(SURVEY_FILTERS_ORIG,cfilt1(1:1))

      if ( jtmp .LE. 0 ) then
         NFILT_ADD = NFILT_ADD + 1
         NONSURVEY_FILTERS_ADD =  & 
           NONSURVEY_FILTERS_ADD(1:NFILT_ADD-1) // cfilt1(1:1)
      endif

601   CONTINUE

    IF ( NFILT_ADD > 0 ) THEN
         NONSURVEY_FILTERS_ADD =  & 
           NONSURVEY_FILTERS_ADD(1:NFILT_ADD) // '  '  ! add padding

       print*,'   Add NON-SURVEY FILTERS = ',  & 
                    NONSURVEY_FILTERS_ADD(1:NFILT_ADD)
       call flush(6)
       SURVEY_FILTERS  = SURVEY_FILTERS(1:NFILTDEF_SURVEY) //  & 
                           NONSURVEY_FILTERS_ADD(1:NFILT_ADD)
       NFILTDEF_SURVEY = NFILTDEF_SURVEY + NFILT_ADD
       L2 = NFILTDEF_SURVEY
    ENDIF

! ----------------------------------------
! create IFILTDEF_MAP_SURVEY


     do i = 1, NFILTDEF_SURVEY
        cfilt1    = SURVEY_FILTERS(i:i)
        IFILT_OBS = FILTINDX(cfilt1)

        IFILTDEF_MAP_SURVEY(i)    = IFILT_OBS
        IFILTDEF_INVMAP_SURVEY(ifilt_obs) = i

        if ( IFILT_OBS .LE. 0 ) then
          CALL PRINT_PREABORT_BANNER("INIT_SURVEY"//char(0),40)
          print*,'   i = ', i,' of NFILTDEF_SURVEY=', NFILTDEF_SURVEY
          print*,'   IFILT_OBS = ', IFILT_OBS
          print*,'   CFILT1    = ', cfilt1

          c1err = 'FILTER = ' // cfilt1 //  & 
               ' is invalid for SURVEY=' // SURVEY_NAME(1:LNAM)
          c2err = 'SURVEY_FILTERS = |' //  & 
                    SURVEY_FILTERS(1:NFILTDEF_SURVEY) // '|'
          CALL MADABORT(FNAM, c1err, c2err )
        endif

     enddo

! ------------
! Now open SURVEY.DEF file to get
!  IDSURVEY
!  IDFIELD

    LL = INDEX ( SNDATA_ROOT , ' ' ) - 1
    SURVEYFILE = SNDATA_ROOT(1:LL) // '/SURVEY.DEF'

    LL        = INDEX(SURVEYFILE,' ')-1
    NAME_forC = SURVEYFILE(1:LL)//char(0)
    NWD  = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE, NAME_forC,  & 
                 FNAM//char(0), LL, 12)

    DO 110 iwd = 1, NWD-3

       CALL get_PARSE_WORD_fortran(iwd+0, cwd,       LWD0 )
       CALL get_PARSE_WORD_fortran(iwd+1, cwd_next,  LWD1 )
       CALL get_PARSE_WORD_fortran(iwd+2, cwd_next2, LWD2 )
       CALL get_PARSE_WORD_fortran(iwd+3, cwd_next3, LWD3 )
       if ( cwd .EQ. ' ' ) goto 110

       if ( cwd(1:7) .EQ. 'SURVEY:' ) then
          L1 = index(cwd_next,' ') - 1
          SURVEY_TMP  = cwd_next(1:L1)
          read(cwd_next2,*) IDTMP

!    store all survey names in case there are subSurveys
!    such as LOWZ(CFA3)
          NSURVEY_LIST = NSURVEY_LIST + 1
          SURVEY_NAME_LIST(NSURVEY_LIST) = SURVEY_TMP
          IDSURVEY_LIST(NSURVEY_LIST)    = IDTMP

          if ( SURVEY_TMP(1:L1) .EQ. SURVEY_NAME(1:LNAM) ) then
!     IDSURVEY should never change, but IDSUBSURVEY can change
             IDSURVEY=IDTMP ; IDSUBSURVEY=IDTMP
          endif
       endif

! ------------
! check SN field

       LKEY   = cwd .EQ. 'FIELD:'
       LSRVY  = .FALSE.
       IF ( LKEY ) THEN
          read(cwd_next3,*) SURVEY_TMP
          LSRVY  = SURVEY_TMP .EQ. SURVEY_NAME ! field for this survey
       ENDIF

       if ( LKEY .and. LSRVY ) then

        call UPCASE(cwd_next, upper)
        cfield = upper  ! FIELD in SURVEY.DEF file

        read(cwd_next2,*) IDFIELD

        IF ( IDFIELD > 0 ) THEN
           NFIELD_SURVEY = NFIELD_SURVEY + 1
           SURVEY_FIELDNAME(NFIELD_SURVEY)= cfield(1:MXCHAR_FIELDNAME)
           SURVEY_IDFIELD(NFIELD_SURVEY)  = IDFIELD
        ENDIF

        IF ( NFIELD_SURVEY > MXIDFIELD ) THEN
           write(c1err,681) MXIDFIELD
 681         format('NSURVEY_FIELD exceeds bound of ', I4)
           c2err = 'Check MXIDFIELD and SURVEY.DEF file.'
           CALL MADABORT(FNAM, c1err, c2err )
        ENDIF

! now compare with user request, and convert everything to upper case

         DO i = 1, NIDFIELD_LIST
! convert user name to uppercase
           CALL UPCASE(SNFIELD_LIST(i), upper)
           SNFIELD_LIST(i) = upper   ! user-specified field
           LL = INDEX( cfield, ' ' ) - 1

           if ( cfield(1:LL) .EQ. SNFIELD_LIST(i)(1:LL) ) then
             SNFIELD_LIST(i) = cfield  ! truncate extra characters
             IDFIELD_LIST(i) = IDFIELD
             write(6,206) CFIELD(1:LL), IDFIELD
206            format(T5,'Found requested FIELD=',A,2x,  & 
                     'with IDFIELD=', I3 )
           endif
          END DO

       endif  ! end of FIELD: if-block

110   CONTINUE

! ------------------------
! print some info to stdout

    IF ( IDSURVEY .LE. 0 .and. SUBSURVEY_NAME .EQ. '' ) then
      c1err = 'Could not find IDSURVEY for SURVEY='  & 
              // SURVEY_NAME(1:LNAM)
      c2err = 'Check ' // SURVEYFILE(1:72)
      CALL MADABORT(FNAM, c1err, c2err )
    endif

    if ( DO_GETINFO ) return ! Apr 2021

    write(6,201) SURVEY_NAME(1:LNAM), IDSURVEY
201   format(T5,'Found SURVEY  = ',A, 3x, '(IDSURVEY=',I3,')' )

    write(6,202) SURVEY_FILTERS(1:L2), NFILTDEF_SURVEY
202   format(T5,'Found FILTERS = ',A, 3x, '(NFILTDEF=',I3,')' )

    write(6,204) NSN_VERS
204   format(T5,'Found ',I7,' SN Candidates to Process.')

    print*,' '
    call flush(6)

! --------------------------------------------------------

    IF ( .NOT. DOALL_SNFIELD ) THEN
      DO i    = 1, NIDFIELD_LIST
        cfield  = SNFIELD_LIST(i)
        idfield = IDFIELD_LIST(i)
        if ( idfield .LT. 0 ) then
           LL = index( cfield, ' ' ) - 1
           c1err = 'Unknown FIELD: ' // cfield(1:LL)
           c2err = 'Check namelist variable SNFIELD_LIST '
           CALL MADABORT(FNAM, c1err, c2err )
        endif
      ENDDO
    ENDIF


! ---------------------------------
! check for FILTER_REPLACE option

    CALL INIT_FILTER_REPLACE()

! Apr 2021: check option to use survey name for TEXTFILE_PREFIX
    if ( TEXTFILE_PREFIX(1:6) == 'SURVEY' ) then
        TEXTFILE_PREFIX = SURVEY_NAME
    endif

    RETURN
  END SUBROUTINE INIT_SURVEY

! ===================================
    SUBROUTINE INIT_FILTER_REPLACE()

! Created Jun 2013 by R.K.
! Initialize IFILTOBS_REPLACE array based on user input
! string FILTER_REPLACE read from &SNLCINP namelist.
! After this init, can use FILTINDX_REPLACE function.
! 
! Example:  FILTER_REPLACE = 'UGRIZ -> ugriz' will
! internally replace UGRIZ with ugriz so that we can
! analyze 'ugriz' (intead of ugrizUGRIZ) and get
! all of the epochs. This is intended for cases in
! which UGRIZ and ugriz bands are very similar.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER  & 
        IFILTOBS, NF1, NF2, NFARROW, IFILT, NWD, LF  & 
       ,IFILTOBS_ORIG

    CHARACTER  & 
         ARROW*60  & 
        ,FILTLIST1*(MXFILT_ALL)  & 
        ,FILTLIST2*(MXFILT_ALL)  & 
        ,NAME_forC*(MXFILT_ALL)  & 
        ,FNAM*20  & 
        ,CFILT_ORIG*2  & 
        ,CFILT_REPLACE*2

! function
    INTEGER FILTINDX
    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------- BEGIN ------------

    NFILT_REPLACE = 0
    DO IFILTOBS = 1, MXFILT_ALL
      IFILTOBS_REPLACE(IFILTOBS) = IFILTOBS  ! default is 1-to-1 map
    ENDDO

    FNAM = 'INIT_FILTER_REPLACE'

    IF ( FILTER_REPLACE .NE. ' ' ) THEN

      LF        = INDEX(FILTER_REPLACE,' ', BACK=.TRUE.) - 1
      NAME_forC = FILTER_REPLACE(1:LF) // char(0)
      NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,NAME_forC,  & 
                   FNAM//char(0), LF, 20)

      CALL get_PARSE_WORD_fortran(1, FILTLIST1, NF1)
      CALL get_PARSE_WORD_fortran(2, ARROW,     NFARROW)
      CALL get_PARSE_WORD_fortran(3, FILTLIST2, NF2)

! idiot checks

      if ( NF1 .NE. NF2 ) THEN
         c1err = 'Cannot use FILTER_REPLACE with ' //  & 
                    'different size lists.'
         c2err = 'N(' // FILTLIST1(1:NF1) // ')' //  & 
               ' != N(' // FILTLIST2(1:NF2) // ')'
         CALL MADABORT(FNAM, c1err, c2err )
      ENDIF
      IF ( ARROW(1:2) .NE. '->' ) THEN
         c1err = 'Second string of FILTER_REPLACE must by  ->'
         c2err = 'Check &SNLCINP namelist.'
         CALL MADABORT(FNAM, c1err, c2err )
      ENDIF

      NFILT_REPLACE = NF1
      DO IFILT = 1, NFILT_REPLACE
         CFILT_ORIG    = FILTLIST1(IFILT:IFILT)
         CFILT_REPLACE = FILTLIST2(IFILT:IFILT)

         IFILTOBS_ORIG                 = FILTINDX(CFILT_ORIG)
         IFILTOBS_REPLACE(IFILTOBS_ORIG) = FILTINDX(CFILT_REPLACE)

         write(6,20) CFILT_ORIG, CFILT_REPLACE, IFILTOBS_ORIG,  & 
                          IFILTOBS_REPLACE(IFILTOBS_ORIG)
 20        format(T10,'FILTER_REPLACE:  ', A, ' -> ', A, 3x,  & 
              '(',I3,' -> ', I3, ')' )
         CALL FLUSH(6)
      ENDDO
      print*,' '

!        print*,' xxx FILTLIST1 = ', FILTLIST1(1:12), NF1
!        print*,' xxx FILTLIST2 = ', FILTLIST2(1:12), NF2
!        print*,' xxx ARROW = ', ARROW
!        print*,' xxx REPLACE ', NFILT_REPLACE,'  filters '


    ENDIF

    RETURN
  END SUBROUTINE INIT_FILTER_REPLACE

! =============================
    SUBROUTINE SET_SURVEY(NAME, NFILTDEF, IFILTDEF, LAMSHIFT )
! 
! May 18, 2008
! 
! Set survey and IFILTDEF_[INV]MAP_SURVEY variables
! by passing them (i.e., from simulation) rather than
! by calling RDSURVEY.
! 
! Jun 20 2017: pass LAMSHIFT as argument.
! 
! ---------


    USE SNPAR
    USE CTRLCOM
    USE FILTCOM
    USE SNCUTS
    USE SNLCINP_NML

    IMPLICIT NONE

! subroutine args

    CHARACTER NAME*(*)         ! (I) name of survey

    INTEGER  & 
          NFILTDEF             &  ! (I) number of survey filters
         ,IFILTDEF(NFILTDEF)  ! (I) absolute filter indices

    REAL LAMSHIFT(NFILTDEF) ! (I) LAMSHIFT for each band

! local args

    INTEGER ifilt, ifilt_obs

! ---------------- BEGIN ------------

    SURVEY_NAME     = NAME
    NFILTDEF_SURVEY = NFILTDEF


    DO ifilt    = 1, NFILTDEF
      ifilt_obs = IFILTDEF(ifilt)
      IFILTDEF_MAP_SURVEY(ifilt)        = ifilt_obs
      IFILTDEF_INVMAP_SURVEY(ifilt_obs) = ifilt
      FILTER_LAMSHIFT_FILT(ifilt_obs)   = LAMSHIFT(IFILT)
    ENDDO


! Apr 2013: make sure to abort on bad kcor redshift.
    ABORT_ON_BADZ = .TRUE.

    RETURN
  END SUBROUTINE SET_SURVEY

! ================================================
    SUBROUTINE ABORT_ON_CUTWIN_TREST()

! Feb 2014
! Utility called by photo-z fitting program to abort
! if any of the SNANA-Trest cuts are set ... because
! for photo-z fit we cannot make an a-priori cut that
! depends on redsfhit. Try CUTWIN_TOBSMIN, CUTWIN_TOBSMAX,
! and other cuts that do NOT depend on the fitted photo-z.
! 
!  To override this abort, set
!  &SNLCINP
!    ABORT_ON_TRESTCUT = F
! -------------------------------------------


    USE SNPAR
    USE SNCUTS
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER NCUT, i, icut, ICUTLIST(20)
    INTEGER NCUT_USE, LEN
    LOGICAL LLO, LHI
    REAL    VLO, VHI

    CHARACTER c1err*76, C2err*76

! --------------- BEGIN --------------
    NCUT = 0

    NCUT = NCUT + 1 ;   ICUTLIST(NCUT) = CUTBIT_TREST
    NCUT = NCUT + 1 ;   ICUTLIST(NCUT) = CUTBIT_TRESTMIN
    NCUT = NCUT + 1 ;   ICUTLIST(NCUT) = CUTBIT_TRESTMAX
    NCUT = NCUT + 1 ;   ICUTLIST(NCUT) = CUTBIT_TGAPMAX
    NCUT = NCUT + 1 ;   ICUTLIST(NCUT) = CUTBIT_T0GAPMAX

    NCUT_USE = 0
    Do 100 i = 1, NCUT
       icut  = ICUTLIST(i)
       VLO   = CUTWIN_VAR(1,icut)
       VHI   = CUTWIN_VAR(2,icut)

       LLO   = ( VLO > -.99*CUTVAL_OPEN )
       LHI   = ( VHI < +.99*CUTVAL_OPEN  )
       if ( LLO .or. LHI ) THEN
          NCUT_USE = NCUT_USE + 1
          LEN = INDEX(CUTVAR_NAME(icut),' ') - 1
          write(6,60) CUTVAR_NAME(icut), VLO, VHI
 60         format(T5,'WARNING: should not use CUTWIN_', A15,  & 
                  2x, F10.3,' - ', F10.3 )
          call flush(6)
       endif
 100  CONTINUE

    IF ( NCUT_USE > 0 .and. ABORT_ON_TRESTCUT ) THEN
       write(c1err,160) NCUT_USE
 160     format(I3,' Trest &SNLCINP cuts above are forbidden.' )
       c2err = 'Can override abort with &SNLCINP ABORT_ON_TRESTCUT=F'
       CALL MADABORT('ABORT_ON_CUTWIN_TREST', c1err, c2err )
    ENDIF

    RETURN
  END SUBROUTINE ABORT_ON_CUTWIN_TREST

! ================================================
    SUBROUTINE READ_EPIGNORE_FILE ( IVERS  )
! 
! Created May 6, 2008 by R.Kessler
! Read list of epochs to ignore.
! 
! Use namelist variable EPOCH_IGNORE_FILE:
! if 'DEFAULT', then read [SNDATA_PREFIX].IGNORE.
! 
! Fill variables in common EPXCOM
! 
! File syntax is
! 
!  IGNORE:  <CID1>  <MJD1>  <FILTER1>
!  IGNORE:  <CID2>  <MJD2>  <FILTER2>
!  etc ...
! 
!  where <FILTER> is 'g' or 'r' or 'i' etc ...
! 
!  and/or
! 
! PHOTFLAG_MSKREJ:  <mask>  ! apply to PHOTFLAG
!    or
! PHOTFLAG_BITLIST_REJECT: 2,3,4  ! convert bits to mask
! 
! PRIORITY for setting PHOTFLAG_MSKREJ
!   + user input via &SNLCINP -> top priority
!   + [VERSION].IGNORE file   -> next priority
! 
! Strategy is that PHOTFLAG_MSKREJ is initialized to -1 before
! reading &SNLCINP ... which is before calling this function.
! If PHOTFLAG_MSKREJ >=0 here, then it picked up a value from &SNLCINP
! and it is NOT allowed to be modified here.
! If PHOTFLAG_MSKREJ < 0 here, then it can be updated by contents
!    of [VERSION].IGNORE file.
! 
! Aug 27 2020: for externally supplied .IGNORE  file, make DOCANA check.
! Feb 06 2021: IGNORE file is now optional ! No longer required !
! May 13 2021: check for PHOTFLAG_MSKREJ in IGNORE file
! 
! ----------------------------------

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER  & 
          IVERS        ! (I) use this version index

! local var

    INTEGER LL, LL1, LL2, LL3, LV, LF, BIT
    INTEGER iwd, NWD, N, LUN, NREAD_IGNORE
    LOGICAL LEXIST, FOUND_EPIGNORE, NOTSET
    REAL*8 MJD8
    CHARACTER  & 
         LOCAL_FILENAME*(MXCHAR_FILENAME)  & 
        ,NAME_forC*(MXCHAR_FILENAME)  & 
        ,VERSION*(MXCHAR_VERSION)  & 
        ,cwd*(MXCHAR_FILEWORD), cflt*4  & 
        ,cwd_tmp*(MXCHAR_FILEWORD)  & 
        ,ccid*(MXCHAR_CCID)  & 
        ,ccid_last*(MXCHAR_CCID)  & 
        ,STRING_BITLIST*60  & 
        ,FNAM*30

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------- BEGIN ----------

    FNAM = 'READ_EPIGNORE_FILE'

    LUN = LUNTMP

    VERSION  = VERSION_PHOTOMETRY(ivers)
    LV       = index(VERSION,' ') - 1

    FOUND_EPIGNORE = .FALSE.
    STRING_BITLIST = ''

    CALL PRBANNER ( "READ_EPIGNORE_FILE:" )

    if ( EPOCH_IGNORE_FILE(1:4) .EQ. 'NONE' ) GOTO 8000
    if ( EPOCH_IGNORE_FILE(1:4) .EQ. 'NULL' ) GOTO 8000

! check if file is default, or if provided by user.

    IF ( EPOCH_IGNORE_FILE(1:7) .EQ. 'DEFAULT' ) THEN
      LL = INDEX( SNDATA_PREFIX, ' ' ) - 1
      LOCAL_FILENAME = SNDATA_PREFIX(1:LL) // '.IGNORE'

!       bail if default file does not exist (Feb 2021)
      inquire (file=LOCAL_FILENAME, exist=LEXIST)
      if ( .NOT. LEXIST ) THEN
         print*,'    Default IGNORE file does not exist -> skip read.'
         call flush(6)
         GOTO 8000
      endif
    ELSE
      LOCAL_FILENAME = EPOCH_IGNORE_FILE
    ENDIF

! for multiple versions, don't bother re-reading  the same
! IGNORE file
    IF ( IVERS > 1 ) then
      LF  = INDEX(LOCAL_FILENAME,' ') - 1
      IF ( LOCAL_FILENAME(1:LF) .EQ.  & 
             EPOCH_IGNORE_LASTFILE(1:LF) ) THEN
           RETURN
      ENDIF
    ENDIF

! -----------
! Now start reading file

    LL = INDEX( LOCAL_FILENAME, ' ' ) - 1
    NAME_forC = LOCAL_FILENAME(1:LL)//char(0)

    NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE, NAME_forC,  & 
               FNAM//char(0), LL, 30)

    write(6,50) LOCAL_FILENAME(1:LL)
50    format(T5,'Read epochs to ignore from: ', /, T5, A)

    NREAD_IGNORE = 0
    CCID_LAST    = EPOCH_IGNORE_CCID(1)

    DO 200 iwd = 1, NWD

       CALL get_PARSE_WORD_fortran(iwd,cwd,LL)

       if ( cwd .EQ. 'IGNORE:' ) then

         NREAD_IGNORE = NREAD_IGNORE + 1
         CALL get_PARSE_WORD_fortran(iwd+1,CCID,LL1)
         if ( IVERS .GT. 1 .and. ccid .EQ. CCID_LAST ) then
           print*,'   IGNORE file already read -> SKIP re-reading. '
           return
         endif

         NEPOCH_IGNORE = NEPOCH_IGNORE + 1
         N = NEPOCH_IGNORE

         if ( N .GT. MXEPOCH_IGNORE ) then
            write(c1err,661) N
661           format('NEPOCH_IGNORE=',I4,' exceeds array bound.')
            c2err = 'Check file: ' // LOCAL_FILENAME(1:LL)
            CALL MADABORT("READ_EPIGNORE", c1err, c2err )
         endif

         CALL get_PARSE_WORD_fortran(iwd+2,CWD_TMP,LL2)
         read(cwd_tmp,*) MJD8

         CALL get_PARSE_WORD_fortran(iwd+3,CWD_TMP,LL3)
         cflt =  CWD_TMP(1:4)
         EPOCH_IGNORE_CCID(N)  = CCID
         EPOCH_IGNORE_FILT(N)  = cflt
         EPOCH_IGNORE_MJD(N)   = MJD8

         LL = INDEX(CCID,' ' ) - 1
         IF ( LL .LE. 0 ) LL = MXCHAR_CCID

!           write(6,55) N, CCID(1:LL), MJD8, cflt
! 5         format(T4,I3,': Add ignore-epoch to list for CID=',A,
!     &         3x,'MJD=',F9.3, 2x,'FILT=',A)

       endif

! May 2021: check for PHOTFLAGs to reject epochs
       NOTSET = (PHOTFLAG_MSKREJ(1) < 0)
       if ( cwd .EQ. 'PHOTFLAG_MSKREJ:' .and. NOTSET ) then
          CALL get_PARSE_WORD_fortran(iwd+1,cwd_tmp,LL1)
          read(cwd_tmp,*) PHOTFLAG_MSKREJ(1)
       endif

       if ( cwd .EQ. 'PHOTFLAG_BITLIST_REJECT:' ) then
          CALL get_PARSE_WORD_fortran(iwd+1,STRING_BITLIST,LL1)
       endif

200   ENDDO

    EPOCH_IGNORE_LASTFILE = LOCAL_FILENAME
    call FLUSH(6)

! -------------------------------------------
8000  CONTINUE

! parse BITLIST after reading IGNORE file to avoid parse_WORD conflicts
    NOTSET = (PHOTFLAG_MSKREJ(1) < 0)
    IF ( STRING_BITLIST .NE. ' ' .and. NOTSET ) THEN
       PHOTFLAG_MSKREJ(1) = 0
       LL = INDEX(STRING_BITLIST,' ') - 1
       NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,  & 
                  STRING_BITLIST(1:LL)//char(0), FNAM//char(0), LL, 30)
       do iwd = 1, NWD
          CALL get_PARSE_WORD_fortran(iwd,cwd_tmp,LL1)
          read(cwd_tmp,*) BIT
          PHOTFLAG_BITLIST_REJECT(iwd) = BIT
          PHOTFLAG_MSKREJ(1) = IBSET(PHOTFLAG_MSKREJ(1),BIT)
       enddo
    ENDIF

! if PHOTFLAG_MSKREJ is still negative, then it has not been set;
! in this case, set it to zero so that it doesn't reject anything.
    IF ( PHOTFLAG_MSKREJ(1) < 0 ) PHOTFLAG_MSKREJ(1) = 0

    if ( NEPOCH_IGNORE > 0 ) then
       FOUND_EPIGNORE = .TRUE.
       write(6,801) NEPOCH_IGNORE
 801     format(T5,'EPOCH_IGNORE: Found ', I5,' epochs to IGNORE')
    endif

    if ( PHOTFLAG_MSKREJ(1) > 0 ) then
       write(6,802) PHOTFLAG_MSKREJ(1)
       FOUND_EPIGNORE = .TRUE.
 802     format(T5,'EPOCH_IGNORE: ignore epochs with PHOTFLAG & ', I8 )
    endif

    if ( .NOT. FOUND_EPIGNORE ) then
       print*,'   EPOCH_IGNORE: ' //  & 
                'No epochs are ignored => use all epochs. '
    endif
    print*,' '

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE READ_EPIGNORE_FILE 

! =============================
    LOGICAL FUNCTION PASS_EPIGNORE_FILE(EP)

! Created Feb 2021
! Returns True if this epoch (EP) is NOT listed in the
! epoch IGNORE file associated with VERSION.IGNORE,
! or with user input EPOCH_IGNORE_FILE.
! Note that cuts are NOT applied here; this function simply
! checks against an already existing list of epochs to ignore.
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER EP ! (I) epoch index

    INTEGER  i, LL, ifilt_obs
    LOGICAL  LFILT, IGNORE, LMJD, LCID
    REAL*8    MJD8
    CHARACTER CFILT*2, CCID*(MXCHAR_CCID)

! --------------- BEGIN ----------

    IGNORE = .FALSE.

    DO 100 i = 1, NEPOCH_IGNORE

      CCID = SNLC_CCID
      LL   = ISNLC_LENCCID
      LCID = EPOCH_IGNORE_CCID(i)(1:LL) .EQ. CCID(1:LL)
      IF ( .NOT. LCID ) GOTO 100

! leave big margin for MJD check (.002 days) in case IGNORE
! file has round-off errors for the MJD.
      MJD8 = SNLC8_MJD(ep)
      LMJD = ( ABS(EPOCH_IGNORE_MJD(i)-MJD8) < 0.002 )
      IF ( .NOT. LMJD ) GOTO 100

      ifilt_obs = ISNLC_IFILT_OBS(ep)
      cfilt     = filtdef_string(ifilt_obs:ifilt_obs)
      LFILT     = EPOCH_IGNORE_FILT(i)(1:1) .EQ. CFILT(1:1)
      IF ( .NOT. LFILT ) GOTO 100

! if we get here, then give message to ignore epoch
! and skip down to over-write part.

      IGNORE = .TRUE.

      IF ( STDOUT_UPDATE ) THEN
        write(6,101) CCID, MJD8, cfilt
      ENDIF
101     format(T8,'IGNORE EPOCH FOR: CID=',A8,  & 
               3x,'MJD=',F9.3, 3x,'FILT=',A )

100   CONTINUE

    PASS_EPIGNORE_FILE = (.NOT. IGNORE)

    RETURN
  END FUNCTION PASS_EPIGNORE_FILE

! ================================================
    SUBROUTINE CHECK_FORMAT()

! Created Jun 15, 2011 by R.Kessler
! Check if format is FITS or TEXT;
! set global logical FORMAT_FITS or FORMAT_TEXT
! 
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE FITSCOM

    IMPLICIT NONE

    INTEGER    IVERS, OPT, ISTAT, LEN_VERS, LEN_PATH
    CHARACTER  VERSION*(MXCHAR_VERSION), PATH*(MXCHAR_PATH)

! ------------ BEGIN -----------

! check for subfolder; e.g., VERSION_PHOTOMETRY = 'JLA2014/JLA2014_CSP'
    CALL PREP_VERSION_SUBFOLDER()

! use first version on list since all versions must have same format
    IVERS    =  1
    VERSION  = VERSION_PHOTOMETRY(ivers)
    LEN_VERS = INDEX(VERSION,' ') - 1
    VERSION  = VERSION(1:LEN_VERS) // char(0)

    LEN_PATH = INDEX(PRIVATE_DATA_PATH,' ') - 1
    PATH     = PRIVATE_DATA_PATH(1:LEN_PATH) // char(0)

    OPT      = 1  ! option to check format and return
    ISTAT    = RD_SNFITSIO_PREP(OPT, PATH, VERSION,  & 
           LEN_PATH, LEN_VERS)

    IF ( ISTAT .GT. 0 ) THEN
       FORMAT_FITS = .TRUE.
    ELSE
       FORMAT_TEXT = .TRUE.
    ENDIF

    RETURN
  END SUBROUTINE CHECK_FORMAT


! ================================================
    SUBROUTINE SELECT_EPOCH_DRIVER(OPT)

! Created Feb 2021 (used by REFAC_READ_DATA)
! 
! OPT += 1 -> Apply all epoch cuts except for Trest
!                 (e..g, before calling SET_PEAKMJD
! 
! OPT += 2 -> Apply only Trest cut (e.g, after SET_PEAKMJD)
! 
! OPT  = 3 -> apply all obs/epoch cuts


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PKMJDCOM

    IMPLICIT NONE

    INTEGER OPT  ! (I)

    REAL*8    MJD8, MJD8_LAST
    INTEGER   NMJD_STORE_NEW, NMJD_STORE_ORIG, NEWMJD
    INTEGER   EPMIN, EPMAX, ep_orig, ep_new, IMJD
    INTEGER   IFILT, IFILT_OBS
    REAL      PEAKMJD, MJD, TOBS, TREST, z1, z
    LOGICAL   APPLY_CUTS, APPLY_TREST
    LOGICAL   PASS_CUTS, PASS_TREST, PASS
    CHARACTER FNAM*20

! function
    LOGICAL PASS_EPOCH_CUTS, PASS_EPOCH_TREST

! ------------------- BEGIN -------------------

    FNAM = 'SELECT_EPOCH_DRIVER'
    NMJD_STORE_NEW  = 0
    NMJD_STORE_ORIG = ISNLC_NEWMJD_STORE
    EP_NEW          = 0
    MJD8_LAST       = -9.0

    APPLY_CUTS  = BTEST(OPT,0)
    APPLY_TREST = BTEST(OPT,1)

    DO 101 NEWMJD = 1, ISNLC_NEWMJD_STORE
       EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
       EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

       PASS_CUTS  = .TRUE.
       PASS_TREST = .TRUE.

       if ( APPLY_CUTS ) then
          PASS_CUTS = PASS_EPOCH_CUTS(1,EPMIN)
       endif

       if ( APPLY_TREST ) then
          PEAKMJD = SNLC_SEARCH_PEAKMJD
          z       = SNLC_REDSHIFT
          PASS_TREST = PASS_EPOCH_TREST(z,PEAKMJD,EPMIN)
       endif

       PASS = PASS_CUTS .and. PASS_TREST

       if ( .NOT. PASS ) goto 101

       NMJD_STORE_NEW = NMJD_STORE_NEW + 1
       IMJD           = NMJD_STORE_NEW
       ISNLC_NFILT_NEWMJD(IMJD) = ISNLC_NFILT_NEWMJD(NEWMJD)

! overwrite EPOCH-dependent arrays to exclude rejected epochs

      DO 102 EP_ORIG = EPMIN, EPMAX

       EP_NEW = EP_NEW + 1

       MJD8  = SNLC8_MJD(EP_ORIG)
       if ( MJD8 > MJD8_LAST ) then
          if ( EP_ORIG == EPMIN ) then
            ISNLC_EPOCH_RANGE_NEWMJD(1,IMJD) = EP_NEW
          endif
            ISNLC_EPOCH_RANGE_NEWMJD(2,IMJD) = EP_NEW
       endif
       MJD8_LAST = MJD8

       if ( EP_NEW > EP_ORIG ) then
          CALL PRINT_PREABORT_BANNER(FNAM//char(0),40)
          print*,'    OPT = ', OPT
          print*,'    EPMIN, EPMAX = ', EPMIN, EPMAX
          print*,'    NEWMJD(ORIG,NEW) = ', NEWMJD, IMJD
          write(c1err,61) EP_NEW, EP_ORIG
61          format('Invalid EP_NEW=',I4, ' > EP_ORIG=',I4)
          c2err = 'Something is messed up.'
          CALL MADABORT(FNAM, c1err, c2err)
       endif

       CALL MOVE_SNLC_ARRAYS(EP_NEW,EP_ORIG)

102     CONTINUE
101     CONTINUE

    ISNLC_NEWMJD_STORE = NMJD_STORE_NEW

    RETURN
  END SUBROUTINE SELECT_EPOCH_DRIVER


! ===========================================
    LOGICAL FUNCTION PASS_EPOCH_CUTS(OPTMASK,EP)

! Created Feb 2021
! Return TRUE if epoch passes explicit epoch cuts such as
! PSF, PHOTFLAG, etc ..
!   (no implicit cuts such as TREST)
! 
! Mar 25 2021: init BADPHOT=F (bug fix)
! -----------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PKMJDCOM

    IMPLICIT NONE

! input args
    INTEGER OPTMASK, EP  ! (I)

! local args
    INTEGER PHOTFLAG, OV, i, IFILT, IFILT_OBS, IDFIELD
    LOGICAL BADPHOT, REJECT, LXMJD, LTMP, NO_DETECT
    REAL    ZP, PSF, MJD, SNR
! functions
    LOGICAL PASS_EPIGNORE_FILE, SELECT_EARLYLC

! ----------- BEGIN ------------

    PASS_EPOCH_CUTS = .TRUE.  ! init return arg
    REJECT          = .FALSE.

    IFILT_OBS = ISNLC_IFILT_OBS(EP)
    IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
    IDFIELD   = ISNLC_IDFIELD(EP)
    PHOTFLAG  = ISNLC_PHOTFLAG(ep)
    ZP        = SNLC_ZEROPT_forCUT(ep)
    PSF       = SNLC_PSF_FWHM_ARCSEC(EP)
    MJD       = SNGL(SNLC8_MJD(ep))

! Jul 2020: check option to write list of rejected epochs
#if defined(SNANA)
    CALL WRITE_EPOCH_IGNORE_EXEC(ep)
#endif

    if ( .NOT. PASS_EPIGNORE_FILE(EP) ) then
       PASS_EPOCH_CUTS = .FALSE.
       return
    endif

! ----------------------------------------------
! now check PHOTFLAG mask.

    BADPHOT = .FALSE.

! check OR-logic
! ==> reject epoch if any PHOTFLAG_MSKREJ bit is set
    OV   = IAND ( PHOTFLAG_MSKREJ(1), PHOTFLAG )
    IF ( OV > 0 ) THEN
       BADPHOT = .TRUE.  ;  REJECT  = .TRUE.
    ENDIF

! Now check AND logic
! ==> reject epoch if all PHOTFLAG_MSKREJ(2:5) bits are set

    do i = 2, 5
      OV   = IAND ( PHOTFLAG_MSKREJ(i), PHOTFLAG )
      IF ( OV > 0 .and. OV == PHOTFLAG_MSKREJ(i) ) THEN
           BADPHOT = .TRUE. ; REJECT  = .TRUE.
      ENDIF
    enddo

! ----------------------
! Check EARLY-epoch select logic (Mar 2015)

    IF ( .not. SELECT_EARLYLC(ep) ) REJECT = .TRUE.

! -------------------------------------
! July 11 2020
!  Check cuts on PSF and ZP ... doesn't change analysis,
!  but allows new OUT_EPOCH_IGNORE_FILE option to work properly.

! xxx      IF( RESTORE_DES3YR ) GOTO 600 ! skip ZP and PSF cut for DES3YR

    IF ( ZP < CUTWIN_ZP(1)   ) REJECT = .TRUE.
    IF ( ZP > CUTWIN_ZP(2)   ) REJECT = .TRUE.

    IF ( PSF < CUTWIN_PSF(1) ) REJECT = .TRUE.
    IF ( PSF > CUTWIN_PSF(2) ) REJECT = .TRUE.

! MJD cut
    LXMJD = ( MJD > CUTWIN_MJD_EXCLUDE(1) .and.  & 
                MJD < CUTWIN_MJD_EXCLUDE(2) )
    IF ( LXMJD ) REJECT = .TRUE.

    IF ( MJD < CUTWIN_MJD(1) ) REJECT = .TRUE.
    IF ( MJD > CUTWIN_MJD(2) ) REJECT = .TRUE.

! - - - - - - - - - - - -
! reject epochs for filter band that has SNRMIN > 99
    if ( IFILT < 900 ) then
      if ( CUTWIN_SNRMIN_FILT(1,IFILT) > 99.9999 ) REJECT = .TRUE.
    endif

! check selected field.
    LTMP = DOALL_SNFIELD .or. (IDFIELD > 0)
    IF ( .NOT. LTMP ) REJECT = .TRUE.

! Sep 2020: apply SNR cut for non-detections - to reject crazy fluxes.
    IF ( PHOTFLAG_DETECT > 0 ) THEN
       NO_DETECT = ( IAND(PHOTFLAG,PHOTFLAG_DETECT) == 0 )
       SNR       = SNLC_FLUXCAL(EP) / SNLC_FLUXCAL_ERRTOT(EP)
       if ( NO_DETECT ) then
          if ( SNR < CUTWIN_SNR_NODETECT(1)  ) REJECT = .TRUE.
          if ( SNR > CUTWIN_SNR_NODETECT(2)  ) REJECT = .TRUE.
       endif
    ENDIF

! ----------------------------
 600  CONTINUE

! set return-function value
    PASS_EPOCH_CUTS = ( .not. REJECT)

! check option to count BAD PHOTOMETRY epochs based on PHOTFLAG
    IF ( BTEST(OPTMASK,0) ) THEN
       NEPOCH_TOT = NEPOCH_TOT + 1
       IF ( .NOT. REJECT)  NEPOCH_USE = NEPOCH_USE + 1
       IF ( BADPHOT ) THEN
          NEPOCH_BADPHOT     = NEPOCH_BADPHOT     + 1
          NEPOCH_BADPHOT_SUM = NEPOCH_BADPHOT_SUM + 1
       ENDIF
    ENDIF


    RETURN
  END FUNCTION PASS_EPOCH_CUTS

! = = = = = = = = = = = = = = = = =
    LOGICAL FUNCTION PASS_EPOCH_TREST(z,PEAKMJD,EP)

! Created Feb 2021
! Return TRUE if epoch passes TREST cut


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PKMJDCOM

    IMPLICIT NONE

! input args
    REAL z, PEAKMJD  ! (I) redshift and peak MJD
    INTEGER EP       ! (I) epoch index to check Trest cut

! local args

    REAL MJD, ZZ, Trest
    LOGICAL VALID_PEAKMJD

! ----------- BEGIN ------------

    PASS_EPOCH_TREST = .TRUE.

    VALID_PEAKMJD = ( PEAKMJD > 40000.  )  ! valid PEAKMJD

! evaluate epoch cut in rest-frame
    IF ( VALID_PEAKMJD ) THEN
       ZZ = 1.0                       ! if no valid z
       IF ( z > 0.0 ) ZZ  = 1.0 + z   ! valid z
       MJD   = SNGL( SNLC8_MJD(EP) )
       Trest = (MJD - PEAKMJD) / ZZ

       PASS_EPOCH_TREST =  & 
              ( Trest .GE. cutwin_Trest(1) .and.  & 
                TRest .LE. cutwin_Trest(2) )
    ENDIF

    RETURN
  END FUNCTION PASS_EPOCH_TREST


! =============================
    LOGICAL FUNCTION PASS_SNTYPE(itype)
! 
! Created Feb 22 2021 [modified from SNTYPESTAT]
! 
! Returns True if ITYPE passes cuts; returns false otherwise.
! Also increments global NACCEPT_TYPE if ITYPE is valid.
! 
! Jul 12 2021: abort if ITYPE < 0
! Mar 08 2022: apply CUTWIN_SNTYPE
! -------------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! function args

    INTEGER  ITYPE  ! SN-integer type from data header

! local args

    CHARACTER FNAM*12
    INTEGER i
    REAL RTYPE

! -------------- BEGIN --------------

    FNAM = 'PASS_SNTYPE'
    PASS_SNTYPE = .FALSE.

    IF ( ITYPE < 0 ) THEN
       write(c1err,61) ITYPE, SNLC_CCID
       write(c2err,62) MXTYPE
61       format('Invalid data header TYPE = ', I6, ' for CID=',A )
62       format('Valid TYPE range is 0 to ', I5)
       CALL MADABORT(FNAM, c1err, c2err)
    ENDIF

    do i = 1, NSNTYPE_IGNORE
      if ( SNTYPE_IGNORE(i) .EQ. ITYPE ) then
         PASS_SNTYPE = .FALSE.
         RETURN
      endif
    enddo

    RTYPE = float(ITYPE)
    IF ( RTYPE < CUTWIN_SNTYPE(1) ) RETURN
    IF ( RTYPE > CUTWIN_SNTYPE(2) ) RETURN

    IF ( NSNTYPE_LIST .LE. 0 ) THEN
       ! don't check  anything
       PASS_SNTYPE = .TRUE.
       RETURN
    ENDIF


    do i = 1, NSNTYPE_LIST
      if ( SNTYPE_LIST(i) .EQ. ITYPE ) PASS_SNTYPE = .TRUE.
    enddo

    IF ( PASS_SNTYPE ) THEN
       NACCEPT_TYPE = NACCEPT_TYPE + 1
    ELSE
       IF ( LDMP_SNFAIL ) then
            print*,'  ** WARNING ** CID= ',snlc_ccid,  & 
             ' (TYPE=',itype,')  FAILS TYPE CUT '
       endif
    endif

    RETURN
  END FUNCTION PASS_SNTYPE


! ======================================
    SUBROUTINE UNPACK_SNFITSIO_STR(NVAL,KEYNAME,STRING)

! 
!  Created Jun 2011
!  Following a call to RD_SNFITSIO_STR(IROW, KEYNAME, STRING ..)
!  this routine unpacks and decodes the retured STRING and fills
!  the following arrays,
! 
!  KEYNAME    filled array
!   FLT       ISNLC_IFILT_OBS
!   FIELD     ISNLC_IDFIELD & SNLC_FIELD
! 
! Feb 11 2021: fix computation of LENS (length of STRING)
! -----------------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! subtroutine args
    INTEGER   NVAL    ! (I) expect number of values to extract from STRING

    CHARACTER  & 
         KEYNAME*(*)   &  ! (I) table par name read from fits file
        ,STRING*(*)   ! (I) string returned from RD_SNFITSIO_STR

! local args

    INTEGER iwd, LWD, NWD, LK, LENS, MSKOPT
    CHARACTER  CWD*80, FNAM*30
    LOGICAL LVAR_FLT, LVAR_FIELD

! functions
    INTEGER  FILTINDX_REPLACE, GET_IDFIELD

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ----------------- BEGIN --------------

    IF ( NVAL .LT. 0 ) RETURN  ! all epochs fail cut (May 16, 2012)

    FNAM = 'UNPACK_SNFITSIO_STR'

! determine variable type from KEYNAME

    LVAR_FLT   = .FALSE.
    LVAR_FIELD = .FALSE.

    IF ( KEYNAME(1:3) .EQ. 'FLT' ) THEN
      LVAR_FLT = .TRUE.
    ELSE IF ( KEYNAME(1:5) .EQ. 'FIELD' ) THEN
      LVAR_FIELD = .TRUE.
    ELSE
      c1err = 'Unrecognized KEYNAME = ' // KEYNAME
      c2err = 'Check FITS table.'
      CALL MADABORT('UNPACK_SNFITSIO_STR', c1err,c2err)
    ENDIF

! break up STRING into substrings
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    LENS   = INDEX(STRING//' ',' ' ) - 1  ! Feb 2021

!       print*,' xxx UNPACK: ', LENS, '|', STRING, '|'
    NWD    = STORE_PARSE_WORDS(MSKOPT, STRING(1:LENS)//char(0),  & 
                  FNAM//char(0), LENS, 30 )

    IF  ( NWD .NE. NVAL ) THEN
      LK = INDEX(KEYNAME//' ' ,' ') -1
      CALL PRINT_PREABORT_BANNER("UNPACK_SNFITSIO_STR"//char(0),40)
      print*,'   NWD, NVAL = ', NWD, NVAL
      print*,'   SNLC_CCID = ', SNLC_CCID
      print*,'   KEYNAME   = ', KEYNAME,'      LENKEY=', LK
	print*,'   ISVAR_[FLT,FIELD] = ',  & 
               LVAR_FLT, LVAR_FIELD
      print*,'   LEN(STRING) = ', LENS
      print*,'   STRING = ', STRING(1:1000)
      call flush(6)

      write(c1err,61) NWD, NVAL
      write(c2err,62)  & 
             SNLC_CCID(1:ISNLC_LENCCID)  & 
            ,KEYNAME(1:LK)

61      format('NWD=',I5,' substrings, but expected NVAL=',I6)
62      format('CID=',A, 2x, 'KEYNAME=',A )
      CALL MADABORT('UNPACK_SNFITSIO_STR', c1err,c2err)

    ENDIF

    DO 100 iwd = 1, NWD

       CALL get_PARSE_WORD_fortran(iwd,CWD,LWD)

       if ( LVAR_FLT ) then
          ISNLC_IFILT_OBS(iwd) = FILTINDX_REPLACE(CWD)
       else if ( LVAR_FIELD ) then
          SNLC_FIELD(iwd)    = CWD(1:LWD)
          ISNLC_IDFIELD(iwd) = GET_IDFIELD(CWD)
       endif

100   CONTINUE

    RETURN
  END SUBROUTINE UNPACK_SNFITSIO_STR
! ========================================
    SUBROUTINE CHECK_FITSABORT(IVERS,NRD,KEY)

! Created Jun 16, 2011
! ABORT if NRD = 0 or of NRD exceeds bound of MXEPOCH
! Call this routine after each fits-read for a REQIURED key.
! 
! Note that NRD = -9 is OK because this means that
! all epochs failed the epoch-MASK  cuts.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE


    INTEGER  & 
        IVERS    &  ! (I) version index
       ,NRD     ! (I) Number of elements returned for fits key

    CHARACTER KEY*(*)  ! (I) name of fits key

! local var

    CHARACTER VERSION*100
    INTEGER LK, LV

! ---------- BEGIN ----------


    IF ( NRD .EQ. 0 .or. NRD .GT. MXEPOCH ) THEN

       print*,' ERROR: NRD = ', NRD

       LK = INDEX(KEY,' ') - 1
       c1err = 'Non-existent FITS parameter: ' // KEY(1:LK)  & 
                 // ' for CID='// SNLC_CCID(1:ISNLC_LENCCID)

       VERSION  = VERSION_PHOTOMETRY(ivers)
       LV = INDEX(VERSION,' ') - 1
       c2err = 'Check fits files for version ' // VERSION(1:LV)

       CALL MADABORT("CHECK_FITSABORT",c1err, c2err)
    ENDIF

    RETURN
  END SUBROUTINE CHECK_FITSABORT



! ==================================
    SUBROUTINE PRINT_RDSN()

! Called before SNANA_DRIVER, print one-line summary of SN that has
! been read/parsed.
! 
! Apr 19 2022: replace REJECT_PRESCALE with MOD function to have
!              more predictable prints.
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER N_UPDATE_PS1, N_UPDATE_PS10, PS_UPDATE

! ------------- BEGIN -----------

! determine STDOUT_UPDATE logical to dump info

    N_UPDATE_PS1  = 50 ! update first 50 events
    N_UPDATE_PS10 = 300
    if ( REDUCE_STDOUT_BATCH ) then
       N_UPDATE_PS1  = 5
       N_UPDATE_PS10 = 100
    endif

    IF ( N_SNLC_PROC < N_UPDATE_PS1 ) THEN
       PS_UPDATE = 1
    ELSE IF ( N_SNLC_PROC < N_UPDATE_PS10 ) THEN
       PS_UPDATE = 10
    ELSE
       PS_UPDATE = 100
    ENDIF

    IF ( FORCE_STDOUT_BATCH ) THEN
       STDOUT_UPDATE = .TRUE.  ! for interactive debug only
    ELSE
       STDOUT_UPDATE = (MOD(N_SNLC_PROC,PS_UPDATE) .EQ. 0)
    ENDIF


    if ( .NOT. STDOUT_UPDATE ) RETURN

    print*,' *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-' //  & 
             '*-*-*-*-*-*-*-*-*-*-*-*-*-*-'

    IF ( SNLC_CASTCID .EQ. 'CHAR' ) then
        write(6,30) SNLC_CID, SNLC_CCID(1:ISNLC_LENCCID)
30        format(T4,'Assign integer CID=',I6,' for ', A)
    endif

    write(6,200)  & 
            SNLC_CCID  & 
          , ISNLC_NEWMJD_FOUND  & 
          , ISNLC_NEWMJD_STORE

200   format(T5,' Done Reading CID ', A16,  & 
             2x, 'NMJD(found,stored)=', I4,1x,I4 )

    CALL FLUSH(6)
    RETURN
  END SUBROUTINE PRINT_RDSN

! =================================
    SUBROUTINE SET_HOSTGAL_PREFIX(IGAL,PREFIX,LENPRE)

! Feb 2019
! For inoput IGAL, return PREFIX = 'HOSTGAL' or 'HOSTGAL2', etc ..,
! along with LENPRE = length of PREFIX string.
! The output is used for parsing and writing output.

    INTEGER IGAL         ! (I) 1,2 ... host-match index
    CHARACTER PREFIX*20  ! (O) key prefix
    INTEGER LENPRE       ! (O) length of prefix string

    IF ( IGAL == 1 ) THEN
       PREFIX = 'HOSTGAL'; LENPRE=7
    ELSE
       write(PREFIX,20) IGAL
 20      format('HOSTGAL',I1.1)
       LENPRE = 8
    ENDIF
    RETURN
  END SUBROUTINE SET_HOSTGAL_PREFIX

! ==================================
    SUBROUTINE SET_EPVAR_MISC(ep)
! Created Jul 11 2020:
! Compute a few misc variables from those read from data files
!   + SNLC_ZEROPT_forCUT(ep) based on cut in ADU or Npe
!   + SNLC_PSF_FWHM_ARCSEC(ep)


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER ep      ! (I) epoch index

    REAL   ZPADD_pe, GAIN
    REAL*8 PSFSIG1, PSFSIG2, PSFRAT, NEA

    EXTERNAL NoiseEquivAperture
    REAL*8   NoiseEquivAperture

! -------------- BEGIN ------------

    SNLC_ZEROPT_forCUT(ep) = SNLC_ZEROPT(ep) ! default is ADU

    if ( cutwin_zpNPE(1) > 0. ) then
       GAIN     = SNLC_GAIN(ep)

       if ( GAIN < 0.01 ) then
          GAIN = 0.001  ! -> results in very low ZP that should fail cut
       endif

       ZPADD_pe = 2.5*log10(GAIN)
       SNLC_ZEROPT_forCUT(ep) = SNLC_ZEROPT(ep) + ZPADD_pe
    endif

! if PSF unit is NEA, convert back to effected PSF-sigma;
! otherwise compute and store NEA
    IF ( UNIT_PSF_NEA ) then
       NEA     = DBLE(SNLC_PSF_NEA(ep))
       PSFSIG1 = SQRT(NEA/(4.0*PI))
       SNLC_PSF_SIG1(ep)  = SNGL(PSFSIG1)
       SNLC_PSF_SIG2(ep)  = 0.0
       SNLC_PSF_RATIO(ep) = 0.0
    else
       PSFSIG1 = DBLE(SNLC_PSF_SIG1(ep))
       PSFSIG2 = DBLE(SNLC_PSF_SIG2(ep))
       PSFRAT  = DBLE(SNLC_PSF_RATIO(ep))
       NEA     = NoiseEquivAperture(PSFSIG1,PSFSIG2,PSFRAT)
       SNLC_PSF_NEA(ep) = SNGL(NEA)
    endif

    SNLC_PSF_FWHM_ARCSEC(ep) =  & 
         SNLC_PSF_SIG1(ep) * (SNLC_PIXSIZE * 2.355)

    RETURN
  END SUBROUTINE SET_EPVAR_MISC

! =================================
    SUBROUTINE CHECKSTRING_CCID(CCID)

!  Jun 2016: abort if illegal char in SNLC_CCID


    CHARACTER CCID*(*)

! local

    INTEGER, PARAMETER  :: NBADCHAR = 14
    INTEGER ichar
    CHARACTER  BADCHAR_LIST*(NBADCHAR)
    CHARACTER  c1err*44, c2err*44, c1*2
! ----------- BEGIN ---------------

    BADCHAR_LIST = ':;!@#$%^&*()=?'

    DO 100 ichar = 1, NBADCHAR
       C1 = BADCHAR_LIST(ichar:ichar)
       IF ( INDEX(CCID,C1) > 0 ) THEN
          c1err = 'Invalid CCID name = ' // CCID
          c2err = 'because it contains invalid char = ' // c1
          CALL MADABORT('CHECKSTRING_CCID', c1err, c2err)
       ENDIF
100   CONTINUE

    RETURN
  END SUBROUTINE CHECKSTRING_CCID

! ==============================
    SUBROUTINE CHECK_FILTERS_KEY(FILTERS)

! Apr 2013  [code moved from PARSE_HEAD]
! Called from PARSE_HEAD for text files to check that
! FILTERS argument is always the same. Abort if the
! FILTERS argument ever changes.
! 
! -------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    CHARACTER FILTERS*(*)  ! (I) argument of "FILTERS:" key

    INTEGER   NFILT, IFILT, itmp, IFILT_OBS
    CHARACTER cfilt1*2

! functions
    INTEGER FILTINDX

! ---------------- BEGIN ----------------

    NFILT     = INDEX ( FILTERS//' ', ' ' ) - 1

    DO itmp       = 1, NFILT
       cfilt1     = FILTERS(itmp:itmp)
       ifilt_obs  = FILTINDX(cfilt1)
       IFILT      = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

! check if sparse index is valid.
! Apr 2013: better error messaging.

       if ( ifilt .LE. 0 ) then
          CALL PRINT_PREABORT_BANNER("PARSE_HEAD"//char(0),40)
          print*,' FILTERS  defined in 1st data file:'
          print*,'      ', SURVEY_FILTERS(1:NFILTDEF_SURVEY)
          print*,' FILTERS  defined for '  & 
                 // SNLC_CCID(1:ISNLC_LENCCID), ' : '
          print*,'      ', FILTERS(1:NFILT)
          print*,'  --> FILTERS argument must be the same ',  & 
                      'in each data file. '

          c1err = 'Filter ' // cfilt1 // ' is unrecognized '  & 
                 // ' for CID=' // SNLC_CCID(1:ISNLC_LENCCID) // ' ; '

          c2err = 'see PRE-ABORT comments above.'
          CALL MADABORT("PARSE_HEAD", c1err, c2err)
       endif

    ENDDO                     ! end of ITMP loop

    RETURN
  END SUBROUTINE CHECK_FILTERS_KEY


! ===================================
    SUBROUTINE PARSE_PARENTHESES(INKEY,OUTSTR)

! Parse inKey to extract varName from ().
! Example:
!   Input : INKEY   = SIMSED_PARAM(COSANGLE)
!  -->
!   Ouptut: OUTSTR = COSANGLE
! 
! ---------------------

! subroutine args
    CHARACTER INKEY*(*)   ! (I) input key with varname in ()
    CHARACTER OUTSTR*(*)  ! (O) string in ()

! local args
    INTEGER J1, J2, LK
    CHARACTER C1ERR*72, C2ERR*72
! -------------- BEGIN ---------------

    J1 = INDEX( INKEY , '(' )
    J2 = INDEX( INKEY , ')' )

    IF ( J1 .EQ. 0 .or. J2 .EQ. 0 ) THEN
       LK = INDEX( INKEY, ' ' ) - 1
       c1err = 'Invalid INKEY = ' // INKEY(1:LK)
       c2err = 'Expecting to find key-string  SIMSED_XXX(PARNAME)'
       CALL MADABORT("PARSE_PARENTHESES", c1err,c2err)
    ENDIF

    OUTSTR = INKEY(j1+1:j2-1)

    RETURN
  END SUBROUTINE PARSE_PARENTHESES


! ====================================
    SUBROUTINE PARSE_CID ( ccid, iauc, cid, USECID )
! 
! Convert character string CCID into integer CID.
! Returns logical USECID=T if this CID should be processed.
! 
! 
! Jul 28 2014; allow ABORT_ON_DUPLCID=F to count duplicates without aborting
! 
! Jul 31 2015: move GET_INTERP_MJDLIST call before LCIDSELECT call
! 
! Dec 2 2015: pass IAUC arg.
! 
! ---------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE INTERPCM

    IMPLICIT NONE

    CHARACTER  CCID*(*)  ! (I) character string for CID
    CHARACTER  IAUC*(*)  ! (I) char string for IAUC name
    INTEGER    CID       ! (O) integer CID
    LOGICAL    USECID    ! (O) T => process this CID

! local var

    INTEGER LCHAR, NMJD, ISTAT
    LOGICAL  USECID_LOCAL, ADDFLAG, LDUPL

    REAL*8  MJD8LIST(MXINTERP)

! function
    INTEGER CIDASSIGN, EXEC_CIDMASK
    LOGICAL LCIDSELECT

! --------------- BEGIN --------------

    CID = 0
    USECID = .FALSE.

! check if this CID should be processed.
    LCHAR = INDEX ( CCID, ' ' ) - 1
    LCHAR = MIN(MXCHAR_CCID,LCHAR) ! check limit on Number of chars

! convert char-string to integer CID
    CID = CIDASSIGN(CCID)

    IF ( CID .LE. 0 ) THEN
      c1err = 'Cannot determine integer CID for: ' // CCID
      write(c2err,602)  LCHAR
602     format('LEN(CCID) = ', I5 )
      CALL MADABORT("PARSE_CID", c1err, c2err)
    ENDIF

    IF ( CID > MXCID ) THEN
      write(c1err,611) CID, MXCID
      C2err = '   '
611     format('CID=',I9,' exceeds MXCID=',I9 )
      CALL MADABORT("PARSE_CID", c1err, c2err)
    ENDIF


! check for CIDs to interpolate (May 18, 2012)
    CALL GET_INTERP_MJDLIST(CCID, NMJD, MJD8LIST, ADDFLAG)
    IF ( ADDFLAG ) THEN
        NCCID_LIST = NCCID_LIST + 1
        SNCCID_LIST(NCCID_LIST) = CCID  ! Dec 4 2012
        USECID_LOCAL = .TRUE.
        GOTO 800
    ENDIF

! ----------------------------------------------------
    IF ( .NOT. LCIDSELECT(cid,ccid,iauc) ) RETURN
    SNLC_CCID      = CCID    ! in case CCID -> IAUC (Dec 3, 2015)
    ISNLC_LENCCID  = INDEX(CCID,' ') - 1  ! Mar 20 2016
    USECID_LOCAL   = .TRUE.
    GOTO 800

! ---------------------------
800   CONTINUE

    IF ( CID < MXCID_CHECK ) THEN
      ISTAT = EXEC_CIDMASK(2,CID)  ! returns 1 if CID bit is set
      LDUPL = (ISTAT > 0 )
    ELSE
      LDUPL = .FALSE. ! cannot check duplicates for very large CID
    ENDIF

    IF ( USECID_LOCAL ) THEN

! if this CID has already been used (found), abort.

      if ( LDUPL  ) THEN

        IF ( ABORT_ON_DUPLCID ) then
          write(c1err,801) CID, CCID(1:LCHAR)
801         format('CID= ',I9, ' (',A, ') appears more than once; ' )
          c2err = 'each CID can appear only once.'
          CALL MADABORT("PARSE_CID", c1err, c2err)
        ELSE
           print*,' *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-' //  & 
                  '*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
           write(6,866) CCID
 866         format(T10,'!!! DUPLICATE WARNING: reject CID=',A)
           call flush(6)
           N_DUPLICATE_CID = N_DUPLICATE_CID + 1
           return  ! return with USECID=F
        endif
      endif

      IF ( CID < MXCID_CHECK ) THEN
         ISTAT = EXEC_CIDMASK(1,CID)  ! set CID bit
      ENDIF

      NACCEPT_CID = NACCEPT_CID + 1
      USECID      = .TRUE.

    ELSE

      if ( LDMP_SNFAIL ) then
         print*,'  ** WARNING ** CID ',cid,' FAILS CID CUT '
         CALL FLUSH(6)
      endif

    ENDIF

    RETURN
  END SUBROUTINE PARSE_CID 


! ============================================
    INTEGER FUNCTION GET_IDSURVEY(SURVEY)


    USE SNPAR
    USE CTRLCOM

    IMPLICIT NONE

    CHARACTER SURVEY*(*)  ! (I) name of survey

    INTEGER ID, j1, j2, i
    LOGICAL LMATCH
! ---------- BEGIN ------------

    ID = -9 ; GET_IDSURVEY=-9

    DO 101 i = 1, NSURVEY_LIST
       j1 = index(SURVEY,' ') - 1
       j2 = index(SURVEY_NAME_LIST(i),' ') - 1
       if ( j1 .NE. j2 ) goto 101
       LMATCH = ( SURVEY(1:j1) .EQ. SURVEY_NAME_LIST(i)(1:j2) )
       if ( LMATCH )  ID = IDSURVEY_LIST(i)
101   CONTINUE

    GET_IDSURVEY = ID

    RETURN
  END FUNCTION GET_IDSURVEY



! =============================================
    INTEGER FUNCTION CIDASSIGN(CCID)
! 
! Created May 28, 2008 by R.Kessler
! Assign integer CID to character string CCID
! If CCID is already an integer, then CID = CCID;
! otherwise set CID = IFILE.
! 
! Note that if a string-named SN appears in
! multiple PHOTOMETRY version, it can have
! a different CID in each case. For integer
! CIDs (like SDSS), the CID is always the same.
! 
! 
! ---------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER CCID*(*)  ! (I)

    INTEGER IERR, CID, LCHAR

! ----------- BEGIN -------------

    CID = -1 ; IERR=0
    read( ccid, 20 , iostat = IERR ) CID
20    format(I10)    ! was I8

    IF ( IERR .EQ. 0 ) THEN  ! CCID is an integer
       SNLC_CASTCID = 'INT'
    ELSE
       LCHAR = INDEX ( CCID, ' ' ) - 1
       CID   = ABSO_INDEX
       SNLC_CASTCID = 'CHAR'
    ENDIF

    CIDASSIGN = CID

    RETURN
  END FUNCTION CIDASSIGN



! ==============================
    SUBROUTINE PARSE_FILTSTRING(OPT,string, NFILT, IFILTDEF, XFILTVAL)

! -------------------------------
! Created May 2008 by R.Kessler
! 
! Translate input STRING into filter indices
! and values. Assume STRING is of the form
! 
!  STRING = 'cfilt1 xval1  cfilt2  xval2 cfilt2 xval2 ...'
! 
! where a float-value follows each single-char
! filter-string.
! 
! Oct 26 2015: replace MXFILT_OBS arg with NFILT
! Dec 27 2015: fix bug initializing IFILTDEF and XFILTVAL
!              Loop 1-NFILT instead of 1-MXFILT_OBS
! 
! May 20 2016: declare output arrays (MXFILT_OBS) instead of (NFILT)
! 
! Feb 20 2017: update to allow multiple bands gluded together.
!              'gri .01'  is equivalent to  'g .01  r .01  i .01'
! 
! Feb 4 2020:
!   + check 'ALL' or 'all' to set same mag for all bands.
!   + fix bug parsing  'abcd 0.02' ; last band value wasn't set.
! 
! --------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

! subroutine args
    CHARACTER string*(*)  ! (I) string to parse

    INTEGER  & 
         OPT                   &  ! (I) 1=> error check obs filter
        ,NFILT                 &  ! (O) number of filters in string
        ,IFILTDEF(MXFILT_OBS) ! (O) absolute filter indices

    REAL  & 
         XFILTVAL(MXFILT_OBS)   ! (O) float values

! local var

    INTEGER iwd, NWD, ifilt_tmp, LEN, NFILT_TMP, NFILT_LAST, i
    INTEGER MSKOPT
    CHARACTER cwd*(MXCHAR_FILEWORD), band*2, FNAM*30
    REAL VAL
    LOGICAL LERR

! function
    INTEGER FILTINDX
    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS
! --------------- BEGIN ------------

    NFILT = 0
    IFILT_TMP = 0

    FNAM = 'PARSE_FILTSTRING'
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    NWD = STORE_PARSE_WORDS(MSKOPT,STRING//char(0),  & 
                FNAM//char(0), 80, 30)

    DO i = 1, MXFILT_OBS
      IFILTDEF(i) = 0
      XFILTVAL(i) = 0.  ! init output to zero
    ENDDO

! - - - - - - - - - - - - -
! Feb 2020: check 'all' option to set same mag error for all filters
    IF ( STRING(1:3) .EQ. 'ALL' .or. STRING(1:3).EQ. 'all' ) then
       iwd = 2
       CALL get_PARSE_WORD_fortran(iwd,cwd,LEN) ! read mag_err
       read(cwd,*) VAL
       do ifilt_tmp = 1, NFILTDEF_SURVEY
         IFILTDEF(ifilt_tmp) = IFILTDEF_MAP_SURVEY(ifilt_tmp)
         XFILTVAL(ifilt_tmp) = VAL
       enddo
       NFILT = NFILTDEF_SURVEY
       RETURN
    ENDIF

! - - - - - - - - -

    NFILT_LAST = 0
    DO 200 iwd = 1, NWD

       CALL get_PARSE_WORD_fortran(iwd,cwd,LEN)

       if ( IFILT_TMP == 0 ) then
          NFILT_LAST = NFILT
          NFILT_TMP  = LEN
          do i = 1, NFILT_TMP
            BAND            = cwd(i:i)
            NFILT           = NFILT + 1
            IFILT_TMP       = FILTINDX(BAND)
            IFILTDEF(NFILT) = IFILT_TMP

! require valid observer-filter if OPT=1
            LERR = IFILTDEF_INVMAP_SURVEY(IFILT_TMP) .LE. 0
            if ( LERR .and. OPT==1 ) then
              print*,'  SURVEY_FILTERS = ', SURVEY_FILTERS
              LEN = index(cwd,' ') - 1
              c1err = 'Invalid obs-frame filter = '''  & 
                   // BAND // ''' from nml string'
              c2err = '= ''' // string(1:60) // '''  '
              CALL MADABORT('PARSE_FILTSTRING', c1err, c2err )
            endif
         enddo
       else
          read(cwd,*) VAL
          do i = NFILT_LAST+1, NFILT_LAST+NFILT_TMP
             XFILTVAL(i) = VAL
          enddo
          ifilt_tmp = 0
       endif
200   CONTINUE


    RETURN
  END SUBROUTINE PARSE_FILTSTRING


! ====================================================
    INTEGER FUNCTION PARSE_INTLIST(STRING_INLIST,OUTLIST)

! Created Jan 25 2018
! Translate comma-separated list of integers in STRING_LIST
! into array of integers (OUTLIST).
! Function returns number of integers in OUTLIST array.
! 
! Aug 25 2025:
!   STRING_LIST length of 60 -> MXCHAR_PATH[=120]
!   MXINT=10 -> 50
! 

    USE SNPAR

    IMPLICIT NONE

    INTEGER, PARAMETER :: MXINT=50

    CHARACTER STRING_INLIST*(MXCHAR_PATH)  ! (I) comma-separated list of integers
    INTEGER   OUTLIST(MXINT)      ! (O) output array of integers

! local arg

    INTEGER iwd, NWD, LEN, MSKOPT
    CHARACTER CWD*(MXCHAR_PATH), C1ERR*72, C2ERR*72, FNAM*16
! function
    INTEGER STORE_PARSE_WORDS

! ---------------- BEGIN ----------------
    FNAM   = 'PARSE_INTLIST'
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    LEN    = INDEX(STRING_INLIST,' ') - 1
    NWD    = STORE_PARSE_WORDS(MSKOPT, STRING_INLIST//char(0),  & 
                  FNAM//char(0), LEN, 16)

    IF ( NWD > MXINT ) THEN
       CALL PRINT_PREABORT_BANNER(FNAM(1:14)//char(0),14)
       print*,'     STRING_INLIST = ', STRING_INLIST
       print*,'     LEN(STRING_INLIST) = ', LEN
       write(c1err,61) NWD, MXINT
 61      format('NWD=',I2,' exceeds bound of MXINT=',I2)
       c2err = 'Check input args'
       CALL MADABORT("PARSE_INTLIST", c1err, c2err )
    ENDIF

    DO 200 iwd = 1, NWD
       CALL get_PARSE_WORD_fortran(iwd,cwd,LEN)
       READ(cwd,*) OUTLIST(iwd)
! c         print*,' xxx ', iwd, OUTLIST(iwd)
 200  CONTINUE

    OUTLIST(NWD+1) = 0

    PARSE_INTLIST = NWD
    RETURN
  END FUNCTION PARSE_INTLIST

! ===========================================
    SUBROUTINE PARSE_COMMASEP_LIST(KEY,LINE)

! Created May 30 2019
! Parse command LINE argument and load global array
! based on KEY value.
! E.g., KEY = 'SNCCID_LIST and  LINE = '2004hq,2006ab' ->
!    SNCCID_LIST = '2004hq', '2006ab'
! 


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER KEY*(*), LINE*(*)  ! (I)

    INTEGER  MSKOPT, NWD, iwd, LL
    INTEGER  STORE_PARSE_WORDS
    CHARACTER WDTMP*60, FNAM*30
    EXTERNAL STORE_PARSE_WORDS

! --------------- BEGIN ------------

    FNAM = 'PARSE_COMMASEP_LIST'
    LL = INDEX(LINE,' ') - 1
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    NWD    = STORE_PARSE_WORDS(MSKOPT,LINE(1:LL)//char(0),  & 
                   FNAM//char(0), LL, 30)

    DO iwd = 1, NWD

       if ( KEY .EQ. 'SNCCID_LIST' ) then
          CALL get_PARSE_WORD_fortran(iwd,  & 
                         SNCCID_LIST(iwd), LL)

       else if ( KEY .EQ. 'SNCCID_IGNORE' ) then
          CALL get_PARSE_WORD_fortran(iwd,  & 
                         SNCCID_IGNORE(iwd), LL)

       else if ( KEY .EQ. 'VERSION_PHOTOMETRY' ) then
          CALL get_PARSE_WORD_fortran(iwd,  & 
                         VERSION_PHOTOMETRY(iwd), LL)

       else if ( KEY .EQ. 'PHOTFLAG_BITLIST' ) then
          CALL get_PARSE_WORD_fortran(iwd, WDTMP, LL)
          read(WDTMP,*) PHOTFLAG_BITLIST_REJECT(iwd)
       endif

!        LL = INDEX(SNCCID_LIST(iwd),' ') - 1
!        print*,' xxx SNCCID = |', SNCCID_LIST(iwd)(1:LL), '|'
    ENDDO

    RETURN
  END SUBROUTINE PARSE_COMMASEP_LIST

! ====================================================
    INTEGER FUNCTION PARSE_NML_STRINGLIST(STRLIST,NCHAR)

! Created Mar 12 2015 by R.Kessler
! Function returns number of non-null elements,
! and parses strings separated by blanks.
! 
! The following inputs all produce the same output
! 
! Input:
!      STRLIST = 'E1', 'E2', 'S1', 'S2'
!          or
!      STRLIST = 'E1 E2 S1 S2'
!          or
!      STRLIST = 'E1 E2',  'S1 S2'

! all produce  Output:
!    STRLIST = 'E1', 'E2', 'S1', 'S2'
!    and FUN=4
! 
! ------------------------------


    USE SNPAR

    IMPLICIT NONE

    INTEGER NCHAR   ! (I)
    CHARACTER STRLIST(MXLISTNML)*(NCHAR) ! (I,O)

! local var

    CHARACTER STRLIST_LOCAL(MXLISTNML)*(NCHAR), FNAM*30
    INTEGER i, j, NLIST, iwd, NWD, LL, MSKOPT

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------- BEGIN ------------

    FNAM = 'PARSE_NML_STRINGLIST'

    NLIST = 0
    MSKOPT = MSKOPT_PARSE_WORDS_STRING

! transfer input string to local string
    DO i = 1, MXLISTNML
       STRLIST_LOCAL(i) = STRLIST(i)
    ENDDO

    i = 1
    DO WHILE ( STRLIST_LOCAL(i) .NE. '' )

! remove optional '+' symbols between strings, and substitute blank.
       DO j = 1, NCHAR
         IF ( STRLIST_LOCAL(i)(j:j) .EQ. '+' ) then
            STRLIST_LOCAL(i)(j:j) = ' '
         endif
       ENDDO

       NWD = STORE_PARSE_WORDS(MSKOPT,STRLIST_LOCAL(i)//char(0),  & 
                  FNAM//char(0), 100, 30 )

       DO iwd = 1, NWD
          NLIST = NLIST + 1
          CALL get_PARSE_WORD_fortran(iwd,STRLIST(NLIST),LL)
       ENDDO

       i = i + 1

    END DO  ! end while

    PARSE_NML_STRINGLIST = NLIST

    RETURN
  END FUNCTION PARSE_NML_STRINGLIST

! ====================================================
    SUBROUTINE CIDSTRING(CID,CCID,LENCCID)
! 
! Created Jan 15, 2011
! If integer CID < 1 million, return 6-char CCID and LENCCID=6;
! otherwise return 8-char CCID and LENCCID=8.
! 

    USE SNPAR

    IMPLICIT NONE

    INTEGER   CID                ! (I) integer cand. id
    CHARACTER CCID*(MXCHAR_CCID)  ! (O) char-string for CCID
    INTEGER   LENCCID            ! (O) length of string

! local
    CHARACTER C1ERR*72, C2ERR*72

! ----------- BEGIN -----------

    IF ( CID .LT. 1000000 ) THEN    ! 1 million
       write(ccid,26) CID
26       format(I6.6)
       LENCCID = 6
    ELSE IF ( CID < 100000000 ) THEN  ! 100 million
       write(ccid,28) CID
28       format(I8.8)
       LENCCID = 8
    ELSE IF ( CID < MXCID ) THEN
       write(ccid,29) CID
29       format(I9.9)
       LENCCID = 9
    ELSE
       write(C1ERR,66) CID, MXCID
66       format('CID=',I9,' exceeds bound (MXCID=', I9,')' )
       C2err = '   '
       CALL MADABORT("CIDSTRING", c1err, c2err )
    ENDIF

    RETURN
  END SUBROUTINE CIDSTRING

! ===============================================
    INTEGER FUNCTION FILTINDX_REPLACE ( cfilt )

! Created Jun 2013 by RK.
! Same as FILTINDX, but apply FILTER_REPLACE (see INIT_FILTER_REPLACE)
! Examples:
!   FILTER_REPLACE = ''       : do exactly the same thing as FILTINDX.
!   FILTER_REPLACE = 'U -> u' : for U filter, return index for 'u'
! 
! ------


    USE SNPAR
    USE FILTCOM

    IMPLICIT NONE

    character cfilt*(*)  ! (I) filter name to parse

    INTEGER   IFILTDEF_ORIG
!  function
    INTEGER FILTINDX

! ---------------- BEGIN ---------------
    IFILTDEF_ORIG    =  FILTINDX(CFILT)
    FILTINDX_REPLACE =  IFILTOBS_REPLACE(IFILTDEF_ORIG)

    RETURN
  END FUNCTION FILTINDX_REPLACE 

! ===============================================
    INTEGER FUNCTION FILTINDX ( cfilt )
! 
! May 2008 R.Kessler
! Returns absolute integer filter-index 1:MXFLT_ALL by
! parsing character name "cfilt".
! Assumes that last character is the 1-char symbol;
! i.e, SDSS-g => g,  CTIO4m-R => R, etc ...


    USE SNPAR
    USE FILTCOM

    IMPLICIT NONE

    character cfilt*(*)  ! (I) filter name to parse

    integer LL, ifilt
    character cfilt1*1, ctest*1


! ----------- BEGIN -----------

    FILTINDX = 0

    LL  = index(CFILT//' ',' ' ) - 1
    cfilt1 = CFILT(LL:LL)

    DO ifilt = 1, MXFILT_ALL
      ctest = filtdef_string(ifilt:ifilt)
      if ( ctest(1:1) .eq. cfilt1(1:1) ) then
         FILTINDX = ifilt
         return
      endif
    ENDDO

    RETURN
  END FUNCTION FILTINDX 

! ==============================
    LOGICAL FUNCTION ISBXFILT(IFILT,NAME,FRAME)

! Created May 2012
! Set global locical EXIST_BXFILT_OBS[REST]=T if this
! X-filter contains 'BX'.
! 


    USE SNPAR
    USE FILTCOM

    IMPLICIT NONE

! function args

    INTEGER IFILT  ! (I) absolute filter index
    CHARACTER  & 
         NAME*(*)     &  ! (I) full name of filter
        ,FRAME*(*)   ! (I) 'OBS' or 'REST'

! local args

    INTEGER JX
    CHARACTER C2*2

! ----------------- BEGIN ----------------

    ISBXFILT = .FALSE.
    IF ( IFILT .NE. IFILT_BESS_BX ) RETURN

    JX = INDEX(NAME,'X')

    IF ( JX .GT. 1 ) THEN
        C2 = NAME(JX-1:JX)
        IF( C2 .EQ. 'BX' ) then
           ISBXFILT = .TRUE.
        else
           return
        endif

! now that ISBXFILT=T  set global logical
        IF ( FRAME .EQ. 'OBS'  ) EXIST_BXFILT_OBS  = .TRUE.
        IF ( FRAME .EQ. 'REST' ) EXIST_BXFILT_REST = .TRUE.

    ELSE

        IF ( FRAME .EQ. 'OBS'   ) ISBXFILT = EXIST_BXFILT_OBS
        IF ( FRAME .EQ. 'REST'  ) ISBXFILT = EXIST_BXFILT_REST

    ENDIF

    RETURN
  END FUNCTION ISBXFILT



! ================================================
    INTEGER FUNCTION GET_IDFIELD ( field_name )
! 
! May 2009
! Return integer ID for this FIELD.
! IDFIELD is defined $SNDATA_ROOT/SURVEY.DEF
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    character field_name*(*)  ! name of telescope

! local
    character field_upper*60
    INTEGER LL, i
! --------------- BEGIN -----------------

    GET_IDFIELD = -9

    CALL UPCASE(field_name, field_upper) ! returns field_upper

    if ( field_upper(1:3) .EQ. 'ALL' ) then
      GET_IDFIELD = 0
      return
    endif

    DO i = 1, NIDFIELD_LIST  ! user-list of fields
      LL = INDEX(SNFIELD_LIST(i),' ' ) - 1
      IF ( FIELD_UPPER(1:LL) .EQ. SNFIELD_LIST(i)(1:LL) ) THEN
        GET_IDFIELD = IDFIELD_LIST(i)
      ENDIF
    END DO

    RETURN
  END FUNCTION GET_IDFIELD 

! ==================================
    SUBROUTINE ENDMJD_PROC(NEWMJD)
! 
! Aug 2007, R.Kessler
! 
! Erase bad mesaurements.
! If all measurements (filters) are erased, then erase
! this MJD.
! 
! May 9 2016: modify IGNORE logic
! 
! ---------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

! input args
    INTEGER NEWMJD  ! index for SN and MJD

! local args

    INTEGER  & 
         EPMIN, EPMAX, ep, ep2  & 
        ,ifilt_obs, ifilt_obs2, ifilt  & 
        ,NEPSKIP, NFILT, EPMAX_OLD, NMOVED

    REAL FLUXCAL, FLUXCAL_ERR

    LOGICAL  & 
         LVALID(MXFILT_ALL)  & 
        ,LMOVED(MXFILT_ALL)  & 
        ,LTMP, IGNORE

    character cfilt1

! ------------- BEGIN -------------

    NFILT = ISNLC_NFILT_NEWMJD(NEWMJD)

! for this NEWMJD, erase NULL measurements to save
! memory for old files with 5 filters per _line even
! when only 1 filter has data.
! Note that valid measurements are re-sorted
! and EPMAX is lowered.

    EPMIN  =  ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
    EPMAX  =  ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

    do ifilt_obs = 1, MXFILT_ALL
        LVALID(ifilt_obs) = .FALSE.
        LMOVED(ifilt_obs) = .FALSE.
    enddo

    NEPSKIP = 0

    DO 100 ep      = EPMIN, EPMAX
       ifilt_obs   = ISNLC_IFILT_OBS(ep)
       fluxcal     = SNLC_FLUXCAL(ep)
       fluxcal_err = SNLC_FLUXCAL_ERRTOT(ep)
!          mag         = SNLC_MAG(ep)
       IGNORE      = (fluxcal_err < 0.0 ) .and.  & 
               (.not. REFORMAT_SAVE_BADEPOCHS )

! make sure filter is valid by checking that
! the flux-error is non-negative.

       IFILT = IFILTDEF_INVMAP_SURVEY(IFILT_OBS)
       if ( IFILT .LT. 1 .or. IFILT .GT. NFILTDEF_SURVEY ) then

          cfilt1 = FILTDEF_STRING(ifilt_obs:ifilt_obs)

          write(c1err,650)  & 
               ifilt_obs, cfilt1, SNLC8_MJD(EPMIN), SNLC_CCID
650         format('Invalid IFILT_OBS=',I3, '(',A,')' ,  & 
                '  at MJD=',F9.3, 3x,'CID=',A)

          write(c2err,651) SURVEY_FILTERS(1:NFILTDEF_SURVEY)
651         format('Valid IFILT_OBS  are ', A )

          CALL MADABORT("ENDMJD_PROC", c1err, c2err )
       endif

       LMOVED(ifilt_obs) = .FALSE.

! check for valid measure.
! Set LSORT=T for any invalid measurement

       LTMP  =  (.not. IGNORE )  & 
           .and.  (CUTWIN_SNRMIN_FILT(1,ifilt) < 100.0)

       if ( LTMP ) then
          LVALID(ifilt_obs) = .TRUE.
       else
         NEPSKIP =  NEPSKIP + 1
         LVALID(ifilt_obs) = .FALSE.
       endif
100   CONTINUE

!       print*,' xxx LVALID = ', LVALID

! -------------------------------------------------------
! increment NEWMJD and NEPOCH if we have not erased
! all measurements

    if ( NEPSKIP .EQ. NFILT ) RETURN ! erase entire MJD

! ----------------

    ISNLC_NEWMJD_STORE = NEWMJD  ! increment new MJD

    ISNLC_NEPOCH_STORE = ISNLC_NEPOCH_STORE +  & 
         ISNLC_NFILT_NEWMJD(NEWMJD)

    if ( NEPSKIP .EQ. 0 ) RETURN  ! all OK => bail

! if we get here, then there is at least one invalid
! measurement (epoch) for this NEWMJD; re-sort
! epochs to skip bad epochs.

    EPMAX_OLD = EPMAX
    EPMAX     = EPMAX_OLD - NEPSKIP

    IF ( EPMAX .LT. EPMIN ) then

      write(c1err,660) EPMIN,EPMAX_OLD, EPMIN, EPMAX,  & 
             NFILT, NEPSKIP
660     format('EPMIN,MAX=',2I4, ' -> ', 2I4, 3x,  & 
             'NEPFILT=',I2, 2x,  'NEPSKIP=',I2  )
      write(c2err,661) SNLC8_MJD(EPMIN), SNLC_CCID
661     format('at MJD = ', F9.3, 4x, 'CID=', A )

      CALL MADABORT("ENDMJD_PROC", c1err, c2err )

    ENDIF


! loop over epochs and re-sort with NULL measurements removed.
! Logic is tricky.

    DO 200 ep = EPMIN, EPMAX_OLD

       ifilt_obs   = ISNLC_IFILT_OBS(ep)
       if ( LVALID(ifilt_obs) ) goto 200

! here we have invalid measure at epoch EP ...
! copy next valid epoch here.

       NMOVED = 0
       DO 202 ep2 = ep+1, EPMAX_OLD
         ifilt_obs2   = ISNLC_IFILT_OBS(ep2)
         LTMP = .NOT. LMOVED(ifilt_obs2) .and. NMOVED.EQ.0
         if ( LVALID(ifilt_obs2) .and. LTMP ) then
            CALL MOVE_SNLC_ARRAYS(ep,ep2)  ! move ep2 -> ep
            LMOVED(ifilt_obs2) = .TRUE.
            LVALID(ifilt_obs2) = .FALSE.
            NMOVED = NMOVED + 1
         endif

202      CONTINUE

200   CONTINUE  ! end of 'ep' loop

! ---------------------------------------
! adjust global variables that depend on NEPSKIP or EPMAX

     ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD) = EPMIN
     ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD) = EPMAX

     ISNLC_NEPOCH_STORE =  & 
       ISNLC_NEPOCH_STORE - NEPSKIP

     ISNLC_NEPOCH_FOUND =  & 
       ISNLC_NEPOCH_FOUND - NEPSKIP

     ISNLC_NFILT_NEWMJD(NEWMJD) =  & 
       ISNLC_NFILT_NEWMJD(NEWMJD) - NEPSKIP

    RETURN
  END SUBROUTINE ENDMJD_PROC

! ===================================
    SUBROUTINE MOVE_SNLC_ARRAYS(ep1,ep2)
! 
! 
! Jun 25 2019: include SIM_XXX arrays
! Nov 12 2019: check LSIM_MAGOBS to set SIM_EPMAGOBS(ep)
! Feb 17 2021:
!   for I/O refactor, include SNLC_ZEROPT_forCUTS and SNLC_PSF_FWHM
! Mar 01 2021: include SNLC_PSF_NEA


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER ep1, ep2  ! (I)

! ------------- BEGIN ------------

    SNLC8_MJD(ep1)           = SNLC8_MJD(ep2)
    ISNLC_IDFIELD(ep1)       = ISNLC_IDFIELD(ep2)
    SNLC_FIELD(ep1)          = SNLC_FIELD(ep2)
    ISNLC_IFILT_OBS(ep1)     = ISNLC_IFILT_OBS(ep2)
    SNLC_FLUXCAL(ep1)        = SNLC_FLUXCAL(ep2)
    SNLC_FLUXCAL_ERRTOT(ep1) = SNLC_FLUXCAL_ERRTOT(ep2)
    ISNLC_PHOTFLAG(ep1)      = ISNLC_PHOTFLAG(ep2)
    SNLC_PHOTPROB(ep1)       = SNLC_PHOTPROB(ep2)
    SNLC_PSF_SIG1(ep1)       = SNLC_PSF_SIG1(ep2)
    SNLC_PSF_SIG2(ep1)       = SNLC_PSF_SIG2(ep2)
    SNLC_PSF_RATIO(ep1)      = SNLC_PSF_RATIO(ep2)
    SNLC_PSF_NEA(ep1)        = SNLC_PSF_NEA(ep2)
    SNLC_SKYSIG(ep1)         = SNLC_SKYSIG(ep2)
    SNLC_SKYSIG_T(ep1)       = SNLC_SKYSIG_T(ep2)
    SNLC_RDNOISE(ep1)        = SNLC_RDNOISE(ep2)
    SNLC_ZEROPT(ep1)         = SNLC_ZEROPT(ep2)
    SNLC_ZEROPT_ERR(ep1)     = SNLC_ZEROPT_ERR(ep2)
    SNLC_TEXPOSE(ep1)        = SNLC_TEXPOSE(ep2)
    SNLC_GAIN(ep1)           = SNLC_GAIN(ep2)
    SNLC_XPIX(ep1)           = SNLC_XPIX(ep2)
    SNLC_YPIX(ep1)           = SNLC_YPIX(ep2)

    SNLC_ZEROPT_forCUT(ep1)   = SNLC_ZEROPT_forCUT(ep2)
    SNLC_PSF_FWHM_ARCSEC(ep1) = SNLC_PSF_FWHM_ARCSEC(ep2)

! beware that DETNUM is read from header, and also from
! PHOT section of TEXT format ... but not read from PHOT-FITS file??
    ISNLC_DETNUM(EP1)        = ISNLC_DETNUM(EP2)
    ISNLC_IMGNUM(EP1)        = ISNLC_IMGNUM(EP2)

    IF ( ISJOB_SIM ) THEN
      SIM_EPMAGOBS(EP1)          = SIM_EPMAGOBS(EP2)
      SIM_EPFLUXCAL(EP1)         = SIM_EPFLUXCAL(EP2)
      SIM_EPCHI2FLUX(EP1)        = SIM_EPCHI2FLUX(EP2)
      SIM_EPPULL(EP1)            = SIM_EPPULL(EP2)
    ENDIF

    IF ( LSIM_SNANA ) THEN
      SIM_EPFLUXCAL_HOSTERR(EP1) = SIM_EPFLUXCAL_HOSTERR(EP2)
      SIM_EPSNRMON(EP1)          = SIM_EPSNRMON(EP2)
    ENDIF

! - - -
    SNLC_MAG(ep1)            = SNLC_MAG(ep2)
    SNLC_MAG_ERRPLUS(ep1)    = SNLC_MAG_ERRPLUS(ep2)
    SNLC_MAG_ERRMINUS(ep1)   = SNLC_MAG_ERRMINUS(ep2)

    RETURN
  END SUBROUTINE MOVE_SNLC_ARRAYS

! ===================================
    SUBROUTINE INIT_SNLC()
! 
! Init SNLC_XXX arrays for this SN
! Called for each SN.
! 
! Nov 22 2017: use NEP_RESET to reduce CPU time for init
! Feb 07 2021: rename INIT_SNDATA -> INIT_SNLC to avoid confusion
!               with INIT_SNDATA in sntools.c
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT
    USE PRIVCOM
    USE EARLYCOM
    USE SPECCOM

    IMPLICIT NONE

! local

    INTEGER iep, ifilt, ipar, i, NEP_RESET, igal

! ------------- BEGIN -------------

    FOUND_SURVEY = .FALSE.

    SNLC_CID    = -9
    SNLC_CCID   = ''
    SNLC_NAME_IAUC       = 'NONE'   ! Aug 31 2017 (don't write blank to FITS file)
    SNLC_NAME_TRANSIENT  = 'NONE'   ! Jul 26 2024

    NSPECTRUM   = 0    ! for reading spectra for plot-table
    NLAMBIN_READ = 0

    ISNLC_NEWMJD_HEAD   =  0
    ISNLC_NEWMJD_FOUND  =  0
    ISNLC_NEWMJD_STORE  =  0
    ISNLC_NEWMJD_CUTS   =  0
    ISNLC_NEPOCH_FOUND  =  0
    ISNLC_NEPOCH_STORE  =  0
    ISNLC_NEPOCH_PHOTPROB = 0  ! NEPOCH with PHOTPROB >= 0
    ISNLC_NFILT_SNRMAX  =  0
    ISNLC_NFILT_SNRMAX2 =  0
    ISNLC_NFIELD_OVP    =  0

    ISNLC_NFILT_TRESTMIN =  0
    ISNLC_NFILT_TRESTMAX =  0
    ISNLC_NFILT_TREST2   =  0

    ISNLC_NOBS_PREDETECT   =  0
    ISNLC_NOBS_DETECT      =  0
    SNLC_TLIVE_DETECT      = -9.0
    SNLC8_MJD_TRIGGER      = -99.0
    SNLC8_MJD_DETECT_FIRST = -99.0
    SNLC8_MJD_DETECT_LAST  = -99.0

    NEPOCH_BADPHOT      =   0
    ISNLC_FAKE          =   0  ! Feb 2021: default is real data
    ISNLC_DETNUM(1)     =  -9
    ISNLC_IMGNUM(1)     =  -9

    SNLC8_RA            =  -9.0
    SNLC8_DEC           =  -9.0

    SNLC_ZHELIO        =  -9.0
    SNLC_ZHELIO_ERR    =   0.0
    SNLC_ZCMB          =  -9.0
    SNLC_ZCMB_ERR      =   0.0
    SNLC_ZSN           =  -9.0
    SNLC_ZSN_ERR       =   0.0
    SNLC_REDSHIFT      =  -9.0
    SNLC_REDSHIFT_ERR  =   0.0
    SNLC_VPEC          =   0.0
    SNLC_VPEC_ERR      =   0.0
    ISNLC_zFLAG        =  -999

    SNLC_PIXSIZE       =  -9.0
    SNLC_NXPIX         =  -9.0
    SNLC_NYPIX         =  -9.0

    SNLC_MWEBV         =   0.0
    SNLC_MWEBV_ERR     =   0.0
    SNLC_SEARCH_PEAKMJD    =  -9.0
    SNLC_SNANAFIT_PEAKMJD = -9.0
    ISNLC_TYPE         =  0

    SNLC_PHOTPROB_MIN = 1.0

    SNLC_AREAFRAC_AVG = -9.0

! init a few SNHOST arrays

    SNHOST_NMATCH            =  0
    SNHOST_NMATCH2           =  0
    SNHOST_CONFUSION         = -9.0
    do igal = 1, MXSNHOST
       SNHOST_FLAG(igal)        = 0
       SNHOST_OBJID(igal)       = -9
      DSNHOST_OBJID(igal)       = -9.0
       SNHOST8_RA(igal)         = -999.0
       SNHOST8_DEC(igal)        = -999.0
       SNHOST_ANGSEP(igal)      = -9.0
       SNHOST_DDLR(igal)        = -9.0
       SNHOST_ZPHOT(igal)       = -9.0
       SNHOST_ZPHOT_ERR(igal)   = -9.0
       SNHOST_ZSPEC(igal)       = -9.0
       SNHOST_ZSPEC_ERR(igal)   = -9.0

       SNHOST_LOGMASS(igal)       = -9999.0
       SNHOST_LOGMASS_ERR(igal)   = -9999.0
       SNHOST_LOGSFR(igal)        = -9999.0
       SNHOST_LOGSFR_ERR(igal)    = -9999.0
       SNHOST_LOGsSFR(igal)       = -9999.0
       SNHOST_LOGsSFR_ERR(igal)   = -9999.0
       SNHOST_COLOR(igal)         = -9999.0
       SNHOST_COLOR_ERR(igal)     = -9999.0
    enddo
! --------

    SNLC_SNRMAX_FILT(0) = -9.
    SNLC_SNRSUM         = 0.0

    DO ifilt = 1, MXFILT_OBS
       ISNLC_NEPOCH_FILT(ifilt)    =   0
       EXIST_FILT(ifilt)          = .FALSE.
       SNLC_FLUXCALMAX(ifilt)     = -9.
       SNLC_FLUXCALMAX_ERR(ifilt) = -9.
       SNLC_MWXT_MAG(ifilt)       = 0.0
       SNLC_MWXT_MAGERR(ifilt)    = 0.0
       SNLC_MWXT_FLUXFRAC(ifilt)  = 0.0

       SNHOST_SBFLUXCAL(ifilt)      = -999.9
       SNHOST_SBFLUXCAL_ERR(ifilt)  = -999.9
       SNHOST_SBMAG(ifilt)          = -999.9
       do igal = 1, MXSNHOST
          SNHOST_MAGOBS(ifilt,igal)         = +99.0
          SNHOST_MAGOBS_ERR(ifilt,igal)     = +99.0
       enddo
       SNLC_SNRMAX_FILT(ifilt)     = -9.
       SNLC_SNRMAX_SORT(ifilt)     = -9.
       SNLC_FLUXCAL_OFF(ifilt)     = 0.0

       SNLC_FIELDLIST = 'VOID'

       DO ipar = 1, NPAR_ANYLC
         SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar) = -99.0
         SNLC_SNANAFIT_PEAKMJD_FITERR(ifilt,ipar) = -99.0
       ENDDO

       SIM_TEMPLATEMAG(ifilt) = 99.0
       SIM_LCWIDTH(ifilt)     = 0.0
    ENDDO

    LSIM_TRUE_SNIa  = .FALSE.
    SIM_MODEL_INDEX = -9
    SIM_COLORPAR    = -9.
    SIM_COLORLAW    = -9.
    SIM_DLMAG       = -9.
    SIM_LENSDMU     =  0.
    SIM_MUSHIFT     =  0.

    SIM_WGT_POPULATION = 1.0
    SIM_SHAPEPAR    = -9.
    SIM_SHAPELAW    = -9.
    SIM_COLORPAR    = -9.
    SIM_COLORLAW    = -9.
    SIM_SALT2x0     = -9.
    SIM_SALT2mb     = -9.
    SIM_AV          = -9.    ! July 2016
    SIM_RV          = -9.    ! idem
    SIM_TEMPLATE_INDEX    = -9
    SIM_SEARCHEFF_MASK = 0
    SIM_MAGSMEAR_COH   = 0.0
    SIM_SALT2gammaDM   = 0.0

    SIM_NGEN_LIBID = 0
    SIM_NOBS_UNDEFINED  = 0
    SIM_SUBSAMPLE_INDEX = -9

    IF ( NCALL_SNANA_DRIVER < 2 ) then
       NEP_RESET = MXEPOCH
    ELSE
       NEP_RESET = ISNLC_NEPOCH_STORE
    ENDIF
!      print*,' xxx NCALL,NEP_RESET=', NCALL_SNANA_DRIVER, NEP_RESET

    DO iep = 1, NEP_RESET

      SNLC_FIELD(iep) = 'VOID'   ! avoid pandas issues with NULL

      ISNLC_NFILT_NEWMJD(iep)          =  0
      ISNLC_EPOCH_RANGE_NEWMJD(1,iep)  =  -9
      ISNLC_EPOCH_RANGE_NEWMJD(2,iep)  =  -9
      ISNLC_IFILT_OBS(iep)      =  -9

      ISNLC_IDFIELD(iep)        =  -9

      SNLC8_MJD(iep)           =  -9.0

      ISNLC_PHOTFLAG(iep)  = 0
      SNLC_PHOTPROB(iep)   = 0.0

      SNLC_SKYSIG(iep)     = 0.0
      SNLC_SKYSIG_T(iep)   = 0.0  ! Mar 2018
      SNLC_PSF_SIG1(iep)   = 0.0
      SNLC_PSF_SIG2(iep)   = 0.0
      SNLC_PSF_RATIO(iep)  = 0.0
      SNLC_PSF_NEA(iep)    = 0.0
      SNLC_GAIN(iep)       = 0.0
      SNLC_RDNOISE(iep)    = 0.0

      SNLC_SNR(iep)             =  0.0
      SNLC_FLUXCAL_ERRCALC(iep) =  0.0
      SNLC_FLUXCAL_ERRTEST(iep) = -9.0
      SNLC_FLUXCAL_HOSTERRCALC(iep) = 0.0
      SNLC_FLUXCAL(iep)         =  0.0
      SNLC_FLUXCAL_ERRTOT(iep)  = -9.0  ! Oct 8 2014
      SNLC_MAG(iep)             =  0.0
      SNLC_MAG_ERRPLUS(iep)     =  0.0
      SNLC_MAG_ERRMINUS(iep)    =  0.0

      SNLC_ZEROPT(iep)          =  0.0
      SNLC_ZEROPT_ERR(iep)      =  0.0
      SNLC_TEXPOSE(iep)         =  0.0

      SNLC_XPIX(iep) = -9.0
      SNLC_YPIX(iep) = -9.0

      SNLC_DTOBS(iep)          = -9.0
      SNLC_DTOBS_SAMEFILT(iep) = -9.0

      SNLC_AREAFRAC(iep)  = -9.0

    END DO  ! end of epoch loop

! zero out the FITVAL_STORE array if fit-option is chosen

    IF ( DO_FIT ) THEN

         FITPROB_ITER1   = 0.0  ! July 2024
	   FITCHI2RED_INI  = 0.0
	   FITCHI2RED_INI2 = 0.0

       do i = 1, 4
         FITCHI2_STORE(i)      = 0.0
         FITPROBCHI2_STORE(i)  = 0.0
         LCCHI2_STORE(i)       = 0.0
         LCPROBCHI2_STORE(i)   = 0.0

         NDOF_STORE(i)  = 0
       enddo

       NDOF_PRIOR = 0

       do ipar = 1, MXFITSTORE
          FITVAL_STORE(ipar) = 0.0
          FITERR_STORE(ipar) = 0.0

          LCVAL_STORE(ipar) = 0.0
          LCERR_STORE(ipar) = 0.0
          LCFRACERRDIF_STORE(ipar) = 0.0
       enddo
    ENDIF


    DO i   = 1, MXERRTYPE
       NERRTYPE(i) = 0
    ENDDO

    DO i = 1, MXVAR_PRIVATE
      PRIVATE_VALUE(i)    = PRIVATE_NULL
    ENDDO

    LSNCUTS       = .FALSE.
    PASS_PRIVCUTS = .TRUE.
    PASS_SIMCUTS  = .TRUE.
    CUTFLAG_SNANA =  0
    ERRFLAG_FIT   = -9  ! init to no fit

    MADE_LCPLOT = .FALSE.

    NSEASON_TOT      = 0
    NSEASON_ACTIVE    = 0

    ISNLC_CUTFLAG_REQEP    = 1
    ISNLC_CUTFLAG_PRIVATE  = 1
    ISNLC_NOBS_TREST       = 0
    ISNLC_CUTFLAG_SIMVAR   = 1
    ISNLC_WRMASK_FLUXCOR_SNANA = 0
    ISNLC_RDMASK_FLUXCOR_SNANA = 0

    NOBS_EARLYLC    = 0
    NNIGHT_EARLYLC  = 0
    MJDLAST_EARLYLC = -999.
    MJDLAST_SELECT  = -999.
    NPHOTMASK_START_EARLYLC = 0
    NSNR_START_EARLYLC      = 0

    DO i = 1, NPAR_SIMSED
       SIMSED_PARVAL(i) = -999.0
    ENDDO
    DO i = 1, NPAR_LCLIB
       LCLIB_PARVAL(i) = -999.0
    ENDDO

    RETURN
  END SUBROUTINE INIT_SNLC


! =======================================
    SUBROUTINE INIT_CUTMASK ( IERR )
! 
! Creatd Jan 31, 2006 by R.Kessler
! Initialize CUTMASK_ALL and CUTWIN_XXX
! 
! Jan 4 2021:
!  + remove obsolete logic with !LSIM_SNANA && ibit==CUTBIT_SEARCH
!  + check OPT_SNCID_LIST
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IERR ! (O)  0=> OK
    INTEGER ibit

! -------------- BEGIN -------------------

    CALL PRBANNER('INIT_CUTMASK')

    IERR = 0  ! init output arg

    if ( NCUTBIT .GE. 64 ) then
       print*,' INIT_CUTBIT FATAL ERROR: NCUTBIT=', NCUTBIT
       print*,' ***** ABORT ***** '
       CALL EXIT(EXIT_ERRCODE)
    endif

    CUTMASK8_SN_ALL    = 0  ! init common block var
    CUTMASK8_MJD_ALL   = 0  ! init common block var

    DO 100 ibit = 1, NCUTBIT
       if ( ibit .LE. CUTBIT_MJD_MARKER ) then
         CUTMASK8_SN_ALL = IBSET ( CUTMASK8_SN_ALL, ibit-1 )
       else
         CUTMASK8_MJD_ALL = IBSET ( CUTMASK8_MJD_ALL, ibit-1 )
       endif
100   CONTINUE


! Jan 2021: check option to require ONLY CID-select.

    if ( USE_SNCID_FILE ) THEN
        CUTMASK8_SN_ALL = 0
        CUTMASK8_SN_ALL = IBSET(CUTMASK8_SN_ALL,CUTBIT_CID-1)
        print*,' '
        print*,'      **** APPLY ONLY CID-SELECT CUT **** '
        print*,' '
        CALL FLUSH(6)
    endif

    print*,' INIT_CUTMASK: NCUTBIT          = ', NCUTBIT
    print*,' INIT_CUTMASK: CUTMASK8_SN_ALL  = ',CUTMASK8_SN_ALL
    print*,' INIT_CUTMASK: CUTMASK8_MJD_ALL = ',CUTMASK8_MJD_ALL
    call flush(6)

! - - - - - - - - -  -
! July 11 2020:
!  check option to write list of rejected epochs
!   (e.g., to pass to classifier)

    CALL WRITE_EPOCH_IGNORE_INIT()

    RETURN
  END SUBROUTINE INIT_CUTMASK 

! =============================================
    SUBROUTINE WRITE_EPOCH_IGNORE_INIT()


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER LEN
    REAL*8  CUTWIN_TMP(2)
    CHARACTER CUTMODE*12, CUTVARNAME*20
    EXTERNAL write_epoch_list_init, write_epoch_list_addvar

! ---------- BEGIN --------

    LEN = INDEX(OUT_EPOCH_IGNORE_FILE,' ') - 1
    if ( LEN <= 0 ) RETURN

    CALL write_epoch_list_init(  & 
           OUT_EPOCH_IGNORE_FILE(1:LEN)//char(0), LEN+1)

    CUTMODE = 'REJECT' // char(0)
!     start with PHOTFLAG cut
    CUTWIN_TMP(1)  = DBLE(PHOTFLAG_MSKREJ(1))
    CUTWIN_TMP(2)  = -9.0
    CUTVARNAME     = 'PHOTFLAG' // char(0)
    CALL write_epoch_list_addvar(CUTVARNAME, CUTWIN_TMP, CUTMODE,  & 
           20, 12)

! next check PHOTPROB
    CUTWIN_TMP(1)  = CUTWIN_PHOTPROB(1)
    CUTWIN_TMP(2)  = CUTWIN_PHOTPROB(2)
    CUTVARNAME     = 'PHOTPROB' // char(0)
    CALL write_epoch_list_addvar(CUTVARNAME, CUTWIN_TMP, CUTMODE,  & 
           20, 12)

! zero point
    CUTWIN_TMP(1)  = CUTWIN_ZP(1)
    CUTWIN_TMP(2)  = CUTWIN_ZP(2)
    CUTVARNAME     = 'ZP' // char(0)
    CALL write_epoch_list_addvar(CUTVARNAME, CUTWIN_TMP, CUTMODE,  & 
           20, 12)

! PSF
    CUTWIN_TMP(1) = CUTWIN_PSF(1)
    CUTWIN_TMP(2) = CUTWIN_PSF(2)
    CUTVARNAME    = 'PSF' // char(0)
    CALL write_epoch_list_addvar(CUTVARNAME, CUTWIN_TMP, CUTMODE,  & 
           20, 12)

!      print*,' '
!      print*,' xxx DEBUG STOP xxx '
!      STOP

    RETURN
  END SUBROUTINE WRITE_EPOCH_IGNORE_INIT

! =============================================
    SUBROUTINE WRITE_EPOCH_IGNORE_EXEC(ep)


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER ep ! index to process

    INTEGER   NVAR, IFILT_OBS
    CHARACTER CCID*(MXCHAR_CCID), BAND*2
    REAL*8  MJD, VALUES(20)
    EXTERNAL write_epoch_list_exec

! -------- BEGIN -----------
    IF ( OUT_EPOCH_IGNORE_FILE .EQ. '' ) RETURN

    NVAR = 0

    IFILT_OBS = ISNLC_IFILT_OBS(ep)
    CCID = SNLC_CCID(1:ISNLC_LENCCID) // char(0)
    MJD  = SNLC8_MJD(ep)
    BAND = FILTDEF_STRING(ifilt_obs:ifilt_obs) // char(0)

    NVAR=NVAR+1;  VALUES(NVAR) = DBLE( ISNLC_PHOTFLAG(ep) )
    NVAR=NVAR+1;  VALUES(NVAR) = DBLE( SNLC_PHOTPROB(ep)  )
    NVAR=NVAR+1;  VALUES(NVAR) = DBLE( SNLC_ZEROPT_forCUT(ep) )
    NVAR=NVAR+1;  VALUES(NVAR) = DBLE( SNLC_PSF_FWHM_ARCSEC(ep) )

    CALL write_epoch_list_exec(CCID, MJD, BAND, VALUES,  & 
            MXCHAR_CCID, 2 )

    RETURN
  END SUBROUTINE WRITE_EPOCH_IGNORE_EXEC

! =======================================
    SUBROUTINE INIT_CUTNAMES ( IERR )
! 
! Creatd Jan 31, 2007 by R.Kessler
! Initialize strings for cut names.
! 
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IERR ! (O)  0=> OK
    INTEGER ibit, ifilt, ifilt_obs
    character cflt*1

! -------------- BEGIN -------------------

    IERR = 0  ! init output arg

    do ibit = 1, NCUTBIT
      cutvar_name(ibit)  = 'UNDEFINED:'
    enddo

! hard-wire names of cut variables (maxchar = 28)

    cutvar_name(CUTBIT_CID)          = 'CID:'
    cutvar_name(CUTBIT_SNTYPE)       = 'TYPE:'
    cutvar_name(CUTBIT_RA)           = 'RA:'
    cutvar_name(CUTBIT_DEC)          = 'DEC:'
    cutvar_name(CUTBIT_HOSTSEP)      = 'HOST-SN sep:'
    cutvar_name(CUTBIT_TRESTMIN)     = 'Trestmin:'
    cutvar_name(CUTBIT_TRESTMAX)     = 'Trestmax:'
    cutvar_name(CUTBIT_TRESTRANGE)   = 'TrestRange:'
    cutvar_name(CUTBIT_TGAPMAX)      = 'TGAPmax:'
    cutvar_name(CUTBIT_T0GAPMAX)     = 'T0GAPmax:'
    cutvar_name(CUTBIT_TobsMIN)      = 'TobsMin:'
    cutvar_name(CUTBIT_TobsMAX)      = 'TobsMax:'
    cutvar_name(CUTBIT_PEAKMJD)      = 'PEAKMJD:'
    cutvar_name(CUTBIT_NOBS_PREDETECT) = 'NOBS_PREDETECT:'
    cutvar_name(CUTBIT_Nepoch)       = 'Nepoch:'
    cutvar_name(CUTBIT_REDSHIFT)     = 'Redshift:'
    cutvar_name(CUTBIT_REDSHIFT_ERR) = 'Redshift-ERROR:'
    cutvar_name(CUTBIT_PSF)          = 'PSF:'
    cutvar_name(CUTBIT_ZP)           = 'ZP:'
    cutvar_name(CUTBIT_ZPERR)        = 'ZPERR:'
    cutvar_name(CUTBIT_PHOTPROB)     = 'PHOTPROB:'
    cutvar_name(CUTBIT_MWEBV)        = 'MWEBV:'
    cutvar_name(CUTBIT_NSEASON_ACTIVE) = 'NSEASON_ACTIVE:'
    cutvar_name(CUTBIT_REQEP)        = 'REQEP_CUTFLAG:'
    cutvar_name(CUTBIT_PRIVATE)      = 'PRIVATEVAR_CUTS:'
    cutvar_name(CUTBIT_SIMVAR)       = 'SIMVAR_CUTS:'
    cutvar_name(CUTBIT_TREST)        = 'Trest:'
    cutvar_name(CUTBIT_TOBS)         = 'Tobs:'
    cutvar_name(CUTBIT_ERRTEST)      = 'ERR(CALC)/ERR(TRUE):'
    cutvar_name(CUTBIT_SIM_PULL)     = '(F-Ftrue)/sigF:'
    cutvar_name(CUTBIT_TREST_TRUEFLUX2) = 'Trest_TRUEFLUX:'
    cutvar_name(CUTBIT_SNRMAX)       = 'SNRmax:'
    cutvar_name(CUTBIT_SNRMAX2)      = 'SNRmax2:'
    cutvar_name(CUTBIT_SNRSUM)       = 'SNRSUM:'
    cutvar_name(CUTBIT_NFILT_SNRMAX) = 'NFILT_SNRmax:'
    cutvar_name(CUTBIT_NFILT_SNRMAX2)= 'NFILT_SNRmax2:'
    cutvar_name(CUTBIT_NFILT_TRESTMIN) = 'NFILT_Trestmin:'
    cutvar_name(CUTBIT_NFILT_TRESTMAX) = 'NFILT_Trestmax:'
    cutvar_name(CUTBIT_NFILT_TREST2)   = 'NFILT_Trest2:'
    cutvar_name(CUTBIT_NFIELD)         = 'NFIELD:'
    cutvar_name(CUTBIT_SEARCH)         = 'SEARCHEFF_MASK:'

    DO ifilt      = 1, NFILT_SNRMAX
        ifilt_obs = IFILT_SNRMAX(ifilt)
        cflt      = FILTDEF_STRING(ifilt_obs:ifilt_obs)
        ibit = CUTBIT_OFFSET_SNRMAX + ifilt
        cutvar_name(ibit)  = 'SNRMAX-' // cflt // ':'
    ENDDO

    DO ifilt      = 1, NFILT_HOST_SBFLUX  ! Aug 2021
        ifilt_obs = ifilt_HOST_SBFLUX(ifilt)
        cflt      = FILTDEF_STRING(ifilt_obs:ifilt_obs)
        ibit      = CUTBIT_OFFSET_SBFLUX + ifilt
        cutvar_name(ibit)  = 'SBFLUX-' // cflt // ':'
    ENDDO

    RETURN
  END SUBROUTINE INIT_CUTNAMES 

! ======================================
    SUBROUTINE MWEBV_FLUXCOR()
! 
! Compute MW extinction in each defined passband,
! and correct SNLC_FLUXCAL[_ERRTOT]
! 
! 
! Mar 06, 2025: make explicit check on <lam> for each SURVEY_FILTER
! --------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT

    IMPLICIT NONE

! local var

    INTEGER IFILT, IFILT_OBS, NEWMJD, EPMIN, EPMAX, ep, OPT, i

    REAL*8  & 
         XTMW8,  XTMW8_PLUSERR  & 
        ,AVMW8,  AVMW8_PLUSERR  & 
        ,MWEBV8, MWEBV8_PLUSERR, PARLIST_MWCOLORLAW8(10)  & 
        ,LAM8, GALextinct, RV8

    REAL    XTMW_cor, arg
    CHARACTER CFILT*2, FNAM*14

    LOGICAL LDMP

! ------------ BEGIN -------------

    LDMP = .FALSE.

    FNAM = 'MWEBV_FLUXCOR'
    IF ( USESIM_TRUEFLUX .and. LSIM_SNANA ) THEN
       SNLC_MWEBV     = SIM_MWEBV
       SNLC_MWEBV_ERR = 0.0
    ENDIF

    if ( .NOT. EXIST_CALIB_FILE ) RETURN ! Oct 2024

! determin MW extinction in each passband

     MWEBV8          = DBLE(SNLC_MWEBV )
     MWEBV8_PLUSERR  = DBLE(SNLC_MWEBV+SNLC_MWEBV_ERR)

     RV8             = DBLE(RV_MWCOLORLAW)
     AVMW8           = RV8 * MWEBV8
     AVMW8_PLUSERR   = RV8 * MWEBV8_PLUSERR

     OPT = OPT_MWCOLORLAW

     DO i = 1, 10
       PARLIST_MWCOLORLAW8(i) = DBLE(PARLIST_MWCOLORLAW(i))
     ENDDO

     if (  LDMP ) THEN
        write(6,122) SNLC_CCID(1:ISNLC_LENCCID),  & 
             SURVEY_FILTERS(1:NFILTDEF_SURVEY), MWEBV8, OPT
122       format(T5,'xxx --- XTMW(',A,'-', A,') for MWEBV=',F5.3,  & 
               2x,'OPT=',I3 )
     endif


    DO ifilt = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        LAM8      = DBLE ( FILTOBS_LAMAVG(ifilt_obs) )
        CFILT     = FILTDEF_STRING(ifilt_obs:ifilt_obs)

        if ( LAM8 < 1.0 ) THEN
           write(c1err,601) cfilt, LAM8
601          format('SURVEY_FILTER ', A1,' has invalid <lam> = ',  & 
                  F8.2, ' A' )
           c2err = 'Check if ' // CFILT(1:1) //  & 
                 '-band is defined in kcor/calib file.'
           CALL MADABORT(FNAM, c1err, c2err )
        endif

        XTMW8     =  & 
                 GALextinct ( RV8, AVMW8, LAM8,  & 
                      OPT, PARLIST_MWCOLORLAW8, FNAM//char(0), 20 )
        XTMW8_PLUSERR  =  & 
                 GALextinct ( RV8, AVMW8_PLUSERR, LAM8,  & 
                      OPT, PARLIST_MWCOLORLAW8, FNAM//char(0), 20 )

        SNLC_MWXT_MAG(ifilt)    = SNGL(XTMW8)
        SNLC_MWXT_MAGERR(ifilt) = SNGL(XTMW8_PLUSERR-XTMW8)

        if ( LDMP ) then
          write(6,123) CFILT, LAM8, XTMW8
123         format(T5,'xxx ', A,' : <LAM>=',F7.0, 2x, 'XTMW=',F6.4)
          call flush(6)
        endif

! store extinction in flux-fraction units.
! FLUXFRAC = Flux(Extincted)/flux(top-of-Galaxy)
        arg = SNGL(-0.4 * XTMW8)
        SNLC_MWXT_FLUXFRAC(ifilt) = 10.0**arg  ! <= 1.0

    END DO  ! end of ifilt loop


! correct data fluxes if USE_MWCOR=T
    IF ( .NOT.  USE_MWCOR )  RETURN

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

        EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
        EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

        DO EP = EPMIN, EPMAX

          IFILT_OBS = ISNLC_IFILT_OBS(ep)
          IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
          XTMW_cor  = 1.0 - SNLC_MWXT_FLUXFRAC(ifilt)

          SNLC_FLUXCAL(ep) =  & 
            SNLC_FLUXCAL(ep) / XTMW_cor

          SNLC_FLUXCAL_ERRTOT(ep) =  & 
            SNLC_FLUXCAL_ERRTOT(ep) / XTMW_cor

        ENDDO

200   CONTINUE  ! NEWMJD

    RETURN
  END SUBROUTINE MWEBV_FLUXCOR

! ===========================================
    SUBROUTINE RDFILE_USERTAGS()
! 
! Created Aug 2011 by R.Kessler
! If namelist USERTAGS_FILE is set, then read integer tag
! for each SN specified. These tags are then included in
! analysis ntuples and fitres files. Motivation is to
! easily tag subsets such as SN confirmed from a particular
! telescope.  For SN not listed in the USERTAGS_FILE,
! the tag is automatically set to zero.
! 
! If an unknown CCID is found in the USERTAGS_FILE,
! it is ignored (i.e., no abort)
! 
! The format of the user-tag file is
! 
!  SN:  <CID1>  <TAG-VALUE1>
!  SN:  <CID2>  <TAG-VALUE2>
!  SN:  <CID3>  <TAG-VALUE3>
!  etc ...
! 
! --------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE USRTAGCM

    IMPLICIT NONE

    INTEGER LEN, ITAG, iwd, NWD, N
    CHARACTER  & 
         CWD*60  & 
        ,CKEY*(MXCHAR_CCID)  & 
        ,CCID*(MXCHAR_CCID)  & 
        ,NAME_ForC*(MXCHAR_FILENAME)  & 
        ,FNAM*30

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ------------- BEGIN ----------------

    N_USERTAGS    = 0

    IF ( USERTAGS_FILE .EQ. ' ' ) RETURN

    FNAM = 'RDFILE_USERTAGS'

    LEN = INDEX(USERTAGS_FILE,' ') - 1

! init all tags to zero
    DO itag = 1, MXUSERTAG
      USERTAG_VALUELIST(itag) = 0
      USERTAG_CCIDLIST(itag)  = ''
      USERTAG_USED(itag)      = 0
    ENDDO

    CALL PRBANNER("RDFILE_USERTAGS")
    print*,'  Read user tags to identify subsets. '

    NAME_ForC = USERTAGS_FILE(1:LEN)//char(0)
    NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,NAME_forC,  & 
                FNAM//char(0), LEN, 30 )

    DO 100 iwd = 1, NWD-1

       CALL get_PARSE_WORD_fortran(iwd+0,CWD,LEN)
       CKEY = CWD(1:MXCHAR_CCID)

       if ( CKEY .EQ. 'SN:' ) then

         N_USERTAGS = N_USERTAGS + 1
         N = N_USERTAGS

         if ( N .GT. MXUSERTAG ) then
           write(c1err,61) MXUSERTAG
61           format('Number of USER TAGS exceeds bound of ',I5)
           c2err = 'Check ' // USERTAGS_FILE(1:60)
           CALL MADABORT("RDFILE_USERTAGS", c1err, c2err )
         endif

         CALL get_PARSE_WORD_fortran(iwd+1,CWD,LEN)
         READ ( CWD, *) CCID

         CALL get_PARSE_WORD_fortran(iwd+2,CWD,LEN)
         READ ( CWD, *) ITAG

         USERTAG_CCIDLIST(N)  = CCID
         USERTAG_VALUELIST(N) = ITAG

       endif

100   CONTINUE

    RETURN
  END SUBROUTINE RDFILE_USERTAGS


! =========================================
    SUBROUTINE GET_USERTAG(CCID)

! 
! Load global USERTAG value for this SN
! 


    USE SNPAR
    USE USRTAGCM

    IMPLICIT NONE

    CHARACTER CCID*(*)  ! (I)

    INTEGER i, LCCID, LTMP
    CHARACTER CCID_TMP*(MXCHAR_CCID)

! -------------- BEGIN --------------

    USERTAG = -999
    IF ( N_USERTAGS .LE. 0 ) RETURN

    LCCID = INDEX(CCID,' ') - 1
    DO 100 i = 1, N_USERTAGS

       if ( USERTAG_USED(i) .NE. 0 ) GOTO 100

       CCID_TMP = USERTAG_CCIDLIST(i)
       LTMP     = INDEX(CCID_TMP,' ') - 1

       IF ( LTMP          .NE. LCCID             ) GOTO 100
       IF ( CCID(1:LCCID) .NE. CCID_TMP(1:LCCID) ) GOTO 100

       USERTAG = USERTAG_VALUELIST(i)
       USERTAG_USED(i) = 1
       RETURN
100   CONTINUE


    RETURN
  END SUBROUTINE GET_USERTAG

! ==================================
    SUBROUTINE MULTISEASON(IFLAG)

! Created Oct 2014
! Driver for quick analysis to check for multi-season variability.
! 
! May 1 2019: convert GET_MULTISEASON args to double.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IFLAG ! (I)

    INTEGER NEWMJD, EPMIN, EPMAX, EP, NOBS, iSeason
    INTEGER ITMP_REJECT(MXEPOCH)
    REAL*8  TMP8_MJD(MXEPOCH)
    REAL*8  TMP8_FLUX(MXEPOCH), TMP8_FLUXERR(MXEPOCH)
    REAL*8  OUT_CHI2RED(MXSEASON), OUT_AVGFLUX(MXSEASON)
    REAL*8  OUT_MJDMIN(MXSEASON),  OUT_MJDMAX(MXSEASON)
    REAL   PARLIST(8), CHI2RED, CUTVAL
    CHARACTER CCID_forC*(MXCHAR_CCID)
    LOGICAL USE_TGAP, USE_MASK

! -------------------- BEGIN --------------

    IF ( IFLAG == IFLAG_INI ) THEN

       USE_TGAP = MULTISEASON_TGAP < 1.0E8
       USE_MASK = MULTISEASON_OPTMASK > 0
       if ( .NOT. USE_TGAP ) THEN
          MULTISEASON_OPTMASK = 0
          RETURN
       endif
       if ( .NOT. USE_MASK ) MULTISEASON_OPTMASK = 1

! copy user-analysis parameters to list
       PARLIST(1) = MULTISEASON_OPTMASK
       PARLIST(2) = MULTISEASON_TGAP
       PARLIST(3) = MULTISEASON_NREJECT_OUTLIER

! call to C function in multiseason.c
       CALL INIT_MULTISEASON(PARLIST)

       RETURN
    ENDIF

    if ( MULTISEASON_OPTMASK == 0 ) return

! ------ if we get here then analyze this CID -------

! create local array of MJD, FLUX, FLUXERR that were selected by SNRECON
    NOBS = 0

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 201 EP = EPMIN, EPMAX
         if ( ISNLC_SNRECON_USE(ep) == 0 ) GOTO 201
         NOBS = NOBS + 1
         TMP8_MJD(NOBS)     = SNLC8_MJD(ep)
         TMP8_FLUX(NOBS)    = DBLE( SNLC_FLUXCAL(ep) )
         TMP8_FLUXERR(NOBS) = DBLE( SNLC_FLUXCAL_ERRTOT(ep) )
 201    CONTINUE
 200  CONTINUE


! --------------
    CCID_forC = SNLC_CCID(1:ISNLC_LENCCID) // char(0)
    NOBS      = ISNLC_NEWMJD_STORE

    CALL GET_MULTISEASON(  & 
            CCID_forC  & 
           ,NOBS, TMP8_MJD, TMP8_FLUX, TMP8_FLUXERR  &  ! (I)
           ,ITMP_REJECT          &  ! (O) mask of reject epochs
           ,NSEASON_TOT          &  ! (O) number of seasons found
           ,OUT_CHI2RED          &  ! (O) reduced chi2 per season
           ,OUT_MJDMIN           &  ! (O) min-MJD for each season
           ,OUT_MJDMAX           &  ! (O) max-MJD for each season
           ,OUT_AVGFLUX          &  ! (O) avg flux for each season
           )

! transfer to float globals
    DO iSeason = 1, NSEASON_TOT
       MULTISEASON_CHI2RED(iSeason) = SNGL(OUT_CHI2RED(iSeason))
       MULTISEASON_MJDMIN(iSeason)  = SNGL(OUT_MJDMIN(iSeason))
       MULTISEASON_MJDMAX(iSeason)  = SNGL(OUT_MJDMAX(iSeason))
       MULTISEASON_AVGFLUX(iSeason) = SNGL(OUT_AVGFLUX(iSeason))
    ENDDO

    CUTVAL = MULTISEASON_CHI2RED_ACTIVE
    DO iSeason = 1, NSEASON_TOT
      CHI2RED = MULTISEASON_CHI2RED(iSeason)
      IF ( CHI2RED > CUTVAL ) THEN
         NSEASON_ACTIVE = NSEASON_ACTIVE + 1
      ENDIF
    ENDDO

! ------------------------------------------
! transfer ITMP_REJECT flag to ISNLC_RECON_USE
! so that plots show rejected epochs

    NOBS = 0

    DO 2000 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 2001 EP = EPMIN, EPMAX
         if ( ISNLC_SNRECON_USE(ep) == 0 ) GOTO 2001
         NOBS = NOBS + 1
         if ( ITMP_REJECT(NOBS) == 1 ) then
            ISNLC_SNRECON_USE(ep) = 0
         endif
 2001   CONTINUE
 2000 CONTINUE

    RETURN
  END SUBROUTINE MULTISEASON


! ============================================
    LOGICAL FUNCTION SELECT_EARLYLC(EP)
! 
! Created March 2015 by R.Kessler
! Return TRUE if EARLYLC-requirements are satsified.
! Logical AND is applied to each selection defined in
! subroutine PARSE_EARLYLC_STRING().
! 
! Note that all epochs BEFORE EARLYLC requirements are accepted;
! epochs are rejected after all of the EARLYLC obs or nights
! are found.
! 
! 9/15/2017: implement SNR_START and PHOTMASK_START (see manual)
! 
! -------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE EARLYCOM
    USE FILTCOM

    IMPLICIT NONE

    INTEGER EP    ! (I) epoch to analyze

! local var

    REAL SNR, PROB, MJD
    INTEGER IFILT_OBS, MASK, OVPMASK, OVPMASK_START, NMJD_DIF
    LOGICAL LCUTS, LSNR, LPHOTPROB, LPHOTMASK, LFILT
    LOGICAL LDMP, SAMENIGHT, LNOSEL
    CHARACTER CFILT*2

! --------------- BEGIN --------------

    SELECT_EARLYLC = .TRUE.
    IF ( EARLYLC_STRING .EQ. '' ) RETURN

!       print*,' xxx HELOO from SELECT_EARLYLC: ep = ', ep

! strip off info for this epoch.
    IFILT_OBS = ISNLC_IFILT_OBS(ep)
    CFILT     = FILTDEF_STRING(ifilt_obs:ifilt_obs)
    PROB      = SNLC_PHOTPROB(ep)
    MASK      = ISNLC_PHOTFLAG(ep)
    MJD       = SNGL( SNLC8_MJD(ep) )

    if ( SNLC_FLUXCAL_ERRTOT(ep) > 0 ) then
       SNR  = SNLC_FLUXCAL(ep) / SNLC_FLUXCAL_ERRTOT(ep)
       if ( SNR > SNR_START_EARLYLC ) then
          NSNR_START_EARLYLC = NSNR_START_EARLYLC + 1
       endif
    else
       SNR = 0.0
    endif

    IF ( PHOTMASK_EARLYLC .NE. 0 ) THEN
      OVPMASK = IAND(PHOTMASK_EARLYLC,MASK)
    ELSE
      OVPMASK = 999  ! any number > 0
    ENDIF

! check PHOTMASK_START required to start counting epochs (9/15/2017)
! --> ignore epochs before PHOTMASK_START is satisfied.
    IF ( PHOTMASK_START_EARLYLC > 0 ) THEN
      OVPMASK_START = IAND(PHOTMASK_START_EARLYLC,MASK)
    ELSE
      OVPMASK_START = 999  ! any number > 0
    ENDIF

    IF ( OVPMASK_START > 0 ) THEN
       NPHOTMASK_START_EARLYLC = NPHOTMASK_START_EARLYLC + 1
    ENDIF

! set cut-logicals.
    LFILT      = (INDEX(FILTERS_EARLYLC,CFILT(1:1)) > 0)

    LSNR       = (SNR .GE. SNRMIN_EARLYLC) .AND.  & 
                   (NSNR_START_EARLYLC>0)

    LPHOTPROB  = PROB .GE. PHOTPROBMIN_EARLYLC

    LPHOTMASK  = ( OVPMASK > 0) .AND.  & 
                   ( NPHOTMASK_START_EARLYLC > 0)

    SAMENIGHT  = ( (MJD-MJDLAST_EARLYLC) < DT_SAMENIGHT )

    LCUTS = (LFILT .and. LSNR .and. LPHOTPROB .and. LPHOTMASK )

! apply all cuts to increment number of early obs.
    IF ( LCUTS ) THEN
       NOBS_EARLYLC = NOBS_EARLYLC + 1
       if ( .NOT. SAMENIGHT ) then
          NNIGHT_EARLYLC = NNIGHT_EARLYLC + 1
          MJDLAST_EARLYLC = MJD
       endif
    ENDIF

! ---------------------------------------------------
! check if we reached max number of OBS or NIGHTS

    IF( NOBS_EARLYLC > MAXOBS_EARLYLC ) THEN
       SELECT_EARLYLC = .FALSE.
    ENDIF

    IF ( NNIGHT_EARLYLC > MAXNIGHT_EARLYLC ) THEN
       SELECT_EARLYLC = .FALSE.
    ENDIF

    IF ( SELECT_EARLYLC ) MJDLAST_SELECT = MJD


! check NDAYADD option to add number of days instead of
! number of obs or nights.

    LNOSEL = (.NOT. SELECT_EARLYLC) .and. (NDAYADD_EARLYLC>0)
    NMJD_DIF = -9
    if ( LNOSEL .and. MJDLAST_SELECT > 10000 ) then
       NMJD_DIF = int(MJD - MJDLAST_SELECT+0.5)
       SELECT_EARLYLC = (NMJD_DIF < NDAYADD_EARLYLC)
    endif

! xxxxxxxxxxxxxxx
!      print*,' xxx LNOSEL=',LNOSEL,
!     &   '  MJD,MJDLAST_SELECT=', MJD,MJDLAST_SELECT
!      print*,' xxx --> NMJD_DIF=', NMJD_DIF,
!     &     '   SELECT=', SELECT_EARLYLC
! xxxxxxxxxxxxxxx

    LDMP = .FALSE.  ! (PROB < 1.0E9)
    if ( LDMP ) THEN
       print*,' xxx ------------------------------------- '
       write(6,666) EP, CFILT, LFILT,  & 
                 MJD, PROB, NOBS_EARLYLC, NNIGHT_EARLYLC
666      format(T2, 'xxx EP=',I3'-',A1, '(',L1,')', 2x,  & 
               'MJD=',F9.3, 2x, 'PROB=',F5.2, 2x,  & 
               'N[OBS,NITE]=',2I3 )

       print*,' xxx SNR, SNRMIN_EARLYLC = ', SNR, SNRMIN_EARLYLC
       print*,' xxx MAXOBS_EARLYLC, MAXNIGHT_EARLYLC = ',  & 
                      MAXOBS_EARLYLC, MAXNIGHT_EARLYLC

       print*,' xxx NPHOTMASK_START_EARLYLC = ',  & 
                      NPHOTMASK_START_EARLYLC

       write(6,667) LFILT, LSNR, LPHOTPROB, LPHOTMASK
667      format(T2, 'xxx LCUT(FILT,SNR,PHOTPROB,PHOTMASK) = ',4L3)
       print*,' xxx LCUTS, SELECT = ', LCUTS, SELECT_EARLYLC
       call flush(6)
    endif

    RETURN
  END FUNCTION SELECT_EARLYLC

! ========================================
    SUBROUTINE PARSE_EARLYLC_STRING()

! Mar 2015 -
!  parse &SNLCINP input EARLYLC_STRING and set global
!  variables in common EARLYCOM.
!  Also print summary to stdout.
! 
! Valid keys contained in the EARLYLC_STRING are
! (with example values):
!   MAXOBS:            4        ! default = 999
!   MAXNIGHT:          4        ! default = 999
!   FILTERS:         riz      ! default = all bands
!   SNRMIN:          5.0      ! default = -999
!   PHOTPROBMIN:     0.5      ! default = -999
!   PHOTMASK        4096      ! default =  0
!   NDAYADD           10      ! default =  0
!   PHOTMASK_START  4096      ! default =  0    (disabled)
!   SNR_START          5      ! default =  999  (disabled)
! 
! 9/15/2017: check for PHOTMASK_START
! 

! -------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE EARLYCOM

    IMPLICIT NONE

    INTEGER iwd, NWD, NOTUSED, L0, L1, MSKOPT
    CHARACTER  & 
         cwd0*(MXCHAR_FILEWORD),  & 
         cwd1*(MXCHAR_FILEWORD), FNAM*22
    LOGICAL USEWD(200)

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! -------------- BEGIN -----------------

    FNAM = 'PARSE_EARLYLC_STRING'
    MSKOPT = MSKOPT_PARSE_WORDS_STRING

! set defaults

    NDAYADD_EARLYLC        =  0
    MAXOBS_EARLYLC         =  999
    MAXNIGHT_EARLYLC       =  999
    SNRMIN_EARLYLC         = -999.
    PHOTPROBMIN_EARLYLC    = -999.
    PHOTMASK_EARLYLC       =  0
    PHOTMASK_START_EARLYLC =  0
    NPHOTMASK_START_EARLYLC=  0
    SNR_START_EARLYLC      =  -9999.
    NSNR_START_EARLYLC     =  0
    FILTERS_EARLYLC        = SURVEY_FILTERS

    IF ( EARLYLC_STRING .EQ. '' ) RETURN

    NWD = STORE_PARSE_WORDS(MSKOPT, EARLYLC_STRING//char(0),  & 
               FNAM//char(0), 100, 30)

    DO iwd = 1, NWD
       USEWD(iwd) = .FALSE.
    ENDDO

    DO 10 iwd = 1, NWD

        cwd0 = '' ; cwd1 = ''

        CALL get_PARSE_WORD_fortran(iwd+0, cwd0, L0)
        if ( iwd < NWD ) then
           CALL get_PARSE_WORD_fortran(iwd+1, cwd1, L1)
        endif

        if ( cwd0(1:6) .EQ. 'MAXOBS' ) then
           read(cwd1,* ) MAXOBS_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:7) .EQ. 'NDAYADD' ) then
           read(cwd1,* ) NDAYADD_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:8) .EQ. 'MAXNIGHT' ) then
           read(cwd1,* ) MAXNIGHT_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:7) .EQ. 'FILTERS' ) then
           FILTERS_EARLYLC = cwd1
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:6) .EQ. 'SNRMIN' ) then
           read(cwd1,* ) SNRMIN_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:11) .EQ. 'PHOTPROBMIN' ) then
           read(cwd1,* ) PHOTPROBMIN_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

!   be careful with two keys that have same 'PHOTMASK' in string
        if ( cwd0(1:14) .EQ. 'PHOTMASK_START' ) then
           read(cwd1,* ) PHOTMASK_START_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        else if ( cwd0(1:8) .EQ. 'PHOTMASK' ) then
           read(cwd1,* ) PHOTMASK_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

        if ( cwd0(1:9) .EQ. 'SNR_START' ) then
           read(cwd1,* ) SNR_START_EARLYLC
           USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
        endif

10    CONTINUE

    GLOBAL_BANNER = 'Init Selection of Early Part of Light Curve'
    CALL PRBANNER (GLOBAL_BANNER(1:60) )

    write(6,50) MAXOBS_EARLYLC, MAXNIGHT_EARLYLC
50    format(T5,'Select first ',I2,' obs or ',I2,' nights satisfying:')
    print*,'        FILTER  in   ', FILTERS_EARLYLC(1:40)
    print*,'   and  PHOTPROB > ', PHOTPROBMIN_EARLYLC
    print*,'   and  SNR      > ', SNRMIN_EARLYLC
    print*,'   and  NDAYADD  = ', NDAYADD_EARLYLC

    if ( PHOTMASK_EARLYLC > 0 ) then
       print*,'   and  PHOTFLAG contains ', PHOTMASK_EARLYLC
    endif

    if ( PHOTMASK_START_EARLYLC > 0 ) then
       print*,'   After epoch where PHOTFLAG contains ',  & 
               PHOTMASK_START_EARLYLC
    endif

    if ( SNR_START_EARLYLC > -900. ) then
       print*,'   After epoch where SNR > ', SNR_START_EARLYLC
    endif

    print*,' '
    call flush(6)

! ------------------------------------
! if any words are not used --> ABORT !!!

    NOTUSED =  0
    DO 66 iwd = 1, NWD
      if ( .NOT. USEWD(iwd) ) then
          NOTUSED = NOTUSED + 1
          CALL get_PARSE_WORD_fortran(iwd, cwd0, L0)
          write(6,666) cwd0
666         format('ERROR: EARLYLC_STRING contains invalid ', A20)
      endif
66    CONTINUE

    IF ( NOTUSED > 0 ) THEN
      write(c1err,61) NOTUSED
      write(c2err,62) EARLYLC_STRING
61      format('Found ',I2,' undefined words in')
62      format('EARLYLC_STRING="', A48, '"' )
      CALL MADABORT(FNAM,C1ERR,C2ERR)
    ENDIF

    RETURN
  END SUBROUTINE PARSE_EARLYLC_STRING

! ========================================================
    INTEGER FUNCTION ISTAT_REQUIRE_EPOCHS()

! Created Sep 2017
! Return 1 (TRUE) if epochs satify user input REQUIRE_EPOCHS_STRING.
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE REQEPCOM

    IMPLICIT NONE

    INTEGER  & 
          NEWMJD, EPMIN, EPMAX, EP, i, NOBS(3)

    REAL*8  & 
         MJD_FIRST, MJD_LAST, MJD, MJD_DIF  & 
        ,MJD_WIN_BEFORE(2), MJD_WIN_AFTER(2)  & 
        ,SNR, FLUX, FLUXERR, z1

    CHARACTER FNAM*24
    LOGICAL  LPASS, LTRANGE(3), LDMP

! ------------------- BEGIN ----------------

    FNAM       = 'ISTAT_REQUIRE_EPOCHS'

    ISTAT_REQUIRE_EPOCHS     = 1
    NDAYS_ABOVE_SNRMIN_REQEP = 0.0
    IF ( NFILT_REQEP .EQ. 0 ) RETURN

    MJD_FIRST  = -9.0
    MJD_LAST   = -9.0

    do i = 1, 3
      LTRANGE(i) = .FALSE.
      NOBS(i)    = 0
    ENDDO

    IF ( ISFRAME_OBS_REQEP ) then
       z1   = 1.0
    else
       z1  = 1.0 + SNLC_REDSHIFT
       if ( SNLC_REDSHIFT < 0.0 ) then
          c1err = 'Negative redshift !'
          c2err = 'Cannot compute TREST'
          CALL MADABORT(FNAM, c1err, c2err )
       endif
    endif

! - - - - - - -  - -

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

    DO 201 EP = EPMIN, EPMAX

         MJD      = SNLC8_MJD(EP)
         FLUX     = DBLE( SNLC_FLUXCAL(EP) )
         FLUXERR  = DBLE( SNLC_FLUXCAL_ERRTOT(EP) )
         SNR      = 0.0
         IF ( FLUXERR > 0.00001 ) THEN
            SNR = FLUX / FLUXERR
         ENDIF

         IF ( SNR > SNRMIN_REQEP ) THEN
            if ( MJD_FIRST < 0.0 ) MJD_FIRST = MJD
            MJD_LAST = MJD
            NOBS(2) = NOBS(2) + 1
         ENDIF
201   CONTINUE
200   CONTINUE


! -----------------------------------
! check range of epochs above SNR
    MJD_DIF = MJD_LAST - MJD_FIRST
    NDAYS_ABOVE_SNRMIN_REQEP = SNGL( MJD_DIF/z1 )
    if ( NDAYS_ABOVE_SNRMIN_REQEP > TRANGE_REQEP(2)  & 
              .OR. NOBS(2)==0 ) THEN
       ISTAT_REQUIRE_EPOCHS = 0
       GOTO 888
    endif

! We have a short-enough transient.
! Check that there there is at least 1 obs before and after
! the live range, so that a longer-duration event could
! have been observed.

    MJD_WIN_BEFORE(1) = MJD_FIRST - TRANGE_REQEP(1)*z1
    MJD_WIN_BEFORE(2) = MJD_FIRST - 0.001

    MJD_WIN_AFTER(1) = MJD_LAST + 0.001
    MJD_WIN_AFTER(2) = MJD_LAST + TRANGE_REQEP(3)*z1

    DO 400 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
    DO 401 EP = EPMIN, EPMAX

         MJD      = SNLC8_MJD(EP)

         if ( MJD > MJD_WIN_BEFORE(1) .and.  & 
                MJD < MJD_WIN_BEFORE(2) ) then
           NOBS(1) = NOBS(1) + 1
         endif

         if ( MJD > MJD_WIN_AFTER(1) .and.  & 
                MJD < MJD_WIN_AFTER(2) ) then
           NOBS(3) = NOBS(3) + 1
         endif

401   CONTINUE
400   CONTINUE

! ------------------------------------------------------

    DO i = 1, 3
      if ( NOBS(i)>0 .OR. TRANGE_REQEP(i) < 1.0E-9 ) THEN
         LTRANGE(i) = .TRUE.
      endif
    ENDDO

    LPASS = ( LTRANGE(1) .and. LTRANGE(2) .and. LTRANGE(3) )
    IF ( LPASS ) THEN
       ISTAT_REQUIRE_EPOCHS = 1
    ELSE
       ISTAT_REQUIRE_EPOCHS = 0
    ENDIF

! - - - - - - - - - - -
888   CONTINUE

    LDMP = ( SNLC_CID < -5 )
    IF ( LDMP ) THEN
      print*,' '
      print*,' xxx ---------- DUMP ISTAT_REQUIRE_EPOCHS ------ '
      print*,' xxx CID = ', SNLC_CID
      print*,' xxx TRANGE  = ', TRANGE_REQEP
      print*,' xxx LTRANGE = ', LTRANGE
      print*,' xxx NOBS    = ', NOBS
      print*,' xxx MJD_FIRST/LAST = ',  & 
            sngl(MJD_FIRST), '/', sngl(MJD_LAST)
      print*,' xxx ISTAT_REQUIRE_EPOCHS=', ISTAT_REQUIRE_EPOCHS
      print*,' '
    ENDIF

    RETURN
  END FUNCTION ISTAT_REQUIRE_EPOCHS


! ========================================================
    SUBROUTINE PARSE_REQUIRE_EPOCHS_STRING()

! Created Sep 2017
! Parse &SNLCINP input REQUIRE_EPOCHS_STRING.
! Main goal is to select fast transients and reject the
! larger source of Supernova.
! 
! Example:
!   REQUIRE_EPOCHS_STRING = 'FILTERS riz  TOBS_RANGES  14 20 25 SNRMIN 5'
!   REQUIRE_EPOCHS_STRING = 'FILTERS riz  TREST_RANGES 14 20 25 SNRMIN 5'
! 
! Legacy Example:
!    REQUIRE_EPOCHS_STRING = 'gri 14 20 21'
! --> For the gri passbands:
!  + Require (MJD_max - MJD_min)_DETECT < 20 days = [TRANGE(2)]
!  + require NOBS>=1 within 14 days before MJD_min  [TRANGE(1)]
!  + require NOBS>=1 within 21 days after  MJD_max  [TRANGE(3)]
! 
!  For TRANGE(2), a "detection" is defined by SNR > CUTWIN_SNRMAX(1).
!  The NOBS>=1 requirements ensure observations before and after
!  the detect-range to avoid edge effects such as a falling
!  light curve at the start of a season where we can't see the
!  entire LC. Such edge cases can be rejected.
! 
!  Setting TRANGE(1)=0 or TRANGE(3)=0 will disable the
!  corresponding requirement.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE REQEPCOM

    IMPLICIT NONE

    INTEGER ifilt, IFILT_OBS, iwd, NWD, i, L0, L1, MSKOPT
    character FNAM*28, CFILT*2, cwd0*80, cwd1*80
    LOGICAL LEGACY, USEWD(40)

    INTEGER FILTINDX

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! --------------- BEGIN ----------------

    FNAM = 'PARSE_REQUIRE_EPOCHS_STRING'
    MSKOPT = MSKOPT_PARSE_WORDS_STRING

! set defaults
    ISFRAME_REST_REQEP = .FALSE.
    ISFRAME_OBS_REQEP  = .FALSE.
    DO i = 1, 3
      TRANGE_REQEP(i)  = 0.0
    ENDDO
    NFILT_REQEP     = 0
    FILTLIST_REQEP  = ''
    SNRMIN_REQEP    = 0.0   ! defines detection

    IF ( REQUIRE_EPOCHS_STRING .EQ. '' ) RETURN

    NWD = STORE_PARSE_WORDS(MSKOPT,REQUIRE_EPOCHS_STRING//char(0),  & 
                FNAM//char(0), 100, 22 )

    LEGACY = ( NWD .EQ. 4 )

    IF ( LEGACY ) THEN
      ISFRAME_OBS_REQEP = .TRUE.
      call get_PARSE_WORD_fortran(1, cwd0, L0)
      read( CWD0,* ) FILTLIST_REQEP

      call get_PARSE_WORD_fortran(2, cwd0, L0)
      read( CWD0,* ) TRANGE_REQEP(1)

      call get_PARSE_WORD_fortran(3, cwd0, L0)
      read( CWD0,* ) TRANGE_REQEP(2)

      call get_PARSE_WORD_fortran(4, cwd0, L0)
      read( CWD0,* ) TRANGE_REQEP(3)

      SNRMIN_REQEP = CUTWIN_SNRMAX(1)
    ELSE
       DO iwd = 1, NWD
          USEWD(iwd) = .FALSE.
       ENDDO

       DO 44 iwd = 1, NWD-1

         call get_PARSE_WORD_fortran(iwd+0, cwd0, L0)
         call get_PARSE_WORD_fortran(iwd+1, cwd1, L1)

         if ( cwd0(1:7) .EQ. 'FILTERS' ) then
            read(cwd1,* ) FILTLIST_REQEP
            USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
         else if ( cwd0(1:6) .EQ. 'SNRMIN' ) then
            read(cwd1,* ) SNRMIN_REQEP
            USEWD(iwd) = .TRUE.;  USEWD(iwd+1) = .TRUE.
         else if ( cwd0(1:11) .EQ. 'TOBS_RANGES' ) then
            USEWD(iwd) = .TRUE.;
            ISFRAME_OBS_REQEP = .TRUE.
            do i = 1, 3
              call get_PARSE_WORD_fortran(iwd+i, cwd1, L1)
              read(cwd1,* ) TRANGE_REQEP(i); USEWD(iwd+i)=.TRUE.
            enddo
         else if ( cwd0(1:12) .EQ. 'TREST_RANGES' ) then
            USEWD(iwd) = .TRUE.;
            ISFRAME_REST_REQEP = .TRUE.
            do i = 1, 3
              call get_PARSE_WORD_fortran(iwd+i, cwd1, L1)
              read(cwd1,* ) TRANGE_REQEP(i); USEWD(iwd+i)=.TRUE.
            enddo
         endif

44       CONTINUE
    ENDIF

! ------------------------------------------
! store absolute filter indices in sparse list
    NFILT_REQEP = INDEX(FILTLIST_REQEP,' ') - 1
    DO 10 ifilt=1, NFILT_REQEP
      CFILT = FILTLIST_REQEP(ifilt:ifilt)
      IFILT_OBS = FILTINDX(cfilt // ' ')
      IFILTLIST_REQEP(IFILT) = IFILT_OBS
10    CONTINUE

    CALL PRBANNER(FNAM)
    print*,'   Require Epochs for filters : ',  & 
             FILTLIST_REQEP(1:NFILT_REQEP)

    write(6,41) 'Detection Requires SNR  > ', SNRMIN_REQEP
    write(6,42) '(Tlast - Tfirst)_Detect < ', TRANGE_REQEP(2)
    write(6,42) 'Pre-detect  veto range  = ', TRANGE_REQEP(1)
    write(6,42) 'Post-detect veto range  = ', TRANGE_REQEP(3)
    print*,' ISFRAME[OBS,REST] = ',  & 
           ISFRAME_OBS_REQEP, ISFRAME_REST_REQEP
41    format(T8,A, F6.1 )
42    format(T8,A, F6.1, ' days' )

    RETURN
  END SUBROUTINE PARSE_REQUIRE_EPOCHS_STRING



! ========================================================
    SUBROUTINE INIT_MAGCOR(MAGCOR_INFILE)

! Created Dec 7 2016 by R.Kessler
! Read optional list of mag-correction vs. epoch from ASCII file.
! The file syntax is a FITRES-formatted file containing
!    VARNAMES:  ROW  MAGCOR EXTRA1 EXTRA2 ...
!    SN: [CID]-[MJD]-[BAND]  [MAGCOR] [EXTRA1] [EXTRA2]
!    etc ..
! and see snana manual section  "Epoch-Dependent Mag Corrections".
! Here we read only ROW and MAGCOR, and ignore all other variables.
! ROW must be first column, but MAGCOR can be in any column.
! 
! If there is a minus sign in front of the file name,
!     MAGCOR_INFILE = '-/BLA/BLA/MAGCOR.DAT'
! then subtract MAGCOR instead of adding.
! 
!  NVAR: 2
!  VARNAMES:  ROW  MAGCOR
!  ROW:  [CID]-[MJD]-[BAND]  [MAGCOR]
!  ROW:  [CID]-[MJD]-[BAND]  [MAGCOR]
!   etc ...
! 
! Additional columns are allowed, but will be ignored.
! 
! Nov 13 2018: pass MAGCOR_INFILE as input argument.
! 
! ----------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER MAGCOR_INFILE*(MXCHAR_FILENAME)  ! (I)

    INTEGER LEN1, LEN2, LEN3, OPTMASK_STORE, ISTAT
    CHARACTER  & 
           cFILE*(MXCHAR_FILENAME)  & 
          ,cTABLE*20, cVARLIST*20, FNAM*12

! functions
    INTEGER  SNTABLE_AUTOSTORE_INIT, EXIST_VARNAME_AUTOSTORE
    EXTERNAL SNTABLE_AUTOSTORE_INIT, EXIST_VARNAME_AUTOSTORE
    INTEGER ISTAT_MAGCOR_ADD

! ------------ BEGIN -------------

    NSTORE_MAGCOR = 0
    NUSE_MAGCOR   = 0
    IF ( MAGCOR_INFILE .EQ. ' '    ) RETURN
    IF ( MAGCOR_INFILE .EQ. 'NULL' ) RETURN
    IF ( MAGCOR_INFILE .EQ. 'NONE' ) RETURN

    SIGN_MAGCOR = +1  ! default is to add
    IF ( MAGCOR_INFILE(1:1) .EQ. '-' ) THEN
       SIGN_MAGCOR = -1  ! subtract instead
       MAGCOR_INFILE = MAGCOR_INFILE(2:)
    ENDIF

    FNAM = 'INIT_MAGCOR'
    CALL PRBANNER('INIT_MAGCOR: read mag-correction vs. epoch')

    CALL ENVreplace(MAGCOR_INFILE)
    LEN1   = INDEX(MAGCOR_INFILE,' ') - 1
    cFILE  = MAGCOR_INFILE(1:LEN1) // char(0)

! check header to see if if MAGCOR has been added in version
! ISTAT=1 if MAGCOR alread added ; ISTAT=0 if not already added.
    ISTAT = ISTAT_MAGCOR_ADD(MAGCOR_INFILE,VERSION_PHOTOMETRY(1))

! avoid double-counting MAGCOR
    IF ( ISTAT==0 .and. SIGN_MAGCOR < 0 ) THEN
       C1ERR = 'MAGCOR not in data fluxes --> '
       C2ERR = 'Will not subtract MAGCOR'
       CALL MADABORT(FNAM,C1ERR,C2ERR)
    ENDIF
    IF ( ISTAT==1 .and. SIGN_MAGCOR > 0 ) THEN
       C1ERR = 'MAGCOR already in data fluxes --> '
       C2ERR = 'Will not double-count MAGCOR.'
       CALL MADABORT(FNAM,C1ERR,C2ERR)
    ENDIF

    cTABLE = 'MAGCOR' // char(0)  ! any table name will work for ASCII
    LEN2   = INDEX(cTABLE,' ') - 1

    cVARLIST = 'MAGCOR' // char(0)
    LEN3     = INDEX(cVARLIST,' ') - 1

    OPTMASK_STORE = 5  ! 1->5 on May 11 2017
    NSTORE_MAGCOR =  & 
             SNTABLE_AUTOSTORE_INIT(cFILE, cTABLE, cVARLIST,  & 
                     OPTMASK_STORE,  LEN1, LEN2, LEN3 )

    write(6,20) NSTORE_MAGCOR
 20   format(/, T4,'Stored ', I6, ' epochs of MAG corrections. ')
    print*,'    Sign of MAGCOR : ', SIGN_MAGCOR
    print*,' '  ;  CALL FLUSH(6)

    RETURN
  END SUBROUTINE INIT_MAGCOR


! ========================================
    INTEGER FUNCTION ISTAT_MAGCOR_ADD(MAGCOR_FILE,VERSION)

! Check header of MAGCOR_FILE to see if VERSION has MAGCOR
! already added.  CHeck MAGCOR_FILE for
! 
!   VERSION_ADD:  <VERSION>
! 
! If MAGCOR already added in data files, return 1.
! If MAGCOR not added, return 0
! 
! This function is used to help avoid double-counting MAGCOR,
! but the MAGCOR-add process (outside SNANA) must add  the
! VERSION_ADD key at the top of the file.
! 
! May 22 2018: skip lines with comment field


    USE SNPAR

    IMPLICIT NONE

! inputs
    CHARACTER MAGCOR_FILE*(*), VERSION*(*)

! local
    INTEGER FOUND_ROW, LUN, NWD, LEN, LENVER, NLINE, MSKOPT
    LOGICAL SAME_LEN, SAME_VER, KEYMATCH
    CHARACTER LINE*200, cwd*80, VERSION_ADD*80, FNAM*22
    CHARACTER C1ERR*72, C2ERR*72

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS
! --------------- BEGIN ---------------

    ISTAT_MAGCOR_ADD = 0 ! init function value
    FNAM   = 'ISTAT_MAGCOR_ADD'
    MSKOPT = MSKOPT_PARSE_WORDS_STRING
    LENVER = INDEX(VERSION,' ' ) - 1

    NLINE=0; FOUND_ROW=0 ; LUN = 24
    OPEN(   UNIT   = LUN  & 
            , FILE   = MAGCOR_FILE  & 
            , STATUS = 'OLD'  & 
                 )


! keep reading until first ROW key is found

    DO 200 WHILE ( FOUND_ROW == 0 )
       READ(LUN,100) LINE
 100     format(A160)

       IF ( LINE(1:1) .EQ. '#' ) GOTO 200

       NLINE = NLINE + 1

       if ( NLINE > 50 ) then
          C1ERR = 'Could not find ROW key after 50 lines.'
          C2ERR = 'Check ' // MAGCOR_FILE
          CALL MADABORT(FNAM,C1ERR,C2ERR)
       endif

       NWD = STORE_PARSE_WORDS(MSKOPT, LINE//char(0),  & 
                  FNAM//char(0), 200, 22)

       if ( NWD < 2 ) goto 200

       CALL get_PARSE_WORD_fortran(1,cwd,LEN)

       IF ( cwd(1:4) .EQ. 'ROW:' ) goto 500

       KEYMATCH = cwd(1:16) .EQ. 'VERSION_ADD:'
       IF ( KEYMATCH ) THEN
          CALL get_PARSE_WORD_fortran(2,VERSION_ADD,LEN)
          SAME_LEN = ( LEN .EQ. LENVER )
          SAME_VER = ( VERSION(1:LEN) .EQ. VERSION_ADD(1:LEN) )
          IF ( SAME_LEN .and. SAME_VER ) THEN
            ISTAT_MAGCOR_ADD = 1
            GOTO 500
          ENDIF
       ENDIF
200   CONTINUE

500   CONTINUE
    CLOSE(UNIT = LUN)

    RETURN
  END FUNCTION ISTAT_MAGCOR_ADD

! ========================================================
    SUBROUTINE EXEC_MAGCOR(ep)

! Dec 2016
! Apply MAGCOR read in INIT_MAGCOR. The corrections are
! applied to the FLUXCAL.
! Modify SNLC_FLUXCAL(ep) and SNLC_MAG(ep)
! 
! May 1 2017: abort on crazy MAGCOR
! ---------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER ep ! epoch index

! local variables

    INTEGER IFILT_OBS, ISTAT, L, L2
    REAL*8 MJD, DVAL
    REAL*4 MAGCOR, FCOR
    CHARACTER cVARNAME*20, BAND*2, cDUM*20
    CHARACTER STR_EPID1*60, STR_EPID2*60

    REAL, PARAMETER :: MAGCOR_CRAZY = 0.2

! function
    EXTERNAL SNTABLE_AUTOSTORE_READ
! --------------- BEGIN -----------------

    IF ( NSTORE_MAGCOR <= 0 ) RETURN

    IFILT_OBS = ISNLC_IFILT_OBS(ep)
    BAND      = FILTDEF_STRING(ifilt_obs:ifilt_obs)
    MJD       = SNLC8_MJD(ep)

! construct unique epoch-identifier by gluing SNID+MJD+BAND
! into a single string.

    L  = INDEX(SNLC_CCID,' ') - 1
    WRITE(STR_EPID1,40) SNLC_CCID(1:L), MJD, BAND, char(0)
40    FORMAT(A,'-',F9.3,'-', A1, A )

    cVARNAME = 'MAGCOR' // char(0)

    CALL SNTABLE_AUTOSTORE_READ(STR_EPID1, cVARNAME, ISTAT,  & 
                       DVAL,cDUM, 60,10,10 )

! if no MAGCOR, then try again with IAUC name
    IF ( ISTAT .NE. 0 ) THEN
      L2  = INDEX(SNLC_NAME_IAUC,' ') - 1
      WRITE(STR_EPID2,40) SNLC_NAME_IAUC(1:L2), MJD, BAND, char(0)
      CALL SNTABLE_AUTOSTORE_READ(STR_EPID2, cVARNAME, ISTAT,  & 
                       DVAL,cDUM, 60,10,10 )
    ENDIF

    IF ( ISTAT .EQ. 0 ) then
      NUSE_MAGCOR  = NUSE_MAGCOR + 1
      MAGCOR = sngl(DVAL) * SIGN_MAGCOR

! trap crazy MAGCOR
      if ( abs(MAGCOR) > MAGCOR_CRAZY ) then
         write(C1ERR,61) MAGCOR
61         format('Crazy MAGCOR = ', G12.4 )
         C2ERR = 'Check ' // STR_EPID1(1:L+12)
         CALL MADABORT("EXEC_MAGCOR",C1ERR,C2ERR)
      endif

      FCOR   = 10**(-0.4*MAGCOR)

      SNLC_FLUXCAL(ep) = SNLC_FLUXCAL(ep) * FCOR
      SNLC_MAG(ep)     = SNLC_MAG(ep) + MAGCOR

!        write(6,66) DVAL, STR_EPID1, ISTAT
! 6      format(' xxx MAGCOR=',F6.3,' for ', A20,'  ISTAT=',I3)
!        call flush(6)
    ENDIF

    CALL SETMASK_FLUXCOR_SNANA(MASK_FLUXCOR_SNANA)

    RETURN
  END SUBROUTINE EXEC_MAGCOR

! ========================================================
    SUBROUTINE END_MAGCOR()

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! --------------- BEGIN --------------

    IF ( NSTORE_MAGCOR == 0 ) RETURN

    print*,' '
    print*,' MAGCOR SUMMARY:'
    write(6,20) NUSE_MAGCOR, NSTORE_MAGCOR
20    format(T8,'Used ',I6,' MAGCOR epochs from ',  & 
                I6,' read from file.' )

    RETURN
  END SUBROUTINE END_MAGCOR


! ========================================================
    SUBROUTINE RD_VPEC_FILE()
! ----------------------------

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! ------------- BEGIN ------------
    IF ( VPEC_FILE .EQ. '' ) RETURN
    C1ERR = '&SNLCIP VPEC_FILE input is obsolete'
    C2ERR = 'Use HEADER_OVERRIDE_FILE'
    CALL MADABORT("RD_VPEC_FILE",C1ERR,C2ERR)   ! Jun 8 2020
    RETURN
  END SUBROUTINE RD_VPEC_FILE


! ================================
    SUBROUTINE REDSHIFT_HD()

! Created Oct 2013 by R.Kessler.
! Compute redshift 'ZHD' for  Hubble diagram, which is
! ZCMB + vPec correction.
! 
! Nov 30 2016:
!  + re-factor to use externally-computed SNLC_VPEC.
!    i.e., move the read-vpec-file stuff into RD_VPEC_FILE.
! 
! 
! Jan 26 2018: set SNLC_ZPEC[_ERR]
! 
! ----------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    LOGICAL  DOCOR_VPEC
    REAL     ZERR1, ZERR2, ZPEC, ZPECERR

! ------------ BEGIN ---------

! init default with ZHD = CMB redshift
    SNLC_ZHD        = SNLC_ZCMB
    SNLC_ZHD_ERR    = SNLC_ZCMB_ERR

    DOCOR_VPEC = .FALSE.
    IF ( SNLC_VPEC     .NE. 0.0 ) DOCOR_VPEC = .TRUE.
    IF ( SNLC_VPEC_ERR .NE. 0.0 ) DOCOR_VPEC = .TRUE.

    SNLC_ZPEC     = SNLC_VPEC     / SNGL(CLIGHT)
    SNLC_ZPEC_ERR = SNLC_VPEC_ERR / SNGL(CLIGHT)

! ------------------------------------------
! re-compute ZHD and its error if SNLC_VPEC has non-null value
! Jan 3 2018: use VPEC (sigma_z) correction from Eq. 1 (A1) in
!             Davis 2012 [https://arxiv.org/pdf/1012.2912.pdf]

    IF ( DOCOR_VPEC ) THEN

       ZPEC     = ( SNLC_VPEC     / SNGL(CLIGHT) )
       ZPECERR  = ( SNLC_VPEC_ERR / SNGL(CLIGHT) )

       if ( RESTORE_WRONG_VPEC ) then
!    use incorrect VPEC convention ... to be deleted someday
          SNLC_ZHD = (1.0+SNLC_ZCMB) * (1.0+ZPEC) - 1.0  ! Jan 2018;
       else
!    use correct VPEC convention
          SNLC_ZHD = (1.0+SNLC_ZCMB) / (1.0+ZPEC) - 1.0  ! Oct 2020
       endif

       ZERR1 = SNLC_ZCMB_ERR
       ZERR2 = ZPECERR * (1.0 + SNLC_ZCMB) ! Eq A1, Davis 2012
       SNLC_ZHD_ERR = sqrt(ZERR1*ZERR1 + ZERR2*ZERR2)
    ENDIF

    RETURN
  END SUBROUTINE REDSHIFT_HD

! ============================================
    SUBROUTINE SNRECON()
! 
! Created Feb 10, 2006 by R.Kessler
! 
! Miscellaneous reconstruction/computation of
! useful variables for analysis.
! 
! Includes:
! 
!  * ISNLC_NFILT_TRESTMIN[MAX]
!  * ISNLC_NFILT_TREST2
!  * SNLC_TREST(epoch)        = MJD - MJDatPEAK(search) / 1+z
!  * SNLC_TOBS(epoch)         = MJD - MJDatPEAK(search)
!  *
!  * SNLC_FLUXCAL_ERRCALC(ifilt,epoch,isn)  ! calc flux error
!  * SNHOST_ZPHOT[_ERR]
!  * SNLC_SNRMAX_FILT(ifilt)
!  * SNLC_SNRMAX_SORT(rank)
!  * SNLC_SNRMAX_IFILTDEF(ifilt)
!  * SNLC_FLUXCALMAX(ifilt)
!  * SNLC_FLUXCALMAX_ERR(ifilt)  ! Nov 2022
! 
!  * correct SNLC_MAG and FLUXCAL for AB mag-offsets
! 
!  * erase filter-epochs with bad photometry flag
! 
! 
! 
! Feb 12 2018: for USE_SNHOST_ZPHOT=T, convert host zphot to host zcmb
! Feb 16 2018: compute SNHOST_SBMAG(ifilt)
! Mar 18 2017: fill SNLC_DTOBS[_SAMEFILT] = time since last obs
! Sep 16 2018: compute NOBS_DETECT and TLIVE_DETECT
! 
! Nov 12 2018: check APPLY_FLUXCOR_BUG to apply bug of computing
!              SNRMAX variables before EXEC_FUDGE_FLUXCAL()
! 
! Feb 08 2020: apply PSF, ZPERR, PHOTPROB cuts here so that USEFLAG is set.
! Feb 12 2020: check MWEBV_FORCE
! Nov 30 2020: compute SIM_EPCHI2FLUX for SNANA+EPOCHS table
! Feb 9  2022: fix bug setting SUBSURVEY_NAME_LIST
! -------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PRIVCOM

    IMPLICIT NONE

! local var

    INTEGER  & 
         newmjd, ifilt, ifilt_obs, LENLIST, LENNAME, INDEX_OVP  & 
        ,NBAND_SNRMAX, NBAND_SNRMAX2  & 
        ,EPMIN, EPMAX, ep, isort, PHOTFLAG  & 
        ,ORDER, INDEX_SORT(MXFILT_OBS)  ! SORT args

    REAL  & 
         Fluxcal, Fluxcal_err, TRUEFLUX, ERRTEST  & 
        ,Xnsig, MAG, SIM_MAG, magoff, ZP, PHOTPROB, TMP  & 
        ,Tobs, z1, z, Trest, scale, arg, PSF, ZPERR  & 
        ,Zhost, Zhosterr, LAMOBS, LAMZ1,  SNRMAX, SUM_AREAFRAC

    REAL*8  MJD8, EBV8, EBVERR8
    REAL*8  MJD8_LASTALL, MJD8_LASTFILT(MXFILT_ALL)
    REAL*8  MU_OLD, MU_NEW, z8
    LOGICAL  & 
         LSIG, LLAM, LZ, LTEST, LXMJD, LINCMJD  & 
        ,LPSF, LZPERR, LPP, LERRTEST  & 
        ,LSNRMAX(MXFILT_ALL)  & 
        ,LSNRMAX2(MXFILT_ALL)  & 
        ,LFLUX, LERR, LTMP, USE4SNRMAX

! function
    REAL*8   DLMAG8_REF
    EXTERNAL SORTFLOAT, modify_MWEBV_SFD
    INTEGER  ISTAT_REQUIRE_EPOCHS

! ------------------ BEGIN -----------------

! check to increment SUBSURVEY_LIST
! WARNING: If CFA follows CFA3, the code will think that CFA already exists as a survey
    IF ( IDSURVEY .ne. IDSUBSURVEY ) THEN
      LENNAME = INDEX(SUBSURVEY_NAME,' ') - 1
      INDEX_OVP = INDEX(SUBSURVEY_NAME_LIST,SUBSURVEY_NAME(1:LENNAME))
      IF ( INDEX_OVP<=0 ) THEN
        LENLIST = INDEX(SUBSURVEY_NAME_LIST,' ') - 1
        if ( LENLIST == 0 ) then
           SUBSURVEY_NAME_LIST = SUBSURVEY_NAME
        else
          SUBSURVEY_NAME_LIST = SUBSURVEY_NAME_LIST(1:LENLIST)  & 
                // ',' // SUBSURVEY_NAME
        endif
      ENDIF
    ENDIF

! - - - - - - - - -
! prepare redshifts; e.g. use photo-z, apply systematic shifts, etc ...
    CALL SET_REDSHIFTS()

! - - - - -
    Z    = SNLC_REDSHIFT

    if ( Z > 0.0 ) then
       Z1 = 1.0  + z
    else
       Z1 = 1.0  ! in case Z = -9 is undefined; Trest -> Tobs
    endif

    ISNLC_NEPOCH_USE = 0


! ------ compute MJDMIN/MJDMAX here before calling SET_PKMJD ------

    SNLC8_MJDMIN  = 1.0E8
    SNLC8_MJDMAX  = 0.0
    DO NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

      DO EP = EPMIN, EPMAX
        SNLC8_FLUXCAL(EP)        = DBLE(SNLC_FLUXCAL(ep))
        SNLC8_FLUXCAL_ERRTOT(EP) = DBLE(SNLC_FLUXCAL_ERRTOT(ep))
        MJD8  = SNLC8_MJD(EPMIN)
        IF( MJD8 .GT. SNLC8_MJDMAX ) SNLC8_MJDMAX = MJD8
        IF( MJD8 .LT. SNLC8_MJDMIN ) SNLC8_MJDMIN = MJD8
     ENDDO
    ENDDO

! Sep 21 2013: check option to modify MWEBV_SFD;
! Note that MWEBV and MWEBV_ERR are modified if OPT_MWEBV > 0
    EBV8 = DBLE(SNLC_MWEBV) ; EBVERR8 = DBLE(SNLC_MWEBV_ERR)

    CALL modify_MWEBV_SFD(OPT_MWEBV, SNLC8_RA, SNLC8_DEC,  & 
                            EBV8, EBVERR8)  ! I -> O

    SNLC_MWEBV     = SNGL(EBV8)     * MWEBV_SCALE + MWEBV_SHIFT
    SNLC_MWEBV_ERR = SNGL(EBVERR8)  * MWEBV_SCALE
    IF ( MWEBV_FORCE > -0.001 ) THEN
       SNLC_MWEBV = MWEBV_FORCE
       SNLC_MWEBV_ERR = 0.0
    ENDIF

! compute Galactic extinction in each band
! (Oct 11 2013: code moved from MAIN to here)
    CALL MWEBV_FLUXCOR()

! determine redshift for Hubble diagram (Oct 2013)
    CALL REDSHIFT_HD()

#if defined(MINUIT)
! get approx peakMJD
    if ( OPT_SETPKMJD > 0 ) THEN
       CALL SET_PEAKMJD()
    endif
#endif

    CALL SELECT_EPOCH_DRIVER(2) ! apply Trest cuts

! -----------------------------------------------
! pick host photoZ

    Zhost    = -9.0
    Zhosterr = -9.0

    if ( SNHOST_ZPHOT(1) .GT. 0.0 ) then
       zhost    = SNHOST_ZPHOT(1)
       zhosterr = SNHOST_ZPHOT_ERR(1)
    else
      SNHOST_ZPHOT(1)     = Zhost
      SNHOST_ZPHOT_ERR(1) = Zhosterr
    endif

     CALL SET_SNHOST_CONFUSION()

! compute lumi-distance using "standard" cosmology
! Useful as initial value for fits.

    IF ( Z .GT. 1.0E-5 ) THEN
       z8 = dble(z)
       SNLC_DLMAG  = SNGL(DLMAG8_REF(z8) )
    ELSE
       SNLC_DLMAG  = -9.0  ! undefined
    ENDIF


    do ifilt = 1, NFILTDEF_SURVEY
       LSNRMAX(ifilt)       = .FALSE.
       LSNRMAX2(ifilt)      = .FALSE.
       MJD8_LASTFILT(ifilt) = -9.0
       CALL SET_SBMAG(ifilt)  ! Feb 17 2018: compute SBMAG from SBFLUX
    enddo

    NBAND_SNRMAX  = 0
    NBAND_SNRMAX2 = 0
    MJD8_LASTALL  = -9.0
    NEP_SIM_MAGOBS =  0
    SUM_AREAFRAC   = 0.0

! fill ISNLC_NFIELD with number of fields;
    CALL COUNTFIELDS()

! init computation of PROB_TRUEFLUX
    CALL PROB_TRUEFLUX_CALC(0)

! = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! = = = = = = = = START LOOP OVER EPOCHS = = = = = = = =
! = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

      MJD8  = SNLC8_MJD(EPMIN)
      Tobs  = SNGL(MJD8) - SNLC_SEARCH_PEAKMJD
      Trest = Tobs / z1

      SNLC_TOBS(EPMIN)  = Tobs
      SNLC_TREST(EPMIN) = Trest

      DO 201 EP = EPMIN, EPMAX

         ISNLC_SNRECON_USE(ep) = 0

         SNLC_TREST(EP) = Trest
         SNLC_TOBS(EP)  = Tobs
         IFILT_OBS = ISNLC_IFILT_OBS(ep)
         IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

         CALL CHECK_DUPLICATE_MJDBAND(ep) ! Jun 13 2017

! Feb 6 2020: fudge errors here before SNRMAX is computed

         CALL EXEC_FUDGE_FLUXCAL(ep)  ! flux and fluxerr fudges

         CALL SET_MAG(ep)  ! Jan 23 2018

         USE4SNRMAX = LFILTDEF_SNRMAX(ifilt_obs) ! Apr 2014

! test filter mean only if KCOR file is read.
         LLAM = .TRUE.

! Mar 22 2023 restore LLAM test that was mistakenly removed when
!             fortran-KCOR code was replaced with C code.
! 
        LZ = ( Z > 0 .and. Z < 998 )
        IF ( LZ .and. EXIST_CALIB_FILE ) THEN
          LAMOBS = FILTOBS_LAMAVG(ifilt_obs)
          LAMZ1  = LAMOBS / Z1
          LLAM   =  & 
                 LAMZ1  .GE. CUTWIN_LAMREST(1) .and.  & 
                 LAMZ1  .LE. CUTWIN_LAMREST(2) .and.  & 
                 LAMOBS .GE. CUTWIN_LAMOBS(1)  .and.  & 
                 LAMOBS .LE. CUTWIN_LAMOBS(2)
        ENDIF
! 
! strip off local variables.
        Fluxcal      = SNLC_FLUXCAL(ep)
        Fluxcal_err  = SNLC_FLUXCAL_ERRTOT(ep)
        mag          = SNLC_MAG(ep)

! Jan 3 2018: skip saturated mags
        if ( FLUXCAL_ERR .GE. 0.9E8 .AND. LSIM_SNANA ) GOTO 201

! convert sim-mag to SIM_FLUXCA
        IF ( ISJOB_SIM ) THEN
           ZP        = ZP_FLUXCAL
           SIM_MAG   = SIM_EPMAGOBS(ep)
           if ( SIM_MAG < 90.0 ) then
             NEP_SIM_MAGOBS = NEP_SIM_MAGOBS + 1 ! Nov 25 2019
             ARG            = 0.4*(ZP - SIM_MAG)
             TRUEFLUX       = 10.0**(ARG)
           else
             TRUEFLUX = 0.0
           endif
           SIM_EPFLUXCAL(ep) = TRUEFLUX
           IF ( USESIM_TRUEFLUX ) SNLC_FLUXCAL(ep)=TRUEFLUX

           TMP  = (TRUEFLUX - FLUXCAL)/ FLUXCAL_ERR
           SIM_EPCHI2FLUX(ep) = TMP * TMP  ! Nov 30 2020
           SIM_EPPULL(ep)     = TMP        ! Mar 16 2021
           CALL PROB_TRUEFLUX_CALC(ep)
        ELSE
           SIM_EPFLUXCAL(ep)  = -9.0
           SIM_EPCHI2FLUX(ep) = -9.0
           SIM_EPPULL(ep)     = -9.0
        ENDIF

#if defined(SNANA) || defined(SNFIT)
        CALL FLUXERRCALC( ep, Fluxcal, Fluxcal_err )
#endif

! compute significance.
        Xnsig  = -999.0
        LERR   = (Fluxcal_err > 0.0 )
        LFLUX  = (abs(Fluxcal+9.0) > 1.0E-5 )  !  FLUXCAL != -9
        if ( LERR .and. LFLUX ) then
           Xnsig  = Fluxcal / Fluxcal_err
        endif

! ----------------
! Correct MAG and FLUXCAL  for  MAGOFF_AB.
! This correction is done AFTER Xsig = S/N is calculated

         magoff = MAGOBS_SHIFT_ZP_FILT(ifilt_OBS)
         SNLC_MAG(ep) = SNLC_MAG(ep) + magoff
         FLUXSCALE_ZP_FILT(IFILT_OBS) = 10.0**(-0.4 * magoff)

         if ( magoff .ne. 0.0 ) then
           scale = FLUXSCALE_ZP_FILT(IFILT_OBS)
           SNLC_FLUXCAL(ep) =  & 
             SNLC_FLUXCAL(ep) * Scale
           SNLC_FLUXCAL_ERRTOT(ep) =  & 
             SNLC_FLUXCAL_ERRTOT(ep) * Scale
         endif


! Mar 2018: store time since last obs: all bands, and same band
        IF ( MJD8_LASTALL > 0.0 ) THEN
           SNLC_DTOBS(ep) = SNGL(MJD8 - MJD8_LASTALL)
        ENDIF
        IF ( MJD8_LASTFILT(ifilt) > 0.0 ) THEN
           SNLC_DTOBS_SAMEFILT(ep) = SNGL(MJD8 - MJD8_LASTFILT(ifilt))
        ENDIF
        MJD8_LASTALL         = MJD8
        MJD8_LASTFILT(ifilt) = MJD8

! ----------------
! signal is accepted if it passes S/noise cut,
! or it is within the Trest window.

        SNLC_SNR(ep) = Xnsig

! xxxxxxx mark delete Mar 2025 xxxxxxxx
! check how many obs are in the INCLUDE_MJD window (Aug 2015)
!          LINCMJD = MJD8 .GE. CUTWIN_MJD_INCLUDE(1) .and.
!     &              MJD8 .LE. CUTWIN_MJD_INCLUDE(2)
!          IF ( LINCMJD ) ISNLC_NMJD_INCLUDE = ISNLC_NMJD_INCLUDE + 1
! xxxxxxxxxxxxxxxx

! check MJD-exclude window to reject epochs (but keep LC)
        LXMJD = MJD8 .GE. CUTWIN_MJD_EXCLUDE(1) .and.  & 
                  MJD8 .LE. CUTWIN_MJD_EXCLUDE(2)

! Feb 2020: check PSF and ZPERR cuts here
        PSF  = SNLC_PSF_FWHM_ARCSEC(EP)
        LPSF = PSF .GE. CUTWIN_PSF(1) .and. PSF.LE.CUTWIN_PSF(2)

        ZPERR  = SNLC_ZEROPT_ERR(EP)
        LZPERR = ZPERR .GE. CUTWIN_ZPERR(1) .and.  & 
                   ZPERR .LE. CUTWIN_ZPERR(2)

        PHOTPROB = SNLC_PHOTPROB(ep)
        LPP  = (PHOTPROB > CUTWIN_PHOTPROB(1) .OR. PHOTPROB<-0.001)

        ERRTEST = SNLC_FLUXCAL_ERRTEST(ep)
        LERRTEST = ERRTEST .GE. CUTWIN_ERRTEST(1) .and.  &  ! Mar 2021
                     ERRTEST .LE. CUTWIN_ERRTEST(2)

        LSIG = Xnsig .GE. cutwin_snrmin_filt(1,ifilt)  & 
           .and. Tobs  .GE. cutwin_Tobs(1)  & 
           .and. Tobs  .LE. cutwin_Tobs(2)  & 
           .and. Trest .GE. cutwin_Trest(1)  & 
           .and. Trest .LE. cutwin_Trest(2)  & 
           .and. MJD8  .GE. cutwin_MJD(1)  & 
           .and. MJD8  .LE. cutwin_MJD(2)  & 
           .and. .not. LXMJD      &  ! Feb 2012
           .and. LLAM             &  ! rest-frame lambda OK (7/08/08)
           .and. LPSF .and. LZPERR .and. LPP     &  ! Feb 2020
           .and. LERRTEST                        &  ! Mar 2021
           .and. fluxcal_err >  0.0    ! flux err is > 0

        IF ( .NOT. LSIG ) GOTO 201 ! Nov 30, 2011

! set use-flag for this epoch
        ISNLC_SNRECON_USE(ep) = 1

        if ( USE4SNRMAX ) then
           SNLC_SNRMAX_FILT(0)     =  & 
                  MAX ( SNLC_SNRMAX_FILT(0), Xnsig )
           SNLC_SNRMAX_FILT(ifilt) =  & 
                  MAX ( SNLC_SNRMAX_FILT(ifilt), Xnsig )
        endif

        if ( Fluxcal > SNLC_FLUXCALMAX(ifilt) ) then
           SNLC_FLUXCALMAX(ifilt)     = Fluxcal
           SNLC_FLUXCALMAX_ERR(ifilt) = Fluxcal_err
        endif

        EXIST_FILT(ifilt) = .TRUE.
        ISNLC_NEPOCH_USE  = ISNLC_NEPOCH_USE + 1

! track number of filters that pass SNRMAX cut
        LTEST = XNSIG .GE. CUTWIN_SNRMAX(1)
        LTMP  = USE4SNRMAX .and. (.not. LSNRMAX(ifilt))

        IF ( LTEST .and. LTMP ) THEN
           NBAND_SNRMAX = NBAND_SNRMAX + 1
           LSNRMAX(ifilt) = .TRUE.
        ENDIF

! repeat for 2nd SNRMAX2 cut
        LTEST = XNSIG .GE. CUTWIN_SNRMAX2(1)
        LTMP  = USE4SNRMAX .and. (.not. LSNRMAX2(ifilt))
        IF ( LTEST .and. LTMP ) THEN
           NBAND_SNRMAX2 = NBAND_SNRMAX2 + 1
           LSNRMAX2(ifilt) = .TRUE.
        ENDIF

! check PHOTPROB (Mar 2018)
        PHOTPROB = SNLC_PHOTPROB(ep)
        if ( PHOTPROB > -0.0001 ) then
           SNLC_PHOTPROB_MIN = MIN(SNLC_PHOTPROB_MIN,PHOTPROB)
           ISNLC_NEPOCH_PHOTPROB = ISNLC_NEPOCH_PHOTPROB + 1
        endif

! if PHOTFLAG detection bit is set by user, compute MJD of first & last detection
        PHOTFLAG = ISNLC_PHOTFLAG(ep)
        if ( IAND(PHOTFLAG,PHOTFLAG_DETECT) > 0 ) then
           ISNLC_NOBS_DETECT = ISNLC_NOBS_DETECT + 1
           if ( ISNLC_NOBS_DETECT == 1 ) then
              SNLC8_MJD_DETECT_FIRST = MJD8
           endif
           SNLC8_MJD_DETECT_LAST = MJD8
        endif

! check for trigger bit set by survey (5.2019)
        if ( IAND(PHOTFLAG,PHOTFLAG_TRIGGER) > 0 ) then
           SNLC8_MJD_TRIGGER = MJD8
        endif

#if defined(SNANA) || defined(SNFIT)
        CALL CCD_AREAFRAC(ep)
        IF ( SNLC_AREAFRAC(ep) > -0.001 ) then
           SUM_AREAFRAC = SUM_AREAFRAC + SNLC_AREAFRAC(ep)
        ENDIF
#endif

! Mar 2025 : count NOBS_PREDETECT
        TMP = MJD8 - SNLC8_MJD_DETECT_FIRST  ! TMP < 0 for pre-detection
        if ( TMP .GE.  CUTWIN_TOBS_PREDETECT(1) .and.  & 
               TMP .LE.  CUTWIN_TOBS_PREDETECT(2)  ) then
           ISNLC_NOBS_PREDETECT =  ISNLC_NOBS_PREDETECT + 1
        endif

201     CONTINUE     ! come here of LSIG fails

200   CONTINUE  ! NEWMJD

! - - - - - -  -

    SNLC_TLIVE_DETECT =  & 
         SNGL(SNLC8_MJD_DETECT_LAST - SNLC8_MJD_DETECT_FIRST)

    ISNLC_NFILT_SNRMAX    = NBAND_SNRMAX
    ISNLC_NFILT_SNRMAX2   = NBAND_SNRMAX2

    ORDER = -1
    CALL SORTFLOAT( NFILTDEF_SURVEY, SNLC_SNRMAX_FILT(1),  & 
                       ORDER, INDEX_SORT(1) )

    DO isort  = 1, NFILTDEF_SURVEY
       ifilt  = INDEX_SORT(isort)
       SNRMAX = SNLC_SNRMAX_FILT(ifilt)
       SNLC_SNRMAX_SORT(isort) = SNRMAX
    ENDDO

! -----------------------------------------
    CALL GET_USERTAG(SNLC_CCID )

! ----------------------------------------
! Oct 2014: check for multi-season activity

    CALL MULTISEASON(IFLAG_ANA)

! check epochs above/below threshold
    ISNLC_CUTFLAG_REQEP = ISTAT_REQUIRE_EPOCHS()

! evalute cuts on PRIVATE & SIMVAR variables
    CALL EVAL_PRIVATEVAR_CUTS()
    CALL EVAL_SIMVAR_CUTS()

! check option to estimate LC width
    CALL GET_SIM_LCWIDTH(1)

! check option to compute PROB_TRUEFLUX_CALC per band
    CALL PROB_TRUEFLUX_CALC(-1)

! May 2017: load ADDCOL arrays
    CALL TABLE_ADDCOL_LOAD()

    IF ( ISNLC_NEPOCH_USE > 0 ) then
      SNLC_AREAFRAC_AVG = SUM_AREAFRAC/float(ISNLC_NEPOCH_USE)
    ENDIF

    RETURN
  END SUBROUTINE SNRECON

! =============================
    SUBROUTINE SET_REDSHIFTS()
! 
! Created Jun 2021
! Called near beginning of SNRECON, this is a driver routine
! to finalize redshifts:
!  + option to use host photo-z
!  + fix for missing zHELIO (on very old surveys)
!  + systematic shifts
!  + etc ...
! 
! Aug 23 2025: check USESIM_REDSHIFT option (moved from snlc_fit to here)
! 
    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE
    INTEGER IERR

! -------------- BEGIN ----------

    if ( USESIM_REDSHIFT ) then  ! Aug 23 2025
       SNLC_REDSHIFT = SIM_REDSHIFT_HELIO
       SNLC_zHELIO   = SIM_REDSHIFT_HELIO
       SNLC_zCMB     = SIM_REDSHIFT_CMB
       return
    endif

! Nov 30 2025: if there are quantiles, compute MEAN and STDDEV 
    CALL SET_SNHOST_QZPHOT(METHOD_SPLINE_QUANTILES_DEFAULT,IERR)

! Check option to use host photo-z as redshift.
    CALL SET_SNHOST_ZPHOT()

! if ZHELIO < 0, convert zCMB -> zHELIO
    CALL SET_ZHELIO()

! check for systematic shifts
    CALL SET_zSHIFT()

    RETURN
  END SUBROUTINE SET_REDSHIFTS

! =============================
    SUBROUTINE SET_ZHELIO()

! Creatd Feb 25 2018
! If SNLC_ZHELIO < 0, transform from zCMB.
! Otherwise, check zHEL <-> zCMB transformation.
! If zCMB < 0 for TEXT format, compute zCMB from zHEL.
! 
!  Jun 8 2020
!   + fix awful bug with sign of OPT sent to translator.
!      was translating zHELIO -> CMB instead of zCMB -> zHELIO
!   + if SNLC_HELIO[CMB] are both given in data file, make
!     crosscheck here using zTOL_HELIO2CMB test on dz/z
! 
! Mar 15 2021: if zCMB < 0 for text format, use computed value.
! Jun 07 2022: do not abort on bad tolerance if REFORMAT_SNANA=T
! Jan 28 2024: do NOT check zTOL if zERR > 1.0 due to sim artifact
!              from GENSIGMA_REDSHIFT option;
!                https://github.com/RickKessler/SNANA/issues/1265
! 
! ---------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT
    CHARACTER EQ*4, FNAM*12
    REAL*8  ZHEL8, ZCMB8, DZ8, zRATIO8
    LOGICAL BADz, SKIP_BADz

    REAL*8 zhelio_zcmb_translator ! function

! ------------- BEGIN ------------

    FNAM = 'SET_ZHELIO'

! now we have undefined zCMB

    EQ     = 'eq' // char(0)

! don't apply test on VERY low redshifts to avoid silly round-off
! or truncation errors. This is a hack to protect DES DIFFIMG that
! has redshifts written as %.4f. Next dataFile update will have %.6f
! precisions and then can apply test all the way down to z ~ 1E-4

    IF ( SNLC_ZHELIO > 0.003 ) THEN  ! 0.003 is temp hack; should be 1E-4
       OPT    = +1  ! --> convert ZHEL to ZCMB
       ZHEL8 = DBLE(SNLC_ZHELIO)
       ZCMB8 =  & 
            zhelio_zcmb_translator(ZHEL8,SNLC8_RA,SNLC8_DEC,EQ,OPT, 4)

! For TEXT format only, set zCMB if not provided (Mar 15 2021);
! allows quick/sloppy makeDataFiles with real data.
       IF ( SNLC_zCMB < 0.0 .and. FORMAT_TEXT ) THEN
         SNLC_zCMB     = SNGL(ZCMB8)
         SNLC_zCMB_ERR = SNLC_ZHELIO_ERR
       ENDIF

       zRATIO8 = zCMB8/DBLE(SNLC_zCMB)
       DZ8     = abs(1.0 - zRATIO8)

       BADz      = ( SNLC_zCMB > 0.0 .and. DZ8 > zTOL_HELIO2CMB )
       SKIP_BADz = REFORMAT_SNANA .or. (SNLC_ZHELIO_ERR .GE. 0.9999)

       if ( BADz .and. .NOT.SKIP_BADz ) THEN
          CALL PRINT_PREABORT_BANNER(FNAM(1:10)//char(0),12)
          print*,' zHEL(dataFile) = ', SNLC_ZHELIO
          print*,' zCMB(dataFile) = ', SNLC_ZCMB
          print*,' zCMB(compute)  = ', ZCMB8
          print*,' zRATIO = zCMB(compute)/zCMB(dataFile) = ',  & 
                       sngl(zRATIO8)
          print*,' zTOL_HELIO2CMB = ', zTOL_HELIO2CMB
          C1ERR='Computed zCMB (from zHEL) outside tolerance w.r.t.'
          C2ERR='REDSHIFT_CMB in data file for CID=' // SNLC_CCID
          CALL MADABORT(FNAM,C1ERR,C2ERR)
       endif

    ELSE if ( SNLC_ZHELIO < 0.0 ) THEN
       OPT    = -1 ! --> convert CMB -> zHEL
!          OPT    = +1 ! --> BUG
       ZCMB8  = DBLE ( SNLC_ZCMB )
       ZHEL8  =  & 
          zhelio_zcmb_translator(ZCMB8,SNLC8_RA,SNLC8_DEC,EQ,OPT, 4)

       SNLC_ZHELIO       = SNGL(ZHEL8)
       SNLC_ZHELIO_ERR   = SNLC_ZCMB_ERR
       SNLC_REDSHIFT     = SNLC_ZHELIO
       SNLC_REDSHIFT_ERR = SNLC_ZHELIO_ERR
    ENDIF

    RETURN
  END SUBROUTINE SET_ZHELIO


! =============================================
    SUBROUTINE SET_SNHOST_QZPHOT(METHOD_SPLINE_QUANTILES,IERR)
!
! Created Nov 30 2025:
! If there are host photot-z quanitiles, compute MEAN and STDDEV
!  [code moved from snlc_fit.F90 to here so that MEAN and STDDEV
!    appear in SNANA table without having to do LC fit]

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER :: METHOD_SPLINE_QUANTILES*(*)  ! (I) method to interpolate PDF
    INTEGER   :: IERR                         ! (O) return error code (0 = no error)

! local var
    INTEGER   :: q, LM, IERR_ZPDF, IPRINT
    REAL*8    :: ZPHOT_Q(MXZPHOT_Q), ZPHOT_PROB(MXZPHOT_Q), MEAN, STD
    CHARACTER :: CCID_forC(MXCHAR_CCID)


! ------------- BEGIN ---------------

    if ( SNHOST_NZPHOT_Q <= 0 ) RETURN

    do q = 1, SNHOST_NZPHOT_Q
       ZPHOT_PROB(q) = DBLE(SNHOST_ZPHOT_PERCENTILE(q)) / 100.  ! 0 <= PROB <= 1
       ZPHOT_Q(q)    = DBLE(SNHOST_ZPHOT_Q(1,q))
    enddo
    
    LM = INDEX(METHOD_SPLINE_QUANTILES,' ') - 1
    CCID_forC = SNLC_CCID(1:ISNLC_LENCCID) // char(0)
    IPRINT    = 0   ! set to 1 for dump

    CALL init_zPDF_spline(SNHOST_NZPHOT_Q, ZPHOT_PROB, ZPHOT_Q,  &
         CCID_forC, METHOD_SPLINE_QUANTILES(1:LM)//char(0),  &
         IPRINT, MEAN, STD, IERR, ISNLC_LENCCID, 20)

    if (IERR .EQ. 0 ) then
       SNHOST_QZPHOT_MEAN(1) = MEAN ! store mean & std in 4 byte global
       SNHOST_QZPHOT_STD(1)  = STD
    endif

    return
    END SUBROUTINE SET_SNHOST_QZPHOT

! =============================================
    SUBROUTINE SET_SNHOST_ZPHOT()

! if USE_SNHOST_ZPHOT=T then set redshift variables to host zphot.
! This ignores accurate specz to enable using photo-z.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT
    CHARACTER EQ*4
    REAL*8  ZCMB8, ZHEL8

    REAL*8 zhelio_zcmb_translator ! function
! -------------- BEGIN ----------------

    IF ( .not. (USE_SNHOST_ZPHOT .or. USE_HOSTGAL_PHOTOZ) ) RETURN  ! same variable meaning ???

    if ( .not. EXIST_SNHOST_ZPHOT ) then
        C1ERR = 'Cannot USE_HOSTGAL_PHOTOZ for CID='//SNLC_CCID
        C2ERR = 'Check data files, or set USE_HOSTGAL_PHOTOZ=F'
        CALL MADABORT("SNRECON", c1err, c2err )
    endif

    OPT    = 1     ! --> convert ZHEL to ZCMB
    EQ     = 'eq' // char(0)
    ZHEL8  = DBLE( SNHOST_ZPHOT(1) )
    ZCMB8  = zhelio_zcmb_translator(ZHEL8, SNLC8_RA, SNLC8_DEC, EQ, OPT, 4)

    SNLC_ZCMB         = SNGL(ZCMB8)
    SNLC_ZHELIO       = SNGL(ZHEL8)
    SNLC_REDSHIFT     = SNGL(ZHEL8)

    SNLC_ZHELIO_ERR     = SNHOST_ZPHOT_ERR(1)
    SNLC_ZCMB_ERR       = SNHOST_ZPHOT_ERR(1)
    SNLC_REDSHIFT_ERR   = SNHOST_ZPHOT_ERR(1)

    RETURN
  END SUBROUTINE SET_SNHOST_ZPHOT

! =============================================
    SUBROUTINE SET_zSHIFT

! Created Jun 5 2021
! Apply systematic z-shifts passed from &SNLCINP
! 
! Dec 29 2024: allow ztol of 1E-5 to allow for slop in round-off in TEXT format
! Apr 22 2025: fix bug and shift ZPHOT_Q (quantiles)


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    LOGICAL LSHIFT_FINAL, LSHIFT_HOST_ZPHOT, LSHIFT_HOST_ZSPEC
    INTEGER OPT, q, igal
    REAL    zshift, ztol, zdif

    REAL*8 zhelio_zcmb_translator ! function
! -------------- BEGIN ----------------

! bail if none of the zSHIFT options are used
    if ( .not. DOzSHIFT ) RETURN

    LSHIFT_FINAL       = ( REDSHIFT_FINAL_SHIFT .NE. 0.0 )
    LSHIFT_HOST_ZSPEC  = ( HOSTGAL_ZSPEC_SHIFT  .NE. 0.0 )
    LSHIFT_HOST_ZPHOT  = ( HOSTGAL_ZPHOT_SHIFT  .NE. 0.0 )

    zshift = -9.0
    IF ( LSHIFT_FINAL ) THEN
       zshift  = REDSHIFT_FINAL_SHIFT
    ENDIF

! Dec 29 2024: allow tolerance for round-off in text format
    ZTOL = 1.0E-5

    IF ( LSHIFT_HOST_ZSPEC ) THEN
       zdif = SNHOST_ZSPEC(1) - SNLC_REDSHIFT
       if ( abs(zdif) < ztol ) then
          zshift = HOSTGAL_ZSPEC_SHIFT + REDSHIFT_FINAL_SHIFT
          SNHOST_ZSPEC = SNHOST_ZSPEC + zshift
       endif
    ENDIF

    IF ( LSHIFT_HOST_ZPHOT ) THEN
       zdif = SNHOST_ZPHOT(1) - SNLC_REDSHIFT
       if ( abs(zdif) < ztol ) then
          zshift = HOSTGAL_ZPHOT_SHIFT + REDSHIFT_FINAL_SHIFT
          SNHOST_ZPHOT = SNHOST_zPHOT + zshift
       endif

       if ( SNHOST_NZPHOT_Q > 0 ) then
          do q    = 1, SNHOST_NZPHOT_Q
          do igal = 1, MXSNHOST
	      SNHOST_ZPHOT_Q(igal,q) = SNHOST_ZPHOT_Q(igal,q) + zshift ! Apr 22 2025
          enddo
          enddo
       endif
    ENDIF

    IF ( zshift .NE. -9.0 ) THEN
       if ( SNLC_REDSHIFT > 0.) SNLC_REDSHIFT = SNLC_REDSHIFT+ zshift
       if ( SNLC_ZHELIO   > 0.) SNLC_ZHELIO   = SNLC_ZHELIO  + zshift
       if ( SNLC_ZCMB     > 0.) SNLC_ZCMB     = SNLC_ZCMB    + zshift
    ENDIF

    RETURN
  END SUBROUTINE SET_zSHIFT

! ============================================
    LOGICAL FUNCTION IS_ZSPEC()

! Created March 20 2024 by R.Kessler and A.Mitra
! Return True if best redshift is a spec-z.
! Optimal logic is to use MASK_REDSHIFT_SOURCE that
! explcitly sets bit0 for host zspec and bit1 for SN zspec.
! However, if this MASK is not provided, one can assume
! specz if redshift error is less than &SNLCINP variable
! QUANTILE_ZERRMIN.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! define local zspec mask;
!   grep MASK_REDSHIFT_SOURCE sndata.h
    INTEGER, PARAMETER :: MASK_ZSPEC_SOURCE = 3   ! 1=zspec_host, 2=zspec(SN)
    LOGICAL  IS_SMALL_ZERR, IS_MASK_ZSPEC

    IS_SMALL_ZERR = (SNLC_REDSHIFT_ERR < QUANTILE_ZERRMIN) .and. (SNLC_REDSHIFT_ERR >= 0.0)
    IS_MASK_ZSPEC = (IAND(ISNLC_zSOURCE,MASK_ZSPEC_SOURCE) > 0)

    IS_ZSPEC = IS_SMALL_ZERR .or. IS_MASK_ZSPEC

    RETURN
  END FUNCTION IS_ZSPEC

! ===============================
    SUBROUTINE SET_SNHOST_CONFUSION()
! 
! Created Sep 2 2022
! Compute HOST_CONFUSION based on Eq 3 in Gupta et al, 2016
! 


    USE SNDATCOM

    IMPLICIT NONE

    REAL*8   HC, DDLR_LIST(MXSNHOST)
    INTEGER  igal, NMATCH2
    LOGICAL LDMP / .FALSE. /
    character cCID*40

    REAL*8   HOST_CONFUSION
    EXTERNAL HOST_CONFUSION   ! see sntools.c

! ----------- BEGIN -------------

    SNHOST_CONFUSION = -9.0

    NMATCH2 = SNHOST_NMATCH2

! only the first MXSNHOST matches are stored
    if ( NMATCH2 > MXSNHOST ) NMATCH2 = MXSNHOST

    if ( NMATCH2 .LE. 1 ) RETURN

! convert single precions to double precision
    DO igal = 1, NMATCH2
      DDLR_LIST(igal) = DBLE(SNHOST_DDLR(igal))
    ENDDO
    cCID = SNLC_CCID(1:ISNLC_LENCCID) // char(0)
    HC = HOST_CONFUSION(cCID,NMATCH2,DDLR_LIST,40)

    IF ( LDMP ) THEN
       print*,' xxx ----------------------------------------- '
       print*, ' xxx DUMP for SET_SNHOST_CONFUSION '
       write(6,60) SNLC_CCID(1:12), NMATCH2,  & 
               DDLR_LIST(1), DDLR_LIST(2)
60       format('  xxx CID=',A, 3x, 'NMATCH=',I2, 3x,'DDLR=', 2F9.3 )
       print*,' xxx HOST_CONFUSION = ', HC
       print*,' xxx '
 	 CALL FLUSH(6)
    ENDIF

    SNHOST_CONFUSION = SNGL(HC)

    RETURN
  END SUBROUTINE SET_SNHOST_CONFUSION

! ===============================
    SUBROUTINE SET_MAG(ep)

! Created Jan 23 2018
! Compute observer mag SNLC_MAG(ep) to reduce
! data volume in files. Note that SNLC_MAGERR
! is not set because it is not used.
! 
! Aug 2022: check LDMP_SATURATE
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER ep  ! (I) epoch index

    REAL FLUXCAL, FLUXCAL_ERR, MAG, ZP
    INTEGER IFILT_OBS
    CHARACTER CFILT*2
! --------- BEGIN ------------

    Fluxcal      = SNLC_FLUXCAL(ep)
    Fluxcal_err  = SNLC_FLUXCAL_ERRTOT(ep)

    IF ( FLUXCAL_ERR > 0.99E8 ) THEN
      MAG = MAG_SATURATE
      if ( LDMP_SATURATE ) then
           IFILT_OBS = ISNLC_IFILT_OBS(ep)
           cfilt     = FILTDEF_STRING(ifilt_obs:ifilt_obs)
           write(6,20) SNLC_CCID(1:8), SNLC_ZHELIO,  & 
                         SNLC8_MJD(ep), CFILT(1:1), SIM_EPMAGOBS(ep)
20           format(T2,'SATURATION: ',  & 
                'CID=',A8, 2x, 'zHEL=',F6.4, 2x,  & 
                'MJD=',F9.3,'-',A, 2x, 'TrueMag=',F8.3 )
      endif
    ELSE IF ( FLUXCAL < 1.0E-9 ) THEN
      MAG = 128.0  ! negative flux
    ELSE
      ZP  = SNGL(ZP_FLUXCAL)
      MAG = ZP - 2.5*LOG10(FLUXCAL)
    ENDIF

    SNLC_MAG(ep) = MAG

    RETURN
  END SUBROUTINE SET_MAG


! ===============================
    SUBROUTINE SET_SBMAG(ifilt)

! Created Feb 17 2018
! Compute surface-brightness mag SNLC_SBMAG(ep) from SB flux.
! 
! Mar 23 2021: SBFLUXMIN -> 1.0 (was 0.1)
! Aug 11 2021: SBFLUX = -99 -> SBMAG=99
! 

    USE SNDATCOM

    IMPLICIT NONE

    INTEGER ifilt  ! (I) sparse filter index
    REAL SBFLUX, SBMAG, SBFLUXMIN, ZP
! --------- BEGIN ------------
    SBFLUX = SNHOST_SBFLUXCAL(ifilt)
    SBFLUXMIN = 1.0

    ZP = SNGL( ZP_FLUXCAL )
    IF ( SBFLUX > SBFLUXMIN ) THEN
       SBMAG = ZP - 2.5*LOG10(SBFLUX)
    ELSE
       SBMAG = ZP - 2.5*LOG10(SBFLUXMIN)
    ENDIF

    IF ( SBFLUX .LE. -999.0 ) SBMAG = 99.0 ! Aug 2021 - undefined SB

    SNHOST_SBMAG(ifilt) = SBMAG
    RETURN
  END SUBROUTINE SET_SBMAG

! ========================
    SUBROUTINE PARSE_SIMVAR_CUTS()

! Dec 2018
! Parse &SNLCINP input SIMVAR_CUTWIN_STRING


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER NTMP, i, NWD, LWD
    CHARACTER ctmp1*60, ctmp2*60, ctmp3*60, FNAM*20

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! --------- BEGIN ---------
    IF ( .NOT. ISJOB_SIM )  RETURN
    IF ( SIMVAR_CUTWIN_STRING .EQ. '' ) RETURN

    FNAM = 'PARSE_SIMVAR_CUTS'
    NTMP = 0
    LWD = MXCHAR_CUTNAME-2
    NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,  & 
           SIMVAR_CUTWIN_STRING(1:LWD)//char(0), FNAM//char(0), LWD, 20)

    DO i = 1, NWD, 3
       NTMP = NTMP + 1
       CALL get_PARSE_WORD_fortran(i+0, ctmp1, LWD)
       CALL get_PARSE_WORD_fortran(i+1, ctmp2, LWD)
       CALL get_PARSE_WORD_fortran(i+2, ctmp3, LWD)
       SIMVAR_CUTWIN_LIST(NTMP) = CTMP1(1:40)
       READ(CTMP2,*)  SIMVAR_CUTWIN(1,NTMP)
       READ(CTMP3,*)  SIMVAR_CUTWIN(2,NTMP)
    ENDDO
    NSIMVAR_CUTWIN = NTMP

    if ( NTMP > MXCUT_PRIVATE ) then
       write(c1err,61) NTMP, MXCUT_PRIVATE
 61      format('NSIMVAR_CUTWIN=',I2,' exceeds bound of ', I2)
       c2err = 'Check &SNLCINP input SIMCAR_CUTWIN_STRING'
       CALL MADABORT(FNAM, c1err, c2err )
    endif

    DO i = 1, NTMP
       LWD = INDEX(SIMVAR_CUTWIN_LIST(i),' ') - 1
       write(6,77) SIMVAR_CUTWIN_LIST(i)(1:LWD),  & 
              SIMVAR_CUTWIN(1,i), SIMVAR_CUTWIN(2,i)
 77      format(T5,'SIMVAR_CUT on ', A20,' : ', F10.4,' to ', F10.4);
    ENDDO

    RETURN
  END SUBROUTINE PARSE_SIMVAR_CUTS


! ================================
    SUBROUTINE EVAL_SIMVAR_CUTS()

! Created Dec 2018
! Evaluate PASS_SIMCUTS based on user-input SIMVAR_CUTWIN_STRING
! Check hard-coded SIM_xxx names, and also check SIMSED_xxx names.
! 
! Nov 25 2019: check NEP_SIM_MAGOBS = number os SIM_MAGOBS < 99
! Jun 23 2020: abort of VARNAME not found
! Mar 07 2025: check SIM_PEAKMAG_[band]


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    LOGICAL FOUND
    INTEGER i, ipar, LENTMP, ifilt, ifilt_obs
    REAL*8  VMIN, VMAX
    CHARACTER VARNAME*40, CTMP*60, VARNAME_TMP*60, FNAM*20, cfilt*1

! --------------- BEGIN ------------

    IF ( NSIMVAR_CUTWIN == 0 ) RETURN
    FNAM = 'EVAL_SIMVAR_CUTS'

    DO 100 i = 1, NSIMVAR_CUTWIN
       FOUND   = .FALSE.
       VMIN    = DBLE( SIMVAR_CUTWIN(1,i) )
       VMAX    = DBLE( SIMVAR_CUTWIN(2,i) )
       VARNAME = SIMVAR_CUTWIN_LIST(i)
       if ( VARNAME .EQ. 'SIM_ZCMB' ) then
          FOUND   = .TRUE.
          if ( SIM_REDSHIFT_CMB < VMIN ) PASS_SIMCUTS = .FALSE.
          if ( SIM_REDSHIFT_CMB > VMAX ) PASS_SIMCUTS = .FALSE.

       else if ( VARNAME .EQ. 'SIM_PEAKMJD' ) then
          FOUND   = .TRUE.
          if ( SIM_PEAKMJD < VMIN ) PASS_SIMCUTS = .FALSE.
          if ( SIM_PEAKMJD > VMAX ) PASS_SIMCUTS = .FALSE.

       else if ( VARNAME .EQ. 'SIM_SALT2c' ) then
          FOUND   = .TRUE.
          if ( SIM_COLORPAR < VMIN ) PASS_SIMCUTS = .FALSE.
          if ( SIM_COLORPAR > VMAX ) PASS_SIMCUTS = .FALSE.

       else if ( VARNAME .EQ. 'SIM_SALT2x1' ) then
          FOUND   = .TRUE.
          if ( SIM_SHAPEPAR < VMIN ) PASS_SIMCUTS = .FALSE.
          if ( SIM_SHAPEPAR > VMAX ) PASS_SIMCUTS = .FALSE.

       else if ( VARNAME .EQ. 'NEP_SIM_MAGOBS' ) then
          FOUND   = .TRUE.
          if ( NEP_SIM_MAGOBS < int(VMIN) ) PASS_SIMCUTS = .FALSE.
          if ( NEP_SIM_MAGOBS > int(VMAX) ) PASS_SIMCUTS = .FALSE.
       endif

! check SIMSED variables
       DO IPAR    = 1, NPAR_SIMSED
          CTMP    = SIMSED_PARNAME(ipar)
          LENTMP  = INDEX(CTMP,' ') - 1
          VARNAME_TMP = 'SIMSED_' // CTMP(1:LENTMP)
          if ( VARNAME .EQ. VARNAME_TMP ) then
          FOUND   = .TRUE.
             if ( SIMSED_PARVAL(ipar) < VMIN ) PASS_SIMCUTS = .FALSE.
             if ( SIMSED_PARVAL(ipar) > VMAX ) PASS_SIMCUTS = .FALSE.
          endif
       ENDDO

!  Mar 2025: check PEAKMAG
       do ifilt = 1, NFILTDEF_SURVEY
          ifilt_obs   = IFILTDEF_MAP_SURVEY(ifilt)
          cfilt       = FILTDEF_STRING(ifilt_obs:ifilt_obs)
          VARNAME_TMP = 'SIM_PEAKMAG_' // cfilt
          if ( VARNAME .EQ. VARNAME_TMP ) then
             FOUND   = .TRUE.
             if ( SIM_PEAKMAG(ifilt) < VMIN ) PASS_SIMCUTS = .FALSE.
             if ( SIM_PEAKMAG(ifilt) > VMAX ) PASS_SIMCUTS = .FALSE.
          endif
       enddo


       if ( .NOT. FOUND ) then
          c1err = 'Could not evaluate cut on ' // VARNAME(1:20)
          c2err = 'Check &SNLCINP  SIMVAR_CUTWIN_STRING'
          CALL MADABORT(FNAM, c1err, c2err )
       endif
 100  CONTINUE


    RETURN
  END SUBROUTINE EVAL_SIMVAR_CUTS

! ======================================
    SUBROUTINE EVAL_PRIVATEVAR_CUTS()

! 
! Created Nov 3 2014
! Evaluate optional cuts on private variables.
! PASS_PRIVCUTS has been initialized to TRUE; set this
! flag to FALSE if any private variable cut-window fails.
! The private cuts are defined by &SNLCINP variable

! PRIVATE_CUTWIN_STRING = 'VARNAME  CUTMIN  CUTMAX'
!   or
! PRIVATE_CUTWIN_STRING = 'VARNAME  != VETO_VALUE'
! 
! Apr 27 2021: refactor to allow veto option
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE PRIVCOM

    IMPLICIT NONE

    INTEGER ivar, APPLY_CUTFLAG
    REAL*8 VAL, CUTMIN, CUTMAX
! -------------- BEGIN ------------

    IF ( NCUT_PRIVATE .EQ. 0 ) RETURN

    DO 100 ivar = 1, NVAR_PRIVATE
       APPLY_CUTFLAG = USE_PRIVATE_CUTWIN(ivar)
       if ( APPLY_CUTFLAG == 0 ) goto 100

       ! strip off info stored from init stage
       VAL    = PRIVATE_VALUE(ivar)
       CUTMIN = PRIVATE_CUTWIN(1,IVAR)
       CUTMAX = PRIVATE_CUTWIN(2,IVAR)

       IF ( APPLY_CUTFLAG < 0 ) then
          ! veto cut (Apr 2021)
          if ( VAL == CUTMIN ) PASS_PRIVCUTS = .FALSE.
       ELSE
          ! regular CUT-window
          if(VAL<CUTMIN .or. VAL>CUTMAX ) PASS_PRIVCUTS = .FALSE.
       ENDIF

       if ( .not. PASS_PRIVCUTS ) RETURN

 100  CONTINUE

    RETURN
  END SUBROUTINE EVAL_PRIVATEVAR_CUTS

! ===============================================
    SUBROUTINE INIT_FUDGE_FLUXCAL(IERR)

! 
! Init fudges on FLUXCAL[ERR]
! Apr 2017: call INIT_NONLIN().
! Dec 2017: use IGNOREFILE function on FUDGE_HOSTNOISE_FILE
!              to allow NONE
! Feb 15 2018: call INIT_FLUXERRMODEL
! Apr 15 2018: mass OPTMASK argument to INIT_FLUXERRMODEL


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IERR  ! (I) return 0 if all ok

    INTEGER LL, OPTMASK, MSKTMP, RDMASK
    CHARACTER cFILE*(MXCHAR_FILENAME), cNONE*6, FNAM*20
    LOGICAL FORCE_FLUXERRCOR, FORCE_FLUXCOR, VALID_FILE
    LOGICAL ALLOW_FLUXERRCOR, ALLOW_FLUXCOR, ISDATA

! function in sntools.c
    INTEGER IGNOREFILE
! --------------- BEGIN ----------------

    IERR = 0
    FNAM = 'INIT_FUDGE_FLUXCAL'

    DOFUDGE_HOSTNOISE    = .FALSE.
    DOFUDGE_NONLIN       = .FALSE.
    DOFUDGE_FLUXERRMODEL = .FALSE.
    FORCEMASK_FLUXCOR    = 0
    ISDATA = (.NOT. LSIM_SNANA)

! --------------------------------------------------------
! check FORCE option in case flux fudges already applied

    RDMASK = ISNLC_RDMASK_FLUXCOR_SNANA

    MSKTMP = MASK_FLUXERRCOR_SNANA
    CALL SETFORCE_FLUXCOR_SNANA(MSKTMP,FLUXERRMODEL_FILE)
    CALL SETFORCE_FLUXCOR_SNANA(MSKTMP,FUDGE_HOSTNOISE_FILE)
    FORCE_FLUXERRCOR = ( IAND(MSKTMP,FORCEMASK_FLUXCOR) > 0 )
    ALLOW_FLUXERRCOR =  & 
        FORCE_FLUXERRCOR .or. (IAND(RDMASK,MSKTMP)==0 )

    MSKTMP = MASK_FLUXCOR_SNANA
    CALL SETFORCE_FLUXCOR_SNANA(MSKTMP,MAGCOR_FILE)
    CALL SETFORCE_FLUXCOR_SNANA(MSKTMP,NONLINEARITY_FILE)
    FORCE_FLUXCOR  = ( IAND(MSKTMP,FORCEMASK_FLUXCOR) > 0 )
    ALLOW_FLUXCOR =  & 
        FORCE_FLUXCOR .or. (IAND(RDMASK,MSKTMP)==0 )

! ------------------------------------------------
    CALL ENVreplace(FUDGE_HOSTNOISE_FILE)
    LL    = INDEX(FUDGE_HOSTNOISE_FILE,' ')-1
    cFILE = FUDGE_HOSTNOISE_FILE(1:LL) // char(0)
    IF ( IGNOREFILE(cFILE,LL) == 0 ) THEN
       DOFUDGE_HOSTNOISE = .TRUE.
       CALL INIT_NOISEMODEL_HOST_LEGACY( cFILE, LL );
    ENDIF

! Apr 2017: check nonlinearity from user-supplied map
    CALL ENVreplace(NONLINEARITY_FILE);
    LL    = INDEX(NONLINEARITY_FILE,' ')-1
    cFILE = NONLINEARITY_FILE(1:LL) // char(0)
    if ( IGNOREFILE(cFILE,LL) == 0 ) then
       DOFUDGE_NONLIN       = .TRUE.
       CALL INIT_NONLIN(cFile,LL)
    endif

! ------------------------------
! Feb 2018: check FLUXERRMAPs for DATA

    OPTMASK = FLUXERRMODEL_OPTMASK
    if ( REQUIRE_DOCANA .NE. 0 ) then
!   note that MASK_REQUIRE_DOCANA_FLUXERRMAP=1024
!   in sntools_fluxErrModels.h
       FLUXERRMODEL_OPTMASK = FLUXERRMODEL_OPTMASK + 1024
    endif

    CALL ENVreplace(FLUXERRMODEL_FILE);
    LL    = INDEX(FLUXERRMODEL_FILE,' ')-1
    cFILE = FLUXERRMODEL_FILE(1:LL) // char(0)
    VALID_FILE = (IGNOREFILE(cFILE,LL) == 0)
    if ( VALID_FILE ) then
      CALL PRINTMSG_FLUXCOR(MASK_FLUXERRCOR_SNANA,"FLUXERRMODEL_FILE")
    endif
    IF ( VALID_FILE .and. ALLOW_FLUXERRCOR .and. ISDATA ) THEN
       cNONE = 'NONE' // char(0)
       DOFUDGE_FLUXERRMODEL = .TRUE.
       CALL INIT_FLUXERRMODEL(OPTMASK,cFile,cNONE,cNONE,LL,4,4)
    ENDIF

! March 2018: check FLUXERRMAPs for SIM
    CALL ENVreplace(SIM_FLUXERRMODEL_FILE);
    LL    = INDEX(SIM_FLUXERRMODEL_FILE,' ')-1
    cFILE = SIM_FLUXERRMODEL_FILE(1:LL) // char(0)
    IF ( IGNOREFILE(cFILE,LL) == 0  .and. LSIM_SNANA ) THEN
       cNONE = 'NONE' // char(0)
       DOFUDGE_FLUXERRMODEL = .TRUE.
       CALL INIT_FLUXERRMODEL(OPTMASK,cFile,cNONE,cNONE,LL,4,4)
    ENDIF

! ------------------------------------------------
! check MAGCOR (e.g. chromatic corrections) for data
! CALL INIT_MAGCOR()

    CALL ENVreplace(MAGCOR_FILE);
    LL    = INDEX(MAGCOR_FILE,' ')-1
    cFILE = MAGCOR_FILE(1:LL) // char(0)
    VALID_FILE = (IGNOREFILE(cFILE,LL) == 0)
    if ( VALID_FILE ) then
      CALL PRINTMSG_FLUXCOR(MASK_FLUXCOR_SNANA,"MAGCOR_FILE")
    endif

    IF ( VALID_FILE .and. ALLOW_FLUXCOR .and. ISDATA ) THEN
       CALL INIT_MAGCOR(MAGCOR_FILE);
    ENDIF

! check SIM_MAGCOR_FILE
    CALL ENVreplace(SIM_MAGCOR_FILE);
    LL    = INDEX(SIM_MAGCOR_FILE,' ')-1
    cFILE = SIM_MAGCOR_FILE(1:LL) // char(0)
    IF ( IGNOREFILE(cFILE,LL) == 0  .and. LSIM_SNANA ) THEN
       CALL INIT_MAGCOR(SIM_MAGCOR_FILE);
    ENDIF

! - - - - - - - - - - - - -
! error checking
    if ( DOFUDGE_HOSTNOISE .and. DOFUDGE_FLUXERRMODEL ) then
       c1err = 'Cannot specify FUDGE_HOSTNOISE_FILE ' //  & 
                 'and FLUXERRMODEL_FILE'
       c2err = 'Use one file or the other.'
       CALL MADABORT(FNAM, c1err, c2err )
    endif

    RETURN
  END SUBROUTINE INIT_FUDGE_FLUXCAL


! ==============================
    SUBROUTINE PRINTMSG_FLUXCOR(INMASK,CORFILE)

! Created Nov 14 2018
! Print appropriate message for FLUXCOR or FLUXERR correction.
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER  INMASK   ! (I) 1=fluxcor, 2=fluxerrcor

    INTEGER  RDMASK, LEN
    LOGICAL  ALLOW, FORCE, ISDATA
    CHARACTER*(*) CORFILE

! ------------------- BEGIN ------------------

    ISDATA = (.NOT. LSIM_SNANA)
    RDMASK = ISNLC_RDMASK_FLUXCOR_SNANA
    LEN    = INDEX(CORFILE,' ') - 1

    FORCE = ( IAND(INMASK,FORCEMASK_FLUXCOR)  > 0   )
    ALLOW = ( FORCE .or. ( IAND(INMASK,RDMASK)==0 ) )

! start with default behavior of ignoring already-applied correction
    IF ( (.NOT. ALLOW) .and. ISDATA ) THEN
       print*,' '
       print*,' !!! '//CORFILE// ' CORRECTION ALREADY APPLIED '
       print*,' !!!  ==> IGNORE '// CORFILE
       print*,' '
    ENDIF

! check option to FORCE correction, even if already applied
    IF ( ALLOW .and. FORCE .and. ISDATA ) THEN
       print*,' '
       print*,' !!! '//CORFILE// ' CORRECTION ALREADY APPLIED '
       print*,' !!! ==> FORCE '//CORFILE//' ANYWAY   '
       print*,' '
    ENDIF

! for sims, give reminder that corrections are not applied to sims.
    IF ( LSIM_SNANA ) THEN
      print*,' '
      print*,'  !!! IGNORE '//CORFILE//' for SIMULATION'
      print*,'  !!!   (or switch to SIM_'//CORFILE//' option)  '
      print*,' '
    ENDIF

    call flush(6)

    RETURN
  END SUBROUTINE PRINTMSG_FLUXCOR


! ===============================
    SUBROUTINE EXEC_FUDGE_FLUXCAL(EP)

! Created Nov 15, 2011 by R.Kessler and J.Marriner
! Fudge FLUXCAL and ERROR based on user namelist option
! FUDGE_FLUXCAL_OFFSET/ERROR
! 
! The following arrays are modified
!   - SNLC_FLUXCAL
!   - SNLC_FLUXCAL_ERRTOT
! 
! Note that the SNLC_MAG array is not modified !!!
! 
! 
! Mar 9 2020:
!   fix awful bug where previous FLUXCAL was used to apply
!   FUDGE_MAG_ERROR. Note RESTORE_DES3YR restores bug.
! 
! Dec 2 2022: add string-length args to get_FLUXERRMODEL call.
! 
! Jun 11 2024: implement SIM_FUDGE_MAG_ERROR_FILT (for DES-SN5YR)
! Jul 30 2024: move 'SNR = FLUXCAL/ERR1' aft3er ERR1 is computed to avoid divide-by-zero
! ------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER EP  ! (I) epoch index and SN index

! local var

    REAL*8  & 
         ERR1, ERR2a, ERR2b, ERR3, ERR3_SIM, ERR_SIM, ERR_DATA  & 
        ,ERR_ORIG, SCALE, FUDGE, SNR, SNR_PROTECT, LOGSNR  & 
        ,SQERR, SQERR1, SQERR2a, SQERR2b, SQERR3, SQERRHOST, SQERRMAP  & 
        ,RDNOISE_pe, FLUXCAL, FLUXCAL_BUG, FLUXADU  & 
        ,AREA, NEA, PSFSIG1, PSFSIG2, PSFRATIO  & 
        ,GALMAG, SBMAG, SBFLUX, SNSEP, NOISEPAR(2)  & 
        ,MJD, GAIN, SKYSIG, ZP, ZPDIF, PARLIST(20), Texpose  & 
        ,SQTMP, Fscale, Fpe_source, Fpe_sky, Fpe_galaxy, GENMAG  & 
        ,Fpe_list(3)

    INTEGER*8 GALID
    INTEGER IFILT_OBS, IFILT, OPT, NPAR, LL
    CHARACTER cfilt*2, FNAM*20, cFIELD*32, cBAND*2
    LOGICAL DOFUDGE

! function
    EXTERNAL NoiseEquivAperture, GET_NONLIN, get_FLUXERRMODEL
    REAL*8   NoiseEquivAperture, GET_NONLIN

! ------------ BEGIN ----------------

    FNAM = 'EXEC_FUDGE_FLUXCAL'
    IFILT_OBS = ISNLC_IFILT_OBS(ep)
    IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
    cfilt     = FILTDEF_STRING(ifilt_obs:ifilt_obs)

    SNLC_FLUXCAL(ep) = SNLC_FLUXCAL(ep)  & 
           + SNLC_FLUXCAL_OFF(ifilt)                &  ! offset in data file
           + FUDGE_FLUXCAL_OFFSET_FILT(IFILT_OBS)  ! user-nml offset

! get C function arguments
    cBAND        = cfilt(1:1) // char(0)
    LL           = INDEX(SNLC_FIELD(ep), ' ') - 1
    cFIELD       = SNLC_FIELD(ep)(1:LL) // char(0)
    GALID        = SNHOST_OBJID(1)
    SNSEP      = DBLE( SNHOST_ANGSEP(1) )
    GALMAG     = DBLE( SNHOST_MAGOBS(ifilt,1) )
    SBFLUX     = DBLE( SNHOST_SBFLUXCAL(ifilt) )
    SBMAG      = DBLE( SNHOST_SBMAG(ifilt) )
    FLUXCAL    = DBLE( SNLC_FLUXCAL(ep) )
    PSFSIG1    = DBLE( SNLC_PSF_SIG1(EP) )  ! Gauss sigma, pixel
    PSFSIG2    = DBLE( SNLC_PSF_SIG2(EP) )  ! idem, 2nd component
    PSFRATIO   = DBLE( SNLC_PSF_RATIO(EP) )
    NEA        = DBLE( SNLC_PSF_NEA(EP) )
    GAIN       = DBLE( SNLC_GAIN(ep) )
    SKYSIG     = DBLE( SNLC_SKYSIG(ep) )          ! sigma, ADU/pixel
    RDNOISE_pe = DBLE(SNLC_RDNOISE(ep) )  ! e-/pix
    ZP         = DBLE( SNLC_ZEROPT(ep) )
    ZPDIF      = ZP - ZP_FLUXCAL
    MJD        = SNLC8_MJD(EP)

! - - - - - - - - - - - - - -
    SCALE    = DBLE(FUDGE_FLUXERR_SCALE_FILT(IFILT_OBS))
    ERR_ORIG = DBLE(SNLC_FLUXCAL_ERRTOT(ep))
    ERR1  = ERR_ORIG * SCALE
    SNR   = FLUXCAL/ERR1

    FUDGE = FUDGE_FLUXCAL_ERROR_FILT(IFILT_OBS)
    ERR2a = DBLE(FUDGE)  ! fudge error

    FUDGE = FUDGE_FLUXCAL_ERRPIX_FILT(IFILT_OBS)
    ERR2b = DBLE(FUDGE) ! fudge err/pix

    FUDGE = FUDGE_MAG_ERROR_FILT(IFILT_OBS)
    ERR3  = DBLE(FUDGE)* FLUXCAL ! fudge mag err

    SQERR1     = ERR1  * ERR1
    SQERR2a    = ERR2a * ERR2a
    IF ( ERR2a .LT. 0.0 ) SQERR2a = -SQERR2a
    SQERR2b    = ERR2b * ERR2b   ! * AREA  ! sum in quad over pixels
    SQERR3     = ERR3  * ERR3
    SQERRHOST  = 0.0
    SQERRMAP   = 0.0

    IF ( LSIM_SNANA ) THEN
      FUDGE     = SIM_FUDGE_MAG_ERROR_FILT(IFILT_OBS)
      ERR3_SIM  = DBLE(FUDGE)* FLUXCAL ! fudge mag err
	SQERR3    = SQERR3 + (ERR3_SIM*ERR3_SIM)
    ENDIF

! --------------------------------------------
!  Check for map corrections

    IF ( DOFUDGE_HOSTNOISE  ) THEN  ! compute ERRHOST
! call function to return
!  noisePar(1) =  FLUXCAL noise per pixel
!  noisePar(2) =  ERRSCALE to multiply skynoise
       CALL GEN_NOISEMODEL_HOST_LEGACY(cBAND,cFIELD,GALID,  & 
                 GALMAG, SBMAG, SNSEP, noisePar, 2,30);
       SQTMP     = NOISEPAR(2)**2 - 1.0 ! don't double-count error
       SQERRHOST = SQERR1 * SQTMP
    ENDIF


    IF ( DOFUDGE_FLUXERRMODEL  ) THEN ! compute ERRMAP
      OPT = 0 ; NPAR=0
      SNR_PROTECT = MAX(SNR,0.126)    ! --> LOGSNR >= -0.9
      LOGSNR = log10(SNR_PROTECT)

!  grep IPAR_FLUXERRMAP sntools_fluxErrModels.h
      NPAR = NPAR+1;  PARLIST(NPAR) = MJD
      NPAR = NPAR+1;  PARLIST(NPAR) = PSFSIG1
      NPAR = NPAR+1;  PARLIST(NPAR) = SKYSIG
      NPAR = NPAR+1;  PARLIST(NPAR) = ZP
      NPAR = NPAR+1;  PARLIST(NPAR) = LOGSNR
      NPAR = NPAR+1;  PARLIST(NPAR) = SBMAG
      NPAR = NPAR+1;  PARLIST(NPAR) = GALMAG
      NPAR = NPAR+1;  PARLIST(NPAR) = SNSEP
      CALL  get_FLUXERRMODEL(OPT,ERR1,cBAND,cFIELD,NPAR, PARLIST,  & 
                  ERR_SIM, ERR_DATA, 2,30 )  ! <== return ERR_SIM, ERR_DATA

      SQERRMAP = (ERR_DATA*ERR_DATA) - SQERR1 ! note: can be negative
    ENDIF


! avoid invalid fluxes.
    IF ( FLUXCAL .EQ. -9.0 ) RETURN

! bail if nothing to do
    DOFUDGE =  & 
           ( ERR1  > 0.0 ) .or.  & 
           ( ERR2a > 0.0 ) .or.  & 
           ( ERR2b > 0.0 ) .or.  & 
           ( ERR3  > 0.0 ) .or.  & 
           ( SQERRHOST  > 0.0 ) .or.  & 
           ( SQERRMAP   > 0.0 )

    if ( .NOT. DOFUDGE ) GOTO 880  ! jump to nonlinearity flux-fudge

    AREA = 0.0
    IF ( ERR2b > 0.0 ) THEN
       ! sum noise over pixels in effective background area
       if ( PSFSIG1 < 0.002 ) then
          c1err = 'Undefined PSF for FUDGE_FLUXCAL_ERRPIX'
          write(c2err,62) SNLC_CCID, MJD, cfilt
 62         format('CID=',A10,' for MJD=', F9.3,'-', A)
          CALL MADABORT(FNAM, c1err, c2err )
       endif
       AREA  = NEA
    ENDIF


    SQERR = SQERR1 + SQERR2a + (SQERR2b*AREA) + SQERR3  & 
            + SQERRHOST + SQERRMAP

    IF ( SQERR .LE. 0.0 ) THEN
       CALL PRINT_PREABORT_BANNER(FNAM(1:18)//char(0),40)
       print*,'   ERR1(from data file): ',          sngl(ERR1)
       print*,'   ERR2a(user-fluxcal fudge): ',     sngl(ERR2a)
       print*,'   ERR2b(user-fluxcal fudge/pix): ', sngl(ERR2b)
       print*,'   ERR3(user magerr fudge): ' ,      sngl(ERR3)
       print*,'   ERRHOST(FUDGE_HOSTNOISE_FILE): ',  & 
                        sngl( sqrt(SQERRHOST) )
       print*,'   ERRMAP(FLUXERRMODEL_FILE): ',  & 
                        sngl( sqrt(SQERRMAP) )
       write(c1err,61)  & 
              SNLC_CCID(1:ISNLC_LENCCID), cfilt, MJD
 61      format('Negative flux error for CID=',A, '-',A,  & 
               3x, 'and MJD=',F9.3)
       c2err = 'Check error fudges in &SNLCINP'
       CALL MADABORT(FNAM, c1err, c2err )
    ENDIF

    SNLC_FLUXCAL_ERRTOT(ep) = SNGL ( SQRT(SQERR) )

    CALL SETMASK_FLUXCOR_SNANA(MASK_FLUXERRCOR_SNANA)

! ---------------------------
880   CONTINUE
! ---------------------------

! check GET_NONLIN option in sntools_nonlinearity.c
    IF ( DOFUDGE_NONLIN ) THEN

      if ( GAIN < 0.01 .or. SKYSIG < 0.0001 .or. ZP < 0.01 ) then
         CALL PRINT_PREABORT_BANNER(FNAM(1:18)//char(0),40)
         print*,'    BAND    = ', CFILT
         print*,'    GAIN    = ', sngl(GAIN), ' e-/pix'
         print*,'    ZP      = ', sngl(ZP),'  (ADU->MAG)'
         print*,'    SKYSIG  = ', sngl(SKYSIG),' ADU/pix '
         print*,'    RDNOISE = ', sngl(RDNOISE_pe),' e-/pix'
         c1err = 'Canot compute NON-LINEARITY'
         c2err = 'See above values'
         CALL MADABORT(FNAM, c1err, c2err )
      endif

    ! sum noise over pixels in effective background area
      FLUXADU = FLUXCAL * 10.0**(0.4*ZPDIF)

      Texpose    = SNLC_TEXPOSE(ep)
      Fpe_source = FLUXADU * GAIN
      IF ( Fpe_source < 0.0 ) Fpe_source = 0.01 ! must be positive
      Fpe_sky    = NEA * (SKYSIG*GAIN)**2
      Fpe_galaxy = 0.0  ! ?? NEA * SBMAG * GAIN ! ??
      GENMAG     = DBLE( SNLC_MAG(ep) )
      Fpe_list(1) = Fpe_source
      Fpe_list(2) = Fpe_sky
      Fpe_list(3) = Fpe_galaxy

      Fscale = GET_NONLIN(SNLC_CCID(1:ISNLC_LENCCID)//char(0),  & 
                        cfilt(1:1)//char(0), Texpose, NEA,  & 
                        Fpe_list, GENMAG, ISNLC_LENCCID, 2 )
      SNLC_FLUXCAL(ep) = SNLC_FLUXCAL(ep) * SNGL(Fscale)
      CALL SETMASK_FLUXCOR_SNANA(MASK_FLUXCOR_SNANA)
    ENDIF

! ---------------------
! mag-corrections; e.g.. chromatic corrections from FGCM
    CALL EXEC_MAGCOR(ep)

    RETURN
  END SUBROUTINE EXEC_FUDGE_FLUXCAL

! ===========================================
    SUBROUTINE SETFORCE_FLUXCOR_SNANA(MASK,INFILE)

! Created Nov 14 2018
! INFILE = '<fileName>'       --> do nothing
! INFILE = '<fileName FORCE>' --> set MASK bit of FORCEMASK_FLUXCOR
!                                 and set INFILE = fileName
! - - - - - -

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER MASK         ! (I) mask to set for FORCEMASK
    CHARACTER INFILE*(*) ! (I) name of input file

    CHARACTER NAME_forC*(MXCHAR_FILENAME), cwd*(MXCHAR_FILENAME), FNAM*30
    INTEGER iwd, NWD, LENFILE,MSKOPT

    LOGICAL  IGNOREFILE_fortran
    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! ---------------- BEGIN --------------

    if ( IGNOREFILE_fortran(INFILE) ) RETURN  ! Apr 7 2022

    FNAM      = 'SETFORCE_FLUXCOR_SNANA'
    MSKOPT    = MSKOPT_PARSE_WORDS_STRING
    LENFILE   = INDEX(INFILE,' ') + 10 ! allow FORCE key after space
    NAME_forC = INFILE(1:LENFILE-1)//char(0)
    NWD = STORE_PARSE_WORDS(MSKOPT, NAME_forC,  & 
               FNAM//char(0), LENFILE, 30)

    IF ( NWD <= 1 ) RETURN

! with more than 1 word, check for FORCE key

    DO 10 iwd = 1, NWD
       CALL get_PARSE_WORD_fortran(iwd,cwd,LENFILE)

       if ( iwd==1 ) then
          INFILE = cwd
       else if ( cwd(1:5) .EQ. 'FORCE' ) then
          FORCEMASK_FLUXCOR = IOR(FORCEMASK_FLUXCOR,MASK)
       endif
10    CONTINUE


    RETURN
  END SUBROUTINE SETFORCE_FLUXCOR_SNANA

! ===========================================
    SUBROUTINE SETMASK_FLUXCOR_SNANA(MASK)


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER MASK  ! (I) mask to set

    INTEGER ALREADY_SET, FORCE
    CHARACTER FNAM*24

! --------------- BEGIN ------------------

    FNAM = 'SETMASK_FLUXCOR_SNANA'

! first check mask read from data file to see if fluxcor fudge
! was already applied

    ALREADY_SET = IAND(ISNLC_RDMASK_FLUXCOR_SNANA,MASK)
    FORCE       = IAND(FORCEMASK_FLUXCOR,MASK)

    if ( ALREADY_SET > 0 .and. FORCE == 0 ) RETURN

! xxxxxxxxxxxxxxxx
!      if ( ALREADY_SET > 0 .and. FORCE == 0 ) then
!         CMASK(MASK_FLUXCOR_SNANA)    = 'FLUXCOR/MAGCOR'
!         CMASK(MASK_FLUXERRCOR_SNANA) = 'FLUXERRCOR'
!         C1ERR = 'Will not apply ' // CMASK(MASK)
!         C2ERR = 'because data files already include '//CMASK(MASK)
!         CALL MADABORT(FNAM, c1err, c2err )
!      endif
! xxxxxxxxxxxxxxxx

    ISNLC_WRMASK_FLUXCOR_SNANA =  & 
           IOR(ISNLC_WRMASK_FLUXCOR_SNANA,MASK)

    RETURN
  END SUBROUTINE SETMASK_FLUXCOR_SNANA

! =================================
    SUBROUTINE COUNTFIELDS()
! 
! Jun 8, 2009 R.Kessler
! 
! Fill ISNLC_NFIELD with number of different fields used by this SN.
! NFIELD = 1 for isolated field, but can be > 1 for overlapping fields.
! 
! 
! --------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER  & 
         ep, NEP, i, ovp  & 
        ,IDF, IDF_OVPLIST(MXFIELD_OVP)  & 
        ,NFIELD, LENTMP, LENTOT, L2

    CHARACTER FIELDTMP*(MXCHAR_FIELDNAME), FNAM*12
    LOGICAL USE, FIRST

! ------------- BEGIN ----------

    NEP    = ISNLC_NEPOCH_STORE
    NFIELD = 0
    FNAM   = 'COUNTFIELDS'

    DO 101 ep = 1, NEP
        IDF   = ISNLC_IDFIELD(ep)

        USE = .FALSE.
        DO 102 i = 1, MIN(NFIELD,MXFIELD_OVP)
          if ( SNLC_FIELD(ep).EQ.SNLC_FIELD_OVPLIST(i) ) USE=.TRUE.
 102      CONTINUE

        if ( .NOT. USE ) THEN
           NFIELD = NFIELD + 1
           IF ( NFIELD .LE. MXFIELD_OVP ) then
              SNLC_FIELD_OVPLIST(NFIELD) = SNLC_FIELD(ep)
              IDF_OVPLIST(NFIELD) = IDF
           ENDIF
        endif  ! end of USE
101   CONTINUE  ! ep loop


    IF ( NFIELD > MXFIELD_OVP ) THEN
       CALL PRINT_PREABORT_BANNER(FNAM(1:11)//char(0),20)

       DO  i = 1, MIN(NFIELD,MXFIELD_OVP)
          write(6,60) i, SNLC_FIELD_OVPLIST(i)
60          format(T5,'FIELD_OVPLIST(',I2,') = ', A)
       ENDDO

       write(c1err,61) NFIELD
       write(c2err,62) MXFIELD_OVP
 61      format('NFIELD_OVP = ', I3,' exceeds bound of')
 62      format('MXFIELD_OVP = ', I3 )
       CALL MADABORT(FNAM, c1err, c2err )
    ENDIF

    ISNLC_NFIELD_OVP = NFIELD

! prepare catenated string SNLC_FIELD_forC that is sorted based
! on survey fields.  e.g., if SURVEY_FIELDS = 'E1', 'E2' then
! SNLC_FIELDS = E2,E1 gets constructed as SNLC_FIELD_forC = 'E1+E2'

    SNLC_FIELDLIST = ''
    LENTOT = 0
    FIRST  = .TRUE.

    DO 300 i = 1, NFIELD_SURVEY
      FIELDTMP = SURVEY_FIELDNAME(i)  ! from SURVEY.DEF file
      LENTMP   = INDEX(FIELDTMP,' ')-1

      DO 301 ovp = 1, NFIELD  ! overlap fields in this SN
        L2 = INDEX(SNLC_FIELD_OVPLIST(ovp),' ')-1
        IF ( LENTMP .NE. L2 ) GOTO 301

        IF ( SNLC_FIELD_OVPLIST(ovp) .EQ. FIELDTMP ) THEN

           if ( FIRST ) then
             SNLC_FIELDLIST = FIELDTMP
             LENTOT = LENTOT + LENTMP
             FIRST  = .FALSE.
           else
              SNLC_FIELDLIST =  & 
                SNLC_FIELDLIST(1:LENTOT) // '+' // FIELDTMP(1:LENTMP)
              LENTOT = LENTOT + LENTMP + 1
           endif

        ENDIF
 301    CONTINUE
 300  CONTINUE

    IF ( LENTOT .EQ. 0  ) THEN
       SNLC_FIELDLIST = 'VOID'
    ENDIF

    RETURN
  END SUBROUTINE COUNTFIELDS


! ========================================
    SUBROUTINE U3BAND()
! 
! Created Sep 13, 2008 by R.Kessler
! 
! Special flag to require exactly three passbands
! that pass the global SNRCUT, and one of the passbands
! must be U-band. If only two passbands exist, the SN
! is rejected. If 4 or more passbands exist, then the
! reddest extra passbands are erased and the SN is used.
! For example, if griz -> UBVR, then all z-band fluxes
! and errors are set to -9.0, and SNRECON is called
! again; this SN is then processed as if z-band had
! never existed.
! 
! 
! Jul 16, 2011: pass ISN as argument and remove IERR arg
! --------------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT

    IMPLICIT NONE

    INTEGER  & 
         NBAND, NTMP  & 
        ,IFILT, IFILT_OBS, IFILT_OBS_LIST(MXFILT_ALL)  & 
        ,ORDER  & 
        ,INDEX_SORT(MXFILT_ALL)  & 
        ,i, isort  & 
        ,NERASE

    REAL  & 
         ZZ, RESTLAM  & 
        ,RESTLAM_LIST(MXFILT_ALL)

    LOGICAL LSNRMAX, LU, LKEEP

    character cfilt*2, ccid*(MXCHAR_CCID), cdiscard*32

    REAL, PARAMETER :: ULAMMAX = 3900.0  ! max U-band lambda in rest-frame

    EXTERNAL SORTFLOAT

! ---------- BEGIN ---------

    NERASE = 0

    if ( .not. LTEST_U3BAND ) RETURN

    CCID = SNLC_CCID
    LKEEP = .TRUE.  ! keep SN by default

! require at least 3 passbands with SNRMAX cut.
! If not, set PEAKMJD = -9.0 and return so that
! this SN is rejected.

      NBAND = ISNLC_NFILT_SNRMAX
      IF ( NBAND .LT. 3 ) THEN
        LKEEP  = .FALSE.
        cdiscard = 'NBAND(SNRMAX) < 3'
        GOTO 150
      ENDIF

      ZZ          = 1.0  + SNLC_REDSHIFT


! filter loop: get list of observer passbands and check if
! U-band is there.

      NTMP = 0
      LU   = .FALSE.

      DO 101 ifilt = 1, NFILTDEF_SURVEY

         if ( .NOT. EXIST_FILT(ifilt) ) goto 101

         ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)

         RESTLAM = FILTOBS_LAMAVG(ifilt_obs) / ZZ
         LSNRMAX = SNLC_SNRMAX_FILT(ifilt) .GT. CUTWIN_SNRMAX(1)

         if ( .NOT. LSNRMAX ) GOTO 101

         IF ( RESTLAM .LT. ULAMMAX .and.  & 
                RESTLAM .GT. CUTWIN_RESTLAM(1) ) LU = .TRUE.

         NTMP = NTMP + 1
         IFILT_OBS_LIST(NTMP) = IFILT_OBS
         RESTLAM_LIST(NTMP)   = RESTLAM

101     CONTINUE  ! UFILT loop


      if ( .NOT. LU ) then
         LKEEP = .FALSE.
         cdiscard = 'No filter maps onto rest-U'
         goto 150
      endif

      if ( LU .and. NBAND .EQ. 3 ) GOTO 150

! here we have more than 3 bands with good SNRMAX.
! Erase reddest bands until we have the 3 bluest.

     ORDER = + 1 ! increasing order
     CALL SORTFLOAT( NBAND, RESTLAM_LIST, ORDER, INDEX_SORT)

! erase reddest filters as if they never existed
    DO i = 4, NTMP
      isort = INDEX_SORT(i)
      IFILT_OBS = IFILT_OBS_LIST(isort)
      cfilt = filtdef_string(ifilt_obs:ifilt_obs)
      RESTLAM   = RESTLAM_LIST(isort)

      write(6,440) CCID, cfilt, RESTLAM
440     format(T8,'U3BAND CID=',A8, ' : erase filter=',A,  & 
             ' with rest-Lambda=',F6.0,' A' )

      NERASE =  NERASE + 1
      CALL ERASE_FILTER(ifilt_obs)
    ENDDO

! --------
150     CONTINUE
      IF ( .NOT. LKEEP ) THEN
        SNLC_SEARCH_PEAKMJD = -9.0
        write(6,441) CCID, cdiscard
441       format(T8,'U3BAND CID=',A8, ' => DISCARD because ', A )
      ENDIF

! ------------------------
! call SNRECON again after erasing filter(s)
    IF ( NERASE .GT. 0 ) THEN
      CALL SNRECON()
    ENDIF


    RETURN
  END SUBROUTINE U3BAND


! ========================================
    SUBROUTINE ERASE_FILTER(ifilt_obs)
! 
! Created Sep 13, 2008 by R.Kessler
! 
! Erase "ifilt_obs" for SN index ISN.
! Erasing is done by setting fluxes and errors to -9.
! 
! --------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IFILT_OBS   !(O) SN index, and filter index

    INTEGER NEWMJD, EPMIN, EPMAX, ep, ifilt
    LOGICAL LTMP

! -------- BEGIN -------

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

      DO EP = EPMIN, EPMAX

        LTMP = IFILT_OBS .EQ. ISNLC_IFILT_OBS(ep)

        IF ( LTMP ) THEN
          SNLC_FLUXCAL(ep)        = -9.0
          SNLC_FLUXCAL_ERRTOT(ep) = -9.0
          SNLC_MAG(ep)        = -9.0
          SNLC_ZEROPT(ep)     = -9.0
          SNLC_ZEROPT_ERR(ep) = -9.0

          IFILT  = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
          EXIST_FILT(ifilt)         = .FALSE.
          ISNLC_NEPOCH_FILT(ifilt)   =  0
          SNLC_SNRMAX_FILT(0)       = -9.
          SNLC_SNRMAX_FILT(ifilt)   = -9.
          SNLC_SNRMAX_SORT(ifilt)   = -9.
          SNLC_FLUXCALMAX(ifilt)    = -9.
          SNLC_FLUXCALMAX_ERR(ifilt)  = -9.

        ENDIF

      ENDDO

200   CONTINUE

    RETURN
  END SUBROUTINE ERASE_FILTER

! =============================================
    DOUBLE PRECISION FUNCTION DLMAG8_REF(Z8)
! Oct 23 2020: refactor DLMAG function to use sntools_cosmology.c


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    REAL*8 Z8             ! input redshift
    REAL*8 zCMB, zHEL, vPEC, H0, OM, OL, w0, wa

    REAL*8   DLMAG_fortC    ! C function in sntools_cosmology.c
    EXTERNAL DLMAG_fortC

! ------- BEGIN -------

    H0   = H0_REF(1)
    OM   = OMAT_REF(1)
    OL   = OLAM_REF(1)
    w0   = W0_REF(1)
    wa   = WA_REF(1)
    zCMB = z8
    zHEL = z8
    vPEC = 0.0
    DLMAG8_REF = DLMAG_fortC(zCMB, zHEL, vPEC, H0, OM, OL, w0, wa)
    RETURN
  END FUNCTION DLMAG8_REF



! =============================
    INTEGER FUNCTION DETNUMSTAT(DETNUM)
! 
! Dec 2017
! Returns 0 if DETNUM is valid; returns -1 otherwise.
! -------------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! function args
    INTEGER DETNUM

! local args
    INTEGER i
    LOGICAL LDETNUM

! -------------- BEGIN --------------

    IF ( NDETNUM_LIST .LE. 0 ) THEN
       ! don't check  anything
       DETNUMSTAT = 0
       RETURN
    ENDIF

    LDETNUM = .FALSE.
    do i = 1, NDETNUM_LIST
      if ( DETNUM_LIST(i) .EQ. DETNUM ) LDETNUM = .TRUE.
    enddo


    IF ( LDETNUM ) THEN
       DETNUMSTAT   = 0
    ELSE
       DETNUMSTAT = -1
    endif

    RETURN
  END FUNCTION DETNUMSTAT


! ================================================
    SUBROUTINE MAKE_SIMLIB_FILE(OPT_FLAG)

! Created Feb 2016
! Create LIBID entry in simlib file, based on info from data file.
! Beware that event distribution is not necessarily uniform since
! it has the data distribution.
! 
! OPT_FLAG = 1 --> open simlib file and write header
! OPT_FLAG = 2 --> add next LIBID
! OPT_FLAG = 3 --> close file
! 
!  OPT_SIMLIB_OUT & 1 -> write SIM_MAGOBS for all obs,
!                        even for true flux = 0 and mag=99.
!  OPT_SIMLIB_OUT & 2 -> write SIM_MAGOBS only for true flux > 0
! 
! 
! Oct 15, 2019: fix format for SKYSIG : 6.2f -> 7.2f
! Nov 22, 2019: adapt to work on fakes and sims, and write MAG column.
! Oct 05, 2020: add ZP requirement to compute metadata; see FOUND_METADATA
! Jan 19, 2021: add pad space between MJD and IDEXPT
! Mar 26, 2021:
!   + suppress GALID since it wrote 12345 (garbage) for DC2
!   + skip LIBID of PEAKMJD < 1000
! 
! --------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SIMLIBCOM
    USE SPECCOM

    IMPLICIT NONE

    INTEGER OPT_FLAG ! 1=> init file, 2=> add simlib entry

! local args
    REAL*8  MJD, GAIN, RDNOISE, SKYSIG, PSF1, PSF2, PSFRAT
    REAL*8  ZP, ZPERR, MAG, FLUX, FLUXERR, PIXSIZE, LAMOBS
    REAL*8  MJD_MIN, TOBS_MIN, PEAKMJD, DVAL
    INTEGER LIBID, IFILT, IFILT_OBS, LASTEP(MXFILT_OBS)
    INTEGER LEN0, LEN1, LEN2, EPMIN, EPMAX, NEWMJD, NOBS, iep, EP
    INTEGER LENcMAG, ISTAT, LENF, NFILTDEF_TMP
    LOGICAL IS_REAL_DATA, WRALL_SIM_MAGOBS, WRSET_SIM_MAGOBS
    LOGICAL WR_SIM_MAGOBS, FOUND_METADATA
    CHARACTER BAND*4, FNAM*20, SEDCMD*200, cDUM*20, CCID*32
    CHARACTER cMAG*12, cNCUTS*12, STR_OLD*40, STR_NEW*60

! ------------------ BEGIN ------------------

    IF( SIMLIB_OUTFILE .EQ. '' ) RETURN

    FNAM = 'MAKE_SIMLIB_FILE'
    IS_REAL_DATA = ( .NOT. ISJOB_SIM )

! check option to write ALL SIM_MAGOBS, even for SIM_MAGOBS=99 where
! no fake or sim event was generated and there is zero flux.
    WRALL_SIM_MAGOBS =  & 
          BTEST(OPT_SIMLIB_OUT,0) .and. (.not.IS_REAL_DATA)

! check option to write SIM_MAGOBS  only for true flux > 0
    WRSET_SIM_MAGOBS =  & 
          BTEST(OPT_SIMLIB_OUT,1) .and. (.not.IS_REAL_DATA)
    WR_SIM_MAGOBS = WRALL_SIM_MAGOBS .or.  WRSET_SIM_MAGOBS

    NFILTDEF_TMP = NFILTDEF_OBS

! -------------------------------------------

    IF ( OPT_FLAG == 1 ) THEN

#if defined(SNANA)
      if ( NFILTDEF_TMP == 0 ) then
         C1ERR = 'NFILTDEF_OBS=0 --> no defined obs-frame filters'
         C2ERR = 'Check &SNLCINP input CALIB_FILE'
         CALL MADABORT(FNAM, c1err, c2err)
      endif
#endif


#if defined(SNFIT) || defined(PSNID)
    C1ERR = 'Cannot use SIMLIB_OUT option with snlc_fit or psnid.'
    C2ERR = 'Use snana.exe instead. '
    CALL MADABORT(FNAM, c1err, c2err)
#endif

      GLOBAL_BANNER =  & 
            'Prepare to write SIMLIB file from data'
      CALL PRBANNER (GLOBAL_BANNER(1:60) )
      LEN1 = INDEX(SIMLIB_OUTFILE,' ') -1
      print*,'   Open output SIMLIB file: ', SIMLIB_OUTFILE(1:LEN1)

      OPEN(   UNIT   = LUNOUT  & 
              , FILE   = SIMLIB_OUTFILE  & 
              , STATUS = 'UNKNOWN'         )

! - - - - - -
! write DOCANA

      WRITE(LUNOUT,20)  'DOCUMENTATION:'

      WRITE(LUNOUT,24)  & 
             'PURPOSE: simulation cadence created from data'

      WRITE(LUNOUT,24)  & 
             'VERSION_PHOTOMETRY: ' // VERSION_PHOTOMETRY(1)(1:60)

      WRITE(LUNOUT,24) 'NOTES: '
      WRITE(LUNOUT,24)  & 
             '- Created from snana.exe using SIMLIB_OUTFILE option.'
      if ( WR_SIM_MAGOBS ) then
         WRITE(LUNOUT,24)  & 
              '- SIM_MAGOBS are copied to MAG column.'
      endif

      WRITE(LUNOUT,20) 'DOCUMENTATION_END:'
      WRITE(LUNOUT,20) ' '

! - - - - - - - -
      LEN1 = INDEX(SURVEY_NAME,' ') - 1
      LEN2 = INDEX(SURVEY_FILTERS,' ') - 1
      WRITE(LUNOUT,10) 'SURVEY:      ',  SURVEY_NAME(1:LEN1)
      WRITE(LUNOUT,10) 'FILTERS:     ',  SURVEY_FILTERS(1:LEN2)
10      FORMAT(A, 1x, A)  ! KEY VAL

20      FORMAT(A)
24      FORMAT(T4,A)

      WRITE(LUNOUT,10) 'NLIBID:     ',  'NLIBID_REPLACE'

      if ( IS_REAL_DATA ) then
        WRITE(LUNOUT,10) 'PSF_UNIT:    ',  'ARCSEC_FWHM'
        WRITE(LUNOUT,10) 'SKYSIG_UNIT: ',  'ADU_PER_SQARCSEC'
        CALL SET_SIMLIB_SKYMAG()
        CALL SET_SIMLIB_ZPERR()
     endif

      WRITE(LUNOUT, '(/,A,/)' ) 'BEGIN LIBGEN'
! 90      FORMAT(/, 'BEGIN LIBGEN',/  )

      CALL FLUSH(LUNOUT)
      CALL FLUSH(6)

      RETURN
    ENDIF

    IF ( OPT_FLAG == 3 ) THEN
       write(LUNOUT,400)
400      format(/, 'END_OF_SIMLIB: ', / )
       close(UNIT = LUNOUT)

! replace N_SNLC_CUTS with comment giving number of SIMLIB entries
       write(cNCUTS,'(I6,2x)') N_SNLC_CUTS
       STR_OLD = 'NLIBID_REPLACE '
       STR_NEW = cNCUTS
       LEN0    = index(STR_OLD,' ') - 1
       LEN1    = 8    ! for cNCUT
       LEN2    = INDEX(SIMLIB_OUTFILE,' ') - 1
       write(SEDCMD,505)  & 
             STR_OLD(1:LEN0), STR_NEW(1:LEN1), SIMLIB_OUTFILE(1:LEN2)
!        print*,' xxx SEDCMD = ', SEDCMD
       CALL SYSTEM(SEDCMD)

505      FORMAT('sed -i ', "'s/", A, '/',A, "/g1'", 2x, A)

! if SUBSURVEY_LIST is not blank, use system 'sed' command
! to insert its comma-separate list at top of SIMLIB file.
! fixed key SUBSURVEY_LIST to be printed after SURVEY key (M. Vincenzi Feb 2022)
       IF ( SUBSURVEY_NAME_LIST .NE. '' ) THEN
          LEN1   = INDEX(SUBSURVEY_NAME_LIST,' ') - 1
          LEN2   = INDEX(SIMLIB_OUTFILE,' ') - 1
          SEDCMD = "sed -i '/SURVEY:/a  SUBSURVEY_LIST: "  & 
                // SUBSURVEY_NAME_LIST(1:LEN1) // "' "  & 
                // SIMLIB_OUTFILE(1:LEN2)
          CALL SYSTEM(SEDCMD)
       ENDIF

       return
    ENDIF

! ------------- OPT_FLAG=2 below -------------------------

    IF ( ISNLC_NEWMJD_STORE < 4 ) RETURN
    IF ( SNLC_SEARCH_PEAKMJD < 1000.0 ) RETURN ! Mar 26 2021

    DO IFILT_OBS = 1, MXFILT_OBS
       LASTEP(IFILT_OBS) = -9
    ENDDO


! ----------
! first loop through and count NOBS
    NOBS = 0
    DO 25 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 21 EP = EPMIN, EPMAX
         if ( WR_SIM_MAGOBS ) then
            MAG = SIM_EPMAGOBS(ep)
            if ( WRSET_SIM_MAGOBS .and. MAG > 98.0 ) goto 21
         endif
         NOBS = NOBS + 1
21     CONTINUE
25    CONTINUE

! -----------
    PIXSIZE = SNLC_PIXSIZE
    if ( PIXSIZE < 0.0 ) PIXSIZE = PIXSIZE_GUESS

    WRITE(LUNOUT,66)
66    FORMAT('#--------------------------------------------' )

! check option to force TMIN
    MJD_MIN  = SNLC8_MJD(1)
    PEAKMJD  = SNLC_SEARCH_PEAKMJD
    IF ( SIMLIB_OUT_TMINFIX < 999 ) THEN
      PEAKMJD = MJD_MIN - DBLE(SIMLIB_OUT_TMINFIX)
    ENDIF
    TOBS_MIN = MJD_MIN - PEAKMJD

    LIBID = SNLC_CID      ! UNIQUE_CIDINT(SNLC_CCID)

    WRITE(LUNOUT,101) LIBID, SNLC_CCID(1:ISNLC_LENCCID), TOBS_MIN
 101  FORMAT('LIBID: ', I8, 5x, '# cadence from SNID=', A, 2x,  & 
          'min(MJD-PEAKMJD)=',F5.1  )

    WRITE(LUNOUT,103) SNLC8_RA, SNLC8_DEC, SNLC_MWEBV
 103  FORMAT('RA: ', F12.6, 8x, 'DEC:', F12.6, 5x, 'MWEBV: ',F7.4)

    WRITE(LUNOUT,105) NOBS, PIXSIZE,  & 
             SNLC_REDSHIFT,  PEAKMJD
 105  FORMAT('NOBS: ', I4, 5x, 'PIXSIZE: ', F6.3, 5x,  & 
             'REDSHIFT: 'F8.5, 5x, 'PEAKMJD: ', F9.3 )

    IF ( SUBSURVEY_NAME .NE. ' ' ) THEN
       WRITE(LUNOUT,10) 'SUBSURVEY:', SUBSURVEY_NAME
    ENDIF

! Nov 23 2019: add more info for simulated data
    IF ( WR_SIM_MAGOBS  ) THEN
!           write(LUNOUT,"('GALID: ',I9)" ) SNHOST_OBJID(1)
       LENF = INDEX(SNLC_FIELDLIST,' ') - 1
       write(LUNOUT,"('FIELD: ',A )" ) SNLC_FIELDLIST(1:LENF) ! includes overlaps
    ENDIF

! - - - - - - - - - -

! - - - - -
    CALL SET_SIMLIB_TEMPLATE_INFO()  ! Feb 5 2020

    IF ( FOUND_TEMPLATE_INFO ) THEN
        write(LUNOUT,"('TEMPLATE_ZPT:',$)")
        DO ifilt = 1, NFILTDEF_SURVEY
           write(LUNOUT,"(F8.3,$)" ) SIMLIB_TEMPLATE_ZPT(ifilt)
        ENDDO
        write(LUNOUT,"(' ')")

        write(LUNOUT,"('TEMPLATE_SKYSIG:',$)")
        DO ifilt = 1, NFILTDEF_SURVEY
           write(LUNOUT,"(F9.2,$)" ) SIMLIB_TEMPLATE_SKYSIG(ifilt)
        ENDDO
        write(LUNOUT,"(' ')")

    ENDIF

! May 29 2020: write TAKE_SPECTRUM keys if there are spectra
    CALL WRITE_SIMLIB_TAKE_SPECTRUM(LUNOUT)

! write column header
    WRITE(LUNOUT,110)
110   FORMAT(/,'#   MJD    IDEXPT  BAND  GAIN RDNOISE  SKYSIG  ',  & 
               'PSF1 PSF2 PSFRAT    ZP  ZPERR  MAG' )

    DO IFILT_OBS = 1, MXFILT_OBS
       LASTEP(IFILT_OBS) = -9
    ENDDO

! --------
    DO 300 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

      DO 301 IEP = EPMIN, EPMAX

         EP        = IEP
         MJD       = SNLC8_MJD(EP)
         IFILT_OBS = ISNLC_IFILT_OBS(ep)
         IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
         BAND      = FILTDEF_STRING(ifilt_obs:ifilt_obs)

         if ( EP < 0 ) goto 301  ! nothing valid to print

         GAIN      = SNGL( SNLC_GAIN(ep) )
         RDNOISE   = SNGL( SNLC_RDNOISE(ep) )
         SKYSIG    = SNGL( SNLC_SKYSIG(ep) )
         PSF1      = SNGL( SNLC_PSF_SIG1(ep) )
         PSF2      = SNGL( SNLC_PSF_SIG2(ep) )
         PSFRAT    = SNGL( SNLC_PSF_RATIO(ep) )
         ZP        = SNGL( SNLC_ZEROPT(ep) )
         ZPERR     = SNGL( SNLC_ZEROPT_ERR(ep) )

         cMAG      = '99' ; LENcMAG=2
         if ( WR_SIM_MAGOBS ) then
            MAG = SIM_EPMAGOBS(ep)
            if ( WRSET_SIM_MAGOBS .and. MAG > 98.0 ) goto 301
            write(cMAG,"(F8.4)") MAG
            LENcMAG = 8
         endif

         IF ( PSF1 < 1.0E-5 .and. SKYSIG > 1.0E-5 ) GOTO 301

! if simlib info is not available, estimate it from FLUX and FLUXERR.
! Mainly for low-z sample which has no meta data.

         FOUND_METADATA =  & 
             ( PSF1 > 1.0E-5 .and. SKYSIG > 1.0E-5 .and. ZP > 1.0)

         IF ( .not. FOUND_METADATA ) THEN
            GAIN    = 1.0 ;  RDNOISE = 1.0
            ZPERR   = DBLE( SIMLIB_ZPERR(ifilt_obs) )
            PSF2    = 0.0 ;  PSFRAT  = 0.0
            FLUX    = SNLC_FLUXCAL(EP)
            FLUXERR = SNLC_FLUXCAL_ERRTOT(EP)
            LAMOBS  = FILTOBS_LAMAVG(IFILT_OBS)

            if ( LAMOBS < 0.0 ) then
               c1err = 'Invalid/undefined band = ' // band
               c2err = 'Make sure all bands are defined in kcor file'
               CALL MADABORT(FNAM, C1ERR, C2ERR )
            endif

            CALL COMPUTE_SIMLIB_INFO(  & 
                     SNLC_CCID,FLUX, FLUXERR, LAMOBS,PIXSIZE,  & 
                     PSF1, SKYSIG, ZP)     ! returned
         ENDIF

         write(LUNOUT,120)  & 
              MJD, EP, BAND,     &  ! EP is a dummy ID
              GAIN, RDNOISE,  & 
              SKYSIG, PSF1, PSF2, PSFRAT,  & 
              ZP, ZPERR, cMAG(1:LENcMAG)

120        format('S: ',  & 
                  F9.3, 1x, I4, 2x, A2, 2x,   &  ! MJD, ID, BAND
                  F5.2, 2x, F5.2, 2x,         &  ! GAIN, RDNOISE
                  F7.2, 1x, 3(F6.2),2x,       &  ! SKYSIG, PSF info
                  F6.3, 1x, F5.3, 1x, A       &  ! ZP, ZPERR, MAG
                      )

301   CONTINUE
300   CONTINUE

    WRITE(LUNOUT,199) SNLC_CID
199   FORMAT('END_LIBID: ', I8, / )
    CALL FLUSH(6)

    RETURN
  END SUBROUTINE MAKE_SIMLIB_FILE


! =========================================
    SUBROUTINE SET_SIMLIB_TEMPLATE_INFO()

! Created Feb 2020
! If there is SKYSIG_T (template) info, set arrays
!    SIMLIB_TEMPLATE_ZPT[SKYSIG](ifilt)
! Search all epochs and extract typical epoch for each band.
! 
! Mar 26 2021: require SKYSIG_T > 0


! local variables

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SIMLIBCOM

    IMPLICIT NONE

    INTEGER ep, EPMIN, EPMAX, NEWMJD, IFILT, IFILT_OBS
    LOGICAL IS_SET(MXFILT_ALL), LZPT, LSKY
    REAL ZPT, SKYSIG

! ------------- BEGIN -------------

    DO IFILT = 1, NFILTDEF_SURVEY
      SIMLIB_TEMPLATE_ZPT(ifilt)    =  0.0
      SIMLIB_TEMPLATE_SKYSIG(ifilt) =  0.0
      IS_SET(IFILT) = .FALSE.
    END DO


    FOUND_TEMPLATE_INFO = .FALSE.

! find first epoch for each band

    DO 300 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 301 EP = EPMIN, EPMAX
         IFILT_OBS = ISNLC_IFILT_OBS(ep)
         IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
         ZPT       = SNLC_ZEROPT(ep)
         SKYSIG    = SNLC_SKYSIG_T(ep)
         LZPT      = (ZPT > 10.0)
         LSKY      = SKYSIG > 0.001
         if ( .not.IS_SET(ifilt) .and. LZPT .and. LSKY) then
           IS_SET(ifilt)                 = .TRUE.
           FOUND_TEMPLATE_INFO           = .TRUE.
           SIMLIB_TEMPLATE_ZPT(ifilt)    = ZPT
           SIMLIB_TEMPLATE_SKYSIG(ifilt) = SKYSIG
         endif
 301    CONTINUE
 300  CONTINUE

    RETURN
  END SUBROUTINE SET_SIMLIB_TEMPLATE_INFO

! ====================================
    SUBROUTINE SET_SIMLIB_SKYMAG()

! Check SURVEY_NAME for space or ground missing;
! set SKYMAG_LIST and PSF_FWHM accordingly.
! Beware that list of space surveys is hard-wired.


    USE SNPAR
    USE CTRLCOM
    USE SIMLIBCOM

    IMPLICIT NONE

    INTEGER i, LENS

! -------------- BEGIN --------------

    ISGROUND = .TRUE.  ! default
    IF ( INDEX(SURVEY_NAME,'HST'     ) > 0 ) ISGROUND = .FALSE.
    IF ( INDEX(SURVEY_NAME,'JWST'    ) > 0 ) ISGROUND = .FALSE.
    IF ( INDEX(SURVEY_NAME,'WFIRST'  ) > 0 ) ISGROUND = .FALSE.
    IF ( INDEX(SURVEY_NAME,'CANDELS' ) > 0 ) ISGROUND = .FALSE.
    IF ( INDEX(SURVEY_NAME,'ROMAN'   ) > 0 ) ISGROUND = .FALSE.  ! Oct 20 2025
    IF ( INDEX(SURVEY_NAME,'EUCLID'  ) > 0 ) ISGROUND = .FALSE.  ! Oct 20 2025
    IF ( INDEX(SURVEY_NAME,'CSST'    ) > 0 ) ISGROUND = .FALSE.  ! Oct 20 2025

    LENS = INDEX(SURVEY_NAME,' ') - 1

    IF ( ISGROUND ) THEN
        write(6,20) SURVEY_NAME(1:LENS), 'GROUND'
        NLIST_SKY = NLIST_SKY_GROUND
        DO i = 1, NLIST_SKY_GROUND
          SKYLAM_LIST(i) = SKYLAM_GROUND_LIST(i)
          SKYMAG_LIST(i) = SKYMAG_GROUND_LIST(i)
        ENDDO
        PSF_FWHM_GUESS   = 1.0   ! arcSec
        PIXSIZE_GUESS    = 0.5   ! arcSec
        ADD_SKYSIG_PIX   = 0.05  ! skySigma per pixel to add in loop

    ELSE
        write(6,20) SURVEY_NAME(1:LENS), 'SPACE'
        NLIST_SKY = NLIST_SKY_SPACE
        DO i = 1, NLIST_SKY_SPACE
          SKYLAM_LIST(i) = SKYLAM_SPACE_LIST(i)
          SKYMAG_LIST(i) = SKYMAG_SPACE_LIST(i)
        ENDDO
        PSF_FWHM_GUESS   = 0.2  ! arcSec
        PIXSIZE_GUESS    = 0.1
        ADD_SKYSIG_PIX   = 0.01
    ENDIF

20    format(T5,'Survey ',A, ' --> ', A, ' instrument params.' )
    write(6,40) 'PSF_FWHM ' , PSF_FWHM_GUESS, 'arcSeconds'
    write(6,40) 'PIXSIZE  ' , PIXSIZE_GUESS, 'arcSeconds'
40    format(T10, 'Set ', A,' = ', F6.3, 2x, A )

    DO i = 1, NLIST_SKY
       write(6,50) SKYLAM_LIST(i), SKYMAG_LIST(i)
50       format(T10, 'Set SKYMAG(',F7.0,' A) = ', F5.2, ' / asec^2'  )
    ENDDO

! - - - - - -
! Write info to SIMLIB file
    if ( ISGROUND ) then
       write(LUNOUT,11) 'GROUND'
    else
       write(LUNOUT,11) 'SPACE'
    endif
11      format(/, '# Assume instrument parameters for ', A)
      DO i = 1, NLIST_SKY
         write(LUNOUT,12) SKYLAM_LIST(i), SKYMAG_LIST(i)
12         format('# Assume SKYMAG(',F7.0,') = ', F5.2,' mag/asec^2')
      ENDDO

    RETURN
  END SUBROUTINE SET_SIMLIB_SKYMAG

! ====================================
    SUBROUTINE SET_SIMLIB_ZPERR()

! Set SIMLIB_ZPERR(ifilt_obs) using &SNLCINP input
!   SIMLIB_ZPERR_LIST = 'abc .01 def .02 hij .014'
! 
! If SIMLIB_ZPERR_LIST='', set all ZPERR = 0.01
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE SIMLIBCOM
    USE FILTCOM

    IMPLICIT NONE

    INTEGER ifilt, IFILTDEF, NF, iafilt(MXFILT_OBS)
    REAL xafilt(MXFILT_OBS), ZPERR

! -------------- BEGIN --------------

    DO IFILT = 1, MXFILT_ALL
      SIMLIB_ZPERR(IFILT) = 0.01  ! default
    ENDDO

    IF ( SIMLIB_ZPERR_LIST .EQ. '' ) RETURN

    CALL PARSE_FILTSTRING(0, SIMLIB_ZPERR_LIST, NF,iafilt,xafilt)
    DO ifilt    = 1, NF
       ifiltdef = IAFILT(ifilt)
       ZPERR    = xafilt(ifilt)
       SIMLIB_ZPERR(ifiltdef) = ZPERR
    ENDDO

    RETURN
  END SUBROUTINE SET_SIMLIB_ZPERR


! =======================================================
    SUBROUTINE COMPUTE_SIMLIB_INFO(CCID,FLUX,FLUXERR,LAMOBS,PIXSIZE,  & 
                                     PSF,SKYSIG,ZP)

! 
! Created Jan 26 2017 by R.Kessler
! When &SNLCINP namelist SIMLIB_OUTFILE is specified (fileName),
! and the data files have no meta data, this function is called
! as a hack to compute SIMLIB row entries for PSF, SKYSIG and ZP.
! 
! Strategy here is numerical solution with assumptions on
! the PSF and sky brightness. Maybe can be done analytically,
! but I didn't try.

! * Fix PSF.
! * assume a sky brightness vs. wavelength estimate from LSST
!     deep-drill cadence. Interpolate SKYMAG_REF based on
!     input LAMOBS = central wavelength of filter.
! * loop over grid of SKYSIG values, increments of 0.1
! * for each SKYSIG, SNR_meas = Fpe/sqrt(Fpe+SQSKYSG) is solved
!    (quadratic equation) for ZP(pe) between FLUXCAL and MAG.
! * Use ZP and SKYSIG to compute "SKYMAG_predict"
! * SKYSIG with smallest |SKYMAG_predict - SKYMAG_REF|
!    is the result.
! 
! Something similar was previously done for the low-z similb
! using only CFA3_KEPLERCAM; the resulting CFA3_KEPLERCAM.SIMLIB
! was used to simulate low-z in
!    http://adsabs.harvard.edu/abs/2013ApJ...764...48K
! and also in the JLA(Betoule 2014) and PS1(Sconlic 2014)
! analyses.
! Unfortunately, I cannot find the script or used to create
! CFA3_KEPLERCAM.SIMLIB .
! 
! May 2024: if LAMOBS is too large, abort with clear message.
! 
! ----------------------------------------


    USE SNPAR
    USE SNLCCOM
    USE SIMLIBCOM
    USE CTRLCOM

    IMPLICIT NONE

! subroutine args
    CHARACTER CCID*(MXCHAR_CCID)  ! (I) for error message only

    DOUBLE PRECISION  & 
         FLUX               &  ! (I) FLUXCAL
        ,FLUXERR            &  ! (I) error on above
        ,LAMOBS             &  ! (I) central wavelength of filter
        ,PIXSIZE            &  ! (I) pixel size, arcsec
        ,PSF, SKYSIG, ZP   ! (O)


! local args
    INTEGER NLOOP, MXLOOP
    REAL*8  & 
         PSF_SIG, SKYMAG, TMP_SKYMAG  & 
        ,TMP_SKYSIG_PIX, TMP_SKYSIG_ASEC, TMP_SQSKYSIG_TOT  & 
        ,TMP_ZP, TMP_Fsky_pix,  TMP_FSKY_ASEC, TMP_Fpe  & 
        ,SNR, SQSNR, NEA, MAG, MAG_MAX,  ARG, AREA_PIXEL  & 
        ,SKYMAGDIF, SKYMAGDIF_MIN, SKYMAGDIF_LAST, SNR_CHECK

    CHARACTER FNAM*20, C1ERR*80, C2ERR*80

    REAL*8   INTERP_1DFUN
    EXTERNAL INTERP_1DFUN

    LOGICAL LDMP

! --------------- BEGIN ---------------

! init output
    PSF = 0.0; SKYSIG=0.0 ; ZP=0.0

!      LDMP = .TRUE.
    LDMP = .FALSE.

    FNAM = 'COMPUTE_SIMLIB_INFO'

! set MAG_MAX based on Flux = FluxErr
    MAG_MAX = ZP_FLUXCAL - 2.5*LOG10(FLUXERR)
    if ( FLUX > FLUXERR ) then
       MAG = ZP_FLUXCAL - 2.5*LOG10(FLUX)
    else
       MAG = MAG_MAX
    endif

! abort on insane FLUXERR
    IF ( FLUXERR < 1.0E-8 ) THEN
       write(c1err,161) 'FLUXERR', FLUXERR, FLUX
161      format('Insane ', A,' value = ', F12.4,  & 
               3x, 'for FLUX=',F12.4)
       c2err = 'Check CID = ' // CCID
       CALL MADABORT(FNAM, c1err, c2err)
    ENDIF

! hard-wire PSF since we can't solve for all 3 items
    PSF      = PSF_FWHM_GUESS         ! set output
    PSF_SIG  = PSF_FWHM_GUESS/2.355   ! PSF in sigma, arcsec
    NEA      = 4.0*PI * (PSF_SIG*PSF_SIG)  ! Noise-equiv area, arcSec^2
    AREA_PIXEL = PIXSIZE**2     ! area of 1 pixel, arcsec^2

    SQSNR = (FLUX**2) / ( FLUXERR**2 + (0.0*FLUX)**2 )
    IF ( SQSNR < 2.0 ) SQSNR = 2.0  ! protect crazy SIMLIB values
    SNR   = sqrt(SQSNR)
    SKYMAGDIF_MIN = 9999.

    if ( LAMOBS > SKYLAM_LIST(NLIST_SKY) ) THEN
       write(C1ERR,601) LAMOBS, SKYLAM_LIST(NLIST_SKY)
601      format('LAMOBS=',F8.0,' exceeds max LAM=', F8.0 )
       C2ERR = 'Need to extend SKYLAM and SKYMAG ' //  & 
                 'arrays in SIMLIBCOM'
       CALL MADABORT(FNAM, c1err, c2err)
    endif

! linear interpolation to get estimate of Sky mag per asec^2
    SKYMAG = interp_1dfun(1, LAMOBS, NLIST_SKY,  & 
                        SKYLAM_LIST, SKYMAG_LIST,  & 
                        FNAM//char(0), 20)

    IF ( LDMP ) THEN
      print*,' xxx --------------------------------- '
      print*,' xxx FLUX = ', sngl(FLUX), ' +- ', sngl(FLUXERR)
      print*,' xxx LAMOBS, SKYMAG = ', sngl(LAMOBS), sngl(SKYMAG)
      print*,' xxx NEA  = ', sngl(NEA)
      print*,' xxx PSF  = ', sngl(PSF_FWHM_GUESS), ' FWHM '
      print*,' xxx SNR  = ', sngl(SNR)
      print*,' xxx AREA_PIXEL = ', sngl(AREA_PIXEL)
    ENDIF

    SKYMAGDIF_LAST = 9999; SKYMAGDIF = 999. ;
    TMP_SKYSIG_PIX = ADD_SKYSIG_PIX
    NLOOP = 0 ; MXLOOP = 100000

    DO WHILE ( SKYMAGDIF < SKYMAGDIF_LAST )
       NLOOP = NLOOP + 1
       SKYMAGDIF_LAST   = SKYMAGDIF

!  increment skysig (pe) per pix
       TMP_SKYSIG_PIX   = TMP_SKYSIG_PIX + ADD_SKYSIG_PIX
       TMP_SKYSIG_ASEC  = TMP_SKYSIG_PIX / PIXSIZE  ! sqrt(AREA_PIXEL)

       TMP_SQSKYSIG_TOT = NEA * (TMP_SKYSIG_ASEC**2) ! total SQSKYSIG
       TMP_FSKY_PIX     = TMP_SKYSIG_PIX**2         ! sky val/pix, pe
       TMP_FSKY_ASEC    = TMP_FSKY_PIX / AREA_PIXEL ! sky val/asec^2

       ARG           = 1.0 + 4.0 * TMP_SQSKYSIG_TOT / SQSNR
       TMP_Fpe       = (SQSNR/2.0) * ( 1.0 + sqrt(ARG) )
       TMP_ZP        = MAG + 2.5*LOG10(TMP_Fpe)

! check how close computed TMP_SKYMAG is to expected SKYMAG
       TMP_SKYMAG    = TMP_ZP - 2.5*log10(TMP_FSKY_ASEC)
       SKYMAGDIF     = abs(TMP_SKYMAG - SKYMAG)

       if ( SKYMAGDIF < SKYMAGDIF_MIN ) then
          SKYMAGDIF_MIN = SKYMAGDIF
          SKYSIG    = TMP_SKYSIG_PIX/PIXSIZE     ! per arcSec^2
          ZP        = TMP_ZP                     ! pe per exposure
       endif

       IF ( NLOOP > MXLOOP ) THEN
          CALL PRINT_PREABORT_BANNER(FNAM(1:18)//char(0),40)
          print*,'   CID = ', CCID
          print*,'   FLUX = ', sngl(FLUX), ' +- ', sngl(FLUXERR)
          print*,'   LAMOBS, SKYMAG = ', sngl(LAMOBS), sngl(SKYMAG)
          print*,'   NEA  = ', sngl(NEA)
          print*,'   PSF  = ', sngl(PSF_FWHM_GUESS), ' FWHM '
          print*,'   SNR  = ', sngl(SNR)
          print*,'   AREA_PIXEL = ', sngl(AREA_PIXEL)
          print*,'   NLOOP = ', NLOOP
          print*,'   last SKYSIG = ', sngl(TMP_SKYSIG_PIX/PIXSIZE)
          print*,'   last ZP     = ', sngl(TMP_ZP)
          C1ERR = 'Cannot determine SKYSIG '
          C2ERR = 'after NLOOP tries '
          CALL MADABORT(FNAM, c1err, c2err)
       ENDIF

       IF ( LDMP ) THEN
          SNR_CHECK = TMP_Fpe / sqrt(TMP_Fpe + TMP_SQSKYSIG_TOT )
         write(6,60) NLOOP, TMP_SKYSIG_PIX, SKYMAGDIF, SNR_CHECK/SNR
60         format(' xxx ', I3,' SkySig=',F5.2, 5x, 'SKYMAGDIF=',F8.5,  & 
               5x,'SNR(test)=',F6.4 )
       ENDIF
    ENDDO

    IF ( LDMP ) THEN
      print*,' xxx FINAL PSF,SKYSIG,ZP = ',  & 
             sngl(PSF), sngl(SKYSIG), sngl(ZP)
        print*,' xxx DEBUG STOP' ; CALL FLUSH(6)
        CALL EXIT(EXIT_ERRCODE)
    ENDIF

    RETURN
  END SUBROUTINE COMPUTE_SIMLIB_INFO

! ==================
    SUBROUTINE WRITE_SIMLIB_TAKE_SPECTRUM(LUN)

! Created May 29 2020
! If there are spectra, write TAKE_SPECTRUM keys to
! SIMLIB header (LUNOUT)
! 


    USE SNDATCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

    INTEGER LUN ! (I) write to this LUN

! discard if FLUXERR/FLUXERR_MEDIAN exceeds this value
! --> avoid crazy large errors that result in bad FLUXUM/FLUXERRSUM.
    REAL, PARAMETER :: TOL_FLUXERR_MEDIAN_RATIO = 50.0

    REAL*8  TOBS, LAMMIN, LAMMAX, SNR
    REAL*8  MJD, LAMMIN_BIN, LAMMAX_BIN
    REAL*8  FLUXSUM, SQERRSUM, FLAM, FLAMERR, DLAM
    REAL*8  FLAMERR_AVG, FLAMERR_RMS, FLAMERR_MEDIAN, RATIO
    INTEGER ispec, NLAMBIN, ilam, NBIN_REJECT
    INTEGER LEN_TOBS, LEN_SNR, LEN_LAMMIN, LEN_LAMMAX
    CHARACTER  & 
          STR_TOBS*20, STR_SNR*20,  & 
          STR_LAMMIN*20, STR_LAMMAX*20

    EXTERNAL ARRAYSTAT
! ----------- BEGIN --------

    DO 222 ispec = 1, NSPECTRUM

       CALL RDSPEC_DRIVER(ispec)  ! apr 5 2021

       MJD       = MJD_SPECTRUM(ispec)
       TOBS      = TOBS_SPECTRUM(ispec)
       NLAMBIN   = NLAMBIN_SPECTRUM(ispec)
       LAMMIN    = LAMMIN_SPECTRUM(1)
       LAMMAX    = LAMMAX_SPECTRUM(NLAMBIN)

       FLUXSUM = 0.0 ; SQERRSUM = 0.0

       CALL ARRAYSTAT(NLAMBIN, FLAMERR_SPECTRUM,  & 
                        FLAMERR_AVG, FLAMERR_RMS, FLAMERR_MEDIAN)

       NBIN_REJECT = 0
       DO ilam = 1, NLAMBIN
           FLAM    = FLAM_SPECTRUM(ilam)
           FLAMERR = FLAMERR_SPECTRUM(ilam)
           LAMMIN_BIN = LAMMIN_SPECTRUM(ilam)
           LAMMAX_BIN = LAMMAX_SPECTRUM(ilam)

           RATIO = FLAMERR/FLAMERR_MEDIAN
           if ( RATIO > TOL_FLUXERR_MEDIAN_RATIO ) then
              NBIN_REJECT = NBIN_REJECT + 1
           else
             DLAM     = 2.0 ! LAMMAX_BIN - LAMMIN_BIN
             FLUXSUM  = FLUXSUM  + (FLAM*DLAM)
             SQERRSUM = SQERRSUM + (FLAMERR*DLAM)**2
           endif
       ENDDO
       SNR = FLUXSUM / sqrt(SQERRSUM)
       if ( NBIN_REJECT > 0 ) then
          write(6,212) SNLC_CCID(1:ISNLC_LENCCID),  & 
                       NBIN_REJECT, NLAMBIN, MJD
 212        format(t5, 'WARNING(',A,'): flag ', I4,' of ', I5, 2x,  & 
              'bins as FLAMERR outlier (MJD=',F9.3,')' )
       endif

       CALL dble2string_nospace(TOBS, STR_TOBS, LEN_TOBS)
       CALL dble2string_nospace(SNR,  STR_SNR,  LEN_SNR)
       CALL dble2string_nospace(LAMMIN,  STR_LAMMIN, LEN_LAMMIN)
       CALL dble2string_nospace(LAMMAX,  STR_LAMMAX, LEN_LAMMAX)

       WRITE(LUN,225)  & 
            STR_TOBS(1:LEN_TOBS), STR_SNR(1:LEN_SNR),  & 
            STR_LAMMIN(1:LEN_LAMMIN), STR_LAMMAX(1:LEN_LAMMAX)
 225     FORMAT('TAKE_SPECTRUM:  TOBS(', A, ')', 3x,  & 
                 'SNR(',A,')', 3x, 'SNR_LAMOBS(',A,':',A, ')' )

 222  CONTINUE

    RETURN
  END SUBROUTINE WRITE_SIMLIB_TAKE_SPECTRUM

! ================================================
    SUBROUTINE dble2string_nospace(VAL, STRING, LEN)
    REAL*8 VAL            ! (I) double value to convert
    CHARACTER STRING*(*)  ! (O) string with float value
    INTEGER  LEN          ! (O) length of string

    write(STRING,240) VAL
 240  format(F10.1)
    STRING = ADJUSTL(STRING)
    LEN = INDEX(STRING,' ') - 1
    RETURN
  END SUBROUTINE dble2string_nospace

    SUBROUTINE float2string_nospace(VAL, STRING, LEN)
    REAL VAL              ! (I) float value to convert
    CHARACTER STRING*(*)  ! (O) string with float value
    INTEGER  LEN          ! (O) length of string

    write(STRING,240) VAL
 240  format(F10.1)

    STRING = ADJUSTL(STRING)
    LEN = INDEX(STRING,' ') - 1

    RETURN
  END SUBROUTINE float2string_nospace
! =======================================================
#if defined(SNANA) || defined(SNFIT)
    SUBROUTINE FLUXERRCALC(iepoch, FLUXCAL, FLUXCAL_ERR )
! 
! Compute error from SKY,PSF and ZPT the same way as the
! simulation does.
! 
! Fill
!   SNLC_SKYFLUXCAL(iepoch)           ! fluxcal per pixel
!   SNLC_FLUXCAL_ERRCALC(iepoch)      ! calculated flux error
!   SNLC_FLUXCAL_HOSTERRCALC(iepoch)  ! idem for host noise contribution
! where
!   SKYFLUX_ERR_PHOTO^2 = TOTAL_ERR^2 - SIGNAL_ERR^2
! 
! and SKYFLUX_ERR_CALC is determined from the
! effective aperture.
! 
! Feb 12,2012 - complete re-write.
!   The ZEROPT is now for the search run instead of for the
!   template run, so no need for obsolete KSUM to translate ZPT.
!   The 'template' noise is ingored since we assume that
!   many template images exist. To include template noise,
!   the SKY_SIG entry in the data file can be increased.
!   Another caveat is that FLUXERR is computed here in ADU,
!   but the reference FLUX can be in ADU or uJy.  To avoid
!   this confusion, the FLUXCAL error is computed from
!     FLUXCAL_ERRCALC = FLUXADU_ERRCALC * ADUSCALE / 10**[0.4*ZPT]
!   Thus comparing the calculated FLUXCAL_ERRCALC to the real
!   FLUXCAL_ERRTOT is always unamgiguous
! 
!   The calibrated FLUXCAL and its error are passed as arguments
!   to allow using fitted quantities that may work better for
!   low SNR.
! 
! 
! Sep 21 2017:
!  + implement CHECK_SNANA_DUMP
!  + check PIXSIZE !!!
!  + allow calc for any format (previously skipped TERSE)
! 
! Mar 2 2018:
!  + REAL -> REAL*8
! 
! Nov 15 2019: load SNLC_FLUXCAL_ERRTEST = ERRCALC/ERRTRUE
! 
! --------------------------------------------

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER iepoch                ! (I) epoch & SN indices
    REAL    FLUXCAL, FLUXCAL_ERR  ! (I) FLUXCAL and its error

! local var

    REAL  & 
         GAIN  & 
        ,RDNOISE_pe, RDNOISE_ADU  & 
        ,SQERR_TOT  & 
        ,SQERR_SIGNAL  & 
        ,SQERR_SKY       &  ! integrated over effective aperture
        ,SQERR_ZP  & 
        ,SQSUM, TMPRAT, TMPERR  & 
        ,SKY_SIG_ADU, SKY_SIG2_ADU  &  ! search and template (ADU/pixel)
        ,SKY_SIGTOT_ADU  & 
        ,SKY_AVG_ADU  & 
        ,SBFLUXCAL, SBFLUXCAL_ERR, SBFLUXADU  & 
        ,FLUXADU, FLUXADU_ERR  & 
        ,ZP,  ZPERR, ZPSCALE, ZPDIF  & 
        ,SQ1, SQRD,  SQERR_HOST, HOST_SKY_RATIO, ERRCALC, ERRTRUE  & 
        ,SKYFLUX_ERR_PHOTO  & 
        ,SKYFLUX_ERR_CALC  & 
        ,PIXSIZE

    INTEGER IFILT_OBS, IFILT, DUMPFLAG, LENS
    CHARACTER cfilt*2, FNAM*12, STRING_INVALID*200

! function
    EXTERNAL NoiseEquivAperture, CHECK_SNANA_DUMP, ABORT_SNANA_DUMP
    REAL*8   NoiseEquivAperture
    INTEGER  CHECK_SNANA_DUMP_forC

    REAL*8  AREA, PSFSIG1, PSFSIG2, PSFRATIO, MJD8  ! args for above

! ------------- BEGIN --------------

    ERRCALC = -9.0
    SNLC_FLUXCAL_ERRCALC(iepoch)     = -9.0
    FNAM = 'FLUXERRCALC'

! initialize output

    SKYFLUX_ERR_PHOTO                = -9.0
    SKYFLUX_ERR_CALC                 = -9.0

! strip off observing conditions into local variables.

    ifilt_obs  = ISNLC_IFILT_OBS(iepoch)
    IFILT      = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

    SKY_SIG_ADU   = SNLC_SKYSIG(iepoch)   ! search sky-sigma, ADU/pix
    SKY_SIG2_ADU  = max(0.0,SNLC_SKYSIG_T(iepoch))   ! idem for template
    PSFSIG1       = DBLE(SNLC_PSF_SIG1(iepoch))    ! Gauss sigma, pixel
    PSFSIG2       = DBLE(SNLC_PSF_SIG2(iepoch))    ! idem, 2nd component
    PSFRATIO      = DBLE(SNLC_PSF_RATIO(iepoch))
    GAIN          = SNLC_GAIN(iepoch)             ! e-/ADU
    RDNOISE_pe    = SNLC_RDNOISE(iepoch)          ! e-/pix
    ZP            = SNLC_ZEROPT(iepoch)
    ZPERR         = SNLC_ZEROPT_ERR(iepoch)
    SBFLUXCAL     = SNHOST_SBFLUXCAL(ifilt)
    SBFLUXCAL_ERR = SNHOST_SBFLUXCAL_ERR(ifilt)
    PIXSIZE       = SNLC_PIXSIZE
    ERRTRUE       = SNLC_FLUXCAL_ERRTOT(iepoch)

    SKY_SIGTOT_ADU = SQRT(SKY_SIG_ADU**2 + SKY_SIG2_ADU**2)

! check a few things.
    IF ( FLUXERRCALC_ZPTERR  .GE. 0.0 ) THEN
        ZPERR = FLUXERRCALC_ZPTERR
    ENDIF

    IF ( GAIN > 0.0 ) THEN
      RDNOISE_ADU  = RDNOISE_pe/GAIN     ! translate RDNOISE from e- to ADU
      SKY_AVG_ADU  = GAIN * (SKY_SIGTOT_ADU**2) ! needed to scale host noise
    ELSE
      SKY_AVG_ADU  = 0.0
      RDNOISE_ADU  = 0.0
    ENDIF

! ----------------------------------------------
! check for invalid/undefined values

    STRING_INVALID = ''
    if ( SKY_SIG_ADU   < 1.0E-4  ) then
       CALL APPEND_STRING(STRING_INVALID, 'SKY_SIG_ADU')
    endif
    if ( RDNOISE_ADU   < 0.0  ) then
       CALL APPEND_STRING(STRING_INVALID, 'RDNOISE')
    endif
    if ( PSFSIG1       < 0.001   ) then
       CALL APPEND_STRING(STRING_INVALID, 'PSFSIG1')
    endif
    if ( PSFSIG2       < 0.0     ) then
       CALL APPEND_STRING(STRING_INVALID, 'PSFSIG2')
    endif
    if ( PSFRATIO      < 0.0     ) then
       CALL APPEND_STRING(STRING_INVALID, 'PSFRATIO')
    endif
    if ( GAIN          < 0.01    ) then
       CALL APPEND_STRING(STRING_INVALID, 'GAIN')
    endif
    IF ( ZP < 0.01 .or. ZP > 60.0    ) then
       CALL APPEND_STRING(STRING_INVALID, 'ZP')
    endif
    IF ( ZPERR < 0.0 .or. ZPERR > 3.99     ) then
       CALL APPEND_STRING(STRING_INVALID, 'ZPERR')
    endif
    IF ( PIXSIZE       < 0.0001  ) then
       CALL APPEND_STRING(STRING_INVALID, 'PIXSIZE')
    endif

    IF ( STRING_INVALID .NE. '' ) GOTO 777  ! skip to dump

! set logical as soon as we get here.
    DO_FLUXERRCALC = .TRUE.

! use ZPT to convert calibated flux back into ADU
! (since FLUX in the data file can be either ADU or uJy)

    ZPDIF       = ZP - SNGL(ZP_FLUXCAL)
    ZPSCALE     = 10**(0.4*ZPDIF)   ! FLUXADU/FLUXCAL ratio
    FLUXADU     = ZPSCALE * FLUXCAL
    FLUXADU_ERR = ZPSCALE * FLUXCAL_ERR

! get signal error, but avoid problems with large negative fluxes
    tmprat = FLUXADU / FLUXADU_ERR
    if ( tmprat .GT. 0.5 ) then
       SQERR_SIGNAL = FLUXADU / GAIN   ! ADU^2
    else
       SQERR_SIGNAL = 0.0
    endif

! get effective area for sky-bkg = 1/[integral PSF^2], pixelsize^2
    AREA    = SNLC_PSF_NEA(iepoch) ! Mar 1 2021

! multiply skyerr^2 (per pixel) by effective area.
! Include optional read noise

    SQ1  = SKY_SIGTOT_ADU  * SKY_SIGTOT_ADU
    SQRD = RDNOISE_ADU  * RDNOISE_ADU
    SQERR_SKY = SNGL( AREA * ( SQ1 + SQRD ) )

! Include host noise if surface brightness (SB) is defined.

    SQERR_HOST = 0.0
    IF ( SBFLUXCAL > 0.001 .and. SKY_AVG_ADU > 0.0 ) THEN
       SBFLUXCAL  = SBFLUXCAL*(PIXSIZE**2)  ! -> FLUXCAL per pixel
       SBFLUXADU  = SBFLUXCAL * ZPSCALE
       SQERR_HOST = SNGL(AREA*SQ1*(SBFLUXADU/SKY_AVG_ADU)) ! scale sky noise^2
       SQERR_SKY  = SQERR_SKY + SQERR_HOST
    ENDIF
    HOST_SKY_RATIO = sqrt ( SQERR_HOST / ( SQERR_SKY-SQERR_HOST) )

    SKYFLUX_ERR_CALC = SQRT(SQERR_SKY)    ! ADU

! compute calibrated skyflux per pixel^2
!  (not for error calc, but to compare with host surface brightess)

    SNLC_SKYFLUXCAL(iepoch) =  & 
         ((SKY_SIGTOT_ADU**2)*GAIN) / ZPSCALE

! compute error from zero-point smearing

    TMPERR    = 10**(0.4*ZPERR) - 1.0
    SQERR_ZP  = (FLUXADU * TMPERR)**2

! now combine calculated sky error and signal poisson noise
! for best estimate of flux error. This is what the simulation does.

    SQSUM = SQERR_SKY + SQERR_SIGNAL + SQERR_ZP
    ERRCALC                          = SQRT(SQSUM)/ZPSCALE
    SNLC_FLUXCAL_ERRCALC(iepoch)     = ERRCALC
    SNLC_FLUXCAL_HOSTERRCALC(iepoch) = SQRT(SQERR_HOST)/ZPSCALE

    IF ( ERRCALC > 0.0 .and. ERRTRUE > 0.0 ) THEN
      SNLC_FLUXCAL_ERRTEST(iepoch)  = ERRCALC/ERRTRUE
    ENDIF

! ------------------------------------
! Get photometry sky-error by subracting
! photo-stat signal error from total error

    SQERR_TOT    = FLUXADU_ERR * FLUXADU_ERR
    SQERR_SKY    = SQERR_TOT - SQERR_SIGNAL

    IF ( SQERR_SKY .GT. 0.0  ) THEN
       SKYFLUX_ERR_PHOTO = SQRT(SQERR_SKY)  ! ADU
    ELSE
       SKYFLUX_ERR_PHOTO = -9.0
    ENDIF

! --------------------------
! check for DUMP

777   CONTINUE
    IF ( DUMP_STRING .EQ. '' ) RETURN
    cfilt    = filtdef_string(ifilt_obs:ifilt_obs)
    MJD8     = SNLC8_MJD(iepoch)
    DUMPFLAG = CHECK_SNANA_DUMP_forC(FNAM,SNLC_CCID,CFILT,MJD8)

    if ( DUMPFLAG > 0 ) then

      print*,' '
      print*,' XXX ---------------------------------------------- '
      print*,' XXX FLUXERRCALC DUMP: '
      print*,' XXX CID=', SNLC_CCID(1:ISNLC_LENCCID),  & 
               '  MJD=', sngl(MJD8), '  FILTER=',CFILT
      print*,' XXX FLUXADU = ', FLUXADU,' +- ', FLUXADU_ERR
      print*,' XXX FLUXCAL = ', FLUXCAL,' +- ', FLUXCAL_ERR
      print*,' XXX ZPT(AVG,ERR,SCALE) = ',  & 
                     ZP, ZPERR, ZPSCALE
      print*,' XXX PSFSIG1, PSFSIG2   = ',  & 
                     sngl(PSFSIG1), sngl(PSFSIG2),' Gauss pix'
      print*,' XXX GAIN               = ', GAIN
      print*,' XXX PIXSIZE, AREA      = ', PIXSIZE, sngl(AREA)
      print*,' XXX RDNOISE_ADU        = ', RDNOISE_ADU
      print*,' XXX SKY_SIG_ADU        = ', SKY_SIG_ADU,  & 
                   ' ADU/pix (search)'
      print*,' XXX SKY_SIG2_ADU       = ', SKY_SIG2_ADU,  & 
                   ' ADU/pix (template)'
      print*,' XXX SIG(HOST)/SIG(SKY) = ', HOST_SKY_RATIO
      print*,' XXX SKY_AVG_ADU        = ', SKY_AVG_ADU
      print*,' XXX SQERR_SIGNAL       = ', SQERR_SIGNAL,' ADU'
      print*,' XXX SQERR[SKY,TOT]     = ',  & 
              SQERR_SKY, SQERR_TOT,' ADU '
      print*,' XXX SIM_MAGOBS         = ', SIM_EPMAGOBS(iepoch)
      print*,' XXX SIM_LIBID          = ', SIM_LIBID

      write(6,60) SKYFLUX_ERR_PHOTO, SKYFLUX_ERR_CALC
60      format(T2,' XXX SKYERR_ADU(PHOTO,CALC) =',F10.3,F10.3 )


      write(6,61) ERRTRUE, ERRCALC
61      format(T3,'XXX FLUXCAL_ERR(PHOTO,CALC)=',F10.3,F10.3 )
      print*,'   '

      IF ( STRING_INVALID .NE. '' ) THEN
        LENS = INDEX(STRING_INVALID,' ' ) - 1
        print*,' XXX Invalid quantities: ', STRING_INVALID(1:LENS)
      ENDIF

      if ( DUMPFLAG == 2 ) CALL ABORT_SNANA_DUMP()

    endif

    RETURN
  END SUBROUTINE FLUXERRCALC
#endif

! =====================================
    SUBROUTINE APPEND_STRING(STRING,APPEND)
! return STRING = STRING // ':' // APPEND
! Note that there are no blank spaces. Each item is separated by colon.

    CHARACTER STRING*(*), APPEND*(*)
    INTEGER LENS
! ----------- BEGIN ------------
    LENS = INDEX(STRING,' ' ) - 1
    IF ( LENS == 0 ) then
      STRING = APPEND
    ELSE
      STRING = STRING(1:LENS) // ':' // APPEND
    ENDIF

    RETURN
  END SUBROUTINE APPEND_STRING

! =====================================
    INTEGER FUNCTION CHECK_SNANA_DUMP_forC(FUN,CCID,BAND,MJD)

! add char(0) terminators to strings and then call
! C function CHECK_SNANA_DUMP(...)


! inputs
    CHARACTER FUN*(*), CCID*(*), BAND*(*)
    REAL*8  MJD

! local args
    INTEGER LEN_FUN, LEN_CCID, LEN_BAND, DUMPFLAG

! Declare C function
    EXTERNAL CHECK_SNANA_DUMP
    INTEGER CHECK_SNANA_DUMP
! ------------- BEGIN --------------

    LEN_FUN  = INDEX(FUN,  ' ' ) - 1
    LEN_CCID = INDEX(CCID, ' ' ) - 1
    LEN_BAND = INDEX(BAND, ' ' ) - 1

!      print*,' xxx FUN = ', FUN, LEN_FUN
!      print*,' xxx CCID = ', CCID, LEN_CCID
!      print*,' xxx BAND = ', BAND, LEN_BAND

    DUMPFLAG = CHECK_SNANA_DUMP(  & 
              FUN(1:LEN_FUN) // char(0) ,  & 
              CCID(1:LEN_CCID) // char(0) ,  & 
              BAND(1:LEN_BAND) // char(0),  & 
              MJD  )

    CHECK_SNANA_DUMP_forC = DUMPFLAG

    RETURN
  END FUNCTION CHECK_SNANA_DUMP_forC

! ==========================================
    SUBROUTINE CCD_AREAFRAC(ep)

! Created aug 2014 by R.Kessler
! Use input globals NXPIX,NYPIX, XPIX(ep), YPIX(ep)
! to compute SNLC_AREAFRAC(ep) = A/AMAX where
! 
! A    = area of largest rectangle containing XPIX, YPIX
! AMAX = area of CCD.
! 
! For uniform CCD illumination, the distribution of AREAFRAC
! is uniform; AREAFRAC=0 at center of CCD, and AREAFRAC=1.00
! at edge.
! 
! Also note that this works only for data because there is
! no pixel-coord input to the simulation ... but maybe
! someday XPIX and YPIX keys will be added to SIMLIB header.
! 
! --------------

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER ep  ! (I) epoch index

    REAL XCEN, YCEN, YSCALE, X, Y, xy, A, AMAX

! ------------- BEGIN -------

    SNLC_AREAFRAC(ep) = -9.0 ;
    if ( SNLC_NXPIX .LE. 0.0 .or. SNLC_NYPIX .LE. 0.0 ) RETURN

!      print*,' xxx NXPIX,NYPIX = ', SNLC_NXPIX, SNLC_NYPIX

    XCEN   = SNLC_NXPIX/2.0
    YCEN   = SNLC_NYPIX/2.0
    YSCALE = SNLC_NYPIX/SNLC_NXPIX

    X  = ABS(SNLC_XPIX(ep) - XCEN)
    Y  = ABS(SNLC_YPIX(ep) - YCEN) / YSCALE  ! --> scale to square

    if ( X > Y ) THEN
       xy = X  ! half-len of square
    else
       xy = Y
    endif

    AMAX = SNLC_NXPIX**2  ! max possible area
    A    = 4. * xy * xy   ! area containing XPIX, YPIX/Yscale

    SNLC_AREAFRAC(ep) = A/AMAX

    RETURN
  END SUBROUTINE CCD_AREAFRAC
! +SELF.

! ==========================================
    SUBROUTINE MON_SNANA(IFLAG)

! Re-Created Feb 2013
! Driver for monitor tables.
! May 2 2020: add call to TABLE_MARZ


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IFLAG   ! (I) see IFLAG_XXX params

! local var
    LOGICAL   LTMP
    CHARACTER FNAM*10

! ------------ BEGIN ---------
    FNAM = 'MON_SNANA'

    IF( .NOT. USE_TABLEFILE ) RETURN

! ------------------------------------------------------
    IF ( MADE_LCPLOT )  N_SNLC_PLOT = N_SNLC_PLOT + 1
! ------------------------------------------------------

! fill OUTLIER table first so that SNANA table can include summary info
    IF ( OPT_TABLE(ITABLE_OUTLIER) > 0 ) THEN
       CALL TABLE_OUTLIER(IDTABLE_OUTLIER,IFLAG)
    ENDIF

    IF ( OPT_TABLE(ITABLE_SNANA) > 0 ) THEN
       CALL TABLE_SNANA(IDTABLE_SNANA,IFLAG)
    ENDIF


    IF ( MARZFILE_OUT .NE. '' ) THEN
       CALL TABLE_MARZ(IDTABLE_MARZ,IFLAG)
    ENDIF

    RETURN
  END SUBROUTINE MON_SNANA

! ========================================
    SUBROUTINE MON_HUBBLEREF ( IERR )
! 
! Plot theory MAG vs. Z for various Omega_{M,LAM,w}
! to have convenient theory references, and to
! crosscheck DLZ8 function.
! 
! May 29, 2007: add 407 with Omega_MAT=1
! Feb 02, 2013: switch to SNHIST_XXX functions in sntools_output.c.
! Jul 08, 2016: plot curves from z=0-2 (instead of 0-1)
! Apr 23, 2019: add 8th plot with Orad = 1.2E-5
! Oct 23, 2020: use DLMAG_fortC to get DLMAG


    USE SNPAR

    IMPLICIT NONE

    INTEGER IERR ! (O) 0=> OK

! local var

    INTEGER iz, NBZ(2), iplot, hid, LENTIT
    REAL*8  Z0, vPEC, Z0arg(2), zbinsize, ZLIM(2), DLMAG
    character chis*80

! define parameters for theory plots

    INTEGER, PARAMETER :: NPL = 8, NDIM=1

    REAL*8  & 
         OM(NPL)  / 0.3, 0.3,  0.3,  0.3,  0.0,  0.0,  1.0,  0.3 /  & 
        ,OL(NPL)  / 0.7, 0.7,  0.7,  0.0,  0.0,  1.0,  0.0,  0.7 /  & 
        ,OR(NPL)  / 0.0, 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.2E-5/  & 
        ,w0(NPL)  /-1.0,-0.95,-0.95, 0.0,  0.0, -1.0, -1.0, -1.0 /  & 
        ,wa(NPL)  / 0.0,-0.15, 0.15, 0.0,  0.0,  0.0,  0.0,  0.0 /

!     &  ,W(NPL)   /-1.0,-0.9, -0.84, 0.0,  0.0, -1.0, -1.0, -1.0 /
!     &  ,DWDZ(NPL)/ 0.0, 0.0,  0.65, 0.0,  0.0,  0.0,  0.0,  0.0 /

! functions

    REAL*8 DLMAG_FORTC

! ---------------- BEGIN ------------
    IERR = 0

    CALL PRBANNER ( "PLOT THEORY HUBBLE DIAGRAMS for REFERENCE")

! define Z bins for each plot

    ZLIM(1) = 0.0
    ZLIM(2) = 2.0
    NBZ(1)  = 200

    ZLIM(2) = 4.0
    NBZ(1)  = 200
! ccc

    Zbinsize = ( ZLIM(2)-ZLIM(1) ) / float(NBZ(1))

! book histos

    LENTIT = LEN(chis)

    DO iplot = 1, NPL
       write(chis,20) OM(iplot), OL(iplot),  & 
                        OR(iplot)/1.0E-5, w0(iplot), wa(iplot)
20       format('MU vs. z (OM,OL,Or*E5 =', 3(F5.3,',') , ')',  & 
                  2x, 'w0,wa = ',2F6.2, ')'  )
       hid = 400 + iplot
       chis = chis(1:70) // char(0)
       CALL SNHIST_INIT(NDIM, hid, chis, NBZ, zlim(1), zlim(2),  & 
                            LENTIT )
    END DO

! comput Mag vs. Z and fill histograms
    vPEC = 0.0

    DO 100 iz    = 1, NBZ(1)
      Z0 = zbinsize * ( float(iz)-0.5 )
    DO 200 iplot = 1, NPL

! use H0_DEFAULT instead of namelist H0_REF so that plot is always the same.

      DLMAG = DLMAG_fortC(z0,z0, vPEC, H0_DEFAULT,  & 
          OM(iplot), OL(iplot), w0(iplot), wa(iplot) )

      HID = 400 + iplot
      Z0arg(1) = z0  ! use array to aboif gcc 10 warning
      CALL SNHIST_FILL( NDIM, HID, Z0arg, DLMAG )

200   CONTINUE
100   CONTINUE   ! iz


    RETURN
  END SUBROUTINE MON_HUBBLEREF 


! ===============================================
    LOGICAL FUNCTION MJDSELECT ( imjd, cutbit )
! 
! Returns TRUE if this IMJD passes cuts
! CUTBIT = 0 => apply all cuts.
! 
! If cutbit < 0, then all cuts are applied EXCEPT
! for |cutbit|; i.e., the cut on |cutbit| is ignored.
! This is useful to plot variables with cuts on all other
! variables except the one of interest. This allows you to
! see what gets cut.
! 
! Note that if |cutbit| <= CUTBIT_MJD_MARKER,
! then this function really just selects from the
! SN-dependent cuts (like RA, redshfit, etc ...)
! and ignores the epoch index.
! 
! Aug 22, 2007: modify to select NEWMJD instead of epoch.
! Nov 12, 2018: real MJD -> real*8 MJD
! -------------------------------------------

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

! declare subroutine args

    INTEGER  & 
         imjd      &  ! (I) epoch index
        ,cutbit   ! (I) cut to ignore (if negative)

! local var

    INTEGER EPMIN, acutbit   ! abs(cutbit)
    REAL*8    MJD
    INTEGER*8 CUTMASK8_ALL, CUTMASK8_TMP, OVP8

    LOGICAL LCUT_MJD, LCUT_SN

! -------------- BEGIN ------------------

    MJDSELECT =  .FALSE.

! -----------------------------------------------
! first make a few idiot checks

    EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,imjd)
    MJD   = SNGL(SNLC8_MJD(EPMIN))
    if ( MJD .LT. 40000. ) return

! ----------------------------------------

    IF ( CUTBIT .GE. 0 ) THEN
       MJDSELECT =  MJDMASK(imjd) .GT. 0
       RETURN
    ENDIF

! --------------------------------------------------
! Now we have a negative cutbit;
! check all cuts except |cutbit|

    acutbit = abs(cutbit)

    if ( acutbit .GT. NCUTBIT ) then
        write(c1err,60) cutbit
60        format('bad value of CUTBIT = ', I4 )
        CALL MADABORT("MJDSELECT", c1err, "")
    endif

! get cut-mask with all bits set.
! Note "SN" and EPOCHS have different cutmask.

    if ( acutbit .LE. CUTBIT_MJD_MARKER ) then
         CUTMASK8_ALL =  CUTMASK8_SN_ALL
         LCUT_SN    = .TRUE.
         LCUT_MJD   = .FALSE.
    else
         CUTMASK8_ALL =  CUTMASK8_MJD_ALL
         LCUT_SN      = .FALSE.
         LCUT_MJD     = .TRUE.
    endif

! define mask with all cuts except IBIT cut

    CUTMASK8_TMP = IBCLR ( CUTMASK8_ALL, acutbit-1 )

    IF ( LCUT_SN ) then
      OVP8  = IAND ( CUTMASK8_TMP, CUTMASK8_SN )
    ELSE
      OVP8  = IAND ( CUTMASK8_TMP, CUTMASK8_MJD(imjd))
    ENDIF

!  MJD passes if overlap mask equals mask with all
!  cuts EXCEPT |cutbit|

    MJDSELECT = ( OVP8 .EQ. CUTMASK8_TMP )

    RETURN
  END FUNCTION MJDSELECT 




! ==========================================
    SUBROUTINE GET_TRESTVAR(  & 
           NTLIST, TLIST, FLIST, SNRLIST       &  ! (I)
         , CUTWIN_TREST                        &  ! (I)
         , CUTWIN_TRESTMIN, CUTWIN_TRESTMAX    &  ! (I)
         , CUTWIN_TREST2                       &  ! (I)
         , NFILT, TMIN, TMAX                   &  ! (O)
         , NFTMIN, NFTMAX, NFT2                &  ! (O)
         , TGAPMAX, T0GAPMAX                   &  ! (O)
         , SNRSUM     )                       ! (O)

! Created May 25, 2009 by R.Kessler
! Return rest-frame variables related to TREST,
! where TREST = 0 at peak.
! NOTES:
!   - assumes that TLIST is ordered from min to max.
! 
! 
! Jul 01 2024 pass and process SNRLIST
! --------------------


    USE SNPAR

    IMPLICIT NONE

    INTEGER  & 
         NTLIST         &  ! (I) size of TLIST array
        ,FLIST(NTLIST) ! (I) obs-filter list (to count NFTMIN[MAX])

    REAL  & 
         TLIST(NTLIST)          &  ! (I) rest-frame list
        ,SNRLIST(NTLIST)        &  ! (I) SNR list
        ,CUTWIN_TREST(2)        &  ! (I) global Trest-window
        ,CUTWIN_TRESTMIN(2)     &  ! (I) window to count NFTMIN
        ,CUTWIN_TRESTMAX(2)     &  ! (I) idem for NFTMAX
        ,CUTWIN_TREST2(2)      ! (I) generic window

    REAL  & 
         TMIN, TMAX       &  ! (O) min, max Trest
        ,TGAPMAX          &  ! (O) max rest-frame gap (days)
        ,T0GAPMAX        ! (O) idem near peak

    INTEGER  & 
         NFILT            &  ! (O) NFILT in FLIST array
        ,NFTMIN           &  ! (O) Nfilt passing CUTWIN_TRESTMIN
        ,NFTMAX           &  ! (O) Nfilt passing CUTWIN_TRESTMAX
        ,NFT2            ! (O) Nfilt passing CUTWIN_TREST2

! local car

    INTEGER  & 
         i, ORDER, isort  & 
        ,INDEX_SORT(MXEPOCH)  & 
        ,IFILT_obs

    REAL  & 
        Trest, Tlast, TGAP, SNR, SNRSUM, SQSNRSUM  & 
       ,T0GAPWIN(2)

    LOGICAL  & 
         LT0GAP, LTMP  & 
        ,LTMIN(MXFILT_ALL)  & 
        ,LTMAX(MXFILT_ALL)  & 
        ,LT2(MXFILT_ALL)  & 
        ,LFILT(MXFILT_ALL)

! function

    EXTERNAL SORTFLOAT

! ------------------- BEGIN ---------------

! init output args
    Tmin      = +999.
    Tmax      = -999.
    TGAPMAX   = -999.
    T0GAPMAX  = -999.
    NFTMIN    = 0
    NFTMAX    = 0
    NFILT     = 0
    NFT2      = 0
    SQSNRSUM  = 0.0
    SNRSUM    = 0.0

! local stuff
    Tlast     =  999.

! construct T0GAP-window from other cuts
    T0GAPWIN(1)  = CUTWIN_TRESTMIN(2)
    T0GAPWIN(2)  = CUTWIN_TRESTMAX(1)

    DO IFILT_OBS = 1, MXFILT_ALL
      LFILT(IFILT_OBS) = .FALSE.
      LTMIN(IFILT_OBS) = .FALSE.
      LTMAX(IFILT_OBS) = .FALSE.
      LT2(IFILT_OBS)   = .FALSE.
    ENDDO


    ORDER = +1
    CALL sortFloat(NTLIST, TLIST, ORDER, INDEX_SORT )

    DO 444 i   = 1, NTLIST

       isort      = INDEX_SORT(i)
       Trest      = TLIST(isort)
       ifilt_obs  = FLIST(isort)
       snr        = SNRLIST(isort)

       if ( Trest .LT. CUTWIN_TREST(1)   ) GOTO 444
       if ( Trest .GT. CUTWIN_TREST(2)   ) GOTO 444

       IF ( .not. LFILT(ifilt_obs) ) then
          NFILT = NFILT + 1
          LFILT(ifilt_obs) = .TRUE.
       ENDIF

       LTMP = Trest .GE. CUTWIN_TRESTMIN(1) .and.  & 
                Trest .LE. CUTWIN_TRESTMIN(2)
       IF ( LTMP .and.  .not. LTMIN(ifilt_obs) ) then
          NFTMIN = NFTMIN + 1
          LTMIN(ifilt_obs) = .TRUE.
       ENDIF

       LTMP = Trest .GE. CUTWIN_TRESTMAX(1) .and.  & 
                Trest .LE. CUTWIN_TRESTMAX(2)
       IF ( LTMP .and.  .not. LTMAX(ifilt_obs) ) then
          NFTMAX = NFTMAX + 1
          LTMAX(ifilt_obs) = .TRUE.
       ENDIF

       LTMP = Trest .GE. CUTWIN_TREST2(1) .and.  & 
                Trest .LE. CUTWIN_TREST2(2)
       IF ( LTMP .and.  .not. LT2(ifilt_obs) ) then
          NFT2 = NFT2 + 1
          LT2(ifilt_obs) = .TRUE.
       ENDIF

       IF ( LTMP ) THEN
         SQSNRSUM = SQSNRSUM + (SNR*SNR)
       ENDIF

       Tmin = MIN ( Tmin, Trest )
       Tmax = MAX ( Tmax, Trest )

       if ( i .EQ. 1 ) then
         Tlast = Trest
         GOTO 444
       endif

       TGAP = Trest - Tlast
       if ( TGAP .GT. TGAPMAX ) then
           TGAPMAX = TGAP
       endif

       LT0GAP = .FALSE.
       if ( Trest > T0GAPWIN(1) .and. Trest < T0GAPWIN(2) )  & 
             LT0GAP = .TRUE.
       if ( Tlast > T0GAPWIN(1) .and. Tlast < T0GAPWIN(2) )  & 
             LT0GAP = .TRUE.
       if ( Trest > T0GAPWIN(2) .and. Tlast < T0GAPWIN(1) )  & 
             LT0GAP = .TRUE.

       if ( LT0GAP .and. TGAP .GT. T0GAPMAX ) then
         T0GAPMAX = TGAP
       endif

       Tlast = Trest

444   CONTINUE

    if ( SQSNRSUM > 1.0E-8 ) SNRSUM = sqrt(SQSNRSUM)

    RETURN
  END SUBROUTINE GET_TRESTVAR




! ====================================================
    SUBROUTINE DMP_CUTFAIL ( IMJD, idebug )
! 
! Created Jan 31, 2006 by R.Kessler
! 
! Dump CUTBITS that fail for this 'IMJD'
! 
! IDEBUG = 1 : general info for users
! IDEBUG = 2 : include detailed cut-masks for development
! 
! Aug 2007: make special note when telescope is not found.
! Nov 2019: remove telescope logic and CUTBIT_IDTEL


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IMJD   ! (I) sparse SN index and epoch
    INTEGER idebug     ! (I) debug level

! local var

    INTEGER  cid, ibit, ii, EPMIN
    character ccid*(MXCHAR_CCID)
    LOGICAL LCUT, LTEST

! ------------------- BEGIN --------------

    cid  = SNLC_CID
    ccid = SNLC_CCID

    EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,IMJD)

! write MJD-cutmask in hex format

    IF ( IDEBUG .GT. 1 ) then
      write(6,20) CUTMASK8_MJD(IMJD), ccid, IMJD
20      format(T5, 'CUTMASK8_MJD = x', z16.16, 2x,  & 
              'for  SN ', A12, 2x, 'NEWMJD ',I3 )
    ENDIF

! list failed bits for epoch-dependent cuts only

    DO ibit = CUTBIT_MJD_MARKER+1, NCUTBIT

      LTEST = BTEST ( CUTMASK8_MJD_ALL,        ibit-1 )
      LCUT  = BTEST ( CUTMASK8_MJD(IMJD),  ibit-1 )

      if ( LTEST .and.  .not. LCUT ) then
         ii = INDEX(cutvar_name(ibit),':') - 1

         write(6,30)  & 
              IMJD, snlc8_MJD(EPMIN),  & 
              ibit, cutvar_name(ibit)(1:ii)

30         format(T5,'NEWMJD ',I2,1x,'(MJD=',F9.3,')',2x,  & 
               'Failed CUTBIT ', I3, ' : ', A  )

      endif

    ENDDO

    RETURN
  END SUBROUTINE DMP_CUTFAIL 

! =============================================
    SUBROUTINE DMP_SNFAIL()
! 
! Dump utility for SN that fails cuts
! Shows cuts that fail for SN, and cuts that fail for each epoch.
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER imjd, ibit, ii
    LOGICAL LCUT, LTEST

! ---------------- BEGIN ------------

    print*,' '
    print*,' ========== DMP_SNFAIL for CID ', SNLC_CID,  & 
              ' ================= '

    DO 100 ibit = 1, CUTBIT_MJD_MARKER

      LTEST = BTEST ( CUTMASK8_SN_ALL,  ibit-1 )
      LCUT  = BTEST ( CUTMASK8_SN,      ibit-1 )

      if ( LTEST .and.  .not. LCUT ) then
         ii = INDEX(cutvar_name(ibit),':') - 1

         write(6,30) ibit, cutvar_name(ibit)(1:ii)  & 
                ,SNLC_CUTVAR(ibit)  & 
                ,CUTWIN_VAR(1,ibit)  & 
                ,CUTWIN_VAR(2,ibit)
30         format(T5, 'Failed SN cut ',I2,' : ', A,  & 
              ' =', F8.2, ' (cutwin=',F8.1,' to ',G10.3,  ')'  )
      endif

100   CONTINUE

    do imjd = 1, ISNLC_NEWMJD_STORE
       CALL DMP_CUTFAIL ( IMJD, 1 )  ! MJD cut-failures
    enddo

    print*,' ===================================== '
    print*,' '
    RETURN
  END SUBROUTINE DMP_SNFAIL


! ==================================================
    SUBROUTINE SET_CUTMASK()

! 
! Created Jan 31, 2006 by R.Kessler
! 
! After SN data are read in, this routine
! sets CUTMASK8_SN and CUTMASK8_MJD(epoch),
! and then determines which SN and epochs pass/fail cuts.
! 
! To test SN
!   if ( CUTMASK8_SN .EQ. CUTMASK8_SN_ALL ) then
!      use this SN
!   endif
! 
! To test epochs,
! 
!    if (  btest(EPOCHMASK(isn),epoch-1 )  ) then
!       use this epoch
!    endif
! 
!       HISTORY
!    ~~~~~~~~~~~~~
! Feb 27 2018: set cut for ZP in Npe; see CUTWIN_ZP
! Feb 06 2020: add cut for Trest_trueflux
! Feb 20 2020: make sure PHOTPROB cut passes if PHOTPROB < 0
! Mar 16 2021: load SIM_EPPULL to enable cut on PULL
! Jul 01 2024: compute SNRLIST and pass to GET_TRESTVAR
! Nov 25 2025: check SNCUT_NOBS_TREST
! 
! -----------------------------


    USE SNDATCOM
    USE FILTCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

! local var

    INTEGER  & 
         NEWMJD, imjd, NEWMJD_CUTS  & 
        ,EPMIN, EPMAX, ep, ifilt, ifilt_obs  & 
        ,ICUT, i, cid,  NTLIST, MASK, VBOSE  & 
        ,FLIST(MXEPOCH)  & 
        ,NFILT, NFTMIN, NFTMAX, NFT2

    INTEGER*8  CUTMASK8_TMP, OVP8_SN, OVP8_ALL

    LOGICAL LTEST, LCUT, LCUT1, LCUT2

    REAL  TMP, Tmin, Tmax, ZZ, TGAPMAX, T0GAPMAX
    REAL  TLIST(MXEPOCH), SNRLIST(MXEPOCH)
    REAL  TREST, TRUEMAG, PHOTPROB, SNRSUM

! functions
    LOGICAL MJDSELECT, LCIDSELECT
    INTEGER  EVAL_SNCUT_NOBS_TREST
    EXTERNAL EVAL_SNCUT_NOBS_TREST

! ---------------- BEGIN ---------------

    NEWMJD       = ISNLC_NEWMJD_STORE
    NEWMJD_CUTS  = 0

! init common block var output

    CUTMASK8_SN   = 0


    Do imjd = 1, NEWMJD
       CUTMASK8_MJD(imjd) = 0
       MJDMASK(imjd)      = 0  ! init
    End Do

    CID = SNLC_CID

    IF ( SNLC_REDSHIFT > 0. ) THEN
       zz  = 1.0 + SNLC_REDSHIFT
    ELSE
       zz = 1.0  ! avoid crazy Tobsmin/max
    ENDIF

! check privat-var cuts and user cuts from private code
    IF ( .not. PASS_PRIVCUTS ) ISNLC_CUTFLAG_PRIVATE = 0
    IF ( .not. PASS_SIMCUTS  ) ISNLC_CUTFLAG_SIMVAR  = 0

! --------------------------------------------------------
! first evaluate SNLC_CUTVAR array for MJD-dependent
! SNLC_XXX arrays.

    DO 100 imjd = 1, NEWMJD

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,imjd)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,imjd)

! Dec 18 2017: replace cloudAVG with PSF cut. Note that cut is in ARCSEC
       SNEP_CUTVAR(CUTBIT_PSF,imjd)  =  & 
             SNLC_PSF_SIG1(EPMIN) * (SNLC_PIXSIZE * 2.355)

       SNEP_CUTVAR(CUTBIT_ZP,imjd) =  & 
             SNLC_ZEROPT_forCUT(EPMIN)

       SNEP_CUTVAR(CUTBIT_ZPERR,imjd) =  & 
             SNLC_ZEROPT_ERR(EPMIN)      ! Feb 2020

       PHOTPROB = SNLC_PHOTPROB(EPMIN)
       IF ( PHOTPROB < 0.0 ) PHOTPROB = 1.0  ! undefined PHOTPROB passes cut
       SNEP_CUTVAR(CUTBIT_PHOTPROB,imjd)  = PHOTPROB

       SNEP_CUTVAR(CUTBIT_ERRTEST,imjd)  =  & 
             SNLC_FLUXCAL_ERRTEST(EPMIN)

       SNEP_CUTVAR(CUTBIT_SIM_PULL,imjd)  =  & 
             SIM_EPPULL(EPMIN)

       SNEP_CUTVAR(CUTBIT_TREST,imjd)  =  & 
             SNLC_TREST(EPMIN)

       SNEP_CUTVAR(CUTBIT_TOBS,imjd)  =  & 
             SNLC_TOBS(EPMIN)

       SNEP_CUTVAR(CUTBIT_MJD,imjd)  =  & 
             SNGL( SNLC8_MJD(EPMIN) )

! Requiring valid true flux is a little tricky.
! Motivation is for fake light curves on images where a few epochs
! didn't get overlaid for unknown reasons.
       SNEP_CUTVAR(CUTBIT_TREST_TRUEFLUX2,imjd) = 1.0
       IF ( ISJOB_SIM ) THEN
          Trest   = SNLC_TREST(EPMIN)
          TRUEMAG = SIM_EPMAGOBS(EPMIN)
          LCUT1   = (TREST > CUTWIN_TREST_TRUEFLUX(1) .and.  & 
                       TREST < CUTWIN_TREST_TRUEFLUX(2) )
          LCUT2   = (TRUEMAG < 40)
          if ( LCUT1 .and. (.not. LCUT2) ) then
            SNEP_CUTVAR(CUTBIT_TREST_TRUEFLUX2,imjd) = 0.0
          endif
       ENDIF

100   CONTINUE

! ----------------------------------------------------
! Evaluate cutmask for each epoch

    DO 200 imjd = 1, NEWMJD

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,imjd)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,imjd)

      CUTMASK8_TMP = 0
      CUTMASK8_MJD(imjd) = CUTMASK8_TMP

    DO 201 icut  = CUTBIT_MJD_MARKER+1, NCUTBIT

        LTEST = BTEST ( CUTMASK8_MJD_ALL, icut-1 )
        if ( .NOT. LTEST ) goto 201

        tmp = snep_cutvar(icut,imjd)
        LCUT = tmp .GE. cutwin_var(1,icut) .and.  & 
                 tmp .LE. cutwin_var(2,icut)

! set cut-bit if cut is satisfied
        if ( LCUT ) then
           CUTMASK8_TMP = IBSET ( CUTMASK8_TMP, icut-1 )
        endif


201   CONTINUE   ! end ICUT loop

! if all cuts are satisfied for this epoch, then set MJDMASK

      CUTMASK8_MJD(imjd) = CUTMASK8_TMP


      if ( CUTMASK8_TMP .EQ. CUTMASK8_MJD_ALL ) then
         MJDMASK(imjd) = 1
         NEWMJD_CUTS   = NEWMJD_CUTS + 1
      endif

200   CONTINUE  ! end of IMJD loop

! -----------------------------------------------
! store number of epochs passing cuts,

    ISNLC_NEWMJD_CUTS = NEWMJD_CUTS

! -------------------------------------------
! fill SNLC_CUTVAR array for SN-dependent variables
! Note that this is done AFTER epoch-dependent cutmask
! so that NEWMJD_CUTS can be evaluated correctly.

       imjd  = 1

       SNLC_CUTVAR(CUTBIT_CID)       =  & 
              float ( SNLC_CID )

       SNLC_CUTVAR(CUTBIT_SNTYPE)     =  & 
              float ( ISNLC_TYPE )

       SNLC_CUTVAR(CUTBIT_REDSHIFT)  =  & 
             SNLC_REDSHIFT

       SNLC_CUTVAR(CUTBIT_REDSHIFT_ERR)  =  & 
             SNLC_REDSHIFT_ERR

       SNLC_CUTVAR(CUTBIT_RA)  =  & 
             SNGL( SNLC8_RA )
       SNLC_CUTVAR(CUTBIT_DEC)  =  & 
             SNGL( SNLC8_DEC )

       SNLC_CUTVAR(CUTBIT_HOSTSEP)  = SNHOST_ANGSEP(1)

       SNLC_CUTVAR(CUTBIT_PEAKMJD)  =  & 
             SNLC_SEARCH_PEAKMJD

       SNLC_CUTVAR(CUTBIT_NOBS_PREDETECT)  =  & 
             float ( ISNLC_NOBS_PREDETECT )

       SNLC_CUTVAR(CUTBIT_SNRMAX)  =  & 
             SNLC_SNRMAX_FILT(0)
       SNLC_CUTVAR(CUTBIT_SNRMAX2)  =  & 
             SNLC_SNRMAX_FILT(0)  ! temp val; changed below

       SNLC_CUTVAR(CUTBIT_NEPOCH)  =  & 
             float ( ISNLC_NEPOCH_STORE )

       SNLC_CUTVAR(CUTBIT_SEARCH)  =  & 
             SIM_SEARCHEFF_MASK

       SNLC_CUTVAR(CUTBIT_NFIELD)  =  & 
             float ( ISNLC_NFIELD_OVP )

       SNLC_CUTVAR(CUTBIT_MWEBV)  =  & 
              SNLC_MWEBV               ! added May 2012

       SNLC_CUTVAR(CUTBIT_NSEASON_ACTIVE)  =  & 
              float(NSEASON_ACTIVE)

       SNLC_CUTVAR(CUTBIT_REQEP)  =  & 
              float( ISNLC_CUTFLAG_REQEP )   ! added Sep 17 2017

       SNLC_CUTVAR(CUTBIT_PRIVATE)  =  & 
              float( ISNLC_CUTFLAG_PRIVATE )   ! added Nov 3 2014

       SNLC_CUTVAR(CUTBIT_SIMVAR)  =  & 
              float( ISNLC_CUTFLAG_SIMVAR )    ! added Dec 2018

! ---------------------------------------
! set TREST-related variables.
! Warning: this code must go here after evaluating
! cutmask from all other cuts.

! first build TLIST array for epochs that pass MJDSELECT()
     NTLIST = 0
     DO 444 imjd = 1, NEWMJD
       if ( .NOT. MJDSELECT(imjd,0)  ) GOTO 444
       EPMIN     = ISNLC_EPOCH_RANGE_NEWMJD(1,imjd)
       EPMAX     = ISNLC_EPOCH_RANGE_NEWMJD(2,imjd)
       DO EP     = EPMIN, EPMAX
         NTLIST  = NTLIST + 1
         TLIST(NTLIST)   = SNLC_Trest(EP)
         FLIST(NTLIST)   = ISNLC_IFILT_OBS(EP)
         SNRLIST(NTLIST) = SNLC_SNR(EP)
       ENDDO
444    CONTINUE

     CALL GET_TRESTVAR(NTLIST, TLIST, FLIST, SNRLIST   &  ! (I)
             ,CUTWIN_TREST                         &  ! (I)
             ,CUTWIN_TRESTMIN, CUTWIN_TRESTMAX     &  ! (I)
             ,CUTWIN_TREST2                        &  ! (I)
             ,NFILT, TMIN, TMAX                    &  ! (O)
             ,NFTMIN, NFTMAX, NFT2                 &  ! (O)
             ,TGAPMAX, T0GAPMAX, SNRSUM )         ! (O)

     SNLC_TGAPMAX  = TGAPMAX
     SNLC_T0GAPMAX = T0GAPMAX
     SNLC_SNRSUM   = SNRSUM

     SNLC_Trestmin = Tmin
     SNLC_Trestmax = Tmax
     SNLC_TrestRange = Tmax - Tmin    ! Dec 2017

     SNLC_Tobsmin = Tmin * ZZ
     SNLC_Tobsmax = Tmax * ZZ

     ISNLC_NFILT_TRESTMIN  = NFTMIN
     ISNLC_NFILT_TRESTMAX  = NFTMAX
     ISNLC_NFILT_TREST2    = NFT2

! ----------------------------------------------------

     imjd  = 1
     EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,imjd)

! - - - - -
! Nov 24 2025: check NOBS cut in multiple Trest ranges
     VBOSE = 0
     MASK = EVAL_SNCUT_NOBS_TREST(NTLIST, TLIST, VBOSE)
     SNLC_CUTVAR(CUTBIT_MASK_NOBS_TREST) = float(MASK)
!  - - - - -

     SNLC_CUTVAR(CUTBIT_TRESTMIN)  =  & 
             SNLC_TRESTMIN
     SNLC_CUTVAR(CUTBIT_TRESTMAX)  =  & 
             SNLC_TRESTMAX
     SNLC_CUTVAR(CUTBIT_TRESTRANGE)  =  & 
             SNLC_TRESTRANGE

     SNLC_CUTVAR(CUTBIT_TobsMIN)  =  & 
             SNLC_TobsMIN
     SNLC_CUTVAR(CUTBIT_TobsMax)  =  & 
             SNLC_TobsMax

     SNLC_CUTVAR(CUTBIT_TGAPMAX)  =  & 
             SNLC_TGAPMAX
     SNLC_CUTVAR(CUTBIT_T0GAPMAX)  =  & 
             SNLC_T0GAPMAX

     SNLC_CUTVAR(CUTBIT_NFILT_SNRMAX)  =  & 
             float ( ISNLC_NFILT_SNRMAX )

     SNLC_CUTVAR(CUTBIT_NFILT_SNRMAX2)  =  & 
             float ( ISNLC_NFILT_SNRMAX2 )

     SNLC_CUTVAR(CUTBIT_NFILT_TRESTMIN)  =  & 
             float ( ISNLC_NFILT_TRESTMIN )

     SNLC_CUTVAR(CUTBIT_NFILT_TRESTMAX)  =  & 
             float ( ISNLC_NFILT_TRESTMAX )

     SNLC_CUTVAR(CUTBIT_NFILT_TREST2)  =  & 
             float ( ISNLC_NFILT_TREST2 )

! always make sure TREST2 cut passes since the only cut of
! interest in in NFILT_TREST2.
     SNLC_CUTVAR(CUTBIT_TREST2)  = CUTWIN_TREST2(1) + 0.1

     SNLC_CUTVAR(CUTBIT_SNRSUM)  = SNLC_SNRSUM

! be careful with the SNRMAX indices.
     DO i         = 1, NFILT_SNRMAX  ! sparse index for SNRMAX cut
        ifilt_obs = IFILT_SNRMAX(i)  ! absolute filtetr index
        ifilt     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)   ! sparse filt index
        icut      = CUTBIT_OFFSET_SNRMAX + i
        SNLC_CUTVAR(icut) = SNLC_SNRMAX_FILT(ifilt)
     ENDDO

! be careful with the HOST_SBFLUX indices.
     DO i         = 1, NFILT_HOST_SBFLUX  ! sparse index for SNRMAX cut
        ifilt_obs = IFILT_HOST_SBFLUX(i)  ! absolute filtetr index
        ifilt     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)   ! sparse filt index
        icut      = CUTBIT_OFFSET_SBFLUX + i
        SNLC_CUTVAR(icut) =  & 
               SNHOST_SBFLUXCAL(ifilt)
     ENDDO

! -----------------------------------------------------
! After epochs, evaluate SN cut-mask
!  (epoch-independent => use 1st epoch of SNLC_CUTVAR)
! This part comes AFTER epochmask so that we can cut
! on NEWMJD_CUTS

       CUTMASK8_TMP = 0

    DO 300 icut  = 1, CUTBIT_MJD_MARKER

        LTEST = BTEST (CUTMASK8_SN_ALL , icut-1 )
        if ( .NOT. LTEST ) GOTO 300

        tmp = snlc_cutvar(icut)
        LCUT = tmp .GE. cutwin_var(1,icut) .and.  & 
                 tmp .LE. cutwin_var(2,icut)

! special check for CID cut. This cut is probably redundant
! with the cut in PARSE_CID.

        if ( icut .EQ. CUTBIT_CID ) then
          LCUT = LCIDSELECT(SNLC_CID,SNLC_CCID,SNLC_NAME_IAUC)
        endif

! set cut-bit if cut is satisfied

        if ( LCUT ) then
           CUTMASK8_TMP = IBSET ( CUTMASK8_TMP, icut-1 )
           NACCEPT_CUT(icut) = NACCEPT_CUT(icut) + 1
        endif

300   CONTINUE   ! end ICUT loop

    CUTMASK8_SN = CUTMASK8_TMP

! set SN cut-logical

    IF ( CUTMASK8_SN .EQ. CUTMASK8_SN_ALL  ) THEN
       LSNCUTS       = .TRUE.
       CUTFLAG_SNANA = IBSET(CUTFLAG_SNANA,0)
       N_SNLC_CUTS   = N_SNLC_CUTS + 1
       if (NSPECTRUM > 0) N_SNLC_SPEC = N_SNLC_SPEC + 1

       if(SNHOST_ZSPEC(1) > 0.) N_SNHOST_ZSPEC = N_SNHOST_ZSPEC+1
       if(SNHOST_ZPHOT(1) > 0.) N_SNHOST_ZPHOT = N_SNHOST_ZPHOT+1

	 N_MASK_zSOURCE_LC_CUTS(ISNLC_zSOURCE) =  & 
         N_MASK_zSOURCE_LC_CUTS(ISNLC_zSOURCE) + 1

    ENDIF

! Jan 3, 2013
! increment Number of SN passing each incremental cut,
! and versus type. This info is dumped at the end of the
! job to trace where SN are lost.

    IF ( ISNLC_TYPE < 0 ) RETURN

    IF ( ISNLC_TYPE > MXTYPE ) THEN
      write(C1err,661) ISNLC_TYPE, MXTYPE
661     format('TYPE=',I5,' exceeds bound: MXTYPE=',I4 )
      C2err = 'Check CID = ' // SNLC_CCID
      CALL MADABORT("SET_CUTMASK", c1err, c2err )
    ENDIF

    CUTMASK8_TMP = 0
    DO icut = 1, NCUTBIT_SNLC
       CUTMASK8_TMP = IBSET(CUTMASK8_TMP,icut-1)
       OVP8_ALL     = IAND(CUTMASK8_TMP,CUTMASK8_SN_ALL)
       OVP8_SN      = IAND(CUTMASK8_SN,OVP8_ALL)
       IF ( OVP8_SN .EQ. OVP8_ALL ) THEN
          NPASSCUT_INCREMENT(ISNLC_TYPE,ICUT) =  & 
            NPASSCUT_INCREMENT(ISNLC_TYPE,ICUT) + 1

! increment all types with ITYPE=-1 (Aug 2013)
          NPASSCUT_INCREMENT(-1,ICUT) =  & 
            NPASSCUT_INCREMENT(-1,ICUT) + 1
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE SET_CUTMASK


! ==============================================
    LOGICAL FUNCTION LCIDSELECT(CID,CCID,NAME_IAUC)
! 
! Created Dec 4, 2012 by R.Kessler
! Returns TRUE if current CID/CCID is accepted by
! any of
!  CUTWIN_CID     (integer window of CID to accept)
!  SNCID_LIST     (integer list of CID to accept)
!  SNCCID_LIST    (char list of CCID to accept)
! 
! and is NOT rejected by
!  SNCID_IGNORE       (integer list of CID to ignore)
!  SNCCID_IGNORE_ALL  (char list of CCID to ignore)
! 
! This function replaces much of the obsolete CIDMASK(CID) array
! to avoid declaring such a large array.
! 
! 
! Oct 10 2020: fix bug checking SNCCID_IGNORE_ALL
! 
! Jan 08 2025: set INDEX_CID_MATCH inside check for SNCID_LIST in case SNCID_LIST
!              has overlap with CID in SNCID_LIST_FILE.
! 
! ---------------------------


    USE SNPAR
    USE SNCUTS
    USE SNLCINP_NML
    USE INTERPCM
    USE SNLCCOM

    IMPLICIT NONE

! input function args
    INTEGER CID
    CHARACTER CCID*(*), NAME_IAUC*(*)

! local var
    INTEGER i, LEN, LENCCID, LENIAUC, ISN_MATCH
    REAL XCID
    CHARACTER CCID_forC*(MXCHAR_CCID), CCIDTMP*(MXCHAR_CCID)

! functions
    LOGICAL   STRINGMATCH
    INTEGER   MATCH_CIDLIST_EXEC
    EXTERNAL  MATCH_CIDLIST_EXEC
! ------------------ BEGIN ---------------

    LCIDSELECT = .FALSE.
    LENCCID    = INDEX(CCID//' ' , ' ') - 1
    CCID_forC  = CCID(1:LENCCID) // char(0)
    INDEX_CID_MATCH = -9
! ----------------------------
! first check IGNORE Lists

    i = 1
    DO 30 WHILE ( SNCID_IGNORE(i) .GE. 1 )
      if ( CID .EQ. SNCID_IGNORE(i) ) GOTO 800
      i = i + 1
30    CONTINUE

    DO 40 i = 1, NCCID_IGNORE
      LEN   = INDEX(SNCCID_IGNORE_ALL(i),' ') - 1
      IF ( LEN .NE. LENCCID ) GOTO 40
      if ( CCID(1:LEN) .EQ. SNCCID_IGNORE_ALL(i)(1:LEN)) GOTO 800
40    CONTINUE

! -----------------------------
! Check integer window

    XCID = FLOAT(CID)
    IF ( XCID .GE. CUTWIN_CID(1) .and.  & 
           XCID .LE. CUTWIN_CID(2) ) then
      LCIDSELECT = .TRUE.
      GOTO 800
    ENDIF

! Since user-defined list is small, use brute-force looping to match

    IF ( NCID_LIST > 0 .OR. NCCID_LIST > 0 ) THEN
      i = 1
      DO 130 WHILE ( SNCID_LIST(i) .GE. 1 )
        if ( CID .EQ. SNCID_LIST(i) ) then
          LCIDSELECT = .TRUE.
          INDEX_CID_MATCH = MATCH_CIDLIST_EXEC(CCID_forC,LENCCID) ! in case we need fit prior from file
! xxx mark delete May 28 2025          GOTO 800
        endif
        i = i + 1
130     CONTINUE

! To do: combine SNCID_LIST and SNCCID_LIST into single SNCCID_LIST to avoid duplicate code logic !!!
!        make sure to start with clean new git tag in case things go really bad.

      i = 1
      DO 131 WHILE ( SNCCID_LIST(i) .NE. ' ' )
        if ( CCID .EQ. SNCCID_LIST(i) ) then
          LCIDSELECT = .TRUE.
          INDEX_CID_MATCH = MATCH_CIDLIST_EXEC(CCID_forC,LENCCID) ! in case we need fit prior from file
        endif
        i = i + 1
131     CONTINUE

      ! May 2025: always bail after checking user list, even if there is a file that defines prior.
      GOTO 800

    ENDIF  ! end user-define lists

! - - - - -
! If we get here, there is no user list and only a list in a file.

    ISN_MATCH = MATCH_CIDLIST_EXEC(CCID_forC,LENCCID)
    if ( ISN_MATCH  < 0 ) then   ! try IAUC name
       LENIAUC    = INDEX(NAME_IAUC,' ') - 1
       CCID_forC  = NAME_IAUC(1:LENIAUC) // char(0)
       ISN_MATCH  = MATCH_CIDLIST_EXEC(CCID_forC,LENIAUC) ! use hash table
       if ( ISN_MATCH .GE. 0 ) CCID = NAME_IAUC ! write IAUC to output table
    endif
    INDEX_CID_MATCH = ISN_MATCH     ! load global
    LCIDSELECT = (ISN_MATCH .GE. 0)
    if ( LCIDSELECT ) GOTO 800

! - - - - -
! finally, check list of CIDs on the MJD-interpolation list;
! see &SNLCINP namelist input SNMJD_LIST_FILE

    DO 240 i = 1, N_INTERP_MJDLIST
      CCIDTMP = INTERP_CCIDLIST(i)

      IF ( STRINGMATCH(CCIDTMP,CCID) ) THEN
         LCIDSELECT = .TRUE.
         GOTO 800
      ENDIF

! Dec 2015: also check match with IAUC name
      IF ( STRINGMATCH(CCIDTMP,NAME_IAUC) ) THEN
         CCID       = NAME_IAUC  ! so that IAUC is written to output
         LCIDSELECT = .TRUE.
         GOTO 800
      ENDIF

240   CONTINUE

800   CONTINUE

    RETURN
  END FUNCTION LCIDSELECT


! ====================================
    LOGICAL FUNCTION STRINGMATCH(STR1,STR2)

! Created Dec 2015
! Return True of STR1=STR2


    CHARACTER STR1*(*), STR2*(*)
    INTEGER   LEN1, LEN2

! ------------------ BEGIN -------------

    STRINGMATCH = .FALSE.
    LEN1 = INDEX(STR1,' ') -1
    LEN2 = INDEX(STR2,' ') -1
    IF ( LEN1 .NE. LEN2 ) RETURN

    IF ( STR1(1:LEN1) .EQ. STR2(1:LEN2) ) THEN
      STRINGMATCH = .TRUE.
    ENDIF

    RETURN
  END FUNCTION STRINGMATCH

! ======================
    DOUBLE PRECISION FUNCTION SMOOTH_STEPFUN8(SEP,SEPMAX)
! 
! Created Sep 8, 2009 by R.Kessler
! 
! Define smooth function that goes from 0 to 1 between
! -SEPMAX and SEPMAX, and returns 0.5 at SEP=0.
! Return value of function at SEP.
! Function is atan.
! 
! NOTE: this should be moved to sntools.c (C code)
! 
! --------------------------


! +CDE,SNPAR.

    REAL*8 SEP, SEPMAX  ! (I)
    REAL*8 TAU, F
! ----------- BEGIN ----------

    IF ( SEP .GT. SEPMAX ) THEN
       SMOOTH_STEPFUN8 = 1.0
       RETURN
    ENDIF
    IF ( SEP .LT. -SEPMAX ) THEN
       SMOOTH_STEPFUN8 = 0.0
       RETURN
    ENDIF

    TAU  = .1 * SEPMAX
    F    = 0.5 * (1. + ATAN(SEP/TAU)/ATAN(SEPMAX/TAU) )
    SMOOTH_STEPFUN8 = F

    RETURN
  END FUNCTION SMOOTH_STEPFUN8


! ======================================
    INTEGER FUNCTION FIRSTBIN_INTERP (  & 
         VALUE       &  ! (I) value to find bin for
        ,RANGE       &  ! (I) allowed range for VALUE
        ,BINSIZE     &  ! (I) binsize for RANGE
          )

! 
! Return first bin to use for interpoliation.
! Asumes interp is done with 3 bins, where 2nd
! bin contains VALUE.
! 
! ------------------------------


    REAL  VALUE, RANGE(2), BINSIZE  ! (I)

! local args

    INTEGER ibin_cen, ibin_start, NBIN
    REAL dif

! ----------------- BEGIN ------------

    FIRSTBIN_INTERP  = -9

    dif      = value - range(1)
    ibin_cen = int ( (dif / binsize) + 0.5 ) + 1

    dif      = range(2) - range(1) + binsize/10000.
    NBIN     = int( dif / binsize ) + 1

! find starting bin for grid interpolation

    if ( ibin_cen .LE. 1 ) then
        ibin_start = 1
    else if ( ibin_cen .GE. NBIN ) then
        ibin_start = NBIN - 2
    else
        ibin_start = ibin_cen - 1   ! nominal case
    endif

    FIRSTBIN_INTERP  = ibin_start

    RETURN
  END FUNCTION FIRSTBIN_INTERP 


! ==================================================
    LOGICAL FUNCTION LCUTVAR ( var, icut, VBOSE )

    USE SNPAR
    USE SNCUTS
    USE SNLCINP_NML

    IMPLICIT NONE

    REAL    VAR   ! (I) variable to  test cutwin
    INTEGER ICUT  ! (I) points to cutwin
    LOGICAL VBOSE ! (I) T=> print PASS or FAIL for each cut

! local var

    INTEGER ITMP
    CHARACTER CTMP*6

! -------------- BEGIN --------------

    LCUTVAR = var .GE. cutwin_var(1,icut) .and.  & 
                var .LE. cutwin_var(2,icut)

    IF ( VBOSE ) THEN
        if ( LCUTVAR ) then
           CTMP = 'PASSES'
        else
           CTMP = 'FAILED'
        endif

        ITMP = INDEX(cutvar_name(icut),':') - 1
        ITMP = MIN(ITMP,14)
        write(6,20) cutvar_name(icut)(1:ITMP), var, CTMP
20        format(T10,A18,' = ', G10.3, 1x, A, ' CUT.' )
    ENDIF

    RETURN
  END FUNCTION LCUTVAR 

! ====================================
    SUBROUTINE DMPTRACE ( banner )
    character banner*(*)
    print*,'  XXXXX ', BANNER
    CALL FLUSH(6)
    RETURN
  END SUBROUTINE DMPTRACE 


! ========================================
    INTEGER FUNCTION IPLOT_JOBSPLIT(NPLOT)
! 
! Created Sep 28, 2012 by R.Kessler
! 
! For single job, returns IPLOT = NPLOT is the
! incremental (sparse) index. For non-trivial JOBSPLIT input,
!      IPLOT = (NPLOT-1) * NSPLIT_TOT + ISPLIT
! so that the IPLOT are unique among each split job.
! Needed so that snana#fitres plot-macro works on merged
! histo files from split_and_fit.pl script.
! 
! -----------------------------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER NPLOT  ! (I) sequential plot number for this job

    INTEGER NJOBTOT, IJOB, IPLOT

! -------------- BEGIN ---------------

! set default for non-split jobs.
    IPLOT = NPLOT

! check for internal split (mainly for FITS format)
    NJOBTOT = JOBSPLIT(2) ! total number of split jobs
    IJOB    = JOBSPLIT(1) ! do IJOB of NJOBTOT
    IF ( NJOBTOT .GT. 1 ) THEN
      IPLOT = (NPLOT-1) * NJOBTOT + IJOB
    ENDIF

! check for text-format where SPLIT has been done externally
    NJOBTOT = JOBSPLIT_EXTERNAL(2) ! total number of split jobs
    IJOB    = JOBSPLIT_EXTERNAL(1) ! do IJOB of NJOBTOT
    IF ( NJOBTOT .GT. 1 ) THEN
      IPLOT = (NPLOT-1) * NJOBTOT + IJOB
    ENDIF

! load output
    IPLOT_JOBSPLIT = IPLOT

    RETURN
  END FUNCTION IPLOT_JOBSPLIT

! =============================================
    SUBROUTINE GET_SIMNAME_TYPE(name)

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    CHARACTER NAME*(*)  ! (I)
    INTEGER LL

! --------- BEGIN ---------
    NAME = ''
    IF ( .NOT. LSIM_SNANA ) RETURN
    LL = INDEX(SIMNAME_TYPE,' ') - 1
    NAME = SIMNAME_TYPE(1:LL) // char(0)  ! Ia, Ib, IIN, etc ...
    RETURN
  END SUBROUTINE GET_SIMNAME_TYPE


! =============================
    LOGICAL FUNCTION DOPLOT_SNLC()

! Created Aug 30 2013
! Returns T to make light curve plot for this SN.
! 
! ----------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, LENTMP
    LOGICAL DOPLOT

! ------------ BEGIN -------------

    DOPLOT = .FALSE.

    if ( OPT_TABLE(ITABLE_SNLCPAK) == 0 ) goto 888  ! Sep 8 2014

    if ( N_SNLC_PLOT < MXLC_PLOT ) THEN
       DOPLOT = .TRUE.
       GOTO 888
    endif

! check list of CCIDs to plot

    DO 100 i = 1, NCCID_PLOT
      LENTMP = INDEX( SNCCID_PLOT(i)//' ', ' ') - 1

      IF ( LENTMP .NE. ISNLC_LENCCID ) GOTO 100

      if ( SNCCID_PLOT(i) .EQ. SNLC_CCID ) then
        DOPLOT = .TRUE.
        GOTO 888
      endif

100   CONTINUE

888   CONTINUE
    DOPLOT_SNLC = DOPLOT

    RETURN
  END FUNCTION DOPLOT_SNLC

! =============================
    LOGICAL FUNCTION DOPLOT_SPEC()

! Created May 30 2020
! Returns T to make spectra plot(s) for this SN
! 
! ----------


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER i, LENTMP
    LOGICAL DOPLOT

! ------------ BEGIN -------------

    DOPLOT = .FALSE.
    if ( OPT_TABLE(ITABLE_SPECPAK) == 0 ) goto 888

    if ( N_SNLC_PLOT < MXLC_PLOT ) THEN
       DOPLOT = .TRUE.
       GOTO 888
    endif

! check list of CCIDs to plot

    DO 100 i = 1, NCCID_PLOT
      LENTMP = INDEX( SNCCID_PLOT(i)//' ', ' ') - 1

      IF ( LENTMP .NE. ISNLC_LENCCID ) GOTO 100

      if ( SNCCID_PLOT(i) .EQ. SNLC_CCID ) then
        DOPLOT = .TRUE.
        GOTO 888
      endif

100   CONTINUE

888   CONTINUE
    DOPLOT_SPEC = DOPLOT

    RETURN
  END FUNCTION DOPLOT_SPEC

! -----------------------------------------
#if defined(SNANA)
    INTEGER FUNCTION SNANA_GET_NLCPLOT()

! Created May 17 2014 by R.Kessler
! Return NFIT_PER_SN = 0 or 1 based on SNLCINP options.
! For snana.exe job only.
! Sep 7 2014: use OPT_TABLE.


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER NLCPLOT, OPT
! --------------- BEGIN ------------

    NLCPLOT = 0  ! default -> no LC plotted
    OPT = OPT_TABLE(ITABLE_SNLCPAK)
    if ( OPT > 0 ) NLCPLOT = 1

! return function value
    SNANA_GET_NLCPLOT = NLCPLOT

    RETURN
  END FUNCTION SNANA_GET_NLCPLOT
#endif

! ==============================
#if defined(SNANA)
    SUBROUTINE SNLCPLOT()

! Created Feb 2013 by R.Kessler
! Called only by snana.exe to
!  * create sub dir
!  * pack light curves and PKMJD fit-curve
!  * cdtopdir
! 
! Call wrappers in sntools_output.c to that there
! no native calls to CERNLIB or ROOT.
! 
! 
! Apr 11 2019:  MXEP_SNLCPAK -> 4*MXEPOCH (was 10*MXEPOCH)
! Jan 16 2020:  pass sim fluxes (see VSIMFLUX)
! -----------------------------------------------------------

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PKMJDCOM

    IMPLICIT NONE

    CHARACTER  CCID*(MXCHAR_CCID), TEXT_forC*80

    LOGICAL  OVMODEL_ANYFUN, LTMP, DOPLOT_FILT(MXFILT_OBS)
    INTEGER  & 
         LENCCID, LENTEXT, NEWMJD, EPMIN, EPMAX, EP  & 
        ,IFILT_OBS, IFILT, ipar, NBT, i, NFILT

    REAL*8  & 
          Z, Z1, Tobs, Trest, MJD, MJD_PLOT, TOBS_PLOT  & 
         ,FLUX_DATA,  FLUXERR_DATA, FLUX_MODEL  & 
         ,SQDIF, SQERR, CHI2, TMIN, TMAX, DT  & 
         ,XVAL8(NPAR_ANYLC)

! SNLCPAK arrays

    INTEGER, PARAMETER :: MXEP_SNLCPAK = 4*MXEPOCH 

    REAL*8  & 
         VMJD(MXEP_SNLCPAK)  & 
        ,VTOBS(MXEP_SNLCPAK)  & 
        ,VFLUX(MXEP_SNLCPAK)  & 
        ,VFLUX_ERR(MXEP_SNLCPAK)  & 
        ,VSIMFLUX(MXEP_SNLCPAK)  & 
        ,VCHI2(MXEP_SNLCPAK)  & 
        ,VREJECT(MXEP_SNLCPAK)  & 
        ,VERRCALC(MXEP_SNLCPAK)  & 
        ,VDUMERR(MXEP_SNLCPAK)   &  ! all zeros
! 
        ,VBAND_NDOF(MXFILT_OBS)  & 
        ,VBAND_CHI2(MXFILT_OBS)  & 
        ,VBAND_PKFLUX(MXFILT_OBS)  & 
        ,VBAND_PKFLUX_ERR(MXFILT_OBS)  & 
        ,VBAND_PKMJD(MXFILT_OBS)  & 
        ,VBAND_PKMJD_ERR(MXFILT_OBS)

    INTEGER  & 
         NOBS  & 
        ,VFILTOBS(MXEP_SNLCPAK)

! functions
    EXTERNAL  ANYLCFUN
    REAL*8    ANYLCFUN

    LOGICAL DOPLOT_SNLC

! ---------------- BEGIN --------------

! check if we should make the plot
    IF ( .NOT. DOPLOT_SNLC() ) RETURN

! xxx mark dele      CALL SNLCPAK_NFIT(1)      ! May 11 2014

    LENCCID = ISNLC_LENCCID
    CCID = SNLC_CCID(1:LENCCID) // char(0)

    OVMODEL_ANYFUN = (SNLC_SNANAFIT_PEAKMJD > 10.0) !show best peak-MJD fit
    NFILT   = NFILTDEF_SURVEY

    write(6,10) SNLC_CCID(1:LENCCID)
 10   format(T8,'SNLCPLOT: pack CID=',A,' for plotting.')
    call flush(6)

! create subdir
    CALL MAKEDIR_OUTPUT(CCID, SNLC_CID, LENCCID)

! prepare text string(s) to display on plot
    LENTEXT = LEN(TEXT_forC)


    write(TEXT_forC,20) ISNLC_TYPE, SNLC_REDSHIFT, char(0)
 20   format('TYPE = ', I3, 3x,'z=',F5.3, 3x, A )
    CALL SNLCPAK_DISPLAYTEXT(CCID, TEXT_forC, LENCCID, LENTEXT )

    IF ( OVMODEL_ANYFUN ) THEN
       TEXT_forC = 'Fit  ANYFUN  to estimate peak MJD' // char(0)
       CALL SNLCPAK_DISPLAYTEXT(CCID, TEXT_forC, LENCCID, LENTEXT )
    ENDIF

! flag which filter(s) to plot
    DO IFILT = 1, NFILT
      VBAND_NDOF(ifilt)   = 0.0           ! init this array
      VBAND_CHI2(ifilt)   = 0.0           ! init this array
      VBAND_PKFLUX(ifilt) = 0.0
      VBAND_PKMJD(ifilt)  = 0.0

      DOPLOT_FILT(IFILT) = .TRUE.  ! default is all
      IFILT_OBS = IFILTDEF_MAP_SURVEY(ifilt)
      LTMP = ( PKMJD_FIT(ifilt_obs) .GT. 100.0 )
      IF ( OVMODEL_ANYFUN .and. .NOT. LTMP ) then
        DOPLOT_FILT(IFILT) = .FALSE.
      ENDIF

! set default PKMJD
      VBAND_PKMJD(IFILT)     = SNLC_SEARCH_PEAKMJD
      VBAND_PKMJD_ERR(IFILT) = 0.0

    ENDDO

! -------------------------------
    Z    = SNLC_REDSHIFT
    Z1   = 1 + Z
! -----------------------------

    NOBS = 0
    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)

      MJD   = SNLC8_MJD(EPMIN)
      Tobs  = MJD - SNLC_SEARCH_PEAKMJD
      Trest = Tobs / Z1

      DO 201 EP = EPMIN, EPMAX

        FLUX_DATA     = SNLC_FLUXCAL(ep)
        FLUXERR_DATA  = SNLC_FLUXCAL_ERRTOT(ep)
        if ( FLUXERR_DATA .LE. 0.0 ) GOTO 201

        IFILT_OBS = ISNLC_IFILT_OBS(ep)
        IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

        IF ( .NOT. DOPLOT_FILT(IFILT) ) GOTO 201

        VBAND_NDOF(ifilt) = VBAND_NDOF(ifilt) + 1.0 ! note this is real*8
        NOBS = NOBS + 1

        IF ( NOBS .GT. MXEP_SNLCPAK ) THEN
          write(C1ERR,661) NOBS, MXEP_SNLCPAK
 661        format('NOBS=',I5,' exceeds bound of MXEP_SNLCMAX=',I5)
          c2err = 'Increase bound for data array.'
          CALL MADABORT('SNLCPLOT(SNANA)', C1ERR, C2ERR)
        ENDIF

        VMJD(NOBS)      = MJD
        VTOBS(NOBS)     = TOBS
        VFLUX(NOBS)     = SNLC_FLUXCAL(ep)
        VFLUX_ERR(NOBS) = SNLC_FLUXCAL_ERRTOT(ep)
        VERRCALC(NOBS)  = SNLC_FLUXCAL_ERRCALC(ep)
        VSIMFLUX(NOBS)  = SIM_EPFLUXCAL(ep)   ! Jan 2020

! Sep 7 2022: check options to fold LC within a single cylce.
        if ( MJDPERIOD_PLOT > .01 .and. NOBS > 1 ) then
           DT = MJD - VMJD(1) - MJDSHIFT_PLOT
           VMJD(NOBS)  = MOD(DT,MJDPERIOD_PLOT) + VMJD(1)

           DT = TOBS - VTOBS(1) - MJDSHIFT_PLOT
           VTOBS(NOBS) = MOD(DT,MJDPERIOD_PLOT) + VTOBS(1)
        endif

        VFILTOBS(NOBS)  = IFILT_OBS
        VREJECT(NOBS)   = DBLE( 1 - ISNLC_SNRECON_USE(ep) )
        VDUMERR(NOBS)   = 0.0  ! dummy arguments

! fetch model value and VCHI2 (if we fit the "ANYLC" function)

        if ( OVMODEL_ANYFUN ) then
          do ipar = 1, NPAR_ANYLC
            XVAL8(ipar) =  & 
              SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar)
          enddo

          FLUX_MODEL  = ANYLCFUN(MJD,XVAL8)
          SQDIF       = (FLUX_DATA - FLUX_MODEL)**2
          SQERR       = FLUXERR_DATA**2
          CHI2        = SQDIF/SQERR
          VCHI2(NOBS) = CHI2
          VBAND_CHI2(IFILT) = VBAND_CHI2(IFILT) + CHI2

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          if ( ifilt_obs == 4 .and. MJD < 56600. ) then
             write(6,621) MJD, TOBS, FLUX_DATA,  & 
                   FLUX_MODEL, FLUXERR_DATA
 621           format(' xxx MJD=',F9.3, '  TOBS=',F6.1,  & 
                   '  FLUX(D,M)=',2G10.3, '  FLUX_ERR=', G8.3 )
          endif
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        endif  ! end of OVMODEL

201     CONTINUE  ! EP
200   CONTINUE   ! NEWMJD

! ------------------------------------------

! store data fluxes
    CALL SNLCPAK_DATA(CCID, NOBS, VMJD, VTOBS, VFLUX, VFLUX_ERR,  & 
           VFILTOBS, SNLCPAK_EPFLAG_FLUXDATA, LENCCID)

! store REJECT flags
    CALL SNLCPAK_DATA(CCID, NOBS, VMJD, VTOBS, VREJECT, VDUMERR,  & 
           VFILTOBS, SNLCPAK_EPFLAG_REJECT, LENCCID)

! store data-FITFUN chi2
    IF ( OVMODEL_ANYFUN ) THEN
      CALL SNLCPAK_DATA(CCID, NOBS, VMJD, VTOBS, VCHI2, VDUMERR,  & 
           VFILTOBS, SNLCPAK_EPFLAG_CHI2, LENCCID)
    ENDIF

! check for sim fluxes (Jan 2020)
    IF ( ISJOB_SIM ) THEN
      CALL SNLCPAK_DATA(CCID, NOBS, VMJD, VTOBS,  & 
           VSIMFLUX, VFLUX_ERR,  & 
           VFILTOBS, SNLCPAK_EPFLAG_FLUXSIM, LENCCID)
    ENDIF

! ------------------------------------------
    IF ( .NOT. OVMODEL_ANYFUN ) GOTO 500
! ------------------------------------------

! plot smooth best-fit ANYFUN function in 2-day bins.

    TMIN = -40.0
    TMAX = 140.0
    NBT  = int( (TMAX-TMIN) / DTOBS_MODEL_PLOT )
    DT   = (TMAX - TMIN) / float(NBT)
    NOBS = 0

    DO 300 ifilt = 1, NFILT
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)

       VBAND_PKFLUX(ifilt)     = DBLE( PKFLUX_FIT(ifilt_obs) )
       VBAND_PKFLUX_ERR(ifilt) = DBLE( PKFLUX_ERR(ifilt_obs) )

       IF ( .not. EXIST_FILT(ifilt) ) GOTO 300

       DO 320 i = 1, NBT
          Tobs  = TMIN + dT * ( float(i) - 0.5 )
          MJD   = Tobs + SNLC_SEARCH_PEAKMJD
          do ipar = 1, NPAR_ANYLC
             XVAL8(ipar) = SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar)
          enddo

          NOBS = NOBS + 1

          IF ( NOBS .GT. MXEP_SNLCPAK ) THEN
             write(C1ERR,661) NOBS, MXEP_SNLCPAK
             c2err = 'Increase bound for best-fit (ANYFUN) array.'
             CALL MADABORT('SNLCPLOT(SNANA)', C1ERR, C2ERR)
          ENDIF

          FLUX_MODEL       = ANYLCFUN(MJD,XVAL8)

! avoid plotting artifacts with very small numbers (Oct 29 2014)
          IF( FLUX_MODEL < 1.0E-20 ) FLUX_MODEL = 0.0

          VMJD(NOBS)       = MJD
          VTOBS(NOBS)      = TOBS
          VFLUX(NOBS)      = FLUX_MODEL
          VFLUX_ERR(NOBS)  = 0.0
          VFILTOBS(NOBS)   = IFILT_OBS

320     CONTINUE  ! epoch-bin

300   CONTINUE  ! ifilt


! -------------------------------------------------

! store best-fit ANYLCFUN

    CALL SNLCPAK_DATA(CCID, NOBS, VMJD, VTOBS, VFLUX, VFLUX_ERR,  & 
           VFILTOBS, SNLCPAK_EPFLAG_FITFUN, LENCCID)

! store filter-dependent quantities.
    VMJD(1)  = -9999.0
    VTOBS(1) = -9999.0  ! dummy for unused TOBS arg

    CALL SNLCPAK_DATA(CCID, NFILT, VMJD, VTOBS, VBAND_NDOF, VDUMERR,  & 
           IFILTDEF_MAP_SURVEY, SNLCPAK_BANDFLAG_NDOF, LENCCID)

    CALL SNLCPAK_DATA(CCID, NFILT, VMJD,VTOBS, VBAND_CHI2, VDUMERR,  & 
           IFILTDEF_MAP_SURVEY, SNLCPAK_BANDFLAG_CHI2, LENCCID)

    CALL SNLCPAK_DATA(CCID, NFILT, VTOBS, VMJD,  & 
         VBAND_PKFLUX, VBAND_PKFLUX_ERR,  & 
         IFILTDEF_MAP_SURVEY, SNLCPAK_BANDFLAG_PKFLUX, LENCCID)

! --------------------------------------------------
! fill output structure with above; here iy goes to a disk file.

 500  CONTINUE

    CALL SNLCPAK_NFIT(1)      ! move from abvove, Oct 23 2025

    CALL SNLCPAK_DATA(CCID, NFILT, VMJD, VTOBS,  & 
          VBAND_PKMJD,VBAND_PKMJD_ERR,  & 
          IFILTDEF_MAP_SURVEY, SNLCPAK_BANDFLAG_PKMJD, LENCCID)

    CALL SNLCPAK_FILL(CCID, LENCCID)

! back to topdir
    CALL CDTOPDIR_OUTPUT(STDOUT_UPDATE)

! set global flag.
    MADE_LCPLOT = .TRUE.

    RETURN
  END SUBROUTINE SNLCPLOT
#endif


! ==============================================
    SUBROUTINE SPECPLOT()

! Apr 2019
! Prepare spectra table for plotting (analog of SNLCPLOT)
! Works only for TEXT output; not for ROOT
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE SPECCOM

    IMPLICIT NONE

    INTEGER LENCCID, LENTEXT, ispec
    CHARACTER  CCID_forC*(MXCHAR_CCID)
    LOGICAL DOPLOT_SPEC

    EXTERNAL SPECPAK_DATA
! ----------------- BEGIN ---------------

! check if we should make the plot
    IF ( .NOT. DOPLOT_SPEC() ) RETURN

    IF ( MARZFILE_OUT .NE. '' ) RETURN

    LENCCID = ISNLC_LENCCID
    CCID_forC = SNLC_CCID(1:LENCCID) // char(0)

    write(6,10) SNLC_CCID(1:LENCCID)
 10   format(T8,'SPECPLOT: pack CID=',A,' for plotting.')
    call flush(6)

! load each spectrum into SPECPAK C-function util

    DO 100 ispec = 1, NSPECTRUM
       CALL RDSPEC_DRIVER(ispec)  ! Apr 5 2021

       CALL SPECPAK_DATA(  & 
               CCID_forC  & 
             , ID_SPECTRUM(ispec)  & 
             , MJD_SPECTRUM(ispec)  & 
             , TOBS_SPECTRUM(ispec)     &  ! May 7 2019
             , TEXPOSE_SPECTRUM(ispec)  & 
             , NLAMBIN_SPECTRUM(ispec)  & 
             , LAMMIN_SPECTRUM  &  ! array of LAMMIN in each bin
             , LAMMAX_SPECTRUM  &  ! array of LAMMAX in each bin
             , FLAM_SPECTRUM  & 
             , FLAMERR_SPECTRUM  & 
             , LENCCID )
 100     CONTINUE

! write out spectrum to TEXT table (ROOT not enabled)
    IF ( NSPECTRUM > 0 ) then
       CALL SPECPAK_FILL(CCID_forC, LENCCID)
    ENDIF

    RETURN
  END SUBROUTINE SPECPLOT


! ==============================================
    SUBROUTINE CHECK_DUPLICATE_MJDBAND(iep)

! Created Jun 13, 2017
! print & count error if current MJD & BAND is the same
! as previous epoch.


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IEP ! (I) epoch index

    REAL*8 MJD, MJD_LAST
    INTEGER IFILTOBS, IFILTOBS_LAST, i
    CHARACTER CFILT*2, cLIBID*40

! ---------------- BEGIN --------------

!       IF ( .NOT. ABORT_ON_DUPLMJD  ) RETURN
    IF ( IEP <= 1 ) RETURN

    MJD      = SNLC8_MJD(IEP)
    MJD_LAST = SNLC8_MJD(IEP-1)
    IFILTOBS = ISNLC_IFILT_OBS(IEP)
    IFILTOBS_LAST = ISNLC_IFILT_OBS(IEP-1)

    IF ( MJD      .NE. MJD_LAST      ) RETURN
    IF ( IFILTOBS .NE. IFILTOBS_LAST ) RETURN

    N_DUPLICATE_MJD = N_DUPLICATE_MJD + 1

! check if this MJD was already flagged as duplicate
    DO i = 1, min(200,NSTORE_DUPLICATE_MJD)
      IF ( MJD .EQ. DUPLICATE_MJDLIST(i) ) RETURN
    ENDDO

    NSTORE_DUPLICATE_MJD = NSTORE_DUPLICATE_MJD + 1

    if ( NSTORE_DUPLICATE_MJD < 200 ) then
       DUPLICATE_MJDLIST(NSTORE_DUPLICATE_MJD) = MJD
       CFILT = filtdef_string(ifiltobs:ifiltobs)
	 CLIBID = ''
	 IF ( LSIM_SNANA ) THEN
          write(CLIBID,600) SIM_LIBID
600	    format('LIBID=',I9)
       ENDIF
       write(6,666) SNLC_CCID(1:ISNLC_LENCCID), MJD, CFILT,  & 
              SNLC_FIELDLIST, cLIBID
666      format(/,T3,'*** DUPLICATE MJD WARNING *** CID=',A, 3x,  & 
              'MJD=', F10.4, 3x, 'BAND=', A, 2x,'FIELD=',A6,2x, A20 /)
       CALL FLUSH(6)
    endif

    RETURN
  END SUBROUTINE CHECK_DUPLICATE_MJDBAND


! =======================
! =======================

! call utilities to call table-interface C-functions

! ===========================================
    SUBROUTINE TABLE_OUTLIER(IDTABLE,IFLAG)
! 
! Created Mar 17 2021
! Create and fill table for OUTLIER fluxes w.r.t, SIM truth or LC model fit.
! Goal is to replace outlier option in sntable_dump so that making
! OUTLIER table is a 1-step process (snana.exe) instead of 2 steps
! (snana.exe + sntable_dump.py)
! 
! Aug 11 2021: check SBMAG cut


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML
    USE FILTCOM
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER IDTABLE  ! (I) table ID
    INTEGER IFLAG    ! (I) see IFLAG_XXX params

! local var

    INTEGER   LENNAME, LENFMT, NEWMJD, EPMIN, EPMAX, EP, BIT
    INTEGER   IFILTOBS, IFILT
    LOGICAL   DOFILL_TABLE, LCUT_NSIG, LCUT_FTRUE, LCUT_SBMAG
    REAL      NSIG, FTRUE, SBMAG, FLUXCAL_FIT, FLUXCAL_ERR_FIT
    CHARACTER NAME_forC*40, TEXTFMT*20, TEXTFMT_forC*20, FNAM*14

    REAL  NSIG_OUTLIER_FIT  ! function
    EXTERNAL  SNTABLE_CREATE, SNTABLE_FILL
! ---------------- BEGIN -----------------

    FNAM = 'TABLE_OUTLIER'

! for snana job, must be sim from either SNANA or FAKE
    IF ( ISJOB_SNANA .and. (.NOT. ISJOB_SIM)  ) THEN
       c1err = 'Cannot make OUTLIER table for ' //  & 
                 'real data using snana.exe.'
       c2err = 'Process SIM or FAKES; or run LC fit on real data'
       CALL MADABORT(FNAM, c1err, c2err)
    ENDIF

! FIT outliers work on real data or sim, so error check not needed.


    IF ( IFLAG .EQ. IFLAG_INI ) THEN

       OUTLIER_TABLE_NAME  = 'OUTLIER'
       LENNAME   = INDEX(OUTLIER_TABLE_NAME, ' ') - 1
       NAME_forC = OUTLIER_TABLE_NAME(1:LENNAME) // char(0)

       TEXTFMT  = TEXTFORMAT_TABLE(ITABLE_OUTLIER)
       LENFMT   = INDEX(TEXTFMT, ' ') - 1
       TEXTFMT_forC = TEXTFMT(1:LENFMT) // char(0)

       CALL SNTABLE_CREATE(IDTABLE, NAME_forC, TEXTFMT_forC,  & 
                LENNAME,LENFMT)  ! C fun

       CALL INIT_TABLE_OUTLIERVAR(IDTABLE, 'OUTLIER')

       CALL COUNT_OUTLIERS(0,0)

       RETURN
    ENDIF


! --------------------------
! fill table.

    if ( ISJOB_SNANA ) then
       DOFILL_TABLE = BTEST(CUTFLAG_SNANA,0)
    else
       DOFILL_TABLE = BTEST(CUTFLAG_SNANA,1)
    endif

    NEP_OUTLIER_PEREVT = 0
    IF ( .NOT. DOFILL_TABLE ) RETURN

!   terminate strings for C table-functions (snana.F90 routine)
    CALL TABLE_STRING_TERMINATION(+1)

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 201 EP = EPMIN, EPMAX
         IFILTOBS   = ISNLC_IFILT_OBS(ep)
         IFILT      = IFILTDEF_INVMAP_SURVEY(IFILTOBS)
         if ( ISNLC_SNRECON_USE(ep) == 0 ) goto 201

         if ( ISJOB_SNANA ) then
            NSIG   = abs(SIM_EPPULL(ep))
            FTRUE  = SIM_EPFLUXCAL(ep)
            SBMAG  = SNHOST_SBMAG(ifilt)
         else
#if defined(SNFIT)
            NSIG   = NSIG_OUTLIER_FIT(ep)
            FTRUE  = 1000.0  ! something guaranteed to pass cut
            SBMAG  = 20.0    ! idem
#endif
         endif

         if ( NSIG  < -0.01            ) goto 201
         if ( FTRUE < FTRUECUT_OUTLIER ) goto 201   ! Mar 24 2021

         LCUT_NSIG  = ( NSIG  >= NSIGCUT_OUTLIER  )
         LCUT_FTRUE = ( FTRUE >= FTRUECUT_OUTLIER )
         LCUT_SBMAG = ( SBMAG <= SBMAGCUT_OUTLIER ) ! Aug 2021

         CALL COUNT_OUTLIERS(1,ep)

         if ( LCUT_NSIG .and. LCUT_FTRUE .and. LCUT_SBMAG ) then
            NEP_OUTLIER_TOT    = NEP_OUTLIER_TOT    + 1
            NEP_OUTLIER_PEREVT = NEP_OUTLIER_PEREVT + 1
            CALL LOAD_OUTLIER_VAR(ep)
            CALL SNTABLE_FILL(IDTABLE) ! generic C function
            CALL COUNT_OUTLIERS(2,ep)
         endif
 201    CONTINUE
 200  CONTINUE

    CALL TABLE_STRING_TERMINATION(-1) ! remove termination

    RETURN
  END SUBROUTINE TABLE_OUTLIER


! ======================================
    SUBROUTINE COUNT_OUTLIERS(OPT,EP)


    USE SNDATCOM
    USE SNANAFIT
    USE FILTCOM
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER OPT  ! (I) 0=init, 1=count all, 2=count outliers, 3=print summary
    INTEGER EP   ! (I) epoch index

    INTEGER IFILTOBS, BIT, PHOTFLAG
    REAL    XNALL(2), XNOUTLIER(2), XNBULK(2)
    REAL    EFF_OUTLIER, EFF_ERR, EFF_BULK, FOM, XX, XXX
    LOGICAL ISBITSET

! ------------ BEGIN -----------

    if ( NSIGCUT_OUTLIER < 1.0 ) RETURN

    IF ( EP > 0 )  THEN
       IFILTOBS   = ISNLC_IFILT_OBS(ep)
       PHOTFLAG   = ISNLC_PHOTFLAG(ep)
    ENDIF

    IF ( OPT == 0 ) THEN
       ! one-time init
       NEP_OUTLIER_TOT = 0
       DO ifiltobs = 1, MXFILT_ALL
          NEP_OUTLIER_CHECK(ifiltobs) = 0  ! total epochs checked per band
          NEP_OUTLIER_FOUND(ifiltobs) = 0  ! Noutlier per band
       ENDDO

       DO BIT = 1, MXBIT_PHOTFLAG
          NEP_OUTLIER_VS_PHOTBIT(1,BIT) = 0
          NEP_OUTLIER_VS_PHOTBIT(2,BIT) = 0
          NEP_ALL_VS_PHOTBIT(1,BIT) = 0
          NEP_ALL_VS_PHOTBIT(2,BIT) = 0
       ENDDO

    ELSE IF ( OPT == 1 ) THEN
       ! passes cuts, but before checking if it's an outlier
       NEP_OUTLIER_CHECK(ifiltobs) = NEP_OUTLIER_CHECK(ifiltobs) + 1

       do bit = 1, MXBIT_PHOTFLAG
          ISBITSET = BTEST(PHOTFLAG,bit-1)
          NEP_ALL_VS_PHOTBIT(1,BIT) = NEP_ALL_VS_PHOTBIT(1,BIT) + 1
          if ( ISBITSET ) then
             NEP_ALL_VS_PHOTBIT(2,BIT) = NEP_ALL_VS_PHOTBIT(2,BIT) + 1
          endif
       ENDDO

    ELSE IF ( OPT == 2 ) THEN
       ! it's an OUTLIER
       NEP_OUTLIER_FOUND(ifiltobs) = NEP_OUTLIER_FOUND(ifiltobs) + 1

       do bit = 1, MXBIT_PHOTFLAG
          ISBITSET = BTEST(PHOTFLAG,bit-1)
          NEP_OUTLIER_VS_PHOTBIT(1,BIT) =  & 
            NEP_OUTLIER_VS_PHOTBIT(1,BIT) + 1
          if ( ISBITSET ) then
             NEP_OUTLIER_VS_PHOTBIT(2,BIT) =  & 
               NEP_OUTLIER_VS_PHOTBIT(2,BIT) + 1
          endif
       enddo

    ELSE IF ( OPT == 3 ) THEN
       ! print summary
       print*,' '
       write(6,10) NSIGCUT_OUTLIER

 10      format('  PHOTFLAG BIT PERFORMANCE for REJECTING ',  & 
                   F5.1, ' sigma OUTLIERS:', / )

       print*,'  BIT   EFF(OUTLIER)   1-EFF(BULK)   FoM'
       print*,'  ----------------------------------------- '
       do bit = 1, MXBIT_PHOTFLAG
          XNALL(1)     = float(NEP_ALL_VS_PHOTBIT(1,BIT))
          XNALL(2)     = float(NEP_ALL_VS_PHOTBIT(2,BIT))
          XNOUTLIER(1) = float(NEP_OUTLIER_VS_PHOTBIT(1,BIT))
          XNOUTLIER(2) = float(NEP_OUTLIER_VS_PHOTBIT(2,BIT))
          XNBULK(1)    = XNALL(1) - XNOUTLIER(1) ! not outliers
          XNBULK(2)    = XNALL(2) - XNOUTLIER(2) ! not outliers
          if ( XNALL(2) > 0.0 .and. XNOUTLIER(1) > 0.0 ) then
             EFF_OUTLIER = XNOUTLIER(2) / XNOUTLIER(1)
             XX      = XNOUTLIER(2) * (XNOUTLIER(1) - XNOUTLIER(2))
             XXX     = XNOUTLIER(1) * XNOUTLIER(1) * XNOUTLIER(1)
             EFF_ERR     = sqrt(XX/XXX)
             EFF_BULK    = 1.0 - XNBULK(2) / XNBULK(1)
             FOM         = EFF_OUTLIER * EFF_BULK
             write(6,30) bit-1, EFF_OUTLIER, EFF_ERR, EFF_BULK, FoM
 30            format(T5,I2, 2x, F6.3,' +_', F6.3, 3x, F6.3, 3x, F6.3 )
          endif
       enddo
       call flush(6)
    ENDIF

    RETURN
  END SUBROUTINE COUNT_OUTLIERS

! =============================
    SUBROUTINE INIT_TABLE_OUTLIERVAR(ID, BLOCK)


    USE SNDATCOM
    USE SNANAFIT
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER   ID          ! (I) table ID
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)

! local var

    INTEGER  & 
         IFILT, IFILT_OBS, ITEXT  & 
        ,LENBLOCK, LENLIST, LENNAME, LENV, ipar, ivar

    LOGICAL LTMP, ADDCOL_SETPKMJD
    CHARACTER varlist*100, CTMP*60, CBLOCK*40
    CHARACTER FNAM*22, CFILT*2
    EXTERNAL  & 
          SNTABLE_ADDCOL  & 
         ,SNTABLE_ADDCOL_int  & 
         ,SNTABLE_ADDCOL_flt  & 
         ,SNTABLE_ADDCOL_dbl  & 
         ,SNTABLE_ADDCOL_str

! ------------ BEGIN ----------

    FNAM = 'INIT_TABLE_OUTLIERVAR'

    write(6,10) BLOCK, ID
 10   format(T6,'Create BLOCK = ',A,'  for TABLE ID = ', I5 )
    call flush(6)

    LENBLOCK     = INDEX(BLOCK//' ',' ') - 1
    CBLOCK       = BLOCK(1:LENBLOCK) // char(0)
    LENLIST      = LEN(VARLIST)

!  for CCID, use special variable SNLC_CCID_forC that has char(0)
!  termination so that parsing CCID with C code is easier.

    VARLIST = 'CCID:C*20' // char(0)
    CALL SNTABLE_ADDCOL_str(ID, CBLOCK, SNLC_CCID, VARLIST,1,  & 
                         LENBLOCK, 20 )

! ROWNUM is easy identifier for CID-MJD-BAND-DETNUM;
! e.g., pass ROWNUM to script that fetches stamps.
    VARLIST = 'ROWNUM:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NEP_OUTLIER_TOT, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'FIELD:C*12' // char(0)
    CALL SNTABLE_ADDCOL_str(ID, CBLOCK, OUTLIER_FIELD, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'MJD:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, OUTLIER_MJD, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'BAND:C*4' // char(0)
    CALL SNTABLE_ADDCOL_str(ID, CBLOCK, OUTLIER_BAND, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'IFILTOBS:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, OUTLIER_IFILTOBS, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'DETNUM:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, OUTLIER_DETNUM, VARLIST,1,  & 
                         LENBLOCK, 20 )

! --
    VARLIST = 'IXPIX:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, OUTLIER_IXPIX, VARLIST,1,  & 
                         LENBLOCK, 20 )
    VARLIST = 'IYPIX:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, OUTLIER_IYPIX, VARLIST,1,  & 
                         LENBLOCK, 20 )
    VARLIST = 'AREAFRAC:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_AREAFRAC, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'PHOTFLAG:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, OUTLIER_PHOTFLAG, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'FLUXCAL_DATA:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_DATA,  & 
                 VARLIST, 1, LENBLOCK, 20 )

    VARLIST = 'FLUXCAL_ERR_DATA:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_ERR_DATA,  & 
              VARLIST, 1, LENBLOCK, 20 )

    VARLIST = 'FLUXCAL_ERR_CALC:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK,OUTLIER_FLUXCAL_ERR_CALC,  & 
                  VARLIST, 1, LENBLOCK, 20 )

    if ( ISJOB_SNANA ) then
       ! outlier w.r.t. true sim/fake flux (SIM_MAGOBS column in data file)
       VARLIST = 'FLUXCAL_TRUE:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_TRUE,  & 
              VARLIST, 1, LENBLOCK, 20 )
    else
       ! outlier w.r.t. LC fit
       VARLIST = 'FLUXCAL_FIT:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_FIT,  & 
              VARLIST, 1, LENBLOCK, 20 )

       VARLIST = 'FLUXCAL_ERR_FIT:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_ERR_FIT,  & 
               VARLIST, 1, LENBLOCK, 20 )

       VARLIST = 'FLUXCAL_ERRTOT_FIT:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_FLUXCAL_ERRTOT_FIT,  & 
               VARLIST, 1, LENBLOCK, 20 )
    endif

    VARLIST = 'NSIG:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_NSIG, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'LOGSNR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_LOGSNR, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'SBMAG:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_SBMAG, VARLIST,1,  & 
                         LENBLOCK, 20 )

! - - - -
    VARLIST = 'ZP:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_ZP, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'PSF:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_PSF, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'SKYSIG:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_SKYSIG, VARLIST,1,  & 
                         LENBLOCK, 20 )

#if defined(SNFIT)
    VARLIST = 'z:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_z, VARLIST,1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'TREST:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, OUTLIER_TREST, VARLIST,1,  & 
                         LENBLOCK, 20 )

#endif
    RETURN
  END SUBROUTINE INIT_TABLE_OUTLIERVAR


! =============================
    SUBROUTINE LOAD_OUTLIER_VAR(ep)

! Created Mar 17 2021
! Load OUTLIER_XXX variables for epoch 'ep'.
! 
    USE SNDATCOM
    USE SNANAFIT
    USE FILTCOM
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER  ep  ! (I) epoch index to load

    INTEGER IFILT_OBS, IFILT
    CHARACTER BAND*4
    REAL SNR

! ---------- BEGIN ---------

    IFILT_OBS   = ISNLC_IFILT_OBS(ep)
    IFILT       = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
    BAND        = filtdef_string(ifilt_obs:ifilt_obs)

    OUTLIER_MJD          = SNLC8_MJD(ep)
    OUTLIER_IFILTOBS     = IFILT_OBS
    OUTLIER_BAND         = BAND
    OUTLIER_FIELD        = SNLC_FIELD(ep)
    OUTLIER_DETNUM       = ISNLC_DETNUM(ep)
    OUTLIER_PHOTFLAG     = ISNLC_PHOTFLAG(ep)
    OUTLIER_FLUXCAL_DATA      = SNLC_FLUXCAL(ep)
    OUTLIER_FLUXCAL_ERR_DATA  = SNLC_FLUXCAL_ERRTOT(ep)
    OUTLIER_FLUXCAL_ERR_CALC  = SNLC_FLUXCAL_ERRCALC(ep)
    OUTLIER_NSIG         = abs(SIM_EPPULL(ep))
    OUTLIER_ZP           = SNLC_ZEROPT(ep)
    OUTLIER_PSF          = SNLC_PSF_SIG1(ep)
    OUTLIER_SKYSIG       = SNLC_SKYSIG(ep)

    SNR  = OUTLIER_FLUXCAL_DATA/OUTLIER_FLUXCAL_ERR_DATA
    OUTLIER_LOGSNR       = LOG10(max(0.1,SNR))

    OUTLIER_SBMAG        = SNHOST_SBMAG(ifilt)

    OUTLIER_IXPIX     = INT(SNLC_XPIX(ep))
    OUTLIER_IYPIX     = INT(SNLC_YPIX(ep))
    OUTLIER_AREAFRAC  = SNLC_AREAFRAC(ep)

    OUTLIER_z         = SNLC_REDSHIFT
    OUTLIER_TREST     = SNLC_TREST(ep)

    IF ( ISJOB_SNANA ) THEN
       ! true flux from SIM_MAGOBS
       OUTLIER_FLUXCAL_TRUE      = SIM_EPFLUXCAL(ep)
    ELSE
       ! get model fit flux
       OUTLIER_FLUXCAL_FIT        = SNLC_FLUXCAL_FIT(ep)
       OUTLIER_FLUXCAL_ERR_FIT    = SNLC_FLUXCAL_ERR_FIT(ep)
       OUTLIER_FLUXCAL_ERRTOT_FIT = SNLC_FLUXCAL_ERRTOT_FIT(ep)
    ENDIF

    IF ( USE_TABLEFILE_ROOT ) THEN
       OUTLIER_BAND = OUTLIER_BAND(1:1) // char(0)
    ENDIF

    RETURN
  END SUBROUTINE LOAD_OUTLIER_VAR


! ===========================================
#if defined(SNFIT)
    REAL FUNCTION NSIG_OUTLIER_FIT(ep)

! Created Mar 17 2021
! Function returns Nsig = | FLUXCAL_DATA / ERRTOT |
! where ERRTOT^2 = FLUXCAL_ERR_DATA^2 + FLUXCAL_ERR_FIT^2
! Initial use is for OUTLIER table.
! 

    USE SNDATCOM
    USE SNANAFIT
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER  ep  ! (I) epoch index

    REAL  FLUXCAL_FIT, FLUXCAL_ERR_FIT
    REAL  FLUXCAL_DATA, FLUXCAL_ERR_DATA, SQERR, ERRTOT, NSIG
    LOGICAL REJECT
! --------- BEGIN ---------

    NSIG = 0.0

    CALL GET_FITFLUX(ep, FLUXCAL_FIT, FLUXCAL_ERR_FIT,REJECT)

    FLUXCAL_DATA     = SNLC_FLUXCAL(ep)
    FLUXCAL_ERR_DATA = SNLC_FLUXCAL_ERRTOT(ep)

    IF ( FLUXCAL_ERR_FIT > 0.00001 ) THEN
       SQERR = (FLUXCAL_ERR_FIT  * FLUXCAL_ERR_FIT) +  & 
                 (FLUXCAL_ERR_DATA * FLUXCAL_ERR_DATA)

       ERRTOT = SQRT(SQERR)
       NSIG   = abs(FLUXCAL_DATA-FLUXCAL_FIT) / ERRTOT
    ELSE
       ERRTOT = FLUXCAL_ERR_DATA
       NSIG   = -9.0  ! epoch not used
    ENDIF

!      print*,' xxx ep, DIF/ERR = ',
!     &    ep, (FLUXCAL_DATA-FLUXCAL_FIT), '/', ERRTOT, NSIG

! load fit flux info to global so that it can be used elsewhere
    SNLC_FLUXCAL_FIT(ep)        = FLUXCAL_FIT
    SNLC_FLUXCAL_ERR_FIT(ep)    = FLUXCAL_ERR_FIT
    SNLC_FLUXCAL_ERRTOT_FIT(ep) = ERRTOT

! set function value
    NSIG_OUTLIER_FIT = NSIG

    RETURN
  END FUNCTION NSIG_OUTLIER_FIT
#endif

! ==================================
    SUBROUTINE PRINT_OUTLIER_SUMMARY()

! Created Mar 17 2021
! print OUTLIER fraction per band to stdout.
! If TEXTFILE_PREFIX is set, prepend top of file with summary
! using clumsy method of writing a separate header file and
! then using system cat command.


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML
    USE FILTCOM
    USE OUTLIERCOM

    IMPLICIT NONE

    INTEGER IFILT, IFILTOBS, N0, N1, LEN_WHAT
    INTEGER LEN1, LEN2, LEN3, LEN_HEAD
    REAL FRAC
    CHARACTER cfilt*2, TEXT_WHAT*40, CMD_CAT*200, CMD_MOVE*200
    CHARACTER TEXTFILE*200, TMP_HEADER_FILE*200, LINE*60, LINE2*60
    LOGICAL   PREPEND_TEXTFILE

! ---------- BEGIN -----------

    call count_outliers(3,0)

    PREPEND_TEXTFILE = ( TEXTFILE_PREFIX .NE. ' ' )
    IF ( PREPEND_TEXTFILE ) THEN
! construct name of TEXT file to prepend outlier summary.
! Beware that TEXTFILE is fragile.
       LEN1 = INDEX(TEXTFILE_PREFIX,   ' ' ) - 1
       LEN2 = INDEX(OUTLIER_TABLE_NAME,' ' ) - 1
       TEXTFILE = TEXTFILE_PREFIX(1:LEN1) // '.' //  & 
             OUTLIER_TABLE_NAME(1:LEN2) // '.TEXT'
       LEN3 = INDEX(TEXTFILE,' ') - 1
       TMP_HEADER_FILE =  & 
             'TMP_HEAD_' // TEXTFILE_PREFIX(1:LEN1) // '.TEXT'
       LEN_HEAD = INDEX(TMP_HEADER_FILE,' ') - 1
       OPEN(UNIT=LUNTMP, FILE=TMP_HEADER_FILE, STATUS='UNKNOWN')
    ENDIF

! - - - - - - -

    IF ( ISJOB_SNANA ) THEN
       TEXT_WHAT  = '|F_data-F_true(SIM)|/Sigma_F'
    ELSE
       TEXT_WHAT  = '|F_data-F_fit|/Sigma_F'
    ENDIF
    LEN_WHAT = INDEX(TEXT_WHAT,' ') - 1

    LINE = '    FLUX-OUTLIER SUMMARY per BAND for '
    write(LINE2,10) TEXT_WHAT(1:LEN_WHAT), NSIGCUT_OUTLIER
 10   format(T5, A, ' > ', F5.2, ' : ' )

    print*, ' '
    print*, LINE
    print*, LINE2
    call flush(6)

    if (PREPEND_TEXTFILE) then
       WRITE(LUNTMP,70) '=========================================='
       WRITE(LUNTMP,70) LINE
       WRITE(LUNTMP,70) LINE2
    endif

 70   format('# ', A)
! - - - - - -

    DO ifilt  = 1, NFILTDEF_SURVEY
       ifiltobs  = IFILTDEF_MAP_SURVEY(ifilt)
       cfilt     = filtdef_string(ifiltobs:ifiltobs)

       N0 = NEP_OUTLIER_CHECK(ifiltobs)
       N1 = NEP_OUTLIER_FOUND(ifiltobs)
       if ( N0 > 0 ) then
          FRAC = FLOAT(N1) / FLOAT(N0)
       else
          FRAC = 0.0
       endif

       write(LINE,20) cfilt, FRAC, N1, N0
 20      format(T8,A1, '-Outlier frac = ', F8.5, 3x,  & 
            '(', I6, '/', I6, ')'  )

       print*, LINE
       if (PREPEND_TEXTFILE) WRITE(LUNTMP,70) LINE

    ENDDO

    if(PREPEND_TEXTFILE) then
       WRITE(LUNTMP,70) ' '
       CLOSE(UNIT = LUNTMP)
! finally, use cat to glue header comments and table.
       CMD_CAT = 'cat ' // TEXTFILE(1:LEN3) //  & 
               ' >> ' // TMP_HEADER_FILE
       CMD_MOVE = 'mv ' // TMP_HEADER_FILE(1:LEN_HEAD) //  ' ' //  & 
                 TEXTFILE(1:LEN3)
       CALL SYSTEM(CMD_CAT // ' ; ' // CMD_MOVE)
    endif

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE PRINT_OUTLIER_SUMMARY

! ===========================================
    SUBROUTINE TABLE_SNANA(IDTABLE,IFLAG)
! 
! Apr 26 2017: check PRESCALE_TABLE
! 
! ------------------------------

    USE SNPAR
    USE SNCUTS
    USE CTRLCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER IDTABLE  ! (I) table ID
    INTEGER IFLAG    ! (I) see IFLAG_XXX params

! local var

    LOGICAL   DO_FILL, LPS
    REAL*8    PS
    INTEGER   LENNAME, LENFMT
    CHARACTER NAME*40, TEXTFMT*20, TEXTFMT_forC*20, FNAM*12
    EXTERNAL SNTABLE_CREATE, SNTABLE_FILL

    LOGICAL REJECT_PRESCALE  ! function

! -------------------- BEGIN -----------------
    FNAM = 'TABLE_SNANA'

    IF ( IFLAG .EQ. IFLAG_INI ) THEN

       NAME     = 'SNANA' // char(0)
       LENNAME  = INDEX(NAME,    ' ') - 1

       TEXTFMT  = TEXTFORMAT_TABLE(ITABLE_SNANA)
       LENFMT   = INDEX(TEXTFMT, ' ') - 1
       TEXTFMT_forC = TEXTFMT(1:LENFMT) // char(0)

       CALL SNTABLE_CREATE(IDTABLE,NAME,TEXTFMT_forC,  & 
                LENNAME,LENFMT)  ! C fun

! make table of scalars
       CALL INIT_TABLE_SNANAVAR(IDTABLE, 'SNANA', 1)

! check option to store variables for each obs (Mar 2015)
       IF ( OPT_TABLE(ITABLE_SNANA) .EQ. 2 ) THEN
         CALL INIT_TABLE_SNOBSVAR(IDTABLE, 'SNOBS')
       ENDIF

       IF ( OPT_TABLE(ITABLE_SNANA) .EQ. 4 ) THEN ! Jan 2019
         CALL INIT_TABLE_SIM_MAGOBS(IDTABLE, 'MODELMAG')
       ENDIF

! for sims ...

       CALL INIT_TABLE_SIMVAR(IDTABLE, 'SIM')
       CALL FLUSH(6)
       RETURN
    ENDIF

! --------------------------
! fill table.

    DO_FILL = .FALSE.

    IF ( CUTMASK_SNANA_TABLE == 63 )  & 
                DO_FILL = .TRUE.   ! default --> include everything

! check subset options
    IF (CUTFLAG_SNANA == 0 .and. CUTMASK_SNANA_TABLE == 0 )  & 
                DO_FILL = .TRUE.
    IF( IAND(CUTFLAG_SNANA,CUTMASK_SNANA_TABLE) > 0 )  & 
                DO_FILL = .TRUE.


    IF ( DO_FILL ) THEN

!       terminate strings for C table-functions (snana.F90 routine)
      CALL TABLE_STRING_TERMINATION(+1)

! check option to store MODEL-MAG on epoch grid (Feb 2019)
      IF ( OPT_TABLE(ITABLE_SNANA) .EQ. 4 ) THEN
         CALL LOAD_TABLE_SIM_MAGOBS()
      ENDIF

      PS  = DBLE(PRESCALE_TABLE(ITABLE_SNANA))
      LPS = REJECT_PRESCALE(NCALL_SNANA_DRIVER,PS)
      IF ( .not. LPS ) THEN
        CALL SNTABLE_FILL(IDTABLE)  ! generic C function
      ENDIF

      CALL TABLE_STRING_TERMINATION(-1)  ! remove termination
    ENDIF

    RETURN
  END SUBROUTINE TABLE_SNANA

! =====================================
    SUBROUTINE INIT_TABLE_SNANAVAR(ID,BLOCK,IFLAG)
! 
! Created Feb 1, 2013
! 
! Initialize variables for snana-analysis table.
! Call generic function SNTABLE_ADDCOL() that loads a table
! for each selected format: CERNLIB, ROOT, HDF5 ...
! 
! IFLAG = 1 => regular SNANA variables
! IFLAG = 2 => use _FIT variables instead (if they exist) since
!              post-fit variables may be different than the
!              original SNANA variables (e.g., SNRMAX3 can change
!              if a filter is excluded from a fit.
! 
! 
! May 31 2016: add zCMB[ERR] to text file and remove legacy "z,zerr"
!              from text file [since z=zHELIO]
! 
! Jan  8 2018: write zHEL[ERR] to text file
! Mar  9 2018: add PHOTPROB info, MIN and Nepochs
! Mar 10 2019: add HOST_OBJID with BLOCK_LL for long long int
! May 23 2019: bug fix: book HOST_RA[DEC] as double instead of float.
! Feb 20 2020: write MWEBV to TEXT table.
! May 14 2020: add zFLAG
! Mar 14 2021: always include PKMJDINI in TEXT output
! Mar 19 2021: add NEP_OUTLIER if OUTLIER table is defined
! Aug 13 2021: include HOST_SBMAG, and HOST_SB -> HOST_SBFLUXCAL
! Sep 09 2023: update to write 2nd HOST match.
! Jul 23 2024: add FITPROB_ITER1, FITCHI2RED_INI[2]
! Feb 05 2025: add LENSDMU[ERR]
! -------------------------------------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE FILTCOM
    USE SNLCINP_NML
    USE USRTAGCM
    USE PRIVCOM
    USE PKMJDCOM
    USE REQEPCOM
    USE OUTLIERCOM
    USE TRUECHI2COM

    IMPLICIT NONE

    INTEGER   ID          ! (I) table ID
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)
    INTEGER   IFLAG       ! (I) SNANA or _FIT variables

! local var

    INTEGER  & 
         IFILT, IFILT_OBS, ITEXT, ITEXT2, ITEXT_NO, ITEXT_YES, IGAL  & 
        ,LENBLOCK, LENLIST, LENNAME, LENV, ipar, ivar

    LOGICAL LTMP, ADDCOL_SETPKMJD
    CHARACTER varlist*1000, CTMP*60, CBLOCK*40
    CHARACTER FNAM*20, VARNAME*40, CFILT*2
    EXTERNAL  & 
         SNTABLE_ADDCOL  & 
        ,SNTABLE_ADDCOL_int  & 
        ,SNTABLE_ADDCOL_flt  & 
        ,SNTABLE_ADDCOL_dbl  & 
        ,SNTABLE_ADDCOL_str

! ---------------- BEGIN ------------

    FNAM = 'INIT_TABLE_SNANAVAR'
    ITEXT_NO  = 0   ! --> for ROOT table only
    ITEXT_YES = 1   ! --> for TEXT and ROOT tables

    write(6,10) BLOCK, ID
 10   format(T6,'Create BLOCK = ',A,'  for TABLE ID = ', I5 )
    call flush(6)

    LENBLOCK     = INDEX(BLOCK//' ',' ') - 1
    CBLOCK       = BLOCK(1:LENBLOCK) // char(0)
    LENLIST      = LEN(VARLIST)

! check filter list for filter-dependent columns (Feb 2017)
    ADDCOL_FILTERS = SURVEY_FILTERS  ! default
    IF ( NFILT_REMAP_TABLE > 0 ) THEN
       ADDCOL_FILTERS = FILTLIST_REMAP_TABLE
    ENDIF

    LTMP = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_SAVEPAR)
    ADDCOL_SETPKMJD = ( LTMP .and. OPT_SETPKMJD > 0 )

! - - - - - - ID variables - - - - - -

    VARLIST = 'CID:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SNLC_CID, VARLIST,0, LENBLOCK, 20 )

!  for CCID, use special variable SNLC_CCID_forC that has char(0)
!  termination so that parsing CCID with C code is easier.

    VARLIST = 'CCID:C*20' // char(0)
    CALL SNTABLE_ADDCOL_str(ID, CBLOCK, SNLC_CCID, VARLIST,1, LENBLOCK, 20 )

    if ( WRTABLEFILE_IAUC ) then
      VARLIST = 'IAUC:C*20' // char(0)
      CALL SNTABLE_ADDCOL_str(ID, CBLOCK, SNLC_NAME_IAUC, VARLIST,1, LENBLOCK, 20 )
    endif

! Aug 2022: check optional user-column; e.g., to tag magshift value
    if ( SNTABLE_APPEND_VARNAME .NE. '' ) then
       LENV = INDEX(SNTABLE_APPEND_VARNAME,' ') - 1
       VARLIST = SNTABLE_APPEND_VARNAME(1:LENV) // ':F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK,SNTABLE_APPEND_VALUE, VARLIST,1,  LENBLOCK, 20 )
    endif
! to do: when combined low-z is ready change IDSURVEY -> IDSUBSURVEY,
! but keep same name 'IDSURVEY'
    VARLIST = 'IDSURVEY:I' //  char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, IDSUBSURVEY, VARLIST,1, LENBLOCK, 20 )

    VARLIST = 'TYPE:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_TYPE, VARLIST,1, LENBLOCK, 20 )

    VARLIST = 'NFIELD_OVP:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFIELD_OVP, VARLIST, 0, LENBLOCK, 20 )

    VARLIST = 'FIELD:C*20' // char(0)
    CALL SNTABLE_ADDCOL_str(ID, CBLOCK, SNLC_FIELDLIST, VARLIST,1, LENBLOCK, 20 )

! Feb 15 2018: remove condition on DETNUM
    ITEXT = ITEXT_NO
    IF(NDETNUM_LIST>0) ITEXT = ITEXT_YES
    VARLIST = 'DETNUM:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_DETNUM(1), VARLIST,ITEXT,  LENBLOCK, 20 )

! -------------------------------------------
! CUTFLAG_SNANA is a mask as follows:
! 0    -> failed SNANA cuts
! bit0 -> passed SNANA cuts
! bit1 -> passed FIT cuts (i.e, valid fit + re-applied snana cuts)

    VARLIST = 'CUTFLAG_SNANA:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, CUTFLAG_SNANA, VARLIST,1,  LENBLOCK, 20)

! -------------------------------------------
! ERRFLAG_FIT = -9 -> no fit, just SNANA job
! ERRFLAG_FIT =  0 -> no fitting errors.
! ERRFLAG_FIT >  0 -> fit failed, see ERRFLAG codes in fitting program.

    ITEXT = ITEXT_NO
#if defined(SNFIT)
    ITEXT=ITEXT_YES
#endif

    VARLIST = 'ERRFLAG_FIT:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ERRFLAG_FIT, VARLIST, ITEXT,   LENBLOCK, 20 )

! --------------------------------------------
    IF ( N_USERTAGS .GT. 0) then
      VARLIST = 'USERTAG:I'  // char(0)
      CALL SNTABLE_ADDCOL_int(ID, CBLOCK, USERTAG, VARLIST,1,  LENBLOCK, 20)
    endif

    VARLIST = 'RA:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_RA, VARLIST,0,  LENBLOCK, 20 )
    VARLIST = 'DEC:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_DEC, VARLIST,0,  LENBLOCK, 20 )

    VARLIST = 'AREAFRAC_AVG:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_AREAFRAC_AVG, VARLIST,0,  LENBLOCK, 20 )

! -- NEPOCH ( Aug 2013)

    VARLIST = 'NEPOCH:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NEPOCH_STORE,  & 
                     VARLIST,0,  LENBLOCK, 20 )

    VARLIST = 'NEPOCH_BADPHOT:I' // char(0)      ! added Jun 2016
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NEPOCH_BADPHOT, VARLIST, 0,  & 
                         LENBLOCK, 20 )


! Mar 2018: add info about PHOTPROB
    VARLIST = 'NEPOCH_PHOTPROB:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NEPOCH_PHOTPROB,  & 
                  VARLIST,0,  LENBLOCK, 20 )

    VARLIST = 'PHOTPROB_MIN:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_PHOTPROB_MIN,  & 
                      VARLIST,0,   LENBLOCK, 20 )

! Sep 2018: add optonal info about detections
! Nov 2020: write to TEXT file
    IF ( PHOTFLAG_DETECT > 0 ) THEN
       VARLIST = 'NOBS_DETECT:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NOBS_DETECT, VARLIST,1, LENBLOCK, 20 )
       VARLIST = 'TLIVE_DETECT:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TLIVE_DETECT, VARLIST,1,   LENBLOCK, 20 )
    ENDIF


    IF ( OPT_TABLE(ITABLE_OUTLIER) > 0 ) THEN
       VARLIST = 'NEPOCH_OUTLIER:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NEP_OUTLIER_PEREVT,  & 
                  VARLIST,1,   LENBLOCK, 20 )
    ENDIF
! - - - - -  redshift - - - - -

    VARLIST = 'zHEL:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZHELIO, VARLIST, 1, LENBLOCK, 20 )
    VARLIST = 'zHELERR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZHELIO_ERR, VARLIST, 1, LENBLOCK, 20 )

    VARLIST = 'zCMB:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZCMB, VARLIST,1,  LENBLOCK, 20 )
    VARLIST = 'zCMBERR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZCMB_ERR, VARLIST,1,  LENBLOCK, 20 )

    VARLIST = 'zHD:F' // char(0)  ! Hubble diagram redshift
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZHD, VARLIST,1, LENBLOCK, 20 )
    VARLIST = 'zHDERR:F' // char(0)  ! Hubble diagram redshift
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_ZHD_ERR, VARLIST,1, LENBLOCK, 20 )

    VARLIST = 'zFLAG:I' // char(0)  ! z-quality flag (May 2020)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_zFLAG, VARLIST, 0, LENBLOCK, 20 )

! xxxxxxxx mark delete Nov 30 2025xxxxxxx
! z = zCMB for spec z, but for photo-z fit z & zERR are the original
! data redshifts to allow plotting zHD-z resids
!    VARLIST = 'z:F,zERR:F'  // char(0)
!    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_REDSHIFT, VARLIST, 0,  & 
!                         LENBLOCK, 20)
! xxxxxxxxxx end mark xxxxxxx

! tack on VPEC
    VARLIST = 'VPEC:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_VPEC, VARLIST, 1, LENBLOCK, 20 )

    VARLIST = 'VPECERR:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_VPEC_ERR, VARLIST, 1, LENBLOCK, 20 )

! tack on LENSDMU (Feb 2025)
    VARLIST = 'LENSDMU:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_LENSDMU, VARLIST, 1,  LENBLOCK, 20 )

    VARLIST = 'LENSDMUERR:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_LENSDMU_ERR, VARLIST, 1,  LENBLOCK, 20 )

! - - - - -  Galactic extinction - - - -  -
! Feb 6 2020: write MWEBV to TEXT (but not the error)

    VARLIST = 'MWEBV:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_MWEBV, VARLIST, 1,  & 
                         LENBLOCK, 20 )

    VARLIST = 'MWEBVERR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_MWEBV, VARLIST, 0,  & 
                         LENBLOCK, 20 )

! - - - - - - HOST GALAXY - - - -
    ITEXT = ITEXT_NO
    if ( WRTABLEFILE_HOST_TEXT .or. WRTABLEFILE_HOST2_TEXT ) ITEXT=1

    ITEXT2 = ITEXT_NO
    if ( WRTABLEFILE_HOST2_TEXT ) ITEXT2=ITEXT_YES  ! for 2nd host only

    VARLIST =  'HOST_NMATCH:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SNHOST_NMATCH,  & 
                    VARLIST,ITEXT,   LENBLOCK, 20 )
    VARLIST =  'HOST_NMATCH2:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SNHOST_NMATCH2,  & 
                    VARLIST,ITEXT,    LENBLOCK, 20 )

    IGAL = 1
    CALL INIT_TABLE_HOSTVAR(ID, IGAL, ITEXT, BLOCK)

    IGAL = 2
    CALL INIT_TABLE_HOSTVAR(ID, IGAL, ITEXT2, BLOCK)


! Apr 2022
! store host surface brightness if SB exists or host-mag exists.
! The host-mag check is a safety feature in case SB is undefined.
! E.g., DES event 1248110 has SBFLUX=-999 even though host mags are
! all in the 21-22 range. Likely a subtle problem measuring SBFLUX
! on the template image.

    IF ( EXIST_SNHOST_SB .OR. EXIST_SNHOST_MAGOBS) THEN
!        HOST_SB -> HOST_SBFLUXCAL
      CALL TABLE_VARLIST_FILTERS('HOST_SBFLUXCAL','F', ADDCOL_FILTERS,  & 
          VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_SNHOST_SBFLUXCAL(1),  & 
          VARLIST(1:LENLIST)//char(0), itext, LENBLOCK, LENLIST)

      CALL TABLE_VARLIST_FILTERS('HOST_SBMAG','F', ADDCOL_FILTERS,  & 
          VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_SNHOST_SBMAG(1),  & 
          VARLIST(1:LENLIST)//char(0), itext, LENBLOCK, LENLIST)
    ENDIF


! - - - - -  pre-fit quantities to show for all events  - -  -

! PKMJDINI (e.g., from max-flux clump method)

    VARLIST = 'PKMJDINI:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_SEARCH_PEAKMJD,  & 
                 VARLIST, ITEXT_YES, LENBLOCK, 20 ) ! needed by SNN classifier

    VARLIST = 'FITPROB_ITER1:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, FITPROB_ITER1,  & 
                    VARLIST, ITEXT_NO,   LENBLOCK, 20 )

    VARLIST = 'FITCHI2RED_INI:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, FITCHI2RED_INI,  & 
                    VARLIST, ITEXT_NO,   LENBLOCK, 20 )

    VARLIST = 'FITCHI2RED_INI2:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, FITCHI2RED_INI2,  & 
                    VARLIST, ITEXT_NO,   LENBLOCK, 20 )

! - - - - time above threshold for REQUIRE_EPOCHS - - - - - -

    IF( NFILT_REQEP > 0 ) THEN
      VARLIST = 'NDAYS_ABOVE_SNRMIN:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, NDAYS_ABOVE_SNRMIN_REQEP,  & 
                          VARLIST, ITEXT_YES,       LENBLOCK, 40 )
    ENDIF

! - - - - - max flux per band (May 2016) - - - - - - - - -

      CALL TABLE_VARLIST_FILTERS('FLUXCALMAX','F', ADDCOL_FILTERS,  & 
          VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_FLUXCALMAX(1),  & 
          VARLIST(1:LENLIST)//char(0), ITEXT_NO, LENBLOCK, LENLIST)

! ... and the uncertainty (Nov 2022)

      CALL TABLE_VARLIST_FILTERS('FLUXCALMAX_ERR','F',ADDCOL_FILTERS,  & 
         VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_FLUXCALMAX_ERR(1),  & 
          VARLIST(1:LENLIST)//char(0), ITEXT_NO, LENBLOCK, LENLIST)

! - - - - - - optional fit-params from PKMJD fit (Oct 20 2014) - - -

    IF ( ADDCOL_SETPKMJD ) THEN
       DO ipar = IPAR_T0, IPAR_A2
          VARNAME = 'FIT_' // PKPARNAME(ipar)
          LENV    = INDEX(VARNAME,' ') - 1
          CALL TABLE_VARLIST_FILTERS(VARNAME(1:LENV),'F',  & 
                 SURVEY_FILTERS, VARLIST, LENLIST)

          CALL SNTABLE_ADDCOL_flt(ID, CBLOCK,  & 
              SNLC_SNANAFIT_PEAKMJD_FITPAR(1,ipar),  & 
              VARLIST(1:LENLIST)//char(0), ITEXT_YES, LENBLOCK, LENLIST)
       ENDDO

! ----------------------------
!  PKMJD per filter

       DO 222 ifilt  = 1, NFILTDEF_SURVEY
          ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
          cfilt     = filtdef_string(ifilt_obs:ifilt_obs)

          IF( IFILTOBS_REPLACE(IFILT_OBS) .NE. IFILT_OBS) GOTO 222

          VARNAME   = 'PKMJD_FIT_' // cfilt(1:1) // ':F'
          CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, PKMJD_FIT(ifilt_obs),  & 
                 VARNAME(1:20)//char(0), ITEXT_YES, LENBLOCK, 20 )

          VARNAME   = 'PKMJD_ERR_' // cfilt(1:1) //  ':F'
          CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, PKMJD_ERR(ifilt_obs),  & 
                 VARNAME(1:20)//char(0), ITEXT_YES, LENBLOCK, 20 )

 222     CONTINUE

! --------------------------
!      fit chi2 per filter
        CALL TABLE_VARLIST_FILTERS('FIT_CHI2','F', ADDCOL_FILTERS,  & 
            VARLIST, LENLIST )
        CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_CHI2_FITPKMJD(1),  & 
            VARLIST(1:LENLIST)//char(0), ITEXT_YES, LENBLOCK, LENLIST)

! --------------------------
!      fit NDOF per filter
        CALL TABLE_VARLIST_FILTERS('FIT_NDOF','I', SURVEY_FILTERS,  & 
            VARLIST, LENLIST)
        CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NDOF_FITPKMJD(1),  & 
            VARLIST(1:LENLIST)//char(0), ITEXT_YES, LENBLOCK, LENLIST)

    ENDIF  ! end of SETPKMJD if-block

! - - - - - -
    IF ( LPROB_TRUEFLUX .and. ISJOB_SIM ) THEN

! global PROB_TRUEFLUX for all bands combined
      VARNAME   = 'PROB_TRUEFLUX:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, PROB_TRUEFLUX(0),  & 
              VARNAME(1:20), ITEXT_YES, LENBLOCK, 20)

! ... and for each individual band
      CALL TABLE_VARLIST_FILTERS('PROB_TRUEFLUX','F', ADDCOL_FILTERS,  & 
           VARLIST, LENLIST )
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_PROB_TRUEFLUX(1),  & 
          VARLIST(1:LENLIST)//char(0), ITEXT_YES, LENBLOCK, LENLIST)

      CALL TABLE_VARLIST_FILTERS('NDOF_TRUEFLUX','I', ADDCOL_FILTERS,  & 
               VARLIST, LENLIST )
      CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ADDCOL_NDOF_TRUEFLUX(1),  & 
          VARLIST(1:LENLIST)//char(0), ITEXT_YES, LENBLOCK, LENLIST)

    ENDIF

! - - - - -  min and max MJD
    VARLIST = 'MJDMIN:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_MJDMIN, VARLIST, ITEXT_NO,   LENBLOCK, 20 )
    VARLIST = 'MJDMAX:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_MJDMAX, VARLIST, ITEXT_NO,   LENBLOCK, 20 )

! - - - - - MJD_DETECT and MJD_TRIGGER (Sep 2023) - - - - -

    VARLIST = 'MJD_TRIGGER:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_MJD_TRIGGER, VARLIST, ITEXT_NO,   LENBLOCK, 20 )

    VARLIST = 'MJD_DETECT_FIRST:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_MJD_DETECT_FIRST, VARLIST, ITEXT_NO,   LENBLOCK, 20 )
    VARLIST = 'MJD_DETECT_LAST:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNLC8_MJD_DETECT_LAST, VARLIST, ITEXT_NO,   LENBLOCK, 20 )


! - - - -  multi-season variability - - - -
    IF( MULTISEASON_OPTMASK > 0 ) THEN

! Aug 15 2018: include NTOT and NLC in output and TEXT output
      VARLIST = 'NSEASON_TOT:I' // char(0)
      CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NSEASON_TOT, VARLIST,  & 
                         ITEXT_YES, LENBLOCK, 40 )

      VARLIST = 'NSEASON_ACTIVE:I' // char(0)
      CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NSEASON_ACTIVE, VARLIST,  & 
                         ITEXT_YES, LENBLOCK, 40 )

      VARLIST = 'MULTISEASON_CHI2RED_1:F,' //  & 
                  'MULTISEASON_CHI2RED_2:F,' //  & 
                  'MULTISEASON_CHI2RED_3:F,' //  & 
                  'MULTISEASON_CHI2RED_4:F,' //  & 
                  'MULTISEASON_CHI2RED_5:F'  //  & 
                char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, MULTISEASON_CHI2RED(1),  & 
                     VARLIST, ITEXT_NO,     LENBLOCK, 160 )

      VARLIST = 'MULTISEASON_AVGFLUX_1:F,' //  & 
                  'MULTISEASON_AVGFLUX_2:F,' //  & 
                  'MULTISEASON_AVGFLUX_3:F,' //  & 
                  'MULTISEASON_AVGFLUX_4:F,' //  & 
                  'MULTISEASON_AVGFLUX_5:F' //  & 
                char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, MULTISEASON_AVGFLUX(1),  & 
                    VARLIST, ITEXT_NO,    LENBLOCK, 160 )
    ENDIF

! ------------------------------------


    IF ( IFLAG .EQ. 1 ) THEN

! snana variables with no fitting.

       VARLIST = 'TrestMIN:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TRESTMIN, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST = 'TrestMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TRESTMAX, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST = 'TrestRange:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TRESTRANGE, VARLIST, ITEXT_NO, LENBLOCK, 40 )

       VARLIST = 'TobsMIN:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TOBSMIN, VARLIST, ITEXT_NO, LENBLOCK,  40 )
       VARLIST = 'TobsMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TOBSMAX, VARLIST, ITEXT_NO, LENBLOCK,  40 )

       VARLIST = 'TGAPMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_TGAPMAX, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST = 'T0GAPMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_T0GAPMAX, VARLIST, ITEXT_NO, LENBLOCK, 40 )

       VARLIST =  'NFILT_Tmin:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFILT_TRESTMIN, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST =  'NFILT_Tmax:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFILT_TRESTMAX, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST =  'NFILT_Trest2:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFILT_TREST2, VARLIST, ITEXT_NO, LENBLOCK, 40 )

       VARLIST = 'NFILT_SNRMAX:I'  // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFILT_SNRMAX, VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST = 'NFILT_SNRMAX2:I'  // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, ISNLC_NFILT_SNRMAX2, VARLIST, ITEXT_NO, LENBLOCK, 40 )

       VARLIST = 'SNRMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_SNRMAX_SORT(1), VARLIST, ITEXT_NO,   LENBLOCK, 20 )

       VARLIST = 'SNRMAX1:F,SNRMAX2:F,SNRMAX3:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_SNRMAX_SORT(1), VARLIST, ITEXT_YES, LENBLOCK, 40 )

       VARLIST = 'SNRSUM:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNLC_SNRSUM,  VARLIST, ITEXT_YES,   LENBLOCK, 20 )
! ----------
! Below are SNANA variables for which the FIT analog is not
! available here -> book FIT-analogs in snlc_fit.F90

       CALL TABLE_VARLIST_FILTERS('SNRMAX', 'F', ADDCOL_FILTERS, VARLIST, LENLIST )
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_SNRMAX(1),  & 
          VARLIST(1:LENLIST)//char(0),ITEXT_NO, LENBLOCK, LENLIST)

       CALL TABLE_VARLIST_FILTERS('XTMW', 'F', ADDCOL_FILTERS, VARLIST, LENLIST )
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, ADDCOL_XTMW(1),  & 
          VARLIST(1:LENLIST)//char(0), ITEXT_NO, LENBLOCK, LENLIST)

    ELSE

! variables re-evaluated after fit

       VARLIST = 'TrestMIN:F'  // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TRESTMIN_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )
       VARLIST = 'TrestMAX:F'  // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TRESTMAX_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )
       VARLIST = 'TrestRange:F'  // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TRESTRANGE_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )

       VARLIST = 'TobsMIN:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TOBSMIN_FIT,  VARLIST, ITEXT_NO, LENBLOCK, 40 )
       VARLIST = 'TobsMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TOBSMAX_FIT,  VARLIST, ITEXT_NO, LENBLOCK, 40 )

       VARLIST = 'TGAPMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TGAPMAX_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )
       VARLIST = 'T0GAPMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, T0GAPMAX_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )

       VARLIST = 'NFILT_Tmin:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NFILT_TRESTMIN_FIT, VARLIST, ITEXT_NO,   LENBLOCK, 40 )
       VARLIST = 'NFILT_Tmax:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NFILT_TRESTMAX_FIT, VARLIST, ITEXT_NO,   LENBLOCK, 40 )


       VARLIST = 'NFILT_SNRMAX:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NFILT_SNRMAX_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )
       VARLIST = 'NFILT_SNRMAX2:I' // char(0)
       CALL SNTABLE_ADDCOL_int(ID, CBLOCK, NFILT_SNRMAX2_FIT, VARLIST, ITEXT_NO,  LENBLOCK, 40 )

       VARLIST = 'SNRMAX:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNRMAX_SORT_FIT(1), VARLIST, ITEXT_NO,    LENBLOCK, 20 )

       VARLIST = 'SNRMAX1:F,SNRMAX2:F,SNRMAX3:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNRMAX_SORT_FIT(1), VARLIST, ITEXT_YES,  LENBLOCK, 40 )

       VARLIST = 'SNRSUM:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNRSUM_FIT, VARLIST, ITEXT_YES,   LENBLOCK, 20 )

    ENDIF

! - - - - - - private vars - - - - -

! for PRIVATE variables in data files,
! glue together variables in comma-separate VARLIST

    IF ( NVAR_PRIVATE .GT. 0 ) THEN
         CTMP    = PRIVATE_VARNAME(1)(1:MXCHAR_FILEWORD-1) // ' '
         LENNAME = INDEX(CTMP,' ') - 1
         VARLIST = CTMP(1:LENNAME) // ':D'
      do ivar = 2, NVAR_PRIVATE
         CTMP    = PRIVATE_VARNAME(ivar)
         LENNAME = INDEX(CTMP,' ') - 1
         LENLIST = INDEX(VARLIST,' ') - 1
         VARLIST = VARLIST(1:LENLIST) // ',' //  & 
              CTMP(1:LENNAME) // ':D'

      enddo

      LENLIST = INDEX(VARLIST,' ') - 1
      VARLIST = VARLIST(1:LENLIST) // char(0)

      CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, PRIVATE_VALUE(1), VARLIST, ITEXT_NO, LENBLOCK, LENLIST)
    ENDIF     ! end of NPAR_PRIVATE if-block

! --- Sep 23 2017: user variables --------
    IF ( NTABLEVAR_USER > 0 ) THEN
      DO ivar = 1, NTABLEVAR_USER
        VARLIST = TABLEVARNAME_USER(ivar)
        LENLIST = INDEX(VARLIST,' ') - 1
        VARLIST = VARLIST(1:LENLIST) // ':F' // char(0)
        CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, TABLEVALUE_USER(1), VARLIST, ITEXT_YES, LENBLOCK, LENLIST+2)
      ENDDO
    ENDIF

    RETURN
  END SUBROUTINE INIT_TABLE_SNANAVAR


! ======================================
    SUBROUTINE INIT_TABLE_HOSTVAR(ID, IGAL, ITEXT, BLOCK)

! Created Sep 2023
! init table column for HOST_XXX (igal=1) or HOST2_XXX (igal=2)
! 
!  Inputs
!   ID   : table ID
!   IGAL : 1 for best match, >1 for other matches
!   ITEXT: 1 -> write to text also
! 
! 
! Feb 23 2024: pass ITEXT for HOST MAGS so that they appear with
!              SNTABLE_LIST = 'SNANA(text:key,text:host)'
! 
! Jul 13 2024: fix bug in which all gal property indices were 1 instead of IGAL.
! Sep 24 2024: fix so that MAGOBS works for IGAL > 1


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! subroutine args
    INTEGER ID, IGAL, ITEXT   ! (I)
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)

! local var
    CHARACTER PREFIX*8, VARLIST*400, CBLOCK*40
    INTEGER   LP, LENLIST, LENBLOCK, ITEXT_LOCAL

! ------------- BEGIN ----------

    if ( IGAL == 1 ) then
      PREFIX = 'HOST_'
      LP     = 5
    else
      write(PREFIX,20) IGAL
      LP = 6
20      format('HOST',I1,'_')
    endif

    LENBLOCK     = INDEX(BLOCK//' ',' ') - 1
    CBLOCK       = BLOCK(1:LENBLOCK) // char(0)
    LENLIST      = LEN(VARLIST)

! store HOST_OBJID as double precision because
!  * root tree supports Long64_t, but does not read back correctly
!  * sntable_dump extraction will write real*8 as integer.

    VARLIST =  PREFIX(1:LP) // 'OBJID:D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, DSNHOST_OBJID(IGAL), VARLIST,ITEXT, LENBLOCK, 20 )

    VARLIST =  PREFIX(1:LP) // 'ZPHOT:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_ZPHOT(IGAL), VARLIST,ITEXT, LENBLOCK, 40 )
    VARLIST =  PREFIX(1:LP) // 'ZPHOTERR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_ZPHOT_ERR(IGAL), VARLIST,ITEXT, LENBLOCK, 40 )

! 
    if(SNHOST_NZPHOT_Q > 0) then
       VARLIST =  PREFIX(1:LP) // 'QZPHOT:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_QZPHOT_MEAN(IGAL), VARLIST,ITEXT, LENBLOCK, 40 )
       VARLIST =  PREFIX(1:LP) // 'QZPHOTSTD:F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_QZPHOT_STD(IGAL), VARLIST,ITEXT, LENBLOCK, 40 )
    end if

    VARLIST =  PREFIX(1:LP) // 'ZSPEC:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_ZSPEC(IGAL),  & 
                   VARLIST,ITEXT,    LENBLOCK, 40 )
    VARLIST =  PREFIX(1:LP) // 'ZSPECERR:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_ZSPEC_ERR(IGAL),  & 
                    VARLIST,ITEXT,  LENBLOCK, 40 )

    VARLIST = PREFIX(1:LP) // 'RA:D'  // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNHOST8_RA(IGAL),  & 
                VARLIST,ITEXT,    LENBLOCK, 20 )
    VARLIST = PREFIX(1:LP) // 'DEC:D'  // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID, CBLOCK, SNHOST8_DEC(IGAL),  & 
                VARLIST,ITEXT,   LENBLOCK, 20 )

    VARLIST = PREFIX(1:LP) // 'ANGSEP:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_ANGSEP(IGAL),  & 
               VARLIST,ITEXT,   LENBLOCK, 20 )

    VARLIST = PREFIX(1:LP) // 'DDLR:F'  // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_DDLR(IGAL),  & 
                VARLIST,ITEXT,    LENBLOCK, 20 )

    if ( igal == 1 ) then
      VARLIST = 'HOST_CONFUSION:F'  // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_CONFUSION, VARLIST,ITEXT, LENBLOCK, 20 )
    endif

! physical properties
    ITEXT_LOCAL = 1  ! always write properties to text format
    if ( IGAL > 1 ) ITEXT_LOCAL = ITEXT  ! optional for other host matches

    IF ( EXIST_SNHOST_LOGMASS ) THEN
      VARLIST =  PREFIX(1:LP) // 'LOGMASS:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGMASS(IGAL), VARLIST, ITEXT_LOCAL, LENBLOCK, 40 )
      VARLIST =  PREFIX(1:LP) // 'LOGMASS_ERR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGMASS_ERR(IGAL), VARLIST, ITEXT_LOCAL, LENBLOCK, 40 )
    ENDIF

    IF ( EXIST_SNHOST_LOGSFR ) THEN
      VARLIST =  PREFIX(1:LP) // 'LOGSFR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGSFR(IGAL), VARLIST, ITEXT_LOCAL,    LENBLOCK, 40 )
      VARLIST =  PREFIX(1:LP) // 'LOGSFR_ERR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGSFR_ERR(IGAL), VARLIST, ITEXT_LOCAL,    LENBLOCK, 40 )
    ENDIF

    IF ( EXIST_SNHOST_LOGsSFR ) THEN
      VARLIST =  PREFIX(1:LP) // 'LOGsSFR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGsSFR(IGAL), VARLIST, ITEXT_LOCAL,   LENBLOCK, 40 )
      VARLIST =  PREFIX(1:LP) // 'LOGsSFR_ERR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_LOGsSFR_ERR(IGAL), VARLIST, ITEXT_LOCAL,    LENBLOCK, 40 )
    ENDIF

    IF ( EXIST_SNHOST_COLOR ) THEN
      VARLIST =  PREFIX(1:LP) // 'COLOR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_COLOR(IGAL), VARLIST, ITEXT_LOCAL,    LENBLOCK, 40 )
      VARLIST =  PREFIX(1:LP) // 'COLOR_ERR:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SNHOST_COLOR_ERR(IGAL), VARLIST, ITEXT_LOCAL,     LENBLOCK, 40 )
    ENDIF

! Sep 24 2024: update to allow writing MAGOBS for both host matches.
    IF ( EXIST_SNHOST_MAGOBS ) THEN    ! host mags
      CALL TABLE_VARLIST_FILTERS( PREFIX(1:LP)//'MAG', 'F', ADDCOL_FILTERS, VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID,CBLOCK,ADDCOL_SNHOST_MAGOBS(1,IGAL), VARLIST(1:LENLIST)//char(0),  & 
            ITEXT_LOCAL, LENBLOCK, LENLIST)
    ENDIF

    RETURN
  END SUBROUTINE INIT_TABLE_HOSTVAR

! =============================================
    SUBROUTINE INIT_TABLE_SNOBSVAR(ID, BLOCK)

! Mar 2 2015:
! Book epoch info in SNANA table.
! 
! Sep 23 2017: reduce MXOBS down to 1000 so that MJD*8 is handled by hbook.
! Oct 16 2017: include CCDNUM if it is there.
! Mar 03 2018: include SIM_FLUXCAL_HOSTERR(NOBS):F
! Mar 18 2018: include DTOBS(NOBS):F and DTOBS_SAMEFILT(NOBS):F
! 
! Nov 30 2020:
!    + change a few VARNAMEs to match names in FITRES+RESIDUALS
!    + add CHI2FLUX_SIM
!    + add several variables: GAIN, SKYSIG_T ...
!    + remove a few obscure variables
!       (DTOBSLAST, DTOBSLAST_SAMEFILT, FLUXCAL_HOSTERRCALC)
! - - - - - - - - -


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER   ID          ! (I) table id
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)

    INTEGER, PARAMETER :: MXOBS_SNTABLE = 900  !  cannot exceed 1000 for MJD*8

! local var

    INTEGER LENV, LENB, LENTMP
    CHARACTER CBLK*40, CNOBS*14, VARNAME*40
    LOGICAL  IGNORE

    EXTERNAL  & 
          SNTABLE_ADDCOL  & 
         ,SNTABLE_ADDCOL_int  & 
         ,SNTABLE_ADDCOL_flt  & 
         ,SNTABLE_ADDCOL_dbl  & 
         ,SNTABLE_ADDCOL_str
    LOGICAL  IGNOREFILE_fortran

! --------------- BEGIN -------------

    write(6,10) BLOCK, ID
 10   format(T6,'Create BLOCK = ',A,'  for TABLE ID = ', I5)
    call flush(6)

    LENV   = LEN(VARNAME)
    LENB   = INDEX(BLOCK//' ',' ') - 1
    CBLK   = BLOCK(1:LENB) // char(0)

    write(CNOBS,15) MXOBS_SNTABLE
15    format('NOBS[',I4.4,']' )  ! no spaces or commas in string
    VARNAME = CNOBS(1:12)  // char(0)  ! Nobs used in fit + NREJECT
    CALL SNTABLE_ADDCOL_int(ID, CBLK, ISNLC_NEWMJD_STORE,  VARNAME, 0, LENB, LENV)

    VARNAME = 'USEFLAG(NOBS):I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID,CBLK, ISNLC_SNRECON_USE(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'PHOTFLAG(NOBS):I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID,CBLK, ISNLC_PHOTFLAG(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'PHOTPROB(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_PHOTPROB(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'IFILTOBS(NOBS):I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID,CBLK, ISNLC_IFILT_OBS(1),  VARNAME, 0, LENB,LENV)

    VARNAME = 'MJD(NOBS):D' // char(0)
    CALL SNTABLE_ADDCOL_dbl(ID,CBLK, SNLC8_MJD(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'FLUXCAL_DATA(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_FLUXCAL(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'FLUXCAL_DATA_ERR(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_FLUXCAL_ERRTOT(1), VARNAME, 0, LENB,LENV)

! -- observing conditions

    VARNAME = 'ZP(NOBS):F' // char(0)    ! match FITRES+RESID table
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_ZEROPT(1),  VARNAME, 0, LENB,LENV)

    VARNAME = 'ZP_ERR(NOBS):F' // char(0)    ! match FITRES+RESIDS table
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_ZEROPT_ERR(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'PSF(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_PSF_SIG1(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'SKYSIG(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_SKYSIG(1),  VARNAME, 0, LENB,LENV)

    VARNAME = 'SKYSIG_T(NOBS):F' // char(0)    ! Nov 30 2020
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_SKYSIG_T(1), VARNAME, 0, LENB,LENV)

    VARNAME = 'GAIN(NOBS):F' // char(0)    ! Nov 30 2020
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_GAIN(1),  VARNAME, 0, LENB,LENV)

! --- error calculated from PSF+SKY+ZP, as the simulation would do

    VARNAME = 'FLUXCAL_ERRCALC(NOBS):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SNLC_FLUXCAL_ERRCALC(1), VARNAME, 0, LENB,LENV)

! - - - - - - - - - - - - - - - -

    IF ( ISJOB_SIM ) THEN
       VARNAME = 'FLUXCAL_SIM(NOBS):F' // char(0)   ! match FITRES+RESIDS
       CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_EPFLUXCAL(1), VARNAME, 0, LENB,LENV)

       VARNAME = 'CHI2FLUX_SIM(NOBS):F' // char(0)  ! match FITRES+RESIDS
       CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_EPCHI2FLUX(1),  VARNAME, 0, LENB,LENV)

       VARNAME = 'MAGOBS_SIM(NOBS):F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_EPMAGOBS(1),  VARNAME, 0, LENB,LENV)
    ENDIF

    IF ( LSIM_SNANA ) THEN
       VARNAME = 'SIM_FLUXCAL_HOSTERR(NOBS):F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_EPFLUXCAL_HOSTERR(1),  VARNAME, 0, LENB,LENV)
    ENDIF

! include SNR_MAG[mag] if it exists.
    LENTMP  = INDEX(SIMNAME_SNRMON,' ')-1
    IGNORE  = ( IGNOREFILE_fortran(SIMNAME_SNRMON) )
    IF ( LSIM_SNANA .and. (.not.IGNORE)  ) THEN
       VARNAME = SIMNAME_SNRMON(1:LENTMP)// '(NOBS):F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_EPSNRMON(1), VARNAME, 0, LENB,LENV)
    ENDIF

    RETURN
  END SUBROUTINE INIT_TABLE_SNOBSVAR


! =============================================
    SUBROUTINE INIT_TABLE_SIM_MAGOBS(ID, BLOCK)

! Feb 2019
! Book only MODEL-MAG vs. EPOCH (much smaller output than SNANA+EPOCHS)
! MODEL-MAG is interpolated on a grid to enable plotting model colors.
! 
!  NEP_SIM_MODELGRID, SIM_MODELGRID_TOBS, SIM_MODELGRID_MAGOBS


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER   ID          ! (I) table id
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)

! local var

    INTEGER LENV, LENB, IFILT, IFILT_OBS
    LOGICAL IGNORE
    CHARACTER CBLK*40, CNOBS*20, VARNAME*40, cfilt*2

    EXTERNAL  & 
          SNTABLE_ADDCOL  & 
         ,SNTABLE_ADDCOL_int  & 
         ,SNTABLE_ADDCOL_flt  & 
         ,SNTABLE_ADDCOL_dbl  & 
         ,SNTABLE_ADDCOL_str

! --------------- BEGIN -------------

    IF ( .not. LSIM_SNANA ) RETURN ! sim only

    write(6,10) BLOCK, ID
 10   format(T6,'Create BLOCK = ',A,'  for TABLE ID = ', I5)
    call flush(6)

    LENV   = LEN(VARNAME)
    LENB   = INDEX(BLOCK//' ',' ') - 1
    CBLK   = BLOCK(1:LENB) // char(0)

    write(CNOBS,15) MXEP_MODELGRID
15    format('NEP_MODEL[',I4.4,']' )  ! no spaces or commas in string
    VARNAME = CNOBS(1:15)  // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLK, NEP_SIM_MODELGRID,  VARNAME, 0, LENB, LENV)

    VARNAME = 'EPOCH_MODEL(NEP_MODEL):F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID,CBLK, SIM_MODELGRID_TOBS(1),  VARNAME, 0, LENB,LENV)

    DO ifilt = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        cfilt     = filtdef_string(ifilt_obs:ifilt_obs)
        VARNAME   = 'MAGOBS_MODEL_' // cfilt(1:1) //  & 
                          '(NEP_MODEL):F' // char(0)
        CALL SNTABLE_ADDCOL_flt(ID,CBLK,  SIM_MODELGRID_MAGOBS(1,ifilt),  VARNAME, 0, LENB,LENV )
    ENDDO

    RETURN
  END SUBROUTINE INIT_TABLE_SIM_MAGOBS

! =============================================
    SUBROUTINE LOAD_TABLE_SIM_MAGOBS()

! Feb 2019:
!   load SIM_MODELGRID arrays for MODELMAG table.
!   This CWNT is designed to enable plotting model color vs. epoch
! 
! NEP_SIM_MODELGRID, SIM_MODELGRID_TOBS, SIM_MODELGRID_MAGOBS


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IGRID, NEP, IFILT, IFILT_OBS
    REAL MAGOBS, TOBS

! function
    REAL SIM_MAGOBS_INTERP

! ---------- BEGIN ----------

    IF ( .not. LSIM_SNANA ) RETURN ! sim only

! define arbitrary epoch grid

    NEP = 0
    DO IGRID = -40, 100, 2
       NEP  = NEP + 1
       TOBS = FLOAT(IGRID)
       SIM_MODELGRID_TOBS(NEP) = TOBS

       DO ifilt = 1, NFILTDEF_SURVEY
          MAGOBS = SIM_MAGOBS_INTERP(ifilt,TOBS)
          SIM_MODELGRID_MAGOBS(NEP,IFILT) = MAGOBS
       ENDDO
    ENDDO
    NEP_SIM_MODELGRID = NEP
    RETURN
  END SUBROUTINE LOAD_TABLE_SIM_MAGOBS

! ===========================
    REAL FUNCTION SIM_MAGOBS_INTERP(IFILT,Tobs)

! Return interpolated SIM_MAGOBS for inputs:
!   * IFILT  = sparse filter index
!   * Tobs   = MJD - SIM_PEAKMJD


    USE SNDATCOM
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IFILT  ! (I)
    REAL Tobs      ! (I)

    INTEGER NEWMJD, EP, EPMIN, EPMAX, IFILT_OBS, EP0, EP1
    REAL MAGOBS, Tobs_tmp, T0, T1, MAG0, MAG1, frac
    REAL*8 MJD8
! ------------ BEGIN -----------

    MAGOBS = 99.0  ! init value

    T0=-9999. ; T1=-9999. ; EP0 = -9; EP1 = -9
    DO 100 NEWMJD = 1, ISNLC_NEWMJD_STORE
      EPMIN  = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX  = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      DO 101 EP = EPMIN, EPMAX
          IFILT_OBS = ISNLC_IFILT_OBS(ep)
          IF ( IFILT .NE. IFILTDEF_INVMAP_SURVEY(ifilt_obs)) GOTO 101
          MJD8      = SNLC8_MJD(EP)
          Tobs_tmp  = sngl(MJD8) - SIM_PEAKMJD
          MAGOBS    = SIM_EPMAGOBS(ep)
          if ( Tobs_tmp < Tobs ) then
             EP0 = EP; T0 = Tobs_tmp; MAG0 = MAGOBS
          endif
          if ( Tobs_tmp > Tobs .and. T1 < -999. ) THEN
             EP1 = ep; T1 = Tobs_tmp ; MAG1 = MAGOBS ; goto 444
          endif
 101     CONTINUE
 100  CONTINUE

 444  CONTINUE

    IF ( EP0 > 0 .and. EP1 > 0 ) THEN
       FRAC   = (Tobs-T0)/(T1-T0)
       MAGOBS = MAG0 + (MAG1-MAG0) * frac
    ENDIF

    SIM_MAGOBS_INTERP = MAGOBS

    RETURN
  END FUNCTION SIM_MAGOBS_INTERP

! ========================================
    SUBROUTINE INIT_TABLE_SIMVAR(ID,BLOCK)

! Created Feb 02, 2013
! Initialize table for simulated variables.
! 
! Feb 27, 2020: add SIM_HOSTLIB_GALID
! Mar 19, 2020: add SIM_RV  to table
! Sep 11, 2020: add PySEDMODEL = BYOSED or SNEMO
! Dec 28, 2020: write SIM_gammaDM to TEXT file if non-zero
! Dec     2022: remove SIM_TYPE_INDEX
! Jan  8  2023: restore SIM_TYPE_INDEX (forgot why it was removed)
! Dec 23  2023: add SIM_WGT_POPULATION
! 

    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER   ID          ! (I) table id
    CHARACTER BLOCK*(*)   ! (I) name of BLOCK (optional)

! local var

    CHARACTER  & 
          CTMP*(MXCHAR_FILEWORD)  & 
        , VARLIST*400, CBLOCK*40
    INTEGER LENLIST, IPAR, LENBLOCK, LENTMP, LENM, ITEXT

    EXTERNAL  & 
          SNTABLE_ADDCOL  & 
         ,SNTABLE_ADDCOL_int  & 
         ,SNTABLE_ADDCOL_flt  & 
         ,SNTABLE_ADDCOL_dbl  & 
         ,SNTABLE_ADDCOL_str

! ---------------- BEGIN ---------------

    IF ( .NOT. LSIM_SNANA         ) RETURN
    IF ( .NOT. WRTABLEFILE_SIMVAR ) RETURN


    write(6,10) BLOCK, ID
 10   format(T6,'Create BLOCK = ',A,'  for TABLE ID = ', I5)
    call flush(6)

    LENBLOCK = INDEX(BLOCK//' ',' ') - 1
    LENLIST  = LEN(VARLIST)
    CBLOCK   = BLOCK(1:LENBLOCK) // char(0)


    IF ( SIM_SUBSAMPLE_INDEX >= 0 ) THEN
      VARLIST = 'SIM_SUBSAMPLE_INDEX:I' // char(0)
      CALL SNTABLE_ADDCOL_int(ID,CBLOCK,SIM_SUBSAMPLE_INDEX,  VARLIST, 0,     LENBLOCK, LENLIST)
    ENDIF

    VARLIST = 'SIM_GENTYPE:I' //char(0)  ! legacy name
    CALL SNTABLE_ADDCOL(ID, CBLOCK, SIM_GENTYPE,  VARLIST, 1,   LENBLOCK, LENLIST)

    VARLIST = 'SIM_TEMPLATE_INDEX:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_TEMPLATE_INDEX,  VARLIST,1,  LENBLOCK, LENLIST)


    VARLIST = 'SIM_LIBID:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_LIBID,  VARLIST, 1,  LENBLOCK, LENLIST)
    VARLIST = 'SIM_NGEN_LIBID:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_NGEN_LIBID,  VARLIST, 1,  LENBLOCK, LENLIST)

    VARLIST = 'SIM_NOBS_UNDEFINED:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_NOBS_UNDEFINED, VARLIST,0, LENBLOCK, LENLIST)

    VARLIST = 'SIM_SEARCHEFF_MASK:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_SEARCHEFF_MASK, VARLIST, 0, LENBLOCK, LENLIST)

    VARLIST = 'SIM_MAGSMEAR_COH:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_MAGSMEAR_COH, VARLIST, 0,  LENBLOCK, LENLIST)

    VARLIST = 'SIM_MWEBV:F' // char(0)  ! added Nov 26 2016
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_MWEBV,  & 
                      VARLIST, 0,   LENBLOCK, LENLIST)

    VARLIST = 'SIM_ZCMB:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_REDSHIFT_CMB,  & 
                    VARLIST, 1,     LENBLOCK, LENLIST)

    VARLIST = 'SIM_ZHELIO:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_REDSHIFT_HELIO,  & 
                      VARLIST, 0,    LENBLOCK, LENLIST)

    VARLIST = 'SIM_ZHOST:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_REDSHIFT_HOST,  & 
                   VARLIST, 0,      LENBLOCK, LENLIST)

    VARLIST = 'SIM_ZHOST_MATCH:F' // char(0)   ! Jun 5 2025
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_REDSHIFT_HOST_MATCH,  & 
                   VARLIST, 1,      LENBLOCK, LENLIST)

    VARLIST = 'SIM_ZFLAG:I' // char(0)
    CALL SNTABLE_ADDCOL_int(ID, CBLOCK, SIM_REDSHIFT_FLAG,  & 
                    VARLIST, 1,     LENBLOCK, LENLIST)

    VARLIST = 'SIM_VPEC:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_VPEC,  & 
                     VARLIST, 1,    LENBLOCK, LENLIST)

    VARLIST = 'SIM_DLMAG:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_DLMAG,  & 
                     VARLIST, 1,    LENBLOCK, LENLIST)

! Feb 4 2025: add SIM_LENSDMU to text fitres file
    VARLIST = 'SIM_LENSDMU:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_LENSDMU,  & 
                     VARLIST, 1,    LENBLOCK, LENLIST)

    IF ( SIM_MUSHIFT .NE. 0.0 ) THEN   ! hope that first events does not have SIM_MUSHIFT=0
      VARLIST = 'SIM_MUSHIFT:F' // char(0)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_MUSHIFT, VARLIST, 1,  LENBLOCK, LENLIST)
    ENDIF

    VARLIST = 'SIM_PEAKMJD:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_PEAKMJD,  & 
                     VARLIST, 1,     LENBLOCK, LENLIST)

! - -
    VARLIST = 'SIM_WGT_POP' // ':F'   // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_WGT_POPULATION,  & 
                     VARLIST, 1,     LENBLOCK, LENLIST)

    LENTMP = INDEX(SIMNAME_SHAPEPAR,' ') - 1  ! x1, delta, dm15
    VARLIST = SIMNAME_SHAPEPAR(1:LENTMP) // ':F'   // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_SHAPEPAR, VARLIST, 1,     LENBLOCK, LENLIST)

    LENTMP = INDEX(SIMNAME_COLORPAR,' ') - 1  ! c, AV
    VARLIST = SIMNAME_COLORPAR(1:LENTMP) // ':F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_COLORPAR, VARLIST, 1,     LENBLOCK, LENLIST)

    LENTMP = INDEX(SIMNAME_SHAPELAW,' ') - 1   ! SIM_alpha
    VARLIST = SIMNAME_SHAPELAW(1:LENTMP) // ':F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_SHAPELAW, VARLIST, 1,     LENBLOCK, LENLIST)

    LENTMP = INDEX(SIMNAME_COLORLAW,' ') - 1   ! SIM_beta
    VARLIST = SIMNAME_COLORLAW(1:LENTMP) // ':F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_COLORLAW,  VARLIST, 1,   LENBLOCK, LENLIST)

! add SIM_gammaDM if not zero.
! Sep 7 2021: always write SIM_gammaDM in case first value is 0
!    followed by non0-zero values.
    VARLIST = 'SIM_gammaDM:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_SALT2gammaDM,  & 
                      VARLIST, 1,   LENBLOCK, LENLIST )

    VARLIST = 'SIM_x0:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_SALT2x0,  & 
                      VARLIST, 1,   LENBLOCK, LENLIST)

    VARLIST = 'SIM_mB:F' // char(0)
    CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_SALT2mb,  & 
                      VARLIST, 1,   LENBLOCK, LENLIST)

! add SIM_AV only if not already defined.
! This is for the SALT2 model with intinsic color (c) and
! external scatter (SIM_AV) from dust.
! Aug 2 2022: perform test only for SALT2 model.

    if ( SIM_MODEL_INDEX .EQ. MODEL_SALT2 ) then
      if ( SIMNAME_COLORPAR(1:6) .NE. 'SIM_AV' ) then
        VARLIST = 'SIM_AV:F' // char(0)
        CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_AV, VARLIST, 1, LENBLOCK, LENLIST)
        VARLIST = 'SIM_RV:F' // char(0)
        CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_RV,  VARLIST, 1,  LENBLOCK, LENLIST)
      endif
    endif

!  SIMSED model params; glue SIMSED_ prefix (Feb 2014)
!  Apr 28 2015: set DOTEXT flag to 1

    DO IPAR    = 1, NPAR_SIMSED
       CTMP    = SIMSED_PARNAME(ipar)
       LENTMP  = INDEX(CTMP,' ') - 1
       VARLIST = 'SIMSED_' // CTMP(1:LENTMP) // ':F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLOCK, SIMSED_PARVAL(ipar),  VARLIST, 1, LENBLOCK, LENLIST)
    ENDDO


! PySEDMODEL model params (BYOSED, SNEMO, ...)
    DO IPAR    = 1, NPAR_PySEDMODEL
       LENM    = INDEX(PySEDMODEL_NAME,' ') - 1
       CTMP    = PySEDMODEL_PARNAME(ipar)
       LENTMP  = INDEX(CTMP,' ') - 1
       VARLIST = PySEDMODEL_NAME(1:LENM) // '_'  & 
                   // CTMP(1:LENTMP) // ':F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLOCK, PySEDMODEL_PARVAL(ipar),  VARLIST, 1, LENBLOCK, LENLIST)
    ENDDO

! Feb 2018: LCLIB model params
    IF ( SIM_MODEL_INDEX == MODEL_LCLIB ) THEN
      CALL TABLE_VARLIST_FILTERS('SIM_TEMPLATEMAG','F',ADDCOL_FILTERS, VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID,CBLOCK, SIM_TEMPLATEMAG(1),  & 
                VARLIST(1:LENLIST)//char(0), 0, LENBLOCK, LENLIST)

      DO IPAR    = 1, NPAR_LCLIB
         CTMP    = LCLIB_PARNAME(ipar)
         LENTMP  = INDEX(CTMP,' ') - 1
         VARLIST = 'LCLIB_' // CTMP(1:LENTMP) // ':F' // char(0)
         CALL SNTABLE_ADDCOL_flt(ID,CBLOCK, LCLIB_PARVAL(ipar),  VARLIST,1,   LENBLOCK, LENLIST)
      ENDDO
    ENDIF

! SIM_HOSTLIB params

    ITEXT = 0
    if ( WRTABLEFILE_HOST_TEXT ) ITEXT=1
    VARLIST = 'SIM_HOSTLIB_GALID' // ':D' // char(0) ! Feb 2020
       CALL SNTABLE_ADDCOL_dbl(ID,CBLOCK, DSIM_HOSTLIB_GALID,  VARLIST, ITEXT,    LENBLOCK, LENLIST)

    DO IPAR    = 1, NPAR_SIM_HOSTLIB
       CTMP    = SIM_HOSTLIB_PARNAME(ipar)
       LENTMP  = INDEX(CTMP,' ') - 1
       VARLIST = 'SIM_HOSTLIB_' // CTMP(1:LENTMP) // ':F' // char(0)
       CALL SNTABLE_ADDCOL_flt(ID,CBLOCK, SIM_HOSTLIB_PARVAL(ipar),  & 
              VARLIST, ITEXT,    LENBLOCK, LENLIST)
    ENDDO

! - - - - - optional LCWIDTH for simulated MAGOBS -----
    if ( OPTSIM_LCWIDTH > 0 ) then
      CALL TABLE_VARLIST_FILTERS('SIM_LCWIDTH','F', ADDCOL_FILTERS,  & 
                 VARLIST, LENLIST)
      CALL SNTABLE_ADDCOL_flt(ID, CBLOCK, SIM_LCWIDTH(1),  & 
            VARLIST(1:LENLIST)//char(0), 0, LENBLOCK, LENLIST)
    endif

    RETURN
  END SUBROUTINE INIT_TABLE_SIMVAR



! ===================================================
    SUBROUTINE TABLE_VARLIST_FILTERS(PREFIX,CAST,FILTERS,  & 
                  VARLIST,LENLIST)

! Created Sep 28, 2012
! 
! For input PREFIX  and defined FILTERS,  return
! comma-separate list of variable names for column-wise ntuple.
! Example for FILTERS=gri,
! 
! VARLIST = '[PREFIX]_g,[PREFIX]_r,[PREFIX]_i'
! 
! If a filter is case-insensitive duplicate (i.e., u & U)
! then append '2' at the end of the variable name.
! The returned VARLIST is used in the HBNAME calls,
! and hence case-insensitive duplicates must be distiguished
! for ntuples.
! 
! Feb  1, 2013: include ':F' after each variable.
! Feb  5, 2013: return LENLIST
! Feb 15, 2013: set max of MXFILT = 30 filters.
! Jun 24, 2013: ignore remapped filter(s)
! Oct 21, 2014: add CAST as input argument
! Jul 13, 2020: MXFILT -> 15
! ------------------------------------


    USE SNPAR
    USE FILTCOM

    IMPLICIT NONE

    CHARACTER  & 
         PREFIX*(*)    &  ! (I) prefix for variable names
        ,CAST*(*)      &  ! (I) 'F', 'I', 'D'
        ,FILTERS*(*)   &  ! (I) list of filters; i.e, 'gri'
        ,VARLIST*(*)  ! (O) output list of varnames.

    INTEGER LENLIST    ! (O) length of VARLIST

! local var

    INTEGER NFILT, ifilt, lenv, MXFILT, LENMAX, LENSUM
    INTEGER IFILTOBS
    LOGICAL LTMP
    character varname*32, cfilt*1, c1err*76, c2err*76

! function
    LOGICAL DUPLICATE_FILTER
    INTEGER FILTINDX

! ---------------- BEGIN ----------

    LENMAX = LEN(VARLIST)  ! abort if larger than this
    LENSUM = 0
!  cxxx      MXFILT = 30
    MXFILT = 15

    NFILT = INDEX(FILTERS,' ' ) - 1
    IF ( NFILT > MXFILT ) NFILT = MXFILT
    VARLIST = ''

! -----------
    DO 100 ifilt = 1, NFILT
       cfilt     = FILTERS(ifilt:ifilt)
       LTMP      = DUPLICATE_FILTER(ifilt,FILTERS)

       IFILTOBS = FILTINDX(cfilt // ' ')
       IF( IFILTOBS_REPLACE(IFILTOBS) .NE. IFILTOBS) GOTO 100

       varname   = PREFIX // '_' //  cfilt
       LENV = INDEX(VARNAME,' ') - 1
       if ( LTMP ) VARNAME = VARNAME(1:LENV) // '2'

       LENV = INDEX(VARNAME,' ') - 1
       LENSUM = LENSUM + LENV + 2  ! for error message only

       IF ( ifilt .eq. 1 ) then
          varlist = VARNAME(1:LENV) // ':' // CAST
       else if ( LENSUM < LENMAX ) then
          lenlist = INDEX(varlist,' ') - 1
          varlist = varlist(1:LENLIST) // ',' //  & 
                  VARNAME(1:LENV) // ':' // CAST
       endif

100   CONTINUE

! ---------------------
! abort if required size (LENSUM) of VARLIST exceed size passed (LENMAX)

    IF ( LENSUM > LENMAX ) THEN
      write(c1err,661) LENSUM, LENMAX
661     format('LENSUM=',I4, ' but LEN(VARLIST)=',I4,  & 
                ' is too short.')
      write(c2err,662) PREFIX, NFILT, FILTERS(1:NFILT)
662     format('PREFIX=',A,2x,'NFILT=',I3,'(',A,'+)'  )
      CALL MADABORT("TABLE_VARLIST_FILTERS", c1err, c2err)
    ENDIF

    LENLIST = INDEX(VARLIST,' ') - 1

    RETURN
  END SUBROUTINE TABLE_VARLIST_FILTERS

! ======================================
    SUBROUTINE TABLE_ADDCOL_LOAD()

! Created Feb 2017
! Load filter-dependent ADDCOL_XXX arrays.
! Use nominal filter indices, or REMAPED filters if
! &SNLCINP namelist input SNTABLE_FILTER_REMAP is set.


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE PKMJDCOM
    USE TRUECHI2COM

    IMPLICIT NONE

    INTEGER ifilt, NFILT, IFILTOBS, IFILTOBS_REMAP, IFILT_REMAP
    INTEGER LEN
    CHARACTER  CFILTOBS*2, CFILTOBS_REMAP*2, FNAM*18
    LOGICAL DO_NOMINAL, DO_REMAP, VALID_BAND
    LOGICAL USE, USEFILT_REMAP(MXFILT_ALL)
    REAL    SNRMAX, VAL_OLD, VAL_NEW
! function
    INTEGER FILTINDX

! ------------- BEGIN --------------

    FNAM = 'TABLE_ADDCOL_LOAD'

    DO_NOMINAL = ( NFILT_REMAP_TABLE == 0 )
    DO_REMAP   = ( .NOT. DO_NOMINAL)

    NFILT   = NFILTDEF_SURVEY

    IF ( DO_REMAP) THEN
      DO  ifilt=1, MXFILT_ALL
        USEFILT_REMAP(ifilt) = .FALSE.
      ENDDO
    ENDIF

    DO IFILT_REMAP = 1, NFILT_REMAP_TABLE
         ADDCOL_SNRMAX(IFILT_REMAP)     = -9.0
         ADDCOL_XTMW(IFILT_REMAP)       = -9.0
         ADDCOL_FLUXCALMAX(IFILT_REMAP)        = -9.0
         ADDCOL_FLUXCALMAX_ERR(IFILT_REMAP)     = -9.0
         ADDCOL_PROB_TRUEFLUX(IFILT_REMAP)      = -9.0
         ADDCOL_NDOF_TRUEFLUX(IFILT_REMAP)      =  0
         ADDCOL_SNHOST_MAGOBS(IFILT_REMAP,1)    = -9.0
         ADDCOL_SNHOST_MAGOBS(IFILT_REMAP,2)    = -9.0
         ADDCOL_SNHOST_SBFLUXCAL(IFILT_REMAP)   = -9.0
         ADDCOL_SNHOST_SBMAG(IFILT_REMAP)       = -9.0
    ENDDO

    DO 100 ifilt = 1, NFILT
       IFILTOBS = IFILTDEF_MAP_SURVEY(ifilt)
       CFILTOBS = filtdef_string(ifiltobs:ifiltobs)

       IF ( DO_NOMINAL ) then
         ADDCOL_SNRMAX(ifilt)           = SNLC_SNRMAX_FILT(ifilt)
         ADDCOL_XTMW(ifilt)             = SNLC_MWXT_FLUXFRAC(ifilt)
         ADDCOL_FLUXCALMAX(ifilt)       = SNLC_FLUXCALMAX(ifilt)
         ADDCOL_FLUXCALMAX_ERR(ifilt)   = SNLC_FLUXCALMAX_ERR(ifilt)
         ADDCOL_PROB_TRUEFLUX(ifilt)    = PROB_TRUEFLUX(ifilt)
         ADDCOL_NDOF_TRUEFLUX(ifilt)    = NDOF_TRUEFLUX(ifilt)
         ADDCOL_CHI2_FITPKMJD(ifilt)    = CHI2_FITPKMJD(ifilt)
         ADDCOL_SNHOST_MAGOBS(ifilt,1)  = SNHOST_MAGOBS(ifilt,1)
         ADDCOL_SNHOST_MAGOBS(ifilt,2)  = SNHOST_MAGOBS(ifilt,2)
         ADDCOL_SNHOST_SBFLUXCAL(ifilt) = SNHOST_SBFLUXCAL(ifilt)
         ADDCOL_SNHOST_SBMAG(ifilt)     = SNHOST_SBMAG(ifilt)
       ELSE
           ! DO_REMAP
         CALL FILTER_REMAP_FETCH(IFILTOBS,   &  ! (I)
                 IFILTOBS_REMAP,IFILT_REMAP)  ! (O)

         SNRMAX         = SNLC_SNRMAX_FILT(ifilt)
         VALID_BAND     = (SNRMAX > -8.999)
         if ( .NOT. VALID_BAND ) GOTO 100

         IFILT_REMAP    = IFILT_REMAP + 1  ! C -> fortran index
         CFILTOBS_REMAP = ADDCOL_FILTERS(ifilt_remap:ifilt_remap)
         USE            = USEFILT_REMAP(IFILTOBS_REMAP)

! check that each band is mapped.
         if ( IFILTOBS_REMAP<0 .or. IFILT_REMAP<0 ) then
           C1ERR = 'Missing REMAP for Band= ' // CFILTOBS(1:1) //  & 
                '  for CID=' // SNLC_CCID(1:ISNLC_LENCCID)
           C2ERR = 'Check SNTABLE_FILTER_REMAP in &SNLCINP'
           CALL MADABORT(FNAM, c1err, c2err)
         endif

         USEFILT_REMAP(IFILTOBS_REMAP) = .TRUE.

         VAL_OLD = ADDCOL_SNRMAX(IFILT_REMAP)
         VAL_NEW = SNLC_SNRMAX_FILT(ifilt)
         ADDCOL_SNRMAX(IFILT_REMAP) = MAX(VAL_OLD,VAL_NEW)

         ADDCOL_XTMW(IFILT_REMAP) = SNLC_MWXT_FLUXFRAC(ifilt)
         ADDCOL_PROB_TRUEFLUX(IFILT_REMAP) = PROB_TRUEFLUX(ifilt)
         ADDCOL_NDOF_TRUEFLUX(IFILT_REMAP) = NDOF_TRUEFLUX(ifilt)

         VAL_OLD = ADDCOL_FLUXCALMAX(IFILT_REMAP)
         VAL_NEW = SNLC_FLUXCALMAX(ifilt)
         ADDCOL_FLUXCALMAX(IFILT_REMAP) = MAX(VAL_OLD,VAL_NEW)

         VAL_OLD = ADDCOL_CHI2_FITPKMJD(IFILT_REMAP)
         VAL_NEW = CHI2_FITPKMJD(ifilt)
         ADDCOL_CHI2_FITPKMJD(IFILT_REMAP) = VAL_OLD + VAL_NEW

         ADDCOL_SNHOST_MAGOBS(IFILT_REMAP,1) =  & 
                  SNHOST_MAGOBS(ifilt,1)      ! any band
         ADDCOL_SNHOST_MAGOBS(IFILT_REMAP,2) =  & 
                  SNHOST_MAGOBS(ifilt,2)      ! any band
         ADDCOL_SNHOST_SBFLUXCAL(IFILT_REMAP) =  & 
                  SNHOST_SBFLUXCAL(ifilt)   ! any band
         ADDCOL_SNHOST_SBMAG(IFILT_REMAP) =  & 
                  SNHOST_SBMAG(ifilt)   ! any band

       ENDIF

100   CONTINUE

    RETURN
  END SUBROUTINE TABLE_ADDCOL_LOAD

! ======================================
    SUBROUTINE TABLE_STRING_TERMINATION(OPT)

! Created Jan 2014 by R.Kessler
! For ROOT, terminate strings with char(0) to be used by C++/root code.
! Note that the list of table string-variables is hard-wired;
! initial set is SNLC_CCID and SNLC_FIELDLIST.
! 
! 
! OPT > 0 --> add termination
! OPT < 0 --> remove termination
! 
! Oct 20 2020: fix index for removing termination char


    USE SNDATCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT  ! (I) integer option
    INTEGER LENCID, LENFLD

! ------------ BEGIN --------------

    IF ( USE_TABLEFILE_ROOT ) THEN
       LENCID  = INDEX(SNLC_CCID,     ' ') - 1
       LENFLD  = INDEX(SNLC_FIELDLIST,' ') - 1

       IF ( OPT > 0 ) THEN
         ! add termination
         SNLC_CCID      = SNLC_CCID(1:LENCID)      // char(0)
         SNLC_FIELDLIST = SNLC_FIELDLIST(1:LENFLD) // char(0)
       ELSE
         ! remove termination char
         SNLC_CCID(LENCID:LENCID)      = ' '
         SNLC_FIELDLIST(LENFLD:LENFLD) = ' '
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE TABLE_STRING_TERMINATION

! ===========================================
    SUBROUTINE TABLE_MARZ(IDTABLE,IFLAG)
! 
! Created May 2020 by R.Kessler
! Driver for MARZ table with host-galaxy spectra.
! 
! May 3 2021: include GALID as part of NAME
! Nov 16 2023: update text to handle 12 digit GALID instead of 9 digits
! ------------------------------

    USE SNPAR
    USE SNCUTS
    USE CTRLCOM
    USE SNLCCOM
    USE SNLCINP_NML
    USE SNFILECOM
    USE SPECCOM
    USE SNSIMCOM

    IMPLICIT NONE

    INTEGER IDTABLE  ! (I) table ID
    INTEGER IFLAG    ! (I) see IFLAG_XXX params

! local var

    LOGICAL   DO_FILL
    INTEGER   LENNAME, LENFMT, ISPEC, LENCCID, LENz, LENGALID
    CHARACTER FNAM*12, NAME*40, TEXTFMT*20, TEXTFMT_forC*20
    CHARACTER CCID_forC*40, cGALID*24, cz*20
    REAL      z
    INTEGER*8 GALID
    EXTERNAL SNTABLE_CREATE, SNTABLE_FILL, SPECPAK_DATA

! -------------------- BEGIN -----------------
    FNAM = 'TABLE_MARZ'

#if defined(SNFIT) || defined(PSNID)
    c1err = 'Cannot write MARZ spectra file with fit code'
    c2err = 'MARZFILE_OUT works only with snana.exe '
    CALL MADABORT(FNAM, c1err, c2err)
#endif

    IF ( IFLAG .EQ. IFLAG_INI ) THEN
       NAME     = 'MARZ' // char(0)
       LENNAME  = INDEX(NAME,    ' ') - 1
       TEXTFMT_forC = 'FITS' // char(0)
       LENFMT       =  4
       CALL SNTABLE_CREATE(IDTABLE, NAME, TEXTFMT_forC,  & 
                LENNAME,LENFMT)  ! C fun

       RETURN
    ENDIF

! --------------------------
! fill table.

    DO_FILL = LSNCUTS
    IF ( DO_FILL ) THEN
       ISPEC = 1
       LENCCID   = ISNLC_LENCCID
       CCID_forC = SNLC_CCID(1:LENCCID)
       if ( LSIM_SNANA ) then
          z      = SIM_REDSHIFT_HELIO
          GALID  = SIM_HOSTLIB_GALID
          cGALID ='GALID=None'

          write(cz,41) z
 41         format('z=',F5.3)
          if ( SIM_HOSTLIB_GALID > 0 ) then
             write(cGALID,44) SIM_HOSTLIB_GALID
 44            format('GALID=', I12.12 )
          endif

          LENz     = INDEX(cz,' ') - 1
          LENGALID = INDEX(cGALID,' ') - 1

          write(CCID_forC,40) SNLC_CCID(1:LENCCID),  & 
                cz(1:LENz), cGALID(1:LENGALID)
40          format(A, '[', A, ',', A, ']' )
          LENCCID = LENCCID + LENz + LENGALID + 3


       endif

       write(6,10) CCID_forC(1:LENCCID)
 10      format(T8,'SPECPLOT: pack CID=',A,' for MARZ table.')

       CCID_forC = CCID_forC(1:LENCCID) // char(0)

       CALL RDSPEC_DRIVER(ISPEC)  ! Apr 5 2021

       CALL SPECPAK_DATA(  & 
               CCID_forC  & 
             , ID_SPECTRUM(ispec)  & 
             , MJD_SPECTRUM(ispec)  & 
             , TOBS_SPECTRUM(ispec)     &  ! May 7 2019
             , TEXPOSE_SPECTRUM(ispec)  & 
             , NLAMBIN_SPECTRUM(ispec)  & 
             , LAMMIN_SPECTRUM  & 
             , LAMMAX_SPECTRUM  & 
             , FLAM_SPECTRUM  & 
             , FLAMERR_SPECTRUM  & 
             , LENCCID )

        CALL SPECPAK_FILL(CCID_forC, LENCCID)
    ENDIF

    RETURN
  END SUBROUTINE TABLE_MARZ


! =============================================
    LOGICAL FUNCTION DUPLICATE_FILTER(i1,FILTERS)

! Created Sep 29, 2012
! Return TRUE if any previous filter in FILTERS
! is defined after using LOCASE on both filters.
! 
! Example 1: i1=3 and FILTERS='ugrizUGRIZ'
!            i1=3 corresponds to r, and previous filters are unique.
!                --> returns FALSE since there are no duplicates.
! 
! Example 2: i1=8 and FILTERS='ugrizUGRIZ'
!            i1=8 corresponds to R, which is the same as 'r' after
!            converting to lower case --> returns TRUE.
! 
! ------------------------------------------

    INTEGER i1             ! (I) check this filter in the list
    CHARACTER FILTERS*(*)  ! (I) list of filster; i.e, 'gri'

! local var

    INTEGER i2
    character cfilt1*2, cfilt2*2, cfilt1_lower*2, cfilt2_lower*2

! ----------------- BEGIN ------------

    DUPLICATE_FILTER = .FALSE.

    if ( i1 .eq. 1 ) return

    cfilt1      = FILTERS(i1:i1)
    CALL LOCASE(cfilt1, cfilt1_lower)

    DO i2 = 1, i1-1
       cfilt2      = FILTERS(i2:i2)
       CALL LOCASE(cfilt2, cfilt2_lower)
       if ( cfilt1_lower .EQ. cfilt2_lower ) then
          DUPLICATE_FILTER = .TRUE.
          return
       endif
    ENDDO

    RETURN
  END FUNCTION DUPLICATE_FILTER

! ================================================
    SUBROUTINE INIT_GALextinct()

! Created Sep 19 2013 by R.Kessler
! 
! - Hierarchy of OPT_MWCOLORLAW (and RV) for data, in order of priority
!    -> OPT_MWCOLORLAW         (if >=0 -> read from &SNLCINP namelist)
!    -> OPT_MWCOLORLAW_DEFAULT
! 
! - Hierarchy of OPT_MWCOLORLAW for sim
!    -> OPT_MWCOLORLAW         (if >=0 -> read from &SNLCINP namelist)
!    -> SIMOPT_MWCOLORLAW      (from sim files)
! 
! 
! - ABORT on invalid option
! 
! - print info for colorLaw and MWEBV-updates
! 
! ------------------------------------------


    USE SNPAR
    USE CTRLCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE SNSIMCOM

    IMPLICIT NONE

    LOGICAL   LDATA, SETMWCL_USER, SETEBV_USER, SETRV_USER
    INTEGER   LENTEXT
    CHARACTER nameOpt*20, CTEXT*60
    EXTERNAL  text_MWoption

! ---------------- BEGIN -----------

    CALL PRBANNER("INIT_GALextinct")

    IF ( RESTORE_DES5YR ) OPT_MWCOLORLAW = -99 ! restore approx Fitzpatrick99 CL

    LDATA         = .NOT. LSIM_SNANA
! xxx mark delete May 28 2025 SETMWCL_USER  = ( OPT_MWCOLORLAW  > -1   )
    SETMWCL_USER  = ( OPT_MWCOLORLAW  > -998 )  ! allow -99 for approx Fitzpatrick (May 28 2025)
    SETEBV_USER   = ( OPT_MWEBV       > -8   )
    SETRV_USER    = ( RV_MWCOLORLAW   > -0.01 )

! ---------------------------------

    IF ( LDATA ) THEN

!     ---------- check RV ---------
       if ( SETRV_USER ) then
          write(6,22) RV_MWCOLORLAW, 'User (&SNLCINP)'
       else
          RV_MWCOLORLAW = RV_MWCOLORLAW_DEFAULT
          write(6,22) RV_MWCOLORLAW, 'Default'
       endif

!     --------- check color law ------------
       if ( SETMWCL_USER ) then
          write(6,20) OPT_MWCOLORLAW, 'User (&SNLCINP)'
       else
          OPT_MWCOLORLAW = OPT_MWCOLORLAW_DEFAULT
          write(6,20) OPT_MWCOLORLAW, 'Default'
       endif

    ENDIF


    IF ( LSIM_SNANA ) THEN

       if ( SETRV_USER ) then
          write(6,22) RV_MWCOLORLAW, 'User (&SNLCINP)'
       else
          ! use simulated RV if not set by user
          RV_MWCOLORLAW = SIM_MWRV
          write(6,22) RV_MWCOLORLAW, 'Simulation'
       endif

       if ( SETMWCL_USER ) then
          write(6,20) OPT_MWCOLORLAW, 'User (&SNLCINP)'
       else
          ! use simulated color law if not set by user
          OPT_MWCOLORLAW = SIMOPT_MWCOLORLAW
          write(6,20) OPT_MWCOLORLAW, 'Simulation'
       endif

    ENDIF

 20   format(T5, 'OPT_MWCOLORLAW -> ',I3   , ' from ', A)
 22   format(T5, 'RV_MWCOLORLAW  -> ',F5.2 , ' from ', A)


    nameOpt = "COLORLAW" // char(0)
    call text_MWoption(nameOpt, OPT_MWCOLORLAW, CTEXT, 8,40)

    LENTEXT = INDEX(CTEXT,char(0)) - 1
    write(6,30) OPT_MWCOLORLAW, CTEXT(1:LENTEXT)
 30   format(T5,'OPT_MWCOLORLAW =  ', I4,' -> ', A )
    call flush(6)

! ---------------------------------------------
! check for option to modify MWEBV;

    print*,' '


    IF ( LDATA ) THEN
       if ( SETEBV_USER ) then
          write(6,21) OPT_MWEBV, 'User (&SNLCINP)'
       else
       ! use default MWEBV map if not set by user
          OPT_MWEBV = OPT_MWEBV_DEFAULT
          write(6,21) OPT_MWEBV, 'Default'
       endif
    ENDIF

 21   format(T5, 'OPT_MWEBV -> ',I3,' from ', A)

    IF ( LSIM_SNANA ) THEN
       if ( SETEBV_USER ) then
          write(6,21) OPT_MWEBV, 'User (&SNLCINP)'
       else
       ! use simulated MWEBV value if option is not set by user
          OPT_MWEBV = 1       ! don't change  MWEBV in sim-data file
          write(6,21) SIMOPT_MWEBV, 'Simulation'
       endif
    ENDIF

    nameOpt = "EBV" // char(0)
    call text_MWoption(nameOpt, OPT_MWEBV, CTEXT, 3,40)
    LENTEXT = INDEX(CTEXT,char(0)) - 1
    write(6,40) OPT_MWEBV, CTEXT(1:LENTEXT)
 40   format(T5,'OPT_MWEBV =  ',I3,' -> ', A )
    call flush(6)

    RETURN
  END SUBROUTINE INIT_GALextinct

! =======================
! =======================

! ===========================================
    INTEGER FUNCTION NEAREST_IFILTDEF_WRAPPER(  & 
         OPT           &  ! (I) option; see below
        ,ifiltdef      &  ! (I) observer filter index
        ,rank_want     &  ! (I) 1=nearest, 2=2nd nearest, 3=3rd nearest
        ,Z             &  ! (I) redshift
        ,LAMDIF_MIN    &  ! (O) min lam-distance (A)
          )

! 
! Created Nov 2022
! Temporary wrapper to select either LEGACY/fortran function,
! or REFACTORED/C function. This function should be removed
! after the transistion to refactored code.


    USE SNDATCOM
    USE SNANAFIT
    USE FILTCOM
    USE SNLCINP_NML

    IMPLICIT NONE

    INTEGER OPT, IFILTDEF, RANK_WANT
    REAL    z,  LAMDIF_MIN
    REAL*8  z8, LAMDIF_MIN8

    CHARACTER FNAM*32
    INTEGER NEAREST_IFILT_REST     ! lefacy fortran func
    INTEGER NEAREST_IFILTDEF_REST  ! refactored C code
! ------------- BEGIN -------------

       z8 = DBLE(z)
       FNAM = 'NEAREST_IFILTDEF_WRAPPER' // char(0)
       NEAREST_IFILTDEF_WRAPPER =  & 
         NEAREST_IFILTDEF_REST(OPT, IFILTDEF,RANK_WANT,  & 
                               z8, FNAM, LAMDIF_MIN8, 32)
       LAMDIF_MIN = SNGL(LAMDIF_MIN8)


    RETURN
  END FUNCTION NEAREST_IFILTDEF_WRAPPER


! =======================================
    SUBROUTINE FILTER_UPDATE_DRIVER()
! 
! Nov 8, 2010 R.Kessler
! 
! --------------


    USE SNPAR
    USE CTRLCOM
    USE SNLCCOM
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

! local var

    INTEGER IFILT, IFILT_OBS, OPT
    LOGICAL LSAME, LUPD, ISBX
    CHARACTER C1ERR*80, C2err*80

! function
    LOGICAL ISBXFILT
! ------------- BEGIN ----------

    OPT = OPT_FILTER_UPDATE
    IF ( OPT .LE. 0 ) RETURN

    write(6,20) SNLC_CCID(1:ISNLC_LENCCID)
 20   format(T5,'***** UPDATE FILTER TRANSMISSION FOR CID = ',  & 
             A, ' ***** ')


    LUPD = ( OPT .EQ. OPT_FILTUPD_EACHSN ) .or.  & 
             ( OPT .EQ. OPT_FILTUPD_MAP    )

! determine path with filter update; set global FILTER_UPDATE_DIR
! Make sure to do this outside the IFILT loop to avoid redundant
! loops through the maps.
    IF ( LUPD ) CALL GET_FILTER_UDPATE_DIR(SNLC_CCID)

! -----------------------------------------
! update filter transmissions based on option

    DO 10 IFILT = 1, NFILTDEF_SURVEY

      if ( ISNLC_NEPOCH_FILT(ifilt) .LE. 0 ) GOTO 10

      ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
      DOFLAG_FILTER_UPDATE(ifilt_obs) = OPT_FILTER_UPDATE

      ISBX      = ISBXFILT(ifilt_obs,'', 'OBS' )

      LSAME = (OPT .EQ. OPT_FILTUPD_SAMEFILT)  .or. ISBX

      if ( LSAME  ) then
          ! fancy way to use same filter response for each SN
         CALL FILTER_UPDATE_SAMEFILT(ifilt)
      else if ( LUPD ) then
         CALL FILTER_UPDATE_EACHSN(ifilt)

      else
         write(c1err,60) OPT
60         format('Invalid  OPT_FILTER_UPDATE = ', I3 )
         c2err = 'Check OPT_FILTUPD_*  parameters in snana.F90'
         CALL MADABORT("FILTER_UPDATE_DRIVER", C1ERR, C2ERR)
      endif

 10   CONTINUE


! ----------------------------------------------------
!  Finally, update ZP offsets for all filters.

    IF ( LUPD ) THEN
      CALL FILTER_UPDATE_ZPOFF()
    ENDIF

! later add stuff for K-corrections ...

    RETURN
  END SUBROUTINE FILTER_UPDATE_DRIVER

! ====================================
    SUBROUTINE FILTER_UPDATE_ZPOFF()

! Created July 2011
! Update MAGOBS_SHIFT_ZP_FILT(ifilt_obs) based on the
! zeropoint (AB) offsets in text-file $FDIR/ZPOFF.DAT .
! If this ZPOFF.DAT file does not exist, then just set
! the offsets to the user offsets.
! 
! -------------


    USE SNPAR
    USE CTRLCOM
    USE SNLCCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE FILTCOM
    USE FILTUPDCM
! CDE,PARSECOM.

    INTEGER IFILT, IFILT_OBS
    REAL ZPOFF_UPD(MXFILT_ALL)

! --------------- BEGIN -------------

     DO ifilt  = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        ZPOFF_UPD(ifilt_obs) = 0.0
     ENDDO

! update ZP-shift for each filter (i.e., add to user-shift)
    DO ifilt  = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       MAGOBS_SHIFT_ZP_FILT(ifilt_obs) =  & 
         MAGOBS_SHIFT_ZP_USER(ifilt_obs) + ZPOFF_UPD(ifilt_obs)
    ENDDO

    CALL DMP_ZPOFF

    RETURN
  END SUBROUTINE FILTER_UPDATE_ZPOFF

! ================================================
    SUBROUTINE GET_FILTER_UDPATE_DIR(CCID)

! Created May 18, 2012
! Find and load global  FILTER_UPDATE_DIR for this CCID.
! Used for filter-update option.


    USE SNPAR
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

    CHARACTER CCID*(*)  ! (I) SN name

    INTEGER LPDIR, LSDIR, LPATH, LCC, LMAP, IMAP1, IMAP2
    CHARACTER MAPDIR*100
! functions
    INTEGER IMAP2_FUPDMATCH

! -------------------- BEGIN --------------

    LCC   = INDEX(CCID,' ') - 1
    LPDIR = INDEX(PREFIX_UPD_FILTDIR,' ') -1
    LSDIR = INDEX(SUFFIX_UPD_FILTDIR,' ') -1
    LPATH = INDEX(FILTER_UPDATE_TOPDIR,' ') - 1


    if ( OPT_FILTER_UPDATE .EQ. OPT_FILTUPD_EACHSN ) then
       FILTER_UPDATE_DIR = FILTER_UPDATE_TOPDIR(1:LPATH) // '/'  & 
              // PREFIX_UPD_FILTDIR(1:LPDIR)  & 
              // CCID(1:LCC)  & 
              // SUFFIX_UPD_FILTDIR(1:LSDIR)

    else if ( OPT_FILTER_UPDATE .EQ. OPT_FILTUPD_MAP ) then

! search 2nd map for CCID match to get IMAP1

       IMAP2   = IMAP2_FUPDMATCH(CCID)
       IMAP1   = MAP2_FILTER_UPDATE_PTRMAP1(IMAP2)
       MAPDIR  = MAP1_FILTER_UPDATE_SUBDIR(IMAP1)
       LMAP    = INDEX(MAPDIR,' ') - 1
       FILTER_UPDATE_DIR = FILTER_UPDATE_TOPDIR(1:LPATH) // '/'  & 
                   // MAPDIR(1:LMAP)

    endif

    RETURN
  END SUBROUTINE GET_FILTER_UDPATE_DIR

! ================================================
    SUBROUTINE FILTER_UPDATE_EACHSN(ifilt)

! Nov 8, 2010
! Read/update filter response for this sparse 'ifilt'
! Store in FILTOBS_TRANS_UPD array using same lambda
! binning as in KCOR file. Updated lambda binning
! can be different since it gets interpolated to the
! original grid.
! 
! Nov 2022: refactor to use new & old kcor/calib code.
! Dec 2023: remove fortran OPEN and rely on rd2column to read gzipped or unzipped files.
! 


    USE SNPAR
    USE CTRLCOM
    USE SNLCCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE FILTCOM
    USE FILTUPDCM
! +CDE,PARSECOM.

! subroutine args
    INTEGER IFILT     ! (I) sparse indices

! local var

    INTEGER  & 
         IFILT_OBS, L0, LCC, LREF, LSN, LDIR, LTMPF, LCOM  & 
        ,LPSN, LPREF, LSSN, LSREF, ilam, NLAM_UPD, NLAM_STORE  & 
        ,istat, MXLAM, OPT_INTERP, OPT_FRAME, OPT_RD2COL

    CHARACTER  & 
         FILTFILE_SN*(MXCHAR_FILENAME)    &  ! SN tranmission vs. lambda
        ,FILTFILE_REF*(MXCHAR_FILENAME)   &  ! idem for primary ref
        ,TMPFILE*(MXCHAR_FILENAME)  & 
        ,CCID*(MXCHAR_CCID)  & 
        ,cfilt*1  & 
        ,c1err*80  & 
        ,c2err*80  & 
        ,comment*40

    REAL*8  & 
         LAM  & 
        ,FLAM_UPD(MXLAMBIN_FILT)     &  ! temp lambda array
        ,FTRANS_UPD(MXLAMBIN_FILT)   &  ! temp filterTrans array
        ,Trans_upd  & 
        ,FLAM_STORE(MXLAMBIN_FILT)  &  ! store trans on this lam array
        ,FTRANSSN_STORE(MXLAMBIN_FILT)  & 
        ,FTRANSREF_STORE(MXLAMBIN_FILT)

    LOGICAL  LEXIST_SN, LEXIST_REF, UPD_DUMP

! define utilities from sntools.c
    INTEGER  RD2COLUMNFILE
    EXTERNAL RD2COLUMNFILE

    REAL*8   INTERP_1DFUN
    EXTERNAL INTERP_1DFUN
    EXTERNAL LOAD_FILTERTRANS_CALIB, GET_CALIB_FILTERLAM

! ---------------- BEGIN --------------

    UPD_DUMP = .FALSE.  ! internal dump logical

    ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
    cfilt     = filtdef_string(ifilt_obs:ifilt_obs)
    LCC       = INDEX(SNLC_CCID, ' '     ) - 1
    CCID      = SNLC_CCID(1:LCC)


     OPT_FRAME = OPT_FILTOBS - 1
     CALL  GET_CALIB_FILTERLAM(OPT_FRAME, IFILT_OBS,  & 
                            NLAM_STORE, FLAM_STORE)


! - - - - -
    if (UPD_DUMP) then
      print*,' '
      print*, 'XXXXXX ----------------------------------------- '
      print*, 'XXXXXX FILTER_UPDATE_EACHSN  DUMP start for ',  & 
             CCID(1:LCC),'-', cfilt
      call flush(6)
    end if

    OPT_INTERP = 1   ! 1=linear;  2=quadratic

    comment = CCID(1:LCC) // '-filter-update' // char(0)
    LCOM    = INDEX(comment,' ') - 1

!     get string lengths for file name

    LPSN  = INDEX(PREFIX_UPD_TRANSSN,' ') -1
    LSSN  = INDEX(SUFFIX_UPD_TRANSSN,' ') -1

    LPREF = INDEX(PREFIX_UPD_TRANSREF,' ') -1
    LSREF = INDEX(SUFFIX_UPD_TRANSREF,' ') -1

    L0    = INDEX(FILTER_UPDATE_PATH, ' ' ) - 1
    LDIR = INDEX(FILTER_UPDATE_DIR,' ') - 1

! construct full filename of SN and REF-transmission files.

    TMPFILE = PREFIX_UPD_TRANSSN(1:LPSN)  & 
                 // cfilt // SUFFIX_UPD_TRANSSN(1:LSSN)
    LTMPF = INDEX(TMPFILE,' ') - 1
    FILTFILE_SN = FILTER_UPDATE_DIR(1:LDIR)  & 
                       // '/' // TMPFILE(1:LTMPF)

    TMPFILE = PREFIX_UPD_TRANSREF(1:LPREF)  & 
                // cfilt // SUFFIX_UPD_TRANSREF(1:LSREF)
    LTMPF = INDEX(TMPFILE,' ') - 1
    FILTFILE_REF = FILTER_UPDATE_DIR(1:LDIR)  & 
                // '/' // TMPFILE(1:LTMPF)


    if ( UPD_DUMP ) then
      print*, '++++ UPDATING FILTER '  & 
                ,'IFILT = ', ifilt_obs,' '  & 
               , cfilt(1:1), ' with NLAM '  & 
               , NLAM_STORE , ' for SN ', CCID
      call flush(6)
     end if


! Use IOSTAT option to figure out which of the files
! FILTFILE_REF  and FILTFILE_SN exist

    OPT_RD2COL = 1  ! rd2columun does NOT abort on missing file

!     check if FILTFILE_REF exists?


    MXLAM   = MXLAMBIN_FILT
    LREF    = INDEX(FILTFILE_REF,' ') - 1
    TMPFILE = FILTFILE_REF(1:LREF) // char(0)
    istat   = rd2columnFile (TMPFILE, MXLAM, NLAM_UPD,  & 
                       FLAM_UPD, FTRANS_UPD, OPT_RD2COL, LREF)  ! (out)
    LEXIST_REF = (istat > 0 )

    if (UPD_DUMP) then
      print*, '  LEXIST_REF = ', LEXIST_REF,' for ', CCID
      call flush(6)
    end if

! - - - - - -
!    now interpolate to the original grid
!    for files which exist

    IF ( LEXIST_REF ) THEN
      DO 20 ilam = 1, NLAM_STORE
          FILTOBS_TRANSREF_UPD(ilam,IFILT_OBS) = 0.0

          LAM = FLAM_STORE(ilam)
          FTRANSREF_STORE(ilam) = 0.0

          if ( LAM .LT. FLAM_UPD(1)        ) GOTO 20
          if ( LAM .GT. FLAM_UPD(NLAM_UPD) ) GOTO 20

          Trans_upd = interp_1dfun(OPT_INTERP, LAM,  & 
                        NLAM_UPD, FLAM_UPD, FTRANS_UPD, comment, LCOM)

          FILTOBS_TRANSREF_UPD(ilam,IFILT_OBS) = sngl(Trans_upd)
          FTRANSREF_STORE(ilam)  = Trans_upd

20      CONTINUE  ! ILAM loop over original lambda binning
      if (UPD_DUMP) then
        print*, "REF interp done"
        call flush(6)
      end if

    ENDIF  ! end LEXIST_REF


! ------------------------------------
!     FILTFILE_SN exists?

    LSN     = INDEX(FILTFILE_SN,' ') - 1
    MXLAM   = MXLAMBIN_FILT
    TMPFILE = FILTFILE_SN(1:LSN) // char(0)
    istat   = rd2columnFile (TMPFILE, MXLAM, NLAM_UPD,  & 
                       FLAM_UPD, FTRANS_UPD, OPT_RD2COL, LSN)  ! (out)
    LEXIST_SN = (istat > 0 )

! - - - - - - -
!   INTERPOLATE to get FILTOBS_TRANSSN_UPD

    IF (LEXIST_SN) THEN
!       print*, ' Interpolate SN file', NLAM_STORE
      DO 70 ilam = 1, NLAM_STORE
          FILTOBS_TRANSSN_UPD(ilam,IFILT_OBS) = 0.0

          LAM = FLAM_STORE(ilam)
          FTRANSSN_STORE(ilam) = 0.0

          if ( LAM .LT. FLAM_UPD(1)        ) GOTO 70
          if ( LAM .GT. FLAM_UPD(NLAM_UPD) ) GOTO 70

          Trans_upd = interp_1DFUN(OPT_INTERP, LAM,  & 
                   NLAM_UPD, FLAM_UPD, FTRANS_UPD, comment, LCOM )

          FILTOBS_TRANSSN_UPD(ilam,IFILT_OBS) = sngl(Trans_upd)
          FTRANSSN_STORE(ilam) = Trans_upd

70      CONTINUE  ! ILAM loop over original lambda binning
      if (UPD_DUMP) then
        print*, 'SN Interp done for', NLAM_STORE, CCID
        call flush(6)
      end if
    ENDIF


!    If neither file exists, abort program


    if (.NOT.( LEXIST_REF .OR. LEXIST_SN) ) then
      c1err = 'Could not update filter-trans for CID = '  & 
             // CCID(1:LCC) // '-' // cfilt
      c2err = 'FILTDIR=' // FILTER_UPDATE_DIR(1:LDIR)
      CALL MADABORT('FILTER_UPDATE_EACHSN', c1err, c2err )
    endif

!     IF either SN or REF file does not exist, copy the other
!     file for the transmission functions


!     SN File does not exist. IF REF file exists, SN trans= Ref Trans
!                             IF REF file does not exist, SN trans undef
!                             but should have been aborted above

    if (.NOT. LEXIST_SN ) then
      DO 170 ilam = 1, NLAM_STORE
          FILTOBS_TRANSSN_UPD(ilam,IFILT_OBS) =  & 
            FILTOBS_TRANSREF_UPD(ilam,IFILT_OBS)

          FTRANSSN_STORE(ilam) = FTRANSREF_STORE(ilam)
170     CONTINUE  ! ILAM loop over original lambda binning
    end if


! - - - - - - - - -
!     REF File does not exist. IF SN file exists, Ref trans= SN Trans
!                             IF SN file does not exist, Ref trans undef
!                             but should have been aborted above
    IF (.NOT. LEXIST_REF ) THEN
      DO 270 ilam = 1, NLAM_STORE
          LAM  = FILTOBS_LAMBDA(ilam,IFILT_OBS)
          FILTOBS_TRANSREF_UPD(ilam,IFILT_OBS) =  & 
            FILTOBS_TRANSSN_UPD(ilam,IFILT_OBS)

          FTRANSREF_STORE(ilam) = FTRANSSN_STORE(ilam)
270      CONTINUE  ! ILAM loop over original lambda binning
    ENDIF

    CALL LOAD_FILTERTRANS_CALIB(OPT_FRAME, IFILT_OBS, NLAM_STORE,  & 
               FLAM_STORE, FTRANSSN_STORE )

! -----------------------------
    if ( UPD_DUMP ) then
      print*, '%%%% ', NLAM_STORE
      if (mod(NLAM_STORE,2).EQ.0) then
         ilam = NLAM_STORE/2
      else
       ilam = (NLAM_STORE+1)/2
      end if
      LAM  = FILTOBS_LAMBDA(ilam,IFILT_OBS)
      print*, 'xxxxx', NLAM_STORE  & 
              , sngl(FTRANSSN_STORE(ilam))  & 
              , sngl(FTRANSREF_STORE(ilam))

      print*, 'xxxx SN =' , CCID
      print*, 'xxxx LEXIST_REF =', LEXIST_REF
      print*, 'xxxx LEXIST_SN  =', LEXIST_SN
      print*, 'xxxx FILE_REF   =', FILTFILE_REF
      print*, 'xxxx FILE_SN    =', FILTFILE_SN

      print*, 'XXXXXX FILTER_UPDATE_EACHSN  DUMP end for ', CCID
      call flush(6)
    endif


    RETURN
  END SUBROUTINE FILTER_UPDATE_EACHSN


! =======================================
    INTEGER FUNCTION IMAP2_FUPDMATCH(CCID)

! Created May 19, 2012
! Search array IMAP2_FILTER_UPDATE_CCID(imap2)
! to find CCID match; return IMAP2 = index of array.
! Used to match CCID to sparse index associated
! with the filter-subdir.


    USE SNPAR
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

    CHARACTER CCID*(*)  ! (I) SN name to match

    INTEGER IMAP2, L0, L1
    CHARACTER CCID_TMP*(MXCHAR_CCID), C1ERR*72, C2ERR*72

! ---------------- BEGIN ------------------

    IMAP2_FUPDMATCH = -9
    L0 = INDEX(CCID,' ') - 1

    DO 100 IMAP2 = 1, NMAP2_FILTER_UPDATE

       CCID_TMP  = MAP2_FILTER_UPDATE_CCID(IMAP2)
       L1 = INDEX(CCID_TMP // ' ',' ') - 1

       IF ( L0 .NE. L1 ) GOTO 100
       IF ( CCID(1:L0) .NE. CCID_TMP(1:L1) ) GOTO 100

! if we get here then return
       IMAP2_FUPDMATCH = IMAP2
       RETURN

100   CONTINUE


! if we get here then abort.
    write(c1err,61) CCID(1:8)
61    format('Could not find MAP2_FILTER_UPDATE match for CID=',A)
    c2err = 'Check FILTER.INFO file.'
    CALL MADABORT('IMAP2_FUPDMATCH', c1err, c2err )


    RETURN
  END FUNCTION IMAP2_FUPDMATCH

! ================================================
    SUBROUTINE FILTER_UPDATE_SAMEFILT(ifilt)

    USE SNPAR
    USE CTRLCOM
    USE SNLCCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE FILTCOM
    USE FILTUPDCM
! +CDE,PARSECOM.

    INTEGER ifilt  ! (I) sparse SN & filter index

! local variables

    INTEGER ilam, NLAM, ifilt_obs

! ------------ BEGIN -------------

    ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
    NLAM      = NLAMBIN_FILTOBS(ifilt_obs)

    DO 100 ilam = 1, NLAM

        FILTOBS_TRANSSN_UPD(ilam,IFILT_OBS) =  & 
          FILTOBS_TRANS(ilam,IFILT_OBS)

        FILTOBS_TRANSREF_UPD(ilam,IFILT_OBS) =  & 
          FILTOBS_TRANS(ilam,IFILT_OBS)
100   CONTINUE

    RETURN
  END SUBROUTINE FILTER_UPDATE_SAMEFILT

! ==========================================
    SUBROUTINE SET_ZPOFF()
! 
! Created July 17, 2011 by R.Kessler
! Called from RDKCOR to update MAGOBS_SHIFT_ZP_FILT(ifilt_obs)
! for each filter. Update based on FILTOBS_ZPOFF_SNPHOT(ifilt_obs)
! that is read either from the KCOR file or from the ZPOFF_FILE.
! 
! Feb 22, 2012:
!  Use MAGOBS_SHIFT_ZP_USER(ifilt_obs)  so that final
!  SHIFT_ZP is OK even if RDKCOR is called multiple times.
! -------------


    USE SNPAR
    USE CTRLCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE SNFILECOM
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IFILT, IFILT_OBS

! ------------- BEGIN -----------

! update ZP-shift for each filter (i.e., add to user-shift)
    DO ifilt  = 1, NFILTDEF_SURVEY
        ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
        MAGOBS_SHIFT_ZP_FILT(ifilt_obs) =  & 
          MAGOBS_SHIFT_ZP_USER(ifilt_obs) +  & 
          FILTOBS_ZPOFF_SNPHOT(ifilt_obs)
    ENDDO

! check for lambda-dependent MAGOBS_SHIFT (Sep 30, 2015)
    CALL MAGOBS_SHIFT_UPDATE()

! dump ZPOFF values vs. filter
    CALL DMP_ZPOFF

    RETURN
  END SUBROUTINE SET_ZPOFF

! ==========================================
    SUBROUTINE DMP_ZPOFF

! --------------------------
! Created Jul 17, 2011
! Dump zeropoint offsets to screen.
! 
! Feb 2016: dump only obs-filters in kcor file; ingore un-defined filters.
! Nov 2022: separate filter loop over OBS and REST frame so that it works
!           for both legacy and refac kcor.
! 
! July 2024: clarify stdout comment that ZP offsets are from ZPOFF.DAT file
!            and from user MAGOBS_SHIFT_ZP ... not from kcor-input file.
! -------------------------

    USE SNPAR
    USE CTRLCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    REAL      ZPOFF_TOT
    INTEGER   ifilt_obs, ifilt_rest, ifilt
    character CFILT*2

! ------------------ BEGIN ------------

    print*,' '
    print*,'    Extra ZP offsets from ZPOFF.DAT + ' //  & 
                  ' SNLCINP MAGOBS_SHIFT_ZP:'

    DO 100 ifilt = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       if ( ifilt_obs .EQ. IFILT_BESS_BX ) goto 100
       cfilt     = filtdef_string(ifilt_obs:ifilt_obs)
       ZPOFF_TOT = MAGOBS_SHIFT_ZP_FILT(ifilt_obs)
       write(6,801) 'MAGOBS', cfilt(1:1), ZPOFF_TOT
100   CONTINUE

    DO 200 ifilt = 1, NFILTDEF_REST
       ifilt_rest = IFILTDEF_MAP_REST(ifilt)
       if ( ifilt_rest .EQ. IFILT_BESS_BX ) goto 200
       cfilt     = filtdef_string(ifilt_rest:ifilt_rest)
       ZPOFF_TOT = MAGREST_SHIFT_ZP_FILT(ifilt_rest)
       write(6,801) 'MAGREST', cfilt(1:1), ZPOFF_TOT
200   CONTINUE

801     format(T10,'Will apply extra ',  & 
           A,'_SHIFT_ZP(',A,') = ', F7.4 )

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE DMP_ZPOFF


! ======================================
    LOGICAL FUNCTION FILTBTEST(MSK8,ibit)
! Mar 29, 2011
! Return true if 1 <= IBIT < 90 is set in MSK8.
! Note that LSB=1 (not 0)
! 

    INTEGER   IBIT    ! (I) bit to set
    INTEGER*8 MSK8(2) ! (I/O) 8-byte mask to set

    INTEGER IBIT2
    INTEGER, PARAMETER ::  MXBIT =  60

    LOGICAL LTMP
! ------------ BEGIN ------------

    if ( IBIT .LE. MXBIT ) then
      LTMP = BTEST(MSK8(1),IBIT-1)
    else if ( IBIT .LT. 2*MXBIT ) then
      IBIT2  = IBIT-MXBIT
      LTMP   = BTEST(MSK8(2),IBIT2-1)
    endif
    FILTBTEST = LTMP
    RETURN
  END FUNCTION FILTBTEST

! ======================================
    SUBROUTINE FILTBSET(MSK8,ibit)
! Mar 29, 2011
! set bit IBIT(1 to 2*MXBIT) for long int MSK8.
! Note that the input argument is also the output
! with IBIT set.

    INTEGER   IBIT    ! (I) bit to set
    INTEGER*8 MSK8(2) ! (I/O) 8-byte mask to set

    INTEGER IBIT2
    INTEGER, PARAMETER :: MXBIT =  60

! ------------ BEGIN ------------

    if ( IBIT .LE. MXBIT ) then
      MSK8(1) = IBSET(MSK8(1),IBIT-1)
    else if ( IBIT .LT. 2*MXBIT ) then
      IBIT2   = IBIT-MXBIT
      MSK8(2) = IBSET(MSK8(2),IBIT2-1)
    endif

    RETURN
  END SUBROUTINE FILTBSET

! =======================================
    SUBROUTINE MAGOBS_SHIFT_UPDATE()

! Created Sep 20 2015
! Loop over filters and apply optional shift-vs-lambda
! for ZP and PRIMARY mag. Allows coherent mag-shifts
! with just a few parameters instead of giving a mag
! shift for each band.
! See &SNLCINP parameters MAGOBS_SHIFT_[ZP,PRIMARY]_PARAMS(3)
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER IFILT, IFILTDEF
    REAL MAGOBS_SHIFT_USRFUN  ! function
! --------------- BEGIN -----------

    DO ifilt = 1, NFILTDEF_SURVEY
       IFILTDEF = IFILTDEF_MAP_SURVEY(ifilt)

       magobs_shift_zp_filt(ifiltdef) =  & 
         magobs_shift_zp_filt(ifiltdef) +  & 
                MAGOBS_SHIFT_USRFUN(ifiltdef, 'ZP') ! Sep 2015

       magobs_shift_primary_filt(ifiltdef) =  & 
         magobs_shift_primary_filt(ifiltdef) +  & 
                MAGOBS_SHIFT_USRFUN(ifiltdef, 'PRIMARY') ! Sep 2015
    ENDDO

    RETURN
  END SUBROUTINE MAGOBS_SHIFT_UPDATE

! =======================================
    REAL FUNCTION MAGOBS_SHIFT_USRFUN(IFILTDEF, WHAT )

! Created Sep 30 2015
! Compute shift for what='ZP' or what='PRIMARY' for this IFILTDEF.
! Use SNLCINP parameters MAGOBS_SHIFT_[ZP,PRIMARY]_PARAMS to
! determine polynomial function of wavelength. This allows computing
! correlated systematics among filters using a function of wavelength
! so that you don't have to enter a shift for each band separately.
! 
! BEWARE: POLYNOMIAL FUNCTION OF LAMBDA IN MICRONS (NOT ANGSTROMS)
! 
! ----------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM

    IMPLICIT NONE

    INTEGER    IFILTDEF   ! (I) absolute filter index
    CHARACTER  WHAT*(*)   ! (I) 'ZP' or 'PRIMARY'

! local var

    REAL LAM_um, SQLAM_um, COEF(0:2)
    CHARACTER FNAM*20

! ------------------- BEGIN ----------------

    FNAM = 'MAGOBS_SHIFT_USRFUN'
    MAGOBS_SHIFT_USRFUN = 0.0

    IF ( WHAT(1:2) .EQ. 'ZP' ) THEN
       COEF(0) = MAGOBS_SHIFT_ZP_PARAMS(1)
       COEF(1) = MAGOBS_SHIFT_ZP_PARAMS(2)
       COEF(2) = MAGOBS_SHIFT_ZP_PARAMS(3)
    ELSE IF ( WHAT(1:7) .EQ. 'PRIMARY' ) THEN
       COEF(0) = MAGOBS_SHIFT_PRIMARY_PARAMS(1)
       COEF(1) = MAGOBS_SHIFT_PRIMARY_PARAMS(2)
       COEF(2) = MAGOBS_SHIFT_PRIMARY_PARAMS(3)
    ELSE
       C1ERR = 'Invalid WHAT=' // WHAT
       C2ERR = 'for MAGOBS_SHIFT_XXX'
       CALL  MADABORT(FNAM, C1ERR, C2ERR)
    ENDIF

    LAM_um   = FILTOBS_LAMAVG(ifiltdef)/1.0E4
    SQLAM_um = LAM_um * LAM_um

    MAGOBS_SHIFT_USRFUN =  & 
            COEF(0)  & 
          + COEF(1) * LAM_um  & 
          + COEF(2) * SQLAM_um

! xxxxxxxxxxxxxxxxx
!      write(6,66) WHAT(1:2), MAGOBS_SHIFT_USRFUN, IFILTDEF, LAM_um
! 6    format(' xxx ', A2, ' = ', F7.4,' for IFILTDEF=',I3,
!     &      3x,'LAM=',F7.3,' um')
!      CALL FLUSH(6)
! xxxxxxxxxxxxxxxxx

    RETURN
  END FUNCTION MAGOBS_SHIFT_USRFUN


! =======================================
    SUBROUTINE FILTER_UPDATE_INIT(IERR)
! 
! If FILTER_UPDATE_PATH is defined, read FILTER.INFO
! file for instructions on updating the filters for
! each SN.
! 
! Jun 15 2016:  CALL ENVreplace(FILTER_UPDATE_PATH)
! 

    USE SNPAR
    USE CTRLCOM
    USE SNCUTS
    USE SNLCINP_NML
    USE SNFILECOM
! +CDE,PARSECOM.
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

    INTEGER IERR

! local var

    INTEGER  & 
         LEN_UPDPATH, LEN_FULLPATH, L1, LEN_ROOT  & 
        ,iwd, NWD, IERROPEN, ifilt, LWD0, LWD1

    character  & 
         FILTINFO_FILE*(MXCHAR_FILENAME)  & 
        ,NAME_forC*(MXCHAR_FILENAME)  & 
        ,cwd*60, cwd_next*60, msg*100, FNAM*30

    INTEGER  STORE_PARSE_WORDS
    EXTERNAL STORE_PARSE_WORDS

! --------------- BEGIN -------------

    IERR = 0
    FNAM = 'FILTER_UPDATE_INIT'

    OPT_FILTER_UPDATE    = 0
    DO ifilt=1, MXFILT_ALL
      DOFLAG_FILTER_UPDATE(ifilt) = 0
    ENDDO

    PREFIX_UPD_TRANSSN  = ' '
    PREFIX_UPD_TRANSREF = ' '
    PREFIX_UPD_FILTDIR  = ' '
    SUFFIX_UPD_TRANSSN  = ' '
    SUFFIX_UPD_TRANSREF = ' '
    SUFFIX_UPD_FILTDIR  = ' '
    FILTINFO_UPD_SN     = .FALSE.
    FILTINFO_UPD_REF    = .FALSE.

    IF ( FILTER_UPDATE_PATH .EQ. ' ' ) RETURN

    CALL ENVreplace(FILTER_UPDATE_PATH)

    msg =  & 
         "FILTER_UPDATE_INIT: Prepare SN-dependent Filter response"
    CALL PRBANNER(msg(1:60))

! -------------------------------------------------------
! first try $SNDATA_ROOT/filters/[FILTER_UPDATE_PATH]/
! if not there, then try absolute path.

    LEN_UPDPATH  = INDEX(FILTER_UPDATE_PATH,' ' ) -  1
    LEN_ROOT     = INDEX(SNDATA_ROOT,' ') - 1

    FILTER_UPDATE_TOPDIR = SNDATA_ROOT(1:LEN_ROOT)  & 
          // '/filters/' // FILTER_UPDATE_PATH(1:LEN_UPDPATH)
    LEN_FULLPATH  = INDEX(FILTER_UPDATE_TOPDIR,' ') - 1


! check if/where FILTER.INFO exists ...
    FILTINFO_FILE = FILTER_UPDATE_TOPDIR(1:LEN_FULLPATH)  & 
                       // '/FILTER.INFO'
    OPEN(UNIT=LUNTMP, FILE = FILTINFO_FILE,  & 
           IOSTAT = IERROPEN, STATUS='OLD')
    IF ( IERROPEN .EQ. 0 ) THEN
      CLOSE(UNIT = LUNTMP)
    ELSE
! try absolute path
      LEN_FULLPATH = LEN_UPDPATH
      FILTER_UPDATE_TOPDIR  = FILTER_UPDATE_PATH
      FILTINFO_FILE = FILTER_UPDATE_PATH(1:LEN_UPDPATH)  & 
                 // '/FILTER.INFO'
    ENDIF

    L1 = INDEX(FILTINFO_FILE,' ') - 1
    print*,' Update Filter Response for each SN from : '
    print*,' ', FILTINFO_FILE(1:L1)

    NAME_forC = FILTINFO_FILE(1:L1)//char(0)
    NWD = STORE_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,NAME_forC,  & 
               FNAM//char(0), L1, 30 )


! parse enough to know which type of FILTER.INFO file we have;
! set global  OPT_FILTER_UPDATE and call PARSE1_XXX or PARSE2_XXX
! to do the appropriate parsing.

    DO 10 iwd = 1, min(1000,NWD)

        CALL get_PARSE_WORD_fortran(iwd+0, cwd,      LWD0 )
        CALL get_PARSE_WORD_fortran(iwd+1, cwd_next, LWD1 )

        if ( cwd .EQ. ' ' ) goto 10

        if ( CWD .EQ. 'FILTER_SUBDIR:' ) THEN
           CALL PARSE1_FILTER_INFO(NWD)
           OPT_FILTER_UPDATE  = OPT_FILTUPD_EACHSN
           RETURN
        ENDIF

        if ( CWD .EQ. 'FILTER_PATH:' ) THEN
           CALL PARSE2_FILTER_INFO(NWD)
           OPT_FILTER_UPDATE  = OPT_FILTUPD_MAP
           RETURN
        ENDIF

10    CONTINUE


! if we get here then abort.

    c1err = 'Cannot parse FILTER.INFO file'
    c2err = FILTER_UPDATE_TOPDIR(1:LEN_FULLPATH)
    CALL MADABORT("FILTER_UPDATE_INIT", c1err, c2err)

    RETURN
  END SUBROUTINE FILTER_UPDATE_INIT

! ======================================
    SUBROUTINE PARSE1_FILTER_INFO(NWD)

! parse FILTER.INFO file that has a different filter set
! for each SN.


    USE SNPAR
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

    INTEGER NWD ! (I) number of words in file
    INTEGER IWD, istar, LEN

    character  & 
         FILTER_SUBDIR*(MXCHAR_PATH)  & 
        ,cwd*(MXCHAR_FILEWORD)  & 
        ,cwd_next*(MXCHAR_FILEWORD)  & 
        ,copt*(MXCHAR_FILEWORD)  & 
        ,c1err*80, c2err*80

    LOGICAL LDUMP

! ----------- BEGIN ----------

    DO 10 iwd = 1, NWD-1

        CALL get_PARSE_WORD_fortran(iwd+0, cwd,      LEN)
        CALL get_PARSE_WORD_fortran(iwd+1, cwd_next, LEN)

        if ( cwd .EQ. ' ' ) goto 10

        if ( cwd .EQ. 'FILTER_SUBDIR:' ) then
          FILTER_SUBDIR = cwd_next
          istar = Index(FILTER_SUBDIR,'{SNID}')
          PREFIX_UPD_FILTDIR = FILTER_SUBDIR(1:istar-1)
          SUFFIX_UPD_FILTDIR = FILTER_SUBDIR(istar+6:99)
        endif

        if ( cwd .EQ. 'FILETRANS_SN:' ) then
          CALL SET_FILTINFO_UPD("SN", cwd_next)
        endif
        if ( cwd .EQ. 'FILETRANS_REF:' ) then
          CALL SET_FILTINFO_UPD("REF", cwd_next)
        endif

        if ( cwd(1:7) .EQ. 'OPTION:' ) then
           copt = cwd_next
           if ( copt(1:7) .EQ. 'DEFAULT' ) then

           else if ( copt(1:8) .EQ. 'SAMEFILT' ) then
             OPT_FILTER_UPDATE  = OPT_FILTUPD_SAMEFILT
           else
             C1ERR = 'Invalid FilterFile.INFO OPTION: ' // copt(1:20)
             C2ERR = 'Valid options: DEFAULT SAMEFILT '
             CALL MADABORT("FILTER_UPDATE_INIT", C1ERR, C2ERR )
           endif
        endif

10    CONTINUE

    LDUMP = .false.
    if ( LDUMP ) then
        print*, ' xxxxxx FILTER_SUBDIR= ', FILTER_SUBDIR(1:20)
        print*, ' xxxxxx   PREFIX_UPD_FILTDIR =',  & 
                             PREFIX_UPD_FILTDIR(1:20)
        print*, ' xxxxxx   SUFFIX_UPD_FILTDIR =',  & 
                             SUFFIX_UPD_FILTDIR(1:20)
    endif

    RETURN
  END SUBROUTINE PARSE1_FILTER_INFO

! ======================================
    SUBROUTINE PARSE2_FILTER_INFO(NWD)

! Created July 16, 2011 by R.Kessler
! parse FILTER.INFO file that has a defined set of filters
! such as one set per SDSS-CCD column.
! IMAP is a sparse index 1-NMAP
! INDX is the absolute filter-path index in the FILTER.INFO file.
! 
! Load MAP1_FILTER_UPDATE (subdir vs. index) and
! MAP2_FILTER_UPDATE (index vs. CCID)
! 
! ------------------------


    USE SNPAR
! CDE,PARSECOM.
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

    INTEGER NWD  ! (I) number of words to check

    INTEGER  & 
         IWD, INDX, LEN  & 
        ,IMAP1, IMAP1_STORE, NMAP1, NMAP2

    character  & 
         cwd0*(MXCHAR_FILEWORD)  & 
        ,cwd1*(MXCHAR_FILEWORD)  & 
        ,cwd2*(MXCHAR_FILEWORD)  & 
        ,ccid*(MXCHAR_CCID)  & 
        ,c1err*80, c2err*80

! ---------------- BEGIN -------------

    NMAP1_FILTER_UPDATE = 0
    NMAP2_FILTER_UPDATE = 0

    NMAP1 = 0
    NMAP2 = 0

    DO 10 iwd = 1, NWD-2

        call get_PARSE_WORD_fortran(iwd+0, cwd0, LEN)
        call get_PARSE_WORD_fortran(iwd+1, cwd1, LEN)
        call get_PARSE_WORD_fortran(iwd+2, cwd2, LEN)

        if ( cwd0 .EQ. ' ' ) goto 10


        if ( cwd0 .EQ. 'FILTER_PATH:' ) then
           read(cwd1,*) INDX

           NMAP1 = NMAP1 + 1
           NMAP1_FILTER_UPDATE = NMAP1

           if ( NMAP1 .GT. MXMAP1_FILTER_UPDATE ) then
             write(c1err,1660) NMAP1, MXMAP1_FILTER_UPDATE
1660           format('NMAP=',I3,  & 
                 ' exceeds bound of MXMAP1_FILTER_UPDATE=',I3)
             c2err = 'Increase bound or reduce NMAP'
             CALL MADABORT('PARSE2_FILTER_INFO', C1ERR, C2ERR)
           endif
           MAP1_FILTER_UPDATE_INDX(NMAP1)    = INDX
           MAP1_FILTER_UPDATE_SUBDIR(NMAP1)  = CWD2
        endif

        if ( cwd0 .EQ. 'FILETRANS_SN:' ) then
          CALL SET_FILTINFO_UPD('SN', cwd1)
        endif
        if ( cwd0 .EQ. 'FILETRANS_REF:' ) then
          CALL SET_FILTINFO_UPD('REF', cwd1)
        endif



        if ( cwd0 .EQ. 'SN:' ) then
           CCID = cwd1(1:MXCHAR_CCID)
           read(cwd2,*) INDX

! translate INDX to sparse IMAP1
           IMAP1_STORE = -9
           DO imap1 = 1, NMAP1_FILTER_UPDATE
              if ( MAP1_FILTER_UPDATE_INDX(imap1) .EQ. INDX) then
                 imap1_store = imap1
              endif
           END DO

           if ( IMAP1_STORE .LE. 0 ) then
             write(c1err,670) INDX, CCID
             c2err = 'Check FILTER.INFO file.'
670            format('Undefined INDX=',I8,'   for  CID = ',A8)
             CALL MADABORT('PARSE2_FILTER_INFO', C1ERR, C2ERR)
           endif

! update 2nd map
           NMAP2 = NMAP2 + 1
           NMAP2_FILTER_UPDATE = NMAP2

           if ( NMAP2 .GT. MXSNLC_FILTUPD ) then
              write(c1err,661) MXSNLC_FILTUPD
661             format('MAP2 size exceeds bound of MXSNLC_FILTUPD=',I6)
              c2err = 'Check FILTER.INFO file.'
              CALL MADABORT("PARSE2_FILTER_INFO", C1ERR, C2ERR)
           endif

           MAP2_FILTER_UPDATE_CCID(NMAP2)    = CCID
           MAP2_FILTER_UPDATE_PTRMAP1(NMAP2) = IMAP1_STORE

        endif  !   'SN:' key


10    CONTINUE

! ---------------------------------------------------------

    write(6,60) NMAP1_FILTER_UPDATE
60    format(T4,'Found ',I2,' filter sets to use for updates.')

    write(6,61) NMAP2_FILTER_UPDATE
61    format(T4,'Filter update map stored for ',I6,' SNe.' )

    print*,' -------------------------------------------------- '

    RETURN
  END SUBROUTINE PARSE2_FILTER_INFO

! =================================
    SUBROUTINE SET_FILTINFO_UPD(KEY,STRING)
! 
! Created Jul 16, 2011
! 
! Parse STRING for asterisk, and store PREFIX and SUFFIX
! before/after the asterisk. Used later to determine the
! name of the filter transmission  file where the  asterisk
! is replaced with the filter-character.
! 
! Jul 5 2013: string-length fixes to avoid segfault in debug mode.
!             See use of LKEY and LSTR.
! 
! ---------


    USE SNPAR
    USE FILTCOM
    USE FILTUPDCM

    IMPLICIT NONE

! input arguments
    CHARACTER KEY*(*)     ! (I) 'SN' or 'REF'
    CHARACTER STRING*(*)  ! (I) string to define filter filename

! local variables

    INTEGER istar, LSTR, LKEY
    LOGICAL LDUMP
! ------------ BEGIN -------------

    LSTR = INDEX(STRING,' ' ) - 1
    LKEY = INDEX(KEY//' ',   ' ' ) - 1

    istar = Index(STRING,'*')

    IF ( KEY(1:LKEY) .EQ. 'REF' ) THEN
        PREFIX_UPD_TRANSREF = STRING(1:istar-1)
        SUFFIX_UPD_TRANSREF = STRING(istar+1:LSTR)
        FILTINFO_UPD_REF = .TRUE.

    ELSE IF ( KEY(1:LKEY) .EQ. 'SN' ) THEN
        PREFIX_UPD_TRANSSN = STRING(1:istar-1)
        SUFFIX_UPD_TRANSSN = STRING(istar+1:LSTR)
        FILTINFO_UPD_SN = .TRUE.

    ELSE

    ENDIF

! print comment about filename syntax

    LDUMP = .TRUE.
    IF ( LDUMP ) THEN
       write(6,20) KEY, STRING(1:LSTR)
20       format(T4,'Filter-Trans filenames(', A, ') : ''',A,  & 
                '''  with * => filter' )
    ENDIF

    RETURN
  END SUBROUTINE SET_FILTINFO_UPD



! ==============================

! utility to fit for the peak-MJD.
! ==============================

#if defined(MINUIT)

! ====================================
    SUBROUTINE SET_PEAKMJD()

! Created Oct 24, 2009 by R.Kessler
! 
! Analyze and/or fit light curve to estimate peak-MJD.
! Fill SNLC_SEARCH_PEAKMJD to be used as the initial
! Tpeak-estimate in the light curve fit.
! 
! This routine will not work properly for low SNR light curves.
! The minimum SNR on a light curve must be at least 4.
! 
! Namelist variable OPT_SETPKMJD determines options;
! see OPTBIT_SETPKMJD_XXX options in global declaration.
! 
!           HISTORY
! 
! May 10 2016: remove PSF>2 epochs to get fluxMax (avoid crazy fluxes)
!              See new parameter PSFMAX_forFLUXMAX
! 
! May 17, 2019: implement new options OPTBIT_SETPKMJD_[TRIGGER,MAXFLUX2]
! May 20, 2019: return if less then 2 epochs.
! Jan 27, 2020: use STDOUT_UPDATE flag to print screen update.
! May 01, 2020: implement SHIFT_SETPKMJD
! May 05, 2020: implement OPT_SETPKMJD_SIM option
! -------------------------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE SNANAFIT
    USE PKMJDCOM

    IMPLICIT NONE

    INTEGER  & 
         EPMAX, EP  & 
        ,EP_atFLUXMAX(0:MXFILT_ALL), EPTEST_atFLUXMAX(0:MXFILT_ALL)  & 
        ,IFILT_OBS, ifilt, ifilt_snrmin  & 
        ,IFILT_OBS2, IFILT_OUTLIER, IFILT_REJECT  & 
        ,NPKMJD, LUN, VBOSE

    REAL  & 
         PKMJDSUM, PKMJDAVG, PKMJDDIF, PKMJDDIF_MAX  & 
        ,PKMJD, DIF, WGTSUM, WGT  & 
        ,FLUX, FLUXERR, SNR, SNRMIN  & 
        ,FLUXMAX(0:MXFILT_ALL)  & 
        ,MJDatFLUXMAX(0:MXFILT_ALL)  & 
        ,SNRatFLUXMAX(0:MXFILT_ALL)

    REAL*8  & 
         MJD8, F8, PKMJD8, PKMJD8ERR, PKFLUX8, PKFLUX8ERR  & 
        ,PARLIST(3), FITWIN_MJD8(2)

    LOGICAL  & 
          LSET, LDMP, LDMP_USER, NOABORT, LABORT  & 
         ,USEFILT(MXFILT_ALL)  & 
         ,USE_MJDatFLUXMAX, USE_MJDatFMAXCLUMP  & 
         ,USE_MJDatTRIGGER  & 
         ,USE_SIMPKMJD  & 
         ,FIRST_INIT, DOFIT, NOFIT

    CHARACTER cfilt*8, FNAM*12, METHOD*24, CCID*(MXCHAR_CCID)

! define hard-wired cut-parameters
    REAL, PARAMETER ::                 & 
         MJDMIN_LC           = 40000.  &  ! min MJD to consider for light curve
        ,PKMJD_OUTLIER       = 30.0    &  ! reject PKMJD off by this much
        ,MAX_SNRMIN_REJECT   = 2.0        ! max SNRMIN to reject filter

! function
    EXTERNAL init_obs_atFLUXMAX, get_obs_atFLUXMAX

! ---------------- BEGIN ----------------

    if ( OPT_SETPKMJD <= 0 ) RETURN

    FNAM = 'SET_PEAKMJD'
    CCID = SNLC_CCID

    if ( ISNLC_NEWMJD_STORE < 2 ) RETURN

! July 12 2013: check option to use MJD at max flux (no fitting)
    USE_MJDatFLUXMAX  = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_FLUXMAX)
    USE_MJDatFMAXCLUMP =  & 
                 BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_FLUXMAX2)  & 
           .or.  BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_FLUXMAX3)
    USE_MJDatTRIGGER  = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_TRIGGER)
    USE_SIMPKMJD      = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_SIM)
    DOFIT = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_ANYFUN) .or.  &  ! Bazin09
              BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_POLYCOR)
    NOFIT = (.NOT. DOFIT)

    if ( USE_MJDatTRIGGER .and. PHOTFLAG_TRIGGER == 0 ) then
      c1err = 'Must set PHOTFLAG_TRIGGER in &SNLCINP'
      c2err = 'to set PEAKMJD = MJD_TRIGGER'
      CALL MADABORT(FNAM, c1err, c2err)
    endif

! May 5 2020: check option to use true SIM_PEAKMJD
    IF ( USE_SIMPKMJD ) THEN
      if ( LSIM_SNANA ) then
          PKMJD = SIM_PEAKMJD
          GOTO 800
      else
         c1err = 'Cannot use OPT_SETPKMJD += 2048 for real data'
         c2err = 'This option (use SIM_PEAKMJD) is for sim only'
         CALL MADABORT(FNAM, c1err, c2err)
      endif
    ENDIF

! bail if the user CUTWIN_PEAKMJD does not contain  any observations:
    IF ( SNLC8_MJDMIN .GT. CUTWIN_PEAKMJD(2) ) RETURN
    IF ( SNLC8_MJDMAX .LT. CUTWIN_PEAKMJD(1) ) RETURN

! first-time init:
!  + print comments to stdout
!  + open separate log file to dump MNFIT stuff for each SN
    FIRST_INIT = (NFIT_PKMJD==0)

    IF  ( FIRST_INIT ) THEN  ! init
       LSET = .FALSE.
       print*,' '
       if ( USE_MJDatFLUXMAX ) then
          write(6,80) 'MJD at max flux' ; LSET=.TRUE.
          METHOD = '(FLuxMax)'
       else if ( USE_MJDatFMAXCLUMP ) then
          write(6,80)  & 
               'MJD at max flux in dense clump of SNR-detections.'
          LSET=.TRUE.
          METHOD = '(FmaxClump)'
       else if ( USE_MJDatTRIGGER ) then
          write(6,80) 'MJD_TRIGGER' ; LSET=.TRUE.
          METHOD = '(MJD_trigger)'
       endif

       if( DOFIT ) then
          write(6,80) 'Bazin fit' ; LSET=.TRUE.
          METHOD = '(FluxMax+BazinFit)'
          if ( USE_MJDatFMAXCLUMP ) METHOD = '(FmaxClump+BazinFit)'
       endif
 80      format(' SET_PEAKMJD method: ', A)
       print*,' '
       call flush(6)

       if ( .not. LSET ) then
          write(6,80) 'UNKNOWN'
          write(c1err,86) OPT_SETPKMJD
 86         format('Unknown OPT_SETPKMJD=',I5,' in &SNLCINP.')
          c2err = 'See manual for OPT_SETPKMJD'
          CALL MADABORT(FNAM, c1err,c2err )
       endif

       VBOSE = 1
       PARLIST(1) = MJDWIN_SETPKMJD    ! user MJDWIN from &SNLCINP
       PARLIST(2) = SNRCUT_SETPKMJD    ! user SNRCUT from &SNLCINP
       PARLIST(3) = SNRMIN_forFLUXMAX  ! backup SNRCUT if above always fails
       call init_obs_atFLUXMAX(OPT_SETPKMJD, PARLIST, VBOSE)

! open outfile to dump MINUIT stdout
       if ( DOFIT ) then
          LUN = LUNPKMJD
          OPEN(UNIT=LUN,FILE=MNFIT_PKMJD_LOGFILE,status='UNKNOWN')
       endif

! open INTERP-OUT file if set
       IF ( SNMJD_OUT_FILE .NE. ' ' ) THEN
          LUN = LUNINTERP
          OPEN(UNIT=LUN,FILE=SNMJD_OUT_FILE,status='UNKNOWN')
          write(LUNINTERP,46)
 46         format('VARNAMES:',2x,'CID', 4x,'FILT', 5x,'MJD',  & 
                 10x,'FLUX      FLUXERR')
       ENDIF

    ENDIF  ! end of FIRST_INIT

    NFIT_PKMJD = NFIT_PKMJD  + 1

    global_banner =  & 
         'SET_PEAKMJD: Find Estimate of PEAK-MJD for CID='  //  & 
         CCID
!       CALL PRBANNER ( global_banner )

    IF ( NFILTDEF_SURVEY .GT. MXFILT_OBS ) THEN
      write(c1err,1661) NFILTDEF_SURVEY
      write(c2err,1662) MXFILT_OBS
1661    format('There are ', I3,' obs-filters,')
1662    format('but SET_PEAKMJD works for up to ',I2,' filters.')
      CALL MADABORT(FNAM, c1err, c2err)
    ENDIF

! initialize

    DO 480 IFILT_OBS = 0, MXFILT_ALL
       EP_atFLUXMAX(IFILT_OBS)  = -9
       FLUXMAX(IFILT_OBS)       =  0.0
       MJDatFLUXMAX(IFILT_OBS)  = -9.0
       SNRatFLUXMAX(IFILT_OBS)  = -9.0

       IF ( IFILT_OBS .EQ. 0 ) goto 480

       PKMJD_FIT(IFILT_OBS)     =  0.0
       PKMJD_ERR(IFILT_OBS)     = -9.0
       PKFLUX_FIT(IFILT_OBS)    = -9.0
       PKFLUX_ERR(IFILT_OBS)    = -9.0

       USEFILT(ifilt_obs) = .FALSE.
 480  CONTINUE


    PKMJD_ERRMIN = 1.0E9
    PKMJD_ERRWGT = 0.0

    PKFLUX_ERRMIN = 1.0E9
    PKFLUX_ERRWGT = 0.0

! find epoch with max flux for each filter using C-function

     CALL get_obs_atFLUXMAX(  & 
            CCID(1:ISNLC_LENCCID)//char(0),  & 
            ISNLC_NEWMJD_STORE,  & 
            SNLC_FLUXCAL, SNLC_FLUXCAL_ERRTOT,  & 
            SNLC8_MJD, ISNLC_IFILT_OBS,  & 
            EP_atFLUXMAX(0),             &  ! <== returned array vs. IFILTOBS
            MXCHAR_CCID )

! use EP at max flux to set FLUXMAX, MJD at max, SNR at max.
    DO 40 ifilt = 0, NFILTDEF_SURVEY
       if ( ifilt == 0 ) then
          ifilt_obs = 0
       else
          ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       endif
       EP = EP_atFLUXMAX(ifilt_obs)+1 ! from C code

       if ( EP < 0 ) then
         SNRatFLUXMAX(ifilt_obs) = 0.01
         MJDatFLUXMAX(ifilt_obs) = 0.01
         goto 40
       endif

       FLUX      = SNLC_FLUXCAL(EP)
       FLUXERR   = SNLC_FLUXCAL_ERRTOT(EP)

       FLUXMAX(ifilt_obs)      = FLUX
       MJDatFLUXMAX(ifilt_obs) = SNGL(SNLC8_MJD(ep))
       if ( FLUXERR > 0.0 ) then
          SNRatFLUXMAX(ifilt_obs) = FLUX/FLUXERR
       else
          SNRatFLUXMAX(ifilt_obs) = 0.0
       endif
40     CONTINUE

! -------------------------------------
! check option to use MJD at max-flux or at TRIGGER
    IF ( NOFIT ) THEN
       IF ( USE_MJDatFLUXMAX .or. USE_MJDatFMAXCLUMP ) THEN
          PKMJD  = MJDatFLUXMAX(0)
       ELSE IF ( USE_MJDatTRIGGER ) THEN
          PKMJD  = SNGL(SNLC8_MJD_TRIGGER)
       ELSE
          C1err = 'Invalid NOFIT option.'
          c2err = 'Expecting MJD-at-maxFlux or MJD-at-Trigger'
          CALL MADABORT(FNAM, c1err, c2err)
       ENDIF

       SNLC_SEARCH_PEAKMJD    = PKMJD
       IF ( STDOUT_UPDATE ) then
          write(6,50) CCID(1:ISNLC_LENCCID), PKMJD, METHOD
          call FLUSH(6)
       ENDIF

       RETURN
    ENDIF


! --------- Below do some kind of Bazin fit -----------

! --------------
! Find valid filter with smallest SNR at epoch with max-FLux,
! and that has SNR below MAX_SNRMIN_REJECT.
! If all filters have SNR > MAX_SNRMIN_REJECT then use all filters;
! otherwise toss the one with smallest SNR at FLUXMAX.

    SNRMIN = 999999999.
    IFILT_SNRMIN = -9
    DO ifilt = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       SNR = SNRatFLUXMAX(ifilt_obs)

       if (    SNR .LT. SNRMIN  & 
           .and. SNR .LT. MAX_SNRMIN_REJECT  & 
           .and. SNR .GT. 0.0  ) then

          SNRMIN = SNR
          ifilt_snrmin = ifilt_obs  ! filter to reject
       endif
    ENDDO

! -----------------------------------------
! Now fit each filter to general function using MINUIT.
! Take weighted PKMJD-average so that badly measured
! PKMJD don't have much effect

    NPKMJD   = 0
    PKMJDSUM = 0.0
    WGTSUM   = 0.0
    PKMJDDIF_MAX = 0.0
    FITWIN_MJD8(1) = MJDatFLUXMAX(0) + CUTWIN_TOBS(1)
    FITWIN_MJD8(2) = MJDatFLUXMAX(0) + CUTWIN_TOBS(2)

    DO 200 ifilt = 1, NFILTDEF_SURVEY
       ifilt_obs = IFILTDEF_MAP_SURVEY(ifilt)
       MJD8      = DBLE( MJDatFLUXMAX(ifilt_obs) )
       F8        = DBLE ( FLUXMAX(IFILT_OBS) )

       if ( MJD8      .LT. MJDMIN_LC    ) GOTO 200
       if ( ifilt_obs .EQ. ifilt_snrmin ) GOTO 200

       USEFILT(ifilt_obs) = .TRUE.

       CALL MNFIT_PKMJD(IFILT_OBS, FITWIN_MJD8, F8, MJD8,    &  ! inputs
                    PKMJD8, PKMJD8ERR, PKFLUX8, PKFLUX8ERR)   ! return args

! store result for each filter.

         PKFLUX_FIT(ifilt_obs) = SNGL(PKFLUX8)
         PKFLUX_ERR(ifilt_obs) = SNGL(PKFLUX8ERR)

         PKMJD_FIT(ifilt_obs) = SNGL(PKMJD8)
         PKMJD_ERR(ifilt_obs) = SNGL(PKMJD8ERR)
         NPKMJD = NPKMJD + 1

! check option to interpolate flux
         CALL SNMJD_INTERP_ANYLC(ifilt_obs)

200   CONTINUE

    IF ( NPKMJD .EQ. 0 ) then

       NOABORT = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_NOABORT)  & 
              .or. (OPT_REFORMAT_FITS > 0)  & 
              .or. (.not. ABORT_ON_NOPKMJD)  ! 7.29.2014

       LABORT  = (.NOT. NOABORT)

       IF ( LABORT ) then
          C1err = 'NMJD=0 => no estimate of PKMJD !!!'
          c2err = 'Check SN ' // CCID
          CALL MADABORT(FNAM, c1err, c2err)
       ELSE
          write(6,240) CCID
 240        format(T8,'Could not estimate PKMJD for CID=',A)
          call flush(6)
          PKMJD = -1.
          GOTO 800
       ENDIF
    ENDIF

! convert sum of 1/err^2 into wgted error.
    IF  ( PKMJD_ERRWGT .GT. 0.0 ) THEN
       PKMJD_ERRWGT  = SQRT(1.0/PKMJD_ERRWGT)
    ENDIF
    IF ( PKFLUX_ERRWGT .GT. 0.0 ) THEN
       PKFLUX_ERRWGT = SQRT(1.0/PKFLUX_ERRWGT)
    ENDIF

! if there are more than 2 filters, check for
! crazy outlier.

    IFILT_REJECT    = -9
    IFILT_OUTLIER   = -9

    IF ( NPKMJD .GE. 3 ) THEN
    DO 252 ifilt_obs2 = 1, MXFILT_ALL  ! ignore this filter in average
       if ( .not. USEFILT(ifilt_obs2) ) goto 252
       PKMJDSUM = 0.0
       WGTSUM = 0.0

    DO 250 ifilt_obs  = 1, MXFILT_ALL
       if ( .not. USEFILT(ifilt_obs) ) goto 250
       if ( ifilt_obs2 .NE. ifilt_obs ) then
         PKMJDSUM = PKMJDSUM + (PKMJD_FIT(ifilt_obs)-SNGL(MJDOFF))
         WGTSUM   = WGTSUM + 1.0
       endif
250   CONTINUE

! test average against filter that was excluded

       PKMJDAVG = (PKMJDSUM/WGTSUM + SNGL(MJDOFF) )
       PKMJDDIF = abs(PKMJDAVG - PKMJD_FIT(ifilt_obs2))

       if ( PKMJDDIF .GT. PKMJDDIF_MAX ) then
          PKMJDDIF_MAX  = PKMJDDIF
          IFILT_OUTLIER = ifilt_obs2
       endif

252   CONTINUE

       if ( PKMJDDIF_MAX .GT. PKMJD_OUTLIER ) then
         USEFILT(ifilt_outlier) = .FALSE.
         IFILT_REJECT = IFILT_OUTLIER

!            print*,' xxx REJECT FILT=',IFILT_OUTLIER,
!     &       '  for CID=',CCID,'  PKMJDDIF=', PKMJDDIF_MAX
       endif

    ENDIF

! ---------------------------------------
! Finally, take PKMJD-wgted average of remaining filters.

    PKMJDSUM = 0.0
    WGTSUM   = 0.0
    NPKMJD   = 0

    DO 300 ifilt_obs = 1, MXFILT_ALL
       if ( .not. USEFILT(ifilt_obs) ) goto 300

         PKMJD8    = PKMJD_FIT(ifilt_obs)
         PKMJD8ERR = PKMJD_ERR(ifilt_obs)

         WGT    = SNGL(1.0/PKMJD8ERR) ! * 10.0/MAX(PKMJDDIF,10.0)

         PKMJDDIF = SNGL(PKMJD8 - MJDOFF)

         PKMJDSUM = PKMJDSUM + WGT*PKMJDDIF
         WGTSUM = WGTSUM + WGT
         NPKMJD = NPKMJD + 1

300   CONTINUE

    PKMJDAVG = PKMJDSUM/WGTSUM + SNGL(MJDOFF)
    PKMJD    = PKMJDAVG

    if ( STDOUT_UPDATE ) THEN
      write(6,50) CCID(1:ISNLC_LENCCID), PKMJD, METHOD
50      format(T8,'SET_PEAKMJD(',A,') --> ', F9.3, 3x, A )
      call flush(6)
    endif

! -----------------------------------------------

    LDMP_USER =  BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_DUMP)
    LDMP      = LDMP_USER

    IF ( LDMP ) THEN
      print*,'  --------------------------------------------- '
      print*,'  BEGIN DUMP of SET_PEAKMJD for SN ', CCID
      print*,'  '
      print*,'                  MJDat    SNRat      PKMJD  PKMJD  '
      print*,'  filt  FLUXMAX   FLUXMAX  FLUXMAX    fit    err    '
      print*,'  --------------------------------------------------'

      DO 660 ifilt_obs = 1, MXFILT_ALL
         if ( .not. USEFILT(ifilt_obs) ) goto 660
         cfilt =  filtdef_string(ifilt_obs:ifilt_obs)
         write(6,661) cfilt  & 
             , FLUXMAX(ifilt_obs)  & 
             , MJDatFLUXMAX(ifilt_obs)  & 
             , SNRatFLUXMAX(ifilt_obs)  & 
             , PKMJD_FIT(ifilt_obs)  & 
             , PKMJD_ERR(ifilt_obs)
661        format( T5,A, T7,F10.1, F11.2, F6.1, F11.2, F7.2 )

660      CONTINUE
       print*, ' '

      if ( ifilt_snrmin .GT. 0 ) then
        cfilt =  filtdef_string(ifilt_snrmin:ifilt_snrmin)
      else
        cfilt = 'none'
      endif
      print*,'    Filter with min SNR at max: ', cfilt

      if ( ifilt_reject .GT. 0 ) then
        cfilt =  filtdef_string(ifilt_reject:ifilt_reject)
      else
        cfilt = 'none'
      endif
      print*,'    PKMJD-outlier filter rejected : ', cfilt

      print*,'    PEAKMJD(CALC,SEARCH)  = ',  & 
                  PKMJD, SNLC_SEARCH_PEAKMJD

      print*,'    CALC-SEARCH PKMJD DIF = ', DIF
      print*,'  END DUMP of SET_PEAKMJD for SN ', CCID
      print*,'  --------------------------------------------- '
    ENDIF

! -----------------------------------------
800   CONTINUE

    PKMJD = PKMJD + SHIFT_SETPKMJD
    SNLC_SNANAFIT_PEAKMJD = PKMJD
    SNLC_SEARCH_PEAKMJD   = PKMJD

    RETURN
  END SUBROUTINE SET_PEAKMJD


! =========================================
    SUBROUTINE MNFIT_PKMJD(ifilt_obs, FITWIN_MJD,  & 
              FLUXatMAX, MJDatMAX, PKMJD,PKMJDERR, PKFLUX, PKFLUXERR )

! -----------------------------------------------
! Created Oct 29, 2009 by R.Kessler
! Fit FCNANYLC (any LC function) to estimate PKMJD.
! 
! METHOD:
! First get  approx PEAKMJD = filter-average of MJDs at max flux.
! Then grid-minimize Eq. 1 of  Bazin et al., (SNLC CC rate)
! from astro-ph 0904.1066.
! 
!               exp[ -(t-t0)/Tfall ]
!  f(t) = A * --------------------------  +  B
!              1 + exp[ -(t-t0)/Trise ]
! 
! where t0, Trise, Tfall, A, & B=0 are the five
! fit-parameters for each filter. We only care about
! getting an approximate function, so we don't care
! about errors.
! 
! Note that the exponent in their denomintor is missing a
! minus sign (with Trise > 0). The t0 parameter is not really
! the time at peak; after the fit has converged, the true Tpeak
! is found by setting the derivative=0.
! 
! Jun 6, 2010: for POLYCOR option, do BTEST on OPT_SETPKMJD
! 
! Oct 26, 2010: remove "SET PRint" command so that we get
!               correlation matrix printed to screen.
! 
! Dec 9, 2011: keep track of min-error and  wgted error for
!              T0 and F0err/F0.
! 
! May 28 2019: add FITWIN_MJD argument to restrict range of epochs fit.
! Sep 22 2020: MNARG -> MNARG(10) to avoid warnings with gcc v10
! ----------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE PKMJDCOM

    IMPLICIT NONE

! function args

    INTEGER IFILT_OBS   ! (I) SN index and obs filter to fit

    REAL*8  & 
         FITWIN_MJD(2)  &  ! (I) MJD window to fit observations
        ,MJDatMAX       &  ! (I) MJD at epoch with max flux
        ,FLUXatMAX      &  ! (I) flux at MJD with max flux
        ,PKMJD       &  ! (O) MJD at peak of function
        ,PKMJDERR    &  ! (O) approx error on PKMJD (used for wgted avg)
        ,PKFLUX      &  ! (O) peak flux from fit
        ,PKFLUXERR  ! (O) error on above

! MINUIT stuff

    INTEGER IERR, IPAR, IVARBL

    INTEGER  & 
         IUNIN, IUOUT, IUWRI  & 
        ,NFIXPAR, IFILT

    EXTERNAL FCNANYLC, ANYLCFUN
    DOUBLE  PRECISION  ANYLCFUN

    LOGICAL USE_ANYFUN, USE_POLYCOR

    REAL*8  & 
         MNARG(10)  & 
        ,FIXLIST(NPAR_ANYLC)  & 
        ,INIBND(2,NPAR_ANYLC)  & 
        ,INISTP(NPAR_ANYLC)  & 
        ,INIVAL(NPAR_ANYLC)  & 
        ,TMPVAL(10), TMPERR(10), BND1, BND2  & 
        ,T0, TRISE, TFALL, ARG, T0ERR  & 
        ,INIVAL_TMP, COV  & 
        ,ERRTMP

! functions
    REAL*8 COVPKFLUXFUN

! ---------------- BEGIN -------------

    USE_ANYFUN  = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_ANYFUN)
    USE_POLYCOR = BTEST(OPT_SETPKMJD,OPTBIT_SETPKMJD_POLYCOR)

    IUNIN = 5
    IUOUT = LUNPKMJD   ! write to file instead of STDOUT=6
    IUWRI = 7

    IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)

! ---------------------------------
! write header to log file so we know what CID and filter

    write(LUNPKMJD,20)  & 
         SNLC_CCID  & 
        ,filtdef_string(ifilt_obs:ifilt_obs)
20    format(/, T2, 60('#'),  & 
             /, T10,'MNFIT_PKMJD for CID=',A8,' FILTER=',A,  & 
             /, T2, 60('#'), / )

! ------------------------

    CALL MNINIT(IUNIN,IUOUT,IUWRI)  ! init MINUIT

    MNARG(1) = dble(0.0)
    CALL MNEXCM(FCNANYLC, 'SET NOWarnings', MNARG,0,IERR, ANYLCFUN)

    MNARG(1) = DBLE(MINUIT_PRINT_LEVEL)
    CALL MNEXCM(FCNANYLC, 'SET PRI', MNARG, 1, IERR, ANYLCFUN )



! init each parameter
    INIBND(1,IPAR_ISN) = 0
    INIBND(2,IPAR_ISN) = 9999999.
    INIVAL(IPAR_ISN)   = DBLE(SNLC_CID)
    INISTP(IPAR_ISN)   = 0.0       ! fixed parameter

    INIBND(1,IPAR_FILT) = 0
    INIBND(2,IPAR_FILT) = 100
    INIVAL(IPAR_FILT)   = DBLE(IFILT_OBS)
    INISTP(IPAR_FILT)   = 0.0     ! fixed parameter

    INIBND(1,IPAR_MJDMIN) = 0
    INIBND(2,IPAR_MJDMIN) = 150000
    INIVAL(IPAR_MJDMIN)   = FITWIN_MJD(1)
    INISTP(IPAR_MJDMIN)   = 0.0     ! fixed parameter

    INIBND(1,IPAR_MJDMAX) = 0
    INIBND(2,IPAR_MJDMAX) = 150000
    INIVAL(IPAR_MJDMAX)   = FITWIN_MJD(2)
    INISTP(IPAR_MJDMAX)   = 0.0     ! fixed parameter

    INIBND(1,IPAR_T0) = MJDatMAX - 50.0
    INIBND(2,IPAR_T0) = MJDatMAX + 50.0
    INIVAL(IPAR_T0)   = MJDatMAX
    INISTP(IPAR_T0)   = 2.0

    INIBND(1,IPAR_TRISE) =  1.0
    INIBND(2,IPAR_TRISE) =  40.
    INIVAL(IPAR_TRISE)   =  5.0
    INISTP(IPAR_TRISE)   =  2.0

    INIBND(1,IPAR_TFALL) = 1.0
    INIBND(2,IPAR_TFALL) = 100.0
    INIVAL(IPAR_TFALL)   = 10.0
    INISTP(IPAR_TFALL)   = 2.0

    INIBND(1,IPAR_A0) = FLUXatMAX * 0.2
    INIBND(2,IPAR_A0) = FLUXatMAX * 5.0
    INIVAL(IPAR_A0)   = FLUXatMAX * 1.2
    INISTP(IPAR_A0)   = FLUXatMAX * 0.1

    IF ( USE_POLYCOR ) THEN
       INIVAL_TMP = 0.1
    ELSE
       INIVAL_TMP = 0.0
    ENDIF

    INIBND(1,IPAR_A1)   = -20.0
    INIBND(2,IPAR_A1)   =  20.0
    INIVAL(IPAR_A1)     =   0.0
    INISTP(IPAR_A1)     = INIVAL_TMP

    INIBND(1,IPAR_A2)   = -20.0
    INIBND(2,IPAR_A2)   =  20.0
    INIVAL(IPAR_A2)     =   0.0
    INISTP(IPAR_A2)     = INIVAL_TMP
!      PKPARNAME(IPAR_A2)  =  'A2'


    CHI2_FITPKMJD(IFILT) = -9.0
    NDOF_FITPKMJD(IFILT) = -9

    NFIXPAR = 0
    DO IPAR = 1, NPAR_ANYLC

       CALL MNPARM(IPAR,    PKPARNAME(IPAR)  & 
             ,INIVAL(IPAR),   INISTP(IPAR)  & 
             ,INIBND(1,IPAR), INIBND(2,IPAR)  & 
             ,IERR   )

       IF ( IERR .NE. 0 ) THEN
          WRITE (6,'(A,A10,A)') '  ERROR initializing ',  & 
                  PKPARNAME(ipar), '  with MNPARM '
          RETURN
       ENDIF

! keep track of fixed parameters
      if ( INISTP(IPAR) .EQ. 0.0 ) then
         NFIXPAR =  NFIXPAR + 1
         FIXLIST(NFIXPAR) = DBLE ( IPAR )
      endif

    ENDDO  ! end of IPAR loop

    CALL FLUSH(6)

    CALL MNEXCM(FCNANYLC, 'FIX', FIXLIST, NFIXPAR, IERR, ANYLCFUN )

    MNARG(1) = dble(30000.0)  ! MAXCALLS
    CALL MNEXCM(FCNANYLC, 'MINIMIZE', MNARG, 0, IERR, ANYLCFUN )

    DO IPAR = 1,  NPAR_ANYLC

      CALL MNPOUT(  IPAR, PKPARNAME(IPAR)  & 
                    , TMPVAL(IPAR), TMPERR(IPAR)  & 
                    , BND1, BND2, IVARBL )

      SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar) = SNGL(TMPVAL(IPAR))
      SNLC_SNANAFIT_PEAKMJD_FITERR(ifilt,ipar) = SNGL(TMPERR(IPAR))

    ENDDO

    MNARG(1) = ZERO8
    CALL MNEXCM ( FCNANYLC, 'EXIT',  MNARG, 0, IERR, ANYLCFUN )
    PRINT *,'  MNFIT_PKMJD: EXIT returns IERR = ', IERR

! compute PKMJD from setting 1st deriviate = 0
    Tfall = TMPVAL(IPAR_TFALL)
    Trise = TMPVAL(IPAR_TRISE)
    T0    = TMPVAL(IPAR_T0)
    T0ERR = TMPERR(IPAR_T0)

    arg = (Tfall - Trise)/Trise
    if ( ARG .GT. 0.1 ) then
       PKMJD = T0 + Trise * DLOG(ARG)
       PKMJDERR = MIN(T0ERR,4.0)  ! max error is 4 days
    else
       PKMJD = MJDatMAX
       PKMJDERR = 4.0
    endif


    if ( OPT_REFORMAT_SALT2 .GT. 0 ) then
       PKMJDERR = T0ERR
    endif

! get error on peak flux from fit

    CALL MNEMAT(FITERRMAT_PKMJD(1,1,ifilt_obs),NPAR_ANYLC)

    COV = COVPKFLUXFUN(ifilt_obs,PKMJD,PKMJD)

    PKFLUX    = ANYLCFUN(PKMJD,TMPVAL)
    PKFLUXERR = SQRT(COV)

! Dec 9, 2011
! keep track of min error and wgted error
    ERRTMP  = T0ERR
    IF ( ERRTMP .LT. PKMJD_ERRMIN .and. ERRTMP .GT. 0.0 ) THEN
       PKMJD_ERRMIN = SNGL(ERRTMP)
       PKMJD_ERRWGT = PKMJD_ERRWGT + 1.0/SNGL(ERRTMP*ERRTMP)
    ENDIF

    ERRTMP  = PKFLUXERR/PKFLUX
    IF ( ERRTMP .LT. PKFLUX_ERRMIN .and.  & 
           ERRTMP .GT. 0.0   .and.  & 
           ERRTMP .NE. 1.0 ) THEN
      PKFLUX_ERRMIN = SNGL(ERRTMP)
      PKFLUX_ERRWGT = PKFLUX_ERRWGT + 1.0/SNGL(ERRTMP*ERRTMP)
    ENDIF

    RETURN
  END SUBROUTINE MNFIT_PKMJD

! =========================================
    SUBROUTINE FCNANYLC(NVAR,GRAD,CHI2,XVAL,IFLAG, ANYLCFUN )
! 
! Created Oct 25, 2009
! MINUIT-callable function to return data-model CHI2
! using the "ANY LC" model from Bazin et al.
! See comments at top of MNFIT_PKMJD.
! 
! May 28, 2019: Use MJDMIN/MJDMAX to restrict epoch range.
! 
! --------------------------------------------


    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE PKMJDCOM

    IMPLICIT NONE

    INTEGER NVAR, IFLAG
    REAL*8  & 
         XVAL(*)   &  ! (I) T0, Trise, Tfall, A, B
        ,GRAD(*)  & 
        ,CHI2  & 
        ,ANYLCFUN

    EXTERNAL ANYLCFUN

! local var

    INTEGER  & 
         ISN, IFILT_OBS, IFILT  & 
        ,NEWMJD, EPMIN, EPMAX, EP, NDOF

    REAL*8  & 
         MJDMIN, MJDMAX, FLUX, FLUXERR, MJD, SQDIF, SQERR, SNR  & 
        ,MODEL_FLUX, MODEL_FLUXERR

    REAL*8, PARAMETER :: MODEL_MAGERR = 0.05

! --------------- BEGIN -------------

    ISN         = INT ( XVAL(1) + 0.0001 )
    IFILT_OBS   = INT ( XVAL(2) + 0.0001 )
    MJDMIN      = XVAL(3)
    MJDMAX      = XVAL(4)

    CHI2 = 0.0
    NDOF = 0

    DO 200 NEWMJD = 1, ISNLC_NEWMJD_STORE

      EPMIN = ISNLC_EPOCH_RANGE_NEWMJD(1,NEWMJD)
      EPMAX = ISNLC_EPOCH_RANGE_NEWMJD(2,NEWMJD)
      MJD   = SNLC8_MJD(EPMIN)

      if ( MJD < MJDMIN ) goto 200
      if ( MJD > MJDMAX ) goto 200

      MODEL_FLUX    = ANYLCFUN(MJD,XVAL)
      MODEL_FLUXERR = MODEL_FLUX * MODEL_MAGERR

      DO 201 EP = EPMIN, EPMAX

          if ( IFILT_OBS .NE. ISNLC_IFILT_OBS(ep) ) GOTO 201

          FLUX      = SNLC_FLUXCAL(ep)
          FLUXERR   = SNLC_FLUXCAL_ERRTOT(ep)

! SNR requirement makes things worse ??
          SNR = FLUX/FLUXERR

          SQDIF = (FLUX - MODEL_FLUX)**2
          SQERR = FLUXERR*FLUXERR + MODEL_FLUXERR*MODEL_FLUXERR
          CHI2  = CHI2 + SQDIF/SQERR
          NDOF = NDOF + 1

201     CONTINUE   ! ep
200   CONTINUE  ! end of NEWMJD

! --------------------------------------
! load chi2 and NDOF in global array

    IF (IFLAG == FCNFLAG_LAST ) THEN
      IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
      CHI2_FITPKMJD(IFILT) = SNGL(CHI2)
      NDOF_FITPKMJD(IFILT) = NDOF
    ENDIF

    RETURN
  END SUBROUTINE FCNANYLC

! ======================================
    DOUBLE PRECISION FUNCTION ANYLCFUN(MJD,XVAL)
! 
! Created Dec 5, 2009 by R.Kessler
! Any light curve function based on Bazin et al,, 2009 (SNLS CC rate)
! 
! May 28 2019: replace hard-code XVAL-indices with IPAR_XXX
! ------------

    USE SNPAR
    USE PKMJDCOM

    IMPLICIT NONE

! function args

    REAL*8 MJD, XVAL(*)

! local args

    REAL*8  & 
         T0, TRISE, TFALL, A0, A1, A2  & 
        ,ARG1, ARG2, DUM1, DUM2, TDIF, POLYCOR  & 
        ,MODEL_FLUX

! ---------- BEGIN ---------

    T0    = XVAL(IPAR_T0)
    TRISE = XVAL(IPAR_TRISE)
    TFALL = XVAL(IPAR_TFALL)
    A0    = XVAL(IPAR_A0)
    A1    = XVAL(IPAR_A1)
    A2    = XVAL(IPAR_A2)

    TDIF  = MJD - T0

! numerator
    ARG1  = TDIF / Tfall
    DUM1  = EXP(-ARG1)

! denominator
    ARG2  = TDIF / Trise
    DUM2  = 1.0 + EXP(-ARG2)

! polynominal correction
    POLYCOR = A0 * ( 1.0 + A1*TDIF + A2*TDIF*TDIF )

! final flux

    MODEL_FLUX = POLYCOR * (DUM1/DUM2)

! Aug 31 2019: protect against crazy values.
    if ( MODEL_FLUX > +1.0E20 ) MODEL_FLUX =  1.0E20
    if ( MODEL_FLUX < -1.0E20 ) MODEL_FLUX = -1.0E20

    ANYLCFUN = MODEL_FLUX

    RETURN
  END FUNCTION ANYLCFUN

! =====================================
    DOUBLE PRECISION FUNCTION COVPKFLUXFUN(ifilt_obs, MJD1, MJD2)

! 
! Created April 12, 2010 by R.Kessler
! 
! Similar to COVARFLUX in snlc_fit.F90, but this is
! for the PKMJD fit to estimate error on peak flux.
! Since the user-function calls are different than
! in snlc_fit.F90, we must have a separate COV function
! here to evaluate errors.
! 
! Using full fitpar covariance, FITERRMAT(ipar1,ipar2)
!  ( see eq. 32.24 in PDG statistics section)
! 
!                           dF1   dF2
! U12       = \sum_{k,l}  ------ ------ * FITERRMAT_PKMJD(k,l)
!                         dPAR_k dPAR_l
! 
!  where PAR_k = fit-parameter with index 'k'.
!  and F1,F2 are fluxes at days with index 1,2.
! 
! If T1 = T2,  then U12 = square of model
! uncertainty at this epoch.
! 
! 
! 
! Jan 9, 2012: fix bug : use IFILT in SNLC_SNANAFIT_PEAKMJD_FITPAR
! -------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML
    USE FILTCOM
    USE PKMJDCOM

    IMPLICIT NONE

! subroutine args:

    INTEGER  & 
         IFILT_OBS  ! (I) ifilt_obs

    REAL*8  MJD1, MJD2  ! (I) MJDs to get covariance

! local var

    INTEGER  & 
         ipar, ipar_k, ipar_l, isp_k, isp_l, IFILT

    REAL*8  & 
         XVAL(NPAR_ANYLC), XERR(NPAR_ANYLC)  & 
        ,FLUX1, FLUX2  & 
        ,FTMP1, FTMP2  & 
        ,SAVEVAL_k, SAVEVAL_l  & 
        ,ERRPAR_k,  ERRPAR_l  & 
        ,dF1dVAL_k, dF2dVAL_l, FF  & 
        ,V_kl, U12

! functions
    REAL*8  ANYLCFUN

! -------------- BEGIN -------------

    COVPKFLUXFUN = 0.0
    U12 = 0.0   ! init output

! get reference fluxes to computer deriviates below.

    IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
    DO ipar = 1, NPAR_ANYLC
      XVAL(ipar) = SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar)
      XERR(ipar) = SNLC_SNANAFIT_PEAKMJD_FITERR(ifilt,ipar)
    ENDDO
    Flux1   =  ANYLCFUN ( MJD1, XVAL )
    Flux2   =  ANYLCFUN ( MJD2, XVAL )

    isp_k = 0
    DO 301 ipar_k = 1, NPAR_ANYLC
       IF ( XERR(ipar_k) .LE. 0.0 ) goto 301
       isp_k = isp_k + 1  ! for COV matrix

    isp_l = 0
    DO 302 ipar_l = 1, NPAR_ANYLC
       IF ( XERR(ipar_l) .LE. 0.0 ) goto 302
       isp_l = isp_l + 1 ! for COV matrix

       SAVEVAL_k = XVAL(ipar_k)
       SAVEVAL_l = XVAL(ipar_l)

       ERRPAR_k  = XERR(ipar_k)
       ERRPAR_l  = XERR(ipar_l)

! change fitpar values by +1 sigma and re-evaluate Flux.

       XVAL(ipar_k) = XVAL(ipar_k) + ERRPAR_k
       Ftmp1  = ANYLCFUN ( MJD1, XVAL )
       XVAL(ipar_k) = SAVEVAL_k

! repeat for 2nd fitpar ...
       XVAL(ipar_l) = XVAL(ipar_l) + ERRPAR_l
       Ftmp2   = ANYLCFUN ( MJD2, XVAL )
       XVAL(ipar_l) = SAVEVAL_l

! make sure that errors are non-zero before dividing.

       c2err = ' '
       if ( ERRPAR_k .EQ. 0.0 ) then
         write(c1err,662) PARNAME_STORE(ipar_k)
         CALL MADABORT("COVPKFLUXFUN", c1err, c2err )
       endif
       if ( ERRPAR_l .EQ. 0.0 ) then
         write(c1err,662) PARNAME_STORE(ipar_l)
         CALL MADABORT("COVPKFLUXFUN", c1err, c2err )
       endif
662      format('Cannot compute FLUXERR with ERR(',A,') = 0' )

! evaluate partial deriviates for flux
       dF1dVAL_k = (Ftmp1 - Flux1) / ERRPAR_k
       dF2dVAL_l = (Ftmp2 - Flux2) / ERRPAR_l

       FF   = dF1dVAL_k * dF2dVAL_l
       V_kl = FITERRMAT_PKMJD(isp_k,isp_l,ifilt_obs)
       U12 = U12 + FF * V_kl

!         print*,'  xxxx V_(',ipar_k,ipar_l, ') = ', V_kl,
!     &      '   RHO=', V_kl/(ERRPAR_k*ERRPAR_l)

302   CONTINUE
301   CONTINUE

    COVPKFLUXFUN = U12

    RETURN
  END FUNCTION COVPKFLUXFUN

! ============================
    SUBROUTINE SNMJD_INTERP_ANYLC(ifilt_obs)
! 
! 
! Jan 9, 2012: write output to LUNINTERP if SNMJD_OUT_FILE is defined.
! 

    USE SNDATCOM
    USE SNLCINP_NML
    USE FILTCOM
    USE INTERPCM

    IMPLICIT NONE

    INTEGER ifilt_obs   ! (I) absolute filter index

! local args

    INTEGER  NMJD, imjd,  LENCCID, ipar, IFILT, LUN

    REAL*8  & 
         FLUX, FLUXERR, MJD, COV, MJDLIST(MXINTERP)  & 
        ,XVAL(NPAR_ANYLC)  & 
        ,XERR(NPAR_ANYLC)

    CHARACTER cfilt*1, ccid*(MXCHAR_CCID)

    LOGICAL ADDFLAG  ! dummy arg

! functions
    REAL*8 COVPKFLUXFUN, ANYLCFUN

! ---------------- BEGIN --------------

    IF ( N_INTERP_MJDLIST .LE. 0 ) RETURN

    CCID = SNLC_CCID
    LENCCID = INDEX(CCID,' ') - 1

! get NMJD and MJDLIST;
    CALL GET_INTERP_MJDLIST(CCID, NMJD, MJDLIST, ADDFLAG)

    IF ( NMJD .EQ. 0 ) RETURN

    IFILT     = IFILTDEF_INVMAP_SURVEY(ifilt_obs)
    cfilt     = filtdef_string(ifilt_obs:ifilt_obs)

    IF ( SNMJD_OUT_FILE .EQ. ' ' ) THEN
       LUN = -1
    ELSE
       LUN = LUNINTERP
    ENDIF

    DO 200 imjd = 1, NMJD

      MJD  = MJDLIST(imjd)

! IF MJD=0, then replace with PEAKMJD (Aug 20, 2012)
      IF ( MJD .LT. 1.0  ) THEN
          MJD = SNLC_SEARCH_PEAKMJD
      ENDIF

      DO ipar = 1, NPAR_ANYLC
        XVAL(ipar) = SNLC_SNANAFIT_PEAKMJD_FITPAR(ifilt,ipar)
        XERR(ipar) = SNLC_SNANAFIT_PEAKMJD_FITERR(ifilt,ipar)
      ENDDO

      COV     = COVPKFLUXFUN(ifilt_obs,MJD,MJD)
      FLUX    = ANYLCFUN(MJD,XVAL)
      FLUXERR = SQRT(COV)

      write(6,60) SNLC_CCID(1:LENCCID),  & 
             cfilt, MJD, FLUX, FLUXERR
60      format(T4,'INTERP-FLUX(CID=',A, '-', A,  & 
               ' MJD=', F9.3, ')= ', G10.4, ' +- ', G10.4 )

      IF ( LUN .GT. 0 ) THEN
         write(LUN,70) SNLC_CCID(1:12),  & 
             cfilt, MJD, FLUX, FLUXERR
70         format('INTERP:',2x, A, 2x,A, 3x, F9.3, 2x, G12.4, G12.4)
      ENDIF

200   CONTINUE

    RETURN
  END SUBROUTINE SNMJD_INTERP_ANYLC

#endif



#if defined(SNFIT)

! ================================================================================
! ================================================================================
! Subroutines below are the base code for light curve fitting (e.g., snlc_fit)
! ================================================================================
! ================================================================================
! 
! 
! 
! ==============================================
    SUBROUTINE MNFIT_DRIVER (  & 
         CCID         &  ! (I) cand id
        ,NFITPAR      &  ! (I) number of fit parameters
        ,INIVAL       &  ! (I) initial parmater values
        ,INISTP       &  ! (I) initial step sizes (0=> fixed parameter)
        ,INIBND       &  ! (I) parameter bounds (0,0 => no bound)
        ,PARNAME      &  ! (I) list of parmater names
        ,USE_MINOS    &  ! (I) T=> use minos
        ,PRINT_LEVEL  &  ! (I) integer print level (-1=none)
        ,FITVAL        &  ! (O) final fit values
        ,FITERR_PLUS   &  ! (O) final fit errors, positive
        ,FITERR_MINUS  &  ! (O) final fit errors, negative
        ,FITCHI2       &  ! (O) min chi2
        ,NFIXPAR       &  ! (O) number of fixed parameters
        ,ERRTYPE       &  ! (O) error types (MINOS vs. PARABOLIC)
        ,MNSTAT_COV    &  ! (O) status of cov matrix (see minuit manual)
        ,IERR          &  ! (O) 0=>OK
            )
! 
! Created 2 Feb 28, 2006 by D.Cinabro and R.Kessler
! 
! Interface for minuit fitting
! 
! This is a SNANA routine ... not a MINUIT routine
! 
! Note that the non-integer/characater inputs are Double Precision!!!!
! For more information on using MINUIT, type minuit in Google and follor
! the URL with cern in it.
! 
! 
! Aug 31,  2009: Add USE_MINOS argument.
! 
! May 21, 2012: set print level using SNLCINP namelist MINUIT_PRINT_LEVEL
! 
! Jan 03 2016: pass new output arg MNSTAT_COV
! Apr 19 2022: set DO_PRINT for printing to suppress STDOUT for batch jobs
! May 08 2024: return IERR !=0  on NaN for any fit par value
! 
! -------------------------------------------------


    USE SNPAR
    USE CTRLCOM

    IMPLICIT NONE

! arguments

    CHARACTER CCID*(*)
    INTEGER   NFITPAR, NFIXPAR, PRINT_LEVEL

    LOGICAL USE_MINOS  ! (I)

    DOUBLE PRECISION  & 
         INIVAL(NFITPAR)  & 
        ,INISTP(NFITPAR)  & 
        ,INIBND(2,NFITPAR)       &  ! (I) fit bounds (zero=> no bound)
        ,FITVAL(NFITPAR)         &  ! (O) final fit values
        ,FITERR_PLUS(NFITPAR)    &  ! (O) final fit errors
        ,FITERR_MINUS(NFITPAR)   &  ! (O) final fit errors
        ,FITCHI2                ! (O) min chi2 value

    CHARACTER PARNAME(NFITPAR)*20  ! (MXCHAR_PARNAME)

    INTEGER ERRTYPE(NFITPAR), MNSTAT_COV
    INTEGER IERR, N_NAN

! local var

    EXTERNAL FCNSNLC, USRFUN
    DOUBLE PRECISION  USRFUN

    INTEGER  & 
         IPAR, IVARBL, i, IFLAG, NARG  & 
        ,NPARI, NPARX, ISTAT

    DOUBLE PRECISION  & 
         FIXLIST(NFITPAR)  & 
        ,BND1, BND2  & 
        ,STRATEGY(4), PARG(4)  & 
        ,MAXCALLS(4)  & 
        ,GRAD(NFITPAR)  & 
        ,CHI2  & 
        ,ARGLIST(20)  & 
        ,EPLUS, EMINUS, EPARAB, GLOBCC  & 
        ,FEDM, ERRDEF  & 
        ,ERRSYM      ! local SYMMMETRIC fiterr

    LOGICAL LFIX, DO_PRINT

! Unit numbers for input and output

    INTEGER IUNIN,IUOUT,IUWRI
    DATA IUNIN/5/
    DATA IUOUT/6/
    DATA IUWRI/7/

! -------------------- BEGIN -------------------------


    IERR = 0
    N_NAN   = 0
    NFIXPAR = 0
    FITCHI2 = 0.0
    DO_PRINT = STDOUT_UPDATE
! 
! Initialize MINUIT, input unit, output unit, save unit
! 

    CALL MNINIT(IUNIN,IUOUT,IUWRI)

    PARG(1) = DBLE(PRINT_LEVEL)
    CALL MNEXCM(FCNSNLC, 'SET PRI', PARG, 1, IERR, USRFUN )

! Define parameters in MINUIT

    DO IPAR = 1, NFITPAR

       CALL MNPARM(IPAR, PARNAME(IPAR)  & 
             ,INIVAL(IPAR), INISTP(IPAR)  & 
             ,INIBND(1,IPAR)    &  ! lo-bound
             ,INIBND(2,IPAR)    &  ! hi-bound
             ,IERR   )

       IF ( IERR .NE. 0 ) THEN
	    IERR = ERRFLAG_MNFIT_INITPAR
          WRITE (6,'(A,A10,A)') '  ERROR initializing ',  & 
                  PARNAME(ipar), '  with MNPARM '
          RETURN
       ENDIF

! init fit params
       FITVAL(ipar)       = 0.0
       FITERR_PLUS(ipar)  = 0.0
       FITERR_MINUS(ipar) = 0.0
       ERRTYPE(ipar)      = 0

! Set list of fixed parameters with step size = 0

       IF ( INISTP(IPAR) .EQ. 0.0 ) THEN
          NFIXPAR = NFIXPAR + 1
          FIXLIST(NFIXPAR) = DFLOAT(IPAR)
       ENDIF

    ENDDO

    CALL FLUSH(6)

! --------------------------------------------
! fix param if its step size is zero

    IF ( NFIXPAR .GT. 0 ) THEN

       CALL MNEXCM(FCNSNLC,'FIX',FIXLIST, NFIXPAR,IERR,USRFUN)
       IF ( IERR .NE. 0 ) THEN
          IERR = ERRFLAG_MNFIT_FIXPAR
          WRITE (6,'(A)') '  ERROR fixing parameters'
          RETURN
       ENDIF
    ENDIF

! call FCN function with init flag: IFLAG=2
    IFLAG = 2
    CALL FCNSNLC ( NFITPAR, GRAD, chi2, FITVAL, IFLAG, USRFUN )

! Set Strategy 2 for more accurate derivative calculation

    NARG = 1
    STRATEGY(1) = DBLE(2.0)   ! 1=default,  2 => better results

    IF ( DO_PRINT ) THEN
      print*,'  Set MINUIT STRATEGY = ', int(STRATEGY(1))
    ENDIF

    CALL MNEXCM(FCNSNLC,'SET STR', STRATEGY, NARG, IERR, USRFUN )

    MAXCALLS(1) = dble(30000.0)

! Actually do the fit

    IF ( DO_PRINT ) THEN
      PRINT *,' '
      PRINT *,' ------------------------------------------------ '
    ENDIF

    IF ( USE_MINOS ) THEN
      CALL MNEXCM(FCNSNLC, 'MINOS', MAXCALLS, NARG, IERR, USRFUN)
      IF ( DO_PRINT ) THEN
        PRINT *,'  MNFIT_DRIVER: MINOS returns IERR = ', IERR
      ENDIF
    ELSE
      CALL MNEXCM(FCNSNLC, 'MINIMIZE', MAXCALLS, 0, IERR, USRFUN)
      IF ( DO_PRINT ) THEN
        PRINT *,'  MNFIT_DRIVER: MIGRAD returns IERR = ', IERR
      ENDIF
    ENDIF

    IF ( DO_PRINT ) THEN
       PRINT *,' ------------------------------------------------ '
    ENDIF

! -------------------------------------------
! get the fit result

    IF ( DO_PRINT ) THEN
      PRINT *,' '
      PRINT *,'   MNFIT_DRIVER: extract fit parameters with MNPOUT'
      CALL FLUSH(6)
    ENDIF

    DO 44 IPAR = 1, NFITPAR

      LFIX = ( INISTP(ipar) .EQ. 0.0 )

      CALL MNPOUT(  IPAR, PARNAME(IPAR)  & 
                    , FITVAL(IPAR), ERRSYM  & 
                    , BND1,BND2, IVARBL )

      IF ( ISNAN(FITVAL(IPAR))  ) THEN
	   N_NAN = N_NAN + 1
	   write(6,644) PARNAME(IPAR), CCID
644	   format(' MNFIT_DRIVER ERROR: ', A,' = NaN for CID = ', A)
         call flush(6)
      ENDIF

      CALL MNERRS( IPAR, EPLUS, EMINUS, EPARAB, GLOBCC )

! xxxxxxxxxxxxxxxxxxxxxxxx
!        print*, IPAR, ': ERR(EPARAB)=', sngl(EPARAB)
!     &            ,'  ERR(+-) = ', sngl(EPLUS), sngl(EMINUS)
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! take average of +- errors (note that EMINUS is negative)
! If MINOS errors are zero, then use parabolic error above.

      if ( EPLUS .GT. 0.0 .and. EMINUS .LT. 0.0 ) THEN
         FITERR_PLUS(IPAR)   = EPLUS
         FITERR_MINUS(IPAR)  = EMINUS
         ERRTYPE(ipar)       = ERRTYPE_MINOS
      else
         FITERR_PLUS(IPAR)   = +ERRSYM
         FITERR_MINUS(IPAR)  = -ERRSYM
         ERRTYPE(ipar)       = ERRTYPE_PARAB
      endif

44    CONTINUE  ! end of loop over IPAR


    IF ( N_NAN > 0 ) THEN
       IERR = ERRFLAG_MNFIT_NAN
       RETURN  ! May 2024
    ENDIF

! -------------------------------------------------------
! get chi2 and status; fill FITCHI2 return arg
! From Minuit manual, here are the ISTAT values from MNSTAT:
! 
! ISTAT  meaning
!  0    Not calculated at all
!  1    Diagonal approximation only, not accurate
!  2    Full matrix, but forced positive-definite
!  3    Full accurate covariance matrix
!          (After MIGRAD, this is the indication of normal convergence.)

    CALL MNSTAT(FITCHI2, FEDM, ERRDEF, NPARI, NPARX, ISTAT)
    MNSTAT_COV = ISTAT  ! 3=OK

! ----------------------------------


! Stop fit and make LAST call to function

    PARG(1) = 0.0
    CALL MNEXCM(FCNSNLC, 'EXIT',  PARG, 0, IERR, USRFUN )

    IF ( DO_PRINT ) THEN
      PRINT *,'   MNFIT_DRIVER: EXIT returns IERR = ', IERR
      CALL FLUSH(6)
    ENDIF

! ---------------------------------
    RETURN
  END SUBROUTINE MNFIT_DRIVER 

! ===========================
    SUBROUTINE MNFIT_STOREPAR(iter,IERR)
! 
! Nov 17, 2011: replace CERNLIB's PROB with snana's PROB_CHI2NDOF
! Jan 04, 2016: protect sqrt(negative SQERR8)
! May 05, 2018: init FITERRMAT=0 instead of -9. Allows USESIM_xxx=T
! May 23, 2018: set LCFRACERRDIF_STORE(ipar)
! Feb 22, 2020: set NDOF_PRIOR=0  here on each iteration.
! Apr 16, 2024: write FITPROB to stdout
! Jul 30, 2024: protect EPLUS/EMINUS when EMINUS=0


    USE SNDATCOM
    USE SNLCINP_NML
    USE SNANAFIT

    IMPLICIT NONE

! subroutine args

    INTEGER  & 
         iter    &  ! (I) iterationn
        ,IERR   ! (O) error flag; 0=> OK

! local variables

    INTEGER ipar, NDOF, ipar1, ipar2, i
    INTEGER J, JPAR(MXFITPAR), J1, J2, etmp
    REAL EPLUS, EMINUS, PCHI2
    LOGICAL LBADERR, FLOATBOTH

    REAL*8  & 
         CHI8, GRAD8(MXFITPAR)  & 
        ,EMAT8, ERR8(2), SQERR8(2), ERR_FINAL, ERR_LAST, ERR_AVG

! function

    REAL*8   PROB_CHI2NDOF
    REAL*8   USRFUN
    EXTERNAL USRFUN

! -------------- BEGIN ------------------

    IERR = 0   ! set output flag to OK

! store chi2

    FITCHI2_STORE(1)     = FITCHI2_MIN ! total chi2

! get prior-chi2 using special flag

    CALL FCNSNLC(NFITPAR_MN, GRAD8, CHI8, FITVAL(1,iter),  & 
             FCNFLAG_PRIOR_ONLY, USRFUN)

! determine effective Ndof for prior.

    if ( CHI8 .GT. 0.01 ) then
      NDOF_PRIOR = 1
    else
      NDOF_PRIOR = 0
    endif

    FITCHI2_STORE(3) = CHI8  ! prior chi2
    FITCHI2_STORE(2) =  & 
           FITCHI2_STORE(1) - FITCHI2_STORE(3)  ! total - prior

! store chi2 from ln(sigma) terms
    CALL FCNSNLC(NFITPAR_MN, GRAD8, CHI8, FITVAL(1,iter),  & 
             FCNFLAG_SIGMA_ONLY, USRFUN)
    FITCHI2_STORE(4) = CHI8

! store number of degrees of freedom

    NDOF          = NEPOCH_FIT(0) - (NFITPAR_MN - NFIXPAR)
    NDOF_STORE(1) = NDOF + NDOF_PRIOR  ! full chi2
    NDOF_STORE(2) = NDOF               ! chi2 excluding prior
    NDOF_STORE(3) = NDOF_PRIOR         ! prior chi2 only
    NDOF_STORE(4) = 1                  ! ln(sigma)

! Convert fitchi2 into fitprobs.

    DO i     = 1, 4
      NDOF   = NDOF_STORE(i)
      CHI8   = DBLE(FITCHI2_STORE(i))

      if ( NDOF .GT. 0 .and. CHI8 .GT. 0.0 ) THEN
         PCHI2  = SNGL( PROB_CHI2NDOF(CHI8,NDOF) )
      ELSE
         PCHI2  = 1.0
      ENDIF

      FITPROBCHI2_STORE(i) = PCHI2
      LCCHI2_STORE(i)      = FITCHI2_STORE(i)
      LCPROBCHI2_STORE(i)  = FITPROBCHI2_STORE(i)
    ENDDO

    IF ( ITER==1 ) FITPROB_ITER1 = FITPROBCHI2_STORE(1)  ! July 2024


! now store fit parameters.

    DO 444 ipar = 1, NFITPAR_MN

       FLOATPAR(ipar) = ( INISTP(ipar) .NE. 0.0 )

! if there are no floated params, then set errors
! to user-initialized values: FITERR(ipar)

       if ( NFIXPAR .EQ. NFITPAR_MN  ) then
         FITERR_PLUS(ipar,iter)   = FITERR(ipar,iter)
         FITERR_MINUS(ipar,iter)  = FITERR(ipar,iter)
       endif

       EPLUS  = FITERR_PLUS(ipar,iter)
       EMINUS = ABS ( FITERR_MINUS(ipar,iter) )

! if fit error is way too small, set error to INISTP
! and set error type to BAD

       LBADERR = Eplus  .LT. ERRMAX_BAD(ipar)  & 
           .and.   Eminus .LT. ERRMAX_BAD(ipar)

! 9/22/2007: allow tiny errors when the exposure time is large.
       IF ( SIM_EXPOSURE_TIME(1) .GT. 10.0 ) then
          LBADERR = .FALSE.
       endif

! set bad-error flag of fit-value has not moved from initial value.
       LBADERR = LBADERR .or.  & 
                  ( INIVAL(ipar) .EQ. FITVAL(ipar,iter) )

       if ( FLOATPAR(ipar) .and. LBADERR ) then
          FITERR_PLUS(ipar,iter)  = +INISTP(ipar)
          FITERR_MINUS(ipar,iter) = -INISTP(ipar)
          EPLUS              = FITERR_PLUS(ipar,iter)
          EMINUS             = ABS ( FITERR_MINUS(ipar,iter) )
          ERRTYPE(ipar)      = ERRTYPE_BAD
       endif

       FITERR_RATIO(ipar,iter)   = EPLUS / (EMINUS+1.0E-10)
       FITERR(ipar,iter)         = 0.5 * ( EPLUS + EMINUS )

       FITVAL_STORE(ipar)    = FITVAL(ipar,iter)
       FITERR_STORE(ipar)    = FITERR(ipar,iter)
       INIVAL_STORE(ipar)    = INIVAL(ipar)
       ERRTYPE_STORE(ipar)   = ERRTYPE(ipar)

! store FITVCAL in 'LC' array (might get over-written later by PDFVAL

       LCVAL_STORE(ipar) = FITVAL_STORE(ipar)
       LCERR_STORE(ipar) = MIN(99.9,FITERR_STORE(ipar))

! print results for floated parameters only

       if ( FLOATPAR(ipar)  & 
           .or. NFIXPAR .EQ. NFITPAR_MN  & 
           .or. IPAR    .EQ. 1  & 
                 ) then
          CALL PRINT_FITPAR(SNLC_CCID, iter, ipar, 'fit' )
       endif

! do stuff on last iteration
       IF ( iter .EQ. NFIT_ITERATION .and. FLOATPAR(ipar) ) then
          etmp = ERRTYPE(ipar)
          NERRTYPE(etmp) = NERRTYPE(etmp) + 1

          if ( iter > 1 ) then
            ERR_FINAL = FITERR(ipar,iter)
            ERR_LAST  = FITERR(ipar,iter-1)
            ERR_AVG   = 0.5*(ERR_FINAL+ERR_LAST)
            if ( ERR_AVG > 1.0E-12 ) then
              LCFRACERRDIF_STORE(ipar) = (ERR_FINAL-ERR_LAST)/ERR_AVG
            endif
          endif
       ENDIF

444   CONTINUE   ! end loop over IPAR


! - - - - - -
    if ( STDOUT_UPDATE ) then
       write(6,446) SNLC_CCID(1:ISNLC_LENCCID),  & 
                      iter, FITPROBCHI2_STORE(1)
446      format(T8, 'CID ', A, 3x, 'ITER=',I1, 3x, 'FITPROB = ', E10.3 )
       CALL FLUSH(6)
    endif

! ----------------------
! 9/29/2007: get error matrix

    IF ( LREPEAT_MINOS ) MNSTAT_COV = MNSTAT_COV + 10
    CALL MNEMAT(FITERRMAT_SPARSE,MXFITPAR)

! compute correlation matrix from error matrix

    J = 0
    DO ipar = 1, NFITPAR_MN
      if ( FLOATPAR(ipar) ) then
         j = j + 1
         JPAR(ipar) = j
      endif
    ENDDO

    DO ipar1 = 1, NFITPAR_MN
    DO ipar2 = 1, NFITPAR_MN

       FITERRMAT(ipar1,ipar2) = 0.0
       FITCORMAT(ipar1,ipar2) = 0.0

       FLOATBOTH = ( FLOATPAR(ipar1) .and. FLOATPAR(ipar2) )
       if ( FLOATBOTH ) then
         J1 = JPAR(ipar1)
         J2 = JPAR(ipar2)
         EMAT8     = FITERRMAT_SPARSE(J1,J2)
         SQERR8(1) = FITERRMAT_SPARSE(J1,J1)
         SQERR8(2) = FITERRMAT_SPARSE(J2,J2)

         IF  ( SQERR8(1) > 0.0 .and. SQERR8(2) > 0.0 ) then
           ERR8(1) = sqrt ( SQERR8(1) )
           ERR8(2) = sqrt ( SQERR8(2) )
           FITCORMAT(ipar1,ipar2) = EMAT8/(ERR8(1) * ERR8(2) )
           FITERRMAT(ipar1,ipar2) = EMAT8
         else

         endif

       endif
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE MNFIT_STOREPAR

! ======================================
    CHARACTER*1 FUNCTION ERRTYPE_STR(ERRTYPE)

    USE SNPAR

    IMPLICIT NONE

    INTEGER ERRTYPE  ! (I) error type

! --------------- BEGIN -------------

    if ( ERRTYPE .EQ. ERRTYPE_MINOS ) then
      ERRTYPE_STR = 'M'
    else if ( ERRTYPE .EQ. ERRTYPE_PARAB ) then
      ERRTYPE_STR = 'P'
    else if ( ERRTYPE .EQ. ERRTYPE_MARG ) then
      ERRTYPE_STR = 'm'
    else if ( ERRTYPE .EQ. ERRTYPE_BAD ) then
      ERRTYPE_STR = '?'
    else if ( ERRTYPE .EQ. 0 ) then
      ERRTYPE_STR = ' '
    else
      ERRTYPE_STR = '?'
    endif

    RETURN
  END FUNCTION ERRTYPE_STR

! =======================================
    SUBROUTINE PRINT_FITPAR(CCID, iter, ipar, type )
! 
! standard print-line for fit-result or pdf-result.
! 
! type = 'fit' : uses current FITVAL and FITERR
! type = 'pdf' : uses current PDFVAL and PDFERR
! 
! Note that ITER is used only for 'fit' option
! 
! Jun 2013: return if we are no longer in VERBOSE mode and iter < NFIT_ITER
!            (to reduce output for BIG jobs)
! 

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! subroutine args

    INTEGER  ITER, IPAR  ! (I) SN cand id and IPAR to print
    CHARACTER  & 
          CCID*(*)    &  ! (I) char string cand id
         ,type*(*)   ! (I) 'fit' or 'pdf' type of result

! local var

    INTEGER LL
    REAL*8  VAL, ERR
    LOGICAL LQUIET
    CHARACTER STRING_VAL*12

! function
    CHARACTER ERRTYPE_STR*1

! --------------- BEGIN -------------

    LQUIET = .NOT. STDOUT_UPDATE

    IF ( LQUIET ) RETURN

    IF ( type .EQ. 'fit' ) then
      VAL = FITVAL(ipar,iter)
      ERR = FITERR(ipar,iter)
    ELSE IF ( type .EQ. 'pdf' ) then
      VAL = PDFVAL(ipar)
      ERR = PDFERR(ipar)
    ELSE

    ENDIF

! write par value into string using format
    WRITE(STRING_VAL,211) VAL
211   FORMAT(G11.5)
    IF ( VAL > 1.0E4 ) THEN  ! this is PEAKMJD
      WRITE(STRING_VAL,212) VAL
212     FORMAT(F11.3)
    ENDIF

    LL = INDEX(CCID,' ') - 1
    IF ( LL .LE. 0 ) LL = MXCHAR_CCID
    WRITE (6,430)  CCID(1:LL), type, IPAR  & 
          , PARNAME_STORE(IPAR), STRING_VAL, ERR  & 
          , ERRTYPE_STR ( ERRTYPE(IPAR) )

 430    FORMAT (T8,'CID ',A,2x, A,'par(', I2, '):',A14,' = ',  & 
                       A11,' +/- ', G10.4, '(',A,')' )


    CALL FLUSH(6)

    RETURN
  END SUBROUTINE PRINT_FITPAR


! ======================================
    SUBROUTINE PDF_INIT()
! 
! Created Nov 24, 2009 by R.Kessler
! 
! Init PDFXXX arrays before marginalizing or before running
! the MCMC option.  The main issue here is to make sure that
! fixed parameters (like ITER) get transfered to
! the PDFVAL array.
! 
! ---------------------------

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! local var
    INTEGER IPAR, IPAR2

! --------------- BEGIN -----------------

    DO 40 ipar = 1, NFITPAR_MN
       PDFVAL(ipar) = FITVAL(ipar,NFIT_ITERATION)

! reset PDF-marginalized values for floated parameters
       IF ( FLOATPAR(ipar) ) then
          PDFERR(ipar)  = -9.0
          PDFVAL(ipar)  = -9.0
       ENDIF

! zero the covariance matrix

       DO ipar2 = 1, NFITPAR_MN
         PDFERRMAT(ipar,ipar2) = 0.0
         PDFCORMAT(ipar,ipar2) = 0.0
       ENDDO
40    CONTINUE

    RETURN
  END SUBROUTINE PDF_INIT
! ======================================
    SUBROUTINE PDF_STORE()
! 
! Created Nov 11, 2009 by R.Kessler
! 
! Utility to store marginalized/MCMC PDF values,
! and re-compute chi2 and fit-probs.
! [Code moved from end of MARG_DRIVER]
! 
! Nov 17, 2011: replace CERNLIB's PROB() with PROB_CHI2NDOF
! 
! -----------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! local var

    INTEGER IPAR, NDOF, i

    DOUBLE PRECISION  & 
          GRAD8(MXFITPAR)  & 
         ,CHI8, USRFUN

    REAL  PCHI2
    EXTERNAL USRFUN

! function
    REAL*8 PROB_CHI2NDOF

! ----------------- BEGIN -----------

! print/store results

    DO 30 ipar = 1, NFITPAR_MN

       PDFVAL_STORE(ipar)   = PDFVAL(ipar)
       PDFERR_STORE(ipar)   = PDFERR(ipar)
       PDFPROB2_STORE(ipar) = PDFPROB2(ipar)

! over-write LCVAL[ERR] array with pdf-avarges

       LCVAL_STORE(ipar) = PDFVAL_STORE(ipar)
       LCERR_STORE(ipar) = MIN(99.9,PDFERR_STORE(ipar))

       if ( .NOT. FLOATPAR(ipar)  ) goto 30

       CALL PRINT_FITPAR( SNLC_CCID, 0, ipar, 'pdf' )
       CALL FLUSH(6)

30    CONTINUE

! ----------------------------------------------------
! re-compute chi2, prior-chi2 and prob with marginalized values.

! start with PRIOR-only chi2.
    CALL FCNSNLC(NFITPAR_MN, GRAD8, CHI8, PDFVAL,  & 
             FCNFLAG_PRIOR_ONLY, USRFUN)
    LCCHI2_STORE(3) = CHI8

! now get contribution from ln(sigma)-terms

    CALL FCNSNLC(NFITPAR_MN, GRAD8, CHI8, PDFVAL,  & 
             FCNFLAG_SIGMA_ONLY, USRFUN)
    LCCHI2_STORE(4) = CHI8

! Now get full chi2 ...

    CALL FCNSNLC(NFITPAR_MN, GRAD8, CHI8, PDFVAL,  & 
            FCNFLAG_LAST, USRFUN)
    LCCHI2_STORE(1) = CHI8

! data-only chi2 if full-prior chi2.
    LCCHI2_STORE(2) = CHI8 - LCCHI2_STORE(3)

! compute fit-probs.

    DO i     = 1, 4
      chi8   = LCCHI2_STORE(i)
      NDOF   = NDOF_STORE(i)

      if ( NDOF .GT. 0 .and. CHI8 .GT. 0.0 ) then
         PCHI2  = SNGL( PROB_CHI2NDOF(CHI8,NDOF) )
      else
         PCHI2 = 1.0
      endif

      LCPROBCHI2_STORE(i) = PCHI2
    ENDDO

    CALL FLUSH(6)

    RETURN
  END SUBROUTINE PDF_STORE


! =======================================
    SUBROUTINE MARG_DRIVER( HOFF_MARG, OPT,  & 
                 MAX_INTEGPDF, NGRID_FINAL, NSIGMA)
! 
! Created Aug 3, 2006 by R.Kessler
! 
! Utility to compute pdf-averaged quantity for each fitpar.
! 
! Call this routine AFTER fit is done ...
! this routine uses current FITVAL and FITERR,
! and fills PDFVAL(ipar) and PDFERR(ipar).
! 
! Histograms book/filled:
!   HOFF:         PDF value for each FCNPDF function call
!   HOFF + ipar : 1-dim PDF distribution for each fitted ipar
! 
! 
!  May 5, 2007: fix dumb bug. Need to init PDVAL(ipar) = FITVAL(ipar)
!               before integration to make sure that fixed parameters
!               are set in PDFVAL
! 
! Oct 16, 2009: call INTEGPDF twice. First time set NGRID=7 to get
!               better estimate of errors. Second time set
!               NGRID = NGRID_FINAL. Change should run faster
!               if there are fewer iterations needed with NGRID_FINAL.
! 
! Nov 24, 2009: call PDF_INIT() to init PDFXXX arrays
! 
! -------------------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! input args

    INTEGER  & 
          HOFF_MARG       &  ! (I) fill plots with hbook offset = HOFF
         ,OPT             &  ! (I) options
         ,MAX_INTEGPDF    &  ! (I) max number of iterations
         ,NGRID_FINAL    ! (I) # bins for each integrated dimension

    REAL  NSIGMA     ! (I) integrate +_ NSIGMA for exact pdf.

! local var

    INTEGER  & 
         ipar, ipar2  & 
        ,JTIME1, JTIME2, JDIFTIME  & 
        ,NEVAL         &  ! number of function evaluations
        ,LL, NDOF  & 
        ,i, IERR  & 
        ,NGRID, HOFF, NHDIM, HID_PDF, NBPDF(2)

    character chis*80, copt*6

    REAL*8 XMIN(2), XMAX(2), TMP
    LOGICAL LSYMERR

! functions
    CHARACTER ERRTYPE_STR*1
    REAL*8  PROB

! FCN args

    DOUBLE PRECISION  & 
          GRAD8(MXFITPAR)  & 
         ,CHI8  & 
         ,USRFUN

    EXTERNAL USRFUN

! ------------------ BEGIN ------------

! require at least a few bins ...

    if ( NGRID_FINAL .LE. 0 ) RETURN

! set chi2 value for FCN function to quit

    IF ( PDFMIN .GT. 0.0 .and.  & 
              OPT .EQ. OPT_INTEGPDF_QUITCHI2 ) THEN
      FITCHI2_QUIT = -2.0*DLOG(PDFMIN) + FITCHI2_MIN
    ELSE
      FITCHI2_QUIT = 1.0E20
    ENDIF

    COPT = 'GRID'   ! only option so far

! -----------------

    write(6,19) SNLC_CCID, NSIGMA, copt
19    format(/,T5,'MARG_DRIVER(CID ',A6,'): ',  & 
            'compute P.D.F(+- ',F3.1,' sigma)', 2x, 'method=',A)

    IF ( FITCHI2_QUIT .LT. 1.0E19 ) then
      print*,'   CPU-saver ON  => FCNSNLC quits when CHI2 > ',  & 
               FITCHI2_QUIT
    ELSE
      print*,'   CPU-saver OFF => FCNSNLC always computes full CHI2.'
    ENDIF

    CALL FLUSH(6)

    USEPDF_MARG  = .TRUE.
    ISTAGE_SNANA = ISTAGE_TEST

! --------------------------------------


! init PDFVAL = FITVAL so that fixed parameters get transfered.

    CALL PDF_INIT()

! ---------------------
! get exact pdf by integrating over other fit-parameters

    JTIME1 = TIME()

! First marginalize with just 7 grid-points per variable.

    OPT   = 1   ! 1st round estimate
    NGRID = 7
    HOFF  = 0   ! skip histograms
    CALL INTEGPDF( OPT, HOFF,  & 
              MAX_INTEGPDF, NGRID, DBLE(NSIGMA), NEVAL, IERR )

! final marginalization; use previous PDF for grid size estimate

    OPT   = 2
    NGRID = NGRID_FINAL
    HOFF  = HOFF_MARG
    CALL INTEGPDF( OPT, HOFF,  & 
              MAX_INTEGPDF, NGRID, DBLE(NSIGMA), NEVAL, IERR )

! compute integration time.

    JTIME2   = TIME()
    JDIFTIME = JTIME2 - JTIME1

    LL = INDEX(SNLC_CCID,' ') - 1
    write(6,80) NEVAL, JDIFTIME, SNLC_CCID(1:LL)
80    format(T5,'MARG_DRIVER: Finished ',I7,' function calls in ',  & 
               I4,' seconds  (SN ', A,')'   )
    print*,' '

! keep track of integration times.

    NCALL_INTEGPDF     = NCALL_INTEGPDF + 1
    TIME_INTEGPDF      = JDIFTIME
    TIMESUM_INTEGPDF   = TIMESUM_INTEGPDF + JDIFTIME
    tmp                = DBLE(TIMESUM_INTEGPDF)/DBLE(NCALL_INTEGPDF)
    TIMEAVG_INTEGPDF   = INT(TMP+0.5)

    if ( mod(NCALL_INTEGPDF,5) .EQ. 0 ) then
       print*,' '
       print*,'   (AVERAGE INTEGPDF TIME: ',  & 
                    TIMEAVG_INTEGPDF,'  seconds)'
       print*,' '
    endif

! call utility to store PDF results
    CALL PDF_STORE()

    RETURN
  END SUBROUTINE MARG_DRIVER


! =======================================
    SUBROUTINE INTEGPDF(OPT, HOFF,  & 
                 MAX_INTEGPDF, NGRID, NSIGMA, NEVAL, IERR )
! ---------------------
!  Retruns p.d.f(DLMAG) integratged over other parameters;
!  integration is from +-NSIGMA * FITERR over each
!  dimension with non-zero INISTP.
! 
!  The calling routine must set FITVAL(IPAR_DLMAG) = DLMAG,
!  and also set INISTP(IPAR_DLMAG) = 0.0 so that FCNPDF
!  knows to ignore the DLMAG-dimension in the integration.
! 
!  OPT=1 => first estimate with small NGRID
!  OPT=2 => final estimate with final NGRID
! 
! 
! histograms are booked / filled for
! 
!  HOFF          : function value for each call
!  HOFF + ipar   : pdf for each floated "ipar"
! 
! 
!  Feb 24, 2007: major upgrade to iterate if problem
!                is detected. See LREDO logic.
! 
!                PDFERR(ipar) is now the RMS of the pdf distribution.
! 
! Apr 28, 2007:
!  on 2nd iteration when prob at edge is too high, make more robust
!  estimate of integration region. Previously, integ-region was extended
!  by three times the shift in PDFVAL. Now, a Gaussian profile is
!  assumed, and an effective SIGMA is computed based on prob(at edge)
!  and current PDFVAL(ipar).  The integration limmit is then
!  changed to PDFVAL + NSIGMA*SIGMA
!  This improvement should help when MINUIT returns an error that
!  is way too small, but is still not flagged by BADERR.
! 
! May 3, 2007: accept MAX_INTEGPDF as argument
! 
! Aug 20, 2008: change MXPAR from 8 to 10
!               (to accomodate IPAR_LUMIPAR2 in STRETCH2 model)
! 
! Oct 16, 2009:
!      use OPT=1,2 to determine which NGRID-iteration
!      Compute covariance & correlations: PDFCORMAT(ipar1,ipar2)
! 
!      Fill PDRPROB2(ipar)
! 
! Jan 4, 2010: add IERR argument. For PDF=0 error, abort only
!              if OPT=2. This gives both NGRID values a chance
!              to  succeed.
! 
! Oct 01, 2012: use LCPLOT utility instead of HBOOK1 and HPAK
! 
! Feb 06, 2013: replace LCPLOT util with SNHIST
! 
! Jun 10 2013: protect ABORT when LPDFZERO=T using user namelist
!              ABORT_ON_MARGPDF0
! 
! -------------------------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! function aarguments

    INTEGER  & 
         OPT       &  ! (I) option
        ,MAX_INTEGPDF   &  ! (I) max # times to integrate
        ,NGRID     &  ! (I) # grid-bins for each dimension
        ,HOFF      &  ! (I) hbook offset
        ,NEVAL     &  ! (O) number of function calls.
        ,IERR     ! (O) 0=>OK

    REAL*8  NSIGMA   ! (I) integrate +- NSIGMA in each dimension

! local args
    INTEGER, PARAMETER :: &
          MXPAR   = 10    &  ! max number of fit paramters/dimensions
         ,MXGRID  = 30  

    INTEGER  & 
         IPAR, IPAR2  & 
        ,NBINTOT  & 
        ,IBIN  & 
        ,IBIN_OFF  & 
        ,IGRID  & 
        ,NDIM, IDIM, IDIM2  & 
        ,IPAR_DIM(MXPAR)   &  ! IPAR for each dimension to integrate
        ,IBIN_DIM(MXPAR)   &  ! local grid-bin for each dimension
        ,NN, i  & 
        ,NPASS  & 
        ,HID, NHDIM, LL, NUM  & 
        ,ITER  & 
        ,NPDF, IBIN_PLOT, NB(2)

    REAL*8  & 
         PARVAL_MIN(MXPAR)  & 
        ,PARVAL_MAX(MXPAR)  & 
        ,PARVAL_BINSIZE(MXPAR)  & 
        ,PARDIF_MIN(MXPAR)  & 
        ,PARDIF_MAX(MXPAR)  & 
        ,PARVAL(MXPAR)  & 
        ,PARVAL_GRID(MXGRID,MXPAR)  & 
        ,TMP, TMPVAL  & 
        ,DVOL  & 
        ,X8(MXPAR)  & 
        ,PDF, XPDF(2), WGT  & 
        ,PDFWSUMCOR(0:MXPAR,0:MXPAR)    &  ! wgted sum for correlations
        ,PDFSUMCOR  & 
        ,PDFWSUM(0:MXPAR)           &  ! wgted sum for each ipar
        ,PDFSUM                     &  ! PDF sum for each ipar
        ,PDFSUM_GRID(MXGRID,MXPAR)   &  ! PDF vs. par to plot 1-d pdf
        ,PDFMAX(MXPAR)              &  ! max over grid for each ipar
        ,TMP_RANGE(2,MXPAR)  & 
        ,PDF1D(MXGRID,MXPAR)  & 
        ,SUM0, SUM1, SUM2, XN, XTMP  & 
        ,SQERR, E12, E1xE2, PDFTMP  & 
        ,PDF_NBR1, PDF_NBR2  & 
        ,XMIN(2), XMAX(2), XVAL(2)


    LOGICAL  & 
         LTMP  & 
        ,LPDFZERO  & 
        ,LREDO  & 
        ,LREDO_ALL  & 
        ,LDMP_DEBUG

    CHARACTER chis*80, choice*12

! function

    REAL*8   FCNPDF
    EXTERNAL FCNPDF

! ----------------- BEGIN ------------

    NEVAL    = 0  ! init output arg
    IERR     = 0

    NPASS    = 0

    ITER = NFIT_ITERATION

    DO ipar = 1, NFITPAR_MN
       PARVAL_MIN(ipar) = 0.0
       PARVAL_MAX(ipar) = 0.0
       PARDIF_MIN(ipar) = 0.0
       PARDIF_MAX(ipar) = 0.0
    ENDDO

    LDMP_DEBUG = .FALSE.

    NHDIM    = 1  ! 1D histo

! =====================================
2     CONTINUE
    NPASS = NPASS + 1

! init some useful things.

    DVOL      = 1.0
    NBINTOT   = 1
    NDIM      = 0  ! number of dimensions to integrate PDF
    PDFSUM    = 0.0
    PDFSUMCOR = 0.0
    LREDO_ALL = .FALSE.

! set up integration limits for each parameter

    DO 30 ipar = 1, NFITPAR_MN

       PDFWSUM(ipar) =  0.0
       PDFMAX(ipar)  =  0.0

       do ipar2 = 1, NFITPAR_MN
          PDFWSUMCOR(ipar,ipar2) = 0.0
       enddo

       do igrid = 1, NGRID
          PDFSUM_GRID(igrid,ipar) = 0.0
       enddo

       if ( .not. FLOATPAR(ipar)  ) goto 30 ! skip fixed params

       NBINTOT = NBINTOT * NGRID

! keep track of which parameters to integrate
! (i.e., to ignore fixed parameters)

       NDIM = NDIM + 1
       IPAR_DIM(NDIM) = IPAR
       TMP_RANGE(1,ipar) = PARVAL_MIN(ipar)
       TMP_RANGE(2,ipar) = PARVAL_MAX(ipar)

       CALL INTEGRANGE(IPAR,NGRID,NSIGMA,PDF1D(1,ipar),   &  ! inputs
                LREDO, TMP_RANGE(1,ipar) ) ! outputs are LREDO   RANGE(1:2)

! set global REDO flag if any parameter-range is adjusted.

       LREDO_ALL = LREDO_ALL .OR. LREDO

30    CONTINUE  ! end of IPAR loop

! --------------------------------------------
! check if integration should proceed.
! Allow no more than three tries ... give "BEWARE" warning
! after 3 tries.

    IF ( NPASS .GT. 1 ) THEN

      IF ( LREDO_ALL ) THEN
         LL = INDEX(SNLC_CCID,' ') - 1

         if ( NPASS .LE. MAX_INTEGPDF ) then
            print*,'  ==> Integrate ',SNLC_CCID(1:LL),  & 
                  ' again with NGRID=', NGRID
            CALL FLUSH(6)
         else if ( OPT .EQ. 1 ) then
            goto 800

         else if ( OPT .EQ. 2 ) then  ! warn on final NGRID only

            print*,'  ==> Store ',SNLC_CCID(1:LL),  & 
                ' result, but BEWARE !!!'
            CALL FLUSH(6)
!              print*,'  ==> Cannot converge for ',SNLC_CCID(1:LL)
!              IERR = -9  ! Dec 16, 2011

            GOTO 800  ! skip integration
         endif

      ELSE
         GOTO 800  ! skip integration
      ENDIF
    ENDIF

      CALL FLUSH(6)

! ----
! if we get here, then adjust integration range

    DO 31 ipar = 1, NFITPAR_MN

      if ( .not. FLOATPAR(ipar)  ) goto 31 ! skip fixed params

      PARVAL_MIN(ipar) = TMP_RANGE(1,ipar)
      PARVAL_MAX(ipar) = TMP_RANGE(2,ipar)

      TMP = FITVAL(ipar,ITER)
      PARDIF_MIN(ipar) = PARVAL_MIN(ipar) - TMP
      PARDIF_MAX(ipar) = PARVAL_MAX(ipar) - TMP

      TMP = PARVAL_MAX(ipar) - PARVAL_MIN(ipar)
      PARVAL_BINSIZE(ipar) = TMP / float(NGRID)

! compute volume element (used in brute-force method)

      DVOL = DVOL *  PARVAL_BINSIZE(ipar)

31    CONTINUE


! Now do the integration

    LPDFZERO = .TRUE.

    DO 770 IBIN = 1, NBINTOT

! determine local grid-bin for each dimension to integrate;
! the load local PARVAL with value at each grid-point.

       ibin_off = 0
       do idim  = 1, NDIM
          NN    = NGRID**(NDIM-idim)
          ibin_dim(idim) = (ibin - 1 - ibin_off)/NN + 1
          ibin_off       = ibin_off + (ibin_dim(idim) - 1) * NN

          ipar     = ipar_dim(idim)  ! fetch fit par index
          igrid    = ibin_dim(idim)
          TMP      = float( igrid ) - 0.5

          XTMP     = PARVAL_MIN(ipar)  & 
                     + PARVAL_BINSIZE(ipar) * TMP

          PARVAL(ipar)            = XTMP
          PARVAL_GRID(igrid,ipar) = XTMP

       enddo  ! end loop of NDIM

! get X8 array that contains only parameters to integrate
! (fixed parameters are weeded out from PARVAL)

       CALL FITVAL_FLOAT(PARVAL, NDIM, X8) ! returns NDIM and X8

       PDF   = FCNPDF(NDIM,X8)   ! evaluate normalized PDF

       NEVAL = NEVAL + 1         ! increment # function calls

       IF ( PDF .EQ. 0.0 ) goto 771

       LPDFZERO = .FALSE.

       PDFSUM   = PDFSUM + PDF  ! increment total integral

       do idim  = 1, NDIM
          ipar  = ipar_dim(idim)  ! fetch fit par index

          PDFWSUM(ipar) = PDFWSUM(ipar) + PDF * PARVAL(ipar)

          do idim2 = 1, NDIM
             ipar2 = ipar_dim(idim2)
             PDFWSUMCOR(ipar,ipar2) = PDFWSUMCOR(ipar,ipar2)  & 
                      + PDF * PARVAL(ipar) * PARVAL(ipar2)
          enddo

! sum pdf separately for each grid point and for each IPAR
! so that 1-dim PDF can be plotted for each IPAR.

          igrid = ibin_dim(idim)
          PDFSUM_GRID(igrid,ipar) =  & 
            PDFSUM_GRID(igrid,ipar) + PDF

          if ( PDFSUM_GRID(igrid,ipar) .GT. PDFMAX(ipar) ) then
            PDFMAX(ipar) = PDFSUM_GRID(igrid,ipar)
          endif
       enddo  ! end of idim loop

        if ( PDF > 0.0 .and. HOFF > 0 ) then
          xpdf(1)   = DLOG10(PDF)
          xpdf(1)   = max ( -19.999, xpdf(1) )
          wgt       = 1.0
          CALL SNHIST_FILL(NHDIM, HOFF, XPDF, WGT )
        endif

771       continue
        if ( MOD(IBIN,10000)  .EQ. 0 ) then
           print*,'      Processing grid-bin ',  & 
                 ibin,'/', NBINTOT
           CALL FLUSH(6)
        endif

770   CONTINUE

    IF ( LPDFZERO ) THEN
      print*,' '
      print*,'  WARNING: INTEGPDF ERROR for CID=', SNLC_CCID
      print*,'  pdf function is zero everywhere with NGRID=',NGRID
      print*,' '
      IERR = -9

      IF ( OPT .EQ. 2 .and. ABORT_ON_MARGPDF0 ) THEN
        print*,' ***** ABORT ***** '
        CALL EXIT(EXIT_ERRCODE)
      ELSE
        RETURN
      ENDIF

    ENDIF

! ---------------------------------------
! fill histograms and PDFVAL_STORE array

    DO idim  = 1, NDIM
       ipar  = ipar_dim(idim)  ! fetch fit par index

       PDFVAL(ipar) = PDFWSUM(ipar) / PDFSUM

       DO igrid = 1, NGRID
          PDF = PDFSUM_GRID(igrid,ipar) / PDFMAX(ipar)

          if ( PDF .LT. 1.0E-30 ) THEN
            PDF1D(igrid,ipar) = 1.0E-20  ! avoid hbook bit problems
          else
            PDF1D(igrid,ipar) = PDF
          endif

       ENDDO  ! end of igrid loop

    ENDDO  ! end of IDIM loop

! determine PDF error (RMS) without using hbook so
! that we can pass HOFF=0 and skip histograms

    DO idim  = 1, NDIM
       ipar  = ipar_dim(idim)  ! fetch fit par index
       SUM0 = 0.0
       SUM1 = 0.0
       SUM2 = 0.0
       NPDF = 0
       PDFPROB2(ipar) = 0.0

    DO igrid = 1, NGRID

       XTMP   = PARVAL_GRID(igrid,ipar)
       PDFTMP = PDF1D(igrid,ipar)
       if ( PDFTMP .GT. 1.0E-6 ) NPDF = NPDF + 1

       SUM0   = SUM0 + PDFTMP
       SUM1   = SUM1 + PDFTMP * XTMP
       SUM2   = SUM2 + PDFTMP * XTMP * XTMP

       PDF_NBR1 = 0.0
       PDF_NBR2 = 0.0
       if ( igrid .GT. 1     )  PDF_NBR1 = PDF1D(igrid-1,ipar)
       if ( igrid .LT. NGRID )  PDF_NBR2 = PDF1D(igrid+1,ipar)

! check for 2nd local maximum
       if (    PDFTMP .LT. .99  & 
           .and. PDFTMP .GT. PDF_NBR1  & 
           .and. PDFTMP .GT. PDF_NBR2 ) then
          PDFPROB2(ipar) = PDFTMP
       endif
    ENDDO  ! igrid

      SQERR = SUM2/SUM0 - (SUM1/SUM0)**2

      if ( SQERR .LT. -0.00001 ) THEN ! allow for numerical rounding
              print*,' '
              print*,' SUM[0,1,2] = ', SUM0, SUM1, SUM2
              print*,' (SUM2/SUM0)     = ', SUM2/SUM0
              print*,' (SUM1/SUM0)**2  = ', (SUM1/SUM0)**2
              print*,' SQERR           = ', SQERR
              print*,' NPASS=',NPASS,'  OPT=',OPT
              c1err = ' SWERR < 0 for ' // PARNAME_STORE(ipar)
              CALL MADABORT("INTEGPDF", c1err, "  ")
       endif

! to compute error from RMS, require more than 1 PDF bin to be non-zero
! (to avoid pathologies from PDFERR -> 0)

       if ( SQERR .GT. 0.0 .and. NPDF .GT. 1 ) then
          PDFERR(ipar) = SQRT( SQERR )
       else
          PDFERR(ipar) = FITERR(ipar,iter)  ! PDF err = fit err for now
       endif

       ERRTYPE(ipar) = ERRTYPE_MARG
    ENDDO  ! idim


! -------------------------------------------------
! evaluate covariances

    DO 400 idim  = 1, NDIM
       ipar      = ipar_dim(idim)
    DO 401 idim2 = 1, NDIM
       ipar2     = ipar_dim(idim2)

       E12   = PDFWSUMCOR(ipar,ipar2)/PDFSUM
       E1xE2 = PDFVAL(ipar)*PDFVAL(ipar2)
       PDFERRMAT(ipar,ipar2) = E12 - E1xE2

       SQERR = PDFERR(ipar)*PDFERR(ipar2)
       PDFCORMAT(ipar,ipar2) = PDFERRMAT(ipar,ipar2)/SQERR

401   CONTINUE
400   CONTINUE


! ------------------------
     IF ( LDMP_DEBUG ) THEN

        print*,' xxx ====================================== '
        PRINT*,' xxx NPASS = ', NPASS

       DO idim  = 1, NDIM
          ipar  = ipar_dim(idim)

          print*,' - - - - - - - - - - - - - - - - - '
          write(6,665) PARNAME_STORE(ipar),  & 
            PDFVAL(ipar), PDFERR(ipar), PDFPROB2(ipar)

          write(6,666) PARNAME_STORE(ipar),  & 
            PARVAL_MIN(ipar), PARVAL_MAX(ipar)

          write(6,667) 'PARVAL',  & 
                 ( PARVAL_GRID(igrid,ipar), igrid=1,11)
          write(6,667) 'PDFVAL',  & 
                 ( PDF1D(igrid,ipar), igrid=1,11)

665        format(T3, 'xxx ',A6, 2x, 'PDFVAL = ',  & 
                G10.3,' +- ', G10.3, 3x, 'PROB2=',F5.3 )
666        format(T3, 'xxx ',A6, 2x, 'PARDIF(MIN,MAX)=',2F9.3)

667        format(T3,'xxx ',A, '=', 11F7.3 )
         CALL FLUSH(6)
       ENDDO
     ENDIF
! ------------------------

!  ------------------------------------------------------
! always go back to start to check if we have converged

         goto 2

! ==================================================
! Check option to plot PDF for each floated parameter

800   CONTINUE
    CALL FLUSH(6)

    IF ( HOFF .LE. 0 ) RETURN

    DO idim  = 1, NDIM
       ipar  = ipar_dim(idim)  ! fetch fit par index
       hid   = HOFF + ipar

       LL = index ( PARNAME_STORE(ipar), ' ' ) - 1
       write(chis,21)  & 
                PARNAME_STORE(IPAR)(1:LL)  & 
              , PARNAME_STORE(IPAR)(1:LL)  & 
              , SNLC_CCID(1:ISNLC_LENCCID), char(0)

21       format(' margin. PDF( ',A,'-',A,'(fit) ) for CID=',A, A)


       TMPVAL    = FITVAL(ipar,NFIT_ITERATION)
       xmin(1)   = PARVAL_MIN(ipar) - TMPVAL
       xmax(1)   = PARVAL_MAX(ipar) - TMPVAL
       NB(1)     = NGRID

       CALL SNHIST_INIT(NHDIM, HID, CHIS//char(0),  & 
                 NB, XMIN, XMAX, LEN(chis) )

       DO igrid = 1, NGRID
         TMP      = DBLE(igrid) - 0.5
         XVAL(1)  = XMIN(1) + PARVAL_BINSIZE(ipar) * TMP
         WGT      = PDF1D(igrid,ipar)
         CALL SNHIST_FILL( NHDIM, HID, XVAL, WGT )
       ENDDO

    ENDDO

    RETURN
  END SUBROUTINE INTEGPDF


! =======================================
    SUBROUTINE INTEGRANGE(IPAR,NGRID,NSIGMA,PDFLAST,  & 
            LREDO,RANGE)

! 
! Created Apr 30, 2007 by R.Kessler
! 
! Determine integration range for IPAR parameter based on
! 1-dim PDF values from previous integration.
! 
! On first pass, RANGE(1)=RANGE(2)=0 and integration range
! is just based on NSIGMA * FITERR from MINUIT.
! 
! Oct 23, 2007: TMP_PROB -> min (0.98, GRIDPROB)
!               to avoid dividing by ln(near 1)
! 
! Dec 3, 2007:  min(0.98,GRIDPROB) -> min(0.95,GRIDPROB)
!               so that integ-range is not opened so much
!               (short-term fix for photoZ fit to 14888)
! 
! Oct 17, 2009: define PDFLAST_MAX = .95 with mutiple PDF bins,
!               or =.5 with just one PDF bin. Avoids severe bin-
!               extensions when there is just one non-zero PDF bin.
! 
!               Makes sure that RANGE(1:2) is withing INIBND
! 
! Dec 16, 2011:
!    - write more info on PDFPROBLEM at endge of 1 dim PDF.
!    - for PHOTOZ, skip upper INIBND check on RANGE(2)
! 
! -----------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! subroutine args

    INTEGER  & 
         IPAR    &  ! (I) parameter to determine integration range
        ,NGRID  ! (I) number of grid points to integrate

    REAL*8  & 
         NSIGMA           &  ! (I)
        ,PDFLAST(NGRID)   &  ! (I) 1-dim PDF at grid points of last integration
        ,RANGE(2)        ! (I,O) integration limits: old -> new

    LOGICAL LREDO      ! (I) flag to redo integration

! local args

    LOGICAL  & 
         LSYMERR  & 
        ,LTMP  & 
        ,LPDFMAX  & 
        ,BADEDGE(2)  & 
        ,LFIRST  & 
        ,USE_PDFERR

    INTEGER  & 
         NREDO  & 
        ,N1D_TOT            &  ! total number of 1D bins = NGRID
        ,N1D_PROB0          &  ! number of 1D bins with negligible prob
        ,NPDF_GOOD  & 
        ,igrid  & 
        ,igrid_min  & 
        ,igrid_max  & 
        ,L1, L2, i  & 
        ,ITER

    REAL*8  & 
         XNSIG, XNGRID, xgrid  & 
        ,RANGE_OLD(2)  & 
        ,RANGE_NEW(2)  & 
        ,BINSIZE_OLD  & 
        ,BINSIZE_NEW  & 
        ,TMP_PLUS  & 
        ,TMP_MINUS  & 
        ,TMP_SIG  & 
        ,TMP_PROB  & 
        ,TMP_EDGE  & 
        ,TMP_PDF  & 
        ,TMP_VAL  & 
        ,SHIFT, BND  & 
        ,PDFMAX_atEDGE  & 
        ,PDFLAST_MAX

    CHARACTER PDFPROBLEM(40)*48

! -------------- BEGIN -------------

    ITER = NFIT_ITERATION

! init a bunch of stuff.

    XNSIG         = NSIGMA
    XNGRID        = FLOAT(NGRID)

    LFIRST = (RANGE(1) .EQ. 0.0 .and. RANGE(2) .EQ. 0.0)

! store old integration range since it gets over-written later.
    RANGE_OLD(1) = RANGE(1)
    RANGE_OLD(2) = RANGE(2)
    BINSIZE_OLD  = (RANGE_OLD(2)-RANGE_OLD(1))/XNGRID

    NREDO        = 0
    LREDO        = .FALSE.  ! output arg !!
    TMP_PLUS     = 0.0
    TMP_MINUS    = 0.0

    IF ( LFIRST ) THEN

       USE_PDFERR = PDFERR(ipar) .GT. 0.0
       LSYMERR = abs(FITERR_RATIO(ipar,ITER)-1.0) .LT. 0.20

       IF ( USE_PDFERR ) THEN  ! PDF error gives better estimate
         TMP_PLUS    =  XNSIG * PDFERR(ipar)
         TMP_MINUS   = -XNSIG * PDFERR(ipar)  ! neg number
       ELSE IF ( LSYMERR ) THEN
         TMP_PLUS    =  XNSIG * FITERR(ipar,ITER)
         TMP_MINUS   = -XNSIG * FITERR(ipar,ITER)  ! neg number
       ELSE
         TMP_PLUS    =  XNSIG * FITERR_PLUS(ipar,ITER)
         TMP_MINUS   =  XNSIG * FITERR_MINUS(ipar,ITER)  ! neg number
       ENDIF

       RANGE(1) = FITVAL(ipar,ITER) + TMP_MINUS
       RANGE(2) = FITVAL(ipar,ITER) + TMP_PLUS

       LREDO = .TRUE.
       RETURN

    ENDIF

! ---------------------------------------------------
! loop over grid points and find bin with first,last
! non-negligible PDF value.
! Also count NPDF_GOOD = number of bins with PDF > PDFMIN

    N1D_TOT    = 0
    N1D_PROB0  = 0
    LPDFMAX    = .FALSE.
    igrid_min  = 1
    igrid_max  = NGRID

    DO 100 igrid = 1, NGRID

       TMP_PDF = PDFLAST(igrid)
       if ( TMP_PDF .GT. 0.90 ) LPDFMAX = .TRUE.

       N1D_TOT = N1D_TOT + 1

       if ( TMP_PDF .GT. PDFMIN_GOOD ) then
           NPDF_GOOD = NPDF_GOOD + 1
       else
           N1D_PROB0 = N1D_PROB0 + 1

           if ( LPDFMAX .and. IGRID_MAX .EQ. NGRID ) then
              igrid_max = igrid  ! first negligible pdf past max
           endif

           if ( .not. LPDFMAX ) then
              igrid_min = igrid  ! last negligible pdf before max
           endif
       endif

100   CONTINUE  ! end of igrid loop

    IF ( NPDF_GOOD .GT. 1 ) THEN
        PDFLAST_MAX = 0.95
    ELSE
        PDFLAST_MAX = 0.80  ! Guass approx is bad with just one bin
    ENDIF


! -----------------------------------------------------------
! if PDF at edge is too big, then compute adjustment to extend edge.
! Adjustent is based on assumption of Gaussian profile.

    SHIFT      = PDFVAL(ipar) - FITVAL(ipar,ITER)
    BADEDGE(1) = .FALSE.
    BADEDGE(2) = .FALSE.
    PDFMAX_atEDGE = MAX ( PDFLAST(1), PDFLAST(NGRID) )

    LTMP = PDFLAST(1) .GT. PDFMAX_EDGE
    IF ( LTMP  ) THEN
       BADEDGE(1) = .TRUE.
       TMP_PROB   = min ( PDFLAST(1), PDFLAST_MAX )
       TMP_EDGE   = RANGE_OLD(1) - PDFVAL(ipar) ! distance to min edge
       TMP_SIG    = TMP_EDGE / SQRT ( -2.0 * DLOG(TMP_PROB) )
       TMP_MINUS  = XNSIG * TMP_SIG - TMP_EDGE

! avoid making integration window smaller
       if ( TMP_MINUS .GT. 0.0 ) then
            TMP_MINUS = -abs(3.*SHIFT)
       endif
    ENDIF



    LTMP = PDFLAST(NGRID) .GT. PDFMAX_EDGE
    IF ( LTMP ) THEN
        BADEDGE(2) = .TRUE.
        TMP_PROB   = min ( PDFLAST(NGRID), PDFLAST_MAX )
        TMP_EDGE   = RANGE_OLD(2) - PDFVAL(ipar) ! distance to edge
        TMP_SIG    = TMP_EDGE / SQRT ( -2.0 * DLOG(TMP_PROB) )
        TMP_PLUS   = XNSIG * TMP_SIG - TMP_EDGE

! avoid making integration window smaller
        if ( TMP_PLUS .LT. 0.0 ) then
             TMP_PLUS = +abs(3.*SHIFT)
        endif
    ENDIF

! ----------------------------
! update new integration range to account for
! extended range(s).

    RANGE_NEW(1) = RANGE_OLD(1) + TMP_MINUS
    RANGE_NEW(2) = RANGE_OLD(2) + TMP_PLUS
    BINSIZE_NEW  = (RANGE_NEW(2)-RANGE_NEW(1))/XNGRID

! -------------------------------------
! check for bins that have ~0 prob, and reduce integration
! range to eliminate these 0-prob bins.

    if ( igrid_min .GT. 1 ) then
       xgrid        = float(igrid_min) - 0.5
       TMP_VAL      = RANGE_OLD(1) + BINSIZE_OLD * xgrid
       RANGE_NEW(1) = TMP_VAL - BINSIZE_NEW/2.0
    endif

    if ( igrid_max .LT. NGRID ) then
       xgrid        = float(igrid_max) - 0.5
       TMP_VAL      = RANGE_OLD(1) + BINSIZE_OLD * xgrid
       RANGE_NEW(2) = TMP_VAL + BINSIZE_NEW/2.0
    endif

! store NEW integration range in output subroutine arg.

    RANGE(1) = RANGE_NEW(1)
    RANGE(2) = RANGE_NEW(2)

! ---------------------------------------------------
! check of there is any reason that integration needs to be redone.
!  - Prob(edge) is too high
!  - Only 1 bin in 1D PDF
!  - Too many Zero-prob bins'
! 
!  Set LREDO flag = T if there is a problem.

    IF ( BADEDGE(1)  ) then
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),441) PDFLAST(1), 'lo'
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),1441) RANGE_OLD, RANGE_NEW
    ENDIF
    IF ( BADEDGE(2)  ) then
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),441) PDFLAST(NGRID), 'hi'
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),1441) RANGE_OLD, RANGE_NEW
    ENDIF

441        format('PROB=',F6.4, ' at ',A,' edge of 1 dim PDF')
1441       format('RANGE=(', F6.2, ',', F6.2, ')' , ' -> ',  & 
                        '(', F6.2, ',', F6.2, ')'   )  ! Dec 2011

    IF ( NPDF_GOOD .LT. 3  ) then
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),442)
442        format('<= 2 non-zero PDF bins.' )
    ENDIF

    IF ( N1D_PROB0 .GT. 3 ) THEN
         NREDO = NREDO + 1
         write(pdfproblem(NREDO),443) N1D_PROB0, N1D_TOT
443        format('PROB=0 for ', I3, '/' , I3, 2x,  & 
                 'PDF bins' )
    ENDIF

! ----------------------------------------------
! print list of PDF-problems to STDOUT

    IF ( NREDO .GT. 0 ) THEN

       LREDO = .TRUE.

       L1 = index ( PARNAME_STORE(ipar), ' ' ) - 1
       L2 = index ( SNLC_CCID, ' ' ) - 1

       do i = 1, NREDO
          write(6,449)  & 
                PARNAME_STORE(ipar)(1:L1)  & 
              , SNLC_CCID(1:L2)  & 
              , PDFPROBLEM(i)
449         format(T8,A,'-PDF PROBLEM(',A,') : ', A)
       enddo

    ENDIF

! don't allow range to go outside of official INIBND range.

    BND = INIBND(1,IPAR)
    IF (RANGE(1) .LT. BND ) RANGE(1) = BND

! special exception for PHOTOZ
    IF ( PARNAME_STORE(IPAR)(1:6) .NE. 'PHOTOZ' ) THEN
      BND = INIBND(2,IPAR)
      IF (RANGE(2) .GT. BND ) RANGE(2) = BND
    ENDIF

    RETURN
  END SUBROUTINE INTEGRANGE

! =======================================
    DOUBLE PRECISION FUNCTION FCNPDF(NDIM,X)
! 
!  Retruns exp[ - (chi2 - chi2min) / 2 ]
! 
!  where CHI2 is returned from a call to FCNSNLC.
!  This function is designed to accept only the
!  floated parameters so that it can be used by
!  more generic integration routines that don't
!  know about fixed parameters.
! 
! Example: fit-params 1,2 are fixed and 3-6 float.
!          Then pass NDIM=4 and X(4) = array of floated parameters.
! 
! -----------------------

    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! function input

    INTEGER NDIM     ! (I) number of dimensions
    REAL*8  X(NDIM)  ! (I) args for chi2-function

! local args

    INTEGER  & 
         ipar_all  & 
        ,ipar_float  & 
        ,ipar

! FCNSNLC args

    INTEGER IFLAG

    DOUBLE PRECISION  & 
          GRAD(MXFITPAR)  & 
         ,FITVAL_LOC(MXFITPAR)  & 
         ,CHI2, DELCHI2  & 
         ,USRFUN

    EXTERNAL USRFUN

! ---------- BEGIN -----------------

    FCNPDF  = 0.0
    IFLAG   = FCNFLAG_FAST

! restore input parameters X to local FITVAL array;
! note that fixed parameters are loaded from FITVAL_STORE array.

    ipar_float = 0

    DO ipar_all = 1, NFITPAR_MN

      if ( FLOATPAR(ipar_all)  ) then  ! is a floated param
         ipar_float = ipar_float + 1
         FITVAL_LOC(ipar_all) = X(ipar_float)
      else
         FITVAL_LOC(ipar_all) = FITVAL(ipar_all,NFIT_ITERATION)
      endif

    ENDDO

!      write(6,20) (FITVAL_LOC(ipar), ipar=1, 6)

    CALL FCNSNLC(NFITPAR_MN,GRAD,CHI2,FITVAL_LOC,IFLAG,USRFUN)

! set FCNPDF only if FCN function did not have its
! chi2 got too big.

    if ( chi2 .LE. FITCHI2_QUIT ) then
       DELCHI2 = chi2 - FITCHI2_MIN
       FCNPDF  = dexp(-DELCHI2/2.0)
    endif

    RETURN
  END FUNCTION FCNPDF



! ==============================================
    SUBROUTINE FITVAL_FLOAT(FITVAL_LOC,NDIM,X)
! 
! FITVAL is the array of NFITPAR_MN parameters.
! Returns NDIM = number of floated parmaters,
! and X(NDIM) = array of floated paramters.
! i.e, the fixed paramters are ignored.
! 
! -----------------------------------


    USE SNDATCOM
    USE SNANAFIT
    USE SNLCINP_NML

    IMPLICIT NONE

! declare function args

    REAL*8  FITVAL_LOC(*)     ! (I) fit parameters
    INTEGER NDIM              ! (O) number of floated params
    REAL*8  X(*)              ! (O) array of floated params

! local args
    INTEGER ipar
! --------------- BEGIN -------------

    NDIM = 0
    DO ipar = 1, NFITPAR_MN
      if ( FLOATPAR(ipar) ) then
         NDIM    = NDIM + 1
         X(NDIM) = FITVAL_LOC(ipar)
      endif
    ENDDO

    RETURN
  END SUBROUTINE FITVAL_FLOAT

#endif

! END:




