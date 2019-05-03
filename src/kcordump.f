CDECK  ID>,  MAIN. 
      PROGRAM KCORDUMP

      IMPLICIT NONE

C +CDE,DUMPCOM. inserted below by ypatchy. 

      REAL
     &   REDSHIFT_DMP_BINSIZE
     &  ,TREST_DMP_MIN
     &  ,TREST_DMP_MAX
     &  ,TREST_DMP_BINSIZE

      PARAMETER (
     &   REDSHIFT_DMP_BINSIZE = 0.01
     &  ,TREST_DMP_MIN        = -20.0
     &  ,TREST_DMP_MAX        = +80.0
     &  ,TREST_DMP_BINSIZE    =   1.0
     &    )


      CHARACTER
     &   STRING_Kxy*8, STRING_REST_COLOR*4  ! Apr 2016

      INTEGER
     &   IFILT_REST_DMP       ! (I) rest-frame filter index
     &  ,IFILT_OBS_DMP        ! (I) observer frame filter index
     &  ,IFILT_RESTCOLOR_DMP(2) ! (I) rest-color (default=0.0) (Apr 2016)

      REAL
     &   Trest_DMP       ! (I) rest-frame epoch, days
     &  ,redshift_DMP    ! (I) redshift
     &  ,REDSHIFT_DMP_RANGE(2)
     &  ,TREST_DMP_RANGE(2)
     &  ,MAGREST_COLOR, MAGREST(2)

      INTEGER
     &   NZBIN_DMP, NTBIN_DMP

      COMMON / DUMPCOM /
     &   STRING_Kxy,     STRING_REST_COLOR, MAGREST_COLOR, MAGREST
     &  ,IFILT_REST_DMP, IFILT_OBS_DMP, IFILT_RESTCOLOR_DMP
     &  ,Trest_DMP,      redshift_DMP
     &  ,TREST_DMP_RANGE,   REDSHIFT_DMP_RANGE
     &  ,NZBIN_DMP,         NTBIN_DMP

C ==============================================

C +CDE,SNDATCOM. inserted below by ypatchy. 


C +CDE,SNPAR. inserted below by ypatchy. 

      CHARACTER  SNTABLE_LIST_DEFAULT*60

c parameters used by snana code

      INTEGER
     &   MXVERS, MXSURVEY, MXSNLC,MXCID, MXCID_CHECK,MXEPOCH, MXITER
     &  ,MXFILT_ALL, MXFILT_OBS, MXFILT_REST, ONE
     &  ,MXFILT_KCOR, MXFILT_SNRMAX, MXEP_MODELGRID
     &  ,MXIDTEL, MXIDFIELD, MXFIELD_OVP, MXSEASON, MXSNHOST
     &  ,MXTABLE1D_KCOR,MXTABLE1D_LCMAG,MXTABLE1D_MWXT,MXTABLE1D_AVWARP
     &  ,MXZBIN_KCOR, MXTBIN_KCOR, MXAVBIN_KCOR, MXCBIN_AVWARP
     &  ,MNTYPE, MXTYPE, MXLISTNML, MXCCID_LIST
     &  ,MXVAR_PRIVATE, MXCUT_PRIVATE, MXVAR_TERSE
     &  ,HOFF_NOCUTS, HOFF_CUTS, HOFF_SIM
     &  ,LUNHIS, LUNNML, LUNDAT, LUNTMP, LUNFIT, LUNPKMJD, LUNDMP
     &  ,LUNCID, LUNOUT
     &  ,LUNRES1, LUNRES2, LUNINTERP, LUNLIST, LUNIGNORE, LUNSALT2
     &  ,ISTAGE_INIT, ISTAGE_RDSN, ISTAGE_CUTS, ISTAGE_USRANA
     &  ,ISTAGE_TEST
     &  ,INTERP_LINEAR, INTERP_SMOOTH, INTERP_ZSMOOTH
     &  ,MXERRTYPE, ERRTYPE_MINOS, ERRTYPE_PARAB, ERRTYPE_MARG
     &  ,ERRTYPE_BAD
     &  ,FCNFLAG_USER, FCNFLAG_FAST, FCNFLAG_LAST, FCNFLAG_USESIM
     &  ,FCNFLAG_PRIOR_ONLY, FCNFLAG_SIGMA_ONLY
     &  ,OPT_INTEGPDF_QUITCHI2, OPT_INTEGPDF_FULL
     &  ,OPT_SNXT_CCM89, OPT_SNXT_SJPAR
     &  ,OPT_MWCOLORLAW_DEFAULT, OPT_MWEBV_DEFAULT
     &  ,OPT_KCORERR_SJ, OPT_KCORERR_SJ5, OPT_KCORERR_SMOOTH
     &  ,MXLINE_ARGS, MXEPOCH_IGNORE
     &  ,NPAR_ANYLC, MXPAR_SIMSED, MXPAR_LCLIB
     &  ,MXLAMBIN_SNSED, MXCUTBIT
     &  ,OPT_FILTUPD_EACHSN, OPT_FILTUPD_MAP, OPT_FILTUPD_SAMEFILT
     &  ,OPT_FILTOBS, OPT_FILTREST
     &  ,ISTAT_READAGAIN, ISTAT_SKIP
     &  ,IDTABLE_MCMC, NFIT_VERBOSE
     &  ,ITABLE_SNANA,  ITABLE_FITRES, ITABLE_SNLCPAK, ITABLE_SPECPAK
     &  ,ITABLE_SPECTRA, MXTABLE
     &  ,IDTABLE_SNANA, IDTABLE_FITRES, IDTABLE_SPECTRA
     &  ,IDTABLE_CIDPTR
     &  ,OPT_PARSTORE_TEXTTABLE
c
     &  ,MODEL_STRETCH
     &  ,MODEL_STRETCH2
     &  ,MODEL_MLCS2k2
     &  ,MODEL_SNOOPY
     &  ,MODEL_SALT2
     &  ,MODEL_SIMSED
     &  ,MODEL_BYOSED
     &  ,MODEL_NON1A
     &  ,MODEL_LCLIB    ! Sep 2017
     &  ,MODEL_FIXMAG   ! force mags to user-value
     &  ,MXMODEL_INDEX
c
     &  ,MXCHAR_CCID, MXCHAR_VERSION, MXCHAR_SURVEY
     &  ,MXCHAR_PATH, MXCHAR_FILENAME, MXCHAR_MODELNAME
     &  ,MXCHAR_FIELDNAME, MXCHAR_PARNAME, MXCHAR_CUTNAME
     &  ,MXCHAR_FILEWORD
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA     ! DATA flux per epoch
     &  ,SNLCPAK_EPFLAG_REJECT   ! REJECT flag (1=> excluded from fit)
     &  ,SNLCPAK_EPFLAG_CHI2     ! data-fit chi2 per epoch
     &  ,SNLCPAK_EPFLAG_FITFUN   ! smooth  fitfun curve
     &  ,SNLCPAK_EPFLAG_FLUXSIM  ! SIM flux per epoch
     &  ,SNLCPAK_EPFLAG_FLUXREST ! rest-frame flux per epoch (optional)
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST ! idem for sim truth
     &  ,SNLCPAK_BANDFLAG_PKFLUX   ! peak flux vs. filter
     &  ,SNLCPAK_BANDFLAG_PKMJD    ! peak MJD  vs. filter
     &  ,SNLCPAK_BANDFLAG_NDOF     ! Ndof vs. filter
     &  ,SNLCPAK_BANDFLAG_CHI2     ! chi2 vs. filter
c
     &  ,IFLAG_INI, IFLAG_ANA, IFLAG_END
     &  ,MAG_SATURATE
     &  ,MASK_FLUXCOR_SNANA, MASK_FLUXERRCOR_SNANA
c
     &  ,EXIT_ERRCODE_SNANA, EXIT_ERRCODE_SNFIT, EXIT_ERRCODE_PSNID

      REAL*8
     &   ZERO8, ONE8, TEN8,  PI, LOGTEN
     &  ,PDFMIN, PDFMAX_EDGE, PDFMIN_GOOD
     &  ,KCORPACK
     &  ,XTMW_FRACERR
     &  ,MJDOFF
     &  ,ZEROPOINT_FLUXCAL_DEFAULT
     &  ,RV_MWCOLORLAW_DEFAULT
     &  ,CUTVAL_OPEN, IGNORE_HEADVAL

      REAL NULLVAL

      PARAMETER (                   ! BEGIN SNANA PARAMS
     &   SNTABLE_LIST_DEFAULT = 'SNANA  FITRES  LCPLOT'
     &  ,MXVERS        = 50         ! max number of versions to read
     &  ,MXSURVEY      = 100        ! max number of SURVEY_NAMEs to store
     &  ,MXITER        = 12         ! max # fit iterations
     &  ,NFIT_VERBOSE  = 500        ! Number of fits for verbose printouts
     &  ,MXSNLC        = 4000000    ! max number of SNe (fits format)
     &  ,MXCCID_LIST   = 10000      ! max size of SNCCID_LIST_ALL
     &  ,MXCID       = 299 999 999  ! max CID = 300 million -1 (9 digits)
     &  ,MXCID_CHECK =  99 999 999  ! MXCID to check duplicates
     &  ,MXEPOCH     = 2000        ! max number of filter-epochs per SN
     &  ,MXEP_MODELGRID = 200     ! max model grid for SNANA+SIM_MAGOBS table
     &  ,MXIDTEL     = 200        ! max telescope ID in SURVEY.DEF
     &  ,MXIDFIELD   = 200        ! max FIELD ID in SURVEY.DEF
     &  ,MXFIELD_OVP = 12         ! max number of overlapping fields
     &  ,MXSEASON    = 100        ! max number of seasons
     &  ,MXSNHOST    = 2          ! max number of host matches to read/write
     &  ,MXVAR_PRIVATE = 40       ! max number of private variables
     &  ,MXCUT_PRIVATE = 10        ! max number of cuts on private var
     &  ,MXVAR_TERSE   = 30       ! max number of text-data  columms
     &  ,MXLISTNML   = 52    ! max list size for some NML lists
     &  ,MXFILT_OBS  = 62    ! max number of used observer filters
     &  ,MXFILT_ALL  = 80    ! max number of all possible filter defs
     &  ,MXFILT_REST = 20    ! max number of rest-frame filter types (was 12)
     &  ,MXFILT_KCOR = MXFILT_REST ! max number of filters used in Kcor table
     &  ,MXFILT_SNRMAX = 10       ! no more than this many SNRMAX(filt) cuts
     &  ,MNTYPE       =   1       ! min sn "type"
     &  ,MXTYPE       = 1000      ! max sn "type"
     &  ,MXZBIN_KCOR  = 100       ! max # Z-bins for KCOR tables
     &  ,MXTBIN_KCOR  = 150       ! max # Epochs for KCOR tables
     &  ,MXAVBIN_KCOR = 100       ! max # AV bins for KCOR tables
     &  ,MXCBIN_AVWARP = 100      ! max color index for AVWARP table
     &  ,MXTABLE1D_KCOR   = 8 000 000  ! max number of KCOR bins
     &  ,MXTABLE1D_AVWARP = 1 000 000  ! max number of bins for AVWARP table
     &  ,MXTABLE1D_LCMAG  = 2 000 000  ! idem for LCMAG  table
     &  ,MXTABLE1D_MWXT   = 2 000 000  ! idem for MWXT table
     &  ,KCORPACK      = 1000.  ! store KCOR * KCORPACK as I*2
     &  ,XTMW_FRACERR  = 0.16   ! error on MW Xtinc is 16% of XTMW
c
     &  ,HOFF_NOCUTS = 100
     &  ,HOFF_CUTS   = 200
     &  ,HOFF_SIM    = 20000
     &  ,LUNHIS    = 21  ! LUN for hbook output file
     &  ,LUNNML    = 22  ! LUN for input namelist
     &  ,LUNDAT    = 23  ! LUN for data
     &  ,LUNTMP    = 24
     &  ,LUNOUT    = 25
     &  ,LUNFIT    = 26
     &  ,LUNLIST   = 27
     &  ,LUNIGNORE = 28
     &  ,LUNDMP    = 29
     &  ,LUNRES1   = 31  ! for DMP_FITRES
     &  ,LUNRES2   = 32
     &  ,LUNINTERP = 33  ! for FLUX,MAGs interpolated at SN,MJD
     &  ,LUNSALT2  = 34  ! for SALT2 dictFile
     &  ,LUNPKMJD  = 35
     &  ,LUNCID    = 36  ! reserved for reading SNCID_LIST_FILE
     &  ,ISTAGE_INIT    = 10       ! init has finished
     &  ,ISTAGE_RDSN    = 20       ! SN have been read
     &  ,ISTAGE_CUTS    = 30       ! cuts have been applied
     &  ,ISTAGE_USRANA  = 40       ! USRANA has been called.
     &  ,ISTAGE_TEST    = 90       ! set for testing
     &  ,ZERO8          = 0.0
     &  ,ONE8           = 1.0
     &  ,ONE            = 1
     &  ,TEN8           = 10.0
     &  ,LOGTEN         = 2.302585092994
     &  ,INTERP_LINEAR  = 1    ! option flag to use linear DFINT
     &  ,INTERP_SMOOTH  = 2    ! option flag to use smoothing
     &  ,INTERP_ZSMOOTH = 3    ! linear DFINT, but smooth along z
c
     &  ,MXERRTYPE      = 10
     &  ,ERRTYPE_MINOS  = 1
     &  ,ERRTYPE_PARAB  = 2
     &  ,ERRTYPE_MARG   = 3  ! marginalized error
     &  ,ERRTYPE_BAD    = 6
C
     &  ,FCNFLAG_LAST     = 3    ! last MINUIT call
     &  ,FCNFLAG_USER     = 90   ! pass this IFLAG for user-FCNSNLC calls
     &  ,FCNFLAG_FAST     = 91   ! go as fast as possible (no LAST if-block)
     &  ,FCNFLAG_PRIOR_ONLY = 92
     &  ,FCNFLAG_SIGMA_ONLY = 93
     &  ,FCNFLAG_USESIM   = 99   ! pass this IFLAG to use SIM params
c
     &  ,PI     = 3.1415926535898
     &  ,PDFMIN      = 1.0E-5   ! used to speed up PDF integration
     &  ,PDFMAX_EDGE = 0.03     ! max allowed PDF value at edges
     &  ,PDFMIN_GOOD = 1.0E-4   ! PDF > PDFMIN_GOOD counts as good point
     &  ,OPT_INTEGPDF_QUITCHI2 = 2  ! abort FCNSNLC if chi2 > quitchi2
     &  ,OPT_INTEGPDF_FULL     = 1  ! do full FCNSNLC evaluation
c
     &  ,OPT_SNXT_CCM89   = 1  ! exact SN etinction using INIT_XTHOST
     &  ,OPT_SNXT_SJPAR   = 2  ! SN extinction with Jha's parameters
     &  ,OPT_KCORERR_SMOOTH = 1 ! use smooth half-Gaussian
     &  ,OPT_KCORERR_SJ     = 2   ! use Saurabh's Kcor error
     &  ,OPT_KCORERR_SJ5    = 5  ! x5 Saurabh's Kcor error
c
     &  ,RV_MWCOLORLAW_DEFAULT  = 3.1   ! A_V/E(B-V)
     &  ,OPT_MWCOLORLAW_DEFAULT = 94    ! ODonnel 94
     &  ,OPT_MWEBV_DEFAULT      =  1    ! whatever is in the data file
c
     &  ,MXLINE_ARGS      = 100
     &  ,ZEROPOINT_FLUXCAL_DEFAULT  = 27.5
     &  ,MXEPOCH_IGNORE   = 1000
     &  ,NPAR_ANYLC       = 8    ! for MNFIT_PKMJD
     &  ,MXPAR_SIMSED     = 100  ! max number of SIMSED parameters
     &  ,MXPAR_LCLIB      = 40   ! should be same as in genmag_LCLIB.h
     &  ,MXLAMBIN_SNSED   = 4000 ! Jan 2017: raised from 3000
     &  ,MXCUTBIT         = 64   ! max number of cut bits
c
     &  ,OPT_FILTUPD_EACHSN   = 1  ! default filter-updates
     &  ,OPT_FILTUPD_MAP      = 2  ! default filter-updates
     &  ,OPT_FILTUPD_SAMEFILT = 3  ! test: use same filter each SN
     &  ,OPT_FILTREST         = 1
     &  ,OPT_FILTOBS          = 2
     &  ,ISTAT_READAGAIN      = 7
     &  ,ISTAT_SKIP           = -1
c
     &  ,ITABLE_SNANA=1, ITABLE_FITRES=2, ITABLE_SNLCPAK=3
     &  ,ITABLE_SPECPAK=4, ITABLE_SPECTRA=5
     &  ,MXTABLE = 10
     &  ,IDTABLE_SNANA   = 7100
     &  ,IDTABLE_FITRES  = 7788
     &  ,IDTABLE_SPECTRA = 8000  ! SALT2 model spectra from LC fit
     &  ,IDTABLE_CIDPTR  = 1600
     &  ,IDTABLE_MCMC    = 7711
     &  ,OPT_PARSTORE_TEXTTABLE = 1  ! tag subset for TEXT table.
c
     &  ,MODEL_STRETCH    = 1
     &  ,MODEL_STRETCH2   = 2
     &  ,MODEL_MLCS2k2    = 3
     &  ,MODEL_SNOOPY     = 4
     &  ,MODEL_SALT2      = 6
     &  ,MODEL_SIMSED     = 7
     &  ,MODEL_BYOSED     = 8
     &  ,MODEL_NON1A      = 10
     &  ,MODEL_LCLIB      = 12
     &  ,MODEL_FIXMAG     = 20
     &  ,MXMODEL_INDEX    = 20
     &  ,NULLVAL          = -99999.
     &  ,IGNORE_HEADVAL   = -5555.
c
     &  ,MXCHAR_CCID       = 20   ! max len of CCID string (i.e, SN name)
     &  ,MXCHAR_VERSION    = 72   ! max len of VERSION_PHOTOMETRY
     &  ,MXCHAR_SURVEY     = 40   ! max len of SURVEY_NAME
     &  ,MXCHAR_PATH       = 160  ! max len of path
     &  ,MXCHAR_FILENAME   = 200  ! max len of filename with full path
     &  ,MXCHAR_MODELNAME  = 72   ! max len of model name
     &  ,MXCHAR_FIELDNAME  = 20   ! max len of field name
     &  ,MXCHAR_PARNAME    = 20   ! max len of parameter name
     &  ,MXCHAR_CUTNAME    = 160  ! to define cut names
     &  ,MXCHAR_FILEWORD   =  60  ! size of FILEWORD_LIST
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA    = 1    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_REJECT      = 2
     &  ,SNLCPAK_EPFLAG_CHI2        = 3
     &  ,SNLCPAK_EPFLAG_FITFUN      = 4
     &  ,SNLCPAK_EPFLAG_FLUXSIM     = 5    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_FLUXREST    = 6
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST = 7
     &  ,SNLCPAK_BANDFLAG_NDOF    = 100
     &  ,SNLCPAK_BANDFLAG_PKFLUX  = 101  ! filter-dependent
     &  ,SNLCPAK_BANDFLAG_PKMJD   = 102
     &  ,SNLCPAK_BANDFLAG_CHI2    = 103
c
     &  ,IFLAG_INI=1, IFLAG_ANA=2, IFLAG_END=3
     &  ,MAG_SATURATE = -7.0   ! for sim only
     &  ,CUTVAL_OPEN = 1.0E12  ! cutwin value to accept everything.
     &	,MASK_FLUXCOR_SNANA    = 1
     &  ,MASK_FLUXERRCOR_SNANA = 2
c
     &  ,EXIT_ERRCODE_SNANA = 21
     &  ,EXIT_ERRCODE_SNFIT = 22
     &  ,EXIT_ERRCODE_PSNID = 23
     &      )

c physical constants

      REAL*8
     &   PARSEC, CLIGHT, PEAKMAG_AT_10PC
     &  ,Zat10pc
     &  ,OMAT_DEFAULT, OMATERR_DEFAULT
     &  ,OLAM_DEFAULT, OLAMERR_DEFAULT
     &  ,ORAD_DEFAULT
     &  ,H0_DEFAULT,   H0ERR_DEFAULT
     &  ,W0_DEFAULT,   W0ERR_DEFAULT
     &  ,DWDA_DEFAULT, DWDAERR_DEFAULT

      PARAMETER (
     &   PARSEC        = 3.085678E13  ! 1 parsec (km)
     &  ,CLIGHT        = 2.998E5      ! c (km/sec)
c
     &  ,H0_DEFAULT        = 70.0 / ( 1.0E6 * Parsec )  ! standard value
     &  ,W0_DEFAULT        = -1.0
     &  ,DWDA_DEFAULT      =  0.0
     &  ,OMAT_DEFAULT      =  0.3
     &  ,OLAM_DEFAULT      =  0.7
     &  ,ORAD_DEFAULT      =  1.2E-5
c
     &  ,H0ERR_DEFAULT        =  7.0 / ( 1.0E6 * Parsec )
     &  ,W0ERR_DEFAULT        =  0.1
     &  ,DWDAERR_DEFAULT      =  0.0
     &  ,OMATERR_DEFAULT      =  0.03
     &  ,OLAMERR_DEFAULT      =  0.03
c
     &  ,PEAKMAG_AT_10PC   = -19.6
     &  ,Zat10pc           = 2.34E-9   ! magic redshift at 10 pc
     &  ,MJDOFF            = 0.0       ! 53000.
     &     )



C +CDE,SNFILECOM. inserted below by ypatchy. 

c Sep 9, 2010: Pulled out of SNDATACOM

c define HEADER MASK BITS for required variables.
      INTEGER
     &   HEADBIT_SNID,   HEADBIT_IAUC, HEADBIT_CIDSEL
     &  ,HEADBIT_SURVEY, HEADBIT_FILTERS
     &  ,HEADBIT_RA,    HEADBIT_DEC
     &  ,HEADBIT_MWEBV, HEADBIT_Z
     &  ,NBIT_HEADMASK
     &  ,HEADMASK_REQUIRED
c
     &  ,HEADMASK  ! reset for reach SN

      PARAMETER (
     &    HEADBIT_SNID    = 0  ! MASK = 1
     &   ,HEADBIT_IAUC    = 1  ! MASK = 2  ! added Dec 2015
     &   ,HEADBIT_CIDSEL  = 2  ! MASK = 4  ! added Dec 2015
     &   ,HEADBIT_SURVEY  = 3
     &   ,HEADBIT_FILTERS = 4
     &   ,HEADBIT_RA      = 5
     &   ,HEADBIT_DEC     = 6
     &   ,HEADBIT_MWEBV   = 7
     &   ,HEADBIT_Z       = 8
     &   ,NBIT_HEADMASK   = 9
     &   ,HEADMASK_REQUIRED  = 2**(NBIT_HEADMASK)-1 - 2 ! don't require IAUC
     &      )

      CHARACTER
     &   SNDATA_ROOT*(MXCHAR_PATH)
     &  ,SNANA_DIR*(MXCHAR_PATH)
     &  ,SNDATA_PATH*(MXCHAR_PATH)      ! subdir with data or sim files
     &  ,SNLIST_FILE*(MXCHAR_FILENAME)  ! input list of SNDATA Files
     &  ,SNREADME_FILE(MXVERS)*(MXCHAR_FILENAME)   ! name of EVERY README file
     &  ,SNDATA_FILE_CURRENT*(MXCHAR_FILENAME)     ! current file being read
     &  ,GLOBAL_BANNER*120
     &  ,SNDATA_PREFIX*(MXCHAR_FILENAME)  ! $SNDATA_ROOT/lcmerge/$VERSION
     &  ,C1ERR*88, C2ERR*88    ! generic error strings

      INTEGER
     &  ABSO_INDEX(MXSNLC) ! absolute (IFILE or IROW) vs. ISN index

      LOGICAL LFLAG_RDHEAD_ONLY

      COMMON / SNFILECOM /
     &   SNDATA_ROOT, SNANA_DIR, SNLIST_FILE
     &  ,SNDATA_FILE_CURRENT
     &  ,GLOBAL_BANNER, SNDATA_PREFIX, SNREADME_FILE, SNDATA_PATH
     &  ,C1ERR, C2ERR, ABSO_INDEX, HEADMASK, LFLAG_RDHEAD_ONLY



C +CDE,CTRLCOM. inserted below by ypatchy. 

c control variables and counters.

      INTEGER
     &   ISTAGE_SNANA       ! current stage of processing
     &  ,NACCEPT_CUT(MXCUTBIT) ! Number of SN that pass each cut
     &  ,NACCEPT_CID           ! # SN with valid CID
     &  ,NACCEPT_TYPE          ! # SN with valid TYPE
     &  ,NACCEPT_Z             ! idem with valid redshift
     &  ,NACCEPT_ZERR          ! idem with valid redshift error
     &  ,NMJD_IDTEL            ! # MJDs with valid telescope id
     &  ,APPLY_HEADER_CUTMASK  ! mask of applied header cuts
c
     &  ,JTIME_START
     &  ,JTIME_LOOPSTART        ! time at start of fits with TIME()
     &  ,JTIME_LOOPEND          ! time and end
     &  ,NCALL_SNANA_DRIVER
c
     &  ,NPASSCUT_INCREMENT(-1:MXTYPE,100)  ! 100 > NCUTBIT_SNLC
     &  ,NPASSCUT_FIT(-1:MXTYPE)
     &  ,NVAR_NEARNBR          ! Number of NEARNBR variables to analyze
     &  ,NSTORE_MAGCOR         ! number of stored MAGCOR values
     &  ,NUSE_MAGCOR           ! number of used MAGCOR values
     &  ,SIGN_MAGCOR           ! add or subtract
     &  ,FORCEMASK_FLUXCOR   ! mask to force fluxCor, even if already applied
     &  ,EXIT_ERRCODE        ! used for abort

      LOGICAL
     &   DO_FIT
     &  ,DO_FLUXERRCALC    ! T => compute error from PSF,SKY & ZPT
     &  ,LSIM_SNANA        ! simulated with SNANA
     &  ,LSIM_MAGOBS       ! data-like, but with SIM_MAGOBS
     &  ,ISJOB_SNANA          ! =T for snana.exe only.
     &  ,ISJOB_PSNID
     &  ,REFORMAT_SAVE_BADEPOCHS  ! set if bit2 is set on any OPT_REFORMAT
     &  ,STDOUT_UPDATE         ! T => update event to screen
     &  ,DOFUDGE_HOSTNOISE     ! T => FUDGE_HOSTNOISE_FILE is set
     &  ,DOFUDGE_NONLIN        ! T => NONLINEARITY_FILE is set
     &  ,DOFUDGE_FLUXERRMODEL  ! T => FLUXERRMODEL_FILE

      CHARACTER SNANA_VERSION*12

c global survey info
      CHARACTER
     &   SURVEY_NAME*(MXCHAR_SURVEY)
     &  ,SURVEY_NAME_LIST(MXSURVEY)*(MXCHAR_SURVEY) ! in SURVEY.DEF file
     &  ,SUBSURVEY_NAME*(MXCHAR_SURVEY)  ! e.g.,  SURVEY:  BLA(SUBSURVEY)
     &  ,SUBSURVEY_NAME_LIST*(MXCHAR_FILENAME) ! comma-sep list
     &  ,SURVEY_FILTERS*(MXFILT_ALL)  ! read from data file
     &  ,SURVEY_FIELDNAME(MXIDFIELD)*(MXCHAR_FIELDNAME) ! from SURVEY.DEF

      INTEGER
     &   NFIELD_SURVEY             ! number of survey fields in SURVEY.DEF
     &  ,SURVEY_IDFIELD(MXIDFIELD) ! integer ID for each field

      REAL
     &   ZEROPOINT_FLUXCAL(MXFILT_OBS)  ! defines calibrated flux

      INTEGER*4
     &   N_VERSION             ! number of photometry version to read
     &  ,N_SNLC                ! number of SN lightcurves read
     &  ,N_SNLC_CUTS           ! Number of SN after cuts (bookkeeping only)
     &  ,N_SNLC_FIT            ! Number of fitted SN
     &  ,N_SNLC_FITCUTS        ! Number of fitted SN after fit cuts
     &  ,N_SNLC_COVFIX         ! Number of SN with fixed COV to be invertible
     &  ,N_DUPLICATE_CID       ! Number of duplicate CIDs
     &  ,N_DUPLICATE_MJD       ! Number of duplicate MJD+BAND (Jun 2017)
     &  ,NSTORE_DUPLICATE_MJD  ! Number stored
     &  ,N_SNFILE              ! # of SN files to read per version
     &  ,N_SNFILE_LAST         ! idem, as of last version
     &  ,NEPOCH_TOT            ! total number of epochs read
     &  ,NEPOCH_CUT            ! total number of epochs passing cuts
     &  ,NEPOCH_BADPHOT        ! # epochs with bad PHOTFLAG (per event)
     &  ,NEPOCH_BADPHOT_SUM    ! # epochs with bad PHOTFLAG (summed)
     &  ,IDSURVEY              ! survey ID from SURVEY.DEF
     &  ,IDSUBSURVEY           ! =IDSURVEY unless subSurvey is different
     &  ,IDSURVEY_LIST(MXSURVEY) ! corresponds to SURVEY_NAME_LIST
     &  ,NSURVEY_LIST          ! size of SURVEY_NAME_LIST & IDSURVEY_LIST

      LOGICAL*1
     &    EXIST_FILT(MXFILT_OBS)  ! T => at least one point per filt
     &   ,FOUND_SURVEY
     &   ,FORMAT_TEXT    ! ascii/txt for input data
     &   ,FORMAT_TERSE   ! terse-text for input data
     &   ,FORMAT_VERBOSE ! verbose-text for ...
     &   ,FORMAT_FITS    ! snfitsio for input data.

      INTEGER
     &   NVAR_TERSE          ! Number variables for terse format
     &  ,ITERSE_MJD          ! index of MJD in TERSE variable list
     &  ,ITERSE_FILTER       ! idem for filter
     &  ,ITERSE_FIELD
     &  ,ITERSE_FLUXCAL
     &  ,ITERSE_FLUXCALERR
     &  ,ITERSE_PHOTFLAG
     &  ,ITERSE_PHOTPROB
     &  ,ITERSE_ZPFLUX     ! mag <-> native flux zeropoint
     &  ,ITERSE_PSFSIG     ! added July 18, 2013
     &  ,ITERSE_SKYSIG     ! idem
     &  ,ITERSE_SKYSIG_T   ! added Aug 7 2014
     &  ,ITERSE_XPIX       ! added Aug 7 2014
     &  ,ITERSE_YPIX       ! added Aug 7 2014
     &  ,ITERSE_GAIN       ! added Oct 2015 to compute FLUXCAL_ERRCALC
     &  ,ITERSE_MAGOBS
     &  ,ITERSE_CCDNUM     ! added Oct 16 2017
     &  ,ITERSE_SIM_EPFILTREST
     &  ,ITERSE_SIM_EPMAGOBS
     &  ,ITERSE_SIM_EPMAGREST

      INTEGER N_SNLC_PLOT
      LOGICAL MADE_LCPLOT  ! SAVE:  T if LC plot was made

      CHARACTER VARLIST_TERSE(MXLISTNML)*(MXCHAR_CCID) ! TERSE variable names

      COMMON / CTRLCOM /
     &     SNANA_VERSION, SURVEY_NAME, SURVEY_NAME_LIST
     &    ,SUBSURVEY_NAME, SUBSURVEY_NAME_LIST
     &    ,NSURVEY_LIST, IDSURVEY, IDSUBSURVEY, IDSURVEY_LIST
     &    ,SURVEY_FILTERS
     &    ,SURVEY_FIELDNAME, SURVEY_IDFIELD, NFIELD_SURVEY
     &    ,ISJOB_SNANA, ISJOB_PSNID, ZEROPOINT_FLUXCAL
     &    ,NACCEPT_CUT, NACCEPT_CID, NACCEPT_TYPE
     &    ,NACCEPT_Z, NACCEPT_ZERR, NMJD_IDTEL
     &    ,DO_FIT, DO_FLUXERRCALC, NVAR_NEARNBR
     &    ,NSTORE_MAGCOR, NUSE_MAGCOR, SIGN_MAGCOR, FORCEMASK_FLUXCOR
     &    ,LSIM_SNANA, LSIM_MAGOBS, APPLY_HEADER_CUTMASK
     &    ,ISTAGE_SNANA
     &    ,N_VERSION, N_SNLC, N_SNLC_CUTS, N_SNLC_FIT, N_SNLC_FITCUTS
     &    ,N_SNLC_COVFIX, N_SNFILE, N_SNFILE_LAST, N_DUPLICATE_CID
     &    ,N_DUPLICATE_MJD, NSTORE_DUPLICATE_MJD
     &    ,NEPOCH_TOT, NEPOCH_CUT, NEPOCH_BADPHOT, NEPOCH_BADPHOT_SUM
     &    ,REFORMAT_SAVE_BADEPOCHS, STDOUT_UPDATE
     &    ,DOFUDGE_HOSTNOISE, DOFUDGE_NONLIN, DOFUDGE_FLUXERRMODEL
     &    ,NVAR_TERSE, VARLIST_TERSE
     &    ,ITERSE_MJD, ITERSE_FILTER, ITERSE_FIELD
     &    ,ITERSE_FLUXCAL, ITERSE_FLUXCALERR
     &    ,ITERSE_PHOTFLAG, ITERSE_PHOTPROB
     &    ,ITERSE_ZPFLUX, ITERSE_PSFSIG, ITERSE_SKYSIG, ITERSE_SKYSIG_T
     &    ,ITERSE_MAGOBS, ITERSE_XPIX, ITERSE_YPIX, ITERSE_GAIN
     &    ,ITERSE_CCDNUM, ITERSE_SIM_EPFILTREST
     &    ,ITERSE_SIM_EPMAGOBS, ITERSE_SIM_EPMAGREST
     &    ,JTIME_START, JTIME_LOOPSTART, JTIME_LOOPEND
     &    ,NCALL_SNANA_DRIVER, NPASSCUT_INCREMENT, NPASSCUT_FIT
     &    ,N_SNLC_PLOT, MADE_LCPLOT
     &    ,EXIT_ERRCODE

c logical *1 stuff

      COMMON / SNDATCOM1 / EXIST_FILT, FOUND_SURVEY
     &    ,FORMAT_TEXT, FORMAT_TERSE, FORMAT_VERBOSE
     &    ,FORMAT_FITS

c SNTABLE control variables (May 2014)
      LOGICAL
     &   USE_TABLEFILE_HBOOK  ! write table in HBOOK format
     &  ,USE_TABLEFILE_ROOT   ! write table in ROOT format
     &  ,USE_TABLEFILE_TEXT   ! write table in TEXT format
     &  ,USE_TABLEFILE        ! T if either any of the above are set
     &  ,WRTABLEFILE_IAUC     ! T-> SNID -> IAUC (for writing tables)
     &  ,WRTABLEFILE_SIMVAR   ! T-> include SIM_XXX vars for simulation
     &  ,WRTABLEFILE_ZPHOT    ! T-> include ZPHOT info for PHOTOZ fit


      INTEGER
     &   OPT_TABLE(MXTABLE)       ! option(s) for each table
     &  ,CUTMASK_SNANA_TABLE      ! select CUTFLAG_SNANA for SNANA table

      REAL*8
     &   PRESCALE_TABLE(MXTABLE)  ! table pre-scale (to reduce size)

      CHARACTER
     &   TEXTFORMAT_TABLE(MXTABLE)*8

c arrays for TABLE_FILTER_REMAP (Feb 2017)
      INTEGER NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE(MXFILT_ALL)
      CHARACTER FILTLIST_REMAP_TABLE*(MXFILT_ALL)

c Sep 23 2017: allow user codes to add private variables to SNTABLE
      INTEGER   NTABLEVAR_USER  ! e.g., set in snana_private.cra
      CHARACTER TABLEVARNAME_USER(MXVAR_PRIVATE)*40
      REAL      TABLEVALUE_USER(MXVAR_PRIVATE)

      COMMON / SNTABLECOM /
     &     USE_TABLEFILE_HBOOK, USE_TABLEFILE_ROOT
     &    ,USE_TABLEFILE_TEXT, USE_TABLEFILE
     &    ,OPT_TABLE, TEXTFORMAT_TABLE
     &    ,CUTMASK_SNANA_TABLE, WRTABLEFILE_IAUC, WRTABLEFILE_SIMVAR
     &    ,WRTABLEFILE_ZPHOT
c
     &    ,NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE
     &    ,FILTLIST_REMAP_TABLE
     &    ,NTABLEVAR_USER, TABLEVARNAME_USER
     &    ,TABLEVALUE_USER

      COMMON / SNTABLECOM8 / PRESCALE_TABLE

c May 6, 2008: define lists for epochs to ignore

      INTEGER NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
      CHARACTER
     &   EPOCH_IGNORE_CCID(MXEPOCH_IGNORE)*(MXCHAR_CCID)
     &  ,EPOCH_IGNORE_FILT(MXEPOCH_IGNORE)*4
     &  ,EPOCH_IGNORE_LASTFILE*(MXCHAR_FILENAME)

      REAL*8
     &   EPOCH_IGNORE_MJD(MXEPOCH_IGNORE)

      COMMON / EPIGNORE_COM / NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
     &  ,EPOCH_IGNORE_CCID, EPOCH_IGNORE_FILT, EPOCH_IGNORE_LASTFILE
      COMMON / EPIGNORE_COM8 / EPOCH_IGNORE_MJD


      REAL*8  DUPLICATE_MJDLIST(200)
      COMMON / DUPMJDCOM / DUPLICATE_MJDLIST


C +CDE,SNLCCOM. inserted below by ypatchy. 

c Main common block variables for light curves and
c analysis-related variables.

      INTEGER
     &   SNLC_CID               ! integer CID
     &  ,ISNLC_SNRECON_USE(MXEPOCH)  ! 1 -> epoch used in SNRECON
     &  ,ISNLC_PHOTFLAG(MXEPOCH) ! photomety flags
     &  ,ISNLC_LENCCID     ! char-len of CCID
     &  ,ISNLC_LENIAUC     ! char-len of IAUC name
     &  ,ISNLC_VERSION     ! photometry version index vs. ISN
     &  ,ISNLC_TYPE        ! type (120=confirmed Ia, etc ...)
     &  ,ISNLC_ABSO_INDEX  ! absolute index (row or file number)
     &  ,ISNLC_IFILE       ! ifile index, text format only
     &  ,ISNLC_IFILT_OBS(MXEPOCH) ! filt-index for each epoch/SN
     &  ,ISNLC_NEWMJD_HEAD   ! Number of NEWMJDs in header
     &  ,ISNLC_NEWMJD_FOUND  ! Number of NEWMJDs found in file
     &  ,ISNLC_NEWMJD_STORE  ! Number of NEWMJDs stored
     &  ,ISNLC_NEWMJD_CUTS   ! Number of NEWMJDs after cuts
     &  ,ISNLC_EPOCH_RANGE_NEWMJD(2,MXEPOCH)
     &  ,ISNLC_NFILT_NEWMJD(MXEPOCH)   ! # observed filters per NEWMJD
     &  ,ISNLC_NFILT_SNRMAX      ! number of filters passing SNR cut
     &  ,ISNLC_NFILT_SNRMAX2     ! idem for 2nd SNRMAX cut
     &  ,ISNLC_NFILT_TRESTMIN     ! Nfilt passing TRESTMIN cut
     &  ,ISNLC_NFILT_TRESTMAX     ! idem for TRESTMAX
     &  ,ISNLC_NFILT_TREST2       ! Nfilt passing CUTWIN_TREST2 cut
     &  ,ISNLC_NFILT_THRESH(MXEPOCH)
     &  ,ISNLC_NEPOCH_FOUND  ! actual number of epochs found
     &  ,ISNLC_NEPOCH_STORE  ! number of epochs stored in memory
     &  ,ISNLC_NEPOCH_USE       ! used after PHOTMASK cuts
     &  ,ISNLC_NEPOCH_PHOTPROB  ! NEPOCH with PHOTPROB >= 0
     &  ,ISNLC_NEPOCH_DETECT    ! NEPOCH with detection
     &  ,ISNLC_NEPOCH_FILT(MXFILT_OBS)  ! NEPOCH vs. filter
     &  ,ISNLC_NEPOCH_PRESN(MXFILT_OBS) ! pre-explosion epochs
     &  ,ISNLC_NMJD_INCLUDE             ! NOBS in CUTWIN_MJD_INCLUDE window
c
     &  ,ISNLC_FAKE             ! => real data, else it's a fake
     &  ,ISNLC_CCDNUM(MXEPOCH)  ! read from header (May 2017)
     &  ,ISNLC_IDTEL(MXEPOCH)   ! integer telescope id
     &  ,ISNLC_IDFIELD(MXEPOCH) ! integer field id
     &  ,ISNLC_NFIELD_OVP       ! number of fields (>=2 for overlap)
     &  ,ISNLC_CUTFLAG_REQEP        ! idem
     &  ,ISNLC_CUTFLAG_PRIVATE      ! idem for private var cuts
     &  ,ISNLC_CUTFLAG_USRCUTS      ! idem for USRCUTS
     &  ,ISNLC_CUTFLAG_SIMVAR       ! idem for SIMVAR cuts
     &  ,ISNLC_WRMASK_FLUXCOR_SNANA   ! write if fudges applied to data
     &  ,ISNLC_RDMASK_FLUXCOR_SNANA   ! read if SNANA fudges applied to data

      REAL*8
     &   SNLC8_RA        ! RA   vs. SN
     &  ,SNLC8_DEC       ! DEC  vs. SN
     &  ,SNLC8_MJD(MXEPOCH)  ! MJD for each epoch
     &  ,SNLC8_MJDMIN        !   min MJD among all measurements
     &  ,SNLC8_MJDMAX        !   max MJD ...

      REAL
     &   SNLC_ZHELIO        ! redshift used for final analysis
     &  ,SNLC_ZHELIO_ERR    ! error on above
     &  ,SNLC_ZCMB          ! redshift used for final analysis
     &  ,SNLC_ZCMB_ERR      ! error on above
     &  ,SNLC_ZHD           ! redshift used for Hubble diagram
     &  ,SNLC_ZHD_ERR       ! error on above
     &  ,SNLC_ZSN           ! redshift of SN only (ignoring host-z)
     &  ,SNLC_ZSN_ERR       ! error on above
     &  ,SNLC_REDSHIFT      ! redshift used for light curve fit
     &  ,SNLC_REDSHIFT_ERR  ! error on above
     &  ,SNLC_VPEC          ! pec. velocity
     &  ,SNLC_VPEC_ERR      ! error on above
     &  ,SNLC_ZPEC
     &  ,SNLC_ZPEC_ERR
     &  ,SNLC_Trestmin  ! earliest epoch, rest frame days since peak
     &  ,SNLC_Trestmax  ! latest   epoch, rest frame days since peak
     &  ,SNLC_TrestRange
     &  ,SNLC_Tobsmin   ! earliest epoch, obs frame days since peak
     &  ,SNLC_Tobsmax   ! latest   epoch, obs frame days since peak
     &  ,SNLC_TGAPMAX     ! max gap within TREST-range
     &  ,SNLC_T0GAPMAX    ! max gap near peak
     &  ,SNLC_SNRMAX_FILT(0:MXFILT_OBS)   ! max S/N per filter/SN
     &  ,SNLC_SNRMAX_SORT(MXFILT_OBS)     ! 1st, 2nd ... SNRMAX by filt
     &  ,SNLC_FLUXCALMAX(MXFILT_OBS)    ! max flux per filter/SN
     &  ,SNLC_SNANACALC_PEAKMJD     ! SNANA-estimate of PEAKKMJD
     &  ,SNLC_SNANACALC_PEAKMJD_FITPAR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_SNANACALC_PEAKMJD_FITERR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_PHOTPROB(MXEPOCH)   ! generic 'fit probability' per epoch
     &  ,SNLC_PHOTPROB_MIN        ! min photprob for PHOTPROB>0
     &  ,SNLC_AIRMASS(MXEPOCH)
     &  ,SNLC_CLOUDCAM_AVG(MXEPOCH)
     &  ,SNLC_CLOUDCAM_SIG(MXEPOCH)
     &  ,SNLC_MOONDIST(MXEPOCH)    ! in deg.
     &  ,SNLC_MOONPHASE(MXEPOCH)   ! 0 - 1
     &  ,SNLC_TREST(MXEPOCH)          !MJD-SearhMJDPEAK / (1+z)
     &  ,SNLC_GAIN(MXEPOCH)      ! e/AUD
     &  ,SNLC_RDNOISE(MXEPOCH)   ! read noise per pix, e-
     &  ,SNLC_PIXSIZE            ! pixel size
     &  ,SNLC_NXPIX              ! total number of X-pixels (Aug 7 2014)
     &  ,SNLC_NYPIX              ! total number of Y-pixels
     &  ,SNLC_XPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_YPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_AREAFRAC(MXEPOCH)  ! area-frac contained by XPIX,YPIX
     &  ,SNLC_MWEBV              ! Milky Way Galactic E(B-V)
     &  ,SNLC_MWEBV_ERR          ! error on above
     &  ,SNLC_SKY_AVG(MXEPOCH)   ! avg sky/pix in field, ADU
     &  ,SNLC_SKYSIG(MXEPOCH)    ! sigma on above
     &  ,SNLC_SKYSIG_T(MXEPOCH)  ! sigma on template run
     &  ,SNLC_PSF_SIG1(MXEPOCH)  ! sigma, pixels
     &  ,SNLC_PSF_SIG2(MXEPOCH)
     &  ,SNLC_PSF_RATIO(MXEPOCH)
     &  ,SNLC_AREA_NOISE(MXEPOCH)  ! Noise-equivalent area (pixels)
     &  ,SNLC_FLUX(MXEPOCH)
     &  ,SNLC_FLUX_ERRTOT(MXEPOCH)
     &  ,SNLC_FLUX_NSIG(MXEPOCH)   ! significance
     &  ,SNLC_FLUXCAL_OFF(MXFILT_OBS)    ! add SN light from template
     &  ,SNLC_FLUXCAL_ERRCALC(MXEPOCH)   ! calc a-la simulation
     &  ,SNLC_FLUXCAL_HOSTERRCALC(MXEPOCH)   ! idem for host error
     &  ,SNLC_FLUXCAL(MXEPOCH)
     &  ,SNLC_FLUXCAL_ERRTOT(MXEPOCH)
     &  ,SNLC_MAG(MXEPOCH)
     &  ,SNLC_MAG_ERRPLUS(MXEPOCH)
     &  ,SNLC_MAG_ERRMINUS(MXEPOCH)
     &  ,SNLC_ZEROPT(MXEPOCH)
     &  ,SNLC_ZEROPT_ERR(MXEPOCH)
     &  ,SNLC_ZEROPT_forCUT(MXEPOCH)  ! depends on CUTWIN_ZPADU or CUTWIN_ZPNPE
     &  ,SNLC_DLMAG                   ! 5*log10(10pc/DL)
     &  ,SNLC_SKYFLUXCAL(MXEPOCH)     ! calculated sky fluxcal/pixel
     &  ,SNLC_MWXT_MAG(MXFILT_OBS)        ! mag stellar extinct
     &  ,SNLC_MWXT_FLUXFRAC(MXFILT_OBS)   ! same for flux (< 1)
     &  ,SNLC_MWXT_MAGERR(MXFILT_OBS)     ! Galactic mag err per filter
     &  ,SNLC_SEARCH_PEAKMJD            ! external PEAKMJD
     &  ,SNLC_DTOBS(MXEPOCH)            ! time since last obs
     &  ,SNLC_DTOBS_SAMEFILT(MXEPOCH)   ! idem, but same filter
     &  ,SNLC_TLIVE_DETECT     ! MJD(last detection) - MJD(1st detection)

      INTEGER
     &   NSEASON_TOT     ! total number of seasons
     &  ,NSEASON_ACTIVE  ! total number of active seasons

      REAL*4
     &   MULTISEASON_CHI2RED(MXSEASON) ! borrow another 'MX' param
     &  ,MULTISEASON_AVGFLUX(MXSEASON)
     &  ,MULTISEASON_MJDMIN(MXSEASON)
     &  ,MULTISEASON_MJDMAX(MXSEASON)

      CHARACTER
     &   SNLC_CCID*(MXCHAR_CCID)      ! char CID
     &  ,SNLC_IAUC*(MXCHAR_CCID)
     &  ,SNLC_CASTCID*8             ! = 'INT' or 'CHAR' to indicate cast
     &  ,SNLC_FIELD(MXEPOCH)*(MXCHAR_FIELDNAME)    ! char field name
     &  ,SNLC_FIELD_OVPLIST(MXFIELD_OVP)*(MXCHAR_FIELDNAME)  ! sparse list of overlap fields
     &  ,SNLC_FIELDLIST*(MXCHAR_FIELDNAME) ! e.g., '82S', '82N+82S'

      COMMON / SNLC8COM /
     &    SNLC8_RA, SNLC8_DEC
     &  , SNLC8_MJD, SNLC8_MJDMIN, SNLC8_MJDMAX

      COMMON / SNLCCOM /
     &     SNLC_CID, SNLC_CCID,  SNLC_CASTCID, SNLC_IAUC
     &    ,SNLC_FIELD, SNLC_FIELD_OVPLIST, SNLC_FIELDLIST
     &    ,SNLC_SNANACALC_PEAKMJD
     &    ,SNLC_SNANACALC_PEAKMJD_FITPAR
     &    ,SNLC_SNANACALC_PEAKMJD_FITERR
     &    ,SNLC_SEARCH_PEAKMJD
     &    ,SNLC_ZHELIO, SNLC_ZHELIO_ERR
     &    ,SNLC_ZCMB,   SNLC_ZCMB_ERR, SNLC_ZHD, SNLC_ZHD_ERR
     &    ,SNLC_ZSN, SNLC_ZSN_ERR
     &    ,SNLC_REDSHIFT, SNLC_REDSHIFT_ERR
     &    ,SNLC_VPEC, SNLC_VPEC_ERR, SNLC_ZPEC, SNLC_ZPEC_ERR
     &    ,SNLC_TRESTMIN, SNLC_TRESTMAX, SNLC_TrestRange, SNLC_TREST
     &    ,SNLC_TGAPMAX,  SNLC_T0GAPMAX
     &    ,SNLC_Tobsmin, SNLC_Tobsmax
     &    ,SNLC_SNRMAX_FILT, SNLC_SNRMAX_SORT
     &    ,SNLC_FLUXCALMAX, SNLC_PHOTPROB, SNLC_PHOTPROB_MIN
     &    ,SNLC_CLOUDCAM_AVG, SNLC_CLOUDCAM_SIG
     &    ,SNLC_MOONDIST, SNLC_MOONPHASE, SNLC_AIRMASS
     &    ,SNLC_GAIN, SNLC_RDNOISE
     &    ,SNLC_MWXT_MAG, SNLC_MWXT_FLUXFRAC, SNLC_MWXT_MAGERR
     &    ,SNLC_MWEBV, SNLC_MWEBV_ERR
     &    ,SNLC_XPIX, SNLC_YPIX, SNLC_AREAFRAC, SNLC_PIXSIZE
     &    ,SNLC_NXPIX, SNLC_NYPIX
     &    ,SNLC_SKY_AVG,  SNLC_SKYSIG,  SNLC_SKYSIG_T
     &    ,SNLC_PSF_SIG1, SNLC_PSF_SIG2, SNLC_PSF_RATIO
     &    ,SNLC_AREA_NOISE
c
     &    ,SNLC_FLUX, SNLC_FLUX_ERRTOT
     &    ,SNLC_FLUX_NSIG, SNLC_FLUXCAL_OFF
     &    ,SNLC_FLUXCAL_ERRCALC, SNLC_FLUXCAL_HOSTERRCALC
     &    ,SNLC_FLUXCAL, SNLC_FLUXCAL_ERRTOT
     &    ,SNLC_MAG, SNLC_MAG_ERRPLUS, SNLC_MAG_ERRMINUS
     &    ,SNLC_ZEROPT, SNLC_ZEROPT_ERR, SNLC_ZEROPT_forCUT
     &    ,SNLC_SKYFLUXCAL, SNLC_DLMAG
     &    ,SNLC_DTOBS, SNLC_DTOBS_SAMEFILT

      COMMON / ISNCOM /
     &     ISNLC_VERSION, ISNLC_SNRECON_USE, ISNLC_PHOTFLAG
     &    ,ISNLC_NFILT_THRESH, ISNLC_TYPE
     &    ,ISNLC_ABSO_INDEX, ISNLC_IFILE,ISNLC_LENCCID,ISNLC_LENIAUC
     &    ,ISNLC_NFILT_NEWMJD, ISNLC_IFILT_OBS
     &    ,ISNLC_NEWMJD_HEAD, ISNLC_NEWMJD_FOUND, ISNLC_NEWMJD_STORE
     &    ,ISNLC_NEWMJD_CUTS, ISNLC_EPOCH_RANGE_NEWMJD
     &    ,ISNLC_NEPOCH_FOUND, ISNLC_NEPOCH_STORE
     &    ,ISNLC_NEPOCH_USE, ISNLC_NEPOCH_PHOTPROB
     &    ,ISNLC_NEPOCH_FILT,  ISNLC_NEPOCH_PRESN, ISNLC_NMJD_INCLUDE
     &    ,ISNLC_NEPOCH_DETECT, SNLC_TLIVE_DETECT
     &    ,ISNLC_FAKE, ISNLC_CCDNUM
     &    ,ISNLC_NFILT_SNRMAX, ISNLC_NFILT_SNRMAX2
     &    ,ISNLC_NFILT_TRESTMIN, ISNLC_NFILT_TRESTMAX
     &    ,ISNLC_NFILT_TREST2
     &    ,ISNLC_IDTEL, ISNLC_IDFIELD, ISNLC_NFIELD_OVP
     &    ,ISNLC_CUTFLAG_REQEP
     &    ,ISNLC_CUTFLAG_PRIVATE, ISNLC_CUTFLAG_USRCUTS
     &    ,ISNLC_CUTFLAG_SIMVAR
     &    ,ISNLC_WRMASK_FLUXCOR_SNANA, ISNLC_RDMASK_FLUXCOR_SNANA

c - - - - -  ADDCOL stuff - - - - - - -
c ADDCOL arrays are loaded with original filter indices,
c or with REMAPed filter indices. Allows using unique
c pointer for calls to SNTABLE_ADDCOL.

      CHARACTER
     &    ADDCOL_FILTERS*(MXFILT_ALL) ! SURVEY_FILTERS

      REAL
     &    ADDCOL_SNHOST_MAGOBS(MXFILT_ALL,MXSNHOST)
     &   ,ADDCOL_SNHOST_SBFLUXCAL(MXFILT_ALL)
     &   ,ADDCOL_FLUXCALMAX(MXFILT_ALL)
     &   ,ADDCOL_CHI2_FITPKMJD(MXFILT_ALL)
     &   ,ADDCOL_SNRMAX(MXFILT_ALL)
     &   ,ADDCOL_XTMW(MXFILT_ALL)

      INTEGER
     &    ADDCOL_NDOF_FITPKMJD(MXFILT_ALL)

      COMMON / ADDCOL_FILTERCOM /
     &   ADDCOL_FILTERS
     &  ,ADDCOL_SNHOST_MAGOBS, ADDCOL_SNHOST_SBFLUXCAL
     &  ,ADDCOL_FLUXCALMAX, ADDCOL_CHI2_FITPKMJD, ADDCOL_NDOF_FITPKMJD
     &  ,ADDCOL_SNRMAX, ADDCOL_XTMW


c - - - -- MULTI-SEASON STUFF - - - -
      COMMON / MULTISEASONCOM / NSEASON_TOT, NSEASON_ACTIVE
     &  ,MULTISEASON_CHI2RED, MULTISEASON_AVGFLUX
     &  ,MULTISEASON_MJDMIN,  MULTISEASON_MJDMAX


C +CDE,SNCUTS. inserted below by ypatchy. 

c Define cut-window for each variables
c and cut-mask for each epoch.
c User must fill cutwin_XXX arrays in main routine.

      INTEGER
     &   CUTBIT_CID, CUTBIT_SNTYPE
     &  ,CUTBIT_REDSHIFT, CUTBIT_REDSHIFT_ERR
     &  ,CUTBIT_RA,  CUTBIT_DEC
     &  ,CUTBIT_HOSTSEP
     &  ,CUTBIT_Trestmin, CUTBIT_Trestmax, CUTBIT_TrestRange
     &  ,CUTBIT_TREST2
     &  ,CUTBIT_Tgapmax,  CUTBIT_T0gapmax
     &  ,CUTBIT_Tobsmin,  CUTBIT_Tobsmax
     &  ,CUTBIT_SEARCH           ! sim only
     &  ,CUTBIT_TEMPLATE_INDEX   ! sim only (Dec 2018)
     &  ,CUTBIT_PEAKMJD
     &  ,CUTBIT_NMJD_INCLUDE    ! Aug 11 2015
     &  ,CUTBIT_NEPOCH
     &  ,CUTBIT_SNRMAX         ! global SNRMAX cut; any filter(s)
     &  ,CUTBIT_SNRMAX2        ! 2nd global SNRMAX cut; any filter(s)
     &  ,CUTBIT_NFIELD         ! cut on number of fields used by SN
     &  ,CUTBIT_MWEBV          ! cut on Galactic extinct (May 8 2012)
     &  ,CUTBIT_NSEASON_ACTIVE ! Nseason with activity
     &  ,CUTBIT_REQEP          ! required epochs (to select short events)
     &  ,CUTBIT_PRIVATE        ! all private-var cuts together (Nov 3 2014)
     &  ,CUTBIT_SIMVAR         ! SIMVAR cuts
     &  ,CUTBIT_USRCUTS        ! user cuts from private code (Aug 2015)
     &  ,CUTBIT_OFFSET_SBFLUX  ! host surface brightness
     &  ,CUTBIT_OFFSET_SNRMAX  ! SNRMAX cuts in each filter
     &  ,CUTBIT_MJD_MARKER     ! things above do NOT depend on epoch
     &  ,CUTBIT_PSF,  CUTBIT_ZP
     &  ,CUTBIT_MOONDIST,  CUTBIT_PHOTPROB, CUTBIT_AIRMASS
     &  ,CUTBIT_TREST, CUTBIT_MJD
     &  ,CUTBIT_IDTEL          ! valid telescope
     &  ,CUTBIT_NBAND_THRESH
     &  ,CUTBIT_NFILT_SNRMAX
     &  ,CUTBIT_NFILT_SNRMAX2
     &  ,CUTBIT_NFILT_TRESTMIN
     &  ,CUTBIT_NFILT_TRESTMAX
     &  ,CUTBIT_NFILT_TREST2
     &  ,CUTBIT_OFFSET
     &  ,NCUTBIT
     &  ,NCUTBIT_SNLC

      PARAMETER (
     &   CUTBIT_CID            = 1
     &  ,CUTBIT_SNTYPE         = 2
     &  ,CUTBIT_REDSHIFT       = 3
     &  ,CUTBIT_REDSHIFT_ERR   = 4
     &  ,CUTBIT_RA             = 5
     &  ,CUTBIT_DEC            = 6
     &  ,CUTBIT_HOSTSEP        = 7
     &  ,CUTBIT_TrestMIN       = 8
     &  ,CUTBIT_TrestMAX       = 9
     &  ,CUTBIT_TrestRange     = 10  ! Dec 2017 (TrestMax-TrestMin)
     &  ,CUTBIT_Trest2         = 11  ! Oct 2012
     &  ,CUTBIT_Tgapmax        = 12
     &  ,CUTBIT_T0gapMAX       = 13
     &  ,CUTBIT_TobsMIN        = 14
     &  ,CUTBIT_TobsMAX        = 15
     &  ,CUTBIT_PEAKMJD        = 16
     &  ,CUTBIT_NMJD_INCLUDE   = 17
     &  ,CUTBIT_NEPOCH         = 18  ! total # of measurements
     &  ,CUTBIT_SEARCH         = 19  ! found by search (SIM only)
     &  ,CUTBIT_NFILT_SNRMAX   = 20
     &  ,CUTBIT_NFILT_SNRMAX2  = 21
     &  ,CUTBIT_NFILT_TRESTMIN = 22
     &  ,CUTBIT_NFILT_TRESTMAX = 23
     &  ,CUTBIT_NFILT_TREST2   = 24
     &  ,CUTBIT_SNRMAX         = 25
     &  ,CUTBIT_SNRMAX2        = 26
     &  ,CUTBIT_NFIELD         = 27
     &  ,CUTBIT_MWEBV          = 28  ! May 2012
     &  ,CUTBIT_NSEASON_ACTIVE = 29  ! May 2019
     &  ,CUTBIT_REQEP          = 30  ! required epochs (9/17/2017)
     &  ,CUTBIT_PRIVATE        = 31  ! all private-var cuts
     &  ,CUTBIT_SIMVAR        = 32  !
     &  ,CUTBIT_USRCUTS        = 33  ! PASS_USRCUTS=T
     &  ,CUTBIT_OFFSET_SBFLUX  = 34  ! Mar 2012
     &  ,CUTBIT_OFFSET_SNRMAX  = 34 + MXFILT_SNRMAX
     &  ,CUTBIT_MJD_MARKER     = CUTBIT_OFFSET_SNRMAX + MXFILT_SNRMAX
     &  ,NCUTBIT_SNLC          = CUTBIT_MJD_MARKER
c SN-dependent cuts above
c MJD-dependent cuts below
     &  ,CUTBIT_PSF          = CUTBIT_MJD_MARKER + 1  ! PSF cut
     &  ,CUTBIT_ZP           = CUTBIT_MJD_MARKER + 2  ! ZP  cut (Feb 2018)
     &  ,CUTBIT_MOONDIST     = CUTBIT_MJD_MARKER + 3  ! degrees
     &  ,CUTBIT_PHOTPROB     = CUTBIT_MJD_MARKER + 4  ! Oct 23 2018
     &  ,CUTBIT_AIRMASS      = CUTBIT_MJD_MARKER + 5  ! 1=overhead
     &  ,CUTBIT_NBAND_THRESH = CUTBIT_MJD_MARKER + 6
     &  ,CUTBIT_Trest        = CUTBIT_MJD_MARKER + 7
     &  ,CUTBIT_MJD          = CUTBIT_MJD_MARKER + 8
     &  ,CUTBIT_IDTEL        = CUTBIT_MJD_MARKER + 9
c
     &  ,CUTBIT_OFFSET       = CUTBIT_IDTEL
     &  ,NCUTBIT             = CUTBIT_OFFSET
     &      )

c Define cut-masks with 64-BIT integers

      INTEGER*8
     &   CUTMASK8_SN             ! cutmask for each SN
     &  ,CUTMASK8_MJD(MXEPOCH)   ! cutmask vs. MJD, isn
     &  ,CUTMASK8_SN_ALL         ! all bits for SN cuts
     &  ,CUTMASK8_MJD_ALL        ! all bits for MJD cuts

      BYTE MJDMASK(MXEPOCH)  ! logical for each epoch,SN

      LOGICAL
     &   LSNCUTS       ! T=> passes cuts 1 to CUTBIT_MJD_MARKER
     &  ,PASS_USRCUTS  ! T=> user cuts with private code (Oct 2014)
     &  ,PASS_PRIVCUTS ! T=> pass cuts on private var in data header
     &  ,PASS_SIMCUTS  ! T=> pass cuts on SIMVAR

      INTEGER
     &   NCCID_LIST      ! size of SNCCID_LIST
     &  ,NCID_LIST       ! size of SNCID_LIST
     &  ,NCCID_IGNORE    ! size of SNCCID_IGNORE
     &  ,NCID_IGNORE     ! size of NCID_IGNORE
     &  ,CUTFLAG_SNANA   ! bits 0,1 -> LSNCUTS,LSNFITOK=T (for ntuple)
     &  ,ERRFLAG_FIT     ! error flag from fit (0=OK)

c stuff for hbook
      CHARACTER CUTVAR_NAME(NCUTBIT)*28

c define namelist cuts

      INTEGER
     &   SNTYPE_LIST(MXLISTNML)    ! (I) user list of types to select
     &  ,SNTYPE_IGNORE(MXLISTNML)  ! (I) types to ignore
     &  ,SNCID_LIST(MXLISTNML)     ! (I) user-list of integer CIDs
     &  ,CCDNUM_LIST(MXLISTNML)    ! (I) user-list of CCDNUMs (Dec 2017)
     &  ,SIM_TEMPLATE_INDEX_LIST(MXLISTNML)
     &  ,SNCID_IGNORE(MXLISTNML)   ! (I) list of CIDs to ignore
     &  ,PHOTFLAG_MSKREJ(5)        ! (I) PHOTFLAG mask-bits to reject
     &                             !     1,2 => logical OR,AND
     &  ,IDTEL_LIST(MXLISTNML)     ! computed from SNTEL_LIST
     &  ,IDFIELD_LIST(MXLISTNML)   !
     &  ,NIDTEL_LIST               ! size of list
     &  ,NIDFIELD_LIST
     &  ,NSNTYPE_LIST           ! size of SNTYPE_LIST
     &  ,NCCDNUM_LIST           ! size of CCDNUM list
     &  ,NSIM_TEMPLATE_INDEX_LIST
     &  ,NSNTYPE_IGNORE
     &  ,NFILT_SNRMAX           ! number of CUTWIN_SNRMAX_FILT cuts
     &  ,IFILT_SNRMAX(MXFILT_SNRMAX) ! store filt index for SNRMAX cuts
     &  ,NFILT_HOST_SBFLUX      ! Nfilt to cut on HOST_SBFLUX
     &  ,IFILT_HOST_SBFLUX(MXFILT_SNRMAX)

      LOGICAL DOALL_SNTEL, DOALL_SNFIELD

      character
     &    SNTEL_LIST(MXLISTNML)*20    ! (I) list of telescopes (or ALL)
     &   ,SNFIELD_LIST(MXLISTNML)*60  ! (I) list of fields to use (or 'ALL')
     &   ,SNCCID_LIST(MXLISTNML)*(MXCHAR_CCID)   ! (I) list of CIDs to process
     &   ,SNCCID_LIST_ALL(MXCCID_LIST)*(MXCHAR_CCID) ! combined list to process
     &   ,SNCCID_IGNORE(MXLISTNML)*(MXCHAR_CCID) ! (I) list of CIDs to ignore
     &   ,SNCID_IGNORE_FILE*(MXCHAR_FILENAME) ! (I) file with list of CIDs to ignore

      REAL
     &   snlc_cutvar(NCUTBIT_SNLC)
     &  ,snep_cutvar(CUTBIT_MJD_MARKER:NCUTBIT,MXEPOCH)
     &  ,cutwin_var(2,NCUTBIT)
     &  ,cutwin_cid(2)           ! candidate id
     &  ,cutwin_redshift(2)
     &  ,cutwin_redshift_err(2)
     &  ,cutwin_ra(2)
     &  ,cutwin_dec(2)
     &  ,cutwin_hostsep(2)      ! cut on host-SN sep, arcsec
     &  ,cutwin_sbflux_filt(2,MXFILT_SNRMAX)
     &  ,cutwin_Nepoch(2)       !
     &  ,cutwin_snrmax_filt(2,MXFILT_SNRMAX)   ! filled from SNCUT_SNRMAX
     &  ,cutwin_snrmax(2)       ! global SNRMAX cut
     &  ,cutwin_snrmax2(2)       ! 2nd global SNRMAX cut
     &  ,cutwin_nfield(2)        ! number of fields (usually 1)
     &  ,cutwin_mwebv(2)         ! Galactic extinction
     &  ,cutwin_nseason_active(2)
     &  ,cutwin_cutflag_reqep(2)       ! for internal use only
     &  ,cutwin_cutflag_private(2)     ! for internal use only
     &  ,cutwin_cutflag_simvar(2)      ! for internal use only
     &  ,cutwin_cutflag_usrcuts(2)     ! for internal use only
     &  ,cutwin_nfilt_snrmax(2)  ! Nfilt passing global SNRMAX cut
     &  ,cutwin_nfilt_snrmax2(2) ! Nfilt passing 2nd-best SNRMAX
     &  ,cutwin_nfilt_trestmin(2)! Nfilt passing Trestmin cut
     &  ,cutwin_nfilt_trestmax(2)! Nfilt passing Trestmax cut
     &  ,cutwin_nfilt_trest2(2)  ! Nfilt passinng Trest2 cut
     &  ,cutwin_nband_thresh(2)  ! number of ugriz bands above $band_THRESH
     &  ,cutwin_Trestmin(2)      ! window for earliest epoch, rest frame days
     &  ,cutwin_Trestmax(2)      ! window for latest epoch, rest frame day
     &  ,cutwin_TrestRange(2)    ! window on Trestmax - Trestmin
     &  ,cutwin_Trest2(2)        ! window for CUTWIN_NFILT_TREST2
     &  ,cutwin_Tgapmax(2)       ! max Trest-gap within cutwin_TREST(2)
     &  ,cutwin_T0gapmax(2)      ! max Trest-gap near peak
     &  ,cutwin_Tobsmin(2)
     &  ,cutwin_Tobsmax(2)
     &  ,cutwin_peakmjd(2)       ! cut on search peakMJD (for rates)
c
c  Below are Epoch-dependent cuts

     &  ,cutwin_psf(2)           ! PSF cut, FWHM, ARCSEC
     &  ,cutwin_zp(2)            ! ZEROPT cut, ADU or NPE
     &  ,cutwin_zpADU(2)         ! ZEROPT cut, ADU
     &  ,cutwin_zpNPE(2)         ! ZEROPT cut, NPE
     &  ,cutwin_moondist(2)      ! dist betwee SN and moon, deg
     &  ,cutwin_photprob(2)
     &  ,cutwin_Trest(2)         ! window for all epochs, rest frame days
     &  ,cutwin_MJD(2)           ! MJD window
     &  ,CUTWIN_MJD_EXCLUDE(2)          ! MJD window to exclude
     &  ,CUTWIN_MJD_INCLUDE(2)          ! require  MJD in this window
     &  ,CUTWIN_NMJD_INCLUDE(2)         ! for internal use only
     &  ,cutwin_SEARCHEFF_MASK(2) ! for SIM only
     &  ,cutwin_snrmin_filt(2,MXFILT_OBS) ! filled from EPCUT_SNRMIN
     &  ,cutwin_restlam(2)                ! cut on <LAMREST>, no cutBit
     &  ,cutwin_lamrest(2)                ! same as above
     &  ,cutwin_lamobs(2)                 ! cut on <LAMOBS>, no cutBit
c
c define character strings to specify cuts and filters
c 'SNCUT' specifies cut on each SN
c 'EPCUT' specifies selection cut on each epoch

      CHARACTER*(MXCHAR_CUTNAME)
     &  SNCUT_SNRMAX  ! max SNR required in each passband
                      ! example: 'u 10.  g 5.0  r 5.0  i 5.0  z -10.'
c
     & ,EPCUT_SNRMIN  ! min SNR accepted for each epoch/filter
                      ! example: 'u 20000.  g -4.  r -4.  i -4.  z 20000.'
c
     & ,EPCUT_SKYSIG  ! max SKY noise accepted each epoch/filter
                      ! example: 'u 20.  g 50.  r 80.  i 120.  z 200.'
c
     & ,EPCUT_PSFSIG  ! max PSF (arcsec) accepted each epoch/filter
                      ! example: 'u 1.8  g 1.5  r 1.6  i 1.8  z 2.0'
c
     & ,SNCUT_HOST_SBFLUX ! max allows surface brightness,
                          ! example 'r 1000.' -> SBFLUX < 1000

c ----------
c systematic tests for calibration:
c 'U 01 B -0.01' => shift U & B mags of primary

     &  ,MAGOBS_SHIFT_PRIMARY  ! shift primary mags  (for syst test)
     &  ,MAGOBS_SHIFT_ZP       ! shift zero points (e.g.,AB off for SDSS)
     &  ,MAGREST_SHIFT_PRIMARY ! idem for rest-frame mags
     &  ,MAGREST_SHIFT_ZP      ! idem for rest-frame filters
     &  ,FILTER_LAMSHIFT       ! e.g., 'r -2.4  i 6.2'  ! in Angstroms
c -------
c Fluxcal fudges for systematic tests (fudge photometry offsets and errors)
c Note that error is added in quadrature; ERROR<0 is subtracted in quadrature.
     &  ,FUDGE_FLUXCAL_OFFSET
     &  ,FUDGE_FLUXCAL_ERROR   ! fudge net FLUXCAL error in each band
     &  ,FUDGE_FLUXCAL_ERRPIX  ! per-pixel error --> FLUXCAL error propto PSF
     &  ,FUDGE_MAG_ERROR       ! Oct 2013: add stat error per band
     &  ,FUDGE_MAG_COVERR      ! Feb 2015: add covariant error per band

      REAL
     &   MAGOBS_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGOBS_SHIFT_ZP_USER(MXFILT_ALL)      ! user ZP shift
     &  ,MAGOBS_SHIFT_ZP_FILT(MXFILT_ALL)      ! user ZP shift + system ZP
     &  ,FLUXSCALE_ZP_FILT(MXFILT_ALL)         ! corresponding flux scale
     &  ,MAGREST_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGREST_SHIFT_ZP_USER(MXFILT_ALL)      ! shift zero points
     &  ,MAGREST_SHIFT_ZP_FILT(MXFILT_ALL)      ! shift zero points
     &  ,FILTER_LAMSHIFT_FILT(MXFILT_ALL)       ! shift filter trans
c
     &  ,MAGOBS_SHIFT_PRIMARY_PARAMS(3)       ! poly-fun of lambda;
     &  ,MAGOBS_SHIFT_ZP_PARAMS(3)            ! A0 + A1*LAM + A2*LAM^2
     &  ,FUDGE_FLUXCAL_OFFSET_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERRPIX_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_COVERR_FILT(MXFILT_ALL)
     &  ,MWEBV_SCALE        ! scale MW extinc for syst test
     &  ,MWEBV_SHIFT        ! shift MW extinc

      REAL
     &   RV_MWCOLORLAW   ! (I) RV for Galactic extinction

      INTEGER
     &   OPT_MWCOLORLAW    ! (I) MW color law opt (89, 94, 99 ...)
     &  ,OPT_MWEBV         ! (I) option to modify SFD maps

      LOGICAL
     &   USE_MWCOR          ! (I) T=> correct data flux for MW extinc;
                            !     F=> leave data, correct fit model for MW.
     &  ,DOBUG_LAMRANGE     ! (I) implement old lambda-range bug

      REAL
     &   REDSHIFT_FINAL_SHIFT    ! artificial shift
     &  ,FLUXERRCALC_ZPTERR

      COMMON / SNCUTCOM / CUTMASK8_SN, CUTMASK8_MJD
     &   ,CUTMASK8_SN_ALL, CUTMASK8_MJD_ALL
     &   ,cutwin_var, cutvar_name
     &   ,snlc_cutvar, snep_cutvar
     &   ,cutwin_snrmin_filt, LSNCUTS
     &   ,PASS_USRCUTS, PASS_PRIVCUTS, PASS_SIMCUTS
     &   ,CUTFLAG_SNANA , ERRFLAG_FIT
     &   ,NCCID_LIST, NCCID_IGNORE, NCID_LIST, NCID_IGNORE
     &   ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST, SNCCID_LIST_ALL
     &   ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE
     &   ,SNTYPE_IGNORE, CCDNUM_LIST, SIM_TEMPLATE_INDEX_LIST
     &   ,SNTEL_LIST,   IDTEL_LIST,   NIDTEL_LIST,   DOALL_SNTEL
     &   ,SNFIELD_LIST, IDFIELD_LIST, NIDFIELD_LIST, DOALL_SNFIELD
     &   ,NSNTYPE_LIST, NCCDNUM_LIST, NSIM_TEMPLATE_INDEX_LIST
     &   ,NSNTYPE_IGNORE, NFILT_SNRMAX, IFILT_SNRMAX
     &   ,NFILT_HOST_SBFLUX, IFILT_HOST_SBFLUX
     &   ,PHOTFLAG_MSKREJ
     &   ,USE_MWCOR, DOBUG_LAMRANGE
     &   ,MAGOBS_SHIFT_PRIMARY, MAGOBS_SHIFT_PRIMARY_FILT
     &   ,MAGOBS_SHIFT_ZP, MAGOBS_SHIFT_ZP_FILT, FLUXSCALE_ZP_FILT
     &   ,MAGOBS_SHIFT_ZP_USER
     &   ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS
     &   ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_PRIMARY_FILT
     &   ,MAGREST_SHIFT_ZP, MAGREST_SHIFT_ZP_FILT
     &   ,MAGREST_SHIFT_ZP_USER
     &   ,FILTER_LAMSHIFT, FILTER_LAMSHIFT_FILT
     &   ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX
     &   ,FUDGE_FLUXCAL_OFFSET_FILT
     &   ,FUDGE_FLUXCAL_ERROR_FILT, FUDGE_FLUXCAL_ERRPIX_FILT
     &   ,FUDGE_MAG_ERROR, FUDGE_MAG_ERROR_FILT
     &   ,FUDGE_MAG_COVERR, FUDGE_MAG_COVERR_FILT
     &   ,RV_MWCOLORLAW, OPT_MWCOLORLAW, OPT_MWEBV
     &   ,MWEBV_SCALE, MWEBV_SHIFT
     &   ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR
     &   ,CUTWIN_MJD_EXCLUDE, CUTWIN_MJD_INCLUDE

      COMMON / BYTEMASKCOM / MJDMASK
ccc mark delete      COMMON / CIDMASKCOM / MJDMASK,CIDMASK

      EQUIVALENCE
     &    ( cutwin_var(1,cutbit_cid),       cutwin_cid )
     &   ,( cutwin_var(1,cutbit_redshift),  cutwin_redshift )
     &   ,( cutwin_var(1,cutbit_redshift_err),cutwin_redshift_err )
     &   ,( cutwin_var(1,cutbit_ra),        cutwin_ra )
     &   ,( cutwin_var(1,cutbit_dec),       cutwin_dec )
     &   ,( cutwin_var(1,cutbit_hostsep),   cutwin_hostsep )
     &   ,( cutwin_var(1,cutbit_Nepoch),    cutwin_Nepoch )
     &   ,( cutwin_var(1,cutbit_psf),       cutwin_psf )
     &   ,( cutwin_var(1,cutbit_zp),        cutwin_zp  )
     &   ,( cutwin_var(1,cutbit_moondist),  cutwin_moondist  )
     &   ,( cutwin_var(1,cutbit_photprob),  cutwin_photprob  )
     &   ,( cutwin_var(1,cutbit_Nband_thresh), cutwin_nband_thresh )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax), cutwin_nfilt_snrmax )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax2),cutwin_nfilt_snrmax2 )
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmin),cutwin_nfilt_trestmin)
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmax),cutwin_nfilt_trestmax)
     &   ,( cutwin_var(1,cutbit_Nfilt_trest2),  cutwin_nfilt_trest2)
     &   ,( cutwin_var(1,cutbit_Trestmin),  cutwin_Trestmin )
     &   ,( cutwin_var(1,cutbit_Trestmax),  cutwin_Trestmax )
     &   ,( cutwin_var(1,cutbit_TrestRange),cutwin_TrestRange )
     &   ,( cutwin_var(1,cutbit_Trest2),    cutwin_Trest2 )
     &   ,( cutwin_var(1,cutbit_Tgapmax),   cutwin_Tgapmax )
     &   ,( cutwin_var(1,cutbit_T0gapmax),  cutwin_T0gapmax )
     &   ,( cutwin_var(1,cutbit_Tobsmin),   cutwin_Tobsmin )
     &   ,( cutwin_var(1,cutbit_Tobsmax),   cutwin_Tobsmax )
     &   ,( cutwin_var(1,cutbit_Trest),     cutwin_Trest )
     &   ,( cutwin_var(1,cutbit_MJD),       cutwin_MJD   )
     &   ,( cutwin_var(1,cutbit_peakmjd),   cutwin_peakmjd )
     &   ,( cutwin_var(1,cutbit_NMJD_INCLUDE),  cutwin_NMJD_INCLUDE )
     &   ,( cutwin_var(1,cutbit_search),        cutwin_searcheff_mask )
     &   ,( cutwin_var(1,cutbit_snrmax),        cutwin_snrmax  )
     &   ,( cutwin_var(1,cutbit_snrmax2),       cutwin_snrmax2  )
     &   ,( cutwin_var(1,cutbit_nfield),        cutwin_nfield  )
     &   ,( cutwin_var(1,cutbit_mwebv),         cutwin_mwebv  )
     &   ,( cutwin_var(1,cutbit_nseason_active),cutwin_nseason_active)
     &   ,( cutwin_var(1,cutbit_reqep),        cutwin_cutflag_reqep)
     &   ,( cutwin_var(1,cutbit_private),      cutwin_cutflag_private)
     &   ,( cutwin_var(1,cutbit_simvar),       cutwin_cutflag_simvar)
     &   ,( cutwin_var(1,cutbit_usrcuts),      cutwin_cutflag_usrcuts)
     &   ,( cutwin_var(1,cutbit_offset_sbflux+1), cutwin_sbflux_filt)
     &   ,( cutwin_var(1,cutbit_offset_snrmax+1), cutwin_snrmax_filt)
c
     &   ,( cutwin_restlam(1), cutwin_lamrest(1) )



C +CDE,KCORCOM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 

C +CDE,KTABLECM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 
      REAL*4  R4KCORTABLE1D(MXTABLE1D_KCOR)
      COMMON / R4KCORCOM / R4KCORTABLE1D


      INTEGER
     &    LUNKCOR, MXKCOR
     &   ,NKDIM, N4DIM
     &   ,IDMAP_KCOR, IDMAP_AVWARP, IDMAP_LCMAG, IDMAP_MWXT
     &   ,KDIM_T, KDIM_Z, KDIM_AV, KDIM_ifiltr, KDIM_ifilto

      PARAMETER (
     &   LUNKCOR      = 98  ! log unit number to read
     &  ,MXKCOR       = 100 ! max number of KCOR tables to read
     &  ,KDIM_T       = 1  ! epoch KCOR index
     &  ,KDIM_Z       = 2  ! redshift KCOR index
     &  ,KDIM_AV      = 3  ! AV kcor index
     &  ,KDIM_ifiltr  = 4  ! rest-filter KCOR index
     &  ,KDIM_ifilto  = 5  ! obs-filter KCOR index
     &  ,NKDIM        = 5  ! number of KCOR dimensions in map
     &  ,IDMAP_KCOR   = 10 ! for INIT_1DINDEX (any integer is OK)
     &  ,IDMAP_AVWARP = 11
     &  ,IDMAP_LCMAG  = 12
     &  ,IDMAP_MWXT   = 13
     &  ,N4DIM        = 4  ! for 4-d tables
     &   )

c K-cor & photometry-template info

      INTEGER
     &   NCALL_RDKCOR ! # times RDKCOR is called
     &  ,NKCOR_STORE  ! number of stored K-corr tables
     &  ,NZBIN_KCOR   ! # Z-bins in each K-cor table
     &  ,NTBIN_KCOR   ! # T-bins "    "
     &  ,NAVBIN_KCOR  ! # AV bins
     &  ,IAV0         ! IAV index with AV=0
     &  ,IZ0          ! IZ with Z=0
     &  ,NLbin_SNSED      ! # lambda bins for SN SED
     &  ,NTbin_SNSED      ! # epoch bins for SN SED
     &  ,NERR_AVWARP      ! # errors trapped in GET_AVWARP
c
     &  ,NBINTOT_KCOR     ! Total KCOR bins
     &  ,NBINTOT_AVWARP
     &  ,NBINTOT_LCMAG
     &  ,NBINTOT_MWXT
     &  ,NBIN_KCOR(NKDIM)
     &  ,NBIN_AVWARP(NKDIM)
     &  ,NBIN_LCMAG(NKDIM)
     &  ,NBIN_MWXT(NKDIM)

      REAL
     &   Zrange_Kcor(2)    ! min,max redshift in KCOR table
     &  ,Zrange_Kcor_LU(2) ! min,max redshift for lookup
     &  ,Trange_Kcor(2)  ! min,max rest-epoch (days) in KCOR table
     &  ,AVrange_Kcor(2) ! min,max AV
     &  ,ZBINSIZE_KCOR   ! redshift binsize
     &  ,TBINSIZE_KCOR   ! epoch bin size (days)
     &  ,AVbinsize_KCOR  !
     &  ,Cbinsize_AVWARP !
     &  ,GRIDVAL_AV(MXAVBIN_KCOR)  ! store AV(iav) for convenience
     &  ,GRIDVAL_Z(MXZBIN_KCOR)    ! store Z(iz)   for convenience
     &  ,GRIDVAL_T(MXTBIN_KCOR)    ! store Trest(it) for convenience
     &  ,GRIDVAL_C(MXCBIN_AVWARP)  ! store color for AVWARP table
     &  ,Lrange_SNSED(2), LBINSIZE_SNSED
     &  ,Trange_SNSED(2), TBINSIZE_SNSED
     &  ,FLUX_SNSED(MXLAMBIN_SNSED*MXTBIN_KCOR)

      REAL*4
     &   AVWARP_TABLE1D(MXTABLE1D_AVWARP)
     &  ,LCMAG_TABLE1D(MXTABLE1D_LCMAG)
     &  ,MWXT_TABLE1D(MXTABLE1D_MWXT)  ! MilkyWay extinction table

      INTEGER
     &   HIDMAG(MXAVBIN_KCOR,MXFILT_ALL)
     &  ,HIDKCOR(MXFILT_ALL,MXFILT_ALL)

      LOGICAL
     &   USE_AVWARPTABLE
     &  ,LRDZPOFF
     &  ,RDKCOR_STANDALONE
     &  ,EXIST_KCOR(MXFILT_ALL,MXFILT_ALL)

      CHARACTER RESTKCOR_FILTERS*20    ! rest-frame filters from KCOR file

      REAL RVMW_RDKCOR ! RV used for MW table extension
      INTEGER OPT_MWCOLORLAW_RDKCOR

      COMMON / KCORCOM / NCALL_RDKCOR, NKCOR_STORE, RDKCOR_STANDALONE
     &      ,NBIN_KCOR,   NBINTOT_KCOR
     &      ,NBIN_AVWARP, NBINTOT_AVWARP
     &      ,NBIN_LCMAG,  NBINTOT_LCMAG
     &      ,NBIN_MWXT,   NBINTOT_MWXT
     &      ,NZBIN_KCOR,    NTBIN_KCOR,    NAVBIN_KCOR
     &      ,Zrange_KCOR,   Trange_KCOR,   AVrange_KCOR
     &      ,ZBINSIZE_KCOR, TBINSIZE_KCOR, AVbinsize_KCOR
     &      ,CBINSIZE_AVWARP, IAV0, IZ0
     &      ,GRIDVAL_Z, GRIDVAL_T, GRIDVAL_AV, GRIDVAL_C
     &      ,HIDKCOR, HIDMAG
     &      ,Zrange_KCOR_LU, RESTKCOR_FILTERS
     &      ,USE_AVWARPTABLE, LRDZPOFF, NERR_AVWARP
     &      ,NLbin_SNSED, Lrange_SNSED, LBINSIZE_SNSED
     &      ,NTbin_SNSED, Trange_SNSED, TBINSIZE_SNSED
     &      ,FLUX_SNSED
     &      ,LCMAG_TABLE1D, MWXT_TABLE1D, AVWARP_TABLE1D
     &      ,EXIST_KCOR
     &      ,RVMW_RDKCOR, OPT_MWCOLORLAW_RDKCOR

c  July 2016 - define optional spectrograph info
      CHARACTER
     &   SPECTROGRAPH_INSTRUMENT*80
     &  ,SPECTROGRAPH_FILTERLIST*(MXFILT_ALL)
      INTEGER
     &   NFILTDEF_SPECTROGRAPH
     &  ,IFILTDEF_SPECTROGRAPH(MXFILT_ALL)

      COMMON / SPECTROGRAPH_COM /
     &   SPECTROGRAPH_INSTRUMENT, SPECTROGRAPH_FILTERLIST
     &  ,NFILTDEF_SPECTROGRAPH,   IFILTDEF_SPECTROGRAPH

c define temporary 'XXX_RDKCOR' arrays that are used only
c for reading and parsing.

      CHARACTER
     &   PRIMARY_NAME_RDKCOR(10)*40          ! e.g, AB, BD17 ..
     &  ,FILTER_NAME_RDKCOR(MXFILT_ALL)*40   ! e.g, SDSS-g, Bessell90-V
     &  ,KCORINFO_STRING_RDKCOR(MXKCOR)*60   ! Kcor K_xy obsfilt restfilt

      INTEGER
     &   IVER_RDKCOR
     &  ,NPRIM_RDKCOR, NKCOR_RDKCOR
     &  ,INDX_PRIMARY_RDKCOR(MXFILT_ALL)
     &  ,NFILTDEF_RDKCOR               ! Number of all filters
     &  ,IFILTDEF_RDKCOR(MXFILT_ALL)   ! vs. ifilt_rdkcor
     &  ,NFILTOBS_RDKCOR               ! number of obs filters
     &  ,IFILTOBS_RDKCOR(MXKCOR)       ! vs. ikcor
     &  ,IFILT2_RDKCOR(2,MXKCOR)       ! 1=rest:2=obs, 1:NKCOR_STORE
     &  ,MASKFILT_RDKCOR(MXFILT_ALL)   ! bit0,1 -> rest,obs
     &  ,NFILT_DUPLICATE_RDKCOR        ! keep count of duplicate filters
     &  ,IKCOR_RDKCOR(MXKCOR)          ! ikcor vs. istore
     &  ,IPRIM_REF_RDKCOR              ! index for primary  ref.

      REAL
     &   MAG_PRIMARY_RDKCOR(MXFILT_ALL)    ! native mags
     &  ,ZPOFF_PRIMARY_RDKCOR(MXFILT_ALL)  ! mag(native) - mag(synth)
     &  ,ZPOFF_SNPHOT_RDKCOR(MXFILT_ALL)   ! apply to data (from ZPOFF.DAT)

      LOGICAL
     &   ISFITS_RDKCOR
     &  ,LVERBOSE_RDKCOR
     &  ,ISLAMSHIFT_RDKCOR(0:MXFILT_ALL)

      REAL      NULLF_RDKCOR
      INTEGER   NULLI_RDKCOR
      CHARACTER NULLS_RDKCOR*20
      LOGICAL   ANYF_RDKCOR

      COMMON / RDKCORCOM /
     &   IVER_RDKCOR, ISFITS_RDKCOR,LVERBOSE_RDKCOR,ISLAMSHIFT_RDKCOR
     &  ,NKCOR_RDKCOR
     &  ,NPRIM_RDKCOR, INDX_PRIMARY_RDKCOR, PRIMARY_NAME_RDKCOR
     &  ,FILTER_NAME_RDKCOR, KCORINFO_STRING_RDKCOR
     &  ,IFILTDEF_RDKCOR,  NFILTDEF_RDKCOR
     &  ,IFILTOBS_RDKCOR,  NFILTOBS_RDKCOR
     &  ,IFILT2_RDKCOR
     &  ,NFILT_DUPLICATE_RDKCOR, IKCOR_RDKCOR, MASKFILT_RDKCOR
     &  ,IPRIM_REF_RDKCOR
c
     &  ,MAG_PRIMARY_RDKCOR
     &  ,ZPOFF_PRIMARY_RDKCOR
     &  ,ZPOFF_SNPHOT_RDKCOR
c
     &  ,NULLF_RDKCOR,NULLI_RDKCOR,NULLS_RDKCOR,ANYF_RDKCOR


C +CDE,SNHOSTCOM. inserted below by ypatchy. 

c Host galaxy parameters.
c Dec 17, 2012 - add SNHOST_MAGOBS


c logical flags for CWNT
      LOGICAL
     &   EXIST_SNHOST_ANGSEP
     &  ,EXIST_SNHOST_DDLR
     &  ,EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT
     &  ,EXIST_SNHOST_LOGMASS
     &  ,EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS
     &  ,EXIST_SNHOST_SB

      INTEGER*8  SNHOST_OBJID(MXSNHOST)    ! int id
      REAL*8    DSNHOST_OBJID(MXSNHOST)   ! for tables only

      REAL
     &   SNHOST_ANGSEP(MXSNHOST)      ! SN-host sep, arcsec
     &  ,SNHOST_DDLR(MXSNHOST)        ! SNSEP/DLR
     &  ,SNHOST_CONFUSION             ! HC analog from Gupta 2016
     &  ,SNHOST_ZPHOT(MXSNHOST), SNHOST_ZPHOT_ERR(MXSNHOST)
     &  ,SNHOST_ZSPEC(MXSNHOST), SNHOST_ZSPEC_ERR(MXSNHOST)
     &  ,SNHOST_LOGMASS(MXSNHOST)
     &  ,SNHOST_LOGMASS_ERR(MXSNHOST)
     &  ,SNHOST_sSFR(MXSNHOST)
     &  ,SNHOST_sSFR_ERR(MXSNHOST)
     &  ,SNHOST_SBFLUXCAL(MXFILT_ALL)  ! surface brightness FLUXCAL /asec^2
     &  ,SNHOST_MAGOBS(MXFILT_ALL,(MXSNHOST))       ! observer-frame mags
     &  ,SNHOST_MAGOBS_ERR(MXFILT_ALL,(MXSNHOST))   ! error on above

      REAL*8
     &    SNHOST8_RA(MXSNHOST)
     &   ,SNHOST8_DEC(MXSNHOST)

c quantites which do not depend on which host
      INTEGER*4  SNHOST_NMATCH, SNHOST_NMATCH2  ! number of host matches
      REAL
     &   SNHOST_SBFLUXCAL_ERR(MXFILT_ALL)
     &  ,SNHOST_SBMAG(MXFILT_ALL)        ! surface brightness mag/asec^2

      COMMON / SNHOSTCOM /
     &   EXIST_SNHOST_ANGSEP, EXIST_SNHOST_DDLR, EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT,  EXIST_SNHOST_LOGMASS, EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS, EXIST_SNHOST_SB
c
     &  ,SNHOST_NMATCH, SNHOST_NMATCH2
     &  ,SNHOST_ANGSEP, SNHOST_DDLR, SNHOST_CONFUSION
     &  ,SNHOST_ZPHOT,   SNHOST_ZPHOT_ERR
     &  ,SNHOST_ZSPEC,   SNHOST_ZSPEC_ERR
     &  ,SNHOST_LOGMASS, SNHOST_LOGMASS_ERR
     &  ,SNHOST_sSFR,    SNHOST_sSFR_ERR
     &  ,SNHOST_SBFLUXCAL, SNHOST_SBFLUXCAL_ERR, SNHOST_SBMAG
     &  ,SNHOST_MAGOBS, SNHOST_MAGOBS_ERR

      COMMON / SNHOSTCOM8 /
     &   SNHOST_OBJID, DSNHOST_OBJID, SNHOST8_RA, SNHOST8_DEC

cc +CDE,PARSECOM.

C +CDE,SNSIMCOM. inserted below by ypatchy. 

c simulation parameters (if FAKE=2)

      REAL*8  SIM8_RA, SIM8_DECL

      REAL
     &   SIM_REDSHIFT_HELIO
     &  ,SIM_REDSHIFT_CMB
     &  ,SIM_REDSHIFT_HOST
     &  ,SIM_VPEC
     &  ,SIM_DLMAG
     &  ,SIM_LENSDMU
     &  ,SIM_MWEBV
     &  ,SIM_MWRV
     &  ,SIM_SALT2x0
     &  ,SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW
     &  ,SIM_PEAKMJD
     &  ,SIM_EXPOSURE_TIME(MXFILT_OBS)  ! relative exposure time
     &  ,SIM_PEAKMAG(MXFILT_OBS)
     &  ,SIM_EPMAGOBS(MXEPOCH)      ! true epoch mag at each filter/epoch
     &  ,SIM_EPFLUXCAL(MXEPOCH)     ! true fluxcal at each filter/epoch
     &  ,SIM_EPSNRMON(MXEPOCH)      ! optional SNR at MAGMONITOR
     &  ,SIM_EPMAGREST(MXEPOCH)     ! true rest-frame mag
     &  ,SIM_EPFLUXCAL_HOSTERR(MXEPOCH) ! true error from host noise
     &  ,SIMSED_PARVAL(MXPAR_SIMSED)
     &  ,BYOSED_PARVAL(MXPAR_SIMSED)    ! Dec 10 2018
     &  ,LCLIB_PARVAL(MXPAR_LCLIB)
     &  ,SIM_HOSTLIB_PARVAL(MXPAR_SIMSED) ! HOSTLIB params
     &  ,SIM_MAGSMEAR_COH
     &  ,SIM_TEMPLATEMAG(MXFILT_ALL)  ! image-sub template, not LCLIB template
     &  ,SIM_LCWIDTH(MXFILT_ALL)      ! computed from SIM_EPFLUXCAL
c
     &  ,SIM_MODELGRID_TOBS(MXEP_MODELGRID)    ! for SNANA+SIM_MAGOBS table
     &  ,SIM_MODELGRID_MAGOBS(MXEP_MODELGRID,MXFILT_OBS/2) ! idem

      INTEGER
     &   SIM_MODEL_INDEX      ! model index (MLCS,SALT2,NON1a ...)
     &  ,SIM_GENTYPE          ! generated SNTYPE; SIM_TYPE_INDEX
     &  ,SIM_TEMPLATE_INDEX   ! template index for NON1ASED, SIMSED, LCLIB ...
     &  ,SIM_SEARCHEFF_MASK   ! bits 1,2 => found by software,humans
     &  ,SIM_LIBID            ! LIBID for each event
     &  ,SIM_NGEN_LIBID       ! NGEN for this LIBID (usually 1)
     &  ,SIM_NOBS_UNDEFINED   ! NOBS where model is undefined
     &  ,SIM_NSUBSAMPLE_MARK  ! Number of marked sub-samples
     &  ,SIM_SUBSAMPLE_INDEX  ! sub-sample index
     &  ,SIM_REDSHIFT_FLAG    ! points to source of redshift
     &  ,SIMOPT_MWCOLORLAW    ! option for MW color law
     &  ,SIMOPT_MWEBV         ! option to modify MWEBV_SFD map
     &  ,NPAR_SIMSED
     &  ,NPAR_BYOSED
     &  ,NPAR_LCLIB
     &  ,NPAR_SIM_HOSTLIB
     &  ,SIM_EPFILTREST(MXEPOCH)
     &  ,SIMLIB_MSKOPT        ! SIMLIB option mask (Dec 2015)
     &  ,NEP_SIM_MODELGRID

      INTEGER*8  SIM_HOSTLIB_GALID

      CHARACTER
     &   SIMNAME_MODEL*40     ! SALT2, mlcs2k2, NON1A, etc ...
     &  ,SIMNAME_TYPE*12      ! Ia, Ib, II, IIN, etc ...
     &  ,SIMNAME_SHAPEPAR*40  ! DELTA, x1, stretch ...
     &  ,SIMNAME_SHAPELAW*40  ! alpha ...
     &  ,SIMNAME_COLORPAR*40  ! AV, c ...
     &  ,SIMNAME_COLORLAW*40  ! RV,  beta ...
     &  ,SIMLIB_FILENAME*200  ! SIMLIB file used to generate sim
     &  ,SIMSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,SIMSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,BYOSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,BYOSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,LCLIB_KEYWORD(MXPAR_LCLIB)*80 ! full keyname in header
     &  ,LCLIB_PARNAME(MXPAR_LCLIB)*20 ! parname for fitres and ntuple
     &  ,SIM_HOSTLIB_KEYWORD(MXPAR_SIMSED)*60 ! full keyname
     &  ,SIM_HOSTLIB_PARNAME(MXPAR_SIMSED)*60 ! parname
     &  ,SIMNAME_SNRMON*40

      COMMON / SNSIMCOM / SIMNAME_MODEL, SIM_MODEL_INDEX
     &  ,SIM_REDSHIFT_HELIO, SIM_REDSHIFT_CMB, SIM_REDSHIFT_HOST
     &  ,SIM_REDSHIFT_FLAG
     &  ,SIM_VPEC, SIM_DLMAG, SIM_LENSDMU, SIM_SALT2x0, SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW, SIM_PEAKMJD, SIM_MWEBV,SIM_MWRV
     &  ,SIM_EXPOSURE_TIME, SIM_TEMPLATE_INDEX, SIM_TEMPLATEMAG
     &  ,SIM_GENTYPE, SIMNAME_TYPE, SIM_LCWIDTH
     &  ,SIM_SEARCHEFF_MASK, SIM_LIBID, SIM_NGEN_LIBID
     &  ,SIM_NOBS_UNDEFINED, SIM_NSUBSAMPLE_MARK, SIM_SUBSAMPLE_INDEX
     &  ,SIM_PEAKMAG, SIM_EPMAGOBS, SIM_EPFLUXCAL
     &  ,SIM_EPSNRMON, SIM_EPFLUXCAL_HOSTERR
     &  ,SIMNAME_SHAPEPAR,  SIMNAME_SHAPELAW
     &  ,SIMNAME_COLORPAR,  SIMNAME_COLORLAW, SIMNAME_SNRMON
     &  ,SIMLIB_FILENAME, SIMLIB_MSKOPT
     &  ,SIMOPT_MWCOLORLAW, SIMOPT_MWEBV
     &  ,SIM_EPFILTREST, SIM_EPMAGREST,  SIM_MAGSMEAR_COH
     &  ,NEP_SIM_MODELGRID, SIM_MODELGRID_TOBS, SIM_MODELGRID_MAGOBS

      COMMON / SIMSEDCOM /
     &  NPAR_SIMSED, SIMSED_PARVAL, SIMSED_PARNAME, SIMSED_KEYWORD

      COMMON / BYOSEDCOM /
     &   NPAR_BYOSED, BYOSED_PARVAL, BYOSED_PARNAME, BYOSED_KEYWORD

      COMMON / LCLIBCOM /
     &   NPAR_LCLIB, LCLIB_PARVAL, LCLIB_KEYWORD, LCLIB_PARNAME

      COMMON / SIMHOSTCOM /
     &    NPAR_SIM_HOSTLIB, SIM_HOSTLIB_PARVAL
     &   ,SIM_HOSTLIB_PARNAME, SIM_HOSTLIB_KEYWORD

      COMMON / SNSIMCOM8 /
     &         SIM8_RA, SIM8_DECL, SIM_HOSTLIB_GALID



C +CDE,SNLCINP. inserted below by ypatchy. 

c define general input to program:
c SN flux version, hbook files, cut-windows, etc ...
c Dec 30 2012: add KCOR_FILE for fits-formatted tables (same name as in sim)

      CHARACTER*(MXCHAR_FILENAME)  ! user input files
     &   NMLFILE          ! name of input namelist file
     &  ,HFILE_OUT        ! output hbook filename
     &  ,ROOTFILE_OUT     ! output filename for root
     &  ,TEXTFILE_PREFIX  ! prefix for ascii tables
     &  ,SNTABLE_LIST     ! list of SNTABLEs to create
     &  ,SNTABLE_FILTER_REMAP ! remap filters in tables
     &  ,KCOR_FILE            ! fits-formatted KCOR tables
     &  ,hfile_kcor           ! OBSOLETE input -> abort
     &  ,EPOCH_IGNORE_FILE    ! file with epochs to ignor
     &  ,SNMJD_LIST_FILE      ! list of "CID MJD" to process
     &  ,SNMJD_OUT_FILE       ! MJD_LIST out file
     &  ,MNFIT_PKMJD_LOGFILE  ! separate log for PKMJD fits
     &  ,ZPOFF_FILE           ! AB (ZP) off for obs-frame mags
     &  ,USERTAGS_FILE        ! optional int user tag per SN
     &  ,VPEC_FILE            ! pec. velocity corrections
     &  ,HEADER_OVERRIDE_FILE ! header-var overrides
     &  ,SNCID_LIST_FILE      ! list of CID's to process
     &  ,SIMLIB_OUT           ! write simlib entry for each event
     &  ,SIMLIB_ZPERR_LIST    ! 'abc .01 def .02 ghi .005'
     &  ,NONLINEARITY_FILE    ! non-linearity map (Apr 2017)
c
     &  ,FUDGE_HOSTNOISE_FILE  ! legacy: inflate error vs.  hostSB (Mar 2015)
     &  ,FLUXERRMODEL_FILE     ! DATA ONLY: err-fudge maps
     &  ,SIM_FLUXERRMODEL_FILE ! idem, SIM only
     &  ,MAGCOR_FILE           ! DATA ONLY: mag-cor for each CID-MJD-band
     &  ,SIM_MAGCOR_FILE       ! idem, SIM only

      CHARACTER  ! versions
     &   VERSION_PHOTOMETRY(MXVERS)*(MXCHAR_VERSION)   ! SN versions to read
     &  ,VERSION_REPLACE*(MXCHAR_VERSION)  ! replace previous cmd-line VERSION
     &  ,VERSION_REFORMAT_FITS*(MXCHAR_VERSION)   ! output version for
                                                  ! for translating -> FITS

      CHARACTER ! paths
     &   PRIVATE_DATA_PATH*(MXCHAR_PATH)    ! private data subdir
     &  ,FILTER_UPDATE_PATH*(MXCHAR_PATH)   ! SN-dependent filter response

      CHARACTER  ! misc
     &   LINE_ARGS(MXLINE_ARGS)*(MXCHAR_FILENAME)    ! command line args
     &  ,REFORMAT_KEYS*(MXCHAR_FILENAME)  ! global reformat info
     &  ,NONSURVEY_FILTERS*(MXFILT_ALL)   ! non-survey filters to add
     &  ,SNRMAX_FILTERS*(MXFILT_ALL)      ! list of filters for SNRMAX cuts
     &  ,FILTER_REPLACE*(MXFILT_ALL)      ! e.g., 'UGRIZ -> ugriz'
     &  ,FILTLIST_LAMSHIFT*(MXFILT_ALL)   ! list of lam-shifted filters
     &  ,PRIVATE_CUTWIN_STRING(MXCUT_PRIVATE)*(MXCHAR_CUTNAME)
     &  ,PRIVATE_REDSHIFT_CMB*60          ! use PRIVATE variable for redshift
     &  ,SIMVAR_CUTWIN_STRING*(MXCHAR_CUTNAME) ! cuts on SIM_XXX
     &  ,EARLYLC_STRING*(MXCHAR_CUTNAME)  ! see manual
     &  ,REQUIRE_EPOCHS_STRING*100  ! e.g., 'riz 10 7 20' ! uses CUTWIN_SNRMAX
     &  ,DUMP_STRING*100             ! 'funName CID-list'

      INTEGER
     &   NFIT_ITERATION       ! (I) number of fit iterations
     &  ,MINUIT_PRINT_LEVEL   ! -1 -> none
     &  ,INTERP_OPT           ! (I)  interp option (see INTERP_XXX params)
     &  ,NLINE_ARGS
     &  ,OPT_SETPKMJD            ! bit-mask options  to determined PKMJD
     &  ,MJDRANGE_SETPKMJD(2)    ! MJD range to use for SETPKMJD
     &  ,OPT_REFORMAT_TEXT       ! T=> re-write SNDATA files in text format
     &  ,OPT_REFORMAT_TERSE      ! T=> legacy key for text format
     &  ,OPT_REFORMAT_FITS       ! T=> re-write  in FITSformat (data only)
     &  ,OPT_REFORMAT_SALT2      ! T=> re-write SNDATA files in SALT2 format
     &  ,OPTSIM_LCWIDTH          ! T=> option to compute LC width
     &  ,OPT_DEBUG               ! T=> for internal debug
     &  ,JOBSPLIT(2)             ! process (1)-range of (2)=TOTAL
     &  ,JOBSPLIT_EXTERNAL(2)    ! passed by split_and_fit for text format
     &  ,MXLC_FIT                ! stop after this many fits passing all cuts
     &  ,PHOTFLAG_DETECT         ! used to count NOBS_DETECT and TLIVE_DETECT
     &  ,APPLY_FLUXCOR_BUG       ! see calling sequence in SNRECON
     &  ,FLUXERRMODEL_OPTMASK

      LOGICAL
     &   LDMP_SNFAIL        ! (I) T => dump reason for each SN failure
     &  ,LSIM_SEARCH_SPEC   ! (I) T => require simulated SPEC-tag
     &  ,LSIM_SEARCH_ZHOST  ! (I) T => require simulated zHOST
     &  ,LTEST_INTERP  ! (I) T => calls TEST_KCOR
     &  ,LTEST_KCOR    ! (I) T => calls TEST_INTERP
     &  ,LTEST_MAG     ! (I) T => test GET_MAGLC  and GET_MWXT
     &  ,LTEST_U3BAND  ! (I) T => require exactly 3 bands that include U
     &  ,USE_LINE_ARGS(MXLINE_ARGS)
     &  ,USE_SNHOST_ZPHOT    ! (I) T=> replace SNLC_REDSHIFT -> SNHOST_ZPHOT
c
     &  ,ABORT_ON_NOEPOCHS   ! T=> abort if there are no epochs to fit
     &  ,ABORT_ON_BADAVWARP  ! T=> abort if AVwarp cannot be determined
     &  ,ABORT_ON_BADZ       ! T=> abort on bad z in GET_KCOR8
     &  ,ABORT_ON_BADKCOR    ! T=> affects only the init stage (RDKCOR)
     &  ,ABORT_ON_BADSURVEY  ! T=> abort if SURVEY changes
     &  ,ABORT_ON_BADFILTER  ! T=> abort if fit band is not in SURVEY filters
     &  ,ABORT_ON_MARGPDF0   ! T=> abort if marginalized pdf=0 everywhere
     &  ,ABORT_ON_TRESTCUT   ! T=> abort if any Trest cut is set (for photoz)
     &  ,ABORT_ON_DUPLCID    ! T=> abort on duplicate CID
     &  ,ABORT_ON_DUPLMJD    ! T=> abort on repeat MJD+band (Jun 2017)
     &  ,ABORT_ON_NOPKMJD    ! T=> abort if no PKMJDINI (see OPT_SETPKMJD)
     &  ,ABORT_ON_BADHEADER  ! T=> abort if missing header info(HEADMASK_CHECK)
     &  ,USE_MINOS           ! T=> use MINOS instead of MINIMIZE
     &  ,LDMP_AVWARP         ! dump GET_AVWARP8 (debug only)
     &  ,LDMP_KCORFUN        ! dump KCORFUN8 (debug only)
     &  ,USESIM_SNIA         ! default True -> process simulated SNIa
     &  ,USESIM_SNCC         ! default True -> process simulated SNCC
     &  ,USESIM_TRUEFLUX     ! SNLC_FLUXCAL -> SIM_FLUXCAL

c LC plot info

      INTEGER
     &  MXLC_PLOT      ! (I) max number of plots to make with SNLCPAK
     & ,NCCID_PLOT     ! (I) size of SNCCID_PLOT array
      CHARACTER
     &  SNCCID_PLOT(MXLISTNML)*(MXCHAR_CCID) ! (I) string-list of CCIDs to plot

      REAL
     &  DTOBS_MODEL_PLOT  ! (I) binning (days) for overlaid best-fit model

c Oct 2014: variables to check for multi-season transient activity

      INTEGER
     &   MULTISEASON_OPTMASK      ! option(s) for GET_MULTISEASON
      REAL
     &   MULTISEASON_TGAP              ! new season if time-gap > this
     &  ,MULTISEASON_NREJECT_OUTLIER   ! num outliers to reject per season
     &  ,MULTISEASON_CHI2RED_ACTIVE    ! min chi2red for active season

      REAL
     &   SIM_PRESCALE        ! pre-scale applied only to SIM
     &  ,VPEC_ERR_OVERRIDE   ! override VPEC_ERR in data header (km/sec)


c define reference cosmological parameters (double precision!)
c Aug 7, 2009: each parameter has dimension 2 for value & error

      REAL*8
     &   H0_REF(2)    ! reference H0  (70 => 70 MPc^-1)
     &  ,OLAM_REF(2)  ! omega_matter
     &  ,OMAT_REF(2)  ! omega_lamgda
     &  ,ORAD_REF(2)  ! omega_rad
     &  ,W0_REF(2)    ! w = p/rho
     &  ,DWDA_REF(2)  ! dw/da [a = 1/(1+z)]

      INTEGER   NSIMVAR_CUTWIN
      REAL      SIMVAR_CUTWIN(2,MXCUT_PRIVATE)
      CHARACTER SIMVAR_CUTWIN_LIST(MXCUT_PRIVATE)*40

      COMMON   / SNLCINP /
     &      NLINE_ARGS, LINE_ARGS, USE_LINE_ARGS, nmlfile
     &    , VERSION_PHOTOMETRY, VERSION_REPLACE, VERSION_REFORMAT_FITS
     &    , JOBSPLIT, JOBSPLIT_EXTERNAL, SIM_PRESCALE, MXLC_FIT
     &    , PRIVATE_DATA_PATH, FILTER_UPDATE_PATH
     &    , NONSURVEY_FILTERS, SNRMAX_FILTERS, VPEC_ERR_OVERRIDE
     &    , FILTER_REPLACE, FILTLIST_LAMSHIFT
     &    , OPTSIM_LCWIDTH, OPT_REFORMAT_TERSE, OPT_REFORMAT_TEXT
     &    , OPT_REFORMAT_SALT2, REFORMAT_KEYS, OPT_REFORMAT_FITS
     &    , SNMJD_LIST_FILE, SNMJD_OUT_FILE, MNFIT_PKMJD_LOGFILE
     &    , EPOCH_IGNORE_FILE, SNCID_LIST_FILE
     &    , SIMLIB_OUT, SIMLIB_ZPERR_LIST, NONLINEARITY_FILE
     &    , hfile_out, rootfile_out, textfile_prefix
     &    , SNTABLE_LIST, SNTABLE_FILTER_REMAP
     &    , hfile_kcor, KCOR_FILE, ZPOFF_FILE, USERTAGS_FILE
     &    , VPEC_FILE, HEADER_OVERRIDE_FILE
     &    , FLUXERRMODEL_FILE,SIM_FLUXERRMODEL_FILE,FLUXERRMODEL_OPTMASK
     &    , MAGCOR_FILE, SIM_MAGCOR_FILE,  FUDGE_HOSTNOISE_FILE
     &    , NFIT_ITERATION, MINUIT_PRINT_LEVEL, INTERP_OPT, USE_MINOS
     &    , OPT_SETPKMJD, MJDRANGE_SETPKMJD, OPT_DEBUG
     &    , LDMP_SNFAIL, LSIM_SEARCH_SPEC, LSIM_SEARCH_zHOST
     &    , LDMP_AVWARP, LDMP_KCORFUN
     &    , LTEST_KCOR, LTEST_INTERP, LTEST_MAG, LTEST_U3BAND
     &    , USESIM_SNIA, USESIM_SNCC, USESIM_TRUEFLUX, USE_SNHOST_ZPHOT
     &    , ABORT_ON_NOEPOCHS, ABORT_ON_BADAVWARP, ABORT_ON_NOPKMJD
     &    , ABORT_ON_BADZ, ABORT_ON_BADKCOR, ABORT_ON_BADSURVEY
     &    , ABORT_ON_MARGPDF0, ABORT_ON_TRESTCUT
     &    , ABORT_ON_DUPLCID, ABORT_ON_DUPLMJD
     &    , ABORT_ON_BADHEADER, ABORT_ON_BADFILTER
     &    , SNCUT_SNRMAX, SNCUT_HOST_SBFLUX
     &    , EPCUT_SNRMIN, EPCUT_SKYSIG, EPCUT_PSFSIG
     &    , CUTWIN_LAMREST, CUTWIN_LAMOBS
     &    , PRIVATE_CUTWIN_STRING, PRIVATE_REDSHIFT_CMB
     &    , SIMVAR_CUTWIN_STRING
     &    , NSIMVAR_CUTWIN, SIMVAR_CUTWIN, SIMVAR_CUTWIN_LIST
     &    , EARLYLC_STRING, REQUIRE_EPOCHS_STRING, DUMP_STRING
     &    , MXLC_PLOT, NCCID_PLOT, SNCCID_PLOT, DTOBS_MODEL_PLOT
     &    , MULTISEASON_OPTMASK, MULTISEASON_TGAP
     &    , MULTISEASON_NREJECT_OUTLIER, MULTISEASON_CHI2RED_ACTIVE
     &    , APPLY_FLUXCOR_BUG, PHOTFLAG_DETECT
     &    , CUTWIN_ZPADU, CUTWIN_ZPNPE

      COMMON / SNLCINP8 /
     &    H0_REF, OLAM_REF, OMAT_REF, W0_REF, DWDA_REF

      NAMELIST / SNLCINP /
     &      VERSION_PHOTOMETRY, VERSION_REPLACE, VERSION_REFORMAT_FITS
     &    , PRIVATE_DATA_PATH, FILTER_UPDATE_PATH
     &    , NONSURVEY_FILTERS, SNRMAX_FILTERS, VPEC_ERR_OVERRIDE
     &    , FILTER_REPLACE, FILTLIST_LAMSHIFT
     &    , JOBSPLIT, JOBSPLIT_EXTERNAL, SIM_PRESCALE, MXLC_FIT
     &    , OPTSIM_LCWIDTH, OPT_REFORMAT_TERSE, OPT_REFORMAT_TEXT
     &    , OPT_REFORMAT_SALT2, REFORMAT_KEYS, OPT_REFORMAT_FITS
     &    , SNMJD_LIST_FILE, SNMJD_OUT_FILE, MNFIT_PKMJD_LOGFILE
     &    , hfile_out, rootfile_out, textfile_prefix
     &    , SNTABLE_LIST, SNTABLE_FILTER_REMAP
     &    , hfile_kcor, KCOR_FILE, ZPOFF_FILE, USERTAGS_FILE
     &    , VPEC_FILE, HEADER_OVERRIDE_FILE
     &    , EPOCH_IGNORE_FILE, SNCID_LIST_FILE
     &    , SIMLIB_OUT, SIMLIB_ZPERR_LIST, NONLINEARITY_FILE
     &    , NFIT_ITERATION, MINUIT_PRINT_LEVEL, INTERP_OPT, USE_MINOS
     &    , OPT_SETPKMJD, MJDRANGE_SETPKMJD, OPT_DEBUG
     &    , LDMP_SNFAIL, LSIM_SEARCH_SPEC, LSIM_SEARCH_ZHOST
     &    , LTEST_KCOR, LTEST_INTERP, LTEST_MAG, LTEST_U3BAND
     &    , USESIM_SNIA, USESIM_SNCC, USESIM_TRUEFLUX, USE_SNHOST_ZPHOT
     &    , ABORT_ON_NOEPOCHS, ABORT_ON_BADAVWARP, ABORT_ON_NOPKMJD
     &    , ABORT_ON_BADZ, ABORT_ON_BADKCOR, ABORT_ON_BADSURVEY
     &    , ABORT_ON_MARGPDF0, ABORT_ON_TRESTCUT
     &    , ABORT_ON_DUPLCID, ABORT_ON_DUPLMJD
     &    , ABORT_ON_BADHEADER, ABORT_ON_BADFILTER
     &    , H0_REF, OLAM_REF, OMAT_REF, W0_REF, DWDA_REF
     &    , USE_MWCOR, DOBUG_LAMRANGE
     &    , MXLC_PLOT, NCCID_PLOT, SNCCID_PLOT, DTOBS_MODEL_PLOT
     &    , MULTISEASON_OPTMASK, MULTISEASON_TGAP
     &    , MULTISEASON_NREJECT_OUTLIER, MULTISEASON_CHI2RED_ACTIVE
     &    , APPLY_FLUXCOR_BUG, PHOTFLAG_DETECT
c
c Below are cut-variables defined in SNCUTS
c
     &   ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST
     &   ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE
     &   ,SNTYPE_IGNORE, CCDNUM_LIST, SIM_TEMPLATE_INDEX_LIST
     &   ,SNFIELD_LIST, SNTEL_LIST
     &   ,PHOTFLAG_MSKREJ, cutwin_cid
     &   ,cutwin_redshift, cutwin_redshift_err
     &   ,cutwin_ra, cutwin_dec
     &   ,cutwin_hostsep,   cutwin_Nepoch
     &   ,cutwin_snrmax,    cutwin_snrmax2, cutwin_nfield
     &   ,cutwin_mwebv,     cutwin_nseason_active
     &   ,cutwin_Trestmin,   cutwin_Trestmax
     &   ,cutwin_TrestRange, cutwin_Trest
     &   ,cutwin_Tgapmax,  cutwin_T0gapmax
     &   ,cutwin_Tobsmin,  cutwin_Tobsmax
     &   ,cutwin_peakmjd,  cutwin_mjd
c
     &   ,cutwin_psf,       cutwin_zp,  CUTWIN_ZPADU, CUTWIN_ZPNPE
     &   ,cutwin_moondist,  cutwin_photprob
     &   ,cutwin_nband_thresh
     &   ,cutwin_nfilt_snrmax,   cutwin_nfilt_snrmax2
     &   ,cutwin_nfilt_trestmin, cutwin_nfilt_trestmax
     &   ,cutwin_trest2, cutwin_nfilt_trest2
     &   ,cutwin_restlam, cutwin_lamrest, cutwin_lamobs
     &   ,PRIVATE_CUTWIN_STRING, PRIVATE_REDSHIFT_CMB
     &   ,SIMVAR_CUTWIN_STRING
     &   ,EARLYLC_STRING, REQUIRE_EPOCHS_STRING, DUMP_STRING
c
     &   ,SNCUT_SNRMAX, SNCUT_HOST_SBFLUX
     &   ,EPCUT_SNRMIN, EPCUT_SKYSIG, EPCUT_PSFSIG
     &   ,MAGOBS_SHIFT_PRIMARY,  MAGOBS_SHIFT_ZP
     &   ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_ZP, FILTER_LAMSHIFT
     &   ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS
     &   ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX
     &   ,FUDGE_MAG_ERROR, FUDGE_MAG_COVERR
     &   ,FLUXERRMODEL_FILE,SIM_FLUXERRMODEL_FILE,FLUXERRMODEL_OPTMASK
     &   ,MAGCOR_FILE, SIM_MAGCOR_FILE, FUDGE_HOSTNOISE_FILE
     &   ,RV_MWCOLORLAW, OPT_MWCOLORLAW, OPT_MWEBV
     &   ,MWEBV_SCALE, MWEBV_SHIFT
     &   ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR
     &   ,CUTWIN_MJD_EXCLUDE, CUTWIN_MJD_INCLUDE, CUTWIN_NMJD_INCLUDE


c ---------- end of SNLCINP ---------


C +CDE,FILTCOM. inserted below by ypatchy. 

c filter bandpasses.
c Jan 16 2017: MXLAMBIN_FILT -> 4000 (was 3000) for JWST test

      INTEGER MXLAMBIN_FILT, MXLAMBIN_PRIM
      PARAMETER (
     &   MXLAMBIN_FILT = 4000   ! max number of lambda bins to store
     &  ,MXLAMBIN_PRIM = 4000  ! max lambda bins for primary
     &    )

c Define filter chars A-Z, a-z and 0-9, but do not change order of
c original FILTDEF characters so that we can still use older Kcor files.
c
      CHARACTER FILTDEF_STRING*100
      PARAMETER (FILTDEF_STRING =
     &  'ugrizYJHK UBVRIXy0123456789 abcdef ' //  !  32
     &  'ACDEFGLMNOPQSTWZ ' //                    ! +16
     &  'hjklmnopqstvwx '   //                    ! +14 = 62
     &  '&'                             ! +1(autoSynFilter for spectrograph)
     &         )

      CHARACTER
     &   FILTOBS_NAME(MXFILT_ALL)*40   ! full filternames
     &  ,FILTREST_NAME(MXFILT_ALL)*40  ! full filternames
     &  ,PRIMARY_NAME*40               ! name of primary
     &  ,NONSURVEY_FILTERS_ADD*80      ! NONSURVEY_FILTERS that were added

      LOGICAL
     &   LFILTDEF_OBS(MXFILT_ALL)  ! filter-trans define in KCOR file
     &  ,LFILTDEF_REST(MXFILT_ALL) ! filter-trans define in KCOR file
     &  ,LFILTDEF_NONSURVEY(MXFILT_ALL) ! part of NONSURVEY_FILTERS
     &  ,LFILTDEF_SNRMAX(MXFILT_ALL)    ! use for SNRMAX cut
     &  ,FREEZE_SURVEY_FILTERS     ! T=> do not re-read on 2nd version
     &  ,EXIST_BXFILT_OBS          ! T=> observer BX exists
     &  ,EXIST_BXFILT_REST         ! T=> idem for rest-frame

c define indices for legacy filters
      INTEGER
     &   IFILT_SDSS_u, IFILT_SDSS_g, IFILT_SDSS_r
     &  ,IFILT_SDSS_i, IFILT_SDSS_z
     &  ,IFILT_BESS_U, IFILT_BESS_B, IFILT_BESS_V
     &  ,IFILT_BESS_R, IFILT_BESS_I, IFILT_BESS_BX
     &  ,IFILT_Y, IFILT_J, IFILT_H, IFILT_K

c define filter set for survey (observer) and for rest-frame.
c Each set goes 1 - NFILTDEF_[SURVEY,REST]

      INTEGER
     &   NFILTDEF_SURVEY          ! no. survey (obs) filters
     &  ,NFILTDEF_READ            ! no. filters to read (excludes BX)
     &  ,IFILTDEF_MAP_SURVEY(MXFILT_OBS)  ! IFILT_OBS vs. sparse filter index
     &  ,IFILTDEF_INVMAP_SURVEY(MXFILT_ALL) ! sparse index vs. IFILT_OBS
c
     &  ,NFILTDEF_REST                 ! no. rest-filts defined in Kcor file
     &  ,IFILTDEF_MAP_REST(MXFILT_OBS) ! IFILT_REST vs. sparse filter index
     &  ,IFILTDEF_INVMAP_REST(MXFILT_ALL) ! sparse index vs. IFILT_REST
c
     &  ,NFILTDEF_IGNORE_REST
     &  ,IFILTDEF_IGNORE_REST(MXFILT_OBS)  ! nearest rest-filters to ignore

c variables used in fitter; these variables are over-written
c in FITPAR_PREP for each SN

      INTEGER
     &   IFILT_REST_MAP(MXFILT_OBS)  ! idem for rest-frame filter
     &  ,NFILT_OBS_USEFIT                ! number of filters used in fit

      INTEGER*8
     &  IFILT_OBS_EVAL_MASK(2,MXFILT_ALL) ! set for each filter

      character FILTLIST_FIT_USE*64  ! filter-list USED each fit

c define filter properties using MXFILT_ALL
c Nov 12, 2010: split FILT_XXX into FILTOBS_XXX and FILTREST_XXX
c
      INTEGER
     &   NLAMBIN_FILTOBS(MXFILT_ALL)   ! number of bins per transmission curve
     &  ,NLAMBIN_FILTREST(MXFILT_ALL)  ! number of bins per transmission curve
     &  ,NLAMBIN_PRIMARY               ! lambda bins for primary spec


      REAL
     &   FILTOBS_TRANS(MXLAMBIN_FILT,MXFILT_ALL)   ! filter transmissions
     &  ,FILTOBS_TRANSMAX(MXFILT_ALL)              ! max trans
     &  ,FILTOBS_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)  ! corresponding lambda's (A)
     &  ,FILTOBS_LAMAVG(MXFILT_ALL)           ! effective central lambda
     &  ,FILTOBS_LAMRMS(MXFILT_ALL)           ! lambda-RMS
     &  ,FILTOBS_LAMRANGE(2,MXFILT_ALL)       ! min,max range in obs-frames
     &  ,FILTOBS_MAG_PRIMARY(MXFILT_ALL)      ! primary mag vs. ifilt_obs
     &  ,FILTOBS_ZPOFF_PRIMARY(MXFILT_ALL)    ! mag(native) - mag(synth)
     &  ,FILTOBS_ZPOFF_SNPHOT(MXFILT_ALL)     ! apply these ZPOFF to SNphot
     &  ,FILTOBS_LAMSHIFT(MXFILT_ALL)         ! lambda shift per filter
c
     &  ,FILTREST_TRANS(MXLAMBIN_FILT,MXFILT_ALL)   ! filter transmissions
     &  ,FILTREST_TRANSMAX(MXFILT_ALL)              ! max trans
     &  ,FILTREST_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)  !  lambda's (A)
     &  ,FILTREST_LAMAVG(MXFILT_ALL)           ! effective central lambda
     &  ,FILTREST_LAMRMS(MXFILT_ALL)           ! lambda-RMS
     &  ,FILTREST_LAMRANGE(2,MXFILT_ALL)       ! min,max range in rest-frames
     &  ,FILTREST_MAG_PRIMARY(MXFILT_ALL)      ! primary mag vs. ifilt
     &  ,FILTREST_ZPOFF_PRIMARY(MXFILT_ALL)    ! mag(native) - mag(synth)
     &  ,FILTREST_LAMSHIFT(MXFILT_ALL)         ! lambda shift per filter
c
     &  ,PRIMARY_FLUX(MXLAMBIN_PRIM)       ! primary spec for obs filters
     &  ,PRIMARY_LAM(MXLAMBIN_PRIM)


c for remap
      INTEGER NFILT_REPLACE, IFILTOBS_REPLACE(MXFILT_ALL)

      COMMON / FILTCOM / FILTOBS_NAME, FILTREST_NAME
     &   ,IFILT_SDSS_u, IFILT_SDSS_g, IFILT_SDSS_r
     &   ,IFILT_SDSS_i, IFILT_SDSS_z
     &   ,IFILT_BESS_U, IFILT_BESS_B, IFILT_BESS_V
     &   ,IFILT_BESS_R, IFILT_BESS_I, IFILT_BESS_BX
     &   ,IFILT_Y, IFILT_J, IFILT_H, IFILT_K
c
     &   ,NFILTDEF_SURVEY,        NFILTDEF_REST, NFILTDEF_READ
     &   ,FREEZE_SURVEY_FILTERS
     &   ,EXIST_BXFILT_OBS,       EXIST_BXFILT_REST
     &   ,IFILTDEF_MAP_SURVEY,    IFILTDEF_MAP_REST
     &   ,IFILTDEF_INVMAP_SURVEY, IFILTDEF_INVMAP_REST
     &   ,NFILTDEF_IGNORE_REST,   IFILTDEF_IGNORE_REST
     &   ,LFILTDEF_OBS, LFILTDEF_REST
     &   ,LFILTDEF_NONSURVEY, NONSURVEY_FILTERS_ADD
     &   ,LFILTDEF_SNRMAX
c
     &   ,NLAMBIN_FILTOBS, FILTOBS_TRANS, FILTOBS_TRANSMAX
     &   ,FILTOBS_LAMBDA, FILTOBS_LAMAVG, FILTOBS_LAMRMS
     &   ,FILTOBS_LAMRANGE
     &   ,FILTOBS_MAG_PRIMARY, FILTOBS_ZPOFF_PRIMARY
     &   ,FILTOBS_ZPOFF_SNPHOT, FILTOBS_LAMSHIFT
c
     &   ,NLAMBIN_FILTREST, FILTREST_TRANS, FILTREST_TRANSMAX
     &   ,FILTREST_LAMBDA,  FILTREST_LAMAVG, FILTREST_LAMRMS
     &   ,FILTREST_MAG_PRIMARY, FILTREST_ZPOFF_PRIMARY
     &   ,FILTREST_LAMSHIFT, FILTREST_LAMRANGE
c
     &   ,IFILT_REST_MAP
     &   ,FILTLIST_FIT_USE, NFILT_OBS_USEFIT
     &   ,PRIMARY_FLUX, PRIMARY_LAM, NLAMBIN_PRIMARY, PRIMARY_NAME
c
     &   ,NFILT_REPLACE, IFILTOBS_REPLACE

      COMMON / FILTCOM8 / IFILT_OBS_EVAL_MASK


      INTEGER
     &   i, LL, IERR
     &  ,IZ, IT, IFILE, NFILTDEF, ISTAT

      REAL*8 Z8, T8, AVWARP8, KCOR8, ERR8

      character cfilt_rest*1, cfilt_obs*1

c function
      REAL*8  GET_KCOR8, GET_KCORERR8, GET_AVWARP8

C --------------- BEGIN -----------------


c init user args.

      AVWARP8 = DBLE(0.0)

      REDSHIFT_DMP   = -9.0
      TREST_DMP      = -99.0
      IFILT_REST_DMP = -9
      IFILT_OBS_DMP  = -9
      KCOR_FILE      = ''
      IFILT_RESTCOLOR_DMP(1) = -9
      IFILT_RESTCOLOR_DMP(2) = -9

c read command line args.

      NLINE_ARGS = IARGC()
      DO i = 1, NLINE_ARGS
        CALL GETARG(i, LINE_ARGS(i) )  ! read command-line args
        LL = INDEX ( LINE_ARGS(i), ' ' )  - 1
        USE_LINE_ARGS(i) = .FALSE.
        print*,'  COMMAND LINE ARG = ', LINE_ARGS(i)(1:LL)
      ENDDO

      CALL INIT_SNVAR(IERR)

c --------------------------------------------------
C parse command line args
      CALL PARSE_KCORDUMP_ARGS

c --------------------------------------------------
c pass filter info to snana before reading kcor file
      NFILTDEF = 1
      CALL SET_SURVEY("DUMMY", NFILTDEF, IFILT_OBS_DMP )

c -------------------------------------
c read filters and K-corrections

      CALL RDKCOR ( KCOR_FILE, IERR )

c --------------------------------
c  if redshift is 'ALL', then use z-range in kcor file

      IF ( NZBIN_DMP .GT.1 ) THEN
        CALL SET_REDSHIFT_ALL
      ENDIF

      IF ( IFILT_RESTCOLOR_DMP(2) > 0 ) THEN
         print*,' '
         print*,'   RestFrame ', STRING_REST_COLOR,' = ', MAGREST_COLOR
      ENDIF
c -----------------------------
      print*,' '
      print*,' BEGIN KCOR DUMP: '
      print*,' '

      cfilt_rest = filtdef_string(ifilt_rest_dmp:ifilt_rest_dmp)
      cfilt_obs  = filtdef_string(ifilt_obs_dmp:ifilt_obs_dmp)

      write(6,15) cfilt_rest, cfilt_obs
15    format(T8, 'REDSHIFT   Trest   AVwarp    K_',A,A,'    ERROR'  )
      print*,'# -----------------------------------------------------'

      DO iz = 1, NZBIN_DMP

         Z8 = REDSHIFT_DMP_RANGE(1) +
     &        REDSHIFT_DMP_BINSIZE * float(iz-1)

      DO it = 1, NTBIN_DMP

         T8 = TREST_DMP_RANGE(1) +
     &        TREST_DMP_BINSIZE * float(it-1)

         IF( IFILT_RESTCOLOR_DMP(1) > 0 ) THEN
           AVWARP8  = GET_AVWARP8(T8, Zat10pc     ! (I)
     &             ,MAGREST(1)
     &             ,MAGREST(2)                   ! (I)
     &             ,IFILT_RESTCOLOR_DMP(1)
     &             ,IFILT_RESTCOLOR_DMP(2)        ! (I)
     &             ,istat )                       ! (O)
        ENDIF

        KCOR8 = GET_KCOR8 ( ifilt_Rest_dmp, ifilt_obs_dmp,
     &         T8, Z8, AVWARP8 )

        ERR8 = GET_KCORERR8 ( OPT_KCORERR_SJ,
     &        ifilt_Rest_dmp, ifilt_obs_dmp, T8, Z8, AVWARP8 )

        write(6,20) Z8, T8, AVwarp8, KCOR8, ERR8
20      format('KCOR: ', F8.3, F10.2, F8.3,  2F10.5 )

      ENDDO
      ENDDO

      STOP
      END

C ===========================
CDECK  ID>,  PRSKARGS. 
      SUBROUTINE PARSE_KCORDUMP_ARGS()
      IMPLICIT NONE

C +CDE,DUMPCOM. inserted below by ypatchy. 

      REAL
     &   REDSHIFT_DMP_BINSIZE
     &  ,TREST_DMP_MIN
     &  ,TREST_DMP_MAX
     &  ,TREST_DMP_BINSIZE

      PARAMETER (
     &   REDSHIFT_DMP_BINSIZE = 0.01
     &  ,TREST_DMP_MIN        = -20.0
     &  ,TREST_DMP_MAX        = +80.0
     &  ,TREST_DMP_BINSIZE    =   1.0
     &    )


      CHARACTER
     &   STRING_Kxy*8, STRING_REST_COLOR*4  ! Apr 2016

      INTEGER
     &   IFILT_REST_DMP       ! (I) rest-frame filter index
     &  ,IFILT_OBS_DMP        ! (I) observer frame filter index
     &  ,IFILT_RESTCOLOR_DMP(2) ! (I) rest-color (default=0.0) (Apr 2016)

      REAL
     &   Trest_DMP       ! (I) rest-frame epoch, days
     &  ,redshift_DMP    ! (I) redshift
     &  ,REDSHIFT_DMP_RANGE(2)
     &  ,TREST_DMP_RANGE(2)
     &  ,MAGREST_COLOR, MAGREST(2)

      INTEGER
     &   NZBIN_DMP, NTBIN_DMP

      COMMON / DUMPCOM /
     &   STRING_Kxy,     STRING_REST_COLOR, MAGREST_COLOR, MAGREST
     &  ,IFILT_REST_DMP, IFILT_OBS_DMP, IFILT_RESTCOLOR_DMP
     &  ,Trest_DMP,      redshift_DMP
     &  ,TREST_DMP_RANGE,   REDSHIFT_DMP_RANGE
     &  ,NZBIN_DMP,         NTBIN_DMP

C ==============================================

C +CDE,SNDATCOM. inserted below by ypatchy. 


C +CDE,SNPAR. inserted below by ypatchy. 

      CHARACTER  SNTABLE_LIST_DEFAULT*60

c parameters used by snana code

      INTEGER
     &   MXVERS, MXSURVEY, MXSNLC,MXCID, MXCID_CHECK,MXEPOCH, MXITER
     &  ,MXFILT_ALL, MXFILT_OBS, MXFILT_REST, ONE
     &  ,MXFILT_KCOR, MXFILT_SNRMAX, MXEP_MODELGRID
     &  ,MXIDTEL, MXIDFIELD, MXFIELD_OVP, MXSEASON, MXSNHOST
     &  ,MXTABLE1D_KCOR,MXTABLE1D_LCMAG,MXTABLE1D_MWXT,MXTABLE1D_AVWARP
     &  ,MXZBIN_KCOR, MXTBIN_KCOR, MXAVBIN_KCOR, MXCBIN_AVWARP
     &  ,MNTYPE, MXTYPE, MXLISTNML, MXCCID_LIST
     &  ,MXVAR_PRIVATE, MXCUT_PRIVATE, MXVAR_TERSE
     &  ,HOFF_NOCUTS, HOFF_CUTS, HOFF_SIM
     &  ,LUNHIS, LUNNML, LUNDAT, LUNTMP, LUNFIT, LUNPKMJD, LUNDMP
     &  ,LUNCID, LUNOUT
     &  ,LUNRES1, LUNRES2, LUNINTERP, LUNLIST, LUNIGNORE, LUNSALT2
     &  ,ISTAGE_INIT, ISTAGE_RDSN, ISTAGE_CUTS, ISTAGE_USRANA
     &  ,ISTAGE_TEST
     &  ,INTERP_LINEAR, INTERP_SMOOTH, INTERP_ZSMOOTH
     &  ,MXERRTYPE, ERRTYPE_MINOS, ERRTYPE_PARAB, ERRTYPE_MARG
     &  ,ERRTYPE_BAD
     &  ,FCNFLAG_USER, FCNFLAG_FAST, FCNFLAG_LAST, FCNFLAG_USESIM
     &  ,FCNFLAG_PRIOR_ONLY, FCNFLAG_SIGMA_ONLY
     &  ,OPT_INTEGPDF_QUITCHI2, OPT_INTEGPDF_FULL
     &  ,OPT_SNXT_CCM89, OPT_SNXT_SJPAR
     &  ,OPT_MWCOLORLAW_DEFAULT, OPT_MWEBV_DEFAULT
     &  ,OPT_KCORERR_SJ, OPT_KCORERR_SJ5, OPT_KCORERR_SMOOTH
     &  ,MXLINE_ARGS, MXEPOCH_IGNORE
     &  ,NPAR_ANYLC, MXPAR_SIMSED, MXPAR_LCLIB
     &  ,MXLAMBIN_SNSED, MXCUTBIT
     &  ,OPT_FILTUPD_EACHSN, OPT_FILTUPD_MAP, OPT_FILTUPD_SAMEFILT
     &  ,OPT_FILTOBS, OPT_FILTREST
     &  ,ISTAT_READAGAIN, ISTAT_SKIP
     &  ,IDTABLE_MCMC, NFIT_VERBOSE
     &  ,ITABLE_SNANA,  ITABLE_FITRES, ITABLE_SNLCPAK, ITABLE_SPECPAK
     &  ,ITABLE_SPECTRA, MXTABLE
     &  ,IDTABLE_SNANA, IDTABLE_FITRES, IDTABLE_SPECTRA
     &  ,IDTABLE_CIDPTR
     &  ,OPT_PARSTORE_TEXTTABLE
c
     &  ,MODEL_STRETCH
     &  ,MODEL_STRETCH2
     &  ,MODEL_MLCS2k2
     &  ,MODEL_SNOOPY
     &  ,MODEL_SALT2
     &  ,MODEL_SIMSED
     &  ,MODEL_BYOSED
     &  ,MODEL_NON1A
     &  ,MODEL_LCLIB    ! Sep 2017
     &  ,MODEL_FIXMAG   ! force mags to user-value
     &  ,MXMODEL_INDEX
c
     &  ,MXCHAR_CCID, MXCHAR_VERSION, MXCHAR_SURVEY
     &  ,MXCHAR_PATH, MXCHAR_FILENAME, MXCHAR_MODELNAME
     &  ,MXCHAR_FIELDNAME, MXCHAR_PARNAME, MXCHAR_CUTNAME
     &  ,MXCHAR_FILEWORD
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA     ! DATA flux per epoch
     &  ,SNLCPAK_EPFLAG_REJECT   ! REJECT flag (1=> excluded from fit)
     &  ,SNLCPAK_EPFLAG_CHI2     ! data-fit chi2 per epoch
     &  ,SNLCPAK_EPFLAG_FITFUN   ! smooth  fitfun curve
     &  ,SNLCPAK_EPFLAG_FLUXSIM  ! SIM flux per epoch
     &  ,SNLCPAK_EPFLAG_FLUXREST ! rest-frame flux per epoch (optional)
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST ! idem for sim truth
     &  ,SNLCPAK_BANDFLAG_PKFLUX   ! peak flux vs. filter
     &  ,SNLCPAK_BANDFLAG_PKMJD    ! peak MJD  vs. filter
     &  ,SNLCPAK_BANDFLAG_NDOF     ! Ndof vs. filter
     &  ,SNLCPAK_BANDFLAG_CHI2     ! chi2 vs. filter
c
     &  ,IFLAG_INI, IFLAG_ANA, IFLAG_END
     &  ,MAG_SATURATE
     &  ,MASK_FLUXCOR_SNANA, MASK_FLUXERRCOR_SNANA
c
     &  ,EXIT_ERRCODE_SNANA, EXIT_ERRCODE_SNFIT, EXIT_ERRCODE_PSNID

      REAL*8
     &   ZERO8, ONE8, TEN8,  PI, LOGTEN
     &  ,PDFMIN, PDFMAX_EDGE, PDFMIN_GOOD
     &  ,KCORPACK
     &  ,XTMW_FRACERR
     &  ,MJDOFF
     &  ,ZEROPOINT_FLUXCAL_DEFAULT
     &  ,RV_MWCOLORLAW_DEFAULT
     &  ,CUTVAL_OPEN, IGNORE_HEADVAL

      REAL NULLVAL

      PARAMETER (                   ! BEGIN SNANA PARAMS
     &   SNTABLE_LIST_DEFAULT = 'SNANA  FITRES  LCPLOT'
     &  ,MXVERS        = 50         ! max number of versions to read
     &  ,MXSURVEY      = 100        ! max number of SURVEY_NAMEs to store
     &  ,MXITER        = 12         ! max # fit iterations
     &  ,NFIT_VERBOSE  = 500        ! Number of fits for verbose printouts
     &  ,MXSNLC        = 4000000    ! max number of SNe (fits format)
     &  ,MXCCID_LIST   = 10000      ! max size of SNCCID_LIST_ALL
     &  ,MXCID       = 299 999 999  ! max CID = 300 million -1 (9 digits)
     &  ,MXCID_CHECK =  99 999 999  ! MXCID to check duplicates
     &  ,MXEPOCH     = 2000        ! max number of filter-epochs per SN
     &  ,MXEP_MODELGRID = 200     ! max model grid for SNANA+SIM_MAGOBS table
     &  ,MXIDTEL     = 200        ! max telescope ID in SURVEY.DEF
     &  ,MXIDFIELD   = 200        ! max FIELD ID in SURVEY.DEF
     &  ,MXFIELD_OVP = 12         ! max number of overlapping fields
     &  ,MXSEASON    = 100        ! max number of seasons
     &  ,MXSNHOST    = 2          ! max number of host matches to read/write
     &  ,MXVAR_PRIVATE = 40       ! max number of private variables
     &  ,MXCUT_PRIVATE = 10        ! max number of cuts on private var
     &  ,MXVAR_TERSE   = 30       ! max number of text-data  columms
     &  ,MXLISTNML   = 52    ! max list size for some NML lists
     &  ,MXFILT_OBS  = 62    ! max number of used observer filters
     &  ,MXFILT_ALL  = 80    ! max number of all possible filter defs
     &  ,MXFILT_REST = 20    ! max number of rest-frame filter types (was 12)
     &  ,MXFILT_KCOR = MXFILT_REST ! max number of filters used in Kcor table
     &  ,MXFILT_SNRMAX = 10       ! no more than this many SNRMAX(filt) cuts
     &  ,MNTYPE       =   1       ! min sn "type"
     &  ,MXTYPE       = 1000      ! max sn "type"
     &  ,MXZBIN_KCOR  = 100       ! max # Z-bins for KCOR tables
     &  ,MXTBIN_KCOR  = 150       ! max # Epochs for KCOR tables
     &  ,MXAVBIN_KCOR = 100       ! max # AV bins for KCOR tables
     &  ,MXCBIN_AVWARP = 100      ! max color index for AVWARP table
     &  ,MXTABLE1D_KCOR   = 8 000 000  ! max number of KCOR bins
     &  ,MXTABLE1D_AVWARP = 1 000 000  ! max number of bins for AVWARP table
     &  ,MXTABLE1D_LCMAG  = 2 000 000  ! idem for LCMAG  table
     &  ,MXTABLE1D_MWXT   = 2 000 000  ! idem for MWXT table
     &  ,KCORPACK      = 1000.  ! store KCOR * KCORPACK as I*2
     &  ,XTMW_FRACERR  = 0.16   ! error on MW Xtinc is 16% of XTMW
c
     &  ,HOFF_NOCUTS = 100
     &  ,HOFF_CUTS   = 200
     &  ,HOFF_SIM    = 20000
     &  ,LUNHIS    = 21  ! LUN for hbook output file
     &  ,LUNNML    = 22  ! LUN for input namelist
     &  ,LUNDAT    = 23  ! LUN for data
     &  ,LUNTMP    = 24
     &  ,LUNOUT    = 25
     &  ,LUNFIT    = 26
     &  ,LUNLIST   = 27
     &  ,LUNIGNORE = 28
     &  ,LUNDMP    = 29
     &  ,LUNRES1   = 31  ! for DMP_FITRES
     &  ,LUNRES2   = 32
     &  ,LUNINTERP = 33  ! for FLUX,MAGs interpolated at SN,MJD
     &  ,LUNSALT2  = 34  ! for SALT2 dictFile
     &  ,LUNPKMJD  = 35
     &  ,LUNCID    = 36  ! reserved for reading SNCID_LIST_FILE
     &  ,ISTAGE_INIT    = 10       ! init has finished
     &  ,ISTAGE_RDSN    = 20       ! SN have been read
     &  ,ISTAGE_CUTS    = 30       ! cuts have been applied
     &  ,ISTAGE_USRANA  = 40       ! USRANA has been called.
     &  ,ISTAGE_TEST    = 90       ! set for testing
     &  ,ZERO8          = 0.0
     &  ,ONE8           = 1.0
     &  ,ONE            = 1
     &  ,TEN8           = 10.0
     &  ,LOGTEN         = 2.302585092994
     &  ,INTERP_LINEAR  = 1    ! option flag to use linear DFINT
     &  ,INTERP_SMOOTH  = 2    ! option flag to use smoothing
     &  ,INTERP_ZSMOOTH = 3    ! linear DFINT, but smooth along z
c
     &  ,MXERRTYPE      = 10
     &  ,ERRTYPE_MINOS  = 1
     &  ,ERRTYPE_PARAB  = 2
     &  ,ERRTYPE_MARG   = 3  ! marginalized error
     &  ,ERRTYPE_BAD    = 6
C
     &  ,FCNFLAG_LAST     = 3    ! last MINUIT call
     &  ,FCNFLAG_USER     = 90   ! pass this IFLAG for user-FCNSNLC calls
     &  ,FCNFLAG_FAST     = 91   ! go as fast as possible (no LAST if-block)
     &  ,FCNFLAG_PRIOR_ONLY = 92
     &  ,FCNFLAG_SIGMA_ONLY = 93
     &  ,FCNFLAG_USESIM   = 99   ! pass this IFLAG to use SIM params
c
     &  ,PI     = 3.1415926535898
     &  ,PDFMIN      = 1.0E-5   ! used to speed up PDF integration
     &  ,PDFMAX_EDGE = 0.03     ! max allowed PDF value at edges
     &  ,PDFMIN_GOOD = 1.0E-4   ! PDF > PDFMIN_GOOD counts as good point
     &  ,OPT_INTEGPDF_QUITCHI2 = 2  ! abort FCNSNLC if chi2 > quitchi2
     &  ,OPT_INTEGPDF_FULL     = 1  ! do full FCNSNLC evaluation
c
     &  ,OPT_SNXT_CCM89   = 1  ! exact SN etinction using INIT_XTHOST
     &  ,OPT_SNXT_SJPAR   = 2  ! SN extinction with Jha's parameters
     &  ,OPT_KCORERR_SMOOTH = 1 ! use smooth half-Gaussian
     &  ,OPT_KCORERR_SJ     = 2   ! use Saurabh's Kcor error
     &  ,OPT_KCORERR_SJ5    = 5  ! x5 Saurabh's Kcor error
c
     &  ,RV_MWCOLORLAW_DEFAULT  = 3.1   ! A_V/E(B-V)
     &  ,OPT_MWCOLORLAW_DEFAULT = 94    ! ODonnel 94
     &  ,OPT_MWEBV_DEFAULT      =  1    ! whatever is in the data file
c
     &  ,MXLINE_ARGS      = 100
     &  ,ZEROPOINT_FLUXCAL_DEFAULT  = 27.5
     &  ,MXEPOCH_IGNORE   = 1000
     &  ,NPAR_ANYLC       = 8    ! for MNFIT_PKMJD
     &  ,MXPAR_SIMSED     = 100  ! max number of SIMSED parameters
     &  ,MXPAR_LCLIB      = 40   ! should be same as in genmag_LCLIB.h
     &  ,MXLAMBIN_SNSED   = 4000 ! Jan 2017: raised from 3000
     &  ,MXCUTBIT         = 64   ! max number of cut bits
c
     &  ,OPT_FILTUPD_EACHSN   = 1  ! default filter-updates
     &  ,OPT_FILTUPD_MAP      = 2  ! default filter-updates
     &  ,OPT_FILTUPD_SAMEFILT = 3  ! test: use same filter each SN
     &  ,OPT_FILTREST         = 1
     &  ,OPT_FILTOBS          = 2
     &  ,ISTAT_READAGAIN      = 7
     &  ,ISTAT_SKIP           = -1
c
     &  ,ITABLE_SNANA=1, ITABLE_FITRES=2, ITABLE_SNLCPAK=3
     &  ,ITABLE_SPECPAK=4, ITABLE_SPECTRA=5
     &  ,MXTABLE = 10
     &  ,IDTABLE_SNANA   = 7100
     &  ,IDTABLE_FITRES  = 7788
     &  ,IDTABLE_SPECTRA = 8000  ! SALT2 model spectra from LC fit
     &  ,IDTABLE_CIDPTR  = 1600
     &  ,IDTABLE_MCMC    = 7711
     &  ,OPT_PARSTORE_TEXTTABLE = 1  ! tag subset for TEXT table.
c
     &  ,MODEL_STRETCH    = 1
     &  ,MODEL_STRETCH2   = 2
     &  ,MODEL_MLCS2k2    = 3
     &  ,MODEL_SNOOPY     = 4
     &  ,MODEL_SALT2      = 6
     &  ,MODEL_SIMSED     = 7
     &  ,MODEL_BYOSED     = 8
     &  ,MODEL_NON1A      = 10
     &  ,MODEL_LCLIB      = 12
     &  ,MODEL_FIXMAG     = 20
     &  ,MXMODEL_INDEX    = 20
     &  ,NULLVAL          = -99999.
     &  ,IGNORE_HEADVAL   = -5555.
c
     &  ,MXCHAR_CCID       = 20   ! max len of CCID string (i.e, SN name)
     &  ,MXCHAR_VERSION    = 72   ! max len of VERSION_PHOTOMETRY
     &  ,MXCHAR_SURVEY     = 40   ! max len of SURVEY_NAME
     &  ,MXCHAR_PATH       = 160  ! max len of path
     &  ,MXCHAR_FILENAME   = 200  ! max len of filename with full path
     &  ,MXCHAR_MODELNAME  = 72   ! max len of model name
     &  ,MXCHAR_FIELDNAME  = 20   ! max len of field name
     &  ,MXCHAR_PARNAME    = 20   ! max len of parameter name
     &  ,MXCHAR_CUTNAME    = 160  ! to define cut names
     &  ,MXCHAR_FILEWORD   =  60  ! size of FILEWORD_LIST
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA    = 1    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_REJECT      = 2
     &  ,SNLCPAK_EPFLAG_CHI2        = 3
     &  ,SNLCPAK_EPFLAG_FITFUN      = 4
     &  ,SNLCPAK_EPFLAG_FLUXSIM     = 5    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_FLUXREST    = 6
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST = 7
     &  ,SNLCPAK_BANDFLAG_NDOF    = 100
     &  ,SNLCPAK_BANDFLAG_PKFLUX  = 101  ! filter-dependent
     &  ,SNLCPAK_BANDFLAG_PKMJD   = 102
     &  ,SNLCPAK_BANDFLAG_CHI2    = 103
c
     &  ,IFLAG_INI=1, IFLAG_ANA=2, IFLAG_END=3
     &  ,MAG_SATURATE = -7.0   ! for sim only
     &  ,CUTVAL_OPEN = 1.0E12  ! cutwin value to accept everything.
     &	,MASK_FLUXCOR_SNANA    = 1
     &  ,MASK_FLUXERRCOR_SNANA = 2
c
     &  ,EXIT_ERRCODE_SNANA = 21
     &  ,EXIT_ERRCODE_SNFIT = 22
     &  ,EXIT_ERRCODE_PSNID = 23
     &      )

c physical constants

      REAL*8
     &   PARSEC, CLIGHT, PEAKMAG_AT_10PC
     &  ,Zat10pc
     &  ,OMAT_DEFAULT, OMATERR_DEFAULT
     &  ,OLAM_DEFAULT, OLAMERR_DEFAULT
     &  ,ORAD_DEFAULT
     &  ,H0_DEFAULT,   H0ERR_DEFAULT
     &  ,W0_DEFAULT,   W0ERR_DEFAULT
     &  ,DWDA_DEFAULT, DWDAERR_DEFAULT

      PARAMETER (
     &   PARSEC        = 3.085678E13  ! 1 parsec (km)
     &  ,CLIGHT        = 2.998E5      ! c (km/sec)
c
     &  ,H0_DEFAULT        = 70.0 / ( 1.0E6 * Parsec )  ! standard value
     &  ,W0_DEFAULT        = -1.0
     &  ,DWDA_DEFAULT      =  0.0
     &  ,OMAT_DEFAULT      =  0.3
     &  ,OLAM_DEFAULT      =  0.7
     &  ,ORAD_DEFAULT      =  1.2E-5
c
     &  ,H0ERR_DEFAULT        =  7.0 / ( 1.0E6 * Parsec )
     &  ,W0ERR_DEFAULT        =  0.1
     &  ,DWDAERR_DEFAULT      =  0.0
     &  ,OMATERR_DEFAULT      =  0.03
     &  ,OLAMERR_DEFAULT      =  0.03
c
     &  ,PEAKMAG_AT_10PC   = -19.6
     &  ,Zat10pc           = 2.34E-9   ! magic redshift at 10 pc
     &  ,MJDOFF            = 0.0       ! 53000.
     &     )



C +CDE,SNFILECOM. inserted below by ypatchy. 

c Sep 9, 2010: Pulled out of SNDATACOM

c define HEADER MASK BITS for required variables.
      INTEGER
     &   HEADBIT_SNID,   HEADBIT_IAUC, HEADBIT_CIDSEL
     &  ,HEADBIT_SURVEY, HEADBIT_FILTERS
     &  ,HEADBIT_RA,    HEADBIT_DEC
     &  ,HEADBIT_MWEBV, HEADBIT_Z
     &  ,NBIT_HEADMASK
     &  ,HEADMASK_REQUIRED
c
     &  ,HEADMASK  ! reset for reach SN

      PARAMETER (
     &    HEADBIT_SNID    = 0  ! MASK = 1
     &   ,HEADBIT_IAUC    = 1  ! MASK = 2  ! added Dec 2015
     &   ,HEADBIT_CIDSEL  = 2  ! MASK = 4  ! added Dec 2015
     &   ,HEADBIT_SURVEY  = 3
     &   ,HEADBIT_FILTERS = 4
     &   ,HEADBIT_RA      = 5
     &   ,HEADBIT_DEC     = 6
     &   ,HEADBIT_MWEBV   = 7
     &   ,HEADBIT_Z       = 8
     &   ,NBIT_HEADMASK   = 9
     &   ,HEADMASK_REQUIRED  = 2**(NBIT_HEADMASK)-1 - 2 ! don't require IAUC
     &      )

      CHARACTER
     &   SNDATA_ROOT*(MXCHAR_PATH)
     &  ,SNANA_DIR*(MXCHAR_PATH)
     &  ,SNDATA_PATH*(MXCHAR_PATH)      ! subdir with data or sim files
     &  ,SNLIST_FILE*(MXCHAR_FILENAME)  ! input list of SNDATA Files
     &  ,SNREADME_FILE(MXVERS)*(MXCHAR_FILENAME)   ! name of EVERY README file
     &  ,SNDATA_FILE_CURRENT*(MXCHAR_FILENAME)     ! current file being read
     &  ,GLOBAL_BANNER*120
     &  ,SNDATA_PREFIX*(MXCHAR_FILENAME)  ! $SNDATA_ROOT/lcmerge/$VERSION
     &  ,C1ERR*88, C2ERR*88    ! generic error strings

      INTEGER
     &  ABSO_INDEX(MXSNLC) ! absolute (IFILE or IROW) vs. ISN index

      LOGICAL LFLAG_RDHEAD_ONLY

      COMMON / SNFILECOM /
     &   SNDATA_ROOT, SNANA_DIR, SNLIST_FILE
     &  ,SNDATA_FILE_CURRENT
     &  ,GLOBAL_BANNER, SNDATA_PREFIX, SNREADME_FILE, SNDATA_PATH
     &  ,C1ERR, C2ERR, ABSO_INDEX, HEADMASK, LFLAG_RDHEAD_ONLY



C +CDE,CTRLCOM. inserted below by ypatchy. 

c control variables and counters.

      INTEGER
     &   ISTAGE_SNANA       ! current stage of processing
     &  ,NACCEPT_CUT(MXCUTBIT) ! Number of SN that pass each cut
     &  ,NACCEPT_CID           ! # SN with valid CID
     &  ,NACCEPT_TYPE          ! # SN with valid TYPE
     &  ,NACCEPT_Z             ! idem with valid redshift
     &  ,NACCEPT_ZERR          ! idem with valid redshift error
     &  ,NMJD_IDTEL            ! # MJDs with valid telescope id
     &  ,APPLY_HEADER_CUTMASK  ! mask of applied header cuts
c
     &  ,JTIME_START
     &  ,JTIME_LOOPSTART        ! time at start of fits with TIME()
     &  ,JTIME_LOOPEND          ! time and end
     &  ,NCALL_SNANA_DRIVER
c
     &  ,NPASSCUT_INCREMENT(-1:MXTYPE,100)  ! 100 > NCUTBIT_SNLC
     &  ,NPASSCUT_FIT(-1:MXTYPE)
     &  ,NVAR_NEARNBR          ! Number of NEARNBR variables to analyze
     &  ,NSTORE_MAGCOR         ! number of stored MAGCOR values
     &  ,NUSE_MAGCOR           ! number of used MAGCOR values
     &  ,SIGN_MAGCOR           ! add or subtract
     &  ,FORCEMASK_FLUXCOR   ! mask to force fluxCor, even if already applied
     &  ,EXIT_ERRCODE        ! used for abort

      LOGICAL
     &   DO_FIT
     &  ,DO_FLUXERRCALC    ! T => compute error from PSF,SKY & ZPT
     &  ,LSIM_SNANA        ! simulated with SNANA
     &  ,LSIM_MAGOBS       ! data-like, but with SIM_MAGOBS
     &  ,ISJOB_SNANA          ! =T for snana.exe only.
     &  ,ISJOB_PSNID
     &  ,REFORMAT_SAVE_BADEPOCHS  ! set if bit2 is set on any OPT_REFORMAT
     &  ,STDOUT_UPDATE         ! T => update event to screen
     &  ,DOFUDGE_HOSTNOISE     ! T => FUDGE_HOSTNOISE_FILE is set
     &  ,DOFUDGE_NONLIN        ! T => NONLINEARITY_FILE is set
     &  ,DOFUDGE_FLUXERRMODEL  ! T => FLUXERRMODEL_FILE

      CHARACTER SNANA_VERSION*12

c global survey info
      CHARACTER
     &   SURVEY_NAME*(MXCHAR_SURVEY)
     &  ,SURVEY_NAME_LIST(MXSURVEY)*(MXCHAR_SURVEY) ! in SURVEY.DEF file
     &  ,SUBSURVEY_NAME*(MXCHAR_SURVEY)  ! e.g.,  SURVEY:  BLA(SUBSURVEY)
     &  ,SUBSURVEY_NAME_LIST*(MXCHAR_FILENAME) ! comma-sep list
     &  ,SURVEY_FILTERS*(MXFILT_ALL)  ! read from data file
     &  ,SURVEY_FIELDNAME(MXIDFIELD)*(MXCHAR_FIELDNAME) ! from SURVEY.DEF

      INTEGER
     &   NFIELD_SURVEY             ! number of survey fields in SURVEY.DEF
     &  ,SURVEY_IDFIELD(MXIDFIELD) ! integer ID for each field

      REAL
     &   ZEROPOINT_FLUXCAL(MXFILT_OBS)  ! defines calibrated flux

      INTEGER*4
     &   N_VERSION             ! number of photometry version to read
     &  ,N_SNLC                ! number of SN lightcurves read
     &  ,N_SNLC_CUTS           ! Number of SN after cuts (bookkeeping only)
     &  ,N_SNLC_FIT            ! Number of fitted SN
     &  ,N_SNLC_FITCUTS        ! Number of fitted SN after fit cuts
     &  ,N_SNLC_COVFIX         ! Number of SN with fixed COV to be invertible
     &  ,N_DUPLICATE_CID       ! Number of duplicate CIDs
     &  ,N_DUPLICATE_MJD       ! Number of duplicate MJD+BAND (Jun 2017)
     &  ,NSTORE_DUPLICATE_MJD  ! Number stored
     &  ,N_SNFILE              ! # of SN files to read per version
     &  ,N_SNFILE_LAST         ! idem, as of last version
     &  ,NEPOCH_TOT            ! total number of epochs read
     &  ,NEPOCH_CUT            ! total number of epochs passing cuts
     &  ,NEPOCH_BADPHOT        ! # epochs with bad PHOTFLAG (per event)
     &  ,NEPOCH_BADPHOT_SUM    ! # epochs with bad PHOTFLAG (summed)
     &  ,IDSURVEY              ! survey ID from SURVEY.DEF
     &  ,IDSUBSURVEY           ! =IDSURVEY unless subSurvey is different
     &  ,IDSURVEY_LIST(MXSURVEY) ! corresponds to SURVEY_NAME_LIST
     &  ,NSURVEY_LIST          ! size of SURVEY_NAME_LIST & IDSURVEY_LIST

      LOGICAL*1
     &    EXIST_FILT(MXFILT_OBS)  ! T => at least one point per filt
     &   ,FOUND_SURVEY
     &   ,FORMAT_TEXT    ! ascii/txt for input data
     &   ,FORMAT_TERSE   ! terse-text for input data
     &   ,FORMAT_VERBOSE ! verbose-text for ...
     &   ,FORMAT_FITS    ! snfitsio for input data.

      INTEGER
     &   NVAR_TERSE          ! Number variables for terse format
     &  ,ITERSE_MJD          ! index of MJD in TERSE variable list
     &  ,ITERSE_FILTER       ! idem for filter
     &  ,ITERSE_FIELD
     &  ,ITERSE_FLUXCAL
     &  ,ITERSE_FLUXCALERR
     &  ,ITERSE_PHOTFLAG
     &  ,ITERSE_PHOTPROB
     &  ,ITERSE_ZPFLUX     ! mag <-> native flux zeropoint
     &  ,ITERSE_PSFSIG     ! added July 18, 2013
     &  ,ITERSE_SKYSIG     ! idem
     &  ,ITERSE_SKYSIG_T   ! added Aug 7 2014
     &  ,ITERSE_XPIX       ! added Aug 7 2014
     &  ,ITERSE_YPIX       ! added Aug 7 2014
     &  ,ITERSE_GAIN       ! added Oct 2015 to compute FLUXCAL_ERRCALC
     &  ,ITERSE_MAGOBS
     &  ,ITERSE_CCDNUM     ! added Oct 16 2017
     &  ,ITERSE_SIM_EPFILTREST
     &  ,ITERSE_SIM_EPMAGOBS
     &  ,ITERSE_SIM_EPMAGREST

      INTEGER N_SNLC_PLOT
      LOGICAL MADE_LCPLOT  ! SAVE:  T if LC plot was made

      CHARACTER VARLIST_TERSE(MXLISTNML)*(MXCHAR_CCID) ! TERSE variable names

      COMMON / CTRLCOM /
     &     SNANA_VERSION, SURVEY_NAME, SURVEY_NAME_LIST
     &    ,SUBSURVEY_NAME, SUBSURVEY_NAME_LIST
     &    ,NSURVEY_LIST, IDSURVEY, IDSUBSURVEY, IDSURVEY_LIST
     &    ,SURVEY_FILTERS
     &    ,SURVEY_FIELDNAME, SURVEY_IDFIELD, NFIELD_SURVEY
     &    ,ISJOB_SNANA, ISJOB_PSNID, ZEROPOINT_FLUXCAL
     &    ,NACCEPT_CUT, NACCEPT_CID, NACCEPT_TYPE
     &    ,NACCEPT_Z, NACCEPT_ZERR, NMJD_IDTEL
     &    ,DO_FIT, DO_FLUXERRCALC, NVAR_NEARNBR
     &    ,NSTORE_MAGCOR, NUSE_MAGCOR, SIGN_MAGCOR, FORCEMASK_FLUXCOR
     &    ,LSIM_SNANA, LSIM_MAGOBS, APPLY_HEADER_CUTMASK
     &    ,ISTAGE_SNANA
     &    ,N_VERSION, N_SNLC, N_SNLC_CUTS, N_SNLC_FIT, N_SNLC_FITCUTS
     &    ,N_SNLC_COVFIX, N_SNFILE, N_SNFILE_LAST, N_DUPLICATE_CID
     &    ,N_DUPLICATE_MJD, NSTORE_DUPLICATE_MJD
     &    ,NEPOCH_TOT, NEPOCH_CUT, NEPOCH_BADPHOT, NEPOCH_BADPHOT_SUM
     &    ,REFORMAT_SAVE_BADEPOCHS, STDOUT_UPDATE
     &    ,DOFUDGE_HOSTNOISE, DOFUDGE_NONLIN, DOFUDGE_FLUXERRMODEL
     &    ,NVAR_TERSE, VARLIST_TERSE
     &    ,ITERSE_MJD, ITERSE_FILTER, ITERSE_FIELD
     &    ,ITERSE_FLUXCAL, ITERSE_FLUXCALERR
     &    ,ITERSE_PHOTFLAG, ITERSE_PHOTPROB
     &    ,ITERSE_ZPFLUX, ITERSE_PSFSIG, ITERSE_SKYSIG, ITERSE_SKYSIG_T
     &    ,ITERSE_MAGOBS, ITERSE_XPIX, ITERSE_YPIX, ITERSE_GAIN
     &    ,ITERSE_CCDNUM, ITERSE_SIM_EPFILTREST
     &    ,ITERSE_SIM_EPMAGOBS, ITERSE_SIM_EPMAGREST
     &    ,JTIME_START, JTIME_LOOPSTART, JTIME_LOOPEND
     &    ,NCALL_SNANA_DRIVER, NPASSCUT_INCREMENT, NPASSCUT_FIT
     &    ,N_SNLC_PLOT, MADE_LCPLOT
     &    ,EXIT_ERRCODE

c logical *1 stuff

      COMMON / SNDATCOM1 / EXIST_FILT, FOUND_SURVEY
     &    ,FORMAT_TEXT, FORMAT_TERSE, FORMAT_VERBOSE
     &    ,FORMAT_FITS

c SNTABLE control variables (May 2014)
      LOGICAL
     &   USE_TABLEFILE_HBOOK  ! write table in HBOOK format
     &  ,USE_TABLEFILE_ROOT   ! write table in ROOT format
     &  ,USE_TABLEFILE_TEXT   ! write table in TEXT format
     &  ,USE_TABLEFILE        ! T if either any of the above are set
     &  ,WRTABLEFILE_IAUC     ! T-> SNID -> IAUC (for writing tables)
     &  ,WRTABLEFILE_SIMVAR   ! T-> include SIM_XXX vars for simulation
     &  ,WRTABLEFILE_ZPHOT    ! T-> include ZPHOT info for PHOTOZ fit


      INTEGER
     &   OPT_TABLE(MXTABLE)       ! option(s) for each table
     &  ,CUTMASK_SNANA_TABLE      ! select CUTFLAG_SNANA for SNANA table

      REAL*8
     &   PRESCALE_TABLE(MXTABLE)  ! table pre-scale (to reduce size)

      CHARACTER
     &   TEXTFORMAT_TABLE(MXTABLE)*8

c arrays for TABLE_FILTER_REMAP (Feb 2017)
      INTEGER NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE(MXFILT_ALL)
      CHARACTER FILTLIST_REMAP_TABLE*(MXFILT_ALL)

c Sep 23 2017: allow user codes to add private variables to SNTABLE
      INTEGER   NTABLEVAR_USER  ! e.g., set in snana_private.cra
      CHARACTER TABLEVARNAME_USER(MXVAR_PRIVATE)*40
      REAL      TABLEVALUE_USER(MXVAR_PRIVATE)

      COMMON / SNTABLECOM /
     &     USE_TABLEFILE_HBOOK, USE_TABLEFILE_ROOT
     &    ,USE_TABLEFILE_TEXT, USE_TABLEFILE
     &    ,OPT_TABLE, TEXTFORMAT_TABLE
     &    ,CUTMASK_SNANA_TABLE, WRTABLEFILE_IAUC, WRTABLEFILE_SIMVAR
     &    ,WRTABLEFILE_ZPHOT
c
     &    ,NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE
     &    ,FILTLIST_REMAP_TABLE
     &    ,NTABLEVAR_USER, TABLEVARNAME_USER
     &    ,TABLEVALUE_USER

      COMMON / SNTABLECOM8 / PRESCALE_TABLE

c May 6, 2008: define lists for epochs to ignore

      INTEGER NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
      CHARACTER
     &   EPOCH_IGNORE_CCID(MXEPOCH_IGNORE)*(MXCHAR_CCID)
     &  ,EPOCH_IGNORE_FILT(MXEPOCH_IGNORE)*4
     &  ,EPOCH_IGNORE_LASTFILE*(MXCHAR_FILENAME)

      REAL*8
     &   EPOCH_IGNORE_MJD(MXEPOCH_IGNORE)

      COMMON / EPIGNORE_COM / NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
     &  ,EPOCH_IGNORE_CCID, EPOCH_IGNORE_FILT, EPOCH_IGNORE_LASTFILE
      COMMON / EPIGNORE_COM8 / EPOCH_IGNORE_MJD


      REAL*8  DUPLICATE_MJDLIST(200)
      COMMON / DUPMJDCOM / DUPLICATE_MJDLIST


C +CDE,SNLCCOM. inserted below by ypatchy. 

c Main common block variables for light curves and
c analysis-related variables.

      INTEGER
     &   SNLC_CID               ! integer CID
     &  ,ISNLC_SNRECON_USE(MXEPOCH)  ! 1 -> epoch used in SNRECON
     &  ,ISNLC_PHOTFLAG(MXEPOCH) ! photomety flags
     &  ,ISNLC_LENCCID     ! char-len of CCID
     &  ,ISNLC_LENIAUC     ! char-len of IAUC name
     &  ,ISNLC_VERSION     ! photometry version index vs. ISN
     &  ,ISNLC_TYPE        ! type (120=confirmed Ia, etc ...)
     &  ,ISNLC_ABSO_INDEX  ! absolute index (row or file number)
     &  ,ISNLC_IFILE       ! ifile index, text format only
     &  ,ISNLC_IFILT_OBS(MXEPOCH) ! filt-index for each epoch/SN
     &  ,ISNLC_NEWMJD_HEAD   ! Number of NEWMJDs in header
     &  ,ISNLC_NEWMJD_FOUND  ! Number of NEWMJDs found in file
     &  ,ISNLC_NEWMJD_STORE  ! Number of NEWMJDs stored
     &  ,ISNLC_NEWMJD_CUTS   ! Number of NEWMJDs after cuts
     &  ,ISNLC_EPOCH_RANGE_NEWMJD(2,MXEPOCH)
     &  ,ISNLC_NFILT_NEWMJD(MXEPOCH)   ! # observed filters per NEWMJD
     &  ,ISNLC_NFILT_SNRMAX      ! number of filters passing SNR cut
     &  ,ISNLC_NFILT_SNRMAX2     ! idem for 2nd SNRMAX cut
     &  ,ISNLC_NFILT_TRESTMIN     ! Nfilt passing TRESTMIN cut
     &  ,ISNLC_NFILT_TRESTMAX     ! idem for TRESTMAX
     &  ,ISNLC_NFILT_TREST2       ! Nfilt passing CUTWIN_TREST2 cut
     &  ,ISNLC_NFILT_THRESH(MXEPOCH)
     &  ,ISNLC_NEPOCH_FOUND  ! actual number of epochs found
     &  ,ISNLC_NEPOCH_STORE  ! number of epochs stored in memory
     &  ,ISNLC_NEPOCH_USE       ! used after PHOTMASK cuts
     &  ,ISNLC_NEPOCH_PHOTPROB  ! NEPOCH with PHOTPROB >= 0
     &  ,ISNLC_NEPOCH_DETECT    ! NEPOCH with detection
     &  ,ISNLC_NEPOCH_FILT(MXFILT_OBS)  ! NEPOCH vs. filter
     &  ,ISNLC_NEPOCH_PRESN(MXFILT_OBS) ! pre-explosion epochs
     &  ,ISNLC_NMJD_INCLUDE             ! NOBS in CUTWIN_MJD_INCLUDE window
c
     &  ,ISNLC_FAKE             ! => real data, else it's a fake
     &  ,ISNLC_CCDNUM(MXEPOCH)  ! read from header (May 2017)
     &  ,ISNLC_IDTEL(MXEPOCH)   ! integer telescope id
     &  ,ISNLC_IDFIELD(MXEPOCH) ! integer field id
     &  ,ISNLC_NFIELD_OVP       ! number of fields (>=2 for overlap)
     &  ,ISNLC_CUTFLAG_REQEP        ! idem
     &  ,ISNLC_CUTFLAG_PRIVATE      ! idem for private var cuts
     &  ,ISNLC_CUTFLAG_USRCUTS      ! idem for USRCUTS
     &  ,ISNLC_CUTFLAG_SIMVAR       ! idem for SIMVAR cuts
     &  ,ISNLC_WRMASK_FLUXCOR_SNANA   ! write if fudges applied to data
     &  ,ISNLC_RDMASK_FLUXCOR_SNANA   ! read if SNANA fudges applied to data

      REAL*8
     &   SNLC8_RA        ! RA   vs. SN
     &  ,SNLC8_DEC       ! DEC  vs. SN
     &  ,SNLC8_MJD(MXEPOCH)  ! MJD for each epoch
     &  ,SNLC8_MJDMIN        !   min MJD among all measurements
     &  ,SNLC8_MJDMAX        !   max MJD ...

      REAL
     &   SNLC_ZHELIO        ! redshift used for final analysis
     &  ,SNLC_ZHELIO_ERR    ! error on above
     &  ,SNLC_ZCMB          ! redshift used for final analysis
     &  ,SNLC_ZCMB_ERR      ! error on above
     &  ,SNLC_ZHD           ! redshift used for Hubble diagram
     &  ,SNLC_ZHD_ERR       ! error on above
     &  ,SNLC_ZSN           ! redshift of SN only (ignoring host-z)
     &  ,SNLC_ZSN_ERR       ! error on above
     &  ,SNLC_REDSHIFT      ! redshift used for light curve fit
     &  ,SNLC_REDSHIFT_ERR  ! error on above
     &  ,SNLC_VPEC          ! pec. velocity
     &  ,SNLC_VPEC_ERR      ! error on above
     &  ,SNLC_ZPEC
     &  ,SNLC_ZPEC_ERR
     &  ,SNLC_Trestmin  ! earliest epoch, rest frame days since peak
     &  ,SNLC_Trestmax  ! latest   epoch, rest frame days since peak
     &  ,SNLC_TrestRange
     &  ,SNLC_Tobsmin   ! earliest epoch, obs frame days since peak
     &  ,SNLC_Tobsmax   ! latest   epoch, obs frame days since peak
     &  ,SNLC_TGAPMAX     ! max gap within TREST-range
     &  ,SNLC_T0GAPMAX    ! max gap near peak
     &  ,SNLC_SNRMAX_FILT(0:MXFILT_OBS)   ! max S/N per filter/SN
     &  ,SNLC_SNRMAX_SORT(MXFILT_OBS)     ! 1st, 2nd ... SNRMAX by filt
     &  ,SNLC_FLUXCALMAX(MXFILT_OBS)    ! max flux per filter/SN
     &  ,SNLC_SNANACALC_PEAKMJD     ! SNANA-estimate of PEAKKMJD
     &  ,SNLC_SNANACALC_PEAKMJD_FITPAR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_SNANACALC_PEAKMJD_FITERR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_PHOTPROB(MXEPOCH)   ! generic 'fit probability' per epoch
     &  ,SNLC_PHOTPROB_MIN        ! min photprob for PHOTPROB>0
     &  ,SNLC_AIRMASS(MXEPOCH)
     &  ,SNLC_CLOUDCAM_AVG(MXEPOCH)
     &  ,SNLC_CLOUDCAM_SIG(MXEPOCH)
     &  ,SNLC_MOONDIST(MXEPOCH)    ! in deg.
     &  ,SNLC_MOONPHASE(MXEPOCH)   ! 0 - 1
     &  ,SNLC_TREST(MXEPOCH)          !MJD-SearhMJDPEAK / (1+z)
     &  ,SNLC_GAIN(MXEPOCH)      ! e/AUD
     &  ,SNLC_RDNOISE(MXEPOCH)   ! read noise per pix, e-
     &  ,SNLC_PIXSIZE            ! pixel size
     &  ,SNLC_NXPIX              ! total number of X-pixels (Aug 7 2014)
     &  ,SNLC_NYPIX              ! total number of Y-pixels
     &  ,SNLC_XPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_YPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_AREAFRAC(MXEPOCH)  ! area-frac contained by XPIX,YPIX
     &  ,SNLC_MWEBV              ! Milky Way Galactic E(B-V)
     &  ,SNLC_MWEBV_ERR          ! error on above
     &  ,SNLC_SKY_AVG(MXEPOCH)   ! avg sky/pix in field, ADU
     &  ,SNLC_SKYSIG(MXEPOCH)    ! sigma on above
     &  ,SNLC_SKYSIG_T(MXEPOCH)  ! sigma on template run
     &  ,SNLC_PSF_SIG1(MXEPOCH)  ! sigma, pixels
     &  ,SNLC_PSF_SIG2(MXEPOCH)
     &  ,SNLC_PSF_RATIO(MXEPOCH)
     &  ,SNLC_AREA_NOISE(MXEPOCH)  ! Noise-equivalent area (pixels)
     &  ,SNLC_FLUX(MXEPOCH)
     &  ,SNLC_FLUX_ERRTOT(MXEPOCH)
     &  ,SNLC_FLUX_NSIG(MXEPOCH)   ! significance
     &  ,SNLC_FLUXCAL_OFF(MXFILT_OBS)    ! add SN light from template
     &  ,SNLC_FLUXCAL_ERRCALC(MXEPOCH)   ! calc a-la simulation
     &  ,SNLC_FLUXCAL_HOSTERRCALC(MXEPOCH)   ! idem for host error
     &  ,SNLC_FLUXCAL(MXEPOCH)
     &  ,SNLC_FLUXCAL_ERRTOT(MXEPOCH)
     &  ,SNLC_MAG(MXEPOCH)
     &  ,SNLC_MAG_ERRPLUS(MXEPOCH)
     &  ,SNLC_MAG_ERRMINUS(MXEPOCH)
     &  ,SNLC_ZEROPT(MXEPOCH)
     &  ,SNLC_ZEROPT_ERR(MXEPOCH)
     &  ,SNLC_ZEROPT_forCUT(MXEPOCH)  ! depends on CUTWIN_ZPADU or CUTWIN_ZPNPE
     &  ,SNLC_DLMAG                   ! 5*log10(10pc/DL)
     &  ,SNLC_SKYFLUXCAL(MXEPOCH)     ! calculated sky fluxcal/pixel
     &  ,SNLC_MWXT_MAG(MXFILT_OBS)        ! mag stellar extinct
     &  ,SNLC_MWXT_FLUXFRAC(MXFILT_OBS)   ! same for flux (< 1)
     &  ,SNLC_MWXT_MAGERR(MXFILT_OBS)     ! Galactic mag err per filter
     &  ,SNLC_SEARCH_PEAKMJD            ! external PEAKMJD
     &  ,SNLC_DTOBS(MXEPOCH)            ! time since last obs
     &  ,SNLC_DTOBS_SAMEFILT(MXEPOCH)   ! idem, but same filter
     &  ,SNLC_TLIVE_DETECT     ! MJD(last detection) - MJD(1st detection)

      INTEGER
     &   NSEASON_TOT     ! total number of seasons
     &  ,NSEASON_ACTIVE  ! total number of active seasons

      REAL*4
     &   MULTISEASON_CHI2RED(MXSEASON) ! borrow another 'MX' param
     &  ,MULTISEASON_AVGFLUX(MXSEASON)
     &  ,MULTISEASON_MJDMIN(MXSEASON)
     &  ,MULTISEASON_MJDMAX(MXSEASON)

      CHARACTER
     &   SNLC_CCID*(MXCHAR_CCID)      ! char CID
     &  ,SNLC_IAUC*(MXCHAR_CCID)
     &  ,SNLC_CASTCID*8             ! = 'INT' or 'CHAR' to indicate cast
     &  ,SNLC_FIELD(MXEPOCH)*(MXCHAR_FIELDNAME)    ! char field name
     &  ,SNLC_FIELD_OVPLIST(MXFIELD_OVP)*(MXCHAR_FIELDNAME)  ! sparse list of overlap fields
     &  ,SNLC_FIELDLIST*(MXCHAR_FIELDNAME) ! e.g., '82S', '82N+82S'

      COMMON / SNLC8COM /
     &    SNLC8_RA, SNLC8_DEC
     &  , SNLC8_MJD, SNLC8_MJDMIN, SNLC8_MJDMAX

      COMMON / SNLCCOM /
     &     SNLC_CID, SNLC_CCID,  SNLC_CASTCID, SNLC_IAUC
     &    ,SNLC_FIELD, SNLC_FIELD_OVPLIST, SNLC_FIELDLIST
     &    ,SNLC_SNANACALC_PEAKMJD
     &    ,SNLC_SNANACALC_PEAKMJD_FITPAR
     &    ,SNLC_SNANACALC_PEAKMJD_FITERR
     &    ,SNLC_SEARCH_PEAKMJD
     &    ,SNLC_ZHELIO, SNLC_ZHELIO_ERR
     &    ,SNLC_ZCMB,   SNLC_ZCMB_ERR, SNLC_ZHD, SNLC_ZHD_ERR
     &    ,SNLC_ZSN, SNLC_ZSN_ERR
     &    ,SNLC_REDSHIFT, SNLC_REDSHIFT_ERR
     &    ,SNLC_VPEC, SNLC_VPEC_ERR, SNLC_ZPEC, SNLC_ZPEC_ERR
     &    ,SNLC_TRESTMIN, SNLC_TRESTMAX, SNLC_TrestRange, SNLC_TREST
     &    ,SNLC_TGAPMAX,  SNLC_T0GAPMAX
     &    ,SNLC_Tobsmin, SNLC_Tobsmax
     &    ,SNLC_SNRMAX_FILT, SNLC_SNRMAX_SORT
     &    ,SNLC_FLUXCALMAX, SNLC_PHOTPROB, SNLC_PHOTPROB_MIN
     &    ,SNLC_CLOUDCAM_AVG, SNLC_CLOUDCAM_SIG
     &    ,SNLC_MOONDIST, SNLC_MOONPHASE, SNLC_AIRMASS
     &    ,SNLC_GAIN, SNLC_RDNOISE
     &    ,SNLC_MWXT_MAG, SNLC_MWXT_FLUXFRAC, SNLC_MWXT_MAGERR
     &    ,SNLC_MWEBV, SNLC_MWEBV_ERR
     &    ,SNLC_XPIX, SNLC_YPIX, SNLC_AREAFRAC, SNLC_PIXSIZE
     &    ,SNLC_NXPIX, SNLC_NYPIX
     &    ,SNLC_SKY_AVG,  SNLC_SKYSIG,  SNLC_SKYSIG_T
     &    ,SNLC_PSF_SIG1, SNLC_PSF_SIG2, SNLC_PSF_RATIO
     &    ,SNLC_AREA_NOISE
c
     &    ,SNLC_FLUX, SNLC_FLUX_ERRTOT
     &    ,SNLC_FLUX_NSIG, SNLC_FLUXCAL_OFF
     &    ,SNLC_FLUXCAL_ERRCALC, SNLC_FLUXCAL_HOSTERRCALC
     &    ,SNLC_FLUXCAL, SNLC_FLUXCAL_ERRTOT
     &    ,SNLC_MAG, SNLC_MAG_ERRPLUS, SNLC_MAG_ERRMINUS
     &    ,SNLC_ZEROPT, SNLC_ZEROPT_ERR, SNLC_ZEROPT_forCUT
     &    ,SNLC_SKYFLUXCAL, SNLC_DLMAG
     &    ,SNLC_DTOBS, SNLC_DTOBS_SAMEFILT

      COMMON / ISNCOM /
     &     ISNLC_VERSION, ISNLC_SNRECON_USE, ISNLC_PHOTFLAG
     &    ,ISNLC_NFILT_THRESH, ISNLC_TYPE
     &    ,ISNLC_ABSO_INDEX, ISNLC_IFILE,ISNLC_LENCCID,ISNLC_LENIAUC
     &    ,ISNLC_NFILT_NEWMJD, ISNLC_IFILT_OBS
     &    ,ISNLC_NEWMJD_HEAD, ISNLC_NEWMJD_FOUND, ISNLC_NEWMJD_STORE
     &    ,ISNLC_NEWMJD_CUTS, ISNLC_EPOCH_RANGE_NEWMJD
     &    ,ISNLC_NEPOCH_FOUND, ISNLC_NEPOCH_STORE
     &    ,ISNLC_NEPOCH_USE, ISNLC_NEPOCH_PHOTPROB
     &    ,ISNLC_NEPOCH_FILT,  ISNLC_NEPOCH_PRESN, ISNLC_NMJD_INCLUDE
     &    ,ISNLC_NEPOCH_DETECT, SNLC_TLIVE_DETECT
     &    ,ISNLC_FAKE, ISNLC_CCDNUM
     &    ,ISNLC_NFILT_SNRMAX, ISNLC_NFILT_SNRMAX2
     &    ,ISNLC_NFILT_TRESTMIN, ISNLC_NFILT_TRESTMAX
     &    ,ISNLC_NFILT_TREST2
     &    ,ISNLC_IDTEL, ISNLC_IDFIELD, ISNLC_NFIELD_OVP
     &    ,ISNLC_CUTFLAG_REQEP
     &    ,ISNLC_CUTFLAG_PRIVATE, ISNLC_CUTFLAG_USRCUTS
     &    ,ISNLC_CUTFLAG_SIMVAR
     &    ,ISNLC_WRMASK_FLUXCOR_SNANA, ISNLC_RDMASK_FLUXCOR_SNANA

c - - - - -  ADDCOL stuff - - - - - - -
c ADDCOL arrays are loaded with original filter indices,
c or with REMAPed filter indices. Allows using unique
c pointer for calls to SNTABLE_ADDCOL.

      CHARACTER
     &    ADDCOL_FILTERS*(MXFILT_ALL) ! SURVEY_FILTERS

      REAL
     &    ADDCOL_SNHOST_MAGOBS(MXFILT_ALL,MXSNHOST)
     &   ,ADDCOL_SNHOST_SBFLUXCAL(MXFILT_ALL)
     &   ,ADDCOL_FLUXCALMAX(MXFILT_ALL)
     &   ,ADDCOL_CHI2_FITPKMJD(MXFILT_ALL)
     &   ,ADDCOL_SNRMAX(MXFILT_ALL)
     &   ,ADDCOL_XTMW(MXFILT_ALL)

      INTEGER
     &    ADDCOL_NDOF_FITPKMJD(MXFILT_ALL)

      COMMON / ADDCOL_FILTERCOM /
     &   ADDCOL_FILTERS
     &  ,ADDCOL_SNHOST_MAGOBS, ADDCOL_SNHOST_SBFLUXCAL
     &  ,ADDCOL_FLUXCALMAX, ADDCOL_CHI2_FITPKMJD, ADDCOL_NDOF_FITPKMJD
     &  ,ADDCOL_SNRMAX, ADDCOL_XTMW


c - - - -- MULTI-SEASON STUFF - - - -
      COMMON / MULTISEASONCOM / NSEASON_TOT, NSEASON_ACTIVE
     &  ,MULTISEASON_CHI2RED, MULTISEASON_AVGFLUX
     &  ,MULTISEASON_MJDMIN,  MULTISEASON_MJDMAX


C +CDE,SNCUTS. inserted below by ypatchy. 

c Define cut-window for each variables
c and cut-mask for each epoch.
c User must fill cutwin_XXX arrays in main routine.

      INTEGER
     &   CUTBIT_CID, CUTBIT_SNTYPE
     &  ,CUTBIT_REDSHIFT, CUTBIT_REDSHIFT_ERR
     &  ,CUTBIT_RA,  CUTBIT_DEC
     &  ,CUTBIT_HOSTSEP
     &  ,CUTBIT_Trestmin, CUTBIT_Trestmax, CUTBIT_TrestRange
     &  ,CUTBIT_TREST2
     &  ,CUTBIT_Tgapmax,  CUTBIT_T0gapmax
     &  ,CUTBIT_Tobsmin,  CUTBIT_Tobsmax
     &  ,CUTBIT_SEARCH           ! sim only
     &  ,CUTBIT_TEMPLATE_INDEX   ! sim only (Dec 2018)
     &  ,CUTBIT_PEAKMJD
     &  ,CUTBIT_NMJD_INCLUDE    ! Aug 11 2015
     &  ,CUTBIT_NEPOCH
     &  ,CUTBIT_SNRMAX         ! global SNRMAX cut; any filter(s)
     &  ,CUTBIT_SNRMAX2        ! 2nd global SNRMAX cut; any filter(s)
     &  ,CUTBIT_NFIELD         ! cut on number of fields used by SN
     &  ,CUTBIT_MWEBV          ! cut on Galactic extinct (May 8 2012)
     &  ,CUTBIT_NSEASON_ACTIVE ! Nseason with activity
     &  ,CUTBIT_REQEP          ! required epochs (to select short events)
     &  ,CUTBIT_PRIVATE        ! all private-var cuts together (Nov 3 2014)
     &  ,CUTBIT_SIMVAR         ! SIMVAR cuts
     &  ,CUTBIT_USRCUTS        ! user cuts from private code (Aug 2015)
     &  ,CUTBIT_OFFSET_SBFLUX  ! host surface brightness
     &  ,CUTBIT_OFFSET_SNRMAX  ! SNRMAX cuts in each filter
     &  ,CUTBIT_MJD_MARKER     ! things above do NOT depend on epoch
     &  ,CUTBIT_PSF,  CUTBIT_ZP
     &  ,CUTBIT_MOONDIST,  CUTBIT_PHOTPROB, CUTBIT_AIRMASS
     &  ,CUTBIT_TREST, CUTBIT_MJD
     &  ,CUTBIT_IDTEL          ! valid telescope
     &  ,CUTBIT_NBAND_THRESH
     &  ,CUTBIT_NFILT_SNRMAX
     &  ,CUTBIT_NFILT_SNRMAX2
     &  ,CUTBIT_NFILT_TRESTMIN
     &  ,CUTBIT_NFILT_TRESTMAX
     &  ,CUTBIT_NFILT_TREST2
     &  ,CUTBIT_OFFSET
     &  ,NCUTBIT
     &  ,NCUTBIT_SNLC

      PARAMETER (
     &   CUTBIT_CID            = 1
     &  ,CUTBIT_SNTYPE         = 2
     &  ,CUTBIT_REDSHIFT       = 3
     &  ,CUTBIT_REDSHIFT_ERR   = 4
     &  ,CUTBIT_RA             = 5
     &  ,CUTBIT_DEC            = 6
     &  ,CUTBIT_HOSTSEP        = 7
     &  ,CUTBIT_TrestMIN       = 8
     &  ,CUTBIT_TrestMAX       = 9
     &  ,CUTBIT_TrestRange     = 10  ! Dec 2017 (TrestMax-TrestMin)
     &  ,CUTBIT_Trest2         = 11  ! Oct 2012
     &  ,CUTBIT_Tgapmax        = 12
     &  ,CUTBIT_T0gapMAX       = 13
     &  ,CUTBIT_TobsMIN        = 14
     &  ,CUTBIT_TobsMAX        = 15
     &  ,CUTBIT_PEAKMJD        = 16
     &  ,CUTBIT_NMJD_INCLUDE   = 17
     &  ,CUTBIT_NEPOCH         = 18  ! total # of measurements
     &  ,CUTBIT_SEARCH         = 19  ! found by search (SIM only)
     &  ,CUTBIT_NFILT_SNRMAX   = 20
     &  ,CUTBIT_NFILT_SNRMAX2  = 21
     &  ,CUTBIT_NFILT_TRESTMIN = 22
     &  ,CUTBIT_NFILT_TRESTMAX = 23
     &  ,CUTBIT_NFILT_TREST2   = 24
     &  ,CUTBIT_SNRMAX         = 25
     &  ,CUTBIT_SNRMAX2        = 26
     &  ,CUTBIT_NFIELD         = 27
     &  ,CUTBIT_MWEBV          = 28  ! May 2012
     &  ,CUTBIT_NSEASON_ACTIVE = 29  ! May 2019
     &  ,CUTBIT_REQEP          = 30  ! required epochs (9/17/2017)
     &  ,CUTBIT_PRIVATE        = 31  ! all private-var cuts
     &  ,CUTBIT_SIMVAR        = 32  !
     &  ,CUTBIT_USRCUTS        = 33  ! PASS_USRCUTS=T
     &  ,CUTBIT_OFFSET_SBFLUX  = 34  ! Mar 2012
     &  ,CUTBIT_OFFSET_SNRMAX  = 34 + MXFILT_SNRMAX
     &  ,CUTBIT_MJD_MARKER     = CUTBIT_OFFSET_SNRMAX + MXFILT_SNRMAX
     &  ,NCUTBIT_SNLC          = CUTBIT_MJD_MARKER
c SN-dependent cuts above
c MJD-dependent cuts below
     &  ,CUTBIT_PSF          = CUTBIT_MJD_MARKER + 1  ! PSF cut
     &  ,CUTBIT_ZP           = CUTBIT_MJD_MARKER + 2  ! ZP  cut (Feb 2018)
     &  ,CUTBIT_MOONDIST     = CUTBIT_MJD_MARKER + 3  ! degrees
     &  ,CUTBIT_PHOTPROB     = CUTBIT_MJD_MARKER + 4  ! Oct 23 2018
     &  ,CUTBIT_AIRMASS      = CUTBIT_MJD_MARKER + 5  ! 1=overhead
     &  ,CUTBIT_NBAND_THRESH = CUTBIT_MJD_MARKER + 6
     &  ,CUTBIT_Trest        = CUTBIT_MJD_MARKER + 7
     &  ,CUTBIT_MJD          = CUTBIT_MJD_MARKER + 8
     &  ,CUTBIT_IDTEL        = CUTBIT_MJD_MARKER + 9
c
     &  ,CUTBIT_OFFSET       = CUTBIT_IDTEL
     &  ,NCUTBIT             = CUTBIT_OFFSET
     &      )

c Define cut-masks with 64-BIT integers

      INTEGER*8
     &   CUTMASK8_SN             ! cutmask for each SN
     &  ,CUTMASK8_MJD(MXEPOCH)   ! cutmask vs. MJD, isn
     &  ,CUTMASK8_SN_ALL         ! all bits for SN cuts
     &  ,CUTMASK8_MJD_ALL        ! all bits for MJD cuts

      BYTE MJDMASK(MXEPOCH)  ! logical for each epoch,SN

      LOGICAL
     &   LSNCUTS       ! T=> passes cuts 1 to CUTBIT_MJD_MARKER
     &  ,PASS_USRCUTS  ! T=> user cuts with private code (Oct 2014)
     &  ,PASS_PRIVCUTS ! T=> pass cuts on private var in data header
     &  ,PASS_SIMCUTS  ! T=> pass cuts on SIMVAR

      INTEGER
     &   NCCID_LIST      ! size of SNCCID_LIST
     &  ,NCID_LIST       ! size of SNCID_LIST
     &  ,NCCID_IGNORE    ! size of SNCCID_IGNORE
     &  ,NCID_IGNORE     ! size of NCID_IGNORE
     &  ,CUTFLAG_SNANA   ! bits 0,1 -> LSNCUTS,LSNFITOK=T (for ntuple)
     &  ,ERRFLAG_FIT     ! error flag from fit (0=OK)

c stuff for hbook
      CHARACTER CUTVAR_NAME(NCUTBIT)*28

c define namelist cuts

      INTEGER
     &   SNTYPE_LIST(MXLISTNML)    ! (I) user list of types to select
     &  ,SNTYPE_IGNORE(MXLISTNML)  ! (I) types to ignore
     &  ,SNCID_LIST(MXLISTNML)     ! (I) user-list of integer CIDs
     &  ,CCDNUM_LIST(MXLISTNML)    ! (I) user-list of CCDNUMs (Dec 2017)
     &  ,SIM_TEMPLATE_INDEX_LIST(MXLISTNML)
     &  ,SNCID_IGNORE(MXLISTNML)   ! (I) list of CIDs to ignore
     &  ,PHOTFLAG_MSKREJ(5)        ! (I) PHOTFLAG mask-bits to reject
     &                             !     1,2 => logical OR,AND
     &  ,IDTEL_LIST(MXLISTNML)     ! computed from SNTEL_LIST
     &  ,IDFIELD_LIST(MXLISTNML)   !
     &  ,NIDTEL_LIST               ! size of list
     &  ,NIDFIELD_LIST
     &  ,NSNTYPE_LIST           ! size of SNTYPE_LIST
     &  ,NCCDNUM_LIST           ! size of CCDNUM list
     &  ,NSIM_TEMPLATE_INDEX_LIST
     &  ,NSNTYPE_IGNORE
     &  ,NFILT_SNRMAX           ! number of CUTWIN_SNRMAX_FILT cuts
     &  ,IFILT_SNRMAX(MXFILT_SNRMAX) ! store filt index for SNRMAX cuts
     &  ,NFILT_HOST_SBFLUX      ! Nfilt to cut on HOST_SBFLUX
     &  ,IFILT_HOST_SBFLUX(MXFILT_SNRMAX)

      LOGICAL DOALL_SNTEL, DOALL_SNFIELD

      character
     &    SNTEL_LIST(MXLISTNML)*20    ! (I) list of telescopes (or ALL)
     &   ,SNFIELD_LIST(MXLISTNML)*60  ! (I) list of fields to use (or 'ALL')
     &   ,SNCCID_LIST(MXLISTNML)*(MXCHAR_CCID)   ! (I) list of CIDs to process
     &   ,SNCCID_LIST_ALL(MXCCID_LIST)*(MXCHAR_CCID) ! combined list to process
     &   ,SNCCID_IGNORE(MXLISTNML)*(MXCHAR_CCID) ! (I) list of CIDs to ignore
     &   ,SNCID_IGNORE_FILE*(MXCHAR_FILENAME) ! (I) file with list of CIDs to ignore

      REAL
     &   snlc_cutvar(NCUTBIT_SNLC)
     &  ,snep_cutvar(CUTBIT_MJD_MARKER:NCUTBIT,MXEPOCH)
     &  ,cutwin_var(2,NCUTBIT)
     &  ,cutwin_cid(2)           ! candidate id
     &  ,cutwin_redshift(2)
     &  ,cutwin_redshift_err(2)
     &  ,cutwin_ra(2)
     &  ,cutwin_dec(2)
     &  ,cutwin_hostsep(2)      ! cut on host-SN sep, arcsec
     &  ,cutwin_sbflux_filt(2,MXFILT_SNRMAX)
     &  ,cutwin_Nepoch(2)       !
     &  ,cutwin_snrmax_filt(2,MXFILT_SNRMAX)   ! filled from SNCUT_SNRMAX
     &  ,cutwin_snrmax(2)       ! global SNRMAX cut
     &  ,cutwin_snrmax2(2)       ! 2nd global SNRMAX cut
     &  ,cutwin_nfield(2)        ! number of fields (usually 1)
     &  ,cutwin_mwebv(2)         ! Galactic extinction
     &  ,cutwin_nseason_active(2)
     &  ,cutwin_cutflag_reqep(2)       ! for internal use only
     &  ,cutwin_cutflag_private(2)     ! for internal use only
     &  ,cutwin_cutflag_simvar(2)      ! for internal use only
     &  ,cutwin_cutflag_usrcuts(2)     ! for internal use only
     &  ,cutwin_nfilt_snrmax(2)  ! Nfilt passing global SNRMAX cut
     &  ,cutwin_nfilt_snrmax2(2) ! Nfilt passing 2nd-best SNRMAX
     &  ,cutwin_nfilt_trestmin(2)! Nfilt passing Trestmin cut
     &  ,cutwin_nfilt_trestmax(2)! Nfilt passing Trestmax cut
     &  ,cutwin_nfilt_trest2(2)  ! Nfilt passinng Trest2 cut
     &  ,cutwin_nband_thresh(2)  ! number of ugriz bands above $band_THRESH
     &  ,cutwin_Trestmin(2)      ! window for earliest epoch, rest frame days
     &  ,cutwin_Trestmax(2)      ! window for latest epoch, rest frame day
     &  ,cutwin_TrestRange(2)    ! window on Trestmax - Trestmin
     &  ,cutwin_Trest2(2)        ! window for CUTWIN_NFILT_TREST2
     &  ,cutwin_Tgapmax(2)       ! max Trest-gap within cutwin_TREST(2)
     &  ,cutwin_T0gapmax(2)      ! max Trest-gap near peak
     &  ,cutwin_Tobsmin(2)
     &  ,cutwin_Tobsmax(2)
     &  ,cutwin_peakmjd(2)       ! cut on search peakMJD (for rates)
c
c  Below are Epoch-dependent cuts

     &  ,cutwin_psf(2)           ! PSF cut, FWHM, ARCSEC
     &  ,cutwin_zp(2)            ! ZEROPT cut, ADU or NPE
     &  ,cutwin_zpADU(2)         ! ZEROPT cut, ADU
     &  ,cutwin_zpNPE(2)         ! ZEROPT cut, NPE
     &  ,cutwin_moondist(2)      ! dist betwee SN and moon, deg
     &  ,cutwin_photprob(2)
     &  ,cutwin_Trest(2)         ! window for all epochs, rest frame days
     &  ,cutwin_MJD(2)           ! MJD window
     &  ,CUTWIN_MJD_EXCLUDE(2)          ! MJD window to exclude
     &  ,CUTWIN_MJD_INCLUDE(2)          ! require  MJD in this window
     &  ,CUTWIN_NMJD_INCLUDE(2)         ! for internal use only
     &  ,cutwin_SEARCHEFF_MASK(2) ! for SIM only
     &  ,cutwin_snrmin_filt(2,MXFILT_OBS) ! filled from EPCUT_SNRMIN
     &  ,cutwin_restlam(2)                ! cut on <LAMREST>, no cutBit
     &  ,cutwin_lamrest(2)                ! same as above
     &  ,cutwin_lamobs(2)                 ! cut on <LAMOBS>, no cutBit
c
c define character strings to specify cuts and filters
c 'SNCUT' specifies cut on each SN
c 'EPCUT' specifies selection cut on each epoch

      CHARACTER*(MXCHAR_CUTNAME)
     &  SNCUT_SNRMAX  ! max SNR required in each passband
                      ! example: 'u 10.  g 5.0  r 5.0  i 5.0  z -10.'
c
     & ,EPCUT_SNRMIN  ! min SNR accepted for each epoch/filter
                      ! example: 'u 20000.  g -4.  r -4.  i -4.  z 20000.'
c
     & ,EPCUT_SKYSIG  ! max SKY noise accepted each epoch/filter
                      ! example: 'u 20.  g 50.  r 80.  i 120.  z 200.'
c
     & ,EPCUT_PSFSIG  ! max PSF (arcsec) accepted each epoch/filter
                      ! example: 'u 1.8  g 1.5  r 1.6  i 1.8  z 2.0'
c
     & ,SNCUT_HOST_SBFLUX ! max allows surface brightness,
                          ! example 'r 1000.' -> SBFLUX < 1000

c ----------
c systematic tests for calibration:
c 'U 01 B -0.01' => shift U & B mags of primary

     &  ,MAGOBS_SHIFT_PRIMARY  ! shift primary mags  (for syst test)
     &  ,MAGOBS_SHIFT_ZP       ! shift zero points (e.g.,AB off for SDSS)
     &  ,MAGREST_SHIFT_PRIMARY ! idem for rest-frame mags
     &  ,MAGREST_SHIFT_ZP      ! idem for rest-frame filters
     &  ,FILTER_LAMSHIFT       ! e.g., 'r -2.4  i 6.2'  ! in Angstroms
c -------
c Fluxcal fudges for systematic tests (fudge photometry offsets and errors)
c Note that error is added in quadrature; ERROR<0 is subtracted in quadrature.
     &  ,FUDGE_FLUXCAL_OFFSET
     &  ,FUDGE_FLUXCAL_ERROR   ! fudge net FLUXCAL error in each band
     &  ,FUDGE_FLUXCAL_ERRPIX  ! per-pixel error --> FLUXCAL error propto PSF
     &  ,FUDGE_MAG_ERROR       ! Oct 2013: add stat error per band
     &  ,FUDGE_MAG_COVERR      ! Feb 2015: add covariant error per band

      REAL
     &   MAGOBS_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGOBS_SHIFT_ZP_USER(MXFILT_ALL)      ! user ZP shift
     &  ,MAGOBS_SHIFT_ZP_FILT(MXFILT_ALL)      ! user ZP shift + system ZP
     &  ,FLUXSCALE_ZP_FILT(MXFILT_ALL)         ! corresponding flux scale
     &  ,MAGREST_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGREST_SHIFT_ZP_USER(MXFILT_ALL)      ! shift zero points
     &  ,MAGREST_SHIFT_ZP_FILT(MXFILT_ALL)      ! shift zero points
     &  ,FILTER_LAMSHIFT_FILT(MXFILT_ALL)       ! shift filter trans
c
     &  ,MAGOBS_SHIFT_PRIMARY_PARAMS(3)       ! poly-fun of lambda;
     &  ,MAGOBS_SHIFT_ZP_PARAMS(3)            ! A0 + A1*LAM + A2*LAM^2
     &  ,FUDGE_FLUXCAL_OFFSET_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERRPIX_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_COVERR_FILT(MXFILT_ALL)
     &  ,MWEBV_SCALE        ! scale MW extinc for syst test
     &  ,MWEBV_SHIFT        ! shift MW extinc

      REAL
     &   RV_MWCOLORLAW   ! (I) RV for Galactic extinction

      INTEGER
     &   OPT_MWCOLORLAW    ! (I) MW color law opt (89, 94, 99 ...)
     &  ,OPT_MWEBV         ! (I) option to modify SFD maps

      LOGICAL
     &   USE_MWCOR          ! (I) T=> correct data flux for MW extinc;
                            !     F=> leave data, correct fit model for MW.
     &  ,DOBUG_LAMRANGE     ! (I) implement old lambda-range bug

      REAL
     &   REDSHIFT_FINAL_SHIFT    ! artificial shift
     &  ,FLUXERRCALC_ZPTERR

      COMMON / SNCUTCOM / CUTMASK8_SN, CUTMASK8_MJD
     &   ,CUTMASK8_SN_ALL, CUTMASK8_MJD_ALL
     &   ,cutwin_var, cutvar_name
     &   ,snlc_cutvar, snep_cutvar
     &   ,cutwin_snrmin_filt, LSNCUTS
     &   ,PASS_USRCUTS, PASS_PRIVCUTS, PASS_SIMCUTS
     &   ,CUTFLAG_SNANA , ERRFLAG_FIT
     &   ,NCCID_LIST, NCCID_IGNORE, NCID_LIST, NCID_IGNORE
     &   ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST, SNCCID_LIST_ALL
     &   ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE
     &   ,SNTYPE_IGNORE, CCDNUM_LIST, SIM_TEMPLATE_INDEX_LIST
     &   ,SNTEL_LIST,   IDTEL_LIST,   NIDTEL_LIST,   DOALL_SNTEL
     &   ,SNFIELD_LIST, IDFIELD_LIST, NIDFIELD_LIST, DOALL_SNFIELD
     &   ,NSNTYPE_LIST, NCCDNUM_LIST, NSIM_TEMPLATE_INDEX_LIST
     &   ,NSNTYPE_IGNORE, NFILT_SNRMAX, IFILT_SNRMAX
     &   ,NFILT_HOST_SBFLUX, IFILT_HOST_SBFLUX
     &   ,PHOTFLAG_MSKREJ
     &   ,USE_MWCOR, DOBUG_LAMRANGE
     &   ,MAGOBS_SHIFT_PRIMARY, MAGOBS_SHIFT_PRIMARY_FILT
     &   ,MAGOBS_SHIFT_ZP, MAGOBS_SHIFT_ZP_FILT, FLUXSCALE_ZP_FILT
     &   ,MAGOBS_SHIFT_ZP_USER
     &   ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS
     &   ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_PRIMARY_FILT
     &   ,MAGREST_SHIFT_ZP, MAGREST_SHIFT_ZP_FILT
     &   ,MAGREST_SHIFT_ZP_USER
     &   ,FILTER_LAMSHIFT, FILTER_LAMSHIFT_FILT
     &   ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX
     &   ,FUDGE_FLUXCAL_OFFSET_FILT
     &   ,FUDGE_FLUXCAL_ERROR_FILT, FUDGE_FLUXCAL_ERRPIX_FILT
     &   ,FUDGE_MAG_ERROR, FUDGE_MAG_ERROR_FILT
     &   ,FUDGE_MAG_COVERR, FUDGE_MAG_COVERR_FILT
     &   ,RV_MWCOLORLAW, OPT_MWCOLORLAW, OPT_MWEBV
     &   ,MWEBV_SCALE, MWEBV_SHIFT
     &   ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR
     &   ,CUTWIN_MJD_EXCLUDE, CUTWIN_MJD_INCLUDE

      COMMON / BYTEMASKCOM / MJDMASK
ccc mark delete      COMMON / CIDMASKCOM / MJDMASK,CIDMASK

      EQUIVALENCE
     &    ( cutwin_var(1,cutbit_cid),       cutwin_cid )
     &   ,( cutwin_var(1,cutbit_redshift),  cutwin_redshift )
     &   ,( cutwin_var(1,cutbit_redshift_err),cutwin_redshift_err )
     &   ,( cutwin_var(1,cutbit_ra),        cutwin_ra )
     &   ,( cutwin_var(1,cutbit_dec),       cutwin_dec )
     &   ,( cutwin_var(1,cutbit_hostsep),   cutwin_hostsep )
     &   ,( cutwin_var(1,cutbit_Nepoch),    cutwin_Nepoch )
     &   ,( cutwin_var(1,cutbit_psf),       cutwin_psf )
     &   ,( cutwin_var(1,cutbit_zp),        cutwin_zp  )
     &   ,( cutwin_var(1,cutbit_moondist),  cutwin_moondist  )
     &   ,( cutwin_var(1,cutbit_photprob),  cutwin_photprob  )
     &   ,( cutwin_var(1,cutbit_Nband_thresh), cutwin_nband_thresh )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax), cutwin_nfilt_snrmax )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax2),cutwin_nfilt_snrmax2 )
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmin),cutwin_nfilt_trestmin)
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmax),cutwin_nfilt_trestmax)
     &   ,( cutwin_var(1,cutbit_Nfilt_trest2),  cutwin_nfilt_trest2)
     &   ,( cutwin_var(1,cutbit_Trestmin),  cutwin_Trestmin )
     &   ,( cutwin_var(1,cutbit_Trestmax),  cutwin_Trestmax )
     &   ,( cutwin_var(1,cutbit_TrestRange),cutwin_TrestRange )
     &   ,( cutwin_var(1,cutbit_Trest2),    cutwin_Trest2 )
     &   ,( cutwin_var(1,cutbit_Tgapmax),   cutwin_Tgapmax )
     &   ,( cutwin_var(1,cutbit_T0gapmax),  cutwin_T0gapmax )
     &   ,( cutwin_var(1,cutbit_Tobsmin),   cutwin_Tobsmin )
     &   ,( cutwin_var(1,cutbit_Tobsmax),   cutwin_Tobsmax )
     &   ,( cutwin_var(1,cutbit_Trest),     cutwin_Trest )
     &   ,( cutwin_var(1,cutbit_MJD),       cutwin_MJD   )
     &   ,( cutwin_var(1,cutbit_peakmjd),   cutwin_peakmjd )
     &   ,( cutwin_var(1,cutbit_NMJD_INCLUDE),  cutwin_NMJD_INCLUDE )
     &   ,( cutwin_var(1,cutbit_search),        cutwin_searcheff_mask )
     &   ,( cutwin_var(1,cutbit_snrmax),        cutwin_snrmax  )
     &   ,( cutwin_var(1,cutbit_snrmax2),       cutwin_snrmax2  )
     &   ,( cutwin_var(1,cutbit_nfield),        cutwin_nfield  )
     &   ,( cutwin_var(1,cutbit_mwebv),         cutwin_mwebv  )
     &   ,( cutwin_var(1,cutbit_nseason_active),cutwin_nseason_active)
     &   ,( cutwin_var(1,cutbit_reqep),        cutwin_cutflag_reqep)
     &   ,( cutwin_var(1,cutbit_private),      cutwin_cutflag_private)
     &   ,( cutwin_var(1,cutbit_simvar),       cutwin_cutflag_simvar)
     &   ,( cutwin_var(1,cutbit_usrcuts),      cutwin_cutflag_usrcuts)
     &   ,( cutwin_var(1,cutbit_offset_sbflux+1), cutwin_sbflux_filt)
     &   ,( cutwin_var(1,cutbit_offset_snrmax+1), cutwin_snrmax_filt)
c
     &   ,( cutwin_restlam(1), cutwin_lamrest(1) )



C +CDE,KCORCOM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 

C +CDE,KTABLECM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 
      REAL*4  R4KCORTABLE1D(MXTABLE1D_KCOR)
      COMMON / R4KCORCOM / R4KCORTABLE1D


      INTEGER
     &    LUNKCOR, MXKCOR
     &   ,NKDIM, N4DIM
     &   ,IDMAP_KCOR, IDMAP_AVWARP, IDMAP_LCMAG, IDMAP_MWXT
     &   ,KDIM_T, KDIM_Z, KDIM_AV, KDIM_ifiltr, KDIM_ifilto

      PARAMETER (
     &   LUNKCOR      = 98  ! log unit number to read
     &  ,MXKCOR       = 100 ! max number of KCOR tables to read
     &  ,KDIM_T       = 1  ! epoch KCOR index
     &  ,KDIM_Z       = 2  ! redshift KCOR index
     &  ,KDIM_AV      = 3  ! AV kcor index
     &  ,KDIM_ifiltr  = 4  ! rest-filter KCOR index
     &  ,KDIM_ifilto  = 5  ! obs-filter KCOR index
     &  ,NKDIM        = 5  ! number of KCOR dimensions in map
     &  ,IDMAP_KCOR   = 10 ! for INIT_1DINDEX (any integer is OK)
     &  ,IDMAP_AVWARP = 11
     &  ,IDMAP_LCMAG  = 12
     &  ,IDMAP_MWXT   = 13
     &  ,N4DIM        = 4  ! for 4-d tables
     &   )

c K-cor & photometry-template info

      INTEGER
     &   NCALL_RDKCOR ! # times RDKCOR is called
     &  ,NKCOR_STORE  ! number of stored K-corr tables
     &  ,NZBIN_KCOR   ! # Z-bins in each K-cor table
     &  ,NTBIN_KCOR   ! # T-bins "    "
     &  ,NAVBIN_KCOR  ! # AV bins
     &  ,IAV0         ! IAV index with AV=0
     &  ,IZ0          ! IZ with Z=0
     &  ,NLbin_SNSED      ! # lambda bins for SN SED
     &  ,NTbin_SNSED      ! # epoch bins for SN SED
     &  ,NERR_AVWARP      ! # errors trapped in GET_AVWARP
c
     &  ,NBINTOT_KCOR     ! Total KCOR bins
     &  ,NBINTOT_AVWARP
     &  ,NBINTOT_LCMAG
     &  ,NBINTOT_MWXT
     &  ,NBIN_KCOR(NKDIM)
     &  ,NBIN_AVWARP(NKDIM)
     &  ,NBIN_LCMAG(NKDIM)
     &  ,NBIN_MWXT(NKDIM)

      REAL
     &   Zrange_Kcor(2)    ! min,max redshift in KCOR table
     &  ,Zrange_Kcor_LU(2) ! min,max redshift for lookup
     &  ,Trange_Kcor(2)  ! min,max rest-epoch (days) in KCOR table
     &  ,AVrange_Kcor(2) ! min,max AV
     &  ,ZBINSIZE_KCOR   ! redshift binsize
     &  ,TBINSIZE_KCOR   ! epoch bin size (days)
     &  ,AVbinsize_KCOR  !
     &  ,Cbinsize_AVWARP !
     &  ,GRIDVAL_AV(MXAVBIN_KCOR)  ! store AV(iav) for convenience
     &  ,GRIDVAL_Z(MXZBIN_KCOR)    ! store Z(iz)   for convenience
     &  ,GRIDVAL_T(MXTBIN_KCOR)    ! store Trest(it) for convenience
     &  ,GRIDVAL_C(MXCBIN_AVWARP)  ! store color for AVWARP table
     &  ,Lrange_SNSED(2), LBINSIZE_SNSED
     &  ,Trange_SNSED(2), TBINSIZE_SNSED
     &  ,FLUX_SNSED(MXLAMBIN_SNSED*MXTBIN_KCOR)

      REAL*4
     &   AVWARP_TABLE1D(MXTABLE1D_AVWARP)
     &  ,LCMAG_TABLE1D(MXTABLE1D_LCMAG)
     &  ,MWXT_TABLE1D(MXTABLE1D_MWXT)  ! MilkyWay extinction table

      INTEGER
     &   HIDMAG(MXAVBIN_KCOR,MXFILT_ALL)
     &  ,HIDKCOR(MXFILT_ALL,MXFILT_ALL)

      LOGICAL
     &   USE_AVWARPTABLE
     &  ,LRDZPOFF
     &  ,RDKCOR_STANDALONE
     &  ,EXIST_KCOR(MXFILT_ALL,MXFILT_ALL)

      CHARACTER RESTKCOR_FILTERS*20    ! rest-frame filters from KCOR file

      REAL RVMW_RDKCOR ! RV used for MW table extension
      INTEGER OPT_MWCOLORLAW_RDKCOR

      COMMON / KCORCOM / NCALL_RDKCOR, NKCOR_STORE, RDKCOR_STANDALONE
     &      ,NBIN_KCOR,   NBINTOT_KCOR
     &      ,NBIN_AVWARP, NBINTOT_AVWARP
     &      ,NBIN_LCMAG,  NBINTOT_LCMAG
     &      ,NBIN_MWXT,   NBINTOT_MWXT
     &      ,NZBIN_KCOR,    NTBIN_KCOR,    NAVBIN_KCOR
     &      ,Zrange_KCOR,   Trange_KCOR,   AVrange_KCOR
     &      ,ZBINSIZE_KCOR, TBINSIZE_KCOR, AVbinsize_KCOR
     &      ,CBINSIZE_AVWARP, IAV0, IZ0
     &      ,GRIDVAL_Z, GRIDVAL_T, GRIDVAL_AV, GRIDVAL_C
     &      ,HIDKCOR, HIDMAG
     &      ,Zrange_KCOR_LU, RESTKCOR_FILTERS
     &      ,USE_AVWARPTABLE, LRDZPOFF, NERR_AVWARP
     &      ,NLbin_SNSED, Lrange_SNSED, LBINSIZE_SNSED
     &      ,NTbin_SNSED, Trange_SNSED, TBINSIZE_SNSED
     &      ,FLUX_SNSED
     &      ,LCMAG_TABLE1D, MWXT_TABLE1D, AVWARP_TABLE1D
     &      ,EXIST_KCOR
     &      ,RVMW_RDKCOR, OPT_MWCOLORLAW_RDKCOR

c  July 2016 - define optional spectrograph info
      CHARACTER
     &   SPECTROGRAPH_INSTRUMENT*80
     &  ,SPECTROGRAPH_FILTERLIST*(MXFILT_ALL)
      INTEGER
     &   NFILTDEF_SPECTROGRAPH
     &  ,IFILTDEF_SPECTROGRAPH(MXFILT_ALL)

      COMMON / SPECTROGRAPH_COM /
     &   SPECTROGRAPH_INSTRUMENT, SPECTROGRAPH_FILTERLIST
     &  ,NFILTDEF_SPECTROGRAPH,   IFILTDEF_SPECTROGRAPH

c define temporary 'XXX_RDKCOR' arrays that are used only
c for reading and parsing.

      CHARACTER
     &   PRIMARY_NAME_RDKCOR(10)*40          ! e.g, AB, BD17 ..
     &  ,FILTER_NAME_RDKCOR(MXFILT_ALL)*40   ! e.g, SDSS-g, Bessell90-V
     &  ,KCORINFO_STRING_RDKCOR(MXKCOR)*60   ! Kcor K_xy obsfilt restfilt

      INTEGER
     &   IVER_RDKCOR
     &  ,NPRIM_RDKCOR, NKCOR_RDKCOR
     &  ,INDX_PRIMARY_RDKCOR(MXFILT_ALL)
     &  ,NFILTDEF_RDKCOR               ! Number of all filters
     &  ,IFILTDEF_RDKCOR(MXFILT_ALL)   ! vs. ifilt_rdkcor
     &  ,NFILTOBS_RDKCOR               ! number of obs filters
     &  ,IFILTOBS_RDKCOR(MXKCOR)       ! vs. ikcor
     &  ,IFILT2_RDKCOR(2,MXKCOR)       ! 1=rest:2=obs, 1:NKCOR_STORE
     &  ,MASKFILT_RDKCOR(MXFILT_ALL)   ! bit0,1 -> rest,obs
     &  ,NFILT_DUPLICATE_RDKCOR        ! keep count of duplicate filters
     &  ,IKCOR_RDKCOR(MXKCOR)          ! ikcor vs. istore
     &  ,IPRIM_REF_RDKCOR              ! index for primary  ref.

      REAL
     &   MAG_PRIMARY_RDKCOR(MXFILT_ALL)    ! native mags
     &  ,ZPOFF_PRIMARY_RDKCOR(MXFILT_ALL)  ! mag(native) - mag(synth)
     &  ,ZPOFF_SNPHOT_RDKCOR(MXFILT_ALL)   ! apply to data (from ZPOFF.DAT)

      LOGICAL
     &   ISFITS_RDKCOR
     &  ,LVERBOSE_RDKCOR
     &  ,ISLAMSHIFT_RDKCOR(0:MXFILT_ALL)

      REAL      NULLF_RDKCOR
      INTEGER   NULLI_RDKCOR
      CHARACTER NULLS_RDKCOR*20
      LOGICAL   ANYF_RDKCOR

      COMMON / RDKCORCOM /
     &   IVER_RDKCOR, ISFITS_RDKCOR,LVERBOSE_RDKCOR,ISLAMSHIFT_RDKCOR
     &  ,NKCOR_RDKCOR
     &  ,NPRIM_RDKCOR, INDX_PRIMARY_RDKCOR, PRIMARY_NAME_RDKCOR
     &  ,FILTER_NAME_RDKCOR, KCORINFO_STRING_RDKCOR
     &  ,IFILTDEF_RDKCOR,  NFILTDEF_RDKCOR
     &  ,IFILTOBS_RDKCOR,  NFILTOBS_RDKCOR
     &  ,IFILT2_RDKCOR
     &  ,NFILT_DUPLICATE_RDKCOR, IKCOR_RDKCOR, MASKFILT_RDKCOR
     &  ,IPRIM_REF_RDKCOR
c
     &  ,MAG_PRIMARY_RDKCOR
     &  ,ZPOFF_PRIMARY_RDKCOR
     &  ,ZPOFF_SNPHOT_RDKCOR
c
     &  ,NULLF_RDKCOR,NULLI_RDKCOR,NULLS_RDKCOR,ANYF_RDKCOR


C +CDE,SNHOSTCOM. inserted below by ypatchy. 

c Host galaxy parameters.
c Dec 17, 2012 - add SNHOST_MAGOBS


c logical flags for CWNT
      LOGICAL
     &   EXIST_SNHOST_ANGSEP
     &  ,EXIST_SNHOST_DDLR
     &  ,EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT
     &  ,EXIST_SNHOST_LOGMASS
     &  ,EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS
     &  ,EXIST_SNHOST_SB

      INTEGER*8  SNHOST_OBJID(MXSNHOST)    ! int id
      REAL*8    DSNHOST_OBJID(MXSNHOST)   ! for tables only

      REAL
     &   SNHOST_ANGSEP(MXSNHOST)      ! SN-host sep, arcsec
     &  ,SNHOST_DDLR(MXSNHOST)        ! SNSEP/DLR
     &  ,SNHOST_CONFUSION             ! HC analog from Gupta 2016
     &  ,SNHOST_ZPHOT(MXSNHOST), SNHOST_ZPHOT_ERR(MXSNHOST)
     &  ,SNHOST_ZSPEC(MXSNHOST), SNHOST_ZSPEC_ERR(MXSNHOST)
     &  ,SNHOST_LOGMASS(MXSNHOST)
     &  ,SNHOST_LOGMASS_ERR(MXSNHOST)
     &  ,SNHOST_sSFR(MXSNHOST)
     &  ,SNHOST_sSFR_ERR(MXSNHOST)
     &  ,SNHOST_SBFLUXCAL(MXFILT_ALL)  ! surface brightness FLUXCAL /asec^2
     &  ,SNHOST_MAGOBS(MXFILT_ALL,(MXSNHOST))       ! observer-frame mags
     &  ,SNHOST_MAGOBS_ERR(MXFILT_ALL,(MXSNHOST))   ! error on above

      REAL*8
     &    SNHOST8_RA(MXSNHOST)
     &   ,SNHOST8_DEC(MXSNHOST)

c quantites which do not depend on which host
      INTEGER*4  SNHOST_NMATCH, SNHOST_NMATCH2  ! number of host matches
      REAL
     &   SNHOST_SBFLUXCAL_ERR(MXFILT_ALL)
     &  ,SNHOST_SBMAG(MXFILT_ALL)        ! surface brightness mag/asec^2

      COMMON / SNHOSTCOM /
     &   EXIST_SNHOST_ANGSEP, EXIST_SNHOST_DDLR, EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT,  EXIST_SNHOST_LOGMASS, EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS, EXIST_SNHOST_SB
c
     &  ,SNHOST_NMATCH, SNHOST_NMATCH2
     &  ,SNHOST_ANGSEP, SNHOST_DDLR, SNHOST_CONFUSION
     &  ,SNHOST_ZPHOT,   SNHOST_ZPHOT_ERR
     &  ,SNHOST_ZSPEC,   SNHOST_ZSPEC_ERR
     &  ,SNHOST_LOGMASS, SNHOST_LOGMASS_ERR
     &  ,SNHOST_sSFR,    SNHOST_sSFR_ERR
     &  ,SNHOST_SBFLUXCAL, SNHOST_SBFLUXCAL_ERR, SNHOST_SBMAG
     &  ,SNHOST_MAGOBS, SNHOST_MAGOBS_ERR

      COMMON / SNHOSTCOM8 /
     &   SNHOST_OBJID, DSNHOST_OBJID, SNHOST8_RA, SNHOST8_DEC

cc +CDE,PARSECOM.

C +CDE,SNSIMCOM. inserted below by ypatchy. 

c simulation parameters (if FAKE=2)

      REAL*8  SIM8_RA, SIM8_DECL

      REAL
     &   SIM_REDSHIFT_HELIO
     &  ,SIM_REDSHIFT_CMB
     &  ,SIM_REDSHIFT_HOST
     &  ,SIM_VPEC
     &  ,SIM_DLMAG
     &  ,SIM_LENSDMU
     &  ,SIM_MWEBV
     &  ,SIM_MWRV
     &  ,SIM_SALT2x0
     &  ,SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW
     &  ,SIM_PEAKMJD
     &  ,SIM_EXPOSURE_TIME(MXFILT_OBS)  ! relative exposure time
     &  ,SIM_PEAKMAG(MXFILT_OBS)
     &  ,SIM_EPMAGOBS(MXEPOCH)      ! true epoch mag at each filter/epoch
     &  ,SIM_EPFLUXCAL(MXEPOCH)     ! true fluxcal at each filter/epoch
     &  ,SIM_EPSNRMON(MXEPOCH)      ! optional SNR at MAGMONITOR
     &  ,SIM_EPMAGREST(MXEPOCH)     ! true rest-frame mag
     &  ,SIM_EPFLUXCAL_HOSTERR(MXEPOCH) ! true error from host noise
     &  ,SIMSED_PARVAL(MXPAR_SIMSED)
     &  ,BYOSED_PARVAL(MXPAR_SIMSED)    ! Dec 10 2018
     &  ,LCLIB_PARVAL(MXPAR_LCLIB)
     &  ,SIM_HOSTLIB_PARVAL(MXPAR_SIMSED) ! HOSTLIB params
     &  ,SIM_MAGSMEAR_COH
     &  ,SIM_TEMPLATEMAG(MXFILT_ALL)  ! image-sub template, not LCLIB template
     &  ,SIM_LCWIDTH(MXFILT_ALL)      ! computed from SIM_EPFLUXCAL
c
     &  ,SIM_MODELGRID_TOBS(MXEP_MODELGRID)    ! for SNANA+SIM_MAGOBS table
     &  ,SIM_MODELGRID_MAGOBS(MXEP_MODELGRID,MXFILT_OBS/2) ! idem

      INTEGER
     &   SIM_MODEL_INDEX      ! model index (MLCS,SALT2,NON1a ...)
     &  ,SIM_GENTYPE          ! generated SNTYPE; SIM_TYPE_INDEX
     &  ,SIM_TEMPLATE_INDEX   ! template index for NON1ASED, SIMSED, LCLIB ...
     &  ,SIM_SEARCHEFF_MASK   ! bits 1,2 => found by software,humans
     &  ,SIM_LIBID            ! LIBID for each event
     &  ,SIM_NGEN_LIBID       ! NGEN for this LIBID (usually 1)
     &  ,SIM_NOBS_UNDEFINED   ! NOBS where model is undefined
     &  ,SIM_NSUBSAMPLE_MARK  ! Number of marked sub-samples
     &  ,SIM_SUBSAMPLE_INDEX  ! sub-sample index
     &  ,SIM_REDSHIFT_FLAG    ! points to source of redshift
     &  ,SIMOPT_MWCOLORLAW    ! option for MW color law
     &  ,SIMOPT_MWEBV         ! option to modify MWEBV_SFD map
     &  ,NPAR_SIMSED
     &  ,NPAR_BYOSED
     &  ,NPAR_LCLIB
     &  ,NPAR_SIM_HOSTLIB
     &  ,SIM_EPFILTREST(MXEPOCH)
     &  ,SIMLIB_MSKOPT        ! SIMLIB option mask (Dec 2015)
     &  ,NEP_SIM_MODELGRID

      INTEGER*8  SIM_HOSTLIB_GALID

      CHARACTER
     &   SIMNAME_MODEL*40     ! SALT2, mlcs2k2, NON1A, etc ...
     &  ,SIMNAME_TYPE*12      ! Ia, Ib, II, IIN, etc ...
     &  ,SIMNAME_SHAPEPAR*40  ! DELTA, x1, stretch ...
     &  ,SIMNAME_SHAPELAW*40  ! alpha ...
     &  ,SIMNAME_COLORPAR*40  ! AV, c ...
     &  ,SIMNAME_COLORLAW*40  ! RV,  beta ...
     &  ,SIMLIB_FILENAME*200  ! SIMLIB file used to generate sim
     &  ,SIMSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,SIMSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,BYOSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,BYOSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,LCLIB_KEYWORD(MXPAR_LCLIB)*80 ! full keyname in header
     &  ,LCLIB_PARNAME(MXPAR_LCLIB)*20 ! parname for fitres and ntuple
     &  ,SIM_HOSTLIB_KEYWORD(MXPAR_SIMSED)*60 ! full keyname
     &  ,SIM_HOSTLIB_PARNAME(MXPAR_SIMSED)*60 ! parname
     &  ,SIMNAME_SNRMON*40

      COMMON / SNSIMCOM / SIMNAME_MODEL, SIM_MODEL_INDEX
     &  ,SIM_REDSHIFT_HELIO, SIM_REDSHIFT_CMB, SIM_REDSHIFT_HOST
     &  ,SIM_REDSHIFT_FLAG
     &  ,SIM_VPEC, SIM_DLMAG, SIM_LENSDMU, SIM_SALT2x0, SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW, SIM_PEAKMJD, SIM_MWEBV,SIM_MWRV
     &  ,SIM_EXPOSURE_TIME, SIM_TEMPLATE_INDEX, SIM_TEMPLATEMAG
     &  ,SIM_GENTYPE, SIMNAME_TYPE, SIM_LCWIDTH
     &  ,SIM_SEARCHEFF_MASK, SIM_LIBID, SIM_NGEN_LIBID
     &  ,SIM_NOBS_UNDEFINED, SIM_NSUBSAMPLE_MARK, SIM_SUBSAMPLE_INDEX
     &  ,SIM_PEAKMAG, SIM_EPMAGOBS, SIM_EPFLUXCAL
     &  ,SIM_EPSNRMON, SIM_EPFLUXCAL_HOSTERR
     &  ,SIMNAME_SHAPEPAR,  SIMNAME_SHAPELAW
     &  ,SIMNAME_COLORPAR,  SIMNAME_COLORLAW, SIMNAME_SNRMON
     &  ,SIMLIB_FILENAME, SIMLIB_MSKOPT
     &  ,SIMOPT_MWCOLORLAW, SIMOPT_MWEBV
     &  ,SIM_EPFILTREST, SIM_EPMAGREST,  SIM_MAGSMEAR_COH
     &  ,NEP_SIM_MODELGRID, SIM_MODELGRID_TOBS, SIM_MODELGRID_MAGOBS

      COMMON / SIMSEDCOM /
     &  NPAR_SIMSED, SIMSED_PARVAL, SIMSED_PARNAME, SIMSED_KEYWORD

      COMMON / BYOSEDCOM /
     &   NPAR_BYOSED, BYOSED_PARVAL, BYOSED_PARNAME, BYOSED_KEYWORD

      COMMON / LCLIBCOM /
     &   NPAR_LCLIB, LCLIB_PARVAL, LCLIB_KEYWORD, LCLIB_PARNAME

      COMMON / SIMHOSTCOM /
     &    NPAR_SIM_HOSTLIB, SIM_HOSTLIB_PARVAL
     &   ,SIM_HOSTLIB_PARNAME, SIM_HOSTLIB_KEYWORD

      COMMON / SNSIMCOM8 /
     &         SIM8_RA, SIM8_DECL, SIM_HOSTLIB_GALID



C +CDE,SNLCINP. inserted below by ypatchy. 

c define general input to program:
c SN flux version, hbook files, cut-windows, etc ...
c Dec 30 2012: add KCOR_FILE for fits-formatted tables (same name as in sim)

      CHARACTER*(MXCHAR_FILENAME)  ! user input files
     &   NMLFILE          ! name of input namelist file
     &  ,HFILE_OUT        ! output hbook filename
     &  ,ROOTFILE_OUT     ! output filename for root
     &  ,TEXTFILE_PREFIX  ! prefix for ascii tables
     &  ,SNTABLE_LIST     ! list of SNTABLEs to create
     &  ,SNTABLE_FILTER_REMAP ! remap filters in tables
     &  ,KCOR_FILE            ! fits-formatted KCOR tables
     &  ,hfile_kcor           ! OBSOLETE input -> abort
     &  ,EPOCH_IGNORE_FILE    ! file with epochs to ignor
     &  ,SNMJD_LIST_FILE      ! list of "CID MJD" to process
     &  ,SNMJD_OUT_FILE       ! MJD_LIST out file
     &  ,MNFIT_PKMJD_LOGFILE  ! separate log for PKMJD fits
     &  ,ZPOFF_FILE           ! AB (ZP) off for obs-frame mags
     &  ,USERTAGS_FILE        ! optional int user tag per SN
     &  ,VPEC_FILE            ! pec. velocity corrections
     &  ,HEADER_OVERRIDE_FILE ! header-var overrides
     &  ,SNCID_LIST_FILE      ! list of CID's to process
     &  ,SIMLIB_OUT           ! write simlib entry for each event
     &  ,SIMLIB_ZPERR_LIST    ! 'abc .01 def .02 ghi .005'
     &  ,NONLINEARITY_FILE    ! non-linearity map (Apr 2017)
c
     &  ,FUDGE_HOSTNOISE_FILE  ! legacy: inflate error vs.  hostSB (Mar 2015)
     &  ,FLUXERRMODEL_FILE     ! DATA ONLY: err-fudge maps
     &  ,SIM_FLUXERRMODEL_FILE ! idem, SIM only
     &  ,MAGCOR_FILE           ! DATA ONLY: mag-cor for each CID-MJD-band
     &  ,SIM_MAGCOR_FILE       ! idem, SIM only

      CHARACTER  ! versions
     &   VERSION_PHOTOMETRY(MXVERS)*(MXCHAR_VERSION)   ! SN versions to read
     &  ,VERSION_REPLACE*(MXCHAR_VERSION)  ! replace previous cmd-line VERSION
     &  ,VERSION_REFORMAT_FITS*(MXCHAR_VERSION)   ! output version for
                                                  ! for translating -> FITS

      CHARACTER ! paths
     &   PRIVATE_DATA_PATH*(MXCHAR_PATH)    ! private data subdir
     &  ,FILTER_UPDATE_PATH*(MXCHAR_PATH)   ! SN-dependent filter response

      CHARACTER  ! misc
     &   LINE_ARGS(MXLINE_ARGS)*(MXCHAR_FILENAME)    ! command line args
     &  ,REFORMAT_KEYS*(MXCHAR_FILENAME)  ! global reformat info
     &  ,NONSURVEY_FILTERS*(MXFILT_ALL)   ! non-survey filters to add
     &  ,SNRMAX_FILTERS*(MXFILT_ALL)      ! list of filters for SNRMAX cuts
     &  ,FILTER_REPLACE*(MXFILT_ALL)      ! e.g., 'UGRIZ -> ugriz'
     &  ,FILTLIST_LAMSHIFT*(MXFILT_ALL)   ! list of lam-shifted filters
     &  ,PRIVATE_CUTWIN_STRING(MXCUT_PRIVATE)*(MXCHAR_CUTNAME)
     &  ,PRIVATE_REDSHIFT_CMB*60          ! use PRIVATE variable for redshift
     &  ,SIMVAR_CUTWIN_STRING*(MXCHAR_CUTNAME) ! cuts on SIM_XXX
     &  ,EARLYLC_STRING*(MXCHAR_CUTNAME)  ! see manual
     &  ,REQUIRE_EPOCHS_STRING*100  ! e.g., 'riz 10 7 20' ! uses CUTWIN_SNRMAX
     &  ,DUMP_STRING*100             ! 'funName CID-list'

      INTEGER
     &   NFIT_ITERATION       ! (I) number of fit iterations
     &  ,MINUIT_PRINT_LEVEL   ! -1 -> none
     &  ,INTERP_OPT           ! (I)  interp option (see INTERP_XXX params)
     &  ,NLINE_ARGS
     &  ,OPT_SETPKMJD            ! bit-mask options  to determined PKMJD
     &  ,MJDRANGE_SETPKMJD(2)    ! MJD range to use for SETPKMJD
     &  ,OPT_REFORMAT_TEXT       ! T=> re-write SNDATA files in text format
     &  ,OPT_REFORMAT_TERSE      ! T=> legacy key for text format
     &  ,OPT_REFORMAT_FITS       ! T=> re-write  in FITSformat (data only)
     &  ,OPT_REFORMAT_SALT2      ! T=> re-write SNDATA files in SALT2 format
     &  ,OPTSIM_LCWIDTH          ! T=> option to compute LC width
     &  ,OPT_DEBUG               ! T=> for internal debug
     &  ,JOBSPLIT(2)             ! process (1)-range of (2)=TOTAL
     &  ,JOBSPLIT_EXTERNAL(2)    ! passed by split_and_fit for text format
     &  ,MXLC_FIT                ! stop after this many fits passing all cuts
     &  ,PHOTFLAG_DETECT         ! used to count NOBS_DETECT and TLIVE_DETECT
     &  ,APPLY_FLUXCOR_BUG       ! see calling sequence in SNRECON
     &  ,FLUXERRMODEL_OPTMASK

      LOGICAL
     &   LDMP_SNFAIL        ! (I) T => dump reason for each SN failure
     &  ,LSIM_SEARCH_SPEC   ! (I) T => require simulated SPEC-tag
     &  ,LSIM_SEARCH_ZHOST  ! (I) T => require simulated zHOST
     &  ,LTEST_INTERP  ! (I) T => calls TEST_KCOR
     &  ,LTEST_KCOR    ! (I) T => calls TEST_INTERP
     &  ,LTEST_MAG     ! (I) T => test GET_MAGLC  and GET_MWXT
     &  ,LTEST_U3BAND  ! (I) T => require exactly 3 bands that include U
     &  ,USE_LINE_ARGS(MXLINE_ARGS)
     &  ,USE_SNHOST_ZPHOT    ! (I) T=> replace SNLC_REDSHIFT -> SNHOST_ZPHOT
c
     &  ,ABORT_ON_NOEPOCHS   ! T=> abort if there are no epochs to fit
     &  ,ABORT_ON_BADAVWARP  ! T=> abort if AVwarp cannot be determined
     &  ,ABORT_ON_BADZ       ! T=> abort on bad z in GET_KCOR8
     &  ,ABORT_ON_BADKCOR    ! T=> affects only the init stage (RDKCOR)
     &  ,ABORT_ON_BADSURVEY  ! T=> abort if SURVEY changes
     &  ,ABORT_ON_BADFILTER  ! T=> abort if fit band is not in SURVEY filters
     &  ,ABORT_ON_MARGPDF0   ! T=> abort if marginalized pdf=0 everywhere
     &  ,ABORT_ON_TRESTCUT   ! T=> abort if any Trest cut is set (for photoz)
     &  ,ABORT_ON_DUPLCID    ! T=> abort on duplicate CID
     &  ,ABORT_ON_DUPLMJD    ! T=> abort on repeat MJD+band (Jun 2017)
     &  ,ABORT_ON_NOPKMJD    ! T=> abort if no PKMJDINI (see OPT_SETPKMJD)
     &  ,ABORT_ON_BADHEADER  ! T=> abort if missing header info(HEADMASK_CHECK)
     &  ,USE_MINOS           ! T=> use MINOS instead of MINIMIZE
     &  ,LDMP_AVWARP         ! dump GET_AVWARP8 (debug only)
     &  ,LDMP_KCORFUN        ! dump KCORFUN8 (debug only)
     &  ,USESIM_SNIA         ! default True -> process simulated SNIa
     &  ,USESIM_SNCC         ! default True -> process simulated SNCC
     &  ,USESIM_TRUEFLUX     ! SNLC_FLUXCAL -> SIM_FLUXCAL

c LC plot info

      INTEGER
     &  MXLC_PLOT      ! (I) max number of plots to make with SNLCPAK
     & ,NCCID_PLOT     ! (I) size of SNCCID_PLOT array
      CHARACTER
     &  SNCCID_PLOT(MXLISTNML)*(MXCHAR_CCID) ! (I) string-list of CCIDs to plot

      REAL
     &  DTOBS_MODEL_PLOT  ! (I) binning (days) for overlaid best-fit model

c Oct 2014: variables to check for multi-season transient activity

      INTEGER
     &   MULTISEASON_OPTMASK      ! option(s) for GET_MULTISEASON
      REAL
     &   MULTISEASON_TGAP              ! new season if time-gap > this
     &  ,MULTISEASON_NREJECT_OUTLIER   ! num outliers to reject per season
     &  ,MULTISEASON_CHI2RED_ACTIVE    ! min chi2red for active season

      REAL
     &   SIM_PRESCALE        ! pre-scale applied only to SIM
     &  ,VPEC_ERR_OVERRIDE   ! override VPEC_ERR in data header (km/sec)


c define reference cosmological parameters (double precision!)
c Aug 7, 2009: each parameter has dimension 2 for value & error

      REAL*8
     &   H0_REF(2)    ! reference H0  (70 => 70 MPc^-1)
     &  ,OLAM_REF(2)  ! omega_matter
     &  ,OMAT_REF(2)  ! omega_lamgda
     &  ,ORAD_REF(2)  ! omega_rad
     &  ,W0_REF(2)    ! w = p/rho
     &  ,DWDA_REF(2)  ! dw/da [a = 1/(1+z)]

      INTEGER   NSIMVAR_CUTWIN
      REAL      SIMVAR_CUTWIN(2,MXCUT_PRIVATE)
      CHARACTER SIMVAR_CUTWIN_LIST(MXCUT_PRIVATE)*40

      COMMON   / SNLCINP /
     &      NLINE_ARGS, LINE_ARGS, USE_LINE_ARGS, nmlfile
     &    , VERSION_PHOTOMETRY, VERSION_REPLACE, VERSION_REFORMAT_FITS
     &    , JOBSPLIT, JOBSPLIT_EXTERNAL, SIM_PRESCALE, MXLC_FIT
     &    , PRIVATE_DATA_PATH, FILTER_UPDATE_PATH
     &    , NONSURVEY_FILTERS, SNRMAX_FILTERS, VPEC_ERR_OVERRIDE
     &    , FILTER_REPLACE, FILTLIST_LAMSHIFT
     &    , OPTSIM_LCWIDTH, OPT_REFORMAT_TERSE, OPT_REFORMAT_TEXT
     &    , OPT_REFORMAT_SALT2, REFORMAT_KEYS, OPT_REFORMAT_FITS
     &    , SNMJD_LIST_FILE, SNMJD_OUT_FILE, MNFIT_PKMJD_LOGFILE
     &    , EPOCH_IGNORE_FILE, SNCID_LIST_FILE
     &    , SIMLIB_OUT, SIMLIB_ZPERR_LIST, NONLINEARITY_FILE
     &    , hfile_out, rootfile_out, textfile_prefix
     &    , SNTABLE_LIST, SNTABLE_FILTER_REMAP
     &    , hfile_kcor, KCOR_FILE, ZPOFF_FILE, USERTAGS_FILE
     &    , VPEC_FILE, HEADER_OVERRIDE_FILE
     &    , FLUXERRMODEL_FILE,SIM_FLUXERRMODEL_FILE,FLUXERRMODEL_OPTMASK
     &    , MAGCOR_FILE, SIM_MAGCOR_FILE,  FUDGE_HOSTNOISE_FILE
     &    , NFIT_ITERATION, MINUIT_PRINT_LEVEL, INTERP_OPT, USE_MINOS
     &    , OPT_SETPKMJD, MJDRANGE_SETPKMJD, OPT_DEBUG
     &    , LDMP_SNFAIL, LSIM_SEARCH_SPEC, LSIM_SEARCH_zHOST
     &    , LDMP_AVWARP, LDMP_KCORFUN
     &    , LTEST_KCOR, LTEST_INTERP, LTEST_MAG, LTEST_U3BAND
     &    , USESIM_SNIA, USESIM_SNCC, USESIM_TRUEFLUX, USE_SNHOST_ZPHOT
     &    , ABORT_ON_NOEPOCHS, ABORT_ON_BADAVWARP, ABORT_ON_NOPKMJD
     &    , ABORT_ON_BADZ, ABORT_ON_BADKCOR, ABORT_ON_BADSURVEY
     &    , ABORT_ON_MARGPDF0, ABORT_ON_TRESTCUT
     &    , ABORT_ON_DUPLCID, ABORT_ON_DUPLMJD
     &    , ABORT_ON_BADHEADER, ABORT_ON_BADFILTER
     &    , SNCUT_SNRMAX, SNCUT_HOST_SBFLUX
     &    , EPCUT_SNRMIN, EPCUT_SKYSIG, EPCUT_PSFSIG
     &    , CUTWIN_LAMREST, CUTWIN_LAMOBS
     &    , PRIVATE_CUTWIN_STRING, PRIVATE_REDSHIFT_CMB
     &    , SIMVAR_CUTWIN_STRING
     &    , NSIMVAR_CUTWIN, SIMVAR_CUTWIN, SIMVAR_CUTWIN_LIST
     &    , EARLYLC_STRING, REQUIRE_EPOCHS_STRING, DUMP_STRING
     &    , MXLC_PLOT, NCCID_PLOT, SNCCID_PLOT, DTOBS_MODEL_PLOT
     &    , MULTISEASON_OPTMASK, MULTISEASON_TGAP
     &    , MULTISEASON_NREJECT_OUTLIER, MULTISEASON_CHI2RED_ACTIVE
     &    , APPLY_FLUXCOR_BUG, PHOTFLAG_DETECT
     &    , CUTWIN_ZPADU, CUTWIN_ZPNPE

      COMMON / SNLCINP8 /
     &    H0_REF, OLAM_REF, OMAT_REF, W0_REF, DWDA_REF

      NAMELIST / SNLCINP /
     &      VERSION_PHOTOMETRY, VERSION_REPLACE, VERSION_REFORMAT_FITS
     &    , PRIVATE_DATA_PATH, FILTER_UPDATE_PATH
     &    , NONSURVEY_FILTERS, SNRMAX_FILTERS, VPEC_ERR_OVERRIDE
     &    , FILTER_REPLACE, FILTLIST_LAMSHIFT
     &    , JOBSPLIT, JOBSPLIT_EXTERNAL, SIM_PRESCALE, MXLC_FIT
     &    , OPTSIM_LCWIDTH, OPT_REFORMAT_TERSE, OPT_REFORMAT_TEXT
     &    , OPT_REFORMAT_SALT2, REFORMAT_KEYS, OPT_REFORMAT_FITS
     &    , SNMJD_LIST_FILE, SNMJD_OUT_FILE, MNFIT_PKMJD_LOGFILE
     &    , hfile_out, rootfile_out, textfile_prefix
     &    , SNTABLE_LIST, SNTABLE_FILTER_REMAP
     &    , hfile_kcor, KCOR_FILE, ZPOFF_FILE, USERTAGS_FILE
     &    , VPEC_FILE, HEADER_OVERRIDE_FILE
     &    , EPOCH_IGNORE_FILE, SNCID_LIST_FILE
     &    , SIMLIB_OUT, SIMLIB_ZPERR_LIST, NONLINEARITY_FILE
     &    , NFIT_ITERATION, MINUIT_PRINT_LEVEL, INTERP_OPT, USE_MINOS
     &    , OPT_SETPKMJD, MJDRANGE_SETPKMJD, OPT_DEBUG
     &    , LDMP_SNFAIL, LSIM_SEARCH_SPEC, LSIM_SEARCH_ZHOST
     &    , LTEST_KCOR, LTEST_INTERP, LTEST_MAG, LTEST_U3BAND
     &    , USESIM_SNIA, USESIM_SNCC, USESIM_TRUEFLUX, USE_SNHOST_ZPHOT
     &    , ABORT_ON_NOEPOCHS, ABORT_ON_BADAVWARP, ABORT_ON_NOPKMJD
     &    , ABORT_ON_BADZ, ABORT_ON_BADKCOR, ABORT_ON_BADSURVEY
     &    , ABORT_ON_MARGPDF0, ABORT_ON_TRESTCUT
     &    , ABORT_ON_DUPLCID, ABORT_ON_DUPLMJD
     &    , ABORT_ON_BADHEADER, ABORT_ON_BADFILTER
     &    , H0_REF, OLAM_REF, OMAT_REF, W0_REF, DWDA_REF
     &    , USE_MWCOR, DOBUG_LAMRANGE
     &    , MXLC_PLOT, NCCID_PLOT, SNCCID_PLOT, DTOBS_MODEL_PLOT
     &    , MULTISEASON_OPTMASK, MULTISEASON_TGAP
     &    , MULTISEASON_NREJECT_OUTLIER, MULTISEASON_CHI2RED_ACTIVE
     &    , APPLY_FLUXCOR_BUG, PHOTFLAG_DETECT
c
c Below are cut-variables defined in SNCUTS
c
     &   ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST
     &   ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE
     &   ,SNTYPE_IGNORE, CCDNUM_LIST, SIM_TEMPLATE_INDEX_LIST
     &   ,SNFIELD_LIST, SNTEL_LIST
     &   ,PHOTFLAG_MSKREJ, cutwin_cid
     &   ,cutwin_redshift, cutwin_redshift_err
     &   ,cutwin_ra, cutwin_dec
     &   ,cutwin_hostsep,   cutwin_Nepoch
     &   ,cutwin_snrmax,    cutwin_snrmax2, cutwin_nfield
     &   ,cutwin_mwebv,     cutwin_nseason_active
     &   ,cutwin_Trestmin,   cutwin_Trestmax
     &   ,cutwin_TrestRange, cutwin_Trest
     &   ,cutwin_Tgapmax,  cutwin_T0gapmax
     &   ,cutwin_Tobsmin,  cutwin_Tobsmax
     &   ,cutwin_peakmjd,  cutwin_mjd
c
     &   ,cutwin_psf,       cutwin_zp,  CUTWIN_ZPADU, CUTWIN_ZPNPE
     &   ,cutwin_moondist,  cutwin_photprob
     &   ,cutwin_nband_thresh
     &   ,cutwin_nfilt_snrmax,   cutwin_nfilt_snrmax2
     &   ,cutwin_nfilt_trestmin, cutwin_nfilt_trestmax
     &   ,cutwin_trest2, cutwin_nfilt_trest2
     &   ,cutwin_restlam, cutwin_lamrest, cutwin_lamobs
     &   ,PRIVATE_CUTWIN_STRING, PRIVATE_REDSHIFT_CMB
     &   ,SIMVAR_CUTWIN_STRING
     &   ,EARLYLC_STRING, REQUIRE_EPOCHS_STRING, DUMP_STRING
c
     &   ,SNCUT_SNRMAX, SNCUT_HOST_SBFLUX
     &   ,EPCUT_SNRMIN, EPCUT_SKYSIG, EPCUT_PSFSIG
     &   ,MAGOBS_SHIFT_PRIMARY,  MAGOBS_SHIFT_ZP
     &   ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_ZP, FILTER_LAMSHIFT
     &   ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS
     &   ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX
     &   ,FUDGE_MAG_ERROR, FUDGE_MAG_COVERR
     &   ,FLUXERRMODEL_FILE,SIM_FLUXERRMODEL_FILE,FLUXERRMODEL_OPTMASK
     &   ,MAGCOR_FILE, SIM_MAGCOR_FILE, FUDGE_HOSTNOISE_FILE
     &   ,RV_MWCOLORLAW, OPT_MWCOLORLAW, OPT_MWEBV
     &   ,MWEBV_SCALE, MWEBV_SHIFT
     &   ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR
     &   ,CUTWIN_MJD_EXCLUDE, CUTWIN_MJD_INCLUDE, CUTWIN_NMJD_INCLUDE


c ---------- end of SNLCINP ---------


C +CDE,FILTCOM. inserted below by ypatchy. 

c filter bandpasses.
c Jan 16 2017: MXLAMBIN_FILT -> 4000 (was 3000) for JWST test

      INTEGER MXLAMBIN_FILT, MXLAMBIN_PRIM
      PARAMETER (
     &   MXLAMBIN_FILT = 4000   ! max number of lambda bins to store
     &  ,MXLAMBIN_PRIM = 4000  ! max lambda bins for primary
     &    )

c Define filter chars A-Z, a-z and 0-9, but do not change order of
c original FILTDEF characters so that we can still use older Kcor files.
c
      CHARACTER FILTDEF_STRING*100
      PARAMETER (FILTDEF_STRING =
     &  'ugrizYJHK UBVRIXy0123456789 abcdef ' //  !  32
     &  'ACDEFGLMNOPQSTWZ ' //                    ! +16
     &  'hjklmnopqstvwx '   //                    ! +14 = 62
     &  '&'                             ! +1(autoSynFilter for spectrograph)
     &         )

      CHARACTER
     &   FILTOBS_NAME(MXFILT_ALL)*40   ! full filternames
     &  ,FILTREST_NAME(MXFILT_ALL)*40  ! full filternames
     &  ,PRIMARY_NAME*40               ! name of primary
     &  ,NONSURVEY_FILTERS_ADD*80      ! NONSURVEY_FILTERS that were added

      LOGICAL
     &   LFILTDEF_OBS(MXFILT_ALL)  ! filter-trans define in KCOR file
     &  ,LFILTDEF_REST(MXFILT_ALL) ! filter-trans define in KCOR file
     &  ,LFILTDEF_NONSURVEY(MXFILT_ALL) ! part of NONSURVEY_FILTERS
     &  ,LFILTDEF_SNRMAX(MXFILT_ALL)    ! use for SNRMAX cut
     &  ,FREEZE_SURVEY_FILTERS     ! T=> do not re-read on 2nd version
     &  ,EXIST_BXFILT_OBS          ! T=> observer BX exists
     &  ,EXIST_BXFILT_REST         ! T=> idem for rest-frame

c define indices for legacy filters
      INTEGER
     &   IFILT_SDSS_u, IFILT_SDSS_g, IFILT_SDSS_r
     &  ,IFILT_SDSS_i, IFILT_SDSS_z
     &  ,IFILT_BESS_U, IFILT_BESS_B, IFILT_BESS_V
     &  ,IFILT_BESS_R, IFILT_BESS_I, IFILT_BESS_BX
     &  ,IFILT_Y, IFILT_J, IFILT_H, IFILT_K

c define filter set for survey (observer) and for rest-frame.
c Each set goes 1 - NFILTDEF_[SURVEY,REST]

      INTEGER
     &   NFILTDEF_SURVEY          ! no. survey (obs) filters
     &  ,NFILTDEF_READ            ! no. filters to read (excludes BX)
     &  ,IFILTDEF_MAP_SURVEY(MXFILT_OBS)  ! IFILT_OBS vs. sparse filter index
     &  ,IFILTDEF_INVMAP_SURVEY(MXFILT_ALL) ! sparse index vs. IFILT_OBS
c
     &  ,NFILTDEF_REST                 ! no. rest-filts defined in Kcor file
     &  ,IFILTDEF_MAP_REST(MXFILT_OBS) ! IFILT_REST vs. sparse filter index
     &  ,IFILTDEF_INVMAP_REST(MXFILT_ALL) ! sparse index vs. IFILT_REST
c
     &  ,NFILTDEF_IGNORE_REST
     &  ,IFILTDEF_IGNORE_REST(MXFILT_OBS)  ! nearest rest-filters to ignore

c variables used in fitter; these variables are over-written
c in FITPAR_PREP for each SN

      INTEGER
     &   IFILT_REST_MAP(MXFILT_OBS)  ! idem for rest-frame filter
     &  ,NFILT_OBS_USEFIT                ! number of filters used in fit

      INTEGER*8
     &  IFILT_OBS_EVAL_MASK(2,MXFILT_ALL) ! set for each filter

      character FILTLIST_FIT_USE*64  ! filter-list USED each fit

c define filter properties using MXFILT_ALL
c Nov 12, 2010: split FILT_XXX into FILTOBS_XXX and FILTREST_XXX
c
      INTEGER
     &   NLAMBIN_FILTOBS(MXFILT_ALL)   ! number of bins per transmission curve
     &  ,NLAMBIN_FILTREST(MXFILT_ALL)  ! number of bins per transmission curve
     &  ,NLAMBIN_PRIMARY               ! lambda bins for primary spec


      REAL
     &   FILTOBS_TRANS(MXLAMBIN_FILT,MXFILT_ALL)   ! filter transmissions
     &  ,FILTOBS_TRANSMAX(MXFILT_ALL)              ! max trans
     &  ,FILTOBS_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)  ! corresponding lambda's (A)
     &  ,FILTOBS_LAMAVG(MXFILT_ALL)           ! effective central lambda
     &  ,FILTOBS_LAMRMS(MXFILT_ALL)           ! lambda-RMS
     &  ,FILTOBS_LAMRANGE(2,MXFILT_ALL)       ! min,max range in obs-frames
     &  ,FILTOBS_MAG_PRIMARY(MXFILT_ALL)      ! primary mag vs. ifilt_obs
     &  ,FILTOBS_ZPOFF_PRIMARY(MXFILT_ALL)    ! mag(native) - mag(synth)
     &  ,FILTOBS_ZPOFF_SNPHOT(MXFILT_ALL)     ! apply these ZPOFF to SNphot
     &  ,FILTOBS_LAMSHIFT(MXFILT_ALL)         ! lambda shift per filter
c
     &  ,FILTREST_TRANS(MXLAMBIN_FILT,MXFILT_ALL)   ! filter transmissions
     &  ,FILTREST_TRANSMAX(MXFILT_ALL)              ! max trans
     &  ,FILTREST_LAMBDA(MXLAMBIN_FILT,MXFILT_ALL)  !  lambda's (A)
     &  ,FILTREST_LAMAVG(MXFILT_ALL)           ! effective central lambda
     &  ,FILTREST_LAMRMS(MXFILT_ALL)           ! lambda-RMS
     &  ,FILTREST_LAMRANGE(2,MXFILT_ALL)       ! min,max range in rest-frames
     &  ,FILTREST_MAG_PRIMARY(MXFILT_ALL)      ! primary mag vs. ifilt
     &  ,FILTREST_ZPOFF_PRIMARY(MXFILT_ALL)    ! mag(native) - mag(synth)
     &  ,FILTREST_LAMSHIFT(MXFILT_ALL)         ! lambda shift per filter
c
     &  ,PRIMARY_FLUX(MXLAMBIN_PRIM)       ! primary spec for obs filters
     &  ,PRIMARY_LAM(MXLAMBIN_PRIM)


c for remap
      INTEGER NFILT_REPLACE, IFILTOBS_REPLACE(MXFILT_ALL)

      COMMON / FILTCOM / FILTOBS_NAME, FILTREST_NAME
     &   ,IFILT_SDSS_u, IFILT_SDSS_g, IFILT_SDSS_r
     &   ,IFILT_SDSS_i, IFILT_SDSS_z
     &   ,IFILT_BESS_U, IFILT_BESS_B, IFILT_BESS_V
     &   ,IFILT_BESS_R, IFILT_BESS_I, IFILT_BESS_BX
     &   ,IFILT_Y, IFILT_J, IFILT_H, IFILT_K
c
     &   ,NFILTDEF_SURVEY,        NFILTDEF_REST, NFILTDEF_READ
     &   ,FREEZE_SURVEY_FILTERS
     &   ,EXIST_BXFILT_OBS,       EXIST_BXFILT_REST
     &   ,IFILTDEF_MAP_SURVEY,    IFILTDEF_MAP_REST
     &   ,IFILTDEF_INVMAP_SURVEY, IFILTDEF_INVMAP_REST
     &   ,NFILTDEF_IGNORE_REST,   IFILTDEF_IGNORE_REST
     &   ,LFILTDEF_OBS, LFILTDEF_REST
     &   ,LFILTDEF_NONSURVEY, NONSURVEY_FILTERS_ADD
     &   ,LFILTDEF_SNRMAX
c
     &   ,NLAMBIN_FILTOBS, FILTOBS_TRANS, FILTOBS_TRANSMAX
     &   ,FILTOBS_LAMBDA, FILTOBS_LAMAVG, FILTOBS_LAMRMS
     &   ,FILTOBS_LAMRANGE
     &   ,FILTOBS_MAG_PRIMARY, FILTOBS_ZPOFF_PRIMARY
     &   ,FILTOBS_ZPOFF_SNPHOT, FILTOBS_LAMSHIFT
c
     &   ,NLAMBIN_FILTREST, FILTREST_TRANS, FILTREST_TRANSMAX
     &   ,FILTREST_LAMBDA,  FILTREST_LAMAVG, FILTREST_LAMRMS
     &   ,FILTREST_MAG_PRIMARY, FILTREST_ZPOFF_PRIMARY
     &   ,FILTREST_LAMSHIFT, FILTREST_LAMRANGE
c
     &   ,IFILT_REST_MAP
     &   ,FILTLIST_FIT_USE, NFILT_OBS_USEFIT
     &   ,PRIMARY_FLUX, PRIMARY_LAM, NLAMBIN_PRIMARY, PRIMARY_NAME
c
     &   ,NFILT_REPLACE, IFILTOBS_REPLACE

      COMMON / FILTCOM8 / IFILT_OBS_EVAL_MASK


c Apr 28 2016: read K_xy and REST_COLOR

      INTEGER i, LL, ilast, iuse
      CHARACTER ctemp*40, cfilt*2, FNAM*20

      INTEGER FILTINDX

C --------------- BEGIN ---------------

      FNAM = 'PARSE_KCORDUMP_ARGS'
      NTBIN_DMP = 1
      NZBIN_DMP = 1

      i = 1
      ilast = 1

      DO WHILE ( i .LE. NLINE_ARGS )

        LL = INDEX ( LINE_ARGS(i), ' ' ) - 1
        print*,'     PROCESS COMMAND LINE ARG: ', LINE_ARGS(i)(1:LL)

         if ( LINE_ARGS(i) .EQ. 'FILT_REST' ) then
            i = i + 1
            READ(LINE_ARGS(i),*) CFILT
            IFILT_REST_DMP = FILTINDX(CFILT)

         else if ( LINE_ARGS(i) .EQ. 'FILT_OBS' ) then
            i = i + 1
            READ(LINE_ARGS(i),*) CFILT
            IFILT_OBS_DMP = FILTINDX(CFILT)

         else if ( LINE_ARGS(i)(1:2) .EQ. 'K_' ) then ! Apr 28 2016
            STRING_Kxy = LINE_ARGS(i)
            CFILT = STRING_Kxy(3:3) ; IFILT_REST_DMP = FILTINDX(CFILT)
            CFILT = STRING_Kxy(4:4) ; IFILT_OBS_DMP  = FILTINDX(CFILT)
            USE_LINE_ARGS(i) = .TRUE.

         else if ( LINE_ARGS(i) .EQ. 'REST_COLOR' ) then ! Apr 28 2016
            i = i + 1
            READ(LINE_ARGS(i),*) STRING_REST_COLOR  ! e.g., 'B-V'
            if ( STRING_REST_COLOR(2:2) .NE. '-' ) then
              c1err = 'Invalid REST_COLOR argument: '
     &                // STRING_REST_COLOR
              c2err = 'Must be of the form X-Y'
              CALL MADABORT(FNAM,C1ERR, C2ERR)
            endif
            CFILT = STRING_REST_COLOR(1:1)
            IFILT_RESTCOLOR_DMP(1) = FILTINDX(CFILT)
            CFILT = STRING_REST_COLOR(3:3)
            IFILT_RESTCOLOR_DMP(2) = FILTINDX(CFILT)

            i = i + 1
            READ(LINE_ARGS(i),*) MAGREST_COLOR
            MAGREST(1) = MAGREST_COLOR
            MAGREST(2) = 0.0

         else if ( LINE_ARGS(i) .EQ. 'TREST' ) then
            i = i + 1

            if ( LINE_ARGS(i) .EQ. 'ALL' ) then
               CALL SET_TREST_ALL
            else
               READ(LINE_ARGS(i),*) TREST_DMP
               TREST_DMP_RANGE(1) = TREST_DMP
               TREST_DMP_RANGE(2) = TREST_DMP
            endif

         else if ( LINE_ARGS(i) .EQ. 'REDSHIFT' .or.
     &             LINE_ARGS(i) .EQ. 'Z'        .or.
     &             LINE_ARGS(i) .EQ. 'z' ) then

            i = i + 1
            if ( LINE_ARGS(i) .EQ. 'ALL' ) then
               REDSHIFT_DMP = 0.0
               NZBIN_DMP    = 100  ! flag to use REDSFHIT range in KCOR file
            else
               READ(LINE_ARGS(i),*) REDSHIFT_DMP
               REDSHIFT_DMP_RANGE(1) = REDSHIFT_DMP
               REDSHIFT_DMP_RANGE(2) = REDSHIFT_DMP
            endif

         else if ( LINE_ARGS(i) .EQ. 'KCOR_FILE' ) then
            i = i + 1
            KCOR_FILE = LINE_ARGS(i)
         endif

c set logical flag for used LINE_ARGS

         IF ( i > ilast ) THEN
           DO iuse = ilast, i
              USE_LINE_ARGS(iuse) = .TRUE.
           ENDDO
         ENDIF

         i = i + 1
         ilast = i

      ENDDO

c make sure that all command line args are valid;
c if not, then abort.

cc      print*,'  USE_LINE_ARGS = ', ( USE_LINE_ARGS(i), i=1, 10)
      CALL CHECK_LINE_ARGS

C ----------------------------------------------
c  prompt user to enter missing arguments that were
c  not given on the command line.

      print*,' '

      IF ( IFILT_REST_DMP .LE. 0 ) THEN
        write(6,50) 'FILT_REST (' // FILTDEF_STRING(1:62) // ')'
        read(5,55) CFILT
        IFILT_REST_DMP = FILTINDX(CFILT)
      ENDIF

      IF ( IFILT_OBS_DMP .LE. 0 ) THEN
        write(6,50) 'IFILT_OBS (' // FILTDEF_STRING(1:62) // ')'
        read(5,55) CFILT
        IFILT_OBS_DMP = FILTINDX(CFILT)
      ENDIF

      IF ( REDSHIFT_DMP .LT. 0.0 ) THEN
        write(6,50) 'REDSHIFT'
        read(5,*) CTEMP
        if ( CTEMP(1:3) .EQ. 'ALL') then
           NZBIN_DMP = 100  ! flag to use REDSFHIT range in KCOR file
        else
           READ(ctemp,*) REDSHIFT_DMP
           REDSHIFT_DMP_RANGE(1) = REDSHIFT_DMP
           REDSHIFT_DMP_RANGE(2) = REDSHIFT_DMP
        endif

      ENDIF

      IF ( TREST_DMP .LT. -90.0 ) THEN
        write(6,50) 'TREST'
        read(5,*) CTEMP
        if ( CTEMP(1:3) .EQ. 'ALL') then
           call SET_TREST_ALL
        else
           READ(ctemp,*) TREST_DMP
           TREST_DMP_RANGE(1) = TREST_DMP
           TREST_DMP_RANGE(2) = TREST_DMP
        endif
      ENDIF

      IF ( KCOR_FILE .EQ. '' ) THEN
        write(6,50) 'KCOR_FILE'
        read(5,55) KCOR_FILE
      ENDIF

50    format(T5,'Enter ',A,' ==> ', $)
51    format(I2)
52    format(F10.5)
55    format(A)

      RETURN
      END

C =============================
CDECK  ID>,  SET_REDSHIFT_ALL. 
      SUBROUTINE SET_REDSHIFT_ALL
      IMPLICIT NONE

C +CDE,SNDATCOM. inserted below by ypatchy. 


C +CDE,SNPAR. inserted below by ypatchy. 

      CHARACTER  SNTABLE_LIST_DEFAULT*60

c parameters used by snana code

      INTEGER
     &   MXVERS, MXSURVEY, MXSNLC,MXCID, MXCID_CHECK,MXEPOCH, MXITER
     &  ,MXFILT_ALL, MXFILT_OBS, MXFILT_REST, ONE
     &  ,MXFILT_KCOR, MXFILT_SNRMAX, MXEP_MODELGRID
     &  ,MXIDTEL, MXIDFIELD, MXFIELD_OVP, MXSEASON, MXSNHOST
     &  ,MXTABLE1D_KCOR,MXTABLE1D_LCMAG,MXTABLE1D_MWXT,MXTABLE1D_AVWARP
     &  ,MXZBIN_KCOR, MXTBIN_KCOR, MXAVBIN_KCOR, MXCBIN_AVWARP
     &  ,MNTYPE, MXTYPE, MXLISTNML, MXCCID_LIST
     &  ,MXVAR_PRIVATE, MXCUT_PRIVATE, MXVAR_TERSE
     &  ,HOFF_NOCUTS, HOFF_CUTS, HOFF_SIM
     &  ,LUNHIS, LUNNML, LUNDAT, LUNTMP, LUNFIT, LUNPKMJD, LUNDMP
     &  ,LUNCID, LUNOUT
     &  ,LUNRES1, LUNRES2, LUNINTERP, LUNLIST, LUNIGNORE, LUNSALT2
     &  ,ISTAGE_INIT, ISTAGE_RDSN, ISTAGE_CUTS, ISTAGE_USRANA
     &  ,ISTAGE_TEST
     &  ,INTERP_LINEAR, INTERP_SMOOTH, INTERP_ZSMOOTH
     &  ,MXERRTYPE, ERRTYPE_MINOS, ERRTYPE_PARAB, ERRTYPE_MARG
     &  ,ERRTYPE_BAD
     &  ,FCNFLAG_USER, FCNFLAG_FAST, FCNFLAG_LAST, FCNFLAG_USESIM
     &  ,FCNFLAG_PRIOR_ONLY, FCNFLAG_SIGMA_ONLY
     &  ,OPT_INTEGPDF_QUITCHI2, OPT_INTEGPDF_FULL
     &  ,OPT_SNXT_CCM89, OPT_SNXT_SJPAR
     &  ,OPT_MWCOLORLAW_DEFAULT, OPT_MWEBV_DEFAULT
     &  ,OPT_KCORERR_SJ, OPT_KCORERR_SJ5, OPT_KCORERR_SMOOTH
     &  ,MXLINE_ARGS, MXEPOCH_IGNORE
     &  ,NPAR_ANYLC, MXPAR_SIMSED, MXPAR_LCLIB
     &  ,MXLAMBIN_SNSED, MXCUTBIT
     &  ,OPT_FILTUPD_EACHSN, OPT_FILTUPD_MAP, OPT_FILTUPD_SAMEFILT
     &  ,OPT_FILTOBS, OPT_FILTREST
     &  ,ISTAT_READAGAIN, ISTAT_SKIP
     &  ,IDTABLE_MCMC, NFIT_VERBOSE
     &  ,ITABLE_SNANA,  ITABLE_FITRES, ITABLE_SNLCPAK, ITABLE_SPECPAK
     &  ,ITABLE_SPECTRA, MXTABLE
     &  ,IDTABLE_SNANA, IDTABLE_FITRES, IDTABLE_SPECTRA
     &  ,IDTABLE_CIDPTR
     &  ,OPT_PARSTORE_TEXTTABLE
c
     &  ,MODEL_STRETCH
     &  ,MODEL_STRETCH2
     &  ,MODEL_MLCS2k2
     &  ,MODEL_SNOOPY
     &  ,MODEL_SALT2
     &  ,MODEL_SIMSED
     &  ,MODEL_BYOSED
     &  ,MODEL_NON1A
     &  ,MODEL_LCLIB    ! Sep 2017
     &  ,MODEL_FIXMAG   ! force mags to user-value
     &  ,MXMODEL_INDEX
c
     &  ,MXCHAR_CCID, MXCHAR_VERSION, MXCHAR_SURVEY
     &  ,MXCHAR_PATH, MXCHAR_FILENAME, MXCHAR_MODELNAME
     &  ,MXCHAR_FIELDNAME, MXCHAR_PARNAME, MXCHAR_CUTNAME
     &  ,MXCHAR_FILEWORD
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA     ! DATA flux per epoch
     &  ,SNLCPAK_EPFLAG_REJECT   ! REJECT flag (1=> excluded from fit)
     &  ,SNLCPAK_EPFLAG_CHI2     ! data-fit chi2 per epoch
     &  ,SNLCPAK_EPFLAG_FITFUN   ! smooth  fitfun curve
     &  ,SNLCPAK_EPFLAG_FLUXSIM  ! SIM flux per epoch
     &  ,SNLCPAK_EPFLAG_FLUXREST ! rest-frame flux per epoch (optional)
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST ! idem for sim truth
     &  ,SNLCPAK_BANDFLAG_PKFLUX   ! peak flux vs. filter
     &  ,SNLCPAK_BANDFLAG_PKMJD    ! peak MJD  vs. filter
     &  ,SNLCPAK_BANDFLAG_NDOF     ! Ndof vs. filter
     &  ,SNLCPAK_BANDFLAG_CHI2     ! chi2 vs. filter
c
     &  ,IFLAG_INI, IFLAG_ANA, IFLAG_END
     &  ,MAG_SATURATE
     &  ,MASK_FLUXCOR_SNANA, MASK_FLUXERRCOR_SNANA
c
     &  ,EXIT_ERRCODE_SNANA, EXIT_ERRCODE_SNFIT, EXIT_ERRCODE_PSNID

      REAL*8
     &   ZERO8, ONE8, TEN8,  PI, LOGTEN
     &  ,PDFMIN, PDFMAX_EDGE, PDFMIN_GOOD
     &  ,KCORPACK
     &  ,XTMW_FRACERR
     &  ,MJDOFF
     &  ,ZEROPOINT_FLUXCAL_DEFAULT
     &  ,RV_MWCOLORLAW_DEFAULT
     &  ,CUTVAL_OPEN, IGNORE_HEADVAL

      REAL NULLVAL

      PARAMETER (                   ! BEGIN SNANA PARAMS
     &   SNTABLE_LIST_DEFAULT = 'SNANA  FITRES  LCPLOT'
     &  ,MXVERS        = 50         ! max number of versions to read
     &  ,MXSURVEY      = 100        ! max number of SURVEY_NAMEs to store
     &  ,MXITER        = 12         ! max # fit iterations
     &  ,NFIT_VERBOSE  = 500        ! Number of fits for verbose printouts
     &  ,MXSNLC        = 4000000    ! max number of SNe (fits format)
     &  ,MXCCID_LIST   = 10000      ! max size of SNCCID_LIST_ALL
     &  ,MXCID       = 299 999 999  ! max CID = 300 million -1 (9 digits)
     &  ,MXCID_CHECK =  99 999 999  ! MXCID to check duplicates
     &  ,MXEPOCH     = 2000        ! max number of filter-epochs per SN
     &  ,MXEP_MODELGRID = 200     ! max model grid for SNANA+SIM_MAGOBS table
     &  ,MXIDTEL     = 200        ! max telescope ID in SURVEY.DEF
     &  ,MXIDFIELD   = 200        ! max FIELD ID in SURVEY.DEF
     &  ,MXFIELD_OVP = 12         ! max number of overlapping fields
     &  ,MXSEASON    = 100        ! max number of seasons
     &  ,MXSNHOST    = 2          ! max number of host matches to read/write
     &  ,MXVAR_PRIVATE = 40       ! max number of private variables
     &  ,MXCUT_PRIVATE = 10        ! max number of cuts on private var
     &  ,MXVAR_TERSE   = 30       ! max number of text-data  columms
     &  ,MXLISTNML   = 52    ! max list size for some NML lists
     &  ,MXFILT_OBS  = 62    ! max number of used observer filters
     &  ,MXFILT_ALL  = 80    ! max number of all possible filter defs
     &  ,MXFILT_REST = 20    ! max number of rest-frame filter types (was 12)
     &  ,MXFILT_KCOR = MXFILT_REST ! max number of filters used in Kcor table
     &  ,MXFILT_SNRMAX = 10       ! no more than this many SNRMAX(filt) cuts
     &  ,MNTYPE       =   1       ! min sn "type"
     &  ,MXTYPE       = 1000      ! max sn "type"
     &  ,MXZBIN_KCOR  = 100       ! max # Z-bins for KCOR tables
     &  ,MXTBIN_KCOR  = 150       ! max # Epochs for KCOR tables
     &  ,MXAVBIN_KCOR = 100       ! max # AV bins for KCOR tables
     &  ,MXCBIN_AVWARP = 100      ! max color index for AVWARP table
     &  ,MXTABLE1D_KCOR   = 8 000 000  ! max number of KCOR bins
     &  ,MXTABLE1D_AVWARP = 1 000 000  ! max number of bins for AVWARP table
     &  ,MXTABLE1D_LCMAG  = 2 000 000  ! idem for LCMAG  table
     &  ,MXTABLE1D_MWXT   = 2 000 000  ! idem for MWXT table
     &  ,KCORPACK      = 1000.  ! store KCOR * KCORPACK as I*2
     &  ,XTMW_FRACERR  = 0.16   ! error on MW Xtinc is 16% of XTMW
c
     &  ,HOFF_NOCUTS = 100
     &  ,HOFF_CUTS   = 200
     &  ,HOFF_SIM    = 20000
     &  ,LUNHIS    = 21  ! LUN for hbook output file
     &  ,LUNNML    = 22  ! LUN for input namelist
     &  ,LUNDAT    = 23  ! LUN for data
     &  ,LUNTMP    = 24
     &  ,LUNOUT    = 25
     &  ,LUNFIT    = 26
     &  ,LUNLIST   = 27
     &  ,LUNIGNORE = 28
     &  ,LUNDMP    = 29
     &  ,LUNRES1   = 31  ! for DMP_FITRES
     &  ,LUNRES2   = 32
     &  ,LUNINTERP = 33  ! for FLUX,MAGs interpolated at SN,MJD
     &  ,LUNSALT2  = 34  ! for SALT2 dictFile
     &  ,LUNPKMJD  = 35
     &  ,LUNCID    = 36  ! reserved for reading SNCID_LIST_FILE
     &  ,ISTAGE_INIT    = 10       ! init has finished
     &  ,ISTAGE_RDSN    = 20       ! SN have been read
     &  ,ISTAGE_CUTS    = 30       ! cuts have been applied
     &  ,ISTAGE_USRANA  = 40       ! USRANA has been called.
     &  ,ISTAGE_TEST    = 90       ! set for testing
     &  ,ZERO8          = 0.0
     &  ,ONE8           = 1.0
     &  ,ONE            = 1
     &  ,TEN8           = 10.0
     &  ,LOGTEN         = 2.302585092994
     &  ,INTERP_LINEAR  = 1    ! option flag to use linear DFINT
     &  ,INTERP_SMOOTH  = 2    ! option flag to use smoothing
     &  ,INTERP_ZSMOOTH = 3    ! linear DFINT, but smooth along z
c
     &  ,MXERRTYPE      = 10
     &  ,ERRTYPE_MINOS  = 1
     &  ,ERRTYPE_PARAB  = 2
     &  ,ERRTYPE_MARG   = 3  ! marginalized error
     &  ,ERRTYPE_BAD    = 6
C
     &  ,FCNFLAG_LAST     = 3    ! last MINUIT call
     &  ,FCNFLAG_USER     = 90   ! pass this IFLAG for user-FCNSNLC calls
     &  ,FCNFLAG_FAST     = 91   ! go as fast as possible (no LAST if-block)
     &  ,FCNFLAG_PRIOR_ONLY = 92
     &  ,FCNFLAG_SIGMA_ONLY = 93
     &  ,FCNFLAG_USESIM   = 99   ! pass this IFLAG to use SIM params
c
     &  ,PI     = 3.1415926535898
     &  ,PDFMIN      = 1.0E-5   ! used to speed up PDF integration
     &  ,PDFMAX_EDGE = 0.03     ! max allowed PDF value at edges
     &  ,PDFMIN_GOOD = 1.0E-4   ! PDF > PDFMIN_GOOD counts as good point
     &  ,OPT_INTEGPDF_QUITCHI2 = 2  ! abort FCNSNLC if chi2 > quitchi2
     &  ,OPT_INTEGPDF_FULL     = 1  ! do full FCNSNLC evaluation
c
     &  ,OPT_SNXT_CCM89   = 1  ! exact SN etinction using INIT_XTHOST
     &  ,OPT_SNXT_SJPAR   = 2  ! SN extinction with Jha's parameters
     &  ,OPT_KCORERR_SMOOTH = 1 ! use smooth half-Gaussian
     &  ,OPT_KCORERR_SJ     = 2   ! use Saurabh's Kcor error
     &  ,OPT_KCORERR_SJ5    = 5  ! x5 Saurabh's Kcor error
c
     &  ,RV_MWCOLORLAW_DEFAULT  = 3.1   ! A_V/E(B-V)
     &  ,OPT_MWCOLORLAW_DEFAULT = 94    ! ODonnel 94
     &  ,OPT_MWEBV_DEFAULT      =  1    ! whatever is in the data file
c
     &  ,MXLINE_ARGS      = 100
     &  ,ZEROPOINT_FLUXCAL_DEFAULT  = 27.5
     &  ,MXEPOCH_IGNORE   = 1000
     &  ,NPAR_ANYLC       = 8    ! for MNFIT_PKMJD
     &  ,MXPAR_SIMSED     = 100  ! max number of SIMSED parameters
     &  ,MXPAR_LCLIB      = 40   ! should be same as in genmag_LCLIB.h
     &  ,MXLAMBIN_SNSED   = 4000 ! Jan 2017: raised from 3000
     &  ,MXCUTBIT         = 64   ! max number of cut bits
c
     &  ,OPT_FILTUPD_EACHSN   = 1  ! default filter-updates
     &  ,OPT_FILTUPD_MAP      = 2  ! default filter-updates
     &  ,OPT_FILTUPD_SAMEFILT = 3  ! test: use same filter each SN
     &  ,OPT_FILTREST         = 1
     &  ,OPT_FILTOBS          = 2
     &  ,ISTAT_READAGAIN      = 7
     &  ,ISTAT_SKIP           = -1
c
     &  ,ITABLE_SNANA=1, ITABLE_FITRES=2, ITABLE_SNLCPAK=3
     &  ,ITABLE_SPECPAK=4, ITABLE_SPECTRA=5
     &  ,MXTABLE = 10
     &  ,IDTABLE_SNANA   = 7100
     &  ,IDTABLE_FITRES  = 7788
     &  ,IDTABLE_SPECTRA = 8000  ! SALT2 model spectra from LC fit
     &  ,IDTABLE_CIDPTR  = 1600
     &  ,IDTABLE_MCMC    = 7711
     &  ,OPT_PARSTORE_TEXTTABLE = 1  ! tag subset for TEXT table.
c
     &  ,MODEL_STRETCH    = 1
     &  ,MODEL_STRETCH2   = 2
     &  ,MODEL_MLCS2k2    = 3
     &  ,MODEL_SNOOPY     = 4
     &  ,MODEL_SALT2      = 6
     &  ,MODEL_SIMSED     = 7
     &  ,MODEL_BYOSED     = 8
     &  ,MODEL_NON1A      = 10
     &  ,MODEL_LCLIB      = 12
     &  ,MODEL_FIXMAG     = 20
     &  ,MXMODEL_INDEX    = 20
     &  ,NULLVAL          = -99999.
     &  ,IGNORE_HEADVAL   = -5555.
c
     &  ,MXCHAR_CCID       = 20   ! max len of CCID string (i.e, SN name)
     &  ,MXCHAR_VERSION    = 72   ! max len of VERSION_PHOTOMETRY
     &  ,MXCHAR_SURVEY     = 40   ! max len of SURVEY_NAME
     &  ,MXCHAR_PATH       = 160  ! max len of path
     &  ,MXCHAR_FILENAME   = 200  ! max len of filename with full path
     &  ,MXCHAR_MODELNAME  = 72   ! max len of model name
     &  ,MXCHAR_FIELDNAME  = 20   ! max len of field name
     &  ,MXCHAR_PARNAME    = 20   ! max len of parameter name
     &  ,MXCHAR_CUTNAME    = 160  ! to define cut names
     &  ,MXCHAR_FILEWORD   =  60  ! size of FILEWORD_LIST
c
     &  ,SNLCPAK_EPFLAG_FLUXDATA    = 1    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_REJECT      = 2
     &  ,SNLCPAK_EPFLAG_CHI2        = 3
     &  ,SNLCPAK_EPFLAG_FITFUN      = 4
     &  ,SNLCPAK_EPFLAG_FLUXSIM     = 5    ! epoch-dependent
     &  ,SNLCPAK_EPFLAG_FLUXREST    = 6
     &  ,SNLCPAK_EPFLAG_SIMFLUXREST = 7
     &  ,SNLCPAK_BANDFLAG_NDOF    = 100
     &  ,SNLCPAK_BANDFLAG_PKFLUX  = 101  ! filter-dependent
     &  ,SNLCPAK_BANDFLAG_PKMJD   = 102
     &  ,SNLCPAK_BANDFLAG_CHI2    = 103
c
     &  ,IFLAG_INI=1, IFLAG_ANA=2, IFLAG_END=3
     &  ,MAG_SATURATE = -7.0   ! for sim only
     &  ,CUTVAL_OPEN = 1.0E12  ! cutwin value to accept everything.
     &	,MASK_FLUXCOR_SNANA    = 1
     &  ,MASK_FLUXERRCOR_SNANA = 2
c
     &  ,EXIT_ERRCODE_SNANA = 21
     &  ,EXIT_ERRCODE_SNFIT = 22
     &  ,EXIT_ERRCODE_PSNID = 23
     &      )

c physical constants

      REAL*8
     &   PARSEC, CLIGHT, PEAKMAG_AT_10PC
     &  ,Zat10pc
     &  ,OMAT_DEFAULT, OMATERR_DEFAULT
     &  ,OLAM_DEFAULT, OLAMERR_DEFAULT
     &  ,ORAD_DEFAULT
     &  ,H0_DEFAULT,   H0ERR_DEFAULT
     &  ,W0_DEFAULT,   W0ERR_DEFAULT
     &  ,DWDA_DEFAULT, DWDAERR_DEFAULT

      PARAMETER (
     &   PARSEC        = 3.085678E13  ! 1 parsec (km)
     &  ,CLIGHT        = 2.998E5      ! c (km/sec)
c
     &  ,H0_DEFAULT        = 70.0 / ( 1.0E6 * Parsec )  ! standard value
     &  ,W0_DEFAULT        = -1.0
     &  ,DWDA_DEFAULT      =  0.0
     &  ,OMAT_DEFAULT      =  0.3
     &  ,OLAM_DEFAULT      =  0.7
     &  ,ORAD_DEFAULT      =  1.2E-5
c
     &  ,H0ERR_DEFAULT        =  7.0 / ( 1.0E6 * Parsec )
     &  ,W0ERR_DEFAULT        =  0.1
     &  ,DWDAERR_DEFAULT      =  0.0
     &  ,OMATERR_DEFAULT      =  0.03
     &  ,OLAMERR_DEFAULT      =  0.03
c
     &  ,PEAKMAG_AT_10PC   = -19.6
     &  ,Zat10pc           = 2.34E-9   ! magic redshift at 10 pc
     &  ,MJDOFF            = 0.0       ! 53000.
     &     )



C +CDE,SNFILECOM. inserted below by ypatchy. 

c Sep 9, 2010: Pulled out of SNDATACOM

c define HEADER MASK BITS for required variables.
      INTEGER
     &   HEADBIT_SNID,   HEADBIT_IAUC, HEADBIT_CIDSEL
     &  ,HEADBIT_SURVEY, HEADBIT_FILTERS
     &  ,HEADBIT_RA,    HEADBIT_DEC
     &  ,HEADBIT_MWEBV, HEADBIT_Z
     &  ,NBIT_HEADMASK
     &  ,HEADMASK_REQUIRED
c
     &  ,HEADMASK  ! reset for reach SN

      PARAMETER (
     &    HEADBIT_SNID    = 0  ! MASK = 1
     &   ,HEADBIT_IAUC    = 1  ! MASK = 2  ! added Dec 2015
     &   ,HEADBIT_CIDSEL  = 2  ! MASK = 4  ! added Dec 2015
     &   ,HEADBIT_SURVEY  = 3
     &   ,HEADBIT_FILTERS = 4
     &   ,HEADBIT_RA      = 5
     &   ,HEADBIT_DEC     = 6
     &   ,HEADBIT_MWEBV   = 7
     &   ,HEADBIT_Z       = 8
     &   ,NBIT_HEADMASK   = 9
     &   ,HEADMASK_REQUIRED  = 2**(NBIT_HEADMASK)-1 - 2 ! don't require IAUC
     &      )

      CHARACTER
     &   SNDATA_ROOT*(MXCHAR_PATH)
     &  ,SNANA_DIR*(MXCHAR_PATH)
     &  ,SNDATA_PATH*(MXCHAR_PATH)      ! subdir with data or sim files
     &  ,SNLIST_FILE*(MXCHAR_FILENAME)  ! input list of SNDATA Files
     &  ,SNREADME_FILE(MXVERS)*(MXCHAR_FILENAME)   ! name of EVERY README file
     &  ,SNDATA_FILE_CURRENT*(MXCHAR_FILENAME)     ! current file being read
     &  ,GLOBAL_BANNER*120
     &  ,SNDATA_PREFIX*(MXCHAR_FILENAME)  ! $SNDATA_ROOT/lcmerge/$VERSION
     &  ,C1ERR*88, C2ERR*88    ! generic error strings

      INTEGER
     &  ABSO_INDEX(MXSNLC) ! absolute (IFILE or IROW) vs. ISN index

      LOGICAL LFLAG_RDHEAD_ONLY

      COMMON / SNFILECOM /
     &   SNDATA_ROOT, SNANA_DIR, SNLIST_FILE
     &  ,SNDATA_FILE_CURRENT
     &  ,GLOBAL_BANNER, SNDATA_PREFIX, SNREADME_FILE, SNDATA_PATH
     &  ,C1ERR, C2ERR, ABSO_INDEX, HEADMASK, LFLAG_RDHEAD_ONLY



C +CDE,CTRLCOM. inserted below by ypatchy. 

c control variables and counters.

      INTEGER
     &   ISTAGE_SNANA       ! current stage of processing
     &  ,NACCEPT_CUT(MXCUTBIT) ! Number of SN that pass each cut
     &  ,NACCEPT_CID           ! # SN with valid CID
     &  ,NACCEPT_TYPE          ! # SN with valid TYPE
     &  ,NACCEPT_Z             ! idem with valid redshift
     &  ,NACCEPT_ZERR          ! idem with valid redshift error
     &  ,NMJD_IDTEL            ! # MJDs with valid telescope id
     &  ,APPLY_HEADER_CUTMASK  ! mask of applied header cuts
c
     &  ,JTIME_START
     &  ,JTIME_LOOPSTART        ! time at start of fits with TIME()
     &  ,JTIME_LOOPEND          ! time and end
     &  ,NCALL_SNANA_DRIVER
c
     &  ,NPASSCUT_INCREMENT(-1:MXTYPE,100)  ! 100 > NCUTBIT_SNLC
     &  ,NPASSCUT_FIT(-1:MXTYPE)
     &  ,NVAR_NEARNBR          ! Number of NEARNBR variables to analyze
     &  ,NSTORE_MAGCOR         ! number of stored MAGCOR values
     &  ,NUSE_MAGCOR           ! number of used MAGCOR values
     &  ,SIGN_MAGCOR           ! add or subtract
     &  ,FORCEMASK_FLUXCOR   ! mask to force fluxCor, even if already applied
     &  ,EXIT_ERRCODE        ! used for abort

      LOGICAL
     &   DO_FIT
     &  ,DO_FLUXERRCALC    ! T => compute error from PSF,SKY & ZPT
     &  ,LSIM_SNANA        ! simulated with SNANA
     &  ,LSIM_MAGOBS       ! data-like, but with SIM_MAGOBS
     &  ,ISJOB_SNANA          ! =T for snana.exe only.
     &  ,ISJOB_PSNID
     &  ,REFORMAT_SAVE_BADEPOCHS  ! set if bit2 is set on any OPT_REFORMAT
     &  ,STDOUT_UPDATE         ! T => update event to screen
     &  ,DOFUDGE_HOSTNOISE     ! T => FUDGE_HOSTNOISE_FILE is set
     &  ,DOFUDGE_NONLIN        ! T => NONLINEARITY_FILE is set
     &  ,DOFUDGE_FLUXERRMODEL  ! T => FLUXERRMODEL_FILE

      CHARACTER SNANA_VERSION*12

c global survey info
      CHARACTER
     &   SURVEY_NAME*(MXCHAR_SURVEY)
     &  ,SURVEY_NAME_LIST(MXSURVEY)*(MXCHAR_SURVEY) ! in SURVEY.DEF file
     &  ,SUBSURVEY_NAME*(MXCHAR_SURVEY)  ! e.g.,  SURVEY:  BLA(SUBSURVEY)
     &  ,SUBSURVEY_NAME_LIST*(MXCHAR_FILENAME) ! comma-sep list
     &  ,SURVEY_FILTERS*(MXFILT_ALL)  ! read from data file
     &  ,SURVEY_FIELDNAME(MXIDFIELD)*(MXCHAR_FIELDNAME) ! from SURVEY.DEF

      INTEGER
     &   NFIELD_SURVEY             ! number of survey fields in SURVEY.DEF
     &  ,SURVEY_IDFIELD(MXIDFIELD) ! integer ID for each field

      REAL
     &   ZEROPOINT_FLUXCAL(MXFILT_OBS)  ! defines calibrated flux

      INTEGER*4
     &   N_VERSION             ! number of photometry version to read
     &  ,N_SNLC                ! number of SN lightcurves read
     &  ,N_SNLC_CUTS           ! Number of SN after cuts (bookkeeping only)
     &  ,N_SNLC_FIT            ! Number of fitted SN
     &  ,N_SNLC_FITCUTS        ! Number of fitted SN after fit cuts
     &  ,N_SNLC_COVFIX         ! Number of SN with fixed COV to be invertible
     &  ,N_DUPLICATE_CID       ! Number of duplicate CIDs
     &  ,N_DUPLICATE_MJD       ! Number of duplicate MJD+BAND (Jun 2017)
     &  ,NSTORE_DUPLICATE_MJD  ! Number stored
     &  ,N_SNFILE              ! # of SN files to read per version
     &  ,N_SNFILE_LAST         ! idem, as of last version
     &  ,NEPOCH_TOT            ! total number of epochs read
     &  ,NEPOCH_CUT            ! total number of epochs passing cuts
     &  ,NEPOCH_BADPHOT        ! # epochs with bad PHOTFLAG (per event)
     &  ,NEPOCH_BADPHOT_SUM    ! # epochs with bad PHOTFLAG (summed)
     &  ,IDSURVEY              ! survey ID from SURVEY.DEF
     &  ,IDSUBSURVEY           ! =IDSURVEY unless subSurvey is different
     &  ,IDSURVEY_LIST(MXSURVEY) ! corresponds to SURVEY_NAME_LIST
     &  ,NSURVEY_LIST          ! size of SURVEY_NAME_LIST & IDSURVEY_LIST

      LOGICAL*1
     &    EXIST_FILT(MXFILT_OBS)  ! T => at least one point per filt
     &   ,FOUND_SURVEY
     &   ,FORMAT_TEXT    ! ascii/txt for input data
     &   ,FORMAT_TERSE   ! terse-text for input data
     &   ,FORMAT_VERBOSE ! verbose-text for ...
     &   ,FORMAT_FITS    ! snfitsio for input data.

      INTEGER
     &   NVAR_TERSE          ! Number variables for terse format
     &  ,ITERSE_MJD          ! index of MJD in TERSE variable list
     &  ,ITERSE_FILTER       ! idem for filter
     &  ,ITERSE_FIELD
     &  ,ITERSE_FLUXCAL
     &  ,ITERSE_FLUXCALERR
     &  ,ITERSE_PHOTFLAG
     &  ,ITERSE_PHOTPROB
     &  ,ITERSE_ZPFLUX     ! mag <-> native flux zeropoint
     &  ,ITERSE_PSFSIG     ! added July 18, 2013
     &  ,ITERSE_SKYSIG     ! idem
     &  ,ITERSE_SKYSIG_T   ! added Aug 7 2014
     &  ,ITERSE_XPIX       ! added Aug 7 2014
     &  ,ITERSE_YPIX       ! added Aug 7 2014
     &  ,ITERSE_GAIN       ! added Oct 2015 to compute FLUXCAL_ERRCALC
     &  ,ITERSE_MAGOBS
     &  ,ITERSE_CCDNUM     ! added Oct 16 2017
     &  ,ITERSE_SIM_EPFILTREST
     &  ,ITERSE_SIM_EPMAGOBS
     &  ,ITERSE_SIM_EPMAGREST

      INTEGER N_SNLC_PLOT
      LOGICAL MADE_LCPLOT  ! SAVE:  T if LC plot was made

      CHARACTER VARLIST_TERSE(MXLISTNML)*(MXCHAR_CCID) ! TERSE variable names

      COMMON / CTRLCOM /
     &     SNANA_VERSION, SURVEY_NAME, SURVEY_NAME_LIST
     &    ,SUBSURVEY_NAME, SUBSURVEY_NAME_LIST
     &    ,NSURVEY_LIST, IDSURVEY, IDSUBSURVEY, IDSURVEY_LIST
     &    ,SURVEY_FILTERS
     &    ,SURVEY_FIELDNAME, SURVEY_IDFIELD, NFIELD_SURVEY
     &    ,ISJOB_SNANA, ISJOB_PSNID, ZEROPOINT_FLUXCAL
     &    ,NACCEPT_CUT, NACCEPT_CID, NACCEPT_TYPE
     &    ,NACCEPT_Z, NACCEPT_ZERR, NMJD_IDTEL
     &    ,DO_FIT, DO_FLUXERRCALC, NVAR_NEARNBR
     &    ,NSTORE_MAGCOR, NUSE_MAGCOR, SIGN_MAGCOR, FORCEMASK_FLUXCOR
     &    ,LSIM_SNANA, LSIM_MAGOBS, APPLY_HEADER_CUTMASK
     &    ,ISTAGE_SNANA
     &    ,N_VERSION, N_SNLC, N_SNLC_CUTS, N_SNLC_FIT, N_SNLC_FITCUTS
     &    ,N_SNLC_COVFIX, N_SNFILE, N_SNFILE_LAST, N_DUPLICATE_CID
     &    ,N_DUPLICATE_MJD, NSTORE_DUPLICATE_MJD
     &    ,NEPOCH_TOT, NEPOCH_CUT, NEPOCH_BADPHOT, NEPOCH_BADPHOT_SUM
     &    ,REFORMAT_SAVE_BADEPOCHS, STDOUT_UPDATE
     &    ,DOFUDGE_HOSTNOISE, DOFUDGE_NONLIN, DOFUDGE_FLUXERRMODEL
     &    ,NVAR_TERSE, VARLIST_TERSE
     &    ,ITERSE_MJD, ITERSE_FILTER, ITERSE_FIELD
     &    ,ITERSE_FLUXCAL, ITERSE_FLUXCALERR
     &    ,ITERSE_PHOTFLAG, ITERSE_PHOTPROB
     &    ,ITERSE_ZPFLUX, ITERSE_PSFSIG, ITERSE_SKYSIG, ITERSE_SKYSIG_T
     &    ,ITERSE_MAGOBS, ITERSE_XPIX, ITERSE_YPIX, ITERSE_GAIN
     &    ,ITERSE_CCDNUM, ITERSE_SIM_EPFILTREST
     &    ,ITERSE_SIM_EPMAGOBS, ITERSE_SIM_EPMAGREST
     &    ,JTIME_START, JTIME_LOOPSTART, JTIME_LOOPEND
     &    ,NCALL_SNANA_DRIVER, NPASSCUT_INCREMENT, NPASSCUT_FIT
     &    ,N_SNLC_PLOT, MADE_LCPLOT
     &    ,EXIT_ERRCODE

c logical *1 stuff

      COMMON / SNDATCOM1 / EXIST_FILT, FOUND_SURVEY
     &    ,FORMAT_TEXT, FORMAT_TERSE, FORMAT_VERBOSE
     &    ,FORMAT_FITS

c SNTABLE control variables (May 2014)
      LOGICAL
     &   USE_TABLEFILE_HBOOK  ! write table in HBOOK format
     &  ,USE_TABLEFILE_ROOT   ! write table in ROOT format
     &  ,USE_TABLEFILE_TEXT   ! write table in TEXT format
     &  ,USE_TABLEFILE        ! T if either any of the above are set
     &  ,WRTABLEFILE_IAUC     ! T-> SNID -> IAUC (for writing tables)
     &  ,WRTABLEFILE_SIMVAR   ! T-> include SIM_XXX vars for simulation
     &  ,WRTABLEFILE_ZPHOT    ! T-> include ZPHOT info for PHOTOZ fit


      INTEGER
     &   OPT_TABLE(MXTABLE)       ! option(s) for each table
     &  ,CUTMASK_SNANA_TABLE      ! select CUTFLAG_SNANA for SNANA table

      REAL*8
     &   PRESCALE_TABLE(MXTABLE)  ! table pre-scale (to reduce size)

      CHARACTER
     &   TEXTFORMAT_TABLE(MXTABLE)*8

c arrays for TABLE_FILTER_REMAP (Feb 2017)
      INTEGER NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE(MXFILT_ALL)
      CHARACTER FILTLIST_REMAP_TABLE*(MXFILT_ALL)

c Sep 23 2017: allow user codes to add private variables to SNTABLE
      INTEGER   NTABLEVAR_USER  ! e.g., set in snana_private.cra
      CHARACTER TABLEVARNAME_USER(MXVAR_PRIVATE)*40
      REAL      TABLEVALUE_USER(MXVAR_PRIVATE)

      COMMON / SNTABLECOM /
     &     USE_TABLEFILE_HBOOK, USE_TABLEFILE_ROOT
     &    ,USE_TABLEFILE_TEXT, USE_TABLEFILE
     &    ,OPT_TABLE, TEXTFORMAT_TABLE
     &    ,CUTMASK_SNANA_TABLE, WRTABLEFILE_IAUC, WRTABLEFILE_SIMVAR
     &    ,WRTABLEFILE_ZPHOT
c
     &    ,NFILT_REMAP_TABLE, IFILTLIST_REMAP_TABLE
     &    ,FILTLIST_REMAP_TABLE
     &    ,NTABLEVAR_USER, TABLEVARNAME_USER
     &    ,TABLEVALUE_USER

      COMMON / SNTABLECOM8 / PRESCALE_TABLE

c May 6, 2008: define lists for epochs to ignore

      INTEGER NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
      CHARACTER
     &   EPOCH_IGNORE_CCID(MXEPOCH_IGNORE)*(MXCHAR_CCID)
     &  ,EPOCH_IGNORE_FILT(MXEPOCH_IGNORE)*4
     &  ,EPOCH_IGNORE_LASTFILE*(MXCHAR_FILENAME)

      REAL*8
     &   EPOCH_IGNORE_MJD(MXEPOCH_IGNORE)

      COMMON / EPIGNORE_COM / NEPOCH_IGNORE, NEPOCH_IGNORE_WRFITS
     &  ,EPOCH_IGNORE_CCID, EPOCH_IGNORE_FILT, EPOCH_IGNORE_LASTFILE
      COMMON / EPIGNORE_COM8 / EPOCH_IGNORE_MJD


      REAL*8  DUPLICATE_MJDLIST(200)
      COMMON / DUPMJDCOM / DUPLICATE_MJDLIST


C +CDE,SNLCCOM. inserted below by ypatchy. 

c Main common block variables for light curves and
c analysis-related variables.

      INTEGER
     &   SNLC_CID               ! integer CID
     &  ,ISNLC_SNRECON_USE(MXEPOCH)  ! 1 -> epoch used in SNRECON
     &  ,ISNLC_PHOTFLAG(MXEPOCH) ! photomety flags
     &  ,ISNLC_LENCCID     ! char-len of CCID
     &  ,ISNLC_LENIAUC     ! char-len of IAUC name
     &  ,ISNLC_VERSION     ! photometry version index vs. ISN
     &  ,ISNLC_TYPE        ! type (120=confirmed Ia, etc ...)
     &  ,ISNLC_ABSO_INDEX  ! absolute index (row or file number)
     &  ,ISNLC_IFILE       ! ifile index, text format only
     &  ,ISNLC_IFILT_OBS(MXEPOCH) ! filt-index for each epoch/SN
     &  ,ISNLC_NEWMJD_HEAD   ! Number of NEWMJDs in header
     &  ,ISNLC_NEWMJD_FOUND  ! Number of NEWMJDs found in file
     &  ,ISNLC_NEWMJD_STORE  ! Number of NEWMJDs stored
     &  ,ISNLC_NEWMJD_CUTS   ! Number of NEWMJDs after cuts
     &  ,ISNLC_EPOCH_RANGE_NEWMJD(2,MXEPOCH)
     &  ,ISNLC_NFILT_NEWMJD(MXEPOCH)   ! # observed filters per NEWMJD
     &  ,ISNLC_NFILT_SNRMAX      ! number of filters passing SNR cut
     &  ,ISNLC_NFILT_SNRMAX2     ! idem for 2nd SNRMAX cut
     &  ,ISNLC_NFILT_TRESTMIN     ! Nfilt passing TRESTMIN cut
     &  ,ISNLC_NFILT_TRESTMAX     ! idem for TRESTMAX
     &  ,ISNLC_NFILT_TREST2       ! Nfilt passing CUTWIN_TREST2 cut
     &  ,ISNLC_NFILT_THRESH(MXEPOCH)
     &  ,ISNLC_NEPOCH_FOUND  ! actual number of epochs found
     &  ,ISNLC_NEPOCH_STORE  ! number of epochs stored in memory
     &  ,ISNLC_NEPOCH_USE       ! used after PHOTMASK cuts
     &  ,ISNLC_NEPOCH_PHOTPROB  ! NEPOCH with PHOTPROB >= 0
     &  ,ISNLC_NEPOCH_DETECT    ! NEPOCH with detection
     &  ,ISNLC_NEPOCH_FILT(MXFILT_OBS)  ! NEPOCH vs. filter
     &  ,ISNLC_NEPOCH_PRESN(MXFILT_OBS) ! pre-explosion epochs
     &  ,ISNLC_NMJD_INCLUDE             ! NOBS in CUTWIN_MJD_INCLUDE window
c
     &  ,ISNLC_FAKE             ! => real data, else it's a fake
     &  ,ISNLC_CCDNUM(MXEPOCH)  ! read from header (May 2017)
     &  ,ISNLC_IDTEL(MXEPOCH)   ! integer telescope id
     &  ,ISNLC_IDFIELD(MXEPOCH) ! integer field id
     &  ,ISNLC_NFIELD_OVP       ! number of fields (>=2 for overlap)
     &  ,ISNLC_CUTFLAG_REQEP        ! idem
     &  ,ISNLC_CUTFLAG_PRIVATE      ! idem for private var cuts
     &  ,ISNLC_CUTFLAG_USRCUTS      ! idem for USRCUTS
     &  ,ISNLC_CUTFLAG_SIMVAR       ! idem for SIMVAR cuts
     &  ,ISNLC_WRMASK_FLUXCOR_SNANA   ! write if fudges applied to data
     &  ,ISNLC_RDMASK_FLUXCOR_SNANA   ! read if SNANA fudges applied to data

      REAL*8
     &   SNLC8_RA        ! RA   vs. SN
     &  ,SNLC8_DEC       ! DEC  vs. SN
     &  ,SNLC8_MJD(MXEPOCH)  ! MJD for each epoch
     &  ,SNLC8_MJDMIN        !   min MJD among all measurements
     &  ,SNLC8_MJDMAX        !   max MJD ...

      REAL
     &   SNLC_ZHELIO        ! redshift used for final analysis
     &  ,SNLC_ZHELIO_ERR    ! error on above
     &  ,SNLC_ZCMB          ! redshift used for final analysis
     &  ,SNLC_ZCMB_ERR      ! error on above
     &  ,SNLC_ZHD           ! redshift used for Hubble diagram
     &  ,SNLC_ZHD_ERR       ! error on above
     &  ,SNLC_ZSN           ! redshift of SN only (ignoring host-z)
     &  ,SNLC_ZSN_ERR       ! error on above
     &  ,SNLC_REDSHIFT      ! redshift used for light curve fit
     &  ,SNLC_REDSHIFT_ERR  ! error on above
     &  ,SNLC_VPEC          ! pec. velocity
     &  ,SNLC_VPEC_ERR      ! error on above
     &  ,SNLC_ZPEC
     &  ,SNLC_ZPEC_ERR
     &  ,SNLC_Trestmin  ! earliest epoch, rest frame days since peak
     &  ,SNLC_Trestmax  ! latest   epoch, rest frame days since peak
     &  ,SNLC_TrestRange
     &  ,SNLC_Tobsmin   ! earliest epoch, obs frame days since peak
     &  ,SNLC_Tobsmax   ! latest   epoch, obs frame days since peak
     &  ,SNLC_TGAPMAX     ! max gap within TREST-range
     &  ,SNLC_T0GAPMAX    ! max gap near peak
     &  ,SNLC_SNRMAX_FILT(0:MXFILT_OBS)   ! max S/N per filter/SN
     &  ,SNLC_SNRMAX_SORT(MXFILT_OBS)     ! 1st, 2nd ... SNRMAX by filt
     &  ,SNLC_FLUXCALMAX(MXFILT_OBS)    ! max flux per filter/SN
     &  ,SNLC_SNANACALC_PEAKMJD     ! SNANA-estimate of PEAKKMJD
     &  ,SNLC_SNANACALC_PEAKMJD_FITPAR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_SNANACALC_PEAKMJD_FITERR(MXFILT_OBS,NPAR_ANYLC)
     &  ,SNLC_PHOTPROB(MXEPOCH)   ! generic 'fit probability' per epoch
     &  ,SNLC_PHOTPROB_MIN        ! min photprob for PHOTPROB>0
     &  ,SNLC_AIRMASS(MXEPOCH)
     &  ,SNLC_CLOUDCAM_AVG(MXEPOCH)
     &  ,SNLC_CLOUDCAM_SIG(MXEPOCH)
     &  ,SNLC_MOONDIST(MXEPOCH)    ! in deg.
     &  ,SNLC_MOONPHASE(MXEPOCH)   ! 0 - 1
     &  ,SNLC_TREST(MXEPOCH)          !MJD-SearhMJDPEAK / (1+z)
     &  ,SNLC_GAIN(MXEPOCH)      ! e/AUD
     &  ,SNLC_RDNOISE(MXEPOCH)   ! read noise per pix, e-
     &  ,SNLC_PIXSIZE            ! pixel size
     &  ,SNLC_NXPIX              ! total number of X-pixels (Aug 7 2014)
     &  ,SNLC_NYPIX              ! total number of Y-pixels
     &  ,SNLC_XPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_YPIX(MXEPOCH)      ! pixel location
     &  ,SNLC_AREAFRAC(MXEPOCH)  ! area-frac contained by XPIX,YPIX
     &  ,SNLC_MWEBV              ! Milky Way Galactic E(B-V)
     &  ,SNLC_MWEBV_ERR          ! error on above
     &  ,SNLC_SKY_AVG(MXEPOCH)   ! avg sky/pix in field, ADU
     &  ,SNLC_SKYSIG(MXEPOCH)    ! sigma on above
     &  ,SNLC_SKYSIG_T(MXEPOCH)  ! sigma on template run
     &  ,SNLC_PSF_SIG1(MXEPOCH)  ! sigma, pixels
     &  ,SNLC_PSF_SIG2(MXEPOCH)
     &  ,SNLC_PSF_RATIO(MXEPOCH)
     &  ,SNLC_AREA_NOISE(MXEPOCH)  ! Noise-equivalent area (pixels)
     &  ,SNLC_FLUX(MXEPOCH)
     &  ,SNLC_FLUX_ERRTOT(MXEPOCH)
     &  ,SNLC_FLUX_NSIG(MXEPOCH)   ! significance
     &  ,SNLC_FLUXCAL_OFF(MXFILT_OBS)    ! add SN light from template
     &  ,SNLC_FLUXCAL_ERRCALC(MXEPOCH)   ! calc a-la simulation
     &  ,SNLC_FLUXCAL_HOSTERRCALC(MXEPOCH)   ! idem for host error
     &  ,SNLC_FLUXCAL(MXEPOCH)
     &  ,SNLC_FLUXCAL_ERRTOT(MXEPOCH)
     &  ,SNLC_MAG(MXEPOCH)
     &  ,SNLC_MAG_ERRPLUS(MXEPOCH)
     &  ,SNLC_MAG_ERRMINUS(MXEPOCH)
     &  ,SNLC_ZEROPT(MXEPOCH)
     &  ,SNLC_ZEROPT_ERR(MXEPOCH)
     &  ,SNLC_ZEROPT_forCUT(MXEPOCH)  ! depends on CUTWIN_ZPADU or CUTWIN_ZPNPE
     &  ,SNLC_DLMAG                   ! 5*log10(10pc/DL)
     &  ,SNLC_SKYFLUXCAL(MXEPOCH)     ! calculated sky fluxcal/pixel
     &  ,SNLC_MWXT_MAG(MXFILT_OBS)        ! mag stellar extinct
     &  ,SNLC_MWXT_FLUXFRAC(MXFILT_OBS)   ! same for flux (< 1)
     &  ,SNLC_MWXT_MAGERR(MXFILT_OBS)     ! Galactic mag err per filter
     &  ,SNLC_SEARCH_PEAKMJD            ! external PEAKMJD
     &  ,SNLC_DTOBS(MXEPOCH)            ! time since last obs
     &  ,SNLC_DTOBS_SAMEFILT(MXEPOCH)   ! idem, but same filter
     &  ,SNLC_TLIVE_DETECT     ! MJD(last detection) - MJD(1st detection)

      INTEGER
     &   NSEASON_TOT     ! total number of seasons
     &  ,NSEASON_ACTIVE  ! total number of active seasons

      REAL*4
     &   MULTISEASON_CHI2RED(MXSEASON) ! borrow another 'MX' param
     &  ,MULTISEASON_AVGFLUX(MXSEASON)
     &  ,MULTISEASON_MJDMIN(MXSEASON)
     &  ,MULTISEASON_MJDMAX(MXSEASON)

      CHARACTER
     &   SNLC_CCID*(MXCHAR_CCID)      ! char CID
     &  ,SNLC_IAUC*(MXCHAR_CCID)
     &  ,SNLC_CASTCID*8             ! = 'INT' or 'CHAR' to indicate cast
     &  ,SNLC_FIELD(MXEPOCH)*(MXCHAR_FIELDNAME)    ! char field name
     &  ,SNLC_FIELD_OVPLIST(MXFIELD_OVP)*(MXCHAR_FIELDNAME)  ! sparse list of overlap fields
     &  ,SNLC_FIELDLIST*(MXCHAR_FIELDNAME) ! e.g., '82S', '82N+82S'

      COMMON / SNLC8COM /
     &    SNLC8_RA, SNLC8_DEC
     &  , SNLC8_MJD, SNLC8_MJDMIN, SNLC8_MJDMAX

      COMMON / SNLCCOM /
     &     SNLC_CID, SNLC_CCID,  SNLC_CASTCID, SNLC_IAUC
     &    ,SNLC_FIELD, SNLC_FIELD_OVPLIST, SNLC_FIELDLIST
     &    ,SNLC_SNANACALC_PEAKMJD
     &    ,SNLC_SNANACALC_PEAKMJD_FITPAR
     &    ,SNLC_SNANACALC_PEAKMJD_FITERR
     &    ,SNLC_SEARCH_PEAKMJD
     &    ,SNLC_ZHELIO, SNLC_ZHELIO_ERR
     &    ,SNLC_ZCMB,   SNLC_ZCMB_ERR, SNLC_ZHD, SNLC_ZHD_ERR
     &    ,SNLC_ZSN, SNLC_ZSN_ERR
     &    ,SNLC_REDSHIFT, SNLC_REDSHIFT_ERR
     &    ,SNLC_VPEC, SNLC_VPEC_ERR, SNLC_ZPEC, SNLC_ZPEC_ERR
     &    ,SNLC_TRESTMIN, SNLC_TRESTMAX, SNLC_TrestRange, SNLC_TREST
     &    ,SNLC_TGAPMAX,  SNLC_T0GAPMAX
     &    ,SNLC_Tobsmin, SNLC_Tobsmax
     &    ,SNLC_SNRMAX_FILT, SNLC_SNRMAX_SORT
     &    ,SNLC_FLUXCALMAX, SNLC_PHOTPROB, SNLC_PHOTPROB_MIN
     &    ,SNLC_CLOUDCAM_AVG, SNLC_CLOUDCAM_SIG
     &    ,SNLC_MOONDIST, SNLC_MOONPHASE, SNLC_AIRMASS
     &    ,SNLC_GAIN, SNLC_RDNOISE
     &    ,SNLC_MWXT_MAG, SNLC_MWXT_FLUXFRAC, SNLC_MWXT_MAGERR
     &    ,SNLC_MWEBV, SNLC_MWEBV_ERR
     &    ,SNLC_XPIX, SNLC_YPIX, SNLC_AREAFRAC, SNLC_PIXSIZE
     &    ,SNLC_NXPIX, SNLC_NYPIX
     &    ,SNLC_SKY_AVG,  SNLC_SKYSIG,  SNLC_SKYSIG_T
     &    ,SNLC_PSF_SIG1, SNLC_PSF_SIG2, SNLC_PSF_RATIO
     &    ,SNLC_AREA_NOISE
c
     &    ,SNLC_FLUX, SNLC_FLUX_ERRTOT
     &    ,SNLC_FLUX_NSIG, SNLC_FLUXCAL_OFF
     &    ,SNLC_FLUXCAL_ERRCALC, SNLC_FLUXCAL_HOSTERRCALC
     &    ,SNLC_FLUXCAL, SNLC_FLUXCAL_ERRTOT
     &    ,SNLC_MAG, SNLC_MAG_ERRPLUS, SNLC_MAG_ERRMINUS
     &    ,SNLC_ZEROPT, SNLC_ZEROPT_ERR, SNLC_ZEROPT_forCUT
     &    ,SNLC_SKYFLUXCAL, SNLC_DLMAG
     &    ,SNLC_DTOBS, SNLC_DTOBS_SAMEFILT

      COMMON / ISNCOM /
     &     ISNLC_VERSION, ISNLC_SNRECON_USE, ISNLC_PHOTFLAG
     &    ,ISNLC_NFILT_THRESH, ISNLC_TYPE
     &    ,ISNLC_ABSO_INDEX, ISNLC_IFILE,ISNLC_LENCCID,ISNLC_LENIAUC
     &    ,ISNLC_NFILT_NEWMJD, ISNLC_IFILT_OBS
     &    ,ISNLC_NEWMJD_HEAD, ISNLC_NEWMJD_FOUND, ISNLC_NEWMJD_STORE
     &    ,ISNLC_NEWMJD_CUTS, ISNLC_EPOCH_RANGE_NEWMJD
     &    ,ISNLC_NEPOCH_FOUND, ISNLC_NEPOCH_STORE
     &    ,ISNLC_NEPOCH_USE, ISNLC_NEPOCH_PHOTPROB
     &    ,ISNLC_NEPOCH_FILT,  ISNLC_NEPOCH_PRESN, ISNLC_NMJD_INCLUDE
     &    ,ISNLC_NEPOCH_DETECT, SNLC_TLIVE_DETECT
     &    ,ISNLC_FAKE, ISNLC_CCDNUM
     &    ,ISNLC_NFILT_SNRMAX, ISNLC_NFILT_SNRMAX2
     &    ,ISNLC_NFILT_TRESTMIN, ISNLC_NFILT_TRESTMAX
     &    ,ISNLC_NFILT_TREST2
     &    ,ISNLC_IDTEL, ISNLC_IDFIELD, ISNLC_NFIELD_OVP
     &    ,ISNLC_CUTFLAG_REQEP
     &    ,ISNLC_CUTFLAG_PRIVATE, ISNLC_CUTFLAG_USRCUTS
     &    ,ISNLC_CUTFLAG_SIMVAR
     &    ,ISNLC_WRMASK_FLUXCOR_SNANA, ISNLC_RDMASK_FLUXCOR_SNANA

c - - - - -  ADDCOL stuff - - - - - - -
c ADDCOL arrays are loaded with original filter indices,
c or with REMAPed filter indices. Allows using unique
c pointer for calls to SNTABLE_ADDCOL.

      CHARACTER
     &    ADDCOL_FILTERS*(MXFILT_ALL) ! SURVEY_FILTERS

      REAL
     &    ADDCOL_SNHOST_MAGOBS(MXFILT_ALL,MXSNHOST)
     &   ,ADDCOL_SNHOST_SBFLUXCAL(MXFILT_ALL)
     &   ,ADDCOL_FLUXCALMAX(MXFILT_ALL)
     &   ,ADDCOL_CHI2_FITPKMJD(MXFILT_ALL)
     &   ,ADDCOL_SNRMAX(MXFILT_ALL)
     &   ,ADDCOL_XTMW(MXFILT_ALL)

      INTEGER
     &    ADDCOL_NDOF_FITPKMJD(MXFILT_ALL)

      COMMON / ADDCOL_FILTERCOM /
     &   ADDCOL_FILTERS
     &  ,ADDCOL_SNHOST_MAGOBS, ADDCOL_SNHOST_SBFLUXCAL
     &  ,ADDCOL_FLUXCALMAX, ADDCOL_CHI2_FITPKMJD, ADDCOL_NDOF_FITPKMJD
     &  ,ADDCOL_SNRMAX, ADDCOL_XTMW


c - - - -- MULTI-SEASON STUFF - - - -
      COMMON / MULTISEASONCOM / NSEASON_TOT, NSEASON_ACTIVE
     &  ,MULTISEASON_CHI2RED, MULTISEASON_AVGFLUX
     &  ,MULTISEASON_MJDMIN,  MULTISEASON_MJDMAX


C +CDE,SNCUTS. inserted below by ypatchy. 

c Define cut-window for each variables
c and cut-mask for each epoch.
c User must fill cutwin_XXX arrays in main routine.

      INTEGER
     &   CUTBIT_CID, CUTBIT_SNTYPE
     &  ,CUTBIT_REDSHIFT, CUTBIT_REDSHIFT_ERR
     &  ,CUTBIT_RA,  CUTBIT_DEC
     &  ,CUTBIT_HOSTSEP
     &  ,CUTBIT_Trestmin, CUTBIT_Trestmax, CUTBIT_TrestRange
     &  ,CUTBIT_TREST2
     &  ,CUTBIT_Tgapmax,  CUTBIT_T0gapmax
     &  ,CUTBIT_Tobsmin,  CUTBIT_Tobsmax
     &  ,CUTBIT_SEARCH           ! sim only
     &  ,CUTBIT_TEMPLATE_INDEX   ! sim only (Dec 2018)
     &  ,CUTBIT_PEAKMJD
     &  ,CUTBIT_NMJD_INCLUDE    ! Aug 11 2015
     &  ,CUTBIT_NEPOCH
     &  ,CUTBIT_SNRMAX         ! global SNRMAX cut; any filter(s)
     &  ,CUTBIT_SNRMAX2        ! 2nd global SNRMAX cut; any filter(s)
     &  ,CUTBIT_NFIELD         ! cut on number of fields used by SN
     &  ,CUTBIT_MWEBV          ! cut on Galactic extinct (May 8 2012)
     &  ,CUTBIT_NSEASON_ACTIVE ! Nseason with activity
     &  ,CUTBIT_REQEP          ! required epochs (to select short events)
     &  ,CUTBIT_PRIVATE        ! all private-var cuts together (Nov 3 2014)
     &  ,CUTBIT_SIMVAR         ! SIMVAR cuts
     &  ,CUTBIT_USRCUTS        ! user cuts from private code (Aug 2015)
     &  ,CUTBIT_OFFSET_SBFLUX  ! host surface brightness
     &  ,CUTBIT_OFFSET_SNRMAX  ! SNRMAX cuts in each filter
     &  ,CUTBIT_MJD_MARKER     ! things above do NOT depend on epoch
     &  ,CUTBIT_PSF,  CUTBIT_ZP
     &  ,CUTBIT_MOONDIST,  CUTBIT_PHOTPROB, CUTBIT_AIRMASS
     &  ,CUTBIT_TREST, CUTBIT_MJD
     &  ,CUTBIT_IDTEL          ! valid telescope
     &  ,CUTBIT_NBAND_THRESH
     &  ,CUTBIT_NFILT_SNRMAX
     &  ,CUTBIT_NFILT_SNRMAX2
     &  ,CUTBIT_NFILT_TRESTMIN
     &  ,CUTBIT_NFILT_TRESTMAX
     &  ,CUTBIT_NFILT_TREST2
     &  ,CUTBIT_OFFSET
     &  ,NCUTBIT
     &  ,NCUTBIT_SNLC

      PARAMETER (
     &   CUTBIT_CID            = 1
     &  ,CUTBIT_SNTYPE         = 2
     &  ,CUTBIT_REDSHIFT       = 3
     &  ,CUTBIT_REDSHIFT_ERR   = 4
     &  ,CUTBIT_RA             = 5
     &  ,CUTBIT_DEC            = 6
     &  ,CUTBIT_HOSTSEP        = 7
     &  ,CUTBIT_TrestMIN       = 8
     &  ,CUTBIT_TrestMAX       = 9
     &  ,CUTBIT_TrestRange     = 10  ! Dec 2017 (TrestMax-TrestMin)
     &  ,CUTBIT_Trest2         = 11  ! Oct 2012
     &  ,CUTBIT_Tgapmax        = 12
     &  ,CUTBIT_T0gapMAX       = 13
     &  ,CUTBIT_TobsMIN        = 14
     &  ,CUTBIT_TobsMAX        = 15
     &  ,CUTBIT_PEAKMJD        = 16
     &  ,CUTBIT_NMJD_INCLUDE   = 17
     &  ,CUTBIT_NEPOCH         = 18  ! total # of measurements
     &  ,CUTBIT_SEARCH         = 19  ! found by search (SIM only)
     &  ,CUTBIT_NFILT_SNRMAX   = 20
     &  ,CUTBIT_NFILT_SNRMAX2  = 21
     &  ,CUTBIT_NFILT_TRESTMIN = 22
     &  ,CUTBIT_NFILT_TRESTMAX = 23
     &  ,CUTBIT_NFILT_TREST2   = 24
     &  ,CUTBIT_SNRMAX         = 25
     &  ,CUTBIT_SNRMAX2        = 26
     &  ,CUTBIT_NFIELD         = 27
     &  ,CUTBIT_MWEBV          = 28  ! May 2012
     &  ,CUTBIT_NSEASON_ACTIVE = 29  ! May 2019
     &  ,CUTBIT_REQEP          = 30  ! required epochs (9/17/2017)
     &  ,CUTBIT_PRIVATE        = 31  ! all private-var cuts
     &  ,CUTBIT_SIMVAR        = 32  !
     &  ,CUTBIT_USRCUTS        = 33  ! PASS_USRCUTS=T
     &  ,CUTBIT_OFFSET_SBFLUX  = 34  ! Mar 2012
     &  ,CUTBIT_OFFSET_SNRMAX  = 34 + MXFILT_SNRMAX
     &  ,CUTBIT_MJD_MARKER     = CUTBIT_OFFSET_SNRMAX + MXFILT_SNRMAX
     &  ,NCUTBIT_SNLC          = CUTBIT_MJD_MARKER
c SN-dependent cuts above
c MJD-dependent cuts below
     &  ,CUTBIT_PSF          = CUTBIT_MJD_MARKER + 1  ! PSF cut
     &  ,CUTBIT_ZP           = CUTBIT_MJD_MARKER + 2  ! ZP  cut (Feb 2018)
     &  ,CUTBIT_MOONDIST     = CUTBIT_MJD_MARKER + 3  ! degrees
     &  ,CUTBIT_PHOTPROB     = CUTBIT_MJD_MARKER + 4  ! Oct 23 2018
     &  ,CUTBIT_AIRMASS      = CUTBIT_MJD_MARKER + 5  ! 1=overhead
     &  ,CUTBIT_NBAND_THRESH = CUTBIT_MJD_MARKER + 6
     &  ,CUTBIT_Trest        = CUTBIT_MJD_MARKER + 7
     &  ,CUTBIT_MJD          = CUTBIT_MJD_MARKER + 8
     &  ,CUTBIT_IDTEL        = CUTBIT_MJD_MARKER + 9
c
     &  ,CUTBIT_OFFSET       = CUTBIT_IDTEL
     &  ,NCUTBIT             = CUTBIT_OFFSET
     &      )

c Define cut-masks with 64-BIT integers

      INTEGER*8
     &   CUTMASK8_SN             ! cutmask for each SN
     &  ,CUTMASK8_MJD(MXEPOCH)   ! cutmask vs. MJD, isn
     &  ,CUTMASK8_SN_ALL         ! all bits for SN cuts
     &  ,CUTMASK8_MJD_ALL        ! all bits for MJD cuts

      BYTE MJDMASK(MXEPOCH)  ! logical for each epoch,SN

      LOGICAL
     &   LSNCUTS       ! T=> passes cuts 1 to CUTBIT_MJD_MARKER
     &  ,PASS_USRCUTS  ! T=> user cuts with private code (Oct 2014)
     &  ,PASS_PRIVCUTS ! T=> pass cuts on private var in data header
     &  ,PASS_SIMCUTS  ! T=> pass cuts on SIMVAR

      INTEGER
     &   NCCID_LIST      ! size of SNCCID_LIST
     &  ,NCID_LIST       ! size of SNCID_LIST
     &  ,NCCID_IGNORE    ! size of SNCCID_IGNORE
     &  ,NCID_IGNORE     ! size of NCID_IGNORE
     &  ,CUTFLAG_SNANA   ! bits 0,1 -> LSNCUTS,LSNFITOK=T (for ntuple)
     &  ,ERRFLAG_FIT     ! error flag from fit (0=OK)

c stuff for hbook
      CHARACTER CUTVAR_NAME(NCUTBIT)*28

c define namelist cuts

      INTEGER
     &   SNTYPE_LIST(MXLISTNML)    ! (I) user list of types to select
     &  ,SNTYPE_IGNORE(MXLISTNML)  ! (I) types to ignore
     &  ,SNCID_LIST(MXLISTNML)     ! (I) user-list of integer CIDs
     &  ,CCDNUM_LIST(MXLISTNML)    ! (I) user-list of CCDNUMs (Dec 2017)
     &  ,SIM_TEMPLATE_INDEX_LIST(MXLISTNML)
     &  ,SNCID_IGNORE(MXLISTNML)   ! (I) list of CIDs to ignore
     &  ,PHOTFLAG_MSKREJ(5)        ! (I) PHOTFLAG mask-bits to reject
     &                             !     1,2 => logical OR,AND
     &  ,IDTEL_LIST(MXLISTNML)     ! computed from SNTEL_LIST
     &  ,IDFIELD_LIST(MXLISTNML)   !
     &  ,NIDTEL_LIST               ! size of list
     &  ,NIDFIELD_LIST
     &  ,NSNTYPE_LIST           ! size of SNTYPE_LIST
     &  ,NCCDNUM_LIST           ! size of CCDNUM list
     &  ,NSIM_TEMPLATE_INDEX_LIST
     &  ,NSNTYPE_IGNORE
     &  ,NFILT_SNRMAX           ! number of CUTWIN_SNRMAX_FILT cuts
     &  ,IFILT_SNRMAX(MXFILT_SNRMAX) ! store filt index for SNRMAX cuts
     &  ,NFILT_HOST_SBFLUX      ! Nfilt to cut on HOST_SBFLUX
     &  ,IFILT_HOST_SBFLUX(MXFILT_SNRMAX)

      LOGICAL DOALL_SNTEL, DOALL_SNFIELD

      character
     &    SNTEL_LIST(MXLISTNML)*20    ! (I) list of telescopes (or ALL)
     &   ,SNFIELD_LIST(MXLISTNML)*60  ! (I) list of fields to use (or 'ALL')
     &   ,SNCCID_LIST(MXLISTNML)*(MXCHAR_CCID)   ! (I) list of CIDs to process
     &   ,SNCCID_LIST_ALL(MXCCID_LIST)*(MXCHAR_CCID) ! combined list to process
     &   ,SNCCID_IGNORE(MXLISTNML)*(MXCHAR_CCID) ! (I) list of CIDs to ignore
     &   ,SNCID_IGNORE_FILE*(MXCHAR_FILENAME) ! (I) file with list of CIDs to ignore

      REAL
     &   snlc_cutvar(NCUTBIT_SNLC)
     &  ,snep_cutvar(CUTBIT_MJD_MARKER:NCUTBIT,MXEPOCH)
     &  ,cutwin_var(2,NCUTBIT)
     &  ,cutwin_cid(2)           ! candidate id
     &  ,cutwin_redshift(2)
     &  ,cutwin_redshift_err(2)
     &  ,cutwin_ra(2)
     &  ,cutwin_dec(2)
     &  ,cutwin_hostsep(2)      ! cut on host-SN sep, arcsec
     &  ,cutwin_sbflux_filt(2,MXFILT_SNRMAX)
     &  ,cutwin_Nepoch(2)       !
     &  ,cutwin_snrmax_filt(2,MXFILT_SNRMAX)   ! filled from SNCUT_SNRMAX
     &  ,cutwin_snrmax(2)       ! global SNRMAX cut
     &  ,cutwin_snrmax2(2)       ! 2nd global SNRMAX cut
     &  ,cutwin_nfield(2)        ! number of fields (usually 1)
     &  ,cutwin_mwebv(2)         ! Galactic extinction
     &  ,cutwin_nseason_active(2)
     &  ,cutwin_cutflag_reqep(2)       ! for internal use only
     &  ,cutwin_cutflag_private(2)     ! for internal use only
     &  ,cutwin_cutflag_simvar(2)      ! for internal use only
     &  ,cutwin_cutflag_usrcuts(2)     ! for internal use only
     &  ,cutwin_nfilt_snrmax(2)  ! Nfilt passing global SNRMAX cut
     &  ,cutwin_nfilt_snrmax2(2) ! Nfilt passing 2nd-best SNRMAX
     &  ,cutwin_nfilt_trestmin(2)! Nfilt passing Trestmin cut
     &  ,cutwin_nfilt_trestmax(2)! Nfilt passing Trestmax cut
     &  ,cutwin_nfilt_trest2(2)  ! Nfilt passinng Trest2 cut
     &  ,cutwin_nband_thresh(2)  ! number of ugriz bands above $band_THRESH
     &  ,cutwin_Trestmin(2)      ! window for earliest epoch, rest frame days
     &  ,cutwin_Trestmax(2)      ! window for latest epoch, rest frame day
     &  ,cutwin_TrestRange(2)    ! window on Trestmax - Trestmin
     &  ,cutwin_Trest2(2)        ! window for CUTWIN_NFILT_TREST2
     &  ,cutwin_Tgapmax(2)       ! max Trest-gap within cutwin_TREST(2)
     &  ,cutwin_T0gapmax(2)      ! max Trest-gap near peak
     &  ,cutwin_Tobsmin(2)
     &  ,cutwin_Tobsmax(2)
     &  ,cutwin_peakmjd(2)       ! cut on search peakMJD (for rates)
c
c  Below are Epoch-dependent cuts

     &  ,cutwin_psf(2)           ! PSF cut, FWHM, ARCSEC
     &  ,cutwin_zp(2)            ! ZEROPT cut, ADU or NPE
     &  ,cutwin_zpADU(2)         ! ZEROPT cut, ADU
     &  ,cutwin_zpNPE(2)         ! ZEROPT cut, NPE
     &  ,cutwin_moondist(2)      ! dist betwee SN and moon, deg
     &  ,cutwin_photprob(2)
     &  ,cutwin_Trest(2)         ! window for all epochs, rest frame days
     &  ,cutwin_MJD(2)           ! MJD window
     &  ,CUTWIN_MJD_EXCLUDE(2)          ! MJD window to exclude
     &  ,CUTWIN_MJD_INCLUDE(2)          ! require  MJD in this window
     &  ,CUTWIN_NMJD_INCLUDE(2)         ! for internal use only
     &  ,cutwin_SEARCHEFF_MASK(2) ! for SIM only
     &  ,cutwin_snrmin_filt(2,MXFILT_OBS) ! filled from EPCUT_SNRMIN
     &  ,cutwin_restlam(2)                ! cut on <LAMREST>, no cutBit
     &  ,cutwin_lamrest(2)                ! same as above
     &  ,cutwin_lamobs(2)                 ! cut on <LAMOBS>, no cutBit
c
c define character strings to specify cuts and filters
c 'SNCUT' specifies cut on each SN
c 'EPCUT' specifies selection cut on each epoch

      CHARACTER*(MXCHAR_CUTNAME)
     &  SNCUT_SNRMAX  ! max SNR required in each passband
                      ! example: 'u 10.  g 5.0  r 5.0  i 5.0  z -10.'
c
     & ,EPCUT_SNRMIN  ! min SNR accepted for each epoch/filter
                      ! example: 'u 20000.  g -4.  r -4.  i -4.  z 20000.'
c
     & ,EPCUT_SKYSIG  ! max SKY noise accepted each epoch/filter
                      ! example: 'u 20.  g 50.  r 80.  i 120.  z 200.'
c
     & ,EPCUT_PSFSIG  ! max PSF (arcsec) accepted each epoch/filter
                      ! example: 'u 1.8  g 1.5  r 1.6  i 1.8  z 2.0'
c
     & ,SNCUT_HOST_SBFLUX ! max allows surface brightness,
                          ! example 'r 1000.' -> SBFLUX < 1000

c ----------
c systematic tests for calibration:
c 'U 01 B -0.01' => shift U & B mags of primary

     &  ,MAGOBS_SHIFT_PRIMARY  ! shift primary mags  (for syst test)
     &  ,MAGOBS_SHIFT_ZP       ! shift zero points (e.g.,AB off for SDSS)
     &  ,MAGREST_SHIFT_PRIMARY ! idem for rest-frame mags
     &  ,MAGREST_SHIFT_ZP      ! idem for rest-frame filters
     &  ,FILTER_LAMSHIFT       ! e.g., 'r -2.4  i 6.2'  ! in Angstroms
c -------
c Fluxcal fudges for systematic tests (fudge photometry offsets and errors)
c Note that error is added in quadrature; ERROR<0 is subtracted in quadrature.
     &  ,FUDGE_FLUXCAL_OFFSET
     &  ,FUDGE_FLUXCAL_ERROR   ! fudge net FLUXCAL error in each band
     &  ,FUDGE_FLUXCAL_ERRPIX  ! per-pixel error --> FLUXCAL error propto PSF
     &  ,FUDGE_MAG_ERROR       ! Oct 2013: add stat error per band
     &  ,FUDGE_MAG_COVERR      ! Feb 2015: add covariant error per band

      REAL
     &   MAGOBS_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGOBS_SHIFT_ZP_USER(MXFILT_ALL)      ! user ZP shift
     &  ,MAGOBS_SHIFT_ZP_FILT(MXFILT_ALL)      ! user ZP shift + system ZP
     &  ,FLUXSCALE_ZP_FILT(MXFILT_ALL)         ! corresponding flux scale
     &  ,MAGREST_SHIFT_PRIMARY_FILT(MXFILT_ALL) ! shift mag of primary ref
     &  ,MAGREST_SHIFT_ZP_USER(MXFILT_ALL)      ! shift zero points
     &  ,MAGREST_SHIFT_ZP_FILT(MXFILT_ALL)      ! shift zero points
     &  ,FILTER_LAMSHIFT_FILT(MXFILT_ALL)       ! shift filter trans
c
     &  ,MAGOBS_SHIFT_PRIMARY_PARAMS(3)       ! poly-fun of lambda;
     &  ,MAGOBS_SHIFT_ZP_PARAMS(3)            ! A0 + A1*LAM + A2*LAM^2
     &  ,FUDGE_FLUXCAL_OFFSET_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_FLUXCAL_ERRPIX_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_ERROR_FILT(MXFILT_ALL)
     &  ,FUDGE_MAG_COVERR_FILT(MXFILT_ALL)
     &  ,MWEBV_SCALE        ! scale MW extinc for syst test
     &  ,MWEBV_SHIFT        ! shift MW extinc

      REAL
     &   RV_MWCOLORLAW   ! (I) RV for Galactic extinction

      INTEGER
     &   OPT_MWCOLORLAW    ! (I) MW color law opt (89, 94, 99 ...)
     &  ,OPT_MWEBV         ! (I) option to modify SFD maps

      LOGICAL
     &   USE_MWCOR          ! (I) T=> correct data flux for MW extinc;
                            !     F=> leave data, correct fit model for MW.
     &  ,DOBUG_LAMRANGE     ! (I) implement old lambda-range bug

      REAL
     &   REDSHIFT_FINAL_SHIFT    ! artificial shift
     &  ,FLUXERRCALC_ZPTERR

      COMMON / SNCUTCOM / CUTMASK8_SN, CUTMASK8_MJD
     &   ,CUTMASK8_SN_ALL, CUTMASK8_MJD_ALL
     &   ,cutwin_var, cutvar_name
     &   ,snlc_cutvar, snep_cutvar
     &   ,cutwin_snrmin_filt, LSNCUTS
     &   ,PASS_USRCUTS, PASS_PRIVCUTS, PASS_SIMCUTS
     &   ,CUTFLAG_SNANA , ERRFLAG_FIT
     &   ,NCCID_LIST, NCCID_IGNORE, NCID_LIST, NCID_IGNORE
     &   ,SNTYPE_LIST, SNCID_LIST, SNCCID_LIST, SNCCID_LIST_ALL
     &   ,SNCCID_IGNORE, SNCID_IGNORE_FILE, SNCID_IGNORE
     &   ,SNTYPE_IGNORE, CCDNUM_LIST, SIM_TEMPLATE_INDEX_LIST
     &   ,SNTEL_LIST,   IDTEL_LIST,   NIDTEL_LIST,   DOALL_SNTEL
     &   ,SNFIELD_LIST, IDFIELD_LIST, NIDFIELD_LIST, DOALL_SNFIELD
     &   ,NSNTYPE_LIST, NCCDNUM_LIST, NSIM_TEMPLATE_INDEX_LIST
     &   ,NSNTYPE_IGNORE, NFILT_SNRMAX, IFILT_SNRMAX
     &   ,NFILT_HOST_SBFLUX, IFILT_HOST_SBFLUX
     &   ,PHOTFLAG_MSKREJ
     &   ,USE_MWCOR, DOBUG_LAMRANGE
     &   ,MAGOBS_SHIFT_PRIMARY, MAGOBS_SHIFT_PRIMARY_FILT
     &   ,MAGOBS_SHIFT_ZP, MAGOBS_SHIFT_ZP_FILT, FLUXSCALE_ZP_FILT
     &   ,MAGOBS_SHIFT_ZP_USER
     &   ,MAGOBS_SHIFT_PRIMARY_PARAMS, MAGOBS_SHIFT_ZP_PARAMS
     &   ,MAGREST_SHIFT_PRIMARY, MAGREST_SHIFT_PRIMARY_FILT
     &   ,MAGREST_SHIFT_ZP, MAGREST_SHIFT_ZP_FILT
     &   ,MAGREST_SHIFT_ZP_USER
     &   ,FILTER_LAMSHIFT, FILTER_LAMSHIFT_FILT
     &   ,FUDGE_FLUXCAL_OFFSET,FUDGE_FLUXCAL_ERROR,FUDGE_FLUXCAL_ERRPIX
     &   ,FUDGE_FLUXCAL_OFFSET_FILT
     &   ,FUDGE_FLUXCAL_ERROR_FILT, FUDGE_FLUXCAL_ERRPIX_FILT
     &   ,FUDGE_MAG_ERROR, FUDGE_MAG_ERROR_FILT
     &   ,FUDGE_MAG_COVERR, FUDGE_MAG_COVERR_FILT
     &   ,RV_MWCOLORLAW, OPT_MWCOLORLAW, OPT_MWEBV
     &   ,MWEBV_SCALE, MWEBV_SHIFT
     &   ,REDSHIFT_FINAL_SHIFT, FLUXERRCALC_ZPTERR
     &   ,CUTWIN_MJD_EXCLUDE, CUTWIN_MJD_INCLUDE

      COMMON / BYTEMASKCOM / MJDMASK
ccc mark delete      COMMON / CIDMASKCOM / MJDMASK,CIDMASK

      EQUIVALENCE
     &    ( cutwin_var(1,cutbit_cid),       cutwin_cid )
     &   ,( cutwin_var(1,cutbit_redshift),  cutwin_redshift )
     &   ,( cutwin_var(1,cutbit_redshift_err),cutwin_redshift_err )
     &   ,( cutwin_var(1,cutbit_ra),        cutwin_ra )
     &   ,( cutwin_var(1,cutbit_dec),       cutwin_dec )
     &   ,( cutwin_var(1,cutbit_hostsep),   cutwin_hostsep )
     &   ,( cutwin_var(1,cutbit_Nepoch),    cutwin_Nepoch )
     &   ,( cutwin_var(1,cutbit_psf),       cutwin_psf )
     &   ,( cutwin_var(1,cutbit_zp),        cutwin_zp  )
     &   ,( cutwin_var(1,cutbit_moondist),  cutwin_moondist  )
     &   ,( cutwin_var(1,cutbit_photprob),  cutwin_photprob  )
     &   ,( cutwin_var(1,cutbit_Nband_thresh), cutwin_nband_thresh )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax), cutwin_nfilt_snrmax )
     &   ,( cutwin_var(1,cutbit_Nfilt_snrmax2),cutwin_nfilt_snrmax2 )
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmin),cutwin_nfilt_trestmin)
     &   ,( cutwin_var(1,cutbit_Nfilt_trestmax),cutwin_nfilt_trestmax)
     &   ,( cutwin_var(1,cutbit_Nfilt_trest2),  cutwin_nfilt_trest2)
     &   ,( cutwin_var(1,cutbit_Trestmin),  cutwin_Trestmin )
     &   ,( cutwin_var(1,cutbit_Trestmax),  cutwin_Trestmax )
     &   ,( cutwin_var(1,cutbit_TrestRange),cutwin_TrestRange )
     &   ,( cutwin_var(1,cutbit_Trest2),    cutwin_Trest2 )
     &   ,( cutwin_var(1,cutbit_Tgapmax),   cutwin_Tgapmax )
     &   ,( cutwin_var(1,cutbit_T0gapmax),  cutwin_T0gapmax )
     &   ,( cutwin_var(1,cutbit_Tobsmin),   cutwin_Tobsmin )
     &   ,( cutwin_var(1,cutbit_Tobsmax),   cutwin_Tobsmax )
     &   ,( cutwin_var(1,cutbit_Trest),     cutwin_Trest )
     &   ,( cutwin_var(1,cutbit_MJD),       cutwin_MJD   )
     &   ,( cutwin_var(1,cutbit_peakmjd),   cutwin_peakmjd )
     &   ,( cutwin_var(1,cutbit_NMJD_INCLUDE),  cutwin_NMJD_INCLUDE )
     &   ,( cutwin_var(1,cutbit_search),        cutwin_searcheff_mask )
     &   ,( cutwin_var(1,cutbit_snrmax),        cutwin_snrmax  )
     &   ,( cutwin_var(1,cutbit_snrmax2),       cutwin_snrmax2  )
     &   ,( cutwin_var(1,cutbit_nfield),        cutwin_nfield  )
     &   ,( cutwin_var(1,cutbit_mwebv),         cutwin_mwebv  )
     &   ,( cutwin_var(1,cutbit_nseason_active),cutwin_nseason_active)
     &   ,( cutwin_var(1,cutbit_reqep),        cutwin_cutflag_reqep)
     &   ,( cutwin_var(1,cutbit_private),      cutwin_cutflag_private)
     &   ,( cutwin_var(1,cutbit_simvar),       cutwin_cutflag_simvar)
     &   ,( cutwin_var(1,cutbit_usrcuts),      cutwin_cutflag_usrcuts)
     &   ,( cutwin_var(1,cutbit_offset_sbflux+1), cutwin_sbflux_filt)
     &   ,( cutwin_var(1,cutbit_offset_snrmax+1), cutwin_snrmax_filt)
c
     &   ,( cutwin_restlam(1), cutwin_lamrest(1) )



C +CDE,KCORCOM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 

C +CDE,KTABLECM,IF=R4KCOR,I2KCOR. inserted below by ypatchy. 
      REAL*4  R4KCORTABLE1D(MXTABLE1D_KCOR)
      COMMON / R4KCORCOM / R4KCORTABLE1D


      INTEGER
     &    LUNKCOR, MXKCOR
     &   ,NKDIM, N4DIM
     &   ,IDMAP_KCOR, IDMAP_AVWARP, IDMAP_LCMAG, IDMAP_MWXT
     &   ,KDIM_T, KDIM_Z, KDIM_AV, KDIM_ifiltr, KDIM_ifilto

      PARAMETER (
     &   LUNKCOR      = 98  ! log unit number to read
     &  ,MXKCOR       = 100 ! max number of KCOR tables to read
     &  ,KDIM_T       = 1  ! epoch KCOR index
     &  ,KDIM_Z       = 2  ! redshift KCOR index
     &  ,KDIM_AV      = 3  ! AV kcor index
     &  ,KDIM_ifiltr  = 4  ! rest-filter KCOR index
     &  ,KDIM_ifilto  = 5  ! obs-filter KCOR index
     &  ,NKDIM        = 5  ! number of KCOR dimensions in map
     &  ,IDMAP_KCOR   = 10 ! for INIT_1DINDEX (any integer is OK)
     &  ,IDMAP_AVWARP = 11
     &  ,IDMAP_LCMAG  = 12
     &  ,IDMAP_MWXT   = 13
     &  ,N4DIM        = 4  ! for 4-d tables
     &   )

c K-cor & photometry-template info

      INTEGER
     &   NCALL_RDKCOR ! # times RDKCOR is called
     &  ,NKCOR_STORE  ! number of stored K-corr tables
     &  ,NZBIN_KCOR   ! # Z-bins in each K-cor table
     &  ,NTBIN_KCOR   ! # T-bins "    "
     &  ,NAVBIN_KCOR  ! # AV bins
     &  ,IAV0         ! IAV index with AV=0
     &  ,IZ0          ! IZ with Z=0
     &  ,NLbin_SNSED      ! # lambda bins for SN SED
     &  ,NTbin_SNSED      ! # epoch bins for SN SED
     &  ,NERR_AVWARP      ! # errors trapped in GET_AVWARP
c
     &  ,NBINTOT_KCOR     ! Total KCOR bins
     &  ,NBINTOT_AVWARP
     &  ,NBINTOT_LCMAG
     &  ,NBINTOT_MWXT
     &  ,NBIN_KCOR(NKDIM)
     &  ,NBIN_AVWARP(NKDIM)
     &  ,NBIN_LCMAG(NKDIM)
     &  ,NBIN_MWXT(NKDIM)

      REAL
     &   Zrange_Kcor(2)    ! min,max redshift in KCOR table
     &  ,Zrange_Kcor_LU(2) ! min,max redshift for lookup
     &  ,Trange_Kcor(2)  ! min,max rest-epoch (days) in KCOR table
     &  ,AVrange_Kcor(2) ! min,max AV
     &  ,ZBINSIZE_KCOR   ! redshift binsize
     &  ,TBINSIZE_KCOR   ! epoch bin size (days)
     &  ,AVbinsize_KCOR  !
     &  ,Cbinsize_AVWARP !
     &  ,GRIDVAL_AV(MXAVBIN_KCOR)  ! store AV(iav) for convenience
     &  ,GRIDVAL_Z(MXZBIN_KCOR)    ! store Z(iz)   for convenience
     &  ,GRIDVAL_T(MXTBIN_KCOR)    ! store Trest(it) for convenience
     &  ,GRIDVAL_C(MXCBIN_AVWARP)  ! store color for AVWARP table
     &  ,Lrange_SNSED(2), LBINSIZE_SNSED
     &  ,Trange_SNSED(2), TBINSIZE_SNSED
     &  ,FLUX_SNSED(MXLAMBIN_SNSED*MXTBIN_KCOR)

      REAL*4
     &   AVWARP_TABLE1D(MXTABLE1D_AVWARP)
     &  ,LCMAG_TABLE1D(MXTABLE1D_LCMAG)
     &  ,MWXT_TABLE1D(MXTABLE1D_MWXT)  ! MilkyWay extinction table

      INTEGER
     &   HIDMAG(MXAVBIN_KCOR,MXFILT_ALL)
     &  ,HIDKCOR(MXFILT_ALL,MXFILT_ALL)

      LOGICAL
     &   USE_AVWARPTABLE
     &  ,LRDZPOFF
     &  ,RDKCOR_STANDALONE
     &  ,EXIST_KCOR(MXFILT_ALL,MXFILT_ALL)

      CHARACTER RESTKCOR_FILTERS*20    ! rest-frame filters from KCOR file

      REAL RVMW_RDKCOR ! RV used for MW table extension
      INTEGER OPT_MWCOLORLAW_RDKCOR

      COMMON / KCORCOM / NCALL_RDKCOR, NKCOR_STORE, RDKCOR_STANDALONE
     &      ,NBIN_KCOR,   NBINTOT_KCOR
     &      ,NBIN_AVWARP, NBINTOT_AVWARP
     &      ,NBIN_LCMAG,  NBINTOT_LCMAG
     &      ,NBIN_MWXT,   NBINTOT_MWXT
     &      ,NZBIN_KCOR,    NTBIN_KCOR,    NAVBIN_KCOR
     &      ,Zrange_KCOR,   Trange_KCOR,   AVrange_KCOR
     &      ,ZBINSIZE_KCOR, TBINSIZE_KCOR, AVbinsize_KCOR
     &      ,CBINSIZE_AVWARP, IAV0, IZ0
     &      ,GRIDVAL_Z, GRIDVAL_T, GRIDVAL_AV, GRIDVAL_C
     &      ,HIDKCOR, HIDMAG
     &      ,Zrange_KCOR_LU, RESTKCOR_FILTERS
     &      ,USE_AVWARPTABLE, LRDZPOFF, NERR_AVWARP
     &      ,NLbin_SNSED, Lrange_SNSED, LBINSIZE_SNSED
     &      ,NTbin_SNSED, Trange_SNSED, TBINSIZE_SNSED
     &      ,FLUX_SNSED
     &      ,LCMAG_TABLE1D, MWXT_TABLE1D, AVWARP_TABLE1D
     &      ,EXIST_KCOR
     &      ,RVMW_RDKCOR, OPT_MWCOLORLAW_RDKCOR

c  July 2016 - define optional spectrograph info
      CHARACTER
     &   SPECTROGRAPH_INSTRUMENT*80
     &  ,SPECTROGRAPH_FILTERLIST*(MXFILT_ALL)
      INTEGER
     &   NFILTDEF_SPECTROGRAPH
     &  ,IFILTDEF_SPECTROGRAPH(MXFILT_ALL)

      COMMON / SPECTROGRAPH_COM /
     &   SPECTROGRAPH_INSTRUMENT, SPECTROGRAPH_FILTERLIST
     &  ,NFILTDEF_SPECTROGRAPH,   IFILTDEF_SPECTROGRAPH

c define temporary 'XXX_RDKCOR' arrays that are used only
c for reading and parsing.

      CHARACTER
     &   PRIMARY_NAME_RDKCOR(10)*40          ! e.g, AB, BD17 ..
     &  ,FILTER_NAME_RDKCOR(MXFILT_ALL)*40   ! e.g, SDSS-g, Bessell90-V
     &  ,KCORINFO_STRING_RDKCOR(MXKCOR)*60   ! Kcor K_xy obsfilt restfilt

      INTEGER
     &   IVER_RDKCOR
     &  ,NPRIM_RDKCOR, NKCOR_RDKCOR
     &  ,INDX_PRIMARY_RDKCOR(MXFILT_ALL)
     &  ,NFILTDEF_RDKCOR               ! Number of all filters
     &  ,IFILTDEF_RDKCOR(MXFILT_ALL)   ! vs. ifilt_rdkcor
     &  ,NFILTOBS_RDKCOR               ! number of obs filters
     &  ,IFILTOBS_RDKCOR(MXKCOR)       ! vs. ikcor
     &  ,IFILT2_RDKCOR(2,MXKCOR)       ! 1=rest:2=obs, 1:NKCOR_STORE
     &  ,MASKFILT_RDKCOR(MXFILT_ALL)   ! bit0,1 -> rest,obs
     &  ,NFILT_DUPLICATE_RDKCOR        ! keep count of duplicate filters
     &  ,IKCOR_RDKCOR(MXKCOR)          ! ikcor vs. istore
     &  ,IPRIM_REF_RDKCOR              ! index for primary  ref.

      REAL
     &   MAG_PRIMARY_RDKCOR(MXFILT_ALL)    ! native mags
     &  ,ZPOFF_PRIMARY_RDKCOR(MXFILT_ALL)  ! mag(native) - mag(synth)
     &  ,ZPOFF_SNPHOT_RDKCOR(MXFILT_ALL)   ! apply to data (from ZPOFF.DAT)

      LOGICAL
     &   ISFITS_RDKCOR
     &  ,LVERBOSE_RDKCOR
     &  ,ISLAMSHIFT_RDKCOR(0:MXFILT_ALL)

      REAL      NULLF_RDKCOR
      INTEGER   NULLI_RDKCOR
      CHARACTER NULLS_RDKCOR*20
      LOGICAL   ANYF_RDKCOR

      COMMON / RDKCORCOM /
     &   IVER_RDKCOR, ISFITS_RDKCOR,LVERBOSE_RDKCOR,ISLAMSHIFT_RDKCOR
     &  ,NKCOR_RDKCOR
     &  ,NPRIM_RDKCOR, INDX_PRIMARY_RDKCOR, PRIMARY_NAME_RDKCOR
     &  ,FILTER_NAME_RDKCOR, KCORINFO_STRING_RDKCOR
     &  ,IFILTDEF_RDKCOR,  NFILTDEF_RDKCOR
     &  ,IFILTOBS_RDKCOR,  NFILTOBS_RDKCOR
     &  ,IFILT2_RDKCOR
     &  ,NFILT_DUPLICATE_RDKCOR, IKCOR_RDKCOR, MASKFILT_RDKCOR
     &  ,IPRIM_REF_RDKCOR
c
     &  ,MAG_PRIMARY_RDKCOR
     &  ,ZPOFF_PRIMARY_RDKCOR
     &  ,ZPOFF_SNPHOT_RDKCOR
c
     &  ,NULLF_RDKCOR,NULLI_RDKCOR,NULLS_RDKCOR,ANYF_RDKCOR


C +CDE,SNHOSTCOM. inserted below by ypatchy. 

c Host galaxy parameters.
c Dec 17, 2012 - add SNHOST_MAGOBS


c logical flags for CWNT
      LOGICAL
     &   EXIST_SNHOST_ANGSEP
     &  ,EXIST_SNHOST_DDLR
     &  ,EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT
     &  ,EXIST_SNHOST_LOGMASS
     &  ,EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS
     &  ,EXIST_SNHOST_SB

      INTEGER*8  SNHOST_OBJID(MXSNHOST)    ! int id
      REAL*8    DSNHOST_OBJID(MXSNHOST)   ! for tables only

      REAL
     &   SNHOST_ANGSEP(MXSNHOST)      ! SN-host sep, arcsec
     &  ,SNHOST_DDLR(MXSNHOST)        ! SNSEP/DLR
     &  ,SNHOST_CONFUSION             ! HC analog from Gupta 2016
     &  ,SNHOST_ZPHOT(MXSNHOST), SNHOST_ZPHOT_ERR(MXSNHOST)
     &  ,SNHOST_ZSPEC(MXSNHOST), SNHOST_ZSPEC_ERR(MXSNHOST)
     &  ,SNHOST_LOGMASS(MXSNHOST)
     &  ,SNHOST_LOGMASS_ERR(MXSNHOST)
     &  ,SNHOST_sSFR(MXSNHOST)
     &  ,SNHOST_sSFR_ERR(MXSNHOST)
     &  ,SNHOST_SBFLUXCAL(MXFILT_ALL)  ! surface brightness FLUXCAL /asec^2
     &  ,SNHOST_MAGOBS(MXFILT_ALL,(MXSNHOST))       ! observer-frame mags
     &  ,SNHOST_MAGOBS_ERR(MXFILT_ALL,(MXSNHOST))   ! error on above

      REAL*8
     &    SNHOST8_RA(MXSNHOST)
     &   ,SNHOST8_DEC(MXSNHOST)

c quantites which do not depend on which host
      INTEGER*4  SNHOST_NMATCH, SNHOST_NMATCH2  ! number of host matches
      REAL
     &   SNHOST_SBFLUXCAL_ERR(MXFILT_ALL)
     &  ,SNHOST_SBMAG(MXFILT_ALL)        ! surface brightness mag/asec^2

      COMMON / SNHOSTCOM /
     &   EXIST_SNHOST_ANGSEP, EXIST_SNHOST_DDLR, EXIST_SNHOST_CONFUSION
     &  ,EXIST_SNHOST_ZPHOT,  EXIST_SNHOST_LOGMASS, EXIST_SNHOST_sSFR
     &  ,EXIST_SNHOST_MAGOBS, EXIST_SNHOST_SB
c
     &  ,SNHOST_NMATCH, SNHOST_NMATCH2
     &  ,SNHOST_ANGSEP, SNHOST_DDLR, SNHOST_CONFUSION
     &  ,SNHOST_ZPHOT,   SNHOST_ZPHOT_ERR
     &  ,SNHOST_ZSPEC,   SNHOST_ZSPEC_ERR
     &  ,SNHOST_LOGMASS, SNHOST_LOGMASS_ERR
     &  ,SNHOST_sSFR,    SNHOST_sSFR_ERR
     &  ,SNHOST_SBFLUXCAL, SNHOST_SBFLUXCAL_ERR, SNHOST_SBMAG
     &  ,SNHOST_MAGOBS, SNHOST_MAGOBS_ERR

      COMMON / SNHOSTCOM8 /
     &   SNHOST_OBJID, DSNHOST_OBJID, SNHOST8_RA, SNHOST8_DEC

cc +CDE,PARSECOM.

C +CDE,SNSIMCOM. inserted below by ypatchy. 

c simulation parameters (if FAKE=2)

      REAL*8  SIM8_RA, SIM8_DECL

      REAL
     &   SIM_REDSHIFT_HELIO
     &  ,SIM_REDSHIFT_CMB
     &  ,SIM_REDSHIFT_HOST
     &  ,SIM_VPEC
     &  ,SIM_DLMAG
     &  ,SIM_LENSDMU
     &  ,SIM_MWEBV
     &  ,SIM_MWRV
     &  ,SIM_SALT2x0
     &  ,SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW
     &  ,SIM_PEAKMJD
     &  ,SIM_EXPOSURE_TIME(MXFILT_OBS)  ! relative exposure time
     &  ,SIM_PEAKMAG(MXFILT_OBS)
     &  ,SIM_EPMAGOBS(MXEPOCH)      ! true epoch mag at each filter/epoch
     &  ,SIM_EPFLUXCAL(MXEPOCH)     ! true fluxcal at each filter/epoch
     &  ,SIM_EPSNRMON(MXEPOCH)      ! optional SNR at MAGMONITOR
     &  ,SIM_EPMAGREST(MXEPOCH)     ! true rest-frame mag
     &  ,SIM_EPFLUXCAL_HOSTERR(MXEPOCH) ! true error from host noise
     &  ,SIMSED_PARVAL(MXPAR_SIMSED)
     &  ,BYOSED_PARVAL(MXPAR_SIMSED)    ! Dec 10 2018
     &  ,LCLIB_PARVAL(MXPAR_LCLIB)
     &  ,SIM_HOSTLIB_PARVAL(MXPAR_SIMSED) ! HOSTLIB params
     &  ,SIM_MAGSMEAR_COH
     &  ,SIM_TEMPLATEMAG(MXFILT_ALL)  ! image-sub template, not LCLIB template
     &  ,SIM_LCWIDTH(MXFILT_ALL)      ! computed from SIM_EPFLUXCAL
c
     &  ,SIM_MODELGRID_TOBS(MXEP_MODELGRID)    ! for SNANA+SIM_MAGOBS table
     &  ,SIM_MODELGRID_MAGOBS(MXEP_MODELGRID,MXFILT_OBS/2) ! idem

      INTEGER
     &   SIM_MODEL_INDEX      ! model index (MLCS,SALT2,NON1a ...)
     &  ,SIM_GENTYPE          ! generated SNTYPE; SIM_TYPE_INDEX
     &  ,SIM_TEMPLATE_INDEX   ! template index for NON1ASED, SIMSED, LCLIB ...
     &  ,SIM_SEARCHEFF_MASK   ! bits 1,2 => found by software,humans
     &  ,SIM_LIBID            ! LIBID for each event
     &  ,SIM_NGEN_LIBID       ! NGEN for this LIBID (usually 1)
     &  ,SIM_NOBS_UNDEFINED   ! NOBS where model is undefined
     &  ,SIM_NSUBSAMPLE_MARK  ! Number of marked sub-samples
     &  ,SIM_SUBSAMPLE_INDEX  ! sub-sample index
     &  ,SIM_REDSHIFT_FLAG    ! points to source of redshift
     &  ,SIMOPT_MWCOLORLAW    ! option for MW color law
     &  ,SIMOPT_MWEBV         ! option to modify MWEBV_SFD map
     &  ,NPAR_SIMSED
     &  ,NPAR_BYOSED
     &  ,NPAR_LCLIB
     &  ,NPAR_SIM_HOSTLIB
     &  ,SIM_EPFILTREST(MXEPOCH)
     &  ,SIMLIB_MSKOPT        ! SIMLIB option mask (Dec 2015)
     &  ,NEP_SIM_MODELGRID

      INTEGER*8  SIM_HOSTLIB_GALID

      CHARACTER
     &   SIMNAME_MODEL*40     ! SALT2, mlcs2k2, NON1A, etc ...
     &  ,SIMNAME_TYPE*12      ! Ia, Ib, II, IIN, etc ...
     &  ,SIMNAME_SHAPEPAR*40  ! DELTA, x1, stretch ...
     &  ,SIMNAME_SHAPELAW*40  ! alpha ...
     &  ,SIMNAME_COLORPAR*40  ! AV, c ...
     &  ,SIMNAME_COLORLAW*40  ! RV,  beta ...
     &  ,SIMLIB_FILENAME*200  ! SIMLIB file used to generate sim
     &  ,SIMSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,SIMSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,BYOSED_KEYWORD(MXPAR_SIMSED)*80 ! full keyname in header
     &  ,BYOSED_PARNAME(MXPAR_SIMSED)*20 ! parname for fitres and ntuple
     &  ,LCLIB_KEYWORD(MXPAR_LCLIB)*80 ! full keyname in header
     &  ,LCLIB_PARNAME(MXPAR_LCLIB)*20 ! parname for fitres and ntuple
     &  ,SIM_HOSTLIB_KEYWORD(MXPAR_SIMSED)*60 ! full keyname
     &  ,SIM_HOSTLIB_PARNAME(MXPAR_SIMSED)*60 ! parname
     &  ,SIMNAME_SNRMON*40

      COMMON / SNSIMCOM / SIMNAME_MODEL, SIM_MODEL_INDEX
     &  ,SIM_REDSHIFT_HELIO, SIM_REDSHIFT_CMB, SIM_REDSHIFT_HOST
     &  ,SIM_REDSHIFT_FLAG
     &  ,SIM_VPEC, SIM_DLMAG, SIM_LENSDMU, SIM_SALT2x0, SIM_SALT2mb
     &  ,SIM_COLORPAR, SIM_COLORLAW, SIM_AV, SIM_RV
     &  ,SIM_SHAPEPAR, SIM_SHAPELAW, SIM_PEAKMJD, SIM_MWEBV,SIM_MWRV
     &  ,SIM_EXPOSURE_TIME, SIM_TEMPLATE_INDEX, SIM_TEMPLATEMAG
     &  ,SIM_GENTYPE, SIMNAME_TYPE, SIM_LCWIDTH
     &  ,SIM_SEARCHEFF_MASK, SIM_LIBID, SIM_NGEN_LIBID
     &  ,SIM_NOBS_UNDEFINED, SIM_NSUBSAMPLE_MARK, SIM_SUBSAMPLE_INDEX
     &  ,SIM_PEAKMAG, SIM_EPMAGOBS, SIM_EPFLUXCAL
     &  ,SIM_EPSNRMON, SIM_EPFLUXCAL_HOSTERR
     &  ,SIMNAME_SHAPEPAR,  SIMNAME_SHAPELAW
     &  ,SIMNAME_COLORPAR,  SIMNAME_COLORLAW, SIMNAME_SNRMON
     &  ,SIMLIB_FILENAME, SIMLIB_MSKOPT
     &  ,SIMOPT_MWCOLORLAW, SIMOPT_MWEBV
     &  ,SIM_EPFILTREST, SIM_EPMAGREST,  SIM_MAGSMEAR_COH
     &  ,NEP_SIM_MODELGRID, SIM_MODELGRID_TOBS, SIM_MODELGRID_MAGOBS

      COMMON / SIMSEDCOM /
     &  NPAR_SIMSED, SIMSED_PARVAL, SIMSED_PARNAME, SIMSED_KEYWORD

      COMMON / BYOSEDCOM /
     &   NPAR_BYOSED, BYOSED_PARVAL, BYOSED_PARNAME, BYOSED_KEYWORD

      COMMON / LCLIBCOM /
     &   NPAR_LCLIB, LCLIB_PARVAL, LCLIB_KEYWORD, LCLIB_PARNAME

      COMMON / SIMHOSTCOM /
     &    NPAR_SIM_HOSTLIB, SIM_HOSTLIB_PARVAL
     &   ,SIM_HOSTLIB_PARNAME, SIM_HOSTLIB_KEYWORD

      COMMON / SNSIMCOM8 /
     &         SIM8_RA, SIM8_DECL, SIM_HOSTLIB_GALID



C +CDE,DUMPCOM. inserted below by ypatchy. 

      REAL
     &   REDSHIFT_DMP_BINSIZE
     &  ,TREST_DMP_MIN
     &  ,TREST_DMP_MAX
     &  ,TREST_DMP_BINSIZE

      PARAMETER (
     &   REDSHIFT_DMP_BINSIZE = 0.01
     &  ,TREST_DMP_MIN        = -20.0
     &  ,TREST_DMP_MAX        = +80.0
     &  ,TREST_DMP_BINSIZE    =   1.0
     &    )


      CHARACTER
     &   STRING_Kxy*8, STRING_REST_COLOR*4  ! Apr 2016

      INTEGER
     &   IFILT_REST_DMP       ! (I) rest-frame filter index
     &  ,IFILT_OBS_DMP        ! (I) observer frame filter index
     &  ,IFILT_RESTCOLOR_DMP(2) ! (I) rest-color (default=0.0) (Apr 2016)

      REAL
     &   Trest_DMP       ! (I) rest-frame epoch, days
     &  ,redshift_DMP    ! (I) redshift
     &  ,REDSHIFT_DMP_RANGE(2)
     &  ,TREST_DMP_RANGE(2)
     &  ,MAGREST_COLOR, MAGREST(2)

      INTEGER
     &   NZBIN_DMP, NTBIN_DMP

      COMMON / DUMPCOM /
     &   STRING_Kxy,     STRING_REST_COLOR, MAGREST_COLOR, MAGREST
     &  ,IFILT_REST_DMP, IFILT_OBS_DMP, IFILT_RESTCOLOR_DMP
     &  ,Trest_DMP,      redshift_DMP
     &  ,TREST_DMP_RANGE,   REDSHIFT_DMP_RANGE
     &  ,NZBIN_DMP,         NTBIN_DMP

C ==============================================
      REAL DIF, BIN
C ---------- BEGIN ----------------

      REDSHIFT_DMP_RANGE(1) = GRIDVAL_Z(1)
      REDSHIFT_DMP_RANGE(2) = GRIDVAL_Z(NZBIN_KCOR)
      dif = REDSHIFT_DMP_RANGE(2) - REDSHIFT_DMP_RANGE(1)
      bin = REDSHIFT_DMP_BINSIZE
      NZBIN_DMP = int(dif/bin + 0.0001) + 1

      RETURN
      END

C =============================
CDECK  ID>,  SET_TREST_ALL. 
      SUBROUTINE SET_TREST_ALL
      IMPLICIT NONE

C +CDE,DUMPCOM. inserted below by ypatchy. 

      REAL
     &   REDSHIFT_DMP_BINSIZE
     &  ,TREST_DMP_MIN
     &  ,TREST_DMP_MAX
     &  ,TREST_DMP_BINSIZE

      PARAMETER (
     &   REDSHIFT_DMP_BINSIZE = 0.01
     &  ,TREST_DMP_MIN        = -20.0
     &  ,TREST_DMP_MAX        = +80.0
     &  ,TREST_DMP_BINSIZE    =   1.0
     &    )


      CHARACTER
     &   STRING_Kxy*8, STRING_REST_COLOR*4  ! Apr 2016

      INTEGER
     &   IFILT_REST_DMP       ! (I) rest-frame filter index
     &  ,IFILT_OBS_DMP        ! (I) observer frame filter index
     &  ,IFILT_RESTCOLOR_DMP(2) ! (I) rest-color (default=0.0) (Apr 2016)

      REAL
     &   Trest_DMP       ! (I) rest-frame epoch, days
     &  ,redshift_DMP    ! (I) redshift
     &  ,REDSHIFT_DMP_RANGE(2)
     &  ,TREST_DMP_RANGE(2)
     &  ,MAGREST_COLOR, MAGREST(2)

      INTEGER
     &   NZBIN_DMP, NTBIN_DMP

      COMMON / DUMPCOM /
     &   STRING_Kxy,     STRING_REST_COLOR, MAGREST_COLOR, MAGREST
     &  ,IFILT_REST_DMP, IFILT_OBS_DMP, IFILT_RESTCOLOR_DMP
     &  ,Trest_DMP,      redshift_DMP
     &  ,TREST_DMP_RANGE,   REDSHIFT_DMP_RANGE
     &  ,NZBIN_DMP,         NTBIN_DMP

C ==============================================
      REAL DIF
C ---------- BEGIN ----------------
      TREST_DMP          = 0.0
      TREST_DMP_RANGE(1) = TREST_DMP_MIN
      TREST_DMP_RANGE(2) = TREST_DMP_MAX
      DIF = (TREST_DMP_MAX - TREST_DMP_MIN)
      NTBIN_DMP = INT(DIF / TREST_DMP_BINSIZE + 0.001) + 1
      RETURN
      END

