*   Jun 24 2017: MNI->80 (was 50) for SALT2mu program
*   Sep 11 2017: MNI->110
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MINUIT(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
C  CPNAM   Parameter name (10 characters)
C  U       External (visible to user in FCN) value of parameter
C  ALIM, BLIM Lower and upper parameter limits. If both zero, no limits.
C  ERP,ERN Positive and negative MINOS errors, if calculated.
C  WERR    External parameter error (standard deviation, defined by UP)
C  GLOBCC  Global Correlation Coefficient
C  NVARL   =-1 if parameter undefined,      =0 if constant,
C          = 1 if variable without limits,  =4 if variable with limits
C   (Note that if parameter has been fixed, NVARL=1 or =4, and NIOFEX=0)
C  NIOFEX  Internal parameter number, or zero if not currently variable
C  NEXOFI  External parameter number for currently variable parameters
C  X, XT   Internal parameter values (X are sometimes saved in XT)
C  DIRIN   (Internal) step sizes for current step
C  variables with names ending in ..S are saved values for fixed params
C  VHMAT   (Internal) error matrix stored as Half MATrix, since
C                it is symmetric
C  VTHMAT  VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
C
C  ISW definitions:
C      ISW(1) =0 normally, =1 means CALL LIMIT EXCEEDED
C      ISW(2) =0 means no error matrix
C             =1 means only approximate error matrix
C             =2 means full error matrix, but forced pos-def.
C             =3 means good normal full error matrix exists
C      ISW(3) =0 if Minuit is calculating the first derivatives
C             =1 if first derivatives calculated inside FCN
C      ISW(4) =-1 if most recent minimization did not converge.
C             = 0 if problem redefined since most recent minimization.
C             =+1 if most recent minimization did converge.
C      ISW(5) is the PRInt level.  See SHO PRIntlevel
C      ISW(6) = 0 for batch mode, =1 for interactive mode
C                      =-1 for originally interactive temporarily batch
C
C  LWARN is true if warning messges are to be put out (default=true)
C            SET WARN turns it on, set NOWarn turns it off
C  LREPOR is true if exceptional conditions are put out (default=false)
C            SET DEBUG turns it on, SET NODebug turns it off
C  LIMSET is true if a parameter is up against limits (for MINOS)
C  LNOLIM is true if there are no limits on any parameters (not yet used)
C  LNEWMN is true if the previous process has unexpectedly improved FCN
C  LPHEAD is true if a heading should be put out for the next parameter
C        definition, false if a parameter has just been defined
C
      EXTERNAL FCN,FUTIL
      CHARACTER*40 CWHYXT
      DATA CWHYXT/'FOR UNKNOWN REASONS                     '/
      DATA JSYSRD,JSYSWR,JSYSSA/5,6,7/
C                                 . . . . . . . . . . initialize minuit
      WRITE (JSYSWR,'(1X,75(1H*))')
      CALL FLUSH(6)
      CALL MNINIT (JSYSRD,JSYSWR,JSYSSA)
C                                      . . . . initialize new data block
  100 CONTINUE
      WRITE (ISYSWR,'(1X,75(1H*))')
      NBLOCK = NBLOCK + 1
      WRITE (ISYSWR,'(26X,A,I4)')  'MINUIT DATA BLOCK NO.',NBLOCK
      WRITE (ISYSWR,'(1X,75(1H*))')
      CALL FLUSH(6)
C               . . . . . . . . . . .   set parameter lists to undefined
      CALL MNCLER
C                                             . . . . . . . . read title
      CALL MNREAD(FCN,1,IFLGUT,FUTIL)
      IF (IFLGUT .EQ. 2)  GO TO 500
      IF (IFLGUT .EQ. 3)  GO TO 600
C                                        . . . . . . . . read parameters
      CALL MNREAD(FCN,2,IFLGUT,FUTIL)
      IF (IFLGUT .EQ. 2)  GO TO 500
      IF (IFLGUT .EQ. 3)  GO TO 600
      IF (IFLGUT .EQ. 4)  GO TO 700
C                              . . . . . . verify FCN not time-dependent
      WRITE (ISYSWR,'(/A,A)') ' MINUIT: FIRST CALL TO USER FUNCTION,',
     +    ' WITH IFLAG=1'
      NPARX = NPAR
      CALL MNINEX(X)
      FZERO = UNDEFI
      CALL FCN(NPARX,GIN,FZERO,U,1,FUTIL)
      FIRST = UNDEFI
      CALL FCN(NPARX,GIN,FIRST,U,4,FUTIL)
      NFCN = 2
      IF (FZERO.EQ.UNDEFI .AND. FIRST.EQ.UNDEFI)  THEN
          CWHYXT = 'BY ERROR IN USER FUNCTION.  '
          WRITE (ISYSWR,'(/A,A/)') ' USER HAS NOT CALCULATED FUNCTION',
     +    ' VALUE WHEN IFLAG=1 OR 4'
          GO TO 800
      ENDIF
      AMIN = FIRST
      IF (FIRST .EQ. UNDEFI) AMIN=FZERO
      CALL MNPRIN(1,AMIN)
      NFCN = 2
      IF (FIRST .EQ. FZERO)  GO TO 300
      FNEW = 0.0
      CALL FCN(NPARX,GIN,FNEW,U,4,FUTIL)
      IF  (FNEW .NE. AMIN) WRITE (ISYSWR,280) AMIN, FNEW
  280 FORMAT (/' MINUIT WARNING: PROBABLE ERROR IN USER FUNCTION.'/
     +         ' FOR FIXED VALUES OF PARAMETERS, FCN IS TIME-DEPENDENT'/
     +         ' F =',E22.14,' FOR FIRST CALL'/
     +         ' F =',E22.14,' FOR SECOND CALL.'/)
      NFCN = 3
  300 FVAL3 = 2.0*AMIN+1.0
C             
      CALL MNREAD(FCN,3,IFLGUT,FUTIL)
      IF (IFLGUT .EQ. 2)  GO TO 500
      IF (IFLGUT .EQ. 3)  GO TO 600
      IF (IFLGUT .EQ. 4)  GO TO 700
      CWHYXT = 'BY MINUIT COMMAND: '//CWORD
      IF (INDEX(CWORD,'STOP').GT. 0)  GO TO 800
      IF (INDEX(CWORD,'EXI') .GT. 0)  GO TO 800
      IF (INDEX(CWORD,'RET') .EQ. 0)  GO TO 100
      CWHYXT = 'AND RETURNS TO USER PROGRAM.    '
      WRITE (ISYSWR,'(A,A)')  ' ..........MINUIT TERMINATED ',CWHYXT
      CALL FLUSH(6)
      RETURN
C                                           . . . . . . stop conditions
  500 CONTINUE
      CWHYXT = 'BY END-OF-DATA ON PRIMARY INPUT FILE.   '
      GO TO 800
  600 CONTINUE
      CWHYXT = 'BY UNRECOVERABLE READ ERROR ON INPUT.   '
      GO TO 800
  700 CONTINUE
      CWHYXT = ': FATAL ERROR IN PARAMETER DEFINITIONS. '
  800 WRITE (ISYSWR,'(A,A)')  ' ..........MINUIT TERMINATED ',CWHYXT
      CALL FLUSH(6)
      STOP
C
C  ......................entry to set unit numbers  - - - - - - - - - -
      ENTRY MINTIO(I1,I2,I3)
      JSYSRD = I1
      JSYSWR = I2
      JSYSSA = I3
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNAMIN(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called  from many places.  Initializes the value of AMIN by
CC        calling the user function. Prints out the function value and
CC        parameter values if Print Flag value is high enough.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      NPARX = NPAR
      IF (ISW(5) .GE. 1) WRITE (ISYSWR,'(/A,A)') ' FIRST CALL TO ',
     + 'USER FUNCTION AT NEW START POINT, WITH IFLAG=4.'
      CALL MNEXIN(X)
      CALL FCN(NPARX,GIN,FNEW,U,4,FUTIL)
      NFCN = NFCN + 1
      AMIN = FNEW
      EDM = BIGEDM
      CALL FLUSH(6)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNBINS(A1,A2,NAA,BL,BH,NB,BWID)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C         SUBROUTINE TO DETERMINE REASONABLE HISTOGRAM INTERVALS
C         GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2
C         AND DESIRED MAXIMUM NUMBER OF BINS NAA
C         PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID
C         F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988
      PARAMETER (ZERO=0.0, ONE=1.0)
      AL = MIN(A1,A2)
      AH = MAX(A1,A2)
      IF (AL.EQ.AH)  AH = AL + 1.
C         IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE
      IF (NAA .EQ. -1)  GO TO 150
   10 NA = NAA - 1
      IF (NA .LT. 1)  NA = 1
C          GET NOMINAL BIN WIDTH IN EXPON FORM
   20 AWID = (AH-AL)/FLOAT(NA)
      LOG = INT(DLOG10(DBLE(AWID)))
      IF (AWID .LE. ONE)  LOG=LOG-1
      SIGFIG = AWID * (10.00 **(-LOG))
C         ROUND MANTISSA UP TO 2, 2.5, 5, OR 10
      IF(SIGFIG .GT. 2.0)  GO TO 40
      SIGRND = 2.0
      GO TO 100
   40 IF (SIGFIG .GT. 2.5)  GO TO 50
      SIGRND = 2.5
      GO TO 100
   50 IF(SIGFIG .GT. 5.0)  GO TO 60
      SIGRND =5.0
      GO TO 100
   60 SIGRND = 1.0
      LOG = LOG + 1
  100 CONTINUE
      BWID = SIGRND*10.0**LOG
      GO TO 200
C         GET NEW BOUNDS FROM NEW WIDTH BWID
  150 IF (BWID .LE. ZERO)  GO TO 10
  200 CONTINUE
      ALB = AL/BWID
      LWID=ALB
      IF (ALB .LT. ZERO)  LWID=LWID-1
      BL = BWID*FLOAT(LWID)
      ALB = AH/BWID + 1.0
      KWID = ALB
      IF (ALB .LT. ZERO)  KWID=KWID-1
      BH = BWID*FLOAT(KWID)
      NB = KWID-LWID
      IF (NAA .GT. 5)  GO TO 240
      IF (NAA .EQ. -1)  RETURN
C          REQUEST FOR ONE BIN IS DIFFICULT CASE
      IF (NAA .GT. 1 .OR. NB .EQ. 1)  RETURN
      BWID =  BWID*2.0
       NB  = 1
       RETURN
  240 IF (2*NB .NE. NAA)  RETURN
      NA = NA + 1
      GO TO 20
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNCALF(FCN,PVEC,YCALF,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called only from MNIMPR.  Transforms the function FCN
CC        by dividing out the quadratic part in order to find further
CC        minima.    Calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin)
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION PVEC(15)
      NPARX = NPAR
      CALL MNINEX(PVEC)
      CALL FCN(NPARX,GIN,F,U,4,FUTIL)
      NFCN = NFCN + 1
      DO 200 I= 1, NPAR
      GRD(I) = 0.
         DO 200 J= 1, NPAR
         M = MAX(I,J)
         N = MIN(I,J)
         NDEX = M*(M-1)/2 + N
  200    GRD(I) = GRD(I) + VTHMAT(NDEX) * (XT(J)-PVEC(J))
      DENOM = 0.
      DO 210 I= 1, NPAR
  210 DENOM = DENOM + GRD(I) * (XT(I)-PVEC(I))
      IF (DENOM .LE. ZERO)  THEN
         DCOVAR = 1.
         ISW(2) = 0
         DENOM = 1.0
      ENDIF
      YCALF = (F-APSI) / DENOM
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNCLER
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MINUIT and by option from MNEXCM
CC        Resets the parameter list to UNDEFINED
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      NPFIX = 0
      NU = 0
      NPAR = 0
      NFCN = 0
      NWRMES(1) = 0
      NWRMES(2) = 0
      DO 10 I= 1, MAXEXT
      U(I) = 0.0
      CPNAM(I) = CUNDEF
      NVARL(I) = -1
   10 NIOFEX(I) = 0
      CALL MNRSET(1)
      CFROM = 'CLEAR   '
      NFCNFR = NFCN
      CSTATU ='UNDEFINED '
      LNOLIM = .TRUE.
      LPHEAD = .TRUE.
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNCNTR(FCN,KE1,KE2,IERRF,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       to print function contours in two variables, on line printer
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      PARAMETER (NUMBCS=20,NXMAX=115)
      DIMENSION CONTUR(NUMBCS), FCNA(NXMAX),FCNB(NXMAX)
      CHARACTER CLABEL*(NUMBCS)
      CHARACTER CHLN*(NXMAX),CHMID*(NXMAX),CHZERO*(NXMAX)
      DATA CLABEL/'0123456789ABCDEFGHIJ'/
C                 input arguments: parx, pary, devs, ngrid
      IF (KE1.LE.0 .OR. KE2.LE.0)  GO TO 1350
      IF (KE1.GT.NU .OR. KE2.GT.NU)  GO TO 1350
      KI1 = NIOFEX(KE1)
      KI2 = NIOFEX(KE2)
      IF (KI1.LE.0 .OR. KI2.LE.0)  GO TO 1350
      IF (KI1 .EQ. KI2)  GO TO 1350
C
      IF (ISW(2) .LT. 1)  THEN
          CALL MNHESS(FCN,FUTIL)
          CALL MNWERR
          ENDIF
      NPARX = NPAR
      XSAV = U(KE1)
      YSAV = U(KE2)
      DEVS = WORD7(3)
      IF (DEVS .LE. ZERO)  DEVS=2.
      XLO = U(KE1) - DEVS*WERR(KI1)
      XUP = U(KE1) + DEVS*WERR(KI1)
      YLO = U(KE2) - DEVS*WERR(KI2)
      YUP = U(KE2) + DEVS*WERR(KI2)
      NGRID = WORD7(4)
      IF (NGRID .LE. 0)  THEN
          NGRID=25
          NX = MIN(NPAGWD-15,NGRID)
          NY = MIN(NPAGLN-7, NGRID)
      ELSE
          NX = NGRID
          NY = NGRID
      ENDIF
      IF (NX .LT. 11) NX=11
      IF (NY .LT. 11) NY=11
      IF (NX .GE. NXMAX)  NX=NXMAX-1
C         ask if parameter outside limits
      IF (NVARL(KE1) .GT. 1)  THEN
         IF (XLO .LT. ALIM(KE1))  XLO = ALIM(KE1)
         IF (XUP .GT. BLIM(KE1))  XUP = BLIM(KE1)
      ENDIF
      IF (NVARL(KE2) .GT. 1)   THEN
         IF (YLO .LT. ALIM(KE2))  YLO = ALIM(KE2)
         IF (YUP .GT. BLIM(KE2))  YUP = BLIM(KE2)
      ENDIF
      BWIDX = (XUP-XLO)/REAL(NX)
      BWIDY = (YUP-YLO)/REAL(NY)
      IXMID = INT((XSAV-XLO)*REAL(NX)/(XUP-XLO)) + 1
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      DO 185 I= 1, NUMBCS
      CONTUR(I) = AMIN + UP*FLOAT(I-1)**2
  185 CONTINUE
      CONTUR(1) = CONTUR(1) + 0.01*UP
C                fill FCNB to prepare first row, and find column zero
      U(KE2) = YUP
      IXZERO = 0
      XB4 = ONE
      DO 200 IX= 1, NX+1
      U(KE1) = XLO + REAL(IX-1)*BWIDX
      CALL FCN(NPARX,GIN,FF,U,4,FUTIL)
      FCNB(IX) = FF
      IF (XB4.LT.ZERO .AND. U(KE1).GT.ZERO)  IXZERO = IX-1
      XB4 = U(KE1)
      CHMID(IX:IX) = '*'
      CHZERO(IX:IX)= '-'
  200 CONTINUE
      WRITE (ISYSWR,'(A,I3,A,A)') ' Y-AXIS: PARAMETER ',
     +      KE2,': ',CPNAM(KE2)
      IF (IXZERO .GT. 0)  THEN
         CHZERO(IXZERO:IXZERO) = '+'
         CHLN = ' '
         WRITE (ISYSWR,'(12X,A,A)') CHLN(1:IXZERO),'X=0'
      ENDIF
      CALL FLUSH(6)
C                 loop over rows
      DO 280 IY= 1, NY
      UNEXT = U(KE2) - BWIDY
C                 prepare this line's background pattern for contour
      CHLN = ' '
      CHLN(IXMID:IXMID) = '*'
      IF (IXZERO .NE. 0) CHLN(IXZERO:IXZERO) = ':'
      IF (U(KE2).GT.YSAV .AND. UNEXT.LT.YSAV) CHLN=CHMID
      IF (U(KE2).GT.ZERO .AND. UNEXT.LT.ZERO) CHLN=CHZERO
      U(KE2) = UNEXT
      YLABEL = U(KE2) + 0.5*BWIDY
C                 move FCNB to FCNA and fill FCNB with next row
      DO 220 IX= 1, NX+1
      FCNA(IX) = FCNB(IX)
      U(KE1) = XLO + REAL(IX-1)*BWIDX
      CALL FCN(NPARX,GIN,FF,U,4,FUTIL)
      FCNB(IX) = FF
  220 CONTINUE
C                 look for contours crossing the FCNxy squares
      DO 250 IX= 1, NX
      FMX = MAX(FCNA(IX),FCNB(IX),FCNA(IX+1),FCNB(IX+1))
      FMN = MIN(FCNA(IX),FCNB(IX),FCNA(IX+1),FCNB(IX+1))
      DO 230 ICS= 1, NUMBCS
      IF (CONTUR(ICS) .GT. FMN)  GO TO 240
  230 CONTINUE
      GO TO 250
  240 IF (CONTUR(ICS) .LT. FMX) CHLN(IX:IX)=CLABEL(ICS:ICS)
  250 CONTINUE
C                 print a row of the contour plot
      WRITE (ISYSWR,'(1X,G12.4,1X,A)') YLABEL,CHLN(1:NX)
  280 CONTINUE
C                 contours printed, label x-axis
      CHLN = ' '
      CHLN( 1: 1) = 'I'
      CHLN(IXMID:IXMID) = 'I'
      CHLN(NX:NX) = 'I'
      WRITE (ISYSWR,'(14X,A)') CHLN(1:NX)
C                the hardest of all: print x-axis scale!
      CHLN = ' '
      IF (NX .LE. 26) THEN
          NL = MAX(NX-12,2)
          NL2 = NL/2
          WRITE (ISYSWR,'(8X,G12.4,A,G12.4)') XLO,CHLN(1:NL),XUP
          WRITE (ISYSWR,'(14X,A,G12.4)')   CHLN(1:NL2),XSAV
      ELSE
          NL = MAX(NX-24,2)/2
          NL2 = NL
          IF (NL .GT. 10) NL2=NL-6
          WRITE (ISYSWR,'(8X,G12.4,A,G12.4,A,G12.4)')  XLO,
     +      CHLN(1:NL),XSAV,CHLN(1:NL2),XUP
      ENDIF
      WRITE (ISYSWR,'(6X,A,I3,A,A,A,G12.4)') ' X-AXIS: PARAMETER',
     +    KE1,': ',CPNAM(KE1),'  ONE COLUMN=',BWIDX
      WRITE (ISYSWR,'(A,G12.4,A,G12.4,A)') ' FUNCTION VALUES: F(I)=',
     +    AMIN,' +',UP,' *I**2'
C                 finished.  reset input values
      U(KE1) = XSAV
      U(KE2) = YSAV
      IERRF = 0
      CALL FLUSH(6)
      RETURN
 1350 WRITE (ISYSWR,1351)
 1351 FORMAT (' INVALID PARAMETER NUMBER(S) REQUESTED.  IGNORED.' /)
      IERRF = 1
      CALL FLUSH(6)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCOMD(FCN,CRDBIN,ICONDN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called by user.  'Reads' a command string and executes.
CC     Equivalent to MNEXCM except that the command is given as a
CC          character string.
CC
CC     ICONDN = 0: command executed normally
CC              1: command is blank, ignored
CC              2: command line unreadable, ignored
CC              3: unknown command, ignored
CC              4: abnormal termination (e.g., MIGRAD not converged)
CC              5: command is a request to read PARAMETER definitions
CC              6: 'SET INPUT' command
CC              7: 'SET TITLE' command
CC              8: 'SET COVAR' command
CC              9: reserved
CC             10: END command
CC             11: EXIT or STOP command
CC             12: RETURN command
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION PLIST(MAXP)
      CHARACTER COMAND*(MAXCWD)
      CHARACTER CLOWER*26, CUPPER*26
      LOGICAL LEADER
C
      EXTERNAL FCN,FUTIL
      CHARACTER*(*) CRDBIN
      CHARACTER*100 CRDBUF
      DATA CLOWER/'abcdefghijklmnopqrstuvwxyz'/
      DATA CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
      LENBUF = LEN(CRDBIN)
      CRDBUF = CRDBIN
      ICONDN = 0
C     record not case-sensitive, get upper case, strip leading blanks
      LEADER = .TRUE.
      IPOS = 1
         DO 110 I= 1, MIN(MAXCWD,LENBUF)
         IF (CRDBUF(I:I) .EQ. '''') GO TO 111
         IF (CRDBUF(I:I) .EQ. ' ')  THEN
           IF (LEADER) IPOS = IPOS + 1
           GO TO 110
         ENDIF
         LEADER = .FALSE.
           DO 108 IC= 1, 26
           IF (CRDBUF(I:I) .EQ. CLOWER(IC:IC)) CRDBUF(I:I)=CUPPER(IC:IC)
  108      CONTINUE
  110    CONTINUE
  111 CONTINUE
C                     blank or null command
      IF (IPOS .GT. LENBUF)  THEN
         WRITE (ISYSWR,'(A)') ' BLANK COMMAND IGNORED.'
         ICONDN = 1
         GO TO 900
      ENDIF
C                                           . .   preemptive commands
C               if command is 'PARAMETER'
      IF (CRDBUF(IPOS:IPOS+2) .EQ. 'PAR')    THEN
         ICONDN = 5
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               if command is 'SET INPUT'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET INP')  THEN
         ICONDN = 6
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C              if command is 'SET TITLE'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET TIT')  THEN
         ICONDN = 7
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               if command is 'SET COVARIANCE'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET COV')   THEN
         ICONDN = 8
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               crack the command . . . . . . . . . . . . . . . .
      CALL MNCRCK(CRDBUF(IPOS:LENBUF),MAXCWD,COMAND,LNC,
     +                            MAXP,  PLIST, LLIST, IERR,ISYSWR)
      IF (IERR .GT. 0) THEN
            WRITE (ISYSWR,'(A)') ' COMMAND CANNOT BE INTERPRETED'
            ICONDN = 2
            GO TO 900
      ENDIF
C
      CALL MNEXCM(FCN,COMAND(1:LNC),PLIST,LLIST,IERR,FUTIL)
      ICONDN = IERR
  900 CALL FLUSH(6)  
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCONT(FCN,KE1,KE2,NPTU,XPTU,YPTU,IERRF,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Find NPTU points along a contour where the function
CC             FMIN (X(KE1),X(KE2)) =  AMIN+UP
CC       where FMIN is the minimum of FCN with respect to all
CC       the other NPAR-2 variable parameters (if any).
CC   IERRF on return will be equal to the number of points found:
CC     NPTU if normal termination with NPTU points found
CC     -1   if errors in the calling sequence (KE1, KE2 not variable)
CC      0   if less than four points can be found (using MNMNOT)
CC     n>3  if only n points can be found (n < NPTU)
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION XPTU(NPTU), YPTU(NPTU), W(MNI),GCC(MNI)
      CHARACTER CHERE*10
      PARAMETER (CHERE='MNContour ')
      LOGICAL LDEBUG
      EXTERNAL FCN,FUTIL
C                 input arguments: parx, pary, devs, ngrid
      LDEBUG = (IDBG(6) .GE. 1)
      IF (KE1.LE.0 .OR. KE2.LE.0)  GO TO 1350
      IF (KE1.GT.NU .OR. KE2.GT.NU)  GO TO 1350
      KI1 = NIOFEX(KE1)
      KI2 = NIOFEX(KE2)
      IF (KI1.LE.0 .OR. KI2.LE.0)  GO TO 1350
      IF (KI1 .EQ. KI2)  GO TO 1350
      IF (NPTU .LT. 4)  GO TO 1400
C
      NFCNCO = NFCN
      NFCNMX = 100*(NPTU+5)*(NPAR+1)
C           The minimum
      CALL MNCUVE(FCN,FUTIL)
      U1MIN = U(KE1)
      U2MIN = U(KE2)
      IERRF = 0
      CFROM = CHERE
      NFCNFR = NFCNCO
      IF (ISW(5) .GE. 0)  THEN
         WRITE (ISYSWR,'(1X,A,I4,A)')
     +   'START MNCONTOUR CALCULATION OF',NPTU,' POINTS ON CONTOUR.'
         IF (NPAR .GT. 2) THEN
            IF (NPAR .EQ. 3) THEN
              KI3 = 6 - KI1 - KI2
              KE3 = NEXOFI(KI3)
              WRITE (ISYSWR,'(1X,A,I3,2X,A)')
     +        'EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER ',
     +        KE3, CPNAM(KE3)
            ELSE
              WRITE (ISYSWR,'(1X,A,I3,A)')
     +        'EACH POINT IS A MINIMUM WITH RESPECT TO THE OTHER',
     +        NPAR-2, ' VARIABLE PARAMETERS.'
            ENDIF
         ENDIF
      ENDIF
      CALL FLUSH(6)
C
C           Find the first four points using MNMNOT
C              ........................ first two points
      CALL MNMNOT(FCN,KE1,KE2,VAL2PL,VAL2MI,FUTIL)
      IF (ERN(KI1) .EQ. UNDEFI)  THEN
         XPTU(1) = ALIM(KE1)
         CALL MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
      ELSE
         IF (ERN(KI1) .GE. ZERO)  GO TO 1500
         XPTU(1) = U1MIN+ERN(KI1)
      ENDIF
      YPTU(1) = VAL2MI
C
      IF (ERP(KI1) .EQ. UNDEFI)  THEN
         XPTU(3) = BLIM(KE1)
         CALL MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
      ELSE
         IF (ERP(KI1) .LE. ZERO)  GO TO 1500
         XPTU(3) = U1MIN+ERP(KI1)
      ENDIF
      YPTU(3) = VAL2PL
      SCALX = 1.0/(XPTU(3) - XPTU(1))
C              ........................... next two points
      CALL MNMNOT(FCN,KE2,KE1,VAL2PL,VAL2MI,FUTIL)
      IF (ERN(KI2) .EQ. UNDEFI)  THEN
         YPTU(2) = ALIM(KE2)
         CALL MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
      ELSE
         IF (ERN(KI2) .GE. ZERO)  GO TO 1500
         YPTU(2) = U2MIN+ERN(KI2)
      ENDIF
      XPTU(2) = VAL2MI
      IF (ERP(KI2) .EQ. UNDEFI)  THEN
         YPTU(4) = BLIM(KE2)
         CALL MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
      ELSE
         IF (ERP(KI2) .LE. ZERO)  GO TO 1500
         YPTU(4) = U2MIN+ERP(KI2)
      ENDIF
      XPTU(4) = VAL2PL
      SCALY = 1.0/(YPTU(4) - YPTU(2))
      NOWPTS = 4
      NEXT = 5
      IF (LDEBUG) THEN
         WRITE (ISYSWR,'(A)') ' Plot of four points found by MINOS'
         XPT(1) = U1MIN
         YPT(1) = U2MIN
         CHPT(1) = ' '
         NALL = MIN(NOWPTS+1,MAXCPT)
         DO 85 I= 2, NALL
           XPT(I) = XPTU(I-1)
           YPT(I) = YPTU(I-1)
   85    CONTINUE
           CHPT(2)= 'A'
           CHPT(3)= 'B'
           CHPT(4)= 'C'
           CHPT(5)= 'D'
         CALL MNPLOT(XPT,YPT,CHPT,NALL,ISYSWR,NPAGWD,NPAGLN)
      ENDIF
C
C               ..................... save some values before fixing
      ISW2 = ISW(2)
      ISW4 = ISW(4)
      SIGSAV = EDM
      ISTRAV = ISTRAT
      DC = DCOVAR
      APSI  = EPSI*0.5
      ABEST=AMIN
      MPAR=NPAR
      NFMXIN = NFCNMX
      DO 125 I= 1, MPAR
  125 XT(I) = X(I)
      DO 130 J= 1, MPAR*(MPAR+1)/2
  130 VTHMAT(J) = VHMAT(J)
      DO 135 I= 1, MPAR
      GCC(I) = GLOBCC(I)
  135 W(I) = WERR(I)
C                           fix the two parameters in question
      KINTS = NIOFEX(KE1)
      CALL MNFIXP (KINTS,IERR)
      KINTS = NIOFEX(KE2)
      CALL MNFIXP (KINTS,IERR)
C               ......................Fill in the rest of the points
      DO 900 INEW= NEXT, NPTU
C            find the two neighbouring points with largest separation
      BIGDIS = 0.
         DO 200  IOLD = 1, INEW-1
         I2 = IOLD + 1
         IF (I2 .EQ. INEW) I2 = 1
         DIST = (SCALX*(XPTU(IOLD)-XPTU(I2)))**2 +
     +          (SCALY*(YPTU(IOLD)-YPTU(I2)))**2
         IF (DIST .GT. BIGDIS) THEN
            BIGDIS = DIST
            IDIST = IOLD
         ENDIF
  200    CONTINUE
      I1 = IDIST
      I2 = I1 + 1
      IF (I2 .EQ. INEW) I2 = 1
C                   next point goes between I1 and I2
      A1 = HALF
      A2 = HALF
  300 XMIDCR = A1*XPTU(I1) + A2*XPTU(I2)
      YMIDCR = A1*YPTU(I1) + A2*YPTU(I2)
      XDIR = YPTU(I2) - YPTU(I1)
      YDIR = XPTU(I1) - XPTU(I2)
      SCLFAC = MAX(ABS(XDIR*SCALX), ABS(YDIR*SCALY))
      XDIRCR = XDIR/SCLFAC
      YDIRCR = YDIR/SCLFAC
      KE1CR = KE1
      KE2CR = KE2
C                Find the contour crossing point along DIR
      AMIN = ABEST
      CALL MNCROS(FCN,AOPT,IERCR,FUTIL)
      IF (IERCR .GT. 1)  THEN
C              If cannot find mid-point, try closer to point 1
         IF (A1 .GT. HALF) THEN
            IF (ISW(5) .GE. 0)
     +      WRITE (ISYSWR,'(A,A,I3,A)') ' MNCONT CANNOT FIND NEXT',
     +           ' POINT ON CONTOUR.  ONLY ',NOWPTS,' POINTS FOUND.'
            GO TO 950
         ENDIF
         CALL MNWARN('W',CHERE,'Cannot find midpoint, try closer.')
         A1 = 0.75
         A2 = 0.25
         GO TO 300
      ENDIF
C                Contour has been located, insert new point in list
         DO 830 MOVE= NOWPTS,I1+1,-1
         XPTU(MOVE+1) = XPTU(MOVE)
         YPTU(MOVE+1) = YPTU(MOVE)
  830    CONTINUE
      NOWPTS = NOWPTS + 1
      XPTU(I1+1) = XMIDCR + XDIRCR*AOPT
      YPTU(I1+1) = YMIDCR + YDIRCR*AOPT
  900 CONTINUE
  950 CONTINUE
C
      IERRF = NOWPTS
      CSTATU = 'SUCCESSFUL'
      IF (NOWPTS .LT. NPTU)  CSTATU = 'INCOMPLETE'
C                make a lineprinter plot of the contour
      IF (ISW(5) .GE. 0) THEN
         XPT(1) = U1MIN
         YPT(1) = U2MIN
         CHPT(1) = ' '
         NALL = MIN(NOWPTS+1,MAXCPT)
         DO 1000 I= 2, NALL
           XPT(I) = XPTU(I-1)
           YPT(I) = YPTU(I-1)
           CHPT(I)= 'X'
 1000    CONTINUE
         WRITE (ISYSWR,'(A,I3,2X,A)') ' Y-AXIS: PARAMETER ',KE2,
     +        CPNAM(KE2)
         CALL MNPLOT(XPT,YPT,CHPT,NALL,ISYSWR,NPAGWD,NPAGLN)
         WRITE (ISYSWR,'(25X,A,I3,2X,A)') 'X-AXIS: PARAMETER ',
     +         KE1,CPNAM(KE1)
      ENDIF
C                 print out the coordinates around the contour
      IF (ISW(5) .GE. 1)  THEN
         NPCOL = (NOWPTS+1)/2
         NFCOL = NOWPTS/2
         WRITE (ISYSWR,'(/I5,A,G13.5,A,G11.3)') NOWPTS,
     +    ' POINTS ON CONTOUR.   FMIN=',ABEST,'   ERRDEF=',UP
         WRITE (ISYSWR,'(9X,A,3X,A,18X,A,3X,A)')
     +         CPNAM(KE1),CPNAM(KE2),CPNAM(KE1),CPNAM(KE2)
         DO 1050 LINE = 1, NFCOL
           LR = LINE + NPCOL
           WRITE (ISYSWR,'(1X,I5,2G13.5,10X,I5,2G13.5)')
     +     LINE,XPTU(LINE),YPTU(LINE),LR,XPTU(LR),YPTU(LR)
 1050    CONTINUE
         IF (NFCOL .LT. NPCOL) WRITE (ISYSWR,'(1X,I5,2G13.5)')
     +                         NPCOL,XPTU(NPCOL),YPTU(NPCOL)
      ENDIF
C                                    . . contour finished. reset v
      ITAUR = 1
      CALL MNFREE(1)
      CALL MNFREE(1)
      DO 1100 J= 1, MPAR*(MPAR+1)/2
 1100 VHMAT(J) = VTHMAT(J)
      DO 1120 I= 1, MPAR
      GLOBCC(I) = GCC(I)
      WERR(I) = W(I)
 1120 X(I) = XT(I)
      CALL MNINEX (X)
      EDM = SIGSAV
      AMIN = ABEST
      ISW(2) = ISW2
      ISW(4) = ISW4
      DCOVAR = DC
      ITAUR = 0
      NFCNMX = NFMXIN
      ISTRAT = ISTRAV
      U(KE1) = U1MIN
      U(KE2) = U2MIN
      GO TO 2000
C                                     Error returns
 1350 WRITE (ISYSWR,'(A)') ' INVALID PARAMETER NUMBERS.'
      GO TO 1450
 1400 WRITE (ISYSWR,'(A)') ' LESS THAN FOUR POINTS REQUESTED.'
 1450 IERRF = -1
      CSTATU = 'USER ERROR'
      GO TO 2000
 1500 WRITE (ISYSWR,'(A)') ' MNCONT UNABLE TO FIND FOUR POINTS.'
      U(KE1) = U1MIN
      U(KE2) = U2MIN
      IERRF = 0
      CSTATU = 'FAILED'
 2000 CALL FLUSH(6)
      CONTINUE
      CFROM = CHERE
      NFCNFR = NFCNCO
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCRCK(CRDBUF,MAXCWD,COMAND,LNC,
     +                         MXP,   PLIST, LLIST,IERR,ISYSWR)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC
CC       Called from MNREAD.
CC       Cracks the free-format input, expecting zero or more
CC         alphanumeric fields (which it joins into COMAND(1:LNC))
CC         followed by one or more numeric fields separated by
CC         blanks and/or one comma.  The numeric fields are put into
CC         the LLIST (but at most MXP) elements of PLIST.
CC      IERR = 0 if no errors,
CC           = 1 if error(s).
CC      Diagnostic messages are written to ISYSWR
CC
      PARAMETER (MAXELM=25, MXLNEL=19)
      CHARACTER*(*) COMAND, CRDBUF
      CHARACTER CNUMER*13, CELMNT(MAXELM)*(MXLNEL), CNULL*15
      DIMENSION LELMNT(MAXELM),PLIST(MXP)
      DATA CNULL /')NULL STRING   '/
      DATA CNUMER/'123456789-.0+'/
      IELMNT = 0
      LEND = LEN(CRDBUF)
      NEXTB = 1
      IERR = 0
C                                   . . . .  loop over words CELMNT
   10 CONTINUE
      DO 100 IPOS= NEXTB,LEND
         IBEGIN = IPOS
         IF (CRDBUF(IPOS:IPOS).EQ.' ')  GO TO 100
         IF (CRDBUF(IPOS:IPOS).EQ.',')  GO TO 250
         GO TO 150
  100 CONTINUE
         GO TO 300
  150 CONTINUE
C               found beginning of word, look for end
         DO 180 IPOS = IBEGIN+1,LEND
         IF (CRDBUF(IPOS:IPOS).EQ.' ')  GO TO 250
         IF (CRDBUF(IPOS:IPOS).EQ.',')  GO TO 250
  180    CONTINUE
      IPOS = LEND+1
  250 IEND = IPOS-1
      IELMNT = IELMNT + 1
      IF (IEND .GE. IBEGIN) THEN
         CELMNT(IELMNT) = CRDBUF(IBEGIN:IEND)
      ELSE
         CELMNT(IELMNT) = CNULL
      ENDIF
      LELMNT(IELMNT) = IEND-IBEGIN+1
      IF (LELMNT(IELMNT) .GT. MXLNEL)  THEN
         WRITE (ISYSWR, 253) CRDBUF(IBEGIN:IEND),CELMNT(IELMNT)
  253    FORMAT (' MINUIT WARNING: INPUT DATA WORD TOO LONG.'
     +   /'     ORIGINAL:',A
     +   /' TRUNCATED TO:',A)
         LELMNT(IELMNT) = MXLNEL
         ENDIF
      IF (IPOS .GE. LEND) GO TO 300
      IF (IELMNT .GE. MAXELM)  GO TO 300
C                     look for comma or beginning of next word
         DO 280 IPOS= IEND+1,LEND
         IF (CRDBUF(IPOS:IPOS) .EQ. ' ') GO TO 280
         NEXTB = IPOS
         IF (CRDBUF(IPOS:IPOS) .EQ. ',') NEXTB = IPOS+1
         GO TO 10
  280    CONTINUE
C                 All elements found, join the alphabetic ones to
C                                form a command
  300 CONTINUE
      NELMNT = IELMNT
      COMAND = ' '
      LNC = 1
      PLIST(1) = 0.
      LLIST = 0
      IF (IELMNT .EQ. 0)  GO TO 900
      KCMND = 0
         DO 400 IELMNT = 1, NELMNT
         IF (CELMNT(IELMNT) .EQ. CNULL)  GO TO 450
            DO 350 IC= 1, 13
            IF (CELMNT(IELMNT)(1:1) .EQ. CNUMER(IC:IC)) GO TO 450
  350       CONTINUE
         IF (KCMND .GE. MAXCWD) GO TO 400
         LEFT = MAXCWD-KCMND
         LTOADD = LELMNT(IELMNT)
         IF (LTOADD .GT. LEFT) LTOADD=LEFT
         COMAND(KCMND+1:KCMND+LTOADD) = CELMNT(IELMNT)(1:LTOADD)
         KCMND = KCMND + LTOADD
         IF (KCMND .EQ. MAXCWD)  GO TO 400
         KCMND = KCMND + 1
         COMAND(KCMND:KCMND) = ' '
  400    CONTINUE
      LNC = KCMND
      GO TO 900
  450 CONTINUE
      LNC = KCMND
C                      . . . .  we have come to a numeric field
      LLIST = 0
      DO 600 IFLD= IELMNT,NELMNT
      LLIST = LLIST + 1
      IF (LLIST .GT. MXP) THEN
         NREQ = NELMNT-IELMNT+1
         WRITE (ISYSWR,511) NREQ,MXP
  511 FORMAT (/' MINUIT WARNING IN MNCRCK: '/ ' COMMAND HAS INPUT',I5,
     + ' NUMERIC FIELDS, BUT MINUIT CAN ACCEPT ONLY',I3)
          CALL FLUSH(6)
         GO TO 900
      ENDIF
      IF (CELMNT(IFLD) .EQ. CNULL)  THEN
          PLIST(LLIST) = 0.
        ELSE
          READ (CELMNT(IFLD), '(BN,F19.0)',ERR=575) PLIST(LLIST)
      ENDIF
      GO TO 600
  575 WRITE (ISYSWR,'(A,A,A)') ' FORMAT ERROR IN NUMERIC FIELD: "',
     + CELMNT(IFLD)(1:LELMNT(IFLD)),'"'
      IERR = 1
      PLIST(LLIST) = 0.
  600 CONTINUE
C                                  end loop over numeric fields
  900 CONTINUE
      IF (LNC .LE. 0)  LNC=1
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCROS(FCN,AOPT,IERCR,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Find point where MNEVAL=AMIN+UP, along the line through
CC       XMIDCR,YMIDCR with direction XDIRCR,YDIRCR,   where X and Y 
CC       are parameters KE1CR and KE2CR.  If KE2CR=0 (from MINOS),
CC       only KE1CR is varied.  From MNCONT, both are varied.
CC       Crossing point is at
CC        (U(KE1),U(KE2)) = (XMID,YMID) + AOPT*(XDIR,YDIR)
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER CHERE*10, CHARAL*28, CHSIGN*4
      PARAMETER (CHERE='MNCROS    ', MLSB=3, MAXITR=15, TLR=0.01)
      DIMENSION FLSB(MLSB),ALSB(MLSB), COEFF(3)
      LOGICAL LDEBUG
      EXTERNAL FCN,FUTIL
      DATA  CHARAL/' .ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      LDEBUG = (IDBG(6) .GE. 1)
      AMINSV = AMIN
C        convergence when F is within TLF of AIM and next prediction
C        of AOPT is within TLA of previous value of AOPT
      AIM = AMIN + UP
      TLF = TLR*UP
      TLA = TLR
      XPT(1) = 0.0
      YPT(1) = AIM
      CHPT(1) = ' '
      IPT = 1
      IF (KE2CR .EQ. 0) THEN
        XPT(2) = -1.0
        YPT(2) = AMIN
        CHPT(2) = '.'
        IPT = 2
      ENDIF
C                    find the largest allowed A
      AULIM = 100.
      DO 100 IK= 1, 2
         IF (IK .EQ. 1)  THEN
            KEX = KE1CR
            ZMID = XMIDCR
            ZDIR = XDIRCR
         ELSE
            IF (KE2CR .EQ. 0)  GO TO 100
            KEX = KE2CR
            ZMID = YMIDCR
            ZDIR = YDIRCR
         ENDIF
         IF (NVARL(KEX) .LE. 1) GO TO 100
         IF (ZDIR .EQ. ZERO)      GO TO 100
         ZLIM = ALIM(KEX)
         IF (ZDIR .GT. ZERO) ZLIM = BLIM(KEX)
         AULIM = MIN(AULIM,(ZLIM-ZMID)/ZDIR)
  100 CONTINUE
C                  LSB = Line Search Buffer
C          first point
      ANEXT = 0.
      AOPT = ANEXT
      LIMSET = .FALSE.
        IF (AULIM .LT. AOPT+TLA)  LIMSET = .TRUE.
      CALL MNEVAL(FCN,ANEXT,FNEXT,IEREV,FUTIL)
C debug printout:
      IF (LDEBUG) WRITE (ISYSWR,'(A,I8,A,F10.5,A,2F10.5)')
     + ' MNCROS: calls=',NFCN,'   AIM=',AIM,'  F,A=',FNEXT,AOPT
      IF (IEREV .GT. 0)  GO TO 900
      IF (LIMSET .AND. FNEXT .LE. AIM)  GO TO 930
      IPT = IPT + 1
      XPT(IPT) = ANEXT
      YPT(IPT) = FNEXT
      CHPT(IPT)= CHARAL(IPT:IPT)
      ALSB(1) = ANEXT
      FLSB(1) = FNEXT
      FNEXT = MAX(FNEXT,AMINSV+0.1*UP)
      AOPT =  SQRT((UP)/(FNEXT-AMINSV)) - 1.0
      IF (ABS(FNEXT-AIM) .LT. TLF)  GO TO 800
C
      IF (AOPT .LT. -HALF)  AOPT = -HALF
      IF (AOPT .GT. ONE)    AOPT = ONE
      LIMSET = .FALSE.
      IF (AOPT .GT. AULIM)  THEN
              AOPT = AULIM
              LIMSET = .TRUE.
      ENDIF
      CALL MNEVAL(FCN,AOPT,FNEXT,IEREV,FUTIL)
C debug printout:
      IF (LDEBUG) WRITE (ISYSWR,'(A,I8,A,F10.5,A,2F10.5)')
     + ' MNCROS: calls=',NFCN,'   AIM=',AIM,'  F,A=',FNEXT,AOPT
      IF (IEREV .GT. 0)  GO TO 900
      IF (LIMSET .AND. FNEXT .LE. AIM)  GO TO 930
      ALSB(2) = AOPT
      IPT = IPT + 1
      XPT(IPT) = ALSB(2)
      YPT(IPT) = FNEXT
      CHPT(IPT)= CHARAL(IPT:IPT)
      FLSB(2) = FNEXT
      DFDA = (FLSB(2)-FLSB(1))/ (ALSB(2)-ALSB(1))
C                   DFDA must be positive on the contour
      IF (DFDA .GT. ZERO)  GO TO 460
  300    CALL MNWARN('D',CHERE,'Looking for slope of the right sign')
         MAXLK = MAXITR - IPT
         DO 400 IT= 1, MAXLK
            ALSB(1) = ALSB(2)
            FLSB(1) = FLSB(2)
            AOPT = ALSB(1) + 0.2*REAL(IT)
            LIMSET = .FALSE.
            IF (AOPT .GT. AULIM)  THEN
              AOPT = AULIM
              LIMSET = .TRUE.
            ENDIF
            CALL MNEVAL(FCN,AOPT,FNEXT,IEREV,FUTIL)
C debug printout:
      IF (LDEBUG) WRITE (ISYSWR,'(A,I8,A,F10.5,A,2F10.5)')
     + ' MNCROS: calls=',NFCN,'   AIM=',AIM,'  F,A=',FNEXT,AOPT
            IF (IEREV .GT. 0)  GO TO 900
            IF (LIMSET .AND. FNEXT .LE. AIM)  GO TO 930
               ALSB(2) = AOPT
               IPT = IPT + 1
               XPT(IPT) = ALSB(2)
               YPT(IPT) = FNEXT
               CHPT(IPT)= CHARAL(IPT:IPT)
            FLSB(2) = FNEXT
            DFDA = (FLSB(2)-FLSB(1))/ (ALSB(2)-ALSB(1))
            IF (DFDA .GT. ZERO)  GO TO 450
  400    CONTINUE
         CALL MNWARN('W',CHERE,'Cannot find slope of the right sign')
         GO TO 950
  450    CONTINUE
C                    we have two points with the right slope
  460 AOPT = ALSB(2) + (AIM-FLSB(2))/DFDA
      FDIST = MIN(ABS(AIM -FLSB(1)),ABS(AIM -FLSB(2)))
      ADIST = MIN(ABS(AOPT-ALSB(1)),ABS(AOPT-ALSB(2)))
      TLA = TLR
      IF (ABS(AOPT) .GT. ONE)  TLA = TLR*ABS(AOPT)
      IF (ADIST .LT. TLA .AND. FDIST .LT. TLF) GO TO 800
      IF (IPT .GE. MAXITR)  GO TO 950
      BMIN = MIN(ALSB(1),ALSB(2)) - 1.0
      IF (AOPT .LT. BMIN)  AOPT = BMIN
      BMAX = MAX(ALSB(1),ALSB(2)) + 1.0
      IF (AOPT .GT. BMAX)  AOPT = BMAX
C                    Try a third point
      LIMSET = .FALSE.
      IF (AOPT .GT. AULIM) THEN
         AOPT = AULIM
         LIMSET = .TRUE.
      ENDIF
      CALL MNEVAL(FCN,AOPT,FNEXT,IEREV,FUTIL)
C debug printout:
      IF (LDEBUG) WRITE (ISYSWR,'(A,I8,A,F10.5,A,2F10.5)')
     + ' MNCROS: calls=',NFCN,'   AIM=',AIM,'  F,A=',FNEXT,AOPT
      IF (IEREV .GT. 0)  GO TO 900
      IF (LIMSET .AND. FNEXT .LE. AIM)  GO TO 930
      ALSB(3) = AOPT
      IPT = IPT + 1
      XPT(IPT) = ALSB(3)
      YPT(IPT) = FNEXT
      CHPT(IPT)= CHARAL(IPT:IPT)
      FLSB(3) = FNEXT
      INEW = 3
C                now we have three points, ask how many <AIM
      ECARMN = ABS(FNEXT-AIM)
      IBEST = 3
      ECARMX = 0.
      NOLESS = 0
      DO 480 I= 1, 3
         ECART = ABS(FLSB(I) - AIM)
         IF (ECART .GT. ECARMX) THEN
            ECARMX = ECART
            IWORST = I
         ENDIF
         IF (ECART .LT. ECARMN) THEN
            ECARMN = ECART
            IBEST = I
         ENDIF
         IF (FLSB(I) .LT. AIM) NOLESS = NOLESS + 1
  480 CONTINUE
      INEW = IBEST
C           if at least one on each side of AIM, fit a parabola
      IF (NOLESS.EQ.1 .OR. NOLESS.EQ.2) GO TO 500
C           if all three are above AIM, third must be closest to AIM
      IF (NOLESS .EQ. 0 .AND. IBEST .NE. 3)  GO TO 950
C           if all three below, and third is not best, then slope
C             has again gone negative, look for positive slope.
      IF (NOLESS .EQ. 3 .AND. IBEST .NE. 3) THEN
          ALSB(2) = ALSB(3)
          FLSB(2) = FLSB(3)
          GO TO 300
      ENDIF
C           in other cases, new straight line thru last two points
      ALSB(IWORST) = ALSB(3)
      FLSB(IWORST) = FLSB(3)
      DFDA = (FLSB(2)-FLSB(1))/ (ALSB(2)-ALSB(1))
      GO TO 460
C                parabola fit
  500 CALL MNPFIT(ALSB,FLSB,3,COEFF,SDEV)
      IF (COEFF(3) .LE. ZERO)  CALL MNWARN ('D',CHERE,
     +             'Curvature is negative near contour line.')
      DETERM =  COEFF(2)**2 - 4.*COEFF(3)*(COEFF(1)-AIM)
      IF (DETERM .LE. ZERO)   THEN
          CALL MNWARN('D',CHERE,'Problem 2, impossible determinant')
          GO TO 950
      ENDIF
C                Find which root is the right one
      RT = SQRT(DETERM)
      X1 = (-COEFF(2) + RT)/(2.*COEFF(3))
      X2 = (-COEFF(2) - RT)/(2.*COEFF(3))
      S1 = COEFF(2) + 2.*X1*COEFF(3)
      S2 = COEFF(2) + 2.*X2*COEFF(3)
      IF (S1*S2 .GT. ZERO) WRITE (ISYSWR,'(A)') ' MNCONTour problem 1'
      AOPT = X1
      SLOPE = S1
      IF (S2 .GT. ZERO)  THEN
         AOPT = X2
         SLOPE = S2
      ENDIF
C         ask if converged
      TLA = TLR
      IF (ABS(AOPT) .GT. ONE)  TLA = TLR*ABS(AOPT)
      IF (ABS(AOPT-ALSB(IBEST)) .LT. TLA  .AND. 
     &    ABS(FLSB(IBEST)-AIM)  .LT. TLF)  GO TO 800
      IF (IPT .GE. MAXITR)  GO TO 950
C         see if proposed point is in acceptable zone between L and R
C         first find ILEFT, IRIGHT, IOUT and IBEST
      ILEFT = 0
      IRIGHT = 0
      IBEST = 1
      ECARMX = 0.
      ECARMN = ABS(AIM-FLSB(1))
      DO 550 I= 1, 3
      ECART = ABS(FLSB(I) - AIM)
      IF (ECART .LT. ECARMN) THEN
         ECARMN = ECART
         IBEST = I
      ENDIF
      IF (ECART .GT. ECARMX) ECARMX = ECART
      IF (FLSB(I) .GT. AIM)  THEN
         IF (IRIGHT .EQ. 0)  THEN
            IRIGHT = I
         ELSE IF (FLSB(I) .GT. FLSB(IRIGHT)) THEN
            IOUT = I
         ELSE
            IOUT = IRIGHT
            IRIGHT = I
         ENDIF
      ELSE IF (ILEFT .EQ. 0)  THEN
         ILEFT = I
      ELSE IF (FLSB(I) .LT. FLSB(ILEFT)) THEN
         IOUT = I
      ELSE     
         IOUT = ILEFT
         ILEFT = I
      ENDIF
  550 CONTINUE 
C       avoid keeping a very bad point next time around
      IF (ECARMX .GT. 10.*ABS(FLSB(IOUT)-AIM))
     &    AOPT = HALF*AOPT + HALF*HALF*(ALSB(IRIGHT)+ALSB(ILEFT))      
C         knowing ILEFT and IRIGHT, get acceptable window
      SMALLA = 0.1*TLA
      IF (SLOPE*SMALLA .GT. TLF)  SMALLA = TLF/SLOPE
      ALEFT  = ALSB(ILEFT)  + SMALLA
      ARIGHT = ALSB(IRIGHT) - SMALLA
C         move proposed point AOPT into window if necessary
      IF (AOPT .LT. ALEFT)  AOPT = ALEFT
      IF (AOPT .GT. ARIGHT) AOPT = ARIGHT
      IF (ALEFT .GT. ARIGHT)AOPT = HALF*(ALEFT + ARIGHT)
C         see if proposed point outside limits (should be impossible!)
      LIMSET = .FALSE.
      IF (AOPT .GT. AULIM)  THEN
              AOPT = AULIM
              LIMSET = .TRUE.
      ENDIF
C                  Evaluate function at new point AOPT
      CALL MNEVAL(FCN,AOPT,FNEXT,IEREV,FUTIL)
C debug printout:
      IF (LDEBUG) WRITE (ISYSWR,'(A,I8,A,F10.5,A,2F10.5)')
     + ' MNCROS: calls=',NFCN,'   AIM=',AIM,'  F,A=',FNEXT,AOPT
      IF (IEREV .GT. 0)  GO TO 900
      IF (LIMSET .AND. FNEXT .LE. AIM)  GO TO 930
      IPT = IPT + 1
      XPT(IPT) = AOPT
      YPT(IPT) = FNEXT
      CHPT(IPT)= CHARAL(IPT:IPT)
C                Replace odd point by new one
      ALSB(IOUT) = AOPT
      FLSB(IOUT) = FNEXT
C          the new point may not be the best, but it is the only one
C          which could be good enough to pass convergence criteria
      IBEST = IOUT
      GO TO 500
C
C       Contour has been located, return point to MNCONT OR MINOS
  800 CONTINUE
      IERCR = 0
      GO TO 1000
C                error in the minimization
  900 IF (IEREV .EQ. 1)  GO TO 940
      GO TO 950
C                parameter up against limit
  930 IERCR = 1
      GO TO 1000
C                too many calls to FCN
  940 IERCR = 2
      GO TO 1000
C                cannot find next point
  950 IERCR = 3
C                in any case
 1000 CONTINUE
      IF (LDEBUG) THEN
         ITOOHI = 0
         DO 1100 I= 1, IPT
         IF (YPT(I) .GT. AIM+UP) THEN
            YPT(I) = AIM+UP
            CHPT(I) = '+'
            ITOOHI = 1
         ENDIF
 1100    CONTINUE
         CHSIGN = 'POSI'
         IF (XDIRCR .LT. ZERO)  CHSIGN = 'NEGA'
         IF (KE2CR .EQ. 0)  WRITE (ISYSWR, '(2X,A,A,I3)')
     +            CHSIGN,'TIVE MINOS ERROR, PARAMETER ',KE1CR
         IF (ITOOHI .EQ. 1)  WRITE (ISYSWR, '(10X,A)')
     +            'POINTS LABELLED "+" WERE TOO HIGH TO PLOT.'
         IF (IERCR .EQ. 1) WRITE (ISYSWR,'(10X,A)')
     +            'RIGHTMOST POINT IS UP AGAINST LIMIT.'
         CALL MNPLOT(XPT,YPT,CHPT,IPT,ISYSWR,NPAGWD,NPAGLN)
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCUVE(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Makes sure that the current point is a local
CC        minimum and that the error matrix exists,
CC        or at least something good enough for MINOS and MNCONT
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      IF (ISW(4) .LT. 1) THEN
          WRITE (ISYSWR,'(/A,A)')
     +    ' FUNCTION MUST BE MINIMIZED BEFORE CALLING ',CFROM
          APSI = EPSI
          CALL MNMIGR(FCN,FUTIL)
      ENDIF
      IF (ISW(2) .LT. 3)  THEN
         CALL MNHESS(FCN,FUTIL)
         IF (ISW(2) .LT. 1)  THEN
            CALL MNWARN('W',CFROM,'NO ERROR MATRIX.  WILL IMPROVISE.')
            DO 555 I=1,NPAR
              NDEX = I*(I-1)/2
              DO 554 J=1,I-1
              NDEX = NDEX + 1
  554         VHMAT(NDEX) = 0.
            NDEX = NDEX + 1
            IF (G2(I) .LE. ZERO)  THEN
              WINT = WERR(I)
              IEXT = NEXOFI(I)
              IF (NVARL(IEXT) .GT. 1) THEN
                 CALL MNDXDI(X(I),I,DXDI)
                 IF (ABS(DXDI) .LT. .001) THEN
                    WINT = .01
                 ELSE
                    WINT = WINT/ABS(DXDI)
                 ENDIF
              ENDIF
              G2(I) = UP/WINT**2
            ENDIF
            VHMAT(NDEX) = 2./G2(I)
  555       CONTINUE
            ISW(2) = 1
            DCOVAR = 1.
         ELSE
           CALL MNWERR
         ENDIF
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:43  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNDERI(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Calculates the first derivatives of FCN (GRD),
CC        either by finite differences or by transforming the user-
CC        supplied derivatives to internal coordinates,
CC        according to whether ISW(3) is zero or one.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      LOGICAL LDEBUG
      CHARACTER CBF1*22
      NPARX = NPAR
      LDEBUG = (IDBG(2) .GE. 1)
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      IF (ISW(3) .EQ. 1)  GO TO 100
      IF (LDEBUG) THEN
C                       make sure starting at the right place
        CALL MNINEX(X)
        NPARX = NPAR
        CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
        NFCN = NFCN + 1
        IF (FS1 .NE. AMIN) THEN
           DF = AMIN - FS1
           WRITE (CBF1(1:12),'(G12.3)') DF
           CALL MNWARN('D','MNDERI',
     +         'function value differs from AMIN by '//CBF1(1:12) )
           AMIN = FS1
        ENDIF
          WRITE
     +   (ISYSWR,'(/''  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI''/
     +   '' PAR    DERIV     STEP      MINSTEP   OPTSTEP '',
     +   '' D1-D2    2ND DRV'')')
      ENDIF
      DFMIN = 8. * EPSMA2*(ABS(AMIN)+UP)
      VRYSML = 8.* EPSMAC**2
      IF (ISTRAT .LE. 0) THEN
         NCYC = 2
         TLRSTP = 0.5
         TLRGRD = 0.1
      ELSE IF (ISTRAT .EQ. 1) THEN
         NCYC = 3
         TLRSTP = 0.3
         TLRGRD = 0.05
      ELSE
         NCYC = 5
         TLRSTP = 0.1
         TLRGRD = 0.02
      ENDIF
C                                loop over variable parameters
      DO 60  I=1,NPAR
      EPSPRI = EPSMA2 + ABS(GRD(I)*EPSMA2)
C         two-point derivatives always assumed necessary
C         maximum number of cycles over step size depends on strategy
      XTF = X(I)
      STEPB4 = 0.
C                               loop as little as possible here!
      DO 45 ICYC= 1, NCYC
C                 ........ theoretically best step
      OPTSTP = SQRT(DFMIN/(ABS(G2(I))+EPSPRI))
C                     step cannot decrease by more than a factor of ten
      STEP = MAX(OPTSTP, ABS(0.1*GSTEP(I)))
C                 but if parameter has limits, max step size = 0.5
      IF (GSTEP(I).LT.ZERO .AND. STEP.GT.0.5)  STEP=0.5
C                 and not more than ten times the previous step
      STPMAX = 10.*ABS(GSTEP(I))
      IF (STEP .GT. STPMAX)  STEP = STPMAX
C                 minimum step size allowed by machine precision
      STPMIN = MAX(VRYSML, 8.*ABS(EPSMA2*X(I)))
      IF (STEP .LT. STPMIN)  STEP = STPMIN
C                 end of iterations if step change less than factor 2
      IF (ABS((STEP-STEPB4)/STEP) .LT. TLRSTP)  GO TO 50
C         take step positive
      GSTEP(I) = SIGN(STEP, GSTEP(I))
      STEPB4 = STEP
      X(I) = XTF + STEP
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN=NFCN+1
C         take step negative
      X(I) = XTF - STEP
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
      NFCN=NFCN+1
      GRBFOR = GRD(I)
      GRD(I) = (FS1-FS2)/(2.0*STEP)
      G2(I) = (FS1+FS2-2.0*AMIN)/(STEP**2)
      X(I) = XTF
      IF (LDEBUG) THEN
         D1D2 = (FS1+FS2-2.0*AMIN)/STEP
         WRITE (ISYSWR,41) I,GRD(I),STEP,STPMIN,OPTSTP,D1D2,G2(I)
   41    FORMAT (I4,2G11.3,5G10.2)
      ENDIF
C         see if another iteration is necessary
      IF (ABS(GRBFOR-GRD(I))/(ABS(GRD(I))+DFMIN/STEP) .LT. TLRGRD)
     +        GO TO 50
   45 CONTINUE
C                           end of ICYC loop. too many iterations
      IF (NCYC .EQ. 1)  GO TO 50
         WRITE (CBF1,'(2E11.3)')  GRD(I),GRBFOR
         CALL MNWARN('D','MNDERI',
     +         'First derivative not converged. '//CBF1)
   50 CONTINUE
C
   60 CONTINUE
      CALL MNINEX(X)
      RETURN
C                                        .  derivatives calc by fcn
  100 DO 150 IINT= 1, NPAR
      IEXT = NEXOFI(IINT)
      IF (NVARL(IEXT) .GT. 1)  GO TO 120
      GRD(IINT) = GIN(IEXT)
      GO TO 150
  120 DD = (BLIM(IEXT)-ALIM(IEXT))*0.5 *COS(X(IINT))
      GRD(IINT) = GIN(IEXT)*DD
  150 CONTINUE
  200 RETURN
      END
      SUBROUTINE MNDUMMY()
C
C This routine is just to make sure that $MINUIT_DIR
C is used in linking. FOr a while, both MINUIT and
C CERNLIB have the minuit code.
C

      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNDXDI(PINT,IPAR,DXDI)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        calculates the transformation factor between external and
CC        internal parameter values.     this factor is one for
CC        parameters which are not limited.     called from MNEMAT.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      I = NEXOFI(IPAR)
      DXDI = 1.0
      IF (NVARL(I) .GT. 1)
     +      DXDI = 0.5 *ABS((BLIM(I)-ALIM(I)) * COS(PINT))
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEIG(A,NDIMA,N,MITS,WORK,PRECIS,IFAULT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0,  ONE=1.0,   TWO=2.0)
      PARAMETER (TOL=1.0E-35)
      DIMENSION A(NDIMA,*),WORK(*)
C          PRECIS is the machine precision EPSMAC
      IFAULT = 1
C
      I = N
      DO 70 I1 = 2,N
      L = I-2
      F = A(I,I-1)
      GL = ZERO
C
      IF(L .LT. 1) GO TO 25
C
      DO 20 K = 1,L
   20 GL = GL+A(I,K)**2
   25 H = GL + F**2
C
      IF(GL .GT. TOL) GO TO 30
C
      WORK(I) = ZERO
      WORK(N+I) = F
      GO TO 65
   30 L = L+1
C
      GL = SQRT(H)
C
      IF(F .GE. ZERO) GL = -GL
C
      WORK(N+I) = GL
      H = H-F*GL
      A(I,I-1) = F-GL
      F = ZERO
      DO 50 J = 1,L
      A(J,I) = A(I,J)/H
      GL = ZERO
      DO 40 K = 1,J
   40 GL = GL+A(J,K)*A(I,K)
C
      IF(J .GE. L) GO TO 47
C
      J1 = J+1
      DO 45 K = J1,L
   45 GL = GL+A(K,J)*A(I,K)
   47 WORK(N+J) = GL/H
      F = F+GL*A(J,I)
   50 CONTINUE
      HH = F/(H+H)
      DO 60 J = 1,L
      F = A(I,J)
      GL = WORK(N+J)-HH*F
      WORK(N+J) = GL
      DO 60 K = 1,J
      A(J,K) = A(J,K)-F*WORK(N+K)-GL*A(I,K)
   60 CONTINUE
      WORK(I) = H
   65 I = I-1
   70 CONTINUE
      WORK(1) = ZERO
      WORK(N+1) = ZERO
      DO 110 I = 1,N
      L = I-1
C
      IF(WORK(I) .EQ. ZERO .OR. L .EQ. 0) GO TO 100
C
      DO 90 J = 1,L
      GL = ZERO
      DO 80 K = 1,L
   80 GL = GL+A(I,K)*A(K,J)
      DO 90 K = 1,L
      A(K,J) = A(K,J)-GL*A(K,I)
   90 CONTINUE
  100 WORK(I) = A(I,I)
      A(I,I) = ONE
C
      IF(L .EQ. 0) GO TO 110
C
      DO 105 J = 1,L
      A(I,J) = ZERO
      A(J,I) = ZERO
  105 CONTINUE
  110 CONTINUE
C
C
      N1 = N-1
      DO 130 I = 2,N
      I0 = N+I-1
  130 WORK(I0) = WORK(I0+1)
      WORK(N+N) = ZERO
      B = ZERO
      F = ZERO
      DO 210 L = 1,N
      J = 0
      H = PRECIS*(ABS(WORK(L))+ABS(WORK(N+L)))
C
      IF(B .LT. H) B = H
C
      DO 140 M1 = L,N
      M = M1
C
      IF(ABS(WORK(N+M)) .LE. B) GO TO 150
C
  140 CONTINUE
C
  150 IF(M .EQ. L) GO TO 205
C
  160 IF(J .EQ. MITS) RETURN
C
      J = J+1
      PT = (WORK(L+1)-WORK(L))/(TWO*WORK(N+L))
      R = SQRT(PT*PT+ONE)
      PR = PT+R
C
      IF(PT .LT. ZERO) PR=PT-R
C
      H = WORK(L)-WORK(N+L)/PR
      DO 170 I=L,N
  170 WORK(I) = WORK(I)-H
      F = F+H
      PT = WORK(M)
      C = ONE
      S = ZERO
      M1 = M-1
      I = M
      DO 200 I1 = L,M1
      J = I
      I = I-1
      GL = C*WORK(N+I)
      H = C*PT
C
      IF(ABS(PT) .GE. ABS(WORK(N+I))) GO TO 180
C
      C = PT/WORK(N+I)
      R = SQRT(C*C+ONE)
      WORK(N+J) = S*WORK(N+I)*R
      S = ONE/R
      C = C/R
      GO TO 190
  180 C = WORK(N+I)/PT
      R = SQRT(C*C+ONE)
      WORK(N+J) = S*PT*R
      S = C/R
      C = ONE/R
  190 PT = C*WORK(I)-S*GL
      WORK(J) = H+S*(C*GL+S*WORK(I))
      DO 200 K = 1,N
      H = A(K,J)
      A(K,J) = S*A(K,I)+C*H
      A(K,I) = C*A(K,I)-S*H
  200 CONTINUE
      WORK(N+L) = S*PT
      WORK(L) = C*PT
C
      IF(ABS(WORK(N+L)) .GT. B) GO TO 160
C
  205 WORK(L) = WORK(L)+F
  210 CONTINUE
      DO 240 I=1,N1
      K = I
      PT = WORK(I)
      I1 = I+1
      DO 220 J = I1,N
C
      IF(WORK(J) .GE. PT) GO TO 220
C
      K = J
      PT = WORK(J)
  220 CONTINUE
C
      IF(K .EQ. I) GO TO 240
C
      WORK(K) = WORK(I)
      WORK(I) = PT
      DO 230 J=1,N
      PT = A(J,I)
      A(J,I) = A(J,K)
      A(J,K) = PT
  230 CONTINUE
  240 CONTINUE
      IFAULT = 0
C
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEMAT(EMAT,NDIM)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EMAT(NDIM,NDIM)
CC        Calculates the external error matrix from the internal
CC        to be called by user, who must dimension EMAT at (NDIM,NDIM)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      IF (ISW(2) .LT. 1)  RETURN
      IF (ISW(5) .GE. 2)  WRITE (ISYSWR,'(/A,I4,A,I3,A,G10.3)')
     +    ' EXTERNAL ERROR MATRIX.    NDIM=',NDIM,'    NPAR=',NPAR,
     +    '    ERR DEF=',UP
C                    size of matrix to be printed
      NPARD = NPAR
      IF (NDIM .LT. NPAR)  THEN
        NPARD = NDIM
        IF (ISW(5) .GE. 0) WRITE (ISYSWR,'(A,A)') ' USER-DIMENSIONED ',
     +      ' ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.'
      ENDIF
C                 NPERLN is the number of elements that fit on one line
      NPERLN = (NPAGWD-5)/10
      NPERLN = MIN(NPERLN,13)
      IF (ISW(5).GE. 1 .AND. NPARD.GT.NPERLN)  WRITE (ISYSWR,'(A)')
     +     ' ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.'
C                 I counts the rows of the matrix
      DO 110 I= 1, NPARD
         CALL MNDXDI(X(I),I,DXDI)
         KGA = I*(I-1)/2
         DO 100 J= 1, I
            CALL MNDXDI(X(J),J,DXDJ)
            KGB = KGA + J
            EMAT(I,J) = DXDI * VHMAT(KGB) * DXDJ * UP
            EMAT(J,I) = EMAT(I,J)
  100    CONTINUE
  110 CONTINUE
C                    IZ is number of columns to be printed in row I
      IF (ISW(5) .GE. 2)  THEN
      DO 160 I= 1, NPARD
         IZ = NPARD
         IF (NPARD .GE. NPERLN)  IZ = I
         DO 150 K= 1, IZ, NPERLN
           K2 = K + NPERLN - 1
           IF (K2 .GT. IZ)  K2=IZ
           WRITE (ISYSWR,'(1X,13E10.3)')  (EMAT(I,KK),KK=K,K2)
  150    CONTINUE
  160 CONTINUE
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNERRS(NUMBER,EPLUS,EMINUS,EPARAB,GCC)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC    Called by user, utility routine to get MINOS errors
CC    If NUMBER is positive, then it is external parameter number,
CC                  if negative, it is -internal number.
CC    values returned by MNERRS:
CC       EPLUS, EMINUS are MINOS errors of parameter NUMBER,
CC       EPARAB is 'parabolic' error (from error matrix).
CC                 (Errors not calculated are set = 0.)
CC       GCC is global correlation coefficient from error matrix
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      IEX = NUMBER
      IF (NUMBER .LT. 0)  THEN
         IIN = -NUMBER
         IF (IIN .GT. NPAR)  GO TO 900
         IEX = NEXOFI(IIN)
      ENDIF
      IF (IEX .GT. NU .OR. IEX .LE. 0)  GO TO 900
      IIN = NIOFEX(IEX)
      IF (IIN .LE. 0)  GO TO 900
C             IEX is external number, IIN is internal number
      EPLUS = ERP(IIN)
        IF (EPLUS.EQ.UNDEFI)  EPLUS=0.
      EMINUS= ERN(IIN)
        IF (EMINUS.EQ.UNDEFI) EMINUS=0.
      CALL MNDXDI(X(IIN),IIN,DXDI)
      NDIAG = IIN*(IIN+1)/2
      EPARAB = ABS(DXDI*SQRT(ABS(UP*VHMAT(NDIAG))))
C              global correlation coefficient
      GCC = 0.
      IF (ISW(2) .LT. 2)  GO TO 990
      GCC = GLOBCC(IIN)
      GO TO 990
C                  ERROR.  parameter number not valid
  900 EPLUS = 0.
      EMINUS = 0.
      EPARAB = 0.
      GCC = 0.
  990 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEVAL(FCN,ANEXT,FNEXT,IEREV,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC      Evaluates the function being analyzed by MNCROS, which is
CC      generally the minimum of FCN with respect to all remaining
CC      variable parameters.  Common block /MN7XCR/ contains the
CC      data necessary to know the values of U(KE1CR) and U(KE2CR)
CC      to be used, namely     U(KE1CR) = XMIDCR + ANEXT*XDIRCR
CC      and (if KE2CR .NE. 0)  U(KE2CR) = YMIDCR + ANEXT*YDIRCR
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
CC
      EXTERNAL FCN,FUTIL
                          U(KE1CR) = XMIDCR + ANEXT*XDIRCR
      IF ( KE2CR .NE. 0)  U(KE2CR) = YMIDCR + ANEXT*YDIRCR
      CALL MNINEX(X)
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FNEXT,U,4,FUTIL)
      NFCN = NFCN + 1
      IEREV = 0
      IF (NPAR .GT. 0)  THEN
         ITAUR = 1
         AMIN = FNEXT
         ISW(1) = 0
         CALL MNMIGR(FCN,FUTIL)
         ITAUR = 0
         FNEXT = AMIN
         IF (ISW(1) .GE. 1)  IEREV = 1
         IF (ISW(4) .LT. 1)  IEREV = 2
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:45  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEXCM(FCN,COMAND,PLIST,LLIST,IERFLG,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Interprets a command and takes appropriate action,
CC        either directly by skipping to the corresponding code in
CC        MNEXCM, or by setting up a call to a subroutine
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      CHARACTER*(*) COMAND
C   Cannot say DIMENSION PLIST(LLIST) since LLIST can be =0.
      DIMENSION PLIST(*)
      PARAMETER (MXPT=101)
      DIMENSION XPTU(MXPT), YPTU(MXPT)
C  alphabetical order of command names!
      CHARACTER*10 CNAME(40), CNEWAY, CHWHY*18, C26*30, CVBLNK*2
      LOGICAL LTOFIX, LFIXED, LFREED
C
      CHARACTER COMD*4
      CHARACTER CLOWER*26, CUPPER*26
      DATA CLOWER/'abcdefghijklmnopqrstuvwxyz'/
      DATA CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
C  recognized MINUIT commands:
      DATA CNAME( 1) / 'MINImize  ' /
      DATA CNAME( 2) / 'SEEk      ' /
      DATA CNAME( 3) / 'SIMplex   ' /
      DATA CNAME( 4) / 'MIGrad    ' /
      DATA CNAME( 5) / 'MINOs     ' /
      DATA CNAME( 6) / 'SET xxx   ' /
      DATA CNAME( 7) / 'SHOw xxx  ' /
      DATA CNAME( 8) / 'TOP of pag' /
      DATA CNAME( 9) / 'FIX       ' /
      DATA CNAME(10) / 'REStore   ' /
      DATA CNAME(11) / 'RELease   ' /
      DATA CNAME(12) / 'SCAn      ' /
      DATA CNAME(13) / 'CONtour   ' /
      DATA CNAME(14) / 'HESse     ' /
      DATA CNAME(15) / 'SAVe      ' /
      DATA CNAME(16) / 'IMProve   ' /
      DATA CNAME(17) / 'CALl fcn  ' /
      DATA CNAME(18) / 'STAndard  ' /
      DATA CNAME(19) / 'END       ' /
      DATA CNAME(20) / 'EXIt      ' /
      DATA CNAME(21) / 'RETurn    ' /
      DATA CNAME(22) / 'CLEar     ' /
      DATA CNAME(23) / 'HELP      ' /
      DATA CNAME(24) / 'MNContour ' /
      DATA CNAME(25) / 'STOp      ' /
      DATA CNAME(26) / 'JUMp      ' /
      DATA CNAME(27) / '          ' /
      DATA CNAME(28) / '          ' /
      DATA CNAME(29) / '          ' /
      DATA CNAME(30) / '          ' /
      DATA CNAME(31) / '          ' /
      DATA CNAME(32) / '          ' /
      DATA CNAME(33) / '          ' /
C  obsolete commands:
      DATA CNAME(34) / 'COVARIANCE' /
      DATA CNAME(35) / 'PRINTOUT  ' /
      DATA CNAME(36) / 'GRADIENT  ' /
      DATA CNAME(37) / 'MATOUT    ' /
      DATA CNAME(38) / 'ERROR DEF ' /
      DATA CNAME(39) / 'LIMITS    ' /
      DATA CNAME(40) / 'PUNCH     ' /
      DATA NNTOT/40/
C      IERFLG is now (94.5) defined the same as ICONDN in MNCOMD
CC            = 0: command executed normally
CC              1: command is blank, ignored
CC              2: command line unreadable, ignored
CC              3: unknown command, ignored
CC              4: abnormal termination (e.g., MIGRAD not converged)
CC              9: reserved
CC             10: END command
CC             11: EXIT or STOP command
CC             12: RETURN command
      LK = LEN(COMAND)
      IF (LK .GT. MAXCWD) LK=MAXCWD
      CWORD = COMAND(1:LK)
C              get upper case
      DO 16 ICOL= 1, LK
        DO 15 LET= 1, 26
        IF (CWORD(ICOL:ICOL) .EQ. CLOWER(LET:LET))
     +      CWORD(ICOL:ICOL) = CUPPER(LET:LET)
   15   CONTINUE
   16 CONTINUE
C           Copy the first MAXP arguments into COMMON (WORD7), making
C           sure that WORD7(1)=0. if LLIST=0
      DO 20 IW= 1, MAXP
      WORD7(IW) = ZERO
      IF (IW .LE. LLIST) WORD7(IW) = PLIST(IW)
   20 CONTINUE
      ICOMND = ICOMND + 1
      NFCNLC = NFCN
      IF (CWORD(1:7).NE.'SET PRI' .OR. WORD7(1).GE.0.)  THEN
        IF (ISW(5) .GE. 0) THEN
         LNOW = LLIST
         IF (LNOW .GT. 4)  LNOW=4
         WRITE (ISYSWR,25) ICOMND,CWORD(1:LK),(PLIST(I),I=1,LNOW)
   25    FORMAT (1H ,10(1H*)/' **',I5,' **',A,4G12.4)
         INONDE = 0
         IF (LLIST .GT. LNOW) THEN
           KLL = LLIST
           IF (LLIST .GT. MAXP) THEN
              INONDE = 1
              KLL = MAXP
           ENDIF
           WRITE (CVBLNK,'(I2)') LK
           C26 = '(11H **********,'//CVBLNK//'X,4G12.4)'
           WRITE (ISYSWR,C26) (PLIST(I),I=LNOW+1,KLL)
         ENDIF
         WRITE (ISYSWR, '(1H ,10(1H*))' )
         IF (INONDE .GT. 0)  WRITE (ISYSWR, '(1H ,10(1H*),A,I3,A)')
     +        '  ERROR: ABOVE CALL TO MNEXCM TRIED TO PASS MORE THAN ',
     +        MAXP,' PARAMETERS.'
        ENDIF
      ENDIF
      NFCNMX = WORD7(1)
      IF (NFCNMX .LE. 0)  NFCNMX = 200 + 100*NPAR + 5*NPAR**2
      EPSI = WORD7(2)
      IF (EPSI .LE. ZERO)  EPSI = 0.1 * UP
      LNEWMN = .FALSE.
      LPHEAD = .TRUE.
      ISW(1) = 0
      IERFLG = 0
C                look for command in list CNAME . . . . . . . . . .
      DO 80 I= 1, NNTOT
      IF (CWORD(1:3) .EQ. CNAME(I)(1:3))  GO TO 90
   80 CONTINUE
      WRITE (ISYSWR,'(11X,''UNKNOWN COMMAND IGNORED:'',A)') COMAND
      IERFLG = 3
      GO TO 5000
C                normal case: recognized MINUIT command . . . . . . .
   90 CONTINUE
      IF (CWORD(1:4) .EQ. 'MINO') I = 5
      IF (I.NE.6 .AND. I.NE.7 .AND. I.NE.8 .AND. I.NE.23)  THEN
         CFROM = CNAME(I)
         NFCNFR = NFCN
      ENDIF
C              1    2    3    4    5    6    7    8    9   10
      GO TO ( 400, 200, 300, 400, 500, 700, 700, 800, 900,1000,
     1       1100,1200,1300,1400,1500,1600,1700,1800,1900,1900,
     2       1900,2200,2300,2400,1900,2600,3300,3300,3300,3300,
     3       3300,3300,3300,3400,3500,3600,3700,3800,3900,4000) , I
C                                        . . . . . . . . . . seek
  200 CALL MNSEEK(FCN,FUTIL)
      GO TO 5000
C                                        . . . . . . . . . . simplex
  300 CALL MNSIMP(FCN,FUTIL)
      IF (ISW(4) .LT. 1)  IERFLG = 4
      GO TO 5000
C                                        . . . . . . migrad, minimize
  400 CONTINUE
      NF = NFCN
      APSI = EPSI
      CALL MNMIGR(FCN,FUTIL)
      CALL MNWERR
      IF (ISW(4) .GE. 1)         GO TO 5000
        IERFLG = 4
      IF (ISW(1) .EQ. 1)         GO TO 5000
      IF (CWORD(1:3) .EQ. 'MIG') GO TO 5000
      NFCNMX = NFCNMX + NF - NFCN
      NF = NFCN
      CALL MNSIMP(FCN,FUTIL)
      IF (ISW(1) .EQ. 1)  GO TO 5000
      NFCNMX = NFCNMX + NF - NFCN
      CALL MNMIGR(FCN,FUTIL)
         IF (ISW(4) .GE. 1)  IERFLG = 0
      CALL MNWERR
      GO TO 5000
C                                        . . . . . . . . . . minos
  500 CONTINUE
      NSUPER = NFCN + 2*(NPAR+1)*NFCNMX
C          possible loop over new minima
      EPSI = 0.1 * UP
  510 CONTINUE
      CALL MNCUVE(FCN,FUTIL)
      CALL MNMNOS(FCN,FUTIL)
      IF (.NOT. LNEWMN)  GO TO 5000
      CALL MNRSET(0)
      CALL MNMIGR(FCN,FUTIL)
      CALL MNWERR
      IF (NFCN .LT. NSUPER)  GO TO 510
      WRITE (ISYSWR,'(/'' TOO MANY FUNCTION CALLS. MINOS GIVES UP''/)')
      IERFLG = 4
      GO TO 5000
C                                        . . . . . . . . . .set, show
  700 CALL MNSET(FCN,FUTIL)
      GO TO 5000
C                                        . . . . . . . . . . top of page
  800 CONTINUE
      WRITE (ISYSWR,'(1H1)')
      GO TO 5000
C                                        . . . . . . . . . . fix
  900 LTOFIX = .TRUE.
C                                        . . (also release) ....
  901 CONTINUE
      LFREED = .FALSE.
      LFIXED = .FALSE.
      IF (LLIST .EQ. 0)  THEN
         WRITE (ISYSWR,'(A,A)') CWORD,':  NO PARAMETERS REQUESTED '
         GO TO 5000
      ENDIF
      DO 950 ILIST= 1, LLIST
      IEXT = PLIST(ILIST)
      CHWHY = ' IS UNDEFINED.'
      IF (IEXT .LE. 0)         GO TO 930
      IF (IEXT .GT. NU)        GO TO 930
      IF (NVARL(IEXT) .LT. 0)  GO TO 930
      CHWHY = ' IS CONSTANT.  '
      IF (NVARL(IEXT) .EQ. 0)  GO TO 930
      IINT = NIOFEX(IEXT)
      IF (LTOFIX) THEN
         CHWHY = ' ALREADY FIXED.'
         IF (IINT .EQ. 0)      GO TO 930
         CALL MNFIXP(IINT,IERR)
         IF (IERR .EQ. 0) THEN
            LFIXED = .TRUE.
         ELSE
            IERFLG = 4
         ENDIF
      ELSE
         CHWHY = ' ALREADY VARIABLE.'
         IF (IINT .GT. 0)      GO TO 930
         KRL = -IABS(IEXT)
         CALL MNFREE(KRL)
         LFREED = .TRUE.
      ENDIF
      GO TO 950
  930 WRITE (ISYSWR,'(A,I4,A,A)') ' PARAMETER',IEXT,CHWHY,' IGNORED.'
  950 CONTINUE
      IF (LFREED .OR. LFIXED)  CALL MNRSET(0)
      IF (LFREED)  THEN
          ISW(2) = 0
          DCOVAR = 1.
          EDM = BIGEDM
          ISW(4) = 0
      ENDIF
      CALL MNWERR
      IF (ISW(5) .GT. 1)  CALL MNPRIN(5,AMIN)
      GO TO 5000
C                                        . . . . . . . . . . restore
 1000 IT = WORD7(1)
      IF (IT.GT.1 .OR. IT.LT.0)  GO TO 1005
      LFREED = (NPFIX .GT. 0)
      CALL MNFREE(IT)
      IF (LFREED) THEN
         CALL MNRSET(0)
         ISW(2) = 0
         DCOVAR = 1.
         EDM = BIGEDM
      ENDIF
      GO TO 5000
 1005 WRITE (ISYSWR,'(A,I4)') ' IGNORED.  UNKNOWN ARGUMENT:',IT
      IERFLG = 3
      GO TO 5000
C                                        . . . . . . . . . . release
 1100 LTOFIX = .FALSE.
      GO TO 901
C                                       . . . . . . . . . . scan . . .
 1200 CONTINUE
      IEXT = WORD7(1)
      IF (IEXT .LE. 0)  GO TO 1210
      IT2 = 0
      IF (IEXT .LE. NU)  IT2 = NIOFEX(IEXT)
      IF (IT2 .LE. 0)  GO TO 1250
 1210 CALL MNSCAN(FCN,FUTIL)
      GO TO 5000
 1250 WRITE (ISYSWR,'(A,I4,A)') ' PARAMETER',IEXT,' NOT VARIABLE.'
      IERFLG = 3
      GO TO 5000
C                                        . . . . . . . . . . contour
 1300 CONTINUE
      KE1 = WORD7(1)
      KE2 = WORD7(2)
      IF (KE1 .EQ. 0)  THEN
         IF (NPAR .EQ. 2)  THEN
            KE1 = NEXOFI(1)
            KE2 = NEXOFI(2)
         ELSE
            WRITE (ISYSWR,'(A,A)') CWORD,':  NO PARAMETERS REQUESTED '
            IERFLG = 3
            GO TO 5000
         ENDIF
      ENDIF
      NFCNMX = 1000
      CALL MNCNTR(FCN,KE1,KE2,IERRF,FUTIL)
      IF (IERRF .GT. 0)  IERFLG = 3
      GO TO 5000
C                                        . . . . . . . . . . hesse
 1400 CONTINUE
      CALL MNHESS(FCN,FUTIL)
      CALL MNWERR
      IF (ISW(5) .GE. 0)  CALL MNPRIN(2, AMIN)
      IF (ISW(5) .GE. 1)  CALL MNMATU(1)
      GO TO 5000
C                                        . . . . . . . . . . save
 1500 CONTINUE
      CALL MNSAVE
      GO TO 5000
C                                        . . . . . . . . . . improve
 1600 CONTINUE
      CALL MNCUVE(FCN,FUTIL)
      CALL MNIMPR(FCN,FUTIL)
      IF (LNEWMN)  GO TO 400
      IERFLG = 4
      GO TO 5000
C                                        . . . . . . . . . . call fcn
 1700 IFLAG = WORD7(1)
      NPARX = NPAR
      F = UNDEFI
      CALL FCN(NPARX,GIN,F,U,IFLAG,FUTIL)
      NFCN = NFCN + 1
      NOWPRT = 0
      IF (F .NE. UNDEFI)  THEN
         IF (AMIN .EQ. UNDEFI)  THEN
             AMIN = F
             NOWPRT = 1
         ELSE IF (F .LT. AMIN)  THEN
             AMIN = F
             NOWPRT = 1
         ENDIF
         IF (ISW(5).GE.0 .AND. IFLAG.LE.5 .AND. NOWPRT.EQ.1)
     +          CALL MNPRIN(5,AMIN)
         IF (IFLAG .EQ. 3)  FVAL3=F
      ENDIF
      IF (IFLAG .GT. 5)  CALL MNRSET(1)
      GO TO 5000
C                                        . . . . . . . . . . standard
 1800 CALL STAND
      GO TO 5000
C                                       . . . return, stop, end, exit
 1900 IT = WORD7(1)
      IF (FVAL3 .NE. AMIN .AND. IT .EQ. 0)  THEN
        IFLAG = 3
        IF (ISW(5) .GE. 0)
     +WRITE (ISYSWR,'(/A/)') ' CALL TO USER FUNCTION WITH IFLAG = 3'
        NPARX = NPAR
        CALL FCN(NPARX,GIN,F,U,IFLAG,FUTIL)
        NFCN = NFCN + 1
        FVAL3 = F
      ENDIF
      IERFLG = 11
      IF (CWORD(1:3) .EQ. 'END')  IERFLG = 10
      IF (CWORD(1:3) .EQ. 'RET')  IERFLG = 12
      GO TO 5000
C                                        . . . . . . . . . . clear
 2200 CONTINUE
      CALL MNCLER
      IF (ISW(5) .GE. 1)  WRITE (ISYSWR,'(A)')
     + ' MINUIT MEMORY CLEARED. NO PARAMETERS NOW DEFINED.'
      GO TO 5000
C                                        . . . . . . . . . . help
 2300 CONTINUE
CCCC      IF (INDEX(CWORD,'SHO') .GT. 0)  GO TO 700
CCCC      IF (INDEX(CWORD,'SET') .GT. 0)  GO TO 700
      KCOL = 0
      DO 2310 ICOL= 5,LK
        IF (CWORD(ICOL:ICOL) .EQ. ' ') GO TO 2310
        KCOL = ICOL
        GO TO 2320
 2310 CONTINUE
 2320 CONTINUE
      IF (KCOL .EQ. 0)  THEN
         COMD = '*   '
      ELSE
         COMD = CWORD(KCOL:LK)
      ENDIF
      CALL MNHELP(COMD,ISYSWR)
      GO TO 5000
C                                       . . . . . . . . . . MNContour
 2400 CONTINUE
      EPSI = 0.05 * UP
      KE1 = WORD7(1)
      KE2 = WORD7(2)
      IF (KE1.EQ.0 .AND. NPAR.EQ.2) THEN
         KE1 = NEXOFI(1)
         KE2 = NEXOFI(2)
         ENDIF
      NPTU = WORD7(3)
      IF (NPTU .LE. 0)  NPTU=20
      IF (NPTU .GT. MXPT)  NPTU = MXPT
      NFCNMX =  100*(NPTU+5)*(NPAR+1)
      CALL MNCONT(FCN,KE1,KE2,NPTU,XPTU,YPTU,IERRF,FUTIL)
      IF (IERRF .LT. NPTU) IERFLG = 4
      IF (IERRF .EQ. -1)   IERFLG = 3
      GO TO 5000
C                                      . . . . . . . . . . jump
 2600 CONTINUE
      STEP = WORD7(1)
      IF (STEP .LE. ZERO)  STEP = 2.
      RNO = 0.
      IZERO = 0
      DO 2620 I= 1, NPAR
        CALL MNRN15(RNO,IZERO)
        RNO = 2.0*RNO - 1.0
 2620   X(I) = X(I) + RNO*STEP*WERR(I)
      CALL MNINEX(X)
      CALL MNAMIN(FCN,FUTIL)
      CALL MNRSET(0)
      GO TO 5000
C                                      . . . . . . . . . . blank line
 3300 CONTINUE
      WRITE (ISYSWR,'(10X,A)') ' BLANK COMMAND IGNORED.'
      IERFLG = 1
      GO TO 5000
C  . . . . . . . . obsolete commands     . . . . . . . . . . . . . .
C                                      . . . . . . . . . . covariance
 3400 CONTINUE
      WRITE (ISYSWR, '(A)') ' THE "COVARIANCE" COMMAND IS OSBSOLETE.',
     + ' THE COVARIANCE MATRIX IS NOW SAVED IN A DIFFERENT FORMAT',
     + ' WITH THE "SAVE" COMMAND AND READ IN WITH:"SET COVARIANCE"'
      IERFLG = 3
      GO TO 5000
C                                        . . . . . . . . . . printout
 3500 CONTINUE
      CNEWAY = 'SET PRInt '
      GO TO 3100
C                                        . . . . . . . . . . gradient
 3600 CONTINUE
      CNEWAY = 'SET GRAd  '
      GO TO 3100
C                                        . . . . . . . . . . matout
 3700 CONTINUE
      CNEWAY = 'SHOW COVar'
      GO TO 3100
C                                        . . . . . . . . . error def
 3800 CONTINUE
      CNEWAY = 'SET ERRdef'
      GO TO 3100
C                                        . . . . . . . . . . limits
 3900 CONTINUE
      CNEWAY = 'SET LIMits'
      GO TO 3100
C                                        . . . . . . . . . . punch
 4000 CONTINUE
      CNEWAY = 'SAVE      '
C                                ....... come from obsolete commands
 3100 WRITE (ISYSWR, 3101) CWORD,CNEWAY
 3101 FORMAT (' OBSOLETE COMMAND:',1X,A10,5X,'PLEASE USE:',1X,A10)
      CWORD = CNEWAY
      IF (CWORD .EQ. 'SAVE      ') GO TO 1500
      GO TO 700
C                                 . . . . . . . . . . . . . . . . . .
 5000 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEXIN(PINT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Transforms the external parameter values U to internal
CC        values in the dense array PINT. Subroutine MNPINT is used.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION PINT(*)
      LIMSET = .FALSE.
      DO 100  IINT= 1, NPAR
      IEXT = NEXOFI(IINT)
      CALL MNPINT(U(IEXT),IEXT,PINTI)
      PINT(IINT) = PINTI
  100 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNFIXP(IINT,IERR)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        removes parameter IINT from the internal (variable) parameter
CC        list, and arranges the rest of the list to fill the hole.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION YY(MNI)
C                           first see if it can be done
      IERR = 0
      IF (IINT.GT.NPAR .OR. IINT.LE.0)  THEN
         IERR = 1
         WRITE (ISYSWR,'(A,I4)')
     +       ' MINUIT ERROR.  ARGUMENT TO MNFIXP=',IINT
         GO TO 300
      ENDIF
      IEXT = NEXOFI(IINT)
      IF (NPFIX .GE. MNI) THEN
         IERR = 1
         WRITE (ISYSWR,'(A,I4,A,I4)') ' MINUIT CANNOT FIX PARAMETER',
     +   IEXT,' MAXIMUM NUMBER THAT CAN BE FIXED IS',MNI
         GO TO 300
      ENDIF
C                           reduce number of variable parameters by one
      NIOFEX(IEXT) = 0
      NOLD = NPAR
      NPAR = NPAR - 1
C                       save values in case parameter is later restored
      NPFIX = NPFIX + 1
      IPFIX(NPFIX) = IEXT
      LC = IINT
      XS(NPFIX) = X(LC)
      XTS(NPFIX) = XT(LC)
      DIRINS(NPFIX) = WERR(LC)
      GRDS(NPFIX) = GRD(LC)
      G2S(NPFIX) = G2(LC)
      GSTEPS(NPFIX) = GSTEP(LC)
C                        shift values for other parameters to fill hole
      DO 100  IK= IEXT+1, NU
         IF  (NIOFEX(IK) .GT. 0)  THEN
         LC = NIOFEX(IK) - 1
         NIOFEX(IK) = LC
         NEXOFI(LC) = IK
         X(LC)     = X(LC+1)
         XT(LC)    = XT(LC+1)
         DIRIN(LC) = DIRIN(LC+1)
         WERR(LC)  = WERR(LC+1)
         GRD(LC)   = GRD(LC+1)
         G2(LC)    = G2(LC+1)
         GSTEP(LC) = GSTEP(LC+1)
         ENDIF
  100 CONTINUE
      IF (ISW(2) .LE. 0)  GO TO 300
C                    remove one row and one column from variance matrix
      IF (NPAR .LE. 0)  GO TO 300
      DO 260 I= 1, NOLD
      M = MAX(I,IINT)
      N = MIN(I,IINT)
      NDEX = M*(M-1)/2 + N
  260 YY(I)=VHMAT(NDEX)
      YYOVER = 1.0/YY(IINT)
      KNEW = 0
      KOLD = 0
      DO 294 I= 1, NOLD
      DO 292 J= 1, I
      KOLD = KOLD + 1
      IF (J.EQ.IINT .OR. I.EQ.IINT)  GO TO 292
      KNEW = KNEW + 1
      VHMAT(KNEW) = VHMAT(KOLD) - YY(J)*YY(I)*YYOVER
  292 CONTINUE
  294 CONTINUE
  300 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNFREE(K)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Restores one or more fixed parameter(s) to variable status
CC        by inserting it into the internal parameter list at the
CC        appropriate place.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C--       K = 0 means restore all parameters
C--       K = 1 means restore the last parameter fixed
C--       K = -I means restore external parameter I (if possible)
C--       IQ = fix-location where internal parameters were stored
C--       IR = external number of parameter being restored
C--       IS = internal number of parameter being restored
      IF (K .GT. 1)  WRITE (ISYSWR,510)
      IF (NPFIX .LT. 1)  WRITE (ISYSWR,500)
      IF (K.EQ.1 .OR. K.EQ.0)  GO TO 40
C                   release parameter with specified external number
      KA = IABS(K)
      IF (NIOFEX(KA) .EQ. 0)  GO TO 15
      WRITE (ISYSWR,540)
  540 FORMAT (' IGNORED.  PARAMETER SPECIFIED IS ALREADY VARIABLE.')
      RETURN
   15 IF (NPFIX .LT. 1)  GO TO 21
      DO 20 IK= 1, NPFIX
      IF (IPFIX(IK) .EQ. KA)  GO TO 24
   20 CONTINUE
   21 WRITE (ISYSWR,530) KA
  530 FORMAT (' PARAMETER',I4,' NOT FIXED.  CANNOT BE RELEASED.')
      RETURN
   24 IF (IK .EQ. NPFIX)  GO TO 40
C                   move specified parameter to end of list
      IPSAV = KA
      XV = XS(IK)
      XTV = XTS(IK)
      DIRINV = DIRINS(IK)
      GRDV = GRDS(IK)
      G2V = G2S(IK)
      GSTEPV = GSTEPS(IK)
         DO 30 I= IK+1,NPFIX
         IPFIX(I-1) = IPFIX(I)
         XS(I-1) = XS(I)
         XTS(I-1) = XTS(I)
         DIRINS(I-1) = DIRINS(I)
         GRDS(I-1) = GRDS(I)
         G2S(I-1) = G2S(I)
         GSTEPS(I-1) = GSTEPS(I)
   30    CONTINUE
      IPFIX(NPFIX) = IPSAV
      XS(NPFIX) = XV
      XTS(NPFIX) = XTV
      DIRINS(NPFIX) = DIRINV
      GRDS(NPFIX) = GRDV
      G2S(NPFIX) = G2V
      GSTEPS(NPFIX) = GSTEPV
C                restore last parameter in fixed list  -- IPFIX(NPFIX)
   40 CONTINUE
      IF (NPFIX .LT. 1)  GO TO 300
      IR = IPFIX(NPFIX)
      IS = 0
      DO 100 IK= NU, IR, -1
        IF (NIOFEX(IK) .GT. 0) THEN
         LC = NIOFEX(IK) + 1
         IS = LC - 1
         NIOFEX(IK) = LC
         NEXOFI(LC) = IK
         X(LC)     = X(LC-1)
         XT(LC)    = XT(LC-1)
         DIRIN(LC) = DIRIN(LC-1)
         WERR(LC)  = WERR(LC-1)
         GRD(LC)   = GRD(LC-1)
         G2(LC)    = G2(LC-1)
         GSTEP(LC) = GSTEP(LC-1)
        ENDIF
  100 CONTINUE
      NPAR = NPAR + 1
      IF (IS .EQ. 0)   IS = NPAR
      NIOFEX(IR) = IS
      NEXOFI(IS) = IR
      IQ = NPFIX
      X(IS) = XS(IQ)
      XT(IS) = XTS(IQ)
      DIRIN(IS) = DIRINS(IQ)
      WERR(IS)  = DIRINS(IQ)
      GRD(IS) = GRDS(IQ)
      G2(IS) = G2S(IQ)
      GSTEP(IS) = GSTEPS(IQ)
      NPFIX = NPFIX - 1
      ISW(2) = 0
      DCOVAR = 1.
      IF (ISW(5)-ITAUR .GE. 1)  WRITE(ISYSWR,520) IR,CPNAM(IR)
      IF (K.EQ.0)  GO TO 40
  300 CONTINUE
C         if different from internal, external values are taken
      CALL MNEXIN(X)
  400 RETURN
  500 FORMAT (' CALL TO MNFREE IGNORED.  THERE ARE NO FIXED PA',
     + 'RAMETERS'/)
  510 FORMAT (' CALL TO MNFREE IGNORED.  ARGUMENT GREATER THAN ONE'/)
  520 FORMAT (20X, 9HPARAMETER,I4,2H, ,A10,' RESTORED TO VARIABLE.')
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNGRAD(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Called from MNSET
CC       Interprets the SET GRAD command, which informs MINUIT whether
CC       the first derivatives of FCN will be calculated by the user
CC       inside FCN.  It can check the user's derivative calculation
CC       by comparing it with a finite difference approximation.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      EXTERNAL FCN,FUTIL
      CHARACTER*4 CGOOD,CBAD,CNONE,CWD
      LOGICAL LNONE
      DIMENSION GF(MNI)
      PARAMETER (CGOOD='GOOD',CBAD=' BAD',CNONE='NONE')
C
      ISW(3) = 1
      NPARX = NPAR
      IF (WORD7(1) .GT. ZERO)  GO TO 2000
C                  get user-calculated first derivatives from FCN
      DO 30 I= 1, NU
   30 GIN(I) = UNDEFI
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FZERO,U,2,FUTIL)
      NFCN = NFCN + 1
      CALL MNDERI(FCN,FUTIL)
      DO 40 I= 1, NPAR
   40 GF(I) = GRD(I)
C                    get MINUIT-calculated first derivatives
      ISW(3) = 0
      ISTSAV = ISTRAT
      ISTRAT = 2
      CALL MNHES1(FCN,FUTIL)
      ISTRAT = ISTSAV
      WRITE (ISYSWR,51)
   51 FORMAT(/' CHECK OF GRADIENT CALCULATION IN FCN'/12X,'PARAMETER',
     + 6X,9HG(IN FCN) ,3X,9HG(MINUIT) ,2X,'DG(MINUIT)',3X,9HAGREEMENT)
      ISW(3) = 1
      LNONE = .FALSE.
      DO 100 LC = 1, NPAR
      I = NEXOFI(LC)
      CWD = CGOOD
      ERR = DGRD(LC)
      IF (ABS(GF(LC)-GRD(LC)) .GT. ERR)  CWD = CBAD
      IF (GIN(I) .EQ. UNDEFI)  THEN
          CWD = CNONE
          LNONE = .TRUE.
          GF(LC) = 0.
          ENDIF
      IF (CWD .NE. CGOOD)  ISW(3) = 0
      WRITE (ISYSWR,99) I,CPNAM(I),GF(LC),GRD(LC),ERR,CWD
   99 FORMAT (7X,I5,2X ,A10,3E12.4,4X ,A4)
  100 CONTINUE
      IF (LNONE) WRITE (ISYSWR,'(A)')
     +  '  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DERIVATIVE'
      IF (ISW(3) .EQ. 0)  WRITE (ISYSWR,1003)
 1003 FORMAT(/' MINUIT DOES NOT ACCEPT DERIVATIVE CALCULATIONS BY FCN'/
     + ' TO FORCE ACCEPTANCE, ENTER "SET GRAD    1"'/)
C
 2000 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1999/09/03 09:17:47  couet
* - \Cind{} removed in the help of minuit. This was a Tex directive which very
*   likely has been forgotten during a Tex to f77 translation. This didn't
*   compile on RH6.
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNHELP(COMD,LOUT)
*.
*.         HELP routine for MINUIT interactive commands.
*.
*.      COMD ='*   '  prints a global help for all commands
*.      COMD =Command_name: print detailed help for one command.
*.          Note that at least 3 characters must be given for the command name.
*.
*.     Author: Rene Brun
*             comments extracted from the MINUIT documentation file.
*.
      CHARACTER*(*) COMD
      CHARACTER*3 CMD3
*.
*-- command name ASSUMED to be in upper case
*__________________________________________________________________
*--
*--  Global HELP: Summary of all commands
*--  ====================================
*--
      IF(COMD(1:1) .EQ. '*')THEN
         WRITE(LOUT,10000)
         WRITE(LOUT,10001)
         GO TO 99
      ENDIF
10000 FORMAT('   ==>List of MINUIT Interactive commands:',/,
     +' CLEar     Reset all parameter names and values undefined',/,
     +' CONtour   Make contour map of the user function',/,
     +' EXIT      Exit from Interactive Minuit',/,
     +' FIX       Cause parameter(s) to remain constant',/,
     +' HESse     Calculate the Hessian or error matrix.',/,
     +' IMPROVE   Search for a new minimum around current minimum',/,
     +' MIGrad    Minimize by the method of Migrad',/,
     +' MINImize  MIGRAD + SIMPLEX method if Migrad fails',/,
     +' MINOs     Exact (non-linear) parameter error analysis')
10001 FORMAT(' MNContour Calculate one MINOS function contour',/,
     +' PARameter Define or redefine new parameters and values',/,
     +' RELease   Make previously FIXed parameters variable again',/,
     +' REStore   Release last parameter fixed',/,
     +' SAVe      Save current parameter values on a file',/,
     +' SCAn      Scan the user function by varying parameters',/,
     +' SEEk      Minimize by the method of Monte Carlo',/,
     +' SET       Set various MINUIT constants or conditions',/,
     +' SHOw      Show values of current constants or conditions',/,
     +' SIMplex   Minimize by the method of Simplex')
*
      CMD3=COMD(1:3)
*__________________________________________________________________
*--
*--  Command CLEAR
*--  =============
*.
      IF(CMD3.EQ.'CLE')THEN
         WRITE(LOUT,10100)
         GO TO 99
      ENDIF
10100 FORMAT(' ***>CLEAR',/,
     +' Resets all parameter names and values to undefined.',/,
     +' Must normally be followed by a PARameters command or ',/,
     +' equivalent, in order to define parameter values.')
*__________________________________________________________________
*--
*--  Command CONTOUR
*--  ===============
*.
      IF(CMD3.EQ.'CON')THEN
         WRITE(LOUT,10200)
         GO TO 99
      ENDIF
10200 FORMAT(' ***>CONTOUR <par1>  <par2>  [devs]  [ngrid]',/,
     +' Instructs Minuit to trace contour lines of the user function',/,
     +' with respect to the two parameters whose external numbers',/,
     +' are <par1> and <par2>.',/,
     +' Other variable parameters of the function, if any, will have',/,
     +' their values fixed at the current values during the contour',/,
     +' tracing. The optional parameter [devs] (default value 2.)',/,
     +' gives the number of standard deviations in each parameter',/,
     +' which should lie entirely within the plotting area.',/,
     +' Optional parameter [ngrid] (default value 25 unless page',/,
     +' size is too small) determines the resolution of the plot,',/,
     +' i.e. the number of rows and columns of the grid at which the',/,
     +' function will be evaluated. [See also MNContour.]')
*__________________________________________________________________
*--
*--  Command END
*--  ===========
*.
      IF(CMD3.EQ.'END')THEN
         WRITE(LOUT,10300)
         GO TO 99
      ENDIF
10300 FORMAT(' ***>END',/,
     +' Signals the end of a data block (i.e., the end of a fit),',/,
     +' and implies that execution should continue, because another',/,
     +' Data Block follows. A Data Block is a set of Minuit data',/,
     +' consisting of',/,
     +'     (1) A Title,',/,
     +'     (2) One or more Parameter Definitions,',/,
     +'     (3) A blank line, and',/,
     +'     (4) A set of Minuit Commands.',/,
     +' The END command is used when more than one Data Block is to',/,
     +' be used with the same FCN function. It first causes Minuit',/,
     +' to issue a CALL FCN with IFLAG=3, in order to allow FCN to',/,
     +' perform any calculations associated with the final fitted',/,
     +' parameter values, unless a CALL FCN 3 command has already',/,
     +' been executed at the current FCN value.')
*__________________________________________________________________
*.
*--
*--  Command EXIT
*--  ============
      IF(CMD3 .EQ.'EXI')THEN
         WRITE(LOUT,10400)
         GO TO 99
      ENDIF
10400 FORMAT(' ***>EXIT',/,
     +' Signals the end of execution.',/,
     +' The EXIT command first causes Minuit to issue a CALL FCN',/,
     +' with IFLAG=3, to allow FCN to perform any calculations',/,
     +' associated with the final fitted parameter values, unless a',/,
     +' CALL FCN 3 command has already been executed.')
*__________________________________________________________________
*--
*--  Command FIX
*--  ===========
*.
      IF(CMD3.EQ.'FIX')THEN
         WRITE(LOUT,10500)
         GO TO 99
      ENDIF
10500 FORMAT(' ***>FIX} <parno> [parno] ... [parno]',/,
     +' Causes parameter(s) <parno> to be removed from the list of',/,
     +' variable parameters, and their value(s) will remain constant',/,
     +' during subsequent minimizations, etc., until another command',/,
     +' changes their value(s) or status.')
*__________________________________________________________________
*--
*--  Command HESSE
*--  =============
*.
      IF(CMD3.EQ.'HES')THEN
         WRITE(LOUT,10600)
         GO TO 99
      ENDIF
10600 FORMAT(' ***>HESse  [maxcalls]',/,
     +' Calculate, by finite differences, the Hessian or error matrix.',
     +/,'  That is, it calculates the full matrix of second derivatives'
     +,/,' of the function with respect to the currently variable',/,
     +' parameters, and inverts it, printing out the resulting error',/,
     +' matrix. The optional argument [maxcalls] specifies the',/,
     +' (approximate) maximum number of function calls after which',/,
     +' the calculation will be stopped.')
*__________________________________________________________________
*--
*--  Command IMPROVE
*--  ===============
*.
      IF(CMD3.EQ.'IMP')THEN
         WRITE(LOUT,10700)
         GO TO 99
      ENDIF
10700 FORMAT(' ***>IMPROVE  [maxcalls]',/,
     +' If a previous minimization has converged, and the current',/,
     +' values of the parameters therefore correspond to a local',/,
     +' minimum of the function, this command requests a search for',/,
     +' additional distinct local minima.',/,
     +' The optional argument [maxcalls] specifies the (approximate)',/,
     +' maximum number of function calls after which the calculation',/,
     +' will be stopped.')
*__________________________________________________________________
*--
*--  Command MIGRAD
*--  ==============
*.
      IF(CMD3.EQ.'MIG')THEN
         WRITE(LOUT,10800)
         GO TO 99
      ENDIF
10800 FORMAT(' ***>MIGrad  [maxcalls]  [tolerance]',/,
     +' Causes minimization of the function by the method of Migrad,',/,
     +' the most efficient and complete single method, recommended',/,
     +' for general functions (see also MINImize).',/,
     +' The minimization produces as a by-product the error matrix',/,
     +' of the parameters, which is usually reliable unless warning',/,
     +' messages are produced.',/,
     +' The optional argument [maxcalls] specifies the (approximate)',/,
     +' maximum number of function calls after which the calculation',/,
     +' will be stopped even if it has not yet converged.',/,
     +' The optional argument [tolerance] specifies required tolerance',
     +/,' on the function value at the minimum.',/,
     +' The default tolerance is 0.1, and the minimization will stop',/,
     +' when the estimated vertical distance to the minimum (EDM) is',/,
     +' less than 0.001*[tolerance]*UP (see [SET ERRordef]).')
*__________________________________________________________________
*--
*--  Command MINIMIZE
*--  ================
*.
      IF(COMD(1:4).EQ.'MINI')THEN
         WRITE(LOUT,10900)
         GO TO 99
      ENDIF
10900 FORMAT(' ***>MINImize  [maxcalls] [tolerance]',/,
     +' Causes minimization of the function by the method of Migrad,',/,
     +' as does the MIGrad command, but switches to the SIMplex method',
     +/,' if Migrad fails to converge. Arguments are as for MIGrad.',/,
     +' Note that command requires four characters to be unambiguous.')
*__________________________________________________________________
*--
*--  Command MINOS
*--  =============
*.
      IF(COMD(1:4).EQ.'MINO')THEN
         WRITE(LOUT,11000)
         GO TO 99
      ENDIF
11000 FORMAT(' ***>MINOs  [maxcalls]  [parno] [parno] ...',/,
     +' Causes a Minos error analysis to be performed on the parameters'
     +,/,' whose numbers [parno] are specified. If none are specified,',
     +/,' Minos errors are calculated for all variable parameters.',/,
     +' Minos errors may be expensive to calculate, but are very',/,
     +' reliable since they take account of non-linearities in the',/,
     +' problem as well as parameter correlations, and are in general',/
     +' asymmetric.',/,
     +' The optional argument [maxcalls] specifies the (approximate)',/,
     +' maximum number of function calls per parameter requested,',/,
     +' after which the calculation will stop for that parameter.')
*__________________________________________________________________
*--
*--  Command MNCONTOUR
*--  =================
*.
      IF(CMD3.EQ.'MNC')THEN
         WRITE(LOUT,11100)
         GO TO 99
      ENDIF
11100 FORMAT(' ***>MNContour  <par1> <par2> [npts]',/,
     +' Calculates one function contour of FCN with respect to',/,
     +' parameters par1 and par2, with FCN minimized always with',/,
     +' respect to all other NPAR-2 variable parameters (if any).',/,
     +' Minuit will try to find npts points on the contour (default 20)'
     +,/,' If only two parameters are variable at the time, it is not',
     +/,' necessary to specify their numbers. To calculate more than',/,
     +' one contour, it is necessary to SET ERRordef to the appropriate'
     +,/,' value and issue the MNContour command for each contour.')
*__________________________________________________________________
*--
*--  Command PARAMETER
*--  =================
*.
      IF(CMD3.EQ.'PAR')THEN
         WRITE(LOUT,11150)
         GO TO 99
      ENDIF
11150 FORMAT(' ***>PARameters',/,
     +' followed by one or more parameter definitions.',/,
     +' Parameter definitions are of the form:',/,
     +'   <number>  ''name''  <value>  <step>  [lolim] [uplim] ',/,
     +' for example:',/,
     +'  3  ''K width''  1.2   0.1' ,/,
     +' the last definition is followed by a blank line or a zero.')
*__________________________________________________________________
*--
*--  Command RELEASE
*--  ===============
*.
      IF(CMD3.EQ.'REL')THEN
         WRITE(LOUT,11200)
         GO TO 99
      ENDIF
11200 FORMAT(' ***>RELease  <parno> [parno] ... [parno]',/,
     +' If <parno> is the number of a previously variable parameter',/,
     +' which has been fixed by a command: FIX <parno>, then that',/,
     +' parameter will return to variable status.  Otherwise a warning'
     +,/,' message is printed and the command is ignored.',/,
     +' Note that this command operates only on parameters which were',/
     +' at one time variable and have been FIXed. It cannot make',/,
     +' constant parameters variable; that must be done by redefining',/
     +' the parameter with a PARameters command.')
*__________________________________________________________________
*--
*--  Command RESTORE
*--  ===============
*.
      IF(CMD3.EQ.'RES')THEN
         WRITE(LOUT,11300)
         GO TO 99
      ENDIF
11300 FORMAT(' ***>REStore  [code]',/,
     +' If no [code] is specified, this command restores all previously'
     +,/,' FIXed parameters to variable status. If [code]=1, then only',
     +/,' the last parameter FIXed is restored to variable status.',/,
     +' If code is neither zero nor one, the command is ignored.')
*__________________________________________________________________
*--
*--  Command RETURN
*--  ==============
*.
      IF(CMD3.EQ.'RET')THEN
         WRITE(LOUT,11400)
         GO TO 99
      ENDIF
11400 FORMAT(' ***>RETURN',/,
     +' Signals the end of a data block, and instructs Minuit to return'
     +,/,' to the program which called it. The RETurn command first',/,
     +' causes Minuit to CALL FCN with IFLAG=3, in order to allow FCN',/
     +,' to perform any calculations associated with the final fitted',/
     +,' parameter values, unless a CALL FCN 3 command has already been'
     +,/,' executed at the current FCN value.')
*__________________________________________________________________
*--
*--  Command SAVE
*--  ============
*.
      IF(CMD3.EQ.'SAV')THEN
         WRITE(LOUT,11500)
         GO TO 99
      ENDIF
11500 FORMAT(' ***>SAVe',/,
     +' Causes the current parameter values to be saved on a file in',/,
     +' such a format that they can be read in again as Minuit',/,
     +' parameter definitions. If the covariance matrix exists, it is',/
     +,' also output in such a format. The unit number is by default 7,'
     +,/,' or that specified by the user in his call to MINTIO or',/,
     +' MNINIT. The user is responsible for opening the file previous'
     +,/,' to issuing the [SAVe] command (except where this can be done'
     +,/,' interactively).')
*__________________________________________________________________
*--
*--  Command SCAN
*--  ============
*.
      IF(CMD3.EQ.'SCA')THEN
         WRITE(LOUT,11600)
         GO TO 99
      ENDIF
11600 FORMAT(' ***>SCAn  [parno]  [numpts] [from]  [to]',/,
     +' Scans the value of the user function by varying parameter',/,
     +' number [parno], leaving all other parameters fixed at the',/,
     +' current value. If [parno] is not specified, all variable',/,
     +' parameters are scanned in sequence.',/,
     +' The number of points [numpts] in the scan is 40 by default,',/,
     +' and cannot exceed 100. The range of the scan is by default',/,
     +' 2 standard deviations on each side of the current best value,',
     +/,' but can be specified as from [from] to [to].',/,
     +' After each scan, if a new minimum is found, the best parameter'
     +,/,' values are retained as start values for future scans or',/,
     +' minimizations. The curve resulting from each scan is plotted',/
     +,' on the output unit in order to show the approximate behaviour'
     +,/,' of the function.',/,
     +' This command is not intended for minimization, but is sometimes'
     +,/,' useful for debugging the user function or finding a',/,
     +' reasonable starting point.')
*__________________________________________________________________
*--
*--  Command SEEK
*--  ============
*.
      IF(CMD3.EQ.'SEE')THEN
         WRITE(LOUT,11700)
         GO TO 99
      ENDIF
11700 FORMAT(' ***>SEEk  [maxcalls]  [devs]',/,
     +' Causes a Monte Carlo minimization of the function, by choosing',
     +/,' random values of the variable parameters, chosen uniformly',/,
     +' over a hypercube centered at the current best value.',/,
     +' The region size is by default 3 standard deviations on each',/,
     +' side, but can be changed by specifying the value of [devs].')
*__________________________________________________________________
*--
*--  Command SET
*--  ===========
*.
      IF(CMD3.EQ.'SET')THEN
         WRITE(LOUT,11800)
         WRITE(LOUT,11801)
         WRITE(LOUT,11802)
         WRITE(LOUT,11803)
         WRITE(LOUT,11804)
         WRITE(LOUT,11805)
         WRITE(LOUT,11806)
         WRITE(LOUT,11807)
         WRITE(LOUT,11808)
         WRITE(LOUT,11809)
         WRITE(LOUT,11810)
         WRITE(LOUT,11811)
         WRITE(LOUT,11812)
         WRITE(LOUT,11813)
         WRITE(LOUT,11814)
         WRITE(LOUT,11815)
         WRITE(LOUT,11816)
         WRITE(LOUT,11817)
         GO TO 99
      ENDIF
11800 FORMAT(' ***>SET <option_name>',/,/,
     +'  SET BATch',/,
     +'    Informs Minuit that it is running in batch mode.',//,
     +'  SET EPSmachine  <accuracy>',/,
     +'    Informs Minuit that the relative floating point arithmetic',/
     +'    precision is <accuracy>. Minuit determines the nominal',/,
     +'    precision itself, but the SET EPSmachine command can be',/,
     +'    used to override Minuit own determination, when the user',/,
     +'    knows that the FCN function value is not calculated to',/,
     +'    the nominal machine accuracy. Typical values of <accuracy>',/
     +'    are between 10**-5 and 10**-14.')
11801 FORMAT(/,'  SET ERRordef  <up>',/,
     +'    Sets the value of UP (default value= 1.), defining parameter'
     +,/,'    errors. Minuit defines parameter errors as the change',/,
     +'    in parameter value required to change the function value',/,
     +'    by UP. Normally, for chisquared fits UP=1, and for negative'
     +,/,'    log likelihood, UP=0.5.')
11802 FORMAT(/,'   SET GRAdient  [force]',/,
     +'    Informs Minuit that the user function is prepared to',/,
     +'    calculate its own first derivatives and return their values'
     +,/,'    in the array GRAD when IFLAG=2 (see specs of FCN).',/,
     +'    If [force] is not specified, Minuit will calculate',/,
     +'    the FCN derivatives by finite differences at the current',/,
     +'    point and compare with the user calculation at that point,'
     +,/,'    accepting the user values only if they agree.',/,
     +'    If [force]=1, Minuit does not do its own derivative',/,
     +'    calculation, and uses the derivatives calculated in FCN.')
11803 FORMAT(/,'   SET INPut  [unitno]  [filename]',/,
     +'    Causes Minuit, in data-driven mode only, to read subsequent',
     +/,'    commands (or parameter definitions) from a different input'
     +,/,'    file. If no [unitno] is specified, reading reverts to the'
     +,/,'    previous input file, assuming that there was one.',/,
     +'    If [unitno] is specified, and that unit has not been opened,'
     +,/,'    then Minuit attempts to open the file [filename]} if a',/,
     +'    name is specified. If running in interactive mode and',/,
     +'    [filename] is not specified and [unitno] is not opened,',/,
     +'    Minuit prompts the user to enter a file name.',/,
     +'    If the word REWIND is added to the command (note:no blanks',/
     +'    between INPUT and REWIND), the file is rewound before',/,
     +'    reading. Note that this command is implemented in standard',/
     +'    Fortran 77 and the results may depend on the  system;',/,
     +'    for example, if a filename is given under VM/CMS, it must',/,
     +'    be preceeded by a slash.')
11804 FORMAT(/,'   SET INTeractive',/,
     +'    Informs Minuit that it is running interactively.')
11805 FORMAT(/,'   SET LIMits  [parno]  [lolim]  [uplim]',/,
     +'    Allows the user to change the limits on one or all',/,
     +'    parameters. If no arguments are specified, all limits are',/,
     +'    removed from all parameters. If [parno] alone is specified,',
     +/,'    limits are removed from parameter [parno].',/,
     +'    If all arguments are specified, then parameter [parno] will',
     +/,'    be bounded between [lolim] and [uplim].',/,
     +'    Limits can be specified in either order, Minuit will take',/,
     +'    the smaller as [lolim] and the larger as [uplim].',/,
     +'    However, if [lolim] is equal to [uplim], an error condition',
     +/,'    results.')
11806 FORMAT(/,'   SET LINesperpage',/,
     +'     Sets the number of lines for one page of output.',/,
     +'     Default value is 24 for interactive mode')
11807 FORMAT(/,'   SET NOGradient',/,
     +'    The inverse of SET GRAdient, instructs Minuit not to',
     +/,'    use the first derivatives calculated by the user in FCN.')
11808 FORMAT(/,'   SET NOWarnings',/,
     +'    Supresses Minuit warning messages.')
11809 FORMAT(/,'   SET OUTputfile  <unitno>',/,
     +'    Instructs Minuit to write further output to unit <unitno>.')
11810 FORMAT(/,'   SET PAGethrow  <integer>',/,
     +'    Sets the carriage control character for ``new page'' to',/,
     +'    <integer>. Thus the value 1 produces a new page, and 0',/,
     +'    produces a blank line, on some devices (see TOPofpage)')
11811 FORMAT(/,'   SET PARameter  <parno>  <value>',/,
     +'    Sets the value of parameter <parno> to <value>.',/,
     +'    The parameter in question may be variable, fixed, or',/,
     +'    constant, but must be defined.')
11812 FORMAT(/,'   SET PRIntout  <level>',/,
     +'    Sets the print level, determining how much output will be',/,
     +'    produced. Allowed values and their meanings are displayed',/,
     +'    after a SHOw PRInt command, and are currently <level>=:',/,
     +'      [-1]  no output except from SHOW commands',/,
     +'       [0]  minimum output',/,
     +'       [1]  default value, normal output',/,
     +'       [2]  additional output giving intermediate results.',/,
     +'       [3]  maximum output, showing progress of minimizations.',/
     +'    Note: See also the SET WARnings command.')
11813 FORMAT(/,'   SET RANdomgenerator  <seed>',/,
     +'    Sets the seed of the random number generator used in SEEk.',/
     +'    This can be any integer between 10000 and 900000000, for',/,
     +'    example one which was output from a SHOw RANdom command of',/
     +'    a previous run.')
11814 FORMAT(/,'   SET STRategy  <level>',/,
     +'    Sets the strategy to be used in calculating first and second'
     +,/,'    derivatives and in certain minimization methods.',/,
     +'    In general, low values of <level> mean fewer function calls',
     +/,'    and high values mean more reliable minimization.',/,
     +'    Currently allowed values are 0, 1 (default), and 2.')
11815 FORMAT(/,'   SET TITle',/,
     +'    Informs Minuit that the next input line is to be considered',
     +/,'    the (new) title for this task or sub-task.  This is for',/,
     +'    the convenience of the user in reading his output.')
11816 FORMAT(/,'   SET WARnings',/,
     +'    Instructs Minuit to output warning messages when suspicious',
     +/,'    conditions arise which may indicate unreliable results.',/
     +'    This is the default.')
11817 FORMAT(/,'    SET WIDthpage',/,
     +'    Informs Minuit of the output page width.',/,
     +'    Default values are 80 for interactive jobs')
*__________________________________________________________________
*--
*--  Command SHOW
*--  ============
*.
      IF(CMD3.EQ.'SHO')THEN
         WRITE(LOUT,11900)
         WRITE(LOUT,11901)
         WRITE(LOUT,11902)
         WRITE(LOUT,11903)
         WRITE(LOUT,11904)
         GO TO 99
      ENDIF
11900 FORMAT(' ***>SHOw  <option_name>',/,
     +'  All SET XXXX commands have a corresponding SHOw XXXX command.',
     +/,'  In addition, the SHOw commands listed starting here have no',
     +/,'  corresponding SET command for obvious reasons.')
11901 FORMAT(/,'   SHOw CORrelations',/,
     +'    Calculates and prints the parameter correlations from the',/,
     +'    error matrix.')
11902 FORMAT(/,'   SHOw COVariance',/,
     +'    Prints the (external) covariance (error) matrix.')
11903 FORMAT(/,'   SHOw EIGenvalues',/,
     +'    Calculates and prints the eigenvalues of the covariance',/,
     +'    matrix.')
11904 FORMAT(/,'   SHOw FCNvalue',/,
     +'    Prints the current value of FCN.')
*__________________________________________________________________
*--
*--  Command SIMPLEX
*--  ===============
*.
      IF(CMD3.EQ.'SIM')THEN
         WRITE(LOUT,12000)
         GO TO 99
      ENDIF
12000 FORMAT(' ***>SIMplex  [maxcalls]  [tolerance]',/,
     +' Performs a function minimization using the simplex method of',/
     +' Nelder and Mead. Minimization terminates either when the',/,
     +' function has been called (approximately) [maxcalls] times,',/,
     +' or when the estimated vertical distance to minimum (EDM) is',/,
     +' less than [tolerance].',/,
     +' The default value of [tolerance] is 0.1*UP(see SET ERRordef).')
*__________________________________________________________________
*--
*--  Command STANDARD
*--  ================
*.
      IF(CMD3.EQ.'STA')THEN
         WRITE(LOUT,12100)
         GO TO 99
      ENDIF
12100 FORMAT(' ***>STAndard',/,
     +' Causes Minuit to execute the Fortran instruction CALL STAND',/,
     +' where STAND is a subroutine supplied by the user.')
*__________________________________________________________________
*--
*--  Command STOP
*--  ============
*.
      IF(CMD3.EQ.'STO')THEN
         WRITE(LOUT,12200)
         GO TO 99
      ENDIF
12200 FORMAT(' ***>STOP',/,
     +' Same as EXIT.')
*__________________________________________________________________
*--
*--  Command TOPOFPAGE
*--  =================
*.
      IF(CMD3.EQ.'TOP')THEN
         WRITE(LOUT,12300)
         GO TO 99
      ENDIF
12300 FORMAT(' ***>TOPofpage',/,
     +' Causes Minuit to write the character specified in a',/,
     +' SET PAGethrow command (default = 1) to column 1 of the output'
     +,/,' file, which may or may not position your output medium to',
     +/,' the top of a page depending on the device and system.')
*__________________________________________________________________
*
      WRITE(LOUT,13000)
13000 FORMAT(' Unknown MINUIT command. Type HELP for list of commands.')
*
  99  RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNHES1(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC      Called from MNHESS and MNGRAD
CC      Calculate first derivatives (GRD) and uncertainties (DGRD)
CC         and appropriate step sizes GSTEP
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      LOGICAL LDEBUG
      CHARACTER CBF1*22
      LDEBUG = (IDBG(5) .GE. 1)
      IF (ISTRAT .LE. 0) NCYC = 1
      IF (ISTRAT .EQ. 1) NCYC = 2
      IF (ISTRAT .GT. 1) NCYC = 6
      IDRV = 1
      NPARX = NPAR
      DFMIN = 4.*EPSMA2*(ABS(AMIN)+UP)
C                                     main loop over parameters
      DO 100 I= 1, NPAR
      XTF = X(I)
      DMIN = 4.*EPSMA2*ABS(XTF)
      EPSPRI = EPSMA2 + ABS(GRD(I)*EPSMA2)
      OPTSTP = SQRT(DFMIN/(ABS(G2(I))+EPSPRI))
      D = 0.2 * ABS(GSTEP(I))
      IF (D .GT. OPTSTP)  D = OPTSTP
      IF (D .LT. DMIN)  D = DMIN
      CHGOLD = 10000.
C                                       iterate reducing step size
      DO 50 ICYC= 1, NCYC
      X(I) = XTF + D
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTF - D
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTF
C                                       check if step sizes appropriate
      SAG = 0.5*(FS1+FS2-2.0*AMIN)
      GRDOLD = GRD(I)
      GRDNEW = (FS1-FS2)/(2.0*D)
      DGMIN = EPSMAC*(ABS(FS1)+ABS(FS2))/D
      IF (LDEBUG) WRITE (ISYSWR,11) I,IDRV,GSTEP(I),D,G2(I),GRDNEW,SAG
   11 FORMAT (I4,I2,6G12.5)
      IF (GRDNEW .EQ. ZERO)  GO TO 60
      CHANGE = ABS((GRDOLD-GRDNEW)/GRDNEW)
      IF (CHANGE.GT.CHGOLD .AND. ICYC.GT.1)  GO TO 60
      CHGOLD = CHANGE
      GRD(I) = GRDNEW
      GSTEP(I) = SIGN(D,GSTEP(I))
C                  decrease step until first derivative changes by <5%
      IF (CHANGE .LT. 0.05) GO TO 60
      IF (ABS(GRDOLD-GRDNEW) .LT. DGMIN)  GO TO 60
      IF (D .LT. DMIN)  THEN
         CALL MNWARN('D','MNHES1','Step size too small for 1st drv.')
         GO TO 60
      ENDIF
      D = 0.2*D
   50 CONTINUE
C                                       loop satisfied = too many iter
      WRITE (CBF1,'(2G11.3)') GRDOLD,GRDNEW
      CALL MNWARN('D','MNHES1','Too many iterations on D1.'//CBF1)
   60 CONTINUE
      DGRD(I) = MAX(DGMIN,ABS(GRDOLD-GRDNEW))
  100 CONTINUE
C                                        end of first deriv. loop
      CALL MNINEX(X)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNHESS(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Calculates the full second-derivative matrix of FCN
CC        by taking finite differences. When calculating diagonal
CC        elements, it may iterate so that step size is nearly that
CC        which gives function change= UP/10. The first derivatives
CC        of course come as a free side effect, but with a smaller
CC        step size in order to obtain a known accuracy.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION YY(MNI)
      LOGICAL LDEBUG
      CHARACTER CBF1*22
C
      LDEBUG = (IDBG(3) .GE. 1)
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      IF (ISTRAT .LE. 0) THEN
         NCYC = 3
         TLRSTP = 0.5
         TLRG2  = 0.1
      ELSE IF (ISTRAT .EQ. 1) THEN
         NCYC = 5
         TLRSTP = 0.3
         TLRG2  = 0.05
      ELSE
         NCYC = 7
         TLRSTP = 0.1
         TLRG2  = 0.02
      ENDIF
      IF (ISW(5).GE.2 .OR. LDEBUG)  WRITE (ISYSWR,'(A)')
     +   '   START COVARIANCE MATRIX CALCULATION.'
      CFROM = 'HESSE   '
      NFCNFR = NFCN
      CSTATU= 'OK        '
      NPARD = NPAR
C                 make sure starting at the right place
      CALL MNINEX(X)
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FS1 .NE. AMIN) THEN
         DF = AMIN - FS1
         WRITE (CBF1(1:12),'(G12.3)') DF
         CALL MNWARN('D','MNHESS',
     +       'function value differs from AMIN by '//CBF1(1:12) )
      ENDIF
      AMIN = FS1
      IF (LDEBUG) WRITE (ISYSWR,'(A,A)') ' PAR D   GSTEP          ',
     +' D          G2         GRD         SAG    '
C                                        . . . . . . diagonal elements .
C         ISW(2) = 1 if approx, 2 if not posdef, 3 if ok
C         AIMSAG is the sagitta we are aiming for in second deriv calc.
      AIMSAG = SQRT(EPSMA2)*(ABS(AMIN)+UP)
C         Zero the second derivative matrix
      NPAR2 = NPAR*(NPAR+1)/2
      DO 10 I= 1,NPAR2
   10 VHMAT(I) = 0.
C
C         Loop over variable parameters for second derivatives
      IDRV = 2
      DO 100 ID= 1, NPARD
      I = ID + NPAR - NPARD
      IEXT = NEXOFI(I)
      IF (G2(I) .EQ. ZERO) THEN
           WRITE (CBF1(1:4),'(I4)') IEXT
           CALL MNWARN('W','HESSE',
     +      'Second derivative enters zero, param '//CBF1(1:4) )
        WINT = WERR(I)
        IF (NVARL(IEXT) .GT. 1) THEN
           CALL MNDXDI(X(I),I,DXDI)
           IF (ABS(DXDI) .LT. .001) THEN
              WINT = .01
           ELSE
              WINT = WINT/ABS(DXDI)
           ENDIF
        ENDIF
        G2(I) = UP/WINT**2
      ENDIF
      XTF = X(I)
      DMIN = 8.*EPSMA2*ABS(XTF)
C
C                               find step which gives sagitta = AIMSAG
      D = ABS(GSTEP(I))
      DO 40 ICYC= 1, NCYC
C                               loop here only if SAG=0.
      DO 25 MULTPY= 1, 5
C           take two steps
         X(I) = XTF + D
         CALL MNINEX(X)
         NPARX = NPAR
         CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
         NFCN = NFCN + 1
         X(I) = XTF - D
         CALL MNINEX(X)
         CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
         NFCN = NFCN + 1
         X(I) = XTF
         SAG = 0.5*(FS1+FS2-2.0*AMIN)
         IF (SAG .NE. ZERO) GO TO 30
         IF (GSTEP(I) .LT. ZERO) THEN
           IF (D .GE. .5)  GO TO 26
           D = 10.*D
           IF (D .GT. 0.5)  D = 0.51
           GO TO 25
         ENDIF
         D = 10.*D
   25 CONTINUE
   26      WRITE (CBF1(1:4),'(I4)') IEXT
           CALL MNWARN('W','HESSE',
     +      'Second derivative zero for parameter'//CBF1(1:4) )
           GO TO 390
C                             SAG is not zero
   30 G2BFOR = G2(I)
      G2(I) = 2.*SAG/D**2
      GRD(I) = (FS1-FS2)/(2.*D)
      IF (LDEBUG) WRITE (ISYSWR,31) I,IDRV,GSTEP(I),D,G2(I),GRD(I),SAG
   31 FORMAT (I4,I2,6G12.5)
      GSTEP(I) = SIGN(D,GSTEP(I))
      DIRIN(I) = D
      YY(I) = FS1
      DLAST = D
      D = SQRT(2.0*AIMSAG/ABS(G2(I)))
C         if parameter has limits, max int step size = 0.5
      STPINM = 0.5
      IF (GSTEP(I) .LT. ZERO)  D = MIN(D,STPINM)
      IF (D .LT. DMIN)  D = DMIN
C           see if converged
      IF (ABS((D-DLAST)/D)          .LT. TLRSTP)  GO TO 50
      IF (ABS((G2(I)-G2BFOR)/G2(I)) .LT. TLRG2 )  GO TO 50
      D = MIN(D, 10.*DLAST)
      D = MAX(D, 0.1*DLAST)
   40 CONTINUE
C                       end of step size loop
      WRITE (CBF1,'(I2,2E10.2)') IEXT,SAG,AIMSAG
      CALL MNWARN('D','MNHESS','Second Deriv. SAG,AIM= '//CBF1)
C
   50 CONTINUE
      NDEX = I*(I+1)/2
      VHMAT(NDEX) = G2(I)
  100 CONTINUE
C                              end of diagonal second derivative loop
      CALL MNINEX(X)
C                                     refine the first derivatives
      IF (ISTRAT .GT. 0) CALL MNHES1(FCN,FUTIL)
      ISW(2) = 3
      DCOVAR = 0.
C                                        . . . .  off-diagonal elements
      IF (NPAR .EQ. 1)  GO TO 214
      DO 200 I= 1, NPAR
      DO 180 J= 1, I-1
      XTI = X(I)
      XTJ = X(J)
      X(I) = XTI + DIRIN(I)
      X(J) = XTJ + DIRIN(J)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTI
      X(J) = XTJ
      ELEM = (FS1+AMIN-YY(I)-YY(J)) / (DIRIN(I)*DIRIN(J))
      NDEX = I*(I-1)/2 + J
      VHMAT(NDEX) = ELEM
  180 CONTINUE
  200 CONTINUE
  214 CALL MNINEX(X)
C                  verify matrix positive-definite
      CALL MNPSDF
      DO 220 I= 1, NPAR
      DO 219 J= 1, I
      NDEX = I*(I-1)/2 + J
      P(I,J) = VHMAT(NDEX)
  219 P(J,I) = P(I,J)
  220 CONTINUE
      CALL MNVERT(P,MAXINT,MAXINT,NPAR,IFAIL)
      IF (IFAIL .GT. 0)  THEN
        CALL MNWARN('W','HESSE', 'Matrix inversion fails.')
        GO TO 390
      ENDIF
C                                        . . . . . . .  calculate  e d m
      EDM = 0.
        DO 230 I= 1, NPAR
C                              off-diagonal elements
        NDEX = I*(I-1)/2
          DO 225 J= 1, I-1
          NDEX = NDEX + 1
          ZTEMP = 2.0 * P(I,J)
          EDM = EDM + GRD(I)*ZTEMP*GRD(J)
  225     VHMAT(NDEX) = ZTEMP
C                              diagonal elements
        NDEX = NDEX + 1
        VHMAT(NDEX) = 2.0 * P(I,I)
        EDM = EDM  + P(I,I) * GRD(I)**2
  230   CONTINUE
      IF (ISW(5).GE.1 .AND. ISW(2).EQ.3 .AND. ITAUR.EQ.0)
     + WRITE(ISYSWR,'(A)')' COVARIANCE MATRIX CALCULATED SUCCESSFULLY'
      GO TO 900
C                              failure to invert 2nd deriv matrix
  390 ISW(2) = 1
      DCOVAR = 1.
      CSTATU = 'FAILED    '
      IF (ISW(5) .GE. 0) WRITE (ISYSWR,'(A)')
     +        '  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. '
      DO 395 I= 1, NPAR
      NDEX = I*(I-1)/2
      DO 394 J= 1, I-1
      NDEX = NDEX + 1
  394 VHMAT(NDEX) = 0.0
      NDEX = NDEX +1
      G2I = G2(I)
      IF (G2I .LE. ZERO)  G2I = 1.0
  395 VHMAT(NDEX) = 2.0/G2I
  900 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNIMPR(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Attempts to improve on a good local minimum by finding a
CC        better one.   The quadratic part of FCN is removed by MNCALF
CC        and this transformed function is minimized using the simplex
CC        method from several random starting points.
CC        ref. -- Goldstein and Price, Math.Comp. 25, 569 (1971)
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION DSAV(MNI), Y(MNI+1)
      PARAMETER (ALPHA=1.,BETA=0.5,GAMMA=2.0)
      DATA RNUM/0./
      IF (NPAR .LE. 0)  RETURN
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CSTATU = 'UNCHANGED '
      ITAUR = 1
      EPSI = 0.1*UP
      NPFN=NFCN
      NLOOP = WORD7(2)
      IF (NLOOP .LE. 0)  NLOOP = NPAR + 4
      NPARX = NPAR
      NPARP1=NPAR+1
      WG = 1.0/FLOAT(NPAR)
      SIGSAV = EDM
      APSI = AMIN
         DO 2 I= 1, NPAR
         XT(I) = X(I)
         DSAV(I) = WERR(I)
           DO 2 J = 1, I
           NDEX = I*(I-1)/2 + J
           P(I,J) = VHMAT(NDEX)
    2      P(J,I) = P(I,J)
      CALL MNVERT(P,MAXINT,MAXINT,NPAR,IFAIL)
      IF (IFAIL .GE. 1)  GO TO 280
C               Save inverted matrix in VT
         DO 12 I= 1, NPAR
         NDEX = I*(I-1)/2
           DO 12 J= 1, I
           NDEX = NDEX + 1
   12      VTHMAT(NDEX) = P(I,J)
      LOOP = 0
C
   20 CONTINUE
         DO 25 I= 1, NPAR
         DIRIN(I) = 2.0*DSAV(I)
         CALL MNRN15(RNUM,ISEED)
   25    X(I) = XT(I) + 2.0*DIRIN(I)*(RNUM-0.5)
      LOOP = LOOP + 1
      REG = 2.0
      IF (ISW(5) .GE. 0)   WRITE (ISYSWR, 1040) LOOP
   30 CALL  MNCALF(FCN,X,YCALF,FUTIL)
      AMIN = YCALF
C                                        . . . . set up  random simplex
      JL = NPARP1
      JH = NPARP1
      Y(NPARP1) = AMIN
      AMAX = AMIN
         DO 45 I= 1, NPAR
         XI = X(I)
         CALL MNRN15(RNUM,ISEED)
         X(I) = XI - DIRIN(I) *(RNUM-0.5)
         CALL MNCALF(FCN,X,YCALF,FUTIL)
         Y(I) = YCALF
         IF (Y(I) .LT. AMIN)  THEN
            AMIN = Y(I)
            JL = I
         ELSE IF (Y(I) .GT. AMAX)  THEN
            AMAX = Y(I)
            JH = I
         ENDIF
            DO 40 J= 1, NPAR
   40       P(J,I) = X(J)
         P(I,NPARP1) = XI
         X(I) = XI
   45    CONTINUE
C
      EDM = AMIN
      SIG2 = EDM
C                                        . . . . . . .  start main loop
   50 CONTINUE
      IF (AMIN .LT. ZERO)  GO TO 95
      IF (ISW(2) .LE. 2)  GO TO 280
      EP = 0.1*AMIN
      IF (SIG2 .LT. EP   .AND. EDM.LT.EP  )     GO TO 100
      SIG2 = EDM
      IF ((NFCN-NPFN) .GT. NFCNMX)  GO TO 300
C         calculate new point * by reflection
      DO 60 I= 1, NPAR
      PB = 0.
      DO 59 J= 1, NPARP1
   59 PB = PB + WG * P(I,J)
      PBAR(I) = PB - WG * P(I,JH)
   60 PSTAR(I)=(1.+ALPHA)*PBAR(I)-ALPHA*P(I,JH)
      CALL MNCALF(FCN,PSTAR,YCALF,FUTIL)
      YSTAR = YCALF
      IF(YSTAR.GE.AMIN) GO TO 70
C         point * better than jl, calculate new point **
      DO 61 I=1,NPAR
   61 PSTST(I)=GAMMA*PSTAR(I)+(1.-GAMMA)*PBAR(I)
      CALL MNCALF(FCN,PSTST,YCALF,FUTIL)
      YSTST = YCALF
   66 IF (YSTST .LT. Y(JL))  GO TO 67
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      GO TO 50
   67 CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
      GO TO 50
C         point * is not as good as jl
   70 IF (YSTAR .GE. Y(JH))  GO TO 73
      JHOLD = JH
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      IF (JHOLD .NE. JH)  GO TO 50
C         calculate new point **
   73 DO 74 I=1,NPAR
   74 PSTST(I)=BETA*P(I,JH)+(1.-BETA)*PBAR(I)
      CALL MNCALF(FCN,PSTST,YCALF,FUTIL)
      YSTST = YCALF
      IF(YSTST.GT.Y(JH)) GO TO 30
C     point ** is better than jh
      IF (YSTST .LT. AMIN)  GO TO 67
      CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
      GO TO 50
C                                        . . . . . .  end main loop
   95 IF (ISW(5) .GE. 0)  WRITE (ISYSWR,1000)
      REG = 0.1
C                                        . . . . . ask if point is new
  100 CALL MNINEX(X)
      CALL FCN(NPARX,GIN,AMIN,U,4,FUTIL)
      NFCN = NFCN + 1
      DO 120 I= 1, NPAR
      DIRIN(I) = REG*DSAV(I)
      IF (ABS(X(I)-XT(I)) .GT. DIRIN(I)) GO TO 150
  120 CONTINUE
      GO TO 230
  150 NFCNMX = NFCNMX + NPFN - NFCN
      NPFN = NFCN
      CALL MNSIMP(FCN,FUTIL)
      IF (AMIN .GE. APSI)  GO TO 325
      DO 220 I= 1, NPAR
      DIRIN(I) = 0.1 *DSAV(I)
      IF (ABS(X(I)-XT(I)) .GT. DIRIN(I)) GO TO 250
  220 CONTINUE
  230 IF (AMIN .LT. APSI)  GO TO 350
      GO TO 325
C                                        . . . . . . truly new minimum
  250 LNEWMN = .TRUE.
      IF (ISW(2) .GE. 1) THEN
          ISW(2) = 1
          DCOVAR = MAX(DCOVAR,HALF)
      ELSE
          DCOVAR = 1.
      ENDIF
      ITAUR = 0
      NFCNMX = NFCNMX + NPFN - NFCN
      CSTATU = 'NEW MINIMU'
      IF (ISW(5) .GE. 0)      WRITE (ISYSWR,1030)
      RETURN
C                                        . . . return to previous region
  280 IF (ISW(5) .GT. 0) WRITE (ISYSWR,1020)
      GO TO 325
  300 ISW(1) = 1
  325 DO 330 I= 1, NPAR
      DIRIN(I) = 0.01*DSAV(I)
  330 X(I) = XT(I)
      AMIN = APSI
      EDM = SIGSAV
  350 CALL MNINEX(X)
      IF (ISW(5) .GT. 0)    WRITE (ISYSWR,1010)
      CSTATU= 'UNCHANGED '
      CALL MNRSET(0)
      IF (ISW(2) .LT. 2)  GO TO 380
      IF (LOOP .LT. NLOOP .AND. ISW(1) .LT. 1)  GO TO 20
  380 CALL MNPRIN (5,AMIN)
      ITAUR = 0
      RETURN
 1000 FORMAT (54H AN IMPROVEMENT ON THE PREVIOUS MINIMUM HAS BEEN FOUND)
 1010 FORMAT (51H IMPROVE HAS RETURNED TO REGION OF ORIGINAL MINIMUM)
 1020 FORMAT (/44H COVARIANCE MATRIX WAS NOT POSITIVE-DEFINITE)
 1030 FORMAT (/38H IMPROVE HAS FOUND A TRULY NEW MINIMUM/1H ,37(1H*)/)
 1040 FORMAT (/18H START ATTEMPT NO.,I2,  20H TO FIND NEW MINIMUM)
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNINEX(PINT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Transforms from internal coordinates (PINT) to external
CC        parameters (U).   The minimizing routines which work in
CC        internal coordinates call this routine before calling FCN.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION PINT(*)
      DO 100 J= 1, NPAR
      I = NEXOFI(J)
      IF (NVARL(I) .EQ. 1) THEN
         U(I) = PINT(J)
      ELSE
         U(I) = ALIM(I) + 0.5*(SIN(PINT(J)) +1.0) * (BLIM(I)-ALIM(I))
      ENDIF
  100 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.4  1997/09/02 15:16:08  mclareni
* WINNT corrections
*
* Revision 1.3  1997/03/14 17:18:00  mclareni
* WNT mods
*
* Revision 1.2.2.1  1997/01/21 11:33:28  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.2  1996/03/15 18:02:47  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNINIT (I1,I2,I3)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        This is the main initialization subroutine for MINUIT
CC     It initializes some constants in common
CC                (including the logical I/O unit nos.),
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      EXTERNAL INTRAC
      LOGICAL  INTRAC
C            I/O unit numbers
      ISYSRD = I1
      ISYSWR = I2
        ISTKWR(1) = ISYSWR
        NSTKWR = 1
      ISYSSA = I3
      NSTKRD = 0
C               version identifier
      CVRSN = '96.03 '
C               some CONSTANT constants in COMMON
      MAXINT=MNI
      MAXEXT=MNE
      UNDEFI = -54321.
      BIGEDM = 123456.
      CUNDEF = ')UNDEFINED'
      COVMES(0) = 'NO ERROR MATRIX       '
      COVMES(1) = 'ERR MATRIX APPROXIMATE'
      COVMES(2) = 'ERR MATRIX NOT POS-DEF'
      COVMES(3) = 'ERROR MATRIX ACCURATE '
C                some starting values in COMMON
      NBLOCK = 0
      ICOMND = 0
      CTITL = CUNDEF
      CFROM = 'INPUT   '
      NFCNFR = NFCN
      CSTATU= 'INITIALIZE'
      ISW(3) = 0
      ISW(4) = 0
      ISW(5) = 1
C         ISW(6)=0 for batch jobs,  =1 for interactive jobs
C                      =-1 for originally interactive temporarily batch
      ISW(6) = 0
      IF (INTRAC(DUMMY))  ISW(6) = 1
C        DEBUG options set to default values
      DO 10 IDB= 0, MAXDBG
   10 IDBG(IDB) = 0
      LREPOR = .FALSE.
      LWARN  = .TRUE.
      LIMSET = .FALSE.
      LNEWMN = .FALSE.
      ISTRAT = 1
      ITAUR = 0
C        default page dimensions and 'new page' carriage control integer
      NPAGWD = 120
      NPAGLN = 56
      NEWPAG = 1
      IF (ISW(6) .GT. 0) THEN
         NPAGWD = 80
         NPAGLN = 30
         NEWPAG = 0
      ENDIF
      UP = 1.0
      UPDFLT = UP
C                   determine machine accuracy epsmac
      EPSTRY = 0.5
      DO 33 I= 1, 100
      EPSTRY = EPSTRY * 0.5
      EPSP1 = ONE + EPSTRY
      CALL MNTINY(EPSP1, EPSBAK)
      IF (EPSBAK .LT. EPSTRY)  GO TO 35
   33 CONTINUE
      EPSTRY = 1.0E-7
      EPSMAC = 4.0*EPSTRY
      WRITE (ISYSWR,'(A,A,E10.2)') ' MNINIT UNABLE TO DETERMINE',
     + ' ARITHMETIC PRECISION. WILL ASSUME:',EPSMAC
   35 EPSMAC = 8.0 * EPSTRY
      EPSMA2 = 2.0 * SQRT(EPSMAC)
C                 the vlims are a non-negligible distance from pi/2
C         used by MNPINT to set variables "near" the physical limits
      PIBY2 = 2.0*ATAN(1.0)
      DISTNN = 8.0*SQRT(EPSMA2)
      VLIMHI =  PIBY2 - DISTNN
      VLIMLO = -PIBY2 + DISTNN
      CALL MNCLER
      WRITE (ISYSWR,'(3A,I3,A,I3,A,E10.2)')  '  MINUIT RELEASE ',CVRSN,
     +' INITIALIZED.   DIMENSIONS ',MNE,'/',MNI,'  EPSMAC=',EPSMAC
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNINPU(IUNIT,IERR)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC      called by the user to SET INPUT to IUNIT,
CC      an alternative to MNSTIN where the user can specify just
CC      a logical unit number and he is not interrogated about
CC      open files and rewinding, all that is the responsibility
CC      of the user and cannot be fixed interactively.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      IERR = 0
C                              IUNIT = 0, revert to previous input file
      IF (IUNIT .EQ. 0) THEN
        IF (NSTKRD .EQ. 0)  THEN
           WRITE (ISYSWR, '(A)') ' CALL TO MNINPU(0) IGNORED'
           WRITE (ISYSWR, '(A)') ' ALREADY READING FROM PRIMARY INPUT'
        ELSE
          ISYSRD = ISTKRD(NSTKRD)
          NSTKRD = NSTKRD - 1
        ENDIF
C
C                               new input file
      ELSE
          IF (NSTKRD .GE. MAXSTK)  THEN
          WRITE (ISYSWR, '(A)') ' INPUT FILE STACK SIZE EXCEEDED.'
          GO TO 800
          ENDIF
        NSTKRD = NSTKRD + 1
        ISTKRD(NSTKRD) = ISYSRD
        ISYSRD = IUNIT
      ENDIF
C
      RETURN
  800 IERR = 1
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNINTR(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Called by user. Interfaces to MNREAD to allow user to change
CC       easily from Fortran-callable to interactive mode.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      IFLGIN = 3
      CALL MNREAD(FCN,IFLGIN,IFLGUT,FUTIL)
      WRITE (ISYSWR,'(2A/)')  ' END OF MINUIT COMMAND INPUT. ',
     +      '   RETURN TO USER PROGRAM.'
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:48  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNLIMS
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Called from MNSET
CC       Interprets the SET LIM command, to reset the parameter limits
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      CFROM = 'SET LIM '
      NFCNFR = NFCN
      CSTATU= 'NO CHANGE '
      I2 = WORD7(1)
      IF (I2 .GT. MAXEXT .OR. I2 .LT. 0)  GO TO 900
      IF (I2 .GT. 0)  GO TO 30
C                                     set limits on all parameters
      NEWCOD = 4
      IF (WORD7(2) .EQ. WORD7(3))  NEWCOD = 1
      DO 20 INU= 1, NU
      IF (NVARL(INU) .LE. 0)  GO TO 20
      IF (NVARL(INU).EQ.1 .AND. NEWCOD.EQ.1)  GO TO 20
      KINT = NIOFEX(INU)
C             see if parameter has been fixed
      IF (KINT .LE. 0)  THEN
         IF (ISW(5) .GE. 0)  WRITE (ISYSWR,'(11X,A,I3)')
     +      ' LIMITS NOT CHANGED FOR FIXED PARAMETER:',INU
         GO TO 20
      ENDIF
      IF (NEWCOD .EQ. 1)  THEN
C            remove limits from parameter
         IF (ISW(5) .GT. 0)     WRITE (ISYSWR,134)  INU
         CSTATU = 'NEW LIMITS'
         CALL MNDXDI(X(KINT),KINT,DXDI)
         SNEW = GSTEP(KINT)*DXDI
         GSTEP(KINT) = ABS(SNEW)
         NVARL(INU) = 1
      ELSE
C             put limits on parameter
         ALIM(INU) = MIN(WORD7(2),WORD7(3))
         BLIM(INU) = MAX(WORD7(2),WORD7(3))
         IF (ISW(5) .GT. 0) WRITE (ISYSWR,237)  INU,ALIM(INU),BLIM(INU)
         NVARL(INU) = 4
         CSTATU = 'NEW LIMITS'
         GSTEP(KINT) = -0.1
      ENDIF
   20 CONTINUE
      GO TO 900
C                                       set limits on one parameter
   30 IF (NVARL(I2) .LE. 0)  THEN
        WRITE (ISYSWR,'(A,I3,A)') ' PARAMETER ',I2,' IS NOT VARIABLE.'
        GO TO 900
      ENDIF
      KINT = NIOFEX(I2)
C                                       see if parameter was fixed
      IF (KINT .EQ. 0)  THEN
         WRITE (ISYSWR,'(A,I3)')
     +     ' REQUEST TO CHANGE LIMITS ON FIXED PARAMETER:',I2
         DO 82 IFX= 1, NPFIX
         IF (I2 .EQ. IPFIX(IFX)) GO TO 92
   82    CONTINUE
         WRITE (ISYSWR,'(A)') ' MINUIT BUG IN MNLIMS. SEE F. JAMES'
   92    CONTINUE
      ENDIF
      IF (WORD7(2) .NE. WORD7(3))  GO TO 235
C                                       remove limits
      IF (NVARL(I2) .NE. 1)  THEN
         IF (ISW(5) .GT. 0)  WRITE (ISYSWR,134)  I2
  134    FORMAT (30H LIMITS REMOVED FROM PARAMETER  ,I4)
         CSTATU = 'NEW LIMITS'
         IF (KINT .LE. 0)  THEN
            GSTEPS(IFX) = ABS(GSTEPS(IFX))
         ELSE
            CALL MNDXDI(X(KINT),KINT,DXDI)
            IF (ABS(DXDI) .LT. 0.01)  DXDI=0.01
            GSTEP(KINT) = ABS(GSTEP(KINT)*DXDI)
            GRD(KINT) = GRD(KINT)*DXDI
         ENDIF
         NVARL(I2) = 1
      ELSE
         WRITE (ISYSWR,'(A,I3)') ' NO LIMITS SPECIFIED.  PARAMETER ',
     +        I2,' IS ALREADY UNLIMITED.  NO CHANGE.'
      ENDIF
      GO TO 900
C                                        put on limits
  235 ALIM(I2) = MIN(WORD7(2),WORD7(3))
      BLIM(I2) = MAX(WORD7(2),WORD7(3))
      NVARL(I2) = 4
      IF (ISW(5) .GT. 0)   WRITE (ISYSWR,237)  I2,ALIM(I2),BLIM(I2)
  237 FORMAT (10H PARAMETER ,I3, 14H LIMITS SET TO  ,2G15.5)
      CSTATU = 'NEW LIMITS'
      IF (KINT .LE. 0)  THEN
         GSTEPS(IFX) = -0.1
      ELSE
         GSTEP(KINT) = -0.1
      ENDIF
C
  900 CONTINUE
      IF (CSTATU .NE. 'NO CHANGE ')  THEN
        CALL MNEXIN(X)
        CALL MNRSET(1)
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNLINE(FCN,START,FSTART,STEP,SLOPE,TOLER,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Perform a line search from position START
CC        along direction STEP, where the length of vector STEP
CC                   gives the expected position of minimum.
CC        FSTART is value of function at START
CC        SLOPE (if non-zero) is df/dx along STEP at START
CC        TOLER is initial tolerance of minimum in direction STEP
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION START(*), STEP(*)
      PARAMETER (MAXPT=12)
      DIMENSION XPQ(MAXPT),YPQ(MAXPT)
      CHARACTER*1 CHPQ(MAXPT)
      DIMENSION XVALS(3),FVALS(3),COEFF(3)
      CHARACTER*26 CHARAL
      CHARACTER*60 CMESS
      PARAMETER (SLAMBG=5.,ALPHA=2.)
C SLAMBG and ALPHA control the maximum individual steps allowed.
C The first step is always =1. The max length of second step is SLAMBG.
C The max size of subsequent steps is the maximum previous successful
C   step multiplied by ALPHA + the size of most recent successful step,
C   but cannot be smaller than SLAMBG.
      LOGICAL LDEBUG
      DATA CHARAL / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      LDEBUG = (IDBG(1).GE.1)
C                  starting values for overall limits on total step SLAM
      OVERAL = 1000.
      UNDRAL = -100.
C                              debug check if start is ok
      IF (LDEBUG)  THEN
         CALL MNINEX(START)
         CALL FCN(NPARX,GIN,F1,U,4,FUTIL)
         NFCN=NFCN+1
         IF (F1 .NE. FSTART) THEN
             WRITE (ISYSWR,'(A/2E14.5/2X,10F10.5)')
     + ' MNLINE start point not consistent, F values, parameters=',
     +  (X(KK),KK=1,NPAR)
         ENDIF
      ENDIF
C                                      . set up linear search along STEP
      FVMIN = FSTART
      XVMIN = ZERO
      NXYPT = 1
      CHPQ(1) = CHARAL(1:1)
      XPQ(1) = 0.
      YPQ(1) = FSTART
C               SLAMIN = smallest possible value of ABS(SLAM)
      SLAMIN = ZERO
      DO 20 I= 1, NPAR
      IF (STEP(I) .EQ. ZERO)  GO TO 20
      RATIO = ABS(START(I)/STEP(I))
      IF (SLAMIN .EQ. ZERO)     SLAMIN = RATIO
      IF (RATIO .LT. SLAMIN)  SLAMIN = RATIO
   20 X(I) = START(I) + STEP(I)
      IF (SLAMIN .EQ. ZERO)  SLAMIN = EPSMAC
      SLAMIN = SLAMIN*EPSMA2
      NPARX = NPAR
C
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,F1,U,4,FUTIL)
      NFCN=NFCN+1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = 1.
      YPQ(NXYPT) = F1
      IF (F1 .LT. FSTART) THEN
         FVMIN = F1
         XVMIN = 1.0
      ENDIF
C                         . quadr interp using slope GDEL and two points
      SLAM = 1.
      TOLER8 = TOLER
      SLAMAX = SLAMBG
      FLAST = F1
C                         can iterate on two-points (cut) if no imprvmnt
   25 CONTINUE
      DENOM = 2.0*(FLAST-FSTART-SLOPE*SLAM)/SLAM**2
C     IF (DENOM .EQ. ZERO)  DENOM = -0.1*SLOPE
                            SLAM  = 1.
      IF (DENOM .NE. ZERO)  SLAM = -SLOPE/DENOM
      IF (SLAM  .LT. ZERO)  SLAM = SLAMAX
      IF (SLAM .GT. SLAMAX)  SLAM = SLAMAX
      IF (SLAM .LT. TOLER8)  SLAM = TOLER8
      IF (SLAM .LT. SLAMIN)  GO TO 80
      IF (ABS(SLAM-1.0).LT.TOLER8 .AND. F1.LT.FSTART)  GO TO 70
      IF (ABS(SLAM-1.0).LT.TOLER8) SLAM = 1.0+TOLER8
      IF (NXYPT .GE. MAXPT) GO TO 65
      DO 30 I= 1, NPAR
   30 X(I) = START(I) + SLAM*STEP(I)
      CALL MNINEX(X)
      CALL FCN(NPAR,GIN,F2,U,4,FUTIL)
      NFCN = NFCN + 1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = SLAM
      YPQ(NXYPT) = F2
      IF (F2 .LT. FVMIN)  THEN
         FVMIN = F2
         XVMIN = SLAM
      ENDIF
      IF (FSTART .EQ. FVMIN) THEN
         FLAST = F2
         TOLER8 = TOLER*SLAM
         OVERAL = SLAM-TOLER8
         SLAMAX = OVERAL
         GO TO 25
      ENDIF
C                                        . quadr interp using 3 points
      XVALS(1) = XPQ(1)
      FVALS(1) = YPQ(1)
      XVALS(2) = XPQ(NXYPT-1)
      FVALS(2) = YPQ(NXYPT-1)
      XVALS(3) = XPQ(NXYPT)
      FVALS(3) = YPQ(NXYPT)
C                             begin iteration, calculate desired step
   50 CONTINUE
      SLAMAX = MAX(SLAMAX,ALPHA*ABS(XVMIN))
      CALL MNPFIT(XVALS,FVALS,3,COEFF,SDEV)
      IF (COEFF(3) .LE. ZERO)  THEN
         SLOPEM = 2.0*COEFF(3)*XVMIN + COEFF(2)
         IF (SLOPEM .LE. ZERO) THEN
            SLAM = XVMIN + SLAMAX
         ELSE
            SLAM = XVMIN - SLAMAX
         ENDIF
      ELSE
         SLAM = -COEFF(2)/(2.0*COEFF(3))
         IF (SLAM .GT. XVMIN+SLAMAX)  SLAM = XVMIN+SLAMAX
         IF (SLAM .LT. XVMIN-SLAMAX)  SLAM = XVMIN-SLAMAX
      ENDIF
      IF (SLAM .GT. ZERO) THEN
          IF (SLAM .GT. OVERAL) SLAM = OVERAL
      ELSE
          IF (SLAM .LT. UNDRAL) SLAM = UNDRAL
      ENDIF
C               come here if step was cut below
   52 CONTINUE
      TOLER9 = MAX(TOLER8,ABS(TOLER8*SLAM))
      DO 55 IPT= 1, 3
      IF (ABS(SLAM-XVALS(IPT)) .LT. TOLER9)  GO TO 70
   55 CONTINUE
C                take the step
      IF (NXYPT .GE. MAXPT) GO TO 65
      DO 60 I= 1, NPAR
   60 X(I) = START(I)+SLAM*STEP(I)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,F3,U,4,FUTIL)
      NFCN = NFCN + 1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = SLAM
      YPQ(NXYPT) = F3
C             find worst previous point out of three
      FVMAX = FVALS(1)
      NVMAX = 1
      IF (FVALS(2) .GT. FVMAX) THEN
         FVMAX = FVALS(2)
         NVMAX = 2
      ENDIF
      IF (FVALS(3) .GT. FVMAX) THEN
         FVMAX = FVALS(3)
         NVMAX = 3
      ENDIF
C              if latest point worse than all three previous, cut step
      IF (F3 .GE. FVMAX)  THEN
          IF (NXYPT .GE. MAXPT) GO TO 65
          IF (SLAM .GT. XVMIN) OVERAL = MIN(OVERAL,SLAM-TOLER8)
          IF (SLAM .LT. XVMIN) UNDRAL = MAX(UNDRAL,SLAM+TOLER8)
          SLAM = 0.5*(SLAM+XVMIN)
          GO TO 52
      ENDIF
C              prepare another iteration, replace worst previous point
      XVALS(NVMAX) = SLAM
      FVALS(NVMAX) = F3
      IF (F3 .LT. FVMIN)  THEN
         FVMIN = F3
         XVMIN = SLAM
      ELSE
         IF (SLAM .GT. XVMIN) OVERAL = MIN(OVERAL,SLAM-TOLER8)
         IF (SLAM .LT. XVMIN) UNDRAL = MAX(UNDRAL,SLAM+TOLER8)
      ENDIF
      IF (NXYPT .LT. MAXPT)  GO TO 50
C                                            . . end of iteration . . .
C            stop because too many iterations
   65 CMESS = ' LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS '
      IF (LDEBUG) THEN
        WRITE (ISYSWR,'(A/(2X,6G12.4))') ' MNLINE DEBUG: steps=',
     +    (STEP(KK),KK=1,NPAR)
      ENDIF
      GO TO 100
C            stop because within tolerance
   70 CONTINUE
      CMESS = ' LINE SEARCH HAS ATTAINED TOLERANCE '
      GO TO 100
   80 CONTINUE
      CMESS = ' STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM'
  100 CONTINUE
      AMIN = FVMIN
      DO 120 I= 1, NPAR
      DIRIN(I) = STEP(I)*XVMIN
  120 X(I) = START(I) + DIRIN(I)
      CALL MNINEX(X)
      IF (XVMIN .LT. 0.)      CALL MNWARN('D','MNLINE',
     +                   ' LINE MINIMUM IN BACKWARDS DIRECTION')
      IF (FVMIN .EQ. FSTART)  CALL MNWARN('D','MNLINE',
     +                     ' LINE SEARCH FINDS NO IMPROVEMENT ')
      IF (LDEBUG)  THEN
         WRITE (ISYSWR,'('' AFTER'',I3,'' POINTS,'',A)') NXYPT,CMESS
         CALL MNPLOT(XPQ,YPQ,CHPQ,NXYPT,ISYSWR,NPAGWD,NPAGLN)
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMATU(KODE)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        prints the covariance matrix v when KODE=1.
CC        always prints the global correlations, and
CC        calculates and prints the individual correlation coefficients
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION VLINE(MNI)
      ISW2 = ISW(2)
      IF (ISW2 .LT. 1)  THEN
          WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
          GO TO 500
      ENDIF
      IF (NPAR .EQ. 0)  THEN
          WRITE (ISYSWR,'('' MNMATU: NPAR=0'')')
          GO TO 500
          ENDIF
C                                       . . . . .external error matrix
      IF (KODE .EQ. 1)  THEN
         ISW5 = ISW(5)
         ISW(5) = 2
         CALL MNEMAT(P,MAXINT)
           IF (ISW2.LT.3)  WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
         ISW(5) = ISW5
      ENDIF
C                                       . . . . . correlation coeffs. .
      IF (NPAR .LE. 1)   GO TO 500
      CALL MNWERR
C     NCOEF is number of coeff. that fit on one line, not to exceed 20
      NCOEF = (NPAGWD-19)/6
      NCOEF = MIN(NCOEF,20)
      NPARM = MIN(NPAR,NCOEF)
      WRITE (ISYSWR, 150) (NEXOFI(ID),ID=1,NPARM)
  150 FORMAT (/36H PARAMETER  CORRELATION COEFFICIENTS  /
     +         18H       NO.  GLOBAL   ,20I6)
      DO 200 I= 1, NPAR
         IX = NEXOFI(I)
         NDI = I*(I+1)/2
           DO 170 J= 1, NPAR
           M = MAX(I,J)
           N = MIN(I,J)
           NDEX = M*(M-1)/2 + N
           NDJ = J*(J+1)/2
  170      VLINE(J) = VHMAT(NDEX)/SQRT(ABS(VHMAT(NDI)*VHMAT(NDJ)))
         NPARM = MIN(NPAR,NCOEF)
         WRITE (ISYSWR,171)   IX, GLOBCC(I), (VLINE(IT),IT=1,NPARM)
  171    FORMAT (6X,I3,2X,F7.5,1X,20F6.3)
         IF (I.LE.NPARM) GO TO 200
            DO 190 ISO= 1, 10
            NSOFAR = NPARM
            NPARM = MIN(NPAR,NSOFAR+NCOEF)
            WRITE (ISYSWR,181)  (VLINE(IT),IT=NSOFAR+1,NPARM)
  181       FORMAT (19X,20F6.3)
            IF (I .LE. NPARM) GO TO 192
  190       CONTINUE
  192    CONTINUE
  200 CONTINUE
      IF (ISW2.LT.3)  WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
  500 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:49  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMIGR(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Performs a local function minimization using basically the
CC        method of Davidon-Fletcher-Powell as modified by Fletcher
CC        ref. -- Fletcher, Comp.J. 13,317 (1970)   "switching method"
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION GS(MNI), STEP(MNI),  XXS(MNI), FLNU(MNI), VG(MNI)
      LOGICAL LDEBUG
      PARAMETER (TOLER=0.05)
      IF (NPAR .LE. 0)  RETURN
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      LDEBUG = (IDBG(4) .GE. 1)
      CFROM = 'MIGRAD  '
      NFCNFR = NFCN
      NFCNMG = NFCN
      CSTATU= 'INITIATE  '
      ISWTR = ISW(5) - 2*ITAUR
      NPFN = NFCN
      NPARX = NPAR
      LENV = NPAR*(NPAR+1)/2
      NRSTRT = 0
      NPSDF = 0
      LINED2 = 0
      ISW(4) = -1
      RHOTOL = 1.0E-3*APSI
      IF (ISWTR .GE. 1)  WRITE (ISYSWR,470) ISTRAT,RHOTOL
  470 FORMAT (' START MIGRAD MINIMIZATION.  STRATEGY',I2,
     +'.  CONVERGENCE WHEN EDM .LT.',E9.2)
C                                           initialization strategy
      IF (ISTRAT.LT.2 .OR. ISW(2).GE.3)  GO TO 2
C                                come (back) here to restart completely
    1 CONTINUE
      IF (NRSTRT .GT. ISTRAT)  THEN
         CSTATU= 'FAILED    '
         ISW(4) = -1
         GO TO 230
         ENDIF
C                                      . get full covariance and gradient
      CALL MNHESS(FCN,FUTIL)
      CALL MNWERR
      NPSDF = 0
      IF (ISW(2) .GE. 1)  GO TO 10
C                                        . get gradient at start point
    2 CONTINUE
      CALL MNINEX(X)
      IF (ISW(3) .EQ. 1) THEN
          CALL FCN(NPARX,GIN,FZERO,U,2,FUTIL)
          NFCN = NFCN + 1
      ENDIF
      CALL MNDERI(FCN,FUTIL)
      IF (ISW(2) .GE. 1)  GO TO 10
C                                   sometimes start with diagonal matrix
      DO 3 I= 1, NPAR
         XXS(I) = X(I)
         STEP(I) = ZERO
    3 CONTINUE
C                           do line search if second derivative negative
      LINED2 = LINED2 + 1
      IF (LINED2 .LT. (ISTRAT+1)*NPAR) THEN
      DO 5 I= 1, NPAR
         IF (G2(I) .GT. ZERO)  GO TO 5
         STEP(I) = -SIGN(GSTEP(I),GRD(I))
         GDEL = STEP(I)*GRD(I)
         FS = AMIN
         CALL MNLINE(FCN,XXS,FS,STEP,GDEL,TOLER,FUTIL)
         CALL MNWARN('D','MNMIGR','Negative G2 line search')
         IEXT = NEXOFI(I)
         IF (LDEBUG) WRITE (ISYSWR,'(A,I3,2G13.3)')
     +    ' Negative G2 line search, param ',IEXT,FS,AMIN
         GO TO 2
    5 CONTINUE
      ENDIF
C                           make diagonal error matrix
      DO 8 I=1,NPAR
         NDEX = I*(I-1)/2
           DO 7 J=1,I-1
           NDEX = NDEX + 1
    7      VHMAT(NDEX) = 0.
         NDEX = NDEX + 1
         IF (G2(I) .LE. ZERO)  G2(I) = 1.
         VHMAT(NDEX) = 2./G2(I)
    8 CONTINUE
      DCOVAR = 1.
      IF (LDEBUG) WRITE (ISYSWR,'(A,A/(1X,10G10.2))') ' DEBUG MNMIGR,',
     +  ' STARTING MATRIX DIAGONAL,  VHMAT=', (VHMAT(KK),KK=1,LENV)
C                                         ready to start first iteration
   10 CONTINUE
      NRSTRT = NRSTRT + 1
      IF (NRSTRT .GT. ISTRAT+1)  THEN
         CSTATU= 'FAILED    '
         GO TO 230
         ENDIF
      FS = AMIN
C                                        . . . get EDM and set up loop
      EDM = 0.
         DO 18 I= 1, NPAR
         GS(I) = GRD(I)
         XXS(I) = X(I)
         NDEX = I*(I-1)/2
           DO 17 J= 1, I-1
           NDEX = NDEX + 1
   17      EDM = EDM + GS(I)*VHMAT(NDEX)*GS(J)
         NDEX = NDEX + 1
   18    EDM = EDM + 0.5 * GS(I)**2 *VHMAT(NDEX)
      EDM = EDM * 0.5 * (1.0+3.0*DCOVAR)
        IF (EDM .LT. ZERO)  THEN
        CALL MNWARN('W','MIGRAD','STARTING MATRIX NOT POS-DEFINITE.')
        ISW(2) = 0
        DCOVAR = 1.
        GO TO 2
        ENDIF
      IF (ISW(2) .EQ. 0)  EDM=BIGEDM
      ITER = 0
      CALL MNINEX(X)
      CALL MNWERR
      IF (ISWTR .GE. 1)  CALL MNPRIN(3,AMIN)
      IF (ISWTR .GE. 2)  CALL MNMATU(0)
C                                        . . . . .  start main loop
   24 CONTINUE
      IF (NFCN-NPFN .GE. NFCNMX)  GO TO 190
      GDEL = 0.
      GSSQ = 0.
         DO 30  I=1,NPAR
         RI = 0.
         GSSQ = GSSQ + GS(I)**2
           DO 25 J=1,NPAR
           M = MAX(I,J)
           N = MIN(I,J)
           NDEX = M*(M-1)/2 + N
   25      RI = RI + VHMAT(NDEX) *GS(J)
         STEP(I) = -0.5*RI
   30    GDEL = GDEL + STEP(I)*GS(I)
      IF (GSSQ .EQ. ZERO)  THEN
          CALL MNWARN('D','MIGRAD',
     +             ' FIRST DERIVATIVES OF FCN ARE ALL ZERO')
          GO TO 300
      ENDIF
C                 if gdel positive, V not posdef
      IF (GDEL .GE. ZERO)  THEN
         CALL MNWARN('D','MIGRAD',' NEWTON STEP NOT DESCENT.')
         IF (NPSDF .EQ. 1)  GO TO 1
         CALL MNPSDF
         NPSDF = 1
         GO TO 24
         ENDIF
C                                        . . . . do line search
      CALL MNLINE(FCN,XXS,FS,STEP,GDEL,TOLER,FUTIL)
      IF (AMIN .EQ. FS) GO TO 200
      CFROM  = 'MIGRAD  '
      NFCNFR = NFCNMG
      CSTATU= 'PROGRESS  '
C                                        . get gradient at new point
      CALL MNINEX(X)
      IF (ISW(3) .EQ. 1) THEN
          CALL FCN(NPARX,GIN,FZERO,U,2,FUTIL)
          NFCN = NFCN + 1
      ENDIF
      CALL MNDERI(FCN,FUTIL)
C                                         . calculate new EDM
      NPSDF = 0
   81 EDM = 0.
      GVG = 0.
      DELGAM = 0.
      GDGSSQ = 0.
         DO 100 I= 1, NPAR
         RI = 0.
         VGI = 0.
           DO 90 J= 1, NPAR
           M = MAX(I,J)
           N = MIN(I,J)
           NDEX = M*(M-1)/2 + N
           VGI = VGI + VHMAT(NDEX)*(GRD(J)-GS(J))
   90      RI  =  RI + VHMAT(NDEX)* GRD(J)
      VG(I) = VGI*0.5
      GAMI = GRD(I) - GS(I)
      GDGSSQ = GDGSSQ + GAMI**2
      GVG = GVG + GAMI*VG(I)
      DELGAM = DELGAM + DIRIN(I)*GAMI
  100 EDM = EDM + GRD(I)*RI*0.5
      EDM = EDM * 0.5 * (1.0 + 3.0*DCOVAR)
C                          . if EDM negative,  not positive-definite
      IF (EDM .LT. ZERO .OR. GVG .LE. ZERO)  THEN
         CALL MNWARN('D','MIGRAD','NOT POS-DEF. EDM OR GVG NEGATIVE.')
         CSTATU = 'NOT POSDEF'
         IF (NPSDF .EQ. 1)  GO TO 230
         CALL MNPSDF
         NPSDF = 1
         GO TO 81
      ENDIF
C                            print information about this iteration
      ITER = ITER + 1
      IF (ISWTR.GE.3 .OR. (ISWTR.EQ.2.AND.MOD(ITER,10).EQ.1)) THEN
         CALL MNWERR
         CALL MNPRIN(3,AMIN)
      ENDIF
      IF (GDGSSQ .EQ. ZERO)  CALL MNWARN('D','MIGRAD',
     +           'NO CHANGE IN FIRST DERIVATIVES OVER LAST STEP')
      IF (DELGAM .LT. ZERO) CALL MNWARN('D','MIGRAD',
     +          'FIRST DERIVATIVES INCREASING ALONG SEARCH LINE')
C                                        .  update covariance matrix
      CSTATU = 'IMPROVEMNT'
        IF (LDEBUG) WRITE (ISYSWR,'(A,(1X,10G10.3))') ' VHMAT 1 =',
     +             (VHMAT(KK),KK=1,10)
      DSUM = 0.
      VSUM = 0.
         DO  120  I=1, NPAR
           DO  120  J=1, I
           D = DIRIN(I)*DIRIN(J)/DELGAM - VG(I)*VG(J)/GVG
           DSUM = DSUM + ABS(D)
           NDEX = I*(I-1)/2 + J
           VHMAT(NDEX) = VHMAT(NDEX) + 2.0*D
           VSUM = VSUM + ABS(VHMAT(NDEX))
  120      CONTINUE
C                smooth local fluctuations by averaging DCOVAR
      DCOVAR = 0.5*(DCOVAR + DSUM/VSUM)
      IF (ISWTR.GE.3 .OR. LDEBUG) WRITE (ISYSWR,'(A,F5.1,A)')
     +      ' RELATIVE CHANGE IN COV. MATRIX=',DCOVAR*100.,'%'
      IF (LDEBUG) WRITE (ISYSWR,'(A,(1X,10G10.3))') ' VHMAT 2 =',
     +             (VHMAT(KK),KK=1,10)
      IF (DELGAM .LE. GVG)  GO TO 135
      DO 125 I= 1, NPAR
  125 FLNU(I) = DIRIN(I)/DELGAM - VG(I)/GVG
      DO 130 I= 1, NPAR
      DO 130 J= 1, I
      NDEX = I*(I-1)/2 + J
  130 VHMAT(NDEX) = VHMAT(NDEX) + 2.0*GVG*FLNU(I)*FLNU(J)
  135 CONTINUE
C                                              and see if converged
      IF (EDM .LT. 0.1*RHOTOL)  GO TO 300
C                                    if not, prepare next iteration
      DO 140 I= 1, NPAR
      XXS(I) = X(I)
      GS(I) = GRD(I)
  140 CONTINUE
      FS = AMIN
      IF (ISW(2) .EQ. 0  .AND. DCOVAR.LT. 0.5 )  ISW(2) = 1
      IF (ISW(2) .EQ. 3  .AND. DCOVAR.GT. 0.1 )  ISW(2) = 1
      IF (ISW(2) .EQ. 1  .AND. DCOVAR.LT. 0.05)  ISW(2) = 3
      GO TO 24
C                                        . . . . .  end main loop
C                                         . . call limit in MNMIGR
  190 ISW(1) = 1
      IF (ISW(5) .GE. 0)
     +     WRITE (ISYSWR,'(A)')  ' CALL LIMIT EXCEEDED IN MIGRAD.'
      CSTATU = 'CALL LIMIT'
      GO TO 230
C                                         . . fails to improve . .
  200 IF (ISWTR .GE. 1)  WRITE (ISYSWR,'(A)')
     +           ' MIGRAD FAILS TO FIND IMPROVEMENT'
      DO 210 I= 1, NPAR
  210 X(I) = XXS(I)
      IF (EDM .LT. RHOTOL)  GO TO 300
      IF (EDM .LT. ABS(EPSMA2*AMIN))  THEN
         IF (ISWTR .GE. 0)  WRITE (ISYSWR, '(A)')
     +      ' MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT.'
         GO TO 300
         ENDIF
      IF (ISTRAT .LT. 1)  THEN
         IF (ISW(5) .GE. 0) WRITE (ISYSWR, '(A)')
     +    ' MIGRAD FAILS WITH STRATEGY=0.   WILL TRY WITH STRATEGY=1.'
         ISTRAT = 1
      ENDIF
         GO TO 1
C                                         . . fails to converge
  230 IF (ISWTR .GE. 0)  WRITE (ISYSWR,'(A)')
     +    ' MIGRAD TERMINATED WITHOUT CONVERGENCE.'
      IF (ISW(2) .EQ. 3)  ISW(2) = 1
      ISW(4) = -1
      GO TO 400
C                                         . . apparent convergence
  300 IF (ISWTR .GE. 0) WRITE(ISYSWR,'(/A)')
     +   ' MIGRAD MINIMIZATION HAS CONVERGED.'
      IF (ITAUR .EQ. 0) THEN
        IF (ISTRAT .GE. 2 .OR. (ISTRAT.EQ.1.AND.ISW(2).LT.3)) THEN
           IF (ISW(5) .GE. 0)  WRITE (ISYSWR, '(/A)')
     +      ' MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.'
           CALL MNHESS(FCN,FUTIL)
           CALL MNWERR
           NPSDF = 0
           IF (EDM .GT. RHOTOL) GO TO 10
        ENDIF
      ENDIF
      CSTATU='CONVERGED '
      ISW(4) = 1
C                                           come here in any case
  400 CONTINUE
      CFROM = 'MIGRAD  '
      NFCNFR = NFCNMG
      CALL  MNINEX(X)
      CALL MNWERR
      IF (ISWTR .GE. 0)  CALL MNPRIN (3,AMIN)
      IF (ISWTR .GE. 1)  CALL MNMATU(1)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMNOS(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Performs a MINOS error analysis on those parameters for
CC        which it is requested on the MINOS command by calling 
CC        MNMNOT for each parameter requested.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      IF (NPAR .LE. 0)  GO TO 700
      NGOOD = 0
      NBAD = 0
      NFCNMI = NFCN
C                                      . loop over parameters requested
      DO 570 KNT= 1, NPAR
      IF (INT(WORD7(2)) .EQ. 0) THEN
          ILAX = NEXOFI(KNT)
      ELSE
          IF (KNT .GE. 7)  GO TO 580
          ILAX = INT(WORD7(KNT+1))
          IF (ILAX .EQ. 0)  GO TO 580
          IF (ILAX .GT. 0 .AND. ILAX .LE. NU) THEN
             IF (NIOFEX(ILAX) .GT. 0)  GO TO 565
          ENDIF
          WRITE (ISYSWR,564) ILAX
  564     FORMAT (' PARAMETER NUMBER ',I5,' NOT VARIABLE. IGNORED.')
          GO TO 570
      ENDIF
  565 CONTINUE
C                                         calculate one pair of M E's
      ILAX2 = 0
      CALL MNMNOT(FCN,ILAX,ILAX2,VAL2PL,VAL2MI,FUTIL)
      IF (LNEWMN)  GO TO 650
C                                          update NGOOD and NBAD
      IIN = NIOFEX(ILAX)
      IF (ERP(IIN) .GT. ZERO) THEN
         NGOOD=NGOOD+1
      ELSE
         NBAD=NBAD+1
      ENDIF
      IF (ERN(IIN) .LT. ZERO) THEN
         NGOOD=NGOOD+1
      ELSE
         NBAD=NBAD+1
      ENDIF
  570 CONTINUE
C                                           end of loop . . . . . . .
  580 CONTINUE
C                                        . . . . printout final values .
      CFROM = 'MINOS   '
      NFCNFR = NFCNMI
      CSTATU= 'UNCHANGED '
      IF (NGOOD.EQ.0.AND.NBAD.EQ.0) GO TO 700
      IF (NGOOD.GT.0.AND.NBAD.EQ.0) CSTATU='SUCCESSFUL'
      IF (NGOOD.EQ.0.AND.NBAD.GT.0) CSTATU='FAILURE   '
      IF (NGOOD.GT.0.AND.NBAD.GT.0) CSTATU='PROBLEMS  '
      IF (ISW(5) .GE. 0) CALL MNPRIN(4,AMIN)
      IF (ISW(5) .GE. 2) CALL MNMATU(0)
      GO TO 900
C                                        . . . new minimum found . . . .
  650 CONTINUE
      CFROM = 'MINOS   '
      NFCNFR = NFCNMI
      CSTATU= 'NEW MINIMU'
      IF (ISW(5) .GE. 0) CALL MNPRIN(4,AMIN)
      WRITE (ISYSWR,675)
  675 FORMAT(/50H NEW MINIMUM FOUND.  GO BACK TO MINIMIZATION STEP./1H ,
     +60(1H=)/60X,1HV/60X,1HV/60X,1HV/57X,7HVVVVVVV/58X,5HVVVVV/59X,
     +3HVVV/60X,1HV//)
      GO TO 900
  700 WRITE (ISYSWR,'(A)') ' THERE ARE NO MINOS ERRORS TO CALCULATE.'
  900 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMNOT(FCN,ILAX,ILAX2,VAL2PL,VAL2MI,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Performs a MINOS error analysis on one parameter.
CC        The parameter ILAX is varied, and the minimum of the
CC        function with respect to the other parameters is followed
CC        until it crosses the value FMIN+UP.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION XDEV(MNI),W(MNI),GCC(MNI)
      CHARACTER*4 CPOS,CNEG,CSIG
      PARAMETER (CPOS='POSI',CNEG='NEGA')
C                                        . . save and prepare start vals
      ISW2 = ISW(2)
      ISW4 = ISW(4)
      SIGSAV = EDM
      ISTRAV = ISTRAT
      DC = DCOVAR
      LNEWMN = .FALSE.
      APSI  = EPSI*0.5
      ABEST=AMIN
      MPAR=NPAR
      NFMXIN = NFCNMX
      DO 125 I= 1, MPAR
  125 XT(I) = X(I)
      DO 130 J= 1, MPAR*(MPAR+1)/2
  130 VTHMAT(J) = VHMAT(J)
      DO 135 I= 1, MPAR
      GCC(I) = GLOBCC(I)
  135 W(I) = WERR(I)
      IT = NIOFEX(ILAX)
      ERP(IT) = 0.
      ERN(IT) = 0.
      CALL MNINEX(XT)
      UT = U(ILAX)
      IF (NVARL(ILAX) .EQ. 1) THEN
         ALIM(ILAX) = UT -100.*W(IT)
         BLIM(ILAX) = UT +100.*W(IT)
         ENDIF
      NDEX = IT*(IT+1)/2
      XUNIT = SQRT(UP/VTHMAT(NDEX))
      MARC = 0
      DO 162 I= 1, MPAR
      IF (I .EQ. IT)  GO TO 162
      MARC = MARC + 1
         IMAX = MAX(IT,I)
         INDX = IMAX*(IMAX-1)/2 + MIN(IT,I)
      XDEV(MARC) = XUNIT*VTHMAT(INDX)
  162 CONTINUE
C                           fix the parameter in question
      CALL MNFIXP (IT,IERR)
      IF (IERR .GT. 0)  THEN
         WRITE (ISYSWR,'(A,I5,A,I5)')
     +    ' MINUIT ERROR. CANNOT FIX PARAMETER',ILAX,'    INTERNAL',IT
         GO TO 700
      ENDIF
C                       . . . . . Nota Bene: from here on, NPAR=MPAR-1
C      Remember: MNFIXP squeezes IT out of X, XT, WERR, and VHMAT,
C                                                    not W, VTHMAT
      DO 500 ISIG= 1,2
      IF (ISIG .EQ. 1) THEN
         SIG = 1.0
         CSIG = CPOS
      ELSE
         SIG = -1.0
         CSIG = CNEG
      ENDIF
C                                        . sig=sign of error being calcd
      IF (ISW(5) .GT. 1) WRITE (ISYSWR,806)  CSIG,ILAX,CPNAM(ILAX)
  806 FORMAT (/' DETERMINATION OF ',A4,'TIVE MINOS ERROR FOR PARAMETER',
     +    I3, 2X ,A)
      IF (ISW(2).LE.0) CALL MNWARN('D','MINOS','NO COVARIANCE MATRIX.')
      NLIMIT = NFCN + NFMXIN
      ISTRAT = MAX(ISTRAV-1,0)
      DU1 = W(IT)
      U(ILAX) = UT + SIG*DU1
      U(ILAX) = MIN(U(ILAX),BLIM(ILAX))
      U(ILAX) = MAX(U(ILAX),ALIM(ILAX))
      DELU = U(ILAX) - UT
C         stop if already at limit with negligible step size
      IF (ABS(DELU)/(ABS(UT)+ABS(U(ILAX))) .LT. EPSMAC)  GO TO 440
      FAC = DELU/W(IT)
         DO 185 I= 1, NPAR
  185    X(I) = XT(I) + FAC*XDEV(I)
      IF (ISW(5) .GT. 1) WRITE (ISYSWR,801)  ILAX,UT,DELU,U(ILAX)
  801 FORMAT (/' PARAMETER',I4,' SET TO',E11.3,' + ',E10.3,' = ',E12.3)
C                                        loop to hit AMIN+UP
      KE1CR = ILAX
      KE2CR = 0
      XMIDCR = U(ILAX)
      XDIRCR = DELU
C
      AMIN = ABEST
      NFCNMX = NLIMIT - NFCN
      CALL MNCROS(FCN,AOPT,IERCR,FUTIL)
      IF (ABEST-AMIN .GT. 0.01*UP)  GO TO 650
      IF (IERCR .EQ. 1)  GO TO 440
      IF (IERCR .EQ. 2)  GO TO 450
      IF (IERCR .EQ. 3)  GO TO 460
C                                        . error successfully calculated
      EROS = XMIDCR-UT + AOPT*XDIRCR
      IF (ISW(5) .GT. 1) WRITE (ISYSWR,808)  CSIG,ILAX,CPNAM(ILAX),EROS
  808 FORMAT (/9X,4HTHE ,A4,  29HTIVE MINOS ERROR OF PARAMETER,I3,   2H
     +, ,A10,      4H, IS ,E12.4)
      GO TO 480
C                                        . . . . . . . . failure returns
  440 IF (ISW(5) .GE. 1) WRITE(ISYSWR,807)  CSIG,ILAX,CPNAM(ILAX)
  807 FORMAT (5X,'THE ',A4,'TIVE MINOS ERROR OF PARAMETER',I3,', ',A,
     +', EXCEEDS ITS LIMIT.'/)
      EROS = UNDEFI
      GO TO 480
  450 IF (ISW(5) .GE. 1) WRITE (ISYSWR, 802)  CSIG,ILAX,NFMXIN
  802 FORMAT (9X,'THE ',A,'TIVE MINOS ERROR',I4,' REQUIRES MORE THAN',
     +   I5,' FUNCTION CALLS.'/)
      EROS = 0.
      GO TO 480
  460 IF (ISW(5) .GE. 1) WRITE (ISYSWR, 805) CSIG,ILAX
  805 FORMAT (25X,A,'TIVE MINOS ERROR NOT CALCULATED FOR PARAMETER',I4/)
      EROS = 0.
C
  480 IF (ISW(5) .GT. 1) WRITE (ISYSWR,'(5X, 74(1H*))')
      IF (SIG .LT. ZERO)  THEN
         ERN(IT) = EROS
         IF (ILAX2.GT.0 .AND. ILAX2.LE.NU)  VAL2MI = U(ILAX2)
      ELSE
         ERP(IT) = EROS
         IF (ILAX2.GT.0 .AND. ILAX2.LE.NU)  VAL2PL = U(ILAX2)
      ENDIF
  500 CONTINUE
C                                        . . parameter finished. reset v
C                       normal termination
      ITAUR = 1
      CALL MNFREE(1)
      DO 550 J= 1, MPAR*(MPAR+1)/2
  550 VHMAT(J) = VTHMAT(J)
      DO 595 I= 1, MPAR
      WERR(I) = W(I)
      GLOBCC(I) = GCC(I)
  595 X(I) = XT(I)
      CALL MNINEX (X)
      EDM = SIGSAV
      AMIN = ABEST
      ISW(2) = ISW2
      ISW(4) = ISW4
      DCOVAR = DC
      GO TO 700
C                       new minimum
  650 LNEWMN = .TRUE.
      ISW(2) = 0
      DCOVAR = 1.
      ISW(4) = 0
      SAV = U(ILAX)
      ITAUR = 1
      CALL MNFREE(1)
      U(ILAX) = SAV
      CALL MNEXIN(X)
      EDM = BIGEDM
C                       in any case
  700 CONTINUE
      ITAUR = 0
      NFCNMX = NFMXIN
      ISTRAT = ISTRAV
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:50  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPARM(K,CNAMJ,UK,WK,A,B,IERFLG)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MNPARS and user-callable
CC    Implements one parameter definition, that is:
CC          K     (external) parameter number
CC          CNAMK parameter name
CC          UK    starting value
CC          WK    starting step size or uncertainty
CC          A, B  lower and upper physical parameter limits
CC    and sets up (updates) the parameter lists.
CC    Output: IERFLG=0 if no problems
CC                  >0 if MNPARM unable to implement definition
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER*(*) CNAMJ
      CHARACTER  CNAMK*10, CHBUFI*4
C
      CNAMK = CNAMJ
      KINT = NPAR
      IF (K.LT.1 .OR. K.GT.MAXEXT) THEN
C                     parameter number exceeds allowed maximum value
        WRITE (ISYSWR,9)  K,MAXEXT
    9   FORMAT (/' MINUIT USER ERROR.  PARAMETER NUMBER IS',I11/
     +         ',  ALLOWED RANGE IS ONE TO',I4/)
        GO TO 800
      ENDIF
C                     normal parameter request
      KTOFIX = 0
      IF (NVARL(K) .LT. 0) GO TO 50
C         previously defined parameter is being redefined
C                                     find if parameter was fixed
      DO 40 IX= 1, NPFIX
      IF (IPFIX(IX) .EQ. K)  KTOFIX = K
   40 CONTINUE
      IF (KTOFIX .GT. 0)  THEN
         CALL MNWARN('W','PARAM DEF','REDEFINING A FIXED PARAMETER.')
         IF (KINT .GE. MAXINT)  THEN
            WRITE (ISYSWR,'(A)') ' CANNOT RELEASE. MAX NPAR EXCEEDED.'
            GO TO 800
            ENDIF
         CALL MNFREE(-K)
         ENDIF
C                       if redefining previously variable parameter
      IF(NIOFEX(K) .GT. 0) KINT = NPAR-1
   50 CONTINUE
C
C                                      . . .print heading
      IF (LPHEAD .AND. ISW(5).GE.0) THEN
        WRITE (ISYSWR,61)
        LPHEAD = .FALSE.
      ENDIF
   61 FORMAT(/' PARAMETER DEFINITIONS:'/
     +        '    NO.   NAME         VALUE      STEP SIZE      LIMITS')
      IF (WK .GT. ZERO)  GO TO 122
C                                        . . .constant parameter . . . .
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 82)  K,CNAMK,UK
   82 FORMAT (1X,I5,1X,1H',A10,1H',1X,G13.5, '  constant')
      NVL = 0
      GO TO 200
  122 IF (A.EQ.ZERO .AND. B.EQ.ZERO) THEN
C                                      variable parameter without limits
      NVL = 1
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 127)  K,CNAMK,UK,WK
  127 FORMAT (1X,I5,1X,1H',A10,1H',1X,2G13.5, '     no limits')
      ELSE
C                                         variable parameter with limits
      NVL = 4
      LNOLIM = .FALSE.
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 132)  K,CNAMK,UK,WK,A,B
  132 FORMAT(1X,I5,1X,1H',A10,1H',1X,2G13.5,2X,2G13.5)
      ENDIF
C                             . . request for another variable parameter
      KINT = KINT + 1
      IF (KINT .GT. MAXINT)  THEN
         WRITE (ISYSWR,135)  MAXINT
  135    FORMAT (/' MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.'/
     +   ' THIS VERSION OF MINUIT DIMENSIONED FOR',I4//)
         GO TO 800
         ENDIF
      IF (NVL .EQ. 1)  GO TO 200
      IF (A .EQ. B)  THEN
        WRITE (ISYSWR,'(/A,A/A/)') ' USER ERROR IN MINUIT PARAMETER',
     +   ' DEFINITION',' UPPER AND LOWER LIMITS EQUAL.'
        GO TO 800
        ENDIF
      IF (B .LT. A) THEN
         SAV = B
         B = A
         A = SAV
         CALL MNWARN('W','PARAM DEF','PARAMETER LIMITS WERE REVERSED.')
         IF (LWARN) LPHEAD=.TRUE.
         ENDIF
      IF ((B-A) .GT. 1.0E7)  THEN
         WRITE (CHBUFI,'(I4)') K
         CALL MNWARN('W','PARAM DEF',
     +               'LIMITS ON PARAM'//CHBUFI//' TOO FAR APART.')
         IF (LWARN) LPHEAD=.TRUE.
      ENDIF
      DANGER = (B-UK)*(UK-A)
      IF (DANGER .LT. 0.)
     +     CALL MNWARN('W','PARAM DEF','STARTING VALUE OUTSIDE LIMITS.')
      IF (DANGER .EQ. 0.)
     +     CALL MNWARN('W','PARAM DEF','STARTING VALUE IS AT LIMIT.')
  200 CONTINUE
C                           . . . input OK, set values, arrange lists,
C                                    calculate step sizes GSTEP, DIRIN
      CFROM = 'PARAMETR'
      NFCNFR = NFCN
      CSTATU= 'NEW VALUES'
      NU = MAX(NU,K)
      CPNAM(K) = CNAMK
      U(K) = UK
      ALIM(K) = A
      BLIM(K) = B
      NVARL(K) = NVL
C                             K is external number of new parameter
C           LASTIN is the number of var. params with ext. param. no.< K
      LASTIN = 0
      DO 240 IX= 1, K-1
      IF (NIOFEX(IX) .GT. 0)  LASTIN=LASTIN+1
  240 CONTINUE
C                 KINT is new number of variable params, NPAR is old
      IF (KINT .EQ. NPAR)  GO TO 280
      IF (KINT .GT. NPAR) THEN
C                          insert new variable parameter in list
         DO 260 IN= NPAR,LASTIN+1,-1
         IX = NEXOFI(IN)
         NIOFEX(IX) = IN+1
         NEXOFI(IN+1)= IX
         X    (IN+1) = X    (IN)
         XT   (IN+1) = XT   (IN)
         DIRIN(IN+1) = DIRIN(IN)
         G2   (IN+1) = G2   (IN)
         GSTEP(IN+1) = GSTEP(IN)
  260    CONTINUE
      ELSE
C                          remove variable parameter from list
         DO 270 IN= LASTIN+1,KINT
         IX = NEXOFI(IN+1)
         NIOFEX(IX) = IN
         NEXOFI(IN)= IX
         X     (IN)= X    (IN+1)
         XT    (IN)= XT   (IN+1)
         DIRIN (IN)= DIRIN(IN+1)
         G2    (IN)= G2   (IN+1)
         GSTEP (IN)= GSTEP(IN+1)
  270    CONTINUE
      ENDIF
  280 CONTINUE
      IX = K
      NIOFEX(IX) = 0
      NPAR = KINT
      CALL MNRSET(1)
C                                       lists are now arranged . . . .
      IF (NVL .GT. 0)  THEN
         IN = LASTIN+1
         NEXOFI(IN) = IX
         NIOFEX(IX) = IN
         SAV = U(IX)
         CALL MNPINT(SAV,IX,PINTI)
         X(IN) = PINTI
         XT(IN) = X(IN)
         WERR(IN) = WK
         SAV2 = SAV + WK
         CALL MNPINT(SAV2,IX,PINTI)
         VPLU = PINTI - X(IN)
         SAV2 = SAV - WK
         CALL MNPINT(SAV2,IX,PINTI)
         VMINU = PINTI - X(IN)
         DIRIN(IN) = 0.5 * (ABS(VPLU) +ABS(VMINU))
         G2(IN) = 2.0*UP / DIRIN(IN)**2
         GSMIN = 8.*EPSMA2*ABS(X(IN))
         GSTEP(IN) = MAX (GSMIN, 0.1*DIRIN(IN))
         IF (AMIN .NE. UNDEFI) THEN
             SMALL = SQRT(EPSMA2*(AMIN+UP)/UP)
             GSTEP(IN) = MAX(GSMIN, SMALL*DIRIN(IN))
         ENDIF
         GRD  (IN) = G2(IN)*DIRIN(IN)
C                   if parameter has limits
         IF (NVARL(K) .GT. 1) THEN
            IF (GSTEP(IN).GT. 0.5)  GSTEP(IN)=0.5
            GSTEP(IN) = -GSTEP(IN)
         ENDIF
      ENDIF
      IF (KTOFIX .GT. 0)  THEN
         KINFIX = NIOFEX(KTOFIX)
         IF (KINFIX .GT. 0)  CALL MNFIXP(KINFIX,IERR)
         IF (IERR .GT. 0)  GO TO 800
      ENDIF
      IERFLG = 0
      RETURN
C                   error on input, unable to implement request  . . . .
  800 CONTINUE
      IERFLG = 1
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPARS(CRDBUF,ICONDN)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MNREAD and user-callable
CC    Implements one parameter definition, that is:
CC       parses the string CRDBUF and calls MNPARM
C
C output conditions:
C        ICONDN = 0    all OK
C        ICONDN = 1    error, attempt to define parameter is ignored
C        ICONDN = 2    end of parameter definitions
C
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      DIMENSION PLIST(MAXP)
      CHARACTER CNAMK*10, CRDBUF*(*) , CELMNT*20 , COMAND*(MAXCWD)
C
      LENBUF = LEN(CRDBUF)
C                     find out whether fixed or free-field format
      KAPO1 = INDEX(CRDBUF,'''')
      IF (KAPO1 .EQ. 0)  GO TO 150
      KAPO2 = INDEX(CRDBUF(KAPO1+1:),'''')
      IF (KAPO2 .EQ. 0)  GO TO 150
C          new (free-field) format
      KAPO2 = KAPO2 + KAPO1
C                             skip leading blanks if any
         DO 115 ISTART=1, KAPO1-1
         IF (CRDBUF(ISTART:ISTART) .NE. ' ')  GO TO 120
  115    CONTINUE
         GO TO 210
  120 CONTINUE
C                               parameter number integer
      CELMNT = CRDBUF(ISTART:KAPO1-1)
      READ (CELMNT,'(BN,F20.0)',ERR=180) FK
      K = FK
      IF (K .LE. 0)  GO TO 210
      CNAMK = 'PARAM '//CELMNT
      IF (KAPO2-KAPO1 .GT. 1) CNAMK = CRDBUF(KAPO1+1:KAPO2-1)
C  special handling if comma or blanks and a comma follow 'name'
        DO 135 ICY= KAPO2+1,LENBUF
        IF (CRDBUF(ICY:ICY) .EQ. ',') GO TO 139
        IF (CRDBUF(ICY:ICY) .NE. ' ') GO TO 140
  135 CONTINUE
        UK = 0.
        WK = 0.
        A  = 0.
        B = 0.
      GO TO 170
  139 CONTINUE
      ICY = ICY+1
  140 CONTINUE
      IBEGIN = ICY
      CALL MNCRCK(CRDBUF(IBEGIN:),MAXCWD,COMAND,LNC,
     +                             MAXP,PLIST,LLIST, IERR,ISYSWR)
      IF (IERR .GT. 0)  GO TO 180
      UK = PLIST(1)
      WK = 0.
      IF (LLIST .GE. 2)  WK = PLIST(2)
      A = 0.
      IF (LLIST .GE. 3)  A = PLIST(3)
      B = 0.
      IF (LLIST .GE. 4)  B = PLIST(4)
      GO TO 170
C          old (fixed-field) format
  150 CONTINUE
      READ (CRDBUF, 158,ERR=180)  XK,CNAMK,UK,WK,A,B
  158 FORMAT (BN,F10.0, A10, 4F10.0)
      K = XK
      IF (K .EQ. 0)  GO TO 210
C          parameter format cracked, implement parameter definition
  170 CALL MNPARM(K,CNAMK,UK,WK,A,B,IERR)
      ICONDN = IERR
      RETURN
C          format or other error
  180 CONTINUE
      ICONDN = 1
      RETURN
C        end of data
  210 CONTINUE
      ICONDN = 2
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPFIT(PARX2P,PARY2P,NPAR2P,COEF2P,SDEV2P)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     to fit a parabola to npar2p points
C
C   npar2p   no. of points
C   parx2p(i)   x value of point i
C   pary2p(i)   y value of point i
C
C   coef2p(1...3)  coefficients of the fitted parabola
C   y=coef2p(1) + coef2p(2)*x + coef2p(3)*x**2
C   sdev2p= variance
C   method : chi**2 = min equation solved explicitly
      DIMENSION PARX2P(NPAR2P),PARY2P(NPAR2P),COEF2P(NPAR2P)
      DIMENSION CZ(3)
C
      DO 3  I=1,3
    3 CZ(I)=0.
      SDEV2P=0.
      IF(NPAR2P.LT.3) GO TO 10
      F=NPAR2P
C--- center x values for reasons of machine precision
      XM=0.
      DO 2  I=1,NPAR2P
    2 XM=XM+PARX2P(I)
      XM=XM/F
      X2=0.
      X3=0.
      X4=0.
      Y=0.
      Y2=0.
      XY=0.
      X2Y=0.
      DO 1  I=1,NPAR2P
      S=PARX2P(I)-XM
      T=PARY2P(I)
      S2=S*S
      X2=X2+S2
      X3=X3+S*S2
      X4=X4+S2*S2
      Y=Y+T
      Y2=Y2+T*T
      XY=XY+S*T
      X2Y=X2Y+S2*T
    1 CONTINUE
      A=(F*X4-X2**2)*X2-F*X3**2
      IF(A.EQ.0.)  GOTO 10
      CZ(3)=(X2*(F*X2Y-X2*Y)-F*X3*XY)/A
      CZ(2)=(XY-X3*CZ(3))/X2
      CZ(1)=(Y-X2*CZ(3))/F
      IF(NPAR2P.EQ.3)  GOTO 6
      SDEV2P=Y2-(CZ(1)*Y+CZ(2)*XY+CZ(3)*X2Y)
      IF(SDEV2P.LT.0.)  SDEV2P=0.
      SDEV2P=SDEV2P/(F-3.)
    6 CZ(1)=CZ(1)+XM*(XM*CZ(3)-CZ(2))
      CZ(2)=CZ(2)-2.*XM*CZ(3)
   10 CONTINUE
      DO 11  I=1,3
   11 COEF2P(I)=CZ(I)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPINT(PEXTI,I,PINTI)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Calculates the internal parameter value PINTI corresponding
CC        to the external value PEXTI for parameter I.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER CHBUFI*4, CHBUF2*30
      PINTI = PEXTI
      IGO = NVARL(I)
      IF (IGO .EQ. 4)  THEN
C--                          there are two limits
        ALIMI = ALIM(I)
        BLIMI = BLIM(I)
        YY=2.0*(PEXTI-ALIMI)/(BLIMI-ALIMI) - 1.0
        YY2 = YY**2
        IF (YY2 .GE. (1.0- EPSMA2))  THEN
           IF (YY .LT. 0.) THEN
               A = VLIMLO
               CHBUF2 = ' IS AT ITS LOWER ALLOWED LIMIT.'
           ELSE
               A = VLIMHI
               CHBUF2 = ' IS AT ITS UPPER ALLOWED LIMIT.'
           ENDIF
           PINTI = A
           PEXTI = ALIMI + 0.5* (BLIMI-ALIMI) *(SIN(A) +1.0)
           LIMSET = .TRUE.
           WRITE (CHBUFI,'(I4)') I
           IF (YY2 .GT. 1.0) CHBUF2 = ' BROUGHT BACK INSIDE LIMITS.'
           CALL MNWARN('W',CFROM,'VARIABLE'//CHBUFI//CHBUF2)
         ELSE
           PINTI = ASIN(YY)
         ENDIF
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPLOT(XPT,YPT,CHPT,NXYPT,NUNIT,NPAGWD,NPAGLN)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        plots points in array xypt onto one page with labelled axes
CC        NXYPT is the number of points to be plotted
CC        XPT(I) = x-coord. of ith point
CC        YPT(I) = y-coord. of ith point
CC        CHPT(I) = character to be plotted at this position
CC        the input point arrays XPT, YPT, CHPT are destroyed.
CC
      DIMENSION   XPT(*), YPT(*)
      CHARACTER*1 CHPT(*) ,  CHSAV,  CHBEST, CDOT, CSLASH, CBLANK
      PARAMETER (MAXWID=100)
      CHARACTER CLINE*100, CHMESS*30
      DIMENSION XVALUS(12)
      LOGICAL OVERPR
      DATA CDOT,CSLASH,CBLANK/ '.' , '/' , ' '/
      MAXNX = MIN(NPAGWD-20,MAXWID)
      IF (MAXNX .LT. 10)  MAXNX = 10
      MAXNY = NPAGLN
      IF (MAXNY .LT. 10)  MAXNY = 10
      IF (NXYPT .LE. 1)  RETURN
      XBEST = XPT(1)
      YBEST = YPT(1)
      CHBEST = CHPT(1)
C         order the points by decreasing y
      KM1 = NXYPT - 1
      DO 150 I= 1, KM1
      IQUIT = 0
      NI = NXYPT - I
      DO 140 J= 1, NI
      IF (YPT(J) .GT. YPT(J+1)) GO TO 140
        SAVX = XPT(J)
        XPT(J) = XPT(J+1)
        XPT(J+1) = SAVX
        SAVY = YPT(J)
        YPT(J) = YPT(J+1)
        YPT(J+1) = SAVY
        CHSAV = CHPT(J)
        CHPT(J) = CHPT(J+1)
        CHPT(J+1) = CHSAV
      IQUIT = 1
  140 CONTINUE
      IF (IQUIT .EQ. 0) GO TO 160
  150 CONTINUE
  160 CONTINUE
C         find extreme values
      XMAX = XPT(1)
      XMIN = XMAX
      DO 200 I= 1, NXYPT
        IF (XPT(I) .GT. XMAX)  XMAX = XPT(I)
        IF (XPT(I) .LT. XMIN)  XMIN = XPT(I)
  200 CONTINUE
      DXX = 0.001*(XMAX-XMIN)
      XMAX = XMAX + DXX
      XMIN = XMIN - DXX
      CALL MNBINS(XMIN,XMAX,MAXNX,XMIN,XMAX,NX,BWIDX)
      YMAX = YPT(1)
      YMIN = YPT(NXYPT)
      IF (YMAX .EQ. YMIN)  YMAX=YMIN+1.0
      DYY = 0.001*(YMAX-YMIN)
      YMAX = YMAX + DYY
      YMIN = YMIN - DYY
      CALL MNBINS(YMIN,YMAX,MAXNY,YMIN,YMAX,NY,BWIDY)
      ANY = NY
C         if first point is blank, it is an 'origin'
      IF (CHBEST .EQ. CBLANK)  GO TO 50
      XBEST = 0.5 * (XMAX+XMIN)
      YBEST = 0.5 * (YMAX+YMIN)
   50 CONTINUE
C         find scale constants
      AX = 1.0/BWIDX
      AY = 1.0/BWIDY
      BX = -AX*XMIN + 2.0
      BY = -AY*YMIN - 2.0
C         convert points to grid positions
      DO 300 I= 1, NXYPT
      XPT(I) = AX*XPT(I) + BX
  300 YPT(I) = ANY-AY*YPT(I) - BY
      NXBEST = AX*XBEST + BX
      NYBEST = ANY  - AY*YBEST - BY
C         print the points
      NY = NY + 2
      NX = NX + 2
      ISP1 = 1
      LINODD = 1
      OVERPR=.FALSE.
      DO 400 I= 1, NY
      DO 310 IBK= 1, NX
  310 CLINE (IBK:IBK) = CBLANK
      CLINE(1:1) = CDOT
      CLINE(NX:NX) = CDOT
      CLINE(NXBEST:NXBEST) = CDOT
      IF (I.NE.1 .AND. I.NE.NYBEST .AND. I.NE.NY)  GO TO 320
      DO 315 J= 1, NX
  315 CLINE(J:J) = CDOT
  320 CONTINUE
      YPRT = YMAX - FLOAT(I-1)*BWIDY
      IF (ISP1 .GT. NXYPT)  GO TO 350
C         find the points to be plotted on this line
        DO 341 K= ISP1,NXYPT
      KS = YPT(K)
      IF (KS .GT. I)  GO TO 345
      IX = XPT(K)
      IF (CLINE(IX:IX) .EQ.   CDOT)  GO TO 340
      IF (CLINE(IX:IX) .EQ. CBLANK)  GO TO 340
      IF (CLINE(IX:IX) .EQ.CHPT(K))  GO TO 341
      OVERPR = .TRUE.
C         OVERPR is true if one or more positions contains more than
C            one point
      CLINE(IX:IX) = '&'
      GO TO 341
  340 CLINE(IX:IX) = CHPT(K)
  341 CONTINUE
        ISP1 = NXYPT + 1
        GO TO 350
  345   ISP1 = K
  350 CONTINUE
      IF (LINODD .EQ. 1 .OR. I .EQ. NY)  GO TO 380
      LINODD = 1
      WRITE (NUNIT, '(18X,A)')       CLINE(:NX)
      GO TO 400
  380 WRITE (NUNIT,'(1X,G14.7,A,A)') YPRT, ' ..', CLINE(:NX)
      LINODD = 0
  400 CONTINUE
C         print labels on x-axis every ten columns
      DO 410 IBK= 1, NX
      CLINE(IBK:IBK) = CBLANK
      IF (MOD(IBK,10) .EQ. 1)  CLINE(IBK:IBK) = CSLASH
  410 CONTINUE
      WRITE (NUNIT, '(18X,A)')       CLINE(:NX)
C
      DO 430 IBK= 1, 12
  430 XVALUS(IBK) = XMIN + FLOAT(IBK-1)*10.*BWIDX
      ITEN = (NX+9) / 10
      WRITE (NUNIT,'(12X,12G10.4)')  (XVALUS(IBK), IBK=1,ITEN)
      CHMESS = ' '
      IF (OVERPR) CHMESS='   Overprint character is &'
      WRITE (NUNIT,'(25X,A,G13.7,A)') 'ONE COLUMN=',BWIDX, CHMESS
  500 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPOUT(IUEXT,CHNAM,VAL,ERR,XLOLIM,XUPLIM,IUINT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC     User-called
CC   Provides the user with information concerning the current status
CC          of parameter number IUEXT. Namely, it returns:
CC        CHNAM: the name of the parameter
CC        VAL: the current (external) value of the parameter
CC        ERR: the current estimate of the parameter uncertainty
CC        XLOLIM: the lower bound (or zero if no limits)
CC        XUPLIM: the upper bound (or zero if no limits)
CC        IUINT: the internal parameter number (or zero if not variable,
CC           or negative if undefined).
CC  Note also:  If IUEXT is negative, then it is -internal parameter
CC           number, and IUINT is returned as the EXTERNAL number.
CC     Except for IUINT, this is exactly the inverse of MNPARM
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER*(*) CHNAM
      XLOLIM = 0.
      XUPLIM = 0.
      ERR = 0.
      IF (IUEXT .EQ. 0)  GO TO 100
      IF (IUEXT .LT. 0)  THEN
C                   internal parameter number specified
         IINT = -IUEXT
         IF (IINT .GT. NPAR) GO TO 100
         IEXT = NEXOFI(IINT)
         IUINT = IEXT
      ELSE
C                    external parameter number specified
         IEXT = IUEXT
         IF (IEXT .EQ. 0)   GO TO 100
         IF (IEXT .GT. NU)  GO TO 100
         IINT = NIOFEX(IEXT)
         IUINT = IINT
      ENDIF
C                     in both cases
         NVL = NVARL(IEXT)
         IF (NVL .LT. 0) GO TO 100
      CHNAM = CPNAM(IEXT)
      VAL = U(IEXT)
      IF (IINT .GT. 0)  ERR = WERR(IINT)
      IF (NVL .EQ. 4) THEN
         XLOLIM = ALIM(IEXT)
         XUPLIM = BLIM(IEXT)
      ENDIF
      RETURN
C                parameter is undefined
  100 IUINT = -1
      CHNAM = 'undefined'
      VAL = 0.
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPRIN  (INKODE,FVAL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Prints the values of the parameters at the time of the call.
CC        also prints other relevant information such as function value,
CC        estimated distance to minimum, parameter errors, step sizes.
CC
C         According to the value of IKODE, the printout is:
C    IKODE=INKODE= 0    only info about function value
C                  1    parameter values, errors, limits
C                  2    values, errors, step sizes, internal values
C                  3    values, errors, step sizes, first derivs.
C                  4    values, parabolic errors, MINOS errors
C    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2)
C
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      CHARACTER*14 COLHDU(6),COLHDL(6), CX2,CX3,CGETX
      CHARACTER*11 CNAMBF, CBLANK
      CHARACTER  CHEDM*10, CHEVAL*15
      PARAMETER (CGETX='PLEASE GET X..')
      DATA CBLANK/'          '/
C
      IF (NU .EQ. 0)  THEN
       WRITE (ISYSWR,'(A)') ' THERE ARE CURRENTLY NO PARAMETERS DEFINED'
       GO TO 700
      ENDIF
C                  get value of IKODE based in INKODE, ISW(2)
      IKODE = INKODE
      IF (INKODE .EQ. 5) THEN
         IKODE = ISW(2)+1
         IF (IKODE .GT. 3)  IKODE=3
      ENDIF
C                  set 'default' column headings
      DO 5 K= 1, 6
      COLHDU(K) = 'UNDEFINED'
    5 COLHDL(K) = 'COLUMN HEAD'
C              print title if Minos errors, and title exists.
      IF (IKODE.EQ.4 .AND. CTITL.NE.CUNDEF)
     +            WRITE (ISYSWR,'(/A,A)')  ' MINUIT TASK: ',CTITL
C              report function value and status
      IF (FVAL .EQ. UNDEFI) THEN
         CHEVAL = ' unknown       '
      ELSE
         WRITE (CHEVAL,'(G15.7)') FVAL
      ENDIF
         IF (EDM .EQ. BIGEDM) THEN
            CHEDM = ' unknown  '
         ELSE
            WRITE (CHEDM, '(E10.2)') EDM
         ENDIF
      NC = NFCN-NFCNFR
      WRITE (ISYSWR,905)  CHEVAL,CFROM,CSTATU,NC,NFCN
  905 FORMAT (/' FCN=',A,' FROM ',A8,'  STATUS=',A10,I6,' CALLS',
     +         I9,' TOTAL')
      M = ISW(2)
      IF (M.EQ.0 .OR. M.EQ.2 .OR. DCOVAR.EQ.ZERO) THEN
        WRITE (ISYSWR,907) CHEDM,ISTRAT,COVMES(M)
  907   FORMAT (21X,'EDM=',A,'    STRATEGY=',I2,6X,A)
      ELSE
        DCMAX = 1.
        DC = MIN(DCOVAR,DCMAX) * 100.
        WRITE (ISYSWR,908) CHEDM,ISTRAT,DC
  908   FORMAT (21X,'EDM=',A,'  STRATEGY=',I1,'  ERROR MATRIX',
     +     ' UNCERTAINTY=',F5.1,'%')
      ENDIF
C
      IF (IKODE .EQ. 0)  GO TO 700
C               find longest name (for Rene!)
      NTRAIL = 10
      DO 20 I= 1, NU
         IF (NVARL(I) .LT. 0)  GO TO 20
         DO 15 IC= 10,1,-1
            IF (CPNAM(I)(IC:IC) .NE. ' ') GO TO 16
   15    CONTINUE
         IC = 1
   16    LBL = 10-IC
         IF (LBL .LT. NTRAIL)  NTRAIL=LBL
   20 CONTINUE
      NADD = NTRAIL/2 + 1
      IF (IKODE .EQ. 1)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '      PHYSICAL'
         COLHDU(3) = ' LIMITS       '
         COLHDL(2) = '    NEGATIVE  '
         COLHDL(3) = '    POSITIVE  '
      ENDIF
      IF (IKODE .EQ. 2)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '    INTERNAL  '
         COLHDL(2) = '    STEP SIZE '
         COLHDU(3) = '    INTERNAL  '
         COLHDL(3) = '      VALUE   '
      ENDIF
      IF (IKODE .EQ. 3)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '       STEP   '
         COLHDL(2) = '       SIZE   '
         COLHDU(3) = '      FIRST   '
         COLHDL(3) = '   DERIVATIVE '
      ENDIF
      IF (IKODE .EQ. 4)  THEN
         COLHDU(1) = '    PARABOLIC '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '        MINOS '
         COLHDU(3) = 'ERRORS        '
         COLHDL(2) = '   NEGATIVE   '
         COLHDL(3) = '   POSITIVE   '
      ENDIF
C
      IF (IKODE .NE. 4)  THEN
         IF (ISW(2) .LT. 3) COLHDU(1)='  APPROXIMATE '
         IF (ISW(2) .LT. 1) COLHDU(1)=' CURRENT GUESS'
      ENDIF
      NCOL = 3
      WRITE (ISYSWR, 910) (COLHDU(KK),KK=1,NCOL)
      WRITE (ISYSWR, 911) (COLHDL(KK),KK=1,NCOL)
  910 FORMAT (/'  EXT PARAMETER ',     13X       ,6A14)
  911 FORMAT ( '  NO.   NAME    ','    VALUE    ',6A14)
C
C                                        . . . loop over parameters . .
      DO 200 I= 1, NU
      IF (NVARL(I) .LT. 0)  GO TO 200
      L = NIOFEX(I)
      CNAMBF = CBLANK(1:NADD)//CPNAM(I)
      IF (L .EQ. 0)  GO TO 55
C              variable parameter.
      X1 = WERR(L)
      CX2 = CGETX
      CX3 = CGETX
      IF (IKODE .EQ. 1) THEN
         IF (NVARL(I) .LE. 1) THEN
            WRITE (ISYSWR, 952)  I,CNAMBF,U(I),X1
            GO TO 200
         ELSE
         X2 = ALIM(I)
         X3 = BLIM(I)
         ENDIF
      ENDIF
      IF (IKODE .EQ. 2) THEN
         X2 = DIRIN(L)
         X3 = X(L)
      ENDIF
      IF (IKODE .EQ. 3) THEN
         X2 = DIRIN(L)
         X3 = GRD(L)
         IF (NVARL(I).GT.1 .AND. ABS(COS(X(L))) .LT. 0.001)
     +      CX3 = '** at limit **'
      ENDIF
      IF (IKODE .EQ. 4) THEN
         X2 = ERN(L)
           IF (X2.EQ.ZERO)   CX2=' '
           IF (X2.EQ.UNDEFI) CX2='   at limit   '
         X3 = ERP(L)
           IF (X3.EQ.ZERO)   CX3=' '
           IF (X3.EQ.UNDEFI) CX3='   at limit   '
      ENDIF
      IF (CX2.EQ.CGETX) WRITE (CX2,'(G14.5)') X2
      IF (CX3.EQ.CGETX) WRITE (CX3,'(G14.5)') X3
      WRITE (ISYSWR,952)   I,CNAMBF,U(I),X1,CX2,CX3
  952 FORMAT (I4,1X,A11,2G14.5,2A)
C               check if parameter is at limit
      IF (NVARL(I) .LE. 1 .OR. IKODE .EQ. 3)  GO TO 200
      IF (ABS(COS(X(L))) .LT. 0.001)  WRITE (ISYSWR,1004)
 1004 FORMAT (1H ,32X,42HWARNING -   - ABOVE PARAMETER IS AT LIMIT.)
      GO TO 200
C
C                                print constant or fixed parameter.
   55 CONTINUE
                          COLHDU(1) = '   constant   '
      IF (NVARL(I).GT.0)  COLHDU(1) = '     fixed    '
      IF (NVARL(I).EQ.4 .AND. IKODE.EQ.1) THEN
        WRITE (ISYSWR,'(I4,1X,A11,G14.5,A,2G14.5)')
     +     I,CNAMBF,U(I),COLHDU(1),ALIM(I),BLIM(I)
      ELSE
        WRITE (ISYSWR,'(I4,1X,A11,G14.5,A)')  I,CNAMBF,U(I),COLHDU(1)
      ENDIF
  200 CONTINUE
C
      IF (UP.NE.UPDFLT)  WRITE (ISYSWR,'(31X,A,G10.3)') 'ERR DEF=',UP
  700 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:50  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPSDF
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        calculates the eigenvalues of v to see if positive-def.
CC        if not, adds constant along diagonal to make positive.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER CHBUFF*12
      DIMENSION S(MNI)
      EPSMIN = 1.0E-6
      EPSPDF = MAX(EPSMIN, EPSMA2)
      DGMIN = VHMAT(1)
C                        Check if negative or zero on diagonal
      DO 200 I= 1, NPAR
      NDEX = I*(I+1)/2
      IF (VHMAT(NDEX) .LE. ZERO) THEN
          WRITE (CHBUFF(1:3),'(I3)') I
          CALL MNWARN('W',CFROM,
     +'Negative diagonal element'//CHBUFF(1:3)//' in Error Matrix')
      ENDIF
      IF (VHMAT(NDEX) .LT. DGMIN)  DGMIN = VHMAT(NDEX)
  200 CONTINUE
      IF (DGMIN .LE. ZERO) THEN
         DG = (ONE+EPSPDF) - DGMIN
         WRITE (CHBUFF,'(E12.2)') DG
         CALL MNWARN('W',CFROM,
     +     CHBUFF//' added to diagonal of error matrix')
      ELSE
         DG = ZERO
      ENDIF
C                    Store VHMAT in P, make sure diagonal pos.
      DO 213 I= 1, NPAR
      NDEX = I*(I-1)/2
      NDEXD = NDEX + I
      VHMAT(NDEXD) = VHMAT(NDEXD) + DG
      IF (VHMAT(NDEXD) .LE. ZERO)   VHMAT(NDEXD) = 1.0
      S(I) = 1.0/SQRT(VHMAT(NDEXD))
      DO 213 J= 1, I
      NDEX =  NDEX + 1
  213 P(I,J) = VHMAT(NDEX) * S(I)*S(J)
C      call eigen (p,p,maxint,npar,pstar,-npar)
      CALL MNEIG(P,MAXINT,NPAR,MAXINT,PSTAR,EPSPDF,IFAULT)
      PMIN = PSTAR(1)
      PMAX = PSTAR(1)
      DO 215 IP= 2, NPAR
      IF (PSTAR(IP) .LT. PMIN)  PMIN = PSTAR(IP)
      IF (PSTAR(IP) .GT. PMAX)  PMAX = PSTAR(IP)
  215 CONTINUE
      PMAX = MAX(ABS(PMAX), ONE)
      IF ((PMIN .LE. ZERO .AND. LWARN) .OR.  ISW(5) .GE. 2) THEN
         WRITE (ISYSWR,550)
         WRITE (ISYSWR,551) (PSTAR(IP),IP=1,NPAR)
      ENDIF
      IF (PMIN .GT. EPSPDF*PMAX)  GO TO 217
      IF (ISW(2) .EQ. 3)  ISW(2)=2
      PADD = 1.0E-3*PMAX - PMIN
      DO 216 IP= 1, NPAR
      NDEX = IP*(IP+1)/2
  216 VHMAT(NDEX) = VHMAT(NDEX) *(1.0 + PADD)
      CSTATU= 'NOT POSDEF'
      WRITE (CHBUFF,'(G12.5)') PADD
      CALL MNWARN('W',CFROM,
     +   'MATRIX FORCED POS-DEF BY ADDING '//CHBUFF//' TO DIAGONAL.')
  217 CONTINUE
C
  550 FORMAT (' EIGENVALUES OF SECOND-DERIVATIVE MATRIX:' )
  551 FORMAT (7X,6E12.4)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNRAZZ(YNEW,PNEW,Y,JH,JL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called only by MNSIMP (and MNIMPR) to add a new point
CC        and remove an old one from the current simplex, and get the
CC        estimated distance to minimum.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION PNEW(*), Y(*)
      DO 10 I=1,NPAR
   10 P(I,JH) = PNEW(I)
      Y(JH)=YNEW
      IF(YNEW .LT. AMIN) THEN
        DO 15 I=1,NPAR
   15   X(I) = PNEW(I)
        CALL MNINEX(X)
        AMIN = YNEW
        CSTATU = 'PROGRESS  '
        JL=JH
      ENDIF
      JH = 1
      NPARP1 = NPAR+1
   20 DO 25 J=2,NPARP1
      IF (Y(J) .GT. Y(JH))  JH = J
   25 CONTINUE
      EDM = Y(JH) - Y(JL)
      IF (EDM .LE. ZERO)  GO TO 45
      DO 35 I= 1, NPAR
      PBIG = P(I,1)
      PLIT = PBIG
      DO 30 J= 2, NPARP1
      IF (P(I,J) .GT. PBIG)  PBIG = P(I,J)
      IF (P(I,J) .LT. PLIT)  PLIT = P(I,J)
   30 CONTINUE
      DIRIN(I) = PBIG - PLIT
   35 CONTINUE
   40 RETURN
   45 WRITE (ISYSWR, 1000)  NPAR
      GO TO 40
 1000 FORMAT ('   FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE',
     +    I3,' VARIABLE PARAMETERS.' /10X,'VERIFY THAT STEP SIZES ARE',
     +    ' BIG ENOUGH AND CHECK FCN LOGIC.'/1X,79(1H*)/1X,79(1H*)/)
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNREAD(FCN,IFLGIN,IFLGUT,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MINUIT.  Reads all user input to MINUIT.
CC     This routine is highly unstructured and defies normal logic.
CC
CC     IFLGIN indicates the function originally requested:
CC           = 1: read one-line title
CC             2: read parameter definitions
CC             3: read MINUIT commands
CC
CC     IFLGUT= 1: reading terminated normally
CC             2: end-of-data on input
CC             3: unrecoverable read error
CC             4: unable to process parameter requests
CC             5: more than 100 incomprehensible commands
CC internally,
CC     IFLGDO indicates the subfunction to be performed on the next
CC         input record: 1: read a one-line title
CC                       2: read a parameter definition
CC                       3: read a command
CC                       4: read in covariance matrix
CC     for example, when IFLGIN=3, but IFLGDO=1, then it should read
CC       a title, but this was requested by a command, not by MINUIT.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      CHARACTER  CRDBUF*80, CUPBUF*10
      CHARACTER CPROMT(3)*40, CLOWER*26, CUPPER*26
      LOGICAL LEOF
      DATA CPROMT/' ENTER MINUIT TITLE, or "SET INPUT n" : ',
     +            ' ENTER MINUIT PARAMETER DEFINITION:     ',
     +            ' ENTER MINUIT COMMAND:                  '/
C
      DATA CLOWER/'abcdefghijklmnopqrstuvwxyz'/
      DATA CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
      IFLGUT = 1
      IFLGDO = IFLGIN
      LEOF = .FALSE.
      INCOMP = 0
C                                           . . . . read next record
   10 CONTINUE
      IF (ISW(6) .EQ. 1) THEN
           WRITE (ISYSWR,'(A)') CPROMT(IFLGDO)
           IF (IFLGDO .EQ. 2)  LPHEAD = .FALSE.
      ENDIF
      CRDBUF = '   '
      READ (ISYSRD,'(A)',ERR=500,END=45)  CRDBUF
C
C                 CUPBUF is the first few characters in upper case
      CUPBUF(1:10) = CRDBUF(1:10)
      DO 12 I= 1, 10
      IF (CRDBUF(I:I) .EQ. '''') GO TO 13
         DO 11 IC= 1, 26
         IF (CRDBUF(I:I) .EQ. CLOWER(IC:IC)) CUPBUF(I:I)=CUPPER(IC:IC)
   11    CONTINUE
   12 CONTINUE
   13 CONTINUE
C                                           . .   preemptive commands
      LEOF = .FALSE.
      IF (INDEX(CUPBUF,'*EOF') .EQ. 1)    THEN
         WRITE (ISYSWR,'(A,I3)') ' *EOF ENCOUNTERED ON UNIT NO.',ISYSRD
         LPHEAD = .TRUE.
         GO TO 50
         ENDIF
      IF (INDEX(CUPBUF,'SET INP') .EQ. 1)    THEN
         ICOMND = ICOMND + 1
         WRITE (ISYSWR, 21) ICOMND,CRDBUF(1:50)
   21    FORMAT (' **********'/' **',I5,' **',A/' **********')
         LPHEAD = .TRUE.
         GO TO 50
         ENDIF
      GO TO 80
C                                    . . hardware EOF on current ISYSRD
   45 CRDBUF = '*EOF '
      WRITE (ISYSWR,'(A,I3)') ' END OF DATA ON UNIT NO.',ISYSRD
C                                     or SET INPUT command
   50 CONTINUE
         CALL MNSTIN(CRDBUF,IERR)
         IF (IERR .EQ. 0)  GO TO 10
         IF (IERR .EQ. 2)  THEN
            IF (.NOT. LEOF) THEN
               WRITE (ISYSWR,'(A,A/)') ' TWO CONSECUTIVE EOFs ON ',
     +              'PRIMARY INPUT FILE WILL TERMINATE EXECUTION.'
               LEOF = .TRUE.
               GO TO 10
            ENDIF
         ENDIF
         IFLGUT = IERR
         GO TO 900
   80 IF (IFLGDO .GT. 1) GO TO 100
C                            read title        . . . . .   IFLGDO = 1
C              if title is 'SET TITLE', skip and read again
      IF (INDEX(CUPBUF,'SET TIT') .EQ. 1)  GO TO 10
      CALL MNSETI(CRDBUF(1:50))
      WRITE (ISYSWR,'(1X,A50)')  CTITL
      WRITE (ISYSWR,'(1X,78(1H*))')
         LPHEAD = .TRUE.
      IF (IFLGIN .EQ. IFLGDO)  GO TO 900
      IFLGDO = IFLGIN
      GO TO 10
C                            data record is not a title.
  100 CONTINUE
      IF (IFLGDO .GT. 2)  GO TO 300
C                          expect parameter definitions.   IFLGDO = 2
C              if parameter def is 'PARAMETER', skip and read again
      IF (INDEX(CUPBUF,'PAR') .EQ. 1)  GO TO 10
C              if line starts with SET TITLE, read a title first
      IF (INDEX(CUPBUF,'SET TIT') .EQ. 1)  THEN
         IFLGDO = 1
         GO TO 10
         ENDIF
C                      we really have parameter definitions now
      CALL MNPARS(CRDBUF,ICONDP)
      IF (ICONDP .EQ. 0)  GO TO 10
C          format error
      IF (ICONDP .EQ. 1)  THEN
         IF (ISW(6) .EQ. 1)  THEN
           WRITE (ISYSWR,'(A)') ' FORMAT ERROR.  IGNORED.  ENTER AGAIN.'
           GO TO 10
         ELSE
           WRITE (ISYSWR,'(A)') ' ERROR IN PARAMETER DEFINITION'
           IFLGUT = 4
           GO TO 900
         ENDIF
      ENDIF
C                     ICONDP = 2            . . . end parameter requests
      IF (ISW(5).GE.0 .AND. ISW(6).LT.1) WRITE (ISYSWR,'(4X,75(1H*))')
      LPHEAD = .TRUE.
      IF (IFLGIN .EQ. IFLGDO)  GO TO 900
      IFLGDO = IFLGIN
      GO TO 10
C                                              . . . . .   IFLGDO = 3
C                                           read commands
  300 CONTINUE
      CALL MNCOMD(FCN,CRDBUF,ICONDN,FUTIL)
CC     ICONDN = 0: command executed normally
CC              1: command is blank, ignored
CC              2: command line unreadable, ignored
CC              3: unknown command, ignored
CC              4: abnormal termination (e.g., MIGRAD not converged)
CC              5: command is a request to read PARAMETER definitions
CC              6: 'SET INPUT' command
CC              7: 'SET TITLE' command
CC              8: 'SET COVAR' command
CC              9: reserved
CC             10: END command
CC             11: EXIT or STOP command
CC             12: RETURN command
      IF (ICONDN .EQ. 2 .OR. ICONDN .EQ. 3) THEN
         INCOMP = INCOMP + 1
         IF (INCOMP .GT. 100) THEN
            IFLGUT = 5
            GO TO 900
            ENDIF
         ENDIF
C                         parameter
      IF (ICONDN .EQ. 5)  IFLGDO = 2
C                         SET INPUT
      IF (ICONDN .EQ. 6)  GO TO 50
C                         SET TITLE
      IF (ICONDN .EQ. 7)  IFLGDO = 1
C                                        . . . . . . . . . . set covar
      IF (ICONDN .EQ. 8) THEN
         ICOMND = ICOMND + 1
         WRITE (ISYSWR,405) ICOMND,CRDBUF(1:50)
  405    FORMAT (1H ,10(1H*)/' **',I5,' **',A)
         WRITE (ISYSWR, '(1H ,10(1H*))' )
         NPAR2 = NPAR*(NPAR+1)/2
         READ (ISYSRD,420,ERR=500,END=45)  (VHMAT(I),I=1,NPAR2)
  420    FORMAT (BN,7E11.4,3X)
         ISW(2) = 3
         DCOVAR = 0.0
         IF (ISW(5) .GE. 0)  CALL MNMATU(1)
         IF (ISW(5) .GE. 1)  CALL MNPRIN(2,AMIN)
         GO TO 10
         ENDIF
      IF (ICONDN .LT. 10) GO TO 10
      GO TO 900
C                                              . . . . error conditions
  500 IFLGUT = 3
  900 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNRN15(VAL,INSEED)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C         This is a super-portable random number generator.
C         It should not overflow on any 32-bit machine.
C         The cycle is only ~10**9, so use with care!
C         Note especially that VAL must not be undefined on input.
C                    Set Default Starting Seed
      PARAMETER (THREE=3.0)
      DATA ISEED/12345/
      IF (VAL .EQ. THREE)  GO TO 100
C
      INSEED = ISEED
      K = ISEED/53668
      ISEED = 40014*(ISEED-K*53668) - K*12211
      IF (ISEED .LT. 0) ISEED = ISEED + 2147483563
      VAL = REAL(ISEED) * 4.656613E-10
      RETURN
C               "entry" to set seed, flag is VAL=3.
  100 ISEED = INSEED
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNRSET(IOPT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MNCLER and whenever problem changes, for example
CC        after SET LIMITS, SET PARAM, CALL FCN 6
CC    If IOPT=1,
CC        Resets function value and errors to UNDEFINED
CC    If IOPT=0, sets only MINOS errors to undefined
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CSTATU = 'RESET     '
      IF (IOPT .GE. 1)  THEN
        AMIN = UNDEFI
        FVAL3 = 2.0*ABS(AMIN) + 1.
        EDM = BIGEDM
        ISW(4) = 0
        ISW(2) = 0
        DCOVAR = 1.
        ISW(1) = 0
      ENDIF
      LNOLIM = .TRUE.
      DO 10 I= 1, NPAR
      IEXT = NEXOFI(I)

c ---------------------------------------
c protect against IEXT=0 (RK, Mar 2015)
      IF ( IEXT > 0 ) THEN
         IF (NVARL(IEXT) .GE. 4) LNOLIM=.FALSE.   
      ENDIF

      ERP(I) = ZERO
      ERN(I) = ZERO
      GLOBCC(I) = ZERO
   10 CONTINUE
      IF (ISW(2) .GE. 1)  THEN
         ISW(2) = 1
         DCOVAR = MAX(DCOVAR,HALF)
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSAVE
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Writes current parameter values and step sizes onto file ISYSSA
CC          in format which can be reread by Minuit for restarting.
CC       The covariance matrix is also output if it exists.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      LOGICAL LOPEN,LNAME
      CHARACTER CGNAME*64, CFNAME*64, CANSWR*1
C
      INQUIRE(UNIT=ISYSSA,OPENED=LOPEN,NAMED=LNAME,NAME=CGNAME)
      IF (LOPEN) THEN
         IF (.NOT.LNAME) CGNAME='UNNAMED FILE'
         WRITE (ISYSWR,32) ISYSSA,CGNAME
   32    FORMAT (' CURRENT VALUES WILL BE SAVED ON UNIT',I3,': ',A/)
      ELSE
C                new file, open it
         WRITE (ISYSWR,35) ISYSSA
   35    FORMAT (' UNIT',I3,' IS NOT OPENED.')
         IF (ISW(6) .EQ. 1) THEN
            WRITE (ISYSWR,'(A)') ' PLEASE GIVE FILE NAME:'
            READ (ISYSRD,'(A)') CFNAME
            OPEN (UNIT=ISYSSA,FILE=CFNAME,STATUS='NEW',ERR=600)
            CGNAME = CFNAME
         ELSE
            GO TO 650
         ENDIF
      ENDIF
C                               file is now correctly opened
      IF (ISW(6) .EQ. 1)  THEN
         WRITE (ISYSWR,37)  ISYSSA
   37    FORMAT (' SHOULD UNIT',I3,' BE REWOUND BEFORE WRITING TO IT?' )
         READ  (ISYSRD,'(A)')  CANSWR
         IF (CANSWR.EQ.'Y' .OR. CANSWR.EQ.'y') REWIND ISYSSA
      ENDIF
C                               and rewound if requested
      WRITE (ISYSSA,'(10HSET TITLE )',ERR=700)
      WRITE (ISYSSA,'(A)')  CTITL
      WRITE (ISYSSA,'(10HPARAMETERS)')
      NLINES = 3
C                                write out parameter values
      DO 200 I= 1, NU
      IF (NVARL(I) .LT. 0)  GO TO 200
      NLINES = NLINES + 1
      IINT = NIOFEX(I)
      IF (NVARL(I) .GT. 1)  GO TO 100
C         parameter without limits
      WRITE (ISYSSA,1001)  I,CPNAM(I),U(I),WERR(IINT)
      GO TO 200
C         parameter with limits
  100 CONTINUE
      WRITE (ISYSSA,1001) I,CPNAM(I),U(I),WERR(IINT),ALIM(I),BLIM(I)
 1001 FORMAT (1X,I5,1H',A10,1H',4E13.5)
  200 CONTINUE
      WRITE (ISYSSA,'(A)')  ' '
      NLINES = NLINES + 1
C                                  write out covariance matrix, if any
      IF (ISW(2) .LT. 1)  GO TO 750
      WRITE (ISYSSA,1003,ERR=700)  NPAR
 1003 FORMAT ('SET COVARIANCE',I6)
      NPAR2 = NPAR*(NPAR+1)/2
      WRITE (ISYSSA,1004) (VHMAT(I),I=1,NPAR2)
 1004 FORMAT (BN,7E11.4,3X)
      NCOVAR = NPAR2/7 + 1
      IF (MOD(NPAR2,7) .GT. 0)  NCOVAR = NCOVAR + 1
      NLINES = NLINES + NCOVAR
      WRITE (ISYSWR, 501) NLINES,ISYSSA,CGNAME(1:45)
  501 FORMAT (1X,I5,' RECORDS WRITTEN TO UNIT',I4,':',A)
      IF (NCOVAR .GT. 0) WRITE (ISYSWR, 502) NCOVAR
  502 FORMAT (' INCLUDING',I5,' RECORDS FOR THE COVARIANCE MATRIX.'/)
      GO TO 900
C                                           some error conditions
  600 WRITE (ISYSWR,'(A,I4)') ' I/O ERROR: UNABLE TO OPEN UNIT',ISYSSA
      GO TO 900
  650 WRITE (ISYSWR,'(A,I4,A)') ' UNIT',ISYSSA,' IS NOT OPENED.'
      GO TO 900
  700 WRITE (ISYSWR,'(A,I4)') ' ERROR: UNABLE TO WRITE TO UNIT',ISYSSA
      GO TO 900
  750 WRITE (ISYSWR,'(A)') ' THERE IS NO COVARIANCE MATRIX TO SAVE.'
C
  900 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:51  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSCAN(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Scans the values of FCN as a function of one parameter
CC        and plots the resulting values as a curve using MNPLOT.
CC        It may be called to scan one parameter or all parameters.
CC        retains the best function and parameter values found.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      XLREQ = MIN(WORD7(3),WORD7(4))
      XHREQ = MAX(WORD7(3),WORD7(4))
      NCALL = WORD7(2) + 0.01
      IF (NCALL .LE. 1)  NCALL = 41
      IF (NCALL .GT. MAXCPT)  NCALL = MAXCPT
      NCCALL = NCALL
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      IPARWD = WORD7(1) + 0.1
      IPAR = MAX(IPARWD, 0)
      IINT = NIOFEX(IPAR)
      CSTATU = 'NO CHANGE'
      IF (IPARWD .GT. 0)  GO TO 200
C
C         equivalent to a loop over parameters requested
  100 IPAR = IPAR + 1
      IF (IPAR .GT. NU)  GO TO 900
      IINT = NIOFEX(IPAR)
      IF (IINT .LE. 0)  GO TO 100
C         set up range for parameter IPAR
  200 CONTINUE
      UBEST = U(IPAR)
      XPT(1) = UBEST
      YPT(1) = AMIN
      CHPT(1)= ' '
      XPT(2) = UBEST
      YPT(2) = AMIN
      CHPT(2)= 'X'
      NXYPT = 2
      IF (NVARL(IPAR) .GT. 1)  GO TO 300
C         no limits on parameter
      IF (XLREQ .EQ. XHREQ)  GO TO 250
      UNEXT = XLREQ
      STEP = (XHREQ-XLREQ)/FLOAT(NCALL-1)
      GO TO 500
  250 CONTINUE
      XL = UBEST - WERR(IINT)
      XH = UBEST+  WERR(IINT)
      CALL MNBINS(XL,XH,NCALL, UNEXT,UHIGH,NBINS,STEP)
      NCCALL = NBINS + 1
      GO TO 500
C         limits on parameter
  300 CONTINUE
      IF (XLREQ .EQ. XHREQ)  GO TO 350
      XL = MAX(XLREQ,ALIM(IPAR))
      XH = MIN(XHREQ,BLIM(IPAR))
      IF (XL .GE. XH)  GO TO 700
      UNEXT = XL
      STEP = (XH-XL)/FLOAT(NCALL-1)
      GO TO 500
  350 CONTINUE
      UNEXT = ALIM(IPAR)
      STEP = (BLIM(IPAR)-ALIM(IPAR))/FLOAT(NCALL-1)
C         main scanning loop over parameter IPAR
  500 CONTINUE
      DO 600 ICALL = 1, NCCALL
      U(IPAR) = UNEXT
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FNEXT,U,4,FUTIL)
      NFCN = NFCN + 1
      NXYPT = NXYPT + 1
      XPT(NXYPT) = UNEXT
      YPT(NXYPT) = FNEXT
      CHPT(NXYPT) = '*'
      IF (FNEXT .LT. AMIN)  THEN
        AMIN = FNEXT
        UBEST = UNEXT
        CSTATU= 'IMPROVED  '
        ENDIF
  530 CONTINUE
      UNEXT = UNEXT + STEP
  600 CONTINUE
C         finished with scan of parameter IPAR
      U(IPAR) = UBEST
      CALL MNEXIN(X)
      IF (ISW(5) .GE. 1)  THEN
        WRITE (ISYSWR,1001)  NEWPAG,IPAR,CPNAM(IPAR)
        NUNIT = ISYSWR
        CALL MNPLOT(XPT,YPT,CHPT,NXYPT,NUNIT,NPAGWD,NPAGLN)
      ENDIF
      GO TO 800
  700 CONTINUE
      WRITE (ISYSWR,1000) IPAR
  800 CONTINUE
      IF (IPARWD .LE. 0)  GO TO 100
C         finished with all parameters
  900 CONTINUE
      IF (ISW(5) .GE. 0) CALL MNPRIN(5,AMIN)
      RETURN
 1000 FORMAT (46H REQUESTED RANGE OUTSIDE LIMITS FOR PARAMETER  ,I3/)
 1001 FORMAT (I1,'SCAN OF PARAMETER NO.',I3,3H,   ,A10)
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSEEK(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC   Performs a rough (but global) minimization by monte carlo search.
CC        Each time a new minimum is found, the search area is shifted
CC        to be centered at the best value.  Random points are chosen
CC        uniformly over a hypercube determined by current step sizes.
CC   The Metropolis algorithm accepts a worse point with probability
CC      exp(-d/UP), where d is the degradation.  Improved points
CC      are of course always accepted.  Actual steps are random
CC      multiples of the nominal steps (DIRIN).
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      PARAMETER (TWOPI=2.0*3.141593)
      DIMENSION  XBEST(MNI), XMID(MNI)
      MXFAIL = WORD7(1)
      IF (MXFAIL .LE. 0)  MXFAIL=100+20*NPAR
      MXSTEP = 10*MXFAIL
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      ALPHA = WORD7(2)
      IF (ALPHA .LE. ZERO)  ALPHA=3.
      IF (ISW(5) .GE. 1)  WRITE (ISYSWR, 3) MXFAIL,MXSTEP,ALPHA
    3 FORMAT (' MNSEEK: MONTE CARLO MINIMIZATION USING METROPOLIS',
     + ' ALGORITHM'/' TO STOP AFTER',I6,' SUCCESSIVE FAILURES, OR',
     + I7,' STEPS'/' MAXIMUM STEP SIZE IS',F9.3,' ERROR BARS.')
      CSTATU= 'INITIAL  '
      IF (ISW(5) .GE. 2)  CALL MNPRIN(2,AMIN)
      CSTATU = 'UNCHANGED '
      IFAIL = 0
      RNUM = ZERO
      RNUM1 = ZERO
      RNUM2 = ZERO
      NPARX = NPAR
      FLAST = AMIN
C              set up step sizes, starting values
      DO 10 IPAR =  1, NPAR
      IEXT = NEXOFI(IPAR)
      DIRIN(IPAR) = 2.0*ALPHA*WERR(IPAR)
      IF (NVARL(IEXT) .GT. 1)  THEN
C              parameter with limits
         CALL MNDXDI(X(IPAR),IPAR,DXDI)
         IF (DXDI .EQ. ZERO)  DXDI=1.
         DIRIN(IPAR) = 2.0*ALPHA*WERR(IPAR)/DXDI
         IF (ABS(DIRIN(IPAR)).GT.TWOPI)  DIRIN(IPAR)=TWOPI
         ENDIF
      XMID(IPAR) = X(IPAR)
   10 XBEST(IPAR) = X(IPAR)
C                              search loop
      DO 500 ISTEP= 1, MXSTEP
      IF (IFAIL .GE. MXFAIL)  GO TO 600
        DO 100 IPAR= 1, NPAR
        CALL MNRN15(RNUM1,ISEED)
        CALL MNRN15(RNUM2,ISEED)
  100   X(IPAR) = XMID(IPAR) + 0.5*(RNUM1+RNUM2-1.)*DIRIN(IPAR)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FTRY,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FTRY .LT. FLAST)  THEN
         IF (FTRY .LT. AMIN)  THEN
            CSTATU = 'IMPROVEMNT'
            AMIN = FTRY
            DO 200 IB= 1, NPAR
  200       XBEST(IB) = X(IB)
            IFAIL = 0
            IF (ISW(5) .GE. 2) CALL MNPRIN(2,AMIN)
            ENDIF
         GO TO 300
      ELSE
         IFAIL = IFAIL + 1
C                   Metropolis algorithm
         BAR = (AMIN-FTRY)/UP
         CALL MNRN15(RNUM,ISEED)
         IF (BAR .LT. LOG(RNUM))  GO TO 500
      ENDIF
C                    Accept new point, move there
  300 CONTINUE
      DO 350 J= 1, NPAR
      XMID(J) = X(J)
  350 CONTINUE
      FLAST = FTRY
  500 CONTINUE
C                               end search loop
  600 CONTINUE
      IF (ISW(5) .GT. 1) WRITE (ISYSWR,601) IFAIL
  601 FORMAT(' MNSEEK:',I5,' SUCCESSIVE UNSUCCESSFUL TRIALS.')
      DO 700 IB= 1, NPAR
  700 X(IB) = XBEST(IB)
      CALL MNINEX(X)
      IF (ISW(5) .GE. 1)  CALL MNPRIN(2,AMIN)
      IF (ISW(5) .EQ. 0)  CALL MNPRIN(0,AMIN)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:52  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNSET(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Called from MNEXCM
CC        Interprets the commands that start with SET and SHOW
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C
      EXTERNAL FCN,FUTIL
C        file characteristics for SET INPUT
      LOGICAL LNAME
      CHARACTER CFNAME*64, CMODE*16
C       'SET ' or 'SHOW',  'ON ' or 'OFF', 'SUPPRESSED' or 'REPORTED  '
      CHARACTER CKIND*4,    COPT*3,         CWARN*10
C        explanation of print level numbers -1:3  and strategies 0:2
      CHARACTER CPRLEV(-1:3)*34 ,CSTRAT(0:2)*44
C        identification of debug options
      PARAMETER (NUMDBG = 6)
      CHARACTER*40 CDBOPT(0:NUMDBG)
C        things that can be set or shown
      CHARACTER*10 CNAME(30)
      DATA CNAME( 1)/'FCN value '/
      DATA CNAME( 2)/'PARameters'/
      DATA CNAME( 3)/'LIMits    '/
      DATA CNAME( 4)/'COVariance'/
      DATA CNAME( 5)/'CORrelatio'/
      DATA CNAME( 6)/'PRInt levl'/
      DATA CNAME( 7)/'NOGradient'/
      DATA CNAME( 8)/'GRAdient  '/
      DATA CNAME( 9)/'ERRor def '/
      DATA CNAME(10)/'INPut file'/
      DATA CNAME(11)/'WIDth page'/
      DATA CNAME(12)/'LINes page'/
      DATA CNAME(13)/'NOWarnings'/
      DATA CNAME(14)/'WARnings  '/
      DATA CNAME(15)/'RANdom gen'/
      DATA CNAME(16)/'TITle     '/
      DATA CNAME(17)/'STRategy  '/
      DATA CNAME(18)/'EIGenvalue'/
      DATA CNAME(19)/'PAGe throw'/
      DATA CNAME(20)/'MINos errs'/
      DATA CNAME(21)/'EPSmachine'/
      DATA CNAME(22)/'OUTputfile'/
      DATA CNAME(23)/'BATch     '/
      DATA CNAME(24)/'INTeractiv'/
      DATA CNAME(25)/'VERsion   '/
          DATA NNAME/25/
C        options not intended for normal users
      DATA CNAME(26)/'reserve   '/
      DATA CNAME(27)/'NODebug   '/
      DATA CNAME(28)/'DEBug     '/
      DATA CNAME(29)/'SHOw      '/
      DATA CNAME(30)/'SET       '/
          DATA NNTOT/30/
C
      DATA CPRLEV(-1)/'-1: NO OUTPUT EXCEPT FROM "SHOW"  '/
      DATA CPRLEV( 0)/' 0: REDUCED OUTPUT                '/
      DATA CPRLEV( 1)/' 1: NORMAL OUTPUT                 '/
      DATA CPRLEV( 2)/' 2: EXTRA OUTPUT FOR PROBLEM CASES'/
      DATA CPRLEV( 3)/' 3: MAXIMUM OUTPUT                '/
C
      DATA CSTRAT( 0)/' 0: MINIMIZE THE NUMBER OF CALLS TO FUNCTION'/
      DATA CSTRAT( 1)/' 1: TRY TO BALANCE SPEED AGAINST RELIABILITY'/
      DATA CSTRAT( 2)/' 2: MAKE SURE MINIMUM TRUE, ERRORS CORRECT  '/
C
      DATA CDBOPT(0)/'REPORT ALL EXCEPTIONAL CONDITIONS      '/
      DATA CDBOPT(1)/'MNLINE: LINE SEARCH MINIMIZATION       '/
      DATA CDBOPT(2)/'MNDERI: FIRST DERIVATIVE CALCULATIONS  '/
      DATA CDBOPT(3)/'MNHESS: SECOND DERIVATIVE CALCULATIONS '/
      DATA CDBOPT(4)/'MNMIGR: COVARIANCE MATRIX UPDATES      '/
      DATA CDBOPT(5)/'MNHES1: FIRST DERIVATIVE UNCERTAINTIES '/
      DATA CDBOPT(6)/'MNCONT: MNCONTOUR PLOT (MNCROS SEARCH) '/
C
C
      DO 2 I= 1, NNTOT
      IF (INDEX(CWORD(4:10),CNAME(I)(1:3)) .GT. 0)  GO TO 5
    2 CONTINUE
      I = 0
    5 KNAME = I
C
C           Command could be SET xxx, SHOW xxx,  HELP SET or HELP SHOW
      IF (INDEX(CWORD(1:4),'HEL') .GT. 0)  GO TO 2000
      IF (INDEX(CWORD(1:4),'SHO') .GT. 0)  GO TO 1000
      IF (INDEX(CWORD(1:4),'SET') .EQ. 0)  GO TO 1900
C                           ---
      CKIND = 'SET '
C                                        . . . . . . . . . . set unknown
      IF (KNAME .LE. 0)  GO TO 1900
C                                        . . . . . . . . . . set known
      GO TO(3000,  20,  30,  40,3000,  60,  70,  80,  90, 100,
     +       110, 120, 130, 140, 150, 160, 170,3000, 190,3000,
     +       210, 220, 230, 240,3000,1900, 270, 280, 290, 300) , KNAME
C
C                                        . . . . . . . . . . set param
   20 CONTINUE
      IPRM = WORD7(1)
      IF (IPRM .GT. NU)  GO TO 25
      IF (IPRM .LE. 0)   GO TO 25
      IF (NVARL(IPRM) .LT. 0)  GO TO 25
      U(IPRM) = WORD7(2)
      CALL MNEXIN(X)
      ISW2 = ISW(2)
      CALL MNRSET(1)
C        Keep approximate covariance matrix, even if new param value
      ISW(2) = MIN(ISW2,1)
      CFROM = 'SET PARM'
      NFCNFR = NFCN
      CSTATU = 'NEW VALUES'
      GO TO 4000
   25 WRITE (ISYSWR,'(A/)') ' UNDEFINED PARAMETER NUMBER.  IGNORED.'
      GO TO 4000
C                                        . . . . . . . . . . set limits
   30 CALL MNLIMS
      GO TO 4000
C                                        . . . . . . . . . . set covar
   40 CONTINUE
C   this command must be handled by MNREAD, and is not Fortran-callable
      GO TO 3000
C                                        . . . . . . . . . . set print
   60 ISW(5) = WORD7(1)
      GO TO 4000
C                                        . . . . . . . . . . set nograd
   70 ISW(3) = 0
      GO TO 4000
C                                        . . . . . . . . . . set grad
   80 CALL MNGRAD(FCN,FUTIL)
      GO TO 4000
C                                        . . . . . . . . . . set errdef
   90 IF (WORD7(1) .EQ. UP)  GO TO 4000
      IF (WORD7(1) .LE. ZERO)  THEN
         IF (UP .EQ. UPDFLT)  GO TO 4000
         UP = UPDFLT
      ELSE
         UP = WORD7(1)
      ENDIF
      DO 95 I= 1, NPAR
      ERN(I) = 0.
   95 ERP(I) = 0.
      CALL MNWERR
      GO TO 4000
C                                        . . . . . . . . . . set input
C This command must be handled by MNREAD. If it gets this far,
C         it is illegal.
  100 CONTINUE
      GO TO 3000
C                                        . . . . . . . . . . set width
  110 NPAGWD = WORD7(1)
      NPAGWD = MAX(NPAGWD,50)
      GO TO 4000
C                                        . . . . . . . . . . set lines
  120 NPAGLN = WORD7(1)
      GO TO 4000
C                                        . . . . . . . . . . set nowarn
  130 LWARN = .FALSE.
      GO TO 4000
C                                        . . . . . . . . . . set warn
  140 LWARN = .TRUE.
      CALL MNWARN('W','SHO','SHO')
      GO TO 4000
C                                        . . . . . . . . . . set random
  150 JSEED = INT(WORD7(1))
      VAL = 3.
      CALL MNRN15(VAL, JSEED)
      IF (ISW(5) .GT. 0) WRITE (ISYSWR, 151) JSEED
  151 FORMAT (' MINUIT RANDOM NUMBER SEED SET TO ',I10)
      GO TO 4000
C                                        . . . . . . . . . . set title
  160 CONTINUE
C   this command must be handled by MNREAD, and is not Fortran-callable
      GO TO 3000
C                                        . . . . . . . . . set strategy
  170 ISTRAT = WORD7(1)
      ISTRAT = MAX(ISTRAT,0)
      ISTRAT = MIN(ISTRAT,2)
      IF (ISW(5) .GT. 0)  GO TO 1172
      GO TO 4000
C                                       . . . . . . . . . set page throw
  190 NEWPAG = WORD7(1)
      GO TO 1190
C                                        . . . . . . . . . . set epsmac
  210 IF (WORD7(1).GT.ZERO .AND. WORD7(1).LT.0.1) EPSMAC = WORD7(1)
      EPSMA2 = SQRT(EPSMAC)
      GO TO 1210
C                                        . . . . . . . . . . set outputfile
  220 CONTINUE
      IUNIT = WORD7(1)
      ISYSWR = IUNIT
      ISTKWR(1) = IUNIT
      IF (ISW(5) .GE. 0) GO TO 1220
      GO TO 4000
C                                        . . . . . . . . . . set batch
  230 ISW(6) = 0
      IF (ISW(5) .GE. 0)  GO TO 1100
      GO TO 4000
C                                        . . . . . . . . . . set interactive
  240 ISW(6) = 1
      IF (ISW(5) .GE. 0)  GO TO 1100
      GO TO 4000
C                                        . . . . . . . . . . set nodebug
  270 ISET = 0
      GO TO 281
C                                        . . . . . . . . . . set debug
  280 ISET = 1
  281 CONTINUE
      IDBOPT = WORD7(1)
      IF (IDBOPT .GT. NUMDBG) GO TO 288
      IF (IDBOPT .GE. 0) THEN
          IDBG(IDBOPT) = ISET
          IF (ISET .EQ. 1)  IDBG(0) = 1
      ELSE
C             SET DEBUG -1  sets all debug options
          DO 285 ID= 0, NUMDBG
  285     IDBG(ID) = ISET
      ENDIF
      LREPOR = (IDBG(0) .GE. 1)
      CALL MNWARN('D','SHO','SHO')
      GO TO 4000
  288 WRITE (ISYSWR,289) IDBOPT
  289 FORMAT (' UNKNOWN DEBUG OPTION',I6,' REQUESTED. IGNORED')
      GO TO 4000
C                                        . . . . . . . . . . set show
  290 CONTINUE
C                                        . . . . . . . . . . set set
  300 CONTINUE
      GO TO 3000
C                -----------------------------------------------------
 1000 CONTINUE
C               at this point, CWORD must be 'SHOW'
      CKIND = 'SHOW'
      IF (KNAME .LE. 0)  GO TO 1900
      GO TO (1010,1020,1030,1040,1050,1060,1070,1070,1090,1100,
     +       1110,1120,1130,1130,1150,1160,1170,1180,1190,1200,
     +       1210,1220,1100,1100,1250,1900,1270,1270,1290,1300),KNAME
C
C                                        . . . . . . . . . . show fcn
 1010 CONTINUE
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CALL MNPRIN (0,AMIN)
      GO TO 4000
C                                        . . . . . . . . . . show param
 1020 CONTINUE
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CALL MNPRIN (5,AMIN)
      GO TO 4000
C                                        . . . . . . . . . . show limits
 1030 CONTINUE
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CALL MNPRIN (1,AMIN)
      GO TO 4000
C                                        . . . . . . . . . . show covar
 1040 CALL MNMATU(1)
      GO TO 4000
C                                        . . . . . . . . . . show corre
 1050 CALL MNMATU(0)
      GO TO 4000
C                                        . . . . . . . . . . show print
 1060 CONTINUE
      IF (ISW(5) .LT.-1)  ISW(5) = -1
      IF (ISW(5) .GT. 3)  ISW(5) = 3
      WRITE (ISYSWR,'(A)') ' ALLOWED PRINT LEVELS ARE:'
      WRITE (ISYSWR,'(27X,A)') CPRLEV
      WRITE (ISYSWR,1061)  CPRLEV(ISW(5))
 1061 FORMAT (/' CURRENT PRINTOUT LEVEL IS ',A)
      GO TO 4000
C                                        . . . . . . . show nograd, grad
 1070 CONTINUE
      IF (ISW(3) .LE. 0) THEN
         WRITE (ISYSWR, 1081)
 1081    FORMAT(' NOGRAD IS SET.  DERIVATIVES NOT COMPUTED IN FCN.')
      ELSE
         WRITE (ISYSWR, 1082)
 1082    FORMAT('   GRAD IS SET.  USER COMPUTES DERIVATIVES IN FCN.')
      ENDIF
      GO TO 4000
C                                       . . . . . . . . . . show errdef
 1090 WRITE (ISYSWR, 1091)  UP
 1091 FORMAT (' ERRORS CORRESPOND TO FUNCTION CHANGE OF',G13.5)
      GO TO 4000
C                                       . . . . . . . . . . show input,
C                                                batch, or interactive
 1100 CONTINUE
      INQUIRE(UNIT=ISYSRD,NAMED=LNAME,NAME=CFNAME)
      CMODE = 'BATCH MODE      '
      IF (ISW(6) .EQ. 1)  CMODE = 'INTERACTIVE MODE'
      IF (.NOT. LNAME)  CFNAME='unknown'
      WRITE (ISYSWR,1002) CMODE,ISYSRD,CFNAME
 1002 FORMAT (' INPUT NOW BEING READ IN ',A,' FROM UNIT NO.',I3/
     + ' FILENAME: ',A)
      GO TO 4000
C                                       . . . . . . . . . . show width
 1110 WRITE (ISYSWR,1111) NPAGWD
 1111 FORMAT (10X,'PAGE WIDTH IS SET TO',I4,' COLUMNS')
      GO TO 4000
C                                       . . . . . . . . . . show lines
 1120 WRITE (ISYSWR,1121) NPAGLN
 1121 FORMAT (10X,'PAGE LENGTH IS SET TO',I4,' LINES')
      GO TO 4000
C                                       . . . . . . .show nowarn, warn
 1130 CONTINUE
                 CWARN = 'SUPPRESSED'
      IF (LWARN) CWARN = 'REPORTED  '
      WRITE (ISYSWR,1141) CWARN
 1141 FORMAT (' MINUIT WARNING MESSAGES ARE ',A)
      IF (.NOT. LWARN) CALL MNWARN('W','SHO','SHO')
      GO TO 4000
C                                      . . . . . . . . . . show random
 1150 VAL = 0.
      CALL MNRN15(VAL,IGRAIN)
      IKSEED = IGRAIN
      WRITE (ISYSWR, 1151)  IKSEED
 1151 FORMAT (' MINUIT RNDM SEED IS CURRENTLY=',I10/)
      VAL = 3.0
      ISEED = IKSEED
      CALL MNRN15(VAL,ISEED)
      GO TO 4000
C                                        . . . . . . . . . show title
 1160 WRITE (ISYSWR,'(A,A)') ' TITLE OF CURRENT TASK IS:',CTITL
      GO TO 4000
C                                        . . . . . . . show strategy
 1170 WRITE (ISYSWR, '(A)') ' ALLOWED STRATEGIES ARE:'
      WRITE (ISYSWR, '(20X,A)') CSTRAT
 1172 WRITE (ISYSWR, 1175) CSTRAT(ISTRAT)
 1175 FORMAT (/' NOW USING STRATEGY ',A/)
      GO TO 4000
C                                          . . . . . show eigenvalues
 1180 CONTINUE
      ISWSAV = ISW(5)
      ISW(5) = 3
      IF (ISW(2) .LT. 1)  THEN
         WRITE (ISYSWR,'(1X,A)') COVMES(0)
      ELSE
         CALL MNPSDF
      ENDIF
      ISW(5) = ISWSAV
      GO TO 4000
C                                            . . . . . show page throw
 1190 WRITE (ISYSWR,'(A,I3)') ' PAGE THROW CARRIAGE CONTROL =',NEWPAG
      IF (NEWPAG .EQ. 0)
     +    WRITE (ISYSWR,'(A)') ' NO PAGE THROWS IN MINUIT OUTPUT'
      GO TO 4000
C                                        . . . . . . show minos errors
 1200 CONTINUE
      DO 1202 II= 1, NPAR
      IF (ERP(II).GT.ZERO .OR. ERN(II).LT.ZERO)  GO TO 1204
 1202 CONTINUE
      WRITE (ISYSWR,'(A)')
     +   '       THERE ARE NO MINOS ERRORS CURRENTLY VALID.'
      GO TO 4000
 1204 CONTINUE
      CALL MNPRIN(4,AMIN)
      GO TO 4000
C                                        . . . . . . . . . show epsmac
 1210 WRITE (ISYSWR,'(A,E12.3)')
     +  ' FLOATING-POINT NUMBERS ASSUMED ACCURATE TO',EPSMAC
      GO TO 4000
C                                        . . . . . . show outputfiles
 1220 CONTINUE
      WRITE (ISYSWR,'(A,I4)') '  MINUIT PRIMARY OUTPUT TO UNIT',ISYSWR
      GO TO 4000
C                                        . . . . . . show version
 1250 CONTINUE
      WRITE (ISYSWR,'(A,A)') ' THIS IS MINUIT VERSION:',CVRSN
      GO TO 4000
C                                        . . . . . . show nodebug, debug
 1270 CONTINUE
      DO 1285 ID= 0, NUMDBG
      COPT = 'OFF'
      IF (IDBG(ID) .GE. 1)  COPT = 'ON '
 1285 WRITE (ISYSWR,1286) ID, COPT, CDBOPT(ID)
 1286 FORMAT (10X,'DEBUG OPTION',I3,' IS ',A3,' :',A)
      IF (.NOT. LREPOR) CALL MNWARN('D','SHO','SHO')
      GO TO 4000
C                                        . . . . . . . . . . show show
 1290 CKIND = 'SHOW'
      GO TO 2100
C                                        . . . . . . . . . . show set
 1300 CKIND = 'SET '
      GO TO 2100
C                -----------------------------------------------------
C                              UNKNOWN COMMAND
 1900 WRITE (ISYSWR, 1901) CWORD
 1901 FORMAT (' THE COMMAND:',A10,' IS UNKNOWN.'/)
      GO TO 2100
C                -----------------------------------------------------
C                    HELP SHOW,  HELP SET,  SHOW SET, or SHOW SHOW
 2000 CKIND = 'SET '
      IF (INDEX(CWORD(4:10),'SHO') .GT. 0)  CKIND = 'SHOW'
 2100 WRITE (ISYSWR, 2101)  CKIND,CKIND, (CNAME(KK),KK=1,NNAME)
 2101 FORMAT (' THE FORMAT OF THE ',A4,' COMMAND IS:'//
     +   1X,A4,' xxx    [numerical arguments if any]'//
     +   ' WHERE xxx MAY BE ONE OF THE FOLLOWING:'/
     +   (7X,6A12))
      GO TO 4000
C                -----------------------------------------------------
C                               ILLEGAL COMMAND
 3000 WRITE (ISYSWR,'('' ABOVE COMMAND IS ILLEGAL.   IGNORED'')')
 4000 RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSETI(TIT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       Called by user to set or change title of current task.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER*(*) TIT
      CTITL = TIT
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:54  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSIMP(FCN,FUTIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Performs a minimization using the simplex method of Nelder
CC        and Mead (ref. -- Comp. J. 7,308 (1965)).
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION Y(MNI+1)
      DATA ALPHA,BETA,GAMMA,RHOMIN,RHOMAX / 1.0, 0.5, 2.0, 4.0, 8.0/
      IF (NPAR .LE. 0)  RETURN
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CFROM = 'SIMPLEX '
      NFCNFR = NFCN
      CSTATU= 'UNCHANGED '
      NPFN=NFCN
      NPARP1=NPAR+1
      NPARX = NPAR
      RHO1 = 1.0 + ALPHA
      RHO2 = RHO1 + ALPHA*GAMMA
      WG = 1.0/FLOAT(NPAR)
      IF (ISW(5) .GE. 0) WRITE(ISYSWR,100) EPSI
  100 FORMAT(' START SIMPLEX MINIMIZATION.    CONVERGENCE WHEN EDM .LT.'
     +,E10.2 )
         DO 2 I= 1, NPAR
         DIRIN(I) = WERR(I)
           CALL MNDXDI(X(I),I,DXDI)
           IF (DXDI .NE. ZERO) DIRIN(I)=WERR(I)/DXDI
         DMIN = EPSMA2*ABS(X(I))
         IF (DIRIN(I) .LT. DMIN)  DIRIN(I)=DMIN
    2    CONTINUE
C**       choose the initial simplex using single-parameter searches
    1 CONTINUE
      YNPP1 = AMIN
      JL = NPARP1
      Y(NPARP1) = AMIN
      ABSMIN = AMIN
      DO 10 I= 1, NPAR
      AMING = AMIN
      PBAR(I) = X(I)
      BESTX = X(I)
      KG = 0
      NS = 0
      NF = 0
    4 X(I) = BESTX + DIRIN(I)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN, F, U, 4, FUTIL)
      NFCN = NFCN + 1
      IF (F .LT. AMING)  GO TO 6
C         failure
      IF (KG .EQ. 1)  GO TO 8
      KG = -1
      NF = NF + 1
      DIRIN(I) = DIRIN(I) * (-0.4)
      IF (NF .LT. 3)  GO TO 4
C         stop after three failures
      BESTX = X(I)
      DIRIN(I) = DIRIN(I) * 3.0
      AMING = F
      GO TO 8
C
C         success
    6 BESTX = X(I)
      DIRIN(I) = DIRIN(I) * 3.0
      AMING = F
      CSTATU= 'PROGRESS  '
      KG = 1
      NS = NS + 1
      IF (NS .LT. 6)  GO TO 4
C
C         3 failures or 6 successes or
C         local minimum found in ith direction
    8 Y(I) = AMING
      IF (AMING .LT. ABSMIN)  JL = I
      IF (AMING .LT. ABSMIN)  ABSMIN = AMING
      X(I) = BESTX
      DO 9 K= 1, NPAR
    9 P(K,I) = X(K)
   10 CONTINUE
      JH = NPARP1
      AMIN=Y(JL)
      CALL MNRAZZ(YNPP1,PBAR,Y,JH,JL)
      DO 20 I= 1, NPAR
   20 X(I) = P(I,JL)
      CALL MNINEX(X)
      IF (ISW(5) .GE. 1)  CALL MNPRIN(5,AMIN)
      EDM = BIGEDM
      SIG2 = EDM
      NCYCL=0
C                                        . . . . .  start main loop
   50 CONTINUE
      IF (SIG2 .LT. EPSI .AND. EDM.LT.EPSI)     GO TO 76
      SIG2 = EDM
      IF ((NFCN-NPFN) .GT. NFCNMX)  GO TO 78
C         calculate new point * by reflection
      DO 60 I= 1, NPAR
      PB = 0.
      DO 59 J= 1, NPARP1
   59 PB = PB + WG * P(I,J)
      PBAR(I) = PB - WG * P(I,JH)
   60 PSTAR(I)=(1.+ALPHA)*PBAR(I)-ALPHA*P(I,JH)
      CALL MNINEX(PSTAR)
      CALL FCN(NPARX,GIN,YSTAR,U,4,FUTIL)
      NFCN=NFCN+1
      IF(YSTAR.GE.AMIN) GO TO 70
C         point * better than jl, calculate new point **
      CSTATU = 'PROGRESS  '
      DO 61 I=1,NPAR
   61 PSTST(I)=GAMMA*PSTAR(I)+(1.-GAMMA)*PBAR(I)
      CALL MNINEX(PSTST)
      CALL FCN(NPARX,GIN,YSTST,U,4,FUTIL)
      NFCN=NFCN+1
C         try a parabola through ph, pstar, pstst.  min = prho
      Y1 = (YSTAR-Y(JH)) * RHO2
      Y2 = (YSTST-Y(JH)) * RHO1
      RHO = 0.5 * (RHO2*Y1 -RHO1*Y2) / (Y1 -Y2)
      IF (RHO .LT. RHOMIN)  GO TO 66
      IF (RHO .GT. RHOMAX)  RHO = RHOMAX
      DO 64 I= 1, NPAR
   64 PRHO(I) = RHO*PBAR(I) + (1.0-RHO)*P(I,JH)
      CALL MNINEX(PRHO)
      CALL FCN(NPARX,GIN,YRHO, U,4,FUTIL)
      NFCN = NFCN + 1
      IF (YRHO .LT. AMIN)     CSTATU = 'PROGRESS  '
      IF (YRHO .LT. Y(JL) .AND. YRHO .LT. YSTST)  GO TO 65
      IF (YSTST .LT. Y(JL))  GO TO 67
      IF (YRHO .GT. Y(JL))  GO TO 66
C         accept minimum point of parabola, PRHO
   65 CALL MNRAZZ (YRHO,PRHO,Y,JH,JL)
      GO TO 68
   66 IF (YSTST .LT. Y(JL))  GO TO 67
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      GO TO 68
   67 CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
   68 NCYCL=NCYCL+1
      IF (ISW(5) .LT. 2)  GO TO 50
      IF (ISW(5) .GE. 3 .OR. MOD(NCYCL, 10) .EQ. 0) CALL MNPRIN(5,AMIN)
      GO TO 50
C         point * is not as good as jl
   70 IF (YSTAR .GE. Y(JH))  GO TO 73
      JHOLD = JH
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      IF (JHOLD .NE. JH)  GO TO 50
C         calculate new point **
   73 DO 74 I=1,NPAR
   74 PSTST(I)=BETA*P(I,JH)+(1.-BETA)*PBAR(I)
      CALL MNINEX (PSTST)
      CALL FCN(NPARX,GIN,YSTST,U,4,FUTIL)
      NFCN=NFCN+1
      IF(YSTST.GT.Y(JH)) GO TO 1
C     point ** is better than jh
      IF (YSTST .LT. AMIN)     CSTATU = 'PROGRESS  '
      IF (YSTST .LT. AMIN)  GO TO 67
      CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
      GO TO 50
C                                        . . . . . .  end main loop
   76 IF (ISW(5) .GE. 0)  WRITE(ISYSWR,'(A)')
     +                    ' SIMPLEX MINIMIZATION HAS CONVERGED.'
      ISW(4) = 1
      GO TO 80
   78 IF (ISW(5) .GE. 0)  WRITE(ISYSWR,'(A)')
     +                    ' SIMPLEX TERMINATES WITHOUT CONVERGENCE.'
      CSTATU= 'CALL LIMIT'
      ISW(4) = -1
      ISW(1) = 1
   80 DO 82 I=1,NPAR
      PB = 0.
      DO 81 J=1,NPARP1
   81 PB = PB + WG * P(I,J)
   82 PBAR(I) = PB - WG * P(I,JH)
      CALL MNINEX(PBAR)
      CALL FCN(NPARX,GIN,YPBAR,U,4,FUTIL)
      NFCN=NFCN+1
      IF (YPBAR .LT. AMIN)  CALL MNRAZZ(YPBAR,PBAR,Y,JH,JL)
      CALL MNINEX(X)
      IF (NFCNMX+NPFN-NFCN .LT. 3*NPAR)  GO TO 90
      IF (EDM .GT. 2.0*EPSI)  GO TO 1
   90 IF (ISW(5) .GE. 0)  CALL MNPRIN(5, AMIN)
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC       User-called
CC       Provides the user with information concerning the current status
CC          of the current minimization. Namely, it returns:
CC        FMIN: the best function value found so far
CC        FEDM: the estimated vertical distance remaining to minimum
CC        ERRDEF: the value of UP defining parameter uncertainties
CC        NPARI: the number of currently variable parameters
CC        NPARX: the highest (external) parameter number defined by user
CC        ISTAT: a status integer indicating how good is the covariance
CC           matrix:  0= not calculated at all
CC                    1= approximation only, not accurate
CC                    2= full matrix, but forced positive-definite
CC                    3= full accurate covariance matrix
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      FMIN = AMIN
      FEDM = EDM
      ERRDEF = UP
      NPARI = NPAR
      NPARX = NU
      ISTAT = ISW(2)
        IF (EDM  .EQ. BIGEDM)  THEN
            FEDM = UP
        ENDIF
        IF (AMIN .EQ. UNDEFI)  THEN
            FMIN = 0.0
            FEDM = UP
            ISTAT= 0
        ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNSTIN(CRDBUF,IERR)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC Called from MNREAD.
CC Implements the SET INPUT command to change input units.
CC If command is: 'SET INPUT'   'SET INPUT 0'   or  '*EOF',
CC                 or 'SET INPUT , ,  ',
CC                reverts to previous input unit number,if any.
CC
CC      If it is: 'SET INPUT n'  or  'SET INPUT n filename',
CC                changes to new input file, added to stack
CC
CC      IERR = 0: reading terminated normally
CC             2: end-of-data on primary input file
CC             3: unrecoverable read error
CC             4: unable to process request
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER CRDBUF*(*),CUNIT*10,CFNAME*64,CGNAME*64,CANSWR*1
      CHARACTER CMODE*16
      LOGICAL LOPEN,LREWIN,NONAME,LNAME,MNUNPT
      NONAME = .TRUE.
      IERR = 0
      IF (INDEX(CRDBUF,'*EOF') .EQ. 1) GO TO 190
      IF (INDEX(CRDBUF,'*eof') .EQ. 1) GO TO 190
      LEND = LEN(CRDBUF)
C                               look for end of SET INPUT command
        DO 20 IC= 8,LEND
        IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 25
        IF (CRDBUF(IC:IC) .EQ. ',') GO TO 53
   20   CONTINUE
      GO TO 200
   25 CONTINUE
C         look for end of separator between command and first argument
      ICOL = IC+1
         DO 50 IC= ICOL,LEND
         IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 50
         IF (CRDBUF(IC:IC) .EQ. ',') GO TO 53
         GO TO 55
   50 CONTINUE
      GO TO 200
   53 IC = IC + 1
   55 IC1 = IC
C                      see if "REWIND" was requested in command
      LREWIN = .FALSE.
      IF (INDEX(CRDBUF(1:IC1),'REW') .GT. 5)  LREWIN=.TRUE.
      IF (INDEX(CRDBUF(1:IC1),'rew') .GT. 5)  LREWIN=.TRUE.
C                      first argument begins in or after col IC1
      DO 75 IC= IC1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 75
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 200
      GO TO 80
   75 CONTINUE
      GO TO 200
   80 IC1 = IC
C                        first argument really begins in col IC1
      DO 100 IC= IC1+1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 108
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 108
  100 CONTINUE
      IC = LEND + 1
  108 IC2 = IC-1
C                            end of first argument is in col IC2
  110 CONTINUE
      CUNIT = CRDBUF(IC1:IC2)
      WRITE (ISYSWR,'(A,A)') ' UNIT NO. :',CUNIT
      READ (CUNIT,'(BN,F10.0)',ERR=500) FUNIT
      IUNIT = FUNIT
      IF (IUNIT .EQ. 0)  GO TO 200
C                             skip blanks and commas, find file name
      DO 120 IC= IC2+1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 120
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 120
      GO TO 130
  120 CONTINUE
      GO TO 131
  130 CONTINUE
      CFNAME = CRDBUF(IC:LEND)
      NONAME = .FALSE.
      WRITE (ISYSWR, '(A,A)') ' FILE NAME IS:',CFNAME
C              ask if file exists, if not ask for name and open it
  131 CONTINUE
      INQUIRE(UNIT=IUNIT,OPENED=LOPEN,NAMED=LNAME,NAME=CGNAME)
      IF (LOPEN) THEN
         IF (NONAME) THEN
             GO TO 136
         ELSE
             IF (.NOT.LNAME) CGNAME='unknown'
             WRITE (ISYSWR,132) IUNIT,CGNAME,CFNAME
  132        FORMAT (' UNIT',I3,' ALREADY OPENED WITH NAME:',A/
     +                  '                 NEW NAME IGNORED:',A)
         ENDIF
      ELSE
C                new file, open it
         WRITE (ISYSWR,135) IUNIT
  135    FORMAT (' UNIT',I3,' IS NOT OPENED.')
         IF (NONAME) THEN
            WRITE (ISYSWR,'(A)') ' NO FILE NAME GIVEN IN COMMAND.'
            IF (ISW(6) .LT. 1)  GO TO 800
            WRITE (ISYSWR,'(A)') ' PLEASE GIVE FILE NAME:'
            READ (ISYSRD,'(A)') CFNAME
         ENDIF
         OPEN (UNIT=IUNIT,FILE=CFNAME,STATUS='OLD',ERR=600)
         WRITE (ISYSWR,'(A)') ' FILE OPENED SUCCESSFULLY.'
      ENDIF
C                                     . .   file is correctly opened
  136 IF (LREWIN) GO TO 150
      IF (ISW(6) .LT. 1)  GO TO 300
      WRITE (ISYSWR,137)  IUNIT
  137 FORMAT (' SHOULD UNIT',I3,' BE REWOUND?' )
      READ  (ISYSRD,'(A)')  CANSWR
      IF (CANSWR.NE.'Y' .AND. CANSWR.NE.'y') GO TO 300
  150 REWIND IUNIT
      GO TO 300
C                      *EOF
  190 CONTINUE
      IF (NSTKRD .EQ. 0)  THEN
         IERR = 2
         GO TO 900
         ENDIF
C                      revert to previous input file
  200 CONTINUE
      IF (NSTKRD .EQ. 0)  THEN
          WRITE (ISYSWR, '(A,A)') ' COMMAND IGNORED:',CRDBUF
          WRITE (ISYSWR, '(A)') ' ALREADY READING FROM PRIMARY INPUT'
      ELSE
        ISYSRD = ISTKRD(NSTKRD)
        NSTKRD = NSTKRD - 1
        IF (NSTKRD .EQ. 0)  ISW(6) = IABS(ISW(6))
        IF (ISW(5) .GE. 0)  THEN
          INQUIRE(UNIT=ISYSRD,NAMED=LNAME,NAME=CFNAME)
          CMODE = 'BATCH MODE      '
          IF (ISW(6) .EQ. 1)  CMODE = 'INTERACTIVE MODE'
          IF (.NOT.LNAME) CFNAME='unknown'
          IF (MNUNPT(CFNAME))  CFNAME='unprintable'
          WRITE (ISYSWR,290) CMODE,ISYSRD,CFNAME
  290     FORMAT (' INPUT WILL NOW BE READ IN ',A,' FROM UNIT NO.',I3/
     +    ' FILENAME: ',A)
        ENDIF
      ENDIF
      GO TO 900
C                      switch to new input file, add to stack
  300 CONTINUE
      IF (NSTKRD .GE. MAXSTK)  THEN
          WRITE (ISYSWR, '(A)') ' INPUT FILE STACK SIZE EXCEEDED.'
          GO TO 800
          ENDIF
      NSTKRD = NSTKRD + 1
      ISTKRD(NSTKRD) = ISYSRD
      ISYSRD = IUNIT
C                   ISW(6) = 0 for batch, =1 for interactive, and
C                      =-1 for originally interactive temporarily batch
      IF (ISW(6) .EQ. 1)  ISW(6) = -1
      GO TO 900
C                      format error
  500 CONTINUE
      WRITE (ISYSWR,'(A,A)') ' CANNOT READ FOLLOWING AS INTEGER:',CUNIT
      GO TO 800
  600 CONTINUE
      WRITE (ISYSWR, 601) CFNAME
  601 FORMAT (' SYSTEM IS UNABLE TO OPEN FILE:',A)
C                      serious error
  800 CONTINUE
      IERR = 3
  900 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNTINY(EPSP1,EPSBAK)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Compares its argument with the value 1.0, and returns
CC        the value .TRUE. if they are equal.  To find EPSMAC
CC        safely by foiling the Fortran optimizer
CC
      PARAMETER (ONE=1.0)
      EPSBAK =  EPSP1  - ONE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      LOGICAL FUNCTION MNUNPT(CFNAME)
C           is .TRUE. if CFNAME contains unprintable characters.
      CHARACTER CFNAME*(*)
      CHARACTER CPT*80, CP1*40,CP2*40
      PARAMETER (CP1=' ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm')
      PARAMETER (CP2='nopqrstuvwxyz1234567890./;:[]$%*_!@#&+()')
      CPT=CP1//CP2
      MNUNPT = .FALSE.
      L = LEN(CFNAME)
      DO 100 I= 1, L
         DO 50 IC= 1, 80
         IF (CFNAME(I:I) .EQ. CPT(IC:IC))  GO TO 100
   50    CONTINUE
      MNUNPT = .TRUE.
      GO TO 150
  100 CONTINUE
  150 CONTINUE
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNVERS(CV)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC         Returns the Minuit version in CV, char*6
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER*(*) CV
      CV = CVRSN
      RETURN
      END
* Jan 2013 R.Kessler : replace arithmetic IF with standard GOTO
* -----------------------
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.2  1996/03/15 18:02:54  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNVERT(A,L,M,N,IFAIL)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        inverts a symmetric matrix.   matrix is first scaled to
CC        have all ones on the diagonal (equivalent to change of units)
CC        but no pivoting is done since matrix is positive-definite.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      DIMENSION A(L,M) ,PP(MNI), Q(MNI),  S(MNI)
      IFAIL=0
      IF (N .LT. 1)  GO TO 100
      IF (N .GT. MAXINT)  GO TO 100
C                   scale matrix by sqrt of diag elements
      DO 8  I=1,N
      SI = A(I,I)

C ------------------------------ replace arithmetic IF
c      IF (SI) 100,100,8
      IF( SI .LT. 0 ) GOTO 100
      IF( SI .EQ. 0 ) GOTO 100
      IF( SI .GT. 0 ) GOTO 8
C ------------------------------


    8 S(I) = 1.0/SQRT(SI)
      DO 20 I= 1, N
      DO 20 J= 1, N
   20 A(I,J) = A(I,J) *S(I)*S(J)
C                                        . . . start main loop . . . .
      DO 65 I=1,N
      K = I
C                   preparation for elimination step1
      IF (A(K,K) .EQ. ZERO)  GO TO 100
      Q(K)=1./A(K,K)
      PP(K) = 1.0
      A(K,K)=0.0
      KP1=K+1
      KM1=K-1

C ------------------------------ replace arithmetic IF
c      IF(KM1)100,50,40
      IF ( KM1 .LT. 0 ) GOTO 100
      IF ( KM1 .EQ. 0 ) GOTO  50
      IF ( KM1 .GT. 0 ) GOTO  40
C ------------------------------

   40 DO 49 J=1,KM1
      PP(J)=A(J,K)
      Q(J)=A(J,K)*Q(K)
   49 A(J,K)=0.

C ------------------------------ replace arithmetic IF
   50 CONTINUE
c     IF(K-N)51,60,100
      IF ( (K-N) .LT. 0 ) GOTO 51
      IF ( (K-N) .EQ. 0 ) GOTO 60
      IF ( (K-N) .GT. 0 ) GOTO 100
C --------------------------------------

   51 DO 59 J=KP1,N
      PP(J)=A(K,J)
      Q(J)=-A(K,J)*Q(K)
   59 A(K,J)=0.0
C                   elimination proper
   60 DO 65 J=1,N
      DO 65 K=J,N
   65 A(J,K)=A(J,K)+PP(J)*Q(K)
C                   elements of left diagonal and unscaling
      DO 70 J= 1, N
      DO 70 K= 1, J
      A(K,J) = A(K,J) *S(K)*S(J)
   70 A(J,K) = A(K,J)
      RETURN
C                   failure return
  100 IFAIL=1
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNWARN(COPT,CORG,CMES)
C     If COPT='W', CMES is a WARning message from CORG.
C     If COPT='D', CMES is a DEBug message from CORG.
C         If SET WARnings is in effect (the default), this routine
C             prints the warning message CMES coming from CORG.
C         If SET NOWarnings is in effect, the warning message is
C             stored in a circular buffer of length MAXMES.
C         If called with CORG=CMES='SHO', it prints the messages in
C             the circular buffer, FIFO, and empties the buffer.
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER COPT*1, CORG*(*), CMES*(*), CTYP*7
      PARAMETER (MAXMES=10)
      CHARACTER     ORIGIN(MAXMES,2)*10, WARMES(MAXMES,2)*60
      COMMON/MN7WRC/ORIGIN,              WARMES
      COMMON/MN7WRI/NFCWAR(MAXMES,2),ICIRC(2)
      CHARACTER ENGLSH*20
C
      IF (CORG(1:3).EQ.'SHO' .AND. CMES(1:3).EQ.'SHO')  GO TO 200
C             Either print warning or put in buffer
      IF (COPT .EQ. 'W')  THEN
        ITYP = 1
        IF (LWARN) THEN
          WRITE (ISYSWR,'(A,A/A,A)') ' MINUIT WARNING IN ',CORG,
     +              ' ============== ',CMES
          CALL FLUSH(6)
          RETURN
        ENDIF
      ELSE
        ITYP = 2
        IF (LREPOR) THEN
          WRITE (ISYSWR,'(A,A/A,A)') ' MINUIT DEBUG FOR  ',CORG,
     +              ' ============== ',CMES
          RETURN
        ENDIF
      ENDIF
C                 if appropriate flag is off, fill circular buffer
         IF (NWRMES(ITYP) .EQ. 0)  ICIRC(ITYP) = 0
         NWRMES(ITYP) = NWRMES(ITYP) + 1
         ICIRC(ITYP) = ICIRC(ITYP) + 1
         IF (ICIRC(ITYP) .GT. MAXMES) ICIRC(ITYP) = 1
         IC = ICIRC(ITYP)
         ORIGIN(IC,ITYP) = CORG
         WARMES(IC,ITYP) = CMES
         NFCWAR(IC,ITYP) = NFCN
      RETURN
C
C             'SHO WARnings', ask if any suppressed mess in buffer
  200 CONTINUE
      IF (COPT .EQ. 'W') THEN
        ITYP = 1
        CTYP = 'WARNING'
      ELSE
        ITYP = 2
        CTYP = '*DEBUG*'
      ENDIF
      IF (NWRMES(ITYP) .GT. 0) THEN
         ENGLSH = ' WAS SUPPRESSED.  '
         IF (NWRMES(ITYP) .GT. 1) ENGLSH = 'S WERE SUPPRESSED.'
         WRITE (ISYSWR,'(/1X,I5,A,A,A,A/)') NWRMES(ITYP),
     +    ' MINUIT ',CTYP,' MESSAGE', ENGLSH
         NM = NWRMES(ITYP)
         IC = 0
         IF (NM .GT. MAXMES) THEN
              WRITE (ISYSWR,'(A,I2,A)')  ' ONLY THE MOST RECENT ',
     +          MAXMES,' WILL BE LISTED BELOW.'
              NM = MAXMES
              IC = ICIRC(ITYP)
         ENDIF
         WRITE (ISYSWR,'(A)') '  CALLS  ORIGIN         MESSAGE'
           DO 300 I= 1, NM
           IC = IC + 1
           IF (IC .GT. MAXMES)  IC = 1
           WRITE (ISYSWR,'(1X,I6,1X,A,1X,A)')
     +           NFCWAR(IC,ITYP),ORIGIN(IC,ITYP),WARMES(IC,ITYP)
 300       CONTINUE
         NWRMES(ITYP) = 0
         WRITE (ISYSWR,'(1H )')
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNWERR
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC          Calculates the WERR, external parameter errors,
CC      and the global correlation coefficients, to be called
CC      whenever a new covariance matrix is available.
CC
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=110)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
C                         calculate external error if v exists
      IF (ISW(2) .GE. 1) THEN
      DO 100 L= 1, NPAR
        NDEX = L*(L+1)/2
        DX = SQRT(ABS(VHMAT(NDEX)*UP))
        I = NEXOFI(L)
        IF (NVARL(I) .GT. 1)  THEN
          AL = ALIM(I)
          BA = BLIM(I) - AL
          DU1 = AL + 0.5 *(SIN(X(L)+DX) +1.0) * BA - U(I)
          DU2 = AL + 0.5 *(SIN(X(L)-DX) +1.0) * BA - U(I)
          IF (DX .GT. 1.0)  DU1 = BA
          DX = 0.5 * (ABS(DU1) + ABS(DU2))
        ENDIF
        WERR(L) = DX
  100 CONTINUE
      ENDIF
C                          global correlation coefficients
      IF (ISW(2) .GE. 1) THEN
         DO 130 I= 1, NPAR
            GLOBCC(I) = 0.
            K1 = I*(I-1)/2
            DO 130 J= 1, I
               K = K1 + J
               P(I,J) = VHMAT(K)
  130          P(J,I) = P(I,J)
         CALL MNVERT(P,MAXINT,MAXINT,NPAR,IERR)
         IF (IERR .EQ. 0)   THEN
            DO 150 IIN= 1, NPAR
               NDIAG = IIN*(IIN+1)/2
               DENOM = P(IIN,IIN)*VHMAT(NDIAG)
               IF (DENOM.LE.ONE .AND. DENOM.GE.ZERO)  THEN
                   GLOBCC(IIN) = 0.
               ELSE
                   GLOBCC(IIN) = SQRT(1.0-1.0/DENOM)
               ENDIF
  150       CONTINUE
         ENDIF
      ENDIF
      RETURN
      END
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE STAND
*
* $Id: minuit.F,v 1.5 2017/09/11 22:25:21 rkessler Exp $
*
* $Log: minuit.F,v $
* Revision 1.5  2017/09/11 22:25:21  rkessler
* *** empty log message ***
*
* Revision 1.4  2017/07/02 21:43:23  rkessler
* *** empty log message ***
*
* Revision 1.3  2017/06/27 12:44:57  rkessler
* *** empty log message ***
*
* Revision 1.2  2015/03/16 03:40:18  rkessler
* *** empty log message ***
*
* Revision 1.1  2015/02/19 19:55:24  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        optional user-supplied subroutine is called whenever the
CC        command "standard" appears.
CC
      RETURN
      END
