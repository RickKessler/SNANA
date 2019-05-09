c This merges together two or more histogram files, using Hbook's
c  HMERGE subroutine.  This combines Ntuples properly too.
c
c Include the HMERGE function and fix OPEN calls to allow for
c capital letters.
c
c Nov 2 2014: double  IQUEST(10 to 2*64000  for outFIle ...
c             still cannot create out file larger than 262 MB.
c
c Dcc 2014: IQUEST back to 64,000 and reduce LQPAW from 40M to 20M
c            (to hopefully get fewer merge failures)
c -------------------

      program merge_hbook

      implicit none

c HBOOK common block
      integer LQPAW
      parameter (LQPAW = 20 000 000)
      integer PAWARR
      common /PAWC/ PAWARR(LQPAW)

c Local variables
      INTEGER MXFILE, ERRCODE
      PARAMETER ( MXFILE = 130, ERRCODE=17 )
      integer nargs
      integer nfiles, ifile
      character*200 infiles(MXFILE), outfile

c Function type declarations
      integer Iargc
c      integer Lenocc

c------ Beginning of code ------

c Check the number of command-line arguments

      nargs = Iargc()

      if (nargs .le. 1) then
        write (6,'(A)') 'Usage:  ' //
     &      'merge_hbook <infile1> <infile2> [<infile3> ... ] <outfile>'
        CALL EXIT(ERRCODE)
      endif

      if (nargs .gt. MXFILE) then
        write (6,'(A)') 'Too many input files!  See MXFILE limit'
        CALL EXIT(ERRCODE)
      endif

c Get the input and output filenames

      nfiles = nargs - 1
      do ifile = 1, nfiles
        call Getarg( ifile, infiles(ifile) )
      enddo
      call Getarg( nargs, outfile )

c Initialize Hbook

      call Hlimit( LQPAW )

c Do the merge

      call Hmerge( nfiles, infiles, outfile )

      call EXIT(0)
      end

* ------------------------------------
*
* $Id: merge_hbook.F,v 1.6 2019/01/16 19:55:08 rkessler Exp $
*
* $Log: merge_hbook.F,v $
* Revision 1.6  2019/01/16 19:55:08  rkessler
* *** empty log message ***
*
* Revision 1.5  2017/01/17 02:18:45  rkessler
* *** empty log message ***
*
* Revision 1.4  2014/12/26 22:48:05  rkessler
* *** empty log message ***
*
* Revision 1.3  2014/11/29 02:54:39  rkessler
* *** empty log message ***
*
* Revision 1.2  2013/05/21 14:49:19  rkessler
* *** empty log message ***
*
* Revision 1.1  2013/05/03 00:57:17  rkessler
* *** empty log message ***
*
* Revision 1.2  2013/04/27 13:45:12  rkessler
* *** empty log message ***
*
* Revision 1.1  2013/04/26 23:06:36  rkessler
* *** empty log message ***
*
* Revision 1.1.1.1  1996/01/16 17:08:10  mclareni
* First import
*
*
*#include "hbook/pilot.h"
*CMZ :  4.22/07 28/06/94  16.53.03  by  Rene Brun
*-- Author :    Rene Brun, J.C. Marin and Fons Rademakers   26/05/94
      SUBROUTINE HMERGE(NFILES,FILES,FILOUT)
*.==========>
*.
*.  Subroutine to merge the NFILES HBOOK files with identical objects
*.  and directories into FILOUT.
*.
*..=========> (Rene Brun, J.C. Marin, A.A.Rademakers )
*
      INTEGER         IQUEST
      COMMON/QUEST/   IQUEST(100)
*
      INTEGER         NFILES
      CHARACTER*(*) FILES(*), FILOUT
*
      INTEGER         LRECL, ISTAT
      INTEGER         LENOCC
*
      EXTERNAL        HMERGE1, HMERGE2
*____________________________________________________________________
*
      IF (NFILES .LE. 0) THEN
         WRITE(6,'('' No INPUT file. Job aborted'')')
         GOTO 99
      ENDIF
*
*         Open 1st file to get record length
*
*         Auto-record length detection does not work without option X
*         on Alpha/OSF/VMS
*

      LRECL = 0
      CALL HROPEN(21,'FILE',FILES(1),'XP ',LRECL,ISTAT)

      IF (ISTAT .NE. 0) THEN
         PRINT *, 'Cannot open file ', FILES(1)(1:LENOCC(FILES(1)))
         GOTO 99
      ENDIF
      CALL HRENDC('FILE')
*
*         Create big file
*
      IQUEST(10) = 64000
      CALL HROPEN(22,'BIGF',FILOUT(1:LENOCC(FILOUT)),'QNP',LRECL,ISTAT)
      IF (ISTAT .NE. 0) THEN
         PRINT *, 'Cannot create file ', FILOUT(1:LENOCC(FILOUT))
         GOTO 99
      ENDIF
*
*         Open 1st file to read characteristics
*
      LRELC = 0
      CALL HROPEN(21,'FILE',FILES(1),'P',LRECL,ISTAT)
      IF (ISTAT .NE. 0) THEN
         PRINT *, 'Cannot open file ', FILES(1)(1:LENOCC(FILES(1)))
         GOTO 99
      ENDIF
      CALL RZSCAN('//FILE', HMERGE1)
      CALL HRENDC('FILE')
*
*         Loop on all files
*
      DO 50 IFILE = 1, NFILES
         CALL HROPEN(21,'FILE',FILES(IFILE),'P',LRECL,ISTAT)
         IF (ISTAT .NE. 0) THEN
            PRINT *, 'Cannot open file ',
     +                FILES(IFILE)(1:LENOCC(FILES(IFILE)))
            GOTO 99
         ENDIF
         PRINT *, 'Processing file ',
     +             FILES(IFILE)(1:LENOCC(FILES(IFILE)))
         CALL RZSCAN(' ', HMERGE2)
         CALL HRENDC('FILE')
 50          CONTINUE
*
*         Save into the big file and close file
*
      CALL HCDIR('//BIGF',' ')
      CALL HROUT(0,ICYCLE,'T')
      IF (IQUEST(1) .NE. 0) THEN
         PRINT *, 'Error saving merged file'
         GOTO 99
      ENDIF
      CALL HRENDC('BIGF')
      PRINT *, 'Merge completed successfully'
      IQUEST(1)=0
      GO TO 999
*
*         Abnormal end
*
 99       CONTINUE
      IQUEST(1)=1
*
 999     END
