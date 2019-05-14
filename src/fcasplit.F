*
* $Id: fcasplit.F,v 1.1 2014/04/18 02:28:08 rkessler Exp $
*
* $Log: fcasplit.F,v $
* Revision 1.1  2014/04/18 02:28:08  rkessler
* *** empty log message ***
*
* Revision 1.3  1999/11/30 17:37:25  cernsoft
* generate compiler option +DA1.1 for HP architecture
*
* Revision 1.2  1999/03/25 09:58:07  cernlib
* Use g77 and gcc as compilers for Linux
*
* Revision 1.1  1996/05/21 10:04:39  cernlib
* copy from kernlib, remove need for pilot.h
*
* Revision 1.1.1.1  1996/02/15 17:54:16  mclareni
*
* Stolen from kerngen/obsolete

#ifdef CERNLIB_HPUX
#define CERNLIB_QMHPX
#endif

#ifdef CERNLIB_SUN
#define CERNLIB_QMSUN
#endif

#ifdef CERNLIB_DECS
#define CERNLIB_QMVMI
#endif

#ifdef CERNLIB_QMVAOS
#ifdef CERNLIB_QMVMI
#undef CERNLIB_QMVMI
#endif
#define CERNLIB_QMVAO
#endif

#ifdef CERNLIB_IBMRT
#define CERNLIB_QMIRT
#endif

#ifdef CERNLIB_LINUX
#define CERNLIB_QMLNX
#endif

#ifdef CERNLIB_SGI
#define CERNLIB_QMSGI
#endif

#if defined(CERNLIB_VAXVMS)
#define CERNLIB_QMVAX
#endif




      PROGRAM FCASPLIT

C CERN PROGLIB#         FCASPLIT        .VERSION KERNFOR  4.40  940929
C ORIG. 01/10/88  JZ

C#if defined(CERNLIB_QF_APO)
C      INTEGER*2    JARGA2, NARGSA2, PGM_$GET_ARG
C#endif
C#if !defined(CERNLIB_QMAMX)
      PARAMETER   (NARADD = 0)
C#endif
      CHARACTER    CHIDF*8, CHIDC*8, CHIDA*8, CHIDD*4
      CHARACTER    CHIDX*7, CHIDY*6, CHOVER*(*)
      CHARACTER    CHPOF*(*), CHPOC*(*), CHPOA*(*)
      CHARACTER    CHEXFOR*(*), CHEXCC*(*), CHEXAS*(*)
      CHARACTER    CHWHF*(*), CHWHC*(*), CHWHA*(*)
      CHARACTER    CHEXOBJ*(*)

      PARAMETER   (CHIDF = 'CDECK  I' )
      PARAMETER   (CHIDC = '/*DECK I' )
      PARAMETER   (CHIDX =  'DECK  I' )
      PARAMETER   (CHIDY =   'DECK I' )
      PARAMETER   (CHIDD =         'D>, ' )

C#if (!defined(CERNLIB_QMAPO10))&&(defined(CERNLIB_QF_APO))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'OVERWRITE')
C
C      PARAMETER   (CHPOF =
C     +     '-bounds_violation -info 1 -indexl -opt 3 -cpu mathlib_sr10')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.ftn')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.bin')
C
C      PARAMETER   (CHWHF = '/com/ftn  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if (defined(CERNLIB_QMAPO10))&&(defined(CERNLIB_QF_APO))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'OVERWRITE')
C
C      PARAMETER   (CHPOF =
C     +     '-bounds_violation -info 1 -indexl -opt 2 -cpu a88k')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.ftn')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.bin')
C
C      PARAMETER   (CHWHF = '/com/ftn  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if (!defined(CERNLIB_QMAPO10))&&(defined(CERNLIB_QF_APO77))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'OVERWRITE')
C
C      PARAMETER   (CHPOF = '-c -O3 -A cpu,mathlib_sr10' //
C     +     ' -W0,-bounds_violation -W0,-info,1 -W0,-indexl')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if (defined(CERNLIB_QMAPO10))&&(defined(CERNLIB_QF_APO77))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'OVERWRITE')
C
C      PARAMETER   (CHPOF = '-c -O2 -A cpu,a88k' //
C     +     ' -W0,-bounds_violation -W0,-info,1 -W0,-indexl')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMALT)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -Og -w')
C      PARAMETER   (CHPOC = '-c -Og')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'fortran  ')
C      PARAMETER   (CHWHC = 'fxc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMAMX)
C      PARAMETER   (NARADD = -1)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -q -O')
C      PARAMETER   (CHPOC = '-c')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'ftn  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMCRY)
C      PARAMETER   (CHIDA = '*DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = ' ')
C      PARAMETER   (CHPOC = '-c')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'cft77 ')
C      PARAMETER   (CHWHC = 'scc ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMCV64)
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O -p8 -72')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'fc  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMCV32)
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O -fi -72')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'fc  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if (defined(CERNLIB_QMDOS))&&(defined(CERNLIB_QF_NDP))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -vm')
C      PARAMETER   (CHPOC = '-c -n2 -n3 -OLM')
C      PARAMETER   (CHPOA = '-c')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.obj')
C
C      PARAMETER   (CHWHF = 'mf486 ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = '386ASM ')
C#endif
C#if (defined(CERNLIB_QMDOS))&&(defined(CERNLIB_QF_F2C))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-w -Nn802 -Nx802 -NC802')
C      PARAMETER   (CHPOC = '-c ')
C      PARAMETER   (CHPOA = '-c')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f2c ')
C      PARAMETER   (CHWHC = 'gcc  ')
C      PARAMETER   (CHWHA = 'as ')
C#endif
C#if (defined(CERNLIB_QMDOS))&&(defined(CERNLIB_QF_DEC))
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-nolink')
C      PARAMETER   (CHPOC = '-c -DWIN32')
C      PARAMETER   (CHPOA = '-c')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.obj')
C
C      PARAMETER   (CHWHF = 'f77 ')
C      PARAMETER   (CHWHC = 'claxp ')
C      PARAMETER   (CHWHA = 'ASM ')
C#endif
C#if defined(CERNLIB_QMHPX)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O -w +DA1.1 +ppu')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMIBX)
C      PARAMETER   (CHIDA = '*DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -f''opt(3)'' -f''CL(8192)'' ')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'fvs  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMIRT)
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O -qextname -qcharlen=8192')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'xlf  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMLNX)
      CHARACTER    CHIDA2*8
      PARAMETER   (CHIDA2= ' #DECK I' )
      PARAMETER   (CHIDA = ';DECK  I' )
      PARAMETER   (CHOVER= 'UNKNOWN')

      PARAMETER   (CHPOF = '-c -O -fno-automatic')
      PARAMETER   (CHPOC = '-c -O2 -m486')
      PARAMETER   (CHPOA = ' ')

      PARAMETER   (CHEXFOR = '.f')
      PARAMETER   (CHEXCC  = '.c')
      PARAMETER   (CHEXAS  = '.s')
      PARAMETER   (CHEXOBJ = '.o')

      PARAMETER   (CHWHF = 'g77  ')
      PARAMETER   (CHWHC = 'gcc  ')
      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMNXT)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C
C      PARAMETER   (CHPOF = '-c -O')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMOS9)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-R -E -!it -A -h')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'xcc ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMSGI)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMSUN)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O -w66')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMTMO)
C      PARAMETER   (CHIDA = '|DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c')
C      PARAMETER   (CHPOC = '-c')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMVAX)
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '/nolist/check=noover')
C      PARAMETER   (CHPOC = ' ')
C
C      PARAMETER   (CHEXFOR = '.for')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXOBJ = '.obj')
C
C      PARAMETER   (CHWHF = 'fortran')
C      PARAMETER   (CHWHC = 'cc')
C#endif
C#if defined(CERNLIB_QMVAX)
C      PARAMETER   (CHPOA  = '/nolist')
C      PARAMETER   (CHEXAS = '.mar')
C      PARAMETER   (CHWHA  = 'macro')
C
C      CHARACTER    BEGRUT*9, ENDLINE*12
C      PARAMETER   (BEGRUT = ' .TITLE  ')
C      PARAMETER   (ENDLINE= '        .END')
C#endif
C#if defined(CERNLIB_QMVMI)
C      CHARACTER    CHIDA2*8
C      PARAMETER   (CHIDA2= ' #DECK I' )
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'f77  ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif
C#if defined(CERNLIB_QMVAO)
C      CHARACTER    CHIDA2*8
C      PARAMETER   (CHIDA2= ' #DECK I' )
C      PARAMETER   (CHIDA = ';DECK  I' )
C      PARAMETER   (CHOVER= 'UNKNOWN')
C
C      PARAMETER   (CHPOF = '-c -O')
C      PARAMETER   (CHPOC = '-c -O')
C      PARAMETER   (CHPOA = ' ')
C
C      PARAMETER   (CHEXFOR = '.f')
C      PARAMETER   (CHEXCC  = '.c')
C      PARAMETER   (CHEXAS  = '.s')
C      PARAMETER   (CHEXOBJ = '.o')
C
C      PARAMETER   (CHWHF = 'f77  ')
C      PARAMETER   (CHWHC = 'cc   ')
C      PARAMETER   (CHWHA = 'as  ')
C#endif

      PARAMETER   (MLMKLN=100)
      PARAMETER   (MXLENG=128, MXMKLN=64)
      CHARACTER    CHHOLD*(MXLENG)
      CHARACTER    CHOPT(7)*(MXLENG)
      DIMENSION    NCHOPT(7)
      EQUIVALENCE (NOPF,NCHOPT(1))
      EQUIVALENCE (NOPC,NCHOPT(2))
      EQUIVALENCE (NOPA,NCHOPT(3))

      CHARACTER    CHCMP(3)*(MXLENG)
      DIMENSION    NCHCMP(3)
      EQUIVALENCE (CHCMP(1),CHOPT(4)), (NCHCMP(1),NCHOPT(4))
      EQUIVALENCE (NNF,NCHCMP(1))
      EQUIVALENCE (NNC,NCHCMP(2))
      EQUIVALENCE (NNA,NCHCMP(3))

      CHARACTER    CHFIN*(MXLENG)
      EQUIVALENCE (CHFIN,CHOPT(7)), (NCHFIN,NCHOPT(7))

      CHARACTER    CHFSH*(MXLENG), CHFMK*(MXLENG)
      CHARACTER    CHMKLN*(MXMKLN+4)

      PARAMETER   (MXCBUF=100)
      DIMENSION    LXCBUF(MXCBUF)
      CHARACTER*80 CHCBUF(MXCBUF)

      CHARACTER    CHLINE*80, CHNAME*80, CHTEXT*511
      CHARACTER    CHUSE*1

      PARAMETER   (NKEYS=11)
      CHARACTER    CHKEYS(NKEYS)*4
      DATA CHKEYS  / '-noh', '-log'
     +,              '-fo ', '-co ', '-ao '
     +,              '-f  ', '-c  ', '-a  '
     +,              '+fo ', '+co ', '+ao ' /

      LUNPR = 6
#if defined(CERNLIB_QMVAX)
      OPEN (LUNPR,FILE='SYS$OUTPUT',RECL=512,STATUS='OLD')
#endif

      WRITE (LUNPR,9001)
 9001 FORMAT (' FCASPLIT executing.')

#if defined(CERNLIB_QMIBX)
      CALL ERRSET (212, 300,-1,1,1)
#endif

 9002 FORMAT (
     F ' FCASPLIT     [-f  nmft] [-c  nmcc] [-a  nmas]   [-noh] [-log]'
     F/'              [+fo incf] [+co incc] [+ao incs]'
     F/'              [-fo optf] [-co optc] [-ao opts]'
     F/'     f.e  [fca_n] [optf      [optc      [opts]]]'/
     F/' splits file  f.e  having a mixture of Fortran / C / assembler'
     F/' routines into separate files n.f or n.c or n.s, "n" being the'
     F/' name of each routine, creating at the same time a Shell script'
     F/' y.shfca  and a Make file  y.mkfca  either of which can be used'
     F/' to compile all routines individually.'/
     F/' Defaults are defined in fcasplit for the names by which the'
     F/' compilers are called; with the -f, -c, -a options they could'
     F/' be changed.'/
     F/' Defaults are also defined for the options with which they are'
     F/' called; with the -fo, -co, -ao options they can be re-defined;'
     F/' with the +fo, +co, +ao options they can be incremented.'/
     F/' To be backward compatible the options can also be specified by'
     F/' the positional parameters after the file-name.')
 9003 FORMAT (1X
     F/' If the  -noh  option is given (or if the first parameter after'
     F/' the file-name is "fca_n") the identifying header line of each'
     F/' routine is not written out. If the -log option is given'
     F/' the name of each routine is printed on standard output.')
 9004 FORMAT (1X
     F/' Each routine must start with an identifying line :'
     F/'  "', A,  'D>, "      in cols.  1-12  for Fortran'
     F/'  "', A,  'D>, "      in cols.  1-12  for C'
     F/'  "', A,  'D>, "      in cols.  1-12  for assembler'
     F/'   "DECK  ID>, "      in cols.  2-12  or'
     F/'    "DECK ID>, "      in cols.  3-12  for anything else'
     F/'              "name"  in cols. 13-40  gives the name'/
     F/'        In the last two cases, or if "name" contains an'
     F/'        extension,  the file created will be "name"'
     F/'        without extension .f, .c or .s added to it'
     F/'        and without an entry into the script.'/
     F/'        A trailing blank terminates the name,'
     F/'        symbol . followed by blank also terminates,'
     F/'        symbols  ; < # !   all terminate,'
     F/'        symbol   */        also  terminates.'
     F/1X)

C------            Acquire the parameters

#if !defined(CERNLIB_QF_APO)
      NARGS = IARGC() + NARADD
#endif
#if defined(CERNLIB_QF_APO)
      CALL PGM_$GET_ARGS (NARGSA2,LPOIARG)
      NARGS = NARGSA2 - 1
#endif

      CHOPT(1) = CHPOF
      CHOPT(2) = CHPOC
      CHOPT(3) = CHPOA

      CHCMP(1) = CHWHF
      CHCMP(2) = CHWHC
      CHCMP(3) = CHWHA
      CHFIN    = 'f.e'

      IFINFI = 0
      IFLNH  = 0
      IFLLOG = 0

C--       Get the keyed options, until the input file name

      JARG = 0
   11 JOPT = 0
   12 IF (JARG.GE.NARGS)           GO TO 21
      JARG = JARG + 1
#if (!defined(CERNLIB_QF_APO))&&(!defined(CERNLIB_QMHPX))
      CALL GETARG (JARG, CHHOLD)
#endif
#if defined(CERNLIB_QMHPX)
      N = IGETARG (JARG, CHHOLD, MXLENG)
#endif
#if defined(CERNLIB_QF_APO)
      JARGA2 = JARG
      N = PGM_$GET_ARG (JARGA2,CHHOLD,ISTAT)
      CHHOLD(N+1:) = ' '
#endif
#if defined(CERNLIB_QMVAX)
      CALL CUTOL(CHHOLD)
#endif

      IF (JOPT.EQ.0)               GO TO 13

C----     Store option values read

C--       positional trailing options
      IF (JOPT.LT.0)  THEN
          IF (JOPT.EQ.-1)  THEN
              IF (CHHOLD(1:5).EQ.'fca_n')  THEN
                  IFLNH = 7
                  GO TO 12
                ENDIF
            ENDIF
          CHOPT(-JOPT) = CHHOLD
          JOPT = JOPT - 1
          GO TO 12
        ENDIF

C--       keyed options
      IF (JOPT.LE.6)  THEN
          CHOPT(JOPT) = CHHOLD
        ELSE
#if !defined(CERNLIB_QMVAX)
          IF (JOPT.EQ.7)  CHOPT(1) = CHPOF // ' ' // CHHOLD
          IF (JOPT.EQ.8)  CHOPT(2) = CHPOC // ' ' // CHHOLD
          IF (JOPT.EQ.9)  CHOPT(3) = CHPOA // ' ' // CHHOLD
#endif
#if defined(CERNLIB_QMVAX)
          IF (JOPT.EQ.7)  CHOPT(1) = CHPOF // CHHOLD
          IF (JOPT.EQ.8)  CHOPT(2) = CHPOC // CHHOLD
          IF (JOPT.EQ.9)  CHOPT(3) = CHPOA // CHHOLD
#endif
        ENDIF
      GO TO 11

C----     Analyse key, if any

C-         -noh -log  -fo -co -ao   -f  -c  -a  +fo +co +ao
C-            1    2    3   4   5    6   7   8    9  10  11
C-                      1   2   3    4   5   6    7   8   9

   13 DO 14  JKEY=1,NKEYS
      IF (CHHOLD(1:4).EQ.CHKEYS(JKEY))  THEN
          IF (JKEY.EQ.1)  THEN
              IFLNH = 7
              GO TO 11
            ENDIF
          IF (JKEY.EQ.2)  THEN
              IFLLOG = 7
              GO TO 11
            ENDIF
          JOPT = JKEY - 2
          GO TO 12
        ENDIF
   14 CONTINUE

C----     File-name read

      CHFIN = CHHOLD
      IFINFI = 7

C----     Get the positional options if given

      JOPT = -1
      GO TO 12

C------            Parameters have all been read

C--       Length of the option strings = true length + 2 (blanks)

   21 DO  24  JJ=1,7
      DO  23  J=MXLENG,1,-1
      IF (CHOPT(JJ)(J:J).NE.' ')     GO TO 24
   23 CONTINUE
      J = 0
   24 NCHOPT(JJ) = J + 2
      NCHFIN     = NCHFIN - 2
#if defined(CERNLIB_QMVAX)
      NCHOPT(4)  = NCHOPT(4) - 2
      NCHOPT(5)  = NCHOPT(5) - 2
      NCHOPT(6)  = NCHOPT(6) - 2
#endif

C--       Derive the name of the script  yyy.shfca
C-        from the input file,  say  dir/yyy.ext

      N  = NCHFIN
      JA = 1
      JE = N
      J  = N
   26 IF (CHFIN(J:J).EQ.'/')  THEN
          JA = J + 1
          GO TO 27
#if defined(CERNLIB_QMDOS)
          ELSEIF (CHFIN(J:J).EQ.CHAR(92))  THEN
          JA = J + 1
          GO TO 27
#endif
        ELSEIF (CHFIN(J:J).EQ.'.')  THEN
          IF (JE.EQ.N)  JE = J - 1
        ENDIF
      J = J - 1
      IF (J.NE.0)                  GO TO 26
   27 CHFSH  = CHFIN(JA:JE) // '.shfca'
      CHFMK  = CHFIN(JA:JE) // '.mkfca'
      NCHFSH = JE+7 - JA

      IF  (IFINFI.EQ.0)  THEN
          WRITE (LUNPR,9002)
          WRITE (LUNPR,9003)
          WRITE (LUNPR,9004) CHIDF, CHIDC, CHIDA
        ENDIF

      WRITE (LUNPR,9026) CHFIN(1:NCHFIN), CHFSH(1:NCHFSH)
     +,                  CHFMK(1:NCHFSH)
      WRITE (LUNPR,9027) CHCMP(1)(1:NNF),CHOPT(1)(1:NOPF)
      WRITE (LUNPR,9028) CHCMP(2)(1:NNC),CHOPT(2)(1:NOPC)
      WRITE (LUNPR,9029) CHCMP(3)(1:NNA),CHOPT(3)(1:NOPA)
 9026 FORMAT (5X,'        Input file : ',A/
     F        5X,'      Shell script : ',A/
     F        5X,'         Make file : ',A)
 9027 FORMAT (5X,'   Fortran    name : ',A/
     F        5X,'   Fortran options : ',A)
 9028 FORMAT (5X,'        CC    name : ',A/
     F        5X,'        CC options : ',A)
 9029 FORMAT (5X,' Assembler    name : ',A/
     F        5X,' Assembler options : ',A)

C----        Stop if no file-name, help information printed

      IF (IFINFI.EQ.0)  THEN
          WRITE (LUNPR,9030)
          STOP
        ENDIF

 9030 FORMAT (' !!! No file-name given, no execution !!!')

C------------      Process the input file       ------------------

C--                Open input and .shfca file

      OPEN  (11, FILE=CHFIN(1:NCHFIN),STATUS='OLD')
      REWIND 11

      OPEN  (21, FILE=CHFSH(1:NCHFSH),STATUS='OLD',ERR=301)
      CLOSE (21, STATUS='DELETE')
  301 OPEN  (22, FILE=CHFMK(1:NCHFSH),STATUS='OLD',ERR=302)
      CLOSE (22, STATUS='DELETE')
#if defined(CERNLIB_QMVAX)
  302 OPEN  (21, FILE=CHFSH(1:NCHFSH),RECL=512,STATUS=CHOVER)
      OPEN  (22, FILE=CHFMK(1:NCHFSH),RECL=512,STATUS=CHOVER)
#endif
#if !defined(CERNLIB_QMVAX)
  302 OPEN  (21, FILE=CHFSH(1:NCHFSH),STATUS=CHOVER)
      OPEN  (22, FILE=CHFMK(1:NCHFSH),STATUS=CHOVER)
#endif

#if defined(CERNLIB_WINNT)
      write(22,8000) "# Some nmake macros for Win32 apps development"
#endif
#if (defined(CERNLIB_WINNT))&&(defined(CERNLIB_QF_DEC))
      write(22,8000) ".SUFFIXES : .f"
      write(22,8000) "!include <f77.mak>"
#endif
#if defined(CERNLIB_WINNT)
      write(22,8000) "!include <ntwin32.mak>"
#endif

      CHMKLN = 'ROUTINES = '
      NXMKLN = 11

      NTEXT  = -12
      NIGNOR = 0
      NROUT  = 0
      NLINES = 0
      NTCBUF = 0
      NXCBUF = 0

      NLMKLN = 0
      NPCSMK = 0
#if defined(CERNLIB_QMVAX)
C-        Write SHFCA header
      WRITE(21,8000) '$! Command file to compile fortran,c,macro source'
      WRITE(21,8000) '$ PROC_VER = F$ENVIRONMENT("VERIFY_PROCEDURE")'
      WRITE(21,8000) '$ IMAG_VER = F$ENVIRONMENT("VERIFY_IMAGE")'
      WRITE(21,8000) '$ SET VERIFY=(PROCEDURE,NOIMAGE)'
      JTYLA = 0
#endif

C-------           Read next line        -------------------------

   31 CONTINUE
#if (!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMVMI))&&(!defined(CERNLIB_QMVAO))
      READ (11,8000,END=83) CHLINE
      N = LEN(CHLINE)
#endif
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
      READ (11,8031,END=83) N,CHLINE(1:MIN(80,N+1))
 8031 FORMAT (Q,A)
#endif

      DO  33  NCHLN=N,1,-1
      IF (CHLINE(NCHLN:NCHLN).NE.' ')  GO TO 34
   33 CONTINUE
      NCHLN = 0

   34 IF (NCHLN.GE.13)  THEN
          IF (CHLINE(9:12).EQ.CHIDD)   GO TO 59
        ENDIF
   35 IF (NTEXT.LT.0)              GO TO 58
   36 IF (JTYPE.NE.1)              GO TO 41

C--       Check new line is a Fortran comment line

      IF (NCHLN.EQ.0)              GO TO 38
      IF (CHLINE(1:1).EQ.'C')      GO TO 38
      IF (CHLINE(1:1).EQ.'c')      GO TO 38
      IF (CHLINE(1:1).NE.'*')      GO TO 41

   38 IF (NXCBUF.EQ.MXCBUF)        GO TO 41
      NXCBUF = NXCBUF + 1
      IF (NCHLN.GT.0)  CHCBUF(NXCBUF)(1:NCHLN) = CHLINE(1:NCHLN)
      LXCBUF(NXCBUF) = NCHLN
      GO TO 31

C-------           Write next line       -------------------------

C--        NTEXT = <0  idle unheaded lines (not coming here)
C-                  0  normal running
C-                 >0  about to write the first line of new routine
C-                     if =1 : no entry to script and make file

C-        start a new routine only when the first true line is ready

   41 IF (NTEXT.EQ.0)              GO TO 46
      IF (NTEXT.NE.1)              GO TO 81

C--       Open the output file if first line ready

   44 OPEN  (27, FILE=CHNAME(1:NCUM),STATUS='OLD',ERR=441)
      CLOSE (27, STATUS='DELETE')
  441 OPEN  (27, FILE=CHNAME(1:NCUM),STATUS=CHOVER)

      NTEXT = 0
      NROUT = NROUT + 1
      IF (IFLLOG.NE.0)  WRITE (LUNPR,9044) NROUT,CHNAME(1:NCUM)
 9044 FORMAT (' make',I4,1X,A)

#if defined(CERNLIB_QMVAX)
C--       Create .TITLE line for assembler
      JTYLA = JTYPE
      IF (JTYPE.EQ.3)  WRITE (27,9045) BEGRUT,CHNAME(1:NNAM)
 9045 FORMAT (A,A)
#endif

C--       Transfer the comment lines from the buffer

   46 IF (NXCBUF.EQ.0)             GO TO 48
      DO  47  J=1,NXCBUF
      N = LXCBUF(J)
      IF (N.EQ.0)  THEN
          WRITE (27,8000)
        ELSE
          WRITE (27,8000) CHCBUF(J)(1:LXCBUF(J))
        ENDIF
   47 CONTINUE
      NLINES = NLINES + NXCBUF
      NXCBUF = 0

C--       Transfer the current line

   48 IF (NCHLN.EQ.0)  THEN
          WRITE (27,8000)
        ELSE
          WRITE (27,8000) CHLINE(1:NCHLN)
        ENDIF
      NLINES = NLINES + 1
      GO TO 31

C--                Ignore leading unheaded lines

   58 NIGNOR = NIGNOR + 1
      IF (NTEXT.EQ.-1)             GO TO 31
      WRITE (LUNPR,9058) CHLINE(1:NCHLN)
      NTEXT = NTEXT + 1
      IF (NTEXT.NE.-1)             GO TO 31
      WRITE (LUNPR,9059)
      GO TO 31
 9058 FORMAT (' ignored: ',A)
 9059 FORMAT (10X,'...')

C-------           Start new routine ?   -------------------------

   59 JTYNX = 0
      IF (CHLINE(1:8).EQ.CHIDF)    GO TO 61
      IF (CHLINE(1:8).EQ.CHIDC)    GO TO 62
      IF (CHLINE(1:8).EQ.CHIDA)    GO TO 63
      IF (CHLINE(2:8).EQ.CHIDX)    GO TO 60
      IF (CHLINE(3:8).EQ.CHIDY)    GO TO 60
C#if defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
C      IF (CHLINE(1:8).EQ.CHIDA2)   GO TO 63
C#endif
      GO TO 35

C--                Yes, start new routine

   63 JTYNX = 1
   62 JTYNX = JTYNX + 1
   61 JTYNX = JTYNX + 1
   60 JPUT = 0
      JDOT = 0
      JDOP = 0
      NUS  = MIN(NCHLN,40)

C--       skip leading blanks before the name

      JGO = 13
   65 IF (CHLINE(JGO:JGO).EQ.' ')  THEN
          JGO = JGO + 1
          GO TO 65
        ENDIF

      IF (JGO.GT.NUS)              GO TO 35

C--       Convert name to lower case and find termination

      DO 66  JTK=JGO,NUS
      CHUSE = CHLINE(JTK:JTK)
      JV    = ICHAR(CHUSE)
      IF (CHUSE.EQ.' ')            GO TO 67
      IF (CHUSE.EQ.'.')  THEN
          IF (CHLINE(JTK+1:JTK+1).EQ.' ')   GO TO 67
          JDOP = JDOT
          JDOT = JPUT+1
        ENDIF
      IF (CHUSE.EQ.';')            GO TO 67
      IF (CHUSE.EQ.'<')            GO TO 67
      IF (CHUSE.EQ.'#')            GO TO 67
      IF (CHUSE.EQ.'!')            GO TO 67
      IF (CHUSE.EQ.'*')  THEN
          IF (CHLINE(JTK+1:JTK+1).EQ.'/')   GO TO 67
        ENDIF
      IF (JV.LT.91)  THEN
          IF (JV.GE.65)  CHUSE = CHAR(JV+32)
        ENDIF
      JPUT = JPUT + 1
      CHNAME(JPUT:JPUT) = CHUSE
   66 CONTINUE

C--                Handle explicit extension

   67 JTYPE = JTYNX
      IF (JDOT.EQ.JPUT)  THEN
          JDOT = JDOP
          JPUT = JPUT - 1
        ENDIF

      NTEXT = 1
      NPUT  = JPUT
      NCUM  = NPUT
      NNAM  = NPUT
      IF (JDOT.NE.0)  THEN
          NNAM  = JDOT - 1
          JTYPE = -JTYPE
        ENDIF

C----              Handle normal case with standard extensions

      IF (JTYPE.LE.0)              GO TO 80

C--       Fortran

      IF (JTYPE.NE.1)              GO TO 74
      N  =  LEN(CHEXFOR)
      CHNAME(NPUT+1:NPUT+N) = CHEXFOR
      NCUM   =  NPUT + N
#if (!defined(CERNLIB_QF_APO))&&(!defined(CERNLIB_QMVAX))
      CHTEXT =  CHCMP(1)(1:NNF) // CHOPT(1)(1:NOPF) // CHNAME(1:NCUM)
      NTEXT  =  NNF + NOPF + NCUM
#endif
#if defined(CERNLIB_QF_APO)
      CHTEXT =  CHCMP(1)(1:NNF) // CHNAME(1:NCUM) // ' ' //
     +          CHOPT(1)(1:NOPF)
      NTEXT  =  NNF + NCUM + NOPF
#endif
#if defined(CERNLIB_QMVAX)
      CHTEXT ='$ '//CHCMP(1)(1:NNF)//CHOPT(1)(1:NOPF)//CHNAME(1:NCUM)
      NTEXT  =  NNF + NOPF + NCUM +2
#endif
      GO TO 80

C--       CC

   74 IF (JTYPE.NE.2)              GO TO 77
      N  =  LEN(CHEXCC)
      CHNAME(NPUT+1:NPUT+N) = CHEXCC
      NCUM   =  NPUT + N
#if !defined(CERNLIB_QMVAX)
      CHTEXT =  CHCMP(2)(1:NNC) // CHOPT(2)(1:NOPC) // CHNAME(1:NCUM)
      NTEXT  =  NNC + NOPC + NCUM
#endif
#if defined(CERNLIB_QMVAX)
      CHTEXT ='$ '//CHCMP(2)(1:NNC)//CHOPT(2)(1:NOPC)//CHNAME(1:NCUM)
      NTEXT  =  NNC + NOPC + NCUM + 2
#endif
      GO TO 80

C--       Assembler

   77 IF (JTYPE.NE.3)              GO TO 80
      N  =  LEN(CHEXAS)
      CHNAME(NPUT+1:NPUT+N) = CHEXAS
      NCUM   =  NPUT + N
#if (!defined(CERNLIB_QMSUN))&&(!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMVMI))&&(!defined(CERNLIB_QMVAO))
      CHTEXT =  CHCMP(3)(1:NNA) // CHOPT(3)(1:NOPA) // CHNAME(1:NCUM)
      NTEXT  =  NNA + NOPA + NCUM
#endif
#if defined(CERNLIB_QMSUN)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
      CHTEXT =  CHCMP(3)(1:NNA) // '-o ' // CHNAME(1:NPUT) // '.o  ' //
     +          CHOPT(3)(1:NOPA) // CHNAME(1:NCUM)
      NTEXT  =  NNA + NPUT+7 + NOPA + NCUM
#endif
#if defined(CERNLIB_QMVAX)
      CHTEXT ='$ '//CHCMP(3)(1:NNA)//CHOPT(3)(1:NOPA)//CHNAME(1:NCUM)
      NTEXT  =  NNA + NOPA + NCUM + 2
#endif

   80 CONTINUE
#if defined(CERNLIB_QMVAX)
C--       write .END line for assembler
      IF (JTYLA.EQ.3)  WRITE (27,8000) ENDLINE
#endif
      CLOSE (27)
      NTCBUF = NTCBUF + NXCBUF
      NXCBUF = 0
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
      IF (JTYPE.EQ.3)              GO TO 31
#endif
      IF (JTYPE.LT.1)              GO TO 31
      IF (IFLNH.NE.0)              GO TO 31
      GO TO 36

C------       Output of the new deck is just starting

C--       compilation command for this routine to the shell script

   81 WRITE (21,8000) CHTEXT(1:NTEXT)

C--       register the routine name for the Make file

      CHTEXT =  CHNAME(1:NPUT) // CHEXOBJ // ' '
      NTEXT  =  NPUT + LEN(CHEXOBJ) + 1
      IF (NXMKLN+NTEXT.GT.MXMKLN) THEN
          NLMKLN = NLMKLN + 1
          IF (NLMKLN.LT.MLMKLN) THEN
              NXMKLN = NXMKLN + 1
              CHMKLN(NXMKLN:NXMKLN) = CHAR(92)
              WRITE (22,8000) CHMKLN(1:NXMKLN)
              NXMKLN = 0
            ELSE
              WRITE (22,8000) CHMKLN(1:NXMKLN)
              NXMKLN = 0
              NLMKLN = 0
              WRITE (22,9042)
              NPCSMK = NPCSMK + 1
              WRITE (CHMKLN,9043) NPCSMK
              NXMKLN = 11
            ENDIF
        ENDIF

 9042 FORMAT('#')
 9043 FORMAT('ROUTINE',I1,' = ')

      CHMKLN(NXMKLN+1:NXMKLN+NTEXT) = CHTEXT(1:NTEXT)
      NXMKLN = NXMKLN + NTEXT
      GO TO 44

C-------           Done                  -------------------------

C----           Complete the Make file

   83 WRITE (22,8000) CHMKLN(1:NXMKLN-1)

C--       Fortran inference

      CHNAME = CHEXFOR // CHEXOBJ // ':'
      NPUT   = LEN(CHEXFOR) + LEN(CHEXOBJ) + 1
#if !defined(CERNLIB_QF_APO)
      CHTEXT = CHAR(9) // CHCMP(1)(1:NNF) // CHOPT(1)(1:NOPF) //
     +            '$*' // CHEXFOR
      NTEXT   = 3 + NNF + NOPF + LEN(CHEXFOR)
#endif
#if defined(CERNLIB_QF_APO)
      CHTEXT = CHAR(9) // CHCMP(1)(1:NNF) // '$*' // CHEXFOR //
     +             ' ' // CHOPT(1)(1:NOPF)
      NTEXT   = 2 + NNF + LEN(CHEXFOR) + NOPF
#endif
#if (defined(CERNLIB_QMDOS))&&(defined(CERNLIB_QF_F2C))
C
C-- For F2C Fortran to C conversion is followed by gcc compilation
C
      WRITE (22,8000) CHTEXT(1:NTEXT)
#endif
      WRITE (22,9083) CHNAME(1:NPUT), CHTEXT(1:NTEXT)

C--       CC inference

      CHNAME = CHEXCC // CHEXOBJ // ':'
      NPUT   = LEN(CHEXCC) + LEN(CHEXOBJ) + 1
      CHTEXT = CHAR(9) // CHCMP(2)(1:NNC) // CHOPT(2)(1:NOPC) //
     +            '$*' // CHEXCC
      NTEXT   = 3 + NNC + NOPC + LEN(CHEXCC)
#if defined(CERNLIB_QF_APO)
      CHTEXT(NTEXT+1:NTEXT+16) = '; mv $*.o $*.bin'
      NTEXT = NTEXT + 16
#endif
      WRITE (22,9083) CHNAME(1:NPUT), CHTEXT(1:NTEXT)

C--       Assembler inference

      CHNAME = CHEXAS // CHEXOBJ // ':'
      NPUT   = LEN(CHEXAS) + LEN(CHEXOBJ) + 1
#if (!defined(CERNLIB_QMSUN))&&(!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMVMI))&&(!defined(CERNLIB_QMVAO))
      CHTEXT = CHAR(9) // CHCMP(3)(1:NNA) // CHOPT(3)(1:NOPA) //
     +            '$*' // CHEXAS
      NTEXT   = 3 + NNA + NOPA + LEN(CHEXAS)
#endif
#if defined(CERNLIB_QMSUN)||defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
      CHTEXT = CHAR(9) // CHCMP(3)(1:NNA) // CHOPT(3)(1:NOPA) //
     +            '-o $*.o  $*' // CHEXAS
      NTEXT   = 12 + NNA + NOPA + LEN(CHEXAS)
#endif
      WRITE (22,9083) CHNAME(1:NPUT), CHTEXT(1:NTEXT)

      IF (NPCSMK.EQ.0) THEN
         WRITE (22,9084) CHFIN(JA:JE)
      ELSE
         WRITE (22,9085) CHFIN(JA:JE),(I,I=1,NPCSMK)
         WRITE (22,9086) (I,I,I=1,NPCSMK)
      ENDIF

C--                Print summary

      WRITE (LUNPR,9087) NLINES,NROUT
      IF (NIGNOR.NE.0)  WRITE (LUNPR,9088) NIGNOR
      IF (NTCBUF.NE.0)  WRITE (LUNPR,9089) NTCBUF
#if defined(CERNLIB_QMVAX)
      WRITE(21,8000) '$!Restore original user settings  '
      WRITE(21,8000) '$ TEMP=F$VERIFY(PROC_VER,IMAG_VER)'
#endif
#if defined(CERNLIB_QMAPO)
      CALL PGM_$EXIT
#endif
 8000 FORMAT (A)
 9083 FORMAT (/A/A)
 9084 FORMAT(/A,'_all: $(ROUTINES)',/)
 9085 FORMAT(/A,'_all:',10(' rout',I1))
 9086 FORMAT(/,'rout: $(ROUTINES)',/,
     +        ('rout',I1,': $(ROUTINE',I1,')',/))
 9087 FORMAT (1X,I6,' lines written for',I6,' decks')
 9088 FORMAT (1X,I6,' leading unheaded lines ignored.')
 9089 FORMAT (1X,I6,' trailing comment lines ignored.')
      END
