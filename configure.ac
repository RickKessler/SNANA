
AC_INIT([snana], [1.0], [rkessler@uchicago.edu])
AC_PREREQ([2.68]) 
AC_CONFIG_SRCDIR([src/sntools.c])
AC_CONFIG_AUX_DIR([build-aux])
#AM_INIT_AUTOMAKE([foreign -Wall -Werror])
dnl Tell automake that this is not following GNU coding standards (no README, CHANGELOG, etc)
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIR([m4])

dnl setting up gfortran and gcc
AC_PROG_CC([gcc])
AC_PROG_FC([gfortran])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [GSLCBLASLIB="-lgslcblas"], [], [$LIBS])
AC_SUBST([GSLCBLASLIB])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [GSLLIB="-lgsl"], [], [$LIBS])
AC_SUBST([GSLLIB])
AC_CHECK_LIB([cfitsio], [ffopen], [CFITSIOLIB="-lcfitsio"], [AC_MSG_ERROR([cfitsio library not found. Please install CFITSIO.])], [$LIBS])
AC_SUBST([CFITSIOLIB])

AC_ARG_VAR(ROOTSYS)
if test -n "$ROOTSYS"; then
  root_avail=yes
  AM_CONDITIONAL([USE_ROOT], [test "x$root_avail" = "xyes"])
  AC_SUBST(USE_ROOT) 
  echo "ROOTSYS is set to: $ROOTSYS"
fi

dnl set up C++ 17 flags
m4_include([m4/m4_ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX(17)


dnl need to include config.sug and config.guess with the distribution
AC_CANONICAL_HOST


dnl define git version so that C code can write full git version+tag to outputs
AC_PATH_PROG([GIT], [git])

if test -n "$GIT"; then
  GIT_VERSION=`$GIT --version | sed 's/git version //'`
  AC_MSG_NOTICE([Found Git version $GIT_VERSION])
else
  AC_MSG_ERROR([Git not found. Please install Git.])
fi




dnl check for MAC OS
macos=no

dnl case "${host_os}" in
dnl darwin*)
dnl macos=yes
dnl ;;
dnl esac

case "$host_os" in
  darwin*) macos=yes ;;
  *)       macos=no ;;
esac


dnl set OSX to yes for macos
AM_CONDITIONAL([OSX], [test "x$macos" = "xyes"])
AC_MSG_NOTICE([macos $macos])

dnl hard wired to 64 bit machine ... if somebody needs 32 bit, they can change hard-wire manually
have64=yes
havei686=no


AM_CONDITIONAL([BIT64], [test "x$have64" = "xyes"])
AM_CONDITIONAL([I686], [test "x$havei686" = "xyes"])

AC_OUTPUT
