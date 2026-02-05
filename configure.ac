
AC_INIT([snana], [1.0], [rkessler@uchicago.edu])
AC_PREREQ([2.68]) 
AC_CONFIG_SRCDIR([src/sntools.c])
AC_CONFIG_AUX_DIR([build-aux])
#AM_INIT_AUTOMAKE([foreign -Wall -Werror])
dnl Tell automake that this is not following GNU coding standards (no README, CHANGELOG, etc)
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIR([m4])

UTIL_DIR='$(top_srcdir)/util'
AC_SUBST([UTIL_DIR])

dnl setting up gfortran and gcc
AC_PROG_CC([gcc])
AC_PROG_F77([gfortran])
AC_PROG_FC([gfortran])

AC_FC_LIBRARY_LDFLAGS

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])


# Try to get Git version
AC_CHECK_PROG([GIT], [git], [git], [no])
if test "$GIT" != "no"; then
  AC_MSG_CHECKING([Git version])
  GIT_VERSION=`$GIT describe --always --tags 2>/dev/null || echo "unknown"`
  AC_MSG_RESULT([$GIT_VERSION])
else
  GIT_VERSION="unknown"
fi
AC_SUBST([GIT_VERSION])


AM_PATH_PYTHON([3.10])

AC_PATH_PROGS([PYTHON_CONFIG], [python3-config python-config], [no])
AS_IF([test "x$PYTHON_CONFIG" = xno], [AC_MSG_ERROR([cannot find python-config])])

PYTHON_INCLUDES=`$PYTHON_CONFIG --includes`
AC_SUBST([PYTHON_INCLUDES])

if test "x$PYTHON_CONFIG" != "xno"; then
	# We seem to require the python3.x libray, so we use embed
  PYTHON_LIBS=`$PYTHON_CONFIG --embed --ldflags`
#  AC_DEFINE([USE_PYTHON], [1], [Define to 1 if Python is available.])
else
  AC_MSG_ERROR([python3-config not found])
fi
AC_SUBST([PYTHON_LIBS])

#AM_CONDITIONAL([USE_PYTHON], [test "x$PYTHON_CONFIG" = "xyes"])
#AC_SUBST(USE_PYTHON)

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
else
  root_avail=no
fi

AM_CONDITIONAL([USE_ROOT], [test "x$root_avail" = "xyes"])
echo "ROOTSYS is ${ROOTSYS:-not set}"

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
