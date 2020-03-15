AC_INIT([mpl], [1.0], [moller@mollerware.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_PROG_CXX_C_O
AC_PROG_LN_S
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE([enable])

PKG_CHECK_MODULES([GSL], [gsl >= 1.15])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

#PKG_CHECK_MODULES([LIBFFI], [libffi >= 3.0.13], ,
#  [echo "";
#  echo ">>>> You need to install the libffi and libffi-devel packages. <<<<";
#  echo ""; exit])
#AC_SUBST(LIBFFI_CFLAGS)
#AC_SUBST(LIBFFI_LIBS)

#PKG_CHECK_MODULES([PLPLOT], [plplot >= 5.11.1], ,
#  [echo "";
#  echo ">>>> You need to install the plplot and plplot-devel packages. <<<<";
#  echo ""; exit])
#AC_SUBST(PLPLOT_CFLAGS)
#AC_SUBST(PLPLOT_LIBS)

#AC_CHECK_PROG([CHECK_DISPLAY], [display], [true])
#if test x"$CHECK_DISPLAY" != x"true" ; then
#    AC_MSG_WARN([Imagemagick display not found.  Plot visualisation will be disabled.])
#else
#    AC_SUBST(HAVE_DISPLAY, "-DHAVE_DISPLAY")
#fi


AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
