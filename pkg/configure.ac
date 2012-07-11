AC_INIT([GUTS], [0.2.8])
AC_PREREQ(2.67)
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

if test "${GXX}" = yes; then
	gxx_version=`${CXX} -v 2>&1 | grep "^.*g.. version" | sed -e 's/^.*g.. version *//'`
	case ${gxx_version} in
		1.*|2.*)
			AC_MSG_WARN([Only g++ version 3.0 or greater can be used with GUTS.])
			AC_MSG_ERROR([Please use a different compiler.])
		;;
		4.6.*|4.7.*|4.8.*)
			gxx_newer_than_45="-fpermissive"
		;;
	esac
fi

AC_DEFUN([AC_PROG_R], [AC_CHECK_PROG(R,R,yes)])
AC_PROG_R

# Next segment by Kurt Hornik
# adapted from RQuantLib
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
#CXX=`${R_HOME}/bin/R CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

## look for Rscript, but use the one found via R_HOME to allow for multiple installations
# Dirk
AC_DEFUN([AC_PROG_RSCRIPT], [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
AC_PROG_RSCRIPT
#AC_MSG_CHECKING([for Rscript])
if test x"${RSCRIPT}" == x"yes" ; then
     rcpp_ldflags=`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`
     AC_SUBST([LDFLAGS],["${LDFLAGS} $rcpp_ldflags"])
else
    echo "
  Your installation does not appear to have Rscript installed.

  Please make sure that you have a working and complete R installation.
"
    exit 1
fi

AC_CONFIG_MACRO_DIR([m4])
AX_BOOST_BASE([1.45.0])
AC_SUBST([CXXFLAGS],["${CXXFLAGS} ${BOOST_CPPFLAGS}"])
AC_SUBST([LDFLAGS],["${LDFLAGS} ${BOOST_LDFLAGS}"])


AC_OUTPUT(src/Makevars)

echo "Completed configuration and ready to build."
