## Set R_HOME, respecting an environment variable if one is set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=""

## Set the fallback, by default it is nope
arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"

## now use all these
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
