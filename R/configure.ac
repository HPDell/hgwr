AC_INIT([hgwrr], 0.2-2)

AC_LANG(C++)
AC_REQUIRE_CPP

AC_PATH_PROG([PKG_CONFIG], [pkg-config])
if test "${PKG_CONFIG}" != ""; then
    MONGO_CFLAGS=`${PKG_CONFIG} --cflags libmongocxx-static`
    MONGO_LIBS=`${PKG_CONFIG} --libs libmongocxx-static`
else
    AC_MSG_ERROR([pkg-config not found, is it installed?])
fi

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test "${GSL_CONFIG}" != ""; then
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_SUBST(MONGO_CFLAGS)
AC_SUBST(MONGO_LIBS)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
