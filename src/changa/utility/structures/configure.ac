#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([TipsyFile.cpp])
AC_CONFIG_HEADER([config.h])

AC_PROG_CXX
AC_CHECK_TOOL(AR, ar) 

# Checks for header files.
AC_CHECK_HEADERS([stdint.h endian.h machine/endian.h sys/bitypes.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([u_quad_t, u_int64_t])

# Checks for library functions.
AC_HEADER_STDC
AC_CHECK_FUNCS([floor memset sqrt])

AC_CHECK_HEADERS([rpc/types.h])
AC_CHECK_HEADER([tirpc/netconfig.h],[HAVE_TIRPC=yes])
AC_CHECK_HEADER([rpc/xdr.h],,[MISSING_XDR=yes],
[#ifdef HAVE_RPC_TYPES_H
#include <rpc/types.h>
#endif
])
if test "${MISSING_XDR}" = "yes" -a "${HAVE_TIRPC}" = yes; then
   echo "XDR missing but TIRPC installed"
   CPPFLAGS="$CPPFLAGS -I/usr/include/tirpc"
   AC_SEARCH_LIBS([xdr_float], [xdr tirpc], [MISSING_XDR=no], [MISSING_XDR=yes])
else
   AC_SEARCH_LIBS([xdr_float], [xdr],,[MISSING_XDR=yes])
fi
if test "${MISSING_XDR}" = "yes" ; then
   echo "XDR missing, attempting to use xdr with package"
   if test -d ../xdr ; then
   	ln -s ../xdr ./rpc
	(cd ../xdr; make)
   else
	echo "Try downloading the xdr package and place it in the parent directory."
	echo "And do a make in ../xdr"
	exit
   fi
   CPPFLAGS="$CPPFLAGS -I."
   LDFLAGS="$LDFLAGS -L../xdr"
fi

# 128 bit keys:
AC_ARG_ENABLE([bigkeys],
        [AS_HELP_STRING([--enable-bigkeys], [enable 128 bit hash keys])],
        [bigkeys=yes], [bigkeys=no])
if test "$bigkeys" = "yes" ; then
    echo "BIGKEYS selected"
    FLAG_BIGKEYS="-DBIGKEYS"
else
    FLAG_BIGKEYS=""
fi

AC_CHECK_LIB([xdr], [xdr_hyper])
AC_CHECK_FUNCS([xdr_hyper])

AC_SUBST([FLAG_BIGKEYS])
AC_SUBST([USER_LDFLAGS], $LDFLAGS)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
