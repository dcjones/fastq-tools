
AC_INIT( [fastq-tools], [0.3], [dcjones@cs.washington.edu] )
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_C_O

m4_ifdef([AC_TYPE_UINT8_T],  [AC_TYPE_UINT8_T])
m4_ifdef([AC_TYPE_UINT16_T], [AC_TYPE_UINT16_T])
m4_ifdef([AC_TYPE_INT32_T],  [AC_TYPE_INT32_T])
m4_ifdef([AC_TYPE_UINT32_T], [AC_TYPE_UINT32_T])
m4_ifdef([AC_TYPE_UINT64_T], [AC_TYPE_UINT64_T])
m4_ifdef([AC_TYPE_SIZE_T],   [AC_TYPE_SIZE_T])
m4_ifdef([AC_HEADER_STDBOOL], [AC_HEADER_STDBOOL])

opt_CFLAGS="--std=c99 -Wall -Wextra -pedantic -g -O3 -D_C99_SOURCE"
dbg_CFLAGS="--std=c99 -Wall -Wextra -pedantic -g -O0 -D_C99_SOURCE"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

AC_DEFINE(_FILE_OFFSET_BITS, 64)
AC_DEFINE(_POSIX_SOURCE) # needed for fileno

AC_PROG_LIBTOOL

AC_CHECK_FUNC(fileno, ,
              AC_MSG_ERROR([The 'fileno' function is missing.]))

# check pcre
AC_CHECK_PROG(HAVE_PCRE, pcre-config, yes, no)
if test "x$HAVE_PCRE" = "xno"
then
    AC_MSG_ERROR([The PCRE library is needed. See http://www.pcre.org.])
fi
CFLAGS="$CFLAGS $(pcre-config --cflags)"
PCRE_LIBS="$(pcre-config --libs)"
AC_SUBST(PCRE_LIBS)

# check getopt
AC_CHECK_HEADER(getopt.h, ,
                AC_MSG_ERROR([The posix getopt.h header is needed.]))


CXXFLAGS=$CFLAGS

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile
                 doc/Makefile
                 src/version.h])

AC_OUTPUT
