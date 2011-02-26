
AC_INIT( [fastq-tools], [0.1], [dcjones@cs.washington.ed] )
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_C_O

opt_CFLAGS="-Wall -g -O3"
dbg_CFLAGS="-Wall -g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

AC_DEFINE(_FILE_OFFSET_BITS, 64)

AC_PROG_LIBTOOL

# check zlib
AX_CHECK_ZLIB

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

AC_CONFIG_FILES( [Makefile
                  src/Makefile
                  src/swsse2/Makefile] )

AC_OUTPUT
