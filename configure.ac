AC_PREREQ(2.59)
M
AC_INIT(pecnv, [0.1.4], [https://github.com/molpopgen/pecnv])
AC_CONFIG_SRCDIR([src/process_readmappings.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config-h.in])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([enable])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile scripts/Makefile])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AC_ARG_ENABLE(profiling,
[  --enable-profiling    Turn on profiling],
[case "${enableval}" in
  yes) profiling=true ;;
  no)  profiling=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
esac],[profiling=false])
AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)

dnl check for things that we need

dnl htslib
AC_CHECK_HEADER(htslib/bgzf.h,FOUND_HTSLIB_HEADERS=1,[AC_MSG_ERROR([htslib/bgzf.h not found. htslib >= 1.1 required.  See http://htslib.org.  Exiting.])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_HEADERS = x1 )

dnl specific libsequence features
dnl AC_CHECK_HEADER(Sequence/samrecord.hpp,,[AC_MSG_ERROR([Sequence/samrecord.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
dnl AC_CHECK_HEADER(Sequence/samflag.hpp,,[AC_MSG_ERROR([Sequence/samflag.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/bamrecord.hpp,,[AC_MSG_ERROR([Sequence/bamrecord.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/bamreader.hpp,,[AC_MSG_ERROR([Sequence/bamreader.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/samflag.hpp,,[AC_MSG_ERROR([Sequence/samflag.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/IOhelp.hpp,,[AC_MSG_ERROR([Sequence/IOhelp.hpp not found.  libsequence >= 1.8.3 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])


AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib development files required])])


AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
	  AC_CHECK_LIB(sequence, main,,[echo "libsequence runtime library does not seem to be here";exit 1])	
AC_LANG_RESTORE
dnl AM_CONFIG_HEADER(config.h)

AC_CHECK_LIB(z, gzread,,[AC_MSG_ERROR([zlib runtime library does not seem to be here.  Exiting.])])
AC_CHECK_LIB(hts, bgzf_read,FOUND_HTSLIB_LIBRARY=1,[AC_MSG_ERROR([htslib runtime library does not seem to be here.  Exiting.])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_LIBRARY = x1 )

AC_OUTPUT
