AC_PREREQ(2.59)

AC_INIT(pecnv, [0.1.8], [https://github.com/molpopgen/pecnv])
AC_CONFIG_SRCDIR([src/process_readmappings.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([disable])

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

dnl Do C++11 threads require -pthread for linking?

AC_MSG_CHECKING(that C++11 threads compile without the -pthread option)
AC_TRY_LINK(
[
#include <thread>
#include <vector>
void foo(int id){}
],
[
std::vector<std::thread> vt(2);
vt[0]=std::thread(foo,0);
vt[1]=std::thread(foo,1);
for( unsigned i = 0 ; i < vt.size() ; ++i ) vt[i].join();
],
[AC_MSG_RESULT(pthread not needed)]
AM_CONDITIONAL([NEED_PTHREAD], test x0 = x1 )
,
[AC_MSG_RESULT(pthread needed)]
dnl AC_DEFINE([NEED_PTHREAD],[],"adding -pthread to compiler/linker flags")
AM_CONDITIONAL([NEED_PTHREAD], test x1 = x1 )
)

dnl htslib
AC_CHECK_HEADER(htslib/bgzf.h,FOUND_HTSLIB_HEADERS=1,[AC_MSG_ERROR([htslib/bgzf.h not found. htslib >= 1.1 required.  See http://htslib.org.  Exiting.])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_HEADERS = x1 )

dnl specific libsequence features
AC_CHECK_HEADER(Sequence/bamrecord.hpp,,[AC_MSG_ERROR([Sequence/bamrecord.hpp not found.  libsequence >= 1.8.4 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/bamreader.hpp,,[AC_MSG_ERROR([Sequence/bamreader.hpp not found.  libsequence >= 1.8.4 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/samflag.hpp,,[AC_MSG_ERROR([Sequence/samflag.hpp not found.  libsequence >= 1.8.4 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(Sequence/IOhelp.hpp,,[AC_MSG_ERROR([Sequence/IOhelp.hpp not found.  libsequence >= 1.8.4 is required.  Please obtain correct version from http://molpopgen.org/software Exiting.])])
AC_CHECK_HEADER(boost/program_options.hpp,,[AC_MSG_ERROR([boost/program_options.hpp not found.  Please install boost, either via your favorite package manager or from http://www.boost.org])])
AC_CHECK_HEADER(boost/iostreams/filter/gzip.hpp,,[AC_MSG_ERROR([boost/iostreams/filter/gzip.hpp not found.  Please install boost, either via your favorite package manager or from http://www.boost.org])])
AC_CHECK_HEADER(boost/iostreams/filtering_stream.hpp,,[AC_MSG_ERROR([boost/iostreams/filter/filtering_stream.hpp not found.  Please install boost, either via your favorite package manager or from http://www.boost.org])])
AC_CHECK_HEADER(boost/iostreams/device/file.hpp,,[AC_MSG_ERROR([boost/iostreams/device/file.hpp not found.  Please install boost, either via your favorite package manager or from http://www.boost.org])])
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib development files required])])


AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
	  AC_CHECK_LIB(sequence, main,,[echo "libsequence runtime library does not seem to be here";exit 1])	
	  AC_CHECK_LIB(boost_program_options, main,,[echo "boost_program_options runtime library does not seem to be here.  Please install.";exit 1])	
	  AC_CHECK_LIB(boost_iostreams, main,,[echo "boost_iostreams runtime library does not seem to be here.  Please install.";exit 1])	
AC_LANG_RESTORE
dnl AM_CONFIG_HEADER(config.h)

AC_CHECK_LIB(z, gzread,,[AC_MSG_ERROR([zlib runtime library does not seem to be here.  Exiting.])])
AC_CHECK_LIB(hts, bgzf_read,FOUND_HTSLIB_LIBRARY=1,[AC_MSG_ERROR([htslib runtime library does not seem to be here.  Exiting.])])
AM_CONDITIONAL([HAVE_HTSLIB],test x$FOUND_HTSLIB_LIBRARY = x1 )

AC_OUTPUT
