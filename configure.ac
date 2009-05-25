dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl ACLOCAL="aclocal -I m4" autoreconf
dnl
dnl You must have recent versions of Autoconf and Automake installed.

# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $


AC_PREREQ(2.61)
AC_INIT([Specfem3D Globe], [4.0.0], [jtromp AT caltech.edu], [Specfem3DGlobe])
AC_CONFIG_SRCDIR([specfem3D.f90])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])


# 'configure' options

AC_ARG_ENABLE([double-precision],
    [AC_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test "$want_double_precision" = no; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])


# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"
AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])

flags_guess="$SHELL $srcdir/flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_PROG_CC



# influential environment variables

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags for non-critical subroutines])
AC_ARG_VAR(FLAGS_NO_CHECK, [Fortran compiler flags for creating fast, production-run code for critical subroutines])
if test x"$MPIFC" = x; then
    MPIFC=mpif90
fi

AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])
if test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x; then
    LOCAL_PATH_IS_ALSO_GLOBAL=true
fi



# Checks for libraries.


# Checks for header files.

CIT_FC_MPI_HEADER([$MPIFC], [$FLAGS_NO_CHECK])


# Checks for typedefs, structures, and compiler characteristics.


# Checks for library functions.


AC_CONFIG_FILES([Makefile constants.h precision.h])
AC_OUTPUT


dnl end of configure.ac