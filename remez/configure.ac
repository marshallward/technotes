dnl configure.ac for exp_repro test/demo
dnl
dnl Run:
dnl   autoreconf -i
dnl   ./configure
dnl   make

AC_INIT([exp_repro_demo], [1.0])
AC_CONFIG_SRCDIR([profile.f90])
AC_CONFIG_MACRO_DIRS([m4])

dnl Find Fortran compiler
AC_PROG_FC

dnl Check for protect_parens support
MOM6_FC_PROTECT_PARENS
AS_IF([test -n "$PROTECT_PARENS_FCFLAGS"],
  [FCFLAGS="$FCFLAGS $PROTECT_PARENS_FCFLAGS"]
)

dnl Output
AC_SUBST([FCFLAGS])
AC_SUBST([DEFS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
