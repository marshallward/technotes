dnl This file is part of MOM6, the Modular Ocean Model version 6.
dnl See the LICENSE file for licensing information.
dnl SPDX-License-Identifier: Apache-2.0
dnl
dnl MOM6_FC_PROTECT_PARENS
dnl ----------------------
dnl Check if the Fortran compiler respects parentheses in floating-point
dnl expressions (i.e., does not reassociate (a + b) - b into a).
dnl
dnl This is required for the round-to-nearest-integer trick:
dnl   K = (x + round_bias) - round_bias
dnl where round_bias = 1.5 * 2^52.
dnl
dnl Sets:
dnl   mom6_cv_fc_protect_parens=yes/no/unknown
dnl   HAVE_FC_PROTECT_PARENS (AC_DEFINE) if supported
dnl
dnl Compiler flags that enable this:
dnl   Intel ifx/ifort: -assume protect_parens
dnl   GCC gfortran:    -fprotect-parens (default with -std=f2008+)
dnl   NVHPC nvfortran: unknown
dnl
AC_DEFUN([MOM6_FC_PROTECT_PARENS], [
  AC_CACHE_CHECK([whether $FC respects floating-point parentheses],
    [mom6_cv_fc_protect_parens],
    [AC_LANG_PUSH([Fortran])
     AC_RUN_IFELSE(
       [AC_LANG_PROGRAM([], [dnl
      real, parameter :: fracwidth = digits(1.) - 1
      real, parameter :: bias = 1.5 * 2**fracwidth
      if (f(1.23) == 1.0) stop 0
      stop 1
      contains
      function round(x) result(y)
        real, intent(in) :: x
        real :: y
        y = (x + bias) - bias
      end function])
       ],
       [mom6_cv_fc_protect_parens=yes],
       [mom6_cv_fc_protect_parens=no],
       [mom6_cv_fc_protect_parens=unknown])
     AC_LANG_POP([Fortran])])

  AS_IF([test "$mom6_cv_fc_protect_parens" = yes],
    [AC_DEFINE([HAVE_FC_PROTECT_PARENS], [1],
      [Define to 1 if Fortran compiler respects parentheses in FP expressions])])
])
