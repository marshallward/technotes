dnl Determine the flag required to honor parentheses in floating-point
dnl expressions (i.e., prevent reassociation of (a + b) - b into a).
dnl
dnl This is required for the round-to-nearest-integer trick:
dnl   K = (x + round_bias) - round_bias
dnl where round_bias = 1.5 * 2^52.
dnl
dnl Compiler flags that enable this:
dnl   GCC gfortran:    -fprotect-parens (default with -std=f2008+)
dnl   Intel ifx/ifort: -assume protect_parens
dnl   NVHPC nvfortran: (appears to be default)
dnl
AC_DEFUN([MOM6_FC_PROTECT_PARENS], [
  AC_ARG_ENABLE([protect-parens], [
    AS_HELP_STRING([--disable-protect-parens], [do not enforce parentheses protection])
  ])
  AC_ARG_ENABLE([strict-protect-parens], [
    AS_HELP_STRING([--enable-strict-protect-parens],
      [abort if parentheses protection is unsupported]
    )
  ])
  AC_ARG_VAR([PROTECT_PARENS_FCFLAGS],
    [User-defined Fortran flag to enforce parentheses protection])
  AS_IF([test "x$enable_protect_parens" != xno], [
    AS_IF([test -n "$PROTECT_PARENS_FCFLAGS"], [
      AC_MSG_NOTICE([Using user-defined PROTECT_PARENS_FCFLAGS: $PROTECT_PARENS_FCFLAGS])
    ], [
      AC_CACHE_CHECK([for $FC option to honor parentheses],
        [mom6_cv_fc_protect_parens], [
          mom6_cv_fc_protect_parens="unsupported"
          ac_fc_pp_FCFLAGS_save=${FCFLAGS}
          AC_LANG_PUSH([Fortran])
          for ac_flag in none \
            -fprotect-parens \
            "-assume protect_parens"
          do
            test "$ac_flag" != none \
              && FCFLAGS="$ac_fc_pp_FCFLAGS_save $ac_flag"
            AC_RUN_IFELSE([
              AC_LANG_PROGRAM([], [
dnl ---
      real, parameter :: bias = 1.5 * 2.**(digits(1.) - 1)
      if (nint_fast(1.23) == 1.) stop 0
      stop 1
      contains
      function nint_fast(x) result(y)
        real, intent(in) :: x
        real :: y
        y = (x + bias) - bias
      end function
dnl ---
              ])
            ], [
              mom6_cv_fc_protect_parens="$ac_flag"
              break
            ])
          done
          AC_LANG_POP([Fortran])
          FCFLAGS=$ac_fc_pp_FCFLAGS_save
        ]
      )
      AS_CASE([$mom6_cv_fc_protect_parens],
        [none], [
          mom6_cv_fc_protect_parens="none needed"
          PROTECT_PARENS_FCFLAGS=""
          AC_DEFINE([PROTECT_PARENS], [1],
            [Define to 1 if parentheses are protected in FP expressions])
        ],
        [unsupported], [
          AS_IF([test "x$enable_strict_protect_parens" = xyes], [
            AC_MSG_ERROR(
              [No known flag found to protect parentheses; aborting]
            )
          ], [
            AC_MSG_WARN(
              [No known flag found to protect parentheses; using defaults]
            )
            PROTECT_PARENS_FCFLAGS=""
          ])
        ],
        [
          PROTECT_PARENS_FCFLAGS=$mom6_cv_fc_protect_parens
          AC_DEFINE([PROTECT_PARENS], [1],
            [Define to 1 if parentheses are protected in FP expressions])
        ]
      )
    ])
  ])
  AC_SUBST([PROTECT_PARENS_FCFLAGS])
])
