! exp_repro_macros.h - Preprocessor macros for reproducible exp()
!
! These macros provide fallbacks for when compiler flags like
! -assume protect_parens (Intel) are not available.
!
! Usage:
!   #include "exp_repro_macros.h"
!
! Required: round_bias must be defined in scope as:
!   real(kind=WP), parameter :: round_bias = 1.5_WP * 2**digits(1.0_WP)
!
! Build flags:
!   -DPROTECT_PARENS  : Use fast bit-trick rounding (requires compiler
!                       to respect parentheses, e.g. -assume protect_parens)
!   (default)         : Use anint() fallback (slower but always correct)

#ifdef PROTECT_PARENS

! Fast round-to-nearest-integer using bit trick
#define NEAREST_INT(x) anint_fast(x)

#else

! Safe fallback using intrinsic (slower but always correct)
#define NEAREST_INT(x) anint(x)

#endif
