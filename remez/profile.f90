use, intrinsic :: iso_fortran_env, only : real64, real128
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf

use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_get_flag
use, intrinsic :: ieee_exceptions, only : ieee_all
use, intrinsic :: ieee_exceptions, only : ieee_invalid
use, intrinsic :: ieee_exceptions, only : ieee_overflow
use, intrinsic :: ieee_exceptions, only : ieee_underflow
use, intrinsic :: ieee_exceptions, only : ieee_inexact
use, intrinsic :: ieee_exceptions, only : ieee_divide_by_zero

use exp_dev
use exp_cr_mod, only : exp_repro_new => exp_repro

implicit none

integer, parameter :: realq = merge(real128, real64, real128 > 0)

integer, parameter :: niter = 20

!real(kind=real64), parameter :: xmax = 0.5 * log(2._real64)
real(kind=real64), parameter :: xmax = 10._real64
!real(kind=real64), parameter :: xmax = 700._real64

real(kind=realq), allocatable :: xq(:), rq(:)
real(kind=real64), allocatable :: x(:), r(:), re(:), r_cr(:), r_fast(:), r_new(:)

integer :: count_rate, count_max, c1, c2
real :: clock_rate

integer :: i, j, n
real(kind=real64) :: y

n = 10000000
allocate(x(n), r(n), re(n), r_cr(n), r_fast(n), r_new(n))
allocate(xq(n), rq(n))

! Input range

! Symmetric test
x = [(-xmax + (i - 1) * (2. * xmax) / (n-1), i=1,n)]

! Just negative numbers
!x = [(-xmax + (i - 1) * (1. * xmax) / (n-1), i=1,n)]

! NOTE: Using real64 to build the real128 x-points.  I think this is right?
xq = real(x, real128)

! Reference realq values
rq = exp(xq)

! Set up clock
call system_clock(count_rate=count_rate, count_max=count_max)
clock_rate = real(count_rate)

!****

! First verify some values

print '(a26,3(1x,ES25.17E3))', "exp(0):", &
    exp(0._real128), exp(0._real64), exp_dev_cr(0._real64)
print '(a26,3(1x,ES25.17E3))', "exp(-0):", &
    exp(-0._real128), exp(-0._real64), exp_dev_cr(-0._real64)
print '(a26,3(1x,ES25.17E3))', "exp(1):", &
    exp(1._real128), exp(1._real64), exp_dev_cr(1._real64)
print '(a26,3(1x,ES25.17E3))', "exp(-1):", &
    exp(-1._real128), exp(-1._real64), exp_dev_cr(-1._real64)
print '(a26,3(1x,ES25.17E3))', "exp(0.33):", &
!print '(a26,1x,Z32.32,2(1x,Z16.16))', "exp(0.33):", &
    exp(0.33_real128), exp(0.33_real64), exp_dev_cr(0.33_real64)

y = .3225414126648429_real64
print '(a26,3(1x,ES25.17E3))', "exp(.322541412664843):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)

! Range reduction edge cases
y = log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "exp(ln(2)):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "exp(-ln(2)):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = 0.5_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "exp(0.5*ln(2)):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)

! Tiny arguments (tests 1 + x accuracy)
y = 1.0e-15_real64
print '(a26,3(1x,ES25.17E3))', "exp(1e-15):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -1.0e-15_real64
print '(a26,3(1x,ES25.17E3))', "exp(-1e-15):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)

! Round-trip (accumulates log + exp error)
y = log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "exp(log(2)):", &
  exp(log(2.0_real128)), exp(y), exp_dev_cr(y)

! Extreme values

y = 1000._real64
print '(a26,3(1x,ES25.17E3))', "overflow: exp(1000):", &
    exp(1000._real128), exp(y), exp_dev_cr(y)
y = -1000._real64
print '(a26,3(1x,ES25.17E3))', "underflow exp(-1000):", &
    exp(-1000._real128), exp(y), exp_dev_cr(y)

! Near overflow/underflow boundaries
y = 709.78_real64
print '(a26,3(1x,ES25.17E3))', "near overflow (709.78):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -708.39_real64
print '(a26,3(1x,ES25.17E3))', "near underflow (-708.39):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -745.13_real64
print '(a26,3(1x,ES25.17E3))', "subnormal region (-745.13):", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)

! Special values

y = log(huge(1._real64))
print '(a26,3(1x,ES25.17E3))', "largest float:", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = log(tiny(1.0_real64))
print '(a26,3(1x,ES25.17E3))', "smallest float:", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -1060.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "median subnormal:", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
!y = log(2.0_real64**(-1074))
y = -1074.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "smallest subnormal:", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)
y = -1075.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "subnormal zero-cutoff", &
  exp(real(y, real128)), exp(y), exp_dev_cr(y)

y = ieee_value(y, ieee_positive_inf)
print '(a26,3(1x,ES25.17E3))', "+Inf:", &
  exp(ieee_value(0._real128, ieee_positive_inf)), exp(y), exp_dev_cr(y)
y = ieee_value(y, ieee_negative_inf)
print '(a26,3(1x,ES25.17E3))', "-Inf:", &
  exp(ieee_value(0._real128, ieee_negative_inf)), exp(y), exp_dev_cr(y)
y = ieee_value(y, ieee_quiet_nan)
print '(a26,3(1x,ES25.17E3))', "NaN:", &
  exp(ieee_value(0._real128, ieee_quiet_nan)), exp(y), exp_dev_cr(y)
y = -ieee_value(y, ieee_quiet_nan)
print '(a26,3(1x,ES25.17E3))', "-NaN:", &
  exp(-ieee_value(0._real128, ieee_quiet_nan)), exp(y), exp_dev_cr(y)

!***

! Exception flag reports.  Columns are: invalid, overflow, underflow,
! inexact, divide-by-zero.
print '(a26,3(1x,a6))', "IEEE flags:", "exp128", "exp()", "exp_dev"
print '(a26,3(1x,a6))', "flag order:", "IOUXZ", "IOUXZ", "IOUXZ"
call print_exception_flags("normal:", 0.33_real64)
call print_exception_flags("overflow:", 1000.0_real64)
call print_exception_flags("underflow:", -1000.0_real64)
call print_exception_flags("largest float:", log(huge(1.0_real64)))
call print_exception_flags("smallest normal:", log(tiny(1.0_real64)))
call print_exception_flags("+Inf:", ieee_value(0.0_real64, ieee_positive_inf))
call print_exception_flags("-Inf:", ieee_value(0.0_real64, ieee_negative_inf))
call print_exception_flags("NaN:", ieee_value(0.0_real64, ieee_quiet_nan))
call print_exception_flags("-NaN:", -ieee_value(0.0_real64, ieee_quiet_nan))

call ieee_set_flag(ieee_all, .false.)

!****

! Intrinsic exp()
! Use the same explicit inner loop shape as exp_dev_cr() so the timing compares
! function cost rather than array-assignment lowering choices.
do i = 1,3
  re = exp(x)
enddo
call system_clock(count=c1)
do i = 1, niter
  re = exp(x)
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
    maxval(abs(re - rq)), x(maxloc(abs(re - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
    maxval(abs((re - rq) / rq)), x(maxloc(abs((re - rq) / rq)))

!****

! Elemental Remez+Estrin version
do i = 1,3
  r_cr = exp_dev_cr(x)
end do
call system_clock(count=c1)
do i = 1, niter
  r_cr = exp_dev_cr(x)
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp_dev_cr() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp_dev_cr() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_cr - rq)), x(maxloc(abs(r_cr - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_cr - rq) / rq)), x(maxloc(abs((r_cr - rq) / rq)))
print '(a26,1x,g0)', "r - exp():", maxval(abs((r_cr - re) / re))

!****

! Fast finite-normal Remez+Estrin version
do j = 1, size(x)
  r_fast(j) = exp_dev_fast(x(j))
enddo
call system_clock(count=c1)
do i = 1, niter
  do j = 1, size(x)
    r_fast(j) = exp_dev_fast(x(j))
  enddo
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp_dev_fast() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp_dev_fast() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_fast - rq)), x(maxloc(abs(r_fast - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_fast - rq) / rq)), x(maxloc(abs((r_fast - rq) / rq)))
print '(a26,1x,g0)', "r - exp_dev_cr():", maxval(abs((r_fast - r_cr) / r_cr))

!****

! New standalone exp_repro from exp_dev.f90
do i = 1,3
  r_new = exp_repro_new(x)
end do
call system_clock(count=c1)
do i = 1, niter
  r_new = exp_repro_new(x)
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp_repro() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp_repro() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_new - rq)), x(maxloc(abs(r_new - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_new - rq) / rq)), x(maxloc(abs((r_new - rq) / rq)))
print '(a26,1x,g0)', "r - exp_dev_cr():", maxval(abs((r_new - r_cr) / r_cr))

!*!!****
!*!
!*!! Internal 1d loop
!*!call exp_1d(x, r)
!*!call system_clock(count=c1)
!*!do i = 1, niter
!*!  call exp_1d(x, r)
!*!end do
!*!call system_clock(count=c2)
!*!
!*!! Report results
!*!print '(a26,1x,g0)', "exp_1d() time:", (c2 - c1) / clock_rate / niter
!*!print '(a26,1x,g0)', "err:", maxval(abs(r - rq))
!*!
!*!!****
!*!
!*!! Internal 1d GPU loop
!*!r = -1
!*!call exp_1d_do_c(x, r)
!*!call system_clock(count=c1)
!*!do i = 1, niter
!*!  call exp_1d_do_c(x, r)
!*!end do
!*!call system_clock(count=c2)
!*!
!*!! Report results
!*!print '(a26,1x,g0)', "exp_1d_do_c() time:", (c2 - c1) / clock_rate / niter
!*!print '(a26,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a26,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))
!*!
!*!!****
!*!
!*!! External GPU exp()
!*!r = -1
!*!call system_clock(count=c1)
!*!do i = 1, niter
!*!  do concurrent (j=1:n)
!*!    r(j) = exp_dev_cr(x(j))
!*!  end do
!*!end do
!*!call system_clock(count=c2)
!*!
!*!! Report results
!*!print '(a26,1x,g0)', "exp_dev loop() time:", (c2 - c1) / clock_rate / niter
!*!print '(a26,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a26,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))
!*!
!*!!****
!*!
!*!! Taylor GPU loop
!*!r = -1
!*!call system_clock(count=c1)
!*!do i = 1, niter
!*!  do concurrent (j=1:n)
!*!    r(j) = exp_taylor(x(j))
!*!  end do
!*!end do
!*!call system_clock(count=c2)
!*!
!*!! Report results
!*!print '(a26,1x,g0)', "exp_taylor loop() time:", (c2 - c1) / clock_rate / niter
!*!print '(a26,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a26,1x,g0)', "r - r[cpu]", maxval(abs(r - r_t))

contains

subroutine print_exception_flags(label, x0)
  character(len=*), intent(in) :: label
  real(kind=real64), intent(in) :: x0

  real(kind=real128), volatile :: zq
  real(kind=real64), volatile :: z
  logical :: invalid_q, overflow_q, underflow_q, inexact_q, divzero_q
  logical :: invalid_ref, overflow_ref, underflow_ref, inexact_ref, divzero_ref
  logical :: invalid_cr, overflow_cr, underflow_cr, inexact_cr, divzero_cr

  call ieee_set_flag(ieee_all, .false.)
  zq = exp(real(x0, real128))
  call ieee_get_flag(ieee_invalid, invalid_q)
  call ieee_get_flag(ieee_overflow, overflow_q)
  call ieee_get_flag(ieee_underflow, underflow_q)
  call ieee_get_flag(ieee_inexact, inexact_q)
  call ieee_get_flag(ieee_divide_by_zero, divzero_q)

  call ieee_set_flag(ieee_all, .false.)
  z = exp(x0)
  call ieee_get_flag(ieee_invalid, invalid_ref)
  call ieee_get_flag(ieee_overflow, overflow_ref)
  call ieee_get_flag(ieee_underflow, underflow_ref)
  call ieee_get_flag(ieee_inexact, inexact_ref)
  call ieee_get_flag(ieee_divide_by_zero, divzero_ref)

  call ieee_set_flag(ieee_all, .false.)
  z = exp_dev_cr(x0)
  call ieee_get_flag(ieee_invalid, invalid_cr)
  call ieee_get_flag(ieee_overflow, overflow_cr)
  call ieee_get_flag(ieee_underflow, underflow_cr)
  call ieee_get_flag(ieee_inexact, inexact_cr)
  call ieee_get_flag(ieee_divide_by_zero, divzero_cr)

  print '(a26,3(1x,a5))', label, &
      flag_char(invalid_q, 'I') // flag_char(overflow_q, 'O') &
        // flag_char(underflow_q, 'U') // flag_char(inexact_q, 'X') &
        // flag_char(divzero_q, 'Z'), &
      flag_char(invalid_ref, 'I') // flag_char(overflow_ref, 'O') &
        // flag_char(underflow_ref, 'U') // flag_char(inexact_ref, 'X') &
        // flag_char(divzero_ref, 'Z'), &
      flag_char(invalid_cr, 'I') // flag_char(overflow_cr, 'O') &
        // flag_char(underflow_cr, 'U') // flag_char(inexact_cr, 'X') &
        // flag_char(divzero_cr, 'Z')
end subroutine print_exception_flags

pure function flag_char(flag, c) result(s)
  logical, intent(in) :: flag
  character(len=1), intent(in) :: c
  character(len=1) :: s
  s = merge(c, '.', flag)
end function flag_char

end
