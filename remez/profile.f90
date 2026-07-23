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

use exp_repro

implicit none

integer, parameter :: realq = merge(real128, real64, real128 > 0)

integer, parameter :: niter = 20

real(kind=real64), parameter :: xmax = 0.5 * log(2._real64)

real(kind=realq), allocatable :: xq(:), rq(:)
real(kind=real64), allocatable :: x(:), r(:), re(:), r_cr(:), r_fast(:), r_t(:)

integer :: count_rate, count_max, c1, c2
real :: clock_rate

integer :: i, j, n
real(kind=real64) :: y, z

! testing
logical :: invalid_ref, overflow_ref, underflow_ref, inexact_ref, divzero_ref
logical :: invalid_cr, overflow_cr, underflow_cr, inexact_cr, divzero_cr

n = 10000000
allocate(x(n), r(n), re(n), r_cr(n), r_fast(n), r_t(n))
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
    exp(0._real128), exp(0._real64), exp_cr(0._real64)
print '(a26,3(1x,ES25.17E3))', "exp(-0):", &
    exp(-0._real128), exp(-0._real64), exp_cr(-0._real64)
print '(a26,3(1x,ES25.17E3))', "exp(1):", &
    exp(1._real128), exp(1._real64), exp_cr(1._real64)
print '(a26,3(1x,ES25.17E3))', "exp(0.33):", &
!print '(a26,1x,Z32.32,2(1x,Z16.16))', "exp(0.33):", &
    exp(0.33_real128), exp(0.33_real64), exp_cr(0.33_real64)

y = .3225414126648429_real64
print '(a26,3(1x,ES25.17E3))', "exp(.322541412664843):", &
  exp(real(y, real128)), exp(y), exp_cr(y)

! Extreme values

print '(a26,3(1x,ES25.17E3))', "overflow: exp(1000):", &
    exp(1000._real128), exp(1000._real64), exp_cr(1000._real64)
print '(a26,3(1x,ES25.17E3))', "underflow exp(-1000):", &
    exp(-1000._real128), exp(-1000._real64), exp_cr(-1000._real64)

! Special values

y = log(huge(1._real64))
print '(a26,3(1x,ES25.17E3))', "largest float:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = log(tiny(1.0_real64))
print '(a26,3(1x,ES25.17E3))', "smallest float:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = -1060.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "median subnormal:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
!y = log(2.0_real64**(-1074))
y = -1074.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "smallest subnormal:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = -1075.0_real64 * log(2.0_real64)
print '(a26,3(1x,ES25.17E3))', "subnormal zero-cutoff", &
  exp(real(y, real128)), exp(y), exp_cr(y)

y = ieee_value(y, ieee_positive_inf)
print '(a26,3(1x,ES25.17E3))', "+Inf:", &
  exp(ieee_value(0._real128, ieee_positive_inf)), exp(y), exp_cr(y)
y = ieee_value(y, ieee_negative_inf)
print '(a26,3(1x,ES25.17E3))', "-Inf:", &
  exp(ieee_value(0._real128, ieee_negative_inf)), exp(y), exp_cr(y)
y = ieee_value(y, ieee_quiet_nan)
print '(a26,3(1x,ES25.17E3))', "NaN:", &
  exp(ieee_value(0._real128, ieee_quiet_nan)), exp(y), exp_cr(y)

!***

! Exception test
! TODO: Function?

!y = ieee_value(y, ieee_positive_inf)
y = -1000._real64

call ieee_set_flag(ieee_all, .false.)
z = exp(y)

call ieee_get_flag(ieee_invalid, invalid_ref)
call ieee_get_flag(ieee_overflow, overflow_ref)
call ieee_get_flag(ieee_underflow, underflow_ref)
call ieee_get_flag(ieee_inexact, inexact_ref)
call ieee_get_flag(ieee_divide_by_zero, divzero_ref)

call ieee_set_flag(ieee_all, .false.)
z = exp_cr(y)

call ieee_get_flag(ieee_invalid, invalid_cr)
call ieee_get_flag(ieee_overflow, overflow_ref)
call ieee_get_flag(ieee_underflow, underflow_ref)
call ieee_get_flag(ieee_inexact, inexact_ref)
call ieee_get_flag(ieee_divide_by_zero, divzero_ref)

print *, "invalid: ", invalid_ref, invalid_cr
print *, "overflow: ", overflow_ref, overflow_cr
print *, "underflow: ", underflow_ref, underflow_cr
print *, "inexact: ", inexact_ref, inexact_cr
print *, "divzero: ", divzero_ref, divzero_cr

call ieee_set_flag(ieee_all, .false.)

!****

! Intrinsic exp()
! Use the same explicit inner loop shape as exp_cr() so the timing compares
! function cost rather than array-assignment lowering choices.
!$omp simd
do j = 1, size(x)
  re(j) = exp(x(j))
enddo
call system_clock(count=c1)
do i = 1, niter
  !$omp simd
  do j = 1, size(x)
    re(j) = exp(x(j))
  enddo
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
!$omp simd
do j = 1, size(x)
  r_cr(j) = exp_cr(x(j))
enddo
call system_clock(count=c1)
do i = 1, niter
  !$omp simd
  do j = 1, size(x)
    r_cr(j) = exp_cr(x(j))
  enddo
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp_cr() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp_cr() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_cr - rq)), x(maxloc(abs(r_cr - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_cr - rq) / rq)), x(maxloc(abs((r_cr - rq) / rq)))
print '(a26,1x,g0)', "r - exp():", maxval(abs((r_cr - re) / re))

!****

! Fast finite-normal Remez+Estrin version
!$omp simd
do j = 1, size(x)
  r_fast(j) = exp_cr_fast(x(j))
enddo
call system_clock(count=c1)
do i = 1, niter
  !$omp simd
  do j = 1, size(x)
    r_fast(j) = exp_cr_fast(x(j))
  enddo
end do
call system_clock(count=c2)

! Report results
print '(a26,1x,g0)', "exp_cr_fast() time:", (c2 - c1) / clock_rate / niter
print '(a26,1x,g0)', "exp_cr_fast() time/elem:", (c2 - c1) / clock_rate / niter / n
print '(a26,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_fast - rq)), x(maxloc(abs(r_fast - rq)))
print '(a26,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_fast - rq) / rq)), x(maxloc(abs((r_fast - rq) / rq)))
print '(a26,1x,g0)', "r - exp_cr():", maxval(abs((r_fast - r_cr) / r_cr))

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
!*!    r(j) = exp_cr(x(j))
!*!  end do
!*!end do
!*!call system_clock(count=c2)
!*!
!*!! Report results
!*!print '(a26,1x,g0)', "exp_cr loop() time:", (c2 - c1) / clock_rate / niter
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

end
