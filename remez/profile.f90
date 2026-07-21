use, intrinsic :: iso_fortran_env, only : real64, real128
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
use exp_repro

implicit none

integer, parameter :: realq = merge(real128, real64, real128 > 0)

integer, parameter :: niter = 20

real(kind=real64), parameter :: xmax = 0.5 * log(2._real64)

real(kind=realq), allocatable :: xq(:), rq(:)
real(kind=real64), allocatable :: x(:), r(:), re(:), r_cr(:), r_t(:)

integer :: count_rate, count_max, c1, c2
real :: clock_rate

integer :: i, j, n
real(kind=real64) :: y, nan

n = 10000000
allocate(x(n), r(n), re(n), r_cr(n), r_t(n))
allocate(xq(n), rq(n))

! NOTE: Using real64 to build the real128 x-points.  I think this is right?
x = [(-xmax + (i - 1) * (2. * xmax) / (n-1), i=1,n)]
!x = [(-xmax + (i - 1) * (1. * xmax) / (n-1), i=1,n)]
xq = real(x, real128)

! Reference realq values
rq = exp(xq)

! Set up clock
call system_clock(count_rate=count_rate, count_max=count_max)
clock_rate = real(count_rate)

!****

! First verify some values

print '(a22,3(1x,ES25.17E3))', "exp(0):", &
    exp(0._real128), exp(0._real64), exp_cr(0._real64)
print '(a22,3(1x,ES25.17E3))', "exp(1):", &
    exp(1._real128), exp(1._real64), exp_cr(1._real64)
print '(a22,3(1x,ES25.17E3))', "exp(0.33):", &
!print '(a22,1x,Z32.32,2(1x,Z16.16))', "exp(0.33):", &
    exp(0.33_real128), exp(0.33_real64), exp_cr(0.33_real64)

y = .3225414126648429_real64
print '(a22,3(1x,ES25.17E3))', "exp(.322541412664843):", &
  exp(real(y, real128)), exp(y), exp_cr(y)

! Extreme values

print '(a22,3(1x,ES25.17E3))', "overflow: exp(1000):", &
    exp(1000._real128), exp(1000._real64), exp_cr(1000._real64)
print '(a22,3(1x,ES25.17E3))', "underflow exp(-1000):", &
    exp(-1000._real128), exp(-1000._real64), exp_cr(-1000._real64)

! Special values

y = log(huge(1._real64))
print '(a22,3(1x,ES25.17E3))', "largest float:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = log(tiny(1.0_real64))
print '(a22,3(1x,ES25.17E3))', "smallest float:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = -1060.0_real64 * log(2.0_real64)
print '(a22,3(1x,ES25.17E3))', "median subnormal:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
!y = log(2.0_real64**(-1074))
y = -1074.0_real64 * log(2.0_real64)
print '(a22,3(1x,ES25.17E3))', "smallest subnormal:", &
  exp(real(y, real128)), exp(y), exp_cr(y)
y = -1075.0_real64 * log(2.0_real64)
print '(a22,3(1x,ES25.17E3))', "subnormal zero-cutoff", &
  exp(real(y, real128)), exp(y), exp_cr(y)

y = ieee_value(nan, ieee_quiet_nan)
print '(a22,3(1x,ES25.17E3))', "NaN:", &
  exp(ieee_value(y, ieee_quiet_nan)), exp(y), exp_cr(y)

!****

! Intrinsic exp()
re = exp(x)
call system_clock(count=c1)
do i = 1, niter
  re = exp(x)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0,1x,"x=",g0)', "err:", &
    maxval(abs(re - rq)), x(maxloc(abs(re - rq)))
print '(a18,1x,g0,1x,"x=",g0)', "rel err:", &
    maxval(abs((re - rq) / rq)), x(maxloc(abs((re - rq) / rq)))

!****

! Elemental Remez+Estrin version
r_cr = exp_cr(x)
call system_clock(count=c1)
!$omp simd
do i = 1, niter
  r_cr = exp_cr(x)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_cr() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0,1x,"x=",g0)', "err:", &
  maxval(abs(r_cr - rq)), x(maxloc(abs(r_cr - rq)))
print '(a18,1x,g0,1x,"x=",g0)', "rel err:", &
  maxval(abs((r_cr - rq) / rq)), x(maxloc(abs((r_cr - rq) / rq)))
print '(a18,1x,g0)', "r - exp():", maxval(abs((r_cr - re) / re))

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
!*!print '(a18,1x,g0)', "exp_1d() time:", (c2 - c1) / clock_rate / niter
!*!print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
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
!*!print '(a18,1x,g0)', "exp_1d_do_c() time:", (c2 - c1) / clock_rate / niter
!*!print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))
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
!*!print '(a18,1x,g0)', "exp_cr loop() time:", (c2 - c1) / clock_rate / niter
!*!print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))
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
!*!print '(a18,1x,g0)', "exp_taylor loop() time:", (c2 - c1) / clock_rate / niter
!*!print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
!*!print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_t))

end
