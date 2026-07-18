use iso_fortran_env, only : real64, real128
use exp_repro

implicit none

integer, parameter :: niter = 20

real(kind=real64), parameter :: xmax = 0.5 * log(2.)
real(kind=real128), allocatable :: xq(:), rq(:)

real(kind=real64), allocatable :: x(:), r(:)
real(kind=real128), parameter :: xqmax = 0.5 * log(2.)

integer :: count_rate, count_max, c1, c2
real :: clock_rate

integer :: i, n

n = 10000000
allocate(x(n), r(n))
allocate(xq(n), rq(n))

x = [(-xmax + (i - 1) * (2. * xmax) / (n-1), i=1,n)]
xq = [(-xqmax + (i - 1) * (2. * xqmax) / (n-1), i=1,n)]

! Reference real128 values
rq = exp(xq)

! Set up clock
call system_clock(count_rate=count_rate, count_max=count_max)
clock_rate = real(count_rate)

! Intrinsic exp()
r = exp(x)
call system_clock(count=c1)
do i = 1, niter
  r = exp(x)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))

! Elemental version
r = exp_cr(x)
call system_clock(count=c1)
!$omp simd
do i = 1, niter
  r = exp_cr(x)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_cr() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))

! Internal 1d loop
call exp_1d(x, r)
call system_clock(count=c1)
do i = 1, niter
  call exp_1d(x, r)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_1d() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))

end
