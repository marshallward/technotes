use iso_fortran_env, only : real64, real128
use exp_repro

implicit none

integer, parameter :: realq = merge(real128, real64, real128 > 0)

integer, parameter :: niter = 20

real(kind=real64), parameter :: xmax = 0.5 * log(2.)
real(kind=realq), allocatable :: xq(:), rq(:)

real(kind=real64), allocatable :: x(:), r(:), re(:), r_cr(:), r_t(:)
real(kind=realq), parameter :: xqmax = 0.5 * log(2.)

integer :: count_rate, count_max, c1, c2
real :: clock_rate

integer :: i, j, n

n = 10000000
allocate(x(n), r(n), re(n), r_cr(n), r_t(n))
allocate(xq(n), rq(n))

x = [(-xmax + (i - 1) * (2. * xmax) / (n-1), i=1,n)]
xq = [(-xqmax + (i - 1) * (2. * xqmax) / (n-1), i=1,n)]

! Reference realq values
rq = exp(xq)

! Set up clock
call system_clock(count_rate=count_rate, count_max=count_max)
clock_rate = real(count_rate)

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
print '(a18,1x,g0)', "err:", maxval(abs(re - rq))

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
print '(a18,1x,g0)', "err:", maxval(abs(r_cr - rq))
print '(a18,1x,g0)', "r - exp():", maxval(abs(r_cr - re))

!****

! Elemental Taylor version
r_t = exp_taylor(x)
call system_clock(count=c1)
!$omp simd
do i = 1, niter
  r_t = exp_taylor(x)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_taylor() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r_t - rq))
print '(a18,1x,g0)', "r - exp():", maxval(abs(r_t - re))

!****

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

!****

! Internal 1d GPU loop
r = -1
call exp_1d_do_c(x, r)
call system_clock(count=c1)
do i = 1, niter
  call exp_1d_do_c(x, r)
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_1d_do_c() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))

!****

! External GPU exp()
r = -1
call system_clock(count=c1)
do i = 1, niter
  do concurrent (j=1:n)
    r(j) = exp_cr(x(j))
  end do
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_cr loop() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_cr))

!****

! Taylor GPU loop
r = -1
call system_clock(count=c1)
do i = 1, niter
  do concurrent (j=1:n)
    r(j) = exp_taylor(x(j))
  end do
end do
call system_clock(count=c2)

! Report results
print '(a18,1x,g0)', "exp_taylor loop() time:", (c2 - c1) / clock_rate / niter
print '(a18,1x,g0)', "err:", maxval(abs(r - rq))
print '(a18,1x,g0)', "r - r[cpu]", maxval(abs(r - r_t))

end
