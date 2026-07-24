!use, intrinsic :: iso_fortran_env, only : int64, real64
!real(kind=real64), parameter :: bias = 1.5_real64 * 2_int64**52
real, parameter :: bias = 1.5 * 2**(digits(1.0) - 1)
if (f(1.23) == 1.) stop 0
stop 1
contains
function f(x) result(y)
  real, intent(in) :: x
  real :: y
  y = (x + bias) - bias
end function
end
