module exp_repro

use iso_fortran_env, only : real64, real128

implicit none

!*!print *, exp(1000._16)
!*!!print *, exp(1000.)
!*!print *, exp_cr(1000.)
!*!
!*!print *, exp(0.)
!*!print *, exp_cr(0.)

contains

elemental function exp_cr(x) result(a)
  real(kind=real64), intent(in) :: x
    !< Input value [nondim]
  real(kind=real64) :: a
    !< Exponential of x [nondim]

  real(kind=real64) :: r
    ! Rescaled value of x: r = x - K * ln 2
  integer :: K
    ! Scaling factor, where exp(x) = 2^K exp(r)
  real(kind=real64) :: e
    ! Scaled result: exp(r)

  ! TODO: Specify as hex to avoid ambiguous rounding
  real(kind=real64), parameter :: INV_LN2 = 1.4426950408889634_real64
  !real, parameter :: INV_LN2 = 1.44269504088896340735992468100189204
  real(kind=real64), parameter :: TWO_INV_LN2 = 2.8853900817779268_real64
  !real, parameter :: TWO_INV_LN2 = 2.88539008177792681471984936200378409
  real(kind=real64), parameter :: LN2_HI = 6.93147180369123816490e-01_real64
  real(kind=real64), parameter :: LN2_LO = 1.90821492927058770002e-10_real64

  ! For Chebyshev evaluation

  ! Scale to [-0.5 ln 2, 0.5 ln 2]
  K = nint(x * INV_LN2)
  r = (x - K * LN2_HI) - K * LN2_LO

  ! Compute exp(r)
  ! In order to force exp(0) = 1, we actually approximate (exp(r) - 1) / r.
  e = 1. + r * expm1_x_remez(r * TWO_INV_LN2)

  ! Descale 
  a = scale(e, K)

end function exp_cr

pure function expm1_x_remez(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: b0, b1, b2
  integer :: n

  !real, parameter :: t(0:10) = [ &
  !  1.0100546303736508e+00, &
  !  1.7459202119340825e-01, &
  !  1.0069712517838964e-02, &
  !  4.3580091845530559e-04, &
  !  1.5092920250821286e-05, &
  !  4.3566706207388448e-07, &
  !  1.0780554497188584e-08, &
  !  2.3343806108994727e-10, &
  !  4.4933851094935047e-12, &
  !  7.7884991972921424e-14, &
  !  1.2375905543948744e-15 &
  !]

  real(kind=real64), parameter :: t(0:11) = [ &
    1.0100546303736506e+00_real64, &
    1.7459202119340830e-01_real64, &
    1.0069712517839037e-02_real64, &
    4.3580091845535839e-04_real64, &
    1.5092920250787356e-05_real64, &
    4.3566706211233559e-07_real64, &
    1.0780554531543968e-08_real64, &
    2.3343813394501564e-10_real64, &
    4.4933658510432855e-12_real64, &
    7.7841956207353226e-14_real64, &
    1.1969180485188496e-15_real64, &
    3.7974225925737191e-17_real64 &
  ]

  b1 = 0.0
  b2 = 0.0

  do n = 11,1,-1
    b0 = t(n) + 2. * x * b1 - b2
    b2 = b1
    b1 = b0
  enddo

  e = t(0) + x * b1 - b2
end function expm1_x_remez

end module exp_repro
