module exp_repro

use iso_fortran_env, only : int64, real64, real128

implicit none

! Scalar molds
integer(kind=int64), parameter :: int64_mod = 0
real(kind=real64), parameter :: real64_mold = 0.

! Floating point model, where bit layout from high to low is (sign, exp, frac)
integer, parameter :: bias = maxexponent(real64_mold) - 1
  !< The double precision exponent offset
integer, parameter :: explen = storage_size(real64_mold) - digits(real64_mold)
  !< Bit size of exponent
integer, parameter :: expbit = digits(real64_mold) - 1
  !< Position of lowest exponent bit

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

  ! Descale
  integer(kind=int64) :: eb
    ! Bit representation of e
  integer(kind=int64) :: ek
    ! Exponent of descaled exponent

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

  ! The Chebyshev polynomial is normalized to [-1,1] so we have to rescale.
  !e = 1. + r * expm1_x_remez(r * TWO_INV_LN2)

  ! The Estrin form has no such restriction
  e = 1. + r * expm1_x_estrin(r)

  !!!! Descale
  ! (TODO: Move to function)

  !*!! Intrinsic is unfortunately not always inlined
  !*!a = scale(e, K)

  ! Get the bitform of e
  eb = transfer(e, 1_int64)

  ! Extract the exponent (with bias)
  ek = ibits(eb, expbit, explen)

  ! Apply the descaled exponent
  call mvbits(ek + K, 0, explen, eb, expbit)

  a = transfer(eb, 1._real64)
end function exp_cr


pure function expm1_x_remez(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: b0, b1, b2
  integer :: n

  ! n=10 polynomial
  real(kind=real64), parameter :: t(0:10) = [ &
    1.0100546303736508e+00_real64, &
    1.7459202119340825e-01_real64, &
    1.0069712517838964e-02_real64, &
    4.3580091845530559e-04_real64, &
    1.5092920250821286e-05_real64, &
    4.3566706207388448e-07_real64, &
    1.0780554497188584e-08_real64, &
    2.3343806108994727e-10_real64, &
    4.4933851094935047e-12_real64, &
    7.7884991972921424e-14_real64, &
    1.2375905543948744e-15_real64 &
  ]

  ! n=11 polynomial
  !real(kind=real64), parameter :: t(0:11) = [ &
  !  1.0100546303736506e+00_real64, &
  !  1.7459202119340830e-01_real64, &
  !  1.0069712517839037e-02_real64, &
  !  4.3580091845535839e-04_real64, &
  !  1.5092920250787356e-05_real64, &
  !  4.3566706211233559e-07_real64, &
  !  1.0780554531543968e-08_real64, &
  !  2.3343813394501564e-10_real64, &
  !  4.4933658510432855e-12_real64, &
  !  7.7841956207353226e-14_real64, &
  !  1.1969180485188496e-15_real64, &
  !  3.7974225925737191e-17_real64 &
  !]

  b1 = 0.0
  b2 = 0.0

  do n = 10,1,-1
    b0 = t(n) + 2. * x * b1 - b2
    b2 = b1
    b1 = b0
  enddo

  e = t(0) + x * b1 - b2
end function expm1_x_remez


pure function expm1_x_estrin(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: x2, x4, x8
    ! Powers of x
  real(kind=real64) :: b0, b1, b2, b3, b4
  real(kind=real64) :: q0, q1
    ! Dont know yet, maybe levels in the tree

  real(kind=real64), parameter :: c(0:10) = [ &
    9.99999999999999889e-01_real64, &
    5.00000000000003109e-01_real64, &
    1.66666666666674762e-01_real64, &
    4.16666666663442417e-02_real64, &
    8.33333333278445630e-03_real64, &
    1.38888889936273166e-03_real64, &
    1.98412708336599035e-04_real64, &
    2.48014408391593750e-05_real64, &
    2.75566321051561213e-06_real64, &
    2.76488381734320880e-07_real64, &
    2.52312075296588332e-08_real64 &
  ]

  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  b0 = c(0) + x * c(1)
  b1 = c(2) + x * c(3)
  b2 = c(4) + x * c(5)
  b3 = c(6) + x * c(7)
  b4 = c(8) + x * c(9)

  q0 = b0 + x2 * b1
  q1 = b2 + x2 * b3

  e = q0 + x4 * q1 + x8 * (b4 + x2 * c(10))
end function expm1_x_estrin

end module exp_repro
