module exp_repro

use, intrinsic :: iso_fortran_env, only : int64, real64, real128

implicit none

! Scalar molds
integer(kind=int64), parameter :: int64_mold = 0
real(kind=real64), parameter :: real64_mold = 0.

! Floating point model, where bit layout from high to low is (sign, exp, frac)
integer, parameter :: bias = maxexponent(real64_mold) - 1
  !< The double precision exponent offset
integer, parameter :: explen = storage_size(real64_mold) - digits(real64_mold)
  !< Bit size of exponent
integer, parameter :: expbit = digits(real64_mold) - 1
  !< Position of lowest exponent bit

contains

pure subroutine exp_1d(x, a)
  real(kind=real64), intent(in) :: x(:)
  real(kind=real64), intent(out) :: a(:)
  integer :: i

  do i = 1,size(x)
    a(i) = exp_cr(x(i))
  enddo
end subroutine exp_1d


pure subroutine exp_1d_do_c(x, a)
  real(kind=real64), intent(in) :: x(:)
  real(kind=real64), intent(out) :: a(:)
  integer :: i

  do concurrent (i = 1:size(x))
    a(i) = exp_cr(x(i))
  enddo
end subroutine exp_1d_do_c


elemental function exp_cr(x) result(a)
  !$acc routine seq
  real(kind=real64), intent(in) :: x
    !< Input value [nondim]
  real(kind=real64) :: a
    !< Exponential of x [nondim]

  real(kind=real64) :: r
    ! Rescaled value of x; r = x - K * ln 2
  real(kind=real64) :: x_ln2
    ! Intermediate scaling value, x / ln2
  real(kind=real64) :: K
    ! Scaling factor, where a = exp(x) = 2**K exp(r)
    ! Stored as a real to prevent unnecessary type conversion
  real(kind=real64) :: e
    ! Scaled result; e = exp(r)

  ! Descale
  integer(kind=int64) :: eb, kb
    ! Bit representation of e and (round-biased) K
  integer(kind=int64) :: ek
    ! Exponent of descaled exponent

  ! TODO: Specify as hex to avoid ambiguous rounding
  real(kind=real64), parameter :: LN2 = 0.6931471805599453_real64
  real(kind=real64), parameter :: INV_LN2 = 1.4426950408889634_real64
  !real, parameter :: INV_LN2 = 1.44269504088896340735992468100189204
  real(kind=real64), parameter :: TWO_INV_LN2 = 2.8853900817779268_real64
  !real, parameter :: TWO_INV_LN2 = 2.88539008177792681471984936200378409
  real(kind=real64), parameter :: LN2_HI = 6.93147180369123816490e-01_real64
  real(kind=real64), parameter :: LN2_LO = 1.90821492927058770002e-10_real64

  ! Experimental
  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52

  real(kind=real64) :: s

  !$omp declare simd

  ! XXX: Use this for profiling without scaling.
  !a = exp_remez_estrin(x)

  ! Scale to [-0.5 ln 2, 0.5 ln 2]

  x_ln2 = x * INV_LN2

  ! Crude implementation of anint(x_ln2).
  ! NOTE: In x86 GCC, this favors vround instructions over round() calls.
  K = aint(x_ln2 + sign(0.5_real64, x_ln2))

  ! TODO: This is faster but may be removed by optimization
  !K = (x_ln2 + round_bias) - round_bias

  ! Cody-Waite split
  ! This decomposition preserves lower bits after integer cancellation.
  ! TODO: Explain this better
  ! NOTE: This may be optimized to `x - K * (LN2_HI + L2_LI)` which is no
  !   better than x - K * LN2
  r = (x - K * LN2_HI) - K * LN2_LO

  ! This is less accurate if abs(K) is large, but is faster
  !r = x - K * LN2

  ! NOTE: Chebyshev polynomial is normalized to [-1,1] so we have to rescale.
  !e = 1. + r * exp_remez_chebyshev(r * TWO_INV_LN2)
  e = exp_remez_estrin(r)
  !e = exp_remez_estrin9(r)
  !e = exp_taylor_horner(r)
  !e = exp_taylor_estrin(r)

  ! Descale the value

	if (K >= -1020.0_real64 .and. K <= 1020.0_real64) then
    ! Get the bitform of e.
	  eb = transfer(e, int64_mold)

    ! Shift K to the significand's least significant bits.
    ! kb is now the integer value of K (excepting special IEEE values).
	  kb = transfer(K + round_bias, int64_mold)

    ! Apply the K exponent to e's exponent.
	  eb = eb + shiftl(kb, expbit)
	  a = transfer(eb, 1.0_real64)
	else
    ! For exceptional values, fallback to intrinsics for rescaling.
		a = scale(e, int(K))
	end if
end function exp_cr


pure function exp_remez_chebyshev(x) result(e)
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

  ! NOTE: In order to force exp(0) = 1, we estimate (exp(r) - 1) / r.

  b1 = 0.0
  b2 = 0.0

  do n = 10,1,-1
    b0 = t(n) + 2. * x * b1 - b2
    b2 = b1
    b1 = b0
  enddo

  e = t(0) + x * b1 - b2
end function exp_remez_chebyshev


pure function exp_remez_estrin(x) result(e)
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

  ! NOTE: In order to force exp(0) = 1, we estimate (exp(r) - 1) / r.

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

  e = 1 + x * (q0 + x4 * q1 + x8 * (b4 + x2 * c(10)))
end function exp_remez_estrin


pure function exp_remez_estrin9(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(real64), parameter :: c(0:9) = [ &
      1.00000000000000133e+00_real64, &
      5.00000000000001776e-01_real64, &
      1.66666666666236918e-01_real64, &
      4.16666666664353008e-02_real64, &
      8.33333336355678196e-03_real64, &
      1.38888889752610879e-03_real64, &
      1.98411961311610672e-04_real64, &
      2.48014555669475914e-05_real64, &
      2.76301733424592505e-06_real64, &
      2.76447983021274086e-07_real64 &
  ]

  real(real64) :: x2, x4, x8, t0_0, t0_1, t0_2, t0_3, t0_4, &
                  t1_0, t1_1, t1_2, t2_0, t2_1

  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  t0_0 = c(0) + x * c(1)
  t0_1 = c(2) + x * c(3)
  t0_2 = c(4) + x * c(5)
  t0_3 = c(6) + x * c(7)
  t0_4 = c(8) + x * c(9)

  t1_0 = t0_0 + x2 * t0_1
  t1_1 = t0_2 + x2 * t0_3
  t1_2 = t0_4

  t2_0 = t1_0 + x4 * t1_1
  t2_1 = t1_2

  e = 1. + x * (t2_0 + x8 * t2_1)
end function exp_remez_estrin9


elemental function exp_taylor_horner(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64), parameter :: c(0:12) = [ &
    1.0_real64, &
    1.0_real64, &
    1.0_real64 / 2.0_real64, &
    1.0_real64 / 6.0_real64, &
    1.0_real64 / 24.0_real64, &
    1.0_real64 / 120.0_real64, &
    1.0_real64 / 720.0_real64, &
    1.0_real64 / 5040.0_real64, &
    1.0_real64 / 40320.0_real64, &
    1.0_real64 / 362880.0_real64, &
    1.0_real64 / 3628800.0_real64, &
    1.0_real64 / 39916800.0_real64, &
    1.0_real64 / 479001600.0_real64 &
  ]
    !< Taylor coefficients 1/n!

  integer :: n

  e = 0
  do n = 12,0,-1
    e = x * e + c(n)
  enddo
end function exp_taylor_horner


elemental function exp_taylor_estrin(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64), parameter :: c(0:12) = [ &
    1.0_real64, &
    1.0_real64, &
    1.0_real64 / 2.0_real64, &
    1.0_real64 / 6.0_real64, &
    1.0_real64 / 24.0_real64, &
    1.0_real64 / 120.0_real64, &
    1.0_real64 / 720.0_real64, &
    1.0_real64 / 5040.0_real64, &
    1.0_real64 / 40320.0_real64, &
    1.0_real64 / 362880.0_real64, &
    1.0_real64 / 3628800.0_real64, &
    1.0_real64 / 39916800.0_real64, &
    1.0_real64 / 479001600.0_real64 &
  ]
    !< Taylor coefficients 1/n!

  real(kind=real64) :: x2, x4, x8
  real(kind=real64) :: p0, p1, p2, p3, p4, p5
  real(kind=real64) :: q0, q1, q2
  real(kind=real64) :: r0, r1

  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  p0 = c(0)  + x * c(1)
  p1 = c(2)  + x * c(3)
  p2 = c(4)  + x * c(5)
  p3 = c(6)  + x * c(7)
  p4 = c(8)  + x * c(9)
  p5 = c(10) + x * c(11)

  q0 = p0 + x2 * p1
  q1 = p2 + x2 * p3
  q2 = p4 + x2 * p5

  r0 = q0 + x4 * q1
  r1 = q2 + x4 * c(12)

  e = r0 + x8 * r1
end function exp_taylor_estrin

end module exp_repro
