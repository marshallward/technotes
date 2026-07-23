module exp_repro

use, intrinsic :: iso_fortran_env, only : int32, int64, real64, real128

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

! TODO: Portable versions?
integer(int64), parameter :: pos_inf_bits = int(z'7FF0000000000000', int64)
  !< IEEE 64-bit +Inf
integer(int64), parameter :: neg_inf_bits = int(z'FFF0000000000000', int64)
  !< IEEE 64-bit -Inf

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
  real(kind=real64), value, intent(in) :: x
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

  real(kind=real64), parameter :: LN2 = 0.6931471805599453_real64
  real(kind=real64), parameter :: TWO_INV_LN2 = 2.8853900817779268_real64

  ! Experimental
  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52

  ! Subnormals?
  integer(int64) :: j, fb

  ! Even more experimental
  real(kind=real64), parameter :: INV_LN2 = 1.4426950408889634_real64
  real(kind=real64), parameter :: LN2_HI = 6.93147180369123816490e-01_real64
  real(kind=real64), parameter :: LN2_LO = 1.90821492927058770002e-10_real64

  ! Further subdivide [-(1/2N) ln2, +(1/2N) ln2], use tables to scale back up.
  integer, parameter :: NTABLE = 1
  real(kind=real64), parameter :: I_NTABLE = 1._real64 / real(NTABLE, real64)
  real(kind=real64), parameter :: TABLE_INV_LN2 = NTABLE * INV_LN2
  real(kind=real64), parameter :: TABLE_LN2_HI = I_NTABLE * LN2_HI
  real(kind=real64), parameter :: TABLE_LN2_LO = I_NTABLE * LN2_LO

  real(kind=real64) :: Z
  integer(kind=int64) :: iz, kz
  integer(int32) :: table_index

  integer :: i
  real(kind=real64), parameter :: exp2_table(0:NTABLE-1) = &
      [(2._real64**(real(i, real64) / real(NTABLE, real64)), i=0,NTABLE-1)]

  integer(int64), parameter :: index_mask = int(NTABLE - 1, int64)

  ! XXX: Use this to test performance without scaling.
  !a = exp_remez_estrin_5(x)
  !return

  ! Scale to [-1/2N ln 2, +1/2N ln 2]

  !**!! Crude implementation of anint(x_ln2).
  !**!! NOTE: In x86 GCC, this favors vround instructions over round() calls.
  !**!K = aint(x_ln2 + sign(0.5_real64, x_ln2))

  !**!! TODO: This is slightly faster but may be removed by optimization
  !**!!K = (x_ln2 + round_bias) - round_bias

  !**!! Cody-Waite split
  !**!! This decomposition preserves lower bits after integer cancellation.
  !**!! TODO: Explain this better
  !**!! NOTE: Compilers may optimize this to `x - K * (LN2_HI + L2_LI)` which is no
  !**!!   better than x - K * LN2
  !**!r = (x - K * LN2_HI) - K * LN2_LO

  !**!! This is less accurate if abs(K) is large, but is faster
  !**!!r = x - K * LN2

  !**!! NOTE: Chebyshev polynomial is normalized to [-1,1] so we must rescale.
  !**!!e = 1. + r * exp_remez_chebyshev(r * TWO_INV_LN2)
  !**!!e = exp_remez_estrin_9(r)
  !**!e = exp_remez_estrin_10(r)
  !**!!e = exp_remez_estrin_11(r)
  !**!!e = exp_taylor_horner(r)
  !**!!e = exp_taylor_estrin(r)

  ! *** Scale x to r = x - nint(N*x/ln2)

  ! Compute N*x/ln2
  x_ln2 = x * TABLE_INV_LN2

  ! Round to nearest integer: Z = nint(x_ln2)

  ! NOTE: This is fast but may be reduced to Z = x_ln2
  Z = (x_ln2 + round_bias) - round_bias

  ! This is a safer alternative
  !Z = aint(x_ln2 + sign(0.5_real64, x_ln2))

  ! Decompose Z = N * K + table_index, where K = nint(x / ln2)

  ! Extract the lower bits of Z to determine the subscale
  iz = transfer(Z + round_bias, int64_mold)
  table_index = iand(iz, index_mask)

  !K = (Z - real(table_index, real64)) * 0.03125_real64
  kz = iand(iz, not(index_mask))
  K = (transfer(kz, real64_mold) - round_bias) * I_NTABLE

  ! Cody-Waite splitting
  ! Residual after subtracting Z*ln(2)/NTABLE
  ! NOTE: May be reduced to `x - K * (TABLE_LN2_HI + TABLE_L2_LO)` which is no
  !   better than x - Z * TABLE_LN2
  r = (x - Z * TABLE_LN2_HI) - Z * TABLE_LN2_LO

  ! *** Compute exp(r) ***!

  ! Approximate exp(r), then restore the tabulated fractional power.

  ! Use these for N=32
  !e = exp2_table(table_index) * exp_remez_estrin_4(r)
  !e = exp2_table(table_index) * exp_remez_estrin_5(r)

  ! Use with N=1
  e = exp2_table(table_index) * exp_remez_estrin_10(r)

  ! Taylor functions are range-agnostic (i.e. only good near zero!)
  !e = exp2_table(table_index) * exp_taylor_estrin_6(r)

  ! *** Descale the value ***!

  ! Over/underflow subnormal offset
  ! TODO: Keep as real?
  j = merge(1022_int64, 0_int64, K < -1020.0_real64) &
      - merge(1022_int64, 0_int64, K > 1020.0_real64)

  ! Get the bitform of e.
  eb = transfer(e, int64_mold)

  ! Shift K to the significand's least significant bits.
  ! kb is now the integer value of K (with over/underflow correction).
  kb = transfer(K + round_bias, int64_mold)

  ! Apply the K exponent to e's exponent.
  eb = eb + shiftl(kb + j, expbit)
  a  = transfer(eb, 1.0_real64)

  ! Undo the over/underflow correction
  fb = shiftl(1023_int64 - j, expbit)

  ! Apply rescaling if needed
  a = a * transfer(fb, 1.0_real64)

  !*** IEEE corrections ***!

  ! TODO: These correct Inf values but do not account for incorrect signals.

  ! Set exp(-Inf) = 0.
  a = merge(0._real64, a, transfer(x, int64_mold) == neg_inf_bits)

  ! Set exp(+Inf) = +Inf
  r = transfer(pos_inf_bits, real64_mold)
  a = merge(r, a, transfer(x, int64_mold) == pos_inf_bits)
end function exp_cr


elemental function exp_cr_fast(x) result(a)
  !$acc routine seq
  real(kind=real64), value, intent(in) :: x
    !< Input value [nondim]
  real(kind=real64) :: a
    !< Exponential of x [nondim]

  real(kind=real64) :: r
  real(kind=real64) :: x_ln2
  real(kind=real64) :: K
  real(kind=real64) :: e

  integer(kind=int64) :: eb, kb

  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52

  real(kind=real64), parameter :: INV_LN2 = 1.4426950408889634_real64
  real(kind=real64), parameter :: LN2_HI = 6.93147180369123816490e-01_real64
  real(kind=real64), parameter :: LN2_LO = 1.90821492927058770002e-10_real64

  integer, parameter :: NTABLE = 1
  real(kind=real64), parameter :: I_NTABLE = 1._real64 / real(NTABLE, real64)
  real(kind=real64), parameter :: TABLE_INV_LN2 = NTABLE * INV_LN2
  real(kind=real64), parameter :: TABLE_LN2_HI = I_NTABLE * LN2_HI
  real(kind=real64), parameter :: TABLE_LN2_LO = I_NTABLE * LN2_LO

  real(kind=real64) :: Z
  integer(kind=int64) :: iz, kz
  integer(int32) :: table_index

  integer :: i
  real(kind=real64), parameter :: exp2_table(0:NTABLE-1) = &
      [(2._real64**(real(i, real64) / real(NTABLE, real64)), i=0,NTABLE-1)]
  integer(int64), parameter :: index_mask = int(NTABLE - 1, int64)

  x_ln2 = x * TABLE_INV_LN2
  Z = (x_ln2 + round_bias) - round_bias

  iz = transfer(Z + round_bias, int64_mold)
  table_index = iand(iz, index_mask)

  kz = iand(iz, not(index_mask))
  K = (transfer(kz, real64_mold) - round_bias) * I_NTABLE

  r = (x - Z * TABLE_LN2_HI) - Z * TABLE_LN2_LO
  e = exp2_table(table_index) * exp_remez_estrin_10(r)

  eb = transfer(e, int64_mold)
  kb = transfer(K + round_bias, int64_mold)
  eb = eb + shiftl(kb, expbit)
  a = transfer(eb, 1.0_real64)
end function exp_cr_fast


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


! NOTE: Reduced range: +/- 1/32 0.5 ln2
pure function exp_remez_estrin_4(x) result(e)
  real(real64), intent(in) :: x
    !< Reduced argument
  real(real64) :: e
    !< Approximation of exp(x)

  real(real64), parameter :: c(0:4) = [ &
    9.99999999999999889e-01_real64, &
    4.99999999995007050e-01_real64, &
    1.66666666671947628e-01_real64, &
    4.16668555121924925e-02_real64, &
    8.33331131811317072e-03_real64 ]

  real(real64) :: x2, x4
  real(real64) :: a0, a1
  real(real64) :: p

  x2 = x * x
  x4 = x2 * x2

  a0 = c(0) + x * c(1)
  a1 = c(2) + x * c(3)

  p = (a0 + x2 * a1) + x4 * c(4)

  e = 1.0_real64 + x * p
end function exp_remez_estrin_4


! NOTE: Reduced range: +/- 1/32 0.5 ln2
pure function exp_remez_estrin_5(x) result(e)
  real(real64), intent(in) :: x
    !< Reduced argument
  real(real64) :: e
    !< Approximation of exp(x)

  real(real64), parameter :: c(0:5) = [ &
    9.99999999999999778e-01_real64, &
    5.00000000000007438e-01_real64, &
    1.66666666662346141e-01_real64, &
    4.16666663113674576e-02_real64, &
    8.33338618897187799e-03_real64, &
    1.39156572960508467e-03_real64 ]

  real(real64) :: x2, x4
  real(real64) :: a0, a1, a2
  real(real64) :: p

  x2 = x * x
  x4 = x2 * x2

  a0 = c(0) + x * c(1)
  a1 = c(2) + x * c(3)
  a2 = c(4) + x * c(5)

  p = (a0 + x2 * a1) + x4 * a2

  e = 1.0_real64 + x * p
end function exp_remez_estrin_5


pure function exp_remez_estrin_10(x) result(e)
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
end function exp_remez_estrin_10


pure function exp_remez_estrin_9(x) result(e)
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
end function exp_remez_estrin_9


pure function exp_remez_estrin_11(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(real64), parameter :: c(0:11) = [ &
      9.99999999999999667e-01_real64, &
      4.99999999999999500e-01_real64, &
      1.66666666666674040e-01_real64, &
      4.16666666667472874e-02_real64, &
      8.33333333337431259e-03_real64, &
      1.38888888534349124e-03_real64, &
      1.98412691045814590e-04_real64, &
      2.48016479357667924e-05_real64, &
      2.75583751072956796e-06_real64, &
      2.75099984928059077e-07_real64, &
      2.46463609302820685e-08_real64, &
      3.51693798509217157e-09_real64 &
  ]

  real(kind=real64) :: x2, x4, x8
  real(kind=real64) :: a0, a1, a2, a3, a4, a5
  real(kind=real64) :: b0, b1, b2
  real(kind=real64) :: p

  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  a0 = c(0)  + c(1) * x
  a1 = c(2)  + c(3) * x
  a2 = c(4)  + c(5) * x
  a3 = c(6)  + c(7) * x
  a4 = c(8)  + c(9) * x
  a5 = c(10) + c(11) * x

  b0 = a0 + a1 * x2
  b1 = a2 + a3 * x2
  b2 = a4 + a5 * x2

  p = (b0 + b1 * x4) + b2 * x8

  e = 1. + x * p
end function exp_remez_estrin_11


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


pure function exp_taylor_estrin_6(x) result(e)
  real(real64), intent(in) :: x
  real(real64) :: e
  real(real64) :: x2, x4
  real(real64) :: p0, p1, p2

  x2 = x * x
  x4 = x2 * x2

  p0 = 1.0_real64 + x
  p1 = 0.5_real64 + x * (1.0_real64 / 6.0_real64)
  p2 = (1.0_real64 / 24.0_real64) + &
       x * (1.0_real64 / 120.0_real64)

  e = (p0 + x2 * p1) + &
      x4 * (p2 + x2 * (1.0_real64 / 720.0_real64))
end function exp_taylor_estrin_6


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
