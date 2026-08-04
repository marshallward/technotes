#include "exp_repro.h"

!> Reproducible exponential function
!!
!! This module provides a reproducible implementation of exp() that produces
!! bitwise identical results across platforms while matching IEEE signal
!! behavior (overflow, underflow, inexact) with the intrinsic exp().
!!
!! Performance is approximately 1.17x Intel SVML with full IEEE compliance.
module exp_repro_mod

use, intrinsic :: iso_fortran_env, only : int32, int64, real64

implicit none
private

public :: exp_repro

! Scalar mold for transfer()
integer(kind=int64), parameter :: int64_mold = 0
real(kind=real64), parameter :: real64_mold = 0.

! IEEE 754 double precision layout
integer, parameter :: expbit = digits(1._real64) - 1
  !< Position of lowest exponent bit (52)
integer, parameter :: expwidth = storage_size(1._real64) - expbit - 1
  !< Number of exponent bits (11)

! IEEE 754 special values
integer(int64), parameter :: pos_inf_bits = ishft(2_int64**expwidth - 1, expbit)
  !< IEEE 64-bit +Inf
integer(int64), parameter :: neg_inf_bits = ior(pos_inf_bits, ishft(-1_int64, 63))
  !< IEEE 64-bit -Inf

contains

!> Reproducible exponential function
!!
!! Computes exp(x) with bitwise reproducibility across platforms.
!! IEEE signal behavior matches intrinsic exp():
!!   - normal: inexact
!!   - overflow: overflow + inexact
!!   - underflow: underflow + inexact
!!   - +/-Inf: no signals
!!   - NaN: no signals (NaN propagates)
elemental function exp_repro(x) result(a)
  !$acc routine seq
  real(kind=real64), value, intent(in) :: x
    !< Input value
  real(kind=real64) :: a
    !< exp(x)

  ! Finite argument limits where exp(x) overflows or rounds to zero
  real(kind=real64), parameter :: max_exp_arg = log(huge(real64_mold))
  real(kind=real64), parameter :: min_exp_arg = log(tiny(real64_mold)) &
      - real(digits(real64_mold), kind=real64) * log(2.0_real64)

  ! Range reduction: exp(x) = 2^K * exp(r), where r in [-ln2/2, ln2/2]
  real(kind=real64) :: r
  real(kind=real64) :: K
  real(kind=real64) :: e

  ! Bit manipulation for fast 2^K scaling
  integer(kind=int64) :: eb, kb
  integer(int64) :: j, fb, xb

  ! Constants for range reduction
  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52
  real(kind=real64), parameter :: INV_LN2 &
      = transfer(int(z'3FF71547652B82FE', int64), real64_mold)
      ! 1.4426950408889634073599...

  ! Cody-Waite constants for accurate range reduction
  real(kind=real64), parameter :: LN2_HI &
      = transfer(int(z'3FE62E42FEE00000', real64), real64_mold)
      ! 6.93147180369123816490e-01
  real(kind=real64), parameter :: LN2_LO &
      = transfer(real(z'3DEA39EF35793C76', int64), real64_mold)
      ! 1.90821492927058770002e-10

  ! Further subdivide [-(1/2N) ln2, +(1/2N) ln2], use tables to scale back up.
  integer, parameter :: NTABLE = 1
  real(kind=real64), parameter :: I_NTABLE = 1._real64 / real(NTABLE, real64)
  real(kind=real64), parameter :: TABLE_INV_LN2 = NTABLE * INV_LN2
  real(kind=real64), parameter :: TABLE_LN2_HI = I_NTABLE * LN2_HI
  real(kind=real64), parameter :: TABLE_LN2_LO = I_NTABLE * LN2_LO

  integer :: i
  real(kind=real64), parameter :: exp2_table(0:NTABLE-1) = &
      [(2._real64**(real(i, real64) / real(NTABLE, real64)), i=0,NTABLE-1)]
  integer(int64), parameter :: index_mask = int(NTABLE - 1, int64)

  real(kind=real64) :: x_ln2

  ! Tables
  real(kind=real64) :: Z
  integer(kind=int64) :: zb
  integer(int32) :: table_index

  logical :: nonfinite
    ! True if finite x is outside the representable range of exp(x)

  xb = transfer(x, int64_mold)
  nonfinite = iand(xb, pos_inf_bits) == pos_inf_bits

  if (nonfinite) then
    ! exp(-Inf) = 0, pass-through for +Inf and +/-NaN.
    ! Use x + x = x to trigger Invalid for signaled NaNs.
    a = merge(0._real64, x + x, xb == neg_inf_bits)
    return
  endif

  ! *** Range reduction ***

  ! Compute K = nint(x / ln2)
  x_ln2 = x * INV_LN2
  K = NEAREST_INT(x_ln2)

  K = min(max(K, -1024._real64), 1024._real64)

  ! Cody-Waite: r = x - K*ln2 with extended precision
  r = (x - K * LN2_HI) - K * LN2_LO

  !! Wrong, just testing performance
  !if (x < log(tiny(real64_mold))) then
  !  a = transfer(pos_inf_bits, real64_mold)
  !  a = a * a
  !  return
  !endif

  !*!! Table method
  !*!x_ln2 = x * TABLE_INV_LN2
  !*!Z = NEAREST_INT(x_ln2)

  !*!zb = transfer(Z + round_bias, int64_mold)
  !*!table_index = iand(zb, index_mask)

  !*!!K = (Z - real(table_index, real64)) / NTABLE
  !*!kb = iand(zb, not(index_mask))
  !*!K = (transfer(kb, real64_mold) - round_bias) * I_NTABLE

  !*!! Tabled Cody-Waite
  !*!r = (x - Z * TABLE_LN2_HI) - Z * TABLE_LN2_LO

  ! *** Polynomial approximation ***

  ! exp(r) for r in [-ln2/2, ln2/2] using degree-10 Remez minimax polynomial
  !e = exp_remez_estrin_9(r)
  !e = exp_remez_estrin_10(r)
  !e = exp_remez_estrin_11(r)
  !e = exp_remez_estrin_12(r)
  !e = exp_remez_horner_9(r)
  e = exp_remez_horner_10(r)
  !e = exp_remez_horner_10_constrained(r)
  !e = exp_remez_horner_11(r)
  !e = exp_remez_horner_12(r)

  ! Tabled version
  !e = exp2_table(table_index) * exp_remez_horner_10(r)

  ! *** Scaling: multiply by 2^K ***

  ! Handle over/underflow by splitting extreme K values
  j = merge(1022_int64, 0_int64, K < -1020.0_real64) &
      + merge(-1022_int64, 0_int64, K > 1020.0_real64)

  eb = transfer(e, int64_mold)
  kb = transfer(K + round_bias, int64_mold)

  eb = eb + ishft(kb + j, expbit)
  a = transfer(eb, real64_mold)

  ! Apply correction factor for extreme K (triggers IEEE over/underflow signal)
  fb = ishft(1023_int64 - j, expbit)
  a = a * transfer(fb, real64_mold)
end function exp_repro


!> Degree-10 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Estrin's scheme for instruction-level parallelism.
!! Coefficients generated by Sollya.
pure function exp_remez_estrin_10(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: x2, x4, x8
  real(kind=real64) :: b0, b1, b2, b3, b4
  real(kind=real64) :: q0, q1

  ! Remez minimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c(0:10) = [ &
    1.0_real64, &
    5.00000000000000555e-1_real64, &
    1.66666666666666074e-1_real64, &
    4.16666666665738844e-2_real64, &
   !8.33333333337716447e-3_real64, &
    ! This reduces relative error but increases absolute error
    8.33333333337716621e-03_real64, &
    1.38888889322647565e-3_real64, &
    1.98412697469850182e-4_real64, &
    2.48015045964426143e-5_real64, &
    2.75573817985163132e-6_real64, &
    2.76262647076892519e-7_real64, &
    2.50621020021886396e-8_real64 &
  ]

  ! Estrin's scheme: evaluate polynomial with maximum parallelism
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

  ! Final assembly: 1 + x * p(x) ensures exp(0) = 1 exactly
  ! Ambiguous
  !e = 1 + x * (q0 + x4 * q1 + x8 * (b4 + x2 * c(10)))
  ! Better but FMA ambiguity
  !e = 1 + x * (q0 + (x4 * q1 + x8 * (b4 + x2 * c(10))))
  ! Best?
  e = 1 + x * ((q0 + x4 * q1) + x8 * (b4 + x2 * c(10)))
end function exp_remez_estrin_10


!> Degree-12 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_12(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0  = 1.0_real64
  real(kind=real64), parameter :: c1  = 0.5_real64
  real(kind=real64), parameter :: c2  = 0.1666666666666666574148081281236954964697360992431640625_real64
  real(kind=real64), parameter :: c3  = 4.1666666666666678231489839845380629412829875946044921875e-2_real64
  real(kind=real64), parameter :: c4  = 8.33333333333420578359351793551468290388584136962890625e-3_real64
  real(kind=real64), parameter :: c5  = 1.38888888888774706016626669935476456885226070880889892578125e-3_real64
  real(kind=real64), parameter :: c6  = 1.984126983838181126820060518056720866297837346792221069335937e-4_real64
  real(kind=real64), parameter :: c7  = 2.480158734051678335566677724433048979335580952465534210205078e-5_real64
  real(kind=real64), parameter :: c8  = 2.755732358282358854077720433650711129303090274333953857421875e-6_real64
  real(kind=real64), parameter :: c9  = 2.755725932837027351290253886217929135682425112463533878326416e-7_real64
  real(kind=real64), parameter :: c10 = 2.504903185531009510483624108320710455188873311271890997886658e-8_real64
  real(kind=real64), parameter :: c11 = 2.091945316636953120202141592222887245267060052356100641191006e-9_real64
  real(kind=real64), parameter :: c12 = 1.689299804577804335305572395803039298378678267908981069922447e-10_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c12
  p = c11 + x * p
  p = c10 + x * p
  p = c9 + x * p
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_12


!> Degree-10 constrained polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1 - x) / x^2, then computes 1 + x + x^2 * p(x).
!! This form keeps the leading terms exact.
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax with c[0]=1 constrained.
pure function exp_remez_horner_10_constrained(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: x2, p

  ! fpminimax coefficients for (exp(x) - 1 - x) / x^2 on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0 = 0.50000000000000011102230246251565404236316680908203125_real64
  real(kind=real64), parameter :: c1 = 0.1666666666666664353702032030923874117434024810791015625_real64
  real(kind=real64), parameter :: c2 = 4.1666666666623275450120900131878443062305450439453125e-2_real64
  real(kind=real64), parameter :: c3 = 8.333333333356306160677462457897490821778774261474609375e-3_real64
  real(kind=real64), parameter :: c4 = 1.38888889174418350171136271598015810013748705387115478515625e-3_real64
  real(kind=real64), parameter :: c5 = 1.984126978475939985253201358617047844745684415102005004882812e-4_real64
  real(kind=real64), parameter :: c6 = 2.480152106658596513222986290614358040329534560441970825195312e-5_real64
  real(kind=real64), parameter :: c7 = 2.755735516039519399219756556895788435213034972548484802246094e-6_real64
  real(kind=real64), parameter :: c8 = 2.762016599747405094317448739915654698506841668859124183654785e-7_real64
  real(kind=real64), parameter :: c9 = 2.506835135920223397333950524600021392274129539146088063716888e-8_real64

  x2 = x * x

  ! Horner's scheme for p(x) = c0 + x*(c1 + x*(c2 + ...))
  p = c9
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x + x^2 * p(x)
  e = 1 + x + x2 * p
end function exp_remez_horner_10_constrained


!> Degree-12 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme to match Sollya fpminimax optimization.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_estrin_12(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0  = 1.0_real64
  real(kind=real64), parameter :: c1  = 0.5_real64
  real(kind=real64), parameter :: c2  = 0.1666666666666666574148081281236954964697360992431640625_real64
  real(kind=real64), parameter :: c3  = 4.1666666666666678231489839845380629412829875946044921875e-2_real64
  real(kind=real64), parameter :: c4  = 8.33333333333420578359351793551468290388584136962890625e-3_real64
  real(kind=real64), parameter :: c5  = 1.38888888888774706016626669935476456885226070880889892578125e-3_real64
  real(kind=real64), parameter :: c6  = 1.984126983838181126820060518056720866297837346792221069335937e-4_real64
  real(kind=real64), parameter :: c7  = 2.480158734051678335566677724433048979335580952465534210205078e-5_real64
  real(kind=real64), parameter :: c8  = 2.755732358282358854077720433650711129303090274333953857421875e-6_real64
  real(kind=real64), parameter :: c9  = 2.755725932837027351290253886217929135682425112463533878326416e-7_real64
  real(kind=real64), parameter :: c10 = 2.504903185531009510483624108320710455188873311271890997886658e-8_real64
  real(kind=real64), parameter :: c11 = 2.091945316636953120202141592222887245267060052356100641191006e-9_real64
  real(kind=real64), parameter :: c12 = 1.689299804577804335305572395803039298378678267908981069922447e-10_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c12
  p = c11 + x * p
  p = c10 + x * p
  p = c9 + x * p
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_estrin_12


!> Degree-11 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme to match Sollya fpminimax optimization.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_11(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0  = 1.0_real64
  real(kind=real64), parameter :: c1  = 0.5_real64
  real(kind=real64), parameter :: c2  = 0.1666666666666667684371105906393495388329029083251953125_real64
  real(kind=real64), parameter :: c3  = 4.1666666666666608842550800773096852935850620269775390625e-2_real64
  real(kind=real64), parameter :: c4  = 8.3333333333214694438328962178275105543434619903564453125e-3_real64
  real(kind=real64), parameter :: c5  = 1.38888888889239134859232560614827889367006719112396240234375e-3_real64
  real(kind=real64), parameter :: c6  = 1.984126988645421615746478050112955315853469073772430419921875e-4_real64
  real(kind=real64), parameter :: c7  = 2.480158723707121432326336285534296166588319465517997741699219e-5_real64
  real(kind=real64), parameter :: c8  = 2.755724393541999451894908626514713034794112900272011756896973e-6_real64
  real(kind=real64), parameter :: c9  = 2.75573535006209635438574308621828556908894825028255581855774e-7_real64
  real(kind=real64), parameter :: c10 = 2.510911779855824855195022427684081733900711697060614824295044e-8_real64
  real(kind=real64), parameter :: c11 = 2.08892451293866406420839973532167094250056038617913145571947e-9_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c11
  p = c10 + x * p
  p = c9 + x * p
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_11


!> Degree-11 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Estrin's scheme for instruction-level parallelism.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_estrin_11(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: x2, x4, x8
  real(kind=real64) :: b0, b1, b2, b3, b4, b5
  real(kind=real64) :: q0, q1, q2

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c(0:11) = [ &
    1.0_real64, &
    0.5_real64, &
    0.1666666666666667684371105906393495388329029083251953125_real64, &
    4.1666666666666608842550800773096852935850620269775390625e-2_real64, &
    8.3333333333214694438328962178275105543434619903564453125e-3_real64, &
    1.38888888889239134859232560614827889367006719112396240234375e-3_real64, &
    1.984126988645421615746478050112955315853469073772430419921875e-4_real64, &
    2.480158723707121432326336285534296166588319465517997741699219e-5_real64, &
    2.755724393541999451894908626514713034794112900272011756896973e-6_real64, &
    2.75573535006209635438574308621828556908894825028255581855774e-7_real64, &
    2.510911779855824855195022427684081733900711697060614824295044e-8_real64, &
    2.08892451293866406420839973532167094250056038617913145571947e-9_real64 ]

  ! Estrin's scheme: evaluate polynomial with maximum parallelism
  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  ! Pairs: c0+c1*x, c2+c3*x, c4+c5*x, c6+c7*x, c8+c9*x, c10+c11*x
  b0 = c(0) + x * c(1)
  b1 = c(2) + x * c(3)
  b2 = c(4) + x * c(5)
  b3 = c(6) + x * c(7)
  b4 = c(8) + x * c(9)
  b5 = c(10) + x * c(11)

  ! Quads
  q0 = b0 + x2 * b1   ! c0 + c1*x + c2*x^2 + c3*x^3
  q1 = b2 + x2 * b3   ! c4 + c5*x + c6*x^2 + c7*x^3
  q2 = b4 + x2 * b5   ! c8 + c9*x + c10*x^2 + c11*x^3

  ! Final assembly: 1 + x * p(x) ensures exp(0) = 1 exactly
  ! p(x) = q0 + x^4*q1 + x^8*q2
  e = 1 + x * (q0 + x4 * q1 + x8 * q2)
end function exp_remez_estrin_11


!> Degree-10 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_10(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0  = 1.0_real64
  real(kind=real64), parameter :: c1  = 0.50000000000000055511151231257827021181583404541015625_real64
  real(kind=real64), parameter :: c2  = 0.1666666666666660745477201999165117740631103515625_real64
  real(kind=real64), parameter :: c3  = 4.166666666657388440331288848028634674847126007080078125e-2_real64
  real(kind=real64), parameter :: c4  = 8.333333333377164475752607586400699801743030548095703125e-3_real64
  real(kind=real64), parameter :: c5  = 1.38888889322647565531532176663631616975180804729461669921875e-3_real64
  real(kind=real64), parameter :: c6  = 1.984126974698501824807828075591942251776345074176788330078125e-4_real64
  real(kind=real64), parameter :: c7  = 2.480150459644261430602885099006016389466822147369384765625e-5_real64
  real(kind=real64), parameter :: c8  = 2.755738179851631320657146320685093598967796424403786659240723e-6_real64
  real(kind=real64), parameter :: c9  = 2.76262647076892519593864332161370356288898619823157787322998e-7_real64
  real(kind=real64), parameter :: c10 = 2.506210200218863960750001593138364119894845316594000905752182e-8_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c10
  p = c9 + x * p
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_10


!> Degree-9 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_9(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c0 = 1.0000000000000011102230246251565404236316680908203125_real64
  real(kind=real64), parameter :: c1 = 0.49999999999999544808559903685818426311016082763671875_real64
  real(kind=real64), parameter :: c2 = 0.1666666666661721640796400834005908109247684478759765625_real64
  real(kind=real64), parameter :: c3 = 4.1666666667142381041966103794038644991815090179443359375e-2_real64
  real(kind=real64), parameter :: c4 = 8.33333336678350466986131550584104843437671661376953125e-3_real64
  real(kind=real64), parameter :: c5 = 1.38888887615885620215039342184581983019597828388214111328125e-3_real64
  real(kind=real64), parameter :: c6 = 1.98411912796172423302520915200375384301878511905670166015625e-4_real64
  real(kind=real64), parameter :: c7 = 2.480169427719168816518988118779986962181283161044120788574219e-5_real64
  real(kind=real64), parameter :: c8 = 2.763239405763823578756569615544336215862131211906671524047852e-6_real64
  real(kind=real64), parameter :: c9 = 2.75560192705068341480933351642090833877318800659850239753723e-7_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c9
  p = c8 + x * p
  p = c7 + x * p
  p = c6 + x * p
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_9


!> Degree-9 Remez minimax polynomial for exp(x) on [-ln2/2, ln2/2]
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Estrin's scheme for instruction-level parallelism.
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_estrin_9(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in [-ln2/2, ln2/2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: x2, x4, x8
  real(kind=real64) :: b0, b1, b2, b3, b4
  real(kind=real64) :: q0, q1

  ! fpminimax coefficients for (exp(x) - 1) / x on [-ln2/2, ln2/2]
  real(kind=real64), parameter :: c(0:9) = [ &
    1.0000000000000011102230246251565404236316680908203125_real64, &
    0.49999999999999544808559903685818426311016082763671875_real64, &
    0.1666666666661721640796400834005908109247684478759765625_real64, &
    4.1666666667142381041966103794038644991815090179443359375e-2_real64, &
    8.33333336678350466986131550584104843437671661376953125e-3_real64, &
    1.38888887615885620215039342184581983019597828388214111328125e-3_real64, &
    1.98411912796172423302520915200375384301878511905670166015625e-4_real64, &
    2.480169427719168816518988118779986962181283161044120788574219e-5_real64, &
    2.763239405763823578756569615544336215862131211906671524047852e-6_real64, &
    2.75560192705068341480933351642090833877318800659850239753723e-7_real64 ]

  ! Estrin's scheme: evaluate polynomial with maximum parallelism
  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  ! Pairs: c0+c1*x, c2+c3*x, c4+c5*x, c6+c7*x, c8+c9*x
  b0 = c(0) + x * c(1)
  b1 = c(2) + x * c(3)
  b2 = c(4) + x * c(5)
  b3 = c(6) + x * c(7)
  b4 = c(8) + x * c(9)

  ! Quads
  q0 = b0 + x2 * b1   ! c0 + c1*x + c2*x^2 + c3*x^3
  q1 = b2 + x2 * b3   ! c4 + c5*x + c6*x^2 + c7*x^3

  ! Final assembly: 1 + x * p(x) ensures exp(0) = 1 exactly
  ! p(x) = q0 + x^4*q1 + x^8*b4
  e = 1 + x * (q0 + x4 * q1 + x8 * b4)
end function exp_remez_estrin_9


!> Degree-6 Remez minimax polynomial for exp(x) on table-reduced interval
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax for NTABLE=32 interval.
pure function exp_remez_horner_6_n32(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in reduced interval [-ln2/(2*NTABLE), ln2/(2*NTABLE)]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on reduced interval
  real(kind=real64), parameter :: c0 = 1.0_real64
  real(kind=real64), parameter :: c1 = 0.5_real64
  real(kind=real64), parameter :: c2 = 0.1666666666666666574148081281236954964697360992431640625_real64
  real(kind=real64), parameter :: c3 = 4.1666666666487515990890955208669765852391719818115234375e-2_real64
  real(kind=real64), parameter :: c4 = 8.33333333355797990782409812027253792621195316314697265625e-3_real64
  real(kind=real64), parameter :: c5 = 1.3888932526515498401542547668441329733468592166900634765625e-3_real64
  real(kind=real64), parameter :: c6 = 1.984117206235610217975734448359048656129743903875350952148437e-4_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c6
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_6_n32


!> Degree-6 Remez minimax polynomial for exp(x) on table-reduced interval
!!
!! Approximates (exp(x) - 1) / x, then computes 1 + x * p(x).
!! Uses Horner's scheme for reproducibility across compilers.
!! Coefficients generated by Sollya fpminimax for NTABLE=64 interval.
pure function exp_remez_horner_6_n64(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value in reduced interval [-ln2/(2*64), ln2/(2*64)]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on reduced interval
  real(kind=real64), parameter :: c0 = 1.0_real64
  real(kind=real64), parameter :: c1 = 0.5_real64
  real(kind=real64), parameter :: c2 = 0.1666666666666666574148081281236954964697360992431640625_real64
  real(kind=real64), parameter :: c3 = 4.1666666666655471917835029671550728380680084228515625e-2_real64
  real(kind=real64), parameter :: c4 = 8.33333333423471635248436228948776260949671268463134765625e-3_real64
  real(kind=real64), parameter :: c5 = 1.3888899797468317633131196231488502235151827335357666015625e-3_real64
  real(kind=real64), parameter :: c6 = 1.983922810051448062072798617094804285443387925624847412109375e-4_real64

  ! Horner's scheme: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c6
  p = c5 + x * p
  p = c4 + x * p
  p = c3 + x * p
  p = c2 + x * p
  p = c1 + x * p
  p = c0 + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_6_n64


pure function anint_fast(x) result(a)
  real(real64), intent(in) :: x
  real(real64) :: a

  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52

  a = (x + round_bias) - round_bias
end function anint_fast

end module exp_repro_mod
