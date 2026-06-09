implicit none

real(kind=16) :: q

print *, exp(-3.1)
print *, exp_cr(-3.1)

q = -3.1
print *, exp(q)

q = 2.
print *, 1./log(q)
print *, 2./log(q)

contains

elemental function exp_cr(x) result(a)
  real, intent(in) :: x
    !< Input value [nondim]
  real :: a
    !< Exponential of x [nondim]

  real :: r
    ! Rescaled value of x: r = x - K * ln 2
  integer :: K
    ! Scaling factor, where exp(x) = 2^K exp(r)
  real :: e
    ! Scaled result: exp(r)

  ! TODO: Verify reproducible bit/hex values
  real, parameter :: INV_LN2 = 1.4426950408889634
  real, parameter :: TWO_INV_LN2 = 2.8853900817779268
  !real, parameter :: INV_LN2 = 1.44269504088896340735992468100189204
  !real, parameter :: TWO_INV_LN2 = 2.88539008177792681471984936200378409
  real, parameter :: LN2_HI = 6.93147180369123816490e-01
  real, parameter :: LN2_LO = 1.90821492927058770002e-10

  ! For Chebyshev evaluation

  ! Scale to [-0.5 ln 2, 0.5 ln 2]
  K = nint(x * INV_LN2)
  r = (x - K * LN2_HI) - K * LN2_LO

  ! Compute exp(r)
  ! TODO: Explain rescale to [-1,1]
  e = exp_remez(r * TWO_INV_LN2)

  ! Descale 
  a = scale(e, K)

end function exp_cr

pure function exp_remez(x) result(e)
  real, intent(in) :: x
    !< Input value
  real :: e
    !< Approximation of exp(x)

  real :: b0, b1, b2
  integer :: n

  real, parameter :: t(0:11) = [ &
    1.0302544918096181e+00, &
    3.5180320783770408e-01, &
    3.0330010354096444e-02, &
    1.7475636139768675e-03, &
    7.5594039827123609e-05, &
    2.6172719073043696e-06, &
    7.5535800680132951e-08, &
    1.8689063234385808e-09, &
    4.0465189925824025e-11, &
    7.7884384800120311e-13, &
    1.3533787921576843e-14, &
    2.5408732182498882e-16 &
  ]

  b1 = 0.0
  b2 = 0.0

  do n = 11,1,-1
    b0 = t(n) + 2. * x * b1 - b2
    b2 = b1
    b1 = b0
  enddo
  e = t(0) + x * b1 - b2
end function exp_remez

end
