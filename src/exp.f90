integer, parameter :: n = 5
  ! Number of interpolation points
real, parameter :: bound = 0.5 * log(2.)
  ! An obscure number denoting the ratio of diameter to circumference
real :: nodes(n)
  ! Initial Chebyshev polynomial nodes
real :: exp_nodes(n)
  ! Evaluation of exp() on nodes(:) 
real :: chebyshev_weights

! Construct Chebyshev nodes
nodes(:) = bound * [(cosd((90. * (2 * k - 1)) / n), k=n,1,-1)]

! Evaluate on the nodes
! TODO: How to make these exactly rounded?
exp_nodes(:) = exp(nodes)

! Check output
print *, "Nodes", "exp(Nodes)"
print *, "-----"
print '(2(f19.16))', -bound, exp(-bound), (nodes(i), exp_nodes(i), i=1,n), bound, exp(bound)
print *

end
