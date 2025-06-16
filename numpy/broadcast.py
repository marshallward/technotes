import numpy as np

# Broadcast rules are frustrating, but this doc tries to summarize them.

np.random.seed(0)

a = np.random.rand(5,3)
b = np.random.rand(1,3)
c = np.random.rand(3)

print("a=")
print(a)
print()
print("b=")
print(b)

# Lesson 1:
#   If all dimension match then operatons are per-element

print("a + a")
print(a + a)
print("shape(a + a) = {}".format(np.shape(a+a)))

# Lesson 2:
#   If number of dimensions match
#   ... and one of the dimension sizes is one
#   ... then it is *broadcast* to all elements in the other element.

print()
print("a + b = ?")
print(a + b)
print("shape(a + b) = {}".format(np.shape(a*b)))

# Lesson 2.5
#   If one has fewer dimensions and its inner dimensions match the inner
#   dimensions of the other, then it will be promoted to (1, [inner]).
#
# eg: (5,3) + (3,) => (5,3) + (1,3) => broadcast over all 5 elements.

print()
print("a + c = ?")
print(a + c)
print("shape(a + c) = {}".format(np.shape(a + c)))


# Row swapping
# TODO
#   but like this:  a[[2,1]] = a[[1,2]]
