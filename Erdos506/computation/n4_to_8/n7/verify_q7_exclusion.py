"""Exact exclusion of the seven-four-circle orbit for seven points.

An affine normalization does not preserve the standard Euclidean quadratic
form.  We therefore retain the common positive-definite quadratic part

    Q(x,y) = x^2 + 2*k*x*y + h*y^2,    h-k^2 > 0.

The program checks the displayed algebraic reductions used in the accompanying
paper; the final contradiction is mathematical, not a numerical search.
"""

import sympy as sp


a, b, h, k = sp.symbols("a b h k")
t = (h - a) / (a - 1)
d_plus = t**2 + 2 * k * t + h
d_minus = t**2 - 2 * k * t + h
y3 = (t + h) / d_plus
y4 = (a - t) / d_minus
x3 = t * y3
x4 = -t * y4


def q(x, y):
    return x**2 + 2 * k * x * y + h * y**2


# After the circle 1256 gives a=b*h, subtracting the equations of 0123
# and 0356, and of 0145 and 0246, gives x3=t*y3 and x4=-t*y4.
# The first two circle equations then give the displayed y3,y4.  Verify the
# remaining two circle determinants after these substitutions.
e1346 = (
    -a*x3*y4 + a*x4*y3 - a*y3 + a*y4
    + h*y3**2*y4 - h*y3*y4**2
    + 2*k*x3*y3*y4 - 2*k*x4*y3*y4
    + x3**2*y4 - x3*y4 - x4**2*y3 + x4*y3
)
e2345 = (
    b*h*x3*y4 - b*h*x3 - b*h*x4*y3 + b*h*x4
    - h*x3*y4**2 + h*x3*y4 + h*x4*y3**2 - h*x4*y3
    + 2*k*x3*x4*y3 - 2*k*x3*x4*y4
    + x3**2*x4 - x3*x4**2
)

common_denominator = (
    (a**2*h - 2*a**2*k + a**2 + 2*a*h*k - 4*a*h
     + 2*a*k + h**2 - 2*h*k + h)
    *
    (a**2*h + 2*a**2*k + a**2 - 2*a*h*k - 4*a*h
     - 2*a*k + h**2 + 2*h*k + h)
)
r1 = a*h - a*k - h*k + h
r2 = a*k - a - h + k

assert sp.factor(sp.together(e1346.subs(b, a/h))) == (
    2*a*(a - 1)**3*(a - h)*r1 / common_denominator
)
assert sp.factor(sp.together(e2345.subs(b, a/h))) == (
    -2*a*(a - 1)*(a - h)**3*r2 / common_denominator
)
assert sp.factor(r1 + r2) == (a - k) * (h - 1)
assert sp.factor(r2.subs(a, k)) == k**2 - h
assert sp.factor(t + h) == a * (h - 1) / (a - 1)

print("EXACT PASS")
print("Positive definiteness makes d_plus and d_minus positive.")
print("Distinctness gives a!=0, a!=1, b!=1, and y3!=0, hence h!=1.")
print("The last two circle equations force a=k and h=k^2, contradicting h-k^2>0.")
