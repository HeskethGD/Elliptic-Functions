from sympy import *
from sympy.core.relational import Equality

(z, g2, g3, n, m) = symbols('z, g2, g3, n, m')
pw = Function('pw') # Weierstrass P function
pwp = Function('pwp') # Derivative of Weierstrass P function


def nth_diff_pw(nth_order: int) -> Equality:
    
    """
    Builds the symbolic equation for the nth_order derivative of the Weirstrass P function.
    Even orders: result is an O(n/2 + 1) polynomial in the Weierstrass P function.
    Odd orders: result is an O((n-1)/2) polynomial in the Weierstrass P function, multiplied by P Prime.
    
    Args: nth_order (int > 0)
    
    Returns: (sympy.core.relational.Equality)
    
    Params: g2, g3, Weierstrass elliptic invariants
    """
    
    if nth_order < 1:
        print("nth_order must not be less than 1")
        raise
        
    dpw_1 = Eq(diff(pw(z,g2,g3),(z,1)), pwp(z,g2,g3))
    
    if nth_order == 1:
        return dpw_1
        
    dpw_2 = Eq(diff(pw(z,g2,g3),(z,2)), -g2/2 + 6*pw(z,g2,g3)**2)
    
    if nth_order == 2:
        return dpw_2
    
    # Calculate all order diffs from the second order diff
    nth_order_eqs = [Eq(diff(dpw_2.lhs,(z,n)), diff(dpw_2.rhs,(z,n))).doit() for n in range(1, nth_order - 1)]
    
    # Calculate higher orders of pwp using the defining diff eqn
    # (Probably calculating more than it needs but does the job)
    dpw_1_sqrd = Eq(diff(pw(z,g2,g3),z)**2, 4*pw(z,g2,g3)**3 - g2*pw(z,g2,g3) - g3)
    
    nth_order_pwp_sqrd_eqs = [
        Eq(pwp(z,g2,g3)**m, dpw_1_sqrd.rhs**(Rational(m - (m % 2),2)) * pwp(z,g2,g3)**(m % 2)) 
        for m in range(2, nth_order - 1)
    ]
    
    # Recursively substitute lower order diffs and powers of pwp
    nth_order_diff_eq = Eq(
        nth_order_eqs[-1].lhs, 
        nth_order_eqs[-1].rhs
        .subs([_eq.args for _eq in nth_order_eqs[::-1]])
        .subs([dpw_2.args, dpw_1.args])
        .expand()
        .subs([_eq.args for _eq in nth_order_pwp_sqrd_eqs])
        .expand()
        .factor()
    )
        
    return nth_order_diff_eq