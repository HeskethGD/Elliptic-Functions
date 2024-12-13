from sympy import *
from sympy.core.relational import Equality

(z, g2, g3, n, m) = symbols('z, g2, g3, n, m')
pw = Function('pw') # Weierstrass P function
pwp = Function('pwp') # Derivative of Weierstrass P function


def nth_diff_pw(nth_order: int) -> Equality:
    
    """
    Builds the symbolic equation for the nth_order derivative of the Weierstrass P function.
    Odd orders: result is an O((n-1)/2) polynomial in the Weierstrass P function, multiplied by P Prime.
    Even orders: result is an O(n/2 + 1) polynomial in the Weierstrass P function.
    
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

def nth_diff_zeta(nth_order: int) -> Equality:
    
    """
    Builds the symbolic equation for the nth_order derivative of the Weierstrass zeta function.
    Odd orders: result is an O((n-1)/2 + 1) polynomial in the Weierstrass P function.
    Even orders: result is an O(n/2 - 1) polynomial in the Weierstrass P function, multiplied by P Prime.
    
    Args: nth_order (int > 0)
    
    Returns: (sympy.core.relational.Equality)
    
    Params: g2, g3, Weierstrass elliptic invariants
    """
    if nth_order < 1:
        print("nth_order must not be less than 1")
        raise
        
    dzw_pw = Eq(diff(zw(z, g2, g3), (z,nth_order)), -diff(pw(z, g2, g3),(z,nth_order-1)))
    if nth_order == 1:
        return dzw_pw
    nth_order_diff_eq = dzw_pw.subs(*nth_diff_pw(nth_order-1).args)
        
    return nth_order_diff_eq

def nth_diff_sigma(nth_order: int) -> Equality:
    
    """
    Builds the symbolic equation for the nth_order derivative of the Weierstrass sigma function.
    The result is a polynomial in the Weierstrass P, P Prime, zeta, and sigma functions.
    
    Args: nth_order (int > 0)
    
    Returns: (sympy.core.relational.Equality)
    
    Params: g2, g3, Weierstrass elliptic invariants
    """
    
    if nth_order < 1:
        print("nth_order must not be less than 1")
        raise

        
    def _nth_diff_sigma_nth_diff_zw(nth_order: int) -> Equality:
        if nth_order < 1:
            print("nth_order must not be less than 1")
            raise
        return Eq(diff(sigma(z, g2, g3), (z,nth_order)), 
                  diff(sigma(z, g2, g3)*zw(z, g2, g3), (z,nth_order-1)).doit()
                 )
    
    if nth_order == 1:
        return _nth_diff_sigma_nth_diff_zw(nth_order)
        
    # Calculate all order diff sigma in terms of diff zw
    all_orders_dsigma_dzw = [_nth_diff_sigma_nth_diff_zw(_n) for _n in range(1, nth_order + 1)]
    nth_order_dsigma_dzw = all_orders_dsigma_dzw[-1]
    for nth_dsigma_dzw in all_orders_dsigma_dzw[0:-1][::-1]:
        nth_order_dsigma_dzw = Eq(nth_order_dsigma_dzw.lhs,
                                  nth_order_dsigma_dzw.rhs.subs(*nth_dsigma_dzw.args))
        
    # Substitute for derivatives of Weierstrass zeta in terms of polynomials in Weierstrass P and P Prime
    all_orders_dzw_dpw_args = [nth_diff_zeta(_n).args for _n in range(1, nth_order)][::-1]
    nth_order_dsigma_dzw = nth_order_dsigma_dzw.subs(all_orders_dzw_dpw_args).expand()
        
    return nth_order_dsigma_dzw