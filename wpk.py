from sympy import *
pw = Function('pw')
pwp = Function('pwp')
(z, g2, g3) = symbols('z, g2, g3')

def wpk(k: int) -> Equality:
    """
    Symbolic equation for kth order derivative of Weierstrass P function in terms of Weierstrass P and Weierstrass P Prime.
    Even orders: result is an O(k/2 + 1) polynomial in the Weierstrass P function.
    Odd orders: result is an O((k-1)/2) polynomial in the Weierstrass P function, multiplied by Weierstrass P Prime.
    
    Params
    - pw (Function): Weierstrass P function.
    - pwp (Function): First order derivative of Weierstrass P function (Weierstrass P Prime).
    - z (Symbol): Argument of function.
    - g2 (Symbol): Elliptic invariant.
    - g3 (Symbol): Elliptic invariant.
    
    Args
    - k (int): Order of derivative (k > 0).
    
    Returns
    - wpk_equation (Equality): Sympy equation for kth order derivative of Weierstrass P function.
    """
    
    if k < 1:
        print('k must be greater than 0')
        raise
        
    wp_1 = Eq(diff(pw(z,g2,g3),z), pwp(z,g2,g3))
    if k == 1:
        return wp_1
   
    wp_2 = Eq(diff(pw(z, g2, g3), (z, 2)), -g2/2 + 6*pw(z, g2, g3)**2)
    if k == 2: 
        return wp_2
    
    # Recursively calculate the higher order derivative equations and express them in pw and pwp
    wp_1_sqrd = Eq((diff(pw(z,g2,g3),z))**2, 4*pw(z,g2,g3)**3 - g2 * pw(z,g2,g3) - g3)
    wps = [wp_1_sqrd, wp_2]
    for j in range(2, k):
        wps.append(Eq(wps[j-1].lhs.diff(z), wps[j-1].rhs.diff(z).subs([eq.args for eq in wps])).expand())
    
    wpk_equation = Eq(wps[-1].lhs, wps[-1].rhs.subs(*wp_1.args).factor())
    
    return wpk_equation


def run_tests():
    # Checks that the first few orders return the expected equations

    defined_wps = [
        Eq(diff(pw(z,g2,g3),z), pwp(z,g2,g3)),
        Eq(diff(pw(z, g2, g3), (z, 2)), -g2/2 + 6*pw(z, g2, g3)**2),
        Eq(diff(pw(z, g2, g3), (z, 3)), 12*pw(z, g2, g3)*pwp(z, g2, g3)),
        Eq(diff(pw(z, g2, g3), (z, 4)), -18*g2*pw(z, g2, g3) - 12*g3 + 120*pw(z, g2, g3)**3),
        Eq(diff(pw(z, g2, g3), (z, 5)), -18*g2*pwp(z, g2, g3) + 360*pw(z, g2, g3)**2*pwp(z, g2, g3)),
        Eq(diff(pw(z, g2, g3), (z, 6)), 9*g2**2 - 1008*g2*pw(z, g2, g3)**2 - 720*g3*pw(z, g2, g3) + 5040*pw(z, g2, g3)**4)
    ]
    
    for k in range(len(defined_wps)):
        print(defined_wps[k] == wpk(k + 1).expand())