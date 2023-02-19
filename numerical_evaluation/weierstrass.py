# copied from https://github.com/stla/pyweierstrass/blob/main/pyweierstrass/weierstrass.py

from mpmath import jtheta, pi, exp, sqrt, agm, qfrom, mpc, elliprf, im, gamma
from math import isinf

def kleinjinv(z):
    """
    Inverse of the klein-j invariant function (the one with the 1728 factor).
    
    Parameters
    ----------
    z : complex
        A complex number.
        
    Returns
    -------
    complex
        The value of the inverse of the klein-j invariant function at `z`.
        
    """
    if (isinf(z.real) or isinf(z.imag)):
        x = 0
    elif z == 0:
        return mpc(0.5, sqrt(3)/2)
    else:
        z2 = z * z
        z3 = z2 * z
        t = (-z3 + 2304 * z2 + 12288
          * sqrt(3 * (1728 * z2 - z3)) - 884736 * z)**(1/3)
        x = 1/768 * t - (1536 * z - z2) / (768 * t) + (1 - z/768)
    lbd = -(-1 - sqrt(1 - 4*x)) / 2
    return 1j * agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd))

#x = 1/768 (-z^3 + 2304 z^2 + 12288 sqrt(3) sqrt(1728 z^2 - z^3) - 884736 z)^(1/3) - (1536 z - z^2)/(768 (-z^3 + 2304 z^2 + 12288 sqrt(3) sqrt(1728 z^2 - z^3) - 884736 z)^(1/3)) + (768 - z)/768

def eisensenteinE4(tau):
    """
    Eisenstein E-series of weight 4.
    
    Parameters
    ----------
    tau : complex
        A complex number.
        
    Returns
    -------
    complex
        The value of the Eisenstein E4 series at `tau`.
        
    """
    q = qfrom(tau = tau)
    j2 = jtheta(2, 0, q)
    j3 = jtheta(3, 0, q)
    j4 = jtheta(4, 0, q)
    x2 = j2 * j2
    x2 = x2 * x2
    x3 = j3 * j3
    x3 = x3 * x3
    x4 = j4 * j4
    x4 = x4 * x4
    return (
        x2 * x2 + x3 * x3 + x4 * x4
    ) / 2


def eisensenteinE6(tau):
    """
    Eisenstein E-series of weight 6.
    
    Parameters
    ----------
    tau : complex
        A complex number.
        
    Returns
    -------
    complex
        The value of the Eisenstein E6 series at `tau`.
        
    """
    q = qfrom(tau = tau)
    j2 = jtheta(2, 0, q)
    j3 = jtheta(3, 0, q)
    j4 = jtheta(4, 0, q)
    x3 = j3**4
    x4 = j4**4
    return (
        x3**3 + x4**3 - 3 * j2**8 * (x3 + x4)
    ) / 2


def eisensteinG4(tau):
    """
    Eisenstein G-series of weight 4.
    
    Parameters
    ----------
    tau : complex
        A complex number.
        
    Returns
    -------
    complex
        The value of the Eisenstein G4 series at `tau`.
        
    """
    return pi*pi*pi*pi / 45 * eisensenteinE4(tau)


def eisensteinG6(tau):
    """
    Eisenstein G-series of weight 6.
    
    Parameters
    ----------
    tau : complex
        A complex number.
        
    Returns
    -------
    complex
        The value of the Eisenstein G6 series at `tau`.
        
    """
    return 2 * pi**6 / 945 * eisensenteinE6(tau)


def _half_periods(g2, g3):
    if g2 == 0:
        omega1 = g3^(-1/6) * gamma(1/3)**3 / 4 / pi
        tau = mpc(0.5, sqrt(3)/2)
    else:
        g2cube = g2 * g2 * g2
        j = 1728 * g2cube / (g2cube - 27*g3*g3)
        if (isinf(j.real) or isinf(j.imag)):
            return -1j*pi/2/sqrt(3), mpc("inf", "inf")
        tau = kleinjinv(j)
        omega1 = sqrt(
            7 * eisensteinG6(tau) * g2 / (12 * eisensteinG4(tau) * g3)
        )
        # t = 15 / 4 / g2 * eisensteinG4(tau)
        # omega1 = 1j * t**(0.25)
    return omega1, tau

def omega_from_g(g2, g3):
    """
    Half-periods from elliptical invariants.
    
    Parameters
    ----------
    g2 : complex
        A complex number.
    g3 : complex
        A complex number.
        
    Returns
    -------
    complex pair
        The two half-periods corresponding to the elliptical invariants `g2` and `g3`.
        
    """
    omega1, tau = _half_periods(g2, g3)
    return omega1, tau * omega1

def _g_from_omega1_and_tau(omega1, tau):
    q = qfrom(tau = tau)
    j2 = jtheta(2, 0, q)
    j3 = jtheta(3, 0, q)
    g2 = 4/3 * (pi/2/omega1)**4 * (j2**8 - (j2*j3)**4 + j3**8)
    g3 = 8/27 * (pi/2/omega1)**6 * (j2**12 - (
            (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8)
        ) + j3**12)
    return g2, g3

def g_from_omega(omega1, omega2):
    """
    Elliptical invariants from half periods.
    
    Parameters
    ----------
    omega1 : complex
        A complex number.
    omega2 : complex
        A complex number.
        
    Returns
    -------
    complex pair
        The elliptical invariants `g2` and `g3`.
        
    """
    g2, g3 = _g_from_omega1_and_tau(omega1, omega2/omega1)
    return g2, g3

# Weierstrass p-function
def _wp_from_tau(z, tau):
    q = qfrom(tau = tau)
    j2 = jtheta(2, 0, q)
    j3 = jtheta(3, 0, q)
    return (
        (pi * j2 * j3 * jtheta(4, pi*z, q) / jtheta(1, pi*z, q))**2 -
          pi**2 * (j2**4 + j3**4) / 3
    )

# Weierstrass p-function
def wp(z, omega):
    """
    Weierstrass p-function.
    
    Parameters
    ----------
    z : complex
        A complex number.
    omega : complex pair
        A pair of complex numbers, the half-periods.
        
    Returns
    -------
    complex 
        The Weierstrass p-function evaluated at `z`.
        
    """
    omega1, omega2 = omega
    if im(omega2/omega1) <= 0:
        raise Exception("The imaginary part of the periods ratio is negative.")
    return _wp_from_tau(z/omega1/2, omega2/omega1) / omega1 / omega1 / 4

def _wpprime_from_omega1_and_tau(z, omega1, tau):
    q = qfrom(tau = tau)
    w1 = 2 * omega1 / pi
    z1 = - z / 2 / omega1 * pi
    j10 = jtheta(1, 0, q, 1)
    j1z1 = jtheta(1, z1, q)
    f = j10*j10*j10 / (j1z1*j1z1*j1z1 * 
        jtheta(2, 0, q) * jtheta(3, 0, q) * jtheta(4, 0, q))
    return (
        2/(w1*w1*w1) * f * 
        jtheta(2, z1, q) * jtheta(3, z1, q) * jtheta(4, z1, q)
    )

# Weierstrass p-function prime
def wpprime(z, omega):
    """
    Derivative of Weierstrass p-function.
    
    Parameters
    ----------
    z : complex
        A complex number.
    omega : complex pair
        A pair of complex numbers, the half-periods.
        
    Returns
    -------
    complex 
        The derivative of the Weierstrass p-function evaluated at `z`.
        
    """
    omega1, omega2 = omega
    if im(omega2/omega1) <= 0:
        raise Exception("The imaginary part of the periods ratio is negative.")
    return _wpprime_from_omega1_and_tau(z, omega1, omega2 / omega1)


# Weierstrass sigma function
def wsigma(z, omega):
    """
    Weierstrass sigma function.
    
    Parameters
    ----------
    z : complex
        A complex number.
    omega : complex pair
        A pair of complex numbers, the half-periods.
        
    Returns
    -------
    complex 
        The Weierstrass sigma function evaluated at `z`.
        
    """
    omega1, omega2 = omega
    w1 = -2 * omega1 / pi
    tau = omega2 / omega1
    if im(tau) <= 0:
        raise Exception("The imaginary part of the periods ratio is negative.")
    q = qfrom(tau = tau)
    f = jtheta(1, 0, q, 1)
    h = -pi / 6 / w1 * jtheta(1, 0, q, 3) / f
    return w1 * exp(h*z*z/w1/pi) * jtheta(1, z/w1, q) / f
  
# Weierstrass zeta function
def wzeta(z, omega):
    """
    Weierstrass zeta function.
    
    Parameters
    ----------
    z : complex
        A complex number.
    omega : complex pair
        A pair of complex numbers, the half-periods.
        
    Returns
    -------
    complex 
        The Weierstrass zeta function evaluated at `z`.
        
    """
    omega1, omega2 = omega
    if im(omega2/omega1) <= 0:
        raise Exception("The imaginary part of the periods ratio is negative.")
    w1 = - omega1 / pi
    tau = omega2 / omega1
    q = qfrom(tau = tau)
    p = 1 / 2 / w1
    eta1 = p / 6 / w1 * jtheta(1, 0, q, 3) / jtheta(1, 0, q, 1)
    return -eta1 * z + p * jtheta(1, p*z, q, 1) / jtheta(1, p*z, q)

# inverse Weierstrass p-function
def invwp(w, omega):
    """
    Inverse Weierstrass p-function.
    
    Parameters
    ----------
    w : complex
        A complex number.
    omega : complex pair
        A pair of complex numbers, the half-periods.
        
    Returns
    -------
    complex 
        The inverse Weierstrass p-function evaluated at `w`.
        
    """
    omega1, omega2 = omega
    if im(omega2/omega1) <= 0:
        raise Exception("The imaginary part of the periods ratio is negative.")
    e1 = wp(omega1, omega)
    e2 = wp(omega2, omega)
    e3 = wp(-omega1-omega2, omega)
    return elliprf(w-e1, w-e2, w-e3)