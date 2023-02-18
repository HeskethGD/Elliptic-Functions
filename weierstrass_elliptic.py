from mpmath import jtheta, pi, exp, sqrt, polyroots, agm, log, elliprf

class WeierstrassElliptic:
    
    # This is a slightly modified version of the package developed by user stla on github
    # https://github.com/stla/pyweierstrass/blob/main/pyweierstrass/weierstrass.py
    
    def __init__(self):
        pass  

    def p_weierstrass_from_g2_g3(self, g2, g3, derivative = 0):
        # It is not clear what this algorithm is doing re the ordering of the roots
        # https://math.stackexchange.com/questions/4475194/half-periods-ratio-for-the-wp-weierstrass-function
        if derivative < 0 or derivative > 3 or not isinstance(derivative, int):
            raise ValueError("`derivative` must be an integer between 0 and 3.")
        r1, r2, r3 = polyroots([4, 0, -g2, -g3])
        e3 = r3
        a = sqrt(r1 - r3)
        b = sqrt(r1 - r2)
        c = sqrt(r2 - r3)
        w1 = None
        if abs(a + b) < abs(a - b):
            b *= -1
        if abs(a + c) < abs(a - c):
            c *= -1
        if abs(c + 1j*b) < abs(c - 1j*b):
            e3 = r1
            a = sqrt(r3 - r1)
            b = sqrt(r3 - r2)
            c = sqrt(r2 - r1)
            w1 = 1 / agm(1j*b, c) 
        else:
            w1 = 1 / agm(a, b)
        w3 = 1j / agm(a, c)
        q = exp(1j * pi * w3/w1)
        if derivative != 1:
            pw0 = lambda z: (
                e3 + (pi * jtheta(2, 0, q) * jtheta(3, 0, q) * jtheta(4, z/w1, q)
                  / (pi * w1 * jtheta(1, z/w1, q)))**2    
            )
            if derivative == 0:
                return pw0
            if derivative == 2:
                return lambda z: 6*pw0(z)**2 - g2/2
        f = (jtheta(1, 0, q, 1)**3 
            / (jtheta(2, 0, q) * jtheta(3, 0, q) * jtheta(4, 0, q)))
        pw1 = lambda z: (
            -2*(1/w1)**3 * jtheta(2, z/w1, q) * jtheta(3, z/w1, q) 
                * jtheta(4, z/w1, q) * f / jtheta(1, z/w1, q)**3
        )
        if derivative == 1:
            return pw1
        if derivative == 3:
            return lambda z: 12 * pw0(z) * pw1(z)  


#     z = 1+1j
#     g2 = 5 - 3j
#     g3 = 2 - 7j
#     #print(p_weierstrass_from_g2_g3(g2, g3)(z))
#     # (-0.117722596733725 + 0.543400126978109j)
#     g2 = -5 + 3j
#     g3 = 2 + 7j
#     #print(p_weierstrass_from_g2_g3(g2, g3)(z))
#     # (-1.36964875245851 - 2.51228593126988j)
#     z = 1+1j
#     g2 = 5 + 3j
#     g3 = 2 + 7j
#     #print(p_weierstrass_from_g2_g3(g2, g3, derivative=1)(z))
#     # (-1.78040154378359 + 0.213468471404325j)
#     #
#     #print(p_weierstrass_from_g2_g3(g2, g3, derivative=1)(z)**2)
#     #pw = p_weierstrass_from_g2_g3(g2, g3)(z)
#     #print(4 * pw**3 - g2 * pw - g3) # should be equal


    def p_weierstrass_from_w1_w2(self, w1, w2):
        if w2.imag*w1.real < w1.imag*w2.real:
            raise ValueError("No solution. Do you want to exchange `w1 and `w2`?")
        ratio = w2 / w1
        if ratio.imag <= 0:
            raise ValueError("The ratio `w2/w1` must have a positive imaginary part.")
        q = exp(1j * pi * ratio)
        j2 = jtheta(2, 0, q)
        j3 = jtheta(3, 0, q)
        g2 = 4/3 * (pi/2/w1)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
        g3 = 8/27 * (pi/2/w1)**6 * (j2**12 - (
            (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
        ) + j3**12)
        #print("*******")
        #print((g2, g3)) 
        #print("*******")
        return p_weierstrass_from_g2_g3(g2, g3)

#     z = 1+1j
#     w1 = 1.0 + 4.0j
#     w2 = 1.0 + 15.0j
#     #print(p_weierstrass_from_w1_w2(w1, w2)(z))
#     # (0.0076936368424553 - 0.498821838483149j)

    def p_weierstrass_from_tau(self, tau):
        return self.p_weierstrass_from_w1_w2(1.0, tau)

    #print(p_weierstrass_from_tau(1 + 4j)(z))
    #(-0.430565782630798 - 3.62705469588323e-16j)


    def zeta_weierstrass(self, g2, g3):
        r1, r2, r3 = polyroots([4, 0, -g2, -g3])
        a = sqrt(r1 - r3)
        b = sqrt(r1 - r2)
        c = sqrt(r2 - r3)
        w1 = None
        if abs(a + b) < abs(a - b):
            b *= -1
        if abs(a + c) < abs(a - c):
            c *= -1
        if abs(c + 1j*b) < abs(c - 1j*b):
            a = sqrt(r3 - r1)
            b = sqrt(r3 - r2)
            c = sqrt(r2 - r1)
            w1 = 1 / agm(1j*b, c) / 2
        else:
            w1 = 1 / agm(a, b) / 2
        w3 = 1j / agm(a, c) / 2
        #print("#######################################")
        #print(p_weierstrass_from_g2_g3(g2, g3)(w3)) # one of the ri
        q = exp(1j * pi * w3/w1)    
        p = 1 / w1 / 2
        eta1 = p / 6 / w1 * jtheta(1, 0, q, 3) / jtheta(1, 0, q, 1)
        return lambda z: (
            - eta1 * z 
            + p * jtheta(1, p*z, q, 1) / jtheta(1, p*z, q)
        )

#     z = 1+1j
#     g2 = 5 + 3j
#     g3 = 5 + 3j
#     print(zeta_weierstrass(g2, g3)(z)) # same as Wolfram and R
#     # (0.802084165492408 - 0.381791358666872j)


    def sigma_weierstrass(self, g2, g3):
        r1, r2, r3 = polyroots([4, 0, -g2, -g3])
        a = sqrt(r1 - r3)
        b = sqrt(r1 - r2)
        c = sqrt(r2 - r3)
        w1 = None
        if abs(a + b) < abs(a - b):
            b *= -1
        if abs(a + c) < abs(a - c):
            c *= -1
        if abs(c + 1j*b) < abs(c - 1j*b):
            a = sqrt(r3 - r1)
            b = sqrt(r3 - r2)
            c = sqrt(r2 - r1)
            w1 = 1 / agm(1j*b, c)
        else:
            w1 = 1 / agm(a, b)
        w3 = 1j / agm(a, c)
        q = exp(1j * pi * w3/w1)    
        f = jtheta(1, 0, q, 1)
        eta1 = - pi / 6 / w1 * jtheta(1, 0, q, 3) / f
        return lambda z: w1 * exp(eta1 * z*z/w1/pi) * jtheta(1, z/w1, q) / f

#     z = 1+1j
#     g2 = 5 + 3j
#     g3 = 5 - 3j
#     # print(sigma_weierstrass(g2, g3)(z))     
#     # (1.01629457740468 + 1.20442080691231j)   

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