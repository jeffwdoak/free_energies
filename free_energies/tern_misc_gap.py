#!/usr/bin/python

# tern_misc_gap.py v0.2 9-23-2012 Jeff Doak jeff.w.doak@gmail.com

import numpy as np

BoltzConst = 0.08617  # meV/K

def sub_reg_bin(L,x):
    """
    Calculates the sub-regular solution model enthalpy of a binary phase with
    parameters L and composition x.
    """
    H = x*(1.-x)*(L[0]+L[1]*(1.-2.*x))
    return H

def sub_reg_tern(L,x,y):
    """
    Calculates the sub-regular solution model enthalpy of a ternary phase with
    parameters L and compositions x,y.
    """
    H =  (1.-x-y)*x*(L[0]+L[1]*(1.-2.*x-y))
    H += (1.-x-y)*y*(L[2]+L[3]*(1.-x-2.*y))
    H += x*y*(L[4]+L[5]*(x-y))
    H += (1.-x-y)*x*y*((1.-x-y)*L[6]+x*L[7]+y*L[8])
    return H

def dGdx(L,x,y,T):
    """
    Calculates the first derivative of a ternary sub-regular free energy with
    respect to composition variable x.
    """
    val = (1.-2.*x-y)*(L[0]+L[1]*(1.-2.*x-y))-2.*L[1]*x*(1.-x-y)
    val -= y*(L[2]+L[3]*(1.-x-2.*y))+L[3]*y*(1.-x-y)
    val += y*(L[4]+L[5]*(x-y))+x*y*L[5]
    val += (1.-2.*x-y)*y*((1.-x-y)*L[6]+x*L[7]+y*L[8])+(1.-x-y)*x*y*(L[7]-L[6])
    val += BoltzConst*T*np.log(x/(1.-x-y))
    return val

def dGdy(L,x,y,T):
    """
    Calculates the first derivative of a ternary sub-regular free energy with
    respect to composition variable y.
    """
    val = (1.-2.*y-x)*(L[2]+L[3]*(1.-2.*y-x))-2.*L[3]*y*(1.-x-y)
    val -= x*(L[0]+L[1]*(1.-y-2.*x))+L[1]*x*(1.-x-y)
    val += x*(L[4]+L[5]*(x-y))-x*y*L[5]
    val += (1.-2.*y-x)*x*((1.-x-y)*L[6]+x*L[7]+y*L[8])+(1.-x-y)*x*y*(L[8]-L[6])
    val += BoltzConst*T*np.log(y/(1.-x-y))
    return val

def d2Gdx2(L,x,y,T):
    """
    Calculates the second derivative of a ternary sub-regular free energy with
    respect to composition variable x.
    """
    val = -2.*(L[0]+3.*L[1]*(1-2.*x-y)-(L[3]+L[5])*y
          +y*(L[6]+(L[8]-L[6])*y-(L[7]-L[6])*(1-y-6.*x))
                -0.5*BoltzConst*T*(1./x+1./(1.-x-y)))
    return val

def d2Gdy2(L,x,y,T):
    """
    Calculates the second derivative of a ternary sub-regular free energy with
    respect to composition variable y.
    """
    val = -2.*(L[2]+3.*L[3]*(1-2.*y-x)-(L[1]+L[5])*x
          +x*(L[6]+(L[7]-L[6])*x-(L[8]-L[6])*(1-x-6.*y))
                -0.5*BoltzConst*T*(1./y+1./(1.-x-y)))
    return val

def d2Gdxdy(L,x,y,T):
    """
    Calculates the mixed second derivative of a ternary sub-regular free energy
    with respect to composition variables x and y.
    """
    val = -L[0]-2.*L[1]*(1.-3.*x-y)-L[2]-2.*L[3]*(1.-3.*y-x)+L[4]+2.*L[5]*(x-y)
    val += (1.-2.*x-2.*y)*L[6] + x*(2.-4.*y-3.*x)*(L[7]-L[6])
    val += y*(2.-4.*x-3.*y)*(L[8]-L[6]) + BoltzConst*T/(1.-x-y)
    return val

def hessian(L,x,y,T):
    """
    Calculates the Hessian matrix of a sub-regular free energy model.
    """
    hess = np.zeros((2,2,))
    hess[0,0] = d2Gdx2(L,x,y,T)
    hess[0,1] = d2Gdxdy(L,x,y,T)
    hess[1,0] = hess[0,1]
    hess[1,1] = d2Gdy2(L,x,y,T)
    return hess

def hess_det(L,x,y,T):
    """
    Calculates the determinant of the Hessian matrix of a sub-regular free
    energy model.
    """
    det = d2Gdx2(L,x,y,T)*d2Gdy2(L,x,y,T) - d2Gdxdy(L,x,y,T)**2
    return det

def jacobian(L,x1,y1,x2,y2,T):
    """
    Calculates the Jacobian matrix of a sub-regular free energy model at a fixed
    temperature T, with two phases in coexistence with compositions {x1,y1} and
    {x2,y2}, respectively. Composition variable x1 is taken to be independent.
    """
    jac = np.zeros((3,3))
    jac[0,0] = -1.*d2Gx2(L,x2,y2,T)
    jac[0,1] = d2Gdxdy(L,x1,y1,T)
    jac[0,2] = -1.*d2Gdxdy(L,x2,y2,T)
    jac[1,0] = jac[0,2]
    jac[1,1] = d2Gdy2(L,x1,y1,T)
    jac[1,2] = -1.*d2Gdy2(x2,y2,T)
    jac[2,0] = (x2-x1)*d2Gdx2(L,x2,y2,T) + (y2-y1)*d2Gdxdy(x2,y2,T)
    jac[2,1] = dGdy(L,x1,y1,T) - dGdy(L,x2,y2,T)
    return jac

def jac_det(L,x1,y1,x2,y2,T):
    """
    Calculates the determinant of the Jacobian matrix of a sub-regular free
    energy model at a fixed temperature T, with two phases in coexistence with
    compositions {x1,y1} and {x2,y2}, respectively. Composition variable x1 is
    taken to be independent.
    """
    val = np.linalg.det(jacobian(L,x1,y1,x2,y2,T))
    return val

def bin_spinode(L,x=None,T=None):
    """
    Calculates a binary spinodal at either a fixed composition or temperature.
    If x is given, the corresponding T_spinodal is returned. If T is given, the
    two corresponding x_spinodal are returned.
    """
    from numpy.lib.scimath import power
    if x != None:
        T_sp = 2*x*(1-x)*(L[0]+3*L[1]*(1-2.*x))/BoltzConst
        return T_sp
    elif T != None:
        # Coefficients of cubic equation a*x**3+b*x**2+c*x+d
        a = 6.*L[1]
        b = -(L[0]+9.*L[1])
        c = L[0]+3.*L[1]
        d = -T*BoltzConst/2.
        discrim = 18.*a*b*c*d-4.*b**3*d+b**2*c**2-4*a*c**3-27.*a**2*d**2
        if discrim >= 0.:
            temp = 2.*b**3-9.*a*b*c+27.*a**2*d
            x1 = (power(0.5*(temp + power(-27.*a**2*discrim,1./2.)),1./3.)
                + power(0.5*(temp - power(-27.*a**2*discrim,1./2.)),1./3.)
                + b)/(-3.*a)
            x2 = (power(0.5*(temp + power(-27.*a**2*discrim,1./2.)),1./3.)
                    *(1.+np.sqrt(3)*j)
                + power(0.5*(temp - power(-27.*a**2*discrim,1./2.)),1./3.)
                    *(1.-np.sqrt(3)*j)
                - 2.*b)/(6.*a)
            x3 = (power(0.5*(temp + power(-27.*a**2*discrim,1./2.)),1./3.)
                    *(1.-np.sqrt(3)*j)
                + power(0.5*(temp - power(-27.*a**2*discrim,1./2.)),1./3.)
                    *(1.+np.sqrt(3)*j)
                - 2.*b)/(6.*a)
            x_sp = sorted([x1,x2,x3])
            return x_sp[0],x_sp[1]
        else:
            # There exists 1 real root and 2 complex roots
            print "There exists only 1 real root!"
            return

def bin_binode(L,T):
    pass

def bin_crit(L):
    """
    Calculates the critical temperature and composition of a binary sub-regular
    solution model.
    """
    # Taking the minus root of the quadratic equation
    xc = 0.5 + (L[0] - np.sqrt(L[0]**2+27.*L[1]**2))/(18.*L[1])
    Tc = 2*xc*(1-xc)*(L[0]+3*L[1]*(1-2.*xc))/BoltzConst
    return xc,Tc
    
def tern_crit(L,x,y,T):
    pass

def main(args):
    # Read in model parameters
    # Order phases by decreasing binary citical temperatures
    pass

if __name__=="__main__":
    import sys
    main(sys.argv[1:])
