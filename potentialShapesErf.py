import numpy as np
from numpy.polynomial.polynomial import Polynomial
# import scipy
from dataclasses import dataclass, field
from typing import Sequence
import matplotlib.pyplot as plt

from scipy.special import erf, erfc


@dataclass
class Error:
    """
    Error function with definable features including the center, transition sharpness, 
    Z offset, scale, and complement terms in terms of r_sq
    """
    xCenter: float
    sharpness: float
    zOffset: float
    scale: float
    complement: bool

    # setup added init functions
    # def __post_init__(self) -> None:
    #     self.polynomial = Polynomial(np.array(self.coeffs) * self.U0)
    #     self.deriv = self.polynomial.deriv()

    def __call__(self, r_sq: np.ndarray):
        r = np.sqrt(r_sq)
        E = erf(self.sharpness*(r-self.xCenter))*self.scale
        # plt.plot(x,(erf(sharpness*(x-xCenter))*scale+zOffset))
        if self.complement:
            E = 1-E
        E += self.zOffset
        
        dE_dx = 2/np.sqrt(np.pi)*np.exp(-(self.sharpness*(r-self.xCenter))**2)
        if self.complement:
            dE_dx *= -1
        return E, dE_dx

        # x = r_sq/sigma_sq
        # gauss = np.exp(-0.5*x)
        # poly =  self.polynomial(x)
        # dpoly_dx = self.deriv(x)
        # E = gauss * poly
        # dE_dx =  gauss * (-0.5 * poly + dpoly_dx)
        # return E, dE_dx / sigma_sq


def main():
    
    x = np.linspace(0, 15)
    # # Type 1
    # plt.plot(x,(erf(2*(x-10))*0.5+0.5))
    # # d/dx
    # plt.plot(x,2/np.sqrt(np.pi)*np.exp(-(2*(x-10))**2))


    xCenter = 10
    sharpness = 2
    zOffset = 0.5
    scale = 0.5
    complement = False
    # plt.plot(x,(erf(sharpness*(x-xCenter))*scale+zOffset))
    # plt.plot(x,2/np.sqrt(np.pi)*np.exp(-(sharpness*(x-xCenter))**2))
    # plt.show()
    test = Error(xCenter, sharpness, zOffset, scale, complement)
    
    E, dE_dx = test(x**2)
    plt.plot(x, E)
    plt.plot(x, dE_dx)
    plt.show()

    # # Type 2
    # plt.plot(x,(1-erf(2*(x-3)))*0.5)
    # # plt.plot(x,(erfc(2*(x-3)))*0.5)
    # # d/dx
    # plt.plot(x,-2/np.sqrt(np.pi)*0.5*np.exp(-(2*(x-3))**2))

    xCenter = 3
    sharpness = 2
    zOffset = -0.5
    scale = 0.5
    complement = True

    # plt.plot(x,(1-erf(sharpness*(x-xCenter))*scale+zOffset))
    # plt.plot(x,-2/np.sqrt(np.pi)*np.exp(-(sharpness*(x-xCenter))**2))


    test = Error(xCenter, sharpness, zOffset, scale, complement)
    
    E, dE_dx = test(x**2)
    plt.plot(x, E)
    plt.plot(x, dE_dx)
    plt.show()
    
    '''

    # Equation was A exp(-x/2) (1 + Bx)
    # Algorithm:
    # Pick B using randn(), would be positive or negative
    # If B < 0, set A = 1
    # If B > 1/2, set A = -1
    # Otherwise, pick A randomly from +/- 1
    # In mdext, only use lambda >= 0
    # Picking a random sigma
    # coeff = [A, A*B]

    L = 10.
    U0 = 1.  # lambda equivalent to scale all
    np.random.seed(1)
    # sigmaScale = np.random.uniform(10,30)
    # sigma=L/10 # /10 to make sure it goes to zero around half of box by observation
    # sigma=L/sigmaScale # /10 to make sure it goes to zero around half of box by observation
    
    sigma = np.random.uniform(0.5, 3)

    B = np.random.randn() 
    if B > 0.5:
        A = -1.
    elif B < 0.:
        A = 1.
    else:
        A = np.random.choice([1,-1]) 

    coeffs = A*np.array([1, B])

    # TODO lambda only +



    # test = GaussianPolynomial(U0,sigma)
    test = Gaussian(U0, sigma, coeffs)


    dr = 0.1
    r = np.arange(start=0.,stop=L,step=dr)
    r_sq = r **2
    E, r_sq_grad = test(r_sq)
    print("Integral(E^2):", ((E**2).sum() - 0.5*E[0]**2) * dr)
    print("Expected:", sigma * np.sqrt(np.pi) * 0.5) # 1/2 box in

    plt.plot(r,E,label='E')
    # plt.plot(r,r_sq_grad,label='grad')
    # plt.plot(r,r_sq_grad*(-2)*r,label='F')
    plt.xlabel('r')
    plt.legend()
    '''

if __name__ == "__main__":
    main()
