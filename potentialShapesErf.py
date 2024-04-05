import numpy as np
from numpy.polynomial.polynomial import Polynomial
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
    xCenter: float #: center of transition in erf
    sharpness: float #: how sharp the transition is
    zOffset: float #: adjust the height offset of the curve
    scale: float #: scaling of the curve
    complement: bool #: if the complement of the erf is to be used

    # setup added init functions
    # def __post_init__(self) -> None:
    #     self.polynomial = Polynomial(np.array(self.coeffs) * self.U0)
    #     self.deriv = self.polynomial.deriv()

    def __call__(self, r_sq: np.ndarray):
        r = np.sqrt(r_sq)
        E = erf(self.sharpness*(r-self.xCenter))*self.scale
        if self.complement:
            E = 1-E
        E += self.zOffset
        
        dE_dx = 2/np.sqrt(np.pi)*np.exp(-(self.sharpness*(r-self.xCenter))**2)
        if self.complement:
            dE_dx *= -1
        return E, dE_dx



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

    test = Error(xCenter, sharpness, zOffset, scale, complement)
    
    E, dE_dx = test(x**2)
    plt.plot(x, E)
    plt.plot(x, dE_dx)
    plt.show()
    
if __name__ == "__main__":
    main()
