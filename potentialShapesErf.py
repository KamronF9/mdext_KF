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
    correctScale: float #: scaling correction of the curve
    # complement: bool #: if the complement of the erf is to be used
    U0: float #: strength of signal

    # setup added init functions
    # def __post_init__(self) -> None:
    #     self.polynomial = Polynomial(np.array(self.coeffs) * self.U0)
    #     self.deriv = self.polynomial.deriv()

    def __call__(self, r_sq: np.ndarray):
        x = r_sq
        E = erf(self.sharpness*(x-self.xCenter**2))*self.correctScale
        E = 1-E  # complement
        E += self.zOffset
        E *= self.U0
        
        dE_dx = 2/np.sqrt(np.pi)*np.exp(-(self.sharpness*(x-self.xCenter**2))**2)*self.correctScale*self.U0*self.sharpness
        dE_dx *= -1  # complement
        return E, dE_dx



def main():
    
    r = np.linspace(0, 15,num=1000)
    r_sq = r**2
    # # Type 1 - well trap

    # plt.plot(x,(erf(2*(x-10))*0.5+0.5))
    # # d/dx
    # plt.plot(x,2/np.sqrt(np.pi)*np.exp(-(2*(x-10))**2))
    '''
    # for scale in range(4):
    xCenter = 10
    sharpness = 2
    zOffset = 0.5
    correctScale = 0.5
    complement = False

    test = Error(xCenter, sharpness, zOffset, correctScale, complement, U0)

    # test as it would be w r_sq input    
    # E, dE_dx = test(r_sq)
    # plt.plot(r_sq, E)
    # plt.plot(r_sq, dE_dx)
    # plt.show()

    E, dE_dx = test(r)
    plt.plot(r, E)
    plt.plot(r, dE_dx)
    plt.show()
    '''

    # # Type 2 - flat top center:

    # plt.plot(x,(1-erf(2*(x-3)))*0.5)
    # # plt.plot(x,(erfc(2*(x-3)))*0.5)
    # # d/dx
    # plt.plot(x,-2/np.sqrt(np.pi)*0.5*np.exp(-(2*(x-3))**2))

    U0s = np.arange(0,2*2+1)-2
    

    # for U0 in U0s:
    for U0 in [1]:
        print(U0)
        # xCenter = 5
        # sharpness = 0.2
        xCenter = 10 # 10
        sharpness = 0.1 # 0.1
        zOffset = -0.5 # fixed
        correctScale = 0.5 # fixed
        # complement = True # always
        # U0 = -1

        test = Error(xCenter, sharpness, zOffset, correctScale, U0)
        
        E, dE_dx = test(r_sq)
        plt.plot(r, E)
        plt.plot(r, dE_dx)
        plt.plot(np.sqrt(0.5*(r_sq[:-1] + r_sq[1:])), np.diff(E) / np.diff(r_sq))
        plt.show()
    
if __name__ == "__main__":
    main()
