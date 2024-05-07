import numpy as np
from dataclasses import dataclass
from scipy.special import sinc
import matplotlib.pyplot as plt


@dataclass
class Cosine:
    """
    cosine function with definable features including the wavelength and strength
    in terms of r_sq
    """
    halfL: float  #: half-length of box, used to set the wavelength
    U0: float #: strength of signal

    def __call__(self, r_sq: np.ndarray):
        x = r_sq
        k = np.pi / self.halfL
        E = self.U0 * np.cos(np.sqrt(x) * k)
        dE_dx = -self.U0 * sinc(np.sqrt(x)/self.halfL)  * 0.5 * k**2
        return E, dE_dx


def main():
    
    halfL = 10.0
    r = np.linspace(0, halfL,num=1000)
    r_sq = r**2

    for U0 in [2.0]:
        print(U0)
        test = Cosine(halfL, U0)       
        E, dE_dx = test(r_sq)
        plt.plot(r, E)
        plt.plot(r, dE_dx)
        plt.plot(np.sqrt(0.5*(r_sq[:-1] + r_sq[1:])), np.diff(E) / np.diff(r_sq))
        plt.show()
    
if __name__ == "__main__":
    main()