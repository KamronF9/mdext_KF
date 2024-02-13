import numpy as np
from numpy.polynomial.polynomial import Polynomial
import scipy
from dataclasses import dataclass, field
from typing import Sequence
import matplotlib.pyplot as plt

@dataclass
class GaussianPolynomial:
    """potential of gaussian with peak `U0` and width `sigma` and 
    polynomial with `polyCoeffs` in terms of r_sq/sigma_sq
    """
    U0: float  #: Strength
    sigma: float  #: width of Gaussian
    coeffs: Sequence[float] = (1.0,)  #: list of floats for polynomial coefficients
    polynomial: Polynomial = field(init=False)  #: polynomial including U0
    deriv: Polynomial = field(init=False)  #: corresponding derivative

    def __post_init__(self) -> None:
        self.polynomial = Polynomial(np.array(self.coeffs) * self.U0)
        self.deriv = self.polynomial.deriv()

    def __call__(self, r_sq: np.ndarray):
        sigma_sq = self.sigma**2
        x = r_sq/sigma_sq
        gauss = np.exp(-0.5*x)
        poly =  self.polynomial(x)
        dpoly_dx = self.deriv(x)
        E = gauss * poly
        dE_dx =  gauss * (-0.5 * poly + dpoly_dx)
        return E, dE_dx / sigma_sq


def main():
    L = 10.
    U0=1.  # lambda equivalent to scale all
    sigmaScale = np.random.uniform(10,30)
    # sigma=L/10 # /10 to make sure it goes to zero around half of box by observation
    sigma=L/sigmaScale # /10 to make sure it goes to zero around half of box by observation
    # np.random.seed(0)
    
    coeffs = np.random.randn(3) # gauss and poly
    powers = np.arange(len(coeffs))
    power_pair_sums = powers[:, None] + powers[None, :]
    norm_fac = np.sqrt(
        np.sqrt(np.pi)
        / (coeffs @ scipy.special.gamma(power_pair_sums + 0.5) @ coeffs)
    )
    coeffs *= norm_fac

    # test = GaussianPolynomial(U0,sigma)
    test = GaussianPolynomial(U0, sigma, coeffs)


    dr = 0.1
    r = np.arange(start=0.,stop=L,step=dr)
    r_sq = r **2
    E, r_sq_grad = test(r_sq)
    print("Integral(E^2):", ((E**2).sum() - 0.5*E[0]**2) * dr)
    print("Expected:", sigma * np.sqrt(np.pi) * 0.5) # 1/2 box in

    plt.plot(r,E,label='E')
    # plt.plot(r,r_sq_grad,label='grad')
    plt.plot(r,r_sq_grad*(-2)*r,label='F')
    plt.xlabel('r')
    plt.legend()


if __name__ == "__main__":
    main()
