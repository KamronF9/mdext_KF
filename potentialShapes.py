import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt

@dataclass
class GaussianPolynomial:
    """potential of gaussian with peak `U0` and width `sigma` and 
    polynomial with `a0-2` in terms of r_sq/sigma_sq
    """
    U0: float  #: Strength
    sigma: float  #: width of Gaussian
    polyCoeffs: list #: list of floats for polynomial coefficients
    # a0: float #: +/- nondimensional polynomial coeff a0
    # a1: float #: +/- nondimensional polynomial coeff a1
    # a2: float #: +/- nondimensional polynomial coeff a2
    

    def __call__(self, r_sq: np.ndarray):
        sigma_sq = self.sigma**2
        E1 = np.exp(-0.5*r_sq/sigma_sq)
        # E2 = a0 + a1*r_sq/sigma_sq + a2*r_sq**2/sigma_sq**2
        
        # HACK
        E_norm = 1
        if polyCoeffs == [0]*len(polyCoeffs):
            # Gaussian only
            E = self.U0 * E1
            # E_norm = np.linalg.norm(E)
            E *= 1/E_norm
            r_sq_grad = -0.5/sigma_sq*E # dE/dr_sq
        else:
            # Gaussian and polynomial
            polyEq = np.polynomial.polynomial.Polynomial(polyCoeffs, domain=None, window=None)
            E2 = polyEq(r_sq/sigma_sq)
            polyEqDeriv = polyEq.deriv()(r_sq/sigma_sq)
            E = self.U0 * E1 * E2
            # E_norm = np.linalg.norm(E)
            E *= 1/E_norm
            r_sq_grad = self.U0/E_norm * ( -0.5/sigma_sq*E1*E2 + E1*polyEqDeriv )   # dE/dr_sq

        return E, r_sq_grad


L = 10.
U0=1.  # lambda equivalent to scale all
sigmaScale = np.random.uniform(10,30)
# sigma=L/10 # /10 to make sure it goes to zero around half of box by observation
sigma=L/sigmaScale # /10 to make sure it goes to zero around half of box by observation

polyCoeffs = np.random.uniform(-10.,10.,size=(3)).tolist() # gauss and poly
# polyCoeffs = [0.] # gaussian only



test = GaussianPolynomial(U0,sigma,polyCoeffs)


r = np.arange(start=0.,stop=L,step=0.1)
r_sq = r **2
E, r_sq_grad = test(r_sq)

plt.plot(r,E,label='E')
plt.plot(r,r_sq_grad,label='grad')
# plt.plot(r_sq,E,label='E')
# plt.plot(r_sq,r_sq_grad,label='grad')
plt.legend()
# E_norm = E/np.linalg.norm(E)
# (E_norm**2).sum() # is 1 OK

# plt.plot(E_norm)

# normalize power spectrum to 1 = int E^2
# such 

# from hardrods
# def _get_random(grid1d: Grid1D, *, sigma: float, seed: int) -> torch.Tensor:
#     # Determine zero and Nyquist frequency weights / real constraints:
#     iGz = grid1d.iGz
#     Nz = grid1d.grid.shape[2]
#     is_real = torch.logical_or(iGz == 0, 2 * iGz == Nz)
#     Gweight = torch.where(is_real, 1.0, 2.0)
#     # Create white noise with above constraints:
#     Gmag = grid1d.Gmag
#     torch.manual_seed(seed)
#     Gnoise = torch.randn_like(Gmag, dtype=torch.complex128)
#     Gnoise[is_real] = Gnoise[is_real].real.to(torch.complex128)
#     # Filter and normalize:
#     Gnoise *= np.exp(-0.5 * (Gmag * sigma).square())
#     Gnoise *= np.sqrt(1.0 / (qp.utils.abs_squared(Gnoise) * Gweight).sum())
#     return Gnoise