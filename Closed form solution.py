import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import norm
from scipy.optimize import newton
import scipy.optimize as opt
import math
from mpl_toolkits.mplot3d import Axes3D




# Heston model parameters
S0 = 100  # Initial stock price
K = 100   # Strike price
r = 0.05  # Risk-free interest rate
T = 1.0   # Time to maturity
v0 = 0.04 # Initial volatility
kappa = 4.0
theta = 0.04
eta = 0.2
rho = -0.7

# Numerical parameters
x_max = 1.5  # Max x value for integration
N = 100      # Number of intervals for numerical integration




# Define the characteristic function of the Heston model
def characteristic_function(u, x, v, tau):
    alpha = -u**2 / 2 - 1j * u / 2
    beta = kappa - rho * eta * 1j * u
    gamma = 0.5 * eta**2
    
    d = np.sqrt(beta**2 - 4 * alpha * gamma)
    r_plus = (beta + d) / eta**2
    r_minus = (beta - d) / eta**2
    g = r_minus / r_plus
    
    C_tau = (r_minus * tau - 2 / eta**2) * np.log((1 - (1 / g) * np.exp(-d * tau)) / (1 - g))
    D_tau = (r_minus / (1 - np.exp(-d * tau))) * (1 - (1 / g) * np.exp(-d * tau))
    
    return np.exp(1j * u * x) * np.exp(C_tau * v + D_tau * v)

# Define the numerical integration for P_j
def integrate_Pj(j, x, v, tau):
    integral = 0.0
    for n in range(N):
        u_n = -x_max + n * 2 * x_max / N
        integrand = (1 / (2 * np.pi)) * np.exp(-1j * u_n * x) * characteristic_function(u_n, x, v, tau) * (1j * u_n / 2)
        integral += integrand.real  # We take the real part
    return integral

# Calculate pseudo-probabilities P0 and P1 for European call option
def calculate_pseudo_probabilities(x, v, tau):
    P0 = 0.5 + integrate_Pj(0, x, v, tau)
    P1 = 0.5 + integrate_Pj(1, x, v, tau)
    return P0, P1

# Example usage
x_val = np.log(S0 / K)
P0_result, P1_result = calculate_pseudo_probabilities(x_val, v0, T)
print("P0:", P0_result)
print("P1:", P1_result)


