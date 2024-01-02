import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
import numpy as np
from scipy.integrate import quad

# Heston SEMI ANALYTICAL SOLUTION

def PIntegrand(u, lambda_, vbar, eta, rho, v0, r, tau, S0, K, j):
    F = S0 * np.exp(r * tau)
    x = np.log(F / K)
    a = lambda_ * vbar

    if j == 1:
        b = lambda_ - rho * eta
        alpha = -u**2 / 2 - u / 2 * 1j + 1j * u
        beta = lambda_ - rho * eta - rho * eta * 1j * u
    else:
        b = lambda_
        alpha = -u**2 / 2 - u / 2 * 1j
        beta = lambda_ - rho * eta * 1j * u

    gamma = eta**2 / 2
    d = np.sqrt(beta**2 - 4 * alpha * gamma)
    rplus = (beta + d) / (2 * gamma)
    rminus = (beta - d) / (2 * gamma)
    g = rminus / rplus

    D = rminus * (1 - np.exp(-d * tau)) / (1 - g * np.exp(-d * tau))
    C = lambda_ * (rminus * tau - (2 / (eta**2)) * np.log((1 - g * np.exp(-d * tau)) / (1 - g)))

    top = np.exp(C * vbar + D * v0 + 1j * u * x)
    bottom = 1j * u
    real_part = (top / bottom).real
    return real_part


def P(lambda_, vbar, eta, rho, v0, r, tau, S0, K, j):
    value, _ = quad(PIntegrand, 0, np.inf, args=(lambda_, vbar, eta, rho, v0, r, tau, S0, K, j))
    prob = 0.5 + (1/np.pi) * value
    return prob

def HestonCallClosedForm(lambda_, vbar, eta, rho, v0, r, tau, S0, K):

    A = S0 * P(lambda_, vbar, eta, rho, v0, r, tau, S0, K, 1)
    B = K * np.exp(-r * tau) * P(lambda_, vbar, eta, rho, v0, r, tau, S0, K, 0)
    
    return A - B

##-------------------------------------------------------------------------------------------------------

# implementation for finite difference explicit scheme
def Heston_explicit_bf(kappa, eta, sigma, rho, V0, r, T, dt, S0, I, J, K):
    # setting up the directions and steps
    nt = int(np.ceil(T / dt))

    # lower bounds of s and v are both 0
    sbound = S0 * 1.2
    vbound = V0 * 1.3
    ds = S0 / I # step length of s
    dv = V0 / J # step length of v
    s_points = np.arange(0, sbound + ds, ds)
    v_points = np.arange(0, vbound + dv, dv)
    ns = len(s_points) - 1
    nv = len(v_points) - 1
    i_points = np.arange(ns + 1)
    j_points = np.arange(nv + 1)

    # Payoff of a European call option at maturity
    u_initial = np.maximum(s_points - K, 0).reshape((ns + 1, 1))
    U = np.repeat(u_initial, nv + 1, axis=1)
    U_n = copy.deepcopy(U)

    A, B, C, D, E, F = np.zeros_like(U), np.zeros_like(U), np.zeros_like(U), np.zeros_like(U), np.zeros_like(U), np.zeros_like(U)

    U_time = []
    for x in range(nt):
        for i in range(len(i_points)):
            for j in range(len(j_points)):
                
                A[i][j] = 1 - (np.square(i)* v_points[j] * dt) - (np.square(sigma)*j*dt/dv )- (r*dt)
                B[i][j] = rho*sigma*i*j*dt /4
                C[i][j] = 0.5 * (i**2 * v_points[j] - r*i) *dt
                D[i][j] = 0.5 * (i**2 * v_points[j] + r*i) *dt
                E[i][j] = (sigma**2 *j - kappa*(eta - v_points[j]))*dt / (2*dv)
                F[i][j] = (sigma**2 *j + kappa*(eta - v_points[j]))*dt / (2*dv)

        for i in range(1, len(i_points)-1):
            for j in range(1, len(j_points)-1):
                U_n[i][j] = A[i][j]*U[i][j] + B[i][j]*(U[i-1][j-1] - U[i-1][j+1] - U[i+1][j-1] + U[i+1][j+1]) + C[i][j]*U[i-1][j] + D[i][j]*U[i+1][j] + E[i][j]*U[i][j-1] + F[i][j]*U[i][j+1]

        for j in range(len(j_points)):
            U_n[ns][j] = ( -U_n[ns-2][j] + 4*U_n[ns-1][j] + 2*ds)/3
            U_n[0][j] = 0 # S=0 
        
        for i in range(len(i_points)-1):
            U_n[i][nv] = s_points[i]
            P = 1 - r*dt - r*i*dt - ((3*kappa*eta*dt)/(2*dv))
            Q = r*i*dt
            R = kappa*eta*dt / (2*dv)

            U_n[i][0] = P*U[i][0] + Q*U[i+1][0] + R*(4*U[i][1] - U[i][2])
        
        U = copy.deepcopy(U_n)
        U_time.append(copy.deepcopy(U))
        if (x%100 == 0): print("100 steps donee", end="   ")

    U = U[0:I + 1, 0:J + 1]
    return U, U_time
