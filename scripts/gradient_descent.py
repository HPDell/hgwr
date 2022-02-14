"""
Gradient descent algorithm for linear regression.
"""
from typing import Callable, Tuple
import numpy as np


def dri_cost_lm(beta: np.ndarray, args: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    x, y = args
    return (np.sum((x @ beta - y) * x, axis=0)).reshape(x.shape[1], 1) / x.shape[0]


def cost_lm(beta: np.ndarray, args: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    x, y = args
    return np.sum((np.matmul(x, beta) - y) ** 2) / x.shape[0] / 2


def dri_cost_gwr(beta: np.ndarray, args: Tuple[np.ndarray, np.ndarray, np.ndarray, float]) -> np.ndarray:
    x, y, u, b = args
    b2 = b * b
    n_data, _ = x.shape
    d_beta = np.zeros(x.shape)
    for i in range(n_data):
        d = np.sum((np.apply_along_axis(lambda u1, u0: u1 - u0, 1, u, u[i,:])) ** 2, axis=1)
        w = np.exp(- d / b2 / 2)
        wxb = w * (np.sum(x * beta[i,:], axis=1) - y)
        f = np.apply_along_axis(lambda xj, w: xj * w, 0, x, wxb)
        d_beta[i,:] = np.sum(f, axis=0)
    return d_beta  / n_data / n_data


def cost_gwr(beta: np.ndarray, args: Tuple[np.ndarray, np.ndarray, np.ndarray, float]) -> float:
    x, y, u, b = args
    b2 = b * b
    n_data, _ = x.shape
    f = 0.0
    for i in range(n_data):
        d = np.sum((np.apply_along_axis(lambda u1, u0: u1 - u0, 1, u, u[i,:])) ** 2, axis=1)
        w = np.exp(- d / b2 / 2)
        f += np.sum(((w * (np.sum(x * beta[i,:], 1) - y)) ** 2))
    return f / n_data / n_data / 2


GradientFunctionType = Callable[[np.ndarray, Tuple], np.ndarray]


def gradient_descent(params: np.ndarray, fun: GradientFunctionType, dri: GradientFunctionType, args: Tuple, 
                     alpha: float=1e-2, eps: float=1e-10, max_iters: int=10 ** 6):
    iter = 0
    diff = np.inf
    beta0 = params
    beta1 = beta0
    j0 = fun(beta0, args)
    while diff > eps and iter <= max_iters:
        beta1 = beta0 - dri(beta0, args) * alpha
        j1 = fun(beta1, args)
        diff = abs(j1 - j0)
        beta0 = beta1
        j0 = j1
        iter += 1
    return beta1
        

def gradient_lm():
    x = np.loadtxt("../data/x.csv", dtype=np.double)
    y = np.loadtxt("../data/y.csv", dtype=np.double).reshape(x.shape[0], 1)
    beta = np.ones((3,1))
    beta_hat = gradient_descent(params=beta, fun=cost_lm, dri=dri_cost_lm, args=(x, y))
    print(beta_hat)


def gradient_gwr():
    x = np.loadtxt("../data/gwr_x.csv", dtype=np.double)
    y = np.loadtxt("../data/gwr_y.csv", dtype=np.double)
    u = np.loadtxt("../data/gwr_u.csv", dtype=np.double)
    b = 270.0
    beta = np.ones(x.shape)
    beta_hat = gradient_descent(params=beta, fun=cost_gwr, dri=dri_cost_gwr, args=(x, y, u, b), alpha=0.1)
    print(beta_hat)


if __name__=="__main__":
    gradient_gwr()