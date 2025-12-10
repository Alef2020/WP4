import math
from buckling_check import *
from scipy.optimize import minimize

def f(params, p = 1, rho = 23):
    R, t, L = params
    volume = 2*R*math.pi*t*L
    mass = volume * rho
    return mass


def check(params, youngs_modulus =):
    R , t, L = params


conditions = [
    {'type':'ineq', 'fun': x_min}
]

for material in materials:
    youngs_modulus = material.young_moduluss

    result = minimize(
    f, 
    [0,1], 
    constraints=conditions,
    bounds=[(-4, 5),(-2,5)], 
    method='SLSQP',
    options={'disp':True}
    )

print(result.x[0], result.x[1])