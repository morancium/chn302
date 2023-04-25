import sympy
import math
import json

def newton_raphson(f, x0, tolerance=1e-6, max_iterations=10000):
    """
    Find the root of the function `f` using the Newton-Raphson method.
    
    Arguments:
    - `f`: the function to find the root of.
    - `x0`: the initial guess for the root.
    - `tolerance`: the desired tolerance for the root (default: 1e-6).
    - `max_iterations`: the maximum number of iterations to perform (default: 100).
    
    Returns:
    - The approximate root of the function `f`.
    
    Raises:
    - ValueError if the derivative of `f` is zero at `x0`.
    """
    x = sympy.Symbol('x')
    f_prime = sympy.diff(f(x), x)
    f_prime_func = sympy.lambdify(x, f_prime)
    f_func = sympy.lambdify(x, f(x))
    
    if f_prime_func(x0) == 0:
        raise ValueError("The derivative of the function is zero at the initial guess.")
    
    for i in range(max_iterations):
        x_new = x0 - f_func(x0) / f_prime_func(x0)
        if abs(x_new - x0) < tolerance:
            return x_new
        x0 = x_new
        # print(x0)
    raise ValueError(f"Failed to converge after {max_iterations} iterations.")

def Ki(comp,T,P):
    a=comp[0]
    b=comp[1]
    c=comp[2]
    ki=math.pow(10,a-(b/(T+c)))/P
    # print("Ki",ki,a,b,c,T,P)
    return ki

def Hji(comp,T):
    c1,c2,c3=comp
    T+=273.0
    return (c1*T+c2*T**2+c3*T**3)/1e6

def my_function(x):
    return x**2-1

# try:
#     root = newton_raphson(my_function, 0.002)
#     print(root)
# except ValueError as e:
#     print(str(e)) # "The derivative of the function is zero at the initial guess."
