# imports
import math
from scipy.integrate import odeint
# Initial conditions
dx0, dy0, vx0, vy0, dx1, dy1, ... = v0 \
= [1.0, -0.51, 1, -0.49, 1.0, 0.5, ...] # initial conditions
# etc. for other neurons up to #5; we break symmetry for n0 to start
# oscillation
def FHND(v, t):
    dx0, dy0, vx0, vy0, dx1, dy1, [...] = v  # variables
    # Parameters
    e0 = e1 = 0.08
    g = 0.8
    b0 = b1 = 0.46
    Dhead = -0.2  # diffusion constants
    Drest = -0.02
    Dgap  = -0.05
    return [
    vx0 - vx0 ** 3 / 3 - vy0 + Dhead * max(dx0 - vx0, 0),
    e0 * (vx0 - g * vy0 + b0),                           
    dx0 - dx0 ** 3 / 3 - dy0 + Dhead * max(-dx0 + vx0, 0),
    e0 * (dx0 - g * dy0 + b0),
    vx1 - vx1 ** 3 / 3 - vy1 + Drest * max(dx1 - vx1, 0),
    e1 * (vx1 - g * vy1 + b1),
    ...
    ] # return time derivatives

sol = odeint(FHND, v0, t)
