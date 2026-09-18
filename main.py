import scipy.signal as sig
import matplotlib.pyplot as plt
import control
import numpy as np

from CartPole import CartPole

# Parameters defining the sys
n = 2000 # Number of samples
dt = 0.01 # System timebase

M = 1 # Mass of Cart
m = 0.1 # mass of pendulum bob
l = 1 # length from pivot to pendulum's center of mass
g = 9.8 # graviational acceleration
b = 10.0 # cart damping (friction between cart and track)
I = None # moment of Inertia of the rod

z0 = [0, 0, np.pi, 0.1] # initial state [ x , x_dot , theta , theta_dot ]

cart_pole = CartPole(M, m, l, g, b, I)


# building the A and B matrices
if I is None:
    I = (1/3) * m * l**2

D = (M+m)*(I + m*l**2) - (m*l)**2

A = [[0, 1, 0, 0],
     [0, -(I + m*l**2)*b/D, m**2*g*l**2/D, 0],
     [0, 0, 0, 1],
     [0, m*l*b/D, m*g*l*(M+m)/D, 0]]

B = [[0],
     [(I+m*l**2)/D],
     [0],
     [m*l/D]] 

# building Q and R

Q = [[0, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 1000, 0],
     [0, 0, 0, 1000]]

R = [[0.1]]

# LQR
K, S, E = control.lqr(A, B, Q, R)

def rk4_step(z, u, dt, cartpole):
    k1 = cartpole.cartpole_dynamics(z, u)
    k2 = cartpole.cartpole_dynamics(z + dt/2 * k1, u)
    k3 = cartpole.cartpole_dynamics(z + dt/2 * k2, u)
    k4 = cartpole.cartpole_dynamics(z + dt * k3, u)
    return z + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)


# Closed Loop Simulation
z = z0.copy()
state_history = [z.copy()]

for i in range(n-1):
    u = -K @ z #controller(z)  # e.g. -K @ (z - z_eq), or 0 for open-loop test
    z = rk4_step(z, u, dt, cart_pole)
    state_history.append(z.copy())

state_history = np.array(state_history)

plt.plot(np.linspace(0, n*dt, n), state_history[:,2])
plt.grid()
plt.savefig('output.png')

