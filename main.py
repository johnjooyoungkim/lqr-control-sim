
import matplotlib.pyplot as plt
import control
import numpy as np

from CartPole import CartPole
from controller import Controller

np.random.seed(0) #

# Parameters defining the sys
n = 4000 # Number of samples
dt = 0.01 # System timebase

M = 5 # Mass of Cart
m = 1 # mass of pendulum bob
l = 2 # length from pivot to pendulum's center of mass
g = 9.8 # graviational acceleration
b = 1.0 # cart damping (friction between cart and track)
I = None # moment of Inertia of the rod

sensor_noise_std = 0.01

z0 = [0, 0, 3.0, 0.1] # initial state [ x , x_dot , theta , theta_dot ]
z_eq = np.array([3, 0, np.pi, 0]) # equilibrium state (upright)

cart_pole = CartPole(M, m, l, g, b, I, sensor_noise_std)


# building A B C matrices
if I is None:
    I = (1/3) * m * l**2

D = (M+m)*(I + m*l**2) - (m*l)**2

A = [[0, 1, 0, 0],
     [0, -(I + m*l**2)*b/D, m**2*g*l**2/D, 0],
     [0, 0, 0, 1],
     [0, -m*l*b/D, m*g*l*(M+m)/D, 0]]

B = [[0],
     [(I+m*l**2)/D],
     [0],
     [m*l/D]] 

C = [[1, 0, 0, 0],
     [0, 0, 1, 0]]

# LQR: building Q and R

Q = [[1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 3, 0],
     [0, 0, 0, 2]]

R = [[10]]


def rk4_step(z, u, dt, cartpole):
    k1 = cartpole.cartpole_dynamics(z, u)
    k2 = cartpole.cartpole_dynamics(z + dt/2 * k1, u)
    k3 = cartpole.cartpole_dynamics(z + dt/2 * k2, u)
    k4 = cartpole.cartpole_dynamics(z + dt * k3, u)
    return z + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)


# Closed Loop Simulation
z = z0.copy()
y = cart_pole.sensor_data(z)
state_history = [z.copy()]
sensor_history = [y.copy()]

# Pole placement
#p = [-.1,-1.3, -1.4, -0.5]
#K = control.place(A, B, p)

# LQR
K, S, E = control.lqr(A, B, Q, R)

u_history = []

# Initialize Kalman Filter (Q, R, P_0) 
Q_k = np.diag([1e-4, 1e-3, 1e-4, 1e-3])

R_k = np.diag([sensor_noise_std**2, sensor_noise_std**2])

P0 = np.diag([0, 0, 0, 0])

controller = Controller(A, B, C, Q_k, R_k, K, z_eq, z0, P0, dt)


# step
u = [0]

for i in range(n-1):
    y = cart_pole.sensor_data(z)
    u = controller.step(y, u)
    u_history.append(float(u.item()))
    z = rk4_step(z, u, dt, cart_pole)
    sensor_history.append(y.copy())
    state_history.append(z.copy())

state_history = np.array(state_history)
sensor_history = np.array(sensor_history)

t = np.linspace(0, (n-1)*dt, n)

fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

labels = ['Cart position $x$ [m]', 
          'Cart velocity $\\dot{x}$ [m/s]', 
          'Pendulum angle $\\theta$ [rad]', 
          'Angular velocity $\\dot{\\theta}$ [rad/s]']

for i in range(4):
    axs[i].plot(t, state_history[:, i])
    axs[i].set_ylabel(labels[i])
    axs[i].grid(True)
    if i == 2:  # angle plot
        axs[i].axhline(np.pi, color='r', linestyle='--', alpha=0.7, label='upright')
        axs[i].legend()

axs[-1].set_xlabel('Time [s]')
plt.tight_layout()

plt.savefig('states.png', dpi=150)
