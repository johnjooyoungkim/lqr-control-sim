
# Cart-Pole class

import numpy as np

class CartPole:

    def __init__(self, 
                M = 5, # Mass of Cart
                m = 0.2, # mass of pendulum bob
                l = 0.5, # length from pivot to pendulum's center of mass
                g = 9.8, # graviational acceleration
                b = 0.0 # cart damping (friction between cart and track)
                ):
        self.M = M
        self.m = m
        self.l = l
        self.g = g
        self.b = b

    def cartpole_dynamics(self, z, u):
        M = self.M
        m = self.m
        l = self.l
        g = self.g
        b = self.b
        
        x, x_dot, theta, theta_dot = z
        sin_t, cos_t = np.sin(theta), np.cos(theta)

        x_ddot = (u - b*x_dot + m *l*theta_dot*sin_t - m*g*sin_t*cos_t)/ \
                 (M + m - m*cos_t**2)

        theta_ddot = (g*sin_t - cos_t*(u - b*x_dot + m*l*theta_dot**2*sin_t)/(M+m)) / \
                     (l*(4/3 - m*cos_t**2/(M+m)))
        
        return np.array([x_dot, x_ddot, theta_dot, theta_ddot])
