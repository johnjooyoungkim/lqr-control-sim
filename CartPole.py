
# Cart-Pole class

import numpy as np

class CartPole:

    def __init__(self, 
                M = 5, # Mass of Cart
                m = 0.2, # mass of pendulum bob
                l = 0.5, # length from pivot to pendulum's center of mass
                g = 9.8, # graviational acceleration
                b = 0.0, # cart damping (friction between cart and track)
                I = None
                ):
        self.M = M
        self.m = m
        self.l = l
        self.g = g
        self.b = b
        self.I = I

    def cartpole_dynamics(self, z, u):
        M = self.M
        m = self.m
        l = self.l
        g = self.g
        b = self.b
        I = self.I
        
        if I is None:
            I = (1/3) * m * l**2  # uniform rod about center of mass
        
        x, x_dot, theta, theta_dot = z
        sin_t, cos_t = np.sin(theta), np.cos(theta)
        
        u = u.item()

        F = u - b*x_dot + m*l*theta_dot**2*sin_t
        D = (M+m)*(I + m*l**2) - (m*l*cos_t)**2
        
        x_ddot = ((I + m*l**2)*F - m**2*g*l**2*sin_t*cos_t) / D
        theta_ddot = (m*g*l*(M+m)*sin_t - m*l*cos_t*F) / D
        
            
        return np.array([x_dot, x_ddot, theta_dot, theta_ddot])
