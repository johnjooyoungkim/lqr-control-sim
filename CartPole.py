
# Cart-Pole class

import numpy as np

class CartPole:

    def __init__(self, M = 5, m = 0.2, l = 0.5, g = 9.8, b = 0.0, I = None, sensor_noise_std = 0):
        self.M = M # Mass of Cart
        self.m = m # mass of pendulum bob
        self.l = l # length from pivot to pendulum's center of mass
        self.g = g # graviational acceleration
        self.b = b # cart damping (friction between cart and track)
        self.I = I
        self.sensor_noise_std = sensor_noise_std   

    def cartpole_dynamics(self, z, u):
        M, m, l, g, b, I = self.M, self.m, self.l, self.g, self.b, self.I
        if I is None:
            I = (1/3) * m * l**2

        x, x_dot, theta, theta_dot = z
        sin_t, cos_t = np.sin(theta), np.cos(theta)

        if isinstance(u, np.ndarray):
            u = u.item()

        F = u - b*x_dot + m*l*theta_dot**2*sin_t
        D = (M+m)*(I + m*l**2) - (m*l*cos_t)**2

        x_ddot = ((I + m*l**2)*F + m**2*g*l**2*sin_t*cos_t) / D
        theta_ddot = (-m*g*l*(M+m)*sin_t - m*l*cos_t*F) / D
        
        return np.array([x_dot, x_ddot, theta_dot, theta_ddot])

    def sensor_data(self, z): # returns [ x, theta ] with sensor noise and 
        x, x_dot, theta, theta_dot = z

        x_measured = x + np.random.normal(0, self.sensor_noise_std)
        theta_measured = theta + np.random.normal(0, self.sensor_noise_std)

        return np.array([x_measured, theta_measured])