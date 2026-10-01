import numpy as np
import scipy.signal as sig

class Controller:

    def __init__(self, A, B, C, Q_kf, R_kf, K_lqr, z_eq, x_hat0, P0, dt):
        A, B = np.array(A), np.array(B)
        Ad, Bd, _, _, _ = sig.cont2discrete((A, B, np.eye(A.shape[0]), np.zeros((A.shape[0], B.shape[1]))), dt, method='zoh')
        self.A = Ad
        self.B = Bd
        self.C = np.array(C)
        self.Q_d = np.array(Q_kf) * dt
        self.R_kf = R_kf
        self.K_lqr = K_lqr
        self.z_eq = z_eq
        self.x_hat = x_hat0
        self.P = P0

    def step(self, y_measured, u_prev):
        # 1. Kalman filter: predict + update using y_measured
        self.x_hat, self.P = self._kalman_update(self.x_hat, self.P, y_measured, u_prev)
        
        # 2. LQR: compute control from the estimated state
        u = -self.K_lqr @ (self.x_hat - self.z_eq)
        return u

    def _kalman_update(self, x_hat, P, y, u):
        A, B, C, Q, R = self.A, self.B, self.C, self.Q_d, self.R_kf
        n = A.shape[0]

        # Prediction
        delta_hat = x_hat - self.z_eq
        delta_priori = A @ delta_hat + B @ u
        x_priori = delta_priori + self.z_eq
        P_priori = A @ P @ A.T + Q

        # Update 
        Kf = (P_priori @ C.T) @ np.linalg.inv(C @ P_priori @ C.T + R)
        x_posteriori = x_priori + Kf @ (y - C @ x_priori)
        P_posteriori = (np.identity(n) - Kf @ C) @ P_priori
        return x_posteriori, P_posteriori
