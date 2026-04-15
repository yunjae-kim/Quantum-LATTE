import numpy as np

class LennardJones:
    def __init__(self, epsilon=1.0, sigma=1.0):
        self.epsilon = epsilon
        self.sigma = sigma

    def compute_forces(self, positions, box_size):
        delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
        
        delta = delta - box_size * np.round(delta / box_size)
        
        r_sq = np.sum(delta**2, axis=2)
        
        np.fill_diagonal(r_sq, np.inf)

        term6 = (self.sigma**2 / r_sq)**3
        term12 = term6**2
        
        f_mag_over_r = 24 * self.epsilon * (2 * term12 - term6) / r_sq

        forces = np.sum(f_mag_over_r[:, :, np.newaxis] * delta, axis=1)
        
        return forces