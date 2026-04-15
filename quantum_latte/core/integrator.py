import numpy as np

class MDEngine:
    def __init__(self, positions, potential, box_size, dt=0.005, mass=1.0):
        self.N = positions.shape[0]
        self.L = box_size
        self.dt = dt
        self.mass = mass
        self.potential = potential
        
        self.positions = positions.copy()
        
        self.velocities = np.random.randn(self.N, 2)
        self.velocities -= np.mean(self.velocities, axis=0) 
        
        self.forces = self.potential.compute_forces(self.positions, self.L)
        
        self.velocity_trajectory = []

    def verlet_step(self):
        self.positions += self.velocities * self.dt + 0.5 * (self.forces / self.mass) * self.dt**2

        self.positions = self.positions % self.L
        
        new_forces = self.potential.compute_forces(self.positions, self.L)
        
        self.velocities += 0.5 * ((self.forces + new_forces) / self.mass) * self.dt
        
        self.forces = new_forces

    def run_simulation(self, steps, eq_steps=1000):
        print(f"Equilibrating for {eq_steps} steps...")
        for _ in range(eq_steps):
            self.verlet_step()
            
        print(f"Starting production run for {steps} steps...")
        for _ in range(steps):
            self.verlet_step()
            self.velocity_trajectory.append(self.velocities.copy())
            
        print("Simulation complete.")
        return np.array(self.velocity_trajectory)