import numpy as np

def calculate_vacf(velocity_trajectory, max_lag=None):
    n_steps, n_particles, _ = velocity_trajectory.shape
    
    if max_lag is None:
        max_lag = n_steps // 2 
        
    vacf = np.zeros(max_lag)
    counts = np.zeros(max_lag)
    
    stride = 10 
    
    for t0 in range(0, n_steps - max_lag, stride):
        v0 = velocity_trajectory[t0]
        norm = np.sum(v0 * v0) / n_particles 
        
        for t in range(max_lag):
            vt = velocity_trajectory[t0 + t]
            dot_product = np.sum(vt * v0) / n_particles
            vacf[t] += dot_product / norm
            counts[t] += 1
            
    return vacf / counts

def calculate_dos(vacf, dt):
    n_steps = len(vacf)

    window = np.hanning(n_steps)
    vacf_windowed = vacf * window
    
    dos_intensity = np.abs(np.fft.rfft(vacf_windowed))
    frequencies = np.fft.rfftfreq(n_steps, d=dt)
    
    return frequencies, dos_intensity