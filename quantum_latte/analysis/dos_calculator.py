import numpy as np

def calculate_dos(vacf, dt):
    """
    Computes the Density of States (DOS) from the VACF.
    
    vacf: 1D numpy array from autocorrelation.py
    dt: Time step used in the simulation
    Returns: (frequencies, intensities)
    """
    n_steps = len(vacf)
    
    window = np.hanning(n_steps)
    vacf_windowed = vacf * window
    
    dos_intensity = np.abs(np.fft.rfft(vacf_windowed))
    
    frequencies = np.fft.rfftfreq(n_steps, d=dt)
    
    return frequencies, dos_intensity

def find_peaks(frequencies, dos_intensity):
    """Simple utility to find the dominant vibrational frequency."""
    max_idx = np.argmax(dos_intensity)
    return frequencies[max_idx]