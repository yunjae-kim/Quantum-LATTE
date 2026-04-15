import numpy as np
from scipy.ndimage import gaussian_filter1d

def calculate_debye_dos(frequencies, speed_of_sound, density):
    """
    Calculates the theoretical Debye DOS: g(w) = (V / 2*pi^2 * v^3) * w^2
    For simplicity in reduced units, we use g(w) = A * w^2
    """
    # In the Debye model, DOS is proportional to frequency squared
    # until it hits the Debye Cutoff frequency.
    debye_dos = frequencies**2
    
    # Normalize so the heights are comparable for visualization
    if np.max(debye_dos) > 0:
        debye_dos = debye_dos / np.max(debye_dos)
        
    return debye_dos

def smooth_dos(dos_intensity, sigma=2):
    """
    Applies a Gaussian filter to smooth out the noise in the raw MD DOS.
    sigma: Standard deviation for Gaussian kernel (higher = smoother)
    """
    return gaussian_filter1d(dos_intensity, sigma=sigma)

# Example Visualization Script
def plot_comparison(freqs, raw_dos, dt):
    # 1. Smooth the raw data
    smoothed = smooth_dos(raw_dos, sigma=3)
    
    # 2. Get Debye Theory (scaled to match the max of smoothed data)
    theory = calculate_debye_dos(freqs, 1.0, 1.0)
    theory *= np.max(smoothed) 
    
    plt.figure(figsize=(8, 5))
    plt.plot(freqs, raw_dos, alpha=0.3, label='Raw MD DOS (Noisy)', color='gray')
    plt.plot(freqs, smoothed, label='Smoothed MD DOS', linewidth=2, color='blue')
    plt.plot(freqs, theory, '--', label='Debye Theory ($omega^2$)', color='red')
    
    plt.title("Quantum LATTE: MD vs. Debye Model")
    plt.xlabel("Frequency ($omega$)")
    plt.ylabel("Density of States $g(omega)$")
    plt.xlim(0, freqs.max() / 4) # Zoom into the physical frequency range
    plt.ylim(0, np.max(smoothed) * 1.2)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()