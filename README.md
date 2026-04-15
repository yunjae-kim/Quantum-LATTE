# Quantum LATTE (Light Approximation of Thermal and Trajectory Evolution)

Quantum LATTE is a 2D molecular dynamics simulation of argon atoms interacting via the Lennard-Jones potential. It computes the velocity autocorrelation function (VACF) and phonon density of states (PDOS) using classical Newtonian mechanics.

---

## Theory

### The Lennard-Jones Potential

Real atoms neither pass through each other nor stick together permanently. The Lennard-Jones (LJ) potential captures both of these behaviors with two terms:

$$V(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$$

where $r$ is the distance between two atoms, $\varepsilon$ is the depth of the potential well (how strongly they attract), and $\sigma$ is the distance at which $V = 0$ (effective atom diameter). The $r^{-12}$ term represents the repulsion between atoms, which strongly resist being pushed through each other (Pauli exclusion). The $r^{-6}$ term represents attraction from weak van der Waals forces that pull atoms together at moderate distances. At $r = 2^{1/6}\sigma$, the potential is minimized: this is the equilibrium bond length.

This simulation uses reduced units where $\varepsilon = 1$, $\sigma = 1$, and $m = 1$, so all quantities are dimensionless ratios relative to the LJ parameters of argon.

### Equations of Motion

Given the potential $V(r)$, the force on atom $i$ due to atom $j$ is:

$$\vec{F}_{ij} = -\nabla_i V(r_{ij})$$

which works out to:

$$\vec{F}_{ij} = \frac{24\varepsilon}{r^2} \left[ 2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right] \vec{r}_{ij}$$

Each atom obeys Newton's second law, $\vec{F} = m\vec{a}$. The simulation uses these equations forward in time for all 64 atoms.simultaneously.

### Velocity Verlet Integration

To move atoms through time, we use the Velocity Verlet algorithm, a numerical integrator that conserves energy much better than simple Euler methods:

$$\vec{r}(t + \Delta t) = \vec{r}(t) + \vec{v}(t)\Delta t + \frac{1}{2}\frac{\vec{F}(t)}{m}\Delta t^2$$

$$\vec{v}(t + \Delta t) = \vec{v}(t) + \frac{1}{2}\frac{\vec{F}(t) + \vec{F}(t+\Delta t)}{m}\Delta t$$

This is a symplectic integrator: it preserves the geometric structure of Hamiltonian mechanics, preventing the total energy from drifting over long simulations.

### Periodic Boundary Conditions

The simulation box (size $L = 8.0$) uses periodic boundary conditions (PBC): atoms that exit one side re-enter from the opposite side. This mimics bulk matter by eliminating surface effects. The minimum image convention ensures each atom only interacts with the nearest image of every other atom:

$$\Delta x_{\text{corrected}} = \Delta x - L \cdot \text{round}\!\left(\frac{\Delta x}{L}\right)$$

### Velocity Autocorrelation Function (VACF)

The VACF measures how correlated an atom's velocity is with its own velocity at a later time:

$$C(t) = \frac{\langle \vec{v}(t_0) \cdot \vec{v}(t_0 + t) \rangle}{\langle \vec{v}(t_0) \cdot \vec{v}(t_0) \rangle}$$

- $C(0) = 1$ always (perfectly correlated with itself)
- In a gas, $C(t)$ decays quickly to zero: atoms forget their direction after a collision
- In a solid, $C(t)$ oscillates: atoms are trapped in a potential well and bounce back

This simulation runs at $T^* = 0.5$ and $\rho^* \approx 1.0$, placing it firmly in the solid regime of the LJ phase diagram.

### Phonon Density of States (PDOS)

The PDOS tells us how many vibrational modes exist at each frequency. It is obtained by taking the Fourier transform of the VACF (Green-Kubo relation):

$$g(\omega) = \int_0^\infty C(t) \cos(\omega t) \, dt$$

In practice, we use a discrete FFT with a Hanning window to reduce spectral leakage:

$$g(\omega) = \left| \mathcal{F}\left[ C(t) \cdot w(t) \right] \right|$$

where $w(t) = \frac{1}{2}\left(1 - \cos\!\left(\frac{2\pi t}{T}\right)\right)$ is the Hanning window.

Peaks in $g(\omega)$ correspond to characteristic vibrational frequencies of the solid — analogous to phonon branches in a crystal.

---

## Time Complexity

| Component | Complexity | Notes |
|-----------|-----------|-------|
| Force calculation (per step) | $O(N^2)$ | Vectorized over all pairs |
| Verlet integration (per step) | $O(N)$ | Trivially parallel |
| Total simulation | $O(S \cdot N^2)$ | $S$ = number of steps |
| VACF calculation | $O(T_0 \cdot \tau)$ | $T_0$ = origins, $\tau$ = max lag |
| FFT (PDOS) | $O(\tau \log \tau)$ | Fast Fourier Transform |

---

## Results

### Velocity Autocorrelation Function

The VACF shows rapid initial decay followed by weak oscillations, consistent with a dense LJ solid:
- Fast decorrelation at short times reflects strong collisional forces
- Oscillatory tail reflects atoms bouncing within their local potential well
- Does not decay to zero, characteristic of solid-like caging

### Phonon Density of States

The raw PDOS shows many sharp peaks (discrete vibrational modes from the finite system). After Gaussian smoothing ($\sigma = 2$), the spectrum reveals broad phonon bands with:
- A peak near $\omega \approx 2$ (transverse acoustic modes)
- A dominant peak near $\omega \approx 3$ (longitudinal/optical modes)  
- A secondary feature near $\omega \approx 5.5$
- Hard cutoff around $\omega \approx 9$–$10$ (Debye-like frequency cutoff)

The absence of a zero-frequency diffusion peak confirms the system is in a **solid phase**, consistent with the LJ phase diagram at $T^* = 0.5$, $\rho^* \approx 1.0$.

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `NUM_PARTICLES` | 64 | Number of LJ atoms |
| `BOX_SIZE` | 8.0 | Simulation box side length (reduced units) |
| `DT` | 0.005 | Time step |
| `T_target` | 0.5 | Target reduced temperature $T^*$ |
| `EQUILIBRATION_STEPS` | 2000 | Steps before recording begins |
| `PRODUCTION_STEPS` | 8000 | Steps used for analysis |
| `max_lag` | 2000 | Maximum lag for VACF |
