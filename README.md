[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]()

# Quantum Harmonic Oscillator (QHO)
In Quantum Mechanics, the harmonic oscillator is an important pradigm because it provides a model for a variety of systems, such as the modes of the __electrodynamic field (photons)__ and __the vibrations of molecules and solids (phonons)__.

The Quantum Harmonic Oscillator (QHO) is a fundamental model in quantum mechanics, describing systems where a particle is subjected to a potential proportional to the square of its displacement from equilibrium. Examples include molecular vibrations and oscillations of atoms in a solid.

## Classical harmonic oscillator
The simple quantum harmonic oscillator is:
<p align= 'center'>
    $$ V(x) = \frac{1}{2} \times m \times \omega \times x^{2}$$
</p>

where,
- _m_ is the mass of the particle
- $\omega$ is the angular frequency of the oscillator
- _x_ is the displacement from equilibrium

## $Schr\ddot{o}dinger$ equation of the quantum harmonic oscillator

In order to solve this problem quantum mechanically, we follow our standard steps. The $schr\ddot{o}dinger$ equation of the harmonic oscillator is given by

<p align='center'>
    $$E \psi (x) = - \frac{\hbar}{2m} \psi^{"}(x) + \frac{1}{2} m \omega^2 x^2 \psi(x)$$
</p>

where,
- $\hbar$ is the reduced Planck's constant
- $\psi$ is the wavefunction
- $E$ is teh energy of the state

This equation is again a linear differential equation of second order, but now one coefficient is position dependent. From our general considerations we already can anticipate the following:

- The solution $\psi(x)$ are continous.
- The deribatives $\psi^{'}(x)$ of the solutions are also continious.
- As $V(x) \rightarrow \infty$ for $|x| \rightarrow \infty$, the particle cannot excape to infinity at finite energy $E$. This only permits bound states, which decay $\psi(x) \rightarrow 0$ as $|x| \rightarrow \infty$.
- The energies of the bound states are discrete (i.e., only at certain energies we can find valid solutions of the schr$\ddot{o}$inger equation).
- Thground state energy $E_{0}$ will be larger than the classical minimal energy: $E_0 > 0$.
