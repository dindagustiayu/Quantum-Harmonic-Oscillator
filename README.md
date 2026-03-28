[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]()

# Quantum Harmonic Oscillator (QHO)
In Quantum Mechanics, the harmonic oscillator is an important pradigm because it provides a model for a variety of systems, such as the modes of the __electrodynamic field (photons)__ and __the vibrations of molecules and solids (phonons)__.

The Quantum Harmonic Oscillator (QHO) is a fundamental model in quantum mechanics, describing systems where a particle is subjected to a potential proportional to the square of its displacement from equilibrium. Examples include molecular vibrations and oscillations of atoms in a solid.

# Prior knowledge
- Classical Mechanics
- Quantum Fundamentals
- Mathematical Integration

# Preliminaries
- `scipy.integrate`: definite integrals.
- `dblquad`: double integrals.
- `tplquad`: triple integrals.

## Classical harmonic oscillator
The classical harmonic oscillator describes a particle subject to a restoring force $F=-m \omega^2 x$ proportional to the distance from an equilibrium position $x=0$. Newton's equation $m\ddot{x} = F$ results in an oscillatory motion $x(t)=x_0 \ cos \omega t + (v_0/ \omega) \ sin \omega t$, where $\omega = \frac{2 \pi}{T}$ and $T$ is the oscillation period. In this solution, $x_0 = x(0)$ is the initial position and $v_0 = \dot{x} (0)$ is the initial velocity of the particle. According to $F = -V^{'}$, the force $F = -m \omega^2 x$ corresponds to a parabolic potential energy.

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
- The deribatives $\psi'(x)$ of the solutions are also continious.
- As $V(x) \rightarrow \infty$ for $|x| \rightarrow \infty$, the particle cannot excape to infinity at finite energy $E$. This only permits bound states, which decay $\psi(x) \rightarrow 0$ as $|x| \rightarrow \infty$.
- The energies of the bound states are discrete (i.e., only at certain energies we can find valid solutions of the $schr\ddot{o}inger$ equation).
- The ground state energy $E_{0}$ will be larger than the classical minimal energy: $E_0 > 0$.

## Ground state and first excited state

The first two solutions correspond to the ground state and the first excited state of the harmonic oscillator:
 
__Ground state__: $n=0, \ f=1$, which gives to $f' = f^{"} = 0$. This is a solution when $\varepsilon \neq 2n + 1$ if $\epsilon =1 $, which corresponds to an energy $E_{0} = \frac{1}{2} \hbar \omega$. The wavefunction is given by:

<p align='center'>
    $$\psi_{0} (x)= c_{0} exp \left (-x^{2} \frac{m \omega}{2 \hbar} \right)$$
</p>

where $c_{0}$ can be determined from the normalisation condition,

<p align='center'>
    $$\int_{- \infty}^{\infty} |\psi (x)^{2} \ dx = 1 $$
</p>

This gives $c_{0} = \left (\frac{m \omega}{\pi \hbar} \right)^{1/4}$.

## Normalisation of the ground state wavefunction

The ground state wave function of the harmonic oscillator provides us with a good occasion to practice once more the normalisation of the wavefunction. Since we want to interpret $| \psi_{0} (x)^{2} | = P (x)$ as the probability density for position, we require,

<p align='center'>
    $$\int_{- \infty}^{\infty} |\psi_{0} (x)^{2} | \ dx = 1$$
</p>

With $\psi_{0} (x) = c_{0} exp \left (-x^{2} \frac{m \omega}{2 \hbar} \right)$ this integral reads,

<p align='center'>
    $$ \int_{- \infty}^{\infty} c_{0}^{2} exp \left (-x^{2} \frac{m \omega}{\hbar} \right) \ dx$$
</p>

We use the standard integral $\int_{- \infty}^{\infty} exp(-ax^{2}) \ dx = \sqrt{\pi / a}$, where we set $a = \frac{m \omega}{\hbar}$:

<p align='center'>
    $$ \int_{- \infty}^{\infty} c_0^{2} exp \left(-x^{2} \frac{m \omega}{\hbar} \right) \ dx = c_0^{2} \sqrt{\frac{\pi \hbar}{m \omega}} = 1$$
</p>

The equations to 1 if $c_0 = \left(\frac{m \omega}{\pi \hbar} \right)^{1/4}$, in agreement with the value given in the previous section. 

This Python script illustrates the one-dimensional harmonic oscillator ground state wavefunction with the following parameters:

```Python
import numpy as np
from scipy.integrate import quad, dblquad, tplquad

def func(x):
    return np.exp(-x**2)

I, err = quad(func, -np.inf, np.inf)
c0 = 1 / np.sqrt(I)
print('I:', I)
print('c0:', c0)
```
```
I: 1.7724538509055159
c0: 0.7511255444649425
```

The probability density of the oscillator's position is given by $P_{0} (x) = | \psi_{0} (x)|^{2}$ and is non zero outside the classical turning points, $\neq \alpha^{- 1/2}$, a phenomenon known as tunneling. We will calculate the probability f tunneling for an oscillator in the state $\psi_{0}$. 

The wavefunction is symmetric about $x=0$, so the probability of tunneling is

<p align='center'>
    $$ \begin{align} P(x < -\alpha) + P(x > \alpha) &= 2P(x > \alpha) = 2 \sqrt{\frac{\alpha}{\pi}} \int_{\alpha^{-1/2}}^{\infty} exp(-\alpha x^{2}) \ dx \\ &= \frac{2}{\sqrt{\pi}} \int_{1}^{\infty} e^{-y^{2}} \ dy \end{align}$$
</p>

```Python
I, err = quad(func, 1, np.inf)
ptun = 2 / np.sqrt(np.pi) * I
print('I:', I)
print('Ptunneling:', ptun)
```
```
I: 0.13940279264033098
Ptunneling: 0.1572992070502851
```

