# Schrödinger's Billiard

Hello! My name is Rafael and in this repository I perform numerical simulations to calculate properties of electrical transport in quantum billiards such as conductance, power shot noise, concurrence, the entanglement of formation and violations of Bell's inequality. The quantum billiards, also known as quantum dots, are extremely small devices about a few nanometers in which the electrons are carried by leads and scattered in a chaotic cavity, as we can see in Figure 1 [[1]](#1).

<p align="center" width="100%">
  Figure 1<br>
    <img width="33%" src="https://www.lorentz.leidenuniv.nl/beenakkr/mesoscopics/topics/chaos/images/billiard.gif">
</p>

## Conductance

The property which quantifies the number of electrons that can pass through the quantum billiard is called quantum conductance. This property is directly related to the probability of an electron going from a lead to the other due to the wave-particle duality of the electron. Moreover, the conductance fluctuates around an average and it has a specific variance, this fact is known as universal fluctuations of conductance. The mean, variance and skewness of the distribution of conductance depend on the number of open channels in each lead, the potential barrier between the lead and the chaotic resonance cavity and especially the class of symmetry of the hamiltonian. In Figure 2, I plot the conductance of the first 200 realizations for each number of open channels (N) in the leads generated by monte carlo simulation. In each step, we generate a random hamiltonian to calculate the scattering matrix "S" by Mahaux-Weidenmuller formula and the conductance "G" by Landauer-Buttiker formula.

<p align="center" width="100%">
  Figure 2<br>
    <img width="75%" src="https://user-images.githubusercontent.com/53999015/146303088-65fac1ac-5f96-42c4-9b59-6736cf3f6de1.png">
</p>

As we can see in this figure, the distribution of realizations of conductance depends on the number of open channels but also on the β which indicates a specific symmetry of the hamiltonian. These symmetries are fundamentals to describe the physics of quantum billiards. In a Schrödinger billiard, we have three fundamental symmetries which are listed below:

#### Gaussian Orthogonal Ensemble (GOE)
The index β = 1 (blue), also indicated in the literature by "AI", represents a system that has time-reversal symmetry (TRS) and spin rotational symmetry (SRS). The Hamiltonian of this system has to be a real-symmetric matrix.
#### Gaussian Unitary Ensemble (GUE)
The β = 2 (red), also indicated by literature as "A", represents a system that the TRS is broken. This happens when there is a magnetic field applied in the quantum dot, the magnetic field breaks the time-reversal symmetry. When the TRS is broken, it doesn't matter if SRS is preserved or broken. In this case, each element of the Hamiltonian is a complex number.
#### Gaussian Symplectic Ensemble (GSE)
The β = 4 (green), indicated by literature as "AII", represents a system that SRS is broken but TRS is preserved. In this case, the Hamiltonian has a quaternionic structure. In these systems, there is a special interaction called spin-orbit coupling that arises due to the relation between the apparent magnetic field seen from the electron perspective and the magnetic moment of the electron associated with its intrinsic spin. 

In this repository, I analyze the influence of these three fundamental symmetries in the conductance of electrons inside a quantum billiard. I obtain the distributions of conductance, as well as the average and variance in terms of the number of open channels for each class of symmetry. I also perform simulations to see how the conductance of each class of symmetry depends on the barrier between the lead and chaotic cavity, which Γ is called the transparency of the barrier. Moreover, we analyze the weak localization correction of conductance as a function of the barrier's transparency Γ. 

The random-matrix theory of quantum transport has a lot of details and here I just give a brief overview of this project. However, reference [[2]](#2) is the standard guide for this topic. 

## Power Shot Noise
One of the founding fathers of mesoscopic physics Rolf Landauer once said: "The noise is the signal". Since then physicists have turned their attention to noise measurements. In mesoscopic physics, the shot noise is the fluctuations in time of the current's measurements due to the discreteness of the electrical charge. In Figure 3 we see a setup that measures the current and its respective noise when a constant bias voltage is applied to a conductor. The current noise is evaluated by using a spectral analyzer [[3]](#3).

<p align="center" width="100%">
  Figure 3<br>
    <img width="75%" src="https://user-images.githubusercontent.com/53999015/146485017-2cc05951-1ad3-4fec-8768-66c1a83a6082.png">
</p>

## References
<a id="1">[1]</a> 
M. J. Berry, J. A. Katine, R. M. Westervelt, and A. C. Gossard, Phys. Rev. B 50, 17721 (1994).

<a id="2">[2]</a> 
C. W. J. Beenakker, Rev. Mod. Phys. 69, 731 – Published 1 July 1997

<a id="3">[3]</a> 
Shot Noise in Mesoscopic Systems: From Single Particles to Quantum Liquids
Kensuke Kobayashi, Masayuki Hashisaka
Journal of the Physical Society of Japan, 90, 102001 (2021) 10.7566/JPSJ.90.102001
