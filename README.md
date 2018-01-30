# Biomedical
The 3 order random program completes a simulation of mass transfer between blood bag and medicine supplies.

The solution to differential equations is first order Euler method.

The inner blood bag uses hybrid calculation.

The main data uses double buffers for data compressing and storing.

GA optimization algorithm is used and an optimization of parameters Qb Qf Vt1 Vt2 is done with parallel calculation.

Qb Qf are related to time t.

The volumn element which enters tunnel is chosen randomly.

The sum of tiny volumn elements is larger in blood bag. Each time choose multiple volumn elements randomly and let them into the tunnel, in order to keep the consentration the same inside and outside the container.
