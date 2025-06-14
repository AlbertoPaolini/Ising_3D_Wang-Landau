# Ising_gauge_3D_Wang-Landau
Monte Carlo simulation based on Wang-Landau algorithm for 2D Ising model

# Wang-Landau Algorithm for the 3D Ising gauge Model

This project implements the Wang-Landau algorithm to compute the density of states and derive thermodynamic quantities of the 3D Ising gauge model (e.g. average energy, specific heat). In the Ising gauge model, differently from the usual Ising model, the variables are the links. In this simulation, I used the energy as the sum of the so called "plaquette", i.e. the smallest loop in the lattice. The plaquette is defined as the product of 4 links, that form the sides of a square of side-lenght 1.

## Description

The Wang-Landau algorithm is a relatively recent Monte Carlo method more efficent than the previous ones.
The idea of the Wang-Landau algorithm is to directly estimates the density of states g(E) of a system by performing a random walk in energy space.
Thus, the partition function Z can be written as Z = ∑<sub>E</sub> g(E)exp(-βE). The main advantage il that g(E) does not depend on the temperature, we can construct canonical 
distributions at any temperature if we can estimate g(E) with high accuracy for all energies. 


This simulation performs the following steps:

1. Initialize a 2D lattice of spins with random configuration.
2. Perform the Wang-Landau algorithm to estimate the density of states.
3. Compute average energy <E>, energy variance <E^2> - < E >^2, and specific heat for a range of temperatures from 0.5 to 4.1.
4. Repeat simulations for different lattice sizes and save the data.

## Wang-Langau algorithm

The Wang-Landau algorithm is as follows:

1. Initialize the lattice with an arbitrary configuration with energy E1.
2. In principle we don't know g(E), so the simplest approach is to set g(E)=1 for all the possible energies.
3. nitialize the histogram H(E)=0 for all possible energies, which records how frequently each energy level is visited.
4. Initialize the modificator factor f, we suggest f=exp(1).
5. Pick a spin randomly and change its state, now the energy of the system is E2.
6. The new state it's accepted is p(E1->E2) = min(g(E1)/g(E2) ; 1).
7. If the new state with energy E2 is accepted, we multiply g(E2)->g(E2)xf and H(E2)->H(E2)+1.
8. If the new state with energy E2 is not accepted, we multiply g(E1)->g(E1)xf and H(E1)->H(E1)+1.
9. We proceed with the random walk in energy space until we obtain a "flat" histogram H(E), then we reduce the modificator factor. We suggest f->√f .
10. The steps 5-9 are repeated until f<f<sub>Threshold</sub>.

Because g(E) becomes very large, in practice it is preferable to work with the logarithm of the density of states, so that all possible ln g(E) will fit into double precision numbers.
## Structure

- Num_simulation: Allow to perform more independent simulations.
- 'Initialization_Lattice': In this function we initialize the lattice. The elements of the vectors stored in Lattice are the link variables. The first element of each vector is the link along the first direction, the second element is the link along the second direction and the third the link along the third direction. The NN_Map contain the positions of the near-neighbour of each site, thus the element NN_Map[i][0] is the position, in Lattice, of the near neighbour of the i-th site along the direction 0.
- 'Initialization_Plaquette': The plaquettes are initializated in this function, in this way do not have to calculate them in the Wang-Landau process. Plaquette will contain the value of the plaquettes, meanwhile in the Lattice_Plaq the position of the plaquettes attached to each link. In this way, Lattice[i][j] is a vector that contain the position, in the vector Plaquette, of the plaquettes connected to the link in the direction j of the i-th site.
- `Wang_Landau`: Main function implementing the Wang-Landau sampling.
- `compute_total_energy`: Computes the total energy of the system.
- `mean_energy_from_log_g`: Computes <E> from log(g(E)).
- `mean_energy2_from_log_g`: Computes <E^2> from log(g(E)).
- Visualization of specific heat and average energy vs. temperature.

## Output

Two plots are generated for each set of simulations:

- Specific Heat per Spin vs Temperature
- Average Energy per Spin vs Temperature

Additionally, two text files are created for post-processing:

- `specific_heat_data.txt`: Contains specific heat data.
- `Variance_energy_data.txt`: Contains energy variance data.
