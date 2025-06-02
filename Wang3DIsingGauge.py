import numpy as np
import matplotlib.pyplot as plt
import random
import math
import os
Num_Simulations = 1
Lenght = 10
Height = 10


class Links:
    def __init__(self, link, Pos_Plaq):
        self.link = link
        self.Pos_Plaq = Pos_Plaq

def Position( x, y, z ):
    
    return ( ( (x + Lenght)%Lenght ) * Lenght ) * Height + ( ( y + Lenght)% Lenght ) * Height + ( z + Height ) % Height

def get_coordinates(index):
    z = index % Height
    y = (index // Height) % Lenght
    x = (index // (Lenght * Height)) % Lenght
    return [x, y, z]

#--------------------------------------------------INITIALIZATION LATTICE----------------------------------------------------------------

def Initialization_Lattice():

    Lattice_Links = []
    Lattice_Plaq = []
    NN_Map = []


    for i in range(Lenght):
        for j in range(Lenght):
            for k in range(Height):
                link1 = random.choice([1, -1])
                link2 = random.choice([1, -1])
                link3 = random.choice([1, -1])
                final_links = [link1, link2, link3]
                Lattice_Links.append(final_links)
                Lattice_Plaq.append([[],[],[]])

                nn_position = []
                nn_position.append( Position( i+1, j, k ) )
                nn_position.append( Position( i, j+1, k ) )
                nn_position.append( Position( i, j, k+1 ) )
                NN_Map.append(nn_position)
    
    return Lattice_Links, Lattice_Plaq, NN_Map


#--------------------------------------------------INITIALIZATION PLAQUETTES----------------------------------------------------------------

 
def Initialization_Plaquette(Lattice_Links, NN_Map, Lattice_Plaq):

    Plaq = []
    plaq_position = 0

    for i in range(len(Lattice_Links)):

        x_p_pos = NN_Map[i][0]
        y_p_pos = NN_Map[i][1]
        z_p_pos = NN_Map[i][2]

        plaq = Lattice_Links[i][0] * Lattice_Links[i][1] * Lattice_Links[y_p_pos][0] * Lattice_Links[x_p_pos][1] 
        Lattice_Plaq[i][0].append(plaq_position)
        Lattice_Plaq[i][1].append(plaq_position)
        Lattice_Plaq[y_p_pos][0].append(plaq_position)
        Lattice_Plaq[x_p_pos][1].append(plaq_position)
        Plaq.append(plaq)
        plaq_position += 1

        plaq = Lattice_Links[i][0] * Lattice_Links[i][2] * Lattice_Links[z_p_pos][0] * Lattice_Links[x_p_pos][2] 
        Lattice_Plaq[i][0].append(plaq_position)
        Lattice_Plaq[i][2].append(plaq_position)
        Lattice_Plaq[z_p_pos][0].append(plaq_position)
        Lattice_Plaq[x_p_pos][2].append(plaq_position)
        Plaq.append(plaq)
        plaq_position += 1

        plaq = Lattice_Links[i][1] * Lattice_Links[i][2] * Lattice_Links[y_p_pos][2] * Lattice_Links[z_p_pos][1] 
        Lattice_Plaq[i][1].append(plaq_position)
        Lattice_Plaq[i][2].append(plaq_position)
        Lattice_Plaq[y_p_pos][2].append(plaq_position)
        Lattice_Plaq[z_p_pos][1].append(plaq_position)
        Plaq.append(plaq)
        plaq_position += 1
    
    return Plaq


#--------------------------------------------------WANG-LANDAU METHOD--------------------------------------------------------------------

# Function to check if the histogram is "flat", i.e. if the 80% of the data-points are above the average
def is_flat(H, threshold=0.85):
    nonzero = H[H > 0]
    if len(nonzero) == 0:
        return False
    min_H = np.min(nonzero)
    avg_H = np.mean(nonzero)
    return min_H >= threshold * avg_H

def Wang_Landau(Lattice_Links, Lattice_Plaq, Plaq, n, nt, f_in, f_threshold, offset):
    Lattice_Size = n * n * nt
    num_energy_levels = 2 * offset

    H = np.zeros(num_energy_levels)
    log_g = np.zeros(num_energy_levels)
    f = f_in

    current_energy = compute_total_energy(Plaq)

    while f > f_threshold:
        H.fill(0)
        print(f"Current modification factor f: {f}")
        while True:
            for _ in range(1000):
                i = random.randint(0, Lattice_Size - 1)

                for direction in range(3):
                    dE = sum(Plaq[p] for p in Lattice_Plaq[i][direction])
                    
                    new_energy = current_energy + 2*dE
                    e_old_idx = int(current_energy + offset)
                    e_new_idx = int(new_energy + offset)
                    
                    if not (0 <= e_new_idx < len(log_g)) or not (0 <= e_old_idx < len(log_g)):
                        continue

                    delta_log_g = log_g[e_old_idx] - log_g[e_new_idx]
                                
                    try:
                        accept_prob = min(1.0, math.exp(delta_log_g))
                    except OverflowError:
                        accept_prob = 1.0
                                
                    if random.random() < accept_prob:
                        current_energy = new_energy
                        current_index = e_new_idx
                        Lattice_Links[i][direction] *= -1
                        for element in range( len(Lattice_Plaq[i][direction]) ):
                            Plaq[ Lattice_Plaq[i][direction][element]] *= -1
                                
                    else:
                        current_index = e_old_idx
                                
                    log_g[current_index] += math.log(f)
                    H[current_index] += 1
                    
            
            if is_flat(H):
                f = math.sqrt(f)
                break
            
            else:
                continue

    return log_g



# Function to perform the energy in the Wang-Landau algorithm. In this function, energy is the sum of the product of the near-neighbour spins.
def compute_total_energy(Plaq):
    
    return -np.sum(Plaq)

#--------------------------------------------------ENERGY------------------------------------------------------------------

def mean_energy_from_log_g(log_g, beta, offset):
    energies = np.arange(len(log_g)) - offset 
    log_weights = log_g - beta * energies

    valid = log_g > 0
    log_weights = log_weights[valid]
    energies = energies[valid]

    log_weights -= np.max(log_weights)

    weights = np.exp(log_weights)
    Z = np.sum(weights)
    E_avg = np.sum(energies * weights) / Z
    return E_avg
 
def mean_energy2_from_log_g(log_g, beta, offset):
    energies = np.arange(len(log_g)) - offset 
    log_weights = log_g - beta * energies

    valid = log_g > 0
    log_weights = log_weights[valid]
    energies = energies[valid]

    log_weights -= np.max(log_weights)

    weights = np.exp(log_weights)
    Z = np.sum(weights)
    E2_avg = np.sum(weights * energies **2) / Z
    return E2_avg

#--------------------------------------------------MAIN--------------------------------------------------------------------
for _ in range(1):
    All_Energies = []
    All_Heat = []
    All_Variance = []

    for simulation in range(Num_Simulations):

        Lattice_Size= Height * Lenght**2
        offset = int(3 * Lattice_Size)

        SpecificHeat = []
        Average_Energy = []
        Average_Energy2 = []
        Variance = []

        t=0.5
        K_b = 1
        beta = 1
        Temperature = []
        Lattice_Links = []
        Lattice_Plaq = []
        NN_Map = []
        Plaq = []
        
        Lattice_Links, Lattice_Plaq, NN_Map = Initialization_Lattice()

        output_file = os.path.join("Geometry_test.txt")
        file_exists = os.path.isfile(output_file)
   
        output_file = "Geometry_test.txt"
        with open(output_file, "w") as f:
            for i in range(len(Lattice_Links)):
                site_coord = get_coordinates(i)
                f.write(f"Site {i} - Coordinate: {site_coord} in a Lattice {Lenght}^2 x {Height}\n")
                f.write("Near neighbour:\n")
                for nn_index in NN_Map[i]:
                    nn_coord = get_coordinates(nn_index)
                    f.write(f"{nn_coord}\n")
                f.write("\n")

        
        Plaq = Initialization_Plaquette(Lattice_Links, NN_Map, Lattice_Plaq)
        f = np.exp(1.0)
        Threshold = 1.0005
        log_g = Wang_Landau(Lattice_Links, Lattice_Plaq, Plaq, Lenght, Height, f, Threshold, offset)
        plt.plot(log_g)
        plt.show()

        for i in range(60):

            beta = 1.0 / ( K_b * t)
            E_avg = mean_energy_from_log_g(log_g, beta, offset)
            
            print(f"Total Energy = {beta}: {E_avg}")
            E2_avg = mean_energy2_from_log_g(log_g, beta, offset) 
            c = E2_avg - E_avg**2
            Variance.append(c)
            c = c / (K_b * t**2)
            Average_Energy2.append(E2_avg)
            Average_Energy.append(E_avg / (Lattice_Size))
            SpecificHeat.append(c / (Lattice_Size))
            Temperature.append(t)
            t += 0.06
        
        All_Energies.append(Average_Energy)
        All_Heat.append(SpecificHeat)
        All_Variance.append(Variance)



    for i, simulation in enumerate(range(Num_Simulations)):
        plt.plot(Temperature, All_Heat[i], label=f"{Lenght}^2x{Height}")
    plt.xlabel('Temperature')
    plt.ylabel('Specific Heat per Spin')
    plt.title('Specific Heat for different lattice sizes')
    plt.legend()
    plt.grid(True)
    plt.show()

    for i, simulation in enumerate(range(Num_Simulations)):
        plt.plot(Temperature, All_Energies[i], label=f"{Lenght}^2x{Height}")
    plt.xlabel('Temperature')
    plt.ylabel('Average Energy per spin')
    plt.title('Average Energy per spin for different lattice sizes')
    plt.legend()
    plt.grid(True)
    plt.show()

    output_file = os.path.join("specific_heat_data.txt")
    file_exists = os.path.isfile(output_file)
    with open(output_file, "a") as f:
        if not file_exists:
            f.write("# Simulation\tTemperature\tSpecificHeat\n")
        for i, simulation in enumerate(range(Num_Simulations)):
            for j in range(len(Temperature)):
                line = f"{simulation}\t{Temperature[j]:.4f}\t{All_Heat[i][j]:.6f}\n"
                f.write(line)

    output_file = os.path.join("Variance_energy_data.txt")
    file_exists = os.path.isfile(output_file)
    with open(output_file, "a") as f:
        if not file_exists:
            f.write("# Simulation\tTemperature\tVariance\n")
        for i, simulation in enumerate(range(Num_Simulations)):
            for j in range(len(Temperature)):
                line = f"{simulation}\t{Temperature[j]:.4f}\t{All_Variance[i][j]:.6f}\n"
                f.write(line)



