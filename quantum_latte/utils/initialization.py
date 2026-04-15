import numpy as np

def initialize_grid(num_particles, box_size):
    n_side = int(np.ceil(np.sqrt(num_particles)))
    spacing = box_size / n_side
    
    positions = []
    for i in range(n_side):
        for j in range(n_side):
            if len(positions) < num_particles:
                x = i * spacing + spacing / 2
                y = j * spacing + spacing / 2
                positions.append([x, y])
                
    return np.array(positions)