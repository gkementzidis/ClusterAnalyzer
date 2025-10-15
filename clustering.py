import os
import json
import numpy as np

def map_atom_to_mol(file):
    """
    Maps atom IDs to their corresponding molecule IDs and types.

    Input:
    - file: Path to the input file containing atom and molecule information.

    Output:
    - mol_id_maps: List mapping atom IDs to their molecule IDs.
    - n_atoms: Total number of atoms.
    - n_mol: Total number of molecules.
    - atom_types: List mapping atom IDs to their types.
    """
    with open(file, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[2].split()[0])

    mol_id_maps = [0 for _ in range(n_atoms)]
    atom_types = [0 for _ in range(n_atoms)]

    in_atoms_section = False
    first = False

    for line in lines:
        # Detect start of Atoms section
        if line.strip().startswith("Atoms # full"):
            in_atoms_section = True
            continue

        # Keep empty line right after "Atoms"
        if in_atoms_section and line.strip() == "" and not first:
            first = True
            continue
        elif in_atoms_section and line.strip() == "" and first:
            in_atoms_section = False

        # If in Atoms section and line looks like atom data
        if in_atoms_section and line.strip() and line[0].isdigit():
            parts = line.split()
            # parts: [atom-ID, mol-ID, atom-type, x, y, z]
            atom_type = int(parts[2])
            atom_id = int(parts[0])
            mol_id = int(parts[1])
            
            # WARNING -- IDs from LAMMPS are 1-indexed
            mol_id_maps[atom_id - 1] = mol_id - 1
            atom_types[atom_id - 1] = atom_type

    n_mol = max(mol_id_maps) + 1
    
    return mol_id_maps, n_atoms, n_mol, atom_types

def map_mol_to_atom(n_mol, mol_id_map):
    """
    Maps molecule IDs to a list of the atom IDs it consists of.

    Inputs:
    - n_mol: number of molecules
    - mol_id_map: a list of size (# atoms) such that mol_id_map[i] is the molecule id if the i-th atom

    Output:
    - mol_to_atom_map: a list of lists n_mol x (# beads/atoms in the molecule)
    """
    mol_to_atom_map = [[] for _ in range(n_mol)]
    for atom_id, mol_id in enumerate(mol_id_map):
        mol_to_atom_map[mol_id].append(atom_id)
    return mol_to_atom_map

def identify_bonds(lines, n_atoms, n_mol, mol_id_maps):
    """
    Identifies bonds between atoms and creates an adjacency matrix.

    Input:
    - lines: List of lines from the input file.
    - n_atoms: Total number of atoms.
    - n_mol: Total number of molecules.
    - mol_id_maps: List mapping atom IDs to their molecule IDs.

    Output:
    - adj_mat: Adjacency matrix representing bonds between molecules.
    """

    adj_mat = np.zeros((n_mol, n_mol), dtype=int)

    in_bonds_section = False
    first = False

    for line in lines:
        # Detect start of Bonds section
        if line.strip().startswith("Bonds"):
            in_bonds_section = True
            continue

        # Keep empty line right after "Bonds"
        if in_bonds_section and line.strip() == "" and not first:
            first = True
            continue
        elif in_bonds_section and line.strip() == "" and first:
            in_bonds_section = False

        # If in Bonds section and line looks like bond data
        if in_bonds_section and line.strip() and line[0].isdigit():
            parts = line.split()
            # parts: [bond-ID, bond-type, atom1-ID, atom2-ID]
            bond_type = int(parts[1])
            atom1 = int(parts[2])
            atom2 = int(parts[3])

            if bond_type == interm_bond_types[0] or bond_type == interm_bond_types[1]:
                # WARNING -- IDs from LAMMPS are 1-indexed
                adj_mat[mol_id_maps[atom1 - 1], mol_id_maps[atom2 - 1]] = 1
                adj_mat[mol_id_maps[atom2 - 1], mol_id_maps[atom1 - 1]] = 1

    return adj_mat

def extract_coordinates(lines, n_atoms):
    """
    Extracts atomic coordinates from the input lines.

    Input:
    - lines: List of lines from the input file.
    - n_atoms: Total number of atoms.

    Output:
    - coords: Array of shape (n_atoms, 3) containing the x, y, z coordinates of each atom.
    """

    coords = np.zeros((n_atoms, 3), dtype=float)

    in_atoms_section = False
    first = False

    for line in lines:
        # Detect start of Atoms section
        if line.strip().startswith("Atoms # full"):
            in_atoms_section = True
            continue

        # Keep empty line right after "Atoms"
        if in_atoms_section and line.strip() == "" and not first:
            first = True
            continue
        elif in_atoms_section and line.strip() == "" and first:
            in_atoms_section = False

        # If in Atoms section and line looks like atom data
        if in_atoms_section and line.strip() and line[0].isdigit():
            parts = line.split()
            # parts: [atom-ID, mol-ID, atom-type, x, y, z]
            atom_id = int(parts[0])
            x = float(parts[3])
            y = float(parts[4])
            z = float(parts[5])

            # WARNING -- IDs from LAMMPS are 1-indexed
            coords[atom_id - 1] = [x, y, z]

    return coords


def get_connected_components(mat, step):
    """
    Identifies connected components in a molecular adjacency matrix.

    Input:
    - mat: Adjacency matrix representing bonds between molecules.
    - step: Current simulation step (used for debugging/visualization).

    Output:
    - num_components: Number of connected components found.
    - component_sizes: List of sizes of each connected component.
    - components: List of lists, where each inner list contains the indices of atoms in that component.
    - diameters: List of diameters for each connected component.
    """
    visited = [False] * len(mat)
    components = []
    diameters = []

    for i in range(len(mat)):
        if not visited[i]:
            component = []
            stack = [i]
            while stack:
                node = stack.pop()
                if not visited[node]:
                    visited[node] = True
                    component.append(node)
                    for neighbor, connected in enumerate(mat[node]):
                        if connected and not visited[neighbor]:
                            stack.append(neighbor)
                            
            # compute diameter (maximum distance between any two ATOMS between molecules within a component)
            if component:
                component_distances = []
                for ind1 in range(len(component)):
                    for ind2 in range(ind1, len(component)):
                        mol1 = component[ind1]
                        mol2 = component[ind2]
                        for ind3 in range(len(mol_to_atom_map[mol1])):
                            for ind4 in range(len(mol_to_atom_map[mol2])):
                                atom1 = mol_to_atom_map[mol1][ind3]
                                atom2 = mol_to_atom_map[mol2][ind4]
                                if atom_types[atom1] not in alphaC_beads and atom_types[atom2] not in alphaC_beads and atom1 > atom2:
                                    # distance = np.linalg.norm(coordinates[step, atom1] - coordinates[step, atom2]) / 10
                                    
                                    # accounting for Periodic Boundary Conditions ???
                                    dx, dy, dz = (coordinates[step, atom1] - coordinates[step, atom2]) / 10
                                    
                                    # FIX: PERIODIC BC ONLY IF THE DISTANCE CHANGED BY A LOT -- SOME DISTANCES CAN BE AS BIG AS 300nm OR MORE
                                    dx = dx if abs(dx) < 300 else abs(dx) - 600
                                    dy = dy if abs(dy) < 150 else abs(dy) - 300
                                    dz = dz if abs(dz) < 150 else abs(dz) - 300
                                    distance = np.sqrt(dx**2 + dy**2 + dz**2)
                                    
                                    component_distances.append(distance)
                if component_distances:
                    diameters.append(max(component_distances))
            components.append(component)

    return len(components), [len(c) for c in components], components, diameters


if __name__ == "__main__":
    specs = json.load(open("./configs/sample_config.json"))
    print(specs["name"])

    # Identfy global variables from specs
    interm_bond_types = specs["intermolecular_bond_types"]
    # in our fibrin monomer model, we want to exclude the alphaC beads from the diameter calculation
    # to ensure a fair comparison with des-alphaC monomers 
    if "extra" in specs and "exclude_beads_from_diameter" in specs["extra"]:
        alphaC_beads = specs["extra"]["exclude_beads_from_diameter"]

    # Input / output directory and filenames
    data_dir = specs["paths"]["data_dir"]
    if data_dir[-1] != '/':
        data_dir += '/'
    output_dir = specs["paths"]["output_dir"]
    if output_dir[-1] != '/':
        output_dir += '/'
    output_file = specs["filenames"]["output"]

    print("Reading data from {}".format(data_dir))
    file_names = []
    for file in os.listdir(data_dir):
        if file.endswith(".data"):
            file_names.append(file)

    # sort the file names according to the numeric value in the filename
    file_names.sort(key=lambda x: int(x.split('-')[1].split('.')[0]))
    n_steps = len(file_names)
    print("Found {} data files, and sorted them.".format(n_steps))

    # ---------------------------------------------------------------- # 

    mol_id_maps, n_atoms, n_mol, atom_types = map_atom_to_mol(os.path.join(data_dir, file_names[0]))
    print("Number of atoms: {}".format(n_atoms))
    print("Number of molecules: {}".format(n_mol))

    mol_to_atom_map = map_mol_to_atom(n_mol, mol_id_maps)

    # overall adjacency matrix (all steps)
    all_adj_mat = np.zeros((n_steps, n_mol, n_mol), dtype=int)
    # overall coordinates (all steps)
    coordinates = np.zeros((n_steps, n_atoms, 3), dtype=float)

    for i, file in enumerate(file_names):
        with open(os.path.join(data_dir, file), "r") as f:
            lines = f.readlines()

        all_adj_mat[i] = identify_bonds(lines, n_atoms, n_mol, mol_id_maps)
        coordinates[i] = extract_coordinates(lines, n_atoms)

    # we iterate through every step
    # we want to get the number of connected components and the number of nodes in each component
    n_components = np.zeros(n_steps, dtype=int)
    component_sizes = [None for _ in range(n_steps)]
    diameters = [None for _ in range(n_steps)]

    for step in range(n_steps):
        adj_mat = all_adj_mat[step]
        n_component, component_size, components, diameter = get_connected_components(adj_mat, step)
        diameters[step] = diameter
        n_components[step] = n_component
        component_sizes[step] = component_size

    # ---------------------------------------------------------------- # 

    avg_component_size = [np.mean(size) for size in component_sizes]
    avg_diameters = [np.mean(d) for d in diameters]
    # instead of the mean of all clusters, get the mean of the highest 20%
    # remember, I want the MEAN of the TOP 20% of the sizes, not the singular value of the 80th percentile
    # remember, np.percentile returns a singular value, not a list of values
    # avg_component_size_20 = [np.mean(size[size > np.percentile(size, 80)]) for size in component_sizes]
    avg_component_size_20 = []
    for arr in component_sizes:
        a = np.asarray(arr)                 # ensure numpy array
        if a.size == 0:
            avg_component_size_20.append(np.nan)
            continue
        thr = np.percentile(a, 80)
        avg_component_size_20.append(a[a >= thr].mean())

    diameters_20 = []
    for arr in diameters:
        a = np.asarray(arr)
        if a.size == 0:
            diameters_20.append(np.nan)
            continue
        thr = np.percentile(a, 80)
        diameters_20.append(a[a >= thr].mean())

    # same for diameters
    # diameters_20 = [np.mean(d[d > np.percentile(d, 80)]) for d in diameters]
    diameters_max = [np.max(d) for d in diameters]

    # SAVE RESULTS TO A FILE
    # not a txt file
    with open(output_dir + output_file, "w") as f:
        json.dump({
            "avg_component_size": avg_component_size,
            "avg_diameters": avg_diameters,
            "avg_component_size_20": avg_component_size_20,
            "diameters_20": diameters_20,
            "diameters_max": diameters_max
        }, f)