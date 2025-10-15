"""
Functions used at the beginning of the analysis to create correspondence maps between molecule and "atom" (CG bead) IDs.
"""

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
