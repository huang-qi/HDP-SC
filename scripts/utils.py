#ignore warnings
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import json
import os


def read_last_timesteps(file_path, n_timesteps=2):
    # Define block size for reading
    block_size = 1024 * 1024  # 1 MB

    # Open file
    with open(file_path, 'rb') as file:
        # Jump to last block of file
        file.seek(0, 2)  # Move to end of file
        file_size = file.tell()
        iterations = 0
        buffer = b''
        atom_attributes = None

        while file.tell() > 0:
            # Calculate next read position and size
            cursor = max(file_size - (iterations + 1) * block_size, 0)
            size_to_read = min(block_size, file.tell() - cursor)
            file.seek(cursor)
            
            # Read data block
            data = file.read(size_to_read) + buffer
            
            # Find TIMESTEP separators and ATOMS attributes
            timesteps_found = data.count(b'TIMESTEP')
            if atom_attributes is None:
                # Check ATOMS line
                atoms_index = data.find(b'ITEM: ATOMS')
                if atoms_index != -1:
                    end_of_line = data.find(b'\n', atoms_index)
                    atom_line = data[atoms_index:end_of_line].decode('utf-8').strip()
                    atom_attributes = atom_line.split(" ")[2:]
            
            # If enough TIMESTEPS found, process these blocks
            if timesteps_found >= n_timesteps:
                # Find TIMESTEP positions in reverse
                positions = []
                last_pos = len(data)
                while len(positions) < n_timesteps and last_pos != -1:
                    last_pos = data.rfind(b'TIMESTEP', 0, last_pos)
                    if last_pos != -1:
                        positions.append(last_pos)
                        last_pos -= 1
                
                # Process each TIMESTEP block from back to front
                results = []
                for start in reversed(positions):
                    # Find end of TIMESTEP line (start of timestep number line)
                    end_of_line = data.find(b'\n', start + 8) + 1
                    # Find end of timestep number line
                    end_of_timestep = data.find(b'\n', end_of_line)
                    timestep = data[end_of_line:end_of_timestep].decode('utf-8').strip()
                    results.append(timestep)
                
                return results, atom_attributes
            
            # Update buffer with unprocessed data
            buffer = data[:size_to_read]
            iterations += 1
        
        return [], atom_attributes

def generate_bonds_array(build_degree, build_number):
    num_bonds = build_number*(build_degree-1)
    bonds_array = np.zeros((num_bonds, 2))
    for i in range(build_number):
        for j in range(build_degree-1):
            bonds_array[i*(build_degree-1)+j] = [i*build_degree+j, i*build_degree+j+1]
    return bonds_array

def generate_bonds_array_AB(build_degree, build_number):
    num_bonds = build_number * (build_degree * 2 - 1)
    
    2+3*(build_degree-2)
    
    bonds_array = np.zeros((num_bonds, 2), dtype=int)
    bond_index = 0
    for i in range(build_number):
        for j in range(build_degree):
            A_id = (i * build_degree + j) * 2
            B_id = A_id + 1
            bonds_array[bond_index] = [A_id, B_id]
            bond_index += 1
            if j < build_degree - 1:
                next_A_id = A_id + 2
                bonds_array[bond_index] = [A_id, next_A_id]
                bond_index += 1
    return bonds_array

def generate_angles_array(build_degree, build_number):
    num_angles = build_number*(build_degree-2)
    angles_array = np.zeros((num_angles, 3))
    for i in range(build_number):
        for j in range(build_degree-2):
            angles_array[i*(build_degree-2)+j] = [i*build_degree+j, i*build_degree+j+1, i*build_degree+j+2]
    return angles_array

def generate_angles_array_AB(build_degree, build_number):
    num_angles = build_number * (build_degree * 3 - 4)
    angles_array = np.zeros((num_angles, 3), dtype=int)
    angle_index = 0
    for i in range(build_number):
        for j in range(build_degree):
            A_id = (i * build_degree + j) * 2
            B_id = A_id + 1
            if j < build_degree - 1:
                next_A_id = A_id + 2
                next_B_id = next_A_id + 1
                # Angle formed by A, next_A, next_next_A
                if j < build_degree - 2:
                    next_next_A_id = next_A_id + 2
                    angles_array[angle_index] = [A_id, next_A_id, next_next_A_id]
                    angle_index += 1
                # Angle formed by B, A, next_A
                angles_array[angle_index] = [B_id, A_id, next_A_id]
                angle_index += 1
                # Angle formed by A, next_A, next_B
                angles_array[angle_index] = [A_id, next_A_id, next_B_id]
                angle_index += 1
    return angles_array

def generate_dihedrals_array(build_degree, build_number):
    num_dihedrals = build_number*(build_degree-3)
    dihedrals_array = np.zeros((num_dihedrals, 4))
    for i in range(build_number):
        for j in range(build_degree-3):
            dihedrals_array[i*(build_degree-3)+j] = [i*build_degree+j, i*build_degree+j+1, i*build_degree+j+2, i*build_degree+j+3]
    return dihedrals_array

def generate_dihedrals_array_AB(build_degree, build_number):
    num_dihedrals = build_number * ((build_degree - 1) + (build_degree - 2) + (build_degree - 3))
    dihedrals_array = np.zeros((num_dihedrals, 4), dtype=int)
    dihedral_index = 0
    for i in range(build_number):
        for j in range(build_degree):
            A_id = (i * build_degree + j) * 2
            B_id = A_id + 1
            if j < build_degree - 1:
                next_A_id = A_id + 2
                next_B_id = next_A_id + 1
                # Dihedral formed by B, A, next_A, next_B
                dihedrals_array[dihedral_index] = [B_id, A_id, next_A_id, next_B_id]
                dihedral_index += 1
                if j < build_degree - 2:
                    next_next_A_id = next_A_id + 2
                    next_next_B_id = next_next_A_id + 1
                    # Dihedral formed by A, next_A, next_next_A, next_next_B
                    dihedrals_array[dihedral_index] = [A_id, next_A_id, next_next_A_id, next_next_B_id]
                    dihedral_index += 1
                    if j < build_degree - 3:
                        next_next_next_A_id = next_next_A_id + 2
                        # Dihedral formed by A, next_A, next_next_A, next_next_next_A
                        dihedrals_array[dihedral_index] = [A_id, next_A_id, next_next_A_id, next_next_next_A_id]
                        dihedral_index += 1
    return dihedrals_array

def parse_aa_data(file_path):
    num_atoms = None
    num_atom_types = None
    box_size = []
    mass_map = {}

    with open(file_path, "r") as data:
        for i, line in enumerate(data):
            if line.endswith("atoms\n"):
                num_atoms = int(line.split()[0])
                print("num_atoms:", num_atoms)

            if line.endswith("atom types\n"):
                num_atom_types = int(line.split()[0])
                print("num_atom_types:", num_atom_types)

            if line.endswith("xlo xhi\n"):
                box_size = [float(x) for x in line.split()[:2]]
                for _ in range(2):
                    line = next(data)
                    box_size += [float(x) for x in line.split()[:2]]
                print("box_size:", box_size)

            if line.startswith("Masses"):
                next(data)  # Skip the "Masses" line
                mass_map = {}
                for _ in range(num_atom_types):
                    line = next(data)
                    mass_map[int(line.split()[0])] = float(line.split()[1])
                print("mass_map:", mass_map)
                break  # No need to continue reading the file

    return num_atoms, num_atom_types, box_size, mass_map

def generate_psf_file(num_chains, residues_per_chain, mass):
    """
    Generates a PSF file for a given number of chains and residues per chain.

    Args:
    num_chains (int): Number of chains.
    residues_per_chain (int): Number of residues per chain.
    mass (float): Mass of each atom.

    Returns:
    str: A string representing the content of the PSF file.
    """

    # Header
    psf_content = "PSF NAMD\n\n"
    
    # Total number of atoms
    total_atoms = num_chains * residues_per_chain
    psf_content += f"{total_atoms:8} !NATOM\n"

    # Atom section
    atom_id = 1
    for chain_id in range(num_chains):
        for residue_id in range(residues_per_chain):
            psf_content += f"{atom_id:10}{' '}{chain_id:<10}{residue_id + 1:<10}{'1':<10}{'0':<10}{'0':<10}0.000000   {mass}           0\n"
            atom_id += 1

    # Bonds section
    total_bonds = total_atoms - num_chains  # One less bond per chain
    psf_content += f"\n{total_bonds:10} !NBOND: bonds\n"
    bond_lines = []
    for i in range(1, total_atoms):
        if i % residues_per_chain != 0:  # Avoid bonding between chains
            bond_lines.append(f"{i:<10}{i + 1:<10}")
    psf_content += '\n'.join(['\t'.join(bond_lines[i:i+4]) for i in range(0, len(bond_lines), 4)]) + "\n"

    # Angles section
    total_angles = max(0, total_atoms - 2 * num_chains)  # Two less angles per chain
    psf_content += f"\n{total_angles:10} !NTHETA: angles\n"
    angle_lines = []
    for i in range(1, total_atoms - 1):
        if i % residues_per_chain != 0 and i % residues_per_chain != residues_per_chain - 1:  # Avoid angles between chains
            angle_lines.append(f"{i:<10}{i + 1:<10}{i + 2:<10}")
    psf_content += '\n'.join(['\t'.join(angle_lines[i:i+4]) for i in range(0, len(angle_lines), 4)]) + "\n"

    # Dihedrals section
    total_dihedrals = max(0, total_atoms - 3 * num_chains)  # Three less dihedrals per chain
    psf_content += f"\n{total_dihedrals:10} !NPHI: dihedrals\n"
    dihedral_lines = []
    for i in range(1, total_atoms - 2):
        if i % residues_per_chain != 0 and i % residues_per_chain != residues_per_chain - 1 and i % residues_per_chain != residues_per_chain - 2:  # Avoid dihedrals between chains
            dihedral_lines.append(f"{i:<10}{i + 1:<10}{i + 2:<10}{i + 3:<10}")
    psf_content += '\n'.join(['\t'.join(dihedral_lines[i:i+4]) for i in range(0, len(dihedral_lines), 4)]) + "\n"

    # Impropers section (can be added similarly if needed)
    total_impropers = 0  # Modify this if you have data for impropers
    psf_content += f"\n{total_impropers:10} !NIMPHI: impropers\n"
    
    # Other sections with zero entries
    psf_content += "\n         0 !NDON: donors\n\n         0 !NACC: acceptors\n\n         0 !NNB\n\n         0          0 !NGRP\n"

    return psf_content

def generate_psf_file_AB(build_degree, build_number, massA, massB):
    """
    Generates a PSF file for a given number of builds with a specified degree.

    Args:
    build_degree (int): Degree of each build.
    build_number (int): Number of builds.
    mass (float): Mass of each atom.

    Returns:
    str: A string representing the content of the PSF file.
    """

    # Header
    psf_content = "PSF NAMD\n\n"
    
    # Total number of atoms
    total_atoms = build_number * build_degree * 2
    psf_content += f"{total_atoms:8} !NATOM\n"

    # Atom section
    atom_id = 1
    for build_id in range(build_number):
        for degree_id in range(build_degree):
            A_id = atom_id
            B_id = atom_id + 1
            psf_content += f"{A_id:10}{' '}{build_id:<10}{degree_id + 1:<10}{'A':<10}{'A':<10}{'A':<10}0.000000   {massA}           0\n"
            psf_content += f"{B_id:10}{' '}{build_id:<10}{degree_id + 1:<10}{'B':<10}{'B':<10}{'B':<10}0.000000   {massB}           0\n"
            atom_id += 2

    # Bonds section
    bonds_array = generate_bonds_array_AB(build_degree, build_number)
    total_bonds = len(bonds_array)
    psf_content += f"\n{total_bonds:10} !NBOND: bonds\n"
    bond_lines = []
    for bond in bonds_array:
        bond_lines.append(f"{bond[0] + 1:<10}{bond[1] + 1:<10}")
    psf_content += '\n'.join(['\t'.join(bond_lines[i:i+4]) for i in range(0, len(bond_lines), 4)]) + "\n"

    # Angles section
    angles_array = generate_angles_array_AB(build_degree, build_number)
    total_angles = len(angles_array)
    psf_content += f"\n{total_angles:10} !NTHETA: angles\n"
    angle_lines = []
    for angle in angles_array:
        angle_lines.append(f"{angle[0] + 1:<10}{angle[1] + 1:<10}{angle[2] + 1:<10}")
    psf_content += '\n'.join(['\t'.join(angle_lines[i:i+4]) for i in range(0, len(angle_lines), 4)]) + "\n"

    # Dihedrals section
    dihedrals_array = generate_dihedrals_array_AB(build_degree, build_number)
    total_dihedrals = len(dihedrals_array)
    psf_content += f"\n{total_dihedrals:10} !NPHI: dihedrals\n"
    dihedral_lines = []
    for dihedral in dihedrals_array:
        dihedral_lines.append(f"{dihedral[0] + 1:<10}{dihedral[1] + 1:<10}{dihedral[2] + 1:<10}{dihedral[3] + 1:<10}")
    psf_content += '\n'.join(['\t'.join(dihedral_lines[i:i+4]) for i in range(0, len(dihedral_lines), 4)]) + "\n"

    # Impropers section (can be added similarly if needed)
    total_impropers = 0  # Modify this if you have data for impropers
    psf_content += f"\n{total_impropers:10} !NIMPHI: impropers\n"
    
    # Other sections with zero entries
    psf_content += "\n         0 !NDON: donors\n\n         0 !NACC: acceptors\n\n         0 !NNB\n\n         0          0 !NGRP\n"

    return psf_content

#print a dict as tree structure
def print_dict(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            print_dict(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))
