from Bio.PDB import *

def obtain_types_and_coordinates(pdb_file_path):
    '''Generates a list containing the residues' types and coordinates'''

    # Generation of a dictionary for the chemical type associated for each amino acid:
    aa_dictionary = {
        "ARG": "pos", "HIS": "pos", "LYS": "pos", "ASP": "neg", "GLU": "neg",
        "SER": "polar", "THR": "polar", "ASN": "polar", "GLN": "polar",
        "CYS": "special", "SEC": "special", "GLY": "special", "PRO": "special",
        "ALA": "hydrophobic", "VAL": "hydrophobic", "ILE": "hydrophobic",
        "LEU": "hydrophobic", "MET": "hydrophobic", "PHE": "hydrophobic",
        "TYR": "hydrophobic", "TRP": "hydrophobic"
    }

    # The goal of this function is to:
    # Generate a list (lista_res) where each element is itself a list containing
    # four elements: i) the type of residue and ii), iii), and iv) the OX, OY, and OZ
    # coordinates of the alpha carbon (CA) atom:
    lista_res = []

    # Preparation for the use of Bio.PDB module:
    parser = PDBParser()

    # Obtain the structure object:
    estructura = parser.get_structure('X', pdb_file_path)

    # Obtain the residues object:
    residuos = estructura.get_residues()

    # For every residue:
    for res in residuos:
        # Get the residue name (e.g., ALA, GLY, etc.)
        residue_name = res.get_resname()

        # Check if the residue name is in the amino acid dictionary
        if residue_name in aa_dictionary:
            # Get the chemical type of the residue
            residue_type = aa_dictionary[residue_name]

            # Try to get the alpha carbon (CA) atom
            try:
                ca_atom = res['CA']
                ca_coordinates = ca_atom.get_coord()
                # Append the residue type and CA coordinates to lista_res
                lista_res.append([residue_type, ca_coordinates[0], ca_coordinates[1], ca_coordinates[2]])
            except KeyError:
                # If there is no CA atom, skip this residue
                continue

    return lista_res
