#This script receives as parameter the path to the pdb_entry_type.txt, which
#is a file that contains in each line a pdb entry id, the type of biomolecule and
#the technique used to obtain the structure, i.e:
#100d	nuc	diffraction
#101d	nuc	diffraction
#101m	prot	diffraction
#102d	nuc	diffraction
#102l	prot	diffraction
#102m	prot	diffraction
#103d	nuc	NMR
#103l	prot	diffraction
#...
#The script proceeds as follows:
#0)Creates a file named "pdb_id_and_obtained_Q.txt".
#1)Reads a line (starting with the first) of the pdb_entry_type.txt file.
#2)If the line is that corresponding to a protein --i.e if the line
#contains the word "prot"-- then the pdb id contained in that line
#is downloaded from the RCSB database.
#3)The obtain_Q function from
#calculate_Q.py if invoked with a try()
#statement, passing as argument the route of the downloaded pdb file.
#4)If the obtain_Q function works, then the pdb id
#and its corresponding Q value are appended to 
#pdb_id_and_obtained_Q.txt .
#5)Proceed in the same manner with the rest of lines of the
#pdb_entry_type.txt file.
#6)The script also generates a file named pdb_id_residue_percentages.txt
#which contains the residue percentage of each type of amino acid for
#each pdb id

import os
import sys
import urllib.request
from calculate_Q_and_len import obtain_Q_and_len
from Bio.PDB import PDBParser

# Auxiliary function
def download_pdb(pdb_id):
    url = "https://files.rcsb.org/view/" + pdb_id + ".pdb"
    file_name = pdb_id + ".pdb"
    urllib.request.urlretrieve(url, file_name)

# Auxiliary function
def remove_pdb_file(pdb_id, suffix=""):
    """Removes a PDB file with the given pdb_id and optional suffix."""
    file_name = pdb_id + suffix + ".pdb"
    try:
        os.remove(file_name)
        print(f"File {file_name} successfully removed.")
    except FileNotFoundError:
        print(f"File {file_name} not found.")
    except Exception as e:
        print(f"An error occurred while trying to remove {file_name}: {e}")

# Dictionary for residue classification
aa_dictionary = {
    "ARG": "pos", "HIS": "pos", "LYS": "pos", "ASP": "neg", "GLU": "neg",
    "SER": "polar", "THR": "polar", "ASN": "polar", "GLN": "polar",
    "CYS": "special", "SEC": "special", "GLY": "special", "PRO": "special",
    "ALA": "hydrophobic", "VAL": "hydrophobic", "ILE": "hydrophobic",
    "LEU": "hydrophobic", "MET": "hydrophobic", "PHE": "hydrophobic",
    "TYR": "hydrophobic", "TRP": "hydrophobic"
}

# Function to extract sequence from PDB file using BioPython
def extract_sequence_from_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    sequence = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in aa_dictionary:
                    sequence.append(residue.get_resname())
              
    return sequence

# Function to calculate residue percentages
def calculate_residue_percentages(sequence):
    total_residues = len(sequence)
    percentages = {}  # Initialize an empty dictionary to store percentages

    counts = {"hydrophobic": 0, "polar": 0, "pos": 0, "neg": 0, "special": 0}

    for residue in sequence:
        residue_type = aa_dictionary.get(residue, None)
        if residue_type:
            counts[residue_type] += 1

    # Calculate the percentages
    for key in counts:
        percentage = (counts[key] / total_residues) * 100
        percentages[key] = percentage

    return percentages

#0) Creates output files and defines auxiliary lists:
f = open("pdb_id_Q_and_len.txt", "w+")
f_percentages = open("pdb_id_residue_percentages.txt", "w+")
f_percentages.write("PDB_ID\tHydrophobic(%)\tPolar(%)\tPos(%)\tNeg(%)\tSpecial(%)\n")

#1) Reads a line (starting with the first) of the pdb_entry_type.txt file:
g = open(sys.argv[1], "r")
lines_of_pdb_entry_type_txt = g.readlines()
bandera_auxiliar = 0
for j in lines_of_pdb_entry_type_txt:
    try:   
        if (j == "6ri3"):
            bandera_auxiliar = 1
        if ((bandera_auxiliar == 1) & (j != "6ri3")):    
            # Print the current line on the terminal:
            print(j)
            #2) If the line is that corresponding to a protein (not to a nucleic
            # acid or a protein-nucleic acid complex) obtained by diffraction
            # then the pdb id contained in that line is downloaded from 
            # the RCSB database:
            words_of_j = j.split()
            if ((words_of_j[1] == "prot") & (words_of_j[2] == "diffraction")):
                # Extract pdb ID from the line:
                pdb_id_of_j = words_of_j[0]
                # Download that pdb from the RCSB database:
                try:
                    download_pdb(pdb_id_of_j)

                    #3) The obtain_Q function from
                    # calculate_Q.py is invoked with a try()
                    # statement, passing as argument the route of the downloaded pdb file:
                    try:
                        Q_and_len_to_write = obtain_Q_and_len(pdb_id_of_j + ".pdb")
                        #4) If the obtain_Q function works, then the pdb id
                        # and its corresponding Q value and number of residues 
                        # are appended to pdb_id_and_obtained_Q.txt .
                        f.write(str(pdb_id_of_j) + "\t" + str(Q_and_len_to_write[0]) + "\t" + str(Q_and_len_to_write[1]) + "\n")
                    
                        # Extract sequence from PDB file using BioPython
                        sequence = extract_sequence_from_pdb(pdb_id_of_j + ".pdb")
                        
                        # Calculate residue percentages and write to the second output file
                        percentages = calculate_residue_percentages(sequence)
                        f_percentages.write(f"{pdb_id_of_j}\t{percentages['hydrophobic']:.2f}\t{percentages['polar']:.2f}\t{percentages['pos']:.2f}\t{percentages['neg']:.2f}\t{percentages['special']:.2f}\n")
                    
                    except:
                        pass
                    
                    # Remove temporary files
                    remove_pdb_file(pdb_id_of_j)
                    remove_pdb_file(pdb_id_of_j, "-no_hetatoms")
                
                except:
                    pass         

    except:
        pass

f.close()
f_percentages.close()
g.close()




