
import os
import sys
from obtain_Q1_imaginary_protein import obtain_Q1_and_len


#0)Creates a file named "pdb_id_and_obtained_Q1.txt", which
#will contain the Q value and number of residues of each
#imaginary protein:
f = open("Q1_and_number_of_residues_imaginary_protein.txt", "w+")

# Number of imaginary proteins to create:
num_imaginary_proteins = 10000

for i in range(0, num_imaginary_proteins):
    print(i)
    Q1_and_len_to_write = obtain_Q1_and_len()
    f.write(str(Q1_and_len_to_write[0]) + "\t" + str(Q1_and_len_to_write[1]) + "\n")
    
f.close()






