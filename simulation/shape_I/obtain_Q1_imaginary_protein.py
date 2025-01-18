#This function receives as an argument the route of a pdb file and
#returns the Q value of that pdb structure.

def obtain_Q1_and_len():
    '''Extracts Q score and number of residues of imaginary protein''' 
    import sys
    import os
    from generate_artificial_proteins_random_uniform_distribution import obtain_types_and_coordinates_imaginary_protein
    from score_for_one_type import score_type_aa 


    #Generate the list containing the type and coordinates of each residue (lista_tipos_y_coordenadas):
    lista_tipos_y_coordenadas = obtain_types_and_coordinates_imaginary_protein()
    
    
    #Initialize the result (Q):
    Q1 = 0.0

    #Auxiliary variable to store the contribution to Q of each type of amino acid:
    aux = 0.0

    #For each type of residue (not including the "special" type):
    for i in ["neg", "pos", "polar", "hydrophobic"]: 
        aux = score_type_aa(lista_tipos_y_coordenadas, i)
        Q1 = Q1 + aux
    
    number_residues = len(lista_tipos_y_coordenadas)
    
    #Return the result (Q and number of residues of the protein):
    return(Q1, number_residues)





