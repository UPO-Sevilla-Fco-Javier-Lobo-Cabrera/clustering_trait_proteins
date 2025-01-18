#This function receives as an argument the route of a pdb file and
#returns the Q value of that pdb structure.

def obtain_Q_and_len(route_of_pdb_file):
    '''Receives the route of a pdb file and returns the Q value and number of residues''' 
    import sys
    import os
    from types_and_coordinates import obtain_types_and_coordinates
    from score_for_one_type import score_type_aa 


    #Generate the list containing the type and coordinates of each residue (lista_tipos_y_coordenadas):
    lista_tipos_y_coordenadas = obtain_types_and_coordinates(route_of_pdb_file)
    
    #Initialize the result (Q):
    Q = 0.0

    #Auxiliary variable to store the contribution to Q of each type of amino acyd:
    aux = 0.0

    #For each type of residue (not including the "special" type):
    for i in ["neg", "pos", "polar", "hydrophobic"]: 
        aux = score_type_aa(lista_tipos_y_coordenadas, i)
        Q = Q + aux
    
    
    number_residues = len(lista_tipos_y_coordenadas)
    
    #Return the result:
    return(Q, number_residues)






