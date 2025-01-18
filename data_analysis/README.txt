Note: Linux is required.

Execute:

a) python main.py pdb_entry.txt

being pdb_entry.txt the path to the pdb_entry_type.txt, which
is a file that contains in each line a pdb entry id, the type of biomolecule and
the technique used to obtain the structure, i.e:
100d	nuc	diffraction
101d	nuc	diffraction
101m	prot	diffraction
102d	nuc	diffraction
102l	prot	diffraction
102m	prot	diffraction
103d	nuc	NMR
103l	prot	diffraction
...

This execution will generate i)the file pdb_id_Q_and_len.txt (which contains the PDB id, Q value and number of residues)
and ii)the file pdb_id_residue_percentages.txt (which contains the PDB id, and residue composition per type).


b) python main_compactness.py pdb_entry.txt

being pdb_entry.txt the path to the pdb_entry_type.txt, which
is a file that contains in each line a pdb entry id, the type of biomolecule and
the technique used to obtain the structure, i.e:
100d	nuc	diffraction
101d	nuc	diffraction
101m	prot	diffraction
102d	nuc	diffraction
102l	prot	diffraction
102m	prot	diffraction
103d	nuc	NMR
103l	prot	diffraction
...

This execution will generate the file pdb_id_and_compactness.txt, which contains the PDB id and associated residue density in a sphere of minimum volume.


c) The analysis_of_results.Rmd file generates the plots and other useful information in this work from the obtained results. Specifically, it makes its calculations based on the results stored in the directory already_computed_results. So, if you would like analysis_of_results.Rmd to make calculations based on your executions you need first to locate your generated files (pdb_id_Q_and_len.txt, pdb_id_residue_percentages.txt and pdb_id_and_compactness.txt) within already_computed_results (replace the current files in that directory).


