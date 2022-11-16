# TAMS
Tools to Analyze Molecular Simulations

gfortran -c covalent_radius.f90 lexical_sort.f90 md_stuff.f90 rdf3d.f90 molecule_recognition.f90 tams_function.f90 tams.f90

gfortran tams.f90 covalent_radius.o lexical_sort.o md_stuff.o rdf3d.o molecule_recognition.o tams_function.o
