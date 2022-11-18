PROGRAM TAMS

        USE MD_STUFF
        USE MOLECULAR_RECOGNITION
        USE RDF_3D
        USE covalent_radius
        USE lexical_sort
        USE TAMS_FUNCTION
        
        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        !INTEGER :: i, j, s
        !--------------------------!
        
        CALL get_traj_name()

        ! From the number of line, guess number or steps
        !** To be improved (maybe check by size of the file)**!        
        CALL get_step_and_natoms()
        
        CALL SORT_ATOM_NAME()

        ! Getting the box parameters infos
        CALL get_box_parameters() 

        CALL choose_function()

        CLOSE(xyz_unit)
        CLOSE(tamsinput_unit)
END PROGRAM TAMS

