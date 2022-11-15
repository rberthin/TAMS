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
        INTEGER :: io
        INTEGER :: i, j, s!, function_choice
        INTEGER :: rdf3_dr
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: res_rdf3d
        !--------------------------!
        
        CALL get_traj_name()

        ! From the number of line, guess number or steps
        !** To be improved (maybe check by size of the file)**!        
        CALL get_step_and_natoms()
        
        CALL SORT_ATOM_NAME()

        ! Getting the box parameters infos
        CALL get_box_parameters() !boxx, boxy, boxz)

        !WRITE(*,*) 'Which function do you want to compute ?'

        !WRITE(*,*) 'Please, choose a number between 1 and 1'
        !WRITE(*,*) '-------- Available functions ----------'
        !WRITE(*,*) '---------------------------------------'
        !WRITE(*,*) '1- Radial distribution function (RDF)'
        !WRITE(*,*) '---------------------------------------'
        !READ(5,*) function_choice
        call choose_function()
        !IF (function_choice == 1) THEN

                !** READ 1 STEP OF THE XYZ TRAJ **!
                !---------------------------------!
         !       CALL read_xyz() 

                !** RDF 3D STUFF **!
                !------------------!
          !      CALL get_infos_rdf()!ATOM_NAME, n_atoms, rdf3_dr)
           !     REWIND(10)

!                ALLOCATE(res_rdf3d(rdf3_dr, 2))
!                res_rdf3d = RDF3D(int(n_steps), rdf3_dr, xyz_unit, ATOM_NAME, &
!                                  POS, n_atoms, boxx, boxy, boxz)
!
!                OPEN(UNIT = 11, FILE = 'rdf.dat')
!                DO i = 1, rdf3_dr
!                        WRITE(11,'(F10.5, F12.8)') res_rdf3d(i,1), res_rdf3d(i,2)
!                END DO
!        END IF            
END PROGRAM TAMS

