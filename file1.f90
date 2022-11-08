PROGRAM TAMS

        USE MD_STUFF
        USE RDF_3D
        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: n_atoms, n_step

        INTEGER :: io, i, j, s, function_choice
        INTEGER :: xyz_unit = 10, n_line, rdf3_dr
        REAL :: boxx, boxy, boxz, n_steps
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: res_rdf3d
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS
        
        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        CHARACTER(LEN = 30) :: xyz_filename, input_name = 'tams_input.dat'
        !--------------------------!
        ! TO DO
        ! ADD POSSIBILITY TO OPEN TAMS_INPUT.DAT

        WRITE(*,*) 'Name of the XYZ trajectorie ?'
        READ(*,*) xyz_filename

        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        READ(10,*) n_atoms
        REWIND(10)
        WRITE(*,*) n_atoms, 'atoms found'
        n_line = num_lines_file(xyz_unit)
        n_steps = (n_line - 2) / n_atoms
        WRITE(*,*) int(n_steps), 'steps found'

        ALLOCATE(POS(n_atoms, 3))
        ALLOCATE(ATOM_NAME(n_atoms))

        ! Getting the box parameters infos
        CALL get_box_parameters(boxx, boxy, boxz)

        WRITE(*,*) 'Which function do you want to compute ?'
        WRITE(*,*) 'Please, choose a number between 1 and 1'
        WRITE(*,*) '-------- Available functions ----------'
        WRITE(*,*) '---------------------------------------'
        WRITE(*,*) '1- Radial distribution function (RDF)'
        WRITE(*,*) '---------------------------------------'
        READ(5,*) function_choice
        IF (function_choice == 1) THEN

                !** READ 1 STEP OF THE XYZ TRAJ **!
                !---------------------------------!
                CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)

                !** RDF 3D STUFF **!
                !------------------!
                CALL get_infos_rdf(ATOM_NAME, n_atoms, rdf3_dr)
                REWIND(10)

                ALLOCATE(res_rdf3d(rdf3_dr, 2))
                res_rdf3d = RDF3D(int(n_steps), rdf3_dr, xyz_unit, ATOM_NAME, &
                                  POS, n_atoms, boxx, boxy, boxz)

                OPEN(UNIT = 11, FILE = 'rdf.dat')
                DO i = 1, rdf3_dr
                        WRITE(11,'(F10.5, F12.8)') res_rdf3d(i,1), res_rdf3d(i,2)
                END DO
        END IF            
END PROGRAM TAMS

