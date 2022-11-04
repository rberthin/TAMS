PROGRAM TAMS

        USE MD_STUFF
        USE RDF_3D
        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: n_atoms, n_step

        INTEGER :: io, i, j, s
        INTEGER :: xyz_unit = 10, n_line
        REAL :: boxx, boxy, boxz, n_steps
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: boundy
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS
        
        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        CHARACTER(LEN = 30) :: xyz_filename
        !--------------------------!
        WRITE(*,*) 'Number of atoms ?'
        READ(*,*) n_atoms
        WRITE(*,*) 'Name of XYZ trajectorie ?'
        READ(*,*) xyz_filename

        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        n_line = num_lines_file(xyz_unit)
        n_steps = (n_line - 2) / n_atoms
        WRITE(*,*) int(n_steps), 'steps found'

        ALLOCATE(POS(n_atoms, 3))
        ALLOCATE(ATOM_NAME(n_atoms))

        ! Getting the box parameters
        call get_box_parameters(boxx, boxy, boxz)

        !** READ 1 STEP OF THE XYZ TRAJ **!
        !---------------------------------!
        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)

        !** RDF 3D STUFF **!
        !------------------!
        CALL get_infos_rdf3d(ATOM_NAME, n_atoms)
        boundy = RDF3D(POS, n_atoms, boxx, boxy, boxz)
        DO s = 1, int(n_steps-1)
                WRITE(*,*) s
                CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)
                write(*,*) 'done read'
                boundy = RDF3D(POS, n_atoms, boxx, boxy, boxz)
        END DO
            
END PROGRAM TAMS

