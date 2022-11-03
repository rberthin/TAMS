PROGRAM TAMS

        USE MD_STUFF
        USE RDF_3D
        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: n_atoms, n_step

        INTEGER :: io, i, j
        INTEGER :: xyz_unit = 10
        REAL :: boxx, boxy, boxz
        DOUBLE PRECISION, DIMENSION(11) :: boundy
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS
        
        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        CHARACTER(LEN = 30) :: xyz_filename
        !--------------------------!
        WRITE(*,*) 'Number of atoms ?'
        READ(*,*) n_atoms
        WRITE(*,*) 'Name of XYZ trajectorie ?'
        READ(*,*) xyz_filename

        ALLOCATE(POS(n_atoms, 3))
        ALLOCATE(ATOM_NAME(n_atoms))

        ! Getting the box parameters
        call get_box_parameters(boxx, boxy, boxz)


        !** READ 1 STEP OF THE XYZ TRAJ **!
        !---------------------------------!
        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)

        !** RDF 3D STUFF **!
        !------------------!
        CALL get_infos_rdf3d(ATOM_NAME, n_atoms)
        
END PROGRAM TAMS

!END FUNCTION RDF3D         
