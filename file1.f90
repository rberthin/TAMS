PROGRAM TAMS

        USE MD_STUFF
        USE RDF_3D
        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: n_atoms, n_step

        INTEGER :: io, i, j
        INTEGER :: xyz_unit = 10
        DOUBLE PRECISION, EXTERNAL :: RDF3D
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

!https://physics.emory.edu/faculty/weeks/idl/gofr2.html
!******************************************************************************!
!DOUBLE PRECISION FUNCTION RDF3D(dr, rmin, rmax, name_p1, name_p2, num_p1, num_p2, Posref, Posobs) !RESULT (bound) 
!
!        USE MD_STUFF
!
!        IMPLICIT NONE
!        INTEGER :: i, j, k
!        INTEGER, INTENT(IN) :: dr, num_p1, num_p2
!        REAL :: rmax, rmin, d_ref_obs
!        DOUBLE PRECISION, DIMENSION(11) :: bound
!        CHARACTER(LEN = 10), INTENT(IN) :: name_p1, name_p2
!
!        !** Definition of bound **!
!        !-------------------------!
!        bound(1) = rmin
!        DO k = 1, dr
!                bound(j+1) = rmin + j*((rmax - rmin)/dr)
!        END DO
!        !-------------------------!
!
!        DO i = 1, num_p1 ! Loop over all the ref. atoms
!                DO j = 1, num_p2 ! Loop over the observed atoms
!
!                        ! call pbc to change X, Y and Z positions
!                        X = dist_PBC(Posobs(j, 1), Posref(i, 1), boxx)
!                        Y = dist_PBC(Posobs(j, 2), Posref(i, 2), boxy)
!                        Z = dist_PBC(Posobs(j, 3), Posref(i, 3), boxz)
!
!                        ! compute distance ref - obs (with pbc)
!                        d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )
!                
!                        !IF (d_ref_obs <= rmax) THEN
!                        !placer dans le bon bin
!
!                END DO
!        END DO
!END FUNCTION RDF3D         
