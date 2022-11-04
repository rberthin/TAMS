MODULE RDF_3D
        IMPLICIT NONE
        INTEGER :: rdf3_dr, num_ref_atom, num_obs_atom
        REAL :: r_max, r_min
        INTEGER, ALLOCATABLE, DIMENSION(:) :: REF_LOC
        INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_LOC
        CONTAINS

!******************************************************************************!
        SUBROUTINE get_infos_rdf3d(ATOM_NAME, n_atoms)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: n_atoms
                CHARACTER(LEN = 10), DIMENSION(n_atoms), INTENT(IN) :: ATOM_NAME
               
                INTEGER :: c, i
                CHARACTER(LEN = 10) :: rdf3_name1, rdf3_name2

                WRITE(*,*) 'COMPUTATION OF 3D RDF'
                WRITE(*,*) '---------------------'
                WRITE(*,*) 'Name of the reference atom?'
                READ(*,*) rdf3_name1
                IF (ANY(ATOM_NAME == rdf3_name1)) THEN
                        num_ref_atom = COUNT(ATOM_NAME == rdf3_name1)
                        ALLOCATE(REF_LOC(num_ref_atom))
                        c = 1
                        DO i = 1, n_atoms
                                IF (ATOM_NAME(i) == rdf3_name1) THEN
                                        REF_LOC(c) = i
                                        c = c + 1
                                END IF
                        END DO
                        
                        WRITE(*,*) 'Name of the observed atom?'
                        READ(*,*) rdf3_name2
                        IF (ANY(ATOM_NAME == rdf3_name2)) THEN
                                num_obs_atom = COUNT(ATOM_NAME == rdf3_name2)
                                ALLOCATE(OBS_LOC(num_obs_atom))
                                c = 1
                                DO i = 1, n_atoms
                                        IF (ATOM_NAME(i) == rdf3_name2) THEN
                                                OBS_LOC(c) = i
                                                c = c + 1
                                        END IF
                                END DO
                                WRITE(*,*) 'Number of bins?'
                                READ(*,*) rdf3_dr
                        ELSE ; WRITE(*,*) 'I do not know this name, try again! :)'
                        END IF
                ELSE ; WRITE(*,*) 'I do not know this name, try again! :)'
                END IF

                WRITE(*,*) ' Enter the value of r_min :'
                READ(*,*) r_min
                WRITE(*,*) ' Enter the value of r_max :'
                READ(*,*) r_max

        END SUBROUTINE get_infos_rdf3d

!******************************************************************************!
        FUNCTION RDF3D(init_pos, n_atoms, boxx, boxy, boxz) RESULT(in_bound) 
        
                USE MD_STUFF
                !https://physics.emory.edu/faculty/weeks/idl/gofr2.html
        
                IMPLICIT NONE
                INTEGER :: i, j, k
                INTEGER, INTENT(IN) :: n_atoms
                REAL :: boxx, boxy, boxz
                DOUBLE PRECISION :: X, Y, Z, d_ref_obs
                DOUBLE PRECISION, DIMENSION(rdf3_dr+1) :: bound
                INTEGER, DIMENSION(rdf3_dr+1) :: in_bound
                DOUBLE PRECISION, DIMENSION(n_atoms,3) :: init_pos

                !** Definition of each bin **!
                !----------------------------!
                bound(1) = r_min
                DO k = 1, rdf3_dr
                        bound(k+1) = r_min + k*((r_max - r_min)/rdf3_dr)
                END DO
                !----------------------------!
                in_bound = 0
                DO i = 1, num_ref_atom ! Loop over all the ref. atoms
                        DO j = 1, num_obs_atom ! Loop over the observed atoms
                                ! call pbc to change X, Y and Z positions
                                X = dist_PBC(init_pos(OBS_LOC(j), 1), init_pos(REF_LOC(i), 1), boxx)
                                Y = dist_PBC(init_pos(OBS_LOC(j), 2), init_pos(REF_LOC(i), 2), boxy)
                                Z = dist_PBC(init_pos(OBS_LOC(j), 3), init_pos(REF_LOC(i), 3), boxz)

                                ! compute distance ref - obs (with pbc)
                                d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )

                                IF (d_ref_obs <= r_max) THEN
                                        DO k = 1, rdf3_dr+1 ! Loop over the bin
                                                IF (d_ref_obs <= bound(i)) THEN
                                                        in_bound(k) = in_bound(k) + 1
                                                END IF
                                        END DO
                                END IF
!        
                        END DO
                END DO
        END FUNCTION RDF3D  
!******************************************************************************!
        
END MODULE RDF_3D
