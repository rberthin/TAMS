MODULE RDF_3D
        IMPLICIT NONE
        INTEGER :: rdf3_dr, num_ref_atom, num_obs_atom
        REAL :: r_max, r_min
        CONTAINS
!******************************************************************************!
        SUBROUTINE get_infos_rdf3d(ATOM_NAME, n_atoms)
                IMPLICIT NONE
                INTEGER :: c, i
                INTEGER, INTENT(IN) :: n_atoms
                INTEGER, ALLOCATABLE, DIMENSION(:) :: RES_TEST, REF_LOC, OBS_LOC
                CHARACTER(LEN = 10) :: rdf3_name1, rdf3_name2
                CHARACTER(LEN = 10), DIMENSION(n_atoms) :: ATOM_NAME

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
        DOUBLE PRECISION FUNCTION RDF3D(name_p1, name_p2, Posref, Posobs) !RESULT (bound) 
        
                USE MD_STUFF
                !https://physics.emory.edu/faculty/weeks/idl/gofr2.html
        
                IMPLICIT NONE
                INTEGER :: i, j, k
                DOUBLE PRECISION, DIMENSION(11) :: bound
                CHARACTER(LEN = 10), INTENT(IN) :: name_p1, name_p2
                DOUBLE PRECISION, DIMENSION(num_ref_atom) :: Posref
                DOUBLE PRECISION, DIMENSION(num_obs_atom) :: Posobs
                !** Definition of bound **!
                !-------------------------!
                bound(1) = r_min
                DO k = 1, rdf3_dr
                        bound(j+1) = r_min + j*((r_max - r_min)/rdf3_dr)
                END DO
                !-------------------------!
        
!                DO i = 1, num_p1 ! Loop over all the ref. atoms
!                        DO j = 1, num_p2 ! Loop over the observed atoms
!        
!                                ! call pbc to change X, Y and Z positions
!                                X = dist_PBC(Posobs(j, 1), Posref(i, 1), boxx)
!                                Y = dist_PBC(Posobs(j, 2), Posref(i, 2), boxy)
!                                Z = dist_PBC(Posobs(j, 3), Posref(i, 3), boxz)
!        
!                                ! compute distance ref - obs (with pbc)
!                                d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )
!                        
!                                !IF (d_ref_obs <= rmax) THEN
!                                !placer dans le bon bin
!        
!                        END DO
!                END DO
        END FUNCTION RDF3D  
!******************************************************************************!
        
END MODULE RDF_3D
