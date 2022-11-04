MODULE RDF_3D
        IMPLICIT NONE
        INTEGER :: num_ref_atom, num_obs_atom
        REAL :: r_max, r_min
        INTEGER, ALLOCATABLE, DIMENSION(:) :: REF_LOC
        INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_LOC
        CONTAINS

!******************************************************************************!
        SUBROUTINE get_infos_rdf3d(ATOM_NAME, n_atoms, rdf3_dr)
                IMPLICIT NONE

                INTEGER, INTENT(IN) :: n_atoms
                INTEGER :: rdf3_dr
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
        FUNCTION RDF3D(n_steps, rdf3_dr, xyz_unit, ATOM_NAME, init_pos, n_atoms, boxx, boxy, boxz) RESULT(res) 
        
                USE MD_STUFF
                !https://physics.emory.edu/faculty/weeks/idl/gofr2.html
        
                IMPLICIT NONE
                INTEGER :: i, j, k, s, rdf3_dr, n_steps, xyz_unit
                INTEGER, INTENT(IN) :: n_atoms
                REAL :: boxx, boxy, boxz
                LOGICAL :: find
                DOUBLE PRECISION :: X, Y, Z, d_ref_obs
                DOUBLE PRECISION, DIMENSION(rdf3_dr) :: bound, in_bound
                DOUBLE PRECISION, DIMENSION(rdf3_dr, 2) :: res
                DOUBLE PRECISION, DIMENSION(n_atoms,3) :: init_pos
                CHARACTER(LEN = 10), DIMENSION(n_atoms) :: ATOM_NAME

                !** Definition of each bin **!
                !----------------------------!
                bound(1) = r_min
                DO k = 1, rdf3_dr
                        bound(k+1) = r_min + k*((r_max - r_min)/rdf3_dr)
                END DO
                write(*,*) bound
                !----------------------------!
                DO s = 1, n_steps
                        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, init_pos)
                        DO i = 1, num_ref_atom ! Loop over all the ref. atoms
                                DO j = 1, num_obs_atom ! Loop over the observed atoms
                                        ! call pbc to change X, Y and Z positions
                                        X = dist_PBC(init_pos(OBS_LOC(j), 1), init_pos(REF_LOC(i), 1), boxx)
                                        Y = dist_PBC(init_pos(OBS_LOC(j), 2), init_pos(REF_LOC(i), 2), boxy)
                                        Z = dist_PBC(init_pos(OBS_LOC(j), 3), init_pos(REF_LOC(i), 3), boxz)

                                        ! compute distance ref - obs (with pbc)
                                        d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )
                                        
                                        IF (d_ref_obs <= r_max) THEN
                                                find = .TRUE.
                                                DO k = 1, rdf3_dr ! Loop over the bin
                                                        IF (d_ref_obs <= bound(k) .AND. find) THEN
                                                                in_bound(k) = in_bound(k) + 1
                                                                find = .FALSE.   
                                                        END IF
                                                END DO
                                        END IF
!        
                                END DO
                        END DO
                END DO
                DO k = 1, rdf3_dr
                        in_bound(k) = in_bound(k) / n_steps
                        in_bound(k) = in_bound(k) / num_ref_atom
                        if (bound(k)/= 0) then
                                in_bound(k) = in_bound(k) / (4 * 3.14 * bound(k)**2 * ((r_max - r_min)/rdf3_dr) )
                        end if
                        in_bound(k) = in_bound(k) / (num_obs_atom/(boxx*boxy*boxz))
                END DO
                DO k = 1, rdf3_dr
                        res(k, 1) = bound(k)
                        res(k, 2) = in_bound(k)
                END DO
        END FUNCTION RDF3D  
!******************************************************************************!
        
END MODULE RDF_3D
