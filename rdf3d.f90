MODULE RDF_3D
        IMPLICIT NONE
        CONTAINS
!******************************************************************************!
        SUBROUTINE get_infos_rdf3d(ATOM_NAME, n_atoms)
                IMPLICIT NONE
                INTEGER :: c, i
                INTEGER, INTENT(IN) :: n_atoms
                INTEGER :: rdf3_dr, num_ref_atom, num_obs_atom
                REAL :: r_max, r_min
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
        
END MODULE RDF_3D
