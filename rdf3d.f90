MODULE RDF_3D
        USE MD_STUFF
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: pi  = 4*atan(1.0)
        INTEGER :: num_ref_atom, num_obs_atom, rdfdim_choice
        INTEGER :: rdf_dr
        REAL :: r_max, r_min

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: bound, in_bound
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: res
        !DOUBLE PRECISION, DIMENSION(n_atoms,3) :: init_pos

        INTEGER, ALLOCATABLE, DIMENSION(:) :: REF_LOC
        INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_LOC
        CHARACTER(LEN = 1) :: rdf2d_choice
        CONTAINS

!******************************************************************************!
        SUBROUTINE get_infos_rdf()
                IMPLICIT NONE

                INTEGER :: c, i
                LOGICAL :: select_dim = .FALSE., select_name1 = .FALSE.
                LOGICAL :: select_name2 = .FALSE., select_2d = .FALSE.
                CHARACTER(LEN = 10) :: rdf_name1, rdf_name2
                CHARACTER(LEN = 1) :: tempdev_choice
                CHARACTER(LEN = 100) :: crap
                IF (file_input .EQV. .TRUE.) THEN
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) rdfdim_choice
                        IF (rdfdim_choice == 3) THEN
                                WRITE(*,*) '* COMPUTATION OF 3D RDF'
                                select_dim = .TRUE.
                        ELSE IF (rdfdim_choice == 2) THEN
                                READ(tamsinput_unit,*) crap 
                                READ(tamsinput_unit,*) rdf2d_choice
                                WRITE(*,*) '* COMPUTATION OF 2D RDF'
                                select_dim = .TRUE.
                        END IF
                ELSE
                        DO WHILE (select_dim .EQV. .FALSE.)
                                WRITE(*,*) 'Do you want to perform 3D RDF (3) or 2D RDF (2)?'
                                WRITE(tamsinput_unit,'(A)') 'Do you want to perform 3D RDF (3) or 2D RDF (2)?' 
                                READ(*,*)  rdfdim_choice
                                WRITE(tamsinput_unit,'(I0)') rdfdim_choice
                                IF (rdfdim_choice == 3) THEN
                                        WRITE(*,*) 'COMPUTATION OF 3D RDF'
                                        select_dim = .TRUE.
                                        rdf2d_choice = 'n'
                                ELSE IF (rdfdim_choice == 2) THEN
                                        WRITE(*,*) 'In which dimension do you want to compute the 2D RDF (x, y, z)?'
                                        READ(*,*) rdf2d_choice
                                        WRITE(*,*) 'COMPUTATION OF 2D RDF'
                                        select_dim = .TRUE.
                                ELSE
                                        WRITE(*,*) 'I did not understand your answer (:'

                                END IF
                        END DO
                END IF

                WRITE(*,*) '------------------------'

                IF (file_input .EQV. .TRUE.) THEN
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) rdf_name1
                ELSE
                        DO WHILE (select_name1 .EQV. .FALSE.)
                                WRITE(*,*) 'Name of the reference atom?'
                                WRITE(tamsinput_unit,'(A)') 'Name of the reference atom?'
                                READ(*,*) rdf_name1
                                IF (ANY(ATOM_NAME == rdf_name1)) THEN
                                        WRITE(tamsinput_unit,'(A)') rdf_name1
                                        select_name1 = .TRUE.
                                ELSE ; WRITE(*,*) 'I do not know this name, try again! :)'
                                END IF
                        END DO
                END IF

                num_ref_atom = COUNT(ATOM_NAME == rdf_name1)
                ALLOCATE(REF_LOC(num_ref_atom))
                c = 1
                DO i = 1, n_atoms
                        IF (ATOM_NAME(i) == rdf_name1) THEN
                                REF_LOC(c) = i
                                c = c + 1
                        END IF
                END DO

                IF (file_input .EQV. .TRUE.) THEN
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) rdf_name2
                ELSE
                        DO WHILE (select_name2 .EQV. .FALSE.)
                                WRITE(*,*) 'Name of the observed atom?'
                                WRITE(tamsinput_unit,'(A)') 'Name of the observed atom?'
                                READ(*,*) rdf_name2
                                IF (ANY(ATOM_NAME == rdf_name2)) THEN
                                        WRITE(tamsinput_unit,'(A)') rdf_name2
                                        select_name2 = .TRUE.
                                ELSE ; WRITE(*,*) 'I do not know this name, try again! :)'
                                END IF
                        END DO
                END IF

                num_obs_atom = COUNT(ATOM_NAME == rdf_name2)
                ALLOCATE(OBS_LOC(num_obs_atom))
                c = 1
                DO i = 1, n_atoms
                        IF (ATOM_NAME(i) == rdf_name2) THEN
                                OBS_LOC(c) = i
                                c = c + 1
                        END IF
                END DO

                IF (file_input .EQV. .TRUE.) THEN
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) rdf_dr
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) r_min
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) r_max
                ELSE
                        WRITE(*,*) 'Number of bins?'
                        WRITE(tamsinput_unit,'(A)') 'Number of bins?'
                        READ(*,*) rdf_dr
                        WRITE(tamsinput_unit,'(I0)') rdf_dr
                        WRITE(*,*) 'Enter the value of r_min :'
                        WRITE(tamsinput_unit,'(A)') 'Enter the value of r_min :'
                        READ(*,*) r_min
                        WRITE(tamsinput_unit,'(F0.6)') r_min
                        WRITE(*,*) 'Enter the value of r_max :'
                        WRITE(tamsinput_unit,'(A)') 'Enter the value of r_max :'
                        READ(*,*) r_max
                        WRITE(tamsinput_unit,'(F0.6)') r_max
                END IF
                ALLOCATE(bound(rdf_dr+1))
                ALLOCATE(in_bound(rdf_dr))
                ALLOCATE(res(rdf_dr, 2))

                !WRITE(*,*) 'Do you want to save a temporal development of this rdf (y/n)?'
                !READ(*,*) tempdev_choice

                ! Continue temporal development part
                !IF (tempdev_choice == 'y') THEN
                !        WRITE(*,*) 'Do you want to save it for all the
               WRITE(*,*) '------------------------'
               WRITE(*,*) 'BEGIN 3D RDF COMPUTATION' 
        END SUBROUTINE get_infos_rdf

!******************************************************************************!
        SUBROUTINE RDF3D() 
        
                !https://physics.emory.edu/faculty/weeks/idl/gofr2.html
        
                IMPLICIT NONE
                INTEGER :: i, j, k, s
                LOGICAL :: find
                DOUBLE PRECISION :: X, Y, Z, d_ref_obs
                !write(*,*) rdf_dr
                !** Definition of each bin **!
                !----------------------------!
                bound(1) = r_min
                DO k = 1, rdf_dr
                        bound(k+1) = r_min + k*((r_max - r_min)/rdf_dr)
                END DO
                !----------------------------!
                DO s = 1, n_steps
                        write(*,*) s
                        CALL read_xyz() !xyz_unit, n_atoms, ATOM_NAME, init_pos)
                        DO i = 1, num_ref_atom ! Loop over all the ref. atoms
                                DO j = 1, num_obs_atom ! Loop over the observed atoms
                                        IF (rdf2d_choice == 'n') THEN
                                                ! call pbc to change X, Y and Z positions
                                                X = dist_PBC(POS(OBS_LOC(j), 1), POS(REF_LOC(i), 1), boxx)
                                                Y = dist_PBC(POS(OBS_LOC(j), 2), POS(REF_LOC(i), 2), boxy)
                                                Z = dist_PBC(POS(OBS_LOC(j), 3), POS(REF_LOC(i), 3), boxz)
                                        ELSE IF (rdf2d_choice == 'x') THEN
                                                X = POS(REF_LOC(i), 1) - POS(OBS_LOC(j), 1)
                                                Y = dist_PBC(POS(OBS_LOC(j), 2), POS(REF_LOC(i), 2), boxy)
                                                Z = dist_PBC(POS(OBS_LOC(j), 3), POS(REF_LOC(i), 3), boxz)
                                        ELSE IF (rdf2d_choice == 'y') THEN
                                                X = dist_PBC(POS(OBS_LOC(j), 1), POS(REF_LOC(i), 1), boxx)
                                                Y = POS(REF_LOC(i), 2) - POS(OBS_LOC(j), 2)
                                                Z = dist_PBC(POS(OBS_LOC(j), 3), POS(REF_LOC(i), 3), boxz)
                                        ELSE IF (rdf2d_choice == 'z') THEN
                                                ! call pbc to change X, Y and Z positions
                                                X = dist_PBC(POS(OBS_LOC(j), 1), POS(REF_LOC(i), 1), boxx)
                                                Y = dist_PBC(POS(OBS_LOC(j), 2), POS(REF_LOC(i), 2), boxy)
                                                Z = POS(REF_LOC(i), 3) - POS(OBS_LOC(j), 3)
                                        END IF
                                        ! compute distance ref - obs (with pbc)
                                        d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )
                                        
                                        IF (d_ref_obs <= r_max) THEN
                                                find = .TRUE.
                                                DO k = 1, rdf_dr ! Loop over the bin
                                                        IF (d_ref_obs <= bound(k) .AND. find) THEN
                                                                in_bound(k) = in_bound(k) + 1
                                                                find = .FALSE.   
                                                        END IF
                                                END DO
                                        END IF
        
                                END DO
                        END DO
                END DO
                OPEN(UNIT = 11, FILE = 'rdf.dat')
                DO k = 1, rdf_dr
                        in_bound(k) = in_bound(k) / n_steps
                        in_bound(k) = in_bound(k) / num_ref_atom
                        if (bound(k)/= 0) then
                                in_bound(k) = in_bound(k) / (4 * pi * bound(k)**2 * ((r_max - r_min)/rdf_dr) )
                        end if
                        in_bound(k) = in_bound(k) / (num_obs_atom/(boxx*boxy*boxz))

                        res(k, 1) = bound(k)
               !         write(*,*) bound(k)
                        res(k, 2) = in_bound(k)
                        WRITE(11,'(F10.5, F12.8)') res(k, 1), res(k, 2)
                END DO
        END SUBROUTINE RDF3D  
!******************************************************************************!
!        FUNCTION RDF2D(n_steps, rdf3_dr, xyz_unit, ATOM_NAME, init_pos, n_atoms, boxx, boxy, boxz) RESULT(res)
!
!                USE MD_STUFF
!                !https://physics.emory.edu/faculty/weeks/idl/gofr2.html
!
!                IMPLICIT NONE
!                INTEGER :: i, j, k, s, rdf3_dr, n_steps, xyz_unit
!                INTEGER, INTENT(IN) :: n_atoms
!                REAL :: boxx, boxy, boxz
!                LOGICAL :: find
!                DOUBLE PRECISION :: X, Y, Z, d_ref_obs
!                DOUBLE PRECISION, DIMENSION(rdf3_dr) :: bound, in_bound
!                DOUBLE PRECISION, DIMENSION(rdf3_dr, 2) :: res
!                DOUBLE PRECISION, DIMENSION(n_atoms,3) :: init_pos
!                CHARACTER(LEN = 10), DIMENSION(n_atoms) :: ATOM_NAME
!
!                !** Definition of each bin **!
!                !----------------------------!
!                bound(1) = r_min
!                DO k = 1, rdf3_dr
!                        bound(k+1) = r_min + k*((r_max - r_min)/rdf3_dr)
!                END DO
!                !----------------------------!
!                DO s = 1, n_steps
!                        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, init_pos)
!                        DO i = 1, num_ref_atom ! Loop over all the ref. atoms
!                                DO j = 1, num_obs_atom ! Loop over the observed atoms
!                                        ! call pbc to change X, Y and Z positions
!                                        X = dist_PBC(init_pos(OBS_LOC(j), 1), init_pos(REF_LOC(i), 1), boxx)
!                                        Y = dist_PBC(init_pos(OBS_LOC(j), 2), init_pos(REF_LOC(i), 2), boxy)
!                                        Z = dist_PBC(init_pos(OBS_LOC(j), 3), init_pos(REF_LOC(i), 3), boxz)
!
!                                        ! compute distance ref - obs (with pbc)
!                                        d_ref_obs = sqrt( X**2 + Y**2 + Z**2 )
!
!                                        IF (d_ref_obs <= r_max) THEN
!                                                find = .TRUE.
!                                                DO k = 1, rdf3_dr ! Loop over the bin
!                                                        IF (d_ref_obs <= bound(k) .AND. find) THEN
!                                                                in_bound(k) = in_bound(k) + 1
!                                                                find = .FALSE.
!                                                        END IF
!                                                END DO
!                                        END IF        
!                                END DO
!                       END DO
!                END DO
!                DO k = 1, rdf3_dr
!                        in_bound(k) = in_bound(k) / n_steps
!                        in_bound(k) = in_bound(k) / num_ref_atom
!                        if (bound(k)/= 0) then
!                                in_bound(k) = in_bound(k) / (4 * 3.14 * bound(k)**2 * ((r_max - r_min)/rdf3_dr) )
!                        end if
!                        in_bound(k) = in_bound(k) / (num_obs_atom/(boxx*boxy*boxz))
!                END DO
!                DO k = 1, rdf3_dr
!                        res(k, 1) = bound(k)
!                        res(k, 2) = in_bound(k)
!                END DO
!        END FUNCTION RDF2D
END MODULE RDF_3D
