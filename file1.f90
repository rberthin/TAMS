PROGRAM TAMS

        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: io, i, res_test
        INTEGER :: n_atoms, n_step, xyz_unit = 10, rdf3_dr, num_ref_atom
        REAL :: boxx, boxy, boxz
        REAL :: r_max, r_min
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS
        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        CHARACTER(LEN = 30) :: xyz_filename
        CHARACTER(LEN = 10) :: rdf3_name1, rdf3_name2
        !--------------------------!
        WRITE(*,*) 'Number of atoms ?'
        READ(*,*) n_atoms
        WRITE(*,*) 'Name of XYZ trajectorie ?'
        READ(*,*) xyz_filename

        ALLOCATE(POS(n_atoms, 3))
        ALLOCATE(ATOM_NAME(n_atoms))

        ! Getting the box parameters
        call get_box_parameters(boxx, boxy, boxz)

        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)

        !** READ XYZ TRAJ **!
        !-------------------!
        DO i = 1, 4
                CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)
        END DO

        !** RDF 3D STUFF **!
        !------------------!
        WRITE(*,*) 'COMPUTATION OF 3D RDF'
        WRITE(*,*) '---------------------'
        WRITE(*,*) 'Name of the reference atom?'
        READ(*,*) rdf3_name1
        IF (ANY(ATOM_NAME == rdf3_name1)) THEN
                WRITE(*,*) 'Name of the observed atom?'
                READ(*,*) rdf3_name2
                IF (ANY(ATOM_NAME == rdf3_name1)) THEN
                        WRITE(*,*) 'Number of bins?'
                        READ(*,*) rdf3_dr
                ELSE ; WRITE(*,*) 'Do not know this name, try again! :)'
                END IF
        ELSE ; WRITE(*,*) 'Do not know this name, try again! :)'
        END IF

        WRITE(*,*) ' Enter the value of r_min :'
        READ(*,*) r_min
        WRITE(*,*) ' Enter the value of r_max :'
        READ(*,*) r_max

        num_ref_atom = COUNT(ATOM_NAME == rdf3_name1)

END PROGRAM TAMS

!******************************************************************************!
SUBROUTINE read_xyz(file_unit, n_atoms, name_array, pos_array)
        IMPLICIT NONE
        INTEGER :: i, n_atoms, file_unit
        DOUBLE PRECISION, DIMENSION(n_atoms, 3), INTENT(INOUT) :: pos_array
        CHARACTER(LEN = 10), DIMENSION(n_atoms), INTENT(INOUT) :: name_array

        READ(file_unit, *)
        READ(file_unit, *)

        DO i = 1, n_atoms
                READ(file_unit, *) name_array(i), pos_array(i, 1), &
                                   pos_array(i, 2), pos_array(i, 3)
        END DO
END SUBROUTINE read_xyz

!******************************************************************************!
SUBROUTINE get_box_parameters(boxx, boxy, boxz)
        IMPLICIT NONE
        LOGICAL :: step_while = .TRUE. 
        REAL :: boxx, boxy, boxz
        CHARACTER(LEN = 5) answer 

        DO WHILE(step_while)
                WRITE(*,*) 'Is the simulation box cubic (y/n)?'
                READ(*,*) answer
                IF (answer == 'y') THEN
                        WRITE(*,*) 'What is the box length?'
                        READ(*,*) boxx
                        boxy = boxx
                        boxz = boxx
                        step_while = .FALSE.
                ELSE IF (answer == 'n') THEN
                        WRITE(*,*) 'What is the box length along x?'
                        READ(*,*) boxx
                        WRITE(*,*) 'What is the box length along y?'
                        READ(*,*) boxy
                        WRITE(*,*) 'What is the box length along z?'
                        READ(*,*) boxz
                        step_while = .FALSE.
                ELSE
                        WRITE(*,*) 'Something went wrong, try again :)'
                END IF
        END DO

END SUBROUTINE get_box_parameters

!******************************************************************************!
!FUNCTION RDF3D(dr, rmin, rmax, name_p1, name_p2, num_particule)
!
!        IMPLICIT NONE
!
!END FUNCTION RDF3D         
