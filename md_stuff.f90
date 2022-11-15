MODULE MD_STUFF
        IMPLICIT NONE
        INTEGER :: xyz_unit = 10, tamsinput_unit = 20
        INTEGER :: n_atoms, n_steps
        REAL :: boxx, boxy, boxz
        CHARACTER(LEN = 30) :: xyz_filename, input_name = 'tams_input.dat'
        LOGICAL :: file_input

        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS

        CONTAINS
!******************************************************************************!
        SUBROUTINE get_traj_name()
                IMPLICIT NONE
                INTEGER :: io
                CHARACTER(LEN = 100) crap
                INQUIRE(file = input_name, exist = file_input)
                
                IF (file_input .EQV. .TRUE.) THEN
                        OPEN(unit = tamsinput_unit, file = input_name, status='old', iostat=io)
                        WRITE(*,*) 'I found an input file'
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) xyz_filename

                ELSE
                        OPEN(unit = tamsinput_unit, file = input_name)
                        WRITE(*,*) 'Name of the XYZ trajectorie ?'
                        WRITE(tamsinput_unit,*) 'Name of the XYZ trajectorie ?'
                        READ(*,*) xyz_filename
                        WRITE(tamsinput_unit,*) xyz_filename
                END IF

        END SUBROUTINE get_traj_name

!******************************************************************************!
        SUBROUTINE get_step_and_natoms()
                IMPLICIT NONE
                INTEGER :: io
                INTEGER :: n_line

                OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
                READ(xyz_unit,*) n_atoms
                REWIND(xyz_unit)
                WRITE(*,*) n_atoms, 'atoms found'
                n_line = num_lines_file(xyz_unit)
                n_steps = (n_line - 2) / n_atoms
                WRITE(*,*) int(n_steps), 'steps found'
                ALLOCATE(POS(n_atoms, 3))
                ALLOCATE(ATOM_NAME(n_atoms))
        END SUBROUTINE get_step_and_natoms

!******************************************************************************!
        SUBROUTINE read_xyz() !name_array, pos_array)
                IMPLICIT NONE
                INTEGER :: i

                READ(xyz_unit, *)
                READ(xyz_unit, *)

                DO i = 1, n_atoms
                        READ(xyz_unit, *) ATOM_NAME(i), POS(i, 1), &
                                           POS(i, 2), POS(i, 3)
                END DO
        END SUBROUTINE read_xyz

!******************************************************************************!
        SUBROUTINE get_box_parameters()
                IMPLICIT NONE
                LOGICAL :: step_while = .TRUE. ! , file_input
                REAL :: boxx, boxy, boxz
                CHARACTER(LEN = 5) answer
                CHARACTER(LEN = 100) crap

                IF (file_input .EQV. .TRUE.) THEN
                        READ(tamsinput_unit,*) crap
                        READ(tamsinput_unit,*) answer
                        IF (answer == 'y') THEN
                                READ(tamsinput_unit,*) crap
                                READ(tamsinput_unit,*) boxx
                                boxy = boxx
                                boxz = boxx
                        ELSE 
                                READ(tamsinput_unit,*) crap
                                READ(tamsinput_unit,*) boxx
                                READ(tamsinput_unit,*) crap
                                READ(tamsinput_unit,*) boxy
                                READ(tamsinput_unit,*) crap
                                READ(tamsinput_unit,*) boxz
                        END IF

                ELSE
                        DO WHILE(step_while)
                        WRITE(*,*) 'Is the simulation box cubic (y/n)?'
                        WRITE(tamsinput_unit,*) 'Is the simulation box cubic (y/n)?'
                        READ(*,*) answer
                        WRITE(tamsinput_unit,*) answer
                        IF (answer == 'y') THEN
                                WRITE(*,*) 'What is the box length?'
                                WRITE(tamsinput_unit,*) 'What is the box length?'
                                READ(*,*) boxx
                                WRITE(tamsinput_unit,*) boxx
                                boxy = boxx
                                boxz = boxx
                                step_while = .FALSE.
                        ELSE IF (answer == 'n') THEN
                                WRITE(*,*) 'What is the box length along x?'
                                WRITE(tamsinput_unit,*) 'What is the box length along x?'
                                READ(*,*) boxx
                                WRITE(tamsinput_unit,*) boxx
                                WRITE(*,*) 'What is the box length along y?'
                                WRITE(tamsinput_unit,*) 'What is the box length along y?'
                                READ(*,*) boxy
                                WRITE(tamsinput_unit,*) boxy
                                WRITE(*,*) 'What is the box length along z?'
                                WRITE(tamsinput_unit,*) 'What is the box length along z?'
                                READ(*,*) boxz
                                WRITE(tamsinput_unit,*) boxz
                                step_while = .FALSE.
                        ELSE
                                WRITE(*,*) 'Something went wrong, try again :)'
                        END IF
                        END DO
                END IF
        
        END SUBROUTINE get_box_parameters
        
!******************************************************************************!
        REAL(8) FUNCTION dist_PBC(coord_a, coord_b, cell_size) RESULT(res)
                IMPLICIT NONE
                DOUBLE PRECISION, INTENT(IN) :: coord_a, coord_b
                REAL, INTENT(IN) :: cell_size
 
                INTEGER :: i
                REAL(8) :: dxcf, halfboxxrec
 
                halfboxxrec = 2.0/cell_size
 
                dxcf = coord_b - coord_a
 
                ! minimal distance convenction
                dxcf = dxcf - cell_size*int(dxcf*halfboxxrec)

                res = dxcf
        END FUNCTION dist_PBC

!******************************************************************************!
        INTEGER FUNCTION num_lines_file(num_file) RESULT(res)
        IMPLICIT NONE
        INTEGER, INTENT(IN) ::num_file
        INTEGER :: io
        CHARACTER*2 :: inputline
        INTEGER :: n_lin
        n_lin = 0
        DO
           READ(num_file,*,IOSTAT=io) inputline
           IF (io /= 0) THEN
              EXIT
           ELSE
              n_lin = n_lin + 1
           ENDIF
        ENDDO
        REWIND(num_file)
        res = n_lin
        END FUNCTION num_lines_file

!******************************************************************************!

END MODULE MD_STUFF
