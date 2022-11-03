MODULE MD_STUFF
        IMPLICIT NONE
        CONTAINS
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
        REAL(8) FUNCTION dist_PBC(coord_a, coord_b, cell_size) RESULT(res)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: coord_a
                REAL(8), INTENT(IN) :: coord_b
                REAL(8), INTENT(IN) :: cell_size
 
                INTEGER :: i
                REAL(8) :: dxcf, halfboxxrec
 
                halfboxxrec = 2.0/cell_size
 
                dxcf = coord_b - coord_a
 
                ! minimal distance convenction
                dxcf = dxcf - cell_size*int(dxcf*halfboxxrec)

                res=dxcf
        END FUNCTION dist_PBC

!******************************************************************************!
        
END MODULE MD_STUFF