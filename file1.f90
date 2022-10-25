PROGRAM TAMS

        IMPLICIT NONE
        !** VARIABLE DECLARATION **!
        !--------------------------!
        INTEGER :: io
        INTEGER :: n_atoms, n_step, xyz_unit = 10
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: POS
        CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_NAME
        CHARACTER(LEN = 30) xyz_filename
        !--------------------------!
        WRITE(*,*) 'Number of atoms ?'
        READ(*,*) n_atoms
        WRITE(*,*) 'Name of XYZ trajectorie ?'
        READ(*,*) xyz_filename

        ALLOCATE(POS(n_atoms, 3))
        ALLOCATE(ATOM_NAME(n_atoms))
        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        write(*,*) POS(1,1)
        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)
        write(*,*) POS(1,1)
END PROGRAM TAMS

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
