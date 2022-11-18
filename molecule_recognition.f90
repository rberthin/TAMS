MODULE MOLECULAR_RECOGNITION

        USE COVALENT_RADIUS
        USE MD_STUFF
        USE LEXICAL_SORT

        INTEGER :: outfile_unit = 12
        CHARACTER(LEN = 50) :: outfile_tams = 'tams.out'
        CONTAINS
        SUBROUTINE SORT_ATOM_NAME() 
                IMPLICIT NONE
                INTEGER :: i, io, c, m, counter_same_name, loc(1) 
                CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_N, TEMP_ARRAY, ATOM_TYPE
                INTEGER, ALLOCATABLE, DIMENSION(:) :: ATOM_INDEX, COUNT_NAME, TMP_ARRAY
                INTEGER, DIMENSION(n_elem) :: COUNT_TYPE = 0
                CHARACTER(LEN = 2) :: answer
                CHARACTER(LEN = 100) :: crap

                CALL init_periodic_table()

                OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
                OPEN(unit = outfile_unit, file = outfile_tams, iostat=io)

                ALLOCATE(ATOM_INDEX(n_atoms))

                ! Read one step to get the atom name       
                CALL read_xyz() 
                
                ! Sort all the atom name 
                ! Input = ATOM_NAME
                ! ATOM_INDEX = well ordered index

                CALL sort(ATOM_NAME, ATOM_INDEX)
        
                ! Remove duplicates
                c = 1
                counter_same_name = 0
                ! ATOM_N = list of all the name without duplicates
                ALLOCATE(ATOM_N(c)) 
                ! COUNT_NAME = number of each type duplicates
                ALLOCATE(COUNT_NAME(c))
                
                ! Maybe re-do this loop bc shitty way
                DO i = 1, n_atoms-1 
                        IF (ATOM_NAME(i) /= ATOM_NAME(i+1)) THEN
        
                                ATOM_N(c) = ATOM_NAME(i)
                                ALLOCATE(TEMP_ARRAY(size(ATOM_N)+1))
                                TEMP_ARRAY(1:size(ATOM_N)) = ATOM_N
                                DEALLOCATE(ATOM_N)
                                ALLOCATE(ATOM_N(size(TEMP_ARRAY)))
                                ATOM_N(1:size(TEMP_ARRAY)) = TEMP_ARRAY
                                DEALLOCATE(TEMP_ARRAY)

                                !-------------------------------------
                                COUNT_NAME(c) = counter_same_name + 1
                                ALLOCATE(TMP_ARRAY(size(COUNT_NAME)+1))
                                TMP_ARRAY(1:size(COUNT_NAME)) = COUNT_NAME
                                DEALLOCATE(COUNT_NAME)
                                ALLOCATE(COUNT_NAME(size(TMP_ARRAY)))
                                COUNT_NAME(1:size(TMP_ARRAY)) = TMP_ARRAY
                                DEALLOCATE(TMP_ARRAY)
                                c = c + 1
                                counter_same_name = 0
                        ELSE
                                counter_same_name = counter_same_name + 1
                        END IF
                END DO

                ! Final value of array (last index value)
                ATOM_N(c) = ATOM_NAME(n_atoms)
                COUNT_NAME(c) = counter_same_name + 1
                
                ALLOCATE(ATOM_TYPE(size(ATOM_N)))

                WRITE(*,*) size(ATOM_N), 'atoms type found ...'
                WRITE(*,*) 'Starting assigning type to each atom ...'
                WRITE(*,*) '----------------------------------------'
                
                ! Give to each atom type a regular type depending
                ! on the name read in the XYZ file
                DO i = 1, size(ATOM_N)
                        IF (ATOM_N(i)(:1) == 'C') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('l', 'L')
                                        ATOM_TYPE(i) = 'Cl'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Cl (Chlorine) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cl (Chlorine) atom'
                                        COUNT_TYPE(17) = COUNT_TYPE(17) + COUNT_NAME(i)
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'Ca'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ca (Calcium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ca (Calcium) atom'
                                        COUNT_TYPE(20) = COUNT_TYPE(20) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Cr'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Cr (Chromium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cr (Chromium) atom'
                                        COUNT_TYPE(24) = COUNT_TYPE(24) + COUNT_NAME(i)
                                CASE('o', 'O')
                                        ATOM_TYPE(i) = 'Co'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Co (Cobalt) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Co (Cobalt) atom'
                                        COUNT_TYPE(27) = COUNT_TYPE(27) + COUNT_NAME(i)
                                CASE('u', 'U')
                                        ATOM_TYPE(i) = 'Cu'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Cu (Copper) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cu (Copper) atom'
                                        COUNT_TYPE(29) = COUNT_TYPE(29) + COUNT_NAME(i)
                                CASE('d', 'D')
                                        ATOM_TYPE(i) = 'Cd'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Cd (Cadmium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cd (Cadmium) atom'
                                        COUNT_TYPE(48) = COUNT_TYPE(48) + COUNT_NAME(i)
                                CASE('s', 'S')
                                        ATOM_TYPE(i) = 'Cs'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Cs (Cesium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cs (Cesium) atom'
                                        COUNT_TYPE(55) = COUNT_TYPE(55) + COUNT_NAME(i)
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Ce'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ce (Cerium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ce (Cerium) atom' 
                                        COUNT_TYPE(58) = COUNT_TYPE(58) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'C '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as C (Carbon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as C (Carbon) atom'
                                        COUNT_TYPE(6) = COUNT_TYPE(6) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'A') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('l', 'L')
                                        ATOM_TYPE(i) = 'Al'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Al (Aluminium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Al (Aluminium) atom'
                                        COUNT_TYPE(13) = COUNT_TYPE(13) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Ar'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ar (Argon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ar (Argon) atom'
                                        COUNT_TYPE(18) = COUNT_TYPE(18) + COUNT_NAME(i)
                                CASE('u', 'U')
                                        ATOM_TYPE(i) = 'Au'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Au (Gold) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Au (Gold) atom'
                                        COUNT_TYPE(79) = COUNT_TYPE(79) + COUNT_NAME(i)
                                CASE('s', 'S')
                                        ATOM_TYPE(i) = 'As'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as As (Arsenic) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as As (Arsenic) atom'
                                        COUNT_TYPE(33) = COUNT_TYPE(33) + COUNT_NAME(i)
                                CASE('g', 'G')
                                        ATOM_TYPE(i) = 'Ag'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ag (Silver) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ag (Silver) atom'
                                        COUNT_TYPE(47) = COUNT_TYPE(47) + COUNT_NAME(i)
                                CASE('t', 'T')
                                        ATOM_TYPE(i) = 'At'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as At (Astatine) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as At (Astatine) atom'
                                        COUNT_TYPE(85) = COUNT_TYPE(85) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = '  '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'B') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Be'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Be (Beryllium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Be (Beryllium) atom'
                                        COUNT_TYPE(4) = COUNT_TYPE(4) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Br'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Br (Br) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Br (Br) atom'
                                        COUNT_TYPE(35) = COUNT_TYPE(35) + COUNT_NAME(i)
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'Ba'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ba (Barium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ba (Barium) atom'
                                        COUNT_TYPE(56) = COUNT_TYPE(56) + COUNT_NAME(i)
                                CASE('i', 'I')
                                        ATOM_TYPE(i) = 'Bi'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Bi (Bismuth) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Bi (Bismuth) atom'
                                        COUNT_TYPE(83) = COUNT_TYPE(83) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'B '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as B (Boron) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as B (Boron) atom'
                                        COUNT_TYPE(5) = COUNT_TYPE(5) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'D') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('y', 'Y')
                                        ATOM_TYPE(i) = 'Dy'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Dy (Dysprosium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Dy (Dysprosium) atom'
                                        COUNT_TYPE(66) = COUNT_TYPE(66) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = ' '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'E') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('u', 'U')
                                        ATOM_TYPE(i) = 'Eu'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Eu (Europium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Eu (Europium) atom'
                                        COUNT_TYPE(63) = COUNT_TYPE(63) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Er'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Eu (Erbium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Eu (Erbium) atom'
                                        COUNT_TYPE(68) = COUNT_TYPE(68) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = ' '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'F') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Fe'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Fe (Iron) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Fe (Iron) atom'
                                        COUNT_TYPE(26) = COUNT_TYPE(26) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'F '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                                        COUNT_TYPE(9) = COUNT_TYPE(9) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'G') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'Ga'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ga (Gallium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                                        COUNT_TYPE(31) = COUNT_TYPE(31) + COUNT_NAME(i)
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Ge'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ge (Germanium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ge (Germanium) atom'
                                        COUNT_TYPE(32) = COUNT_TYPE(32) + COUNT_NAME(i)
                                CASE('d', 'D') 
                                        ATOM_TYPE(i) = 'Gd'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Gd (Gadolinium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Gd (Gadolinium) atom'
                                        COUNT_TYPE(64) = COUNT_TYPE(64) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = ' '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'H') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'He'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as He (Helium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as He (Helium) atom'
                                        COUNT_TYPE(2) = COUNT_TYPE(2) + COUNT_NAME(i)
                                CASE('o', 'O')
                                        ATOM_TYPE(i) = 'Ho'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ho (Holmium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ho (Holmium) atom'
                                        COUNT_TYPE(67) = COUNT_TYPE(67) + COUNT_NAME(i)
                                CASE('f', 'F')
                                        ATOM_TYPE(i) = 'Hf'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Hf (Hafnium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Hf (Hafnium) atom'
                                        COUNT_TYPE(72) = COUNT_TYPE(72) + COUNT_NAME(i)
                                CASE('g', 'G')
                                        ATOM_TYPE(i) = 'Hg'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Hg (Mercury) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Hg (Mercury) atom'
                                        COUNT_TYPE(80) = COUNT_TYPE(80) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'H '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as H (Hydrogen) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as H (Hydrogen) atom'
                                        COUNT_TYPE(1) = COUNT_TYPE(1) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'I') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('n', 'N')
                                        ATOM_TYPE(i) = 'In'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as In (Indium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as In (Indium) atom'
                                        COUNT_TYPE(49) = COUNT_TYPE(49) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Ir'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ir (Iridium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ir (Iridium) atom'
                                        COUNT_TYPE(77) = COUNT_TYPE(77) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'I '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as I (Iodine) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as I (Iodine) atom'
                                        COUNT_TYPE(53) = COUNT_TYPE(53) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'K') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Kr'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Kr (Krypton) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Kr (Krypton) atom'
                                        COUNT_TYPE(36) = COUNT_TYPE(36) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'K '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as K (Potassium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as K (Potassium) atom'
                                        COUNT_TYPE(19) = COUNT_TYPE(19) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'L') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('i', 'I')
                                        ATOM_TYPE(i) = 'Li'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Li (Lithium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Li (Lithium) atom'
                                        COUNT_TYPE(3) = COUNT_TYPE(3) + COUNT_NAME(i)
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'La'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as La (Lantanum) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as La (Lantanum) atom'
                                        COUNT_TYPE(57) = COUNT_TYPE(57) + COUNT_NAME(i)
                                CASE('u', 'U')
                                        ATOM_TYPE(i) = 'Lu'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Lu (Lutetium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Lu (Lutetium) atom'
                                        COUNT_TYPE(71) = COUNT_TYPE(71) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = ' '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'M') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('g', 'G')
                                        ATOM_TYPE(i) = 'Mg'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Mg (Magnesium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Mg (Magnesium) atom'
                                        COUNT_TYPE(12) = COUNT_TYPE(12) + COUNT_NAME(i)
                                CASE('n', 'N')
                                        ATOM_TYPE(i) = 'Mn'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Mn (Manganese) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Mn (Manganese) atom'
                                        COUNT_TYPE(25) = COUNT_TYPE(25) + COUNT_NAME(i)
                                CASE('o', 'O')
                                        ATOM_TYPE(i) = 'Mo'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Mo (Molybdenum) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Mo (Molybdenum) atom'
                                        COUNT_TYPE(42) = COUNT_TYPE(42) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = ' '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
                               
                        ELSE IF (ATOM_N(i)(:1) == 'N') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Ne'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ne (Neon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ne (Neon) atom'
                                        COUNT_TYPE(10) = COUNT_TYPE(10) + COUNT_NAME(i)
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'Na'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Na (Sodium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Na (Sodium) atom'
                                        COUNT_TYPE(11) = COUNT_TYPE(11) + COUNT_NAME(i)
                                CASE('i', 'I')
                                        ATOM_TYPE(i) = 'Ni'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ni (Nickel) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ni (Nickel) atom'
                                        COUNT_TYPE(28) = COUNT_TYPE(28) + COUNT_NAME(i)
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Nb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Nb (Neobium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Nb (Neobium) atom'
                                        COUNT_TYPE(41) = COUNT_TYPE(41) + COUNT_NAME(i)
                                CASE('d', 'D')
                                        ATOM_TYPE(i) = 'Nd'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Nd (Neodymium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Nd (Neodymium) atom'
                                        COUNT_TYPE(60) = COUNT_TYPE(60) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'N '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as N (Nitrogen) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as N (Nitrogen) atom'
                                        COUNT_TYPE(7) = COUNT_TYPE(7) + COUNT_NAME(i)
                                END SELECT                        
        
                        ELSE IF (ATOM_N(i)(:1) == 'O') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('s', 'S')
                                        ATOM_TYPE(i) = 'Os'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Os (Osmium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Os (Osmium) atom'
                                        COUNT_TYPE(76) = COUNT_TYPE(76) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'O '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as O (Oxygen) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as O (Oxygen) atom'
                                        COUNT_TYPE(8) = COUNT_TYPE(8) + COUNT_NAME(i)
                                END SELECT                         
        
                        ELSE IF (ATOM_N(i)(:1) == 'P') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('o', 'O')
                                        ATOM_TYPE(i) = 'Po'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Po (Polonium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Po (Polonium) atom'
                                        COUNT_TYPE(84) = COUNT_TYPE(84) + COUNT_NAME(i)
                                CASE('m', 'M')
                                        ATOM_TYPE(i) = 'Pm'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Pm (Promethium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Pm (Promethium) atom'
                                        COUNT_TYPE(61) = COUNT_TYPE(61) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Pr'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Pr (Praseodymium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Pr (Praseodymium) atom'
                                        COUNT_TYPE(59) = COUNT_TYPE(59) + COUNT_NAME(i)
                                CASE('t', 'T')
                                        ATOM_TYPE(i) = 'Pt'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Pt (Platinum) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Pt (Platinum) atom'
                                        COUNT_TYPE(78) = COUNT_TYPE(78) + COUNT_NAME(i)
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Pb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Pb (Lead) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Pb (Lead) atom'
                                        COUNT_TYPE(82) = COUNT_TYPE(82) + COUNT_NAME(i)
                                CASE('d', 'D')
                                        ATOM_TYPE(i) = 'Pd'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Pd (Palladium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Pd (Palladium) atom'
                                        COUNT_TYPE(46) = COUNT_TYPE(46) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'P '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as P (Phosphorus) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as P (Phosphorus) atom'
                                        COUNT_TYPE(15) = COUNT_TYPE(15) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'R') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Rb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Rb (Rubidium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Rb (Rubidium) atom'
                                        COUNT_TYPE(37) = COUNT_TYPE(37) + COUNT_NAME(i)
                                CASE('u', 'U')
                                        ATOM_TYPE(i) = 'Ru'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ru (Ruthenium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ru (Ruthenium) atom'
                                        COUNT_TYPE(44) = COUNT_TYPE(44) + COUNT_NAME(i)
                                CASE('h', 'H')
                                        ATOM_TYPE(i) = 'Rh'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Rh (Rhodium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Rh (Rhodium) atom'
                                        COUNT_TYPE(45) = COUNT_TYPE(45) + COUNT_NAME(i)
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Re'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Re (Rhenium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Re (Rhenium) atom'
                                        COUNT_TYPE(75) = COUNT_TYPE(75) + COUNT_NAME(i)
                                CASE('n', 'N')
                                        ATOM_TYPE(i) = 'Rn'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Rn (Radon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Rn (Radon) atom'
                                        COUNT_TYPE(86) = COUNT_TYPE(86) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = '  '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                        
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'S') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('i', 'I')
                                        ATOM_TYPE(i) = 'Si'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Si (Silicon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Si (Silicon) atom'
                                        COUNT_TYPE(14) = COUNT_TYPE(14) + COUNT_NAME(i)
                                CASE('c', 'C')
                                        ATOM_TYPE(i) = 'Sc'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Sc (Scandium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Sc (Scandium) atom'
                                        COUNT_TYPE(21) = COUNT_TYPE(21) + COUNT_NAME(i)
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Se'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Se (Selenium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Se (Selenium) atom'
                                        COUNT_TYPE(34) = COUNT_TYPE(34) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Sr'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Sr (Strontium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Sr (Strontium) atom'
                                        COUNT_TYPE(38) = COUNT_TYPE(38) + COUNT_NAME(i)
                                CASE('n', 'N')
                                        ATOM_TYPE(i) = 'Sn'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Sn (Tin) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Sn (Tin) atom'
                                        COUNT_TYPE(50) = COUNT_TYPE(50) + COUNT_NAME(i)
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Sb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Sb (Antimony) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Sb (Antimony) atom'
                                        COUNT_TYPE(51) = COUNT_TYPE(51) + COUNT_NAME(i)
                                CASE('m', 'M')
                                        ATOM_TYPE(i) = 'Sm'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Sm (Samarium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Sm (Samarium) atom'
                                        COUNT_TYPE(62) = COUNT_TYPE(62) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'S '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as S (Sulfur) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as S (Sulfur) atom'
                                        COUNT_TYPE(16) = COUNT_TYPE(16) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'T') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('i', 'I')
                                        ATOM_TYPE(i) = 'Ti'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ti (Titanium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ti (Titanium) atom'
                                        COUNT_TYPE(22) = COUNT_TYPE(22) + COUNT_NAME(i)
                                CASE('c', 'C')
                                        ATOM_TYPE(i) = 'Tc'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Tc (Technetium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Tc (Technetium) atom'
                                        COUNT_TYPE(43) = COUNT_TYPE(43) + COUNT_NAME(i)
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Te'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Te (Tellurium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Te (Tellurium) atom'
                                        COUNT_TYPE(52) = COUNT_TYPE(52) + COUNT_NAME(i)
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Tb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Tb (Terbium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Tb (Terbium) atom'
                                        COUNT_TYPE(65) = COUNT_TYPE(65) + COUNT_NAME(i)
                                CASE('m', 'M')
                                        ATOM_TYPE(i) = 'Tm'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Tm (Thulium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Tm (Thulium) atom'
                                        COUNT_TYPE(69) = COUNT_TYPE(69) + COUNT_NAME(i)
                                CASE('a', 'A')
                                        ATOM_TYPE(i) = 'Ta'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Ta (Tantalum) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Ta (Tantalum) atom'
                                        COUNT_TYPE(73) = COUNT_TYPE(73) + COUNT_NAME(i)
                                CASE('l', 'L')
                                        ATOM_TYPE(i) = 'Tl'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Tl (Thallium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Tl (Thallium) atom'
                                        COUNT_TYPE(81) = COUNT_TYPE(81) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = '  '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'V') THEN
                                ATOM_TYPE(i) = 'V  '
                                WRITE(outfile_unit,*) ATOM_N(i), 'assigned as V (Vanadium) atom'
                                WRITE(*,*) ATOM_N(i), 'assigned as V (Vanadium) atom'
                                COUNT_TYPE(23) = COUNT_TYPE(23) + COUNT_NAME(i)
        
                        ELSE IF (ATOM_N(i)(:1) == 'W') THEN
                                ATOM_TYPE(i) = 'W  '
                                WRITE(outfile_unit,*) ATOM_N(i), 'assigned as W (Tungsten) atom'
                                WRITE(*,*) ATOM_N(i), 'assigned as W (Tungsten) atom'
                                COUNT_TYPE(74) = COUNT_TYPE(74) + COUNT_NAME(i)
        
                        ELSE IF (ATOM_N(i)(:1) == 'X') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('e', 'E')
                                        ATOM_TYPE(i) = 'Xe'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Xe (Xenon) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Xe (Xenon) atom'
                                        COUNT_TYPE(54) = COUNT_TYPE(54) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = '  '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'Y') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('b', 'B')
                                        ATOM_TYPE(i) = 'Yb'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Yb (Ytterbium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Yb (Ytterbium) atom'
                                        COUNT_TYPE(70) = COUNT_TYPE(70) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = 'Y '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Y (Yttrium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Y (Yttrium) atom'
                                        COUNT_TYPE(39) = COUNT_TYPE(39) + COUNT_NAME(i)
                                END SELECT
        
                        ELSE IF (ATOM_N(i)(:1) == 'Z') THEN
                                SELECT CASE (ATOM_N(i)(2:2))
                                CASE('n', 'N')
                                        ATOM_TYPE(i) = 'Zn'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Zn (Zinc) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Zn (Zinc) atom'
                                        COUNT_TYPE(30) = COUNT_TYPE(30) + COUNT_NAME(i)
                                CASE('r', 'R')
                                        ATOM_TYPE(i) = 'Zr'
                                        WRITE(outfile_unit,*) ATOM_N(i), 'assigned as Zr (Zirconium) atom'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Zr (Zirconium) atom'
                                        COUNT_TYPE(40) = COUNT_TYPE(40) + COUNT_NAME(i)
                                CASE DEFAULT
                                        ATOM_TYPE(i) = '  '
                                        WRITE(outfile_unit,*) ATOM_N(i), 'not assigned'
                                        WRITE(*,*) ATOM_N(i), 'not assigned'
                                END SELECT
                        ELSE
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END IF
                END DO

                WRITE(*,*) 'Done ...................................'
                WRITE(*,*) '----------------------------------------'

                ! Count of atom type before any change
                WRITE(*,*) 'Found:'
                WRITE(*,*) '------'

                DO i = 1, size(COUNT_TYPE)
                        IF (COUNT_TYPE(i) /= 0) THEN
                                WRITE(*,*) COUNT_TYPE(i), ptable(i) % symbol
                        END IF
                END DO
                WRITE(*,*) '  ' 

                ! Gestion des not assigned

                IF (SUM(COUNT_TYPE) /= n_atoms) THEN
                        WRITE(*,'(I0, A)') n_atoms - SUM(COUNT_TYPE), ' atoms are not assigned'
                
                        DO i = 1, size(ATOM_N)
                                IF (ATOM_TYPE(i) == '  ') THEN
                                        IF (file_input .EQV. .TRUE.) THEN
                                                READ(tamsinput_unit,*) crap
                                                READ(tamsinput_unit,*) answer
                                        ELSE
                                                WRITE(*,'(A, A, A)') 'What is the atom type of ', ATOM_N(i), '?'
                                                WRITE(tamsinput_unit,'(A, A, A)') 'What is the atom type of ', ATOM_N(i), '?'
                                                READ(5,*) answer
                                                WRITE(tamsinput_unit,*) answer
                                        END IF
                                        ATOM_TYPE(i) = answer
                                        loc = FINDLOC(ptable % symbol, value=answer)
                                        COUNT_TYPE(loc) = COUNT_TYPE(loc) + COUNT_NAME(i)
                                END IF
                        END DO
                ENDIF

                WRITE(*,*) 'Update..'
                WRITE(*,*) 'Found:'
                WRITE(*,*) '------'

                DO i = 1, size(COUNT_TYPE)
                        IF (COUNT_TYPE(i) /= 0) THEN
                                WRITE(*,*) COUNT_TYPE(i), ptable(i) % symbol
                        END IF
                END DO
                
                WRITE(*,*) '  '

                ! 1 - Gestion des types d'atomes type Cr etc ..
                !   - Update du tableau COUNT_TYPE
                WRITE(*,*) 'Do you want to change something (y/n)?'
                READ(5,*) answer
                IF (answer == 'y' .OR. answer == 'Y') THEN
                        WRITE(*,*) 'okok'
                ENDIF
                !DO i = 1, size(COUNT_TYPE)
                ! 3 - Compte final ( a faire aussi avant pt. 1)

                ! 4 - WRITE(*,*) 'Do you want to change something ?'
                ! 5 - exclude from bond recognition (li etc..)

        !WRITE(*,*) ATOM_TYPE
        !IF (ANY(ATOM_TYPE == 'Cr')) THEN
        !        WRITE(*,*) 'Any Cr atom in the system ?'
        !        READ(*,*) answer
        !        IF (answer == 'n') THEN
        !                WRITE(*,*) 'Do you want to replace Cr assigned atoms by carbon atoms'
        !                READ(*,*) answer
        !                IF (answer == 'y') THEN
        !                        WHERE (ATOM_TYPE == 'Cr')
        !                                ATOM_TYPE = 'C'
        !                        END WHERE
        !                ELSE
        !                        WRITE(*,*) 'Cr atoms replace as not assigned'
        !                        WHERE (ATOM_TYPE == 'Cr')
        !                                ATOM_TYPE = '  '
        !                        END WHERE
        !                END IF
        !        END IF
        !END IF

        !WRITE(*,*) ATOM_TYPE
        !CALL sort(ATOM_NAME, ATOM_INDEX)
        !OPEN(unit =12, file = 'res2.dat', iostat=io)
        !DO i = 1, n_atoms
        !   WRITE(12, *) ATOM_NAME(i)
        !DO i = 1, n_atoms
                !IF (ANY(ATOM_NAME(i) == ANAME)) THEN
        !if ([(any(ATOM_NAME(i) == ptable % symbol), i = 1, n_atoms)]) then
        !                print *, "YES"
        !print *, ATOM_NAME(23)
                !WRITE(*,*) 'match'
                !END IF

        END SUBROUTINE SORT_ATOM_NAME
END MODULE
