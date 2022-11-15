MODULE MOLECULAR_RECOGNITION

        USE COVALENT_RADIUS
        USE MD_STUFF
        USE LEXICAL_SORT

        CONTAINS
        SUBROUTINE SORT_ATOM_NAME() !xyz_unit, xyz_filename, n_atoms)       
                IMPLICIT NONE
                INTEGER :: i, io, c, m, counter_same_name 
                !INTEGER, INTENT(IN) :: xyz_unit, n_atoms
                !DOUBLE PRECISION, DIMENSION(n_atoms,3) :: POS
                !CHARACTER(LEN = 10), DIMENSION(n_atoms) :: ATOM_NAME
                CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_N, TEMP_ARRAY, ATOM_TYPE
                INTEGER, ALLOCATABLE, DIMENSION(:) :: ATOM_INDEX, COUNT_NAME, TMP_ARRAY
                INTEGER, DIMENSION(n_elem) :: COUNT_TYPE
                CHARACTER(LEN = 1) :: answer
                !CHARACTER(LEN = 30), INTENT(IN) :: xyz_filename
                CALL init_periodic_table()

        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        ALLOCATE(ATOM_INDEX(n_atoms))
        write(*,*) 'coucou 0'
        ! Read one step to get the atom name       
        CALL read_xyz() !xyz_unit, n_atoms, ATOM_NAME, POS)
        write(*,*) 'coucou'
        ! Sort all the atom name 
        CALL sort(ATOM_NAME, ATOM_INDEX)

        ! Remove duplicates
        c = 1
        counter_same_name = 0
        ALLOCATE(ATOM_N(c))
        ALLOCATE(COUNT_NAME(c))

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

        ATOM_N(c) = ATOM_NAME(n_atoms)
        COUNT_NAME(c) = counter_same_name + 1
        
        ALLOCATE(ATOM_TYPE(size(ATOM_N)))

        WRITE(*,*) size(ATOM_N), 'atoms type found ...'
        WRITE(*,*) 'Starting assigning type to each atom ...'
        WRITE(*,*) '----------------------------------------'
        
        write(*,*) size(COUNT_NAME)
        write(*,*) size(ATOM_N)
        
        DO i = 1, size(ATOM_N)
                IF (ATOM_N(i)(:1) == 'C') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('l', 'L')
                                ATOM_TYPE(i) = 'Cl'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cl (Chlorine) atom'
                                COUNT_TYPE(17) = COUNT_TYPE(17) + 1
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'Ca'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ca (Calcium) atom'
                                COUNT_TYPE(20) = COUNT_TYPE(20) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Cr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cr (Chromium) atom'
                                COUNT_TYPE(24) = COUNT_TYPE(24) + 1
                        CASE('o', 'O')
                                ATOM_TYPE(i) = 'Co'
                                WRITE(*,*) ATOM_N(i), 'assigned as Co (Cobalt) atom'
                                COUNT_TYPE(27) = COUNT_TYPE(27) + 1
                        CASE('u', 'U')
                                ATOM_TYPE(i) = 'Cu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cu (Copper) atom'
                                COUNT_TYPE(29) = COUNT_TYPE(29) + 1
                        CASE('d', 'D')
                                ATOM_TYPE(i) = 'Cd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cd (Cadmium) atom'
                                COUNT_TYPE(48) = COUNT_TYPE(48) + 1
                        CASE('s', 'S')
                                ATOM_TYPE(i) = 'Cs'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cs (Cesium) atom'
                                COUNT_TYPE(55) = COUNT_TYPE(55) + 1
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Ce'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ce (Cerium) atom'
                                COUNT_TYPE(58) = COUNT_TYPE(58) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'C '
                                WRITE(*,*) ATOM_N(i), 'assigned as C (Carbon) atom'
                                COUNT_TYPE(6) = COUNT_TYPE(6) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'A') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('l', 'L')
                                ATOM_TYPE(i) = 'Al'
                                WRITE(*,*) ATOM_N(i), 'assigned as Al (Aluminium) atom'
                                COUNT_TYPE(13) = COUNT_TYPE(13) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Ar'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ar (Argon) atom'
                                COUNT_TYPE(18) = COUNT_TYPE(18) + 1
                        CASE('u', 'U')
                                ATOM_TYPE(i) = 'Au'
                                WRITE(*,*) ATOM_N(i), 'assigned as Au (Gold) atom'
                                COUNT_TYPE(79) = COUNT_TYPE(79) + 1
                        CASE('s', 'S')
                                ATOM_TYPE(i) = 'As'
                                WRITE(*,*) ATOM_N(i), 'assigned as As (Arsenic) atom'
                                COUNT_TYPE(33) = COUNT_TYPE(33) + 1
                        CASE('g', 'G')
                                ATOM_TYPE(i) = 'Ag'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ag (Silver) atom'
                                COUNT_TYPE(47) = COUNT_TYPE(47) + 1
                        CASE('t', 'T')
                                ATOM_TYPE(i) = 'At'
                                WRITE(*,*) ATOM_N(i), 'assigned as At (Astatine) atom'
                                COUNT_TYPE(85) = COUNT_TYPE(85) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'B') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Be'
                                WRITE(*,*) ATOM_N(i), 'assigned as Be (Beryllium) atom'
                                COUNT_TYPE(4) = COUNT_TYPE(4) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Br'
                                WRITE(*,*) ATOM_N(i), 'assigned as Br (Br) atom'
                                COUNT_TYPE(35) = COUNT_TYPE(35) + 1
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'Ba'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ba (Barium) atom'
                                COUNT_TYPE(56) = COUNT_TYPE(56) + 1
                        CASE('i', 'I')
                                ATOM_TYPE(i) = 'Bi'
                                WRITE(*,*) ATOM_N(i), 'assigned as Bi (Bismuth) atom'
                                COUNT_TYPE(83) = COUNT_TYPE(83) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'B '
                                WRITE(*,*) ATOM_N(i), 'assigned as B (Boron) atom'
                                COUNT_TYPE(5) = COUNT_TYPE(5) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'D') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('y', 'Y')
                                ATOM_TYPE(i) = 'Dy'
                                WRITE(*,*) ATOM_N(i), 'assigned as Dy (Dysprosium) atom'
                                COUNT_TYPE(66) = COUNT_TYPE(66) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'E') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('u', 'U')
                                ATOM_TYPE(i) = 'Eu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Eu (Europium) atom'
                                COUNT_TYPE(63) = COUNT_TYPE(63) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Er'
                                WRITE(*,*) ATOM_N(i), 'assigned as Eu (Erbium) atom'
                                COUNT_TYPE(68) = COUNT_TYPE(68) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'F') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Fe'
                                WRITE(*,*) ATOM_N(i), 'assigned as Fe (Iron) atom'
                                COUNT_TYPE(26) = COUNT_TYPE(26) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'F '
                                WRITE(*,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                                COUNT_TYPE(9) = COUNT_TYPE(9) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'G') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'Ga'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ga (Gallium) atom'
                                COUNT_TYPE(31) = COUNT_TYPE(31) + 1
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Ge'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ge (Germanium) atom'
                                COUNT_TYPE(32) = COUNT_TYPE(32) + 1
                        CASE('d', 'D') 
                                ATOM_TYPE(i) = 'Gd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Gd (Gadolinium) atom'
                                COUNT_TYPE(64) = COUNT_TYPE(64) + 1 
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'H') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'He'
                                WRITE(*,*) ATOM_N(i), 'assigned as He (Helium) atom'
                                COUNT_TYPE(2) = COUNT_TYPE(2) + 1
                        CASE('o', 'O')
                                ATOM_TYPE(i) = 'Ho'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ho (Holmium) atom'
                                COUNT_TYPE(67) = COUNT_TYPE(67) + 1
                        CASE('f', 'F')
                                ATOM_TYPE(i) = 'Hf'
                                WRITE(*,*) ATOM_N(i), 'assigned as Hf (Hafnium) atom'
                                COUNT_TYPE(72) = COUNT_TYPE(72) + 1
                        CASE('g', 'G')
                                ATOM_TYPE(i) = 'Hg'
                                WRITE(*,*) ATOM_N(i), 'assigned as Hg (Mercury) atom'
                                COUNT_TYPE(80) = COUNT_TYPE(80) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'H '
                                WRITE(*,*) ATOM_N(i), 'assigned as H (Hydrogen) atom'
                                COUNT_TYPE(1) = COUNT_TYPE(1) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'I') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('n', 'N')
                                ATOM_TYPE(i) = 'In'
                                WRITE(*,*) ATOM_N(i), 'assigned as In (Indium) atom'
                                COUNT_TYPE(49) = COUNT_TYPE(49) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Ir'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ir (Iridium) atom'
                                COUNT_TYPE(77) = COUNT_TYPE(77) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'I '
                                WRITE(*,*) ATOM_N(i), 'assigned as I (Iodine) atom'
                                COUNT_TYPE(53) = COUNT_TYPE(53) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'K') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Kr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Kr (Krypton) atom'
                                COUNT_TYPE(36) = COUNT_TYPE(36) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'K '
                                WRITE(*,*) ATOM_N(i), 'assigned as K (Potassium) atom'
                                COUNT_TYPE(19) = COUNT_TYPE(19) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'L') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i', 'I')
                                ATOM_TYPE(i) = 'Li'
                                WRITE(*,*) ATOM_N(i), 'assigned as Li (Lithium) atom'
                                COUNT_TYPE(3) = COUNT_TYPE(3) + 1
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'La'
                                WRITE(*,*) ATOM_N(i), 'assigned as La (Lantanum) atom'
                                COUNT_TYPE(57) = COUNT_TYPE(57) + 1
                        CASE('u', 'U')
                                ATOM_TYPE(i) = 'Lu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Lu (Lutetium) atom'
                                COUNT_TYPE(71) = COUNT_TYPE(71) + 1 
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'M') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('g', 'G')
                                ATOM_TYPE(i) = 'Mg'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mg (Magnesium) atom'
                                COUNT_TYPE(12) = COUNT_TYPE(12) + 1
                        CASE('n', 'N')
                                ATOM_TYPE(i) = 'Mn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mn (Manganese) atom'
                                COUNT_TYPE(25) = COUNT_TYPE(25) + 1
                        CASE('o', 'O')
                                ATOM_TYPE(i) = 'Mo'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mo (Molybdenum) atom'
                                COUNT_TYPE(42) = COUNT_TYPE(42) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT
                       
                ELSE IF (ATOM_N(i)(:1) == 'N') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Ne'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ne (Neon) atom'
                                COUNT_TYPE(10) = COUNT_TYPE(10) + 1
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'Na'
                                WRITE(*,*) ATOM_N(i), 'assigned as Na (Sodium) atom'
                                COUNT_TYPE(11) = COUNT_TYPE(11) + 1
                        CASE('i', 'I')
                                ATOM_TYPE(i) = 'Ni'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ni (Nickel) atom'
                                COUNT_TYPE(28) = COUNT_TYPE(28) + 1
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Nb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Nb (Neobium) atom'
                                COUNT_TYPE(41) = COUNT_TYPE(41) + 1
                        CASE('d', 'D')
                                ATOM_TYPE(i) = 'Nd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Nd (Neodymium) atom'
                                COUNT_TYPE(60) = COUNT_TYPE(60) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'N '
                                WRITE(*,*) ATOM_N(i), 'assigned as N (Nitrogen) atom'
                                COUNT_TYPE(7) = COUNT_TYPE(7) + 1
                        END SELECT                        

                ELSE IF (ATOM_N(i)(:1) == 'O') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('s', 'S')
                                ATOM_TYPE(i) = 'Os'
                                WRITE(*,*) ATOM_N(i), 'assigned as Os (Osmium) atom'
                                COUNT_TYPE(76) = COUNT_TYPE(76) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'O '
                                WRITE(*,*) ATOM_N(i), 'assigned as O (Oxygen) atom'
                                COUNT_TYPE(8) = COUNT_TYPE(8) + 1
                        END SELECT                         

                ELSE IF (ATOM_N(i)(:1) == 'P') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('o', 'O')
                                ATOM_TYPE(i) = 'Po'
                                WRITE(*,*) ATOM_N(i), 'assigned as Po (Polonium) atom'
                                COUNT_TYPE(84) = COUNT_TYPE(84) + 1
                        CASE('m', 'M')
                                ATOM_TYPE(i) = 'Pm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pm (Promethium) atom'
                                COUNT_TYPE(61) = COUNT_TYPE(61) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Pr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pr (Praseodymium) atom'
                                COUNT_TYPE(59) = COUNT_TYPE(59) + 1
                        CASE('t', 'T')
                                ATOM_TYPE(i) = 'Pt'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pt (Platinum) atom'
                                COUNT_TYPE(78) = COUNT_TYPE(78) + 1
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Pb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pb (Lead) atom'
                                COUNT_TYPE(82) = COUNT_TYPE(82) + 1
                        CASE('d', 'D')
                                ATOM_TYPE(i) = 'Pd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pd (Palladium) atom'
                                COUNT_TYPE(46) = COUNT_TYPE(46) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'P '
                                WRITE(*,*) ATOM_N(i), 'assigned as P (Phosphorus) atom'
                                COUNT_TYPE(15) = COUNT_TYPE(15) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'R') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Rb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rb (Rubidium) atom'
                                COUNT_TYPE(37) = COUNT_TYPE(37) + 1
                        CASE('u', 'U')
                                ATOM_TYPE(i) = 'Ru'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ru (Ruthenium) atom'
                                COUNT_TYPE(44) = COUNT_TYPE(44) + 1
                        CASE('h', 'H')
                                ATOM_TYPE(i) = 'Rh'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rh (Rhodium) atom'
                                COUNT_TYPE(45) = COUNT_TYPE(45) + 1
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Re'
                                WRITE(*,*) ATOM_N(i), 'assigned as Re (Rhenium) atom'
                                COUNT_TYPE(75) = COUNT_TYPE(75) + 1
                        CASE('n', 'N')
                                ATOM_TYPE(i) = 'Rn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rn (Radon) atom'
                                COUNT_TYPE(86) = COUNT_TYPE(86) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'S') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i', 'I')
                                ATOM_TYPE(i) = 'Si'
                                WRITE(*,*) ATOM_N(i), 'assigned as Si (Silicon) atom'
                                COUNT_TYPE(14) = COUNT_TYPE(14) + 1
                        CASE('c', 'C')
                                ATOM_TYPE(i) = 'Sc'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sc (Scandium) atom'
                                COUNT_TYPE(21) = COUNT_TYPE(21) + 1
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Se'
                                WRITE(*,*) ATOM_N(i), 'assigned as Se (Selenium) atom'
                                COUNT_TYPE(34) = COUNT_TYPE(34) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Sr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sr (Strontium) atom'
                                COUNT_TYPE(38) = COUNT_TYPE(38) + 1
                        CASE('n', 'N')
                                ATOM_TYPE(i) = 'Sn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sn (Tin) atom'
                                COUNT_TYPE(50) = COUNT_TYPE(50) + 1
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Sb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sb (Antimony) atom'
                                COUNT_TYPE(51) = COUNT_TYPE(51) + 1
                        CASE('m', 'M')
                                ATOM_TYPE(i) = 'Sm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sm (Samarium) atom'
                                COUNT_TYPE(62) = COUNT_TYPE(62) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'S '
                                WRITE(*,*) ATOM_N(i), 'assigned as S (Sulfur) atom'
                                COUNT_TYPE(16) = COUNT_TYPE(16) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'T') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i', 'I')
                                ATOM_TYPE(i) = 'Ti'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ti (Titanium) atom'
                                COUNT_TYPE(22) = COUNT_TYPE(22) + 1
                        CASE('c', 'C')
                                ATOM_TYPE(i) = 'Tc'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tc (Technetium) atom'
                                COUNT_TYPE(43) = COUNT_TYPE(43) + 1
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Te'
                                WRITE(*,*) ATOM_N(i), 'assigned as Te (Tellurium) atom'
                                COUNT_TYPE(52) = COUNT_TYPE(52) + 1
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Tb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tb (Terbium) atom'
                                COUNT_TYPE(65) = COUNT_TYPE(65) + 1
                        CASE('m', 'M')
                                ATOM_TYPE(i) = 'Tm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tm (Thulium) atom'
                                COUNT_TYPE(69) = COUNT_TYPE(69) + 1
                        CASE('a', 'A')
                                ATOM_TYPE(i) = 'Ta'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ta (Tantalum) atom'
                                COUNT_TYPE(73) = COUNT_TYPE(73) + 1
                        CASE('l', 'L')
                                ATOM_TYPE(i) = 'Tl'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tl (Thallium) atom'
                                COUNT_TYPE(81) = COUNT_TYPE(81) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'V') THEN
                        ATOM_TYPE(i) = 'V  '
                        WRITE(*,*) ATOM_N(i), 'assigned as V (Vanadium) atom'
                        COUNT_TYPE(23) = COUNT_TYPE(23) + 1

                ELSE IF (ATOM_N(i)(:1) == 'W') THEN
                        ATOM_TYPE(i) = 'W  '
                        WRITE(*,*) ATOM_N(i), 'assigned as W (Tungsten) atom'
                        COUNT_TYPE(74) = COUNT_TYPE(74) + 1

                ELSE IF (ATOM_N(i)(:1) == 'X') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e', 'E')
                                ATOM_TYPE(i) = 'Xe'
                                WRITE(*,*) ATOM_N(i), 'assigned as Xe (Xenon) atom'
                                COUNT_TYPE(54) = COUNT_TYPE(54) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'Y') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('b', 'B')
                                ATOM_TYPE(i) = 'Yb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Yb (Ytterbium) atom'
                                COUNT_TYPE(70) = COUNT_TYPE(70) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'Y '
                                WRITE(*,*) ATOM_N(i), 'assigned as Y (Yttrium) atom'
                                COUNT_TYPE(39) = COUNT_TYPE(39) + 1
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'Z') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('n', 'N')
                                ATOM_TYPE(i) = 'Zn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Zn (Zinc) atom'
                                COUNT_TYPE(30) = COUNT_TYPE(30) + 1
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Zr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Zr (Zirconium) atom'
                                COUNT_TYPE(40) = COUNT_TYPE(40) + 1
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT
                ELSE
                        ATOM_TYPE(i) = '  '
                        WRITE(*,*) ATOM_N(i), 'not assigned'
                END IF
        END DO
        ! Final count of atom type before any change


        ! 1 - Gestion des types d'atomes type Cr etc ..
        !   - Update du tableau COUNT_TYPE

        ! 2 - Gestion des not assigned 
        !   - Update du tableau COUNT_TYPE

        ! 3 - Compte final ( a faire aussi avant pt. 1)
        !   - Relier COUNT_TYPE au tableau ATOM_NAME
        !   - Lors du sort compter le nombre d'atomes du même type et multiplier

        ! 4 - WRITE(*,*) 'Do you want to change something ?'


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
