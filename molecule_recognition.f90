MODULE MOLECULAR_RECOGNITION

        USE COVALENT_RADIUS
        USE MD_STUFF
        USE LEXICAL_SORT

        CONTAINS
        SUBROUTINE SORT_ATOM_NAME(xyz_unit, xyz_filename, n_atoms)       
                IMPLICIT NONE
                INTEGER :: i, io, c, m 
                INTEGER, INTENT(IN) :: xyz_unit, n_atoms
                DOUBLE PRECISION, DIMENSION(n_atoms,3) :: POS
                CHARACTER(LEN = 10), DIMENSION(n_atoms) :: ATOM_NAME
                CHARACTER(LEN = 10), ALLOCATABLE, DIMENSION(:) :: ATOM_N, TEMP_ARRAY, ATOM_TYPE
                INTEGER, ALLOCATABLE, DIMENSION(:) :: ATOM_INDEX
                CHARACTER(LEN = 1) :: answer
                CHARACTER(LEN = 30), INTENT(IN) :: xyz_filename
                CALL init_periodic_table()

        OPEN(unit = xyz_unit, file = xyz_filename, status='old', iostat=io)
        ALLOCATE(ATOM_INDEX(n_atoms))
        write(*,*) 'allocate done!'
        ! Read one step to get the atom name       
        CALL read_xyz(xyz_unit, n_atoms, ATOM_NAME, POS)

        ! Sort all the atom name 
        CALL sort(ATOM_NAME, ATOM_INDEX)

        ! Remove duplicates
        c = 1
        ALLOCATE(ATOM_N(c))

        DO i = 1, n_atoms-1 
                IF (ATOM_NAME(i) /= ATOM_NAME(i+1)) THEN
                        ATOM_N(c) = ATOM_NAME(i)
                        ALLOCATE(TEMP_ARRAY(size(ATOM_N)+1))
                        TEMP_ARRAY(1:size(ATOM_N)) = ATOM_N
                        DEALLOCATE(ATOM_N)
                        ALLOCATE(ATOM_N(size(TEMP_ARRAY)))
                        ATOM_N(1:size(TEMP_ARRAY)) = TEMP_ARRAY
                        DEALLOCATE(TEMP_ARRAY)
                        c = c + 1 
                END IF
        END DO
        ATOM_N(c) = ATOM_NAME(n_atoms)
        ALLOCATE(ATOM_TYPE(size(ATOM_N)))
        WRITE(*,*) size(ATOM_N), 'atoms type found ...'
        WRITE(*,*) 'Starting assigning type to each atom ...'
        WRITE(*,*) '----------------------------------------'

        DO i = 1, size(ATOM_N)
                IF (ATOM_N(i)(:1) == 'C') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('l')
                                ATOM_TYPE(i) = 'Cl'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cl (Chlorine) atom'
                        CASE('a')
                                ATOM_TYPE(i) = 'Ca'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ca (Calcium) atom'
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Cr'
                                        WRITE(*,*) ATOM_N(i), 'assigned as Cr (Chromium) atom'
                        CASE('o')
                                ATOM_TYPE(i) = 'Co'
                                WRITE(*,*) ATOM_N(i), 'assigned as Co (Cobalt) atom'
                        CASE('u')
                                ATOM_TYPE(i) = 'Cu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cu (Copper) atom'
                        CASE('d')
                                ATOM_TYPE(i) = 'Cd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cd (Cadmium) atom'
                        CASE('s')
                                ATOM_TYPE(i) = 'Cs'
                                WRITE(*,*) ATOM_N(i), 'assigned as Cs (Cesium) atom'
                        CASE('e')
                                ATOM_TYPE(i) = 'Ce'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ce (Cerium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'C '
                                WRITE(*,*) ATOM_N(i), 'assigned as C (Carbon) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'A') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('l')
                                ATOM_TYPE(i) = 'Al'
                                WRITE(*,*) ATOM_N(i), 'assigned as Al (Aluminium) atom'
                        CASE('r', 'R')
                                ATOM_TYPE(i) = 'Ar'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ar (Argon) atom'
                        CASE('u')
                                ATOM_TYPE(i) = 'Au'
                                WRITE(*,*) ATOM_N(i), 'assigned as Au (Gold) atom'
                        CASE('s')
                                ATOM_TYPE(i) = 'As'
                                WRITE(*,*) ATOM_N(i), 'assigned as As (Arsenic) atom'
                        CASE('g')
                                ATOM_TYPE(i) = 'Ag'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ag (Silver) atom'
                        CASE('t')
                                ATOM_TYPE(i) = 'At'
                                WRITE(*,*) ATOM_N(i), 'assigned as At (Astatine) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'B') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'Be'
                                WRITE(*,*) ATOM_N(i), 'assigned as Be (Beryllium) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Br'
                                WRITE(*,*) ATOM_N(i), 'assigned as Br (Br) atom'
                        CASE('a')
                                ATOM_TYPE(i) = 'Ba'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ba (Barium) atom'
                        CASE('i')
                                ATOM_TYPE(i) = 'Bi'
                                WRITE(*,*) ATOM_N(i), 'assigned as Bi (Bismuth) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'B '
                                WRITE(*,*) ATOM_N(i), 'assigned as B (Boron) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'D') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('y')
                                ATOM_TYPE(i) = 'Dy'
                                WRITE(*,*) ATOM_N(i), 'assigned as Dy (Dysprosium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'E') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('u')
                                ATOM_TYPE(i) = 'Eu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Eu (Europium) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Er'
                                WRITE(*,*) ATOM_N(i), 'assigned as Eu (Erbium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'F') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'Fe'
                                WRITE(*,*) ATOM_N(i), 'assigned as Fe (Iron) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'F '
                                WRITE(*,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'F') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'Fe'
                                WRITE(*,*) ATOM_N(i), 'assigned as Fe (Iron) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'F '
                                WRITE(*,*) ATOM_N(i), 'assigned as F (Fluorine) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'G') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('a')
                                ATOM_TYPE(i) = 'Ga'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ga (Gallium) atom'
                        CASE('e')
                                ATOM_TYPE(i) = 'Ge'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ge (Germanium) atom'
                        CASE('d') 
                                ATOM_TYPE(i) = 'Gd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Gd (Gadolinium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'H') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'He'
                                WRITE(*,*) ATOM_N(i), 'assigned as He (Helium) atom'
                        CASE('o')
                                ATOM_TYPE(i) = 'Ho'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ho (Holmium) atom'
                        CASE('f')
                                ATOM_TYPE(i) = 'Hf'
                                WRITE(*,*) ATOM_N(i), 'assigned as Hf (Hafnium) atom'
                        CASE('g')
                                ATOM_TYPE(i) = 'Hg'
                                WRITE(*,*) ATOM_N(i), 'assigned as Hg (Mercury) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'H '
                                WRITE(*,*) ATOM_N(i), 'assigned as H (Hydrogen) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'I') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('n')
                                ATOM_TYPE(i) = 'In'
                                WRITE(*,*) ATOM_N(i), 'assigned as In (Indium) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Ir'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ir (Iridium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'I '
                                WRITE(*,*) ATOM_N(i), 'assigned as I (Iodine) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'K') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('r')
                                ATOM_TYPE(i) = 'Kr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Kr (Krypton) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'I '
                                WRITE(*,*) ATOM_N(i), 'assigned as K (Potassium) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'L') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i')
                                ATOM_TYPE(i) = 'Li'
                                WRITE(*,*) ATOM_N(i), 'assigned as Li (Lithium) atom'
                        CASE('a')
                                ATOM_TYPE(i) = 'La'
                                WRITE(*,*) ATOM_N(i), 'assigned as La (Lantanum) atom'
                        CASE('u')
                                ATOM_TYPE(i) = 'Lu'
                                WRITE(*,*) ATOM_N(i), 'assigned as Lu (Lutetium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'M') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('g')
                                ATOM_TYPE(i) = 'Mg'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mg (Magnesium) atom'
                        CASE('n')
                                ATOM_TYPE(i) = 'Mn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mn (Manganese) atom'
                        CASE('o')
                                ATOM_TYPE(i) = 'Mo'
                                WRITE(*,*) ATOM_N(i), 'assigned as Mo (Molybdenum) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = ' '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT
                       
                ELSE IF (ATOM_N(i)(:1) == 'N') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'Ne'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ne (Neon) atom'
                        CASE('a')
                                ATOM_TYPE(i) = 'Na'
                                WRITE(*,*) ATOM_N(i), 'assigned as Na (Sodium) atom'
                        CASE('i')
                                ATOM_TYPE(i) = 'Ni'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ni (Nickel) atom'
                        CASE('b')
                                ATOM_TYPE(i) = 'Nb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Nb (Neobium) atom'
                        CASE('d')
                                ATOM_TYPE(i) = 'Nd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Nd (Neodymium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'N '
                                WRITE(*,*) ATOM_N(i), 'assigned as N (Nitrogen) atom'
                        END SELECT                        

                ELSE IF (ATOM_N(i)(:1) == 'O') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('s')
                                ATOM_TYPE(i) = 'Os'
                                WRITE(*,*) ATOM_N(i), 'assigned as Os (Osmium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'O '
                                WRITE(*,*) ATOM_N(i), 'assigned as O (Oxygen) atom'
                        END SELECT                         

                ELSE IF (ATOM_N(i)(:1) == 'P') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('o')
                                ATOM_TYPE(i) = 'Po'
                                WRITE(*,*) ATOM_N(i), 'assigned as Po (Polonium) atom'
                        CASE('m')
                                ATOM_TYPE(i) = 'Pm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pm (Promethium) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Pr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pr (Praseodymium) atom'
                        CASE('t')
                                ATOM_TYPE(i) = 'Pt'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pt (Platinum) atom'
                        CASE('b')
                                ATOM_TYPE(i) = 'Pb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pb (Lead) atom'
                        CASE('d')
                                ATOM_TYPE(i) = 'Pd'
                                WRITE(*,*) ATOM_N(i), 'assigned as Pd (Palladium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'P '
                                WRITE(*,*) ATOM_N(i), 'assigned as P (Phosphorus) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'R') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('b')
                                ATOM_TYPE(i) = 'Rb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rb (Rubidium) atom'
                        CASE('u')
                                ATOM_TYPE(i) = 'Ru'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ru (Ruthenium) atom'
                        CASE('h')
                                ATOM_TYPE(i) = 'Rh'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rh (Rhodium) atom'
                        CASE('e')
                                ATOM_TYPE(i) = 'Re'
                                WRITE(*,*) ATOM_N(i), 'assigned as Re (Rhenium) atom'
                        CASE('n')
                                ATOM_TYPE(i) = 'Rn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Rn (Radon) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'S') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i')
                                ATOM_TYPE(i) = 'Si'
                                WRITE(*,*) ATOM_N(i), 'assigned as Si (Silicon) atom'
                        CASE('c')
                                ATOM_TYPE(i) = 'Sc'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sc (Scandium) atom'
                        CASE('e')
                                ATOM_TYPE(i) = 'Se'
                                WRITE(*,*) ATOM_N(i), 'assigned as Se (Selenium) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Sr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sr (Strontium) atom'
                        CASE('n')
                                ATOM_TYPE(i) = 'Sn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sn (Tin) atom'
                        CASE('b')
                                ATOM_TYPE(i) = 'Sb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sb (Antimony) atom'
                        CASE('m')
                                ATOM_TYPE(i) = 'Sm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Sm (Samarium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'S '
                                WRITE(*,*) ATOM_N(i), 'assigned as S (Sulfur) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'T') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('i')
                                ATOM_TYPE(i) = 'Ti'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ti (Titanium) atom'
                        CASE('c')
                                ATOM_TYPE(i) = 'Tc'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tc (Technetium) atom'
                        CASE('e')
                                ATOM_TYPE(i) = 'Te'
                                WRITE(*,*) ATOM_N(i), 'assigned as Te (Tellurium) atom'
                        CASE('b')
                                ATOM_TYPE(i) = 'Tb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tb (Terbium) atom'
                        CASE('m')
                                ATOM_TYPE(i) = 'Tm'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tm (Thulium) atom'
                        CASE('a')
                                ATOM_TYPE(i) = 'Ta'
                                WRITE(*,*) ATOM_N(i), 'assigned as Ta (Tantalum) atom'
                        CASE('l')
                                ATOM_TYPE(i) = 'Tl'
                                WRITE(*,*) ATOM_N(i), 'assigned as Tl (Thallium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'V') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'V  '
                                WRITE(*,*) ATOM_N(i), 'assigned as V (Vanadium) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'W') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'W  '
                                WRITE(*,*) ATOM_N(i), 'assigned as W (Tungsten) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'X') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('e')
                                ATOM_TYPE(i) = 'Xe'
                                WRITE(*,*) ATOM_N(i), 'assigned as Xe (Xenon) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'Y') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('b')
                                ATOM_TYPE(i) = 'Yb'
                                WRITE(*,*) ATOM_N(i), 'assigned as Yb (Ytterbium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = 'Y '
                                WRITE(*,*) ATOM_N(i), 'assigned as Y (Yttrium) atom'
                        END SELECT

                ELSE IF (ATOM_N(i)(:1) == 'Z') THEN
                        SELECT CASE (ATOM_N(i)(2:2))
                        CASE('n')
                                ATOM_TYPE(i) = 'Zn'
                                WRITE(*,*) ATOM_N(i), 'assigned as Zn (Zinc) atom'
                        CASE('r')
                                ATOM_TYPE(i) = 'Zr'
                                WRITE(*,*) ATOM_N(i), 'assigned as Zr (Zirconium) atom'
                        CASE DEFAULT
                                ATOM_TYPE(i) = '  '
                                WRITE(*,*) ATOM_N(i), 'not assigned'
                        END SELECT
                ELSE
                        ATOM_TYPE(i) = '  '
                        WRITE(*,*) ATOM_N(i), 'not assigned'
                END IF
        END DO
        WRITE(*,*) ATOM_TYPE
        IF (ANY(ATOM_TYPE == 'Cr')) THEN
                WRITE(*,*) 'Any Cr atom in the system ?'
                READ(*,*) answer
                IF (answer == 'n') THEN
                        WHERE (ATOM_TYPE == 'Cr')
                                ATOM_TYPE = 'C'
                        end where
                end if
        END IF
        WRITE(*,*) ATOM_TYPE
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
