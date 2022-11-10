MODULE COVALENT_RADIUS
        IMPLICIT NONE
        PUBLIC :: init_periodic_table, atom, ptable, n_elem
!************************************************************        
        TYPE atom
                CHARACTER (LEN = 2) :: symbol
                INTEGER :: symbol_number
                REAL :: mass
                REAL :: covalent_radius
        END TYPE atom

        INTEGER, PARAMETER :: n_elem = 86

        TYPE (atom) :: ptable(n_elem)
        CONTAINS

!************************************************************
        SUBROUTINE init_periodic_table()

        ! Hydrogen
        ptable(1) % symbol = 'H '
        ptable(1) % symbol_number = 1
        ptable(1) % mass = 1.007825 
        ptable(1) % covalent_radius = 0.32 
        
        ! Helium
        ptable(2) % symbol = 'He'
        ptable(2) % symbol_number = 2
        ptable(2) % mass = 4.00260 
        ptable(2) % covalent_radius = 0.9300 
        
        ! Lithium
        ptable(3) % symbol = 'Li'
        ptable(3) % symbol_number = 3
        ptable(3) % mass = 7.01600 
        ptable(3) % covalent_radius = 1.2300 
        
        ! Beryllium
        ptable(4) % symbol = 'Be'
        ptable(4) % symbol_number = 4
        ptable(4) % mass = 9.01218 
        ptable(4) % covalent_radius = 0.9000 
        
        ! Boron
        ptable(5) % symbol = 'B '
        ptable(5) % symbol_number = 5
        ptable(5) % mass = 11.00931 
        ptable(5) % covalent_radius = 0.8200  
        
        ! Carbon
        ptable(6) % symbol = 'C '
        ptable(6) % symbol_number = 6
        ptable(6) % mass = 12.0000 
        ptable(6) % covalent_radius = 0.7700 
        
        ! Nitrogen
        ptable(7) % symbol = 'N '
        ptable(7) % symbol_number = 7
        ptable(7) % mass = 14.00307 
        ptable(7) % covalent_radius = 0.7500 
        
        ! Oxygen
        ptable(8) % symbol = 'O '
        ptable(8) % symbol_number = 8
        ptable(8) % mass = 15.99491 
        ptable(8) % covalent_radius = 0.7300 
        
        ! Fluorine
        ptable(9) % symbol = 'F '
        ptable(9) % symbol_number = 9
        ptable(9) % mass = 18.99840 
        ptable(9) % covalent_radius = 0.7200 
        
        ! Neon
        ptable(10) % symbol = 'Ne'
        ptable(10) % symbol_number = 10
        ptable(10) % mass = 19.99244 
        ptable(10) % covalent_radius = 0.7100 
        
        ! Sodium
        ptable(11) % symbol = 'Na'
        ptable(11) % symbol_number = 11
        ptable(11) % mass = 22.9898 
        ptable(11) % covalent_radius = 1.5400 

        ! Magnesium
        ptable(12) % symbol = 'Mg'
        ptable(12) % symbol_number = 12
        ptable(12) % mass = 23.98504 
        ptable(12) % covalent_radius = 1.3600 

        ! Aluminium
        ptable(13) % symbol = 'Al'
        ptable(13) % symbol_number = 13
        ptable(13) % mass = 26.98153 
        ptable(13) % covalent_radius = 1.1800 

        ! Silicon
        ptable(14) % symbol = 'Si'
        ptable(14) % symbol_number = 14
        ptable(14) % mass = 27.97693 
        ptable(14) % covalent_radius = 1.1100 

        ! Phosphorus
        ptable(15) % symbol = 'P '
        ptable(15) % symbol_number = 15
        ptable(15) % mass = 30.97376 
        ptable(15) % covalent_radius = 1.0600 

        !Sulfur
        ptable(16) % symbol = 'S '
        ptable(16) % symbol_number = 16
        ptable(16) % mass = 31.97207 
        ptable(16) % covalent_radius = 1.0200 

        ! Chlorine
        ptable(17) % symbol = 'Cl'
        ptable(17) % symbol_number = 17
        ptable(17) % mass = 34.96885 
        ptable(17) % covalent_radius = 0.9900 

        ! Argon
        ptable(18) % symbol = 'Ar'
        ptable(18) % symbol_number = 18
        ptable(18) % mass = 39.94800 
        ptable(18) % covalent_radius = 0.9800 

        !Potassium
        ptable(19) % symbol = 'K '
        ptable(19) % symbol_number = 19
        ptable(19) % mass = 38.96371 
        ptable(19) % covalent_radius = 2.0300 

        ! Calcium
        ptable(20) % symbol = 'Ca'
        ptable(20) % symbol_number = 20
        ptable(20) % mass = 39.96259 
        ptable(20) % covalent_radius = 1.7400 

        ! Scandium
        ptable(21) % symbol = 'Sc'
        ptable(21) % symbol_number = 21
        ptable(21) % mass = 44.95592 
        ptable(21) % covalent_radius = 1.4400 

        ! Titanium
        ptable(22) % symbol = 'Ti'
        ptable(22) % symbol_number = 22
        ptable(22) % mass = 48.00000 
        ptable(22) % covalent_radius = 1.3200 

        ! Vanadium 
        ptable(23) % symbol = 'V '
        ptable(23) % symbol_number = 23
        ptable(23) % mass = 50.94400 
        ptable(23) % covalent_radius = 1.2200 

        ! Chromium
        ptable(24) % symbol = 'Cr'
        ptable(24) % symbol_number = 24
        ptable(24) % mass = 51.94050 
        ptable(24) % covalent_radius = 1.1800 
        
        ! Manganese
        ptable(25) % symbol = 'Mn'
        ptable(25) % symbol_number = 25
        ptable(25) % mass = 54.93810 
        ptable(25) % covalent_radius = 1.1700 

        ! Iron
        ptable(26) % symbol = 'Fe'
        ptable(26) % symbol_number = 26
        ptable(26) % mass = 55.93490 
        ptable(26) % covalent_radius = 1.1700 

        ! Cobalt
        ptable(27) % symbol = 'Co'
        ptable(27) % symbol_number = 27
        ptable(27) % mass = 58.93320 
        ptable(27) % covalent_radius = 1.1600 
 
        ! Nickel
        ptable(28) % symbol = 'Ni'
        ptable(28) % symbol_number = 28
        ptable(28) % mass = 57.93530 
        ptable(28) % covalent_radius = 1.1500 

        ! Copper
        ptable(29) % symbol = 'Cu'
        ptable(29) % symbol_number = 29
        ptable(29) % mass = 63.20000 
        ptable(29) % covalent_radius = 1.1700 

        ! Zinc
        ptable(30) % symbol = 'Zn'
        ptable(30) % symbol_number = 30
        ptable(30) % mass = 63.92910 
        ptable(30) % covalent_radius = 1.2500 

        ! Gallium
        ptable(31) % symbol = 'Ga'
        ptable(31) % symbol_number = 31
        ptable(31) % mass = 68.92570 
        ptable(31) % covalent_radius = 1.2600 

        ! Germanium
        ptable(32) % symbol = 'Ge'
        ptable(32) % symbol_number = 32
        ptable(32) % mass = 73.92140 
        ptable(32) % covalent_radius = 1.2200 

        ! Arsenic
        ptable(33) % symbol = 'As'
        ptable(33) % symbol_number = 33
        ptable(33) % mass = 74.92160 
        ptable(33) % covalent_radius = 1.2000 

        ! Selenium
        ptable(34) % symbol = 'Se'
        ptable(34) % symbol_number = 34
        ptable(34) % mass = 79.91650 
        ptable(34) % covalent_radius = 1.1600 

        ! Bromine
        ptable(35) % symbol = 'Br'
        ptable(35) % symbol_number = 35
        ptable(35) % mass = 78.91830 
        ptable(35) % covalent_radius = 1.1400 

        ! Krypton
        ptable(36) % symbol = 'Kr'
        ptable(36) % symbol_number = 36
        ptable(36) % mass = 84.00000 
        ptable(36) % covalent_radius = 1.1200 

        ! Rubidium
        ptable(37) % symbol = 'Rb'
        ptable(37) % symbol_number = 37
        ptable(37) % mass = 84.91170 
        ptable(37) % covalent_radius = 2.1600 
        
        ! Strontium
        ptable(38) % symbol = 'Sr'
        ptable(38) % symbol_number = 38
        ptable(38) % mass = 87.90560 
        ptable(38) % covalent_radius = 1.9100 

        ! Yttrium
        ptable(39) % symbol = 'Y '
        ptable(39) % symbol_number = 39
        ptable(39) % mass = 88.90590 
        ptable(39) % covalent_radius = 1.6200 
        
        ! Zirconium
        ptable(40) % symbol = 'Zr'
        ptable(40) % symbol_number = 40
        ptable(40) % mass = 89.90430 
        ptable(40) % covalent_radius = 1.4500 
        
        ! Niobium
        ptable(41) % symbol = 'Nb'
        ptable(41) % symbol_number = 41
        ptable(41) % mass = 92.90600 
        ptable(41) % covalent_radius = 1.3400 

        ! Molybdenum
        ptable(42) % symbol = 'Mo'
        ptable(41) % symbol_number = 41
        ptable(42) % mass = 97.90550 
        ptable(42) % covalent_radius = 1.3000 
        
        ! Technetium
        ptable(43) % symbol = 'Tc'
        ptable(42) % symbol_number = 42
        ptable(43) % mass = 98.90600 
        ptable(43) % covalent_radius = 1.2700 
        
        ! Ruthenium
        ptable(44) % symbol = 'Ru'
        ptable(44) % symbol_number = 44
        ptable(44) % mass = 101.90370 
        ptable(44) % covalent_radius = 1.2500 
        
        ! Rhodium
        ptable(45) % symbol = 'Rh'
        ptable(45) % symbol_number = 45
        ptable(45) % mass = 102.90480 
        ptable(45) % covalent_radius = 1.2500 

        ! Palladium
        ptable(46) % symbol = 'Pd'
        ptable(46) % symbol_number = 46
        ptable(46) % mass = 105.90320 
        ptable(46) % covalent_radius = 1.2800 

        ! Silver
        ptable(47) % symbol = 'Ag'
        ptable(47) % symbol_number = 47
        ptable(47) % mass = 106.90509 
        ptable(47) % covalent_radius = 1.3400 
        
        ! Cadmium
        ptable(48) % symbol = 'Cd'
        ptable(48) % symbol_number = 48
        ptable(48) % mass = 113.90360 
        ptable(48) % covalent_radius = 1.4800 
        
        ! Indium
        ptable(49) % symbol = 'In'
        ptable(49) % symbol_number = 49
        ptable(49) % mass = 114.90410 
        ptable(49) % covalent_radius = 1.4400 

        ! Tin
        ptable(50) % symbol = 'Sn'
        ptable(50) % symbol_number = 50
        ptable(50) % mass = 120.00000 
        ptable(50) % covalent_radius = 1.4100 
        
        ! Antimony
        ptable(51) % symbol = 'Sb'
        ptable(51) % symbol_number = 51
        ptable(51) % mass = 120.90380 
        ptable(51) % covalent_radius = 1.4000 

        ! Tellurium
        ptable(52) % symbol = 'Te'
        ptable(52) % symbol_number = 52
        ptable(52) % mass = 129.90670 
        ptable(52) % covalent_radius = 1.3600 
        
        ! Iodine
        ptable(53) % symbol = 'I '
        ptable(53) % symbol_number = 53
        ptable(53) % mass = 126.90440 
        ptable(53) % covalent_radius = 1.3300 
        
        ! Xenon
        ptable(54) % symbol = 'Xe'
        ptable(54) % symbol_number = 54
        ptable(54) % mass = 131.90420 
        ptable(54) % covalent_radius = 1.3100 

        ! Cesium
        ptable(55) % symbol = 'Cs'
        ptable(55) % symbol_number = 55
        ptable(55) % mass = 132.90510 
        ptable(55) % covalent_radius = 2.3500 
        
        ! Barium
        ptable(56) % symbol = 'Ba'
        ptable(56) % symbol_number = 56
        ptable(56) % mass = 137.90500 
        ptable(56) % covalent_radius = 1.9800 
        
        ! Lantanum
        ptable(57) % symbol = 'La'
        ptable(57) % symbol_number = 57
        ptable(57) % mass = 138.90610 
        ptable(57) % covalent_radius = 1.6900 
        
        ! Cerium
        ptable(58) % symbol = 'Ce'
        ptable(58) % symbol_number = 58
        ptable(58) % mass = 139.90530 
        ptable(58) % covalent_radius = 1.6500 
        
        ! Praseodymium
        ptable(59) % symbol = 'Pr'
        ptable(59) % symbol_number = 59
        ptable(59) % mass = 140.90740 
        ptable(59) % covalent_radius = 1.6500 

        ! Neodymium
        ptable(60) % symbol = 'Nd'
        ptable(60) % symbol_number = 60
        ptable(60) % mass = 141.90750 
        ptable(60) % covalent_radius = 1.6400 
        
        ! Promethium
        ptable(61) % symbol = 'Pm'
        ptable(61) % symbol_number = 61
        ptable(61) % mass = 144.91300 
        ptable(61) % covalent_radius = 1.6300 
        
        ! Samarium
        ptable(62) % symbol = 'Sm'
        ptable(62) % symbol_number = 62
        ptable(62) % mass = 151.91950 
        ptable(62) % covalent_radius = 1.6200 
        
        ! Europium
        ptable(63) % symbol = 'Eu'
        ptable(63) % symbol_number = 63
        ptable(63) % mass = 152.92090 
        ptable(63) % covalent_radius = 1.8500 

        ! Gadolinium
        ptable(64) % symbol = 'Gd'
        ptable(64) % symbol_number = 64
        ptable(64) % mass = 157.92410 
        ptable(64) % covalent_radius = 1.6100 

        ! Terbium
        ptable(65) % symbol = 'Tb'
        ptable(65) % symbol_number = 65
        ptable(65) % mass = 158.92500 
        ptable(65) % covalent_radius = 1.5900 

        ! Dysprosium
        ptable(66) % symbol = 'Dy'
        ptable(66) % symbol_number = 66
        ptable(66) % mass = 163.92880 
        ptable(66) % covalent_radius = 1.5900 
        
        ! Holmium
        ptable(67) % symbol = 'Ho'
        ptable(67) % symbol_number = 67
        ptable(67) % mass = 164.93000 
        ptable(67) % covalent_radius = 1.5800 
        
        ! Erbium
        ptable(68) % symbol = 'Er'
        ptable(68) % symbol_number = 68
        ptable(68) % mass = 165.93040 
        ptable(68) % covalent_radius = 1.5700 
        
        ! Thulium
        ptable(69) % symbol = 'Tm'
        ptable(69) % symbol_number = 69
        ptable(69) % mass = 168.93440 
        ptable(69) % covalent_radius = 1.5600 

        ! Ytterbium
        ptable(70) % symbol = 'Yb'
        ptable(70) % symbol_number = 70
        ptable(70) % mass = 173.93900 
        ptable(70) % covalent_radius = 1.5600 

        ! Lutetium
        ptable(71) % symbol = 'Lu'
        ptable(71) % symbol_number = 71
        ptable(71) % mass = 174.94090 
        ptable(71) % covalent_radius = 1.5600 
        
        ! Hafnium
        ptable(72) % symbol = 'Hf'
        ptable(72) % symbol_number = 72
        ptable(72) % mass = 179.94680 
        ptable(72) % covalent_radius = 1.4400 
        
        ! Tantalum
        ptable(73) % symbol = 'Ta'
        ptable(73) % symbol_number = 73
        ptable(73) % mass = 180.94800 
        ptable(73) % covalent_radius = 1.3400 
        
        ! Tungsten
        ptable(74) % symbol = 'W '
        ptable(74) % symbol_number = 74
        ptable(74) % mass = 183.95100 
        ptable(74) % covalent_radius = 1.3000 

        ! Rhenium
        ptable(75) % symbol = 'Re'
        ptable(75) % symbol_number = 75
        ptable(75) % mass = 186.95600 
        ptable(75) % covalent_radius = 1.2800 
        
        ! Osmium
        ptable(76) % symbol = 'Os'
        ptable(76) % symbol_number = 76
        ptable(76) % mass = 192.00000 
        ptable(76) % covalent_radius = 1.2600 
        
        ! Iridium
        ptable(77) % symbol = 'Ir'
        ptable(77) % symbol_number = 77
        ptable(77) % mass = 192.96330 
        ptable(77) % covalent_radius = 1.2700 
        
        ! Platinum
        ptable(78) % symbol = 'Pt'
        ptable(78) % symbol_number = 78
        ptable(78) % mass = 194.96480 
        ptable(78) % covalent_radius = 1.3000 

        ! Gold
        ptable(79) % symbol = 'Au'
        ptable(79) % symbol_number = 79
        ptable(79) % mass = 196.96660 
        ptable(79) % covalent_radius = 1.3400 

        ! Mercury
        ptable(80) % symbol = 'Hg'
        ptable(80) % symbol_number = 80
        ptable(80) % mass = 201.97060 
        ptable(80) % covalent_radius = 1.4900 

        ! Thallium
        ptable(81) % symbol = 'Tl'
        ptable(81) % symbol_number = 81
        ptable(81) % mass = 204.97450 
        ptable(81) % covalent_radius = 1.4800 

        ! Lead
        ptable(82) % symbol = 'Pb'
        ptable(82) % symbol_number = 82
        ptable(82) % mass = 207.97660 
        ptable(82) % covalent_radius = 1.4700 

        ! Bismuth
        ptable(83) % symbol = 'Bi'
        ptable(83) % symbol_number = 83
        ptable(83) % mass = 208.98040 
        ptable(83) % covalent_radius = 1.4600 

        ! Polonium
        ptable(84) % symbol = 'Po'
        ptable(84) % symbol_number = 84
        ptable(84) % mass = 209.98290 
        ptable(84) % covalent_radius = 1.4600 

        ! Astatine
        ptable(85) % symbol = 'At'
        ptable(85) % symbol_number = 85
        ptable(85) % mass = 209.98700 
        ptable(85) % covalent_radius = 1.4500 

        ! Radon
        ptable(86) % symbol = 'Rn'
        ptable(86) % symbol_number = 86
        ptable(86) % mass = 222.01750 
        ptable(86) % covalent_radius = 0.00

        END SUBROUTINE init_periodic_table
!************************************************************

END MODULE
