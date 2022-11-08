MODULE COVALENT_RADII
        IMPLICIT NONE
        PUBLIC :: atom
!************************************************************        
        TYPE atom
                CHARACTER (LEN = 2) :: symbol
                REAL :: mass
                REAL :: covalent_radius
        END TYPE atom

        INTEGER, PARAMETER :: n_elem = 86
        TYPE (element) :: ptable (0:nelem)
        CONTAINS

!************************************************************
        SUBROUTINE init_periodic_table()

        ! Hydrogen
        ptable(1) % symbol = 'H '
        ptable(1) % mass = 1.007825_dp
        ptable(1) % covalent_radius = 0.32_dp
        
        ! Helium
        ptable(2) % symbol = 'He'
        ptable(2) % mass = 4.00260_dp
        ptable(2) % covalent_radius = 0.9300_dp
        
        ! Lithium
        ptable(3) % symbol = 'Li'
        ptable(3) % mass = 7.01600_dp
        ptable(3) % covalent_radius = 1.2300_dp
        
        ! Beryllium
        ptable(4) % symbol = 'Be'
        ptable(4) % mass = 9.01218_dp
        ptable(4) % covalent_radius = 0.9000_dp
        
        ! Boron
        ptable(5) % symbol = 'B '
        ptable(5) % mass = 11.00931_dp
        ptable(5) % covalent_radius = 0.8200_dp 
        
        ! Carbon
        ptable(6) % symbol = 'C '
        ptable(6) % mass = 12.0000_dp
        ptable(6) % covalent_radius = 0.7700_dp
        
        ! Nitrogen
        ptable(7) % symbol = 'N '
        ptable(7) % mass = 14.00307_dp
        ptable(7) % covalent_radius = 0.7500_dp
        
        ! Oxygen
        ptable(8) % symbol = 'O '
        ptable(8) % mass = 15.99491_dp
        ptable(8) % covalent_radius = 0.7300_dp
        
        ! Fluorine
        ptable(9) % symbol = 'F '
        ptable(9) % mass = 18.99840_dp
        ptable(9) % covalent_radius = 0.7200_dp
        
        ! Neon
        ptable(10) % symbol = 'Ne'
        ptable(10) % mass = 19.99244_dp
        ptable(10) % covalent_radius = 0.7100_dp
        
        ! Sodium
        ptable(11) % symbol = 'Na'
        ptable(11) % mass = 22.9898_dp
        ptable(11) % covalent_radius = 1.5400_dp

        ! Magnesium
        ptable(12) % symbol = 'Mg'
        ptable(12) % mass = 23.98504_dp
        ptable(12) % covalent_radius = 1.3600_dp

        ! Aluminium
        ptable(13) % symbol = 'Al'
        ptable(13) % mass = 26.98153_dp
        ptable(13) % covalent_radius = 1.1800_dp

        ! Silicon
        ptable(14) % symbol = 'Si'
        ptable(14) % mass = 27.97693_dp
        ptable(14) % covalent_radius = 1.1100_dp

        ! Phosphorus
        ptable(15) % symbol = 'P '
        ptable(15) % mass = 30.97376_dp
        ptable(15) % covalent_radius = 1.0600_dp

        !Sulfur
        ptable(16) % symbol = 'S '
        ptable(16) % mass = 31.97207_dp
        ptable(16) % covalent_radius = 1.0200_dp

        ! Chlorine
        ptable(17) % symbol = 'Cl'
        ptable(17) % mass = 34.96885_dp
        ptable(17) % covalent_radius = 0.9900_dp

        ! Argon
        ptable(18) % symbol = 'Ar'
        ptable(18) % mass = 39.94800_dp
        ptable(18) % covalent_radius = 0.9800_dp

        !Potassium
        ptable(19) % symbol = 'K '
        ptable(19) % mass = 38.96371_dp
        ptable(19) % covalent_radius = 2.0300_dp

        ! Calcium
        ptable(20) % symbol = 'Ca'
        ptable(20) % mass = 39.96259_dp
        ptable(20) % covalent_radius = 1.7400_dp

        ! Scandium
        ptable(21) % symbol = 'Sc'
        ptable(21) % mass = 44.95592_dp
        ptable(21) % covalent_radius = 1.4400_dp

        ! Titanium
        ptable(22) % symbol = 'Ti'
        ptable(22) % mass = 48.00000_dp
        ptable(22) % covalent_radius = 1.3200_dp

        ! Vanadium 
        ptable(23) % symbol = 'V '
        ptable(23) % mass = 50.94400_dp
        ptable(23) % covalent_radius = 1.2200_dp

        ! Chromium
        ptable(24) % symbol = 'Cr'
        ptable(24) % mass = 51.94050_dp
        ptable(24) % covalent_radius = 1.1800_dp
        
        ! Manganese
        ptable(25) % symbol = 'Mn'
        ptable(25) % mass = 54.93810_dp
        ptable(25) % covalent_radius = 1.1700_dp

        ! Iron
        ptable(26) % symbol = 'Fe'
        ptable(26) % mass = 55.93490_dp
        ptable(26) % covalent_radius = 1.1700_dp

        ! Cobalt
        ptable(27) % symbol = 'Co'
        ptable(27) % mass = 58.93320_dp
        ptable(27) % covalent_radius = 1.1600_dp
 
        ! Nickel
        ptable(28) % symbol = 'Ni'
        ptable(28) % mass = 57.93530_dp
        ptable(28) % covalent_radius = 1.1500_dp

        ! Copper
        ptable(29) % symbol = 'Cu'
        ptable(29) % mass = 63.20000_dp
        ptable(29) % covalent_radius = 1.1700_dp

        ! Zinc
        ptable(30) % symbol = 'Zn'
        ptable(30) % mass = 63.92910_dp
        ptable(30) % covalent_radius = 1.2500_dp

        ! Gallium
        ptable(31) % symbol = 'Ga'
        ptable(31) % mass = 68.92570_dp
        ptable(31) % covalent_radius = 1.2600_dp

        ! Germanium
        ptable(32) % symbol = 'Ge'
        ptable(32) % mass = 73.92140_dp
        ptable(32) % covalent_radius = 1.2200_dp

        ! Arsenic
        ptable(33) % symbol = 'As'
        ptable(33) % mass = 74.92160_dp
        ptable(33) % covalent_radius = 1.2000_dp

        ! Selenium
        ptable(34) % symbol = 'Se'
        ptable(34) % mass = 79.91650_dp
        ptable(34) % covalent_radius = 1.1600_dp

        ! Bromine
        ptable(35) % symbol = 'Br'
        ptable(35) % mass = 78.91830_dp
        ptable(35) % covalent_radius = 1.1400_dp

        ! Krypton
        ptable(36) % symbol = 'Kr'
        ptable(36) % mass = 84.00000_dp
        ptable(36) % covalent_radius = 1.1200_dp

        ! Rubidium
        ptable(37) % symbol = 'Rb'
        ptable(37) % mass = 84.91170_dp
        ptable(37) % covalent_radius = 2.1600_dp
        
        ! Strontium
        ptable(38) % symbol = 'Sr'
        ptable(38) % mass = 87.90560_dp
        ptable(38) % covalent_radius = 1.9100_dp

        ! Yttrium
        ptable(39) % symbol = 'Y '
        ptable(39) % mass = 88.90590_dp
        ptable(39) % covalent_radius = 1.6200_dp
        
        ! Zirconium
        ptable(40) % symbol = 'Zr'
        ptable(40) % mass = 89.90430_dp
        ptable(40) % covalent_radius = 1.4500_dp
        
        ! Niobium
        ptable(41) % symbol = 'Nb'
        ptable(41) % mass = 92.90600_dp
        ptable(41) % covalent_radius = 1.3400_dp

        ! Molybdenum
        ptable(42) % symbol = 'Mo'
        ptable(42) % mass = 97.90550_dp
        ptable(42) % covalent_radius = 1.3000_dp
        
        ! Technetium
        ptable(43) % symbol = 'Tc'
        ptable(43) % mass = 98.90600_dp
        ptable(43) % covalent_radius = 1.2700_dp
        
        ! Ruthenium
        ptable(44) % symbol = 'Ru'
        ptable(44) % mass = 101.90370_dp
        ptable(44) % covalent_radius = 1.2500_dp
        
        ! Rhodium
        ptable(45) % symbol = 'Rh'
        ptable(45) % mass = 102.90480_dp
        ptable(45) % covalent_radius = 1.2500_dp

        ! Palladium
        ptable(46) % symbol = 'Pd'
        ptable(46) % mass = 105.90320_dp
        ptable(46) % covalent_radius = 1.2800_dp

        ! Silver
        ptable(47) % symbol = 'Ag'
        ptable(47) % mass = 106.90509_dp
        ptable(47) % covalent_radius = 1.3400_dp
        
        ! Cadmium
        ptable(48) % symbol = 'Cd'
        ptable(48) % mass = 113.90360_dp
        ptable(48) % covalent_radius = 1.4800_dp
        
        ! Indium
        ptable(49) % symbol = 'In'
        ptable(49) % mass = 114.90410_dp
        ptable(49) % covalent_radius = 1.4400_dp

        ! Tin
        ptable(50) % symbol = 'Sn'
        ptable(50) % mass = 120.00000_dp
        ptable(50) % covalent_radius = 1.4100_dp
        
        ! Antimony
        ptable(51) % symbol = 'Sb'
        ptable(51) % mass = 120.90380_dp
        ptable(51) % covalent_radius = 1.4000_dp

        ! Tellurium
        ptable(52) % symbol = 'Te'
        ptable(52) % mass = 129.90670_dp
        ptable(52) % covalent_radius = 1.3600_dp
        
        ! Iodine
        ptable(53) % symbol = 'I '
        ptable(53) % mass = 126.90440_dp
        ptable(53) % covalent_radius = 1.3300_dp
        
        ! Xenon
        ptable(54) % symbol = 'Xe'
        ptable(54) % mass = 131.90420_dp
        ptable(54) % covalent_radius = 1.3100_dp

        ! Cesium
        ptable(55) % symbol = 'Cs'
        ptable(55) % mass = 132.90510_dp
        ptable(55) % covalent_radius = 2.3500_dp
        
        ! Barium
        ptable(56) % symbol = 'Ba'
        ptable(56) % mass = 137.90500_dp
        ptable(56) % covalent_radius = 1.9800_dp
        
        ! Lantanum
        ptable(57) % symbol = 'La'
        ptable(57) % mass = 138.90610_dp
        ptable(57) % covalent_radius = 1.6900_dp
        
        ! Cerium
        ptable(58) % symbol = 'Ce'
        ptable(58) % mass = 139.90530_dp
        ptable(58) % covalent_radius = 1.6500_dp
        
        ! Praseodymium
        ptable(59) % symbol = 'Pr'
        ptable(59) % mass = 140.90740_dp
        ptable(59) % covalent_radius = 1.6500_dp

        ! Neodymium
        ptable(60) % symbol = 'Nd'
        ptable(60) % mass = 141.90750_dp
        ptable(60) % covalent_radius = 1.6400_dp
        
        ! Promethium
        ptable(61) % symbol = 'Pm'
        ptable(61) % mass = 144.91300_dp
        ptable(61) % covalent_radius = 1.6300_dp
        
        ! Samarium
        ptable(62) % symbol = 'Sm'
        ptable(62) % mass = 151.91950_dp
        ptable(62) % covalent_radius = 1.6200_dp
        
        ! Europium
        ptable(63) % symbol = 'Eu'
        ptable(63) % mass = 152.92090_dp
        ptable(63) % covalent_radius = 1.8500_dp

        ! Gadolinium
        ptable(64) % symbol = 'Gd'
        ptable(64) % mass = 157.92410_dp
        ptable(64) % covalent_radius = 1.6100_dp

        ! Terbium
        ptable(65) % symbol = 'Tb'
        ptable(65) % mass = 158.92500_dp
        ptable(65) % covalent_radius = 1.5900_dp

        ! Dysprosium
        ptable(66) % symbol = 'Dy'
        ptable(66) % mass = 163.92880_dp
        ptable(66) % covalent_radius = 1.5900_dp
        
        ! Holmium
        ptable(67) % symbol = 'Ho'
        ptable(67) % mass = 164.93000_dp
        ptable(67) % covalent_radius = 1.5800_dp
        
        ! Erbium
        ptable(68) % symbol = 'Er'
        ptable(68) % mass = 165.93040_dp
        ptable(68) % covalent_radius = 1.5700_dp
        
        ! Thulium
        ptable(69) % symbol = 'Tm'
        ptable(69) % mass = 168.93440_dp
        ptable(69) % covalent_radius = 1.5600_dp

        ! Ytterbium
        ptable(70) % symbol = 'Yb'
        ptable(70) % mass = 173.93900_dp
        ptable(70) % covalent_radius = 1.5600_dp

        ! Lutetium
        ptable(71) % symbol = 'Lu'
        ptable(71) % mass = 174.94090_dp
        ptable(71) % covalent_radius = 1.5600_dp
        
        ! Hafnium
        ptable(72) % symbol = 'Hf'
        ptable(72) % mass = 179.94680_dp
        ptable(72) % covalent_radius = 1.4400_dp
        
        ! Tantalum
        ptable(73) % symbol = 'Ta'
        ptable(73) % mass = 180.94800_dp
        ptable(73) % covalent_radius = 1.3400_dp
        
        ! Tungsten
        ptable(74) % symbol = 'W '
        ptable(74) % mass = 183.95100_dp
        ptable(74) % covalent_radius = 1.3000_dp

        ! Rhenium
        ptable(75) % symbol = 'Re'
        ptable(75) % mass = 186.95600_dp
        ptable(75) % covalent_radius = 1.2800_dp
        
        ! Osmium
        ptable(76) % symbol = 'Os'
        ptable(76) % mass = 192.00000_dp
        ptable(76) % covalent_radius = 1.2600_dp
        
        ! Iridium
        ptable(77) % symbol = 'Ir'
        ptable(77) % mass = 192.96330_dp
        ptable(77) % covalent_radius = 1.2700_dp
        
        ! Platinum
        ptable(78) % symbol = 'Pt'
        ptable(78) % mass = 194.96480_dp
        ptable(78) % covalent_radius = 1.3000_dp

        ! Gold
        ptable(79) % symbol = 'Au'
        ptable(79) % mass = 196.96660_dp
        ptable(79) % covalent_radius = 1.3400_dp

        ! Mercury
        ptable(80) % symbol = 'Hg'
        ptable(80) % mass = 201.97060_dp
        ptable(80) % covalent_radius = 1.4900_dp

        ! Thallium
        ptable(81) % symbol = 'Tl'
        ptable(81) % mass = 204.97450_dp
        ptable(81) % covalent_radius = 1.4800_dp

        ! Lead
        ptable(82) % symbol = 'Pb'
        ptable(82) % mass = 207.97660_dp
        ptable(82) % covalent_radius = 1.4700_dp

        ! Bismuth
        ptable(83) % symbol = 'Bi'
        ptable(83) % mass = 208.98040_dp
        ptable(83) % covalent_radius = 1.4600_dp

        ! Polonium
        ptable(84) % symbol = 'Po'
        ptable(84) % mass = 209.98290_dp
        ptable(84) % covalent_radius = 1.4600_dp

        ! Astatine
        ptable(85) % symbol = 'At'
        ptable(85) % mass = 209.98700_dp
        ptable(85) % covalent_radius = 1.4500_dp

        ! Radon
        ptable(86) % symbol = 'Rn'
        ptable(86) % mass = 222.01750_dp
        ptable(86) % covalent_radius = z

        END SUBROUTINE init_periodic_table
!************************************************************

END MODULE
