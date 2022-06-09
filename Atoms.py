"""Infrastructure for storing periodic table data for Slater shielding code.
"""
principalLabels = [1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,
                   9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13,
                   14, 14, 14, 15, 15, 16]  # This is principal quantum number N. combine with magnetic labels to get "1s" "2sp", "3sp", etc

azimuthalLabels = ["s", "sp", "sp", "d", "sp", "d", "f", "sp", "d", "f", "g", "sp", "d", "f", "g", "h", "sp", "d", "f",
                   "g", "h", "i", "sp", "d", "f", "g", "h", "i", "j", "sp", "d", "f", "g", "h", "i", "j", "k", "sp",
                   "d", "f", "g", "h", "i", "j", "sp", "d", "f", "g", "h", "i", "sp", "d", "f", "g", "h", "sp", "d",
                   "f", "g", "sp", "d", "f", "sp", "d", "sp"]  # This is the magnetic quantum number l, combined into Slater groups.


class NewAtom:
    """A new atom object to better suit my needs."""
    def __init__(self, atomicNumber, name, occupancy):
        self.atomicNumber = atomicNumber
        self.name = name
        self.occupancy = occupancy
        self.principalQuantumNumberLabelList = []
        self.azimuthalQuantumNumberLabelList = []
        for i in range(len(self.occupancy)):
            self.principalQuantumNumberLabelList.append(principalLabels[i])
            self.azimuthalQuantumNumberLabelList.append(azimuthalLabels[i])


class Atom:
    """quick and dirty class for structured atomic data.
    Use with Periodic Table to get a dictionary of elements
    """

    def __init__(self, Z, label, occupancy, config="Normal", filename=None):
        self.Z = Z
        self.label = label
        self.occupancy = occupancy
        self.config = "Normal"
        self.filename = filename
        self.densintegral = None

    def __repr__(self):
        line1 = "#" + self.label + " Atom" + "\n"
        if self.config != "Normal":
            linex = "#" + self.config + "\n"  # This will always be "Normal" because it is defined explicitly instead of
                                    # being read in from the argument passed to it in __init__. Is this intentional?
                                    # Looking at posIon(), it looks like this might just be a placeholder?
        else:
            linex = ""
        line2 = "#Nuclear charge = " + str(self.Z) + "\n"
        line3 = "#Occupancy list = " + str(self.occupancy) + "\n"
        if self.filename is not None:
            line4 = "#Filename for density = " + self.filename + "\n"
        else:
            line4 = "#Uses Slater model for density"
        if self.densintegral is not None:
            line5 = "#Density Integral = " + str(self.densintegral) + "\n"
        else:
            line5 = ""
        return line1 + linex + line2 + line3 + line4 + line5


#    def initialize(self):

class Grid:
    """Holder class for grid on which to display atom"""
    name = "Exponential"
    amesh = 1.05
    hmin = 0.006256
    rmax = 10.0

    def __init__(self, amesh=None, hmin=None, rmax=None):
        if amesh is not None:
            self.amesh = amesh
        if hmin is not None:
            self.hmin = hmin
        if rmax is not None:  # again, these are all defined explicitly, so will ALWAYS be None. Is there a reason?
            self.rmax = rmax

    def __str__(self):
        str1 = self.name + " Grid:\n"
        str1 += "amesh = %f\n" % self.amesh
        str1 += "hmin  = %f\n" % self.hmin
        str1 += "rmax  = %f\n" % self.rmax
        return str1


# Create a Dictionary of Neutral Atoms.
powers = [1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
          9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15,
          16]  # this is principal quantum number N
# Why not defined through R (or leave out 5d's, 4f's)
occupancylabels = ["1s", "2sp", "3sp", "3d", "4sp", "4d", "4f", "5sp", "5d", "5f", "5g", "6sp", "6d", "6f", "6g", "6h",
                   "7sp", "7d", "7f", "7g", "7h", "7i", "8sp", "8d", "8f", "8g", "8h", "8i", "8j", "9sp", "9d", "9f",
                   "9g", "9h", "9i", "9j", "9k", "10sp", "10d", "10f", "10g", "10h", "10i", "10j", "11sp", "11d", "11f",
                   "11g", "11h", "11i", "12sp", "12d", "12f", "12g", "12h", "13sp", "13d", "13f", "13g", "14sp", "14d",
                   "14f", "15sp", "15d", "16sp"]
maxoccupancy = [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 14, 18, 22, 8, 10, 14, 18, 22, 26, 8, 10, 14, 18, 22, 26,
                30, 8, 10, 14, 18, 22, 26, 30, 34, 8, 10, 14, 18, 22, 26, 30, 8, 10, 14, 18, 22, 26, 8, 10, 14, 18, 22,
                8, 10, 14, 18, 8, 10, 14, 8, 10, 8]
periodictable = {}
# Normal
periodictable['H'] = Atom(1, "H", [1])  # 1s
periodictable['He'] = Atom(2, "He", [2])  # 1s2
periodictable['Li'] = Atom(3, "Li", [2, 1])  # [He] 2sN
periodictable['Be'] = Atom(4, "Be", [2, 2])  #
periodictable['B'] = Atom(5, "B", [2, 3])  # [He] 2s2 2pN
periodictable['C'] = Atom(6, "C", [2, 4])  #
periodictable['N'] = Atom(7, "N", [2, 5])  #
periodictable['O'] = Atom(8, "O", [2, 6])  #
periodictable['F'] = Atom(9, "F", [2, 7])  #
periodictable['Ne'] = Atom(10, "Ne", [2, 8])  # [He] 2s2 2p6
periodictable['Na'] = Atom(11, "Na", [2, 8, 1])  # [Ne] 3s1
periodictable['Mg'] = Atom(12, "Mg", [2, 8, 2])  # [Ne] 3s2
periodictable['Al'] = Atom(13, "Al", [2, 8, 3])  # [Ne] 3s2 3p1
periodictable['Si'] = Atom(14, "Si", [2, 8, 4])  # [Ne] 3s2 3p2
periodictable['P'] = Atom(15, "P", [2, 8, 5])  # [Ne] 3s2 3p3
periodictable['S'] = Atom(16, "S", [2, 8, 6])  # [Ne] 3s2 3p4
periodictable['Cl'] = Atom(17, "Cl", [2, 8, 7])  # [Ne] 3s2 3p5
periodictable['Ar'] = Atom(18, "Ar", [2, 8, 8])  # [Ne] 3s2 3p6
periodictable['K'] = Atom(19, "K", [2, 8, 8, 0, 1])  # [Ne] 3s2 3p6
periodictable['Ca'] = Atom(20, "Ca", [2, 8, 8, 0, 2])  # [Ne] 3s2 3p6
periodictable['Sc'] = Atom(21, "Sc", [2, 8, 8, 1, 2])  # [Ne] 3s2 3p6 3d1
periodictable['Ti'] = Atom(22, "Ti", [2, 8, 8, 2, 2])  # [Ne] 3s2 3p6
periodictable['V'] = Atom(23, "V", [2, 8, 8, 3, 2])  # [Ne] 3s2 3p6
periodictable['Cr'] = Atom(24, "Cr", [2, 8, 8, 5, 1])  # [Ne] 3s2 3p6
periodictable['Mn'] = Atom(25, "Mn", [2, 8, 8, 5, 2])  # [Ne] 3s2 3p6
periodictable['Fe'] = Atom(26, "Fe", [2, 8, 8, 6, 2])  # [Ar] 4s2 3d6
periodictable['Co'] = Atom(27, "Co", [2, 8, 8, 7, 2])  # [Ne] 3s2 3p6
periodictable['Ni'] = Atom(28, "Ni", [2, 8, 8, 8, 2])  # [Ne] 3s2 3p6
periodictable['Cu'] = Atom(29, "Cu", [2, 8, 8, 10, 1])  # [Ar] 4s1 3d10
periodictable['Zn'] = Atom(30, "Zn", [2, 8, 8, 10, 2])  # [Ne] 3s2 3p6
periodictable['Ga'] = Atom(31, "Ga", [2, 8, 8, 10, 3])  # [Ne] 3s2 3p6
periodictable['Ge'] = Atom(32, "Ge", [2, 8, 8, 10, 4])
periodictable['As'] = Atom(33, "As", [2, 8, 8, 10, 5])
periodictable['Se'] = Atom(34, "Se", [2, 8, 8, 10, 6])
periodictable['Br'] = Atom(35, "Br", [2, 8, 8, 10, 7])
periodictable['Kr'] = Atom(36, "Kr", [2, 8, 8, 10, 8])  # [Ar] 3d10 4s2 4p6
periodictable['Rb'] = Atom(37, "Rb", [2, 8, 8, 10, 8, 0, 0, 1])  # [Ar] 3d10 4s2 4p6
periodictable['Sr'] = Atom(38, "Sr", [2, 8, 8, 10, 8, 0, 0, 2])  # [Ar] 3d10 4s2 4p6
periodictable['Y'] = Atom(39, "Y", [2, 8, 8, 10, 8, 1, 0, 2])
periodictable['Zr'] = Atom(40, "Zr", [2, 8, 8, 10, 8, 2, 0, 2])
periodictable['Nb'] = Atom(41, "Nb", [2, 8, 8, 10, 8, 4, 0, 1])
periodictable['Mo'] = Atom(42, "Mo", [2, 8, 8, 10, 8, 5, 0, 1])
periodictable['Tc'] = Atom(43, "Tc", [2, 8, 8, 10, 8, 6, 0, 1])
periodictable['Ru'] = Atom(44, "Ru", [2, 8, 8, 10, 8, 7, 0, 1])
periodictable['Rh'] = Atom(45, "Rh", [2, 8, 8, 10, 8, 8, 0, 1])
periodictable['Pd'] = Atom(46, "Pd", [2, 8, 8, 10, 8, 10])
periodictable['Ag'] = Atom(47, "Ag", [2, 8, 8, 10, 8, 10, 0, 1])
periodictable['Cd'] = Atom(48, "Cd", [2, 8, 8, 10, 8, 10, 0, 2])
periodictable['In'] = Atom(49, "In", [2, 8, 8, 10, 8, 10, 0, 3])
periodictable['Sn'] = Atom(50, "Sn", [2, 8, 8, 10, 8, 10, 0, 4])
periodictable['Sb'] = Atom(51, "Sb", [2, 8, 8, 10, 8, 10, 0, 5])
periodictable['Te'] = Atom(52, "Te", [2, 8, 8, 10, 8, 10, 0, 6])
periodictable['I'] = Atom(53, "I", [2, 8, 8, 10, 8, 10, 0, 7])
periodictable['Xe'] = Atom(54, "Xe", [2, 8, 8, 10, 8, 10, 0, 8])  # [Kr] 4d10 5s2 5p6
periodictable['Cs'] = Atom(55, "Cs", [2, 8, 8, 10, 8, 10, 0, 8, 0, 0, 0, 1])
periodictable['Ba'] = Atom(56, "Ba", [2, 8, 8, 10, 8, 10, 0, 8, 0, 0, 0, 2])  # [Xe] 6s2
periodictable['La'] = Atom(57, "La", [2, 8, 8, 10, 8, 10, 0, 8, 1, 0, 0, 2])
periodictable['Ce'] = Atom(58, "Ce", [2, 8, 8, 10, 8, 10, 2, 8, 0, 0, 0, 2])
periodictable['Pr'] = Atom(59, "Pr", [2, 8, 8, 10, 8, 10, 3, 8, 0, 0, 0, 2])
periodictable['Nd'] = Atom(60, "Nd", [2, 8, 8, 10, 8, 10, 4, 8, 0, 0, 0, 2])
periodictable['Pm'] = Atom(61, "Pm", [2, 8, 8, 10, 8, 10, 5, 8, 0, 0, 0, 2])
periodictable['Sm'] = Atom(62, "Sm", [2, 8, 8, 10, 8, 10, 6, 8, 0, 0, 0, 2])
periodictable['Eu'] = Atom(63, "Eu", [2, 8, 8, 10, 8, 10, 7, 8, 0, 0, 0, 2])
periodictable['Gd'] = Atom(64, "Gd", [2, 8, 8, 10, 8, 10, 7, 8, 1, 0, 0, 2])
periodictable['Tb'] = Atom(65, "Tb", [2, 8, 8, 10, 8, 10, 8, 8, 1, 0, 0, 2])
periodictable['Dy'] = Atom(66, "Dy", [2, 8, 8, 10, 8, 10, 10, 8, 0, 0, 0, 2])
periodictable['Ho'] = Atom(67, "Ho", [2, 8, 8, 10, 8, 10, 11, 8, 0, 0, 0, 2])
periodictable['Er'] = Atom(68, "Er", [2, 8, 8, 10, 8, 10, 12, 8, 0, 0, 0, 2])
periodictable['Tm'] = Atom(69, "Tm", [2, 8, 8, 10, 8, 10, 13, 8, 0, 0, 0, 2])
periodictable['Yb'] = Atom(70, "Dy", [2, 8, 8, 10, 8, 10, 14, 8, 0, 0, 0, 2])
periodictable['Lu'] = Atom(71, "Dy", [2, 8, 8, 10, 8, 10, 14, 8, 1, 0, 0, 2])
periodictable['Hf'] = Atom(72, "Hf", [2, 8, 8, 10, 8, 10, 14, 8, 2, 0, 0, 2])
periodictable['Ta'] = Atom(73, "Ta", [2, 8, 8, 10, 8, 10, 14, 8, 3, 0, 0, 2])
periodictable['W'] = Atom(74, "W", [2, 8, 8, 10, 8, 10, 14, 8, 4, 0, 0, 2])
periodictable['Re'] = Atom(75, "Re", [2, 8, 8, 10, 8, 10, 14, 8, 5, 0, 0, 2])
periodictable['Os'] = Atom(76, "Os", [2, 8, 8, 10, 8, 10, 14, 8, 6, 0, 0, 2])
periodictable['Ir'] = Atom(77, "Ir", [2, 8, 8, 10, 8, 10, 14, 8, 9])
periodictable['Pt'] = Atom(78, "Pt", [2, 8, 8, 10, 8, 10, 14, 8, 9, 0, 0, 1])
periodictable['Au'] = Atom(79, "Au", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 1])
periodictable['Hg'] = Atom(80, "Hg", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 2])
periodictable['Tl'] = Atom(81, "Tl", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 3])
periodictable['Pb'] = Atom(82, "Pb", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 4])
periodictable['Bi'] = Atom(83, "Bi", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 5])
periodictable['Po'] = Atom(84, "Po", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 6])
periodictable['At'] = Atom(85, "At", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 7])
periodictable['Rn'] = Atom(86, "Rn", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 8])  # [Xe] 4f14 5d10 6s2 6p6
periodictable['Fr'] = Atom(87, "Fr", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 8, 1])
periodictable['Ra'] = Atom(88, "Ra", [2, 8, 8, 10, 8, 10, 14, 8, 10, 0, 0, 8, 2])  # [Rn] 7s2

periodictable['H2p'] = Atom(1, "H2p", [2, 3])
periodictable['H3d'] = Atom(1, "H3d", [0, 0, 0, 1])
periodictable['Uuo'] = Atom(118, "Uuo", [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 0, 8, 10, 0, 0, 0, 8])

# uber atoms
periodictable['U20'] = Atom(120, "U20", [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 0, 8, 10, 0, 0, 0, 8, 0, 0, 0, 0, 0, 2])
periodictable['U38'] = Atom(138, "U38", [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 0, 0, 0, 8, 0, 0, 0, 0, 0, 2])
periodictable['U68'] = Atom(168, "U68", [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 14, 0, 0, 8, 10, 0, 0, 0, 0, 8])
periodictable['U70'] = Atom(170, "U70",
                            [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 14, 0, 0, 8, 10, 0, 0, 0, 0, 8, 0, 0, 0, 0,
                             0, 0, 2])
periodictable['T62'] = Atom(362, "T62",
                            [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 14, 18, 22, 8, 10, 14, 18, 22, 0, 8, 10, 14,
                             18, 0, 0, 0, 8, 10, 14, 0, 0, 0, 0, 0, 8, 10, 0, 0, 0, 0, 0, 8])
periodictable['N76'] = Atom(976, "N76",
                            [2, 8, 8, 10, 8, 10, 14, 8, 10, 14, 18, 8, 10, 14, 18, 22, 8, 10, 14, 18, 22, 26, 8, 10, 14,
                             18, 22, 26, 30, 8, 10, 14, 18, 22, 26, 30, 34, 8, 10, 14, 18, 22, 26, 30, 8, 10, 14, 18,
                             22, 26, 8, 10, 14, 18, 22, 8, 10, 14, 18, 8, 10, 14, 8, 10, 8])

# Finite Jellium drop system from Cyrus Umrigar
periodictable['Jellium3'] = Atom(18, "Jellium3", [2, 6, 0, 10])  # Jellium -- 1s2 2p6 3d10 -- psp!!
periodictable['Jellium3Smooth'] = Atom(18, "Jellium3Smooth", [2, 6, 0, 10])  # Jellium 1s2 2p6 3d10

periodictable_ions = {}
# Negative ions:  ##### Add Fluorides
periodictable_ions['H -'] = Atom(1, "H -", [2])  # 1s
periodictable_ions['He -'] = Atom(2, "He -", [2, 1])  # 1s2
periodictable_ions['Li -'] = Atom(3, "Li -", [2, 2])  # [He] 2sN
periodictable_ions['Be -'] = Atom(4, "Be -", [2, 3])  #
periodictable_ions['B -'] = Atom(5, "B -", [2, 4])  # [He] 2s2 2pN
periodictable_ions['C -'] = Atom(6, "C -", [2, 5])  #
periodictable_ions['N -'] = Atom(7, "N -", [2, 6])  #
periodictable_ions['O -'] = Atom(8, "O -", [2, 7])  #
periodictable_ions['Fl -'] = Atom(9, "Fl -", [2, 8])  #
periodictable_ions['Ne -'] = Atom(10, "Ne -", [2, 8, 1])  # [He] 2s2 2p6
periodictable_ions['Na -'] = Atom(11, "Na -", [2, 8, 2])  # [Ne] 3s1
periodictable_ions['Mg -'] = Atom(12, "Mg -", [2, 8, 3])  # [Ne] 3s2
periodictable_ions['Al -'] = Atom(13, "Al -", [2, 8, 4])  # [Ne] 3s2 3p1
periodictable_ions['Si -'] = Atom(14, "Si -", [2, 8, 5])  # [Ne] 3s2 3p2
periodictable_ions['P -'] = Atom(15, "P -", [2, 8, 6])  # [Ne] 3s2 3p3
periodictable_ions['S -'] = Atom(16, "S -", [2, 8, 7])  # [Ne] 3s2 3p4
periodictable_ions['Cl -'] = Atom(17, "Cl -", [2, 8, 8])  # [Ne] 3s2 3p5
periodictable_ions['Ar -'] = Atom(18, "Ar -", [2, 8, 8, 1])  # [Ne] 3s2 3p6
periodictable_ions['Fe -'] = Atom(26, "Fe -", [2, 8, 8, 6, 3])  # [Ar] 4s2 3d6
periodictable_ions['Kr -'] = Atom(36, "Kr -", [2, 8, 8, 10, 8, 1])  # [Ar] 3d10 4s2 4p6
periodictable_ions['Xe -'] = Atom(54, "Xe -", [2, 8, 8, 10, 8, 10, 0, 8, 1])  # [Kr] 4d10 5s2 5p6
# periodictable['Rn -'] = Atom(86, "Rn -",  [2,8,8,10,8,10,14,8,10,8,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) # [Xe] 4f14 5d10 6s2 6p6
# Positive ion   ###### Add simple metals
periodictable_ions['H +'] = Atom(1, "H +", [0])  # 1s
periodictable_ions['He +'] = Atom(2, "He +", [1])  # 1s2
periodictable_ions['Li +'] = Atom(3, "Li +", [2])  # [He] 2sN
periodictable_ions['Be +'] = Atom(4, "Be +", [2, 1])  #
periodictable_ions['B +'] = Atom(5, "B +", [2, 2])  # [He] 2s2 2pN
periodictable_ions['C +'] = Atom(6, "C +", [2, 3])  #
periodictable_ions['N +'] = Atom(7, "N +", [2, 4])  #
periodictable_ions['O +'] = Atom(8, "O +", [2, 5])  #
periodictable_ions['Fl +'] = Atom(9, "Fl +", [2, 6])  #
periodictable_ions['Ne +'] = Atom(10, "Ne +", [2, 7])  # [He] 2s2 2p6
periodictable_ions['Na +'] = Atom(11, "Na +", [2, 8])  # [Ne] 3s1
periodictable_ions['Mg +'] = Atom(12, "Mg +", [2, 8, 1])  # [Ne] 3s2
periodictable_ions['Al +'] = Atom(13, "Al +", [2, 8, 2])  # [Ne] 3s2 3p1
periodictable_ions['Si +'] = Atom(14, "Si +", [2, 8, 3])  # [Ne] 3s2 3p2
periodictable_ions['P +'] = Atom(15, "P +", [2, 8, 4])  # [Ne] 3s2 3p3
periodictable_ions['S +'] = Atom(16, "S +", [2, 8, 5])  # [Ne] 3s2 3p4
periodictable_ions['Cl +'] = Atom(17, "Cl +", [2, 8, 6])  # [Ne] 3s2 3p5
periodictable_ions['Ar +'] = Atom(18, "Ar +", [2, 8, 7])  # [Ne] 3s2 3p6
periodictable_ions['K +'] = Atom(19, "K +", [2, 8, 8])  # [Ne] 3s2 3p6
periodictable_ions['Fe +'] = Atom(26, "Fe +", [2, 8, 8, 6, 1])  # [Ar] 4s2 3d6
periodictable_ions['Kr +'] = Atom(36, "Kr +", [2, 8, 8, 10, 7])  # [Ar] 3d10 4s2 4p6
periodictable_ions['Rb +'] = Atom(37, "Rb +", [2, 8, 8, 10, 8])  # [Ar] 3d10 4s2 4p6
periodictable_ions['Xe +'] = Atom(54, "Xe +", [2, 8, 8, 10, 8, 10, 0, 7])  # [Kr] 4d10 5s2 5p6
periodictable_ions['Cs +'] = Atom(55, "Cs +", [2, 8, 8, 10, 8, 10, 0, 8])  # [Kr] 4d10 5s2 5p6
periodictable_ions['Rn +'] = Atom(86, "Rn +", [2, 8, 8, 10, 8, 10, 14, 8, 10, 7])  # [Xe] 4f14 5d10 6s2 6p6


# atom configuration definitions
def PosIon(atom):
    label = atom.label + " +"
    atom.occupancy = periodictable_ions[label].occupancy
    atom.config = "Positive Ion"
    return atom


def NegIon(atom):
    label = atom.label + " -"
    atom.occupancy = periodictable_ions[label].occupancy
    atom.config = "Negative Ion"
    return atom


def InnerShell(atom):
    atom.config = "1s^2 Shell"
    atom.occupancy = periodictable["He"].occupancy
    return atom


def Valence(atom):
    # Keep last occupied shell only in normal configuration
    have_valence = False
    for i in range(len(atom.occupancy), 0, -1):
        if atom.occupancy[i - 1] != 0:
            if have_valence == False:
                have_valence = True
            else:
                atom.occupancy[i - 1] = 0
    atom.config = "Valence Shell"
    # FIX Shielding!!!
    return atom


def Normal(atom):  # This literally does nothing? (Daniel 2-6-2022)
    return atom


configTable = {"Normal": Normal, "PosIon": PosIon, "NegIon": NegIon, "InnerShell": InnerShell, "Valence": Valence}


def setAtomConfig(atom, atomconfig):
    return configTable[atomconfig](atom)
