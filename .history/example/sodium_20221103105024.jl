using DFTK
using Unitful
using UnitfulAtomic
using Plots

# Setup a sodium lattice, see the following link for more information
# http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/bcc/bcc.php
a = 0.429u"nm"

# each column is a lattice vector
lattice = a / 2 * [[1.0 -1.0 1.0] [1.0 1.0 -1.0] [-1.0 1.0 1.0]]

Na = ElementPsp(:Na, psp=load_psp("Na_PBE.upf"))