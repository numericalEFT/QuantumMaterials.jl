using DFTK
using Unitful
using UnitfulAtomic
using Plots

# Setup a sodium lattice, see the following link for more information
# http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/bcc/bcc.php
a = 0.429u"nm"

# each column is a lattice vector
lattice = a / 2 * [[1.0 -1.0 1.0] [1.0 1.0 -1.0] [-1.0 1.0 1.0]]
pseudo_file = joinpath(@__DIR__, "Na_PBE.upf")
Na = ElementPsp(:Na, psp=load_psp(pseudo_file))
atoms = [Na,]
positions = [[0.0, 0.0, 0.0],]

smearing = DFTK.Smearing.Gaussian()
temperature = 0.005

model = 
