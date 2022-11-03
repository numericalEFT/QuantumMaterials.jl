using DFTK
using Unitful
using UnitfulAtomic
using Plots

# Setup a sodium lattice, see the following link for more information
# http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/bcc/bcc.php
a = 0.429u"nm"

# each column is a lattice vector
lattice = a / 2 * [[1.0 -1.0 1.0]; [1.0 1.0 -1.0]; [-1.0 1.0 1.0]]
pseudo_file = joinpath(@__DIR__, "Na_PBE.upf")
Na = ElementPsp(:Na, psp=load_psp(pseudo_file))
atoms = [Na,]
positions = [[0.0, 0.0, 0.0],]

smearing = DFTK.Smearing.Gaussian()
temperature = 0.01u"Ry"

model = model_PBE(lattice, atoms, positions; temperature, smearing)

Ecut = 48.0u"Ry"
kgrid = [16, 16, 16]

basis = PlaneWaveBasis(model; Ecut, kgrid=kgrid)

scfres = self_consistent_field(basis, tol=1e-8);

# plot
plot_dos(scfres) # in hartree unit
plot_bandstructure(scfres, kline_density=10, unit=u"eV", ρ=scfres.ρ)

# calculate band width
