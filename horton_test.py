

from horton import *

import numpy as np

#([[1.48123726, -0.93019123, -0.        ], [-0. ,         0.11720081, -0.        ], [-1.48123726, -0.93019123, -0.        ]])
coord= np.array([[ 0.783837, -0.492236, 0], [0.,  0.06202 , 0.],
       [-0.783837, -0.492236, 0.]])
numb= np.squeeze(np.array([[1, 8, 1]]))

mol= Molecule(coordinates = coord,  numbers = numb)

obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')

# Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)

# Create alpha orbitals
exp_alpha = lf.create_expansion()

# Initial guess
guess_core_hamiltonian(olp, kin, na, exp_alpha)

# Setup integration grids with default settings
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

# Construction of Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
libxc_term = RLibXCHybridGGA('xc_o3lyp')
terms = [
    ROneBodyTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RGridGroup(obasis, grid, [
        libxc_term,
    ]),
    RExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
    ROneBodyTerm(na, 'ne'),
]
ham = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons)
occ_model = AufbauOccModel(5)
occ_model.assign(exp_alpha)

# Optimal damping SCF cycle
# - construct the initial density matrix
dm_alpha = exp_alpha.to_dm()
# - scf solver
scf_solver = ODASCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, dm_alpha)
