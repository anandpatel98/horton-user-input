import user_input_parsing as uip
import atom_attributes as aa
import default_selection as ds
import input_conversion as ic
from horton import *


import numpy as np

def exchange_term(restriction = ic.mol.restriction):
    if restriction == 'restricted'and ic.hybrid:
        return RExchangeTerm(er, 'x_hf', ic.hybrid.get_exx_fraction())
    if restriction == 'unrestricted'and ic.hybrid:
        return UExchangeTerm(er, 'x_hf', ic.hybrid.get_exx_fraction())



obasis = get_gobasis(ic.mol.coordinates, ic.mol.numbers, ic.mol.basis)
lf= DenseLinalgFactory(obasis.nbasis)
    
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(ic.mol.coordinates, ic.mol.pseudo_numbers, lf) # mandatory pseudo_numbers??
er = obasis.compute_electron_repulsion(lf)

if ds.user_input_defaults["calculation"] == "hf":
    
    if ds.user_input_defaults['restriction'] == 'restricted':
        exp_alpha = lf.create_expansion()
        guess_core_hamiltonian(olp, kin, na, exp_alpha)
        external = {'nn': compute_nucnuc(ic.mol.coordinates, ic.mol.pseudo_numbers)}
        terms = [
            ROneBodyTerm(kin, 'kin'),
            RDirectTerm(er, 'hartree'),
            RExchangeTerm(er, 'x_hf'),
            ROneBodyTerm(na, 'ne'),
        ]
        ham = REffHam(terms, external)
        occ_model = AufbauOccModel(ic.nalpha)
        occ_model.assign(exp_alpha) # Taken from hfs run file. Ask if this is okay
        dm_alpha = exp_alpha.to_dm()
        scf_solver= ic.scf_solver
        if ic.expansion:
            scf_solver(ham, lf, olp, occ_model, exp_alpha)
        else:
            scf_solver(ham, lf, olp, occ_model, dm_alpha)
            
    if ds.user_input_defaults['restriction'] == 'unrestricted':
        exp_beta = lf.create_expansion()
        exp_alpha = lf.create_expansion()
        guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)
        external = {'nn': compute_nucnuc(ic.mol.coordinates, ic.mol.pseudo_numbers)}
        terms = [
            UOneBodyTerm(kin, 'kin'),
            UDirectTerm(er, 'hartree'),
            UExchangeTerm(er, 'x_hf'),
            UOneBodyTerm(na, 'ne'),
        ]
        ham = UEffHam(terms, external)
        occ_model = AufbauOccModel(ic.nalpha,ic.nbeta)
        occ_model.assign(exp_alpha, exp_beta) # Taken from hfs run file. Ask if this is okay
        dm_alpha = exp_alpha.to_dm()
        dm_beta= exp_beta.to_dm()
        scf_solver= ic.scf_solver
        if ic.expansion:
            scf_solver(ham, lf, olp, occ_model, exp_alpha, exp_beta)
        else:
            scf_solver(ham, lf, olp, occ_model, dm_alpha, dm_beta)
    

    
if ds.user_input_defaults["calculation"] == "dft":

    if ds.user_input_defaults['restriction'] == 'restricted':
        exp_alpha = lf.create_expansion()
        guess_core_hamiltonian(olp, kin, na, exp_alpha)
        grid = ic.grid
        
        external = {'nn': compute_nucnuc(ic.mol.coordinates, ic.mol.pseudo_numbers)}
        terms = [
            ROneBodyTerm(kin, 'kin'),
            RDirectTerm(er, 'hartree'),
            RGridGroup(obasis, grid, ic.lib),
            ROneBodyTerm(na, 'ne'),
        ]
        if exchange_term() != None:
            terms.append(exchange_term())
        ham = REffHam(terms, external)
        
        occ_model = AufbauOccModel(ic.nalpha)
        occ_model.assign(exp_alpha) # Taken from hfs run file. Ask if this is okay
        dm_alpha = exp_alpha.to_dm()
        
        scf_solver = ic.scf_solver
        if ic.expansion:
            scf_solver(ham, lf, olp, occ_model, exp_alpha)
        else:
            scf_solver(ham, lf, olp, occ_model, dm_alpha)
            
    if ds.user_input_defaults['restriction'] == 'unrestricted':
        exp_beta = lf.create_expansion()
        exp_alpha = lf.create_expansion()
        guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)
        grid = ic.grid  
        
        external = {'nn': compute_nucnuc(ic.mol.coordinates, ic.mol.pseudo_numbers)}

        terms = [
            UOneBodyTerm(kin, 'kin'),
            UDirectTerm(er, 'hartree'),
            UGridGroup(obasis, grid, ic.lib),
            UOneBodyTerm(na, 'ne'),
        ]
        if exchange_term()  != None:
            terms.append(exchange_term())
        ham = UEffHam(terms, external)

        occ_model = AufbauOccModel(ic.nalpha,ic.nbeta)
        occ_model.assign(exp_alpha, exp_beta) 
        dm_alpha = exp_alpha.to_dm()
        dm_beta= exp_beta.to_dm()
        
        scf_solver = ic.scf_solver

        if ic.expansion:
            scf_solver(ham, lf, olp, occ_model, exp_alpha, exp_beta)
        else:
            scf_solver(ham, lf, olp, occ_model, dm_alpha, dm_beta)
           
